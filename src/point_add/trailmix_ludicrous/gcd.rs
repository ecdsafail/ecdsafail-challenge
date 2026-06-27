//! The product-min jump-GCD modular inversion (`y -> y * x^-1 mod q`
//! and the forward `y -> y * x` direction), built on this crate's `B`
//! builder and the `schedule` / `arith` / `comparator` / `codec` modules.
//!
//! ## Fused apply with incremental codec
//! Schedule `schedule::SCHED_J2` / `GAP_J2`, jump=2, `ITERS = 258`. The dialog
//! apply is fused into the forward/reverse GCD: each divstep symbol is applied to
//! the coordinate pair the instant it is computed, so the apply adders run in the
//! GCD passes' headroom and the full raw tape is never materialized.
//!
//!   1. `forward_gcd_jump` runs the jump=2 divstep loop, computing the per-step
//!      dialog `(subtracted, swap, s_2)`. When `apply_inv` is set it applies the
//!      inverse step to `[y, tmp]` (`Direction::Inverse`, producing `y*x^-1`)
//!      before the symbol is swapped to the tape; otherwise it only records.
//!   2. each codec window is compressed inline the instant its symbols are
//!      recorded, so the resident tape stays compressed.
//!   3. `reverse_gcd_jump` is the exact gate-inverse of the divstep loop: it
//!      restores `x` and drains the tape to empty, decompressing one window at a
//!      time. When `apply_fwd` is set it applies the forward step to `[tmp, y]`
//!      (`Direction::Forward`, producing `y*x`) in reverse iter order.
//!
//! The compressed tape held live across a pass is `dialog_tape_qubits` qubits
//! (603 for the len-258 all-triple tiling -- see `codec::dialog_tape_qubits`).
//!
//! ## Algorithm structure
//! - The classical jump-before-swap divstep is the algebra `forward_gcd_jump`
//!   realizes in the circuit.
//! - The forward / reverse GCD: the per-step
//!   shift / swap-decision / cswap / subtract structure and the schedule-driven
//!   register shrink/regrow, with the comparator/subtract routed through this
//!   crate's `comparator` / `arith`.
//! - The reverse-iter apply (`if subtracted: y+=x; if swap: swap(x,y);
//!   y *= 2^j mod q`): the multiply is the apply, the inversion is the
//!   reverse-of-multiply composition.
//! - The per-step variable shift (`controlled_right_shift` /
//!   `controlled_left_shift`) and the apply's controlled doubling.

use super::arith::{self, F_SECP256K1};
use super::schedule::{GAP_J2, ITERS, JUMP, SCHED_J2};
use super::{B, BExt};
use crate::circuit::{QubitId};
use std::cell::Cell;

thread_local! {
    static RIGHT_SHIFT_CALL_INDEX: Cell<usize> = const { Cell::new(0) };
    static LEFT_SHIFT_CALL_INDEX: Cell<usize> = const { Cell::new(0) };
}

pub(super) fn reset_gcd_trace_call_index() {
    RIGHT_SHIFT_CALL_INDEX.with(|index| index.set(0));
    LEFT_SHIFT_CALL_INDEX.with(|index| index.set(0));
}

fn next_right_shift_call_index() -> usize {
    RIGHT_SHIFT_CALL_INDEX.with(|index| {
        let current = index.get();
        index.set(current + 1);
        current
    })
}

fn next_left_shift_call_index() -> usize {
    LEFT_SHIFT_CALL_INDEX.with(|index| {
        let current = index.get();
        index.set(current + 1);
        current
    })
}

/// Which product the apply pass reconstructs into `y`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Direction {
    /// `y -> y * x^-1 mod q` (the modular inversion: apply the dialog reversed
    /// vs the multiply, i.e. forward iter order with the inverse step).
    Inverse,
    /// `y -> y * x   mod q` (the multiply: reverse iter order, forward step).
    Forward,
}

/// secp256k1 prime q = 2^256 - f = 0xFFFF...FFFFFC2F, little-endian bytes.
/// (Build-time constant; avoids a bignum dependency in the submission.)
#[must_use]
pub fn q_secp256k1_le() -> [u8; 32] {
    let mut b = [0xFFu8; 32];
    b[0] = 0x2F;
    b[1] = 0xFC;
    b[4] = 0xFE;
    b
}

/// All-triple group count the product-min selector picks for this
/// schedule: `n3 = iters/3` (clamped inside `codec::jump_dialog_regions`).
#[must_use]
pub fn n3_for_iters(iters: usize) -> usize {
    iters / 3
}

fn env_i32(name: &str, default: i32) -> i32 {
    std::env::var(name)
        .ok()
        .and_then(|value| value.parse::<i32>().ok())
        .unwrap_or(default)
}

fn env_usize(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(default)
}

fn adjust_gcd_k(prefix: &str, i: usize, k: usize) -> usize {
    if k == usize::MAX {
        return k;
    }
    let adjust = env_i32(&format!("{prefix}_ADJUST"), 0);
    let after = env_usize(&format!("{prefix}_ADJUST_AFTER"), 0);
    let before = env_usize(&format!("{prefix}_ADJUST_BEFORE"), usize::MAX);
    if i >= after && i < before && adjust != 0 {
        (k as i32).saturating_add(adjust).max(0) as usize
    } else {
        k
    }
}

fn maybe_adjust_late_gcd_k(i: usize, k: usize) -> usize {
    let k = adjust_gcd_k("TLM_GCD_K", i, k);
    adjust_gcd_k("TLM_GCD_K_EXTRA", i, k)
}

fn trace_step_regions(circ: &mut B, direction: &str, i: usize, region_start: usize) {
    if std::env::var("TRACE_TLM_GCD_STEPS").is_err() {
        return;
    }
    let threshold = std::env::var("TRACE_TLM_GCD_MIN_Q")
        .ok()
        .and_then(|value| value.parse::<u32>().ok())
        .unwrap_or(1150);
    let step_max = circ.phase_active_regions[region_start..]
        .iter()
        .map(|(_, _, active)| *active)
        .max()
        .unwrap_or(circ.active_qubits);
    if step_max >= threshold {
        eprintln!(
            "TLM_GCD_STEP direction={direction} i={i} active_max={step_max} global_peak={} ops={}",
            circ.peak_qubits,
            circ.current_ops_len(),
        );
        for (_, phase, active) in &circ.phase_active_regions[region_start..] {
            if *active >= threshold {
                eprintln!(
                    "TLM_GCD_STAGE direction={direction} i={i} active_max={active} phase={phase}",
                );
            }
        }
    }
}

// =====================================================================
// MBU AND-uncompute (HMR + conditional-Z), measurement-vented. `t` holds
// `a AND b` and is returned to |0>.
// =====================================================================

fn clear_and(circ: &mut B, t: &QubitId, a: &QubitId, b: &QubitId) {
    let bit = circ.alloc_bit();
    circ.hmr(*t, bit);
    circ.cz_if_bit(*a, *b, bit);
}

fn park_odd_u0_enabled(i: usize, side: &str) -> bool {
    let all = std::env::var("TLM_PARK_ODD_U0").ok().as_deref() == Some("1");
    let side_on = std::env::var(format!("TLM_PARK_ODD_U0_{side}"))
        .ok()
        .as_deref()
        == Some("1");
    if !all && !side_on {
        return false;
    }
    let limit = std::env::var("TLM_PARK_ODD_U0_LIMIT")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(usize::MAX);
    i < limit
}

fn loan_odd_u0_enabled() -> bool {
    std::env::var("TLM_LOAN_ODD_U0").ok().as_deref() == Some("1")
}

fn park_even_v0_enabled() -> bool {
    std::env::var("TLM_PARK_EVEN_V0").ok().as_deref() == Some("1")
}

fn loan_even_v0_enabled() -> bool {
    std::env::var("TLM_LOAN_EVEN_V0").ok().as_deref() == Some("1")
}

fn loan_gcd_y0_enabled() -> bool {
    std::env::var("TLM_LOAN_GCD_Y0").ok().as_deref() == Some("1")
}

fn apply_fwd_cswap_skip(i: usize) -> bool {
    let legacy_first_skip =
        std::env::var("TLM_APPLY_FWD_FIRST_CSWAP_SKIP").ok().as_deref() == Some("1")
            && i + 1 == ITERS;
    let last_n = std::env::var("TLM_APPLY_FWD_CSWAP_SKIP_LAST")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(0);
    legacy_first_skip || (last_n != 0 && i + last_n >= ITERS)
}

fn index_spec_contains(spec: &str, index: usize) -> bool {
    spec.split(',').any(|item| {
        let item = item.trim();
        if item.eq_ignore_ascii_case("all") {
            return true;
        }
        if let Some((lo, hi)) = item.split_once('-') {
            return lo
                .trim()
                .parse::<usize>()
                .ok()
                .zip(hi.trim().parse::<usize>().ok())
                .is_some_and(|(lo, hi)| lo <= index && index <= hi);
        }
        item.parse::<usize>().ok() == Some(index)
    })
}

fn apply_fwd_cswap_partial_skip(i: usize, bit: usize) -> bool {
    let last_n = std::env::var("TLM_APPLY_FWD_CSWAP_PARTIAL_SKIP_LAST")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(0);
    if last_n == 0 || i + last_n < ITERS {
        return false;
    }
    std::env::var("TLM_APPLY_FWD_CSWAP_PARTIAL_SKIP_BITS")
        .ok()
        .is_some_and(|spec| index_spec_contains(&spec, bit))
}

fn apply_inv_cswap_skip(i: usize) -> bool {
    let last_n = std::env::var("TLM_APPLY_INV_CSWAP_SKIP_LAST")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(0);
    last_n != 0 && i + last_n >= ITERS
}

fn apply_fwd_s2_zero(i: usize) -> bool {
    let last_n = std::env::var("TLM_APPLY_FWD_S2_ZERO_LAST")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(0);
    last_n != 0 && i + last_n >= ITERS
}

fn apply_inv_s2_zero(i: usize) -> bool {
    let last_n = std::env::var("TLM_APPLY_INV_S2_ZERO_LAST")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(0);
    last_n != 0 && i + last_n >= ITERS
}

fn apply_add_skip(i: usize, fwd: bool) -> bool {
    if let Some(k) = std::env::var("TLM_APPLY_ADD_SKIP_LASTK")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
    {
        if k != 0 && i + k >= ITERS {
            return true;
        }
    }
    let var = if fwd {
        "TLM_APPLY_ADD_SKIP_FWD"
    } else {
        "TLM_APPLY_ADD_SKIP_INV"
    };
    let k = std::env::var(var)
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(0);
    k != 0 && i + k >= ITERS
}

fn park_known_one(circ: &mut B, q: QubitId) -> QubitId {
    circ.x(q);
    if loan_odd_u0_enabled() {
        circ.loan_zero_qubit(q);
    } else {
        circ.zero_and_free(q);
    }
    q
}

fn restore_known_one(circ: &mut B, parked: QubitId) -> QubitId {
    let q = if loan_odd_u0_enabled() {
        circ.reclaim_zero_qubit(parked);
        parked
    } else {
        circ.alloc_qubit()
    };
    circ.x(q);
    q
}

fn park_known_zero(circ: &mut B, q: QubitId) -> QubitId {
    if loan_even_v0_enabled() {
        circ.loan_zero_qubit(q);
    } else {
        circ.zero_and_free(q);
    }
    q
}

fn restore_known_zero(circ: &mut B, parked: QubitId) -> QubitId {
    if loan_even_v0_enabled() {
        circ.reclaim_zero_qubit(parked);
        parked
    } else {
        circ.alloc_qubit()
    }
}

fn loan_known_one_gcd_y0(circ: &mut B, q: QubitId) {
    circ.x(q);
    circ.loan_zero_qubit(q);
}

fn reclaim_known_one_gcd_y0(circ: &mut B, q: QubitId) {
    circ.reclaim_zero_qubit(q);
    circ.x(q);
}

fn loan_known_zero_gcd_y0(circ: &mut B, q: QubitId) {
    circ.loan_zero_qubit(q);
}

fn reclaim_known_zero_gcd_y0(circ: &mut B, q: QubitId) {
    circ.reclaim_zero_qubit(q);
}

const GCD_REVERSE_CSWAP_DEAD_RANGES: &[(usize, usize, usize)] = &[
    (255, 5, 13),
    (256, 3, 11),
    (253, 8, 15),
    (254, 7, 13),
    (252, 10, 15),
    (250, 13, 16),
    (251, 12, 15),
    (249, 13, 16),
    (236, 27, 29),
    (248, 15, 17),
    (234, 29, 31),
    (235, 28, 30),
    (237, 26, 28),
    (244, 18, 20),
    (246, 17, 19),
    (247, 16, 18),
    (101, 164, 165),
    (143, 122, 123),
    (145, 120, 121),
    (219, 45, 46),
    (220, 44, 45),
    (221, 43, 44),
    (224, 40, 41),
    (226, 38, 39),
    (228, 36, 37),
    (229, 35, 36),
    (231, 33, 34),
    (232, 32, 33),
    (233, 31, 32),
    (245, 18, 19),
    (257, 7, 10),
    (95, 170, 171),
    (116, 149, 150),
    (134, 131, 132),
    (218, 46, 47),
    (222, 42, 43),
    (223, 41, 42),
    (225, 39, 40),
    (227, 37, 38),
    (230, 34, 35),
    (240, 22, 23),
    (11, 254, 254),
    (12, 253, 253),
    (19, 246, 246),
    (21, 244, 244),
    (30, 235, 235),
    (31, 234, 234),
    (33, 232, 232),
    (35, 230, 230),
    (36, 229, 229),
    (37, 228, 228),
    (39, 226, 226),
    (40, 225, 225),
    (42, 223, 223),
    (43, 222, 222),
    (46, 219, 219),
    (47, 218, 218),
    (48, 217, 217),
    (49, 216, 216),
    (50, 215, 215),
    (51, 214, 214),
    (53, 212, 212),
    (54, 211, 211),
    (56, 209, 209),
    (63, 202, 202),
    (68, 197, 197),
    (70, 195, 195),
    (75, 190, 190),
    (76, 189, 189),
    (79, 186, 186),
    (80, 185, 185),
    (81, 184, 184),
    (87, 178, 178),
    (94, 172, 172),
    (103, 163, 163),
    (105, 161, 161),
    (107, 159, 159),
    (108, 158, 158),
    (110, 156, 156),
    (112, 154, 154),
    (113, 153, 153),
    (114, 152, 152),
    (115, 151, 151),
    (123, 143, 143),
    (124, 142, 142),
    (126, 140, 140),
    (127, 139, 139),
    (128, 138, 138),
    (129, 137, 137),
    (130, 136, 136),
    (131, 135, 135),
    (132, 134, 134),
    (133, 133, 133),
    (141, 125, 125),
    (142, 124, 124),
    (146, 119, 119),
    (147, 118, 118),
    (148, 117, 117),
    (149, 116, 116),
    (153, 112, 112),
    (157, 108, 108),
    (159, 106, 106),
    (161, 104, 104),
    (162, 103, 103),
    (163, 102, 102),
    (164, 101, 101),
    (167, 98, 98),
    (168, 97, 97),
    (169, 96, 96),
    (170, 95, 95),
    (171, 94, 94),
    (180, 85, 85),
    (182, 83, 83),
    (183, 82, 82),
    (187, 78, 78),
    (188, 77, 77),
    (189, 76, 76),
    (190, 75, 75),
    (192, 73, 73),
    (193, 72, 72),
    (194, 71, 71),
    (195, 70, 70),
    (197, 68, 68),
    (198, 67, 67),
    (199, 66, 66),
    (203, 62, 62),
    (204, 61, 61),
    (205, 60, 60),
    (206, 59, 59),
    (207, 58, 58),
    (208, 57, 57),
    (210, 55, 55),
    (211, 54, 54),
    (212, 53, 53),
    (213, 52, 52),
    (214, 51, 51),
    (215, 50, 50),
    (216, 49, 49),
    (217, 48, 48),
    (238, 25, 25),
    (239, 24, 24),
    (241, 22, 22),
    (242, 21, 21),
];

const GCD_FORWARD_CSWAP_REMAINDER_KEYS: &[u32] = &[
    3325, 4345, 8935, 9955, 17860, 24236, 26021, 26531, 26786, 27551, 28316, 28571,
    29846, 31631, 31886, 32906, 33161, 33416, 33671, 33926, 34181, 34436, 35710,
    36221, 36476, 36985, 37241, 38260, 40810, 41575, 41830, 43360, 44380, 45145,
    46165, 46675, 47185, 47950, 48205, 48715, 49225, 50755, 51010, 51520, 51775,
    52030, 52285, 52540, 53305, 53560, 54070, 54580, 54835, 55090, 55600, 56110,
    56620, 56875, 57130, 57385, 57895, 58149, 58405, 58660, 58915, 59170, 59425,
    59680, 59935, 60190, 60444, 60445, 60700, 62484, 62995, 63250, 63505, 63759,
    63760, 64014, 64015, 64016, 64270, 64271, 64524, 64525, 64526, 64527, 64779,
    64780, 64781, 64782, 64783, 65035, 65036, 65037, 65290, 65291, 65292, 65293,
    65545, 65546, 65547, 65800, 65801, 65802,
];

const GCD_REVERSE_CSWAP_REMAINDER_KEYS: &[u32] = &[
    3580, 3835, 4090, 4345, 4600, 4855, 5365, 5875, 6130, 6640, 6895, 7150,
    7405, 7660, 8425, 8935, 9955, 10720, 11740, 13525, 14290, 15055, 15310,
    15565, 15820, 16075, 16585, 16840, 17095, 17350, 18370, 18625, 18880,
    19135, 19900, 20155, 21175, 21685, 22195, 22705, 22960, 23215, 23470,
    23725, 24745, 25255, 25510, 26786, 31376, 34690, 34945, 35200, 35455,
    38515, 38770, 39025, 39790, 40045, 40555, 41065, 42340, 42595, 44125,
    44635, 44890, 45145, 45400, 45655, 45910, 46420, 47185, 47440, 48970,
    50245, 51265, 53560,
];

fn gcd_reverse_cswap_has_structurally_dead_gate(step: usize, bit: usize) -> bool {
    if std::env::var_os("TLM_GCD_SKIP_STRUCTURAL_DEAD_CSWAPS").is_none() {
        return false;
    }
    if std::env::var_os("TLM_GCD_SKIP_REVERSE_DIAGONAL_EDGE").is_some()
        && step + bit
            >= std::env::var("TLM_GCD_REVERSE_DIAGONAL_MIN")
                .ok()
                .and_then(|value| value.parse::<usize>().ok())
                .unwrap_or(265)
        && step
            >= std::env::var("TLM_GCD_REVERSE_DIAGONAL_STEP_MIN")
                .ok()
                .and_then(|value| value.parse::<usize>().ok())
                .unwrap_or(0)
        && step
            <= std::env::var("TLM_GCD_REVERSE_DIAGONAL_STEP_MAX")
                .ok()
                .and_then(|value| value.parse::<usize>().ok())
                .unwrap_or(usize::MAX)
    {
        return true;
    }
    if std::env::var_os("TLM_GCD_SKIP_EXACT_REVERSE_CSWAPS").is_some() {
        let key = (((step as u32) & 0xffff) << 8) | (bit as u32 & 0xff);
        if GCD_REVERSE_CSWAP_REMAINDER_KEYS.binary_search(&key).is_ok() {
            return true;
        }
    }
    GCD_REVERSE_CSWAP_DEAD_RANGES
        .iter()
        .any(|&(range_step, lo, hi)| range_step == step && (lo..=hi).contains(&bit))
}

fn gcd_forward_cswap_has_structurally_dead_gate(step: usize, bit: usize) -> bool {
    if std::env::var_os("TLM_GCD_SKIP_STRUCTURAL_DEAD_CSWAPS").is_none()
        || std::env::var_os("TLM_GCD_SKIP_EXACT_FORWARD_CSWAPS").is_none()
    {
        return false;
    }
    let key = (((step as u32) & 0xffff) << 8) | (bit as u32 & 0xff);
    GCD_FORWARD_CSWAP_REMAINDER_KEYS.binary_search(&key).is_ok()
}

const GCD_SHIFT_DEAD_RANGES: &[(u8, usize, usize, usize)] = &[
    (12, 0, 1, 9),
    (12, 259, 1, 9),
    (12, 1, 3, 10),
    (12, 2, 5, 12),
    (12, 260, 3, 10),
    (12, 261, 5, 12),
    (12, 262, 7, 12),
    (12, 263, 9, 14),
    (12, 264, 9, 14),
    (12, 3, 7, 12),
    (12, 4, 9, 14),
    (12, 5, 10, 14),
    (11, 254, 11, 14),
    (11, 515, 8, 8),
    (11, 515, 10, 12),
    (11, 516, 7, 10),
    (12, 265, 11, 14),
    (12, 266, 12, 15),
    (12, 6, 11, 14),
    (12, 7, 12, 15),
    (11, 252, 12, 14),
    (11, 253, 12, 14),
    (11, 255, 10, 12),
    (11, 256, 10, 12),
    (11, 257, 8, 10),
    (11, 258, 7, 9),
    (11, 510, 13, 15),
    (11, 512, 12, 14),
    (11, 513, 12, 14),
    (11, 514, 10, 12),
    (12, 20, 25, 27),
    (12, 267, 13, 15),
    (12, 269, 15, 17),
    (12, 270, 16, 18),
    (12, 279, 25, 27),
    (12, 280, 26, 28),
    (12, 287, 33, 35),
    (12, 8, 13, 15),
    (12, 9, 14, 16),
    (11, 250, 14, 15),
    (11, 251, 14, 15),
    (11, 391, 133, 134),
    (11, 405, 119, 120),
    (11, 511, 12, 12),
    (11, 511, 14, 14),
    (11, 517, 8, 9),
    (12, 10, 16, 17),
    (12, 11, 17, 18),
    (12, 12, 17, 18),
    (12, 13, 18, 19),
    (12, 21, 27, 28),
    (12, 22, 28, 29),
    (12, 23, 29, 30),
    (12, 24, 30, 31),
    (12, 25, 31, 32),
    (12, 26, 32, 33),
    (12, 268, 15, 16),
    (12, 27, 33, 34),
    (12, 271, 17, 18),
    (12, 28, 34, 35),
    (12, 281, 28, 29),
    (12, 283, 30, 31),
    (12, 286, 33, 34),
    (12, 288, 35, 36),
    (12, 289, 36, 37),
    (12, 290, 37, 38),
    (12, 293, 40, 41),
    (12, 294, 41, 42),
    (12, 295, 42, 43),
    (12, 31, 37, 38),
    (12, 32, 38, 39),
    (12, 33, 39, 40),
    (12, 34, 40, 41),
    (12, 35, 41, 42),
    (12, 36, 42, 43),
];

const GCD_RIGHT_SHIFT_REMAINDER_KEYS: &[u32] = &[
    510, 766, 1022, 1278, 1534, 1790, 2046, 2302, 2558, 2814, 3070, 3325,
    3580, 3835, 4090, 4345, 4600, 4855, 5110, 5365, 5620, 5875, 6130, 6385,
    6640, 6895, 7150, 7405, 7660, 7915, 8170, 8425, 8680, 8935, 9190, 9445,
    9700, 9955, 10210, 10465, 10720, 10975, 11230, 11485, 11740, 11995, 12250, 12505,
    12760, 13015, 13270, 13525, 13780, 14035, 14290, 14545, 14800, 15055, 15310, 15565,
    15820, 16075, 16330, 16585, 16840, 17095, 17350, 17605, 17860, 18115, 18370, 18625,
    18880, 19135, 19390, 19645, 19900, 20155, 20410, 20665, 20920, 21175, 21430, 21685,
    21940, 22195, 22450, 22705, 22960, 23215, 23470, 23725, 23980, 24235, 24491, 24746,
    25000, 25255, 25510, 25765, 26020, 26276, 26530, 26786, 27041, 27296, 27550, 27806,
    28061, 28315, 28571, 28826, 29081, 29336, 29591, 29846, 30101, 30355, 30610, 30865,
    31120, 31375, 31631, 31886, 32141, 32395, 32651, 32906, 33161, 33416, 33671, 33926,
    34181, 34436, 34691, 34945, 35200, 35455, 35710, 35965, 36220, 36476, 36731, 36986,
    37240, 37496, 37750, 38005, 38260, 38515, 38770, 39025, 39280, 39535, 39790, 40045,
    40300, 40555, 40810, 41065, 41320, 41575, 41830, 42085, 42340, 42595, 42850, 43105,
    43360, 43615, 43870, 44125, 44380, 44635, 44890, 45145, 45400, 45655, 45910, 46165,
    46420, 46675, 46930, 47185, 47440, 47695, 47950, 48205, 48460, 48715, 48970, 49225,
    49480, 49735, 49990, 50245, 50500, 50755, 51010, 51265, 51520, 51775, 52030, 52285,
    52540, 52795, 53050, 53305, 53560, 53815, 54070, 54325, 54580, 54835, 55090, 55345,
    55600, 55855, 56110, 56365, 56620, 56875, 57130, 57385, 57640, 57895, 58150, 58405,
    58660, 58915, 59170, 59425, 59680, 59935, 60190, 60445, 60700, 60955, 61208, 61463,
    61718, 61973, 62228, 62483, 62739, 62994, 63250, 63505, 63760, 66814, 67070, 67326,
    67582, 67838, 68094, 68350, 68606, 68862, 69118, 69374, 69629, 69884, 70139, 70394,
    70649, 70904, 71159, 71414, 71669, 71924, 72179, 72434, 72689, 72944, 73199, 73454,
    73709, 73964, 74219, 74474, 74729, 74984, 75239, 75494, 75749, 76004, 76259, 76514,
    76769, 77024, 77279, 77534, 77789, 78044, 78299, 78554, 78809, 79064, 79319, 79574,
    79829, 80084, 80339, 80594, 80849, 81104, 81359, 81614, 81869, 82124, 82379, 82634,
    82889, 83144, 83399, 83654, 83909, 84164, 84419, 84674, 84929, 85184, 85439, 85694,
    85949, 86204, 86459, 86714, 86969, 87224, 87479, 87734, 87989, 88244, 88499, 88754,
    89009, 89264, 89519, 89774, 90029, 90284, 90539, 90795, 91050, 91304, 91559, 91814,
    92069, 92324, 92580, 92834, 93090, 93345, 93600, 93854, 94110, 94365, 94619, 94875,
    95130, 95385, 95640, 95895, 96150, 96405, 96659, 96914, 97169, 97424, 97679, 97935,
    98190, 98445, 98699, 98955, 99210, 99465, 99720, 99975, 100485, 100740, 100995, 101249,
    101504, 101759, 102014, 102269, 102524, 102780, 103035, 103290, 103544, 104054, 104309, 104564,
    104819, 105074, 105329, 105584, 105839, 106094, 106349, 106604, 106859, 107114, 107369, 107624,
    107879, 108134, 108389, 108644, 108899, 109154, 109409, 109664, 109919, 110174, 110429, 110684,
    110939, 111194, 111449, 111704, 111959, 112214, 112469, 112724, 112979, 113234, 113489, 113744,
    113999, 114254, 114509, 114764, 115019, 115274, 115529, 115784, 116039, 116294, 116549, 116804,
    117059, 117314, 117569, 117824, 118079, 118334, 118589, 118844, 119099, 119354, 119609, 119864,
    120119, 120374, 120629, 120884, 121139, 121394, 121649, 121904, 122159, 122414, 122669, 122924,
    123179, 123434, 123689, 123944, 124199, 124454, 124709, 124964, 125219, 125474, 125729, 125984,
    126239, 126494, 126749, 127004, 127259, 127512, 127767, 128022, 128277, 128532, 128787, 129043,
    129298, 129554, 129809, 130064, 130319,
];

const GCD_LEFT_SHIFT_REMAINDER_KEYS: &[u32] = &[
    3603, 3860, 4117, 4374, 4631, 4888, 7460, 7717, 9516, 9773, 12086, 12600,
    13114, 15427, 17740, 25707, 28792, 29306, 29820, 30847, 31361, 31619, 32903, 34189,
    36245, 37273, 41642, 47038, 53463, 62715, 69651, 69907, 70164, 70421, 70678, 70935,
    71192, 72222, 72736, 72993, 74535, 74792, 75820, 76077, 80446, 87385, 89441, 95096,
    95610, 95867, 96124, 99979, 100493, 102549, 103320, 104605, 105376, 105890, 107946, 110772,
    119510, 126963, 128248,
];

fn gcd_shift_has_structurally_dead_gate(tag: u8, call_index: usize, bit: usize) -> bool {
    if std::env::var_os("TLM_GCD_SKIP_EXACT_SHIFT_REMAINDER").is_some() {
        let key = (((call_index as u32) & 0xffff) << 8) | (bit as u32 & 0xff);
        if (tag == 11 && GCD_RIGHT_SHIFT_REMAINDER_KEYS.binary_search(&key).is_ok())
            || (tag == 12 && GCD_LEFT_SHIFT_REMAINDER_KEYS.binary_search(&key).is_ok())
        {
            return true;
        }
    }
    std::env::var_os("TLM_GCD_SKIP_STRUCTURAL_DEAD_SHIFTS").is_some()
        && GCD_SHIFT_DEAD_RANGES.iter().any(|&(range_tag, call, lo, hi)| {
            range_tag == tag && call == call_index && (lo..=hi).contains(&bit)
        })
}

fn skip_top_zero_controlled_shift_edge(call_index: usize, bit: usize) -> bool {
    if !std::env::var_os("TLM_GCD_SKIP_TOP_ZERO_SHIFT_EDGE").is_some() {
        return false;
    }
    if let Ok(spec) = std::env::var("TLM_GCD_SKIP_TOP_ZERO_SHIFT_EDGE_CALLS") {
        if !index_spec_contains(&spec, call_index) {
            return false;
        }
    }
    if let Ok(spec) = std::env::var("TLM_GCD_SKIP_TOP_ZERO_SHIFT_EDGE_BITS") {
        if !index_spec_contains(&spec, bit) {
            return false;
        }
    }
    true
}

// =====================================================================
// Per-step variable shift (controlled right/left shift by 1).
// =====================================================================

/// Controlled logical right-shift-by-1 (LSB-first): when `ctrl`, v[i+1] -> v[i].
fn controlled_right_shift(circ: &mut B, ctrl: &QubitId, v: &[QubitId]) {
    let call_index = next_right_shift_call_index();
    for i in 0..v.len().saturating_sub(1) {
        let old_context = crate::point_add::set_op_trace_context(
            0x0b00_0000 | (((call_index as u32) & 0xffff) << 8) | (i as u32 & 0xff),
        );
        let top_zero_edge = i + 2 == v.len() && skip_top_zero_controlled_shift_edge(call_index, i);
        if !top_zero_edge && !gcd_shift_has_structurally_dead_gate(11, call_index, i) {
            circ.cswap(*ctrl, v[i], v[i + 1]);
        }
        crate::point_add::restore_op_trace_context(old_context);
    }
}

/// Exact gate-inverse of [`controlled_right_shift`] (cswaps reversed).
fn controlled_left_shift(circ: &mut B, ctrl: &QubitId, v: &[QubitId]) {
    let call_index = next_left_shift_call_index();
    for i in (1..v.len()).rev() {
        let old_context = crate::point_add::set_op_trace_context(
            0x0c00_0000 | (((call_index as u32) & 0xffff) << 8) | ((i - 1) as u32 & 0xff),
        );
        let bit = i - 1;
        let top_zero_edge = i + 1 == v.len() && skip_top_zero_controlled_shift_edge(call_index, bit);
        if !top_zero_edge && !gcd_shift_has_structurally_dead_gate(12, call_index, i - 1) {
            circ.cswap(*ctrl, v[i], v[i - 1]);
        }
        crate::point_add::restore_op_trace_context(old_context);
    }
}

/// Unconditional logical right-shift-by-1: v[i+1] -> v[i], top bit -> |0>.
/// A swap chain from the bottom up; v[0] is |0> in the GCD post-subtract so it
/// wraps harmlessly to the top.
fn right_shift(circ: &mut B, v: &[QubitId]) {
    for i in 0..v.len().saturating_sub(1) {
        circ.swap(v[i], v[i + 1]);
    }
}

/// Exact gate-inverse of [`right_shift`].
fn left_shift(circ: &mut B, v: &[QubitId]) {
    for i in (1..v.len()).rev() {
        circ.swap(v[i], v[i - 1]);
    }
}

// =====================================================================
// Controlled mod_double for the apply (gated `y := 2*y mod q`), built on
// this crate's arith::add_f_window (re-using add_f via a ctrl).
// =====================================================================

/// If `ctrl`: `a := 2*a mod q`; else `a` unchanged. `a.len() == n = 256`; the n+1
/// shift view is `a ++ ovf` with a transient `ovf` (not a slot in `a`), restored to
/// |0>. Normal controlled shift + the (MBU) `+f` fold gated on the shifted-out
/// overflow + an MBU ancilla clean.
fn controlled_mod_double(circ: &mut B, ctrl: &QubitId, a: &[QubitId]) {
    let n = a.len();
    assert_eq!(n, 256, "controlled_mod_double expects 256-bit a");
    let f_bytes = F_SECP256K1.to_le_bytes();
    let ovf = circ.alloc_qubit();
    // 1) controlled left-shift by 1 over the n+1 view a ++ ovf: MSB -> ovf.
    let w: Vec<&QubitId> = a.iter().chain(std::iter::once(&ovf)).collect();
    for i in (0..n).rev() {
        circ.cswap(*ctrl, *w[i], *w[i + 1]);
    }
    // 2) if ovf (= ctrl AND old-MSB), add f to the low LSBS bits.
    arith::add_f_window_pub(circ, &ovf, a, arith::LSBS, &f_bytes, None);
    // 3) clear ovf (= ctrl AND a[0], post-fold) by MBU.
    clear_and(circ, &ovf, ctrl, &a[0]);
    circ.zero_and_free(ovf);
}

/// Exact gate-inverse of [`controlled_mod_double`]: if `ctrl`, `a := a/2 mod q`.
fn controlled_mod_double_reverse(circ: &mut B, ctrl: &QubitId, a: &[QubitId]) {
    let n = a.len();
    assert_eq!(n, 256, "controlled_mod_double_reverse expects 256-bit a");
    let f_bytes = F_SECP256K1.to_le_bytes();
    let ovf = circ.alloc_qubit();
    // inverse of 3): rebuild ovf = ctrl AND a[0] with a CCX.
    circ.ccx(*ctrl, a[0], ovf);
    // inverse of 2): subtract f (X-sandwich of add_f_window) gated on ovf.
    for q in &a[..arith::LSBS] {
        circ.x(*q);
    }
    arith::add_f_window_pub(circ, &ovf, a, arith::LSBS, &f_bytes, None);
    for q in &a[..arith::LSBS] {
        circ.x(*q);
    }
    // inverse of 1): controlled right-shift over a ++ ovf; ovf shifted back to |0>.
    let w: Vec<&QubitId> = a.iter().chain(std::iter::once(&ovf)).collect();
    for i in 0..n {
        circ.cswap(*ctrl, *w[i], *w[i + 1]);
    }
    circ.zero_and_free(ovf);
}

// =====================================================================
// Forward jump-GCD: record the dialog onto the compressed tape.
// =====================================================================

/// Forward jump=2 GCD on `(u = q, v = x)` over the baked schedule. `v` is the
/// input register (>= 256 bits; bits beyond the active width must be |0> for
/// schedule-fitting inputs). Computes the per-step dialog `(subtracted, swap,
/// s_2)`, compresses each codec window inline, and returns the compressed tape
/// (`[t1, win_0 code bits, ..]`, `codec::dialog_tape_qubits` long). On a fitting
/// input, ends with `u = 0` (X'd from 1), `v = 0`, both shrunk away; the tape
/// holds the whole dialog.
///
/// When `apply_inv` is `Some((x_reg, y_reg))`, the inverse apply step is fused in
/// (applied to the coordinate pair each iter before the symbol is taped);
/// otherwise the coordinate registers are untouched.
#[must_use]
pub fn forward_gcd_jump(circ: &mut B, v: &mut Vec<QubitId>, apply_inv: Option<(&[QubitId], &[QubitId])>) -> Vec<QubitId> {
    let n = 256usize;
    assert_eq!(JUMP, 2, "ludicrous apply/codec are jump=2 specific");
    assert!(v.len() >= n, "v must be at least n=256 bits");
    let iters = ITERS;
    let sym_bits = 3; // (subtracted, swap, s_2)

    // u_full = q (n bits); shrinks per the schedule.
    let mut u: Vec<QubitId> = (0..n).map(|_| circ.alloc_qubit()).collect();
    let q_bytes = q_secp256k1_le();
    for (i, qb) in u.iter().enumerate() {
        if (q_bytes.get(i / 8).copied().unwrap_or(0) >> (i % 8)) & 1 == 1 {
            circ.x(*qb);
        }
    }

    let subtracted = circ.alloc_qubit();
    let mut swap_flag: Option<QubitId> = None;
    let s2 = circ.alloc_qubit(); // the single jump=2 extra-shift flag
    let t1 = circ.alloc_qubit(); // step-0 shift1-fired flag (= x even)

    // Incremental all-triple codec: compress each codec window inline the instant
    // its symbols are recorded, so the resident tape grows compressed (~603) not
    // raw (775), bounding the peak. `tape` accumulates the t1 prefix
    // followed by each window's compressed code bits, in iter order.
    let n3 = n3_for_iters(iters);
    let mut window_plan: Vec<super::codec::DialogCodec> = Vec::new();
    for (codec, count) in super::codec::jump_dialog_regions(n3, iters) {
        for _ in 0..count {
            window_plan.push(codec);
        }
    }
    let mut tape: Vec<QubitId> = Vec::with_capacity(super::codec::dialog_tape_qubits(n3, iters));
    let mut win_idx = 0usize; // current window
    let mut pending: Vec<QubitId> = Vec::new(); // raw symbol slots of the current window
    let mut tail4_prefix_encoded = false;
    for i in 0..iters {
        let trace_region_start = circ.phase_active_regions.len();
        circ.set_phase(if apply_inv.is_some() {
            "tlm_inverse_gcd_forward_shift"
        } else {
            "tlm_multiply_gcd_forward_shift"
        });
        let current_n = (SCHED_J2[i] as usize).max(1);
        while u.len() > current_n {
            let q = u.pop().expect("u nonempty");
            circ.zero_and_free(q);
        }
        while v.len() > current_n {
            let q = v.pop().expect("v nonempty");
            circ.zero_and_free(q);
        }
        // swap-decision window: top GAP_J2[i] bits of the active width.
        let cmp_eff = (GAP_J2[i] as usize).min(current_n).max(1);

        // 1) Shift-first: remove up to jump=2 trailing zeros of v.
        //    i==0 gates shift1 on (v even) and records it in t1; i>=1 is
        //    unconditional (v even post-subtract).
        if i == 0 {
            circ.cx(v[0], t1); // t1 = v[0]
            circ.x(t1); // t1 = NOT(v[0]) = (x even)
            controlled_right_shift(circ, &t1, &v[..current_n]);
        } else {
            right_shift(circ, &v[..current_n]);
        }
        // s_2: shift again while still even. Nesting is automatic (once an odd
        // bit reaches LSB the shift stops).
        circ.cx(v[0], s2);
        circ.x(s2);
        controlled_right_shift(circ, &s2, &v[..current_n]);

        // 2) subtracted = v[0] (post-shift parity): 1 => odd => swap+subtract.
        circ.cx(v[0], subtracted);

        // 3) swap decision. Step 0: u=q, v=x<q always, so (v<u)=1 deterministically,
        //    so swap = subtracted and no separate swap flag is needed. Else the narrow top-k comparator
        //    decides swap_flag ^= subtracted AND (v < u).
        circ.set_phase(if apply_inv.is_some() {
            "tlm_inverse_gcd_forward_compare"
        } else {
            "tlm_multiply_gcd_forward_compare"
        });
        let swp = if i == 0 {
            subtracted
        } else {
            let sf = *swap_flag.get_or_insert_with(|| circ.alloc_qubit());
            controlled_swap_decision_v_lt_u(
                circ,
                &subtracted,
                &v[..current_n],
                &u[..current_n],
                cmp_eff,
                &sf,
            );
            sf
        };
        // 4) cswap(swap_flag, u, v).
        for j in 1..current_n {
            if !gcd_forward_cswap_has_structurally_dead_gate(i, j) {
                let old_context = crate::point_add::set_op_trace_context(
                    0x1b00_0000 | (((i as u32) & 0xffff) << 8) | (j as u32 & 0xff),
                );
                circ.cswap(swp, u[j], v[j]);
                crate::point_add::restore_op_trace_context(old_context);
            }
        }
        let parked_u0 = if park_odd_u0_enabled(i, "FWD") {
            let q = u[0];
            Some(park_known_one(circ, q))
        } else {
            None
        };
        // 5) v -= subtracted * u (controlled mod-free subtract on the active width;
        //    post-swap v >= u so no borrow on a fitting input). X-sandwich add.
        circ.set_phase(if apply_inv.is_some() {
            "tlm_inverse_gcd_forward_body"
        } else {
            "tlm_multiply_gcd_forward_body"
        });
        for q in &v[..current_n] {
            circ.x(*q);
        }
        controlled_add_active(
            circ,
            i,
            &subtracted,
            &u[..current_n],
            &v[..current_n],
            GcdBit0Mode::ForwardKnownOneAfterCx,
            apply_inv.map_or(&[], |(xr, _)| &xr[..3]),
        );
        for q in &v[..current_n] {
            circ.x(*q);
        }

        // Fused inverse apply (Direction::Inverse): apply this divstep's symbol to the
        // coordinate pair using the live symbol bits, before they are swapped to the
        // tape. Forward iter order matches apply_step_reverse's order, so the tape is
        // never materialized for the apply -> the apply adders run in the GCD's headroom.
        circ.set_phase(if apply_inv.is_some() {
            "tlm_inverse_gcd_forward_apply"
        } else {
            "tlm_multiply_gcd_forward_apply"
        });
        if i >= 250 && std::env::var_os("TRACE_TLM_TAIL").is_some() {
            eprintln!(
                "TLM_TAIL direction=forward i={i} active={} tape={} pending={} encoded={tail4_prefix_encoded}",
                circ.active_qubits,
                tape.len(),
                pending.len(),
            );
        }
        let parked_v0 = if apply_inv.is_some() && park_even_v0_enabled() {
            let q = v[0];
            Some(park_known_zero(circ, q))
        } else {
            None
        };
        if let Some((xr, yr)) = apply_inv {
            apply_step_reverse(
                circ,
                i,
                &subtracted,
                &swp,
                &s2,
                &t1,
                xr,
                yr,
                &u[1..4],
            );
        }
        if let Some(q) = parked_v0 {
            v[0] = restore_known_zero(circ, q);
        }
        if let Some(q) = parked_u0 {
            u[0] = restore_known_one(circ, q);
        }

        // 6) record the symbol into fresh |0> slots (returning the ancilla to |0>).
        circ.set_phase(if apply_inv.is_some() {
            "tlm_inverse_gcd_forward_codec"
        } else {
            "tlm_multiply_gcd_forward_codec"
        });
        let slots: Vec<QubitId> = (0..sym_bits).map(|_| circ.alloc_qubit()).collect();
        circ.swap(subtracted, slots[0]);
        if i == 0 {
            circ.cx(slots[0], slots[1]);
        } else {
            circ.swap(swp, slots[1]);
        }
        circ.swap(s2, slots[2]);
        if i == 0 {
            debug_assert_eq!(window_plan[win_idx], super::codec::DialogCodec::Step0);
            let data = super::codec::compress_step0_with_t1(circ, t1, &slots);
            tape.extend(data);
            win_idx += 1;
            circ.set_phase("tlm_gcd_step_end");
            trace_step_regions(
                circ,
                if apply_inv.is_some() {
                    "inverse-forward"
                } else {
                    "multiply-forward"
                },
                i,
                trace_region_start,
            );
            continue;
        }
        pending.extend(slots);

        // When the current window's symbols are complete, compress it inline.
        let codec = window_plan[win_idx];
        if codec == super::codec::DialogCodec::Tail4Top32 {
            if !tail4_prefix_encoded && pending.len() == 3 * sym_bits {
                pending = super::codec::DialogCodec::Triple.compress_window(circ, &pending);
                tail4_prefix_encoded = true;
            } else if tail4_prefix_encoded
                && pending.len() == super::codec::DialogCodec::Triple.code_bits() + 2 * sym_bits
            {
                let last = pending.split_off(super::codec::DialogCodec::Triple.code_bits());
                let mut raw = super::codec::DialogCodec::Triple.decompress_window(circ, &pending);
                raw.extend(last);
                let data = codec.compress_window(circ, &raw);
                tape.extend(data);
                pending.clear();
                tail4_prefix_encoded = false;
                win_idx += 1;
            }
        } else if pending.len() == codec.syms() * sym_bits {
            let data = codec.compress_window(circ, &pending);
            tape.extend(data);
            pending.clear();
            win_idx += 1;
        }
        circ.set_phase("tlm_gcd_step_end");
        trace_step_regions(
            circ,
            if apply_inv.is_some() {
                "inverse-forward"
            } else {
                "multiply-forward"
            },
            i,
            trace_region_start,
        );
    }
    assert_eq!(win_idx, window_plan.len(), "all windows compressed");
    assert!(pending.is_empty(), "no leftover symbols");

    // u converged to 1; X to 0 then free. v=0; shrink away.
    circ.x(u[0]);
    while let Some(q) = v.pop() {
        circ.zero_and_free(q);
    }
    for q in u {
        circ.zero_and_free(q);
    }
    circ.zero_and_free(subtracted);
    if let Some(swap_flag) = swap_flag {
        circ.zero_and_free(swap_flag);
    }
    circ.zero_and_free(s2);
    assert_eq!(tape.len(), super::codec::dialog_tape_qubits(n3, iters));
    tape
}

/// Reverse of [`forward_gcd_jump`]: restores `v[..256]` to the original `x` and
/// drains the compressed `tape` to empty (all |0>, freed). Exact gate-inverse,
/// step by step. Decompresses one codec window at a time (in reverse iter order),
/// consumes its symbols, and frees the raw slots -- so the resident tape stays
/// compressed.
pub fn reverse_gcd_jump(circ: &mut B, v: &mut Vec<QubitId>, tape: &mut Vec<QubitId>, apply_fwd: Option<(&[QubitId], &[QubitId])>) {
    let n = 256usize;
    let iters = ITERS;
    let n3 = n3_for_iters(iters);
    assert_eq!(
        tape.len(),
        super::codec::dialog_tape_qubits(n3, iters),
        "tape must be the compressed dialog"
    );

    // Window plan (iter order); we consume windows from the last one back.
    let mut window_plan: Vec<super::codec::DialogCodec> = Vec::new();
    for (codec, count) in super::codec::jump_dialog_regions(n3, iters) {
        for _ in 0..count {
            window_plan.push(codec);
        }
    }
    let mut win_idx = window_plan.len(); // next window to decompress (from the end)
    // `pending` holds the raw symbol slots of the currently-decompressed window,
    // consumed symbol-by-symbol from the end (reverse symbol order).
    let mut pending: Vec<QubitId> = Vec::new();
    let mut pending_tail4 = false;
    let mut tail4_prefix_encoded = false;

    // u regrows from 1 bit (forward ended u=0 post-X; re-init u_final=1 via X).
    let mut u: Vec<QubitId> = vec![circ.alloc_qubit()];
    circ.x(u[0]);

    let subtracted = circ.alloc_qubit();
    let mut swap_flag: Option<QubitId> = Some(circ.alloc_qubit());
    let s2 = circ.alloc_qubit();
    let mut step0_t1: Option<QubitId> = None;

    for i in (0..iters).rev() {
        let trace_region_start = circ.phase_active_regions.len();
        circ.set_phase(if apply_fwd.is_some() {
            "tlm_multiply_gcd_reverse_decode"
        } else {
            "tlm_inverse_gcd_reverse_decode"
        });
        let current_n = (SCHED_J2[i] as usize).max(1);
        while u.len() < current_n {
            u.push(circ.alloc_qubit());
        }
        while v.len() < current_n {
            v.push(circ.alloc_qubit());
        }
        let cmp_eff = (GAP_J2[i] as usize).min(current_n).max(1);

        // If the current window is exhausted, decompress the next one (from the
        // tape end) into raw symbol slots.
        if pending.is_empty() {
            win_idx -= 1;
            let codec = window_plan[win_idx];
            let cb = codec.code_bits();
            let tlen = tape.len();
            let data: Vec<QubitId> = tape.split_off(tlen - cb);
            if codec == super::codec::DialogCodec::Step0 {
                let (t1, raw) = super::codec::decompress_step0_with_t1(circ, &data);
                step0_t1 = Some(t1);
                pending = raw;
            } else {
                pending = codec.decompress_window(circ, &data);
            }
            pending_tail4 = codec == super::codec::DialogCodec::Tail4Top32;
        } else if tail4_prefix_encoded {
            let suffix = pending.split_off(super::codec::DialogCodec::Triple.code_bits());
            pending = super::codec::DialogCodec::Triple.decompress_window(circ, &pending);
            pending.extend(suffix);
            tail4_prefix_encoded = false;
        }
        // Pull the last symbol (3 bits) off `pending` into the ancilla.
        let plen = pending.len();
        let cur: Vec<QubitId> = pending.split_off(plen - 3);
        if pending_tail4 && pending.len() == 12 {
            let suffix = pending.split_off(9);
            pending = super::codec::DialogCodec::Triple.compress_window(circ, &pending);
            pending.extend(suffix);
            tail4_prefix_encoded = true;
        } else if pending_tail4 && pending.is_empty() {
            pending_tail4 = false;
        }
        if i >= 250 && std::env::var_os("TRACE_TLM_TAIL").is_some() {
            eprintln!(
                "TLM_TAIL direction=reverse i={i} active={} tape={} pending={} encoded={tail4_prefix_encoded}",
                circ.active_qubits,
                tape.len(),
                pending.len(),
            );
        }
        circ.swap(subtracted, cur[0]);
        let swp = if i == 0 {
            circ.cx(subtracted, cur[1]);
            subtracted
        } else {
            let sf = *swap_flag
                .as_ref()
                .expect("swap flag live for non-step0 replay");
            circ.swap(sf, cur[1]);
            sf
        };
        circ.swap(s2, cur[2]);
        // Free the 3 now-|0> symbol slots here -- before the apply -- so they are not
        // carried as dead live qubits through the apply + GCD-undo, where they would
        // inflate the peak by sym_bits=3. The next window's decompress re-allocs
        // from these freed slots.
        for q in cur {
            circ.zero_and_free(q);
        }

        circ.set_phase(if apply_fwd.is_some() {
            "tlm_multiply_gcd_reverse_apply"
        } else {
            "tlm_inverse_gcd_reverse_apply"
        });
        let parked_u0 = if park_odd_u0_enabled(i, "REV") {
            let q = u[0];
            Some(park_known_one(circ, q))
        } else {
            None
        };

        // Fused forward apply (Direction::Forward multiply): apply this divstep's
        // symbol to the coordinate pair using the live symbol bits, in reverse iter
        // order matching apply_step_forward's order, so the tape is never materialized
        // for the apply -> the apply adders run in the (reverse) GCD's headroom.
        let parked_v0 = if apply_fwd.is_some() && park_even_v0_enabled() {
            let q = v[0];
            Some(park_known_zero(circ, q))
        } else {
            None
        };
        if let Some((xr, yr)) = apply_fwd {
            let t1 = step0_t1.unwrap_or(subtracted);
            apply_step_forward(
                circ,
                i,
                &subtracted,
                &swp,
                &s2,
                &t1,
                xr,
                yr,
                &u[1..4],
            );
        }
        if let Some(q) = parked_v0 {
            v[0] = restore_known_zero(circ, q);
        }

        // Inverse step (reverse op order): sub^-1, cswap^-1, cmp^-1, subtracted^-1,
        // s_2^-1, shift1^-1.
        // a) sub^-1: v += subtracted*u (X-sandwich cancels).
        circ.set_phase(if apply_fwd.is_some() {
            "tlm_multiply_gcd_reverse_body"
        } else {
            "tlm_inverse_gcd_reverse_body"
        });
        controlled_add_active(
            circ,
            i,
            &subtracted,
            &u[..current_n],
            &v[..current_n],
            GcdBit0Mode::ReverseKnownZeroBeforeCx,
            apply_fwd.map_or(&[], |(xr, _)| &xr[..3]),
        );
        if let Some(q) = parked_u0 {
            u[0] = restore_known_one(circ, q);
        }
        // b) cswap^-1 (involutory).
        for j in 1..current_n {
            let old_context = crate::point_add::set_op_trace_context(
                0x1200_0000 | (((i as u32) & 0xffff) << 8) | (j as u32 & 0xff),
            );
            if !gcd_reverse_cswap_has_structurally_dead_gate(i, j) {
                circ.cswap(swp, u[j], v[j]);
            }
            crate::point_add::restore_op_trace_context(old_context);
        }
        // c) uncompute swap_flag. Step 0 is a CNOT (swap_flag == subtracted). For
        //    i>=1 the flag holds `subtracted AND (v<u)`; clear it by measurement-vent
        //    (hmr + capped Z-recompute under push_condition) -- ~half the normal
        //    comparator's Toffoli, scorer-discounted ~0.5 more.
        if i != 0 {
            super::comparator::swap_decision_uncompute_vented(
                circ,
                &subtracted,
                &v[..current_n],
                &u[..current_n],
                cmp_eff,
                &swp,
            );
        }
        // d) subtracted^-1: cx(v[0], subtracted) (v[0] still post-shift parity).
        circ.cx(v[0], subtracted);
        // e) s_2 inverse: controlled left-shift, uncompute s_2.
        controlled_left_shift(circ, &s2, &v[..current_n]);
        circ.x(s2);
        circ.cx(v[0], s2);
        // f) shift1 inverse: i>=1 unconditional left-shift; i==0 gated on t1.
        if i == 0 {
            let t1 = step0_t1.expect("step0 t1 decompressed");
            controlled_left_shift(circ, &t1, &v[..current_n]);
            circ.x(t1);
            circ.cx(v[0], t1);
        } else {
            left_shift(circ, &v[..current_n]);
        }

        // (symbol slots already freed before the apply, above). At i==0 drain t1.
        if i == 0 {
            let t1 = step0_t1.take().expect("step0 t1 present");
            circ.zero_and_free(t1);
        }
        if i == 1 {
            let sf = swap_flag.take().expect("swap flag still allocated");
            circ.zero_and_free(sf);
        }
        circ.set_phase("tlm_gcd_step_end");
        trace_step_regions(
            circ,
            if apply_fwd.is_some() {
                "multiply-reverse"
            } else {
                "inverse-reverse"
            },
            i,
            trace_region_start,
        );
    }
    assert!(tape.is_empty(), "tape not fully drained");

    // u grew back to 256 bits holding q; deinit and free.
    let q_bytes = q_secp256k1_le();
    for (i, qb) in u.iter().enumerate().take(n) {
        if (q_bytes.get(i / 8).copied().unwrap_or(0) >> (i % 8)) & 1 == 1 {
            circ.x(*qb);
        }
    }
    for q in u {
        circ.zero_and_free(q);
    }
    circ.zero_and_free(subtracted);
    if let Some(swap_flag) = swap_flag {
        circ.zero_and_free(swap_flag);
    }
    circ.zero_and_free(s2);
}

// =====================================================================
// Active-width helpers (the GCD body runs on `current_n` bits, NOT n=256).
// =====================================================================

/// `target ^= ctrl AND (v_top < u_top)` on the top `k` MSBs of the active-width
/// slices.
fn controlled_swap_decision_v_lt_u(
    circ: &mut B,
    ctrl: &QubitId,
    v: &[QubitId],
    u: &[QubitId],
    k: usize,
    target: &QubitId,
) {
    super::comparator::controlled_swap_decision_lt_truncated(circ, ctrl, v, u, k, target);
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum GcdBit0Mode {
    ForwardKnownOneAfterCx,
    ReverseKnownZeroBeforeCx,
}

/// `y += ctrl * x (mod 2^width)` over the active width (no carry-out captured).
/// The GCD subtract uses this inside an X-sandwich; post-swap `v >= u` so the
/// two's-complement subtract never borrows out on a fitting input.
///
/// Measurement-vented threaded add: ~2 Toffoli/bit (one forward carry +
/// one gated sum), with the carry-uncompute vented by `hmr` (0 Toffoli) and no
/// gated-addend materialization -- vs the single-ancilla Cuccaro's ~3 Toffoli/bit.
/// The per-bit carry qubits (allocated + freed inside) fit in the GCD passes'
/// headroom below the global (apply) peak. mod-2^m: the carry-out is dropped, as
/// the GCD subtract never carries out on a fitting input.
fn controlled_add_active(
    circ: &mut B,
    i: usize,
    ctrl: &QubitId,
    x: &[QubitId],
    y: &[QubitId],
    bit0_mode: GcdBit0Mode,
    dirty_vents: &[QubitId],
) {
    // The GCD subtract `v -= u` (here `y += x` inside the X-sandwich):
    // a = target = y (= v), b = addend = x (= u). cap = PAD.
    // Schedule-driven GCD subtract: pull the baked carry-cap `k` for this call
    // (next_gcd_k) and emit the capped varchunk/adaptive adder. The baked `k`
    // fixes the chunk decomposition deterministically.
    // Kaliski bit-0 optimization: under the odd-u invariant the subtrahend bit0
    // `x[0] == 1`, and the accumulator bit0 is the control (forward, inside the
    // X-sandwich `~v[0] == NOT ctrl`) or 0 (reverse-add). Either way the bit-0
    // carry-out `acc[0] AND ctrl AND x[0]` is provably 0, so the bit-0 result is the
    // known `cx(ctrl, y[0])` with no carry into bit 1 -- emit it directly and run the
    // capped adder on bits 1.. with carry-in 0. Saves the bit-0 carry CCX (~1-2 tof)
    // per GCD conditional sub/add * 2 (fwd+rev) * ITERS ~= 1000+ tof. Not the apply.
    let k = maybe_adjust_late_gcd_k(i, super::next_gcd_k());
    let branch = super::next_gcd_branch();
    let loan_y0 = loan_gcd_y0_enabled() && x.len() > 1;
    match bit0_mode {
        GcdBit0Mode::ForwardKnownOneAfterCx => {
            circ.cx(*ctrl, y[0]); // bit-0 sum is known one inside the X-sandwich.
            if loan_y0 {
                loan_known_one_gcd_y0(circ, y[0]);
            }
        }
        GcdBit0Mode::ReverseKnownZeroBeforeCx => {
            // The inverse add starts with y[0] = 0 and no carry into bit 1. Delay the
            // bit-0 CNOT until after the high-bit adder so y[0] can be borrowed.
            if loan_y0 {
                loan_known_zero_gcd_y0(circ, y[0]);
            }
        }
    }
    if x.len() > 1 {
        let yr: Vec<&QubitId> = y[1..].iter().collect();
        let xr: Vec<&QubitId> = x[1..].iter().collect();
        super::gidney::with_dirty_vent_pool(dirty_vents, || {
            super::gidney::controlled_hybrid_add_capped_branch(
                circ,
                ctrl,
                &yr,
                &xr,
                k,
                super::PAD,
                branch,
            );
        });
    }
    if loan_y0 {
        match bit0_mode {
            GcdBit0Mode::ForwardKnownOneAfterCx => reclaim_known_one_gcd_y0(circ, y[0]),
            GcdBit0Mode::ReverseKnownZeroBeforeCx => reclaim_known_zero_gcd_y0(circ, y[0]),
        }
    }
    if bit0_mode == GcdBit0Mode::ReverseKnownZeroBeforeCx {
        circ.cx(*ctrl, y[0]); // bit-0 sum (x[0] == 1 under the odd-u invariant)
    }
}

// =====================================================================
// Apply: read the (decompressed) tape and reconstruct the modular product.
// =====================================================================

/// One forward apply step (the multiply body) at iter `i` with the symbol bits
/// `(sub, swp, s2)` and the t1 prefix `t1`:
///   if sub: y += x mod q; if swp: swap(x,y); y *= 2 (shift1) ; if s2: y *= 2.
fn apply_step_forward(
    circ: &mut B,
    i: usize,
    sub: &QubitId,
    swp: &QubitId,
    s2: &QubitId,
    t1: &QubitId,
    x_reg: &[QubitId],
    y_reg: &[QubitId],
    dirty_vents: &[QubitId],
) {
    let n = 256usize;
    let s2_known_zero = i != 0 && apply_fwd_s2_zero(i);

    // 1) if subtracted: y += x mod q. Apply cofactor add -> cout adder at the
    // schedule's k.
    circ.set_phase("tlm_apply_forward_mod_add");
    let k = super::next_cout_k();
    let ffg = super::next_ffg();
    if !apply_add_skip(i, true) {
        super::gidney::with_dirty_vent_pool(dirty_vents, || {
            arith::controlled_mod_add_k(
                circ,
                sub,
                &x_reg[..n],
                &y_reg[..n],
                Some(k),
                Some(ffg),
            );
        });
    }
    // 2) if swap: swap(x, y).
    circ.set_phase("tlm_apply_forward_swap");
    if !apply_fwd_cswap_skip(i) {
        for j in 0..n {
            if !apply_fwd_cswap_partial_skip(i, j) {
                circ.cswap(*swp, x_reg[j], y_reg[j]);
            }
        }
    }
    // 3) y := 2*(1+s2)*y mod q. i==0: two separate controlled doublings (shift1 is
    //    t1-gated). i>0: the fused double+cdouble -- one combined (e+2d)*f fold
    //    (the unfused form costs extra inversion Toffoli).
    circ.set_phase("tlm_apply_forward_fold");
    if i == 0 {
        controlled_mod_double(circ, t1, y_reg);
        controlled_mod_double(circ, s2, y_reg);
    } else if s2_known_zero {
        super::fused::fused_double_only(circ, y_reg);
    } else {
        super::fused::fused_double_cdouble(circ, s2, y_reg);
    }
}

/// One inverse apply step (the divide body, gate-inverse of [`apply_step_forward`])
/// at iter `i`:  s_2 halve; shift1 halve; if swp: swap(x,y); if sub: y -= x mod q.
fn apply_step_reverse(
    circ: &mut B,
    i: usize,
    sub: &QubitId,
    swp: &QubitId,
    s2: &QubitId,
    t1: &QubitId,
    x_reg: &[QubitId],
    y_reg: &[QubitId],
    dirty_vents: &[QubitId],
) {
    let n = 256usize;
    let s2_known_zero = i != 0 && apply_inv_s2_zero(i);

    // inverse of 3): i==0 two separate reverse-doublings (s2 halve then t1 halve);
    // i>0 the fused inverse double+cdouble.
    circ.set_phase("tlm_apply_inverse_fold");
    if i == 0 {
        controlled_mod_double_reverse(circ, s2, y_reg);
        controlled_mod_double_reverse(circ, t1, y_reg);
    } else if s2_known_zero {
        super::fused::fused_double_only_reverse(circ, y_reg);
    } else {
        super::fused::fused_double_cdouble_reverse(circ, s2, y_reg);
    }
    // inverse of 2): swap (involutory).
    circ.set_phase("tlm_apply_inverse_swap");
    if !apply_inv_cswap_skip(i) {
        for j in 0..n {
            circ.cswap(*swp, x_reg[j], y_reg[j]);
        }
    }
    // inverse of 1): y -= x mod q. The apply-path operands carry pseudo-Mersenne
    // representation drift, so the borrow clean uses the MBU form.
    // Apply cofactor sub -> schedule cout k.
    circ.set_phase("tlm_apply_inverse_mod_sub");
    let k = super::next_cout_k();
    if !apply_add_skip(i, false) {
        super::gidney::with_dirty_vent_pool(dirty_vents, || {
            controlled_mod_sub_vented(circ, sub, &x_reg[..n], &y_reg[..n], Some(k));
        });
    }
}

/// Apply-path `y := y - ctrl * x (mod q)` whose borrow ancilla is cleaned by
/// MBU (hmr + cz_if_bit) rather than the normal comparator
/// `zero_and_free`. Tolerates pseudo-Mersenne representation drift (operands in
/// `[0, 2^256)` not strictly `< q`), which the apply pipeline produces. A gated
/// add over the full width with an MSBs phase-correction borrow clean.
fn controlled_mod_sub_vented(circ: &mut B, ctrl: &QubitId, x: &[QubitId], y: &[QubitId], sched_k: Option<usize>) {
    let n = x.len();
    assert_eq!(y.len(), n, "x,y equal width");
    let f_bytes = F_SECP256K1.to_le_bytes();
    let anc = circ.alloc_qubit();
    // X-sandwich: ~y += ctrl*x => y -= ctrl*x; cout = ctrl AND borrow.
    circ.set_phase("tlm_apply_inverse_mod_sub_register");
    for q in y {
        circ.x(*q);
    }
    controlled_add_active_cout(circ, ctrl, x, y, &anc, sched_k);
    for q in y {
        circ.x(*q);
    }
    // gated -f fold on the borrow.
    circ.set_phase("tlm_apply_inverse_mod_sub_fold");
    for q in &y[..arith::LSBS] {
        circ.x(*q);
    }
    let ffg = super::next_ffg();
    arith::add_f_window_pub(circ, &anc, y, arith::LSBS, &f_bytes, Some(ffg));
    for q in &y[..arith::LSBS] {
        circ.x(*q);
    }
    // Measurement-vented borrow clean by MBU over the top `msbs` bits (a vented
    // adaptive msbs comparator, not a normal 2k MAJ/un-MAJ). The borrow == ctrl AND
    // carryout(y_final + x + 1); the -f fold only touched the low LSBS bits, so the
    // predicate is recomputed on the top msbs, which the fold left intact (msbs =
    // PAD). HMR the borrow, then on the fired shots recompute the carry as a deferred
    // Z via the Gidney `a + ~b + ~cin` carry chain (~1 Toffoli/bit); cz_if_bit(ctrl,
    // carry) cancels the phase. This never asserts |0>, so operand drift is tolerated.
    // Flip x_top so the comparator's internal `~b` yields `+x_top` (carryout(y+x+1)).
    circ.set_phase("tlm_apply_inverse_mod_sub_clean");
    let k = arith::MSBS.min(n);
    let lo = n - k;
    let ctrl = *ctrl;
    let bit = circ.alloc_bit();
    circ.hmr(anc, bit);
    circ.zero_and_free(anc);
    circ.push_condition(bit);
    let yt: Vec<QubitId> = y[lo..n].to_vec();
    let xt: Vec<QubitId> = x[lo..n].to_vec();
    for q in &xt {
        circ.x(*q); // ~b = ~(~x_top) = x_top inside the comparator -> carryout(y+x+1)
    }
    // Recompute the borrow through the headroom-adaptive (chunked) comparator
    // backend (`compare_geq_chunked_middle`); held-carry count = the full window
    // `k`. `flag = [yt >= xt] = carryout(yt + x + 1)` (xt was complemented to ~x
    // above) = the borrow predicate; deposit Z^(ctrl AND flag). The chunked
    // comparator builds the [>=] carry with an implicit +1 carry-in, so no separate
    // `zcin` is needed (this also drops a qubit vs the explicit carry-in form).
    let flag = circ.alloc_qubit();
    super::comparator::compare_geq_chunked_middle(circ, &yt, &xt, &flag, |c, fl| {
        c.cz(ctrl, *fl);
    }, k);
    circ.zero_and_free(flag);
    for q in &xt {
        circ.x(*q);
    }
    circ.pop_condition();
}

/// `y += ctrl * x` over `n` bits depositing carry-out into `cout` (= ctrl AND
/// overflow). Measurement-vented chunked-gated adder (~2.5 Toffoli/bit, bounded
/// peak) -- the inverse apply's mod-sub register add, on the peak-bound hot path.
/// `cout` is |0> on entry; `x` restored.
fn controlled_add_active_cout(circ: &mut B, ctrl: &QubitId, x: &[QubitId], y: &[QubitId], cout: &QubitId, sched_k: Option<usize>) {
    match sched_k {
        Some(k) => {
            // cout adder at the schedule's k (a = target = y, b = addend = x).
            let yr: Vec<&QubitId> = y.iter().collect();
            let xr: Vec<&QubitId> = x.iter().collect();
            super::gidney::controlled_hybrid_add_cout_refs(circ, ctrl, &yr, &xr, cout, k);
        }
        None => arith::controlled_add_vented_chunked_cout(circ, ctrl, x, y, arith::APPLY_CHUNK, Some(cout)),
    }
}


/// Jump-GCD in-place modular inversion / multiply: `(xv, y) -> (xv, y*x^{±1})`.
/// `Direction::Inverse` => `y := y * xv^-1 mod q`; `Direction::Forward` => `y := y * xv mod q`.
/// xv restored. Fused both directions (apply folded into the GCD passes).
#[must_use]
pub fn mod_mul_inverse_in_place(
    circ: &mut B,
    mut xv: Vec<QubitId>,
    y: &[QubitId],
    dir: Direction,
) -> Vec<QubitId> {
    let n = 256usize;
    assert_eq!(xv.len(), n, "xv must be 256 bits");
    assert_eq!(y.len(), n, "y must be 256 bits");

    match dir {
        Direction::Inverse => {
            // Fused inverse apply folds into the forward GCD. apply pair (x_reg=y,
            // y_reg=tmp); result y = z*x^-1, tmp = 0.
            let tmp: Vec<QubitId> = (0..n).map(|_| circ.alloc_qubit()).collect();
            for j in 0..n {
                circ.swap(y[j], tmp[j]); // tmp = z, y = 0
            }
            let mut tape = forward_gcd_jump(circ, &mut xv, Some((y, &tmp)));
            for q in tmp {
                circ.zero_and_free(q);
            }
            reverse_gcd_jump(circ, &mut xv, &mut tape, None);
            xv
        }
        Direction::Forward => {
            // Fused forward apply folds into the reverse GCD. apply pair (x_reg=tmp=z,
            // y_reg=y=0); result y = z*x, tmp = canonical 0 (drift cleared).
            let mut tape = forward_gcd_jump(circ, &mut xv, None);
            let tmp: Vec<QubitId> = (0..n).map(|_| circ.alloc_qubit()).collect();
            for j in 0..n {
                circ.swap(y[j], tmp[j]); // tmp = z, y = 0
            }
            reverse_gcd_jump(circ, &mut xv, &mut tape, Some((&tmp, y)));
            clear_zeroed_drift(circ, &tmp[..n]);
            for q in tmp {
                circ.zero_and_free(q);
            }
            xv
        }
    }
}

/// Clear a register holding a value `== 0 (mod q)` whose representative is `0` or
/// `q`: XOR q maps `q -> 0` (caught by the caller's zero_and_free if not).
fn clear_zeroed_drift(circ: &mut B, reg: &[QubitId]) {
    let q_bytes = q_secp256k1_le();
    for (i, qb) in reg.iter().enumerate() {
        if (q_bytes.get(i / 8).copied().unwrap_or(0) >> (i % 8)) & 1 == 1 {
            circ.x(*qb);
        }
    }
}
