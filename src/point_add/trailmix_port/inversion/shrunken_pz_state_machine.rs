//! Reversible unpacked PZ inversion as a bit-by-bit pipelined state machine
//! (design reference: `scripts/kaliski_test.py` `pz_big_step`). This supersedes
//! the full-division `shrunken_pz_primitives` module, whose coarser granularity
//! needed a fat quotient pad and did not handle large termination quotients.
//!
//! Per iteration (fixed count ~= sum of quotient bitlengths), gated on the state
//! flags so termination is intrinsic (no separate counter):
//!   DIVISION substep:  s = bitlen(A)-bitlen(B); align B<<s; if A>=B { A-=B;
//!                      `q_div` ^= 1<<s }; restore B>>s. A<B => `div_active=0`.
//!   MULTIPLY substep (pipelined): s = `ctz(q_mul)`; clear it; a += b<<s; restore.
//!                      `q_mul==0` => swap a,b; flip parity; `mul_active=0`.
//!   TRANSITION: q_div->q_mul; swap A,B; divide builds the NEXT quotient while
//!               the multiply drains the PREVIOUS. q pads are TINY (one quotient).
//! All shifts are `controlled_cyclic_rotate` (rotate-in-place, fixed width).
//! Up front: normalize x -> min(x, P-x) (sgn); final a corrected by parity ^ sgn.

#![allow(dead_code)]

use crate::point_add::trailmix_port::circuit::{Circuit, QReg};
use crate::point_add::trailmix_port::inversion::shrunken_pz_primitives::{
    borrow_compare_gated_refs, borrow_compare_gated_refs_with_carry, borrow_compare_refs,
    borrow_compare_refs_with_carry,
};

fn env_usize(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(default)
}

fn trailmix_srot_width() -> usize {
    // The generated schedule's shift bounds need six bits on valid samples.
    // Keep an env override for experiments.
    env_usize("TRAILMIX_SROT_W", 6).max(1)
}

fn q954_srot_counter7_requested() -> bool {
    std::env::var("LOWQ_Q954_SROT_COUNTER7").ok().as_deref() == Some("1")
}

fn trailmix_logical_srot_width() -> usize {
    trailmix_srot_width() + usize::from(q954_srot_counter7_requested())
}

fn trailmix_counter_width() -> usize {
    if std::env::var("TRAILMIX_NO_COUNTER").ok().as_deref() == Some("1") {
        0
    } else {
        env_usize("TRAILMIX_COUNTER_W", 10)
    }
}

fn trailmix_q_width(wq: usize) -> usize {
    let w = wq.max(1);
    std::env::var("TRAILMIX_Q_CAP")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .map_or(w, |cap| w.min(cap.max(1)))
}

/// Per-step quotient width with SELECTIVE peak-targeting.
///
/// The global qubit peak at a `shrunken_pz` step is
///   2*max(wa,wb) + 2*max(wca,wcb) + q_width + FIXED.
/// A blunt global `TRAILMIX_Q_CAP` clamps q on ALL ~490 steps (most have
/// universal q in 23..38), but only the peak-binding step(s) need a smaller q
/// to lower the global peak. Clamping the rest just manufactures classical
/// misses (overflowed quotients) without helping the peak.
///
/// `TRAILMIX_Q_TARGET=T` instead gives each step a budget so that its working
/// width never exceeds T: `q <= T - 2*max(wa,wb) - 2*max(wca,wcb)`. Steps whose
/// other registers are small keep their full natural q (no miss); only the
/// wide-carry peak step(s) get q trimmed, and only by the minimum needed.
/// Falls back to `trailmix_q_width` (global cap) when `TRAILMIX_Q_TARGET` unset.
/// Cap the shared A/B register width (both A and B are resized to max(wa,wb)).
/// `TRAILMIX_AB_CAP` trims it on the steps where it would otherwise bind the peak.
fn trailmix_ab_width(wab: usize) -> usize {
    let w = wab.max(1);
    std::env::var("TRAILMIX_AB_CAP")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .map_or(w, |c| w.min(c.max(1)))
}

/// Cap the shared ca/cb cofactor register width (both resized to max(wca,wcb)).
/// `TRAILMIX_CACB_CAP` trims the dominant 2*245 carry pair at the peak step.
fn trailmix_cacb_width(wcacb: usize) -> usize {
    let w = wcacb.max(1);
    std::env::var("TRAILMIX_CACB_CAP")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .map_or(w, |c| w.min(c.max(1)))
}

/// Fuse the immutable input-sign bit into the EEA parity bit and reclaim the
/// released persistent wire. If `s` is the
/// input sign and `p` is the original EEA parity, the fused state is `s XOR p`.
/// The slope-correction control is therefore its negation, and reverse EEA
/// restores `s XOR 1`, from which `s` is recovered and uncomputed exactly.
fn sign_parity_q_reuse_enabled() -> bool {
    if std::env::var("TRAILMIX_SIGN_PARITY_Q_REUSE")
        .ok()
        .as_deref()
        != Some("1")
    {
        return false;
    }
    assert!(
        matches!(
            std::env::var("TRAILMIX_Q_TARGET").ok().as_deref(),
            Some("683" | "684")
        ),
        "TRAILMIX_SIGN_PARITY_Q_REUSE is sealed to Q_TARGET=683/684"
    );
    true
}

fn trailmix_q_width_step(wq: usize, wa: usize, wb: usize, wca: usize, wcb: usize) -> usize {
    let natural = wq.max(1);
    let target = std::env::var("TRAILMIX_Q_TARGET")
        .ok()
        .and_then(|s| s.parse::<usize>().ok());
    let Some(target) = target else {
        return trailmix_q_width(wq);
    };
    // q budget is computed from the (possibly capped) A/B and ca/cb widths so the
    // working width 2*ab + 2*cacb + q meets `target` consistently with the resizes.
    let other = 2 * trailmix_ab_width(wa.max(wb)) + 2 * trailmix_cacb_width(wca.max(wcb));
    // Q_TARGET=684 retains the audited quotient widths. The two reclaimed
    // persistent lanes are physical savings and are not added back to q.
    let budget = target.saturating_sub(other).max(1);
    // Still honor a global Q_CAP if both are set (take the tighter bound).
    let capped = natural.min(budget);
    std::env::var("TRAILMIX_Q_CAP")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .map_or(capped, |cap| capped.min(cap.max(1)))
        .max(1)
}

fn compute_active(c: &mut Circuit, counter: &[QReg], candidates: &[&QReg]) -> QReg {
    let active = c.alloc_qreg("active");
    if counter.is_empty() {
        c.x(&active);
    } else if lowq_q959_selective_borrow_enabled() {
        toggle_zero_dirty(c, counter, &active, candidates, &[&active]);
    } else {
        or_is_zero(c, counter, &active);
    }
    active
}

fn uncompute_active(c: &mut Circuit, counter: &[QReg], active: &QReg, candidates: &[&QReg]) {
    if counter.is_empty() {
        c.x(active);
    } else if lowq_q959_selective_borrow_enabled() {
        toggle_zero_dirty(c, counter, active, candidates, &[active]);
    } else {
        or_is_zero(c, counter, active);
    }
}

/// `p + 1` (secp256k1 base field prime) as 33 LE bytes.
fn p_plus_1_bytes() -> Vec<u8> {
    vec![
        0x30, 0xfc, 0xff, 0xff, 0xfe, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x00,
    ]
}

/// Controlled field-negate `a := (p - a) mod p` IFF `g` (a in [0,p), 257-bit).
/// Self-inverse. `~a + (p+1) ≡ p - a (mod 2^257)`; canonical for a in [1,p).
/// (Relocated from `kaliski_spooky::unpacked` so `shrunken_pz` has no spooky-Kaliski dep.)
pub fn controlled_field_neg(c: &mut Circuit, g: &QReg, a: &[QReg]) {
    use crate::point_add::trailmix_port::arith::const_add::controlled_add_const;
    for q in a {
        c.cx(g, q);
    }
    controlled_add_const(c, g, a, &p_plus_1_bytes());
}

/// Canonical controlled field negation. Unlike `controlled_field_neg`, this
/// leaves zero at zero when the control is set instead of producing the
/// congruent but noncanonical representative `p`.
fn controlled_field_neg_canonical(c: &mut Circuit, g: &QReg, a: &[QReg]) {
    assert_eq!(a.len(), 257, "canonical field negation requires 257 lanes");
    let nonzero = c.alloc_qreg("field-neg.nonzero");
    let apply = c.alloc_qreg("field-neg.apply");
    or_nonzero(c, a, &nonzero);
    c.ccx(g, &nonzero, &apply);
    controlled_field_neg(c, &apply, a);
    c.ccx(g, &nonzero, &apply);
    or_nonzero(c, a, &nonzero);
    c.zero_and_free(apply);
    c.zero_and_free(nonzero);
}

/// `s += bitlen(a) - bitlen(b)` (clz diff), bound by `bound`. After alignment in
/// the division substep, s is the shift to apply. Inverse: swap a,b.
/// LEAN `bit_length`: `s += bitlen(src)` (or `-=` if dec), via a reversible
/// prefix-AND ladder + gray-code deposit -- ~2n ccx (ladder build+unbuild) with
/// NO per-row position-equality. Supersedes the first-hit scan (~38 tof/row from
/// the per-row `toggle_on_cursor_eq_const` uncompute of `is_hit`).
///
/// Construction (MSB-first running flag `f_i` = "no 1 bit strictly above i"):
///   - prefix-AND ladder over ~src (X-bracketed) gives every `f_i` as a ladder
///     qubit, fully reversibly (fwd builds, rev unbuilds).
///   - deposit pos (init = n) ^= (i ^ (i+1)) gated on `f_i`, for i = n-1..0. The
///     gray differences telescope: pos collapses to the MSB index p (= bitlen-1).
///   - s += (pos + 1)  [bitlen]; then uncompute pos (re-run deposit) + ladder.
///
/// PRE: src nonzero (EEA gcd / nonzero quotient pad). For src==0 this returns
/// bitlen=1 (pos stays 0, +1); callers must not pass an all-zero src.
/// _middle core. Builds the prefix-AND ladder over ~src, deposits the MSB index
/// (= bitlen-1) into the caller's `pos` register (PRE: pos = |n>) in the FORWARD
/// sweep, runs `body` (which sees pos = MSB index), then unbuilds.
///
/// `body` returns whether the deposit should be UNDONE on the reverse sweep:
///   - `false` (DEFAULT, 3n): pos is KEPT at the MSB index -- the caller owns it
///     and must clear it later (e.g. via the SM's reverse). One consume = 3n.
///   - `true` (4n): the deposit is re-run on the reverse, returning pos to |n>.
///     Use when pos is a throwaway temp whose value was folded elsewhere in body.
///
/// The gray-code deposit is pure XOR (CX gated on a single flag materialized from
/// the prefix-AND with one ccx, then HMR-freed) -- so each consume is 1 toffoli
/// per position. Prefix build+unbuild = 2n; consume = n/sweep.
fn bit_length_lean_middle(
    circ: &mut Circuit,
    src: &[&QReg],
    pos: &[QReg],
    body: impl FnOnce(&mut Circuit) -> bool,
) {
    use crate::point_add::trailmix_port::arith::khattar_gidney::{kg_prefix_ancilla_count, KgPrefixAnd};
    let n = src.len();
    if n == 0 {
        body(circ);
        return;
    }
    // ~src (X-bracket); the prefix-AND reads the complemented bits.
    for q in src {
        circ.x(q);
    }
    // q = ~src MSB-first: q[j] = ~src[n-1-j]. The log*-ancilla KG streaming
    // prefix-AND gives, at layer i, AND(ctrls) = AND(q[0..i]) = "no 1 in top i
    // positions" = f_k ("no 1 strictly above k") for k = n-1-i. ctrls is 1-2 qubits
    // (KG conditionally-clean form), so the deposit is the KG prefix-controlled-X
    // consumer directly: CX (1 ctrl, zero toffoli) or CCX (2 ctrls) per gray bit --
    // NO mcx materialize. Total ~3n-4n (2n prefix compute + n-2n consume).
    let qbits: Vec<&QReg> = src.iter().rev().copied().collect();
    let nanc = kg_prefix_ancilla_count(n);
    let anc_owned = circ.alloc_qreg_bits("bll.kganc", nanc);
    let anc: Vec<&QReg> = anc_owned.iter().collect();
    let flag = circ.alloc_qreg("bll.flag");

    // Deposit at layer i (position k = n-1-i): gray-XOR (k ^ (k+1)) into pos gated
    // on f_k = AND(ctrls). For a 2-qubit ctrls, materialize f_k onto `flag` with ONE
    // ccx, CX the gray bits (free), then free `flag` via clear_and (HMR + cz_if_bit,
    // ZERO toffoli) -- so the consume is 1 toffoli/position. For <=1 ctrl the gray
    // bits are a direct CX/X (zero toffoli). pos starts at |n>; the gray differences
    // telescope it to the MSB index p. Self-inverse, so reverse undoes pos to |n>.
    fn deposit_step(
        circ: &mut Circuit,
        i: usize,
        ctrls: &[&QReg],
        pos: &[QReg],
        flag: &QReg,
        n: usize,
    ) {
        if i >= n {
            return; // i == n is the empty (k = -1) layer
        }
        let k = n - 1 - i;
        let gd = k ^ (k + 1);
        let bits: Vec<usize> = (0..pos.len()).filter(|&b| (gd >> b) & 1 == 1).collect();
        if bits.is_empty() {
            return;
        }
        match ctrls {
            [] => {
                for &b in &bits {
                    circ.x(&pos[b]);
                }
            }
            [c] => {
                for &b in &bits {
                    circ.cx(c, &pos[b]);
                }
            }
            [a, b2] => {
                circ.ccx(a, b2, flag); // flag = f_k (1 toffoli)
                for &b in &bits {
                    circ.cx(flag, &pos[b]); // free
                }
                circ.clear_and(flag, a, b2); // free flag via HMR+CZ (0 toffoli)
            }
            _ => unreachable!("KG prefix ctrls is <=2 qubits"),
        }
    }

    let kg = KgPrefixAnd::new(&qbits, &anc);
    let done = kg.forward(circ, |c, i, ctrls| deposit_step(c, i, ctrls, pos, &flag, n)); // pos -> p
    let clean = body(circ);
    if clean {
        // 4n: re-run the deposit on the reverse, returning pos to |n>.
        done.reverse(circ, |c, i, ctrls| deposit_step(c, i, ctrls, pos, &flag, n));
    } else {
        // 3n: unbuild the prefix only; pos stays at the MSB index (caller-owned).
        done.reverse(circ, |_, _, _| {});
    }
    circ.zero_and_free(flag);
    drop(anc);
    for q in anc_owned {
        circ.zero_and_free(q);
    }
    for q in src {
        circ.x(q);
    }
}

/// `s += bitlen(src)` (or `-=` if dec). Built from [`bit_length_lean_middle`]:
/// pos = MSB index in the middle, then `s ±= (pos + 1)`. With `dec` this clears a
/// register `s` that already holds `bitlen(src)` (the "same method" both ways).
fn bit_length_lean(circ: &mut Circuit, src: &[&QReg], s: &[QReg], dec: bool) {
    let n = src.len();
    if n == 0 {
        return;
    }
    let pbl = circ.push_section("p.bitlen");
    // pos holds transient gray values up to (n-1)^n < 2n; reuse s's width (equal-
    // width so the Cuccaro add s += pos is clean).
    let pos_w = s.len();
    debug_assert!(
        (n as u64) <= (1u64 << (pos_w - 1)),
        "bit_length_lean: s width {pos_w} too small for n={n}"
    );
    let pos = circ.alloc_qreg_bits("bll.pos", pos_w);
    xor_const(circ, &pos, n); // pos = n  (PRE for the middle)
    bit_length_lean_middle(circ, src, &pos, |circ| {
        // pos = MSB index = bitlen-1; s ±= (pos + 1).
        if dec {
            for q in s {
                circ.x(q);
            }
        }
        let pref: Vec<&QReg> = pos.iter().collect();
        let sref: Vec<&QReg> = s.iter().collect();
        add_refs(circ, &sref, &pref); // s += pos
        let one = circ.alloc_qreg("bll.one");
        circ.x(&one);
        ctrl_inc(circ, &one, s); // s += 1  (bitlen = p + 1)
        circ.x(&one);
        circ.zero_and_free(one);
        if dec {
            for q in s {
                circ.x(q);
            }
        }
        true // pos is a throwaway temp -> clean on reverse (4n)
    });
    xor_const(circ, &pos, n); // pos back to |0>
    for q in pos {
        circ.zero_and_free(q);
    }
    circ.pop_section(&pbl);
}

fn lowq_clz_diff_const_fold_enabled() -> bool {
    if std::env::var("LOWQ_CLZ_DIFF_CONST_FOLD").ok().as_deref() != Some("1") {
        return false;
    }
    let target = std::env::var("TRAILMIX_Q_TARGET")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .expect("LOWQ_CLZ_DIFF_CONST_FOLD requires an integer TRAILMIX_Q_TARGET");
    assert!(
        matches!(target, 683 | 684 | 685),
        "LOWQ_CLZ_DIFF_CONST_FOLD repair audit permits Q_TARGET 683/684/685"
    );
    true
}

fn lowq_hybrid_clz_enabled() -> bool {
    if std::env::var("LOWQ_HYBRID_CLZ").ok().as_deref() != Some("1") {
        return false;
    }
    assert_eq!(
        trailmix_logical_srot_width(),
        5,
        "LOWQ_HYBRID_CLZ requires the five-bit shift register"
    );
    assert_eq!(
        env_usize("TRAILMIX_THIN_CLZ_WINDOW", 0),
        78,
        "LOWQ_HYBRID_CLZ is sealed to the audited 78-bit windows"
    );
    assert!(
        matches!(env_usize("TRAILMIX_Q_TARGET", 0), 683 | 684 | 685),
        "LOWQ_HYBRID_CLZ repair audit permits Q_TARGET 683/684/685"
    );
    true
}

fn lowq_exact_ctz_enabled() -> bool {
    if std::env::var("LOWQ_EXACT_CTZ").ok().as_deref() != Some("1") {
        return false;
    }
    assert_eq!(
        trailmix_logical_srot_width(),
        5,
        "LOWQ_EXACT_CTZ requires the five-bit shift register"
    );
    assert!(
        matches!(env_usize("TRAILMIX_Q_TARGET", 0), 683 | 684 | 685),
        "LOWQ_EXACT_CTZ repair audit permits Q_TARGET 683/684/685"
    );
    true
}

fn lowq_hybrid_clz_kg_mcx_enabled() -> bool {
    std::env::var("LOWQ_HYBRID_CLZ_KG_MCX").ok().as_deref() == Some("1")
}

fn lowq_hybrid_clz_prefix_parity_enabled() -> bool {
    std::env::var("LOWQ_HYBRID_CLZ_PREFIX_PARITY").ok().as_deref() == Some("1")
}

fn lowq_hybrid_clz_noalloc_add_enabled() -> bool {
    std::env::var("LOWQ_HYBRID_CLZ_NOALLOC_ADD").ok().as_deref() == Some("1")
}

fn lowq_q957_target683_enabled() -> bool {
    std::env::var("LOWQ_Q957_TARGET683").ok().as_deref() == Some("1")
}

fn lowq_q959_selective_borrow_enabled() -> bool {
    if std::env::var("LOWQ_Q959_SELECTIVE_BORROW").ok().as_deref() != Some("1") {
        return false;
    }
    assert_eq!(
        trailmix_logical_srot_width(),
        5,
        "LOWQ_Q959_SELECTIVE_BORROW requires five shift lanes"
    );
    assert_eq!(
        env_usize("TRAILMIX_THIN_CLZ_WINDOW", 0),
        78,
        "LOWQ_Q959_SELECTIVE_BORROW is sealed to the 78-bit schedule"
    );
    let q_target = env_usize("TRAILMIX_Q_TARGET", 0);
    assert!(
        q_target == 684 || (q_target == 683 && lowq_q957_target683_enabled()),
        "LOWQ_Q959_SELECTIVE_BORROW requires Q_TARGET=684 or the Q957 target683 route"
    );
    assert_eq!(
        std::env::var("TRAILMIX_SIGN_PARITY_Q_REUSE").ok().as_deref(),
        Some("1"),
        "LOWQ_Q959_SELECTIVE_BORROW requires sign/parity fusion"
    );
    assert_eq!(
        std::env::var("LOWQ_EXACT_CTZ").ok().as_deref(),
        Some("1"),
        "LOWQ_Q959_SELECTIVE_BORROW requires exact in-place CTZ"
    );
    true
}

fn lowq_q958_gated_compare_enabled() -> bool {
    if std::env::var("LOWQ_Q958_GATED_COMPARE").ok().as_deref() != Some("1") {
        return false;
    }
    assert!(
        lowq_q959_selective_borrow_enabled(),
        "LOWQ_Q958_GATED_COMPARE requires the sealed selective-borrow route"
    );
    true
}

fn lowq_q956_off_borrow_enabled() -> bool {
    if std::env::var("LOWQ_Q956_OFF_BORROW").ok().as_deref() != Some("1") {
        return false;
    }
    assert!(
        lowq_q958_gated_compare_enabled(),
        "LOWQ_Q956_OFF_BORROW requires the sealed Q958 gated-comparator route"
    );
    assert!(
        lowq_q957_target683_enabled(),
        "LOWQ_Q956_OFF_BORROW requires the Q957 target683 route"
    );
    assert_eq!(
        env_usize("TRAILMIX_Q_TARGET", 0),
        683,
        "LOWQ_Q956_OFF_BORROW is sealed to Q_TARGET=683"
    );
    assert_eq!(
        env_usize("TRAILMIX_Q_CAP", 0),
        99,
        "LOWQ_Q956_OFF_BORROW is sealed to Q_CAP=99"
    );
    assert_eq!(
        env_usize("TRAILMIX_COUNTER_W", 0),
        8,
        "LOWQ_Q956_OFF_BORROW is sealed to the eight-bit counter"
    );
    assert_eq!(
        trailmix_logical_srot_width(),
        5,
        "LOWQ_Q956_OFF_BORROW requires five logical arithmetic shift lanes"
    );
    assert!(
        std::env::var_os("TRAILMIX_PASSENGER_TOP_Q_REUSE").is_none(),
        "LOWQ_Q956_OFF_BORROW forbids passenger-top reuse"
    );
    true
}

fn lowq_q955_off_canonical_enabled() -> bool {
    if std::env::var("LOWQ_Q955_OFF_CANONICAL").ok().as_deref() != Some("1") {
        return false;
    }
    assert!(
        lowq_q956_off_borrow_enabled(),
        "LOWQ_Q955_OFF_CANONICAL requires the sealed Q956 off-borrow route"
    );
    assert_eq!(
        env_usize("TRAILMIX_Q_TARGET", 0),
        683,
        "LOWQ_Q955_OFF_CANONICAL is sealed to Q_TARGET=683"
    );
    assert_eq!(
        env_usize("TRAILMIX_Q_CAP", 0),
        99,
        "LOWQ_Q955_OFF_CANONICAL preserves the Q_CAP=99 support widths"
    );
    assert_eq!(
        env_usize("TRAILMIX_COUNTER_W", 0),
        8,
        "LOWQ_Q955_OFF_CANONICAL requires the eight-bit counter"
    );
    assert!(
        std::env::var_os("TRAILMIX_PASSENGER_TOP_Q_REUSE").is_none(),
        "LOWQ_Q955_OFF_CANONICAL forbids passenger-top reuse"
    );
    assert!(
        std::env::var_os("TRAILMIX_Q_MODEL_GUARD").is_none(),
        "LOWQ_Q955_OFF_CANONICAL forbids TRAILMIX_Q_MODEL_GUARD"
    );
    true
}

const Q954_FIRST_TERMINAL_ROW: usize = 371;
const Q954_LAST_CTZ_BIT4_ROW: usize = 477;
const Q954_LAST_RAW_BIT4_BARREL_ROW: usize = 495;
const Q954_MAX_PRE_BODY_COUNTER: usize = 124;
const Q954_MAX_FINAL_COUNTER: usize = 159;
const Q954_SCHEDULE_FINGERPRINT: u64 = 0xf128_4a16_5e9c_235d;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Q954ScheduleCertificate {
    pub rows: usize,
    pub first_terminal_capable_row: usize,
    pub last_ctz_bit4_row: usize,
    pub last_raw_bit4_barrel_row: usize,
    pub max_pre_body_counter: usize,
    pub max_final_counter: usize,
    pub fingerprint: u64,
}

/// Bind the counter[7] alias proof to the exact generated 530-row schedule.
/// Terminal timing is authoritative support metadata; the raw barrel cutoff is
/// also re-derived from the checked shift-bound rows. Any schedule-array change
/// must be reviewed and issued a new fingerprint before this route can run.
#[doc(hidden)]
pub fn q954_srot_counter7_schedule_certificate() -> Q954ScheduleCertificate {
    use crate::point_add::trailmix_port::inversion::shrunken_pz_schedule::{
        SHRUNKEN_PZ_A, SHRUNKEN_PZ_A_LO, SHRUNKEN_PZ_B, SHRUNKEN_PZ_B_LO,
        SHRUNKEN_PZ_CA, SHRUNKEN_PZ_CA_LO, SHRUNKEN_PZ_CB, SHRUNKEN_PZ_CB_LO,
        SHRUNKEN_PZ_NSTEPS, SHRUNKEN_PZ_Q, SHRUNKEN_PZ_Q_LO, SHRUNKEN_PZ_S2,
        SHRUNKEN_PZ_SDIV,
    };

    static CERTIFICATE: std::sync::OnceLock<Q954ScheduleCertificate> =
        std::sync::OnceLock::new();
    *CERTIFICATE.get_or_init(|| {
        fn mix(hash: &mut u64, value: u16) {
            for byte in value.to_le_bytes() {
                *hash ^= u64::from(byte);
                *hash = hash.wrapping_mul(0x100_0000_01b3);
            }
        }

        assert_eq!(SHRUNKEN_PZ_NSTEPS, 530, "Q954 schedule row-count drift");
        let mut fingerprint = 0xcbf2_9ce4_8422_2325u64;
        mix(&mut fingerprint, SHRUNKEN_PZ_NSTEPS as u16);
        for rows in [
            &SHRUNKEN_PZ_A,
            &SHRUNKEN_PZ_B,
            &SHRUNKEN_PZ_CA,
            &SHRUNKEN_PZ_CB,
            &SHRUNKEN_PZ_Q,
            &SHRUNKEN_PZ_A_LO,
            &SHRUNKEN_PZ_B_LO,
            &SHRUNKEN_PZ_CA_LO,
            &SHRUNKEN_PZ_CB_LO,
            &SHRUNKEN_PZ_Q_LO,
            &SHRUNKEN_PZ_SDIV,
            &SHRUNKEN_PZ_S2,
        ] {
            for &value in rows {
                mix(&mut fingerprint, value);
            }
        }
        assert_eq!(
            fingerprint, Q954_SCHEDULE_FINGERPRINT,
            "Q954 schedule fingerprint drift; reject counter[7] alias"
        );

        let last_raw_bit4_barrel_row = (0..SHRUNKEN_PZ_NSTEPS)
            .rev()
            .find(|&row| SHRUNKEN_PZ_SDIV[row] >= 16 || SHRUNKEN_PZ_S2[row] >= 16)
            .expect("Q954 schedule must exercise barrel bit4");
        assert_eq!(
            last_raw_bit4_barrel_row, Q954_LAST_RAW_BIT4_BARREL_ROW,
            "Q954 raw bit4 barrel cutoff drift"
        );
        assert_eq!(
            Q954_LAST_RAW_BIT4_BARREL_ROW - Q954_FIRST_TERMINAL_ROW,
            Q954_MAX_PRE_BODY_COUNTER,
            "Q954 pre-body counter bound drift"
        );
        assert_eq!(
            SHRUNKEN_PZ_NSTEPS - Q954_FIRST_TERMINAL_ROW,
            Q954_MAX_FINAL_COUNTER,
            "Q954 final counter bound drift"
        );
        assert!(
            Q954_MAX_PRE_BODY_COUNTER < (1 << 7),
            "Q954 counter[7] is not clean at the final bit4 barrel use"
        );
        assert!(
            Q954_MAX_FINAL_COUNTER < (1 << 8),
            "Q954 terminal counter exceeds its eight-bit register"
        );

        Q954ScheduleCertificate {
            rows: SHRUNKEN_PZ_NSTEPS,
            first_terminal_capable_row: Q954_FIRST_TERMINAL_ROW,
            last_ctz_bit4_row: Q954_LAST_CTZ_BIT4_ROW,
            last_raw_bit4_barrel_row,
            max_pre_body_counter: Q954_MAX_PRE_BODY_COUNTER,
            max_final_counter: Q954_MAX_FINAL_COUNTER,
            fingerprint,
        }
    })
}

fn lowq_q954_srot_counter7_enabled() -> bool {
    if !q954_srot_counter7_requested() {
        return false;
    }
    assert!(
        lowq_q955_off_canonical_enabled(),
        "LOWQ_Q954_SROT_COUNTER7 requires the composed Q955 canonical route"
    );
    assert_eq!(
        trailmix_srot_width(),
        4,
        "LOWQ_Q954_SROT_COUNTER7 allocates exactly four owned shift lanes"
    );
    assert_eq!(
        trailmix_counter_width(),
        8,
        "LOWQ_Q954_SROT_COUNTER7 requires counter[7]"
    );
    let certificate = q954_srot_counter7_schedule_certificate();
    assert_eq!(certificate.rows, 530);
    true
}

fn q954_ctz_width(row: usize) -> usize {
    if lowq_q954_srot_counter7_enabled() && row > Q954_LAST_CTZ_BIT4_ROW {
        4
    } else {
        5
    }
}

fn with_arithmetic_srot_view<'a, R>(
    owned: &'a [QReg],
    counter: &'a [QReg],
    body: impl FnOnce(&[&'a QReg]) -> R,
) -> R {
    if lowq_q954_srot_counter7_enabled() {
        assert_eq!(owned.len(), 4, "Q954 owned shift-lane count drift");
        assert_eq!(counter.len(), 8, "Q954 counter width drift");
        // This borrowed view exists only inside an already-gated arithmetic
        // body. The body is an exact cleanup block, so counter[7] is restored
        // before gate-holder or done logic evaluates the full counter again.
        let split = [
            &owned[0],
            &owned[1],
            &owned[2],
            &owned[3],
            &counter[7],
        ];
        body(&split)
    } else {
        let refs: Vec<&QReg> = owned.iter().collect();
        body(&refs)
    }
}

/// Undo the route-specific representation used to reconstruct a field product.
/// Q955 keeps every Horner state canonical; earlier routes retain the original
/// rfold representation and its matched cleanup.
pub(crate) fn shrunken_pz_product_undo(
    c: &mut Circuit,
    result: &[QReg],
    a: &[QReg],
    b: &[QReg],
) {
    if lowq_q955_off_canonical_enabled() {
        crate::point_add::trailmix_port::arith::rfold_mbu::mod_mul_canonical_mbu_undo(
            c, result, a, b,
        );
    } else {
        crate::point_add::trailmix_port::arith::rfold_mbu::mod_mul_rfold_mbu_undo(
            c, result, a, b,
        );
    }
}

fn assert_q956_off_alias(
    off: &QReg,
    counter: &[QReg],
    s_rot: &[QReg],
) {
    assert!(!counter.is_empty(), "Q956 off borrow requires a counter lane");
    assert!(
        std::ptr::eq(off, &counter[0]),
        "Q956 off must alias counter[0] exactly"
    );
    assert!(s_rot.len() >= 3, "Q956 boundary predicates require s_rot[0..3]");
    assert!(
        s_rot.iter().all(|lane| !std::ptr::eq(off, lane))
            && counter[1..].iter().all(|lane| !std::ptr::eq(off, lane)),
        "Q956 off alias overlaps a protected state lane"
    );
}

/// One controlled fixed-distance shift layer. The forward direction is a
/// logical left shift on the promised branch because its top `distance` lanes
/// are zero. Reversing the pair order is the exact inverse for arbitrary data.
fn controlled_fixed_shift(
    circ: &mut Circuit,
    reg: &[QReg],
    control: &QReg,
    distance: usize,
    forward: bool,
) {
    if distance == 0 || distance >= reg.len() {
        return;
    }
    if forward {
        for hi in (distance..reg.len()).rev() {
            circ.cswap(control, &reg[hi], &reg[hi - distance]);
        }
    } else {
        for hi in distance..reg.len() {
            circ.cswap(control, &reg[hi], &reg[hi - distance]);
        }
    }
}

/// Toggle `out` iff the highest `prefix` lanes of `src` are all zero. The
/// peer register supplies restored dirty lenders and is unchanged.
fn toggle_zero_prefix_dirty(
    circ: &mut Circuit,
    src: &[QReg],
    prefix: usize,
    out: &QReg,
    peer: &[QReg],
    clean_scratch: &[QReg],
) {
    use crate::point_add::trailmix_port::arith::khattar_gidney::{
        kg_prefix_ancilla_count, xor_and_of_khattar_gidney_refs_with_anc,
    };
    use crate::point_add::trailmix_port::arith::mcx::mcx_dirty_ladder;

    assert!(prefix > 0 && prefix < src.len());
    let controls_owned = &src[src.len() - prefix..];
    for q in controls_owned {
        circ.x(q);
    }
    let controls: Vec<&QReg> = controls_owned.iter().collect();
    let clean_refs: Vec<&QReg> = clean_scratch.iter().collect();
    if lowq_hybrid_clz_kg_mcx_enabled()
        && prefix >= 6
        && clean_refs.len() >= kg_prefix_ancilla_count(prefix)
    {
        xor_and_of_khattar_gidney_refs_with_anc(circ, &controls, out, &clean_refs);
    } else {
        let dirty: Vec<&QReg> = peer.iter().take(prefix.saturating_sub(2)).collect();
        assert_eq!(
            dirty.len(),
            prefix.saturating_sub(2),
            "LOWQ_HYBRID_CLZ peer lender shortage"
        );
        mcx_dirty_ladder(circ, &controls, out, &dirty);
    }
    for q in controls_owned.iter().rev() {
        circ.x(q);
    }
}

/// Toggle `out` iff `active` is set and the lowest `prefix` lanes of `src` are
/// all zero. Lenders may contain arbitrary data and are restored exactly.
fn toggle_active_zero_low_dirty(
    circ: &mut Circuit,
    src: &[QReg],
    prefix: usize,
    active: &QReg,
    out: &QReg,
    lenders: &[&QReg],
) {
    use crate::point_add::trailmix_port::arith::mcx::mcx_dirty_ladder;

    assert!(prefix > 0 && prefix < src.len());
    let controls_owned = &src[..prefix];
    for q in controls_owned {
        circ.x(q);
    }
    let mut controls: Vec<&QReg> = Vec::with_capacity(prefix + 1);
    controls.push(active);
    controls.extend(controls_owned.iter());
    let need = controls.len().saturating_sub(2);
    assert!(
        lenders.len() >= need,
        "LOWQ_EXACT_CTZ lender shortage: need={need} have={}",
        lenders.len()
    );
    mcx_dirty_ladder(circ, &controls, out, &lenders[..need]);
    for q in controls_owned.iter().rev() {
        circ.x(q);
    }
}

/// Compute `transcript = clz(src)` and normalize `src` to an MSB-one word.
/// Each branch bit controls one power-of-two shift and is retained until the
/// inverse restores `src`, so the map is bijective on the full basis space.
fn binary_clz_compute(
    circ: &mut Circuit,
    src: &[QReg],
    peer: &[QReg],
    transcript: &[QReg],
) {
    assert!(!src.is_empty() && src.len() <= (1usize << transcript.len()));
    for bit in (0..transcript.len()).rev() {
        let distance = 1usize << bit;
        if distance >= src.len() {
            continue;
        }
        toggle_zero_prefix_dirty(circ, src, distance, &transcript[bit], peer, &transcript[..bit]);
        controlled_fixed_shift(circ, src, &transcript[bit], distance, true);
    }
}

fn binary_clz_uncompute(
    circ: &mut Circuit,
    src: &[QReg],
    peer: &[QReg],
    transcript: &[QReg],
) {
    for bit in 0..transcript.len() {
        let distance = 1usize << bit;
        if distance >= src.len() {
            continue;
        }
        controlled_fixed_shift(circ, src, &transcript[bit], distance, false);
        toggle_zero_prefix_dirty(circ, src, distance, &transcript[bit], peer, &transcript[..bit]);
    }
}

fn toggle_prefix_controlled_by_active(
    circ: &mut Circuit,
    ctrls: &[&QReg],
    active: &QReg,
    out: &QReg,
    flag: &QReg,
) {
    match ctrls {
        [] => circ.cx(active, out),
        [c] => circ.ccx(active, c, out),
        [a, b] => {
            circ.ccx(a, b, flag);
            circ.ccx(active, flag, out);
            circ.clear_and(flag, a, b);
        }
        _ => panic!(
            "toggle_prefix_controlled_by_active: expected <=2 KG controls, got {}",
            ctrls.len()
        ),
    }
}

fn toggle_clz_parity_prefix_stream(
    circ: &mut Circuit,
    src: &[QReg],
    active: &QReg,
    out: &QReg,
    scratch: &[QReg],
) -> bool {
    use crate::point_add::trailmix_port::arith::khattar_gidney::{
        kg_prefix_ancilla_count, KgPrefixAnd,
    };

    if src.len() <= 1 {
        return true;
    }
    let qbits: Vec<&QReg> = src.iter().rev().take(src.len() - 1).collect();
    let nanc = kg_prefix_ancilla_count(qbits.len());
    if scratch.len() < nanc + 1 {
        return false;
    }
    let anc: Vec<&QReg> = scratch[..nanc].iter().collect();
    let flag = &scratch[nanc];

    for &q in &qbits {
        circ.x(q);
    }
    KgPrefixAnd::new(&qbits, &anc)
        .forward(circ, |_, _, _| {})
        .reverse(circ, |c, i, ctrls| {
            if i > 0 {
                toggle_prefix_controlled_by_active(c, ctrls, active, out, flag);
            }
        });
    for &q in qbits.iter().rev() {
        circ.x(q);
    }
    true
}

/// PRE: `s=0`. Deposit `active*ctz(q)` directly into `s`, using `s` itself as
/// the branch transcript. The final left-shift sweep restores multi-hot q while
/// intentionally retaining s.
fn exact_multihot_ctz_deposit(
    circ: &mut Circuit,
    q: &[QReg],
    s: &[&QReg],
    active: &QReg,
    lenders: &[&QReg],
) {
    assert!(!s.is_empty() && s.len() <= 5, "LOWQ exact CTZ output width");
    let prev = circ.push_section("p.hctz.deposit");
    for bit in (0..s.len()).rev() {
        let distance = 1usize << bit;
        if distance >= q.len() {
            continue;
        }
        toggle_active_zero_low_dirty(circ, q, distance, active, s[bit], lenders);
        controlled_fixed_shift(circ, q, s[bit], distance, false);
    }
    for bit in 0..s.len() {
        let distance = 1usize << bit;
        if distance < q.len() {
            controlled_fixed_shift(circ, q, s[bit], distance, true);
        }
    }
    circ.pop_section(&prev);
}

/// Exact gate inverse of `exact_multihot_ctz_deposit`.
/// PRE: `s=active*ctz(q)`. Restores q after the temporary normalization and
/// clears s to zero.
fn exact_multihot_ctz_erase(
    circ: &mut Circuit,
    q: &[QReg],
    s: &[&QReg],
    active: &QReg,
    lenders: &[&QReg],
) {
    assert!(!s.is_empty() && s.len() <= 5, "LOWQ exact CTZ output width");
    let prev = circ.push_section("p.hctz.erase");
    for bit in (0..s.len()).rev() {
        let distance = 1usize << bit;
        if distance < q.len() {
            controlled_fixed_shift(circ, q, s[bit], distance, false);
        }
    }
    for bit in 0..s.len() {
        let distance = 1usize << bit;
        if distance >= q.len() {
            continue;
        }
        controlled_fixed_shift(circ, q, s[bit], distance, true);
        toggle_active_zero_low_dirty(circ, q, distance, active, s[bit], lenders);
    }
    circ.pop_section(&prev);
}

fn collect_dirty_lenders<'a>(
    candidates: impl IntoIterator<Item = &'a QReg>,
    controls: &[&QReg],
    action: &[&QReg],
) -> Vec<&'a QReg> {
    let mut out: Vec<&'a QReg> = Vec::new();
    for q in candidates {
        if controls.iter().any(|c| std::ptr::eq(*c, q))
            || action.iter().any(|a| std::ptr::eq(*a, q))
            || out.iter().any(|d| std::ptr::eq(*d, q))
        {
            continue;
        }
        out.push(q);
    }
    out
}

/// Exact multi-control toggle using arbitrary dirty lenders. Every lender is
/// restored before return; no clean quantum lane is allocated.
fn dirty_controlled_x(
    circ: &mut Circuit,
    controls: &[&QReg],
    target: &QReg,
    candidates: &[&QReg],
    action: &[&QReg],
) {
    use crate::point_add::trailmix_port::arith::mcx::mcx_dirty_ladder;

    let dirty = collect_dirty_lenders(candidates.iter().copied(), controls, action);
    let need = controls.len().saturating_sub(2);
    assert!(
        dirty.len() >= need,
        "Q959 selective-borrow lender shortage: controls={} need={} have={}",
        controls.len(),
        need,
        dirty.len()
    );
    mcx_dirty_ladder(circ, controls, target, &dirty[..need]);
}

fn dirty_controlled_inc_suffix(
    circ: &mut Circuit,
    selector: &[&QReg],
    target: &[&QReg],
    lo: usize,
    subtract: bool,
    candidates: &[&QReg],
) {
    let action = target.to_vec();
    for i in (lo + 1..target.len()).rev() {
        let lower = target[lo..i].to_vec();
        if subtract {
            for q in &lower {
                circ.x(q);
            }
        }
        let mut controls = selector.to_vec();
        controls.extend(lower.iter().copied());
        dirty_controlled_x(circ, &controls, target[i], candidates, &action);
        if subtract {
            for q in lower.iter().rev() {
                circ.x(q);
            }
        }
    }
    dirty_controlled_x(circ, selector, target[lo], candidates, &action);
}

fn dirty_controlled_add_const(
    circ: &mut Circuit,
    selector: &[&QReg],
    target: &[&QReg],
    value: usize,
    subtract: bool,
    candidates: &[&QReg],
) {
    let mask = (1usize << target.len()) - 1;
    let value = value & mask;
    for bit in 0..target.len() {
        if (value >> bit) & 1 == 1 {
            dirty_controlled_inc_suffix(circ, selector, target, bit, subtract, candidates);
        }
    }
}

/// Add or subtract the promised nonzero source bit length directly into the
/// existing five-bit target. Scanning from high to low, the borrowed selector
/// gate latches exactly once at the source MSB. One unit update per remaining
/// position then deposits `msb - lo + 1`; a final controlled constant update
/// supplies `lo`. The high suffix stays complemented across adjacent selectors
/// instead of being rebuilt for every candidate MSB.
fn direct_bitlen_update(
    circ: &mut Circuit,
    src: &[QReg],
    peer: &[QReg],
    lo: usize,
    target: &[&QReg],
    active: &QReg,
    selector_gate: &QReg,
    subtract: bool,
    extra_lenders: &[&QReg],
) {
    let lo = lo.min(src.len().saturating_sub(1));
    let candidates: Vec<&QReg> = peer
        .iter()
        .chain(src.iter())
        .chain(extra_lenders.iter().copied())
        .collect();
    let mut action = target.to_vec();
    action.push(selector_gate);
    for k in (lo..src.len()).rev() {
        let mut selector = Vec::with_capacity(src.len() - k + 1);
        selector.push(active);
        selector.push(&src[k]);
        selector.extend(src[k + 1..].iter());
        // The selector is one-hot over k. Once it fires, every lower-k selector
        // is false because the complemented suffix contains the true MSB.
        dirty_controlled_x(circ, &selector, selector_gate, &candidates, &action);
        if lowq_q956_off_borrow_enabled() {
            // `selector_gate` aliases counter[0]. It is guaranteed clean only
            // when active, so every read must carry the active predicate too.
            dirty_controlled_inc_suffix(
                circ,
                &[active, selector_gate],
                target,
                0,
                subtract,
                &candidates,
            );
        } else {
            dirty_controlled_inc_suffix(
                circ,
                &[selector_gate],
                target,
                0,
                subtract,
                &candidates,
            );
        }
        if k > lo {
            circ.x(&src[k]);
        }
    }
    for q in src[lo + 1..].iter().rev() {
        circ.x(q);
    }

    // On the promised support, active implies that the selected window is
    // nonzero, so exactly one MSB selector fired and selector_gate == active.
    circ.cx(active, selector_gate);
    if lo != 0 {
        dirty_controlled_add_const(circ, &[active], target, lo, subtract, &candidates);
    }
}

fn direct_bitlen_diff_update(
    circ: &mut Circuit,
    a: &[QReg],
    b: &[QReg],
    lo_a: usize,
    lo_b: usize,
    target: &[&QReg],
    active: &QReg,
    selector_gate: &QReg,
    subtract_diff: bool,
    extra_lenders: &[&QReg],
) {
    let prev = circ.push_section("p.dbitlen");
    direct_bitlen_update(
        circ,
        a,
        b,
        lo_a,
        target,
        active,
        selector_gate,
        subtract_diff,
        extra_lenders,
    );
    direct_bitlen_update(
        circ,
        b,
        a,
        lo_b,
        target,
        active,
        selector_gate,
        !subtract_diff,
        extra_lenders,
    );
    circ.pop_section(&prev);
}

fn direct_bitlen_diff_parity(
    circ: &mut Circuit,
    a: &[QReg],
    b: &[QReg],
    lo_a: usize,
    lo_b: usize,
    out: &QReg,
    active: &QReg,
    extra_lenders: &[&QReg],
) {
    let prev = circ.push_section("p.dbitlen.parity");
    let action = [out];
    for (src, peer, lo) in [(a, b, lo_a), (b, a, lo_b)] {
        let lo = lo.min(src.len().saturating_sub(1));
        let candidates: Vec<&QReg> = peer
            .iter()
            .chain(src.iter())
            .chain(extra_lenders.iter().copied())
            .collect();
        for k in (lo..src.len()).rev() {
            if (k + 1) & 1 == 1 {
                let mut selector = Vec::with_capacity(src.len() - k + 1);
                selector.push(active);
                selector.push(&src[k]);
                selector.extend(src[k + 1..].iter());
                dirty_controlled_x(circ, &selector, out, &candidates, &action);
            }
            if k > lo {
                circ.x(&src[k]);
            }
        }
        for q in src[lo + 1..].iter().rev() {
            circ.x(q);
        }
    }
    circ.pop_section(&prev);
}

fn hybrid_transcript_width(max_window_len: usize) -> usize {
    let branch_bits = if max_window_len <= 1 {
        0
    } else {
        usize::BITS as usize - (max_window_len - 1).leading_zeros() as usize
    };
    branch_bits.max(5)
}

fn selective_direct_bitlen_needed(
    circ: &mut Circuit,
    a: &[QReg],
    b: &[QReg],
    lo_a: usize,
    lo_b: usize,
) -> bool {
    if !lowq_q959_selective_borrow_enabled() {
        return false;
    }
    circ.flush_pending_frees();
    let aw = a.len().saturating_sub(lo_a.min(a.len().saturating_sub(1)));
    let bw = b.len().saturating_sub(lo_b.min(b.len().saturating_sub(1)));
    let peak_target = if lowq_q954_srot_counter7_enabled() {
        954
    } else if lowq_q956_off_borrow_enabled() {
        956
    } else if lowq_q957_target683_enabled() {
        957
    } else if lowq_q958_gated_compare_enabled() {
        958
    } else {
        959
    };
    circ.b.active_qubits as usize + hybrid_transcript_width(aw.max(bw)) > peak_target
}

/// Deposit `active*(bitlen(a)-bitlen(b))` into the existing five-bit shift
/// register. Equal full register widths imply
/// `bitlen(a)-bitlen(b) = clz(b)-clz(a)` even when the audited low windows
/// differ. A single seven-bit transcript is reused sequentially.
fn hybrid_bitlen_diff_update(
    circ: &mut Circuit,
    a: &[QReg],
    b: &[QReg],
    lo_a: usize,
    lo_b: usize,
    target: &[&QReg],
    active: &QReg,
    subtract_diff: bool,
) {
    use crate::point_add::trailmix_port::inversion::shrunken_pz_primitives::{
        ctrl_add, ctrl_add_dirty_lenders, ctrl_sub, ctrl_sub_dirty_lenders,
    };

    assert_eq!(a.len(), b.len(), "LOWQ_HYBRID_CLZ requires equal full widths");
    assert_eq!(target.len(), 5, "LOWQ_HYBRID_CLZ target width");
    let prev = circ.push_section("p.hclz");
    let a_window = &a[lo_a.min(a.len() - 1)..];
    let b_window = &b[lo_b.min(b.len() - 1)..];
    let transcript = circ.alloc_qreg_bits(
        "hybrid.clz",
        hybrid_transcript_width(a_window.len().max(b_window.len())),
    );
    let target_refs = target.to_vec();
    let low_refs: Vec<&QReg> = transcript[..target.len()].iter().collect();
    let noalloc_add = lowq_hybrid_clz_noalloc_add_enabled();

    binary_clz_compute(circ, a_window, b, &transcript);
    if subtract_diff {
        if noalloc_add {
            ctrl_add_dirty_lenders(circ, active, &target_refs, &low_refs);
        } else {
            ctrl_add(circ, active, &target_refs, &low_refs);
        }
    } else {
        if noalloc_add {
            ctrl_sub_dirty_lenders(circ, active, &target_refs, &low_refs);
        } else {
            ctrl_sub(circ, active, &target_refs, &low_refs);
        }
    }
    binary_clz_uncompute(circ, a_window, b, &transcript);

    binary_clz_compute(circ, b_window, a, &transcript);
    if subtract_diff {
        if noalloc_add {
            ctrl_sub_dirty_lenders(circ, active, &target_refs, &low_refs);
        } else {
            ctrl_sub(circ, active, &target_refs, &low_refs);
        }
    } else {
        if noalloc_add {
            ctrl_add_dirty_lenders(circ, active, &target_refs, &low_refs);
        } else {
            ctrl_add(circ, active, &target_refs, &low_refs);
        }
    }
    binary_clz_uncompute(circ, b_window, a, &transcript);

    for q in transcript {
        circ.zero_and_free(q);
    }
    circ.pop_section(&prev);
}

fn hybrid_bitlen_diff_parity(
    circ: &mut Circuit,
    a: &[QReg],
    b: &[QReg],
    lo_a: usize,
    lo_b: usize,
    out: &QReg,
    active: &QReg,
) {
    assert_eq!(a.len(), b.len(), "LOWQ_HYBRID_CLZ requires equal full widths");
    let prev = circ.push_section("p.hclz.parity");
    let a_window = &a[lo_a.min(a.len() - 1)..];
    let b_window = &b[lo_b.min(b.len() - 1)..];
    let transcript = circ.alloc_qreg_bits(
        "hybrid.clz",
        hybrid_transcript_width(a_window.len().max(b_window.len())),
    );

    if lowq_hybrid_clz_prefix_parity_enabled()
        && toggle_clz_parity_prefix_stream(circ, a_window, active, out, &transcript)
        && toggle_clz_parity_prefix_stream(circ, b_window, active, out, &transcript)
    {
        // Fast exact parity path: clz(x) mod 2 is the XOR of all non-empty
        // top-zero prefix flags of x. No controlled shifts are needed.
    } else {
        binary_clz_compute(circ, a_window, b, &transcript);
        circ.ccx(active, &transcript[0], out);
        binary_clz_uncompute(circ, a_window, b, &transcript);
        binary_clz_compute(circ, b_window, a, &transcript);
        circ.ccx(active, &transcript[0], out);
        binary_clz_uncompute(circ, b_window, a, &transcript);
    }

    for q in transcript {
        circ.zero_and_free(q);
    }
    circ.pop_section(&prev);
}

/// `_middle` form of the clz-diff compute-USE-uncompute pattern: deposits the two
/// bitlen positions into the internal `pa`/`pb` ancillae, FOLDS the diff
/// d = bitlen(a)-bitlen(b) (windowed) INTO `pa`, runs `body(circ, &pa)` with `pa`
/// holding the diff, then restores `pa` and un-deposits to |0>. No caller-supplied
/// diff register -- `pa` IS the diff, so nothing extra is live at the peak (this is
/// the `shrunken_pz_divide_forward` peak section). `w` sizes pa/pb (must hold the window MSB
/// index and the signed diff). Scans un-nested (one KG ancilla set live at a time).
fn clz_diff_body_middle(
    circ: &mut Circuit,
    a: &[QReg],
    b: &[QReg],
    w: usize,
    lo_a: usize,
    lo_b: usize,
    body: impl FnOnce(&mut Circuit, &[QReg]),
) {
    use crate::point_add::trailmix_port::arith::ripple_add::add_const;
    let pbl = circ.push_section("p.bitlen");
    let aw: Vec<&QReg> = a[lo_a..a.len()].iter().collect();
    let bw: Vec<&QReg> = b[lo_b..b.len()].iter().collect();
    let pa = circ.alloc_qreg_bits("clzm.pa", w);
    let pb = circ.alloc_qreg_bits("clzm.pb", w);
    let add_pa = |circ: &mut Circuit, pa: &[QReg], v: i64| {
        let val = i128::from(v).rem_euclid(1i128 << w) as u128;
        let bytes: Vec<u8> = (0..w.div_ceil(8)).map(|i| (val >> (8 * i)) as u8).collect();
        add_const(circ, pa, &bytes);
    };
    let (na, nb) = (aw.len(), bw.len());
    // UN-NESTED scans: deposit pos_a then pos_b SEQUENTIALLY (one KG ancilla set
    // live at a time, not both nested). `bit_length_lean_middle` with a `|_| false`
    // body deposits pos (na -> MSB index) and leaves it; the pos-telescoping is a
    // fixed XOR-set gated on `src` only (independent of pos's value), hence
    // self-inverse -- the SAME call run again returns pos (MSB index -> na), so it
    // doubles as the un-deposit phase.
    xor_const(circ, &pa, na);
    bit_length_lean_middle(circ, &aw, &pa, |_| false); // pa = pos_a
    xor_const(circ, &pb, nb);
    bit_length_lean_middle(circ, &bw, &pb, |_| false); // pb = pos_b

    let const_fold = lowq_clz_diff_const_fold_enabled();
    if const_fold {
        // Constants commute across the subtract. This is the q980 reduction:
        // one modular constant add instead of two, with no extra live wires.
        {
            let par: Vec<&QReg> = pa.iter().collect();
            let pbr: Vec<&QReg> = pb.iter().collect();
            sub_refs(circ, &par, &pbr);
        }
        add_pa(circ, &pa, lo_a as i64 - lo_b as i64);
    } else {
        {
            let par: Vec<&QReg> = pa.iter().collect();
            let pbr: Vec<&QReg> = pb.iter().collect();
            add_pa(circ, &pa, 1 + lo_a as i64);
            sub_refs(circ, &par, &pbr);
        }
        add_pa(circ, &pa, -(1 + lo_b as i64));
    }

    body(circ, &pa); // USE pa (= diff)

    if const_fold {
        {
            let par: Vec<&QReg> = pa.iter().collect();
            let pbr: Vec<&QReg> = pb.iter().collect();
            add_refs(circ, &par, &pbr);
        }
        add_pa(circ, &pa, lo_b as i64 - lo_a as i64);
    } else {
        add_pa(circ, &pa, 1 + lo_b as i64);
        {
            let par: Vec<&QReg> = pa.iter().collect();
            let pbr: Vec<&QReg> = pb.iter().collect();
            add_refs(circ, &par, &pbr);
        }
        add_pa(circ, &pa, -(1 + lo_a as i64));
    }

    // un-deposit (self-inverse clean=false calls, reverse order).
    bit_length_lean_middle(circ, &bw, &pb, |_| false); // pb -> nb
    xor_const(circ, &pb, nb); // pb -> 0
    bit_length_lean_middle(circ, &aw, &pa, |_| false); // pa -> na
    xor_const(circ, &pa, na); // pa -> 0
    for q in pa {
        circ.zero_and_free(q);
    }
    for q in pb {
        circ.zero_and_free(q);
    }
    circ.pop_section(&pbl);
}

/// Rotate-LEFT `reg` in place by the quantum amount `s` (= reg << s, since the
/// aligned value's bitlen <= reg width so no nonzero bit wraps). Uses the ACYCLIC
/// `barrel_shift_inplace` (exactly `s.len()` layers, no wrap) rather than
/// `controlled_cyclic_rotate` (s.len()+1 full-width layers incl. a spurious
/// offset layer, + cyclic wrap churn): ~1.28x fewer cswaps. The no-wrap
/// precondition (top s bits of reg are |0>) is exactly the existing one.
/// forward=true is `<< s`; forward=false (restore) is `>> s`, Fredkin self-inverse.
fn barrel_shift_refs(circ: &mut Circuit, reg: &[QReg], s: &[&QReg], forward: bool) {
    let n = reg.len();
    if n == 0 || s.is_empty() {
        return;
    }
    let prev = circ.push_section("p.shift");
    let layers: Box<dyn Iterator<Item = usize>> = if forward {
        Box::new(0..s.len())
    } else {
        Box::new((0..s.len()).rev())
    };
    for bit in layers {
        let distance = 1usize << bit;
        if distance >= n {
            continue;
        }
        let pairs: Box<dyn Iterator<Item = usize>> = if forward {
            Box::new((distance..n).rev())
        } else {
            Box::new(distance..n)
        };
        for hi in pairs {
            let lo = hi - distance;
            circ.cx(&reg[lo], &reg[hi]);
            circ.ccx(s[bit], &reg[hi], &reg[lo]);
            circ.cx(&reg[lo], &reg[hi]);
        }
    }
    circ.pop_section(&prev);
}

fn rotate_left(circ: &mut Circuit, reg: &[QReg], s: &[&QReg]) {
    barrel_shift_refs(circ, reg, s, true);
}
fn rotate_right(circ: &mut Circuit, reg: &[QReg], s: &[&QReg]) {
    barrel_shift_refs(circ, reg, s, false);
}

/// Shift by one under `active AND off`. On the Q956 route `off` is a counter
/// lane that may be one on inactive branches, so using it as a lone Fredkin
/// control would corrupt those branches. The three-control swap is emitted
/// directly with restored dirty lenders and allocates no conjunction lane.
fn rotate_one_by_off(
    circ: &mut Circuit,
    reg: &[QReg],
    active: &QReg,
    off: &QReg,
    forward: bool,
    candidates: &[&QReg],
) {
    if !lowq_q956_off_borrow_enabled() {
        let control = [off];
        if forward {
            rotate_left(circ, reg, &control);
        } else {
            rotate_right(circ, reg, &control);
        }
        return;
    }
    if reg.len() < 2 {
        return;
    }

    let prev = circ.push_section("p.shift.off-borrow");
    let action: Vec<&QReg> = reg.iter().collect();
    let pairs: Box<dyn Iterator<Item = usize>> = if forward {
        Box::new((1..reg.len()).rev())
    } else {
        Box::new(1..reg.len())
    };
    for hi in pairs {
        let lo = hi - 1;
        circ.cx(&reg[lo], &reg[hi]);
        let controls = [active, off, &reg[hi]];
        dirty_controlled_x(circ, &controls, &reg[lo], candidates, &action);
        circ.cx(&reg[lo], &reg[hi]);
    }
    circ.pop_section(&prev);
}

/// `q[i] ^= active AND (s == i)` = `q ^= active·(1<<s)` -- the q-demux via KG
/// `unary_iterate_log_star` (~2 ccx/step) instead of a per-bit `eq_const_inplace` loop
/// (~58 tof/bit, ~30x more). active=0 => s masked to 0 => only i=0 gate fires,
/// `ANDed` with active=0 -> no-op. Self-inverse; `s` restored on exit.
fn set_bit_at_s_gated(
    circ: &mut Circuit,
    q_div: &[QReg],
    s: &[&QReg],
    active: &QReg,
    borrowed_gate: &QReg,
    lenders: &[&QReg],
) {
    let n_pad = q_div.len();
    if n_pad == 0 {
        return;
    }
    let prev = circ.push_section("p.demux");
    if lowq_q959_selective_borrow_enabled() {
        let mask_borrowed_reads = lowq_q956_off_borrow_enabled();
        let mut action: Vec<&QReg> = q_div.iter().collect();
        action.push(borrowed_gate);
        for (i, target) in q_div.iter().enumerate() {
            for (bit, q) in s.iter().enumerate() {
                if (i >> bit) & 1 == 0 {
                    circ.x(q);
                }
            }
            let mut controls = Vec::with_capacity(s.len() + 1);
            controls.push(active);
            controls.extend(s.iter().copied());
            dirty_controlled_x(circ, &controls, borrowed_gate, lenders, &action);
            if mask_borrowed_reads {
                circ.ccx(active, borrowed_gate, target);
            } else {
                circ.cx(borrowed_gate, target);
            }
            dirty_controlled_x(circ, &controls, borrowed_gate, lenders, &action);
            for (bit, q) in s.iter().enumerate().rev() {
                if (i >> bit) & 1 == 0 {
                    circ.x(q);
                }
            }
        }
        circ.pop_section(&prev);
        return;
    }

    use crate::point_add::trailmix_port::arith::khattar_gidney::unary_iterate_log_star;
    unary_iterate_log_star(circ, s, n_pad, |c, i, gate| {
        c.ccx(active, gate, &q_div[i]);
    });
    circ.pop_section(&prev);
}

/// Unconditional `a -= b` (mod 2^len) via two's complement (X-bracket + add).
fn sub_refs(circ: &mut Circuit, a: &[&QReg], b: &[&QReg]) {
    use crate::point_add::trailmix_port::inversion::shrunken_pz_primitives::ctrl_sub;
    let one = circ.alloc_qreg("sm.one");
    circ.x(&one);
    ctrl_sub(circ, &one, a, b); // gated on |1> = unconditional
    circ.x(&one);
    circ.zero_and_free(one);
}

/// Controlled decrement `s -= 1` iff `g` (X-bracket + controlled increment).
fn ctrl_dec(circ: &mut Circuit, g: &QReg, s: &[QReg]) {
    use crate::point_add::trailmix_port::arith::khattar_gidney::cinc_khattar_gidney;
    for q in s {
        circ.x(q);
    }
    cinc_khattar_gidney(circ, s, g); // a=s, ctrl=g
    for q in s {
        circ.x(q);
    }
}

/// Controlled increment `s += 1` iff `g`.
fn ctrl_inc(circ: &mut Circuit, g: &QReg, s: &[QReg]) {
    use crate::point_add::trailmix_port::arith::khattar_gidney::cinc_khattar_gidney;
    cinc_khattar_gidney(circ, s, g);
}

fn ctrl_inc_refs(circ: &mut Circuit, g: &QReg, s: &[&QReg]) {
    use crate::point_add::trailmix_port::arith::khattar_gidney::cinc_khattar_gidney_refs;
    cinc_khattar_gidney_refs(circ, s, g);
}

fn ctrl_dec_refs(circ: &mut Circuit, g: &QReg, s: &[&QReg]) {
    for q in s {
        circ.x(q);
    }
    ctrl_inc_refs(circ, g, s);
    for q in s {
        circ.x(q);
    }
}

fn ctrl_inc_by_off(
    circ: &mut Circuit,
    active: &QReg,
    off: &QReg,
    s: &[&QReg],
    candidates: &[&QReg],
) {
    if lowq_q956_off_borrow_enabled() {
        dirty_controlled_inc_suffix(circ, &[active, off], s, 0, false, candidates);
    } else {
        ctrl_inc_refs(circ, off, s);
    }
}

fn ctrl_dec_by_off(
    circ: &mut Circuit,
    active: &QReg,
    off: &QReg,
    s: &[&QReg],
    candidates: &[&QReg],
) {
    if lowq_q956_off_borrow_enabled() {
        dirty_controlled_inc_suffix(circ, &[active, off], s, 0, true, candidates);
    } else {
        ctrl_dec_refs(circ, off, s);
    }
}

/// Unconditional `a += b` (mod 2^len) via a |1>-gated controlled add.
fn add_refs(circ: &mut Circuit, a: &[&QReg], b: &[&QReg]) {
    use crate::point_add::trailmix_port::inversion::shrunken_pz_primitives::ctrl_add;
    let one = circ.alloc_qreg("sm.one_a");
    circ.x(&one);
    ctrl_add(circ, &one, a, b);
    circ.x(&one);
    circ.zero_and_free(one);
}

/// Unpacked PZ state-machine registers. gcd pair (`a_gcd=A`, `b_gcd=B`) shrinks;
/// cofactor pair (ca=|a|, cb=|b|) grows. `q_div/q_mul` are the quotient pads
/// (~one quotient, ~26 bits each): `q_div` is built by the division (`q_div^=1`<<s),
/// swapped to `q_mul`, and DRAINED by the multiply (a += b<<`ctz(q_mul)`, clearing
/// it) -- the pipelined drain is what keeps the quotient record at one-quotient
/// size instead of a full ~256-bit tape. NOT removable (scripts/
/// `pz_fused_nopad_proto.py`: fusing gives the right inverse but s-recovery from
/// the cofactors mismatches ~30%, and an undrained pad accumulates a full tape).
pub struct PzSmRegs {
    pub a_gcd: Vec<QReg>,
    pub b_gcd: Vec<QReg>,
    pub ca: Vec<QReg>,
    pub cb: Vec<QReg>,
    pub q_div: Vec<QReg>,
    pub q_mul: Vec<QReg>,
}

/// Single-qubit state flags + sign. Invariant matches `pz_big_step`.
pub struct PzSmFlags {
    pub div_active: QReg,
    pub mul_active: QReg,
    pub offset: QReg,
    pub parity: QReg,
    pub sgn: QReg,
}

/// Load/unload the classical constant `c` into `reg` via X gates (self-inverse).
fn xor_const(circ: &mut Circuit, reg: &[QReg], c: usize) {
    for (j, q) in reg.iter().enumerate() {
        if (c >> j) & 1 == 1 {
            circ.x(q);
        }
    }
}

/// Magnitude compare `out ^= (a < b)` narrowed to the schedule window
/// `[lo, min(a.len, b.len))`. Used for the ALIGNED offset/o compares where a and
/// b share a bitlen (MSB guaranteed in [lo, hi) by the schedule), so the top bits
/// decide the order; a tie below `lo` (prob ~2^-(hi-lo) per the window width)
/// flips the result -- within the whole-pass tail tolerance. Forward and inverse
/// substeps call this with the same `lo`, so the (possibly-wrong) flag is
/// computed identically both ways and round-trips cleanly. Restores a,b.
/// NOT for the magnitude GATES (`g_mul/g_div)`: there A,B get arbitrarily close at
/// the div<->mul transition, so a deep tie is common, not a 2^-w tail.
fn narrow_lt(circ: &mut Circuit, a: &[QReg], b: &[QReg], out: &QReg, lo: usize) {
    let hi = a.len().min(b.len());
    let lo = lo.min(hi.saturating_sub(1));
    let ar: Vec<&QReg> = a[lo..hi].iter().collect();
    let br: Vec<&QReg> = b[lo..hi].iter().collect();
    borrow_compare_refs(circ, &ar, &br, out);
}

/// Toggle `out` by `active AND (a < b)` without materializing a separate
/// comparison result. The comparator still restores its one clean carry lane,
/// so this saves one peak-live qubit and one complete comparator replay.
fn narrow_lt_controlled(
    circ: &mut Circuit,
    a: &[QReg],
    b: &[QReg],
    out: &QReg,
    active: &QReg,
    lo: usize,
) {
    let hi = a.len().min(b.len());
    let lo = lo.min(hi.saturating_sub(1));
    let ar: Vec<&QReg> = a[lo..hi].iter().collect();
    let br: Vec<&QReg> = b[lo..hi].iter().collect();
    borrow_compare_gated_refs(circ, &ar, &br, active, out);
}

/// WINDOWED division substep: same as `division_substep_act` but the two clz
/// computations scan only the schedule's clz windows (`lo_a`/`lo_b` = window low
/// bounds for A/B) and the B<<s / restore rotates use `rot_bits` shift bits
/// (shift bound) instead of the full `s_rot` width. The offset-clean clz operates
/// on (A, `B_aligned`), both ~bitlen(A), so it reuses the A window (`lo_a`). For
/// in-schedule inputs this is gate-identical to `division_substep_act`.
#[allow(clippy::too_many_arguments)]
pub fn division_substep_windowed(
    circ: &mut Circuit,
    a: &[QReg],
    b: &[QReg],
    q_div: &[QReg],
    s_rot: &[&QReg],
    offset: &QReg,
    active: &QReg,
    extra_lenders: &[&QReg],
    lo_a: usize,
    lo_b: usize,
    rot_bits: usize,
    ctz_bits: usize,
) {
    use crate::point_add::trailmix_port::inversion::shrunken_pz_primitives::ctrl_sub;
    let aref: Vec<&QReg> = a.iter().collect();
    let bref: Vec<&QReg> = b.iter().collect();
    let n_pad = q_div.len();
    let rb = rot_bits.min(s_rot.len());
    let w = s_rot.len();
    let off_lenders: Vec<&QReg> = a
        .iter()
        .chain(b.iter())
        .chain(q_div.iter())
        .chain(s_rot.iter().copied())
        .chain(extra_lenders.iter().copied())
        .collect();

    // diff = bitlen(A)-bitlen(B) (windowed _middle, folded into the clz's own pa);
    // mask s_rot = diff AND active.
    if selective_direct_bitlen_needed(circ, a, b, lo_a, lo_b) {
        direct_bitlen_diff_update(
            circ,
            a,
            b,
            lo_a,
            lo_b,
            s_rot,
            active,
            offset,
            false,
            extra_lenders,
        );
    } else if lowq_hybrid_clz_enabled() {
        hybrid_bitlen_diff_update(circ, a, b, lo_a, lo_b, s_rot, active, false);
    } else {
        clz_diff_body_middle(circ, a, b, w, lo_a, lo_b, |circ, diff| {
            for j in 0..w {
                circ.ccx(active, &diff[j], &s_rot[j]);
            }
        });
    }

    rotate_left(circ, b, &s_rot[0..rb]); // B <<= s if active (bounded rotator)

    // offset = active AND (A < B_aligned) -- narrowed (A,B_aligned share bitlen).
    if lowq_q958_gated_compare_enabled() {
        narrow_lt_controlled(circ, a, b, offset, active, lo_a);
    } else {
        let or = circ.alloc_qreg("dg.offr");
        narrow_lt(circ, a, b, &or, lo_a);
        circ.ccx(active, &or, offset);
        narrow_lt(circ, a, b, &or, lo_a);
        circ.zero_and_free(or);
    }
    rotate_one_by_off(circ, b, active, offset, false, &off_lenders); // B >>= 1 if offset
    ctrl_dec_by_off(circ, active, offset, s_rot, &off_lenders); // s_rot -= 1 if offset => s_eff

    // clean offset via windowed _middle clz on (A, B_aligned) -> A window. The diff
    // lives in the clz's pa (this clz is the shrunken_pz_divide_forward peak section).
    if selective_direct_bitlen_needed(circ, a, b, lo_a, lo_a) {
        direct_bitlen_diff_parity(circ, a, b, lo_a, lo_a, offset, active, extra_lenders);
    } else if lowq_hybrid_clz_enabled() {
        hybrid_bitlen_diff_parity(circ, a, b, lo_a, lo_a, offset, active);
    } else {
        clz_diff_body_middle(circ, a, b, w, lo_a, lo_a, |circ, diff| {
            circ.ccx(active, &diff[0], offset);
        });
    }

    ctrl_sub(circ, active, &aref, &bref); // A -= B_aligned if active

    let demux_lenders: Vec<&QReg> = a
        .iter()
        .chain(b.iter())
        .chain(extra_lenders.iter().copied())
        .collect();
    set_bit_at_s_gated(circ, q_div, s_rot, active, offset, &demux_lenders);

    rotate_right(circ, b, &s_rot[0..rb]); // restore B >>= s_eff (bounded rotator)

    if lowq_exact_ctz_enabled() {
        let lenders: Vec<&QReg> = a
            .iter()
            .chain(b.iter())
            .chain(extra_lenders.iter().copied())
            .collect();
        exact_multihot_ctz_erase(
            circ,
            q_div,
            &s_rot[..ctz_bits.min(s_rot.len())],
            active,
            &lenders,
        );
    } else {
        let t = circ.alloc_qreg_bits("dg.ctz", w);
        xor_const(circ, &t, n_pad);
        let rev: Vec<&QReg> = q_div.iter().rev().collect();
        bit_length_lean(circ, &rev, &t, true);
        let srr = s_rot.to_vec();
        let tr: Vec<&QReg> = t.iter().collect();
        ctrl_sub(circ, active, &srr, &tr);
        bit_length_lean(circ, &rev, &t, false);
        xor_const(circ, &t, n_pad);
        for lane in t {
            circ.zero_and_free(lane);
        }
    }
}

/// Gate-by-gate INVERSE of `division_substep_windowed` (for the backward pass).
/// Reverses the op sequence; the compute-use-uncompute blocks (clz-mask, offset,
/// offset-clean, q-demux) are self-inverse and run as-is; `rotate_left`<->right,
/// ctrl_sub->ctrl_add, ctrl_dec->ctrl_inc flip. Restores A += B<<`s_eff`, clears
/// the `q_div` bit, leaving `A/B/q_div/s/s_rot/offset` as before the forward step.
#[allow(clippy::too_many_arguments)]
pub fn division_substep_windowed_inv(
    circ: &mut Circuit,
    a: &[QReg],
    b: &[QReg],
    q_div: &[QReg],
    s_rot: &[&QReg],
    offset: &QReg,
    active: &QReg,
    extra_lenders: &[&QReg],
    lo_a: usize,
    lo_b: usize,
    rot_bits: usize,
    ctz_bits: usize,
) {
    use crate::point_add::trailmix_port::inversion::shrunken_pz_primitives::ctrl_add;
    let aref: Vec<&QReg> = a.iter().collect();
    let bref: Vec<&QReg> = b.iter().collect();
    let n_pad = q_div.len();
    let rb = rot_bits.min(s_rot.len());
    let w = s_rot.len();
    let off_lenders: Vec<&QReg> = a
        .iter()
        .chain(b.iter())
        .chain(q_div.iter())
        .chain(s_rot.iter().copied())
        .chain(extra_lenders.iter().copied())
        .collect();

    // 12' reconstruct s_rot from the least-significant set quotient bit.
    if lowq_exact_ctz_enabled() {
        let lenders: Vec<&QReg> = a
            .iter()
            .chain(b.iter())
            .chain(extra_lenders.iter().copied())
            .collect();
        exact_multihot_ctz_deposit(
            circ,
            q_div,
            &s_rot[..ctz_bits.min(s_rot.len())],
            active,
            &lenders,
        );
    } else {
        let t = circ.alloc_qreg_bits("dg.ctz", w);
        xor_const(circ, &t, n_pad);
        let rev: Vec<&QReg> = q_div.iter().rev().collect();
        bit_length_lean(circ, &rev, &t, true);
        let srr = s_rot.to_vec();
        let tr: Vec<&QReg> = t.iter().collect();
        ctrl_add(circ, active, &srr, &tr);
        bit_length_lean(circ, &rev, &t, false);
        xor_const(circ, &t, n_pad);
        for lane in t {
            circ.zero_and_free(lane);
        }
    }
    // 11' rotate_left (was rotate_right restore).
    rotate_left(circ, b, &s_rot[0..rb]);
    // 10' q_div demux (self-inverse XOR).
    let demux_lenders: Vec<&QReg> = a
        .iter()
        .chain(b.iter())
        .chain(extra_lenders.iter().copied())
        .collect();
    set_bit_at_s_gated(circ, q_div, s_rot, active, offset, &demux_lenders);
                                                    // 9' ctrl_sub -> ctrl_add (restore A += B_aligned).
    ctrl_add(circ, active, &aref, &bref);
    // 8' offset clean (self-inverse, _middle); diff in the clz's pa.
    if selective_direct_bitlen_needed(circ, a, b, lo_a, lo_a) {
        direct_bitlen_diff_parity(circ, a, b, lo_a, lo_a, offset, active, extra_lenders);
    } else if lowq_hybrid_clz_enabled() {
        hybrid_bitlen_diff_parity(circ, a, b, lo_a, lo_a, offset, active);
    } else {
        clz_diff_body_middle(circ, a, b, w, lo_a, lo_a, |circ, diff| {
            circ.ccx(active, &diff[0], offset);
        });
    }
    // 7' ctrl_dec -> ctrl_inc.
    ctrl_inc_by_off(circ, active, offset, s_rot, &off_lenders);
    // 6' rotate_left (was rotate_right by offset).
    rotate_one_by_off(circ, b, active, offset, true, &off_lenders);
    // 5' offset compute (self-inverse) -- narrowed, same window as forward.
    if lowq_q958_gated_compare_enabled() {
        narrow_lt_controlled(circ, a, b, offset, active, lo_a);
    } else {
        let or = circ.alloc_qreg("dg.offr");
        narrow_lt(circ, a, b, &or, lo_a);
        circ.ccx(active, &or, offset);
        narrow_lt(circ, a, b, &or, lo_a);
        circ.zero_and_free(or);
    }
    // 4' rotate_right (was rotate_left B<<s).
    rotate_right(circ, b, &s_rot[0..rb]);
    // 3',2',1' clz-mask block (self-inverse, _middle) -- clears s_rot to |0>.
    if selective_direct_bitlen_needed(circ, a, b, lo_a, lo_b) {
        direct_bitlen_diff_update(
            circ,
            a,
            b,
            lo_a,
            lo_b,
            s_rot,
            active,
            offset,
            true,
            extra_lenders,
        );
    } else if lowq_hybrid_clz_enabled() {
        hybrid_bitlen_diff_update(circ, a, b, lo_a, lo_b, s_rot, active, true);
    } else {
        clz_diff_body_middle(circ, a, b, w, lo_a, lo_b, |circ, diff| {
            for j in 0..w {
                circ.ccx(active, &diff[j], &s_rot[j]);
            }
        });
    }
}

/// `out ^= (reg != 0)` (restores reg).
fn or_nonzero(circ: &mut Circuit, reg: &[QReg], out: &QReg) {
    use crate::point_add::trailmix_port::arith::mcx::mcx_clean_k;
    let prev = circ.push_section("p.ornz");
    for q in reg {
        circ.x(q);
    }
    let refs: Vec<&QReg> = reg.iter().collect();
    mcx_clean_k(circ, &refs, out); // out ^= (reg == 0)
    for q in reg {
        circ.x(q);
    }
    circ.x(out); // out ^= (reg != 0)
    circ.pop_section(&prev);
}

/// `out ^= (reg == 0)` via X-bracket + mcx (clean, self-inverse, restores reg).
fn or_is_zero(circ: &mut Circuit, reg: &[QReg], out: &QReg) {
    use crate::point_add::trailmix_port::arith::mcx::mcx_clean_k;
    let prev = circ.push_section("p.orz");
    for q in reg {
        circ.x(q);
    }
    let refs: Vec<&QReg> = reg.iter().collect();
    mcx_clean_k(circ, &refs, out); // out ^= (reg == 0)
    for q in reg {
        circ.x(q);
    }
    circ.pop_section(&prev);
}

fn toggle_zero_dirty(
    circ: &mut Circuit,
    reg: &[QReg],
    out: &QReg,
    candidates: &[&QReg],
    action: &[&QReg],
) {
    for q in reg {
        circ.x(q);
    }
    let controls: Vec<&QReg> = reg.iter().collect();
    dirty_controlled_x(circ, &controls, out, candidates, action);
    for q in reg.iter().rev() {
        circ.x(q);
    }
}

fn toggle_nonzero_dirty(
    circ: &mut Circuit,
    reg: &[QReg],
    out: &QReg,
    candidates: &[&QReg],
    action: &[&QReg],
) {
    toggle_zero_dirty(circ, reg, out, candidates, action);
    circ.x(out);
}

#[allow(clippy::too_many_arguments)]
fn borrowed_swap_in_place(
    circ: &mut Circuit,
    aa: &[QReg],
    bb: &[QReg],
    cca: &[QReg],
    ccb: &[QReg],
    qq: &[QReg],
    counter: &[QReg],
    parity: &QReg,
    s_rot: &[QReg],
    off: &QReg,
) {
    assert!(s_rot.len() >= 2, "Q959 swap predicate lanes");
    let gate = if lowq_q956_off_borrow_enabled() {
        assert_q956_off_alias(off, counter, s_rot);
        assert!(!std::ptr::eq(off, parity), "Q956 off aliases parity");
        &s_rot[2]
    } else {
        off
    };
    let qz = &s_rot[0];
    let anz = &s_rot[1];
    let candidates: Vec<&QReg> = aa
        .iter()
        .chain(bb.iter())
        .chain(cca.iter())
        .chain(ccb.iter())
        .chain(qq.iter())
        .chain(counter.iter())
        .chain(s_rot.iter())
        .chain(std::iter::once(parity))
        .chain(std::iter::once(off))
        .collect();
    let action = [qz, anz, gate];
    let prev = circ.push_section("p.swap.borrowed");

    // At every step boundary s_rot is clean. Retain the two predicates in its
    // first lanes and materialize their active conjunction in a third lane on
    // Q956, leaving the conditionally-clean counter alias untouched.
    toggle_zero_dirty(circ, qq, qz, &candidates, &action);
    toggle_nonzero_dirty(circ, aa, anz, &candidates, &action);
    let toggle_gate = |circ: &mut Circuit| {
        for q in counter {
            circ.x(q);
        }
        let mut controls: Vec<&QReg> = counter.iter().collect();
        controls.push(qz);
        controls.push(anz);
        dirty_controlled_x(circ, &controls, gate, &candidates, &action);
        for q in counter.iter().rev() {
            circ.x(q);
        }
    };
    toggle_gate(circ);
    for j in 0..aa.len() {
        circ.cswap(gate, &aa[j], &bb[j]);
    }
    for j in 0..cca.len() {
        circ.cswap(gate, &cca[j], &ccb[j]);
    }
    circ.cx(gate, parity);
    toggle_gate(circ);
    toggle_nonzero_dirty(circ, aa, anz, &candidates, &action);
    toggle_zero_dirty(circ, qq, qz, &candidates, &action);
    circ.pop_section(&prev);
}

/// WINDOWED multiply substep: same as `multiply_substep_act` but the two clz
/// computations scan the schedule's cofactor clz windows. The `o` clz is on
/// (ca, cb<<s2), both ~bitlen(ca) -> ca window (`ca_window`). The s_rot-clean clz is
/// on (cb, ca) -> cb/ca windows. The cb<<s2 / restore rotates use `rot_bits`.
/// q (ctz) is small -> not windowed. Gate-identical for in-schedule inputs.
#[allow(clippy::too_many_arguments)]
pub fn multiply_substep_windowed(
    circ: &mut Circuit,
    a: &[QReg],
    b: &[QReg],
    q_mul: &[QReg],
    s_rot: &[&QReg],
    off: &QReg,
    active: &QReg,
    extra_lenders: &[&QReg],
    ca_window: usize,
    cb_window: usize,
    rot_bits: usize,
    ctz_bits: usize,
) {
    use crate::point_add::trailmix_port::inversion::shrunken_pz_primitives::ctrl_add;
    let aref: Vec<&QReg> = a.iter().collect();
    let bref: Vec<&QReg> = b.iter().collect();
    let n_pad = q_mul.len();
    let rb = rot_bits.min(s_rot.len());
    let w = s_rot.len();
    let off_lenders: Vec<&QReg> = a
        .iter()
        .chain(b.iter())
        .chain(q_mul.iter())
        .chain(s_rot.iter().copied())
        .chain(extra_lenders.iter().copied())
        .collect();

    if lowq_exact_ctz_enabled() {
        let lenders: Vec<&QReg> = a
            .iter()
            .chain(b.iter())
            .chain(extra_lenders.iter().copied())
            .collect();
        exact_multihot_ctz_deposit(
            circ,
            q_mul,
            &s_rot[..ctz_bits.min(s_rot.len())],
            active,
            &lenders,
        );
    } else {
        let t = circ.alloc_qreg_bits("mg.ctz", w);
        let rev: Vec<&QReg> = q_mul.iter().rev().collect();
        xor_const(circ, &t, n_pad);
        bit_length_lean(circ, &rev, &t, true);
        for j in 0..w {
            circ.ccx(active, &t[j], &s_rot[j]);
        }
        bit_length_lean(circ, &rev, &t, false);
        xor_const(circ, &t, n_pad);
        for lane in t {
            circ.zero_and_free(lane);
        }
    }

    let demux_lenders: Vec<&QReg> = a
        .iter()
        .chain(b.iter())
        .chain(extra_lenders.iter().copied())
        .collect();
    set_bit_at_s_gated(circ, q_mul, s_rot, active, off, &demux_lenders);

    rotate_left(circ, b, &s_rot[0..rb]); // b <<= s if active (bounded rotator)
    ctrl_add(circ, active, &aref, &bref); // a += b<<s if active

    // o = active AND (bitlen(ca) != bitlen(cb<<s2)) -- ca window, _middle; diff in
    // the clz's pa. This clz is the shrunken_pz_divide_forward peak section.
    if selective_direct_bitlen_needed(circ, a, b, ca_window, ca_window) {
        direct_bitlen_diff_parity(
            circ,
            a,
            b,
            ca_window,
            ca_window,
            off,
            active,
            extra_lenders,
        );
    } else if lowq_hybrid_clz_enabled() {
        hybrid_bitlen_diff_parity(circ, a, b, ca_window, ca_window, off, active);
    } else {
        clz_diff_body_middle(circ, a, b, w, ca_window, ca_window, |circ, diff| {
            circ.ccx(active, &diff[0], off);
        });
    }
    rotate_one_by_off(circ, b, active, off, true, &off_lenders); // b <<= 1 if o
    ctrl_inc_by_off(circ, active, off, s_rot, &off_lenders);
    if lowq_q958_gated_compare_enabled() {
        narrow_lt_controlled(circ, a, b, off, active, ca_window);
    } else {
        let lt = circ.alloc_qreg("mg.cleanlt");
        narrow_lt(circ, a, b, &lt, ca_window);
        circ.ccx(active, &lt, off);
        narrow_lt(circ, a, b, &lt, ca_window);
        circ.zero_and_free(lt);
    }
    rotate_right(circ, b, &s_rot[0..rb]); // restore b >>= s_eff (bounded rotator)

    // clean s_rot via _middle clz on (cb, ca): s_rot += (bitlen(cb)-bitlen(ca)).
    if selective_direct_bitlen_needed(circ, b, a, cb_window, ca_window) {
        direct_bitlen_diff_update(
            circ,
            b,
            a,
            cb_window,
            ca_window,
            s_rot,
            active,
            off,
            false,
            extra_lenders,
        );
    } else if lowq_hybrid_clz_enabled() {
        hybrid_bitlen_diff_update(circ, b, a, cb_window, ca_window, s_rot, active, false);
    } else {
        clz_diff_body_middle(circ, b, a, w, cb_window, ca_window, |circ, diff| {
            let srr = s_rot.to_vec();
            let ter: Vec<&QReg> = diff.iter().collect();
            ctrl_add(circ, active, &srr, &ter);
        });
    }
}

/// Gate-by-gate INVERSE of `multiply_substep_windowed` (backward pass). Reverses
/// the sequence; clz/o/q-demux blocks are self-inverse; `rotate_left`<->right,
/// ctrl_add->ctrl_sub, ctrl_inc->ctrl_dec flip. Restores ca -= cb<<s2, re-sets
/// the `q_mul` bit.
#[allow(clippy::too_many_arguments)]
pub fn multiply_substep_windowed_inv(
    circ: &mut Circuit,
    a: &[QReg],
    b: &[QReg],
    q_mul: &[QReg],
    s_rot: &[&QReg],
    off: &QReg,
    active: &QReg,
    extra_lenders: &[&QReg],
    ca_window: usize,
    cb_window: usize,
    rot_bits: usize,
    ctz_bits: usize,
) {
    use crate::point_add::trailmix_port::inversion::shrunken_pz_primitives::{ctrl_add, ctrl_sub};
    let aref: Vec<&QReg> = a.iter().collect();
    let bref: Vec<&QReg> = b.iter().collect();
    let n_pad = q_mul.len();
    let rb = rot_bits.min(s_rot.len());
    let w = s_rot.len();
    let _ = ctrl_add;
    let off_lenders: Vec<&QReg> = a
        .iter()
        .chain(b.iter())
        .chain(q_mul.iter())
        .chain(s_rot.iter().copied())
        .chain(extra_lenders.iter().copied())
        .collect();

    // 10' s_rot clean inverse: ctrl_add -> ctrl_sub (_middle); diff in the clz's pa.
    if selective_direct_bitlen_needed(circ, b, a, cb_window, ca_window) {
        direct_bitlen_diff_update(
            circ,
            b,
            a,
            cb_window,
            ca_window,
            s_rot,
            active,
            off,
            true,
            extra_lenders,
        );
    } else if lowq_hybrid_clz_enabled() {
        hybrid_bitlen_diff_update(circ, b, a, cb_window, ca_window, s_rot, active, true);
    } else {
        clz_diff_body_middle(circ, b, a, w, cb_window, ca_window, |circ, diff| {
            let srr = s_rot.to_vec();
            let ter: Vec<&QReg> = diff.iter().collect();
            ctrl_sub(circ, active, &srr, &ter);
        });
    }
    // 9' rotate_left (was rotate_right restore).
    rotate_left(circ, b, &s_rot[0..rb]);
    // 8' clean-o block (self-inverse) -- narrowed, same window as forward.
    if lowq_q958_gated_compare_enabled() {
        narrow_lt_controlled(circ, a, b, off, active, ca_window);
    } else {
        let lt = circ.alloc_qreg("mg.cleanlt");
        narrow_lt(circ, a, b, &lt, ca_window);
        circ.ccx(active, &lt, off);
        narrow_lt(circ, a, b, &lt, ca_window);
        circ.zero_and_free(lt);
    }
    // 7' ctrl_inc -> ctrl_dec.
    ctrl_dec_by_off(circ, active, off, s_rot, &off_lenders);
    // 6' rotate_right (was rotate_left by o).
    rotate_one_by_off(circ, b, active, off, false, &off_lenders);
    // 5' o clz block (self-inverse, _middle); diff in the clz's pa.
    if selective_direct_bitlen_needed(circ, a, b, ca_window, ca_window) {
        direct_bitlen_diff_parity(
            circ,
            a,
            b,
            ca_window,
            ca_window,
            off,
            active,
            extra_lenders,
        );
    } else if lowq_hybrid_clz_enabled() {
        hybrid_bitlen_diff_parity(circ, a, b, ca_window, ca_window, off, active);
    } else {
        clz_diff_body_middle(circ, a, b, w, ca_window, ca_window, |circ, diff| {
            circ.ccx(active, &diff[0], off);
        });
    }
    // 4' ctrl_add -> ctrl_sub (undo ca += cb<<s2).
    ctrl_sub(circ, active, &aref, &bref);
    // 3' rotate_right (was rotate_left cb<<s2).
    rotate_right(circ, b, &s_rot[0..rb]);
    // 2' q_mul clear demux (self-inverse).
    let demux_lenders: Vec<&QReg> = a
        .iter()
        .chain(b.iter())
        .chain(extra_lenders.iter().copied())
        .collect();
    set_bit_at_s_gated(circ, q_mul, s_rot, active, off, &demux_lenders);
    // 1' clear the least-significant-set-bit index from s_rot.
    if lowq_exact_ctz_enabled() {
        let lenders: Vec<&QReg> = a
            .iter()
            .chain(b.iter())
            .chain(extra_lenders.iter().copied())
            .collect();
        exact_multihot_ctz_erase(
            circ,
            q_mul,
            &s_rot[..ctz_bits.min(s_rot.len())],
            active,
            &lenders,
        );
    } else {
        let t = circ.alloc_qreg_bits("mg.ctz", w);
        let rev: Vec<&QReg> = q_mul.iter().rev().collect();
        xor_const(circ, &t, n_pad);
        bit_length_lean(circ, &rev, &t, true);
        for j in 0..w {
            circ.ccx(active, &t[j], &s_rot[j]);
        }
        bit_length_lean(circ, &rev, &t, false);
        xor_const(circ, &t, n_pad);
        for lane in t {
            circ.zero_and_free(lane);
        }
    }
}

// NEXT (reversible_pz_notes.md has the primitive mapping):
//   fn normalize_input(circ, x, sgn)               -- x -> min(x,P-x), set sgn
//   fn division_substep(circ, regs, flags, s, bound)
//   fn multiply_substep(circ, regs, flags, s, bound)
//   fn transition(circ, regs, flags)
//   fn iterate(circ, regs, flags, n_iters)         -- the fixed-count driver
//   fn recover_inverse(circ, regs, flags)          -- parity^sgn sign fix
//   test pz_sm_faithful  -- per-iter contract vs a Rust port of pz_big_step

// ===== shrunken_pz reversible inversion step driver (shared fwd/back, used by
// the round-trip test AND the EC-add) =====

// ---- shared forward/backward step helpers (used by the round-trip) ----

/// Compute `g = active AND (x < y)` directly from the comparator carry, run the
/// gated body, then clear `g` with the same restored-input comparison. No
/// separate `(x < y)` lane is retained across the body.
pub(crate) fn gate_hold(
    c: &mut Circuit,
    x: &[QReg],
    y: &[QReg],
    active: &QReg,
    g: &QReg,
    borrowed_carry: Option<&QReg>,
    body: impl FnOnce(&mut Circuit, &QReg),
) {
    let xr: Vec<&QReg> = x.iter().collect();
    let yr: Vec<&QReg> = y.iter().collect();
    let compare = |c: &mut Circuit| {
        if let Some(carry) = borrowed_carry {
            borrow_compare_gated_refs_with_carry(c, &xr, &yr, active, g, carry);
        } else {
            borrow_compare_gated_refs(c, &xr, &yr, active, g);
        }
    };
    compare(c);
    body(c, g);
    compare(c);
}

/// Run a gated body with `g = (counter == 0) AND (x < y)` while allocating
/// only `g`. The existing parity lane temporarily hosts the active predicate;
/// its prior value is parked in `g`, swapped back before the body, and restored
/// exactly after the second comparison. Two clean shift lanes host the
/// comparator output and carry only while the body is not running.
#[allow(clippy::too_many_arguments)]
fn gate_hold_counter_zero(
    c: &mut Circuit,
    x: &[QReg],
    y: &[QReg],
    counter: &[QReg],
    parity: &QReg,
    s_rot: &[QReg],
    g: &QReg,
    candidates: &[&QReg],
    body: impl FnOnce(&mut Circuit, &QReg),
) {
    assert!(s_rot.len() >= 2, "Q959 comparator borrow lanes");
    let lt = &s_rot[0];
    let carry = &s_rot[1];
    let action = [g, parity, lt, carry];
    let xr: Vec<&QReg> = x.iter().collect();
    let yr: Vec<&QReg> = y.iter().collect();

    let swap_parity_gate = |c: &mut Circuit| {
        c.cx(parity, g);
        c.cx(g, parity);
        c.cx(parity, g);
    };
    let toggle_active = |c: &mut Circuit| {
        if counter.is_empty() {
            c.x(parity);
        } else {
            toggle_zero_dirty(c, counter, parity, candidates, &action);
        }
    };
    let compare = |c: &mut Circuit| {
        borrow_compare_refs_with_carry(c, &xr, &yr, lt, carry);
    };
    let remove_nonless_active = |c: &mut Circuit| {
        for q in counter {
            c.x(q);
        }
        c.x(lt);
        let mut controls: Vec<&QReg> = counter.iter().collect();
        controls.push(lt);
        dirty_controlled_x(c, &controls, g, candidates, &action);
        c.x(lt);
        for q in counter.iter().rev() {
            c.x(q);
        }
    };

    // Park P in g, clear parity, compute active in parity, then swap the two:
    // parity=P and g=active. Remove the active-and-not-less branch to obtain
    // g=active-and-less.
    c.cx(parity, g);
    c.cx(g, parity);
    toggle_active(c);
    swap_parity_gate(c);
    compare(c);
    remove_nonless_active(c);
    compare(c);

    body(c, g);

    // Exact reverse of the preparation above.
    compare(c);
    remove_nonless_active(c);
    compare(c);
    swap_parity_gate(c);
    toggle_active(c);
    c.cx(g, parity);
    c.cx(parity, g);
}

/// done-counter (forward: counter += conv) / its inverse (counter -= conv),
/// conv = (A==0 & q==0). `done` is clean scratch (|0> at exit). User's recipe.
pub(crate) fn done_counter_fn(
    c: &mut Circuit,
    aa: &[QReg],
    qq: &[QReg],
    counter: &[QReg],
    s_rot: &[QReg],
    off: &QReg,
    candidates: &[&QReg],
    inverse: bool,
) {
    if counter.is_empty() {
        return;
    }
    if lowq_q959_selective_borrow_enabled() {
        assert!(s_rot.len() >= 2, "Q959 done predicate lanes");
        let done = if lowq_q956_off_borrow_enabled() {
            // Boundary logic gets a truly clean lane; counter[0] is reserved
            // for conditionally-clean borrowing only inside active bodies.
            assert_q956_off_alias(off, counter, s_rot);
            &s_rot[2]
        } else {
            off
        };
        let az = &s_rot[0];
        let qz = &s_rot[1];
        let counter_refs: Vec<&QReg> = counter.iter().collect();
        let action = [done, az, qz];
        let conv = |c: &mut Circuit| {
            toggle_zero_dirty(c, aa, az, candidates, &action);
            toggle_zero_dirty(c, qq, qz, candidates, &action);
            c.ccx(az, qz, done);
            toggle_zero_dirty(c, qq, qz, candidates, &action);
            toggle_zero_dirty(c, aa, az, candidates, &action);
        };
        let cnz = |c: &mut Circuit| {
            toggle_nonzero_dirty(c, counter, az, candidates, &action);
            c.cx(az, done);
            toggle_nonzero_dirty(c, counter, az, candidates, &action);
        };
        if inverse {
            cnz(c);
            dirty_controlled_inc_suffix(c, &[done], &counter_refs, 0, true, candidates);
            conv(c);
        } else {
            conv(c);
            dirty_controlled_inc_suffix(c, &[done], &counter_refs, 0, false, candidates);
            cnz(c);
        }
        return;
    }

    let done = c.alloc_qreg("done");
    let conv = |c: &mut Circuit, done: &QReg| {
        let az = c.alloc_qreg("d.az");
        let qz = c.alloc_qreg("d.qz");
        or_is_zero(c, aa, &az);
        or_is_zero(c, qq, &qz);
        c.ccx(&az, &qz, done); // done ^= (A==0 & q==0)
        or_is_zero(c, qq, &qz);
        or_is_zero(c, aa, &az);
        c.zero_and_free(qz);
        c.zero_and_free(az);
    };
    let cnz = |c: &mut Circuit, done: &QReg| {
        let z = c.alloc_qreg("d.cnz");
        or_nonzero(c, counter, &z);
        c.cx(&z, done); // done ^= (counter != 0)
        or_nonzero(c, counter, &z);
        c.zero_and_free(z);
    };
    if inverse {
        cnz(c, &done);
        ctrl_dec(c, &done, counter);
        conv(c, &done);
    } else {
        conv(c, &done);
        ctrl_inc(c, &done, counter);
        cnz(c, &done);
    }
    c.zero_and_free(done);
}

/// One forward (inverse=false) or backward (inverse=true) `shrunken_pz` step on the
/// dynamic-W registers at their current width. Resize is done by the caller.
#[allow(clippy::too_many_arguments)]
pub(crate) fn shrunken_pz_pass_step(
    c: &mut Circuit,
    aa: &[QReg],
    bb: &[QReg],
    cca: &[QReg],
    ccb: &[QReg],
    qq: &[QReg],
    counter: &[QReg],
    parity: &QReg,
    s_rot: &[QReg],
    off: &QReg,
    i: usize,
    inverse: bool,
) {
    use crate::point_add::trailmix_port::inversion::shrunken_pz_schedule::{reg_los, shift_bounds};
    fn rb(b: usize) -> usize {
        if b == 0 {
            1
        } else {
            64 - (b as u64).leading_zeros() as usize
        }
    }
    let (lo_a, lo_b, ca_window, cb_window, _) = reg_los(i);
    let (sdb, s2b) = shift_bounds(i);
    let ctz_bits = q954_ctz_width(i);
    if lowq_q954_srot_counter7_enabled() {
        assert_eq!(s_rot.len(), 4, "Q954 boundary must see four owned shift lanes");
        assert_eq!(counter.len(), 8, "Q954 boundary counter width");
        assert!(
            s_rot.iter().all(|lane| !std::ptr::eq(lane, &counter[7])),
            "Q954 boundary received the arithmetic-only counter[7] alias"
        );
        assert_eq!(
            ctz_bits,
            if i <= Q954_LAST_CTZ_BIT4_ROW { 5 } else { 4 },
            "Q954 CTZ bit4 cutoff drift"
        );
    }
    if lowq_q956_off_borrow_enabled() {
        assert_q956_off_alias(off, counter, s_rot);
        assert!(!std::ptr::eq(off, parity), "Q956 off aliases parity");
        assert!(
            aa.iter()
                .chain(bb.iter())
                .chain(cca.iter())
                .chain(ccb.iter())
                .chain(qq.iter())
                .all(|lane| !std::ptr::eq(off, lane)),
            "Q956 off alias overlaps a dynamic state register"
        );
    }
    // Swap, gated g_swap=(q==0 & A!=0 & active). HOLD the (q==0)/(A!=0) flags
    // across the cswaps so or_nonzero(A)/or_is_zero(q) run 2x not 4x per step
    // (the swap preserves both predicates: q untouched, A_new=B_old!=0).
    let swap = |c: &mut Circuit, active: &QReg| {
        let qz = c.alloc_qreg("sw.qz");
        let anz = c.alloc_qreg("sw.anz");
        or_is_zero(c, qq, &qz);
        or_nonzero(c, aa, &anz);
        let t = c.alloc_qreg("sw.t");
        let g = c.alloc_qreg("g_swap");
        c.ccx(&qz, &anz, &t); // t = (q==0 & A!=0)
        c.ccx(&t, active, &g); // g_swap = t AND active
        for j in 0..aa.len() {
            c.cswap(&g, &aa[j], &bb[j]);
        }
        for j in 0..cca.len() {
            c.cswap(&g, &cca[j], &ccb[j]);
        }
        c.cx(&g, parity);
        c.ccx(&t, active, &g); // uncompute g (t,active preserved)
        c.ccx(&qz, &anz, &t); // uncompute t (qz held; anz=A_old!=0)
        c.zero_and_free(g);
        c.zero_and_free(t);
        or_nonzero(c, aa, &anz); // post-swap A=B_old!=0 -> clears anz
        or_is_zero(c, qq, &qz);
        c.zero_and_free(anz);
        c.zero_and_free(qz);
    };

    let boundary_candidates: Vec<&QReg> = aa
        .iter()
        .chain(bb.iter())
        .chain(cca.iter())
        .chain(ccb.iter())
        .chain(qq.iter())
        .chain(counter.iter())
        .chain(s_rot.iter())
        .chain(std::iter::once(parity))
        .chain(std::iter::once(off))
        .collect();

    if lowq_q959_selective_borrow_enabled() {
        if inverse {
            done_counter_fn(
                c,
                aa,
                qq,
                counter,
                s_rot,
                off,
                &boundary_candidates,
                true,
            );
            borrowed_swap_in_place(c, aa, bb, cca, ccb, qq, counter, parity, s_rot, off);

            let g_div = c.alloc_qreg("g_div");
            gate_hold_counter_zero(
                c,
                cca,
                ccb,
                counter,
                parity,
                s_rot,
                &g_div,
                &boundary_candidates,
                |c, g| {
                    let lenders: Vec<&QReg> = cca.iter().chain(ccb.iter()).collect();
                    with_arithmetic_srot_view(s_rot, counter, |s_rot_arith| {
                        division_substep_windowed_inv(
                            c,
                            aa,
                            bb,
                            qq,
                            s_rot_arith,
                            off,
                            g,
                            &lenders,
                            lo_a,
                            lo_b,
                            rb(sdb),
                            ctz_bits,
                        );
                    });
                },
            );
            c.zero_and_free(g_div);

            let g_mul = c.alloc_qreg("g_mul");
            gate_hold_counter_zero(
                c,
                aa,
                bb,
                counter,
                parity,
                s_rot,
                &g_mul,
                &boundary_candidates,
                |c, g| {
                    let lenders: Vec<&QReg> = aa.iter().chain(bb.iter()).collect();
                    with_arithmetic_srot_view(s_rot, counter, |s_rot_arith| {
                        multiply_substep_windowed_inv(
                            c,
                            cca,
                            ccb,
                            qq,
                            s_rot_arith,
                            off,
                            g,
                            &lenders,
                            ca_window,
                            cb_window,
                            rb(s2b),
                            ctz_bits,
                        );
                    });
                },
            );
            c.zero_and_free(g_mul);
        } else {
            let g_mul = c.alloc_qreg("g_mul");
            gate_hold_counter_zero(
                c,
                aa,
                bb,
                counter,
                parity,
                s_rot,
                &g_mul,
                &boundary_candidates,
                |c, g| {
                    let lenders: Vec<&QReg> = aa.iter().chain(bb.iter()).collect();
                    with_arithmetic_srot_view(s_rot, counter, |s_rot_arith| {
                        multiply_substep_windowed(
                            c,
                            cca,
                            ccb,
                            qq,
                            s_rot_arith,
                            off,
                            g,
                            &lenders,
                            ca_window,
                            cb_window,
                            rb(s2b),
                            ctz_bits,
                        );
                    });
                },
            );
            c.zero_and_free(g_mul);

            let g_div = c.alloc_qreg("g_div");
            gate_hold_counter_zero(
                c,
                cca,
                ccb,
                counter,
                parity,
                s_rot,
                &g_div,
                &boundary_candidates,
                |c, g| {
                    let lenders: Vec<&QReg> = cca.iter().chain(ccb.iter()).collect();
                    with_arithmetic_srot_view(s_rot, counter, |s_rot_arith| {
                        division_substep_windowed(
                            c,
                            aa,
                            bb,
                            qq,
                            s_rot_arith,
                            off,
                            g,
                            &lenders,
                            lo_a,
                            lo_b,
                            rb(sdb),
                            ctz_bits,
                        );
                    });
                },
            );
            c.zero_and_free(g_div);

            borrowed_swap_in_place(c, aa, bb, cca, ccb, qq, counter, parity, s_rot, off);
            done_counter_fn(
                c,
                aa,
                qq,
                counter,
                s_rot,
                off,
                &boundary_candidates,
                false,
            );
        }
        return;
    }

    if inverse {
        done_counter_fn(c, aa, qq, counter, s_rot, off, &boundary_candidates, true);
        let active = compute_active(c, counter, &boundary_candidates);
        swap(c, &active); // self-inverse
        let g_div = c.alloc_qreg("g_div");
        gate_hold(
            c,
            cca,
            ccb,
            &active,
            &g_div,
            lowq_q959_selective_borrow_enabled().then_some(&s_rot[0]),
            |c, g| {
            let lenders: Vec<&QReg> = cca.iter().chain(ccb.iter()).collect();
            with_arithmetic_srot_view(s_rot, counter, |s_rot_arith| {
                division_substep_windowed_inv(
                    c, aa, bb, qq, s_rot_arith, off, g, &lenders, lo_a, lo_b, rb(sdb),
                    ctz_bits,
                );
            });
            },
        );
        c.zero_and_free(g_div);
        let g_mul = c.alloc_qreg("g_mul");
        gate_hold(
            c,
            aa,
            bb,
            &active,
            &g_mul,
            lowq_q959_selective_borrow_enabled().then_some(&s_rot[0]),
            |c, g| {
            let lenders: Vec<&QReg> = aa.iter().chain(bb.iter()).collect();
            with_arithmetic_srot_view(s_rot, counter, |s_rot_arith| {
                multiply_substep_windowed_inv(
                    c,
                    cca,
                    ccb,
                    qq,
                    s_rot_arith,
                    off,
                    g,
                    &lenders,
                    ca_window,
                    cb_window,
                    rb(s2b),
                    ctz_bits,
                );
            });
            },
        );
        c.zero_and_free(g_mul);
        uncompute_active(c, counter, &active, &boundary_candidates);
        c.zero_and_free(active);
    } else {
        let active = compute_active(c, counter, &boundary_candidates);
        let g_mul = c.alloc_qreg("g_mul");
        gate_hold(
            c,
            aa,
            bb,
            &active,
            &g_mul,
            lowq_q959_selective_borrow_enabled().then_some(&s_rot[0]),
            |c, g| {
            let lenders: Vec<&QReg> = aa.iter().chain(bb.iter()).collect();
            with_arithmetic_srot_view(s_rot, counter, |s_rot_arith| {
                multiply_substep_windowed(
                    c,
                    cca,
                    ccb,
                    qq,
                    s_rot_arith,
                    off,
                    g,
                    &lenders,
                    ca_window,
                    cb_window,
                    rb(s2b),
                    ctz_bits,
                );
            });
            },
        );
        c.zero_and_free(g_mul);
        let g_div = c.alloc_qreg("g_div");
        gate_hold(
            c,
            cca,
            ccb,
            &active,
            &g_div,
            lowq_q959_selective_borrow_enabled().then_some(&s_rot[0]),
            |c, g| {
            let lenders: Vec<&QReg> = cca.iter().chain(ccb.iter()).collect();
            with_arithmetic_srot_view(s_rot, counter, |s_rot_arith| {
                division_substep_windowed(
                    c, aa, bb, qq, s_rot_arith, off, g, &lenders, lo_a, lo_b, rb(sdb),
                    ctz_bits,
                );
            });
            },
        );
        c.zero_and_free(g_div);
        swap(c, &active);
        uncompute_active(c, counter, &active, &boundary_candidates);
        c.zero_and_free(active);
        done_counter_fn(c, aa, qq, counter, s_rot, off, &boundary_candidates, false);
    }
}

/// Resize a dynamic-W register to `target` bits: free high qubits (must be |0>)
/// or alloc fresh |0> ones, in place.
pub(crate) fn shrunken_pz_resize(c: &mut Circuit, reg: &mut Vec<QReg>, target: usize, name: &str) {
    while reg.len() > target {
        let q = reg.pop().unwrap();
        c.zero_and_free(q);
    }
    while reg.len() < target {
        let k = reg.len();
        reg.push(c.alloc_qreg(&format!("{name}[{k}]")));
    }
}

/// Drop the proven-clean top lane of the canonical persistent slope while the
/// backward EEA owns the peak. The canonical multiplier and field negation
/// guarantee the precondition `lambda[256] = |0>`.
fn release_q955_canonical_lambda_top(c: &mut Circuit, lambda: &mut Vec<QReg>) {
    assert_eq!(
        lambda.len(),
        257,
        "Q955 canonical lambda must enter the reverse EEA with 257 lanes"
    );
    let active_before = c.b.active_qubits;
    let top = lambda.pop().expect("Q955 canonical lambda top lane");
    c.zero_and_free(top);
    assert_eq!(lambda.len(), 256, "Q955 reverse EEA keeps 256 lambda lanes");
    assert_eq!(
        c.b.active_qubits + 1,
        active_before,
        "Q955 canonical lambda release must save exactly one live qubit"
    );
}

/// Restore the public 257-bit slope shape with a newly allocated clean lane.
fn restore_q955_canonical_lambda_top(c: &mut Circuit, lambda: &mut Vec<QReg>) {
    assert_eq!(
        lambda.len(),
        256,
        "Q955 canonical lambda must leave the reverse EEA with 256 lanes"
    );
    let active_before = c.b.active_qubits;
    lambda.push(c.alloc_qreg("shpzdiv.lambda[256].restored"));
    assert_eq!(lambda.len(), 257, "Q955 lambda API requires 257 lanes");
    assert_eq!(
        c.b.active_qubits,
        active_before + 1,
        "Q955 lambda top restoration must allocate exactly one clean qubit"
    );
}

/// Remove a canonical field passenger's proven-zero 257th lane while an EEA
/// traversal owns the peak. The lane is deallocated, not borrowed as workspace.
fn release_q954_canonical_passenger_top(
    c: &mut Circuit,
    passenger: &mut Vec<QReg>,
    context: &str,
) {
    assert!(lowq_q954_srot_counter7_enabled());
    assert_eq!(
        passenger.len(),
        257,
        "Q954 {context} passenger must enter EEA with 257 lanes"
    );
    let live_before = c.b.active_qubits as usize;
    let top = passenger.pop().expect("Q954 canonical passenger top lane");
    c.zero_and_free(top);
    assert_eq!(passenger.len(), 256, "Q954 {context} passenger width");
    c.flush_pending_frees();
    assert_eq!(
        c.b.active_qubits as usize + 1,
        live_before,
        "Q954 {context} passenger release must save one live qubit"
    );
}

fn restore_q954_canonical_passenger_top(
    c: &mut Circuit,
    passenger: &mut Vec<QReg>,
    name: &str,
) {
    assert!(lowq_q954_srot_counter7_enabled());
    assert_eq!(
        passenger.len(),
        256,
        "Q954 canonical passenger must leave EEA with 256 lanes"
    );
    let live_before = c.b.active_qubits as usize;
    passenger.push(c.alloc_qreg(name));
    assert_eq!(passenger.len(), 257, "Q954 passenger API requires 257 lanes");
    assert_eq!(
        c.b.active_qubits as usize,
        live_before + 1,
        "Q954 passenger top restoration must allocate one clean qubit"
    );
}

fn shrunken_pz_shrink(c: &mut Circuit, reg: &mut Vec<QReg>, target: usize) {
    while reg.len() > target {
        let q = reg.pop().unwrap();
        c.zero_and_free(q);
    }
}

/// FORWARD `shrunken_pz` inversion driver. PRE: the registers hold the `S_0` state at width
/// `reg_widths(0)` -- A=p, B=|x| (sign-adjusted, < p/2), ca=0, cb=1, q=0,
/// counter=0, parity=1. Runs all `SHRUNKEN_PZ_NSTEPS` forward steps (resizing per step),
/// leaving the modular inverse of |x| in `ccb` (up to the `parity` bit: the true
/// value is `parity ? cb : p-cb`), with A=p, B=|x| at the EEA terminal. `s`,
/// `s_rot` (9 bits each), `off`, `parity`, `counter` (10 bits) are fixed-width.
#[allow(clippy::too_many_arguments)]
pub(crate) fn shrunken_pz_invert_forward(
    c: &mut Circuit,
    aa: &mut Vec<QReg>,
    bb: &mut Vec<QReg>,
    cca: &mut Vec<QReg>,
    ccb: &mut Vec<QReg>,
    qq: &mut Vec<QReg>,
    counter: &[QReg],
    parity: &QReg,
    s_rot: &[QReg],
    off: &QReg,
) {
    use crate::point_add::trailmix_port::inversion::shrunken_pz_schedule::{reg_widths, SHRUNKEN_PZ_NSTEPS};
    for i in 0..SHRUNKEN_PZ_NSTEPS {
        let (wa, wb, wca, wcb, wq) = reg_widths(i);
        let wab = trailmix_ab_width(wa.max(wb));
        let wcacb = trailmix_cacb_width(wca.max(wcb));
        let wq = trailmix_q_width_step(wq, wa, wb, wca, wcb);
        // Rebalance the pack transactionally: release every shrinking high
        // lane before allocating any growing lane. The previous interleaving
        // created a one-qubit transient above the final schedule width.
        shrunken_pz_shrink(c, aa, wab);
        shrunken_pz_shrink(c, bb, wab);
        shrunken_pz_shrink(c, cca, wcacb);
        shrunken_pz_shrink(c, ccb, wcacb);
        shrunken_pz_shrink(c, qq, wq);
        shrunken_pz_resize(c, aa, wab, "A");
        shrunken_pz_resize(c, bb, wab, "B");
        shrunken_pz_resize(c, cca, wcacb, "ca");
        shrunken_pz_resize(c, ccb, wcacb, "cb");
        shrunken_pz_resize(c, qq, wq, "q");
        shrunken_pz_pass_step(
            c, aa, bb, cca, ccb, qq, counter, parity, s_rot, off, i, false,
        );
    }
}

/// BACKWARD `shrunken_pz` inversion driver (gate-for-gate inverse of `shrunken_pz_invert_forward`).
/// Restores the `S_0` state (A=p, B=|x|, ca=0, cb=1, q=0, counter=0, parity=1) and
/// uncomputes the inverse from `ccb`. Resizes back down per step.
#[allow(clippy::too_many_arguments)]
pub(crate) fn shrunken_pz_invert_backward(
    c: &mut Circuit,
    aa: &mut Vec<QReg>,
    bb: &mut Vec<QReg>,
    cca: &mut Vec<QReg>,
    ccb: &mut Vec<QReg>,
    qq: &mut Vec<QReg>,
    counter: &[QReg],
    parity: &QReg,
    s_rot: &[QReg],
    off: &QReg,
) {
    use crate::point_add::trailmix_port::inversion::shrunken_pz_schedule::{reg_widths, SHRUNKEN_PZ_NSTEPS};
    for i in (0..SHRUNKEN_PZ_NSTEPS).rev() {
        shrunken_pz_pass_step(
            c, aa, bb, cca, ccb, qq, counter, parity, s_rot, off, i, true,
        );
        if i > 0 {
            let (wa, wb, wca, wcb, wq) = reg_widths(i - 1);
            let wab = trailmix_ab_width(wa.max(wb));
            let wcacb = trailmix_cacb_width(wca.max(wcb));
            let wq = trailmix_q_width_step(wq, wa, wb, wca, wcb);
            shrunken_pz_shrink(c, aa, wab);
            shrunken_pz_shrink(c, bb, wab);
            shrunken_pz_shrink(c, cca, wcacb);
            shrunken_pz_shrink(c, ccb, wcacb);
            shrunken_pz_shrink(c, qq, wq);
            shrunken_pz_resize(c, aa, wab, "A");
            shrunken_pz_resize(c, bb, wab, "B");
            shrunken_pz_resize(c, cca, wcacb, "ca");
            shrunken_pz_resize(c, ccb, wcacb, "cb");
            shrunken_pz_resize(c, qq, wq, "q");
        }
    }
}

/// `lambda = dy / dx mod p`, with `dx` and `dy` PRESERVED. `dx`, `dy` are 257-bit
/// registers holding field elements in [0, p). Returns `(dx, dy, lambda)` -- dx
/// and dy unchanged (dy reconstructed via the HMR-ghost trick), lambda = dy·dx^-1.
/// With `LOWQ_Q955_OFF_CANONICAL=1`, lambda is produced by the exact canonical
/// multiplier and its clean top lane is absent during reverse EEA. The API still
/// returns 257 lanes by appending a new clean top lane afterward.
/// With `LOWQ_Q954_SROT_COUNTER7=1`, canonical dy[256] is likewise absent during
/// the initial forward EEA and restored only after its constant pack is removed.
pub fn shrunken_pz_divide_forward(
    c: &mut Circuit,
    mut dx: Vec<QReg>,
    mut dy: Vec<QReg>,
) -> (Vec<QReg>, Vec<QReg>, Vec<QReg>) {
    use crate::point_add::trailmix_port::arith::compare::compare_geq_const;
    use crate::point_add::trailmix_port::arith::rfold_mbu::{
        mod_mul_canonical_mbu, mod_mul_rfold_mbu,
    };
    use crate::point_add::trailmix_port::inversion::shrunken_pz_schedule::reg_widths;
    use crate::point_add::trailmix_port::num_bigint::BigUint;
    assert_eq!(dx.len(), 257);
    assert_eq!(dy.len(), 257);
    let canonical_lambda_top_off = lowq_q955_off_canonical_enabled();
    let canonical_passenger_top_off = lowq_q954_srot_counter7_enabled();
    if canonical_passenger_top_off {
        release_q954_canonical_passenger_top(c, &mut dy, "divide-forward dy");
    }
    // sgn = dx > p/2  <=>  dx >= (p+1)/2.
    let half_bytes = vec![
        0x18, 0xfe, 0xff, 0x7f, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f, 0x00,
    ];
    let p_bytes = crate::point_add::trailmix_port::mod_arith::SECP256K1_P_LE;

    // --- sign-adjust dx -> |dx| < p/2 (the schedule assumes |x| < p/2) ---
    let reuse_sign_wire = sign_parity_q_reuse_enabled();
    let mut fused_parity = reuse_sign_wire.then(|| c.alloc_qreg("shpzdiv.par_sgn"));
    let sgn = (!reuse_sign_wire).then(|| c.alloc_qreg("shpzdiv.sgn"));
    let sign_control = fused_parity.as_ref().or(sgn.as_ref()).unwrap();
    compare_geq_const(c, &dx, &half_bytes, sign_control);
    controlled_field_neg(c, sign_control, &dx); // dx := (sgn ? p-dx : dx) = |dx|

    // --- set up the inversion S_0 state (B = |dx|, A = p, cb = 1, parity = 1) ---
    let (a0, b0, ca0, cb0, q0) = reg_widths(0);
    let (wg0, wc0) = (a0.max(b0), ca0.max(cb0));
    shrunken_pz_resize(c, &mut dx, wg0, "B"); // |dx| becomes the EEA B register
    let mut aa = c.alloc_qreg_bits("shpzdiv.A", wg0);
    let mut cca = c.alloc_qreg_bits("shpzdiv.ca", wc0);
    let mut ccb = c.alloc_qreg_bits("shpzdiv.cb", wc0);
    let mut qq = c.alloc_qreg_bits("shpzdiv.q", q0.max(1));
    let s_rot = c.alloc_qreg_bits("shpzdiv.srot", trailmix_srot_width());
    let parity = fused_parity
        .take()
        .unwrap_or_else(|| c.alloc_qreg("shpzdiv.par"));
    let counter = c.alloc_qreg_bits("shpzdiv.ctr", trailmix_counter_width());
    let off_owned = (!lowq_q956_off_borrow_enabled()).then(|| c.alloc_qreg("shpzdiv.off"));
    let off = off_owned.as_ref().unwrap_or_else(|| {
        assert_q956_off_alias(&counter[0], &counter, &s_rot);
        &counter[0]
    });
    let load_p = |c: &mut Circuit, reg: &[QReg]| {
        for (j, q) in reg.iter().enumerate() {
            if j < 256 && (p_bytes[j / 8] >> (j % 8)) & 1 == 1 {
                c.x(q);
            }
        }
    };
    load_p(c, &aa); // A = p
    c.x(&ccb[0]); // cb = 1
    c.x(&parity); // parity = 1, or fused parity = 1 XOR sign

    // --- forward inversion: 1/|dx| in cb (up to the parity bit) ---
    shrunken_pz_invert_forward(
        c, &mut aa, &mut dx, &mut cca, &mut ccb, &mut qq, &counter, &parity, &s_rot, off,
    );

    // --- TEAR DOWN the EEA pack before creating lambda. At convergence the PZ
    // state is A=0, B=1, ca=p, q=0 (all CONSTANTS) and cb=1/|dx| (the only data).
    // Free the constant registers (0-Toffoli uncompute) so only cb is live during
    // the multiply -- saves ~ca(258) qubits at the peak. Re-create them (cheap)
    // before the backward. ---
    let (ta, tb, tca, tq) = (aa.len(), dx.len(), cca.len(), qq.len());
    load_p(c, &cca); // ca: p -> 0
    c.x(&dx[0]); // B: 1 -> 0
    for q in std::mem::take(&mut aa) {
        c.zero_and_free(q); // A = 0
    }
    for q in std::mem::take(&mut dx) {
        c.zero_and_free(q); // B = 0
    }
    for q in std::mem::take(&mut cca) {
        c.zero_and_free(q); // ca = 0
    }
    for q in std::mem::take(&mut qq) {
        c.zero_and_free(q); // q = 0
    }
    if canonical_passenger_top_off {
        restore_q954_canonical_passenger_top(c, &mut dy, "shpzdiv.dy[256]");
    }
    // --- lambda = dy * (1/|dx|), parity/sign corrected (only cb live in the pack) ---
    let cb_w = ccb.len();
    shrunken_pz_resize(c, &mut ccb, 257, "cb"); // pad the inverse to 257 for mod_mul
    let mut lambda = c.alloc_qreg_bits("shpzdiv.lambda", 257);
    if canonical_lambda_top_off {
        mod_mul_canonical_mbu(c, &lambda, &ccb[..257], &dy);
    } else {
        mod_mul_rfold_mbu(c, &lambda, &ccb[..257], &dy); // lambda_raw = dy * cb
    }
    shrunken_pz_resize(c, &mut ccb, cb_w, "cb"); // restore width for the backward
    // 1/dx = (-1)^{sgn + (1-parity)} * cb. With fusion, the live parity
    // lane already equals sgn XOR parity, so its X-bracket is exactly f.
    if reuse_sign_wire {
        c.x(&parity);
        if canonical_lambda_top_off {
            controlled_field_neg_canonical(c, &parity, &lambda);
        } else {
            controlled_field_neg(c, &parity, &lambda);
        }
        c.x(&parity);
    } else {
        let sgn = sgn.as_ref().unwrap();
        let f = c.alloc_qreg("shpzdiv.negf");
        c.cx(sgn, &f);
        c.cx(&parity, &f);
        c.x(&f); // f = NOT(sgn XOR parity)
        if canonical_lambda_top_off {
            controlled_field_neg_canonical(c, &f, &lambda);
        } else {
            controlled_field_neg(c, &f, &lambda);
        }
        c.x(&f);
        c.cx(&parity, &f);
        c.cx(sgn, &f); // uncompute f
        c.zero_and_free(f);
    }

    if canonical_lambda_top_off {
        release_q955_canonical_lambda_top(c, &mut lambda);
    }

    // --- GHOST dy (HMR each bit) so the reverse runs dy-free ---
    let mut ghosts = Vec::with_capacity(dy.len());
    for q in &dy {
        ghosts.push(c.hmr_ghost(q));
    }
    for q in dy {
        c.zero_and_free(q);
    }
    // --- RE-CREATE the constant pack (A=0, B=1, ca=p, q=0) for the backward ---
    aa = c.alloc_qreg_bits("shpzdiv.A", ta); // A = 0
    dx = c.alloc_qreg_bits("shpzdiv.B", tb);
    c.x(&dx[0]); // B = 1
    cca = c.alloc_qreg_bits("shpzdiv.ca", tca);
    load_p(c, &cca); // ca = p
    qq = c.alloc_qreg_bits("shpzdiv.q", tq); // q = 0

    // --- backward inversion: restore B = |dx|, uncompute cb/parity ---
    shrunken_pz_invert_backward(
        c, &mut aa, &mut dx, &mut cca, &mut ccb, &mut qq, &counter, &parity, &s_rot, off,
    );

    // --- free the clean inversion ancillas (S_0: A=p, ca=0, cb=1, q=0) ---
    if !reuse_sign_wire {
        c.x(&parity);
    }
    c.x(&ccb[0]); // cb: 1 -> 0
    load_p(c, &aa); // A: p -> 0
    for q in aa.into_iter().chain(cca).chain(ccb).chain(qq) {
        c.zero_and_free(q);
    }
    if let Some(off) = off_owned {
        c.zero_and_free(off);
    }
    for q in s_rot.into_iter().chain(counter) {
        c.zero_and_free(q);
    }

    // --- un-sign-adjust: |dx| -> dx, uncompute sign state ---
    shrunken_pz_resize(c, &mut dx, 257, "dx");
    if reuse_sign_wire {
        // Reverse EEA restored fused parity = 1 XOR sign.
        c.x(&parity);
        controlled_field_neg(c, &parity, &dx);
        compare_geq_const(c, &dx, &half_bytes, &parity);
        c.zero_and_free(parity);
    } else {
        let sgn = sgn.unwrap();
        controlled_field_neg(c, &sgn, &dx);
        compare_geq_const(c, &dx, &half_bytes, &sgn);
        c.zero_and_free(sgn);
        c.zero_and_free(parity);
    }

    // --- reconstruct dy = lambda * dx and EXORCIZE the ghosts ---
    if canonical_lambda_top_off {
        restore_q955_canonical_lambda_top(c, &mut lambda);
    }
    assert_eq!(lambda.len(), 257, "slope API and raw dy roundtrip require 257 lanes");
    let dy_new = c.alloc_qreg_bits("shpzdiv.dy", 257);
    if canonical_lambda_top_off {
        mod_mul_canonical_mbu(c, &dy_new, &lambda[..257], &dx);
    } else {
        mod_mul_rfold_mbu(c, &dy_new, &lambda[..257], &dx);
    }
    for (g, q) in ghosts.into_iter().zip(dy_new.iter()) {
        c.resolve_ghost(g, q);
    }

    (dx, dy_new, lambda)
}

/// CANCEL the `shrunken_pz` slope: given `lambda` = `new_dy` / `new_dx` (live, 257), drive it to
/// |0> and FREE it, with `new_dx` (dx) and `new_dy` (dy) PRESERVED. Returns
/// (`new_dx`, `new_dy`). By EC linearity `new_dy/new_dx` == lambda, so this is the
/// alt-witness cleanup that removes the slope ancilla after the point coordinates
/// are computed.
///
/// Mirror of `shrunken_pz_divide_forward`, but it GHOSTS lambda (not dy) up front so only
/// `new_dy` rides through the inversion as the passenger (peak = EEA-peak + 256, same
/// as forward). After inverting `new_dx` -> cb = `1/|new_dx`|, it recomputes
/// temp = `new_dy` * cb (parity/sign corrected) = `new_dy/new_dx` == lambda's original
/// value, resolves the lambda-ghost against temp (exorcizing it), uncomputes temp
/// through the matching route-specific multiplier inverse, then reverse-inverts
/// to restore `new_dx`. On Q954, canonical new_dy[256] is released independently
/// around both EEA traversals and restored for the intervening multiply and final
/// API result.
pub fn shrunken_pz_divide_cancel(
    c: &mut Circuit,
    mut dx: Vec<QReg>,
    mut dy: Vec<QReg>,
    lambda: Vec<QReg>,
) -> (Vec<QReg>, Vec<QReg>) {
    use crate::point_add::trailmix_port::arith::compare::compare_geq_const;
    use crate::point_add::trailmix_port::arith::rfold_mbu::{
        mod_mul_canonical_mbu, mod_mul_canonical_mbu_undo, mod_mul_rfold_mbu,
        mod_mul_rfold_mbu_undo,
    };
    use crate::point_add::trailmix_port::inversion::shrunken_pz_schedule::reg_widths;
    use crate::point_add::trailmix_port::num_bigint::BigUint;
    assert_eq!(dx.len(), 257);
    assert_eq!(dy.len(), 257);
    assert_eq!(lambda.len(), 257);
    let canonical_lambda = lowq_q955_off_canonical_enabled();
    let canonical_passenger_top_off = lowq_q954_srot_counter7_enabled();
    if canonical_passenger_top_off {
        release_q954_canonical_passenger_top(c, &mut dy, "cancel-forward new_dy");
    }
    let half_bytes = vec![
        0x18, 0xfe, 0xff, 0x7f, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f, 0x00,
    ];
    let p_bytes = crate::point_add::trailmix_port::mod_arith::SECP256K1_P_LE;

    // --- sign-adjust new_dx -> |new_dx| < p/2 ---
    let reuse_sign_wire = sign_parity_q_reuse_enabled();
    let mut fused_parity = reuse_sign_wire.then(|| c.alloc_qreg("shpzcan.par_sgn"));
    let sgn = (!reuse_sign_wire).then(|| c.alloc_qreg("shpzcan.sgn"));
    let sign_control = fused_parity.as_ref().or(sgn.as_ref()).unwrap();
    compare_geq_const(c, &dx, &half_bytes, sign_control);
    controlled_field_neg(c, sign_control, &dx);

    // --- GHOST lambda (HMR each bit, free 257q) so the inversion runs lambda-free;
    // new_dy is the sole 256-bit passenger (peak = EEA-peak + 256). ---
    let mut lam_ghosts = Vec::with_capacity(lambda.len());
    for q in &lambda {
        lam_ghosts.push(c.hmr_ghost(q));
    }
    for q in lambda {
        c.zero_and_free(q);
    }

    // --- set up the inversion S_0 (B = |new_dx|, A = p, cb = 1, parity = 1) ---
    let (a0, b0, ca0, cb0, q0) = reg_widths(0);
    let (wg0, wc0) = (a0.max(b0), ca0.max(cb0));
    shrunken_pz_resize(c, &mut dx, wg0, "B");
    let mut aa = c.alloc_qreg_bits("shpzcan.A", wg0);
    let mut cca = c.alloc_qreg_bits("shpzcan.ca", wc0);
    let mut ccb = c.alloc_qreg_bits("shpzcan.cb", wc0);
    let mut qq = c.alloc_qreg_bits("shpzcan.q", q0.max(1));
    let s_rot = c.alloc_qreg_bits("shpzcan.srot", trailmix_srot_width());
    let parity = fused_parity
        .take()
        .unwrap_or_else(|| c.alloc_qreg("shpzcan.par"));
    let counter = c.alloc_qreg_bits("shpzcan.ctr", trailmix_counter_width());
    let off_owned = (!lowq_q956_off_borrow_enabled()).then(|| c.alloc_qreg("shpzcan.off"));
    let off = off_owned.as_ref().unwrap_or_else(|| {
        assert_q956_off_alias(&counter[0], &counter, &s_rot);
        &counter[0]
    });
    let load_p = |c: &mut Circuit, reg: &[QReg]| {
        for (j, q) in reg.iter().enumerate() {
            if j < 256 && (p_bytes[j / 8] >> (j % 8)) & 1 == 1 {
                c.x(q);
            }
        }
    };
    load_p(c, &aa);
    c.x(&ccb[0]);
    c.x(&parity); // parity = 1, or fused parity = 1 XOR sign

    // --- forward inversion: 1/|new_dx| in cb (passenger: new_dy) ---
    shrunken_pz_invert_forward(
        c, &mut aa, &mut dx, &mut cca, &mut ccb, &mut qq, &counter, &parity, &s_rot, off,
    );

    // --- tear down the constant pack (A=0,B=1,ca=p,q=0); keep cb=1/|new_dx| ---
    let (ta, tb, tca, tq) = (aa.len(), dx.len(), cca.len(), qq.len());
    load_p(c, &cca);
    c.x(&dx[0]);
    for q in std::mem::take(&mut aa) {
        c.zero_and_free(q);
    }
    for q in std::mem::take(&mut dx) {
        c.zero_and_free(q);
    }
    for q in std::mem::take(&mut cca) {
        c.zero_and_free(q);
    }
    for q in std::mem::take(&mut qq) {
        c.zero_and_free(q);
    }
    if canonical_passenger_top_off {
        restore_q954_canonical_passenger_top(c, &mut dy, "shpzcan.dy[256]");
    }
    // --- temp = new_dy * (1/|new_dx|), parity/sign corrected = new_dy/new_dx, the
    // original value of lambda. Resolve the lambda-ghost against it, then uncompute
    // temp. ---
    let cb_w = ccb.len();
    shrunken_pz_resize(c, &mut ccb, 257, "cb");
    let temp = c.alloc_qreg_bits("shpzcan.temp", 257);
    if canonical_lambda {
        mod_mul_canonical_mbu(c, &temp, &ccb[..257], &dy);
    } else {
        mod_mul_rfold_mbu(c, &temp, &ccb[..257], &dy);
    }
    if reuse_sign_wire {
        c.x(&parity); // fused parity -> f = NOT(sgn XOR parity)
        if canonical_lambda {
            controlled_field_neg_canonical(c, &parity, &temp);
        } else {
            controlled_field_neg(c, &parity, &temp);
        }
        for (g, q) in lam_ghosts.into_iter().zip(temp.iter()) {
            c.resolve_ghost(g, q);
        }
        if canonical_lambda {
            controlled_field_neg_canonical(c, &parity, &temp);
        } else {
            controlled_field_neg(c, &parity, &temp);
        }
        c.x(&parity);
    } else {
        let sgn = sgn.as_ref().unwrap();
        let f = c.alloc_qreg("shpzcan.negf");
        c.cx(sgn, &f);
        c.cx(&parity, &f);
        c.x(&f); // f = NOT(sgn XOR parity)
        if canonical_lambda {
            controlled_field_neg_canonical(c, &f, &temp);
        } else {
            controlled_field_neg(c, &f, &temp);
        }
        for (g, q) in lam_ghosts.into_iter().zip(temp.iter()) {
            c.resolve_ghost(g, q); // exorcize lambda (temp == lambda's value)
        }
        if canonical_lambda {
            controlled_field_neg_canonical(c, &f, &temp);
        } else {
            controlled_field_neg(c, &f, &temp);
        }
        c.x(&f);
        c.cx(&parity, &f);
        c.cx(sgn, &f); // uncompute f
        c.zero_and_free(f);
    }
    if canonical_lambda {
        mod_mul_canonical_mbu_undo(c, &temp, &ccb[..257], &dy);
    } else {
        mod_mul_rfold_mbu_undo(c, &temp, &ccb[..257], &dy);
    }
    for q in temp {
        c.zero_and_free(q);
    }
    shrunken_pz_resize(c, &mut ccb, cb_w, "cb");

    if canonical_passenger_top_off {
        release_q954_canonical_passenger_top(c, &mut dy, "cancel-reverse new_dy");
    }

    // --- re-create the pack, backward inversion (restore B=|new_dx|) ---
    aa = c.alloc_qreg_bits("shpzcan.A", ta);
    dx = c.alloc_qreg_bits("shpzcan.B", tb);
    c.x(&dx[0]);
    cca = c.alloc_qreg_bits("shpzcan.ca", tca);
    load_p(c, &cca);
    qq = c.alloc_qreg_bits("shpzcan.q", tq);
    shrunken_pz_invert_backward(
        c, &mut aa, &mut dx, &mut cca, &mut ccb, &mut qq, &counter, &parity, &s_rot, off,
    );

    // --- free the clean inversion ancillas (S_0: A=p, ca=0, cb=1, q=0) ---
    if !reuse_sign_wire {
        c.x(&parity);
    }
    c.x(&ccb[0]);
    load_p(c, &aa);
    for q in aa.into_iter().chain(cca).chain(ccb).chain(qq) {
        c.zero_and_free(q);
    }
    if let Some(off) = off_owned {
        c.zero_and_free(off);
    }
    for q in s_rot.into_iter().chain(counter) {
        c.zero_and_free(q);
    }
    if canonical_passenger_top_off {
        restore_q954_canonical_passenger_top(c, &mut dy, "shpzcan.dy[256]");
    }

    // --- un-sign-adjust: |new_dx| -> new_dx, uncompute sign state ---
    shrunken_pz_resize(c, &mut dx, 257, "dx");
    if reuse_sign_wire {
        c.x(&parity);
        controlled_field_neg(c, &parity, &dx);
        compare_geq_const(c, &dx, &half_bytes, &parity);
        c.zero_and_free(parity);
    } else {
        let sgn = sgn.unwrap();
        controlled_field_neg(c, &sgn, &dx);
        compare_geq_const(c, &dx, &half_bytes, &sgn);
        c.zero_and_free(sgn);
        c.zero_and_free(parity);
    }

    (dx, dy)
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct GatedCompareExhaustiveReport {
    pub widths_checked: usize,
    pub comparator_states_checked: usize,
    pub gate_hold_states_checked: usize,
    pub borrowed_comparator_states_checked: usize,
    pub borrowed_gate_hold_states_checked: usize,
    pub max_comparator_extra_qubits: usize,
    pub max_gate_hold_extra_qubits: usize,
    pub max_borrowed_comparator_extra_qubits: usize,
    pub max_borrowed_gate_hold_extra_qubits: usize,
}

/// Exhaustively verify the active-gated comparator and `gate_hold` skeleton for
/// widths one through five over every basis state.
#[doc(hidden)]
pub fn exhaustive_gated_compare_check() -> GatedCompareExhaustiveReport {
    use crate::circuit::{Op, OperationType};

    fn apply(ops: &[Op], mut state: u64) -> u64 {
        let bit = |state: u64, id: u64| ((state >> id) & 1) != 0;
        for op in ops {
            match op.kind {
                OperationType::X => state ^= 1u64 << op.q_target.0,
                OperationType::CX => {
                    if bit(state, op.q_control1.0) {
                        state ^= 1u64 << op.q_target.0;
                    }
                }
                OperationType::CCX => {
                    if bit(state, op.q_control1.0) && bit(state, op.q_control2.0) {
                        state ^= 1u64 << op.q_target.0;
                    }
                }
                OperationType::R => {
                    assert!(!bit(state, op.q_target.0), "freed comparator carry was not zero");
                }
                other => panic!("gated comparator emitted unexpected gate {other:?}"),
            }
        }
        state
    }

    fn word(state: u64, start: usize, width: usize) -> u64 {
        (state >> start) & ((1u64 << width) - 1)
    }

    let mut comparator_states_checked = 0usize;
    let mut gate_hold_states_checked = 0usize;
    let mut borrowed_comparator_states_checked = 0usize;
    let mut borrowed_gate_hold_states_checked = 0usize;
    let mut max_comparator_extra_qubits = 0usize;
    let mut max_gate_hold_extra_qubits = 0usize;
    let mut max_borrowed_comparator_extra_qubits = 0usize;
    let mut max_borrowed_gate_hold_extra_qubits = 0usize;

    for width in 1..=5usize {
        let mut c = Circuit::new();
        let active = c.alloc_qreg("gated-check.active");
        let out = c.alloc_qreg("gated-check.out");
        let v = c.alloc_qreg_bits("gated-check.v", width);
        let u = c.alloc_qreg_bits("gated-check.u", width);
        let vr: Vec<&QReg> = v.iter().collect();
        let ur: Vec<&QReg> = u.iter().collect();
        borrow_compare_gated_refs(&mut c, &vr, &ur, &active, &out);
        let external = 2 * width + 2;
        let builder = c.into_builder();
        let extra = builder.peak_qubits as usize - external;
        max_comparator_extra_qubits = max_comparator_extra_qubits.max(extra);
        assert_eq!(extra, 1, "width={width}: gated comparator peak changed");

        for input in 0..(1u64 << external) {
            comparator_states_checked += 1;
            let active_pre = input & 1;
            let out_pre = (input >> 1) & 1;
            let v_pre = word(input, 2, width);
            let u_pre = word(input, 2 + width, width);
            let got = apply(&builder.ops, input);
            let expected_out = out_pre ^ (active_pre & u64::from(v_pre < u_pre));
            assert_eq!(got & 1, active_pre, "width={width}: active changed");
            assert_eq!((got >> 1) & 1, expected_out, "width={width} input={input}");
            assert_eq!(word(got, 2, width), v_pre, "width={width}: v changed");
            assert_eq!(word(got, 2 + width, width), u_pre, "width={width}: u changed");
        }

        let mut c = Circuit::new();
        let active = c.alloc_qreg("borrowed-check.active");
        let out = c.alloc_qreg("borrowed-check.out");
        let carry = c.alloc_qreg("borrowed-check.carry");
        let v = c.alloc_qreg_bits("borrowed-check.v", width);
        let u = c.alloc_qreg_bits("borrowed-check.u", width);
        let vr: Vec<&QReg> = v.iter().collect();
        let ur: Vec<&QReg> = u.iter().collect();
        borrow_compare_gated_refs_with_carry(&mut c, &vr, &ur, &active, &out, &carry);
        let external = 2 * width + 3;
        let builder = c.into_builder();
        let extra = builder.peak_qubits as usize - external;
        max_borrowed_comparator_extra_qubits =
            max_borrowed_comparator_extra_qubits.max(extra);
        assert_eq!(extra, 0, "width={width}: borrowed comparator allocated");

        for input in 0..(1u64 << external) {
            if (input >> 2) & 1 != 0 {
                continue;
            }
            borrowed_comparator_states_checked += 1;
            let active_pre = input & 1;
            let out_pre = (input >> 1) & 1;
            let v_pre = word(input, 3, width);
            let u_pre = word(input, 3 + width, width);
            let got = apply(&builder.ops, input);
            let expected_out = out_pre ^ (active_pre & u64::from(v_pre < u_pre));
            assert_eq!(got & 1, active_pre, "width={width}: active changed");
            assert_eq!((got >> 1) & 1, expected_out, "width={width} input={input}");
            assert_eq!((got >> 2) & 1, 0, "width={width}: carry not restored");
            assert_eq!(word(got, 3, width), v_pre, "width={width}: v changed");
            assert_eq!(word(got, 3 + width, width), u_pre, "width={width}: u changed");
        }

        let mut c = Circuit::new();
        let active = c.alloc_qreg("gate-hold-check.active");
        let g = c.alloc_qreg("gate-hold-check.g");
        let body = c.alloc_qreg("gate-hold-check.body");
        let x = c.alloc_qreg_bits("gate-hold-check.x", width);
        let y = c.alloc_qreg_bits("gate-hold-check.y", width);
        gate_hold(&mut c, &x, &y, &active, &g, None, |c, gate| {
            c.cx(gate, &body)
        });
        let external = 2 * width + 3;
        let builder = c.into_builder();
        let extra = builder.peak_qubits as usize - external;
        max_gate_hold_extra_qubits = max_gate_hold_extra_qubits.max(extra);
        assert_eq!(extra, 1, "width={width}: gate_hold peak changed");

        for input in 0..(1u64 << external) {
            gate_hold_states_checked += 1;
            let active_pre = input & 1;
            let g_pre = (input >> 1) & 1;
            let body_pre = (input >> 2) & 1;
            let x_pre = word(input, 3, width);
            let y_pre = word(input, 3 + width, width);
            let got = apply(&builder.ops, input);
            let gate = g_pre ^ (active_pre & u64::from(x_pre < y_pre));
            assert_eq!(got & 1, active_pre, "width={width}: active changed");
            assert_eq!((got >> 1) & 1, g_pre, "width={width}: g not restored");
            assert_eq!((got >> 2) & 1, body_pre ^ gate, "width={width}: body mismatch");
            assert_eq!(word(got, 3, width), x_pre, "width={width}: x changed");
            assert_eq!(word(got, 3 + width, width), y_pre, "width={width}: y changed");
        }

        let mut c = Circuit::new();
        let active = c.alloc_qreg("borrowed-hold.active");
        let g = c.alloc_qreg("borrowed-hold.g");
        let body = c.alloc_qreg("borrowed-hold.body");
        let carry = c.alloc_qreg("borrowed-hold.carry");
        let x = c.alloc_qreg_bits("borrowed-hold.x", width);
        let y = c.alloc_qreg_bits("borrowed-hold.y", width);
        gate_hold(
            &mut c,
            &x,
            &y,
            &active,
            &g,
            Some(&carry),
            |c, gate| c.cx(gate, &body),
        );
        let external = 2 * width + 4;
        let builder = c.into_builder();
        let extra = builder.peak_qubits as usize - external;
        max_borrowed_gate_hold_extra_qubits = max_borrowed_gate_hold_extra_qubits.max(extra);
        assert_eq!(extra, 0, "width={width}: borrowed gate_hold allocated");

        for input in 0..(1u64 << external) {
            if (input >> 3) & 1 != 0 {
                continue;
            }
            borrowed_gate_hold_states_checked += 1;
            let active_pre = input & 1;
            let g_pre = (input >> 1) & 1;
            let body_pre = (input >> 2) & 1;
            let x_pre = word(input, 4, width);
            let y_pre = word(input, 4 + width, width);
            let got = apply(&builder.ops, input);
            let gate = g_pre ^ (active_pre & u64::from(x_pre < y_pre));
            assert_eq!(got & 1, active_pre, "width={width}: active changed");
            assert_eq!((got >> 1) & 1, g_pre, "width={width}: g not restored");
            assert_eq!((got >> 2) & 1, body_pre ^ gate, "width={width}: body mismatch");
            assert_eq!((got >> 3) & 1, 0, "width={width}: carry not restored");
            assert_eq!(word(got, 4, width), x_pre, "width={width}: x changed");
            assert_eq!(word(got, 4 + width, width), y_pre, "width={width}: y changed");
        }
    }

    GatedCompareExhaustiveReport {
        widths_checked: 5,
        comparator_states_checked,
        gate_hold_states_checked,
        borrowed_comparator_states_checked,
        borrowed_gate_hold_states_checked,
        max_comparator_extra_qubits,
        max_gate_hold_extra_qubits,
        max_borrowed_comparator_extra_qubits,
        max_borrowed_gate_hold_extra_qubits,
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SelectiveBorrowExhaustiveReport {
    pub bitlen_widths_checked: usize,
    pub bitlen_states_checked: usize,
    pub bitlen_parity_states_checked: usize,
    pub counter_gate_widths_checked: usize,
    pub counter_gate_states_checked: usize,
    pub done_widths_checked: usize,
    pub done_states_checked: usize,
    pub demux_widths_checked: usize,
    pub demux_states_checked: usize,
    pub swap_widths_checked: usize,
    pub swap_states_checked: usize,
    pub max_bitlen_extra_qubits: usize,
    pub max_bitlen_parity_extra_qubits: usize,
    pub max_counter_gate_extra_qubits: usize,
    pub max_done_extra_qubits: usize,
    pub max_demux_extra_qubits: usize,
    pub max_swap_extra_qubits: usize,
}

/// Exhaustively verify the actual borrowed demultiplexer and promised-support
/// swap circuits on small basis spaces. The caller must enable the sealed Q959
/// route so this exercises the production branches.
#[doc(hidden)]
pub fn exhaustive_selective_borrow_check() -> SelectiveBorrowExhaustiveReport {
    use crate::circuit::{Op, OperationType};

    assert!(lowq_q959_selective_borrow_enabled());

    fn apply(ops: &[Op], mut state: u64) -> u64 {
        let bit = |state: u64, id: u64| ((state >> id) & 1) != 0;
        for op in ops {
            match op.kind {
                OperationType::X => state ^= 1u64 << op.q_target.0,
                OperationType::CX => {
                    if bit(state, op.q_control1.0) {
                        state ^= 1u64 << op.q_target.0;
                    }
                }
                OperationType::CCX => {
                    if bit(state, op.q_control1.0) && bit(state, op.q_control2.0) {
                        state ^= 1u64 << op.q_target.0;
                    }
                }
                OperationType::R => {
                    assert!(!bit(state, op.q_target.0), "borrowed lane was not zero");
                }
                other => panic!("selective-borrow proof emitted unexpected gate {other:?}"),
            }
        }
        state
    }

    fn word(state: u64, start: usize, width: usize) -> u64 {
        (state >> start) & ((1u64 << width) - 1)
    }

    let mut demux_states_checked = 0usize;
    let mut swap_states_checked = 0usize;
    let mut bitlen_states_checked = 0usize;
    let mut bitlen_parity_states_checked = 0usize;
    let mut counter_gate_states_checked = 0usize;
    let mut done_states_checked = 0usize;
    let mut max_bitlen_extra_qubits = 0usize;
    let mut max_bitlen_parity_extra_qubits = 0usize;
    let mut max_counter_gate_extra_qubits = 0usize;
    let mut max_done_extra_qubits = 0usize;
    let mut max_demux_extra_qubits = 0usize;
    let mut max_swap_extra_qubits = 0usize;

    for width in 1..=3usize {
        for lo_a in 0..width {
            for lo_b in 0..width {
                for subtract_diff in [false, true] {
                    let mut c = Circuit::new();
                    let active = c.alloc_qreg("bitlen-check.active");
                    let gate = c.alloc_qreg("bitlen-check.gate");
                    let target = c.alloc_qreg_bits("bitlen-check.target", 3);
                    let target_refs: Vec<&QReg> = target.iter().collect();
                    let a = c.alloc_qreg_bits("bitlen-check.a", width);
                    let b = c.alloc_qreg_bits("bitlen-check.b", width);
                    let extra = c.alloc_qreg_bits("bitlen-check.extra", 2);
                    let extra_refs: Vec<&QReg> = extra.iter().collect();
                    direct_bitlen_diff_update(
                        &mut c,
                        &a,
                        &b,
                        lo_a,
                        lo_b,
                        &target_refs,
                        &active,
                        &gate,
                        subtract_diff,
                        &extra_refs,
                    );
                    let external = 7 + 2 * width;
                    let builder = c.into_builder();
                    let added = builder.peak_qubits as usize - external;
                    max_bitlen_extra_qubits = max_bitlen_extra_qubits.max(added);
                    assert_eq!(
                        added, 0,
                        "width={width} lo_a={lo_a} lo_b={lo_b}: direct bitlen allocated"
                    );

                    let target_start = 2;
                    let a_start = target_start + 3;
                    let b_start = a_start + width;
                    let target_mask = 0b111u64;
                    for input in 0..(1u64 << external) {
                        if (input >> 1) & 1 != 0 {
                            continue;
                        }
                        let a_pre = word(input, a_start, width);
                        let b_pre = word(input, b_start, width);
                        if (a_pre >> lo_a) == 0 || (b_pre >> lo_b) == 0 {
                            continue;
                        }
                        bitlen_states_checked += 1;
                        let active_pre = input & 1;
                        let target_pre = word(input, target_start, 3);
                        let a_len = 64 - a_pre.leading_zeros() as u64;
                        let b_len = 64 - b_pre.leading_zeros() as u64;
                        let delta = if subtract_diff {
                            b_len.wrapping_sub(a_len)
                        } else {
                            a_len.wrapping_sub(b_len)
                        };
                        let target_post =
                            target_pre.wrapping_add(active_pre * delta) & target_mask;
                        let expected = (input & !(target_mask << target_start))
                            | (target_post << target_start);
                        let got = apply(&builder.ops, input);
                        assert_eq!(
                            got, expected,
                            "bitlen width={width} lo_a={lo_a} lo_b={lo_b} input={input}"
                        );
                    }
                }
            }
        }
    }

    for width in 1..=3usize {
        for lo_a in 0..width {
            for lo_b in 0..width {
                let mut c = Circuit::new();
                let active = c.alloc_qreg("bitlen-parity-check.active");
                let out = c.alloc_qreg("bitlen-parity-check.out");
                let a = c.alloc_qreg_bits("bitlen-parity-check.a", width);
                let b = c.alloc_qreg_bits("bitlen-parity-check.b", width);
                let extra = c.alloc_qreg_bits("bitlen-parity-check.extra", 2);
                let extra_refs: Vec<&QReg> = extra.iter().collect();
                direct_bitlen_diff_parity(
                    &mut c,
                    &a,
                    &b,
                    lo_a,
                    lo_b,
                    &out,
                    &active,
                    &extra_refs,
                );
                let external = 4 + 2 * width;
                let builder = c.into_builder();
                let added = builder.peak_qubits as usize - external;
                max_bitlen_parity_extra_qubits =
                    max_bitlen_parity_extra_qubits.max(added);
                assert_eq!(
                    added, 0,
                    "width={width} lo_a={lo_a} lo_b={lo_b}: direct parity allocated"
                );

                let a_start = 2;
                let b_start = a_start + width;
                for input in 0..(1u64 << external) {
                    let active_pre = input & 1;
                    let a_pre = word(input, a_start, width);
                    let b_pre = word(input, b_start, width);
                    if active_pre == 1 && ((a_pre >> lo_a) == 0 || (b_pre >> lo_b) == 0) {
                        continue;
                    }
                    bitlen_parity_states_checked += 1;
                    let a_len = 64 - a_pre.leading_zeros() as u64;
                    let b_len = 64 - b_pre.leading_zeros() as u64;
                    let expected = input ^ (active_pre * ((a_len ^ b_len) & 1) << 1);
                    let got = apply(&builder.ops, input);
                    assert_eq!(
                        got, expected,
                        "bitlen parity width={width} lo_a={lo_a} lo_b={lo_b} input={input}"
                    );
                }
            }
        }
    }

    for width in 1..=3usize {
        let mut c = Circuit::new();
        let parity = c.alloc_qreg("counter-gate.parity");
        let g = c.alloc_qreg("counter-gate.g");
        let s_rot = c.alloc_qreg_bits("counter-gate.s", 2);
        let body = c.alloc_qreg("counter-gate.body");
        let x = c.alloc_qreg_bits("counter-gate.x", width);
        let y = c.alloc_qreg_bits("counter-gate.y", width);
        let counter = c.alloc_qreg_bits("counter-gate.counter", 2);
        let extra = c.alloc_qreg_bits("counter-gate.extra", 2);
        let candidates: Vec<&QReg> = x
            .iter()
            .chain(y.iter())
            .chain(counter.iter())
            .chain(extra.iter())
            .chain(std::iter::once(&body))
            .chain(std::iter::once(&parity))
            .chain(s_rot.iter())
            .chain(std::iter::once(&g))
            .collect();
        gate_hold_counter_zero(
            &mut c,
            &x,
            &y,
            &counter,
            &parity,
            &s_rot,
            &g,
            &candidates,
            |c, gate| c.cx(gate, &body),
        );
        let external = 9 + 2 * width;
        let builder = c.into_builder();
        let added = builder.peak_qubits as usize - external;
        max_counter_gate_extra_qubits = max_counter_gate_extra_qubits.max(added);
        assert_eq!(added, 0, "width={width}: counter gate allocated");

        let x_start = 5;
        let y_start = x_start + width;
        let counter_start = y_start + width;
        for input in 0..(1u64 << external) {
            if input & 0b1110 != 0 {
                continue;
            }
            counter_gate_states_checked += 1;
            let x_pre = word(input, x_start, width);
            let y_pre = word(input, y_start, width);
            let counter_pre = word(input, counter_start, 2);
            let gate_pre = u64::from(counter_pre == 0 && x_pre < y_pre);
            let expected = input ^ (gate_pre << 4);
            let got = apply(&builder.ops, input);
            assert_eq!(got, expected, "counter gate width={width} input={input}");
        }
    }

    for width in 1..=3usize {
        for inverse in [false, true] {
            let mut c = Circuit::new();
            let off = c.alloc_qreg("done-check.off");
            let s_rot = c.alloc_qreg_bits("done-check.s", 2);
            let aa = c.alloc_qreg_bits("done-check.a", width);
            let qq = c.alloc_qreg_bits("done-check.q", 2);
            let counter = c.alloc_qreg_bits("done-check.counter", 2);
            let extra = c.alloc_qreg_bits("done-check.extra", 2);
            let candidates: Vec<&QReg> = aa
                .iter()
                .chain(qq.iter())
                .chain(counter.iter())
                .chain(extra.iter())
                .chain(s_rot.iter())
                .chain(std::iter::once(&off))
                .collect();
            done_counter_fn(
                &mut c,
                &aa,
                &qq,
                &counter,
                &s_rot,
                &off,
                &candidates,
                inverse,
            );
            let external = 9 + width;
            let builder = c.into_builder();
            let added = builder.peak_qubits as usize - external;
            max_done_extra_qubits = max_done_extra_qubits.max(added);
            assert_eq!(added, 0, "width={width}: done counter allocated");

            let a_start = 3;
            let q_start = a_start + width;
            let counter_start = q_start + 2;
            for input in 0..(1u64 << external) {
                if input & 0b111 != 0 {
                    continue;
                }
                let a_pre = word(input, a_start, width);
                let q_pre = word(input, q_start, 2);
                let counter_pre = word(input, counter_start, 2);
                let conv = a_pre == 0 && q_pre == 0;
                if inverse {
                    if (counter_pre == 0) != !conv {
                        continue;
                    }
                } else if (counter_pre > 0 && !conv) || (counter_pre == 3 && conv) {
                    continue;
                }
                done_states_checked += 1;
                let counter_post = if inverse {
                    counter_pre.saturating_sub(u64::from(counter_pre > 0))
                } else {
                    counter_pre + u64::from(conv)
                };
                let mask = 0b11u64 << counter_start;
                let expected = (input & !mask) | (counter_post << counter_start);
                let got = apply(&builder.ops, input);
                assert_eq!(got, expected, "done width={width} inverse={inverse} input={input}");
            }
        }
    }

    for width in 1..=3usize {
        let mut c = Circuit::new();
        let active = c.alloc_qreg("demux-check.active");
        let gate = c.alloc_qreg("demux-check.gate");
        let s = c.alloc_qreg_bits("demux-check.s", width);
        let q = c.alloc_qreg_bits("demux-check.q", 1usize << width);
        let lender_regs = c.alloc_qreg_bits("demux-check.lender", width.saturating_sub(1));
        let lenders: Vec<&QReg> = lender_regs.iter().collect();
        let s_refs: Vec<&QReg> = s.iter().collect();
        set_bit_at_s_gated(&mut c, &q, &s_refs, &active, &gate, &lenders);
        let external = 2 + width + (1usize << width) + lender_regs.len();
        let builder = c.into_builder();
        let extra = builder.peak_qubits as usize - external;
        max_demux_extra_qubits = max_demux_extra_qubits.max(extra);
        assert_eq!(extra, 0, "width={width}: borrowed demux allocated");

        let s_start = 2;
        let q_start = s_start + width;
        for input in 0..(1u64 << external) {
            if (input >> 1) & 1 != 0 {
                continue;
            }
            demux_states_checked += 1;
            let active_pre = input & 1;
            let selected = word(input, s_start, width) as usize;
            let expected = input ^ (active_pre << (q_start + selected));
            let got = apply(&builder.ops, input);
            assert_eq!(got, expected, "demux width={width} input={input}");
        }
    }

    for width in 1..=2usize {
        let mut c = Circuit::new();
        let _active = c.alloc_qreg("swap-check.spectator");
        let parity = c.alloc_qreg("swap-check.parity");
        let off = c.alloc_qreg("swap-check.off");
        let s_rot = c.alloc_qreg_bits("swap-check.s", 2);
        let aa = c.alloc_qreg_bits("swap-check.a", width);
        let bb = c.alloc_qreg_bits("swap-check.b", width);
        let cca = c.alloc_qreg_bits("swap-check.ca", width);
        let ccb = c.alloc_qreg_bits("swap-check.cb", width);
        let qq = c.alloc_qreg_bits("swap-check.q", width);
        let counter = c.alloc_qreg_bits("swap-check.counter", 1);
        borrowed_swap_in_place(
            &mut c, &aa, &bb, &cca, &ccb, &qq, &counter, &parity, &s_rot, &off,
        );
        let external = 6 + 5 * width;
        let builder = c.into_builder();
        let extra = builder.peak_qubits as usize - external;
        max_swap_extra_qubits = max_swap_extra_qubits.max(extra);
        assert_eq!(extra, 0, "width={width}: borrowed swap allocated");

        let a_start = 5;
        let b_start = a_start + width;
        let ca_start = b_start + width;
        let cb_start = ca_start + width;
        let q_start = cb_start + width;
        let counter_start = q_start + width;
        let mask = (1u64 << width) - 1;
        for input in 0..(1u64 << external) {
            if input & 0b1_1100 != 0 {
                continue;
            }
            let a_pre = word(input, a_start, width);
            let b_pre = word(input, b_start, width);
            let q_pre = word(input, q_start, width);
            let counter_pre = word(input, counter_start, 1);
            let gate = counter_pre == 0 && q_pre == 0 && a_pre != 0;
            if gate && b_pre == 0 {
                continue;
            }
            swap_states_checked += 1;
            let mut expected = input;
            if gate {
                let ca_pre = word(input, ca_start, width);
                let cb_pre = word(input, cb_start, width);
                expected &= !((mask << a_start)
                    | (mask << b_start)
                    | (mask << ca_start)
                    | (mask << cb_start));
                expected |= b_pre << a_start;
                expected |= a_pre << b_start;
                expected |= cb_pre << ca_start;
                expected |= ca_pre << cb_start;
                expected ^= 1 << 1;
            }
            let got = apply(&builder.ops, input);
            assert_eq!(got, expected, "swap width={width} input={input}");
        }
    }

    SelectiveBorrowExhaustiveReport {
        bitlen_widths_checked: 3,
        bitlen_states_checked,
        bitlen_parity_states_checked,
        counter_gate_widths_checked: 3,
        counter_gate_states_checked,
        done_widths_checked: 3,
        done_states_checked,
        demux_widths_checked: 3,
        demux_states_checked,
        swap_widths_checked: 2,
        swap_states_checked,
        max_bitlen_extra_qubits,
        max_bitlen_parity_extra_qubits,
        max_counter_gate_extra_qubits,
        max_done_extra_qubits,
        max_demux_extra_qubits,
        max_swap_extra_qubits,
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct OffBorrowExhaustiveReport {
    pub body_widths_checked: usize,
    pub forward_states_checked: usize,
    pub reverse_states_checked: usize,
    pub roundtrip_states_checked: usize,
    pub composed_widths_checked: usize,
    pub composed_states_checked: usize,
    pub composed_roundtrip_states_checked: usize,
    pub demux_widths_checked: usize,
    pub demux_states_checked: usize,
    pub done_widths_checked: usize,
    pub done_states_checked: usize,
    pub swap_widths_checked: usize,
    pub swap_states_checked: usize,
    pub max_body_extra_qubits: usize,
    pub max_composed_extra_qubits: usize,
    pub max_demux_extra_qubits: usize,
    pub max_done_extra_qubits: usize,
    pub max_swap_extra_qubits: usize,
}

/// Exhaustively verify the Q956 support contract on small basis spaces.
/// Active branches require the borrowed counter lane to be zero; inactive
/// branches deliberately cover both lane values and must remain unchanged.
/// The body checks exercise the production masked shift/increment primitives
/// in both directions and as a complete forward/reverse cleanup pair.
#[doc(hidden)]
pub fn exhaustive_off_borrow_check() -> OffBorrowExhaustiveReport {
    use crate::circuit::{Op, OperationType};

    assert!(lowq_q956_off_borrow_enabled());

    fn apply(ops: &[Op], mut state: u64) -> u64 {
        let bit = |state: u64, id: u64| ((state >> id) & 1) != 0;
        for op in ops {
            match op.kind {
                OperationType::X => state ^= 1u64 << op.q_target.0,
                OperationType::CX => {
                    if bit(state, op.q_control1.0) {
                        state ^= 1u64 << op.q_target.0;
                    }
                }
                OperationType::CCX => {
                    if bit(state, op.q_control1.0) && bit(state, op.q_control2.0) {
                        state ^= 1u64 << op.q_target.0;
                    }
                }
                OperationType::R => {
                    assert!(!bit(state, op.q_target.0), "Q956 proof freed a nonzero lane");
                }
                other => panic!("Q956 off-borrow proof emitted unexpected gate {other:?}"),
            }
        }
        state
    }

    fn word(state: u64, start: usize, width: usize) -> u64 {
        (state >> start) & ((1u64 << width) - 1)
    }

    let mut forward_states_checked = 0usize;
    let mut reverse_states_checked = 0usize;
    let mut roundtrip_states_checked = 0usize;
    let mut composed_states_checked = 0usize;
    let mut composed_roundtrip_states_checked = 0usize;
    let mut demux_states_checked = 0usize;
    let mut done_states_checked = 0usize;
    let mut swap_states_checked = 0usize;
    let mut max_body_extra_qubits = 0usize;
    let mut max_composed_extra_qubits = 0usize;
    let mut max_demux_extra_qubits = 0usize;
    let mut max_done_extra_qubits = 0usize;
    let mut max_swap_extra_qubits = 0usize;

    for width in 2..=4usize {
        let build = |forward: bool, roundtrip: bool| {
            let mut c = Circuit::new();
            let active = c.alloc_qreg("off-body.active");
            let predicate = c.alloc_qreg("off-body.predicate");
            let off = c.alloc_qreg("off-body.borrowed");
            let s = c.alloc_qreg_bits("off-body.s", 3);
            let value = c.alloc_qreg_bits("off-body.value", width);
            let lender_regs = c.alloc_qreg_bits("off-body.lender", 4);
            let s_refs: Vec<&QReg> = s.iter().collect();
            let candidates: Vec<&QReg> = value
                .iter()
                .chain(s.iter())
                .chain(lender_regs.iter())
                .chain(std::iter::once(&active))
                .chain(std::iter::once(&predicate))
                .chain(std::iter::once(&off))
                .collect();
            let emit_forward = |c: &mut Circuit| {
                c.ccx(&active, &predicate, &off);
                rotate_one_by_off(c, &value, &active, &off, true, &candidates);
                ctrl_inc_by_off(c, &active, &off, &s_refs, &candidates);
                c.ccx(&active, &predicate, &off);
            };
            let emit_reverse = |c: &mut Circuit| {
                c.ccx(&active, &predicate, &off);
                ctrl_dec_by_off(c, &active, &off, &s_refs, &candidates);
                rotate_one_by_off(c, &value, &active, &off, false, &candidates);
                c.ccx(&active, &predicate, &off);
            };
            if forward {
                emit_forward(&mut c);
                if roundtrip {
                    emit_reverse(&mut c);
                }
            } else {
                emit_reverse(&mut c);
            }
            c.into_builder()
        };

        let external = 10 + width;
        let forward = build(true, false);
        let reverse = build(false, false);
        let roundtrip = build(true, true);
        for builder in [&forward, &reverse, &roundtrip] {
            let extra = builder.peak_qubits as usize - external;
            max_body_extra_qubits = max_body_extra_qubits.max(extra);
            assert_eq!(extra, 0, "width={width}: Q956 body allocated a lane");
        }

        let s_start = 3;
        let value_start = 6;
        let s_mask = 0b111u64;
        let value_mask = (1u64 << width) - 1;
        for input in 0..(1u64 << external) {
            let active = input & 1;
            let predicate = (input >> 1) & 1;
            let off = (input >> 2) & 1;
            if active == 1 && off != 0 {
                continue;
            }
            let gate = active & predicate;
            let s_pre = word(input, s_start, 3);
            let value_pre = word(input, value_start, width);

            if gate == 0 || value_pre >> (width - 1) == 0 {
                forward_states_checked += 1;
                let s_post = s_pre.wrapping_add(gate) & s_mask;
                let value_post = (value_pre << gate) & value_mask;
                let expected = (input
                    & !((s_mask << s_start) | (value_mask << value_start)))
                    | (s_post << s_start)
                    | (value_post << value_start);
                assert_eq!(
                    apply(&forward.ops, input),
                    expected,
                    "Q956 forward width={width} input={input}"
                );
                roundtrip_states_checked += 1;
                assert_eq!(
                    apply(&roundtrip.ops, input),
                    input,
                    "Q956 roundtrip width={width} input={input}"
                );
            }

            if gate == 0 || value_pre & 1 == 0 {
                reverse_states_checked += 1;
                let s_post = s_pre.wrapping_sub(gate) & s_mask;
                let value_post = value_pre >> gate;
                let expected = (input
                    & !((s_mask << s_start) | (value_mask << value_start)))
                    | (s_post << s_start)
                    | (value_post << value_start);
                assert_eq!(
                    apply(&reverse.ops, input),
                    expected,
                    "Q956 reverse width={width} input={input}"
                );
            }
        }
    }

    for width in 1..=2usize {
        let build = |roundtrip: bool| {
            let mut c = Circuit::new();
            let parity = c.alloc_qreg("off-composed.parity");
            let g = c.alloc_qreg("off-composed.g");
            let s_rot = c.alloc_qreg_bits("off-composed.srot", 2);
            let predicate = c.alloc_qreg("off-composed.predicate");
            let target = c.alloc_qreg_bits("off-composed.target", 2);
            let value = c.alloc_qreg_bits("off-composed.value", 2);
            let x = c.alloc_qreg_bits("off-composed.x", width);
            let y = c.alloc_qreg_bits("off-composed.y", width);
            let counter = c.alloc_qreg_bits("off-composed.counter", 2);
            let target_refs: Vec<&QReg> = target.iter().collect();
            let candidates: Vec<&QReg> = x
                .iter()
                .chain(y.iter())
                .chain(counter.iter())
                .chain(value.iter())
                .chain(target.iter())
                .chain(std::iter::once(&predicate))
                .chain(std::iter::once(&parity))
                .chain(s_rot.iter())
                .chain(std::iter::once(&g))
                .collect();
            let emit = |c: &mut Circuit, forward: bool| {
                gate_hold_counter_zero(
                    c,
                    &x,
                    &y,
                    &counter,
                    &parity,
                    &s_rot,
                    &g,
                    &candidates,
                    |c, active| {
                        let off = &counter[0];
                        c.ccx(active, &predicate, off);
                        if forward {
                            rotate_one_by_off(c, &value, active, off, true, &candidates);
                            ctrl_inc_by_off(c, active, off, &target_refs, &candidates);
                        } else {
                            ctrl_dec_by_off(c, active, off, &target_refs, &candidates);
                            rotate_one_by_off(c, &value, active, off, false, &candidates);
                        }
                        c.ccx(active, &predicate, off);
                    },
                );
            };
            emit(&mut c, true);
            if roundtrip {
                emit(&mut c, false);
            }
            c.into_builder()
        };

        let external = 11 + 2 * width;
        let forward = build(false);
        let roundtrip = build(true);
        for builder in [&forward, &roundtrip] {
            let extra = builder.peak_qubits as usize - external;
            max_composed_extra_qubits = max_composed_extra_qubits.max(extra);
            assert_eq!(
                extra, 0,
                "width={width}: composed Q956 gate holder allocated a lane"
            );
        }

        let target_start = 5;
        let value_start = 7;
        let x_start = 9;
        let y_start = x_start + width;
        let counter_start = y_start + width;
        for input in 0..(1u64 << external) {
            if input & 0b1110 != 0 {
                continue;
            }
            let predicate = (input >> 4) & 1;
            let target_pre = word(input, target_start, 2);
            let value_pre = word(input, value_start, 2);
            let x_pre = word(input, x_start, width);
            let y_pre = word(input, y_start, width);
            let counter_pre = word(input, counter_start, 2);
            let gate = u64::from(counter_pre == 0 && x_pre < y_pre) * predicate;
            if gate == 1 && value_pre >> 1 != 0 {
                continue;
            }
            composed_states_checked += 1;
            let target_post = target_pre.wrapping_add(gate) & 0b11;
            let value_post = (value_pre << gate) & 0b11;
            let expected = (input & !((0b11 << target_start) | (0b11 << value_start)))
                | (target_post << target_start)
                | (value_post << value_start);
            assert_eq!(
                apply(&forward.ops, input),
                expected,
                "Q956 composed width={width} input={input}"
            );
            composed_roundtrip_states_checked += 1;
            assert_eq!(
                apply(&roundtrip.ops, input),
                input,
                "Q956 composed roundtrip width={width} input={input}"
            );
        }
    }

    for width in 1..=3usize {
        let mut c = Circuit::new();
        let active = c.alloc_qreg("off-demux.active");
        let off = c.alloc_qreg("off-demux.borrowed");
        let s = c.alloc_qreg_bits("off-demux.s", width);
        let q = c.alloc_qreg_bits("off-demux.q", 1usize << width);
        let lender_regs = c.alloc_qreg_bits("off-demux.lender", width);
        let lenders: Vec<&QReg> = lender_regs.iter().collect();
        let s_refs: Vec<&QReg> = s.iter().collect();
        set_bit_at_s_gated(&mut c, &q, &s_refs, &active, &off, &lenders);
        let external = 2 + 2 * width + (1usize << width);
        let builder = c.into_builder();
        let extra = builder.peak_qubits as usize - external;
        max_demux_extra_qubits = max_demux_extra_qubits.max(extra);
        assert_eq!(extra, 0, "width={width}: Q956 demux allocated a lane");

        let s_start = 2;
        let q_start = s_start + width;
        for input in 0..(1u64 << external) {
            let active_pre = input & 1;
            let off_pre = (input >> 1) & 1;
            if active_pre == 1 && off_pre != 0 {
                continue;
            }
            demux_states_checked += 1;
            let selected = word(input, s_start, width) as usize;
            let expected = input ^ (active_pre << (q_start + selected));
            assert_eq!(
                apply(&builder.ops, input),
                expected,
                "Q956 demux width={width} input={input}"
            );
        }
    }

    for width in 1..=2usize {
        for inverse in [false, true] {
            let mut c = Circuit::new();
            let s_rot = c.alloc_qreg_bits("off-done.s", 3);
            let aa = c.alloc_qreg_bits("off-done.a", width);
            let qq = c.alloc_qreg_bits("off-done.q", 2);
            let counter = c.alloc_qreg_bits("off-done.counter", 2);
            let extra = c.alloc_qreg_bits("off-done.extra", 3);
            let off = &counter[0];
            let candidates: Vec<&QReg> = aa
                .iter()
                .chain(qq.iter())
                .chain(counter.iter())
                .chain(extra.iter())
                .chain(s_rot.iter())
                .collect();
            done_counter_fn(
                &mut c,
                &aa,
                &qq,
                &counter,
                &s_rot,
                off,
                &candidates,
                inverse,
            );
            let external = 10 + width;
            let builder = c.into_builder();
            let extra = builder.peak_qubits as usize - external;
            max_done_extra_qubits = max_done_extra_qubits.max(extra);
            assert_eq!(extra, 0, "width={width}: Q956 done allocated a lane");

            let a_start = 3;
            let q_start = a_start + width;
            let counter_start = q_start + 2;
            for input in 0..(1u64 << external) {
                if input & 0b111 != 0 {
                    continue;
                }
                let a_pre = word(input, a_start, width);
                let q_pre = word(input, q_start, 2);
                let counter_pre = word(input, counter_start, 2);
                let conv = a_pre == 0 && q_pre == 0;
                if inverse {
                    if (counter_pre == 0) != !conv {
                        continue;
                    }
                } else if (counter_pre > 0 && !conv) || (counter_pre == 3 && conv) {
                    continue;
                }
                done_states_checked += 1;
                let counter_post = if inverse {
                    counter_pre.saturating_sub(u64::from(counter_pre > 0))
                } else {
                    counter_pre + u64::from(conv)
                };
                let mask = 0b11u64 << counter_start;
                let expected = (input & !mask) | (counter_post << counter_start);
                assert_eq!(
                    apply(&builder.ops, input),
                    expected,
                    "Q956 done width={width} inverse={inverse} input={input}"
                );
            }
        }
    }

    for width in 1..=2usize {
        let mut c = Circuit::new();
        let parity = c.alloc_qreg("off-swap.parity");
        let s_rot = c.alloc_qreg_bits("off-swap.s", 3);
        let aa = c.alloc_qreg_bits("off-swap.a", width);
        let bb = c.alloc_qreg_bits("off-swap.b", width);
        let cca = c.alloc_qreg_bits("off-swap.ca", width);
        let ccb = c.alloc_qreg_bits("off-swap.cb", width);
        let qq = c.alloc_qreg_bits("off-swap.q", width);
        let counter = c.alloc_qreg_bits("off-swap.counter", 1);
        let off = &counter[0];
        borrowed_swap_in_place(
            &mut c, &aa, &bb, &cca, &ccb, &qq, &counter, &parity, &s_rot, off,
        );
        let external = 5 + 5 * width;
        let builder = c.into_builder();
        let extra = builder.peak_qubits as usize - external;
        max_swap_extra_qubits = max_swap_extra_qubits.max(extra);
        assert_eq!(extra, 0, "width={width}: Q956 swap allocated a lane");

        let a_start = 4;
        let b_start = a_start + width;
        let ca_start = b_start + width;
        let cb_start = ca_start + width;
        let q_start = cb_start + width;
        let counter_start = q_start + width;
        let mask = (1u64 << width) - 1;
        for input in 0..(1u64 << external) {
            if input & 0b1110 != 0 {
                continue;
            }
            let a_pre = word(input, a_start, width);
            let b_pre = word(input, b_start, width);
            let q_pre = word(input, q_start, width);
            let counter_pre = word(input, counter_start, 1);
            let gate = counter_pre == 0 && q_pre == 0 && a_pre != 0;
            if gate && b_pre == 0 {
                continue;
            }
            swap_states_checked += 1;
            let mut expected = input;
            if gate {
                let ca_pre = word(input, ca_start, width);
                let cb_pre = word(input, cb_start, width);
                expected &= !((mask << a_start)
                    | (mask << b_start)
                    | (mask << ca_start)
                    | (mask << cb_start));
                expected |= b_pre << a_start;
                expected |= a_pre << b_start;
                expected |= cb_pre << ca_start;
                expected |= ca_pre << cb_start;
                expected ^= 1;
            }
            assert_eq!(
                apply(&builder.ops, input),
                expected,
                "Q956 swap width={width} input={input}"
            );
        }
    }

    OffBorrowExhaustiveReport {
        body_widths_checked: 3,
        forward_states_checked,
        reverse_states_checked,
        roundtrip_states_checked,
        composed_widths_checked: 2,
        composed_states_checked,
        composed_roundtrip_states_checked,
        demux_widths_checked: 3,
        demux_states_checked,
        done_widths_checked: 2,
        done_states_checked,
        swap_widths_checked: 2,
        swap_states_checked,
        max_body_extra_qubits,
        max_composed_extra_qubits,
        max_demux_extra_qubits,
        max_done_extra_qubits,
        max_swap_extra_qubits,
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Q954SrotCounter7Report {
    pub high_lane_differential_cases_checked: usize,
    pub body_forward_cases_checked: usize,
    pub body_reverse_cases_checked: usize,
    pub body_roundtrip_cases_checked: usize,
    pub late_inactive_cases_checked: usize,
    pub passenger_cases_checked: usize,
    pub passenger_release_intervals_checked: usize,
    pub srot_qubits_saved: usize,
    pub passenger_qubits_saved: usize,
    pub max_helper_toffoli_increase: usize,
}

/// Differentially check the borrowed high shift lane against a five-owned-lane
/// reference, including production division/multiply bodies and the late
/// inactive branch where counter[7] is one. Also exercise all canonical
/// passenger-top release intervals used by the Q954 route.
#[doc(hidden)]
pub fn q954_srot_counter7_roundtrip_check() -> Q954SrotCounter7Report {
    use crate::circuit::{OperationType, QubitId};
    use crate::point_add::B;
    use crate::sim::Simulator;
    use sha3::{
        digest::{ExtendableOutput, Update},
        Shake128,
    };

    assert!(lowq_q954_srot_counter7_enabled());
    let certificate = q954_srot_counter7_schedule_certificate();
    assert_eq!(certificate.first_terminal_capable_row, 371);
    assert_eq!(certificate.last_ctz_bit4_row, 477);
    assert_eq!(certificate.last_raw_bit4_barrel_row, 495);
    assert_eq!(certificate.max_pre_body_counter, 124);
    assert_eq!(certificate.max_final_counter, 159);

    struct Harness {
        builder: B,
        registers: Vec<Vec<u32>>,
        external: Vec<u32>,
    }

    #[derive(Debug, Eq, PartialEq)]
    struct Snapshot {
        registers: Vec<u64>,
        phase: u64,
        internal_clean: bool,
    }

    fn ids(reg: &[QReg]) -> Vec<u32> {
        reg.iter().map(QReg::id).collect()
    }

    fn external_ids(registers: &[Vec<u32>]) -> Vec<u32> {
        let mut out = Vec::new();
        for &id in registers.iter().flatten() {
            if !out.contains(&id) {
                out.push(id);
            }
        }
        out
    }

    fn simulate(harness: &Harness, values: &[u64]) -> Snapshot {
        assert_eq!(harness.registers.len(), values.len());
        let mut seed = Shake128::default();
        seed.update(b"q954-srot-counter7-differential");
        let mut xof = seed.finalize_xof();
        let mut sim = Simulator::new(
            harness.builder.next_qubit as usize,
            harness.builder.next_bit as usize,
            &mut xof,
        );
        for (register, &value) in harness.registers.iter().zip(values) {
            assert!(register.len() <= 64);
            for (bit, &id) in register.iter().enumerate() {
                if (value >> bit) & 1 == 1 {
                    *sim.qubit_mut(QubitId(u64::from(id))) |= 1;
                }
            }
        }
        sim.apply_iter(harness.builder.ops.iter());
        let registers = harness
            .registers
            .iter()
            .map(|register| {
                register.iter().enumerate().fold(0u64, |value, (bit, &id)| {
                    value | ((sim.qubit(QubitId(u64::from(id))) & 1) << bit)
                })
            })
            .collect();
        let internal_clean = (0..harness.builder.next_qubit).all(|id| {
            harness.external.contains(&id)
                || sim.qubit(QubitId(u64::from(id))) & 1 == 0
        });
        Snapshot {
            registers,
            phase: sim.phase & 1,
            internal_clean,
        }
    }

    fn toffoli(builder: &B) -> usize {
        builder
            .ops
            .iter()
            .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
            .count()
    }

    // mode: 0=forward, 1=reverse, 2=forward+reverse.
    fn build_body(split: bool, multiply: bool, mode: u8) -> Harness {
        let mut c = Circuit::new();
        let a = c.alloc_qreg_bits("q954-body.a", 20);
        let b = c.alloc_qreg_bits("q954-body.b", 20);
        let q = c.alloc_qreg_bits("q954-body.q", 17);
        let owned = c.alloc_qreg_bits("q954-body.srot", if split { 4 } else { 5 });
        let counter = c.alloc_qreg_bits("q954-body.counter", 8);
        let active = c.alloc_qreg("q954-body.active");
        let lender_regs = c.alloc_qreg_bits("q954-body.lenders", 20);
        let lenders: Vec<&QReg> = lender_regs.iter().collect();
        let s_rot: Vec<&QReg> = if split {
            vec![
                &owned[0],
                &owned[1],
                &owned[2],
                &owned[3],
                &counter[7],
            ]
        } else {
            owned.iter().collect()
        };
        let off = &counter[0];
        let emit = |c: &mut Circuit, inverse: bool| {
            if multiply {
                if inverse {
                    multiply_substep_windowed_inv(
                        c, &a, &b, &q, &s_rot, off, &active, &lenders, 0, 0, 5, 5,
                    );
                } else {
                    multiply_substep_windowed(
                        c, &a, &b, &q, &s_rot, off, &active, &lenders, 0, 0, 5, 5,
                    );
                }
            } else if inverse {
                division_substep_windowed_inv(
                    c, &a, &b, &q, &s_rot, off, &active, &lenders, 0, 0, 5, 5,
                );
            } else {
                division_substep_windowed(
                    c, &a, &b, &q, &s_rot, off, &active, &lenders, 0, 0, 5, 5,
                );
            }
        };
        match mode {
            0 => emit(&mut c, false),
            1 => emit(&mut c, true),
            2 => {
                emit(&mut c, false);
                emit(&mut c, true);
            }
            _ => unreachable!(),
        }

        let mut s_ids = ids(&owned);
        if split {
            s_ids.push(counter[7].id());
        }
        let registers = vec![
            ids(&a),
            ids(&b),
            ids(&q),
            s_ids,
            ids(&counter),
            vec![active.id()],
            ids(&lender_regs),
        ];
        let external = external_ids(&registers);
        Harness {
            builder: c.into_builder(),
            registers,
            external,
        }
    }

    fn build_high_carry(split: bool, mode: u8) -> Harness {
        let mut c = Circuit::new();
        let owned = c.alloc_qreg_bits("q954-carry.srot", if split { 4 } else { 5 });
        let counter = c.alloc_qreg_bits("q954-carry.counter", 8);
        let active = c.alloc_qreg("q954-carry.active");
        let lenders = c.alloc_qreg_bits("q954-carry.lenders", 8);
        let candidates: Vec<&QReg> = lenders
            .iter()
            .chain(owned.iter())
            .chain(counter.iter())
            .chain(std::iter::once(&active))
            .collect();
        let s_rot: Vec<&QReg> = if split {
            vec![
                &owned[0],
                &owned[1],
                &owned[2],
                &owned[3],
                &counter[7],
            ]
        } else {
            owned.iter().collect()
        };
        let off = &counter[0];
        match mode {
            0 => ctrl_inc_by_off(&mut c, &active, off, &s_rot, &candidates),
            1 => ctrl_dec_by_off(&mut c, &active, off, &s_rot, &candidates),
            2 => {
                ctrl_inc_by_off(&mut c, &active, off, &s_rot, &candidates);
                ctrl_dec_by_off(&mut c, &active, off, &s_rot, &candidates);
            }
            _ => unreachable!(),
        }
        let mut s_ids = ids(&owned);
        if split {
            s_ids.push(counter[7].id());
        }
        let registers = vec![s_ids, ids(&counter), vec![active.id()], ids(&lenders)];
        let external = external_ids(&registers);
        Harness {
            builder: c.into_builder(),
            registers,
            external,
        }
    }

    fn build_late_inactive(roundtrip: bool) -> Harness {
        let mut c = Circuit::new();
        let owned = c.alloc_qreg_bits("q954-late.srot", 4);
        let counter = c.alloc_qreg_bits("q954-late.counter", 8);
        let parity = c.alloc_qreg("q954-late.parity");
        let gate = c.alloc_qreg("q954-late.gate");
        let x = c.alloc_qreg_bits("q954-late.x", 2);
        let y = c.alloc_qreg_bits("q954-late.y", 2);
        let a = c.alloc_qreg_bits("q954-late.a", 20);
        let b = c.alloc_qreg_bits("q954-late.b", 20);
        let q = c.alloc_qreg_bits("q954-late.q", 15);
        let lenders = c.alloc_qreg_bits("q954-late.lenders", 20);
        let candidates: Vec<&QReg> = owned
            .iter()
            .chain(counter.iter())
            .chain(std::iter::once(&parity))
            .chain(std::iter::once(&gate))
            .chain(x.iter())
            .chain(y.iter())
            .chain(a.iter())
            .chain(b.iter())
            .chain(q.iter())
            .chain(lenders.iter())
            .collect();
        let lender_refs: Vec<&QReg> = lenders.iter().collect();
        let emit = |c: &mut Circuit, inverse: bool| {
            gate_hold_counter_zero(
                c,
                &x,
                &y,
                &counter,
                &parity,
                &owned,
                &gate,
                &candidates,
                |c, active| {
                    with_arithmetic_srot_view(&owned, &counter, |s_rot| {
                        if inverse {
                            multiply_substep_windowed_inv(
                                c,
                                &a,
                                &b,
                                &q,
                                s_rot,
                                &counter[0],
                                active,
                                &lender_refs,
                                0,
                                0,
                                4,
                                4,
                            );
                        } else {
                            multiply_substep_windowed(
                                c,
                                &a,
                                &b,
                                &q,
                                s_rot,
                                &counter[0],
                                active,
                                &lender_refs,
                                0,
                                0,
                                4,
                                4,
                            );
                        }
                    });
                },
            );
        };
        emit(&mut c, false);
        if roundtrip {
            emit(&mut c, true);
        }
        let mut s_ids = ids(&owned);
        s_ids.push(counter[7].id());
        let registers = vec![
            s_ids,
            ids(&counter),
            vec![parity.id()],
            vec![gate.id()],
            ids(&x),
            ids(&y),
            ids(&a),
            ids(&b),
            ids(&q),
            ids(&lenders),
        ];
        let external = external_ids(&registers);
        Harness {
            builder: c.into_builder(),
            registers,
            external,
        }
    }

    let mut high_lane_differential_cases_checked = 0usize;
    let mut body_forward_cases_checked = 0usize;
    let mut body_reverse_cases_checked = 0usize;
    let mut body_roundtrip_cases_checked = 0usize;
    let mut max_helper_toffoli_increase = 0usize;

    for multiply in [false, true] {
        let pre = if multiply {
            [0, 1, 1 << 16, 0, 0, 1, 0]
        } else {
            [1 << 16, 1, 0, 0, 0, 1, 0]
        };
        let post = if multiply {
            [1 << 16, 1, 0, 0, 0, 1, 0]
        } else {
            [0, 1, 1 << 16, 0, 0, 1, 0]
        };
        for mode in 0..=2u8 {
            let canonical = build_body(false, multiply, mode);
            let split = build_body(true, multiply, mode);
            let input = if mode == 1 { &post } else { &pre };
            let expected = if mode == 0 { &post } else { &pre };
            let canonical_out = simulate(&canonical, input);
            let split_out = simulate(&split, input);
            assert_eq!(canonical_out.registers, expected);
            assert_eq!(split_out.registers, expected);
            assert_eq!(split_out.registers, canonical_out.registers);
            assert!(canonical_out.internal_clean && split_out.internal_clean);
            if mode == 2 {
                assert_eq!(canonical_out.phase, 0);
                assert_eq!(split_out.phase, 0);
            }
            let canonical_t = toffoli(&canonical.builder);
            let split_t = toffoli(&split.builder);
            assert_eq!(
                split.builder.peak_qubits + 1,
                canonical.builder.peak_qubits,
                "Q954 split body must remove exactly one physical shift lane"
            );
            max_helper_toffoli_increase =
                max_helper_toffoli_increase.max(split_t.saturating_sub(canonical_t));
            assert!(
                split_t <= canonical_t,
                "Q954 split body increased helper Toffoli count"
            );
            match mode {
                0 => body_forward_cases_checked += 1,
                1 => body_reverse_cases_checked += 1,
                2 => body_roundtrip_cases_checked += 1,
                _ => unreachable!(),
            }
        }
    }

    for mode in 0..=2u8 {
        let canonical = build_high_carry(false, mode);
        let split = build_high_carry(true, mode);
        let pre = [15, 1, 1, 0];
        let post_canonical = [16, 1, 1, 0];
        let post_split = [16, 129, 1, 0];
        let input_canonical = if mode == 1 { &post_canonical } else { &pre };
        let input_split = if mode == 1 { &post_split } else { &pre };
        let canonical_out = simulate(&canonical, input_canonical);
        let split_out = simulate(&split, input_split);
        let expected_canonical = if mode == 0 { &post_canonical } else { &pre };
        let expected_split = if mode == 0 { &post_split } else { &pre };
        assert_eq!(canonical_out.registers, expected_canonical);
        assert_eq!(split_out.registers, expected_split);
        assert_eq!(canonical_out.registers[0], split_out.registers[0]);
        assert_eq!(canonical_out.registers[1] & 0x7f, split_out.registers[1] & 0x7f);
        assert!(canonical_out.internal_clean && split_out.internal_clean);
        let canonical_t = toffoli(&canonical.builder);
        let split_t = toffoli(&split.builder);
        assert_eq!(
            split.builder.peak_qubits + 1,
            canonical.builder.peak_qubits,
            "Q954 split carry must remove exactly one physical shift lane"
        );
        max_helper_toffoli_increase =
            max_helper_toffoli_increase.max(split_t.saturating_sub(canonical_t));
        assert!(split_t <= canonical_t, "Q954 split carry increased Toffoli count");
        high_lane_differential_cases_checked += 1;
    }

    let late_values = [16, 128, 1, 0, 0, 1, 3, 1, 8, 0];
    let late_forward = build_late_inactive(false);
    let late_roundtrip = build_late_inactive(true);
    for (case, harness) in [&late_forward, &late_roundtrip].into_iter().enumerate() {
        let out = simulate(harness, &late_values);
        assert_eq!(out.registers, late_values);
        assert_eq!(out.registers[0], 16, "late counter[7] alias was not restored");
        assert_eq!(out.registers[1], 128, "late counter changed");
        assert!(out.internal_clean, "late inactive body retained an ancilla");
        if case == 1 {
            assert_eq!(out.phase, 0);
        }
    }

    // Canonical passenger lifetime proof: three production intervals (divide
    // forward, cancel forward, cancel reverse), each removing exactly one lane.
    let mut c = Circuit::new();
    let mut passenger = c.alloc_qreg_bits("q954-passenger", 257);
    let lower_ids = ids(&passenger[..64]);
    let all_lower_ids = ids(&passenger[..256]);
    let witness = c.alloc_qreg("q954-passenger.top-witness");
    c.cx(&passenger[256], &witness);
    let live_before = c.b.active_qubits as usize;
    for interval in 0..3 {
        release_q954_canonical_passenger_top(
            &mut c,
            &mut passenger,
            match interval {
                0 => "proof divide-forward",
                1 => "proof cancel-forward",
                _ => "proof cancel-reverse",
            },
        );
        let workspace = c.alloc_qreg_bits("q954-passenger.workspace", 7);
        assert_eq!(
            c.b.active_qubits as usize,
            live_before - 1 + workspace.len(),
            "Q954 passenger interval did not save exactly one lane"
        );
        for lane in workspace {
            c.zero_and_free(lane);
        }
        restore_q954_canonical_passenger_top(
            &mut c,
            &mut passenger,
            "q954-passenger[256]",
        );
        assert_eq!(c.b.active_qubits as usize, live_before);
    }
    let final_top_id = passenger[256].id();
    let registers = vec![lower_ids, vec![final_top_id], vec![witness.id()]];
    let mut external = all_lower_ids;
    external.push(final_top_id);
    external.push(witness.id());
    let passenger_harness = Harness {
        builder: c.into_builder(),
        registers,
        external,
    };
    let mut passenger_cases_checked = 0usize;
    for value in [0, 1, 2, 3, u64::MAX, 0x0123_4567_89ab_cdef] {
        let out = simulate(&passenger_harness, &[value, 0, 0]);
        assert_eq!(out.registers, [value, 0, 0]);
        assert!(out.internal_clean);
        passenger_cases_checked += 1;
    }

    assert_eq!(max_helper_toffoli_increase, 0);
    Q954SrotCounter7Report {
        high_lane_differential_cases_checked,
        body_forward_cases_checked,
        body_reverse_cases_checked,
        body_roundtrip_cases_checked,
        late_inactive_cases_checked: 2,
        passenger_cases_checked,
        passenger_release_intervals_checked: 3,
        srot_qubits_saved: 1,
        passenger_qubits_saved: 1,
        max_helper_toffoli_increase,
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct CanonicalLambdaLifetimeReport {
    pub cases_checked: usize,
    pub controlled_zero_cases_checked: usize,
    pub canonical_roundtrip_cases_checked: usize,
    pub lambda_lanes_before_reverse: usize,
    pub lambda_lanes_during_reverse: usize,
    pub lambda_lanes_after_reverse: usize,
    pub reverse_workspace_lanes: usize,
    pub reverse_live_qubits: usize,
    pub unreleased_reverse_live_qubits: usize,
    pub reverse_qubits_saved: usize,
    pub emitted_ops: usize,
    pub emitted_hmr: usize,
    pub emitted_resets: usize,
    pub max_internal_extra_qubits: usize,
}

/// Exercise the Q955 lambda ownership boundary with the production arithmetic.
/// The canonical product is sign-corrected, its top lane is observed and
/// released, a reverse-EEA workspace is allocated while only 256 lambda lanes
/// remain, and a fresh clean lane restores the 257-bit API. A separate
/// canonical product is then run forward and backward to prove that the
/// short-lived cleanup remains an exact matched pair.
#[doc(hidden)]
pub fn q955_off_canonical_lifetime_roundtrip_check() -> CanonicalLambdaLifetimeReport {
    use crate::circuit::{OperationType, QubitId};
    use crate::point_add::trailmix_port::arith::rfold_mbu::{
        mod_mul_canonical_mbu, mod_mul_canonical_mbu_undo,
    };
    use crate::point_add::SECP256K1_P;
    use crate::sim::Simulator;
    use ruint::aliases::U256;
    use sha3::digest::{ExtendableOutput, Update, XofReader};

    const LANES: usize = 257;
    const REVERSE_WORKSPACE_LANES: usize = 7;

    fn ids(reg: &[QReg]) -> Vec<u32> {
        reg.iter().map(QReg::id).collect()
    }

    fn load_u256<R: XofReader>(
        sim: &mut Simulator<'_, R>,
        reg: &[u32],
        value: U256,
        shot: usize,
    ) {
        for (i, &id) in reg.iter().take(256).enumerate() {
            if value.bit(i) {
                *sim.qubit_mut(QubitId(u64::from(id))) |= 1u64 << shot;
            }
        }
    }

    fn read_u256<R: XofReader>(sim: &Simulator<'_, R>, reg: &[u32], shot: usize) -> U256 {
        let mut value = U256::ZERO;
        for (i, &id) in reg.iter().take(256).enumerate() {
            if ((sim.qubit(QubitId(u64::from(id))) >> shot) & 1) != 0 {
                value.set_bit(i, true);
            }
        }
        value
    }

    assert!(lowq_q955_off_canonical_enabled());
    let mut c = Circuit::new();
    assert!(!c.b.count_only, "Q955 lifetime proof requires emitted operations");

    let a = c.alloc_qreg_bits("q955-life.a", LANES);
    let b = c.alloc_qreg_bits("q955-life.b", LANES);
    let negate = c.alloc_qreg("q955-life.negate");
    let mut lambda = c.alloc_qreg_bits("q955-life.lambda", LANES);
    mod_mul_canonical_mbu(&mut c, &lambda, &a, &b);
    controlled_field_neg_canonical(&mut c, &negate, &lambda);

    // Capture the production precondition before the lane is reset and reused.
    let top_witness = c.alloc_qreg("q955-life.top-witness");
    c.cx(&lambda[256], &top_witness);
    let lambda_lanes_before_reverse = lambda.len();
    let live_before_release = c.b.active_qubits as usize;
    release_q955_canonical_lambda_top(&mut c, &mut lambda);
    let lambda_lanes_during_reverse = lambda.len();

    let reverse_workspace =
        c.alloc_qreg_bits("q955-life.reverse-workspace", REVERSE_WORKSPACE_LANES);
    let reverse_live_qubits = c.b.active_qubits as usize;
    let unreleased_reverse_live_qubits = live_before_release + REVERSE_WORKSPACE_LANES;
    assert_eq!(
        reverse_live_qubits + 1,
        unreleased_reverse_live_qubits,
        "Q955 lambda lifetime must remove exactly one reverse-EEA lane"
    );
    for lane in &reverse_workspace {
        c.x(lane);
        c.x(lane);
    }
    for lane in reverse_workspace {
        c.zero_and_free(lane);
    }

    restore_q955_canonical_lambda_top(&mut c, &mut lambda);
    let lambda_lanes_after_reverse = lambda.len();
    assert_eq!(
        c.b.active_qubits as usize,
        live_before_release,
        "Q955 lambda lifetime must restore the pre-release live width"
    );

    // Leave this diagnostic output live so the simulator can observe exact
    // zero after the canonical forward/undo pair and the intervening sign
    // round trip used by cancellation.
    let roundtrip_temp = c.alloc_qreg_bits("q955-life.canonical-temp", LANES);
    mod_mul_canonical_mbu(&mut c, &roundtrip_temp, &lambda, &a);
    controlled_field_neg_canonical(&mut c, &negate, &roundtrip_temp);
    controlled_field_neg_canonical(&mut c, &negate, &roundtrip_temp);
    mod_mul_canonical_mbu_undo(&mut c, &roundtrip_temp, &lambda, &a);

    let a_ids = ids(&a);
    let b_ids = ids(&b);
    let lambda_ids = ids(&lambda);
    let roundtrip_temp_ids = ids(&roundtrip_temp);
    let negate_id = negate.id();
    let top_witness_id = top_witness.id();
    let external = a_ids.len()
        + b_ids.len()
        + lambda_ids.len()
        + roundtrip_temp_ids.len()
        + 2;
    assert_eq!(
        c.b.active_qubits as usize,
        external,
        "Q955 lifetime proof retained an internal quantum ancilla"
    );
    let builder = c.into_builder();

    let mut cases = Vec::with_capacity(64);
    for shot in 0..64usize {
        let low_a = U256::from((shot + 1) as u64);
        let low_b = U256::from((3 * shot + 5) as u64);
        let a_value = if shot == 1 {
            U256::ZERO
        } else if shot & 1 == 0 {
            low_a
        } else {
            SECP256K1_P - low_a
        };
        let b_value = if shot % 3 == 0 {
            SECP256K1_P - low_b
        } else {
            low_b
        };
        assert!(a_value < SECP256K1_P);
        assert!(b_value != U256::ZERO && b_value < SECP256K1_P);
        cases.push((a_value, b_value, shot & 1 != 0));
    }

    let mut seed = sha3::Shake128::default();
    seed.update(b"q955-off-canonical-lifetime-roundtrip");
    let mut xof = seed.finalize_xof();
    let mut sim = Simulator::new(
        builder.next_qubit as usize,
        builder.next_bit as usize,
        &mut xof,
    );
    sim.clear_for_shot();
    for (shot, &(a_value, b_value, negate_value)) in cases.iter().enumerate() {
        load_u256(&mut sim, &a_ids, a_value, shot);
        load_u256(&mut sim, &b_ids, b_value, shot);
        if negate_value {
            *sim.qubit_mut(QubitId(u64::from(negate_id))) |= 1u64 << shot;
        }
    }
    sim.apply_iter(builder.ops.iter());

    assert_eq!(
        sim.qubit(QubitId(u64::from(top_witness_id))),
        0,
        "canonical lambda top lane was not zero before release"
    );
    assert_eq!(
        sim.qubit(QubitId(u64::from(lambda_ids[256]))),
        0,
        "restored lambda top lane was not clean"
    );
    for (shot, &(a_value, b_value, negate_value)) in cases.iter().enumerate() {
        let product = a_value.mul_mod(b_value, SECP256K1_P);
        let expected = if negate_value && product != U256::ZERO {
            SECP256K1_P - product
        } else {
            product
        };
        assert_eq!(read_u256(&sim, &a_ids, shot), a_value, "a changed at shot {shot}");
        assert_eq!(read_u256(&sim, &b_ids, shot), b_value, "b changed at shot {shot}");
        assert_eq!(
            read_u256(&sim, &lambda_ids, shot),
            expected,
            "canonical lambda changed across its lifetime at shot {shot}"
        );
        assert_eq!(
            (sim.qubit(QubitId(u64::from(negate_id))) >> shot) & 1,
            negate_value as u64,
            "negation control changed at shot {shot}"
        );
        assert_eq!(
            read_u256(&sim, &roundtrip_temp_ids, shot),
            U256::ZERO,
            "matched canonical product/undo did not roundtrip at shot {shot}"
        );
        assert_eq!(
            (sim.qubit(QubitId(u64::from(roundtrip_temp_ids[256]))) >> shot) & 1,
            0,
            "canonical roundtrip left its top lane set at shot {shot}"
        );
    }
    assert_eq!(sim.phase, 0, "Q955 lifetime proof left phase garbage");
    let controlled_zero_cases_checked = cases
        .iter()
        .filter(|(a, b, negate)| {
            *negate && (*a).mul_mod(*b, SECP256K1_P) == U256::ZERO
        })
        .count();
    assert!(
        controlled_zero_cases_checked > 0,
        "Q955 lifetime proof must cover controlled zero"
    );

    for id in a_ids
        .iter()
        .chain(b_ids.iter())
        .chain(lambda_ids.iter())
        .chain(roundtrip_temp_ids.iter())
        .copied()
        .chain([negate_id, top_witness_id])
    {
        *sim.qubit_mut(QubitId(u64::from(id))) = 0;
    }
    for q in 0..builder.next_qubit as usize {
        assert_eq!(
            sim.qubit(QubitId(q as u64)),
            0,
            "Q955 lifetime proof left quantum ancilla q{q} dirty"
        );
    }

    CanonicalLambdaLifetimeReport {
        cases_checked: cases.len(),
        controlled_zero_cases_checked,
        canonical_roundtrip_cases_checked: cases.len(),
        lambda_lanes_before_reverse,
        lambda_lanes_during_reverse,
        lambda_lanes_after_reverse,
        reverse_workspace_lanes: REVERSE_WORKSPACE_LANES,
        reverse_live_qubits,
        unreleased_reverse_live_qubits,
        reverse_qubits_saved: unreleased_reverse_live_qubits - reverse_live_qubits,
        emitted_ops: builder.ops.len(),
        emitted_hmr: builder
            .ops
            .iter()
            .filter(|op| op.kind == OperationType::Hmr)
            .count(),
        emitted_resets: builder
            .ops
            .iter()
            .filter(|op| op.kind == OperationType::R)
            .count(),
        max_internal_extra_qubits: builder.peak_qubits as usize - external,
    }
}

#[cfg(test)]
mod tests;
