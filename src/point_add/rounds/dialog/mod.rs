//! Dialog-GCD modular inversion. This `mod.rs` holds the raw-log path (config
//! levers, per-step comparators, controlled add/sub, tobitvector / ipmul /
//! quotient / apply emitters, and the `emit_dialog_gcd_raw_pa` driver). The
//! `compressed` sidecar (round763 compressor + runway/composite scratch + the
//! `emit_dialog_gcd_compressed_sidecar_*` block-lifecycle emitters) lives in the
//! sibling module.
use super::*;

mod compressed;
mod config;
pub(crate) use compressed::*;
pub(crate) use config::*;

pub(crate) fn round84_emit_fused_square_xtail(
    b: &mut B,
    tx: &[QubitId],
    lam: &[QubitId],
    ox: &[BitId],
    p: U256,
) {
    b.set_phase("round84_fused_square_xtail_dx_sub_lam_square_lowq");
    if std::env::var("ROUND84_XTAIL_KARATSUBA").ok().as_deref() == Some("1") {
        // Squaring-aware 1-level Karatsuba square (default OFF). Overrides the
        // ROUND84_XTAIL_SCHOOLBOOK default set in configure_ecdsafail_submission_route.
        squaring_sub_from_acc_karatsuba(b, tx, lam, p);
    } else if std::env::var("ROUND84_XTAIL_WALK_SQUARE").ok().as_deref() == Some("1") {
        squaring_sub_from_acc_walk_controls_lowq(b, tx, lam, p);
    } else if std::env::var("ROUND84_XTAIL_SCHOOLBOOK").ok().as_deref() == Some("1") {
        squaring_sub_from_acc_schoolbook(b, tx, lam, p);
    } else {
        squaring_sub_from_acc_schoolbook_lowq_shift22(b, tx, lam, p);
    }
    if dialog_fuse_c_form_enabled() {
        return;
    }
    b.set_phase("round84_fused_square_xtail_add_double_ox");
    mod_add_double_qb(b, tx, ox, p);
    b.set_phase("round84_fused_square_xtail_negate_to_x3");
    mod_neg_inplace_fast(b, tx, p);
}


pub(crate) fn dialog_gcd_cmp_gt_truncated_into_width(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    flag: QubitId,
    compare_bits: usize,
) {
    assert_eq!(u.len(), v.len());
    assert!(!u.is_empty());
    let compare_bits = compare_bits.min(u.len()).max(1);
    let start = u.len() - compare_bits;
    cmp_lt_into_fast(b, &v[start..], &u[start..], flag);
}

pub(crate) fn dialog_gcd_ccx_cmp_gt_truncated_into_width(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    ctrl: QubitId,
    target: QubitId,
    compare_bits: usize,
) {
    assert_eq!(u.len(), v.len());
    assert!(!u.is_empty());
    let compare_bits = compare_bits.min(u.len()).max(1);
    let start = u.len() - compare_bits;
    ccx_cmp_lt_into_fast(b, &v[start..], &u[start..], ctrl, target);
}

pub(crate) fn dialog_gcd_branch_bits_host_comparator_enabled() -> bool {
    std::env::var("DIALOG_GCD_BRANCH_BITS_HOST_COMPARATOR")
        .ok()
        .as_deref()
        == Some("1")
}

/// Truncated controlled branch-bit comparator that hosts its borrow `c_in` +
/// `carries` transient on a borrowed clean slice (the idle future-log region)
/// when one of sufficient length is supplied, freeing the peak qubit the fresh
/// allocation would otherwise consume at the branch_bits instant. Falls back to
/// the self-allocating comparator when no slice (or a too-short one) is given, so
/// behaviour is identical to `dialog_gcd_ccx_cmp_gt_truncated_into_width` in that
/// case. Value-exact either way.
pub(crate) fn dialog_gcd_ccx_cmp_gt_truncated_into_width_hosted(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    ctrl: QubitId,
    target: QubitId,
    compare_bits: usize,
    borrowed: Option<&[QubitId]>,
) {
    assert_eq!(u.len(), v.len());
    assert!(!u.is_empty());
    let compare_bits = compare_bits.min(u.len()).max(1);
    let start = u.len() - compare_bits;
    let cmp_u = &v[start..];
    let cmp_v = &u[start..];
    let n = cmp_u.len();
    // Need c_in (1) + carries (n) = n+1 clean lanes. PARTIAL hosting: borrow the
    // future-log prefix that fits and allocate only the deficit, instead of
    // all-or-nothing (which fully self-allocs n+1 at the late GCD steps where the
    // slice runs short, pinning the branch_bits peak at 1446). The borrowed-carries
    // comparator indexes c_in and each carries[i] independently, so a gathered
    // [borrowed_prefix ++ owned] vec is value-identical; borrowed lanes are restored
    // to |0> by the measured backward inv-MAJ sweep, owned lanes are freed.
    let need = n + 1;
    let avail = borrowed.map(|s| s.len()).unwrap_or(0);
    if dialog_gcd_partial_host_comparator_enabled() && avail > 0 && avail < need {
        let slice = borrowed.expect("avail>0");
        let owned = b.alloc_qubits(need - avail);
        let mut clean: Vec<QubitId> = Vec::with_capacity(need);
        clean.extend_from_slice(slice);
        clean.extend_from_slice(&owned);
        let (c_in, carries) = clean.split_first().expect("need >= 1");
        ccx_cmp_lt_into_fast_borrowed_carries(b, cmp_u, cmp_v, ctrl, target, *c_in, &carries[..n]);
        b.free_vec(&owned);
    } else if let Some(slice) = borrowed.filter(|s| s.len() >= need) {
        let (c_in, carries) = slice.split_first().expect("slice len >= n+1 > 0");
        ccx_cmp_lt_into_fast_borrowed_carries(b, cmp_u, cmp_v, ctrl, target, *c_in, &carries[..n]);
    } else {
        ccx_cmp_lt_into_fast(b, cmp_u, cmp_v, ctrl, target);
    }
}

pub(crate) fn dialog_gcd_cmp_gt_truncated_phase_conditioned_hosted(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    ctrl: QubitId,
    phase: BitId,
    compare_bits: usize,
    borrowed: Option<&[QubitId]>,
) {
    assert_eq!(u.len(), v.len());
    assert!(!u.is_empty());
    let compare_bits = compare_bits.min(u.len()).max(1);
    let start = u.len() - compare_bits;
    let cmp_u = &v[start..];
    let cmp_v = &u[start..];
    let n = cmp_u.len();
    let need = n + 1;
    let avail = borrowed.map(|s| s.len()).unwrap_or(0);
    if dialog_gcd_partial_host_comparator_enabled() && avail > 0 && avail < need {
        let slice = borrowed.expect("avail>0");
        let owned = b.alloc_qubits(need - avail);
        let mut clean: Vec<QubitId> = Vec::with_capacity(need);
        clean.extend_from_slice(slice);
        clean.extend_from_slice(&owned);
        let (c_in, carries) = clean.split_first().expect("need >= 1");
        cmp_lt_phase_conditioned_borrowed_carries(
            b,
            cmp_u,
            cmp_v,
            *c_in,
            &carries[..n],
            ctrl,
            phase,
        );
        b.free_vec(&owned);
    } else if let Some(slice) = borrowed.filter(|s| s.len() >= need) {
        let (c_in, carries) = slice.split_first().expect("slice len >= n+1 > 0");
        cmp_lt_phase_conditioned_borrowed_carries(
            b,
            cmp_u,
            cmp_v,
            *c_in,
            &carries[..n],
            ctrl,
            phase,
        );
    } else {
        let c_in = b.alloc_qubit();
        cmp_lt_phase_conditioned_with_cin(b, cmp_u, cmp_v, c_in, ctrl, phase);
        b.free(c_in);
    }
}

fn dialog_gcd_cmp_lt_phase_conditioned_hosted(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    ctrl: QubitId,
    phase: BitId,
    borrowed: Option<&[QubitId]>,
) {
    let n = u.len();
    assert_eq!(v.len(), n);
    assert!(n > 0);
    let need = n + 1;
    let avail = borrowed.map_or(0, <[QubitId]>::len);
    if dialog_gcd_partial_host_comparator_enabled() && avail > 0 && avail < need {
        let borrowed = borrowed.expect("avail > 0");
        let owned = b.alloc_qubits(need - avail);
        let mut clean = Vec::with_capacity(need);
        clean.extend_from_slice(borrowed);
        clean.extend_from_slice(&owned);
        let (c_in, carries) = clean.split_first().expect("need >= 1");
        cmp_lt_phase_conditioned_borrowed_carries(
            b,
            u,
            v,
            *c_in,
            &carries[..n],
            ctrl,
            phase,
        );
        b.free_vec(&owned);
    } else if let Some(borrowed) = borrowed.filter(|slice| slice.len() >= need) {
        let (c_in, carries) = borrowed.split_first().expect("borrowed len >= n + 1");
        cmp_lt_phase_conditioned_borrowed_carries(
            b,
            u,
            v,
            *c_in,
            &carries[..n],
            ctrl,
            phase,
        );
    } else {
        let c_in = b.alloc_qubit();
        cmp_lt_phase_conditioned_with_cin(b, u, v, c_in, ctrl, phase);
        b.free(c_in);
    }
}

pub(crate) fn dialog_gcd_partial_host_comparator_enabled() -> bool {
    std::env::var("DIALOG_GCD_PARTIAL_HOST_COMPARATOR")
        .ok()
        .as_deref()
        != Some("0")
}


pub(crate) fn dialog_gcd_shift_right_assuming_even(b: &mut B, v: &[QubitId]) {
    assert!(!v.is_empty());
    for i in 0..v.len() - 1 {
        b.swap(v[i], v[i + 1]);
    }
}

pub(crate) fn dialog_gcd_unshift_right_assuming_even(b: &mut B, v: &[QubitId]) {
    assert!(!v.is_empty());
    for i in (0..v.len() - 1).rev() {
        b.swap(v[i], v[i + 1]);
    }
}

pub(crate) fn dialog_gcd_width_margin() -> f64 {
    // W-TRUNC safety margin added to the empirical bit-length envelope.
    // Default 37.0 reproduces pldallairedemers' baseline byte-for-byte.
    // Lowering it tightens every GCD-body width (cswap/sub/add) -> fewer
    // Toffoli, peak-neutral (early steps clamp at N). Co-tune with reroll.
    std::env::var("DIALOG_GCD_WIDTH_MARGIN")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .filter(|m| m.is_finite() && *m >= 0.0 && *m <= N as f64)
        .unwrap_or(37.0)
}

pub(crate) fn dialog_gcd_width_slope() -> f64 {
    // Per-step shrink rate of the realizable max(bitlen(u),bitlen(v)).
    // Default 0.5*1.415 = 0.7075 reproduces the baseline.
    std::env::var("DIALOG_GCD_WIDTH_SLOPE_X1000")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .filter(|s| s.is_finite() && *s > 0.0 && *s <= 4000.0)
        .map(|s| s / 1000.0)
        .unwrap_or(0.5 * 1.415)
}

pub(crate) fn dialog_gcd_tobitvector_active_width(step: usize) -> usize {
    if !dialog_gcd_raw_tobitvector_variable_width_enabled() {
        return N;
    }
    let ideal = N as f64 - (step as f64) * dialog_gcd_width_slope() + dialog_gcd_width_margin();
    let rounded = ((ideal.max(1.0) / 2.0).ceil() as usize) * 2;
    rounded
        .saturating_add(dialog_gcd_width_step_bump(step))
        .clamp(1, N)
}

fn dialog_gcd_step_map_value(env: &str, step: usize) -> usize {
    let Ok(map) = std::env::var(env) else {
        return 0;
    };
    map.split(',')
        .filter_map(|entry| {
            let (s, value) = entry.trim().split_once(':')?;
            Some((
                s.trim().parse::<usize>().ok()?,
                value.trim().parse::<usize>().ok()?,
            ))
        })
        .filter_map(|(s, value)| (s == step).then_some(value))
        .sum()
}

fn dialog_gcd_step_map_override(env: &str, step: usize) -> Option<usize> {
    let map = std::env::var(env).ok()?;
    map.split(',').rev().find_map(|entry| {
        let (raw_step, raw_value) = entry.trim().split_once(':')?;
        if raw_step.trim().parse::<usize>().ok()? != step {
            return None;
        }
        raw_value.trim().parse::<usize>().ok()
    })
}

pub(crate) fn dialog_gcd_width_step_bump(step: usize) -> usize {
    dialog_gcd_step_map_value("DIALOG_GCD_WIDTH_STEP_BUMPS", step)
}

pub(crate) fn dialog_gcd_body_step_giveback(step: usize) -> usize {
    dialog_gcd_step_map_value("DIALOG_GCD_BODY_STEP_GIVEBACKS", step)
}

pub(crate) fn dialog_gcd_fused_fold_carry_trunc_window(
    step: Option<usize>,
) -> Option<usize> {
    step.and_then(|step| {
        dialog_gcd_step_map_override("DIALOG_GCD_FOLD_CARRY_TRUNC_STEP_WINDOWS", step)
    })
    .filter(|&window| window > 0)
    .or_else(fold_only_carry_trunc_window)
    .or_else(double_carry_trunc_window)
}

/// Carry-tail truncation window for the materialized controlled sub/add BODY
/// (and its gated LOAD). Default 0 (OFF). When `w > 0`, the controlled
/// `acc -= ctrl·subtrahend` / `acc += ctrl·addend` only loads + ripples the
/// low `active_width - w` bits. The GCD work registers u/v are bounded by the
/// realizable bitlen, which sits `WIDTH_MARGIN` (=28) bits below `active_width`,
/// so the top `w <= margin` bits of both operands are 0 in the no-truncation
/// regime: the gated LOAD there is `ctrl & 0 = 0` and the body's top carries
/// are 0, so neither the load nor the carry ripple above `active_width - w`
/// affects the result. Failure mode (a step whose realizable bitlen actually
/// reaches into the truncated window) is selected away by the co-tuned reroll,
/// exactly like the global WIDTH_MARGIN — but applied to the sub/add ONLY,
/// leaving the cswap and comparator at full active_width. Returns the truncated
/// body width, clamped to >= 2.
pub(crate) fn dialog_gcd_body_carry_band_trim(step: usize) -> Option<usize> {
    let trims = std::env::var("DIALOG_GCD_BODY_CARRY_BAND_TRIMS").ok()?;
    if trims.trim().is_empty() {
        return None;
    }
    let trims: Vec<usize> = trims
        .split(',')
        .filter_map(|s| s.trim().parse::<usize>().ok())
        .collect();
    if trims.is_empty() {
        return None;
    }
    let iters = dialog_gcd_active_iterations().max(1);
    let band_size = ((iters + trims.len() - 1) / trims.len()).max(1);
    let band = (step / band_size).min(trims.len() - 1);
    Some(trims[band])
}

pub(crate) fn dialog_gcd_tobitvector_cswap_width(active_width: usize, step: usize) -> usize {
    if std::env::var("DIALOG_GCD_TOBITVECTOR_CSWAP_BODY_TRIM")
        .ok()
        .as_deref()
        == Some("1")
    {
        dialog_gcd_body_carry_trunc_width(active_width, step).min(active_width)
    } else {
        active_width
    }
}

pub(crate) fn dialog_gcd_tobitvector_shift_width(active_width: usize, step: usize) -> usize {
    if std::env::var("DIALOG_GCD_TOBITVECTOR_SHIFT_BODY_TRIM")
        .ok()
        .as_deref()
        == Some("1")
    {
        dialog_gcd_body_carry_trunc_width(active_width, step).min(active_width)
    } else {
        active_width
    }
}

pub(crate) fn dialog_gcd_body_carry_trunc_width(active_width: usize, step: usize) -> usize {
    let mut w = dialog_gcd_body_carry_band_trim(step).unwrap_or_else(|| {
        std::env::var("DIALOG_GCD_BODY_CARRY_TRUNC_W")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(0)
    });
    if dialog_gcd_trio_width_notch_enabled() && step == dialog_gcd_trio_width_notch_step() {
        w = w.saturating_add(dialog_gcd_trio_width_notch_extra());
    }
    // Multi-step binder notch (gated, default OFF). When
    // DIALOG_GCD_BINDER_NOTCH_STEPS lists `step`, trim an extra
    // DIALOG_GCD_BINDER_NOTCH_EXTRA (default 2) high bits off the materialized
    // sub/add body at THIS step too. Under the active nocin body the composite
    // scratch ask is want = 2*body_len-1, so trimming body_w by k drops the
    // owned deficit (and thus the compressed-block trio peak) by k at each
    // listed binder step. Value-exact on the reachable GCD support: at the
    // width-clamped binder steps the realizable bitlen sits WIDTH_MARGIN below
    // active_width, so the trimmed top bits of both operands are |0> (the gated
    // load there is ctrl & 0 = 0 and the carry ripple above the cut is 0).
    // Absent the env this is a no-op -> byte-identical to the accepted stream.
    if dialog_gcd_binder_notch_steps().contains(&step) {
        w = w.saturating_add(dialog_gcd_binder_notch_extra());
    }
    w = w.saturating_add(dialog_gcd_binder_notch_map_extra(step));
    w = w.saturating_sub(dialog_gcd_body_step_giveback(step));
    active_width.saturating_sub(w).max(2)
}

pub(crate) fn dialog_gcd_binder_notch_steps() -> Vec<usize> {
    std::env::var("DIALOG_GCD_BINDER_NOTCH_STEPS")
        .ok()
        .map(|s| {
            s.split(',')
                .filter_map(|t| t.trim().parse::<usize>().ok())
                .collect()
        })
        .unwrap_or_default()
}

pub(crate) fn dialog_gcd_binder_notch_extra() -> usize {
    std::env::var("DIALOG_GCD_BINDER_NOTCH_EXTRA")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(2)
}

pub(crate) fn dialog_gcd_binder_notch_map_extra(step: usize) -> usize {
    let Ok(map) = std::env::var("DIALOG_GCD_BINDER_NOTCH_MAP") else {
        return 0;
    };
    map.split(',')
        .filter_map(|entry| {
            let (s, extra) = entry.trim().split_once(':')?;
            Some((
                s.trim().parse::<usize>().ok()?,
                extra.trim().parse::<usize>().ok()?,
            ))
        })
        .filter_map(|(s, extra)| (s == step).then_some(extra))
        .sum()
}

pub(crate) fn dialog_gcd_trio_width_notch_enabled() -> bool {
    // Default-on successor from aaf9616: the current route inherited its body
    // geometry, and this one-step notch is needed to reclaim the 1306q tier.
    std::env::var("DIALOG_GCD_TRIO_WIDTH_NOTCH").ok().as_deref() != Some("0")
}

pub(crate) fn dialog_gcd_trio_width_notch_step() -> usize {
    std::env::var("DIALOG_GCD_TRIO_WIDTH_NOTCH_STEP")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(11)
}

pub(crate) fn dialog_gcd_trio_width_notch_extra() -> usize {
    std::env::var("DIALOG_GCD_TRIO_WIDTH_NOTCH_EXTRA")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(2)
}


pub(crate) fn dialog_gcd_host_gated_enabled() -> bool {
    // Port of our KAL_GZ_EARLY_RECOVER carry-pool relocation: host the
    // materialized `gated` register (width = active_width, up to 256 at peak)
    // on the provably-|0> future-log slots that already host the ripple carry,
    // instead of allocating fresh ancilla. The borrowed slice (when long enough
    // for carry + gated = 2n-1) is split: [..n-1] = carry, [n-1..2n-1] = gated.
    // Both are restored to |0> (carry by the adder, gated by measurement-clear),
    // so the future-log slots are clean for the future blocks that own them.
    // Peak-neutral->down: removes the +256 fresh ancilla at the GCD-body peak.
    // Default off = byte-identical baseline.
    std::env::var("DIALOG_GCD_HOST_GATED").ok().as_deref() == Some("1")
}

pub(crate) fn dialog_gcd_body_host_cin_enabled() -> bool {
    // When the odd-u low-bit fastpath is active (body_start>=1), the low gated
    // slot gated[0] is never loaded or cleared, so it stays |0> across the body
    // and is distinct from the operands and the borrowed carry lane. Hosting the
    // Cuccaro carry-in there instead of a fresh alloc removes the single qubit
    // that pinned the materialized add/sub BODY one slot above the marker tier.
    // Value-exact (c_in=0 is the carry-in either way; returned to |0>).
    std::env::var("DIALOG_GCD_BODY_HOST_CIN").ok().as_deref() == Some("1")
}

pub(crate) fn dialog_gcd_selected_body_nocin_enabled() -> bool {
    // Successor to BODY_HOST_CIN for the odd-lowbit fastpath (body_start>=1):
    // the materialized selected add/sub body consumes NO physical incoming-carry
    // lane at all. The carry/borrow into body_start=1 is semantically zero on the
    // reachable GCD support (subtrahend[0]=1, acc[0]=ctrl), so the Cuccaro chain
    // is seeded from the known-zero with the c_in register folded out entirely
    // (see cuccaro_{add,sub}_fast_borrowed_carries_no_cin). This drops the
    // selected-body host demand from 2*body_w-1 to 2*body_w-3 (one structural gap
    // lane + the former c_in lane both vanish), moving the three GCD tobitvector
    // siblings off the 1320 tier without reusing the wrapper-unsafe gap-as-c_in
    // slice that the COMPACT probe (closed) tried. Default off until traced.
    matches!(
        std::env::var("DIALOG_GCD_SELECTED_BODY_NOCIN").ok().as_deref(),
        Some("1") | Some("2")
    )
}

/// Diagnostic mode 2: use the no-c_in BODY but keep the legacy `2n-1` composite
/// pool and BODY_HOST_CIN slice offsets (gated = c[n-1..2n-1], its [0] left
/// unused/clean). This isolates the body arithmetic from the host repack — it
/// yields no peak win (pool unchanged) but, if eval is 0/0/0, proves the no-c_in
/// body is route-correct and any failure under mode 1 is in the host compaction.
pub(crate) fn dialog_gcd_selected_body_nocin_keep_pool() -> bool {
    std::env::var("DIALOG_GCD_SELECTED_BODY_NOCIN").ok().as_deref() == Some("2")
}

/// Per-step count of high source bits streamed through the controlled low-q
/// suffix instead of being materialized. A value of one uses the cheaper
/// top-bit specialization. Default off when the map is absent.
pub(crate) fn dialog_gcd_selected_body_stream_suffix_bits(step: usize, body_len: usize) -> usize {
    let Ok(map) = std::env::var("DIALOG_GCD_SELECTED_BODY_STREAM_SUFFIX_MAP") else {
        return 0;
    };
    map.split(',')
        .find_map(|entry| {
            let (entry_step, entry_bits) = entry.trim().split_once(':')?;
            let entry_step = entry_step.parse::<usize>().ok()?;
            let entry_bits = entry_bits.parse::<usize>().ok()?;
            (entry_step == step).then_some(entry_bits)
        })
        .unwrap_or(0)
        .min(body_len.saturating_sub(1))
}

pub(crate) fn dialog_gcd_selected_body_stream_top_enabled(step: usize, body_len: usize) -> bool {
    dialog_gcd_selected_body_stream_suffix_bits(step, body_len) == 1
}

pub(crate) fn dialog_gcd_selected_body_stream_topclean_bits(
    step: usize,
    prefix_len: usize,
) -> usize {
    if prefix_len <= 1 {
        return 0;
    }
    let global = std::env::var("DIALOG_GCD_SELECTED_BODY_STREAM_TOPCLEAN")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0);
    let mapped = dialog_gcd_step_map_value("DIALOG_GCD_SELECTED_BODY_STREAM_TOPCLEAN_MAP", step);
    global.max(mapped).min(prefix_len - 1)
}

pub(crate) fn dialog_gcd_late_borrow_uv_high_enabled() -> bool {
    std::env::var("DIALOG_GCD_LATE_BORROW_UV_HIGH")
        .ok()
        .as_deref()
        == Some("1")
}

/// Pick the carry/gated borrow slice for a GCD step. Prefer the compressed
/// future-log; when it is too short to host the full gated(n)+carry(n-1) lane
/// (late steps, where the compressed future region has shrunk), fall back to the
/// high zero bits of `u`. By the same premise the width truncation relies on,
/// `u < 2^active_width` here so `u[active_width..]` is |0>; it is already
/// allocated, so borrowing it as scratch is peak-neutral and adds no failure
/// modes (any input with nonzero u-high already fails the truncation). The
/// returned slice is disjoint from `u[..active_width]` and the `v` accumulator.
pub(crate) fn dialog_gcd_pick_borrow_slice<'a>(
    future: Option<&'a [QubitId]>,
    u: &'a [QubitId],
    active_width: usize,
) -> Option<&'a [QubitId]> {
    if dialog_gcd_late_borrow_uv_high_enabled() && active_width >= 1 {
        let want = 2 * active_width - 1;
        let short = future.map_or(true, |s| s.len() < want);
        if short && u.len() >= active_width + want {
            return Some(&u[active_width..active_width + want]);
        }
    }
    future
}

pub(crate) fn dialog_gcd_controlled_sub_selected(
    b: &mut B,
    subtrahend: &[QubitId],
    acc: &[QubitId],
    ctrl: QubitId,
    borrowed_carries: Option<&[QubitId]>,
    step: usize,
) {
    assert_eq!(subtrahend.len(), acc.len());
    assert!(!subtrahend.is_empty());
    if dialog_gcd_raw_tobitvector_materialized_sub_enabled() {
        let n = subtrahend.len();
        let body_w = dialog_gcd_body_carry_trunc_width(n, step);
        let odd_lowbit_fast = dialog_gcd_odd_u_lowbit_fastpath_enabled();
        let body_start = if odd_lowbit_fast { 1 } else { 0 };
        let body_len = body_w.saturating_sub(body_start);
        let stream_suffix = dialog_gcd_selected_body_stream_suffix_bits(step, body_len);
        let nocin_need = if stream_suffix >= 2
            && !dialog_gcd_selected_body_nocin_keep_pool()
        {
            2 * (body_len - stream_suffix) + 1
        } else if dialog_gcd_selected_body_stream_top_enabled(step, body_len)
            && !dialog_gcd_selected_body_nocin_keep_pool()
            && body_len >= 2
        {
            2 * (body_len - 1)
        } else if dialog_gcd_selected_body_nocin_keep_pool() {
            // Legacy gated offset c[n..n+body_len] needs the full 2n-1 pool.
            (n + body_len).max(2 * body_len - 1)
        } else {
            2 * body_len - 1
        };
        let nocin = dialog_gcd_selected_body_nocin_enabled()
            && body_start >= 1
            && body_len >= 1
            && borrowed_carries.map_or(false, |c| c.len() >= nocin_need);
        if nocin {
            if stream_suffix >= 2 {
                let prefix_len = body_len - stream_suffix;
                let c = borrowed_carries.expect("nocin requires borrowed carries");
                let (carries, rest) = c.split_at(prefix_len);
                let (gated, rest) = rest.split_at(prefix_len);
                let scratch = rest[0];
                b.set_phase("dialog_gcd_raw_tobitvector_materialized_sub_load");
                for j in 0..prefix_len {
                    b.ccx(ctrl, subtrahend[body_start + j], gated[j]);
                }
                b.cx(ctrl, acc[0]);
                b.set_phase("dialog_gcd_raw_tobitvector_materialized_sub_body");
                cuccaro_sub_fast_prefix_ctrl_suffix_no_cin(
                    b,
                    gated,
                    &subtrahend[body_start + prefix_len..body_w],
                    &acc[body_start..body_w],
                    ctrl,
                    carries,
                    scratch,
                );
                b.set_phase("dialog_gcd_raw_tobitvector_materialized_sub_clear");
                for j in 0..prefix_len {
                    let m = b.alloc_bit();
                    b.hmr(gated[j], m);
                    b.cz_if(ctrl, subtrahend[body_start + j], m);
                }
                return;
            }
            if dialog_gcd_selected_body_stream_top_enabled(step, body_len) && body_len >= 2 {
                let lower_len = body_len - 1;
                let c = borrowed_carries.expect("nocin requires borrowed carries");
                let (carries, gated) = c.split_at(lower_len);
                let gated = &gated[..lower_len];
                b.set_phase("dialog_gcd_raw_tobitvector_materialized_sub_load");
                for j in 0..lower_len {
                    b.ccx(ctrl, subtrahend[body_start + j], gated[j]);
                }
                b.cx(ctrl, acc[0]);
                b.set_phase("dialog_gcd_raw_tobitvector_materialized_sub_body");
                b.ccx(ctrl, subtrahend[body_w - 1], acc[body_w - 1]);
                cuccaro_sub_fast_low_to_ext_borrowed_carries_no_cin(
                    b,
                    gated,
                    &acc[body_start..body_w],
                    carries,
                );
                b.set_phase("dialog_gcd_raw_tobitvector_materialized_sub_clear");
                for j in 0..lower_len {
                    let m = b.alloc_bit();
                    b.hmr(gated[j], m);
                    b.cz_if(ctrl, subtrahend[body_start + j], m);
                }
                return;
            }
            // No-physical-c_in body: host demand 2*body_len-1 (== 2*body_w-3).
            // carries = borrowed[..body_len-1], gated = borrowed[body_len-1..2*body_len-1].
            // Diagnostic keep-pool (mode 2) instead uses the BODY_HOST_CIN offsets
            // (carries low, gated = c[n-1+1..] on the legacy 2n-1 pool) to isolate
            // the body arithmetic from the host repack.
            let c = borrowed_carries.expect("nocin requires borrowed carries");
            let (carries, gated): (&[QubitId], &[QubitId]) =
                if dialog_gcd_selected_body_nocin_keep_pool() {
                    // Legacy gated = c[n-1..2n-1]; gated[0]=c[n-1] is the unused
                    // (clean) former c_in slot, so operand lands on c[n..2n-1].
                    let carry_need = body_len - 1;
                    (&c[..carry_need], &c[n..n + body_len])
                } else {
                    let carry_need = body_len - 1;
                    (&c[..carry_need], &c[carry_need..carry_need + body_len])
                };
            b.set_phase("dialog_gcd_raw_tobitvector_materialized_sub_load");
            for j in 0..body_len {
                b.ccx(ctrl, subtrahend[body_start + j], gated[j]);
            }
            // Reachable GCD states have subtrahend[0]=1 and acc[0]=ctrl here:
            // ctrl - ctrl has result bit 0 and no borrow into bit 1 (the omitted
            // c_in). This is exactly the premise the no-c_in body relies on.
            b.cx(ctrl, acc[0]);
            b.set_phase("dialog_gcd_raw_tobitvector_materialized_sub_body");
            cuccaro_sub_fast_borrowed_carries_no_cin(
                b,
                gated,
                &acc[body_start..body_w],
                carries,
            );
            b.set_phase("dialog_gcd_raw_tobitvector_materialized_sub_clear");
            for j in 0..body_len {
                let m = b.alloc_bit();
                b.hmr(gated[j], m);
                b.cz_if(ctrl, subtrahend[body_start + j], m);
            }
            return;
        }
        // Host the gated register on the tail of the borrowed clean slice when
        // it is long enough for both carry (n-1) and gated (n).
        let gated_host: Option<&[QubitId]> = if dialog_gcd_host_gated_enabled() {
            borrowed_carries.and_then(|c| {
                if c.len() >= 2 * n - 1 {
                    Some(&c[n - 1..2 * n - 1])
                } else {
                    None
                }
            })
        } else {
            None
        };
        let mut gated_owned: Vec<QubitId> = Vec::new();
        let gated: &[QubitId] = match gated_host {
            Some(h) => h,
            None => {
                gated_owned = b.alloc_qubits(n);
                gated_owned.as_slice()
            }
        };
        b.set_phase("dialog_gcd_raw_tobitvector_materialized_sub_load");
        for i in body_start..body_w {
            b.ccx(ctrl, subtrahend[i], gated[i]);
        }
        if odd_lowbit_fast {
            // Reachable GCD states have subtrahend[0]=1 and acc[0]=ctrl here:
            // ctrl - ctrl has result bit 0 and no borrow into bit 1.
            b.cx(ctrl, acc[0]);
        }
        b.set_phase("dialog_gcd_raw_tobitvector_materialized_sub_body");
        if body_start < body_w {
            if let Some(carries) =
                borrowed_carries.filter(|carries| carries.len() >= body_len.saturating_sub(1))
            {
                if dialog_gcd_body_host_cin_enabled() && body_start >= 1 {
                    // gated[0] is unused (load/clear start at body_start) and |0>:
                    // use it as the Cuccaro carry-in, dropping the fresh c_in alloc.
                    cuccaro_sub_fast_borrowed_carries(
                        b,
                        &gated[body_start..body_w],
                        &acc[body_start..body_w],
                        gated[0],
                        &carries[..body_len.saturating_sub(1)],
                    );
                } else {
                    sub_nbit_qq_fast_borrowed_carries(
                        b,
                        &gated[body_start..body_w],
                        &acc[body_start..body_w],
                        &carries[..body_len.saturating_sub(1)],
                    );
                }
            } else {
                sub_nbit_qq_fast(b, &gated[body_start..body_w], &acc[body_start..body_w]);
            }
        }
        b.set_phase("dialog_gcd_raw_tobitvector_materialized_sub_clear");
        for i in body_start..body_w {
            let m = b.alloc_bit();
            b.hmr(gated[i], m);
            b.cz_if(ctrl, subtrahend[i], m);
        }
        if gated_host.is_none() {
            b.free_vec(&gated_owned);
        }
    } else {
        let n = subtrahend.len();
        if dialog_gcd_ctrl_body_vented_enabled() {
            if let Some(vents) =
                borrowed_carries.filter(|c| n >= 2 && c.len() >= n - 1)
            {
                cuccaro_sub_ctrl_vented(b, subtrahend, acc, ctrl, &vents[..n - 1]);
                return;
            }
        }
        cucc_sub_ctrl_lowq(b, subtrahend, acc, ctrl);
    }
}

pub(crate) fn dialog_gcd_controlled_add_selected(
    b: &mut B,
    addend: &[QubitId],
    acc: &[QubitId],
    ctrl: QubitId,
    borrowed_carries: Option<&[QubitId]>,
    step: usize,
) {
    assert_eq!(addend.len(), acc.len());
    assert!(!addend.is_empty());
    if dialog_gcd_raw_tobitvector_materialized_sub_enabled() {
        let n = addend.len();
        let body_w = dialog_gcd_body_carry_trunc_width(n, step);
        let odd_lowbit_fast = dialog_gcd_odd_u_lowbit_fastpath_enabled();
        let body_start = if odd_lowbit_fast { 1 } else { 0 };
        let body_len = body_w.saturating_sub(body_start);
        let stream_suffix = dialog_gcd_selected_body_stream_suffix_bits(step, body_len);
        let nocin_need = if stream_suffix >= 2
            && !dialog_gcd_selected_body_nocin_keep_pool()
        {
            2 * (body_len - stream_suffix) + 1
        } else if dialog_gcd_selected_body_stream_top_enabled(step, body_len)
            && !dialog_gcd_selected_body_nocin_keep_pool()
            && body_len >= 2
        {
            2 * (body_len - 1)
        } else if dialog_gcd_selected_body_nocin_keep_pool() {
            // Legacy gated offset c[n..n+body_len] needs the full 2n-1 pool.
            (n + body_len).max(2 * body_len - 1)
        } else {
            2 * body_len - 1
        };
        let nocin = dialog_gcd_selected_body_nocin_enabled()
            && body_start >= 1
            && body_len >= 1
            && borrowed_carries.map_or(false, |c| c.len() >= nocin_need);
        if nocin {
            if stream_suffix >= 2 {
                let prefix_len = body_len - stream_suffix;
                let c = borrowed_carries.expect("nocin requires borrowed carries");
                let (carries, rest) = c.split_at(prefix_len);
                let (gated, rest) = rest.split_at(prefix_len);
                let scratch = rest[0];
                b.set_phase("dialog_gcd_raw_tobitvector_materialized_add_load");
                for j in 0..prefix_len {
                    b.ccx(ctrl, addend[body_start + j], gated[j]);
                }
                b.cx(ctrl, acc[0]);
                b.set_phase("dialog_gcd_raw_tobitvector_materialized_add_body");
                cuccaro_add_fast_prefix_ctrl_suffix_no_cin(
                    b,
                    gated,
                    &addend[body_start + prefix_len..body_w],
                    &acc[body_start..body_w],
                    ctrl,
                    carries,
                    scratch,
                );
                b.set_phase("dialog_gcd_raw_tobitvector_materialized_add_clear");
                for j in 0..prefix_len {
                    let m = b.alloc_bit();
                    b.hmr(gated[j], m);
                    b.cz_if(ctrl, addend[body_start + j], m);
                }
                return;
            }
            if dialog_gcd_selected_body_stream_top_enabled(step, body_len) && body_len >= 2 {
                let lower_len = body_len - 1;
                let c = borrowed_carries.expect("nocin requires borrowed carries");
                let (carries, gated) = c.split_at(lower_len);
                let gated = &gated[..lower_len];
                b.set_phase("dialog_gcd_raw_tobitvector_materialized_add_load");
                for j in 0..lower_len {
                    b.ccx(ctrl, addend[body_start + j], gated[j]);
                }
                b.cx(ctrl, acc[0]);
                b.set_phase("dialog_gcd_raw_tobitvector_materialized_add_body");
                cuccaro_add_fast_low_to_ext_borrowed_carries_no_cin(
                    b,
                    gated,
                    &acc[body_start..body_w],
                    carries,
                );
                b.ccx(ctrl, addend[body_w - 1], acc[body_w - 1]);
                b.set_phase("dialog_gcd_raw_tobitvector_materialized_add_clear");
                for j in 0..lower_len {
                    let m = b.alloc_bit();
                    b.hmr(gated[j], m);
                    b.cz_if(ctrl, addend[body_start + j], m);
                }
                return;
            }
            // No-physical-c_in inverse body: host demand 2*body_len-1 (==2*body_w-3).
            let c = borrowed_carries.expect("nocin requires borrowed carries");
            let (carries, gated): (&[QubitId], &[QubitId]) =
                if dialog_gcd_selected_body_nocin_keep_pool() {
                    let carry_need = body_len - 1;
                    (&c[..carry_need], &c[n..n + body_len])
                } else {
                    let carry_need = body_len - 1;
                    (&c[..carry_need], &c[carry_need..carry_need + body_len])
                };
            b.set_phase("dialog_gcd_raw_tobitvector_materialized_add_load");
            for j in 0..body_len {
                b.ccx(ctrl, addend[body_start + j], gated[j]);
            }
            // In reverse, acc[0] is zero after unshift and addend[0]=1: adding
            // ctrl sets the low result bit with no carry into bit 1 (omitted c_in).
            b.cx(ctrl, acc[0]);
            b.set_phase("dialog_gcd_raw_tobitvector_materialized_add_body");
            cuccaro_add_fast_borrowed_carries_no_cin(
                b,
                gated,
                &acc[body_start..body_w],
                carries,
            );
            b.set_phase("dialog_gcd_raw_tobitvector_materialized_add_clear");
            for j in 0..body_len {
                let m = b.alloc_bit();
                b.hmr(gated[j], m);
                b.cz_if(ctrl, addend[body_start + j], m);
            }
            return;
        }
        let gated_host: Option<&[QubitId]> = if dialog_gcd_host_gated_enabled() {
            borrowed_carries.and_then(|c| {
                if c.len() >= 2 * n - 1 {
                    Some(&c[n - 1..2 * n - 1])
                } else {
                    None
                }
            })
        } else {
            None
        };
        let mut gated_owned: Vec<QubitId> = Vec::new();
        let gated: &[QubitId] = match gated_host {
            Some(h) => h,
            None => {
                gated_owned = b.alloc_qubits(n);
                gated_owned.as_slice()
            }
        };
        b.set_phase("dialog_gcd_raw_tobitvector_materialized_add_load");
        for i in body_start..body_w {
            b.ccx(ctrl, addend[i], gated[i]);
        }
        if odd_lowbit_fast {
            // In reverse, acc[0] is zero after unshift and addend[0]=1:
            // adding ctrl sets the low result bit with no carry into bit 1.
            b.cx(ctrl, acc[0]);
        }
        b.set_phase("dialog_gcd_raw_tobitvector_materialized_add_body");
        if body_start < body_w {
            if let Some(carries) =
                borrowed_carries.filter(|carries| carries.len() >= body_len.saturating_sub(1))
            {
                if dialog_gcd_body_host_cin_enabled() && body_start >= 1 {
                    // gated[0] is unused (load/clear start at body_start) and |0>:
                    // use it as the Cuccaro carry-in, dropping the fresh c_in alloc.
                    cuccaro_add_fast_borrowed_carries(
                        b,
                        &gated[body_start..body_w],
                        &acc[body_start..body_w],
                        gated[0],
                        &carries[..body_len.saturating_sub(1)],
                    );
                } else {
                    add_nbit_qq_fast_borrowed_carries(
                        b,
                        &gated[body_start..body_w],
                        &acc[body_start..body_w],
                        &carries[..body_len.saturating_sub(1)],
                    );
                }
            } else {
                add_nbit_qq_fast(b, &gated[body_start..body_w], &acc[body_start..body_w]);
            }
        }
        b.set_phase("dialog_gcd_raw_tobitvector_materialized_add_clear");
        for i in body_start..body_w {
            let m = b.alloc_bit();
            b.hmr(gated[i], m);
            b.cz_if(ctrl, addend[i], m);
        }
        if gated_host.is_none() {
            b.free_vec(&gated_owned);
        }
    } else {
        let n = addend.len();
        if dialog_gcd_ctrl_body_vented_enabled() {
            if let Some(vents) =
                borrowed_carries.filter(|c| n >= 2 && c.len() >= n - 1)
            {
                cuccaro_add_ctrl_vented(b, addend, acc, ctrl, &vents[..n - 1]);
                return;
            }
        }
        cucc_add_ctrl_lowq(b, addend, acc, ctrl);
    }
}

pub(crate) fn dialog_gcd_future_log_carry_slice(
    dialog_log: &[QubitId],
    step: usize,
    active_width: usize,
) -> Option<&[QubitId]> {
    if !dialog_gcd_raw_tobitvector_borrow_future_log_carries_enabled() {
        return None;
    }
    let carry_need = active_width.saturating_sub(1);
    let want = if dialog_gcd_host_gated_enabled() {
        2 * active_width - 1
    } else {
        carry_need
    };
    let start = 2 * (step + 1);
    dialog_log
        .get(start..)
        .filter(|future| future.len() >= carry_need)
        .map(|future| &future[..future.len().min(want)])
}

pub(crate) fn emit_dialog_gcd_raw_tobitvector_steps(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    dialog_log: &[QubitId],
) {
    assert_eq!(u.len(), N);
    assert_eq!(v.len(), N);
    assert!(dialog_log.len() >= 2 * dialog_gcd_active_iterations());

    for step in 0..dialog_gcd_active_iterations() {
        let b0 = dialog_log[2 * step];
        let b0_and_b1 = dialog_log[2 * step + 1];
        let cmp = b.alloc_qubit();
        let active_width = dialog_gcd_tobitvector_active_width(step);
        let u_active = &u[..active_width];
        let v_active = &v[..active_width];
        let compare_bits = dialog_gcd_compare_bits_for_step(step, active_width);

        b.set_phase("dialog_gcd_raw_tobitvector_branch_bits");
        b.cx(v[0], b0);
        if dialog_gcd_fused_branch_bits_enabled() {
            dialog_gcd_ccx_cmp_gt_truncated_into_width(
                b,
                u_active,
                v_active,
                b0,
                b0_and_b1,
                compare_bits,
            );
        } else {
            dialog_gcd_cmp_gt_truncated_into_width(b, u_active, v_active, cmp, compare_bits);
            b.ccx(b0, cmp, b0_and_b1);
            dialog_gcd_cmp_gt_truncated_into_width(b, u_active, v_active, cmp, compare_bits);
        }
        b.free(cmp);

        b.set_phase("dialog_gcd_raw_tobitvector_cswap");
        for (i, (&ui, &vi)) in u_active.iter().zip(v_active.iter()).enumerate() {
            if i == 0 && dialog_gcd_odd_u_lowbit_fastpath_enabled() {
                continue;
            }
            cswap(b, b0_and_b1, ui, vi);
        }

        b.set_phase("dialog_gcd_raw_tobitvector_subtract");
        let borrowed_carries = dialog_gcd_future_log_carry_slice(dialog_log, step, active_width);
        dialog_gcd_controlled_sub_selected(b, u_active, v_active, b0, borrowed_carries, step);

        b.set_phase("dialog_gcd_raw_tobitvector_shift");
        dialog_gcd_shift_right_assuming_even(b, v_active);
    }
}

pub(crate) fn emit_dialog_gcd_raw_tobitvector_steps_reverse(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    dialog_log: &[QubitId],
) {
    assert_eq!(u.len(), N);
    assert_eq!(v.len(), N);
    assert!(dialog_log.len() >= 2 * dialog_gcd_active_iterations());

    for step in (0..dialog_gcd_active_iterations()).rev() {
        let b0 = dialog_log[2 * step];
        let b0_and_b1 = dialog_log[2 * step + 1];
        let cmp = b.alloc_qubit();
        let active_width = dialog_gcd_tobitvector_active_width(step);
        let u_active = &u[..active_width];
        let v_active = &v[..active_width];
        let compare_bits = dialog_gcd_compare_bits_for_step(step, active_width);

        b.set_phase("dialog_gcd_raw_tobitvector_reverse_unshift");
        dialog_gcd_unshift_right_assuming_even(b, v_active);

        b.set_phase("dialog_gcd_raw_tobitvector_reverse_add");
        let borrowed_carries = dialog_gcd_future_log_carry_slice(dialog_log, step, active_width);
        dialog_gcd_controlled_add_selected(b, u_active, v_active, b0, borrowed_carries, step);

        b.set_phase("dialog_gcd_raw_tobitvector_reverse_cswap");
        for (i, (&ui, &vi)) in u_active.iter().zip(v_active.iter()).enumerate() {
            if i == 0 && dialog_gcd_odd_u_lowbit_fastpath_enabled() {
                continue;
            }
            cswap(b, b0_and_b1, ui, vi);
        }

        b.set_phase("dialog_gcd_raw_tobitvector_reverse_branch_bits");
        if dialog_gcd_fused_branch_bits_enabled() {
            dialog_gcd_ccx_cmp_gt_truncated_into_width(
                b,
                u_active,
                v_active,
                b0,
                b0_and_b1,
                compare_bits,
            );
        } else {
            dialog_gcd_cmp_gt_truncated_into_width(b, u_active, v_active, cmp, compare_bits);
            b.ccx(b0, cmp, b0_and_b1);
            dialog_gcd_cmp_gt_truncated_into_width(b, u_active, v_active, cmp, compare_bits);
        }
        b.free(cmp);
        b.cx(v[0], b0);
    }
}


pub(crate) fn dialog_gcd_cmod_add_pseudomersenne_lowq(
    b: &mut B,
    acc: &[QubitId],
    a: &[QubitId],
    ctrl: QubitId,
    p: U256,
) {
    assert_eq!(acc.len(), N);
    assert_eq!(a.len(), N);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1u64));

    let (acc_ext, acc_ovf) = ext_reg(b, acc);
    let a_ovf = b.alloc_qubit();
    let mut a_ext = a.to_vec();
    a_ext.push(a_ovf);
    let c_in = b.alloc_qubit();
    let scratch = b.alloc_qubit();

    b.set_phase("dialog_gcd_direct_special_cadd_raw_sum");
    cuccaro_add_ctrl_lowq(b, &a_ext, &acc_ext, ctrl, c_in, scratch);
    b.free(scratch);
    b.free(c_in);
    b.free(a_ovf);

    // If the controlled 256-bit add overflowed, subtract p by adding
    // c = 2^256 - p to the low word.  The low slice is the explicit
    // approximation knob: carry beyond this window is treated as a rare
    // arithmetic failure branch, not as phase dirt.
    b.set_phase("dialog_gcd_direct_special_overflow_fold");
    cadd_nbit_const_fast(b, &acc[..DIALOG_GCD_SPECIAL_ADD_LSBS], c, acc_ovf);

    // For successful branches this is the exact overflow cleanup identity:
    // after subtracting p, the final low word is smaller than the addend iff
    // the overflow branch happened.  The omitted no-overflow sum>=p case is
    // the approximation budgeted by the caller.
    b.set_phase("dialog_gcd_direct_special_overflow_clean");
    cmp_lt_into(b, acc, a, acc_ovf);
    unext_reg(b, acc_ovf);
}

pub(crate) fn dialog_gcd_cmod_add_materialized_pseudomersenne(
    b: &mut B,
    acc: &[QubitId],
    a: &[QubitId],
    ctrl: QubitId,
    p: U256,
) {
    dialog_gcd_cmod_add_materialized_pseudomersenne_at_step(b, acc, a, ctrl, p, None);
}

pub(crate) fn dialog_gcd_cmod_add_materialized_pseudomersenne_at_step(
    b: &mut B,
    acc: &[QubitId],
    a: &[QubitId],
    ctrl: QubitId,
    p: U256,
    step: Option<usize>,
) {
    dialog_gcd_cmod_add_materialized_pseudomersenne_with_clean_scratch_at_step(
        b,
        acc,
        a,
        ctrl,
        p,
        &[],
        step,
    );
}

pub(crate) fn dialog_gcd_cmod_add_materialized_pseudomersenne_with_clean_scratch(
    b: &mut B,
    acc: &[QubitId],
    a: &[QubitId],
    ctrl: QubitId,
    p: U256,
    clean_scratch: &[QubitId],
) {
    dialog_gcd_cmod_add_materialized_pseudomersenne_with_clean_scratch_at_step(
        b,
        acc,
        a,
        ctrl,
        p,
        clean_scratch,
        None,
    );
}

pub(crate) fn dialog_gcd_cmod_add_materialized_pseudomersenne_with_clean_scratch_at_step(
    b: &mut B,
    acc: &[QubitId],
    a: &[QubitId],
    ctrl: QubitId,
    p: U256,
    clean_scratch: &[QubitId],
    step: Option<usize>,
) {
    assert_eq!(acc.len(), N);
    assert_eq!(a.len(), N);
    if let Some(blocks) = dialog_gcd_apply_chunked_f_blocks()
        .filter(|_| dialog_gcd_raw_apply_truncated_clean_enabled())
    {
        dialog_gcd_cmod_add_materialized_pseudomersenne_chunked(
            b,
            acc,
            a,
            ctrl,
            p,
            blocks,
            clean_scratch,
            step,
        );
        return;
    }
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1u64));

    let f = b.alloc_qubits(N);
    b.set_phase("dialog_gcd_materialized_special_load_addend");
    for i in 0..N {
        b.ccx(ctrl, a[i], f[i]);
    }

    let (acc_ext, acc_ovf) = ext_reg(b, acc);
    let c_in = b.alloc_qubit();

    b.set_phase("dialog_gcd_materialized_special_raw_sum");
    if let Some(w) = dialog_gcd_apply_window_blocks() {
        cuccaro_add_fast_windowed_low_to_ext(b, &f, &acc_ext, c_in, w);
    } else {
        let f_ovf = b.alloc_qubit();
        let mut f_ext = f.clone();
        f_ext.push(f_ovf);
        cuccaro_add_fast(b, &f_ext, &acc_ext, c_in);
        b.free(f_ovf);
    }
    b.free(c_in);

    b.set_phase("dialog_gcd_materialized_special_overflow_fold");
    if let Some(w) = dialog_gcd_special_fold_carry_trunc_window(step) {
        cadd_nbit_const_direct_trunc_fast(b, &acc[..DIALOG_GCD_SPECIAL_ADD_LSBS], c, acc_ovf, w);
    } else {
        cadd_nbit_const_fast(b, &acc[..DIALOG_GCD_SPECIAL_ADD_LSBS], c, acc_ovf);
    }

    b.set_phase("dialog_gcd_materialized_special_overflow_clean");
    if dialog_gcd_raw_apply_truncated_clean_enabled() {
        let compare_start = N - dialog_gcd_special_overflow_clean_compare_bits(step);
        cmp_lt_into_fast(b, &acc[compare_start..], &f[compare_start..], acc_ovf);
    } else {
        cmp_lt_into(b, acc, &f, acc_ovf);
    }
    unext_reg(b, acc_ovf);

    b.set_phase("dialog_gcd_materialized_special_clear_addend");
    for i in 0..N {
        let m = b.alloc_bit();
        b.hmr(f[i], m);
        b.cz_if(ctrl, a[i], m);
    }
    b.free_vec(&f);
}

pub(crate) fn dialog_gcd_measured_apply_sub_enabled() -> bool {
    std::env::var("DIALOG_GCD_MEASURED_APPLY_SUB")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_apply_window_blocks() -> Option<usize> {
    std::env::var("DIALOG_GCD_APPLY_WINDOW_BLOCKS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|&w| w >= 2)
}

fn dialog_gcd_clean_truncated_underflow_with_borrowed(
    b: &mut B,
    acc: &[QubitId],
    a: &[QubitId],
    ctrl: QubitId,
    acc_ovf: QubitId,
    step: Option<usize>,
    borrowed: Option<&[QubitId]>,
) {
    let compare_start = N - dialog_gcd_special_underflow_clean_compare_bits(step);
    for &q in &a[compare_start..] {
        b.x(q);
    }
    if dialog_gcd_special_clean_conditional_replay_enabled() {
        let phase = b.alloc_bit();
        b.hmr(acc_ovf, phase);
        b.z_if(ctrl, phase);
        dialog_gcd_cmp_lt_phase_conditioned_hosted(
            b,
            &acc[compare_start..],
            &a[compare_start..],
            ctrl,
            phase,
            borrowed,
        );
    } else {
        b.cx(ctrl, acc_ovf);
        ccx_cmp_lt_into_fast(b, &acc[compare_start..], &a[compare_start..], ctrl, acc_ovf);
    }
    for &q in &a[compare_start..] {
        b.x(q);
    }
}

pub(crate) fn dialog_gcd_clean_truncated_underflow(
    b: &mut B,
    acc: &[QubitId],
    a: &[QubitId],
    ctrl: QubitId,
    acc_ovf: QubitId,
    step: Option<usize>,
) {
    dialog_gcd_clean_truncated_underflow_with_borrowed(
        b,
        acc,
        a,
        ctrl,
        acc_ovf,
        step,
        None,
    );
}

pub(crate) fn dialog_gcd_special_underflow_clean_compare_bits(step: Option<usize>) -> usize {
    dialog_gcd_special_clean_compare_bits_from_env(
        step,
        "DIALOG_GCD_SPECIAL_UNDERFLOW_CLEAN_STEP_BITS",
    )
}

pub(crate) fn dialog_gcd_special_overflow_clean_compare_bits(step: Option<usize>) -> usize {
    dialog_gcd_special_clean_compare_bits_from_env(
        step,
        "DIALOG_GCD_SPECIAL_OVERFLOW_CLEAN_STEP_BITS",
    )
}

pub(crate) fn dialog_gcd_special_fold_carry_trunc_window(
    step: Option<usize>,
) -> Option<usize> {
    step.and_then(|step| {
        dialog_gcd_step_map_override(
            "DIALOG_GCD_SPECIAL_FOLD_CARRY_TRUNC_STEP_WINDOWS",
            step,
        )
    })
    .filter(|&window| window > 0)
    .or_else(fold_carry_trunc_window)
}

pub(crate) fn dialog_gcd_special_clean_compare_bits_from_env(
    step: Option<usize>,
    env_name: &str,
) -> usize {
    let default_bits = dialog_gcd_apply_clean_compare_bits();
    let Some(step) = step else {
        return default_bits;
    };
    let Ok(spec) = std::env::var(env_name) else {
        return default_bits;
    };
    for item in spec.split(',') {
        let Some((raw_step, raw_bits)) = item.trim().split_once(':') else {
            continue;
        };
        if raw_step.trim().parse::<usize>().ok() != Some(step) {
            continue;
        }
        if let Ok(bits) = raw_bits.trim().parse::<usize>() {
            if (1..=N).contains(&bits) {
                return bits;
            }
        }
    }
    default_bits
}

pub(crate) fn dialog_gcd_load_controlled_slice(
    b: &mut B,
    ctrl: QubitId,
    source: &[QubitId],
    lo: usize,
    hi: usize,
) -> Vec<QubitId> {
    assert!(lo <= hi);
    assert!(hi <= source.len());
    let out = b.alloc_qubits(hi - lo);
    for (i, &q) in source[lo..hi].iter().enumerate() {
        b.ccx(ctrl, q, out[i]);
    }
    out
}

pub(crate) fn dialog_gcd_clear_controlled_slice_hmr(
    b: &mut B,
    ctrl: QubitId,
    source: &[QubitId],
    lo: usize,
    loaded: &[QubitId],
) {
    assert!(lo + loaded.len() <= source.len());
    for (i, &q) in loaded.iter().enumerate() {
        let m = b.alloc_bit();
        b.hmr(q, m);
        b.cz_if(ctrl, source[lo + i], m);
    }
}

pub(crate) fn dialog_gcd_chunk_hi(blocks: usize, block: usize, ext_n: usize) -> usize {
    if let Some(cuts) = dialog_gcd_apply_chunked_f_cuts() {
        assert_eq!(
            cuts.len() + 1,
            blocks,
            "DIALOG_GCD_APPLY_CHUNKED_F_CUTS must contain blocks-1 cuts"
        );
        assert!(
            cuts.first().is_some_and(|&cut| cut > 0)
                && cuts.windows(2).all(|pair| pair[0] < pair[1])
                && cuts.last().is_some_and(|&cut| cut < ext_n),
            "DIALOG_GCD_APPLY_CHUNKED_F_CUTS must be strictly increasing in 1..{ext_n}: {cuts:?}"
        );
        if block < cuts.len() {
            return cuts[block];
        }
    }
    if blocks == 4 && dialog_gcd_apply_chunked_f_custom4_enabled() {
        let cuts = [
            dialog_gcd_apply_chunked_f_cut().unwrap_or(ext_n / 4),
            dialog_gcd_apply_chunked_f_cut2().unwrap_or(ext_n / 2),
            dialog_gcd_apply_chunked_f_cut3().unwrap_or(3 * ext_n / 4),
        ];
        assert!(
            cuts[0] < cuts[1] && cuts[1] < cuts[2] && cuts[2] < ext_n,
            "custom four-chunk apply boundaries must be strictly increasing and below {ext_n}: {cuts:?}"
        );
        if block < cuts.len() {
            return cuts[block];
        }
    }
    if blocks == 5 && dialog_gcd_apply_chunked_f_custom5_enabled() {
        let cuts = [
            dialog_gcd_apply_chunked_f_cut().unwrap_or(ext_n / 5),
            dialog_gcd_apply_chunked_f_cut2().unwrap_or((2 * ext_n) / 5),
            dialog_gcd_apply_chunked_f_cut3().unwrap_or((3 * ext_n) / 5),
            dialog_gcd_apply_chunked_f_cut4().unwrap_or((4 * ext_n) / 5),
        ];
        assert!(
            cuts[0] < cuts[1] && cuts[1] < cuts[2] && cuts[2] < cuts[3] && cuts[3] < ext_n,
            "custom five-chunk apply boundaries must be strictly increasing and below {ext_n}: {cuts:?}"
        );
        if block < cuts.len() {
            return cuts[block];
        }
    }
    if block == 0 && blocks <= 3 {
        return dialog_gcd_apply_chunked_f_cut()
            .unwrap_or(ext_n / 2)
            .min(ext_n - 1);
    }
    if blocks == 3 && block == 1 {
        return dialog_gcd_apply_chunked_f_cut2()
            .unwrap_or(2 * ext_n / 3)
            .min(ext_n - 1);
    }
    ((block + 1) * ext_n) / blocks
}

fn dialog_gcd_add_fast_with_borrowed_carries(
    b: &mut B,
    a: &[QubitId],
    acc: &[QubitId],
    c_in: QubitId,
    borrowed: &[QubitId],
) {
    let needed = a.len().saturating_sub(1);
    let borrowed = &borrowed[..borrowed.len().min(needed)];
    let owned = b.alloc_qubits(needed - borrowed.len());
    let mut carries = Vec::with_capacity(needed);
    carries.extend_from_slice(borrowed);
    carries.extend_from_slice(&owned);
    cuccaro_add_fast_borrowed_carries(b, a, acc, c_in, &carries);
    b.free_vec(&owned);
}

fn dialog_gcd_sub_fast_with_borrowed_carries(
    b: &mut B,
    a: &[QubitId],
    acc: &[QubitId],
    c_in: QubitId,
    borrowed: &[QubitId],
) {
    let needed = a.len().saturating_sub(1);
    let borrowed = &borrowed[..borrowed.len().min(needed)];
    let owned = b.alloc_qubits(needed - borrowed.len());
    let mut carries = Vec::with_capacity(needed);
    carries.extend_from_slice(borrowed);
    carries.extend_from_slice(&owned);
    cuccaro_sub_fast_borrowed_carries(b, a, acc, c_in, &carries);
    b.free_vec(&owned);
}

fn dialog_gcd_add_fast_low_to_ext_with_borrowed_carries(
    b: &mut B,
    a: &[QubitId],
    acc_ext: &[QubitId],
    c_in: QubitId,
    borrowed: &[QubitId],
) {
    let needed = a.len();
    let borrowed = &borrowed[..borrowed.len().min(needed)];
    let owned = b.alloc_qubits(needed - borrowed.len());
    let mut carries = Vec::with_capacity(needed);
    carries.extend_from_slice(borrowed);
    carries.extend_from_slice(&owned);
    cuccaro_add_fast_low_to_ext_borrowed_carries(b, a, acc_ext, c_in, &carries);
    b.free_vec(&owned);
}

fn dialog_gcd_add_fast_low_to_ext_with_borrowed_carries_topclean(
    b: &mut B,
    a: &[QubitId],
    acc_ext: &[QubitId],
    c_in: QubitId,
    borrowed_carries: &[QubitId],
    clean_top: usize,
) {
    let clean_top = clean_top.min(a.len().saturating_sub(1));
    if clean_top == 0 {
        return dialog_gcd_add_fast_low_to_ext_with_borrowed_carries(
            b,
            a,
            acc_ext,
            c_in,
            borrowed_carries,
        );
    }
    let needed_carries = a.len() - clean_top;
    let borrowed = borrowed_carries.len().min(needed_carries);
    let owned = b.alloc_qubits(needed_carries - borrowed);
    let mut carries = Vec::with_capacity(needed_carries);
    carries.extend_from_slice(&borrowed_carries[..borrowed]);
    carries.extend_from_slice(&owned);
    cuccaro_add_fast_low_to_ext_borrowed_carries_topclean(
        b,
        a,
        acc_ext,
        c_in,
        &carries,
        clean_top,
    );
    b.free_vec(&owned);
}

fn dialog_gcd_sub_fast_low_to_ext_with_borrowed_carries(
    b: &mut B,
    a: &[QubitId],
    acc_ext: &[QubitId],
    c_in: QubitId,
    borrowed: &[QubitId],
) {
    let needed = a.len();
    let borrowed = &borrowed[..borrowed.len().min(needed)];
    let owned = b.alloc_qubits(needed - borrowed.len());
    let mut carries = Vec::with_capacity(needed);
    carries.extend_from_slice(borrowed);
    carries.extend_from_slice(&owned);
    cuccaro_sub_fast_low_to_ext_borrowed_carries(b, a, acc_ext, c_in, &carries);
    b.free_vec(&owned);
}

fn dialog_gcd_sub_fast_low_to_ext_with_borrowed_carries_topclean(
    b: &mut B,
    a: &[QubitId],
    acc_ext: &[QubitId],
    c_in: QubitId,
    borrowed_carries: &[QubitId],
    clean_top: usize,
) {
    let clean_top = clean_top.min(a.len().saturating_sub(1));
    if clean_top == 0 {
        return dialog_gcd_sub_fast_low_to_ext_with_borrowed_carries(
            b,
            a,
            acc_ext,
            c_in,
            borrowed_carries,
        );
    }
    let needed_carries = a.len() - clean_top;
    let borrowed = borrowed_carries.len().min(needed_carries);
    let owned = b.alloc_qubits(needed_carries - borrowed);
    let mut carries = Vec::with_capacity(needed_carries);
    carries.extend_from_slice(&borrowed_carries[..borrowed]);
    carries.extend_from_slice(&owned);
    cuccaro_sub_fast_low_to_ext_borrowed_carries_topclean(
        b,
        a,
        acc_ext,
        c_in,
        &carries,
        clean_top,
    );
    b.free_vec(&owned);
}

fn dialog_gcd_conditional_boundary_replay(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    ctrl: QubitId,
    c_in: QubitId,
    targets: &[(QubitId, usize)],
) {
    assert!(!targets.is_empty());
    assert!(targets.windows(2).all(|w| w[0].1 < w[1].1));
    for index in (0..targets.len()).rev() {
        let (target, p) = targets[index];
        let (start, carry_in) = if index == 0 {
            (0, c_in)
        } else {
            (targets[index - 1].1, targets[index - 1].0)
        };
        let phase = b.alloc_bit();
        b.hmr(target, phase);
        cmp_lt_phase_conditioned_with_cin(
            b,
            &u[start..p],
            &v[start..p],
            carry_in,
            ctrl,
            phase,
        );
    }
}

fn dialog_gcd_conditional_boundary_replay_free_owned(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    ctrl: QubitId,
    c_in: QubitId,
    targets: &[(QubitId, usize, bool)],
) {
    assert!(!targets.is_empty());
    assert!(targets.windows(2).all(|w| w[0].1 < w[1].1));
    for index in (0..targets.len()).rev() {
        let (target, p, owned_target) = targets[index];
        let (start, carry_in) = if index == 0 {
            (0, c_in)
        } else {
            (targets[index - 1].1, targets[index - 1].0)
        };
        let phase = b.alloc_bit();
        b.hmr(target, phase);
        if owned_target {
            b.free(target);
        }
        cmp_lt_phase_conditioned_with_cin(
            b,
            &u[start..p],
            &v[start..p],
            carry_in,
            ctrl,
            phase,
        );
    }
}

fn dialog_gcd_apply_auto_topclean_bits(
    active_before_ripple: u32,
    source_len: usize,
    future_boundary_carries: &[QubitId],
) -> usize {
    let Some(target) = dialog_gcd_apply_chunked_f_auto_topclean_target() else {
        return 0;
    };
    if source_len <= 1 {
        return 0;
    }
    let future_borrowed = future_boundary_carries.len().min(source_len);
    let owned_carries_without_topclean = source_len - future_borrowed;
    let projected_peak = active_before_ripple as usize + owned_carries_without_topclean;
    let needed = projected_peak.saturating_sub(target as usize);
    needed
        .min(dialog_gcd_apply_chunked_f_auto_topclean_max_bits())
        .min(source_len - 1)
}

pub(crate) fn dialog_gcd_add_ctrl_chunked_low_to_ext(
    b: &mut B,
    source: &[QubitId],
    acc_ext: &[QubitId],
    ctrl: QubitId,
    c_in: QubitId,
    blocks: usize,
    clean_scratch: &[QubitId],
) {
    let n = source.len();
    assert_eq!(acc_ext.len(), n + 1);
    for (i, &q) in clean_scratch.iter().enumerate() {
        assert!(!clean_scratch[..i].contains(&q));
        assert!(!source.contains(&q));
        assert!(!acc_ext.contains(&q));
        assert_ne!(q, ctrl);
        assert_ne!(q, c_in);
    }
    let ext_n = acc_ext.len();
    let blocks = blocks.max(2).min(ext_n);
    let mut carry = c_in;
    let mut lo = 0usize;
    // The low-to-extended-register primitive represents the source high zero
    // implicitly. Otherwise reserve one borrowed cell for that transient lane.
    let implicit_high_zero = dialog_gcd_apply_implicit_high_zero_enabled();
    let zero_host = (!implicit_high_zero)
        .then(|| clean_scratch.first().copied())
        .flatten();
    let boundary_hosts = &clean_scratch
        [usize::from(!implicit_high_zero && zero_host.is_some())..];
    let mut couts: Vec<(QubitId, usize, bool)> = Vec::new();

    for blk in 0..blocks {
        let hi = dialog_gcd_chunk_hi(blocks, blk, ext_n);
        if hi <= lo {
            continue;
        }
        if blk == blocks - 1 || hi == ext_n {
            b.set_phase("dialog_gcd_apply_chunk_add_final_load");
            let f = dialog_gcd_load_controlled_slice(b, ctrl, source, lo.min(n), n);
            b.set_phase("dialog_gcd_apply_chunk_add_final_ripple");
            let final_topclean = dialog_gcd_apply_final_topclean_bits()
                .max(dialog_gcd_apply_auto_topclean_bits(b.active_qubits, f.len(), &[]));
            if final_topclean > 0 {
                cuccaro_add_fast_low_to_ext_topclean(b, &f, &acc_ext[lo..hi], carry, final_topclean);
            } else if let Some(window_blocks) = dialog_gcd_apply_final_windowed_fast_blocks() {
                cuccaro_add_fast_windowed_low_to_ext(
                    b,
                    &f,
                    &acc_ext[lo..hi],
                    carry,
                    window_blocks,
                );
            } else if dialog_gcd_apply_final_lowq_enabled() {
                let zero = b.alloc_qubit();
                let mut f_ext = f.clone();
                f_ext.push(zero);
                cuccaro_add(b, &f_ext, &acc_ext[lo..hi], carry);
                b.free(zero);
            } else {
                cuccaro_add_fast_low_to_ext(b, &f, &acc_ext[lo..hi], carry);
            }
            b.set_phase("dialog_gcd_apply_chunk_add_final_clear");
            dialog_gcd_clear_controlled_slice_hmr(b, ctrl, source, lo.min(n), &f);
            b.free_vec(&f);
            break;
        }

        assert!(hi <= n);
        b.set_phase("dialog_gcd_apply_chunk_add_load");
        let f = dialog_gcd_load_controlled_slice(b, ctrl, source, lo, hi);
        let needs_distinct_zero =
            carry == c_in || !dialog_gcd_apply_chunked_f_reuse_cin_zero_enabled();
        let (zero, owned_zero) = if needs_distinct_zero {
            zero_host.map_or_else(|| (b.alloc_qubit(), true), |q| (q, false))
        } else {
            (c_in, false)
        };
        let (cout, owned_cout) = boundary_hosts
            .get(couts.len())
            .copied()
            .map_or_else(|| (b.alloc_qubit(), true), |q| (q, false));
        let mut acc_block = acc_ext[lo..hi].to_vec();
        acc_block.push(cout);
        let future_boundary_carries = if dialog_gcd_apply_borrow_future_boundary_carries_enabled() {
            boundary_hosts.get(couts.len() + 1..).unwrap_or(&[])
        } else {
            &[]
        };
        let topclean_bits = if implicit_high_zero {
            dialog_gcd_apply_auto_topclean_bits(b.active_qubits, f.len(), future_boundary_carries)
        } else {
            0
        };
        b.set_phase("dialog_gcd_apply_chunk_add_ripple");
        if implicit_high_zero {
            if topclean_bits > 0 {
                dialog_gcd_add_fast_low_to_ext_with_borrowed_carries_topclean(
                    b,
                    &f,
                    &acc_block,
                    carry,
                    future_boundary_carries,
                    topclean_bits,
                );
            } else {
                dialog_gcd_add_fast_low_to_ext_with_borrowed_carries(
                    b,
                    &f,
                    &acc_block,
                    carry,
                    future_boundary_carries,
                );
            }
        } else {
            let mut a_block = f.clone();
            a_block.push(zero);
            dialog_gcd_add_fast_with_borrowed_carries(
                b,
                &a_block,
                &acc_block,
                carry,
                future_boundary_carries,
            );
        }
        if owned_zero {
            b.free(zero);
        }
        b.set_phase("dialog_gcd_apply_chunk_add_clear");
        dialog_gcd_clear_controlled_slice_hmr(b, ctrl, source, lo, &f);
        b.free_vec(&f);
        couts.push((cout, hi, owned_cout));
        carry = cout;
        lo = hi;
    }

    let mut boundary_replay_freed_owned = false;
    if dialog_gcd_apply_chunked_f_fuse_boundary_clears_enabled() {
        if let Some(&(_, p, _)) = couts.last() {
            b.set_phase("dialog_gcd_apply_chunk_add_boundary_clear");
            let targets = couts
                .iter()
                .map(|&(cout, p, _)| (cout, p))
                .collect::<Vec<_>>();
            if dialog_gcd_apply_boundary_conditional_replay_enabled() {
                if dialog_gcd_apply_boundary_free_owned_during_replay_enabled() {
                    dialog_gcd_conditional_boundary_replay_free_owned(
                        b,
                        &acc_ext[..p],
                        &source[..p],
                        ctrl,
                        c_in,
                        &couts,
                    );
                    boundary_replay_freed_owned = true;
                } else {
                    dialog_gcd_conditional_boundary_replay(
                        b,
                        &acc_ext[..p],
                        &source[..p],
                        ctrl,
                        c_in,
                        &targets,
                    );
                }
            } else if let Some(split) = dialog_gcd_apply_boundary_split() {
                ccx_cmp_lt_into_fast_prefix_targets_split(
                    b,
                    &acc_ext[..p],
                    &source[..p],
                    ctrl,
                    &targets,
                    split.min(p.saturating_sub(1)),
                );
            } else {
                ccx_cmp_lt_into_fast_prefix_targets(b, &acc_ext[..p], &source[..p], ctrl, &targets);
            }
        }
    } else {
        for &(cout, p, _) in couts.iter().rev() {
            b.set_phase("dialog_gcd_apply_chunk_add_boundary_clear");
            ccx_cmp_lt_into_fast(b, &acc_ext[..p], &source[..p], ctrl, cout);
        }
    }
    for &(cout, _, owned_cout) in couts.iter().rev() {
        if owned_cout && !boundary_replay_freed_owned {
            b.free(cout);
        }
    }
}

pub(crate) fn dialog_gcd_sub_ctrl_chunked_low_to_ext(
    b: &mut B,
    source: &[QubitId],
    acc_ext: &[QubitId],
    ctrl: QubitId,
    c_in: QubitId,
    blocks: usize,
    clean_scratch: &[QubitId],
) {
    let n = source.len();
    assert_eq!(acc_ext.len(), n + 1);
    for (i, &q) in clean_scratch.iter().enumerate() {
        assert!(!clean_scratch[..i].contains(&q));
        assert!(!source.contains(&q));
        assert!(!acc_ext.contains(&q));
        assert_ne!(q, ctrl);
        assert_ne!(q, c_in);
    }
    let ext_n = acc_ext.len();
    let blocks = blocks.max(2).min(ext_n);
    let mut borrow = c_in;
    let mut lo = 0usize;
    let implicit_high_zero = dialog_gcd_apply_implicit_high_zero_enabled();
    let zero_host = (!implicit_high_zero)
        .then(|| clean_scratch.first().copied())
        .flatten();
    let boundary_hosts = &clean_scratch
        [usize::from(!implicit_high_zero && zero_host.is_some())..];
    let mut bouts: Vec<(QubitId, usize, bool)> = Vec::new();

    for blk in 0..blocks {
        let hi = dialog_gcd_chunk_hi(blocks, blk, ext_n);
        if hi <= lo {
            continue;
        }
        if blk == blocks - 1 || hi == ext_n {
            b.set_phase("dialog_gcd_apply_chunk_sub_final_load");
            let f = dialog_gcd_load_controlled_slice(b, ctrl, source, lo.min(n), n);
            b.set_phase("dialog_gcd_apply_chunk_sub_final_ripple");
            let final_topclean = dialog_gcd_apply_final_topclean_bits()
                .max(dialog_gcd_apply_auto_topclean_bits(b.active_qubits, f.len(), &[]));
            if final_topclean > 0 {
                cuccaro_sub_fast_low_to_ext_topclean(b, &f, &acc_ext[lo..hi], borrow, final_topclean);
            } else if let Some(window_blocks) = dialog_gcd_apply_final_windowed_fast_blocks() {
                cuccaro_sub_fast_windowed_low_to_ext(
                    b,
                    &f,
                    &acc_ext[lo..hi],
                    borrow,
                    window_blocks,
                );
            } else if dialog_gcd_apply_final_lowq_enabled() {
                let zero = b.alloc_qubit();
                let mut f_ext = f.clone();
                f_ext.push(zero);
                cuccaro_sub(b, &f_ext, &acc_ext[lo..hi], borrow);
                b.free(zero);
            } else {
                cuccaro_sub_fast_low_to_ext(b, &f, &acc_ext[lo..hi], borrow);
            }
            b.set_phase("dialog_gcd_apply_chunk_sub_final_clear");
            dialog_gcd_clear_controlled_slice_hmr(b, ctrl, source, lo.min(n), &f);
            b.free_vec(&f);
            break;
        }

        assert!(hi <= n);
        b.set_phase("dialog_gcd_apply_chunk_sub_load");
        let f = dialog_gcd_load_controlled_slice(b, ctrl, source, lo, hi);
        let needs_distinct_zero =
            borrow == c_in || !dialog_gcd_apply_chunked_f_reuse_cin_zero_enabled();
        let (zero, owned_zero) = if needs_distinct_zero {
            zero_host.map_or_else(|| (b.alloc_qubit(), true), |q| (q, false))
        } else {
            (c_in, false)
        };
        let (bout, owned_bout) = boundary_hosts
            .get(bouts.len())
            .copied()
            .map_or_else(|| (b.alloc_qubit(), true), |q| (q, false));
        let mut acc_block = acc_ext[lo..hi].to_vec();
        acc_block.push(bout);
        let future_boundary_carries = if dialog_gcd_apply_borrow_future_boundary_carries_enabled() {
            boundary_hosts.get(bouts.len() + 1..).unwrap_or(&[])
        } else {
            &[]
        };
        let topclean_bits = if implicit_high_zero {
            dialog_gcd_apply_auto_topclean_bits(b.active_qubits, f.len(), future_boundary_carries)
        } else {
            0
        };
        b.set_phase("dialog_gcd_apply_chunk_sub_ripple");
        if implicit_high_zero {
            if topclean_bits > 0 {
                dialog_gcd_sub_fast_low_to_ext_with_borrowed_carries_topclean(
                    b,
                    &f,
                    &acc_block,
                    borrow,
                    future_boundary_carries,
                    topclean_bits,
                );
            } else {
                dialog_gcd_sub_fast_low_to_ext_with_borrowed_carries(
                    b,
                    &f,
                    &acc_block,
                    borrow,
                    future_boundary_carries,
                );
            }
        } else {
            let mut a_block = f.clone();
            a_block.push(zero);
            dialog_gcd_sub_fast_with_borrowed_carries(
                b,
                &a_block,
                &acc_block,
                borrow,
                future_boundary_carries,
            );
        }
        if owned_zero {
            b.free(zero);
        }
        b.set_phase("dialog_gcd_apply_chunk_sub_clear");
        dialog_gcd_clear_controlled_slice_hmr(b, ctrl, source, lo, &f);
        b.free_vec(&f);
        bouts.push((bout, hi, owned_bout));
        borrow = bout;
        lo = hi;
    }

    let mut boundary_replay_freed_owned = false;
    if dialog_gcd_apply_chunked_f_fuse_boundary_clears_enabled() {
        if let Some(&(_, p, _)) = bouts.last() {
            b.set_phase("dialog_gcd_apply_chunk_sub_boundary_clear");
            for i in 0..p {
                b.x(source[i]);
            }
            let targets = bouts
                .iter()
                .map(|&(bout, p, _)| (bout, p))
                .collect::<Vec<_>>();
            if dialog_gcd_apply_boundary_conditional_replay_enabled() {
                if dialog_gcd_apply_boundary_free_owned_during_replay_enabled() {
                    dialog_gcd_conditional_boundary_replay_free_owned(
                        b,
                        &source[..p],
                        &acc_ext[..p],
                        ctrl,
                        c_in,
                        &bouts,
                    );
                    boundary_replay_freed_owned = true;
                } else {
                    dialog_gcd_conditional_boundary_replay(
                        b,
                        &source[..p],
                        &acc_ext[..p],
                        ctrl,
                        c_in,
                        &targets,
                    );
                }
            } else if let Some(split) = dialog_gcd_apply_boundary_split() {
                ccx_cmp_lt_into_fast_prefix_targets_split(
                    b,
                    &source[..p],
                    &acc_ext[..p],
                    ctrl,
                    &targets,
                    split.min(p.saturating_sub(1)),
                );
            } else {
                ccx_cmp_lt_into_fast_prefix_targets(b, &source[..p], &acc_ext[..p], ctrl, &targets);
            }
            for i in 0..p {
                b.x(source[i]);
            }
        }
    } else {
        for &(bout, p, _) in bouts.iter().rev() {
            b.set_phase("dialog_gcd_apply_chunk_sub_boundary_clear");
            for i in 0..p {
                b.x(source[i]);
            }
            ccx_cmp_lt_into_fast(b, &source[..p], &acc_ext[..p], ctrl, bout);
            for i in 0..p {
                b.x(source[i]);
            }
        }
    }
    for &(bout, _, owned_bout) in bouts.iter().rev() {
        if owned_bout && !boundary_replay_freed_owned {
            b.free(bout);
        }
    }
}

pub(crate) fn dialog_gcd_cmod_add_materialized_pseudomersenne_chunked(
    b: &mut B,
    acc: &[QubitId],
    a: &[QubitId],
    ctrl: QubitId,
    p: U256,
    blocks: usize,
    clean_scratch: &[QubitId],
    step: Option<usize>,
) {
    assert_eq!(acc.len(), N);
    assert_eq!(a.len(), N);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1u64));

    let (acc_ext, acc_ovf) = ext_reg(b, acc);
    for (i, &q) in clean_scratch.iter().enumerate() {
        assert!(!clean_scratch[..i].contains(&q));
        assert!(!acc_ext.contains(&q));
        assert!(!a.contains(&q));
        assert_ne!(q, ctrl);
    }
    let (c_in, owned_c_in, inner_scratch) = clean_scratch.split_first().map_or_else(
        || (b.alloc_qubit(), true, &[][..]),
        |(&q, rest)| (q, false, rest),
    );

    b.set_phase("dialog_gcd_materialized_special_chunked_raw_sum");
    dialog_gcd_add_ctrl_chunked_low_to_ext(b, a, &acc_ext, ctrl, c_in, blocks, inner_scratch);
    if owned_c_in {
        b.free(c_in);
    }

    b.set_phase("dialog_gcd_materialized_special_overflow_fold");
    if let Some(w) = dialog_gcd_special_fold_carry_trunc_window(step) {
        let borrowed_carries = if std::env::var("DIALOG_GCD_SPECIAL_FOLD_BORROW_CARRIES")
            .ok()
            .as_deref()
            == Some("1")
        {
            // The chunk carry-in is back to |0> after the raw sum and is idle
            // during the fold, so it is a valid carry host alongside the
            // remaining clean scratch.
            clean_scratch
        } else {
            &[]
        };
        if std::env::var("DIALOG_GCD_SPECIAL_FOLD_RELEASE_SCRATCH")
            .ok()
            .as_deref()
            == Some("1")
            && !borrowed_carries.is_empty()
        {
            assert_eq!(
                std::env::var("DIALOG_GCD_K2_APPLY_INPLACE_RAW_BLOCK")
                    .ok()
                    .as_deref(),
                Some("1"),
                "special-fold scratch release requires owned in-place apply scratch"
            );
            cadd_nbit_const_direct_trunc_fast_releasing_scratch_at_step(
                b,
                &acc[..DIALOG_GCD_SPECIAL_ADD_LSBS],
                c,
                acc_ovf,
                w,
                borrowed_carries,
                step,
            );
        } else {
            cadd_nbit_const_direct_trunc_fast_borrowed_carries(
                b,
                &acc[..DIALOG_GCD_SPECIAL_ADD_LSBS],
                c,
                acc_ovf,
                w,
                borrowed_carries,
            );
        }
    } else {
        cadd_nbit_const_fast(b, &acc[..DIALOG_GCD_SPECIAL_ADD_LSBS], c, acc_ovf);
    }

    b.set_phase("dialog_gcd_materialized_special_overflow_clean");
    let compare_start = N - dialog_gcd_special_overflow_clean_compare_bits(step);
    if dialog_gcd_special_clean_conditional_replay_enabled() {
        let phase = b.alloc_bit();
        b.hmr(acc_ovf, phase);
        dialog_gcd_cmp_lt_phase_conditioned_hosted(
            b,
            &acc[compare_start..],
            &a[compare_start..],
            ctrl,
            phase,
            Some(clean_scratch),
        );
    } else {
        ccx_cmp_lt_into_fast(b, &acc[compare_start..], &a[compare_start..], ctrl, acc_ovf);
    }
    unext_reg(b, acc_ovf);
}

pub(crate) fn dialog_gcd_cmod_sub_materialized_pseudomersenne_chunked(
    b: &mut B,
    acc: &[QubitId],
    a: &[QubitId],
    ctrl: QubitId,
    p: U256,
    blocks: usize,
    clean_scratch: &[QubitId],
    step: Option<usize>,
) {
    assert_eq!(acc.len(), N);
    assert_eq!(a.len(), N);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1u64));

    let (acc_ext, acc_ovf) = ext_reg(b, acc);
    for (i, &q) in clean_scratch.iter().enumerate() {
        assert!(!clean_scratch[..i].contains(&q));
        assert!(!acc_ext.contains(&q));
        assert!(!a.contains(&q));
        assert_ne!(q, ctrl);
    }
    let (c_in, owned_c_in, inner_scratch) = clean_scratch.split_first().map_or_else(
        || (b.alloc_qubit(), true, &[][..]),
        |(&q, rest)| (q, false, rest),
    );

    b.set_phase("dialog_gcd_materialized_special_chunked_raw_difference");
    dialog_gcd_sub_ctrl_chunked_low_to_ext(b, a, &acc_ext, ctrl, c_in, blocks, inner_scratch);
    if owned_c_in {
        b.free(c_in);
    }

    b.set_phase("dialog_gcd_materialized_special_underflow_fold");
    if let Some(w) = dialog_gcd_special_fold_carry_trunc_window(step) {
        let borrowed_carries = if std::env::var("DIALOG_GCD_SPECIAL_FOLD_BORROW_CARRIES")
            .ok()
            .as_deref()
            == Some("1")
        {
            // The chunk borrow-in is back to |0> after the raw difference and
            // can host one fold borrow without increasing the live set.
            clean_scratch
        } else {
            &[]
        };
        if std::env::var("DIALOG_GCD_SPECIAL_FOLD_RELEASE_SCRATCH")
            .ok()
            .as_deref()
            == Some("1")
            && !borrowed_carries.is_empty()
        {
            assert_eq!(
                std::env::var("DIALOG_GCD_K2_APPLY_INPLACE_RAW_BLOCK")
                    .ok()
                    .as_deref(),
                Some("1"),
                "special-fold scratch release requires owned in-place apply scratch"
            );
            csub_nbit_const_direct_trunc_fast_releasing_scratch_at_step(
                b,
                &acc[..DIALOG_GCD_SPECIAL_ADD_LSBS],
                c,
                acc_ovf,
                w,
                borrowed_carries,
                step,
            );
        } else {
            csub_nbit_const_direct_trunc_fast_borrowed_carries(
                b,
                &acc[..DIALOG_GCD_SPECIAL_ADD_LSBS],
                c,
                acc_ovf,
                w,
                borrowed_carries,
            );
        }
    } else {
        csub_nbit_const_fast(b, &acc[..DIALOG_GCD_SPECIAL_ADD_LSBS], c, acc_ovf);
    }

    b.set_phase("dialog_gcd_materialized_special_underflow_clean");
    dialog_gcd_clean_truncated_underflow_with_borrowed(
        b,
        acc,
        a,
        ctrl,
        acc_ovf,
        step,
        Some(clean_scratch),
    );
    unext_reg(b, acc_ovf);
}

pub(crate) fn dialog_gcd_cmod_sub_materialized_pseudomersenne(
    b: &mut B,
    acc: &[QubitId],
    a: &[QubitId],
    ctrl: QubitId,
    p: U256,
) {
    dialog_gcd_cmod_sub_materialized_pseudomersenne_at_step(b, acc, a, ctrl, p, None);
}

pub(crate) fn dialog_gcd_cmod_sub_materialized_pseudomersenne_at_step(
    b: &mut B,
    acc: &[QubitId],
    a: &[QubitId],
    ctrl: QubitId,
    p: U256,
    step: Option<usize>,
) {
    dialog_gcd_cmod_sub_materialized_pseudomersenne_with_clean_scratch_at_step(
        b,
        acc,
        a,
        ctrl,
        p,
        &[],
        step,
    );
}

pub(crate) fn dialog_gcd_cmod_sub_materialized_pseudomersenne_with_clean_scratch(
    b: &mut B,
    acc: &[QubitId],
    a: &[QubitId],
    ctrl: QubitId,
    p: U256,
    clean_scratch: &[QubitId],
) {
    dialog_gcd_cmod_sub_materialized_pseudomersenne_with_clean_scratch_at_step(
        b,
        acc,
        a,
        ctrl,
        p,
        clean_scratch,
        None,
    );
}

pub(crate) fn dialog_gcd_cmod_sub_materialized_pseudomersenne_with_clean_scratch_at_step(
    b: &mut B,
    acc: &[QubitId],
    a: &[QubitId],
    ctrl: QubitId,
    p: U256,
    clean_scratch: &[QubitId],
    step: Option<usize>,
) {
    assert_eq!(acc.len(), N);
    assert_eq!(a.len(), N);
    if let Some(blocks) = dialog_gcd_apply_chunked_f_blocks()
        .filter(|_| dialog_gcd_raw_apply_truncated_clean_enabled())
        .filter(|_| dialog_gcd_measured_apply_sub_enabled())
    {
        dialog_gcd_cmod_sub_materialized_pseudomersenne_chunked(
            b,
            acc,
            a,
            ctrl,
            p,
            blocks,
            clean_scratch,
            step,
        );
        return;
    }
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1u64));

    let f = b.alloc_qubits(N);
    b.set_phase("dialog_gcd_materialized_special_load_subtrahend");
    for i in 0..N {
        b.ccx(ctrl, a[i], f[i]);
    }

    let (acc_ext, acc_ovf) = ext_reg(b, acc);

    b.set_phase("dialog_gcd_materialized_special_raw_difference");
    if dialog_gcd_measured_apply_sub_enabled() {
        // Measured (Gidney) difference: ~n Toffoli instead of the ~2n of the
        // non-fast cuccaro_sub uncompute. Peak-safe: the symmetric apply ADD
        // already runs cuccaro_add_fast with its carry lane in this same phase.
        let c_in = b.alloc_qubit();
        if let Some(w) = dialog_gcd_apply_window_blocks() {
            cuccaro_sub_fast_windowed_low_to_ext(b, &f, &acc_ext, c_in, w);
        } else {
            let f_ovf = b.alloc_qubit();
            let mut f_ext = f.clone();
            f_ext.push(f_ovf);
            cuccaro_sub_fast(b, &f_ext, &acc_ext, c_in);
            b.free(f_ovf);
        }
        b.free(c_in);
    } else {
        let f_ovf = b.alloc_qubit();
        let mut f_ext = f.clone();
        f_ext.push(f_ovf);
        sub_nbit_qq(b, &f_ext, &acc_ext);
        b.free(f_ovf);
    }

    b.set_phase("dialog_gcd_materialized_special_underflow_fold");
    if let Some(w) = dialog_gcd_special_fold_carry_trunc_window(step) {
        csub_nbit_const_direct_trunc_fast(b, &acc[..DIALOG_GCD_SPECIAL_ADD_LSBS], c, acc_ovf, w);
    } else {
        csub_nbit_const_fast(b, &acc[..DIALOG_GCD_SPECIAL_ADD_LSBS], c, acc_ovf);
    }

    b.set_phase("dialog_gcd_materialized_special_underflow_clean");
    if dialog_gcd_raw_apply_truncated_clean_enabled() {
        dialog_gcd_clean_truncated_underflow(b, acc, a, ctrl, acc_ovf, step);
    } else {
        b.x(acc_ovf);
        mod_neg_inplace_fast(b, &f, p);
        cmp_lt_into_fast(b, acc, &f, acc_ovf);
        mod_neg_inplace_fast(b, &f, p);
    }
    unext_reg(b, acc_ovf);

    b.set_phase("dialog_gcd_materialized_special_clear_subtrahend");
    for i in 0..N {
        let m = b.alloc_bit();
        b.hmr(f[i], m);
        b.cz_if(ctrl, a[i], m);
    }
    b.free_vec(&f);
}

pub(crate) fn emit_dialog_gcd_raw_apply_bitvector(
    b: &mut B,
    dialog_log: &[QubitId],
    x: &[QubitId],
    y: &[QubitId],
    p: U256,
) {
    assert!(dialog_log.len() >= 2 * dialog_gcd_active_iterations());
    assert_eq!(x.len(), N);
    assert_eq!(y.len(), N);

    for step in (0..dialog_gcd_active_iterations()).rev() {
        let b0 = dialog_log[2 * step];
        let b0_and_b1 = dialog_log[2 * step + 1];

        b.set_phase("dialog_gcd_raw_apply_double_y");
        mod_double_inplace_fast(b, y, p);

        b.set_phase("dialog_gcd_raw_apply_cadd");
        if dialog_gcd_raw_apply_materialized_special_add_enabled() {
            dialog_gcd_cmod_add_materialized_pseudomersenne_at_step(b, y, x, b0, p, Some(step));
        } else if dialog_gcd_raw_apply_direct_special_add_enabled() {
            dialog_gcd_cmod_add_pseudomersenne_lowq(b, y, x, b0, p);
        } else {
            cmod_add_qq_lowq(b, y, x, b0, p);
        }

        b.set_phase("dialog_gcd_raw_apply_cswap");
        for (&xi, &yi) in x.iter().zip(y.iter()) {
            cswap(b, b0_and_b1, xi, yi);
        }
    }
}

pub(crate) fn emit_dialog_gcd_raw_apply_bitvector_reverse_exact(
    b: &mut B,
    dialog_log: &[QubitId],
    x: &[QubitId],
    y: &[QubitId],
    p: U256,
) {
    assert!(dialog_log.len() >= 2 * dialog_gcd_active_iterations());
    assert_eq!(x.len(), N);
    assert_eq!(y.len(), N);

    for step in 0..dialog_gcd_active_iterations() {
        let b0 = dialog_log[2 * step];
        let b0_and_b1 = dialog_log[2 * step + 1];

        b.set_phase("dialog_gcd_raw_apply_reverse_cswap");
        for (&xi, &yi) in x.iter().zip(y.iter()) {
            cswap(b, b0_and_b1, xi, yi);
        }

        b.set_phase("dialog_gcd_raw_apply_reverse_csub");
        if dialog_gcd_raw_apply_reverse_materialized_special_sub_enabled() {
            dialog_gcd_cmod_sub_materialized_pseudomersenne_at_step(b, y, x, b0, p, Some(step));
        } else if dialog_gcd_raw_apply_reverse_fast_sub_enabled() {
            cmod_sub_qq(b, y, x, b0, p);
        } else {
            cmod_sub_qq_lowq(b, y, x, b0, p);
        }

        b.set_phase("dialog_gcd_raw_apply_reverse_halve_y");
        mod_halve_inplace_fast(b, y, p);
    }
}

pub(crate) fn cmod_sub_qq_lowq_borrowed_subtrahend(
    b: &mut B,
    acc: &[QubitId],
    a: &[QubitId],
    ctrl: QubitId,
    p: U256,
    f: &[QubitId],
) {
    assert_eq!(acc.len(), N);
    assert_eq!(a.len(), N);
    assert_eq!(f.len(), N);

    for i in 0..N {
        b.ccx(ctrl, a[i], f[i]);
    }
    mod_sub_qq(b, acc, f, p);
    for i in (0..N).rev() {
        b.ccx(ctrl, a[i], f[i]);
    }
}

pub(crate) fn dialog_gcd_cmod_sub_materialized_pseudomersenne_borrowed_subtrahend(
    b: &mut B,
    acc: &[QubitId],
    a: &[QubitId],
    ctrl: QubitId,
    p: U256,
    f: &[QubitId],
) {
    dialog_gcd_cmod_sub_materialized_pseudomersenne_borrowed_subtrahend_at_step(
        b,
        acc,
        a,
        ctrl,
        p,
        f,
        None,
    );
}

pub(crate) fn dialog_gcd_cmod_sub_materialized_pseudomersenne_borrowed_subtrahend_at_step(
    b: &mut B,
    acc: &[QubitId],
    a: &[QubitId],
    ctrl: QubitId,
    p: U256,
    f: &[QubitId],
    step: Option<usize>,
) {
    assert_eq!(acc.len(), N);
    assert_eq!(a.len(), N);
    assert_eq!(f.len(), N);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1u64));

    b.set_phase("dialog_gcd_materialized_special_borrowed_load_subtrahend");
    for i in 0..N {
        b.ccx(ctrl, a[i], f[i]);
    }

    let (acc_ext, acc_ovf) = ext_reg(b, acc);
    let f_ovf = b.alloc_qubit();
    let mut f_ext = f.to_vec();
    f_ext.push(f_ovf);

    b.set_phase("dialog_gcd_materialized_special_borrowed_raw_difference");
    sub_nbit_qq(b, &f_ext, &acc_ext);
    b.free(f_ovf);

    b.set_phase("dialog_gcd_materialized_special_borrowed_underflow_fold");
    if let Some(w) = dialog_gcd_special_fold_carry_trunc_window(step) {
        csub_nbit_const_direct_trunc_fast(b, &acc[..DIALOG_GCD_SPECIAL_ADD_LSBS], c, acc_ovf, w);
    } else {
        csub_nbit_const_fast(b, &acc[..DIALOG_GCD_SPECIAL_ADD_LSBS], c, acc_ovf);
    }

    b.set_phase("dialog_gcd_materialized_special_borrowed_underflow_clean");
    if dialog_gcd_raw_apply_truncated_clean_enabled() {
        dialog_gcd_clean_truncated_underflow(b, acc, a, ctrl, acc_ovf, step);
    } else {
        b.x(acc_ovf);
        mod_neg_inplace_fast(b, f, p);
        cmp_lt_into_fast(b, acc, f, acc_ovf);
        mod_neg_inplace_fast(b, f, p);
    }
    unext_reg(b, acc_ovf);

    b.set_phase("dialog_gcd_materialized_special_borrowed_clear_subtrahend");
    for i in (0..N).rev() {
        b.ccx(ctrl, a[i], f[i]);
    }
}

pub(crate) fn emit_dialog_gcd_raw_apply_bitvector_reverse_borrowed_subtrahend(
    b: &mut B,
    dialog_log: &[QubitId],
    x: &[QubitId],
    y: &[QubitId],
    p: U256,
    f: &[QubitId],
) {
    assert!(dialog_log.len() >= 2 * dialog_gcd_active_iterations());
    assert_eq!(x.len(), N);
    assert_eq!(y.len(), N);
    assert_eq!(f.len(), N);

    for step in 0..dialog_gcd_active_iterations() {
        let b0 = dialog_log[2 * step];
        let b0_and_b1 = dialog_log[2 * step + 1];

        b.set_phase("dialog_gcd_raw_apply_reverse_borrowed_cswap");
        for (&xi, &yi) in x.iter().zip(y.iter()) {
            cswap(b, b0_and_b1, xi, yi);
        }

        b.set_phase("dialog_gcd_raw_apply_reverse_borrowed_csub");
        if dialog_gcd_raw_apply_reverse_materialized_special_sub_enabled() {
            dialog_gcd_cmod_sub_materialized_pseudomersenne_borrowed_subtrahend_at_step(
                b,
                y,
                x,
                b0,
                p,
                f,
                Some(step),
            );
        } else {
            cmod_sub_qq_lowq_borrowed_subtrahend(b, y, x, b0, p, f);
        }

        b.set_phase("dialog_gcd_raw_apply_reverse_borrowed_halve_y");
        mod_halve_inplace_fast(b, y, p);
    }
}


pub(crate) fn emit_dialog_gcd_raw_ipmul(b: &mut B, factor: &[QubitId], target: &[QubitId], p: U256) {
    assert_eq!(factor.len(), N);
    assert_eq!(target.len(), N);

    if dialog_gcd_compressed_sidecar_log_enabled() {
        emit_dialog_gcd_compressed_sidecar_ipmul(b, factor, target, p);
        return;
    }

    let dialog_log = b.alloc_qubits(DIALOG_GCD_RAW_LOG_BITS);
    let u = b.alloc_qubits(N);
    b.set_phase("dialog_gcd_raw_ipmul_load_p");
    for i in 0..N {
        if bit(p, i) {
            b.x(u[i]);
        }
    }

    b.set_phase("dialog_gcd_raw_ipmul_tobitvector");
    emit_dialog_gcd_raw_tobitvector_steps(b, &u, factor, &dialog_log);

    if dialog_gcd_raw_ipmul_terminal_reuse_enabled() {
        b.set_phase("dialog_gcd_raw_ipmul_release_terminal_u");
        b.x(u[0]);
        b.free_vec(&u);

        b.set_phase("dialog_gcd_raw_ipmul_apply_bitvector_reuse_factor_zero");
        emit_dialog_gcd_raw_apply_bitvector(b, &dialog_log, target, factor, p);

        if dialog_gcd_raw_ipmul_clear_p_residual_enabled() {
            b.set_phase("dialog_gcd_raw_ipmul_clear_p_residual_source_lane");
            for i in 0..N {
                if bit(p, i) {
                    b.x(target[i]);
                }
            }
        }

        b.set_phase("dialog_gcd_raw_ipmul_swap_product_into_target");
        for i in 0..N {
            b.swap(target[i], factor[i]);
        }

        b.set_phase("dialog_gcd_raw_ipmul_reacquire_terminal_u");
        b.reacquire_vec(&u);
        b.set_phase("dialog_gcd_raw_ipmul_seed_terminal_u");
        b.x(u[0]);

        b.set_phase("dialog_gcd_raw_ipmul_uncompute_tobitvector");
        emit_dialog_gcd_raw_tobitvector_steps_reverse(b, &u, factor, &dialog_log);

        b.set_phase("dialog_gcd_raw_ipmul_unload_p");
        for i in 0..N {
            if bit(p, i) {
                b.x(u[i]);
            }
        }
        b.free_vec(&u);
        b.free_vec(&dialog_log);
        return;
    }

    let tmp = b.alloc_qubits(N);
    b.set_phase("dialog_gcd_raw_ipmul_apply_bitvector");
    emit_dialog_gcd_raw_apply_bitvector(b, &dialog_log, target, &tmp, p);

    b.set_phase("dialog_gcd_raw_ipmul_swap_product_into_target");
    for i in 0..N {
        b.swap(target[i], tmp[i]);
    }

    b.set_phase("dialog_gcd_raw_ipmul_free_zero_tmp");
    b.free_vec(&tmp);

    b.set_phase("dialog_gcd_raw_ipmul_uncompute_tobitvector");
    emit_dialog_gcd_raw_tobitvector_steps_reverse(b, &u, factor, &dialog_log);

    b.set_phase("dialog_gcd_raw_ipmul_unload_p");
    for i in 0..N {
        if bit(p, i) {
            b.x(u[i]);
        }
    }
    b.free_vec(&u);
    b.free_vec(&dialog_log);
}

pub(crate) fn emit_dialog_gcd_raw_quotient(b: &mut B, factor: &[QubitId], target: &[QubitId], p: U256) {
    assert_eq!(factor.len(), N);
    assert_eq!(target.len(), N);

    if dialog_gcd_compressed_sidecar_log_enabled() {
        emit_dialog_gcd_compressed_sidecar_quotient(b, factor, target, p);
        return;
    }

    let dialog_log = b.alloc_qubits(DIALOG_GCD_RAW_LOG_BITS);
    let u = b.alloc_qubits(N);
    b.set_phase("dialog_gcd_raw_quotient_load_p");
    for i in 0..N {
        if bit(p, i) {
            b.x(u[i]);
        }
    }

    b.set_phase("dialog_gcd_raw_quotient_tobitvector");
    emit_dialog_gcd_raw_tobitvector_steps(b, &u, factor, &dialog_log);

    if dialog_gcd_raw_quotient_keep_terminal_u_enabled() {
        b.set_phase("dialog_gcd_raw_quotient_zero_terminal_u_for_borrow");
        b.x(u[0]);

        b.set_phase("dialog_gcd_raw_quotient_apply_reverse_reuse_factor_zero_keep_u");
        emit_dialog_gcd_raw_apply_bitvector_reverse_borrowed_subtrahend(
            b,
            &dialog_log,
            factor,
            target,
            p,
            &u,
        );

        b.set_phase("dialog_gcd_raw_quotient_swap_quotient_into_target_keep_u");
        for i in 0..N {
            b.swap(target[i], factor[i]);
        }

        b.set_phase("dialog_gcd_raw_quotient_restore_terminal_u_after_borrow");
        b.x(u[0]);

        b.set_phase("dialog_gcd_raw_quotient_uncompute_tobitvector_keep_u");
        emit_dialog_gcd_raw_tobitvector_steps_reverse(b, &u, factor, &dialog_log);

        b.set_phase("dialog_gcd_raw_quotient_unload_p_keep_u");
        for i in 0..N {
            if bit(p, i) {
                b.x(u[i]);
            }
        }
        b.free_vec(&u);
        b.free_vec(&dialog_log);
        return;
    }

    if dialog_gcd_raw_quotient_terminal_reuse_enabled() {
        b.set_phase("dialog_gcd_raw_quotient_release_terminal_u");
        b.x(u[0]);
        b.free_vec(&u);

        b.set_phase("dialog_gcd_raw_quotient_apply_reverse_reuse_factor_zero");
        emit_dialog_gcd_raw_apply_bitvector_reverse_exact(b, &dialog_log, factor, target, p);

        b.set_phase("dialog_gcd_raw_quotient_swap_quotient_into_target");
        for i in 0..N {
            b.swap(target[i], factor[i]);
        }

        b.set_phase("dialog_gcd_raw_quotient_reacquire_terminal_u");
        b.reacquire_vec(&u);
        b.set_phase("dialog_gcd_raw_quotient_seed_terminal_u");
        b.x(u[0]);

        b.set_phase("dialog_gcd_raw_quotient_uncompute_tobitvector");
        emit_dialog_gcd_raw_tobitvector_steps_reverse(b, &u, factor, &dialog_log);

        b.set_phase("dialog_gcd_raw_quotient_unload_p");
        for i in 0..N {
            if bit(p, i) {
                b.x(u[i]);
            }
        }
        b.free_vec(&u);
        b.free_vec(&dialog_log);
        return;
    }

    let tmp = b.alloc_qubits(N);
    b.set_phase("dialog_gcd_raw_quotient_apply_reverse");
    emit_dialog_gcd_raw_apply_bitvector_reverse_exact(b, &dialog_log, &tmp, target, p);

    b.set_phase("dialog_gcd_raw_quotient_swap_quotient_into_target");
    for i in 0..N {
        b.swap(target[i], tmp[i]);
    }

    b.set_phase("dialog_gcd_raw_quotient_free_zero_tmp");
    b.free_vec(&tmp);

    b.set_phase("dialog_gcd_raw_quotient_uncompute_tobitvector");
    emit_dialog_gcd_raw_tobitvector_steps_reverse(b, &u, factor, &dialog_log);

    b.set_phase("dialog_gcd_raw_quotient_unload_p");
    for i in 0..N {
        if bit(p, i) {
            b.x(u[i]);
        }
    }
    b.free_vec(&u);
    b.free_vec(&dialog_log);
}

pub(crate) fn emit_dialog_gcd_raw_pa(
    b: &mut B,
    tx: &[QubitId],
    ty: &[QubitId],
    ox: &[BitId],
    oy: &[BitId],
    p: U256,
) {
    assert_eq!(tx.len(), N);
    assert_eq!(ty.len(), N);
    assert_eq!(ox.len(), N);
    assert_eq!(oy.len(), N);

    b.set_phase("dialog_gcd_raw_pa_pair1_quotient");
    emit_dialog_gcd_raw_quotient(b, tx, ty, p);
    if dialog_gcd_raw_pa_stop_after_quotient_enabled() {
        return;
    }

    round84_emit_fused_square_xtail(b, tx, ty, ox, p);
    if dialog_gcd_raw_pa_stop_after_xtail_enabled() {
        return;
    }

    b.set_phase("dialog_gcd_raw_pa_c_ox_minus_rx");
    if dialog_fuse_c_form_enabled() {
        mod_add_triple_qb(b, tx, ox, p);
    } else {
        mod_sub_qb(b, tx, ox, p);
        mod_neg_inplace_fast(b, tx, p);
    }
    if dialog_gcd_raw_pa_stop_after_c_enabled() {
        return;
    }

    b.set_phase("dialog_gcd_raw_pa_pair2_product");
    emit_dialog_gcd_raw_ipmul(b, tx, ty, p);
    if dialog_gcd_raw_pa_stop_after_pair2_enabled() {
        return;
    }

    b.set_phase("dialog_gcd_raw_pa_y_output");
    mod_sub_qb(b, ty, oy, p);

    b.set_phase("dialog_gcd_raw_pa_x_restore");
    if dialog_fuse_x_restore_enabled() {
        mod_const_minus_reg_qb(b, tx, ox, p);
    } else {
        mod_neg_inplace_fast(b, tx, p);
        mod_add_qb(b, tx, ox, p);
    }
}

