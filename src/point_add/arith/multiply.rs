//! Multiplication and squaring: schoolbook + Karatsuba multiply, symmetric
//! squaring (incl. self-hosted / hosted variants), the controlled add/subtract
//! used by the schoolbook walk, and the `squaring_sub_from_acc_*` reducers.
use super::*;

/// Low-peak variant of `mod_mul_write_into_zero_acc_schoolbook`: uses
/// `schoolbook_mul_into_addsub_lowq` + `_inverse_lowq` instead of the fast
/// variants, saving ~n qubits at peak at the cost of ~n extra Toffolis per
/// row.
///
/// NOTE: microbench (n=256) shows this DOES NOT reduce the local peak
/// (schoolbook_fast 1797 = schoolbook_lowq 1797); the Solinas reduction +
/// acc lifetimes already dominate, and the lowq carry saving is hidden
/// underneath. We also observed a deterministic phase-garbage batch when
/// wiring this in at pair1_mul1 (1/20480 shots, ALT_SEED tag=5, across
/// two runs), so this helper is currently DEAD CODE kept only as a paper
/// trail for the negative result. See `autoresearch.ideas.md`.
#[allow(dead_code)]
pub(crate) fn mod_mul_write_into_zero_acc_schoolbook_lowq(
    b: &mut B,
    acc: &[QubitId],
    x: &[QubitId],
    y: &[QubitId],
    p: U256,
) {
    let n = acc.len();
    debug_assert_eq!(n, 256);

    let tmp_ext = b.alloc_qubits(2 * n);
    schoolbook_mul_into_addsub_lowq(b, x, y, &tmp_ext);

    let lo: Vec<QubitId> = tmp_ext[0..n].to_vec();
    let hi: Vec<QubitId> = tmp_ext[n..2 * n].to_vec();
    mod_add_qq_fast_from_zero(b, acc, &lo, p);
    mod_add_qq_fast(b, acc, &hi, p);
    for _ in 0..4 {
        mod_double_inplace_fast(b, &hi, p);
    }
    mod_add_qq_fast(b, acc, &hi, p);
    for _ in 0..2 {
        mod_double_inplace_fast(b, &hi, p);
    }
    mod_sub_qq_fast(b, acc, &hi, p);
    for _ in 0..4 {
        mod_double_inplace_fast(b, &hi, p);
    }
    mod_add_qq_fast(b, acc, &hi, p);
    let (spill, flag_inv, ovf) = mod_shift_left_by_k(b, &hi, p, 22);
    mod_add_qq(b, acc, &hi, p);
    mod_shift_right_by_k(b, &hi, p, 22, spill, flag_inv, ovf);
    for _ in 0..10 {
        mod_halve_inplace_fast(b, &hi, p);
    }

    schoolbook_mul_into_addsub_lowq_inverse(b, x, y, &tmp_ext);
    b.free_vec(&tmp_ext);
}

// ─────────────────────────────────────────────────────────────────────────────────────
// Litinski add-subtract (arXiv:2410.00899) primitives
// ─────────────────────────────────────────────────────────────────────────────────────

/// Low-peak variant of `controlled_add_subtract_fast` using non-fast
/// Cuccaro (no carry ancillae). Saves ~n qubits of transient peak at the
/// cost of ~n extra Toffolis per call. Useful when called inside the
/// Kaliski-body mul sites where peak is tight.
pub(crate) fn controlled_add_subtract_lowq(
    b: &mut B,
    x: &[QubitId],
    acc: &[QubitId],
    ctrl: QubitId,
) {
    let n = x.len();
    debug_assert_eq!(acc.len(), n + 1);

    let pad = b.alloc_qubit();
    let mut x_ext = x.to_vec();
    x_ext.push(pad);

    let c_in = b.alloc_qubit();

    b.x(ctrl);
    for i in 0..n {
        b.cx(ctrl, x_ext[i]);
    }
    b.cx(ctrl, c_in);

    cuccaro_add(b, &x_ext, acc, c_in);

    b.cx(ctrl, c_in);
    for i in 0..n {
        b.cx(ctrl, x_ext[i]);
    }
    b.x(ctrl);

    b.free(c_in);
    b.free(pad);
}

/// Inverse of `controlled_add_subtract_lowq`.
pub(crate) fn controlled_add_subtract_lowq_inverse(
    b: &mut B,
    x: &[QubitId],
    acc: &[QubitId],
    ctrl: QubitId,
) {
    let n = x.len();
    debug_assert_eq!(acc.len(), n + 1);

    let pad = b.alloc_qubit();
    let mut x_ext = x.to_vec();
    x_ext.push(pad);

    let c_in = b.alloc_qubit();

    b.x(ctrl);
    for i in 0..n {
        b.cx(ctrl, x_ext[i]);
    }
    b.cx(ctrl, c_in);

    cuccaro_sub(b, &x_ext, acc, c_in);

    b.cx(ctrl, c_in);
    for i in 0..n {
        b.cx(ctrl, x_ext[i]);
    }
    b.x(ctrl);

    b.free(c_in);
    b.free(pad);
}

/// Low-peak variant of `schoolbook_mul_into_addsub`: uses non-fast Cuccaro
/// (`cuccaro_add`) inside the `controlled_add_subtract` core and in the
/// correction adders. Saves roughly `n` transient qubits at peak vs. the
/// `_fast` variant at the cost of ~n extra Toffolis per row. Top-level
/// semantics identical to `schoolbook_mul_into_addsub`.
pub(crate) fn schoolbook_mul_into_addsub_lowq(
    b: &mut B,
    x: &[QubitId],
    y: &[QubitId],
    tmp_ext: &[QubitId],
) {
    let n = x.len();
    debug_assert_eq!(y.len(), n);
    debug_assert_eq!(tmp_ext.len(), 2 * n);

    let low = b.alloc_qubit();
    let mut wide: Vec<QubitId> = Vec::with_capacity(2 * n + 1);
    wide.push(low);
    wide.extend_from_slice(tmp_ext);

    for k in 0..n {
        let slice: Vec<QubitId> = wide[k..k + n + 1].to_vec();
        controlled_add_subtract_lowq(b, x, &slice, y[k]);
    }

    // +2^n * (y + 1)
    {
        let pad = b.alloc_qubit();
        let mut y_ext = y.to_vec();
        y_ext.push(pad);
        let slice: Vec<QubitId> = wide[n..2 * n + 1].to_vec();
        let c_in = b.alloc_qubit();
        b.x(c_in);
        cuccaro_add(b, &y_ext, &slice, c_in);
        b.x(c_in);
        b.free(c_in);
        b.free(pad);
    }

    // -2^{2n}
    b.x(wide[2 * n]);

    // -x full (2n+1)-bit sub
    {
        let mut x_ext: Vec<QubitId> = x.to_vec();
        while x_ext.len() < 2 * n + 1 {
            x_ext.push(b.alloc_qubit());
        }
        let c_in = b.alloc_qubit();
        cuccaro_sub(b, &x_ext, &wide, c_in);
        b.free(c_in);
        for _ in n..2 * n + 1 {
            let q = x_ext.pop().unwrap();
            b.free(q);
        }
    }

    // +2^n * x
    {
        let pad = b.alloc_qubit();
        let mut x_ext = x.to_vec();
        x_ext.push(pad);
        let slice: Vec<QubitId> = wide[n..2 * n + 1].to_vec();
        let c_in = b.alloc_qubit();
        cuccaro_add(b, &x_ext, &slice, c_in);
        b.free(c_in);
        b.free(pad);
    }

    b.free(low);
}

/// Exact gate-level inverse of `schoolbook_mul_into_addsub_lowq`.
pub(crate) fn schoolbook_mul_into_addsub_lowq_inverse(
    b: &mut B,
    x: &[QubitId],
    y: &[QubitId],
    tmp_ext: &[QubitId],
) {
    let n = x.len();
    debug_assert_eq!(y.len(), n);
    debug_assert_eq!(tmp_ext.len(), 2 * n);

    let low = b.alloc_qubit();
    let mut wide: Vec<QubitId> = Vec::with_capacity(2 * n + 1);
    wide.push(low);
    wide.extend_from_slice(tmp_ext);

    // Reverse correction 4: sub x at bit n.
    {
        let pad = b.alloc_qubit();
        let mut x_ext = x.to_vec();
        x_ext.push(pad);
        let slice: Vec<QubitId> = wide[n..2 * n + 1].to_vec();
        let c_in = b.alloc_qubit();
        cuccaro_sub(b, &x_ext, &slice, c_in);
        b.free(c_in);
        b.free(pad);
    }
    // Reverse correction 3.
    {
        let mut x_ext: Vec<QubitId> = x.to_vec();
        while x_ext.len() < 2 * n + 1 {
            x_ext.push(b.alloc_qubit());
        }
        let c_in = b.alloc_qubit();
        cuccaro_add(b, &x_ext, &wide, c_in);
        b.free(c_in);
        for _ in n..2 * n + 1 {
            let q = x_ext.pop().unwrap();
            b.free(q);
        }
    }
    // Reverse correction 2.
    b.x(wide[2 * n]);
    // Reverse correction 1.
    {
        let pad = b.alloc_qubit();
        let mut y_ext = y.to_vec();
        y_ext.push(pad);
        let slice: Vec<QubitId> = wide[n..2 * n + 1].to_vec();
        let c_in = b.alloc_qubit();
        b.x(c_in);
        cuccaro_sub(b, &y_ext, &slice, c_in);
        b.x(c_in);
        b.free(c_in);
        b.free(pad);
    }
    for k in (0..n).rev() {
        let slice: Vec<QubitId> = wide[k..k + n + 1].to_vec();
        controlled_add_subtract_lowq_inverse(b, x, &slice, y[k]);
    }

    b.free(low);
}

// ═══════════════════════════════════════════════════════════════════════════
//  1-level Karatsuba multiplication
// ═══════════════════════════════════════════════════════════════════════════

pub(crate) fn karatsuba_half_sum_compute(
    b: &mut B,
    lo: &[QubitId],
    hi: &[QubitId],
    acc: &[QubitId],
) {
    let h = lo.len();
    debug_assert_eq!(h, hi.len());
    debug_assert_eq!(acc.len(), h + 1);
    for i in 0..h {
        b.cx(lo[i], acc[i]);
    }
    let hi_pad = b.alloc_qubit();
    let mut hi_ext = hi.to_vec();
    hi_ext.push(hi_pad);
    add_nbit_qq_fast(b, &hi_ext, acc);
    b.free(hi_pad);
}

pub(crate) fn karatsuba_half_sum_uncompute(
    b: &mut B,
    lo: &[QubitId],
    hi: &[QubitId],
    acc: &[QubitId],
) {
    let h = lo.len();
    let hi_pad = b.alloc_qubit();
    let mut hi_ext = hi.to_vec();
    hi_ext.push(hi_pad);
    sub_nbit_qq_fast(b, &hi_ext, acc);
    b.free(hi_pad);
    for i in 0..h {
        b.cx(lo[i], acc[i]);
    }
}

// ─── 2-level Karatsuba variants (recursive on inner half-mults) ───
// Costs 2 extra z1_inner registers of ~2*(n/4+1) qubits each (~260 total for n=256).
// Higher peak qubits; use only at low-peak mul sites.

/// Symmetric schoolbook for squaring: x² = sum_i x[i]·2^(2i) + sum_{i<j} 2·x[i]·x[j]·2^(i+j).
/// Each cross-product is computed ONCE (instead of twice in full schoolbook),
/// halving the AND count + Cuccaro_add length. Saves ~130k CCX per squaring.
///
/// Row i layout (width n-i): bit 0 = diagonal x[i] at position 2i, bit 1 = 0
/// (gap), bit k+2 = cross-product (x[i] AND x[i+1+k]) at position i+(i+1+k)+1.
pub(crate) fn schoolbook_square_symmetric(b: &mut B, x: &[QubitId], tmp_ext: &[QubitId]) {
    let n = x.len();
    debug_assert_eq!(tmp_ext.len(), 2 * n);
    for i in 0..n {
        // Width: bit 0 = diag at pos 2i, bit 1 = gap, bits 2..(n-i) = cross-
        // products at positions 2i+2..i+n. Last bit index = n-i, so width = n-i+1.
        // Edge case: i = n-1 has only the diagonal, width = 1.
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let num_cross = if i + 1 < n { n - i - 1 } else { 0 };
        // num_cross = number of cross-products in this row = width - 2 when width >= 2.
        let row = b.alloc_qubits(width);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            b.ccx(x[i], x[i + 1 + k], row[k + 2]);
        }
        let pad = b.alloc_qubit();
        let mut row_padded = row.clone();
        row_padded.push(pad);
        let slice: Vec<QubitId> = tmp_ext[2 * i..2 * i + width + 1].to_vec();
        let c_in = b.alloc_qubit();
        cuccaro_add_fast(b, &row_padded, &slice, c_in);
        b.free(c_in);
        b.free(pad);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            let m = b.alloc_bit();
            b.hmr(row[k + 2], m);
            b.cz_if(x[i], x[i + 1 + k], m);
        }
        b.free_vec(&row);
    }
}

pub(crate) fn schoolbook_square_symmetric_inverse(b: &mut B, x: &[QubitId], tmp_ext: &[QubitId]) {
    let n = x.len();
    for i in (0..n).rev() {
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let num_cross = if i + 1 < n { n - i - 1 } else { 0 };
        let row = b.alloc_qubits(width);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            b.ccx(x[i], x[i + 1 + k], row[k + 2]);
        }
        let pad = b.alloc_qubit();
        let mut row_padded = row.clone();
        row_padded.push(pad);
        let slice: Vec<QubitId> = tmp_ext[2 * i..2 * i + width + 1].to_vec();
        let c_in = b.alloc_qubit();
        cuccaro_sub_fast(b, &row_padded, &slice, c_in);
        b.free(c_in);
        b.free(pad);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            let m = b.alloc_bit();
            b.hmr(row[k + 2], m);
            b.cz_if(x[i], x[i + 1 + k], m);
        }
        b.free_vec(&row);
    }
}

pub(crate) fn schoolbook_square_symmetric_lowq(b: &mut B, x: &[QubitId], tmp_ext: &[QubitId]) {
    let n = x.len();
    debug_assert_eq!(tmp_ext.len(), 2 * n);
    for i in 0..n {
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let num_cross = if i + 1 < n { n - i - 1 } else { 0 };
        let row = b.alloc_qubits(width);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            b.ccx(x[i], x[i + 1 + k], row[k + 2]);
        }
        let pad = b.alloc_qubit();
        let mut row_padded = row.clone();
        row_padded.push(pad);
        let slice: Vec<QubitId> = tmp_ext[2 * i..2 * i + width + 1].to_vec();
        let c_in = b.alloc_qubit();
        cuccaro_add(b, &row_padded, &slice, c_in);
        b.free(c_in);
        b.free(pad);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            let m = b.alloc_bit();
            b.hmr(row[k + 2], m);
            b.cz_if(x[i], x[i + 1 + k], m);
        }
        b.free_vec(&row);
    }
}

pub(crate) fn schoolbook_square_symmetric_lowq_inverse(
    b: &mut B,
    x: &[QubitId],
    tmp_ext: &[QubitId],
) {
    let n = x.len();
    for i in (0..n).rev() {
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let num_cross = if i + 1 < n { n - i - 1 } else { 0 };
        let row = b.alloc_qubits(width);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            b.ccx(x[i], x[i + 1 + k], row[k + 2]);
        }
        let pad = b.alloc_qubit();
        let mut row_padded = row.clone();
        row_padded.push(pad);
        let slice: Vec<QubitId> = tmp_ext[2 * i..2 * i + width + 1].to_vec();
        let c_in = b.alloc_qubit();
        cuccaro_sub(b, &row_padded, &slice, c_in);
        b.free(c_in);
        b.free(pad);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            let m = b.alloc_bit();
            b.hmr(row[k + 2], m);
            b.cz_if(x[i], x[i + 1 + k], m);
        }
        b.free_vec(&row);
    }
}

/// Like `schoolbook_square_symmetric` (fast, measurement UMA) but the per-row
/// Cuccaro carry lane is hosted on a caller-supplied clean register `host`
/// (returned clean) instead of a fresh allocation. Toffoli-identical to the
/// fast square, peak-identical to the lowq square — used for the z0 lobe of the
/// round84 Karatsuba square, where the not-yet-written z2 slice is clean scratch.
pub(crate) fn schoolbook_square_symmetric_hosted(
    b: &mut B,
    x: &[QubitId],
    tmp_ext: &[QubitId],
    host: &[QubitId],
) {
    let n = x.len();
    debug_assert_eq!(tmp_ext.len(), 2 * n);
    if square_selfhost_safe_lane_reuse_enabled() {
        assert_qubit_slices_disjoint(&[x, tmp_ext, host]);
    }
    for i in 0..n {
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let num_cross = if i + 1 < n { n - i - 1 } else { 0 };
        let row = b.alloc_qubits(width);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            b.ccx(x[i], x[i + 1 + k], row[k + 2]);
        }
        let slice: Vec<QubitId> = tmp_ext[2 * i..2 * i + width + 1].to_vec();
        if square_selfhost_safe_lane_reuse_enabled() {
            // The z2 sibling host is clean and disjoint from x and z0.  It has
            // ample room for both the width carry lanes and one clean c_in.
            assert!(host.len() > width);
            cuccaro_add_fast_low_to_ext_borrowed_carries(
                b,
                &row,
                &slice,
                host[width],
                &host[..width],
            );
        } else {
            let pad = b.alloc_qubit();
            let mut row_padded = row.clone();
            row_padded.push(pad);
            let c_in = b.alloc_qubit();
            cuccaro_add_fast_borrowed_carries(
                b,
                &row_padded,
                &slice,
                c_in,
                &host[..row_padded.len() - 1],
            );
            b.free(c_in);
            b.free(pad);
        }
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            let m = b.alloc_bit();
            b.hmr(row[k + 2], m);
            b.cz_if(x[i], x[i + 1 + k], m);
        }
        b.free_vec(&row);
    }
}

pub(crate) fn schoolbook_square_symmetric_hosted_inverse(
    b: &mut B,
    x: &[QubitId],
    tmp_ext: &[QubitId],
    host: &[QubitId],
) {
    let n = x.len();
    if square_selfhost_safe_lane_reuse_enabled() {
        assert_qubit_slices_disjoint(&[x, tmp_ext, host]);
    }
    for i in (0..n).rev() {
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let num_cross = if i + 1 < n { n - i - 1 } else { 0 };
        let row = b.alloc_qubits(width);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            b.ccx(x[i], x[i + 1 + k], row[k + 2]);
        }
        let slice: Vec<QubitId> = tmp_ext[2 * i..2 * i + width + 1].to_vec();
        if square_selfhost_safe_lane_reuse_enabled() {
            assert!(host.len() > width);
            cuccaro_sub_fast_low_to_ext_borrowed_carries(
                b,
                &row,
                &slice,
                host[width],
                &host[..width],
            );
        } else {
            let pad = b.alloc_qubit();
            let mut row_padded = row.clone();
            row_padded.push(pad);
            let c_in = b.alloc_qubit();
            cuccaro_sub_fast_borrowed_carries(
                b,
                &row_padded,
                &slice,
                c_in,
                &host[..row_padded.len() - 1],
            );
            b.free(c_in);
            b.free(pad);
        }
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            let m = b.alloc_bit();
            b.hmr(row[k + 2], m);
            b.cz_if(x[i], x[i + 1 + k], m);
        }
        b.free_vec(&row);
    }
}

/// Experimental square-only reclaim.  This is deliberately opt-in: every lane
/// borrowed by the prototype is either an untouched high tail of the square
/// accumulator, a caller-proved square bit that is exactly zero, or a clean
/// sibling square destination.  Dirty-but-idle data and operand aliases are not
/// eligible.
pub(crate) fn square_selfhost_safe_lane_reuse_enabled() -> bool {
    std::env::var("SQUARE_SELFHOST_SAFE_LANE_REUSE")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn assert_qubit_slices_disjoint(slices: &[&[QubitId]]) {
    let mut seen = std::collections::BTreeSet::new();
    for slice in slices {
        for &q in *slice {
            assert!(seen.insert(q), "scratch lane q{} aliases an operand", q.0);
        }
    }
}

pub(crate) fn square_selfhost_gate_suffix_carries(n: usize) -> usize {
    std::env::var("SQUARE_SELFHOST_GATE_SUFFIX_CARRIES")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0)
        .min(n.saturating_sub(1))
}

/// Like `schoolbook_square_symmetric_lowq` but converts the per-row Cuccaro
/// UMA-uncompute (CCX, executed every shot) into measurement-based (fast)
/// uncompute, WITHOUT a separate clean host register. The fast carry lane is
/// hosted on the slice's OWN not-yet-written high zeros
/// (`tmp_ext[2i+width+1 ..]`, which rows 0..=i never touch) topped up with a
/// small global remainder (<=3 qubits, since the lane width exceeds the clean
/// tail by exactly the 3-bit diagonal/gap/pad overhead). Unlike
/// `schoolbook_square_symmetric_hosted` this needs no sibling clean register,
/// so it applies where the sibling slice is occupied (the Karatsuba z2 square).
/// Peak rises only by the global remainder (<=3); Toffoli drops by the whole
/// UMA-uncompute. Under `SQUARE_SELFHOST_SAFE_LANE_REUSE=1`, the source-high
/// zero is represented structurally (no allocated `pad`) and an optional
/// caller-proved clean supplement is consumed before the global remainder. The
/// borrowed carries are returned clean by the HMR uncompute.
/// Peak-bounded row window for the selfhosted square. When set (>=2), each
/// schoolbook square row's add into `tmp_ext` is sliced into this many windows;
/// the transient row register holds only one window's worth of cross-term
/// qubits at a time (peak ~= 1024 + width/windows + boundary carries) instead of
/// the full row (peak ~= 1024 + 257). Value-exact: the same product lands in
/// `tmp_ext`. Cost: a per-boundary carry-clean comparator that rebuilds the row
/// prefix (extra CCX), traded for the dropped peak qubits.
pub(crate) fn square_row_windows() -> usize {
    std::env::var("SQUARE_ROW_WINDOWS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0)
}

/// Minimum row width below which a row is built monolithically (windowing a
/// narrow row buys no peak but still pays the comparator tax).
fn square_row_window_min_width() -> usize {
    std::env::var("SQUARE_ROW_WINDOW_MIN_WIDTH")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(96)
}

/// When >0, each row is windowed into the *minimum* number of windows that
/// keeps every window's source segment <= this width. Rows narrow enough to fit
/// in one segment are built monolithically (no comparator tax). This minimizes
/// the carry-recovery comparator overhead: only the rows wide enough to break
/// the peak budget get windowed, and only into as many windows as needed. When
/// set, it overrides the fixed SQUARE_ROW_WINDOWS count.
fn square_row_max_seg() -> usize {
    std::env::var("SQUARE_ROW_MAX_SEG")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0)
}

/// Optional truncation for the row-window boundary-carry cleanup comparator.
/// Default 0 means exact/full-width.  When set below the segment width, cleanup
/// compares only the high suffix of the segment and final partial sum.  This is
/// a deliberate island-hunt knob: it keeps the same low peak and saves Toffoli,
/// but wrong suffix ties leave the boundary carry dirty.
fn square_cleanup_direction(raw: &str) -> Option<bool> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "f" | "forward" | "false" | "0" => Some(false),
        "r" | "reverse" | "true" | "1" => Some(true),
        _ => None,
    }
}

fn square_row_window_clean_compare_bits(row: usize, window: usize, reverse: bool) -> usize {
    let default_bits = std::env::var("SQUARE_ROW_WINDOW_CLEAN_COMPARE_BITS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0);
    let row_bits = std::env::var("SQUARE_ROW_WINDOW_CLEAN_ROW_BITS")
        .ok()
        .and_then(|spec| {
            spec.split(',').rev().find_map(|item| {
                let (raw_row, raw_bits) = item.trim().split_once(':')?;
                if raw_row.trim().parse::<usize>().ok()? != row {
                    return None;
                }
                raw_bits
                    .trim()
                    .parse::<usize>()
                    .ok()
                    .filter(|bits| (1..=N).contains(bits))
            })
        })
        .unwrap_or(default_bits);
    let Ok(spec) = std::env::var("SQUARE_ROW_WINDOW_CLEAN_SITE_BITS") else {
        return row_bits;
    };
    for item in spec.split(',').rev() {
        let fields: Vec<_> = item.trim().split(':').map(str::trim).collect();
        if fields.len() != 4 {
            continue;
        }
        if fields[0].parse::<usize>().ok() != Some(row)
            || fields[1].parse::<usize>().ok() != Some(window)
            || square_cleanup_direction(fields[2]) != Some(reverse)
        {
            continue;
        }
        if let Ok(bits) = fields[3].parse::<usize>() {
            if (1..=N).contains(&bits) {
                return bits;
            }
        }
    }
    row_bits
}

fn square_row_window_measured_carry_clear_enabled() -> bool {
    std::env::var("SQUARE_ROW_WINDOW_MEASURED_CARRY_CLEAR")
        .ok()
        .as_deref()
        == Some("1")
}

/// Set row bit `j` of square row `i` into `t`. Bit 0 = x_i (diagonal low),
/// bit 1 = 0 (gap), bit 2+k = x_i & x_{i+1+k} (doubled cross term).
fn square_row_bit_set(b: &mut B, x: &[QubitId], i: usize, j: usize, t: QubitId) {
    if j == 0 {
        b.cx(x[i], t);
    } else if j == 1 {
        // gap bit: zero, nothing to do
    } else {
        b.ccx(x[i], x[i + 1 + (j - 2)], t);
    }
}

/// Measurement-based clear of a row bit set by `square_row_bit_set` (the bit is
/// known to equal its set expression at clear time).
fn square_row_bit_clear_hmr(b: &mut B, x: &[QubitId], i: usize, j: usize, t: QubitId) {
    if j == 0 {
        b.cx(x[i], t);
    } else if j == 1 {
        // gap bit: nothing
    } else {
        let m = b.alloc_bit();
        b.hmr(t, m);
        b.cz_if(x[i], x[i + 1 + (j - 2)], m);
    }
}

/// Windowed selfhosted square row add: `tmp_ext[2i ..] += row_i` where
/// `row_i` has `width` bits, built one window at a time. `forward=true` adds,
/// `forward=false` subtracts (the inverse). Value-identical to a single
/// `cuccaro_{add,sub}` of the full row into `tmp_ext[2i..2i+width+1]`.
///
/// The full-width add is split into a chain of low-to-ext adds. Window `w`
/// covers row bits `[lo..hi)` and writes `tmp_ext[base+lo .. base+hi+1]` (the
/// extra high cell absorbs the window carry). Because windows are contiguous in
/// `tmp_ext`, window `w`'s carry lands in `tmp_ext[base+hi]`, which is the low
/// cell of window `w+1` — so the carry chains *through* `tmp_ext` with no
/// separate carry-out ancilla and no boundary comparators. The per-window
/// Cuccaro carry lane is borrowed from `tmp_ext`'s not-yet-written high zeros
/// (rows `0..=i` never touch `tmp_ext[2i+width+1 ..]`), topped up by a small
/// global remainder, so the transient overhead is only the `seg_w`-wide source
/// window. Forward order low→high; inverse must mirror it high→low.
fn square_row_windowed_apply(
    b: &mut B,
    x: &[QubitId],
    tmp_ext: &[QubitId],
    i: usize,
    width: usize,
    windows: usize,
    forward: bool,
) {
    let base = 2 * i;
    let windows = windows.max(1).min(width);

    // Window boundaries over the row bit range [0, width).
    let bounds: Vec<(usize, usize)> = (0..windows)
        .map(|w| {
            let lo = (w * width) / windows;
            let hi = ((w + 1) * width) / windows;
            (lo, hi)
        })
        .filter(|&(lo, hi)| hi > lo)
        .collect();
    let nwin = bounds.len();

    // Each interior window's carry-out (forward) / borrow-out (inverse) is
    // captured in a fresh clean ancilla `cout`, fed as the carry/borrow-IN of
    // the next window so the carry ripples across the boundary into the
    // accumulator. The final window captures its carry into tmp_ext[base+width].
    // Interior couts are NOT clean after being consumed as the next c_in
    // (Cuccaro restores c_in to the carry value), so they are uncomputed by a
    // *local* width-bounded comparator that recovers the carry from the final
    // partial sum and the rebuilt source window — peak stays ~1024 + 2*seg_w.
    //
    // The inverse (sub) is built as the structural mirror of the forward (add):
    // same window order and carry chaining, add->sub, with the borrow-recovery
    // comparator X-wrapped per the Cuccaro sub convention. It SUBTRACTS the same
    // row value the forward ADDED, so tmp_ext returns to its pre-square state.

    let build_seg = |b: &mut B, lo: usize, hi: usize| -> Vec<QubitId> {
        let seg = b.alloc_qubits(hi - lo);
        for (k, &q) in seg.iter().enumerate() {
            square_row_bit_set(b, x, i, lo + k, q);
        }
        seg
    };
    let clear_seg = |b: &mut B, lo: usize, seg: &[QubitId]| {
        for (k, &q) in seg.iter().enumerate() {
            square_row_bit_clear_hmr(b, x, i, lo + k, q);
        }
        b.free_vec(seg);
    };

    // The Cuccaro carry lane for each window add/sub is borrowed from tmp_ext's
    // clean high tail (positions beyond this row's footprint base+width+1, which
    // no row 0..=i touches), so the per-window transient overhead is only the
    // seg_w source bits + a 0-pad + the cout ancilla (~seg_w+2), never an
    // allocated carry array. The interior carry-out cleanup uses the *slow*
    // (carry-array-free) comparator, so cleanup is peak-flat (+0 beyond seg).
    let row_top = base + width + 1; // first clean tmp_ext cell above the row.
    let borrow_lane = |b: &mut B, _need: usize| -> Vec<QubitId> {
        // Always available: tmp_ext beyond row_top is clean and >= seg_w wide
        // for every window (seg_w <= width and the high tail is wide enough).
        tmp_ext[row_top..row_top + _need].to_vec()
    };

    // carry/borrow-in for window 0 is a clean zero.
    let mut carry_in = b.alloc_qubit();
    let first_carry = carry_in;
    let mut couts: Vec<(QubitId, usize, usize, QubitId, usize)> = Vec::new();
    for (wi, &(lo, hi)) in bounds.iter().enumerate() {
        let last = wi == nwin - 1;
        let seg = build_seg(b, lo, hi);
        let seg_w = hi - lo;
        // Build a_block = seg ++ 0pad, acc_block = tmp[lo..hi] ++ high, n = seg_w+1.
        let pad = b.alloc_qubit();
        let mut a_block = seg.clone();
        a_block.push(pad);
        let high = if last {
            // Final window: high carry lands in the (clean) tmp_ext[base+width].
            tmp_ext[base + hi]
        } else {
            b.alloc_qubit()
        };
        let mut acc_block: Vec<QubitId> = tmp_ext[base + lo..base + hi].to_vec();
        acc_block.push(high);
        let nblk = a_block.len();
        let carries = borrow_lane(b, nblk - 1);
        if forward {
            cuccaro_add_fast_borrowed_carries(b, &a_block, &acc_block, carry_in, &carries);
        } else {
            cuccaro_sub_fast_borrowed_carries(b, &a_block, &acc_block, carry_in, &carries);
        }
        b.free(pad);
        if last {
            // nothing extra: carry already in tmp_ext[base+width].
        } else {
            couts.push((high, lo, hi, carry_in, wi));
            carry_in = high;
        }
        clear_seg(b, lo, &seg);
    }
    // Reverse sweep: clean each interior cout with a local comparator. The
    // measured-uncompute fast comparator (~n CCX) borrows its n-wide carry lane
    // from tmp_ext's clean high tail, so cleanup adds no peak qubits. Setting
    // SQUARE_ROW_WINDOW_SLOW_CMP=1 falls back to the carry-array-free slow
    // comparator (~2n CCX, also peak-flat) for cross-checking.
    let slow_cmp = std::env::var("SQUARE_ROW_WINDOW_SLOW_CMP").ok().as_deref() == Some("1");
    let measured_clear = square_row_window_measured_carry_clear_enabled();
    for &(cout, lo, hi, cin, window) in couts.iter().rev() {
        let clean_cmp_bits = square_row_window_clean_compare_bits(i, window, !forward);
        let seg_w = hi - lo;
        let trunc_w = if clean_cmp_bits == 0 {
            seg_w
        } else {
            clean_cmp_bits.min(seg_w)
        };
        if trunc_w < seg_w {
            let suffix_lo = hi - trunc_w;
            let seg = build_seg(b, suffix_lo, hi);
            let carries = tmp_ext[row_top..row_top + trunc_w].to_vec();
            let cmp_cin = b.alloc_qubit();
            if forward {
                if measured_clear {
                    let phase = b.alloc_bit();
                    b.hmr(cout, phase);
                    cmp_lt_phase_conditioned_with_cin_borrowed_carries(
                        b,
                        &tmp_ext[base + suffix_lo..base + hi],
                        &seg,
                        cmp_cin,
                        &carries,
                        phase,
                    );
                } else {
                    cmp_lt_into_fast_with_cin_borrowed_carries(
                        b,
                        &tmp_ext[base + suffix_lo..base + hi],
                        &seg,
                        cmp_cin,
                        cout,
                        &carries,
                    );
                }
            } else {
                for &q in &seg {
                    b.x(q);
                }
                if measured_clear {
                    let phase = b.alloc_bit();
                    b.hmr(cout, phase);
                    cmp_lt_phase_conditioned_with_cin_borrowed_carries(
                        b,
                        &seg,
                        &tmp_ext[base + suffix_lo..base + hi],
                        cmp_cin,
                        &carries,
                        phase,
                    );
                } else {
                    cmp_lt_into_fast_with_cin_borrowed_carries(
                        b,
                        &seg,
                        &tmp_ext[base + suffix_lo..base + hi],
                        cmp_cin,
                        cout,
                        &carries,
                    );
                }
                for &q in &seg {
                    b.x(q);
                }
            }
            b.free(cmp_cin);
            clear_seg(b, suffix_lo, &seg);
        } else {
            let seg = build_seg(b, lo, hi);
            let carries = tmp_ext[row_top..row_top + seg_w].to_vec();
            if forward {
                // carry_out = (partial_sum < seg + cin)
                if measured_clear {
                    let phase = b.alloc_bit();
                    b.hmr(cout, phase);
                    cmp_lt_phase_conditioned_with_cin_borrowed_carries(
                        b,
                        &tmp_ext[base + lo..base + hi],
                        &seg,
                        cin,
                        &carries,
                        phase,
                    );
                } else if slow_cmp {
                    cmp_lt_into_with_cin_slow(b, &tmp_ext[base + lo..base + hi], &seg, cin, cout);
                } else {
                    cmp_lt_into_fast_with_cin_borrowed_carries(
                        b,
                        &tmp_ext[base + lo..base + hi],
                        &seg,
                        cin,
                        cout,
                        &carries,
                    );
                }
            } else {
                // borrow_out = (seg + cin > partial_diff)
                for k in 0..seg_w {
                    b.x(seg[k]);
                }
                if measured_clear {
                    let phase = b.alloc_bit();
                    b.hmr(cout, phase);
                    cmp_lt_phase_conditioned_with_cin_borrowed_carries(
                        b,
                        &seg,
                        &tmp_ext[base + lo..base + hi],
                        cin,
                        &carries,
                        phase,
                    );
                } else if slow_cmp {
                    cmp_lt_into_with_cin_slow(b, &seg, &tmp_ext[base + lo..base + hi], cin, cout);
                } else {
                    cmp_lt_into_fast_with_cin_borrowed_carries(
                        b,
                        &seg,
                        &tmp_ext[base + lo..base + hi],
                        cin,
                        cout,
                        &carries,
                    );
                }
                for k in 0..seg_w {
                    b.x(seg[k]);
                }
            }
            clear_seg(b, lo, &seg);
        }
        b.free(cout);
    }
    b.free(first_carry);
}

pub(crate) fn schoolbook_square_symmetric_lowq_selfhosted(
    b: &mut B,
    x: &[QubitId],
    tmp_ext: &[QubitId],
) {
    schoolbook_square_symmetric_lowq_selfhosted_with_clean_supplement(b, x, tmp_ext, &[]);
}

pub(crate) fn schoolbook_square_symmetric_lowq_selfhosted_with_clean_supplement(
    b: &mut B,
    x: &[QubitId],
    tmp_ext: &[QubitId],
    clean_supplement: &[QubitId],
) {
    let n = x.len();
    debug_assert_eq!(tmp_ext.len(), 2 * n);
    let safe_reuse = square_selfhost_safe_lane_reuse_enabled();
    if safe_reuse {
        assert_qubit_slices_disjoint(&[x, tmp_ext, clean_supplement]);
    }
    let gate_prefix_rows = std::env::var("SQUARE_SELFHOST_GATE_PREFIX_ROWS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0);
    let row_windows = square_row_windows();
    let row_window_min = square_row_window_min_width();
    let max_seg = square_row_max_seg();
    for i in 0..n {
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let num_cross = if i + 1 < n { n - i - 1 } else { 0 };
        if max_seg > 0 && i >= gate_prefix_rows && width > max_seg {
            let w = width.div_ceil(max_seg);
            square_row_windowed_apply(b, x, tmp_ext, i, width, w, true);
            continue;
        }
        if max_seg == 0 && row_windows >= 1 && i >= gate_prefix_rows && width >= row_window_min {
            square_row_windowed_apply(b, x, tmp_ext, i, width, row_windows, true);
            continue;
        }
        let row = b.alloc_qubits(width);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            b.ccx(x[i], x[i + 1 + k], row[k + 2]);
        }
        let hi = 2 * i + width + 1;
        let slice: Vec<QubitId> = tmp_ext[2 * i..hi].to_vec();
        if i < gate_prefix_rows {
            let pad = b.alloc_qubit();
            let mut row_padded = row.clone();
            row_padded.push(pad);
            let c_in = b.alloc_qubit();
            cuccaro_add(b, &row_padded, &slice, c_in);
            b.free(c_in);
            b.free(pad);
        } else if safe_reuse {
            let need = row.len() - square_selfhost_gate_suffix_carries(row.len());
            let avail = tmp_ext.len() - hi;
            let from_tmp = need.min(avail);
            let from_supplement = (need - from_tmp).min(clean_supplement.len());
            let from_global = need - from_tmp - from_supplement;
            let gpool = b.alloc_qubits(from_global);
            let mut carries: Vec<QubitId> = tmp_ext[hi..hi + from_tmp].to_vec();
            carries.extend_from_slice(&clean_supplement[..from_supplement]);
            carries.extend_from_slice(&gpool);
            cuccaro_add_fast_low_to_ext_borrowed_carries_no_cin(b, &row, &slice, &carries);
            b.free_vec(&gpool);
        } else {
            let pad = b.alloc_qubit();
            let mut row_padded = row.clone();
            row_padded.push(pad);
            let c_in = b.alloc_qubit();
            let need = row_padded.len() - 1;
            let avail = tmp_ext.len() - hi;
            let from_tmp = need.min(avail);
            let from_global = need - from_tmp;
            let gpool = b.alloc_qubits(from_global);
            let mut carries: Vec<QubitId> = tmp_ext[hi..hi + from_tmp].to_vec();
            carries.extend_from_slice(&gpool);
            cuccaro_add_fast_borrowed_carries(b, &row_padded, &slice, c_in, &carries);
            b.free(c_in);
            b.free_vec(&gpool);
            b.free(pad);
        }
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            let m = b.alloc_bit();
            b.hmr(row[k + 2], m);
            b.cz_if(x[i], x[i + 1 + k], m);
        }
        b.free_vec(&row);
    }
}

pub(crate) fn schoolbook_square_symmetric_lowq_selfhosted_inverse(
    b: &mut B,
    x: &[QubitId],
    tmp_ext: &[QubitId],
) {
    schoolbook_square_symmetric_lowq_selfhosted_inverse_with_clean_supplement(b, x, tmp_ext, &[]);
}

pub(crate) fn schoolbook_square_symmetric_lowq_selfhosted_inverse_with_clean_supplement(
    b: &mut B,
    x: &[QubitId],
    tmp_ext: &[QubitId],
    clean_supplement: &[QubitId],
) {
    let n = x.len();
    debug_assert_eq!(tmp_ext.len(), 2 * n);
    let safe_reuse = square_selfhost_safe_lane_reuse_enabled();
    if safe_reuse {
        assert_qubit_slices_disjoint(&[x, tmp_ext, clean_supplement]);
    }
    let gate_prefix_rows = std::env::var("SQUARE_SELFHOST_GATE_PREFIX_ROWS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0);
    let row_windows = square_row_windows();
    let row_window_min = square_row_window_min_width();
    let max_seg = square_row_max_seg();
    for i in (0..n).rev() {
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let num_cross = if i + 1 < n { n - i - 1 } else { 0 };
        if max_seg > 0 && i >= gate_prefix_rows && width > max_seg {
            let w = width.div_ceil(max_seg);
            square_row_windowed_apply(b, x, tmp_ext, i, width, w, false);
            continue;
        }
        if max_seg == 0 && row_windows >= 1 && i >= gate_prefix_rows && width >= row_window_min {
            square_row_windowed_apply(b, x, tmp_ext, i, width, row_windows, false);
            continue;
        }
        let row = b.alloc_qubits(width);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            b.ccx(x[i], x[i + 1 + k], row[k + 2]);
        }
        let hi = 2 * i + width + 1;
        let slice: Vec<QubitId> = tmp_ext[2 * i..hi].to_vec();
        if i < gate_prefix_rows {
            let pad = b.alloc_qubit();
            let mut row_padded = row.clone();
            row_padded.push(pad);
            let c_in = b.alloc_qubit();
            cuccaro_sub(b, &row_padded, &slice, c_in);
            b.free(c_in);
            b.free(pad);
        } else if safe_reuse {
            let need = row.len() - square_selfhost_gate_suffix_carries(row.len());
            let avail = tmp_ext.len() - hi;
            let from_tmp = need.min(avail);
            let from_supplement = (need - from_tmp).min(clean_supplement.len());
            let from_global = need - from_tmp - from_supplement;
            let gpool = b.alloc_qubits(from_global);
            let mut carries: Vec<QubitId> = tmp_ext[hi..hi + from_tmp].to_vec();
            carries.extend_from_slice(&clean_supplement[..from_supplement]);
            carries.extend_from_slice(&gpool);
            cuccaro_sub_fast_low_to_ext_borrowed_carries_no_cin(b, &row, &slice, &carries);
            b.free_vec(&gpool);
        } else {
            let pad = b.alloc_qubit();
            let mut row_padded = row.clone();
            row_padded.push(pad);
            let c_in = b.alloc_qubit();
            let need = row_padded.len() - 1;
            let avail = tmp_ext.len() - hi;
            let from_tmp = need.min(avail);
            let from_global = need - from_tmp;
            let gpool = b.alloc_qubits(from_global);
            let mut carries: Vec<QubitId> = tmp_ext[hi..hi + from_tmp].to_vec();
            carries.extend_from_slice(&gpool);
            cuccaro_sub_fast_borrowed_carries(b, &row_padded, &slice, c_in, &carries);
            b.free(c_in);
            b.free_vec(&gpool);
            b.free(pad);
        }
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            let m = b.alloc_bit();
            b.hmr(row[k + 2], m);
            b.cz_if(x[i], x[i + 1 + k], m);
        }
        b.free_vec(&row);
    }
}

/// Gate for the measured-uncompute (self-hosted) Karatsuba z2 square. Defaults
/// ON; set KARA_Z2_SELFHOST=0 to fall back to the plain ancilla-free lowq z2
/// square (CCX UMA-uncompute).
pub(crate) fn kara_z2_selfhost_enabled() -> bool {
    std::env::var("KARA_Z2_SELFHOST").ok().as_deref() != Some("0")
}

/// Gate for the measured-uncompute (self-hosted) round84 x-tail full-width
/// lam^2 square. Defaults ON; set XTAIL_SQ_SELFHOST=0 to fall back to the plain
/// ancilla-free lowq square (CCX UMA-uncompute).
pub(crate) fn xtail_sq_selfhost_enabled() -> bool {
    std::env::var("XTAIL_SQ_SELFHOST").ok().as_deref() != Some("0")
}

fn round84_inplace_solinas_fold_enabled() -> bool {
    std::env::var("ROUND84_INPLACE_SOLINAS_FOLD")
        .ok()
        .as_deref()
        == Some("1")
}

// tofprof CAT-4 lever: the in-place Solinas fold/unfold build their adders from
// the COHERENT cuccaro_add/sub (maj/uma, ~2 CCX/bit, 0 carry ancilla). The
// fold/unfold phases run at active=1160, i.e. 137 qubits below the 1297 peak, so
// the SMALL adders (quotient*c product = 33-bit, narrow correction = 66-bit,
// quotient-update spill <=34-bit) can use the MEASURED cuccaro_*_fast (~1 CCX/bit
// + Hmr-uncompute) peak-neutrally => ~1 CCX/bit saved on those. The BIG fold-step
// adders (224..256-bit) are left coherent (a fast version would need ~256 carry
// lanes -> 1160+256=1416 > 1297 = peak-positive). Default OFF (byte-identical).
fn round84_fold_fast_add_enabled() -> bool {
    std::env::var("ROUND84_FOLD_FAST_ADD").ok().as_deref() == Some("1")
}
#[inline]
fn round84_add_small(b: &mut B, a: &[QubitId], acc: &[QubitId]) {
    if round84_fold_fast_add_enabled() {
        add_nbit_qq_fast(b, a, acc);
    } else {
        add_nbit_qq(b, a, acc);
    }
}
#[inline]
fn round84_sub_small(b: &mut B, a: &[QubitId], acc: &[QubitId]) {
    if round84_fold_fast_add_enabled() {
        sub_nbit_qq_fast(b, a, acc);
    } else {
        sub_nbit_qq(b, a, acc);
    }
}

fn round84_inplace_quotient_carry_trunc_window() -> usize {
    std::env::var("ROUND84_INPLACE_QUOTIENT_CARRY_TRUNC_W")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(21)
        .max(1)
}

fn round84_inplace_vent_carry_enabled() -> bool {
    std::env::var("ROUND84_INPLACE_VENT_CARRY").ok().as_deref() == Some("1")
}

fn round84_correction_wrap_borrow_quotient_top_enabled() -> bool {
    std::env::var("ROUND84_CORRECTION_WRAP_BORROW_QUOTIENT_TOP")
        .ok()
        .as_deref()
        == Some("1")
}

fn round84_keep_quotient_product_enabled() -> bool {
    std::env::var("ROUND84_KEEP_QUOTIENT_PRODUCT")
        .ok()
        .as_deref()
        == Some("1")
}

fn round84_qprod_naf_enabled() -> bool {
    std::env::var("R84_QPROD_NAF").ok().as_deref() == Some("1")
}

fn round84_qprod_short_enabled() -> bool {
    std::env::var("ROUND84_QPROD_SHORT").ok().as_deref() == Some("1")
}

struct Round84FoldStep {
    shift: usize,
    add: bool,
    wrap: QubitId,
}

struct Round84AggregateFold {
    steps: Vec<Round84FoldStep>,
    quotient: Vec<QubitId>,
    correction_wrap: QubitId,
    correction_wrap_owned: bool,
    product: Option<Vec<QubitId>>,
}

fn round84_update_fold_quotient(
    b: &mut B,
    quotient: &[QubitId],
    hi: &[QubitId],
    step: &Round84FoldStep,
    inverse: bool,
) {
    let add = step.add != inverse;
    let update_wrap = |b: &mut B| {
        if add {
            cadd_nbit_const_direct_fast(b, quotient, U256::from(1), step.wrap);
        } else {
            csub_nbit_const_direct_fast(b, quotient, U256::from(1), step.wrap);
        }
    };
    let update_spill = |b: &mut B| {
        if step.shift == 0 {
            return;
        }
        let pad = b.alloc_qubits(quotient.len() - step.shift);
        let mut spill = hi[hi.len() - step.shift..].to_vec();
        spill.extend_from_slice(&pad);
        if add {
            add_nbit_qq(b, &spill, quotient);
        } else {
            sub_nbit_qq(b, &spill, quotient);
        }
        b.free_vec(&pad);
    };

    if inverse {
        update_wrap(b);
        update_spill(b);
    } else {
        update_spill(b);
        update_wrap(b);
    }
}

fn round84_compute_quotient_c_product(
    b: &mut B,
    quotient: &[QubitId],
    dirty: &[QubitId],
) -> Vec<QubitId> {
    // quotient <= c, so its low 33 bits suffice and quotient*c fits in 66 bits.
    let q = &quotient[..33];
    let product = b.alloc_qubits(66);
    for i in 0..q.len() {
        b.cx(q[i], product[i]);
    }
    if round84_qprod_naf_enabled() {
        for (shift, add) in [(10usize, true), (32, true), (5, false), (4, false)] {
            if round84_qprod_vent_pad_enabled()
                && (product.len() - shift - q.len()) >= round84_qprod_vent_pad_min_width()
            {
                round84_qprod_shifted_addsub_vented(b, q, &product, shift, add, dirty);
                continue;
            }
            let target = &product[shift..];
            if round84_qprod_short_enabled() {
                if add {
                    add_short_to_long_qq_fast_no_cin(b, q, target);
                } else {
                    sub_short_to_long_qq_fast_no_cin(b, q, target);
                }
            } else {
                let pad = b.alloc_qubits(target.len() - q.len());
                let mut source = q.to_vec();
                source.extend_from_slice(&pad);
                if add {
                    round84_add_small(b, &source, target);
                } else {
                    round84_sub_small(b, &source, target);
                }
                b.free_vec(&pad);
            }
        }
    } else {
        for shift in [4usize, 6, 7, 8, 9, 32] {
            let target = &product[shift..];
            if round84_qprod_short_enabled() {
                add_short_to_long_qq_fast_no_cin(b, q, target);
            } else {
                let pad = b.alloc_qubits(target.len() - q.len());
                let mut source = q.to_vec();
                source.extend_from_slice(&pad);
                round84_add_small(b, &source, target);
                b.free_vec(&pad);
            }
        }
    }
    product
}

fn round84_uncompute_quotient_c_product(
    b: &mut B,
    quotient: &[QubitId],
    product: &[QubitId],
    dirty: &[QubitId],
) {
    let q = &quotient[..33];
    if round84_qprod_naf_enabled() {
        for (shift, add) in [(10usize, true), (32, true), (5, false), (4, false)]
            .into_iter()
            .rev()
        {
            if round84_qprod_vent_pad_enabled()
                && (product.len() - shift - q.len()) >= round84_qprod_vent_pad_min_width()
            {
                // uncompute = inverse op: add->sub, sub->add.
                round84_qprod_shifted_addsub_vented(b, q, product, shift, !add, dirty);
                continue;
            }
            let target = &product[shift..];
            if round84_qprod_short_enabled() {
                if add {
                    sub_short_to_long_qq_fast_no_cin(b, q, target);
                } else {
                    add_short_to_long_qq_fast_no_cin(b, q, target);
                }
            } else {
                let pad = b.alloc_qubits(target.len() - q.len());
                let mut source = q.to_vec();
                source.extend_from_slice(&pad);
                if add {
                    round84_sub_small(b, &source, target);
                } else {
                    round84_add_small(b, &source, target);
                }
                b.free_vec(&pad);
            }
        }
    } else {
        for shift in [4usize, 6, 7, 8, 9, 32].into_iter().rev() {
            let target = &product[shift..];
            if round84_qprod_short_enabled() {
                sub_short_to_long_qq_fast_no_cin(b, q, target);
            } else {
                let pad = b.alloc_qubits(target.len() - q.len());
                let mut source = q.to_vec();
                source.extend_from_slice(&pad);
                round84_sub_small(b, &source, target);
                b.free_vec(&pad);
            }
        }
    }
    for i in 0..q.len() {
        b.cx(q[i], product[i]);
    }
    b.free_vec(product);
}

fn round84_add_narrow_correction(
    b: &mut B,
    lo: &[QubitId],
    product: &[QubitId],
    dirty: &[QubitId],
    borrowed_wrap: Option<QubitId>,
) -> (QubitId, bool) {
    let (wrap, owned_wrap) = borrowed_wrap.map_or_else(|| (b.alloc_qubit(), true), |q| (q, false));
    let source_top = b.alloc_qubit();
    let mut target_ext = lo[..product.len()].to_vec();
    target_ext.push(wrap);
    let mut source_ext = product.to_vec();
    source_ext.push(source_top);
    round84_add_small(b, &source_ext, &target_ext);
    b.free(source_top);
    let high = &lo[product.len()..];
    if round84_inplace_vent_carry_enabled() {
        let clean2 = [b.alloc_qubit(), b.alloc_qubit()];
        venting::ciadd_dirty_2clean_classical(
            b,
            high,
            &dirty[..high.len() - 2],
            &clean2,
            1,
            wrap,
            false,
        );
        b.free(clean2[0]);
        b.free(clean2[1]);
    } else {
        cadd_nbit_const_direct_trunc_fast(
            b,
            high,
            U256::from(1),
            wrap,
            round84_inplace_quotient_carry_trunc_window(),
        );
    }
    (wrap, owned_wrap)
}

fn round84_sub_narrow_correction(
    b: &mut B,
    lo: &[QubitId],
    product: &[QubitId],
    wrap: QubitId,
    dirty: &[QubitId],
    owned_wrap: bool,
) {
    let high = &lo[product.len()..];
    if round84_inplace_vent_carry_enabled() {
        let clean2 = [b.alloc_qubit(), b.alloc_qubit()];
        venting::cisub_dirty_2clean_classical(b, high, &dirty[..high.len() - 2], &clean2, 1, wrap);
        b.free(clean2[0]);
        b.free(clean2[1]);
    } else {
        csub_nbit_const_direct_trunc_fast(
            b,
            high,
            U256::from(1),
            wrap,
            round84_inplace_quotient_carry_trunc_window(),
        );
    }
    let source_top = b.alloc_qubit();
    let mut target_ext = lo[..product.len()].to_vec();
    target_ext.push(wrap);
    let mut source_ext = product.to_vec();
    source_ext.push(source_top);
    round84_sub_small(b, &source_ext, &target_ext);
    b.free(source_top);
    if owned_wrap {
        b.free(wrap);
    }
}

/// Reversibly fold `hi*c` into `lo`, where `c = 2^256-p`.
///
/// Each signed shifted add/sub retains its 2^256 quotient contribution. The
/// five contributions are accumulated into a 34-bit register, multiplied by
/// sparse `c` once, and added to `lo`. The returned state is sufficient to
/// restore the original square after `lo` has been consumed.
fn round84_fold_hi_into_lo_aggregate(
    b: &mut B,
    lo: &[QubitId],
    hi: &[QubitId],
    dirty: &[QubitId],
) -> Round84AggregateFold {
    let n = lo.len();
    let quotient = b.alloc_qubits(34);
    // c = 2^32 + 977 = 2^32 + 2^10 - 2^5 - 2^4 + 1.
    let terms = [
        (0usize, true),
        (4, false),
        (5, false),
        (10, true),
        (32, true),
    ];
    let mut steps = Vec::with_capacity(terms.len());

    for (shift, add) in terms {
        let width = n - shift;
        let wrap = b.alloc_qubit();
        let source_top = b.alloc_qubit();
        let mut target_ext = lo[shift..].to_vec();
        target_ext.push(wrap);
        let mut source_ext = hi[..width].to_vec();
        source_ext.push(source_top);
        if add {
            add_nbit_qq(b, &source_ext, &target_ext);
        } else {
            sub_nbit_qq(b, &source_ext, &target_ext);
        }
        b.free(source_top);

        let step = Round84FoldStep { shift, add, wrap };
        round84_update_fold_quotient(b, &quotient, hi, &step, false);
        steps.push(step);
    }

    let product = round84_compute_quotient_c_product(b, &quotient, dirty);
    let borrowed_correction_wrap =
        round84_correction_wrap_borrow_quotient_top_enabled().then_some(quotient[33]);
    let (correction_wrap, correction_wrap_owned) =
        round84_add_narrow_correction(b, lo, &product, dirty, borrowed_correction_wrap);
    let product = if round84_keep_quotient_product_enabled() {
        Some(product)
    } else {
        round84_uncompute_quotient_c_product(b, &quotient, &product, dirty);
        None
    };
    Round84AggregateFold {
        steps,
        quotient,
        correction_wrap,
        correction_wrap_owned,
        product,
    }
}

fn round84_unfold_hi_from_lo_aggregate(
    b: &mut B,
    lo: &[QubitId],
    hi: &[QubitId],
    dirty: &[QubitId],
    state: Round84AggregateFold,
) {
    let product = state
        .product
        .unwrap_or_else(|| round84_compute_quotient_c_product(b, &state.quotient, dirty));
    round84_sub_narrow_correction(
        b,
        lo,
        &product,
        state.correction_wrap,
        dirty,
        state.correction_wrap_owned,
    );
    round84_uncompute_quotient_c_product(b, &state.quotient, &product, dirty);

    for step in state.steps.into_iter().rev() {
        round84_update_fold_quotient(b, &state.quotient, hi, &step, true);
        let width = lo.len() - step.shift;
        let source_top = b.alloc_qubit();
        let mut target_ext = lo[step.shift..].to_vec();
        target_ext.push(step.wrap);
        let mut source_ext = hi[..width].to_vec();
        source_ext.push(source_top);
        if step.add {
            sub_nbit_qq(b, &source_ext, &target_ext);
        } else {
            add_nbit_qq(b, &source_ext, &target_ext);
        }
        b.free(source_top);
        b.free(step.wrap);
    }
    b.free_vec(&state.quotient);
}

/// Schoolbook squarer with Bennett uncompute. For squaring `tmp_ext = x*x`
/// (2n bits, no mod reduction), then sub from acc with on-the-fly Solinas
/// reduction, then uncompute tmp_ext via gate-level inverse. Saves ~170k
/// CCX vs walk-x squaring (459k → 289k) by avoiding 256 expensive
/// cmod_add_qq calls (each 5n) in favor of 2n²=131k of cheap AND+Cuccaro.
pub(crate) fn squaring_sub_from_acc_schoolbook(b: &mut B, acc: &[QubitId], x: &[QubitId], p: U256) {
    let n = acc.len();
    debug_assert_eq!(n, 256);
    debug_assert_eq!(x.len(), n);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));

    // Wide accumulator (2n bits) starts at 0.
    let tmp_ext = b.alloc_qubits(2 * n);

    // Phase 1: symmetric schoolbook tmp_ext = x*x (~half the CCX of full).
    schoolbook_square_symmetric(b, x, &tmp_ext);

    // Phase 2: subtract (lo + hi*c mod p) from acc.
    // For each set bit k of c, sub (hi shifted by k mod p) from acc, by
    // walking hi via mod_double in place. Sub lo first.
    let lo: Vec<QubitId> = tmp_ext[0..n].to_vec();
    let hi: Vec<QubitId> = tmp_ext[n..2 * n].to_vec();
    mod_sub_qq_fast(b, acc, &lo, p);
    let _ = c;
    // 977 consolidation: c = {+2^0, +2^4, -2^6, +2^10, +2^32}. For acc-=hi·c, signs flip:
    // acc -= hi·2^0, acc -= hi·2^4, acc += hi·2^6, acc -= hi·2^10, acc -= hi·2^32.
    mod_sub_qq_fast(b, acc, &hi, p);
    for _ in 0..4 {
        mod_double_inplace_fast(b, &hi, p);
    }
    mod_sub_qq_fast(b, acc, &hi, p);
    for _ in 0..2 {
        mod_double_inplace_fast(b, &hi, p);
    }
    mod_add_qq_fast(b, acc, &hi, p); // sign flipped
    for _ in 0..4 {
        mod_double_inplace_fast(b, &hi, p);
    }
    mod_sub_qq_fast(b, acc, &hi, p);
    let (spill, flag_inv, ovf) = mod_shift_left_by_k(b, &hi, p, 22);
    mod_sub_qq(b, acc, &hi, p);
    mod_shift_right_by_k(b, &hi, p, 22, spill, flag_inv, ovf);
    for _ in 0..10 {
        mod_halve_inplace_fast(b, &hi, p);
    }

    // Phase 3: uncompute tmp_ext via symmetric schoolbook inverse.
    schoolbook_square_symmetric_inverse(b, x, &tmp_ext);

    b.free_vec(&tmp_ext);
}

/// Squaring-aware 1-level Karatsuba variant of [`squaring_sub_from_acc_schoolbook`].
///
/// Computes `acc -= x^2 mod p` (Solinas-reduced) via a 1-level Karatsuba
/// SQUARE. Split `x = hi‖lo` (`h = n/2` bits each) and form the three
/// SYMMETRIC sub-squares
///   z0 = lo^2,  z2 = hi^2,  z1 = (lo+hi)^2,
/// then combine `z1 -= z0 + z2` (= 2·lo·hi) and add the middle term:
///   x^2 = z0 + (z1 - z0 - z2)·2^h + z2·2^{2h}.
/// Each sub-square is the existing symmetric square (`schoolbook_square_symmetric`,
/// cross-products counted once via Gidney-uncomputed AND lanes), so the dominant
/// cross-product AND budget drops ~25 % vs the symmetric 256-bit schoolbook
/// square: 3·(n/2)(n/2-1)/2 cross ANDs instead of n(n-1)/2. Using a plain
/// Karatsuba MUL with x=y would re-introduce the cross terms and be strictly
/// worse — the symmetry of the SQUARE is what buys the win.
///
/// Peak control: the (lo+hi)^2 square is emitted FIRST, before the 2n-bit
/// `tmp_ext` result register is allocated, and its `x_sum` operand is freed
/// before `tmp_ext` is taken — so the z1 step (z1_reg + x_sum + row) and the
/// z0/z2 step (tmp_ext + z1_reg + row) never coexist. The combine carries use
/// the non-fast (ancilla-free) Cuccaro, and the Solinas lanes default to the
/// low-peak set (non-fast add/sub, direct-const double/halve, lowq shift) so the
/// extra z1_reg register (2(h+1) q) is absorbed without pushing the affine
/// square phase over the global GCD-body peak binder (~1567 < 1698).
pub(crate) fn squaring_sub_from_acc_karatsuba(b: &mut B, acc: &[QubitId], x: &[QubitId], p: U256) {
    let n = acc.len();
    debug_assert_eq!(n, 256);
    debug_assert_eq!(x.len(), n);
    let h = n / 2;
    let x_lo: Vec<QubitId> = x[0..h].to_vec();
    let x_hi: Vec<QubitId> = x[h..n].to_vec();

    // z1_reg holds z1 = (lo+hi)^2, width 2*(h+1).
    let mut z1_reg = b.alloc_qubits(2 * (h + 1));
    // KARA_FREE_Z1_TOPBIT: after z1 -= z0; z1 -= z2, z1_reg holds 2*lo*hi < 2^257,
    // so its top bit (index 2(h+1)-1 = 257) is provably 0 throughout the Solinas
    // peak. Free it for that window; re-grab a fresh zero before z1 += z2 restores
    // (lo+hi)^2 for the inverse uncompute. Bennett-clean (free zero, alloc zero).
    let free_z1_top = std::env::var("KARA_FREE_Z1_TOPBIT").ok().as_deref() == Some("1");
    // The z0=lo^2 / z2=hi^2 squares coexist with tmp_ext(2n)+z1_reg, and the
    // _fast symmetric square allocates a ~(h)-wide cuccaro carry lane on top of
    // its ~(h)-wide row — that lane is the round84 peak binder. The ancilla-free
    // _lowq square drops the carry lane (peak −~h) at a higher Toffoli cost.
    // z1=(lo+hi)^2 is computed before tmp_ext (low peak), so it stays _fast.
    let z02_lowq = std::env::var("KARA_Z02_LOWQ").ok().as_deref() == Some("1");

    // ── Forward z1 = (lo+hi)^2 FIRST (tmp_ext not yet allocated → low peak). ──
    {
        let x_sum = b.alloc_qubits(h + 1);
        karatsuba_half_sum_compute(b, &x_lo, &x_hi, &x_sum);
        schoolbook_square_symmetric(b, &x_sum, &z1_reg);
        karatsuba_half_sum_uncompute(b, &x_lo, &x_hi, &x_sum);
        b.free_vec(&x_sum);
    }

    // 2n-bit result accumulator for x^2 (allocated after the z1 square so its
    // 2n qubits never coexist with the z1 operand/row registers).
    let tmp_ext = b.alloc_qubits(2 * n);

    // z0 = lo^2 → tmp_ext[0..2h], z2 = hi^2 → tmp_ext[2h..4h].
    {
        let slice: Vec<QubitId> = tmp_ext[0..2 * h].to_vec();
        if z02_lowq {
            // z2 slice (tmp_ext[2h..4h]) is still clean here → host z0's fast
            // carry there (Toffoli-free peak drop) instead of paying lowq.
            let host: Vec<QubitId> = tmp_ext[2 * h..4 * h].to_vec();
            schoolbook_square_symmetric_hosted(b, &x_lo, &slice, &host);
        } else {
            schoolbook_square_symmetric(b, &x_lo, &slice);
        }
    }
    {
        let slice: Vec<QubitId> = tmp_ext[2 * h..4 * h].to_vec();
        if z02_lowq {
            if kara_z2_selfhost_enabled() {
                if square_selfhost_safe_lane_reuse_enabled() {
                    // z1=(lo+hi)^2 and z0=lo^2 are exact integer squares here.
                    // Every square is 0 or 1 mod 4, so bit 1 of each register is
                    // provably |0>.  Both lanes are disjoint from x_hi, z2, and
                    // z2's own untouched-tail carry lanes.
                    let clean_square_bits = [z1_reg[1], tmp_ext[1]];
                    schoolbook_square_symmetric_lowq_selfhosted_with_clean_supplement(
                        b,
                        &x_hi,
                        &slice,
                        &clean_square_bits,
                    );
                } else {
                    schoolbook_square_symmetric_lowq_selfhosted(b, &x_hi, &slice);
                }
            } else {
                schoolbook_square_symmetric_lowq(b, &x_hi, &slice);
            }
        } else {
            schoolbook_square_symmetric(b, &x_hi, &slice);
        }
    }

    // Combine: z1 -= z0; z1 -= z2; mid (tmp_ext[h..4h]) += z1. Non-fast Cuccaro
    // (no carry ancilla) keeps the peak flat while tmp_ext + z1_reg are live.
    {
        let pad = b.alloc_qubits(2);
        let mut z0_ext: Vec<QubitId> = tmp_ext[0..2 * h].to_vec();
        z0_ext.extend_from_slice(&pad);
        sub_nbit_qq(b, &z0_ext, &z1_reg);
        b.free_vec(&pad);
    }
    {
        let pad = b.alloc_qubits(2);
        let mut z2_ext: Vec<QubitId> = tmp_ext[2 * h..4 * h].to_vec();
        z2_ext.extend_from_slice(&pad);
        sub_nbit_qq(b, &z2_ext, &z1_reg);
        b.free_vec(&pad);
    }
    // z1_reg == 2*lo*hi < 2^257 here ⇒ bit 257 is 0. Release it for the peak window.
    if free_z1_top {
        let top = z1_reg.pop().expect("z1_reg width 2*(h+1) >= 2");
        b.free(top);
    }
    {
        let pad = b.alloc_qubits(3 * h - z1_reg.len());
        let mut z1_ext: Vec<QubitId> = z1_reg.to_vec();
        z1_ext.extend_from_slice(&pad);
        let acc_slice: Vec<QubitId> = tmp_ext[h..4 * h].to_vec();
        add_nbit_qq(b, &z1_ext, &acc_slice);
        b.free_vec(&pad);
    }

    // ── Solinas reduction: acc -= (lo + hi·c) mod p. ──
    // z1_reg (2(h+1) q) is still live through this whole block, so the lanes
    // that allocate a full-width carry ancilla (fast Cuccaro add/sub, fast
    // shift) bind the affine-square phase peak. Each lane defaults to its
    // low-peak (ancilla-free) variant so the phase peak stays below the global
    // GCD-body binder; per-lane env knobs select the higher-peak fast variants
    // for measurement (each computes the SAME value on `acc`, so any mix is
    // value-correct):
    //   KARA_SOL_MOD_FAST=1   → fast mod add/sub          (else non-fast)
    //   KARA_SOL_DBL_FAST=1   → fast in-place double/halve (else direct-const)
    //   KARA_SOL_SHIFT_FAST=1 → fast shift-by-22          (else lowq shift)
    let mod_fast = std::env::var("KARA_SOL_MOD_FAST").ok().as_deref() == Some("1");
    let dbl_fast = std::env::var("KARA_SOL_DBL_FAST").ok().as_deref() == Some("1");
    let shift_fast = std::env::var("KARA_SOL_SHIFT_FAST").ok().as_deref() == Some("1");
    let lo: Vec<QubitId> = tmp_ext[0..n].to_vec();
    let hi: Vec<QubitId> = tmp_ext[n..2 * n].to_vec();
    // The non-fast mod_add/sub materialize a 256-q load_const for the Solinas
    // `c` correction, which coexists with tmp_ext + z1_reg and binds the phase
    // peak. The vent form hosts that correction on the operand `a_ext` (dirty,
    // value-preserved) for 2 clean qubits, dropping the transient ~n.
    let mod_vent = std::env::var("KARA_SOL_MOD_VENT").ok().as_deref() == Some("1");
    let mod_sub = |b: &mut B, acc: &[QubitId], a: &[QubitId]| {
        if mod_vent {
            mod_sub_qq_vent(b, acc, a, p);
        } else if mod_fast {
            mod_sub_qq_fast(b, acc, a, p);
        } else {
            mod_sub_qq(b, acc, a, p);
        }
    };
    let mod_add = |b: &mut B, acc: &[QubitId], a: &[QubitId]| {
        if mod_vent {
            mod_add_qq_vent(b, acc, a, p);
        } else if mod_fast {
            mod_add_qq_fast(b, acc, a, p);
        } else {
            mod_add_qq(b, acc, a, p);
        }
    };
    let mod_dbl = |b: &mut B, v: &[QubitId]| {
        if dbl_fast {
            mod_double_inplace_fast(b, v, p);
        } else {
            mod_double_inplace_direct_const_fast(b, v, p);
        }
    };
    let mod_hlv = |b: &mut B, v: &[QubitId]| {
        if dbl_fast {
            mod_halve_inplace_fast(b, v, p);
        } else {
            mod_halve_inplace_direct_const_fast(b, v, p);
        }
    };
    b.set_phase("r84k_sol_subadd");
    mod_sub(b, acc, &lo);
    mod_sub(b, acc, &hi);
    for _ in 0..4 {
        mod_dbl(b, &hi);
    }
    mod_sub(b, acc, &hi);
    for _ in 0..2 {
        mod_dbl(b, &hi);
    }
    mod_add(b, acc, &hi); // sign flipped
    for _ in 0..4 {
        mod_dbl(b, &hi);
    }
    mod_sub(b, acc, &hi);
    b.set_phase("r84k_sol_shift");
    // The shift-by-22 lane binds the affine-square phase peak: its lowq form
    // allocates a ~(n+1)-wide `padded` scratch on top of the live z1_reg+tmp_ext,
    // overflowing the free pool. `acc` (tx) is idle and value-preserved during the
    // shift itself, so the dirty-borrow form hosts that scratch on `acc` (venting
    // 2-clean), dropping the phase peak well under the GCD-apply binder. Same value
    // on `acc`; gated so it can be A/B compared.
    let shift_dirty = std::env::var("ROUND84_XTAIL_BORROW_CARRIES")
        .ok()
        .as_deref()
        == Some("1");
    if shift_dirty {
        // Dirty-doubles form of `acc -= hi * 2^22 mod p`: 22 in-place doubles
        // (each borrows `acc` via Gidney venting) avoid the shift's persistent
        // k-wide `spill` lane that — stacked on the live z1_reg+tmp_ext base —
        // pushed the shift/mid-sub over the GCD-apply binder. `acc` is idle and
        // value-preserved during each double/halve, so the phase peak drops well
        // under 1558. Mirrors the schoolbook_peak_lowq D1 reduction lane.
        b.set_phase("r84k_sol_dbl22");
        for _ in 0..22 {
            mod_dbl(b, &hi);
        }
        b.set_phase("r84k_sol_midsub");
        mod_sub(b, acc, &hi);
        b.set_phase("r84k_sol_hlv22");
        for _ in 0..22 {
            mod_hlv(b, &hi);
        }
    } else {
        b.set_phase("r84k_sol_shiftL");
        let (spill, flag_inv, ovf) = if shift_fast {
            mod_shift_left_by_k(b, &hi, p, 22)
        } else {
            mod_shift_left_by_k_lowq(b, &hi, p, 22)
        };
        b.set_phase("r84k_sol_midsub");
        mod_sub(b, acc, &hi);
        b.set_phase("r84k_sol_shiftR");
        if shift_fast {
            mod_shift_right_by_k(b, &hi, p, 22, spill, flag_inv, ovf);
        } else {
            mod_shift_right_by_k_lowq(b, &hi, p, 22, spill, flag_inv, ovf);
        }
    }
    b.set_phase("r84k_sol_halve");
    for _ in 0..10 {
        mod_hlv(b, &hi);
    }

    // ── Inverse combine: mid -= z1; z1 += z2; z1 += z0. ──
    b.set_phase("r84k_inv_combine");
    {
        let pad = b.alloc_qubits(3 * h - z1_reg.len());
        let mut z1_ext: Vec<QubitId> = z1_reg.to_vec();
        z1_ext.extend_from_slice(&pad);
        let acc_slice: Vec<QubitId> = tmp_ext[h..4 * h].to_vec();
        sub_nbit_qq(b, &z1_ext, &acc_slice);
        b.free_vec(&pad);
    }
    // Restore z1_reg top bit (fresh zero) before z1 += z2 can re-set it.
    if free_z1_top {
        let top = b.alloc_qubit();
        z1_reg.push(top);
    }
    {
        let pad = b.alloc_qubits(2);
        let mut z2_ext: Vec<QubitId> = tmp_ext[2 * h..4 * h].to_vec();
        z2_ext.extend_from_slice(&pad);
        add_nbit_qq(b, &z2_ext, &z1_reg);
        b.free_vec(&pad);
    }
    {
        let pad = b.alloc_qubits(2);
        let mut z0_ext: Vec<QubitId> = tmp_ext[0..2 * h].to_vec();
        z0_ext.extend_from_slice(&pad);
        add_nbit_qq(b, &z0_ext, &z1_reg);
        b.free_vec(&pad);
    }

    // Uncompute z2, z0 (reverse of forward compute order), then free tmp_ext.
    b.set_phase("r84k_z_inv_squares");
    {
        let slice: Vec<QubitId> = tmp_ext[2 * h..4 * h].to_vec();
        if z02_lowq {
            if kara_z2_selfhost_enabled() {
                if square_selfhost_safe_lane_reuse_enabled() {
                    // Inverse-combine restored the exact z1 and z0 squares
                    // before this block, so their square-bit-1 lanes are clean
                    // scratch again (the mirror of the forward z2 proof).
                    let clean_square_bits = [z1_reg[1], tmp_ext[1]];
                    schoolbook_square_symmetric_lowq_selfhosted_inverse_with_clean_supplement(
                        b,
                        &x_hi,
                        &slice,
                        &clean_square_bits,
                    );
                } else {
                    schoolbook_square_symmetric_lowq_selfhosted_inverse(b, &x_hi, &slice);
                }
            } else {
                schoolbook_square_symmetric_lowq_inverse(b, &x_hi, &slice);
            }
        } else {
            schoolbook_square_symmetric_inverse(b, &x_hi, &slice);
        }
    }
    {
        let slice: Vec<QubitId> = tmp_ext[0..2 * h].to_vec();
        if z02_lowq {
            // z2 slice was just uncomputed above → clean again, host inv-z0's
            // borrow there (mirror of the forward z0 hosting).
            let host: Vec<QubitId> = tmp_ext[2 * h..4 * h].to_vec();
            schoolbook_square_symmetric_hosted_inverse(b, &x_lo, &slice, &host);
        } else {
            schoolbook_square_symmetric_inverse(b, &x_lo, &slice);
        }
    }
    b.free_vec(&tmp_ext);

    // Uncompute z1 last (mirrors the forward z1-first ordering, tmp_ext freed).
    {
        let x_sum = b.alloc_qubits(h + 1);
        karatsuba_half_sum_compute(b, &x_lo, &x_hi, &x_sum);
        schoolbook_square_symmetric_inverse(b, &x_sum, &z1_reg);
        karatsuba_half_sum_uncompute(b, &x_lo, &x_hi, &x_sum);
        b.free_vec(&x_sum);
    }

    b.free_vec(&z1_reg);
}

pub(crate) fn squaring_sub_from_acc_schoolbook_lowq_shift22(
    b: &mut B,
    acc: &[QubitId],
    x: &[QubitId],
    p: U256,
) {
    let n = acc.len();
    debug_assert_eq!(n, 256);
    debug_assert_eq!(x.len(), n);

    let tmp_ext = b.alloc_qubits(2 * n);
    b.set_phase("round84_inplace_solinas_square_forward");
    if xtail_sq_selfhost_enabled() {
        schoolbook_square_symmetric_lowq_selfhosted(b, x, &tmp_ext);
    } else {
        schoolbook_square_symmetric_lowq(b, x, &tmp_ext);
    }

    let lo: Vec<QubitId> = tmp_ext[0..n].to_vec();
    let hi: Vec<QubitId> = tmp_ext[n..2 * n].to_vec();
    if round84_inplace_solinas_fold_enabled() {
        b.set_phase("round84_inplace_solinas_fold");
        let state = round84_fold_hi_into_lo_aggregate(b, &lo, &hi, acc);
        b.set_phase("round84_inplace_solinas_sub");
        mod_sub_qq_vent(b, acc, &lo, p);
        b.set_phase("round84_inplace_solinas_unfold");
        round84_unfold_hi_from_lo_aggregate(b, &lo, &hi, acc, state);
    } else {
        mod_sub_qq(b, acc, &lo, p);
        mod_sub_qq(b, acc, &hi, p);
        for _ in 0..4 {
            mod_double_inplace_direct_const_fast(b, &hi, p);
        }
        mod_sub_qq(b, acc, &hi, p);
        for _ in 0..2 {
            mod_double_inplace_direct_const_fast(b, &hi, p);
        }
        mod_add_qq(b, acc, &hi, p);
        for _ in 0..4 {
            mod_double_inplace_direct_const_fast(b, &hi, p);
        }
        mod_sub_qq(b, acc, &hi, p);
        let (spill, flag_inv, ovf) = mod_shift_left_by_k_lowq(b, &hi, p, 22);
        if r84_lowq_enabled() {
            mod_sub_qq_lowq(b, acc, &hi, p);
        } else {
            mod_sub_qq(b, acc, &hi, p);
        }
        mod_shift_right_by_k_lowq(b, &hi, p, 22, spill, flag_inv, ovf);
        for _ in 0..10 {
            mod_halve_inplace_direct_const_fast(b, &hi, p);
        }
    }

    b.set_phase("round84_inplace_solinas_square_inverse");
    if xtail_sq_selfhost_enabled() {
        schoolbook_square_symmetric_lowq_selfhosted_inverse(b, x, &tmp_ext);
    } else {
        schoolbook_square_symmetric_lowq_inverse(b, x, &tmp_ext);
    }
    b.free_vec(&tmp_ext);
}

pub(crate) fn squaring_sub_from_acc_walk_controls_lowq(
    b: &mut B,
    acc: &[QubitId],
    x: &[QubitId],
    p: U256,
) {
    let n = acc.len();
    debug_assert_eq!(n, 256);
    debug_assert_eq!(x.len(), n);

    let ctrl_copy = b.alloc_qubits(n);
    for i in 0..n {
        b.cx(x[i], ctrl_copy[i]);
    }

    mod_neg_inplace_fast(b, x, p);
    for i in 0..n {
        cmod_add_qq(b, acc, x, ctrl_copy[i], p);
        if i < n - 1 {
            mod_double_inplace_fast(b, x, p);
        }
    }
    for _ in 0..(n - 1) {
        mod_halve_inplace_fast(b, x, p);
    }
    mod_neg_inplace_fast(b, x, p);

    for i in 0..n {
        b.cx(x[i], ctrl_copy[i]);
    }
    b.free_vec(&ctrl_copy);
}

// HYP-12 lever (the round84 Solinas fold/unfold wall @1221). The fold's
// quotient*c product is built by shifted adds of the 33-bit quotient `q` into
// the 66-bit `product`. To match widths the caller zero-extends `q` with a
// `pad` whose width is `product.len()-shift-33` (=29 at shift=4). MEASURED:
// the shift=4 pad (29 transient |0> lanes) is the SOLE binder that pins the
// fold phase at 1221 (UNQ_sh4=1221, the next is UNQ_sh5=1219). The high `pad`
// bits are all 0, so the work on `product[shift+33..]` is a pure carry ripple,
// not a real add. This lever replaces the 29-lane pad with a single carry
// `wrap` + a Gidney measure-vented carry ripple (`ciadd/cisub_dirty_2clean`,
// borrowing the idle `acc`/`dirty` lanes + 2 clean, uncompute=0) and an
// ancilla-free `cmp_lt_into` wrap-uncompute (1 c_in, ~n CCX, no carry array).
// Net: the qprod transient drops from product+~30 to product+~3 => the fold
// peak falls below 1221, exposing the global drop to 1220 (with SEG<=193 the
// square is already <=1220 there). Value-exact (a permutation that round-
// trips); default OFF (byte-identical base).
fn round84_qprod_vent_pad_enabled() -> bool {
    std::env::var("ROUND84_QPROD_VENT_PAD").ok().as_deref() == Some("1")
}

// Only the WIDEST-pad shifted add binds the fold peak (MEASURED: shift=4 pins
// 1221, shift=5 sits at 1219). Venting the narrower-pad shifts adds Toffoli
// (a cmp_lt wrap-uncompute + a vented ripple) for no peak gain, so by default
// the lever only vents shifts whose pad width exceeds this threshold. The
// shift=4 pad is `66-4-33 = 29`; shift=5 is 28; set the cutoff at 29 so only
// shift=4 vents. Override with ROUND84_QPROD_VENT_PAD_MINW to vent more.
fn round84_qprod_vent_pad_min_width() -> usize {
    std::env::var("ROUND84_QPROD_VENT_PAD_MINW")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(29)
}

// NOTE (measured): venting must cover BOTH the fold's qprod-uncompute AND the
// unfold's qprod-compute — each builds the 66-bit product with the 29-lane
// shift=4 pad and reaches 1221 (the fold-compute @1220 and unfold-uncompute
// @1220 are 1 below, but the round-trip needs all four product builds vented
// to clear the wall). The lever therefore vents every shift>=MINW build in
// both compute and uncompute.

/// One shifted small-add of `q` (33-bit) into `product[shift..]`, with the
/// high zero-extension realized as a vented carry ripple instead of a `pad`.
/// `add=true` => `product[shift..] += q`; `add=false` => `-= q`. The carry/
/// borrow `wrap` is recomputed-and-freed in place (no residue), so this is the
/// exact width-matched equivalent of the padded `cuccaro_add/sub` it replaces.
fn round84_qprod_shifted_addsub_vented(
    b: &mut B,
    q: &[QubitId],
    product: &[QubitId],
    shift: usize,
    add: bool,
    dirty: &[QubitId],
) {
    let m = q.len(); // 33
    let total = product.len() - shift;
    debug_assert!(total >= m);
    let high_w = total - m; // width of the zero-extension (carry ripple region)

    // The vent helper needs n>4 dirty/clean lanes; for tiny tails fall back to
    // the padded coherent add (these shifts never bind the peak).
    if high_w < 5 || dirty.len() < high_w.saturating_sub(2) {
        let target = &product[shift..];
        let pad = b.alloc_qubits(high_w);
        let mut source = q.to_vec();
        source.extend_from_slice(&pad);
        if add {
            round84_add_small(b, &source, target);
        } else {
            round84_sub_small(b, &source, target);
        }
        b.free_vec(&pad);
        return;
    }

    let wrap = b.alloc_qubit();
    let mut low_ext = product[shift..shift + m].to_vec();
    low_ext.push(wrap);
    let high = &product[shift + m..];

    if add {
        // product[shift..shift+m] += q, carry-out -> wrap.
        let c_in = b.alloc_qubit();
        cuccaro_add_low_to_ext_clean(b, q, &low_ext, c_in);
        b.free(c_in);
        // Ripple the carry: product[shift+m..] += wrap (vented, dirty-borrowed).
        let clean2 = [b.alloc_qubit(), b.alloc_qubit()];
        venting::ciadd_dirty_2clean_classical(
            b,
            high,
            &dirty[..high_w - 2],
            &clean2,
            1,
            wrap,
            false,
        );
        b.free(clean2[1]);
        b.free(clean2[0]);
        // Uncompute wrap: carry == (new_low < q). cmp_lt_into uses 1 c_in only.
        cmp_lt_into(b, &product[shift..shift + m], q, wrap);
    } else {
        // product[shift..shift+m] -= q, borrow-out -> wrap.
        let c_in = b.alloc_qubit();
        cuccaro_sub_low_to_ext_clean(b, q, &low_ext, c_in);
        b.free(c_in);
        // Ripple the borrow: product[shift+m..] -= wrap (vented).
        let clean2 = [b.alloc_qubit(), b.alloc_qubit()];
        venting::cisub_dirty_2clean_classical(b, high, &dirty[..high_w - 2], &clean2, 1, wrap);
        b.free(clean2[1]);
        b.free(clean2[0]);
        // Uncompute wrap: borrow == carry_out(new_low + q) == (~q < new_low).
        for &qb in q {
            b.x(qb);
        }
        cmp_lt_into(b, q, &product[shift..shift + m], wrap);
        for &qb in q {
            b.x(qb);
        }
    }
    b.free(wrap);
}
