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
pub(crate) fn controlled_add_subtract_lowq(b: &mut B, x: &[QubitId], acc: &[QubitId], ctrl: QubitId) {
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
pub(crate) fn controlled_add_subtract_lowq_inverse(b: &mut B, x: &[QubitId], acc: &[QubitId], ctrl: QubitId) {
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
pub(crate) fn schoolbook_mul_into_addsub_lowq(b: &mut B, x: &[QubitId], y: &[QubitId], tmp_ext: &[QubitId]) {
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

pub(crate) fn karatsuba_half_sum_compute(b: &mut B, lo: &[QubitId], hi: &[QubitId], acc: &[QubitId]) {
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

pub(crate) fn karatsuba_half_sum_uncompute(b: &mut B, lo: &[QubitId], hi: &[QubitId], acc: &[QubitId]) {
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

pub(crate) fn schoolbook_square_symmetric_lowq_inverse(b: &mut B, x: &[QubitId], tmp_ext: &[QubitId]) {
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
pub(crate) fn schoolbook_square_symmetric_lowq_selfhosted(b: &mut B, x: &[QubitId], tmp_ext: &[QubitId]) {
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
    for i in 0..n {
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let num_cross = if i + 1 < n { n - i - 1 } else { 0 };
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
    for i in (0..n).rev() {
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let num_cross = if i + 1 < n { n - i - 1 } else { 0 };
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
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));

    let tmp_ext = b.alloc_qubits(2 * n);
    if xtail_sq_selfhost_enabled() {
        schoolbook_square_symmetric_lowq_selfhosted(b, x, &tmp_ext);
    } else {
        schoolbook_square_symmetric_lowq(b, x, &tmp_ext);
    }

    let lo: Vec<QubitId> = tmp_ext[0..n].to_vec();
    let hi: Vec<QubitId> = tmp_ext[n..2 * n].to_vec();
    mod_sub_qq(b, acc, &lo, p);
    let _ = c;
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

    if xtail_sq_selfhost_enabled() {
        schoolbook_square_symmetric_lowq_selfhosted_inverse(b, x, &tmp_ext);
    } else {
        schoolbook_square_symmetric_lowq_inverse(b, x, &tmp_ext);
    }
    b.free_vec(&tmp_ext);
}

pub(crate) fn squaring_sub_from_acc_walk_controls_lowq(b: &mut B, acc: &[QubitId], x: &[QubitId], p: U256) {
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

