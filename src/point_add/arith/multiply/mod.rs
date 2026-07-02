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


mod squaring;
pub(crate) use squaring::*;
