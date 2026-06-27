use super::*;

pub(crate) fn add_nbit_qq_fast(b: &mut B, a: &[QubitId], acc: &[QubitId]) {
    assert_eq!(a.len(), acc.len());
    let c_in = b.alloc_qubit();
    cuccaro_add_fast(b, a, acc, c_in);
    b.free(c_in);
}

/// Fast `acc -= a mod 2^n` using measurement-based Cuccaro.
pub(crate) fn sub_nbit_qq_fast(b: &mut B, a: &[QubitId], acc: &[QubitId]) {
    assert_eq!(a.len(), acc.len());
    let c_in = b.alloc_qubit();
    cuccaro_sub_fast(b, a, acc, c_in);
    b.free(c_in);
}

pub(crate) fn add_nbit_qq_fast_borrowed_carries(
    b: &mut B,
    a: &[QubitId],
    acc: &[QubitId],
    carries: &[QubitId],
) {
    assert_eq!(a.len(), acc.len());
    let c_in = b.alloc_qubit();
    cuccaro_add_fast_borrowed_carries(b, a, acc, c_in, carries);
    b.free(c_in);
}

pub(crate) fn sub_nbit_qq_fast_borrowed_carries(
    b: &mut B,
    a: &[QubitId],
    acc: &[QubitId],
    carries: &[QubitId],
) {
    assert_eq!(a.len(), acc.len());
    let c_in = b.alloc_qubit();
    cuccaro_sub_fast_borrowed_carries(b, a, acc, c_in, carries);
    b.free(c_in);
}

#[inline]
fn maj3_into_clean_2ccx(b: &mut B, x: QubitId, y: QubitId, z: QubitId, target: QubitId) {
    debug_assert!(x != y && x != z && x != target && y != z && y != target && z != target);
    b.ccx(x, z, target);
    b.cx(x, z);
    b.ccx(y, z, target);
    b.cx(x, z);
}

/// Exact measured add of a short source into a longer accumulator, without
/// materializing the zero-valued high suffix of the source.
pub(crate) fn add_short_to_long_qq_fast_no_cin(b: &mut B, a: &[QubitId], acc: &[QubitId]) {
    let m = a.len();
    let n = acc.len();
    assert!(m > 0);
    assert!(m <= n);
    if n == 1 {
        b.cx(a[0], acc[0]);
        return;
    }

    let carries = b.alloc_qubits(n - 1);
    for i in 0..n - 1 {
        if i < m {
            if i == 0 {
                b.ccx(acc[i], a[i], carries[i]);
            } else {
                maj3_into_clean_2ccx(b, acc[i], a[i], carries[i - 1], carries[i]);
            }
        } else {
            b.ccx(acc[i], carries[i - 1], carries[i]);
        }
    }

    for i in 0..n {
        if i < m {
            b.cx(a[i], acc[i]);
        }
        if i > 0 {
            b.cx(carries[i - 1], acc[i]);
        }
    }

    for i in (0..n - 1).rev() {
        let bit = b.alloc_bit();
        b.hmr(carries[i], bit);
        if i < m {
            b.x(acc[i]);
            b.cz_if(acc[i], a[i], bit);
            if i > 0 {
                b.cz_if(acc[i], carries[i - 1], bit);
                b.x(acc[i]);
                b.cz_if(a[i], carries[i - 1], bit);
            } else {
                b.x(acc[i]);
            }
        } else {
            b.x(acc[i]);
            b.cz_if(acc[i], carries[i - 1], bit);
            b.x(acc[i]);
        }
    }
    b.free_vec(&carries);
}

/// Exact measured subtract of a short source from a longer accumulator, without
/// materializing the zero-valued high suffix of the source.
pub(crate) fn sub_short_to_long_qq_fast_no_cin(b: &mut B, a: &[QubitId], acc: &[QubitId]) {
    let m = a.len();
    let n = acc.len();
    assert!(m > 0);
    assert!(m <= n);
    if n == 1 {
        b.cx(a[0], acc[0]);
        return;
    }

    let borrows = b.alloc_qubits(n - 1);
    for i in 0..n - 1 {
        if i < m {
            b.x(acc[i]);
            if i == 0 {
                b.ccx(acc[i], a[i], borrows[i]);
            } else {
                maj3_into_clean_2ccx(b, acc[i], a[i], borrows[i - 1], borrows[i]);
            }
            b.x(acc[i]);
        } else {
            b.x(acc[i]);
            b.ccx(acc[i], borrows[i - 1], borrows[i]);
            b.x(acc[i]);
        }
    }

    for i in 0..n {
        if i < m {
            b.cx(a[i], acc[i]);
        }
        if i > 0 {
            b.cx(borrows[i - 1], acc[i]);
        }
    }

    for i in (0..n - 1).rev() {
        let bit = b.alloc_bit();
        b.hmr(borrows[i], bit);
        if i < m {
            b.cz_if(acc[i], a[i], bit);
            if i > 0 {
                b.cz_if(acc[i], borrows[i - 1], bit);
                b.cz_if(a[i], borrows[i - 1], bit);
            }
        } else {
            b.cz_if(acc[i], borrows[i - 1], bit);
        }
    }
    b.free_vec(&borrows);
}

/// `acc += a mod 2^n`. Caller must pre-extend both slices if they want the
/// top carry absorbed into the accumulator (i.e. pass n+1-bit slices with
/// top bits 0 to get a full n+1-bit add). The carry-out beyond the slice
/// is discarded via `R` on the `z` ancilla — safe when both inputs fit
/// in n-1 bits (as in our mod-p layer where both < 2p < 2^{n+1}).
pub(crate) fn add_nbit_qq(b: &mut B, a: &[QubitId], acc: &[QubitId]) {
    assert_eq!(a.len(), acc.len());
    let c_in = b.alloc_qubit();
    cuccaro_add(b, a, acc, c_in);
    b.free(c_in);
}

pub(crate) fn sub_nbit_qq(b: &mut B, a: &[QubitId], acc: &[QubitId]) {
    assert_eq!(a.len(), acc.len());
    let c_in = b.alloc_qubit();
    cuccaro_sub(b, a, acc, c_in);
    b.free(c_in);
}

pub(crate) fn add_nbit_const(b: &mut B, acc: &[QubitId], c: U256) {
    let n = acc.len();
    let a = load_const(b, n, c);
    add_nbit_qq(b, &a, acc);
    unload_const(b, &a, c);
}

pub(crate) fn sub_nbit_const(b: &mut B, acc: &[QubitId], c: U256) {
    let n = acc.len();
    let a = load_const(b, n, c);
    sub_nbit_qq(b, &a, acc);
    unload_const(b, &a, c);
}
