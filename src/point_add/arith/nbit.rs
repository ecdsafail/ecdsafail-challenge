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

#[allow(dead_code)] // retained reference/alternative impl; not on active build path
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

#[allow(dead_code)] // retained reference/alternative impl; not on active build path
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


