//! Bit-shift primitives: logical left/right shift by one position.
//! Extracted from the former `mbu_primitives` / `poc_arith` files.

use crate::point_add::trailmix_port::circuit::{Circuit, QReg};

/// Left shift by 1: a <<= 1. Implemented as rotation (MSB wraps to
/// LSB), so the caller MUST ensure a[n-1] == |0> before calling.
/// The nw = nb+1 invariant guarantees this for all `mod_mul` operands.
pub fn left_shift(circ: &mut Circuit, a: &[QReg]) {
    let n = a.len();
    for i in (1..n).rev() {
        circ.swap(&a[i], &a[i - 1]);
    }
}

pub fn right_shift(circ: &mut Circuit, a: &[QReg]) {
    let n = a.len();
    for i in 0..n - 1 {
        circ.swap(&a[i], &a[i + 1]);
    }
}
