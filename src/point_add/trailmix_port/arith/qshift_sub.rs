//! In-place log-depth barrel shifter: shift a quantum register `b` by a
//! quantum amount `s`, one cswap layer per bit of `s` (layer i shifts by 2^i
//! when `s[i] = 1`). Used by the shrunken-PZ divstep to align the cofactor
//! registers.
//!
//! Precondition: the top `s_max` bits of `b` must be |0> on entry (where
//! `s_max = 2^len(s) - 1`); otherwise high bits shift off the top of the
//! in-place register and `b` is not restored by the reverse shifter.

use crate::point_add::trailmix_port::circuit::{Circuit, QReg};

fn controlled_shift_layer(
    circ: &mut Circuit,
    reg: &[QReg],
    control: &QReg,
    distance: usize,
    forward: bool,
) {
    let n = reg.len();
    if distance == 0 || distance >= n {
        return;
    }
    let mut pairs: Vec<(usize, usize)> =
        ((distance..n).rev()).map(|j| (j, j - distance)).collect();
    if !forward {
        pairs.reverse();
    }
    for (hi, lo) in pairs {
        circ.cx(&reg[lo], &reg[hi]);
        circ.ccx(control, &reg[hi], &reg[lo]);
        circ.cx(&reg[lo], &reg[hi]);
    }
}

/// Shift by one fixed classical distance under a quantum control. `forward`
/// moves bits toward higher indices; the opposite orientation is its inverse.
pub fn controlled_shift_fixed_inplace(
    circ: &mut Circuit,
    reg: &[QReg],
    control: &QReg,
    distance: usize,
    forward: bool,
) {
    let prev = circ.push_section("p.shift");
    controlled_shift_layer(circ, reg, control, distance, forward);
    circ.pop_section(&prev);
}

/// In-place barrel shift `b <<= s` (toward higher indices) when `forward`, or
/// the exact inverse when `!forward` (for uncomputation). Each layer is a
/// Fredkin (CX-CCX-CX) cswap per affected position; no ancillae.
pub fn barrel_shift_inplace(circ: &mut Circuit, b: &[QReg], s: &[QReg], forward: bool) {
    let n = b.len();
    if n == 0 || s.is_empty() {
        return;
    }
    let prev = circ.push_section("p.shift");
    let layer_order: Vec<usize> = if forward {
        (0..s.len()).collect()
    } else {
        (0..s.len()).rev().collect()
    };
    for &i in &layer_order {
        let k = 1usize << i;
        controlled_shift_layer(circ, b, &s[i], k, forward);
    }
    circ.pop_section(&prev);
}
