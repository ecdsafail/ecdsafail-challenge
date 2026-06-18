//! In-place log-depth barrel shifter: shift a quantum register `b` by a
//! quantum amount `s`, one cswap layer per bit of `s` (layer i shifts by 2^i
//! when `s[i] = 1`). Used by the shrunken-PZ divstep to align the cofactor
//! registers.
//!
//! Precondition: the top `s_max` bits of `b` must be |0> on entry (where
//! `s_max = 2^len(s) - 1`); otherwise high bits shift off the top of the
//! in-place register and `b` is not restored by the reverse shifter.

use crate::point_add::trailmix_port::circuit::{Circuit, QReg};

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
        if k >= n {
            // Whole register would shift off-end; nothing to do
            // (precondition guarantees those bits are 0).
            continue;
        }
        // cswap pairs (j, j-k) for j = n-1 down to k.
        // Forward: top-to-bottom; reverse: bottom-to-top.
        let mut pairs: Vec<(usize, usize)> = ((k..n).rev()).map(|j| (j, j - k)).collect();
        if !forward {
            pairs.reverse();
        }
        for (hi, lo) in pairs {
            // cswap(s[i], b[hi], b[lo]) via Fredkin = CX-CCX-CX.
            circ.cx(&b[lo], &b[hi]);
            circ.ccx(&s[i], &b[hi], &b[lo]);
            circ.cx(&b[lo], &b[hi]);
        }
    }
    circ.pop_section(&prev);
}
