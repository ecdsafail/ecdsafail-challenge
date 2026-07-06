//! Top-level secp256k1 affine point-add driver. It composes the `arith` /
//! `gcd` / `square` primitives into the full in-place point addition, built on
//! this crate's `B` builder.
//!
//! It adds a quantum point P to a classical point Q in place, computing
//! `(x2, y2) -> P + Q` via the a-independent secp add formula:
//!
//!   3:  x2 -= ox                 (coordinate const-subtract)
//!   4:  y2 -= oy                 (coordinate const-subtract)
//!   6:  y2 *= x2^-1 mod q        (gcd inversion: y2 becomes the slope lambda)
//!   7:  x2 += 3*ox               (coordinate const-add of 3*ox)
//!  10:  x2 -= lambda^2 mod q     (symmetric mod-square-subtract)
//!  11:  y2 *= x2   mod q         (gcd forward multiply)
//!  14:  y2 -= oy                 (coordinate const-subtract)
//!  15:  x2 := ox - x2  mod q     (mod-negate then const-add ox)
//!
//! Post: `(x2, y2)` holds the affine `P + Q` (mod q, on the common path).
//!
//! Steps 6 and 11 are the two GCD passes that share a single Schrottenloher
//! jump-GCD inversion.
//!
//! ## Register convention
//! - `x2`: the quantum point's x-coordinate working register. 256 data bits.
//!   `mod_mul_inverse_in_place` consumes and returns it (restored to dx after
//!   the inversion, to x2_new after the forward multiply).
//! - `y2`: 256 bits holding the quantum point's y. Coordinate ops + the square
//!   touch `y2[..256]`.
//! - `ox`, `oy`: 256-bit CLASSICAL (`BitId`) input registers holding the other
//!   point Q's coordinates (one value per shot -- the fuzzer's runtime control).
//!   Each coordinate step LOADS them into a transient quantum temp (`x_if_bit`),
//!   applies the q-q mod-add/sub, then UNLOADS -- so they are never resident at
//!   the GCD peak (the product-min choice: -512q vs holding them as quantum).
//! The gcd apply's 256-bit scratch is allocated + freed INSIDE
//! `mod_mul_inverse_in_place`, so it is never resident at the GCD peak nor
//! across the square step.

use super::arith::{mod_add, mod_double, mod_double_reverse, mod_neg, mod_sub_vented};
use super::gcd::{mod_mul_inverse_in_place, Direction};
use super::square::mod_square_sub_pm_secp256k1_symmetric;
use super::{B, BExt};
use crate::circuit::{BitId, QubitId};

const N: usize = 256;

/// `dst := dst (+|-) coord (mod q)`, where `coord` is the other point's 256-bit
/// coordinate held in a CLASSICAL `BitId` register (value < q, one per shot). We
/// LOAD it into a transient 256-bit quantum temp (`x_if_bit`), do the
/// unconditional q-q pseudo-Mersenne mod-add/sub (`mod_add`/`mod_sub_vented`),
/// then UNLOAD the temp back to |0>. Keeping ox/oy classical -- loaded only at these
/// (off-peak) coordinate steps, never resident during the GCD -- is the
/// product-min choice (-512q vs holding both as quantum registers). The temp is
/// freed inside the step.
fn coord_addsub(circ: &mut B, dst: &[QubitId], coord: &[BitId], subtract: bool) {
    debug_assert_eq!(dst.len(), N);
    debug_assert_eq!(coord.len(), N);
    let temp = circ.alloc_qubits(N);
    for i in 0..N {
        circ.x_if_bit(temp[i], coord[i]); // load: temp := coord (per-shot classical)
    }
    // UNCONTROLLED vented q-q mod-add/sub (NOT |1>-ctrl normal Cuccaro). dst is
    // modified; `temp` (= coord) is UNTOUCHED, so the unload below is clean on
    // every input.
    if subtract {
        mod_sub_vented(circ, &temp, dst);
    } else {
        mod_add(circ, &temp, dst);
    }
    for i in 0..N {
        circ.x_if_bit(temp[i], coord[i]); // unload: temp := coord XOR coord == 0
    }
    for q in temp {
        circ.zero_and_free(q);
    }
}

/// `dst += 3*coord (mod q)` via `coord + 2*coord`: one mod-add of `coord`, a
/// mod-double of a loaded copy, a mod-add of `2*coord`, then the matched reverse
/// mod-double to restore the copy to `coord` (one fewer mod-add than three). The
/// double needs a 257-bit operand with a |0> overflow slot, so we copy `coord`
/// into a 257-bit temp (top bit |0>) and clear it afterwards.
fn coord_add3x(circ: &mut B, dst: &[QubitId], coord: &[BitId]) {
    debug_assert_eq!(dst.len(), N);
    debug_assert_eq!(coord.len(), N);
    // LOAD coord into a 257-bit temp (low 256 = coord, bit 256 = |0> double
    // overflow slot). 3*coord is done as coord + 2*coord -- two representative
    // q-q mod-adds + a mod-double, NOT a classical `3*ox mod q` precompute.
    let temp: Vec<QubitId> = (0..=N).map(|_| circ.alloc_qubit()).collect();
    for i in 0..N {
        circ.x_if_bit(temp[i], coord[i]); // temp[..N] := coord (per-shot classical)
    }
    // UNCONTROLLED vented; dst modified, temp[..N] restored to coord by
    // mod_double_reverse before unload.
    mod_add(circ, &temp[..N], dst); // dst += coord
    mod_double(circ, &temp); // temp := 2*coord mod q
    mod_add(circ, &temp[..N], dst); // dst += 2*coord
    mod_double_reverse(circ, &temp); // temp := coord (value-inverse of the double)
    for i in 0..N {
        circ.x_if_bit(temp[i], coord[i]); // unload: temp := 0
    }
    for q in temp {
        circ.zero_and_free(q);
    }
}

/// In-place reverse mod-subtract: `x := coord - x (mod q)` over 256 bits, where
/// `coord` is a 256-bit CLASSICAL register holding a generic value < q, via the
/// identity `coord - x = -(x - coord)`:
///   t := coord (load classical into a quantum temp);
///   x := x - coord    (UNCONTROLLED q-q sub; dst = x, t = coord UNTOUCHED);
///   unload t (= coord XOR coord = 0); free t -- clean on EVERY input (t never
///     modified, so NO temp-restore round-trip, unlike a `coord - x` into the temp);
///   x := -x = mod_neg(x) = coord - x.
/// The two mod ops (sub + negate-via-const-add) are representative quantum
/// arithmetic against the quantum `x`; `coord` being classical only changes the
/// load/unload from `CX` to `x_if_bit` (0 Toffoli). Boundary: `mod_neg` lands on
/// `q` only when `x - coord == 0` (i.e. x == coord, a degenerate input),
/// excluded with the other generic-add preconditions.
fn coord_rsub(circ: &mut B, x: &[QubitId], coord: &[BitId]) {
    debug_assert_eq!(x.len(), N);
    debug_assert_eq!(coord.len(), N);
    let t: Vec<QubitId> = (0..N).map(|_| circ.alloc_qubit()).collect();
    for i in 0..N {
        circ.x_if_bit(t[i], coord[i]); // load: t := coord (per-shot classical)
    }
    mod_sub_vented(circ, &t, x); // x := x - coord  (dst = x; t = coord UNTOUCHED)
    for i in 0..N {
        circ.x_if_bit(t[i], coord[i]); // unload: t := coord XOR coord == 0 (clean)
    }
    for q in t {
        circ.zero_and_free(q);
    }
    mod_neg(circ, x); // x := -(x - coord) = coord - x
}

/// Product-min secp256k1 in-place EC point addition: `(x2, y2) -> P + Q`.
///
/// `x2`: 256-bit register, pre = P.x in [0,q), post = (P+Q).x mod q.
/// `y2`: 256-bit register, pre = P.y in [0,q), post = (P+Q).y mod q.
/// `ox`, `oy`: 256-bit CLASSICAL (`BitId`) registers holding Q.x, Q.y in [0,q)
///   (the other point, one value per shot -- loaded into transient quantum
///   temps only at the off-peak coordinate steps, never resident at the GCD peak).
/// The gcd apply's 256-bit scratch is internal to `mod_mul_inverse_in_place`.
///
/// Preconditions: P != Q and P != -Q (generic add, no doubling / identity).
/// `x2` (= the inversion's GCD input `dx = P.x - Q.x`) must be a schedule
/// FITTING input -- a width-truncating `dx` makes the forward GCD's
/// register-shrink `zero_and_free` panic.
pub fn ec_add(
    circ: &mut B,
    x2: &mut Vec<QubitId>,
    y2: &[QubitId],
    ox: &[BitId],
    oy: &[BitId],
) {
    assert_eq!(x2.len(), N, "x2 is 256 bits");
    assert_eq!(y2.len(), N, "y2 is 256 bits");
    assert_eq!(ox.len(), N, "ox is 256 classical bits");
    assert_eq!(oy.len(), N, "oy is 256 classical bits");

    // Step 3/4: x2 -= ox ; y2 -= oy.  => (dx, dy).
    coord_addsub(circ, x2, ox, true);
    coord_addsub(circ, &y2[..N], oy, true);

    // Step 6: y2 *= x2^-1 (gcd inversion). The multiplicand starts in y2[..N]
    // (= dy); after the inverse apply y2 holds lambda = dy * dx^-1 mod q, and x2
    // is restored to dx. `mod_mul_inverse_in_place` takes/returns the 256-bit x
    // register; y2 and the internal tmp scratch are 256-bit.
    let xv = std::mem::take(x2);
    *x2 = mod_mul_inverse_in_place(circ, xv, y2, Direction::Inverse);

    // Step 7: x2 += 3*ox.  => (P.x + 2*Q.x, lambda)  [x2 currently = dx = P.x-Q.x].
    coord_add3x(circ, x2, ox);

    // Step 10: x2 -= lambda^2 mod q.  (lambda = y2[..N]).
    mod_square_sub_pm_secp256k1_symmetric(circ, &y2[..N], x2);

    // Step 11: y2 *= x2 (gcd forward multiply). x2 restored to the post-square
    // value; y2 = lambda * x2 mod q.
    let xv = std::mem::take(x2);
    *x2 = mod_mul_inverse_in_place(circ, xv, y2, Direction::Forward);

    // Step 14: y2 -= oy.   Step 15: x2 := ox - x2.  => (P+Q).x.
    coord_addsub(circ, &y2[..N], oy, true);
    coord_rsub(circ, x2, ox);
}

