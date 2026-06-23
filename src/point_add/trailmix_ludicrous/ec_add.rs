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

use super::arith::{
    mod_add, mod_add_exact, mod_neg, mod_sub_classical_low3, mod_sub_shifted_low, mod_sub_vented,
};
use super::gcd::{mod_mul_inverse_in_place, Direction};
use super::square::mod_square_sub_pm_secp256k1_symmetric;
use super::{B, BExt};
use crate::point_add::{arith::mod_const_minus_reg_qb, SECP256K1_P};
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
    let split_low3 = subtract
        && std::env::var("TLM_COORD_SPLIT_LOW3")
            .ok()
            .as_deref()
            .unwrap_or("0")
            != "0";
    if split_low3 {
        let temp = circ.alloc_qubits(N - 3);
        for i in 3..N {
            circ.x_if_bit(temp[i - 3], coord[i]);
        }
        mod_sub_shifted_low(circ, &temp, dst, 3);
        for i in 3..N {
            circ.x_if_bit(temp[i - 3], coord[i]);
        }
        for q in temp {
            circ.zero_and_free(q);
        }
        mod_sub_classical_low3(circ, dst, &coord[..3]);
        return;
    }
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

/// `dst += 3*coord (mod q)`.
///
/// `coord` (= ox) is CLASSICAL, so `3*coord mod q` is itself a classical value
/// known per shot. We derive it in the classical control system (BitId
/// arithmetic -- ZERO Toffoli, ZERO Clifford in the scorer) and then do a SINGLE
/// generic-classical mod-add of that derived register.
///
/// The previous form did `coord + 2*coord` with two q-q mod-adds plus a
/// mod-double + reverse (~754 Toffoli). Folding the `*3 mod q` into the free
/// classical domain leaves just one mod-add (~326 Toffoli) -- the doubling, the
/// second mod-add, and the 257-bit temp all vanish.
fn coord_add3x(circ: &mut B, dst: &[QubitId], coord: &[BitId]) {
    debug_assert_eq!(dst.len(), N);
    debug_assert_eq!(coord.len(), N);
    // Derive t = 3*coord mod q  (classical, free).
    let three_coord = classical_times3_mod_q(circ, coord);
    // Single generic-classical mod-add: dst += t (mod q). `t < q` by construction.
    // Load t into a transient quantum temp, do an EXACT (full-width) mod-add, unload.
    let temp = circ.alloc_qubits(N);
    for i in 0..N {
        circ.x_if_bit(temp[i], three_coord[i]); // temp := t (per-shot classical)
    }
    mod_add_exact(circ, &temp, dst); // dst += t (mod q), exact -- temp untouched
    for i in 0..N {
        circ.x_if_bit(temp[i], three_coord[i]); // unload: temp := 0
    }
    for q in temp {
        circ.zero_and_free(q);
    }
    // Release the derived classical bits back to |0> (clean: store 0).
    for &b in &three_coord {
        circ.bit_store0(b);
    }
}

/// Compute `t = 3*coord mod q` (q = 2^256 - C, C = 2^32 + 977) entirely in the
/// classical (BitId) domain and return the 256 result bits, value in [0, q).
///
/// EVERY op here is BitStore0/BitStore1/BitInvert (some condition-gated): they
/// touch only the per-shot classical control bits, so the scorer counts ZERO
/// Toffoli and ZERO Clifford for the whole derivation. `coord < q` (a canonical
/// field element) on every shot.
///
/// All arithmetic bodies run UNCONDITIONALLY (every scratch bit is cleared and
/// updated on all shots). Per-shot data choices are encoded by building data
/// registers (e.g. `C*hi`, `q`-or-`0`) via AND-gated copies, never by gating an
/// adder. This keeps scratch bits clean on all shots.
///
/// Math: s = 3*coord < 3*2^256, so s = hi*2^256 + lo with hi in {0,1,2}, lo the
/// low 256 bits. 2^256 ≡ C (mod q), so s ≡ hi*C + lo (mod q). r = lo + hi*C <
/// 2^256 + 2*C < 2q, so a single conditional `-q` lands in [0, q).
fn classical_times3_mod_q(circ: &mut B, coord: &[BitId]) -> Vec<BitId> {
    debug_assert_eq!(coord.len(), N);
    const C: u128 = (1u128 << 32) + 977; // q = 2^256 - C

    // 1) s = 3*coord as a 258-bit classical value (lo = s[0..256], hi = s[256..258]).
    let s: Vec<BitId> = circ.alloc_bits(N + 2);
    for &b in &s {
        circ.bit_store0(b);
    }
    classical_add_into(circ, &s, coord); // s  = coord
    classical_add_into(circ, &s, coord); // s += coord
    classical_add_into(circ, &s, coord); // s += coord  => 3*coord

    // 2) hival = hi*C, where hi = s[256] + 2*s[257] in {0,1,2}. Build the addend
    //    register hival (35 bits is enough: 2*C < 2^34) by AND-gated constant
    //    copies, unconditionally, then ripple-add into the low part.
    //    hi*C = C*(s256) + 2C*(s257).
    let r: Vec<BitId> = circ.alloc_bits(N + 1); // r holds lo + hi*C, < 2^257
    for i in 0..N {
        circ.bit_copy(r[i], s[i]); // r[0..256] = lo
    }
    circ.bit_store0(r[N]); // overflow slot
    // Build addend register av = C*s256 + 2C*s257 (a small classical value), then
    // ripple-add av into r. av fits in 35 bits.
    let av_bits = 35usize;
    let av: Vec<BitId> = circ.alloc_bits(av_bits);
    classical_set_const_times_bit(circ, &av, C, s[N], false); // av  = C   if s256
    classical_add_const_times_bit(circ, &av, 2 * C, s[N + 1]); // av += 2C  if s257
    classical_add_into(circ, &r, &av); // r += av  => r = lo + hi*C

    // 3) Conditional subtract of q. r - q = r - 2^256 + C = (r + C) mod 2^256 with
    //    the carry telling us r >= q. Compute tmp = r + C (258-bit). r >= q iff
    //    tmp has any bit >= 256 set. If so result = low256(tmp); else result = low256(r).
    let tmp: Vec<BitId> = circ.alloc_bits(N + 2);
    for i in 0..(N + 1) {
        circ.bit_copy(tmp[i], r[i]);
    }
    circ.bit_store0(tmp[N + 1]);
    {
        // unconditional add of the constant C into tmp
        let cbits: Vec<BitId> = circ.alloc_bits(av_bits);
        classical_set_const(circ, &cbits, C); // cbits = C
        classical_add_into(circ, &tmp, &cbits);
        for &b in &cbits {
            circ.bit_store0(b);
        }
    }
    // geflag = tmp[256] | tmp[257]   (r >= q)
    let geflag = circ.alloc_bit();
    circ.bit_store0(geflag);
    circ.push_condition(tmp[N]);
    circ.bit_store1(geflag);
    circ.pop_condition();
    circ.push_condition(tmp[N + 1]);
    circ.bit_store1(geflag);
    circ.pop_condition();

    // result[i] = geflag ? tmp[i] : r[i]
    let result: Vec<BitId> = circ.alloc_bits(N);
    for i in 0..N {
        circ.bit_store0(result[i]);
        circ.push_condition(geflag);
        circ.push_condition(tmp[i]);
        circ.bit_store1(result[i]); // geflag & tmp[i]
        circ.pop_condition();
        circ.pop_condition();
        circ.bit_invert(geflag); // !geflag
        circ.push_condition(geflag);
        circ.push_condition(r[i]);
        circ.bit_store1(result[i]); // !geflag & r[i]
        circ.pop_condition();
        circ.pop_condition();
        circ.bit_invert(geflag); // restore
    }

    // Release scratch back to 0 (free, clean).
    circ.bit_store0(geflag);
    for &b in tmp.iter().chain(av.iter()).chain(r.iter()).chain(s.iter()) {
        circ.bit_store0(b);
    }
    result
}

/// Set classical register `dst := k` (compile-time constant, mod 2^dst.len()).
fn classical_set_const(circ: &mut B, dst: &[BitId], k: u128) {
    for (i, &b) in dst.iter().enumerate() {
        let bit = i < 128 && ((k >> i) & 1) == 1;
        if bit {
            circ.bit_store0(b);
            circ.bit_invert(b); // unconditional set to 1
        } else {
            circ.bit_store0(b);
        }
    }
}

/// `dst := k AND gate` per bit, i.e. dst = (gate ? k : 0). If `accumulate` is
/// false this overwrites dst; the helper here always overwrites.
fn classical_set_const_times_bit(circ: &mut B, dst: &[BitId], k: u128, gate: BitId, _accumulate: bool) {
    for (i, &b) in dst.iter().enumerate() {
        circ.bit_store0(b);
        let bit = i < 128 && ((k >> i) & 1) == 1;
        if bit {
            circ.push_condition(gate);
            circ.bit_store1(b); // b = gate
            circ.pop_condition();
        }
    }
}

/// `dst += (k AND gate)` via building the gated constant register and ripple-add.
fn classical_add_const_times_bit(circ: &mut B, dst: &[BitId], k: u128, gate: BitId) {
    let w = dst.len();
    let addend: Vec<BitId> = circ.alloc_bits(w);
    classical_set_const_times_bit(circ, &addend, k, gate, false);
    classical_add_into(circ, dst, &addend);
    for &b in &addend {
        circ.bit_store0(b);
    }
}

/// UNCONDITIONAL classical ripple add `acc += addend (mod 2^acc.len())`, all over
/// BitIds. `addend` may be shorter than `acc` (zero-extended). Free in the scorer.
/// Every scratch bit is cleared and updated on ALL shots (no condition gating),
/// so it is clean regardless of per-shot data.
fn classical_add_into(circ: &mut B, acc: &[BitId], addend: &[BitId]) {
    let carry = circ.alloc_bit();
    circ.bit_store0(carry);
    let newcarry = circ.alloc_bit();
    for i in 0..acc.len() {
        let a_i = addend.get(i).copied();
        // newcarry = majority(acc[i], a_i, carry), computed before overwriting acc[i].
        circ.bit_store0(newcarry);
        if let Some(a) = a_i {
            circ.bit_and_xor_into(newcarry, acc[i], a);
            circ.bit_and_xor_into(newcarry, acc[i], carry);
            circ.bit_and_xor_into(newcarry, a, carry);
        } else {
            circ.bit_and_xor_into(newcarry, acc[i], carry);
        }
        // sum bit -> acc[i]
        if let Some(a) = a_i {
            circ.bit_xor_into(acc[i], a);
        }
        circ.bit_xor_into(acc[i], carry);
        // carry := newcarry
        circ.bit_copy(carry, newcarry);
    }
    circ.bit_store0(newcarry);
    circ.bit_store0(carry);
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
    if std::env::var("TLM_FUSE_X_RESTORE")
        .ok()
        .as_deref()
        == Some("1")
    {
        mod_const_minus_reg_qb(circ, x, coord, SECP256K1_P);
        return;
    }
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
    circ.set_phase("tlm_coord_x_sub");
    coord_addsub(circ, x2, ox, true);
    circ.set_phase("tlm_coord_y_sub");
    coord_addsub(circ, &y2[..N], oy, true);

    // Step 6: y2 *= x2^-1 (gcd inversion). The multiplicand starts in y2[..N]
    // (= dy); after the inverse apply y2 holds lambda = dy * dx^-1 mod q, and x2
    // is restored to dx. `mod_mul_inverse_in_place` takes/returns the 256-bit x
    // register; y2 and the internal tmp scratch are 256-bit.
    circ.set_phase("tlm_inverse");
    let xv = std::mem::take(x2);
    *x2 = mod_mul_inverse_in_place(circ, xv, y2, Direction::Inverse);

    // Step 7: x2 += 3*ox.  => (P.x + 2*Q.x, lambda)  [x2 currently = dx = P.x-Q.x].
    circ.set_phase("tlm_coord_add3x");
    coord_add3x(circ, x2, ox);

    // Step 10: x2 -= lambda^2 mod q.  (lambda = y2[..N]).
    circ.set_phase("tlm_square");
    mod_square_sub_pm_secp256k1_symmetric(circ, &y2[..N], x2);

    // Step 11: y2 *= x2 (gcd forward multiply). x2 restored to the post-square
    // value; y2 = lambda * x2 mod q.
    circ.set_phase("tlm_forward_multiply");
    let xv = std::mem::take(x2);
    *x2 = mod_mul_inverse_in_place(circ, xv, y2, Direction::Forward);

    // Step 14: y2 -= oy.   Step 15: x2 := ox - x2.  => (P+Q).x.
    circ.set_phase("tlm_coord_y_sub_final");
    coord_addsub(circ, &y2[..N], oy, true);
    circ.set_phase("tlm_coord_rsub_final");
    coord_rsub(circ, x2, ox);
}


/// TEST-ONLY: build a tiny circuit that loads classical ox into reg0 (256 bits),
/// computes t = 3*ox mod q classically, and writes t into reg1 (256 classical
/// bits) so a harness can read it back. Returns (ops, ox_reg_bits, t_reg_bits).
pub fn build_times3_test() -> (Vec<crate::circuit::Op>, Vec<BitId>, Vec<BitId>) {
    let mut circ = B::new_for_test();
    let ox = circ.alloc_bits(N);
    let t = classical_times3_mod_q(&mut circ, &ox);
    // copy t into a stable output register tout
    let tout = circ.alloc_bits(N);
    for i in 0..N {
        circ.bit_copy(tout[i], t[i]);
    }
    circ.declare_bit_register(&ox);
    circ.declare_bit_register(&tout);
    (circ.take_ops(), ox, tout)
}

/// TEST-ONLY: build a circuit that loads classical ox into reg-ox, a quantum
/// dst (P.x) into reg-dst, runs coord_add3x (dst += 3*ox mod q), and exposes
/// reg-dst (quantum) + reg-ox so a harness can verify dst' == (dst + 3*ox) mod q.
pub fn build_add3x_test() -> (Vec<crate::circuit::Op>, Vec<BitId>, Vec<QubitId>) {
    let mut circ = B::new_for_test();
    let dst: Vec<QubitId> = circ.alloc_qubits(N);
    let ox = circ.alloc_bits(N);
    coord_add3x(&mut circ, &dst, &ox);
    circ.declare_qubit_register(&dst);
    circ.declare_bit_register(&ox);
    (circ.take_ops(), ox, dst)
}

/// ORIGINAL add3x body (for representative comparison).
fn coord_add3x_orig(circ: &mut B, dst: &[QubitId], coord: &[BitId]) {
    let temp: Vec<QubitId> = (0..=N).map(|_| circ.alloc_qubit()).collect();
    for i in 0..N {
        circ.x_if_bit(temp[i], coord[i]);
    }
    mod_add(circ, &temp[..N], dst);
    super::arith::mod_double(circ, &temp);
    mod_add(circ, &temp[..N], dst);
    super::arith::mod_double_reverse(circ, &temp);
    for i in 0..N {
        circ.x_if_bit(temp[i], coord[i]);
    }
    for q in temp {
        circ.zero_and_free(q);
    }
}

pub fn build_add3x_test_orig() -> (Vec<crate::circuit::Op>, Vec<BitId>, Vec<QubitId>) {
    let mut circ = B::new_for_test();
    let dst: Vec<QubitId> = circ.alloc_qubits(N);
    let ox = circ.alloc_bits(N);
    coord_add3x_orig(&mut circ, &dst, &ox);
    circ.declare_qubit_register(&dst);
    circ.declare_bit_register(&ox);
    (circ.take_ops(), ox, dst)
}
