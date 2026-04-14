//! Reversible secp256k1 point addition circuit.
//!
//! THE editable file for the research loop. Everything else in `src/` is
//! stable harness; all circuit construction lives here.
//!
//! This circuit is specialized to secp256k1. The curve parameters
//!   p = 2^256 - 2^32 - 977
//!   a = 0, b = 7
//! are hard-coded. Specialization lets later optimization passes exploit
//! the Solinas structure of p (sparse low word, mostly-ones upper words)
//! for faster modular reduction. Generalizing is an explicit non-goal.
//!
//! # Interface
//! `build(b)` allocates four 256-wide registers in declaration order —
//! target_x (qubits), target_y (qubits), offset_x (bits), offset_y (bits)
//! — and emits gates that mutate the target registers into (P + Q) where
//! P is the quantum point in targets and Q is the classical point in
//! offsets. The harness validates against `WeierstrassEllipticCurve::add`.
//!
//! # Algorithm
//! Standard affine addition with Roetteler-style two-Kaliski uncomputation:
//!
//!   1. Px -= Qx,  Py -= Qy          (register now holds dx, dy)
//!   2. kaliski_inv_inplace(Px)       (Px ← dx^{-1})
//!   3. lam += Py * Px                (lam ← (dy)(dx^{-1}) = λ)
//!   4. kaliski_inv_inplace(Px)       (Px ← dx)
//!   5. Py -= lam * Px                (Py ← 0)
//!   6. Px -= lam*lam                 (Px ← dx - λ²)
//!   7. Px ← -Px                      (Px ← λ² - dx)
//!   8. Px -= 2*Qx                    (Px ← λ² - Px_orig - Qx = Rx)
//!   9. Py += lam * Qx                (Py ← λ·Qx)
//!  10. Py -= lam * Px                (Py ← λ·Qx - λ·Rx)
//!  11. Py -= Qy                      (Py ← Ry, via the identity
//!                                      Ry = λ(Qx - Rx) - Qy)
//!  12. Uncompute lam via the inverse path using the (Rx, Ry) state.
//!
//! Step 12 in detail (uses the identity λ = (Qy + Ry) / (Qx - Rx)):
//!     a. Px -= Qx; Px ← -Px            (Px ← Qx - Rx)
//!     b. kaliski_inv_inplace(Px)       (Px ← (Qx - Rx)^{-1})
//!     c. lam -= Py * Px                (lam -= Ry / (Qx - Rx))
//!     d. lam -= Qy * Px                (lam -= Qy / (Qx - Rx))
//!                                        → lam = 0
//!     e. kaliski_inv_inplace(Px)       (Px ← Qx - Rx)
//!     f. Px ← -Px; Px += Qx            (Px ← Rx)
//!
//! # Primitive layer
//! All modular arithmetic is built on a single Cuccaro ripple-carry
//! adder operating on `(n+1)`-wide extended registers. Subtract =
//! forward complement + add + back complement. Modular reduction
//! after add/sub is: (cond-sub p) + (cond-add p) controlled by the
//! resulting sign bit.
//!
//! # Current status
//! First-pass baseline: correctness-first, no optimization. Kaliski is
//! implemented as the textbook binary almost-inverse (2n iterations).
//! Expected gate counts far exceed zenodo's targets; the research loop
//! reduces them.

use alloy_primitives::U256;

use crate::builder::{Builder, Layout};
use crate::circuit::{BitId, QubitId};

pub const N: usize = 256;

/// secp256k1 prime:  p = 2^256 - 2^32 - 977.
pub const SECP256K1_P: U256 = U256::from_limbs([
    0xFFFFFFFEFFFFFC2F,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
]);

/// secp256k1 curve coefficient a = 0.
pub const SECP256K1_A: U256 = U256::ZERO;

/// secp256k1 curve coefficient b = 7.
pub const SECP256K1_B: U256 = U256::from_limbs([7, 0, 0, 0]);

// ─── helpers: bit access on U256 ────────────────────────────────────────────

fn bit(c: U256, i: usize) -> bool {
    // alloy's U256::bit returns bool for index < 256.
    c.bit(i)
}

// ═══════════════════════════════════════════════════════════════════════════
//  Cuccaro ripple-carry adder
// ═══════════════════════════════════════════════════════════════════════════
//
// Operates on two n-wide qubit registers `a` (addend, unchanged) and
// `acc` (accumulator, becomes a + acc mod 2^n). Also takes:
//   * c_in: one ancilla qubit, = 0 on entry, = 0 on exit (unchanged)
//   * z   : one ancilla qubit, = 0 on entry, = carry_out ⊕ z_in on exit
//           (i.e., the output carry is XORed into z; pass a fresh 0 bit
//           to receive the high bit)
//
// Based on Cuccaro et al. 2004 (arXiv:quant-ph/0410184), Figure 3.
//
// `MAJ(x, y, w)` triple:
//     CX(w, y)        # y ← y ⊕ w
//     CX(w, x)        # x ← x ⊕ w
//     CCX(x, y, w)    # w ← w ⊕ (x·y)        w becomes MAJ(w_old, y_old, x_old)
//
// `UMA(x, y, w)` triple (undoes MAJ, leaves sum bit in y):
//     CCX(x, y, w)
//     CX(w, x)
//     CX(x, y)

fn maj(b: &mut Builder, x: QubitId, y: QubitId, w: QubitId) {
    b.cx(w, y);
    b.cx(w, x);
    b.ccx(x, y, w);
}

fn uma(b: &mut Builder, x: QubitId, y: QubitId, w: QubitId) {
    b.ccx(x, y, w);
    b.cx(w, x);
    b.cx(x, y);
}

/// In-place addition `acc += a mod 2^n` on quantum n-bit registers.
/// * `c_in` is a fresh ancilla qubit at 0 on entry and returns to 0.
/// * `a` unchanged; `acc` becomes (a + acc) mod 2^n.
/// Pure mod-2^n: the high carry is discarded (no `z` ancilla). This is
/// honestly reversible because the last MAJ/UMA pair cancel out the
/// carry information on `a[n-1]`.
fn cuccaro_add(b: &mut Builder, a: &[QubitId], acc: &[QubitId], c_in: QubitId) {
    let n = a.len();
    assert_eq!(n, acc.len());
    if n == 0 { return; }
    if n == 1 {
        // acc[0] += a[0] + c_in  mod 2 ; c_in → 0
        b.cx(c_in, acc[0]);
        b.cx(a[0], acc[0]);
        return;
    }

    // Forward MAJ sweep.
    maj(b, c_in, acc[0], a[0]);
    for i in 1..n - 1 {
        maj(b, a[i - 1], acc[i], a[i]);
    }

    // Final sum bit: sum[n-1] = acc[n-1] XOR a[n-1] XOR carry_in_to_n-1,
    // where carry_in_to_n-1 is in a[n-2] after the MAJ sweep.
    b.cx(a[n - 2], acc[n - 1]);
    b.cx(a[n - 1], acc[n - 1]);

    // Reverse UMA sweep (skips the final MAJ since we didn't do it).
    for i in (1..n - 1).rev() {
        uma(b, a[i - 1], acc[i], a[i]);
    }
    uma(b, c_in, acc[0], a[0]);
}

/// Reverse of `cuccaro_add`: performs `acc -= a mod 2^n`.
/// Implemented as the exact inverse gate sequence of `cuccaro_add`.
fn cuccaro_sub(b: &mut Builder, a: &[QubitId], acc: &[QubitId], c_in: QubitId) {
    let n = a.len();
    assert_eq!(n, acc.len());
    if n == 0 { return; }
    if n == 1 {
        // Inverse of (cx c_in acc; cx a acc) is the same two gates in reverse.
        b.cx(a[0], acc[0]);
        b.cx(c_in, acc[0]);
        return;
    }

    // Inverse of `uma(c_in, acc[0], a[0])`, then the rest of UMA sweep
    // in reverse order.
    inv_uma(b, c_in, acc[0], a[0]);
    for i in 1..n - 1 {
        inv_uma(b, a[i - 1], acc[i], a[i]);
    }

    // Inverse of the final sum writes (both CX self-inverse; reverse order).
    b.cx(a[n - 1], acc[n - 1]);
    b.cx(a[n - 2], acc[n - 1]);

    // Inverse of the forward MAJ sweep.
    for i in (1..n - 1).rev() {
        inv_maj(b, a[i - 1], acc[i], a[i]);
    }
    inv_maj(b, c_in, acc[0], a[0]);
}

fn inv_maj(b: &mut Builder, x: QubitId, y: QubitId, w: QubitId) {
    // maj = CX(w,y); CX(w,x); CCX(x,y,w)
    // inv = CCX(x,y,w); CX(w,x); CX(w,y)
    b.ccx(x, y, w);
    b.cx(w, x);
    b.cx(w, y);
}

fn inv_uma(b: &mut Builder, x: QubitId, y: QubitId, w: QubitId) {
    // uma = CCX(x,y,w); CX(w,x); CX(x,y)
    // inv = CX(x,y); CX(w,x); CCX(x,y,w)
    b.cx(x, y);
    b.cx(w, x);
    b.ccx(x, y, w);
}

// ═══════════════════════════════════════════════════════════════════════════
//  Loading classical operands into a fresh qubit register
// ═══════════════════════════════════════════════════════════════════════════
//
// Cuccaro needs two qubit registers. To add a classical constant or a
// classical bit register to a quantum register, we allocate a fresh
// qubit register, load the classical value into it, run Cuccaro, then
// unload. The load/unload is not counted against Toffolis.

fn load_const(b: &mut Builder, n: usize, c: U256) -> Vec<QubitId> {
    let qs = b.alloc_qubits(n);
    for i in 0..n {
        if bit(c, i) {
            b.x(qs[i]);
        }
    }
    qs
}

fn unload_const(b: &mut Builder, qs: &[QubitId], c: U256) {
    for i in 0..qs.len() {
        if bit(c, i) {
            b.x(qs[i]);
        }
    }
    b.free_qubits_vec(qs);
}

fn load_bits(b: &mut Builder, bits: &[BitId]) -> Vec<QubitId> {
    let n = bits.len();
    let qs = b.alloc_qubits(n);
    for i in 0..n {
        // qs[i] ← bits[i] via conditional X
        b.x_if(qs[i], bits[i]);
    }
    qs
}

fn unload_bits(b: &mut Builder, qs: &[QubitId], bits: &[BitId]) {
    for i in 0..qs.len() {
        b.x_if(qs[i], bits[i]);
    }
    b.free_qubits_vec(qs);
}

// ═══════════════════════════════════════════════════════════════════════════
//  Extended registers and modular reduction
// ═══════════════════════════════════════════════════════════════════════════
//
// All modular arithmetic operates on "extended" registers of width n+1
// where bit n is an overflow/sign ancilla. The primitive quantum
// registers handed to us (Px, Py) are exactly n=256 wide; the extension
// bit is a transient ancilla allocated for the duration of a mod-op.

/// Build an (n+1)-bit view by attaching a freshly-allocated 0 ancilla.
fn ext_reg(b: &mut Builder, reg: &[QubitId]) -> (Vec<QubitId>, QubitId) {
    let ovf = b.alloc_qubit();
    let mut r = reg.to_vec();
    r.push(ovf);
    (r, ovf)
}

/// Release the overflow ancilla (which must be 0 on exit).
fn unext_reg(b: &mut Builder, ovf: QubitId) {
    b.free_qubit(ovf);
}

/// `acc := (acc + a) mod p`. Both `acc` and `a` are n-bit quantum registers
/// with value in [0, p). Extends both to n+1 bits, adds, subtracts p,
/// conditionally adds p back based on the resulting sign bit.
fn mod_add_qq(b: &mut Builder, acc: &[QubitId], a: &[QubitId], p: U256) {
    let n = acc.len();
    assert_eq!(n, a.len());

    let (acc_ext, acc_ovf) = ext_reg(b, acc);
    let (a_ext, a_ovf) = ext_reg(b, a);

    // (n+1)-bit add.
    add_nbit_qq(b, &a_ext, &acc_ext);
    // acc_ext is in [0, 2p) ⊂ [0, 2^{n+1}). Subtract p.
    sub_nbit_const(b, &acc_ext, p);
    // Sign bit (acc_ovf) = 1 iff we went negative. Conditionally add p.
    cadd_nbit_const(b, &acc_ext, p, acc_ovf);
    // acc_ovf is now 0.
    unext_reg(b, a_ovf);
    unext_reg(b, acc_ovf);
    let _ = (acc_ext, a_ext);
}

fn mod_sub_qq(b: &mut Builder, acc: &[QubitId], a: &[QubitId], p: U256) {
    // acc := (acc - a) mod p
    //      = (acc + (p - a)) mod p
    // Use the "add then conditional add" pattern with the subtraction
    // realized by running cuccaro_add on a_complemented then adjusting.
    //
    // Simpler: mod_add_qq(acc, -a) where -a is computed in a temporary.
    // We compute tmp = (p - a) mod p using mod_neg_into_tmp, then
    // mod_add_qq(acc, tmp), then uncompute tmp.
    let n = acc.len();
    let tmp = b.alloc_qubits(n);
    // tmp ← a (copy via CX)
    for i in 0..n {
        b.cx(a[i], tmp[i]);
    }
    // tmp ← -tmp mod p
    mod_neg_inplace(b, &tmp, p);
    // acc += tmp
    mod_add_qq(b, acc, &tmp, p);
    // Uncompute tmp: undo mod_neg, then undo the copy.
    mod_neg_inplace(b, &tmp, p);
    for i in 0..n {
        b.cx(a[i], tmp[i]);
    }
    b.free_qubits_vec(&tmp);
}

fn mod_add_qc(b: &mut Builder, acc: &[QubitId], c: U256, p: U256) {
    // acc := (acc + c) mod p. c is a compile-time constant.
    let n = acc.len();
    let a = load_const(b, n, c);
    mod_add_qq(b, acc, &a, p);
    unload_const(b, &a, c);
}

fn mod_sub_qc(b: &mut Builder, acc: &[QubitId], c: U256, p: U256) {
    // acc := (acc - c) mod p = acc + (p - c) mod p.
    let n = acc.len();
    let c_neg = (p - (c % p)) % p;
    let a = load_const(b, n, c_neg);
    mod_add_qq(b, acc, &a, p);
    unload_const(b, &a, c_neg);
}

fn mod_add_qb(b: &mut Builder, acc: &[QubitId], bits: &[BitId], p: U256) {
    // acc := (acc + bits) mod p. `bits` is a classical bit register.
    let a = load_bits(b, bits);
    mod_add_qq(b, acc, &a, p);
    unload_bits(b, &a, bits);
}

fn mod_sub_qb(b: &mut Builder, acc: &[QubitId], bits: &[BitId], p: U256) {
    // acc := (acc - bits) mod p
    //     = (acc + (p - bits)) mod p
    // Implement as load → neg → add → uneg → unload.
    let n = acc.len();
    let a = load_bits(b, bits);
    mod_neg_inplace(b, &a, p);
    mod_add_qq(b, acc, &a, p);
    mod_neg_inplace(b, &a, p);
    unload_bits(b, &a, bits);
    let _ = n;
}

/// `v := (p - v) mod p`. Operates on an n-bit register in [0, p).
///
/// Implementation uses the reversible identity:
///     p - v = NOT(v) + (p + 1)         (all arithmetic mod 2^n)
/// which holds because NOT(v) = 2^n - 1 - v, so NOT(v) + p + 1 = 2^n + (p - v).
///
/// For v = 0 the result is p, not 0 (non-canonical but ≡ 0 mod p).
/// EC preconditions (dx, dy nonzero) avoid this case in practice.
fn mod_neg_inplace(b: &mut Builder, v: &[QubitId], p: U256) {
    for &q in v {
        b.x(q);
    }
    add_nbit_const(b, v, p.wrapping_add(U256::from(1)));
}

// ═══════════════════════════════════════════════════════════════════════════
//  Non-modular n-bit primitives
// ═══════════════════════════════════════════════════════════════════════════

/// `acc += a mod 2^n`. Caller must pre-extend both slices if they want the
/// top carry absorbed into the accumulator (i.e. pass n+1-bit slices with
/// top bits 0 to get a full n+1-bit add). The carry-out beyond the slice
/// is discarded via `R` on the `z` ancilla — safe when both inputs fit
/// in n-1 bits (as in our mod-p layer where both < 2p < 2^{n+1}).
fn add_nbit_qq(b: &mut Builder, a: &[QubitId], acc: &[QubitId]) {
    assert_eq!(a.len(), acc.len());
    let c_in = b.alloc_qubit();
    cuccaro_add(b, a, acc, c_in);
    b.free_qubit(c_in);
}

fn sub_nbit_qq(b: &mut Builder, a: &[QubitId], acc: &[QubitId]) {
    assert_eq!(a.len(), acc.len());
    let c_in = b.alloc_qubit();
    cuccaro_sub(b, a, acc, c_in);
    b.free_qubit(c_in);
}

fn add_nbit_const(b: &mut Builder, acc: &[QubitId], c: U256) {
    let n = acc.len();
    let a = load_const(b, n, c);
    add_nbit_qq(b, &a, acc);
    unload_const(b, &a, c);
}

fn sub_nbit_const(b: &mut Builder, acc: &[QubitId], c: U256) {
    let n = acc.len();
    let a = load_const(b, n, c);
    sub_nbit_qq(b, &a, acc);
    unload_const(b, &a, c);
}

fn csub_nbit_const(b: &mut Builder, acc: &[QubitId], c: U256, ctrl: QubitId) {
    // acc -= (ctrl ? c : 0). Mirror of cadd_nbit_const.
    let n = acc.len();
    let a = b.alloc_qubits(n);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, a[i]);
        }
    }
    sub_nbit_qq(b, &a, acc);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, a[i]);
        }
    }
    b.free_qubits_vec(&a);
}

fn cadd_nbit_const(b: &mut Builder, acc: &[QubitId], c: U256, ctrl: QubitId) {
    // Conditional add of constant c, controlled by qubit ctrl.
    // Trick: load c into a qubit register via CX-from-ctrl gates
    // (so the loaded value is (ctrl ? c : 0)), then unconditional add,
    // then unload.
    let n = acc.len();
    let a = b.alloc_qubits(n);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, a[i]);
        }
    }
    add_nbit_qq(b, &a, acc);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, a[i]);
        }
    }
    b.free_qubits_vec(&a);
}


// ═══════════════════════════════════════════════════════════════════════════
//  Modular multiplication
// ═══════════════════════════════════════════════════════════════════════════
//
// Shift-and-add, MSB-to-LSB. `acc += x*y mod p`. Iteration:
//
//     for i from n-1 down to 0:
//         acc := 2*acc mod p
//         if y[i]:  acc := acc + x mod p
//
// For q*q mul, y[i] is a qubit; we implement the conditional add by
// CCX-copying x (gated on y[i]) into a temporary, adding, and
// uncopying. For q*b mul, y[i] is a classical bit and the copy is
// done with CX_if gates.

/// `v := 2*v mod p`. In-place via shift-left (swap cascade) then mod-reduce.
///
/// The swap cascade is a hard logical shift: after it, v[0]=0,
/// v[i] = old_v[i-1] for i∈[1,n), and an ovf ancilla holds old_v[n-1].
/// The (n+1)-bit value is exactly 2*old_v ∈ [0, 2p) ⊂ [0, 2^{n+1}).
/// We then subtract p and conditionally add it back using ovf as the
/// sign bit — identical to the mod_add_qq reduction tail.
fn mod_double_inplace(b: &mut Builder, v: &[QubitId], p: U256) {
    let n = v.len();
    let ovf = b.alloc_qubit();

    // Shift left by 1 via swaps: introduces a 0 into v[0], pushes v[n-1] → ovf.
    b.swap(v[n - 1], ovf);
    for i in (0..n - 1).rev() {
        b.swap(v[i], v[i + 1]);
    }

    // Build the (n+1)-bit view and run the standard sub-p / cond-add-p tail.
    let mut v_ext: Vec<QubitId> = v.to_vec();
    v_ext.push(ovf);
    sub_nbit_const(b, &v_ext, p);
    cadd_nbit_const(b, &v_ext, p, ovf);

    // ovf should now be 0.
    b.free_qubit(ovf);
}

/// Inverse of `mod_double_inplace`: `v := v/2 mod p` (where the "halving"
/// is w.r.t. the group (Z/p)*). Since double+halve round-trip, we just
/// run the `mod_double_inplace` gate sequence backwards.
fn mod_halve_inplace(b: &mut Builder, v: &[QubitId], p: U256) {
    let n = v.len();
    let ovf = b.alloc_qubit();
    let mut v_ext: Vec<QubitId> = v.to_vec();
    v_ext.push(ovf);
    // Inverse of cadd_nbit_const(v_ext, p, ovf):
    csub_nbit_const(b, &v_ext, p, ovf);
    // Inverse of sub_nbit_const(v_ext, p):
    add_nbit_const(b, &v_ext, p);
    // Inverse of the swap cascade (swaps self-inverse, order reversed):
    for i in 0..n - 1 {
        b.swap(v[i], v[i + 1]);
    }
    b.swap(v[n - 1], ovf);
    b.free_qubit(ovf);
}

// ═══════════════════════════════════════════════════════════════════════════
//  Conditional modular add/sub helpers
// ═══════════════════════════════════════════════════════════════════════════
//
// Used by the multipliers. Each variant loads `(ctrl ? a : 0)` into a
// fresh temporary via CCX or CX_if, runs the unconditional mod_add_qq /
// mod_sub_qq, then unloads.

fn cmod_add_qq(b: &mut Builder, acc: &[QubitId], a: &[QubitId], ctrl: QubitId, p: U256) {
    let n = acc.len();
    let f = b.alloc_qubits(n);
    for i in 0..n {
        b.ccx(ctrl, a[i], f[i]);
    }
    mod_add_qq(b, acc, &f, p);
    for i in 0..n {
        b.ccx(ctrl, a[i], f[i]);
    }
    b.free_qubits_vec(&f);
}

fn cmod_sub_qq(b: &mut Builder, acc: &[QubitId], a: &[QubitId], ctrl: QubitId, p: U256) {
    let n = acc.len();
    let f = b.alloc_qubits(n);
    for i in 0..n {
        b.ccx(ctrl, a[i], f[i]);
    }
    mod_sub_qq(b, acc, &f, p);
    for i in 0..n {
        b.ccx(ctrl, a[i], f[i]);
    }
    b.free_qubits_vec(&f);
}

fn cmod_add_qq_bit(b: &mut Builder, acc: &[QubitId], a: &[QubitId], ctrl: BitId, p: U256) {
    let n = acc.len();
    let f = b.alloc_qubits(n);
    for i in 0..n {
        b.cx_if(a[i], f[i], ctrl);
    }
    mod_add_qq(b, acc, &f, p);
    for i in 0..n {
        b.cx_if(a[i], f[i], ctrl);
    }
    b.free_qubits_vec(&f);
}

fn cmod_sub_qq_bit(b: &mut Builder, acc: &[QubitId], a: &[QubitId], ctrl: BitId, p: U256) {
    let n = acc.len();
    let f = b.alloc_qubits(n);
    for i in 0..n {
        b.cx_if(a[i], f[i], ctrl);
    }
    mod_sub_qq(b, acc, &f, p);
    for i in 0..n {
        b.cx_if(a[i], f[i], ctrl);
    }
    b.free_qubits_vec(&f);
}

fn mod_mul_add_qq(
    b: &mut Builder,
    acc: &[QubitId],
    x: &[QubitId],
    y: &[QubitId],
    p: U256,
) {
    let n = acc.len();
    let tmp = b.alloc_qubits(n);
    for i in 0..n { b.cx(x[i], tmp[i]); }
    for i in 0..n {
        cmod_add_qq(b, acc, &tmp, y[i], p);
        if i < n - 1 { mod_double_inplace(b, &tmp, p); }
    }
    for _ in 0..(n - 1) { mod_halve_inplace(b, &tmp, p); }
    for i in 0..n { b.cx(x[i], tmp[i]); }
    b.free_qubits_vec(&tmp);
}

fn mod_mul_sub_qq(
    b: &mut Builder,
    acc: &[QubitId],
    x: &[QubitId],
    y: &[QubitId],
    p: U256,
) {
    let n = acc.len();
    let tmp = b.alloc_qubits(n);
    for i in 0..n { b.cx(x[i], tmp[i]); }
    for i in 0..n {
        cmod_sub_qq(b, acc, &tmp, y[i], p);
        if i < n - 1 { mod_double_inplace(b, &tmp, p); }
    }
    for _ in 0..(n - 1) { mod_halve_inplace(b, &tmp, p); }
    for i in 0..n { b.cx(x[i], tmp[i]); }
    b.free_qubits_vec(&tmp);
}

fn mod_mul_add_qb(
    b: &mut Builder,
    acc: &[QubitId],
    x: &[QubitId],
    y: &[BitId],
    p: U256,
) {
    let n = acc.len();
    let tmp = b.alloc_qubits(n);
    for i in 0..n { b.cx(x[i], tmp[i]); }
    for i in 0..n {
        cmod_add_qq_bit(b, acc, &tmp, y[i], p);
        if i < n - 1 { mod_double_inplace(b, &tmp, p); }
    }
    for _ in 0..(n - 1) { mod_halve_inplace(b, &tmp, p); }
    for i in 0..n { b.cx(x[i], tmp[i]); }
    b.free_qubits_vec(&tmp);
}

fn mod_mul_sub_qb(
    b: &mut Builder,
    acc: &[QubitId],
    x: &[QubitId],
    y: &[BitId],
    p: U256,
) {
    let n = acc.len();
    let tmp = b.alloc_qubits(n);
    for i in 0..n { b.cx(x[i], tmp[i]); }
    for i in 0..n {
        cmod_sub_qq_bit(b, acc, &tmp, y[i], p);
        if i < n - 1 { mod_double_inplace(b, &tmp, p); }
    }
    for _ in 0..(n - 1) { mod_halve_inplace(b, &tmp, p); }
    for i in 0..n { b.cx(x[i], tmp[i]); }
    b.free_qubits_vec(&tmp);
}

// ═══════════════════════════════════════════════════════════════════════════
//  Kaliski almost-inverse
// ═══════════════════════════════════════════════════════════════════════════

/// Custom modular halving:  v ← v/2 mod p.
/// If v odd: v := (v + p)/2; else v := v/2.
fn mod_halve_v2(b: &mut Builder, v: &[QubitId], p: U256) {
    let n = v.len();
    let odd = b.alloc_qubit();
    b.cx(v[0], odd);
    let ovf = b.alloc_qubit();
    let mut v_ext: Vec<QubitId> = v.to_vec(); v_ext.push(ovf);
    cadd_nbit_const(b, &v_ext, p, odd);
    // Now v_ext is even and in [0, 2p) ⊂ [0, 2^{n+1}).
    // Shift right: drop bit 0, fill new top bit from ovf.
    for i in 0..n - 1 {
        b.swap(v[i], v[i + 1]);
    }
    b.swap(v[n - 1], ovf);
    b.free_qubit(ovf);
    b.free_qubit(odd);
}

fn cmod_double_inplace(b: &mut Builder, v: &[QubitId], p: U256, ctrl: QubitId) {
    let n = v.len();
    let tmp = b.alloc_qubits(n);
    for i in 0..n { b.ccx(ctrl, v[i], tmp[i]); }
    mod_add_qq(b, v, &tmp, p);
    b.free_qubits_vec(&tmp);
}

fn cmod_halve_inplace(b: &mut Builder, v: &[QubitId], p: U256, ctrl: QubitId) {
    mod_halve_v2(b, v, p);
    b.x(ctrl);
    cmod_double_inplace(b, v, p, ctrl);
    b.x(ctrl);
}

/// flag ^= (u < v).  u and v are n-wide qubit registers, both holding values
/// in [0, 2^{n-1}) so a difference fits in n bits with sign in the top bit.
/// We extend by one ancilla for safety.
fn cmp_lt_into(b: &mut Builder, u: &[QubitId], v: &[QubitId], flag: QubitId) {
    let n = u.len();
    assert_eq!(n, v.len());
    let u_top = b.alloc_qubit();
    let v_top = b.alloc_qubit();
    let mut ue: Vec<QubitId> = u.to_vec(); ue.push(u_top);
    let mut ve: Vec<QubitId> = v.to_vec(); ve.push(v_top);
    sub_nbit_qq(b, &ve, &ue);
    b.cx(u_top, flag);
    add_nbit_qq(b, &ve, &ue);
    b.free_qubit(v_top);
    b.free_qubit(u_top);
    let _ = (ue, ve);
}

/// flag ^= (v != 0). Computes OR of all bits of v into a scratch ancilla,
/// CXs into flag, then properly uncomputes the scratch.
///
/// We use the simple chain: `or[0] = v[0]`, `or[i] = or[i-1] OR v[i]`.
/// OR via de Morgan: `or[i] = NOT((NOT or[i-1]) AND (NOT v[i]))`, i.e.
///   x(or[i-1]); x(v[i]); ccx(or[i-1], v[i], or[i]); x(or[i]);
///   x(v[i]); x(or[i-1]);
/// Each `or[i]` is a fresh ancilla. We compute the chain, CX `or[n-1]`
/// into `flag`, then reverse the chain to return every ancilla to |0⟩.
fn cmp_neq_zero_into(b: &mut Builder, v: &[QubitId], flag: QubitId) {
    let n = v.len();
    assert!(n > 0);
    if n == 1 {
        b.cx(v[0], flag);
        return;
    }

    let or_chain: Vec<QubitId> = b.alloc_qubits(n - 1);
    // or_chain[0] = v[0] OR v[1]
    or_step(b, v[0], v[1], or_chain[0]);
    for i in 1..n - 1 {
        or_step(b, or_chain[i - 1], v[i + 1], or_chain[i]);
    }

    // flag ^= or_chain[n-2]
    b.cx(or_chain[n - 2], flag);

    // Uncompute.
    for i in (1..n - 1).rev() {
        or_step(b, or_chain[i - 1], v[i + 1], or_chain[i]);
    }
    or_step(b, v[0], v[1], or_chain[0]);

    b.free_qubits_vec(&or_chain);
}

/// out ^= (x OR y). `out` starts 0. Uses the de-Morgan form:
///   x(x); x(y); ccx(x, y, out); x(out); x(y); x(x);
/// After this, out = x OR y (assuming out started at 0). Its inverse is
/// the same gate sequence run in reverse — since it's symmetric (all gates
/// involutions, palindromic structure), running the exact same helper
/// again uncomputes it.
fn or_step(b: &mut Builder, x: QubitId, y: QubitId, out: QubitId) {
    b.x(x);
    b.x(y);
    b.ccx(x, y, out);
    b.x(out);
    b.x(y);
    b.x(x);
}

fn kaliski_inv_inplace(b: &mut Builder, v_in: &[QubitId], p: U256) {
    let n = v_in.len();
    let m = n;

    let u: Vec<QubitId> = b.alloc_qubits(m);
    let v_w: Vec<QubitId> = b.alloc_qubits(m);
    let r: Vec<QubitId> = b.alloc_qubits(m);
    let s: Vec<QubitId> = b.alloc_qubits(m);

    // Init u = p (only low n bits significant; bit n is 0 since p < 2^n).
    for i in 0..n {
        if bit(p, i) { b.x(u[i]); }
    }
    // Copy v_in into v_w (low n bits).
    for i in 0..n { b.cx(v_in[i], v_w[i]); }
    // Clear v_in by XOR-ing v_w (which equals v_in here).
    for i in 0..n { b.cx(v_w[i], v_in[i]); }
    // Init s = 1.
    b.x(s[0]);
    // r = 0 already.

    let max_iter = 2 * n + 2;
    for _ in 0..max_iter {
        // active := (v_w != 0).  Compute fresh each iteration.
        let active = b.alloc_qubit();
        cmp_neq_zero_into(b, &v_w, active);

        // u_even = NOT u[0]; v_even = NOT v_w[0].
        let u_even = b.alloc_qubit();
        b.cx(u[0], u_even);
        b.x(u_even);
        let v_even = b.alloc_qubit();
        b.cx(v_w[0], v_even);
        b.x(v_even);

        // Branch flags (mutually exclusive):
        //   A: active & v_even                          → halve v_w; s *= 2
        //   B: active & !v_even & u_even                → halve u; r *= 2
        //   C: active & !v_even & !u_even & (u > v_w)   → u-=v_w; halve u; r += s; s *= 2 ... NO
        //
        // Standard binary EGCD updates when both odd:
        //   if u >= v: u := u - v; r := r - s;  (now u even) then halve u, halve r
        //   if u <  v: v := v - u; s := s - r;  (now v even) then halve v, halve s
        //
        // We'll do the branches in two phases: first the subtract+swap-style work,
        // then a unified halving.

        // gt = (v_w < u) i.e. u > v_w → use cmp_lt_into(v_w, u, gt)
        let gt = b.alloc_qubit();
        cmp_lt_into(b, &v_w, &u, gt);
        // Actually we only need this when both odd. Keep gt always; we gate it later.

        // both_odd = active & !u_even & !v_even
        let both_odd = b.alloc_qubit();
        // both_odd := active AND NOT u_even AND NOT v_even
        // Compute via temporary: t1 = active AND NOT u_even
        let t1 = b.alloc_qubit();
        b.x(u_even);
        b.ccx(active, u_even, t1);
        b.x(u_even);
        b.x(v_even);
        b.ccx(t1, v_even, both_odd);
        b.x(v_even);
        b.free_qubit(t1);

        // case_C = both_odd & gt   (u > v_w)
        let case_c = b.alloc_qubit();
        b.ccx(both_odd, gt, case_c);
        // case_D = both_odd & NOT gt
        let case_d = b.alloc_qubit();
        b.x(gt);
        b.ccx(both_odd, gt, case_d);
        b.x(gt);

        // case_A = active & v_even
        let case_a = b.alloc_qubit();
        b.ccx(active, v_even, case_a);
        // case_B = active & !v_even & u_even
        let case_b = b.alloc_qubit();
        let t2 = b.alloc_qubit();
        b.x(v_even);
        b.ccx(active, v_even, t2);
        b.x(v_even);
        b.ccx(t2, u_even, case_b);
        b.free_qubit(t2);

        // ─ Apply C: u -= v_w; r -= s; (both controlled by case_c) ─
        cmod_sub_qq(b, &u, &v_w, case_c, p);
        cmod_sub_qq(b, &r, &s, case_c, p);
        // ─ Apply D: v_w -= u; s -= r; (controlled by case_d) ─
        cmod_sub_qq(b, &v_w, &u, case_d, p);
        cmod_sub_qq(b, &s, &r, case_d, p);

        // Now under case_C, u is even; under case_D, v_w is even.
        // Combined "u becomes halvable" = case_B OR case_C
        // Combined "v_w becomes halvable" = case_A OR case_D
        let halve_u = b.alloc_qubit();
        b.cx(case_b, halve_u);
        b.cx(case_c, halve_u);
        let halve_v = b.alloc_qubit();
        b.cx(case_a, halve_v);
        b.cx(case_d, halve_v);
        // For r/s halving:
        // After case A: s should be doubled (not halved). Actually, in the
        // standard binary EGCD with INVARIANT r*v_orig ≡ u, s*v_orig ≡ v_w,
        // when v_w is halved we need s to be halved too (so s*v_orig = v_w/2 → s/=2 mod p).
        // Same: when u halved, r halved.
        // So: r halves alongside u, s halves alongside v_w.
        let halve_r = halve_u;  // alias
        let halve_s = halve_v;  // alias

        cmod_halve_inplace(b, &u, p, halve_u);
        cmod_halve_inplace(b, &v_w, p, halve_v);
        cmod_halve_inplace(b, &r, p, halve_r);
        cmod_halve_inplace(b, &s, p, halve_s);

        b.free_qubit(halve_v);
        b.free_qubit(halve_u);
        b.free_qubit(case_b);
        b.free_qubit(case_a);
        b.free_qubit(case_d);
        b.free_qubit(case_c);
        b.free_qubit(both_odd);
        b.free_qubit(gt);
        b.free_qubit(v_even);
        b.free_qubit(u_even);
        b.free_qubit(active);
    }

    // After the loop, u should be 1 and r ≡ v_orig^{-1} (mod p) (or its negation;
    // the exact sign depends on convention). We have r*v_orig ≡ u ≡ 1 (mod p),
    // so r is the inverse.  But r may be in [0, p) — that's fine.

    // Reduce r mod p just in case (should already be in [0, p)).
    // Copy r into v_in.
    for i in 0..n { b.cx(r[i], v_in[i]); }

    b.free_qubits_vec(&s);
    b.free_qubits_vec(&r);
    b.free_qubits_vec(&v_w);
    b.free_qubits_vec(&u);
}

// ═══════════════════════════════════════════════════════════════════════════
//  Top-level point addition
// ═══════════════════════════════════════════════════════════════════════════

pub fn build(b: &mut Builder) -> Layout {
    // Register 0: target_x (quantum)
    let tx = b.alloc_qubits(N);
    let target_x = b.declare_qubit_register(&tx);
    // Register 1: target_y (quantum)
    let ty = b.alloc_qubits(N);
    let target_y = b.declare_qubit_register(&ty);
    // Register 2: offset_x (classical bits)
    let ox = b.alloc_bits(N);
    let offset_x = b.declare_bit_register(&ox);
    // Register 3: offset_y (classical bits)
    let oy = b.alloc_bits(N);
    let offset_y = b.declare_bit_register(&oy);

    // === Point add ===
    //
    // NOTE: the subroutines `mod_mul_*` and `kaliski_inv_inplace` are
    // currently stubbed with `unimplemented!`. Calling `build` will
    // panic at circuit-construction time until those are filled in.
    // This scaffold compiles and exercises the Cuccaro adder layer +
    // the register declarations so the harness interface is validated.

    let p = SECP256K1_P;

    // Step 1-2: Px -= Qx, Py -= Qy
    mod_sub_qb(b, &tx, &ox, p);
    mod_sub_qb(b, &ty, &oy, p);


    let lam = b.alloc_qubits(N);

    kaliski_inv_inplace(b, &tx, p);              // Px ← dx^{-1}
    mod_mul_add_qq(b, &lam, &ty, &tx, p);        // lam += dy · dx^{-1} = λ
    kaliski_inv_inplace(b, &tx, p);              // Px ← dx
    mod_mul_sub_qq(b, &ty, &lam, &tx, p);        // Py -= λ·dx = 0

    // Px := λ² - Px_orig - Qx
    mod_mul_sub_qq(b, &tx, &lam, &lam, p);       // Px ← dx - λ²
    mod_neg_inplace(b, &tx, p);                  // Px ← λ² - dx
    mod_sub_qb(b, &tx, &ox, p);
    mod_sub_qb(b, &tx, &ox, p);                  // Px ← Rx

    // Py := λ·Qx − λ·Rx − Qy
    mod_mul_add_qb(b, &ty, &lam, &ox, p);
    mod_mul_sub_qq(b, &ty, &lam, &tx, p);
    mod_sub_qb(b, &ty, &oy, p);

    // Uncompute lam using λ = (Qy + Ry) / (Qx - Rx).
    mod_sub_qb(b, &tx, &ox, p);
    mod_neg_inplace(b, &tx, p);
    kaliski_inv_inplace(b, &tx, p);
    mod_mul_sub_qq(b, &lam, &ty, &tx, p);
    mod_mul_sub_qb(b, &lam, &tx, &oy, p);
    kaliski_inv_inplace(b, &tx, p);
    mod_neg_inplace(b, &tx, p);
    mod_add_qb(b, &tx, &ox, p);

    b.free_qubits_vec(&lam);

    Layout { target_x, target_y, offset_x, offset_y }
}

