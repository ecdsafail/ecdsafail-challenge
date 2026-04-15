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
use crate::circuit::{BitId, OperationType, QubitId};

// ═══════════════════════════════════════════════════════════════════════════
//  emit_inverse: run a closure, pop the ops it emitted, and re-emit them
//  reversed.
//
//  The closure may contain `alloc_qubit` / `assert_zero_and_free` calls;
//  the R ops that `assert_zero_and_free` produces are SKIPPED during
//  reverse replay. This relies on the forward being "clean" — i.e. each
//  free lands on a qubit that the forward gates already drove to |0⟩
//  before the R. Under that invariant, the reverse gate sequence brings
//  the same qubit back to |0⟩ at the "alloc" point (pre-forward-allocation),
//  and the R we skipped is unnecessary.
//
//  The forward's internal alloc/free bookkeeping in the Builder's free
//  pool is NOT undone by the reverse — the pool state at reverse exit
//  equals the pool state at forward exit. Subsequent allocations in the
//  parent scope reuse those qubit IDs, seeing them at |0⟩ (as zeroed by
//  the reverse gate sequence).
// ═══════════════════════════════════════════════════════════════════════════
fn emit_inverse<F: FnOnce(&mut Builder)>(b: &mut Builder, f: F) {
    let start = b.ops.len();
    f(b);
    let end = b.ops.len();
    // Extract the forward slice and drop it from the builder.
    let fwd: Vec<_> = b.ops[start..end].to_vec();
    b.ops.truncate(start);
    for op in fwd.into_iter().rev() {
        match op.kind {
            OperationType::X
            | OperationType::Z
            | OperationType::CX
            | OperationType::CZ
            | OperationType::CCX
            | OperationType::CCZ
            | OperationType::Swap => b.ops.push(op),
            // R ops are the free markers. They're not directly reversible
            // as gates, but in a clean forward they're preceded by gates
            // that already zero the qubit. We skip them in reverse.
            OperationType::R => {}
            // Metadata ops (register declarations, debug prints) don't
            // affect state and shouldn't appear inside an emit_inverse
            // closure anyway, but skip them if they do.
            OperationType::Register
            | OperationType::AppendToRegister
            | OperationType::DebugPrint => {}
            _ => panic!(
                "emit_inverse: non-invertible op kind {:?} inside forward block",
                op.kind
            ),
        }
    }
}

/// Runs `compute`, then `body`, then the inverse of `compute` — the
/// "with conjugate" pattern from qrisp. `compute` must emit only
/// reversible gates (no alloc/free/R).
fn conjugate<F, G>(b: &mut Builder, compute: F, body: G)
where
    F: Fn(&mut Builder),
    G: FnOnce(&mut Builder),
{
    compute(b);
    body(b);
    emit_inverse(b, compute);
}

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
    b.assert_zero_and_free_vec(qs);
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
    b.assert_zero_and_free_vec(qs);
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
    b.assert_zero_and_free(ovf);
}

/// `acc := (acc + a) mod p`. Both `acc` and `a` are n-bit quantum registers
/// with value in [0, p). Solinas reduction using c = 2^n - p: sum ∈ [0, 2p),
/// then add c, branch on top bit to either clear it (reduction) or undo
/// the add (no reduction). Saves one full (n+1)-wide Cuccaro compared to
/// the sub-p/add-p/csub-p pattern.
fn mod_add_qq(b: &mut Builder, acc: &[QubitId], a: &[QubitId], p: U256) {
    let n = acc.len();
    assert_eq!(n, a.len());
    debug_assert_eq!(n, 256);

    let (acc_ext, acc_ovf) = ext_reg(b, acc);
    let (a_ext, a_ovf) = ext_reg(b, a);

    // Step 1: (n+1)-bit add. acc_ext ∈ [0, 2p).
    add_nbit_qq(b, &a_ext, &acc_ext);

    // Step 2: add c. If sum was >= p, the top bit of (sum + c) becomes 1.
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));
    add_nbit_const(b, &acc_ext, c);

    // Step 3: flag := acc_ovf (= top bit of sum + c).
    let flag = b.alloc_qubit();
    b.cx(acc_ovf, flag);

    // Step 4: if flag=0 (no reduction needed), undo the add of c.
    b.x(flag);
    csub_nbit_const(b, &acc_ext, c, flag);
    b.x(flag);

    // Step 5: if flag=1, clear the top bit (drops 2^n → yields sum - p).
    b.cx(flag, acc_ovf);

    // Step 6: uncompute flag. Same identity as the old version:
    //   flag == (acc_final < a_orig)
    // because in the flag=1 case acc_final = acc_orig + a - p < a (since acc_orig < p),
    // and in the flag=0 case acc_final = acc_orig + a ≥ a.
    cmp_lt_into(b, &acc_ext[..n], &a_ext[..n], flag);
    b.assert_zero_and_free(flag);

    unext_reg(b, a_ovf);
    unext_reg(b, acc_ovf);
    let _ = (acc_ext, a_ext);
}

fn mod_sub_qq(b: &mut Builder, acc: &[QubitId], a: &[QubitId], p: U256) {
    // mod_add_qq is a bijection on (acc, a): (acc, a) ↦ (acc + a mod p, a).
    // Its gate-level inverse therefore acts as (acc, a) ↦ (acc - a mod p, a),
    // which is exactly what we want. emit_inverse replays the forward's gates
    // reversed, skipping R markers — valid because mod_add_qq is clean
    // (every ancilla is driven to |0⟩ before its R).
    let a_copy: Vec<QubitId> = a.to_vec();
    emit_inverse(b, move |b| mod_add_qq(b, acc, &a_copy, p));
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
    // Gate-inverse of mod_add_qb, by the same bijection argument as mod_sub_qq.
    let bits_copy: Vec<BitId> = bits.to_vec();
    emit_inverse(b, move |b| mod_add_qb(b, acc, &bits_copy, p));
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
    b.assert_zero_and_free(c_in);
}

fn sub_nbit_qq(b: &mut Builder, a: &[QubitId], acc: &[QubitId]) {
    assert_eq!(a.len(), acc.len());
    let c_in = b.alloc_qubit();
    cuccaro_sub(b, a, acc, c_in);
    b.assert_zero_and_free(c_in);
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
    b.assert_zero_and_free_vec(&a);
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
    b.assert_zero_and_free_vec(&a);
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

/// `v := 2*v mod p`. In-place via shift-left (swap cascade) + Solinas-style
/// mod reduction. For secp256k1, p = 2^n - c with c = 2^32 + 977, so
/// `T - p = T + c - 2^n`. The reduction becomes: add c, branch on the top
/// bit of the (n+1)-wide shifted register — if set, clear it; else undo
/// the add. Costs two full (n+1)-wide Cuccaro adds instead of three.
fn mod_double_inplace(b: &mut Builder, v: &[QubitId], p: U256) {
    let n = v.len();
    let ovf = b.alloc_qubit();

    // Shift left by 1 via swaps: introduces a 0 into v[0], pushes v[n-1] → ovf.
    b.swap(v[n - 1], ovf);
    for i in (0..n - 1).rev() {
        b.swap(v[i], v[i + 1]);
    }

    let mut v_ext: Vec<QubitId> = v.to_vec();
    v_ext.push(ovf);

    // c = 2^n - p (= 2^32 + 977 for secp256k1). Assumes n == 256 so that
    // 2^n wraps cleanly in U256::MAX + 1 arithmetic.
    debug_assert_eq!(n, 256);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));

    // S := T + c. Fits in n+1 bits.
    add_nbit_const(b, &v_ext, c);

    // flag := (S >= 2^n) = S[n]. S[n]==1 iff we need the reduction.
    let flag = b.alloc_qubit();
    b.cx(ovf, flag);

    // If flag=0, undo the add (we didn't need to reduce).
    b.x(flag);
    csub_nbit_const(b, &v_ext, c, flag);
    b.x(flag);

    // If flag=1, clear the top bit (drops the 2^n from S, giving T - p).
    b.cx(flag, ovf);

    // Uncompute flag via parity: flag == v[0] after the operation.
    // Case flag=0: v = T = 2*v_orig (even) → v[0]=0.
    // Case flag=1: v = T - p. T even, p odd → v is odd → v[0]=1.
    b.cx(v[0], flag);
    b.assert_zero_and_free(flag);
    b.assert_zero_and_free(ovf);
}

/// `v := v/2 mod p`. Gate-inverse of `mod_double_inplace`.
fn mod_halve_inplace(b: &mut Builder, v: &[QubitId], p: U256) {
    let v_copy: Vec<QubitId> = v.to_vec();
    emit_inverse(b, move |b| mod_double_inplace(b, &v_copy, p));
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
    b.assert_zero_and_free_vec(&f);
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
    b.assert_zero_and_free_vec(&f);
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
    b.assert_zero_and_free_vec(&f);
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
    b.assert_zero_and_free_vec(&f);
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
    b.assert_zero_and_free_vec(&tmp);
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
    b.assert_zero_and_free_vec(&tmp);
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
        // Mask the whole conditional-add body by y[i]: on shots where
        // y[i]=0 nothing needs to happen AND nothing should be counted.
        b.push_condition(y[i]);
        cmod_add_qq_bit(b, acc, &tmp, y[i], p);
        b.pop_condition();
        if i < n - 1 { mod_double_inplace(b, &tmp, p); }
    }
    for _ in 0..(n - 1) { mod_halve_inplace(b, &tmp, p); }
    for i in 0..n { b.cx(x[i], tmp[i]); }
    b.assert_zero_and_free_vec(&tmp);
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
        b.push_condition(y[i]);
        cmod_sub_qq_bit(b, acc, &tmp, y[i], p);
        b.pop_condition();
        if i < n - 1 { mod_double_inplace(b, &tmp, p); }
    }
    for _ in 0..(n - 1) { mod_halve_inplace(b, &tmp, p); }
    for i in 0..n { b.cx(x[i], tmp[i]); }
    b.assert_zero_and_free_vec(&tmp);
}

// ═══════════════════════════════════════════════════════════════════════════
//  Kaliski almost-inverse
// ═══════════════════════════════════════════════════════════════════════════

/// Fredkin (controlled swap): swap (a, t) if ctrl. Decomposed as CX/CCX/CX.
fn cswap(b: &mut Builder, ctrl: QubitId, a: QubitId, t: QubitId) {
    b.cx(t, a);
    b.ccx(ctrl, a, t);
    b.cx(t, a);
}

fn cmod_double_inplace(b: &mut Builder, v: &[QubitId], p: U256, ctrl: QubitId) {
    let n = v.len();
    let ovf = b.alloc_qubit();
    let mut v_ext: Vec<QubitId> = v.to_vec();
    v_ext.push(ovf);

    // Conditional left-shift: if ctrl=1, v[n-1] → ovf; v[i] → v[i+1].
    cswap(b, ctrl, v[n - 1], ovf);
    for i in (0..n - 1).rev() {
        cswap(b, ctrl, v[i], v[i + 1]);
    }

    csub_nbit_const(b, &v_ext, p, ctrl);
    cadd_nbit_const(b, &v_ext, p, ovf);
    // ovf ends at 0 by the same argument as mod_double_inplace.
    b.assert_zero_and_free(ovf);
}

/// `cmod_halve_inplace` = exact inverse of `cmod_double_inplace`.
fn cmod_halve_inplace(b: &mut Builder, v: &[QubitId], p: U256, ctrl: QubitId) {
    let n = v.len();
    let ovf = b.alloc_qubit();
    let mut v_ext: Vec<QubitId> = v.to_vec();
    v_ext.push(ovf);

    // Inverse of: cadd(v_ext, p, ovf).
    csub_nbit_const(b, &v_ext, p, ovf);
    // Inverse of: csub(v_ext, p, ctrl).
    cadd_nbit_const(b, &v_ext, p, ctrl);
    // Inverse of cswap cascade (self-inverse; reversed order).
    for i in 0..n - 1 {
        cswap(b, ctrl, v[i], v[i + 1]);
    }
    cswap(b, ctrl, v[n - 1], ovf);

    b.assert_zero_and_free(ovf);
}

/// flag ^= (u < v).  Non-destructive on u and v.
///
/// Uses a MAJ-only carry chain instead of the full sub+add pattern.
/// Identity: u < v iff carry-out of (~u + v) = 1, since
///   ~u + v = (2^n - 1 - u) + v = (v - u) + (2^n - 1)
/// which overflows 2^n iff v - u ≥ 1 iff v > u. We negate u in place,
/// run a forward MAJ sweep over (~u, v, c_in=0), capture u[n-1] (which
/// holds the high carry after the chain), then run the inverse MAJ
/// sweep + un-negate to restore u and v. Cost ≈ 2n CCX, half of the
/// previous sub+add (≈ 4n CCX).
fn cmp_lt_into(b: &mut Builder, u: &[QubitId], v: &[QubitId], flag: QubitId) {
    let n = u.len();
    assert_eq!(n, v.len());

    let c_in = b.alloc_qubit();

    // ~u in place (X is free in the metric).
    for i in 0..n { b.x(u[i]); }

    // Forward MAJ sweep — n MAJs (one more than cuccaro_add, which omits
    // the top one because it doesn't need the carry-out).
    maj(b, c_in, v[0], u[0]);
    for i in 1..n {
        maj(b, u[i - 1], v[i], u[i]);
    }
    // u[n-1] now holds the high carry = (u < v).
    b.cx(u[n - 1], flag);

    // Inverse sweep restores u and v to their (negated u) state.
    for i in (1..n).rev() {
        inv_maj(b, u[i - 1], v[i], u[i]);
    }
    inv_maj(b, c_in, v[0], u[0]);

    // Un-negate u.
    for i in 0..n { b.x(u[i]); }

    b.assert_zero_and_free(c_in);
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

    b.assert_zero_and_free_vec(&or_chain);
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

// ═══════════════════════════════════════════════════════════════════════════
//  Primitives for the Kaliski port (qrisp-style)
// ═══════════════════════════════════════════════════════════════════════════

/// 2-controlled X with per-control polarity. `polarity=true` means positive
/// control; `false` means anti-control (ctrl=0 triggers).
fn mcx2_polar(
    b: &mut Builder,
    c1: QubitId, p1: bool,
    c2: QubitId, p2: bool,
    target: QubitId,
) {
    if !p1 { b.x(c1); }
    if !p2 { b.x(c2); }
    b.ccx(c1, c2, target);
    if !p2 { b.x(c2); }
    if !p1 { b.x(c1); }
}

/// 3-controlled X with per-control polarity. Uses a borrowed scratch qubit
/// (must be supplied clean, returns clean).
fn mcx3_polar(
    b: &mut Builder,
    c1: QubitId, p1: bool,
    c2: QubitId, p2: bool,
    c3: QubitId, p3: bool,
    target: QubitId,
    scratch: QubitId,
) {
    if !p1 { b.x(c1); }
    if !p2 { b.x(c2); }
    if !p3 { b.x(c3); }
    b.ccx(c1, c2, scratch);
    b.ccx(scratch, c3, target);
    b.ccx(c1, c2, scratch);
    if !p3 { b.x(c3); }
    if !p2 { b.x(c2); }
    if !p1 { b.x(c1); }
}

/// flag ^= (v == 0).  Uses cmp_neq_zero_into internally.
fn cmp_eq_zero_into(b: &mut Builder, v: &[QubitId], flag: QubitId) {
    b.x(flag);
    cmp_neq_zero_into(b, v, flag);
}

/// flag ^= (u > v).  Symmetric to cmp_lt_into(v, u, flag).
fn cmp_gt_into(b: &mut Builder, u: &[QubitId], v: &[QubitId], flag: QubitId) {
    cmp_lt_into(b, v, u, flag);
}

/// Controlled n-bit subtract mod 2^n: if ctrl, acc -= a. Both are n-wide
/// qubit slices. Not a mod-p operation.
fn cucc_sub_ctrl(b: &mut Builder, a: &[QubitId], acc: &[QubitId], ctrl: QubitId) {
    let n = a.len();
    let tmp = b.alloc_qubits(n);
    for i in 0..n { b.ccx(ctrl, a[i], tmp[i]); }
    sub_nbit_qq(b, &tmp, acc);
    for i in 0..n { b.ccx(ctrl, a[i], tmp[i]); }
    b.assert_zero_and_free_vec(&tmp);
}

/// Controlled n-bit add mod 2^n: if ctrl, acc += a.
fn cucc_add_ctrl(b: &mut Builder, a: &[QubitId], acc: &[QubitId], ctrl: QubitId) {
    let n = a.len();
    let tmp = b.alloc_qubits(n);
    for i in 0..n { b.ccx(ctrl, a[i], tmp[i]); }
    add_nbit_qq(b, &tmp, acc);
    for i in 0..n { b.ccx(ctrl, a[i], tmp[i]); }
    b.assert_zero_and_free_vec(&tmp);
}

/// Controlled shift-right by 1 of an n-bit register. ASSUMES v[0]=0 when
/// ctrl=1 (so no information is lost). Implemented as a controlled swap
/// cascade: if ctrl=1, new v[i] = old v[i+1] for i < n-1, new v[n-1] = 0.
fn c_shift_right_1(b: &mut Builder, v: &[QubitId], ctrl: QubitId) {
    let n = v.len();
    for i in 0..(n - 1) {
        cswap(b, ctrl, v[i], v[i + 1]);
    }
}

/// Unconditional shift-left by 1 of an (n+1)-bit register. ASSUMES r[n]=0
/// before the shift. After the shift: r[0]=0, r[i] = old r[i-1] for i ∈ [1, n].
fn shift_left_1(b: &mut Builder, r: &[QubitId]) {
    let n1 = r.len();  // n+1
    // Swap r[n] ↔ r[0] first: r[0] gets the known-0 top bit.
    b.swap(r[n1 - 1], r[0]);
    // Then propagate: swap r[n] ↔ r[n-1], r[n-1] ↔ r[n-2], ..., r[2] ↔ r[1].
    for i in (2..n1).rev() {
        b.swap(r[i], r[i - 1]);
    }
}

/// Inverse of `shift_left_1`: shifts an (n+1)-bit register right by 1.
/// ASSUMES r[0]=0 before the shift (i.e., was even).
#[allow(dead_code)]
fn shift_right_1(b: &mut Builder, r: &[QubitId]) {
    let n1 = r.len();
    for i in 2..n1 {
        b.swap(r[i], r[i - 1]);
    }
    b.swap(r[n1 - 1], r[0]);
}

/// flag ^= (r > c).  r is (n+1)-wide; c is a compile-time constant.
/// Non-destructive: r is restored at the end.
fn cmp_gt_const_n1(b: &mut Builder, r: &[QubitId], c: U256, flag: QubitId) {
    let n1 = r.len();
    let c_plus_1 = c.wrapping_add(U256::from(1));
    sub_nbit_const(b, r, c_plus_1);
    // If r - (c+1) >= 0 (top bit 0), then r > c.
    b.x(r[n1 - 1]);
    b.cx(r[n1 - 1], flag);
    b.x(r[n1 - 1]);
    add_nbit_const(b, r, c_plus_1);
}

/// Classical modular inverse via Fermat's little theorem. Used ONLY at
/// circuit-construction time to compute correction constants.
#[allow(dead_code)]
fn classical_modinv(a: U256, p: U256) -> U256 {
    // a^(p-2) mod p via square-and-multiply.
    let exponent = p.wrapping_sub(U256::from(2));
    let mut result = U256::from(1);
    let mut base = a % p;
    for i in 0..256 {
        if exponent.bit(i) {
            result = mulmod(result, base, p);
        }
        base = mulmod(base, base, p);
    }
    result
}

/// Classical modular multiplication used to compute correction constants
/// at build time.
fn mulmod(a: U256, b: U256, p: U256) -> U256 {
    // Naive (a * b) mod p — both < p < 2^256, so the product may overflow
    // 256 bits. Use U256's widening mul if available; else do it in u512
    // via chunks. alloy's U256 has `mul_mod`.
    a.mul_mod(b, p)
}

// ═══════════════════════════════════════════════════════════════════════════
//  Kaliski binary almost-inverse (qrisp-style, standard form)
// ═══════════════════════════════════════════════════════════════════════════
//
// Faithful port of `kaliski_mod_inv` from the qrisp reference at
// `quantum-elliptic-curve-logarithm/src/quantum/ec_arithmetic.py`.
//
// The function computes `v_in := v_in^{-1} mod p` in place, using a
// self-contained scratch region that is zeroed at function exit. Every
// per-iteration ancilla is uncomputed via the `conjugate` pattern or via
// classical invariants (e.g. `a ^= NOT s[0]` at the end of each iteration).
//
// Difference from qrisp: we work in STANDARD form, no Montgomery
// conversion. The final r register holds `-v_orig^{-1} * 2^{2n} mod p`
// instead of the Montgomery version. We compensate via a single in-place
// classical-constant multiplication by K = (2^{-2n}) mod p at function
// end, which gets us back to v_orig^{-1}.
//
// Assumption: v_in is a nonzero element of (Z/p)*. The test harness
// filters out the v_orig = 0 case before calling `build`, so we skip the
// two phase-fix blocks that qrisp needs for v_orig = 0.

/// Emit the inner iteration body. Takes the persistent state as parameters.
/// Per-iteration transients (`is_zero`, `l_gt`) are allocated and freed
/// WITHIN this function, via the conjugate pattern. The persistent flags
/// `a_f, b_f, add_f` carry no data across iterations (each iteration resets
/// them via classical uncomputation).
fn kaliski_iteration(
    b: &mut Builder,
    p: U256,
    u: &[QubitId],
    v_w: &[QubitId],
    r: &[QubitId],
    s: &[QubitId],
    m_i: QubitId,
    f: QubitId,
    a_f: QubitId,
    b_f: QubitId,
    add_f: QubitId,
) {
    let n = u.len();
    let n1 = r.len();  // n+1

    // ─── STEP 0: is_zero = (v_w == 0);  m[i] ^= (f AND is_zero);  f ^= m[i] ───
    let is_zero = b.alloc_qubit();
    // compute is_zero ^= (v_w == 0)
    cmp_eq_zero_into(b, v_w, is_zero);
    b.ccx(f, is_zero, m_i);
    // uncompute is_zero (run the same function; cmp_eq_zero is self-inverse
    // as a gate sequence because each step is involutive).
    cmp_eq_zero_into(b, v_w, is_zero);
    b.assert_zero_and_free(is_zero);
    b.cx(m_i, f);

    // ─── STEP 1 ───
    //   a ^= (f=1 AND u[0]=0)
    //   m[i] ^= (f=1 AND a=0 AND v_w[0]=0)
    //   b ^= a; b ^= m[i]
    mcx2_polar(b, f, true, u[0], false, a_f);
    {
        // Borrow a scratch qubit for the 3-control mcx.
        let scratch = b.alloc_qubit();
        mcx3_polar(b, f, true, a_f, false, v_w[0], false, m_i, scratch);
        b.assert_zero_and_free(scratch);
    }
    b.cx(a_f, b_f);
    b.cx(m_i, b_f);

    // ─── STEP 2: with conjugate(l = u > v_w): ... ───
    //   a ^= (f=1 AND l=1 AND b=0)
    //   m[i] ^= (f=1 AND l=1 AND b=0)
    let l_gt = b.alloc_qubit();
    cmp_gt_into(b, u, v_w, l_gt);
    {
        let scratch = b.alloc_qubit();
        mcx3_polar(b, f, true, l_gt, true, b_f, false, a_f, scratch);
        mcx3_polar(b, f, true, l_gt, true, b_f, false, m_i, scratch);
        b.assert_zero_and_free(scratch);
    }
    // Uncompute l_gt (cmp_gt is self-inverse as a gate sequence).
    cmp_gt_into(b, u, v_w, l_gt);
    b.assert_zero_and_free(l_gt);

    // ─── STEP 3: with control(a): swap(u, v_w); swap(r, s) ───
    for j in 0..n { cswap(b, a_f, u[j], v_w[j]); }
    for j in 0..n1 { cswap(b, a_f, r[j], s[j]); }

    // ─── STEP 4 ───
    //   add ^= (f=1 AND b=0)
    //   with control(add): v_w -= u; s += r
    mcx2_polar(b, f, true, b_f, false, add_f);
    // v_w -= u mod 2^n, controlled by add_f.
    cucc_sub_ctrl(b, u, v_w, add_f);
    // s += r mod 2^(n+1), controlled by add_f.
    cucc_add_ctrl(b, r, s, add_f);

    // ─── STEP 5: uncompute add; uncompute b ───
    mcx2_polar(b, f, true, b_f, false, add_f);
    b.cx(m_i, b_f);
    b.cx(a_f, b_f);

    // ─── STEP 6: v_w := v_w / 2 (shift right by 1), controlled by f ───
    // At this point, if f=1 then v_w is even (low bit 0).
    c_shift_right_1(b, v_w, f);

    // ─── STEP 7: r := 2*r (shift left by 1 in (n+1) bits) ───
    shift_left_1(b, r);

    // ─── STEP 8: if r ≥ p: r -= p (Solinas fold). r is even after shift_left
    //   and p is odd, so `r = p` never occurs and `r ≥ p` ≡ `r > p`. Cheaper
    //   than the prior cmp_gt_const_n1 + csub_nbit_const(p) pair.
    {
        let n1 = r.len();
        let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));
        // r' = r + c (fits in n+1 bits since r < 2p and c = 2^n - p).
        add_nbit_const(b, r, c);
        let flag = b.alloc_qubit();
        b.cx(r[n1 - 1], flag);   // flag := top bit of r' = (r ≥ p)
        // If flag=0: undo the add of c (we don't reduce).
        b.x(flag);
        csub_nbit_const(b, r, c, flag);
        b.x(flag);
        // If flag=1: clear top bit, giving r + c - 2^n = r - p.
        b.cx(flag, r[n1 - 1]);
        // Uncompute flag via parity (r was even pre-step; r-p odd if reduced).
        b.cx(r[0], flag);
        b.assert_zero_and_free(flag);
    }

    // ─── STEP 9: with control(a): swap(u, v_w); swap(r, s) (again) ───
    for j in 0..n { cswap(b, a_f, u[j], v_w[j]); }
    for j in 0..n1 { cswap(b, a_f, r[j], s[j]); }

    // ─── STEP 10: uncompute a via `a ^= NOT s[0]` ───
    // After STEP 9's swap, the invariant (from qrisp) is that
    //   a == NOT s[0]
    // Hence `cx(NOT s[0], a)` zeros a.
    b.x(s[0]);
    b.cx(s[0], a_f);
    b.x(s[0]);
}

/// In-place classical-constant multiplication: v := v * c mod p.
///
/// Uses the standard compute-in-fresh-then-uncompute pattern:
///   1. tmp = 0
///   2. tmp += v * c                         (shift-and-add, classical c)
///   3. v -= tmp * c^{-1} = v - v*c*c^{-1} = 0  (classical c^{-1})
///   4. swap v, tmp
///   5. free tmp
fn in_place_mul_const(b: &mut Builder, v: &[QubitId], c: U256, p: U256) {
    let n = v.len();
    let tmp = b.alloc_qubits(n);
    mul_by_const_acc(b, v, c, &tmp, p, false);       // tmp += v * c
    let c_inv = classical_modinv(c, p);
    mul_by_const_acc(b, &tmp, c_inv, v, p, true);    // v -= tmp * c_inv
    for i in 0..n { b.swap(v[i], tmp[i]); }
    b.assert_zero_and_free_vec(&tmp);
}

/// `acc ±= x * c mod p`. `c` is a classical constant. Does NOT fold acc.
/// Maintains a doubling copy of x in a temp register; adds it to acc at
/// positions where c has a bit set.
fn mul_by_const_acc(
    b: &mut Builder,
    x: &[QubitId],
    c: U256,
    acc: &[QubitId],
    p: U256,
    subtract: bool,
) {
    let n = x.len();
    if c == U256::ZERO { return; }

    // tmp := x  (via CX copy)
    let tmp = b.alloc_qubits(n);
    for i in 0..n { b.cx(x[i], tmp[i]); }

    // Iterate bits of c from LSB to MSB. At step i, tmp holds x * 2^i mod p.
    // Add tmp to acc if bit i of c is set. Then double tmp for the next step.
    //
    // We iterate up through the highest set bit of c, plus any trailing zero
    // bits (we must double enough times to make uncomputation clean).
    let mut top = 0usize;
    for i in 0..256 {
        if bit(c, i) { top = i; }
    }

    for i in 0..=top {
        if bit(c, i) {
            if subtract {
                mod_sub_qq(b, acc, &tmp, p);
            } else {
                mod_add_qq(b, acc, &tmp, p);
            }
        }
        if i < top {
            mod_double_inplace(b, &tmp, p);
        }
    }

    // At this point tmp = x * 2^top mod p. Halve it back `top` times to
    // recover x, then uncompute via cx.
    for _ in 0..top {
        mod_halve_inplace(b, &tmp, p);
    }
    for i in 0..n { b.cx(x[i], tmp[i]); }
    b.assert_zero_and_free_vec(&tmp);
}

/// Persistent state for the Kaliski forward computation. Transients are
/// allocated inside the iteration body; `emit_inverse` will correctly
/// reverse them because it skips R ops (the free markers) in the reverse
/// stream, and our forward guarantees each free lands on a |0⟩ qubit.
struct KaliskiState {
    u: Vec<QubitId>,       // n qubits
    v_w: Vec<QubitId>,     // n qubits
    r: Vec<QubitId>,       // n+1 qubits
    s: Vec<QubitId>,       // n+1 qubits
    m_hist: Vec<QubitId>,  // 2n qubits
    f_flag: QubitId,
    a_flag: QubitId,
    b_flag: QubitId,
    add_flag: QubitId,
}

fn alloc_kaliski_state(b: &mut Builder, n: usize) -> KaliskiState {
    KaliskiState {
        u: b.alloc_qubits(n),
        v_w: b.alloc_qubits(n),
        r: b.alloc_qubits(n + 1),
        s: b.alloc_qubits(n + 1),
        m_hist: b.alloc_qubits(2 * n),
        f_flag: b.alloc_qubit(),
        a_flag: b.alloc_qubit(),
        b_flag: b.alloc_qubit(),
        add_flag: b.alloc_qubit(),
    }
}

fn free_kaliski_state(b: &mut Builder, st: KaliskiState) {
    b.assert_zero_and_free(st.add_flag);
    b.assert_zero_and_free(st.b_flag);
    b.assert_zero_and_free(st.a_flag);
    b.assert_zero_and_free(st.f_flag);
    b.assert_zero_and_free_vec(&st.m_hist);
    b.assert_zero_and_free_vec(&st.s);
    b.assert_zero_and_free_vec(&st.r);
    b.assert_zero_and_free_vec(&st.v_w);
    b.assert_zero_and_free_vec(&st.u);
}

/// Forward-only Kaliski computation. Reads `v_in` (never writes), populates
/// `st.*` with the algorithm's intermediate state. After this returns:
///   - `v_in` is unchanged
///   - `st.r[..n]` holds the RAW Kaliski inverse `v^{-1} * 2^{2n} mod p`
///   - everything else in `st` is populated with deterministic iteration history
///
/// The caller is responsible for applying the classical correction factor
/// `K = 2^{-2n} mod p` and for calling `emit_inverse(kaliski_forward)` to
/// restore `st.*` to all zero.
fn kaliski_forward(b: &mut Builder, v_in: &[QubitId], st: &KaliskiState, p: U256) {
    let n = v_in.len();

    // ─── Init ───
    // u := p (classical load)
    for i in 0..n { if bit(p, i) { b.x(st.u[i]); } }
    // v_w := v_in  (CX-copy; v_in unchanged)
    for i in 0..n { b.cx(v_in[i], st.v_w[i]); }
    // s := 1
    b.x(st.s[0]);
    // f := 1
    b.x(st.f_flag);

    // ─── 2n iterations ───
    for i in 0..(2 * n) {
        kaliski_iteration(
            b, p, &st.u, &st.v_w, &st.r, &st.s,
            st.m_hist[i],
            st.f_flag, st.a_flag, st.b_flag, st.add_flag,
        );
    }

    // After the loop for nonzero v_in, classical invariants give:
    //   u = 1, v_w = 0, f = 0, a = b = add = 0
    //   r = raw coefficient (related to -v^{-1} * 2^{2n})
    //   s = some coefficient
    // Apply inpl_rsub to r so st.r contains the POSITIVE raw inverse form:
    //   r := (p - r) mod 2^(n+1)
    // via `x(r); add_nbit_const(r, p+1)` on the full (n+1)-bit register.
    for &q in &st.r { b.x(q); }
    add_nbit_const(b, &st.r, p.wrapping_add(U256::from(1)));
}

/// Compute `output ^= v_in^{-1} mod p` without mutating `v_in`. `st` and
/// `output` are caller-provided scratch (both |0⟩ on entry) and are left
/// in their algorithm-intermediate state: `st` returns to |0⟩, `output`
/// holds the inverse. The caller uncomputes `output` by running this
/// same function under `emit_inverse`.
fn kal_compute_into(
    b: &mut Builder,
    v_in: &[QubitId],
    output: &[QubitId],
    st: &KaliskiState,
    p: U256,
) {
    let n = v_in.len();
    // Forward pass: st.r[..n] holds the raw form = inverse * 2^(2n) mod p.
    kaliski_forward(b, v_in, st, p);
    // Copy the raw form into output (output ^= raw).
    for i in 0..n { b.cx(st.r[i], output[i]); }
    // Uncompute st entirely.
    emit_inverse(b, |b| kaliski_forward(b, v_in, st, p));
    // Correction: output *= 2^(-2n) mod p ≡ halve output 2n times. Cheaper
    // than the generic in_place_mul_const(output, k_const) by a factor of
    // ~2× because mod_halve_inplace is just one Solinas reduction.
    for _ in 0..(2 * n) { mod_halve_inplace(b, output, p); }
}

fn kaliski_inv_inplace(b: &mut Builder, v_in: &[QubitId], p: U256) {
    let n = v_in.len();

    // Bennett compute-copy-uncompute pattern. Each call of
    // `kaliski_inv_inplace` maps v_in ↔ v_in^{-1} (involution), with
    // internal scratch fully zeroed by function end.
    let st = alloc_kaliski_state(b, n);
    let output = b.alloc_qubits(n);

    // ─── Phase 1: compute inverse of v_in into output ───
    kaliski_forward(b, v_in, &st, p);
    // st.r[..n] now holds raw inverse (in mod 2p, low n bits).
    // Apply classical correction: st.r[..n] *= K mod p, where K = 2^{-2n} mod p.
    let two_2n = pow_mod_2_k(p, 2 * n);
    let k_const = classical_modinv(two_2n, p);
    in_place_mul_const(b, &st.r[..n], k_const, p);
    // Copy exact inverse into output.
    for i in 0..n { b.cx(st.r[i], output[i]); }
    // Undo the correction: st.r[..n] *= K^{-1} mod p.
    in_place_mul_const(b, &st.r[..n], two_2n, p);
    // Now st is back to its post-kaliski_forward state. Reverse the forward.
    emit_inverse(b, |b| kaliski_forward(b, v_in, &st, p));
    // st is all 0 again. v_in unchanged. output = v_in^{-1}.

    // Swap v_in and output.
    for i in 0..n { b.swap(v_in[i], output[i]); }
    // v_in = inverse, output = v_orig.

    // ─── Phase 2: zero output via a second Bennett pass ───
    // Compute inverse of current v_in (which is v_orig^{-1}), = v_orig,
    // and XOR it into output. Since output currently = v_orig, the XOR
    // zeroes output.
    kaliski_forward(b, v_in, &st, p);
    in_place_mul_const(b, &st.r[..n], k_const, p);
    for i in 0..n { b.cx(st.r[i], output[i]); }   // output ^= v_orig = 0
    in_place_mul_const(b, &st.r[..n], two_2n, p);
    emit_inverse(b, |b| kaliski_forward(b, v_in, &st, p));
    // st all 0, output all 0 (hopefully), v_in = inverse.

    b.assert_zero_and_free_vec(&output);
    free_kaliski_state(b, st);
}

/// Classical: compute `2^k mod p`.
fn pow_mod_2_k(p: U256, k: usize) -> U256 {
    let mut r = U256::from(1);
    let two = U256::from(2);
    for _ in 0..k {
        r = mulmod(r, two, p);
    }
    r
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

    // Pair 1 (folded): keep tx holding dx throughout, compute dx^{-1} into
    // `inv` ancilla, use it, then uncompute. Replaces two kaliski_inv_inplace
    // involutions (≈4 Bennett passes) with one kal_compute_into and its
    // emit_inverse (≈2 Bennett passes).
    {
        let st1 = alloc_kaliski_state(b, N);
        let inv = b.alloc_qubits(N);
        kal_compute_into(b, &tx, &inv, &st1, p);     // inv = dx^{-1}
        mod_mul_add_qq(b, &lam, &ty, &inv, p);       // lam += dy · dx^{-1} = λ
        mod_mul_sub_qq(b, &ty, &lam, &tx, p);        // Py -= λ·dx = 0
        emit_inverse(b, |b| kal_compute_into(b, &tx, &inv, &st1, p));
        b.assert_zero_and_free_vec(&inv);
        free_kaliski_state(b, st1);
    }

    // Px := λ² - Px_orig - Qx
    mod_mul_sub_qq(b, &tx, &lam, &lam, p);
    mod_neg_inplace(b, &tx, p);
    mod_sub_qb(b, &tx, &ox, p);
    mod_sub_qb(b, &tx, &ox, p);

    // Py := λ·Qx − λ·Rx − Qy
    mod_mul_add_qb(b, &ty, &lam, &ox, p);
    mod_mul_sub_qq(b, &ty, &lam, &tx, p);
    mod_sub_qb(b, &ty, &oy, p);

    // Uncompute lam using λ = (Qy + Ry) / (Qx - Rx).
    mod_sub_qb(b, &tx, &ox, p);
    mod_neg_inplace(b, &tx, p);                   // tx = Qx - Rx
    // Pair 2 (folded): keep tx = Qx-Rx throughout, put its inverse in `inv`.
    {
        let st2 = alloc_kaliski_state(b, N);
        let inv = b.alloc_qubits(N);
        kal_compute_into(b, &tx, &inv, &st2, p);
        mod_mul_sub_qq(b, &lam, &ty, &inv, p);
        mod_mul_sub_qb(b, &lam, &inv, &oy, p);
        emit_inverse(b, |b| kal_compute_into(b, &tx, &inv, &st2, p));
        b.assert_zero_and_free_vec(&inv);
        free_kaliski_state(b, st2);
    }
    mod_neg_inplace(b, &tx, p);                   // tx = -(Qx-Rx) = Rx - Qx
    mod_add_qb(b, &tx, &ox, p);                   // tx = Rx

    b.assert_zero_and_free_vec(&lam);

    Layout { target_x, target_y, offset_x, offset_y }
}


