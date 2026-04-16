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

use crate::circuit::{BitId, Op, OperationType, QubitId, RegisterId};

struct B {
    pub ops: Vec<Op>,
    pub next_qubit: u32,
    pub next_bit: u32,
    pub next_register: u32,
    pub free_qubits: Vec<u32>,
}

impl B {
    fn new() -> Self {
        Self { ops: Vec::new(), next_qubit: 0, next_bit: 0, next_register: 0, free_qubits: Vec::new() }
    }
    fn alloc_qubit(&mut self) -> QubitId {
        if let Some(q) = self.free_qubits.pop() { QubitId(q) }
        else { let q = self.next_qubit; self.next_qubit += 1; QubitId(q) }
    }
    fn alloc_qubits(&mut self, n: usize) -> Vec<QubitId> { (0..n).map(|_| self.alloc_qubit()).collect() }
    fn alloc_bit(&mut self) -> BitId { let b = self.next_bit; self.next_bit += 1; BitId(b) }
    fn alloc_bits(&mut self, n: usize) -> Vec<BitId> { (0..n).map(|_| self.alloc_bit()).collect() }
    fn free(&mut self, q: QubitId) { self.r(q); self.free_qubits.push(q.0); }
    fn free_vec(&mut self, qs: &[QubitId]) { for &q in qs { self.free(q); } }
    fn declare_qubit_register(&mut self, qs: &[QubitId]) {
        let r = RegisterId(self.next_register); self.next_register += 1;
        for &q in qs { let mut op = Op::empty(); op.kind = OperationType::AppendToRegister; op.q_target = q; op.r_target = r; self.ops.push(op); }
        let mut op = Op::empty(); op.kind = OperationType::Register; op.r_target = r; self.ops.push(op);
    }
    fn declare_bit_register(&mut self, bs: &[BitId]) {
        let r = RegisterId(self.next_register); self.next_register += 1;
        for &b in bs { let mut op = Op::empty(); op.kind = OperationType::AppendToRegister; op.c_target = b; op.r_target = r; self.ops.push(op); }
        let mut op = Op::empty(); op.kind = OperationType::Register; op.r_target = r; self.ops.push(op);
    }
    fn x(&mut self, q: QubitId) { let mut op = Op::empty(); op.kind = OperationType::X; op.q_target = q; self.ops.push(op); }
    fn z(&mut self, q: QubitId) { let mut op = Op::empty(); op.kind = OperationType::Z; op.q_target = q; self.ops.push(op); }
    fn cx(&mut self, ctrl: QubitId, tgt: QubitId) { let mut op = Op::empty(); op.kind = OperationType::CX; op.q_control1 = ctrl; op.q_target = tgt; self.ops.push(op); }
    fn cz(&mut self, a: QubitId, b: QubitId) { let mut op = Op::empty(); op.kind = OperationType::CZ; op.q_control1 = a; op.q_target = b; self.ops.push(op); }
    fn ccx(&mut self, c1: QubitId, c2: QubitId, tgt: QubitId) { let mut op = Op::empty(); op.kind = OperationType::CCX; op.q_control2 = c1; op.q_control1 = c2; op.q_target = tgt; self.ops.push(op); }
    fn ccz(&mut self, c1: QubitId, c2: QubitId, tgt: QubitId) { let mut op = Op::empty(); op.kind = OperationType::CCZ; op.q_control2 = c1; op.q_control1 = c2; op.q_target = tgt; self.ops.push(op); }
    fn swap(&mut self, a: QubitId, b: QubitId) { let mut op = Op::empty(); op.kind = OperationType::Swap; op.q_control1 = a; op.q_target = b; self.ops.push(op); }
    fn r(&mut self, q: QubitId) { let mut op = Op::empty(); op.kind = OperationType::R; op.q_target = q; self.ops.push(op); }
    fn x_if(&mut self, q: QubitId, cond: BitId) { let mut op = Op::empty(); op.kind = OperationType::X; op.q_target = q; op.c_condition = cond; self.ops.push(op); }
    fn cx_if(&mut self, ctrl: QubitId, tgt: QubitId, cond: BitId) { let mut op = Op::empty(); op.kind = OperationType::CX; op.q_control1 = ctrl; op.q_target = tgt; op.c_condition = cond; self.ops.push(op); }
    fn ccx_if(&mut self, c1: QubitId, c2: QubitId, tgt: QubitId, cond: BitId) { let mut op = Op::empty(); op.kind = OperationType::CCX; op.q_control2 = c1; op.q_control1 = c2; op.q_target = tgt; op.c_condition = cond; self.ops.push(op); }
    fn push_condition(&mut self, cond: BitId) { let mut op = Op::empty(); op.kind = OperationType::PushCondition; op.c_condition = cond; self.ops.push(op); }
    fn pop_condition(&mut self) { let mut op = Op::empty(); op.kind = OperationType::PopCondition; self.ops.push(op); }
    // ── Measurement / phase / classical bit ops ──
    fn hmr(&mut self, q: QubitId, c: BitId) { let mut op = Op::empty(); op.kind = OperationType::Hmr; op.q_target = q; op.c_target = c; self.ops.push(op); }
    fn neg(&mut self) { let mut op = Op::empty(); op.kind = OperationType::Neg; self.ops.push(op); }
    fn bit_invert(&mut self, c: BitId) { let mut op = Op::empty(); op.kind = OperationType::BitInvert; op.c_target = c; self.ops.push(op); }
    fn bit_store0(&mut self, c: BitId) { let mut op = Op::empty(); op.kind = OperationType::BitStore0; op.c_target = c; self.ops.push(op); }
    fn bit_store1(&mut self, c: BitId) { let mut op = Op::empty(); op.kind = OperationType::BitStore1; op.c_target = c; self.ops.push(op); }
    // ── Classically-conditioned variants for all remaining gates ──
    fn z_if(&mut self, q: QubitId, cond: BitId) { let mut op = Op::empty(); op.kind = OperationType::Z; op.q_target = q; op.c_condition = cond; self.ops.push(op); }
    fn cz_if(&mut self, a: QubitId, b: QubitId, cond: BitId) { let mut op = Op::empty(); op.kind = OperationType::CZ; op.q_control1 = a; op.q_target = b; op.c_condition = cond; self.ops.push(op); }
    fn ccz_if(&mut self, c1: QubitId, c2: QubitId, tgt: QubitId, cond: BitId) { let mut op = Op::empty(); op.kind = OperationType::CCZ; op.q_control2 = c1; op.q_control1 = c2; op.q_target = tgt; op.c_condition = cond; self.ops.push(op); }
    fn swap_if(&mut self, a: QubitId, b: QubitId, cond: BitId) { let mut op = Op::empty(); op.kind = OperationType::Swap; op.q_control1 = a; op.q_target = b; op.c_condition = cond; self.ops.push(op); }
    fn neg_if(&mut self, cond: BitId) { let mut op = Op::empty(); op.kind = OperationType::Neg; op.c_condition = cond; self.ops.push(op); }
    fn hmr_if(&mut self, q: QubitId, c: BitId, cond: BitId) { let mut op = Op::empty(); op.kind = OperationType::Hmr; op.q_target = q; op.c_target = c; op.c_condition = cond; self.ops.push(op); }
    fn bit_invert_if(&mut self, c: BitId, cond: BitId) { let mut op = Op::empty(); op.kind = OperationType::BitInvert; op.c_target = c; op.c_condition = cond; self.ops.push(op); }
    fn bit_store0_if(&mut self, c: BitId, cond: BitId) { let mut op = Op::empty(); op.kind = OperationType::BitStore0; op.c_target = c; op.c_condition = cond; self.ops.push(op); }
    fn bit_store1_if(&mut self, c: BitId, cond: BitId) { let mut op = Op::empty(); op.kind = OperationType::BitStore1; op.c_target = c; op.c_condition = cond; self.ops.push(op); }
    fn r_if(&mut self, q: QubitId, cond: BitId) { let mut op = Op::empty(); op.kind = OperationType::R; op.q_target = q; op.c_condition = cond; self.ops.push(op); }
    // ── Gidney measurement-based AND uncomputation (convenience) ──
    // Uncomputes `tgt = c1 AND c2` using HMR + phase feedback.
    // Cost: 0 Toffoli (1 HMR + 1 classically-conditioned CZ).
    // Precondition: tgt holds (c1 AND c2) computed by a prior CCX.
    fn uncompute_and(&mut self, c1: QubitId, c2: QubitId, tgt: QubitId) {
        let m = self.alloc_bit();
        self.hmr(tgt, m);
        self.cz_if(c1, c2, m);
        self.neg_if(m);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  emit_inverse: run a closure, pop the ops it emitted, and re-emit them
//  reversed.
//
//  The closure may contain `alloc_qubit` / `free` calls;
//  the R ops that `free` produces are SKIPPED during
//  reverse replay. This relies on the forward being "clean" — i.e. each
//  free lands on a qubit that the forward gates already drove to |0⟩
//  before the R. Under that invariant, the reverse gate sequence brings
//  the same qubit back to |0⟩ at the "alloc" point (pre-forward-allocation),
//  and the R we skipped is unnecessary.
//
//  The forward's internal alloc/free bookkeeping in the B's free
//  pool is NOT undone by the reverse — the pool state at reverse exit
//  equals the pool state at forward exit. Subsequent allocations in the
//  parent scope reuse those qubit IDs, seeing them at |0⟩ (as zeroed by
//  the reverse gate sequence).
// ═══════════════════════════════════════════════════════════════════════════
fn emit_inverse<F: FnOnce(&mut B)>(b: &mut B, f: F) {
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
fn conjugate<F, G>(b: &mut B, compute: F, body: G)
where
    F: Fn(&mut B),
    G: FnOnce(&mut B),
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

fn maj(b: &mut B, x: QubitId, y: QubitId, w: QubitId) {
    b.cx(w, y);
    b.cx(w, x);
    b.ccx(x, y, w);
}

fn uma(b: &mut B, x: QubitId, y: QubitId, w: QubitId) {
    b.ccx(x, y, w);
    b.cx(w, x);
    b.cx(x, y);
}

/// Fast Cuccaro add using carry ancillae + measurement-based UMA.
/// Same interface as `cuccaro_add` but uses n-1 carry ancillae so the
/// UMA sweep costs 0 Toffoli (measurement only). NOT emit_inverse-safe.
fn cuccaro_add_fast(b: &mut B, a: &[QubitId], acc: &[QubitId], c_in: QubitId) {
    let n = a.len();
    assert_eq!(n, acc.len());
    if n == 0 { return; }
    if n == 1 {
        b.cx(c_in, acc[0]);
        b.cx(a[0], acc[0]);
        return;
    }

    let carries = b.alloc_qubits(n - 1);

    // Forward MAJ sweep with carry ancillae.
    // Step 0: MAJ(c_in, acc[0], a[0]) → carry into carries[0]
    b.cx(a[0], acc[0]);
    b.cx(a[0], c_in);
    b.ccx(c_in, acc[0], carries[0]);
    b.cx(carries[0], a[0]);
    // Steps 1..n-2: MAJ(a[i-1], acc[i], a[i]) → carry into carries[i]
    for i in 1..n - 1 {
        b.cx(a[i], acc[i]);
        b.cx(a[i], a[i - 1]);
        b.ccx(a[i - 1], acc[i], carries[i]);
        b.cx(carries[i], a[i]);
    }

    // Final sum bit (same as original cuccaro_add)
    b.cx(a[n - 2], acc[n - 1]);
    b.cx(a[n - 1], acc[n - 1]);

    // Backward UMA sweep with measurement-based carry uncompute (0 Toffoli).
    for i in (1..n - 1).rev() {
        b.cx(carries[i], a[i]);
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        b.cz_if(a[i - 1], acc[i], m);
        b.cx(a[i], a[i - 1]);
        b.cx(a[i - 1], acc[i]);
    }
    // Step 0 UMA:
    b.cx(carries[0], a[0]);
    let m0 = b.alloc_bit();
    b.hmr(carries[0], m0);
    b.cz_if(c_in, acc[0], m0);
    b.cx(a[0], c_in);
    b.cx(c_in, acc[0]);

    b.free_vec(&carries);
}

/// In-place addition `acc += a mod 2^n` on quantum n-bit registers.
/// * `c_in` is a fresh ancilla qubit at 0 on entry and returns to 0.
/// * `a` unchanged; `acc` becomes (a + acc) mod 2^n.
/// Pure mod-2^n: the high carry is discarded (no `z` ancilla). This is
/// honestly reversible because the last MAJ/UMA pair cancel out the
/// carry information on `a[n-1]`.
fn cuccaro_add(b: &mut B, a: &[QubitId], acc: &[QubitId], c_in: QubitId) {
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
fn cuccaro_sub(b: &mut B, a: &[QubitId], acc: &[QubitId], c_in: QubitId) {
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

fn inv_maj(b: &mut B, x: QubitId, y: QubitId, w: QubitId) {
    // maj = CX(w,y); CX(w,x); CCX(x,y,w)
    // inv = CCX(x,y,w); CX(w,x); CX(w,y)
    b.ccx(x, y, w);
    b.cx(w, x);
    b.cx(w, y);
}

fn inv_uma(b: &mut B, x: QubitId, y: QubitId, w: QubitId) {
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

fn load_const(b: &mut B, n: usize, c: U256) -> Vec<QubitId> {
    let qs = b.alloc_qubits(n);
    for i in 0..n {
        if bit(c, i) {
            b.x(qs[i]);
        }
    }
    qs
}

fn unload_const(b: &mut B, qs: &[QubitId], c: U256) {
    for i in 0..qs.len() {
        if bit(c, i) {
            b.x(qs[i]);
        }
    }
    b.free_vec(qs);
}

fn load_bits(b: &mut B, bits: &[BitId]) -> Vec<QubitId> {
    let n = bits.len();
    let qs = b.alloc_qubits(n);
    for i in 0..n {
        // qs[i] ← bits[i] via conditional X
        b.x_if(qs[i], bits[i]);
    }
    qs
}

fn unload_bits(b: &mut B, qs: &[QubitId], bits: &[BitId]) {
    for i in 0..qs.len() {
        b.x_if(qs[i], bits[i]);
    }
    b.free_vec(qs);
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
fn ext_reg(b: &mut B, reg: &[QubitId]) -> (Vec<QubitId>, QubitId) {
    let ovf = b.alloc_qubit();
    let mut r = reg.to_vec();
    r.push(ovf);
    (r, ovf)
}

/// Release the overflow ancilla (which must be 0 on exit).
fn unext_reg(b: &mut B, ovf: QubitId) {
    b.free(ovf);
}

/// `acc := (acc + a) mod p`. Both `acc` and `a` are n-bit quantum registers
/// with value in [0, p). Solinas reduction using c = 2^n - p: sum ∈ [0, 2p),
/// then add c, branch on top bit to either clear it (reduction) or undo
/// the add (no reduction). Saves one full (n+1)-wide Cuccaro compared to
/// the sub-p/add-p/csub-p pattern.
fn mod_add_qq(b: &mut B, acc: &[QubitId], a: &[QubitId], p: U256) {
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
    b.free(flag);

    unext_reg(b, a_ovf);
    unext_reg(b, acc_ovf);
    let _ = (acc_ext, a_ext);
}

fn mod_sub_qq(b: &mut B, acc: &[QubitId], a: &[QubitId], p: U256) {
    // mod_add_qq is a bijection on (acc, a): (acc, a) ↦ (acc + a mod p, a).
    // Its gate-level inverse therefore acts as (acc, a) ↦ (acc - a mod p, a),
    // which is exactly what we want. emit_inverse replays the forward's gates
    // reversed, skipping R markers — valid because mod_add_qq is clean
    // (every ancilla is driven to |0⟩ before its R).
    let a_copy: Vec<QubitId> = a.to_vec();
    emit_inverse(b, move |b| mod_add_qq(b, acc, &a_copy, p));
}

/// Fast `acc := (acc - a) mod p`. Direct sub + conditional add-p + flag
/// uncompute via neg+cmp_lt+neg. All ops use measurement-based Cuccaro.
fn mod_sub_qq_fast(b: &mut B, acc: &[QubitId], a: &[QubitId], p: U256) {
    let n = acc.len();
    assert_eq!(n, a.len());
    debug_assert_eq!(n, 256);

    let (acc_ext, acc_ovf) = ext_reg(b, acc);
    let (a_ext, a_ovf) = ext_reg(b, a);

    // Step 1: (n+1)-bit sub.
    sub_nbit_qq_fast(b, &a_ext, &acc_ext);

    // Step 2: flag = acc_ovf (=1 iff underflow, i.e. acc < a).
    let flag = b.alloc_qubit();
    b.cx(acc_ovf, flag);
    // We only need the borrow as a separate flag; the low register is
    // corrected modulo 2^n, so clear the extension bit immediately.
    b.cx(flag, acc_ovf);

    // Step 3: underflow correction. With p = 2^n - c, the wrapped 256-bit
    // subtraction needs only a conditional subtract of c on the low register.
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));
    csub_nbit_const_fast(b, &acc_ext[..n], c, flag);

    // Step 4: uncompute flag. Identity: flag = NOT(acc_final < (p - a)).
    // Negate a in place, compare, un-negate.
    b.x(flag);
    mod_neg_inplace_fast(b, &a_ext[..n], p);
    cmp_lt_into_fast(b, &acc_ext[..n], &a_ext[..n], flag);
    mod_neg_inplace_fast(b, &a_ext[..n], p);
    b.free(flag);

    unext_reg(b, a_ovf);
    unext_reg(b, acc_ovf);
    let _ = (acc_ext, a_ext);
}

/// Fast mod_neg using measurement-based Cuccaro for the addition.
fn mod_neg_inplace_fast(b: &mut B, v: &[QubitId], p: U256) {
    for &q in v { b.x(q); }
    let n = v.len();
    let ca = load_const(b, n, p.wrapping_add(U256::from(1)));
    add_nbit_qq_fast(b, &ca, v);
    unload_const(b, &ca, p.wrapping_add(U256::from(1)));
}

fn mod_add_qc(b: &mut B, acc: &[QubitId], c: U256, p: U256) {
    // acc := (acc + c) mod p. c is a compile-time constant.
    let n = acc.len();
    let a = load_const(b, n, c);
    mod_add_qq_fast(b, acc, &a, p);
    unload_const(b, &a, c);
}

fn mod_sub_qc(b: &mut B, acc: &[QubitId], c: U256, p: U256) {
    // acc := (acc - c) mod p = acc + (p - c) mod p.
    let n = acc.len();
    let c_neg = (p - (c % p)) % p;
    let a = load_const(b, n, c_neg);
    mod_add_qq_fast(b, acc, &a, p);
    unload_const(b, &a, c_neg);
}

fn mod_add_qb(b: &mut B, acc: &[QubitId], bits: &[BitId], p: U256) {
    // acc := (acc + bits) mod p. `bits` is a classical bit register.
    let a = load_bits(b, bits);
    mod_add_qq_fast(b, acc, &a, p);
    unload_bits(b, &a, bits);
}

fn mod_sub_qb(b: &mut B, acc: &[QubitId], bits: &[BitId], p: U256) {
    // acc -= bits mod p. Uses fast mod_sub_qq via neg+add+neg.
    let a = load_bits(b, bits);
    mod_sub_qq_fast(b, acc, &a, p);
    unload_bits(b, &a, bits);
}

/// `v := (p - v) mod p`. Operates on an n-bit register in [0, p).
///
/// Implementation uses the reversible identity:
///     p - v = NOT(v) + (p + 1)         (all arithmetic mod 2^n)
/// which holds because NOT(v) = 2^n - 1 - v, so NOT(v) + p + 1 = 2^n + (p - v).
///
/// For v = 0 the result is p, not 0 (non-canonical but ≡ 0 mod p).
/// EC preconditions (dx, dy nonzero) avoid this case in practice.
fn mod_neg_inplace(b: &mut B, v: &[QubitId], p: U256) {
    for &q in v {
        b.x(q);
    }
    add_nbit_const(b, v, p.wrapping_add(U256::from(1)));
}

// ═══════════════════════════════════════════════════════════════════════════
//  Non-modular n-bit primitives
// ═══════════════════════════════════════════════════════════════════════════

/// Fast Cuccaro sub: `acc -= a mod 2^n` with measurement UMA (0 Toffoli
/// for UMA sweep). Exact gate-level inverse of `cuccaro_add_fast`.
fn cuccaro_sub_fast(b: &mut B, a: &[QubitId], acc: &[QubitId], c_in: QubitId) {
    let n = a.len();
    assert_eq!(n, acc.len());
    if n == 0 { return; }
    if n == 1 {
        b.cx(a[0], acc[0]);
        b.cx(c_in, acc[0]);
        return;
    }

    let carries = b.alloc_qubits(n - 1);

    // Forward inv_UMA sweep with carry ancillae (reversed UMA from cuccaro_sub).
    // Step 0:
    b.cx(c_in, acc[0]);
    b.cx(a[0], c_in);
    b.ccx(c_in, acc[0], carries[0]);
    b.cx(carries[0], a[0]);
    // Steps 1..n-2:
    for i in 1..n - 1 {
        b.cx(a[i - 1], acc[i]);
        b.cx(a[i], a[i - 1]);
        b.ccx(a[i - 1], acc[i], carries[i]);
        b.cx(carries[i], a[i]);
    }

    // Final sum bit (reversed from cuccaro_add)
    b.cx(a[n - 1], acc[n - 1]);
    b.cx(a[n - 2], acc[n - 1]);

    // Backward inv_MAJ sweep with measurement.
    for i in (1..n - 1).rev() {
        b.cx(carries[i], a[i]);
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        b.cz_if(a[i - 1], acc[i], m);
        b.cx(a[i], a[i - 1]);
        b.cx(a[i], acc[i]);
    }
    b.cx(carries[0], a[0]);
    let m0 = b.alloc_bit();
    b.hmr(carries[0], m0);
    b.cz_if(c_in, acc[0], m0);
    b.cx(a[0], c_in);
    b.cx(a[0], acc[0]);

    b.free_vec(&carries);
}

/// Fast `acc += a mod 2^n` using measurement-based Cuccaro.
fn add_nbit_qq_fast(b: &mut B, a: &[QubitId], acc: &[QubitId]) {
    assert_eq!(a.len(), acc.len());
    let c_in = b.alloc_qubit();
    cuccaro_add_fast(b, a, acc, c_in);
    b.free(c_in);
}

/// Fast `acc -= a mod 2^n` using measurement-based Cuccaro.
fn sub_nbit_qq_fast(b: &mut B, a: &[QubitId], acc: &[QubitId]) {
    assert_eq!(a.len(), acc.len());
    let c_in = b.alloc_qubit();
    cuccaro_sub_fast(b, a, acc, c_in);
    b.free(c_in);
}

/// `acc += a mod 2^n`. Caller must pre-extend both slices if they want the
/// top carry absorbed into the accumulator (i.e. pass n+1-bit slices with
/// top bits 0 to get a full n+1-bit add). The carry-out beyond the slice
/// is discarded via `R` on the `z` ancilla — safe when both inputs fit
/// in n-1 bits (as in our mod-p layer where both < 2p < 2^{n+1}).
fn add_nbit_qq(b: &mut B, a: &[QubitId], acc: &[QubitId]) {
    assert_eq!(a.len(), acc.len());
    let c_in = b.alloc_qubit();
    cuccaro_add(b, a, acc, c_in);
    b.free(c_in);
}

fn sub_nbit_qq(b: &mut B, a: &[QubitId], acc: &[QubitId]) {
    assert_eq!(a.len(), acc.len());
    let c_in = b.alloc_qubit();
    cuccaro_sub(b, a, acc, c_in);
    b.free(c_in);
}

fn add_nbit_const(b: &mut B, acc: &[QubitId], c: U256) {
    let n = acc.len();
    let a = load_const(b, n, c);
    add_nbit_qq(b, &a, acc);
    unload_const(b, &a, c);
}

fn sub_nbit_const(b: &mut B, acc: &[QubitId], c: U256) {
    let n = acc.len();
    let a = load_const(b, n, c);
    sub_nbit_qq(b, &a, acc);
    unload_const(b, &a, c);
}

fn csub_nbit_const(b: &mut B, acc: &[QubitId], c: U256, ctrl: QubitId) {
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
    b.free_vec(&a);
}

fn cadd_nbit_const(b: &mut B, acc: &[QubitId], c: U256, ctrl: QubitId) {
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
    b.free_vec(&a);
}

fn csub_nbit_const_fast(b: &mut B, acc: &[QubitId], c: U256, ctrl: QubitId) {
    let n = acc.len();
    let a = b.alloc_qubits(n);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, a[i]);
        }
    }
    sub_nbit_qq_fast(b, &a, acc);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, a[i]);
        }
    }
    b.free_vec(&a);
}

fn cadd_nbit_const_fast(b: &mut B, acc: &[QubitId], c: U256, ctrl: QubitId) {
    let n = acc.len();
    let a = b.alloc_qubits(n);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, a[i]);
        }
    }
    add_nbit_qq_fast(b, &a, acc);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, a[i]);
        }
    }
    b.free_vec(&a);
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
fn mod_double_inplace(b: &mut B, v: &[QubitId], p: U256) {
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
    b.free(flag);
    b.free(ovf);
}

/// Fast `v := 2*v mod p` using measurement-based Cuccaro.
fn mod_double_inplace_fast(b: &mut B, v: &[QubitId], p: U256) {
    let n = v.len();
    let ovf = b.alloc_qubit();
    b.swap(v[n - 1], ovf);
    for i in (0..n - 1).rev() { b.swap(v[i], v[i + 1]); }
    debug_assert_eq!(n, 256);
    // For secp256k1, p = 2^n - c. After the shift, the old top bit is in
    // `ovf` and the low register holds T mod 2^n for T = 2*v. If ovf=1 then
    // T = 2^n + low and T mod p = low + c; otherwise T mod p = low.
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));
    cadd_nbit_const_fast(b, v, c, ovf);
    // Result parity equals the old top bit: even if ovf=0, odd if ovf=1.
    b.cx(v[0], ovf);
    b.free(ovf);
}

/// Fast `v := v/2 mod p`. Explicit reverse of `mod_double_inplace` with
/// measurement-based Cuccaro (not emit_inverse).
fn mod_halve_inplace_fast(b: &mut B, v: &[QubitId], p: U256) {
    let n = v.len();
    let ovf = b.alloc_qubit();
    debug_assert_eq!(n, 256);
    // If v is odd, then v = low + c for some even `low`; subtract c before
    // shifting and reinsert the parity bit at the top. If v is even, this is
    // just an ordinary right shift.
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));
    b.cx(v[0], ovf);
    csub_nbit_const_fast(b, v, c, ovf);
    for i in 0..n - 1 { b.swap(v[i], v[i + 1]); }
    b.swap(v[n - 1], ovf);
    b.free(ovf);
}

/// `v := v/2 mod p`. Gate-inverse of `mod_double_inplace`.
fn mod_halve_inplace(b: &mut B, v: &[QubitId], p: U256) {
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

/// Like `cmp_lt_into` but uses carry-ancilla + measurement-based uncompute
/// for the inv_MAJ sweep. Saves n CCX. NOT emit_inverse-safe.
fn cmp_lt_into_fast(b: &mut B, u: &[QubitId], v: &[QubitId], flag: QubitId) {
    let n = u.len();
    assert_eq!(n, v.len());
    let c_in = b.alloc_qubit();
    let carries = b.alloc_qubits(n);
    for i in 0..n { b.x(u[i]); }

    // Forward MAJ sweep with carry ancillae
    b.cx(u[0], v[0]);
    b.cx(u[0], c_in);
    b.ccx(c_in, v[0], carries[0]);
    b.cx(carries[0], u[0]);
    for i in 1..n {
        b.cx(u[i], v[i]);
        b.cx(u[i], u[i - 1]);
        b.ccx(u[i - 1], v[i], carries[i]);
        b.cx(carries[i], u[i]);
    }

    b.cx(u[n - 1], flag);

    // Backward inv_MAJ with measurement
    for i in (1..n).rev() {
        b.cx(carries[i], u[i]);
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        b.cz_if(u[i - 1], v[i], m);
        b.cx(u[i], u[i - 1]);
        b.cx(u[i], v[i]);
    }
    b.cx(carries[0], u[0]);
    let m0 = b.alloc_bit();
    b.hmr(carries[0], m0);
    b.cz_if(c_in, v[0], m0);
    b.cx(u[0], c_in);
    b.cx(u[0], v[0]);

    for i in 0..n { b.x(u[i]); }
    b.free_vec(&carries);
    b.free(c_in);
}

/// Like `mod_add_qq` but uses `cmp_lt_into_fast` for the flag uncompute.
/// NOT safe inside emit_inverse blocks.
fn mod_add_qq_fast(b: &mut B, acc: &[QubitId], a: &[QubitId], p: U256) {
    let n = acc.len();
    assert_eq!(n, a.len());
    debug_assert_eq!(n, 256);

    let (acc_ext, acc_ovf) = ext_reg(b, acc);
    let (a_ext, a_ovf) = ext_reg(b, a);

    // Use fast (measurement-based) Cuccaro everywhere.
    add_nbit_qq_fast(b, &a_ext, &acc_ext);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));
    // add_nbit_const with fast Cuccaro
    {
        let n1 = acc_ext.len();
        let ca = load_const(b, n1, c);
        add_nbit_qq_fast(b, &ca, &acc_ext);
        unload_const(b, &ca, c);
    }
    let flag = b.alloc_qubit();
    b.cx(acc_ovf, flag);
    b.x(flag);
    // csub_nbit_const with fast Cuccaro
    {
        let n1 = acc_ext.len();
        let ca = b.alloc_qubits(n1);
        for i in 0..n1 { if bit(c, i) { b.cx(flag, ca[i]); } }
        sub_nbit_qq_fast(b, &ca, &acc_ext);
        for i in 0..n1 { if bit(c, i) { b.cx(flag, ca[i]); } }
        b.free_vec(&ca);
    }
    b.x(flag);
    b.cx(flag, acc_ovf);
    cmp_lt_into_fast(b, &acc_ext[..n], &a_ext[..n], flag);
    b.free(flag);

    unext_reg(b, a_ovf);
    unext_reg(b, acc_ovf);
    let _ = (acc_ext, a_ext);
}

fn cmod_add_qq(b: &mut B, acc: &[QubitId], a: &[QubitId], ctrl: QubitId, p: U256) {
    let n = acc.len();
    let f = b.alloc_qubits(n);
    for i in 0..n {
        b.ccx(ctrl, a[i], f[i]);
    }
    mod_add_qq_fast(b, acc, &f, p);
    // Gidney measurement-based AND uncomputation: f[i] = ctrl AND a[i],
    // which is unchanged by mod_add_qq (Cuccaro restores the addend).
    // HMR + classically-conditioned CZ costs 0 Toffoli vs 256 CCX.
    for i in 0..n {
        let m = b.alloc_bit();
        b.hmr(f[i], m);
        b.cz_if(ctrl, a[i], m);
    }
    b.free_vec(&f);
}

fn cmod_sub_qq(b: &mut B, acc: &[QubitId], a: &[QubitId], ctrl: QubitId, p: U256) {
    let n = acc.len();
    let f = b.alloc_qubits(n);
    for i in 0..n {
        b.ccx(ctrl, a[i], f[i]);
    }
    mod_sub_qq_fast(b, acc, &f, p);
    for i in 0..n {
        let m = b.alloc_bit();
        b.hmr(f[i], m);
        b.cz_if(ctrl, a[i], m);
    }
    b.free_vec(&f);
}

fn cmod_add_qq_bit(b: &mut B, acc: &[QubitId], a: &[QubitId], ctrl: BitId, p: U256) {
    let n = acc.len();
    let f = b.alloc_qubits(n);
    for i in 0..n {
        b.cx_if(a[i], f[i], ctrl);
    }
    mod_add_qq_fast(b, acc, &f, p);
    for i in 0..n {
        b.cx_if(a[i], f[i], ctrl);
    }
    b.free_vec(&f);
}

fn cmod_sub_qq_bit(b: &mut B, acc: &[QubitId], a: &[QubitId], ctrl: BitId, p: U256) {
    let n = acc.len();
    let f = b.alloc_qubits(n);
    for i in 0..n {
        b.cx_if(a[i], f[i], ctrl);
    }
    mod_sub_qq_fast(b, acc, &f, p);
    for i in 0..n {
        b.cx_if(a[i], f[i], ctrl);
    }
    b.free_vec(&f);
}

fn mod_mul_add_qq(
    b: &mut B,
    acc: &[QubitId],
    x: &[QubitId],
    y: &[QubitId],
    p: U256,
) {
    // acc += x * y mod p. Walk the multiplicand in place to avoid the
    // doubled tmp register and its qubit cost. For squaring, snapshot the
    // original control bits once before the in-place doubling walk.
    let n = acc.len();
    let is_squaring = x[0] == y[0];
    if is_squaring {
        let ctrl_copy = b.alloc_qubits(n);
        for i in 0..n { b.cx(x[i], ctrl_copy[i]); }
        for i in 0..n {
            cmod_add_qq(b, acc, x, ctrl_copy[i], p);
            if i < n - 1 { mod_double_inplace_fast(b, x, p); }
        }
        for _ in 0..(n - 1) { mod_halve_inplace_fast(b, x, p); }
        for i in 0..n { b.cx(x[i], ctrl_copy[i]); }
        b.free_vec(&ctrl_copy);
    } else {
        for i in 0..n {
            cmod_add_qq(b, acc, x, y[i], p);
            if i < n - 1 { mod_double_inplace_fast(b, x, p); }
        }
        for _ in 0..(n - 1) { mod_halve_inplace_fast(b, x, p); }
    }
}

/// Horner-method multiplication: acc += x * y mod p.
/// REQUIRES acc = 0 on entry. Doubles the accumulator (MSB-first),
/// avoiding the tmp register and 255 halvings entirely.
fn mod_mul_horner_add_qq(
    b: &mut B,
    acc: &[QubitId],
    x: &[QubitId],
    y: &[QubitId],
    p: U256,
) {
    let n = acc.len();
    for i in (0..n).rev() {
        if i < n - 1 { mod_double_inplace_fast(b, acc, p); }
        cmod_add_qq(b, acc, x, y[i], p);
    }
}

/// Exact inverse of `mod_mul_horner_add_qq` on the accumulator:
/// if `acc` currently holds `x * y mod p`, this maps it back to 0 while
/// leaving `x` and `y` unchanged.
fn mod_mul_horner_unadd_qq(
    b: &mut B,
    acc: &[QubitId],
    x: &[QubitId],
    y: &[QubitId],
    p: U256,
) {
    let n = acc.len();
    let is_squaring = x[0] == y[0];
    if is_squaring {
        for i in 0..n {
            cmod_sub_qq(b, acc, x, y[i], p);
            if i < n - 1 { mod_halve_inplace_fast(b, acc, p); }
        }
    } else {
        mod_neg_inplace_fast(b, x, p);
        for i in 0..n {
            cmod_add_qq(b, acc, x, y[i], p);
            if i < n - 1 { mod_halve_inplace_fast(b, acc, p); }
        }
        mod_neg_inplace_fast(b, x, p);
    }
}

/// Horner-method multiplication: acc -= x * y mod p (= acc += (p-x)*y).
/// REQUIRES acc = 0 on entry.
fn mod_mul_horner_sub_qq(
    b: &mut B,
    acc: &[QubitId],
    x: &[QubitId],
    y: &[QubitId],
    p: U256,
) {
    let n = acc.len();
    let is_squaring = x[0] == y[0];
    // Negate x, then Horner-add. For squaring: x=y, negating x also
    // negates y, giving (-x)*(-x)=x² (ADDITION, not subtraction).
    // So squaring can't use 2-neg trick.
    if is_squaring {
        mod_neg_inplace_fast(b, x, p);
        for i in (0..n).rev() {
            if i < n - 1 { mod_double_inplace_fast(b, acc, p); }
            cmod_add_qq(b, acc, x, y[i], p);
        }
        mod_neg_inplace_fast(b, x, p);
    } else {
        mod_neg_inplace_fast(b, x, p);
        for i in (0..n).rev() {
            if i < n - 1 { mod_double_inplace_fast(b, acc, p); }
            cmod_add_qq(b, acc, x, y[i], p);
        }
        mod_neg_inplace_fast(b, x, p);
    }
}

fn mod_mul_sub_qq(
    b: &mut B,
    acc: &[QubitId],
    x: &[QubitId],
    y: &[QubitId],
    p: U256,
) {
    // acc -= x * y mod p. Negate x, run schoolbook ADD (cheaper than sub),
    // then restore x. For x≠y we can walk the negated multiplicand in place
    // and halve it back afterwards, avoiding the doubled tmp register. For
    // squaring we snapshot the original control bits once into `ctrl_copy`,
    // then reuse the same in-place walk on the negated x.
    let n = acc.len();
    let is_squaring = x[0] == y[0]; // same register → squaring
    if is_squaring {
        // Hold the original x bits fixed for control while x itself walks
        // through (-x)*2^i mod p.
        let ctrl_copy = b.alloc_qubits(n);
        for i in 0..n { b.cx(x[i], ctrl_copy[i]); }
        mod_neg_inplace_fast(b, x, p);
        for i in 0..n {
            cmod_add_qq(b, acc, x, ctrl_copy[i], p);
            if i < n - 1 { mod_double_inplace_fast(b, x, p); }
        }
        for _ in 0..(n - 1) { mod_halve_inplace_fast(b, x, p); }
        mod_neg_inplace_fast(b, x, p);
        for i in 0..n { b.cx(x[i], ctrl_copy[i]); }
        b.free_vec(&ctrl_copy);
    } else {
        // Keep x negated during the loop and walk it in place.
        mod_neg_inplace_fast(b, x, p);
        for i in 0..n {
            cmod_add_qq(b, acc, x, y[i], p);
            if i < n - 1 { mod_double_inplace_fast(b, x, p); }
        }
        for _ in 0..(n - 1) { mod_halve_inplace_fast(b, x, p); }
        mod_neg_inplace_fast(b, x, p);
    }
}

fn mod_mul_add_qb(
    b: &mut B,
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
        if i < n - 1 { mod_double_inplace_fast(b, &tmp, p); }
    }
    for _ in 0..(n - 1) { mod_halve_inplace_fast(b, &tmp, p); }
    for i in 0..n { b.cx(x[i], tmp[i]); }
    b.free_vec(&tmp);
}

fn mod_mul_sub_qb(
    b: &mut B,
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
        if i < n - 1 { mod_double_inplace_fast(b, &tmp, p); }
    }
    for _ in 0..(n - 1) { mod_halve_inplace_fast(b, &tmp, p); }
    for i in 0..n { b.cx(x[i], tmp[i]); }
    b.free_vec(&tmp);
}

// ═══════════════════════════════════════════════════════════════════════════
//  Kaliski almost-inverse
// ═══════════════════════════════════════════════════════════════════════════

/// Fredkin (controlled swap): swap (a, t) if ctrl. Decomposed as CX/CCX/CX.
fn cswap(b: &mut B, ctrl: QubitId, a: QubitId, t: QubitId) {
    b.cx(t, a);
    b.ccx(ctrl, a, t);
    b.cx(t, a);
}

fn cmod_double_inplace(b: &mut B, v: &[QubitId], p: U256, ctrl: QubitId) {
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
    b.free(ovf);
}

/// `cmod_halve_inplace` = exact inverse of `cmod_double_inplace`.
fn cmod_halve_inplace(b: &mut B, v: &[QubitId], p: U256, ctrl: QubitId) {
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

    b.free(ovf);
}

/// Run `body` with `flag` holding (u < v), then uncompute the flag and
/// restore u, v. Uses carry-ancilla + measurement-based uncomputation
/// for the inv_MAJ sweep (0 Toffoli instead of n CCX).
/// Cost ≈ n CCX (forward MAJ) + body + 0 CCX (measurement inv_MAJ).
fn with_lt<F: FnOnce(&mut B)>(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    flag: QubitId,
    body: F,
) {
    let n = u.len();
    assert_eq!(n, v.len());
    let c_in = b.alloc_qubit();
    let carries = b.alloc_qubits(n);
    for i in 0..n { b.x(u[i]); }

    // Forward MAJ sweep with separate carry ancillae.
    // maj_with_carry: CX(w,y); CX(w,x); CCX(x_new,y_new,carry); CX(carry,w)
    // Step 0: (x=c_in, y=v[0], w=u[0])
    b.cx(u[0], v[0]);
    b.cx(u[0], c_in);
    b.ccx(c_in, v[0], carries[0]);
    b.cx(carries[0], u[0]);
    // Steps 1..n-1: (x=u[i-1], y=v[i], w=u[i])
    for i in 1..n {
        b.cx(u[i], v[i]);
        b.cx(u[i], u[i - 1]);
        b.ccx(u[i - 1], v[i], carries[i]);
        b.cx(carries[i], u[i]);
    }

    b.cx(u[n - 1], flag);
    body(b);
    b.cx(u[n - 1], flag);

    // Backward inv_MAJ sweep with measurement-based carry uncompute (0 Toffoli).
    // inv_maj_with_carry: CX(carry,w); HMR+CZ(carry,x,y); CX(w,x); CX(w,y)
    for i in (1..n).rev() {
        b.cx(carries[i], u[i]);             // restore w = u[i]
        let m = b.alloc_bit();
        b.hmr(carries[i], m);               // measure carry
        b.cz_if(u[i - 1], v[i], m);         // phase correction
        b.cx(u[i], u[i - 1]);               // restore x = u[i-1]
        b.cx(u[i], v[i]);                   // restore y = v[i]
    }
    // Step 0: (x=c_in, y=v[0], w=u[0])
    b.cx(carries[0], u[0]);
    let m0 = b.alloc_bit();
    b.hmr(carries[0], m0);
    b.cz_if(c_in, v[0], m0);
    b.cx(u[0], c_in);
    b.cx(u[0], v[0]);

    for i in 0..n { b.x(u[i]); }
    b.free_vec(&carries);
    b.free(c_in);
}

/// Symmetric helper: runs `body` with `flag` holding (u > v).
fn with_gt<F: FnOnce(&mut B)>(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    flag: QubitId,
    body: F,
) {
    with_lt(b, v, u, flag, body)
}

/// Run `body` with `flag` holding (v == 0), then uncompute. Single forward
/// OR chain + body + single inverse OR chain — half the cost of two
/// `cmp_eq_zero_into` calls.
fn with_eq_zero<F: FnOnce(&mut B)>(
    b: &mut B,
    v: &[QubitId],
    flag: QubitId,
    body: F,
) {
    let n = v.len();
    assert!(n > 0);
    if n == 1 {
        b.x(v[0]);
        b.cx(v[0], flag);
        body(b);
        b.cx(v[0], flag);
        b.x(v[0]);
        return;
    }
    let or_chain: Vec<QubitId> = b.alloc_qubits(n - 1);
    or_step(b, v[0], v[1], or_chain[0]);
    for i in 1..n - 1 {
        or_step(b, or_chain[i - 1], v[i + 1], or_chain[i]);
    }
    // or_chain[n-2] = (v != 0). Take complement for "== 0".
    b.x(or_chain[n - 2]);
    b.cx(or_chain[n - 2], flag);
    b.x(or_chain[n - 2]);
    body(b);
    b.x(or_chain[n - 2]);
    b.cx(or_chain[n - 2], flag);
    b.x(or_chain[n - 2]);
    for i in (1..n - 1).rev() {
        or_step(b, or_chain[i - 1], v[i + 1], or_chain[i]);
    }
    or_step(b, v[0], v[1], or_chain[0]);
    b.free_vec(&or_chain);
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
fn cmp_lt_into(b: &mut B, u: &[QubitId], v: &[QubitId], flag: QubitId) {
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

    b.free(c_in);
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
fn cmp_neq_zero_into(b: &mut B, v: &[QubitId], flag: QubitId) {
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

    b.free_vec(&or_chain);
}

/// out ^= (x OR y). `out` starts 0. Uses the de-Morgan form:
///   x(x); x(y); ccx(x, y, out); x(out); x(y); x(x);
/// After this, out = x OR y (assuming out started at 0). Its inverse is
/// the same gate sequence run in reverse — since it's symmetric (all gates
/// involutions, palindromic structure), running the exact same helper
/// again uncomputes it.
fn or_step(b: &mut B, x: QubitId, y: QubitId, out: QubitId) {
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
    b: &mut B,
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
    b: &mut B,
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
fn cmp_eq_zero_into(b: &mut B, v: &[QubitId], flag: QubitId) {
    b.x(flag);
    cmp_neq_zero_into(b, v, flag);
}

/// flag ^= (u > v).  Symmetric to cmp_lt_into(v, u, flag).
fn cmp_gt_into(b: &mut B, u: &[QubitId], v: &[QubitId], flag: QubitId) {
    cmp_lt_into(b, v, u, flag);
}

/// Controlled n-bit subtract mod 2^n: if ctrl, acc -= a. Both are n-wide
/// qubit slices. Not a mod-p operation.
fn cucc_sub_ctrl(b: &mut B, a: &[QubitId], acc: &[QubitId], ctrl: QubitId) {
    let n = a.len();
    let tmp = b.alloc_qubits(n);
    for i in 0..n { b.ccx(ctrl, a[i], tmp[i]); }
    sub_nbit_qq(b, &tmp, acc);
    for i in 0..n { b.ccx(ctrl, a[i], tmp[i]); }
    b.free_vec(&tmp);
}

/// Controlled n-bit add mod 2^n: if ctrl, acc += a.
fn cucc_add_ctrl(b: &mut B, a: &[QubitId], acc: &[QubitId], ctrl: QubitId) {
    let n = a.len();
    let tmp = b.alloc_qubits(n);
    for i in 0..n { b.ccx(ctrl, a[i], tmp[i]); }
    add_nbit_qq(b, &tmp, acc);
    for i in 0..n { b.ccx(ctrl, a[i], tmp[i]); }
    b.free_vec(&tmp);
}

/// Controlled shift-right by 1 of an n-bit register. ASSUMES v[0]=0 when
/// ctrl=1 (so no information is lost). Implemented as a controlled swap
/// cascade: if ctrl=1, new v[i] = old v[i+1] for i < n-1, new v[n-1] = 0.
fn c_shift_right_1(b: &mut B, v: &[QubitId], ctrl: QubitId) {
    let n = v.len();
    for i in 0..(n - 1) {
        cswap(b, ctrl, v[i], v[i + 1]);
    }
}

/// Unconditional shift-left by 1 of an (n+1)-bit register. ASSUMES r[n]=0
/// before the shift. After the shift: r[0]=0, r[i] = old r[i-1] for i ∈ [1, n].
fn shift_left_1(b: &mut B, r: &[QubitId]) {
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
fn shift_right_1(b: &mut B, r: &[QubitId]) {
    let n1 = r.len();
    for i in 2..n1 {
        b.swap(r[i], r[i - 1]);
    }
    b.swap(r[n1 - 1], r[0]);
}

/// flag ^= (r > c).  r is (n+1)-wide; c is a compile-time constant.
/// Non-destructive: r is restored at the end.
fn cmp_gt_const_n1(b: &mut B, r: &[QubitId], c: U256, flag: QubitId) {
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
    b: &mut B,
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
    // add_f is zero here (uncomputed in step 5 of the prior iter, or 0 at
    // init). Borrow it as the is_zero scratch instead of allocating.
    // Uses with_eq_zero_fast: measurement-based OR chain uncompute (saves 255 CCX).
    with_eq_zero_fast(b, v_w, add_f, |b| {
        b.ccx(f, add_f, m_i);
    });
    b.cx(m_i, f);

    // ─── STEP 1 ───
    //   a ^= (f=1 AND u[0]=0)
    //   m[i] ^= (f=1 AND a=0 AND v_w[0]=0)  [= f AND u[0] AND NOT v_w[0]]
    //   b ^= a; b ^= m[i]
    //
    // Shared-intermediate trick: compute z = f AND u[0] once into b_f
    // (known 0 here), then derive a_f = f XOR z = f AND NOT u[0] via CX,
    // and update m_i via ccx(z, NOT v_w[0], m_i). Uncompute z, then set
    // b_f to a_f XOR m_i as before. Saves 1 CCX per iter vs mcx2+mcx3.
    b.ccx(f, u[0], b_f);                  // b_f = f AND u[0] (z)
    b.cx(f, a_f);
    b.cx(b_f, a_f);                       // a_f = f XOR z = f AND NOT u[0]
    b.x(v_w[0]);
    b.ccx(b_f, v_w[0], m_i);              // m_i ^= z AND NOT v_w[0]
    b.x(v_w[0]);
    // Measurement-uncompute z (= f AND u[0]) from b_f: 0 CCX.
    {
        let zm = b.alloc_bit();
        b.hmr(b_f, zm);
        b.cz_if(f, u[0], zm);
    }
    b.cx(a_f, b_f);
    b.cx(m_i, b_f);                       // b_f = a_f XOR m_i

    // ─── STEP 2: with l = u > v_w: a ^= (f AND l AND ¬b); m_i ^= same.
    // Both targets share the same (f, l_gt, ¬b_f) control. Compute the AND
    // of f and l_gt once, apply to both targets, then uncompute.
    let l_gt = b.alloc_qubit();
    with_gt(b, u, v_w, l_gt, |b| {
        // add_f is 0 at this point (not yet computed in step 4). Borrow it
        // as the inner scratch for the (f AND l_gt) helper.
        b.x(b_f);                          // negate polarity of b_f
        b.ccx(f, l_gt, add_f);             // add_f = f AND l_gt
        b.ccx(add_f, b_f, a_f);            // a_f ^= add_f AND ¬b_f_orig
        b.ccx(add_f, b_f, m_i);            // m_i ^= same
        // Measurement-uncompute add_f (= f AND l_gt): 0 CCX.
        {
            let am = b.alloc_bit();
            b.hmr(add_f, am);
            b.cz_if(f, l_gt, am);
        }
        b.x(b_f);
    });
    b.free(l_gt);

    // ─── STEP 3: with control(a): swap(u, v_w); swap(r, s) ───
    for j in 0..n { cswap(b, a_f, u[j], v_w[j]); }
    for j in 0..n1 { cswap(b, a_f, r[j], s[j]); }

    // ─── STEP 4 ───
    //   add ^= (f=1 AND b=0)
    //   with control(add): v_w -= u; s += r
    //
    // Fused dual controlled sub+add: reuse one tmp register across both ops.
    // Load tmp[0..n] = add_f AND u, do sub on v_w, then transform tmp to
    // add_f AND r in place (without unloading + reloading) by temporarily
    // XOR'ing r into u and re-applying ccx(add_f, u, tmp). Then load
    // tmp[n] = add_f AND r[n], do add on s, unload. Saves n CCX/iter.
    mcx2_polar(b, f, true, b_f, false, add_f);
    {
        let tmp = b.alloc_qubits(n1);
        // Load tmp[0..n] = add_f AND u. tmp[n] stays 0.
        for i in 0..n { b.ccx(add_f, u[i], tmp[i]); }
        sub_nbit_qq_fast(b, &tmp[..n], v_w);
        // Transform tmp[0..n] from "add_f AND u" to "add_f AND r".
        for i in 0..n { b.cx(r[i], u[i]); }
        for i in 0..n { b.ccx(add_f, u[i], tmp[i]); }
        for i in 0..n { b.cx(r[i], u[i]); }
        // Cuccaro add on (n+1)-bit tmp+s with tmp[n]=0.
        add_nbit_qq_fast(b, &tmp, s);
        // XOR the missing (add_f AND r[n]) bit into s[n] directly. Saves the
        // pair of load/unload CCX that would otherwise target tmp[n].
        b.ccx(add_f, r[n], s[n]);
        // Unload tmp[0..n] via measurement-based AND uncompute (0 Toffoli).
        // tmp[i] = add_f AND r[i], unchanged by Cuccaro add (addend preserved).
        for i in 0..n {
            let m = b.alloc_bit();
            b.hmr(tmp[i], m);
            b.cz_if(add_f, r[i], m);
        }
        // tmp[n] is already 0.
        b.free_vec(&tmp);
    }

    // ─── STEP 5: uncompute add; uncompute b ───
    // Measurement-uncompute add_f = f AND (NOT b_f): 0 CCX.
    b.x(b_f);
    {
        let sm = b.alloc_bit();
        b.hmr(add_f, sm);
        b.cz_if(f, b_f, sm);
    }
    b.x(b_f);
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
        // Fast Cuccaro for add_nbit_const
        {
            let ca = load_const(b, n1, c);
            add_nbit_qq_fast(b, &ca, r);
            unload_const(b, &ca, c);
        }
        let flag = b.alloc_qubit();
        b.cx(r[n1 - 1], flag);   // flag := top bit of r' = (r ≥ p)
        // If flag=0: undo the add of c (we don't reduce).
        b.x(flag);
        // Fast Cuccaro for csub_nbit_const
        {
            let ca = b.alloc_qubits(n1);
            for i in 0..n1 { if bit(c, i) { b.cx(flag, ca[i]); } }
            sub_nbit_qq_fast(b, &ca, r);
            for i in 0..n1 { if bit(c, i) { b.cx(flag, ca[i]); } }
            b.free_vec(&ca);
        }
        b.x(flag);
        // If flag=1: clear top bit, giving r + c - 2^n = r - p.
        b.cx(flag, r[n1 - 1]);
        // Uncompute flag via parity (r was even pre-step; r-p odd if reduced).
        b.cx(r[0], flag);
        b.free(flag);
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
fn in_place_mul_const(b: &mut B, v: &[QubitId], c: U256, p: U256) {
    let n = v.len();
    let tmp = b.alloc_qubits(n);
    mul_by_const_acc(b, v, c, &tmp, p, false);       // tmp += v * c
    let c_inv = classical_modinv(c, p);
    mul_by_const_acc(b, &tmp, c_inv, v, p, true);    // v -= tmp * c_inv
    for i in 0..n { b.swap(v[i], tmp[i]); }
    b.free_vec(&tmp);
}

/// `acc ±= x * c mod p`. `c` is a classical constant. Does NOT fold acc.
/// Maintains a doubling copy of x in a temp register; adds it to acc at
/// positions where c has a bit set.
fn mul_by_const_acc(
    b: &mut B,
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
                mod_sub_qq_fast(b, acc, &tmp, p);
            } else {
                mod_add_qq_fast(b, acc, &tmp, p);
            }
        }
        if i < top {
            mod_double_inplace_fast(b, &tmp, p);
        }
    }

    // At this point tmp = x * 2^top mod p. Halve it back `top` times to
    // recover x, then uncompute via cx.
    for _ in 0..top {
        mod_halve_inplace_fast(b, &tmp, p);
    }
    for i in 0..n { b.cx(x[i], tmp[i]); }
    b.free_vec(&tmp);
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

fn alloc_kaliski_state(b: &mut B, n: usize) -> KaliskiState {
    KaliskiState {
        u: b.alloc_qubits(n),
        v_w: b.alloc_qubits(n),
        r: b.alloc_qubits(n + 1),
        s: b.alloc_qubits(n + 1),
        m_hist: b.alloc_qubits(2 * n - 1),
        f_flag: b.alloc_qubit(),
        a_flag: b.alloc_qubit(),
        b_flag: b.alloc_qubit(),
        add_flag: b.alloc_qubit(),
    }
}

fn free_kaliski_state(b: &mut B, st: KaliskiState) {
    b.free(st.add_flag);
    b.free(st.b_flag);
    b.free(st.a_flag);
    b.free(st.f_flag);
    b.free_vec(&st.m_hist);
    b.free_vec(&st.s);
    b.free_vec(&st.r);
    b.free_vec(&st.v_w);
    b.free_vec(&st.u);
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
fn kaliski_forward(b: &mut B, v_in: &[QubitId], st: &KaliskiState, p: U256) {
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
    for i in 0..(2 * n - 1) {
        kaliski_iteration(
            b, p, &st.u, &st.v_w, &st.r, &st.s,
            st.m_hist[i],
            st.f_flag, st.a_flag, st.b_flag, st.add_flag,
        );
    }

    // After the loop for nonzero v_in, classical invariants give:
    //   u = 1, v_w = 0, f = 0, a = b = add = 0
    //   r = raw coefficient (the NEGATIVE form: r = -v^{-1} * 2^{2n} mod p)
    //   s = some coefficient
    // We skip the `x(r); add_nbit_const(r, p+1)` negation (~2n CCX per call,
    // 4 calls total ≈ 8n Toffoli saved). Callers compensate by using the
    // negated inv: body multiplications that would normally `mul_add` with
    // +inv become `mul_sub` with -inv, and vice versa.
}

/// Like `with_eq_zero` but uses measurement-based uncomputation for the
/// backward OR chain (0 Toffoli instead of n-1 CCX). NOT safe inside
/// emit_inverse blocks (uses HMR ops).
fn with_eq_zero_fast<F: FnOnce(&mut B)>(
    b: &mut B,
    v: &[QubitId],
    flag: QubitId,
    body: F,
) {
    let n = v.len();
    assert!(n > 0);
    if n == 1 {
        b.x(v[0]);
        b.cx(v[0], flag);
        body(b);
        b.cx(v[0], flag);
        b.x(v[0]);
        return;
    }
    let or_chain: Vec<QubitId> = b.alloc_qubits(n - 1);
    // Forward OR chain (n-1 CCX)
    or_step(b, v[0], v[1], or_chain[0]);
    for i in 1..n - 1 {
        or_step(b, or_chain[i - 1], v[i + 1], or_chain[i]);
    }
    b.x(or_chain[n - 2]);
    b.cx(or_chain[n - 2], flag);
    b.x(or_chain[n - 2]);
    body(b);
    b.x(or_chain[n - 2]);
    b.cx(or_chain[n - 2], flag);
    b.x(or_chain[n - 2]);
    // Measurement-based uncompute (0 Toffoli)
    for i in (1..n - 1).rev() {
        or_step_uncompute(b, or_chain[i - 1], v[i + 1], or_chain[i]);
    }
    or_step_uncompute(b, v[0], v[1], or_chain[0]);
    b.free_vec(&or_chain);
}

/// Measurement-based uncompute of one or_step: uncomputes
/// `out = x OR y` using HMR + CZ (0 Toffoli).
/// Precondition: out = x OR y (was computed by or_step(x, y, out)).
/// After this: out = 0.
fn or_step_uncompute(b: &mut B, x: QubitId, y: QubitId, out: QubitId) {
    // out currently holds NOT((NOT x) AND (NOT y)) = x OR y.
    // Flip to get the AND value: (NOT x) AND (NOT y).
    b.x(out);
    // Now match the AND controls: flip x and y.
    b.x(x);
    b.x(y);
    let m = b.alloc_bit();
    b.hmr(out, m);        // measure; out → 0
    b.cz_if(x, y, m);    // phase correction with (NOT x_orig, NOT y_orig) controls
    b.x(y);
    b.x(x);
}

/// Reverse of a single kaliski_iteration. Uses measurement-based
/// uncomputation for the OR chain (with_eq_zero) and the step-4 tmp
/// unload, saving ~511 CCX per iteration vs the gate-reversed version.
fn kaliski_iteration_backward(
    b: &mut B,
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
    let n1 = r.len();

    // ── Reverse STEP 10 ─────────────────────────────────────────────────
    b.x(s[0]);
    b.cx(s[0], a_f);
    b.x(s[0]);

    // ── Reverse STEP 9 ─────────────────────────────────────────────────
    for j in (0..n1).rev() { cswap(b, a_f, r[j], s[j]); }
    for j in (0..n).rev() { cswap(b, a_f, u[j], v_w[j]); }

    // ── Reverse STEP 8 (Solinas fold) with fast Cuccaro ────────────────
    {
        let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));
        let flag = b.alloc_qubit();
        b.cx(r[0], flag);
        b.cx(flag, r[n1 - 1]);
        b.x(flag);
        // cadd_nbit_const with fast Cuccaro
        {
            let ca = b.alloc_qubits(n1);
            for i in 0..n1 { if bit(c, i) { b.cx(flag, ca[i]); } }
            add_nbit_qq_fast(b, &ca, r);
            for i in 0..n1 { if bit(c, i) { b.cx(flag, ca[i]); } }
            b.free_vec(&ca);
        }
        b.x(flag);
        b.cx(r[n1 - 1], flag);
        b.free(flag);
        // sub_nbit_const with fast Cuccaro
        {
            let ca = load_const(b, n1, c);
            sub_nbit_qq_fast(b, &ca, r);
            unload_const(b, &ca, c);
        }
    }

    // ── Reverse STEP 7 ─────────────────────────────────────────────────
    shift_right_1(b, r);

    // ── Reverse STEP 6 (conditional shift-right → shift-left) ───────────
    for i in (0..(n - 1)).rev() {
        cswap(b, f, v_w[i], v_w[i + 1]);
    }

    // ── Reverse STEP 5 ─────────────────────────────────────────────────
    b.cx(a_f, b_f);
    b.cx(m_i, b_f);
    mcx2_polar(b, f, true, b_f, false, add_f);

    // ── Reverse STEP 4 (with measurement uncompute for unload) ─────────
    {
        let tmp = b.alloc_qubits(n1);
        // Load tmp = AND(add_f, r[0..n]) for the reversed add
        for i in 0..n { b.ccx(add_f, r[i], tmp[i]); }
        // Reversed (G): CCX(add_f, r[n], s[n])
        b.ccx(add_f, r[n], s[n]);
        // Reversed (F): sub_nbit_qq_fast(tmp, s)
        sub_nbit_qq_fast(b, &tmp, s);
        // Reversed (E): transform tmp from AND(add_f,r) → AND(add_f,u)
        for i in 0..n { b.cx(r[i], u[i]); }
        for i in 0..n { b.ccx(add_f, u[i], tmp[i]); }
        for i in 0..n { b.cx(r[i], u[i]); }
        // Reversed (D): add_nbit_qq_fast(tmp[..n], v_w)
        add_nbit_qq_fast(b, &tmp[..n], v_w);
        // Reversed (C): unload AND(add_f,u) via measurement (0 Toffoli)
        for i in 0..n {
            let m = b.alloc_bit();
            b.hmr(tmp[i], m);
            b.cz_if(add_f, u[i], m);
        }
        b.free_vec(&tmp);
    }
    // Reversed (A): measurement-uncompute add_f = f AND (NOT b_f)
    b.x(b_f);
    {
        let sm = b.alloc_bit();
        b.hmr(add_f, sm);
        b.cz_if(f, b_f, sm);
    }
    b.x(b_f);

    // ── Reverse STEP 3 ─────────────────────────────────────────────────
    for j in (0..n1).rev() { cswap(b, a_f, r[j], s[j]); }
    for j in (0..n).rev() { cswap(b, a_f, u[j], v_w[j]); }

    // ── Reverse STEP 2 (with_gt body is self-inverse) ──────────────────
    let l_gt = b.alloc_qubit();
    with_gt(b, u, v_w, l_gt, |b| {
        b.x(b_f);
        b.ccx(f, l_gt, add_f);
        b.ccx(add_f, b_f, m_i);
        b.ccx(add_f, b_f, a_f);
        // Measurement-uncompute add_f = f AND l_gt: 0 CCX.
        {
            let am = b.alloc_bit();
            b.hmr(add_f, am);
            b.cz_if(f, l_gt, am);
        }
        b.x(b_f);
    });
    b.free(l_gt);

    // ── Reverse STEP 1 ─────────────────────────────────────────────────
    b.cx(m_i, b_f);
    b.cx(a_f, b_f);
    b.ccx(f, u[0], b_f);
    b.x(v_w[0]);
    b.ccx(b_f, v_w[0], m_i);
    b.x(v_w[0]);
    b.cx(b_f, a_f);
    b.cx(f, a_f);
    // Measurement-uncompute z = f AND u[0] from b_f: 0 CCX.
    {
        let zm = b.alloc_bit();
        b.hmr(b_f, zm);
        b.cz_if(f, u[0], zm);
    }

    // ── Reverse STEP 0 (with measurement uncompute of OR chain) ────────
    b.cx(m_i, f);
    {
        // Recompute OR chain (forward direction), run body, then
        // measurement-uncompute the chain.
        let nv = v_w.len();
        if nv == 1 {
            b.x(v_w[0]);
            b.cx(v_w[0], add_f);
            b.ccx(f, add_f, m_i);
            b.cx(v_w[0], add_f);
            b.x(v_w[0]);
        } else {
            let or_chain: Vec<QubitId> = b.alloc_qubits(nv - 1);
            or_step(b, v_w[0], v_w[1], or_chain[0]);
            for i in 1..nv - 1 {
                or_step(b, or_chain[i - 1], v_w[i + 1], or_chain[i]);
            }
            b.x(or_chain[nv - 2]);
            b.cx(or_chain[nv - 2], add_f);
            b.x(or_chain[nv - 2]);
            // Body
            b.ccx(f, add_f, m_i);
            // Uncompute flag
            b.x(or_chain[nv - 2]);
            b.cx(or_chain[nv - 2], add_f);
            b.x(or_chain[nv - 2]);
            // Measurement-based uncompute of OR chain (0 Toffoli)
            for i in (1..nv - 1).rev() {
                or_step_uncompute(b, or_chain[i - 1], v_w[i + 1], or_chain[i]);
            }
            or_step_uncompute(b, v_w[0], v_w[1], or_chain[0]);
            b.free_vec(&or_chain);
        }
    }
}

/// Explicit backward pass for kaliski_forward. Uses measurement-based
/// uncomputation to save ~511 CCX per iteration vs emit_inverse.
fn kaliski_backward(b: &mut B, v_in: &[QubitId], st: &KaliskiState, p: U256) {
    let n = v_in.len();

    // ─── Reverse 2n iterations (in reverse order) ───
    for i in (0..(2 * n - 1)).rev() {
        kaliski_iteration_backward(
            b, p, &st.u, &st.v_w, &st.r, &st.s,
            st.m_hist[i],
            st.f_flag, st.a_flag, st.b_flag, st.add_flag,
        );
    }

    // ─── Reverse Init ───
    b.x(st.f_flag);
    b.x(st.s[0]);
    for i in 0..n { b.cx(v_in[i], st.v_w[i]); }
    for i in 0..n { if bit(p, i) { b.x(st.u[i]); } }
}

/// Run `body` with `inv` holding `v_in^{-1} mod p`, leaving `v_in`
/// unchanged. Allocates the kaliski state and `inv` register itself, then
/// frees them at the end. The body must NOT touch `st` or `v_in`.
///
/// Implementation keeps `st` live across the body, so we only run
/// `kaliski_forward` ONCE (and its emit_inverse once), instead of the
/// 4-call structure of the previous Bennett-cleaned `kal_compute_into`.
/// Halves the dominant kaliski cost.
fn with_kal_inv_raw<F: FnOnce(&mut B, &[QubitId])>(
    b: &mut B,
    v_in: &[QubitId],
    p: U256,
    body: F,
) {
    let n = v_in.len();
    let st = alloc_kaliski_state(b, n);

    // Forward kaliski. st.r[..n] holds raw = v_in^{-1} * 2^(2n) mod p.
    kaliski_forward(b, v_in, &st, p);

    let r_low: Vec<QubitId> = st.r[..n].to_vec();
    body(b, &r_low);
    // Explicit backward pass (uses measurement-based uncompute, saves
    // ~511 CCX per iteration vs the emit_inverse version).
    kaliski_backward(b, v_in, &st, p);

    free_kaliski_state(b, st);
}

fn with_kal_inv<F: FnOnce(&mut B, &[QubitId])>(
    b: &mut B,
    v_in: &[QubitId],
    p: U256,
    body: F,
) {
    let n = v_in.len();
    with_kal_inv_raw(b, v_in, p, |b, inv_raw| {
        // Kaliski's raw output carries a 2^(2n-1) factor. Apply the
        // correction in place when callers need the exact inverse.
        for _ in 0..(2 * n - 1) { mod_halve_inplace_fast(b, inv_raw, p); }
        body(b, inv_raw);
        for _ in 0..(2 * n - 1) { mod_double_inplace_fast(b, inv_raw, p); }
    });
}

fn kaliski_inv_inplace(b: &mut B, v_in: &[QubitId], p: U256) {
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

    b.free_vec(&output);
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

pub fn build() -> Vec<Op> {
    let b = &mut B::new();
    // Register 0: target_x (quantum)
    let tx = b.alloc_qubits(N);
    b.declare_qubit_register(&tx);
    // Register 1: target_y (quantum)
    let ty = b.alloc_qubits(N);
    b.declare_qubit_register(&ty);
    // Register 2: offset_x (classical bits)
    let ox = b.alloc_bits(N);
    b.declare_bit_register(&ox);
    // Register 3: offset_y (classical bits)
    let oy = b.alloc_bits(N);
    b.declare_bit_register(&oy);

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

    // Pair 1: Kaliski's raw inverse carries a 2^(2N-1) factor. Fold that
    // scale onto lam itself, then halve lam down once. This avoids the
    // inverse-register restore pass entirely.
    with_kal_inv_raw(b, &tx, p, |b, inv_raw| {
        mod_mul_horner_add_qq(b, &lam, &ty, inv_raw, p); // lam = (-λ) * 2^(2N-1)
        for _ in 0..(2 * N - 1) { mod_halve_inplace_fast(b, &lam, p); }
        // ty = dy = λ·dx. Flip lam to +λ, run the reverse-Horner unadd,
        // then restore lam to the negative convention used below.
        mod_neg_inplace_fast(b, &lam, p);
        mod_mul_horner_unadd_qq(b, &ty, &lam, &tx, p);
        mod_neg_inplace_fast(b, &lam, p);
    });

    // Px := λ² - Px_orig - Qx. Rearranged: tx = dx - λ². Add 2Qx, then
    // negate: -(dx - λ² + 2Qx) = λ² - dx - 2Qx = Rx. mod_add_qb is
    // cheaper than mod_sub_qb (1024 vs 1280 per call, saves 512 total).
    mod_mul_sub_qq(b, &tx, &lam, &lam, p);
    mod_add_qb(b, &tx, &ox, p);
    mod_add_qb(b, &tx, &ox, p);
    mod_neg_inplace_fast(b, &tx, p);

    // With lam = -λ, this multiply produces ty = -(Ry + Qy).
    {
        let dxr = b.alloc_qubits(N);
        for i in 0..N { b.x_if(dxr[i], ox[i]); }     // dxr = Qx
        mod_sub_qq_fast(b, &dxr, &tx, p);              // dxr = Qx − Rx
        mod_mul_horner_add_qq(b, &ty, &lam, &dxr, p);  // ty += (-λ)·(Qx − Rx) = -(Ry + Qy)
        mod_add_qq_fast(b, &dxr, &tx, p);              // dxr += Rx → Qx
        for i in 0..N { b.x_if(dxr[i], ox[i]); }     // dxr = 0
        b.free_vec(&dxr);
    }
    // Pair 2 uses ty = -(Ry+Qy), which matches the negative lam convention.
    //
    // Uncompute lam. Use (Rx - Qx) as kaliski input; inv = -(Rx-Qx)^{-1}
    // = (Qx-Rx)^{-1}, so inv·ty = -λ and the same `mod_mul_sub_qq` zeroes lam.
    mod_sub_qb(b, &tx, &ox, p);                   // tx = Rx - Qx
    with_kal_inv_raw(b, &tx, p, |b, inv_raw| {
        // Scale lam up to match Kaliski's raw 2^(2N-1)-scaled inverse.
        // The product then zeroes lam directly, so no down-scaling restore
        // is needed afterward.
        for _ in 0..(2 * N - 1) { mod_double_inplace_fast(b, &lam, p); }
        mod_mul_horner_unadd_qq(b, &lam, inv_raw, &ty, p);
        mod_add_qb(b, &ty, &oy, p);                   // ty = -Ry
        mod_neg_inplace_fast(b, &ty, p);              // ty = Ry
    });
    mod_add_qb(b, &tx, &ox, p);                   // tx = Rx

    b.free_vec(&lam);

    b.ops.clone()
}
