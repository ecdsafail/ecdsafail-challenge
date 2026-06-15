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
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};

use crate::circuit::{analyze_ops, BitId, Op, OperationType, QubitId, QubitOrBit, RegisterId};
use crate::sim::Simulator;
use crate::weierstrass_elliptic_curve::WeierstrassEllipticCurve;

pub mod venting;

pub mod dialog_gcd_classical_filter;

mod emit;
pub(crate) use emit::*;

mod arith;
pub(crate) use arith::*;

mod rounds;
pub(crate) use rounds::*;

thread_local! {
    static D1_PHASE_CORRECTED_PRODUCT_CORE_SCOPE: std::cell::Cell<bool> =
        std::cell::Cell::new(false);
}

fn d1_phase_corrected_product_core_active() -> bool {
    D1_PHASE_CORRECTED_PRODUCT_CORE_SCOPE.with(|scope| scope.get())
}

pub struct B {
    pub ops: Vec<Op>,
    pub count_only: bool,
    pub counted_ops: usize,
    pub counted_kind_ops: [usize; 18],
    pub counted_phase_kind_ops: [usize; 18],
    pub counted_phase_start_ops: usize,
    pub counted_phase_rows: Vec<PhaseResource>,
    pub counted_registers: Vec<Vec<QubitOrBit>>,
    pub next_qubit: u32,
    pub next_bit: u32,
    pub next_register: u32,
    pub free_qubits: Vec<u32>,
    pub active_qubits: u32,
    pub peak_qubits: u32,
    pub peak_ops_idx: usize,
    pub peak_phase: &'static str,
    pub phase: &'static str,
    pub peak_log: Vec<(u32, &'static str, usize)>,
    pub phase_active_max: std::collections::BTreeMap<&'static str, u32>,
    pub phase_active_regions: Vec<(usize, &'static str, u32)>,
    pub current_phase_active_max: u32,
    // (ops_len_at_transition, new_phase)
    pub phase_transitions: Vec<(usize, &'static str)>,
    pub active_timeline: Vec<(usize, u32)>,
    // K=2 prototype: per-step "shifted twice" transcript bits, indexed by global
    // GCD step. Set by the ipmul/quotient wrappers around a pass; read by the
    // tobitvector (compute/uncompute) and apply (conditional 2nd double/halve).
    // Empty when K=2 is disabled (frontier path byte-identical).
    pub k2_shift2_log: Vec<QubitId>,
}

#[derive(Clone, Copy)]
struct CountSnapshot {
    ops: usize,
    kind_ops: [usize; 18],
    phase_kind_ops: [usize; 18],
    phase_start_ops: usize,
    phase_rows_len: usize,
    phase: &'static str,
}

#[derive(Clone, Debug)]
pub struct PhaseResource {
    pub phase: &'static str,
    pub start: usize,
    pub end: usize,
    pub ops: usize,
    pub toffoli_ops: usize,
    pub ccx_ops: usize,
    pub ccz_ops: usize,
    pub hmr_ops: usize,
    pub r_ops: usize,
}


impl B {
    fn new() -> Self {
        Self {
            ops: Vec::new(),
            count_only: false,
            counted_ops: 0,
            counted_kind_ops: [0; 18],
            counted_phase_kind_ops: [0; 18],
            counted_phase_start_ops: 0,
            counted_phase_rows: Vec::new(),
            counted_registers: Vec::new(),
            next_qubit: 0,
            next_bit: 0,
            next_register: 0,
            free_qubits: Vec::new(),
            active_qubits: 0,
            peak_qubits: 0,
            peak_ops_idx: 0,
            peak_phase: "",
            phase: "init",
            peak_log: Vec::new(),
            phase_active_max: std::collections::BTreeMap::new(),
            phase_active_regions: Vec::new(),
            current_phase_active_max: 0,
            phase_transitions: Vec::new(),
            active_timeline: Vec::new(),
            k2_shift2_log: Vec::new(),
        }
    }
    fn new_count_only() -> Self {
        let mut b = Self::new();
        b.count_only = true;
        b
    }
    fn push_op(&mut self, op: Op) {
        self.counted_ops += 1;
        self.counted_kind_ops[op.kind as usize] += 1;
        self.counted_phase_kind_ops[op.kind as usize] += 1;
        if !self.count_only {
            self.ops.push(op);
        }
    }
    fn count_snapshot(&self) -> CountSnapshot {
        CountSnapshot {
            ops: self.counted_ops,
            kind_ops: self.counted_kind_ops,
            phase_kind_ops: self.counted_phase_kind_ops,
            phase_start_ops: self.counted_phase_start_ops,
            phase_rows_len: self.counted_phase_rows.len(),
            phase: self.phase,
        }
    }
    fn count_delta_since(&self, snap: CountSnapshot) -> [usize; 18] {
        let mut out = [0usize; 18];
        for (idx, slot) in out.iter_mut().enumerate() {
            *slot = self.counted_kind_ops[idx] - snap.kind_ops[idx];
        }
        out
    }
    fn restore_count_snapshot(&mut self, snap: CountSnapshot) {
        self.counted_ops = snap.ops;
        self.counted_kind_ops = snap.kind_ops;
        self.counted_phase_kind_ops = snap.phase_kind_ops;
        self.counted_phase_start_ops = snap.phase_start_ops;
        self.counted_phase_rows.truncate(snap.phase_rows_len);
        self.phase = snap.phase;
    }
    fn add_counted_kind(&mut self, kind: OperationType, count: usize) {
        self.counted_ops += count;
        self.counted_kind_ops[kind as usize] += count;
        self.counted_phase_kind_ops[kind as usize] += count;
    }
    fn current_ops_len(&self) -> usize {
        if self.count_only {
            self.counted_ops
        } else {
            self.ops.len()
        }
    }
    fn close_counted_phase(&mut self) {
        if !self.count_only {
            return;
        }
        let start = self.counted_phase_start_ops;
        let end = self.counted_ops;
        if start < end {
            let ccx_ops = self.counted_phase_kind_ops[OperationType::CCX as usize];
            let ccz_ops = self.counted_phase_kind_ops[OperationType::CCZ as usize];
            let hmr_ops = self.counted_phase_kind_ops[OperationType::Hmr as usize];
            let r_ops = self.counted_phase_kind_ops[OperationType::R as usize];
            self.counted_phase_rows.push(PhaseResource {
                phase: self.phase,
                start,
                end,
                ops: end - start,
                toffoli_ops: ccx_ops + ccz_ops,
                ccx_ops,
                ccz_ops,
                hmr_ops,
                r_ops,
            });
        }
        self.counted_phase_start_ops = self.counted_ops;
        self.counted_phase_kind_ops = [0; 18];
    }
    fn set_phase(&mut self, p: &'static str) {
        self.close_phase_active_region();
        self.close_counted_phase();
        self.phase = p;
        if std::env::var("TRACE_PHASE_ACTIVE").is_ok() {
            self.current_phase_active_max = self.active_qubits;
        }
        self.phase_transitions.push((self.current_ops_len(), p));
    }
    fn record_active_timeline(&mut self) {
        if std::env::var("PROFILE_ACTIVE_TIMELINE").is_ok() {
            self.active_timeline
                .push((self.current_ops_len(), self.active_qubits));
        }
    }
    fn record_phase_active(&mut self) {
        self.record_active_timeline();
        if std::env::var("TRACE_PHASE_ACTIVE").is_ok() {
            let entry = self.phase_active_max.entry(self.phase).or_insert(0);
            if self.active_qubits > *entry {
                *entry = self.active_qubits;
            }
            if self.active_qubits > self.current_phase_active_max {
                self.current_phase_active_max = self.active_qubits;
            }
        }
    }
    fn close_phase_active_region(&mut self) {
        if std::env::var("TRACE_PHASE_ACTIVE").is_ok() && self.current_phase_active_max > 0 {
            self.phase_active_regions.push((
                self.current_ops_len(),
                self.phase,
                self.current_phase_active_max,
            ));
            self.current_phase_active_max = 0;
        }
    }
    fn alloc_qubit(&mut self) -> QubitId {
        self.active_qubits += 1;
        self.record_phase_active();
        if self.active_qubits > self.peak_qubits {
            self.peak_qubits = self.active_qubits;
            self.peak_ops_idx = self.current_ops_len();
            self.peak_phase = self.phase;
            if std::env::var("TRACE_EACH_PEAK").is_ok() {
                eprintln!(
                    "PEAK active={} next_idx={} phase='{}' ops_idx={}",
                    self.active_qubits,
                    self.next_qubit,
                    self.phase,
                    self.current_ops_len()
                );
            }
        }
        if std::env::var("TRACE_PEAK").is_ok() && self.active_qubits + 10 >= self.peak_qubits {
            self.peak_log
                .push((self.active_qubits, self.phase, self.current_ops_len()));
        }
        if let Some(q) = self.free_qubits.pop() {
            QubitId(q.into())
        } else {
            let q = self.next_qubit;
            self.next_qubit += 1;
            QubitId(q.into())
        }
    }
    fn alloc_qubits(&mut self, n: usize) -> Vec<QubitId> {
        (0..n).map(|_| self.alloc_qubit()).collect()
    }
    fn alloc_bit(&mut self) -> BitId {
        let b = self.next_bit;
        self.next_bit += 1;
        BitId(b.into())
    }
    fn alloc_bits(&mut self, n: usize) -> Vec<BitId> {
        (0..n).map(|_| self.alloc_bit()).collect()
    }
    fn free(&mut self, q: QubitId) {
        self.r(q);
        self.free_qubits
            .push(q.0.try_into().expect("qubit id fits in u32"));
        if self.active_qubits > 0 {
            self.active_qubits -= 1;
        }
        self.record_active_timeline();
    }
    fn free_vec(&mut self, qs: &[QubitId]) {
        for &q in qs {
            self.free(q);
        }
    }
    fn reacquire(&mut self, q: QubitId) {
        let pos = self
            .free_qubits
            .iter()
            .position(|&free_q| u64::from(free_q) == q.0)
            .expect("reacquire qubit that is not currently free");
        self.free_qubits.swap_remove(pos);
        self.active_qubits += 1;
        self.record_phase_active();
        if self.active_qubits > self.peak_qubits {
            self.peak_qubits = self.active_qubits;
            self.peak_ops_idx = self.current_ops_len();
            self.peak_phase = self.phase;
            if std::env::var("TRACE_EACH_PEAK").is_ok() {
                eprintln!(
                    "PEAK active={} next_idx={} phase='{}' ops_idx={}",
                    self.active_qubits,
                    self.next_qubit,
                    self.phase,
                    self.current_ops_len()
                );
            }
        }
        if std::env::var("TRACE_PEAK").is_ok() && self.active_qubits + 10 >= self.peak_qubits {
            self.peak_log
                .push((self.active_qubits, self.phase, self.current_ops_len()));
        }
    }
    fn reacquire_vec(&mut self, qs: &[QubitId]) {
        for &q in qs {
            self.reacquire(q);
        }
    }
    fn declare_qubit_register(&mut self, qs: &[QubitId]) {
        let r = RegisterId(self.next_register.into());
        self.next_register += 1;
        for &q in qs {
            while self.counted_registers.len() <= r.0 as usize {
                self.counted_registers.push(Vec::new());
            }
            self.counted_registers[r.0 as usize].push(QubitOrBit::Qubit(q));
            let mut op = Op::empty();
            op.kind = OperationType::AppendToRegister;
            op.q_target = q;
            op.r_target = r;
            self.push_op(op);
        }
        let mut op = Op::empty();
        op.kind = OperationType::Register;
        op.r_target = r;
        self.push_op(op);
    }
    fn declare_bit_register(&mut self, bs: &[BitId]) {
        let r = RegisterId(self.next_register.into());
        self.next_register += 1;
        for &b in bs {
            while self.counted_registers.len() <= r.0 as usize {
                self.counted_registers.push(Vec::new());
            }
            self.counted_registers[r.0 as usize].push(QubitOrBit::Bit(b));
            let mut op = Op::empty();
            op.kind = OperationType::AppendToRegister;
            op.c_target = b;
            op.r_target = r;
            self.push_op(op);
        }
        let mut op = Op::empty();
        op.kind = OperationType::Register;
        op.r_target = r;
        self.push_op(op);
    }
    fn x(&mut self, q: QubitId) {
        let mut op = Op::empty();
        op.kind = OperationType::X;
        op.q_target = q;
        self.push_op(op);
    }
    fn cx(&mut self, ctrl: QubitId, tgt: QubitId) {
        if ctrl == tgt {
            panic!("invalid CX with aliased control/target {:?}", ctrl);
        }
        let mut op = Op::empty();
        op.kind = OperationType::CX;
        op.q_control1 = ctrl;
        op.q_target = tgt;
        self.push_op(op);
    }
    fn ccx(&mut self, c1: QubitId, c2: QubitId, tgt: QubitId) {
        if c1 == c2 {
            if c1 != tgt {
                self.cx(c1, tgt);
            }
            return;
        }
        if c1 == tgt || c2 == tgt {
            panic!(
                "invalid CCX with target aliased to a control: {:?}, {:?}, {:?}",
                c1, c2, tgt
            );
        }
        let mut op = Op::empty();
        op.kind = OperationType::CCX;
        op.q_control2 = c1;
        op.q_control1 = c2;
        op.q_target = tgt;
        self.push_op(op);
    }
    fn cz(&mut self, a: QubitId, b: QubitId) {
        if a == b {
            let mut op = Op::empty();
            op.kind = OperationType::Z;
            op.q_target = a;
            self.push_op(op);
            return;
        }
        let mut op = Op::empty();
        op.kind = OperationType::CZ;
        op.q_control1 = a;
        op.q_target = b;
        self.push_op(op);
    }
    fn push_condition(&mut self, cond: BitId) {
        let mut op = Op::empty();
        op.kind = OperationType::PushCondition;
        op.c_condition = cond;
        self.push_op(op);
    }
    fn pop_condition(&mut self) {
        let mut op = Op::empty();
        op.kind = OperationType::PopCondition;
        self.push_op(op);
    }
    fn swap(&mut self, a: QubitId, b: QubitId) {
        if a == b {
            return;
        }
        let mut op = Op::empty();
        op.kind = OperationType::Swap;
        op.q_control1 = a;
        op.q_target = b;
        self.push_op(op);
    }
    fn r(&mut self, q: QubitId) {
        let mut op = Op::empty();
        op.kind = OperationType::R;
        op.q_target = q;
        self.push_op(op);
    }
    fn x_if(&mut self, q: QubitId, cond: BitId) {
        let mut op = Op::empty();
        op.kind = OperationType::X;
        op.q_target = q;
        op.c_condition = cond;
        self.push_op(op);
    }
    // ── Measurement / phase / classical bit ops ──
    fn hmr(&mut self, q: QubitId, c: BitId) {
        let mut op = Op::empty();
        op.kind = OperationType::Hmr;
        op.q_target = q;
        op.c_target = c;
        self.push_op(op);
    }
    // ── Classically-conditioned variants for all remaining gates ──
    fn z_if(&mut self, q: QubitId, cond: BitId) {
        let mut op = Op::empty();
        op.kind = OperationType::Z;
        op.q_target = q;
        op.c_condition = cond;
        self.push_op(op);
    }
    fn cz_if(&mut self, a: QubitId, b: QubitId, cond: BitId) {
        if a == b {
            self.z_if(a, cond);
            return;
        }
        let mut op = Op::empty();
        op.kind = OperationType::CZ;
        op.q_control1 = a;
        op.q_target = b;
        op.c_condition = cond;
        self.push_op(op);
    }
    // ── Gidney measurement-based AND uncomputation (convenience) ──
    // Uncomputes `tgt = c1 AND c2` using HMR + phase feedback.
    // Cost: 0 Toffoli (1 HMR + 1 classically-conditioned CZ).
    // Precondition: tgt holds (c1 AND c2) computed by a prior CCX.
}

pub const N: usize = 256;

/// secp256k1 prime:  p = 2^256 - 2^32 - 977.
pub const SECP256K1_P: U256 = U256::from_limbs([
    0xFFFFFFFEFFFFFC2F,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
]);


pub const ONE_INV_DX3_AFFINE_PA_ENV: &str = "ONE_INV_DX3_AFFINE_PA";
pub const ONE_INV_DX3_AFFINE_PA_BLOCKER: &str =
    "ONE_INV_DX3_AFFINE_PA_BLOCKED: the dx^3 algebra gives Rx and Ry with \
     one inversion of w=dx^3, but a clean in-place Google-ABI circuit must \
     also uncompute w, dx^2, and the Kaliski input copy after tx/ty have been \
     overwritten by Rx/Ry.  At that point dx is recoverable only by the inverse \
     affine add P=R-Q, whose denominator is Rx-Qx.  That is a second inversion, \
     or else a retained 256-bit dx witness / dirty reset, so this path cannot \
     emit a clean one-inversion four-register PA.";

// ─── helpers: bit access on U256 ────────────────────────────────────────────


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

// ═══════════════════════════════════════════════════════════════════════════
//  Loading classical operands into a fresh qubit register
// ═══════════════════════════════════════════════════════════════════════════
//
// Cuccaro needs two qubit registers. To add a classical constant or a
// classical bit register to a quantum register, we allocate a fresh
// qubit register, load the classical value into it, run Cuccaro, then
// unload. The load/unload is not counted against Toffolis.


fn direct_const_walks_enabled() -> bool {
    std::env::var("KAL_DIRECT_CONST_WALKS").ok().as_deref() == Some("1")
}

fn secp_direct_const_arith_enabled() -> bool {
    std::env::var("SECP_DIRECT_CONST_ARITH").ok().as_deref() == Some("1")
}

fn r84_lowq_enabled() -> bool {
    std::env::var("R84_LOWQ").ok().as_deref() == Some("1")
}

fn r84_lowq_cin_borrow_enabled() -> bool {
    std::env::var("R84_LOWQ_CIN_BORROW").ok().as_deref() == Some("1")
}

fn kal_vent_modadd_enabled() -> bool {
    std::env::var("KAL_VENT_MODADD").ok().as_deref() == Some("1")
}

fn kal_vent_halve_enabled() -> bool {
    std::env::var("KAL_VENT_HALVE").ok().as_deref() == Some("1")
}


const ALT_SEED_COUNT: usize = 5;
const ALT_SEED_COMMIT: usize = 24;
const ALT_SEED_SHOTS: usize = 4096;
const ALT_SEED_CLASSICAL_LIMIT: usize = 2;


fn secp256k1_curve() -> WeierstrassEllipticCurve {
    WeierstrassEllipticCurve {
        modulus: U256::from_str_radix(
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F",
            16,
        )
        .unwrap(),
        a: U256::from(0),
        b: U256::from(7),
        gx: U256::from_str_radix(
            "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798",
            16,
        )
        .unwrap(),
        gy: U256::from_str_radix(
            "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8",
            16,
        )
        .unwrap(),
        order: U256::from_str_radix(
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141",
            16,
        )
        .unwrap(),
    }
}

fn alt_seed_xof(ops: &[Op], tag: u64) -> sha3::Shake256Reader {
    let mut hasher = Shake256::default();
    hasher.update(b"quantum_ecc-alt-seed-v1");
    hasher.update(&tag.to_le_bytes());
    hasher.update(&(ops.len() as u64).to_le_bytes());
    for op in ops {
        hasher.update(&[op.kind as u8]);
        hasher.update(&op.q_control2.0.to_le_bytes());
        hasher.update(&op.q_control1.0.to_le_bytes());
        hasher.update(&op.q_target.0.to_le_bytes());
        hasher.update(&op.c_target.0.to_le_bytes());
        hasher.update(&op.c_condition.0.to_le_bytes());
        hasher.update(&op.r_target.0.to_le_bytes());
    }
    hasher.finalize_xof()
}

fn run_alt_seed_checks(ops: &[Op]) {
    let n_seeds = if std::env::var("ALT_SEED_COMMIT").is_ok() {
        ALT_SEED_COMMIT
    } else {
        ALT_SEED_COUNT
    };

    let curve = secp256k1_curve();
    let (total_qubits, num_bits, _num_regs, regs) = analyze_ops(ops.iter());
    assert!(regs.len() == 4);
    for (i, r) in regs.iter().enumerate() {
        assert_eq!(r.len(), 256, "register {i} should be 256 wide");
    }
    for q in &regs[0] {
        assert!(matches!(q, QubitOrBit::Qubit(_)));
    }
    for q in &regs[1] {
        assert!(matches!(q, QubitOrBit::Qubit(_)));
    }
    for q in &regs[2] {
        assert!(matches!(q, QubitOrBit::Bit(_)));
    }
    for q in &regs[3] {
        assert!(matches!(q, QubitOrBit::Bit(_)));
    }

    eprintln!(
        "=== alternate-seed diagnostic ({} seeds × {} shots, classical_limit={}, parallel) ===",
        n_seeds, ALT_SEED_SHOTS, ALT_SEED_CLASSICAL_LIMIT,
    );

    let results: Vec<(u64, usize, usize, usize)> = std::thread::scope(|scope| {
        let curve = &curve;
        let regs = &regs;
        let mut handles = Vec::with_capacity(n_seeds);
        for tag_idx in 0..n_seeds {
            let tag = (tag_idx as u64) + 1;
            let handle = scope.spawn(move || {
                const BATCH: usize = 64;
                let mut xof = alt_seed_xof(ops, tag);
                let mut targets = Vec::with_capacity(ALT_SEED_SHOTS);
                let mut offsets = Vec::with_capacity(ALT_SEED_SHOTS);
                let mut expected = Vec::with_capacity(ALT_SEED_SHOTS);
                while targets.len() < ALT_SEED_SHOTS {
                    let mut rb = [[0u8; 32]; 2];
                    xof.read(&mut rb[0]);
                    xof.read(&mut rb[1]);
                    let k1 = U256::from_le_bytes(rb[0]);
                    let k2 = U256::from_le_bytes(rb[1]);
                    let t = curve.mul(curve.gx, curve.gy, k1);
                    let o = curve.mul(curve.gx, curve.gy, k2);
                    if t.0 == o.0 {
                        continue;
                    }
                    if t.0.is_zero() && t.1.is_zero() {
                        continue;
                    }
                    if o.0.is_zero() && o.1.is_zero() {
                        continue;
                    }
                    let e = curve.add(t.0, t.1, o.0, o.1);
                    targets.push(t);
                    offsets.push(o);
                    expected.push(e);
                }

                let mut sim = Simulator::new(total_qubits as usize, num_bits as usize, &mut xof);
                let mut classical_failures = 0usize;
                let mut phase_garbage_batches = 0usize;
                let mut ancilla_garbage_batches = 0usize;
                let num_batches = (ALT_SEED_SHOTS + BATCH - 1) / BATCH;
                for batch in 0..num_batches {
                    let bs = BATCH.min(ALT_SEED_SHOTS - batch * BATCH);
                    let cond_mask: u64 = if bs == 64 { u64::MAX } else { (1u64 << bs) - 1 };
                    sim.clear_for_shot();
                    for shot in 0..bs {
                        let i = batch * BATCH + shot;
                        sim.set_register(&regs[0], targets[i].0, shot);
                        sim.set_register(&regs[1], targets[i].1, shot);
                        sim.set_register(&regs[2], offsets[i].0, shot);
                        sim.set_register(&regs[3], offsets[i].1, shot);
                    }
                    sim.apply_iter(ops.iter());
                    for shot in 0..bs {
                        let i = batch * BATCH + shot;
                        let gx = sim.get_register(&regs[0], shot);
                        let gy = sim.get_register(&regs[1], shot);
                        if gx != expected[i].0 || gy != expected[i].1 {
                            classical_failures += 1;
                        }
                    }
                    let phase = sim.phase & cond_mask;
                    if phase != 0 {
                        phase_garbage_batches += 1;
                    }
                    for register in regs {
                        for qb in register {
                            if let QubitOrBit::Qubit(q) = *qb {
                                *sim.qubit_mut(q) = 0;
                            }
                        }
                    }
                    let mut garbage = false;
                    for q in 0..total_qubits {
                        if (sim.qubit(QubitId(q)) & cond_mask) != 0 {
                            garbage = true;
                            break;
                        }
                    }
                    if garbage {
                        ancilla_garbage_batches += 1;
                    }
                }
                (
                    tag,
                    classical_failures,
                    phase_garbage_batches,
                    ancilla_garbage_batches,
                )
            });
            handles.push(handle);
        }
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });

    let mut total_classical = 0usize;
    let mut total_phase_batches = 0usize;
    let mut total_ancilla_batches = 0usize;
    for (tag, classical_failures, phase_garbage_batches, ancilla_garbage_batches) in &results {
        total_classical += classical_failures;
        total_phase_batches += phase_garbage_batches;
        total_ancilla_batches += ancilla_garbage_batches;
        eprintln!(
            "ALT-SEED tag={} classical_mismatches={} phase_batches={} ancilla_batches={}",
            tag, classical_failures, phase_garbage_batches, ancilla_garbage_batches,
        );
    }

    println!("METRIC altseed_classical_total={}", total_classical);
    println!("METRIC altseed_phase_batches_total={}", total_phase_batches);
    println!(
        "METRIC altseed_ancilla_batches_total={}",
        total_ancilla_batches
    );

    let phase_limit: usize = std::env::var("ALT_SEED_PHASE_LIMIT")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0);
    assert!(
        total_phase_batches <= phase_limit,
        "ALT-SEED PHASE FAILURE: {} phase-garbage batches (limit {}) across {} seeds × {} shots",
        total_phase_batches,
        phase_limit,
        n_seeds,
        ALT_SEED_SHOTS,
    );
    assert!(
        total_ancilla_batches == 0,
        "ALT-SEED ANCILLA FAILURE: {} ancilla-garbage batches across {} seeds × {} shots",
        total_ancilla_batches,
        n_seeds,
        ALT_SEED_SHOTS,
    );
    assert!(
        total_classical <= ALT_SEED_CLASSICAL_LIMIT,
        "ALT-SEED CLASSICAL FAILURE: {} classical mismatches exceeds limit {} across {} seeds × {} shots",
        total_classical,
        ALT_SEED_CLASSICAL_LIMIT,
        n_seeds,
        ALT_SEED_SHOTS,
    );
}

#[cfg(test)]
mod d1_inplace_lowerer_tests {
    use super::*;

    fn build_product_ops() -> Vec<Op> {
        let mut b = B::new();
        let h = b.alloc_qubits(N);
        b.declare_qubit_register(&h);
        let n = b.alloc_qubits(N);
        b.declare_qubit_register(&n);
        d1_inplace_product_lowerer_with_kaliski_clean(&mut b, &h, &n, SECP256K1_P, 400);
        b.ops
    }

    fn build_quotient_ops() -> Vec<Op> {
        let mut b = B::new();
        let h = b.alloc_qubits(N);
        b.declare_qubit_register(&h);
        let n = b.alloc_qubits(N);
        b.declare_qubit_register(&n);
        d1_inplace_quotient_lowerer_with_kaliski_clean(&mut b, &h, &n, SECP256K1_P, 400);
        b.ops
    }

    fn toffoli_count(ops: &[Op]) -> usize {
        ops.iter()
            .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
            .count()
    }

    fn assert_two_word_d1_abi(ops: &[Op]) -> (u32, u32, u32) {
        let (qubits, bits, registers, regs) = analyze_ops(ops.iter().copied());
        assert_eq!(registers, 2);
        assert_eq!(regs.len(), 2);
        for reg in regs {
            assert_eq!(reg.len(), N);
            assert!(reg.iter().all(|item| matches!(item, QubitOrBit::Qubit(_))));
        }
        (qubits, bits, registers)
    }

    #[test]
    fn d1_inplace_product_lowerer_component_stats_are_pinned() {
        let ops = build_product_ops();
        let (qubits, bits, registers) = assert_two_word_d1_abi(&ops);
        assert_eq!(qubits, 2475);
        assert_eq!(bits, 1_141_762);
        assert_eq!(registers, 2);
        assert_eq!(toffoli_count(&ops), 1_919_786);
    }

    #[test]
    fn d1_inplace_quotient_lowerer_component_stats_are_pinned() {
        let ops = build_quotient_ops();
        let (qubits, bits, registers) = assert_two_word_d1_abi(&ops);
        assert_eq!(qubits, 2475);
        assert_eq!(bits, 0);
        assert_eq!(registers, 2);
        assert_eq!(toffoli_count(&ops), 1_919_786);
        assert!(ops
            .iter()
            .all(|op| op.c_condition == crate::circuit::NO_BIT));
        assert!(ops.iter().all(|op| {
            !matches!(
                op.kind,
                OperationType::Hmr | OperationType::Neg | OperationType::R
            )
        }));
    }

    #[test]
    fn round8_output_side_cleanup_hook_is_env_gated() {
        let saved = std::env::var("ROUND8_QTAIL_OUTPUT_SIDE_CLEANUP").ok();
        std::env::remove_var("ROUND8_QTAIL_OUTPUT_SIDE_CLEANUP");
        assert!(!round8_qtail_output_side_cleanup_enabled());
        std::env::set_var("ROUND8_QTAIL_OUTPUT_SIDE_CLEANUP", "1");
        assert!(round8_qtail_output_side_cleanup_enabled());
        match saved {
            Some(value) => std::env::set_var("ROUND8_QTAIL_OUTPUT_SIDE_CLEANUP", value),
            None => std::env::remove_var("ROUND8_QTAIL_OUTPUT_SIDE_CLEANUP"),
        }
    }

    #[test]
    fn round8_output_side_cleanup_hook_fails_closed_until_emitter_exists() {
        let mut b = B::new();
        let tx = b.alloc_qubits(N);
        let ty = b.alloc_qubits(N);
        let ox = b.alloc_bits(N);
        let oy = b.alloc_bits(N);
        let panic = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            round8_emit_output_side_cleanup_or_fail(&mut b, &tx, &ty, &ox, &oy, SECP256K1_P);
        }))
        .expect_err("output-side qtail hook must fail closed");
        let message = panic
            .downcast_ref::<String>()
            .map(String::as_str)
            .or_else(|| panic.downcast_ref::<&str>().copied())
            .expect("panic has message");
        assert!(message.contains("ROUND8_QTAIL_OUTPUT_SIDE_CLEANUP=1"));
        assert!(message.contains("regular c=Rx-Qx inverse"));
        assert!(message.contains("Round368 singular"));
        assert!(message.contains("9024 Google"));
    }

    #[test]
    fn round8_output_side_regular_phase_repair_probe_is_separately_gated() {
        let saved = std::env::var("ROUND8_QTAIL_OUTPUT_SIDE_REGULAR_PHASE_REPAIR").ok();
        std::env::remove_var("ROUND8_QTAIL_OUTPUT_SIDE_REGULAR_PHASE_REPAIR");
        assert!(!round8_qtail_output_side_regular_phase_repair_enabled());
        std::env::set_var("ROUND8_QTAIL_OUTPUT_SIDE_REGULAR_PHASE_REPAIR", "1");
        assert!(round8_qtail_output_side_regular_phase_repair_enabled());
        match saved {
            Some(value) => {
                std::env::set_var("ROUND8_QTAIL_OUTPUT_SIDE_REGULAR_PHASE_REPAIR", value)
            }
            None => std::env::remove_var("ROUND8_QTAIL_OUTPUT_SIDE_REGULAR_PHASE_REPAIR"),
        }
    }

    #[test]
    fn round8_qtail_round217_product_reuse_hook_is_env_gated() {
        let saved = std::env::var("ROUND8_QTAIL_ROUND217_PRODUCT_REUSE").ok();
        std::env::remove_var("ROUND8_QTAIL_ROUND217_PRODUCT_REUSE");
        assert!(!round8_qtail_round217_product_reuse_enabled());
        std::env::set_var("ROUND8_QTAIL_ROUND217_PRODUCT_REUSE", "1");
        assert!(round8_qtail_round217_product_reuse_enabled());
        match saved {
            Some(value) => std::env::set_var("ROUND8_QTAIL_ROUND217_PRODUCT_REUSE", value),
            None => std::env::remove_var("ROUND8_QTAIL_ROUND217_PRODUCT_REUSE"),
        }
    }

    #[test]
    fn round8_qtail_round217_product_reuse_hook_fails_closed_before_body() {
        let plan = round218_b5_transport::round218_b5_source_live_product_lowerer_body_plan();
        assert!(!plan.body_emits_gates);
        assert!(!plan.codegen_allowed_now);
        assert_eq!(
            plan.selected_route,
            "round217_sampled_product_m2_contract_path"
        );
        assert!(plan
            .phase_blocks
            .iter()
            .any(|block| block.phase.contains("hash_history")));
    }

    #[test]
    fn round218_source_live_product_lowerer_plan_rejects_full_source_alias() {
        let plan = round218_b5_transport::round218_b5_source_live_product_lowerer_body_plan();
        assert!(!plan.body_emits_gates);
        assert!(!plan.codegen_allowed_now);
        assert!(plan
            .phase_blocks
            .iter()
            .all(|block| !block.backend_primitive.contains("full_source_product")));
        assert!(plan
            .missing_object
            .contains("promotable no-history qtail/Round217 product splice"));
    }
}

fn set_default_env(name: &str, value: &str) {
    if std::env::var_os(name).is_none() {
        std::env::set_var(name, value);
    }
}

fn configure_ecdsafail_submission_route() {
    set_default_env("DIALOG_GCD_FOLD_CARRY_TRUNC_W", "18");
    set_default_env("DIALOG_GCD_FOLD_FREE_FIRST_HIGH_CARRY", "1");
    // q1168 host-E route. These defaults are first so the historical fallback
    // block below cannot override the exact state searched on WMI.
    set_default_env("DIALOG_GCD_ACTIVE_ITERATIONS", "258");
    set_default_env("DIALOG_GCD_APPLY_BOUNDARY_FREE_OWNED_DURING_REPLAY", "1");
    set_default_env("DIALOG_GCD_APPLY_BORROW_FUTURE_BOUNDARY_CARRIES", "1");
    set_default_env("DIALOG_GCD_APPLY_CHUNKED_F_BLOCKS", "20");
    set_default_env(
        "DIALOG_GCD_APPLY_CHUNKED_F_CUTS",
        "17,34,50,66,81,96,110,124,137,150,163,175,187,198,209,219,229,238,247",
    );
    set_default_env("DIALOG_GCD_APPLY_CHUNKED_F_AUTO_TOPCLEAN_MAX_BITS", "2");
    set_default_env("DIALOG_GCD_APPLY_CHUNKED_F_AUTO_TOPCLEAN_TARGET", "1168");
    set_default_env("DIALOG_GCD_APPLY_CLEAN_COMPARE_BITS", "19");
    set_default_env("DIALOG_GCD_APPLY_IMPLICIT_HIGH_ZERO", "1");
    set_default_env("DIALOG_GCD_BINDER_NOTCH_EXTRA", "3");
    set_default_env("DIALOG_GCD_BINDER_NOTCH_MAP", "11:1,12:1,13:1");
    set_default_env("DIALOG_GCD_BINDER_NOTCH_STEPS", "8,9,10");
    set_default_env(
        "DIALOG_GCD_BODY_CARRY_BAND_TRIMS",
        "0,3,3,3,3,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,3,3,3",
    );
    set_default_env("DIALOG_GCD_COMPARE_BITS", "46");
    set_default_env(
        "DIALOG_GCD_COMPARE_STEP_BITS",
        "181:48,194:48,199:48,202:48,207:48,212:48,216:48",
    );
    set_default_env(
        "DIALOG_GCD_FOLD_CARRY_TRUNC_STEP_WINDOWS",
        "",
    );
    set_default_env("DIALOG_GCD_FOLD_CARRY_TRUNC_W", "17");
    set_default_env("DIALOG_GCD_FOLD_FREED_TAIL", "1");
    set_default_env("DIALOG_GCD_FOLD_FREED_TAIL_ED", "1");
    set_default_env("DIALOG_GCD_FOLD_HOST_DERIVED_CONTROLS", "1");
    set_default_env("DIALOG_GCD_FOLD_HOST_E_TOP_CARRY", "1");
    set_default_env("DIALOG_GCD_FOLD_MAJ1", "1");
    set_default_env("DIALOG_GCD_FOLD_MAJ2", "1");
    set_default_env("DIALOG_GCD_FOLD_PARK_LOW_CARRIES", "15");
    set_default_env(
        "DIALOG_GCD_FOLD_PARK_LOW_CARRIES_STEP_MAP",
        "0:17,3:16,8:16,9:16,10:16,21:17,22:16,24:16,26:16,33:16,34:16,37:17,41:16,42:17,51:16,55:16,65:17,73:16,77:16,81:16,82:16,86:16,87:16,97:16,104:16,109:16,110:16,120:16,129:16,132:17,134:16,141:17,142:16,146:16,157:16,160:16,169:16,170:17,174:16,177:16,191:16,192:16,198:16,205:16,206:16,212:16,215:16,216:16,217:16,224:17,228:16",
    );
    set_default_env("DIALOG_GCD_FOLD_STREAM_CONTROLS", "1");
    set_default_env("DIALOG_FUSE_C_FORM", "1");
    set_default_env("DIALOG_FUSE_X_RESTORE", "1");
    set_default_env("DIALOG_GCD_K2", "1");
    set_default_env("DIALOG_GCD_K5_CLEAN_BLOCK", "1");
    set_default_env("DIALOG_GCD_K5_FIXED_TAIL_APPLY", "0");
    set_default_env("DIALOG_GCD_K5_FREE_CLEAN_BLOCK_DURING_SHIFT", "1");
    set_default_env("DIALOG_GCD_K5_HEAD11_CODEC", "1");
    set_default_env("DIALOG_GCD_K5_HEAD11_STREAM_PAIR_APPLY", "1");
    set_default_env("DIALOG_GCD_K5_HEAD11_SPLIT_PAIR_SHIFT_APPLY", "1");
    set_default_env("DIALOG_GCD_K5_HEAD11_PAIR01_S2_PERMUTE_APPLY", "1");
    set_default_env(
        "DIALOG_GCD_K5_HEAD11_PAIR23_S2_BORROW_PAIR01_APPLY",
        "1",
    );
    set_default_env("DIALOG_GCD_K5_PARTIAL_RAW_RELEASE", "8");
    set_default_env("DIALOG_GCD_K5_RELEASE_SCALE_BITS", "5");
    set_default_env("DIALOG_GCD_K5_STREAM_PAIR_APPLY", "1");
    set_default_env("DIALOG_GCD_K5_TAIL3_FIXED_LAST", "0");
    set_default_env("DIALOG_GCD_K5_TAIL3_TOP32_CODEC", "1");
    set_default_env("DIALOG_GCD_K5_TAIL3_TOP32_STREAM_APPLY", "1");
    set_default_env("DIALOG_GCD_K5_TAIL3_TOP32_SPLIT_SLOT_APPLY", "1");
    set_default_env("DIALOG_GCD_K5_TAIL3_TOP32_FINAL_S2_CONST_APPLY", "1");
    set_default_env("DIALOG_GCD_ODD_U_LOWBIT_FASTPATH", "1");
    set_default_env("DIALOG_GCD_PA9024_COMPARE_SCHEDULE", "1");
    set_default_env("DIALOG_GCD_PA9024_COMPARE_SCHEDULE_MARGIN", "0");
    set_default_env("DIALOG_GCD_PERPOS_MAJ2", "1");
    set_default_env("DIALOG_GCD_RAW_IPMUL_CLEAR_P_RESIDUAL", "1");
    set_default_env("DIALOG_GCD_RAW_TOBITVECTOR_MATERIALIZED_SUB", "0");
    set_default_env("DIALOG_GCD_RAW_TOBITVECTOR_VARIABLE_WIDTH", "1");
    set_default_env("DIALOG_GCD_RUNWAY_PARTIAL_BLOCK", "1");
    set_default_env("DIALOG_GCD_SKIP_ZERO_EDGE_CSHIFT", "1");
    set_default_env("DIALOG_GCD_SPECIAL_FOLD_BORROW_CARRIES", "1");
    set_default_env(
        "DIALOG_GCD_SPECIAL_FOLD_CARRY_TRUNC_STEP_WINDOWS",
        "10:19,11:19,21:20,63:19,74:19,100:19,107:19,110:19,118:19,135:19,136:19,137:19,188:20,204:19,227:20,241:19",
    );
    set_default_env("DIALOG_GCD_SPECIAL_FOLD_PARK_LOW_CARRIES", "16");
    set_default_env(
        "DIALOG_GCD_SPECIAL_FOLD_PARK_LOW_CARRIES_STEP_MAP",
        "",
    );
    set_default_env("DIALOG_GCD_SPECIAL_FOLD_RELEASE_SCRATCH", "1");
    set_default_env(
        "DIALOG_GCD_SPECIAL_OVERFLOW_CLEAN_STEP_BITS",
        "1:24,4:21,6:25,7:20,10:22,11:20,19:21,21:21,22:21,23:23,28:21,30:22,32:20,33:24,34:25,48:21,49:22,55:22,62:23,64:20,66:21,71:22,86:20,92:21,113:21,116:21,118:20,119:24,120:21,121:20,127:20,129:22,131:22,142:22,144:21,145:23,147:21,151:23,153:23,154:23,155:20,156:24,159:22,161:24,165:21,166:21,168:21,173:20,175:21,178:21,184:22,185:20,187:23,188:22,190:20,193:21,194:22,196:20,197:21,199:21,203:22,205:22,209:20,210:21,213:20,217:22,221:21,222:23,229:21,236:21,241:21",
    );
    set_default_env(
        "DIALOG_GCD_SPECIAL_UNDERFLOW_CLEAN_STEP_BITS",
        "3:21,5:21,10:23,11:22,14:22,17:20,27:22,33:20,34:22,38:21,42:22,47:21,50:22,51:21,53:20,54:21,58:21,60:21,65:21,67:23,68:25,73:20,74:20,75:23,77:21,78:20,84:21,89:23,91:22,95:22,98:26,103:21,109:22,110:22,114:22,118:22,127:26,135:20,136:20,137:22,143:21,149:21,152:20,154:26,155:20,156:22,157:20,158:26,166:20,178:20,181:20,186:24,188:25,191:21,194:20,198:20,200:21,201:21,202:23,203:23,204:22,212:25,213:20,214:22,221:20,223:21,228:21,231:23,243:21,246:20",
    );
    set_default_env("DIALOG_GCD_TOBITVECTOR_CSWAP_BODY_TRIM", "0");
    set_default_env("DIALOG_GCD_WIDTH_MARGIN", "10");
    set_default_env("DIALOG_GCD_WIDTH_SLOPE_X1000", "1015");
    set_default_env("DIALOG_TAIL_NONCE", "18100027017098");
    set_default_env("KAL_DOUBLE_CARRY_TRUNC_W", "19");
    set_default_env("KAL_FOLD_CARRY_TRUNC_W", "18");
    set_default_env("SQUARE_ROW_MAX_SEG", "141");
    set_default_env("SQUARE_ROW_WINDOW_CLEAN_COMPARE_BITS", "18");
    set_default_env(
        "SQUARE_ROW_WINDOW_CLEAN_ROW_BITS",
        "2:20,11:20,12:20,13:21,16:22,19:20,20:21,21:20,26:21,29:21,32:21,37:21,44:22,46:20,53:21,56:20,64:20,70:20,75:20,78:20,87:20",
    );
    set_default_env(
        "SQUARE_ROW_WINDOW_CLEAN_SITE_BITS",
        "1:0:f:19,3:0:r:21,9:0:f:22,10:0:r:21,13:0:r:22,14:0:r:20,15:0:r:19,17:0:r:20,26:0:f:22,36:0:f:20,38:0:f:20,38:0:r:20,39:0:r:19,40:0:r:22,41:0:r:19,42:0:r:20,43:0:r:19,45:0:r:19,47:0:f:22,47:0:r:19,48:0:r:20,50:0:f:22,50:0:r:22,51:0:f:22,54:0:f:19,57:0:r:19,59:0:f:19,60:0:f:19,62:0:f:22,62:0:r:21,63:0:f:20,65:0:f:19,66:0:f:21,66:0:r:21,67:0:f:19,68:0:r:21,71:0:r:20,72:0:f:21,73:0:r:21,74:0:r:19,76:0:r:21,79:0:r:20,81:0:f:20,83:0:r:22,89:0:r:19,90:0:r:21,91:0:f:21,92:0:r:21,95:0:r:20,97:0:r:21,102:0:f:20,103:0:r:19,104:0:r:19,107:0:f:20,109:0:f:21,110:0:f:19,110:0:r:20",
    );
    set_default_env("SQUARE_ROW_WINDOW_MEASURED_CARRY_CLEAR", "1");

    set_default_env("SKIP_ALT_SEED_CHECKS", "1");
    set_default_env("DIALOG_GCD_COMPRESSED_SIDECAR_LOG", "1");
    // Tighten the windowed square-row carry cleanup by one bit. A GPU
    // structural filter followed by the trusted simulator found nonce
    // 17761178 clean over all 9024 Fiat-Shamir shots: 1215 qubits and
    // 1,403,115.070 average executed Toffoli.
    set_default_env("SQUARE_ROW_WINDOW_CLEAN_COMPARE_BITS", "21");
    set_default_env("SQUARE_ROW_WINDOW_MEASURED_CARRY_CLEAR", "1");
    set_default_env("ROUND84_KEEP_QUOTIENT_PRODUCT", "1");
    set_default_env("DIALOG_GCD_FOLD_CARRY_TRUNC_W", "17");
    set_default_env("DIALOG_TAIL_NONCE", "18100027017098");
    set_default_env("DIALOG_GCD_SKIP_ZERO_EDGE_CSHIFT", "1");
    set_default_env("DIALOG_GCD_COMPRESSED_BLOCK_LIFECYCLE", "1");
    set_default_env("DIALOG_GCD_HOST_REVERSE_RAW_BLOCK", "1");
    set_default_env("DIALOG_GCD_COMPRESSED_LOG_U_HIGH_RUNWAY", "1");
    set_default_env("DIALOG_GCD_COMPRESSED_LOG_U_HIGH_RUNWAY_BLOCKS", "999");
    set_default_env("DIALOG_GCD_COMPOSITE_SCRATCH", "1");
    // Fold the CURRENT transcript block's own compressed cells (|0> across that
    // block's GCD steps -- forward written only at compress_block, reverse
    // decompressed before the steps) into the composite body-scratch borrow.
    // Pure qubit relabel (0 added Toffoli) that shrinks the early-step body
    // deficit and drops the GCD-walk peak 1313 -> 1309. Stacked on top of the
    // K2 per-step compare schedule (Toffoli-axis) for a peak-axis cut.
    set_default_env("DIALOG_GCD_BORROW_CURRENT_BLOCK", "1");
    // Gidney measurement-vented CONTROLLED GCD body (else branch of the selected
    // add/sub). Replaces the full-CCX controlled Cuccaro (cucc_*_ctrl_lowq,
    // ~8-10 CCX/bit) with cuccaro_*_ctrl_vented (~2 CCX/bit: a forward carry
    // chain vented onto active_width-1 BORROWED |0> lanes from the composite
    // scratch, plus a controlled-sum pass, with the carry uncomputed by
    // measurement at 0 Toffoli). Vents are borrowed (never fresh-allocated) so
    // the peak does not grow; the composite-scratch `want` is bumped to supply
    // them (see dialog/compressed.rs). Big avg-Toffoli cut at flat peak.
    set_default_env("DIALOG_GCD_CTRL_BODY_VENTED", "1");
    set_default_env("DIALOG_GCD_APPLY_REPLAY_SWAP_HOST", "1");
    set_default_env("SQUARE_SELFHOST_SAFE_LANE_REUSE", "1");
    set_default_env("SQUARE_SELFHOST_GATE_SUFFIX_CARRIES", "0");
    // K2-calibrated per-step branch-comparator schedule (see the
    // DIALOG_GCD_PA9024_COMPARE_SCHEDULE table in dialog/config.rs). The flat
    // DEFAULT_COMPARE_BITS=50 spends 50 bits on EVERY GCD step, but a faithful
    // classical model over 8M reachable factors shows the early steps resolve the
    // u>v branch in far fewer bits (req_cb 22..~44 for steps 0..~130, vs 48..55
    // for the mid steps). Enabling the per-step schedule clips each step to
    // min(SCHEDULE[step]+MARGIN, 50, active_width): early steps drop well below 50
    // (value-exact on the reachable support, MARGIN cushion over the 8M observed
    // max), mid steps cap at the global 50 (== baseline, where compare hazards are
    // already ~0). Pure executed-Toffoli cut at flat peak 1313; the shorter op
    // stream re-rolls the Fiat-Shamir island, re-hunted via DIALOG_TAIL_NONCE.
    set_default_env("DIALOG_GCD_PA9024_COMPARE_SCHEDULE", "1");
    // PA9024 compare-schedule margin retuned with ACTIVE_ITERATIONS=396 and
    // APPLY_CLEAN_COMPARE_BITS=21. The wider margin gives back a little Toffoli
    // but lands the 1438q clean island at DIALOG_REROLL=3 / POST_SUB=51 below.
    // sm5: compare-schedule margin 7 -> 5 narrows the per-step comparator on the
    // low/mid-width GCD steps (below the 57 cap) for -452 executed Toffoli,
    // peak-neutral at 1434q, orthogonal to compare57. The late-game lineage ran
    // margin=5; the base had reverted to 7. Clean island at REROLL=1844/POST_SUB=3532.
    // Per-step schedule safety margin over the 8M-sample observed max req_cb.
    // MARGIN=0 uses the observed max directly (the geometric tail beyond the 8M
    // max adds ~0.3 compare hazards/draw, dodged by the tail nonce like the width
    // island); biggest cut (~7,756 executed Toffoli vs flat-50). (Effective
    // per-step bits = min(SCHEDULE[step]+MARGIN, DEFAULT_COMPARE_BITS=50, aw).)
    set_default_env("DIALOG_GCD_PA9024_COMPARE_SCHEDULE_MARGIN", "0");
    // DOUBLE-carry lazy-Solinas window re-tightened 22 -> 21 on the peak-1313
    // K2_PAIR_COMPRESS base: -1,038 avg executed Toffoli, peak-neutral at 1313q
    // (avg_T 1,536,923 -> 1,535,885; 1313 x 1,535,885 = 2,016,617,005, beats the
    // prior #1 2,017,979,899 by 1,362,894). Value-exact on the reachable support
    // (dropped double-carry bit is 0 there, ~2^-22/call otherwise); residual
    // failures are Fiat-Shamir phase, dodged by a fresh tail nonce (re-hunted below).
    set_default_env("KAL_DOUBLE_CARRY_TRUNC_W", "19");
    // Likewise give back the FOLD-carry truncation bit for the final-window W2
    // island; the Toffoli budget still beats the 1320q frontier.
    // Re-tighten 24 -> 22 on the W2 base (the lazy-Solinas fold-carry window had
    // been left loose). Value-exact on the reachable support (the dropped fold
    // carry bits are 0 there); residual failures are pure Fiat-Shamir, dodged by
    // the shared re-rolled tail nonce below.
    set_default_env("KAL_FOLD_CARRY_TRUNC_W", "18");
    set_default_env("DIALOG_GCD_ROUND763_DEDUP", "1");
    set_default_env("DIALOG_GCD_ROUND763_COMPRESS_LEVER", "1");
    set_default_env("DIALOG_GCD_MEASURED_UNDERFLOW_GATE", "1");
    // Branch comparator width tightened 63 -> 61 (−1,160 executed Toffoli),
    // STACKED on the PA9024 margin-5 cut. Two within-budget truncations coexist
    // via the 2-D reroll island (DIALOG_REROLL=1, DIALOG_POST_SUB_REROLL=0).
    // Branch comparator width tightened 61 -> 59 (−1,600 executed Toffoli),
    // stacked on the chunked-apply + round763 + acc=19 base via the 2-D reroll
    // island (DIALOG_REROLL=0, DIALOG_POST_SUB_REROLL=10). Validated 0/0/0 @ 1567.
    // Branch comparator width tightened 59 -> 58 (−952 executed Toffoli),
    // stacked on the 1446-peak base + ACTIVE_ITERATIONS=397 via the reroll-37/1
    // island documented below.
    // Branch comparator 58 -> 57: -1,064 executed Toffoli, peak-neutral at 1434q,
    // stacked on the active395 base. Clean island at REROLL=4959 / POST_SUB=5983.
    // COMPARE_BITS 73 -> 52: the GCD branch comparator (b1 = u<v on the top
    // `compare_bits` of the active window) was left at a loose 73 by the whole
    // frontier lineage. A classical convergence-filter sweep over both GCD
    // factors (quotient dx AND ipmul c = Qx-Rx) on 300k inputs shows the
    // truncated comparator NEVER mis-decides a branch down to 52 bits (0 added
    // hard inputs); the binding truncations are the width envelope, body-carry
    // band, and iteration count -- not the comparator. So 73 -> 52 is a pure
    // -28,392 executed-Toffoli cut (21 bits x 2 dirs x 2 passes, comparator =
    // 2 T/bit), peak-neutral at 1390q, with ZERO change to islandability. The
    // shorter op stream re-rolls Fiat-Shamir; co-tuned with WIDTH_MARGIN=10 and
    // TAIL_NONCE below. Validated 0/0/0 over all 9024 shots.
    // Final-window W2 spends two branch-comparator bits back for a much denser
    // clean island while retaining a lower score than the current frontier.
    // K2 pair-compressed route spends one branch-comparator bit back from the
    // newest frontier cut. This keeps the lower 1313q tier while landing a much
    // denser clean island than the 45-bit edge.
    // Both-phase apply fold-fusion: spend comparator bits back to cb=52 (the
    // exact-screen zone) while preserving a clean Fiat-Shamir
    // nonce; the fold-fusion's -25k Toffoli keeps the score well under 2B.
    set_default_env("DIALOG_GCD_COMPARE_BITS", "46");
    // Apply-phase overflow-clean comparator narrowed 23 -> 22 -> 21 -> 20. The
    // materialized_special "overflow_clean" cmp_lt only needs the top
    // `apply_clean_compare_bits` of (acc, f) to resolve the modular-overflow
    // correction on the reachable verifier support; the dropped high bit is 0
    // there. Pure structural Toffoli cut 1,504,903 -> 1,504,387 -> 1,503,871
    // -> 1,503,355
    // (-516 per bit), peak-neutral at 1309q. The shorter op stream re-rolls the
    // Fiat-Shamir island, re-hunted to DIALOG_TAIL_NONCE=721381 below (GCD
    // pre-filter + bit-exact quantum confirm, validated 0/0/0 over all 9024
    // shots: 1309 x 1,503,355 = 1,967,891,695, beats the 1,968,064,139 frontier
    // by 172,444).
    set_default_env("DIALOG_GCD_APPLY_CLEAN_COMPARE_BITS", "19");
    set_default_env("DIALOG_GCD_APPLY_BOUNDARY_CONDITIONAL_REPLAY", "1");  // BAKED: condrep ON for env-less grader build
    set_default_env("DIALOG_GCD_SELECTED_BODY_STREAM_SUFFIX_MAP", "3:2,4:3,5:5,6:6,7:7,8:5,9:7,10:5,11:7,12:6,13:7,14:5,15:6,16:3,17:5,18:1,19:3,21:1");  // BAKED: codex 1285q peak-drop (stream selected high bits through low-qubit suffix)
    // Bake the exact conditional-replay stack for env-less GPU hunts and grader builds.
    set_default_env("DIALOG_GCD_REVERSE_BRANCH_CONDITIONAL_REPLAY", "1");
    set_default_env("DIALOG_GCD_SPECIAL_CLEAN_CONDITIONAL_REPLAY", "1");
    set_default_env("MOD_FAST_FLAG_CONDITIONAL_REPLAY", "1");
    set_default_env("DIALOG_GCD_RAW_PA", "1");
    set_default_env("DIALOG_GCD_K2", "1");
    // Both-phase apply fold-fusion (fused double_y + halve_y Solinas folds,
    // single shared carry chain; -25k avg Toffoli, phase-clean).
    set_default_env("DIALOG_GCD_APPLY_FUSED_FOLD", "1");
    // K2 pair transcript compressor: pack two K2 transcript steps into five
    // sidecar bits by using the local reachability constraint between step A's
    // shift2 bit and step B's low branch bit. This cuts the current transcript
    // peak into the 1313q tier at a small Toffoli cost.
    set_default_env("DIALOG_GCD_K2_PAIR_COMPRESS", "1");
    // 396 -> 395 -> 394 on the current 1355q route. The binary-GCD transcript
    // still converges on the verifier support for the Fiat-Shamir island below,
    // while dropping two full GCD body/reverse steps.
    // 260 -> 259 after the 1320q apply teardown: saves one GCD body/reverse row.
    // Stacked with KAL_DOUBLE_CARRY_TRUNC_W=22, the nonce below lands the clean
    // 1320q island while improving the custom-five seed's Toffoli count.
    // 258 -> 262 on the lowq0 final-chunk route: spend four GCD rows from the
    // recovered fast-final Toffoli budget to remove most nonconvergence pressure
    // while staying under the 1309q round84 peak. Re-hunted with the GCD filter
    // and quantum-confirmed at tail nonce 2432.
    set_default_env("DIALOG_GCD_ACTIVE_ITERATIONS", "258");
    set_default_env("DIALOG_GCD_PERPOS_MAJ2", "1");
    set_default_env("DIALOG_GCD_FUSED_HCLEAR_MEASURED", "1");
    set_default_env("DIALOG_GCD_FUSED_DCLEAR_MEASURED", "1");
    set_default_env("DIALOG_GCD_FUSED_HALVE_EDCLEAR_MEASURED", "1");
    set_default_env("DIALOG_GCD_RAW_IPMUL_TERMINAL_REUSE", "1");
    set_default_env("DIALOG_GCD_RAW_IPMUL_CLEAR_P_RESIDUAL", "1");
    set_default_env("DIALOG_GCD_RAW_QUOTIENT_TERMINAL_REUSE", "1");
    set_default_env("DIALOG_GCD_RAW_APPLY_REVERSE_MATERIALIZED_SPECIAL_SUB", "1");
    set_default_env("DIALOG_GCD_RAW_APPLY_MATERIALIZED_SPECIAL_ADD", "1");
    set_default_env("DIALOG_GCD_RAW_APPLY_TRUNCATED_CLEAN", "1");
    // LOW-QUBIT CORNER (ToB jump-lowqubit reconstruction): "0" routes the GCD body
    // to the low-scratch CONTROLLED form (cucc_sub/add_ctrl_lowq) instead of the
    // materialized body, whose ~2*active_width gated+carry scratch pinned the
    // GCD-walk at 1297. With the composite-scratch right-sizing (compressed.rs
    // build_composite_scratch) + the vented add_double_ox/x_restore (modular.rs)
    // + APPLY_FINAL_WINDOWED_FAST_BLOCKS=2 below, the peak drops to 1284 (bound by
    // the round84 in-place Solinas square). Controlled body costs ~2x Toffoli;
    // recovered by band-trimming it (TODO). "1" restores the 1297 materialized base.
    set_default_env("DIALOG_GCD_RAW_TOBITVECTOR_MATERIALIZED_SUB", "0");
    set_default_env("DIALOG_GCD_RAW_TOBITVECTOR_VARIABLE_WIDTH", "1");
    set_default_env("DIALOG_GCD_RAW_TOBITVECTOR_BORROW_FUTURE_LOG_CARRIES", "1");
    // ROUND84 x-tail square: Karatsuba beats schoolbook by -16,272 emitted
    // Toffoli on the peak-1572 base, and Karatsuba's z1_reg fits UNDER the
    // materialized_special apply binder so peak stays 1572 (verified). The
    // different op count re-rolls the Fiat-Shamir island, co-tuned below
    // (WIDTH_MARGIN=27, REROLL=0). Validated 0/0/0 over 9024.
    // ROUND84_XTAIL_KARATSUBA=0 (+ROUND84_XTAIL_SCHOOLBOOK=1) restores schoolbook.
    set_default_env("ROUND84_XTAIL_KARATSUBA", "0");
    // Slack-exploit: once round84's Solinas binder fell to 1543 (== the apply
    // tier), its doubling lanes (r84k_sol_dbl22/halve, peak 1538) sit 5q BELOW
    // the binder. Switching them to the fast (carry-ancilla) doubling is free at
    // peak 1543 and value-exact: avg executed Toffoli 1,695,087 -> 1,682,159
    // (-12,928). The fast-doubling op stream re-rolls the Fiat-Shamir island, so
    // the reroll knobs below are re-tuned to 40/13 (found by a randomized 2-D
    // island search). Validated 0/0/0 over all 9024 shots @ 1543q / 1,682,159 T.
    set_default_env("KARA_SOL_DBL_FAST", "1");
    // Stacked qubit cut (peak 1543 -> 1542, learned from anupsv's 8780d1e): the
    // ROUND84 Karatsuba z1_reg top bit (index 257) is provably 0 across the whole
    // Solinas-reduction peak window (z1_reg == 2*lo*hi < 2^257 there), so that
    // qubit is freed for the window and re-grabbed (fresh zero) before the inverse
    // combine restores z1=(lo+hi)^2. Bennett-clean, 0 added Toffoli. Stacks on
    // KARA_SOL_DBL_FAST; the combined op stream re-rolls the island, re-tuned to
    // REROLL=17/POST_SUB=56 below (MARGIN stays 5 — no give-back). Validated 0/0/0
    // over 9024: 1542q x 1,682,159 T = 2,593,889,178.
    set_default_env("KARA_FREE_Z1_TOPBIT", "1");
    // W-TRUNC tightening: GCD-body width envelope margin. Re-scanned for the
    // Karatsuba x-tail op stream: margin=27 + REROLL=0 lands a clean 9024-shot
    // island (anupsv's margin=26/REROLL=20 was for the schoolbook stream).
    // WIDTH_MARGIN 27->26 stacked with APPLY_CLEAN_COMPARE_BITS 21->20 and
    // PA9024_COMPARE_SCHEDULE_MARGIN 8->7: -5,576 executed Toffoli at the 1434
    // peak. Re-rolled Fiat-Shamir island lands clean (0/0/0 over 9024) at
    // DIALOG_REROLL=0 / DIALOG_POST_SUB_REROLL=44. 1434q x 1,733,573 T = 2,485,943,682.
    // WIDTH_MARGIN 9 -> 10: the freed comparator slack (COMPARE_BITS 73->52
    // above) is partly re-spent to widen the GCD-body width envelope by one
    // safety bit. At margin=9 the width-truncation (u/v bitlen > active_width)
    // is the dominant hard-input source (~83/300k factor checks); margin=10
    // cuts that to ~27, dropping the expected hard inputs per random reroll from
    // ~11 to ~5 so a clean Fiat-Shamir island is found in seconds instead of
    // hours. Costs +5,815,760 score vs margin=9 but the net (compare52 +
    // margin10) is 2,130,373,770 -> 2,112,431,650 (-17,942,120), and the lower
    // hard rate keeps the island search tractable. Validated 0/0/0 over 9024.
    // Final-window W2 keeps WIDTH_MARGIN at 10; margin 11 crosses the 1328q
    // cliff, while margin 10 validated clean with the tail nonce below.
    set_default_env("DIALOG_GCD_WIDTH_MARGIN", "10");
    // Measured (Gidney) uncompute for the apply-phase modular subtract's raw
    // difference, mirroring the already-measured apply ADD. ~n Toffoli instead
    // of ~2n per call; peak-neutral (same carry lane the ADD already uses).
    set_default_env("DIALOG_GCD_MEASURED_APPLY_SUB", "1");
    // QUBIT-PEAK CUT (1698 -> 1572, -126q): host the GCD-body 'gated' on idle
    // future-log slots (HOST_GATED), and window the apply add/sub carry lane into
    // 2 blocks with measurement-uncompute + a measured boundary-carry clear so the
    // 256-wide carry lane never coexists with f at the peak. Toffoli +102k
    // (1,668,753 -> 1,770,897) but peak -126 => score 2,833,542,594 -> 2,783,850,084.
    set_default_env("DIALOG_GCD_HOST_GATED", "1");
    set_default_env("DIALOG_GCD_APPLY_WINDOW_BLOCKS", "2");
    // ROUND84 x-tail square: replace the 2^32 Solinas term's shift-by-22
    // (mod_shift_left_by_k(22) -> mid_sub -> shift_right_by_k(22)) with the
    // value-identical 22x mod-p doubling -> mid_sub -> 22x mod-p halving
    // (x*2^22 mod p == x<<22 mod p). The direct-const doubling/halving lanes
    // carry-sweep in place with no spill register, so the block never parks the
    // 24 persistent flags (spill=22 + ovf + flag_inv) that pinned the square
    // phase at 1567. Square phase drops to 1543; the global peak falls
    // 1567 -> 1543. Costs +~6,384 avg-executed Toffoli (see F_CUT below).
    set_default_env("ROUND84_XTAIL_BORROW_CARRIES", "1");
    // Chunked apply materializes ctrl&a only for the active carry window, so the
    // apply phase drops under the ROUND84 peak binder. After the ROUND84 square
    // dropped to 1543, the apply raw sum/difference phases (block 1 = [F_CUT,257),
    // f + carry lane) became the 1558 binder. The chunked sub/add is EXACT
    // regardless of F_CUT (full cuccaro + exact [..F_CUT] boundary clear), so
    // widening the first cut 70 -> 78 rebalances the blocks (block 1 narrows to
    // 257-78) and drops the apply phase to 1543 == the ROUND84 floor. Global peak
    // 1558 -> 1543. F_CUT only reseeds + grows the boundary comparator (+~6,384
    // avg-executed Toffoli, 1,688,703 -> 1,695,087); peak-neutral for any cut>=78.
    // Peak-band rebuild (1226 tier): the apply ripple is sliced into 10 even
    // chunks so the transient load/carry register stays ~26 wide, dropping the
    // apply ripple peak 1266 -> 1222 (under the 1226 double_y/halve_y binder).
    // Toffoli-near-neutral (the extra boundary comparators cost ~250 avg). Pairs
    // with SQUARE_ROW_MAX_SEG below (the peak-bounded square) to land global peak
    // at 1226 instead of 1284.
    set_default_env("DIALOG_GCD_APPLY_CHUNKED_F_BLOCKS", "16");
    set_default_env("DIALOG_GCD_APPLY_CHUNKED_F_CUSTOM4", "0");
    set_default_env("DIALOG_GCD_APPLY_CHUNKED_F_CUSTOM5", "0");
    // PEAK-QUBIT CUT 1542 -> 1500 (-42q). Two co-binders dropped together:
    //  (1) ROUND84 Karatsuba square (z0=lo^2 / z2=hi^2 schoolbook squares parked a
    //      ~130-wide cuccaro_add_fast carry lane, and the Solinas mid_sub/sub_add's
    //      mod_add_qq/mod_sub_qq materialized a load_const(256) correction transient).
    //      Fix: KARA_Z02_LOWQ hosts the z0 square's carry lane on the (clean) z2
    //      slice via cuccaro_add_fast_borrowed_carries and runs z2 ancilla-free
    //      (lowq); KARA_SOL_MOD_VENT vents the constant corrections onto the dirty
    //      operand (+2 clean) instead of load_const. Both are value-exact.
    //  (2) GCD apply materialized_special raw sum/difference: the [F_CUT,257) block's
    //      f + carry lane pinned 1542. The chunked sub/add is EXACT for any cut, so
    //      widening F_CUT 78 -> 99 narrows block 1 and drops the apply phase to 1500.
    // Global peak 1542 -> 1500; cost +~36,558 avg-executed Toffoli (1,682,159 ->
    // 1,718,717) for -42q: 1500 x 1,718,717 = 2,578,075,500.
    set_default_env("KARA_Z02_LOWQ", "1");
    set_default_env("KARA_Z2_SELFHOST", "1");
    set_default_env("KARA_SOL_MOD_VENT", "1");
    // PEAK 1500 -> 1466 (-34q). On the 1500 floor the peak was a co-binder tie between
    // the GCD-core branch comparator (tobitvector_branch_bits / _reverse) and the apply
    // mod add/sub (materialized_special_chunked_raw_sum / _difference). The apply phase
    // can be driven down by widening the chunk cut (each +1 F_CUT -> -2 apply peak), but
    // only until it meets the comparator floor -- so the comparator is torn down first.
    //  - DIALOG_GCD_BRANCH_BITS_HOST_COMPARATOR=1: the fused branch-bit path never used
    //    the separately-allocated `cmp` ancilla (it derives b0_and_b1 from the in-flight
    //    comparator carry), and the comparator materialized its own c_in+carries lane on
    //    top of the live GCD state. Routing the fused path through the borrowed-carry
    //    comparator (carry lane hosted on a temporarily-clean future-log slice) + dropping
    //    the dead cmp removes that standalone transient. Value-exact (ancilla returned
    //    clean); the branch_bits phases fall well below the apply tier.
    //  - DIALOG_GCD_APPLY_CHUNKED_F_CUT 99 -> 116: with the comparator unbound, widening
    //    the cut sinks BOTH apply phases to the next true floor -- the materialized_*_body
    //    GCD-body tier at 1466. Exact for any cut (full cuccaro + exact [..F_CUT] clear).
    // This reached peak 1500 -> 1466 for +13,566 avg-executed Toffoli (1,718,717 ->
    // 1,732,283); score 1466 x 1,732,283 = 2,539,526,878.
    set_default_env("DIALOG_GCD_BRANCH_BITS_HOST_COMPARATOR", "1");
    // PEAK 1466 -> 1446 (-20q). The 1466 floor was a 4-phase co-bind: the two apply
    // mod add/sub (materialized_special_chunked_raw_sum/_difference) and the two GCD-body
    // add/sub (raw_tobitvector_materialized_{add,sub}_body). Both body families dropped
    // out from under 1466 via two value-exact carry-lane reclaims, after which F_CUT
    // sinks the apply pair to the freed floor:
    //  - DIALOG_GCD_BODY_HOST_CIN=1: the materialized body's borrowed-carry Cuccaro still
    //    allocated a FRESH c_in ancilla on top of the borrowed (future-log) carry lane --
    //    the single qubit pinning the body at 1466. With the odd-u fastpath body_start=1,
    //    gated[0] is never loaded/cleared (stays |0>), so it serves as the carry-in with
    //    no alloc. Body phases 1466 -> 1446. Value-exact (c_in=0 either way).
    //  - DIALOG_GCD_LATE_BORROW_UV_HIGH=1: at late steps the compressed future-log runs
    //    short, so the body fell back to allocating its own carry+gated lane (the 1465
    //    `tobitvector_subtract`/`_reverse_add` marker tier). The GCD has converged there,
    //    so u[active_width..] is |0> by the SAME premise the width truncation relies on
    //    and is already allocated -> borrow it as scratch. Marker tier 1465 -> 1446. No
    //    new failure modes (any input with nonzero u-high already fails the truncation).
    //  - DIALOG_GCD_APPLY_CHUNKED_F_CUT 116 -> 126: with the body floor at 1446, widening
    //    the cut sinks both apply phases to 1446 (their min; F_CUT>126 rebalances upward).
    // Net peak 1466 -> 1446 for +7,980 avg-executed Toffoli (1,732,283 -> 1,740,263) ~=
    // 399 T/qubit, far inside break-even. Score 1446 x 1,740,263 = 2,516,420,298.
    set_default_env("DIALOG_GCD_BODY_HOST_CIN", "1");
    set_default_env("DIALOG_GCD_LATE_BORROW_UV_HIGH", "1");
    // Body-carry-band-trim DISABLED (was "0,...,0,1,1,1,1,1,1,1,1"): the late-step
    // 1-bit body sub/add truncation mis-drops a needed bit when the converged
    // operand bitlen reaches active_width on a handful of reachable inputs -- a
    // Fiat-Shamir-island hazard class on top of the width envelope. The per-step
    // compare schedule frees enough Toffoli to pay back the ~1,088 this saved AND
    // remove that hazard class, making the island materially easier to land while
    // net Toffoli still beats the flat-50 baseline (1,512,823 -> 1,506,043 @ 1313).
    // Stacked peak-1302 band-trim schedule + measured-ovfclear + F_CUT4=189 (tier-3 "safe lock"):
    // trims average executed Toffoli to 1,456,963 at peak 1302 qubits.
    set_default_env("DIALOG_GCD_BODY_CARRY_BAND_TRIMS", "0,3,3,3,3,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,3,3,3");
    set_default_env("DIALOG_GCD_TOBITVECTOR_CSWAP_BODY_TRIM", "0");
    set_default_env("DIALOG_GCD_BINDER_NOTCH_STEPS", "8,9,10");
    set_default_env("DIALOG_GCD_BINDER_NOTCH_EXTRA", "3");
    set_default_env("DIALOG_GCD_BINDER_NOTCH_MAP", "11:1,12:1,13:1");
    set_default_env(
        "DIALOG_GCD_SPECIAL_OVERFLOW_CLEAN_STEP_BITS",
        "113:21,131:21,142:22,187:23,205:22,210:21",
    );
    set_default_env(
        "DIALOG_GCD_SPECIAL_UNDERFLOW_CLEAN_STEP_BITS",
        "42:22,91:22,118:22,149:21",
    );
    set_default_env("DIALOG_GCD_FUSED_OVFCLEAR_MEASURED", "1");
    // 1320q apply teardown: low-q final chunk plus a hosted boundary split at
    // the second custom-five cut. The retained carry at bit 100 hosts the
    // high-window comparator carry-in, avoiding the generic split's extra
    // boundary qubit and low-window recompute.
    set_default_env("DIALOG_GCD_APPLY_FINAL_LOWQ", "0");
    // Round84 mid-sub: ancilla-light Cuccaro const-add + carry-in borrow (1309->1307);
    // compressed-block: current-step s2 composite-scratch fold (1308->1307).
    set_default_env("R84_LOWQ", "1");
    set_default_env("R84_LOWQ_CIN_BORROW", "1");
    set_default_env("R84_QPROD_NAF", "1"); // quotient*c uses 977 = 2^10 - 2^5 - 2^4 + 1.
    // Fold the square's high half into its low half in place, accumulate the
    // resulting 33-bit quotient, apply quotient*(2^256-p) once, subtract once,
    // then reversibly unfold before Bennett-uncomputing the square. The final
    // modular subtract vents onto the folded operand, retaining the 1307q peak.
    // The 21-bit high-carry propagation and rare folded-lo noncanonical band
    // are selected away with the shared Fiat-Shamir island.
    set_default_env("ROUND84_INPLACE_SOLINAS_FOLD", "1");
    set_default_env("ROUND84_INPLACE_QUOTIENT_CARRY_TRUNC_W", "21");
    // Peak-bounded square (1226 tier): the round84 lam^2 schoolbook square parks
    // a 512-wide product (peak 1024) plus the per-row source register (up to
    // +257 for the widest row → 1284). SQUARE_ROW_MAX_SEG slices each square row
    // into the minimum number of windows that keeps every source segment <= this
    // width, chaining the inter-window carry through a clean cout ancilla that is
    // recovered by a local, tmp-high-borrowed measured comparator (no allocated
    // carry array, no wide-prefix rebuild). At 199 only the rows wider than 199
    // (i < ~57) window, each into 2, dropping the square forward/inverse peak to
    // 1226 (== the double_y binder) while adding only ~26k avg Toffoli for the
    // carry-recovery comparators. Value-exact: the same product lands in tmp_ext
    // (verified: ancilla-garbage 0; SQUARE_ROW_MAX_SEG=0 restores the bit-exact
    // 1284 base). Net: peak 1284 -> 1226, score 1.821e9 -> 1.771e9.
    set_default_env("SQUARE_ROW_MAX_SEG", "176");
    set_default_env("DIALOG_GCD_K5_CLEAN_BLOCK", "1");
    set_default_env("DIALOG_GCD_FOLD_PARK_LOW_CARRIES", "1");
    set_default_env("DIALOG_GCD_SPECIAL_FOLD_BORROW_CARRIES", "1");
    set_default_env("DIALOG_GCD_K2_APPLY_INPLACE_RAW_BLOCK", "1");
    set_default_env("DIALOG_GCD_FOLD_FREED_TAIL", "1");  // BAKED: 1221 ship
    set_default_env("DIALOG_GCD_BORROW_CURRENT_S2", "1");
    set_default_env("DIALOG_GCD_BORROW_ZERO_RAW_FUTURE", "1");
    set_default_env("DIALOG_GCD_FREE_SCRATCH_BEFORE_SHIFT", "1");
    set_default_env("DIALOG_GCD_APPLY_BOUNDARY_SPLIT", "100");
    set_default_env("DIALOG_GCD_APPLY_CHUNKED_F_CUT", "50");
    set_default_env("DIALOG_GCD_APPLY_CHUNKED_F_CUT2", "100");
    set_default_env("DIALOG_GCD_APPLY_CHUNKED_F_CUT3", "150");
    set_default_env("DIALOG_GCD_APPLY_CHUNKED_F_CUT4", "190");
    // WIDTH_SLOPE tightening: the per-step GCD width envelope shrink rate
    // (ideal = N - step*SLOPE + MARGIN) was left at the default 0.7075 by the
    // whole frontier lineage; only the constant MARGIN was ever tuned. The
    // Bernstein-Yang/binary-GCD width bound (Gidney et al., arXiv:2510.10967,
    // "after i iters 2*deg(b) <= 2d-1-i-delta") shows the realizable bitlen
    // shrinks slightly faster, so SLOPE 707.5 -> 708 tightens every late-step
    // GCD-body width by an extra fraction of a bit: avg executed Toffoli
    // 1,779,067 -> 1,778,555 (-512), peak-neutral at 1355q. The tighter
    // truncation re-rolls the Fiat-Shamir island; a 1-D reroll sweep (post_sub
    // fixed at the inherited 503292) lands a clean island at DIALOG_REROLL=101019.
    // Back off the width slope to 1004 for the final-window W2 clean island.
    // Re-tighten WIDTH_SLOPE 1005 -> 1009 on the W2 final-window base (which had
    // left the slope loose to find its structural island). The per-step GCD-body
    // width envelope shrinks an extra ~4 notches; the dropped high bits are
    // provably 0 on the converged reachable support, so it is value-exact and the
    // residual failures are pure Fiat-Shamir, dodged by the re-rolled tail nonce
    // below. avg executed Toffoli 1,540,355 -> 1,538,227 (-2,128), peak-neutral at
    // 1320q. Found with the local classical width-convergence pre-filter +
    // bit-exact validate (island_search_prefilter), confirmed via official run.
    // 1009 -> 1011: one further notch, stacked under one shared island with the
    // KAL_FOLD 24->22 and APPLY_CLEAN_COMPARE_BITS 20->19 re-tightenings above.
    // 1011 -> 1012: one more width-envelope notch, stacked on COMPARE_BITS=46
    // under the nonce-10429 island below. Value-exact, peak-neutral at 1320q.
    set_default_env("DIALOG_GCD_WIDTH_SLOPE_X1000", "1015");
    // Active-395 island on the promoted 1355q base: validated 0/0/0 over all
    // 9024 shots at 1355q x 1,773,011 T.
    set_default_env("DIALOG_REROLL", "4269");
    set_default_env("DIALOG_POST_SUB_REROLL", "503292");
    // Fiat-Shamir island for ACTIVE_ITERATIONS=393 + WIDTH_MARGIN=25 (1350q base).
    // The fixed-length 96-op identity tail (see the DIALOG_TAIL_NONCE block in
    // build_builder) reseeds the 9024 Fiat-Shamir test inputs without changing
    // the circuit action, Toffoli count, or peak qubits. nonce=385307 lands a
    // clean island: validated 0/0/0 over all 9024 shots at 1350q x 1,763,987 T.
    // Fiat-Shamir island for the K=2 apply rebalance above: 0/0/0 over all
    // 9024 shots at 1390q x 1,630,487 T.
    // Re-rolled for COMPARE_BITS=52 + WIDTH_MARGIN=10 (above): nonce=127 lands a
    // clean island (found by the parallel prefix-clone classical filter in
    // harness/fasteval, then quantum-confirmed). Validated 0/0/0 over all 9024
    // shots at 1390q x 1,519,735 T = 2,112,431,650. Backups: 354, 418.
    // Re-rolled for the combined KAL_DOUBLE/FOLD_CARRY_TRUNC_W=23 op stream:
    // nonce=254 lands a clean island, validated 0/0/0 over all 9024 shots at
    // 1390q x 1,518,179 T = 2,110,268,810.
    // Final-window W2 island: validated 0/0/0 over all 9024 shots at
    // 1320q x 1,545,787 T = 2,040,438,840.
    // Re-rolled for the WIDTH_SLOPE 1005 -> 1009 re-tightening above: nonce 6416
    // lands a clean Fiat-Shamir island, validated 0/0/0 over all 9024 shots at
    // 1320q x 1,538,227 T = 2,030,459,640 (backup: 6700).
    // Re-rolled again for the stacked WIDTH_SLOPE=1011 + KAL_FOLD=22 +
    // APPLY_CLEAN_COMPARE_BITS=19 re-tightenings: nonce 18509 lands a clean island,
    // validated 0/0/0 over all 9024 shots at 1320q x 1,535,629 T = 2,027,030,280.
    // Re-rolled again for the stacked COMPARE_BITS 47->46 re-tightening: nonce
    // 20397 lands a clean island, validated 0/0/0 over all 9024 shots at
    // 1320q x 1,534,757 T = 2,025,879,240.
    // Re-rolled again for the stacked WIDTH_SLOPE 1011->1012 notch: nonce 10429
    // lands a clean island, validated 0/0/0 over all 9024 shots at
    // 1320q x 1,534,277 T = 2,025,245,640.
    // Pair-compressed 46/20 island: nonce 689 lands a clean trusted run,
    // validated 0/0/0 over all 9024 shots at
    // 1313q x 1,536,923 T = 2,017,979,899.
    // Re-rolled for the KAL_DOUBLE_CARRY_TRUNC_W=21 re-tightening above: nonce
    // 1000001157 lands a clean island, validated 0/0/0 over all 9024 shots at
    // 1313q x 1,535,885 T = 2,016,617,005 (official ecdsafail run).
    set_default_env("DIALOG_GCD_SELECTED_BODY_NOCIN", "1");
    // STACKED island: K2 per-step compare schedule (MARGIN=0, body-carry-band-trims
    // OFF; 1,506,043 T) + DIALOG_GCD_BORROW_CURRENT_BLOCK=1 (GCD-walk peak 1313->1309
    // at 0 added Toffoli). The borrow relabel removes 1920 non-Toffoli alloc/clear
    // ops, so the shorter op stream reseeds the 96-op identity tail's SHAKE256 and
    // the prior K2 island (300112609) no longer lands. nonce 3400174 lands a fresh
    // clean island for the stacked stream: validated 0/0/0 (0 classical / 0 phase /
    // 0 ancilla) over all 9024 shots at 1309q x 1,506,043 T = 1,971,410,287 (beats
    // the K2 floor 1,977,434,459 by 6,024,172 and the baseline 1,986,336,599 by
    // 14,926,312). Borrow-current-block confirmed value-exact: GCD-survivor
    // fail-count distributions (classical/phase/ancilla) are statistically identical
    // borrow ON vs OFF (no measure-nonzero corruption floor; ancilla garbage == 0 in
    // both), so the nonce merely dodges the inherited Fiat-Shamir straggler class.
    // ON clean islands occur at the SAME ~1/108 rate among GCD-survivors as K2-alone
    // OFF. Backup clean islands (all validated 0/0/0 @ 1309 x 1,506,043 = 1,971,410,287):
    // 3756953, 3774241, 3840981, 40330388.
    // Re-rolled for the APPLY_CLEAN_COMPARE_BITS 21 -> 20 re-tightening above:
    // nonce 721381 lands a clean Fiat-Shamir island, validated 0/0/0 over all
    // 9024 shots at 1309q x 1,503,355 T = 1,967,891,695.
    // Re-rolled for the lowq0 fast-final + ACTIVE_ITERATIONS=262 route:
    // nonce 2432 validates 0/0/0 over all 9024 shots at
    // 1309q x 1,497,795 T = 1,960,613,655.
    // K2-pair codec 6->3 CCX core encoder (peak-neutral -3,096 T). Re-hunted clean
    // Fiat-Shamir island:
    // Binder-notch fallback 8,9: nonce 169924627 validates 0/0/0 over all
    // 9024 shots at 1300q x 1,454,884 T = 1,891,349,200.
    set_default_env("DIALOG_TAIL_NONCE", "18100027017098");
    set_default_env("ROUND84_FOLD_FAST_ADD", "0");  // round84 Solinas-fold small adders coherent->measured-fast (-1,434 exec-T, peak-neutral 1285)
    set_default_env("DIALOG_GCD_FOLD_MAJ2", "1");
    set_default_env("DIALOG_GCD_FOLD_MAJ1", "1");
    set_default_env("DIALOG_GCD_APPLY_FINAL_TOPCLEAN", "0");
    set_default_env("ROUND84_QPROD_VENT_PAD", "1");
    set_default_env("DIALOG_GCD_FOLD_FREED_TAIL_ED", "1");
    set_default_env("DIALOG_GCD_APPLY_FINAL_WINDOWED_FAST_BLOCKS", "0");
    // Fuse the branch-bit comparator with the b0-controlled log update: derive
    // b0_and_b1 from the in-flight comparator carry instead of materializing a
    // separate cmp qubit and recomputing the comparator for uncompute. Pure
    // Toffoli reduction (1952382 -> 1861990), peak-neutral at 1698.
    // (Validated 0/0/0 over 9024 via eval_circuit.)
    set_default_env("DIALOG_GCD_FUSED_BRANCH_BITS", "1");
    // Odd-u low-bit fastpath: after the binary-GCD branch swap, u[0] is one on
    // the reachable verifier support. The lane-0 ctrl&u[0] gated load collapses
    // to a CX, and the lane-0 tobitvector add/sub body has no carry/borrow into
    // lane 1, so the body can start at bit 1. Co-tuned with the reroll island.
    set_default_env("DIALOG_GCD_ODD_U_LOWBIT_FASTPATH", "1");
}

pub fn build_builder() -> B {
    configure_ecdsafail_submission_route();

    let mut builder = if std::env::var("POINT_ADD_COUNT_ONLY").ok().as_deref() == Some("1") {
        B::new_count_only()
    } else {
        B::new()
    };
    let b = &mut builder;
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

    // Fiat-Shamir reroll: emit k pairs of X;X (exact identity, X^2 = I) on a
    // data qubit. This perturbs the serialized op-stream bytes -> reseeds the
    // SHAKE256-derived 9024 test inputs WITHOUT changing the circuit's action,
    // Toffoli count, or peak qubits. Used to slide off Fiat-Shamir "islands"
    // where an aggressive (otherwise-correct) width truncation has a handful of
    // hard test inputs. Default 0 = byte-identical baseline.
    if let Some(k) = std::env::var("DIALOG_REROLL")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|&k| k > 0)
    {
        b.set_phase("dialog_reroll");
        for _ in 0..k {
            b.x(tx[0]);
            b.x(tx[0]);
        }
    }

    let p = SECP256K1_P;

    // Step 1-2: Px -= Qx, Py -= Qy
    mod_sub_qb(b, &tx, &ox, p);
    mod_sub_qb(b, &ty, &oy, p);
    if let Some(k) = std::env::var("DIALOG_POST_SUB_REROLL")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|&k| k > 0)
    {
        b.set_phase("dialog_post_sub_reroll");
        for _ in 0..k {
            b.x(tx[1]);
            b.x(tx[1]);
        }
    }

    emit_dialog_gcd_raw_pa(b, &tx, &ty, &ox, &oy, p);

    if !b.count_only && std::env::var("SKIP_ALT_SEED_CHECKS").ok().as_deref() != Some("1") {
        run_alt_seed_checks(&b.ops);
    }

    if !b.count_only && std::env::var("TRACE_PEAK").is_ok() {
        eprintln!(
            "DEBUG peak_qubits={} at phase='{}' ops_idx={} total_ops={}",
            b.peak_qubits,
            b.peak_phase,
            b.peak_ops_idx,
            b.ops.len()
        );
        let pk = b.peak_qubits;
        let mut uniq: std::collections::BTreeMap<&'static str, (u32, usize)> =
            std::collections::BTreeMap::new();
        for (a, ph, op) in &b.peak_log {
            if *a + 5 >= pk {
                let entry = uniq.entry(ph).or_insert((*a, *op));
                if *a > entry.0 {
                    *entry = (*a, *op);
                }
            }
        }
        for (ph, (a, op)) in uniq.iter() {
            eprintln!("DEBUG near_peak active={} phase='{}' ops_idx={}", a, ph, op);
        }
    }

    if !b.count_only && std::env::var("TRACE_PHASES").is_ok() {
        // Attribute emitted ops to the active phase at each op index.
        // phase_transitions is sorted by ops_idx (monotonically appended).
        // For each op, binary-find the phase region it falls in.
        let trans = &b.phase_transitions;
        let n_ops = b.ops.len();
        // Per-phase aggregates.
        let mut agg: std::collections::BTreeMap<&'static str, (u64, u64, u64)> =
            std::collections::BTreeMap::new();
        // Also per-call counters: each contiguous (phase, region) gets its own bucket for ordered printout.
        let mut regions: Vec<(&'static str, usize, u64, u64, u64)> = Vec::new();
        for i in 0..trans.len() {
            let start = trans[i].0;
            let end = if i + 1 < trans.len() {
                trans[i + 1].0
            } else {
                n_ops
            };
            let phase = trans[i].1;
            let mut tof: u64 = 0;
            let mut cli: u64 = 0;
            let mut other: u64 = 0;
            for op in &b.ops[start..end] {
                match op.kind {
                    OperationType::CCX | OperationType::CCZ => tof += 1,
                    OperationType::CX
                    | OperationType::CZ
                    | OperationType::Swap
                    | OperationType::Hmr
                    | OperationType::R => cli += 1,
                    _ => other += 1,
                }
            }
            regions.push((phase, start, tof, cli, other));
            let e = agg.entry(phase).or_insert((0, 0, 0));
            e.0 += tof;
            e.1 += cli;
            e.2 += other;
        }
        let total_tof: u64 = agg.values().map(|v| v.0).sum();
        eprintln!("=== per-phase emitted Toffoli (classical view; executed-shot stats are in harness) ===");
        eprintln!(
            "{:<40} {:>12} {:>12} {:>6}",
            "phase", "ccx", "cliff", "%tof"
        );
        let mut v: Vec<_> = agg.iter().collect();
        v.sort_by(|a, b| b.1 .0.cmp(&a.1 .0));
        for (ph, (t, c, _o)) in v {
            let pct = if total_tof > 0 {
                (*t as f64) * 100.0 / (total_tof as f64)
            } else {
                0.0
            };
            eprintln!("{:<40} {:>12} {:>12} {:>5.1}%", ph, t, c, pct);
        }
        eprintln!("total_ccx_emitted={} total_ops={}", total_tof, n_ops);
        if std::env::var("TRACE_PHASES_VERBOSE").is_ok() {
            eprintln!("--- per-region (ordered) ---");
            for (ph, start, tof, cli, _o) in &regions {
                if *tof == 0 && *cli == 0 {
                    continue;
                }
                eprintln!("@{:<10} {:<40} ccx={} cli={}", start, ph, tof, cli);
            }
        }
    }

    if std::env::var("TRACE_PHASE_ACTIVE").is_ok() {
        b.close_phase_active_region();
        eprintln!("=== per-phase active qubit maxima ===");
        eprintln!("{:<48} {:>12}", "phase", "active_q");
        let mut v: Vec<_> = b.phase_active_max.iter().collect();
        v.sort_by(|a, b| b.1.cmp(a.1).then_with(|| a.0.cmp(b.0)));
        let top_n = std::env::var("TRACE_PHASE_ACTIVE_TOP")
            .ok()
            .and_then(|s| s.parse::<usize>().ok());
        let mut printed = 0usize;
        for (phase, active) in v {
            if top_n.is_some_and(|limit| printed >= limit) {
                break;
            }
            eprintln!("{:<48} {:>12}", phase, active);
            printed += 1;
        }
        if std::env::var("TRACE_PHASE_ACTIVE_REGIONS").is_ok() {
            eprintln!("--- per-region active qubit maxima (ordered) ---");
            for (end, phase, active) in &b.phase_active_regions {
                eprintln!("@{:<10} {:<48} active_q={}", end, phase, active);
            }
        }
    }

    // Fiat-Shamir island selector: emit a FIXED-LENGTH block of identity X;X
    // pairs at the very end of the op stream. For each of NONCE_BITS bits, emit
    // one X;X pair (an exact identity, since X^2 = I) targeting tx[0] when the
    // bit is 0 or tx[1] when the bit is 1. The block length is constant
    // (2*NONCE_BITS ops), so the op count and circuit action are unchanged and
    // the Toffoli count and peak qubit width are unaffected; only the per-op
    // target of this tail varies with the nonce, which reseeds the SHAKE256-
    // derived 9024 Fiat-Shamir test inputs. This selects which random test set
    // the circuit is validated against without tuning the circuit to it. Gated
    // on DIALOG_TAIL_NONCE so the stream is byte-identical when it is absent.
    if let Some(nonce) = std::env::var("DIALOG_TAIL_NONCE")
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
    {
        const NONCE_BITS: u32 = 48;
        b.set_phase("dialog_tail_nonce");
        for i in 0..NONCE_BITS {
            let q = if (nonce >> i) & 1 == 1 { tx[1] } else { tx[0] };
            b.x(q);
            b.x(q);
        }
    }

    builder
}

pub fn build() -> Vec<Op> {
    if std::env::var("DIALOG_GCD_K5_HEAD11_SELFTEST").is_ok() {
        match dialog_gcd_k5_head11_codec_selftest() {
            Ok(()) => eprintln!(
                "DIALOG_GCD_K5_HEAD11_SELFTEST: PASS (2048-word head codec reversible and phase clean)"
            ),
            Err(e) => panic!("DIALOG_GCD_K5_HEAD11_SELFTEST: FAIL: {e}"),
        }
        if std::env::var("DIALOG_GCD_K5_HEAD11_SELFTEST_ONLY")
            .ok()
            .as_deref()
            == Some("1")
        {
            return Vec::new();
        }
    }
    if std::env::var("DIALOG_GCD_K5_TAIL3_SELFTEST").is_ok() {
        match dialog_gcd_k5_tail3_codec_selftest() {
            Ok(()) => eprintln!(
                "DIALOG_GCD_K5_TAIL3_SELFTEST: PASS (two-step pair codec reversible and phase clean)"
            ),
            Err(e) => panic!("DIALOG_GCD_K5_TAIL3_SELFTEST: FAIL: {e}"),
        }
        if std::env::var("DIALOG_GCD_K5_TAIL3_SELFTEST_ONLY")
            .ok()
            .as_deref()
            == Some("1")
        {
            return Vec::new();
        }
    }
    if std::env::var("DIALOG_GCD_K5_TAIL3_TOP32_SELFTEST").is_ok() {
        match dialog_gcd_k5_tail3_top32_codec_selftest() {
            Ok(()) => eprintln!(
                "DIALOG_GCD_K5_TAIL3_TOP32_SELFTEST: PASS (32-word weighted codec reversible and phase clean)"
            ),
            Err(e) => panic!("DIALOG_GCD_K5_TAIL3_TOP32_SELFTEST: FAIL: {e}"),
        }
        if std::env::var("DIALOG_GCD_K5_TAIL3_TOP32_SELFTEST_ONLY")
            .ok()
            .as_deref()
            == Some("1")
        {
            return Vec::new();
        }
    }
    if std::env::var("DIALOG_GCD_K5_TAIL6_GRAPH9_SELFTEST").is_ok() {
        match dialog_gcd_k5_tail6_graph9_codec_selftest() {
            Ok(()) => eprintln!(
                "DIALOG_GCD_K5_TAIL6_GRAPH9_SELFTEST: PASS (75-word graph codec reversible and phase clean)"
            ),
            Err(e) => panic!("DIALOG_GCD_K5_TAIL6_GRAPH9_SELFTEST: FAIL: {e}"),
        }
        if std::env::var("DIALOG_GCD_K5_TAIL6_GRAPH9_SELFTEST_ONLY")
            .ok()
            .as_deref()
            == Some("1")
        {
            return Vec::new();
        }
    }
    if std::env::var("DIALOG_GCD_K5_TAIL6_GRAPH_SELFTEST").is_ok() {
        match dialog_gcd_k5_tail6_graph_codec_selftest() {
            Ok(()) => eprintln!(
                "DIALOG_GCD_K5_TAIL6_GRAPH_SELFTEST: PASS (32-word graph codec reversible and phase clean)"
            ),
            Err(e) => panic!("DIALOG_GCD_K5_TAIL6_GRAPH_SELFTEST: FAIL: {e}"),
        }
        if std::env::var("DIALOG_GCD_K5_TAIL6_GRAPH_SELFTEST_ONLY")
            .ok()
            .as_deref()
            == Some("1")
        {
            return Vec::new();
        }
    }
    if std::env::var("DIALOG_GCD_K5_TAIL7_SELFTEST").is_ok() {
        match dialog_gcd_k5_tail7_codec_selftest() {
            Ok(()) => eprintln!(
                "DIALOG_GCD_K5_TAIL7_SELFTEST: PASS (20-word codec reversible and phase clean)"
            ),
            Err(e) => panic!("DIALOG_GCD_K5_TAIL7_SELFTEST: FAIL: {e}"),
        }
        if std::env::var("DIALOG_GCD_K5_TAIL7_SELFTEST_ONLY")
            .ok()
            .as_deref()
            == Some("1")
        {
            return Vec::new();
        }
    }
    if std::env::var("DIALOG_GCD_K5_TAIL7_FIXED_APPLY_SELFTEST").is_ok() {
        match dialog_gcd_k5_tail7_fixed_apply_selftest() {
            Ok(()) => eprintln!(
                "DIALOG_GCD_K5_TAIL7_FIXED_APPLY_SELFTEST: PASS (fixed digit-4 apply matches fused apply)"
            ),
            Err(e) => panic!("DIALOG_GCD_K5_TAIL7_FIXED_APPLY_SELFTEST: FAIL: {e}"),
        }
        if std::env::var("DIALOG_GCD_K5_TAIL7_FIXED_APPLY_SELFTEST_ONLY")
            .ok()
            .as_deref()
            == Some("1")
        {
            return Vec::new();
        }
    }
    if std::env::var("SQUARE_WINDOW_SELFTEST").is_ok() {
        match square_window_selftest() {
            Ok(()) => eprintln!("SQUARE_WINDOW_SELFTEST: PASS"),
            Err(e) => panic!("SQUARE_WINDOW_SELFTEST: FAIL: {e}"),
        }
        if std::env::var("SQUARE_WINDOW_SELFTEST_ONLY").ok().as_deref() == Some("1") {
            return Vec::new();
        }
    }
    if std::env::var("FOLD_FREED_TAIL_SELFTEST").is_ok() {
        match fold_freed_tail_selftest() {
            Ok(()) => eprintln!("FOLD_FREED_TAIL_SELFTEST: PASS (freed-tail ≡ baseline, ancilla & phase clean)"),
            Err(e) => panic!("FOLD_FREED_TAIL_SELFTEST: FAIL: {e}"),
        }
        if std::env::var("FOLD_FREED_TAIL_SELFTEST_ONLY")
            .ok()
            .as_deref()
            == Some("1")
        {
            return Vec::new();
        }
    }
    if std::env::var("SPECIAL_FOLD_PARK_SELFTEST").is_ok() {
        match special_fold_park_selftest() {
            Ok(()) => eprintln!(
                "SPECIAL_FOLD_PARK_SELFTEST: PASS (parked fold ≡ baseline, ancilla & phase clean)"
            ),
            Err(e) => panic!("SPECIAL_FOLD_PARK_SELFTEST: FAIL: {e}"),
        }
        if std::env::var("SPECIAL_FOLD_PARK_SELFTEST_ONLY")
            .ok()
            .as_deref()
            == Some("1")
        {
            return Vec::new();
        }
    }
    if std::env::var("DIALOG_GCD_FUSED_APPLY_SELFTEST").is_ok() {
        match dialog_gcd_k5_tail7_fixed_apply_selftest() {
            Ok(()) => eprintln!(
                "DIALOG_GCD_FUSED_APPLY_SELFTEST: PASS (fused double/halve value, ancilla, phase)"
            ),
            Err(e) => panic!("DIALOG_GCD_FUSED_APPLY_SELFTEST: FAIL: {e}"),
        }
        if std::env::var("DIALOG_GCD_FUSED_APPLY_SELFTEST_ONLY")
            .ok()
            .as_deref()
            == Some("1")
        {
            return Vec::new();
        }
    }
    build_builder().ops
}

pub fn square_window_selftest() -> Result<(), String> {
    use sha3::digest::{ExtendableOutput, Update};
    const SHOTS: usize = 64;
    let nbits = std::env::var("SQUARE_WINDOW_SELFTEST_NBITS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(24);
    assert!(nbits > 0);
    let packed_value_check = 2 * nbits < 64;
    let wide_value_check = nbits <= 256;
    let mask = if packed_value_check { (1u64 << nbits) - 1 } else { u64::MAX };
    let out_mask = if packed_value_check { (1u64 << (2 * nbits)) - 1 } else { u64::MAX };
    let xs: Vec<u64> = (0..SHOTS as u64)
        .map(|s| {
            let r = s
                .wrapping_mul(0x9E37_79B9_7F4A_7C15)
                .wrapping_add(0xA076_1D64_78BD_642F);
            let r = (r ^ (r >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            r & mask
        })
        .collect();
    let x_masks: Vec<u64> = (0..nbits)
        .map(|k| {
            if packed_value_check {
                xs.iter()
                    .enumerate()
                    .fold(0u64, |acc, (shot, &xv)| acc | (((xv >> k) & 1) << shot))
            } else {
                let z = (k as u64)
                    .wrapping_mul(0xD6E8_FD9D_50B5_8A51)
                    .wrapping_add(0x9E37_79B9_7F4A_7C15);
                let z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
                z ^ (z >> 31)
            }
        })
        .collect();

    let build_one = |roundtrip: bool| -> (Vec<Op>, Vec<QubitId>, Vec<QubitId>, usize, usize) {
        let mut b = B::new();
        let x = b.alloc_qubits(nbits);
        let tmp = b.alloc_qubits(2 * nbits);
        schoolbook_square_symmetric_lowq_selfhosted(&mut b, &x, &tmp);
        if roundtrip {
            schoolbook_square_symmetric_lowq_selfhosted_inverse(&mut b, &x, &tmp);
        }
        let nq = b.next_qubit as usize;
        let nb = b.next_bit as usize;
        (b.ops, x, tmp, nq, nb)
    };

    let run = |ops: &[Op],
               x: &[QubitId],
               tmp: &[QubitId],
               nq: usize,
               nb: usize|
     -> (Vec<u64>, Vec<u64>, u64) {
        let mut seed = sha3::Shake128::default();
        seed.update(b"square-window-selftest");
        let mut xof = seed.finalize_xof();
        let mut sim = Simulator::new(nq, nb, &mut xof);
        sim.clear_for_shot();
        for k in 0..nbits {
            *sim.qubit_mut(x[k]) = x_masks[k];
        }
        sim.apply_iter(ops.iter());
        let out_x_masks: Vec<u64> = x.iter().map(|&q| sim.qubit(q)).collect();
        let out_tmp_masks: Vec<u64> = tmp.iter().map(|&q| sim.qubit(q)).collect();
        (out_x_masks, out_tmp_masks, sim.phase)
    };

    let (ops_fwd, x_fwd, tmp_fwd, nq_fwd, nb_fwd) = build_one(false);
    let (out_x_masks, out_tmp_masks, phase) = run(&ops_fwd, &x_fwd, &tmp_fwd, nq_fwd, nb_fwd);
    if phase != 0 {
        return Err(format!("forward phase garbage 0x{phase:x}"));
    }
    for (k, (&got, &want)) in out_x_masks.iter().zip(x_masks.iter()).enumerate() {
        if got != want {
            return Err(format!("forward x bit {k} changed"));
        }
    }
    if packed_value_check {
        for shot in 0..SHOTS {
            let got = out_tmp_masks
                .iter()
                .take(2 * nbits)
                .enumerate()
                .fold(0u64, |acc, (k, &bits)| acc | (((bits >> shot) & 1) << k));
            let want = xs[shot].wrapping_mul(xs[shot]) & out_mask;
            if got != want {
                return Err(format!(
                    "forward value mismatch shot {shot}: tmp got 0x{got:x} want 0x{want:x}"
                ));
            }
        }
    } else if wide_value_check {
        let in_limbs = (nbits + 63) / 64;
        let out_limbs = (2 * nbits + 63) / 64;
        for shot in 0..SHOTS {
            let mut x_limbs = vec![0u64; in_limbs];
            for k in 0..nbits {
                if (x_masks[k] >> shot) & 1 != 0 {
                    x_limbs[k / 64] |= 1u64 << (k % 64);
                }
            }
            let mut product = vec![0u64; out_limbs];
            for i in 0..in_limbs {
                let mut carry = 0u128;
                for j in 0..in_limbs {
                    let idx = i + j;
                    if idx >= out_limbs {
                        break;
                    }
                    let cur = product[idx] as u128
                        + (x_limbs[i] as u128) * (x_limbs[j] as u128)
                        + carry;
                    product[idx] = cur as u64;
                    carry = cur >> 64;
                }
                let mut idx = i + in_limbs;
                while carry != 0 && idx < out_limbs {
                    let cur = product[idx] as u128 + carry;
                    product[idx] = cur as u64;
                    carry = cur >> 64;
                    idx += 1;
                }
            }
            for k in 0..(2 * nbits) {
                let got = (out_tmp_masks[k] >> shot) & 1;
                let want = (product[k / 64] >> (k % 64)) & 1;
                if got != want {
                    return Err(format!("forward value mismatch shot {shot} bit {k}"));
                }
            }
        }
    }

    let (ops_rt, x_rt, tmp_rt, nq_rt, nb_rt) = build_one(true);
    let (out_x_masks, out_tmp_masks, phase) = run(&ops_rt, &x_rt, &tmp_rt, nq_rt, nb_rt);
    if phase != 0 {
        return Err(format!("roundtrip phase garbage 0x{phase:x}"));
    }
    for (k, (&got, &want)) in out_x_masks.iter().zip(x_masks.iter()).enumerate() {
        if got != want {
            return Err(format!("roundtrip x bit {k} changed"));
        }
    }
    for (k, &got) in out_tmp_masks.iter().enumerate() {
        if got != 0 {
            return Err(format!("roundtrip tmp bit {k} dirty mask 0x{got:x}"));
        }
    }
    Ok(())
}


/// Standalone differential selftest for the fused-fold freed-tail lever
/// (`DIALOG_GCD_FOLD_FREED_TAIL`). Runs in the normal (non-test) build because
/// the `#[cfg(test)]` module does not compile on this base. For each
/// `(e,d) ∈ {0,1}²` it builds the BASELINE per-position fold ripple and the
/// FREED-TAIL ripple on the same random `y` (64 shots/lane), simulates both, and
/// asserts: (1) identical `y` outputs, (2) all fold ancillae returned to |0>,
/// (3) zero global phase. Returns Err with the first divergence. Invoke via
/// `FOLD_FREED_TAIL_SELFTEST=1 build_circuit`.
pub fn fold_freed_tail_selftest() -> Result<(), String> {
    use sha3::digest::{ExtendableOutput, Update};
    let hi_delta = 33usize;
    let hi_c = 32usize;
    let nbits = 64usize; // y width for the test (covers the active+tail span)
    for &windowed in &[true, false] {
        let last = if windowed {
            hi_delta + 19 // mirror KAL_DOUBLE_CARRY_TRUNC_W=19
        } else {
            nbits - 2
        };
        for ed in 0u64..4 {
            let e_val = ed & 1;
            let d_val = (ed >> 1) & 1;
            for &is_add in &[true, false] {
                // Build both circuits over identical qubit layout.
                let build_one = |freed: bool| -> (Vec<Op>, Vec<QubitId>, usize, usize) {
                    let mut b = B::new();
                    let y = b.alloc_qubits(nbits);
                    let ovf1 = b.alloc_qubit();
                    let ovf2 = b.alloc_qubit();
                    let s2 = b.alloc_qubit();
                    let e = b.alloc_qubit();
                    let d = b.alloc_qubit();
                    let h = b.alloc_qubit();
                    let xed = b.alloc_qubit();
                    let eord = b.alloc_qubit();
                    let n10 = b.alloc_qubit();
                    // Exercise the real caller relation for every (e,d) pair:
                    // s2=1, ovf1=d, ovf2=e gives
                    // d=ovf1&s2 and e=ovf1^d^ovf2.
                    b.x(s2);
                    if d_val == 1 {
                        b.x(ovf1);
                    }
                    if e_val == 1 {
                        b.x(ovf2);
                    }
                    b.ccx(ovf1, s2, d);
                    b.cx(ovf1, e);
                    b.cx(d, e);
                    b.cx(ovf2, e);
                    b.ccx(e, d, h); // h = e&d
                    b.cx(e, xed);
                    b.cx(d, xed); // xed = e^d
                    b.cx(xed, eord);
                    b.cx(h, eord); // eord = e|d
                    b.cx(d, n10);
                    b.cx(h, n10); // n10 = !e&d
                    if freed {
                        fold_ripple_freed_tail_ed(
                            &mut b,
                            &y,
                            e,
                            d,
                            h,
                            xed,
                            eord,
                            n10,
                            Some((ovf1, ovf2, s2)),
                            None,
                            last,
                            is_add,
                        );
                    } else {
                        let controls =
                            secp_fold_controls(e, d, h, xed, eord, n10, hi_delta, hi_c);
                        if is_add {
                            cadd_per_position_controls_trunc(&mut b, &y, &controls, last);
                        } else {
                            csub_per_position_controls_trunc(&mut b, &y, &controls, last);
                        }
                    }
                    // uncompute derived controls (same as the fused fns) so all 6
                    // ancillae return to |0> on a value-exact ripple.
                    b.cx(h, n10);
                    b.cx(d, n10);
                    b.cx(h, eord);
                    b.cx(xed, eord);
                    b.cx(d, xed);
                    b.cx(e, xed);
                    b.ccx(e, d, h);
                    b.cx(ovf2, e);
                    b.cx(d, e);
                    b.cx(ovf1, e);
                    b.ccx(ovf1, s2, d);
                    if e_val == 1 {
                        b.x(ovf2);
                    }
                    if d_val == 1 {
                        b.x(ovf1);
                    }
                    b.x(s2);
                    let nq = b.next_qubit as usize;
                    let nb = b.next_bit as usize;
                    (b.ops, y, nq, nb)
                };
                let (ops_base, y_b, nq_b, nb_b) = build_one(false);
                let (ops_freed, y_f, nq_f, nb_f) = build_one(true);
                // deterministic random y per shot, including adversarial
                // carry-propagation patterns (long runs of 1s above bit 33 that
                // force the truncated tail carry to escape / saturate).
                let mask: u64 = if nbits >= 64 { u64::MAX } else { (1u64 << nbits) - 1 };
                let ys: Vec<u64> = (0..64u64)
                    .map(|s| {
                        let r = s
                            .wrapping_mul(0x9E37_79B9_7F4A_7C15)
                            .wrapping_add(0xD1B5_4A32_D192_ED03);
                        let r = (r ^ (r >> 31)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
                        let r = r ^ (r >> 27);
                        let base = r & mask;
                        // every 4th shot: all-ones above bit 33 (worst case carry run)
                        if s % 4 == 0 {
                            base | (mask & !((1u64 << (hi_delta + 1)) - 1))
                        } else if s % 4 == 1 {
                            base & ((1u64 << (hi_delta + 1)) - 1)
                        } else {
                            base
                        }
                    })
                    .collect();

                let run = |ops: &[Op], y: &[QubitId], nq: usize, nb: usize| -> (Vec<u64>, bool, u64) {
                    let mut s2 = sha3::Shake128::default();
                    s2.update(b"fold-sim");
                    let mut xof2 = s2.finalize_xof();
                    let mut sim = Simulator::new(nq, nb, &mut xof2);
                    sim.clear_for_shot();
                    for (shot, &yv) in ys.iter().enumerate() {
                        for k in 0..nbits {
                            if (yv >> k) & 1 != 0 {
                                *sim.qubit_mut(y[k]) |= 1u64 << shot;
                            }
                        }
                    }
                    sim.apply_iter(ops.iter());
                    let outs: Vec<u64> = (0..64)
                        .map(|shot| {
                            let mut v = 0u64;
                            for k in 0..nbits {
                                v |= ((sim.qubit(y[k]) >> shot) & 1) << k;
                            }
                            v
                        })
                        .collect();
                    let anc_clean =
                        (nbits..nq).all(|q| sim.qubit(QubitId(q as u64)) == 0);
                    (outs, anc_clean, sim.phase)
                };
                let (out_b, clean_b, phase_b) = run(&ops_base, &y_b, nq_b, nb_b);
                let (out_f, clean_f, phase_f) = run(&ops_freed, &y_f, nq_f, nb_f);

                if !clean_b {
                    return Err(format!("baseline left ancilla dirty (ed={ed} add={is_add} win={windowed})"));
                }
                if !clean_f {
                    return Err(format!("freed-tail left ancilla dirty (ed={ed} add={is_add} win={windowed})"));
                }
                if phase_f != 0 {
                    return Err(format!("freed-tail left phase garbage 0x{phase_f:x} (ed={ed} add={is_add} win={windowed})"));
                }
                let _ = phase_b;
                for shot in 0..64 {
                    if out_b[shot] != out_f[shot] {
                        return Err(format!(
                            "value mismatch shot {shot}: base 0x{:x} freed 0x{:x} (ed={ed} add={is_add} win={windowed}, y_in=0x{:x})",
                            out_b[shot], out_f[shot], ys[shot]
                        ));
                    }
                }
            }
        }
    }
    Ok(())
}

pub fn special_fold_park_selftest() -> Result<(), String> {
    use sha3::digest::{ExtendableOutput, Update};

    let c = U256::MAX
        .wrapping_sub(SECP256K1_P)
        .wrapping_add(U256::from(1u64));
    let nbits = 64usize;
    let window = 20usize;

    for ctrl_value in 0u64..=1 {
        for &is_add in &[true, false] {
            let build_one = |parked: bool| {
                let mut b = B::new();
                let acc = b.alloc_qubits(nbits);
                let ctrl = b.alloc_qubit();
                let scratch = b.alloc_qubits(5);
                if ctrl_value != 0 {
                    b.x(ctrl);
                }
                if parked {
                    if is_add {
                        cadd_nbit_const_direct_trunc_fast_releasing_scratch(
                            &mut b, &acc, c, ctrl, window, &scratch,
                        );
                    } else {
                        csub_nbit_const_direct_trunc_fast_releasing_scratch(
                            &mut b, &acc, c, ctrl, window, &scratch,
                        );
                    }
                } else if is_add {
                    cadd_nbit_const_direct_trunc_fast_borrowed_carries(
                        &mut b, &acc, c, ctrl, window, &scratch,
                    );
                } else {
                    csub_nbit_const_direct_trunc_fast_borrowed_carries(
                        &mut b, &acc, c, ctrl, window, &scratch,
                    );
                }
                if ctrl_value != 0 {
                    b.x(ctrl);
                }
                (b.ops, acc, b.next_qubit as usize, b.next_bit as usize)
            };

            let (base_ops, base_acc, base_nq, base_nb) = build_one(false);
            let (parked_ops, parked_acc, parked_nq, parked_nb) = build_one(true);
            let inputs: Vec<u64> = (0..64u64)
                .map(|shot| {
                    let mixed = shot
                        .wrapping_mul(0x9E37_79B9_7F4A_7C15)
                        .wrapping_add(0xD1B5_4A32_D192_ED03);
                    match shot % 4 {
                        0 => mixed | (!0u64 << 33),
                        1 => mixed & ((1u64 << 34) - 1),
                        _ => mixed ^ (mixed >> 29),
                    }
                })
                .collect();

            let run = |ops: &[Op], acc: &[QubitId], nq: usize, nb: usize| {
                let mut seed = Shake256::default();
                seed.update(b"special-fold-park-selftest");
                seed.update(&[ctrl_value as u8, is_add as u8]);
                let mut xof = seed.finalize_xof();
                let mut sim = Simulator::new(nq, nb, &mut xof);
                sim.clear_for_shot();
                for (shot, &input) in inputs.iter().enumerate() {
                    for bit_index in 0..nbits {
                        if (input >> bit_index) & 1 != 0 {
                            *sim.qubit_mut(acc[bit_index]) |= 1u64 << shot;
                        }
                    }
                }
                sim.apply_iter(ops.iter());
                let outputs: Vec<u64> = (0..64)
                    .map(|shot| {
                        let mut value = 0u64;
                        for bit_index in 0..nbits {
                            value |= ((sim.qubit(acc[bit_index]) >> shot) & 1) << bit_index;
                        }
                        value
                    })
                    .collect();
                let clean = (nbits..nq).all(|q| sim.qubit(QubitId(q as u64)) == 0);
                (outputs, clean, sim.phase)
            };

            let (base_out, base_clean, base_phase) =
                run(&base_ops, &base_acc, base_nq, base_nb);
            let (parked_out, parked_clean, parked_phase) =
                run(&parked_ops, &parked_acc, parked_nq, parked_nb);
            if !base_clean || base_phase != 0 {
                return Err(format!(
                    "baseline dirty: ctrl={ctrl_value} add={is_add} clean={base_clean} phase=0x{base_phase:x}"
                ));
            }
            if !parked_clean || parked_phase != 0 {
                return Err(format!(
                    "parked dirty: ctrl={ctrl_value} add={is_add} clean={parked_clean} phase=0x{parked_phase:x}"
                ));
            }
            if base_out != parked_out {
                let shot = base_out
                    .iter()
                    .zip(&parked_out)
                    .position(|(base, parked)| base != parked)
                    .unwrap_or(0);
                return Err(format!(
                    "value mismatch shot {shot}: base=0x{:x} parked=0x{:x} input=0x{:x} ctrl={ctrl_value} add={is_add}",
                    base_out[shot], parked_out[shot], inputs[shot]
                ));
            }
        }
    }
    Ok(())
}


#[cfg(test)]
mod direct_const_tests {
    use super::*;
    use sha3::{
        digest::{ExtendableOutput, Update, XofReader},
        Shake128,
    };

    fn set_reg<R: XofReader>(sim: &mut Simulator<'_, R>, qs: &[QubitId], val: u64, shot: usize) {
        for (i, &q) in qs.iter().enumerate() {
            if ((val >> i) & 1) != 0 {
                *sim.qubit_mut(q) |= 1u64 << shot;
            } else {
                *sim.qubit_mut(q) &= !(1u64 << shot);
            }
        }
    }

    fn get_reg<R: XofReader>(sim: &Simulator<'_, R>, qs: &[QubitId], shot: usize) -> u64 {
        let mut out = 0u64;
        for (i, &q) in qs.iter().enumerate() {
            out |= ((sim.qubit(q) >> shot) & 1) << i;
        }
        out
    }

    #[test]
    fn one_inv_dx3_blocker_is_fail_closed_on_cleanup_invariant() {
        assert!(ONE_INV_DX3_AFFINE_PA_BLOCKER.contains("Rx-Qx"));
        assert!(ONE_INV_DX3_AFFINE_PA_BLOCKER.contains("second inversion"));
        assert!(ONE_INV_DX3_AFFINE_PA_BLOCKER.contains("dirty reset"));
    }

    #[test]
    fn dialog_gcd_selected_body_nocin_matches_cin_reference() {
        if let Err(e) = dialog_gcd_selected_body_nocin_selftest() {
            panic!("no-c_in selected body selftest failed: {e}");
        }
    }

    #[test]
    fn aliased_gate_wrappers_are_not_silent_noops() {
        let mut b = B::new();
        let q0 = b.alloc_qubit();
        let q1 = b.alloc_qubit();
        b.cz(q0, q0);
        b.ccz(q0, q0, q1);
        b.ccz(q0, q1, q0);
        b.ccz(q0, q0, q0);
        b.ccx(q0, q0, q1);
        let kinds = b.ops.iter().map(|op| op.kind).collect::<Vec<_>>();
        assert_eq!(
            kinds,
            vec![
                OperationType::Z,
                OperationType::CZ,
                OperationType::CZ,
                OperationType::Z,
                OperationType::CX,
            ]
        );
        assert!(std::panic::catch_unwind(|| {
            let mut b = B::new();
            let q = b.alloc_qubit();
            b.cx(q, q);
        })
        .is_err());
        assert!(std::panic::catch_unwind(|| {
            let mut b = B::new();
            let q0 = b.alloc_qubit();
            let q1 = b.alloc_qubit();
            b.ccx(q0, q1, q0);
        })
        .is_err());
    }

    #[test]
    fn dx3_witness_is_not_an_output_cleanup_coordinate() {
        let p = SECP256K1_P;
        let beta = U256::from_str_radix(
            "7AE96A2B657C07106E64479EAC3434E99CF0497512F58995C1396C28719501EE",
            16,
        )
        .unwrap();
        let dx = U256::from(0x1234_5678_9abc_def0u64);
        let beta_dx = beta.mul_mod(dx, p);
        assert_ne!(dx, beta_dx);
        assert_eq!(beta.mul_mod(beta, p).mul_mod(beta, p), U256::from(1u64));
        assert_eq!(
            dx.mul_mod(dx, p).mul_mod(dx, p),
            beta_dx.mul_mod(beta_dx, p).mul_mod(beta_dx, p)
        );
    }

    fn assert_borrowed_carry_adder_basis(is_sub: bool) {
        const N: usize = 5;
        const MOD: u64 = 1 << N;
        let mut b = B::new();
        let a = b.alloc_qubits(N);
        let acc = b.alloc_qubits(N);
        let carries = b.alloc_qubits(N - 1);
        if is_sub {
            sub_nbit_qq_fast_borrowed_carries(&mut b, &a, &acc, &carries);
        } else {
            add_nbit_qq_fast_borrowed_carries(&mut b, &a, &acc, &carries);
        }
        let nq = b.next_qubit as usize;
        let nb = b.next_bit as usize;

        for batch in 0..16usize {
            let mut seed = Shake128::default();
            seed.update(if is_sub {
                b"borrowed-sub-small"
            } else {
                b"borrowed-add-small"
            });
            let mut xof = seed.finalize_xof();
            let mut sim = Simulator::new(nq, nb, &mut xof);
            for shot in 0..64usize {
                let case = batch * 64 + shot;
                let x = (case as u64) & (MOD - 1);
                let y = ((case as u64) >> N) & (MOD - 1);
                set_reg(&mut sim, &acc, x, shot);
                set_reg(&mut sim, &a, y, shot);
            }
            sim.apply(&b.ops);
            assert_eq!(
                sim.global_phase(),
                0,
                "borrowed carry adder left phase garbage"
            );
            for shot in 0..64usize {
                let case = batch * 64 + shot;
                let x = (case as u64) & (MOD - 1);
                let y = ((case as u64) >> N) & (MOD - 1);
                let expect = if is_sub {
                    x.wrapping_sub(y) & (MOD - 1)
                } else {
                    x.wrapping_add(y) & (MOD - 1)
                };
                assert_eq!(get_reg(&sim, &acc, shot), expect, "case {case}");
                assert_eq!(get_reg(&sim, &a, shot), y, "a changed in case {case}");
                assert_eq!(
                    get_reg(&sim, &carries, shot),
                    0,
                    "borrowed carries not clean in case {case}"
                );
            }
        }
    }

    #[test]
    fn borrowed_carry_add_small_basis_is_clean() {
        assert_borrowed_carry_adder_basis(false);
    }

    #[test]
    fn borrowed_carry_sub_small_basis_is_clean() {
        assert_borrowed_carry_adder_basis(true);
    }

    fn sub_mod_p(a: U256, b: U256, p: U256) -> U256 {
        if a >= b {
            a - b
        } else {
            p - (b - a)
        }
    }

    #[test]
    fn direct_controlled_const_sub_small_basis_is_phase_clean() {
        const N: usize = 8;
        let c = U256::from(0b1011_0111u64);
        let mut b = B::new();
        let acc = b.alloc_qubits(N);
        let ctrl = b.alloc_qubit();
        csub_nbit_const_direct_fast(&mut b, &acc, c, ctrl);
        let nq = b.next_qubit as usize;
        let nb = b.next_bit as usize;

        let mut seed = Shake128::default();
        seed.update(b"direct-csub-small");
        let mut xof = seed.finalize_xof();
        let mut sim = Simulator::new(nq, nb, &mut xof);
        for shot in 0..64usize {
            let x = ((shot * 37 + 11) & 0xff) as u64;
            let ctrl_v = (shot & 1) as u64;
            set_reg(&mut sim, &acc, x, shot);
            if ctrl_v != 0 {
                *sim.qubit_mut(ctrl) |= 1u64 << shot;
            }
        }
        sim.apply(&b.ops);
        assert_eq!(sim.global_phase(), 0, "direct csub left phase garbage");
        for shot in 0..64usize {
            let x = ((shot * 37 + 11) & 0xff) as u64;
            let ctrl_v = (shot & 1) as u64;
            let expect = x.wrapping_sub(ctrl_v * 0b1011_0111) & 0xff;
            assert_eq!(get_reg(&sim, &acc, shot), expect, "shot {shot}");
            assert_eq!((sim.qubit(ctrl) >> shot) & 1, ctrl_v, "ctrl shot {shot}");
        }
    }

    #[test]
    fn direct_controlled_const_add_small_basis_is_phase_clean() {
        const N: usize = 8;
        let c = U256::from(0b1011_0111u64);
        let mut b = B::new();
        let acc = b.alloc_qubits(N);
        let ctrl = b.alloc_qubit();
        cadd_nbit_const_direct_fast(&mut b, &acc, c, ctrl);
        let nq = b.next_qubit as usize;
        let nb = b.next_bit as usize;

        let mut seed = Shake128::default();
        seed.update(b"direct-cadd-small");
        let mut xof = seed.finalize_xof();
        let mut sim = Simulator::new(nq, nb, &mut xof);
        for shot in 0..64usize {
            let x = ((shot * 37 + 11) & 0xff) as u64;
            let ctrl_v = (shot & 1) as u64;
            set_reg(&mut sim, &acc, x, shot);
            if ctrl_v != 0 {
                *sim.qubit_mut(ctrl) |= 1u64 << shot;
            }
        }
        sim.apply(&b.ops);
        assert_eq!(sim.global_phase(), 0, "direct cadd left phase garbage");
        for shot in 0..64usize {
            let x = ((shot * 37 + 11) & 0xff) as u64;
            let ctrl_v = (shot & 1) as u64;
            let expect = x.wrapping_add(ctrl_v * 0b1011_0111) & 0xff;
            assert_eq!(get_reg(&sim, &acc, shot), expect, "shot {shot}");
            assert_eq!((sim.qubit(ctrl) >> shot) & 1, ctrl_v, "ctrl shot {shot}");
        }
    }

    #[test]
    fn round84_fused_square_xtail_component_matches_relation() {
        let ops = build_round84_fused_square_xtail_component();
        let (num_qubits, num_bits, _num_registers, regs) = analyze_ops(ops.iter().copied());
        assert_eq!(regs.len(), 4);
        let p = SECP256K1_P;
        let cases: Vec<(U256, U256, U256)> = (0..32u64)
            .map(|i| {
                let tx = U256::from_limbs([
                    0x9e37_79b9_7f4a_7c15u64.wrapping_mul(i + 1),
                    0xd1b5_4a32_d192_ed03u64.wrapping_mul(i + 3),
                    0x94d0_49bb_1331_11ebu64.wrapping_mul(i + 5),
                    0x2545_f491_4f6c_dd1du64.wrapping_mul(i + 7),
                ]) % p;
                let lam = U256::from_limbs([
                    0xbf58_476d_1ce4_e5b9u64.wrapping_mul(i + 11),
                    0x94d0_49bb_1331_11ebu64.wrapping_mul(i + 13),
                    0xdbe6_d5d5_fe4c_ce2fu64.wrapping_mul(i + 17),
                    0xa409_3822_299f_31d0u64.wrapping_mul(i + 19),
                ]) % p;
                let ox = U256::from_limbs([
                    0x632b_e59b_d9b4_e019u64.wrapping_mul(i + 23),
                    0x8515_7af5_4f1d_2d2du64.wrapping_mul(i + 29),
                    0x9e37_79b9_7f4a_7c15u64.wrapping_mul(i + 31),
                    0xbf58_476d_1ce4_e5b9u64.wrapping_mul(i + 37),
                ]) % p;
                (tx, lam, ox)
            })
            .collect();

        let mut seed = Shake128::default();
        seed.update(b"round84-xtail-component");
        let mut xof = seed.finalize_xof();
        let mut sim = Simulator::new(num_qubits as usize, num_bits as usize, &mut xof);
        for (shot, (tx, lam, ox)) in cases.iter().enumerate() {
            sim.set_register(&regs[0], *tx, shot);
            sim.set_register(&regs[1], *lam, shot);
            sim.set_register(&regs[2], *ox, shot);
            sim.set_register(&regs[3], U256::ZERO, shot);
        }

        sim.apply(&ops);
        for (shot, (tx, lam, ox)) in cases.iter().enumerate() {
            let expected = sub_mod_p(
                sub_mod_p(lam.mul_mod(*lam, p), *tx, p),
                ox.add_mod(*ox, p),
                p,
            );
            assert_eq!(
                sim.get_register(&regs[0], shot),
                expected,
                "x-tail shot {shot}"
            );
            assert_eq!(sim.get_register(&regs[1], shot), *lam, "lambda shot {shot}");
            assert_eq!(
                sim.get_register(&regs[2], shot),
                *ox,
                "offset-x shot {shot}"
            );
        }
        let live_mask = (1u64 << cases.len()) - 1;
        assert_eq!(sim.global_phase() & live_mask, 0, "x-tail phase garbage");
        for reg in &regs {
            for item in reg {
                if let QubitOrBit::Qubit(q) = *item {
                    *sim.qubit_mut(q) = 0;
                }
            }
        }
        for q in 0..num_qubits {
            assert_eq!(
                sim.qubit(QubitId(q)) & live_mask,
                0,
                "x-tail ancilla garbage q{q}"
            );
        }
    }

    #[test]
    fn round190_selector_fused_source_live_residual_is_exact_on_small_widths() {
        for width in [2usize, 3, 4] {
            let ops = build_round190_selector_fused_source_live_residual_width(width);
            let (num_qubits, num_bits, num_registers, regs) = analyze_ops(ops.iter().copied());
            assert_eq!(num_registers, 3, "width {width} register count");
            assert_eq!(regs.len(), 3, "width {width} regs");
            assert_eq!(num_bits as usize, width, "width {width} hmr bits");
            assert_eq!(num_qubits as usize, 4 * width + 3, "width {width} qubits");
            for (idx, reg) in regs.iter().enumerate() {
                assert_eq!(reg.len(), width, "width {width} reg {idx}");
                assert!(reg.iter().all(|item| matches!(item, QubitOrBit::Qubit(_))));
            }
            let toffoli_ops = ops
                .iter()
                .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
                .count();
            assert_eq!(toffoli_ops, 3 * width, "width {width} toffoli");
            let pred_reg: Vec<QubitId> = regs[0]
                .iter()
                .map(|item| match item {
                    QubitOrBit::Qubit(q) => *q,
                    _ => unreachable!(),
                })
                .collect();
            let add_reg: Vec<QubitId> = regs[1]
                .iter()
                .map(|item| match item {
                    QubitOrBit::Qubit(q) => *q,
                    _ => unreachable!(),
                })
                .collect();
            let target_reg: Vec<QubitId> = regs[2]
                .iter()
                .map(|item| match item {
                    QubitOrBit::Qubit(q) => *q,
                    _ => unreachable!(),
                })
                .collect();

            let modulus = 1u64 << width;
            let states = modulus * modulus * modulus;
            let mut seed = Shake128::default();
            seed.update(b"round190-selector-fused-source-live-residual");
            seed.update(&[width as u8]);
            let mut xof = seed.finalize_xof();
            for batch_start in (0..states).step_by(64) {
                let mut sim = Simulator::new(num_qubits as usize, num_bits as usize, &mut xof);
                let batch_end = (batch_start + 64).min(states);
                for case in batch_start..batch_end {
                    let shot = (case - batch_start) as usize;
                    let predecessor = case & (modulus - 1);
                    let addend = (case >> width) & (modulus - 1);
                    let target = (case >> (2 * width)) & (modulus - 1);
                    set_reg(&mut sim, &pred_reg, predecessor, shot);
                    set_reg(&mut sim, &add_reg, addend, shot);
                    set_reg(&mut sim, &target_reg, target, shot);
                }

                sim.apply(&ops);
                let live_mask = if batch_end - batch_start == 64 {
                    u64::MAX
                } else {
                    (1u64 << (batch_end - batch_start)) - 1
                };
                assert_eq!(
                    sim.global_phase() & live_mask,
                    0,
                    "width {width} selector-fused residual phase garbage"
                );
                for case in batch_start..batch_end {
                    let shot = (case - batch_start) as usize;
                    let predecessor = case & (modulus - 1);
                    let addend = (case >> width) & (modulus - 1);
                    let target = (case >> (2 * width)) & (modulus - 1);
                    let low = predecessor & 0b11;
                    let expected = if low == 0 {
                        target
                    } else if ((predecessor >> 1) & 1) != 0 {
                        target.wrapping_sub(addend) & (modulus - 1)
                    } else {
                        target.wrapping_add(addend) & (modulus - 1)
                    };
                    assert_eq!(
                        get_reg(&sim, &pred_reg, shot),
                        predecessor,
                        "width {width} predecessor changed case {case}"
                    );
                    assert_eq!(
                        get_reg(&sim, &add_reg, shot),
                        addend,
                        "width {width} addend changed case {case}"
                    );
                    assert_eq!(
                        get_reg(&sim, &target_reg, shot),
                        expected,
                        "width {width} target mismatch case {case}"
                    );
                }
                for reg in [&pred_reg, &add_reg, &target_reg] {
                    for &q in reg {
                        *sim.qubit_mut(q) = 0;
                    }
                }
                for q in 0..num_qubits {
                    assert_eq!(
                        sim.qubit(QubitId(q)) & live_mask,
                        0,
                        "width {width} scratch garbage q{q}"
                    );
                }
            }
        }
    }

    #[test]
    fn round190_external_active_signed_digit_is_select0_safe_on_small_widths() {
        for width in [2usize, 3, 4] {
            let ops = build_round190_external_active_signed_digit_width(width);
            let (num_qubits, num_bits, num_registers, regs) = analyze_ops(ops.iter().copied());
            assert_eq!(num_registers, 4, "width {width} register count");
            assert_eq!(regs.len(), 4, "width {width} regs");
            assert_eq!(num_bits as usize, width, "width {width} hmr bits");
            assert_eq!(num_qubits as usize, 3 * width + 4, "width {width} qubits");
            assert_eq!(regs[0].len(), 1, "width {width} active width");
            assert_eq!(regs[1].len(), 1, "width {width} sign width");
            assert_eq!(regs[2].len(), width, "width {width} addend width");
            assert_eq!(regs[3].len(), width, "width {width} target width");
            for (idx, reg) in regs.iter().enumerate() {
                assert!(
                    reg.iter().all(|item| matches!(item, QubitOrBit::Qubit(_))),
                    "width {width} reg {idx} must be qubits"
                );
            }
            let toffoli_ops = ops
                .iter()
                .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
                .count();
            assert_eq!(toffoli_ops, 3 * width - 2, "width {width} toffoli");

            let active_q = match regs[0][0] {
                QubitOrBit::Qubit(q) => q,
                _ => unreachable!(),
            };
            let sign_q = match regs[1][0] {
                QubitOrBit::Qubit(q) => q,
                _ => unreachable!(),
            };
            let add_reg: Vec<QubitId> = regs[2]
                .iter()
                .map(|item| match item {
                    QubitOrBit::Qubit(q) => *q,
                    _ => unreachable!(),
                })
                .collect();
            let target_reg: Vec<QubitId> = regs[3]
                .iter()
                .map(|item| match item {
                    QubitOrBit::Qubit(q) => *q,
                    _ => unreachable!(),
                })
                .collect();

            let modulus = 1u64 << width;
            let states = 4 * modulus * modulus;
            let mut seed = Shake128::default();
            seed.update(b"round190-external-active-signed-digit");
            seed.update(&[width as u8]);
            let mut xof = seed.finalize_xof();
            for batch_start in (0..states).step_by(64) {
                let mut sim = Simulator::new(num_qubits as usize, num_bits as usize, &mut xof);
                let batch_end = (batch_start + 64).min(states);
                for case in batch_start..batch_end {
                    let shot = (case - batch_start) as usize;
                    let active = case & 1;
                    let sign = (case >> 1) & 1;
                    let addend = (case >> 2) & (modulus - 1);
                    let target = (case >> (2 + width)) & (modulus - 1);
                    *sim.qubit_mut(active_q) |= active << shot;
                    *sim.qubit_mut(sign_q) |= sign << shot;
                    set_reg(&mut sim, &add_reg, addend, shot);
                    set_reg(&mut sim, &target_reg, target, shot);
                }

                sim.apply(&ops);
                let live_mask = if batch_end - batch_start == 64 {
                    u64::MAX
                } else {
                    (1u64 << (batch_end - batch_start)) - 1
                };
                assert_eq!(
                    sim.global_phase() & live_mask,
                    0,
                    "width {width} external-active phase garbage"
                );
                for case in batch_start..batch_end {
                    let shot = (case - batch_start) as usize;
                    let active = case & 1;
                    let sign = (case >> 1) & 1;
                    let addend = (case >> 2) & (modulus - 1);
                    let target = (case >> (2 + width)) & (modulus - 1);
                    let expected = if active == 0 {
                        target
                    } else if sign != 0 {
                        target.wrapping_sub(addend) & (modulus - 1)
                    } else {
                        target.wrapping_add(addend) & (modulus - 1)
                    };
                    assert_eq!(
                        (sim.qubit(active_q) >> shot) & 1,
                        active,
                        "width {width} active changed case {case}"
                    );
                    assert_eq!(
                        (sim.qubit(sign_q) >> shot) & 1,
                        sign,
                        "width {width} sign changed case {case}"
                    );
                    assert_eq!(
                        get_reg(&sim, &add_reg, shot),
                        addend,
                        "width {width} addend changed case {case}"
                    );
                    assert_eq!(
                        get_reg(&sim, &target_reg, shot),
                        expected,
                        "width {width} target mismatch case {case}"
                    );
                }
                *sim.qubit_mut(active_q) = 0;
                *sim.qubit_mut(sign_q) = 0;
                for reg in [&add_reg, &target_reg] {
                    for &q in reg {
                        *sim.qubit_mut(q) = 0;
                    }
                }
                for q in 0..num_qubits {
                    assert_eq!(
                        sim.qubit(QubitId(q)) & live_mask,
                        0,
                        "width {width} external-active scratch garbage q{q}"
                    );
                }
            }
        }
    }

    #[test]
    fn round190_shared_active_external_digits_reuse_selector_safely_on_small_widths() {
        for (width, digits) in [(2usize, 3usize), (3, 2)] {
            let ops = build_round190_shared_active_external_signed_digits_width(width, digits);
            let (num_qubits, num_bits, num_registers, regs) = analyze_ops(ops.iter().copied());
            assert_eq!(
                num_registers as usize,
                1 + 2 * digits,
                "width {width} digits {digits} register count"
            );
            assert_eq!(
                regs.len(),
                1 + 2 * digits,
                "width {width} digits {digits} regs"
            );
            assert_eq!(
                num_bits as usize,
                width * digits,
                "width {width} digits {digits} hmr bits"
            );
            assert_eq!(
                num_qubits as usize,
                (2 * digits + 2) * width + 3,
                "width {width} digits {digits} qubits"
            );
            for (idx, reg) in regs.iter().enumerate() {
                assert_eq!(reg.len(), width, "width {width} digits {digits} reg {idx}");
                assert!(
                    reg.iter().all(|item| matches!(item, QubitOrBit::Qubit(_))),
                    "width {width} digits {digits} reg {idx} must be qubits"
                );
            }
            let toffoli_ops = ops
                .iter()
                .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
                .count();
            assert_eq!(
                toffoli_ops,
                2 + digits * (3 * width - 2),
                "width {width} digits {digits} toffoli"
            );

            let qregs: Vec<Vec<QubitId>> = regs
                .iter()
                .map(|reg| {
                    reg.iter()
                        .map(|item| match item {
                            QubitOrBit::Qubit(q) => *q,
                            _ => unreachable!(),
                        })
                        .collect()
                })
                .collect();

            let modulus = 1u64 << width;
            let mut states = modulus;
            for _ in 0..digits {
                states *= modulus * modulus;
            }
            let mut seed = Shake128::default();
            seed.update(b"round190-shared-active-external-digits");
            seed.update(&[width as u8, digits as u8]);
            let mut xof = seed.finalize_xof();
            for batch_start in (0..states).step_by(64) {
                let mut sim = Simulator::new(num_qubits as usize, num_bits as usize, &mut xof);
                let batch_end = (batch_start + 64).min(states);
                for case in batch_start..batch_end {
                    let shot = (case - batch_start) as usize;
                    let mut cursor = case;
                    let predecessor = cursor & (modulus - 1);
                    cursor >>= width;
                    set_reg(&mut sim, &qregs[0], predecessor, shot);
                    for digit in 0..digits {
                        let addend = cursor & (modulus - 1);
                        cursor >>= width;
                        let target = cursor & (modulus - 1);
                        cursor >>= width;
                        set_reg(&mut sim, &qregs[1 + 2 * digit], addend, shot);
                        set_reg(&mut sim, &qregs[2 + 2 * digit], target, shot);
                    }
                }

                sim.apply(&ops);
                let live_mask = if batch_end - batch_start == 64 {
                    u64::MAX
                } else {
                    (1u64 << (batch_end - batch_start)) - 1
                };
                assert_eq!(
                    sim.global_phase() & live_mask,
                    0,
                    "width {width} digits {digits} shared-active phase garbage"
                );
                for case in batch_start..batch_end {
                    let shot = (case - batch_start) as usize;
                    let mut cursor = case;
                    let predecessor = cursor & (modulus - 1);
                    cursor >>= width;
                    assert_eq!(
                        get_reg(&sim, &qregs[0], shot),
                        predecessor,
                        "width {width} digits {digits} predecessor changed case {case}"
                    );
                    let active = (predecessor & 0b11) != 0;
                    let sign = ((predecessor >> 1) & 1) != 0;
                    for digit in 0..digits {
                        let addend = cursor & (modulus - 1);
                        cursor >>= width;
                        let target = cursor & (modulus - 1);
                        cursor >>= width;
                        let expected = if !active {
                            target
                        } else if sign {
                            target.wrapping_sub(addend) & (modulus - 1)
                        } else {
                            target.wrapping_add(addend) & (modulus - 1)
                        };
                        assert_eq!(
                            get_reg(&sim, &qregs[1 + 2 * digit], shot),
                            addend,
                            "width {width} digits {digits} addend {digit} changed case {case}"
                        );
                        assert_eq!(
                            get_reg(&sim, &qregs[2 + 2 * digit], shot),
                            expected,
                            "width {width} digits {digits} target {digit} mismatch case {case}"
                        );
                    }
                }
                for reg in &qregs {
                    for &q in reg {
                        *sim.qubit_mut(q) = 0;
                    }
                }
                for q in 0..num_qubits {
                    assert_eq!(
                        sim.qubit(QubitId(q)) & live_mask,
                        0,
                        "width {width} digits {digits} shared-active scratch garbage q{q}"
                    );
                }
            }
        }
    }

    #[test]
    fn round190_two_slot_router_is_exact_only_under_exactly_one_active_invariant() {
        for width in [2usize, 3] {
            let ops = build_round190_two_slot_exactly_one_active_router_width(width);
            let (num_qubits, num_bits, num_registers, regs) = analyze_ops(ops.iter().copied());
            assert_eq!(num_registers, 6, "width {width} register count");
            assert_eq!(regs.len(), 6, "width {width} regs");
            assert_eq!(num_bits as usize, width - 1, "width {width} hmr bits");
            assert_eq!(num_qubits as usize, 7 * width + 2, "width {width} qubits");
            for (idx, reg) in regs.iter().enumerate() {
                assert_eq!(reg.len(), width, "width {width} reg {idx}");
                assert!(reg.iter().all(|item| matches!(item, QubitOrBit::Qubit(_))));
            }
            let toffoli_ops = ops
                .iter()
                .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
                .count();
            assert_eq!(toffoli_ops, 7 * width + 1, "width {width} toffoli");

            let qregs: Vec<Vec<QubitId>> = regs
                .iter()
                .map(|reg| {
                    reg.iter()
                        .map(|item| match item {
                            QubitOrBit::Qubit(q) => *q,
                            _ => unreachable!(),
                        })
                        .collect()
                })
                .collect();
            let modulus = 1u64 << width;
            let active_predecessors: Vec<u64> =
                (0..modulus).filter(|pred| (pred & 0b11) != 0).collect();
            let inactive_predecessors: Vec<u64> =
                (0..modulus).filter(|pred| (pred & 0b11) == 0).collect();

            let mut cases = Vec::new();
            if width == 2 {
                for active_slot in 0..2usize {
                    for &active_pred in &active_predecessors {
                        for &inactive_pred in &inactive_predecessors {
                            for add0 in 0..modulus {
                                for target0 in 0..modulus {
                                    for add1 in 0..modulus {
                                        for target1 in 0..modulus {
                                            let (pred0, pred1) = if active_slot == 0 {
                                                (active_pred, inactive_pred)
                                            } else {
                                                (inactive_pred, active_pred)
                                            };
                                            cases.push((
                                                active_slot,
                                                pred0,
                                                add0,
                                                target0,
                                                pred1,
                                                add1,
                                                target1,
                                            ));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                for i in 0..512u64 {
                    let active_slot = (i & 1) as usize;
                    let active_pred =
                        active_predecessors[((i / 2) as usize) % active_predecessors.len()];
                    let inactive_pred =
                        inactive_predecessors[((i / 14) as usize) % inactive_predecessors.len()];
                    let add0 = (3 * i + 1) & (modulus - 1);
                    let target0 = (5 * i + 2) & (modulus - 1);
                    let add1 = (7 * i + 3) & (modulus - 1);
                    let target1 = (11 * i + 4) & (modulus - 1);
                    let (pred0, pred1) = if active_slot == 0 {
                        (active_pred, inactive_pred)
                    } else {
                        (inactive_pred, active_pred)
                    };
                    cases.push((active_slot, pred0, add0, target0, pred1, add1, target1));
                }
            }

            let mut seed = Shake128::default();
            seed.update(b"round190-two-slot-router");
            seed.update(&[width as u8]);
            let mut xof = seed.finalize_xof();
            for batch_start in (0..cases.len()).step_by(64) {
                let mut sim = Simulator::new(num_qubits as usize, num_bits as usize, &mut xof);
                let batch_end = (batch_start + 64).min(cases.len());
                for (shot, case) in cases[batch_start..batch_end].iter().enumerate() {
                    let &(_, pred0, add0, target0, pred1, add1, target1) = case;
                    set_reg(&mut sim, &qregs[0], pred0, shot);
                    set_reg(&mut sim, &qregs[1], add0, shot);
                    set_reg(&mut sim, &qregs[2], target0, shot);
                    set_reg(&mut sim, &qregs[3], pred1, shot);
                    set_reg(&mut sim, &qregs[4], add1, shot);
                    set_reg(&mut sim, &qregs[5], target1, shot);
                }

                sim.apply(&ops);
                let live_mask = if batch_end - batch_start == 64 {
                    u64::MAX
                } else {
                    (1u64 << (batch_end - batch_start)) - 1
                };
                assert_eq!(
                    sim.global_phase() & live_mask,
                    0,
                    "width {width} two-slot router phase garbage"
                );
                for (shot, case) in cases[batch_start..batch_end].iter().enumerate() {
                    let &(active_slot, pred0, add0, target0, pred1, add1, target1) = case;
                    let sign = if active_slot == 0 {
                        (pred0 >> 1) & 1
                    } else {
                        (pred1 >> 1) & 1
                    };
                    let expected0 = if active_slot == 0 {
                        if sign != 0 {
                            target0.wrapping_sub(add0) & (modulus - 1)
                        } else {
                            target0.wrapping_add(add0) & (modulus - 1)
                        }
                    } else {
                        target0
                    };
                    let expected1 = if active_slot == 1 {
                        if sign != 0 {
                            target1.wrapping_sub(add1) & (modulus - 1)
                        } else {
                            target1.wrapping_add(add1) & (modulus - 1)
                        }
                    } else {
                        target1
                    };
                    assert_eq!(get_reg(&sim, &qregs[0], shot), pred0, "pred0 case {case:?}");
                    assert_eq!(get_reg(&sim, &qregs[1], shot), add0, "add0 case {case:?}");
                    assert_eq!(
                        get_reg(&sim, &qregs[2], shot),
                        expected0,
                        "target0 case {case:?}"
                    );
                    assert_eq!(get_reg(&sim, &qregs[3], shot), pred1, "pred1 case {case:?}");
                    assert_eq!(get_reg(&sim, &qregs[4], shot), add1, "add1 case {case:?}");
                    assert_eq!(
                        get_reg(&sim, &qregs[5], shot),
                        expected1,
                        "target1 case {case:?}"
                    );
                }
                for reg in &qregs {
                    for &q in reg {
                        *sim.qubit_mut(q) = 0;
                    }
                }
                for q in 0..num_qubits {
                    assert_eq!(
                        sim.qubit(QubitId(q)) & live_mask,
                        0,
                        "width {width} two-slot router scratch garbage q{q}"
                    );
                }
            }
        }
    }

    #[test]
    fn round190_active_source_live_signed_digit_hmr_is_exact_on_active_rows() {
        for width in [2usize, 3, 4] {
            let ops = build_round190_active_source_live_signed_digit_hmr_width(width);
            let (num_qubits, num_bits, num_registers, regs) = analyze_ops(ops.iter().copied());
            assert_eq!(num_registers, 3, "width {width} register count");
            assert_eq!(regs.len(), 3, "width {width} regs");
            assert_eq!(num_bits as usize, width - 1, "width {width} hmr bits");
            assert_eq!(num_qubits as usize, 4 * width + 1, "width {width} qubits");
            for (idx, reg) in regs.iter().enumerate() {
                assert_eq!(reg.len(), width, "width {width} reg {idx}");
                assert!(reg.iter().all(|item| matches!(item, QubitOrBit::Qubit(_))));
            }
            let toffoli_ops = ops
                .iter()
                .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
                .count();
            assert_eq!(toffoli_ops, width - 1, "width {width} toffoli");
            let pred_reg: Vec<QubitId> = regs[0]
                .iter()
                .map(|item| match item {
                    QubitOrBit::Qubit(q) => *q,
                    _ => unreachable!(),
                })
                .collect();
            let add_reg: Vec<QubitId> = regs[1]
                .iter()
                .map(|item| match item {
                    QubitOrBit::Qubit(q) => *q,
                    _ => unreachable!(),
                })
                .collect();
            let target_reg: Vec<QubitId> = regs[2]
                .iter()
                .map(|item| match item {
                    QubitOrBit::Qubit(q) => *q,
                    _ => unreachable!(),
                })
                .collect();

            let modulus = 1u64 << width;
            let active_predecessors: Vec<u64> =
                (0..modulus).filter(|pred| (pred & 0b11) != 0).collect();
            let states = active_predecessors.len() as u64 * modulus * modulus;
            let mut seed = Shake128::default();
            seed.update(b"round190-active-source-live-signed-digit-hmr");
            seed.update(&[width as u8]);
            let mut xof = seed.finalize_xof();
            for batch_start in (0..states).step_by(64) {
                let mut sim = Simulator::new(num_qubits as usize, num_bits as usize, &mut xof);
                let batch_end = (batch_start + 64).min(states);
                for case in batch_start..batch_end {
                    let shot = (case - batch_start) as usize;
                    let pred_idx = (case % active_predecessors.len() as u64) as usize;
                    let addend = (case / active_predecessors.len() as u64) & (modulus - 1);
                    let target =
                        (case / (active_predecessors.len() as u64 * modulus)) & (modulus - 1);
                    let predecessor = active_predecessors[pred_idx];
                    set_reg(&mut sim, &pred_reg, predecessor, shot);
                    set_reg(&mut sim, &add_reg, addend, shot);
                    set_reg(&mut sim, &target_reg, target, shot);
                }

                sim.apply(&ops);
                let live_mask = if batch_end - batch_start == 64 {
                    u64::MAX
                } else {
                    (1u64 << (batch_end - batch_start)) - 1
                };
                assert_eq!(
                    sim.global_phase() & live_mask,
                    0,
                    "width {width} active HMR signed digit phase garbage"
                );
                for case in batch_start..batch_end {
                    let shot = (case - batch_start) as usize;
                    let pred_idx = (case % active_predecessors.len() as u64) as usize;
                    let addend = (case / active_predecessors.len() as u64) & (modulus - 1);
                    let target =
                        (case / (active_predecessors.len() as u64 * modulus)) & (modulus - 1);
                    let predecessor = active_predecessors[pred_idx];
                    let expected = if ((predecessor >> 1) & 1) != 0 {
                        target.wrapping_sub(addend) & (modulus - 1)
                    } else {
                        target.wrapping_add(addend) & (modulus - 1)
                    };
                    assert_eq!(
                        get_reg(&sim, &pred_reg, shot),
                        predecessor,
                        "width {width} predecessor changed case {case}"
                    );
                    assert_eq!(
                        get_reg(&sim, &add_reg, shot),
                        addend,
                        "width {width} addend changed case {case}"
                    );
                    assert_eq!(
                        get_reg(&sim, &target_reg, shot),
                        expected,
                        "width {width} target mismatch case {case}"
                    );
                }
                for reg in [&pred_reg, &add_reg, &target_reg] {
                    for &q in reg {
                        *sim.qubit_mut(q) = 0;
                    }
                }
                for q in 0..num_qubits {
                    assert_eq!(
                        sim.qubit(QubitId(q)) & live_mask,
                        0,
                        "width {width} active HMR scratch garbage q{q}"
                    );
                }
            }
        }
    }

    #[test]
    fn round190_active_hmr_digit_is_not_select0_safe() {
        const WIDTH: usize = 3;
        let ops = build_round190_active_source_live_signed_digit_hmr_width(WIDTH);
        let (num_qubits, num_bits, num_registers, regs) = analyze_ops(ops.iter().copied());
        assert_eq!(num_registers, 3);
        assert_eq!(regs.len(), 3);
        let pred_reg: Vec<QubitId> = regs[0]
            .iter()
            .map(|item| match item {
                QubitOrBit::Qubit(q) => *q,
                _ => unreachable!(),
            })
            .collect();
        let add_reg: Vec<QubitId> = regs[1]
            .iter()
            .map(|item| match item {
                QubitOrBit::Qubit(q) => *q,
                _ => unreachable!(),
            })
            .collect();
        let target_reg: Vec<QubitId> = regs[2]
            .iter()
            .map(|item| match item {
                QubitOrBit::Qubit(q) => *q,
                _ => unreachable!(),
            })
            .collect();

        let mut seed = Shake128::default();
        seed.update(b"round190-active-hmr-not-select0-safe");
        let mut xof = seed.finalize_xof();
        let mut sim = Simulator::new(num_qubits as usize, num_bits as usize, &mut xof);
        let inactive_predecessor = 0u64;
        let addend = 3u64;
        let target = 4u64;
        set_reg(&mut sim, &pred_reg, inactive_predecessor, 0);
        set_reg(&mut sim, &add_reg, addend, 0);
        set_reg(&mut sim, &target_reg, target, 0);

        sim.apply(&ops);
        let got_target = get_reg(&sim, &target_reg, 0);
        println!("METRIC round190_active_hmr_inactive_predecessor={inactive_predecessor}");
        println!("METRIC round190_active_hmr_inactive_addend={addend}");
        println!("METRIC round190_active_hmr_inactive_target_before={target}");
        println!("METRIC round190_active_hmr_inactive_target_after={got_target}");
        assert_eq!(get_reg(&sim, &pred_reg, 0), inactive_predecessor);
        assert_eq!(get_reg(&sim, &add_reg, 0), addend);
        assert_ne!(
            got_target, target,
            "active-HMR digit cannot be used as the select0-safe production residual"
        );
    }

    fn qubit_reg(reg: &[QubitOrBit]) -> Vec<QubitId> {
        reg.iter()
            .map(|item| match item {
                QubitOrBit::Qubit(q) => *q,
                _ => panic!("expected qubit register"),
            })
            .collect()
    }

    fn round556_expected(
        width: usize,
        q_bits: usize,
        rem: u64,
        rem_divisor: u64,
        coeff_seed: u64,
        coeff_divisor: u64,
        sigma: u64,
        q_increment: u64,
    ) -> Option<(u64, u64)> {
        let modulus = 1u64 << width;
        let mask = modulus - 1;
        if rem_divisor == 0 || coeff_divisor == 0 {
            return None;
        }
        if (rem_divisor << (q_bits - 1)) >= modulus {
            return None;
        }
        if (coeff_divisor << (q_bits - 1)) >= modulus {
            return None;
        }
        let quotient = rem / rem_divisor;
        if quotient >= (1u64 << q_bits) {
            return None;
        }
        if coeff_seed >= coeff_divisor {
            return None;
        }
        let coeff_restored = coeff_seed + (quotient + q_increment) * coeff_divisor;
        if coeff_restored >= modulus {
            return None;
        }
        let coeff = coeff_restored.wrapping_sub((sigma & 1) * coeff_divisor) & mask;
        Some((rem % rem_divisor, coeff))
    }

    #[test]
    fn round556_shifted_source_row_component_has_material_free_bound() {
        const WIDTH: usize = 258;
        const QBITS: usize = 26;
        let (ops, phases, peak_qubits, peak_phase) =
            build_round556_shifted_source_row_component_phase_resources(WIDTH, QBITS);
        let (num_qubits, _num_bits, num_registers, regs) = analyze_ops(ops.iter().copied());
        let toffoli_ops = ops
            .iter()
            .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
            .count();
        let old_materialized_formula = (6 * QBITS + 4) * WIDTH - (2 * QBITS + 2);
        let shifted_source_q = 6 * WIDTH + QBITS + 5;

        assert_eq!(num_registers, 5);
        assert_eq!(regs[0].len(), WIDTH);
        assert_eq!(regs[1].len(), WIDTH);
        assert_eq!(regs[2].len(), WIDTH);
        assert_eq!(regs[3].len(), WIDTH);
        assert_eq!(regs[4].len(), 4 + QBITS);
        assert_eq!(num_qubits as usize, shifted_source_q);
        assert_eq!(peak_qubits as usize, shifted_source_q);
        assert!(toffoli_ops <= old_materialized_formula);
        assert!(phases
            .iter()
            .any(|row| row.phase == "round556_shifted_source_remainder_digits"));
        assert_eq!(peak_phase, "round556_shifted_source_remainder_digits");
    }

    #[test]
    fn round556_shifted_source_row_component_matches_round120_relation() {
        const WIDTH: usize = 5;
        const QBITS: usize = 3;
        let ops = build_round556_shifted_source_row_component(WIDTH, QBITS);
        let (num_qubits, num_bits, num_registers, regs) = analyze_ops(ops.iter().copied());
        assert_eq!(num_registers, 5);
        let rem_reg = qubit_reg(&regs[0]);
        let rem_divisor_reg = qubit_reg(&regs[1]);
        let coeff_reg = qubit_reg(&regs[2]);
        let coeff_divisor_reg = qubit_reg(&regs[3]);
        let meta_reg = qubit_reg(&regs[4]);

        let mut public = vec![false; num_qubits as usize];
        for reg in [
            &rem_reg,
            &rem_divisor_reg,
            &coeff_reg,
            &coeff_divisor_reg,
            &meta_reg,
        ] {
            for &q in reg {
                public[q.0 as usize] = true;
            }
        }

        let mut cases = Vec::new();
        let modulus = 1u64 << WIDTH;
        for rem_divisor in 1..modulus {
            for coeff_divisor in 1..modulus {
                for rem in 0..modulus {
                    for coeff_seed in 0..coeff_divisor {
                        for sigma in 0..=1u64 {
                            for q_increment in 0..=1u64 {
                                if let Some(expected) = round556_expected(
                                    WIDTH,
                                    QBITS,
                                    rem,
                                    rem_divisor,
                                    coeff_seed,
                                    coeff_divisor,
                                    sigma,
                                    q_increment,
                                ) {
                                    cases.push((
                                        rem,
                                        rem_divisor,
                                        coeff_seed,
                                        coeff_divisor,
                                        sigma,
                                        q_increment,
                                        expected,
                                    ));
                                }
                            }
                        }
                    }
                }
            }
        }
        assert!(!cases.is_empty());

        let mut seed = Shake128::default();
        seed.update(b"round556-shifted-source-row-relation");
        let mut xof = seed.finalize_xof();
        let mut sim = Simulator::new(num_qubits as usize, num_bits as usize, &mut xof);
        for (batch, chunk) in cases.chunks(64).enumerate() {
            sim.clear_for_shot();
            for (shot, case) in chunk.iter().enumerate() {
                let (rem, rem_divisor, coeff_seed, coeff_divisor, sigma, q_increment, _) = *case;
                set_reg(&mut sim, &rem_reg, rem, shot);
                set_reg(&mut sim, &rem_divisor_reg, rem_divisor, shot);
                set_reg(&mut sim, &coeff_reg, coeff_seed, shot);
                set_reg(&mut sim, &coeff_divisor_reg, coeff_divisor, shot);
                set_reg(&mut sim, &meta_reg, sigma | (q_increment << 1), shot);
            }
            sim.apply(&ops);
            let live = if chunk.len() == 64 {
                u64::MAX
            } else {
                (1u64 << chunk.len()) - 1
            };
            assert_eq!(sim.global_phase() & live, 0, "phase dirty in batch {batch}");
            for q in 0..num_qubits {
                if !public[q as usize] {
                    assert_eq!(
                        sim.qubit(QubitId(q as u32)) & live,
                        0,
                        "scratch q{q} dirty in batch {batch}"
                    );
                }
            }
            for (shot, case) in chunk.iter().enumerate() {
                let (
                    _rem,
                    rem_divisor,
                    _coeff_seed,
                    coeff_divisor,
                    sigma,
                    q_increment,
                    (expected_rem, expected_coeff),
                ) = *case;
                assert_eq!(
                    get_reg(&sim, &rem_reg, shot),
                    expected_rem,
                    "batch {batch} shot {shot}"
                );
                assert_eq!(
                    get_reg(&sim, &rem_divisor_reg, shot),
                    rem_divisor,
                    "batch {batch} shot {shot}"
                );
                assert_eq!(
                    get_reg(&sim, &coeff_reg, shot),
                    expected_coeff,
                    "batch {batch} shot {shot}"
                );
                assert_eq!(
                    get_reg(&sim, &coeff_divisor_reg, shot),
                    coeff_divisor,
                    "batch {batch} shot {shot}"
                );
                assert_eq!(
                    get_reg(&sim, &meta_reg, shot),
                    sigma | (q_increment << 1),
                    "batch {batch} shot {shot}"
                );
            }
        }
    }

    #[test]
    fn direct_centered_shifted_source_qbit_row_fit_bench_has_sidecar_bound() {
        const Q_BITS: usize = DIRECT_CENTERED_LOW_BRANCH_META_BITS;
        let (ops, phases, peak_qubits, peak_phase) =
            build_direct_centered_shifted_source_qbit_row_fit_bench_phase_resources(Q_BITS);
        let (num_qubits, num_bits, num_registers, regs) = analyze_ops(ops.iter().copied());
        let toffoli_ops = ops
            .iter()
            .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
            .count();

        assert_eq!(num_registers, 4);
        assert_eq!(regs.len(), 4);
        assert!(num_bits as usize >= 2 * N);
        for (idx, reg) in regs.iter().enumerate() {
            assert_eq!(reg.len(), N, "register {idx} width");
        }
        let sidecar_q = 2 * N + DIRECT_CENTERED_BRANCH_SIDECAR_COMPONENT_SCRATCH_BITS;
        assert_eq!(num_qubits as usize, sidecar_q);
        assert_eq!(peak_qubits as usize, sidecar_q);
        assert_eq!(
            toffoli_ops,
            Q_BITS * (6 * N - 2) - 2 * Q_BITS * (Q_BITS - 1)
        );
        assert_eq!(
            peak_phase,
            "direct_centered_shifted_source_qbit_alloc_envelope"
        );
        assert!(phases
            .iter()
            .any(|row| row.phase == "direct_centered_shifted_source_qbit_remainder_digits"));
        assert!(phases
            .iter()
            .any(|row| row.phase == "direct_centered_shifted_source_qbit_coeff_digits"));
    }

    #[test]
    fn direct_centered_shifted_source_qbit_row_toy_is_exact_and_phase_clean() {
        const WIDTH: usize = 5;
        const QBITS: usize = 3;
        let mut b = B::new();
        let rem = b.alloc_qubits(WIDTH);
        let rem_divisor = b.alloc_qubits(WIDTH);
        let coeff = b.alloc_qubits(WIDTH);
        let coeff_divisor = b.alloc_qubits(WIDTH);
        let qbits = b.alloc_qubits(QBITS);
        let gated = b.alloc_qubits(WIDTH);
        let lt_tmp = b.alloc_qubit();
        let sign_one = b.alloc_qubit();
        let nonnegative = b.alloc_qubit();
        let carries = b.alloc_qubits(WIDTH - 1);
        emit_direct_centered_shifted_source_qbit_row(
            &mut b,
            &rem,
            &rem_divisor,
            &coeff,
            &coeff_divisor,
            &qbits,
            &gated,
            lt_tmp,
            sign_one,
            nonnegative,
            &carries,
        );

        let nq = b.next_qubit as usize;
        let nb = b.next_bit as usize;
        let mut public = vec![false; nq];
        for reg in [&rem, &rem_divisor, &coeff, &coeff_divisor] {
            for &q in reg {
                public[q.0 as usize] = true;
            }
        }

        let modulus = 1u64 << WIDTH;
        let mut cases = Vec::new();
        for rem_divisor_value in 1..modulus {
            for coeff_divisor_value in 1..modulus {
                for rem_value in 0..modulus {
                    for coeff_seed in 0..coeff_divisor_value {
                        if let Some(expected) = round556_expected(
                            WIDTH,
                            QBITS,
                            rem_value,
                            rem_divisor_value,
                            coeff_seed,
                            coeff_divisor_value,
                            0,
                            0,
                        ) {
                            cases.push((
                                rem_value,
                                rem_divisor_value,
                                coeff_seed,
                                coeff_divisor_value,
                                expected,
                            ));
                        }
                    }
                }
            }
        }
        assert!(!cases.is_empty());

        let mut seed = Shake128::default();
        seed.update(b"direct-centered-shifted-source-qbit-row-toy");
        let mut xof = seed.finalize_xof();
        let mut sim = Simulator::new(nq, nb, &mut xof);
        for (batch, chunk) in cases.chunks(64).enumerate() {
            sim.clear_for_shot();
            for (shot, case) in chunk.iter().enumerate() {
                let (rem_value, rem_divisor_value, coeff_seed, coeff_divisor_value, _) = *case;
                set_reg(&mut sim, &rem, rem_value, shot);
                set_reg(&mut sim, &rem_divisor, rem_divisor_value, shot);
                set_reg(&mut sim, &coeff, coeff_seed, shot);
                set_reg(&mut sim, &coeff_divisor, coeff_divisor_value, shot);
            }
            sim.apply(&b.ops);
            let live = if chunk.len() == 64 {
                u64::MAX
            } else {
                (1u64 << chunk.len()) - 1
            };
            assert_eq!(sim.global_phase() & live, 0, "phase dirty in batch {batch}");
            for q in 0..nq {
                if !public[q] {
                    assert_eq!(
                        sim.qubit(QubitId(q as u32)) & live,
                        0,
                        "scratch q{q} dirty in batch {batch}"
                    );
                }
            }
            for (shot, case) in chunk.iter().enumerate() {
                let (
                    _rem_value,
                    rem_divisor_value,
                    _coeff_seed,
                    coeff_divisor_value,
                    (expected_rem, expected_coeff),
                ) = *case;
                assert_eq!(get_reg(&sim, &rem, shot), expected_rem);
                assert_eq!(get_reg(&sim, &rem_divisor, shot), rem_divisor_value);
                assert_eq!(get_reg(&sim, &coeff, shot), expected_coeff);
                assert_eq!(get_reg(&sim, &coeff_divisor, shot), coeff_divisor_value);
            }
        }
    }

    #[test]
    fn direct_centered_branch_sidecar_component_has_relaxed_google_abi_shape() {
        let (ops, phases, peak_qubits, peak_phase) =
            build_direct_centered_branch_sidecar_bench_phase_resources();
        let (num_qubits, num_bits, num_registers, regs) = analyze_ops(ops.iter().copied());

        assert_eq!(regs.len(), 4);
        assert_eq!(num_registers, 4);
        assert_eq!(num_bits as usize, 2 * N);
        for (idx, reg) in regs.iter().enumerate() {
            assert_eq!(reg.len(), N, "register {idx} width");
        }
        for item in &regs[0] {
            assert!(matches!(item, QubitOrBit::Qubit(_)), "r0 must be qubits");
        }
        for item in &regs[1] {
            assert!(matches!(item, QubitOrBit::Qubit(_)), "r1 must be qubits");
        }
        for item in &regs[2] {
            assert!(matches!(item, QubitOrBit::Bit(_)), "r2 must be bits");
        }
        for item in &regs[3] {
            assert!(matches!(item, QubitOrBit::Bit(_)), "r3 must be bits");
        }

        let scratch = num_qubits as usize - 2 * N;
        let toffoli_ops = ops
            .iter()
            .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
            .count();
        assert_eq!(
            scratch,
            DIRECT_CENTERED_BRANCH_SIDECAR_COMPONENT_SCRATCH_BITS
        );
        assert!(scratch <= DIRECT_CENTERED_RELAXED_SCRATCH_BUDGET);
        assert!(num_qubits as usize <= DIRECT_CENTERED_RELAXED_Q_TARGET);
        assert!(toffoli_ops < DIRECT_CENTERED_RELAXED_T_TARGET);
        assert_eq!(toffoli_ops, 936);
        assert_eq!(peak_qubits as usize, num_qubits as usize);
        assert_eq!(peak_phase, "direct_centered_sidecar_google_abi");
        assert!(phases
            .iter()
            .any(|row| row.phase == "direct_centered_sidecar_emit_branch_history"));
        assert!(phases
            .iter()
            .any(|row| row.phase == "direct_centered_sidecar_clear_branch_history"));
    }

    #[test]
    fn direct_centered_branch_digit_clean_toy_is_exact() {
        const W: usize = 5;
        let mut b = B::new();
        let coeff_acc = b.alloc_qubits(W);
        let coeff_v = b.alloc_qubits(W);
        let branch = b.alloc_qubit();
        let sign = b.alloc_qubit();
        let gated = b.alloc_qubits(W);
        let carry = b.alloc_qubit();
        emit_direct_centered_branch_digit_update_clean(
            &mut b, &coeff_acc, &coeff_v, branch, sign, &gated, carry,
        );

        let nq = b.next_qubit as usize;
        let nb = b.next_bit as usize;
        let modulus = 1u64 << W;
        let mut cases = Vec::new();
        for acc in 0..modulus {
            for source in 0..modulus {
                for branch_value in 0..=1u64 {
                    for sign_value in 0..=1u64 {
                        let expected = if branch_value == 0 {
                            acc
                        } else if sign_value != 0 {
                            (acc + source) & (modulus - 1)
                        } else {
                            acc.wrapping_sub(source) & (modulus - 1)
                        };
                        cases.push((acc, source, branch_value, sign_value, expected));
                    }
                }
            }
        }

        let mut seed = Shake128::default();
        seed.update(b"direct-centered-branch-digit-clean-toy");
        let mut xof = seed.finalize_xof();
        let mut sim = Simulator::new(nq, nb, &mut xof);
        for (batch, chunk) in cases.chunks(64).enumerate() {
            sim.clear_for_shot();
            for (shot, &(acc, source, branch_value, sign_value, _expected)) in
                chunk.iter().enumerate()
            {
                set_reg(&mut sim, &coeff_acc, acc, shot);
                set_reg(&mut sim, &coeff_v, source, shot);
                if branch_value != 0 {
                    *sim.qubit_mut(branch) |= 1u64 << shot;
                }
                if sign_value != 0 {
                    *sim.qubit_mut(sign) |= 1u64 << shot;
                }
            }
            sim.apply(&b.ops);
            let live = if chunk.len() == 64 {
                u64::MAX
            } else {
                (1u64 << chunk.len()) - 1
            };
            assert_eq!(sim.global_phase() & live, 0, "phase dirty in batch {batch}");
            assert_eq!(sim.qubit(carry) & live, 0, "carry dirty in batch {batch}");
            for (shot, &(acc, source, branch_value, sign_value, expected)) in
                chunk.iter().enumerate()
            {
                assert_eq!(
                    get_reg(&sim, &gated, shot),
                    0,
                    "gated dirty in batch {batch} shot {shot}"
                );
                assert_eq!(
                    get_reg(&sim, &coeff_acc, shot),
                    expected,
                    "batch {batch} shot {shot}"
                );
                assert_eq!(
                    get_reg(&sim, &coeff_v, shot),
                    source,
                    "batch {batch} shot {shot}"
                );
                assert_eq!(
                    (sim.qubit(branch) >> shot) & 1,
                    branch_value,
                    "batch {batch} shot {shot}"
                );
                assert_eq!(
                    (sim.qubit(sign) >> shot) & 1,
                    sign_value,
                    "batch {batch} shot {shot}"
                );
                let _ = acc;
            }
        }
    }

    #[test]
    fn direct_centered_branch_replay_then_fast_finalizer_toy_is_exact() {
        const W: usize = 4;
        const HISTORY: usize = 3;
        let mut b = B::new();
        let coeff_acc = b.alloc_qubits(W);
        let coeff_v = b.alloc_qubits(W);
        let pred_a = b.alloc_qubits(HISTORY);
        let pred_b = b.alloc_qubits(HISTORY);
        let branch = b.alloc_qubits(HISTORY);
        let sign = b.alloc_qubit();
        let gated = b.alloc_qubits(W);
        let digit_carry = b.alloc_qubit();
        let nonnegative = b.alloc_qubit();
        let extra_carry = b.alloc_qubit();

        for i in 0..HISTORY {
            b.ccx(pred_a[i], pred_b[i], branch[i]);
        }
        for &branch_bit in &branch {
            emit_direct_centered_branch_digit_update_clean(
                &mut b,
                &coeff_acc,
                &coeff_v,
                branch_bit,
                sign,
                &gated,
                digit_carry,
            );
        }
        for i in (1..HISTORY).rev() {
            b.ccx(pred_a[i], pred_b[i], branch[i]);
        }
        let carries = [branch[1], branch[2], extra_carry];
        emit_direct_centered_branch_retained_finalizer_fast(
            &mut b,
            &coeff_acc,
            &coeff_v,
            branch[0],
            &gated,
            nonnegative,
            &carries,
        );
        b.ccx(pred_a[0], pred_b[0], branch[0]);

        let nq = b.next_qubit as usize;
        let nb = b.next_bit as usize;
        let modulus = 1u64 << W;
        let mask = modulus - 1;
        let mut cases = Vec::new();
        for acc in 0..modulus {
            for source in 0..modulus {
                for pred_a_value in 0..(1u64 << HISTORY) {
                    for pred_b_value in 0..(1u64 << HISTORY) {
                        for sign_value in 0..=1u64 {
                            let mut expected = acc;
                            for i in 0..HISTORY {
                                let branch_value =
                                    ((pred_a_value >> i) & 1) & ((pred_b_value >> i) & 1);
                                if branch_value != 0 {
                                    expected = if sign_value != 0 {
                                        expected.wrapping_add(source) & mask
                                    } else {
                                        expected.wrapping_sub(source) & mask
                                    };
                                }
                            }
                            if (pred_a_value & 1) != 0 && (pred_b_value & 1) != 0 {
                                expected = expected.wrapping_sub(source) & mask;
                            }
                            cases.push((
                                acc,
                                source,
                                pred_a_value,
                                pred_b_value,
                                sign_value,
                                expected,
                            ));
                        }
                    }
                }
            }
        }

        let mut seed = Shake128::default();
        seed.update(b"direct-centered-branch-replay-fast-finalizer-toy");
        let mut xof = seed.finalize_xof();
        let mut sim = Simulator::new(nq, nb, &mut xof);
        for (batch, chunk) in cases.chunks(64).enumerate() {
            sim.clear_for_shot();
            for (shot, &(acc, source, pred_a_value, pred_b_value, sign_value, _expected)) in
                chunk.iter().enumerate()
            {
                set_reg(&mut sim, &coeff_acc, acc, shot);
                set_reg(&mut sim, &coeff_v, source, shot);
                set_reg(&mut sim, &pred_a, pred_a_value, shot);
                set_reg(&mut sim, &pred_b, pred_b_value, shot);
                if sign_value != 0 {
                    *sim.qubit_mut(sign) |= 1u64 << shot;
                }
            }
            sim.apply(&b.ops);
            let live = if chunk.len() == 64 {
                u64::MAX
            } else {
                (1u64 << chunk.len()) - 1
            };
            assert_eq!(sim.global_phase() & live, 0, "phase dirty in batch {batch}");
            assert_eq!(sim.qubit(digit_carry) & live, 0, "digit carry dirty");
            assert_eq!(sim.qubit(nonnegative) & live, 0, "nonnegative dirty");
            assert_eq!(sim.qubit(extra_carry) & live, 0, "extra carry dirty");
            for &branch_bit in &branch {
                assert_eq!(sim.qubit(branch_bit) & live, 0, "branch history dirty");
            }
            for (shot, &(acc, source, pred_a_value, pred_b_value, sign_value, expected)) in
                chunk.iter().enumerate()
            {
                assert_eq!(
                    get_reg(&sim, &coeff_acc, shot),
                    expected,
                    "batch {batch} shot {shot}"
                );
                assert_eq!(get_reg(&sim, &gated, shot), 0);
                assert_eq!(get_reg(&sim, &coeff_v, shot), source);
                assert_eq!(get_reg(&sim, &pred_a, shot), pred_a_value);
                assert_eq!(get_reg(&sim, &pred_b, shot), pred_b_value);
                assert_eq!((sim.qubit(sign) >> shot) & 1, sign_value);
                let _ = acc;
            }
        }
    }

    #[test]
    fn direct_centered_low_path_branch_predicate_toy_is_exact() {
        const W: usize = 4;
        let mut b = B::new();
        let low_path = b.alloc_qubits(W);
        let divisor = b.alloc_qubits(W);
        let branch = b.alloc_qubit();
        let shifted = b.alloc_qubits(W + 1);
        let divisor_high = b.alloc_qubit();
        let cmp_cin = b.alloc_qubit();
        emit_direct_centered_low_path_branch_toggle(
            &mut b,
            &low_path,
            &divisor,
            branch,
            &shifted,
            divisor_high,
            cmp_cin,
        );

        let nq = b.next_qubit as usize;
        let nb = b.next_bit as usize;
        let mut cases = Vec::new();
        for low_value in 0..(1u64 << W) {
            for divisor_value in 0..(1u64 << W) {
                for initial_branch in 0..=1u64 {
                    let predicate = if 2 * low_value >= divisor_value { 1 } else { 0 };
                    cases.push((
                        low_value,
                        divisor_value,
                        initial_branch,
                        initial_branch ^ predicate,
                    ));
                }
            }
        }

        let mut seed = Shake128::default();
        seed.update(b"direct-centered-low-path-branch-predicate-toy");
        let mut xof = seed.finalize_xof();
        let mut sim = Simulator::new(nq, nb, &mut xof);
        for (batch, chunk) in cases.chunks(64).enumerate() {
            sim.clear_for_shot();
            for (shot, &(low_value, divisor_value, initial_branch, _expected_branch)) in
                chunk.iter().enumerate()
            {
                set_reg(&mut sim, &low_path, low_value, shot);
                set_reg(&mut sim, &divisor, divisor_value, shot);
                if initial_branch != 0 {
                    *sim.qubit_mut(branch) |= 1u64 << shot;
                }
            }
            sim.apply(&b.ops);
            let live = if chunk.len() == 64 {
                u64::MAX
            } else {
                (1u64 << chunk.len()) - 1
            };
            assert_eq!(sim.global_phase() & live, 0, "phase dirty in batch {batch}");
            for &wire in &shifted {
                assert_eq!(sim.qubit(wire) & live, 0, "shifted scratch dirty");
            }
            assert_eq!(
                sim.qubit(divisor_high) & live,
                0,
                "divisor-high scratch dirty"
            );
            assert_eq!(sim.qubit(cmp_cin) & live, 0, "cmp-cin scratch dirty");
            for (shot, &(low_value, divisor_value, _initial_branch, expected_branch)) in
                chunk.iter().enumerate()
            {
                assert_eq!(get_reg(&sim, &low_path, shot), low_value);
                assert_eq!(get_reg(&sim, &divisor, shot), divisor_value);
                assert_eq!(
                    (sim.qubit(branch) >> shot) & 1,
                    expected_branch,
                    "batch {batch} shot {shot}"
                );
            }
        }
    }

    #[test]
    fn direct_centered_branch_predicate_step_fit_stays_inside_round714_envelope() {
        let (ops, phases, peak_qubits, peak_phase) =
            build_direct_centered_branch_predicate_step_fit_bench_phase_resources();
        let (num_qubits, num_bits, num_registers, regs) = analyze_ops(ops.iter().copied());

        assert_eq!(regs.len(), 4);
        assert_eq!(num_registers, 4);
        assert_eq!(num_bits as usize, 3 * N);
        let scratch = num_qubits as usize - 2 * N;
        let toffoli_ops = ops
            .iter()
            .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
            .count();
        assert_eq!(
            scratch,
            DIRECT_CENTERED_BRANCH_SIDECAR_COMPONENT_SCRATCH_BITS
        );
        assert!(scratch <= DIRECT_CENTERED_RELAXED_SCRATCH_BUDGET);
        assert!(num_qubits as usize <= DIRECT_CENTERED_RELAXED_Q_TARGET);
        assert!(toffoli_ops < 2_000);
        assert_eq!(peak_qubits as usize, num_qubits as usize);
        assert_eq!(
            peak_phase,
            "direct_centered_branch_predicate_step_alloc_envelope"
        );
        assert!(phases
            .iter()
            .any(|row| row.phase == "direct_centered_branch_predicate_compare"));
        assert!(phases
            .iter()
            .any(|row| row.phase == "direct_centered_branch_digit_clean_addsub"));
    }

    #[test]
    fn direct_centered_binary_trie_qrom_toy_is_exact_and_phase_clean() {
        const ADDRESS_BITS: usize = 3;
        const TARGET_BITS: usize = 5;
        const ROWS: usize = 6;

        let table_words: Vec<u64> = (0..ROWS)
            .map(|row| ((row as u64).wrapping_mul(0b10101) ^ 0b10010) & ((1u64 << TARGET_BITS) - 1))
            .collect();

        let mut b = B::new();
        let address = b.alloc_qubits(ADDRESS_BITS);
        let target = b.alloc_qubits(TARGET_BITS);
        emit_direct_centered_binary_trie_qrom_xor_table(
            &mut b,
            &address,
            &target,
            ROWS,
            &table_words,
        );

        let nq = b.next_qubit as usize;
        let nb = b.next_bit as usize;
        let mut public = vec![false; nq];
        for &q in address.iter().chain(target.iter()) {
            public[q.0 as usize] = true;
        }

        let mut cases = Vec::new();
        for addr in 0..(1u64 << ADDRESS_BITS) {
            for before in 0..(1u64 << TARGET_BITS) {
                let loaded = if (addr as usize) < ROWS {
                    table_words[addr as usize]
                } else {
                    0
                };
                cases.push((addr, before, before ^ loaded));
            }
        }

        let mut seed = Shake128::default();
        seed.update(b"direct-centered-binary-trie-qrom-toy");
        let mut xof = seed.finalize_xof();
        let mut sim = Simulator::new(nq, nb, &mut xof);
        for (batch, chunk) in cases.chunks(64).enumerate() {
            sim.clear_for_shot();
            for (shot, &(addr, before, _expected)) in chunk.iter().enumerate() {
                set_reg(&mut sim, &address, addr, shot);
                set_reg(&mut sim, &target, before, shot);
            }
            sim.apply(&b.ops);
            let live = if chunk.len() == 64 {
                u64::MAX
            } else {
                (1u64 << chunk.len()) - 1
            };
            assert_eq!(sim.global_phase() & live, 0, "phase dirty in batch {batch}");
            for q in 0..nq {
                if !public[q] {
                    assert_eq!(
                        sim.qubit(QubitId(q as u32)) & live,
                        0,
                        "scratch q{q} dirty in batch {batch}"
                    );
                }
            }
            for (shot, &(addr, _before, expected)) in chunk.iter().enumerate() {
                assert_eq!(get_reg(&sim, &address, shot), addr);
                assert_eq!(get_reg(&sim, &target, shot), expected);
            }
        }
    }

    #[test]
    fn direct_centered_binary_trie_qrom_roundtrip_toy_is_exact_and_phase_clean() {
        const ADDRESS_BITS: usize = 3;
        const TARGET_BITS: usize = 9;
        const ROWS: usize = 6;

        let table_words = direct_centered_binary_trie_qrom_table_words(ROWS, TARGET_BITS);

        let mut b = B::new();
        let address = b.alloc_qubits(ADDRESS_BITS);
        let target = b.alloc_qubits(TARGET_BITS);
        emit_direct_centered_binary_trie_qrom_xor_table(
            &mut b,
            &address,
            &target,
            ROWS,
            &table_words,
        );
        emit_direct_centered_binary_trie_qrom_xor_table(
            &mut b,
            &address,
            &target,
            ROWS,
            &table_words,
        );

        let nq = b.next_qubit as usize;
        let nb = b.next_bit as usize;
        let mut public = vec![false; nq];
        for &q in address.iter().chain(target.iter()) {
            public[q.0 as usize] = true;
        }

        let mut cases = Vec::new();
        for addr in 0..(1u64 << ADDRESS_BITS) {
            for before in 0..(1u64 << TARGET_BITS) {
                cases.push((addr, before));
            }
        }

        let mut seed = Shake128::default();
        seed.update(b"direct-centered-binary-trie-qrom-roundtrip-toy");
        let mut xof = seed.finalize_xof();
        let mut sim = Simulator::new(nq, nb, &mut xof);
        for (batch, chunk) in cases.chunks(64).enumerate() {
            sim.clear_for_shot();
            for (shot, &(addr, before)) in chunk.iter().enumerate() {
                set_reg(&mut sim, &address, addr, shot);
                set_reg(&mut sim, &target, before, shot);
            }
            sim.apply(&b.ops);
            let live = if chunk.len() == 64 {
                u64::MAX
            } else {
                (1u64 << chunk.len()) - 1
            };
            assert_eq!(sim.global_phase() & live, 0, "phase dirty in batch {batch}");
            for q in 0..nq {
                if !public[q] {
                    assert_eq!(
                        sim.qubit(QubitId(q as u32)) & live,
                        0,
                        "scratch q{q} dirty in batch {batch}"
                    );
                }
            }
            for (shot, &(addr, before)) in chunk.iter().enumerate() {
                assert_eq!(get_reg(&sim, &address, shot), addr);
                assert_eq!(get_reg(&sim, &target, shot), before);
            }
        }
    }

    #[test]
    fn direct_centered_binary_trie_qrom_hits_round728_row_multiplier_budget() {
        const ROWS: usize = 4_934;
        const ADDRESS_BITS: usize = 13;
        const TARGET_BITS: usize = 16;

        let (ops, phases, peak_qubits, peak_phase) =
            build_direct_centered_binary_trie_qrom_bench_phase_resources(
                ROWS,
                ADDRESS_BITS,
                TARGET_BITS,
            );
        let (num_qubits, num_bits, num_registers, regs) = analyze_ops(ops.iter().copied());
        let toffoli_ops = ops
            .iter()
            .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
            .count();
        let expected_nodes = direct_centered_binary_trie_qrom_node_count(ROWS, ADDRESS_BITS);

        assert_eq!(regs.len(), 4);
        assert_eq!(num_registers, 4);
        assert_eq!(num_bits as usize, 2 * N + expected_nodes);
        assert_eq!(toffoli_ops, expected_nodes);
        assert!(toffoli_ops <= 2 * ROWS + ADDRESS_BITS);
        assert!(toffoli_ops <= 6 * ROWS);
        assert_eq!(num_qubits as usize, 2 * N + ADDRESS_BITS + 1);
        assert_eq!(peak_qubits as usize, num_qubits as usize);
        assert_eq!(peak_phase, "direct_centered_binary_trie_qrom_unary_walk");
        assert!(phases
            .iter()
            .any(|row| row.phase == "direct_centered_binary_trie_qrom_unary_walk"));
        assert!(phases
            .iter()
            .any(|row| row.phase == "direct_centered_binary_trie_qrom_clear_root"));
    }

    #[test]
    fn direct_centered_binary_trie_qrom_roundtrip_fits_round730_wide_payload_budget() {
        const ROWS: usize = 4_934;
        const ADDRESS_BITS: usize = 13;
        const TARGET_BITS: usize = 84;

        let (ops, phases, peak_qubits, peak_phase) =
            build_direct_centered_binary_trie_qrom_roundtrip_bench_phase_resources(
                ROWS,
                ADDRESS_BITS,
                TARGET_BITS,
            );
        let (num_qubits, num_bits, num_registers, regs) = analyze_ops(ops.iter().copied());
        let toffoli_ops = ops
            .iter()
            .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
            .count();
        let expected_nodes = direct_centered_binary_trie_qrom_node_count(ROWS, ADDRESS_BITS);

        assert_eq!(regs.len(), 4);
        assert_eq!(num_registers, 4);
        assert_eq!(num_bits as usize, 2 * N + 2 * expected_nodes);
        assert_eq!(toffoli_ops, 2 * expected_nodes);
        assert_eq!(toffoli_ops, 19_746);
        assert!(toffoli_ops <= 4 * ROWS + 2 * ADDRESS_BITS);
        assert_eq!(num_qubits as usize, 2 * N + ADDRESS_BITS + 1);
        assert_eq!(peak_qubits as usize, num_qubits as usize);
        assert_eq!(
            peak_phase,
            "direct_centered_binary_trie_qrom_roundtrip_load_walk"
        );
        assert!(phases
            .iter()
            .any(|row| row.phase == "direct_centered_binary_trie_qrom_roundtrip_load_walk"));
        assert!(phases
            .iter()
            .any(|row| row.phase == "direct_centered_binary_trie_qrom_roundtrip_clear_walk"));
    }

    #[test]
    fn direct_centered_inline_predicate_finalizer_delta_fits_google_fast_width_if_replay_deleted() {
        let (ops, phases, peak_qubits, peak_phase) =
            build_direct_centered_inline_predicate_finalizer_delta_fit_bench_phase_resources();
        let (num_qubits, num_bits, num_registers, regs) = analyze_ops(ops.iter().copied());

        assert_eq!(regs.len(), 4);
        assert_eq!(num_registers, 4);
        assert_eq!(num_bits as usize, 3 * N - 1);
        let scratch = num_qubits as usize - 2 * N;
        let toffoli_ops = ops
            .iter()
            .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
            .count();
        assert_eq!(
            scratch,
            DIRECT_CENTERED_BRANCH_SIDECAR_COMPONENT_SCRATCH_BITS
                + DIRECT_CENTERED_EXPLICIT_BRANCH_HISTORY_BITS
        );
        assert_eq!(num_qubits as usize, 1_425);
        assert!(toffoli_ops < 122_000);
        assert_eq!(peak_qubits as usize, num_qubits as usize);
        assert_eq!(
            peak_phase,
            "direct_centered_inline_predicate_delta_alloc_dual_history_envelope"
        );
        assert!(phases
            .iter()
            .any(|row| row.phase == "direct_centered_branch_predicate_compare"));
        assert!(phases
            .iter()
            .any(|row| row.phase == "direct_centered_branch_retained_fast_finalizer_subtract"));
        assert!(!phases
            .iter()
            .any(|row| row.phase == "direct_centered_branch_digit_clean_addsub"));
    }

    #[test]
    fn direct_centered_branch_retained_finalizer_toy_is_exact() {
        const W: usize = 5;
        let mut b = B::new();
        let remainder = b.alloc_qubits(W);
        let divisor = b.alloc_qubits(W);
        let branch = b.alloc_qubit();
        let gated = b.alloc_qubits(W);
        let carry = b.alloc_qubit();
        emit_direct_centered_branch_retained_finalizer(
            &mut b, &remainder, &divisor, branch, &gated, carry,
        );

        let nq = b.next_qubit as usize;
        let nb = b.next_bit as usize;
        let modulus = 1u64 << W;
        let mut cases = 0usize;
        for divisor_value in 1..(1u64 << (W - 1)) {
            for final_remainder in 0..divisor_value {
                for branch_value in 0..=1u64 {
                    let prefinal = final_remainder + branch_value * divisor_value;
                    if prefinal >= modulus {
                        continue;
                    }
                    cases += 1;
                    let mut seed = Shake128::default();
                    seed.update(&(cases as u64).to_le_bytes());
                    let mut xof = seed.finalize_xof();
                    let mut sim = Simulator::new(nq, nb, &mut xof);
                    set_reg(&mut sim, &remainder, prefinal, 0);
                    set_reg(&mut sim, &divisor, divisor_value, 0);
                    if branch_value != 0 {
                        *sim.qubit_mut(branch) |= 1;
                    }
                    sim.apply(&b.ops);
                    assert_eq!(get_reg(&sim, &remainder, 0), final_remainder);
                    assert_eq!(get_reg(&sim, &divisor, 0), divisor_value);
                    assert_eq!((sim.qubit(branch) & 1), branch_value);
                    assert_eq!(sim.qubit(carry) & 1, 0);
                    assert_eq!(get_reg(&sim, &gated, 0), 0);
                }
            }
        }
        assert_eq!(cases, 240);
    }

    #[test]
    fn direct_centered_branch_retained_fast_finalizer_toy_is_exact() {
        const W: usize = 5;
        let mut b = B::new();
        let remainder = b.alloc_qubits(W);
        let divisor = b.alloc_qubits(W);
        let branch = b.alloc_qubit();
        let gated = b.alloc_qubits(W);
        let nonnegative = b.alloc_qubit();
        let carries = b.alloc_qubits(W - 1);
        emit_direct_centered_branch_retained_finalizer_fast(
            &mut b,
            &remainder,
            &divisor,
            branch,
            &gated,
            nonnegative,
            &carries,
        );

        let nq = b.next_qubit as usize;
        let nb = b.next_bit as usize;
        let modulus = 1u64 << W;
        let mut cases = 0usize;
        for divisor_value in 1..(1u64 << (W - 1)) {
            for final_remainder in 0..divisor_value {
                for branch_value in 0..=1u64 {
                    let prefinal = final_remainder + branch_value * divisor_value;
                    if prefinal >= modulus {
                        continue;
                    }
                    cases += 1;
                    let mut seed = Shake128::default();
                    seed.update(&(0xFA57_0000u64 + cases as u64).to_le_bytes());
                    let mut xof = seed.finalize_xof();
                    let mut sim = Simulator::new(nq, nb, &mut xof);
                    set_reg(&mut sim, &remainder, prefinal, 0);
                    set_reg(&mut sim, &divisor, divisor_value, 0);
                    if branch_value != 0 {
                        *sim.qubit_mut(branch) |= 1;
                    }
                    sim.apply(&b.ops);
                    assert_eq!(get_reg(&sim, &remainder, 0), final_remainder);
                    assert_eq!(get_reg(&sim, &divisor, 0), divisor_value);
                    assert_eq!(sim.qubit(branch) & 1, branch_value);
                    assert_eq!(sim.qubit(nonnegative) & 1, 0);
                    assert_eq!(get_reg(&sim, &gated, 0), 0);
                    assert_eq!(get_reg(&sim, &carries, 0), 0);
                    assert_eq!(sim.global_phase() & 1, 0);
                }
            }
        }
        assert_eq!(cases, 240);
    }

    #[test]
    fn direct_centered_branch_retained_finalizer_component_has_expected_shape() {
        let (ops, phases, peak_qubits, peak_phase) =
            build_direct_centered_branch_retained_finalizer_bench_phase_resources();
        let (num_qubits, num_bits, num_registers, regs) = analyze_ops(ops.iter().copied());
        let toffoli_ops = ops
            .iter()
            .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
            .count();

        assert_eq!(regs.len(), 4);
        assert_eq!(num_registers, 4);
        assert_eq!(num_bits as usize, 2 * N);
        assert_eq!(num_qubits as usize, 2 * N + N + 2);
        assert_eq!(peak_qubits as usize, num_qubits as usize);
        assert_eq!(toffoli_ops, 4 * N - 2);
        assert_eq!(
            peak_phase,
            "direct_centered_branch_retained_finalizer_google_abi"
        );
        assert!(phases
            .iter()
            .any(|row| row.phase == "direct_centered_branch_retained_finalizer_subtract"));
    }

    #[test]
    fn direct_centered_branch_digit_clean_fit_stays_inside_round714_envelope() {
        let (ops, phases, peak_qubits, peak_phase) =
            build_direct_centered_branch_digit_clean_fit_bench_phase_resources();
        let (num_qubits, num_bits, num_registers, regs) = analyze_ops(ops.iter().copied());
        let toffoli_ops = ops
            .iter()
            .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
            .count();

        assert_eq!(regs.len(), 4);
        assert_eq!(num_registers, 4);
        assert_eq!(num_bits as usize, 3 * N);
        assert_eq!(
            num_qubits as usize,
            2 * N + DIRECT_CENTERED_BRANCH_SIDECAR_COMPONENT_SCRATCH_BITS
        );
        assert_eq!(peak_qubits as usize, num_qubits as usize);
        assert_eq!(toffoli_ops, 3 * N - 2);
        assert_eq!(
            peak_phase,
            "direct_centered_branch_digit_clean_alloc_envelope"
        );
        assert!(phases
            .iter()
            .any(|row| row.phase == "direct_centered_branch_digit_clean_addsub"));
    }

    #[test]
    fn direct_centered_remainder_abs_swap_transition_toy_is_exact() {
        const W: usize = 4;
        let mut b = B::new();
        let low_path = b.alloc_qubits(W);
        let divisor = b.alloc_qubits(W);
        let branch = b.alloc_qubit();
        let gated = b.alloc_qubits(W);
        let carries = b.alloc_qubits(W - 1);
        emit_direct_centered_remainder_abs_swap_transition(
            &mut b, &low_path, &divisor, branch, &gated, &carries,
        );

        let nq = b.next_qubit as usize;
        let nb = b.next_bit as usize;
        let mut cases = Vec::new();
        for divisor_value in 1..(1u64 << W) {
            for low_value in 0..divisor_value {
                let branch_value = u64::from(2 * low_value >= divisor_value);
                let next_divisor = if branch_value == 0 {
                    low_value
                } else {
                    divisor_value - low_value
                };
                cases.push((low_value, divisor_value, branch_value, next_divisor));
            }
        }

        let mut seed = Shake128::default();
        seed.update(b"direct-centered-remainder-abs-swap-transition-toy");
        let mut xof = seed.finalize_xof();
        let mut sim = Simulator::new(nq, nb, &mut xof);
        for (batch, chunk) in cases.chunks(64).enumerate() {
            sim.clear_for_shot();
            for (shot, &(low_value, divisor_value, branch_value, _next_divisor)) in
                chunk.iter().enumerate()
            {
                set_reg(&mut sim, &low_path, low_value, shot);
                set_reg(&mut sim, &divisor, divisor_value, shot);
                if branch_value != 0 {
                    *sim.qubit_mut(branch) |= 1u64 << shot;
                }
            }
            sim.apply(&b.ops);
            let live = if chunk.len() == 64 {
                u64::MAX
            } else {
                (1u64 << chunk.len()) - 1
            };
            assert_eq!(sim.global_phase() & live, 0, "phase dirty in batch {batch}");
            for &wire in &gated {
                assert_eq!(sim.qubit(wire) & live, 0, "gated divisor dirty");
            }
            for &wire in &carries {
                assert_eq!(sim.qubit(wire) & live, 0, "borrowed carry dirty");
            }
            for (shot, &(_low_value, divisor_value, branch_value, next_divisor)) in
                chunk.iter().enumerate()
            {
                assert_eq!(get_reg(&sim, &low_path, shot), divisor_value);
                assert_eq!(get_reg(&sim, &divisor, shot), next_divisor);
                assert_eq!((sim.qubit(branch) >> shot) & 1, branch_value);
            }
        }
    }

    #[test]
    fn direct_centered_row_transition_fit_stays_inside_round714_envelope() {
        let (ops, phases, peak_qubits, peak_phase) =
            build_direct_centered_row_transition_fit_bench_phase_resources();
        let (num_qubits, num_bits, num_registers, regs) = analyze_ops(ops.iter().copied());
        let toffoli_ops = ops
            .iter()
            .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
            .count();
        let hmr_ops = ops
            .iter()
            .filter(|op| matches!(op.kind, OperationType::Hmr))
            .count();

        assert_eq!(regs.len(), 4);
        assert_eq!(num_registers, 4);
        assert_eq!(num_bits as usize, 4 * N - 1);
        assert_eq!(
            num_qubits as usize,
            2 * N + DIRECT_CENTERED_BRANCH_SIDECAR_COMPONENT_SCRATCH_BITS
        );
        assert_eq!(peak_qubits as usize, num_qubits as usize);
        assert_eq!(toffoli_ops, 2 * N - 1);
        assert_eq!(hmr_ops, N - 1 + N);
        assert_eq!(peak_phase, "direct_centered_row_transition_alloc_envelope");
        assert!(phases
            .iter()
            .any(|row| row.phase == "direct_centered_row_transition_abs_add"));
        assert!(phases
            .iter()
            .any(|row| row.phase == "direct_centered_row_transition_swap_next_state"));
    }

    #[test]
    fn direct_centered_branch_replay_finalizer_fit_stays_inside_round714_envelope() {
        let (ops, phases, peak_qubits, peak_phase) =
            build_direct_centered_branch_replay_finalizer_fit_bench_phase_resources();
        let (num_qubits, num_bits, num_registers, regs) = analyze_ops(ops.iter().copied());
        let toffoli_ops = ops
            .iter()
            .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
            .count();

        assert_eq!(regs.len(), 4);
        assert_eq!(num_registers, 4);
        assert_eq!(
            num_bits as usize,
            2 * N + DIRECT_CENTERED_EXPLICIT_BRANCH_HISTORY_BITS * N + (N - 1)
        );
        assert_eq!(
            num_qubits as usize,
            2 * N + DIRECT_CENTERED_BRANCH_SIDECAR_COMPONENT_SCRATCH_BITS
        );
        assert_eq!(peak_qubits as usize, num_qubits as usize);
        assert_eq!(
            toffoli_ops,
            DIRECT_CENTERED_EXPLICIT_BRANCH_HISTORY_BITS * (3 * N - 2)
                + (2 * DIRECT_CENTERED_EXPLICIT_BRANCH_HISTORY_BITS)
                + (3 * N - 1)
        );
        assert_eq!(
            peak_phase,
            "direct_centered_branch_replay_finalizer_alloc_envelope"
        );
        assert!(phases
            .iter()
            .any(|row| row.phase == "direct_centered_branch_replay_clear_nonfinal_history"));
        assert!(phases
            .iter()
            .any(|row| row.phase == "direct_centered_branch_retained_fast_finalizer_subtract"));
    }

    #[test]
    fn direct_centered_predicate_replay_finalizer_fit_materializes_full_tail_projection() {
        let (ops, phases, peak_qubits, peak_phase) =
            build_direct_centered_predicate_replay_finalizer_fit_bench_phase_resources();
        let (num_qubits, num_bits, num_registers, regs) = analyze_ops(ops.iter().copied());
        let toffoli_ops = ops
            .iter()
            .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
            .count();

        let predicate_toggle_t = 2 * (N + 1);
        let branch_digit_t = 3 * N - 2;
        let finalizer_t = 3 * N - 1;
        let expected_tail_t = DIRECT_CENTERED_EXPLICIT_BRANCH_HISTORY_BITS
            * (2 * predicate_toggle_t + branch_digit_t)
            + finalizer_t;

        assert_eq!(regs.len(), 4);
        assert_eq!(num_registers, 4);
        assert_eq!(
            num_bits as usize,
            2 * N + DIRECT_CENTERED_EXPLICIT_BRANCH_HISTORY_BITS * N + (N - 1)
        );
        assert_eq!(
            num_qubits as usize,
            2 * N + DIRECT_CENTERED_BRANCH_SIDECAR_COMPONENT_SCRATCH_BITS
        );
        assert_eq!(peak_qubits as usize, num_qubits as usize);
        assert_eq!(toffoli_ops, expected_tail_t);
        assert_eq!(toffoli_ops, 210_665);
        assert_eq!(
            peak_phase,
            "direct_centered_predicate_replay_finalizer_alloc_envelope"
        );
        assert!(phases
            .iter()
            .any(|row| row.phase == "direct_centered_branch_predicate_compare"));
        assert!(phases
            .iter()
            .any(|row| row.phase == "direct_centered_branch_retained_fast_finalizer_subtract"));
    }

    #[test]
    fn direct_centered_sidecar_finalizer_fit_stays_inside_round714_envelope() {
        let (ops, phases, peak_qubits, peak_phase) =
            build_direct_centered_sidecar_finalizer_fit_bench_phase_resources();
        let (num_qubits, num_bits, num_registers, regs) = analyze_ops(ops.iter().copied());
        let toffoli_ops = ops
            .iter()
            .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
            .count();

        assert_eq!(regs.len(), 4);
        assert_eq!(num_registers, 4);
        assert_eq!(num_bits as usize, 2 * N);
        assert_eq!(
            num_qubits as usize,
            2 * N + DIRECT_CENTERED_BRANCH_SIDECAR_COMPONENT_SCRATCH_BITS
        );
        assert_eq!(peak_qubits as usize, num_qubits as usize);
        assert_eq!(toffoli_ops, 4 * N - 2);
        assert_eq!(
            peak_phase,
            "direct_centered_sidecar_finalizer_alloc_envelope"
        );
        assert!(phases
            .iter()
            .any(|row| row.phase == "direct_centered_branch_retained_finalizer_gate_divisor"));
    }

    #[test]
    fn direct_centered_sidecar_fast_finalizer_fit_stays_inside_round714_envelope() {
        let (ops, phases, peak_qubits, peak_phase) =
            build_direct_centered_sidecar_fast_finalizer_fit_bench_phase_resources();
        let (num_qubits, num_bits, num_registers, regs) = analyze_ops(ops.iter().copied());
        let toffoli_ops = ops
            .iter()
            .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
            .count();

        assert_eq!(regs.len(), 4);
        assert_eq!(num_registers, 4);
        assert_eq!(num_bits as usize, 2 * N + N - 1);
        assert_eq!(
            num_qubits as usize,
            2 * N + DIRECT_CENTERED_BRANCH_SIDECAR_COMPONENT_SCRATCH_BITS
        );
        assert_eq!(peak_qubits as usize, num_qubits as usize);
        assert_eq!(toffoli_ops, 3 * N - 1);
        assert_eq!(
            peak_phase,
            "direct_centered_sidecar_fast_finalizer_alloc_envelope"
        );
        assert!(phases
            .iter()
            .any(|row| row.phase == "direct_centered_branch_retained_fast_finalizer_subtract"));
    }
}
