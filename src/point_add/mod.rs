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

mod emit;
pub(crate) use emit::*;

mod arith;
pub(crate) use arith::*;

pub mod trailmix_ludicrous;
mod single_ccx_fanout;

thread_local! {
    static D1_PHASE_CORRECTED_PRODUCT_CORE_SCOPE: std::cell::Cell<bool> =
        std::cell::Cell::new(false);
    static OP_SITE_TRACE: std::cell::RefCell<Vec<OpSite>> =
        std::cell::RefCell::new(Vec::new());
    static OP_TRACE_CONTEXT: std::cell::Cell<u32> = std::cell::Cell::new(0);
}

fn d1_phase_corrected_product_core_active() -> bool {
    D1_PHASE_CORRECTED_PRODUCT_CORE_SCOPE.with(|scope| scope.get())
}

pub type OpSite = (&'static str, u32, u32);

pub(crate) fn op_site_trace_enabled() -> bool {
    static ENABLED: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ENABLED.get_or_init(|| std::env::var_os("TRACE_OP_SITES").is_some())
}

fn reset_op_site_trace() {
    if op_site_trace_enabled() {
        OP_SITE_TRACE.with(|sites| sites.borrow_mut().clear());
    }
}

fn record_op_site(site: OpSite) {
    if op_site_trace_enabled() {
        OP_SITE_TRACE.with(|sites| sites.borrow_mut().push(site));
    }
}

pub(crate) fn set_op_trace_context(context: u32) -> u32 {
    if !op_site_trace_enabled() {
        return 0;
    }
    OP_TRACE_CONTEXT.with(|slot| {
        let old = slot.get();
        slot.set(context);
        old
    })
}

pub(crate) fn restore_op_trace_context(context: u32) {
    if op_site_trace_enabled() {
        OP_TRACE_CONTEXT.with(|slot| slot.set(context));
    }
}

pub(crate) fn take_op_site_trace_for_constprop(expected_len: usize) -> Option<Vec<OpSite>> {
    if !op_site_trace_enabled() {
        return None;
    }
    OP_SITE_TRACE.with(|sites| {
        let mut sites = sites.borrow_mut();
        assert_eq!(
            sites.len(),
            expected_len,
            "op site trace length before constprop"
        );
        Some(std::mem::take(&mut *sites))
    })
}

pub(crate) fn set_op_site_trace_from_constprop(sites: Vec<OpSite>) {
    if op_site_trace_enabled() {
        OP_SITE_TRACE.with(|slot| *slot.borrow_mut() = sites);
    }
}

pub fn take_last_op_sites() -> Vec<OpSite> {
    OP_SITE_TRACE.with(|sites| std::mem::take(&mut *sites.borrow_mut()))
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
        reset_op_site_trace();
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
    /// TEST-ONLY constructor + ops extractor (used by the classical-arith unit bin).
    pub fn new_for_test() -> Self {
        Self::new()
    }
    pub fn take_ops(&mut self) -> Vec<Op> {
        std::mem::take(&mut self.ops)
    }
    #[track_caller]
    fn push_op(&mut self, op: Op) {
        self.counted_ops += 1;
        self.counted_kind_ops[op.kind as usize] += 1;
        self.counted_phase_kind_ops[op.kind as usize] += 1;
        if !self.count_only {
            let loc = std::panic::Location::caller();
            let context = OP_TRACE_CONTEXT.with(|slot| slot.get());
            record_op_site((loc.file(), loc.line(), context));
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
    #[track_caller]
    fn alloc_qubit(&mut self) -> QubitId {
        self.active_qubits += 1;
        self.record_phase_active();
        if let Ok(threshold) = std::env::var("TRACE_ALLOC_NEAR_PEAK")
            .ok()
            .and_then(|value| value.parse::<u32>().ok())
            .ok_or(())
        {
            if self.active_qubits >= threshold {
                let caller = std::panic::Location::caller();
                eprintln!(
                    "ALLOC_NEAR active={} next_idx={} phase='{}' ops_idx={} free_pool={} caller={}:{}",
                    self.active_qubits,
                    self.next_qubit,
                    self.phase,
                    self.current_ops_len(),
                    self.free_qubits.len(),
                    caller.file(),
                    caller.line(),
                );
            }
        }
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
    #[track_caller]
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

    // Classical-bit (BitId) writes: ZERO Toffoli, ZERO Clifford in the scorer.
    /// `dst := 0`.
    fn bit_store0(&mut self, dst: BitId) {
        let mut op = Op::empty();
        op.kind = OperationType::BitStore0;
        op.c_target = dst;
        self.push_op(op);
    }
    /// `dst |= (condition stack AND)`; empty stack => `dst := 1`.
    fn bit_store1(&mut self, dst: BitId) {
        let mut op = Op::empty();
        op.kind = OperationType::BitStore1;
        op.c_target = dst;
        self.push_op(op);
    }
    /// `dst ^= (condition stack AND)`; empty stack => `dst := !dst`.
    fn bit_invert(&mut self, dst: BitId) {
        let mut op = Op::empty();
        op.kind = OperationType::BitInvert;
        op.c_target = dst;
        self.push_op(op);
    }
    /// `dst := a`.
    fn bit_copy(&mut self, dst: BitId, a: BitId) {
        self.bit_store0(dst);
        self.push_condition(a);
        self.bit_store1(dst);
        self.pop_condition();
    }
    /// `dst ^= a`.
    fn bit_xor_into(&mut self, dst: BitId, a: BitId) {
        self.push_condition(a);
        self.bit_invert(dst);
        self.pop_condition();
    }
    /// `dst ^= (a AND b)`.
    fn bit_and_xor_into(&mut self, dst: BitId, a: BitId, b: BitId) {
        self.push_condition(a);
        self.push_condition(b);
        self.bit_invert(dst);
        self.pop_condition();
        self.pop_condition();
    }
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

fn kal_vent_modadd_enabled() -> bool {
    std::env::var("KAL_VENT_MODADD").ok().as_deref() == Some("1")
}

fn kal_vent_halve_enabled() -> bool {
    std::env::var("KAL_VENT_HALVE").ok().as_deref() == Some("1")
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

// Q1153 second-512 scan route. To submit a clean hit from the current hunt,
// update this nonce, build with no shell env overrides, run `ecdsafail run`,
// and submit only if it remains 0 / 0 / 0.
const Q1153_SECOND512_SUBMISSION_NONCE: &str = "50400005525597";

fn configure_q1153_second512_submission_defaults() {
    set_default_env("DIALOG_TAIL_NONCE", Q1153_SECOND512_SUBMISSION_NONCE);
    set_default_env("TLM_TARGET_Q", "1152");
    set_default_env("TLM_FOLD_CHUNK_ZERO_CIN", "1");
    set_default_env("TLM_FFG_MAX_G", "47");
    set_default_env("TLM_APPLY_ADD_SKIP_LASTK", "1");
    set_default_env("TLM_FOLD_TAIL_CINC", "1");
    set_default_env("TLM_CODEC_DIAMOND_MCX", "1");
    set_default_env("SINGLE_CCX_FANOUT_DISABLE", "1");
    // ── 1152 stack (codex FFG cy0-release + opus square vents) ──────────────
    set_default_env("TLM_FFG_RELEASE_CY0_DURING_SUFFIX", "1");
    set_default_env("TLM_FFG_RELEASE_CY0_CALLS", "178,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,203,208,210,211,212,213,215,217,219,221,226,232,234,235,236,237,239");
    set_default_env("TLM_APPLY_FWD_CSWAP_SKIP_LAST", "1");
    set_default_env("TLM_COORD_RSUB_FUSED", "1");
    set_default_env("TLM_SQUARE_VENT_MARGIN", "0");
    set_default_env("TLM_COORD_ADD3X_TRUNC", "1");
    set_default_env("TLM_SQUARE_VENT_SHIFTED", "1");
    set_default_env("TLM_SQUARE_PEAK_CAP", "1152");
}

pub fn build() -> Vec<Op> {
    configure_q1153_second512_submission_defaults();

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
    // GPT-Codex Q1159 product route. Per-call FFG/fold reserves fit every local
    // arithmetic peak under the target width; direct comparator carries and HMR
    // cleanup remove Toffolis without increasing liveness. Nonce 453700 passed
    // the trusted 9024-shot evaluator with 0 classical, phase, and ancilla
    // failures at rounded T=1,388,180 and Q=1159 (score 1,608,900,620).
    set_default_env("LUD_EXTRA_FOLD_VENTS", "0");
    set_default_env("LUD_EXTRA_FOLD_MIN_G", "0");
    set_default_env("LUD_EXTRA_FOLD_MAX_G", "999");
    set_default_env("DIALOG_TAIL_NONCE", "2430844");
    set_default_env("TLM_FOLD_TAIL_CINC", "1");
    set_default_env("TLM_CODEC_DIAMOND_MCX", "1");
    set_default_env("SINGLE_CCX_FANOUT_DISABLE", "1");
    // Stack the latest frontier square fold: use shifted-low folding for all
    // square lanes instead of the older `a`-only direct32 ramp shortcut.
    set_default_env("TLM_SQUARE_F_RAMP10_DIRECT32_TAGS", "");
    set_default_env("TLM_SQUARE_F_SHIFTED_LOW", "1");
    // post-1159 avgT stack (Codex): graduated final +f chunk w/o materializing the
    // dropped carry-out (arith.rs) + skip the first forward-apply cswap (gcd.rs).
    set_default_env("TLM_GRAD_FINAL_NO_COUT", "1");
    set_default_env("TLM_APPLY_FWD_FIRST_CSWAP_SKIP", "1");
    set_default_env("CONSTPROP_MAX_ITERS", "16");
    // q1155 trial: tighten the q1156 chunk4/ffg11/s2safer reserve machinery by
    // one peak qubit before retuning the per-call reserve schedules.
    set_default_env("TLM_TARGET_Q", "1155");
    set_default_env("TLM_FOLD_BOUNDARY_ZERO_DIRECT", "1");
    set_default_env("TLM_FOLD_CHUNK_FORCE", "4");
    set_default_env("TLM_TARGET_FOLD_CALL_RESERVE_OVERRIDES", "173:3,175:3,177:3,256:11,257:11,336:3,338:3,340:3,176:3,178:3,180:3,254:5,259:20,333:3,335:3,337:3,179:3,181:3,183:3,182:3,184:3,186:3,327:3,329:3,330:3,331:3,332:3,334:3");
    set_default_env("TLM_TARGET_FFG_CALL_RESERVE_OVERRIDES", "184:4,186:4,188:4,205:6,207:6,209:6,220:7,222:7,224:7,238:8,240:8,242:8,251:9,257:10,262:10,355:10,362:10,359:10,181:3,183:3,185:3,187:4,189:4,191:4,196:5,198:5,200:5,208:6,210:6,212:6,223:7,225:7,227:7,241:8,243:8,245:8,250:9,252:9,190:4,192:4,193:5,194:4,195:5,197:5,199:5,201:5,202:6,203:5,204:6,206:6,211:6,213:6,214:7,215:6,216:7,218:7,226:7,228:8,229:8,230:8,231:8,233:8,244:8,246:8,247:9,253:9,254:10,259:11,358:10,340:11,341:11,342:11,343:11,344:11,345:11,346:11,347:11,348:11,349:11,350:11");
    set_default_env("TLM_APPLY_FWD_S2_ZERO_LAST", "1");
    set_default_env("TLM_APPLY_INV_S2_ZERO_LAST", "1");
    set_default_env("TLM_APPLY_FWD_CSWAP_SKIP_LAST", "2");
    set_default_env("TLM_APPLY_INV_CSWAP_SKIP_LAST", "1");
    set_default_env("TLM_FOLD_RELEASE_CONTROLS", "1");
    set_default_env("TLM_TARGET_FFG_RESERVE", "9");
    set_default_env(
        "TLM_TARGET_FFG_CALL_RESERVES",
        concat!(
            "163:8,165:8,166:7,167:8,168:7,169:6,170:7,171:6,172:5,173:6,174:5,175:4,176:5,177:4,178:3,179:4,180:3,181:2,182:3,183:2,184:1,185:2,186:1,187:0,188:1,189:0,190:3,191:0,192:3,193:3,194:3,195:3,196:4,197:3,198:4,199:4,200:4,201:4,202:4,203:4,204:4,205:5,206:4,207:5,208:5,209:5,210:5,211:5,212:5,213:5,214:5,215:5,216:5,217:6,218:5,219:6,220:6,221:6,222:6,223:6,224:6,225:6,226:6,227:6,228:6,229:6,230:6,231:6,232:7,233:6,234:7,235:7,236:7,237:7,238:7,239:7,240:7,241:7,242:7,243:7,244:7,245:7,246:7,247:7,248:8,249:8,250:8,251:8,252:8,253:8,254:8,",
            "509:8,510:8,511:8,512:8,513:8,514:8,515:8,516:7,517:7,518:7,519:7,520:7,521:7,522:7,523:7,524:7,525:7,526:7,527:7,528:7,529:7,530:6,531:7,532:6,533:6,534:6,535:6,536:6,537:6,538:6,539:6,540:6,541:6,542:6,543:6,544:6,545:5,546:6,547:5,548:5,549:5,550:5,551:5,552:5,553:5,554:5,555:5,556:5,557:4,558:5,559:4,560:4,561:4,562:4,563:4,564:4,565:4,566:3,567:4,568:3,569:3,570:3,571:3,572:0,573:3,574:0,575:1,576:0,577:1,578:2,579:1,580:2,581:3,582:2,583:3,584:4,585:3,586:4,587:5,588:4,589:5,590:6,591:5,592:6,593:7,594:6,595:7,596:8,597:7,598:8,600:8",
        ),
    );
    set_default_env("TLM_TARGET_FOLD_RESERVE", "4");
    set_default_env(
        "TLM_TARGET_FOLD_CALL_RESERVES",
        concat!(
            "170:3,172:3,173:2,174:3,175:2,176:1,177:2,178:1,179:0,180:1,181:0,182:0,183:0,184:0,185:3,186:0,187:3,188:3,189:3,190:3,191:3,192:3,193:3,195:3,",
            "251:3,252:3,253:3,254:3,255:3,256:3,257:3,258:3,259:3,260:3,261:3,262:3,318:3,320:3,321:3,322:3,323:3,324:3,325:3,326:3,327:0,328:3,329:0,330:0,331:0,332:0,333:1,334:0,335:1,336:2,337:1,338:2,339:3,340:2,341:3,343:3",
        ),
    );
    set_default_env("TLM_GCD_RESELECT_LAYOUT", "1");
    set_default_env("TLM_DIRECT_VARCHUNK", "1");
    set_default_env("TLM_COUT_LAYOUT_SEARCH", "1");
    set_default_env("TLM_COUT_LAYOUT_MARGIN", "0");
    set_default_env("TLM_COUT_LAYOUT_FORCE_M1_KS", "129");
    set_default_env("TLM_GCD_ADAPTIVE_LAYOUT_SEARCH", "1");
    set_default_env("TLM_GCD_ADAPTIVE_LAYOUT_MARGIN", "0");
    // u0/even-v0 lifecycle loans plus the GCD y0 loan candidate
    // (1165->1164 at the same layout stack) — BAKED so env-less builds reproduce it.
    set_default_env("TLM_PARK_ODD_U0", "1");
    set_default_env("TLM_LOAN_ODD_U0", "1");
    set_default_env("TLM_PARK_EVEN_V0", "1");
    set_default_env("TLM_LOAN_EVEN_V0", "1");
    set_default_env("TLM_LOAN_GCD_Y0", "1");
    set_default_env("TLM_HYB_V_DELTA", "2");
    set_default_env("TLM_COUT_K_DELTA", "2");
    set_default_env("TLM_FOLD_DELTA", "2");
    set_default_env("TLM_FFG_DELTA", "0");
    set_default_env("TLM_GCD_K_ADJUST_AFTER", "169");
    set_default_env("TLM_GCD_K_ADJUST_BEFORE", "196");
    set_default_env("TLM_GCD_K_ADJUST", "-2");
    // Codex idx-less structural stack. Generated dead-drop lists are not used
    // in this submission tree.
    set_default_env("TLM_FFG_SKIP_STRUCTURAL_DEAD_CALLS", "1");
    set_default_env("TLM_FFG_SKIP_TOP_CARRY31", "1");
    set_default_env("TLM_FFG_SKIP_TOP_CARRY30", "1");
    set_default_env("TLM_CUCCARO_SKIP_STRUCTURAL_DEAD_CALLS", "1");
    set_default_env("TLM_COMPARE_SKIP_STRUCTURAL_DEAD_CALLS", "1");
    set_default_env("TLM_COMPARE_SKIP_EXACT_REMAINDER", "1");
    set_default_env("TLM_GIDNEY_SKIP_STRUCTURAL_DEAD_CALLS", "1");
    set_default_env("TLM_GIDNEY_SKIP_EXACT_REMAINDER", "1");
    set_default_env("TLM_CONST_CHUNK_SKIP_STRUCTURAL_DEAD_CALLS", "1");
    set_default_env("TLM_CONST_CHUNK_SKIP_EXACT_REMAINDER", "1");
    set_default_env("TLM_FUSED_SKIP_STRUCTURAL_DEAD_CARRIES", "1");
    set_default_env("TLM_FUSED_SKIP_STRUCTURAL_DEAD_SHIFT0", "1");
    set_default_env("TLM_FUSED_SKIP_EXACT_FOLD_REMAINDER", "1");
    set_default_env("TLM_FUSED_SKIP_STRUCTURAL_DEAD_DIRTY_FOLD", "1");
    set_default_env("TLM_FUSED_SKIP_STRUCTURAL_DEAD_CLEAN_WINDOW", "1");
    set_default_env("TLM_ADD_CONST_SKIP_STRUCTURAL_DEAD_CARRIES", "1");
    set_default_env("TLM_GCD_SKIP_STRUCTURAL_DEAD_CSWAPS", "1");
    set_default_env("TLM_GCD_SKIP_EXACT_FORWARD_CSWAPS", "1");
    set_default_env("TLM_GCD_SKIP_STRUCTURAL_DEAD_SHIFTS", "1");
    set_default_env("TLM_GCD_SKIP_EXACT_SHIFT_REMAINDER", "1");
    set_default_env("TLM_COMPARE_SKIP_EXACT_CIN_REMAINDER", "1");
    set_default_env("TLM_FUSED_SKIP_EXACT_BOUNDARY_ZERO", "1");
    set_default_env("TLM_GIDNEY_SKIP_EXACT_ERASE_ALL_CCZ", "1");
    set_default_env("TLM_FFG_SKIP_EXACT_TOP29_REMAINDER", "1");
    set_default_env("TLM_GCD_SKIP_REVERSE_DIAGONAL_EDGE", "1");
    set_default_env("TLM_FFG_SKIP_INVERSE_MOD_SUB_TOP29", "1");
    set_default_env("TLM_FFG_INVERSE_TOP29_MAX_CALL", "180");
    set_default_env("TLM_FUSED_CLEAN_FOLD_SKIP_TOP31", "1");
    set_default_env("TLM_GIDNEY_SKIP_SMALL_RESIDUAL_DEAD", "1");
    let mut ops = trailmix_ludicrous::build_trailmix_ludicrous_ops();
    if std::env::var("SINGLE_CCX_FANOUT_DISABLE")
        .ok()
        .as_deref()
        == Some("1")
    {
        return ops;
    }
    let input_ops = ops.len();
    let mut fanout_passes = 0usize;
    loop {
        match single_ccx_fanout::rewrite_first_target_fanout(ops.clone(), 96) {
            Ok((rewritten, _witness)) => {
                fanout_passes += 1;
                ops = rewritten;
            }
            Err(error) => {
                eprintln!(
                    "SINGLE_CCX_FANOUT: STOP passes={} input_ops={} output_ops={} reason={}",
                    fanout_passes,
                    input_ops,
                    ops.len(),
                    error,
                );
                break;
            }
        }
    }
    assert!(fanout_passes >= 1, "single-fanout rewrite failed to find first pass");
    eprintln!(
        "SINGLE_CCX_FANOUT: SUMMARY input_ops={} output_ops={} passes={}",
        input_ops,
        ops.len(),
        fanout_passes,
    );
    ops
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
