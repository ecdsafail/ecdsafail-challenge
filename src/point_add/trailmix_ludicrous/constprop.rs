//! Sound classical constant-propagation / peephole pass over the emitted
//! op-stream.
//!
//! ## What it does
//! Forward abstract interpretation tracking, for every qubit and classical bit,
//! a provable constant value in the lattice {Zero, One, Unknown}. The seeds are
//! the *genuinely* classical initial states of the circuit:
//!   * every qubit starts at |0> in the simulator's `clear_for_shot`, EXCEPT the
//!     IO input registers reg0 (P.x) and reg1 (P.y) which `set_register` fills
//!     with per-shot input data -> those qubits are seeded Unknown;
//!   * every classical bit starts at 0.
//!
//! Propagation is sound by construction: a value is marked `One`/`Zero` only
//! when it provably holds for EVERY basis input on EVERY shot. The moment any
//! doubt exists (an unknown control, a conditional write that is not a no-op, a
//! measurement, a phase op, an input qubit) the value collapses to Unknown.
//!
//! ## The peephole, applied only when PROVABLY known
//!   * CCX with a known-0 quantum control  -> DROP (no-op, still scored).
//!   * CCX with a known-1 quantum control   -> FOLD to CX on the other control.
//!   * CCX with both controls known-1       -> FOLD to X on the target.
//! Dropping a counted CCX removes one executed-Toffoli; folding to CX/X moves it
//! to the uncounted (Clifford / classically-trackable) tier. Both preserve the
//! computed unitary on all basis states.
//!
//! ## Soundness corroboration (step 2)
//! `CONSTPROP_VERIFY=<n>` runs the *unmodified* op-stream through the real
//! `Simulator` over `n` diverse nonces x 9024 shots and records, for every gate
//! the static analysis flagged, whether its quantum control(s) were ALWAYS the
//! claimed constant across all those inputs. Any flagged gate that fails this
//! empirical check is reported (and would indicate a static-analysis bug); the
//! transform is then restricted to the intersection.

use crate::circuit::{BitId, NO_BIT, NO_QUBIT, Op, OperationType, QubitId};

/// Sentinel "never touched" timestamp for the inverse-pair last-touch arrays.
/// Treated as strictly before every op index (i.e. touched before all ops).
const NEVER: usize = usize::MAX;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Val {
    Zero,
    One,
    Unknown,
}

use Val::*;

/// Result of the static analysis: per index in the op-stream, what transform (if
/// any) applies, plus the bookkeeping for reporting.
#[derive(Clone, Copy, Debug, Default)]
pub struct ConstPropStats {
    pub ccx_total: usize,
    pub dropped: usize,
    pub folded_cx: usize,
    pub folded_x: usize,
}

/// Per-gate transform decision, recorded so the empirical verifier can check the
/// claimed control constancy on the *original* stream.
#[derive(Clone, Copy, Debug)]
enum Decision {
    /// Keep unchanged.
    Keep,
    /// Drop: the named control qubit is claimed provably Zero.
    DropZeroCtrl { ctrl: QubitId },
    /// Fold CCX -> CX on `keep_ctrl`; `one_ctrl` is claimed provably One.
    FoldCx { one_ctrl: QubitId, keep_ctrl: QubitId },
    /// Fold CCX -> X on target; both controls claimed provably One.
    FoldX { c1: QubitId, c2: QubitId },
    /// Drop: the two controls are provably COMPLEMENTARY (a = NOT b on every
    /// shot), so a AND b = 0 and the CCX is a no-op. `a`,`b` recorded for the
    /// empirical verifier.
    DropComplementCtrls { a: QubitId, b: QubitId },
    /// Fold CCX -> CX(keep_ctrl, target): the two controls are provably EQUAL
    /// (a = b on every shot), so a AND b = a. `a`,`b` recorded for verification.
    FoldEqualCtrls { a: QubitId, b: QubitId, keep_ctrl: QubitId },
}

struct Analyzer {
    q: Vec<Val>,
    b: Vec<Val>,
    /// condition stack of bits (pushed by PushCondition).
    cond_stack: Vec<BitId>,
}

impl Analyzer {
    fn qv(&self, id: QubitId) -> Val {
        if id == NO_QUBIT { Unknown } else { self.q[id.0 as usize] }
    }
    fn bv(&self, id: BitId) -> Val {
        if id == NO_BIT { Unknown } else { self.b[id.0 as usize] }
    }
    fn set_q(&mut self, id: QubitId, v: Val) {
        self.q[id.0 as usize] = v;
    }
    fn set_b(&mut self, id: BitId, v: Val) {
        self.b[id.0 as usize] = v;
    }

    /// Is the *effective* classical condition for this op provably always-true on
    /// every shot? That means: every bit on the condition stack is provably One,
    /// and the op's own `c_condition` (if any) is provably One.
    fn cond_always_true(&self, op: &Op) -> bool {
        for &c in &self.cond_stack {
            if self.bv(c) != One {
                return false;
            }
        }
        if op.c_condition != NO_BIT && self.bv(op.c_condition) != One {
            return false;
        }
        true
    }

    /// Could the op possibly NOT execute on some shot? (condition not provably
    /// always-true). If so, a write must be merged with the prior value.
    fn cond_maybe_false(&self, op: &Op) -> bool {
        !self.cond_always_true(op)
    }
}

/// XOR of two abstract values when both known.
fn xor_val(a: Val, b: Val) -> Val {
    match (a, b) {
        (Zero, x) | (x, Zero) => x,
        (One, One) => Zero,
        _ => Unknown,
    }
}

/// AND of two abstract values.
fn and_val(a: Val, b: Val) -> Val {
    match (a, b) {
        (Zero, _) | (_, Zero) => Zero,
        (One, One) => One,
        _ => Unknown,
    }
}

/// Merge: the target is written with `new` only on *some* shots (condition may be
/// false), keeping `old` on the rest. Provably constant only if old==new.
fn merge(old: Val, new: Val) -> Val {
    if old == new { old } else { Unknown }
}

/// Run the static forward analysis. `input_qubits` lists the qubit ids that are
/// seeded Unknown (the IO data registers); all other qubit ids are seeded Zero.
fn analyze(ops: &[Op], num_q: usize, num_b: usize, input_qubits: &[QubitId]) -> (Vec<Decision>, ConstPropStats) {
    let mut a = Analyzer {
        q: vec![Zero; num_q],
        b: vec![Zero; num_b],
        cond_stack: Vec::new(),
    };
    for &q in input_qubits {
        a.q[q.0 as usize] = Unknown;
    }

    let mut decisions = vec![Decision::Keep; ops.len()];
    let mut stats = ConstPropStats::default();

    for (i, op) in ops.iter().enumerate() {
        match op.kind {
            OperationType::PushCondition => {
                a.cond_stack.push(op.c_condition);
            }
            OperationType::PopCondition => {
                a.cond_stack.pop();
            }
            OperationType::CCX => {
                stats.ccx_total += 1;
                let c1 = a.qv(op.q_control1);
                let c2 = a.qv(op.q_control2);
                // Decide the transform from the controls (independent of the
                // classical condition: if a control is provably Zero the gate is
                // a no-op on every shot regardless of the condition mask).
                if c1 == Zero {
                    decisions[i] = Decision::DropZeroCtrl { ctrl: op.q_control1 };
                    stats.dropped += 1;
                    // No state change (no-op).
                } else if c2 == Zero {
                    decisions[i] = Decision::DropZeroCtrl { ctrl: op.q_control2 };
                    stats.dropped += 1;
                    // No state change.
                } else if c1 == One && c2 == One {
                    decisions[i] = Decision::FoldX { c1: op.q_control1, c2: op.q_control2 };
                    stats.folded_x += 1;
                    // Effect: target ^= cond. New target value:
                    let tgt = a.qv(op.q_target);
                    let nv = xor_val(tgt, One);
                    let res = if a.cond_maybe_false(op) { merge(tgt, nv) } else { nv };
                    a.set_q(op.q_target, res);
                } else if c1 == One {
                    decisions[i] = Decision::FoldCx { one_ctrl: op.q_control1, keep_ctrl: op.q_control2 };
                    stats.folded_cx += 1;
                    // Effect: target ^= cond & c2(unknown/one) -> unknown-ish.
                    let tgt = a.qv(op.q_target);
                    let delta = c2; // c2 is One or Unknown here
                    let nv = xor_val(tgt, delta);
                    let res = if a.cond_maybe_false(op) { merge(tgt, nv) } else { nv };
                    a.set_q(op.q_target, res);
                } else if c2 == One {
                    decisions[i] = Decision::FoldCx { one_ctrl: op.q_control2, keep_ctrl: op.q_control1 };
                    stats.folded_cx += 1;
                    let tgt = a.qv(op.q_target);
                    let delta = c1;
                    let nv = xor_val(tgt, delta);
                    let res = if a.cond_maybe_false(op) { merge(tgt, nv) } else { nv };
                    a.set_q(op.q_target, res);
                } else {
                    // Both controls unknown: target ^= cond & c1 & c2 -> if the
                    // and is provably Zero target unchanged, else Unknown.
                    let delta = and_val(c1, c2);
                    let tgt = a.qv(op.q_target);
                    let nv = xor_val(tgt, delta);
                    let res = if a.cond_maybe_false(op) { merge(tgt, nv) } else { nv };
                    a.set_q(op.q_target, res);
                }
            }
            OperationType::CX => {
                let ctrl = a.qv(op.q_control1);
                let tgt = a.qv(op.q_target);
                let nv = xor_val(tgt, ctrl);
                let res = if a.cond_maybe_false(op) { merge(tgt, nv) } else { nv };
                a.set_q(op.q_target, res);
            }
            OperationType::X => {
                let tgt = a.qv(op.q_target);
                let nv = xor_val(tgt, One);
                let res = if a.cond_maybe_false(op) { merge(tgt, nv) } else { nv };
                a.set_q(op.q_target, res);
            }
            OperationType::Swap => {
                // q_control1 <-> q_target, executed where cond is set.
                let va = a.qv(op.q_control1);
                let vt = a.qv(op.q_target);
                if a.cond_maybe_false(op) {
                    // partial swap: each side becomes merge(self, other).
                    a.set_q(op.q_control1, merge(va, vt));
                    a.set_q(op.q_target, merge(vt, va));
                } else {
                    a.set_q(op.q_control1, vt);
                    a.set_q(op.q_target, va);
                }
            }
            OperationType::R => {
                // Reset to |0> where cond set; else unchanged.
                let tgt = a.qv(op.q_target);
                let res = if a.cond_maybe_false(op) { merge(tgt, Zero) } else { Zero };
                a.set_q(op.q_target, res);
            }
            OperationType::Hmr => {
                // Measurement: classical bit becomes a fresh random value
                // (Unknown), qubit demolished to |0> where cond set.
                let res = if a.cond_maybe_false(op) {
                    merge(a.bv(op.c_target), Unknown)
                } else {
                    Unknown
                };
                a.set_b(op.c_target, res);
                let tgt = a.qv(op.q_target);
                let qres = if a.cond_maybe_false(op) { merge(tgt, Zero) } else { Zero };
                a.set_q(op.q_target, qres);
            }
            OperationType::BitStore0 => {
                let cur = a.bv(op.c_target);
                let res = if a.cond_maybe_false(op) { merge(cur, Zero) } else { Zero };
                a.set_b(op.c_target, res);
            }
            OperationType::BitStore1 => {
                let cur = a.bv(op.c_target);
                let res = if a.cond_maybe_false(op) { merge(cur, One) } else { One };
                a.set_b(op.c_target, res);
            }
            OperationType::BitInvert => {
                let cur = a.bv(op.c_target);
                let nv = xor_val(cur, One);
                let res = if a.cond_maybe_false(op) { merge(cur, nv) } else { nv };
                a.set_b(op.c_target, res);
            }
            // Phase-only ops, register metadata, debug: no effect on the
            // computational basis values we track.
            OperationType::Z
            | OperationType::CZ
            | OperationType::CCZ
            | OperationType::Neg
            | OperationType::Register
            | OperationType::AppendToRegister
            | OperationType::DebugPrint => {}
        }
    }

    (decisions, stats)
}

// ── Affine GF(2) relation analysis ───────────────────────────────────────────
//
// A complementary refinement to the const-prop lattice. It tracks, for every
// qubit, an EXACT affine form over GF(2):
//     value(q) = const(q)  XOR  (XOR of the symbolic source variables in set(q))
// Sources: each input data qubit (reg0/reg1) is its own independent variable;
// every |0>-seeded ancilla starts as the empty set with const 0. A measurement
// result (Hmr) and a CCX target (nonlinear) become a FRESH independent variable
// (a well-defined but unconstrained value) so the relation algebra stays exact.
//
// With this we can PROVE two CCX controls a,b are related on EVERY shot:
//   * set(a)==set(b) AND const(a)==const(b)  ->  a == b always  -> a&b = a
//        => CCX(a,b,t) == CX(a,t)  (fold; removes a counted Toffoli).
//   * set(a)==set(b) AND const(a)!=const(b)  ->  a == NOT b      -> a&b = 0
//        => CCX(a,b,t) is a no-op  (drop; removes a counted Toffoli).
// Both hold regardless of the CCX's own classical condition (the relationship is
// proved over all shots), exactly like the constant-control drops.
//
// Soundness:
//   * By construction — a value's affine form is updated only by the exact GF(2)
//     transfer of each gate, and any write that might NOT happen (a possibly-
//     false condition) collapses the target to a fresh variable (precision loss
//     only, never an unsound claim). A set exceeding CAP_SET also collapses to a
//     fresh variable.
//   * Empirically — the CONSTPROP_VERIFY harness checks, for every flagged CCX,
//     that the controls were ALWAYS equal (or always complementary) across many
//     nonces x 9024 shots, reverting any that ever fail.

/// Cap on affine-set size; beyond it a qubit collapses to a fresh variable.
const CAP_SET: usize = 2048;

struct Affine {
    cst: Vec<bool>,
    set: Vec<Vec<u32>>,
    nextvar: u32,
    cond_stack: Vec<BitId>,
    // bit values reuse the const-prop {Zero,One,Unknown} lattice so we know when
    // a condition is provably always-true.
    b: Vec<Val>,
}

/// Symmetric difference (XOR) of two sorted, de-duplicated u32 sets.
fn xor_set(a: &[u32], b: &[u32]) -> Vec<u32> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0usize, 0usize);
    while i < a.len() && j < b.len() {
        if a[i] < b[j] {
            out.push(a[i]);
            i += 1;
        } else if a[i] > b[j] {
            out.push(b[j]);
            j += 1;
        } else {
            i += 1;
            j += 1;
        }
    }
    out.extend_from_slice(&a[i..]);
    out.extend_from_slice(&b[j..]);
    out
}

impl Affine {
    fn fresh(&mut self) -> Vec<u32> {
        let v = self.nextvar;
        self.nextvar += 1;
        vec![v]
    }
    fn bv(&self, id: BitId) -> Val {
        if id == NO_BIT { Unknown } else { self.b[id.0 as usize] }
    }
    /// Could this op fail to execute on some shot?
    fn cond_maybe_false(&self, op: &Op) -> bool {
        for &c in &self.cond_stack {
            if self.bv(c) != One {
                return true;
            }
        }
        if op.c_condition != NO_BIT && self.bv(op.c_condition) != One {
            return true;
        }
        false
    }
}

/// Forward affine-relation analysis. Returns per-op Decisions only for CCX gates
/// whose controls are provably equal/complementary (Keep otherwise) plus counts.
fn analyze_affine(
    ops: &[Op],
    num_q: usize,
    num_b: usize,
    input_qubits: &[QubitId],
) -> (Vec<Decision>, usize, usize) {
    // Bits are seeded Unknown (the classical input registers reg2/reg3 carry
    // per-shot data, and we make no assumption about any bit). This only ever
    // makes a condition "maybe false", collapsing the affected qubit write to a
    // fresh variable — strictly conservative, never an unsound relation claim.
    let mut af = Affine {
        cst: vec![false; num_q],
        set: vec![Vec::new(); num_q],
        nextvar: 0,
        cond_stack: Vec::new(),
        b: vec![Unknown; num_b],
    };
    for &q in input_qubits {
        let v = af.fresh();
        af.set[q.0 as usize] = v;
    }

    let mut decisions = vec![Decision::Keep; ops.len()];
    let mut fold_eq = 0usize;
    let mut drop_comp = 0usize;

    for (i, op) in ops.iter().enumerate() {
        match op.kind {
            OperationType::PushCondition => af.cond_stack.push(op.c_condition),
            OperationType::PopCondition => {
                af.cond_stack.pop();
            }
            OperationType::X => {
                let t = op.q_target.0 as usize;
                if af.cond_maybe_false(op) {
                    af.set[t] = af.fresh();
                    af.cst[t] = false;
                } else {
                    af.cst[t] ^= true;
                }
            }
            OperationType::CX => {
                let c = op.q_control1.0 as usize;
                let t = op.q_target.0 as usize;
                if af.cond_maybe_false(op) {
                    af.set[t] = af.fresh();
                    af.cst[t] = false;
                } else {
                    let ns = xor_set(&af.set[t], &af.set[c]);
                    af.cst[t] ^= af.cst[c];
                    if ns.len() > CAP_SET {
                        af.set[t] = af.fresh();
                        af.cst[t] = false;
                    } else {
                        af.set[t] = ns;
                    }
                }
            }
            OperationType::CCX => {
                let a = op.q_control1.0 as usize;
                let b = op.q_control2.0 as usize;
                let t = op.q_target.0 as usize;
                // Prove the control relationship (independent of this CCX's own
                // condition: the relation holds on every shot).
                if af.set[a] == af.set[b] {
                    if af.cst[a] == af.cst[b] {
                        decisions[i] = Decision::FoldEqualCtrls {
                            a: op.q_control1,
                            b: op.q_control2,
                            keep_ctrl: op.q_control1,
                        };
                        fold_eq += 1;
                    } else {
                        decisions[i] = Decision::DropComplementCtrls {
                            a: op.q_control1,
                            b: op.q_control2,
                        };
                        drop_comp += 1;
                    }
                }
                // Target becomes a fresh (nonlinear) variable. If we proved the
                // controls equal, the gate behaves as CX(a,t) -> t ^= a; but a
                // fresh var is a sound over-approximation of the resulting value
                // and keeps the analysis simple/exact for relations.
                af.set[t] = af.fresh();
                af.cst[t] = false;
            }
            OperationType::Swap => {
                let x = op.q_control1.0 as usize;
                let y = op.q_target.0 as usize;
                if af.cond_maybe_false(op) {
                    af.set[x] = af.fresh();
                    af.cst[x] = false;
                    af.set[y] = af.fresh();
                    af.cst[y] = false;
                } else {
                    af.set.swap(x, y);
                    af.cst.swap(x, y);
                }
            }
            OperationType::R => {
                let t = op.q_target.0 as usize;
                if af.cond_maybe_false(op) {
                    af.set[t] = af.fresh();
                    af.cst[t] = false;
                } else {
                    af.set[t] = Vec::new();
                    af.cst[t] = false;
                }
            }
            OperationType::Hmr => {
                let t = op.q_target.0 as usize;
                af.set[t] = af.fresh();
                af.cst[t] = false;
                if op.c_target != NO_BIT {
                    af.b[op.c_target.0 as usize] = Unknown;
                }
            }
            // Bit lattice tracking so condition-always-true is known. Mirror the
            // const-prop bit transfer (conservative).
            OperationType::BitStore0 => {
                if op.c_target != NO_BIT {
                    let cur = af.bv(op.c_target);
                    af.b[op.c_target.0 as usize] =
                        if af.cond_maybe_false(op) { merge(cur, Zero) } else { Zero };
                }
            }
            OperationType::BitStore1 => {
                if op.c_target != NO_BIT {
                    let cur = af.bv(op.c_target);
                    af.b[op.c_target.0 as usize] =
                        if af.cond_maybe_false(op) { merge(cur, One) } else { One };
                }
            }
            OperationType::BitInvert => {
                if op.c_target != NO_BIT {
                    let cur = af.bv(op.c_target);
                    let nv = xor_val(cur, One);
                    af.b[op.c_target.0 as usize] =
                        if af.cond_maybe_false(op) { merge(cur, nv) } else { nv };
                }
            }
            // Phase-only / metadata: no effect on computational-basis values.
            OperationType::Z
            | OperationType::CZ
            | OperationType::CCZ
            | OperationType::Neg
            | OperationType::Register
            | OperationType::AppendToRegister
            | OperationType::DebugPrint => {}
        }
    }

    (decisions, fold_eq, drop_comp)
}

/// Apply the static decisions, rewriting the op-stream. Returns the new stream.
fn apply_decisions(ops: &[Op], decisions: &[Decision]) -> Vec<Op> {
    let mut out = Vec::with_capacity(ops.len());
    for (i, op) in ops.iter().enumerate() {
        match decisions[i] {
            Decision::Keep => out.push(*op),
            Decision::DropZeroCtrl { .. } => { /* drop entirely */ }
            Decision::FoldCx { keep_ctrl, .. } => {
                let mut nop = Op::empty();
                nop.kind = OperationType::CX;
                nop.q_control1 = keep_ctrl;
                nop.q_target = op.q_target;
                nop.c_condition = op.c_condition;
                out.push(nop);
            }
            Decision::FoldX { .. } => {
                let mut nop = Op::empty();
                nop.kind = OperationType::X;
                nop.q_target = op.q_target;
                nop.c_condition = op.c_condition;
                out.push(nop);
            }
            Decision::DropComplementCtrls { .. } => { /* a&b==0: drop entirely */ }
            Decision::FoldEqualCtrls { keep_ctrl, .. } => {
                // a==b so a&b==a: CCX(a,b,t) == CX(a,t).
                let mut nop = Op::empty();
                nop.kind = OperationType::CX;
                nop.q_control1 = keep_ctrl;
                nop.q_target = op.q_target;
                nop.c_condition = op.c_condition;
                out.push(nop);
            }
        }
    }
    out
}

// ── Inverse-pair (self-inverse) CCX cancellation ─────────────────────────────
//
// CCX is self-inverse. Two CCX gates G1@i, G2@j (i<j) with the SAME target `t`
// and the SAME (unordered) control pair {a,b} compose to identity on every shot
// PROVIDED that, at the instant each executes, the contribution
//     v = cond & q[a] & q[b]
// they XOR into the target is identical AND nothing in between observes the
// half-applied target. The simulator applies CCX as `q[t] ^= cond & q[a] & q[b]`
// (sim.rs), so the pair is a true identity when, for every op strictly between i
// and j:
//   (1) `t` is neither written NOR read  — removing both XORs must not change any
//       intermediate value that depends on q[t];
//   (2) `a` and `b` are not WRITTEN       — so q[a],q[b] reaching G2 equal those
//       reaching G1 (reads of a/b are fine: we don't alter them);
//   (3) the effective condition mask is identical at i and j — guaranteed by
//       requiring no PushCondition/PopCondition between them (condition-stack
//       epoch unchanged), the identical `c_condition` field, and no write to that
//       `c_condition` bit in between.
// Under (1)-(3) the two gates cancel exactly and BOTH are removed, eliminating
// two counted Toffolis (the scorer counts each CCX by its classical condition
// mask, sim.rs:86, regardless of the quantum control values).
//
// Soundness is by construction; it is additionally corroborated empirically by
// the CONSTPROP_VERIFY harness, which replays the ORIGINAL stream and asserts
// that for every cancelled pair the live target register is bit-identical just
// before G1 and just before G2 (so the XOR'd contribution is provably equal) and
// that no intervening op touched the shared qubits.

/// One self-inverse CCX-pair cancellation: both ops at the given indices are
/// dropped.
#[derive(Clone, Copy, Debug)]
struct PairKill {
    first: usize,
    second: usize,
}

/// Forward scan that finds sound adjacent-in-dependency self-inverse CCX pairs.
/// Returns the list of (first,second) index pairs to drop and the count.
fn find_inverse_pairs(ops: &[Op], num_q: usize, num_b: usize) -> Vec<PairKill> {
    // Monotonic op index doubles as a timestamp. For each qubit we record the
    // last index that WROTE it and the last that READ it; for each bit the last
    // index that WROTE it. A pending CCX recorded at index `p` survives as a
    // cancellation candidate only while its operands' last-touch stamps stay <=p.
    let mut wlast_q = vec![usize::MAX; num_q]; // MAX sentinel = "written at init"? no
    let mut rlast_q = vec![usize::MAX; num_q];
    let mut wlast_b = vec![usize::MAX; num_b];
    // Use 0 as "never touched"; shift indices by +1 so a real touch is >=1 and
    // the MAX sentinel is unreachable. Simpler: init to a value that is never a
    // valid "between" — use a separate "never" marker of usize::MAX meaning the
    // qubit has not been touched yet, treated as <= any p (touched before all).
    for v in wlast_q.iter_mut() { *v = NEVER; }
    for v in rlast_q.iter_mut() { *v = NEVER; }
    for v in wlast_b.iter_mut() { *v = NEVER; }

    // Pending CCX candidate per target qubit: the most recent CCX whose target is
    // this qubit and which has not yet been touched/invalidated. Keyed by target.
    // Stored: (index, ctrl_a, ctrl_b sorted, c_condition, cond_epoch).
    #[derive(Clone, Copy)]
    struct Pending {
        idx: usize,
        a: u64,
        b: u64,
        cb: u64,
        epoch: u64,
    }
    let mut pending: Vec<Option<Pending>> = vec![None; num_q];

    let mut cond_epoch: u64 = 0;
    // Bit ids currently on the condition stack. The effective condition mask is
    // (AND of these bits' live values) & c_condition. Even with the epoch
    // unchanged, a write to any of these bits BETWEEN the pair would change the
    // mask, so we additionally require none of them was written after `p`.
    let mut cond_stack: Vec<u64> = Vec::new();
    let mut killed = vec![false; ops.len()];
    let mut pairs = Vec::new();

    // helper: is stamp s strictly AFTER p? (touched between p and now)
    #[inline]
    fn touched_after(s: usize, p: usize) -> bool {
        s != NEVER && s > p
    }

    for (i, op) in ops.iter().enumerate() {
        match op.kind {
            OperationType::PushCondition => {
                cond_epoch += 1;
                cond_stack.push(op.c_condition.0);
            }
            OperationType::PopCondition => {
                cond_epoch += 1;
                cond_stack.pop();
            }
            OperationType::CCX => {
                let c1 = op.q_control1.0;
                let c2 = op.q_control2.0;
                let t = op.q_target.0;
                let (a, b) = if c1 <= c2 { (c1, c2) } else { (c2, c1) };
                let cb = op.c_condition.0;

                // Try to cancel against a pending CCX on the same target.
                let mut cancelled = false;
                if let Some(p) = pending[t as usize] {
                    let same_gate = p.a == a && p.b == b && p.cb == cb;
                    let same_epoch = p.epoch == cond_epoch;
                    // controls not written between; target not touched between;
                    // condition bit not written between.
                    let ctrls_clean = !touched_after(wlast_q[a as usize], p.idx)
                        && !touched_after(wlast_q[b as usize], p.idx);
                    let tgt_clean = !touched_after(wlast_q[t as usize], p.idx)
                        && !touched_after(rlast_q[t as usize], p.idx);
                    let cond_clean = cb == u64::MAX
                        || !touched_after(wlast_b[cb as usize], p.idx);
                    // No condition-stack bit written between p and i. (same_epoch
                    // already guarantees the stack membership is identical.)
                    let stack_clean = same_epoch
                        && cond_stack
                            .iter()
                            .all(|&sb| sb == u64::MAX || !touched_after(wlast_b[sb as usize], p.idx));
                    if same_gate && same_epoch && ctrls_clean && tgt_clean && cond_clean && stack_clean {
                        killed[p.idx] = true;
                        killed[i] = true;
                        pairs.push(PairKill { first: p.idx, second: i });
                        pending[t as usize] = None;
                        cancelled = true;
                    }
                }

                if !cancelled {
                    // This CCX reads a,b and writes t. Record touches AND make it
                    // the new pending candidate for target t.
                    // (Record touches for the operands as of this op.)
                    rlast_q[a as usize] = i;
                    rlast_q[b as usize] = i;
                    wlast_q[t as usize] = i;
                    if cb != u64::MAX {
                        // condition bit is read, not written — no wlast update.
                    }
                    pending[t as usize] = Some(Pending {
                        idx: i,
                        a,
                        b,
                        cb,
                        epoch: cond_epoch,
                    });
                } else {
                    // Both gates removed. The operands were untouched between, so
                    // their last-touch stamps remain whatever they were before p
                    // (the removed pair contributes no surviving touch). Leave the
                    // stamps as-is; conservative correctness is preserved because
                    // any future op still sees stamps <= p which is fine.
                }
            }
            OperationType::CX => {
                rlast_q[op.q_control1.0 as usize] = i;
                wlast_q[op.q_target.0 as usize] = i;
                pending[op.q_target.0 as usize] = None; // target overwritten
            }
            OperationType::X => {
                wlast_q[op.q_target.0 as usize] = i;
                pending[op.q_target.0 as usize] = None;
            }
            OperationType::Swap => {
                let x = op.q_control1.0 as usize;
                let y = op.q_target.0 as usize;
                rlast_q[x] = i; rlast_q[y] = i;
                wlast_q[x] = i; wlast_q[y] = i;
                pending[x] = None;
                pending[y] = None;
            }
            OperationType::R => {
                wlast_q[op.q_target.0 as usize] = i;
                pending[op.q_target.0 as usize] = None;
            }
            OperationType::Hmr => {
                wlast_q[op.q_target.0 as usize] = i;
                if op.c_target.0 != u64::MAX { wlast_b[op.c_target.0 as usize] = i; }
                pending[op.q_target.0 as usize] = None;
            }
            OperationType::CCZ => {
                // Reads three qubits (phase only, no write). Reads invalidate the
                // target-untouched requirement, so mark all three as read.
                rlast_q[op.q_control1.0 as usize] = i;
                rlast_q[op.q_control2.0 as usize] = i;
                rlast_q[op.q_target.0 as usize] = i;
            }
            OperationType::CZ => {
                rlast_q[op.q_control1.0 as usize] = i;
                rlast_q[op.q_target.0 as usize] = i;
            }
            OperationType::Z => {
                rlast_q[op.q_target.0 as usize] = i;
            }
            OperationType::BitInvert
            | OperationType::BitStore0
            | OperationType::BitStore1 => {
                if op.c_target.0 != u64::MAX { wlast_b[op.c_target.0 as usize] = i; }
            }
            OperationType::Neg
            | OperationType::Register
            | OperationType::AppendToRegister
            | OperationType::DebugPrint => {}
        }
    }

    pairs
}

/// Public entry: run the sound const-prop peephole over `ops`.
///
/// `input_qubits` = the qubit ids that hold per-shot input data (reg0 + reg1);
/// every other qubit id is a |0>-seeded ancilla.
pub fn run(ops: Vec<Op>, input_qubits: &[QubitId]) -> Vec<Op> {
    let (num_q, num_b) = dims(&ops);
    let nonces_verify = std::env::var("CONSTPROP_VERIFY")
        .ok()
        .and_then(|s| s.parse::<usize>().ok());
    // When set, skip the (already board-verified) first-wave const-prop
    // verification and only empirically check the NEW transforms: cascade
    // const-prop drops (iter>=2) and all affine relations. Lets a focused
    // verification run quickly.
    let verify_new_only = std::env::var("CONSTPROP_VERIFY_NEW_ONLY").ok().as_deref() == Some("1");

    // Fixpoint loop: const-prop, then inverse-pair cancellation. Either pass can
    // expose new opportunities for the other (a dropped/folded gate can make a
    // downstream qubit provably constant; removing a CCX can bring an earlier and
    // later CCX into cancellation adjacency). Iterate until a full sweep makes no
    // transform. Each individual pass is sound on its current input by
    // construction, so the composition is sound.
    let mut cur = ops;
    let mut iter = 0usize;
    let mut tot_dropped = 0usize;
    let mut tot_folded_cx = 0usize;
    let mut tot_folded_x = 0usize;
    let mut tot_pairs = 0usize;
    let mut tot_aff_drop = 0usize;
    let mut tot_aff_fold = 0usize;
    let affine_disabled = std::env::var("CONSTPROP_AFFINE_DISABLE").ok().as_deref() == Some("1");
    let max_iters = std::env::var("CONSTPROP_MAX_ITERS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(16);

    loop {
        iter += 1;

        // ── const-prop pass ──
        let (mut decisions, stats) = analyze(&cur, num_q, num_b, input_qubits);

        // Empirical corroboration of the const-prop decisions on the CURRENT
        // (this-iteration input) stream, run on every iteration that produces any
        // transform — so cascade drops exposed by later iterations are verified
        // too, not just the first wave.
        if let Some(nonces) = nonces_verify {
            if stats.dropped + stats.folded_cx + stats.folded_x > 0
                && !(verify_new_only && iter == 1)
            {
                let surviving = verify_control_constancy(&cur, &decisions, num_q, num_b, nonces);
                let mut kept = 0usize;
                let mut killed = 0usize;
                for (i, ok) in surviving.iter().enumerate() {
                    if !matches!(decisions[i], Decision::Keep) {
                        if *ok {
                            kept += 1;
                        } else {
                            killed += 1;
                            decisions[i] = Decision::Keep;
                        }
                    }
                }
                eprintln!(
                    "CONSTPROP_VERIFY iter={} nonces={} shots_each=9024 transforms_static={} passed_empirical={} REVERTED_unsound={}",
                    iter,
                    nonces,
                    stats.dropped + stats.folded_cx + stats.folded_x,
                    kept,
                    killed
                );
            }
        }

        let cp_transforms = stats.dropped + stats.folded_cx + stats.folded_x;
        tot_dropped += stats.dropped;
        tot_folded_cx += stats.folded_cx;
        tot_folded_x += stats.folded_x;
        cur = apply_decisions(&cur, &decisions);

        // ── inverse-pair cancellation pass ──
        let (nq2, nb2) = dims(&cur);
        let pairs = find_inverse_pairs(&cur, nq2, nb2);

        // Empirical corroboration of the pair cancellations on this stream.
        if let Some(nonces) = nonces_verify {
            if !pairs.is_empty() {
                let bad = verify_inverse_pairs(&cur, &pairs, nq2, nb2, nonces);
                eprintln!(
                    "CONSTPROP_PAIR_VERIFY iter={} nonces={} pairs={} UNSOUND_pairs={}",
                    iter,
                    nonces,
                    pairs.len(),
                    bad
                );
                if bad != 0 {
                    panic!(
                        "INVERSE-PAIR CANCELLATION UNSOUND: {} of {} pairs failed empirical check",
                        bad,
                        pairs.len()
                    );
                }
            }
        }

        let pair_transforms = pairs.len();
        tot_pairs += pair_transforms;
        if pair_transforms > 0 {
            let mut kill = vec![false; cur.len()];
            for p in &pairs {
                kill[p.first] = true;
                kill[p.second] = true;
            }
            let mut out = Vec::with_capacity(cur.len() - 2 * pair_transforms);
            for (i, op) in cur.iter().enumerate() {
                if !kill[i] {
                    out.push(*op);
                }
            }
            cur = out;
        }

        // ── affine-relation pass (equal/complementary controls) ──
        let (mut aff_drop, mut aff_fold) = (0usize, 0usize);
        if !affine_disabled {
            let (nq3, nb3) = dims(&cur);
            let (mut adec, fold_eq, drop_comp) =
                analyze_affine(&cur, nq3, nb3, input_qubits);

            // Empirical corroboration of the affine relationship claims.
            if let Some(nonces) = nonces_verify {
                if fold_eq + drop_comp > 0 {
                    let surviving =
                        verify_affine_relations(&cur, &adec, nq3, nb3, nonces);
                    let mut killed = 0usize;
                    for (i, ok) in surviving.iter().enumerate() {
                        if matches!(
                            adec[i],
                            Decision::DropComplementCtrls { .. }
                                | Decision::FoldEqualCtrls { .. }
                        ) && !*ok
                        {
                            killed += 1;
                            adec[i] = Decision::Keep;
                        }
                    }
                    eprintln!(
                        "CONSTPROP_AFFINE_VERIFY iter={} nonces={} fold_eq={} drop_comp={} REVERTED_unsound={}",
                        iter, nonces, fold_eq, drop_comp, killed
                    );
                    if killed != 0 {
                        panic!(
                            "AFFINE RELATION CLAIM UNSOUND: {} flagged CCX failed empirical check",
                            killed
                        );
                    }
                }
            }

            // Recount after any reversion and rewrite.
            for d in &adec {
                match d {
                    Decision::DropComplementCtrls { .. } => aff_drop += 1,
                    Decision::FoldEqualCtrls { .. } => aff_fold += 1,
                    _ => {}
                }
            }
            if aff_drop + aff_fold > 0 {
                cur = apply_decisions(&cur, &adec);
            }
            let _ = (fold_eq, drop_comp);
        }
        tot_aff_drop += aff_drop;
        tot_aff_fold += aff_fold;

        eprintln!(
            "CONSTPROP iter={} ccx_total={} dropped={} folded_cx={} folded_x={} inverse_pairs={} aff_drop={} aff_fold={} (this-iter toffoli removed = {})",
            iter,
            stats.ccx_total,
            stats.dropped,
            stats.folded_cx,
            stats.folded_x,
            pair_transforms,
            aff_drop,
            aff_fold,
            cp_transforms + 2 * pair_transforms + aff_drop + aff_fold,
        );

        if cp_transforms == 0 && pair_transforms == 0 && aff_drop + aff_fold == 0 {
            break;
        }
        if iter >= max_iters {
            eprintln!("CONSTPROP reached max_iters={}, stopping", max_iters);
            break;
        }
    }

    eprintln!(
        "CONSTPROP TOTAL iters={} dropped={} folded_cx={} folded_x={} inverse_pairs={} aff_drop={} aff_fold={} (toffoli removed = {})",
        iter,
        tot_dropped,
        tot_folded_cx,
        tot_folded_x,
        tot_pairs,
        tot_aff_drop,
        tot_aff_fold,
        tot_dropped + tot_folded_cx + tot_folded_x + 2 * tot_pairs + tot_aff_drop + tot_aff_fold,
    );

    cur
}

fn dims(ops: &[Op]) -> (usize, usize) {
    let mut nq = 0u64;
    let mut nb = 0u64;
    for op in ops {
        for q in [op.q_control2, op.q_control1, op.q_target] {
            if q != NO_QUBIT {
                nq = nq.max(q.0 + 1);
            }
        }
        for b in [op.c_target, op.c_condition] {
            if b != NO_BIT {
                nb = nb.max(b.0 + 1);
            }
        }
    }
    (nq as usize, nb as usize)
}

/// Empirically verify, over `nonces` diverse Fiat-Shamir seeds x 9024 shots,
/// that every transformed gate's claimed-constant control(s) are ALWAYS that
/// constant. Returns a per-op bool: true = the static claim held on all inputs
/// (or the op was Keep), false = a counter-example was observed (UNSOUND claim).
///
/// We re-derive the verifier's exact input distribution and run the *original*
/// op-stream through the real `Simulator`, snapshotting the relevant qubit
/// register just before each flagged op.
fn verify_control_constancy(
    ops: &[Op],
    decisions: &[Decision],
    num_q: usize,
    num_b: usize,
    nonces: usize,
) -> Vec<bool> {
    use crate::circuit::{analyze_ops, QubitOrBit};
    use crate::sim::Simulator;
    use crate::weierstrass_elliptic_curve::WeierstrassEllipticCurve;
    use alloy_primitives::U256;
    use sha3::{digest::{ExtendableOutput, Update, XofReader}, Shake256};

    let curve = WeierstrassEllipticCurve {
        modulus: U256::from_str_radix("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16).unwrap(),
        a: U256::from(0u64),
        b: U256::from(7u64),
        gx: U256::from_str_radix("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16).unwrap(),
        gy: U256::from_str_radix("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16).unwrap(),
        order: U256::from_str_radix("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141", 16).unwrap(),
    };

    let (_tq, _tb, _nr, regs) = analyze_ops(ops.iter());
    assert_eq!(regs.len(), 4, "expected 4 IO registers");

    // Indices of flagged ops and the control qubit(s) + claimed value to verify.
    // claim: Vec<(qubit, expected_val_bit)>  (expected 0 or 1)
    let mut flagged: Vec<(usize, Vec<(QubitId, u64)>)> = Vec::new();
    for (i, d) in decisions.iter().enumerate() {
        match *d {
            Decision::Keep => {}
            Decision::DropZeroCtrl { ctrl } => flagged.push((i, vec![(ctrl, 0)])),
            Decision::FoldCx { one_ctrl, .. } => flagged.push((i, vec![(one_ctrl, 1)])),
            Decision::FoldX { c1, c2 } => flagged.push((i, vec![(c1, 1), (c2, 1)])),
            // Affine-relation decisions are verified separately by
            // verify_affine_relations; not produced for this verifier.
            Decision::DropComplementCtrls { .. } | Decision::FoldEqualCtrls { .. } => {}
        }
    }
    let mut ok = vec![true; ops.len()];
    if flagged.is_empty() {
        return ok;
    }

    // Flat op-index -> position-in-`flagged` (u32::MAX = not flagged). A flat
    // array index is far cheaper than a HashMap on the 10.7M-op hot path.
    let mut flag_pos = vec![u32::MAX; ops.len()];
    for (p, (i, _)) in flagged.iter().enumerate() {
        flag_pos[*i] = p as u32;
    }

    const NUM_TESTS: usize = 9024;
    const BATCH: usize = 64;

    for nonce in 0..nonces {
        // Diverse seed per nonce: mix the official domain tag with the nonce.
        let mut hasher = Shake256::default();
        hasher.update(b"quantum_ecc-fiat-shamir-v2");
        hasher.update(&(ops.len() as u64).to_le_bytes());
        // Perturb with nonce so each run draws a different input population.
        hasher.update(b"CONSTPROP_VERIFY");
        hasher.update(&(nonce as u64).to_le_bytes());
        let mut xof = hasher.finalize_xof();

        let mut targets = Vec::new();
        let mut offsets = Vec::new();
        for _ in 0..NUM_TESTS {
            let mut rb = [[0u8; 32]; 2];
            xof.read(&mut rb[0]);
            xof.read(&mut rb[1]);
            let k1 = U256::from_le_bytes(rb[0]);
            let k2 = U256::from_le_bytes(rb[1]);
            let t = curve.mul(curve.gx, curve.gy, k1);
            let o = curve.mul(curve.gx, curve.gy, k2);
            if t.0 == o.0 { continue; }
            if t.0.is_zero() && t.1.is_zero() { continue; }
            if o.0.is_zero() && o.1.is_zero() { continue; }
            targets.push(t);
            offsets.push(o);
        }
        let n = targets.len();
        let num_batches = (n + BATCH - 1) / BATCH;

        let mut sim = Simulator::new(num_q, num_b, &mut xof);
        for batch in 0..num_batches {
            let bs = BATCH.min(n - batch * BATCH);
            sim.clear_for_shot();
            for shot in 0..bs {
                let i = batch * BATCH + shot;
                sim.set_register(&regs[0], targets[i].0, shot);
                sim.set_register(&regs[1], targets[i].1, shot);
                sim.set_register(&regs[2], offsets[i].0, shot);
                sim.set_register(&regs[3], offsets[i].1, shot);
            }
            let cond_mask: u64 = if bs == 64 { u64::MAX } else { (1u64 << bs) - 1 };

            // Step the simulator op by op, checking claimed controls just before
            // each flagged op executes.
            step_and_check(&mut sim, ops, &flag_pos, &flagged, &mut ok, cond_mask);
        }
        let bad = ok.iter().filter(|b| !**b).count();
        eprintln!(
            "CONSTPROP_PROGRESS nonce={}/{} shots={} cumulative_failed_claims={}",
            nonce + 1, nonces, n, bad
        );
    }
    let _ = QubitOrBit::Bit; // silence potential unused import lints
    ok
}

/// Single-batch op-by-op driver that mirrors `Simulator::apply_iter` exactly,
/// but snapshots the claimed control qubit values just before each flagged op.
fn step_and_check<R: sha3::digest::XofReader>(
    sim: &mut crate::sim::Simulator<R>,
    ops: &[Op],
    flag_pos: &[u32],
    flagged: &[(usize, Vec<(QubitId, u64)>)],
    ok: &mut [bool],
    cond_mask: u64,
) {
    // We cannot reuse apply_iter (it consumes the whole iter), so replicate the
    // condition-stack machinery and check before each op.
    let mut condition_stack: Vec<u64> = Vec::new();
    let mut current_base_condition = u64::MAX;

    for (idx, op) in ops.iter().enumerate() {
        // Check claim BEFORE executing this op.
        let fp = flag_pos[idx];
        if fp != u32::MAX {
            let p = fp as usize;
            for &(qid, expected) in &flagged[p].1 {
                let live = sim.qubit(qid) & cond_mask;
                let claim_ok = if expected == 0 {
                    // claimed Zero: no live shot has the bit set.
                    live == 0
                } else {
                    // claimed One: every live shot has the bit set.
                    live == cond_mask
                };
                if !claim_ok {
                    ok[idx] = false;
                }
            }
        }

        // Now replicate the simulator step for this single op.
        let mut cond = current_base_condition;
        if op.c_condition != NO_BIT {
            cond &= sim.bit(op.c_condition);
        }
        match op.kind {
            OperationType::CCX => {
                let v = cond & sim.qubit(op.q_control1) & sim.qubit(op.q_control2);
                *sim.qubit_mut(op.q_target) ^= v;
            }
            OperationType::CX => {
                let v = cond & sim.qubit(op.q_control1);
                *sim.qubit_mut(op.q_target) ^= v;
            }
            OperationType::Swap => {
                let mut q_c1 = sim.qubit(op.q_control1);
                let mut q_t = sim.qubit(op.q_target);
                q_c1 ^= q_t;
                q_t ^= cond & q_c1;
                q_c1 ^= q_t;
                *sim.qubit_mut(op.q_control1) = q_c1;
                *sim.qubit_mut(op.q_target) = q_t;
            }
            OperationType::X => {
                *sim.qubit_mut(op.q_target) ^= cond;
            }
            OperationType::CCZ => {
                let v = cond & sim.qubit(op.q_target) & sim.qubit(op.q_control1) & sim.qubit(op.q_control2);
                sim.phase ^= v;
            }
            OperationType::CZ => {
                let v = cond & sim.qubit(op.q_target) & sim.qubit(op.q_control1);
                sim.phase ^= v;
            }
            OperationType::Z => {
                let v = cond & sim.qubit(op.q_target);
                sim.phase ^= v;
            }
            OperationType::Neg => {
                sim.phase ^= cond;
            }
            OperationType::Hmr => {
                let mut buf = [0u8; 8];
                sim.xof.read(&mut buf);
                let rng_val = u64::from_le_bytes(buf);
                *sim.bit_mut(op.c_target) &= !cond;
                *sim.bit_mut(op.c_target) ^= rng_val & cond;
                sim.phase ^= sim.qubit(op.q_target) & rng_val & cond;
                *sim.qubit_mut(op.q_target) &= !cond;
            }
            OperationType::R => {
                let mut buf = [0u8; 8];
                sim.xof.read(&mut buf);
                let rng_val = u64::from_le_bytes(buf);
                sim.phase ^= sim.qubit(op.q_target) & rng_val & cond;
                *sim.qubit_mut(op.q_target) &= !cond;
            }
            OperationType::BitInvert => {
                *sim.bit_mut(op.c_target) ^= cond;
            }
            OperationType::BitStore0 => {
                *sim.bit_mut(op.c_target) &= !cond;
            }
            OperationType::BitStore1 => {
                *sim.bit_mut(op.c_target) |= cond;
            }
            OperationType::AppendToRegister
            | OperationType::Register
            | OperationType::DebugPrint => {}
            OperationType::PushCondition => {
                condition_stack.push(current_base_condition);
                current_base_condition &= sim.bit(op.c_condition);
            }
            OperationType::PopCondition => {
                if let Some(val) = condition_stack.pop() {
                    current_base_condition = val;
                }
            }
        }
    }
}

/// Empirically corroborate inverse-pair cancellations over `nonces` x 9024 shots.
/// For each cancelled pair (G1@first, G2@second) we replay the ORIGINAL stream
/// and assert, across all live shots, that:
///   (a) the XOR contribution `cond & q[a] & q[b]` is bit-identical at G1 and G2
///       (so the two CCX truly XOR the same value into the target), AND
///   (b) the target register q[t] is bit-identical just before G1 and just before
///       G2 (so no intervening op observed a half-applied target).
/// Either failing on ANY shot marks the pair unsound. Returns the count of
/// unsound pairs (must be 0).
fn verify_inverse_pairs(
    ops: &[Op],
    pairs: &[PairKill],
    num_q: usize,
    num_b: usize,
    nonces: usize,
) -> usize {
    use crate::circuit::analyze_ops;
    use crate::sim::Simulator;
    use crate::weierstrass_elliptic_curve::WeierstrassEllipticCurve;
    use alloy_primitives::U256;
    use sha3::{digest::{ExtendableOutput, Update, XofReader}, Shake256};

    let curve = WeierstrassEllipticCurve {
        modulus: U256::from_str_radix("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16).unwrap(),
        a: U256::from(0u64),
        b: U256::from(7u64),
        gx: U256::from_str_radix("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16).unwrap(),
        gy: U256::from_str_radix("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16).unwrap(),
        order: U256::from_str_radix("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141", 16).unwrap(),
    };

    let (_tq, _tb, _nr, regs) = analyze_ops(ops.iter());
    assert_eq!(regs.len(), 4, "expected 4 IO registers");

    // For each op index that is a pair endpoint, record (pair_pos, is_first).
    // We snapshot the pair's contribution & target at each endpoint.
    let mut endpoint: Vec<u32> = vec![u32::MAX; ops.len()]; // pair index
    let mut is_first_at: Vec<bool> = vec![false; ops.len()];
    for (p, pk) in pairs.iter().enumerate() {
        endpoint[pk.first] = p as u32;
        is_first_at[pk.first] = true;
        endpoint[pk.second] = p as u32;
        is_first_at[pk.second] = false;
    }

    // Per-pair accumulated mismatch flag across all batches/nonces.
    let mut bad_pair = vec![false; pairs.len()];

    const NUM_TESTS: usize = 9024;
    const BATCH: usize = 64;

    for nonce in 0..nonces {
        let mut hasher = Shake256::default();
        hasher.update(b"quantum_ecc-fiat-shamir-v2");
        hasher.update(&(ops.len() as u64).to_le_bytes());
        hasher.update(b"CONSTPROP_PAIR_VERIFY");
        hasher.update(&(nonce as u64).to_le_bytes());
        let mut xof = hasher.finalize_xof();

        let mut targets = Vec::new();
        let mut offsets = Vec::new();
        for _ in 0..NUM_TESTS {
            let mut rb = [[0u8; 32]; 2];
            xof.read(&mut rb[0]);
            xof.read(&mut rb[1]);
            let k1 = U256::from_le_bytes(rb[0]);
            let k2 = U256::from_le_bytes(rb[1]);
            let t = curve.mul(curve.gx, curve.gy, k1);
            let o = curve.mul(curve.gx, curve.gy, k2);
            if t.0 == o.0 { continue; }
            if t.0.is_zero() && t.1.is_zero() { continue; }
            if o.0.is_zero() && o.1.is_zero() { continue; }
            targets.push(t);
            offsets.push(o);
        }
        let n = targets.len();
        let num_batches = (n + BATCH - 1) / BATCH;

        let mut sim = Simulator::new(num_q, num_b, &mut xof);
        // Per-pair snapshot of the FIRST endpoint, valid within the current
        // batch: (contribution mask, target value). Indexed by pair position.
        let mut snap_contrib = vec![0u64; pairs.len()];
        let mut snap_tgt = vec![0u64; pairs.len()];
        let mut snap_seen = vec![false; pairs.len()];

        for batch in 0..num_batches {
            let bs = BATCH.min(n - batch * BATCH);
            sim.clear_for_shot();
            for shot in 0..bs {
                let i = batch * BATCH + shot;
                sim.set_register(&regs[0], targets[i].0, shot);
                sim.set_register(&regs[1], targets[i].1, shot);
                sim.set_register(&regs[2], offsets[i].0, shot);
                sim.set_register(&regs[3], offsets[i].1, shot);
            }
            let cond_mask: u64 = if bs == 64 { u64::MAX } else { (1u64 << bs) - 1 };
            for s in snap_seen.iter_mut() { *s = false; }

            step_and_check_pairs(
                &mut sim,
                ops,
                pairs,
                &endpoint,
                &is_first_at,
                &mut snap_contrib,
                &mut snap_tgt,
                &mut snap_seen,
                &mut bad_pair,
                cond_mask,
            );
        }
        let bad = bad_pair.iter().filter(|b| **b).count();
        eprintln!(
            "CONSTPROP_PAIR_PROGRESS nonce={}/{} shots={} cumulative_unsound_pairs={}",
            nonce + 1, nonces, n, bad
        );
    }

    bad_pair.iter().filter(|b| **b).count()
}

/// Single-batch driver mirroring `Simulator::apply_iter`, snapshotting each
/// inverse pair's contribution `cond & q[a] & q[b]` and target value at its two
/// endpoints and flagging any mismatch (across live shots) into `bad_pair`.
fn step_and_check_pairs<R: sha3::digest::XofReader>(
    sim: &mut crate::sim::Simulator<R>,
    ops: &[Op],
    pairs: &[PairKill],
    endpoint: &[u32],
    is_first_at: &[bool],
    snap_contrib: &mut [u64],
    snap_tgt: &mut [u64],
    snap_seen: &mut [bool],
    bad_pair: &mut [bool],
    cond_mask: u64,
) {
    let mut condition_stack: Vec<u64> = Vec::new();
    let mut current_base_condition = u64::MAX;

    for (idx, op) in ops.iter().enumerate() {
        // Snapshot/compare BEFORE executing this op if it is a pair endpoint.
        let pp = endpoint[idx];
        if pp != u32::MAX {
            let p = pp as usize;
            // Recompute the effective condition mask for this op.
            let mut cond = current_base_condition;
            if op.c_condition != NO_BIT {
                cond &= sim.bit(op.c_condition);
            }
            let a = op.q_control1;
            let b = op.q_control2;
            let t = op.q_target;
            let contrib = (cond & sim.qubit(a) & sim.qubit(b)) & cond_mask;
            let tgt = sim.qubit(t) & cond_mask;
            if is_first_at[idx] {
                snap_contrib[p] = contrib;
                snap_tgt[p] = tgt;
                snap_seen[p] = true;
            } else if snap_seen[p] {
                if contrib != snap_contrib[p] || tgt != snap_tgt[p] {
                    bad_pair[p] = true;
                }
            } else {
                // Second endpoint seen without a first in this batch: should not
                // happen (both are in the same stream); flag conservatively.
                bad_pair[p] = true;
            }
        }

        // Replicate the simulator step for this op.
        let mut cond = current_base_condition;
        if op.c_condition != NO_BIT {
            cond &= sim.bit(op.c_condition);
        }
        match op.kind {
            OperationType::CCX => {
                let v = cond & sim.qubit(op.q_control1) & sim.qubit(op.q_control2);
                *sim.qubit_mut(op.q_target) ^= v;
            }
            OperationType::CX => {
                let v = cond & sim.qubit(op.q_control1);
                *sim.qubit_mut(op.q_target) ^= v;
            }
            OperationType::Swap => {
                let mut q_c1 = sim.qubit(op.q_control1);
                let mut q_t = sim.qubit(op.q_target);
                q_c1 ^= q_t;
                q_t ^= cond & q_c1;
                q_c1 ^= q_t;
                *sim.qubit_mut(op.q_control1) = q_c1;
                *sim.qubit_mut(op.q_target) = q_t;
            }
            OperationType::X => {
                *sim.qubit_mut(op.q_target) ^= cond;
            }
            OperationType::CCZ => {
                let v = cond & sim.qubit(op.q_target) & sim.qubit(op.q_control1) & sim.qubit(op.q_control2);
                sim.phase ^= v;
            }
            OperationType::CZ => {
                let v = cond & sim.qubit(op.q_target) & sim.qubit(op.q_control1);
                sim.phase ^= v;
            }
            OperationType::Z => {
                let v = cond & sim.qubit(op.q_target);
                sim.phase ^= v;
            }
            OperationType::Neg => {
                sim.phase ^= cond;
            }
            OperationType::Hmr => {
                let mut buf = [0u8; 8];
                sim.xof.read(&mut buf);
                let rng_val = u64::from_le_bytes(buf);
                *sim.bit_mut(op.c_target) &= !cond;
                *sim.bit_mut(op.c_target) ^= rng_val & cond;
                sim.phase ^= sim.qubit(op.q_target) & rng_val & cond;
                *sim.qubit_mut(op.q_target) &= !cond;
            }
            OperationType::R => {
                let mut buf = [0u8; 8];
                sim.xof.read(&mut buf);
                let rng_val = u64::from_le_bytes(buf);
                sim.phase ^= sim.qubit(op.q_target) & rng_val & cond;
                *sim.qubit_mut(op.q_target) &= !cond;
            }
            OperationType::BitInvert => {
                *sim.bit_mut(op.c_target) ^= cond;
            }
            OperationType::BitStore0 => {
                *sim.bit_mut(op.c_target) &= !cond;
            }
            OperationType::BitStore1 => {
                *sim.bit_mut(op.c_target) |= cond;
            }
            OperationType::AppendToRegister
            | OperationType::Register
            | OperationType::DebugPrint => {}
            OperationType::PushCondition => {
                condition_stack.push(current_base_condition);
                current_base_condition &= sim.bit(op.c_condition);
            }
            OperationType::PopCondition => {
                if let Some(val) = condition_stack.pop() {
                    current_base_condition = val;
                }
            }
        }
    }
    let _ = pairs;
}

/// Empirically corroborate affine-relation claims over `nonces` x 9024 shots.
/// For each flagged CCX, replays the ORIGINAL stream and asserts that, across all
/// live shots, the two controls are ALWAYS equal (FoldEqualCtrls) or ALWAYS
/// complementary (DropComplementCtrls). Returns a per-op bool (true = claim held
/// on all inputs, or Keep). Any false marks an UNSOUND claim.
fn verify_affine_relations(
    ops: &[Op],
    decisions: &[Decision],
    num_q: usize,
    num_b: usize,
    nonces: usize,
) -> Vec<bool> {
    use crate::circuit::analyze_ops;
    use crate::sim::Simulator;
    use crate::weierstrass_elliptic_curve::WeierstrassEllipticCurve;
    use alloy_primitives::U256;
    use sha3::{digest::{ExtendableOutput, Update, XofReader}, Shake256};

    let curve = WeierstrassEllipticCurve {
        modulus: U256::from_str_radix("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16).unwrap(),
        a: U256::from(0u64),
        b: U256::from(7u64),
        gx: U256::from_str_radix("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16).unwrap(),
        gy: U256::from_str_radix("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16).unwrap(),
        order: U256::from_str_radix("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141", 16).unwrap(),
    };

    let (_tq, _tb, _nr, regs) = analyze_ops(ops.iter());
    assert_eq!(regs.len(), 4, "expected 4 IO registers");

    // For each flagged op: (a, b, want_equal). want_equal=true -> assert a==b on
    // every live shot; false -> assert a==NOT b (a^b all-ones) on every live shot.
    let mut want_equal = vec![false; ops.len()];
    let mut flagged_idx: Vec<usize> = Vec::new();
    let mut is_flagged = vec![false; ops.len()];
    for (i, d) in decisions.iter().enumerate() {
        match *d {
            Decision::FoldEqualCtrls { .. } => {
                want_equal[i] = true;
                is_flagged[i] = true;
                flagged_idx.push(i);
            }
            Decision::DropComplementCtrls { .. } => {
                want_equal[i] = false;
                is_flagged[i] = true;
                flagged_idx.push(i);
            }
            _ => {}
        }
    }
    let mut ok = vec![true; ops.len()];
    if flagged_idx.is_empty() {
        return ok;
    }

    const NUM_TESTS: usize = 9024;
    const BATCH: usize = 64;

    for nonce in 0..nonces {
        let mut hasher = Shake256::default();
        hasher.update(b"quantum_ecc-fiat-shamir-v2");
        hasher.update(&(ops.len() as u64).to_le_bytes());
        hasher.update(b"CONSTPROP_AFFINE_VERIFY");
        hasher.update(&(nonce as u64).to_le_bytes());
        let mut xof = hasher.finalize_xof();

        let mut targets = Vec::new();
        let mut offsets = Vec::new();
        for _ in 0..NUM_TESTS {
            let mut rb = [[0u8; 32]; 2];
            xof.read(&mut rb[0]);
            xof.read(&mut rb[1]);
            let k1 = U256::from_le_bytes(rb[0]);
            let k2 = U256::from_le_bytes(rb[1]);
            let t = curve.mul(curve.gx, curve.gy, k1);
            let o = curve.mul(curve.gx, curve.gy, k2);
            if t.0 == o.0 { continue; }
            if t.0.is_zero() && t.1.is_zero() { continue; }
            if o.0.is_zero() && o.1.is_zero() { continue; }
            targets.push(t);
            offsets.push(o);
        }
        let n = targets.len();
        let num_batches = (n + BATCH - 1) / BATCH;

        let mut sim = Simulator::new(num_q, num_b, &mut xof);
        for batch in 0..num_batches {
            let bs = BATCH.min(n - batch * BATCH);
            sim.clear_for_shot();
            for shot in 0..bs {
                let i = batch * BATCH + shot;
                sim.set_register(&regs[0], targets[i].0, shot);
                sim.set_register(&regs[1], targets[i].1, shot);
                sim.set_register(&regs[2], offsets[i].0, shot);
                sim.set_register(&regs[3], offsets[i].1, shot);
            }
            let cond_mask: u64 = if bs == 64 { u64::MAX } else { (1u64 << bs) - 1 };
            step_and_check_affine(
                &mut sim,
                ops,
                &is_flagged,
                &want_equal,
                &mut ok,
                cond_mask,
            );
        }
        let bad = flagged_idx.iter().filter(|&&i| !ok[i]).count();
        eprintln!(
            "CONSTPROP_AFFINE_PROGRESS nonce={}/{} shots={} cumulative_failed_claims={}",
            nonce + 1, nonces, n, bad
        );
    }
    ok
}

/// Single-batch driver: checks affine relationship claims just before each
/// flagged CCX executes, then replicates the simulator step.
fn step_and_check_affine<R: sha3::digest::XofReader>(
    sim: &mut crate::sim::Simulator<R>,
    ops: &[Op],
    is_flagged: &[bool],
    want_equal: &[bool],
    ok: &mut [bool],
    cond_mask: u64,
) {
    let mut condition_stack: Vec<u64> = Vec::new();
    let mut current_base_condition = u64::MAX;

    for (idx, op) in ops.iter().enumerate() {
        if is_flagged[idx] {
            // Check across ALL shots (the relationship is claimed to hold on every
            // shot regardless of the gate's own condition).
            let va = sim.qubit(op.q_control1) & cond_mask;
            let vb = sim.qubit(op.q_control2) & cond_mask;
            let claim_ok = if want_equal[idx] {
                va == vb
            } else {
                (va ^ vb) == cond_mask
            };
            if !claim_ok {
                ok[idx] = false;
            }
        }

        let mut cond = current_base_condition;
        if op.c_condition != NO_BIT {
            cond &= sim.bit(op.c_condition);
        }
        match op.kind {
            OperationType::CCX => {
                let v = cond & sim.qubit(op.q_control1) & sim.qubit(op.q_control2);
                *sim.qubit_mut(op.q_target) ^= v;
            }
            OperationType::CX => {
                let v = cond & sim.qubit(op.q_control1);
                *sim.qubit_mut(op.q_target) ^= v;
            }
            OperationType::Swap => {
                let mut q_c1 = sim.qubit(op.q_control1);
                let mut q_t = sim.qubit(op.q_target);
                q_c1 ^= q_t;
                q_t ^= cond & q_c1;
                q_c1 ^= q_t;
                *sim.qubit_mut(op.q_control1) = q_c1;
                *sim.qubit_mut(op.q_target) = q_t;
            }
            OperationType::X => {
                *sim.qubit_mut(op.q_target) ^= cond;
            }
            OperationType::CCZ => {
                let v = cond & sim.qubit(op.q_target) & sim.qubit(op.q_control1) & sim.qubit(op.q_control2);
                sim.phase ^= v;
            }
            OperationType::CZ => {
                let v = cond & sim.qubit(op.q_target) & sim.qubit(op.q_control1);
                sim.phase ^= v;
            }
            OperationType::Z => {
                let v = cond & sim.qubit(op.q_target);
                sim.phase ^= v;
            }
            OperationType::Neg => {
                sim.phase ^= cond;
            }
            OperationType::Hmr => {
                let mut buf = [0u8; 8];
                sim.xof.read(&mut buf);
                let rng_val = u64::from_le_bytes(buf);
                *sim.bit_mut(op.c_target) &= !cond;
                *sim.bit_mut(op.c_target) ^= rng_val & cond;
                sim.phase ^= sim.qubit(op.q_target) & rng_val & cond;
                *sim.qubit_mut(op.q_target) &= !cond;
            }
            OperationType::R => {
                let mut buf = [0u8; 8];
                sim.xof.read(&mut buf);
                let rng_val = u64::from_le_bytes(buf);
                sim.phase ^= sim.qubit(op.q_target) & rng_val & cond;
                *sim.qubit_mut(op.q_target) &= !cond;
            }
            OperationType::BitInvert => {
                *sim.bit_mut(op.c_target) ^= cond;
            }
            OperationType::BitStore0 => {
                *sim.bit_mut(op.c_target) &= !cond;
            }
            OperationType::BitStore1 => {
                *sim.bit_mut(op.c_target) |= cond;
            }
            OperationType::AppendToRegister
            | OperationType::Register
            | OperationType::DebugPrint => {}
            OperationType::PushCondition => {
                condition_stack.push(current_base_condition);
                current_base_condition &= sim.bit(op.c_condition);
            }
            OperationType::PopCondition => {
                if let Some(val) = condition_stack.pop() {
                    current_base_condition = val;
                }
            }
        }
    }
}
