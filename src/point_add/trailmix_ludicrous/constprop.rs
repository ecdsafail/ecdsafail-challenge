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
        }
    }
    out
}

/// Public entry: run the sound const-prop peephole over `ops`.
///
/// `input_qubits` = the qubit ids that hold per-shot input data (reg0 + reg1);
/// every other qubit id is a |0>-seeded ancilla.
pub fn run(ops: Vec<Op>, input_qubits: &[QubitId]) -> Vec<Op> {
    let (num_q, num_b) = dims(&ops);
    let (mut decisions, stats) = analyze(&ops, num_q, num_b, input_qubits);

    // Optional empirical corroboration over many diverse nonces.
    if let Ok(nv) = std::env::var("CONSTPROP_VERIFY") {
        if let Ok(nonces) = nv.parse::<usize>() {
            let surviving = verify_control_constancy(&ops, &decisions, num_q, num_b, nonces);
            // Restrict to gates that pass the empirical check.
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
                "CONSTPROP_VERIFY nonces={} shots_each=9024 transforms_static={} passed_empirical={} REVERTED_unsound={}",
                nonces,
                stats.dropped + stats.folded_cx + stats.folded_x,
                kept,
                killed
            );
        }
    }

    eprintln!(
        "CONSTPROP ccx_total={} dropped={} folded_cx={} folded_x={} (toffoli removed = {})",
        stats.ccx_total,
        stats.dropped,
        stats.folded_cx,
        stats.folded_x,
        stats.dropped + stats.folded_cx + stats.folded_x
    );

    apply_decisions(&ops, &decisions)
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
