
use crate::circuit::{BitId, NO_BIT, NO_QUBIT, Op, OperationType, QubitId};
use crate::point_add::OpSite;

const NEVER: usize = usize::MAX;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Val {
    Zero,
    One,
    Unknown,
}

use Val::*;

#[derive(Clone, Copy, Debug, Default)]
pub struct ConstPropStats {
    pub ccx_total: usize,
    pub dropped: usize,
    pub folded_cx: usize,
    pub folded_x: usize,
}

#[derive(Clone, Copy, Debug)]
enum Decision {

    Keep,

    DropZeroCtrl { ctrl: QubitId },

    FoldCx { one_ctrl: QubitId, keep_ctrl: QubitId },

    FoldX { c1: QubitId, c2: QubitId },

    DropComplementCtrls { a: QubitId, b: QubitId },

    FoldEqualCtrls { a: QubitId, b: QubitId, keep_ctrl: QubitId },
}

struct Analyzer {
    q: Vec<Val>,
    b: Vec<Val>,

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

    fn cond_maybe_false(&self, op: &Op) -> bool {
        !self.cond_always_true(op)
    }
}

fn xor_val(a: Val, b: Val) -> Val {
    match (a, b) {
        (Zero, x) | (x, Zero) => x,
        (One, One) => Zero,
        _ => Unknown,
    }
}

fn and_val(a: Val, b: Val) -> Val {
    match (a, b) {
        (Zero, _) | (_, Zero) => Zero,
        (One, One) => One,
        _ => Unknown,
    }
}

fn merge(old: Val, new: Val) -> Val {
    if old == new { old } else { Unknown }
}

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

                if c1 == Zero {
                    decisions[i] = Decision::DropZeroCtrl { ctrl: op.q_control1 };
                    stats.dropped += 1;

                } else if c2 == Zero {
                    decisions[i] = Decision::DropZeroCtrl { ctrl: op.q_control2 };
                    stats.dropped += 1;

                } else if c1 == One && c2 == One {
                    decisions[i] = Decision::FoldX { c1: op.q_control1, c2: op.q_control2 };
                    stats.folded_x += 1;

                    let tgt = a.qv(op.q_target);
                    let nv = xor_val(tgt, One);
                    let res = if a.cond_maybe_false(op) { merge(tgt, nv) } else { nv };
                    a.set_q(op.q_target, res);
                } else if c1 == One {
                    decisions[i] = Decision::FoldCx { one_ctrl: op.q_control1, keep_ctrl: op.q_control2 };
                    stats.folded_cx += 1;

                    let tgt = a.qv(op.q_target);
                    let delta = c2;
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

                let va = a.qv(op.q_control1);
                let vt = a.qv(op.q_target);
                if a.cond_maybe_false(op) {

                    a.set_q(op.q_control1, merge(va, vt));
                    a.set_q(op.q_target, merge(vt, va));
                } else {
                    a.set_q(op.q_control1, vt);
                    a.set_q(op.q_target, va);
                }
            }
            OperationType::R => {

                let tgt = a.qv(op.q_target);
                let res = if a.cond_maybe_false(op) { merge(tgt, Zero) } else { Zero };
                a.set_q(op.q_target, res);
            }
            OperationType::Hmr => {

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

const CAP_SET: usize = 2048;

struct Affine {
    cst: Vec<bool>,
    set: Vec<Vec<u32>>,
    nextvar: u32,
    cond_stack: Vec<BitId>,

    b: Vec<Val>,
}

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

fn analyze_affine(
    ops: &[Op],
    num_q: usize,
    num_b: usize,
    input_qubits: &[QubitId],
) -> (Vec<Decision>, usize, usize) {

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

                if af.set[a] == af.set[b] {
                    if af.cst[a] == af.cst[b] {
                        decisions[i] = Decision::FoldEqualCtrls {
                            a: op.q_control1,
                            b: op.q_control2,
                            keep_ctrl: op.q_control1,
                        };
                        fold_eq += 1;

                        if af.cond_maybe_false(op) {
                            af.set[t] = af.fresh();
                            af.cst[t] = false;
                        } else {
                            let ns = xor_set(&af.set[t], &af.set[a]);
                            af.cst[t] ^= af.cst[a];
                            if ns.len() > CAP_SET {
                                af.set[t] = af.fresh();
                                af.cst[t] = false;
                            } else {
                                af.set[t] = ns;
                            }
                        }
                    } else {
                        decisions[i] = Decision::DropComplementCtrls {
                            a: op.q_control1,
                            b: op.q_control2,
                        };
                        drop_comp += 1;

                    }
                } else {

                    af.set[t] = af.fresh();
                    af.cst[t] = false;
                }
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

fn apply_decisions(ops: &[Op], decisions: &[Decision]) -> Vec<Op> {
    let mut out = Vec::with_capacity(ops.len());
    for (i, op) in ops.iter().enumerate() {
        match decisions[i] {
            Decision::Keep => out.push(*op),
            Decision::DropZeroCtrl { .. } => {  }
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
            Decision::DropComplementCtrls { .. } => {  }
            Decision::FoldEqualCtrls { keep_ctrl, .. } => {

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

fn apply_site_decisions(sites: &[OpSite], decisions: &[Decision]) -> Vec<OpSite> {
    let mut out = Vec::with_capacity(sites.len());
    for (i, site) in sites.iter().copied().enumerate() {
        match decisions[i] {
            Decision::Keep
            | Decision::FoldCx { .. }
            | Decision::FoldX { .. }
            | Decision::FoldEqualCtrls { .. } => out.push(site),
            Decision::DropZeroCtrl { .. } | Decision::DropComplementCtrls { .. } => {}
        }
    }
    out
}

fn filter_sites(sites: &[OpSite], kill: &[bool]) -> Vec<OpSite> {
    sites
        .iter()
        .copied()
        .enumerate()
        .filter_map(|(i, site)| (!kill[i]).then_some(site))
        .collect()
}

#[derive(Clone, Copy, Debug)]
struct PairKill {
    first: usize,
    second: usize,
}

#[derive(Clone, Copy)]
struct WEvent {
    idx: u32,
    src: u32,
    cond: u32,
    epoch: u32,
}

#[inline]
fn wev_written_between(ev: &[WEvent], lo: u32, hi: u32) -> bool {
    if hi <= lo + 1 {
        return false;
    }
    let start = ev.partition_point(|e| e.idx <= lo);
    start < ev.len() && ev[start].idx < hi
}

#[inline]
fn bit_written_between(ev: &[u32], lo: u32, hi: u32) -> bool {
    if hi <= lo + 1 {
        return false;
    }
    let start = ev.partition_point(|&x| x <= lo);
    start < ev.len() && ev[start] < hi
}

fn control_net_restored(
    ctrl: u64,
    p_idx: usize,
    cur_epoch: u64,
    cond_stack: &[u64],
    wev_q: &[Vec<WEvent>],
    wev_b: &[Vec<u32>],
) -> bool {
    let events = &wev_q[ctrl as usize];
    let p = p_idx as u32;
    let start = events.partition_point(|e| e.idx <= p);
    let suffix = &events[start..];
    if suffix.is_empty() {
        return true;
    }
    let mut stack: Vec<WEvent> = Vec::new();
    for &e in suffix {
        if e.src != u32::MAX {
            if let Some(&top) = stack.last() {

                let same = top.src == e.src
                    && top.cond == e.cond
                    && top.epoch == e.epoch
                    && e.epoch as u64 == cur_epoch;
                if same {
                    let src_ok =
                        !wev_written_between(&wev_q[e.src as usize], top.idx, e.idx);
                    let cond_ok = e.cond == u32::MAX
                        || !bit_written_between(&wev_b[e.cond as usize], top.idx, e.idx);
                    let stack_ok = cond_stack.iter().all(|&sb| {
                        sb == u64::MAX
                            || !bit_written_between(&wev_b[sb as usize], top.idx, e.idx)
                    });
                    if src_ok && cond_ok && stack_ok {
                        stack.pop();
                        continue;
                    }
                }
            }
        }
        stack.push(e);
    }
    stack.is_empty()
}

fn find_inverse_pairs(
    ops: &[Op],
    num_q: usize,
    num_b: usize,
    straddle: bool,
) -> (Vec<PairKill>, usize) {

    let mut wlast_q = vec![usize::MAX; num_q];
    let mut rlast_q = vec![usize::MAX; num_q];
    let mut wlast_b = vec![usize::MAX; num_b];

    for v in wlast_q.iter_mut() { *v = NEVER; }
    for v in rlast_q.iter_mut() { *v = NEVER; }
    for v in wlast_b.iter_mut() { *v = NEVER; }

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

    let mut cond_stack: Vec<u64> = Vec::new();
    let mut killed = vec![false; ops.len()];
    let mut pairs = Vec::new();

    let mut wev_q: Vec<Vec<WEvent>> = if straddle {
        vec![Vec::new(); num_q]
    } else {
        Vec::new()
    };
    let mut wev_b: Vec<Vec<u32>> = if straddle {
        vec![Vec::new(); num_b]
    } else {
        Vec::new()
    };

    let mut straddle_extra = 0usize;

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

                let mut cancelled = false;
                if let Some(p) = pending[t as usize] {
                    let same_gate = p.a == a && p.b == b && p.cb == cb;
                    let same_epoch = p.epoch == cond_epoch;

                    let ctrls_clean = !touched_after(wlast_q[a as usize], p.idx)
                        && !touched_after(wlast_q[b as usize], p.idx);

                    let ctrls_ok = if ctrls_clean {
                        true
                    } else if straddle {
                        control_net_restored(a, p.idx, cond_epoch, &cond_stack, &wev_q, &wev_b)
                            && control_net_restored(
                                b, p.idx, cond_epoch, &cond_stack, &wev_q, &wev_b,
                            )
                    } else {
                        false
                    };
                    let tgt_clean = !touched_after(wlast_q[t as usize], p.idx)
                        && !touched_after(rlast_q[t as usize], p.idx);
                    let cond_clean = cb == u64::MAX
                        || !touched_after(wlast_b[cb as usize], p.idx);

                    let stack_clean = same_epoch
                        && cond_stack
                            .iter()
                            .all(|&sb| sb == u64::MAX || !touched_after(wlast_b[sb as usize], p.idx));
                    if same_gate && same_epoch && ctrls_ok && tgt_clean && cond_clean && stack_clean {
                        killed[p.idx] = true;
                        killed[i] = true;
                        pairs.push(PairKill { first: p.idx, second: i });
                        pending[t as usize] = None;
                        cancelled = true;
                        if !ctrls_clean {
                            straddle_extra += 1;
                        }
                    }
                }

                if !cancelled {

                    rlast_q[a as usize] = i;
                    rlast_q[b as usize] = i;
                    wlast_q[t as usize] = i;
                    if cb != u64::MAX {

                    }
                    if straddle {

                        wev_q[t as usize].push(WEvent {
                            idx: i as u32,
                            src: u32::MAX,
                            cond: u32::MAX,
                            epoch: cond_epoch as u32,
                        });
                    }
                    pending[t as usize] = Some(Pending {
                        idx: i,
                        a,
                        b,
                        cb,
                        epoch: cond_epoch,
                    });
                } else {

                }
            }
            OperationType::CX => {
                rlast_q[op.q_control1.0 as usize] = i;
                wlast_q[op.q_target.0 as usize] = i;
                pending[op.q_target.0 as usize] = None;
                if straddle {

                    wev_q[op.q_target.0 as usize].push(WEvent {
                        idx: i as u32,
                        src: op.q_control1.0 as u32,
                        cond: op.c_condition.0 as u32,
                        epoch: cond_epoch as u32,
                    });
                }
            }
            OperationType::X => {
                wlast_q[op.q_target.0 as usize] = i;
                pending[op.q_target.0 as usize] = None;
                if straddle {

                    wev_q[op.q_target.0 as usize].push(WEvent {
                        idx: i as u32,
                        src: u32::MAX,
                        cond: u32::MAX,
                        epoch: cond_epoch as u32,
                    });
                }
            }
            OperationType::Swap => {
                let x = op.q_control1.0 as usize;
                let y = op.q_target.0 as usize;
                rlast_q[x] = i; rlast_q[y] = i;
                wlast_q[x] = i; wlast_q[y] = i;
                pending[x] = None;
                pending[y] = None;
                if straddle {
                    wev_q[x].push(WEvent { idx: i as u32, src: u32::MAX, cond: u32::MAX, epoch: cond_epoch as u32 });
                    wev_q[y].push(WEvent { idx: i as u32, src: u32::MAX, cond: u32::MAX, epoch: cond_epoch as u32 });
                }
            }
            OperationType::R => {
                wlast_q[op.q_target.0 as usize] = i;
                pending[op.q_target.0 as usize] = None;
                if straddle {
                    wev_q[op.q_target.0 as usize].push(WEvent { idx: i as u32, src: u32::MAX, cond: u32::MAX, epoch: cond_epoch as u32 });
                }
            }
            OperationType::Hmr => {
                wlast_q[op.q_target.0 as usize] = i;
                if op.c_target.0 != u64::MAX { wlast_b[op.c_target.0 as usize] = i; }
                pending[op.q_target.0 as usize] = None;
                if straddle {
                    wev_q[op.q_target.0 as usize].push(WEvent { idx: i as u32, src: u32::MAX, cond: u32::MAX, epoch: cond_epoch as u32 });
                    if op.c_target.0 != u64::MAX { wev_b[op.c_target.0 as usize].push(i as u32); }
                }
            }
            OperationType::CCZ => {

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
                if straddle && op.c_target.0 != u64::MAX {
                    wev_b[op.c_target.0 as usize].push(i as u32);
                }
            }
            OperationType::Neg
            | OperationType::Register
            | OperationType::AppendToRegister
            | OperationType::DebugPrint => {}
        }
    }

    (pairs, straddle_extra)
}

/// W018 / W044: straddle-aware CCZ self-inverse cancellation.
///
/// CCZ is diagonal and fully symmetric in its three qubits. Two CCZ on the same
/// unordered triple {a,b,c} compose to identity **iff**, between them, all three
/// qubits are NET-RESTORED (their per-branch computational-basis values at the 2nd
/// CCZ equal those at the 1st) and the condition context is unchanged. Under that
/// premise CCZ.U.CCZ = U exactly -- an identity in value AND phase -- so removing the
/// pair is bit-exact by construction (no phase census needed; contrast the M-60
/// never-fire census, which had no identity and broke on the phase channel).
///
/// Net-restore is decided by the SAME sound analysis (`control_net_restored`) the CCX
/// straddle path uses, applied to each of the three qubits. This pass does NOT cancel
/// CCX -- it treats every CCX as an opaque write -- so it is a conservative lower
/// bound on the straddle-restorable CCZ pairs, but every cancellation it makes is
/// sound. Intended to run on the FINAL post-`apply_m60_dead_t10` stream, so it never
/// perturbs the dead_t10 absolute-index skip-set.
pub(crate) fn ccz_straddle_cancel(ops: Vec<Op>) -> Vec<Op> {
    let (num_q, num_b) = dims(&ops);
    const OPAQUE: u32 = u32::MAX;

    let mut wev_q: Vec<Vec<WEvent>> = vec![Vec::new(); num_q];
    let mut wev_b: Vec<Vec<u32>> = vec![Vec::new(); num_b];
    let mut cond_epoch: u64 = 0;
    let mut cond_stack: Vec<u64> = Vec::new();

    #[derive(Clone, Copy)]
    struct PendCcz {
        idx: usize,
        cb: u64,
        epoch: u64,
    }
    let mut pending: std::collections::HashMap<(u64, u64, u64), PendCcz> =
        std::collections::HashMap::new();
    let mut killed = vec![false; ops.len()];

    let mut total_ccz = 0usize; // real (3-distinct-qubit) CCZ seen
    let mut candidates = 0usize; // same-triple, same cond/epoch (pre net-restore)
    let mut cancelled = 0usize;

    let push_q = |wev_q: &mut Vec<Vec<WEvent>>, q: u64, ev: WEvent| {
        if (q as usize) < wev_q.len() {
            wev_q[q as usize].push(ev);
        }
    };
    let push_b = |wev_b: &mut Vec<Vec<u32>>, b: u64, i: u32| {
        if (b as usize) < wev_b.len() {
            wev_b[b as usize].push(i);
        }
    };

    for (i, op) in ops.iter().enumerate() {
        let iu = i as u32;
        let ep = cond_epoch as u32;
        match op.kind {
            OperationType::PushCondition => {
                cond_epoch += 1;
                cond_stack.push(op.c_condition.0);
            }
            OperationType::PopCondition => {
                cond_epoch += 1;
                cond_stack.pop();
            }
            OperationType::CCZ => {
                let mut tri = [op.q_control1.0, op.q_control2.0, op.q_target.0];
                tri.sort_unstable();
                if tri[2] != u64::MAX && tri[0] != tri[1] && tri[1] != tri[2] {
                    total_ccz += 1;
                    let key = (tri[0], tri[1], tri[2]);
                    let cb = op.c_condition.0;
                    let pend = pending.get(&key).copied();
                    let mut did_cancel = false;
                    if let Some(p) = pend {
                        if p.cb == cb && p.epoch == cond_epoch {
                            candidates += 1;
                            let lo = p.idx as u32;
                            let qs_restored = tri.iter().all(|&q| {
                                control_net_restored(
                                    q, p.idx, cond_epoch, &cond_stack, &wev_q, &wev_b,
                                )
                            });
                            let cond_ok = cb == u64::MAX
                                || !bit_written_between(&wev_b[cb as usize], lo, iu);
                            let stack_ok = cond_stack.iter().all(|&sb| {
                                sb == u64::MAX
                                    || !bit_written_between(&wev_b[sb as usize], lo, iu)
                            });
                            if qs_restored && cond_ok && stack_ok {
                                killed[p.idx] = true;
                                killed[i] = true;
                                did_cancel = true;
                                cancelled += 1;
                            }
                        }
                    }
                    if did_cancel {
                        pending.remove(&key);
                    } else {
                        pending.insert(
                            key,
                            PendCcz {
                                idx: i,
                                cb,
                                epoch: cond_epoch,
                            },
                        );
                    }
                }
                // CCZ is diagonal: writes nothing, records no write-event.
            }
            OperationType::CX => {
                push_q(
                    &mut wev_q,
                    op.q_target.0,
                    WEvent {
                        idx: iu,
                        src: op.q_control1.0 as u32,
                        cond: op.c_condition.0 as u32,
                        epoch: ep,
                    },
                );
            }
            OperationType::CCX
            | OperationType::X
            | OperationType::R => {
                push_q(
                    &mut wev_q,
                    op.q_target.0,
                    WEvent { idx: iu, src: OPAQUE, cond: OPAQUE, epoch: ep },
                );
            }
            OperationType::Swap => {
                push_q(
                    &mut wev_q,
                    op.q_control1.0,
                    WEvent { idx: iu, src: OPAQUE, cond: OPAQUE, epoch: ep },
                );
                push_q(
                    &mut wev_q,
                    op.q_target.0,
                    WEvent { idx: iu, src: OPAQUE, cond: OPAQUE, epoch: ep },
                );
            }
            OperationType::Hmr => {
                push_q(
                    &mut wev_q,
                    op.q_target.0,
                    WEvent { idx: iu, src: OPAQUE, cond: OPAQUE, epoch: ep },
                );
                push_b(&mut wev_b, op.c_target.0, iu);
            }
            OperationType::BitInvert
            | OperationType::BitStore0
            | OperationType::BitStore1 => {
                push_b(&mut wev_b, op.c_target.0, iu);
            }
            OperationType::CZ
            | OperationType::Z
            | OperationType::Neg
            | OperationType::Register
            | OperationType::AppendToRegister
            | OperationType::DebugPrint => {}
        }
    }

    let n_before = ops.len();
    let kept: Vec<Op> = ops
        .into_iter()
        .enumerate()
        .filter_map(|(i, op)| if killed[i] { None } else { Some(op) })
        .collect();
    eprintln!(
        "  [W018 CCZ straddle] total_ccz={} same_triple_candidates={} cancelled_pairs={} removed_ccz={} -> {} ops",
        total_ccz,
        candidates,
        cancelled,
        n_before - kept.len(),
        kept.len()
    );
    kept
}

/// DIRECT MEASUREMENT (corpus-independent): run the shipped, sound CCX self-inverse
/// matcher on the FINAL post-fanout / post-dead_t10 stream. The production constprop
/// pass runs BEFORE `single_ccx_fanout` and `apply_m60_dead_t10`, both of which rewrite
/// the stream afterward -- so any self-inverse CCX adjacencies those two passes create
/// have never been seen by a canceller. Every pair `find_inverse_pairs` returns is a
/// proven self-inverse (same controls/target, clean or net-restored between) -> removing
/// it is bit-exact. Gated OFF by default (`TLM_CCX_FINAL_CANCEL=1` to enable) so the
/// baseline op-stream is unchanged for differential comparison. `straddle=false` by
/// default = strict clean case only (definitely bit-exact); `TLM_CCX_FINAL_STRADDLE=1`
/// widens to net-restore (reuses the CCX straddle path).
pub(crate) fn ccx_final_cancel(ops: Vec<Op>) -> Vec<Op> {
    if std::env::var("TLM_CCX_FINAL_CANCEL").ok().as_deref() != Some("1") {
        return ops;
    }
    let (nq, nb) = dims(&ops);
    let straddle = std::env::var("TLM_CCX_FINAL_STRADDLE").ok().as_deref() == Some("1");
    let (pairs, straddle_extra) = find_inverse_pairs(&ops, nq, nb, straddle);
    let mut killed = vec![false; ops.len()];
    for p in &pairs {
        killed[p.first] = true;
        killed[p.second] = true;
    }
    let kept: Vec<Op> = ops
        .into_iter()
        .enumerate()
        .filter_map(|(i, o)| if killed[i] { None } else { Some(o) })
        .collect();
    eprintln!(
        "  [FINAL CCX cancel] straddle={} pairs={} straddle_extra={} removed_ccx={} -> {} ops",
        straddle,
        pairs.len(),
        straddle_extra,
        pairs.len() * 2,
        kept.len()
    );
    kept
}

pub fn run(ops: Vec<Op>, input_qubits: &[QubitId]) -> Vec<Op> {
    let (num_q, num_b) = dims(&ops);
    let nonces_verify = std::env::var("CONSTPROP_VERIFY")
        .ok()
        .and_then(|s| s.parse::<usize>().ok());

    let verify_new_only = std::env::var("CONSTPROP_VERIFY_NEW_ONLY").ok().as_deref() == Some("1");

    let straddle = std::env::var("TLM_CONSTPROP_STRADDLE").ok().as_deref() == Some("1");

    let mut cur_sites = crate::point_add::take_op_site_trace_for_constprop(ops.len());
    let mut cur = ops;
    let mut iter = 0usize;
    let mut tot_dropped = 0usize;
    let mut tot_folded_cx = 0usize;
    let mut tot_folded_x = 0usize;
    let mut tot_pairs = 0usize;
    let mut tot_aff_drop = 0usize;
    let mut tot_aff_fold = 0usize;
    let mut tot_straddle_extra = 0usize;
    let affine_disabled = std::env::var("CONSTPROP_AFFINE_DISABLE").ok().as_deref() == Some("1");
    let max_iters = std::env::var("CONSTPROP_MAX_ITERS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(16);

    loop {
        iter += 1;

        let (mut decisions, stats) = analyze(&cur, num_q, num_b, input_qubits);

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
        if let Some(sites) = cur_sites.as_mut() {
            *sites = apply_site_decisions(sites, &decisions);
        }
        cur = apply_decisions(&cur, &decisions);

        let (nq2, nb2) = dims(&cur);
        let (pairs, straddle_extra) = find_inverse_pairs(&cur, nq2, nb2, straddle);
        tot_straddle_extra += straddle_extra;
        if straddle && straddle_extra > 0 {
            eprintln!(
                "CONSTPROP_STRADDLE iter={} extra_pairs={} (extra toffoli removed = {})",
                iter, straddle_extra, 2 * straddle_extra
            );
        }

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
            if let Some(sites) = cur_sites.as_mut() {
                *sites = filter_sites(sites, &kill);
            }
            cur = out;
        }

        let (mut aff_drop, mut aff_fold) = (0usize, 0usize);
        if !affine_disabled {
            let (nq3, nb3) = dims(&cur);
            let (mut adec, fold_eq, drop_comp) =
                analyze_affine(&cur, nq3, nb3, input_qubits);

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

            for d in &adec {
                match d {
                    Decision::DropComplementCtrls { .. } => aff_drop += 1,
                    Decision::FoldEqualCtrls { .. } => aff_fold += 1,
                    _ => {}
                }
            }
            if aff_drop + aff_fold > 0 {
                if let Some(sites) = cur_sites.as_mut() {
                    *sites = apply_site_decisions(sites, &adec);
                }
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
    if straddle {
        eprintln!(
            "CONSTPROP_STRADDLE TOTAL straddle_extra_pairs={} (of inverse_pairs={})",
            tot_straddle_extra, tot_pairs
        );
    }

    if let Some(sites) = cur_sites {
        crate::point_add::set_op_site_trace_from_constprop(sites);
    }

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

    let mut flagged: Vec<(usize, Vec<(QubitId, u64)>)> = Vec::new();
    for (i, d) in decisions.iter().enumerate() {
        match *d {
            Decision::Keep => {}
            Decision::DropZeroCtrl { ctrl } => flagged.push((i, vec![(ctrl, 0)])),
            Decision::FoldCx { one_ctrl, .. } => flagged.push((i, vec![(one_ctrl, 1)])),
            Decision::FoldX { c1, c2 } => flagged.push((i, vec![(c1, 1), (c2, 1)])),

            Decision::DropComplementCtrls { .. } | Decision::FoldEqualCtrls { .. } => {}
        }
    }
    let mut ok = vec![true; ops.len()];
    if flagged.is_empty() {
        return ok;
    }

    let mut flag_pos = vec![u32::MAX; ops.len()];
    for (p, (i, _)) in flagged.iter().enumerate() {
        flag_pos[*i] = p as u32;
    }

    const NUM_TESTS: usize = 9024;
    const BATCH: usize = 64;

    for nonce in 0..nonces {

        let mut hasher = Shake256::default();
        hasher.update(b"quantum_ecc-fiat-shamir-v2");
        hasher.update(&(ops.len() as u64).to_le_bytes());

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

            step_and_check(&mut sim, ops, &flag_pos, &flagged, &mut ok, cond_mask);
        }
        let bad = ok.iter().filter(|b| !**b).count();
        eprintln!(
            "CONSTPROP_PROGRESS nonce={}/{} shots={} cumulative_failed_claims={}",
            nonce + 1, nonces, n, bad
        );
    }
    let _ = QubitOrBit::Bit;
    ok
}

fn step_and_check<R: sha3::digest::XofReader>(
    sim: &mut crate::sim::Simulator<R>,
    ops: &[Op],
    flag_pos: &[u32],
    flagged: &[(usize, Vec<(QubitId, u64)>)],
    ok: &mut [bool],
    cond_mask: u64,
) {

    let mut condition_stack: Vec<u64> = Vec::new();
    let mut current_base_condition = u64::MAX;

    for (idx, op) in ops.iter().enumerate() {

        let fp = flag_pos[idx];
        if fp != u32::MAX {
            let p = fp as usize;
            for &(qid, expected) in &flagged[p].1 {
                let live = sim.qubit(qid) & cond_mask;
                let claim_ok = if expected == 0 {

                    live == 0
                } else {

                    live == cond_mask
                };
                if !claim_ok {
                    ok[idx] = false;
                }
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

#[cfg(test)]
mod affine_transfer_tests {
    use super::*;

    fn gate(kind: OperationType, c1: u64, c2: u64, target: u64) -> Op {
        let mut op = Op::empty();
        op.kind = kind;
        op.q_control1 = QubitId(c1);
        op.q_control2 = QubitId(c2);
        op.q_target = QubitId(target);
        op
    }

    #[test]
    fn closes_equal_and_complementary_chains_conservatively() {
        let equal_chain = vec![
            gate(OperationType::CX, 0, u64::MAX, 1),
            gate(OperationType::CCX, 0, 1, 2),
            gate(OperationType::CCX, 2, 0, 3),
        ];
        let (decisions, fold_eq, drop_comp) =
            analyze_affine(&equal_chain, 4, 0, &[QubitId(0)]);
        assert_eq!((fold_eq, drop_comp), (2, 0));
        assert!(matches!(decisions[1], Decision::FoldEqualCtrls { .. }));
        assert!(matches!(decisions[2], Decision::FoldEqualCtrls { .. }));

        let complement_chain = vec![
            gate(OperationType::CX, 0, u64::MAX, 1),
            gate(OperationType::X, u64::MAX, u64::MAX, 1),
            gate(OperationType::CX, 0, u64::MAX, 2),
            gate(OperationType::CCX, 0, 1, 2),
            gate(OperationType::CCX, 2, 0, 3),
        ];
        let (decisions, fold_eq, drop_comp) =
            analyze_affine(&complement_chain, 4, 0, &[QubitId(0)]);
        assert_eq!((fold_eq, drop_comp), (1, 1));
        assert!(matches!(decisions[3], Decision::DropComplementCtrls { .. }));
        assert!(matches!(decisions[4], Decision::FoldEqualCtrls { .. }));

        let mut conditional = gate(OperationType::CCX, 0, 1, 2);
        conditional.c_condition = BitId(0);
        let conditional_chain = vec![
            gate(OperationType::CX, 0, u64::MAX, 1),
            conditional,
            gate(OperationType::CCX, 2, 0, 3),
        ];
        let (decisions, fold_eq, drop_comp) =
            analyze_affine(&conditional_chain, 4, 1, &[QubitId(0)]);
        assert_eq!((fold_eq, drop_comp), (1, 0));
        assert!(matches!(decisions[1], Decision::FoldEqualCtrls { .. }));
        assert!(matches!(decisions[2], Decision::Keep));
    }
}

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

    let mut endpoint: Vec<u32> = vec![u32::MAX; ops.len()];
    let mut is_first_at: Vec<bool> = vec![false; ops.len()];
    for (p, pk) in pairs.iter().enumerate() {
        endpoint[pk.first] = p as u32;
        is_first_at[pk.first] = true;
        endpoint[pk.second] = p as u32;
        is_first_at[pk.second] = false;
    }

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

        let pp = endpoint[idx];
        if pp != u32::MAX {
            let p = pp as usize;

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

                snap_tgt[p] = tgt ^ contrib;
                snap_seen[p] = true;
            } else if snap_seen[p] {
                if contrib != snap_contrib[p] || tgt != snap_tgt[p] {
                    bad_pair[p] = true;
                }
            } else {

                bad_pair[p] = true;
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
    let _ = pairs;
}

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
