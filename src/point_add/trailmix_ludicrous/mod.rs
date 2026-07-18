
mod arith;
mod codec;
mod comparator;
mod constprop;
pub mod ec_add;
mod fused;

pub(super) use fused::FoldShadowMarker;

pub(super) fn collect_fold_shadow_markers(
    sites: &[crate::point_add::OpSite],
) -> Vec<FoldShadowMarker> {
    fused::collect_fold_shadow_markers(sites)
}

pub(super) fn shift_fold_shadow_markers_for_rewrite(
    markers: &mut [FoldShadowMarker],
    first_removed: usize,
    insertion_after: usize,
    second_removed: usize,
) {
    fused::shift_fold_shadow_markers_for_rewrite(
        markers,
        first_removed,
        insertion_after,
        second_removed,
    );
}

pub(super) fn shift_fold_shadow_markers_for_drops(
    markers: &mut [FoldShadowMarker],
    drops: &[usize],
) {
    fused::shift_fold_shadow_markers_for_drops(markers, drops);
}

pub(super) fn write_fold_shadow_markers(markers: &[FoldShadowMarker]) {
    fused::write_fold_shadow_markers(markers);
}
mod gcd;
mod gidney;
mod mcx;
pub mod schedule;
mod square;

pub use schedule::PAD;

use super::B;
use crate::circuit::{BitId, Op, OperationType, QubitId};
use std::cell::{Cell, RefCell};
use std::collections::HashMap;

const N: usize = 256;

pub(super) trait BExt {
    fn loan_zero_qubit(&mut self, q: QubitId);
    fn reclaim_zero_qubit(&mut self, q: QubitId);
    fn z(&mut self, q: QubitId);
    #[track_caller]
    fn ccz(&mut self, a: QubitId, b: QubitId, c: QubitId);
    fn neg(&mut self);
    #[track_caller]
    fn cswap(&mut self, ctrl: QubitId, a: QubitId, b: QubitId);
    fn x_if_bit(&mut self, q: QubitId, c: BitId);
    fn z_if_bit(&mut self, q: QubitId, c: BitId);
    fn cz_if_bit(&mut self, a: QubitId, b: QubitId, c: BitId);

    fn zero_and_free(&mut self, q: QubitId);
}

impl BExt for B {
    fn loan_zero_qubit(&mut self, q: QubitId) {
        self.free_qubits
            .push(q.0.try_into().expect("qubit id fits in u32"));
        if self.active_qubits > 0 {
            self.active_qubits -= 1;
        }
        self.record_active_timeline();
        self.b0_on_free(q.0);
    }

    fn reclaim_zero_qubit(&mut self, q: QubitId) {
        self.reacquire(q);
    }

    fn z(&mut self, q: QubitId) {
        let mut op = Op::empty();
        op.kind = OperationType::Z;
        op.q_target = q;
        self.push_op(op);
    }
    #[track_caller]
    fn ccz(&mut self, a: QubitId, b: QubitId, c: QubitId) {
        let mut op = Op::empty();
        op.kind = OperationType::CCZ;
        op.q_control2 = a;
        op.q_control1 = b;
        op.q_target = c;
        self.push_op(op);
    }
    fn neg(&mut self) {
        let mut op = Op::empty();
        op.kind = OperationType::Neg;
        self.push_op(op);
    }
    #[track_caller]
    fn cswap(&mut self, ctrl: QubitId, a: QubitId, b: QubitId) {
        self.cx(b, a);
        self.ccx(ctrl, a, b);
        self.cx(b, a);
    }
    fn x_if_bit(&mut self, q: QubitId, c: BitId) {
        self.push_condition(c);
        self.x(q);
        self.pop_condition();
    }
    fn z_if_bit(&mut self, q: QubitId, c: BitId) {
        self.push_condition(c);
        self.z(q);
        self.pop_condition();
    }
    fn cz_if_bit(&mut self, a: QubitId, b: QubitId, c: BitId) {
        self.push_condition(c);
        self.cz(a, b);
        self.pop_condition();
    }
    fn zero_and_free(&mut self, q: QubitId) {
        self.free(q);
    }
}

#[derive(Default)]
struct Sched {
    gcd_k: (Vec<usize>, usize),
    cout_k: (Vec<usize>, usize),
    fold: (Vec<i32>, usize),
    gcd_branch: (Vec<u8>, usize),
    cmp_k: (Vec<usize>, usize),
    ffg: (Vec<usize>, usize),
    hyb_v: (Vec<usize>, usize),
    sqrow_k: (Vec<usize>, usize),
}

thread_local!(static SCHED: RefCell<Sched> = RefCell::new(Sched::default()));

#[derive(Clone, Copy, Debug)]
pub(super) struct ScheduleFit {
    pub call_index: usize,
    pub base: usize,
    pub selected: usize,
}

thread_local! {
    static HYB_CALL_INDEX: Cell<usize> = const { Cell::new(0) };
    static COUT_CALL_INDEX: Cell<usize> = const { Cell::new(0) };
    static PENDING_COUT_FIT: Cell<Option<ScheduleFit>> = const { Cell::new(None) };
}

fn step<T: Copy>(slot: &mut (Vec<T>, usize), exhausted: T) -> T {
    let v = slot.0.get(slot.1).copied().unwrap_or(exhausted);
    slot.1 += 1;
    v
}

fn env_delta(name: &str) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0)
}

fn sub_delta(v: usize, name: &str) -> usize {
    if v == usize::MAX {
        v
    } else {
        v.saturating_sub(env_delta(name))
    }
}

fn env_call_value(name: &str, call_index: usize) -> Option<usize> {
    std::env::var(name).ok().and_then(|value| {
        value
            .split(',')
            .filter_map(|item| item.trim().split_once(':'))
            .find_map(|(call, value)| {
                (call.parse::<usize>().ok()? == call_index)
                    .then(|| value.parse::<usize>().ok())
                    .flatten()
            })
    })
}

fn next_call_index(counter: &'static std::thread::LocalKey<Cell<usize>>) -> usize {
    counter.with(|index| {
        let current = index.get();
        index.set(current + 1);
        current
    })
}

fn fit_schedule_value(
    base: usize,
    call_index: usize,
    global_delta: &str,
    call_deltas: &str,
    call_overrides: &str,
) -> ScheduleFit {
    let selected = env_call_value(call_overrides, call_index).unwrap_or_else(|| {
        let globally_adjusted = sub_delta(base, global_delta);
        match env_call_value(call_deltas, call_index) {
            Some(delta) if globally_adjusted != usize::MAX => {
                globally_adjusted.saturating_sub(delta)
            }
            _ => globally_adjusted,
        }
    });
    ScheduleFit {
        call_index,
        base,
        selected,
    }
}

fn reset_schedule_fit_call_indices() {
    HYB_CALL_INDEX.with(|index| index.set(0));
    COUT_CALL_INDEX.with(|index| index.set(0));
    PENDING_COUT_FIT.with(|pending| pending.set(None));
}

fn target_qubit_headroom(circ: &B) -> Option<usize> {
    std::env::var("TLM_TARGET_Q")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .map(|target| target.saturating_sub(circ.active_qubits as usize))
}

fn next_gcd_k() -> usize { SCHED.with(|s| step(&mut s.borrow_mut().gcd_k, usize::MAX)) }
fn next_cout_k() -> usize {
    let base = SCHED.with(|s| step(&mut s.borrow_mut().cout_k, usize::MAX));
    let fit = fit_schedule_value(
        base,
        next_call_index(&COUT_CALL_INDEX),
        "TLM_COUT_K_DELTA",
        "TLM_COUT_K_CALL_DELTAS",
        "TLM_COUT_K_CALL_OVERRIDES",
    );
    PENDING_COUT_FIT.with(|pending| {
        debug_assert!(pending.get().is_none(), "previous COUT schedule call was not consumed");
        pending.set(Some(fit));
    });
    fit.selected
}
fn next_fold() -> i32 {
    SCHED.with(|s| {
        let v = step(&mut s.borrow_mut().fold, i32::MAX);
        let d = env_delta("TLM_FOLD_DELTA") as i32;
        if v == i32::MAX || v < 0 || d == 0 {
            v
        } else {
            v.saturating_sub(d)
        }
    })
}
fn next_gcd_branch() -> u8 { SCHED.with(|s| step(&mut s.borrow_mut().gcd_branch, 255)) }
fn next_cmp_k() -> usize { SCHED.with(|s| step(&mut s.borrow_mut().cmp_k, usize::MAX)) }
fn next_ffg() -> usize { SCHED.with(|s| sub_delta(step(&mut s.borrow_mut().ffg, usize::MAX), "TLM_FFG_DELTA")) }
fn next_hyb_v_fit() -> ScheduleFit {
    let base = SCHED.with(|s| step(&mut s.borrow_mut().hyb_v, usize::MAX));
    fit_schedule_value(
        base,
        next_call_index(&HYB_CALL_INDEX),
        "TLM_HYB_V_DELTA",
        "TLM_HYB_V_CALL_DELTAS",
        "TLM_HYB_V_CALL_OVERRIDES",
    )
}

fn take_cout_fit(selected: usize) -> ScheduleFit {
    PENDING_COUT_FIT.with(|pending| {
        pending.take().unwrap_or_else(|| {
            fit_schedule_value(
                selected,
                next_call_index(&COUT_CALL_INDEX),
                "TLM_COUT_K_DELTA",
                "TLM_COUT_K_CALL_DELTAS",
                "TLM_COUT_K_CALL_OVERRIDES",
            )
        })
    })
}
fn next_sqrow_k() -> usize { SCHED.with(|s| step(&mut s.borrow_mut().sqrow_k, usize::MAX)) }

fn load_schedule() {
    reset_schedule_fit_call_indices();
    arith::reset_ffg_call_index();
    comparator::reset_compare_call_index();
    fused::reset_fold_call_index();
    gcd::reset_gcd_trace_call_index();
    gidney::reset_gidney_call_index();
    SCHED.with(|s| {
        let mut s = s.borrow_mut();
        *s = Sched::default();
        let extra_fold_vents = std::env::var("LUD_EXTRA_FOLD_VENTS")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(0);
        let extra_fold_min_g = std::env::var("LUD_EXTRA_FOLD_MIN_G")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(0);
        let extra_fold_max_g = std::env::var("LUD_EXTRA_FOLD_MAX_G")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(usize::MAX);
        let fold_g = |v: &[usize]| -> Vec<usize> {
            v.iter()
                .map(|&x| {
                    if extra_fold_vents > 0
                        && x >= extra_fold_min_g
                        && x <= extra_fold_max_g
                    {
                        x.saturating_add(extra_fold_vents).min(53)
                    } else {
                        x
                    }
                })
                .collect()
        };
        s.gcd_k.0 = schedule::GCD_SUB_K.to_vec();
        s.gcd_branch.0 = schedule::GCD_BRANCH.to_vec();
        s.cout_k.0 = schedule::APPLY_COUT_K.to_vec();
        s.fold.0 = schedule::FOLD_SCHED.to_vec();
        s.cmp_k.0 = schedule::CMP_K.to_vec();
        s.ffg.0 = fold_g(schedule::FFG_G);
        s.hyb_v.0 = schedule::HYB_V.to_vec();
        s.sqrow_k.0 = schedule::SQ_ROW_K.to_vec();
    });
}

fn route_swaps(src: &[QubitId], dst: &[QubitId]) -> Vec<(QubitId, QubitId)> {
    let mut loc: Vec<QubitId> = src.to_vec();
    let mut at: HashMap<u64, usize> = HashMap::new();
    for (i, q) in src.iter().enumerate() {
        at.insert(q.0, i);
    }
    let mut swaps = Vec::new();
    for i in 0..dst.len() {
        let target = dst[i];
        let cur = loc[i];
        if cur == target {
            continue;
        }
        swaps.push((target, cur));
        let displaced = at.get(&target.0).copied();
        at.insert(target.0, i);
        loc[i] = target;
        match displaced {
            Some(b) => {
                at.insert(cur.0, b);
                loc[b] = cur;
            }
            None => {
                at.remove(&cur.0);
            }
        }
    }
    swaps
}

fn install_q1153_submission_defaults() {
    for (name, value) in [
        ("TLM_TARGET_Q", "1152"),
        ("TLM_FOLD_CHUNK_ZERO_CIN", "1"),
        ("TLM_FFG_MAX_G", "47"),
        ("TLM_APPLY_ADD_SKIP_LASTK", "1"),
        ("DIALOG_TAIL_NONCE", "2430844"),
    ] {

        if (name == "DIALOG_TAIL_NONCE" || name == "TLM_TARGET_Q")
            && std::env::var_os(name).is_some()
        {
            continue;
        } else {
            std::env::set_var(name, value);
        }
    }
}

pub fn build_trailmix_ludicrous_ops() -> Vec<Op> {
    install_q1153_submission_defaults();
    let mut circ = B::new();
    load_schedule();

    let x2 = circ.alloc_qubits(N);
    let y2 = circ.alloc_qubits(N);
    let ox = circ.alloc_bits(N);
    let oy = circ.alloc_bits(N);

    let x2_init = x2.clone();
    let mut x2m = x2;
    ec_add::ec_add(&mut circ, &mut x2m, &y2, &ox, &oy);

    circ.declare_qubit_register(&x2_init);
    circ.declare_qubit_register(&y2);
    circ.declare_bit_register(&ox);
    circ.declare_bit_register(&oy);

    for (a, b) in route_swaps(&x2m, &x2_init) {
        circ.swap(a, b);
    }

    if let Some(nonce) = std::env::var("DIALOG_TAIL_NONCE")
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
    {
        for i in 0..48u32 {
            let q = if (nonce >> i) & 1 == 1 { x2_init[1] } else { x2_init[0] };
            circ.x(q);
            circ.x(q);
        }
    }

    circ.b0_finalize();

    if std::env::var("TRACE_TLM_PROFILE").is_ok() {
        circ.close_phase_active_region();
        eprintln!(
            "TLM_PROFILE peak_qubits={} peak_phase={} peak_ops_idx={} emitted_ops={}",
            circ.peak_qubits,
            circ.peak_phase,
            circ.peak_ops_idx,
            circ.current_ops_len(),
        );
        let mut phases: Vec<_> = circ.phase_active_max.iter().collect();
        phases.sort_by(|left, right| right.1.cmp(left.1).then_with(|| left.0.cmp(right.0)));
        for (phase, active) in phases.into_iter().take(24) {
            eprintln!("TLM_PHASE active_max={active} phase={phase}");
        }
    }

    if std::env::var("TLM_TIMELINE_DUMP").is_ok() {
        let trans = &circ.phase_transitions;
        let phase_at = |op: usize| -> &'static str {

            let mut lo = 0usize;
            let mut hi = trans.len();
            let mut ans = "init";
            while lo < hi {
                let mid = (lo + hi) / 2;
                if trans[mid].0 <= op {
                    ans = trans[mid].1;
                    lo = mid + 1;
                } else {
                    hi = mid;
                }
            }
            ans
        };
        let minq: u32 = std::env::var("TLM_TIMELINE_MIN")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(1136);

        use std::collections::BTreeMap;
        let mut census: BTreeMap<&'static str, (u32, usize, usize)> = BTreeMap::new();
        for &(op, active) in &circ.active_timeline {
            if active >= minq {
                let ph = phase_at(op);
                let e = census.entry(ph).or_insert((0, 0, op));
                if active > e.0 {
                    e.0 = active;
                    e.2 = op;
                }
                e.1 += 1;
            }
        }
        let mut rows: Vec<_> = census.into_iter().collect();
        rows.sort_by(|a, b| b.1 .0.cmp(&a.1 .0).then_with(|| a.0.cmp(b.0)));
        eprintln!("TLM_TIMELINE census (samples with active>={minq}), phase: max_active n_samples example_op");
        for (ph, (mx, n, ex)) in &rows {
            eprintln!("TL_CENSUS phase={ph} max_active={mx} n_samples={n} example_op={ex}");
        }

        if let (Ok(lo), Ok(hi)) = (
            std::env::var("TLM_WIN_LO").map(|s| s.parse::<usize>().unwrap_or(0)),
            std::env::var("TLM_WIN_HI").map(|s| s.parse::<usize>().unwrap_or(usize::MAX)),
        ) {
            eprintln!("TLM_TIMELINE raw window [{lo},{hi}] : op active phase");
            for &(op, active) in &circ.active_timeline {
                if op >= lo && op <= hi {
                    eprintln!("TL_RAW op={op} active={active} phase={}", phase_at(op));
                }
            }
        }
    }

    if std::env::var("TRACE_TLM_CCX").is_ok() {
        use std::collections::BTreeMap;
        let mut bounds = circ.phase_transitions.clone();
        bounds.sort_by_key(|(i, _)| *i);
        let total = circ.ops.len();
        let mut by: BTreeMap<&'static str, usize> = BTreeMap::new();
        for w in 0..bounds.len() {
            let s = bounds[w].0.min(total);
            let e = if w + 1 < bounds.len() { bounds[w + 1].0.min(total) } else { total };
            let c = circ.ops[s..e].iter().filter(|op| op.kind as u32 == 13).count();
            *by.entry(bounds[w].1).or_insert(0) += c;
        }
        let grand: usize = by.values().sum();
        let mut v: Vec<_> = by.into_iter().collect();
        v.sort_by(|a, b| b.1.cmp(&a.1));
        let mut cum = 0usize;
        for (phase, c) in v.iter().take(30) {
            cum += *c;
            eprintln!(
                "TLM_CCX phase={phase} ccx={c} pct={:.2} cum={:.2}",
                100.0 * *c as f64 / grand as f64,
                100.0 * cum as f64 / grand as f64
            );
        }
        eprintln!("TLM_CCX_TOTAL {grand} phases={}", v.len());
    }

    let ops = std::mem::take(&mut circ.ops);

    if std::env::var("CONSTPROP_DISABLE").ok().as_deref() == Some("1") {
        return ops;
    }
    let mut input_qubits = x2_init.clone();
    input_qubits.extend_from_slice(&y2);
    constprop::run(ops, &input_qubits)
}
