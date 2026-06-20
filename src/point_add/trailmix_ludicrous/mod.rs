//! Product-min secp256k1 EC point-add for `point_add`, built directly on
//! quantum_ecc's own `B` builder. Emitted as the default circuit by
//! [`super::build`].
//!
//! The circuit modules (`ec_add`/`arith`/`gcd`/`square`/`comparator`/`codec`/
//! `gidney`/`fused`/`mcx`) call `B` directly. Two things `B` does not itself
//! provide live here: the [`BExt`] trait adds the gates the modules need that
//! `B` lacks (`z`/`ccz`/`neg`/`cswap` and the bit-conditioned forms), and a
//! thread-local [`Sched`] holds the per-call replay cursors that drive the
//! product-min operating point. So the emitted op-stream is native quantum_ecc
//! `Op`s produced through `B`'s own allocation + recycling.

mod arith;
mod codec;
mod comparator;
mod constprop;
mod ec_add;
mod fused;
mod gcd;
mod gidney;
mod mcx;
pub mod schedule;
mod square;

pub use schedule::PAD;

use super::B;
use crate::circuit::{BitId, Op, OperationType, QubitId};
use std::cell::RefCell;
use std::collections::HashMap;

const N: usize = 256;

// ── gates B does not expose ──────────────────────────────────────────────
// `B` provides x/cx/ccx/cz/swap/hmr and the alloc/condition primitives. The
// circuit modules also need a phase Z, a doubly-controlled Z, the free Neg, a
// Fredkin swap, and bit-conditioned forms; this trait supplies them on `B`.
pub(super) trait BExt {
    fn loan_zero_qubit(&mut self, q: QubitId);
    fn reclaim_zero_qubit(&mut self, q: QubitId);
    fn z(&mut self, q: QubitId);
    fn ccz(&mut self, a: QubitId, b: QubitId, c: QubitId);
    fn neg(&mut self);
    fn cswap(&mut self, ctrl: QubitId, a: QubitId, b: QubitId);
    fn x_if_bit(&mut self, q: QubitId, c: BitId);
    fn z_if_bit(&mut self, q: QubitId, c: BitId);
    fn cz_if_bit(&mut self, a: QubitId, b: QubitId, c: BitId);
    /// Return a qubit to the free-list (emits an R reset; the caller must have
    /// already uncomputed it to |0>).
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

// ── per-call schedule replay ─────────────────────────────────────────────
// The product-min operating point is reached by replaying baked per-call
// choices (carry caps, vent counts, branch and fold selections). They are set
// once at the start of the build and read in sequence as the circuit is
// emitted; a thread-local holds the cursors so the circuit fns stay `&mut B`.
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

/// Read the next entry of a `(values, cursor)` slot, returning `exhausted` past
/// the end (a sentinel meaning "no schedule constraint": full headroom).
fn step<T: Copy>(slot: &mut (Vec<T>, usize), exhausted: T) -> T {
    let v = slot.0.get(slot.1).copied().unwrap_or(exhausted);
    slot.1 += 1;
    v
}

fn env_i32(name: &str, default: i32) -> i32 {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse::<i32>().ok())
        .unwrap_or(default)
}

fn env_usize(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(default)
}

fn sub_delta(v: usize, name: &str) -> usize {
    if v == usize::MAX {
        v
    } else {
        v.saturating_sub(env_usize(name, 0))
    }
}

fn range_adjust_usize(idx: usize, v: usize, prefix: &str) -> usize {
    if v == usize::MAX {
        return v;
    }
    let adj = env_i32(&format!("{prefix}_ADJUST"), 0);
    if adj == 0 {
        return v;
    }
    let after = env_usize(&format!("{prefix}_ADJUST_AFTER"), 0);
    let before = env_usize(&format!("{prefix}_ADJUST_BEFORE"), usize::MAX);
    if idx >= after && idx < before {
        (v as i32).saturating_add(adj).max(0) as usize
    } else {
        v
    }
}

fn range_adjust_i32(idx: usize, v: i32, prefix: &str) -> i32 {
    if v == i32::MAX || v < 0 {
        return v;
    }
    let adj = env_i32(&format!("{prefix}_ADJUST"), 0);
    if adj == 0 {
        return v;
    }
    let after = env_usize(&format!("{prefix}_ADJUST_AFTER"), 0);
    let before = env_usize(&format!("{prefix}_ADJUST_BEFORE"), usize::MAX);
    if idx >= after && idx < before {
        v.saturating_add(adj).max(0)
    } else {
        v
    }
}

fn next_gcd_k() -> usize { SCHED.with(|s| step(&mut s.borrow_mut().gcd_k, usize::MAX)) }
fn next_cout_k() -> usize {
    SCHED.with(|s| {
        let mut s = s.borrow_mut();
        let idx = s.cout_k.1;
        let v = sub_delta(step(&mut s.cout_k, usize::MAX), "TLM_COUT_K_DELTA");
        range_adjust_usize(idx, v, "TLM_COUT_K")
    })
}
fn next_fold() -> i32 {
    SCHED.with(|s| {
        let mut s = s.borrow_mut();
        let idx = s.fold.1;
        let mut v = step(&mut s.fold, i32::MAX);
        let delta = env_usize("TLM_FOLD_DELTA", 0) as i32;
        if v != i32::MAX && v >= 0 && delta != 0 {
            v = v.saturating_sub(delta);
        }
        range_adjust_i32(idx, v, "TLM_FOLD")
    })
}
fn next_gcd_branch() -> u8 { SCHED.with(|s| step(&mut s.borrow_mut().gcd_branch, 255)) }
fn next_cmp_k() -> usize { SCHED.with(|s| step(&mut s.borrow_mut().cmp_k, usize::MAX)) }
fn next_ffg() -> usize {
    SCHED.with(|s| {
        let mut s = s.borrow_mut();
        let idx = s.ffg.1;
        let v = sub_delta(step(&mut s.ffg, usize::MAX), "TLM_FFG_DELTA");
        range_adjust_usize(idx, v, "TLM_FFG")
    })
}
fn next_hyb_v() -> usize {
    SCHED.with(|s| {
        let mut s = s.borrow_mut();
        let idx = s.hyb_v.1;
        let v = sub_delta(step(&mut s.hyb_v, usize::MAX), "TLM_HYB_V_DELTA");
        range_adjust_usize(idx, v, "TLM_HYB_V")
    })
}
fn next_sqrow_k() -> usize { SCHED.with(|s| step(&mut s.borrow_mut().sqrow_k, usize::MAX)) }

/// Load the product-min jump schedule onto the thread-local cursors.
fn load_schedule() {
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

/// Swaps that route the value at qubit `src[i]` to qubit `dst[i]` (placing R.x bit
/// i into reg0 slot dst[i]). Qubits in `dst` not in `src` are |0> ancilla, so the
/// routing is value-preserving.
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

/// Build the product-min EC-add op-stream natively via `B`, with the 4
/// evaluator registers (reg0=R.x qubits, reg1=R.y qubits, reg2=Q.x bits,
/// reg3=Q.y bits) and the grinding tail nonce appended.
pub fn build_trailmix_ludicrous_ops() -> Vec<Op> {
    let mut circ = B::new();
    load_schedule();

    // Allocation order fixes the ids that become the IO registers.
    let x2 = circ.alloc_qubits(N); // reg0: P.x -> R.x
    let y2 = circ.alloc_qubits(N); // reg1: P.y -> R.y
    let ox = circ.alloc_bits(N); // reg2: Q.x (classical)
    let oy = circ.alloc_bits(N); // reg3: Q.y (classical)

    let x2_init = x2.clone();
    let mut x2m = x2;
    ec_add::ec_add(&mut circ, &mut x2m, &y2, &ox, &oy);

    // ── register declarations + result routing + tail nonce ──
    circ.declare_qubit_register(&x2_init);
    circ.declare_qubit_register(&y2);
    circ.declare_bit_register(&ox);
    circ.declare_bit_register(&oy);

    // Route R.x (scattered ids x2m) back into reg0 ids (x2_init).
    for (a, b) in route_swaps(&x2m, &x2_init) {
        circ.swap(a, b);
    }

    // Grinding tail nonce: 48 X;X identity pairs (96 X gates) on reg0[0]/reg0[1].
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

    if std::env::var_os("TRACE_PEAK").is_some()
        || std::env::var_os("TRACE_PHASE_ACTIVE").is_some()
    {
        eprintln!(
            "TLM_TRACE peak_qubits={} phase='{}' ops_idx={} total_ops={}",
            circ.peak_qubits,
            circ.peak_phase,
            circ.peak_ops_idx,
            circ.ops.len()
        );
    }
    if std::env::var_os("TRACE_PEAK_VERBOSE").is_some() {
        let pk = circ.peak_qubits;
        for (active, phase, op) in &circ.peak_log {
            if *active + 10 >= pk {
                eprintln!("TLM_TRACE near_peak active={} phase='{}' ops_idx={}", active, phase, op);
            }
        }
    }
    if std::env::var_os("TRACE_PHASE_ACTIVE").is_some() {
        let mut v: Vec<_> = circ.phase_active_max.iter().collect();
        v.sort_by_key(|(_, active)| std::cmp::Reverse(**active));
        let top_n = std::env::var("TRACE_PHASE_ACTIVE_TOP")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(20);
        for (phase, active) in v.into_iter().take(top_n) {
            eprintln!("TLM_TRACE active_by_phase active={} phase='{}'", active, phase);
        }
    }

    let ops = std::mem::take(&mut circ.ops);

    // Sound classical constant-propagation peephole: drop CCX with a provably
    // |0> quantum control (no-op but still scored) and fold CCX with a provably
    // |1> control to CX/X. reg0 (x2_init) and reg1 (y2) hold per-shot input data
    // -> seeded Unknown; every other qubit id is a |0> ancilla. Can be disabled
    // with CONSTPROP_DISABLE=1.
    if std::env::var("CONSTPROP_DISABLE").ok().as_deref() == Some("1") {
        return ops;
    }
    let mut input_qubits = x2_init.clone();
    input_qubits.extend_from_slice(&y2);
    constprop::run(ops, &input_qubits)
}
