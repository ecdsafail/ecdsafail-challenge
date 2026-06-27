//! The GCD swap-decision comparator for the product-min secp256k1 EC-add,
//! built on this crate's `B` builder.
//!
//! ## What the schedule needs
//! At GCD jump step `i` the swap decision compares the GCD cofactors `u`,`v`
//! over a narrow window, the per-step `GAP_J2[i]`. The comparator scans the
//! top `k` bits of the operands; it mis-decides the `u <-> v` swap iff the
//! highest differing bit of `u`,`v` sits below the window, i.e. it returns
//! "equal -> no swap" when the top-`k` MSBs of the two operands agree. On
//! uniform operands this happens with probability ~2^-k per call.
//!
//! ## Carry handling: chunked Cuccaro + Gidney held carries
//! The backend is `compare_geq_chunked_middle`. The bottom `[0, n-k)` bits run an
//! in-place Cuccaro `a >= b` MAJ chain (one live carry, uncomputed exactly by the
//! self-inverse UMA). The top `[n-k, n)` bits hold `k` Gidney carries that are
//! measure-erased (MBU) on the reverse, so only `k+1` carries are live when the
//! caller body runs. The held-carry count `k` is supplied per call from the
//! schedule (`next_cmp_k`): `k = 0` is pure in-place Cuccaro (peak-safe), `k >= n`
//! is full Gidney.

use super::{BExt, B};
use crate::circuit::QubitId;
use std::cell::Cell;

thread_local! {
    static CMP_CALL_INDEX: Cell<usize> = const { Cell::new(0) };
    static CMP_HOST_RESERVE_CALL_INDEX: Cell<Option<usize>> = const { Cell::new(None) };
}

fn env_cmp_call_i32(name: &str, call_index: usize) -> Option<i32> {
    std::env::var(name).ok().and_then(|value| {
        value
            .split(',')
            .filter_map(|item| item.trim().split_once(':'))
            .find_map(|(call, value)| {
                (call.parse::<usize>().ok()? == call_index)
                    .then(|| value.parse::<i32>().ok())
                    .flatten()
            })
    })
}

fn adjusted_active_k(call_index: usize, k: usize, max_k: usize) -> usize {
    if let Some(override_k) = env_cmp_call_i32("TLM_CMP_ACTIVE_K_CALL_OVERRIDES", call_index) {
        return override_k.clamp(1, max_k as i32) as usize;
    }
    let global_delta = env_i32("TLM_CMP_ACTIVE_K_DELTA");
    let call_delta = env_cmp_call_i32("TLM_CMP_ACTIVE_K_CALL_DELTAS", call_index).unwrap_or(0);
    (k as i32)
        .saturating_add(global_delta)
        .saturating_add(call_delta)
        .clamp(1, max_k as i32) as usize
}

fn env_i32(name: &str) -> i32 {
    std::env::var(name)
        .ok()
        .and_then(|value| value.parse::<i32>().ok())
        .unwrap_or(0)
}

fn adjusted_carry_k(call_index: usize, k: usize, max_k: usize) -> usize {
    if let Some(override_k) = env_cmp_call_i32("TLM_CMP_CARRY_K_CALL_OVERRIDES", call_index) {
        return override_k.clamp(0, max_k as i32) as usize;
    }
    let global_delta = env_i32("TLM_CMP_CARRY_K_DELTA");
    let call_delta = env_cmp_call_i32("TLM_CMP_CARRY_K_CALL_DELTAS", call_index).unwrap_or(0);
    (k as i32)
        .saturating_add(global_delta)
        .saturating_add(call_delta)
        .clamp(0, max_k as i32) as usize
}

fn adjusted_host_reserve(call_index: usize, default: usize) -> usize {
    if let Some(override_k) = env_cmp_call_i32("TLM_CMP_HOST_RESERVE_CALL_OVERRIDES", call_index) {
        return override_k.max(0) as usize;
    }
    let call_delta = env_cmp_call_i32("TLM_CMP_HOST_RESERVE_CALL_DELTAS", call_index).unwrap_or(0);
    (default as i32).saturating_add(call_delta).max(0) as usize
}

/// Direct-final-carry form of the chunked `a >= b` comparator. `k` is the held
/// carry count (clamped to `n`). `body(circ, carry)` sees
/// `carry = (a >= b)` and must restore it before returning; `a`/`b` and all
/// carries are then cleaned.
fn compare_geq_chunked_middle_direct<F: FnOnce(&mut B, &QubitId)>(
    circ: &mut B,
    a: &[QubitId],
    b: &[QubitId],
    body: F,
    k: usize,
) {
    let n = a.len();
    assert_eq!(
        b.len(),
        n,
        "compare_geq_chunked_middle_direct: a,b equal width"
    );
    assert!(
        n > 0,
        "compare_geq_chunked_middle_direct: nonempty operands"
    );
    // Force `host_reserve` extra low bits onto the in-place Cuccaro bottom (which
    // reuses the single live `c` lane) instead of holding them as Gidney carries.
    // This is value-EXACT (the predicate is the full n-bit compare either way); it
    // just shrinks the simultaneously-held carry set by `host_reserve`, relieving the
    // GCD swap-decision comparator wall by that many qubits at the peak.
    let global_host_reserve = std::env::var("TLM_CMP_HOST_RESERVE")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(0);
    let host_reserve = CMP_HOST_RESERVE_CALL_INDEX.with(|index| {
        index
            .get()
            .map_or(global_host_reserve, |call| adjusted_host_reserve(call, global_host_reserve))
    });
    let k = super::target_qubit_headroom(circ)
        .map_or(k, |headroom| {
            k.min(headroom.saturating_sub(1 + host_reserve))
        })
        .min(n);
    let split = n - k; // bottom [0, split) in-place; top [split, n) held.
    let mut cy: Vec<Option<QubitId>> = (0..=n).map(|_| None).collect();
    let c = circ.alloc_qubit();
    circ.x(c); // c_0 = 1
               // Forward bottom: in-place Cuccaro MAJ (only `c` live).
    for i in 0..split {
        circ.x(b[i]);
        circ.cx(c, b[i]);
        circ.cx(c, a[i]);
        circ.ccx(a[i], b[i], c); // c = c_{i+1}
    }
    cy[split] = Some(c);
    // Forward top: Gidney held carries.
    for i in split..n {
        let next = circ.alloc_qubit();
        {
            let ci = cy[i].as_ref().unwrap();
            circ.x(b[i]);
            circ.cx(*ci, b[i]);
            circ.cx(*ci, a[i]);
            circ.ccx(a[i], b[i], next);
            circ.cx(*ci, next); // next = c_{i+1}
        }
        cy[i + 1] = Some(next);
    }
    body(circ, cy[n].as_ref().unwrap());
    // Reverse top: measure-erase the held carries.
    for i in (split..n).rev() {
        let next = cy[i + 1].take().unwrap();
        circ.cx(*cy[i].as_ref().unwrap(), next); // c_{i+1} -> ta&tb
        let bit = circ.alloc_bit();
        circ.hmr(next, bit);
        circ.zero_and_free(next);
        circ.cz_if_bit(a[i], b[i], bit);
        circ.cx(*cy[i].as_ref().unwrap(), a[i]);
        circ.cx(*cy[i].as_ref().unwrap(), b[i]);
        circ.x(b[i]);
    }
    // Reverse bottom: in-place UMA.
    let c = cy[split].take().unwrap();
    for i in (0..split).rev() {
        circ.ccx(a[i], b[i], c);
        circ.cx(c, a[i]);
        circ.cx(c, b[i]);
        circ.x(b[i]);
    }
    circ.x(c); // c_0 = 1 -> 0
    circ.zero_and_free(c);
}

/// Compatibility form for callers that need a persistent predicate qubit.
/// New middle-only users should use the direct-final-carry form above.
pub fn compare_geq_chunked_middle<F: FnOnce(&mut B, &QubitId)>(
    circ: &mut B,
    a: &[QubitId],
    b: &[QubitId],
    flag: &QubitId,
    body: F,
    k: usize,
) {
    assert_eq!(
        b.len(),
        a.len(),
        "compare_geq_chunked_middle: a,b equal width"
    );
    if a.is_empty() {
        circ.x(*flag);
        body(circ, flag);
        circ.x(*flag);
        return;
    }
    compare_geq_chunked_middle_direct(
        circ,
        a,
        b,
        |c, carry| {
            c.cx(*carry, *flag);
            body(c, flag);
            c.cx(*carry, *flag);
        },
        k,
    );
}

/// Controlled form: `target ^= ctrl AND (u_top < v_top)` on the top `k` MSBs.
/// Used where the swap decision is itself gated.
pub fn controlled_swap_decision_lt_truncated(
    circ: &mut B,
    ctrl: &QubitId,
    u: &[QubitId],
    v: &[QubitId],
    k: usize,
    target: &QubitId,
) {
    assert!(
        k > 0 && k <= u.len() && k <= v.len(),
        "k must fit in both operands"
    );
    let call = CMP_CALL_INDEX.with(|index| {
        let call = index.get();
        index.set(call + 1);
        call
    });
    let effective_k = adjusted_active_k(call, k, u.len().min(v.len()));
    let u_top: Vec<QubitId> = u[u.len() - effective_k..].to_vec();
    let v_top: Vec<QubitId> = v[v.len() - effective_k..].to_vec();
    // The callback consumes the final carry directly, eliminating the old
    // predicate-copy lane. Spend exactly that lane on one additional held carry;
    // compare_geq_chunked_middle_direct's target-Q clamp still bounds the same peak.
    let base_ck = super::next_cmp_k();
    let ck = adjusted_carry_k(call, base_ck.saturating_add(1), u_top.len());
    let trace_call = std::env::var_os("TRACE_TLM_CMP").is_some().then_some(call);
    let timeline_start = circ.active_timeline.len();
    let entry_active = circ.active_qubits;
    let entry_ops = circ.current_ops_len();
    CMP_HOST_RESERVE_CALL_INDEX.with(|index| index.set(Some(call)));
    compare_geq_chunked_middle_direct(
        circ,
        &u_top,
        &v_top,
        |c, carry| {
            c.x(*carry);
            c.ccx(*ctrl, *carry, *target);
            c.x(*carry);
        },
        ck,
    );
    CMP_HOST_RESERVE_CALL_INDEX.with(|index| index.set(None));
    if let Some(call) = trace_call {
        let local_peak = circ.active_timeline[timeline_start..]
            .iter()
            .map(|(_, active)| *active)
            .max()
            .unwrap_or(circ.active_qubits);
        eprintln!(
            "TLM_CMP call={call} phase={} n={} k={} base_ck={} ck={} entry_active={} local_peak={} global_peak={} entry_ops={} exit_ops={}",
            circ.phase,
            u.len(),
            effective_k,
            base_ck,
            ck,
            entry_active,
            local_peak,
            circ.peak_qubits,
            entry_ops,
            circ.current_ops_len(),
        );
    }
}

/// Measurement-vented `a + b + cin` carry chain with a middle callback. Computes
/// the chain `cy[i+1] = carry of (a + ~b + ~cin)` bit-by-bit, then at the top bit
/// hands `(ta = a_top ^ c, tb = ~b_top ^ c, c_{n-1})` to `body` -- the final carry
/// `cy[n] = (ta AND tb) XOR c_{n-1}` is NOT built, so `body` deposits its phase via
/// a bare Z/CZ/CCZ on those three (no value flip), riding through the reverse
/// measure-uncompute. Reverse vents each internal carry by `hmr` + `cz_if_bit`.
/// `a`,`b`,`cin` restored. Equal-width `a`,`b` (the chunked-erase caller).
///
/// `carry-out(a + b + cin) = NOT carry-out(a + ~b + ~cin)`, so a caller wanting
/// to test `[a + b + cin >= 2^n]` reads the complement of the built predicate.
pub fn compare_geq_cin_middle<F: FnOnce(&mut B, &QubitId, &QubitId, &QubitId)>(
    circ: &mut B,
    a: &[QubitId],
    b: &[QubitId],
    cin: &QubitId,
    body: F,
) {
    let n = a.len();
    assert_eq!(b.len(), n, "compare_geq_cin_middle: a,b equal width");
    assert!(n >= 1, "needs >= 1 bit");
    if std::env::var("TLM_COMPARE_CIN_CUCCARO").ok().as_deref() == Some("1") {
        let c = circ.alloc_qubit();
        circ.x(c);
        circ.cx(*cin, c); // c_0 = ~cin
        for i in 0..n - 1 {
            circ.x(b[i]);
            circ.cx(c, b[i]);
            circ.cx(c, a[i]);
            circ.ccx(a[i], b[i], c); // c = c_{i+1}
        }
        {
            let i = n - 1;
            circ.x(b[i]);
            circ.cx(c, b[i]);
            circ.cx(c, a[i]);
            body(circ, &a[i], &b[i], &c);
            circ.cx(c, a[i]);
            circ.cx(c, b[i]);
            circ.x(b[i]);
        }
        for i in (0..n - 1).rev() {
            circ.ccx(a[i], b[i], c);
            circ.cx(c, a[i]);
            circ.cx(c, b[i]);
            circ.x(b[i]);
        }
        circ.cx(*cin, c);
        circ.x(c);
        circ.zero_and_free(c);
        return;
    }
    let borrow_cin_c0 = std::env::var("TLM_COMPARE_BORROW_CIN_C0").ok().as_deref() == Some("1");
    let mut cy: Vec<Option<QubitId>> = Vec::with_capacity(n);
    if borrow_cin_c0 {
        circ.x(*cin); // cy[0] = ~cin, restored before return.
        cy.push(Some(*cin));
    } else {
        let c0 = circ.alloc_qubit();
        circ.x(c0);
        circ.cx(*cin, c0); // cy[0] = 1 ^ cin = ~cin (carry-in of a + ~b + ~cin)
        cy.push(Some(c0));
    }
    for i in 0..n - 1 {
        let next = circ.alloc_qubit();
        let ci = cy[i].as_ref().unwrap();
        circ.x(b[i]);
        circ.cx(*ci, b[i]);
        circ.cx(*ci, a[i]);
        circ.ccx(a[i], b[i], next);
        circ.cx(*ci, next);
        cy.push(Some(next));
    }
    // Top bit: fold only, hand (ta, tb, c_{n-1}) to body.
    {
        let i = n - 1;
        let ci = cy[i].as_ref().unwrap();
        circ.x(b[i]);
        circ.cx(*ci, b[i]);
        circ.cx(*ci, a[i]);
        body(circ, &a[i], &b[i], ci);
        circ.cx(*ci, a[i]);
        circ.cx(*ci, b[i]);
        circ.x(b[i]);
    }
    // Reverse: vent cy[1..n-1] via hmr, restore a/b.
    for i in (0..n - 1).rev() {
        let next = cy[i + 1].take().unwrap();
        let ci_raw = cy[i].as_ref().unwrap();
        circ.cx(*ci_raw, next); // next = ta_i & tb_i
        let bit = circ.alloc_bit();
        circ.hmr(next, bit);
        circ.zero_and_free(next);
        circ.cz_if_bit(a[i], b[i], bit);
        circ.cx(*cy[i].as_ref().unwrap(), a[i]);
        circ.cx(*cy[i].as_ref().unwrap(), b[i]);
        circ.x(b[i]);
    }
    let c0 = cy[0].take().unwrap();
    if borrow_cin_c0 {
        debug_assert_eq!(c0, *cin);
        circ.x(*cin);
    } else {
        circ.cx(*cin, c0); // ~cin -> 1
        circ.x(c0); // 1 -> 0
        circ.zero_and_free(c0);
    }
}

/// Vented uncompute of a GCD swap-decision flag that holds `ctrl AND (v_top <
/// u_top)` (the forward decision). HMR the flag to |0> (0 Toffoli), then
/// under the HMR `push_condition` recompute the predicate as a deferred Z, using
/// the direct-final-carry comparator on `(v_top, u_top)`. `v`,`u` restored. The forward
/// computes the flag normally (a value); only this reverse clear vents.
///
/// `[v >= u] = carryout(v + ~u + 1)` over the top-`k` window; the predicate is
/// `ctrl AND (v < u) = ctrl AND NOT[v>=u]`. Deposit the phase directly through
/// the final carry via `X(carry); CZ(ctrl, carry); X(carry)`.
pub fn swap_decision_uncompute_vented(
    circ: &mut B,
    ctrl: &QubitId,
    v: &[QubitId],
    u: &[QubitId],
    k: usize,
    flag: &QubitId,
) {
    assert!(
        k > 0 && k <= v.len() && k <= u.len(),
        "k must fit in both operands"
    );
    let v_top: Vec<QubitId> = v[v.len() - k..].to_vec();
    let u_top: Vec<QubitId> = u[u.len() - k..].to_vec();
    // Match the forward decision: removing the predicate-copy lane funds one
    // additional held carry without changing peak liveness.
    let ck = super::next_cmp_k().saturating_add(1);
    let bit = circ.alloc_bit();
    circ.hmr(*flag, bit);
    circ.push_condition(bit);
    compare_geq_chunked_middle_direct(
        circ,
        &v_top,
        &u_top,
        |c, carry| {
            // Deposit Z^(ctrl AND NOT carry) = Z^(ctrl AND [v < u]), gated by the HMR
            // condition (push_condition). Same phase as the cin (ta,tb,c_prev) form.
            c.x(*carry);
            c.cz(*ctrl, *carry);
            c.x(*carry);
        },
        ck,
    );
    circ.pop_condition();
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::circuit::OperationType;
    use crate::sim::Simulator;
    use sha3::{
        digest::{ExtendableOutput, Update},
        Shake256,
    };

    fn alloc_case(circ: &mut B, n: usize) -> (Vec<QubitId>, Vec<QubitId>, QubitId, QubitId) {
        let a = (0..n).map(|_| circ.alloc_qubit()).collect();
        let b = (0..n).map(|_| circ.alloc_qubit()).collect();
        let ctrl = circ.alloc_qubit();
        let target = circ.alloc_qubit();
        (a, b, ctrl, target)
    }

    fn xor_value(circ: &mut B, qs: &[QubitId], value: usize) {
        for (i, &q) in qs.iter().enumerate() {
            if (value >> i) & 1 != 0 {
                circ.x(q);
            }
        }
    }

    fn simulate(circ: &B) -> (Vec<u64>, u64) {
        let mut shake = Shake256::default();
        shake.update(b"comparator-direct-final-carry-test");
        let mut xof = shake.finalize_xof();
        let mut sim = Simulator::new(circ.next_qubit as usize, circ.next_bit as usize, &mut xof);
        sim.apply_iter(circ.ops.iter());
        (sim.qubits, sim.phase)
    }

    fn read_uniform(qs: &[QubitId], qubits: &[u64]) -> usize {
        qs.iter().enumerate().fold(0usize, |value, (i, q)| {
            let lane = qubits[q.0 as usize];
            assert!(
                lane == 0 || lane == u64::MAX,
                "nonuniform data lane q{}",
                q.0
            );
            value | (usize::from(lane == u64::MAX) << i)
        })
    }

    fn toffoli_count(circ: &B) -> usize {
        circ.ops
            .iter()
            .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
            .count()
    }

    #[test]
    fn direct_final_carry_is_exhaustive_for_small_widths() {
        for n in 1..=4 {
            let limit = 1usize << n;
            for held in 0..=n {
                for a_value in 0..limit {
                    for b_value in 0..limit {
                        for ctrl_value in 0..=1usize {
                            for target_value in 0..=1usize {
                                let mut circ = B::new();
                                let (a, b, ctrl, target) = alloc_case(&mut circ, n);
                                xor_value(&mut circ, &a, a_value);
                                xor_value(&mut circ, &b, b_value);
                                if ctrl_value != 0 {
                                    circ.x(ctrl);
                                }
                                if target_value != 0 {
                                    circ.x(target);
                                }
                                compare_geq_chunked_middle_direct(
                                    &mut circ,
                                    &a,
                                    &b,
                                    |c, carry| {
                                        c.x(*carry);
                                        c.ccx(ctrl, *carry, target);
                                        c.cz(ctrl, *carry);
                                        c.x(*carry);
                                    },
                                    held,
                                );

                                assert_eq!(circ.active_qubits as usize, 2 * n + 2);
                                let (qubits, phase) = simulate(&circ);
                                let predicate = ctrl_value != 0 && a_value < b_value;
                                assert_eq!(read_uniform(&a, &qubits), a_value);
                                assert_eq!(read_uniform(&b, &qubits), b_value);
                                assert_eq!(
                                    qubits[ctrl.0 as usize],
                                    if ctrl_value != 0 { u64::MAX } else { 0 }
                                );
                                assert_eq!(
                                    qubits[target.0 as usize],
                                    if (target_value != 0) ^ predicate {
                                        u64::MAX
                                    } else {
                                        0
                                    },
                                );
                                assert_eq!(phase, if predicate { u64::MAX } else { 0 });
                                assert!(qubits[2 * n + 2..].iter().all(|&q| q == 0));
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn freed_predicate_lane_funds_one_held_carry() {
        for n in 1..=8 {
            for held in 0..n {
                let mut legacy = B::new();
                let (a, b, ctrl, target) = alloc_case(&mut legacy, n);
                let flag = legacy.alloc_qubit();
                compare_geq_chunked_middle(
                    &mut legacy,
                    &a,
                    &b,
                    &flag,
                    |c, flag| {
                        c.x(*flag);
                        c.ccx(ctrl, *flag, target);
                        c.x(*flag);
                    },
                    held,
                );
                legacy.zero_and_free(flag);

                let mut direct = B::new();
                let (a, b, ctrl, target) = alloc_case(&mut direct, n);
                compare_geq_chunked_middle_direct(
                    &mut direct,
                    &a,
                    &b,
                    |c, carry| {
                        c.x(*carry);
                        c.ccx(ctrl, *carry, target);
                        c.x(*carry);
                    },
                    held + 1,
                );

                assert_eq!(direct.peak_qubits, legacy.peak_qubits, "n={n} held={held}");
                assert_eq!(
                    toffoli_count(&direct) + 1,
                    toffoli_count(&legacy),
                    "n={n} held={held}",
                );
            }
        }
    }
}
