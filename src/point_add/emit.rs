use super::*;

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
pub(crate) fn emit_inverse<F: FnOnce(&mut B)>(b: &mut B, f: F) {
    if b.count_only {
        let snap = b.count_snapshot();
        if b.fiat_hash.is_some() {
            let hash_before = b.fiat_hash.clone();
            b.count_only_capture_stack.push(Vec::new());
            f(b);
            let forward = b
                .count_only_capture_stack
                .pop()
                .expect("count-only inverse capture stack");
            b.restore_count_snapshot(snap);
            b.fiat_hash = hash_before;
            emit_inverse_ops_allowing_clean_resets(
                b,
                &forward,
                "count-only hashed emit_inverse",
            );
            return;
        }
        f(b);
        let delta = b.count_delta_since(snap);
        b.restore_count_snapshot(snap);
        add_inverse_count_delta(b, &delta);
        return;
    }
    let start = b.ops.len();
    f(b);
    let end = b.ops.len();
    // Extract the forward slice and drop it from the builder.
    let fwd: Vec<_> = b.ops[start..end].to_vec();
    b.ops.truncate(start);
    emit_inverse_ops_allowing_clean_resets(b, &fwd, "emit_inverse");
}

pub(crate) fn add_inverse_count_delta(b: &mut B, delta: &[usize; 18]) {
    for kind in [
        OperationType::X,
        OperationType::Z,
        OperationType::CX,
        OperationType::CZ,
        OperationType::CCX,
        OperationType::CCZ,
        OperationType::Swap,
    ] {
        b.add_counted_kind(kind, delta[kind as usize]);
    }
}

pub(crate) fn emit_inverse_ops_allowing_clean_resets(b: &mut B, fwd: &[Op], context: &'static str) {
    for op in fwd.iter().rev().copied() {
        match op.kind {
            OperationType::X
            | OperationType::Z
            | OperationType::CX
            | OperationType::CZ
            | OperationType::CCX
            | OperationType::CCZ
            | OperationType::Swap => b.push_op(op),
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
                "{context}: non-invertible op kind {:?} inside forward block",
                op.kind
            ),
        }
    }
}
