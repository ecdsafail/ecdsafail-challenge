//! Tier B (issue #4), first increment: MEASURE how point additions compose into
//! the ECDLP double-and-add ladder, instead of asserting it.
//!
//! The scored circuit is one optimized point addition. Shor's ladder runs the
//! same primitive many times against a running accumulator; because the addend
//! is a classical register (and a qubit-recycled QFT feeds scalar bits forward
//! as classical conditions), a controlled addition has the SAME gate count as an
//! uncontrolled one, and successive additions reuse the same accumulator +
//! ancilla. We measure that by chaining the built op stream `k` times (same
//! qubit ids ⇒ real cross-copy dependencies + ancilla reuse) through the harness
//! analyzers — no pipeline re-run, no materialized megastream.
//!
//! Confirms three previously-asserted claims (ADR 0007):
//!   - Toffoli is additive: total = k · PA_Toffoli.
//!   - peak width is FLAT in k (ancilla reused) ⇒ `ECDLP_Qubits = PA_Qubits + w`.
//!   - toffoli-depth is serial: composed ≈ k · PA_depth (accumulator serializes).
//!
//! `#[cfg(test)]` only; never in the scored circuit. (Does NOT cover the windowed
//! QROM lookup or the QFT — those need the quantum-addend build, deferred in #4.)

use super::build;
use crate::circuit::{analyze_depth, analyze_ops, OperationType};

fn toffoli_count<'a>(ops: impl Iterator<Item = &'a crate::circuit::Op>) -> u64 {
    ops.filter(|o| matches!(o.kind, OperationType::CCX | OperationType::CCZ))
        .count() as u64
}

// Heavy: builds the full point-add circuit and re-scans the op stream for k=2..4.
// Opt-in so it never slows the default `cargo test`; run with `--ignored`.
#[test]
#[ignore = "heavy measurement harness; run explicitly with `cargo test -- --ignored`"]
fn ladder_composition_is_additive_flat_width_serial_depth() {
    let ops = build();

    let pa_tof = toffoli_count(ops.iter());
    let (pa_qubits, pa_bits, _r, _regs) = analyze_ops(ops.iter());
    let nq = pa_qubits as usize;
    let nb = pa_bits as usize;
    let pa_tdepth = analyze_depth(ops.iter(), nq, nb).toffoli_depth;

    // A well-formed point-add has non-zero Toffoli count and depth; assert it so
    // the ratio reporting below can never divide by zero.
    assert!(
        pa_tof > 0 && pa_tdepth > 0,
        "degenerate build: pa_tof={pa_tof}, pa_tdepth={pa_tdepth}"
    );

    eprintln!("\n=== issue #4 measured ladder composition (chained k additions) ===");
    eprintln!(
        "  k=1 (per addition): toffoli={pa_tof}  peak_qubits={pa_qubits}  toffoli_depth={pa_tdepth}"
    );

    for k in [2u64, 3, 4] {
        let ku = k as usize;
        let tof = toffoli_count((0..ku).flat_map(|_| ops.iter()));
        let (qubits, _b, _r, _regs) = analyze_ops((0..ku).flat_map(|_| ops.iter()));
        let tdepth = analyze_depth((0..ku).flat_map(|_| ops.iter()), nq, nb).toffoli_depth;
        eprintln!(
            "  k={k}: toffoli={tof} (={}x)  peak_qubits={qubits} (Δ={})  toffoli_depth={tdepth} (={:.3}x)",
            tof / pa_tof,
            qubits as i64 - pa_qubits as i64,
            tdepth as f64 / pa_tdepth as f64,
        );

        // Toffoli additive.
        assert_eq!(tof, k * pa_tof, "k={k}: toffoli not additive");
        // Peak width flat: chaining reuses the same qubit ids (ancilla returned
        // to |0>), so the composed width equals a single addition's.
        assert_eq!(
            qubits, pa_qubits,
            "k={k}: peak width grew (ancilla not reused)"
        );
        // Depth serial: the accumulator written by copy i is read by copy i+1, so
        // the non-Clifford critical path is ~k times a single addition's.
        assert_eq!(tdepth, k * pa_tdepth, "k={k}: toffoli-depth not serial");
    }

    // Project to the paper's optimal window (w=16 ⇒ 28 windowed additions): the
    // measured composition says the addition term is exactly 28·PA on a flat
    // (PA_qubits)-wide register file, serialized to 28·PA_depth. The lookup and
    // QFT terms remain derived (deferred; see ADR 0007).
    let n_add = 28u64;
    eprintln!(
        "  => measured addition term @w=16: {}·PA = {} Toffoli, peak {} qubits, depth {}",
        n_add,
        n_add * pa_tof,
        pa_qubits,
        n_add * pa_tdepth,
    );
}
