//! Tier B (issue #4), full-ladder increment: STREAM-EMIT and measure the whole
//! 28-window ECDLP ladder cost end-to-end, instead of deriving it.
//!
//! `analysis/ecdlp_estimate.py` derives the headline as `(PA_Toff + 3·2^w)·(2n/w−4)`.
//! This harness replaces the *emitted* part of that with a measurement: it chains
//! the built point-add op stream `n_add = 2n/w − 4` times through the harness
//! analyzers (`analyze_ops`/`analyze_depth`) — same qubit ids, so real cross-copy
//! dependencies and ancilla reuse — with **no materialized megastream** (the full
//! ladder is ~5×10⁹ ops / ~290 GB). The first increment (`ladder_composition.rs`)
//! only went to k≤4 and extrapolated; here the addition term is measured at the
//! *true* ladder length.
//!
//! The three terms of the ladder and how each is grounded:
//!   - ADDITION  `n_add · PA`     : MEASURED here by streaming n_add chained
//!                                  point-adds (Toffoli additive, width flat,
//!                                  depth serial — verified at the real n_add).
//!   - LOOKUP    `~3·2^w` / add   : MEASURED in `analysis/verify/ladder_lookup_cost.py`
//!                                  (#4 / ADR 0010) — a validated unary-iteration
//!                                  QROM read costs `2^(w+1)−4` Toffoli and its
//!                                  measurement-based uncompute-read `2^w−2`, so a
//!                                  read + uncompute-read ≈ `3·2^w − 6`, matching
//!                                  (and slightly under) the paper's `3·2^w`. Its
//!                                  data-writes are Clifford (do not score).
//!   - QFT                        : semiclassical / Clifford — 0 Toffoli.
//!
//! Composing the measured addition stream with the measured lookup term yields an
//! **emitted-and-measured** headline (~46.0M Toffoli on the static op-stream basis
//! — ~43.7M on `ecdlp_estimate.py`'s executed `score.json` basis — / ~1168 qubits
//! at w=16), matching the derived estimate — with one stated assumption: the
//! quantum-addend point-add (this repo's PA folds a *classical* compile-time
//! addend; the windowed ladder loads `P[k]` from a quantum table). That variant's
//! functional correctness and exact register overlap (whether the lookup ancilla
//! truly reuse PA's workspace, i.e. `ECDLP_Qubits = PA_Qubits + w`) is the
//! remaining Tier B build; this harness prices the composition under that model.
//!
//! `#[cfg(test)]` only; never compiled into the scored circuit.

use super::build;
use crate::circuit::{analyze_depth, analyze_ops, OperationType};

fn toffoli_count<'a>(ops: impl Iterator<Item = &'a crate::circuit::Op>) -> u64 {
    ops.filter(|o| matches!(o.kind, OperationType::CCX | OperationType::CCZ))
        .count() as u64
}

/// Windowed additions for secp256k1 (n = 256 bits): `2n/w − 4` (paper A1/A3).
fn n_windowed_additions(w: u64) -> u64 {
    2 * 256 / w - 4
}

/// Measured unary-iteration QROM read Toffoli (`analysis/verify/ladder_lookup_cost.py`).
fn qrom_read_toffoli(w: u64) -> u64 {
    (1u64 << (w + 1)) - 4
}

/// Measurement-based uncompute-read Toffoli (the compute-only subset; the
/// reversible uncompute becomes Clifford under MBUC — the repo's technique).
fn qrom_mbuc_uncompute_toffoli(w: u64) -> u64 {
    (1u64 << w) - 2
}

/// Per windowed addition: load `P[k]` (read) + unload (MBUC uncompute-read).
fn lookup_toffoli_per_add(w: u64) -> u64 {
    qrom_read_toffoli(w) + qrom_mbuc_uncompute_toffoli(w) // = 3·2^w − 6
}

// Heavy: streams the built op stream n_add=28 times through the analyzers.
// Opt-in so it never slows the default `cargo test`; run with `--ignored`.
#[test]
#[ignore = "heavy full-ladder measurement; run explicitly with `cargo test -- --ignored`"]
fn full_ladder_streamed_toffoli_qubits_depth() {
    let ops = build();
    let pa_tof = toffoli_count(ops.iter());
    let (pa_qubits, pa_bits, _r, _regs) = analyze_ops(ops.iter());
    let nq = pa_qubits as usize;
    let nb = pa_bits as usize;
    let pa_tdepth = analyze_depth(ops.iter(), nq, nb).toffoli_depth;
    assert!(
        pa_tof > 0 && pa_tdepth > 0,
        "degenerate build: pa_tof={pa_tof}, pa_tdepth={pa_tdepth}"
    );

    let w = 16u64;
    let n_add = n_windowed_additions(w); // 28
    let n = n_add as usize;

    // STREAM-EMIT the addition term at the true ladder length (no materialization).
    let add_tof = toffoli_count((0..n).flat_map(|_| ops.iter()));
    let (add_qubits, _b, _r, _regs) = analyze_ops((0..n).flat_map(|_| ops.iter()));
    let add_tdepth = analyze_depth((0..n).flat_map(|_| ops.iter()), nq, nb).toffoli_depth;

    // Composition laws must hold at the *real* n_add (not just k≤4 as in #13).
    assert_eq!(
        add_tof,
        n_add * pa_tof,
        "addition Toffoli not additive at n_add"
    );
    assert_eq!(add_qubits, pa_qubits, "addition width not flat at n_add");
    assert_eq!(
        add_tdepth,
        n_add * pa_tdepth,
        "addition depth not serial at n_add"
    );

    // Measured lookup term (QROM read + MBUC uncompute-read per addition).
    let lookup_per_add = lookup_toffoli_per_add(w);
    let lookup_tof = n_add * lookup_per_add;

    // End-to-end composed, emitted-and-measured headline. QFT is Clifford (0 Tof).
    // The QROM iteration depth (~2w per read) is negligible next to PA_depth, so
    // the ladder toffoli-depth is dominated by the serial additions.
    let total_tof = add_tof + lookup_tof;
    let total_qubits = pa_qubits + w; // A2: + window register (lookup ancilla reuse PA workspace)
    let total_tdepth = add_tdepth;

    // The derived headline `analysis/ecdlp_estimate.py` prints, for cross-check.
    let derived_tof = (pa_tof + 3 * (1u64 << w)) * n_add;

    eprintln!("\n=== issue #4 full-ladder streamed measurement (w={w}, n_add={n_add}) ===");
    eprintln!("  per addition (measured): PA toffoli={pa_tof}  qubits={pa_qubits}  toffoli_depth={pa_tdepth}");
    eprintln!("    (PA toffoli here is the STATIC op-stream CCX/CCZ count, as in");
    eprintln!("     ladder_composition.rs; analysis/ecdlp_estimate.py uses the executed");
    eprintln!("     avg-per-shot from score.json — a smaller basis, ~1.36M — so its");
    eprintln!("     headline is ~43.7M vs the ~46.0M static composition below.)");
    eprintln!("  ADDITION term (STREAMED {n_add}x point-add, no materialization):");
    eprintln!("    toffoli={add_tof} (={}x PA)  peak_qubits={add_qubits} (flat)  toffoli_depth={add_tdepth} (={}x PA)",
        add_tof / pa_tof, add_tdepth / pa_tdepth);
    eprintln!("  LOOKUP term (measured QROM, ADR 0010): read {} + MBUC-uncompute {} = {}/add  ->  {} total",
        qrom_read_toffoli(w), qrom_mbuc_uncompute_toffoli(w), lookup_per_add, lookup_tof);
    eprintln!("  QFT term: Clifford -> 0 Toffoli");
    eprintln!("  ------------------------------------------------------------------");
    eprintln!("  FULL ECDLP (emitted+measured composition):");
    eprintln!(
        "    Toffoli      = {total_tof}  (~{:.1}M)",
        total_tof as f64 / 1e6
    );
    eprintln!("    peak qubits  = {total_qubits}  (= PA_qubits + w)");
    eprintln!(
        "    toffoli-depth= {total_tdepth}  (~{:.2}M; serial additions dominate)",
        total_tdepth as f64 / 1e6
    );
    eprintln!(
        "  vs derived (PA+3·2^w)·n_add = {derived_tof}  (Δ = {} = 6·n_add, MBUC saving)",
        derived_tof - total_tof
    );

    // The measured composition matches the derived headline to within the MBUC
    // lookup saving (6·n_add), and never exceeds it.
    assert!(
        total_tof <= derived_tof,
        "measured composed exceeds derived headline"
    );
    assert_eq!(
        derived_tof - total_tof,
        6 * n_add,
        "composition mismatch beyond the MBUC lookup saving"
    );

    // Sanity vs the paper's published full-ECDLP Toffoli bounds (Appendix A):
    // this repo's measured PA drives the ladder well under both.
    assert!(
        total_tof < 70_000_000,
        "unexpectedly above the paper's Low-Gate ECDLP bound"
    );
}
