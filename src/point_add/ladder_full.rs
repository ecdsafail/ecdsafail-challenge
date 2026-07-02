//! Tier B (issue #4), full-ladder increment: STREAM-EMIT and count the whole
//! 28-window ECDLP ladder — additions AND table lookups — end-to-end, instead of
//! deriving it.
//!
//! `analysis/ecdlp_estimate.py` derives the headline as `(PA_Toff + 3·2^w)·(2n/w−4)`.
//! This harness replaces the *emitted* part with a measurement: per window it
//! chains an emitted **QROM-read op stream** (the lookup) with the built
//! **point-add op stream** (the addition) and another read (the unload), and
//! counts the whole thing through `analyze_ops` / `analyze_depth` — same qubit
//! ids across windows (real cross-copy hazards + ancilla reuse), with **no
//! materialized megastream** (a full ladder is ~5×10⁹ ops / ~290 GB).
//!
//! The lookup is emitted as an optimized **unary-iteration** QROM selector
//! (Gidney 2018 §III.C), identical in structure to
//! `analysis/verify/ladder_lookup_cost.py` (ADR 0010) — so this Rust emission
//! independently reproduces that Python measurement: a read is `2^(w+1)−4`
//! Toffoli. Only the Toffoli-bearing **selector** is emitted; a lookup's
//! data-writes are `2^w·d` **Clifford** `CX` (~290 GB of ops at w=16, d=512) that
//! do not affect the Toffoli/score metric — they are the reason the full stream
//! cannot be materialized, and are counted analytically as Clifford (0 Toffoli).
//!
//! The three terms and how each is grounded here:
//!   - ADDITION  `n_add · PA`  : EMITTED (chained point-add op stream at the true
//!                               n_add; Toffoli additive, width flat, depth serial).
//!   - LOOKUP    per addition  : EMITTED (unary-iteration QROM read + unload). A
//!                               reversible unload is another read; with the
//!                               repo's measurement-based uncompute the unload is
//!                               Clifford, so per-addition lookup Toffoli is
//!                               `read + read/2 = 3·2^w − 6 ≤ 3·2^w`.
//!   - QFT                     : semiclassical / Clifford — 0 Toffoli.
//!
//! Reported (static op-stream basis, w=16): full-ladder Toffoli is emitted+counted
//! (~47.8M reversible / ~46.0M with MBUC unload), matching the derived
//! `(PA+3·2^w)·n_add` to within the MBUC saving. `ecdlp_estimate.py`'s headline
//! uses the *executed* avg-per-shot PA (~1.36M → ~43.7M), a smaller basis; this
//! harness cross-validates that closed form rather than replacing it.
//!
//! **Not measured here** (the one remaining Tier B build): the **peak qubits** and
//! **exact depth** of a *functionally composed* ladder need the **quantum-addend
//! point-add** — this repo's PA folds a *classical* compile-time addend, whereas
//! the ladder loads `P[k]` from a *quantum* table that the addition then consumes.
//! Only that variant fixes the real register overlap behind `ECDLP_Qubits =
//! PA_Qubits + w` (A2) and the read→add data dependency. Emitting the lookup on
//! *disjoint* ids (as here) would OVER-count width and UNDER-count the serial
//! depth, so width is reported per A2 and depth as the measured add-dominated
//! critical path, both flagged. Issue #5's mid-ladder ∞/`dx=0` residual lands in
//! that same quantum-addend testbed and is out of scope for this cost harness.
//!
//! `#[cfg(test)]` only; never compiled into the scored circuit.

use super::build;
use crate::circuit::{analyze_depth, analyze_ops, Op, OperationType, QubitId};

fn toffoli_count<'a>(ops: impl Iterator<Item = &'a Op>) -> u64 {
    ops.filter(|o| matches!(o.kind, OperationType::CCX | OperationType::CCZ))
        .count() as u64
}

fn ccx(c1: u64, c2: u64, t: u64) -> Op {
    let mut o = Op::empty();
    o.kind = OperationType::CCX;
    o.q_control1 = QubitId(c1);
    o.q_control2 = QubitId(c2);
    o.q_target = QubitId(t);
    o
}
fn cx(c1: u64, t: u64) -> Op {
    let mut o = Op::empty();
    o.kind = OperationType::CX;
    o.q_control1 = QubitId(c1);
    o.q_target = QubitId(t);
    o
}
fn x(t: u64) -> Op {
    let mut o = Op::empty();
    o.kind = OperationType::X;
    o.q_target = QubitId(t);
    o
}

/// Emit the Toffoli-bearing **selector** of a unary-iteration QROM read over `w`
/// window qubits (`win_base..`) using a reused single-ancilla-per-level spine
/// (`anc_base..`). Data-writes (Clifford `CX`, `2^w·d` of them) are intentionally
/// NOT emitted — they do not affect the Toffoli metric and would be ~290 GB.
/// Produces `2^(w+1) − 4` CCX (matching `ladder_lookup_cost.py`).
fn emit_qrom_read_selector(w: u64, win_base: u64, anc_base: u64) -> Vec<Op> {
    fn rec(level: u64, ctrl: Option<u64>, w: u64, win_base: u64, anc_base: u64, out: &mut Vec<Op>) {
        if level == w {
            return; // leaf: data-writes are Clifford, not emitted (see above)
        }
        let a = anc_base + level;
        let ab = win_base + level;
        // a = ctrl AND addr[level]   (top level: ctrl always-on -> a = addr)
        match ctrl {
            None => out.push(cx(ab, a)),
            Some(c) => out.push(ccx(c, ab, a)),
        }
        rec(level + 1, Some(a), w, win_base, anc_base, out); // addr[level] == 1
                                                             // a ^= ctrl  =>  a = ctrl AND NOT addr[level]
        match ctrl {
            None => out.push(x(a)),
            Some(c) => out.push(cx(c, a)),
        }
        rec(level + 1, Some(a), w, win_base, anc_base, out); // addr[level] == 0
                                                             // restore a = ctrl AND addr[level], then uncompute a -> 0
        match ctrl {
            None => {
                out.push(x(a));
                out.push(cx(ab, a));
            }
            Some(c) => {
                out.push(cx(c, a));
                out.push(ccx(c, ab, a));
            }
        }
    }
    let mut out = Vec::new();
    rec(0, None, w, win_base, anc_base, &mut out);
    out
}

// Heavy: streams (point-add + emitted QROM) x n_add through the analyzers.
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
    let n_add = 2 * 256 / w - 4; // 28 windowed additions (paper A1/A3)
    let n = n_add as usize;

    // EMIT the lookup: a unary-iteration QROM read selector, on ids disjoint from
    // (and reused across windows above) the point-add's register file.
    let win_base = pa_qubits;
    let anc_base = pa_qubits + w;
    let top = anc_base + w; // one past the highest emitted id
    let qrom = emit_qrom_read_selector(w, win_base, anc_base);
    let qrom_read_tof = toffoli_count(qrom.iter());
    // The Rust emission reproduces the Python QROM measurement (ADR 0010).
    assert_eq!(
        qrom_read_tof,
        (1u64 << (w + 1)) - 4,
        "emitted QROM read Toffoli != 2^(w+1)-4"
    );

    // STREAM-EMIT the full reversible ladder [read, point-add, read] x n_add and
    // COUNT Toffoli end-to-end (addition + lookup, no materialization).
    let ladder = || (0..n).flat_map(|_| qrom.iter().chain(ops.iter()).chain(qrom.iter()));
    let ladder_tof_rev = toffoli_count(ladder());
    // Toffoli-depth of the emitted composed stream. QROM runs on disjoint ids, so
    // it composes in parallel with the additions -> the serial accumulator chain
    // dominates; the real read->add dependency (quantum-addend PA) would add only
    // the negligible ~n_add*2w QROM depth.
    let ladder_tdepth = analyze_depth(ladder(), top as usize, nb).toffoli_depth;

    // Addition term on its own, at the true n_add (composition laws).
    let add_tof = toffoli_count((0..n).flat_map(|_| ops.iter()));
    let (add_qubits, _b, _r, _regs) = analyze_ops((0..n).flat_map(|_| ops.iter()));
    let add_tdepth = analyze_depth((0..n).flat_map(|_| ops.iter()), nq, nb).toffoli_depth;
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

    // Emitted reversible total = per window (read + PA + read).
    assert_eq!(
        ladder_tof_rev,
        n_add * (pa_tof + 2 * qrom_read_tof),
        "emitted reversible ladder Toffoli inconsistent with its parts"
    );

    // MBUC-optimized: the unload uses measurement-based uncompute (Clifford), so
    // its Toffoli is the compute-only subset of a read (= read/2 = 2^w-2). Each
    // subcount is emitted; this is the paper's `3*2^w` lookup, decomposed.
    let lookup_mbuc_per_add = qrom_read_tof + qrom_read_tof / 2; // = 3*2^w - 6
    let total_tof_mbuc = add_tof + n_add * lookup_mbuc_per_add;

    // Peak qubits: reported per A2 (register reuse). A naive disjoint emission
    // would over-count; the true overlap needs the quantum-addend PA.
    let total_qubits_a2 = pa_qubits + w;
    let disjoint_peak = top; // what emitting the lookup on separate ids costs

    let derived_tof = (pa_tof + 3 * (1u64 << w)) * n_add;
    let delta = derived_tof as i128 - total_tof_mbuc as i128; // signed: never underflows

    eprintln!("\n=== issue #4 full-ladder streamed measurement (w={w}, n_add={n_add}) ===");
    eprintln!("  per addition (measured): PA toffoli={pa_tof}  qubits={pa_qubits}  toffoli_depth={pa_tdepth}");
    eprintln!("    (static op-stream CCX/CCZ count, as ladder_composition.rs; ecdlp_estimate.py");
    eprintln!("     uses the executed avg-per-shot from score.json, ~1.36M -> ~43.7M headline.)");
    eprintln!("  LOOKUP (EMITTED unary-iteration QROM selector): read = {qrom_read_tof} Toffoli");
    eprintln!(
        "    (reproduces ladder_lookup_cost.py's 2^(w+1)-4; data-writes are Clifford, not emitted)"
    );
    eprintln!("  ADDITION (EMITTED, streamed {n_add}x): toffoli={add_tof} (={}x PA)  width={add_qubits} (flat)  depth={add_tdepth} (={}x PA)",
        add_tof / pa_tof, add_tdepth / pa_tdepth);
    eprintln!("  ------------------------------------------------------------------");
    eprintln!("  FULL LADDER, emitted+counted end-to-end ([read, add, read] x {n_add}):");
    eprintln!(
        "    Toffoli (reversible unload) = {ladder_tof_rev}  (~{:.1}M)",
        ladder_tof_rev as f64 / 1e6
    );
    eprintln!(
        "    Toffoli (MBUC unload)       = {total_tof_mbuc}  (~{:.1}M)",
        total_tof_mbuc as f64 / 1e6
    );
    eprintln!(
        "    toffoli-depth (measured)    = {ladder_tdepth}  (~{:.2}M; add-dominated)",
        ladder_tdepth as f64 / 1e6
    );
    eprintln!("    peak qubits (A2, reuse)     = {total_qubits_a2}   [disjoint-emit peak = {disjoint_peak}, over-counts]");
    eprintln!(
        "  vs derived (PA+3·2^w)·n_add = {derived_tof}  (Δ = {delta} = 6·n_add, MBUC saving)"
    );

    // Underflow-safe cross-checks: measured MBUC total is at/under the derived
    // headline, off only by the MBUC lookup saving.
    assert!(
        total_tof_mbuc <= derived_tof,
        "measured MBUC total exceeds derived headline"
    );
    assert_eq!(
        delta,
        (6 * n_add) as i128,
        "composition mismatch beyond the MBUC saving"
    );
    // Measured depth is add-dominated (QROM parallel on disjoint ids).
    assert_eq!(
        ladder_tdepth, add_tdepth,
        "composed depth not add-dominated as expected"
    );
    // Sanity vs the paper's published full-ECDLP Toffoli bounds (Appendix A).
    assert!(
        total_tof_mbuc < 70_000_000,
        "unexpectedly above the paper's Low-Gate ECDLP bound"
    );
}
