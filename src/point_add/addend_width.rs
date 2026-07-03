//! Tier B (issue #27), second increment: MEASURE the *width* half of the
//! classical-vs-quantum-addend gap — the `ECDLP_Qubits = PA_Qubits + w` (A2)
//! assumption that `ladder_full.rs` reports but flags as "the true overlap needs
//! the quantum-addend PA" (lines 38–47, 195–198), and that ADR 0012 left open
//! after closing the *Toffoli* half (that gap is ≤ 0.05% of PA).
//!
//! ## The finding: the width gap is REAL (opposite of the Toffoli gap)
//!
//! The scored PA keeps the other point `Q = (ox, oy)` **classical**. `coord_addsub`
//! loads a coordinate into a fresh 256-qubit temp *at each coordinate step*, runs
//! the q-q add, then **frees it within the step** — so the addend is resident only
//! at the (off-peak) coordinate phases, never at the GCD peak. Measured profile
//! (`TRACE_PHASE_ACTIVE=1 TRACE_TLM_PROFILE=1 cargo run --release --bin build_circuit`):
//!
//! ```text
//!   peak  active_max = 1152  phase = tlm_multiply_gcd_reverse_body   (the GCD apply)
//!   coord active_max = 1026  phase = tlm_coord_{x,y}_sub / add3x / rsub_final
//!                            (1026 already INCLUDES the 256-qubit addend temp)
//! ```
//!
//! A QROM-fed **quantum** addend cannot do this. `ox` is consumed at steps 3, 7, 15
//! and `oy` at steps 4, 14 (`ec_add.rs`) — both **straddle both GCD passes and the
//! square**, i.e. the peak. A single QROM read of `P[k]` must therefore stay
//! resident *across the peak*, where — crucially — it cannot overlap the GCD
//! scratch, because the preserved addend (needed *after* the GCD) and the GCD
//! scratch (in use *during* it) are **both live at the peak**. And the scored
//! build's peak is *tight*: its max qubit id (`analyze_ops` = 1152 = `score.json`
//! qubits) equals the profiler's peak active (1152), so at the peak the free list
//! is empty — there is no freed slot for the addend to reuse.
//!
//! => a faithful quantum-addend port adds the resident addend's **full width** on
//! top of PA_Qubits: `+256` if one coordinate is held while the other is re-read
//! per side (lower bound — one coord must survive the peak), up to `+512` for a
//! single `P[k]=(x,y)` read held across the whole addition. Plus the window
//! register (`+w`, the A2 term) and the QROM's `w` unary-iteration ancilla.
//!
//! So for THIS product-min PA, `ECDLP_Qubits = PA_Qubits + w` (= 1168) **undercounts**:
//! the quantum-addend width is `PA_Qubits + (256..512) + w` ≈ 1408..1664. A2 holds
//! for the *paper's* PA because its ZK bound (1175 low-qubit / 1425 low-gate)
//! already prices a resident quantum addend into a tighter arithmetic core; this
//! repo spent its width budget on a product-min GCD and stayed under bound by
//! keeping the addend classical — an advantage that a faithful port would erase.
//!
//! ## What this harness checks — and what carries the conclusion
//!
//! `analyze_ops(..).0` is the **allocation span** (`max qubit id + 1`), NOT a
//! live/peak-active count. So this harness does not *independently* measure a
//! peak: it builds the scored circuit and CONSTRUCTS a resident-addend port by
//! shifting every qubit id up by the addend width `Δ` (placing a held `Δ`-qubit
//! addend register at ids `[0, Δ)`, loaded/uncomputed at the ends). Its span is
//! then `PA_span + Δ` essentially by construction — the harness only confirms the
//! span *bookkeeping* (a held `Δ`-qubit register that cannot reuse a freed id adds
//! exactly `Δ` ids; the `NO_QUBIT` sentinel is preserved under the shift).
//!
//! The substantive **peak** conclusion — that a resident quantum addend raises the
//! *peak active* qubit count by its full width — rests on the separately-**measured
//! profiler fact** quoted above: the scored build's peak active (`1152`) equals its
//! max qubit id (`analyze_ops` = `score.json` qubits = `1152`). That equality means
//! the peak is *tight* (the free list is empty at the peak), so a register that
//! must stay live across the peak has no freed slot to reuse and therefore adds new
//! ids to the peak, not just to the span. The harness is the reproducible span-side
//! of ladder_full.rs's flagged "disjoint-emit over-counts" caveat; the tight-peak
//! profiler measurement is what turns a span delta into a peak delta.
//!
//! `#[cfg(test)]` only; never compiled into the scored circuit.

use super::build;
use crate::circuit::{analyze_ops, Op, OperationType, QubitId, NO_QUBIT};

/// Shift a single qubit id by `delta`, leaving the `NO_QUBIT` sentinel untouched.
fn shift_q(q: QubitId, delta: u64) -> QubitId {
    if q == NO_QUBIT {
        q
    } else {
        QubitId(q.0 + delta)
    }
}

/// Model a quantum-addend port that holds a `delta`-qubit addend register resident
/// across the WHOLE addition: place it at ids `[0, delta)` and run the scored
/// circuit on ids shifted up by `delta`. An `X` on each addend qubit at the start
/// (QROM load) and end (uncompute) makes the register live end-to-end. The port's
/// allocation span is then `analyze_ops(port).0` (= max qubit id + 1), which is
/// `PA_span + delta` by construction — see the module note on why, given the
/// measured tight peak, this span delta is also a peak-active delta.
fn resident_addend_port(scored: &[Op], delta: u64) -> Vec<Op> {
    let mut out: Vec<Op> = Vec::with_capacity(scored.len() + 2 * delta as usize);
    let x = |t: u64| {
        let mut o = Op::empty();
        o.kind = OperationType::X;
        o.q_target = QubitId(t);
        o
    };
    for id in 0..delta {
        out.push(x(id)); // QROM load of the resident addend register
    }
    for op in scored {
        let mut o = *op;
        o.q_control1 = shift_q(o.q_control1, delta);
        o.q_control2 = shift_q(o.q_control2, delta);
        o.q_target = shift_q(o.q_target, delta);
        out.push(o);
    }
    for id in 0..delta {
        out.push(x(id)); // uncompute the addend register
    }
    out
}

#[test]
#[ignore = "heavy; builds the scored circuit, run with --ignored --exact"]
fn quantum_addend_width_gap() {
    const N: u64 = 256; // one secp256k1 coordinate
    const W: u64 = 16; // paper's optimal window (A2's `+w` term)

    let scored = build();
    let (pa_qubits, _b, _r, _regs) = analyze_ops(scored.iter());
    assert!(pa_qubits > 0, "degenerate build");

    // Construct + measure the resident-addend port width for the two anchor cases.
    let port_one = resident_addend_port(&scored, N); // hold one coordinate (lower bound)
    let (q_one, _b, _r, _regs) = analyze_ops(port_one.iter());
    let port_both = resident_addend_port(&scored, 2 * N); // hold P[k] = (x, y)
    let (q_both, _b, _r, _regs) = analyze_ops(port_both.iter());

    // A2 headline vs the measured quantum-addend port.
    let a2 = pa_qubits + W; // ECDLP_Qubits = PA_Qubits + w (paper A2)
    let port_one_full = q_one + W; // + window register
    let port_both_full = q_both + W;

    eprintln!("\n=== issue #27 quantum-addend WIDTH gap (A2) ===");
    eprintln!("  scored PA span (analyze_ops = max id + 1) : qubits={pa_qubits}");
    eprintln!("    (= score.json qubits, and separately == profiler peak ACTIVE, so the");
    eprintln!("     peak is tight: free list empty at the peak -> a resident addend held");
    eprintln!("     across the peak cannot reuse a freed slot, so a span delta is a peak delta)");
    eprintln!("  classical addend: resident only off-peak (coord steps ~1026 < {pa_qubits} peak),");
    eprintln!("    re-loaded+freed per step -> never coexists with the GCD scratch.");
    eprintln!("  quantum addend (QROM P[k]) must straddle the peak (ox@3/7/15, oy@4/14):");
    eprintln!("    hold one coord  : port span={q_one}  (= PA + {N})");
    eprintln!(
        "    hold P[k]=(x,y) : port span={q_both}  (= PA + {})",
        2 * N
    );
    eprintln!("  A2 headline PA+w = {a2};  quantum-addend port span = {port_one_full}..{port_both_full} (+w)");
    eprintln!(
        "  => A2's `+w only` UNDERCOUNTS this classical-addend PA by {N}..{} qubits.",
        2 * N
    );
    eprintln!("     (the paper's PA_Qubits bound already prices a resident addend; this repo");
    eprintln!(
        "      stayed under bound by keeping the addend classical — a port would erase that.)"
    );

    // Checked findings (span bookkeeping; the peak conclusion is argued from the
    // tight-peak profiler equality documented in the module note above):
    // (1) A held Δ-qubit addend register adds its full width to the allocation span
    //     PA + Δ. Given the measured tight peak (peak active == max id == PA), this
    //     span delta is also a peak-active delta — the addend cannot overlap the
    //     GCD scratch, since both must be live across the (already full) peak.
    assert_eq!(
        q_one,
        pa_qubits + N,
        "holding one coordinate did not add its full {N} qubits to the span"
    );
    assert_eq!(
        q_both,
        pa_qubits + 2 * N,
        "holding P[k]=(x,y) did not add its full {} qubits to the span",
        2 * N
    );
    // (2) The quantum-addend width therefore exceeds A2's PA+w by at least one
    //     resident coordinate — the width caveat is real, not negligible.
    assert!(
        port_one_full > a2 + N - 1,
        "measured quantum-addend port ({port_one_full}) not materially above A2 ({a2})"
    );
}
