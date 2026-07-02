# ADR 0007 — Tier B: measuring the full-ECDLP ladder

**Status:** Accepted (strategy); first increment landed, quantum-addend/QROM/QFT deferred
**Date:** 2026-07-02

## Context

The full-ECDLP cost (`analysis/ecdlp_estimate.py`, [ADR 0003](0003-ground-ecdlp-estimate-in-source-paper.md))
is *derived*: `ECDLP_Toff = (PA_Toff + 3·2^w)(2n/w − 4)`, `ECDLP_Qubits =
PA_Qubits + w`, with `PA_*` the measured per-addition primitive and the rest
analytic. Issue #4 asks to replace the derived composition with emitted-and-
measured numbers.

Findings that shape what is feasible:

- **The additions compose faithfully as classical-feedforward-controlled copies.**
  The scored primitive takes the addend as a *classical* `BitId` register
  (`ec_add(circ, x2, y2, ox, oy)`), and Shor's qubit-recycled (semiclassical) QFT
  measures each scalar/window bit and feeds it forward as a classical condition —
  exactly the repo's `c_condition` mechanism. So a controlled addition has the
  **same Toffoli count** as an uncontrolled one, and the ladder's addition term
  is `n_add · PA_Toff` by additivity (classical addends don't change gate
  structure). This is not an approximation; it is exact for the addition term.
- **What is genuinely un-measured** is therefore only: (i) the windowed **QROM
  lookup** that selects `P[k]` — the scored primitive uses a classical addend, so
  a faithful windowed ladder needs a *quantum-addend* point-add variant (or the
  paper's measurement-based unary-iteration lookup writing into the addend
  register), which is new circuitry; and (ii) the **QFT** (Clifford, ~0 Toffoli).
  Both are small relative to the dominant `28 · PA_Toff` term.

## Decision

Measure what composes faithfully now; build the rest as a tracked follow-on.

- **First increment (landed):** measure the ADDITION composition directly on the
  already-optimized single-addition op stream by chaining its iterator `k` times
  through `circuit::analyze_ops` / `analyze_depth` (no pipeline re-run, no memory
  blow-up). This turns three previously-asserted claims into measured facts:
  additivity of Toffoli (`k·PA`), **flat peak width** across `k` (ancilla reused
  ⇒ confirms `ECDLP_Qubits = PA_Qubits + w`), and **serial toffoli-depth**
  (`≈ k·PA_depth`, accumulator serializes). Lives in a `#[cfg(test)]` harness;
  circuit byte-identical.
- **Deferred (the real remaining build):** a quantum-addend point-add variant +
  the paper's optimized unary-iteration QROM + the semiclassical QFT, composed
  into the actual 28-window ladder and measured end-to-end. This is where the
  lookup term stops being derived. Large; tracked in #4.

## Consequences

- The dominant term (`28 · PA`) is now measured/justified, not asserted; the
  estimate's headline (~43.7M Toffoli, ~1168 qubits) rests on measured addition
  composition + measured per-addition + a derived lookup/QFT tail.
- The remaining derived pieces (QROM `3·2^w`, QFT) are explicitly the smaller
  terms; refining them needs the quantum-addend circuitry, not more analysis.
- Depends on nothing new for the first increment; the deferred build is the
  natural home for the #5 ∞-accumulator structural fix (the "first window =
  direct lookup" that removes the ∞ start) — so #4 and #5 converge there.
- Consistent with [ADR 0001](0001-analysis-layer-isolated-from-score.md):
  measurement harness is `#[cfg(test)]`, score-neutral.
