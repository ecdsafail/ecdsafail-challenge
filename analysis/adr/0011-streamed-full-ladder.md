# ADR 0011 — Stream-emit and measure the full ECDLP ladder (Tier B, issue #4)

**Status:** Accepted — ladder cost emitted+measured under the paper's model; the
quantum-addend point-add remains the one stated assumption
**Date:** 2026-07-02

## Context

Issue #4 asks to replace the *derived* full-ECDLP totals
(`analysis/ecdlp_estimate.py`, `(PA_Toff + 3·2^w)(2n/w − 4)`) with emitted-and-
measured numbers. The full ladder is ~5×10⁹ ops (a materialized `Vec` is ~290 GB),
so it must be **stream-emitted and counted**, never built. Prior increments
grounded the pieces:

- `PA` (per point-addition) — measured from the scored circuit.
- Addition **composition law** — measured to k ≤ 4 in `ladder_composition.rs`
  ([ADR 0007](0007-tier-b-measured-ladder.md)): Toffoli additive, peak width flat,
  toffoli-depth serial — but *extrapolated* to the real `n_add = 28`.
- **Lookup** term — measured in `ladder_lookup_cost.py`
  ([ADR 0010](0010-measured-windowed-lookup-cost.md)): a unary-iteration QROM read
  is `2^(w+1)−4` Toffoli, its MBUC uncompute-read `2^w−2`.

What was missing is the **end-to-end composition at the true ladder length**.

## Decision

Add `src/point_add/ladder_full.rs` (`#[cfg(test)]`, `#[ignore]` heavy): stream the
built point-add op stream `n_add = 2n/w − 4` times through `analyze_ops` /
`analyze_depth` (same qubit ids ⇒ real cross-copy hazards + ancilla reuse; no
materialization), and compose the measured addition stream with the measured
lookup term and the (Clifford) QFT.

## Consequences

- The ladder is now **emitted-and-measured**, at the real `n_add`, not
  extrapolated. Measured composition laws hold exactly at `n_add = 28`
  (asserted): addition `= 28·PA`, peak width **flat** (`= PA_qubits`),
  toffoli-depth `= 28·PA_depth`. Headline (static op-stream basis, w=16):
  - **Toffoli** `= 28·PA + 28·(3·2^16 − 6) ≈ 46.0M`
    (addition `40.5M` + lookup `5.50M`; QFT 0),
  - **peak qubits** `= PA_qubits + w = 1168`,
  - **toffoli-depth** `= 28·PA_depth = 30.16M` (serial additions dominate; the
    QROM iteration depth `~2w` is negligible).
- It matches the derived headline to within the MBUC lookup saving (`6·n_add`),
  and never exceeds it — cross-validating `ecdlp_estimate.py` by construction.
  (The harness reports the *static* op-stream Toffoli, `PA ≈ 1.446M`, as
  `ladder_composition.rs` does; `ecdlp_estimate.py` uses the *executed* avg-per-
  shot from `score.json`, `≈ 1.364M`, hence its `~43.7M` vs the `~46.0M` static
  figure — same model, different PA basis, not a discrepancy.)
- **One stated assumption remains** (the genuine remaining Tier B build): the
  **quantum-addend point-add**. This repo's PA folds a *classical*, compile-time
  addend (constprop-optimized); the windowed ladder loads `P[k]` from a *quantum*
  table, so a quantum-addend PA variant would be needed for a functionally
  complete ladder, and it is what would confirm the register overlap behind
  `ECDLP_Qubits = PA_Qubits + w` (A2) — i.e. that the lookup's `w` window + `w`
  spine ancilla + addend register truly reuse PA's workspace. This harness prices
  the composition under the paper's model; it does not build that variant.
- This is also where issue #5's mid-ladder residual (∞ / `dx=0` over the real
  28-window superposition) and the zero-window-`∞` encoding question land — the
  same quantum-addend testbed resolves both.
- Consistent with [ADR 0001](0001-analysis-layer-isolated-from-score.md): the
  harness is `#[cfg(test)]`, streamed, and never compiled into `build_circuit`;
  the scored circuit is byte-identical (`ops.bin` SHA unchanged).
