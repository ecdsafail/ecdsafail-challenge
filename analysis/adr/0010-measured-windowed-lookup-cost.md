# ADR 0010 — Measure the windowed-lookup (QROM) cost (Tier B, issue #4)

**Status:** Accepted — lookup term grounded by a validated unary-iteration QROM;
end-to-end ladder composition still deferred
**Date:** 2026-07-02

## Context

The full-ECDLP estimate (`analysis/ecdlp_estimate.py`,
[ADR 0003](0003-ground-ecdlp-estimate-in-source-paper.md)) is
`ECDLP_Toff = (PA_Toff + 3·2^w)(2n/w − 4)`. Two of the three factors are already
grounded: `PA_Toff` is **measured** from the scored circuit, and the
addition-composition law (`n_add · PA`, flat width, serial depth) was **measured**
in [ADR 0007](0007-tier-b-measured-ladder.md) / issue #4's first increment. The
remaining derived piece is the **`3·2^w` table-lookup term** — the windowed QROM
that loads the precomputed multiple `P[k]` from a `2^w`-entry table indexed by the
`w`-qubit window register. The scored primitive takes a *classical* addend, so it
contains **no** QROM; nothing in the repo measured this term. It was cited from the
paper.

## Decision

Ground the lookup term with a validated construction rather than leaving it cited,
without building the full quantum-addend ladder:
`analysis/verify/ladder_lookup_cost.py` (suite stage 5/8).

- Build the read as an **optimized unary-iteration QROM** (Gidney 2018 §III.C, the
  primitive the paper cites) — a single reused-ancilla "spine" recursion over the
  `w` address bits, `out ^= T[addr]` — as an actual kickmix circuit.
- **Validate** it exhaustively with the reference-faithful `kickmix_sim.py`
  ([ADR 0004](0004-cross-validate-against-reference-circuits.md)): correct read for
  every address × random tables, address/table registers unchanged, all `w`
  iteration ancilla returned to `|0⟩`, global phase `+1`.
- **Measure** its Toffoli count and ancilla width, and compare to `3·2^w`.

## Consequences

- The lookup term now has a validated circuit and a **measured** cost:
  `Toffoli(read) = 2^(w+1) − 4` with `w` ancilla — e.g. `131,068` Toffoli / `16`
  ancilla at the optimal `w = 16`, versus the paper's `3·2^16 = 196,608`. So the
  read is `~0.67×` the paper's per-addition budget, and the ancilla scaling `w`
  matches `ECDLP_Qubits = PA_Qubits + w` (A2). `ecdlp_estimate.py` keeps `3·2^w`
  as the **conservative headline** but now prints the measured grounding.
- The compute half of the read (`2^w − 2` Toffoli) is a subset of the same
  validated circuit; with measurement-based uncomputation (the repo's / paper's
  technique, as the adders and `controlled_lookup.py` use) the uncompute becomes
  Clifford, so a windowed addition's read + uncompute-read stays near `2·2^w`,
  within the `3·2^w` budget.
- **Still deferred** (the remaining Tier B work, issue #4): composing the
  quantum-addend point-add + this QROM + the semiclassical QFT into the actual
  28-window ladder and stream-measuring the end-to-end Toffoli/qubit/depth (the
  full `Vec` is ~290 GB, so it must be stream-emitted). This ADR grounds the
  lookup *term*; the end-to-end emission is the next increment. It is also where
  issue #5's mid-ladder residual and zero-window-`∞` question land.
- Consistent with [ADR 0001](0001-analysis-layer-isolated-from-score.md):
  analysis-only, pure-Python, deterministic, reuses the validated kickmix
  simulator; zero effect on the scored circuit.
