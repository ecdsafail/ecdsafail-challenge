# ADR 0005 — Validate the ladder lookup primitive by construction

**Status:** Accepted
**Date:** 2026-07-02

## Context

[ADR 0004](0004-cross-validate-against-reference-circuits.md) validated the
source paper's reference *adder* circuits against an independent kickmix
simulator. The other ECDLP-ladder primitive is the windowed table lookup — the
`3·2^w` term in the cost formula ([ADR 0003](0003-ground-ecdlp-estimate-in-source-paper.md)).
The paper ships `table_lookup_3x3.kmx` for it, but investigation (issue #3)
established:

- Our `verify/kickmix_sim.py` is equivalent instruction-for-instruction to the
  reference simulator `original/zkp_ecc_zenodo_v2/lib/src/sim.rs`, so a no-op
  result is not a simulator bug.
- `table_lookup_3x3.kmx` ships only as `.kmx` + `.svg` (no test-case / fuzzer /
  proof, unlike every iadd circuit): it is an **illustrative extract**. Its
  unary-iteration selector accumulator is `R`-reset to `|0⟩` and would be driven
  by an outer control absent from the standalone snippet; no accumulator/register
  /layout combination recovers correctness (< 4/8). It is therefore not
  fuzz-validatable as-is.

## Decision

Validate the lookup primitive by **constructing** a self-contained controlled
lookup rather than relying on the non-runnable extract. `verify/controlled_lookup.py`
generates a kickmix circuit computing

    r0 ^= (ctrl ? r2[r1] : 0)

for parameterized address width `a`, data width `d`, a `2^a`-entry classical
table `r2`, and an explicit control qubit — in two forms: a reversible
CCX-ladder uncompute, and a measurement-based-uncomputation (HMR + CZ phase
correction) uncompute matching the paper's technique. It fuzz-validates both
(exhaustive addresses × both control values × random tables) for: correct
output, a no-op when `ctrl = 0`, all selector ancilla returned to `|0⟩`, and
global phase `+1`. Wired into `analysis/run.sh` (stage 4/6).

## Consequences

- The ladder's lookup primitive is now validated in a form we control, closing
  the gap left by the non-runnable reference extract — both the reversible and
  the MBUC (paper-technique) uncomputations are exercised for lookups, not just
  adders.
- The test is shown to have teeth: removing the MBUC `CZ` phase corrections makes
  the phase check fail 63/64.
- The constructed circuit is a *correct* lookup, not the paper's minimal-Toffoli
  one; it validates the function and the uncomputation discipline, not the exact
  `3·2^w` gate count. The reference extract remains a diagram (issue #3 stays open
  for reproducing it in full controlled context if ever needed).
- Consistent with [ADR 0001](0001-analysis-layer-isolated-from-score.md):
  deterministic (seeded), analysis-only, no effect on the scored circuit.
