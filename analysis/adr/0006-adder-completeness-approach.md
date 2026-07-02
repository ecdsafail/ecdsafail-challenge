# ADR 0006 — Approach to adder completeness (cost estimate → verified attack)

**Status:** Accepted (strategy); the A-vs-B path choice is gated by an experiment
**Date:** 2026-07-02

## Context

The ECDLP numbers (`analysis/ecdlp_estimate.py`, [ADR 0003](0003-ground-ecdlp-estimate-in-source-paper.md))
carry `completeness_overhead = 1.0` — i.e. they assume the exceptional cases of
elliptic-curve addition are handled at negligible cost. That assumption is what
keeps the result a *cost estimate* rather than a *verified attack* (issue #5).

Grounding in the code:

- `trailmix_ludicrous/ec_add.rs` implements **only the generic chord formula**
  (it forms `dx = x2 − ox` and inverts it). The reference `add()`
  (`weierstrass_elliptic_curve.rs:67–121`) branches on four cases; only the
  generic one is built. Missing: P = ∞, Q = ∞, doubling (P = Q), and P = −Q → ∞.
- Two distinct failure modes result:
  - **∞ operand** (a point is the (0,0) sentinel): `dx ≠ 0`, no crash, but the
    generic formula returns a **wrong affine point**.
  - **`dx = 0`** (x1 == x2: doubling or P = −Q): modular inverse of 0 is
    undefined → the GCD pass returns garbage. **Open:** does the circuit leave
    *clean* ancilla + correct phase (⇒ merely wrong output) or *dirty* state
    (⇒ reversibility-breaking)?
- The **scored** circuit is unaffected: `eval_circuit.rs:245–266` draws P, Q as
  random multiples of G and skips ∞, so `dx = 0` has prob ~2⁻²⁵⁵ and never
  occurs in 9024 shots. Completeness only matters in a full Shor ladder over
  superposition (where the accumulator *starts* at ∞), so this work depends on
  #4 (Tier B ladder) as a testbed.

## Decision

Pursue completeness **analysis-first (Path A)**, matching the source paper, and
fall back to complete formulas (Path B) only if A cannot be established:

- **Path A — negligibility (Roetteler–Naehrig–Svore–Lauter 2017).** Show the
  wrong-amplitude fraction stays below Shor's failure tolerance (the paper only
  proves ~99% correctness via Fiat–Shamir fuzz; Shor tolerates a small bad
  fraction). The ∞-accumulator is not rare initially and must be removed
  structurally — the paper's "first windowed addition = direct lookup" trick
  (Appendix A) appears to do exactly this; confirm in the #4 ladder. Then bound
  the residual `dx = 0` amplitude.
- **Path B — complete formulas (Renes–Costello–Batina 2016).** Exception-free
  adder; unconditionally correct, but new circuitry (~1.5–2× cost) that changes
  the cost model and perturbs the optimized primitive. Reserved as fallback.

**Gating experiment (do first, before committing to A vs B):** run the *existing*
circuit on crafted `dx = 0` and ∞ inputs via a `#[cfg(test)]` harness in
`point_add`, and record whether the ancilla/phase stay clean. Clean ⇒ exceptions
are "wrong output" only (Path A viable); dirty ⇒ reversibility-breaking (pushes
toward Path B or an explicit guard).

## Consequences

- Gives a decision procedure rather than a premature commitment; the experiment
  is cheap (no #4 dependency, no score impact — tests don't run in
  `build_circuit`) and highly informative.
- If Path A succeeds, the ECDLP result upgrades from cost estimate to verified
  attack with no change to the optimized circuit. If Path B is forced, the cost
  model must be re-derived.
- Consistent with [ADR 0001](0001-analysis-layer-isolated-from-score.md): the
  experiment and any analysis are score-neutral; only a Path-B rebuild would
  touch the scored circuit, which would be a separate, explicit decision.
- Full definition of done: a rigorous (proof- or simulation-backed) argument that
  full-ladder success meets Shor's threshold, or a validated exception-free
  adder. Tracked in issue #5.
