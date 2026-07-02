# ADR 0009 — Circuit-demonstrate the ∞-start removal (direct-lookup first window)

**Status:** Accepted — amplitude-1 ∞ start shown removed by an actual reversible
lookup; mid-ladder residual deferred to the Tier B ladder (#4)
**Date:** 2026-07-02

## Context

Completeness (issue #5) splits into two exceptional-input concerns, and they are
*not* of equal weight ([ADR 0006](0006-adder-completeness-approach.md),
[ADR 0008](0008-empirical-completeness-collision-rate.md),
`completeness_argument.md`):

- **(a) The ∞-accumulator start.** The running accumulator begins at ∞ with
  **amplitude 1** — certain, not rare. Negligibility cannot touch it; if the
  adder is fed ∞ at t=0 the attack fails outright. It must be removed
  *structurally*. The source paper (Babbush et al. 2026, Appendix A) does so by
  replacing the **first windowed addition with a direct table lookup** that
  *writes* the accumulator instead of adding into ∞.
- **(b) The zero-window / mid-ladder residual.** `addend = ∞` (a zero window) and
  `dx = 0` collisions, both far below Shor's ~1% tolerance
  ([ADR 0008](0008-empirical-completeness-collision-rate.md): `~2⁻¹¹` and
  `~2⁻²⁵⁰` respectively).

Part (a) is the load-bearing correctness claim, and until now it rested only on
the paper's prose and on being *assumed* by `completeness_argument.md §3` /
`completeness_collision_rate.py`. Nothing in the repo exercised the init as an
actual circuit.

## Decision

Demonstrate (a) at circuit level with a self-contained, score-neutral verifier
(`analysis/verify/direct_lookup_init.py`, suite stage 6/8), reusing primitives
already validated rather than building new circuitry:

- Use the **controlled table-lookup QROM** validated in
  [ADR 0005](0005-validate-lookup-by-construction.md) /
  `controlled_lookup.py` (issue #3/#9) as the init write `acc ^= T[w]` on an
  accumulator register that starts at `|0…0⟩`, with the classical table
  `T[w] = affine coords of [w]·P`.
- Verify, via the reference-faithful `kickmix_sim.py`, that for **every** window
  value `w`: the register holds exactly the coordinates of `[w]·P` (write
  correct), the selector ancilla return to `|0⟩`, the global phase is `+1`, and —
  the (a) property — the register is the `(0,0)` ∞ sentinel **iff `w = 0`**.
- Include a `ctrl = 0` **negative control** (no write ⇒ accumulator stays ∞) so
  it is demonstrably the *write* that removes ∞, and run both the reversible and
  MBUC uncompute forms.
- Scope the exhaustive check to a small prime-order toy curve (the property
  `[w]P = ∞ ⇔ w ≡ 0` is curve-independent), plus a **secp256k1 256-bit spot
  check** to show the write composes at attack scale.

## Consequences

- The amplitude-1 ∞ start is now **removed by construction, circuit-demonstrated**
  — not assumed. For the definite initial window value the accumulator is a real
  affine point, so the adder is never fed ∞ at t=0. The only residual is the
  `w = 0` zero-window case, which is part (b) and already negligible.
- This closes the **load-bearing half** of issue #5. What remains is the
  *mid-ladder* residual (∞ / `dx=0`) across the real 28-window two-scalar
  superposition, which is already bounded analytically
  ([ADR 0008](0008-empirical-completeness-collision-rate.md)) but whose
  end-to-end circuit confirmation needs the **Tier B ladder**
  ([ADR 0007](0007-tier-b-measured-ladder.md), issue #4) — the quantum-addend
  point-add + unary-iteration QROM + semiclassical QFT. Part (b)'s open question
  (does the windowed lookup emit ∞ at `w = 0`, or does the encoding avoid it?)
  falls out of that same testbed.
- Consistent with [ADR 0001](0001-analysis-layer-isolated-from-score.md):
  analysis-only, pure-Python, deterministic, reuses validated primitives; zero
  effect on the scored circuit.
