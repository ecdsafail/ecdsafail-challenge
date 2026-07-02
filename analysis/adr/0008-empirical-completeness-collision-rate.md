# ADR 0008 — Empirically validate (and sharpen) the completeness collision rate

**Status:** Accepted — equidistribution heuristic validated on toy curves; dominant
exceptional term re-identified as the zero-window ∞ case
**Date:** 2026-07-02

## Context

The completeness negligibility argument ([issue #5](https://github.com/CaptainEmpower/ecdsafail-challenge/issues/5),
`analysis/completeness_argument.md`, Path A of [ADR 0006](0006-adder-completeness-approach.md))
bounds the incomplete affine adder's exceptional amplitude at `≈ 2⁻²⁵⁰` — a union
bound of `28 · 2/n` over the windowed ladder, where the per-addition `2/n` is the
probability the running accumulator `A` collides in x-coordinate with the fixed
addend `M` (`A ∈ {M, −M}`, the `dx=0` branch). The argument's own caveat flags
the weak link:

> "Equidistribution is heuristic. The `~1/n` collision rate assumes the
> accumulator's x-coordinate is approximately uniform over the superposition."

That assumption was asserted, not measured. This ADR records the decision to
*measure* it, and what the measurement changed.

## Decision

Add an exact, self-contained measurement (`analysis/verify/completeness_collision_rate.py`,
suite stage 5/7) instead of leaving the rate as a heuristic:

- **Scalar model, validated against a real curve.** For a prime-order group,
  `[s]P` and `[t]P` share an x-coordinate **iff** `t ≡ ±s (mod n)`, and `[0]P = ∞`.
  Part A confirms this exhaustively on a small prime-order Weierstrass curve, so
  the whole exceptional-input structure (`dx=0 ⇔ A ∈ {M,−M}`; ∞ ⇔ the `(0,0)`
  sentinel) is captured **exactly** by scalar arithmetic mod `n` — no
  floating-point curve arithmetic, no Monte-Carlo error.
- **Exact per-addition rates by convolution.** Part B computes the exact
  distribution of the running accumulator at each windowed addition of the
  combined `[a]P + [b]Q` ladder (per-window precomputed multiples, direct-lookup
  first window per #5 §3), over the full power-of-two Shor register. It reads off
  the exact probability of each exceptional branch: `dx=0`, `addend = ∞` (a zero
  window selects the `[0]·P = ∞` table entry), and `accumulator = ∞`.
- **Validate, then extrapolate.** Confirm the mechanism on toy `n`, then extend
  the validated per-addition rates to attack parameters (`n ≈ 2²⁵⁶`, `w = 16`,
  28 additions).

## Consequences

- **The `2/n` heuristic holds — and is robust.** Measured `dx=0` rate is
  `0.47–0.81 × 2/n` across `n ∈ {1009, 2003, 4093}` and `w ∈ {2,4,5}`. Notably it
  stays `O(1)·2/n` **even when the accumulator is far from uniform** (concentration
  up to `250×` the uniform mass): the addend multiple `M = v·c` sweeps the group
  as the window value `v` varies, so the collision rate is insensitive to the
  accumulator's shape. This is a *stronger* footing than the original
  "accumulator ≈ uniform" argument. Extrapolated, `dx=0` total `≈ 28·2/n ≈ 2⁻²⁵⁰`,
  matching §4.
- **The dominant exceptional term is not `dx=0`.** The `addend = ∞` (zero-window)
  branch has probability **exactly `1/2^w`** per addition — `n`-independent —
  giving a ladder total `≈ 28/2¹⁶ ≈ 2⁻¹¹`, over 200 bits larger than the `dx=0`
  term. It is still `≪` Shor's `~1%` tolerance, so **Path A's conclusion is
  unchanged** (`completeness_overhead = 1.0` stands). But it means §4's `2⁻²⁵⁰`
  headline silently assumed the windowed lookup never emits the `∞` table entry.
  The honest bound is `≈ 2⁻¹¹`, and it is conditional on the lookup encoding
  (signed-digit / offset representation avoiding zero windows, or a controlled
  skip of zero windows). Stating that condition is the concrete next step toward a
  machine-checked completeness result — sharper than "equidistribution is
  heuristic".
- Consistent with [ADR 0001](0001-analysis-layer-isolated-from-score.md): the
  measurement is analysis-only Python, deterministic, no numpy/z3 dependency, and
  cannot affect the scored circuit.
- Open follow-on (tracked in #5): confirm the actual lookup encoding of the #4
  Tier-B ladder either avoids zero-window `∞` addends or handles them, which
  turns the `≈ 2⁻¹¹` conditional bound into an unconditional one.
