# ADR 0015 — Offset window encoding removes the zero-window ∞ exceptional term

**Status:** Accepted — dominant exceptional term (zero-window ∞) removed
structurally; completeness bound sharpened from ~2⁻¹¹ back to the dx=0-limited ~2⁻²⁵⁰
**Date:** 2026-07-02

## Context

[ADR 0008](0008-empirical-completeness-collision-rate.md) / #15 measured the
incomplete affine adder's exceptional-input rate across the combined `[a]P+[b]Q`
windowed ladder ([issue #5](https://github.com/CaptainEmpower/ecdsafail-challenge/issues/5))
and split it into two very different terms:

- **`dx=0` collisions** (`A ∈ {M, −M}`, both finite): `~2/n` per addition →
  ladder total `~28·2/n ≈ 2⁻²⁵⁰`.
- **zero-window ∞** (a window digit is `0`, so the lookup selects the `[0]·P = ∞`
  table entry as the addend): exactly `1/2^w` per addition → ladder total
  `~28/2¹⁶ ≈ 2⁻¹¹` — ~240 bits larger, i.e. the **dominant** exceptional term.

Both sit ≪ Shor's ~1% tolerance, so Path A's conclusion is unchanged either way,
but the `2⁻²⁵⁰` headline in `completeness_argument.md §4` silently assumes the
lookup never emits the ∞ entry. #15 named the fix as the concrete next step:
"a signed-digit / offset encoding, or a controlled skip of zero windows".

## Decision

Add `analysis/verify/offset_window_encoding.py` (suite stage 8/10) implementing
and validating the **offset encoding**. Standard base-`2^w` windowing adds
`[digit · 2^{w i}]·P` per window, which is `∞` when `digit = 0`. Shift every digit
by one and add `[(digit+1) · 2^{w i}]·P` instead: the emitted index
`digit+1 ∈ [1, 2^w]` is never `0`, so for `2^w < n` (the real `w=16 ≪ n≈2²⁵⁶`) the
addend is a finite point at every window — the ∞ entry is never selected. Over the
`t` windows of a base this adds a fixed classical constant `S = Σ_{i} 2^{w i}`, so
the offset ladder computes `[a]P+[b]Q + [(1+d)S]P`; subtracting the single
compile-time point `[(1+d)S]P` (folded into precomputation) recovers `[a]P+[b]Q`.

- **Part A** validates on a real prime-order toy curve, exhaustively: standard
  windowing emits ∞ exactly at a zero digit, the offset encoding **never** emits ∞,
  and the corrected offset ladder yields `[a]P+[b]Q` for every `(a, b)`.
- **Part B** re-runs #15's exact distribution measurement with the offset digit
  set `v ∈ [1, 2^w]` and reads off the exceptional rates.

The alternative (a controlled skip of zero windows) reaches the same end but costs
a per-window comparator and data-dependent control; the offset encoding is a
compile-time reindexing plus one classical correction, so it is preferred.

## Consequences

- The `addend = ∞` (zero-window) rate is **exactly 0** on every measured config —
  the dominant exceptional term is removed *structurally*, not by negligibility.
- The `dx=0` rate is unchanged (`O(1)·2/n`; the offset multiple still sweeps the
  group), so the ladder total returns to the `dx=0`-limited `~2⁻²⁵⁰` — the
  `completeness_argument.md §4` headline now holds under an explicit, validated
  encoding condition rather than a silent assumption.
- Extrapolated (w=16, 28 additions): standard zero-window ∞ `≈ 2⁻¹¹` → offset `0`;
  offset total exceptional `≈ 2⁻²⁴⁸` (dx=0 only). Still ≪ Shor's tolerance, but
  now the *dominant* term is eliminated rather than merely bounded.
- Analysis-only, deterministic, pure-Python; reuses the scalar model + exact
  convolution validated in #15. Consistent with
  [ADR 0001](0001-analysis-layer-isolated-from-score.md): no effect on the scored
  circuit. What remains for #5 is the circuit-level **mid-ladder** verification
  over the real superposition, which needs the quantum-addend testbed shared with
  [issue #4](https://github.com/CaptainEmpower/ecdsafail-challenge/issues/4) /
  [ADR 0011](0011-streamed-full-ladder.md).
