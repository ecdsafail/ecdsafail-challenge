# Completeness negligibility argument (issue #5, Path A)

This turns the ECDLP result from a *cost estimate* toward a *verified attack* by
arguing — quantitatively — that the incomplete affine adder this repo implements
is sufficient for a working Shor-ECDLP run, i.e. the exceptional cases it
mishandles occur with amplitude far below Shor's tolerance. This is the
Roetteler–Naehrig–Svore–Lauter 2017 style argument, made concrete with this
repo's measured behaviour ([ADR 0006](adr/0006-adder-completeness-approach.md)).

It is an **argument, not a machine-checked proof**; the caveats at the end state
exactly what is heuristic.

## 1. Shor's tolerance to a small wrong fraction

Shor's period-finding succeeds as long as all but a small fraction of the
computational basis states are computed correctly (in value **and** phase). The
source paper (Babbush et al. 2026, Appendix A.5) proves only ~99% correctness via
Fiat–Shamir fuzz and relies on exactly this: "a superposition with 1% of points
in the wrong place will cause the algorithm to fail at most 1% of the time." So
it suffices to show the total amplitude landing on an exceptional adder input is
≪ 1% (and can be repeated a few times to overcome any residual failure).

## 2. What the adder does on exceptional inputs (measured)

`src/point_add/completeness_probe.rs` ran the built circuit on crafted
exceptional inputs (16 RNG seeds). The signature is uniform:

| input | ancilla | output | phase |
|---|---|---|---|
| doubling (dx=0) | clean 16/16 | wrong | corrupted on ~9/16 |
| P=−Q (dx=0) | clean 16/16 | wrong | corrupted on ~4/16 |
| ∞ accumulator | clean 16/16 | wrong | corrupted on ~8/16 |

Key facts used below: exceptional inputs **do not leak the ancilla** (no
register-basis corruption of *other* basis states — the error is confined to the
offending state), but they **do** inject wrong output and probabilistic phase
garbage. So each exceptional basis state contributes **at most its own
amplitude** to Shor's failure probability. It is therefore enough to bound the
summed amplitude of exceptional inputs.

## 3. The ∞-accumulator must be removed structurally (not by negligibility)

The running accumulator starts at ∞ with **amplitude 1** (before any addition),
so ∞ is *not* rare at the start and cannot be waved away. The paper's ladder
removes it structurally: the **first windowed addition is replaced by a direct
table lookup** (Appendix A) that *writes* the initial accumulator instead of
adding into ∞. Hence the adder is never fed ∞ as the accumulator at t=0.

After the first window, the accumulator equals `[a']P + [b']Q` for the
partial scalars accumulated so far; it is ∞ only when that partial scalar is
`≡ 0 (mod n)`, i.e. **one** value out of the group order `n ≈ 2²⁵⁶` — amplitude
`~1/n ≈ 2⁻²⁵⁶` per addition, which falls under the §4 bound. So the ∞ case is
handled: amplitude-1 start removed by construction, residual occurrences
negligible.

## 4. The `dx=0` collisions are negligible

An addition adds a *fixed* precomputed classical multiple `M = P[k]`. It hits the
`dx=0` branch (x-coordinates equal) exactly when the accumulator `A` satisfies
`A.x == M.x`, i.e. `A ∈ {M, −M}` (the two points sharing that x). Over the
superposition of scalars the accumulator ranges over ~`n` group points; treating
it as approximately equidistributed, the amplitude with `A ∈ {M, −M}` is

```
P[dx=0 at one addition]  ≈  2 / n  ≈  2⁻²⁵⁵
```

The windowed ladder performs `2n/w − 4 = 28` windowed additions (w=16), so by a
union bound the **total** exceptional amplitude across the whole run is

```
28 · 2/n  ≈  56 / n  ≈  2⁻²⁵⁰   ≪   10⁻²  (Shor's tolerance)
```

— over 240 bits of margin. Doubling (`A == M`) and `P=−Q` (`A == −M`) are the two
sub-cases and are already included in the `A ∈ {M,−M}` count.

## 5. Conclusion

Combining §3 and §4: after the direct-lookup first window removes the amplitude-1
∞ start, the total amplitude on any exceptional adder input across the full
28-addition ladder is `≈ 2⁻²⁵⁰`, over 240 bits below Shor's ~1% tolerance.
The incomplete affine adder is therefore sufficient for a working attack — no
complete formulas (Path B) are required — matching the standard argument in the
literature. This is what justifies `completeness_overhead = 1.0` in
`analysis/ecdlp_estimate.py`.

## Caveats (what keeps this an argument, not a proof)

- **Equidistribution is heuristic.** The `~1/n` collision rate assumes the
  accumulator's x-coordinate is approximately uniform over the superposition. A
  rigorous proof would bound the actual distribution of partial-scalar multiples
  (or invoke a specific ladder ordering that provably avoids `{M,−M}`), as
  Roetteler et al. discuss. The 240-bit margin is large enough that even a very
  non-uniform distribution stays negligible, but this is not machine-checked.
- **The ∞-removal depends on the ladder using direct-lookup initialisation.**
  This repo builds one addition, not the ladder; §3 relies on the paper's
  structure. Confirming it requires the Tier B ladder build ([issue #4](https://github.com/CaptainEmpower/ecdsafail-challenge/issues/4)),
  which is where this fix lands.
- **Phase, not just value.** §2 shows exceptions also corrupt phase; the bound in
  §4 covers this because a mis-phased state at amplitude ε contributes ≤ ε to the
  failure probability, exactly as a wrong-value state does.
