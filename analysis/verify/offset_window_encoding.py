#!/usr/bin/env python3
"""Remove the *dominant* exceptional term of the windowed ladder — the
zero-window ∞ addend — with an offset window encoding, and re-measure the
sharpened completeness bound (issue #5, part (b)).

Background. `completeness_collision_rate.py` (#15, ADR 0008) measured the affine
adder's exceptional-input rate across the combined `[a]P + [b]Q` windowed ladder
and found two very different terms:

  - `dx=0` collisions (A ∈ {M, −M}, both finite):  ~`2/n` per addition
        -> ladder total ~`28 · 2/n ≈ 2⁻²⁵⁰`, and
  - **zero-window ∞** (a window digit is 0, so the table lookup selects the
        `[0]·P = ∞` entry as the addend):  exactly `1/2^w` per addition
        -> ladder total ~`28 / 2¹⁶ ≈ 2⁻¹¹`  <-  ~240 bits LARGER, the dominant term.

Both sit ≪ Shor's ~1% tolerance, so Path A's conclusion holds either way — but the
`2⁻²⁵⁰` headline in `completeness_argument.md §4` implicitly assumes the windowed
lookup never emits the ∞ table entry. #15 named the fix as the concrete next step:
"a signed-digit / offset encoding, or a controlled skip of zero windows". This
module implements and validates the **offset encoding**, which removes the
zero-window ∞ term *structurally* (not by negligibility), so the honest bound
returns to the `dx=0`-limited `~2⁻²⁵⁰`.

The offset encoding. Standard base-`2^w` windowing adds, for window `i` of base
`P`, the multiple `[digit · 2^{w i}]·P`; when `digit = 0` that multiple is `∞`
(the zero-window addend). Shift every digit by one:

    add  [(digit + 1) · 2^{w i}]·P     instead of   [digit · 2^{w i}]·P

Now the emitted index `digit + 1 ∈ [1, 2^w]` is never `0`, so (for `2^w < n`,
i.e. the real `w=16 ≪ n≈2²⁵⁶`) the addend is a *finite* point at every window —
the ∞ table entry is never selected. Summed over the `t` windows of a base this
adds a fixed classical constant `S = Σ_{i=0}^{t−1} 2^{w i}`, so the offset ladder
computes `[a + S]P + [b + S]Q = [a]P + [b]Q + [(1+d)S]P`; subtracting the single
compile-time point `[(1+d)S]P` (folded into precomputation) recovers `[a]P+[b]Q`.

Part A validates the encoding on a real prime-order toy curve, exhaustively:
correct combined point for every `(a, b)`, and **no window ever emits ∞** (while
standard windowing does, exactly when a digit is 0). Part B re-runs the exact
distribution measurement of #15 with the offset digit set and confirms the
`addend=∞` rate is now **exactly 0**, `dx=0` is unchanged, and the extrapolated
w=16 total drops from `~2⁻¹¹` to `~2⁻²⁵⁰`.

Analysis-only, deterministic, pure-Python (no numpy / z3); reuses the scalar
model + exact-convolution machinery validated in #15. Never touches the scored
circuit.
"""
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from completeness_collision_rate import (  # noqa: E402
    INF,
    convolve,
    find_prime_order_curve,
    measure_ladder,
)


# --------------------------------------------------------------------------- #
# Part A — the offset encoding is correct and ∞-free (exhaustive, real curve).
# --------------------------------------------------------------------------- #

def window_count(n, w):
    """Smallest t with 2^{w t} >= n, so every scalar in [0, n) is representable
    by t base-2^w windows."""
    if w <= 0:
        raise ValueError(f"window width w must be positive, got {w}")
    t = 0
    while (1 << (w * t)) < n:
        t += 1
    return t


def validate_offset_encoding(w=2):
    """Exhaustively check, on a real prime-order toy curve, that the offset
    windowed decomposition (a) never emits the ∞ table entry, while standard
    windowing does exactly at a zero digit, and (b) computes [a]P + [b]Q for
    every (a, b) after the single classical correction. Returns (ok, detail)."""
    c, P, n = find_prime_order_curve()
    # dlog table over the REAL curve: mult[s] = [s]P. Ties every scalar claim
    # below to actual curve geometry (mult[s] is INF iff s ≡ 0, prime order).
    mult = [c.mul(s, P) for s in range(n)]
    assert mult[0] is INF and all(pt is not INF for pt in mult[1:])

    t = window_count(n, w)
    if (1 << w) >= n:
        return False, f"toy w={w}: need 2^w < n for a structurally ∞-free offset"
    mask = (1 << w) - 1
    d = 2  # a nontrivial secret; Q = [d]P
    S = sum(1 << (w * i) for i in range(t))       # per-base offset constant
    corr = ((1 + d) * S) % n                       # combined classical correction

    # (a) ∞-emission per window value, curve-level: standard hits ∞ exactly at
    #     digit 0; offset never does (index digit+1 ∈ [1, 2^w], nonzero mod n).
    for i in range(t):
        base_c = (1 << (w * i)) % n
        for g in range(1 << w):
            std_scalar = (g * base_c) % n          # standard addend multiple
            off_scalar = ((g + 1) * base_c) % n    # offset addend multiple
            if (mult[std_scalar] is INF) != (g == 0):
                return False, f"standard window {i} digit {g}: ∞ not iff digit=0"
            if mult[off_scalar] is INF:
                return False, f"offset window {i} digit {g}: emitted ∞ (index {g+1})"

    # (b) exhaustive combined correctness over every (a, b) in [0, n).
    for a in range(n):
        for b in range(n):
            acc = 0                                # accumulated dlog (offset ladder)
            for i in range(t):                     # base-P windows
                g = (a >> (w * i)) & mask
                acc = (acc + (g + 1) * (1 << (w * i))) % n
            for j in range(t):                     # base-Q windows (Q = [d]P)
                g = (b >> (w * j)) & mask
                acc = (acc + (g + 1) * (1 << (w * j)) * d) % n
            acc = (acc - corr) % n                  # subtract the classical offset
            if mult[acc] != mult[(a + d * b) % n]:  # points must match
                return False, f"offset ladder wrong at a={a}, b={b}"

    return True, (f"curve y²=x³+{c.a}x+{c.b}/F_{c.p}, n={n}, w={w}, {t} windows/base: "
                  f"offset never emits ∞; [a]P+[b]Q correct for all {n}² (a,b)")


# --------------------------------------------------------------------------- #
# Part B — re-measure the exceptional rates with the offset digit set.
# --------------------------------------------------------------------------- #
#
# Same exact-distribution method as measure_ladder() in #15, but each window
# adds v*c for v in [1, 2^w] (the offset digit set) instead of [0, 2^w). Because
# c is a unit mod the prime n, v*c ≡ 0 iff v ≡ 0 (mod n); for 2^w < n no v in
# [1, 2^w] hits that, so the addend=∞ branch has probability exactly 0.

def measure_ladder_offset(n, w, d, label):
    """Exact per-addition exceptional probabilities for the OFFSET ladder.

    Mirrors completeness_collision_rate.measure_ladder but with the offset digit
    set v ∈ [1, 2^w]; returns the same dict shape (Fractions)."""
    t = window_count(n, w)
    windows = [pow(2, w * i, n) for i in range(t)]           # base P
    windows += [(pow(2, w * j, n) * d) % n for j in range(t)]  # base Q = [d]P
    # offset masses: v in [1, 2^w] (never the zero digit)
    masses = [[(v * c) % n for v in range(1, (1 << w) + 1)] for c in windows]

    dist = [0] * n                                   # accumulator after init write
    for s in masses[0]:
        dist[s] += 1
    denom = 1 << w

    sums = {"dx0": Fraction(0), "addend_inf": Fraction(0),
            "acc_inf": Fraction(0), "any": Fraction(0)}
    n_adds = 0
    for k in range(1, len(windows)):
        n_adds += 1
        c = windows[k]
        # addend = ∞ iff v*c ≡ 0 (mod n), v ∈ [1, 2^w]: exactly 0 for 2^w < n.
        cnt_add_inf = sum(1 for v in range(1, (1 << w) + 1) if (v * c) % n == 0)
        p_add_inf = Fraction(cnt_add_inf, 1 << w)
        p_acc_inf = Fraction(dist[0], denom)
        cnt_dx0 = 0
        for v in range(1, (1 << w) + 1):
            mv = (v * c) % n
            if mv == 0:                              # that is the addend=∞ branch
                continue
            cnt_dx0 += dist[mv] + dist[(n - mv) % n]
        p_dx0 = Fraction(cnt_dx0, denom * (1 << w))
        p_any = p_add_inf + p_acc_inf - p_add_inf * p_acc_inf + p_dx0

        sums["dx0"] += p_dx0
        sums["addend_inf"] += p_add_inf
        sums["acc_inf"] += p_acc_inf
        sums["any"] += p_any

        dist = convolve(dist, masses[k], n)
        denom *= (1 << w)

    sums.update(n=n, w=w, t=t, n_adds=n_adds, label=label,
                dx0_per_add=sums["dx0"] / n_adds, twoovern=Fraction(2, n))
    return sums


CONFIGS = [
    # (n prime, window w, secret d) — matches completeness_collision_rate.py
    (1009, 2, 613),
    (1009, 5, 613),
    (2003, 4, 877),
    (4093, 4, 2531),
]


def main():
    print("=" * 74)
    print(" Offset window encoding: remove the dominant zero-window ∞ exceptional")
    print(" term of the ladder, sharpen the completeness bound (issue #5b)")
    print("=" * 74)
    print()

    ok = True

    # Part A ---------------------------------------------------------------- #
    good, detail = validate_offset_encoding(w=2)
    ok &= good
    print("Part A — offset encoding is ∞-free and correct (exhaustive, real curve)")
    print("-" * 74)
    print(f"  [{'ok' if good else 'XX'}] {detail}\n")

    # Part B ---------------------------------------------------------------- #
    print("Part B — exact exceptional rates: standard vs offset digit set")
    print("-" * 74)
    print("  'addend∞/add' = per-addition zero-window ∞ rate (the dominant term);")
    print("  offset drives it to EXACTLY 0.  'dx0/add' is unchanged (~O(1)·2/n).")
    print()
    hdr = (f"  {'n':>5} {'w':>2} {'adds':>4} | {'std ∞/add':>10} {'off ∞/add':>10} "
           f"| {'std dx0/add':>11} {'off dx0/add':>11} {'2/n':>9}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))

    results = []
    for n, w, d in CONFIGS:
        std = measure_ladder(n, w, d, f"n={n},w={w}")
        off = measure_ladder_offset(n, w, d, f"n={n},w={w}")
        results.append((std, off))
        std_inf = std["addend_inf"] / std["n_adds"]
        off_inf = off["addend_inf"] / off["n_adds"]
        print(f"  {n:>5} {w:>2} {off['n_adds']:>4} | "
              f"{float(std_inf):>10.3e} {float(off_inf):>10.3e} | "
              f"{float(std['dx0_per_add']):>11.3e} {float(off['dx0_per_add']):>11.3e} "
              f"{float(off['twoovern']):>9.2e}")
    print()

    # Extrapolation to attack parameters ------------------------------------ #
    import math
    N_REAL, W_REAL, ADDS_REAL = 2 ** 256, 16, 28
    pref = float(max(o["dx0_per_add"] / o["twoovern"] for _, o in results))
    dx0_real = ADDS_REAL * pref * (2.0 / N_REAL)
    std_zerowin = ADDS_REAL * (1.0 / (1 << W_REAL))       # standard: dominant term
    print("Extrapolation to attack parameters  (n≈2²⁵⁶, w=16, 28 additions)")
    print("-" * 74)
    print(f"  standard: zero-window ∞ total ~ 28/2¹⁶      = {std_zerowin:.2e}  "
          f"(≈ 2^{math.log2(std_zerowin):.0f})  <- was dominant")
    print(f"  offset:   zero-window ∞ total               = 0          "
          f"(structurally ∞-free)")
    print(f"  offset:   total exceptional  ~ dx=0 only    = {dx0_real:.2e}  "
          f"(≈ 2^{math.log2(dx0_real):.0f})")
    print(f"  -> offset sharpens the headline from ~2⁻¹¹ back to the dx=0-limited "
          f"~2⁻²⁵⁰.")
    print()

    # Regression / finding locks -------------------------------------------- #
    notes = []
    off_inf_zero = all(o["addend_inf"] == 0 for _, o in results)
    ok &= off_inf_zero
    notes.append(f"[{'ok' if off_inf_zero else 'XX'}] offset addend=∞ rate == 0 "
                 f"exactly on every config (zero-window term removed)")

    std_inf_pos = all(s["addend_inf"] > 0 for s, _ in results)
    ok &= std_inf_pos
    notes.append(f"[{'ok' if std_inf_pos else 'XX'}] standard addend=∞ rate > 0 "
                 f"(the term the offset removes was really there)")

    ratios = [float(o["dx0_per_add"] / o["twoovern"]) for _, o in results]
    dx0_ok = all(0.2 <= r <= 6.0 for r in ratios)
    ok &= dx0_ok
    notes.append(f"[{'ok' if dx0_ok else 'XX'}] offset dx=0 rate still in [0.2, 6]× of "
                 f"2/n ({min(ratios):.2f}–{max(ratios):.2f}×) — unaffected by the shift")

    sharper = std_zerowin > dx0_real and dx0_real < 1e-60
    ok &= sharper
    notes.append(f"[{'ok' if sharper else 'XX'}] extrapolated offset total "
                 f"{dx0_real:.2e} ≪ standard {std_zerowin:.2e} — bound sharpened "
                 f"to ≈2^{math.log2(dx0_real):.0f}")

    print("Findings")
    print("-" * 74)
    for ln in notes:
        print("  " + ln)
    print()
    print("=" * 74)
    if ok:
        print(" RESULT: the offset window encoding removes the dominant zero-window ∞")
        print(" exceptional term STRUCTURALLY (addend never ∞), verified ∞-free and")
        print(" correct on a real toy curve. The completeness headline sharpens from")
        print(" ~2⁻¹¹ back to the dx=0-limited ~2⁻²⁵⁰ (issue #5b; ADR 0015).")
        print("=" * 74)
        return 0
    print(" RESULT: FAILURE — see [XX] above.")
    print("=" * 74)
    return 1


if __name__ == "__main__":
    sys.exit(main())
