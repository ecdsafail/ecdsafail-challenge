#!/usr/bin/env python3
"""Empirically measure the exceptional-input rate of the incomplete affine adder
across a faithful windowed ECDLP ladder — the missing evidence for issue #5.

The completeness *argument* (`analysis/completeness_argument.md`, PR #14, Path A)
bounds the `dx=0` exceptional amplitude at `~2/n` per addition by assuming the
running accumulator equidistributes over the group. Its own caveat flags this as
heuristic:

    "Equidistribution is heuristic. The ~1/n collision rate assumes the
     accumulator's x-coordinate is approximately uniform over the superposition."

This module turns that heuristic into a *measured* fact on toy curves. It does two
things:

  Part A — validate the scalar model against a real curve.
    On a small prime-order Weierstrass curve it confirms, exhaustively, the fact
    the whole measurement rests on: for a prime-order group, `[s]P` and `[t]P`
    share an x-coordinate **iff** `t ≡ ±s (mod n)`, and `[0]P = ∞`. So the entire
    exceptional-input structure of the affine adder (`dx=0 ⇔ A ∈ {M,−M}`; ∞ ⇔
    the (0,0) sentinel) is captured *exactly* by scalar arithmetic mod n — no
    floating-point curve arithmetic, no sampling error.

  Part B — measure the exceptional rate on the windowed ladder (exact).
    Using that model, it computes the *exact* distribution of the running
    accumulator at each windowed addition of the combined `[a]P + [b]Q` ladder
    (the paper's per-window precomputed-multiple form, with the direct-lookup
    first window that removes the ∞ start — issue #5 §3), over the full
    power-of-two Shor register `a,b ∈ [0, 2^m)`. From that distribution it reads
    off, with no Monte-Carlo error, the per-addition probability of each
    exceptional branch:
      - `dx=0`  collisions  (A ∈ {M,−M}, both finite)  — the §4 quantity,
      - addend = ∞          (a zero window selects the [0]·P = ∞ table entry),
      - accumulator = ∞     (a partial scalar ≡ 0 mod n).
    Summing per-addition probabilities gives the same union bound §4 uses, so the
    measured total is directly comparable to §4's `28 · 2/n`.

What the measurement shows (see the printed table):
  1. The `dx=0` rate tracks `2/n` up to a small constant — the accumulator does
     NOT pathologically concentrate. This is the empirical backing the argument
     lacked; the 240-bit margin in §4 easily absorbs the measured constant.
  2. The union-bound total is instead dominated by the **zero-window / ∞** terms
     at `~1/2^w` per addition (not `2/n`). Still far below Shor's ~1% tolerance,
     so the Path-A conclusion holds — but it means §4's headline bound implicitly
     assumes the windowed lookup never emits the ∞ table entry (a signed-digit /
     offset encoding, or a controlled skip). That is a sharper, testable
     condition than "equidistribution is heuristic", and is the concrete next
     item for a machine-checked completeness result.

Analysis-only, deterministic, pure-Python (no numpy / no z3). Never touches the
scored circuit.
"""
import sys

# --------------------------------------------------------------------------- #
# Part A — a real prime-order elliptic curve, to validate the scalar model.
# --------------------------------------------------------------------------- #

INF = None  # the point at infinity / (0,0) sentinel


def _legendre(a, p):
    return pow(a % p, (p - 1) // 2, p)


class Curve:
    """y^2 = x^3 + a x + b over F_p (complete, reference-correct addition)."""

    def __init__(self, p, a, b):
        self.p, self.a, self.b = p, a, b

    def is_on(self, pt):
        if pt is INF:
            return True
        x, y = pt
        return (y * y - (x * x * x + self.a * x + self.b)) % self.p == 0

    def add(self, P, Q):
        p = self.p
        if P is INF:
            return Q
        if Q is INF:
            return P
        x1, y1 = P
        x2, y2 = Q
        if x1 == x2 and (y1 + y2) % p == 0:
            return INF  # P == -Q
        if P == Q:
            m = (3 * x1 * x1 + self.a) * pow(2 * y1, -1, p) % p
        else:
            m = (y2 - y1) * pow(x2 - x1, -1, p) % p
        x3 = (m * m - x1 - x2) % p
        y3 = (m * (x1 - x3) - y1) % p
        return (x3, y3)

    def mul(self, k, P):
        R, Q, k = INF, P, k % (self.order if hasattr(self, "order") else 10**9)
        while k:
            if k & 1:
                R = self.add(R, Q)
            Q = self.add(Q, Q)
            k >>= 1
        return R

    def points(self):
        pts = [INF]
        for x in range(self.p):
            rhs = (x * x * x + self.a * x + self.b) % self.p
            if rhs == 0:
                pts.append((x, 0))
            elif _legendre(rhs, self.p) == 1:
                y = pow(rhs, (self.p + 1) // 4, self.p) if self.p % 4 == 3 else _tonelli(rhs, self.p)
                pts.append((x, y))
                pts.append((x, self.p - y))
        return pts


def _tonelli(n, p):
    # Minimal Tonelli-Shanks (only hit for p % 4 == 1 curves).
    if _legendre(n, p) != 1:
        raise ValueError("no sqrt")
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    if s == 1:
        return pow(n, (p + 1) // 4, p)
    z = 2
    while _legendre(z, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        bexp = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, bexp * bexp % p, t * bexp * bexp % p, r * bexp % p
    return r


def find_prime_order_curve():
    """Return (curve, generator, order) for a small prime-order curve."""
    def is_prime(x):
        if x < 2:
            return False
        for d in range(2, int(x**0.5) + 1):
            if x % d == 0:
                return False
        return True

    for p in [199, 211, 223, 227, 229, 233, 239, 251, 263, 271]:
        if p % 4 != 3:
            continue  # keep sqrt trivial
        for a in range(0, 4):
            for b in range(1, 6):
                if (4 * a**3 + 27 * b**2) % p == 0:
                    continue
                c = Curve(p, a, b)
                n = len(c.points())  # includes INF -> full group order
                if is_prime(n) and n > 20:
                    c.order = n
                    gen = next(pt for pt in c.points() if pt is not INF)
                    return c, gen, n
    raise RuntimeError("no small prime-order curve found")


def validate_scalar_model():
    """Confirm: [s]P and [t]P share x  <=>  t == +-s (mod n); [0]P = INF."""
    c, P, n = find_prime_order_curve()
    mult = [c.mul(s, P) for s in range(n)]  # dlog table: index = scalar
    assert mult[0] is INF, "scalar 0 must be the point at infinity"
    assert all(pt is not INF for pt in mult[1:]), "prime order: only 0 -> INF"

    fails = 0
    for s in range(n):
        for t in range(n):
            same_x = (
                mult[s] is not INF
                and mult[t] is not INF
                and mult[s][0] == mult[t][0]
            )
            model = (s != 0) and (t != 0) and ((t - s) % n == 0 or (t + s) % n == 0)
            if same_x != model:
                fails += 1
    ok = fails == 0
    print("Part A — scalar model vs. a real prime-order curve")
    print("-" * 74)
    print(f"  curve y^2 = x^3 + {c.a}x + {c.b}  over F_{c.p},  prime order n = {n}")
    print(f"  checked all {n}x{n} scalar pairs:  x-collision  <=>  t == +-s (mod n)")
    print(f"  [{'ok' if ok else 'XX'}] model exact "
          f"({'0 mismatches' if ok else str(fails) + ' MISMATCHES'})\n")
    return ok


# --------------------------------------------------------------------------- #
# Part B — exact exceptional-rate measurement on the windowed ladder.
# --------------------------------------------------------------------------- #
#
# Group = cyclic of prime order n (Part A shows this models the curve exactly).
# A point is its discrete log s in Z_n; INF is s == 0; [s]P and [t]P share an
# x-coordinate iff t == +-s (mod n).
#
# The combined windowed ladder computes [a]P + [b]Q = [a + d*b]P by adding, one
# per window, a precomputed classical multiple selected by that window's value:
#   window (constant c, value v)  contributes scalar  v * c (mod n),
# where c = 2^{w*i}          for the i-th window of a (base P),
#       c = 2^{w*j} * d      for the j-th window of b (base Q = [d]P).
# There are NO doublings (they are folded into the per-window tables) — exactly
# the "adds a fixed precomputed classical multiple M = P[k]" model of §4.
#
# The Shor register is a power of two: a, b in [0, 2^m), m = w*t, 2^m >= n, so
# window values are exactly uniform on [0, 2^w). The first window is a DIRECT
# LOOKUP that writes the accumulator (issue #5 §3), so ∞ is not the certain start.


def contribution_dist(c, n, w):
    """Distribution over Z_n of  v*c mod n  for v uniform on [0, 2^w).

    Returned as a list of 2^w point masses (scalar values); each has weight
    1/2^w. Kept as the raw list so convolutions stay O(2^w * n), not O(n^2)."""
    return [(v * c) % n for v in range(1 << w)]


def convolve(dist, masses, n, w):
    """dist (length-n prob vector) (+) the uniform mass list -> new length-n vec."""
    out = [0.0] * n
    scale = 1.0 / (1 << w)
    for m in masses:
        # add the m-shift of dist, weighted by 1/2^w
        for y in range(n):
            out[(y + m) % n] += dist[y] * scale
    return out


def measure_ladder(n, w, d, label):
    """Exact per-addition exceptional probabilities for the combined ladder.

    Returns a dict of summed (union-bound) probabilities per exceptional class,
    plus the worst-case accumulator non-uniformity observed."""
    t = 0
    while (1 << (w * t)) < n:  # smallest t with 2^{w t} >= n
        t += 1
    m_bits = w * t

    # window list: P windows i=0..t-1, then Q windows j=0..t-1
    windows = [(pow(2, w * i, n)) for i in range(t)]  # base P
    windows += [(pow(2, w * j, n) * d) % n for j in range(t)]  # base Q = [d]P
    masses = [contribution_dist(c, n, w) for c in windows]

    # accumulator distribution AFTER the direct-lookup first window (window 0),
    # which writes the accumulator instead of adding into ∞.
    dist = [0.0] * n
    for s in masses[0]:
        dist[s] += 1.0 / (1 << w)

    sums = {"dx0": 0.0, "addend_inf": 0.0, "acc_inf": 0.0, "any": 0.0}
    max_nonuniform = 0.0  # max_y | dist[y] - 1/n | * n, over pre-add distributions
    n_adds = 0

    for k in range(1, len(windows)):  # each is one adder addition
        n_adds += 1
        c = windows[k]
        # non-uniformity of the accumulator seen by this addition
        dev = max(abs(dist[y] - 1.0 / n) for y in range(n)) * n
        max_nonuniform = max(max_nonuniform, dev)

        p_add_inf = 1.0 / (1 << w)              # v == 0  -> addend M = ∞
        p_acc_inf = dist[0]                     # A == ∞  (partial scalar ≡ 0)
        # dx=0 (both finite): A ∈ {M, −M}, M = v*c, v != 0. A != 0 automatically.
        p_dx0 = 0.0
        for v in range(1, 1 << w):
            mv = (v * c) % n
            p_dx0 += (dist[mv] + dist[(n - mv) % n])
        p_dx0 /= (1 << w)

        # exact "exceptional at this addition" = union of the three (disjoint:
        # dx=0 requires both finite, so no overlap with the ∞ cases; A=∞ and
        # M=∞ can co-occur, so subtract that overlap).
        p_any = p_add_inf + p_acc_inf - p_add_inf * p_acc_inf + p_dx0

        sums["dx0"] += p_dx0
        sums["addend_inf"] += p_add_inf
        sums["acc_inf"] += p_acc_inf
        sums["any"] += p_any

        dist = convolve(dist, masses[k], n, w)

    sums.update(n=n, w=w, t=t, m_bits=m_bits, n_adds=n_adds,
                max_nonuniform=max_nonuniform, label=label,
                dx0_per_add=sums["dx0"] / n_adds, twoovern=2.0 / n)
    return sums


def main():
    print("=" * 74)
    print(" Empirical completeness: exceptional-input rate of the affine adder")
    print(" across a windowed ECDLP ladder   (issue #5, backs completeness §4)")
    print("=" * 74)
    print()

    model_ok = validate_scalar_model()

    print("Part B — exact per-addition dx=0 collision rate vs. the 2/n heuristic")
    print("-" * 74)
    print("  'nonunif' = max accumulator concentration seen (x the uniform mass);")
    print("  the KEY column is 'ratio' = measured dx=0 rate / (2/n).  It stays O(1)")
    print("  even when nonunif is huge, because the addend multiple sweeps the")
    print("  group -> the collision rate is insensitive to accumulator shape.")
    print()
    hdr = (f"  {'n':>5} {'w':>2} {'adds':>4} {'nonunif':>8}  "
           f"{'dx0/add':>9} {'2/n':>9} {'ratio':>6}  {'Sum dx0':>9}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))

    configs = [
        # (n prime, window w, secret d)
        (1009, 2, 613),
        (1009, 5, 613),
        (2003, 4, 877),
        (4093, 4, 2531),
    ]
    results = []
    for n, w, d in configs:
        r = measure_ladder(n, w, d, f"n={n},w={w}")
        results.append(r)
        print(f"  {r['n']:>5} {r['w']:>2} {r['n_adds']:>4} {r['max_nonuniform']:>8.1f}  "
              f"{r['dx0_per_add']:>9.2e} {r['twoovern']:>9.2e} "
              f"{r['dx0_per_add'] / r['twoovern']:>6.2f}  {r['dx0']:>9.2e}")
    print()
    print("  The other exceptional branches are n-independent per addition:")
    print("    addend = ∞ (zero window):  exactly 1/2^w   (measured, matches)")
    print("    accumulator = ∞ residual:  ~1/2^w early, ~1/n once spread")
    print("  So the ladder's total exceptional amplitude splits into a dx=0 term")
    print("  ~ (#adds)·2/n and a zero-window ∞ term ~ (#adds)/2^w.")
    print()

    # ------------------------------------------------------------------ #
    # Extrapolate the validated per-addition rates to attack parameters.
    # secp256k1: n ~ 2^256, paper's windowed ladder w=16, 28 additions.
    # dx=0 uses the measured O(1) prefactor (<=1 across toy n); zero-window
    # is the exact analytic 1/2^w.  Both are compared to Shor's ~1% budget.
    # ------------------------------------------------------------------ #
    import math

    N_REAL = 2 ** 256          # ~ secp256k1 group order
    W_REAL, ADDS_REAL = 16, 28
    pref = max(r["dx0_per_add"] / r["twoovern"] for r in results)  # conservative
    dx0_real = ADDS_REAL * pref * (2.0 / N_REAL)
    zerowin_real = ADDS_REAL * (1.0 / (1 << W_REAL))
    total_real = dx0_real + zerowin_real

    print("Extrapolation to attack parameters  (n≈2^256, w=16, 28 additions)")
    print("-" * 74)
    print(f"  dx=0 total        ~ 28 · {pref:.2f} · 2/n   = {dx0_real:.2e}  "
          f"(≈ 2^{math.log2(dx0_real):.0f})")
    print(f"  zero-window ∞     ~ 28 / 2^16            = {zerowin_real:.2e}  "
          f"(≈ 2^{math.log2(zerowin_real):.0f})  <- dominant")
    print(f"  total exceptional                        = {total_real:.2e}  "
          f"vs Shor tolerance ~1e-2")
    print()

    # ---- assertions: locked as a regression ---- #
    ok = model_ok
    notes = []

    # (1) equidistribution: dx=0 stays within a small constant of 2/n on every
    #     toy config, independent of (large) accumulator non-uniformity.
    ratios = [r["dx0_per_add"] / r["twoovern"] for r in results]
    c1 = all(0.2 <= x <= 6.0 for x in ratios)
    ok &= c1
    notes.append(f"[{'ok' if c1 else 'XX'}] dx=0 rate in [0.2, 6]x of 2/n on all "
                 f"configs (measured {min(ratios):.2f}-{max(ratios):.2f}x) "
                 f"-> equidistribution validated")

    # (2) addend=∞ probability is exactly 1/2^w per addition (n-independent).
    c2 = all(abs(r["addend_inf"] / r["n_adds"] - 1.0 / (1 << r["w"])) < 1e-12
             for r in results)
    ok &= c2
    notes.append(f"[{'ok' if c2 else 'XX'}] addend=∞ rate == 1/2^w exactly "
                 f"-> zero-window term is analytic")

    # (3) at attack parameters both terms are far below Shor's ~1% budget.
    c3 = total_real < 1e-2 and dx0_real < 1e-60
    ok &= c3
    notes.append(f"[{'ok' if c3 else 'XX'}] extrapolated total {total_real:.2e} < 1e-2 "
                 f"-> Path A holds at n≈2^256 (dx=0 ≈ 2^{math.log2(dx0_real):.0f})")

    # (4) informational: the zero-window ∞ term dominates dx=0 at real params.
    notes.append(f"[--] dominant exceptional term at w=16 is the zero-window ∞ "
                 f"case ({zerowin_real:.2e}), not dx=0 -> §4 must state the "
                 f"lookup-encoding condition (see ADR 0008)")

    print("Findings")
    print("-" * 74)
    for line in notes:
        print("  " + line)
    print()
    print("=" * 74)
    if ok:
        print(" RESULT: the 2/n equidistribution heuristic is validated (dx=0 rate")
        print(" tracks 2/n within a small constant, robustly). Extrapolated to")
        print(" attack parameters the total exceptional amplitude is ~2^-11 << 1%,")
        print(" so completeness Path A holds — but the dominant term is the")
        print(" zero-window ∞ case, a condition §4 must state (ADR 0008).")
        print("=" * 74)
        return 0
    print(" RESULT: FAILURE — a locked expectation regressed (see [XX] above).")
    print("=" * 74)
    return 1


if __name__ == "__main__":
    sys.exit(main())
