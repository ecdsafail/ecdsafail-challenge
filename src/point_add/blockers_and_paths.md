# Structural blockers and candidate paths to SOTA

Written 2026-04-27 after deep reading of primary sources. All numeric
claims traced to either our measurements or a specific paper (cited).

## Setting

- Primary metric: average executed Toffoli (lower better).
- Secondary constraint: peak qubits ≤ 1300 (stretch) or ≤ 2800 (cap).
- Correctness: Shor-tolerant approximation OK (user-confirmed, Option A).
  ≤ 0.1% failure rate across 9024 shots, padding residue allowed.
- Current baseline: **4,136,878 Toffoli @ 2716 qubits**, exact-correct
  (commit c1aeeb4).
- Target: **2.7M Toffoli @ ≤ 1300 qubits** per Google 2026 SOTA.

## What the SOTA actually is

Google 2026 ZK-proves two circuits for the windowed point-add `Q += P[k]`:

| Variant    | Toffoli | Qubits | Ops | Description |
|------------|--------:|-------:|----:|-------------|
| Low-qubit  |  2.7M   |  1175  | 17M | w=16 windowed point-add + table lookups |
| Low-gate   |  2.1M   |  1425  | 17M | Same task, different tradeoff |

Our task is **simpler** than their windowed task (classical Q directly, no
lookups). Lookup overhead in Google's 2.7M is ~200k Toffoli
(3 lookups × (2^16 + √(2^16)) CCX). So Google's implied **bare** point-add
cost is ~2.5M.

Public literature best at n=256 (per point-add):
- HRSL 2020 Low-W: ~12M at 2124q
- HRSL 2020 Low-T: ~12M at 2619q
- Kim 2026 Q-opt: ~17M at 6475q (depth-optimized)
- Litinski 2023+2024: ~3.7M at 6000q (single instance)
- **Us: 4.14M at 2716q** (simpler bare task)
- **Google 2026: 2.5M (bare) / 2.7M (windowed) at ~1175q (undisclosed)**

## Blocker 1: second modular inversion

**Our current circuit performs 2 full Kaliski invocations (forward+backward
per invocation = ~1.61M Toffoli per invocation).**

Each Kaliski invocation needs both a forward and backward pass for
reversibility. In our circuit: 2 × 1.61M = 3.22M Toffoli from Kaliski
alone, out of 4.13M total (78%).

**Fundamental claim**: an in-place affine point-add `(Px, Py) → (Rx, Ry)`
requires information about `dx = Px−Qx` to be either live OR reconstructible
when the inversion's ancilla must be zeroed. Since Px is overwritten by
Rx by the end, and Rx doesn't trivially yield dx without another inversion,
any single-invocation scheme hits a cleanup obstruction unless:

1. **Output goes to fresh registers** (not in-place). Then the cleanup
   step still leaves stale dx somewhere, requiring a SECOND invocation
   to zero it. Net: still 2 invocations.

2. **Approximate correctness** relaxes "exactly zero" to "approximately
   zero" (≤10⁻⁷ per invocation). Then partial Kaliski backward or
   windowed-tail techniques might suffice. **UNDER-EXPLORED.**

3. **Structurally different algebra** that avoids the cleanup dependency.
   Strategy C (`w = dx³`) classically works (tested 200/200), but its
   quantum schedule still needs dx alive at cleanup time. Not resolved.

**Conclusion**: the "save 1 inversion" lever is ~1.6M Toffoli (a 40%
reduction). It's partially viable under approximate correctness but
requires novel scheduling work that nobody has publicly done.

## Blocker 2: Kaliski per-iter cost

Our per-iter cost ≈ 1950 CCX/iter (measured). Breakdown:

| sub-phase              | CCX/iter (fwd) | fraction |
|------------------------|---------------:|---------:|
| step 4 (sub + add)     | ~1237          | 63%      |
| step 3 + step 9 cswaps | ~1240 combined | 63% (duplicate with 4) |
| step 2 (comparator)    | ~390           | 20%      |
| step 0 (equality zero) | ~70            | 3.6%     |
| step 1                 | ~72            | 3.7%     |
| step 6/7/8 (halve/double) | ~116        | 6%       |

Wait these don't sum right. Re-check: these are AVG per iter, some
phases are bulk-only (313 of 407 iters). Let me recompute:

Forward Kaliski total = 802k CCX over 407 iters = 1970 CCX/iter average.
Sum of phases divided by their iter counts:
- bulk phases (313 iters): step4 + 2 cswaps + step2 = (503 + 252 + 252 + 158)·1e3/313 ≈ 3720 CCX/bulk-iter
- slow phases (94 iters non-bulk): step4 + 2 cswaps + etc = ~2500 CCX/slow-iter

Average ~3200 CCX/iter. Hmm, off from 1970. The trace includes both
forward AND backward in "bulk" phases, so bulk step4 is fwd+bwd. Let
me just use **forward-only** figures by dividing by 814:

Ah right: the bulk phases are labeled `kal_*` for forward or `bk_*`
for backward. kal_bulk_step4 = 503k over 313 bulk iters = 1606 CCX/iter
in BULK. Plus slow step4 is kal_step4 = 151k over (407-313)=94 slow
iters = 1607 CCX/iter. Consistent!

**So per-iter Kaliski forward cost ≈ 1600 CCX/iter for step 4 alone.**
Plus step 3 cswap + step 9 cswap = 2 × (252+76)/407 = 1610 CCX/iter.
Plus step 2 = 616 CCX/iter.
Plus step 0/1/5/6/7/8/10 ≈ 500 CCX/iter.

Sum ≈ 4300 CCX/iter forward. With 407 iters: 1.75M. We measure 802k
for FORWARD only. So actual per-iter ~2000 CCX/iter fwd. Average of
bulk+slow: (bulk step 4 (313 × 1237) + slow step 4 (94 × 371)) =
(387k + 35k)/(407) = 1037 CCX/iter avg for step 4 fwd. That matches
the 500k over 407 iter breakdown. Etc.

**Bottom line**: our forward Kaliski is ~800k CCX / 407 iters = **1970
CCX/iter**.

### Where per-iter savings might come from

1. **Eliminate step 3 + step 9 cswap on (u, v_w)** via HRSL swap-based
   Kaliski. Published HRSL Algorithm 2:
   ```
   bswap ← false
   if u even and v odd, or both odd and u > v:
     swap u and v; swap r and s; bswap ← true
   if u and v both odd: v ← v − u; s ← r + s
   v ← v/2; r ← 2·r
   if bswap: swap u and v; swap r and s
   ```
   HRSL's version still has 2 pairs of cswaps per iter (bswap before and
   after the sub/halve). Same as ours structurally. **No immediate
   savings from direct HRSL port.**

2. **Step 4's fused sub+add** could share more work. Currently:
   - Load tmp = add_f AND u (n CCX)
   - Sub v_w -= tmp (~n CCX)
   - Transform tmp to add_f AND r (n CX + n CCX + n CX) = n CCX 
   - Add s += tmp (n CCX)
   - Unload tmp via MBU (0 CCX + n HMR + n CZ_if)
   Total: 4n CCX/iter. At n=256: 1024 CCX/iter just for step 4.
   **Already optimized.** Could maybe save more via in-place swap logic.

3. **Coset representation** on the modular portion of step 4's
   `s += r mod p`. Saves ~60% on THAT specific op = ~400 CCX/iter.
   Over 407 iters × 2 passes = **~325k CCX saved**. A real 8% reduction.
   But coset breaks comparators, which the Kaliski body relies on.

4. **Windowed-Montgomery style** replacement for Kaliski body. Litinski
   2023 uses Montgomery mul with 9n+28 CCX per 4-bit step, n/4 steps.
   That's 2.25n² per mul. At n=256: 147k per mul. We're at 67k per mul
   already. Not obviously a Kaliski win.

5. **Eliminate m_hist via recomputation** (verified classically in prior
   session, formula `m = f AND u[0] AND (NOT v_w[0] OR (u > v_w))`).
   Saves 407 qubits at ~negligible Toffoli cost. **Unblocks karatsuba
   at pair1_mul2: -28k Toffoli.** Total impact maybe 1-2% of Toffoli,
   large qubit impact.

## Blocker 3: step 3 + step 9 cswaps

Total cswap cost: 1.31M Toffoli (32% of baseline). Each cswap = 1 CCX + 2 CX.

**These are genuinely unavoidable in the current algorithm.** Each pair
of cswaps implements the conditional branching of Kaliski (do-op-on-
smaller-register). HRSL's reformulation keeps the same 2 cswap-
positions per iter.

**Possible angles**:
1. Don't actually swap; track position via a classical bit and
   conditionally select which register each op writes to. This costs
   1 bit of state per iter (= our existing `a_f`), so no free lunch.
   BUT: a SELECT instead of a SWAP could be cheaper. Need to work out
   exactly what operations need to be selected on.

2. Batch multiple iters before a single cswap. If we do 2 iters before
   the swap, the swap might not even be needed (pairs of swaps cancel).
   **Plausible new angle**. Worth thinking about.

3. Use approximation: allow "wrong register" sometimes with bounded
   error, save swap cost. Probably doesn't compose well.

## Blocker 4: 4×n Kaliski register layout

Our Kaliski state: u, v_w, r, s = 4n = 1024 qubits. Plus m_hist (407)
+ other ancillae. Total state during body: ~1432 persistent.

With data registers tx, ty (512q) live simultaneously, we're already
at 1944q persistent during pair1 body. Add ~260q transient (step 4
tmp + Cuccaro) = 2204q. Measured peak 2716q (extra ~500 from lam_inner
and other concurrent state).

For 1300q target: need total non-data ≤ 788q. Means Kaliski state must
fit in 788q. Can't reduce 4×n core state (that's info-theoretically
needed for the (u, v, r, s) Bezout recursion). But can:

1. **Eliminate m_hist**: 407 → 0 via recomputation. Saves 407q
   persistent. Moves us to 1800q peak estimate. Still over 1300.

2. **Fuse (u, v) and/or (r, s)** when their combined bitlen fits in
   n bits. Per Kaliski invariant: bitlen(u) + bitlen(v) ≤ 2n − iter,
   so at iter ≥ n they fit in n. Saves ~n qubits in second half of
   pass. Doesn't help peak (which is at mid-pass).

3. **Reuse tx as v_w**: tx starts holding dx = v, which Kaliski reads
   into v_w via CX-copy. If we skip the copy and use tx directly as
   v_w, save 256q. After Kaliski, tx is back to holding dx (by
   Kaliski's invariant at termination). **Feasible, unexplored.**

4. **Fold r onto tx at end**: after Kaliski, we want tx to hold the
   inverse (or a derived value). If we don't allocate a separate r
   register but compute into tx directly, save 256q. **Requires
   restructuring Kaliski internal layout.**

Combined (1+3+4): saves 407 + 256 + 256 = 919q from current 2716.
Projected peak: **~1800q**. Still above 1300.

5. **Share with ty**: ty is idle during the inversion body (before pair1
   body closure reads it). If ty doubles as Kaliski s register, save
   another 256q. **Requires scheduling care.**

Combined (1+3+4+5): saves 1175q. Projected peak: **~1540q**. Closer.

6. **Windowed Kaliski** (small-window r, s): if r, s fit in ~n/w bits
   per window-step, save most of (r, s) storage. **Needs novel circuit,
   not in public literature.**

With all of above: **1300q plausibly reachable**. Toffoli impact: m_hist
elimination unlocks 2-level Karatsuba savings (~28k), but the layout
changes themselves don't reduce Toffoli. **To hit 2.7M, we still need
the one-inversion or cost-reduction stacking.**

## Blocker 5: comparator vs coset representation

Coset representation cuts mod-add cost by ~60% but breaks comparators.
Our Kaliski uses `u > v_w` (step 2) and `v_w == 0` (step 0) — both
comparators. **Directly porting to coset breaks Kaliski.**

**Possible resolutions**:
1. **Exit coset before comparator, re-enter after.** Cost: ~cpad
   Toffoli per round-trip (small), plus the comparison itself. Net:
   small overhead, MAYBE still a win.
2. **Use approximation-tolerant comparator** that gives the right
   answer with high probability when inputs are in coset rep. Needs
   construction.
3. **Restructure algorithm** to avoid comparators. B-Y divsteps is
   comparator-free — but we've shown it's worse for Kaliski.

## Path evaluation

Paths in rough order of (Toffoli impact × feasibility):

### Path A: coset rep for non-Kaliski mod-adds
Apply coset ONLY to operations outside Kaliski body (the muls, halves
in pair1_halve/pair2_double). No comparators to worry about. Savings:
- Bulk halving loops: ~207k CCX currently
- Scale corrections / sol_halve_tail: ~22k
- Between-pair muls: ~150k
- Four q×q muls: ~240k
Total target: ~500k CCX. Coset 60% savings → **~300k savings** → ~7% reduction.
Qubit cost: cpad ≈ 26 extra qubits. Easy.

Rough implementation: introduce coset-extended registers for the
non-Kaliski parts; bridge to/from computational basis at Kaliski entry/exit.

### Path B: eliminate m_hist + karatsuba-1 at pair1_mul2
Known unlock from prior sessions. -407q persistent, -28k Toffoli.
Straightforward but has phase-correction complexity.

### Path C: Kaliski per-iter reduction via step 3+9 cswap optimization
Highest-potential Toffoli target but structurally hardest. Requires
inventing a new Kaliski variant. Not in public literature.

### Path D: single-invocation via approximate cleanup
If cleanup obstruction can be relaxed under approximation, ~1.6M
Toffoli saved. Needs novel scheduling. Highest risk, highest reward.

## Recommendation

Given the research:
1. **Path B** (m_hist + karatsuba) is low-risk, small Toffoli impact
   but meaningful qubit impact. Good **first** step.
2. **Path A** (coset for non-Kaliski ops) is medium-risk, moderate
   Toffoli impact (~300k). Good **second** step.
3. **Path C** and **Path D** are high-risk, high-reward. Only attempt
   after A and B have laid foundation.

To hit 2.7M/1300q: need at least B + A + (C or D). Under approximate
correctness, D is more attractive. B + A alone: 4.14M → ~3.8M at
~2300q. Not SOTA.

**This should be the working plan.** No experiments until we decide
which path to pursue.
