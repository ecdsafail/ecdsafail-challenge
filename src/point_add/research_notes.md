# Research notes — inversion moonshots inside `src/point_add/`

Session: 2026-04-22 (continued, moonshot-only work).

This file keeps all moonshot literature / classical-analysis work under
`src/point_add/`, per the current scope rules.

## Deliverable 1 (classical B-Y on secp256k1) — confirmed

Implemented classical `divstep2` reference and modular-inverse recovery in
`src/point_add/by.rs`, then ran a 10,000-input secp256k1 survey.

Results:

| metric | value |
|---|---|
| theoretical bound `⌈(49·256 + 57)/17⌉` | 742 |
| observed minimum iters | 502 |
| observed maximum iters | 567 |
| observed mean iters | 531.01 |
| max `|δ|` observed | 20 |
| modinv matches (vs Fermat) | 10,000 / 10,000 |

Interpretation:
- The BY safegcd upper bound is pessimistic by ~24% on secp256k1 inputs.
- However, this is **not enough** to save plain B-Y: the per-iter reversible
  cost is still too high relative to Kaliski.

## Deliverable 2 (algorithm-space survey) — corrected final version

### 1. Kaliski almost-inverse (baseline)
- Classical ref: Burton S. Kaliski Jr., “The Montgomery inverse and its
  applications,” IEEE Trans. Computers 44(8), 1995.
- Quantum / reversible refs:
  - Roetteler–Naehrig–Svore–Lauter 2017, arXiv:1706.06752.
  - Häner–Roetteler–Soeken 2020, arXiv:2001.09580 / ePrint 2020/077.
- Iterations in our tuned circuit: 399.
- Measured per-iter reversible cost: ~2180 CCX.
- Per-pass cost: ~1.81M CCX.

### 2. Bernstein–Yang divstep2 (w = 1)
- Ref: Bernstein–Yang 2019, ePrint 2019/266.
- Reversible implementation: unpublished / would be novel.
- Empirical iterations on secp256k1: max 567, mean 531.
- Per-iter reversible estimate: 10–12n CCX.
- Conclusion: still worse than Kaliski.

### 3. Bernstein–Yang jumpdivsteps2 (w > 1)
- Ref: Bernstein–Yang 2019, Figure 10.2 / §10.
- Reversible implementation: unpublished / would be novel.

#### 3a. Corrected matrix-growth result
A previous version of the jump survey undercounted the scaled transition
matrix. After fixing it, the 100,000-sample survey now shows the **full
scaled** transition matrices do hit the theoretical `2^w` growth.

Corrected survey over 100,000 random low-word states:

| w | max observed `|entry|` | max log2 | mean log2 | theoretical max log2 |
|---|---:|---:|---:|---:|
| 4  | 16    | 4.00  | 2.03 | 4  |
| 8  | 256   | 8.00  | 4.28 | 8  |
| 12 | 4096  | 12.00 | 6.34 | 12 |
| 16 | 65536 | 16.00 | 8.19 | 16 |

Interpretation:
- The **maximum** entry size really does hit the full `2^w` growth.
- So a faithful reversible matrix-apply must still handle `w`-bit classical
  coefficients.
- That restores the pessimistic reversible cost model: batching by `w` does
  not automatically beat Kaliski.

#### 3b. New follow-up: matrix *histogram* / compression structure
Even if entries hit `2^w`, a quantum QROM implementation might still benefit
if the number of **distinct** jump matrices is tiny compared to the total
state space. I measured this exactly for all low-word states with
`delta ∈ [-20, 20]`, odd `f_low`, and arbitrary `g_low`.

State-space and matrix-count results:

| w | total states | distinct matrices | compression factor |
|---|---:|---:|---:|
| 4 | 5,248 | 656 | 8× |
| 6 | 83,968 | 2,624 | 32× |
| 8 | 1,343,488 | 10,496 | 128× |

Strong pattern:
- distinct matrices scale as `2^(2w+1) / 2^(w-1) = 2^(w+2)`
- equivalently, **compression factor = 2^(w−1)**

Most common matrix pattern observed:
- w=4:  `[[-4,  2], [ -2, -3]]`
- w=6:  `[[-16, 2], [ -8, -3]]`
- w=8:  `[[-64, 2], [-32, -3]]`

This means:
- a naive QROM over all `(delta, f_low, g_low)` states is wasteful;
- but a compressed QROM over matrix classes is possible.

However, the compressed class count is still:
- w=8  → 10,496 entries
- w=12 → likely ~167,936 entries
- w=16 → likely ~2,686,976 entries

This is *much* better than the raw state space, but still large enough that
we would need a serious select-swap/QROM design to exploit it.

#### 3c. Updated verdict on jumped B-Y
The moonshot is not completely dead, but it is no longer a literature-backed
"obvious next step." The corrected state is:
- matrix coefficients are still worst-case `w`-bit,
- but the matrix family collapses by a factor `2^(w−1)`.

So jumped B-Y is now a **QROM-compression problem**. If we can exploit the
small number of distinct matrices *and* make matrix-apply cheap enough, it
might still work. This is now clearly novel research territory.

### 4. Montgomery inverse (Savaş–Koç)
- Classical ref: Savaş–Koç 2000, “The Montgomery modular inverse revisited.”
- Quantum / reversible refs: effectively same family as RNSL/HRSL Kaliski.
- Conclusion: not a distinct win over Kaliski in our setting.

### 5. Lehmer-style GCDs
- Classical refs: Lehmer 1938; Jebelean 1993.
- Reversible implementation: unpublished / novel.
- Main issue: runtime matrix selection depends on quantum data, so a faithful
  reversible implementation needs a QROM keyed by top bits. No concrete,
  literature-backed reversible cost win established yet.
- Still potentially interesting as novel research, but much less grounded than
  a compressed jumped-BY route, because we now have actual empirical structure
  for B-Y matrices and none for Lehmer.

### 6. Fermat / addition-chain inversion
- Standard classical method; discussed in cryptographic resource estimates.
- Prime-field reversible cost is far too large (hundreds of multiplications).
- Not competitive.

### 7. Itoh–Tsujii
- Only for GF(2^n), not GF(p).
- Not applicable to secp256k1.

## Deliverable 3 — final conclusion (updated after histogram analysis)

**Conclusion: `hybrid Kaliski-jump is the bet.`**

This is a refinement of the earlier “no known path remains” conclusion.
Here's why.

### Why full replacement B-Y is still not the best bet
Full BY jumpdivsteps2 still has two major problems:
1. matrix entries hit the full `2^w` growth;
2. coefficient tracking and cleanup are all-new machinery.

So a *full* B-Y replacement remains very high-risk.

### Why the histogram result changes the hybrid story
The jump histogram shows that there are far fewer distinct local transition
matrices than raw low-word states. This suggests a new hybrid direction:

> keep Kaliski's global state machine and cleanup structure,
> but replace some local per-iteration parity/cswap/sub/halve stretches
> with **pre-batched micro-transitions** chosen from a compressed matrix/QROM.

This avoids the hardest part of BY (full coefficient-system replacement)
while still exploiting the empirical low-word structure.

### The actual research bet now
The best remaining moonshot is:

## **Hybrid Kaliski-jump**

Specifically:
- batch 4–8 Kaliski micro-iterations at a time,
- use a compressed lookup over the small family of observed local transition
  classes,
- keep the existing `(r, s, m_hist)` cleanup logic as much as possible,
- only replace the expensive `(u, v_w)` update path.

### Why this beats the alternatives
- Better grounded than full B-Y replacement.
- Better grounded than Lehmer (which still lacks any empirical structure here).
- Avoids Montgomery / Jacobian cleanup obstruction.
- Targets the actual hot path: Kaliski is ~81% of total cost.

## Proposals for future sessions

### Proposal P1: enumerate Kaliski 4-step local transition classes
Take the current Kaliski update on `(u, v_w)` and enumerate the exact state
transition induced by 4 low-bit steps, keyed by the same style of low-word
state used in BY. Measure:
- number of distinct transition classes,
- coefficient growth,
- whether a compressed QROM could represent them cheaply.

This is the next classical moonshot task.

### Proposal P2: compressed QROM design study
Given the jump-matrix class count, derive a reversible cost model for:
- raw-state QROM,
- matrix-class QROM,
- select-swap QROM,
- unary-encoded small-matrix QROM.

### Proposal P3: full BY jump prototype only if P1 fails
If Kaliski-local batching has no class compression, then fall back to the
full jumped-BY moonshot.

## Bottom line

After correcting the jumpdivstep survey and adding the exact histogram study,
my best current research judgement is:

> The best moonshot is no longer “replace Kaliski with B-Y.”
> The best moonshot is **compress and batch Kaliski using B-Y-style local
> transition classes**.

That's novel research, but it is the most focused path that still directly
attacks the 81% inversion budget.
