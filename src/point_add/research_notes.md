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

#### 3b. Exact matrix-family compression result
Even if entries hit `2^w`, a quantum QROM implementation might still benefit
if the number of **distinct** transition matrices is tiny compared to the raw
state space. I measured this exactly for all low-word states with
`delta ∈ [-20, 20]`, odd `f_low`, and arbitrary `g_low`.

Results:

| w | total states | distinct matrices | compression factor |
|---|---:|---:|---:|
| 4 | 5,248 | 656 | 8× |
| 6 | 83,968 | 2,624 | 32× |
| 8 | 1,343,488 | 10,496 | 128× |

Pattern:
- compression factor = `2^(w−1)` exactly on the observed range.
- equivalently, distinct matrix count appears to scale like `2^(w+2)`.

This does **not** rescue full jumped B-Y by itself, but it is a strong sign
that *compressed local transition classes* are real and exploitable.

#### 3c. Updated verdict on jumped B-Y
Full jumped B-Y still looks too expensive as a drop-in replacement, because:
- matrix entries hit the full `2^w` growth,
- full coefficient tracking would still need to carry those `w`-bit entries,
- cleanup is all-new machinery.

But the compression result changes the local-batching story.

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
- Still potentially interesting as novel research, but now less grounded than
  a compressed Kaliski-local batching route, because we have exact empirical
  class-compression data for the latter.

### 6. Fermat / addition-chain inversion
- Standard classical method; discussed in cryptographic resource estimates.
- Prime-field reversible cost is far too large (hundreds of multiplications).
- Not competitive.

### 7. Itoh–Tsujii
- Only for GF(2^n), not GF(p).
- Not applicable to secp256k1.

## Stronger result: coefficient-side compression matches (u, v) compression

A remaining risk in the hybrid Kaliski-jump idea was that even if the `(u, v)`
window transition family compressed well, the coefficient-side `(r, s)`
transforms might explode and ruin the QROM story.

I derived the per-case coefficient matrices directly from the implemented
`kaliski_iteration` logic:

- UEven: `(r, s) -> (r, 2s)`
- VEven: `(r, s) -> (2r, s)`
- UGtV : `(r, s) -> (r+s, 2s)`
- VGtU : `(r, s) -> (2r, r+s)`

Then I ran the same exact 10,000-input window survey for those coefficient-side
matrices.

**Result:** the `(r, s)` side compresses **identically** to the `(u, v)` side.

| w | t | distinct `(u,v)` mats | distinct `(r,s)` mats | max `|entry|` | mean mats/class |
|---|---:|---:|---:|---:|---:|
| 6 | 4 | 125 | 125 | 16 | 4.506 |
| 8 | 4 | 125 | 125 | 16 | 4.493 |
| 8 | 6 | 1133 | 1133 | 64 | 9.461 |

This removed the biggest remaining objection to the hybrid Kaliski-jump
moonshot.

## Strongest result so far: the **joint** transition family also stays tiny

I pushed the classical analysis one step further and measured the *joint* local
transition object that a reversible batched primitive would actually need to
know: the pair `(uv_mat, rs_mat)`, not just each side separately.

Result on the same 10,000 secp256k1 trajectories:

| w | t | distinct `(u,v)` mats | distinct `(r,s)` mats | distinct joint pairs |
|---|---:|---:|---:|---:|
| 6 | 4 | 125 | 125 | **125** |
| 8 | 4 | 125 | 125 | **125** |
| 8 | 6 | 1133 | 1133 | **1133** |

This is the strongest empirical result in the project so far.

Interpretation:
- The coefficient-side transform is not merely similarly compressible — in the
  sampled data it is effectively **functionally locked** to the `(u, v)` side.
- So a hybrid batched primitive may need only **one compressed lookup** for the
  whole local Kaliski window.

## Strongest result so far, refined again: modest side information collapses ambiguity

The remaining practical question is whether the raw key `(u mod 2^w, v mod 2^w)`
is already enough to select the local transition class, or whether we need extra
metadata (which would cost qubits / logic in the eventual quantum version).

I added `src/point_add/kaliski_jump_extra.rs` and measured how much the branch-
sequence ambiguity drops as we augment the key.

For `w = 8`, `t = 4` on 10,000 real secp256k1 trajectories:

| key | mean sequences/class | max sequences/class | singleton classes |
|---|---:|---:|---:|
| `low = (u_low, v_low)` | 4.492 | 16 | 4,102 |
| `low + cmp0` | 2.570 | 8 | 28,731 |
| `low + cmp0 + cmp1` | 1.742 | 4 | 78,817 |
| `low + cmp0 + cmp1 + low1` | 1.696 | 4 | 163,675 |

Interpretation:
- Just adding the **initial compare bit** nearly halves the ambiguity.
- Adding the **compare bit after the first micro-step** cuts the average class
  ambiguity to ~1.74 and the maximum to 4.
- Even the strongest tested key only gets down to ~1.70 average, so there is
  still some residual ambiguity. But it is *tiny*.

This is a huge deal:
- it suggests a practical hybrid batched primitive does **not** need a full
  branch history or a massive QROM key,
- and that a small amount of dynamically-computed side information may be enough
  to select from a very small family of local transition classes.

## New strongest result: brute-force key search finds a nearly-deterministic 3-bit side channel

I added `src/point_add/kaliski_key_search.rs` and brute-forced feature subsets
of size up to 4 over a reasonable feature family built from:
- compare bits `cmp0, cmp1, cmp2`,
- a few low bits of `(u1, v1)` and `(u2, v2)`.

On 300 random secp256k1 trajectories (~108k 4-step windows), the best key found was:

> **`(u_low, v_low, cmp0, cmp1, cmp2)`**

with statistics:

| key | mean sequences/class | max sequences/class | singleton classes |
|---|---:|---:|---:|
| `(u_low, v_low, cmp0, cmp1, cmp2)` | **1.034** | **2** | **80,193** |

This is dramatically better than the hand-picked `(cmp0, cmp1)` key.
Interpretation:
- three compare bits almost completely determine the 4-step branch sequence,
- the residual ambiguity is only 2-way in the worst observed class,
- and most classes are effectively deterministic.

This is the first evidence that a hybrid Kaliski-jump primitive could plausibly
be driven by:
- low 8 bits of `(u, v)`
- plus just **three compare bits**

instead of a large branch-history register.

## Current best moonshot conclusion

**Conclusion: `hybrid Kaliski-jump is the bet.`**

This is now stronger than before.

### Why full B-Y replacement is not the best bet
Full BY jumpdivsteps2 still has two major problems:
1. matrix entries hit the full `2^w` growth;
2. coefficient tracking and cleanup are all-new machinery.

So a *full* B-Y replacement remains very high-risk.

### Why the new key-search result matters
The exact histogram showed the local transition family is small.
The ambiguity survey showed compare bits collapse the family.
The brute-force key search now shows that **three compare bits are almost
sufficient to identify the exact 4-step branch sequence**.

That means a plausible hybrid primitive could work like:
1. read `(u_low, v_low)`,
2. compute `cmp0, cmp1, cmp2`,
3. lookup 1 of a tiny set of candidate 4-step transforms,
4. if needed, resolve a residual 2-way ambiguity with one extra cheap bit.

That is the most concrete reversible interface I have found so far.

## New classical proposal: hybrid Kaliski-jump

### Model
Standard Kaliski / binary almost-inverse update on `(u, v)` has four branch
cases:

```text
if u even:                   (u, v) ← (u/2, v)
elif v even:                 (u, v) ← (u, v/2)
elif u > v:                  (u, v) ← ((u-v)/2, v)
else:                        (u, v) ← (u, (v-u)/2)
```

Each step is a linear map with a shared `1/2` factor. Over `t` steps we get
an integer 2×2 matrix `P_t` with

```text
(u_t, v_t)^T = (1 / 2^t) · P_t · (u_0, v_0)^T.
```

The classical question is: along actual secp256k1 trajectories, keyed by low
`w` bits of `(u, v)` and a tiny amount of extra branch metadata, how many
possible `P_t` arise?

### Best current empirical lead
For `w = 8`, `t = 4`:
- only **125** joint `(uv, rs)` transition classes globally,
- key `(u_low, v_low, cmp0, cmp1, cmp2)` gives mean ambiguity **1.034**,
- worst observed ambiguity only **2**,
- matrices bounded by `|entry| ≤ 16`.

This is currently the most actionable structural lead toward reducing the 81%
inversion budget.

## Proposed next sessions

### P1. Enumerate the exact 125 four-step joint classes
For `t = 4`, produce:
- canonical representative branch sequences,
- the exact `(uv_mat, rs_mat)` pair,
- the low-bit preconditions / compare-bit conditions under which they occur.

This is the final classical step before a real reversible design sketch.

### P2. Build a concrete reversible selector cost model
Now that the best key is effectively
`(u_low, v_low, cmp0, cmp1, cmp2)` with max ambiguity 2,
we should estimate the true reversible cost of:
- forming the three compare bits,
- indexing the candidate class set,
- resolving the residual 2-way ambiguity.

### P3. Compare `t = 4` vs `t = 6`
`t = 4` has tiny matrices and nearly-deterministic keys.
`t = 6` has larger matrices and larger class ambiguity, but fewer windows.
Current evidence strongly favors `t = 4` as the first real prototype.

## Bottom line

The strongest current research judgement is:

> The best moonshot is **hybrid Kaliski-jump batching** over 4-step windows,
> keyed by low bits plus three compare bits, because the exact local transition
> family is tiny on both the state side `(u, v_w)` and the coefficient side
> `(r, s)`, and the branch-sequence ambiguity is almost gone.

That's still novel research, but unlike the other moonshots, it now has
clear empirical support directly tied to the 81%-of-budget hot path.
