# λ — the hidden third axis

## The setup

`apply_tail_nonce` (mod.rs:1714-1726) asserts the last 96 ops are all `X` and rewrites **only `q_target`** on 48
adjacent `X;X` identity pairs. So the circuit FUNCTION is provably identical for all 2^48 nonces. Only the SHAKE256
Fiat–Shamir seed moves, and with it the 9024 test inputs.

That makes the nonce a clean experimental handle: vary it and you resample the test set from the same circuit.

## The current figure — `6909d15`, n=199

Measured 2026-08-04 on `upstream/main` `6909d15`, score `1,486,468,554`, whose `src/` tree is identical to accepted
submission `ed4b529`. 202 trials, each a full build plus a 9,024-shot `eval_circuit` — **no custom screen**
(`04-traps.md` §4). 99 contiguous nonces, 100 at `2^40` stride, 3 controls.

| statistic | classical | phase-garbage batches |
|---|---|---|
| mean | **17.126** | **11.628** |
| sd (sem) | 4.149 (±0.294) | 3.473 (±0.246) |
| var/mean | **1.005** | **1.037** |
| range | 7..29 | 3..22 |
| runs with zero | **0 / 199** | **0 / 199** |

var/mean ≈ 1 on **both** channels is textbook Poisson with zero overdispersion, which independently proves the
per-shot failure probability is identical at every nonce — i.e. the circuit really is nonce-invariant, and this is
intrinsic error, not overfitting. Ancilla garbage is identically 0 by construction (see Traps below).

Poisson-overlap estimator: `classical ~ Pois(λ_c + λ_both)`, `phase ~ Pois(λ_p + λ_both)` with independent
components, so `Cov(c,p) = λ_both` and `λ_total = mean_c + mean_p − Cov`. Pearson ρ(c,p) = **0.569**.

| | λ_total | P(clean) | trials/seed |
|---|---|---|---|
| lower bound, `max(means)` | 17.126 ± 0.294 | 3.7e-8 | 2.7e7 |
| **covariance estimate** | **20.560** (95% CI 18.007–23.016) | **1.2e-9** | **8.5e8** |
| upper bound, `sum(means)` | 28.754 ± 0.384 | 3.3e-13 | 3.1e12 |

$$\lambda_{\text{total}} = 20.560 \quad\Rightarrow\quad P(\text{clean seed}) = e^{-20.560} = 1.2\times10^{-9}$$

Decomposition **λ_classical_only 8.932, λ_both 8.193, λ_phase_only 3.435**. The CI is a 4,000-resample bootstrap over
whole rows (seed 20260804). Direct observation alone bounds `P(clean)` only at < 1.5e-2 — rule of three on 0/199.

**The phase channel alone costs 31×** — `e^(20.560 − 17.126)` — and almost every estimate in circulation quotes
`e^-(classical mean)`.

The estimator assumes the classical-only and phase-only components are mutually independent. If they are positively
correlated beyond the shared term, λ_both is overstated and **λ_total is overestimated**, moving down toward the 17.13
lower bound. The error is one-way, so the covariance estimate stays conservative for planning.

### Method gates — all four passed before the sweep was trusted

1. **Pristine baseline.** A clean `6909d15` tree builds to `ops.bin` md5 `ef30945f3afcb369192ea32897232d2f` at
   **0/0/0**, avgT 1,288,101.386 × 1154 qubits = `1,486,468,554`, matching upstream exactly.
2. **The control nonce is this head's own** — `200321420125`, baked in at `mod.rs:2384` — not `801dd20`'s
   `62000008397024`. All three control rows returned 0/0/0 at the baseline md5.
3. **The knob is live: 199/199 non-control nonces produced 199 distinct md5 values.** This is the check that catches
   `benchmark.sh` preferring `sudo -n bwrap`, where sudo's `env_reset` strips `SUB4_TAIL_NONCE` before
   `build_circuit` sees it and every trial silently re-measures the shipped nonce.
4. **The analysis code was validated first**, by re-deriving every published `801dd20` figure from that sweep's raw
   data before being pointed at this one.

Raw data and full method live in
[`austinamissah/ecdsa-circuit-optimization`](https://github.com/austinamissah/ecdsa-circuit-optimization):
`docs/lambda-6909d15.md`, `docs/data/lambda-sweep-6909d15.tsv`.

## λ has not moved

| head | source | λ_total | classical | phase |
|---|---|---|---|---|
| `02146ca` | the n=700 fit below | 23.29 | 18.13 | 12.64 |
| `801dd20` | `docs/lambda-measurement.md` | 20.04 | 16.231 | 10.915 |
| **`6909d15`** | **this note, above** | **20.560** | **17.126** | **11.628** |

Bootstrap on the difference `6909d15 − 801dd20` gives **+0.525, 95% CI −2.626 to +3.632, with 37% of resamples at or
below zero.** On this evidence the two heads have the same λ_total.

**That stability is the finding.** Eight accepted submissions and a better score moved λ by nothing measurable. Score
work is neither paying for itself in λ nor costing anything in λ — it is not touching λ at all. Anyone planning
against a moving λ, in either direction, is planning against noise.

Both per-channel means did rise marginally — classical +0.894 (t = 2.24, p = 0.025), phase +0.714 (t = 2.12,
p = 0.034) — but λ_both rose in step (7.111 → 8.193), which is why the channel movement does not propagate to the
total. Two uncorrected comparisons at p ≈ 0.03 is weak evidence: worth noting, not worth acting on.

`6909d15` scores **better** than `801dd20` (1,486,468,554 vs 1,487,590,242) while its channels run slightly
**dirtier**. That is the score-versus-λ tension this note describes, showing up in upstream's own progress.

**23.29 is stale.** It was measured on `02146ca` and the figure above supersedes it.

### Clean seeds remain isolated, not clustered

| block | n | classical | phase |
|---|---|---|---|
| A, contiguous `base+1 … base+99` | 99 | 17.354 ± 0.441 | 11.606 ± 0.328 |
| B, `2^40` stride | 100 | 16.900 ± 0.391 | 11.650 ± 0.369 |

Welch **t = +0.77** (classical) and **−0.09** (phase) — indistinguishable. A contiguous neighbourhood of the shipped
clean nonce is no cleaner than nonces scattered across the space, reproducing the `801dd20` isolation result on a
different head with a different base nonce. The closest any trial came to clean was `base+29` at 8 classical / 4
phase. **Grinding near a known-good nonce buys nothing.**

## The generative model — n=700 on `02146ca`

Superseded as a λ figure by the section above, and retained because it is what established that phase-only failures
exist at all.

| statistic | classical | phase-garbage batches |
|---|---|---|
| mean | 18.127 | 12.636 |
| variance | 18.094 | 11.054 |
| var/mean | **0.998** | — |
| range | 8..30 | 4..23 |
| runs with zero | **0 / 700** | **0 / 700** |

Pearson ρ(cm,pg) = 0.5205. Fitting on conditional means `E[pg|cm]` in bins cm=11..23 (20–69 nonces per bin, no
extrapolation) discriminates decisively between two generative models:

- **A** "phase ⊂ classical" (forces pg=0 when cm=0): SSE **13.37**, residuals systematically curved, and cannot reach
  the observed ρ at any parameter (best fit 0.835 vs observed 0.5205).
- **B** "phase-only failures exist": SSE **2.44**. Fitted λ_classical_only 10.05, λ_both 8.08, λ_phase_only 5.16, so
  λ_total 23.29.

Model B is what the current figure refits: 8.932 / 8.193 / 3.435 is the same three-component shape on n=199.

## What this means

The old head computes a **wrong point addition roughly once per 1,100 inversions**. It ships because a lucky seed was
found once and carried forward, with each subsequent submission accepted only if it kept that seed clean.

So the real objective is:

> **minimise score subject to λ small enough to grind**

and λ is exponentially leveraged: every 1.0 removed multiplies grind yield by *e*.

## Where λ comes from (classical channel, modelled to 88%)

Exact classical emulation of the whole ModDiv incl. the Bezout apply, 6e6 samples, per 9024 shots:

| source | mm |
|---|---|
| divstep convergence tail (ITERS=258 vs ~270 needed) | 5.73 |
| i=257 apply skips (ADD_SKIP_LASTK / S2_ZERO / FWD_CSWAP) | 5.30 |
| SCHED_J2 drops a nonzero bit, walk still terminates | 2.80 |
| LSBS=53 fold-window carry escapes | 2.18 |
| **model total** | **16.01** |
| observed (n=700) | 18.13 |

Residual ~2.1 is the square / non-ModDiv point arithmetic. The source decomposition has not been re-run on
`6909d15`, whose classical mean is 17.126; against that the same model total would leave a residual of ~1.1.

ITERS tail curve (1e6-sample convergence distribution), mm per 9024:
`258→5.228, 259→2.453, 260→1.114, 261→0.483, 262→0.200, 265→0.014`. Steep — the first extra iteration is worth a lot
and the seventh is worth nothing. Cost ≈ 2,930 emitted CCX per iteration, dominated by the apply side (256-bit,
width-independent, so it does NOT get cheaper at the tail).

## Traps

- **`ancilla-garbage = 0` is guaranteed by construction, not evidence.** `B::free` (mod.rs:495) emits an unconditional
  `R`; per sim.rs:149-154 an `R` on a non-|0⟩ qubit flips that shot's phase with p=½ and force-zeroes the qubit with the
  outcome DISCARDED. So no qubit can be dirty at the end and that channel cannot fire. Every would-be ancilla failure
  is laundered into half a phase failure.
- **"Every phase failure is a dirty free" is FALSE.** A census-dropped CCZ that no longer cancels gives phase garbage
  on every batch with ZERO dirty resets. Audit the phase word directly.
- **Don't price a truncation site by `2^-w` alone.** MSBS=19 looks like `9024 × 516 × 2^-19 = 8.9` mismatches; measured
  effect of switching the site fully off (w=48) is **zero**. Three factor-of-two discounts: a top-w tie only means the
  low bits decide (½), the correction is gated on `subtracted` (¾), and the block sits inside `push_condition(hmr_bit)`
  (½). It is also an hmr-uncompute feeding a CZ, so it can only ever produce a *phase* error.

## Triage rule (use this constantly)

| full-9024 result | meaning |
|---|---|
| ~9024 classical | positional desync — a sequentially-addressed table shifted |
| thousands but not 9024 | a repointed gate-DROP table |
| low tens (10–30) | **the intrinsic band. Expected. Not a bug.** |
| saturated 141/141 phase, normal classical | bad phase-correction predicate, or a deleted live gate |
| 0/0/0 | you are on a ground seed |

## Statistics discipline

Per-nonce sd is 4.25 (4.15 at n=199 on `6909d15` — unchanged in kind). **n=1 cannot distinguish Δλ=+7 from Δλ=0.**
A reserve retune measured at n=1 as "cm 19, intrinsic, safe" was **+7.24 λ at n=12** (individual draws 19,21,22,23,23,24,25,27,28,29,31,32 — the first two
sit inside the baseline range). Use n≥12, paired on the same nonce set, and quote a sigma.

Also: avg-executed-Toffoli varies across nonces with sd 13.4 (n=700, span 86). So a single-nonce Toffoli comparison
gates at ~40, not 20. This does NOT gate qubit work (1 qubit = 1152 ppm ≈ 2600× the noise) nor deterministic gate
deletion (verify those by gate count).
