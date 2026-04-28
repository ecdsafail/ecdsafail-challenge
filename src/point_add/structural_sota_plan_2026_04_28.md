# Structural SOTA plan — 2026-04-28

User directive: stop local tuning. This file is the working structural model for
matching Google's secp256k1 point-add frontier, not a micro-optimization list.

## 0. Current measured point

Current committed baseline after the last in-flight retune:

- **4,132,750 Toffoli**
- **2716 qubits**
- exact / phase-clean under the 24-seed gate and checks

Google/ZKP targets:

- low-qubit: **2.7M Toffoli @ 1175 qubits**
- low-gate: **2.1M Toffoli @ 1425 qubits**

So the real gaps are:

- **−1.43M Toffoli** to low-qubit
- **−2.03M Toffoli** to low-gate
- **−1541 qubits** to low-qubit

## 1. The Toffoli gap is one inversion-sized object

Measured decomposition of the current design:

| component | cost |
|---|---:|
| Kaliski invocation #1 (`with_kal_inv_raw`, fwd+body+bwd scale excluded) | ~1.60M |
| Kaliski invocation #2 | ~1.59M |
| non-Kaliski scaffold (muls, scale correction, Solinas, constants) | ~0.94M |
| total | ~4.13M |

Therefore a SOTA-grade design must do one of exactly two things:

1. **Delete one full inversion invocation.** Current primitives then land at
   roughly `4.13M - 1.60M = 2.53M`, already matching Google's low-qubit
   Toffoli target.
2. **Keep two inversions but make each ~45% cheaper.** Need per-invocation
   cost `~1.60M -> ~0.90M`, i.e. save ~700k per Kaliski invocation.

Anything that does not attack one of these two objects cannot match SOTA.
A 50k improvement is useful only if it is a stepping stone toward one of
these two structural changes.

## 2. Why the one-inversion route keeps failing

The map is in-place:

```text
(Px, Py) -> (Rx, Ry) = (Px, Py) + (Qx, Qy)
```

A one-inversion affine formula exists classically (Strategy C, `w = dx^3`),
but reversible cleanup is the hard part. The issue is not the algebraic
formula for `Rx,Ry`; it is zeroing the old input information.

The obstruction in every attempted one-inversion schedule is:

- To uncompute the inversion input (`dx`, `dx^3`, or a product containing it),
  the circuit needs a live source of `dx` after `tx` has become `Rx`.
- Reconstructing `dx = Px-Qx` from `(Rx,Ry,Qx,Qy)` is exactly point subtraction
  by the classical point `Q`.
- Point subtraction has the same denominator `(Rx-Qx)` that the current second
  Kaliski invocation inverts.

So the naive one-inversion design just moves the second inversion into the
cleanup path.

### Strategy C re-estimate at the current baseline

Classically correct formula:

```text
dx = Px-Qx, dy = Py-Qy, w=dx^3
v  = dy^2 - dx^2(Px+Qx)
Rx = v * dx * w^-1
Ry = (dy(dx^2 Qx - v) - w Qy) * w^-1
```

Cost at current 407/403 Kaliski settings:

| block | estimated Toffoli |
|---|---:|
| compute `dx^2, w=dx^3` | 250-300k |
| one Kaliski invocation on `w` | 1.60M |
| compute `dy^2, v, Rx, Ry` with Bennett-clean temps | 1.7-2.0M |
| uncompute `w, dx^2` | 250-300k |
| misc/scale | 100-200k |
| **total** | **3.9-4.4M** |

This is not a SOTA route unless the cleanup can be made triangular/in-place
instead of Bennett-clean. That triangular schedule is the unresolved research
problem.

## 3. Therefore the most credible Toffoli route is a cheaper Kaliski

Current Kaliski cost sources, across both forward and backward and both
invocations:

| substructure | total cost | SOTA relevance |
|---|---:|---|
| step 3 + step 9 cswaps | ~1.0M+ | biggest public lever |
| step 4 controlled `v-=u; s+=r` | ~1.0M | second biggest lever |
| comparator / eqzero / flags | ~0.4M | moderate |
| scale correction loops | ~0.2M | not enough alone |

To reach 2.7M with two inversions, we need to remove/replace **most of the
cswap + step4 cost**. This points to a jumped/windowed/divstep-style Kaliski,
not to local adder swaps.

## 4. Candidate structural programs

### Program A — jumped/windowed Kaliski (highest Toffoli relevance)

Batch `t` binary-GCD microsteps into one matrix update:

```text
[u']   1/2^t [a b] [u]
[v'] =       [c d] [v]

[r']          [A B] [r]
[s'] =        [C D] [s]
```

For `t=2..4`, coefficients are small. The hope is that one matrix-selected
update costs less than `t` copies of:

- two full cswap layers, and
- one step4 sub/add layer.

Target economics:

| window | current cost for t steps | needed jumped cost | result |
|---:|---:|---:|---|
| t=2 | ~2× current step | <1.2× current step | ~40% Kaliski win |
| t=3 | ~3× | <1.8× | ~40% win |
| t=4 | ~4× | <2.4× | ~40% win |

This is exactly the size of win required for SOTA with two inversions.

**Fast invalidation criterion:** if reversible matrix application needs more
than ~2 q-q adds/subs plus one controlled shift per microstep, it cannot beat
current Kaliski. If a low-coefficient t=2 or t=3 schedule can be synthesized
with one scratch n-register and <=~1.5n Toffoli/step equivalent, it is live.

Fresh survey from `kaliski_jump.rs` / scratch replay (10k inputs for the
low-bit table, 2k inputs for quick t sweep):

| t | distinct matrices | max coeff | mean log2 coeff | mean mixed rows |
|---:|---:|---:|---:|---:|
| 2 | 13 | 4 | 1.79 | 0.85 |
| 3 | 41 | 8 | 2.58 | 1.12 |
| 4 | 125 | 16 | 3.29 | 1.34 |
| 6 | 1133 | 64 | 4.71 | 1.63 |

Low-bit lookup is **not** enough: at `w=8,t=4`, each low-bit class still has
`mean 4.49` possible matrices and up to `16`; at `w=8,t=6`, mean `9.46`, max
`62`. A strengthened executable invalidation
`initial_gt_window_classifier_not_approx_good_enough` adds one full comparator
bit (`u>v`) to the low-bit key and still sees a disjoint-sample majority error
of about **60%** for `w=8,t=4`. So a one-comparator window is not even close to
1% approximate correctness.

A real jumped Kaliski must either compute the whole comparator sequence
coherently, use a Bernstein-Yang/divstep variable that avoids full comparisons,
or synthesize a matrix application whose cost beats the per-step loop despite
those predicates.

Positive qubit-side result: `window_hint_bits_can_compress_history_but_not_select_matrix_alone`
records the actual matrix choice as a small per-window hint instead of per-step
history. On 5k sampled trajectories with key `(low8(u), low8(v), u>v)`:

| t | max matrices/key | hint bits/window | total hint bits |
|---:|---:|---:|---:|
| 4  | 8  | 3 | 306 |
| 8  | 23 | 5 | 255 |
| 16 | 34 | 6 | 156 |

So window hints can plausibly save 100-250 history qubits versus `m_hist`, but
only if a selected matrix can be applied cheaper than replaying microsteps. This
is qubit-structural, not yet a Toffoli route.

`selected_matrix_application_arithmetic_intensity_model` measures a simple
row-popcount add/sub model for selected matrices. It ignores QROM, multiplexing,
and reversible cleanup, so real cost is higher:

| t | mean matrix row-add terms | mean raw odd-step add/sub count | max |
|---:|---:|---:|---:|
| 4  | 5.30 | 3.99 | 14 |
| 8  | 13.89 | 7.97 | 44 |
| 16 | 34.94 | 15.73 | 74 |

This means selected-matrix windowing cannot win by reducing arithmetic row
terms; it must win by deleting many cswaps/comparators/control scaffolds. That
focuses the synthesis target sharply: a viable implementation needs a coherent
matrix application that avoids generic controlled-cswap replay and does not pay
QROM/control overhead proportional to all candidate matrices.

Another executable check, `global_window_matrix_indices_do_not_compress_history`,
separates the hint idea from lookup cost. If we store a **global** matrix id
instead of a low-state-keyed short hint, sampled distinct matrices explode:

| t | observed global matrices | global id bits/window | total bits |
|---:|---:|---:|---:|
| 4  | 125    | 7  | 714 |
| 8  | 9478   | 14 | 714 |
| 16 | 111696 | 17 | 442 |

So the qubit compression requires a low-state-keyed QROM/table. It is not a
free history encoding. The next synthesis must include that QROM/control cost.

Next concrete work: synthesize/lower-bound selected matrix application for
`t=4..16` with QROM/control costs included. If it cannot exploit cswap deletion
strongly enough, move to BY/divstep or a different DIV transform.

### BY/divstep jump update (new stronger candidate)

The Bernstein-Yang route deserves renewed attention because branch selection is
local to `(delta, low f, low g)` rather than full `u>v` comparisons. New tests
in `by.rs` add two relevant facts:

1. `jumpdivstep_matrix_arithmetic_intensity_model` row-popcount model for one
   full-width pair under the exact 742-step bound:

| w | mean row-add terms/window | exact windows | mean terms/pair |
|---:|---:|---:|---:|
| 4  | 2.04 | 186 | 379 |
| 8  | 4.51 | 93  | 419 |
| 12 | 7.66 | 62  | 475 |
| 16 | 11.56| 47  | 543 |

2. With approximate tolerance, `approximate_divstep_cutoff_survey` on 20k
   samples gives `q99=549`, `q999=555`, `fail>550≈0.0062`, `fail>560≈0.0001`.
   So a 550-step approximate BY inversion is within the user's 1% failure
   allowance empirically, reducing `w=16` windows from 47 to 35.

`jumpdivstep_budget_model_suggests_live_prototype` turns that into an optimistic
lower-bound budget for applying each selected matrix to three full-width pairs
`(f,g)` plus the two coefficient columns, charging one n-bit add/sub per
row-popcount term:

```text
w=16 exact 742-step bound: 47 windows, ≈416,782 Toffoli lower bound
w=16 approx 550-step cap : 35 windows, ≈310,370 Toffoli lower bound
```

This ignores matrix synthesis, sign handling, reversible cleanup, and modular
normalization, so it is not a forecast. But it is far below the current ~1.6M
Kaliski invocation cost. BY jump inversion is therefore the most concrete live
prototype candidate now.

A first circuit-level calibration, `constant_jump_matrix_apply_cost_probe`,
applies sampled constant `w=16` BY matrices to one full-width pair using the
real add/sub primitives (row formation only, not full reversible update):

```text
mean_ccx      ≈ 3,908 per 274-bit pair
mean_terms    ≈ 12.58
ccx/term      ≈ 310.6
row peak      ≈ 1370q for f,g,out0,out1 + carries
```

Scaling naively to three pairs and 35-47 windows gives roughly 0.4-0.55M
Toffoli for row formation before matrix synthesis/cleanup. That is still far
below current Kaliski's ~1.6M/invocation, so a live BY prototype is justified.
The peak number also shows why register scheduling matters: row formation wants
four full-width registers plus one carry strip; doing three pairs in parallel is
not viable, but sequential coefficient-column updates may fit.

`jumpdivstep_full_state_budget_model` combines the row former with a sequential
six-register BY state model:

```text
width              = 274 bits (256 + w + sign/slack)
state              = (f,g) + two coefficient columns = 6 wide regs
shared row outputs = 2 wide regs
carry strip        = 1 wide reg
modeled peak       ≈ 2514q
exact row cost     ≈ 534k Toffoli
approx row cost    ≈ 398k Toffoli
```

This is the first BY model that simultaneously fits the current 2800q cap and
has a row-formation cost far below current Kaliski. The missing pieces are now
concrete: reversible low-word matrix synthesis, row-output cleanup/swap, sign
normalization, and modular reduction/recovery of the inverse.

This is not yet a full inversion circuit, but it is a better Toffoli-structural
lead than Kaliski low-bit windows: no full comparator sequence, moderate matrix
row intensity, and approximate iteration count is plausible.

#### BY correction after deeper circuit modeling

The next round made the BY picture more precise:

- `fixed_by_coeff_channel_is_tagged_div_when_converged` proves the same
  `y+x` tagged-DIV algebra works for fixed-cap BY. After `K` divsteps, if
  `f=±1,g=0`, then `V*x = sign(f) 2^K` and `R=0 mod p`, so carrying
  `y+x` gives `sign(f)*V*(y+x)*2^-K - 1 = y/x`. At `K=550`, sampled failure
  was `29/5000 = 0.0058`, within the user's 1% allowance.
- `jump_matrix_depends_on_delta_and_g_over_f_ratio` shows the selected matrix
  is determined by `(delta, h=g/f mod 2^w)`, not by both low words. Exact
  enumeration gives `41*2^w` keys for `w=4,6,8`, matching the histogram law.
  For `w=16` this is a 22-bit key (`~2.7M` matrices), not a 33-bit key.
- `scaled_pair_update_cleanup_cost_probe` measures the integer denominator
  jumped replacement with scaled-adjugate cleanup: `≈7744 CCX/window/pair`,
  peak `≈1402q`.
- But the modular coefficient/tag channel is harder than the integer
  denominator. `modular_jump_inverse_cleanup_is_dense_dead_end` shows that
  unscaled modular inverse cleanup uses `2^-w adj(P) mod p`, whose four
  constants have mean popcount `≈814`; this kills naive sparse cleanup.
- `naive_variable_coefficient_jump_apply_is_too_expensive` shows synthesizing
  quantum coefficient bits and applying all possible bits would cost
  `≈5.2M` Toffoli for the 2-pair 35-window tagged DIV alone.
- `by_microstep_inplace_cost_model_is_not_the_jump_win` measures raw coherent
  BY microsteps at `≈5989 CCX/step`, i.e. `≈3.29M` for 550 steps.
- `hybrid_jump_denominator_with_microstep_tag_channel_still_too_costly` tries
  the valid hybrid (jumped integer denominator + raw modular tag channel) and
  gets `≈2.66M` for one tagged DIV.
- `scaled_modular_jump_sparse_cleanup_is_too_expensive_with_current_primitives`
  tries the scaled coefficient convention (`2^-w P` forward + sparse adjugate
  cleanup) with a shared-doubling small-constant modular row former. It still
  costs `≈58.4k CCX/window`, or `≈2.05M` for 35 windows for the modular pair.

The very next lead is the batched `2^-16` shift. For a canonical row value
`T`, choose `m=-T*p^{-1} mod 2^16`, add `m*p`, and shift right by 16. Because
`p=2^256-(2^32+977)`, adding `m*p` is sparse: add `m` at bit 256 and subtract
`m*(2^32+977)` at bits `{0,4,6,7,8,9,32}`. The correction `m` is recovered
from the top 16 output bits except for the negligible set `T < m*(2^32+977)`.

New tests:

- `batched_halve16_top_bits_recover_correction_with_negligible_exception`:
  `0/20000` sampled failures; explicit rare exception `T=1` has `m=13617`,
  top bits `13616`. This is an approximate primitive with failure probability
  around `2^48/p`, far below 1%.
- `highfold_then_batched_halve16_matches_row_distribution`: for sampled BY
  row values `T=a*x+b*y`, first folding `k=T>>256` copies of `p` and then the
  batched halve had `0/40000` failures.
- `approximate_batched_shift_reopens_scaled_by_jump_budget`: high-fold cost
  `≈1862 CCX`, batched shift cost `≈1915 CCX`; integer row+cleanup floor
  `≈6976 CCX`; scaled modular pair window `≈18254 CCX` after also high-folding
  the two old-row cleanup residuals; 35 windows `≈639k` for the modular pair.

`approximate_batched_halve16_canonical_circuit_matches_classical` then builds
and simulates the actual canonical batched-shift circuit on 64 random basis
states, matching the classical `(T+m*p)/2^16` result. Finally,
`windowed_scaled_by_tagged_division_matches_microstep_algebra` validates the
full classical `w=16`, 35-window scaled BY tagged-DIV algebra: `0/3000`
failures at 560 steps, bottom channel zero, and output `sign(f)*r-1 = y/x`.

Caveat: `noncanonical_batched_shift_needs_quotient_uncompute` shows the
highfold quotient is not recoverable from the scaled output alone: `T` and
`T+p` produce the same scaled residue but different low-word corrections. A
real reversible row primitive must therefore keep the quotient, recover it from
row sources, or fuse row reduction with cleanup. The canonical batched shift is
real; the noncanonical row highfold is still an integration problem.

Another tempting branch-history compression was tested and mostly killed:
`low_ratio_window_state_needs_large_rank_history` keeps only
`h=g/f mod 2^16` plus `delta` to select matrices. The window map
`(delta,h)->(delta',h')` is many-to-one; on sampled secp256k1 35-window
trajectories, reversing it needed rank up to `71769` (`17` bits/window), and a
16-bit/window rank would fail about `10.95%` of inversions. So low-ratio-only
state is not the 600-scratch escape.

Carry-slack correction: the earlier shifted-row cost helper truncated carries
when summing multi-bit coefficients. After fixing it to extend addends to the
remaining row width, the 3-pair full BY cleanup model becomes `≈1.03M` Toffoli
but `≈2852q`, just over the current cap. The 2-pair optimistic integer-cleanup
lower bound is `≈575k` at 35 windows but `≈2304q` (`≈1792q` beyond two field
registers), so it is not a 600-scratch primitive by itself.

Positive forward-row progress:

- `noncanonical_scaled_pair_map_is_injective_on_canonical_domain` shows the
  two-row scaled matrix map can be injective on canonical input pairs even
  though one row alone loses representative quotient. This keeps fixed-matrix
  pair replacement algebraically possible.
- `fixed_positive_matrix_forward_rows_clean_m_and_match_classical` builds and
  simulates the first actual noncanonical forward row circuit for the positive
  matrix `[[65536,0],[65535,1]]`: correction `m` is computed from the original
  sources and uncomputed from those same sources after the shift. It matches
  the classical rows on 32 random basis states at `8772 CCX`, peak `1624q` for
  forward rows only.
- `adjugate_m_correction_is_integral_for_sampled_by_matrices` proves the
  general cleanup algebra on samples: if `2^w y = P x + p m` and
  `det(P)=s 2^w`, then `s adj(P) y = x + p*(s adj(P)m/2^w)`, and the
  correction vector is integral.
- `positive_triangular_fixed_matrix_replacement_cleans_old_rows` uses that
  formula to build the first complete fixed-matrix replacement for the
  triangular positive matrix. It computes both scaled rows, recomputes `m` from
  the old sources, zeros the old rows using the noncanonical adjugate residual,
  uncomputes `m` from residual high bits, and uncomputes the residual. It
  simulates correctly on 32 random basis states at `20146 CCX`, peak `1898q`.

This reopens BY as a live SOTA-shaped route but with precise remaining
obstacles: signed/general matrix replacement and branch/matrix history
compression. The next implementation target is to generalize the fixed-matrix
replacement from the triangular positive matrix to arbitrary signed BY matrices;
then attempt a BY tagged-DIV scaffold.

### Program B — triangular one-inversion schedule (highest payoff, highest risk)

Goal: use Strategy C or B2 but avoid Bennett-clean fresh outputs. A successful
schedule must satisfy:

1. Kaliski input (`dx` or `dx^3`) remains uncomputable after output mutation.
2. `tx,ty` are transformed in-place, not via fresh `(rx,ry)` registers.
3. Any copied slope/inverse is phase-uncomputed from live state without
   inverting `(Rx-Qx)`.

Fast invalidation criterion: if the schedule ever contains both
`old dx` and `new Rx` as independent live n-bit values after Kaliski backward,
it has already lost; zeroing one from the other is point-subtraction.

### Program C — Kaliski representation rewrite for qubits

Even a cheaper Kaliski will not hit 1175q with the current state layout.
Required qubit reductions:

| source | potential saving |
|---|---:|
| `m_hist` compression/elimination | 407q |
| fold Kaliski input copy into `tx`/scratch | 256q |
| fold `r` output into next multiplier/output scratch | 256q |
| low-workspace step4 / venting only at peak sites | 200-260q |
| register sharing / length-tracked tail | 100-300q |

The easy `m_i` start-state formula is **not sufficient** because iter-end
state does not expose that fingerprint. A real qubit breakthrough needs either:

- a self-cleaning Kaliski body whose inverse branch is recoverable from end
  state, or
- pebbling/checkpoint recomputation, or
- Luo-style length/location registers.

This is qubit-structural, but by itself does not solve Toffoli.

### Program D — coset/padded arithmetic only after a long region exists

Exact padded/coset add chains cross over at ~12 repeated additions and save
~44% by 256 additions, but cost +500-800 qubits in the current prototype.
Current affine scaffold has no long non-Kaliski add region. Coset becomes
relevant only if Program A creates windowed/batched arithmetic regions or if
we accept a larger representation rewrite.

## 5. Decision rule going forward

Do not pursue a code change unless it plausibly satisfies at least one:

- **Toffoli structural:** can save >=0.5M by deleting an inversion or cutting
  Kaliski per-step cost.
- **Qubit structural:** can save >=400q without >0.5M Toffoli regression.
- **Fast falsification:** conclusively kills a tempting structural path so we
  do not waste another session on it.

Immediate next target: **Program A**, because it is the only public-ish path
whose economics can plausibly produce the missing 1.4M Toffoli while staying
inside the exact harness.
