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
- `signed_matrix_forward_rows_clean_m_and_match_twos_complement` extends the
  forward-row circuit to a signed sampled matrix
  `[[-8192,24576],[-3,1]]`, using arithmetic right shift after adding `m*p`.
  It matches two's-complement classical rows on 32 random basis states at
  `5563 CCX`, peak `1624q`.
- `adjugate_m_correction_is_integral_for_sampled_by_matrices` proves the
  general cleanup algebra on samples: if `2^w y = P x + p m` and
  `det(P)=s 2^w`, then `s adj(P) y = x + p*(s adj(P)m/2^w)`, and the
  correction vector is integral.
- `qcorr_roundtrip_recovers_m_for_sampled_by_matrices` proves the next
  reversibility hook: with `q=s adj(P)m/2^w`, we have `P q = m`. Thus after
  the old source rows are zeroed, the `m` registers can in principle be
  uncomputed from the small `q` registers, and then `q` can be uncomputed from
  residual high bits.
- `positive_triangular_fixed_matrix_replacement_cleans_old_rows` uses that
  formula to build the first complete fixed-matrix replacement for the
  triangular positive matrix. It computes both scaled rows, recomputes `m` from
  the old sources, zeros the old rows using the noncanonical adjugate residual,
  uncomputes `m` from residual high bits, and uncomputes the residual. It
  simulates correctly on 32 random basis states at `20146 CCX`, peak `1898q`.
- `signed_sample_fixed_matrix_replacement_cleans_old_rows` completes the same
  replacement for the signed matrix `[[-8192,24576],[-3,1]]`. It computes
  signed scaled rows, computes `m`, computes signed `q=s adj(P)m/2^16`, zeros
  the old rows, clears `m` via `Pq=m`, clears `q` from residual high bits, and
  uncomputes residuals. It passes 32 random basis states at `13110 CCX`, peak
  `2224q` after freeing unused q sign-extension bits.
- `fixed_matrix_replacement_sample_cost_distribution` generalizes the circuit
  generator to arbitrary signed sampled BY matrices. On 32 sampled `w=16`
  matrices: mean `20991 CCX`, p90 `24234`, max `28099`, peak `2224q` for the
  full pair replacement.

`actual_matrix_sequence_entropy_supports_sub600_history_target` shows the raw
770-bit selector history is not information-theoretically necessary. Over 10k
sampled secp256k1 denominators, the 35-window matrix sequence has empirical
per-window entropy sum `≈449` bits; an independent per-window entropy code gives
`p99≈463` bits, `p999≈465` bits, and `fail>550=0`. This is not a reversible
circuit, but it says a sub-600-bit matrix-history target is plausible.

`by_tagged_div_stored_matrix_upper_bound_model` separates arithmetic from
selection/history. With per-window matrices already known, one tagged-DIV
window costs an integer denominator replacement plus one modular fixed-matrix
replacement. On 32 sampled windows:

```text
mean/window ≈ 28,607 CCX
p90/window  ≈ 35,087 CCX
max/window  ≈ 37,609 CCX
35 windows  ≈ 1,001,258 CCX
scheduled peak model ≈ 2772q
selector history ≈ 770 bits (35 × 22-bit delta,h key)
```

`branch_bits_reconstruct_by_jump_matrix` proves a simpler exact selector:
each `w=16` BY matrix is reconstructed from the 16 odd/even divstep branch bits
plus starting delta. Thus 35 windows need exactly `560` selector bits, no
large matrix IDs/QROM table. `branch_bit_history_by_tagged_div_budget_model`
combines this with the modular replacement peak: `2224 + 560 + 16 = 2800q`.
This is exact for matrix reconstruction but does not solve how to generate the
branch bits from `x` without a denominator pass.

`h_only_compressed_history_by_tagged_div_budget_model` sketches the next,
more aggressive architecture: delete the full integer denominator pair and keep
only low-ratio state `(delta,h)`, plus compressed matrix history. Using the
measured modular fixed-matrix replacement cost and a conservative 480-bit
history budget from the entropy experiment:

```text
mean modular window ≈ 19,219 CCX
35 windows          ≈ 672,650 CCX
modular peak        ≈ 2224q
history budget      = 480q
h/delta/control     ≈ 32q
modeled peak        ≈ 2736q
```

This is the first BY model simultaneously under 1M Toffoli for the DIV-like
component and under the current 2800q cap. It is not a circuit: it requires a
reversible compressed-history selector and an h-only state update/reverse.

`smith_factorization_reduces_by_window_to_inplace_shifts_and_unimodular_maps`
checked the obvious Smith-normal-form route. It proves the diagonal is always
`diag(1,65536)` for sampled `w=16` BY windows, but the naive SNF factors can be
huge (`~3.9e13`), so plain SNF is not a low-cost in-place implementation.
`hermite_factorization_keeps_scaled_by_window_in_place_with_small_coefficients`
fixes this: for 4096 sampled windows it finds small Hermite factors

```text
U P V = [[1,e],[0,65536]], |e| <= 32768,
max coefficient in U,V,U^-1,V^-1,e <= 65536.
```

Thus a scaled window can, algebraically, be done in-place as:

```text
(x0,x1) -> V^-1(x0,x1)
z0      -> (z0 + e*z1) / 2^16      // one batched Solinas shift
(z0,z1) -> U^-1(z0,z1)
```

This is the concrete route from the current double-buffer row replacement to a
600-scratch implementation: no simultaneous old+new pair, only the two live DIV
registers plus carry/shift/control workspace. `fixed_hermite_inplace_modular_window_matches_scaled_by_matrix`
then builds the first actual fixed-window circuit for the sample signed matrix
`[[-8192,24576],[-3,1]]`. It applies `V^-1`, one row shear by `e=21845`, 16
exact modular halvings, and `U^-1`; 32 random basis states match
`2^-16 P(x0,x1)` exactly. Cost/shape:

```text
sample fixed window: 34,489 CCX, peak 1,285q, factor_ops=10
24-sample distribution: mean 33,715 CCX, p90 43,942, max 44,179
35 windows (naive fixed factors): ≈1,180,034 CCX, peak 1,285q
```

This confirms the scratch breakthrough but also shows that naive Euclidean
shear synthesis is more expensive than double-buffer fixed rows. A better route
is to use the 16 branch bits directly as a numerator microprogram and postpone
the common scaling. `fixed_branch_numerator_window_matches_scaled_by_matrix`
implements the fixed-control circuit:

```text
for each branch bit: apply A/B/C numerator matrix
then halve both rows 16 times
sample window: 18,890 CCX, peak 1,029q
64-sample distribution: mean 22,883 CCX, p90 27,588, max 30,913
35 windows: ≈800,900 CCX, peak 1,029q
```

This is now both lower-Toffoli than the double-buffer fixed rows and has the
right scratch shape. But `quantum_controlled_branch_numerator_replay_is_too_expensive_naively`
shows the control tax: implementing every step with generic quantum-controlled
modular adds costs `77,728 CCX/window`, or `≈2.72M` for 35 windows. Therefore
the remaining SOTA blocker is precise: keep the branch-numerator arithmetic,
but avoid paying generic controlled full-width modular adds for the branch
selection. `low_ratio_microstep_update_generates_branch_bits_without_full_denominator`
shows branch generation itself is small: with `h=g/f mod 2^t`, the next branch
bit is `h&1` and h updates 2-adically by

```text
C: h' = h/2
B: h' = (h+1)/2
A: h' = (h-1)/(2h) mod 2^(t-1)
```

So the selector generator can use only a 16-bit h register plus small delta;
the reversibility payload is the branch history. A sparse-correction variant was checked and mostly killed:
`actual_branch_cases_are_not_sparse_enough_for_a_correction_list` finds actual
560-step secp256k1 branch counts

```text
mean(A,B,C) = (133.5, 133.0, 293.5)
p99_A = 154, p999_A = 162
naive A-position list p99 ≈ 1540 bits
```

So A-cases are not a rare payload; a simple A-position correction list is worse
than raw branch history. `selected_replay_budget_requires_more_than_a_signed_mux`
then quantifies the remaining target using measured primitives:

```text
cmod_add = 1280 CCX, mod_add = 1024, double = halve = 255
naive generic controls        ≈ 2.72M
ideal signed mux + static A   ≈ 1.86M
ideal signed mux + value-A LB ≈ 1.28M
fixed-control lower bound     ≈ 0.80M
```

Thus a signed add/sub mux alone is insufficient if the A-only update is still
paid at all 560 possible positions. A SOTA-grade selected replay needs either
value-proportional/block-specialized A handling near the 1.28M lower bound or a
completely different fixed-control-block mechanism. The target is now numeric:
close the `1.28M -> 0.80M` gap without exceeding ~600 scratch.

`enumerated_branch_block_select_explodes_beyond_single_step` kills the naive
block-SELECT version of that idea. Even ignoring equality-control and QROM
overhead, summing all fixed case-sequence bodies gives lower bounds:

```text
b=1: 3 sequences,  ≈2.576M including scaling
b=2: 8 sequences,  ≈5.725M
b=3: 22 sequences, ≈15.105M
b=4: 58 sequences, ≈38.436M
```

So block specialization cannot mean enumerating all branch case sequences and
SELECTing one. It must exploit algebraic sharing between cases or a new
controlled-add primitive. `signed_mux_controlled_modular_add_works_but_not_enough`
implements the obvious shared primitive for the first odd update:

```text
acc += odd ? (neg ? -a : a) : 0
cost = 1790 CCX, peak 1287q
separate cmod_add+cmod_sub = 2560 CCX
```

It is correct on random basis states and saves ~30% for the A/B first update,
but a full selected replay with this mux and static A still costs `≈2.15M`.
Thus the primitive is useful but not sufficient; the A-only update still needs a
non-static treatment or a deeper algebraic refactor.

That refactor now exists at the microstep level. Instead of numerator replay
plus a final `2^-16` scaling, use the scaled BY step directly:

```text
C: (r,s) -> (r, s/2)
B: (r,s) -> (r, (s+r)/2)
A: (r,s) -> (s, (s-r)/2)
```

For A, controlled-swap `(r,s)` first, then compute `s <- -s + r`, then halve
`s`. This removes the A-only `r += s` correction entirely. The implemented
coherent primitive `scaled_by_controlled_microstep_matches_all_cases_and_hits_target_cost`
uses controls `(odd, A)` and matches all three cases on random basis states:

```text
one scaled controlled microstep = 2046 CCX, peak 1287q
560 steps                       ≈ 1,145,760 CCX
```

`scaled_by_controlled_window_matches_jump_matrix` composes 16 such controlled
microsteps for the sample window and verifies the circuit equals
`2^-16 P(r,s)` for the sampled jump matrix:

```text
16-step controlled window = 32,736 CCX, peak 1,317q
matrix = [[-8192,24576],[-3,1]]
```

`scaled_by_controlled_560_scaffold_cost_model_fits_current_cap` then instantiates
all 560 controlled microsteps with raw `(odd,A)` controls:

```text
560-step scaffold = 1,145,760 CCX, peak 2,405q
raw controls      = 1120 qubits
```

So the full arithmetic scaffold fits the current 2800q cap even before history
compression. It is not yet the user's 600-scratch design, but it is now an
actual costed 560-step circuit skeleton, not only an extrapolated one-step
number. `scaled_by_controlled_560_tagged_div_basis_simulation` goes one step
further: for a sampled converged denominator it sets all 560 controls, starts
`(r,s)=(0,y+x)`, simulates the full 1.145M-CCX circuit on basis states, verifies
bottom channel zero, and recovers `y/x` as `sign(f)*r-1`.

The raw `(odd,A)` controls can be compressed further because A is not an
independent bit. `window_pattern_and_delta_reconstruct_a_controls` proves that a
16-bit odd-pattern plus the starting delta reconstructs all 16 A-controls and
the next delta. Thus the history payload can be branch patterns only; A-controls
are decoder scratch. A pessimistic reversible decoder estimate in
`pattern_decoder_budget_fits_branch_decode_margin` uses a 10-bit signed delta
register and still costs only:

```text
~41 CCX/step, ~22,960 CCX total
```

This fits comfortably inside the 150k branch/decode margin used in the whole
point-add budget.

This is the first coherent selected BY replay model in the right Toffoli band.
It is not yet a complete DIV: branch-history compression/cleanup still need
production handling. The controlled-neg zero representative was explicitly
exercised in the one-step test with `(r,s)=(0,0),(0,s),(r,0)` for all valid
A/B/C controls; the following controlled add canonicalizes the temporary `p`
representative in valid A cases. Algebraically this closes the previous 2.72M
selected-replay blocker without QROM or block SELECT.

`branch_pattern_entropy_supports_compressed_history_target` then checks the
history format needed by this scaled microprogram directly. Instead of storing
raw 560 branch bits or matrix IDs, encode each 16-step window as its branch
pattern. On 10k secp256k1 trajectories:

```text
entropy ≈ 440.2 bits
p99 code length ≈ 458.5 bits
p999 code length ≈ 462.1 bits
fixed per-window distinct-pattern IDs = 481 bits
fail > 520 bits = 0
```

So the branch microprogram itself has a sub-500-bit empirical representation.
Combined with a no-clean-temp / dirty-workspace controlled modular add, this is
the concrete scratch path: the arithmetic should use the history bank as dirty
workspace or avoid the 256-bit AND addend, so peak scratch is history-dominated
rather than `history + adder-temp`.

The first whole-point budget with this primitive is now explicit in
`scaled_by_div_point_add_budget_has_sota_margin_if_history_workspace_solved`:

```text
current total                   = 4,132,750
remove two Kaliski invocations  ≈ -3,190,000
keep non-inversion scaffold     ≈    942,750
scaled BY DIV (2046*560)        ≈  1,145,760
branch/decode margin            ≈    150,000
projected point-add             ≈  2,238,510
```

That is below the 2.7M Google low-qubit target and close to the 2.1M low-gate
target. The remaining work is implementation risk, not arithmetic economics:
complete tagged-DIV integration and solve the scratch overlap/decoder.

This reopens BY as a live SOTA-shaped route but with precise remaining
obstacles: branch/matrix history compression, selected Hermite-factor
application, and integration into a 35-window BY tagged-DIV scaffold. The
fixed-matrix replacement itself is now no longer a one-off; sampled arithmetic
is around 1.0M Toffoli for stored-matrix tagged DIV or ~0.67M for the h-only
compressed-history model, plausibly cheaper than Kaliski but not yet a complete
600-scratch primitive.

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
