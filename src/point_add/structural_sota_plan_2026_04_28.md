# Structural SOTA plan — 2026-04-28

User directive: stop local tuning. This file is the working structural model for
matching Google's secp256k1 point-add frontier, not a micro-optimization list.

## 0. Current measured point

Current committed exact baseline:

- **4,089,274 Toffoli**
- **2713 qubits**
- **24,292,787.153 Clifford**, **30,575,944 emitted ops**
- exact / phase-clean under the hardened alt-seed gate and checks

Google/ZKP targets:

- low-qubit: **2.7M Toffoli @ 1175 qubits**
- low-gate: **2.1M Toffoli @ 1425 qubits**

So the real gaps are:

- **−1.389M Toffoli** to low-qubit
- **−1.989M Toffoli** to low-gate
- **−1538 qubits** to low-qubit

## 1. The Toffoli gap is one inversion-sized object

Measured decomposition of the current design:

| component | cost |
|---|---:|
| Kaliski invocation #1 (`with_kal_inv_raw`, fwd+body+bwd scale excluded) | ~1.60M |
| Kaliski invocation #2 | ~1.59M |
| non-Kaliski scaffold (muls, scale correction, Solinas, constants) | ~0.94M |
| total | ~4.11M |

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

A tempting Montgomery-trick/batched-inversion variant is also now closed.  The
two denominators are `dx=Px-Qx` and `bx=Rx-Qx`.  The chord cubic gives the exact
identity

```text
dx * bx = 3 Qx^2 - 2 λ Qy
```

(`chord_product_identity_does_not_batch_the_two_affine_inversions` verifies this
on 200 secp samples).  This product is easy once `λ` is live, but it does not
remove the second inverse: `1/bx = dx / (dx*bx)`.  Replacing `inv(bx)` by
`inv(dx*bx)` merely moves the second Kaliski to a different denominator and adds
at least one variable multiply (`149,889` CCX floor) to recover `1/bx`.

A slope-carried coordinate channel is also closed for the windowed setting.
Representing a point as `(x,λ_Q)` with `λ_Q=(y-Q_y)/(x-Q_x)` would avoid erasing
one slope for a single fixed `Q`, but changing the next selected table point to
`Q'` requires

```text
λ_Q' = (Q_y - Q'_y + λ_Q (x-Q_x))/(x-Q'_x)
```

which is exactly another variable division.  The executable probe
`slope_carried_coordinate_retargeting_is_dense_division` gets full-domain toy
`n=10` degree `19/20`, density `522204/1048576`, and curve-support `n=12`
minimum degree `4`.  A single slope tag therefore does not provide a universal
2n-bit accumulator representation for Google's windowed point-add circuit.

The user's scratch budget tightens the coefficient-DIV story: only ~600 qubits
are available above the two 256-bit input coordinates.  A folded one-pair
Kaliski DIV already spends `u+r=512` scratch, leaving 88 bits for branch cleanup.
`scratch600_sidecar_tag_bits_do_not_fix_kaliski_branch_recovery` gives a small
independent coefficient sidecar the best possible chance and still finds exact
toy branch recovery needs `n-1` sidecar bits (`n=8 -> 7`, with 4 bits leaving
9500 conflicts).  Extrapolated to secp256k1, the sidecar wants ~255 bits, not
88.  The companion test
`unreduced_coefficient_kaliski_self_cleans_but_width_kills_scratch600` shows the
other side of the tradeoff: if coefficient registers are not reduced modulo `p`,
poststate branch recovery is exact on exhaustive toys, but coefficient width
grows as `3n-1` (toy `n=8 -> 23`; secp ≈663-bit coefficients), implying about
1326 scratch in the folded layout.  Thus a 600-scratch design cannot be
"current Kaliski + small tag" or "unreduced Kaliski".  Measuring away the high
quotient bits is not a cheap kickmix fix either:
`measuring_unreduced_coefficient_high_bits_has_dense_phase` gets toy n=10 degree
`20/20` and density `516314/1048576` for a representative high bit in the input
frame.  The viable transform must make the inverse branch local without
field-sized quotient side information in the first place.

New curve-support check: restricting the actual image set to toy curve inputs
makes poststate collisions rare (`curve_restricted_tagged_kaliski_poststate_ambiguity_is_small_but_not_exact`:
ambiguous occurrence fractions `1.12%, 0.886%, 0.578%, 0.476%` for
`n=8,10,12,14`; exact sidecar bits `3,5,6,6`).  That is an information-theoretic
opening for approximate no-history Kaliski **only if** the reverse step can test
curve support cheaply.  It cannot do so from the Kaliski poststate alone:
`secp_curve_support_does_not_make_kaliski_branch_choice_locally_free` runs on
actual secp256k1 curve-supported samples and still sees local predecessor
ambiguity in `23016/26048 = 88.36%` of steps.  The obvious discriminator is the
shifted curve equation for candidate predecessor `(dx,dy)`, i.e. a cubic field
check per Kaliski iteration, which is not SOTA-shaped.  Do not count the rare
curve-support collision rate as a free branch oracle.  The simple compact-tag
variant is also too weak for the 88-bit slack if scaled linearly:
`implementable_curve_sidecar_still_extrapolates_over_88q_slack` keeps an actually
updatable two-lane sidecar modulo `2^b` and finds best toy lane widths
`2,3,5,6,6` for `n=8,10,12,14,16`; pair bits are already 12 at n=16, linearly
extrapolating to 192 bits at secp256k1.  A sublinear/entropy-coded sidecar is
still a possible research object, but the naive two-lane 2-adic tag is not a
600-scratch solution.  The natural one-lane linear/eigenvector compression is
also closed: `one_lane_linear_sidecar_has_no_closed_update_for_all_kaliski_branches`
finds zero common one-dimensional linear sidecar updates modulo `2^1..2^8` for
all four reachable branch matrices.  Sparse rank coding has the same issue in a
different form: `curve_collision_event_oracle_is_dense_in_natural_input_frame`
shows the natural full-domain predicate for whether an input trajectory hits a
curve-supported ambiguous poststate is already degree `19/20`, density
`487716/1048576` at toy `n=10`.  Redundant centered coefficient representatives
also hit the same width wall: `redundant_centered_coefficients_still_need_growing_range_on_curve_support`
needs extra range bits `6,9,11` for `n=8,10,12`, and `n=14` still has 301
conflicts at extra `11`.

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

The obvious branchless alternative, Fermat inversion by a fixed addition chain,
is not in the same cost class with current multiplication primitives.  The
executable floor `fermat_fixed_chain_inversion_floor_misses_sota_by_order`
measures a modular-square layer at `148,342` CCX; any chain reaching an exponent
near `2^256` needs at least 255 square/multiply layers, so even an optimistic
single inverse costs at least `37,827,210` CCX before scratch and Bennett
cleanup.  Avoiding Euclidean branch history by exponentiation is therefore not
a SOTA route.

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

Follow-up hard gate: `selected_matrix_variable_coeff_lower_bound_kills_hybrid_kaliski_windows`
adds the missing selected-coefficient lower bound.  Even omitting QROM/equality
controls, signs, output cleanup, and old-register cleanup, variable-coefficient
row formation for the two Kaliski pairs costs:

```text
t=4:  window_ccx=41,040   windows=102  invocation_lower=4,186,080
t=8:  window_ccx=73,872   windows=51   invocation_lower=3,767,472
t=16: window_ccx=139,536  windows=26   invocation_lower=3,627,936
```

This is above the entire current Kaliski invocation, not near the needed
~0.9M/invocation.  Therefore generic selected-matrix Kaliski windows are dead
with current arithmetic.  Any Kaliski-window revival must avoid quantum variable
coefficients entirely (e.g. fixed-control microprograms or a fundamentally
shared algebra), not just QROM-select a matrix.

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
independent bit. `scaled_by_pattern_history_560_tagged_div_scaffold_reduces_peak`
replaces the 1120 raw control qubits by 560 raw odd-pattern bits plus one
16-bit A scratch window and still simulates tagged-DIV correctly:

```text
raw-pattern 560-step scaffold = 1,145,760 CCX, peak 1,861q
```

`window_pattern_and_delta_reconstruct_a_controls` proves that a 16-bit
odd-pattern plus the starting delta reconstructs all 16 A-controls and the next
delta. Thus the history payload can be branch patterns only; A-controls are
decoder scratch. The decoder is now an actual reversible circuit, not only a
budget line: `reversible_pattern_delta_decoder_matches_and_cleans` implements

```text
(pattern, delta_start, A=0) -> (pattern, delta_next, A_mask)
```

and its inverse. It matches sampled patterns/deltas, restores `delta_start`,
cleans A, and is phase clean:

```text
one 16-step window decoder: 1,776 CCX forward, 3,552 CCX roundtrip, peak 53q
35 windows forward-only:    ≈62,160 CCX
```

This is inside the 150k branch/decode margin used in the whole point-add budget.
`scaled_by_pattern_decoder_560_tagged_div_scaffold_is_clean` wires the decoder
around the full 560-step tagged-DIV replay: expand all A controls from raw
patterns, run the replay, then reverse the decoders to clean A and restore
`delta=1`. It is phase-clean and exact:

```text
decode forward       =    62,160 CCX
replay + decode      = 1,207,920 CCX
roundtrip clean      = 1,270,080 CCX
peak                 = 2,415q
```

This all-A-history schedule is not the final low-scratch implementation, but it
turns the pattern-history replay into a clean reversible circuit. A tempting
window-local schedule was tested in
`window_local_a_clear_fails_phase_with_mbu_microsteps`: keep only 16 A bits plus
350 boundary-delta bits, decode one window, use A immediately, then clear A
before later windows. It reduces peak but fails the phase contract with the
current measurement-based modular microsteps:

```text
window-local A clear: 1,270,080 CCX, peak 2,221q, phase=1
```

Classical tagged-DIV data is correct, but clearing controls early invalidates
phase corrections left by the MBU modular add/halve primitives. Exact modular
microsteps were calibrated in
`exact_scaled_microstep_is_phase_safe_but_too_expensive_for_window_local_clear`:

```text
fast MBU scaled microstep   = 2,046 CCX
exact scaled microstep      = 4,350 CCX
exact 560-step replay       ≈ 2,436,000 CCX
```

So the obvious exact/phase-safe replay would restore window-local control
clearing but lose the SOTA Toffoli shape. A more surgical phase lead now exists:
`live_reduction_flag_microstep_hits_replay_target_but_needs_cleanup` keeps the
modular-add reduction flag live instead of immediately uncomputing it with the
MBU comparator:

```text
live-reduction-flag microstep = 1,790 CCX
560-step replay body          ≈ 1,002,400 CCX
```

`live_reduction_flags_make_window_local_a_clear_phase_safe_candidate` then wraps
this in the window-local A decoder schedule. It is data-correct and phase-clean
while clearing each 16-bit A window immediately:

```text
live flags + window-local A clear: 1,126,720 CCX, peak 2,780q, phase=0
```

The remaining garbage is 560 modular-add reduction flags. The direct recovery
relation is known: `live_reduction_flag_is_recoverable_from_doubled_output_but_cleanup_is_costly`
checks that, away from the fast-negation zero representative edge, the flag is
`odd && (2*out_s mod p) < r_out`. But direct cleanup would need a doubled-output
copy plus a comparator and uncompute, about `766 CCX/flag`, which is more than
the `cmp_lt` cost we skipped.

The flags are also not a sparse side list. `live_reduction_flag_history_is_dense_and_high_entropy`
measures actual tagged-DIV trajectories:

```text
mean true flags              ≈ 133.8 / 560
p99 true flags               = 155
independent entropy upper    ≈ 436.1 bits
```

So a position-list escape is dead; at best they are another compressed history
bank.

A stronger representation rewrite has now been validated algebraically:
`redundant_signed_scaled_by_replay_avoids_reduction_flags_algebraically` skips
modular reduction before every halve. For any signed integer representative `T`,

```text
T/2 mod p  is represented by  (T + (T&1)*p)/2
```

This removes the red/reduction flag entirely. On 2,000 secp256k1 tagged-DIV
samples the redundant signed replay is exact modulo `p`, has zero convergence
failures, and representatives stayed small:

```text
max observed representative magnitude <= 2p
parity_mean                         ≈ 276.5 / 560
parity entropy upper                ≈ 559 bits
```

`redundant_signed_microstep_is_cheap_if_parity_history_can_be_cleaned` builds a
260-bit signed circuit for one step, leaving only the pre-halve parity bit:

```text
redundant signed live-parity microstep = 1,297 CCX
560-step replay body                   ≈ 726,320 CCX
peak                                   = 1,043q
```

A better centered representative discipline also works:
`centered_signed_redundant_replay_stays_within_half_modulus` chooses `+p` or
`-p` on odd pre-halve values according to the sign, keeping the whole tagged
channel in a narrow signed range:

```text
max observed magnitude bits = 255  (< p)
parity_mean                 ≈ 264.9 / 560
```

The centered circuit version has now also been synthesized. A single
centered signed microstep costs `1,560 CCX`, and a full 560-step tagged-DIV
scaffold with raw odd/A controls and dirty parity history is exact and phase
clean:

```text
centered signed live-parity microstep = 1,560 CCX
560-step centered scaffold            = 873,600 CCX
peak                                  = 2,723q
```

The parity bits are dense, but they are not information-theoretically
irrecoverable. `centered_parity_is_recoverable_from_poststate_range_for_add_cases`
shows in a small exact model that, with centered inputs, parity is recovered by
testing whether the even preimage is centered:

```text
B: parity = !(2*s_out - r_out is centered)
A: parity = !(r_out - 2*s_out is centered)
C: parity = !(2*s_out is centered)
```

The inverse/product-clean direction has also been synthesized with the same
cost:

```text
centered signed inverse 560 scaffold = 873,600 CCX
peak                                 = 2,722q
```

The centered replay body is now wired into the benchmark path behind
`BY_CENTERED_REPLAY_BODY_BENCH=1` as a clean no-op smoke hook (one zero odd/A/
parity control reused to avoid adding the unsolved history bank). The real
harness accepts it:

```text
BY_CENTERED_REPLAY_BODY_BENCH=1
avg_toffoli = 4,985,518
qubits      = 2,716
emitted_ops = 38,874,541
altseed/classical/phase/ancilla failures = 0
```

Default path remains unchanged at `4,111,918 Toffoli @ 2,716q`. A stronger
nonzero forward+inverse centered roundtrip hook was initially classically and
ancilla clean but phase dirty. The root cause was not the signed MBU adders: it
was the centered **unhalve** cleanup. The inverse correction chose add/sub
controls from the doubled value's sign, but tried to uncompute them from the
post-correction sign; when parity=1 that sign flips. Keeping a one-qubit
`sign_hist` through the correction fixes the dirty controls and their R-phase.

After the unhalve sign-history fix, the phase isolation changed to:

```text
96-step centered forward+inverse+parity-clear variants
fast MBU controls:          299,616 CCX, phase=0
exact signed add/sub only:  399,264 CCX, phase=0
exact parity ±p controls:   448,800 CCX, phase=0
all exact controls:         548,448 CCX, phase=0

560-step exact-parity/fast-signed: saw only phase 0 across 12 samples
fixed-trace all-exact hook unit:   phase=0, qubits=2,464
```

The requested all-exact fallback is now wired into the real benchmark path as
`BY_CENTERED_CLEAN_ROUNDTRIP_BENCH=1`. It carries a fixed real BY control trace
(raw odd/A plus parity), runs 560 all-exact centered steps, reverses them, and
cleans parity from restored rows. Because the hook changes the circuit hash, it
uses the earlier safe Kaliski bulk prefix `370` under this env flag only. The
real harness accepts it:

```text
BY_CENTERED_CLEAN_ROUNDTRIP_BENCH=1
avg_toffoli = 7,311,738
qubits      = 2,976
emitted_ops = 42,092,380
altseed/classical/phase/ancilla failures = 0
```

This is intentionally not SOTA-shaped, but it proves the clean centered
forward+inverse/parity-clear object can live in the benchmark harness. More
importantly, the sign-history fix reopens the fast MBU centered cleanup path:
the previous data-dependent phase diagnosis was polluted by a dirty unhalve
control, not necessarily by a fundamental signed-adder phase witness.

The fast cleaned body is now also wired behind
`BY_CENTERED_FAST_CLEAN_ROUNDTRIP_BENCH=1`. It uses the same fixed real BY
control trace and corrected unhalve, but uses the fast MBU centered replay body
instead of all-exact arithmetic:

```text
BY_CENTERED_FAST_CLEAN_ROUNDTRIP_BENCH=1
avg_toffoli = 5,860,218
qubits      = 3,235
emitted_ops = 47,024,860
altseed/classical/phase/ancilla failures = 0
```

This is still Kaliski + appended BY roundtrip, not a replacement. The useful
increment is `≈1.748M` Toffoli for a clean forward+inverse centered replay with
raw odd/A/parity controls, matching `2×873,600` plus tiny parity-clear overhead.
The price is peak width (`3,235q`) because raw controls and fast-adder workspace
are all live. It is nevertheless the first production-harness proof that the
SOTA-shaped fast centered replay can be cleaned once the correct unhalve
sign-history is used.

A first denominator-derived integration hook is now wired as
`BY_CENTERED_DENOM_CONTROLS_BENCH=1`. This still appends a no-op after the
normal Kaliski point-add, but it no longer uses a fixed odd/A trace: it copies
the live quantum output-x register into a 560-bit 2-adic denominator state,
reversibly generates BY odd/A controls, runs the clean fast centered replay
roundtrip on scratch, reverses the denominator generator to clear controls, and
uncopies the denominator. The real harness accepts it:

```text
BY_CENTERED_DENOM_CONTROLS_BENCH=1
avg_toffoli = 7,868,378
qubits      = 4,964
emitted_ops = 55,824,052
altseed/classical/phase/ancilla failures = 0
```

This is far too wide and expensive, but it is the first production-harness
component where BY controls are genuinely denominator-derived from live quantum
data and then self-cleaned.

`BY_CENTERED_LIVE_NUM_BENCH=1` takes one more wiring step: it uses the same
quantum-denominator-derived controls, but the centered replay scratch is now a
centered copy of the live quantum `y` output register rather than zero. The
forward+inverse replay is still a no-op, and the centered copy is uncentered and
uncopied afterwards:

```text
BY_CENTERED_LIVE_NUM_BENCH=1
avg_toffoli = 7,870,438
qubits      = 4,965
emitted_ops = 55,834,807
altseed/classical/phase/ancilla failures = 0
```

So the production harness now accepts live denominator-derived controls acting
on live numerator-derived signed data, with all controls/parity/numerator copy
cleaned.

The first actual replacement hook is now wired as
`BY_CENTERED_PAIR1_REPLACE=1`. It replaces pair1's Kaliski inverse and first
schoolbook multiplication by a centered BY tagged-division Bennett computation:

1. copy `dx` into a full-width signed BY denominator state,
2. generate 576 BY odd/A controls (full width, not tapered, so final `sign(f)`
   is available),
3. build a tagged numerator copy `dy+dx`, center it, and run the fast centered
   replay forward,
4. copy `lam=-dy/dx` out of the tagged frame using `lam = 1 - sign(f)*r`,
5. inverse replay to clean parity and numerator scratch,
6. reverse denominator generation to clean controls,
7. use the existing pair1 `mul2` to zero `ty`.

The real point-add harness accepts it:

```text
BY_CENTERED_PAIR1_REPLACE=1
avg_toffoli = 8,172,700
qubits      = 5,589
emitted_ops = 51,801,485
altseed/classical/phase/ancilla failures = 0
```

This is much worse than baseline because it uses the direct full-width 576-step
denominator generator and Bennett replay, and it still keeps pair1 `mul2` plus
all of pair2 Kaliski. But it is the first env-gated path where a Kaliski-sized
piece is actually replaced by a BY tagged-division computation inside the point
add. The 560-step version had one main-harness classical miss; 576 steps passed
all 9,024 checked shots.

`BY_CENTERED_PAIR2_REPLACE=1` now replaces the pair2 Kaliski cleanup as well.
It computes `-(ty/tx)` with the same 576-step centered BY tagged-DIV machinery
and adds it directly into the live `lam` accumulator, avoiding an uncomputed
quotient output register. Pair2's `ty` product is preserved, then the ordinary
`ty -= Qy` cleanup finishes the affine `y` output:

```text
BY_CENTERED_PAIR2_REPLACE=1
avg_toffoli = 8,189,908
qubits      = 5,589
emitted_ops = 51,930,771
altseed/classical/phase/ancilla failures = 0
```

A follow-up changed the quotient-add/copy into the live accumulator from exact
controlled modular add/sub to the fast MBU versions; this stayed phase-clean and
slightly reduced the replacement probes:

```text
BY_CENTERED_PAIR1_REPLACE=1
avg_toffoli = 8,171,417

BY_CENTERED_PAIR2_REPLACE=1
avg_toffoli = 8,187,601
```

Both replacements can be enabled together and pass:

```text
BY_CENTERED_PAIR1_REPLACE=1 BY_CENTERED_PAIR2_REPLACE=1
avg_toffoli = 12,246,560
qubits      = 5,589
emitted_ops = 72,997,701
altseed/classical/phase/ancilla failures = 0
```

Trying to also make the signed-to-modular copy's sign-controlled ±c correction
MBU-fast failed with one phase-garbage batch, so that copy remains exact for now.

`BY_CENTERED_WINDOW_DENOM_REPLACE=1` is the first production-harness wiring of
the lowword/window denominator-control source. For each 16-step window it copies
low signed `(f,g)` plus `delta` into a 34-bit local 2-adic simulator, generates
odd/A controls, applies those selected controls to the full-width signed
denominator state so `sign(f)` is still available, and later reruns the same
window oracle after reversing the full state to clear the histories. This still
uses 576 full-width selected microsteps, so it is deliberately worse than the
direct generator, but it proves the selected/window interface works in the real
quotient-producing replacements:

```text
BY_CENTERED_PAIR1_REPLACE=1 BY_CENTERED_WINDOW_DENOM_REPLACE=1
avg_toffoli = 8,790,041
qubits      = 5,589
emitted_ops = 55,077,495
altseed/classical/phase/ancilla failures = 0

BY_CENTERED_PAIR2_REPLACE=1 BY_CENTERED_WINDOW_DENOM_REPLACE=1
avg_toffoli = 8,806,225
qubits      = 5,589
emitted_ops = 55,195,212
altseed/classical/phase/ancilla failures = 0

BY_CENTERED_PAIR1_REPLACE=1 BY_CENTERED_PAIR2_REPLACE=1 BY_CENTERED_WINDOW_DENOM_REPLACE=1
avg_toffoli = 13,483,808
qubits      = 5,589
emitted_ops = 79,499,589
altseed/classical/phase/ancilla failures = 0
```

`BY_CENTERED_WINDOW_Q_DENOM_REPLACE=1` extends that hook with the actual
`pattern+q` payload: after the local 16-step signed simulator runs, it xors the
34-bit signed `(q0,q1)` rows into raw per-window histories. The histories live
across the real quotient replay/copy and are cleared by rerunning the same
lowword oracle during denominator reverse. The q payload is not consumed yet,
but it is now threaded through both real replacement sites:

```text
BY_CENTERED_PAIR1_REPLACE=1 BY_CENTERED_WINDOW_Q_DENOM_REPLACE=1
avg_toffoli = 8,790,041
qubits      = 8,037
emitted_ops = 55,084,839
altseed/classical/phase/ancilla failures = 0

BY_CENTERED_PAIR2_REPLACE=1 BY_CENTERED_WINDOW_Q_DENOM_REPLACE=1
avg_toffoli = 8,806,225
qubits      = 8,037
emitted_ops = 55,202,556
altseed/classical/phase/ancilla failures = 0

BY_CENTERED_PAIR1_REPLACE=1 BY_CENTERED_PAIR2_REPLACE=1 BY_CENTERED_WINDOW_Q_DENOM_REPLACE=1
avg_toffoli = 13,483,808
qubits      = 8,037
emitted_ops = 79,514,277
altseed/classical/phase/ancilla failures = 0
```

The Toffoli count is unchanged from the control-only window hook because q
copyout is CNOT-only. The +2,448 qubit peak is exactly the raw
`2 × 36 windows × 34 qbits` history bank; final architecture must consume q
window-locally instead of persisting it.

Next deletion: replace the 16 full-width selected-control microsteps with the
`pattern+q` fixed-matrix/window update. The new window hooks prove control
selection and q-payload threading compose; they do not yet consume the q
corrections.

`BY_SCALED_PAIR2_PRODUCT_REPLACE=1` is the ground-up pair2 variant that should
have been tried earlier: do not compute `ty=lam*tx` and then divide it back out.
Instead, use scaled BY inverse/product-clean directly. Generate denominator
controls, frame the live slope as `u=-sign(f)*lam`, and run inverse scaled BY so
`(u,0)->(0,lam*tx)`. This avoids centered parity history, deletes the ordinary
`mul3_between_pair`, and cleans `lam` as the product is written:

```text
BY_SCALED_PAIR2_PRODUCT_REPLACE=1
avg_toffoli = 8,038,619
qubits      = 4,236
emitted_ops = 49,073,481
altseed/classical/phase/ancilla failures = 0
```

This is significantly cleaner than the centered pair2 quotient-add replacement
(`8,187,601 @ 5,589q`) and validates the inverse/product-clean architecture in
the real affine schedule. It is still far worse than default because it still
uses direct 576-step full-width denominator compute/uncompute. This confirms the
hard boundary: **BY only has a SOTA chance if denominator generation becomes a
selected/window/consumed update; product-clean replay itself is viable.**

The first ground-up consumed-denominator idea was invalidated early in
`consumed_denominator_window_branchless_recovery_is_exponentially_ambiguous`.
If we try to consume the BY denominator and later recover/clear a 16-step
branch window from the post-window denominator alone, the inverse relation is
massively many-to-one. The tiny poststate `(delta,f,g)=(0,1,0)` has:

```text
depth 16 predecessor states/patterns = 589,824
```

So a branchless poststate cleanup cannot be the SOTA denominator solution. A
surviving BY design must either carry compressed branch/pattern history, or
consume pattern+q inside a reversible fixed-matrix window update with an
explicit local inverse. This invalidates the tempting “no history, recover from
consumed denominator” route before more point-add wiring.

The remaining constructive fixed-matrix/q window route was then costed as a
single-window object in `last_shot_fixed_matrix_window_consumption_misses_sota_budget`.
It uses the actual reversible fixed-matrix circuit we know how to build: form the
scaled rows, clean the old rows, and clear `m/q/z` locally. This is the object
that would replace 16 full-width denominator microsteps if it were cheap enough.
The measured distribution is:

```text
mean window = 20,323.5 CCX
p90 window  = 24,382 CCX
max window  = 27,641 CCX
max peak    = 2,224q
two 576-step denominators ≈ 1,463,292 CCX
SOTA-shaped target/window ≈ 10,000 CCX
```

This misses the required per-window target by about 2× before replay/scaffold
costs. Therefore the fixed-matrix/q window implementation with current arithmetic
primitives is also not SOTA-shaped. A BY SOTA path would now require a genuinely
new arithmetic primitive for selected window updates, not just better scheduling
of the existing fixed-matrix circuit.

Again, this is not a performance result. It proves both Kaliski inverse-sized
objects in the real affine point-add can be replaced by BY tagged-DIV/product
objects and still satisfy the full exact checker.

Naively synthesizing the range test is too expensive:
`naive_centered_parity_recovery_cost_would_erase_redundant_replay_win` measures
about `1,296 CCX/flag`, or `≈725,760` CCX just to clean 560 parity bits. A cheap
high-bit-only approximation is not acceptable either:
`centered_parity_highbits_recovery_is_too_approximate_without_boundary_fix`
finds `65,709/1,120,000 = 5.87%` mismatches, above the user's 1% tolerance.
It is cheaper (`≈524 CCX/B-flag`) but needs an exact boundary correction.

With the current centered replay body and selected/tapered denominator target,
the folded-cleanup budget is SOTA-shaped:

```text
non-BY scaffold after deletions      ≈   285,766
2× centered signed replay            ≈ 1,747,200
2× selected/tapered denominator      ≈   607,656
projected if parity folded           ≈ 2,640,622
naive parity cleanup projection      ≈ 4,092,142
```

So the parity problem is now precise: **fold centered-range recovery into the
inverse/window arithmetic instead of doing it as a separate post-hoc cleanup**.
The moonshot replay problem has therefore split into two concrete options:
clean/absorb 560 red flags at ~1.00M replay, or move to centered signed
representatives and solve centered-range parity cleanup around an `873,600` CCX
replay body.
A production replay therefore needs either (a) keep A controls/live flags until
the end and clean them globally, or (b) fuse the modular average so the same
flag is recoverable from later state.
`clean_two_replay_by_budget_requires_replay_or_phase_breakthrough` updates the
whole-point economics with the clean decoder:

```text
non-BY scaffold after deleting old phases ≈   285,766 CCX
2× forward decoded replay                ≈ 2,701,606 CCX
2× clean decoded replay                  ≈ 2,825,926 CCX
2× fixed-control replay + 300k branch    ≈ 2,187,566 CCX
```

So merely plumbing denominator generation into the current clean decoded replay
is not enough; it is already slightly over 2.7M before branch generation. The
moonshot target is now sharper: get the selected replay closer to the fixed-
control ~0.8M/replay band, or solve early A clearing without paying exact-
arithmetic costs. `quantum_branch_values_do_not_reduce_replay_toffoli_accounting`
confirms there is no hidden average-case relief from sparse branch values: the
benchmark counts a quantum-controlled CCX as issued on every live shot, so the
fast replay's average Toffoli per shot equals the static `1,145,760` count.

`compressed_pattern_history_scratch_model_is_600q_if_add_workspace_is_removed`
spells out the remaining scratch equation:

```text
current local microstep workspace over (r,s): 775q
compressed pattern history:                  481q
one-window A scratch + delta:                 26q
current additive scratch:                   1282q
if no-clean-temp/dirty controlled add:       ~597q
```

So the 600-scratch target is not compatible with the current clean-temp
`cmod_add_qq` implementation, but it is compatible with the branch-pattern
history if the controlled modular add uses no clean 256-bit addend (or borrows
the history bank as dirty workspace).

The existing venting module already contains the right substrate for this:
`dirty_quantum_offset_adder_is_plausible_cmod_add_substrate` measures
`iadd_dirty_2clean_qoffset`:

```text
quantum offset add: 762 CCX
workspace: 2 clean + 254 dirty, no hidden clean n-register
```

A first controlled version, `ciadd_dirty_3clean_qoffset`, was added and checked
on small basis states. It has the right scratch shape but is too Toffoli-heavy
when built by the naive controlled-qoffset transformation:

```text
n=8 basis check: passes, dirty restored, phase clean
n=256 controlled qoffset: 3557 CCX, peak 770q
scaled BY step with this substrate ≈ 4323 CCX
560-step DIV ≈ 2.42M just for replay
```

So dirty qoffset addition solves the scratch direction but not yet the Toffoli
direction. A newer mask-and-borrow variant gets the opposite tradeoff:
`masked_controlled_qoffset_borrows_offset_as_dirty_gate_good_scratch_short`
computes a clean mask `m_i = ctrl & offset_i`, uses the original offset row as
the dirty workspace for the uncontrolled qoffset adder, then uncomputes the mask.
It is phase-clean and gate-shaped:

```text
n=256 masked-borrow controlled qoffset: 1274 CCX, peak 771q
scaled BY replay with this add:          ≈1,142,400 CCX
compressed history + mask scratch:       ≈766q
```

Thus control-efficient dirty qoffset is possible, but the clean mask misses the
user's ~600-scratch cap by ~166q when added to compressed pattern history. A
linear partial-mask interpolation looked slightly short of the target:
`partial_mask_controlled_qoffset_linear_tradeoff_just_misses_600q_target` leaves
90 mask bits under the 600q cap, estimates add≈2755 CCX, replay560≈1,971,760,
and point-add≈2,764,476 with scaffold+branch margin.  The real one-bit streaming
implementation is better than that model once simple broadcasts are emitted as
direct controlled toggles instead of materializing the mask:
`streamed_mask_controlled_qoffset_fits_scratch_and_hits_lowqubit_target` keeps
one mask bit, passes basis/phase checks, fits scratch (`scratch_with_history≈510q`),
and measures `2542` CCX with `div560≈1,852,480`, projecting point-add≈2,645,196.
So the controlled-add scratch/gate blocker is no longer the main issue for the
low-qubit BY route.  Early invalidation still prevents an integration attempt:
`streamed_mask_qoffset_still_has_no_selector_margin_for_integration` leaves only
87,840 CCX for denominator-derived pattern generation after the measured decoder,
but the cheap lowword pattern oracle is already 208,320 CCX for one denominator
(projected point-add≈2,765,676) and the known exact tapered generator projects
≈4,565,516.  Branch-pattern production/cleanup must be solved below that ~88k
margin before any affine hook-up is justified.

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
current total                   = 4,111,918
remove two Kaliski invocations  ≈ -3,169,168
keep non-inversion scaffold     ≈    942,750
scaled BY DIV (2046*560)        ≈  1,145,760
branch/decode margin            ≈    150,000
projected point-add             ≈  2,238,510
```

`low_scratch_scaled_by_budget_still_beats_27m_after_pair1_mul_deletion` adds an
important correction: true tagged DIV also deletes pair1's two schoolbook
multiplications (`149,889 + 150,145 CCX`) because it directly maps
`(0,dy+dx)->(lambda+1,0)`. With that saving:

```text
scaffold after DIV              ≈    642,716
fast BY projected               ≈  1,938,476
low-scratch vented BY projected ≈  2,650,796
```

Thus even the higher-Toffoli low-scratch/vented variant can beat the 2.7M
Google low-qubit target, while the fast variant is below 2.1M on paper. The
remaining work is implementation risk, not arithmetic economics: complete
tagged-DIV integration and solve the scratch overlap/decoder.

`inverse_scaled_by_560_cleans_lam_and_writes_product` resolves the earlier
pair2-cleanup conceptual gap. The forward scaled BY map sends
`(0, q*x) -> (sign*q, 0)`, so the inverse map sends `(sign*q,0) -> (logical 0,
q*x)`. This is exactly the pair2 operation needed after `tx=Rx-Qx` and
`lam=q=-lambda`: it cleans `lam` and writes `Ry+Qy=q*tx` into `ty`, after which
one classical subtract by `Qy` gives the output. Circuit evidence:

```text
inverse 560-step scaled BY product-clean scaffold: 1,287,440 CCX, peak 2,403q
input  (r,s)=(sign*q,0)
output (r,s)=(logical zero, q*x)
```

A better sign-frame removes the inverse replay's extra subtract cost.
`inverse_scaled_by_560_negr_frame_recovers_fast_cost` stores `u=-r` during the
inverse product-clean pass. The inverse cases become

```text
C: (u,s) -> (u, 2s)
B: (u,s) -> (u, 2s+u)
A: (u,s) -> (u+2s, -u)
```

so the implementation uses `cmod_add` instead of `cmod_sub` and matches the
forward replay cost:

```text
neg-r inverse product-clean: 1,145,760 CCX, peak 2,405q
```

The zero is sometimes the noncanonical representative `p` because of the fast
controlled negation; for the non-exceptional affine path it can be unloaded as a
known logical-zero representative or fixed with a canonical negation variant.
This means the full replacement is two BY replays: forward tagged DIV for pair1
and neg-r inverse product-clean BY for pair2, deleting both Kaliski objects and
the ordinary multiplication cleanup around them. With measured decoder margins,
`low_scratch_scaled_by_budget_still_beats_27m_after_pair1_mul_deletion` now
estimates the two-fast-replay schedule at `≈2,637,286` Toffoli, below the 2.7M
low-qubit target.

This reopens BY as a live SOTA-shaped route but with precise remaining
obstacles: branch/matrix history compression, selected Hermite-factor
application, and integration into a 35-window BY tagged-DIV scaffold. The
fixed-matrix replacement itself is now no longer a one-off; sampled arithmetic
is around 1.0M Toffoli for stored-matrix tagged DIV or ~0.67M for the h-only
compressed-history model, plausibly cheaper than Kaliski but not yet a complete
600-scratch primitive.

A direct quantum branch generator is now explicitly ruled out as the first
integration. `two_adic_branch_generator_matches_classical_prefix_on_small_width`
proves the natural 2-adic generator is algebraically correct: keep `f,g` modulo
`2^W`, emit `odd=g0`, compute `A=(delta>0)&odd`, and update the denominator with
the same swap/neg/add/halve skeleton as the scaled numerator replay. But
`naive_quantum_branch_generator_would_erase_scaled_by_savings` shows why this
must not be wired into the real point-add path:

```text
W=560 direct branch generator step      ≈ 2,880 CCX (optimistic)
one generator compute+uncompute         ≈ 3,225,600 CCX
two generators                           ≈ 6,451,200 CCX
projected point-add with naive generator ≈ 9,028,486 CCX
```

So the immediate integration target is not "BY replay plus naive reversible
branch generation". That would validate correctness but not savings.

`pattern_augmented_low_ratio_state_still_not_forward_complete` also kills the
simplest windowed escape. The 16 branch bits of the current window are
controlled by `(delta, h=g/f mod 2^16)`, but the next window's `h` is not. On
2k sampled secp256k1 trajectories:

```text
h16+pattern next-window mismatch = 66,905 / 70,000 = 95.58%
```

Thus a generator with only a 16-bit h register plus pattern history is not a
forward-complete state. A viable windowed generator needs a sliding
higher-precision 2-adic state, rank/high-bit payload, or a consumed-denominator
schedule. The principled tapered high-precision version was also measured in
`tapered_2adic_branch_generator_cost_is_still_too_high`: keep 560 2-adic bits
initially and drop one active bit per microstep. It halves the uniform-width
cost but is still far outside the SOTA margin:

```text
tapered generator compute             = 1,004,080 CCX
compute+uncompute per denominator      = 2,008,160 CCX
two denominators                       = 4,016,320 CCX
peak                                  = 3372q
```

So neither h16-only nor direct tapered high precision is the integration path.
The missing object is a genuinely windowed/self-cleaning division, not a more
careful version of per-bit branch generation. The self-cleaning mechanism itself
is now validated in miniature: `by_denominator_branch_history_self_cleans_on_reverse`
runs a 64-step, 96-bit 2-adic denominator pass, stores odd/A history, then
reverses the denominator and clears both history bits from the restored
pre-step state. It returns `f,g,delta`, odd history, A history, and phase clean:

```text
64-step 96-bit denominator forward+reverse: 87,808 CCX, peak 524q
```

This is constructive: branch history need not be generated by a separate
compute-copy-uncompute oracle. It can be produced while the denominator is
consumed and then erased while restoring the denominator. `lowword_pattern_oracle_is_cheap_and_clean`
then isolates the window branch oracle: copy low 16 bits of `(f,g)` into a
scratch 2-adic simulator, emit the 16 odd bits, CNOT them to persistent pattern
history, and reverse the simulator. It is exact and phase clean:

```text
16-step lowword pattern oracle: 5,952 CCX, peak 122q
```

So pattern generation for a window is not itself the expensive part. The
expensive part is the full-width selected denominator update. The straightforward
full-width controlled microstep window is measured in
`full_width_denominator_microstep_window_replay_is_not_enough`:

```text
16-step full-width denominator replay window: 26,256 CCX
35-window compute only:                         ≈918,960 CCX
compute+uncompute denominator:                ≈1,837,920 CCX
peak:                                           1128q
```

This is far above the fixed-matrix denominator target. The positive target is
now quantified by `tapered_fixed_matrix_denominator_budget_is_sota_shaped_if_selection_solved`:
assuming each 16-step matrix is already selected, and the active 2-adic width
shrinks from 560 by 16 bits per window, the existing fixed scaled-matrix
replacement costs

```text
35-window denominator compute       ≈ 303,828 CCX mean
compute+uncompute                   ≈ 607,657 CCX mean
max sampled compute                 = 329,595 CCX
```

The peak of the naive standalone cost circuit is high (`3424q`) because it
allocates old/new row buffers, but the Toffoli economics are SOTA-shaped.
`consumed_lowword_window_has_exact_quotient_update_and_pattern_inverse` pins down
the self-cleaning algebra for such a window. With `f=f_low+2^16 F` and
`g=g_low+2^16 G`, a selected matrix `P` gives

```text
(f',g') = P·(F,G) + (P·(f_low,g_low))/2^16
pattern + full quotient state reconstructs old state by sign(det(P))·adj(P)
max sampled lowword quotient correction |q| = 32,757
```

So the consumed low word does not need a separate 32-bit payload; the 16 branch
bits are exactly the determinant-history payload if the quotient state is kept
with enough precision. `lowword_pattern_and_q_oracle_is_still_cheap_and_clean`
turns this into a local reversible circuit by sign-extending the low words,
running 16 signed divsteps, copying out pattern and `(q0,q1)`, then reversing:

```text
pattern+q oracle: 9,408 CCX, peak 262q, qbits=34
```

This is cheap enough to sit next to a selected quotient update. A tempting
fixed-precision shortcut is dead, though:
`fixed_precision_2adic_denominator_branch_curve` shows that 64/96/128/.../256-bit
truncated denominator states predict roughly that many initial steps and then
mismatch essentially every 560-step trajectory:

```text
bits=64  mismatch_rate=1.0000 mean_first≈65
bits=128 mismatch_rate=1.0000 mean_first≈129
bits=256 mismatch_rate=1.0000 mean_first≈257
```

Thus field-width truncated branch generation is not a 1% approximate escape;
the required next object is still very precise: a **selected/fixed-matrix
window update** with algebraic sharing and production scheduling, not per-bit
full-width replay or fixed-width truncation.

The first savings-capable implementation must either:

1. window the denominator/control generator with enough high 2-adic state that
   it avoids 560 full-width microsteps, or
2. use a triangular point-add schedule where the denominator state is consumed
   and later cleaned from the output, avoiding a full compute-copy-uncompute
   branch generator.

This is the line between a benchmark smoke scaffold and a real SOTA candidate.

Default exact path note: `BULK_PREFIX_SAFE_ITERS` was raised from 313 to 375
after env sweeps found 375 phase-clean and 377/380/390 phase-unsafe. This gives
a real but small exact improvement:

```text
old default: 4,132,750 Toffoli @2716q
new default: 4,111,918 Toffoli @2716q
```

It is folded into the baseline but is not a SOTA route.

Approximate Kaliski cutoff is not a shortcut under the current phase-clean
contract. Env-gated threshold probes (`KAL_PAIR1_ITERS`, `KAL_PAIR2_ITERS`) show
that even tiny reductions fail the alternate-seed phase guard:

```text
pair1=403,pair2=403: classical=0, phase_batches=1
pair1=407,pair2=400: classical=0, phase_batches=2
pair1=400,pair2=403: classical=1, phase_batches=1
pair1=390,pair2=390: classical=16, phase_batches=49
```

So current Kaliski iteration counts are effectively tight unless the phase
contract itself is changed. The SOTA path must be a different division
architecture, not truncating the existing one.

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

Immediate next target update after BY invalidation: **Program B only if it
produces a new primitive-level win before point-add wiring**.  The old Program A
BY denominator route is no longer a production target with current primitives.

### Relaxed-2800q quotient-stream checkpoint

The centered Euclid quotient stream is still not a Google-low-qubit candidate:
its parser/boundary state is above the `~600--663q` scratch target.  It remains
useful as a current-cap (`<=2800q`) fallback only if quotient extraction is
packed and phase-clean.

The previous non-restoring ledger was demoted after charging two hidden
corrections: a final negative-remainder fix for the floor extractor and a
nearest-quotient correction.  `direct_centered_nonrestoring_rounding_numerator_reopens_relaxed_margin`
tests the direct-centered variant

```text
q = floor((|u| + |v|/2) / |v|)
```

so the nearest correction is folded into the numerator before extraction.  On
4096 secp samples:

```text
p99 signed-digit payload       = 397 bits
p99 centered binary payload    = 336 bits
p99 quotient count             = 118
p99 final-negative corrections = 69
one-way extraction ledger      = 403,712 CCX
projected relaxed point-add    = 2,804,648 CCX
gap to 3M                      = -195,352 CCX
```

This reopens only the relaxed-2800q route, not the Google-low-qubit route.  The
first hard-piece charge then demotes the current primitive set:
`direct_centered_nonrestoring_current_signed_digit_primitive_kills_margin`
implements the obvious sign-selected add/sub digit as two generic controlled
q-q adders.  It measures

```text
ccx64                         = 516
scaled 257-bit digit           = 2,041 CCX
projected relaxed point-add gap = +2,639,228 CCX
```

That generic implementation is not a viable path.  The follow-up
`direct_centered_nonrestoring_fused_signed_digit_revives_margin` folds the
choice into one two's-complement add:

```text
sign=0: rem += ~shifted_v + 1 = rem - shifted_v
sign=1: rem +=  shifted_v
```

by toggling `shifted_v` under `!sign` and using `!sign` as the fast-Cuccaro
carry-in.  The small-basis simulator check confirms data, addend cleanup, and
phase cleanliness.  Measured cost:

```text
33-bit digit                  = 32 CCX
65-bit digit                  = 64 CCX
257-bit digit                 = 256 CCX
257-bit digit peak            = 772 qubits
projected relaxed point-add gap = -195,352 CCX
```

`direct_centered_nonrestoring_toy_packed_extractor_is_phase_clean` then builds
the first extractor-level toy: fixed boundary, aligned divisor scratch,
shift-down reuse, packed signed-digit history, final negative cleanup, divisor
uncopy, and basis-simulator phase check.  The toy measures 81 CCX for six
11-bit digits:

```text
six 11-bit toy extractor       = 81 CCX
toy peak                       = 55 qubits
full ledger with 2n-1 final fix = 2,925,008 CCX
gap to 3M                      = -74,992 CCX
```

This promotes the direct-centered route from "digit primitive only" to
"toy-extractor viable" under the current 2800q cap.  It is still not production
ready: the missing hard piece is a real 256-bit packed extractor with dynamic
boundary/parser cleanup and integration with coefficient replay.

The margin is extremely narrow.  `direct_centered_nonrestoring_inactive_digit_tax_kills_margin`
charges just one CCX for every skipped static digit position:

```text
p99 static digit positions     = 30,208
p99 active digit payload       = 397
p99 inactive digit positions   = 29,811
one-CCX inactive tax gap       = +44,252 CCX
```

So a production extractor must avoid Toffoli work on inactive quotient-digit
positions.  A controlled full-width scan over all 256 positions per quotient is
not acceptable even with a one-CCX inactive body.

Follow-up: the full 8-level alignment barrel in that ledger was pessimistic for
the sampled direct-centered traces.  `direct_centered_nonrestoring_bounded_barrel_revives_inactive_tax_margin`
uses 32,768 samples and observes maximum non-restoring digit width 24, so the
initial quotient alignment can be modeled with a 5-bit barrel rather than an
8-bit full-width barrel:

```text
bounded barrel bits            = 5
full-barrel relaxed gap        = -74,992 CCX
bounded-barrel relaxed gap     = -437,488 CCX
bounded-barrel saving          = 362,496 CCX
one-CCX inactive tax gap       = -318,244 CCX
```

This reopens the direct-centered relaxed-2800q route as an extractor-engineering
problem.  Fixed K-bit quotient slots are still invalid (K=5 fails 999,633 ppm,
K=12 still fails 49,804 ppm), so the next hard piece is a real packed active
digit loop with a bounded binary initial alignment, not a static 256-position
scan and not fixed-width quotient slots.

Caveat: the 5-bit barrel is a sampled/average-trace model, not a standalone
exact full-domain circuit.  `direct_centered_bounded_barrel_needs_full_domain_fallback`
checks adversarial small denominators and finds `x=1` requires 256
non-restoring digits, hence the full 8 alignment bits.  Charging those three
missing layers back into the inactive-tax ledger adds 362,496 point-add CCX and
flips the one-CCX inactive gap from `-318,244` to `+44,252`.  The direct-centered
route therefore needs either a phase-clean wide-quotient fallback or a way to
make the high alignment layers classically/rarely executed; a bare 5-bit exact
barrel is not correct for the whole field.

Follow-up: the barrel mechanics are not the blocker if the parser can expose
alignment metadata as phase-clean classical bits.  `direct_centered_classical_alignment_metadata_would_remove_barrel_blocker`
uses bit-controlled swaps and measures `0` barrel Toffoli in the toy shifter.
In the exact 8-bit inactive-tax ledger this saves 966,656 point-add CCX, moving
the gap from `+44,252` to `-922,404`.  The remaining hard question is whether
those alignment bits can be produced classically without a dense measurement
phase correction or a static quantum scan.

## 6. Post-BY ground-up attempt: Strategy E slope-coordinate map

New non-BY candidate: make the slope the live y-coordinate and avoid an
independent `lam` register entirely.

```text
dx = Px-Qx
dy = Py-Qy
m  = dy/dx
Rx = m^2 - dx - 2Qx
Ry = -m*(Rx-Qx) - Qy
```

This gives a triangular affine permutation:

```text
(tx,ty)=(dx,dy)
DIV(tx,ty)        -> (dx,m)
tx := m^2-tx-2Qx -> (Rx,m)
IMUL(tx,ty)       -> (Rx,Ry)
```

`strategy_e_slope_coordinate_formula_passes_200` validates the formula on 200
random secp256k1 point pairs.  It is the cleanest one-inversion-looking affine
schedule so far because `dx -> Rx` is an involution once `m` is live, and no
fresh output registers are required.

Fast invalidation: with current reversible primitives, `IMUL(tx,ty)` is not a
schoolbook multiply.  The map `(c,m)->(c,m*c)` is not globally reversible when
`c=0`, and known exact/product-clean constructions need an inverse-sized cleanup
object.  That collapses Strategy E back to the measured pair2
scaled-product-clean world.  The executable budget guard
`strategy_e_slope_coordinate_budget_requires_new_inplace_variable_multiply`
records:

```text
current known product-clean route ≈ 2,988,510 Toffoli before safety margin
if IMUL were schoolbook-like       ≈ 2,022,750 Toffoli
needed IMUL saving                 ≈   965,760 Toffoli
```

Decision: **Strategy E is algebraically validated but circuit-invalidated with
current primitives**.  Do not wire it using existing Kaliski/BY product-clean
machinery; it only becomes SOTA-shaped if we invent a genuinely new phase-clean
in-place variable multiply/divide primitive at roughly schoolbook cost and
without raw inverse/history storage.

First primitive attempt after this decision: destructive Montgomery IMUL.  The
idea was to scan the old multiplier bits into a Montgomery add-and-halve
accumulator, zeroing the multiplier as we go and leaving a product accumulator.
Forward algebra works (`t=x*y*2^-n mod p`), but local reversibility fails.
`destructive_montgomery_product_is_algebraically_promising_but_not_locally_reversible`
finds a concrete small-prime reachable window with 512 valid predecessors for
the same post-accumulator.  Therefore consumed multiplier bits cannot be cleared
from the accumulator alone; the circuit would need history/checkpoints or a
nonlocal inverse, i.e. the same product-clean obstruction.  The possible rescue
"reverse-decode from the final product with a small ambiguity frontier" is also
killed: `destructive_montgomery_reverse_trellis_needs_field_sized_state` expands
to essentially the full `[0,2p)` frontier (`n=8,10,12 -> 502,2042,8186`).
Kill this primitive.

Second primitive attempt: MBUC product cleanup with a phase-only quotient oracle.
If we compute `z=x*y`, X-measure the old multiplier, and only correct phase, the
correction is a known-mask bit of `z/x mod p`.  The toy-field ANF test
`mbuc_product_cleanup_phase_oracle_is_not_low_degree_on_toy_field` shows degree
15/16 and 32518/65536 monomial density for `p=251`.  This kills the cheap
sparse/low-degree phase-correction hope; MBUC cleanup still hides a quotient.

Adjacent kickmix check: maybe we do not measure the old multiplier, only a raw
schoolbook product scratch `t=x*y` while preserving both inputs.  This would
replace the gate-level inverse of the temporary multiplier with phases
`(-1)^(mask·(x*y))`.  `raw_product_measurement_phase_is_dense_not_free_kickmix`
shows that typical product-bit masks are already carry-dense: at toy width
`n=10`, the all-bit mask has degree `19/20` and density `427812/1048576`, and a
single high-bit mask has degree `19/20` and density `120581/1048576`.  So raw
product-scratch MBUC is not a hidden free IMUL/uncompute; high product bits are
full carry logic, not just quadratic `x_i y_j` phases.  A redundant carry-save
version also fails the same smell test: `carry_save_product_scratch_mbu_still_has_dense_phases`
compresses partial products to two carry-save rows and still sees degree `16/16`,
density `20440/65536` for an all-row mask at `n=8` (high-column mask degree
`15/16`, density `3602/65536`).

Third primitive attempt: MBUC the Montgomery loop's internal quotient/history
bits instead of the old multiplier.  In the `(x,old_y)` frame this history is
structured (`n=8` density 2440/65536; `n=10` density 31684/1048576), but an
actual in-place multiply must phase-correct in the surviving `(x,z)` frame.
`montgomery_q_history_phase_in_output_frame_is_dense_dead` gets degree 16/16 and
31032/65536 density for `p=251`, so the structured history becomes quotient-like
again after old-y is gone.  Kill this primitive too.

Fourth primitive attempt: preserve enough denominator information to recover
`λ` from `λ² = Rx + dx + 2Qx` instead of a second division.  This just asks for
square root.  `lambda_square_cleanup_would_require_dense_sqrt_phase` shows the
canonical square-root phase is dense/high-degree on toy primes (`n=8,10,12`
densities `126/256`, `502/1024`, `2072/4096`, degree up to full).  For
secp256k1 this is an exponentiation, not a cleanup primitive.  Kill it.

Forward-only coefficient-transform Kaliski also got a stronger exactness check.
The earlier secp collision sample suggested full post-state `(u,v,r,s,f)` might
recover branches, but `secp_local_poststate_predecessor_branch_is_ambiguous`
now directly enumerates locally consistent inverse branches on actually reached
secp tagged poststates and finds ambiguity in `7205/8140 = 88.5%` of sampled
steps (`hist=[0,935,1833,1779,3593]` for 0..4 predecessor branch counts).
`exhaustive_toy_full_poststate_does_not_recover_forward_branch` also finds exact
tagged-DIV collisions even when the reverse iteration index and full post-state
are known (`n=4,5,6` conflicts `108,1200,5760`).  The approximate exception
version is not viable either: `tagged_full_poststate_branch_ambiguity_is_not_a_rare_exception`
finds ambiguous full-poststate occurrences at about 19%, 22%, 23%, 23%, and 24%
for toy `n=4..8` with nonzero tag `s0=x+y`, so the ambiguity is not a tiny
patchable tail.  `changing_linear_x_tag_does_not_fix_poststate_branch_ambiguity`
checks `s0=y+beta*x` and finds identical ambiguity fractions for nonzero beta;
`any_x_only_additive_tag_is_only_a_bijection_not_branch_history` extends this to
nonlinear additive tags `s0=y+h(x)` and gets the same `n=8` fraction `0.2432`.
For fixed `x`, such tags are only a bijective relabeling of the coefficient
scalar, not branch-history encoding.  Direct regeneration from the preserved initial `x` is not cheap either:
`initial_x_to_branch_history_oracle_is_dense_on_toy_kaliski` sees branch-history
parity as full-degree/near-half-density in `x` (`n=8: degree 8, density 116/256;
`n=12: degree 11, density 1976/4096`).  `exact_branch_history_has_field_entropy_lower_bound`
adds the information floor: complete branch sequences are injective in toy
fields and collision-free on 4096 secp samples, so exact history encoding needs
~256 bits.  `coefficient_transform_history_floor_misses_low_qubit_budget` then
combines this with the remaining layouts: r-as-output coefficient DIV has a
768q scratch floor and second-channel DIV has a 1280q floor, before flags or
adder workspace, versus Google's low-qubit allowance of 663q beyond `tx,ty`.
So branch cleanup cannot simply inspect live post-state, use a tiny x-oracle,
compress to ~64 bits, fit exact coefficient-transform history into low-qubit,
or ignore a negligible poststate-collision tail; it needs a different invariant.
Kim-style unconditional Kaliski was checked through the same lens:
`kim_unconditional_poststate_does_not_recover_branch_flags` keeps the wide
postponed-reduction arithmetic but enumerates exact local inverse branches for
`swap`/`both_odd`, finding `5178/10240 = 50.6%` ambiguous reached secp
poststates.  Thus Kim removes the terminal flag but does not by itself remove
branch history.  Luo/PZ register sharing is also not a hidden low-gate SOTA
route: `luo_pz_gate_slope_is_not_point_add_sota_shaped` converts the published
~976 n^3 whole-ECDLP slope into `31,981,568` Toffoli per point-add proxy, about
11.85x/15.23x above Google low-qubit/low-gate point-add budgets.

Coefficient-transform refinement checked: a single Kaliski coefficient pair
cannot both preserve `x` and expose `y/x` by just using a constant tag.  If
`r0=ρ` and `s0=y+β`, the outputs are `s=ρx` and
`r=k*y+(ρa+βk)`.  `single_coefficient_pair_cannot_preserve_x_and_expose_quotient_by_constant_tag`
shows sampled `(a,k,1)` are affine-non-collinear, so no constants `ρ,β,C` make
`ρa+βk=C`.  The data-dependent cancellation escape is also ugly:
`a_coefficient_cancellation_is_dense_on_toy_kaliski` gets full-degree,
near-half-density ANFs for mask bits of `a(x)` (`n=8: 132/256`,
`n=12: 2094/4096`).  A one-pair DIV still needs either a second channel or a
different triangular Euclidean transform.

secp256k1-specific j=0 automorphism check: `(x,y)->(βx,y)` gives the identity
`(x-Qx)(βx-Qx)(β²x-Qx)=(y-Qy)(y+Qy)`, hence
`λ=(x²+Qx*x+Qx²)/(y+Qy)`.  `secp_j0_endomorphism_slope_denominator_swap_identity_passes`
validates the denominator swap on 200/200 secp samples with zero y-sum
exceptions.  But it is not a cleanup primitive: from the output frame,
`λ=(Rx²+Qx*Rx+Qx²)/(Qy-Ry)`, and
`endomorphism_slope_swap_cleanup_phase_is_still_dense` gets toy n=10 phase ANF
`degree=20/20`, `density=297278/1048576`.  Restricting interpolation to valid
j=0 curve outputs still gives growing min degrees `1,3,3,3,4` through n=12 in
`endomorphism_slope_support_degree_still_grows`.  The automorphism is useful
algebra, but still quotient-like for reversible cleanup.

Coordinate-model escape check: `efficient_curve_model_transforms_need_missing_torsion`
records that secp256k1's prime-order group has `order mod 2 = 1` and
`order mod 3 = 1`.  Base-field Montgomery/Edwards models need rational
2-torsion and Hessian/twisted-Hessian models need rational 3-torsion; birational
maps preserve torsion.  So cheap complete-law model changes are unavailable for
this exact affine benchmark.  Projective/isogenous variants still need affine
conversion cleanup.

Fifth primitive/global attempt: compute the affine output out-of-place, X-measure
the old point, and correct phase using `P = R - Q`.  The toy test
`top_level_mbuc_of_old_point_requires_dense_point_subtraction_phase` shows the
phase oracle for a mask of old point bits as a function of `(R_x,R_y)` has degree
15/16 and density 19540/65536 on `y^2=x^3+7 mod 251`.  The obvious caveat is
that the output is curve-supported, not full-domain.  The follow-up
`curve_support_mbuc_phase_still_scales_not_constant_degree` solves the
support-restricted interpolation problem and still sees growing minimum degree:
`n=4,6,8,10,12 -> 1,3,3,4,4`.  Directly measuring the affine slope is no better:
`measuring_lambda_after_affine_add_still_needs_growing_degree_phase` corrects
`λ=-(R_y+Q_y)/(R_x-Q_x)` from curve-supported outputs and gets minimum degrees
`2,3,3,4,4`.  A slope-pencil coordinate rewrite also hits the same wall:
`pencil_slope_coordinate_needs_dense_root_choice_phase` stores the slope of a
line through fixed `Q`, but the missing root-choice bit is a dense sqrt-like
phase (`n=8: degree 7, density 114/256; n=12: degree 11, density 2006/4096`).
This matches the dimension-threshold/division/sqrt story (`sum_i<=d C(2n,i) >=
~2^n`, giving `d≈0.22n≈56` at secp256k1), not a constant-degree identity.
Sequential one-coordinate MBUC was checked too:
`sequential_old_coordinate_mbuc_still_has_growing_phase_degree` keeps old `x`
live while measuring old `y`, and the support-restricted degree still grows
`n=4,6,8,10 -> 1,2,2,3` (dimension extrapolation for `3n` variables is still
~49 at secp256k1).  Thus generic top-level/sequential MBUC is just
dense/high-degree point-subtraction phase correction, not a cheap
single-inversion cleanup.

## 7. Post-E attempt: Kaliski scale absorption by denominator pre-scaling

This is a smaller structural lever than deleting an inversion, but it attacks a
measured ~206k correction-loop component.  Since raw Kaliski returns
`-v^-1 * 2^iters`, feed it `v' = 2^iters v`; the exposed inverse becomes exact
and the slope halving / cleanup doubling loops disappear.

Pair1 and pair2 were wired behind `KAL_PRESCALE_PAIR1_SAFE=1` /
`KAL_PRESCALE_PAIR2_SAFE=1` using an exact Cuccaro constant prescaler and passed
the full harness:

```text
pair1: avg_toffoli=4,786,373  qubits=2,972  clean
pair2: avg_toffoli=4,771,009  qubits=2,969  clean
```

The all-exact prescaler makes the path much worse than default, but the
algebraic lever is real.  A fast version using the existing measurement-based
constant multiplier was classically correct but phase-unsafe
(`altseed_phase_batches=1`).  The mixed diagnostic variant
`KAL_PRESCALE_PAIR1_MIXED=1` / `KAL_PRESCALE_PAIR2_MIXED=1` use exact q-q
add/sub but fast modular shifts and are clean:

```text
pair1 mixed:  avg_toffoli=4,223,465  qubits=2,972
pair2 mixed:  avg_toffoli=4,220,405  qubits=2,969
both mixed:   avg_toffoli=4,331,952  qubits=2,972
pair1 folded: avg_toffoli=4,223,465  qubits=2,969
pair2 folded: avg_toffoli=4,220,405  qubits=2,965
```

The folded variants push the mixed prescale into Kaliski `v_w` init, removing
the separate scaled denominator register.  They confirm the right integration
point but only save 3-4 qubits in the current peak profile.

So the phase culprit is the fast q-q add/sub in constant multiplication, not the
fast double/halve scale walk.  The mixed prescaler is still ~108-112k above the
true default path, but it is close enough to sharpen the target.

Decision: keep this as a narrow future primitive target only.  A useful
production version requires a secp256k1-specific phase-clean shifted-add
prescaler for `2^iters mod p` (e.g. `2^407 = 2^151(2^32+977) mod p`) costing
well under the current mixed compute+uncompute overhead and ideally folded into
the Kaliski `v_w` initialization so it does not add a persistent n-bit
denominator copy.  Without that primitive, scale absorption is not SOTA-relevant.

## 8. Direct controlled-constant Solinas correction (qubit tradeoff)

The default peak trace after the prescale work showed the global maximum at
Kaliski backward `bk_step6_7_8`, specifically inside modular halve's controlled
Solinas correction.  The old `csub_nbit_const_fast` materializes the classical
addend `ctrl ? (2^32+977) : 0` as a full 256-qubit register before running the
fast ripple.  A direct controlled-constant borrow ripple removes that loaded
addend register and keeps only the carry/borrow chain.

Implemented env-gated helpers:

```text
csub_nbit_const_direct_fast
cadd_nbit_const_direct_fast   # present, but not production-clean yet
KAL_DIRECT_CONST_HALVE=1
KAL_DIRECT_CONST_DOUBLE=1
```

Measured result:

```text
KAL_DIRECT_CONST_HALVE=1
avg_toffoli = 4,121,014
qubits      = 2,715
emitted_ops = 29,985,672
clean       = yes
```

This is not a default Toffoli win (`+9,096` Toffoli versus 4,111,918), but it is
an important qubit-tradeoff result: one peak qubit is removed at very small
Toffoli cost, much cheaper than the earlier dirty-venting halve path.  A peak
trace with the flag enabled moves the bottleneck to Kaliski backward step4:

```text
DEBUG peak_qubits=2715 at phase='bk_step4'
near peaks: bk_bulk_step4, kal_bulk_step4, kal_bulk_step6_7_8
```

The analogous direct controlled-constant add path is logically valid but
phase-cliff sensitive.  Standalone small basis/phase tests for both direct cadd
and csub pass.  In the full point-add harness, direct cadd fails with the tuned
375 Kaliski bulk prefix but passes with the conservative 370 prefix:

```text
KAL_DIRECT_CONST_DOUBLE=1
altseed_classical_total=1
altseed_phase_batches_total=2

KAL_DIRECT_CONST_DOUBLE=1 KAL_BULK3_ITERS=370
avg_toffoli = 4,121,506
qubits      = 2,716
clean       = yes
```

So keep only the halve flag as the useful validated low-qubit tradeoff.  Direct
cadd is an env-gated tested primitive, but it gives no qubit reduction and needs
a less aggressive bulk prefix, so it is not a default-path optimization.  Both
direct cadd+csub together under `KAL_BULK3_ITERS=370` are also clean:

```text
KAL_DIRECT_CONST_HALVE=1 KAL_DIRECT_CONST_DOUBLE=1 KAL_BULK3_ITERS=370
avg_toffoli = 4,130,602
qubits      = 2,715
emitted_ops = 29,250,534
clean       = yes
```

This is useful only for low-emitted-op / low-qubit tradeoff work (`-1` qubit,
`-1.5M` emitted ops, `+18,684` Toffoli versus default).

Default exact path after adding the env-gated primitive initially remained
unchanged at `4,111,918 @2716q`.  A separate local threshold probe then found a
small exact win in the Kaliski coefficient range bound.  The old
`R_SMALL_THRESHOLD=255` controlled when `r := 2r` / reverse halve could skip the
Solinas correction because `r` was guaranteed small.  Env probes showed:

```text
KAL_R_SMALL_THRESHOLD=256  avg_toffoli=4,110,898 clean
KAL_R_SMALL_THRESHOLD=257  avg_toffoli=4,109,878 clean
KAL_R_SMALL_THRESHOLD=258  phase failure (global_phase=0x0000000008000000)
KAL_R_SMALL_THRESHOLD=258 KAL_BULK3_ITERS=370 altseed phase failure=1
```

Default is now `R_SMALL_THRESHOLD=257`:

```text
avg_toffoli = 4,109,878
qubits      = 2,716
clean       = yes
```

This is only a `-2,040` Toffoli local gain, but it is exact and default-clean.
