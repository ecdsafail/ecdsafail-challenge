# Ground-up architecture under ~600 scratch qubits

User framing: Google low-qubit means roughly **600-660 non-data qubits** beyond
`tx,ty` (512 quantum input/output qubits). This document ignores local tuning
and asks what can possibly fit.

## 1. Budget arithmetic

At `n=256`:

```text
data registers: tx, ty = 512q
Google low-qubit total: 1175q
scratch beyond tx,ty: 663q
user mental model: ~600q
```

So a viable low-qubit point-add can have at most:

- two full n-bit scratch registers (`2n = 512`) plus ~90-150 small bits, or
- one full n-bit scratch register plus a compact inversion state, or
- heavy reuse of `tx,ty` as algorithmic work registers.

It **cannot** have three extra n-bit registers, and it definitely cannot have
current Kaliski's `u,v,r,s,m_hist` state.

Current peak has, beyond `tx,ty`:

| live object | qubits |
|---|---:|
| slope `lam` | 256 |
| Kaliski `u,v,r,s` | 1024 |
| `m_hist` | 403-407 |
| transients | 250-520 |
| **non-data total** | **~2200** |

We need to remove roughly **1500 non-data qubits**, not 50.

## 2. Consequence: treat data registers as part of the algorithm

A 600-scratch design cannot say "keep `tx,ty` pristine, compute everything in
ancilla, swap outputs, then uncompute". That Bennett pattern leaves old
`Px,Py` or `dx,dy` in fresh registers and immediately exceeds the budget.

The only plausible low-qubit pattern is:

1. Mutate `tx,ty` into useful intermediates (`dx,dy`, coefficient registers,
   accumulators).
2. Use at most two additional n-bit work registers.
3. Arrange the final inverse/cleanup so that running a reverse transform writes
   the desired output into `tx,ty` rather than restoring the input.

This is the right abstraction: **we need a reversible data transform, not a
Bennett-clean subroutine call.**

## 3. Inversion-state lower bound

Any Euclidean inverse needs, in some representation:

- a denominator state (`u/v` or equivalent), and
- coefficient information connecting the denominator to the inverse.

Current Kaliski stores this as `(u,v,r,s)` plus history. In 600 scratch, the
only way to keep Kaliski-like inversion alive is to fold at least two of those
four n-bit roles into `tx,ty`.

A minimal Kaliski-like layout would have to look like:

| role | storage |
|---|---|
| denominator input `v=dx` | `tx` or one scratch copy |
| other gcd register `u` | scratch A |
| coefficient/output register | scratch B or `ty` |
| second coefficient register | `ty` or eliminated |
| history | not stored, or <=~100 bits |

This is exactly why `m_hist` elimination alone is insufficient: even without
history, the four n-bit Kaliski roles already exceed 600 scratch unless they
are folded into the data registers.

## 4. New structural idea: use Kaliski's coefficient transform on `ty`

Instead of treating Kaliski as an ancilla subroutine, seed its coefficient
register with the data value `dy`.

Use a canonical-mod-p coefficient version of Kaliski. For a fixed denominator
`dx`, the coefficient-side update is a linear transform:

```text
(r_final, s_final)^T = T(dx) (r_initial, s_initial)^T
```

The test module `kaliski_linear_transform.rs` verifies empirically for the
current 407-iteration branch sequence that:

```text
T(dx) = [[ a(dx), k(dx) ],
         [ dx,    0     ]]

k(dx) * dx = -2^407  (mod p)
```

Therefore:

```text
T(dx) * (0, 1)  = (k, 0)          raw inverse
T(dx) * (0, dy) = (k*dy, 0)       scaled slope, ty consumed to zero
T(dx) * (1, 0)  = (a(dx), dx)     exposes dx in the second coefficient
```

This is the first genuinely low-qubit-looking Kaliski algebra found in this
repo: `ty` can be consumed into the coefficient transform instead of being kept
as an external data register plus a separate multiplication `dy * inv(dx)`.

### Why this matters

If we could finish the point-add while the coefficient transform is live, then
run Kaliski backward, `ty` could be written to an arbitrary target value.
Specifically, to finish backward with:

```text
r_initial = 0
s_initial = Ry
```

the state *before* backward must be:

```text
T(dx) * (0, Ry) = (k*Ry, 0)
```

But the dy-seeded forward naturally gives:

```text
(k*dy, 0)
```

So the exact structural subproblem is:

```text
add  k * (Ry - dy)  into r, with s=0, without a second inversion.
```

This is crisp. It replaces the vague "one-inversion cleanup obstruction" with
one algebraic target.

## 5. Current obstruction in the coefficient-transform frame

We know:

```text
k = -2^407 / dx
L = k*dy = scaled(lambda)
Ry = -lambda*(Rx-Qx) - Qy
```

Then:

```text
k*(Ry-dy)
  = -k*lambda*(Rx-Qx) - k*Qy - k*dy
```

The live dy-seeded state gives `L = k*dy`, but not `k` itself. The `k*Qy`
term is the sticking point: multiplying a classical `Qy` by raw `k` requires
access to `k`, i.e. the raw inverse, not just the scaled slope.

This explains why the usual one-inversion schedules leak a slope copy: they
have enough information for `lambda`, but not enough to rewrite the Kaliski
coefficient pair to make backward output `Ry`.

## 6. What would make this a breakthrough?

The coefficient-transform idea becomes a 600-scratch / SOTA route if we can do
one of the following:

1. **Expose both `k` and `k*dy` using the two coefficient registers.**
   Since `T(dx)*(1,0)=(a,dx)` and `T(dx)*(0,1)=(k,0)`, maybe a different
   initialization of `(r,s)` plus the already-live `tx=dx` can recover `k`
   or `k*Qy` without another full inverse.

2. **Choose a different y-coordinate convention so the `k*Qy` term vanishes.**
   Work with shifted `Y` coordinates, e.g. store `Y+Qy` or `Y-Qy`, so that the
   final backward target is `Ry+Qy` instead of `Ry`. If the benchmark output
   can be recovered by a final classical add/sub, this may remove the raw-`k`
   constant term.

3. **Use the `r_initial` channel deliberately.**
   We do not necessarily need backward to end with `r_initial=0`; it could end
   with a known classical constant and then be X-freed. This changes the target
   from `T*(0,Ry)` to `T*(C,Ry) = C*(a,dx)+Ry*(k,0)`, giving an additional
   live `dx*C` in `s_final` and maybe a way to absorb the constant term.

4. **Run a tiny second coefficient transform, not a second full inversion.**
   If only the `k*Qy` term is missing and `Qy` is classical, maybe a
   classical-seeded coefficient pass can be folded into the same branch history
   or a short replay. This would be far cheaper than a full second Kaliski if
   it reuses the branch sequence.

These are structural, not micro. Any one of them could delete the second
inversion and land near 2.5M Toffoli. If all fail, two-inversion SOTA must come
from jumped/windowed Kaliski instead.

## 7. The real primitive we need: in-place modular division

The low-qubit point-add can be phrased around one primitive:

```text
DIV:  (x, y) -> (x, y/x mod p)
```

with all scratch cleaned and `x` preserved. If `DIV` costs roughly one current
Kaliski invocation and fits in ~600 scratch, then point-add becomes:

```text
tx = Px-Qx = dx
ty = Py-Qy = dy
DIV(tx, ty)                    // ty = λ
// tx = λ² - dx - 2Qx = Rx
// ty = λ(Qx-Rx) - Qy = Ry, as an in-place multiply-by-(Qx-Rx)
```

This is conceptually **one inversion**, but it avoids the slope-copy cleanup
obstruction by never materializing `x^-1` as an independent output. It is the
clean abstraction that matches the 600-scratch target.

Current code does **not** have this primitive. `with_kal_inv_raw` computes a
raw inverse into an ancilla and then has to Bennett-clean the inverse state.
The coefficient-transform probe above is a first attempt to derive `DIV` from
Kaliski by seeding the coefficient register with `y`.

### Why a quotient-copy DIV does not fit 600 scratch

A tempting DIV implementation is:

1. Run Kaliski forward with `tx` as the denominator state and `ty` as the
   coefficient seed; this can fit with scratch `u,r` if history is eliminated.
2. Extract/copy the quotient to a separate n-bit register.
3. Run Kaliski backward to restore/clean the Euclidean state.
4. Clear old `ty` and swap in the quotient.

But during backward this needs simultaneously:

```text
tx as v-state, ty as s-state, scratch u, scratch r, quotient_copy
```

That is **three n-bit scratch registers** (`u,r,quotient_copy = 768q`) beyond
`tx,ty`, before flags/history/transients. It already violates the ~600-scratch
budget. Therefore a low-qubit DIV cannot copy the quotient across backward.
The backward transform itself must write the desired output into `ty`.

This is why the coefficient-transform target `(k*Ry,0)` matters: it is not an
optional elegance issue; it is the only way to avoid the third n-bit scratch
register.

## 8. Shifted-Y algebra: first fast invalidation

Try to save the coefficient-transform path by changing the y-coordinate
convention. Let the seed be `S0 = Py + a·Qy = dy + (a+1)Qy`, and the desired
backward output be `S1 = Ry + b·Qy`. The required Kaliski-coefficient update is

```text
k*(S1-S0)
```

where `k = raw_scale/dx` and `L = k*dy = raw_scale*λ` is available.

Compute:

```text
Ry - dy = λ(3Qx - λ²) - Qy
S1 - S0 = λ(3Qx - λ²) + (b-a-2)Qy
```

Choosing `b=a+2` removes the raw `k*Qy` term, but leaves

```text
k * λ * (3Qx - λ²)
  = L * (3Qx - λ²) / dx
  = L * (Qx - Rx - dx) / dx
```

which still requires division by `dx`, i.e. raw `k` or a second inverse. Thus
**affine shifts of Y do not solve the coefficient-transform obstruction**.
They move the missing term from `k*Qy` to `k*λ*(...)`.

## 9. Two-channel coefficient search: partial reduction

Let initial and desired post-backward coefficient pairs be:

```text
initial:       (r0, s0)
desired final: (rF, sF)
```

Kaliski's coefficient transform gives:

```text
current before body = T(dx)(r0, s0)
target before back  = T(dx)(rF, sF)
```

Difference to implement inside the body:

```text
Δr = a(dx)*(rF-r0) + k(dx)*(sF-s0)
Δs = dx*(rF-r0)
```

The unknown coefficient `a(dx)` is just as hard to expose as `k(dx)`. Therefore:

- If the scratch `r` register must be freeable at the end, `rF` must be a known
  constant (usually 0).
- To avoid an `a(dx)` term, we need `rF = r0`.
- Combining these forces `r0 = rF = constant`.

Under those conditions the two-channel problem collapses exactly to the
shifted-Y search in §8, which is already invalidated.

If `r0` is data-dependent, then either:

1. `rF=r0` and the scratch `r` register exits with data-dependent garbage, or
2. `rF` is constant and the body must cancel `a(dx)*(rF-r0)`, requiring access
   to the other unknown transform coefficient `a(dx)` in addition to `k(dx)`.

So the simple two-channel affine family does **not** rescue coefficient-
transform Kaliski. A viable design must be more radical: make `r` itself one
of the final output registers, or use a different Euclidean transform whose
matrix has a triangular form better suited to DIV.

## 10. Output-register use of `r`: reduces to self-cleaning forward Kaliski

If Kaliski coefficient `r` is allowed to become final `ty`, we can avoid the
quotient-copy lower bound:

```text
start:   tx=x, ty=y, scratch r=0, scratch u=p
forward coefficient Kaliski using ty as s:
         r = k*y  (scaled quotient), ty/s = 0, u=1, tx/v=0-ish
scale r by a known constant -> y/x
compute remaining point-add arithmetic using r as the slope/output channel
swap r into ty at the end, leaving r=0
```

This would fit the scratch budget (`u` + `r` = 512q plus small flags) **if and
only if** forward Kaliski can be made self-cleaning, i.e. no persistent
`m_hist` and no backward pass.

The new test `end_state_needs_coefficient_registers_to_recover_branch` shows:

- `(u,v,f)` at iteration end does **not** determine the branch; denominator-only
  recovery has collisions.
- `(u,v,r,s,f)` at iteration end **does** determine the branch on 200×407
  sampled canonical coefficient trajectories.

So a self-cleaning forward Kaliski is not information-theoretically dead for
nonzero coefficient seeds, but its branch-recovery predicate must inspect the
coefficient registers. It is not the cheap 4-bit start-state formula.

The follow-up test `zero_coefficient_seed_loses_branch_information` shows a
critical exactness problem: if the coefficient seed is zero, then even full
`(u,v,r,s,f)` has collisions because `r=s=0` carries no trajectory signal.
Approximate tolerance makes this rare exceptional set acceptable in principle
(`dy=0` is negligible for random points), but it is not the main obstacle.

The stronger follow-up test
`low_bit_end_state_branch_classifier_is_not_approx_good_enough` trains the
best majority lookup from the low 3 bits of `(u,v,r,s,f)` and tests it on
disjoint samples. Error is >50%. So the needed end-state branch predicate is
not a small low-bit heuristic. The branch information is globally present in
full coefficient state, but extracting it cheaply appears hard.

### Approximate-tolerant tag breakthrough: seed with `y+x`

User clarified that **~1% total failure is tolerable**. This makes a tag viable,
but the tag must not introduce an unremovable raw-`k` term. The right tag is the
denominator itself:

```text
s0 = y + x
T(x)*(0, y+x) = (k*y + k*x, 0) = (k*y - 2^ITERS, 0)
```

because `k*x = -2^ITERS` is a known constant. Therefore:

```text
k*y = r + 2^ITERS
 y/x = -(r + 2^ITERS) / 2^ITERS
```

The only zero-tag exceptional set is `y=-x (mod p)`, which is negligible for
random field inputs and fits the approximate-error model. Test
`dx_tagged_seed_recovers_division_with_negligible_exception` verifies this on
random samples.

Circuit validation is also wired behind env var `KAL_TAGGED_DIV_VALIDATE=1`:

- before pair1 Kaliski, set `ty := dy + dx`,
- compute tagged slope `-(λ+1)`,
- consume tagged `ty` to zero with the existing `pair1_mul2`,
- add known constant `1` to recover the ordinary `-λ` used by the remaining
  scaffold.

This default-off integration passes the full 9024-shot harness and 5 alt seeds
with:

```text
KAL_TAGGED_DIV_VALIDATE=1 cargo run --release -- --note tagged-div-validate
avg_toffoli = 4,138,926
qubits      = 2716
classical/phase/ancilla failures = 0
```

It is intentionally a validation path, not an optimization: it adds ~6k Toffoli
because it still uses the old Bennett-clean Kaliski and m_hist. Its importance
is that the tag algebra works in the real circuit.

A stronger default-off implementation was then wired behind
`KAL_TAGGED_DIV_COEFF_CHANNEL=1`: while ordinary Kaliski is running, it carries
an external coefficient pair `(lam, ty)` through the same branch controls. This
computes the tagged quotient directly and consumes `ty` to zero, removing
pair1's two schoolbook multiplications from the scaffold. It passes the real
harness, but it **invalidates the naive side-channel version as a SOTA path**:

```text
KAL_TAGGED_DIV_COEFF_CHANNEL=1 cargo run --release -- --note coeff-channel-div
avg_toffoli = 4,672,021
qubits      = 2,977
classical/phase/ancilla failures = 0
```

Why it loses:

- the external coefficient channel is live during Kaliski, adding a full
  256-qubit `lam` register at the Kaliski peak;
- each Kaliski iteration needs data-channel cswaps, a controlled modular add,
  and a modular double;
- the old inverse-state `(r,s)` and `m_hist` are still present just to clean
  qrisp branch flags and Bennett-uncompute the denominator state.

So the tag algebra is good, but the naive “parallel coefficient side channel”
is too wide and too expensive. The next reduction must remove the ordinary
inverse coefficient registers/history rather than run beside them.

A small positive result for that next scaffold is captured by
`stored_a_and_m_bits_recover_branch_pair`: if the final swap bit `a` is stored
alongside the existing `m_hist`, then `add = f & !(a xor m)` recovers the full
branch pair. This suggests a branch-only Kaliski generator can replace the full
ordinary `(r,s)` inverse sentinel with an `a_hist` bitstream. That saves only
~105 qubits net (`r,s` = 512 removed, `a_hist` = 407 added) if run interleaved,
and still stores history.

The more useful wired version is `KAL_TAGGED_DIV_BRANCH_STREAM=1`:

1. run denominator Kaliski while recording `m_hist`, `a_hist`, and `add_hist`,
2. free the known final denominator state `(u,v,f)=(1,0,0)`,
3. replay the branch histories into the tagged coefficient channel `(lam,ty)`,
4. uncompute the denominator histories.

This passes the real harness and is the first low-qubit tagged-DIV scaffold
under the current 2800q cap:

```text
KAL_TAGGED_DIV_BRANCH_STREAM=1 cargo run --release -- --note branch-stream-div
avg_toffoli = 4,729,076
qubits      = 2,763
classical/phase/ancilla failures = 0
peak phase  = br_rec_step2 / br_stream_coeff_add
```

A compressed variant is wired as `KAL_TAGGED_DIV_BRANCH_TERM=1`. It replaces the
full `add_hist` stream with a 9-bit terminal-iteration index. Coefficient replay
reconstructs active VG adds using `term_idx > i`, leaving only `m_hist+a_hist`
plus the tiny terminal register:

```text
KAL_TAGGED_DIV_BRANCH_TERM=1 cargo run --release -- --note branch-term-div
avg_toffoli = 5,267,537
qubits      = 2,714
classical/phase/ancilla failures = 0
peak phase  = pair2 Kaliski, not branch-DIV
```

This is a useful qubit-shape result: the tagged-DIV scaffold itself is now below
the current baseline peak, and pair2 Kaliski again dominates peak qubits. It is
also an invalidation of naive terminal-index replay as a Toffoli path: a 9-bit
comparator inside every coefficient iteration costs far too much.

The improved compressed variant `KAL_TAGGED_DIV_BRANCH_TERM_ROLL=1` keeps the
same `m_hist+a_hist+term_idx` qubit shape, but carries a rolling active flag
through coefficient replay. Each iteration only tests `term_idx == i` to toggle
the active flag, then uses one add-control `active & !(a xor m)`. This avoids
the double cmod-add and per-iteration greater-than comparator:

```text
KAL_TAGGED_DIV_BRANCH_TERM_ROLL=1 cargo run --release -- --note branch-term-roll-div
avg_toffoli = 4,733,146
qubits      = 2,714
classical/phase/ancilla failures = 0
peak phase  = pair2 Kaliski, not branch-DIV
```

This dominates both previous compressed versions in qubits and Toffoli, but it
is still not SOTA: branch-history recording plus coefficient replay remains
~600k Toffoli worse than the current default. The structural implication is
clear: the branch-DIV qubit shape is plausible, but history replay must be
removed or fused with the point-add body.

A follow-up tried to reuse the same branch-history idea as an exact compact
inversion for pair2 cleanup:

```text
KAL_PAIR2_BRANCH_INV_ROLL=1 cargo run --release -- --note pair2-branch-inv-roll
avg_toffoli = 5,957,442
qubits      = 3,147
classical/phase/ancilla failures = 0
peak phase  = shift22/cmod-add inside coefficient replay
```

This is a useful invalidation. Computing `inv_raw` by branch-record +
coefficient replay + inverse coefficient replay is both wider and much more
expensive than the existing full Kaliski state. The extra live objects are
`m_hist+a_hist+term_idx`, an explicit `(inv_raw, coeff_s)` pair, and the
Solinas cmod-add transients during replay. Therefore **branch-history replay is
not a general compact replacement for Kaliski inversion**; it only makes sense
when the coefficient channel becomes the output and is not reversed.

So these scaffolds prove clean reversible tagged DIV below 2800q, but not SOTA
Toffoli. The remaining gap is to eliminate/compress branch histories without
full replay and/or make the branch predicate self-cleaning.

A further algebraic invalidation is in
`bilinear_invariant_does_not_recover_inverse_branch`. The obvious preserved
relation

```text
r*v + s*u = 0 mod p
```

holds for almost all locally valid inverse-branch candidates, so it cannot be
the cheap self-cleaning predicate. This kills the simplest “try all inverse
branches and keep the one satisfying the invariant” route.

Therefore a self-cleaning DIV now needs:

- a **derived exact/near-exact predicate** over full `(u,v,r,s)` beyond the
  bilinear invariant, much cheaper than storing `m_hist`, and failing only on
  negligible tag-zero / collision subspaces; or
- a different update convention whose inverse branch is local.

Acceptance of a crude local classifier is not enough: >50% per-step blows up.
This is the next hard synthesis problem.

## 11. Fast invalidation tasks still open

1. **End-state branch predicate synthesis**: derive a reversible predicate for
   the previous branch from `(u',v',r',s',f', iter_idx)` cheap enough to replace
   `m_hist`. If it costs ~one comparator + a few modular half/add candidate
   checks per iter, forward-only DIV may still beat the 2-Kaliski scaffold.
   If it costs a full inverse-step replay or many n-bit comparisons, kill it.

2. **Direct DIV synthesis**: ignore current Kaliski structure and design a
   reversible Euclidean map for `(x,y)->(x,y/x)` where `y` is the coefficient
   register throughout and no independent quotient copy is made.

3. **Alternative Euclidean transform search**: seek an update convention whose
   coefficient matrix has form `[[*, k],[0, 1]]` or similar, so that backward
   naturally preserves the quotient rather than requiring another division.

4. **Cost if successful**:
   - one DIV/Kaliski-like invocation: target ~1.6M or less
   - delete `pair1_mul1`, `pair1_mul2`, second Kaliski: save ~1.7M
   - add coefficient modularity overhead in step4: likely +200-400k
   - add final coefficient rewrite: target <=300k
   - expected total if solved: **2.4-2.8M Toffoli**
   - qubits if folded into `ty` and history compressed: plausibly **1100-1500q**

This is now the main ground-up research direction alongside jumped Kaliski.

## 12. Strategy E — slope-coordinate affine permutation

After the BY denominator route missed budget, the cleanest non-BY ground-up
candidate is to make the line slope itself the live coordinate.  The algebra is
now captured in `single_inv_numeric.rs` as `replay_strategy_e_slope_coordinate`.
For

```text
dx = Px-Qx
dy = Py-Qy
m  = dy/dx
```

the point-add map can be written as the in-place-looking triangular schedule

```text
tx: dx -> Rx = m^2 - dx - 2Qx
ty: m  -> Ry = -m*(Rx-Qx) - Qy
```

`strategy_e_slope_coordinate_formula_passes_200` validates this exactly on 200
random secp256k1 point pairs.  This is a real algebraic reduction: `dx -> Rx`
is an involution once `m` is live, and no separate `lam` register is needed in
the final state.

The fast invalidation is equally sharp.  To make this a circuit we need both:

```text
DIV:   (x,y) -> (x,y/x)
IMUL:  (c,m) -> (c,-m*c-Qy)
```

with scratch cleaned and no copied quotient/product.  Known reversible ways to
implement `IMUL` are product-clean multiplication, which is equivalent to the
already-measured pair2 scaled inverse/product-clean primitive.  The budget test
`strategy_e_slope_coordinate_budget_requires_new_inplace_variable_multiply`
records the economics:

```text
current known product-clean route ≈ 2,988,510 Toffoli before safety margin
if IMUL were schoolbook-like       ≈ 2,022,750 Toffoli
needed IMUL saving                 ≈   965,760 Toffoli
```

So Strategy E is **validated algebraically** but **invalidated with current
primitives**.  It becomes SOTA-shaped only if a genuinely new in-place variable
multiply/divide primitive exists, roughly schoolbook-cost, phase-clean, and
without a raw inverse/history bank.  That primitive is more general than BY and
would also solve the earlier pair2 cleanup obstruction; without it, there is no
point wiring another affine scaffold around existing product-clean machinery.

### Attempt E1: destructive Montgomery IMUL

Next candidate for the missing primitive: destructively scan the multiplier bits
through a Montgomery add-and-halve accumulator:

```text
t = 0
for bit b_i of y:
    t += b_i * x
    if t odd: t += p
    t >>= 1
```

Forward algebra is good: after `n` steps, `t = x*y*2^-n mod p` up to one final
canonical subtraction.  If the consumed `y` bits were recoverable locally from
the post-window accumulator, this would give a product-clean multiply with only
one accumulator register and schoolbook-like cost.

`destructive_montgomery_product_is_algebraically_promising_but_not_locally_reversible`
now kills that hope on a small exact instance.  For `p=251`, `a=173`, and an
8-bit consumed window, the reachable poststate `t=223` has **512** valid
`(old_t, consumed_bits)` predecessors.  Thus the window cannot clear the
consumed multiplier bits from the accumulator alone.  A reversible circuit must
keep history/checkpoints or compute a nonlocal inverse of the multiplication,
which is exactly the product-clean obstruction we were trying to avoid.

Decision: destructive Montgomery is a useful failed primitive, not a route to
Strategy E.  A viable IMUL must use a different idea than local recovery from a
Montgomery accumulator.

### Attempt E2: MBUC product cleanup by phase-only quotient oracle

Another possible IMUL rescue is to compute `z=x*y`, measure the old `y` register
in the X basis, and apply only the MBUC phase correction instead of reversibly
recovering `y`.  For a measurement mask `s`, the needed phase is

```text
(-1)^(s · (z/x mod p))
```

as a boolean function of the preserved registers `(x,z)`.  If that phase oracle
were low-degree or sparse, product-clean multiplication might be much cheaper
than division.

`mbuc_product_cleanup_phase_oracle_is_not_low_degree_on_toy_field` kills the
cheap version.  On the 8-bit toy field `p=251`, a fixed mask already gives an
ANF with degree **15 of 16 variables** and density **32518 / 65536** monomials.
This does not prove every possible phase-oracle implementation is expensive,
but it rules out the hoped-for sparse/low-degree correction.  The phase-only
cleanup is just the quotient problem in disguise.

### Attempt E3: MBUC Montgomery quotient-history cleanup

A more structured variant keeps the Montgomery loop's internal quotient bits
`q_i` as the measured history instead of measuring the old multiplier directly.
In the `(x, old_y)` frame, this history is surprisingly sparse on small fields:

```text
n=8, p=251: degree=16, density=2440 / 65536
n=10, p=1021: degree=20, density=31684 / 1048576
```

But MBUC phase correction for an in-place multiply must be expressed in the
**output** frame `(x,z)`, because `old_y` has been replaced by the product.  The
test `montgomery_q_history_phase_in_output_frame_is_dense_dead` computes exactly
that frame transformation (`z = x*y*R^-1`) and gets:

```text
n=8, p=251: degree=16, density=31032 / 65536
```

So the structured q-history collapses back to a quotient-like dense phase when
viewed from the surviving registers.  This invalidates the Montgomery-history
MBUC rescue for Strategy E with current primitives.

### Attempt E4: top-level MBUC of the old affine point

A more global escape from the single-inversion wall would compute the new point
`R` out-of-place, X-measure the old input point `P`, and phase-correct from the
surviving output using `P = R - Q`.  This would avoid explicitly running the
second affine inverse if the phase oracle for point subtraction were cheap.

`top_level_mbuc_of_old_point_requires_dense_point_subtraction_phase` kills the
cheap version on the 8-bit toy curve `y^2=x^3+7 mod 251`: a fixed mask of the
old point bits, viewed as a function of `(R_x,R_y)`, has

```text
degree  = 15 / 16 variables
density = 19540 / 65536 monomials
```

So generic top-level MBUC just turns the affine reversibility wall into a dense
phase version of point subtraction.  It is not a free single-inversion cleanup.

Caveat checked: the full-domain ANF is pessimistic because a correct point-add
only supports curve points.  `curve_support_mbuc_phase_still_scales_not_constant_degree`
solves the support-restricted interpolation problem on toy curves.  The minimum
degree for the same old-point phase grows with `n`:

```text
n=4  p=13    min_degree=1
n=6  p=61    min_degree=3
n=8  p=251   min_degree=3
n=10 p=1021  min_degree=4
n=12 p=4093  min_degree=4
```

This follows the coding-theory dimension threshold, not a constant-degree
identity.  Extrapolating `sum_i<=d C(2n,i) >= ~2^n` puts a generic real-curve
extension near `d≈0.22n` (`≈56` for secp256k1), before sparsity/synthesis cost.
So support restriction does not resurrect generic top-level MBUC; only a highly
specialized sparse kickmix phase would be worth revisiting.

## 12. Attempt F: absorb Kaliski's scale by pre-scaling the denominator

Kaliski exposes a raw coefficient of the form

```text
inv_raw = -v^-1 * 2^iters  (mod p)
```

The ordinary point-add corrects this by applying `iters` halvings to the pair1
slope and `iters` doublings before the pair2 cleanup.  Those two correction
loops cost about **206k Toffoli** total, so a natural idea is to feed Kaliski a
scaled denominator

```text
v' = 2^iters * v
```

which makes the exposed raw inverse exact:

```text
-(v')^-1 * 2^iters = -v^-1 .
```

This validated algebraically and at full-circuit level for both inversion
sites when the prescaler used exact Cuccaro arithmetic:

```text
KAL_PRESCALE_PAIR1_SAFE=1
avg_toffoli = 4,786,373
qubits      = 2,972
clean       = yes

KAL_PRESCALE_PAIR2_SAFE=1
avg_toffoli = 4,771,009
qubits      = 2,969
clean       = yes
```

The all-exact result is much worse than default because the generic
phase-clean constant multiplier computes `2^407 mod p` by 183 modular
 doublings/halvings plus shifted adds, and it also keeps an extra scaled
denominator register live.  A mixed variant then isolated the phase culprit:
use exact q-q add/sub at the sparse constant bits, but keep fast modular
double/halve for the scale walk.

```text
KAL_PRESCALE_PAIR1_MIXED=1
avg_toffoli = 4,223,465
qubits      = 2,972
clean       = yes

KAL_PRESCALE_PAIR2_MIXED=1
avg_toffoli = 4,220,405
qubits      = 2,969
clean       = yes

KAL_PRESCALE_PAIR1_MIXED=1 KAL_PRESCALE_PAIR2_MIXED=1
avg_toffoli = 4,331,952
qubits      = 2,972
clean       = yes

KAL_PRESCALE_PAIR1_FOLDED=1
avg_toffoli = 4,223,465
qubits      = 2,969
clean       = yes

KAL_PRESCALE_PAIR2_FOLDED=1
avg_toffoli = 4,220,405
qubits      = 2,965
clean       = yes
```

The folded variants write the scaled denominator directly into Kaliski's `v_w`
initialization instead of carrying a separate scaled-denominator register; they
save only 3-4 peak qubits in the current allocator/phase profile but prove the
integration is algebraically equivalent.

Single-site mixed probes are only ~108-112k above the default exact path and ~550k below the
all-exact prescaler, proving the fast modular shifts are phase-safe here and the
earlier fast version failed because of the fast q-q add/sub in the constant
multiplier:

```text
KAL_PRESCALE_PAIR1=1
altseed_phase_batches_total = 1
```

Decision: scale absorption is a real algebraic lever, but not with the current
constant-multiply primitive.  It only becomes interesting if we implement a
secp256k1-specific phase-clean shifted-add prescaler for sparse constants like

```text
2^407 mod p = 2^151 * (2^32 + 977)
```

with total cost below roughly half a correction loop per side (≈50k Toffoli for
compute+uncompute), and preferably without an extra persistent n-bit copy.

## Direct controlled-constant Solinas correction note

New env-gated qubit tradeoff:

```text
KAL_DIRECT_CONST_HALVE=1
avg_toffoli = 4,121,014
qubits      = 2,715
clean       = yes
```

The helper `csub_nbit_const_direct_fast` avoids loading the sparse Solinas
constant `2^32+977` into a full 256-qubit `a` register for modular halve's
controlled subtract.  It uses a direct controlled-constant borrow ripple and
measurement-uncomputes the borrow chain.  This removes the previous
`bk_step6_7_8` peak and moves the peak to `bk_step4`.

It is not a Toffoli win (default remains `4,111,918 @ 2716q`), but it is a much
cheaper low-qubit lever than the older dirty-venting halve attempt.  The add
analogue is a valid primitive in standalone tests, but phase-cliff sensitive in
the full harness:

```text
KAL_DIRECT_CONST_DOUBLE=1
altseed_classical_total = 1
altseed_phase_batches_total = 2

KAL_DIRECT_CONST_DOUBLE=1 KAL_BULK3_ITERS=370
avg_toffoli = 4,121,506
qubits      = 2,716
clean       = yes
```

So direct cadd is useful only as an env-gated tested tool; it does not improve
the default Toffoli/qubit point.  Combining both direct cadd and csub under the
conservative prefix is clean at `4,130,602 @ 2715q` with `29,250,534` emitted
ops; this is a low-emitted-op/qubit tradeoff, not a SOTA path.

## R-small no-correction threshold update

The Kaliski `r` smallness shortcut now defaults to 257 iterations instead of
255:

```text
KAL_R_SMALL_THRESHOLD=256  -> 4,110,898 @2716q clean
KAL_R_SMALL_THRESHOLD=257  -> 4,109,878 @2716q clean
KAL_R_SMALL_THRESHOLD=258  -> phase failure
```

`258` also fails with conservative `KAL_BULK3_ITERS=370`, so the clean cliff is
real in this scaffold.  Default exact is now `4,109,878 @2716q`; this is only a
small local gain, not a structural route.
