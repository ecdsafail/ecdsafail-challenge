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

It is still **not** a SOTA path as-is: Toffoli is worse than baseline because it
stores 3 history streams and replays a full modular coefficient pass, with no
Kaliski truncation in the branch generator. But it proves the y+x tagged DIV can
be made into a clean reversible point-add sub-scaffold below 2800q. The
remaining gap is to eliminate/compress branch histories and/or make the branch
predicate self-cleaning.

Therefore a self-cleaning DIV now needs:

- a **derived exact/near-exact predicate** over full `(u,v,r,s)` that is much
  cheaper than storing `m_hist`, and can fail only on the negligible tag-zero
  subspace; or
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
