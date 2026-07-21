# Measured frontier leads (2026-06-07, build_circuit traces)

Supersedes the speculative parts of `2026-06-07-structural-breakthrough-leads.md`.
**Correction:** that note claimed round84 was the peak binder. It is not — measured
peak is the GCD walk (`reverse_add`/`shift`); round84 sits 18 q below peak.

## How these numbers were taken (measured, not estimated)

Two `build_circuit` runs on the configured tier-3 route (no extra env beyond
`configure_ecdsafail_submission_route`):
- `TRACE_PHASES=1` → emitted CCX per phase. **Total emitted CCX = 1,456,963**,
  total ops = 9,767,086. This equals the scored avg executed Toffoli exactly
  (no classically-conditioned CCX in this route), so emitted CCX = score numerator.
- `POINT_ADD_COUNT_ONLY=1 TRACE_PHASE_ACTIVE=1` → per-phase live-qubit maxima.
  **Peak = 1302.**

Per the user's request, no further runs were taken. Anything I could not derive
from these two traces + the source is marked **[needs run]** with the exact probe.

## Measured peak floor map (live qubits per phase)

| phase | active_q |
|---|---|
| `compressed_block_tobitvector_reverse_add` | **1302** ← binder |
| `compressed_block_tobitvector_shift` | **1302** ← binder (idle scratch, see B) |
| `compressed_block_tobitvector_compress_block` | 1301 |
| `compressed_block_tobitvector_reverse_cswap` (impl.) / `apply_chunk_{add,sub}_final_ripple` | 1299 |
| `raw_pa_x_restore`, `round84_fused_square_xtail_add_double_ox` | 1285 |
| `round84_inplace_solinas_square_{forward,inverse}` | 1284 |

Consequence: the next 3 q of peak (1302→1299) are **GCD-walk-only**. Below 1299 you
must *also* cut `compress_block` (1301) and the chunked apply (1299). round84 (1284)
is irrelevant to score until peak drops below 1284.

## Measured Toffoli by category (emitted CCX, = score numerator)

| category | CCX | % | width basis |
|---|---|---|---|
| apply mod add/sub (chunked, both GCDs) | 363,780 | 25.0 | full 256 + chunk boundary clears |
| GCD body sub/add (`materialized_*_{load,body}`) | 275,016 | 18.9 | active_width (band-trimmed) |
| **cswap total** | **271,744** | **18.7** | tobitvector 139,648 (active_width) + apply 132,096 (**full 256**) |
| apply `double_y`+`halve_y` (K2 2nd double/halve, mod p) | 168,216 | 11.5 | full mod-p Solinas |
| tobitvector `shift`+`unshift` (K2 2nd shift) | 139,648 | 9.6 | active_width |
| round84 square fwd+inv | 131,582 | 9.0 | (peak 1284, slack) |
| branch_bits fwd+rev | 40,752 | 2.8 | compare_bits schedule |

Within the apply mod add/sub: `boundary_clear` = 99,588 (6.8%) is pure chunking
overhead (value-exact for any cut), the rest (264k) is the real per-step y+=x mod p.

---

# Leads, in the four requested areas

Format per lead: **files/fns · env · ΔT · Δpeak/phase · correctness · island**.

## A. `DIALOG_GCD_SHIFT_BAND_TRIMS` — small Toffoli knob, **peak-neutral**

1. **Files/fns:** `dialog_gcd_shift_band_trim` (rounds/dialog/mod.rs:208),
   `dialog_gcd_k2_shift_active_width` (mod.rs:229); consumers
   compressed.rs:809-824 (forward 2nd shift) and :938-953 (reverse un-shift).
2. **Env:** `DIALOG_GCD_SHIFT_BAND_TRIMS` (currently unset = OFF). Accepts a
   per-band list OR the literal `body` (reuses `BODY_CARRY_BAND_TRIMS`).
3. **ΔT:** the trimmed phases total 139,648 CCX. Each band trims `w` bits off the
   `(k2_shift_active_width-1)` cswap cascade → saves `Σ_step w(step)` per
   (GCD,direction), ×2 directions ×2 GCDs. With `=body` (schedule
   `0,3,3,3,3,3,1×17,3,3,3`, band_size 10): Σw ≈ 404 ⇒ **≈ −1,600 CCX (−0.11 %)**.
   Beyond-body costs ~1,032 CCX per extra bit-of-trim across all 258×4 slots.
4. **Δpeak:** **0.** The shift phase peaks (1302) because the per-step composite
   scratch is still live (freed only at step end), *not* because of shift width.
   Narrowing the shift does **not** touch peak. → This knob is the wrong tool for
   the "reduce reverse_add/shift peak" goal; it is a pure (small) T lever.
5. **Correctness:** untested. Value-exact requires realizable bitlen ≤
   `aw − w − 1` at each trimmed step (one bit tighter than the body trim, because
   the truncated cascade leaves `v[aw-w-1]` unshifted where the true shift would
   zero it). So `=body` is **1 bit too aggressive** at the boundary.
6. **Island:** `=body` adds a thin hazard class (inputs with realizable bitlen
   exactly `aw−w` at a trimmed step) → current nonce 11201395269 may not survive;
   cheap re-hunt. **Island-free variant:** schedule = `body − 1` per band (floored
   at 0), e.g. `0,2,2,2,2,2,0×17,2,2,2` ⇒ shares the body premise exactly, keeps
   the nonce, but only ≈ −850 CCX. Honest verdict: real but ≤0.1 %.

## B. Reduce the GCD-walk peak (`reverse_add`/`shift`) — the only score-multiplier lever

The peak is `u(256)+v(256)+compressed_log+raw_block+owned`, where
`owned = b.alloc_qubits(want − borrowed)` in `dialog_gcd_build_composite_scratch`
(compressed.rs:352-468), `want = 2·body_len − 1`. The binder is the step whose
`owned` is largest. round84/apply are below, so −1 q here = −1 q global =
−1,456,963 score (≈ 0.077 %/q; one q ≈ 13 of the recent nonce submissions).

**B1 — binder notch (the proven mechanism, already wired):**
1. **Files/fns:** `dialog_gcd_binder_notch_steps`/`_extra` (mod.rs:262,273) feed
   `dialog_gcd_body_carry_trunc_width` (mod.rs:256-258) → shrinks `body_w` →
   `body_len` → `want` → `owned` at the listed steps. (The existing
   `trio_width_notch` step 11 extra 2 is the same trick, already on.)
2. **Env:** `DIALOG_GCD_BINDER_NOTCH_STEPS=<binder step>`, `DIALOG_GCD_BINDER_NOTCH_EXTRA=1`.
3. **ΔT:** ≈ −2 CCX per notched step per GCD-pass it touches (body+cswap-if-also-trimmed). Negligible.
4. **Δpeak:** −`EXTRA` at the binder step **iff** that step is the unique
   `max(owned)` step. **[needs run]** `PROBE_SCRATCH=1` (compressed.rs:454, prints
   `owned` for active_width ≥ 254) → take the `step` with max `owned`; that is the
   binder. Then notch it by 1 and re-`TRACE_PHASE_ACTIVE` to confirm 1302→1301.
   If two steps tie at max, notch both.
5. **Correctness:** untested. Value-exact on reachable support (top `EXTRA` extra
   bits of u,v are |0> by the same realizable-bitlen bound the body trim uses),
   identical hazard *kind* to the accepted route.
6. **Island:** deeper trim at one step ⇒ new straggler inputs ⇒ **needs a fresh
   `DIALOG_TAIL_NONCE`** (re-hunt with the GPU/CPU GCD prefilter, then 9024 eval).
   Density comparable to the existing band-trim islands (~1/108 of GCD-survivors
   per the 2026-06-06 note), so tractable.

**B2 — one more borrow lane (island-free if it lands):**
1. **Files/fns:** the borrow sources in `build_composite_scratch` (future-log,
   current-block cells, `v[aw..]`, `u[aw..]`, current `s2`, sibling `s2`). At the
   binder step (early, wide `aw`) `u[aw..]`/`v[aw..]` are nearly empty, so `owned`
   is the deficit vs the future-log runway.
2. **Env:** none new — needs code: an additional provably-|0> idle source folded
   into `push(...)`. Candidates to check at the binder step: the *next* block's
   not-yet-written compressed cells (beyond current-block), or `b0`/`b0_and_b1`
   raw cells of already-compressed earlier slots in the same block.
3. **ΔT:** 0 (pure relabel).
4. **Δpeak:** −1 if it converts one `owned` lane to borrowed at the binder step. **[needs run]** PROBE_SCRATCH to confirm a clean |0> lane exists there.
5. **Correctness:** value-exact-always if the borrowed cell is provably |0> across
   the step window and restored (it is, by the measured uncompute) — then **no FS
   hazard at all**, nonce 11201395269 survives.
6. **Island:** unchanged (island-free) — this is the preferred peak cut if a lane exists.

**B3 — note:** freeing `composite_scratch.owned` *before* the forward `shift`
(it is idle there, compressed.rs:826) drops the forward `shift` phase off 1302 but
**not** the global peak, because `reverse_add` (compressed.rs:956-985) genuinely
needs the scratch. So B3 alone = 0 score; only B1/B2 move the global peak.

## C. Partial cswap reduction (271,744 CCX = 18.7 %)

**C1 — tobitvector cswap band-trim (island-free, small):**
1. **Files/fns:** the cswap loops compressed.rs:796-802 (fwd) and :987-993 (rev)
   run at **full active_width**; the body sub/add beside them already trims to
   `aw − body_trim` (mod.rs:188 explicitly leaves "cswap and comparator at full
   active_width"). Add a width clamp = `dialog_gcd_body_carry_trunc_width(aw,step)`
   to the cswap loop bound.
2. **Env:** reuse `DIALOG_GCD_BODY_CARRY_BAND_TRIMS` (no new flag) so the cut ≤
   the body's own assumption.
3. **ΔT:** −Σ body_trim per (GCD,dir) on the 139,648 tobitvector-cswap CCX ≈
   **−1,600 CCX (−0.11 %)**.
4. **Δpeak:** 0 (cswap is in-place Fredkin, no scratch).
5. **Correctness:** untested but **value-exact by the body trim's own premise**
   (the swapped high bits are the bits the body already assumes are |0>; swapping
   |0>↔|0> is identity).
6. **Island:** **island-free** — same premise as the accepted body trim, nonce
   11201395269 survives. This is the one cswap cut that is genuinely free; take it.

**C2 — apply cswap (132,096 CCX) is NOT trimmable:** compressed.rs:1106-1108 and
:1146-1148 swap the full 256-bit residues x↔y (Montgomery accumulator pair); they
are not bit-length-bounded, so there is no value-exact truncation. Any reduction
here needs the swap *fused into* the `cadd`/`csub` (a real redesign of
`apply_bitvector`), which is **research, not a measured frontier knob** — flagged,
not claimed.

## D. Low-qubit round84 square — **no score lever at current peak**

1. **Files/fns:** `round84_emit_fused_square_xtail` (rounds/dialog/mod.rs:14) →
   `squaring_sub_from_acc_schoolbook_lowq_shift22` under `ROUND84_INPLACE_SOLINAS_FOLD`.
2. **Env:** `ROUND84_INPLACE_SOLINAS_FOLD=1` (on), `ROUND84_XTAIL_KARATSUBA`,
   `ROUND84_XTAIL_WALK_SQUARE`.
3. **ΔT:** the known faster square (Karatsuba) is −16,272 CCX **but** needs the
   `z1_reg` (~258 q).
4. **Δpeak:** round84 is at **1284 (18 q slack)**. Karatsuba's +258 q ⇒ 1284+258 ≫
   1302 ⇒ becomes the binder. Does **not** fit the 18 q headroom. There is no known
   square variant that trades **≤18 q** for a Toffoli cut (the in-place fold is
   already the low-q schoolbook; the symmetric/fast variants are T-identical and
   only differ in carry-lane hosting).
5. **Correctness:** n/a — nothing to change.
6. **Verdict:** round84 cannot improve score until the GCD-walk peak drops below
   1284. Deprioritize, exactly as the user's "only if peak stays under 1302" gate
   implies. (If a sub-18-q-overhead square speedup is ever found it would be a free
   −T, since round84 has the headroom — but none is known.)

---

## Honest bottom line

On the current frontier the measurable, value-exact knobs are all small:
- **C1 (tobitvector cswap trim to body width):** ≈ −1,600 CCX, peak-neutral,
  **island-free, keep nonce.** Take it first — zero risk.
- **A island-free shift trim:** ≈ −850 CCX, peak-neutral, island-free.
- **B1 binder notch −1 q:** ≈ −1,456,963 score (the biggest single move), but
  **[needs run]**: PROBE_SCRATCH to find the binder step, then a nonce re-hunt.
- **B2 extra borrow lane −1 q:** same score move, **island-free if a |0> lane
  exists** at the binder step — strictly better than B1 if it lands.

Everything ≥1 % (apply mod-add 25 %, GCD body 18.9 %, apply cswap, double/halve)
is bound to a real redesign, not an env knob, and is out of scope for
"measured current-frontier." The combined safe set (C1 + A-island-free + one
peak q via B2) is ≈ −2,450 CCX **and** −1 q ⇒
1301 × 1,454,500 ≈ 1,892,304,500, beating 1,896,965,826 by ≈ 4.66 M (0.25 %),
**without changing the Fiat-Shamir island** if B2 lands clean. Confirm each Δpeak
with `TRACE_PHASE_ACTIVE` and each ΔT with `TRACE_PHASES` before submitting; run
the 9024 eval for the final stack.

---

# Redesign assessment: is there a ≥1 % structural move?

Honest read after walking every big category. Toffoli is ~80 % the two binary-GCD
inversions (tobitvector + apply); the square is 9 %, the rest small. So a ≥1 % move
must make the **inversion** cheaper or cut the **GCD-walk peak**. Knob-level is
exhausted; the candidates below are real rewrites.

## Bet 1 (highest upside, highest risk): implicit-shift δ divstep (Bernstein–Yang)

**Target:** the *physical shift* tax. Measured: tobitvector K2 2nd-shift
(`shift`+`unshift`) = 139,648 (9.6 %); apply `double_y`+`halve_y` = 168,216
(11.5 %) — and `fused_double_y` (compressed.rs:2073) is mostly the two shift
cascades (the conditional 2nd shift is ~256 cswaps, lines 2089-2091) + one fold.
So **~21 % of all Toffoli is spent physically shifting v and re-doubling y every
step.** A Bernstein–Yang `divstep` tracks the relative shift in a small `δ` counter
and never physically shifts — the halving is implicit; the apply mirrors it with a
δ-indexed access instead of a mod-p doubling.

**Why it's the right reference:** `configure_ecdsafail_submission_route` already
cites **Gidney et al. arXiv:2510.10967** for its *width bound only*
(mod.rs:1274). That is a reversible safegcd/BY paper; its **divstep + apply
construction** is exactly the implicit-shift machinery. Read the paper's circuit,
not its inequality.

**Why it might fail / honest caveat:** BY divsteps still do a per-step full-width
conditional add/sub on (u,v) — that work does **not** vanish, only the shifts do.
And BY still has a conditional swap (the δ>0 branch swaps f,g), so this does **not**
remove the cswap (18.7 %); it removes the shift/double layers (~21 %). Net win only
if the δ bookkeeping + implicit-shift apply is cheaper than the shift cascades it
deletes. Plausibly ≥10 %, but unproven until built. It is a ground-up rewrite of
`tobitvector_steps`, `apply_bitvector`, the transcript format, AND the width
envelope. Prototype the divstep in isolation against
`d1_inplace_*_lowerer_component_stats_are_pinned` (mod.rs:838) before touching the
PA. Correctness: untested. Island: new transcript ⇒ full nonce re-hunt.

## Bet 2 (moderate, lower risk): branch bit from the body-subtract carry

**Target:** `branch_bits` = 40,752 (2.8 %). The divstep computes `b1 = (u>v)` with a
standalone truncated comparator (compressed.rs:755-794), then cswaps, then subtracts
`v-=u`. The subtract's borrow-out **is** `(v<u)`. Reorder to take the branch from
the subtract's own carry instead of a separate comparator.

**Why it might fail / honest caveat:** the standalone comparator is *truncated* to
`compare_bits` (avg ~40, scheduled down to 5), which is far narrower than the
`active_width` (~130) subtract. So the existing comparator is already cheaper per
step than reading a full-width subtract carry, and the swap must precede the sub
(needs the bit first). You only win if you can derive the branch from a
`compare_bits`-truncated speculative difference and fold the fix-up into the cswap.
Net likely ~1 %, not 3 %. Medium fiddle; same island if the truncation widths match.

## Ruled out — don't spend time here (with the reason)

- **1 inversion / dx³**: blocked, proven (ONE_INV_DX3 blocker; recovering dx after
  (tx,ty)→(Rx,Ry) needs inverting Rx−Qx = a 2nd inversion). Don't retry.
- **Montgomery batch-invert the two GCDs**: data-dependent — `c=Qx−Rx` only exists
  *after* Rx, which needs the 1st inverse. Can't form the product to batch. Blocked.
- **Eliminate the transcript log (fuse GCD+apply, no storage)**: the log is what
  carries info from the GCD to the target so the GCD registers can be uncomputed
  while the target keeps the result. Fusing would un-compute the result. The log is
  near its entropy floor already (round763 6→5 + K2-pair). Structural, not removable.
- **K=3+ bounded shift**: each higher radix adds a full conditional-shift layer
  (tobitvector) + conditional-double layer (apply) that costs ≈ what the fewer steps
  save (only ~12 % of steps strip a 3rd zero). Net ~neutral. K2 is near-optimal.
- **Batched "jump" divstep (apply the k-step transition matrix at once)**: a Toffoli
  k-bit × n-bit multiply is k·n, done n/k times = n² — same as per-step. The matrix
  entries grow with k, so no Toffoli win in this cost model (it's a software
  locality win only).
- **Fuse the apply double's reduction into the cadd's reduction**: the conditional
  Solinas subtracts are conserved (2y+x ∈ [0,3p) still needs 2 reductions = double's
  1 + add's 1). No saving; `FUSED_FOLD` already captured the carry-chain share.
- **Windowed/table apply (Häner et al.)**: the apply multiplicands (x,y) are quantum,
  not classical, so there is no classical table to look up. N/A here.
- **Lazy/deferred mod-p reduction in the apply**: registers are 256-bit and p≈2²⁵⁶,
  so there is ~0 spare headroom; each per-step double can add a bit, forcing a
  reduction every step. Deferring needs wider y ⇒ raises the apply phase toward the
  1302 binder. Peak-blocked.

**Bottom line:** the only ≥1 % that isn't blocked is Bet 1 (implicit-shift BY
divstep per 2510.10967), and it is a real research rewrite with an unproven net.
Bet 2 is a safer ~1 %. Everything else is sub-1 % knobs (section A–C) or blocked.
