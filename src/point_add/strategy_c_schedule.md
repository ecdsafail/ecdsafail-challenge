# Strategy C schedule: can we fit it?

## Discipline
No speculation. Every register allocation + free is accounted for. Every
claim is checked against `src/point_add/single_inv_numeric.rs::replay_strategy_c`
and against actual qubit tallies from `mod.rs` primitives. If a claim
fails the arithmetic, the row is marked **BLOCKED** — no hand-waving.

## Target
Classical formula (from `replay_strategy_c`, 200/200 passing):
```
dx   = Px - Qx                              (in tx)
dy   = Py - Qy                              (in ty)
dx²  = dx·dx
dx³  = dx²·dx             (= w)
dy²  = dy·dy
v    = dy² - dx²·(Px+Qx)
w⁻¹  = Kaliski(w)

Rx = v · (dx·w⁻¹)                          (→ tx final)
Ry = (dy·(dx²·Qx - v) - w·Qy) · w⁻¹        (→ ty final)
```

Classical fact: one inversion only. Two multiplication "clusters":
Rx-cluster (needs `v, dx, w⁻¹`) and Ry-cluster (needs `dy, dx², v,
w, w⁻¹, Qx, Qy`).

## Primitive qubit peaks (from microbench.md)

| primitive                                 | transient peak |
|-------------------------------------------|---------------:|
| schoolbook mul (n×n→2n)                   | 1797           |
| karatsuba-1 mul                            | 2055           |
| karatsuba-2 mul                            | 2315           |
| schoolbook_addsub fwd (fast Cuccaro)       | 1283           |
| mod_add_qq_fast                            | ~n + n = 512 (?) — need to verify |
| Kaliski forward body peak                   | ~2204 (our baseline incl tx,ty,lam) |

Each mul subroutine allocates its own temps and frees them.
So the "size on top of current" for a mul depends on what's already
alive.

## Inventory during Strategy C

Persistent registers that DO NOT get freed across the body:
- tx (256)
- ty (256)
- m_hist (407)

Registers allocated for the strategy:
- `w_reg` (256): holds dx³, passed as Kaliski `v_in`. Alive from
  pre-Kaliski through post-Kaliski (mostly, until uncomputed). This is
  an extra register not present in current design (current design uses
  `tx` directly as `v_in`).
- `r_low` ⊆ Kaliski state (256): Kaliski output. Alive after forward,
  until backward consumes it.
- Kaliski internal state during fwd/bwd: u, v_w, s, f_flag = 1025q
  (but FREED during body closure, re-alloc'd for backward).

Registers we'll need to allocate AT SOME POINT during strategy C body:
- `dx²`, `dx³`, `dy²`, `v`, and intermediates for Rx/Ry.

## Schedule proposal #1: "build everything, use everything, uncompute everything"

Phases labeled by rough peak qubit count ("p"). Arrows show where a
register is allocated ("+") vs freed ("-").

### Phase 1: pre-Kaliski setup

Goal: compute `w = dx³` into a fresh register while keeping `dx, dy`
live in tx, ty.

  Step | action                                     | live registers                    | peak
  -----|--------------------------------------------|-----------------------------------|-----
   1.1 | alloc `dx2_reg` (256)                      | tx, ty, dx2_reg                   | 768+
   1.2 | compute `dx2_reg := tx·tx`                 | tx, ty, dx2_reg + sb mul internal | 768+1797=2565 peak transient, 768 persistent
   1.3 | alloc `w_reg` (256)                        | tx, ty, dx2_reg, w_reg            | 1024
   1.4 | compute `w_reg := dx2_reg·tx` (→ dx³)      | same + mul internal               | 1024+1797=2821  **BLOCKED** over 2800

**Step 1.4 is right at the cap.** We'd have 1024 persistent + mul's
1797-qubit transient = 2821. Close but just over.

**Mitigation**: use a lower-peak multiplier. If `mul_add_into_acc` only
peaks at ~1283 (schoolbook_addsub fast variant) *relative to the base*,
then 1024 + 1283 = 2307. Need to check: does our `mod_mul_write_into_zero_acc_schoolbook` actually peak at 1797
GLOBALLY, or does it peak at current+1797? The microbench measured
starting from a fresh `B`, so that's global peak-from-empty, which
equals our current state (0) + transient. In our live circuit, existing
persistent qubits ADD to the mul's transient.

Let me be more careful: schoolbook_mul's "1797 peak" is the mul called
in isolation on fresh registers. The mul allocates:
- tmp_ext (2n = 512) — output accumulator
- some Cuccaro carries + ancillas (~n = 256)
- input copies (x_ext, y_ext, pads) (~n)
Total transient during mul: ~4n = 1024 qubits alloc'd by the mul.
Plus the 3 inputs (x, y, acc) at 3n = 768.
Grand total: ~5n = 1280 qubits. Matches microbench 1283 for schoolbook_addsub_fwd.

So mul ADDS ~4n = 1024 qubits on top of its inputs. That's consistent.

In Phase 1.4 context: inputs are dx2_reg (already allocated), tx (already
allocated), w_reg (destination, already allocated). Mul adds 4n = 1024
transient. Pre-mul live: tx(256)+ty(256)+dx2_reg(256)+w_reg(256) = 1024.
Mul transient: +1024. Peak: 2048. **UNDER cap by 750q.** ✓

### Phase 2: Kaliski on w

  Step | action                                     | live                              | peak
  -----|--------------------------------------------|-----------------------------------|-----
   2.1 | Kaliski forward on w_reg                   | tx, ty, dx2_reg, w_reg + kal_state | ~1024+1432+fwd transient
   2.2 | free u, v_w, s, f_flag (KAL_FREE_S)        | tx, ty, dx2_reg, w_reg, r_low, m_hist |  ~1430
   2.3 | body runs here (Phase 3)                   |                                    |
   2.4 | re-alloc u, v_w, s, f_flag                 |                                    |
   2.5 | Kaliski backward on w_reg                   | tx, ty, dx2_reg, w_reg + kal_state | ~1024+1432+bwd transient

**Peak during Kaliski forward**: 1024 (pre-Kaliski) + 1432 (full kal
state incl m_hist) + ~260 (step 4 transient) = ~2716.
Same as current! Because current circuit has similar structure (tx+ty+lam=768+lam_256=1024 persistent inside body + kal state). Actually current has lam NOT dx2 during Kaliski. Either way, peak = ~2716.

### Phase 3: body between forward and backward

Goal inside body (u, v_w, s, f_flag freed, r_low and m_hist alive):
- Compute Rx into tx, Ry into ty, uncompute all intermediates except
  Kaliski state.

At body entry:
- persistent: tx(256), ty(256), dx2_reg(256), w_reg(256), r_low(256), m_hist(407)
  = **1687 qubits alive**
- budget under 2800: **1113 qubits** for transients + intermediates.

#### Phase 3 sub-plan

  Step | action                                     | new alloc    | live delta
  -----|--------------------------------------------|--------------|----------
   3.1 | alloc `dy2_reg` (256)                      | +256         | 1687+256=1943
   3.2 | compute `dy2_reg := ty·ty` (squaring)      | ~1024 mul transient | 1943+1024=2967  **BLOCKED** over 2800

Same pattern. Mul transient is 1024, pre-mul is 1943 → 2967, over cap.

**Mitigation**: use `KAL_FREE_S=1` style tricks to free MORE Kaliski
state during body. `m_hist` is 407q. Is there a way to free it during
body? `m_hist` is needed during backward Kaliski. So no.

What about `r_low`? It holds `w⁻¹`, which is what we need to multiply v
against later. Can't free it.

What about `dx2_reg`? We need it for the Ry formula (core = dx²·Qx − v).
Can we re-derive dx² from live state after Kaliski backward without
keeping it alive? If we need it inside the body AND dx is live in tx,
yes: we could COMPUTE dx² later instead of keeping it alive. But then
we need dx² for BOTH the Rx formula (indirectly through v) and the Ry
formula.

Wait — does the Rx computation need dx²? Rx = v · (dx · w⁻¹). So Rx
only needs v, dx, w⁻¹. Not dx² directly.

But `v = dy² - dx²·(Px+Qx)`, computed earlier. Once v is in a register,
dx² can be freed.

Let me restructure the schedule so dx² is freed AS SOON AS possible.

#### Phase 3 schedule revised

Goal: minimize simultaneous live intermediates.

  Step | action                                     | peak after         | live after
  -----|--------------------------------------------|--------------------|----------
   3.1 | alloc `v_reg` (256), `dy2_reg` (256)       |                    | 1687+512=2199
   3.2 | compute `dy2_reg := ty·ty`                 | 2199+1024=3223  **BLOCKED** | 2199
  
Still fails due to mul peak of +1024 on top of ~2000 persistent.

The fundamental issue is: **mul transients + 1432 Kaliski state +
tx/ty/output-cluster = over 2800**.

### Mul transient analysis

What does a mul ACTUALLY alloc? Let me look:

mod_mul_write_into_zero_acc_schoolbook (n=256):
- tmp_ext (2n = 512)
- inside schoolbook_mul_into_addsub:
  - low (1)
  - pad for y_ext (1)
  - c_in (1)
  - pad for x_ext ext (1)
  - some carry ancillas (depending on add primitive) (~n)
Total: ~2n + (4) + n = ~3n.

Actually looking at the real schoolbook_mul_into_addsub code, it uses
multiple temporary pads and carry ancillas of ~n bits during adders.
Net transient over the mul call: **~2n to 3n qubits**.

For n=256: **512 to 768 transient**.

### Revised peak accounting for Phase 3 sub-plan

If mul transient is ~3n = 768, then Phase 3.1 peak:
- pre-mul live: tx+ty+dx2+w_reg+r_low+m_hist+dy2_reg = 1687+256=1943
- +mul transient 768: peak = 2711.  **UNDER cap.** ✓

If mul transient is 2n = 512:
- peak = 2455. ✓ big margin.

The microbench 1797 is tighter than 768+inputs. Let me reconcile.

**Microbench 1797** = peak observed when schoolbook mul is called FROM
a fresh B, so the 1797 INCLUDES the three input registers (3n=768) plus
its own alloc (2n=512 tmp_ext + some small carries). 768+512+~20 = ~1300.

Hmm but microbench says 1797. Let me verify more carefully.
Actually the microbench likely includes tmp_ext (2n=512) + addend (n=256) +
x (n=256) + y (n=256) + internal pads/carries (~20) = 1300.
1797 - 1300 = 497 more... maybe internal carry regs are longer than I think.

Either way: **mul transient in addition to its 3 input registers is
roughly 2n to 3n (512–768 qubits)**. I'll use 768 as upper bound.

### Returning to Phase 3 budget

At body entry: tx(256)+ty(256)+dx2_reg(256)+w_reg(256)+r_low(256)+m_hist(407) = 1687.

Sub-budget 3.1: alloc dy2_reg(256) → 1943. Peak during mul to compute
dy2_reg = 1943+768 = 2711. ✓ UNDER 2800 cap.

Sub-budget 3.2: compute `v_reg = dy2_reg - dx2_reg·(Px+Qx)`. This is:
  - in-place `dy2_reg -= dx2·(Px+Qx)` to save one register. 
  - dx2·(Px+Qx) is quantum × classical. Cost ~n² schoolbook but with
    classical operand → fewer Toffolis (each bit of classical c just
    triggers a controlled add). 
  - But still needs a mul-like subroutine.
  
After mul-sub: dy2_reg holds v. Peak during mul ≈ 1943+768=2711. ✓

Sub-budget 3.3: now v is in dy2_reg. Can we FREE dx2_reg yet? No, it's
needed for Ry's `dx²·Qx - v`.

live: tx+ty+dx2+w_reg+r_low+m_hist+v_reg = 1687+0 (v_reg overwrote dy2_reg, same slot) = 1687.

Sub-budget 3.4: alloc `Rx_tmp = dx · w⁻¹`. That's tx·r_low → new n-bit.
  - alloc Rx_tmp (+256) → 1943
  - schoolbook mul Rx_tmp := tx · r_low. Peak: 1943+768=2711. ✓
  
live: 1943.

Sub-budget 3.5: alloc `Rx_reg`. Compute `Rx_reg := v · Rx_tmp`.
  - alloc Rx_reg (+256) → 2199
  - mul peak: 2199+768=2967.  **BLOCKED** over 2800

**Mitigation**: free `Rx_tmp` before allocating Rx_reg? No, Rx_tmp is
the operand for Rx_reg's computation. Must be alive during the mul.

Free dx2_reg first? dx2_reg is still needed for Ry's formula. So no.

**Alternative**: put Rx into tx directly. tx currently holds dx. If we
swap or overwrite tx with Rx, we lose dx. But we'll need dx to
un-compute the entire chain at the end (since tx was the input
register).

ACTUALLY — wait. Let me re-think. After Phase 3 completes, we want tx
to hold Rx and ty to hold Ry. So overwriting tx with Rx is DESIRED.
The issue is that Rx's formula includes dx (via `dx·w⁻¹`), so dx must
be live UNTIL the moment we write Rx.

Plan: compute Rx_tmp = dx·w⁻¹ into a fresh register (as above, uses tx=dx), then write `Rx` directly INTO tx: first ZERO tx (uncomputing dx), then XOR v·Rx_tmp into tx. Can we un-compute dx from tx cleanly? We computed dx = tx - ox (mod p) at circuit start. We can un-compute by adding ox back: `tx += ox mod p` would turn tx from dx to Px. But Px is quantum and is "the input" to the circuit — we don't have Px, we had tx = Px at entry which became dx.

Actually classical-bit register ox holds Qx. So `tx += ox mod p` turns tx from dx=Px-Qx back to Px. We DO have that.

OK the zeroing path is actually: if we can't just overwrite, we can do:
1. After computing Rx_tmp, un-compute dx in tx: `tx += ox` → tx = Px.
2. Now tx is NOT dx anymore, so subsequent ops using "dx" must use a
   different source (we can keep dx_reg alive separately).

Hmm this is getting hairy. Let me step back.

### Key realization

The fundamental tension: **the peak inside the Kaliski body is
bounded by (Kaliski state + tx + ty + other-alive-intermediates + mul
transient)**. Currently our 2-Kaliski scaffold has (r+m_hist + tx + ty +
lam + mul transient) = 256+407+256+256+256+260 ≈ 1700. Peak ~2716 includes
step-4 transient during backward Kaliski.

Strategy C needs more intermediates alive simultaneously (dx², v, w_reg,
Rx_tmp, plus Rx_reg, core_reg, numer_reg at various times). Each
additional live n-bit register adds 256 to peak.

If we keep dx² + v + w_reg + Rx_tmp + r_low + m_hist + tx + ty alive
at body peak, that's 8n + 407 = 2455. Add step-4 transient 260 (for
Kaliski backward) = 2715. Right at current peak.

**So Strategy C has SIMILAR peak to current.** The question is whether
we save Toffoli.

### Toffoli budget for Strategy C (Phase 3 body)

Inside body:
- `dy2 := ty·ty`                 130k (sq)
- `v := dy² - dx²·(Px+Qx)`       150k (mul-sub with classical constant) + maybe 2k adder
- `Rx_tmp := tx·r_low`            150k (mul)
- `tx (:= Rx) := v · Rx_tmp`      150k (mul) + reverse path (uncompute old tx?)

For the Ry side:
- `core := dx²·Qx - v`            ~80k (q·classical)
- `dy·core`                       150k
- `w·Qy`                          ~80k
- `numer := dy·core - w·Qy`       1.5k (sub)
- `ty (:= Ry) := numer · r_low`   150k (mul)

Subtotal forward Phase 3: ~1.04M

Uncompute (Bennett-clean) for each intermediate:
- dy2 uncomputed: 130k
- v back to 0 via recomputing dx²·(Px+Qx) and adding: 150k 
- Rx_tmp uncomputed: 150k
- dx²·Qx recompute + core uncompute: 80k
- dy·core uncompute: 150k
- w·Qy uncompute: 80k
- numer uncompute: 150k

Subtotal Phase 3 uncomputes: ~890k

Phase 1 (dx² and dx³ compute): 130k + 150k = 280k

Phase 2 (Kaliski on w): ~800k fwd + ~800k bwd = 1.6M (one Kaliski pass).

Phase 4 (uncompute w_reg = dx³ and dx2_reg after Kaliski backward):
Kaliski backward already uncomputed to leave only r_low changed? No,
Kaliski uncomputes all its state but tx (its input) is preserved. So
w_reg still holds dx³ after backward. We need to uncompute w_reg =
dx³ → 0 by reversing the dx³ mul. Cost: ~150k.
Similarly uncompute dx2_reg (from dx² → 0): 130k.

Plus final classical ops to turn tx back from something to Rx, ty to Ry.

**Total Strategy C estimate**:
- Phase 1: 280k
- Phase 2 (Kaliski): 1.6M
- Phase 3 fwd: 1.04M
- Phase 3 unc: 890k
- Phase 4 (uncompute dx², dx³): 280k
- misc: 100k
- **Total: ~4.2M**

Current baseline: **4.14M**.

**Strategy C estimate is NOT an improvement.** It's a wash, maybe
slightly worse.

### Why the single_inv_plan.md said 3.7-4.3M

The plan used an older baseline (4.91M @ iters=511) and estimated 3.7M
for Strategy C. At iters=407 the old baseline is now 4.14M, and the
same 3.7M became less of an improvement. My re-estimate above (4.2M)
agrees with the upper end of the old range.

**Conclusion**: Strategy C is not a clear win at iters=407. It would
be a CLEAR win if iters were 511 (baseline 4.91M → Strategy C ~4.2M
would be -0.7M = -14%).

But at iters=407 (where we are), **Strategy C is NOT a Toffoli lever**.

## What went wrong with the excitement

I was excited about Strategy C because:
1. It's classically verified (200/200).
2. The single_inv_plan projected -0.6M to -1.2M savings.

But the projection used iters=511 as the baseline. At iters=407, our
CURRENT Kaliski is much cheaper than "half of 2 × Kaliski @ 511" — the
iters reduction was a per-pass optimization that already realized part
of the Strategy C savings.

## What remains as candidate paths

- **Path A (coset on long add chains)**: INVALIDATED — our circuit
  doesn't have 12+ consecutive adds in any one place. Crossover not
  reached. Skip.
- **Path B (m_hist elim + karatsuba pair1_mul2)**: PARTIALLY INVALIDATED
  — the naive m_hist elimination fails because iter-END fingerprint
  doesn't determine m_i (only iter-START does, but iter-start state is
  gone at iter-end). Sophisticated protocols might work but are novel
  research.
- **Path C (step 3+9 cswap restructuring)**: open; novel research.
- **Path D (single invocation via approximation)**: open; requires
  harness change + novel cleanup protocol.
- **Strategy C (this doc)**: NOT A WIN at iters=407.

## What to actually do

Stop chasing structural levers that have already been invalidated.
The tactical path forward is **per-iter Kaliski cost reduction** —
each 1% saved is ~40k Toffoli. Small wins, but they don't require
novel research.

Specific ideas worth thinking about (before picking one):
- Narrow the step-4 tmp load from `min(n, 2n-i)` to something tighter
  based on live-bit-extents during bulk iterations (the Kaliski
  invariant says bitlen(u) ≤ 2n-i, so at iter_idx=100 the load could
  be 156 bits, not 256).
- Replace `mod_double_inplace_fast` inside step 7+8 for big-iter range
  with `mod_double_no_corr` when a tighter invariant can prove r[255]=0.
  (Already done for iter_idx < R_SMALL_THRESHOLD = 255 — check if the
  threshold can be raised.)
- Check if step 2's `with_gt` comparator can use a narrower bit-width
  than current for bulk iters. (May have been done via truncation
  already.)

None of these are huge Toffoli wins individually, but they're the
remaining realistic levers under exact correctness + 2800q cap.
