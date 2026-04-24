# Algebraic & scheduling tricks from the literature we have NOT used

Mined from `/tmp/*.pdf` + arxiv queries, 2026-04-24 session. Ordered by
expected impact for our "n=256 secp256k1 point-add, ~600 qubits over the
512-qubit input register" target.

## 1. Luo 2025 — register sharing via the EEA invariant  (BIG qubit win)

arXiv (retrieved to `/tmp/luo_ec.pdf`, `.txt`). Reports full Shor's
ECDLP on secp256k1 in **1333 total qubits** (vs HRSL 2020's 2124, vs
RNSL 2017's 2338, vs our 2716 for one point-add).

Central trick: the Extended Euclidean Algorithm produces two sequences
`{r_i}` (monotonically decreasing) and `{t_i}` (monotonically
increasing) with the bilinear invariant

    r_{i-1} · t_i + r_i · t_{i-1} = p       (over the whole run)

This means `bitlen(r_{i-1}) + bitlen(t_i) ≤ n+1`, so the two values
**share a single n+2 qubit register** with a known split. Luo pushes
this further: `(r_{i-1}, t_i, q_i)` can all share one n+2 register.

I added `src/point_add/luo_proto.rs` to budget this numerically against our
current scaffold. Results at n=256:
- conservative swap-in estimate: **2084q** total,
- optimistic overlap estimate: **1828q** total,
- savings vs current 2716q: **632–888 qubits**.

So Luo is definitely a *real* qubit lever.
However, even the optimistic estimate is still far above the user's rough
budget of `512 + 600 = 1112` total qubits.

Cost: Luo uses a long-division-style EEA (Proos-Zalka shape), not
Kaliski's halving-based one, so the per-round Toffoli profile is
different. Their full Shor-algorithm Toffoli count is ~976 n³ =
~1.6·10¹⁰ for n=256 whole-algorithm; extracting per-point-add requires
dividing by the 2n=512 point-adds, giving ~3.1·10⁷ Toffoli per point-add
— worse than our 4.18M on raw count. So this is a qubit-vs-Toffoli
tradeoff: Luo wins qubits by 50%, loses Toffoli by ~7x at the same n.

**Takeaway:** register sharing is the single biggest-ticket qubit lever
we're not using. BUT implementing it means switching from our
Kaliski/swap formulation to a PZ/long-division formulation — a total
rewrite of `with_kal_inv_raw` and friends. Multi-week effort.

## 2. Kim 2026 — unconditional execution + postponed modular reduction

**Update from this session:** I added `src/point_add/kim_proto.rs` and split the question in two:

1. **Naive narrow-r port** (wrong model): FAILS. If we keep `r` truncated to 256 bits, then
   - unconditional tail ≠ fixed doubling of the conditional result,
   - final `r` after 512 rounds ≠ `±x^{-1}·2^{2n}`,
   - so the naive claim is false under our current narrow-r classical model.

2. **Kim-style wide-r port** (the right model): PASSES. If we keep `r,s` as 512-bit values and only reduce mod p at the end, then over 200 random secp256k1 inputs:
   - conditional Kaliski termination step `k` ranges in **[339, 384]**, avg **363.06**,
   - unconditional execution to 512 rounds equals the conditional result followed by a deterministic `2^(512-k)` doubling tail,
   - and the low 256-bit residue of the final wide `r` is indeed `±x^{-1}·2^{2n}`.

Interpretation:
- Kim 2a **survives** as a real formulation path, but only with a widened `r,s` model and postponed reduction.
- Therefore it is **not** a low-risk tweak to our present Kaliski scaffold.
- Real import path is probably "Kim 2a + Kim 2b together" (unconditional execution plus postponed reduction), not 2a alone.

Paper `/tmp/kim_2026.pdf`. Two linked tricks, Section 3.3.1:

### 2a. Unconditional execution
Run **all 2n Kaliski rounds unconditionally** — no `v==0` check in any
round. In the rounds after termination, r keeps being doubled into a
2n-bit r register (the arithmetic "overshoots" cleanly). The result is
**already in Montgomery form `x⁻¹ · 2ⁿ mod p` on the nose, every shot.**

Consequence: **no `pair1_halve` or `pair2_double` correction loops**.
Our current build pays
- pair1_halve: 103,785 CCX
- pair2_double: 103,020 CCX

**Direct saving: ~207k CCX (~5.0% of total)** just from eliminating
these two phases. Wireable into our existing scaffold with minimal
structural changes; the `with_kal_inv_raw` body would just skip the
halve/double loops.

### 2b. Postponed modular reduction
Inside each Kaliski round, the v←v-u, r←2r, etc. ops can defer their
mod-p reduction to the final stage. Accumulate into a 2n-bit register
across all rounds; do one big reduction at the end using a QCSA
(quantum carry-save adder) or sequential controlled additions.

Our current `kaliski_iteration` already does modular reduction inside
each round (for r via `mod_double_inplace_fast`/`mod_double_no_corr`
and for s similarly). Eliminating that per-round reduction saves
roughly `O(n) CCX per round × 2n rounds = O(n²)` CCX — napkin estimate
for n=256 puts this in the same ballpark as the inversion cost itself.

**Takeaway:** 2a alone is low-risk, ~5% win, probably 1-day
implementation. 2b is deeper and higher-value but wants a register
size change (2n-bit r) that we currently don't accommodate.

## 3. HRSL 2020 — register reuse across inversion and multiplication

Paper `/tmp/hrsl_2020.pdf`, Fig. 8b. Modular **division** `y/x` is
composed as: invert x (via Kaliski into aux registers), copy-out, undo
the Bennett wrapper. At that point three of the aux registers hold the
known values `{0, 1, p}`. HRSL: **clear them with n X-gates and reuse
them as the workspace for the multiply `y · x⁻¹`**.

Net peak qubit count for divide: ~8n instead of 8n+3n. Saves 3n = 768
qubits during a divide, exactly the workspace size we need to stay
under our cap.

Our current `with_kal_inv_raw` frees `u`, `v_w`, `f_flag` after the
forward pass (already a reuse). We do NOT currently fold the mul
workspace into that reuse: pair1_mul2 allocates its own tmp_ext (2n)
and carries (n). Wiring the mul to reuse the freed Kaliski registers
would save the mul's tmp_ext allocation, potentially ~512 qubits of
transient peak.

**Takeaway:** medium-risk scheduling rewrite. Requires threading the
mul's allocator through `with_kal_inv_raw` so the freed u/v_w slots
get reused instead of new alloc. Could turn a 258q transient into a
0q transient at the pair1_mul2 peak, finally making karatsuba-1
fit in budget at pair1_mul2.

## 4. HRSL 2020 — swap-based Kaliski round

Same paper, Fig. 7b / Algorithm 7b. Instead of our current 4-branch
structure (one of 4 sub-operations per round), use:
- 2 controlled swaps (u↔v, r↔s) at the beginning
- 1 subtraction `v ← v - u` conditioned on `u odd AND v odd`
- 1 addition `s ← r + s` same condition
- 1 unconditional `v ← v/2`
- 1 unconditional `r ← 2r`
- 2 controlled swaps at the end to undo

One subtract + one add per round, vs our four. We partially implement
this via `kaliski_iteration_bulk_prefix3` for the first 315 rounds;
the generic `kaliski_iteration` still carries the 4-branch machinery.

**Takeaway:** fold the swap-based formulation into the full 407/404
iters (not just the bulk prefix). Estimated Toffoli saving: ~500k CCX
across both Kaliskis (roughly 2 of our 4 step4-style ops per round
disappear). Needs very careful phase-correctness work — the partial
bulk-prefix3 is the template.

## 5. Gidney 2019 — windowed classical-quantum addition

Paper `/tmp/gidney_windowed_2019.pdf`. Replace k back-to-back
`quantum_reg += classical_const_i` operations with a single QROM
lookup of the sum-of-k-constants followed by one add. Cost:
`O(W·L/k + k)` Toffoli where `W = n` and `L = 2^k`.

Applies to our `mod_add_qb`, `mod_sub_qb`, `mod_add_double_qb` calls.
In current build we have ~6 such calls × ~1280 CCX = ~7680 CCX total.
Windowing would cut this by a factor of ~4-8 on the per-op cost,
saving ~6000 CCX. Small absolute win.

**Takeaway:** low-impact, skip for now.

## 6. Boneh's Montgomery batch inversion trick

Actually a Montgomery/Peter-Montgomery technique but widely
attributed to Boneh (Stanford crypto course material, BLS signature
paper context). Bundle two inversions into one:

    a⁻¹ = b · (ab)⁻¹
    b⁻¹ = a · (ab)⁻¹

Already explored in Strategy A of `single_inv_plan.md`. Result:
**DEAD** because the Py term gets trapped in the ty register as `ty=Py+Ry`
and we have no classical handle on Py to subtract it out.

**Takeaway:** Boneh's trick is not a new avenue for us — it was
Strategy A, confirmed dead by the 200-trial falsification harness.

## 7. Chevignard 2026 — RNS + projective coords

Paper `/tmp/chevignard_2026.pdf`. Uses Residue Number System
representation + projective coordinates, avoiding the modular
inversion entirely inside the point-add loop by maintaining `(X:Y:Z)`
and only dividing by Z at the end.

Not applicable to our session's contract: our `src/main.rs` harness
requires `(Rx, Ry)` to be written into the SAME quantum registers as
the input `(Px, Py)`, in affine form. A projective representation
would require a full RNS-to-affine conversion which itself is an
inversion — no saving.

**Takeaway:** not applicable, ruled out previously.

---

## Prioritised moves for next session

Ranked by impact × ease, updated after the falsification work in
`src/point_add/kim_proto.rs` and `src/point_add/luo_proto.rs`:

| # | trick | impact (CCX save) | qubit ∆ | implementation risk |
|---|-------|------------------:|--------:|----------------------|
| 1 | HRSL swap-based full Kaliski | -400k to -600k | 0 | MEDIUM (phase-correctness) |
| 2 | HRSL register reuse for pair1_mul2 / mul workspace | 0 direct, but enables -28k/site | -258 peak | MEDIUM (allocator threading) |
| 3 | **Kim 2a+2b together** (unconditional exec + postponed reduction) | maybe -500k to -700k | +n in inversion core, -207k loops removed | HIGH |
| 4 | Luo register sharing | 0 or negative Toffoli | **-632 to -768 peak** | VERY HIGH (rewrite Kaliski as PZ long-division) |
| 5 | Gidney windowed classical-quantum add | ~-6k | 0 | LOW |

Revised combined hope with 1+2:
- **~430k to 630k CCX saved**,
- **~250q peak saved**,
- likely best stable target around **3.55M to 3.75M Toffoli @ ~2460q**.

If Kim 2a+2b lands on top of that, the floor could move toward
**3.1M–3.3M**, but with a much heavier inversion rewrite.

If Luo register sharing also lands, qubits plausibly drop to
**~1.83k–2.08k** in our affine scaffold (measured in `luo_proto.rs`),
still not enough to meet the user's aspirational ~1100–1200 total
without a broader scaffold collapse.

### Hard conclusion from the budget work
To hit the user's qubit target, **inversion tricks alone are insufficient**.
Even a hypothetical *free* inversion block would leave the current affine
scaffold at:
- `tx,ty` live: 512q
- non-inversion affine scaffold: 772q
- total: **1284q**

So after Luo-style inversion sharing, we still also need to collapse at least
~172 qubits from the surrounding affine scaffold — practically, this means:
- kill or fold `lam`, and/or
- eliminate one wide mul transient, and/or
- change the add choreography so one whole n-wide register never exists.
