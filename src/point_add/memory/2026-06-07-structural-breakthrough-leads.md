> **CORRECTION (see `2026-06-07-measured-frontier-leads.md`):** build_circuit traces
> show the peak binder is the **GCD walk** (`compressed_block_tobitvector_reverse_add`/
> `_shift` @ 1302), **not round84** (which is at 1284, 18 q of slack). Lead B below
> (round84 square) therefore does **not** move score at the current peak. The cswap
> figure in Lead A was also low: measured cswap = 271,744 CCX (18.7 %), but the apply
> half (132 k) is full-256 and not truncatable. Read the measured note first.

# Structural breakthrough leads (analysis-only, 2026-06-07)

Model: Claude Opus 4.8. Method: static read of the whole `src/point_add/` tree.
**No benchmark was run for this note** — every number below is either read from
the source/tests or an order-of-magnitude estimate that you MUST confirm with
`TRACE_PHASES` / `TRACE_PHASE_ACTIVE` before spending implementation time. Treat
this as a map of where the big score is hiding, not as a validated result.

Written because the recent loop (see `2026-06-06-tony-anton-audit-loop.md`) has
collapsed into bit-truncation + Fiat-Shamir nonce hunting: every win for the last
many submissions is 0.03 %–0.3 % of score and needs a fresh `DIALOG_TAIL_NONCE`.
That well is dry. The leads here are structural (new divstep / new uncompute /
new scratch layout), which is what actually moves a mature circuit.

---

## 0. Frontier re-anchor (the inherited note is stale)

`configure_ecdsafail_submission_route()` in `mod.rs` is currently wired to the
**tier-3 "safe lock"** route, not the route the 2026-06-06 memory describes:

- `DIALOG_GCD_BODY_CARRY_BAND_TRIMS = "0,3,3,3,...,3,3,3"`, `FUSED_OVFCLEAR_MEASURED=1`,
  `APPLY_CHUNKED_F_CUT4=189`, `ROUND84_INPLACE_SOLINAS_FOLD=1`,
  `ROUND84_INPLACE_QUOTIENT_CARRY_TRUNC_W=21`, `DIALOG_TAIL_NONCE=11201395269`.
- The in-code comment (mod.rs ~line 1356) claims this validates **1302 q ×
  1,456,963 T = 1,896,965,826**.

So the live baseline in the tree is **~1.897e9, peak 1302**, *better* than the
2026-06-06 note's 1,960,613,655 / 1309 q. Re-validate `./benchmark.sh` once to
confirm which one your checkout actually reproduces before comparing against it.
(Score = avg executed Toffoli × **peak** qubits; lower is better.)

---

## 1. Where the cost actually goes (the cost map)

`emit_dialog_gcd_raw_pa` (rounds/dialog/mod.rs:1820) is the whole PA. It is
**exactly two GCD modular inversions** wrapped around one square:

1. `pair1_quotient` — GCD-invert `dx = x1−Qx`, divide `dy` by it → `ty = λ`,
   `tx` kept `= dx`.
2. `round84_emit_fused_square_xtail` — `tx ← λ² − dx − 2·Qx = Rx`. **This is the
   peak binder.** It is `tx (256) + λ (256) + a 2N = 512-qubit product register
   `tmp_ext` + per-row carry scratch`. The in-place Solinas fold
   (`ROUND84_INPLACE_SOLINAS_FOLD`) folds hi→lo *after* the full square to claw
   peak down to ~1302–1307.
3. `c = Qx − Rx` into `tx`.
4. `pair2_product` — GCD-invert `c`, use it to uncompute `λ` → `ty` becomes `Ry`.
5. `ty −= Qy`; restore `tx → Rx`.

Each GCD (`emit_dialog_gcd_raw_{quotient,ipmul}`) is:
`tobitvector_steps` (forward divsteps, build the transcript `dialog_log`) →
`apply_bitvector` (replay transcript onto the full-width target = the actual
multiply/divide) → `tobitvector_steps_reverse` (Bennett-uncompute u and the log).

`tobitvector` step body (rounds/dialog/mod.rs:670-711), per active step:
- branch bits: `cx(v0,b0)` + a truncated comparator (`compare_bits`, scheduled
  down to avg ~40, min 5) → 2 log bits `b0`, `b0_and_b1`.
- **`cswap` of `u_active,v_active`** — `cswap` (adder.rs:951) is `cx;ccx;cx` =
  **exactly 1 Toffoli per bit**, run over the full active width.
- **`controlled_sub` of `v` from `u`** — another full-active-width Cuccaro pass.
- `shift_right_assuming_even(v)` — free (relabel).

So **every divstep does TWO full-active-width Toffoli passes (cswap + sub)**, and
this body runs **4 times** total (forward+reverse × 2 GCDs), once more in each
`apply`. Active width runs ~256 → ~4 over 258 steps (slope ≈ 1.015/step).

### Toffoli budget (estimate — verify with `TRACE_PHASES`)
- Σ active_width over 258 steps ≈ 34 k per pass.
- `cswap` alone ≈ 34 k × 4 passes ≈ **~135 k Toffoli ≈ 9 % of the 1.46 M total**,
  and it uses **no extra scratch** (in-place Fredkin) → removing it is
  **peak-neutral, pure Toffoli**.
- The two `apply` passes (full 256-wide modular add/sub per step) are the other
  large block; already heavily worked (`MEASURED_APPLY_SUB`, chunked-F, fused-fold).

---

## 2. What is already exhausted — do NOT re-spend cycles here

- Fiat-Shamir nonce search (`DIALOG_TAIL_NONCE`, `DIALOG_REROLL`,
  `DIALOG_POST_SUB_REROLL`). The whole 2026-06-06 loop is island-limited; more
  blind sweeps will not find a *structural* win.
- One-bit truncation knobs: `COMPARE_BITS`, `APPLY_CLEAN_COMPARE_BITS`,
  `WIDTH_MARGIN`, `WIDTH_SLOPE`, `KAL_DOUBLE/FOLD_CARRY_TRUNC_W`, the per-step
  compare schedule + margin. Each is ≤0.3 % and re-rolls the island.
- `ACTIVE_ITERATIONS` micro-tuning — sits on a nonconvergence floor.
- **One-inversion PA is provably blocked** for this clean in-place ABI
  (`ONE_INV_DX3_AFFINE_PA_BLOCKER`, mod.rs:507; test at mod.rs:1609). Recovering
  `dx` after `(tx,ty)` are overwritten by `(Rx,Ry)` needs inverting `Rx−Qx` =
  a second inversion. **Two inversions is a hard floor. Stop anyone chasing 1.**

---

## 3. Lead A — kill the per-step `cswap` (highest leverage, peak-neutral)

**Claim:** the divstep spends a separate full-active-width Toffoli pass on
`cswap(u,v)` *in addition to* the controlled-subtract. Literature reversible
binary-GCD / safegcd divsteps fold the swap into the arithmetic (one
conditional ±-subtract steered by a sign/`delta` counter), so the swap pass
disappears. Estimated **~8–10 % Toffoli, 0 qubit cost → ~8–10 % score**, with no
new Fiat-Shamir hazard class (it changes *how* you subtract, not *which bits* you
truncate).

**Strong tell:** `configure_ecdsafail_submission_route` already cites
**Gidney et al., arXiv:2510.10967** — but only for its *width bound*
("after i iters 2·deg(b) ≤ 2d−1−i−δ", mod.rs:1274). That paper is a reversible
safegcd/inversion construction; its **divstep circuit** almost certainly fuses
the swap and is the thing to port, not just its inequality. Read the paper's
divstep, not its appendix bound.

**Concretely:**
1. First *measure* the prize: run with `TRACE_PHASES=1` and read the Toffoli
   attributed to phases `dialog_gcd_raw_tobitvector_cswap` and
   `..._reverse_cswap` across both GCDs. If it is ≥100 k, this lead is real.
2. Replace the `cswap` + `controlled_sub_selected` pair with a single
   sign-steered conditional add/subtract (Bernstein–Yang `divstep`: track a small
   `delta` counter ~9 bits; the branch that currently swaps becomes
   `g ← (g − f)/2` with `f,g` roles selected by `sign(delta)` instead of a data
   swap). The transcript stays ~2 bits/step (log the `delta>0 ∧ g_odd` branch).
3. Keep the existing width-envelope / active-width truncation — it is orthogonal
   and carries over to the BY recurrence (the cited bound is already a BY bound).

**Risk:** this is a genuine re-implementation of the GCD core (forward, reverse,
and the matching `apply` that consumes the new transcript). High effort, but it is
the single largest peak-neutral Toffoli block in the circuit and the one place the
codebase has a paper it half-used. Verify the swap-free divstep is actually
swap-free in the Toffoli model first (some BY formulations still hide a conditional
swap — confirm the paper's does not before committing).

---

## 4. Lead B — shrink round84's 512-qubit square transient (only true peak lever)

**Strategic fact agents keep missing:** the global peak is a *co-bind* between
round84 (the square) and the GCD-walk, both pushed to ~1302. Therefore:

> **Cutting GCD-side qubits below the peak does nothing to score.** The 2026-06-06
> loop's `1285q` restacks were chasing below-peak slack. Score only moves if you
> cut **both** the round84 transient **and** the GCD-walk peak.

round84's square (`schoolbook_square_symmetric*`, multiply.rs:320+) materializes a
**`tmp_ext` of 2N = 512 qubits** for `λ²` before reducing. That 512-wide block is
the largest single scratch in the circuit and it sits exactly at peak. The
in-place Solinas fold already reclaims part of it *after* the fact.

**Idea:** interleave Solinas reduction *into* the accumulation so the product
never fully materializes to 512 bits — stream each high cross-product back through
`2^256 ≡ 2^32 + 977` as it is produced, keeping the accumulator at ~256 + a small
carry band instead of 512. Target: drop the square transient by ~100+ qubits.
Pair it with whatever simultaneously trims the GCD-walk co-peak (e.g. the
transcript-block borrow levers already in the tree) so the *global* peak actually
moves. Each qubit off the global peak is ~1.46 M / 1302 ≈ **1,120 score per qubit**
— i.e. one qubit ≈ two of the recent nonce-grind submissions.

**Risk:** the symmetric square writes cross-products to shifted positions
2i..i+n; streaming reduction must fold high words while later rows still write
into them. Medium-high. But this is the only axis that beats the score *without*
touching the Fiat-Shamir island at all.

---

## 5. Lead C — cheaper transcript-log uncompute (smaller, "uncompute idea")

`tobitvector_steps_reverse` restores `u→p` and `v→factor` (genuinely needed) but
its **only** redundant work is recomputing the truncated comparator each step to
clear the 2-bit log (`b0`, `b0_and_b1`); the cswap/sub there are driven by the
already-present log bits. `b0` is cleared by one `cx(v0,b0)` (free). `b0_and_b1`
still pays a comparator recompute.

Idea: clear `b0_and_b1` by **measurement-based uncompute** (Hmr + phase feedback)
the way the apply-phase AND-clears already do (`FUSED_*CLEAR_MEASURED`,
`MEASURED_APPLY_SUB`). The blocker is that the Gidney phase correction needs the
two set-time controls (`b0` and `cmp = u>v`) live at measure time; `b0` is live
(`= v0`) but `cmp` is a freed ancilla. If `u>v` can be re-expressed cheaply from
currently-live bits (it often resolves on a handful of top bits, exactly what the
per-step compare schedule already exploits), the comparator recompute collapses to
a phase-only correction. Lower leverage than A/B (comparator is already scheduled
small) but it is real, low-risk, and island-neutral. Good warm-up before Lead A.

---

## 6. How to verify any of this BEFORE writing the circuit

- `TRACE_PHASES=1 cargo run --release` → per-phase emitted Toffoli. Confirms the
  cswap / square / apply split and sizes Lead A and B.
- `TRACE_PHASE_ACTIVE=1` (+ `TRACE_PHASE_ACTIVE_REGIONS=1`) → per-phase live-qubit
  maxima. Confirms the round84 ↔ GCD-walk co-bind and which phase is the true peak.
- `DIALOG_GCD_RAW_PA_STOP_AFTER_{QUOTIENT,XTAIL,C,PAIR2}=1` → bisect the PA to
  attribute Toffoli/peak to each stage in isolation.
- Component tests already exist: `round84_fused_square_xtail_component_matches_relation`
  and the `d1_inplace_*_lowerer` pins (mod.rs:838+) — use them as fast oracles for
  a new divstep/square without the full 9024-shot run.

Ranking by expected score impact: **A (~8–10 %) > B (~7 %, harder) > C (small,
safe)**. A and B are independent and stack. None of them touch the nonce search —
that is the point.
