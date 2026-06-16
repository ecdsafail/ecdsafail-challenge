# 2026-06-16 — Sub-800 structural panel: corrected 999 model (4-agent measured)

A 4-way parallel worktree panel (levers B1/B2/B3/C) attacked the 999→<800 qubit
peak. **No lever reached sub-800.** Best measured: **986** (lever C, Track-A only).
But the four agents independently produced a *corrected, measured* model that
**overturns `SUB_800_PLAN.md`'s core assumption.** Trust the numbers below over the
original plan; everything here was measured in count-only mode against
`peak_phase`/section traces.

## What the original plan got wrong

`SUB_800_PLAN.md` modeled 999 as ONE step: `pack 741 + dy passenger 256 + flags 2`,
located in the **division**, and claimed **Track B alone reaches ~745 (sub-800)**.
Both claims are false:

1. **The dominant peak is the MULTIPLY substep's clz, not the division.** Confirmed
   by tagging the division's `clz_diff_body_middle` ("DIVTAG") and seeing the 999
   `peak_phase='ec3.inv_fwd/p.bitlen'` fall *outside* the tag. In the multiply
   substep, `ca` is the LIVE accumulator (`ca += cb<<s` right before the peak clz),
   so it is data-bearing — not free — at the exact peak.
2. **There are multiple structurally-forced 999 regions**, each carrying its own
   256-bit passenger that is sandwiched between ops needing it:
   - `ec3.inv_fwd` forward inversion — passenger = **dy**
   - `ec3.inv_fwd` backward inversion — passenger = **lambda** (the divide output;
     must stay live to uncompute cb)
   - `ec3.alt.cancel` — passenger = **new_dy**
3. **Track B alone CANNOT reach sub-800.** B3 measured the strictly-stronger bound:
   freeing **all four** 256-bit passengers leaves the floor at **808**, now bound by
   the slope reconstruction `mod_mul_rfold_mbu` (`shrunken_pz_state_machine.rs` ~1199:
   `lambda + dx + dy_new + temps ≈ 808`). Passenger elimination is *necessary but
   not sufficient*.

## Measured decomposition of 999

```
999  ≈  pack(~702)  +  one 256-bit passenger  +  clz pa/pb transient(~41)
```
- Peak step ≈ 348–361; `reg_widths` row `[A,B,ca,cb,q] ≈ [~86,~86,254,254,~22]`.
- **ca and cb are BOTH pinned at ~254** (cofactor×remainder≈p forces `ca+cb≈512`).
  No asymmetric slack at the peak — both cofactors are full-width.
- A,B are only ~86–98 bits at the peak ⇒ **no 256-bit zero-hole exists** to relocate
  or pack a passenger into. Peak is information-theoretic; freeing slots elsewhere
  does not lower it.

## Why each lever died (so nobody re-tries them blindly)

- **B1 freed-space reuse** — no 256-bit zero region at the peak (A,B≈86; ca,cb hold
  worst-case data). Backward inversion enters at active 789 carrying lambda
  unremovably. dy and lambda are mutually defining through the consumed dx, so one
  256-bit quantity must always cross the 702 pack. Floor ~960.
- **B2 late dy reconstruction** — **impossible**: `dy = oy − ty_orig` is full-entropy
  input with no surviving source; `ty_orig`/P.y is overwritten in place at
  `ec/point_add.rs:199` and never copied. At the multiply, dx is already torn down
  (B=1) and lambda doesn't exist yet, so `dy = lambda·dx` can't be formed either.
  Also fixes only 1 of 3 regions. (Did measure: forward inv *without* dy peaks at 742.)
- **B3 reversible packing** — nothing to pack into (no zero region); and even total
  passenger elimination floors at 808 (see above). Wrong lever for this circuit.
- **C one-cofactor pack** — `ca` is the live multiply accumulator at the peak, so
  "drop ca, recompute from Bézout" requires a **nested reversible division inside
  every multiply substep** — a multi-day rewrite, extreme correctness risk (one
  mismatch fails all 9024 shots). Track-A schedule tightening
  (`VALIDATE=0, CLZ_WINDOW=64, TRAIN=20000`) narrows only A/B/q → **986 floor**,
  never touches ca/cb.

## The real path to sub-800 (for the next attempt)

Sub-800 needs a **structural rewrite combining TWO reductions**, not a single lever:

1. **One-cofactor EEA** (proper Lever C): run the inversion keeping only `cb` live and
   reconstruct `ca` via the Bézout relation `ca·a + cb·b = gcd` *at the multiply
   where ca is the accumulator* — i.e. a nested-EEA-per-step restructure, removing
   ~254 from the pack across all shared inversion regions.
2. **In-place / windowed slope multiply**: replace `mod_mul_rfold_mbu` (~1199) so it
   does not materialize a fresh 257-bit `dy_new` alongside live `lambda + dx` — this
   is the 808 floor that bites *after* the cofactor is dropped.

Either alone lands ~808–986; **both together** are required to clear 800. This is a
multi-day reversible-arithmetic effort, not a knob/allocation tweak. Budget for it.

## Track A is DEAD for the qubit peak (measured frontier, 2026-06-16)

Follow-up to the panel: tried to "bank" the agent's count-only 986 by passing the
9024-shot correctness gate. It cannot be done. Measured peak-vs-miss frontier (misses =
`miss_factors` from `TRAILMIX_SUPPORT_CHECK`, the harness's own correctness oracle —
`miss_factors>0` ⟺ some real draw's factor exceeds the narrowed width ⟺ that shot
computes WRONG ⟺ eval FAILS):

| schedule (CLZ=128, TRAIN=65536) | peak | miss_factors |
|---|---|---|
| VALIDATE=500000 (baseline)      | 999  | **0** |
| VALIDATE=100000–300000          | 999  | 5–13 |
| VALIDATE=50000                  | 996  | 22 |
| VALIDATE≤9024                   | 986–989 | 34–91 |

- `CLZ_WINDOW=64` and `TRAIN=20000` add misses but give **zero** peak reduction (peak
  stays 999) — pure liabilities; don't touch them.
- The ONLY peak lever is lowering `VALIDATE`, which narrows pack widths — but the FIRST
  qubit of reduction (999→996) already costs 22 misses. A clean tail nonce needs 0
  misses; at ≥22 expected independent heavy-tail misses no feasible nonce search finds
  one (P≈e^-22). Confirmed empirically: a 1,000,000-nonce search at the 996 config found
  **no clean nonce**.
- Root cause: covering the worst of 9024 near-uniform mod-p factors forces ~full 254-bit
  widths (extreme-value statistics) → pack ~702 → peak 999. Reducing the peak *requires*
  under-covering some draws, which breaks them. **There is no valid sub-999 schedule.**

Conclusion: the realistic VERIFIED floor with schedule/knob tuning alone is **999**.
Sub-800 is ONLY reachable via the two-part structural rewrite above (one-cofactor EEA +
in-place slope multiply). Do not chase schedule tuning for the qubit peak again.

## Fast tooling reminders
- Peak (count-only, seconds, 16-byte ops.bin):
  `TRAILMIX_THIN_TRACE=1 TRACE_PEAK=1 POINT_ADD_COUNT_ONLY=1 target/release/build_circuit | grep peak_qubits`
- Section-active maxima via named trace sections + `TRACE_PHASE_ACTIVE`/`TRACE_EACH_PEAK`
  is how the regions above were separated — use it to attribute any new peak.
- Authoritative gate (only on a count-only sub-800 candidate): full
  `build_circuit && eval_circuit` → needs `qubits<800` AND `9024/9024 OK`; then
  `rm -f ops.bin` (7.6 GB).
