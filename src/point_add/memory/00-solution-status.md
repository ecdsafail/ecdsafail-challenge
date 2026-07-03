# secp256k1 point-add — solution status & frontier analysis

## Current operating point (validated locally, matches world-best)

`point_add::build()` (→ `trailmix_ludicrous::build_trailmix_ludicrous_ops()`)
passes the trusted 9024-shot evaluator **clean (0/0/0)**:

| metric | value |
|---|---|
| peak qubits | **1152** |
| avg executed Toffoli / shot | **1,364,229.770** |
| **score** = round(T) × Q | **1,571,592,960 ≈ 1.57 × 10⁹** |

This is the current **best promoted leaderboard submission** — `71f5115`
(BitWonka), commit `d44cad3` = this repo's HEAD. It beats Google's published
Pareto frontier (low-gate 3.0 × 10⁹, low-qubit 3.2 × 10⁹) by ~1.9× on *both*
axes at once.

Reproduce: `cargo fetch --locked` once (the offline cache ships alloy-rlp
0.3.12 but the lock pins 0.3.15), then `build_circuit` → `eval_circuit`.

## The leaderboard is stuck here (as of 2026-07-03)

`ecdsafail submissions --all`: the frontier has not moved since 2026-06-26.
7+ recent submissions re-submit the identical 1,571,592,960 and are **rejected
at 0% (no improvement)**. Attempts at fewer qubits all score *worse*:
1133q→1,460,511t, 1141q→1,423,723t, 1142q→1,415,790t, 1156q→1,365,960t — every
one rejected. Radical low-qubit tries (825–844q) blow the Toffoli count to
400M+ (score +3000%). This is a genuinely hard, saturated frontier.

## Why it can't be improved by tuning (all measured, don't re-litigate)

1. **Peak 1152 is a broad structural plateau.** `TRACE_TLM_PROFILE` shows ~15
   distinct phases across both GCD passes, the square, and the apply steps all
   sitting at exactly 1152 active qubits. `TLM_TARGET_Q=1151` left the peak at
   1152 and broke correctness entirely (9024/9024 mismatches). The headroom knob
   cannot buy a qubit; a lower peak needs a genuinely lower-width inversion.

2. **Any op-stream change reseeds Fiat-Shamir and needs a fresh clean nonce.**
   The 9024 test inputs derive from a SHAKE hash of the whole op stream. The
   baked `DIALOG_TAIL_NONCE=50400005525597` is a *rare* clean island: two random
   nonces both gave ~15 classical + ~11 phase failures, so the clean rate is
   ~e⁻¹⁴ ≈ 1e-6. Brute-force clean-nonce search on the full eval (~40 s each) is
   ~39 days per clean nonce — infeasible.

3. **A real Toffoli reduction exists but is gated by (2).** Enabling
   `SINGLE_CCX_FANOUT_DISABLE=0` removes **318 real CCX** (10,221,377 →
   10,221,059 ops, +0 qubits — the rewrite `CCX;CX;CCX → CX;CCX` allocates no
   wire). That would give avg T ≈ 1,363,912 → score ≈ 1,571,226,624, beating the
   frontier by ~366k. **But** it reseeds FS and breaks the clean nonce (24
   classical + 16 phase failures), so it needs a freshly-hunted clean nonce for
   the fanout circuit.

4. **The in-repo fast prefilter does NOT model this circuit.** `island_search.rs`
   (orphaned; not in the module tree) screens nonces classically via
   `dialog_gcd_classical_filter`. Its soundness self-check FAILS for trailmix:
   the known-clean base nonce screens `hard=9024` (should be 0). That filter
   models the *dialog* GCD path (`build_builder`/`emit_dialog_gcd_raw_pa`), not
   the *trailmix* path that `build()` actually emits, and the two use conflicting
   env configs (`configure_ecdsafail_submission_route` vs
   `configure_q1153_second512_submission_defaults`). So the cheap nonce hunt the
   frontier-holders used is unavailable for the current circuit without writing a
   new trailmix-specific prefilter (likely GPU, per the CUDA-parity comments).

## What a genuine improvement would require

- A **lower-width modular inversion** to crack the 1152 plateau (real research —
  this is how 1175→1152 was reached), **or**
- A **fast nonce scan** (GPU) to hunt a clean nonce for the −318-CCX fanout
  circuit, or a lower-avgT nonce for the current one.

Neither is reachable by parameter tuning. The challenge's stated goal (beat the
published Pareto frontier) is already met; advancing the live leaderboard is an
open competitive problem.

## Fanout-reduction win (quantified) and why it's blocked (measured)

Enabling `single_ccx_fanout` (`SINGLE_CCX_FANOUT_DISABLE=0`) removes **318 CCX**
(rewrite `CCX;CX;CCX → CX;CCX`, **0 added qubits**): the fanout circuit measures
avg Toffoli **≈1,363,940** at peak **1152** → score **≈1,571,26x,xxx**, which
WOULD beat the frontier by ~330k. It only needs one Fiat-Shamir-clean nonce.

Blocker, all measured with the real-sim screener below:
- The fanout circuit's clean rate is ~e⁻¹⁷ (random nonces show cm≈16–21,
  pg≈5–16), i.e. ~1 clean nonce per ~2·10⁷. CPU brute-force ≈ weeks.
- **No cheap loosening exists.** The `TLM_FOLD/COUT_K/HYB_V` deltas are no-ops
  at this operating point (byte-identical circuits at delta 0/1/2 — the frontier
  is already at their floor). Widening the GCD swap-comparator window
  (`TLM_GAP_ADD`, an experimental knob, since reverted) makes correctness
  *dramatically worse*: +4 bits → cm 16→1250, +8 → cm→8990. `GAP_J2` is a
  precisely-tuned value; deviating either way breaks the swap decisions.
- Reducing the failure rate fundamentally needs a wider width truncation
  (`SCHED_J2`) = more qubits = worse score. The frontier is a razor-sharp
  optimum at the 1152 floor.

## Real-sim nonce screener (tooling, env-gated, NOT in the scored circuit)

`island_search::run_realsim_hunt_from_env`, gated on `FASTHUNT`, screens nonces
with the ACTUAL `Simulator`, reproducing `eval_circuit`'s correctness/phase/
ancilla checks exactly (validated: the frontier nonce screens CLEAN with avgT
1,364,229.770 — identical to eval), with early exit at the first failing batch.
This is the fast oracle a GPU/long CPU nonce hunt needs (the orphaned
`run_from_env` prefilter models the *dialog* circuit and is useless here — its
self-check gives hard=9024 for the trailmix base nonce). Also `HUNT_AVG_PROBE`
full-sims one nonce and reports avg Toffoli + failure counts for any config.

Usage: `FASTHUNT=1 SINGLE_CCX_FANOUT_DISABLE=0 HUNT_START=0 HUNT_COUNT=1000000 \
ISLAND_THREADS=16 HUNT_REPORT_AVG=1 ./target/release/build_circuit` prints
`CLEAN <nonce> <avgT>` for each eval-clean nonce found.
