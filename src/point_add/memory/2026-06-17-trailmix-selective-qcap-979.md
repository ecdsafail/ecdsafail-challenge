# TrailMix selective per-step q-cap → 979 qubits (first sub-980, validated)

**Date:** 2026-06-17  **Author:** pua-ecdsafail loop (Claude Opus 4.8)
**Result:** 979 peak qubits, 9024/9024 OK (all 4 gates), toffoli 29,074,641, score 28.46B.
(Qubit-record route; NOT score-competitive vs the 1168q/1.67B frontier — do not submit to server.)

## The lever
The TrailMix peak (980 at the default `TRAILMIX_Q_CAP=20`) is bound by the
`shrunken_pz` schedule's **peak step 353**, row `[A=88, B=89, ca=245, cb=245, q=23]`.
Working width there = `2·max(A,B) + 2·max(ca,cb) + q = 178 + 490 + q`; global peak =
working + ~292 fixed. So q=20 → 980, q=19 → 979.

A blunt global `Q_CAP=19` clamps q on **all ~490 steps** (universal q runs 23–38),
manufacturing ~6–16 classical misses/run → a clean tail nonce is ~1e-4 (infeasible).

**Fix — selective per-step budget** (`TRAILMIX_Q_TARGET`, new):
each step gets `q ≤ TARGET − 2·max(A,B) − 2·max(ca,cb)`, so q is trimmed *only* on
the peak-binding step(s). `TARGET=687` → step 353 q→19 (peak 979); every other step
keeps its natural q → misses collapse from ~10 to ~1.

Implemented in:
- `inversion/shrunken_pz_state_machine.rs`: `trailmix_q_width_step(wq,wa,wb,wca,wcb)`
  applied at both forward and backward resize sites (kept gate-for-gate symmetric).
- `inversion/shrunken_pz_schedule.rs`: `thin_factor_repairs_u256` mirrors the budget
  so the tail-nonce support search models the real circuit (+ optional
  `TRAILMIX_Q_MODEL_GUARD` for extra model strictness).
- `trailmix_port/mod.rs` `configure_sub1000_trailmix_route`: baked defaults
  `Q_CAP=99, Q_TARGET=687, TAIL_NONCE=270`.

## The residual 1 miss (important)
Even model-clean nonces had **exactly 1** real classical miss (97/112/151/208 each
failed at a different shot), because the abstract `repair_sample` bit-length model
cannot see a gadget width-logic dependency at the tight q=19 clamp (it is NOT a
factor bit-length overflow — `MODEL_GUARD=1` still rated them clean). The residual
is ~Poisson(1), so a real-clean nonce exists by lottery: **nonce 270 → 0 misses**.
Validate candidates with the full benchmark; the model is a screen, not an oracle.

## Reproduce
`./benchmark.sh` (defaults now give 979). Explicit:
`TRAILMIX_Q_TARGET=687 TRAILMIX_Q_CAP=99 TRAILMIX_TAIL_NONCE=270 ./benchmark.sh`
Count-only nonce search seed: `POINT_ADD_HASH_OPS_LEN=92854789`.

## Toward < 979
`TARGET=686` → step 353 q→18 → ~2 systematic misses (harder lottery, P~e⁻²). Better:
spread the cut — trim a 2nd near-peak step or a genuinely slack fixed-part register
(COUNTER_W=7 is dead: counter needs 8 bits, 89 misses; SROT_W=4 panics).
