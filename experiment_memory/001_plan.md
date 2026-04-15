# Experiment 001 — Session Plan (2026-04-15)

Baseline confirmed: **18,506,238 Toffoli @ 3083 qubits** (commit 07b05f5).

## Budget (from hot-paths memory)
- Kaliski: ~11M (~60%) — 4 passes × 512 iters × ~5400 CCX/iter
- Mod muls: ~6M (~32%) — 5 mul calls, each ~1.17M
- Kaliski halve/double sandwich: ~1M (~5%)

## Dead ends (don't retry)
- `mod_add_qq` bit-0 parity flag uncompute (2-unknowns-1-eq)
- `with_lt/cmp_lt_into` first-maj ccx skip
- Batch-halve pair in mul unwind (flag uncompute costs as much as save)
- Skipping step 4 first CCX exploiting u[0]=1 invariant

## Ranked next experiments

### Tier 1 — moonshots (big structural wins)
1. **Gidney 2018 carry-lookahead n-Toffoli adder** (vs current Cuccaro 2n). Potentially 5-6M save across all adds. Substantial rewrite of `cuccaro_add` replacement.
2. **4→3 Kaliski passes**: rewrite pair 2 to reuse pair 1's state without forward/reverse. ~2.75M save. Needs algebraic rethink of `lam` uncompute path (not via second kaliski on -dx).
3. **2-bit windowed Kaliski iteration**: halve iteration count 512→256. Needs step-by-step re-derivation of invariants.

### Tier 2 — medium-risk
4. **Windowed mod_mul (Gidney 2019 / Ragavan-Gidney 2025)** with precomputed `{x, 2x, 3x}` table. ~260k per mul × 5 = 1.3M potential.
5. **Lazy modular reduction** in mul: let acc grow to 2p-3p across multiple adds, reduce less often. Need wider register and smarter final reduction.
6. **Fuse step 3+9 cswap cascades** in kaliski_iteration by absorbing swap into step 4's conditional sub/add. ~1030 CCX/iter = ~4M potential if clean.

### Tier 3 — small but definite
7. **Share ancillas across consecutive mod_double calls** in with_kal_inv's halve sandwich. Saves Clifford+qubits, likely 0 Toffoli.
8. **Specialized `add_const_sparse`** for c=2^32+977 via Gidney's in-place ripple. Investigate whether it beats the 2n Cuccaro.

## Attack plan
Start with Tier 1 moonshot #1 (Gidney lookahead adder) — highest-EV. If stuck, try #2 (4→3 Kaliski). Tier 2 fallbacks if Tier 1 fails. Small tiers only as warmups between moonshots.
