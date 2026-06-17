# T14: GPU-accelerated Nonce-Suche (2026-06-17)

## Status
- Classical prefilter updated for current trims (active=258, body trims with 3s, pa9024=1 margin=0, width_margin=10 slope=1017, tobit cswap+shift=1 etc.)
- Base now ~1 (early abort) or 22 full.
- CPU classical_mismatch_count now ~2s for bad nonces (early exit on first mismatch), full ~20s for potential cleans.
- 500 random+ dense samples: 0 clean classical. Density appears << 1/500 under current trims.
- CPU hunt impractical for large scale; GPU essential.

## Deliverables for T14
- Updated `src/bin/classical_mismatch_count.rs` (T14 trims + early-abort for search speed)
- `nonce_hunt_t14.py` : parallel CPU prefilter + full eval, explicit GPU workflow notes
- `tools/cuda/island.cu` : CUDA search skeleton + interface doc (parity with Rust filter)
- `gpu-src/CudaBrainSecp/EcdsaFailFilter.cu` : placeholder for per-step width / trim logic
- Rebuilt bins

## GPU acceleration notes
- Target: run filter on GPU at 3k-4k nonces/s (full 9024-shot zero-phase)
- On host with GPU: build CUDA, shard nonces (1M chunks), run disjoint searches
- Survivors (mm==0 classical) transferred back for trusted `cargo run --release` + `eval_circuit`
- Must keep CUDA in sync with Rust trims / schedule / HardReason logic (use parity jobs over small nonce/shot ranges)
- See prior memory/2026-06-13-q1192-wmi-cuda-search.md for throughput, job arrays, verification

## Next steps (for GPU host or follow-up)
- Run large scale on GPU (millions to billions nonces)
- Collect clean classical candidates
- Local full verification + score update if better than current ~1.67e9
- If clean found, bake nonce or submit

## Current best known (from score.json)
Toffoli ~1,428,708 @ 1168 qubits → score 1,668,730,944

No clean nonce found locally in T14 samples. GPU scale required.
