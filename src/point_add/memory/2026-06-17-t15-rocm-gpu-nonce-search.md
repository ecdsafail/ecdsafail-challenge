# T15: GPU Nonce-Suche auf ROCm (RX 6600 XT) (2026-06-17)

## Goal
Port / set up GPU-accelerated classical GCD prefilter nonce search for the AMD RX 6600 XT using ROCm/HIP (counterpart to T14 CUDA).

## Status after T15
- Created ROCm infrastructure mirroring T14 CUDA:
  - `tools/rocm/island.hip` : HIP kernel skeleton + usage
  - `tools/rocm/build.sh` : compile helper (gfx1032 for Navi 23 RX 6600 XT)
  - `nonce_hunt_t15.py` : updated hunt script with ROCm GPU mode + prepare shards
  - `find_clean_nonce_t15.py` : fast sampler
- Bins are current (classical_mismatch_count with T14/T15 trims)
- CPU prefilter works (16 workers possible on Ryzen 7 1700X)
- ROCm runtime partially present (rocm-smi works, GPU detected as Navi 23), but hipcc not installed yet (needs `sudo apt install hipcc ...`)
- Vulkan available via RADV but not used for compute here.
- No clean classical nonces found in prior T14 samples (~600+); density low under current trims.

## How to use on this machine
1. Install dev tools (requires sudo once):
   ```
   sudo apt-get install hipcc libhipblas-dev libhipcub-dev
   ```
2. Build:
   ```
   cd tools/rocm
   ./build.sh
   ```
3. Run large shard search:
   ```
   ./rocm_island_search --start 2150000000000000 --count 100000000 > survivors.txt
   ```
4. Verify survivors locally:
   ```
   python3 nonce_hunt_t15.py --max-tests 100 ...
   ```
   or directly:
   ```
   DIALOG_TAIL_NONCE=... ./target/release/classical_mismatch_count
   DIALOG_TAIL_NONCE=... cargo run --release -- --note "T15 candidate"
   ```

## Notes on port
- HIP is source-compatible with CUDA (small changes: hip/hip_runtime.h, hip* APIs, --offload-arch=gfx1032)
- Full kernel still needs the heavy lifting: device secp256k1 arith + exact GCD transcript replay with all HardReason checks.
- Keep parity with Rust reference (use small test ranges first).
- Target: millions of nonces tested quickly to find the rare clean islands for current trim set.

## Current metrics (base)
- score.json: 1,668,730,944 (T=1,428,708 Q=1168)
- Base nonce mm ~1 (early abort in filter)

## Next
- Install hipcc + compile + run large scale ROCm search
- Collect 0-mismatch candidates
- Full verify + potential submit if score improves
- If successful, document the found nonce and any score improvement

## Artifacts
- tools/rocm/island.hip
- tools/rocm/build.sh
- nonce_hunt_t15.py
- find_clean_nonce_t15.py
- This memory doc
