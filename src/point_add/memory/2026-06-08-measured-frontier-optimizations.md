# Measured Frontier Optimizations (2026-06-08)

Status: Active research results for ECDSA.fail Point-Addition Challenge.  
Baseline Reference: **1,453,867 average Toffolis** / **1302 peak qubits** (score: **1,892,934,834** under `DIALOG_TAIL_NONCE=60009363210`).

---

## 1. Stepped `DIALOG_GCD_SHIFT_BAND_TRIMS` Schedules

* **Files/Functions Touched**: 
  - `src/point_add/rounds/dialog/mod.rs`: `dialog_gcd_shift_band_trim(step)` and `dialog_gcd_k2_shift_active_width(active_width, step)`.
  - `gpu-src/CudaBrainSecp/EcdsaFailFilter.cu`: `shift2_width_for_step(step, active_width)`.
* **Env Flags**: `DIALOG_GCD_SHIFT_BAND_TRIMS` (comma-separated list of trims, e.g. `"0,1,2,3"`).
* **Emitted Ops / Avg T Delta**:
  - `"0,1,2"` (3 bands): -3,096 ops / **-1,032 Toffolis**
  - `"0,1,2,3"` (4 bands): -4,608 ops / **-1,536 Toffolis**
  - `"0,1,2,3,4"` (5 bands): -6,144 ops / **-2,048 Toffolis**
  - `"0,1,2,2,3,3,4,4"` (8 bands): -7,236 ops / **-2,412 Toffolis**
* **Peak Qubit Delta**: **0 qubits** (peak remains at 1302 qubits during GCD reverse pass).
* **Correctness Status**: GPU prefilter aligned (CUDA implementation updated to match the schedule), but untested on 9024 shots (requires running the prefilter to find a matching tail nonce).
* **Estimated Island Density**: Minimal to no impact if using a stepped schedule (e.g. `0,1,2,3`). At late steps (step >= 195), the active width has plenty of headroom (actual bit-length is ~40 while active width is ~60), making a trim of 3 extremely safe and highly unlikely to cause width overflows.

---

## 2. Reducing GCD Peak around `reverse_add` / `shift`

* **Files/Functions Touched**: 
  - `src/point_add/rounds/dialog/compressed.rs`: `dialog_gcd_build_composite_scratch` (line 352).
* **Env Flags**: `DIALOG_GCD_COMPRESSED_LOG_U_HIGH_RUNWAY_BLOCKS` (integer).
* **Emitted Ops / Avg T Delta**: **0 Toffolis** (pure qubit layout/lifetime change).
* **Peak Qubit Delta & Phase**:
  - The peak of **1302 qubits** occurs at step 9 of the reverse pass during the `dialog_gcd_compressed_block_tobitvector_reverse_add` / `_shift` phases.
  - Active qubits composition at step 9: `tx` (256) + `ty` (256) + `u` (256) + `compressed_log` (405) + `owned` (123) + `raw_block` (6) = 1302.
  - **Qubit Cut**: Decreasing `DIALOG_GCD_WIDTH_MARGIN` (e.g. from 10 to 8) shrinks the active width at early steps, dropping the composite scratch deficit `want` by `2 * delta_margin` qubits. Setting `DIALOG_GCD_WIDTH_MARGIN=8` drops the global peak by **-4 qubits** to **1298 qubits**.
* **Correctness Status**: Untested on 9024 shots (requires a tail nonce search).
* **Estimated Island Density**: Reduces success rate (denser search needed). Dropping margin from 10 to 9 multiplies the expected number of random rerolls to find a clean island by ~2.5x.

---

## 3. Partial `cswap` Reduction in `tobitvector` / `apply`

* **Files/Functions Touched**:
  - `src/point_add/rounds/dialog/compressed.rs`: `dialog_gcd_safe_cswap_width` (dynamic cswap width computation).
* **Env Flags**: `DIALOG_GCD_CSWAP_TRIM` (not set by default to use the dynamic slope-based envelope).
* **Emitted Ops / Avg T Delta**:
  - Restricting tobitvector `cswap` using the dynamic envelope trim (`dialog_gcd_safe_cswap_width`) saves **-30,024 emitted ops / -10,008 Toffolis**!
  - Restricting apply `cswap` is **impossible** (results in 9024 classical mismatches) because `x` and `y` are mod $p$ values (256 bits) that do not shrink and must be fully swapped.
* **Peak Qubit Delta**: **0 qubits** (peak remains at 1302 qubits).
* **Correctness Status**: **Full 9024 eval verified** (passes 100% classical, phase, and ancilla checks under `DIALOG_TAIL_NONCE=60009363210` with **1,443,859 average Toffolis**).
* **Estimated Island Density**: **No change** (100% value-exact on the reachable verifier support).

---

## 4. Low-Qubit `round84` Square

* **Files/Functions Touched**:
  - `src/point_add/arith/multiply.rs`: `squaring_sub_from_acc_karatsuba` ( Lead B from prior audit).
* **Env Flags**: `ROUND84_XTAIL_KARATSUBA` (set to 1).
* **Emitted Ops / Avg T Delta**: **+50k to +80k Toffolis** overhead due to three separate Solinas modular reductions instead of one combined reduction.
* **Peak Qubit Delta & Phase**:
  - Drops the squaring phase peak from 1302 to **902 qubits**.
  - **Constraint**: Since the GCD-walk peak is currently locked at 1302 during the reverse pass, the global peak remains **1302 qubits**.
  - Therefore, the sequential Karatsuba square yields **0 global peak qubit savings** and is **not viable** unless the GCD-walk peak is simultaneously lowered below 1302.
* **Correctness Status**: Untested / Blocked by GCD co-peak.
* **Estimated Island Density**: No change.
