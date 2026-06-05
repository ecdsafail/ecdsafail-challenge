# 2026-06-05 1320q apply teardown

Current promoted base:

- `1382q x 1,514,221T = 2,092,653,422`
- Peak phases are the apply materialized-special chunk raw add/sub pair.
- The next active floor after apply is `1320q` in the compressed-block
  tobitvector phases; ROUND84 is `1309q`.

Peak anatomy after phase-label instrumentation:

- `dialog_gcd_apply_chunk_{add,sub}_final_ripple`: `1382q`
- `dialog_gcd_apply_chunk_{add,sub}_boundary_clear`: `1381q`
- `dialog_gcd_apply_chunk_{add,sub}_ripple`: `1326q`
- `dialog_gcd_compressed_block_tobitvector_*`: `1320q`

Structural teardown prototype:

- `DIALOG_GCD_APPLY_FINAL_LOWQ=1`
  - Replaces only the final apply chunk's fast low-to-ext ripple with the
    no-carry Cuccaro form.
  - Drops final ripple from `1382q`; costs about `+44,376T`.
- `DIALOG_GCD_APPLY_BOUNDARY_SPLIT=<k>`
  - Splits the apply boundary-prefix comparator into two fast windows and
    recomputes the lower window to clear the boundary carry.
  - With cuts `54/108/162`, the minimum 1320-safe split is `53`.
- Rebalanced cuts:
  - `DIALOG_GCD_APPLY_CHUNKED_F_CUT=54`
  - `DIALOG_GCD_APPLY_CHUNKED_F_CUT2=108`
  - `DIALOG_GCD_APPLY_CHUNKED_F_CUT3=162`
  - These push the ordinary non-final chunk ripple under the `1320q` floor.

Best traced under-budget seed:

```text
DIALOG_GCD_APPLY_FINAL_LOWQ=1
DIALOG_GCD_APPLY_BOUNDARY_SPLIT=53
DIALOG_GCD_APPLY_CHUNKED_F_CUT=54
DIALOG_GCD_APPLY_CHUNKED_F_CUT2=108
DIALOG_GCD_APPLY_CHUNKED_F_CUT3=162
KAL_DOUBLE_CARRY_TRUNC_W=22

active=258
peak=1320
emitted_toffoli=1,584,911
score_if_clean=2,092,082,520
```

Validation status on inherited nonce:

- `active=258`, split 53, cuts `54/108/162`, no extra T knob:
  `8 classical`, `4 phase` failures; `1,585,945T` (602T over 1320 budget).
- Add `KAL_DOUBLE_CARRY_TRUNC_W=22`:
  `5 classical`, `2 phase` failures; `1,584,911T` (under budget).
- Nonce prefix `0..17` on the double22 seed found no clean island; best observed
  was nonce `6` with `4 classical`, `4 phase`.

Interpretation:

- The current promising route is no longer brute force over arbitrary
  parameters. It is a structural 1320q apply teardown plus a small island
  problem.
- A large nonce search should use a faster classical/phase filter. Full trusted
  eval is too slow for wide nonce scans and appends many `results.tsv` rows.
- The next non-bruteforce direction is to reduce the residual phase/classical
  failures at the double22 seed, or recover another exact T shave from the
  split comparator so the lower-risk no-double22 active-258 seed can be used.

Follow-up structural improvement:

- Hosted boundary split:
  - When `DIALOG_GCD_APPLY_BOUNDARY_SPLIT` lands exactly on a retained apply
    boundary carry, that carry can host the high-window comparator carry-in.
  - This avoids the extra low-window prepass and boundary qubit used by the
    generic split helper.
  - With four custom chunks and split at the first cut, e.g. cuts `54/108/162`
    and split `54`, the candidate is `1320q` and `1,558,597T`.
  - Cut rebalance down to `49/98/147` stays `1320q`; `48/96/144` crosses the
    final low-q ripple wall at `1323q`.

Best current structural seed:

```text
DIALOG_GCD_APPLY_FINAL_LOWQ=1
DIALOG_GCD_APPLY_BOUNDARY_SPLIT=100
DIALOG_GCD_APPLY_CHUNKED_F_BLOCKS=5
DIALOG_GCD_APPLY_CHUNKED_F_CUSTOM5=1
DIALOG_GCD_APPLY_CHUNKED_F_CUT=50
DIALOG_GCD_APPLY_CHUNKED_F_CUT2=100
DIALOG_GCD_APPLY_CHUNKED_F_CUT3=150
DIALOG_GCD_APPLY_CHUNKED_F_CUT4=200
DIALOG_GCD_ACTIVE_ITERATIONS=260

peak=1320
emitted_toffoli=1,565,417
score_if_clean=2,066,350,440
```

Validation status:

- Inherited nonce `3577`: `2 classical`, `0 phase`, `0 ancilla`.
- Bounded nonce sweep:
  - `0..19` found near misses but no clean island.
  - `60..108` found clean nonce `108`.
- Clean default route:
  - nonce `108`
  - `0 classical`, `0 phase`, `0 ancilla`
  - score `2,066,350,440 = 1320q * 1,565,417T`

Interpretation update:

- The 1320q target is structurally solved with significant T headroom.
- The old phase problem is gone at active 260 with nonce `108`.
- Active 261 still fits `1320q` but regresses to `2 classical`, `2 phase`.
- Active 262 crosses to `1327q`; the current active-row frontier is 260/261.
