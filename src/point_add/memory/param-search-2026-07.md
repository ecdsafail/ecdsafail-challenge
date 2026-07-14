# Parameter-search leads (2026-07)

Notes from a single/double-knob sweep over the current committed circuit,
looking for a cheap score win. **None found** — recorded here so the next
contributor doesn't repeat it. Treat as a lead, not gospel: re-run
`./benchmark.sh` before relying on any line.

## Anchor (current committed circuit)

- **score 1,519,170,048** = avg-executed Toffoli **1,318,724** x peak qubits **1,152**
- 9,024/9,024 shots OK, phase- and ancilla-clean.

Method: every `TLM_*`/`DIALOG_*` default is set via `set_default_env`, so a
candidate is just an environment override of the *same* binaries — no rebuild.
Each run is the full official `./benchmark.sh` (validation + score).

## Results

`FAIL` = harness rejected it (breaks correctness/reversibility/phase, panic
rc=101). `no-op` = valid but byte-identical score.

| knob | value | result |
|------|-------|--------|
| `TLM_APPLY_INV_CSWAP_SKIP_LAST` | 2 (from 1) | FAIL |
| `TLM_APPLY_FWD_CSWAP_SKIP_LAST` | 3 (from 2) | FAIL |
| `TLM_APPLY_INV_S2_ZERO_LAST` | 2 | FAIL |
| `TLM_APPLY_FWD_S2_ZERO_LAST` | 2 | FAIL |
| `TLM_APPLY_ADD_SKIP_INV` | 1 | no-op |
| `TLM_APPLY_ADD_SKIP_FWD` | 1 | no-op |
| `TLM_TARGET_Q` | 1151 (from 1155) | FAIL |
| `TLM_FOLD_DELTA` | 1 (from 2) | FAIL |
| `TLM_COUT_K_DELTA` | 1 (from 2) | no-op |
| `TLM_HYB_V_DELTA` | 1 (from 2) | no-op |
| `TLM_TARGET_FFG_RESERVE` | 8, 7 (from 9) | FAIL |
| `TLM_TARGET_FOLD_RESERVE` | 3, 2 (from 4) | FAIL |
| `TLM_TARGET_FFG_RESERVE` + `TLM_TARGET_FOLD_RESERVE` | 8 + 3 | FAIL |
| `TLM_GCD_K_ADJUST` | -1 | no-op |
| `TLM_FFG_DELTA` | 1 (from 0) | no-op |
| `TLM_COUT_K_DELTA`+`TLM_HYB_V_DELTA`+`TLM_FFG_DELTA` | 1+1+1 | no-op |

## Takeaways

- **Qubit-reserve floor is hard.** Every downward move on `TLM_TARGET_Q`,
  `TLM_TARGET_FFG_RESERVE`, `TLM_TARGET_FOLD_RESERVE` fails validation — the
  scratch budget is already at the minimum this layout can run on. Shaving a
  qubit (worth ~1.3M score) needs a *structural* change, not a smaller global
  reserve. The per-call reserve override lists
  (`TLM_TARGET_FFG_CALL_RESERVES`, `TLM_TARGET_FOLD_CALL_RESERVES`) are the
  finer lever, but only at the specific call that binds the 1,152 peak — worth
  tracing (`TRACE_TLM_GCD_STEPS`) before touching.
- **Skip/margin knobs are non-binding at these values** — bumping them either
  changes nothing or breaks the last-iteration invariants they guard.
- Net: the current point is a tight local optimum for one-knob edits. Progress
  likely requires an algorithmic/structural change (fewer Toffolis in the
  modular-inverse inner loop, or a genuinely narrower qubit layout).

## Where the qubit peak actually lives

Traced with `TRACE_TLM_GCD_STEPS=1 TRACE_TLM_GCD_MIN_Q=1050 build_circuit`:
the GCD (Kaliski) steps top out at **active_max = 1,140**, yet the reported
**global_peak = 1,152**. So the ~12 qubits that bind the score's peak are
allocated *outside* the traced GCD loop — in a non-GCD phase
(fold/square/affine setup). Two consequences for anyone chasing the -1,152
lever:

- Retuning `TLM_TARGET_FFG_CALL_RESERVES` / GCD-side reserves cannot lower the
  global peak — the GCD never reaches it.
- The fold-side reserve *is* on the peak path, but `TLM_TARGET_FOLD_RESERVE`
  down from 4 already fails validation (it's at its floor too).

Shaving even one qubit therefore needs a structural change in that non-GCD
phase, not a reserve knob. A useful next step would be to add a peak-site trace
(the harness already tracks `circ.peak_qubits` / `phase_active_regions`) to name
the exact phase that holds 1,152 live qubits.
