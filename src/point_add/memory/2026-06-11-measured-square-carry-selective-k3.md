# Measured Square Carry and Selective K3

Date: 2026-06-11

## Implemented Levers

- `SQUARE_ROW_WINDOW_MEASURED_CARRY_CLEAR=1`
- `DIALOG_GCD_SELECTIVE_K3_STEP=<step>`
- optional prototype `DIALOG_GCD_SELECTIVE_K3_STEP2=<step>`

All levers remain default-off.

## Verified Structural Facts

- The 256-bit square self-test passes with the measured carry cleanup.
- The measured cleanup preserves the 1,221-qubit peak.
- Selective K3 forward/reverse evaluations left zero ancilla garbage.

## Best Measurements

- Measured cleanup only: 1,404,169.744 average Toffoli, 1,221 qubits, failed
  GCD island.
- Measured cleanup plus K3 step 240: 1,404,876.493 average Toffoli, 1,221
  qubits, failed GCD island.
- K3 step 0, nonce 175488: zero classical failures, one phase batch, zero
  ancilla garbage.

## Negative Evidence

- No K3-step-240 filter-clean nonce in the first 300,000 candidates.
- A second K3 shift raised the peak to 1,222 and was noncompetitive.
- Strict comparator filtering has known false negatives and found no survivor
  in large diagnostic sweeps.

## Next Step

Use a phase-aware filter or distributed nonce search for the lower-cost
step-240 route. Do not submit without a full 9,024-shot clean run.
