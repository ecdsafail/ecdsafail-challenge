# HMR diagnostic finding

A direct HMR-sequence comparison between:
- `kaliski_iteration(..., iter_idx = 0)`
- and `kaliski_iteration_bulk_prefix3(..., iter_idx = 0)`

shows a large mismatch.

## Counts
- generic step-0 HMR count: **1027**
- specialized step-0 HMR count: **768**

## First operand divergence
- generic first HMR target qubit: `QubitId(1283)`
- specialized first HMR target qubit: `QubitId(1285)`

## Interpretation
This is exactly the kind of structural mismatch that can produce a coherent
phase bug under measurement-based uncompute:
- the specialized step is not emitting the same HMR sequence as the generic
  step,
- so even if the classical state update is equivalent, the phase history is not.

This is now the strongest concrete low-level explanation for the remaining phase
bug.
