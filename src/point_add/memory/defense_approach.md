# Defense Module Approach Notes

## What this does

The Muqatta-Z defense module provides the complementary side of the
ecdsa.fail equation. While the challenge optimizes attack cost (fewer
Toffoli gates and qubits = cheaper attack), this module quantifies the
defense response that each circuit improvement implies for Bitcoin's
post-quantum migration planning.

## Core relationship

Every improvement to the leaderboard's circuit score increases the
modeled value of PQC migration. The relationship used here is:

```
defense urgency ~ (improvement ratio) ^ gamma,  gamma = 0.5
```

This is not adversarial framing — better attack-cost estimates motivate
better-calibrated defense investment. The module quantifies how much
defense urgency each circuit improvement implies.

## Integration with circuit optimization

The defense module does not modify the circuit code in
`src/point_add/`. It reads `score.json` after a `cargo run --release`
(or `ecdsafail run`) and recalibrates:

- `w*(PQC)`: modeled share of migration effort allocated to PQC
- `D*`: SegWit-PQ witness discount factor
- `K`: number of fixed-point iterations to convergence

## Reproducibility

Given a fixed score input, all outputs are deterministic:

- numpy seed: 42
- Z-matrix: hardcoded 5x5 (see `z_calibrate.py`)
- Eigenvalues reproduce to 6 decimal places across numpy/scipy versions
  tested (numpy 1.26.2, scipy 1.12.0)

## Caveats

This is a parametric modeling exercise, not a peer-reviewed empirical
result. The Z-matrix entries and `gamma` are chosen, not derived from
external data; treat the urgency outputs as illustrative sensitivity
analysis rather than calibrated forecasts.
