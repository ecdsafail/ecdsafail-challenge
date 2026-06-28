# Muqatta-Z Defense Module — ecdsa.fail Integration

**Author:** Peter Anari Otuke, University of Nairobi
**Framework:** Muqatta-Z Sovereign (Z-parameter formalism)

## Purpose

This defense module consumes ecdsa.fail leaderboard data as input to the
Muqatta-Z spectral framework, converting offensive circuit improvements
into calibrated defensive migration urgency.

The leaderboard score `C` (Toffoli x qubits) feeds the framework via:

```
Z15(t) = Z15(0) * (C(0) / C(t)) ^ gamma
```

where `gamma = 0.5` is a calibrated duality-elasticity constant. As the
leaderboard score drops (circuits improve), the ECDLP-PQC coupling term
increases, which raises the modeled optimal allocation toward PQC
migration, `w*(PQC)`.

## Key outputs

| Metric | Description |
|--------|-------------|
| `gamma` | Duality elasticity exponent applied to the score-improvement ratio |
| `w*(PQC)` | Modeled optimal share of migration effort allocated to PQC |
| `D*` | SegWit-PQ witness discount factor derived from witness-size ratios |
| `K` | Number of fixed-point iterations to converge under the Banach bound |
| `CR` | Spectral stability margin of the calibrated Z-matrix |

## How to reproduce

```bash
pip install numpy scipy
python src/point_add/memory/z_calibrate.py --score <score_from_score.json>
```

See `z_calibrate.py` in this directory for the full calibration
implementation, and `defense_approach.md` for the rationale behind the
approach.
