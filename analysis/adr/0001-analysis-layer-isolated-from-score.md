# ADR 0001 — Analysis layer is isolated from the scored circuit

**Status:** Accepted
**Date:** 2026-07-02

## Context

The challenge score is `round(avg_toffoli_per_shot) * qubits`, produced by a
deterministic run of the evaluator (`src/bin/eval_circuit.rs`) over the circuit
built from `src/point_add/` (the only `editablePaths`). The `analysis/` layer
adds verification (z3, Kani) and physical cost modelling on top of that result.

Any tooling that shared code paths, output files, or build flags with the scored
circuit could silently perturb the score — or worse, tempt us to back-fill a
metric from an analysis script rather than a measured run. The global standard
"don't generate any result; only trust results from deterministic runs" applies
directly here.

## Decision

The analysis layer is strictly downstream and side-effect-free with respect to
scoring:

- It **reads** `score.json` and `depth.json`; it never writes them. `depth.json`
  is produced by a standalone binary (`src/bin/depth_report.rs`) that measures
  `ops.bin` via `circuit::analyze_depth` and does not run the simulator or touch
  `score.json`.
- Analysis lives outside `editablePaths`, so nothing under `analysis/` can change
  the circuit or the score.
- Kani harnesses sit behind `#[cfg(kani)]`; the normal build and `benchmark.sh`
  never compile them.
- Every published number is emitted by a deterministic run (measurement, z3, or
  Kani) or is an explicitly printed assumption — never hand-asserted.

## Consequences

- The score remains reproducible from the circuit alone; analysis can evolve
  freely without leaderboard risk.
- Analysis outputs may lag the circuit if `score.json`/`depth.json` are stale;
  regenerating them is a prerequisite, tracked as an operational step.
- Derived (non-measured) estimates must be labelled as such — see
  [ADR 0002](0002-derived-ecdlp-ladder-factor.md).
