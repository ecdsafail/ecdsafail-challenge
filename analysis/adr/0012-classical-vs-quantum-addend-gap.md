# ADR 0012 — The classical-vs-quantum-addend gap is negligible (Tier B, issue #27)

**Status:** Accepted — the Toffoli gap `ecdlp_estimate.py` flagged is measured
≤ 0.05% of PA; the register-overlap (width) part of #27 remains
**Date:** 2026-07-02

## Context

`analysis/ecdlp_estimate.py` composes the scored per-addition Toffoli (`PA_Toff`)
with the paper's windowed ladder, and flags one assumption:

> the paper's PA table lookup loads `P[k]` from a **quantum** window register; this
> repo's measured PA adds a **classical, compile-time** point, so its constprop
> pass (which folds addend-dependent gates) may not fully transfer to the windowed
> setting.

That "classical-vs-quantum-addend gap" is the keystone of issue #27, and it
propagated into `ladder_full.rs`'s `28·PA` (ADR 0011) as a stated assumption.

## Decision

Measure the gap directly rather than leave it asserted
(`src/point_add/constprop_gap.rs`, `#[cfg(test)]`, `#[ignore]` heavy): build the
circuit with the classical-addend optimizations toggled and count static
`CCX`/`CCZ`.

## Consequences

- **The gap is negligible — ≤ 0.05% of PA.** Structural reason (`coord_addsub`,
  `trailmix_ludicrous/ec_add.rs`): to consume the classical addend the circuit
  allocates a fresh **qubit** register, loads the addend into it with `x_if_bit`
  (Clifford `X`), runs an **uncontrolled quantum-quantum** vented mod-add/sub (a
  full Cuccaro-class adder — the same Toffoli cost a *quantum* addend pays), then
  unloads. So the arithmetic Toffoli count is addend-**value**-independent; the
  classical value only drives the (Clifford) load/unload.
- **Measured:** the only addend-value-dependent optimization is the peephole
  constant-propagation pass (`CONSTPROP_DISABLE`), worth **770 Toffoli
  (0.05%)** of the 1,446,685 static PA. The direct-constant-arithmetic knobs
  (`SECP_DIRECT_CONST_ARITH`, `KAL_DIRECT_CONST_WALKS`) are **inert** for the
  trailmix build path (identical Toffoli) — confirming the scored circuit already
  takes the load-into-qubits path.
- **So `ladder_full.rs`'s `28·PA` and `ecdlp_estimate.py`'s headline already
  reflect the quantum-addend *arithmetic* cost.** The caveat is quantified as
  essentially moot for Toffoli.
- **What #27 still needs** (deferred): the width / register-overlap question —
  whether the QROM-provided addend register and lookup ancilla reuse PA's
  workspace behind `ECDLP_Qubits = PA_Qubits + w` (A2) — plus functional
  correctness of a QROM-fed add. `coord_addsub` *frees* its addend register within
  the addition; a quantum-addend ladder keeps the addend live (from the QROM)
  across the add, so the width analysis differs even though the Toffoli does not.
  That, and issue #28's mid-ladder residual, remain the Tier B build.
- Consistent with [ADR 0001](0001-analysis-layer-isolated-from-score.md): the
  harness is `#[cfg(test)]`, never compiled into `build_circuit`; the scored
  circuit is byte-identical (`ops.bin` SHA unchanged).
