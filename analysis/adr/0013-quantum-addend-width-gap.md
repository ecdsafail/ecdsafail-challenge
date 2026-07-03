# ADR 0013 — The quantum-addend WIDTH gap is real: A2's `+w` undercounts this PA (Tier B, issue #27)

**Status:** Accepted — resolves the width/register-overlap half of #27 that
[ADR 0012](0012-classical-vs-quantum-addend-gap.md) deferred. Unlike the Toffoli
gap (≤ 0.05%, negligible), the **width** gap is material: a faithful
quantum-addend port of this repo's PA needs `PA_Qubits + (256..512) + w`, not
`PA_Qubits + w`.
**Date:** 2026-07-03

## Context

ADR 0012 closed the *Toffoli* half of #27 (the arithmetic is addend-value-
independent) and explicitly deferred the *width* half:

> the width / register-overlap question — whether the QROM-provided addend
> register and lookup ancilla reuse PA's workspace behind `ECDLP_Qubits =
> PA_Qubits + w` (A2).

`ladder_full.rs` (ADR 0011) reports `total_qubits_a2 = pa_qubits + w` but flags
it: *"the true overlap needs the quantum-addend PA … emitting the lookup on
disjoint ids would OVER-count width."* This ADR turns that flag into a
measurement.

## Decision

Establish the quantum-addend width by construction plus one measured fact
(`src/point_add/addend_width.rs`, `#[cfg(test)]`, `#[ignore]`): build the scored
circuit and construct the resident-addend port by shifting every qubit id up by
the addend width `Δ` (placing a held `Δ`-qubit addend register at ids `[0, Δ)`),
then read the port's **allocation span** via `analyze_ops` (= `max qubit id + 1`).
The span is `PA_span + Δ` by construction; the *peak* consequence follows from the
separately-measured tight-peak fact below, not from `analyze_ops` alone (which is
a span, not a live-peak, metric).

## Consequences

- **Measured profile** (`TRACE_PHASE_ACTIVE=1 TRACE_TLM_PROFILE=1 cargo run
  --release --bin build_circuit`): the scored peak is **1152** at the GCD apply
  (`tlm_multiply_gcd_reverse_body`); the coordinate steps run at **1026**, and
  that 1026 *already includes* the 256-qubit classical-addend temp. So the
  classical addend is resident only **off-peak**: `coord_addsub` loads it into a
  fresh temp per coordinate step and frees it within the step, so it never
  coexists with the GCD scratch.

- **A quantum (QROM) addend cannot do this.** `ox` is consumed at steps 3/7/15
  and `oy` at 4/14 (`ec_add.rs`) — both straddle both GCD passes and the square,
  i.e. the peak. A single QROM read of `P[k]` must stay resident *across the
  peak*, where it **cannot overlap the GCD scratch**: the preserved addend
  (needed *after* the GCD) and the GCD scratch (in use *during* it) are both live
  at the peak. The scored peak is *tight* — max qubit id (`analyze_ops` = 1152 =
  `score.json` qubits) equals the profiler's peak active (1152), so the free list
  is empty at the peak and there is no slot to reuse.

- **Port allocation spans** (`quantum_addend_width_gap`; the span delta is a
  peak delta given the tight peak above):
  | model | port span | vs PA |
  |---|---|---|
  | scored PA (classical addend) | 1152 | — |
  | hold one coordinate (lower bound) | **1408** | + 256 |
  | hold `P[k]=(x,y)` | **1664** | + 512 |
  A held addend register adds its **full** width to the span; because the peak is
  tight (free list empty), it cannot overlap the GCD scratch, so this is also the
  peak growth.

- **So A2's `+w` undercounts THIS PA.** `ECDLP_Qubits = PA_Qubits + w` = 1168;
  the measured quantum-addend port is `PA_Qubits + (256..512) + w` ≈ **1408..1664
  (+w)**. The width caveat is the **opposite disposition** from the Toffoli gap:
  real, not negligible.

- **Why the paper's A2 still holds — for the paper's PA.** The paper's ZK
  `PA_Qubits` bound (1175 low-qubit / 1425 low-gate) already prices a *resident
  quantum addend* into a tighter arithmetic core. This repo's product-min PA
  spent its width budget on the GCD and stayed under bound precisely by keeping
  the addend **classical** — an advantage a faithful quantum-addend port would
  erase (indeed reverse: 1408..1664 > 1175/1425). The `+w`-only claim transfers
  to the paper's co-designed PA, not to this repo's.

- **Estimate impact.** `ecdlp_estimate.py`'s qubit headline uses A2 (`PA_Qubits +
  w`) — correct as a restatement of the *paper's* bound, now annotated that a
  verified quantum-addend port of *this* PA carries a `+256..512` width caveat
  (the Toffoli headline is unaffected — ADR 0012). This sharpens, not overturns,
  the estimate: the attack's qubit cost is dominated by whether the addend is
  co-designed into the core (paper) or bolted on (naive port).

- **What #27 still needs** (deferred): functional correctness of a QROM-fed
  quantum-addend add, and a co-designed core that hosts the resident addend
  within budget (the paper's route). Issue #28's mid-ladder ∞/`dx=0` residual
  lands in that same quantum-addend testbed.

- Consistent with [ADR 0001](0001-analysis-layer-isolated-from-score.md): the
  harness is `#[cfg(test)]`, never compiled into `build_circuit`; the scored
  circuit is byte-identical (`ops.bin` SHA unchanged).
