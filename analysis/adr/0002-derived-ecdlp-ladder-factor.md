# ADR 0002 — Derive the full-ECDLP ladder factor instead of hand-picking it

**Status:** Accepted — ladder model superseded by [ADR 0003](0003-ground-ecdlp-estimate-in-source-paper.md)
**Date:** 2026-07-02

> **Superseded (2026-07-02):** the source paper (arXiv:2603.28846v2, `docs/`) was
> located and gives the exact closed-form ladder cost. The decision to *derive*
> rather than hand-pick the factor still stands; the specific `2(n+1)` model and
> `~2^w` lookup here are replaced by the paper's `(PA+3·2^w)(2n/w−4)`, w=16. See
> ADR 0003.

## Context

This repo's circuit is exactly **one** elliptic-curve point addition (quantum
point `+=` a classical, compile-time point). A full Shor attack on the
secp256k1 ECDLP is many such additions arranged as a double-and-add ladder
(Roetteler, Naehrig, Svore, Lauter 2017, arXiv:1706.06752).

To extrapolate from the measured single-addition cost to a full-attack cost,
`cost_model.py` originally multiplied by a hand-picked constant,
`ecdlp_point_additions = 1600`. That number was an unexplained placeholder: it
was not tied to the algorithm's structure and, as it turned out, overestimated
the addition count by roughly 3×.

We considered three tiers of remediation:

- **Tier A — derived estimate:** compute the ladder factor from the algorithm's
  structure and compose it with the measured per-addition metrics.
- **Tier B — emitted + measured full circuit:** build the ladder in the op
  stream and measure it (streaming emission; ~5×10⁹ ops, full simulation
  infeasible at 9024 shots).
- **Tier C — verified attack:** additionally solve completeness (exceptional
  cases: `P==Q`, `P==-Q`, `∞`) with complete addition formulas.

## Decision

Adopt **Tier A**. Replace the hand-picked multiplier with a structurally derived
count and compose it with the measured primitive, in a dedicated script
`analysis/ecdlp_estimate.py`:

- Basic double-and-add uses `2(n+1) = 514` controlled additions for `n = 256`
  (two scalar registers of `n+1` bits). `cost_model.py` now derives the same
  `2(n+1)` factor rather than asserting `1600`.
- Full-attack Toffoli = `514 × 1.36M ≈ 7.0×10⁸`; composed Toffoli-depth
  ≈ `514 × 1.08M` (the accumulator serializes the additions). Windowed variants
  (`w = 4, 8`) are reported as a sensitivity, with `w = 1` as the headline.
- The per-addition inputs (Toffoli, Toffoli-depth, qubits) are read from
  `score.json` + `depth.json`; every other quantity is either derived from the
  structure or a printed assumption, consistent with
  [ADR 0001](0001-analysis-layer-isolated-from-score.md).

## Consequences

- The full-attack estimate is now traceable to the algorithm and the measured
  primitive; changing `n`, the window, or an assumption re-derives it.
- The headline drops from the old `~2.2×10⁹` to `~7.0×10⁸` Toffoli — the correct
  order relative to Gidney–Ekerå's RSA-2048 estimate (`~3×10⁹`).
- **This is a cost estimate, not a verified attack.** Two things remain
  explicitly unmodelled and are the boundary to Tier B/C:
  - **Completeness:** the adder is the incomplete affine formula; a correct
    attack needs exception handling. Priced via `completeness_overhead`
    (default `1.0` = exceptions assumed negligible, per Roetteler).
  - **Register-file width:** the full algorithm needs more than one addition's
    width; the reported peak reuses one addition's ancilla.
- Numbers are derived, not emitted+measured; a measured full circuit is Tier B.
