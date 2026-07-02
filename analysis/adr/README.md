# Architecture Decision Records

Records of the significant decisions behind the `analysis/` scientific-rigor
layer. Each ADR is immutable once **Accepted**; a later decision that changes
course supersedes it with a new record rather than editing history.

Format: Status · Context · Decision · Consequences (lightweight MADR).

| ADR | Title | Status |
|---|---|---|
| [0001](0001-analysis-layer-isolated-from-score.md) | Analysis layer is isolated from the scored circuit | Accepted |
| [0002](0002-derived-ecdlp-ladder-factor.md) | Derive the full-ECDLP ladder factor instead of hand-picking it | Accepted (ladder model superseded by 0003) |
| [0003](0003-ground-ecdlp-estimate-in-source-paper.md) | Ground the ECDLP estimate in the source paper's closed form | Accepted |
| [0004](0004-cross-validate-against-reference-circuits.md) | Cross-validate against the source paper's reference circuits | Accepted |
| [0005](0005-validate-lookup-by-construction.md) | Validate the ladder lookup primitive by construction | Accepted |
| [0006](0006-adder-completeness-approach.md) | Approach to adder completeness (cost estimate → verified attack) | Accepted (Path A viable per gating experiment) |
| [0007](0007-tier-b-measured-ladder.md) | Tier B: measuring the full-ECDLP ladder | Accepted (first increment; QROM/QFT deferred) |
| [0008](0008-empirical-completeness-collision-rate.md) | Empirically validate (and sharpen) the completeness collision rate | Accepted (equidistribution validated; zero-window ∞ term dominant) |
| [0010](0010-measured-windowed-lookup-cost.md) | Measure the windowed-lookup (QROM) cost (Tier B, issue #4) | Accepted (lookup term grounded; end-to-end ladder deferred) |
