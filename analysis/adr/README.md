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
