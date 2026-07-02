# ADR 0004 — Cross-validate against the source paper's reference circuits

**Status:** Accepted
**Date:** 2026-07-02

## Context

The analysis layer already proves the repo's arithmetic at two levels: abstract
z3 lemmas ([§1a–b](../scientific-value.md)) and bit-precise Kani proofs bound to
the real `alloy` U256 type (§1c). Both are *internal* — they prove this repo's
own implementation against a ground truth we also state.

The source paper (Babbush et al. 2026, arXiv:2603.28846v2; see
[ADR 0003](0003-ground-ecdlp-estimate-in-source-paper.md)) ships **reference**
in-place adder circuits in the kickmix `.kmx` format, including negative
controls. `iadd8.kmx`/`iadd64.kmx` are explicitly "a variant of the adder from
quant-ph/0410184" (Cuccaro et al.) — the same primitive this repo's arithmetic
core uses — and `iadd8_with_classical_offset_and_dirty_ancillae.kmx` is a
classical-addend, dirty-borrowed-ancilla adder structurally identical to this
repo's "quantum point += classical point" shape. These are an external oracle we
had not exploited.

## Decision

Add an independent, external cross-validation layer:

- `analysis/verify/kickmix_sim.py` — a spec-faithful `.kmx` simulator, a Python
  re-derivation of the kickmix semantics `src/sim.rs` implements (never sharing
  code with it, so agreement is meaningful).
- `analysis/verify/validate_reference_adders.py` — fuzzes the paper's reference
  adders through that simulator and asserts, deterministically (seeded): correct
  output (`r0 += r1`, or `r0 += r1 + carry`), clean workspace returned to `|0⟩`,
  dirty borrowed ancilla restored to their random input, and global phase `+1`.
  It also asserts the paper's three negative controls
  (`inc3_wrong_{order,phase,garbage}`) are **rejected**.
- Vendor only the tiny, human-readable reference artifacts under
  `analysis/verify/reference_circuits/` (with attribution) so the check is
  self-contained; gitignore the full multi-megabyte Zenodo dump (`/original`).
- Wire it into `analysis/run.sh` as stage 3/5 and document it as
  `scientific-value.md` §1d.

## Consequences

- The kickmix semantics the whole repo depends on (and that `eval_circuit` uses
  to score it) are now shown to reproduce the source paper's own published
  artifacts, and the phase-/ancilla-aware fuzz methodology is shown to actually
  reject buggy circuits — the paper's Appendix A.5 correctness argument,
  reproduced end-to-end.
- Determinism is preserved (fixed seed), consistent with
  [ADR 0001](0001-analysis-layer-isolated-from-score.md): the check reads only
  vendored/untracked reference data and never touches `score.json` or the
  editable circuit.
- Reverse-engineering the classical-offset variant established two facts worth
  recording: its `r2`/`q8` is a **carry-in** (`r0 += r1 + carry`), and its dirty
  borrowed ancilla are exactly the non-register qubits never reset by an `R`
  instruction (clean workspace = the `R`-targets).
- Cost: a dependency on the reference dump being present. Handled by vendoring
  the needed files; the validator degrades to a clean SKIP if neither the
  vendored copy nor the Zenodo dump is found.
