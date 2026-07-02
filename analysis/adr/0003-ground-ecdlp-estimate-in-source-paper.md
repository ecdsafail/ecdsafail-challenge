# ADR 0003 — Ground the ECDLP estimate in the source paper's closed form

**Status:** Accepted (supersedes the ladder model of [ADR 0002](0002-derived-ecdlp-ladder-factor.md))
**Date:** 2026-07-02

## Context

[ADR 0002](0002-derived-ecdlp-ladder-factor.md) replaced a hand-picked
multiplier with a structurally *derived* ladder factor, `2(n+1) = 514` basic
double-and-add additions, and priced the (then-hypothetical) windowed lookup at
`~2^w` Toffoli. That was a reasonable first-principles model in the absence of
the source.

The source paper has since been located in `docs/`:

> Babbush, Zalcman, Gidney, Broughton, Khattar, Neven, Bergamaschi, Drake, Boneh,
> "Securing Elliptic Curve Cryptocurrencies against Quantum Vulnerabilities:
> Resource Estimates and Mitigations", Google Quantum AI, 2026
> (arXiv:2603.28846v2).

Its Appendix A gives the *exact* architecture this challenge is built on — a
kickmix point-addition circuit using measurement-based uncomputation — plus the
closed-form ECDLP cost and the paper's own zero-knowledge-proven bounds. This
supersedes the assumptions in ADR 0002 with authoritative facts.

## Decision

Rebuild `analysis/ecdlp_estimate.py` and the `cost_model.py` extrapolation on the
paper's closed form (Appendix A, eqs. A1–A3):

- `ECDLP_Toff = (PA_Toff + 3·2^w)·(2n/w − 4)` and `ECDLP_Qubits = PA_Qubits + w`,
  optimal window **w = 16** → **28** windowed point additions (not 514; the
  lookup is **3·2^w**, not `2^w`). Register width does **not** grow with the
  number of additions — the ladder reuses the accumulator/ancilla.
- Substitute this repo's **measured** point-addition metrics (`PA_Toff`,
  `PA_Qubits`, `PA_Toffoli-depth`) into the formula. Result: ~43.7M Toffoli at
  1,168 qubits, reaction-limited to ~5 minutes.
- Report the primitive-level comparison against the paper's ZK-proven bounds:
  this repo's PA (1,364,230 Toffoli · 1,152 qubits · 10.2M ops) is under the
  paper's **Low-Qubit** bound (≤ 2.7M / ≤ 1,175 / ≤ 17M) on all three axes, and
  the composed ECDLP undercuts both the paper's Low-Qubit (≤ 90M) and Low-Gate
  (≤ 70M) published Toffoli counts (~2.06× and ~1.60× fewer).
- Cite the paper as the primary reference in both scripts; keep Roetteler et al.
  2017 only as background for the general ladder structure.

## Consequences

- The full-attack estimate is now anchored to the challenge's own source and its
  published numbers, and the paper's formula reproduces its own bounds as a
  cross-check — a much stronger footing than the ADR 0002 model.
- The "beats the frontier" claim is now precise and primitive-level: an improved
  instance of the paper's ZK-proven point addition, not a vague comparison to
  RSA-2048.
- New explicit caveat beyond ADR 0002's completeness point: this repo's measured
  PA folds a **classical, compile-time** addend (constprop), whereas the paper's
  windowed PA loads `P[k]` from a **quantum** table. The `3·2^w` term prices the
  lookup, but residual constprop gains may not fully transfer — the composition
  assumes an addend-independent arithmetic core.
- The paper PDF itself is left untracked in `docs/` (large binary / third-party);
  it is cited by arXiv id rather than committed.
