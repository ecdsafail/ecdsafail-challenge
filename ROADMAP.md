# Roadmap

Tracked, actionable work for this repo. Each item links to its GitHub issue and
to the in-repo docs where the detail and rationale live. This file is an index,
not a second source of truth — decisions live in `analysis/adr/`, and the honest
list of what the analysis does **not** yet cover lives in
`analysis/scientific-value.md` (§2 "Key limitations", §Scope/honesty, §1d).

## Verification / analysis

- [ ] **Validate the reference QROM (unary-iteration table lookup)** —
  [#3](https://github.com/CaptainEmpower/ecdsafail-challenge/issues/3).
  Model measurement-based unary iteration so `table_lookup_3x3.kmx` (the `3·2^w`
  ladder lookup) validates alongside the adders.
  Detail: `analysis/scientific-value.md` §1d, ADR 0004.
- [ ] **Tier B: emit + measure the full ECDLP ladder** —
  [#4](https://github.com/CaptainEmpower/ecdsafail-challenge/issues/4).
  Replace the *derived* ECDLP cost (~43.7M Toffoli / 1168 qubits) with
  emitted-and-measured totals via streaming emission.
  Detail: `analysis/ecdlp_estimate.py`, ADR 0003.
- [ ] **Adder completeness / exceptional cases (Tier C)** —
  [#5](https://github.com/CaptainEmpower/ecdsafail-challenge/issues/5).
  Handle P==Q, P==−Q, ∞ (complete formulas or a negligibility proof) to turn the
  cost estimate into a verified attack.
  Detail: `analysis/scientific-value.md` §2, ADR 0003.

## Challenge score

- [ ] **Reduce Toffoli / peak-qubit liveness** —
  [#6](https://github.com/CaptainEmpower/ecdsafail-challenge/issues/6).
  Current ~1.57e9 (1,364,230 Toffoli × 1,152 qubits) already beats the published
  frontier; both factors are heavily hand-tuned. Editable path: `src/point_add/`.

## Code health

- [ ] **Zero out build warnings (~38)** —
  [#7](https://github.com/CaptainEmpower/ecdsafail-challenge/issues/7).
  Byte-identical to the score; treat `_`-prefixed/dead code as a missing-
  functionality hint (link or implement) rather than blindly deleting.

---

**Provenance.** This repo is a solution to the challenge from Babbush et al.
2026, *Securing Elliptic Curve Cryptocurrencies against Quantum Vulnerabilities*
(arXiv:2603.28846v2); see `analysis/adr/0003-*` and `0004-*`. Remotes: `origin`
is the working fork, `upstream` is the canonical challenge repo.
