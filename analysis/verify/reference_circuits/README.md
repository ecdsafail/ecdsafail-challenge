# Reference kickmix circuits (vendored)

These `.kmx` circuit files and `inc3_test_cases.txt` are **unmodified** reference
artifacts from the source paper's public Zenodo release:

> Babbush, Zalcman, Gidney, Broughton, Khattar, Neven, Bergamaschi, Drake, Boneh,
> "Securing Elliptic Curve Cryptocurrencies against Quantum Vulnerabilities:
> Resource Estimates and Mitigations", Google Quantum AI, 2026
> (arXiv:2603.28846v2), Appendix A. Originals live in the Zenodo upload under
> `docs/example_data/`.

A small subset is vendored here (all tiny, human-readable ASCII) so that
`analysis/verify/validate_reference_adders.py` is self-contained and
reproducible without the full multi-megabyte Zenodo dump (proofs, ELF, SVGs).

- `inc3.kmx` + `inc3_test_cases.txt` — 3-qubit `+1 (mod 8)` incrementer with golden vectors.
- `iadd8.kmx`, `iadd64.kmx` — in-place `r0 += r1` adders (a variant of Cuccaro et al., quant-ph/0410184 — the same primitive this repo uses).
- `iadd8_with_ancillae.kmx` — same, using clean workspace ancilla.
- `iadd8_with_classical_offset_and_dirty_ancillae.kmx` — `r0 += r1 + carry` with a **classical** addend register and **dirty borrowed** ancilla (arXiv:2507.23079) — structurally like this repo's "quantum point += classical point" primitive.
- `inc3_wrong_{order,phase,garbage}.kmx` — the paper's **negative controls**: incorrect incrementers the validator must reject.
