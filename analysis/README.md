# analysis/ — scientific-rigor layer

Turns the challenge circuit from a leaderboard number into a verifiable,
physically-grounded result. Lives outside `editablePaths` (`src/point_add`), so
nothing here can affect the circuit or the score.

| File | What it does |
|---|---|
| `verify/solinas_reduction.py` | z3 proof: `mod_add_qq` computes `(acc+a) mod p` for **all** `acc,a ∈ [0,p)`, and its overflow ancilla uncomputes to \|0⟩. |
| `verify/peephole_identities.py` | z3 proofs of the constprop CCX identities, the ripple-carry adder recurrence, and the borrow-chain comparator (22 lemmas). |
| `verify/run_kani.sh` | Runs the Kani (bit-precise BMC) harnesses in `src/kani_proofs.rs` that bind to the **real Rust `alloy` U256 type** (not an abstract model). |
| `verify/kickmix_sim.py` | Independent, spec-faithful simulator for kickmix `.kmx` circuits (the source paper's format) — re-derives the semantics `src/sim.rs` implements. |
| `verify/validate_reference_adders.py` | Fuzz-validates the **source paper's** reference in-place adders (`verify/reference_circuits/`, from arXiv:2603.28846v2) — correct output, clean/dirty ancilla restored, phase +1 — and confirms its three negative-control circuits are **rejected**. |
| `verify/controlled_lookup.py` | Constructs and validates a self-contained **controlled** table lookup `r0 ^= ctrl ? r2[r1] : 0` (the ladder's `3·2^w` QROM primitive), in both reversible and measurement-based-uncomputation forms — the reference `table_lookup_3x3.kmx` is only an illustrative extract (issue #3). |
| `cost_model.py` | Maps the real `score.json` + `depth.json` metrics to surface-code physical resources (incl. measured runtime + spacetime volume) under explicit, editable assumptions. |
| `ecdlp_estimate.py` | Derives the **full Shor-ECDLP** cost by composing the measured per-addition primitive with the double-and-add ladder structure (`2(n+1)` additions, windowed variants); replaces the old hand-picked multiplier. Analysis-only, no `score.json` impact. |
| `../src/bin/depth_report.rs` | Standalone binary: measures toffoli-depth / gate-depth of `ops.bin` via `circuit::analyze_depth`, writes `depth.json`. Does **not** run the simulator or touch `score.json`. |
| `scientific-value.md` | Synthesis: what is proven, the cost mapping, and the generalizable vs. curve-specific techniques. |
| `adr/` | Architecture decision records for the analysis layer (isolation from scoring, derived ECDLP estimate). |

## Run everything

```bash
cargo run --release --bin depth_report   # measure depth -> depth.json (needs ops.bin)
bash analysis/run.sh                      # z3 proofs + cost model (needs python3 + z3)
bash analysis/verify/run_kani.sh          # Kani proofs on real Rust types (needs cargo kani)
```

`analysis/run.sh` requires `z3` with Python bindings (`python3 -c "import z3"`).
The Kani harnesses live behind `#[cfg(kani)]` in `src/kani_proofs.rs`, so the
normal build and `benchmark.sh` never compile them — zero effect on the score.
Every number is produced by a deterministic run; none are hand-asserted.

## Two-layer verification (why both z3 and Kani)

- **z3** (`verify/*.py`) proves the width-256 arithmetic over abstract
  bitvectors — full field coverage, fast, but a *model* of the algorithm.
- **Kani** (`src/kani_proofs.rs`) proves the exact Rust control flow of the
  Solinas reduction using the real `alloy_primitives::U256` type
  (`solinas_add_u256`, verified against the real secp256k1 prime) and a fast
  small-width twin (`solinas_add_u64`). This binds the proof to the *implementation
  types*, not just the math.
- The division-based `sub_mod` is **not** BMC-tractable (ruint's 256-bit `%` is
  Knuth long division with unbounded loops) — which is itself the argument for
  the division-free Solinas design. That path is covered by the z3 layer.
