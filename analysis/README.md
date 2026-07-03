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
| `verify/ladder_lookup_cost.py` | **Measures** the windowed-lookup (QROM) cost that was derived in the estimate (issue #4, ADR 0010). Builds an optimized **unary-iteration** table read `out ^= T[addr]` as a kickmix circuit in **both** uncompute forms — reversible and measurement-based (MBUC: `HMR` + `CZ`/`Z` phase fixup, issue #29) — validates each exhaustively (correct read, registers unchanged, ancilla cleared, phase `+1`; teeth check that deleting the MBUC fixups fails the phase test), and measures **`2^(w+1)−4` Toffoli reversible / `2^w−2` MBUC, `w` ancilla per read** — below the paper's `3·2^w`. |
| `verify/completeness_collision_rate.py` | **Measures** the affine adder's exceptional-input rate across a faithful windowed ECDLP ladder (issue #5). Exact (distribution convolution, no sampling). Validates the `2/n` equidistribution heuristic behind the completeness argument — dx=0 collisions track `2/n` within a small constant, even under large accumulator non-uniformity — and surfaces that the *dominant* exceptional term at `w=16` is the zero-window `∞` case (`~2^-11`), not dx=0: a lookup-encoding condition §4 must state (ADR 0008). Cross-checked against a real prime-order curve. |
| `verify/direct_lookup_init.py` | **Circuit-level** demonstration that the ladder's amplitude-1 `∞`-accumulator start is removed structurally (issue #5 part (a), ADR 0009). Reuses the validated controlled-lookup QROM to write `acc ^= T[w]` (`T[w] = [w]·P`) into a `\|0⟩` accumulator and shows the register holds a real affine point for every window, is the `(0,0)` `∞` sentinel **iff `w=0`**, keeps ancilla clean / phase `+1`, and stays `∞` under the `ctrl=0` negative control — so the adder is never fed `∞` at t=0. Exhaustive on a toy prime-order curve (both uncompute modes) + a secp256k1 256-bit spot-check. |
| `verify/offset_window_encoding.py` | **Removes** the dominant exceptional term of the ladder — the zero-window `∞` addend (issue #5 part (b), ADR 0015). Implements the **offset window encoding** (each digit `g → g+1`, one classical correction), proves exhaustively on a real toy curve that it never emits the `∞` table entry yet computes `[a]P+[b]Q` for every `(a,b)`, and re-runs #15's exact measurement to show the `addend=∞` rate is now **exactly 0** while `dx=0` is unchanged — sharpening the completeness headline from `~2⁻¹¹` back to the `dx=0`-limited `~2⁻²⁵⁰`. |
| `cost_model.py` | Maps the real `score.json` + `depth.json` metrics to surface-code physical resources (incl. measured runtime + spacetime volume) under explicit, editable assumptions. |
| `ecdlp_estimate.py` | Derives the **full Shor-ECDLP** cost by composing the measured per-addition primitive with the double-and-add ladder structure (`2(n+1)` additions, windowed variants); replaces the old hand-picked multiplier. Analysis-only, no `score.json` impact. |
| `../src/bin/depth_report.rs` | Standalone binary: measures toffoli-depth / gate-depth of `ops.bin` via `circuit::analyze_depth`, writes `depth.json`. Does **not** run the simulator or touch `score.json`. |
| `scientific-value.md` | Synthesis: what is proven, the cost mapping, and the generalizable vs. curve-specific techniques. |
| `completeness_argument.md` | Quantitative negligibility argument (issue #5) that the incomplete affine adder suffices for a working Shor run: exceptional-input amplitude `≈ 2⁻²⁵⁰`, >240 bits below Shor's tolerance. |
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
