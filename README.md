# The secp256k1 Point-Addition Challenge

> **Goal.** Build the cheapest reversible quantum circuit that performs one
> elliptic-curve point addition on **secp256k1**, scored by the product of
> **Toffoli count × peak qubit width**.

---

## Why this matters

Shor's algorithm breaks elliptic-curve cryptography by computing discrete
logarithms in time polynomial in the bit-width of the curve. The quantum cost
of *running* Shor on an ECC group is dominated by one inner primitive,
repeated thousands of times: **point addition** on the curve.

Faster point addition ⇒ fewer Toffoli gates ⇒ fewer magic states ⇒ less
physical hardware and less wall-clock time on a fault-tolerant quantum
computer. Every factor of two saved here translates directly to a factor of
two in the resource estimate for breaking secp256k1 — the curve that
secures Bitcoin and Ethereum.

---

## The benchmark, precisely

You are given a Rust harness that:

1. **Builds** a reversible circuit by calling `point_add::build()`.
   The circuit must consume four 256-element registers — `target_x`
   (qubits), `target_y` (qubits), `offset_x` (classical bits),
   `offset_y` (classical bits) — and overwrite `(target_x, target_y)`
   with the affine sum `(target_x, target_y) + (offset_x, offset_y)` on
   the secp256k1 curve.
2. **Validates** the circuit by simulating it on 9024 random test points.
   Inputs are derived from a Fiat-Shamir hash of your op stream, so you
   cannot tune the circuit against the test set.
3. **Counts** every Toffoli, every Clifford, and the peak number of live
   qubits.
4. **Scores** the run as

   $$\text{score} \;=\; \overline{\text{Toffoli}} \;\times\; \text{peak qubits}$$

   where $\overline{\text{Toffoli}}$ is the average executed Toffoli count
   per shot. **Lower is better.** The score is written to `score.json`.

### What "valid" means

A run is rejected if any of the following fails:

- **Classical correctness.** All 9024 shots must produce the right
  `(R_x, R_y)`.
- **Reversibility.** Every ancilla qubit must be uncomputed to $|0\rangle$
  before being freed. `sim.rs` enforces this on every freed qubit. After
  the forward pass, every non-output qubit must again be $|0\rangle$.
- **Phase cleanliness.** The global phase across all live shots must be
  zero — no leftover phase kickback from a sloppy uncomputation.
- **Forward∘reverse identity.** Running the circuit and then its gate-
  reversed inverse must restore the original state on every qubit.

There are no loopholes. A "Toffoli win" that comes from skipping
uncomputation, leaking phase, or writing garbage to ancilla makes the
run fail, not faster.

### Reference numbers

| | Toffoli (avg/shot) | Peak qubits | Score |
|---|---|---|---|
| Current `main` | 3,942,753 | 2,715 | 1.07 × 10¹⁰ |
| Google's private low-qubit Pareto point | 2,700,000 | 1,175 | 3.2 × 10⁹ |
| Google's private low-gate Pareto point | 2,100,000 | 1,425 | 3.0 × 10⁹ |

We've run a research loop that has cut the score by ~33× from the textbook baseline.
The published Pareto frontier sits roughly **3× lower still**. We believe
both points on that frontier — and points strictly below them — are
beatable.

---

## How to play

```bash
cargo run --release -- --note "what I tried"
```

That single command builds the circuit, validates it, scores it, and
appends one row to `results.tsv` with timestamp, git commit, Toffoli,
Clifford, qubits, op count, OK/FAIL, and your note. The score is also
written to `score.json` in the format

```json
{ "score": 10704574395, "metrics": { "toffoli": 3942753, "qubits": 2715 } }
```

### What you can edit

You may modify **anything inside `src/point_add/`** — split it into
submodules, rewrite primitives, swap algorithms, refactor freely.

You may **not** touch the harness:

- `src/main.rs`, `src/circuit.rs`, `src/sim.rs`,
  `src/weierstrass_elliptic_curve.rs` — these are the contract.
- `Cargo.toml`, `Cargo.lock`, `rust-toolchain` — no new dependencies.
- `results.tsv` directly (the harness appends to it for you).


