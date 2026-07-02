# Scientific value of the ecdsafail-challenge circuit

This document turns the repository from a competitive-optimization artifact into
something with defensible scientific standing. It does three things:

1. **Formally verifies** the algebraic invariants the optimizations rely on
   (previously checked only by sampled simulation) — `analysis/verify/`.
2. **Maps the abstract score to a physical fault-tolerant cost** under stated
   assumptions — `analysis/cost_model.py`.
3. **Extracts the generalizable techniques** from the codebase and separates
   what is reusable from what is harness/curve-specific — this document.

All numbers here come from deterministic runs (`z3`, `score.json`); none are
hand-asserted. Re-run: `python3 analysis/verify/solinas_reduction.py`,
`python3 analysis/verify/peephole_identities.py`, `python3 analysis/cost_model.py`.

---

## 0. What the artifact is

A reversible circuit for **secp256k1 elliptic-curve point addition** — the inner
loop of Shor's algorithm applied to the elliptic-curve discrete-log problem
(ECDLP), i.e. the computation that breaks ECDSA (Bitcoin/Ethereum keys). It is
scored by `round(avg_toffoli_per_shot) × qubits` (`src/bin/eval_circuit.rs:434`),
where "Toffoli" counts CCX+CCZ executions (`src/sim.rs:86`) and "qubits" is the
maximum allocated qubit id + 1 (`src/circuit.rs:356`). Current metrics
(`score.json`): **1,364,230 Toffoli × 1,152 qubits = 1,571,592,960**.

This places the work in **quantum resource estimation**, a legitimate and
cryptographically policy-relevant research area. The improvement is real *if*
(a) the circuit is provably correct, and (b) the score maps to a physical cost.
Sections 1–2 supply exactly those two missing pieces.

---

## 1. Formal correctness (was: empirical only)

The harness validates correctness by *sampled simulation*: 9024 random point
pairs (`benchmark.sh`) plus `CONSTPROP_VERIFY` / `ALT_SEED_*` shot replays. That
establishes correctness on the sampled inputs, not all of them — a subtle bug on
an unsampled input would silently invalidate a "frontier-beating" claim. We
discharge the underlying claims as **theorems over all inputs** (z3 returns
`unsat` on every negation).

### 1a. Solinas modular reduction — the load-bearing arithmetic identity

`mod_add_qq` (`src/point_add/arith/modular.rs:12-49`) computes `(acc + a) mod p`
on `p = 2^256 − 2^32 − 977` using the Solinas trick: add, add `c = 2^256 − p`,
branch on the overflow bit, conditionally undo. The comment asserts this "saves
one full (n+1)-wide Cuccaro" but never proves it. `solinas_reduction.py` models
the algorithm step-for-step as 257-bit vectors and proves, **for all
acc, a ∈ [0, p)**:

```
[PROVED] mod_add_qq: low256 == (acc + a) mod p        for all acc,a in [0,p)
[PROVED] mod_add_qq: overflow flag uncomputes to |0>  (flag == (acc_final < a))
```

The second theorem is a **reversibility** guarantee: the transient overflow
ancilla returns to |0⟩, so the sub-circuit is clean and `emit_inverse`-safe —
exactly the property the challenge's ancilla-uncompute check enforces, now proven
rather than sampled.

### 1b. Peephole, adder, and comparator invariants

`peephole_identities.py` proves the boolean claims behind the gate-level
optimizations (`22/22 lemmas PROVED`):

| Claim | Source | Theorem |
|---|---|---|
| DropZeroCtrl | `constprop.rs` | `a=0 ⇒ CCX(a,b,t)=t` |
| FoldCx | `constprop.rs` | `a=1 ⇒ CCX(a,b,t)=t⊕b` |
| FoldX | `constprop.rs` | `a=1,b=1 ⇒ CCX=¬t` |
| FoldEqualCtrls | `constprop.rs` | `a=b ⇒ CCX(a,b,t)=t⊕a` |
| DropComplementCtrls | `constprop.rs` | `a=¬b ⇒ CCX(a,b,t)=t` |
| InversePairCancellation | `constprop.rs` | `CCX;CCX (controls/target unchanged) = I` |
| Ripple-carry recurrence | `venting.rs`, `arith/adder.rs` | carry chain `= (a+b) mod 2^w`, w∈{1..64} |
| Borrow-chain comparator | `comparator.rs` | final borrow `= (a <ᵤ b)`, w∈{1..64} |

The affine-form analysis in `constprop.rs` (`FoldEqualCtrls`/`DropComplement`)
proves two controls are *always* equal/opposite over GF(2); the z3 lemma
confirms the peephole is sound *given* that premise. The premise itself — that
GF(2) affine equality implies equality on every basis state — is the standard
linearity argument and is what the empirical `CONSTPROP_VERIFY` pass corroborates.

### 1c. Kani bridge — binding the proof to the real Rust types

The z3 lemmas above are a *model* of the arithmetic. `src/kani_proofs.rs`
(compiled only under `cargo kani`, behind `#[cfg(kani)]`) closes the model→code
gap with bit-precise bounded model checking on the **actual `alloy_primitives::U256`
type**:

```
cargo kani --harness solinas_add_u64    VERIFICATION: SUCCESSFUL  (0 of 3 failed,  0.33 s)
cargo kani --harness solinas_add_u256   VERIFICATION: SUCCESSFUL  (0 of 139 failed, 2.2 s)
```

- `solinas_add_u256` reproduces `mod_add_qq`'s extended-register control flow on
  real U256 values, against the **real secp256k1 prime** (`SECP256K1_P`), and
  proves it equals a division-free ground truth for all `a,b ∈ [0,p)`.
- `solinas_add_u64` is a fast small-width twin proving the control flow itself.

A useful negative result: a harness over the real `sub_mod` (which calls ruint's
256-bit `%`) does **not** converge — Knuth long division has data-dependent loops
BMC cannot unwind. Division-based modular arithmetic is not BMC-tractable, which
is precisely the argument for the division-free Solinas design; that path stays
covered by the z3 layer (§1a).

### 1d. Cross-validation against the source paper's reference circuits

The source paper (Babbush et al. 2026, arXiv:2603.28846v2) publishes reference
`iadd` circuits in the kickmix format — and `iadd8.kmx`/`iadd64.kmx` are
explicitly "a variant of the adder from quant-ph/0410184" (Cuccaro et al.), the
**same primitive** this repo's arithmetic core uses.
`analysis/verify/validate_reference_adders.py` runs them through an independent,
spec-faithful kickmix simulator (`verify/kickmix_sim.py`, a Python re-derivation
of the semantics `src/sim.rs` implements) and fuzz-checks, deterministically
(seeded), that:

- **positive controls** compute `r0 += r1` (and `r0 += r1 + carry` for the
  classical-offset variant) with the addend unchanged, all clean workspace
  ancilla returned to `|0⟩`, all **dirty borrowed** ancilla restored to their
  random input, and global phase `+1`;
- **negative controls** — the paper's `inc3_wrong_{order,phase,garbage}.kmx` —
  are **rejected** (wrong output, uncorrected phase kickback, and un-restored
  ancilla respectively).

```
POSITIVE: inc3, iadd8, iadd8_with_ancillae, iadd64,
          iadd8_with_classical_offset_and_dirty_ancillae   -> all PASS
NEGATIVE: inc3_wrong_order / _phase / _garbage              -> all REJECTED
```

The `classical_offset_and_dirty_ancillae` case matters most: it is a *classical*
addend with *dirty* borrowed ancilla and measurement-based-uncomputation phase
correction — structurally the same shape as this repo's "quantum point +=
classical point" primitive. Passing it (and rejecting the negatives) confirms the
kickmix semantics this whole repo relies on reproduce the paper's own artifacts,
and that the phase-/ancilla-aware fuzz methodology (the paper's Appendix A.5
correctness argument, mirrored by `eval_circuit`'s garbage checks) actually
catches bugs.

**Not yet covered:** the reference `table_lookup_3x3.kmx` (a measurement-based
*unary-iteration* QROM, Gidney 2018 §III.C — the `3·2^w` lookup primitive of the
windowed ladder, [ADR 0003](adr/0003-ground-ecdlp-estimate-in-source-paper.md)).
Under the current classical-trajectory simulator its decode accumulator stays
`|0⟩` (the standalone extract lacks the outer control the full circuit would
drive), so it reduces to a no-op and is **not** claimed as validated. Modelling
unary iteration's accumulator re-priming across the address walk is future work;
the adder primitives above are validated, the QROM is not.

### Scope / honesty

This verifies the **algebraic lemmas each optimization class depends on** and
binds the Solinas reduction to the real U256 type — but not a symbolic execution
of the full 28M-gate emitted circuit against the reference point-add (that does
not scale in either solver). The lemmas are the parts where bugs would hide; the
composition into a full point-add is still guarded by the sampled end-to-end
check.

---

## 2. Physical cost model (was: an abstract product)

`Toffoli × qubits` is a proxy; alone it says nothing physical. `cost_model.py`
turns the two real metrics into surface-code resources under **explicit, editable
assumptions** (physical error `1e-3`, threshold `1e-2`, `t_react = 10 µs`,
patch `= 2d²`, measurement-based Toffoli `= 4 T`). Real output for the current
circuit (one point addition):

- **Non-Clifford volume:** 5.46M T @ 4 T/Toffoli (measurement-based, the repo's
  technique) — 9.55M T @ 7 T/Toffoli (Clifford+T textbook upper bound).
- **Per-addition physical qubits:** ≈ 2.0M (d=21) to 3.4M (d=27), including a
  2× factory/routing overhead over the 1,152 logical patches.
- **Runtime:** now **measured**, not just bounded. `depth_report`
  (`src/bin/depth_report.rs` → `depth.json`) computes the non-Clifford critical
  path via `circuit::analyze_depth`: **toffoli-depth 1,077,263** (vs 1,364,230
  Toffoli gates → only **1.27× non-Clifford parallelism** — the circuit is nearly
  serial in its magic-state layer, as expected for ripple-carry modular
  arithmetic). Reaction-limited runtime = 10.77 s (vs the 13.6 s sequential
  upper bound), giving a **spacetime volume ≈ 3.6×10⁷ physical-qubit-seconds**
  at d=27.
- **This circuit vs the source paper's published bounds** (Babbush et al. 2026,
  arXiv:2603.28846v2, `docs/`). The paper zero-knowledge-proves two *point-addition*
  circuits: Low-Qubit (≤ 2.7M Toffoli, ≤ 1,175 qubits, ≤ 17M ops) and Low-Gate
  (≤ 2.1M Toffoli, ≤ 1,425 qubits, ≤ 17M ops). This repo's **measured** point
  addition — **1,364,230 Toffoli · 1,152 qubits · 10.2M ops** — is under the
  Low-Qubit bound on **all three axes**. That is the precise meaning of "beats the
  frontier": it is an improved instance of the paper's own primitive.
- **Full ECDLP extrapolation (paper's closed form, `analysis/ecdlp_estimate.py`):**
  the paper's Appendix A gives `ECDLP_Toff = (PA_Toff + 3·2^w)(2n/w − 4)` and
  `ECDLP_Qubits = PA_Qubits + w`, optimal window **w=16** → `2n/w − 4 = 28`
  windowed point additions. Substituting this repo's measured PA gives
  **(1.36M + 3·2¹⁶)·28 ≈ 43.7M Toffoli at 1,168 qubits**, reaction-limited to
  **~5 minutes** — roughly **2.06× fewer Toffoli** than the paper's published
  Low-Qubit ECDLP (≤ 90M) and 1.60× fewer than Low-Gate (≤ 70M), because the
  improved PA propagates through the ladder. (My earlier `2(n+1)=514` /
  `~7×10⁸` figure used the wrong ladder model and a `2^w` lookup; the paper's
  `28`-addition / `3·2^w`-lookup form supersedes it.)

**Key limitations this surfaces** (all real, all worth fixing):
- The scored "qubits" is `max_id + 1` (total allocated ids), **not peak
  simultaneous width** — the README's "peak qubits" is inaccurate
  (`circuit.rs:356`). A metric that rewarded true peak width would better track
  physical qubit count.
- ~~No depth / T-depth is tracked~~ **RESOLVED**: `circuit::analyze_depth` +
  `depth_report` now measure toffoli-depth and gate-depth (critical path over
  read/write hazards), feeding measured runtime and spacetime volume into the
  cost model.
- The full-attack ladder cost now uses the source paper's exact closed form
  (`(PA+3·2^w)(2n/w−4)`, w=16), but adder completeness (exceptional cases: P==Q,
  P==−Q, ∞) and the classical-vs-quantum-addend gap (this repo's PA folds a
  compile-time classical addend; the windowed ladder loads P[k] from a quantum
  table) remain assumptions pending a full-circuit build; only the per-addition
  figures are measured.

---

## 3. Generalizable techniques (the transferable science)

Catalogued from `src/point_add/`. Provenance strings in the code are real and
were verified (`venting.rs:1,311`, `mod.rs:709,21`, `gcd.rs`).

### Reusable across any modular-arithmetic quantum circuit
- **Cuccaro ripple-carry adder** (`arith/adder.rs`, `mod.rs:709`) — Cuccaro et al.
  2004 (arXiv:quant-ph/0410184). Foundation; 1 carry ancilla.
- **Measurement-based (vented) uncomputation** (`venting.rs`) — Gidney 2025
  (arXiv:2507.23079) + Häner–Roetteler–Soeken 2017 (arXiv:1709.06648). Replaces
  the ~n-Toffoli UMA uncompute with H-measure-reset + deferred conditional-CZ
  phase corrections ⇒ **zero Toffoli** in uncompute, at the cost of classical
  bookkeeping bits. This is the single largest structural saving and transfers to
  any circuit that needs to zero a carry/flag qubit.
- **2-clean-ancilla streaming adder** (`venting.rs:124`) — Gidney 2025 Fig. 2/4.
  Peak O(1) clean ancilla instead of O(n); central to the low-qubit-width score.
- **Kaliski / two-inverse conjugate uncomputation** (`mod.rs:21`, `gcd.rs`) —
  Roetteler et al. 2017 (arXiv:1706.06752) + Bernstein–Yang jump-GCD
  (arXiv:2510.10967). Field inversion reused to uncompute scratch, saving ~2×256
  ancilla qubits.
- **Sound constant-propagation peephole** (`constprop.rs`) — abstract
  interpretation over {0,1,⊥} + GF(2) affine forms drops/folds provably-constant
  CCX gates. General to any reversible circuit with initialized ancillae. Verified
  in §1b.

### Curve/harness-specific (still instructive, less portable)
- **Solinas folding** (`arith/modular.rs`, verified in §1a) — exploits the sparse
  `c = 2^32 + 977`; bespoke per Solinas prime, not general.
- **Fused double / controlled-double and symmetric square-subtract**
  (`trailmix_ludicrous/fused.rs`, `ec_add.rs`) — amortize shared folds/carries;
  depend on the a=0, b=7 group law.
- **PAD-truncated comparator recomputation** (`comparator.rs`, `arith.rs`) — trade
  a `2^-PAD` phase-miss probability for ~n→PAD recomputation width. Tunable but
  problem-specific.
- **Baked schedule / design-space search** (`trailmix_ludicrous/mod.rs`,
  `schedule.rs`, `TLM_*` env knobs) — a Pareto frontier of (carry-cap, vent-count,
  fold-width) operating points, replayed at build time. This is automated
  design-space exploration; the *method* is general, the baked tables are not.

---

## 4. Bottom line

With §1 and §2 in place, the circuit is no longer just a leaderboard number:
its arithmetic core is **proven correct over the whole field** (not just 9024
samples) — at two levels, an abstract-bitvector z3 model (§1a–b) and a
bit-precise Kani proof bound to the real `alloy` U256 type (§1c) — and its score
is **anchored to a physical cost model** with explicit assumptions and now a
**measured** toffoli-depth → runtime → spacetime volume (§2). The remaining gap
to full scientific rigor is concrete: build the full ECDLP circuit to replace the
extrapolation multiplier with a measured count. A stretch goal is symbolic
execution of the emitted op-stream on computational-basis inputs to prove the
*composed* point-add end-to-end.
