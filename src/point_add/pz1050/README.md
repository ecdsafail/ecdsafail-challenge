# shrunken-PZ 1050-qubit secp256k1 point-add

A reproduction of **Trail of Bits TrailMix**'s `shrunken-PZ` reversible inversion
EC point-add (**1050 qubits**), translated into this challenge's op-stream format.
Verified **clean (0 / 0 / 0, all 9024 shots)** through the official `benchmark.sh`
(sandboxed `build_circuit` + trusted `eval_circuit`).

Target objective: **minimize qubit count** (Toffoli is free). This route is
structurally exact `P += Q`, so it is clean on the first try — **no island hunting**.

## Metrics (official harness)

| metric | value |
|---|---|
| qubits | **1050** |
| emitted ops | 99,669,267 |
| avg executed Toffoli | 32,252,071 |
| classical / phase / ancilla | 0 / 0 / 0 |
| shots | 9024 / 9024 OK |

(For comparison, the dialog-GCD SOTA on this repo is 1168 qubits.)

## How it works

**Fully self-contained and submission-clean** — everything lives inside
`src/point_add` (the only editable path) with **no new dependency**: `Cargo.toml`
is unchanged from the base. `ec_shrunken_pz.kmx.lz` (the 1.4 GB kmx compressed to
~38.5 MB with a custom large-window LZSS) is **embedded into the binary**
(`include_bytes!`), and `point_add::build()` (`src/point_add/mod.rs`) decodes it
by **default** with a ~40-line pure-`std` LZSS decoder (`pz1050_lz_decode`),
returning the PZ op stream instead of the dialog circuit. So a bare benchmark run
yields 1050:

```bash
ecdsafail run
# => emitted ops 99,669,267 ; qubits 1050 ; all 9024 shots OK
```

Overrides:
- `POINT_ADD_FROM_KMX=<path> ecdsafail run` — load an external `.kmx` instead.
- `POINT_ADD_DIALOG=1 ecdsafail run` — fall back to the original 1168q dialog circuit.

The two stacks share the same op IR and register layout (`reg0=tx` P.x,
`reg1=ty` P.y quantum; `reg2=ox` Q.x, `reg3=oy` Q.y classical), so translation is
a 1:1 op-name remap (via the repo's own `Circuit::from_kmx` / `Op::from_text`).

Dev tools (outside `src/point_add`, not part of a submission):
`src/bin/kmx_to_ops.rs` (`.kmx -> ops.bin`), `src/bin/pz_lz_encode.rs` (the LZSS
encoder that produced the blob), `src/bin/pz_lz_verify.rs` (bit-exact round-trip
check).

## Regenerating the kmx from TrailMix (provenance)

```bash
cd <trailmix>/trailmix
cargo build --release --bin emit_test_ec_add_shrunken_pz
CIRC_OPS_CAP=150000000 N_CASES=1 \
  ./target/release/emit_test_ec_add_shrunken_pz > /tmp/ec_shrunken_pz.kmx
# stderr: "built: 103436175 ops, peak 1050q, 32937581 tof"
```

The three qubit-saving mechanisms (pipelined Proos–Zalka divstep with no stored
transcript; Monte-Carlo variable-width GCD pack peaking at 741; whole-register
HMR-ghosting of one coordinate passenger) are TrailMix's; see their
`src/inversion/shrunken_pz_*`. Source commit: TrailMix `cd961ff`.
