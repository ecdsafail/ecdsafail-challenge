# Certified Luo et al. EEA primitive stream

This directory contains the executable fixed-schedule EEA used by
`paper2607_eea.rs`. It is derived from Luo et al.,
*Quantum Algorithm for Elliptic Curve Discrete Logarithms with
Space-Efficient Point Addition* (arXiv:2607.13816v1), but it does not rely on
the paper repository's resource-only inverse placeholder or its unsafe
fixed-step and active-window claims.

Upstream source:

- repository: `https://github.com/ZeroWang030221/Space-Efficient-Quantum-Algorithm-for-Elliptic-Curve-Discrete-Logarithms-with-Resource-Estimation`
- commit: `ac1ecffee14b5a977421b75669c52db6b4033646`
- license: MIT; retained in `UPSTREAM_LICENSE`

## Exact repairs

The emitted circuit uses:

- 1,616 microsteps, from a checked universal
  `sum(bit_length(q_i)) <= 404` continuant bound;
- a 10-bit shift pointer, required by the 592 padding rotations for `x=1`;
- a pinned 1,616-row secp256k1 active-window certificate;
- exact inactive-cell controls, endpoint decode, quotient direction, lower
  borrow, folded terminal R guard, and full length updates;
- clean-`C^3X` measurement uncomputation lowered to two executed CCX gates
  plus structural phase repair;
- exact reversed primitive replay for cleanup.

The local stream width is 581 wires: two phase bits, iteration, sign, two
259-bit work registers, two 9-bit lengths, one 10-bit shift, one 9-bit
remainder length, and 22 auxiliary wires. The terminal guard is folded into
the existing R control on the valid Algorithm-3 domain, then uncomputed before
the next phase; no retained guard lane is needed.

Pinned identities:

- active-window table SHA-256:
  `3e1961f5550249604bf044edb65f1d1bc403ed75bd7178e283685ddb4f3cb880`
- generator module SHA-256:
  `5fa0d1bdeb9e8ca76733c913e06f5ff997dd45a44267c144c51583c5ab47a560`
- stream generator SHA-256:
  `37593b625a60b7d255f39d0e704804ab264ad13634d0dd1a24185dea76106741`
- schedule certificate SHA-256:
  `5ed80df7a2a34abdf7ecc0cf2a3d0245af20fe483ea15ff6ffa53f9d466c06cf`
- aggregate manifest SHA-256:
  `bf5924fccc6236f9d50b4fdda7bd6182795e2c8a8543c6f90847ed204876c693`
- independent bit-sliced probe source SHA-256:
  `8fc19c170b59c9e376f2ddfdda04f4a800352d3851a42baa654fd5de3de57003`
- independent probe output SHA-256:
  `ba4bc85013437788aaa49bdaa5525036e2812b6ddec459db4baefb6b67cb3a18`

## Binary format

Each zstd file starts with:

```text
8 bytes  magic = P26EEA2\0
u32 LE   field width = 256
u32 LE   local width = 581
u32 LE   first schedule step (inclusive)
u32 LE   last schedule step (inclusive)
```

The payload is a stream of little-endian `u64` records. Bits `0..3`
encode the primitive kind (`1=X`, `2=CX`, `3=CCX`,
`7=clean-C^3X-MBU`), bits `4..7` encode arity, and five 10-bit local-wire
indices begin at bits 8, 18, 28, 38, and 48. The adjacent JSON file records
per-step counts and the SHA-256 of the uncompressed payload.

There are 36 chunks. The first 35 contain 45 steps each; the last contains
steps 1576 through 1616. The Rust backend checks contiguous headers and emits
each chunk independently to avoid retaining the decoded stream in memory.

For resource accounting, one kind-7 record lowers to two CCX, one HMR, and one
conditional CZ operation. Therefore:

```text
executed T per traversal = ordinary_ccx + 2 * kind7
emitted ops per traversal = records + 3 * kind7
```

The point-add integration executes four traversals: forward and reverse for
the initial quotient, then forward and reverse for exact quotient cleanup.
The verified aggregate is:

```text
records per traversal       = 150,668,315
emitted ops per traversal   = 161,442,371
executed T per traversal    = 59,599,489
four-traversal emitted ops  = 645,769,484
four-traversal executed T   = 238,397,956
```

## Regeneration

With Qiskit installed and the pinned upstream checkout available:

```sh
python generate_eea_blob.py \
  --paper /path/to/pinned-upstream \
  --out chunk-0001-0045.zst \
  --start 1 --end 45 \
  --schedule-end 1616 \
  --module eea_circuit_s835_fastdual_aux22 \
  --aux-size 22 --expected-qubits 581 --level 12
```

Repeat for the ranges encoded in `paper2607_eea.rs`. Promotion requires all
chunk hashes, the independent serialized endpoint/reverse probe, exact
count-only composition, and the official 9,024-shot benchmark.
