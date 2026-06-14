# q1192 WMI CUDA search

## Candidate

The current 1192-qubit candidate uses:

```text
DIALOG_GCD_FOLD_CARRY_TRUNC_W=18
DIALOG_GCD_FOLD_PARK_LOW_CARRIES=13
DIALOG_GCD_FOLD_HOST_N10=1
SQUARE_ROW_MAX_SEG=165
KAL_FOLD_CARRY_TRUNC_W=20
DIALOG_GCD_SPECIAL_FOLD_RELEASE_SCRATCH=1
DIALOG_GCD_SPECIAL_FOLD_PARK_LOW_CARRIES=1
```

The trusted nonce-0 run measured average executed Toffoli `1,419,907.236`,
zero ancilla-garbage batches, 20 classical mismatches, and 15 phase-garbage
batches. A clean nonce at the rounded Toffoli count would score
`1,692,529,144`.

## Exact GPU filter

`tools/cuda/island.cu` was ported from the earlier CUDA searcher and corrected
to match the Rust and Metal models:

- MSB suffix comparisons for apply cleanup;
- per-step overflow and underflow cleanup widths from the state trailer;
- separate square and apply phase-risk counts;
- exact shot windows;
- zero-phase early rejection for production search.

WMI parity job `57587` passed the serialized SHAKE probe and matched all 65
Metal survivors over nonces `[0,100)` and shots `[0,256)`, including phase
breakdown.

## Throughput snapshot

With full 9024-shot zero-phase early rejection and comb-20:

| GPU | Nonces/s |
|---|---:|
| RTX 4090 | 3,664 |
| RTX 3090 | about 3,000 |
| A100 80 GB | 1,800 |
| L40S | 4,454 |

Nonce fan-out did not improve throughput. Comb-20 was only about 1.3% faster
than comb-8, showing that the remaining bottleneck is field and transcript
arithmetic rather than table lookup.

## Campaign

WMI array `57599` searches disjoint one-million-nonce shards starting at
nonce `100,000`, with completion markers under `checkpoints/`. The initial
range is 100 million nonces. A `CLEAN` result is only a filter survivor and
must pass the trusted local evaluator before any submission decision.

## Verification audit

The repository history contains 319 server-accepted snapshots, including the
current `833642f` record at 1,203 qubits. They are verified historical
fallbacks, but they are already submitted and are not new candidates.

The q1192, q1191, q1189, q1188, and q1187 routes must not be described as
submission-ready until an exact nonce passes the trusted evaluator over all
9,024 shots with zero classical, phase, and ancilla failures. Exact arithmetic
tests and cross-backend filter parity are necessary but not sufficient.

## q1187 route

The stream-carry31 host-d route reaches 1,187 qubits with trusted nonce-0
calibration `1,429,540.083` average Toffoli. A clean run would score about
`1,696,863,980`, improving on `833642f` by about `534,133`.

Evidence completed:

- exact freed-tail self-test over add/subtract, full/windowed tails, all
  `(e,d)` combinations, and all 64 packed lanes;
- serialized state SHA-256
  `ae5cf33c53ef72480fc1834cbd61b7bea8d8f022a81273e185e232c0b10a33bd`;
- Rust/Metal/CUDA parity over nonces `[0,100)` and shots `[0,256)`;
- local full-shot search over 50,000 nonces.

No trusted full-shot clean nonce has been found. WMI array `57682` searches
disjoint one-million-nonce shards with two concurrent GPU tasks. Its first two
million nonces completed with zero filter-clean results.

## q1188 parity and scheduling

WMI CUDA parity job `57635` passed all 67 Rust/Metal reference rows for the
q1188 state, including the serialized SHAKE probe and phase-risk counts. Search
array `57636` then started with two concurrent one-million-nonce shards.

The older q1192 array `57599` and q1189 array `57627` were released after q1188
parity completed. q1191 array `57620` remains held because it is superseded.

## Later verification and lower-qubit routes

The q1192 search later produced three nonces that passed the independent Rust
full-shot audit and trusted 9,024-shot evaluator: `36,909,818`, `49,017,993`,
and `77,101,583`. They are three submission-ready artifacts for one distinct
q1192 circuit configuration. Nonce `49,017,993` is strongest at score
`1,692,524,376`. The reproducible package is under
`optimizer/verified/q1192/`.

Four newer routes reached q1190, q1189, q1188, and q1187. Their exact
self-tests, profiles, serialized states, and CUDA parity jobs all pass. Full
9,024-shot WMI searches run as jobs `58213`, `58222`, `58214`, and `58219`,
respectively. None is submission-ready until a filter-clean nonce also passes
the independent Rust audit and trusted evaluator.

## q1186 balanced apply schedule

The accepted frontier moved to submission `ad4cf86` at q1193, average Toffoli
`1,412,391`, and score `1,684,982,463`.

A nonuniform 16-block apply schedule
`16,32,48,65,82,99,116,133,150,167,183,199,214,229,243`, combined with park
19 and temporary release of the clean K5 transcript block during apply shifts,
profiles at q1186 with emitted Toffoli `1,509,838`. Fifteen blocks cannot fit:
the required final block reaches q1210.

The q1186 route passed the fold, special-fold, square-window, and fused-apply
differential self-tests. Its serialized state SHA-256 is
`52b12fa20fdaf8c4a999497be4e4b94b386333dc6dd0053f4388be959754d70e`,
and WMI job `58237` passed Rust/CUDA parity for all 59 stage-256 survivors.
Full 9,024-shot search array `58246` is queued.

Trusted nonce-0 calibration measured average Toffoli `1,430,340.535`, 12
classical mismatches, 12 phase-garbage batches, and zero ancilla-garbage
batches. Even a clean nonce at the rounded cost would score `1,696,384,426`,
so the route still needs about 9,614 fewer average Toffolis to beat the current
record.
