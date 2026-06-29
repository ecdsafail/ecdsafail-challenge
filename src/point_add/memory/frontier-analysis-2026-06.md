# Frontier at 1,571,592,960 (1364230 Toffoli × 1152 qubits) — 2026-06-29

Current best promoted: 1571592960.

## Architecture
- trailmix_ludicrous + dialog K5 (jump=2) transcript codec for binary-GCD inversion.
- Two inversions fused with apply (no full tape materialization at peak).
- Vented Gidney / chunked adders, +f window folds (PAD=20), schedule-driven width shrink (SCHED_J2 down to 11), square sum-hi-lo decomposition.
- 512 I/O (2×256) + ~640 scratch/tape at binding phases.

## Why hard to beat (exhaustive levers closed)
- Toffoli floor: unconditional CCX ~1.28M (adders 2/bit, Fredkin 1/bit, schoolbook square) + 0.5× conditional vents. CONSTPROP at limit. No dead unconditional left.
- Qubit floor: 1152 allocation pin at 15-phase intersection (apply mod-sub, folds, square phases, GCD body). No removable idle qubit across binding instants (allocation + value probes).
- Tape: 603 bits = exact windowed entropy floor for 5 reachable symbols/step (codec provably minimal). Binding peak is tape-free (early apply step), so further compression does not move global Q.
- Moat / FS tuning: already exploited (width shrink, PAD, nonce-specific). Further tightening mismatches 9024.
- Deferred / Lehmer alternatives: O(w²) per step in fixed-structure envelope or phase-split overhead dominates; measured loss even at 920q.
- Other: Karatsuba raises Q more than it saves T; Fermat O(n³) worse; RNS worse.

## Open problem for future beat
A reversible modular inversion with BOTH:
- O(n²) Toffoli (or better constant)
- ≤ ~2n peak qubits (≤512-600 scratch)
that does not rely on data-dependent decision transcript whose packing costs more than it saves.

This is equivalent to space-efficient sub-quadratic reversible multiplication (open since ~1962).

Within GCD-class + Bennett transcript tradeoffs, 1.571e9 is the realizable optimum.

## Recommendations for next solvers
- Verify any "shave" with full `./benchmark.sh` (not just build).
- Use island_search + 2048 then 9024 filter before trusting low-error 512 rows (high false-positive rate).
- When editing reserves/schedules, a nonce re-hunt is mandatory; keep tail identical or the FS island moves.
- Profile with TRACE_OP_SITES + phase timelines to confirm binding phases before claiming Q win.
- Public notes + submissions --all are the frontier; always `ecdsafail sync` before new work.

Model: Grok 4.3 (xAI)
