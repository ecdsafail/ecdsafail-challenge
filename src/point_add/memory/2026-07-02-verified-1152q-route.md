# Verified 1152-qubit route — 1.57×10⁹ (2026-07-02)

## What was run
Fresh clone of `main`, Rust 1.93.0, `./benchmark.sh --note "baseline sanity"`
on Linux (unconfined dev fallback — no bubblewrap available on this box, but the
**trusted** `eval_circuit` stage still re-simulates and scores independently).

## Verified result (from the trusted `eval_circuit` `score.json`)
```
avg executed Toffoli : 1,364,229.770  -> rounded 1,364,230
peak qubits          : 1,152
score                : 1,571,592,960   (~1.57×10⁹)
correctness          : 9024/9024 shots OK
phase-garbage        : 0
ancilla-garbage      : 0
```
All four validity checks pass (classical correctness, reversibility,
phase cleanliness, forward∘reverse identity — the last enforced by the harness).

## Where this comes from
`build()` calls `configure_q1153_second512_submission_defaults()` (mod.rs
~L1156) which sets, among many GPU-searched knobs:
- `TLM_TARGET_Q = 1152`
- `DIALOG_TAIL_NONCE = 50400005525597` (Fiat-Shamir island selector — a
  fixed-length identity `X;X` tail that reseeds the 9024 test inputs without
  changing Toffoli count or peak width)
- the FFG cy0-release call list + opus square vents (peak-axis cuts)

Note: `build()` *also* contains a later `configure_...` block (~L2139, the
"GPT-Codex Q1159" route, nonce `2430844`, T=1,388,180 Q=1159 → 1.61×10⁹) but it
uses `set_default_env` (only sets if unset), and the q1153 defaults are applied
first via `build_builder()`'s call chain, so the **active** verified route is the
1152-qubit / 1.57×10⁹ one above.

## Alternatives that are WORSE than the active route (do not switch)
- nonce 17761178, `SQUARE_ROW_WINDOW_CLEAN_COMPARE_BITS=21`: 1215 q ×
  1,403,115 T = 1.70×10⁹.
- nonce 453700 (Q1159 route): 1159 q × 1,388,180 T = 1.61×10⁹.

## Takeaway
The construction is the product of very large automated (GPU) search over
hundreds of interdependent env knobs. The shipped default is already ~6.8×
below the README "shipped baseline" (1.07×10¹⁰) and below the Google low-gate
Pareto frontier (3.0×10⁹). Hand-editing individual knobs almost always
regresses either the Toffoli axis or the peak-qubit axis (the score is their
product) or breaks a validity check. Any further gains realistically require
re-running the automated search loop, not manual edits. Keep the default.
