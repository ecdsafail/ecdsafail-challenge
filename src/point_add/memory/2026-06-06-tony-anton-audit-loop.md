# Tony + Anton Audit Loop

Status: active solver process for ECDSA.fail.

Purpose: keep optimization work from turning into blind brute force or polished-but-unsupported submission prose. The loop adapts the local Obsidian Tony/RCI pattern (`inspect -> diagnose -> cite evidence -> explain impact -> suggest smallest fix`) and the Anton positioning pattern (`claim stack -> role safety -> product/technical claim hygiene -> positioning fit -> prose/actionability`) to this benchmark.

## Frontier Snapshot

Local CLI checks on 2026-06-06 showed the promoted frontier moved to score `1,969,242,583`, from average executed Toffoli `1,504,387` and peak qubits `1,309` (`c5fada1`). The public note for that frontier described narrowing `DIALOG_GCD_APPLY_CLEAN_COMPARE_BITS` from `23` to `22` and re-hunting `DIALOG_TAIL_NONCE=251235`. Treat that as the baseline until `ecdsafail submissions --all` or `ecdsafail sync` proves otherwise.

## Required Loop

1. Sync frontier:
   - Run `ecdsafail submissions --all`.
   - Read the latest winning submission note.
   - Run `ecdsafail sync` if the promoted best moved.
2. Tony pre-change audit:
   - Problem: the exact waste, risk, or contradiction.
   - Evidence: file/function/env knob, current metric, prior note, or benchmark result.
   - Why it matters: expected Toffoli, qubit, correctness, phase, or cleanup impact.
   - Source check: compare against harness invariants and current promoted best.
   - Smallest useful fix: one bounded change only.
3. Implement the smallest useful fix.
4. Validate:
   - Full candidate: `./benchmark.sh`.
   - Fast probe: direct `build_circuit` / `eval_circuit`, with exact environment and nonce recorded.
   - Always record score, average Toffoli, peak qubits, classical mismatches, phase failures, and ancilla failures.
5. Tony post-run audit:
   - Confirm or reject the hypothesis with metrics.
   - Classify failures as structural, Fiat-Shamir/tail-search-sensitive, or measurement noise.
   - Stop brute force when failures repeat without a source-backed reason.
6. Anton submission audit:
   - Claim stack: exact change, exact score, exact validation status, exact caveat.
   - Role safety: keep ECDSA.fail, Eigen/Google, StarkWare, Starknet, and SNF roles distinct.
   - Claim hygiene: do not claim ECDSA is practically broken today or that a system is fully post-quantum safe.
   - Positioning fit: this is a quantum-circuit optimization benchmark and durability-measurement signal.
   - Prose/actionability: public note must help future solvers reproduce the result or avoid the dead end.
7. Submit only after the Anton gate passes and the audited score beats the current frontier.

## Current Tony Finding

`DIALOG_GCD_COMPARE_BITS=48` looked attractive because it reduced average executed Toffoli from `1,504,903` to `1,504,759` in local failed probes, but repeated known-clean nonce probes still produced classical mismatches and phase failures. That makes it an unproven structural or cleanup-sensitive candidate, not a tail-nonce-only win.

Smallest useful next fix: inspect the compare-screen correctness boundary and supporting cleanup assumptions before running more nonce brute force. If there is no source-backed reason why `48` can be made safe, return to the `49`-bit frontier and search a different bounded hypothesis.

## Current Validated Improvement

Tony pre-change audit selected `DIALOG_GCD_APPLY_CLEAN_COMPARE_BITS=21` because the latest shared note and local trace showed a pure `-516` average executed Toffoli cut at unchanged `1,309` peak qubits. The inherited nonce `251235` failed (`9` classical mismatches, `5` phase batches), so the local GCD pre-filter was used to hunt survivors. Candidate nonce `58422` was GCD-clean but failed full quantum validation with `1` classical mismatch. Candidate nonce `280321` was GCD-clean but failed with `2` classical mismatches and `1` phase batch. Candidate nonce `431581` validated clean over all `9,024` shots.

Validated result: `1,503,871` average executed Toffoli × `1,309` qubits = score `1,968,064,139`, with `0` classical / `0` phase / `0` ancilla failures.

## Public Note Checklist

- Include model and agent context.
- Include exact files or knobs changed.
- Include exact benchmark command and score.
- Include validation counts and caveats.
- Include one useful next lead or one dead end to avoid.
- Do not include API keys, private Obsidian prose, local-only account details, or unsupported strategic claims.
