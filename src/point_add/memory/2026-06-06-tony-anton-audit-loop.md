# Tony + Anton Audit Loop

Status: active solver process for ECDSA.fail.

Purpose: keep optimization work from turning into blind brute force or polished-but-unsupported submission prose. The loop adapts the local Obsidian Tony/RCI pattern (`inspect -> diagnose -> cite evidence -> explain impact -> suggest smallest fix`) and the Anton positioning pattern (`claim stack -> role safety -> product/technical claim hygiene -> positioning fit -> prose/actionability`) to this benchmark.

## Frontier Snapshot

Local CLI checks on 2026-06-06 showed submission `a66b042` promoted as the current frontier with score `1,967,891,695`, from average executed Toffoli `1,503,355` and peak qubits `1,309`. The accepted route narrows `DIALOG_GCD_APPLY_CLEAN_COMPARE_BITS` to `20` and uses `DIALOG_TAIL_NONCE=721381`. Treat that as the baseline until `ecdsafail submissions --all` or `ecdsafail sync` proves otherwise.

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

Validated result: `1,503,871` average executed Toffoli × `1,309` qubits = score `1,968,567,139`, with `0` classical / `0` phase / `0` ancilla failures.

Submission `436b516` promoted at 2026-06-06 08:44 local time. It beat the previous observed promoted frontier `83e3b66` (`1,968,793,475`) by `226,336` score points.

Public correction note `5ec74c1` records that the original submission prose had arithmetic typos in the displayed score and frontier delta; the CLI claimed score, metrics, validation result, and promoted leaderboard result were correct.

## Promoted Successor Frontier

External submission `a66b042` by `jackylee0424` promoted after `436b516`. Public note: apply-clean comparator tightened to `20` with refreshed tail nonce `721381`, validated `0` classical / `0` phase / `0` ancilla over all `9,024` shots at `1,309` qubits × `1,503,355` average Toffoli = score `1,967,891,695`. Local `./benchmark.sh --note 'validate synced a66 frontier'` reproduced the same `0/0/0` result.

## Current Search Audit

2026-06-06 continuation tested three bounded follow-up hypotheses. None produced a submit-ready improvement yet:

- `1285q` restack from submission `83e3b66` plus `DIALOG_GCD_APPLY_CLEAN_COMPARE_BITS=21`: structural probe was `1,531,619` average Toffoli × `1,285` qubits = score `1,968,130,415`, which would beat `436b516` by `436,724` if clean. Nonce `0` failed full eval with `18` classical mismatches and `7` phase batches. Staged GCD search found `320` candidates that passed `2,048` shots, but `0` passed the full `9,024`-shot GCD filter.
- `DIALOG_GCD_APPLY_CLEAN_COMPARE_BITS=20` with `DIALOG_GCD_COMPARE_BITS=50`: structural probe was `1,503,463` average Toffoli × `1,309` qubits = score `1,968,033,067`, which would beat `436b516` by `534,072` if clean. Nonce `0` failed full eval with `11` classical mismatches and `6` phase batches. Staged GCD search found `278` candidates that passed `2,048` shots, but `0` passed the full `9,024`-shot GCD filter.
- `KAL_FOLD_CARRY_TRUNC_W=20`: structural probe was `1,503,355` average Toffoli × `1,309` qubits = score `1,967,891,695`, which would beat `436b516` by `675,444` if clean. Inherited nonce failed full eval with `15` classical mismatches and `5` phase batches. Staged GCD search found `260` candidates that passed `2,048` shots, but `0` passed the full `9,024`-shot GCD filter.

Tony post-run classification: these remain structurally attractive but nonce-island-limited. The sampled failures are full-GCD width/nonconvergence rejections, not branch-comparator mismatches. A next pass should either cover much more nonce space with a faster full-shot filter or find a source-backed way to reduce width/nonconvergence pressure without crossing a qubit break-even cliff.

Public standalone note `6b2eea8f` shares this negative evidence for the collaborative solver pool.

## Second Continuation Audit

2026-06-06 later continuation added four more bounded checks:

- Fast filter tooling: built `/tmp/ecdsafail-fast-filter` using native `secp256k1` generator multiplication and point addition. Cross-check against the original `k256` filter on the same `1285q`/`COMPARE_BITS=48` route and nonce range produced identical `512`-shot hit lists, so the faster tool is acceptable for search triage.
- `KAL_FOLD_CARRY_TRUNC_W=20` with `DIALOG_GCD_WIDTH_SLOPE_X1000=1013`: structural probe was `1,503,835` average Toffoli × `1,309` qubits = score `1,968,520,015`, barely under the frontier by `47,124`. Nonce `0` failed full eval with `275` classical mismatches and `79` phase batches. The `512`-shot scout found `32` early candidates, but all failed at `2,048` shots; `1012` and `1011` were structurally too expensive.
- Current `1309q` route with `DIALOG_GCD_COMPARE_BITS=48`: structural probe was `1,503,727` average Toffoli × `1,309` qubits = score `1,968,378,643`. Nonce `0` failed full eval with `17` classical mismatches and `10` phase batches. Staged search found `300` candidates that passed `2,048` shots, but `0` passed the full `9,024`-shot GCD filter; known clean nonces from adjacent routes did not transfer.
- `1285q` restack from `83e3b66` with `DIALOG_GCD_COMPARE_BITS=48`: structural probe was `1,531,883` average Toffoli × `1,285` qubits = score `1,968,469,655`. Nonce `0` failed full eval with `9` classical mismatches and `3` phase batches. Staged search found `290` candidates that passed `2,048` shots, but `0` passed the full `9,024`-shot GCD filter. Known clean nonces did not transfer; a short direct full-shot search with the fast filter over spaced ranges found no clean nonce before being stopped.
- `DIALOG_GCD_APPLY_FINAL_WINDOWED_FAST_BLOCKS=3/4`: exact but structurally worse. Blocks `3` gave `1,520,899` average Toffoli; blocks `4` gave `1,537,927`, both at `1,309` qubits, so this is not a viable near-frontier path.

Tony post-run classification: every attractive near-frontier route is still bottlenecked by full-shot GCD width/nonconvergence, with comparator mismatches not showing up in the sampled filter rejects. The next useful work is either a genuinely faster full-shot nonce search or a structural qubit-floor change; small Toffoli cuts are now mostly island-limited.

## Post-a66 Search Audit

After syncing to `a66b042`, three immediate one-bit successors were probed from the new frontier. None produced a submit-ready improvement yet:

- `DIALOG_GCD_APPLY_CLEAN_COMPARE_BITS=19`: structural probe was `1,502,839` average Toffoli × `1,309` qubits = score `1,967,216,251`, which would beat `a66b042` by `675,444` if clean. Nonce `0` failed full eval with `17` classical mismatches and `14` phase batches. Staged search found `305` candidates that passed `2,048` shots, but `0` passed the full `9,024`-shot GCD filter.
- `KAL_FOLD_CARRY_TRUNC_W=20`: structural probe was also `1,502,839` average Toffoli × `1,309` qubits = score `1,967,216,251`. Nonce `0` failed full eval with `18` classical mismatches and `10` phase batches. Staged search found `314` candidates that passed `2,048` shots, but `0` passed the full `9,024`-shot GCD filter.
- `DIALOG_GCD_COMPARE_BITS=48`: structural probe was `1,503,211` average Toffoli × `1,309` qubits = score `1,967,703,199`, which would beat `a66b042` by `188,496` if clean. Nonce `0` failed full eval with `11` classical mismatches and `7` phase batches. Staged search found `323` candidates that passed `2,048` shots, but `0` passed the full `9,024`-shot GCD filter.

Tony post-run classification: after `a66b042`, the next one-bit cuts are again structurally attractive but full-shot GCD-island-limited. The sampled rejects remain width/nonconvergence, with no comparator mismatches in these filter passes.

## Post-a66 Extended Audit

2026-06-06 follow-up checked whether the old `1285q` qubit-floor route or a shorter active-iteration schedule could pair with the `a66b042` apply-clean frontier. Neither produced a submit-ready improvement:

- `1285q` restack from `83e3b66` plus `DIALOG_GCD_APPLY_CLEAN_COMPARE_BITS=20`: structural probe with `DIALOG_GCD_COMPARE_BITS=50` was `1,531,103` average Toffoli × `1,285` qubits = score `1,967,467,355`, beating `a66b042` by `424,340` if clean. Nonce `0` failed full eval with `13` classical mismatches and `5` phase batches. Staged GCD search found `285` candidates that passed `2,048` shots, but `0` passed the full `9,024`-shot filter; known adjacent clean nonces did not transfer.
- `1285q` restack plus apply-clean `20` and `DIALOG_GCD_COMPARE_BITS=48`: structural probe was `1,530,851` average Toffoli × `1,285` qubits = score `1,967,143,535`, beating `a66b042` by `748,160` if clean. Nonce `0` failed full eval with `12` classical mismatches and `6` phase batches. Staged GCD search again found `285` candidates passing `2,048` shots, but `0` passed all `9,024` shots; known nonces still did not transfer.
- Split `1285q` levers were not independently viable: `shiftOnly` gave `1,507,999` average Toffoli × `1,308` qubits = score `1,972,462,692`, `suffixOnly` kept `1,309` qubits with `1,526,063` average Toffoli, and disabling both returned to the current `1309q` control. This suggests the old `1285q` win needs the coupled restack, not a single transplantable lever.
- `DIALOG_GCD_ACTIVE_ITERATIONS=257`: structural probe was `1,500,368` average Toffoli × `1,309` qubits = score `1,963,981,712`, beating `a66b042` by about `3.89M` if clean. Nonce `0` failed full eval with `13` classical mismatches and `9` phase batches. Staged GCD search found `132` candidates passing `2,048` shots, but `0` passed all `9,024` shots; full rejects were dominated by nonconvergence (`102`) plus width (`30`).

Tony post-run classification: `ACTIVE_ITERATIONS=257` is the largest structural prize but appears to create a nonconvergence floor, while the `1285q` + apply-clean route remains width/nonconvergence island-limited. Future work should prioritize source-backed convergence or width relief before wider blind nonce sweeps.

## Current Validated Successor

Tony pre-change audit found an exact slack-spend route: `DIALOG_GCD_APPLY_FINAL_LOWQ=0` with `DIALOG_GCD_APPLY_FINAL_WINDOWED_FAST_BLOCKS=0` removes the final apply chunk's low-q/windowed carry overhead while the global peak remains bound by `round84_fused_square_xtail_dx_sub_lam_square_lowq` at `1,309` qubits. The raw fast-final route at active `258` had structural target `1,486,327` average Toffoli × `1,309` qubits = score `1,945,602,043`, but the inherited nonce failed with `19` classical mismatches and `8` phase batches, and the first `500` two-thousand-shot GCD survivors produced no full `9,024`-shot GCD hit.

Smallest useful fix: spend part of that recovered Toffoli budget on convergence by setting `DIALOG_GCD_ACTIVE_ITERATIONS=262`, keeping `DIALOG_GCD_WIDTH_MARGIN=10` and `DIALOG_GCD_WIDTH_SLOPE_X1000=1014`. Active `262` stays at `1,309` peak qubits and structural target `1,497,795` average Toffoli. Current nonce `721381` still failed (`7` classical mismatches and `4` phase batches), but the GCD prefilter became much denser:

- `500` candidates passed the `2,048`-shot filter by nonce `2620`.
- `93` of those passed `4,096` shots.
- `6` passed `8,192` shots: `614`, `1328`, `1718`, `2148`, `2432`, `2499`.
- `4` passed all `9,024` GCD shots: `1328`, `2148`, `2432`, `2499`.

Quantum confirmation results:

- `1328`: GCD-clean but failed with `1` phase-garbage batch.
- `2148`: GCD-clean but failed with `1` classical mismatch and `2` phase-garbage batches.
- `2432`: validated clean over all `9,024` shots with `0` classical / `0` phase / `0` ancilla failures.
- `2499`: GCD-clean but failed with `1` classical mismatch and `2` phase-garbage batches.

Validated result: `1,497,795` average executed Toffoli × `1,309` qubits = score `1,960,613,655`, beating `a66b042` by `7,278,040` score points. Local official path `./benchmark.sh --note 'validate lowq0 active262 nonce2432'` reproduced the clean result and wrote `score.json` with the same score.

## Public Note Checklist

- Include model and agent context.
- Include exact files or knobs changed.
- Include exact benchmark command and score.
- Include validation counts and caveats.
- Include one useful next lead or one dead end to avoid.
- Do not include API keys, private Obsidian prose, local-only account details, or unsupported strategic claims.
