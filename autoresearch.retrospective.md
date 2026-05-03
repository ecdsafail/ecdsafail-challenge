# Autoresearch Retrospective and Future Filter

_Last updated: 2026-04-29 after BY selector/plumbing, plus-minus denominator-normalization failures, and work-style retrospective._

## Why this file exists

I let two attractive routes run too long before stopping to ask whether the remaining uncharged pieces could still fit the global budget:

1. **BY / safegcd replay**: the replay body and local window oracles looked SOTA-shaped, but denominator selection/plumbing/cleanup consumed the margin.
2. **Plus-minus scaled DIV**: the step/history/slack shell became genuinely circuit-validated, but denominator shifting/normalization destroyed the Toffoli budget.

Going forward, a route is not allowed to stay "active" merely because some subcomponent is elegant. It must continue to pass a full-system budget gate with selector/parser/history/cleanup costs either measured or bounded by an explicit remaining budget.

## Process failure modes and work-style changes

### Failure modes I need to avoid

1. **Component-success anchoring.** I over-weighted real local wins: BY replay/window oracles and plus-minus step/history/slack circuits. Those wins were genuine, but they did not imply the full architecture could afford the still-missing denominator selector/parser/normalizer.
2. **Deferred accounting for the hardest piece.** I repeatedly let a route proceed with phrases like "selector later", "parser later", "normalization later", or "scale correction later". In both BY and plus-minus, that deferred piece turned out to be the route-killer.
3. **Proxy/model creep.** A cost model that was initially labelled optimistic kept informing later work after new evidence had invalidated its assumptions. In particular, plus-minus looked excellent while denominator shifts were treated as cheap/relabel-like; once physical denominator motion was charged, the route changed category.
4. **No explicit stop threshold before experiments.** I often measured a promising piece first and only afterwards compared it to the global budget. That makes it too easy to continue because the latest sub-result is good, even if the remaining budget is already impossible.
5. **Sunk-cost integration bias.** After investing in wiring, I tended to ask "what can repair this?" instead of first asking "does even the most generous repair still fit?" This caused extra BY plumbing and plus-minus normalization exploration after the critical cost class was already suspect.
6. **Insufficient separation between active, speculative, and archival ideas.** The ideas file accumulated many true statements, but did not always make the current decision state obvious. That made it easier to accidentally revive a killed route without a new premise.
7. **Cleanup/control underestimation.** I treated exact reversible cleanup, phase cleanliness, direction recovery, top-bit predicates, and boundary decoders as secondary until they became blockers. For this problem they are often the main cost, not bookkeeping.
8. **Not surfacing confidence changes early enough.** I kept working autonomously after my confidence should have dropped from "promising route" to "only a structural breakthrough can save this". That state transition should be written down immediately.

### Concrete changes from now on

1. **Pre-mortem before deep work.** Before spending more than two experiments on a route, write a short go/no-go note with: live state, scratch, all-in Toffoli budget, missing hard pieces, and exact kill conditions.
2. **Hardest-piece-first ordering.** Do not optimize or integrate the easy body while the selector/parser/normalizer/cleanup is uncharged. First build or bound the piece most likely to kill the route.
3. **Lower-bound-first accounting.** For every proposed repair, first compute a generous lower bound. If the lower bound misses the global budget before charging controls/cleanup, stop immediately and demote the route.
4. **Budget ledger update after every kept route experiment.** A `keep` that affects architecture viability must update either `autoresearch.retrospective.md`, `autoresearch.ideas.md`, or `scratch600_frontier.rs` with the remaining gap and the current uncharged blockers. No route may remain implicitly active through old optimistic text.
5. **Two-experiment checkpoint.** After at most two exploratory experiments on the same route, pause for a written decision: continue, demote, or pivot. Continuing requires a numeric budget path, not just a new subcomponent idea.
6. **Promotion requires charged hard pieces.** A candidate can be called "active" only if selector/parser/history/normalization/cleanup are measured or have an explicit affordable budget. Otherwise it is "speculative" or "archival" even if its state shape fits.
7. **Adversarial-accountant checklist.** Before integration, ask: where is the hidden history? who provides the branch/control bit? how is it uncomputed? is phase clean? does the live state fit 600--663 scratch? what is the worst-case rather than p99 cost? would this still pass if all optimistic oracles were charged?
8. **Sunk-cost interrupt.** If a new result shows the route misses by more than the remaining plausible savings, stop coding that route. The next action must be documentation/demotion unless the premise changes.
9. **Backlog hygiene.** Move failed subpaths into the explicit stop list with the premise that killed them. Revive only by naming the new premise that invalidates the old kill reason.
10. **User-visible confidence changes.** When confidence drops materially, say so in the session summary/documentation instead of continuing the loop as if the route remains primary.

## Hard budget gate

Current exact default remains clean at roughly:

```text
avg_toffoli = 4,082,910
qubits      = 2,713
scratch beyond tx,ty = 2,201
```

The Google low-qubit target is approximately:

```text
total target = 2,700,000 Toffoli @ about 1,175q
tx,ty data   = 512q
scratch target beyond tx,ty ~= 600--663q
```

For a low-scratch DIV replacement architecture, the current non-DIV scaffold floor used in the scratch-600 models is about:

```text
scaffold / non-DIV budget ~= 642,716 CCX
remaining for two DIVs + scale + selector/parser/cleanup <= 2,057,284 CCX
```

If the route still has about 404 total update steps across two DIVs, the absolute per-update ceiling is:

```text
2,057,284 / 404 ~= 5,092 CCX/update
```

That ceiling assumes **zero** cost for scale correction, parser cleanup, normalization predicates, and production packing. A believable route should therefore aim for:

```text
<= 4,000--4,600 CCX/update all-in
```

Any route whose optimistic lower bound is already above the ceiling before charging a known-hard parser/selector/cleanup piece must be demoted immediately.

## What BY taught us

BY was attractive because several pieces were real and phase-clean:

- denominator branch history can self-clean in reverse,
- lowword pattern/q oracles are cheap locally,
- selected/window interfaces compose in the real affine path,
- centered signed product-clean replay can be made phase-clean.

But the full route failed because the hard missing piece was not the replay body; it was denominator generation/selection/plumbing:

```text
best fully charged scratch-600 BY row ~= 2,765,676 CCX
remaining gap to 2.7M ~= +65,676 CCX
```

More detailed blockers:

```text
w=4 selector-only projection ~= 2,685,036, but uncharged plumbing kills it
naive full-pair plumbing excess ~= 903,996 CCX
fixed-matrix plumbing excess ~= 306,178 CCX
full-ratio A-inverse projection ~= 9,952,686 CCX
consumed-denominator reverse ambiguity: w4=48 patterns, w16=589,824 patterns
last fixed-window denominator object ~= 20,323 CCX/window mean, about 2x 2.7M target
relaxed 3M fixed-window budget ~= 13,431 CCX/window; free-m/q body still +338,944, last-shot body +496,258
```

**BY is not dead in principle**, but it should only be revived for one specific breakthrough:

```text
a selected/window-local denominator primitive below roughly 10k CCX/window,
with explicit reversible cleanup and no hidden field-sized selector/parser state.
```

2026-04-29 update: a partial-prefix qoffset-mask replay primitive is the first post-retrospective result that reopened a narrow low-scratch BY **one-DIV** budget in a hardest-piece-first way rather than by integration optimism. With 36 lowword windows and 564 harness-scale steps, prefix32/48/64/80/90 scratch is `542/558/574/590/600` and one-DIV projected gaps are `-2,476/-47,596/-92,716/-137,836/-166,036` after charging lowword selector and decoder. Broader local validation passed at n=`8,10,12,16` with phase/dirty restoration.

Important correction after user skepticism: the adversarial two-denominator ledger blocks naive BY promotion. If pair1 tagged-DIV and pair2 product-clean replay each need their own 564-step replay+selector/decode, the total is `4,068,262` (`+1,368,262` over 2.7M). This matches the previous BY blow-up failure mode: a good local replay primitive is not a full point-add architecture. Partial-prefix qoffset is therefore only a useful local primitive unless a separate charged algebra deletes the second denominator/replay.

Do not continue BY integration/plumbing experiments unless that primitive (or a successor) survives those charged hard-piece gates first, including the two-denominator objection. For the relaxed 3M question, the same accountant rule still matters: Strategy E deletes the second denominator algebraically only if its product-clean multiply is a new non-DIV primitive. With the current product-clean replay, the single-DIV side can afford only `911,490` CCX; granting the best fixed-control replay `873,600` leaves just `37,890` for selector/parser/cleanup, while measured decoder alone is `63,936` and lowword selector+decoder is `278,208`. Worse, the known product-clean multiply is itself denominator-controlled: forgetting its second selector makes an optimistic centered ledger look `31,842` under 3M, but charging that second selector/parser adds `278,208`, leaving `+246,366` even with centered product-clean replay and `+518,526` with the current product-clean replay. A direct secp branch-sharing probe found the BY branch streams for `dx` and `Rx-Qx` essentially independent (`odd_mi≈2.44` millibits, `case_mi≈4.85` millibits), so control reuse has no simple empirical support. A follow-up endomorphism denominator-sharing probe on `Py+Qy` vs `Ry-Qy` is also independent (`odd_mi≈3.09` millibits, `case_mi≈5.52` millibits). Thus a <3M low-scratch BY/Strategy-E path still needs a real non-DIV product-clean multiply or a much deeper control-sharing invariant, not just a replay-body number.

## What plus-minus taught us

Plus-minus was attractive for a different reason: it solved the state-shape problem better than most candidates.

Validated pieces now include:

- productive in-place step and inverse/roundtrip at toy widths,
- multi-step composition,
- local direction recovery from coefficient divisibility, so no persistent direction flags,
- active-chain unary history,
- high-bit slack slots used as history storage,
- active-aware terminal no-op fixed loops,
- fixed-bound packed active loop,
- unsigned/signed barrel shift primitives,
- Clifford-only unary-history to binary-k extraction.

The scratch/history model is still the best evidence that this family can fit the Google scratch regime:

```text
scaled plus-minus slack scratch ~= 517 bits in the model
```

But the denominator arithmetic killed the Toffoli path:

```text
repeated physical shifts: W^2, extrapolated 257-bit forward step ~= 150,578 CCX
per-step barrel shifts after Clifford k extraction: ~= 10,243 CCX/update
barrel two-DIV step body ~= 4,138,172 CCX
coefficient offsets + denominator barrels ~= 3,171,400 CCX two-DIV step body
```

Denominator offsets also failed as a simple escape:

```text
denominator offset raw width p99/max = 382/395 bits
periodic normalization p99 count = 89 per DIV
simple public normalization conflicts by step 2
exponent-only normalization mismatch rate = 10,037 ppm
```

The latest generous lower-bound model gives plus-minus a magic exact denominator-normalization oracle and still misses:

```text
base update after coefficient offsets          = 5,794 CCX/update
optimistic p99 denominator normalization cost = 89 * 1,285 CCX per DIV
one DIV step+normalization                    = 1,284,753 CCX
two DIVs                                     = 2,569,506 CCX
total before scale/oracle cleanup             = 3,212,222 CCX
gap before scale to 2.7M                      = +512,222 CCX
```

2026-05-03 update: the offset-normalization path is even worse than that
headline miss. With the current scaffold, the two-DIV budget is `2,057,284`
CCX, so each DIV gets only `1,028,642` CCX. The non-denominator update core
alone is:

```text
202 steps * 5,794 CCX/update = 1,170,388 CCX per DIV
base-core excess             =   141,746 CCX per DIV
required update ceiling      =     5,092 CCX/update
required base-step cut       =       702 CCX/update (~12.1%)
```

So denominator-normalization scheduling cannot rescue this subpath by itself.
It would first need a separate double-digit-percent reduction in the already
optimistic non-denominator update core, before paying any normalization oracle,
scale correction, parser cleanup, or production packing.

Therefore the current plus-minus physical-shift / barrel-shift / offset-normalization subpath is **gate-dead**. Plus-minus should only be revived if a new denominator recurrence eliminates physical denominator shifting/normalization, rather than optimizing the current normalization machinery.

## Introspection cadence for future work

Every new route gets at most **two exploratory experiments** before a go/no-go note is written. The note must answer:

1. **State gate**: what is the persistent live state, and does it fit <=600--663 scratch beyond tx,ty?
2. **Global Toffoli gate**: after adding known scaffold, what is the all-in target budget for the missing piece?
3. **Hard-piece accounting**: selector/parser/history/normalization/cleanup costs must be charged or assigned a maximum affordable budget.
4. **Lower-bound kill test**: if an optimistic lower bound already misses, stop.
5. **Circuit reality check**: before integration, validate a toy reversible circuit for the nontrivial control/cleanup mechanism.
6. **Promotion rule**: no route may be called SOTA-shaped if the only reason it fits is an uncharged parser/selector/normalizer.

## Future approaches with a real chance

See `autoresearch.literature.md` for the 2026-04-29 online sweep. Public low-qubit ECDLP papers currently found either withhold the relevant point-add circuit (Google) or buy qubits with enormous Toffoli counts (Luo/PZ-style register-sharing EEA, CFS/RNS-style low-space lines). So the future focus remains custom structural primitives rather than importing a public inversion circuit.

Ranked by current plausibility:

### 0. Centered / ordinary Euclid quotient stream for the relaxed 3M/current-qubit target

This is **not** a Google-low-qubit candidate, but it is the first post-BY result that looked numerically relevant to the user's relaxed “3M while under qubit budget” question. The old quotient-stream route was killed by the ~600-scratch parser requirement; if the cap is the current project cap (`<=2800q`), explicit quotient boundaries may fit.

Ordinary Euclid lower-bound ledger: payload p99/max `349/355` bits, count p99 `173`, one-boundary scratch p99 `777`, conservative peak with 512q workspace `1801q`; with per-qbit coefficient replay `587` CCX and long-division trial unit `8` CCX, one DIV projects `932,047` and point-add projects `2,506,810` (`-493,190` to 3M). Immediate adversarial correction: this relies on a dynamic/packed extractor. A fixed reversible scan over all 256 shifts per quotient has p99 static bit-trials `11,337,728` (`249.5×` weighted), gap `+43,403,354` even at `1` CCX/bit-trial, and a unit budget of only `0.043` CCX. A packed quotient-bit extractor has a narrow target: one-way extraction budget `486,889`, compare/sub floor `268,032`, leading scans `44,288`, leaving `~1,009` CCX per quotient for shifted-divisor alignment; a generic log barrel would miss by `+718,940` point-add CCX.

Centered Euclid improves the relaxed ledger: payload p99/max `336/341`, count p99 `118`, one-boundary scratch p99 `710`, weighted extraction p99 `43,935`, projected point-add `2,443,100` (`-556,900` to 3M). Fixed scan is still dead (`+28,970,172` at 1 CCX/static trial). A first packed-extractor note overestimated alignment room by forgetting the forward+reverse denominator pass; corrected one-way extraction budget is `490,705`, leaving `~1,716` CCX/quotient for alignment after compare/sub and leading-scan floors. A generic `n log n` barrel at `2048`/quotient would miss by `+156,860` point-add CCX under the 3n compare+masked-sub floor. Fixed-K public-shift slots also fail: `K=4` is barely under budget but fails all samples, `K=5` misses by `+357,704` and still fails `999,633 ppm`, and `K=12` still fails `49,804 ppm` while missing by `+3,983,788`; sampled max quotient bitlength is `23`. New narrow opening: a restoring-subtract extractor (`u -= v<<s`, quotient bit from borrow, add back on borrow) has a `2n` q-bit floor; with a generic barrel this projects `-174,766` to 3M and leaves only `43,691` one-way margin. The payload-bit primitive budget is only `641.65` CCX (`~2.51n`): ideal `2n` fits, `2.5n` barely fits (`-2,222`), but `3n` misses by `+170,322`, and current-style restoring/compare-sub primitives at `4n/5n` miss by `+515,410/+860,498`. A concrete current trial-subtract + masked-addback circuit measured `ccx64=258` (`4.03n`) and scaled gap `+526,194`, confirming existing primitives are not enough. A better signed-digit/non-restoring ledger avoids add-back entirely: one controlled add/sub per quotient payload bit plus a generic barrel gives p99 payload/count `336/118`, extraction one-way `357,888`, margin `132,817`, and projected gap `-531,268` to 3M. Charging nearest-quotient correction still leaves margin: floor payload p99 `325`, centered payload p99 `336`, correction bits p99/max `59/66`, extraction one-way `415,488`, projected gap `-244,516`. So centered Euclid is the better relaxed-2800q candidate, but it remains speculative until a phase-clean non-restoring signed-digit extractor/replayer is built and measured, including exact boundary cleanup, signed final correction, and denominator reverse.

### 1. New denominator-shift-free DIV recurrence

This is the best way to salvage the lessons from plus-minus without carrying its dead denominator cost. First explicit probe after the BY demotion is negative: bounded-quotient subtractive Euclid (`u <- u - qv`, `q <= 15`) avoids physical denominator shifts but explodes the reversible history/parser channel. Even with quotient computation free, q-history alone gives p90/p99 scratch `1608/6276` bits and max-step cap hits `20000`; `q <= 7` has p99 scratch `9514`. So “no shifts by using tiny quotients” is not the needed recurrence. A viable recurrence must avoid both physical denominator shifts and long per-step quotient history.

Requirements:

```text
persistent scratch <= 600--663 bits
all-in two-DIV + scale <= 2.06M CCX
per-update target <= 4.0k--4.6k CCX if ~404 updates remain
no dense bitlength/top-bit normalization predicate
no per-step 257-bit physical/barrel denominator shift
local reversible direction/control recovery
```

Examples worth probing only if they meet the budget gate up front:

- represent denominator scale purely as metadata and never normalize by data-dependent width,
- find a recurrence where the shifted operand is always a coefficient lane, not the denominator lane,
- fuse denominator scaling into the final product-clean channel so no explicit denominator normalization is needed.

### 2. BY selected/window-local denominator primitive

BY remains the best fully charged near-miss, but the missing primitive is precise. The first16/tail streaming-selector low-gate detour was checked and demoted: it fit 1425q scratch at `736` bits but a forward-only fresh tail carry update already exceeded the remaining low-gate Toffoli margin before cleanup/reversibility.

The remaining BY revival condition is:

```text
<= ~10k CCX per 16-step denominator window
explicit inverse/cleanup
no no-history consumed-denominator recovery
no full-ratio A inverse
no selected variable-coefficient row formation that costs field multiplication prices
```

If such a primitive exists, BY can re-enter active status. Otherwise more BY plumbing is not useful.

### 3. A genuine phase-clean in-place variable multiply/DIV primitive

Several algebraic point-add rearrangements become attractive if we can do something like:

```text
(x, y) -> (x, y/x)
```

or a product-clean multiply at near schoolbook cost with no field-sized history. Prior generic MBUC attempts were dense, so this needs a new structural idea, not another generic measure-old-multiplier attempt.

Payoff is high: deleting one Kaliski-scale inversion is still the largest lever. Risk is also high because many toy ANF probes already say generic cleanup phases are dense.

### 4. Solinas history-carry scale correction as a supporting optimization

The history-carry multihalve model can save meaningful cost and may help any route that produces scaled outputs, but it is not enough alone. Build it only when paired with a DIV route whose denominator/update body already passes the global budget gate.

### 5. Half-GCD / quotient-stream fusion

Raw payloads can be close to the scratch target, but every parser/tail attempt has failed so far. Revive only with a fused parser that consumes live denominator state without separate boundary/rank/live-prefix recomputation.

## Explicit stop list

Do not spend main-loop iterations on these unless a genuinely new primitive changes the premise:

- plus-minus per-step repeated shifts,
- plus-minus per-step physical barrel shifts,
- plus-minus denominator offset normalization by public schedule,
- plus-minus exponent-only normalization controls,
- BY full-ratio A-inverse selector,
- BY no-history consumed-denominator cleanup,
- centered/ordinary Euclid raw quotient streams without a parser breakthrough,
- curve-support or top-level MBUC cleanup as a free branch/reciprocal oracle,
- generic in-place multiply cleanup by measuring the old multiplier.
