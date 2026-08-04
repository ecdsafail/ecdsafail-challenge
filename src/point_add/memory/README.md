# Memory

Start with `06-research-status.md` for the current frontier, proved scope, failed approaches, open problems, and exact
re-entry commands. These notes live under `src/point_add` so they travel with submissions.

| file | contents |
|---|---|
| `01-architecture.md` | first-principles decomposition: the algorithm, the inversion, the qubit budget, the Toffoli budget |
| `02-lambda.md` | the intrinsic error rate — the hidden third score axis, and the reason the leaderboard stalls |
| `03-proven-floors.md` | where the headroom is NOT, with proofs (rank bound, multiplicative complexity, exact codec enumeration) |
| `04-traps.md` | four ways an env knob silently no-ops, positional addressing, validation gates |
| `05-qubit-reduction.md` | the measured qubit programme, including the exchange-rate trap |
| `06-research-status.md` | latest research handoff: certified baseline, scoped results, counterexamples, unresolved work, re-entry conditions |
| `07-syntactic-certification.md` | three closed routes to a stream-agnostic dead-gate certificate, each with a passing control; why the 46,286 never-firing gates need a semantic argument |
| `repro/` | compact tested programs retained to reproduce or extend the durable claims |
| [`repro/world_model.py`](repro/world_model.py) | executable evidence, invalidation, history-replay, and promotion-gate model |
| [`repro/hotness.rs`](repro/hotness.rs) | per-gate charge/fire census; self-gating — aborts unless it reconstructs the scorer's own Toffoli total |

The single most important operational fact: **a persistent-set reduction only pays if you lower `TLM_TARGET_Q` by the
same amount**, because the vent pool expands to fill whatever you free. The second most important: **only a
byte-identical `ops.bin` or a full 9024-shot run is evidence.** A healthy peak/Toffoli probe proves nothing.
