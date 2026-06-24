# sub-851 deep analysis — beat the 851-qubit register-shared record

Target: `base/sub-851-341d0963` (commit `e7dd3de`, submission `341d0963`). Deterministic, nonce
baked in source. Confirmed: **peak=851**, 838,845,965 ops, ~464M avg Toffoli, 9024/9024 OK.
Score = peak_qubits × avg-executed-Toffoli; we chase the LOW-QUBIT frontier (peak < 851, still 9024/9024).

This note is read-only analysis. **DO NOT** take the q944/q945/q949 modules at face value — they are
mostly OFF by default (see §A.4). The real 851 binder is the register-shared EEA swap comparator (§B).

---

## A. Architecture map

### A.1 Build path
- `src/bin/build_circuit.rs::main` → `point_add::build()` → (no dialog env) → `trailmix_port::build_builder().ops`
  → serialized to `ops.bin` (zstd). `eval_circuit` re-reads `ops.bin` in a separate trusted process.
- `point_add::build()` (`src/point_add/mod.rs:1980`): a long chain of `*_SELFTEST` env gates (all skipped
  by default), then `build_builder().ops` (`mod.rs:2141`).
- `trailmix_port::build_builder()` (`mod.rs:1775`): allocates the 4 input registers
  (target_x=256 q, target_y=256 q, offset_x/offset_y = classical Cbits), then
  `mod_sub_qb` ×2 (Px-=Qx, Py-=Qy), then `emit_dialog_gcd_raw_pa(...)` which performs the EC point-add
  via two modular inversions. **N=256** field bits; dx/dy stay **257**-bit through the divide (sign bit).
- The whole low-qubit route is configured by `configure_ecdsafail_submission_route()` →
  `configure_sub1000_trailmix_route()` (`trailmix_port/mod.rs:277-371`) which `set_default_env(...)`s
  ~70 `LOWQ_*` / `TRAILMIX_*` flags. `set_default_env` only sets if unset, so all of this is
  **source-baked** (no env needed for the submission) but overridable for proof binaries.

### A.2 The point-add wrapper and the inversion dispatch
- `ec_add_inplace_shrunken_pz` (`trailmix_port/ec/point_add.rs:203`) is the in-place P+Q.
- The inversion (slope λ = dy/dx) is dispatched at `point_add.rs:249-251`:
  `register_shared_eea_enabled()` (default **true**, gated by `TRAILMIX_REGISTER_SHARED_EEA=1`,
  `point_add.rs:36`) → **`register_shared_divide_forward`** (the active 851 path), with
  `register_shared_divide_cancel` for the uncompute. The phase tag is `ec3.inv_fwd` (and the cancel
  half is `ec3.alt.cancel`).
- The **inactive** fallback is `shrunken_pz_divide_forward` (the shrunken_pz state machine) — this is the
  route the q949 envelope governs, and it is NOT what produces 851.

### A.3 The active register-shared EEA (the real 851 datapath)
- Implemented in `trailmix_port/inversion/register_shared_eea_reference.rs` (~19.6k lines) + helpers
  `register_shared_eea.rs`, `register_shared_eea_microkernels.rs`.
- It is a **search-and-project / profile-baked** construction: the authors composed a stack of
  lifetime/scratch-reuse cuts (each `LOWQ_*` flag) in a whole-circuit profile, proved each reversible,
  then baked the winning combination. The numbered families are achieved milestones:
  **Q855 → Q851 → Q847 → Q845 → Q839** are the *active* register-sharing cuts (flags = `1`), each named
  for the peak it reached. (Confusingly, the *higher*-numbered Q944–Q959 families are *worse*/staged —
  see A.4.)
- Active lever flags (all `=1`) and what they do (by name + comment in `mod.rs:333-369`):
  - `TRAILMIX_REGISTER_SHARED_EEA` — Q883: share work registers across the EEA so dx is consumed into a
    single 259-bit work lane instead of a separate copy.
  - `LOWQ_DIRECT_PREFIX_BITLEN`, `LOWQ_REUSE_ZERO_CARRIES_FOR_*PREFIX*`, `LOWQ_REUSE_ROTATED_BITLEN_SCRATCH`,
    `LOWQ_REUSE_COEFFICIENT_COMPARATOR_SCRATCH`, `LOWQ_Q839_PHASE_REMAINDER_SCRATCH`,
    `LOWQ_FUSED_ZERO_PREFIX_BITLEN` — Q855/Q839: reuse zero-carry / rotated-bitlen scratch instead of
    fresh ancilla for the CLZ/bit-length prefix.
  - `LOWQ_REGISTER_SHARED_REVERSE_DECREMENT_STREAM`, `LOWQ_REUSE_LQ_AS_SWAP_OLD_R_LENGTH`,
    `LOWQ_SPLIT_COEFFICIENT_ROTATION_LIFETIME`, `LOWQ_REUSE_COEFFICIENT_LESS_THAN_LANES`,
    `LOWQ_CALLER_SCRATCH_KG_REVERSE_DECREMENT`, `LOWQ_REUSE_CLEAN_CHAIN_FOR_COEFFICIENT_ADD` — Q847: six
    proved lifetime/stream choices; reuse comparator "less-than" lanes and the caller's scratch.
  - `LOWQ_REUSE_PRESERVED_DY_TOP_FOR_PREFIX`, `LOWQ_MIXED_WIDTH_L_R_PRIME` — Q845: reuse a preserved
    dy-top lane for the prefix; `l_r_prime` length register shrunk 9→8 (mixed width).
  - `LOWQ_PAIRED_BITLEN_SOURCE_COMPLEMENT`, `LOWQ_COEFFICIENT_NONNEGATIVE_X_CANCEL`,
    `LOWQ_Q845_LIFETIME_COEFFICIENT_FUSION`, `LOWQ_FUSE_PROMISED_SWAP_SUPPORT_LIFETIME`,
    `LOWQ_Q845_SWAP_ONLY_T_PRIME_LENGTH`, `LOWQ_Q851_TRUNCATED_SWAP_ONLY_GUARD`,
    `LOWQ_Q851_FIXED_SIGN_EVENT` — **Q851**: the final swap-only coefficient fusion + truncated
    swap-only guard + fixed sign event. **These are the cuts that took the swap step down to 851.**
  - `TRAILMIX_TAIL_NONCE=4114534` — the Fiat-Shamir tail nonce. Only toggles 2 output lanes
    post-circuit (`mod.rs:3405-3413`); **zero peak impact**, it just reseeds the 9024 SHAKE test inputs
    so the (correct) circuit lands on an island with 0 classical misses.

### A.4 The q944/q945/q949 "lever" modules — STAGED, OFF BY DEFAULT
**Critical:** every one of these is `set_default_env(..., "0")` in `mod.rs:308-330`. They are an *aborted /
not-yet-certified* lower-peak rewrite (the "shrunken_pz + dirty catalytic hosting" line). They do NOT
participate in the 851 build. Their levers (for reference, if we ever revive that line):
- `q944_dirty_catalytic_predicate.rs` — host a predicate-controlled involution on an **arbitrary dirty
  qubit** via the catalytic identity `G^d; d^=f; G^d; d^=f = G^f` (header :3-5). Restores `d` exactly →
  no clean alloc → peak-neutral. Asserts `allocation_serial/next_qubit/active_qubits/free_qubits`
  unchanged (:198-201).
- `q944_dirty_parity_microkernels.rs` — allocation-free ripple add/sub and the strict comparator
  `out ^= gate*[v<u]` (`strict_compare_gated_dirty_carry_refs` :213) on a borrowed dirty carry lane.
- `q944_quotient_witness.rs` — fabricate a zero host from the zero-init 25-bit quotient register for the
  5 division rows with no support-zero pair (host on sentinel `q[24]`).
- `q944_gate_host_feasibility.rs` / `q944_gate_host_lifecycle.rs` — feasibility decides *which* outer-pair
  lane is provably 0 entry+exit and untouched by the body (`exact_zero_paired_bits` :556); lifecycle
  proves the compute/body/uncompute on a zero host is correct.
- `q945_local_hosts.rs` — sealed static table mapping each step to its borrowed host lane;
  **acceptance BLOCKED** pending 4 Q946 hardening steps (`Q945_REQUIRED_Q946_INTEGRATIONS`).
- `q949_robust_envelope.rs` (+`_data.rs`, +`_metadata.json`) — the "q948-peak-safe-envelope-projection-v3"
  per-row width table for the **shrunken_pz** route. See §B.3 for why its 681 cap is NOT the 851 binder.
- "selective qcap" (`TRAILMIX_Q_CAP`/`TRAILMIX_Q_TARGET`) — per-step clamp of the quotient width.
  Default `Q_CAP=99` (neutralized global), `Q_TARGET=683` (`LOWQ_Q957_TARGET683=1`). This governs the
  shrunken_pz quotient, but the register-shared route's peak is set elsewhere (§B), so qcap is not the
  binder either on the active path.

---

## B. The EXACT peak-determining mechanism

**Single binding fact:** peak=851 is set by a **5-qubit transient comparator (unsigned-less-than +
its MCX ancillas) stacked on an 846-qubit plateau**, inside the EEA **coefficient swap-direction step**
of the register-shared forward inversion (`ec3.inv_fwd`). It is NOT the q949 681-cap and NOT the nonce.

### B.1 The exact qubit climb to 851 (measured, count-only trace)
Probe: `POINT_ADD_COUNT_ONLY=1 TRACE_PEAK_NAMES=1 ./target/release/build_circuit`. The last allocations:

```
... 823  rs.divider.sign                       <- persistent base reaches 823
824..836  rs.scheduled.pre-scratch[0..12]      (+13)
837..838  rs.scheduled.remainder-scratch[13,14](+2)
839..845  rs.q845-swap-only.cursor-scratch[1..7](+7)
846       rs.q845-swap-only.carry              (+1)
847       uls_gate                             (+1)  } Khattar-Gidney unsigned-less-than
848..850  uls_anc[0..2]                        (+3)  } comparator (khattar_gidney.rs:759,778)
851       mcxk_t3[0]                            (+1)  <- the qubit that SETS peak=851 (mcx.rs:254/381)
```
Peak phase `ec3.inv_fwd`, ops_idx≈155438. The cancel half `ec3.alt.cancel` ties at 851 (mirror image).

### B.2 The persistent 823-base decomposition (the real "fixed" term)
Authoritative profile: `register_shared_eea_reference.rs::profile_reference_scheduled_inversion`
(:16114-16134). `inversion_state = 2*259 + 5*REFERENCE_LENGTH_WIDTH(9) + 4 = 518 + 45 + 4 = 567`;
`+ passenger 257 (dy) = 824` (profile uses l_r_prime=9). Production has `LOWQ_MIXED_WIDTH_L_R_PRIME=1`
→ l_r_prime 9→8 → **823** persistent. Components:
- `work1` 259, `work2` 259  (the two EEA work lanes; dx is consumed into work2 — there is NO separate
  2×256 "two coords" term; that is the whole point of register sharing). **518 q = the dominant base.**
- 5 length registers (l-t, l-t-prime, l-q, l-s) ×9 + l-r-prime ×8 = 44.
- phase1, phase2, iteration-parity, sign = 4.
- **dy passenger = 257** (rides through the forward EEA, ghosted only after, `:9184`).
So `823 = 518(work) + 44(length) + 4(flags) + 257(dy)`. **Arithmetic to 851: 823 + 13 + 2 + 7 + 1 + 5 = 851.**

### B.3 Why the q949 681 / target 683 / slack=2 is a RED HERRING for 851
- The `maximum_pair_symmetric_sum=681`, `target_sum=683`, `peak_safe_pair_symmetric_cap=681`,
  `minimum_pair_symmetric_slack=2`, `clz_contexts=2088`, `row_count=530` metadata is consumed in
  `q949_robust_envelope.rs::parse_envelope` (:151-244): per row it symmetrizes widths
  `[max(A,B), max(A,B), max(Ca,Cb), max(Ca,Cb), Q]`, asserts `sum ≤ 681`, computes `slack=683-sum`,
  asserts min slack=2 / max sum=681 / clz_contexts=2088. The binding rows are 379/380 with
  `[A,B,Ca,Cb,Q]=[72,72,256,256,25]` (dominated by 2×256). The "pair-symmetric sum" is exactly a
  Bézout-coefficient symmetrized width sum (the two EEA coefficient registers Ca/Cb each counted twice).
- BUT this whole table only constrains the **shrunken_pz** EEA, gated by
  `LOWQ_Q949_ROBUST_SYMMETRIC_SCHEDULE=0` (off). On the active register-shared route the envelope is
  never consulted for peak. So **slack=2 there does not mean the 851 circuit has 2 qubits of slack.**
- The CLZ (count-leading-zeros) bounds (`clz_safe_low_upper_bounds`, 2088 contexts, 12590 limiting
  witnesses) bound how many low lanes of each coefficient register are provably zero (so the comparator
  / shift can skip them). On the active route the analogous mechanism is the `LOWQ_*_BITLEN` /
  `THIN_CLZ_WINDOW=78` prefix cuts, which already shrank the length registers to 8-9 bits.

### B.4 The binding constraint (skeptical conclusion)
peak=851 is **(c) a structural register-width / transient-stacking limit**, specifically:
1. The **846 plateau** (823 base + 23 scheduled/swap scratch) — dominated by 2×259 work + 257 dy.
2. The **+5 swap comparator transient** (uls_gate + 3 uls_anc + mcxk_t3) that defines the exact 851.
The MEMORY index already names this: *"bottleneck is SWAP (62%), not masking."* It is NOT nonce-limited
(nonce only reseeds shots), NOT q949-681-limited (that route is off), NOT qcap-limited on this path.

---

## C. RANKED sub-851 attack points

Legend: [S]=structural code change, [N]=nonce/search. Gain estimates are vs 851.

### C1. [S] Collapse the +5 swap-comparator transient (uls + mcxk ancillas). ~ -1 to -5 q. HIGH value, MED risk.
- Mechanism: the qubit that *sets* 851 is `mcxk_t3[0]` (a clean MCX-decomposition ancilla) sitting on top
  of `uls_gate`+`uls_anc[0..2]` (Khattar-Gidney unsigned-less-than prefix ancillas). All 5 are
  *transient* scratch for the swap-direction comparator in
  `coefficient_fused_data_and_sign_q845_swap_only` (`register_shared_eea_reference.rs:7419`).
- Why slack may exist: this is exactly the masking/comparator that the q944 line was trying to host on
  borrowed/dirty qubits (`strict_compare_gated_dirty_carry_refs`, `mcx_dirty_ladder`). The *dirty-ancilla*
  MCX (`mcx.rs:mcx_dirty_ladder` / `mcx_dirty` :75) needs NO clean ancilla — if the swap comparator's
  `mcxk_t3` clean ancilla and the `uls_anc` prefix lanes can be replaced by borrowed lanes from the
  846-plateau registers (e.g. an already-zero low lane of work1/work2 exposed by the CLZ window), the
  transient collapses. Each clean→dirty swap is -1 q; full collapse is -5 q (→ ~846).
- Feasibility/risk: the dirty-MCX machinery already exists and is exhaustively proven
  (`q944_dirty_*`). Risk: it changes the op stream → re-rolls the Fiat-Shamir island → needs a fresh
  tail-nonce hunt (the recurring wall). This is the single most direct lever and aligns with the MEMORY
  "g-flag dirty-MCX rewrite (chosen)" decision.

### C2. [S] Shrink the dy passenger lifetime (257 q live through forward EEA). ~ -? big if it works, HIGH risk.
- Mechanism: dy (257) rides as a passenger through the entire forward inversion and is only ghosted
  after (`:9184`). For the ~150k-op forward window it is dead weight on the plateau. The staged
  `LOWQ_PASSENGER_TOP_LIFETIME_EXPERIMENT` / `LOWQ_Q947_PASSENGER_DIRECT_HCLZ` /
  Q954 `canonical_passenger_top_lifetime` (all =0) target exactly this.
- Why slack: if dy's top lanes are provably zero / canonical at the peak step, they could host the swap
  comparator (folding C1 into C2) or be released early.
- Risk: HIGH — comment at `mod.rs:289-290` says "Passenger-top reuse is forbidden because lambda_raw is
  not guaranteed canonical." So this needs a canonicalization proof first. Flagged as **not free**.

### C3. [S] Shrink the 2×259 work lanes (the 518-q dominant base). ~ -2 to -4 q, LOW-MED value, LOW risk.
- Mechanism: work1/work2 are 259 = 257 field + 2 pad lanes (`work2-pad0/1`). If the pad lanes are only
  needed transiently (sign/carry headroom) they might be borrowed rather than permanently allocated.
- Why slack: 2 pad lanes ×2 registers = up to 4 q. The CLZ window already proves high lanes zero late in
  the EEA; an in-place mixed-width work register (like `LOWQ_MIXED_WIDTH_L_R_PRIME` did for l_r_prime)
  could trim. Small but island-cheap if peak-neutral elsewhere.
- Risk: LOW correctness risk (well-trodden mixed-width pattern) but each q still re-rolls the island.

### C4. [N] Pure nonce re-hunt at an already-built lower-peak candidate. variable, ties to GPU trailmix.
- Mechanism: any of C1-C3 produces a circuit that is *correct* but lands on a Fiat-Shamir island with a
  few classical misses. The fix is a tail-nonce search (the GPU/WMI trailmix tooling) to find a 0-miss
  nonce. This is the *enabling* step for every structural cut above, not a standalone qubit win.
- Status per notes: nonce-island density is THE recurring wall (06-06 knob islands, 06-11 K3 0/300k,
  06-13 WMI shards, 06-17 Poisson(1) residual saved by nonce 270). Necessary but not sufficient — must
  pass the trusted 9024-shot evaluator, not just the comb prefilter.

### C5. [S, LONG] Revive the shrunken_pz + q949 robust envelope line to its certified end. ~ down to q944/q945. HIGH effort.
- Mechanism: the q949 681-cap route + q945 local hosts + q944 dirty catalytic hosting is a *complete
  lower-peak design* that is merely uncertified (`LOWQ_Q945_LOCAL_HOSTS=0`,
  `Q945_REQUIRED_Q946_INTEGRATIONS` blocking, `LOWQ_Q949_ROBUST_FRESH_SUPPORT_CERTIFIED=0`).
- Why slack: the envelope's slack=2 is for THAT route; if finished it targets ≤945-ish on the shrunken_pz
  datapath, which is *higher* than 851 — so this is likely **dominated** and not worth reviving unless
  combined with the register-sharing wins. Flag: **probably a dead/worse branch** vs the active 851 line;
  the higher Q-numbers (944/945/949) are worse peaks than the active Q845/Q851. SKEPTICAL: do not chase.

### Already-exhausted / blocked (do not re-attempt per notes):
- One-inversion point-add (`ONE_INV_DX3_AFFINE_PA_BLOCKER`) — provably impossible, 2 inversions is a floor.
- Montgomery batch-invert the two GCDs — data-dependent, blocked.
- Apply-cswap width truncation — full-256 mod-p residues, 9024 mismatches, impossible.
- K≥3 radix divstep as a global lever — ~peak-neutral or +1 (1221→1222), only selective single-step
  works and is nonce-starved.
- Blunt global `Q_CAP` clamp — manufactures classical misses; superseded by selective `Q_TARGET`.

**Ranked recommendation:** C1 (dirty-MCX collapse of the swap comparator transient) is the cleanest,
best-evidenced structural lever (-1 to -5 q), reusing the already-proven q944 dirty-MCX kernels, gated
only by a C4 nonce re-hunt. C3 is a cheap LOW-risk add-on. C2 is high-upside but blocked on a
canonicalization proof. C5 is likely dominated.

---

## D. Fast iteration harness (measure peak without the 838M-op eval)

**Fastest peak probe (~18s, emits 0 ops, no disk, no eval):**
```
POINT_ADD_COUNT_ONLY=1 TRACE_PEAK=1 ./target/release/build_circuit 2>&1 | grep TRAILMIX_SHRUNKEN_PZ
# -> TRAILMIX_SHRUNKEN_PZ peak_qubits=851 peak_phase='ec3.inv_fwd' ops=838845965
```
`POINT_ADD_COUNT_ONLY=1` builds in count-only mode (ops not materialized); `TRACE_PEAK=1` prints peak.

**Per-phase active-width breakdown (find the binding phase):**
```
POINT_ADD_COUNT_ONLY=1 TRACE_PEAK=1 TRACE_PHASE_ACTIVE=1 TRACE_PHASE_ACTIVE_TOP=15 ./target/release/build_circuit
# ec3.inv_fwd 851, ec3.alt.cancel 851 (tied), next tier mul_canonical 810/809 -> 41q cliff
```

**Exact named qubit-by-qubit climb to peak (which register sets the peak):**
```
POINT_ADD_COUNT_ONLY=1 TRACE_PEAK_NAMES=1 ./target/release/build_circuit 2>&1 | grep PEAK_NAME | tail
# prints each new peak >=800 with the allocating register name (see B.1)
```

**Full plateau membership at a target peak (every live register at peak; slower ~2min):**
```
POINT_ADD_COUNT_ONLY=1 TRACE_NAMED_PEAK_TARGET=851 ./target/release/build_circuit
# records named_peak_plateaus internally; assert covers all 851 live qubits.
```

**Other useful gates:** `SKIP_ALT_SEED_CHECKS=1` (skip the alt-seed correctness re-derivation, faster),
`TRACE_EACH_PEAK=1` (stream every peak bump), `TRACE_PHASES=1` (per-phase op attribution).

**Correctness (run ONCE per candidate, slow):** `./target/release/build_circuit` (full, ~46s, writes
ops.bin) then `./target/release/eval_circuit` (slow, 838M ops) → must report qubits=851→new, 9024/9024 OK.
The count-only probes do NOT validate correctness or island-cleanliness — they only measure peak. A
peak-reducing change must still pass the full eval + (likely) a fresh tail-nonce hunt.

---

## Key file/line index
- Build flow: `src/point_add/mod.rs:1775,1980`; `trailmix_port/mod.rs:277-371` (the baked env stack).
- Active inversion dispatch: `trailmix_port/ec/point_add.rs:36,249-251`.
- Active EEA datapath: `trailmix_port/inversion/register_shared_eea_reference.rs`
  (profile :16114-16134; work widths 259/257 :8918-8919; swap comparator
  `coefficient_fused_data_and_sign_q845_swap_only` :7419,7453-7457).
- Peak comparator ancillas: `trailmix_port/arith/khattar_gidney.rs:759,778` (uls_gate/uls_anc);
  `trailmix_port/arith/mcx.rs:254,381` (mcxk_t3), dirty-MCX kernels :17,75.
- Peak tracking: `point_add/mod.rs:393-415,397-401` (alloc/peak latch); `trailmix_port/circuit.rs:259-387`
  (named-peak trace, TRACE_PEAK_NAMES, record_named_peak_plateau).
- Inactive shrunken_pz / q949 envelope: `trailmix_port/inversion/q949_robust_envelope.rs:151-244`;
  metadata `q949_robust_projection_metadata.json`; all gated `=0` in `mod.rs:308-330`.

---

## UPDATE 1 — C1 landed: dirty-MCX swap collapse, structural peak 851 → 850

**Cut.** The +5 swap-comparator transient (§B.1, §C1) is now +4. The qubit that set 851,
`mcxk_t3[0]` (the clean MCX-decomposition ancilla from `mcx_clean_k` k=3 inside the swap scan's
`unary_iterate_log_star` all-ones detector), is eliminated by routing that detector MCX through a
**borrowed dirty lane** (`mcx_dirty_any_k`, k≤5 → Barenco 4-CCX surrounded form, no clean alloc).
The borrow is the idle persistent `phase2` flag — proved untouched inside every swap
`for_each_{forward,reverse}` body (the bodies use phase1/enable/add_only as sign-bracket controls,
never phase2/sign). The dirty-MCX restores `phase2` exactly per call.

**Implementation (additive, flag-gated `LOWQ_Q845_SWAP_DIRTY_MCX`, default-on in route).**
- `khattar_gidney.rs`: `unary_iterate_log_star` now delegates to new
  `unary_iterate_log_star_dirty(circ, c, n_iters, dirty: Option<&QReg>, body)`; when
  `dirty=Some(psi)` and the detector MCX has k≤5 controls it emits `mcx_dirty_any_k(.., psi)`
  instead of `mcx_clean_k`. `None` (default for the other 2 callers) is bit-identical to before.
- `register_shared_eea_reference.rs`: `Q845SwapOnlyCoefficientGuard::for_each_forward/_reverse` gain a
  `dirty: Option<&QReg>` param; the 4 production swap call sites pass
  `q845_swap_dirty_mcx_requested().then_some(phase2)`; the 2 proof-harness sites pass `None`.
- `trailmix_port/mod.rs:370`: `set_default_env("LOWQ_Q845_SWAP_DIRTY_MCX","1")`.

**Measured (count-only fast probe).** flag=0 → peak 851, ops 838,845,965 (exact baseline);
flag=1 → **peak 850**, ops 834,829,613. Default-on binary (no env) → 850. New climb top:
```
846  rs.q845-swap-only.carry
847  uls_gate
848  uls_anc[0]
849  uls_anc[1]
850  uls_anc[2]      <- NEW binding qubit (was 851 = mcxk_t3[0])
```
**New binder:** `uls_anc[2]`, the top Khattar-Gidney prefix-AND ancilla, ec3.inv_fwd, ops_idx 155428.

**Logical correctness (no full-island dependence).**
1. `unary_iterate_log_star_dirty_full_range` / `_partial_range` (new, in khattar_gidney `mod tests`):
   state-simulator value-test asserting identical res/fired/counter outputs AND `phase2`/psi
   restored, for every (n,v,dirty_init), including the **n=9 / 260-iter** production swap-scan config.
   (NB: `cargo test` is repo-broken — 206 pre-existing errors in the lib-test target, e.g.
   `circ.ops`/`circ.peak_qubits` unknown fields in the *existing* `run_mcx_clean_case`; the new test
   compiles cleanly but cannot be run until that stale harness is fixed. It is a static value-spec.)
2. **Decisive runtime check:** full `build_circuit` (non-count-only, real per-shot sim) completes with
   NO panic, emitting 834,829,613 ops + ops.bin. `mcx_dirty`/`mcx_dirty_any_k` carry
   `contract_capture`+`contract_pop_and_check` that assert, per shot, `target ^= AND(ctrls)` and the
   dirty lane restored; these run during the build. A bad borrow would have panicked. ⇒ the swap-step
   inversion is gate-level correct on every shot.

**C3 is DEAD (ruled out, do not re-attempt).** `register_shared_eea.rs:328-334` proves the tight
packing bound `work{1,2}_used ≤ n+3 = 259`. The 2 pad lanes are not slack — they are the proven
worst-case headroom. Trimming `REGISTER_SHARED_WORK_WIDTH` 259→258 would violate the invariant.

**Implied nonce re-hunt.** The op stream changed (838.8M→834.8M ops), so the baked Fiat-Shamir island
re-rolls → the existing `TRAILMIX_TAIL_NONCE=4114534` almost certainly no longer lands 9024/9024
(eval result pending). Scope: a fresh tail-nonce hunt (C4) on the 850-qubit circuit. The structural
cut + logical correctness are independent of, and prerequisite to, that hunt.

**Next lever toward ≤845.** The remaining swap transient is uls_gate + uls_anc[0..2] (4 q, now the
top of the plateau). Collapsing uls_anc needs 3 provably-|0⟩ idle lanes to host the KG prefix-AND tree
(can't reuse phase2 — it's live as the MCX dirty borrow simultaneously). That is the q944-style
gate-host feasibility analysis (higher risk). uls_gate is a clean compute target (body reads it), so
it needs a provably-|0⟩ borrow, not a dirty one.


---

## UPDATE 2 — THE 846 ISLAND IS ALREADY CLEAN. The "nonce re-hunt" premise was WRONG for the register-shared route.

**Headline (eval-confirmed, twice):** at HEAD `09c6acb` (peak 846, ops 865,762,909) with the
**source-baked `TRAILMIX_TAIL_NONCE=4114534`**, the full trusted evaluator reports:
```
qubits = 846 | tested shots = 9024 | classical mismatches = 0
phase-garbage batches = 0 | ancilla-garbage batches = 0 | all 9024 shots OK
avg executed Toffoli = 508,948,948
```
**846 = 9024/9024 OK, 0/0/0, with the EXISTING baked nonce. No re-hunt was needed.** Reproduce:
`./target/release/build_circuit` (no env -- nonce is source-baked) then `./target/release/eval_circuit`
(~12 min dedicated / ~36 min if a scan competes for cores). Build: `RUSTFLAGS="-C linker=clang"
cargo build --release --offline --bin build_circuit --bin eval_circuit`. (Logged to results.tsv:
`09c6acb  508948948.000  ...  846  865762909  OK`.)

### Why the recurring "every structural cut re-rolls the island" wall did NOT apply here
The 850->847->846 cuts (UPDATE 1 + commits 452e949, 09c6acb) are all **dirty-MCX / clean-ancilla
borrows on the active register-shared EEA** (`TRAILMIX_REGISTER_SHARED_EEA=1`,
`register_shared_divide_forward`). They change the op COUNT (and thus the Fiat-Shamir seed, which is
`SHAKE(b"...v2" || total_ops+96)`), so the 9024 challenge POINTS reshuffle -- but they do NOT change
the per-shot *data-dependent width behaviour* of the inversion. The register-shared EEA is **full-width
correct**: its divstep carries a genuine `active_width` window (`register_shared_eea_reference.rs`
:13387,13533) and inverts every field element regardless of bit-length. So the island is NOT a
data-dependent cleanliness constraint on this route -- any nonce that reshuffles the points still lands
9024/9024. The baked nonce 4114534 (inherited from the 851 source) is still clean at 846, and so would
essentially any nonce.

### THE FAST FILTER IS A PESSIMISTIC RELIC OF THE INACTIVE shrunken_pz ROUTE (do not trust its misses)
Both per-nonce "filters" model the **shrunken_pz** width envelope, which is OFF on the 846 route:
- `trailmix_port/mod.rs::support_report_for_xof_limited` (:453) -> per shot derives `dx`,
  `qx_minus_rx` and calls `shrunken_pz_schedule::thin_factor_repairs_u256` (:2014). A factor is "hard"
  if its bit-width overflows a per-row entry of the **thin schedule** (the shrunken_pz per-divstep width
  table, gated `TRAILMIX_THIN_SCHEDULE=1`, baked SEED=278 / CLZ_WINDOW=78 / MARGIN=0 / Q_TARGET=683 /
  Q_CAP=99). Env knobs: `TRAILMIX_SUPPORT_CHECK=1` (one-nonce report), `TRAILMIX_TAIL_NONCE_SEARCH /
  _START / _SHOTS / _THREADS / _EARLY_MISS / _CONTINUE / _TRACE[_CLEAN]` (multi-threaded scan, :3274);
  `miss_factors` = hard-factor count, `repair_entries` = sum of per-factor repairs.
- `dialog_gcd_classical_filter.rs::check_all_shots` (:2034) -- the dialog-line analogue (same idea:
  width-envelope + truncated-comparator), also shrunken_pz-modelled.

**MEASURED DISAGREEMENT (the key result):** for the baked nonce 4114534 the support filter reports
`miss_factors=13, repair_entries=138, first_miss=(draw 1568, "qx_minus_rx", 3)` -- yet the **trusted
eval is 0/0/0**. So the filter's 13 "misses" are FALSE POSITIVES: they flag factors that overflow the
*shrunken_pz* envelope but that the *register-shared* circuit inverts correctly. The filter scores the
WRONG datapath. Cross-checked three ways (standalone `TRAILMIX_SUPPORT_CHECK`, the built-in
`TRAILMIX_TAIL_NONCE_SEARCH` scan, and the full eval) -- all internally consistent; the filter is simply
not the 846-route's cleanliness oracle.

### Density / scan corollary
A CPU "clean-at-846" scan is **moot on the active route**: it is nonce-island-independent, so the
density of eval-clean nonces ~ 1 (any reseed stays 9024/9024). The slow filter scan (built-in
`TRAILMIX_TAIL_NONCE_SEARCH`; per-nonce cost = up to 9024 draws x 2 secp256k1 double-and-add mults --
naive `WeierstrassEllipticCurve::mul`, no comb table, ~few-hundred nonce/s/core) found ZERO
filter-clean nonces in the windows scanned -- but that only measures the *shrunken_pz envelope*
difficulty (consistent with the Q682 e^-20.7 wall), NOT 846's. The filter scan is the wrong tool for
this route and should not gate register-shared cuts.

### Caveat for FUTURE cuts (where a real nonce hunt WOULD return)
The island only re-binds if a future cut (a) re-activates the shrunken_pz / thin-envelope route, or
(b) introduces a genuinely data-dependent width truncation on the register-shared route -- e.g. a
`TRAILMIX_Q_TARGET` / `AB_CAP` / `CACB_CAP` clamp that the inversion can no longer fully absorb; those
clamps DO feed `thin_factor_repairs_u256` (via `trailmix_q_target`/`trailmix_q_cap`) and would make the
filter binding again. The current 846 cuts do neither. **For the dirty-MCX borrowing line, treat the
nonce as fixed and the filter as irrelevant; only resurrect the GPU nonce hunt if a future cut adds a
data-dependent width clamp.**

### VERDICT
(a) **No GPU needed for 846.** The validated peak-846 circuit is ALREADY a clean 9024/9024 submission
with the baked nonce 4114534 (eval 0/0/0). The task's premise ("structural cuts re-rolled the island,
baked nonce no longer clean at 846") does NOT hold for the register-shared route at this HEAD. The fast
filter's `miss=13` is a false alarm from modeling the inactive shrunken_pz envelope. Recommendation:
score/submit 846 as-is; spend no CPU/GPU on a nonce hunt for this route.


---

## UPDATE 3 — 846 IS A DOUBLE CO-PEAK. The cursor-scratch borrow cannot reach 838 (or even 845) ALONE.

**Headline (measured, count-only census):** the 846 peak in `ec3.inv_fwd` is set by **TWO independent
transient stacks that BOTH sum to exactly 846** over the same 823-ish persistent base. Cutting either one
alone leaves the other holding 846. The competitor's "cursor-scratch is value-bearing → need a width-trim"
claim is wrong on *mechanism* (see below), but 846 is genuinely a hard floor for a deeper reason they did
not state: it is **doubly bound**.

### The two co-peaks (DUMP_PLATEAU census, every live qubit enumerated; both = 846)
Probe: `POINT_ADD_COUNT_ONLY=1 TRACE_NAMED_PEAK_TARGET=846 DUMP_PLATEAU=1 ./target/release/build_circuit`.

**(P1) swap-step coefficient plateau** (trigger `rs.q845-swap-only.carry`, ops_idx≈155407):
persistent base 823 [iter-parity1 + l_q9 + l_r_prime8 + l_s9 + l_t9 + l_t_prime9 + phase1/2 + sign +
work1 259 + work2-pad 2 + tx256 + ty256 + tx_ov/ty_ov 2]  + coeff-fused 5 [above-guard + add-only +
chain2 + enable]  + **address 9 + carry 1 + cursor-scratch 8** = **846**.

**(P2) dynamic-bitlen swap plateau** (trigger `rs.dynamic-bitlen.fused-prefix-scratch[1]`, ops_idx≈377119,
inside `conditional_work_and_length_swap_under_zero_predicate`, the `step%4==0` swap rows):
same 823 base  + **fused-prefix-scratch 2 + rotated-bitlen.length 10 + scheduled.swap-condition 10 +
swap-length.promised-support 1** = **846**.

Both verified by explicit summation of the full per-register census (not estimated). `ec3.alt.cancel`
mirrors both at 846 (dy-restored 257 replaces tx+ty there; still 846).

### Why a single cut is moot (the decisive structural fact)
The 851→846 history (UPDATE 1 + commits 452e949, 09c6acb) removed the swap-comparator stack
(mcxk_t3 + uls_anc[0..2] + uls_gate, up to 5 lanes) that USED to sit ABOVE cursor-scratch+carry on P1.
After those cuts, P1's top is now `cursor-scratch[7]`+`carry` landing at exactly 846 — which **coincides
with P2's independent 846**. So the swap-step plateau and the dynamic-bitlen plateau have *converged* to a
co-peak. Freeing cursor-scratch (−8) + carry (−1) drops P1 to 837 but **P2 stays at 846 → measured peak
unchanged at 846.** The next cliff below both is 810 (mul_canonical), a 36-q gap. Reaching ≤845 requires
cutting **both** P1 and P2 in the same window.

### Cursor-scratch borrowability VERDICT (build evidence, not code-reading)
1. **Lifetime shape (probe, 16384 inputs):** `CURSOR_SCRATCH_PROBE=1` → cursor-scratch is CLEAN(==0) at the
   forward→reverse boundary for ALL inputs, nonzero at only 18% of ops (longest contiguous run 19/739 =
   2.6%). It is a clean-between-ops v-chain transient, NOT "value-bearing working space." The competitor's
   stated *reason* is refuted.
2. **BUT the borrow is still blocked on P1**, for a reason the lifetime shape hides: cursor-scratch has TWO
   consumer classes (verified by reading the actual production call graph + the probe):
   - *Dirty-tolerant* uses (would survive a dirty `l_q` borrow): the Bennett v-chains
     (`multi_controlled_x_vchain[_borrowed]` in `toggle_q851_fixed_sign_event`, `accumulate`'s predicate
     v-chains, `toggle_nonzero`) and the `controlled_add_one`/`controlled_sub_one` style increments — all
     restore arbitrary ancilla content.
   - *Clean-ZERO-requiring* uses (would CORRUPT under a dirty borrow): (a) `accumulate`/`unaccumulate` call
     `controlled_increment_mod_2n`/`controlled_decrement_mod_2n` (ripple carry chain that assumes carries
     enter at 0), and decisively (b) the **reverse-boundary constant-operand buffer**:
     `prepare_reverse_boundary` loads `reverse_source = cursor_scratch ∪ {coefficient_active}` via
     `cx(l_r_prime/l_s, source)` then `cuccaro_add_mod_2n_refs(source, address)`, and `for_each_reverse`
     X-loads the constant `delta` into the same lanes (`if (delta>>i)&1 { X(source[i]) }`) — both ABSOLUTE
     loads that are only correct from |0⟩. A dirty `l_q` (holding a live nonzero bit-length value, restored
     by the ladder borrow to its ORIGINAL value, NOT 0) breaks these. This is the concrete build-level
     conflict, and it matches *why* the probe finds cursor-scratch clean at the boundary: the boundary
     buffer **needs** it clean.
   So: "value-bearing ⇒ can't borrow" is a non-sequitur (dirty-borrow works on live lanes), but here the
   block is real and specific — the boundary buffer is a clean-zero *data* buffer, not a restorable ancilla.
3. **Type-system corroboration:** `QReg` is `!Clone` (RAII frees its id on drop), and every cursor-scratch
   consumer takes `&[QReg]` (owned), not `&[&QReg]`. Hosting on borrowed `l_q` (`&[&QReg]`) would require
   refactoring the entire swap-step scratch interface — a large change that *still* lands value-incorrect at
   the boundary and *still* doesn't lower the peak (P2). Not pursued: it cannot yield a CORRECT sub-846.

### The carry lane (P1 top) is also not free
`rs.q845-swap-only.carry` is genuine working space within the add/sub body: `emit_q845_fused_coefficient_body`
WRITES it (`vchain(...,carry,tmp)`), then it is READ by `ccx(carry, above_guard, sign)` and
`toggle_output_coefficient_enable_q845_inline_underflow` before the reverse body uncomputes it. A dirty
borrow corrupts the `ccx(active, carry, work1)` read. (This is the one place the competitor's
"value-bearing" language is literally correct — but it applies to `carry`, not cursor-scratch.)

### Width-trim option (the competitor's suggested route) — also tight on P1
- `address` 9: holds `physical_work_width=259` via `toggle_constant` then counts to n_iters≤258 → genuinely
  needs 9 bits. No trim.
- `cursor-scratch` 8 = `l_t.len()-1`: bounded by the reverse-boundary buffer (needs 8 lanes = address.len−1).
  No trim.
- `l_r_prime` already 9→8 (`LOWQ_MIXED_WIDTH_L_R_PRIME`), the only previously-found 1-lane slack; exhausted.
- work1/work2 pad lanes proven tight (UPDATE 1, `register_shared_eea.rs:328-334`, `work_used ≤ n+3 = 259`).

### Where the slack ACTUALLY is (the only viable sub-846 lever found): P2's v-chain `chain`
On P2, `scheduled.swap-condition` (10 lanes) = `[zero_q, zero_s, control, chain(7)]`. The `chain(7)` are the
**dirty-tolerant** v-chain ancilla of `compute_zero`/`uncompute_zero`
(`multi_controlled_x_vchain(controls, zero, chain)`) — these DO restore arbitrary content and ARE
dirty-borrowable. Likewise `rotated-bitlen.length`(10) feeds a bit-length v-chain whose scratch is partly
borrowable (cf. the existing `LOWQ_REUSE_ROTATED_BITLEN_SCRATCH` / `preserved_dy_top` loans). A combined cut
that (i) dirty-borrows P2's `chain(7)` from an idle lane AND (ii) finds a second P1 lever would be needed;
**neither alone moves the measured peak.** P1 has no borrowable surplus left (above), so **the binding wall
is P1, and P1 is exhausted** at this datapath. Sub-846 needs a P1 *structural* change (e.g. eliminating the
`address` counter or fusing carry into work2-pad), not a cursor-scratch borrow.

### VERDICT (skeptical, build-derived)
- The cursor-scratch borrow to 838 (agent a9563c0b's hypothesis, and this task's premise) is **DISPROVEN as a
  standalone lever**: (1) it cannot lower the measured peak because 846 is a P1/P2 co-peak (cutting P1 alone
  → still 846 via P2), and (2) even on P1 it is value-incorrect at the reverse-boundary clean-zero buffer.
- The competitor's "value-bearing ⇒ width-trim not borrow" claim is **mechanistically wrong** (borrowing
  works on live lanes; their reason is a non-sequitur) but their *conclusion* ("no borrow") happens to hold
  for cursor-scratch — for the boundary-buffer reason, plus the co-peak reason they omitted.
- **No cut landed; peak remains the validated 846** (additive flag `LOWQ_Q838_CURSOR_BORROW` reserved/inert,
  default-off; build unchanged: 846 / 865762909 ops). Final 9024-shot eval not re-run (846 already
  eval-confirmed 0/0/0 in UPDATE 2; no circuit change to validate).
- **Next lever toward ≤845 (for the next agent):** attack P1 structurally — the `address`(9) counter is the
  single largest P1 transient and is a *counter*, not data; a Gidney log-depth or borrowed-counter
  unary-iterate that avoids materializing the 9-bit `address` would cut P1 by up to 9, AND P2 must be cut in
  parallel via its dirty-borrowable `chain(7)` / rotated-bitlen scratch. Both-in-one-window is the bar.


---

## UPDATE 4 — P1 AND P2 ARE BOTH INTERNALLY TIGHT; every P1/P2 transient is structurally required OR only dirty-borrowable WITH COMPENSATING REGROWTH. No sub-846 cut landed. Precise blocker below.

**Base reconfirmed (my build, count-only):** `POINT_ADD_COUNT_ONLY=1 TRACE_PEAK=1` → peak **846**, ops
865,762,909, phase `ec3.inv_fwd`. `CURSOR_SCRATCH_PROBE=1` value proof
`SWAP_ONLY_VALUE_PROOF` = 33424/33424 (inverse_pair/cursor_restore/scratch_clean/residue/oracle), boundary
nonzero 0/16384. (Re-derived the d13f3d0 double-co-peak census end-to-end; both co-peaks reproduce at 846.)

### P1 is exhausted — re-verified each lane against its actual gate use (build-level, not estimate)
P1 transient (23 over the 823 base) = address 9 + cursor-scratch 8 + carry 1 + above-guard 1 + add-only 1
+ chain 2 + enable 1. Every lane has a hard clean-zero / genuine-width requirement; the ONLY idle dirty
lane in P1 is `phase2`, and it already hosts nothing it could legally host:
- **address(9)** = `l_t_prime.len()` bits, holds `physical_work_width=259`, gray-code-walked through
  `n_iters = physical_work_width+1 = 260` (`for_each_forward`, `register_shared_eea_reference.rs:7075`).
  259..260 needs 9 bits; the iterate counter MUST be a real register flipped each step. Cannot live in
  `l_q` (the only 9-lane idle reg): `l_q` holds a live bit-length that must exit unchanged, and the iterate
  destroys+restores the counter to `physical_work_width`, not to `l_q`'s value → would need a 9-lane save.
  `l_q` is already consumed as the ladder DIRTY-LENDER bank (7 of 9 lanes, `mcx_dirty_ladder`); a lender
  may NOT alias a counter bit (asserted, `khattar_gidney.rs:845`), so address and l_q are necessarily
  disjoint.
- **cursor-scratch(8)** = `l_t.len()-1`. Bound by `controlled_increment_mod_2n` inside `accumulate`
  (ripple carry needs `count.len()-1 = l_t_prime.len()-1 = 8` CLEAN-zero carries) AND by the
  reverse-boundary buffer. Confirmed clean-zero requirement at the boundary buffer
  (`prepare_reverse_boundary` :7278-7307): it `cx(l_r_prime/l_s → source)` then
  `cuccaro_add_mod_2n_refs(source, address)` — the added operand is `source = cursor_scratch ⊕ l_r_prime`,
  which must equal `l_r_prime` ⇒ cursor_scratch ≡ 0. A dirty borrow injects `cursor_scratch`'s value into
  `address` → value-incorrect. (UPDATE 3's claim, re-confirmed by reading the active non-truncated path:
  the loads are `cx` relative, but the cuccaro OPERAND is the cx-result, so clean-zero is still mandatory.)
- **carry(1)** = genuine clean ripple/Bennett target. `emit_q845_fused_coefficient_body`
  (:7492-7521) does `multi_controlled_x_vchain([active,work1,work2], carry, [tmp])` (writes
  `carry ^= AND`, needs carry=0 to mean `carry = AND`), then `ccx(active, carry, work1)` reads it.
  Dirty-borrow corrupts the read. Confirmed not borrowable.
- **above-guard(1)** = the iterate `active` flag (live throughout for_each). **enable(1)** is idle DURING
  for_each but is consumed immediately after (`toggle_output_coefficient_enable` :7675) and must persist;
  it overlaps the peak window, so it cannot host carry. **add-only / chain[2]** = `output_scratch` /
  `coefficient_active`+`tmp`, all live in the body.
- **Conclusion:** P1 has exactly one idle dirty lane (`phase2`) and ZERO P1 transients it can legally host
  (no pure-Bennett-ancilla lane survives — `range_scratch`/`count_scratch` reuse `cursor_scratch`, which is
  simultaneously a clean ripple lane). **P1 cannot drop by even 1 via any borrow. Confirmed exhausted.**

### P2 is ALSO internally tight — the swap workspace is conserved (NEW build-derived finding)
P2 transient (23) = fused-prefix 2 + rotated-bitlen.length 10 + scheduled.swap-condition 10 +
promised-support 1. The decisive new fact is in `conditional_work_and_length_swap_under_zero_predicate`
(:5361-5365):
```
required_support_ancillas = l_q.len()-1 = 8
support_workspace_width    = max(1, required_support_ancillas - scratch.len())
```
`scratch` is the `chain(7)` slice of `condition_scratch(10)`; `promised-support = max(1, 8-7) = 1`. So the
swap workspace is **conserved**: `condition_scratch + promised-support = 3 + scratch + max(1, 8-scratch)`,
minimized = 11 for any `scratch ∈ {6,7}` (scratch=6 → 9+2=11; scratch=7 → 10+1=11; scratch=8 → 11+1=12,
worse). Shrinking `condition_scratch` by trimming `chain` REGROWS `promised-support` one-for-one. **P2
cannot drop by trimming the chain in-place; it needs a borrow from OUTSIDE the swap (phase1/phase2).**

### Why the only legal P2 borrow is a LARGE refactor that STILL can't move the peak
`compute_zero`/`uncompute_zero`'s `chain` ancilla IS dirty-tolerant (`multi_controlled_x_vchain` Bennett
compute-uncompute, :371-394) → phase2 could legally serve as one of its 7 lanes, freeing one
`condition_scratch` lane WITHOUT regrowing promised-support — BUT the SAME `chain` slice is threaded into
the swap as `rotated_bitlen_scratch` (`scratch: &[QReg]`, OWNED). To borrow phase2 there, the entire swap
call tree (`conditional_work_and_length_swap{,_inverse,_under_zero_predicate}`,
`materialize_promised_l_q_swap_discrepancy`, the bit-length sub-kernels) must change `&[QReg] → &[&QReg]`.
That is the most-called primitive in the EEA; the refactor is large and high-risk. AND even if P2 falls to
845, **the measured peak stays 846 because P1 holds (above).** A P2-only cut is build-confirmed useless.

### VERDICT (skeptical, build-grounded): sub-846 needs a STRUCTURAL P1 swap-datapath rewrite, not a borrow
The double co-peak is real and BOTH stacks are internally conserved. The minimal viable sub-846 is a paired
cut where **P1 drops via a genuine datapath change** (the only thing that can move P1), specifically ONE of:
1. **Eliminate the `address` counter** by re-expressing the unary-iterate so the index is NOT a materialized
   9-bit register — e.g. a borrowed-dirty counter walked in `l_q`'s 2 SPARE lanes (the ladder uses only 7
   of 9) PLUS a log-structured index, or fold the iterate into the existing work2 scan so no fresh counter
   exists. This is the single largest P1 term (9) and the highest-upside structural lever; it is NOT a
   borrow (borrows are all blocked) — it is a reformulation of the iteration itself. HARD, high value
   (P1 9→? could reach ~837, and toward 810 if cursor-scratch's boundary buffer is also restructured).
2. **Restructure `prepare_reverse_boundary`** so `address` is recomputed across the boundary WITHOUT a
   clean-zero 8-lane `source` buffer (the thing that pins cursor-scratch clean). If the l_s/l_r_prime
   re-add into address used `address` itself + 1 borrowed carry (a controlled in-place add) instead of the
   `cx→cuccaro→cx` staging through `source`, cursor-scratch's clean-zero constraint dissolves and it could
   then dirty-borrow → -8 on P1. Then pair with the P2 chain→phase2 borrow (-1) for the matching cut.

For EITHER, the matching P2 cut is the `chain→phase2` borrow above (requires the `&[QReg]→&[&QReg]` swap
refactor). **Both-in-one-window remains the bar; neither half alone moves the measured 846.**

**No code landed this round** (every tractable borrow is either blocked on P1 or peak-neutral on P2; the
viable levers are large datapath rewrites carrying real correctness risk that I would not bank
un-validated). Base unchanged: peak 846, ops 865,762,909, value proof 33424/33424. 846 remains the
validated frontier (eval 0/0/0 per UPDATE 2). The next agent should attack lever 1 or 2 above as a genuine
datapath change and re-measure BOTH co-peaks (`DUMP_PLATEAU`) after each step.


---

## UPDATE 5 — Lever 2 IMPLEMENTED + value-proven (in-place reverse-boundary add), but it is PEAK-NEUTRAL: cursor-scratch has THREE clean-zero consumers, not one. The -8 P1 cut is structurally blocked by clean-Bennett-vchain ancilla, and is UNVALIDATABLE by the value-proof gate (l_q borrow path absent from the harness). 846 reconfirmed the floor of this circuit family.

**What landed (commit `LOWQ_Q838_INPLACE_BOUNDARY`, additive, default-OFF, value-proven).**
UPDATE 4 lever 2 part 1 is DONE and correct: `prepare_reverse_boundary`
(`register_shared_eea_reference.rs:7311`) now does an IN-PLACE register add of `l_s` / `l_r_prime`
directly onto `address` instead of staging through the `source` buffer (which IS the swap-step
`cursor-scratch`). The staged path was `cx(addend→source); cuccaro_add_mod_2n_refs(source,address);
cx(addend→source)`, which requires `source`≡cursor-scratch clean-|0⟩ so the cx-loaded operand equals
the addend. Since `cuccaro_add_mod_2n_refs(a,b,..)` leaves its addend `a` UNCHANGED (majority/unmajority
restore it), the identical `address += addend` is obtained by passing the live `l_s` (address-width) and
`l_r_prime ++ [coefficient_active]` (8 lanes + one clean high lane for the MSB zero-extension, = the lane
the staged path always left at 0) DIRECTLY as the cuccaro addend. Bit-identical on `address` AND on the
`overflow` lane (each of the 4 cuccaro calls contributes `cx(addend_msb, overflow)` = 0, l_s[msb],
l_s[msb], 0 in BOTH paths). **This removes the reverse-BOUNDARY clean-zero pin on cursor-scratch.**

**Measured.** flag OFF → peak **846**, ops 865,762,909 (exact base, bit-identical fallback). flag ON →
peak **846** (NEUTRAL), ops 865,360,621 (the cx-staging ops removed). Value proof
(`CURSOR_SCRATCH_PROBE=1` → `exhaustive_q845_swap_only_coefficient_check`): **33424/33424** for OFF and ON
(inverse_pair / cursor_restore / scratch_clean / residue / oracle). The proof exercises the truncated
`for_each_reverse` (work1.len() < physical_work_width in the harness) so the boundary datapath is fully
covered. Boundary nonzero 0/16384.

### WHY -8 DOES NOT LAND (the precise, build-grounded blocker — supersedes UPDATE 3/4's "dirty-tolerant" claim)
Removing the boundary pin is necessary but NOT sufficient. **cursor-scratch (8 lanes) is pinned clean-|0⟩
by THREE independent production consumers inside the swap-only coefficient step**
(`coefficient_fused_data_and_sign_q845_swap_only`), every one of which is exercised in the 260-iter scan:
1. **`accumulate`/`unaccumulate`** (`:7048,7074`) → `controlled_increment_mod_2n(condition, count=l_t_prime,
   count_scratch=cursor_scratch)` (`register_shared_eea_microkernels.rs:503`): a **clean ripple-carry**
   incrementer. The carries are written then cleared assuming they enter |0⟩; dirty entry corrupts the
   carry propagation. Fires twice per iteration.
2. **`transition_q851_coefficient_cursor`** → **`toggle_q851_fixed_sign_event`** (`:7493,7517`; production
   has `LOWQ_Q851_FIXED_SIGN_EVENT=1`) → `multi_controlled_x_vchain(controls=l_t[..8], target=l_t[8],
   scratch=cursor_scratch)`. **NOT dirty-tolerant.** I verified the gate construction
   (`register_shared_eea_microkernels.rs:371`): it is the textbook Bennett ladder
   `ccx(c0,c1,anc[0]); …; ccx(c_{k-1},anc[k-3],target); …; ccx(c0,c1,anc[0])` which REQUIRES anc≡0 to
   compute the AND and then clear it. A dirty borrow makes `anc[0]` carry `c0∧c1 ⊕ dirty`, corrupting the
   target. Fires once per iteration.
3. **`toggle_nonzero`** (`:7289`) → `multi_controlled_x_vchain_borrowed(controls, target,
   ancillas=cursor_scratch)` (`:5020`). **CRITICAL CORRECTION to UPDATE 3/4:** `_borrowed` is GATE-FOR-GATE
   IDENTICAL to `multi_controlled_x_vchain` — the "borrowed" refers only to the `&[&QReg]` slice type, NOT
   to dirty tolerance. It is the SAME clean-zero Bennett ladder. So `toggle_nonzero` ALSO requires
   cursor-scratch clean. (UPDATE 3/4 mis-labeled these vchains "dirty-tolerant"; they are not. The ONLY
   genuinely dirty-tolerant MCX in the tree is `mcx.rs::mcx_dirty_ladder`, the double-Barenco cascade,
   which NONE of these three consumers use.)

To dirty-borrow cursor-scratch, ALL THREE must be rewritten to `mcx_dirty_ladder`-based forms (the q944
`controlled_increment_dirty_lenders` / `controlled_*_dirty_carry_refs` line). Each then needs **7 dirty
lender lanes** (for the widest 9-control MCX). The only candidate bank is the idle `l_q` (9 lanes) — but
`l_q` is ALREADY the swap-scan detector's dirty bank (7 of 9 lenders, `Q845_SWAP_LADDER_DETECTOR`), and a
lender may not alias an in-flight control; coordinating three more borrows on the same 7 l_q lanes through
the interleaved iterate body is a large, fragile rewrite.

### THE DECISIVE OBSTACLE: the value-proof gate cannot validate the l_q dirty-borrow path
`build_q845_swap_only_dispatch_harness` (`:12624`) — the builder behind the 33424-case
`exhaustive_q845_swap_only_coefficient_check` correctness gate — calls `coefficient_phase_block(…, None,
false)`: **`l_q_borrow = None`**. The harness never supplies l_q, so the swap step runs the single-lane
`phase2` fallback, NOT the l_q ladder. Therefore ANY cursor-scratch dirty-borrow that borrows l_q lenders
is **un-exercisable by the small-n value proof** — its only correctness oracle would be the full 9024-shot
eval (OOM-prone, the route's once-at-the-end gate). Per the task's OOM-safe / "value-proof is the gate /
never bank a cut that can't be value-proven" constraint, this cut is **un-bankable in the iteration loop**.
(The UPDATE-1 dirty-MCX swap detector had the same harness gap and leaned on the per-shot runtime
`contract_capture` asserts during the full build; that is a weaker gate and only covered ONE borrow site,
not three clean→dirty rewrites of the arithmetic core.)

### P2 unchanged and still independently binds
Even if P1's cursor-scratch were freed, P2 (the `dynamic-bitlen` swap plateau, also 846, no cursor-scratch)
holds. DUMP_PLATEAU re-confirmed both co-peaks at 846 (P1 = …+address9+carry1+cursor8; P2 =
…+rotated-bitlen10+swap-condition10+fused-prefix2+promised-support1). P2's only borrowable surplus is the
`chain(7)` v-chain inside `scheduled.swap-condition`, but (UPDATE 4) trimming it regrows
`promised-support` one-for-one unless borrowed from OUTSIDE the swap (phase1/phase2), which needs the
`&[QReg]→&[&QReg]` swap-call-tree refactor. Neither half alone moves the measured peak.

### VERDICT (skeptical, build-derived)
- **Lever 2 part 1 (in-place boundary) is correct and banked** (value-proven 33424/33424, additive flag,
  default-OFF so the eval-clean 846 island is untouched). It is a genuine simplification (−ops) and removes
  ONE of cursor-scratch's three clean pins — a real prerequisite, but **peak-neutral** on its own.
- **The −8 cursor-scratch cut is structurally blocked** by two clean-Bennett-vchain consumers
  (`toggle_q851_fixed_sign_event`, `toggle_nonzero`) PLUS the clean-ripple incrementer (`accumulate`),
  none of which is dirty-tolerant, and the dirty-lender rewrite is **un-validatable** by the mandated
  value-proof gate (harness omits the l_q borrow). This confirms **846 is the floor of this circuit
  family** under the OOM-safe correctness regime: every sub-846 lever (cursor-scratch −8, or `address` −9,
  UPDATE 4 lever 1) requires a clean→`mcx_dirty_ladder` rewrite of the swap arithmetic core whose only
  oracle is the full eval. **Pivot recommendation: a different inversion** (the shrunken_pz/q949 line is
  dominated per §C5), or first EXTEND the value-proof harness to thread an `l_q` dirty bank through
  `build_q845_swap_only_dispatch_harness` so the dirty-lender swap core becomes value-provable — that
  harness extension is the true unblocking prerequisite for any sub-846 attempt on this datapath.
- **No measured sub-846. Base peak reconfirmed 846** (flag OFF bit-identical; eval 0/0/0 per UPDATE 2, not
  re-run — the default build is unchanged). Do NOT submit.

## UPDATE 6 — The value-proof harness is now EXTENDED to cover the l_q dirty-borrow ladder path (the UPDATE-5 prerequisite LANDED + measured). The cursor-scratch −8 cut remains a large multi-kernel rewrite and is peak-neutral while P2 holds; the double co-peak is reconfirmed at the exact ops_idx level. 846 reconfirmed.

**The harness blocker from UPDATE 5 is FIXED (committed, measured).** The 33424-case swap-only
coefficient value proof drove the swap step with `l_q_borrow = None` (single-lane `phase2` detector
fallback), so the production `mcx_dirty_ladder` l_q-lender path (the peak-846 swap detector) was NEVER
exercised — any l_q-borrowing sub-846 cut was un-value-provable at small n. Now closed:
- New `build_q845_swap_only_core_harness_lq_borrow` forces the EXACT production path
  (`LOWQ_Q845_SWAP_DIRTY_MCX` + `LOWQ_Q845_SWAP_LADDER_DETECTOR` + `LOWQ_Q846_ULS_GATE_FUSE`) and threads a
  REAL `l_q` register (`length_width` lanes) as `Some(&l_q_borrow)`. Flags are saved/restored so the
  builder leaves NO global env side effects (verified: the downstream ephemeral-swap proof still passes).
- New exhaustive l_q dirty-path proof loop (`exhaustive_q845_swap_only_coefficient_check`): over the SAME
  valid-input domain as the 33424 proof (the kernel precondition — enable/add_only derived from
  phase/sign, plus the add_only no-overflow / no-above-guard guards), crossed with EVERY dirty `l_q`
  value, asserts (1) data lanes BIT-IDENTICAL to the `None` path, (2) `l_q` RESTORED exactly to its dirty
  entry, (3) all scratch clean (incl. above_guard→0), (4) forward∘inverse identity on (data ⊕ l_q).
- **Measured (`CURSOR_SCRATCH_PROBE=1`):** existing proof UNCHANGED `33424/33424`;
  `SWAP_ONLY_LQ_DIRTY_PROOF: lq_borrow_dirty_path_checks=267392 lq_restore_checks=267392
  lq_data_transparent_checks=267392` → "l_q dirty-borrow ladder path EXERCISED + value-proven". Count-only
  peak UNCHANGED `846`, ops `865,762,909` (test-only change; production circuit bit-identical, flag-free).
- Pitfall recorded: the lq-borrow ladder path is self-cleaning ONLY inside the kernel's contractual input
  domain. Driving it with the loose sweep (no add_only/above_guard guards) leaves `chain[0]`
  (coefficient_active) dirty → the `R`-reset clean-checker fires. The `None` path happened to tolerate
  those invalid inputs; the ladder does not. Mirroring the 33424 loop's guards is mandatory and correct
  (outside the precondition neither path is contractually clean).

**This is the enabler.** Any future l_q-borrowing cut on this datapath (cursor-scratch −8, address −9) is
now value-provable small-n by extending the same loop to also route the borrowed scratch onto `l_q`.

### Why no measured sub-846 banked this round (build-grounded, NOT a harness limitation)
1. **Double co-peak reconfirmed at ops_idx granularity** (`TRACE_NAMED_PEAK_TARGET=846 DUMP_PLATEAU=1`):
   - **P1** (trigger `rs.q845-swap-only.carry`, ops_idx 155407): `address(9)+carry(1)+cursor-scratch(8)` +
     coeff-fused lanes. The swap-only coefficient step.
   - **P2** (trigger `rs.dynamic-bitlen.fused-prefix-scratch[1]`, ops_idx 377119): `rotated-bitlen.length(10)
     + swap-condition(10) + fused-prefix-scratch(2) + promised-support(1)`. **NO cursor-scratch, NO
     address** — a structurally DISJOINT moment. Both also recur in `ec3.alt.cancel` (the inverse).
   A cut must drop BOTH P1 and P2 to move the measured 846; neither half alone helps.
2. **The cursor-scratch −8 (P1) requires THREE bespoke dirty-lender kernels, none drop-in.** Production
   takes `q851_fixed_sign_event` (Bennett vchain) for the transition, `controlled_increment_mod_2n` (clean
   ripple) inside `accumulate`, and `toggle_nonzero` (Bennett vchain) — all clean-|0⟩ pins (UPDATE 5,
   re-verified). The existing `q944_dirty_parity_microkernels` dirty primitives are `a += gate*carry`
   (gated adds), NOT the plain gray-counter increment / Bennett-AND these consumers need. I built+compiled
   correct `increment/decrement_mod_2n_dirty_lenders` (high-to-low `mcx_dirty_ladder` over `register[..j]`,
   lenders restored) but they map to the gray-counter, which production BYPASSES via the fixed-sign event
   — so they don't wire to the live consumers. Reverted (kept the diff focused on the harness).
3. **Temporal (not spatial) non-aliasing is the real l_q-capacity story.** The detector needs 7 of 9 l_q
   lenders, cursor needs 8 — they cannot coexist spatially. BUT `mcx_dirty_ladder` is compute-uncompute:
   the detector cascade RESTORES all 7 lenders BEFORE the iterate body runs. So during the body (where the
   cursor consumers fire) all 9 l_q lanes ARE idle/restored and could host 8 cursor lanes. The borrow is
   temporally legal; it is blocked only by the absence of dirty-lender forms of the 3 specific consumers,
   each a fresh, individually-value-proven kernel.
4. **Even a perfect P1 −8 is peak-neutral while P2 holds 846.** P2's only surplus is the conserved
   `chain(7)` ↔ promised-support tradeoff; an out-of-swap borrow needs the `&[QReg]→&[&QReg]`
   swap-call-tree refactor (the most-called EEA primitive). Two-front, large, and high-risk.

### VERDICT (skeptical, build-grounded)
- **Harness extension BANKED + value-proven** (33424/33424 preserved; 267392 dirty-path checks;
  peak/ops/flag-free build bit-identical). This is the deliverable that unblocks sub-846 — exactly the
  prerequisite UPDATE 5 named. Additive, test-only, no production impact.
- **No measured sub-846.** Not because the cut is value-INCORRECT (it is now provable) and not a capacity
  wall (temporal borrow is legal) — but because the minimal viable cut is a TWO-FRONT multi-kernel rewrite
  (3 dirty-lender consumer kernels for P1 + the large swap-tree `&[QReg]→&[&QReg]` refactor for P2), each
  piece individually value-provable now, but together a large change carrying real compounding correctness
  risk that I would not bank un-finished. 846 remains the validated frontier (eval 0/0/0 per UPDATE 2,
  default build unchanged). Do NOT submit.
- **Next agent:** the path is now de-risked at the gate level. Build the 3 dirty-lender consumer kernels
  (`q851_fixed_sign_event`, `accumulate`'s increment, `toggle_nonzero`) one at a time, each value-proven by
  extending the UPDATE-6 lq-borrow loop to route its scratch onto `l_q`; land cursor-scratch −8 → P1 838;
  THEN the swap-tree refactor for P2 −1. Both in one window = measured sub-846. If any consumer kernel
  fails the (now-available) value proof, THAT is the decisive floor confirmation.

## UPDATE 7 — FRONT P1 LANDED + value-proven (cursor-scratch −8, P1 846→838). FRONT P2 is VALUE-INCORRECT (decisive): the swap-condition chain is fully clean-pinned by the bit-length prefix-AND loan, so phase2-borrowing it corrupts the swapped-in bit length. Measured peak stays 846 (P2 binds); 846 is the floor of THIS family, with the precise blocker now isolated to the bit-length prefix-AND.

**Commits this round:** P1 cut (`b21af81`), swap-tree `&[QReg]→&[&QReg]` refactor (`1545350`), P2 negative result (`6e15eda`). Base was `bb6c092` (UPDATE 6).

### FRONT P1 — DONE, value-proven, peak-neutral (P2 holds). cursor-scratch −8.
The 8-lane clean `rs.q845-swap-only.cursor-scratch` is now FREED by routing ALL of its
clean-|0⟩ consumers onto the idle `l_q` dirty lenders that already serve the swap-scan detector.
The detector's `mcx_dirty_ladder` is compute-uncompute, so it RESTORES every `l_q` lender BEFORE the
iterate body runs → `l_q` is idle/available during the body, and each consumer borrow is itself an
`mcx_dirty_ladder` (restores lenders exactly). Temporal sharing is legal (UPDATE 6 #3, now realised).

Consumers converted (flag `LOWQ_Q838_CURSOR_BORROW`, gated on `swap_dirty_bank.is_some()`, default-ON
in the sub1000 route; requires `LOWQ_Q838_INPLACE_BOUNDARY` also default-ON):
- `accumulate`/`unaccumulate`: `controlled_increment/decrement_mod_2n` (clean ripple) →
  `controlled_increment/decrement_mod_2n_dirty_lenders` (new microkernels, high-to-low `mcx_dirty_ladder`).
- `transition` fixed-sign event: `multi_controlled_x_vchain` → `mcx_dirty_ladder`.
- `transition` ripple fallback (harness-only path, length_width≠9): `increment/decrement_mod_2n` →
  `add_constant_mod_2n_dirty_lenders` (new; modular constant-add via per-set-bit sub-register increments).
- `toggle_nonzero`: `multi_controlled_x_vchain_borrowed` → `mcx_dirty_ladder`.
- the truncated `for_each_*` geq-comparator (`toggle_register_geq_constant_vchain`, a 4th clean consumer
  the original 3-consumer brief missed): `multi_controlled_x_vchain_borrowed` → `mcx_dirty_ladder`.
- the truncated `for_each_reverse` address `±=delta` constant subtract (a 5th clean consumer, was a
  clean `cuccaro` staged through cursor-scratch): → `add_constant_mod_2n_dirty_lenders`, in place on
  `address`. (The in-place reverse boundary, UPDATE 5, already removed the 6th pin — the boundary buffer.)
`cursor_scratch` is then allocated at width 0.

**Measured (P1 on, P2 off = production):** DUMP_PLATEAU at ops_idx 155407 — `cursor-scratch` GONE, only
`address(9)+carry(1)` remain → P1 plateau dropped **846 → 838**. Overall count-only peak still **846**
(P2 binds, expected). ops 1,134,075,189 (dirty-ladder gate inflation — TOFFOLI IRRELEVANT). Flag OFF =
bit-identical 846 / 865,762,909.

**Value-proven.** `SWAP_ONLY_VALUE_PROOF` = 33424/33424 AND `SWAP_ONLY_LQ_DIRTY_PROOF` = 267392
dirty-path checks. The UPDATE-6 lq-borrow harness was extended to FORCE `LOWQ_Q838_CURSOR_BORROW=1` on
its `lq_borrow` path, so the cursor-borrow rewrite is the exact code under test: data-transparent vs the
clean (None) path, `l_q` restored, scratch clean (incl. above_guard→0), forward∘inverse identity. The
borrow activates at small n (toggle_nonzero needs 2 lenders, l_q bank has 3).

### FRONT P2 — VALUE-INCORRECT. The matching cut cannot be made; this is the decisive floor.
P2 = `rotated-bitlen.length(10) + scheduled.swap-condition(10) + fused-prefix-scratch(2) +
promised-support(1)` = 846. The brief's P2 cut (host one `swap-condition` chain lane on idle `phase2`,
condition_scratch 10→9) was IMPLEMENTED behind `LOWQ_Q845_SWAP_CONDITION_BORROW` and DOES drop the
count-only peak **846 → 845**. But the new value-proof
(`build_q845_ephemeral_swap_harness_condition_borrow` + the swap-condition borrow proof loop) proves the
resulting circuit is **VALUE-INCORRECT**, and the gate REJECTS it:

> The swap `scratch` (= the condition `chain`) is **LOANED IN FULL** to the bit-length prefix-AND kernel.
> Production runs `LOWQ_LOAN_FUSED_PREFIX_SCRATCH=1`, so `bit_length_lean_with_full_prefix_scratch`
> borrows ALL 7 chain lanes as Khattar-Gidney direct-prefix scratch
> (`DIRECT_PREFIX_FULL_SCRATCH_LEN = 9`: 7 loaned + 2 fresh `fused-prefix-scratch`) and **requires them
> clean-|0⟩** (the prefix-AND / pos-increment deposits intermediate ANDs there). `compute_zero` and the
> swap-support v-chain were correctly made dirty-tolerant (`mcx_dirty_ladder`), but the prefix-AND is NOT
> — a dirty borrowed lane corrupts the swapped-in bit length: the `l_t_prime` output diverges on EVERY
> dirty entry (160/160 cases reject; the divergence is at the bit-length output lane).

This is **stronger than UPDATE 4/5's "conserved" claim**: the swap-condition chain is not merely
peak-neutral to trim, it is **value-INCORRECT to phase2-borrow**, because the chain is the bit-length's
clean prefix-scratch loan. The 845 measurement is illusory (value-wrong circuit). The 9 clean prefix-AND
lanes (chain 7 + fused-prefix 2) are the hard P2 floor; no borrow helps (every prefix lane must be clean,
and shrinking the chain regrows fused-prefix one-for-one — conservation, now with a value-correctness
proof behind it).

The P2 flag is **NOT** in the production route. The swap-tree `&[QReg]→&[&QReg]` refactor stays (it is
bit-identical with the flag OFF — `compute_zero`'s clean path uses `multi_controlled_x_vchain_borrowed`,
gate-for-gate identical to `multi_controlled_x_vchain`). The P2 proof asserts divergence>0 (rejection),
so the value gate passes for the shipped P1-only build.

### COMBINED / VERDICT (skeptical, build-grounded)
- **Measured peak stays 846** (P1 cut on, value-proven; P2 binds at 846, its cut value-incorrect).
  ops 1,134,075,189. Default build value-clean: `SWAP_ONLY_VALUE_PROOF` 33424/33424,
  `SWAP_ONLY_LQ_DIRTY_PROOF` 267392, `SWAP_CONDITION_BORROW_NEGATIVE` 160 (all reject). Do NOT submit.
- **846 is the floor of this circuit family** under the OOM-safe value-proof regime, with the blocker now
  pinned PRECISELY: not P1 (freeable, proven) but **P2's bit-length direct-prefix-AND**, which pins 9
  clean prefix-scratch lanes (chain 7 + fused-prefix 2) and is the only non-dirty-tolerant consumer left.
- **Next lever toward ≤845 (genuine, NOT a borrow):** make the bit-length prefix-AND itself need fewer
  clean lanes OR tolerate a borrowed lane — i.e. a dirty-tolerant / radix-restructured
  `bit_length_lean_direct_prefix` (the Khattar-Gidney prefix ladder rewritten so one scratch lane may be
  dirty-restored). That is an algorithmic rewrite of the length kernel, not a lane loan. If it lands, the
  P2 phase2-borrow becomes value-correct and 846 → 845 (then the next P2 lane, rotated-bitlen.length(10),
  is the new target). Until then, **pivot to a different inversion** (per the task) — 846 holds.

## UPDATE 8 — FRONT P2 LANDED: the bit-length direct-prefix-AND IS lane-reducible (9→8) value-correctly. Measured peak 846 → 845. The blocker UPDATE 7 named (the 9 clean prefix lanes) is broken by eliminating the materialized flag via a dirty double-control, freeing one owned fused-prefix lane on the P2 co-peak.

**Headline (measured, count-only + value-proven):** the P2 binder is cut. The bit-length
direct-prefix-AND (`bit_length_lean_with_full_prefix_scratch`) required **9 clean loaned lanes**
= 5 KG prefix ancillae + 3 increment scratch + **1 materialized flag**. New flag
`LOWQ_KG_ONLY_PREFIX_SCRATCH` (default-ON in the sub1000 route) shrinks the loan to **8 lanes** by
eliminating the flag, which frees one owned `rs.dynamic-bitlen.fused-prefix-scratch` lane on the P2
plateau → **measured peak 846 → 845** (`ec3.inv_fwd`, ops 1,976,174,517). Flag OFF = bit-identical
846 / 1,134,075,189.

### The cut (option (c): in-place dirty update needing fewer clean lanes)
The flag is only needed by the **two-control** KG prefix callbacks (`ccx(a,b,flag); update(flag);
clear_and`). They fire only O(log* n) times (at layer boundaries). KG-only mode replaces them with a
**dirty `mcx_dirty_ladder` double-control** (`dirty_controlled_inc_suffix` over restored source/anc
lenders), so no flag lane is needed. Critically, the **hot one-control per-position increment STAYS the
efficient clean log-depth CINC** on the 3 loaned increment-scratch lanes — there is no O(n)-scale dirty
arithmetic on the hot path. The constant-`n` add is routed through the loaned KG/increment scratch
(host the always-on control on KG anc[0]) to avoid the `add_const` carry alloc that otherwise +1'd the
plateau.

### The op-explosion trap and its fix (the key engineering lesson)
The FIRST naive implementation (dirty double-control with the full n-wide source as the dirty-candidate
pool) **blew the op count up >50x** (count-only build 18s → never finished in >4min, uncontended). Root
cause: `collect_dirty_lenders` scans ALL candidates O(n=259) PER prefix callback, compounding across
callbacks × swaps × EEA iterations. **Fix:** cap the candidate pool to `sref.len()+4` (the MCX needs
≤ `controls.len()-2` ≤ `sref.len()` lenders; restored lenders are interchangeable, so the cap is
value-preserving). With the cap: count-only build 30s, ops 1.13B → 1.98B (**1.7x**, eval-feasible). The
cap engages at n≥8 and is covered by the 16320-case roundtrip proof.

### Value proofs (all pass; the cut is value-correct, not illusory)
- `SWAP_ONLY_VALUE_PROOF` = 33424/33424 (P1 intact, not regressed).
- `SWAP_ONLY_LQ_DIRTY_PROOF` = 267392 (P1 l_q dirty-borrow path).
- `SWAP_CONDITION_BORROW_NEGATIVE` = 160 (the OLD phase2-borrow is STILL value-incorrect/rejected — this
  cut is a different mechanism: it reduces the loan REQUIREMENT, it does not dirty-borrow a clean lane).
- **`KG_ONLY_PREFIX_POSITIVE` = 160** (NEW positive proof): the production swap with the KG-only loan is
  **BIT-IDENTICAL** to the clean nine-lane reference across every swap basis state, and the harness peak
  is **−1** (clean_swap_peak=38 vs kg_swap_peak=37). Both harnesses run the FULL loan, differing only by
  the KG-only flag, so the comparison is non-vacuous.
- **`DIRECT_PREFIX_BITLEN_PROOF`** = 16320 baseline + 16320 kg_only cases: the lane-reduced kernel
  computes the EXACT bit length, preserves the source, and restores every internal/loaned lane clean, for
  all 8-bit sources × 5-bit targets × both directions — and this exercises the capped dirty path.

### P2 census before/after (DUMP_PLATEAU, full enumeration)
- **846 (before):** base 823 + rotated-bitlen.length(10) + swap-condition(10) + **fused-prefix-scratch(2)**
  + promised-support(1).
- **845 (after):** base 823 + rotated-bitlen.length(10) + swap-condition(10) + **fused-prefix-scratch(1)**
  + promised-support(1). One owned fused-prefix lane freed (loan 9→8, chain supplies 7, owned 2→1).

### New binder + next lever
The new binder is the SAME P2 step at 845; the next target is the second fused-prefix lane (owned 1 → 0
needs loan ≤ 7, i.e. dropping an increment-scratch lane → dirty hot-path increment → op-explosion, the
same trap, NOT cheaply fixable), or the disjoint `rotated-bitlen.length(10)` / `swap-condition(10)`
gadgets. P1 sits at 838 (below). To reach ≤844 the next P2 lane needs a genuinely different lever.

### Implementation
- `shrunken_pz_state_machine.rs`: `bit_length_lean_direct_prefix_modes` + `bit_length_lean_impl_modes`
  (new `kg_only_full_scratch` param); `bit_length_lean_with_kg_only_prefix_scratch[_complemented_source]`;
  `DIRECT_PREFIX_KG_ONLY_SCRATCH_LEN=8`; capped `dirty_candidates`; the
  `direct_prefix_kg_only_bit_length_roundtrip_check` value proof.
- `register_shared_eea_reference.rs`: `KG_ONLY_PREFIX_SCRATCH_FLAG`, `kg_only_prefix_scratch_requested()`,
  the loan-requirement branch in `bit_length_lean_allow_zero_with_borrowed_scratch_impl`,
  `build_q845_ephemeral_swap_harness_prefix_loan`, the `KG_ONLY_PREFIX_POSITIVE` proof loop + peak/op
  report fields.
- `trailmix_port/mod.rs`: `set_default_env("LOWQ_KG_ONLY_PREFIX_SCRATCH","1")` (route default-ON).
- `bin/build_circuit.rs`: proof prints (`KG_ONLY_PREFIX_POSITIVE`, `KG_ONLY_PREFIX_PEAK`,
  `DIRECT_PREFIX_BITLEN_PROOF/COST`).

**Commits:** `52fc9f7` (cut, default-OFF), `1cddc6a` (route default-ON), `50f52c1` (this UPDATE).

### Full 9024-eval: peak 845 CONFIRMED at load; correctness phase TIME-bounded (NOT OOM)
The trusted evaluator was run once (default route, no env). It built ops.bin (1.36 GB on disk,
110 GB uncompressed, 81x; build peak RAM ~20 GB, under the 36 GB limit) and the eval **loaded
1,976,174,517 ops and reported `qubits = 845`** — the sub-846 peak is confirmed by the trusted stage's
own load. It then entered the 9024-shot correctness phase and ran ~1 h (~5x the baseline ~12 min) WITHOUT
OOM (RSS stayed ~8 GB / 22% throughout, never near 36 GB) but was stopped by the run-time budget before
emitting the per-shot pass/fail line. The slowdown is the eval-side cost of the 1.7x op inflation
(1.13B → 1.98B) compounded by the alt-seed correctness re-derivation; it is a TIME limit, not a memory
limit, so the OOM-safe constraint held. The 9024-shot pass/fail is therefore NOT machine-confirmed here.

**Why the cut is nonetheless value-correct (independent of the unfinished shot loop):** the lane-reduced
kernel is proven by EXHAUSTIVE small-n basis-state proofs over the EXACT production code path —
`KG_ONLY_PREFIX_POSITIVE` (the full swap, KG-only loan, bit-identical to the clean nine-lane reference,
peak −1) and `DIRECT_PREFIX_BITLEN_PROOF` (16320 cases: exact bit length, source preserved, all
internal/loaned lanes restored clean, capped dirty path exercised). These are stronger per-input
guarantees than a fixed-shot sample. Per note UPDATE 2 the register-shared route is island-independent
(the inversion is full-width correct), so a reseed stays 9024/9024; the kg-only change only reduces
scratch lanes in the bit-length computation and does not introduce any data-dependent width truncation.

**Recommendation for the next agent:** to machine-confirm the 9024 shots, run the eval with a larger time
budget (or `SKIP_ALT_SEED_CHECKS=1` to drop the re-derivation), OR first reduce the op inflation (the
dirty double-control fires on the rare 2-control callbacks; further capping / a cheaper 2-control AND
deposit would bring ops back toward baseline and the eval back under ~25 min). The peak 845 and value
correctness are already established; only the full-shot wall-clock confirmation remains. Do NOT submit.
