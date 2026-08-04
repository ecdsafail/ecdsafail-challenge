# Syntactic certification — three closed approaches, three passing controls

Sections 1 and 3 are measured on the certified frontier `cf5aa02` — the promoted source of submission
`0c5b1b7b-561a-48a0-abc6-5fefaffdc0ad`, score **`1,490,805,286`**. Section 2 was measured on `6909d15`
(score `1,486,468,554`) and has **not** been re-measured here; it says so where its numbers appear. Every
conclusion below was checked on both heads and none differs — only counts.

The §1 instrument independently reproduces four entries of `06-research-status.md`'s certified table, which is
what licenses the rest:

| quantity | `06-research-status.md` | [`repro/hotness.rs`](repro/hotness.rs) on `cf5aa02` |
|---|---:|---:|
| average executed Toffoli | `1,291,859.302` | `1,291,859.302` |
| total executed Toffoli | `11,657,738,337` | `11,657,738,337` |
| emitted operations | `9,062,420` | `9,062,420` |
| qubits | `1,154` | `1,154` |

Status vocabulary is `06-research-status.md`'s. All three results below are **Established** within their stated
class: each is a structural argument or an exact enumeration, not a search that ran out of time.

## The target

**46,286 CCX/CCZ — 3.35% of the score, 43,237.5 avgT — fire on none of the 9,024 official shots.** The
strict-beat bar is **0.802 avgT** (the frontier is `1,291,859 × 1154`, so `round(avgT) <= 1,291,858` scores
`1,490,804,132` and wins by 1,154 points), so *one* certified gate is a submission.

The obstacle is `04-traps.md` §1, in its sharpest form: the 9,024 test inputs are a SHAKE256 hash of the whole op
stream (`eval_circuit.rs:204`). Deleting a gate re-rolls every shot. "Never fired on this draw" is a certificate
that cannot survive its own use. Only a stream-agnostic certificate is shippable.

Three ways to get one were tried. All three are closed, and they close for the same reason.

## Why each negative carries a control

A search that finds nothing is indistinguishable from a search that cannot see. Every negative below is paired
with a known-answer or planted-signal control. **Two of the three controls fired before the negative was
believed**, which is the only reason the negatives are stated as results rather than as absences.

---

## 1. Cooling — making gates conditionally cold — CLOSED, class empty

The scorer bills a Toffoli by `popcount` of its condition mask and ignores the control qubits (`src/sim.rs:77-86`).
So a gate charged on fewer shots costs less **whether or not it ever fires**. Tighten condition stacks, harvest the
difference. The mechanism is real and correctly stated in `CEILING.md`.

**Control (known-answer).** [`repro/hotness.rs`](repro/hotness.rs) replays the official measurement — `build()`,
`analyze_ops`, the Fiat–Shamir XOF over the whole stream, the same 9,024-shot draw, the same 64-shot batching —
with `apply_iter` copied dispatch-for-dispatch and one line added to attribute `executed_shots` to the op that
incurred it. It asserts reconstruction of the scorer's own total:

```
GATE ok: attributed=11657738337 == toffoli_gates=11657738337
avgT=1291859.302        (06-research-status.md certifies 1,291,859.302)
```

**This control fired.** A first run drew 9,216 shots instead of 9,024 and gave `avgT=1288114.585` — plausible to
the eye, and wrong. Nothing downstream is meaningful without the `GATE ok` line.

**Result.** 1,347,438 CCX/CCZ, 9,024 shots, total charge 11,657,738,337.

| band | gates | share of gates | share of charge |
|---|---:|---:|---:|
| cold (hotness 0) | **0** | 0.000% | 0.000% |
| partial, all in `[0.4798, 0.5186]` | 111,162 | 8.250% | 4.303% |
| full (hotness exactly 1.0) | 1,236,276 | 91.750% | 95.697% |

**No gate is charged on fewer than 4,330 of 9,024 shots.** The hotness distribution is bimodal at 1.0 and 0.5 with
nothing anywhere else — the decile histogram is empty outside buckets 4, 5 and 10. The partial band's mean is
0.5000 and its ±0.02 spread is sampling noise on 9,024 draws.

The 0.5 band is a fair coin, and that is structural. **Every CCX/CCZ in the stream has `c_condition = NO_BIT`** —
one distinct condition-bit value across all 1.35 M gates. All conditioning therefore comes from the enclosing
`PushCondition` blocks, and the bits those push are `Hmr` measurement outcomes, which are independent of the
quantum controls. Firing is a function of the controls; the condition is a classical bit that cannot see them.
The probability that a given gate's ~2,256 firing shots all land inside a given fair coin's 4,512 true shots is
**2⁻²²⁵⁶ ≈ 10⁻⁶⁷⁹**. Manufacturing the bit instead means measuring the controls — an `Hmr`, which destroys the
qubit — to save at most 0.75 of one Toffoli.

Pricing the lever with a perfect oracle condition granted to all 1,347,438 gates at once: `charge − fire` summed
over the stream is **8,946,656,859 shot-charges = 991,429 avgT = 76.74% of the score** (only 23.26% of charge is
on a shot where the gate fires). Enormous, and unreachable for the reason above. The fire rate concentrates at
**0.25** — exactly the rate at which a Toffoli with two independent uniform controls fires — with **60.2%** of
unconditional gates in `[0.20, 0.30)`. Those gates are not wasteful. Three quarters of a Toffoli's charge is the
price of the shots where its controls are not both 1.

**The candidate class is empty for a structural reason, not a failed search.** Do not re-mine it.

---

## 2. Per-shot census sampling — CLOSED, cannot certify the target

Sample inputs, retain the per-shot effect mask per gate, and certify dead any gate that never fires. This is what
`census.rs` and the shipped `deep_strip_keys.rs` do.

**Control (known-answer).** Replaying the occupancy tripwire (`04-traps.md` §2) against a per-gate dump reproduces
`build_circuit`'s own accept/discard counts, which is what makes the keying — as opposed to the predicates —
provably correct.

**This control also fired**, in the useful direction: it isolated the discrepancy to the certification predicates
rather than to the addressing. On the table it was run against, of the shipped keys the tripwire accepts —
applied in a circuit that passes `9,024/9,024` — the census claims **3,076 dead keys fire (25.02%)** and
**1,674 downgrades violated (42.67%)**. The shipped table demonstrably yields 0/0/0. The census over-observes.

> **Provenance — this section was measured on `6909d15`, not on `cf5aa02`.** Its counts (dead 12,292 accepted /
> 251 stale, downgrade 3,923 / 0) are from a re-mined table of 16,466 keys, and are a third table again: `cf5aa02`
> emits **13,831 keys** and reports
> `[deep-strip-identity] removed 10743 / 10743 dead; downgraded 3088 / 3088 to CX/CZ; 0 stale keys skipped`,
> while `6909d15` emits 17,278. Neither frontier build reports a single stale key. The 25.02% / 42.67% figures
> therefore do **not** reproduce from `cf5aa02` without re-running the miner, and should be read as the magnitude
> of a measured disagreement, not as a property of the shipped table. The mechanism below does not depend on
> which table is used.

**Result.** The approach cannot certify any of the 46,286, and the reason is definitional rather than
quantitative. A sampler observes *firing*. It has no access to *why* a gate does not fire. If a gate is quiet
because of a data invariant — a theorem about what the registers can hold — the census cannot represent that and
can only report an empirical rate at whatever depth it ran. Two censuses at different depths, or over different
input distributions, then disagree about exactly the rarest-firing gates, which is the population in dispute.
"Never fired on this draw" is a statement about one draw, and any edit re-rolls all 9,024 inputs.

This is the same monotonicity `03-proven-floors.md` records under *Dead-gate census — OVER-drawn, not dry*: at
1e9 inputs only 1,290 of 1,442 shipped keys were still never-firing and 153 shipped keys fire. That section states
the depth-dependence as an observation; this section states its mechanism. **A 25%/43% disagreement between two
censuses is not a miner bug and not a defect in the shipped table — it is the signature of certifying by
observation something that is true by invariant.**

> **Scope.** `06-research-status.md`'s open problems are exact-eight joint synthesis, the controlled-addition
> factor-two gap, the streaming dialog ranker, new source implications, and low-qubit representation pricing. It
> does not list a census over-observation gap, and the shipped table works. The overlap is that its *Deep-strip
> localization* row concludes transfer failures originate upstream rather than in the final table — consistent
> with this, and a different claim. The gap resolved here is this fork's, carried in its own census-miner notes.

---

## 3. Affine-relation analysis over GF(2) — CLOSED, no relation exists

A `CCX` never fires if `q(c1) & q(c2) = 0` identically. Two syntactic cases decide that without knowing values: a
control is constant 0, or the controls are complementary (`c1 = ¬c2`). Complementary flag pairs are what binary-GCD
sign and branch logic ought to produce.

Two rungs. `constzero.rs` is a three-valued constant lattice. `affine.rs` carries a full affine form per qubit
(`constant XOR (XOR of atoms)`, atoms XOR-hashed into a `u128`), with `X`/`CX`/`Swap`/`R`/`Hmr` propagating exactly
under provably-`AllOnes` conditions, and `CCX` targets taking a **hash-consed AND term** keyed on the control forms
— an XOR-of-AND graph in which identical subexpressions collapse rather than degrading to opaque unknowns.

**Control (planted signal).** `affine.rs --positive-control` builds a stream where `CX(q0→q1); X(q1)` makes
`q1 = ¬q0`, places a `CCX(q1,q0,q2)` and a `CCZ(q1,q0,q5)` on that complementary pair, and a `CCX(q3,q0,q4)` on
unrelated controls:

```
CERTIFIED never-firing gates : 2
  op 2 kind 13 reason c1=!c2
  op 4 kind 14 reason c1=!c2
POSITIVE CONTROL PASS
```

Both planted pairs detected, the unrelated gate correctly not certified. `constzero.rs` carries matching
non-vacuity evidence: **6,298,889** ops tracked under provably-`AllOnes` conditions (against 2,757,267 Mixed), and
a Zero population moving 67 → 54 → 341 → 598 across four checkpoints as ancillas are allocated and uncomputed.

**Result: zero gates certified, by either rung.**

| diagnostic | value |
|---|---:|
| CCX total | 1,342,695 |
| atoms minted | 2,999,976 |
| CCX with a constant-1 control | **0** |
| CCX with equal controls | **0** |
| distinct AND terms | 1,229,522 |
| — recovered by hash-consing | 6,754 |
| CCX/CCZ whose controls share a single atom | **0** |
| certified | **0** |

**Not one gate in the circuit has controls that are affinely related at all** — not equal, not complementary, not
even sharing one atom. Hash-consing collapses 6,754 of 1,342,695 terms; the remaining 1.23 M nonlinear terms are
genuinely distinct subexpressions, so the zero is a property of the circuit and not an analyzer that gave up. The
constant rung agrees from the other side: all 1,347,438 CCX/CCZ have **both** controls Unknown, and not one has a
provably-One control.

Both rungs `--check` clean against the `cf5aa02` hotness dump — 0 certified, 0 violations — leaving
**46,286 never-firing gates, 0 explained structurally, 46,286 still unexplained.**

The constant rung's zero has its own explanation: `build()` already runs CONSTPROP (`dropped=282, folded_cx=23,
aff_drop=9`), so that certificate class was emptied before the stream was final.

**The complementary-flag intuition is not borne out.**

---

## The unifying finding

All three approaches reason about the **form** of a value.

| approach | what it inspects | why it fails here |
|---|---|---|
| cooling | which classical bit gates the charge | firing is quantum; the classical bits are independent coins |
| census | the empirical fire rate at depth N | an invariant is not visible in a sample at any depth |
| affine | the algebraic expression for each wire | 1.23 M distinct nonlinear terms, no two controls related |

**This circuit computes modular inversion and modular multiplication.** Nearly every value is a nonlinear function
of the inputs, so there is no exploitable form to inspect. Affine structure survives only through `CX`/`X` chains
uninterrupted by any `CCX`, and no such chain reaches any gate's control pair. The 46,286 gates are quiet because
of **what their controls can be**, not because of how those controls are written or how often they were watched.

One sentence covering three separate negatives, and it is the substantive result: the cheap certification routes
are not merely unproductive, they are **the wrong kind of argument**. A sampler cannot see an invariant. The
46,286 are not low-hanging.

## What would certify one

A **semantic** argument over the data invariants of the binary-GCD engine — register bit ranges, mutual exclusion
of branch flags, the loop invariant relating `u`, `v` and the schedule width. The tractable shape:

1. encode **one divstep** as a transition relation (the engine is a loop, so one step is small);
2. discharge the candidate invariant over that step with a bounded model checker or an SMT solver;
3. lift by induction over the 261 divsteps;
4. a gate whose control pair the invariant excludes is then certified for **all** inputs — stream-agnostic,
   λ-free, safe to delete.

This is also the only route that closes the census disagreement from the other side, since the same invariant is
what the sampler cannot see. Research-scale rather than overnight, and now the only identified route from the
46,286 to a submission.

## Standing

Nothing was removed and no configuration beating `1,490,805,286` was produced. Deleting any of the 46,286 on
9,024-shot evidence is exactly the stream-specificity error `04-traps.md` §1 documents. **Prove before removing** —
none of these three approaches proved anything.

## Reproducing

`hotness.rs` is retained as [`repro/hotness.rs`](repro/hotness.rs) because it carries §1 and the 46,286 count.
There is no `[[example]]` entry in `Cargo.toml`, so it is built the way `hyperplane_mitm.cpp` is — copied out,
compiled, scratch directory removed. From the repository root:

```sh
mkdir -p examples && cp src/point_add/memory/repro/hotness.rs examples/
cargo build --release --offline --example hotness
./target/release/examples/hotness -            # summary only
./target/release/examples/hotness /tmp/head    # also writes /tmp/head.hot.tsv, one row per gate
rm -rf examples
```

**Check the `GATE ok` line and that `avgT` matches `eval_circuit` before reading anything else** — the instrument
is meaningless unless it reconstructs the scorer's total, and the assertion is what caught the 9,216-shot draw.
It honours `SUB4_APPLY_STRIP` / `TLM_SCHED_J2_DELTA`, so run it under the same environment as the build being
priced. Verified on `cf5aa02` at 29 s build / 53 s run, reproducing every §1 figure above and four entries of
`06-research-status.md`'s certified table.

The two §3 rungs are not retained here, but they need no patching: `constzero.rs` and `affine.rs` compile and run
against an unmodified `cf5aa02` checkout exactly as `hotness.rs` does. Fetch them from
[`austinamissah/ecdsa-circuit-optimization`](https://github.com/austinamissah/ecdsa-circuit-optimization)
`tools/census/` at commit `ff5d6e0`, drop them in the same `examples/`, and:

```sh
cargo build --release --offline --example affine --example constzero
./target/release/examples/affine    --positive-control     # must print POSITIVE CONTROL PASS
./target/release/examples/constzero --check /tmp/head.hot.tsv
./target/release/examples/affine    --check /tmp/head.hot.tsv
```

Run the positive control first. Every `--check` requires that no certified gate fired on any of the 9,024 shots
and exits non-zero otherwise; both exit 0 on `cf5aa02` with 0 certified and 0 violations. The §2 census is the one
piece not reproducible from this tree — it needs fork-local mining infrastructure.
