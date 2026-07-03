# ADR 0014 — Quantum-addend point-add testbed: a QROM-fed add, verified by simulation (Tier B, issue #27/#28)

**Status:** Accepted — decides the approach for the deferred quantum-addend build
that ADRs 0011/0012/0013 all point to; first increment implemented in
`src/point_add/qaddend_testbed.rs`.
**Date:** 2026-07-03

## Context

Three ADRs converge on one deferred build — a **functionally-correct QROM-fed
quantum addend**:

- [ADR 0011](0011-streamed-full-ladder.md) emits the full ladder but keeps the
  QROM lookup on *disjoint* ids, flagging that the real read→add data dependency
  and register overlap need the quantum-addend PA.
- [ADR 0012](0012-classical-vs-quantum-addend-gap.md) measured the *Toffoli* gap
  negligible: the scored PA already loads its classical addend into a qubit
  register and runs an uncontrolled q-q add (`coord_addsub`), so the arithmetic is
  addend-value-independent.
- [ADR 0013](0013-quantum-addend-width-gap.md) measured the *width* gap real: a
  QROM addend held resident across the GCD peak adds +256..512 qubits.

What none of them do is **exhibit the composition** — read `P[k]` from a quantum
table *into* an addend register, have an adder consume it, and unread — and prove
by simulation that it computes the right sum with all ancilla returned to |0>.
That composition is the keystone of turning the cost estimate into a verified
attack, and issue #28's mid-ladder ∞/`dx=0` residual needs the same testbed.

## Decision

Build a **self-contained, width-parametric** QROM-fed quantum-addend adder as a
`#[cfg(test)]` harness and verify it by **simulation**, rather than mutating the
scored 256-bit point-add:

1. **Reduced, parametric width** (`n`-bit registers, `w`-bit window, prime-free
   integer add mod `2^n` in v1). Justified because (a) the arithmetic is
   width-parametric and addend-value-independent (ADR 0012), so a small instance
   faithfully exercises the *new* surface — a quantum table feeding an adder —
   without the scored circuit's product-min tuning; (b) full-width simulation of
   a QROM-fed 256-bit modular point-add is impractical for exhaustive functional
   checks; (c) it isolates the read→add→unread composition and the register
   overlap from everything already measured.

2. **Real QROM construction.** The read is the same **unary-iteration** selector
   (Gidney 2018 §III.C) as `ladder_lookup_cost.py` (ADR 0010) and
   `ladder_full.rs`, but now WITH the leaf **data-writes** (`CX` of the classical
   table constants into the addend register) that those cost-only harnesses omit.
   So this testbed *functionally validates* that QROM construction, not just its
   Toffoli count.

3. **Real adder shape.** The addend register is consumed by an **uncontrolled
   quantum-quantum ripple-carry (Cuccaro) add**, addend preserved — the exact
   shape `coord_addsub` uses (ADR 0012). The `unread` is a second application of
   the reversible read, returning the addend to |0>.

4. **Verify by simulation** over *all* `2^w` window values and multiple
   accumulator inputs (masked multi-shot, as the existing `*_selftest`
   harnesses): assert `acc' == (acc + P[k]) mod 2^n`, and that the addend, the
   selector ancilla, and the window register are all clean/preserved afterward.

5. **Measure the overlap.** Report the peak width and its breakdown
   (window + addend + accumulator + selector ancilla + carry), the small-scale
   analogue of ADR 0013's full-width +256..512.

## Consequences

- **Functional correctness of the QROM-fed add is demonstrated by construction**,
  closing the "does the composition even work" question that ADR 0011 deferred —
  and independently validating the unary-iteration QROM *with data writes*.
- **The register-overlap picture becomes concrete and executable**: the addend and
  selector ancilla are explicit registers whose peak contribution is measured, not
  argued, complementing ADR 0013's full-width construction.
- **Scope of v1** (honest): integer add mod `2^n` (not yet a field-modular
  reduction) and a single table read→add→unread (not the full ladder). The
  modular-reduction tail is a fixed, addend-independent extension (as
  `coord_addsub`'s `mod_add` shows); issue #28's EC exceptional cases (`P==Q`,
  `dx=0`, ∞) need the *group law* on top and are the next increments on this same
  testbed. These are called out, not silently omitted.
- Consistent with [ADR 0001](0001-analysis-layer-isolated-from-score.md): the
  harness is `#[cfg(test)]`, never compiled into `build_circuit`; the scored
  circuit is byte-identical (`ops.bin` SHA unchanged).
