#!/usr/bin/env python3
"""Validate the source paper's reference in-place adder circuits, and confirm its
negative-control circuits are rejected, using an independent kickmix simulator.

WHY THIS MATTERS FOR THIS REPO
  This repo's scored circuit is a kickmix elliptic-curve point addition whose
  arithmetic core is a Cuccaro ripple-carry adder (arith/adder.rs). The source
  paper (Babbush et al. 2026, arXiv:2603.28846v2) ships reference *iadd* circuits
  in original/zkp_ecc_zenodo_v2/docs/example_data -- and iadd8.kmx/iadd64.kmx are
  explicitly "a variant of the adder from quant-ph/0410184" (Cuccaro et al.), the
  SAME primitive this repo uses. iadd8_with_classical_offset_and_dirty_ancillae
  even adds a *classical* addend register with *dirty* ancilla and MBUC phase
  correction -- structurally the same shape as this repo's "quantum point +=
  classical point" primitive.

  We simulate those reference circuits with analysis/verify/kickmix_sim.py (an
  independent re-derivation of the kickmix semantics that src/sim.rs implements)
  and check, across fuzzed inputs, that:
    (a) they compute r0 += r1 (mod 2^n) with the addend register unchanged,
    (b) all workspace ancilla return to |0> and declared dirty ancilla are
        restored to their input, and
    (c) the tracked global phase stays +1 (MBUC phase kickback fully corrected).
  We also confirm the paper's three *incorrect* increment circuits FAIL these
  checks -- validating that the phase-/ancilla-aware methodology detects bugs,
  which is the paper's whole correctness argument (Appendix A.5) and mirrors
  eval_circuit's classical/phase/ancilla-garbage checks.

Everything is deterministic (fixed RNG seed); no number is hand-asserted.
"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kickmix_sim import Circuit, State, simulate  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
# Prefer the vendored copy (self-contained); fall back to the full Zenodo dump.
_VENDORED = os.path.join(HERE, "reference_circuits")
_ZENODO = os.path.normpath(os.path.join(
    HERE, "..", "..", "original", "zkp_ecc_zenodo_v2", "docs", "example_data"))
DATA = _VENDORED if os.path.isdir(_VENDORED) else _ZENODO
SEED = 0xECD5A
SHOTS = 2000
NEG_SHOTS = 400


def load(name):
    with open(os.path.join(DATA, name)) as f:
        return Circuit.parse(f.read())


def reg_width(circ, r):
    return len(circ.registers[r])


def ancilla_split(circ):
    """Classify non-register qubits into (clean_workspace, dirty_borrowed).

    A clean workspace qubit is reset by an R instruction (so it must start |0>);
    a dirty borrowed qubit is never reset (so it starts arbitrary and must be
    RESTORED to its input) -- this is the dirty-ancilla technique the paper's
    classical-offset variant exercises."""
    reg_q = {i for elems in circ.registers.values() for (k, i) in elems if k == "q"}
    r_targets = {i for (name, tg, _c) in circ.instructions if name == "R"
                 for (k, i) in tg if k == "q"}
    non_reg = circ.qubit_ids - reg_q
    return (non_reg & r_targets), (non_reg - r_targets)


def run_shot(circ, reg_in, reg_expected, rng):
    """Init registers + ancilla, run, and return (ok, reason).

    Clean workspace starts |0> and must return to |0>; dirty borrowed ancilla
    start RANDOM and must be restored to that value; the global phase must be +1."""
    nq = (max(circ.qubit_ids) + 1) if circ.qubit_ids else 0
    nb = (max(circ.bit_ids) + 1) if circ.bit_ids else 0
    st = State(nq, nb)
    for r, v in reg_in.items():
        circ.load_register(st, r, v)
    clean, dirty = ancilla_split(circ)
    dirty_init = {q: (1 if rng.random() < 0.5 else 0) for q in dirty}
    for q, v in dirty_init.items():
        st.qubits[q] = v
    simulate(circ, st, rng)
    for r, exp in reg_expected.items():
        got = circ.read_register(st, r)
        if got != exp:
            return False, f"register r{r} = {got}, expected {exp}"
    for qid in clean:
        if st.qubits.get(qid, 0) != 0:
            return False, f"clean workspace q{qid} left dirty ({st.qubits[qid]})"
    for qid, v0 in dirty_init.items():
        if st.qubits.get(qid, 0) != v0:
            return False, f"dirty ancilla q{qid} not restored ({st.qubits[qid]} != {v0})"
    if st.phase != 1:
        return False, "global phase ended -1 (uncorrected phase kickback)"
    return True, ""


def fuzz_adder(circ, computed, addend, carry, passthrough, shots, rng):
    """r{computed} = (r{computed} + r{addend} + r{carry}) mod 2^width.

    `carry` may be None (plain adder). `passthrough` regs (incl. addend/carry)
    are restored/unchanged. Returns (n_fail, first_reason)."""
    mod = 1 << reg_width(circ, computed)
    fails, reason = 0, ""
    for _ in range(shots):
        reg_in = {computed: rng.randrange(mod), addend: rng.randrange(1 << reg_width(circ, addend))}
        if carry is not None:
            reg_in[carry] = rng.randrange(1 << reg_width(circ, carry))
        for r in passthrough:
            reg_in.setdefault(r, rng.randrange(1 << reg_width(circ, r)))
        total = reg_in[computed] + reg_in[addend] + (reg_in[carry] if carry is not None else 0)
        expected = {computed: total % mod, addend: reg_in[addend]}
        if carry is not None:
            expected[carry] = reg_in[carry]
        for r in passthrough:
            expected[r] = reg_in[r]
        ok, why = run_shot(circ, reg_in, expected, rng)
        if not ok:
            fails += 1
            reason = reason or why
    return fails, reason


def parse_golden(name):
    cases = []
    with open(os.path.join(DATA, name)) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            lhs, rhs = line.split("->")
            cases.append((int(lhs.strip()), int(rhs.strip())))
    return cases


def golden_inc(circ, cases, shots_each, rng):
    """Single register r0; each case is (in -> out). Returns (n_fail, reason)."""
    fails, reason = 0, ""
    for a, out in cases:
        for _ in range(shots_each):
            ok, why = run_shot(circ, {0: a}, {0: out}, rng)
            if not ok:
                fails += 1
                reason = reason or f"input {a}: {why}"
    return fails, reason


def main():
    if not os.path.isdir(DATA):
        print(f"SKIP: reference data not found at {DATA}")
        print("(the original/ Zenodo dump is untracked; nothing to validate)")
        return 0

    rng = random.Random(SEED)
    print("=" * 74)
    print(" Reference kickmix adder validation  (source: arXiv:2603.28846v2, Appendix A)")
    print(f" data: original/.../example_data   seed=0x{SEED:X}   shots/circuit={SHOTS}")
    print("=" * 74)

    inc_cases = parse_golden("inc3_test_cases.txt")

    positives = [
        ("inc3.kmx", "golden", lambda c: golden_inc(c, inc_cases, 64, rng)),
        ("iadd8.kmx", "fuzz r0+=r1", lambda c: fuzz_adder(c, 0, 1, None, [], SHOTS, rng)),
        ("iadd8_with_ancillae.kmx", "fuzz r0+=r1", lambda c: fuzz_adder(c, 0, 1, None, [], SHOTS, rng)),
        ("iadd64.kmx", "fuzz r0+=r1", lambda c: fuzz_adder(c, 0, 1, None, [], SHOTS // 4, rng)),
        ("iadd8_with_classical_offset_and_dirty_ancillae.kmx", "fuzz r0+=r1+cin",
         lambda c: fuzz_adder(c, 0, 1, 2, [], SHOTS, rng)),
    ]
    negatives = [
        ("inc3_wrong_order.kmx", lambda c: golden_inc(c, inc_cases, NEG_SHOTS // 8, rng)),
        ("inc3_wrong_phase.kmx", lambda c: golden_inc(c, inc_cases, NEG_SHOTS // 8, rng)),
        ("inc3_wrong_garbage.kmx", lambda c: golden_inc(c, inc_cases, NEG_SHOTS // 8, rng)),
    ]

    all_ok = True

    print("\nPOSITIVE controls (must PASS: correct output, clean ancilla, phase +1)")
    print("-" * 74)
    for name, kind, fn in positives:
        circ = load(name)
        fails, reason = fn(circ)
        status = "PASS" if fails == 0 else f"FAIL ({fails} bad shots: {reason})"
        all_ok &= fails == 0
        print(f"  [{'ok' if fails == 0 else 'XX'}] {name:<52} {kind:<16} -> {status}")

    print("\nNEGATIVE controls (must be REJECTED: at least one shot fails)")
    print("-" * 74)
    for name, fn in negatives:
        circ = load(name)
        fails, reason = fn(circ)
        detected = fails > 0
        status = f"rejected ({fails} bad shots)" if detected else "NOT DETECTED -- validator too weak!"
        all_ok &= detected
        print(f"  [{'ok' if detected else 'XX'}] {name:<52} {"neg":<16} -> {status}")

    print("\n" + "=" * 74)
    if all_ok:
        print(" RESULT: all reference adders verified; all negative controls rejected.")
        print(" The kickmix semantics this repo relies on reproduce the paper's artifacts.")
        print("=" * 74)
        return 0
    print(" RESULT: FAILURE -- a positive control failed or a negative slipped through.")
    print("=" * 74)
    return 1


if __name__ == "__main__":
    sys.exit(main())
