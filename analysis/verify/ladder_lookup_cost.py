#!/usr/bin/env python3
"""Measure the windowed-lookup (QROM) Toffoli/ancilla cost of the ECDLP ladder
— the largest term in the full-attack estimate that was still *derived*, not
measured (issue #4).

`ecdlp_estimate.py` uses the source paper's closed form
`ECDLP_Toff = (PA_Toff + 3·2^w)(2n/w − 4)`. `PA_Toff` (the per-addition point-add)
is measured from the scored circuit; the `3·2^w` **table-lookup** term — which
loads the precomputed multiple `P[k]` from a `2^w`-entry table indexed by the
`w`-qubit window register — is taken from the paper (ADR 0003). The scored circuit
uses a *classical* addend, so it contains **no** QROM; this term had nothing
measured behind it.

This module builds the lookup as an **optimized unary-iteration QROM**
(Gidney 2018 §III.C — the same primitive the paper cites) as an actual kickmix
circuit, **validates it** exhaustively with the reference-faithful `kickmix_sim.py`
(correct read, address/table unchanged, all iteration ancilla returned to |0⟩,
global phase +1), and **measures** its Toffoli count and ancilla width. A single
QROM read of `N = 2^w` entries costs:

    Toffoli(reversible) = 2^(w+1) − 4        (compute + reversible uncompute)
    Toffoli(compute-only, MBUC bound) = 2^w − 2   (uncompute → Clifford, the
                                                    paper's measurement technique)
    ancilla = w

both **below** the paper's `3·2^w` per-addition budget — so the estimate's lookup
term is conservative, and now has a validated construction behind it. (A windowed
addition uses the read plus an uncompute-read; even counting both, MBUC keeps the
per-addition lookup Toffoli near `2·2^w ≤ 3·2^w`.)

Unary iteration recursion (reused single-ancilla-per-level "spine", w ancilla):

    iterate(level L, control c, prefix):
      a := scratch[L]
      a ^= c AND addr[L]            # CCX (or CX at the always-on top level)
      iterate(L+1, a, prefix|bitL)  # addr[L] == 1 subtree
      a ^= c                        # now a = c AND ¬addr[L]
      iterate(L+1, a, prefix)       # addr[L] == 0 subtree
      a ^= c ;  a ^= c AND addr[L]  # restore + uncompute a -> 0
    leaf: out[j] ^= a  for each data bit j set in table[prefix]   # CX (Clifford)

Analysis-only, deterministic, pure-Python (reuses `kickmix_sim.py`). Never touches
the scored circuit.
"""
import io
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kickmix_sim import Circuit, State, simulate  # noqa: E402


def build_qrom_kmx(w, d):
    """Return kmx text for an unconditional unary-iteration read `out ^= T[addr]`.

    Layout: out = q[0..d), addr = q[d..d+w), scratch = q[d+w..d+2w). Registers
    r0=out, r1=addr, r2=table (classical bits: T[k] bit j is b[k*d+j])."""
    out = list(range(d))
    addr = [d + i for i in range(w)]
    scratch = [d + w + i for i in range(w)]
    L = []
    for q in out:
        L.append(f"APPEND_TO_REGISTER q{q} r0")
    for q in addr:
        L.append(f"APPEND_TO_REGISTER q{q} r1")
    for k in range(1 << w):
        for j in range(d):
            L.append(f"APPEND_TO_REGISTER b{k*d+j} r2")

    def emit(level, ctrl, prefix):
        if level == w:  # leaf: write T[prefix] into out, controlled on `ctrl`
            for j in range(d):
                L.append(f"CX q{ctrl} q{out[j]} if b{prefix*d + j}")
            return
        a, ab = scratch[level], addr[level]
        # a = ctrl AND addr[level]   (top level: ctrl is always-on -> a = addr)
        L.append(f"CX q{ab} q{a}" if ctrl is None else f"CCX q{ctrl} q{ab} q{a}")
        emit(level + 1, a, prefix | (1 << level))          # addr[level] == 1
        # a ^= ctrl  =>  a = ctrl AND NOT addr[level]
        L.append(f"X q{a}" if ctrl is None else f"CX q{ctrl} q{a}")
        emit(level + 1, a, prefix)                         # addr[level] == 0
        # restore a = ctrl AND addr[level], then uncompute a -> 0
        if ctrl is None:
            L.append(f"X q{a}")
            L.append(f"CX q{ab} q{a}")
        else:
            L.append(f"CX q{ctrl} q{a}")
            L.append(f"CCX q{ctrl} q{ab} q{a}")

    emit(0, None, 0)
    return "\n".join(L) + "\n", scratch


def validate(w, d, rng, table_shots=24):
    """Exhaustive over addresses x random tables: correct read, addr/table
    unchanged, iteration ancilla cleared, phase +1."""
    text, scratch = build_qrom_kmx(w, d)
    circ = Circuit.parse(text)
    scratch = set(scratch)
    mod_d = 1 << d
    for _ in range(table_shots):
        table = [rng.randrange(mod_d) for _ in range(1 << w)]
        tval = 0
        for k in range(1 << w):
            tval |= table[k] << (k * d)
        for addr in range(1 << w):
            nq = max(circ.qubit_ids) + 1
            nb = max(circ.bit_ids) + 1
            st = State(nq, nb)
            circ.load_register(st, 1, addr)   # r1 addr  (r0 out defaults to 0)
            circ.load_register(st, 2, tval)   # r2 table
            simulate(circ, st, rng)
            if circ.read_register(st, 0) != table[addr]:
                return False, f"w={w} d={d} addr={addr}: read {circ.read_register(st,0)} != {table[addr]}"
            if circ.read_register(st, 1) != addr:
                return False, f"w={w} d={d} addr={addr}: address register changed"
            if circ.read_register(st, 2) != tval:
                return False, f"w={w} d={d} addr={addr}: table changed"
            if any(st.qubits.get(s, 0) != 0 for s in scratch):
                return False, f"w={w} d={d} addr={addr}: iteration ancilla not cleared"
            if st.phase != 1:
                return False, f"w={w} d={d} addr={addr}: phase -1"
    return True, ""


def measure(w):
    """Count Toffoli (CCX), the compute-only subset, and ancilla width.

    The Toffoli count and ancilla width depend only on `w` (the unary-iteration
    spine), not on the data width `d` or the table contents. So build the spine
    with `d=0`: that omits the O(2^w * d) table/output lines entirely — avoiding a
    large transient for big `w` (e.g. w=16) — while leaving the CCX spine and the
    `w` scratch ancilla identical. Scan the text line-by-line via a StringIO
    iterator rather than `splitlines()`, so no list of all lines is materialized."""
    text, scratch = build_qrom_kmx(w, 0)
    ccx = sum(1 for ln in io.StringIO(text) if ln.startswith("CCX "))
    return {"w": w, "toffoli": ccx, "compute_only": ccx // 2, "ancilla": len(scratch)}


def main():
    rng = random.Random(0x1A4D)
    print("=" * 74)
    print(" Windowed-lookup (QROM) cost: measured unary-iteration table read")
    print(" out ^= T[addr]   —   grounds the derived 3*2^w lookup term (issue #4)")
    print("=" * 74)
    print()

    ok = True

    # (1) Validate the construction is a correct, clean, phase-+1 QROM read.
    print("Validation (exhaustive addresses x random tables, kickmix sim)")
    print("-" * 74)
    for w, d in [(2, 3), (3, 3), (4, 2), (5, 2), (6, 2)]:
        good, detail = validate(w, d, rng)
        ok &= good
        print(f"  [{'ok' if good else 'XX'}] w={w} ({1<<w:>3} entries), d={d} data bits"
              f"{'' if good else '  ' + detail}")
    print()

    # (2) Measure Toffoli(w) / ancilla(w) and compare to the paper's 3*2^w.
    print("Measured cost vs. the paper's 3*2^w per-addition lookup budget")
    print("-" * 74)
    hdr = (f"  {'w':>3} {'entries':>8} {'Toffoli':>10} {'Tof/2^w':>8} "
           f"{'MBUC(=2^w-2)':>12} {'ancilla':>8} {'paper 3*2^w':>12} {'ratio':>6}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    rows = []
    for w in [2, 3, 4, 6, 8, 10, 12, 16]:
        m = measure(w)
        paper = 3 * (1 << w)
        rows.append((m, paper))
        print(f"  {w:>3} {1<<w:>8} {m['toffoli']:>10} {m['toffoli']/(1<<w):>8.2f} "
              f"{m['compute_only']:>12} {m['ancilla']:>8} {paper:>12} "
              f"{m['toffoli']/paper:>6.2f}")
    print()

    # Regression + finding locks.
    notes = []
    closed_form = all(m["toffoli"] == (1 << (m["w"] + 1)) - 4 for m, _ in rows)
    notes.append(f"[{'ok' if closed_form else 'XX'}] Toffoli(w) == 2^(w+1) - 4 "
                 f"(deterministic closed form)")
    ok &= closed_form
    beats = all(m["toffoli"] <= paper for m, paper in rows)
    notes.append(f"[{'ok' if beats else 'XX'}] measured read Toffoli <= paper's 3*2^w "
                 f"on every w (~0.67x)")
    ok &= beats
    anc = all(m["ancilla"] == m["w"] for m, _ in rows)
    notes.append(f"[{'ok' if anc else 'XX'}] ancilla == w (O(w) spine, matches the "
                 f"paper's PA_Qubits + w)")
    ok &= anc
    print("Findings")
    print("-" * 74)
    for ln in notes:
        print("  " + ln)
    m16 = next(m for m, _ in rows if m["w"] == 16)  # reuse the row already measured
    print(f"\n  At the optimal window w=16: measured read = {m16['toffoli']:,} Toffoli "
          f"({m16['ancilla']} ancilla),")
    print(f"  vs the paper's 3*2^16 = {3*(1<<16):,}. A windowed addition needs a read")
    print(f"  plus an uncompute-read; with measurement-based uncompute (the repo's")
    print(f"  technique) the per-addition lookup stays ~2*2^w, within the 3*2^w budget.")
    print()

    print("=" * 74)
    if ok:
        print(" RESULT: the windowed-lookup term now has a validated construction and")
        print(" a measured cost (2^(w+1)-4 Toffoli, w ancilla per read) — below the")
        print(" paper's derived 3*2^w. The estimate's lookup term is grounded, not")
        print(" just cited (issue #4; ADR 0010).")
        print("=" * 74)
        return 0
    print(" RESULT: FAILURE — see [XX] above.")
    print("=" * 74)
    return 1


if __name__ == "__main__":
    sys.exit(main())
