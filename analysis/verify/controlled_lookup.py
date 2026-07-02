#!/usr/bin/env python3
"""Construct and validate a CONTROLLED table-lookup kickmix circuit.

The source paper's shipped `table_lookup_3x3.kmx` is an illustrative extract: its
unary-iteration selector accumulator is driven by an outer control that is absent
from the standalone snippet, and it ships with no test vectors (only a `.svg`), so
it cannot be fuzz-validated as-is (see issue #3 and scientific-value.md §1d). This
module instead *constructs* a self-contained controlled lookup

    r0 ^= (ctrl ? table[addr] : 0)          # r0 = d-bit output, addr = a-bit,
                                            # table = 2^a x d classical bits, ctrl = 1 qubit

as a kickmix circuit, and validates it with the (reference-faithful)
`kickmix_sim.py`. This is the `3·2^w` windowed-lookup primitive of the ECDLP
ladder (ADR 0003) in a form we can actually check.

Two uncomputation strategies are built and validated, so both the plain reversible
and the measurement-based-uncomputation (MBUC, the paper's technique) forms of the
selector ancilla are exercised:

  - "reversible": the AND ladder is uncomputed by replaying the CCX gates.
  - "mbuc":       each selector qubit is cleared with HMR + a CZ phase correction
                  (a Z-basis function of its inputs), like the adders.

Validation (deterministic, seeded) requires, over exhaustive addresses x both
control values x random tables: correct output, address/table/control unchanged,
all selector ancilla returned to |0>, and global phase +1. The ctrl=0 case must
be a genuine no-op.
"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kickmix_sim import Circuit, State, simulate  # noqa: E402


def build_kmx(a, d, mode):
    """Return kmx text for a controlled a-bit-address, d-bit-data table lookup.

    Layout: out = q[0..d), addr = q[d..d+a), ctrl = q[d+a], scratch = q[d+a+1..].
    Registers r0=out, r1=addr, r2=table bits, r3=ctrl. Table bit table[k][j] is
    classical bit b[k*d+j]; the meas bit (mbuc) is b[2^a * d]."""
    out = list(range(d))
    addr = [d + i for i in range(a)]
    ctrl = d + a
    scratch = [d + a + 1 + t for t in range(a)]
    meas = (1 << a) * d
    L = []
    for q in out:
        L.append(f"APPEND_TO_REGISTER q{q} r0")
    for q in addr:
        L.append(f"APPEND_TO_REGISTER q{q} r1")
    for k in range(1 << a):
        for j in range(d):
            L.append(f"APPEND_TO_REGISTER b{k*d+j} r2")
    L.append(f"APPEND_TO_REGISTER q{ctrl} r3")

    for k in range(1 << a):
        # Select address == k: flip the address bits that are 0 in k, so that the
        # AND of all address bits is 1 exactly when addr == k.
        zeros = [addr[i] for i in range(a) if not ((k >> i) & 1)]
        for q in zeros:
            L.append(f"X q{q}")
        # Compute selector = ctrl AND addr[0] AND ... AND addr[a-1] via a CCX ladder.
        L.append(f"R q{scratch[0]}")
        L.append(f"CCX q{ctrl} q{addr[0]} q{scratch[0]}")
        for i in range(1, a):
            L.append(f"R q{scratch[i]}")
            L.append(f"CCX q{scratch[i-1]} q{addr[i]} q{scratch[i]}")
        sel = scratch[a - 1]
        # Write table[k] into the output where the selector is on.
        for j in range(d):
            L.append(f"CX q{sel} q{out[j]} if b{k*d+j}")
        # Uncompute the ladder.
        if mode == "reversible":
            for i in range(a - 1, 0, -1):
                L.append(f"CCX q{scratch[i-1]} q{addr[i]} q{scratch[i]}")
                L.append(f"R q{scratch[i]}")
            L.append(f"CCX q{ctrl} q{addr[0]} q{scratch[0]}")
            L.append(f"R q{scratch[0]}")
        elif mode == "mbuc":
            # scratch[i] = scratch[i-1] AND addr[i]; clear via measurement +
            # phase correction on its Z-basis inputs.
            for i in range(a - 1, 0, -1):
                L.append(f"HMR q{scratch[i]} b{meas}")
                L.append(f"CZ q{scratch[i-1]} q{addr[i]} if b{meas}")
            L.append(f"HMR q{scratch[0]} b{meas}")
            L.append(f"CZ q{ctrl} q{addr[0]} if b{meas}")
        else:
            raise ValueError(mode)
        for q in zeros:
            L.append(f"X q{q}")
    return "\n".join(L) + "\n"


def validate(a, d, mode, rng, table_shots=64):
    circ = Circuit.parse(build_kmx(a, d, mode))
    scratch = {d + a + 1 + t for t in range(a)}
    mod_d = 1 << d
    fails = 0
    first = ""
    for _ in range(table_shots):
        table = [rng.randrange(mod_d) for _ in range(1 << a)]
        tval = 0
        for k in range(1 << a):
            for j in range(d):
                if (table[k] >> j) & 1:
                    tval |= 1 << (k * d + j)
        for ctrl in (0, 1):
            for addr in range(1 << a):
                out_in = rng.randrange(mod_d)
                nq = max(circ.qubit_ids) + 1
                nb = max(circ.bit_ids) + 1
                st = State(nq, nb)
                circ.load_register(st, 0, out_in)   # r0 out
                circ.load_register(st, 1, addr)     # r1 addr
                circ.load_register(st, 2, tval)     # r2 table
                circ.load_register(st, 3, ctrl)     # r3 ctrl
                simulate(circ, st, rng)
                exp = out_in ^ (table[addr] if ctrl else 0)
                bad = None
                if circ.read_register(st, 0) != exp:
                    bad = f"out={circ.read_register(st,0)} exp={exp}"
                elif circ.read_register(st, 1) != addr:
                    bad = "addr changed"
                elif circ.read_register(st, 2) != tval:
                    bad = "table changed"
                elif circ.read_register(st, 3) != ctrl:
                    bad = "ctrl changed"
                elif any(st.qubits.get(s, 0) != 0 for s in scratch):
                    bad = "scratch not cleared"
                elif st.phase != 1:
                    bad = "phase -1"
                if bad:
                    fails += 1
                    first = first or f"a={a} d={d} {mode} ctrl={ctrl} addr={addr}: {bad}"
    return fails, first


def main():
    rng = random.Random(0x100C)
    print("=" * 74)
    print(" Constructed CONTROLLED table-lookup validation  (r0 ^= ctrl ? r2[r1] : 0)")
    print(" the 3*2^w ECDLP ladder lookup primitive, in a self-contained form")
    print("=" * 74)
    cases = [(3, 3), (2, 4), (4, 2)]
    all_ok = True
    for mode in ("reversible", "mbuc"):
        print(f"\n{mode.upper()} uncomputation")
        print("-" * 74)
        for a, d in cases:
            fails, first = validate(a, d, mode, rng)
            ok = fails == 0
            all_ok &= ok
            status = "PASS" if ok else f"FAIL ({fails}: {first})"
            print(f"  [{'ok' if ok else 'XX'}] a={a} addr bits, d={d} data bits "
                  f"({1<<a} entries) -> {status}")
    print("\n" + "=" * 74)
    if all_ok:
        print(" RESULT: controlled lookup verified (both reversible and MBUC forms):")
        print(" correct output, ctrl=0 is a no-op, ancilla cleared, phase +1.")
        print("=" * 74)
        return 0
    print(" RESULT: FAILURE")
    print("=" * 74)
    return 1


if __name__ == "__main__":
    sys.exit(main())
