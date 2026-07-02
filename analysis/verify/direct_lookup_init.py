#!/usr/bin/env python3
"""Circuit-level demonstration that the ladder's ∞-accumulator start is removed
by the direct-lookup first window (issue #5, part (a)).

Background. The running accumulator of the Shor-ECDLP ladder starts at ∞ with
**amplitude 1** (before any addition). Unlike the `dx=0` collisions, this cannot
be waved away by negligibility — it is certain, not rare — so it must be removed
*structurally*. The source paper (Babbush et al. 2026, Appendix A) does this by
replacing the **first windowed addition with a direct table lookup** that *writes*
the accumulator, instead of calling the adder on ∞. `completeness_argument.md §3`
and `completeness_collision_rate.py` both *assume* this; this module binds that
assumption to an actual reversible circuit.

What is demonstrated. Using the already-validated controlled table-lookup
primitive (`controlled_lookup.py`, issue #3/#9) as the init write
`acc ^= T[w]` on an accumulator register that starts at |0…0⟩, with the table
`T[w] = affine coords of [w]·P`:

  - the register ends holding exactly the coordinates of `[w]·P` (write is
    correct), the selector ancilla return to |0⟩, and the global phase is +1;
  - **the register is the `(0,0)` ∞ sentinel iff the window value `w = 0`.**

So for the *definite* (amplitude-1) initial window value the accumulator is a real
affine point — the adder is never fed ∞ at t=0. The residual `w=0` case is the
zero-window ∞ term (`~1/2^w` per addition, issue #5 part (b)), which the
negligibility bound already covers; the amplitude-1 start is gone by construction.

The property (`[w]P = ∞ ⇔ w ≡ 0`) is pure elliptic-curve arithmetic and is
curve-independent, so the exhaustive check runs on a small prime-order toy curve
(reusing `completeness_collision_rate.py`), with a secp256k1 spot-check at the real
256-bit coordinate width to show the write primitive composes at attack scale.

Analysis-only, deterministic, pure-Python (reuses `kickmix_sim.py`). Never touches
the scored circuit.
"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from completeness_collision_rate import INF, Curve, find_prime_order_curve  # noqa: E402
from controlled_lookup import build_kmx  # noqa: E402
from kickmix_sim import Circuit, State, simulate  # noqa: E402


def encode_point(pt, coord_bits):
    """Pack an affine point (or ∞) into a 2*coord_bits integer: x | y<<coord_bits.

    ∞ is the (0,0) sentinel -> 0. A finite point never encodes to 0, because on a
    prime-order (odd) curve no non-identity point has y = 0, so at least the y
    field is nonzero."""
    if pt is INF:
        return 0
    x, y = pt
    return (x & ((1 << coord_bits) - 1)) | (y << coord_bits)


def run_init_lookup(a, coord_bits, table_points, mode, rng, ctrl=1):
    """Build `acc ^= T[addr]` and run it for every window value on a |0⟩
    accumulator. `ctrl=1` is the direct-lookup init (the first window always
    writes); `ctrl=0` is the negative control (no write). Returns per-window
    (written_value, ancilla_clean, phase_ok, regs_preserved), where
    regs_preserved asserts the address (r1), table (r2), and ctrl (r3) registers
    are left untouched — the same invariant the upstream controlled_lookup.py
    validator enforces, so a lookup that silently corrupts an input register
    cannot pass here."""
    d = 2 * coord_bits
    circ = Circuit.parse(build_kmx(a, d, mode))
    scratch = {d + a + 1 + t for t in range(a)}

    # Pack the classical table T[w] (d bits each) into the r2 register integer.
    tval = 0
    for w in range(1 << a):
        tval |= encode_point(table_points[w], coord_bits) << (w * d)

    results = []
    for addr in range(1 << a):
        nq = max(circ.qubit_ids) + 1
        nb = max(circ.bit_ids) + 1
        st = State(nq, nb)
        circ.load_register(st, 0, 0)      # r0 acc = |0...0>  (fresh register)
        circ.load_register(st, 1, addr)   # r1 window value
        circ.load_register(st, 2, tval)   # r2 classical table
        circ.load_register(st, 3, ctrl)   # r3 ctrl
        simulate(circ, st, rng)
        out = circ.read_register(st, 0)
        anc_clean = all(st.qubits.get(s, 0) == 0 for s in scratch)
        regs_preserved = (
            circ.read_register(st, 1) == addr
            and circ.read_register(st, 2) == tval
            and circ.read_register(st, 3) == ctrl
        )
        results.append((out, anc_clean, st.phase == 1, regs_preserved))
    return results


def check_curve(label, c, gen, order, a, coord_bits, rng):
    """Exhaustively verify the direct-lookup init on one curve for both uncompute
    modes. Returns (ok, detail)."""
    assert (1 << a) <= order, "window range must not wrap the group order"
    table = [c.mul(w, gen) for w in range(1 << a)]  # T[w] = [w]P ; T[0] = ∞
    expect = [encode_point(pt, coord_bits) for pt in table]
    # Sanity on the EC side: ∞ (encode 0) occurs at exactly w == 0.
    if [w for w in range(1 << a) if expect[w] == 0] != [0]:
        return False, "toy-curve table: ∞ not uniquely at w=0"

    for mode in ("reversible", "mbuc"):
        res = run_init_lookup(a, coord_bits, table, mode, rng)
        for w, (out, anc, ph, regs) in enumerate(res):
            if out != expect[w]:
                return False, f"{label}/{mode} w={w}: wrote {out}, expected {expect[w]}"
            if not anc:
                return False, f"{label}/{mode} w={w}: selector ancilla not cleared"
            if not ph:
                return False, f"{label}/{mode} w={w}: phase -1"
            if not regs:
                return False, f"{label}/{mode} w={w}: address/table/ctrl register corrupted"
            # The load-bearing property: register == ∞ sentinel iff w == 0.
            if (out == 0) != (w == 0):
                return False, f"{label}/{mode} w={w}: ∞-sentinel/window mismatch"
        # Negative control: with no write (ctrl=0) the accumulator stays ∞ (all
        # zero) for every window — so it is the lookup WRITE, not an artifact,
        # that turns the ∞ start into a real point. The ancilla, phase, and input
        # registers must be preserved on this no-op path too.
        neg = run_init_lookup(a, coord_bits, table, mode, rng, ctrl=0)
        for w, (out, anc, ph, regs) in enumerate(neg):
            if out != 0:
                return False, f"{label}/{mode}: ctrl=0 did not leave the accumulator at ∞"
            if not (anc and ph and regs):
                return False, (f"{label}/{mode} w={w}: ctrl=0 no-op corrupted "
                               "ancilla/phase/registers")
    return True, f"{label}: [w]P written for all {1<<a} windows; ∞ iff w=0; ctrl=0 stays ∞"


def secp256k1():
    p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
    c = Curve(p, 0, 7)
    c.order = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
    g = (
        0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
        0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8,
    )
    return c, g, c.order


def main():
    rng = random.Random(0x0A11)
    print("=" * 74)
    print(" Direct-lookup first window: circuit-level ∞-start removal (issue #5a)")
    print(" acc |0> ^= T[w],  T[w] = [w]P  ->  accumulator is a real point (w!=0)")
    print("=" * 74)
    print()

    ok = True

    # (1) Exhaustive on a small prime-order toy curve, both uncompute modes.
    c, gen, n = find_prime_order_curve()
    coord_bits = c.p.bit_length()
    a = 3  # 8-window table; 2^3 < n so no window value is ≡ 0 except 0 itself
    print("Toy prime-order curve (exhaustive, reversible + MBUC uncompute)")
    print("-" * 74)
    print(f"  curve y^2 = x^3 + {c.a}x + {c.b} over F_{c.p}, prime order n={n}, "
          f"{coord_bits}-bit coords, {1<<a}-entry table")
    good, detail = check_curve("toy", c, gen, n, a, coord_bits, rng)
    ok &= good
    print(f"  [{'ok' if good else 'XX'}] {detail}\n")

    # (2) secp256k1 spot-check at real 256-bit coordinate width (MBUC only,
    #     4-entry table) — shows the write composes at attack scale.
    cs, gs, ns = secp256k1()
    scb = 256
    aa = 2
    table = [cs.mul(w, gs) for w in range(1 << aa)]
    expect = [encode_point(pt, scb) for pt in table]
    res = run_init_lookup(aa, scb, table, "mbuc", rng)
    good2 = all(res[w][0] == expect[w] and res[w][1] and res[w][2] and res[w][3]
                for w in range(1 << aa)) \
        and [w for w in range(1 << aa) if res[w][0] == 0] == [0]
    ok &= good2
    print("secp256k1 spot-check (real 256-bit coords, MBUC, 4-entry table)")
    print("-" * 74)
    for w in range(1 << aa):
        pt = "∞" if table[w] is INF else f"({hex(table[w][0])[:10]}.., ..)"
        print(f"  w={w}: [{w}]G = {pt:<20} written_ok={res[w][0]==expect[w]} "
              f"is_inf={res[w][0]==0}")
    print(f"  [{'ok' if good2 else 'XX'}] write correct at 256-bit width; ∞ iff w=0\n")

    # (3) The contrast that makes this load-bearing: what the init removes.
    print("What the direct-lookup init removes (vs. a naive ∞ start)")
    print("-" * 74)
    print(f"  naive (acc:=∞, then ec_add(∞, T0[w0])):  P[first add sees ∞] = 1.000  "
          f"(amplitude-1, fatal)")
    print(f"  direct-lookup init (acc := T0[w0]):       P[first add sees ∞] = 1/2^w "
          f"= P[w0=0]  (the negligible zero-window term, issue #5b)")
    print("  => the amplitude-1 ∞ start is removed structurally, not by "
          "negligibility.\n")

    print("=" * 74)
    if ok:
        print(" RESULT: the direct-lookup first window writes the accumulator to a")
        print(" real affine point for every nonzero window (∞ only at w=0), so the")
        print(" adder is never fed the amplitude-1 ∞ start. Issue #5 (a) is now")
        print(" circuit-demonstrated on the repo's validated QROM primitive.")
        print("=" * 74)
        return 0
    print(" RESULT: FAILURE — see [XX] above.")
    print("=" * 74)
    return 1


if __name__ == "__main__":
    sys.exit(main())
