#!/usr/bin/env python3
"""Formal (z3/SMT) proof that the Solinas modular-reduction algorithm used by
`mod_add_qq` (src/point_add/arith/modular.rs:12-49) computes (acc + a) mod p
for ALL acc, a in [0, p) on the secp256k1 field -- and that its overflow ancilla
(`flag`) uncomputes to |0>, i.e. the sub-circuit is clean/reversible.

This is the invariant the repo currently only argues in a comment ("Saves one
full (n+1)-wide Cuccaro compared to the sub-p/add-p/csub-p pattern") and checks
empirically via ECC signature traces. Here it is discharged as a theorem over
the full 256-bit field: z3 returns `unsat` on the negation.

The model mirrors the code step-for-step:

  ext_reg:  acc,a  -> 257-bit (top bit = transient overflow ancilla)
  Step 1:   acc_ext += a_ext                     # s = acc + a in [0, 2p)
  Step 2:   acc_ext += c   where c = 2^256 - p    # sets bit256 iff s >= p
  Step 3:   flag  = bit256(acc_ext)
  Step 4:   if flag==0: acc_ext -= c              # (X;csub;X) undoes step 2
  Step 5:   if flag==1: clear bit256              # drops 2^256 -> yields s - p
  Step 6:   flag ^= (acc_final < a)               # uncompute; must return to 0
"""
import sys
from z3 import (BitVec, BitVecVal, Solver, Extract, If, ULT, ULE, unsat)

# secp256k1 prime, exactly as src/point_add/mod.rs:677-682 (limbs, little-endian)
P = (0xFFFFFFFEFFFFFC2F
     | (0xFFFFFFFFFFFFFFFF << 64)
     | (0xFFFFFFFFFFFFFFFF << 128)
     | (0xFFFFFFFFFFFFFFFF << 192))
assert P == (1 << 256) - (1 << 32) - 977, "prime mismatch vs 2^256-2^32-977"

W = 257                       # extended register width (n + 1)
C = (1 << 256) - P            # Solinas constant c = 2^256 - p = 2^32 + 977
MASK256 = (1 << 256) - 1
Cv = BitVecVal(C, W)


def run():
    acc = BitVec('acc', W)
    a = BitVec('a', W)

    pre = [ULT(acc, BitVecVal(P, W)), ULT(a, BitVecVal(P, W))]  # acc, a in [0, p)

    # Step 1: 257-bit add. acc,a < p < 2^256 => s = acc+a < 2p < 2^257 (exact).
    s_val = acc + a
    # Step 2: add c.
    ext2 = s_val + Cv
    # Step 3: flag = top bit.
    flag = Extract(256, 256, ext2)
    # Step 4: undo the +c when flag==0 (the X;csub;X wrapper).
    ext4 = If(flag == BitVecVal(0, 1), ext2 - Cv, ext2)
    # Step 5: clear bit 256 when flag==1.
    ext5 = If(flag == BitVecVal(1, 1), ext4 & BitVecVal(MASK256, W), ext4)

    acc_final = Extract(255, 0, ext5)
    expected = Extract(255, 0, (acc + a) - If(ULE(BitVecVal(P, W), acc + a),
                                              BitVecVal(P, W), BitVecVal(0, W)))
    # (acc+a) mod p, since acc+a in [0,2p): subtract p iff acc+a >= p.

    # --- Theorem 1: value correctness ---
    s = Solver()
    s.add(*pre)
    s.add(acc_final != expected)
    r1 = s.check()
    ok1 = (r1 == unsat)
    print(f"  [{'PROVED ' if ok1 else 'FAILED '}] mod_add_qq: low256 == (acc + a) mod p  for all acc,a in [0,p)")
    if not ok1:
        print(f"      z3={r1}; counterexample={s.model()}")

    # --- Theorem 2: overflow ancilla uncomputes to 0 (step 6 clean) ---
    # Step 6 XORs (acc_final < a) into flag. Clean iff flag == (acc_final < a).
    a256 = Extract(255, 0, a)
    lt = If(ULT(acc_final, a256), BitVecVal(1, 1), BitVecVal(0, 1))
    s2 = Solver()
    s2.add(*pre)
    s2.add(flag != lt)
    r2 = s2.check()
    ok2 = (r2 == unsat)
    print(f"  [{'PROVED ' if ok2 else 'FAILED '}] mod_add_qq: overflow flag uncomputes to |0>  (flag == (acc_final < a))")
    if not ok2:
        print(f"      z3={r2}; counterexample={s2.model()}")

    return ok1 and ok2


if __name__ == "__main__":
    print("== Solinas modular reduction, secp256k1 p = 2^256 - 2^32 - 977 ==")
    print(f"   c = 2^256 - p = {C} = 0x{C:x}")
    ok = run()
    print(f"\n=== {'ALL PROVED' if ok else 'PROOF FAILED'} ===")
    sys.exit(0 if ok else 1)
