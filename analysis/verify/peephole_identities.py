#!/usr/bin/env python3
"""Formal (z3/SMT) proofs of the boolean-logic invariants the peephole and
uncompute optimizations rely on.

These are the claims currently checked only *empirically* (CONSTPROP_VERIFY,
ALT_SEED_*): re-running the op-stream over sampled shots and asserting the
transform preserved the state. Here we discharge them as universally-quantified
theorems: z3 returns `unsat` on the negation, i.e. NO assignment violates the
identity -> proven for all inputs, not just the sampled ones.

Sources in the repo:
  - constprop.rs (DropZeroCtrl / FoldCx / FoldX / FoldEqualCtrls /
    DropComplementCtrls / InversePairCancellation)
  - venting.rs, arith/adder.rs (Cuccaro / HRS carry-xor recurrence)
  - trailmix_ludicrous/comparator.rs (a<b flag via carry chain)
"""
import sys
from z3 import (BitVec, BitVecVal, Solver, Extract, Concat, If, Not, And, ULT,
                unsat, sat)

results = []


def prove(name, claim_solver_setup):
    """claim_solver_setup(s) must add the NEGATION of the claim to solver s.
    unsat => claim holds for all inputs."""
    s = Solver()
    claim_solver_setup(s)
    r = s.check()
    ok = (r == unsat)
    results.append((name, ok, r))
    status = "PROVED " if ok else "FAILED "
    print(f"  [{status}] {name}" + ("" if ok else f"  (z3: {r}; counterexample: {s.model()})"))
    return ok


bit = lambda nm: BitVec(nm, 1)
ONE = BitVecVal(1, 1)
ZERO = BitVecVal(0, 1)
# CCX action on target: t' = t XOR (a AND b)
ccx = lambda a, b, t: t ^ (a & b)

print("== CCX peephole identities (constprop.rs) ==")

# 1.1 DropZeroCtrl: control a==0  =>  CCX is identity on target
prove("1.1 DropZeroCtrl:  a=0  =>  CCX(a,b,t) == t", lambda s: (
    s.add(bit('a') == ZERO),
    s.add(ccx(bit('a'), bit('b'), bit('t')) != bit('t')),
))

# 1.2 FoldCx: control a==1  =>  CCX(a,b,t) == CX(b,t) == t XOR b
prove("1.2 FoldCx:       a=1  =>  CCX(a,b,t) == t XOR b", lambda s: (
    s.add(bit('a') == ONE),
    s.add(ccx(bit('a'), bit('b'), bit('t')) != (bit('t') ^ bit('b'))),
))

# 1.3 FoldX: a==1 and b==1  =>  CCX == X(t) == t XOR 1
prove("1.3 FoldX:        a=1,b=1 => CCX(a,b,t) == NOT t", lambda s: (
    s.add(bit('a') == ONE, bit('b') == ONE),
    s.add(ccx(bit('a'), bit('b'), bit('t')) != (bit('t') ^ ONE)),
))

# 1.4 FoldEqualCtrls: affine analysis proves a==b on every shot => CCX == CX(a,t)
prove("1.4 FoldEqualCtrls:      a==b  =>  CCX(a,b,t) == t XOR a", lambda s: (
    s.add(bit('a') == bit('b')),
    s.add(ccx(bit('a'), bit('b'), bit('t')) != (bit('t') ^ bit('a'))),
))

# 1.5 DropComplementCtrls: a == NOT b on every shot => CCX is a no-op
prove("1.5 DropComplementCtrls: a==~b =>  CCX(a,b,t) == t", lambda s: (
    s.add(bit('a') == (bit('b') ^ ONE)),
    s.add(ccx(bit('a'), bit('b'), bit('t')) != bit('t')),
))

# 2.1 InversePairCancellation: two identical CCX with controls+target untouched
#     between them cancel. v = a AND b applied twice: (t^v)^v == t.
prove("2.1 InversePairCancellation: CCX;CCX (ctrls/target unchanged) == identity", lambda s: (
    s.add(ccx(bit('a'), bit('b'), ccx(bit('a'), bit('b'), bit('t'))) != bit('t')),
))

print("\n== Ripple-carry adder recurrence (Cuccaro / HRS carry-xor, venting.rs) ==")


def prove_adder(w):
    """Prove the carry recurrence used by the vented/HRS/Cuccaro adders computes
    exact integer addition mod 2^w, for a symbolic w-bit a,b."""
    def setup(s):
        a = BitVec(f'a{w}', w)
        b = BitVec(f'b{w}', w)
        # bit-serial ripple: carry c[0]=0; sum[i]=a[i]^b[i]^c[i];
        # c[i+1] = majority(a[i],b[i],c[i]) = (a&b)|(c&(a^b))
        carry = ZERO
        sum_bits = []
        for i in range(w):
            ai = Extract(i, i, a)
            bi = Extract(i, i, b)
            si = ai ^ bi ^ carry
            sum_bits.append(si)
            carry = (ai & bi) | (carry & (ai ^ bi))
        # assemble little-endian
        ssum = sum_bits[0]
        for i in range(1, w):
            ssum = Concat(sum_bits[i], ssum)
        s.add(ssum != (a + b))  # a+b already truncates to w bits (mod 2^w)
    return prove(f"ripple-carry recurrence == (a+b) mod 2^{w}", setup)


for w in (1, 2, 3, 4, 8, 16, 32, 64):
    prove_adder(w)

print("\n== Less-than comparator via borrow chain (comparator.rs) ==")


def prove_cmp(w):
    """flag := (a < b) computed by the subtract-borrow chain equals ULT(a,b)."""
    def setup(s):
        a = BitVec(f'ca{w}', w)
        b = BitVec(f'cb{w}', w)
        # a < b  iff  a - b borrows out of the top bit.
        # borrow[0]=0; borrow[i+1] = (~a[i] & b[i]) | (borrow[i] & ~(a[i]^b[i]))
        borrow = ZERO
        for i in range(w):
            ai = Extract(i, i, a)
            bi = Extract(i, i, b)
            borrow = ((~ai) & bi) | (borrow & (~(ai ^ bi)))
        flag_lt = borrow  # final borrow-out == (a < b)
        s.add(flag_lt != If(ULT(a, b), ONE, ZERO))
    return prove(f"borrow-chain flag == (a <_u b), width {w}", setup)


for w in (1, 2, 3, 4, 8, 16, 32, 64):
    prove_cmp(w)

# ---- summary ----
n_ok = sum(1 for _, ok, _ in results if ok)
print(f"\n=== {n_ok}/{len(results)} lemmas PROVED (unsat on negation) ===")
sys.exit(0 if n_ok == len(results) else 1)
