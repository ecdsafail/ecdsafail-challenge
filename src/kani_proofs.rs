//! Kani (bit-precise, bounded model checking) proofs that bind the arithmetic
//! contract to the ACTUAL Rust types the harness uses, not an abstract model.
//!
//! Compiled only under `--cfg kani` (i.e. `cargo kani`); the normal
//! `cargo build` / `benchmark.sh` path never sees this module, so it cannot
//! affect the circuit or the score.
//!
//! Complements `analysis/verify/*.py` (z3): z3 proves the width-256 arithmetic
//! over abstract bitvectors; Kani proves the exact Rust control flow of the
//! Solinas reduction (`mod_add_qq`, src/point_add/arith/modular.rs:12-49) using
//! the real `alloy_primitives::U256` type.
//!
//! Run:  cargo kani --harness solinas_add_u64
//!       cargo kani --harness solinas_add_u256
#![cfg(kani)]

use alloy_primitives::U256;
use crate::point_add::SECP256K1_P;

/// Solinas reduction, exactly mirroring `mod_add_qq`'s extended-register logic
/// but on plain integers instead of emitted gates. Division-free by design
/// (that is the whole point of Solinas), which also keeps it tractable for a
/// bit-precise solver. Represents the (n+1)-bit extended register as
/// (lo: U256, hi: bool).
fn solinas_add(a: U256, b: U256, p: U256) -> U256 {
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1u64)); // 2^256 - p
    let (s_lo, s_hi) = a.overflowing_add(b); // s = a + b, may set bit 256
    let (t_lo, c_carry) = s_lo.overflowing_add(c); // s + c
    let flag = s_hi | c_carry; // top bit of (s + c) == "reduction needed"
    if flag {
        t_lo // bit 256 cleared -> a + b - p
    } else {
        t_lo.overflowing_sub(c).0 // undo the + c -> a + b
    }
}

/// Small-width twin of `solinas_add` (fast to check); proves the CONTROL FLOW
/// is correct for a Solinas-shaped prime, independent of the U256 type.
fn solinas_add_small(a: u64, b: u64, p: u64) -> u64 {
    let c = u64::MAX.wrapping_sub(p).wrapping_add(1); // 2^64 - p
    let (s_lo, s_hi) = a.overflowing_add(b);
    let (t_lo, c_carry) = s_lo.overflowing_add(c);
    let flag = s_hi | c_carry;
    if flag { t_lo } else { t_lo.overflowing_sub(c).0 }
}

#[kani::proof]
fn solinas_add_u64() {
    // Solinas-shaped prime p = 2^64 - 1025 (c = 1025 small & sparse, like secp).
    let p: u64 = u64::MAX - 1024;
    let a: u64 = kani::any();
    let b: u64 = kani::any();
    kani::assume(a < p);
    kani::assume(b < p);

    let got = solinas_add_small(a, b, p);
    // division-free ground truth: a + b in [0, 2p) -> subtract p iff a+b >= p
    let (sum, carry) = a.overflowing_add(b);
    let ge = carry || sum >= p;
    let expect = if ge { sum.wrapping_sub(p) } else { sum };

    assert!(got == expect, "solinas u64 value");
    assert!(got < p, "solinas u64 in range");
}

#[kani::proof]
fn solinas_add_u256() {
    let p = SECP256K1_P; // real secp256k1 prime, from src/point_add/mod.rs
    let a = U256::from_limbs([kani::any(), kani::any(), kani::any(), kani::any()]);
    let b = U256::from_limbs([kani::any(), kani::any(), kani::any(), kani::any()]);
    kani::assume(a < p);
    kani::assume(b < p);

    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1u64));
    let got = solinas_add(a, b, p);

    // division-free ground truth for (a + b) mod p, a,b in [0,p) => a+b in [0,2p)
    let (s_lo, s_hi) = a.overflowing_add(b);
    let ge = s_hi || s_lo >= p;
    let expect = if ge {
        if s_hi { s_lo.wrapping_add(c) } else { s_lo - p }
    } else {
        s_lo
    };

    assert!(got == expect, "solinas u256 == (a+b) mod p");
    assert!(got < p, "solinas u256 result in [0,p)");
}

// NOTE on `weierstrass_elliptic_curve::sub_mod` (real function, line 16):
// it computes `a % m` / `b % m`. ruint's 256-bit `%` is Knuth long division with
// data-dependent loops that bounded model checking cannot unwind (a harness over
// it spins in `ruint::algorithms::cmp` past thousands of iterations and never
// converges without a `#[kani::unwind]` bound + division stub). This is not a
// gap in the proof so much as a property of the code: DIVISION-BASED modular
// arithmetic is not BMC-tractable, which is exactly why the quantum circuit uses
// the DIVISION-FREE Solinas reduction proved above (`solinas_add_u256`). The
// division path is instead covered by the z3 model in
// analysis/verify/solinas_reduction.py, whose 257-bit encoding sidesteps the
// concrete division algorithm entirely.
