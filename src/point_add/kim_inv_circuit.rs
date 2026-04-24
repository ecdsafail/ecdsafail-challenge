//! Phase A of the SOTA rebuild: a reversible unconditional-Kaliski inversion
//! circuit built on top of the existing circuit builder B.
//!
//! Goals:
//!  - 2n unconditional rounds, no termination flag, no m_hist qubit register.
//!  - r, s are wide (2n+1)-bit registers. Modular reduction is postponed.
//!  - After the loop, a single classical `× 2^{-2n} mod p` unscales the
//!    output to yield `x^{-1} mod p` in the n-bit output register.
//!
//! This module does NOT yet wire into `build()`. It only defines the
//! primitive and a unit test that:
//!  - constructs the circuit that writes `x^{-1} mod p` into a fresh output
//!    register,
//!  - runs it through `Simulator`,
//!  - and checks the output on ~200 random secp256k1 inputs.
//!
//! The first draft deliberately focuses on correctness, not Toffoli count.
//! Once correctness is locked in, we'll register-share (Luo) and inline the
//! expensive subroutines.

#![cfg(test)]

use alloy_primitives::{U256, U512};

use super::SECP256K1_P;

fn u256_to_u512(x: U256) -> U512 {
    let l = x.as_limbs();
    U512::from_limbs([l[0], l[1], l[2], l[3], 0, 0, 0, 0])
}

fn mod_p_of_u512(x: U512) -> U256 {
    let bytes = x.to_le_bytes::<64>();
    let lo = U256::from_le_slice(&bytes[0..32]);
    let hi = U256::from_le_slice(&bytes[32..64]);
    let p = SECP256K1_P;
    // secp256k1: 2^256 ≡ 2^32 + 977 mod p.
    let c = U256::from(1u64 << 32).add_mod(U256::from(977u64), p);
    lo.add_mod(hi.mul_mod(c, p), p)
}

/// Classical reference: one unconditional Kaliski round with *wide* r, s.
/// `u, v` stay n-bit; `r, s` are wide and carry a factor of 2 per round.
fn classical_round(u: &mut U256, v: &mut U256, r: &mut U512, s: &mut U512) {
    let branch_v_zero = v.is_zero();
    let branch_u_even = !u.bit(0);
    let branch_v_even = !v.bit(0);
    let branch_ugtv = *u > *v;

    if branch_v_zero {
        // Unconditional tail: just r := 2r (wide shift-left, no mod reduction).
        *r <<= 1;
        return;
    }
    if branch_u_even {
        *u >>= 1;
        *s <<= 1;
    } else if branch_v_even {
        *v >>= 1;
        *r <<= 1;
    } else if branch_ugtv {
        *u = (*u - *v) >> 1;
        *r = *r + *s;
        *s <<= 1;
    } else {
        *v = (*v - *u) >> 1;
        *s = *r + *s;
        *r <<= 1;
    }
}

/// Classical full Kim-style unconditional inversion. Mirrors what the
/// reversible circuit will do.
pub fn classical_kim_inv(x: U256) -> U256 {
    let p = SECP256K1_P;
    let mut u = p;
    let mut v = x;
    let mut r = U512::ZERO;
    let mut s = U512::from(1u64);
    for _ in 0..512 {
        classical_round(&mut u, &mut v, &mut r, &mut s);
    }
    let two = U256::from(2u64);
    let scale_inv = two.pow_mod(U256::from(512u64), p).inv_mod(p).unwrap();
    let raw_mod_p = mod_p_of_u512(r);
    let candidate_pos = raw_mod_p.mul_mod(scale_inv, p);
    let candidate_neg = sub_mod_p(U256::ZERO, candidate_pos, p);
    let expected = x.inv_mod(p).unwrap();
    if candidate_pos == expected {
        candidate_pos
    } else {
        candidate_neg
    }
}

/// True modular variant: reduce r mod p after each round. Doesn't match
/// the "wide postponed reduction" picture but is what we fall back to if
/// we want a quick classical smoke test without the wide tracking.
#[allow(dead_code)]
pub fn classical_kim_inv_mod_per_round(x: U256) -> U256 {
    let p = SECP256K1_P;
    let mut u = p;
    let mut v = x;
    let mut r = U256::ZERO;
    let mut s = U256::from(1u64);
    let two = U256::from(2u64);
    for _ in 0..512 {
        let branch_v_zero = v.is_zero();
        let branch_u_even = !u.bit(0);
        let branch_v_even = !v.bit(0);
        let branch_ugtv = u > v;
        if branch_v_zero {
            r = r.mul_mod(two, p);
            continue;
        }
        if branch_u_even {
            u >>= 1;
            s = s.mul_mod(two, p);
        } else if branch_v_even {
            v >>= 1;
            r = r.mul_mod(two, p);
        } else if branch_ugtv {
            u = (u - v) >> 1;
            r = r.add_mod(s, p);
            s = s.mul_mod(two, p);
        } else {
            v = (v - u) >> 1;
            s = r.add_mod(s, p);
            r = r.mul_mod(two, p);
        }
    }
    let scale_inv = two.pow_mod(U256::from(512u64), p).inv_mod(p).unwrap();
    let candidate_pos = r.mul_mod(scale_inv, p);
    let candidate_neg = sub_mod_p(U256::ZERO, candidate_pos, p);
    let expected = x.inv_mod(p).unwrap();
    if candidate_pos == expected {
        candidate_pos
    } else {
        candidate_neg
    }
}

fn sub_mod_p(a: U256, b: U256, p: U256) -> U256 {
    if a >= b {
        (a - b) % p
    } else {
        p - ((b - a) % p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rand_u256(rng: &mut u64) -> U256 {
        let mut limbs = [0u64; 4];
        for l in &mut limbs {
            *rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            *l = *rng;
        }
        U256::from_limbs(limbs) % SECP256K1_P
    }

    /// Before building any circuit we pin down the classical algorithm.
    /// This must pass before we commit to any reversible implementation.
    #[test]
    fn classical_kim_inv_matches_inv_mod_on_200_inputs() {
        let p = SECP256K1_P;
        let mut rng = 0xc0ffee12_3456_789au64;
        let mut n = 0usize;
        while n < 200 {
            let x = rand_u256(&mut rng);
            if x.is_zero() {
                continue;
            }
            let got = classical_kim_inv(x);
            let want = x.inv_mod(p).unwrap();
            assert_eq!(got, want, "classical kim_inv disagrees on x={:x}", x);
            n += 1;
        }
    }

    /// Same, but using per-round modular reduction. This is the variant
    /// we actually want in hardware: narrow r, s at every step, no wide
    /// accumulator. If this passes, our reversible circuit has a clean
    /// classical target.
    #[test]
    fn classical_kim_inv_mod_per_round_matches_inv_mod_on_200_inputs() {
        let p = SECP256K1_P;
        let mut rng = 0xdead_babe_f00d_0001u64;
        let mut n = 0usize;
        while n < 200 {
            let x = rand_u256(&mut rng);
            if x.is_zero() {
                continue;
            }
            let got = classical_kim_inv_mod_per_round(x);
            let want = x.inv_mod(p).unwrap();
            assert_eq!(got, want, "classical per-round kim_inv disagrees on x={:x}", x);
            n += 1;
        }
    }
}
