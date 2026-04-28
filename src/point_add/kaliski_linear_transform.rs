//! Ground-up structural probe: use Kaliski's coefficient update as a linear
//! transform on the *data* y-register instead of treating it as disposable
//! ancilla.
//!
//! This is analysis-only (`#[cfg(test)]` module imported from `mod.rs`). It
//! tests a possible 600-scratch architecture:
//!
//! - keep `tx = dx` as the preserved x-difference,
//! - use `ty` as Kaliski's coefficient register `s`, initialized to `dy`,
//! - run a canonical-mod-p coefficient version of Kaliski.
//!
//! If this worked naively, the forward Kaliski would turn `ty=dy` into
//! `s=0` while `r = raw_inv(dx) * dy`, i.e. the scaled slope. Then Kaliski's
//! backward coefficient transform might be used to write the final `Ry` into
//! `ty` without a second inversion. The tests below verify the linear algebra
//! and isolate the remaining obstruction.

#![cfg(test)]
#![allow(dead_code)]

use alloy_primitives::U256;
use sha3::{digest::{ExtendableOutput, Update, XofReader}, Shake128};

use super::SECP256K1_P;

const ITERS: usize = 407;

fn random_element(seed: u64) -> U256 {
    let mut h = Shake128::default();
    h.update(&seed.to_le_bytes());
    let mut reader = h.finalize_xof();
    loop {
        let mut buf = [0u8; 32];
        reader.read(&mut buf);
        let v = U256::from_be_bytes(buf);
        if v != U256::ZERO && v < SECP256K1_P {
            return v;
        }
    }
}

#[inline]
fn sub_mod(a: U256, b: U256, p: U256) -> U256 {
    let (r, borrow) = a.overflowing_sub(b);
    if borrow { r.wrapping_add(p) } else { r }
}

#[inline]
fn neg_mod(a: U256, p: U256) -> U256 {
    if a.is_zero() { a } else { p.wrapping_sub(a) }
}

#[inline]
fn add_mod(a: U256, b: U256, p: U256) -> U256 {
    a.add_mod(b, p)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Branch {
    a_swap: bool,
    add: bool,
}

#[derive(Clone, Copy, Debug)]
struct LinState {
    u: U256,
    v: U256,
    r: U256,
    s: U256,
    f: u8,
}

fn limbs(x: U256) -> [u64; 4] {
    *x.as_limbs()
}

/// The branch sequence depends only on `(u,v,f)`, not on the coefficient
/// values, so it can be separated from the coefficient linear transform.
fn branch_sequence(dx: U256, iters: usize) -> Vec<Branch> {
    let p = SECP256K1_P;
    let mut u = p;
    let mut v = dx;
    let mut f = 1u8;
    let mut out = Vec::with_capacity(iters);
    for _ in 0..iters {
        let mut m = 0u8;
        if f == 1 && v == U256::ZERO { m ^= 1; }
        f ^= m;

        let u0 = if u.bit(0) { 1u8 } else { 0u8 };
        let v0 = if v.bit(0) { 1u8 } else { 0u8 };
        let mut a = 0u8;
        if f == 1 && u0 == 0 { a ^= 1; }
        if f == 1 && u0 == 1 && v0 == 0 { m ^= 1; }
        let b = a ^ m;
        let gt = if u > v { 1u8 } else { 0u8 };
        let delta = (f & gt) & (1 ^ b);
        a ^= delta;
        m ^= delta;
        let add = (f & (1 ^ b)) == 1;
        let a_swap = a == 1;
        out.push(Branch { a_swap, add });

        if a_swap { core::mem::swap(&mut u, &mut v); }
        if add { v = v.wrapping_sub(u); }
        v >>= 1;
        if a_swap { core::mem::swap(&mut u, &mut v); }
        let _ = m;
    }
    out
}

/// Apply the coefficient-side transform with canonical mod-p arithmetic.
/// This is *not* exactly the current circuit's noncanonical `s=p` sentinel;
/// it is the modified architecture needed if `s` is a data register like `dy`.
fn apply_coeffs(seq: &[Branch], mut r: U256, mut s: U256) -> (U256, U256) {
    let p = SECP256K1_P;
    for br in seq {
        if br.a_swap { core::mem::swap(&mut r, &mut s); }
        if br.add { s = add_mod(s, r, p); }
        r = add_mod(r, r, p);
        if br.a_swap { core::mem::swap(&mut r, &mut s); }
    }
    (r, s)
}

fn pow2_mod(e: usize) -> U256 {
    let mut r = U256::from(1u64);
    for _ in 0..e {
        r = add_mod(r, r, SECP256K1_P);
    }
    r
}

fn step_linear_canonical(st: &mut LinState) -> Branch {
    let mut m = 0u8;
    if st.f == 1 && st.v == U256::ZERO { m ^= 1; }
    st.f ^= m;

    let u0 = if st.u.bit(0) { 1u8 } else { 0u8 };
    let v0 = if st.v.bit(0) { 1u8 } else { 0u8 };
    let mut a = 0u8;
    if st.f == 1 && u0 == 0 { a ^= 1; }
    if st.f == 1 && u0 == 1 && v0 == 0 { m ^= 1; }
    let b = a ^ m;
    let gt = if st.u > st.v { 1u8 } else { 0u8 };
    let delta = (st.f & gt) & (1 ^ b);
    a ^= delta;
    m ^= delta;
    let br = Branch { a_swap: a == 1, add: (st.f & (1 ^ b)) == 1 };

    if br.a_swap {
        core::mem::swap(&mut st.u, &mut st.v);
        core::mem::swap(&mut st.r, &mut st.s);
    }
    if br.add {
        st.v = st.v.wrapping_sub(st.u);
        st.s = add_mod(st.s, st.r, SECP256K1_P);
    }
    st.v >>= 1;
    st.r = add_mod(st.r, st.r, SECP256K1_P);
    if br.a_swap {
        core::mem::swap(&mut st.u, &mut st.v);
        core::mem::swap(&mut st.r, &mut st.s);
    }
    br
}

#[test]
fn coefficient_transform_shape() {
    let p = SECP256K1_P;
    let scale = pow2_mod(ITERS);
    for seed in 1..50u64 {
        let dx = random_element(seed);
        let seq = branch_sequence(dx, ITERS);
        let (a, c) = apply_coeffs(&seq, U256::from(1u64), U256::ZERO);
        let (k, d) = apply_coeffs(&seq, U256::ZERO, U256::from(1u64));

        // Empirical theorem for the canonical coefficient transform T(dx):
        //      T = [[a(dx), k(dx)], [dx, 0]]
        // with k(dx) * dx = -2^ITERS mod p.
        assert_eq!(c, dx, "lower-left coefficient is exactly dx");
        assert_eq!(d, U256::ZERO, "lower-right coefficient is zero");
        assert_eq!(k.mul_mod(dx, p), neg_mod(scale, p), "k is the raw inverse scale");
        assert_eq!(k.mul_mod(c, p), neg_mod(scale, p), "determinant relation");
        let _ = a;
    }
}

#[test]
fn dy_seeded_forward_computes_scaled_slope_and_zeroes_s() {
    let p = SECP256K1_P;
    let scale = pow2_mod(ITERS);
    for seed in 1..50u64 {
        let dx = random_element(seed);
        let dy = random_element(seed + 10_000);
        let seq = branch_sequence(dx, ITERS);
        let (r, s) = apply_coeffs(&seq, U256::ZERO, dy);
        let expect = neg_mod(scale, p)
            .mul_mod(dy, p)
            .mul_mod(dx.inv_mod(p).unwrap(), p);
        assert_eq!(r, expect, "r = raw_inv(dx) * dy = scaled slope");
        assert_eq!(s, U256::ZERO, "s/ty is consumed to zero in canonical form");
    }
}

#[test]
fn end_state_needs_coefficient_registers_to_recover_branch() {
    // A forward-only low-qubit DIV would like to run Kaliski without storing
    // m_hist. That requires each iteration's branch bit to be uncomputed from
    // the updated live state. This diagnostic separates two facts:
    //   1. denominator state alone (u,v,f) is NOT enough; many collisions occur.
    //   2. full coefficient state (u,v,r,s,f) WAS enough on this sample set.
    // So a self-cleaning DIV, if it exists, must use the coefficient registers
    // in the branch-recovery predicate; a tiny parity/comparator fingerprint is
    // not enough.
    use std::collections::HashMap;

    let mut denom_seen: HashMap<([u64; 4], [u64; 4], u8), Branch> = HashMap::new();
    let mut full_seen: HashMap<([u64; 4], [u64; 4], [u64; 4], [u64; 4], u8), Branch> = HashMap::new();
    let mut denom_conflicts = 0usize;
    let mut full_conflicts = 0usize;

    for seed in 1..=200u64 {
        let mut st = LinState {
            u: SECP256K1_P,
            v: random_element(seed),
            r: U256::ZERO,
            s: random_element(seed + 10_000),
            f: 1,
        };
        for _ in 0..ITERS {
            let br = step_linear_canonical(&mut st);
            let dk = (limbs(st.u), limbs(st.v), st.f);
            if let Some(prev) = denom_seen.insert(dk, br) {
                if prev != br { denom_conflicts += 1; }
            }
            let fk = (limbs(st.u), limbs(st.v), limbs(st.r), limbs(st.s), st.f);
            if let Some(prev) = full_seen.insert(fk, br) {
                if prev != br { full_conflicts += 1; }
            }
        }
    }

    assert!(denom_conflicts > 0, "denominator-only end-state unexpectedly recovered branches");
    assert_eq!(full_conflicts, 0, "full end-state branch recovery collided in samples");
}

#[test]
fn backward_write_condition_for_ry() {
    // If the coefficient transform is T=[[a,k],[dx,0]], then to have the
    // backward pass finish with `(r_initial=0, s_initial=Ry)`, the final
    // coefficient pair before backward MUST be T*(0,Ry) = (k*Ry, 0).
    // Starting from dy-seeded forward gives (k*dy, 0). So the structural
    // task is exactly to add k*(Ry-dy) into r, while s remains zero.
    // This test records the identity on random field values. It is not a
    // proof of impossibility; it is the crisp algebraic subproblem.
    let p = SECP256K1_P;
    for seed in 1..50u64 {
        let dx = random_element(seed);
        let dy = random_element(seed + 10_000);
        let ry = random_element(seed + 20_000);
        let seq = branch_sequence(dx, ITERS);
        let (k, _) = apply_coeffs(&seq, U256::ZERO, U256::from(1u64));
        let (r_dy, s_dy) = apply_coeffs(&seq, U256::ZERO, dy);
        let (r_ry, s_ry) = apply_coeffs(&seq, U256::ZERO, ry);
        assert_eq!(s_dy, U256::ZERO);
        assert_eq!(s_ry, U256::ZERO);
        assert_eq!(sub_mod(r_ry, r_dy, p), k.mul_mod(sub_mod(ry, dy, p), p));
    }
}
