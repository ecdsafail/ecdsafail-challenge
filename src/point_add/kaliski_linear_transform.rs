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
    step_linear_canonical_with_flags(st).0
}

fn step_linear_canonical_with_flags(st: &mut LinState) -> (Branch, u8, u8) {
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
    (br, a, m)
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
fn single_coefficient_pair_cannot_preserve_x_and_expose_quotient_by_constant_tag() {
    // Try the most tempting one-pair DIV rescue.  Set r0=ρ (nonzero constant)
    // so the lower output s=ρ*x preserves the denominator while seed
    // s0=y+β.  The upper output is
    //     r = k*y + (ρ*a + β*k).
    // If the parenthesized contaminant were a known constant, one coefficient
    // pair would simultaneously expose y/x and keep x, fitting the ~600q
    // target.  This requires an affine relation ρ*a(x)+β*k(x)=C across all x.
    // Three sampled transforms already make (a,k,1) non-collinear, killing all
    // constant-tag/constant-r0 variants of this rescue.
    let p = SECP256K1_P;
    let mut pts = Vec::new();
    for seed in 1..=3u64 {
        let x = random_element(seed);
        let seq = branch_sequence(x, ITERS);
        let (a, lower) = apply_coeffs(&seq, U256::from(1u64), U256::ZERO);
        let (k, zero) = apply_coeffs(&seq, U256::ZERO, U256::from(1u64));
        assert_eq!(lower, x);
        assert_eq!(zero, U256::ZERO);
        pts.push((a, k));
    }
    let (a0, k0) = pts[0];
    let (a1, k1) = pts[1];
    let (a2, k2) = pts[2];
    let da10 = sub_mod(a1, a0, p);
    let dk10 = sub_mod(k1, k0, p);
    let da20 = sub_mod(a2, a0, p);
    let dk20 = sub_mod(k2, k0, p);
    let det = sub_mod(da10.mul_mod(dk20, p), da20.mul_mod(dk10, p), p);
    eprintln!("constant-tag coefficient-pair relation determinant = {det:#x}");
    assert!(!det.is_zero(), "sampled (a,k) were affine-collinear; constant-tag DIV rescue may exist");
}

fn toy_branch_sequence_for_a_coeff(x: u64, p: u64, iters: usize) -> Vec<Branch> {
    let mut u = p;
    let mut v = x;
    let mut f = 1u8;
    let mut out = Vec::with_capacity(iters);
    for _ in 0..iters {
        let mut m = 0u8;
        if f == 1 && v == 0 { m ^= 1; }
        f ^= m;
        let u0 = (u & 1) as u8;
        let v0 = (v & 1) as u8;
        let mut a = 0u8;
        if f == 1 && u0 == 0 { a ^= 1; }
        if f == 1 && u0 == 1 && v0 == 0 { m ^= 1; }
        let b = a ^ m;
        let gt = if u > v { 1u8 } else { 0u8 };
        let delta = (f & gt) & (1 ^ b);
        a ^= delta;
        m ^= delta;
        let br = Branch { a_swap: a == 1, add: (f & (1 ^ b)) == 1 };
        out.push(br);
        if br.a_swap { core::mem::swap(&mut u, &mut v); }
        if br.add {
            assert!(v >= u, "Kaliski branch should subtract smaller from larger");
            v -= u;
        }
        v >>= 1;
        if br.a_swap { core::mem::swap(&mut u, &mut v); }
    }
    out
}

fn toy_apply_coeffs_for_a_coeff(seq: &[Branch], mut r: u64, mut s: u64, p: u64) -> (u64, u64) {
    for br in seq {
        if br.a_swap { core::mem::swap(&mut r, &mut s); }
        if br.add { s = (s + r) % p; }
        r = (2 * r) % p;
        if br.a_swap { core::mem::swap(&mut r, &mut s); }
    }
    (r, s)
}

fn toy_a_coefficient_phase_anf_stats(n: usize, p: u64, mask: u64) -> (usize, usize) {
    let size = 1usize << n;
    let mut anf = vec![0u8; size];
    for x in 0..size {
        let a = if x > 0 && (x as u64) < p {
            let seq = toy_branch_sequence_for_a_coeff(x as u64, p, 2 * n - 1);
            let (a, lower) = toy_apply_coeffs_for_a_coeff(&seq, 1, 0, p);
            assert_eq!(lower, x as u64);
            a
        } else {
            0
        };
        anf[x] = ((a & mask).count_ones() & 1) as u8;
    }
    for bit in 0..n {
        for idx in 0..size {
            if (idx & (1usize << bit)) != 0 {
                anf[idx] ^= anf[idx ^ (1usize << bit)];
            }
        }
    }
    let density = anf.iter().filter(|&&c| c != 0).count();
    let degree = anf
        .iter()
        .enumerate()
        .filter_map(|(i, &c)| if c != 0 { Some(i.count_ones() as usize) } else { None })
        .max()
        .unwrap_or(0);
    (degree, density)
}

fn toy_branch_history_phase_anf_stats(n: usize, p: u64) -> (usize, usize) {
    let size = 1usize << n;
    let mut anf = vec![0u8; size];
    for x in 0..size {
        let val = if x > 0 && (x as u64) < p {
            toy_branch_sequence_for_a_coeff(x as u64, p, 2 * n - 1)
                .iter()
                .enumerate()
                .filter(|(i, _)| i % 3 == 0)
                .fold(0u8, |acc, (_, br)| acc ^ (br.a_swap as u8))
        } else {
            0
        };
        anf[x] = val;
    }
    for bit in 0..n {
        for idx in 0..size {
            if (idx & (1usize << bit)) != 0 {
                anf[idx] ^= anf[idx ^ (1usize << bit)];
            }
        }
    }
    let density = anf.iter().filter(|&&c| c != 0).count();
    let degree = anf
        .iter()
        .enumerate()
        .filter_map(|(i, &c)| if c != 0 { Some(i.count_ones() as usize) } else { None })
        .max()
        .unwrap_or(0);
    (degree, density)
}

fn encode_branch_sequence_for_entropy_test(seq: &[Branch]) -> Vec<u8> {
    let mut out = vec![0u8; (2 * seq.len() + 7) / 8];
    for (i, br) in seq.iter().enumerate() {
        let val = (br.a_swap as u8) | ((br.add as u8) << 1);
        let bit = 2 * i;
        out[bit / 8] |= (val & 1) << (bit % 8);
        out[(bit + 1) / 8] |= ((val >> 1) & 1) << ((bit + 1) % 8);
    }
    out
}

#[test]
fn exact_branch_history_has_field_entropy_lower_bound() {
    // Even if we ignore circuit cost and just ask for a compressed exact
    // branch-history payload, the branch sequence appears to encode the whole
    // denominator.  Exhaustive toy fields are injective, and secp samples are
    // collision-free.  Therefore exact history compression cannot hope for
    // O(log n) or ~100 bits; it has a ~log2(p)≈256-bit information floor.
    use std::collections::BTreeSet;
    for &(n, p) in &[(4usize, 13u64), (6, 61), (8, 251)] {
        let mut seen = BTreeSet::new();
        for x in 1..p {
            let seq = toy_branch_sequence_for_a_coeff(x, p, 2 * n - 1);
            seen.insert(encode_branch_sequence_for_entropy_test(&seq));
        }
        eprintln!(
            "toy branch-history injectivity: n={n}, p={p}, distinct_sequences={} / {}",
            seen.len(),
            p - 1
        );
        assert_eq!(seen.len() as u64, p - 1);
    }

    let mut xs = BTreeSet::new();
    let mut seqs = BTreeSet::new();
    let mut seed = 1u64;
    while xs.len() < 4096 {
        let x = random_element(seed);
        seed += 1;
        if xs.insert(x) {
            let seq = branch_sequence(x, ITERS);
            seqs.insert(encode_branch_sequence_for_entropy_test(&seq));
        }
    }
    eprintln!(
        "secp branch-history sampled injectivity: distinct_sequences={} / {}",
        seqs.len(),
        xs.len()
    );
    assert_eq!(seqs.len(), xs.len());
}

#[test]
fn coefficient_transform_history_floor_misses_low_qubit_budget() {
    // Combine the 256-bit branch-history information floor with the minimum
    // live registers of the remaining coefficient-transform DIV shapes.  This
    // is a qubit lower bound, before flags, carries, comparators, or modular
    // arithmetic workspace.
    const N: usize = 256;
    const GOOGLE_LOW_QUBIT_TOTAL: usize = 1175;
    const DATA_REGS: usize = 2 * N;
    const SCRATCH_ALLOWANCE: usize = GOOGLE_LOW_QUBIT_TOTAL - DATA_REGS;
    let history_floor = N;
    let r_as_output_floor = N /* u */ + N /* r/output channel */ + history_floor;
    let second_channel_floor = N /* u */
        + N /* quotient/output channel */
        + 2 * N /* x-preserving coefficient pair leaves a(x),x */
        + history_floor;
    eprintln!(
        "coefficient-transform scratch floors: allowance={SCRATCH_ALLOWANCE}, r_as_output={r_as_output_floor}, second_channel={second_channel_floor}"
    );
    assert_eq!(SCRATCH_ALLOWANCE, 663);
    assert!(r_as_output_floor > SCRATCH_ALLOWANCE);
    assert!(second_channel_floor > SCRATCH_ALLOWANCE);
}

#[test]
fn initial_x_to_branch_history_oracle_is_dense_on_toy_kaliski() {
    // Compressing m/a history down to the initial denominator x is
    // information-theoretically possible (history is deterministic from x), but
    // it is not a cheap direct oracle.  A sparse parity of branch bits, viewed
    // as a function of the initial x, is already full-degree/dense on toy
    // Kaliski.  So on-the-fly branch regeneration from x is equivalent to
    // rerunning a Kaliski-like computation, not a tiny lookup/phase gadget.
    let cases = [(4usize, 13u64), (6, 61), (8, 251), (10, 1021), (12, 4093)];
    for &(n, p) in &cases {
        let (degree, density) = toy_branch_history_phase_anf_stats(n, p);
        let table = 1usize << n;
        eprintln!(
            "toy Kaliski branch-history phase from x: n={n}, p={p}, degree={degree}, density={density}/{table}"
        );
        assert!(degree + 1 >= n);
        assert!(density > table / 4);
    }
}

#[test]
fn a_coefficient_cancellation_is_dense_on_toy_kaliski() {
    // The constant-tag test above leaves one theoretical escape: preserve x in
    // the lower coefficient output and subtract the contaminant a(x) with a
    // data-dependent circuit.  On toy Kaliski transforms, mask bits of a(x) are
    // already full-degree and near-half-density ANFs.  So cancelling a(x) is not
    // a tiny phase/kickmix correction; it is effectively another Kaliski-like
    // branch computation.
    let cases = [
        (4usize, 13u64, 0b1010u64),
        (6usize, 61u64, 0b10_1010u64),
        (8usize, 251u64, 0b1010_0101u64),
        (10usize, 1021u64, 0b10_1001_0101u64),
        (12usize, 4093u64, 0b1010_0101_0101u64),
    ];
    for &(n, p, mask) in &cases {
        let (degree, density) = toy_a_coefficient_phase_anf_stats(n, p, mask);
        let table = 1usize << n;
        eprintln!(
            "toy Kaliski a(x) phase: n={n}, p={p}, degree={degree}, density={density}/{table}"
        );
        assert!(degree >= n - 1);
        assert!(density > table / 3);
    }
}

#[test]
fn dx_tagged_seed_recovers_division_with_negligible_exception() {
    // Approximate tolerance reopens the self-cleaning DIV route. Seed the
    // coefficient with (y + x) instead of y. Then
    //   T(x)*(0, y+x) = (k*y + k*x, 0) = (k*y - 2^ITERS, 0)
    // because k*x = -2^ITERS. Adding the known scale recovers k*y, and a
    // known rescale gives y/x. The only zero-coefficient exceptional set is
    // y = -x, probability ≈ 1/p for random field inputs.
    let p = SECP256K1_P;
    let scale = pow2_mod(ITERS);
    let scale_inv = scale.inv_mod(p).unwrap();
    for seed in 1..100u64 {
        let x = random_element(seed);
        let y = random_element(seed + 10_000);
        let tagged = add_mod(y, x, p);
        assert_ne!(tagged, U256::ZERO, "random sample hit y=-x exceptional set");
        let seq = branch_sequence(x, ITERS);
        let (r_tagged, s_out) = apply_coeffs(&seq, U256::ZERO, tagged);
        assert_eq!(s_out, U256::ZERO);
        let k_y = add_mod(r_tagged, scale, p); // r + 2^ITERS = k*y
        let quotient = neg_mod(k_y, p).mul_mod(scale_inv, p);
        assert_eq!(quotient, y.mul_mod(x.inv_mod(p).unwrap(), p));
    }
}

#[test]
fn stored_a_and_m_bits_recover_branch_pair() {
    // If we abandon qrisp's full inverse coefficient `(r,s)` sentinel, one
    // plausible branch-only cleanup stores the final swap bit `a` in addition
    // to the existing `m_hist`. The per-step add bit is then not independent:
    // for active steps, add = !(a xor m); after termination f=0 forces add=0.
    // This does not solve the 600-scratch target by itself (it still stores
    // history), but it validates the next branch-only circuit scaffold.
    for seed in 1..200u64 {
        let mut st = LinState {
            u: SECP256K1_P,
            v: random_element(seed),
            r: U256::ZERO,
            s: add_mod(random_element(seed + 10_000), random_element(seed), SECP256K1_P),
            f: 1,
        };
        for _ in 0..ITERS {
            let (br, a, m) = step_linear_canonical_with_flags(&mut st);
            assert_eq!(br.a_swap, a == 1);
            let recovered_add = st.f == 1 && ((a ^ m) == 0);
            assert_eq!(br.add, recovered_add, "add should be recoverable from stored a,m and post f");
        }
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

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct ToyLinKey {
    iter: usize,
    u: u64,
    v: u64,
    r: u64,
    s: u64,
    f: u8,
}

#[derive(Clone, Copy, Debug)]
struct ToyLinState {
    u: u64,
    v: u64,
    r: u64,
    s: u64,
    f: u8,
}

#[derive(Clone, Copy, Debug)]
struct ToyLinStateWithSidecar {
    u: u64,
    v: u64,
    r: u64,
    s: u64,
    tag_r: u64,
    tag_s: u64,
    f: u8,
}

fn toy_step_linear_canonical_with_sidecar(st: &mut ToyLinStateWithSidecar, p: u64) -> Branch {
    let mut m = 0u8;
    if st.f == 1 && st.v == 0 { m ^= 1; }
    st.f ^= m;
    let u0 = (st.u & 1) as u8;
    let v0 = (st.v & 1) as u8;
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
        core::mem::swap(&mut st.tag_r, &mut st.tag_s);
    }
    if br.add {
        assert!(st.v >= st.u);
        st.v -= st.u;
        st.s = (st.s + st.r) % p;
        st.tag_s = (st.tag_s + st.tag_r) % p;
    }
    st.v >>= 1;
    st.r = (2 * st.r) % p;
    st.tag_r = (2 * st.tag_r) % p;
    if br.a_swap {
        core::mem::swap(&mut st.u, &mut st.v);
        core::mem::swap(&mut st.r, &mut st.s);
        core::mem::swap(&mut st.tag_r, &mut st.tag_s);
    }
    br
}

#[derive(Clone, Copy, Debug)]
struct ToyUnreducedCoeffState {
    u: u64,
    v: u64,
    r: u128,
    s: u128,
    f: u8,
}

fn toy_step_unreduced_coeff(st: &mut ToyUnreducedCoeffState) -> Branch {
    let mut m = 0u8;
    if st.f == 1 && st.v == 0 { m ^= 1; }
    st.f ^= m;
    let u0 = (st.u & 1) as u8;
    let v0 = (st.v & 1) as u8;
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
        assert!(st.v >= st.u);
        st.v -= st.u;
        st.s += st.r;
    }
    st.v >>= 1;
    st.r <<= 1;
    if br.a_swap {
        core::mem::swap(&mut st.u, &mut st.v);
        core::mem::swap(&mut st.r, &mut st.s);
    }
    br
}

fn toy_unreduced_coeff_final_r(n: usize, p: u64, x: u64, y: u64) -> Option<u128> {
    if x == 0 || x >= p || y >= p {
        return None;
    }
    let tag = (x + y) % p;
    if tag == 0 {
        return None;
    }
    let mut st = ToyUnreducedCoeffState { u: p, v: x, r: 0, s: tag as u128, f: 1 };
    for _ in 0..(2 * n - 1) {
        toy_step_unreduced_coeff(&mut st);
    }
    Some(st.r)
}

fn toy_unreduced_coeff_highbit_phase_anf_stats(n: usize, p: u64, bit_shift: usize) -> (usize, usize) {
    assert!(n <= 10, "truth table kept small");
    let vars = 2 * n;
    let size = 1usize << vars;
    let limb_mask = (1u64 << n) - 1;
    let mut anf = vec![0u8; size];
    for idx in 0..size {
        let x = (idx as u64) & limb_mask;
        let y = ((idx >> n) as u64) & limb_mask;
        anf[idx] = toy_unreduced_coeff_final_r(n, p, x, y)
            .map(|r| ((r >> bit_shift) & 1) as u8)
            .unwrap_or(0);
    }
    for bit in 0..vars {
        for idx in 0..size {
            if (idx & (1usize << bit)) != 0 {
                anf[idx] ^= anf[idx ^ (1usize << bit)];
            }
        }
    }
    let density = anf.iter().filter(|&&c| c != 0).count();
    let degree = anf
        .iter()
        .enumerate()
        .filter_map(|(i, &c)| if c != 0 { Some(i.count_ones() as usize) } else { None })
        .max()
        .unwrap_or(0);
    (degree, density)
}

fn toy_unreduced_coeff_branch_stats(n: usize, p: u64) -> (usize, usize, usize, usize) {
    use std::collections::BTreeMap;
    type Key = (usize, u64, u64, u128, u128, u8);
    let mut seen: BTreeMap<Key, Branch> = BTreeMap::new();
    let mut conflicts = 0usize;
    let mut total = 0usize;
    let mut max_bits = 0usize;
    for x in 1..p {
        for y in 0..p {
            let tag = (x + y) % p;
            if tag == 0 { continue; }
            let mut st = ToyUnreducedCoeffState { u: p, v: x, r: 0, s: tag as u128, f: 1 };
            for iter in 0..(2 * n - 1) {
                let br = toy_step_unreduced_coeff(&mut st);
                max_bits = max_bits.max((128 - st.r.leading_zeros()) as usize);
                max_bits = max_bits.max((128 - st.s.leading_zeros()) as usize);
                let key = (iter, st.u, st.v, st.r, st.s, st.f);
                if let Some(prev) = seen.insert(key, br) {
                    if prev != br { conflicts += 1; }
                }
                total += 1;
            }
        }
    }
    (conflicts, total, seen.len(), max_bits)
}

fn toy_sidecar_branch_conflicts(n: usize, p: u64, sidecar_bits: usize) -> (usize, usize, usize) {
    use std::collections::BTreeMap;
    type Key = (usize, u64, u64, u64, u64, u8, u64, u64);
    let mask = if sidecar_bits >= 64 { u64::MAX } else { (1u64 << sidecar_bits) - 1 };
    let mut seen: BTreeMap<Key, Branch> = BTreeMap::new();
    let mut conflicts = 0usize;
    let mut total = 0usize;
    for x in 1..p {
        for y in 0..p {
            let tag = (x + y) % p;
            if tag == 0 { continue; }
            let mut st = ToyLinStateWithSidecar {
                u: p,
                v: x,
                r: 0,
                s: tag,
                // Independent known coefficient column.  If a small sidecar
                // could replace branch history, low bits of this evolving tag
                // should disambiguate the poststate branch.
                tag_r: 1,
                tag_s: 0,
                f: 1,
            };
            for iter in 0..(2 * n - 1) {
                let br = toy_step_linear_canonical_with_sidecar(&mut st, p);
                let key = (
                    iter,
                    st.u,
                    st.v,
                    st.r,
                    st.s,
                    st.f,
                    st.tag_r & mask,
                    st.tag_s & mask,
                );
                if let Some(prev) = seen.insert(key, br) {
                    if prev != br { conflicts += 1; }
                }
                total += 1;
            }
        }
    }
    (conflicts, total, seen.len())
}

fn toy_step_linear_canonical(st: &mut ToyLinState, p: u64) -> Branch {
    let mut m = 0u8;
    if st.f == 1 && st.v == 0 { m ^= 1; }
    st.f ^= m;
    let u0 = (st.u & 1) as u8;
    let v0 = (st.v & 1) as u8;
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
        assert!(st.v >= st.u);
        st.v -= st.u;
        st.s = (st.s + st.r) % p;
    }
    st.v >>= 1;
    st.r = (2 * st.r) % p;
    if br.a_swap {
        core::mem::swap(&mut st.u, &mut st.v);
        core::mem::swap(&mut st.r, &mut st.s);
    }
    br
}

#[test]
fn measuring_unreduced_coefficient_high_bits_has_dense_phase() {
    // Natural follow-up: keep unreduced coefficients just long enough to make
    // the Kaliski step locally reversible, then X-measure the high quotient bits
    // to get back under 600 scratch.  The phase correction would need those high
    // bits as boolean functions of the surviving data.  Even in the generous
    // input frame (x,y) the high bits are dense/full-degree on toys, so this is
    // not a cheap MBUC/kickmix escape from the quotient-bit width.
    let cases = [
        (6usize, 61u64, 10usize),
        (8usize, 251u64, 14usize),
        (10usize, 1021u64, 18usize),
    ];
    for &(n, p, bit_shift) in &cases {
        let (degree, density) = toy_unreduced_coeff_highbit_phase_anf_stats(n, p, bit_shift);
        let table = 1usize << (2 * n);
        eprintln!(
            "unreduced coefficient high-bit phase: n={n}, p={p}, bit={bit_shift}, degree={degree}, density={density}/{table}"
        );
        if n == 10 {
            println!("METRIC unreduced_coeff_highbit_degree_n10={degree}");
            println!("METRIC unreduced_coeff_highbit_density_n10={density}");
        }
        assert!(degree + 1 >= 2 * n);
        assert!(density > table / 4);
    }
}

#[test]
fn unreduced_coefficient_kaliski_self_cleans_but_width_kills_scratch600() {
    // First-principles check: the branch ambiguity is created by modularly
    // reducing the coefficient pair.  If we never reduce r,s modulo p, the high
    // quotient bits retain enough information to make the full poststate images
    // disjoint; no m_hist is needed on exhaustive toys.  But the price is that
    // the coefficient registers grow by the iteration count.  With ~407 secp
    // iterations, one data coefficient register would need about 256+407=663
    // bits.  In a folded 600-scratch layout this means scratch u (256) + wide r
    // (663) + extending input s by 407 bits = 1326 scratch, before carries.
    // This is the cleanest way to see where the missing branch history went:
    // it is sitting in the high quotient bits of unreduced coefficients.
    let cases = [(4usize, 13u64), (5, 31), (6, 61), (7, 127), (8, 251)];
    for &(n, p) in &cases {
        let (conflicts, total, states, max_bits) = toy_unreduced_coeff_branch_stats(n, p);
        eprintln!(
            "unreduced coefficient branch recovery: n={n}, p={p}, conflicts={conflicts}, states={states}, total={total}, max_coeff_bits={max_bits}"
        );
        assert_eq!(conflicts, 0);
        assert_eq!(max_bits, 3 * n - 1);
        if n == 8 {
            println!("METRIC unreduced_kaliski_max_coeff_bits_n8={max_bits}");
            println!("METRIC scratch600_unreduced_coeff_bits_secp=663");
            println!("METRIC scratch600_unreduced_scratch_floor=1326");
        }
    }
}

#[test]
fn scratch600_sidecar_tag_bits_do_not_fix_kaliski_branch_recovery() {
    // A 600-scratch DIV over the two 256-bit input registers has a brutal
    // accounting constraint.  Using tx as v and ty as the data coefficient s,
    // a Kaliski-like one-pair DIV already needs scratch for u and r: 2n=512
    // qubits.  Only 88 qubits remain for any branch-cleaning sidecar.  Could a
    // tiny independent tag channel disambiguate each iteration's branch from
    // the poststate, avoiding m_hist?  Exhaustive toy checks say no: even when
    // we evolve a full known coefficient column in parallel and reveal only its
    // low sidecar bits, the minimum exact disambiguating sidecar grows as n-1.
    // This is not a proof for all transforms, but it kills the hoped-for
    // "one coefficient pair + small tag" Kaliski layout under a 600q scratch cap.
    let cases = [(4usize, 13u64, 3usize), (5, 31, 4), (6, 61, 5), (7, 127, 6), (8, 251, 7)];
    let mut last_min = 0usize;
    for &(n, p, expected_min) in &cases {
        let mut min_bits = None;
        for bits in 0..=n {
            let (conflicts, total, states) = toy_sidecar_branch_conflicts(n, p, bits);
            eprintln!(
                "sidecar tag branch recovery: n={n}, p={p}, bits={bits}, conflicts={conflicts}, states={states}, total={total}"
            );
            if n == 8 && bits == 4 {
                println!("METRIC scratch600_sidecar_conflicts_n8_b4={conflicts}");
            }
            if conflicts == 0 {
                min_bits = Some(bits);
                break;
            }
        }
        let min_bits = min_bits.expect("full n-bit sidecar should disambiguate toy state");
        if n == 8 {
            println!("METRIC scratch600_sidecar_min_bits_n8={min_bits}");
            println!("METRIC scratch600_bare_kaliski_div_scratch=512");
            println!("METRIC scratch600_remaining_sidecar=88");
            println!("METRIC scratch600_extrapolated_sidecar_need=255");
        }
        assert_eq!(min_bits, expected_min);
        assert!(min_bits >= last_min);
        last_min = min_bits;
    }
}

#[test]
fn exhaustive_toy_full_poststate_does_not_recover_forward_branch() {
    // The secp256k1 sample above found no collisions when full coefficient
    // state was included, but that is not an information-theoretic guarantee.
    // Exhaustive tiny fields with the tagged nonzero seed s0=x+y show exact
    // collisions even when the reverse iteration index and full post-state are
    // known.  A forward-only self-cleaning Kaliski therefore needs more than
    // "inspect the live post-state"; it needs extra history, a different state
    // invariant, or a deliberately approximate exceptional set.
    use std::collections::BTreeMap;
    let cases = [(4usize, 13u64), (5usize, 31u64), (6usize, 61u64)];
    for &(n, p) in &cases {
        let mut seen: BTreeMap<ToyLinKey, Branch> = BTreeMap::new();
        let mut conflicts = 0usize;
        let mut total = 0usize;
        for x in 1..p {
            for y in 0..p {
                let tag = (x + y) % p;
                if tag == 0 { continue; }
                let mut st = ToyLinState { u: p, v: x, r: 0, s: tag, f: 1 };
                for iter in 0..(2 * n - 1) {
                    let br = toy_step_linear_canonical(&mut st, p);
                    let key = ToyLinKey { iter, u: st.u, v: st.v, r: st.r, s: st.s, f: st.f };
                    if let Some(prev) = seen.insert(key, br) {
                        if prev != br { conflicts += 1; }
                    }
                    total += 1;
                }
            }
        }
        eprintln!(
            "toy full-poststate branch recovery: n={n}, p={p}, total={total}, states={}, conflicts={conflicts}",
            seen.len()
        );
        assert!(conflicts > 0, "toy full post-state unexpectedly determined every branch");
    }
}

fn same_lin_state(a: &LinState, b: &LinState) -> bool {
    a.u == b.u && a.v == b.v && a.r == b.r && a.s == b.s && a.f == b.f
}

fn exact_local_predecessor_branch_count(post: LinState) -> usize {
    use std::collections::BTreeSet;
    let p = SECP256K1_P;
    let inv2 = U256::from(2u64).inv_mod(p).unwrap();
    let mut branches = BTreeSet::new();
    for a_swap in [false, true] {
        for add in [false, true] {
            let (u_after_shift, v_after_shift, r_after_double, s_after_add) = if a_swap {
                (post.v, post.u, post.s, post.r)
            } else {
                (post.u, post.v, post.r, post.s)
            };
            let r_before_double = r_after_double.mul_mod(inv2, p);
            let v_before_shift = v_after_shift << 1usize;
            let (v_before_add, s_before_add) = if add {
                (
                    v_before_shift.wrapping_add(u_after_shift),
                    sub_mod(s_after_add, r_before_double, p),
                )
            } else {
                (v_before_shift, s_after_add)
            };
            let (u0, v0, r0, s0) = if a_swap {
                (v_before_add, u_after_shift, s_before_add, r_before_double)
            } else {
                (u_after_shift, v_before_add, r_before_double, s_before_add)
            };
            for f0 in [0u8, 1u8] {
                let terminal_m = if f0 == 1 && v0 == U256::ZERO { 1u8 } else { 0u8 };
                if (f0 ^ terminal_m) != post.f {
                    continue;
                }
                let mut cand = LinState { u: u0, v: v0, r: r0, s: s0, f: f0 };
                let br = step_linear_canonical(&mut cand);
                if br == (Branch { a_swap, add }) && same_lin_state(&cand, &post) {
                    branches.insert((a_swap, add));
                }
            }
        }
    }
    branches.len()
}

#[test]
fn secp_local_poststate_predecessor_branch_is_ambiguous() {
    // Stronger than collision sampling: for each actually reached secp poststate,
    // enumerate all locally consistent inverse branches and re-run the step to
    // verify them.  Most tagged poststates still have multiple exact predecessor
    // branches.  Therefore there is no exact local poststate predicate hiding in
    // the arithmetic; branch cleanup needs history or a different transform.
    let mut hist = [0usize; 5];
    let mut ambiguous = 0usize;
    let mut total = 0usize;
    for seed in 1..=20u64 {
        let x = random_element(seed);
        let y = random_element(seed + 10_000);
        let mut st = LinState {
            u: SECP256K1_P,
            v: x,
            r: U256::ZERO,
            s: add_mod(x, y, SECP256K1_P),
            f: 1,
        };
        for _ in 0..ITERS {
            step_linear_canonical(&mut st);
            let count = exact_local_predecessor_branch_count(st);
            hist[count] += 1;
            if count > 1 {
                ambiguous += 1;
            }
            total += 1;
        }
    }
    let frac = ambiguous as f64 / total as f64;
    eprintln!(
        "secp exact local poststate predecessor branch counts: hist={hist:?}, ambiguous={ambiguous}/{total}, frac={frac:.6}"
    );
    assert!(frac > 0.60, "local poststate ambiguity unexpectedly rare: frac={frac}");
}

#[test]
fn tagged_full_poststate_branch_ambiguity_is_not_a_rare_exception() {
    // The approximate escape hatch would be to ignore the branch-recovery
    // collisions from `exhaustive_toy_full_poststate_does_not_recover_forward_branch`
    // as a tiny exceptional set.  They are not tiny on exhaustive toy fields.
    // With the nonzero tagged seed s0=x+y (excluding only y=-x), roughly a
    // quarter of step occurrences land in post-states that admit more than one
    // predecessor branch.  That is structural ambiguity, not a negligible
    // exceptional tail that can be patched cheaply.
    use std::collections::BTreeMap;
    let cases = [(4usize, 13u64), (5, 31), (6, 61), (7, 127), (8, 251)];
    for &(n, p) in &cases {
        let mut seen: BTreeMap<ToyLinKey, [usize; 4]> = BTreeMap::new();
        let mut total = 0usize;
        for x in 1..p {
            for y in 0..p {
                let tag = (x + y) % p;
                if tag == 0 { continue; }
                let mut st = ToyLinState { u: p, v: x, r: 0, s: tag, f: 1 };
                for iter in 0..(2 * n - 1) {
                    let br = toy_step_linear_canonical(&mut st, p);
                    let key = ToyLinKey { iter, u: st.u, v: st.v, r: st.r, s: st.s, f: st.f };
                    let idx = (br.a_swap as usize) * 2 + (br.add as usize);
                    seen.entry(key).or_insert([0; 4])[idx] += 1;
                    total += 1;
                }
            }
        }
        let mut ambiguous_keys = 0usize;
        let mut ambiguous_occurrences = 0usize;
        for counts in seen.values() {
            if counts.iter().filter(|&&c| c != 0).count() > 1 {
                ambiguous_keys += 1;
                ambiguous_occurrences += counts.iter().sum::<usize>();
            }
        }
        let frac = ambiguous_occurrences as f64 / total as f64;
        eprintln!(
            "toy tagged full-poststate branch ambiguity: n={n}, p={p}, total={total}, states={}, ambiguous_keys={ambiguous_keys}, ambiguous_occurrences={ambiguous_occurrences}, frac={frac:.6}",
            seen.len()
        );
        assert!(frac > 0.18, "branch ambiguity unexpectedly looked rare: frac={frac}");
    }
}

fn toy_additive_x_tagged_poststate_ambiguous_fraction<F>(n: usize, p: u64, tag_of_x: F) -> f64
where
    F: Fn(u64, u64) -> u64,
{
    use std::collections::BTreeMap;
    let mut seen: BTreeMap<ToyLinKey, [usize; 4]> = BTreeMap::new();
    let mut total = 0usize;
    for x in 1..p {
        let hx = tag_of_x(x, p) % p;
        for y in 0..p {
            let tag = (y + hx) % p;
            if tag == 0 { continue; }
            let mut st = ToyLinState { u: p, v: x, r: 0, s: tag, f: 1 };
            for iter in 0..(2 * n - 1) {
                let br = toy_step_linear_canonical(&mut st, p);
                let key = ToyLinKey { iter, u: st.u, v: st.v, r: st.r, s: st.s, f: st.f };
                let idx = (br.a_swap as usize) * 2 + (br.add as usize);
                seen.entry(key).or_insert([0; 4])[idx] += 1;
                total += 1;
            }
        }
    }
    let ambiguous_occurrences: usize = seen
        .values()
        .filter(|counts| counts.iter().filter(|&&c| c != 0).count() > 1)
        .map(|counts| counts.iter().sum::<usize>())
        .sum();
    ambiguous_occurrences as f64 / total as f64
}

fn toy_tagged_poststate_ambiguous_fraction(n: usize, p: u64, beta: u64) -> f64 {
    toy_additive_x_tagged_poststate_ambiguous_fraction(n, p, |x, p| (beta * x) % p)
}

#[test]
fn changing_linear_x_tag_does_not_fix_poststate_branch_ambiguity() {
    // Tagged DIV can seed s0 = y + beta*x for any known nonzero beta; the
    // output is still k*y plus a known constant because k*x is the fixed scale.
    // If poststate cleanup were only failing for the beta=1 tag, another beta
    // might rescue the idea.  Exhaustive toy fields show the ambiguity fraction
    // is invariant under nonzero beta: scaling/translating the coefficient
    // scalar does not add branch information.
    for &(n, p) in &[(4usize, 13u64), (6, 61), (8, 251)] {
        let base = toy_tagged_poststate_ambiguous_fraction(n, p, 1);
        for &beta in &[2u64, 7, 17] {
            if beta >= p { continue; }
            let got = toy_tagged_poststate_ambiguous_fraction(n, p, beta);
            eprintln!(
                "linear-x tag ambiguity: n={n}, p={p}, beta={beta}, frac={got:.6} (base={base:.6})"
            );
            assert!((got - base).abs() < 1e-12);
        }
    }
}

#[test]
fn any_x_only_additive_tag_is_only_a_bijection_not_branch_history() {
    // Generalize the linear retagging failure.  For fixed denominator x, any
    // tag of the form s0 = y + h(x) is just a bijective relabeling of the
    // coefficient scalar y (with one zero-scalar exception removed).  It cannot
    // encode the Kaliski branch history.  Exhaustive toy checks with nonlinear
    // h(x) confirm the ambiguity fraction is exactly unchanged.
    for &(n, p) in &[(4usize, 13u64), (6, 61), (8, 251)] {
        let base = toy_tagged_poststate_ambiguous_fraction(n, p, 1);
        let nonlinear_tags: &[fn(u64, u64) -> u64] = &[
            |x, p| (x * x + 3 * x + 5) % p,
            |x, p| (x * x % p * x + 7 * x + 11) % p,
        ];
        for (idx, tag) in nonlinear_tags.iter().enumerate() {
            let got = toy_additive_x_tagged_poststate_ambiguous_fraction(n, p, *tag);
            eprintln!(
                "x-only additive tag ambiguity: n={n}, p={p}, tag#{idx}, frac={got:.6} (base={base:.6})"
            );
            assert!((got - base).abs() < 1e-12);
            if n == 8 && idx == 1 {
                println!("METRIC x_only_additive_tag_ambiguity_frac_n8={got:.6}");
                println!("METRIC x_only_additive_tag_ambiguity_invariant=1");
            }
        }
    }
}

#[test]
fn bilinear_invariant_does_not_recover_inverse_branch() {
    // The obvious algebraic invariant of the coefficient transform is
    //     r*v + s*u = 0 (mod p)
    // starting from (u,v,r,s)=(p,x,0,tag). Unfortunately it is preserved by
    // almost all locally valid inverse candidates, so it does not provide the
    // cheap self-cleaning branch predicate we need.
    let p = SECP256K1_P;
    let inv2 = U256::from(2u64).inv_mod(p).unwrap();
    let mut ambiguous = 0usize;
    let mut total = 0usize;

    for seed in 1..=200u64 {
        let x = random_element(seed);
        let y = random_element(seed + 10_000);
        let mut st = LinState { u: p, v: x, r: U256::ZERO, s: add_mod(y, x, p), f: 1 };
        for _ in 0..ITERS {
            let br = step_linear_canonical(&mut st);
            if st.f == 0 { continue; }
            let mut survivors = 0usize;
            let candidates = [
                // (case_is_true, pre_u, pre_v, pre_r, pre_s)
                (!br.a_swap && !br.add, st.u, st.v << 1, st.r.mul_mod(inv2, p), st.s),
                ( br.a_swap && !br.add, st.u << 1usize, st.v, st.r, st.s.mul_mod(inv2, p)),
                ( br.a_swap &&  br.add, (st.u << 1usize).wrapping_add(st.v), st.v, sub_mod(st.r, st.s.mul_mod(inv2, p), p), st.s.mul_mod(inv2, p)),
                (!br.a_swap &&  br.add, st.u, (st.v << 1usize).wrapping_add(st.u), st.r.mul_mod(inv2, p), sub_mod(st.s, st.r.mul_mod(inv2, p), p)),
            ];
            for (_is_true, pu, pv, pr, ps) in candidates {
                let branch_valid = if pu.bit(0) == false {
                    // U-even candidate.
                    true
                } else if pv.bit(0) == false {
                    // V-even candidate.
                    true
                } else {
                    // Odd/odd candidate; either ordering is locally valid.
                    true
                };
                let invariant = add_mod(pr.mul_mod(pv % p, p), ps.mul_mod(pu % p, p), p) == U256::ZERO;
                if branch_valid && invariant {
                    survivors += 1;
                }
            }
            if survivors > 1 { ambiguous += 1; }
            total += 1;
        }
    }
    let frac = ambiguous as f64 / total as f64;
    assert!(frac > 0.90, "bilinear invariant unexpectedly disambiguated branches: ambiguous={frac}");
}

#[test]
fn low_bit_end_state_branch_classifier_is_not_approx_good_enough() {
    // Approximate incorrectness reopens rare exceptional sets, but it does not
    // make a crude local branch predicate viable. Train a best-majority lookup
    // table from low bits of the end-state registers, then test on disjoint
    // samples. Even with coefficient registers included, the error is huge.
    use std::collections::HashMap;

    type Key = (u16, u16, u16, u16, u8);
    const LOW_BITS: u32 = 3;
    let mask = (1u64 << LOW_BITS) - 1;
    let key_of = |st: &LinState| -> Key {
        (
            (st.u.as_limbs()[0] & mask) as u16,
            (st.v.as_limbs()[0] & mask) as u16,
            (st.r.as_limbs()[0] & mask) as u16,
            (st.s.as_limbs()[0] & mask) as u16,
            st.f,
        )
    };

    let mut counts: HashMap<Key, [usize; 4]> = HashMap::new();
    let idx = |br: Branch| -> usize { (br.a_swap as usize) * 2 + (br.add as usize) };

    for seed in 1..=120u64 {
        let mut st = LinState { u: SECP256K1_P, v: random_element(seed), r: U256::ZERO, s: random_element(seed + 10_000), f: 1 };
        for _ in 0..ITERS {
            let br = step_linear_canonical(&mut st);
            let k = key_of(&st);
            counts.entry(k).or_insert([0; 4])[idx(br)] += 1;
        }
    }

    let mut table: HashMap<Key, usize> = HashMap::new();
    for (k, c) in counts {
        let mut best_i = 0usize;
        let mut best_c = 0usize;
        for (i, &v) in c.iter().enumerate() {
            if v > best_c { best_c = v; best_i = i; }
        }
        table.insert(k, best_i);
    }

    let mut wrong = 0usize;
    let mut total = 0usize;
    for seed in 10_001..=10_120u64 {
        let mut st = LinState { u: SECP256K1_P, v: random_element(seed), r: U256::ZERO, s: random_element(seed + 10_000), f: 1 };
        for _ in 0..ITERS {
            let br = step_linear_canonical(&mut st);
            let k = key_of(&st);
            // All 3-bit keys are present in the train set; fallback is arbitrary.
            let pred = table.get(&k).copied().unwrap_or(0);
            if pred != idx(br) { wrong += 1; }
            total += 1;
        }
    }
    let err_rate = wrong as f64 / total as f64;
    assert!(err_rate > 0.50, "low-bit branch classifier unexpectedly good: err={err_rate}");
}

#[test]
fn zero_coefficient_seed_loses_branch_information() {
    // Exact DIV must also handle y=0 (or any value making the coefficient
    // channel uninformative). With r=s=0, full state collapses to the
    // denominator state, and branch recovery collides. Therefore any
    // self-cleaning forward-only Kaliski needs either an additional nonzero
    // tag mixed into the coefficient state or a branch predicate independent
    // of the coefficient scalar.
    use std::collections::HashMap;

    let mut seen: HashMap<([u64; 4], [u64; 4], [u64; 4], [u64; 4], u8), Branch> = HashMap::new();
    let mut conflicts = 0usize;
    for seed in 1..=200u64 {
        let mut st = LinState {
            u: SECP256K1_P,
            v: random_element(seed),
            r: U256::ZERO,
            s: U256::ZERO,
            f: 1,
        };
        for _ in 0..ITERS {
            let br = step_linear_canonical(&mut st);
            let key = (limbs(st.u), limbs(st.v), limbs(st.r), limbs(st.s), st.f);
            if let Some(prev) = seen.insert(key, br) {
                if prev != br { conflicts += 1; }
            }
        }
    }
    assert!(conflicts > 0, "zero coefficient seed unexpectedly preserved branch information");
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
