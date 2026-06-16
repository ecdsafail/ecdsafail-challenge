//! Verified curve parameter sets for the per-curve EC point-addition port:
//! secp256k1, Curve25519 (birational short-Weierstrass form), SM2, and
//! Brainpool P256r1.
//!
//! Everything here is CLASSICAL `num_bigint::BigUint` — no quantum
//! circuits. It provides:
//!   - `CurveParams` + per-curve constructors with constants that are
//!     verified by the test module (on-curve generator, 2G on-curve,
//!     order*G == identity).
//!   - A self-contained short-Weierstrass affine implementation
//!     (`point_add`, `scalar_mul`) that HONORS the curve `a` parameter
//!     (we do not delegate to `zkp_ecc_lib`).
//!   - `ec_add_classical` — the generic affine add (R != +-P), the
//!     reference the quantum circuit is checked against.
//!   - `random_pair` — per-curve random on-curve (P, Q, R=P+Q) generation
//!     mirroring the secp256k1 circuit test's `rand_case`.
//!
//! Curve25519 is the Montgomery curve B*v^2 = u^3 + A*u^2 + u with
//! A = 486662, B = 1 over F_(2^255-19). Its birationally-equivalent
//! short-Weierstrass form y^2 = x^3 + a*x + b is obtained by the standard
//! map (u, v) -> (x, y):
//!   x = u/B + A/(3B),   y = v/B,
//!   a = (3 - A^2) / (3 B^2),
//!   b = (2 A^3 - 9 A B^2) / (27 B^3),
//! all mod p. We measure the group/DLOG attack cost, so a Weierstrass
//! affine add over the same field is the faithful port of the design.

use num_bigint::BigUint;
use num_traits::{One, Zero};
use rand::RngCore;

/// Per-curve modular-reduction strategy metadata. (Not used by the
/// classical reference here; consumed by the per-curve circuit code.)
#[derive(Clone, Debug, PartialEq)]
pub enum Reduction {
    /// q = 2^n - f (f small): pseudo-Mersenne fold.
    PseudoMersenne { f: u64 },
    /// NIST-style Solinas prime (sparse 2^32-word structure), e.g. SM2.
    Solinas,
    /// No special prime structure: generic shift/add + compare-subtract.
    Generic,
}

/// Verified short-Weierstrass curve parameters over `F_p`.
/// A point (x, y) is on the curve iff y^2 == x^3 + a*x + b (mod p).
#[derive(Clone, Debug)]
pub struct CurveParams {
    pub name: &'static str,
    pub p: BigUint,
    pub n: usize,
    pub a: BigUint,
    pub b: BigUint,
    pub gx: BigUint,
    pub gy: BigUint,
    pub order: BigUint,
    pub reduction: Reduction,
}

fn hx(s: &str) -> BigUint {
    BigUint::parse_bytes(s.as_bytes(), 16).expect("valid hex constant")
}

/// Modular inverse via Fermat: a^(p-2) mod p (p prime).
fn mod_inv(a: &BigUint, p: &BigUint) -> BigUint {
    a.modpow(&(p - BigUint::from(2u32)), p)
}

// --------------------------------------------------------------------------
// Curve constructors (constants verified in the test module below).
// --------------------------------------------------------------------------

/// secp256k1: p = 2^256 - 2^32 - 977, a = 0, b = 7.
#[must_use]
pub fn secp256k1() -> CurveParams {
    let p = (BigUint::from(1u32) << 256u32) - BigUint::from((1u64 << 32) + 977);
    CurveParams {
        name: "secp256k1",
        p,
        n: 256,
        a: BigUint::zero(),
        b: BigUint::from(7u32),
        gx: hx("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798"),
        gy: hx("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8"),
        order: hx("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141"),
        reduction: Reduction::PseudoMersenne {
            f: (1u64 << 32) + 977,
        },
    }
}

/// Curve25519 in its birationally-equivalent short-Weierstrass form over
/// F_(2^255-19). a, b, gx, gy are DERIVED from the Montgomery parameters
/// (A = 486662, B = 1, base u = 9) by the standard map (computed here, not
/// hard-coded). order = 8 * L (cofactor 8), L = 2^252 + 277...8493.
#[must_use]
pub fn curve25519() -> CurveParams {
    let p = (BigUint::from(1u32) << 255u32) - BigUint::from(19u32);
    let big_a = BigUint::from(486_662u32);
    let big_b = BigUint::one();

    // Weierstrass a = (3 - A^2) / (3 B^2) mod p. Compute (3 - A^2) mod p
    // additively (A^2 can exceed 3).
    let a_sq = (&big_a * &big_a) % &p;
    let three = BigUint::from(3u32);
    let num_a = (&three + &p - (&a_sq % &p)) % &p; // (3 - A^2) mod p
    let den_a = (&three * &big_b * &big_b) % &p; // 3 B^2
    let a_w = (&num_a * mod_inv(&den_a, &p)) % &p;

    // Weierstrass b = (2 A^3 - 9 A B^2) / (27 B^3) mod p.
    let a_cubed = (&a_sq * &big_a) % &p;
    let two_a3 = (BigUint::from(2u32) * &a_cubed) % &p;
    let nine_ab2 = (BigUint::from(9u32) * &big_a * &big_b * &big_b) % &p;
    let num_b = (&two_a3 + &p - &nine_ab2) % &p; // 2A^3 - 9AB^2 mod p
    let den_b = (BigUint::from(27u32) * &big_b * &big_b * &big_b) % &p; // 27 B^3
    let b_w = (&num_b * mod_inv(&den_b, &p)) % &p;

    // Map the Montgomery base point u = 9 (standard X25519 base). Recover
    // v from B v^2 = u^3 + A u^2 + u, i.e. v^2 = (u^3 + A u^2 + u)/B.
    let u = BigUint::from(9u32);
    let u2 = (&u * &u) % &p;
    let u3 = (&u2 * &u) % &p;
    let rhs = ((&u3 + (&big_a * &u2) % &p + &u) * mod_inv(&big_b, &p)) % &p;
    let v = sqrt_mod_5mod8(&rhs, &p);
    debug_assert_eq!((&v * &v) % &p, rhs, "Curve25519 v sqrt");

    // x_W = u/B + A/(3B), y_W = v/B (mod p).
    let inv_b = mod_inv(&big_b, &p);
    let inv_three_b = mod_inv(&((&three * &big_b) % &p), &p);
    let gx = ((&u * &inv_b) % &p + (&big_a * &inv_three_b) % &p) % &p;
    let gy = (&v * &inv_b) % &p;

    // order = 8 * L.
    let l = (BigUint::from(1u32) << 252u32)
        + BigUint::parse_bytes(b"27742317777372353535851937790883648493", 10).unwrap();
    let order = BigUint::from(8u32) * l;

    CurveParams {
        name: "curve25519",
        p,
        n: 255,
        a: a_w,
        b: b_w,
        gx,
        gy,
        order,
        reduction: Reduction::PseudoMersenne { f: 19 },
    }
}

/// SM2 (GM/T 0003 / GB/T 32918 recommended curve). a = p - 3.
#[must_use]
pub fn sm2() -> CurveParams {
    let p = hx("FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF");
    CurveParams {
        name: "sm2",
        a: hx("FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFC"),
        b: hx("28E9FA9E9D9F5E344D5A9E4BCF6509A7F39789F515AB8F92DDBCBD414D940E93"),
        gx: hx("32C4AE2C1F1981195F9904466A39C9948FE30BBFF2660BE1715A4589334C74C7"),
        gy: hx("BC3736A2F4F6779C59BDCEE36B692153D0A9877CC62A474002DF32E52139F0A0"),
        order: hx("FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123"),
        n: 256,
        p,
        reduction: Reduction::Solinas,
    }
}

/// Brainpool P256r1 (RFC 5639).
#[must_use]
pub fn brainpoolp256r1() -> CurveParams {
    CurveParams {
        name: "brainpoolp256r1",
        p: hx("A9FB57DBA1EEA9BC3E660A909D838D726E3BF623D52620282013481D1F6E5377"),
        a: hx("7D5A0975FC2C3057EEF67530417AFFE7FB8055C126DC5C6CE94A4B44F330B5D9"),
        b: hx("26DC5C6CE94A4B44F330B5D9BBD77CBF958416295CF7E1CE6BCCDC18FF8C07B6"),
        gx: hx("8BD2AEB9CB7E57CB2C4B482FFC81B7AFB9DE27E1E3BD23C23A4453BD9ACE3262"),
        gy: hx("547EF835C3DAC4FD97F8461A14611DC9C27745132DED8E545C1D54C72F046997"),
        order: hx("A9FB57DBA1EEA9BC3E660A909D838D718C397AA3B561A6F7901E0E82974856A7"),
        n: 256,
        reduction: Reduction::Generic,
    }
}

/// Square root mod p for p == 5 (mod 8) (Atkin's algorithm). Returns a
/// root r with r^2 == n (mod p) when n is a QR. Curve25519's p satisfies
/// p == 5 (mod 8); the generic curves never call this.
fn sqrt_mod_5mod8(n: &BigUint, p: &BigUint) -> BigUint {
    let n = n % p;
    // v = (2n)^((p-5)/8)
    let two_n = (BigUint::from(2u32) * &n) % p;
    let exp = (p - BigUint::from(5u32)) / BigUint::from(8u32);
    let v = two_n.modpow(&exp, p);
    // i = 2 n v^2  (i^2 == -1 mod p)
    let i = (BigUint::from(2u32) * &n % p * (&v * &v % p)) % p;
    // r = n v (i - 1)
    let i_minus_1 = (i + p - BigUint::one()) % p;
    (&n * &v % p * i_minus_1) % p
}

// --------------------------------------------------------------------------
// Self-contained short-Weierstrass affine arithmetic (honors `a` and `p`).
// Points are Option<(x, y)>; None is the point at infinity (identity).
// --------------------------------------------------------------------------

fn submod(a: &BigUint, b: &BigUint, p: &BigUint) -> BigUint {
    let a = a % p;
    let b = b % p;
    if a >= b {
        a - b
    } else {
        (a + p) - b
    }
}

/// Affine point addition (handles doubling and the point at infinity).
#[must_use]
pub fn point_add(
    p: &CurveParams,
    x1: &BigUint,
    y1: &BigUint,
    x2: &BigUint,
    y2: &BigUint,
) -> Option<(BigUint, BigUint)> {
    let m = &p.p;
    let lambda = if x1 == x2 {
        // Same x: either doubling or P + (-P) = identity.
        if (y1 + y2) % m == BigUint::zero() || y1.is_zero() {
            return None;
        }
        // lambda = (3 x1^2 + a) / (2 y1)
        let num = (BigUint::from(3u32) * x1 * x1 + &p.a) % m;
        let den = (BigUint::from(2u32) * y1) % m;
        (num * mod_inv(&den, m)) % m
    } else {
        // lambda = (y2 - y1) / (x2 - x1)
        let num = submod(y2, y1, m);
        let den = submod(x2, x1, m);
        (num * mod_inv(&den, m)) % m
    };
    let lam_sq = (&lambda * &lambda) % m;
    let x3 = submod(&submod(&lam_sq, x1, m), x2, m);
    let y3 = submod(&((&lambda * submod(x1, &x3, m)) % m), y1, m);
    Some((x3, y3))
}

/// Scalar multiplication k*(x, y) via double-and-add.
#[must_use]
pub fn scalar_mul(
    p: &CurveParams,
    k: &BigUint,
    x: &BigUint,
    y: &BigUint,
) -> Option<(BigUint, BigUint)> {
    let mut res: Option<(BigUint, BigUint)> = None;
    let mut base = (x.clone(), y.clone());
    let bits = k.bits();
    for i in 0..bits {
        if (k >> i) & BigUint::one() == BigUint::one() {
            res = match res {
                None => Some(base.clone()),
                Some((rx, ry)) => point_add(p, &rx, &ry, &base.0, &base.1),
            };
        }
        if i + 1 < bits {
            base = match point_add(p, &base.0, &base.1, &base.0, &base.1) {
                Some(pt) => pt,
                None => break, // doubling hit identity; no higher bits matter
            };
        }
    }
    res
}

/// Generate a random on-curve pair (P, Q) with distinct x, and R = P + Q.
/// P = `k_P` * G, Q = `k_Q` * G with `k_P`, `k_Q` in [1, order). Mirrors the
/// secp256k1 circuit test's `rand_case`. Returns (px, py, qx, qy, rx, ry).
pub fn random_pair(
    p: &CurveParams,
    rng: &mut dyn RngCore,
) -> (BigUint, BigUint, BigUint, BigUint, BigUint, BigUint) {
    let scalar = |rng: &mut dyn RngCore| -> BigUint {
        let mut bytes = [0u8; 32];
        rng.fill_bytes(&mut bytes);
        let mut k = BigUint::from_bytes_le(&bytes) % &p.order;
        if k.is_zero() {
            k = BigUint::one();
        }
        k
    };
    loop {
        let kp = scalar(rng);
        let kq = scalar(rng);
        if kp == kq {
            continue;
        }
        let Some(pp) = scalar_mul(p, &kp, &p.gx, &p.gy) else {
            continue;
        };
        let Some(qq) = scalar_mul(p, &kq, &p.gx, &p.gy) else {
            continue;
        };
        if pp.0 == qq.0 {
            // Same x => doubling or inverse; the generic-add reference
            // assumes R != +-P, so reject.
            continue;
        }
        let Some(r) = point_add(p, &pp.0, &pp.1, &qq.0, &qq.1) else {
            continue;
        };
        return (pp.0, pp.1, qq.0, qq.1, r.0, r.1);
    }
}

/// Generic affine point addition reference (the quantum circuit is checked
/// against this). Lifted verbatim from `ec_add_classical_secp256k1` with
/// the modulus `q` taken as a parameter. Assumes R != +-P (no doubling /
/// identity special cases).
#[must_use]
pub fn ec_add_classical(
    rx: &BigUint,
    ry: &BigUint,
    p_x: &BigUint,
    p_y: &BigUint,
    q: &BigUint,
) -> (BigUint, BigUint) {
    let dx: BigUint = if rx >= p_x { rx - p_x } else { (rx + q) - p_x };
    let dy: BigUint = if ry >= p_y { ry - p_y } else { (ry + q) - p_y };
    let dx_inv = dx.modpow(&(q - BigUint::from(2u32)), q);
    let lambda = (&dy * &dx_inv) % q;
    let lambda_sq = (&lambda * &lambda) % q;
    let x_new = if lambda_sq >= (rx + p_x) {
        (&lambda_sq - rx - p_x) % q
    } else {
        (&lambda_sq + q + q - rx - p_x) % q
    };
    let p_minus_xnew = if p_x >= &x_new {
        p_x - &x_new
    } else {
        (p_x + q) - &x_new
    };
    let lambda_pmxnew = (&lambda * &p_minus_xnew) % q;
    let y_new = if lambda_pmxnew >= *p_y {
        (&lambda_pmxnew - p_y) % q
    } else {
        (&lambda_pmxnew + q - p_y) % q
    };
    (x_new, y_new)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn all_curves() -> Vec<CurveParams> {
        vec![secp256k1(), curve25519(), sm2(), brainpoolp256r1()]
    }

    /// y^2 == x^3 + a*x + b (mod p).
    fn on_curve(c: &CurveParams, x: &BigUint, y: &BigUint) -> bool {
        let lhs = (y * y) % &c.p;
        let rhs = (((x * x % &c.p) * x) % &c.p + (&c.a * x) % &c.p + &c.b) % &c.p;
        lhs == rhs
    }

    #[test]
    fn generators_on_curve() {
        for c in all_curves() {
            assert!(
                on_curve(&c, &c.gx, &c.gy),
                "{}: generator not on curve",
                c.name
            );
        }
    }

    #[test]
    fn two_g_on_curve() {
        for c in all_curves() {
            let two_g = scalar_mul(&c, &BigUint::from(2u32), &c.gx, &c.gy)
                .unwrap_or_else(|| panic!("{}: 2G is identity", c.name));
            assert!(
                on_curve(&c, &two_g.0, &two_g.1),
                "{}: 2G not on curve",
                c.name
            );
        }
    }

    /// secp256k1 known 2G value (sanity that point_add/scalar_mul is right).
    #[test]
    fn secp256k1_two_g_known_vector() {
        let c = secp256k1();
        let two_g = scalar_mul(&c, &BigUint::from(2u32), &c.gx, &c.gy).unwrap();
        let g2x = hx("C6047F9441ED7D6D3045406E95C07CD85C778E4B8CEF3CA7ABAC09B95C709EE5");
        let g2y = hx("1AE168FEA63DC339A3C58419466CEAEEF7F632653266D0E1236431A950CFE52A");
        assert_eq!(two_g.0, g2x, "secp256k1 2G x");
        assert_eq!(two_g.1, g2y, "secp256k1 2G y");
    }

    /// order * G == identity for every curve.
    #[test]
    fn order_times_g_is_identity() {
        for c in all_curves() {
            assert!(
                scalar_mul(&c, &c.order, &c.gx, &c.gy).is_none(),
                "{}: order*G != identity",
                c.name
            );
        }
    }

    /// Curve25519: the DERIVED Weierstrass generator is on the Weierstrass
    /// curve, and order*G == identity (confirms the birational params).
    #[test]
    fn curve25519_birational_params() {
        let c = curve25519();
        assert!(
            on_curve(&c, &c.gx, &c.gy),
            "curve25519 Weierstrass generator not on curve"
        );
        assert!(
            scalar_mul(&c, &c.order, &c.gx, &c.gy).is_none(),
            "curve25519 order*G != identity"
        );
    }

    /// For every curve: 50 random pairs, the generic reference
    /// `ec_add_classical` agrees with the self-contained `point_add`, and
    /// both equal the scalar_mul-derived R from `random_pair`.
    #[test]
    fn ec_add_classical_matches_point_add_random() {
        let mut rng = rand::thread_rng();
        for c in all_curves() {
            for _ in 0..50 {
                let (px, py, qx, qy, rx, ry) = random_pair(&c, &mut rng);
                // Reference generic add (R != +-P, guaranteed by random_pair).
                let (ax, ay) = ec_add_classical(&px, &py, &qx, &qy, &c.p);
                // Self-contained affine add.
                let (bx, by) =
                    point_add(&c, &px, &py, &qx, &qy).expect("non-identity sum (distinct x)");
                assert_eq!(ax, bx, "{}: ec_add_classical x != point_add x", c.name);
                assert_eq!(ay, by, "{}: ec_add_classical y != point_add y", c.name);
                // Both equal the scalar_mul-derived R.
                assert_eq!(ax, rx, "{}: ec_add_classical x != random_pair R.x", c.name);
                assert_eq!(ay, ry, "{}: ec_add_classical y != random_pair R.y", c.name);
                // Sanity: R is on the curve.
                assert!(on_curve(&c, &rx, &ry), "{}: R not on curve", c.name);
            }
        }
    }
}
