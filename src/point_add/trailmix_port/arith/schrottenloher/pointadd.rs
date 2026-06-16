//! Schrottenloher 2026 Alg 1: reversible in-place EC point addition.
//!
//! Computes `(x2, y2) := R + P` in place, where `R = (x2, y2)` is the
//! QUANTUM point and `P = (ox, oy)` is the OTHER point, supplied in
//! CLASSICAL input registers (one value per shot — the fuzzer's runtime
//! controls). The five field-addend steps that touch P load `ox`/`oy`
//! into a fresh 257-bit temp (`add_creg`), apply the UNCONDITIONAL q-q
//! pseudo-Mersenne mod-add/sub, then unload the temp (`sub_creg`).
//!
//! Alg 1 (c=1 throughout for the non-windowed case):
//!   3:  x2 -= ox           (q-q mod-sub)
//!   4:  y2 -= oy           (q-q mod-sub)
//!   6:  y2 *= x2^-1        (`mod_mul_reverse` — Schrottenloher's
//!                           "inverse direction": the same circuit run
//!                           backwards does division)
//!   7-8: x2 += 3*ox        (3x q-q mod-add; leaves x2 = R.x + 2*ox)
//!   10: x2 -= y2^2         (controlled mod-square-sub; ctrl=|1>)
//!   11: y2 *= x2           (`mod_mul` forward)
//!   12: x2 := -x2          (controlled mod-negation; ctrl=|1>)
//!   14: y2 -= oy           (q-q mod-sub)
//!   15: x2 += ox           (q-q mod-add)
//!
//! Post: (x2, y2) holds R + P.

use crate::point_add::trailmix_port::arith::schrottenloher::gcd_pack::DIALOG_M;
use crate::point_add::trailmix_port::arith::schrottenloher::mod_mul_eea::{
    mod_mul_in_place_eea_secp256k1_m, mod_mul_in_place_eea_secp256k1_reverse_m,
    mod_mul_in_place_jump_reverse_secp256k1, mod_mul_in_place_jump_secp256k1,
};
use crate::point_add::trailmix_port::arith::schrottenloher::pm_prims::{
    controlled_mod_neg_pm_secp256k1, controlled_mod_square_sub_pm_secp256k1, mod_add_pm_secp256k1,
    mod_sub_pm_secp256k1,
};
use crate::point_add::trailmix_port::circuit::{Cbit, Circuit, QReg};
use num_bigint::BigUint;

/// Schrottenloher EC point addition: `(x2, y2) -> R + P` in place.
///
/// `x2_full`: caller-allocated `n + u_padding(n) = 293` bits. Pre: low
///   `n=256` bits hold `R.x`; high bits 0. Post: low 256 bits = (R+P).x
///   (mod q, non-canonical drift).
/// `y2_full`: caller-allocated `n + 1 = 257` bits. Pre: low `n` bits =
///   `R.y`; bit 256 = 0. Post: low 256 bits = (R+P).y.
/// `ox`, `oy`: 256-bit CLASSICAL input registers holding the other
///   point's coordinates in [0, q) (one value per shot).
///
/// Preconditions: R != ±P (no doubling/identity special cases).
pub fn ec_add_inplace_schrottenloher_secp256k1(
    circ: &mut Circuit,
    x2_full: &mut Vec<QReg>,
    y2_full: &[QReg],
    ox: &[Cbit],
    oy: &[Cbit],
) {
    ec_add_inplace_schrottenloher_secp256k1_m(circ, x2_full, y2_full, ox, oy, DIALOG_M, 0);
}

/// Window-size- and vent-parameterized EC point addition.
/// `dialog_m`: 5 = base-3-packed SP1-validated config (1173q / 2.477M tof);
/// 1 = RAW dialog (wider tape, no `compress_5iter`). `vents`: a shared
/// measurement-vent ancilla budget spent at the `apply_bv` peak phase to turn
/// ~2n controlled adds into cheaper vented adds (each vent -1 Toffoli /
/// +1 peak; the pool is reused across ops so peak rises by `vents`, not per
/// op). Configs: (5,0)=1173q/2.477M; (1,0)=1329q/2.272M; (1,90)~=1419q/~2.1M.
pub fn ec_add_inplace_schrottenloher_secp256k1_m(
    circ: &mut Circuit,
    x2_full: &mut Vec<QReg>,
    y2_full: &[QReg],
    ox: &[Cbit],
    oy: &[Cbit],
    dialog_m: usize,
    vents: usize,
) {
    let n = 256usize;
    assert!(x2_full.len() > n);
    assert_eq!(y2_full.len(), n + 1);
    assert_eq!(ox.len(), n, "ox is a 256-bit classical register");
    assert_eq!(oy.len(), n, "oy is a 256-bit classical register");

    let prev = circ.push_section("ec_add_schr");

    // Step 3: x2 -= ox.  => (x, y) = (R.x - P.x, R.y) = (dx, R.y).
    coord_addsub(circ, &x2_full[..=n], ox, true, 1);
    // Step 4: y2 -= oy.  => (dx, dy).
    coord_addsub(circ, y2_full, oy, true, 1);
    // Step 6: y2 *= x2^-1 (mod_mul reverse = division). => (dx, lambda).
    mod_mul_in_place_eea_secp256k1_reverse_m(circ, x2_full, y2_full, dialog_m, vents);
    // Step 7-8: x2 += 3*ox. => (R.x + 2*ox, lambda).
    coord_addsub(circ, &x2_full[..=n], ox, false, 3);
    // Step 10: x2 -= y2^2 (ctrl=|1>). => (R.x + 2*ox - lambda^2, lambda).
    let c = circ.alloc_qreg("alg1.c");
    circ.x(&c);
    controlled_mod_square_sub_pm_secp256k1(circ, &c, y2_full, &x2_full[..=n]);
    // Step 11: y2 *= x2.
    mod_mul_in_place_eea_secp256k1_m(circ, x2_full, y2_full, dialog_m, vents);
    // Step 12: x2 := -x2 (ctrl=|1>). => x = x_new - ox.
    controlled_mod_neg_pm_secp256k1(circ, &c, &x2_full[..=n]);
    circ.x(&c);
    circ.zero_and_free(c);
    // Step 14: y2 -= oy.
    coord_addsub(circ, y2_full, oy, true, 1);
    // Step 15: x2 += ox. => (x_new, y_new) = ec_add_classical(R, P).
    coord_addsub(circ, &x2_full[..=n], ox, false, 1);

    circ.pop_section(&prev);
}

/// Jump-GCD EC point-addition: identical Alg-1 driver to
/// `ec_add_inplace_schrottenloher_secp256k1_m` but the two in-place mod-muls
/// use the jump-before-swap dialog (fewer GCD/apply steps -> lower Toffoli).
/// Low-tof venting budget for the jump apply (coupled = materialize + vent the
/// +f window). Tuned to keep the coupled-apply peak at/just under the ~1420q
/// low-tof cap; raising it lowers Toffoli but raises the apply peak ~1:1.
pub const JUMP_LOWTOF_VENTS: usize = 256;

/// Jump=2 EC-add, LOW-QUBIT config (~1169q / ~2.08M tof): decoupled register
/// venting only (reg 12 fwd / 14 inv, +f borrowed-dirty -> peak-safe).
pub fn ec_add_inplace_schrottenloher_jump_secp256k1(
    circ: &mut Circuit,
    x2_full: &mut Vec<QReg>,
    y2_full: &[QReg],
    ox: &[Cbit],
    oy: &[Cbit],
    jump: usize,
) {
    ec_add_inplace_schrottenloher_jump_cfg(
        circ, x2_full, y2_full, ox, oy, jump, 12, 14, false, true,
    );
}

/// Jump=2 EC-add, LOW-TOF config (cap ~1420q): coupled venting (materialize +
/// vent the +f window, full register vents) spends the qubit headroom to drive
/// the Toffoli count below the low-qubit config.
pub fn ec_add_inplace_schrottenloher_jump_lowtof_secp256k1(
    circ: &mut Circuit,
    x2_full: &mut Vec<QReg>,
    y2_full: &[QReg],
    ox: &[Cbit],
    oy: &[Cbit],
    jump: usize,
) {
    ec_add_inplace_schrottenloher_jump_cfg(
        circ,
        x2_full,
        y2_full,
        ox,
        oy,
        jump,
        JUMP_LOWTOF_VENTS,
        JUMP_LOWTOF_VENTS,
        true,
        true,
    );
}

/// Shared jump EC-add driver. `fwd_vents`/`inv_vents` size the forward/inverse
/// apply register vents; `coupled` chooses borrowed-dirty (false) vs materialized
/// (true) +f folds. See the low-qubit / low-tof wrappers above.
#[allow(clippy::too_many_arguments)]
pub fn ec_add_inplace_schrottenloher_jump_cfg(
    circ: &mut Circuit,
    x2_full: &mut Vec<QReg>,
    y2_full: &[QReg],
    ox: &[Cbit],
    oy: &[Cbit],
    jump: usize,
    fwd_vents: usize,
    inv_vents: usize,
    coupled: bool,
    packed: bool,
) {
    let n = 256usize;
    assert!(x2_full.len() > n);
    assert_eq!(y2_full.len(), n + 1);
    assert_eq!(ox.len(), n);
    assert_eq!(oy.len(), n);

    let prev = circ.push_section("ec_add_schr_jump");
    coord_addsub(circ, &x2_full[..=n], ox, true, 1); // x2 -= ox
    coord_addsub(circ, y2_full, oy, true, 1); // y2 -= oy
    mod_mul_in_place_jump_reverse_secp256k1(
        circ, x2_full, y2_full, jump, inv_vents, coupled, packed,
    ); // y2 *= x2^-1
    coord_addsub(circ, &x2_full[..=n], ox, false, 3); // x2 += 3*ox
    let c = circ.alloc_qreg("alg1.c");
    circ.x(&c);
    controlled_mod_square_sub_pm_secp256k1(circ, &c, y2_full, &x2_full[..=n]); // x2 -= lambda^2
    mod_mul_in_place_jump_secp256k1(circ, x2_full, y2_full, jump, fwd_vents, coupled, packed); // y2 *= x2
    controlled_mod_neg_pm_secp256k1(circ, &c, &x2_full[..=n]); // x2 := -x2
    circ.x(&c);
    circ.zero_and_free(c);
    coord_addsub(circ, y2_full, oy, true, 1); // y2 -= oy
    coord_addsub(circ, &x2_full[..=n], ox, false, 1); // x2 += ox
    circ.pop_section(&prev);
}

/// `dst (257-bit, mod q) (+|-)= times * coord`, where `coord` is a 256-bit
/// CLASSICAL input register. Loads `coord` into a fresh 257-bit temp by
/// XOR-ing each classical bit into the (|0>-initialized) temp via
/// `x_if_bit` — a classically-controlled X, 0 Toffoli and self-inverse —
/// applies `times` UNCONDITIONAL q-q PM mod-add/sub, then unloads the temp
/// with the same XOR (clearing it to |0>) and frees it. (A creg ADD would
/// emit the full ripple-adder Toffolis even into a zero register; the XOR
/// load is free.) These addend steps run at low qubit occupancy (the
/// `mod_mul` workspace is freed around them), so the +257 temp does not
/// raise the `apply_bv` peak.
fn coord_addsub(circ: &mut Circuit, dst: &[QReg], coord: &[Cbit], subtract: bool, times: u32) {
    let n = 256usize;
    let temp = circ.alloc_qreg_bits("ec_coord_tmp", n + 1);
    // Load: temp[i] = coord[i] (temp starts |0>).
    for i in 0..n {
        circ.x_if_bit(&temp[i], coord[i]);
    }
    for _ in 0..times {
        if subtract {
            mod_sub_pm_secp256k1(circ, &temp, dst);
        } else {
            mod_add_pm_secp256k1(circ, &temp, dst);
        }
    }
    // Unload: same XOR clears temp back to |0>.
    for i in 0..n {
        circ.x_if_bit(&temp[i], coord[i]);
    }
    for q in temp {
        circ.zero_and_free(q);
    }
}

/// Classical reference for secp256k1 point addition. (x, y, `P_x`, `P_y`)
/// all in [0, q). Computes R + P in affine, assuming R != ±P.
#[must_use]
pub fn ec_add_classical_secp256k1(
    rx: &BigUint,
    ry: &BigUint,
    p_x: &BigUint,
    p_y: &BigUint,
) -> (BigUint, BigUint) {
    let q: BigUint = (BigUint::from(1u32) << 256u32) - BigUint::from(super::F_SECP256K1);

    let dx: BigUint = if rx >= p_x { rx - p_x } else { (rx + &q) - p_x };
    let dy: BigUint = if ry >= p_y { ry - p_y } else { (ry + &q) - p_y };
    let dx_inv = dx.modpow(&(&q - BigUint::from(2u32)), &q);
    let lambda = (&dy * &dx_inv) % &q;
    let lambda_sq = (&lambda * &lambda) % &q;
    let x_new = if lambda_sq >= (rx + p_x) {
        (&lambda_sq - rx - p_x) % &q
    } else {
        (&lambda_sq + &q + &q - rx - p_x) % &q
    };
    let p_minus_xnew = if p_x >= &x_new {
        p_x - &x_new
    } else {
        (p_x + &q) - &x_new
    };
    let lambda_pmxnew = (&lambda * &p_minus_xnew) % &q;
    let y_new = if &lambda_pmxnew >= p_y {
        (&lambda_pmxnew - p_y) % &q
    } else {
        (&lambda_pmxnew + &q - p_y) % &q
    };
    (x_new, y_new)
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigUint;
    use num_traits::Zero;

    /// secp256k1 generator point G (RFC reference).
    fn secp256k1_gen() -> (BigUint, BigUint) {
        let gx_hex = "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798";
        let gy_hex = "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8";
        (
            BigUint::parse_bytes(gx_hex.as_bytes(), 16).unwrap(),
            BigUint::parse_bytes(gy_hex.as_bytes(), 16).unwrap(),
        )
    }

    /// Verify the classical reference matches sage / hex point computation
    /// for G + G = 2G on secp256k1.
    #[test]
    fn classical_secp256k1_doubles_g_correctly() {
        // 2G on secp256k1 (well-known).
        let g2_x_hex = "C6047F9441ED7D6D3045406E95C07CD85C778E4B8CEF3CA7ABAC09B95C709EE5";
        let g2_y_hex = "1AE168FEA63DC339A3C58419466CEAEEF7F632653266D0E1236431A950CFE52A";
        let (gx, gy) = secp256k1_gen();
        let (g2x_want, g2y_want) = (
            BigUint::parse_bytes(g2_x_hex.as_bytes(), 16).unwrap(),
            BigUint::parse_bytes(g2_y_hex.as_bytes(), 16).unwrap(),
        );
        // Doubling differs from addition; this test only checks the
        // ADD formula on (G, P) where P != G.
        //
        // Instead, just verify the formula gives a point on the curve
        // for (G, 2G).
        let (sum_x, sum_y) = ec_add_classical_secp256k1(&gx, &gy, &g2x_want, &g2y_want);
        // 3G on secp256k1.
        let g3_x_hex = "F9308A019258C31049344F85F89D5229B531C845836F99B08601F113BCE036F9";
        let g3_y_hex = "388F7B0F632DE8140FE337E62A37F3566500A99934C2231B6CB9FD7584B8E672";
        let g3x = BigUint::parse_bytes(g3_x_hex.as_bytes(), 16).unwrap();
        let g3y = BigUint::parse_bytes(g3_y_hex.as_bytes(), 16).unwrap();
        assert_eq!(sum_x, g3x, "G + 2G x");
        assert_eq!(sum_y, g3y, "G + 2G y");
    }

    /// secp256k1 WeierstrassEllipticCurve (the zenodo reference curve),
    /// used to generate random on-curve point pairs for the e2e test.
    fn secp256k1_curve() -> zkp_ecc_lib::WeierstrassEllipticCurve {
        use alloy_primitives::U256;
        zkp_ecc_lib::WeierstrassEllipticCurve {
            modulus: U256::from_str_radix(
                "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F",
                16,
            )
            .unwrap(),
            a: U256::from(0u64),
            b: U256::from(7u64),
            gx: U256::from_str_radix(
                "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798",
                16,
            )
            .unwrap(),
            gy: U256::from_str_radix(
                "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8",
                16,
            )
            .unwrap(),
            order: U256::from_str_radix(
                "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141",
                16,
            )
            .unwrap(),
        }
    }

    /// End-to-end Schrottenloher EC point addition on 64 RANDOM on-curve
    /// pairs (P, Q) with distinct x: quantum P = (x2, y2), classical
    /// Q = (ox, oy) input registers (one value per shot — the fuzzer's
    /// runtime controls). Verifies (x2, y2) -> P + Q (mod q) on every
    /// shot for the given dialog window size. DEBUG_ON_FAIL drops into
    /// the time-travel debugger at the first failing shot.
    fn run_random_pairs(dialog_m: usize, vents: usize, max_peak: u32) {
        use crate::point_add::trailmix_port::arith::schrottenloher::gcd_pack::u_padding;
        use crate::point_add::trailmix_port::circuit::Cbit;
        use alloy_primitives::U256;
        use rand::RngCore;

        let n = 256usize;
        let total = n + u_padding(n);
        let q: BigUint = (BigUint::from(1u32) << 256u32) - BigUint::from(super::super::F_SECP256K1);

        let curve = secp256k1_curve();
        let mut rng = rand::thread_rng();
        let scalar = |rng: &mut dyn RngCore| -> U256 {
            U256::from(rng.next_u64())
                ^ (U256::from(rng.next_u64()) << 64)
                ^ (U256::from(rng.next_u64()) << 128)
                ^ (U256::from(rng.next_u64()) << 192)
        };
        // One random (P, Q, R = P + Q) with P.x != Q.x, distinct points.
        let rand_case = |rng: &mut dyn RngCore| loop {
            let s_p = scalar(rng);
            let s_q = scalar(rng);
            if s_p == U256::ZERO || s_q == U256::ZERO || s_p == s_q {
                continue;
            }
            let p = curve.mul(curve.gx, curve.gy, s_p);
            let qq = curve.mul(curve.gx, curve.gy, s_q);
            if p.0 == qq.0 {
                continue;
            }
            let r = curve.add(p.0, p.1, qq.0, qq.1);
            return (p.0, p.1, qq.0, qq.1, r.0, r.1);
        };

        let mut circ = crate::point_add::trailmix_port::circuit::Circuit::new();
        circ.set_max_qubit_peak(max_peak);
        let mut x2_full = circ.alloc_qreg_bits("x2", total);
        let y2_full = circ.alloc_qreg_bits("y2", n + 1);
        let ox: Vec<Cbit> = (0..n).map(|_| circ.alloc_input_bit()).collect();
        let oy: Vec<Cbit> = (0..n).map(|_| circ.alloc_input_bit()).collect();

        let u256_big = |v: U256| BigUint::from_bytes_le(&v.to_le_bytes::<32>());
        let mut want: Vec<(BigUint, BigUint)> = Vec::with_capacity(64);
        for shot in 0..64 {
            let (px, py, qx, qy, rx, ry) = rand_case(&mut rng);
            circ.sim_load_reg_bytes_shot(&x2_full[..n], &px.to_le_bytes::<32>(), shot);
            circ.sim_load_reg_bytes_shot(&y2_full[..n], &py.to_le_bytes::<32>(), shot);
            circ.sim_load_bits_bytes_shot(&ox, &qx.to_le_bytes::<32>(), shot);
            circ.sim_load_bits_bytes_shot(&oy, &qy.to_le_bytes::<32>(), shot);
            want.push((u256_big(rx), u256_big(ry)));
        }

        let ops0 = circ.ops.len();
        let ccx0 = circ.ccx_emitted;
        let ccz0 = circ.ccz_emitted;
        ec_add_inplace_schrottenloher_secp256k1_m(
            &mut circ,
            &mut x2_full,
            &y2_full,
            &ox,
            &oy,
            dialog_m,
            vents,
        );
        let ops = circ.ops.len() - ops0;
        let tof = (circ.ccx_emitted - ccx0) + (circ.ccz_emitted - ccz0);
        let peak = circ.peak_qubits;
        eprintln!(
            "  cost(ec_add_inplace_schrottenloher_secp256k1, dialog_m={dialog_m}, \
             vents={vents}, random pairs): ops={ops} tof={tof} peak+={peak}"
        );

        let mut outs: Vec<QReg> = Vec::new();
        outs.extend(x2_full);
        outs.extend(y2_full);
        let (sim, detached) = circ.destroy_sim(outs);

        for shot in 0..64 {
            let mut got_x = BigUint::zero();
            for i in 0..n {
                if sim.read_bit_shot(&detached[i], shot) == 1 {
                    got_x |= BigUint::from(1u32) << i;
                }
            }
            let mut got_y = BigUint::zero();
            for i in 0..n {
                if sim.read_bit_shot(&detached[total + i], shot) == 1 {
                    got_y |= BigUint::from(1u32) << i;
                }
            }
            let (want_x, want_y) = &want[shot];
            assert_eq!(&got_x % &q, *want_x, "shot {shot}: P+Q x");
            assert_eq!(&got_y % &q, *want_y, "shot {shot}: P+Q y");
        }
    }

    /// Low-qubit (base-3-packed, dialog_m=5, no vents) config: peak 1173,
    /// below Google's 1175. The SP1-validated headline circuit.
    #[test]
    fn ec_add_inplace_schrottenloher_secp256k1_random_pairs() {
        run_random_pairs(5, 0, 1178);
    }

    /// Higher-qubit / lower-Toffoli (raw, dialog_m=1, no vents) config: the
    /// dialog tape is stored unpacked (2 bits/pair), dropping every base-3
    /// compress_5iter (-204k Toffoli) at a wider tape (peak 1329, +156).
    #[test]
    fn ec_add_inplace_schrottenloher_secp256k1_random_pairs_raw() {
        run_random_pairs(1, 0, 1335);
    }

    /// Higher-qubit / lower-Toffoli RAW + apply_bv-peak venting (dialog_m=1,
    /// vents=90): spends a shared ancilla pool to turn the peak-phase
    /// controlled mod-add/sub (and ~63-bit +f reduction) from the Cuccaro
    /// ~3n path into measurement-vented adders. Measured 1451q / 2.126M tof
    /// (the +f reduction venting materializes the 63-bit constant, so the
    /// peak rises ~122 not ~90).
    #[test]
    fn ec_add_inplace_schrottenloher_secp256k1_random_pairs_raw_vented() {
        run_random_pairs(1, 90, 1455);
    }

    /// LOW-TOFFOLI config (~1416q / ~2.03M Toffoli, below the 1425q / 2.1M):
    /// M=3 SAT-permutation pack (cheap tape shrink 810->675, frees ~132
    /// qubits) + coupled apply_bv-peak venting (register add + materialized
    /// +f reduction). Measured ~1416q / ~2.025M tof.
    #[test]
    fn ec_add_inplace_schrottenloher_secp256k1_random_pairs_m3_vented() {
        run_random_pairs(3, 222, 1420);
    }

    /// Jump-GCD EC-add on 64 random on-curve pairs. Same harness as
    /// `run_random_pairs` but the jump-before-swap driver (fewer steps).
    fn run_random_pairs_jump(jump: usize, max_peak: u32, lowtof: bool) {
        use crate::point_add::trailmix_port::arith::schrottenloher::gcd_pack::u_padding;
        use crate::point_add::trailmix_port::circuit::Cbit;
        use alloy_primitives::U256;
        use rand::RngCore;

        let n = 256usize;
        let total = n + u_padding(n);
        let q: BigUint = (BigUint::from(1u32) << 256u32) - BigUint::from(super::super::F_SECP256K1);
        let curve = secp256k1_curve();
        let mut rng = rand::thread_rng();
        let scalar = |rng: &mut dyn RngCore| -> U256 {
            U256::from(rng.next_u64())
                ^ (U256::from(rng.next_u64()) << 64)
                ^ (U256::from(rng.next_u64()) << 128)
                ^ (U256::from(rng.next_u64()) << 192)
        };
        let rand_case = |rng: &mut dyn RngCore| loop {
            let s_p = scalar(rng);
            let s_q = scalar(rng);
            if s_p == U256::ZERO || s_q == U256::ZERO || s_p == s_q {
                continue;
            }
            let p = curve.mul(curve.gx, curve.gy, s_p);
            let qq = curve.mul(curve.gx, curve.gy, s_q);
            if p.0 == qq.0 {
                continue;
            }
            let r = curve.add(p.0, p.1, qq.0, qq.1);
            return (p.0, p.1, qq.0, qq.1, r.0, r.1);
        };

        let mut circ = crate::point_add::trailmix_port::circuit::Circuit::new();
        circ.set_max_qubit_peak(max_peak);
        let mut x2_full = circ.alloc_qreg_bits("x2", total);
        let y2_full = circ.alloc_qreg_bits("y2", n + 1);
        let ox: Vec<Cbit> = (0..n).map(|_| circ.alloc_input_bit()).collect();
        let oy: Vec<Cbit> = (0..n).map(|_| circ.alloc_input_bit()).collect();

        let u256_big = |v: U256| BigUint::from_bytes_le(&v.to_le_bytes::<32>());
        let mut want: Vec<(BigUint, BigUint)> = Vec::with_capacity(64);
        for shot in 0..64 {
            let (px, py, qx, qy, rx, ry) = rand_case(&mut rng);
            circ.sim_load_reg_bytes_shot(&x2_full[..n], &px.to_le_bytes::<32>(), shot);
            circ.sim_load_reg_bytes_shot(&y2_full[..n], &py.to_le_bytes::<32>(), shot);
            circ.sim_load_bits_bytes_shot(&ox, &qx.to_le_bytes::<32>(), shot);
            circ.sim_load_bits_bytes_shot(&oy, &qy.to_le_bytes::<32>(), shot);
            want.push((u256_big(rx), u256_big(ry)));
        }

        let ops0 = circ.ops.len();
        let ccx0 = circ.ccx_emitted;
        let ccz0 = circ.ccz_emitted;
        if lowtof {
            ec_add_inplace_schrottenloher_jump_lowtof_secp256k1(
                &mut circ,
                &mut x2_full,
                &y2_full,
                &ox,
                &oy,
                jump,
            );
        } else {
            // packed path is jump=2-specific (base-5); other jumps run unpacked
            ec_add_inplace_schrottenloher_jump_cfg(
                &mut circ,
                &mut x2_full,
                &y2_full,
                &ox,
                &oy,
                jump,
                12,
                14,
                false,
                jump == 2,
            );
        }
        let ops = circ.ops.len() - ops0;
        let tof = (circ.ccx_emitted - ccx0) + (circ.ccz_emitted - ccz0);
        let peak = circ.peak_qubits;
        eprintln!("  cost(ec_add jump={jump}, random pairs): ops={ops} tof={tof} peak+={peak}");

        let mut outs: Vec<QReg> = Vec::new();
        outs.extend(x2_full);
        outs.extend(y2_full);
        let (sim, detached) = circ.destroy_sim(outs);
        for shot in 0..64 {
            let mut got_x = BigUint::zero();
            for i in 0..n {
                if sim.read_bit_shot(&detached[i], shot) == 1 {
                    got_x |= BigUint::from(1u32) << i;
                }
            }
            let mut got_y = BigUint::zero();
            for i in 0..n {
                if sim.read_bit_shot(&detached[total + i], shot) == 1 {
                    got_y |= BigUint::from(1u32) << i;
                }
            }
            let (want_x, want_y) = &want[shot];
            assert_eq!(&got_x % &q, *want_x, "jump shot {shot}: P+Q x");
            assert_eq!(&got_y % &q, *want_y, "jump shot {shot}: P+Q y");
        }
    }

    /// Jump=2 EC-add (raw dialog). Peak is the unpacked tape (~814q) + 2 mul
    /// regs; base-5 packing brings it under 1175. tof should beat Google 2.1M.
    #[test]
    fn ec_add_inplace_schrottenloher_jump_random_pairs() {
        run_random_pairs_jump(2, 1200, false);
    }

    /// Jump=2 EC-add, LOW-TOF config (coupled venting, cap ~1420q). tof should
    /// beat the low-qubit config (~2.08M) AND the M=3 low-tof (~2.03M).
    #[test]
    fn ec_add_inplace_schrottenloher_jump_lowtof_random_pairs() {
        run_random_pairs_jump(2, 1420, true);
    }

    /// jump=3 EC-add, UNPACKED (no jump=3 packer yet). Measures the tof/peak of
    /// the higher-jump schedule (237 iters, sum(cn)=31706 vs jump=2's 273/36861)
    /// to estimate the packed product q*tof. Peak is the unpacked tape + 2 mul regs.
    #[test]
    fn ec_add_inplace_schrottenloher_jump3_random_pairs() {
        run_random_pairs_jump(3, 1600, false);
    }
}
