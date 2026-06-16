//! Curve-generic Schrottenloher 2026 Alg 1: reversible in-place EC point
//! addition `(x2, y2) := R + P`.
//!
//! Curve-generic copy of
//! `arith::schrottenloher::pointadd::ec_add_inplace_schrottenloher_secp256k1`
//! (and its `coord_addsub`), routing the field-modular steps through the
//! `ModRed` trait. The affine add formula is curve-`a`-independent, so the
//! driver structure is shared across curves.
//!
//! Alg 1 (c=1 throughout for the non-windowed case):
//!   3:  x2 -= ox           (q-q mod-sub)
//!   4:  y2 -= oy           (q-q mod-sub)
//!   6:  y2 *= x2^-1        (`mod_mul_reverse` = division)
//!   7-8: x2 += 3*ox        (3x q-q mod-add)
//!   10: x2 -= y2^2         (controlled mod-square-sub; ctrl=|1>)
//!   11: y2 *= x2           (`mod_mul` forward)
//!   12: x2 := -x2          (controlled mod-negation; ctrl=|1>)
//!   14: y2 -= oy           (q-q mod-sub)
//!   15: x2 += ox           (q-q mod-add)

use crate::point_add::trailmix_port::circuit::{Cbit, Circuit, QReg};

use super::dialog::{mod_mul_in_place_eea, mod_mul_in_place_eea_reverse};
use super::mod_arith::ModRed;

/// Schrottenloher EC point addition: `(x2, y2) -> R + P` in place.
///
/// `x2_full`: caller-allocated `>= n + 1` bits (typically `n + u_padding(n)`).
///   Pre: low `n` bits hold `R.x`; high bits 0. Post: low n bits = (R+P).x
///   (mod q, non-canonical drift).
/// `y2_full`: caller-allocated `n + 1` bits. Pre: low `n` bits = `R.y`;
///   bit `n` = 0. Post: low n bits = (R+P).y.
/// `ox`, `oy`: `n`-bit CLASSICAL input registers holding the other point's
///   coordinates in [0, q) (one value per shot).
///
/// Preconditions: R != ±P (no doubling/identity special cases).
pub fn ec_add_inplace<R: ModRed>(
    circ: &mut Circuit,
    x2_full: &mut Vec<QReg>,
    y2_full: &[QReg],
    ox: &[Cbit],
    oy: &[Cbit],
    mr: &R,
) {
    let n = mr.n();
    assert!(x2_full.len() > n);
    assert_eq!(y2_full.len(), n + 1);
    assert_eq!(ox.len(), n, "ox is an n-bit classical register");
    assert_eq!(oy.len(), n, "oy is an n-bit classical register");

    let prev = circ.push_section("ec_add_schr");

    // Step 3: x2 -= ox.  => (dx, R.y).
    coord_addsub(circ, &x2_full[..=n], ox, true, 1, mr);
    // Step 4: y2 -= oy.  => (dx, dy).
    coord_addsub(circ, y2_full, oy, true, 1, mr);
    // Step 6: y2 *= x2^-1 (mod_mul reverse = division). => (dx, lambda).
    mod_mul_in_place_eea_reverse(circ, x2_full, y2_full, mr);
    // Step 7-8: x2 += 3*ox. => (R.x + 2*ox, lambda).
    coord_addsub(circ, &x2_full[..=n], ox, false, 3, mr);
    // Step 10: x2 -= y2^2 (ctrl=|1>). => (R.x + 2*ox - lambda^2, lambda).
    let c = circ.alloc_qreg("alg1.c");
    circ.x(&c);
    mr.controlled_mod_square_sub(circ, &c, y2_full, &x2_full[..=n]);
    // Step 11: y2 *= x2.
    mod_mul_in_place_eea(circ, x2_full, y2_full, mr);
    // Step 12: x2 := -x2 (ctrl=|1>). => x = x_new - ox.
    mr.controlled_mod_neg(circ, &c, &x2_full[..=n]);
    circ.x(&c);
    circ.zero_and_free(c);
    // Step 14: y2 -= oy.
    coord_addsub(circ, y2_full, oy, true, 1, mr);
    // Step 15: x2 += ox. => (x_new, y_new) = ec_add_classical(R, P).
    coord_addsub(circ, &x2_full[..=n], ox, false, 1, mr);

    circ.pop_section(&prev);
}

/// `dst (n+1-bit, mod q) (+|-)= times * coord`, where `coord` is an n-bit
/// CLASSICAL input register. Loads `coord` into a fresh `n+1`-bit temp via
/// `x_if_bit` (0 Toffoli, self-inverse), applies `times` UNCONDITIONAL q-q
/// PM mod-add/sub via `mr`, then unloads the temp and frees it.
///
/// Curve-generic copy of `pointadd::coord_addsub`.
fn coord_addsub<R: ModRed>(
    circ: &mut Circuit,
    dst: &[QReg],
    coord: &[Cbit],
    subtract: bool,
    times: u32,
    mr: &R,
) {
    let n = mr.n();
    let temp = circ.alloc_qreg_bits("ec_coord_tmp", n + 1);
    // Load: temp[i] = coord[i] (temp starts |0>).
    for i in 0..n {
        circ.x_if_bit(&temp[i], coord[i]);
    }
    for _ in 0..times {
        if subtract {
            mr.mod_sub_uncond(circ, &temp, dst);
        } else {
            mr.mod_add_uncond(circ, &temp, dst);
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::point_add::trailmix_port::ec::curves::mod_arith::PseudoMersenne;
    use crate::point_add::trailmix_port::ec::curves::params;
    use crate::point_add::trailmix_port::compat::num_bigint::BigUint;
    use crate::point_add::trailmix_port::compat::num_traits::Zero;

    /// Cross-check: the generic driver reproduces the secp256k1 baseline
    /// (~2.48M tof, peak ~1173). If the numbers differ wildly, the
    /// parameterization is wrong.
    #[test]
    fn ec_add_secp256k1_via_generic_driver_random_pairs() {
        use crate::point_add::trailmix_port::arith::schrottenloher::gcd_pack::u_padding;
        use crate::point_add::trailmix_port::circuit::Cbit;

        let n = 256usize;
        let total = n + u_padding(n);
        let cp = params::secp256k1();
        let q = cp.p.clone();
        let mr = PseudoMersenne::new(256, (1u64 << 32) + 977);

        let mut rng = rand::thread_rng();

        let mut circ = crate::point_add::trailmix_port::circuit::Circuit::new();
        circ.set_max_qubit_peak(1185);
        let mut x2_full = circ.alloc_qreg_bits("x2", total);
        let y2_full = circ.alloc_qreg_bits("y2", n + 1);
        let ox: Vec<Cbit> = (0..n).map(|_| circ.alloc_input_bit()).collect();
        let oy: Vec<Cbit> = (0..n).map(|_| circ.alloc_input_bit()).collect();

        let to32 = |v: &BigUint| {
            let mut b = [0u8; 32];
            for (i, x) in v.to_bytes_le().iter().take(32).enumerate() {
                b[i] = *x;
            }
            b
        };
        let mut want: Vec<(BigUint, BigUint)> = Vec::with_capacity(64);
        for shot in 0..64 {
            let (px, py, qx, qy, rx, ry) = params::random_pair(&cp, &mut rng);
            circ.sim_load_reg_bytes_shot(&x2_full[..n], &to32(&px), shot);
            circ.sim_load_reg_bytes_shot(&y2_full[..n], &to32(&py), shot);
            circ.sim_load_bits_bytes_shot(&ox, &to32(&qx), shot);
            circ.sim_load_bits_bytes_shot(&oy, &to32(&qy), shot);
            want.push((rx, ry));
        }

        let ops0 = circ.ops.len();
        let ccx0 = circ.ccx_emitted;
        let ccz0 = circ.ccz_emitted;
        ec_add_inplace(&mut circ, &mut x2_full, &y2_full, &ox, &oy, &mr);
        let ops = circ.ops.len() - ops0;
        let tof = (circ.ccx_emitted - ccx0) + (circ.ccz_emitted - ccz0);
        let peak = circ.peak_qubits;
        eprintln!(
            "  cost(ec_add secp256k1-via-generic, random pairs): ops={ops} tof={tof} peak+={peak}"
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
}
