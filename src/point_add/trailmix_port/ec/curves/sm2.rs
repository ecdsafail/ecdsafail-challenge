//! SM2 (short-Weierstrass over the SM2 prime, q = 2^256 - 2^224 - 2^96 +
//! 2^64 - 1) end-to-end EC-add test through the curve-generic driver. Uses
//! the Solinas (sparse structured-prime fold) reduction strategy: q = 2^256 -
//! `f_s` with `f_s` = 2^224 + 2^96 - 2^64 + 1 sparse, so the conditional-subtract
//! reduction folds the sparse `f_s` instead of the dense q.

#[cfg(test)]
mod tests {
    use crate::point_add::trailmix_port::arith::schrottenloher::gcd_pack::u_padding;
    use crate::point_add::trailmix_port::circuit::{Cbit, QReg};
    use crate::point_add::trailmix_port::ec::curves::driver::ec_add_inplace;
    use crate::point_add::trailmix_port::ec::curves::mod_arith::Solinas;
    use crate::point_add::trailmix_port::ec::curves::params;
    use crate::point_add::trailmix_port::compat::num_bigint::BigUint;
    use crate::point_add::trailmix_port::compat::num_traits::Zero;

    /// End-to-end SM2 EC point addition on 64 RANDOM on-curve pairs
    /// (P, Q) with distinct x: quantum P = (x2, y2), classical Q = (ox, oy)
    /// input registers (one value per shot). Verifies (x2, y2) -> P + Q
    /// (mod q) on every shot.
    #[test]
    fn ec_add_sm2_random_pairs() {
        let n = 256usize;
        let total = n + u_padding(n);
        let cp = params::sm2();
        let q = cp.p.clone();
        let mr = Solinas::new(256, cp.p.clone());

        let mut rng = rand::thread_rng();

        let mut circ = crate::point_add::trailmix_port::circuit::Circuit::new();
        // measured peak 1179 (Solinas sparse-f_s exact reduce; same peak as the
        // generic-prime reduce, 17.48M tof vs 24.32M); +6 margin.
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
        eprintln!("  cost(ec_add sm2, random pairs): ops={ops} tof={tof} peak+={peak}");

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
