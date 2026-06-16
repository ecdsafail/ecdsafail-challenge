//! Standalone cost probe for the dialog in-place multiply (`mod_mul_in_place_eea`)
//! per curve — measures the tof/peak of ONE multiply, to compare against the
//! Montgomery-mult prototype (which measured ~11.3M tof / ~1293 peak on a dense
//! random prime, ~1.8M on secp). This isolates the per-multiply dialog cost
//! that a Montgomery mult would replace for the driver's step-10/11 multiplies.

#[cfg(test)]
mod tests {
    use crate::point_add::trailmix_port::arith::schrottenloher::gcd_pack::u_padding;
    use crate::point_add::trailmix_port::ec::curves::dialog::mod_mul_in_place_eea;
    use crate::point_add::trailmix_port::ec::curves::mod_arith::{GenericPrime, ModRed, PseudoMersenne, Solinas};
    use crate::point_add::trailmix_port::ec::curves::params;
    use crate::point_add::trailmix_port::compat::num_bigint::BigUint;
    use rand::RngCore;

    fn measure<R: ModRed>(name: &str, mr: &R, p: &BigUint, n: usize) {
        let total = n + u_padding(n);
        let mut circ = crate::point_add::trailmix_port::circuit::Circuit::new();
        circ.set_max_qubit_peak(1300);
        let mut x_full = circ.alloc_qreg_bits("x", total);
        let y_full = circ.alloc_qreg_bits("y", n + 1);

        let mut rng = rand::thread_rng();
        let rand_lt = |rng: &mut dyn RngCore| -> BigUint {
            let mut b = [0u8; 32];
            loop {
                rng.fill_bytes(&mut b);
                let v = BigUint::from_bytes_le(&b);
                if &v < p {
                    return v;
                }
            }
        };
        let to32 = |v: &BigUint| {
            let mut b = [0u8; 32];
            for (i, x) in v.to_bytes_le().iter().take(32).enumerate() {
                b[i] = *x;
            }
            b
        };
        for shot in 0..64 {
            let x = rand_lt(&mut rng);
            let y = rand_lt(&mut rng);
            circ.sim_load_reg_bytes_shot(&x_full[..n], &to32(&x), shot);
            circ.sim_load_reg_bytes_shot(&y_full[..n], &to32(&y), shot);
        }

        let ops0 = circ.ops.len();
        let ccx0 = circ.ccx_emitted;
        let ccz0 = circ.ccz_emitted;
        mod_mul_in_place_eea(&mut circ, &mut x_full, &y_full, mr);
        let ops = circ.ops.len() - ops0;
        let tof = (circ.ccx_emitted - ccx0) + (circ.ccz_emitted - ccz0);
        let peak = circ.peak_qubits;
        eprintln!("  DIALOG-MUL cost [{name}]: ops={ops} tof={tof} peak+={peak}");

        // teardown must be clean (ancillae back to |0>); x restored, y = x*y mod p.
        let mut outs: Vec<crate::point_add::trailmix_port::circuit::QReg> = Vec::new();
        outs.extend(x_full);
        outs.extend(y_full);
        let _ = circ.destroy_sim(outs);
    }

    #[test]
    fn dialog_mul_cost_brainpool() {
        let cp = params::brainpoolp256r1();
        measure(
            "brainpool",
            &GenericPrime::new(256, cp.p.clone()),
            &cp.p,
            256,
        );
    }

    #[test]
    fn dialog_mul_cost_sm2() {
        let cp = params::sm2();
        measure("sm2", &Solinas::new(256, cp.p.clone()), &cp.p, 256);
    }

    #[test]
    fn dialog_mul_cost_secp256k1() {
        let cp = params::secp256k1();
        measure(
            "secp256k1",
            &PseudoMersenne::new(256, (1u64 << 32) + 977),
            &cp.p,
            256,
        );
    }
}
