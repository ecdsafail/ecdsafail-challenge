//! Emit the LOW-TOFFOLI Schrottenloher EC point-add as a `.kmx` for the
//! zenodo fuzz tool. Config: `dialog_m = 3` (M=3 SAT-permutation pack, a
//! cheap tape shrink 810->675 freeing ~132 qubits) + `vents = 222`
//! (coupled apply_bv-peak venting: the 257-bit register adds AND the
//! materialized ~63-bit +f reductions). Measured ~1416q / ~2.025M Toffoli
//! -- below Google's low-tof config (1425q / 2.1M) on BOTH axes.
//!
//! Register layout matches the other emitters: reg0/1 = x2/y2 (256 qubits
//! each), reg2/3 = ox/oy (256 classical bits each).

use alloy_primitives::U256;
use rand::RngCore;
use std::io::Write;
use trailmix::arith::schrottenloher::gcd_pack::u_padding;
use trailmix::arith::schrottenloher::pointadd::ec_add_inplace_schrottenloher_secp256k1_m;
use trailmix::circuit::{Cbit, Circuit, QReg};
use zkp_ecc_lib::WeierstrassEllipticCurve;

const DIALOG_M: usize = 3;
const VENTS: usize = 222;

fn secp256k1() -> WeierstrassEllipticCurve {
    WeierstrassEllipticCurve {
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

fn rand_case(
    curve: &WeierstrassEllipticCurve,
    rng: &mut impl RngCore,
) -> (U256, U256, U256, U256, U256, U256) {
    loop {
        let scalar = |rng: &mut dyn RngCore| {
            U256::from(rng.next_u64())
                ^ (U256::from(rng.next_u64()) << 64)
                ^ (U256::from(rng.next_u64()) << 128)
                ^ (U256::from(rng.next_u64()) << 192)
        };
        let (s_p, s_q) = (scalar(rng), scalar(rng));
        if s_p == U256::ZERO || s_q == U256::ZERO || s_p == s_q {
            continue;
        }
        let p = curve.mul(curve.gx, curve.gy, s_p);
        let q = curve.mul(curve.gx, curve.gy, s_q);
        if p.0 == q.0 {
            continue;
        }
        let r = curve.add(p.0, p.1, q.0, q.1);
        return (p.0, p.1, q.0, q.1, r.0, r.1);
    }
}

fn main() {
    let n_cases: usize = std::env::var("N_CASES")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(64);
    let cases_out =
        std::env::var("CASES_OUT").unwrap_or_else(|_| "/tmp/ecs_cases_lowtof.txt".into());

    let curve = secp256k1();
    let mut rng = rand::thread_rng();
    let n = 256usize;
    let total = n + u_padding(n);

    let mut circ = Circuit::new();
    circ.set_max_qubit_peak(1420);
    let mut x2: Vec<QReg> = (0..total)
        .map(|i| circ.alloc_qreg(&format!("x2[{i}]")))
        .collect();
    let y2: Vec<QReg> = {
        let mut v: Vec<QReg> = (0..n)
            .map(|i| circ.alloc_qreg(&format!("y2[{i}]")))
            .collect();
        v.push(circ.alloc_qreg("y2_anc"));
        v
    };
    let ox: Vec<Cbit> = (0..n).map(|_| circ.alloc_input_bit()).collect();
    let oy: Vec<Cbit> = (0..n).map(|_| circ.alloc_input_bit()).collect();

    for shot in 0..64 {
        let (px, py, qx, qy, _rx, _ry) = rand_case(&curve, &mut rng);
        circ.sim_load_reg_bytes_shot(&x2[..n], &px.to_le_bytes::<32>(), shot);
        circ.sim_load_reg_bytes_shot(&y2[..n], &py.to_le_bytes::<32>(), shot);
        circ.sim_load_bits_bytes_shot(&ox, &qx.to_le_bytes::<32>(), shot);
        circ.sim_load_bits_bytes_shot(&oy, &qy.to_le_bytes::<32>(), shot);
    }

    eprintln!("[emit-lowtof] building (dialog_m={DIALOG_M}, vents={VENTS})...");
    ec_add_inplace_schrottenloher_secp256k1_m(&mut circ, &mut x2, &y2, &ox, &oy, DIALOG_M, VENTS);
    eprintln!(
        "[emit-lowtof] built: {} ops, peak {}q, {} tof",
        circ.total_ops(),
        circ.peak_qubits,
        circ.executed_toffoli_shots / 64
    );

    let mut all: Vec<QReg> = std::mem::take(&mut x2);
    all.extend(y2);
    let out = circ.defragment(all);

    circ.register(0);
    for q in &out[..n] {
        circ.append_qreg(q, 0);
    }
    circ.register(1);
    for q in &out[total..total + n] {
        circ.append_qreg(q, 1);
    }
    circ.register(2);
    for b in &ox {
        circ.append_bit(*b, 2);
    }
    circ.register(3);
    for b in &oy {
        circ.append_bit(*b, 3);
    }

    let mut f = std::fs::File::create(&cases_out).expect("create cases file");
    for _ in 0..n_cases {
        let (px, py, qx, qy, rx, ry) = rand_case(&curve, &mut rng);
        writeln!(f, "{px} {py} {qx} {qy} -> {rx} {ry} {qx} {qy}").unwrap();
    }
    eprintln!("[emit-lowtof] wrote {n_cases} cases -> {cases_out}");

    eprintln!("[emit-lowtof] serializing kmx...");
    print!("{}", circ.to_kmx());
    eprintln!("[emit-lowtof] done");

    let _ = circ.destroy_sim(out);
}
