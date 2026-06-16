//! Emit the shrunken-PZ EC point-add as a `.kmx` for the zenodo fuzz tool, with
//! the register layout the fuzzer expects:
//!   reg 0 = tx[0..256]  (QUANTUM)   -- input P.x, output R.x = (P+Q).x
//!   reg 1 = ty[0..256]  (QUANTUM)   -- input P.y, output R.y = (P+Q).y
//!   reg 2 = ox[0..256]  (CLASSICAL) -- input Q.x, unchanged
//!   reg 3 = oy[0..256]  (CLASSICAL) -- input Q.y, unchanged
//!
//! kmx -> stdout. A `P.x P.y Q.x Q.y -> R.x R.y Q.x Q.y` case file -> $CASES_OUT
//! (default /tmp/ec_shrunken_pz_cases.txt), N = $N_CASES (default 64) random secp
//! points (NO schedule prefilter -- the schedule must handle random values).
//!
//! The op stream is input-independent (ox/oy are runtime classical controls), so we
//! load one random 64-shot block only to satisfy the construction-time slope
//! contracts; the fuzz tool re-loads its own cases and re-validates.

use alloy_primitives::U256;
use rand::RngCore;
use std::io::Write;
use trailmix::circuit::{Cbit, Circuit, QReg};
use trailmix::ec::point_add::ec_add_inplace_shrunken_pz;
use zkp_ecc_lib::WeierstrassEllipticCurve;

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

/// One random (P, Q, R=P+Q) with P.x != Q.x and distinct points.
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
        let s_p = scalar(rng);
        let s_q = scalar(rng);
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
        std::env::var("CASES_OUT").unwrap_or_else(|_| "/tmp/ec_shrunken_pz_cases.txt".into());

    let curve = secp256k1();
    let mut rng = rand::thread_rng();
    let n = 256usize;

    let mut circ = Circuit::new();
    circ.set_max_qubit_peak(1300); // shrunken-PZ EC-add peak
    let mut tx: Vec<QReg> = (0..n)
        .map(|i| circ.alloc_qreg(&format!("tx[{i}]")))
        .collect();
    let mut ty: Vec<QReg> = (0..n)
        .map(|i| circ.alloc_qreg(&format!("ty[{i}]")))
        .collect();
    let ox: Vec<Cbit> = (0..n).map(|_| circ.alloc_input_bit()).collect();
    let oy: Vec<Cbit> = (0..n).map(|_| circ.alloc_input_bit()).collect();

    // Load one random 64-shot block so the construction-time slope contracts pass.
    for shot in 0..64 {
        let (px, py, qx, qy, _rx, _ry) = rand_case(&curve, &mut rng);
        circ.sim_load_reg_bytes_shot(&tx[..n], &px.to_le_bytes::<32>(), shot);
        circ.sim_load_reg_bytes_shot(&ty[..n], &py.to_le_bytes::<32>(), shot);
        circ.sim_load_bits_bytes_shot(&ox, &qx.to_le_bytes::<32>(), shot);
        circ.sim_load_bits_bytes_shot(&oy, &qy.to_le_bytes::<32>(), shot);
    }

    eprintln!("[emit] building ec_add_inplace_shrunken_pz...");
    ec_add_inplace_shrunken_pz(&mut circ, &mut tx, &mut ty, &ox, &oy);
    eprintln!(
        "[emit] built: {} ops, peak {}q, {} tof",
        circ.total_ops(),
        circ.peak_qubits,
        circ.executed_toffoli_shots / 64
    );

    // DEFRAGMENT: the in-place divide (gcd-pair resize + dy ghost/reconstruct)
    // migrates tx/ty's qubit ids to scattered high ids, so the registered output
    // ids would no longer match where the input was loaded at the start. The fuzzer
    // writes input AND reads output from the SAME ids, so restore tx/ty to their
    // canonical contiguous ids 0..2n first. out[..n] = tx, out[n..2n] = ty.
    let mut all: Vec<QReg> = std::mem::take(&mut tx);
    all.extend(std::mem::take(&mut ty));
    let out = circ.defragment(all);

    circ.register(0);
    for q in &out[..n] {
        circ.append_qreg(q, 0);
    }
    circ.register(1);
    for q in &out[n..2 * n] {
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

    // Case file: "P.x P.y Q.x Q.y -> R.x R.y Q.x Q.y" (decimal).
    let mut f = std::fs::File::create(&cases_out).expect("create cases file");
    for _ in 0..n_cases {
        let (px, py, qx, qy, rx, ry) = rand_case(&curve, &mut rng);
        writeln!(f, "{px} {py} {qx} {qy} -> {rx} {ry} {qx} {qy}").unwrap();
    }
    eprintln!("[emit] wrote {n_cases} cases -> {cases_out}");

    eprintln!("[emit] serializing kmx...");
    print!("{}", circ.to_kmx());
    eprintln!("[emit] done");

    let _ = circ.destroy_sim(out);
}
