//! Bisection probe: fuzz an individual schrottenloher primitive against
//! trailmix's OWN sim output (self-reference), so canonicalization /
//! PM-drift is factored out and any zenodo-vs-trailmix divergence is a
//! genuine kmx-fidelity bug in that primitive.
//!
//! WHICH=double|mul|mulrev|square|neg|add|sub (default double). The cases
//! file's "expected" column is what trailmix's sim produced, NOT a math
//! reference. So `pass` == zenodo agrees with trailmix on this op.
//!
//! kmx -> stdout; cases -> $CASES_OUT (default /tmp/schr_probe_cases.txt).

use crate::point_add::trailmix_port::compat::num_bigint::BigUint;
use rand::{thread_rng, Rng};
use std::io::Write;
use trailmix::arith::schrottenloher::gcd_pack::{
    forward_gcd_pack_quantum_secp256k1, forward_gcd_pack_quantum_secp256k1_reverse, u_padding,
};
use trailmix::arith::schrottenloher::mod_mul_eea::{
    mod_mul_in_place_eea_secp256k1, mod_mul_in_place_eea_secp256k1_reverse,
};
use trailmix::arith::schrottenloher::pm_prims::{
    controlled_mod_neg_pm_secp256k1, controlled_mod_square_sub_pm_secp256k1, mod_add_pm_secp256k1,
    mod_double_pm_secp256k1, mod_sub_pm_secp256k1,
};
use trailmix::circuit::{Circuit, QReg};

fn q() -> BigUint {
    (BigUint::from(1u32) << 256u32) - BigUint::from(2u32).pow(32) - BigUint::from(977u32)
}

fn load(circ: &mut Circuit, reg: &[QReg], v: &BigUint, shot: usize) {
    let mut b = [0u8; 32];
    for (i, x) in v.to_bytes_le().iter().take(32).enumerate() {
        b[i] = *x;
    }
    circ.sim_load_reg_bytes_shot(&reg[..256], &b, shot);
}

fn read256(sim: &trailmix::circuit::DestroyedSimState, reg: &[QReg], shot: usize) -> BigUint {
    let mut v = BigUint::from(0u32);
    for i in 0..256 {
        if sim.read_bit_shot(&reg[i], shot) == 1 {
            v |= BigUint::from(1u32) << i;
        }
    }
    v
}

fn main() {
    let which = std::env::var("WHICH").unwrap_or_else(|_| "double".into());
    let cases_out =
        std::env::var("CASES_OUT").unwrap_or_else(|_| "/tmp/schr_probe_cases.txt".into());
    let qq = q();
    let mut rng = thread_rng();

    let n = 256usize;
    let total = n + u_padding(n);

    let mut circ = Circuit::new();
    circ.set_max_qubit_peak(1300);

    // Two operands: `a` (wide, for mul workspace) and `b` (257).
    let mut a: Vec<QReg> = circ.alloc_qreg_bits("a", total);
    let mut b: Vec<QReg> = circ.alloc_qreg_bits("b", n + 1);

    // Record inputs per shot.
    let mut ins: Vec<(BigUint, BigUint)> = Vec::with_capacity(64);
    for shot in 0..64 {
        let av = if which == "compress" {
            // Valid radix-3 decompressed: 5 pairs, each in {00,01,11}.
            let mut dv: u32 = 0;
            for pair in 0..5 {
                let bits = match rng.gen_range(0..3u32) {
                    0 => 0b00u32,
                    1 => 0b01,
                    _ => 0b11,
                };
                dv |= bits << (2 * pair);
            }
            BigUint::from(dv)
        } else {
            BigUint::from_bytes_le(&rng.gen::<[u8; 32]>()) % &qq
        };
        let bv = BigUint::from_bytes_le(&rng.gen::<[u8; 32]>()) % &qq;
        load(&mut circ, &a, &av, shot);
        load(&mut circ, &b, &bv, shot);
        ins.push((av, bv));
    }

    eprintln!("[probe] WHICH={which}");
    match which.as_str() {
        "double" => mod_double_pm_secp256k1(&mut circ, &b),
        "add" => mod_add_pm_secp256k1(&mut circ, &a[..n + 1], &b),
        "sub" => mod_sub_pm_secp256k1(&mut circ, &a[..n + 1], &b),
        "compress" => {
            // gcd_compress5 round-trip on a[0..10] (valid radix-3 input).
            // Should restore a[0..10]; embeds compare_geq_const_gidney.
            use trailmix::arith::schrottenloher::gcd_compress5::{
                compress_5iter_refs, compress_5iter_reverse_refs,
            };
            let ds: Vec<QReg> = (0..8).map(|_| circ.alloc_qreg("c5d")).collect();
            {
                let w: [&QReg; 10] = std::array::from_fn(|k| &a[k]);
                let dr: Vec<&QReg> = ds.iter().collect();
                compress_5iter_refs(&mut circ, &w, &dr);
                compress_5iter_reverse_refs(&mut circ, &w, &dr);
            }
            for q in ds {
                circ.zero_and_free(q);
            }
        }
        "ltg" => {
            // Vented top-k less-than (the GCD's comparator), in isolation.
            // ctrl=a[8], target=a[9], compare a[0..8] vs b[0..8] on 8 bits.
            // a/b unchanged except a[9] ^= ctrl AND (a_top < b_top).
            trailmix::arith::schrottenloher::msb_compare::controlled_lt_msbs_gidney(
                &mut circ,
                &a[8],
                &a[0..8],
                &b[0..8],
                8,
                &a[9],
            );
        }
        "gcdrt" => {
            // GCD round-trip: forward pack (a -> garbage, a shrinks) then
            // reverse (garbage -> a restored). a should return to a_in; b
            // is untouched. Isolates the GCD pack/unpack from apply_bv.
            let mut garbage = forward_gcd_pack_quantum_secp256k1(&mut circ, &mut a);
            forward_gcd_pack_quantum_secp256k1_reverse(&mut circ, &mut a, &mut garbage);
            while a.len() < total {
                a.push(circ.alloc_qreg("a_pad_restore"));
            }
        }
        "mul" => mod_mul_in_place_eea_secp256k1(&mut circ, &mut a, &b),
        "mulrev" => mod_mul_in_place_eea_secp256k1_reverse(&mut circ, &mut a, &b),
        "square" => {
            let c = circ.alloc_qreg("c");
            circ.x(&c);
            controlled_mod_square_sub_pm_secp256k1(&mut circ, &c, &b, &a[..n + 1]);
            circ.x(&c);
            circ.zero_and_free(c);
        }
        "neg" => {
            let c = circ.alloc_qreg("c");
            circ.x(&c);
            controlled_mod_neg_pm_secp256k1(&mut circ, &c, &a[..n + 1]);
            circ.x(&c);
            circ.zero_and_free(c);
        }
        other => panic!("unknown WHICH={other}"),
    }
    eprintln!(
        "[probe] built: {} ops, peak {}q, {} tof",
        circ.total_ops(),
        circ.peak_qubits,
        circ.executed_toffoli_shots / 64
    );

    // DEFRAGMENT: internal free+realloc (GCD regrow) migrates a/b's
    // qubit ids; restore them to canonical contiguous ids so the
    // emitted register holds input AND output in the SAME ids (the
    // fuzzer's in-place convention). a (total) then b (n+1) -> ids
    // 0..total-1 and total..total+n.
    let mut all: Vec<QReg> = std::mem::take(&mut a);
    all.extend(std::mem::take(&mut b));
    let out = circ.defragment(all);

    // Register layout: reg0 = a[..256] (out[0..256]), reg1 = b[..256]
    // (out[total..total+256]).
    circ.register(0);
    for q in &out[..n] {
        circ.append_qreg(q, 0);
    }
    circ.register(1);
    for q in &out[total..total + n] {
        circ.append_qreg(q, 1);
    }

    let kmx = circ.to_kmx();

    // trailmix's own sim output -> self-reference expected.
    let (sim, detached) = circ.destroy_sim(out);
    let a_d = &detached[..total];
    let b_d = &detached[total..total + n + 1];

    let mut f = std::fs::File::create(&cases_out).expect("create cases");
    for shot in 0..64 {
        let ao = read256(&sim, a_d, shot);
        let bo = read256(&sim, b_d, shot);
        // input: a_in b_in   ->   output: a_out b_out  (trailmix sim)
        writeln!(f, "{} {} -> {} {}", ins[shot].0, ins[shot].1, ao, bo).unwrap();
    }
    eprintln!("[probe] wrote 64 self-reference cases -> {cases_out}");

    print!("{kmx}");
}
