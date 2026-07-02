//! Gating experiment for issue #5 (adder completeness).
//!
//! Runs the *built* point-add circuit on crafted EXCEPTIONAL inputs — cases the
//! scored fuzzer explicitly skips (`eval_circuit.rs:253` drops `dx=0`, and lines
//! 260-266 drop ∞) — and records whether the circuit leaves **clean ancilla +
//! correct phase**. This decides the completeness strategy (ADR 0006):
//!   - clean ancilla/phase, only the output is wrong  ⇒ Path A (negligibility) is
//!     viable: exceptions are tolerable at low amplitude.
//!   - dirty ancilla or phase garbage  ⇒ reversibility-breaking ⇒ Path B
//!     (complete formulas) or an explicit guard is required.
//!
//! `#[cfg(test)]` only — never compiled into the scored circuit.

use super::build;
use crate::circuit::{analyze_ops, QubitId, QubitOrBit};
use crate::sim::Simulator;
use crate::weierstrass_elliptic_curve::WeierstrassEllipticCurve;
use alloy_primitives::U256;
use sha3::{
    digest::{ExtendableOutput, Update},
    Shake128,
};

fn secp256k1() -> WeierstrassEllipticCurve {
    WeierstrassEllipticCurve {
        modulus: U256::from_str_radix(
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F",
            16,
        )
        .unwrap(),
        a: U256::from(0),
        b: U256::from(7),
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

/// Aggregate over `seeds` independent RNG draws (R/HMR randomness), since a
/// single seed's phase can be a coin-flip. Reports how many draws had a clean
/// phase / clean ancilla and whether the output ever matched the reference.
struct Agg {
    seeds: usize,
    phase_clean: usize,
    ancilla_clean: usize,
    out_matches_ref: usize,
}

fn run_case(
    ops: &[crate::circuit::Op],
    regs: &[Vec<QubitOrBit>],
    total_qubits: u64,
    num_bits: u64,
    seeds: usize,
    px: U256,
    py: U256,
    qx: U256,
    qy: U256,
    expected: (U256, U256),
) -> Agg {
    let mut reg_qubits = std::collections::HashSet::new();
    for r in regs {
        for qb in r {
            if let QubitOrBit::Qubit(q) = qb {
                reg_qubits.insert(q.0);
            }
        }
    }
    let mut agg = Agg { seeds, phase_clean: 0, ancilla_clean: 0, out_matches_ref: 0 };
    for s in 0..seeds {
        let mut hasher = Shake128::default();
        hasher.update(b"completeness-probe-v2");
        hasher.update(&(s as u64).to_le_bytes());
        let mut xof = hasher.finalize_xof();
        let mut sim = Simulator::new(total_qubits as usize, num_bits as usize, &mut xof);
        sim.clear_for_shot();
        sim.set_register(&regs[0], px, 0);
        sim.set_register(&regs[1], py, 0);
        sim.set_register(&regs[2], qx, 0);
        sim.set_register(&regs[3], qy, 0);
        sim.apply_iter(ops.iter());
        let out = (sim.get_register(&regs[0], 0), sim.get_register(&regs[1], 0));
        if (sim.phase >> 0) & 1 == 0 {
            agg.phase_clean += 1;
        }
        let dirty = (0..total_qubits)
            .filter(|q| !reg_qubits.contains(q))
            .any(|q| (sim.qubit(QubitId(q)) >> 0) & 1 != 0);
        if !dirty {
            agg.ancilla_clean += 1;
        }
        if out == expected {
            agg.out_matches_ref += 1;
        }
    }
    let _ = num_bits;
    agg
}

// Heavy: builds the full circuit and simulates it across 16 seeds × 3 cases.
// Opt-in so it never slows the default `cargo test`; run with `--ignored`.
#[test]
#[ignore = "heavy gating probe; run explicitly with `cargo test -- --ignored`"]
fn completeness_probe_exceptional_inputs() {
    let c = secp256k1();
    let g = (c.gx, c.gy);
    let two_g = c.add(g.0, g.1, g.0, g.1);
    let neg_g = (c.gx, sub_mod_p(U256::ZERO, c.gy, c.modulus));
    let q = c.mul(c.gx, c.gy, U256::from(7u64));
    let gen_exp = c.add(g.0, g.1, q.0, q.1);

    // Build the scored circuit once; reuse for every case/seed.
    let ops = build();
    let (tq, nb, _n, regs) = analyze_ops(ops.iter());
    assert_eq!(regs.len(), 4);
    let run = |seeds, px, py, qx, qy, exp| {
        run_case(&ops, &regs, tq, nb, seeds, px, py, qx, qy, exp)
    };

    const K: usize = 16;
    let gen = run(1, g.0, g.1, q.0, q.1, gen_exp); // generic control
    let dbl = run(K, g.0, g.1, g.0, g.1, two_g); // doubling, dx=0
    let inv = run(K, g.0, g.1, neg_g.0, neg_g.1, (U256::ZERO, U256::ZERO)); // P=-Q, dx=0
    let inf = run(K, U256::ZERO, U256::ZERO, q.0, q.1, q); // ∞ accumulator

    let show = |name: &str, a: &Agg| {
        eprintln!(
            "  {name:<22} phase_clean {:>2}/{:<2}  ancilla_clean {:>2}/{:<2}  out==ref {:>2}/{:<2}",
            a.phase_clean, a.seeds, a.ancilla_clean, a.seeds, a.out_matches_ref, a.seeds
        );
    };
    eprintln!("\n=== issue #5 completeness gating probe (K={K} seeds) ===");
    show("generic (dx!=0)", &gen);
    show("doubling (dx=0)", &dbl);
    show("P=-Q ->inf (dx=0)", &inv);
    show("inf accumulator", &inf);

    // Control: the generic case must be correct + clean on every seed (matches
    // the scored fuzzer). If not, the probe is bogus.
    assert_eq!(gen.out_matches_ref, gen.seeds, "generic control wrong output");
    assert_eq!(gen.phase_clean, gen.seeds, "generic control phase dirty");
    assert_eq!(gen.ancilla_clean, gen.seeds, "generic control ancilla dirty");

    // FINDINGS (locked as a regression). All three exceptional inputs share the
    // SAME signature: the ancilla always return to |0> (no register-basis leak),
    // but the OUTPUT is always wrong AND the global PHASE is corrupted on a
    // large fraction of RNG draws — i.e. probabilistic PHASE GARBAGE from
    // uncorrected HMR/R kickback (the kickmix "misused MBUC" failure mode). This
    // is Shor-sensitive: it is NOT merely "wrong output", so completeness must
    // bound the exceptional AMPLITUDE (Path A) and structurally remove the
    // ∞-accumulator (which starts at amplitude 1), not rely on output-only args.
    for (name, a) in [("doubling", &dbl), ("P=-Q", &inv), ("inf", &inf)] {
        assert_eq!(a.ancilla_clean, a.seeds, "{name}: ancilla leaked (expected clean)");
        assert_eq!(a.out_matches_ref, 0, "{name}: output unexpectedly correct");
        assert!(
            a.phase_clean < a.seeds,
            "{name}: phase was clean on all {} seeds — re-examine the phase-garbage claim",
            a.seeds
        );
    }
}

fn sub_mod_p(a: U256, b: U256, p: U256) -> U256 {
    let a = a % p;
    let b = b % p;
    if a >= b {
        a - b
    } else {
        p - (b - a)
    }
}
