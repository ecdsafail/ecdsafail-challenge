//! Experiment harness: builds the point-addition circuit defined in
//! the `point_add/` module, runs it against the zenodo `Simulator` with
//! random test shots, and reports Toffoli / Clifford / qubit counts.
//!
//! Research-loop contract: ONLY files under `src/point_add/` are edited
//! by the loop. This file, `circuit.rs`, `sim.rs`, and
//! `weierstrass_elliptic_curve.rs` are harness and must not be touched.
//!
//! Attribution: `circuit.rs`, `sim.rs`, and `weierstrass_elliptic_curve.rs`
//! are reused verbatim from the `zkp_ecc` Zenodo project under CC BY 4.0.
//! See `NOTICE` at the repository root for details.

#[allow(dead_code)]
mod circuit;
mod point_add;
#[allow(dead_code)]
mod sim;
#[allow(dead_code)]
mod weierstrass_elliptic_curve;

use alloy_primitives::U256;
use circuit::{analyze_ops, Op, QubitOrBit};
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};
use sim::Simulator;
use std::fs::OpenOptions;
use std::io::Write;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};
use weierstrass_elliptic_curve::WeierstrassEllipticCurve;

// ─── secp256k1 parameters ──────────────────────────────────────────────────

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

// ─── Test runner ───────────────────────────────────────────────────────────

const NUM_TESTS: usize = 9024;

/// Hash the circuit's op stream into the seed XOF (Fiat-Shamir).
///
/// Mirrors upstream zenodo's protocol: test inputs are derived from
/// `Shake256(circuit_bytes)`, which makes them deterministic per circuit
/// and impossible to "tune" the circuit against. We don't have rkyv-
/// archived bytes handy here, so we feed each `Op`'s fields into the
/// hasher directly. Any two distinct op streams produce distinct seeds
/// (we hash the count + every field of every op).
fn fiat_shamir_seed(ops: &[Op]) -> sha3::Shake256Reader {
    let mut hasher = Shake256::default();
    hasher.update(b"quantum_ecc-fiat-shamir-v1");
    hasher.update(&(ops.len() as u64).to_le_bytes());
    for op in ops {
        hasher.update(&[op.kind as u8]);
        hasher.update(&op.q_control2.0.to_le_bytes());
        hasher.update(&op.q_control1.0.to_le_bytes());
        hasher.update(&op.q_target.0.to_le_bytes());
        hasher.update(&op.c_target.0.to_le_bytes());
        hasher.update(&op.c_condition.0.to_le_bytes());
        hasher.update(&op.r_target.0.to_le_bytes());
    }
    hasher.finalize_xof()
}

fn run_tests(
    ops: &[Op],
    layout_regs: &[Vec<QubitOrBit>],
    total_qubits: u32,
    num_bits: u32,
) -> (bool, f64, f64, u64, u64, usize, Option<String>) {
    let curve = secp256k1();
    // Fiat-Shamir: a single XOF seeded from the circuit op stream feeds
    // both test-input generation AND the simulator's RNG, exactly as
    // upstream zenodo's program does.
    let mut xof = fiat_shamir_seed(ops);

    // Generate random target/offset points as k*G.
    let mut targets = Vec::with_capacity(NUM_TESTS);
    let mut offsets = Vec::with_capacity(NUM_TESTS);
    let mut expected = Vec::with_capacity(NUM_TESTS);
    for _ in 0..NUM_TESTS {
        let mut rb = [[0u8; 32]; 2];
        xof.read(&mut rb[0]);
        xof.read(&mut rb[1]);
        let k1 = U256::from_le_bytes(rb[0]);
        let k2 = U256::from_le_bytes(rb[1]);
        let t = curve.mul(curve.gx, curve.gy, k1);
        let o = curve.mul(curve.gx, curve.gy, k2);
        // Avoid the doubling / inverse-pair cases the baseline doesn't handle.
        if t.0 == o.0 {
            continue;
        }
        if t.0.is_zero() && t.1.is_zero() {
            continue;
        }
        if o.0.is_zero() && o.1.is_zero() {
            continue;
        }
        let e = curve.add(t.0, t.1, o.0, o.1);
        targets.push(t);
        offsets.push(o);
        expected.push(e);
    }
    let n = targets.len();

    let mut sim = Simulator::new(total_qubits as usize, num_bits as usize, &mut xof);
    let mut ok = true;
    let mut fail_reason: Option<String> = None;
    let mut classical_failures = 0usize;
    let mut phase_garbage_batches = 0usize;
    let mut ancilla_garbage_batches = 0usize;

    const BATCH: usize = 64;
    let num_batches = (n + BATCH - 1) / BATCH;
    for batch in 0..num_batches {
        let bs = BATCH.min(n - batch * BATCH);
        // `cond_mask`: bit i is 1 iff shot i is "live" in this batch.
        // Used for the phase + end-state garbage checks below.
        let cond_mask: u64 = if bs == 64 { u64::MAX } else { (1u64 << bs) - 1 };

        sim.clear_for_shot();
        for shot in 0..bs {
            let i = batch * BATCH + shot;
            sim.set_register(&layout_regs[0], targets[i].0, shot);
            sim.set_register(&layout_regs[1], targets[i].1, shot);
            sim.set_register(&layout_regs[2], offsets[i].0, shot);
            sim.set_register(&layout_regs[3], offsets[i].1, shot);
        }

        // ─── Forward pass ────────────────────────────────────────────────
        // sim.rs is upstream verbatim. R randomizes the phase on dirty
        // frees; the phase check below catches it probabilistically.
        sim.apply(ops);

        // ─── Correctness check ──────────────────────────────────────────
        for shot in 0..bs {
            let i = batch * BATCH + shot;
            let gx = sim.get_register(&layout_regs[0], shot);
            let gy = sim.get_register(&layout_regs[1], shot);
            if gx != expected[i].0 || gy != expected[i].1 {
                classical_failures += 1;
                if fail_reason.is_none() {
                    fail_reason = Some(format!(
                        "CLASSICAL MISMATCH shot {i}: got ({:#x},{:#x}) exp ({:#x},{:#x})",
                        gx, gy, expected[i].0, expected[i].1
                    ));
                }
                ok = false;
            }
        }

        // ─── Phase garbage check ────────────────────────────────────────
        // Upstream zenodo's protocol: after forward, the global phase must
        // be 0 across all live shots. Catches misuses of phase-flipping
        // gates (Z/CZ/CCZ) or bad Hmr uncomputation.
        let phase = sim.global_phase() & cond_mask;
        if phase != 0 {
            phase_garbage_batches += 1;
            let msg = format!(
                "PHASE GARBAGE: global_phase = {:#018x} across {} live shots (must be 0)",
                phase, bs
            );
            if fail_reason.is_none() {
                fail_reason = Some(msg);
            }
            ok = false;
        }

        // ─── End-state ancillary garbage check ──────────────────────────
        // Upstream zenodo's reversibility contract: at end of forward, every
        // qubit OUTSIDE the declared output registers must be |0⟩ on every
        // live shot. We zero the register qubits first, then sweep.
        for register in layout_regs {
            for qb in register {
                if let QubitOrBit::Qubit(q) = *qb {
                    *sim.qubit_mut(q) = 0;
                }
            }
        }
        let mut garbage_q: Option<u32> = None;
        for q in 0..total_qubits {
            let v = sim.qubit(circuit::QubitId(q)) & cond_mask;
            if v != 0 {
                garbage_q = Some(q);
                break;
            }
        }
        if let Some(q) = garbage_q {
            ancilla_garbage_batches += 1;
            let v = sim.qubit(circuit::QubitId(q)) & cond_mask;
            let msg = format!(
                "ANCILLA GARBAGE: qubit {} = {:#018x} (live shots) at end of forward; \
                 every non-register qubit must be |0⟩ on every live shot",
                q, v
            );
            if fail_reason.is_none() {
                fail_reason = Some(msg);
            }
            ok = false;
        }
    }

    println!("  tested shots            : {}", n);
    println!("  classical mismatches    : {}", classical_failures);
    println!("  phase-garbage batches   : {}", phase_garbage_batches);
    println!("  ancilla-garbage batches : {}", ancilla_garbage_batches);

    let denom = n.max(1) as f64;
    let avg_cliff = sim.stats.clifford_gates as f64 / denom;
    let avg_tof = sim.stats.toffoli_gates as f64 / denom;
    (
        ok,
        avg_cliff,
        avg_tof,
        sim.stats.toffoli_gates,
        sim.stats.clifford_gates,
        n,
        fail_reason,
    )
}

fn parse_note() -> String {
    let mut args = std::env::args().skip(1);
    let mut note = String::new();
    while let Some(a) = args.next() {
        if a == "--note" {
            if let Some(v) = args.next() {
                note = v;
            }
        } else if let Some(rest) = a.strip_prefix("--note=") {
            note = rest.to_string();
        }
    }
    note.replace('\t', " ").replace('\n', " ")
}

fn git_commit_short() -> String {
    Command::new("git")
        .args(["rev-parse", "--short", "HEAD"])
        .output()
        .ok()
        .and_then(|o| {
            if o.status.success() {
                Some(String::from_utf8_lossy(&o.stdout).trim().to_string())
            } else {
                None
            }
        })
        .unwrap_or_else(|| "nogit".to_string())
}

fn append_results_row(
    correct: &str,
    avg_tof: f64,
    avg_cliff: f64,
    qubits: u32,
    ops_len: usize,
    note: &str,
) {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let commit = git_commit_short();
    let safe_note = note.replace('\t', " ").replace('\n', " ");
    let row = format!(
        "{ts}\t{commit}\t{avg_tof:.3}\t{avg_cliff:.3}\t{qubits}\t{ops_len}\t{correct}\t{safe_note}\n"
    );
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/results.tsv");
    match OpenOptions::new().create(true).append(true).open(path) {
        Ok(mut f) => {
            if let Err(e) = f.write_all(row.as_bytes()) {
                eprintln!("warning: failed to write results.tsv: {e}");
            }
        }
        Err(e) => eprintln!("warning: failed to open results.tsv: {e}"),
    }
}

fn write_score(avg_tof: f64, qubits: u32) {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/score.txt");
    let score = avg_tof * f64::from(qubits);
    if let Err(e) = std::fs::write(path, format!("{score:.3}\n")) {
        eprintln!("warning: failed to write score.txt: {e}");
    }
}

fn main() {
    let note = parse_note();
    println!("=== quantum_ecc: secp256k1 point addition baseline ===\n");

    println!("-- building circuit --");
    let ops = point_add::build();

    let (total_qubits, num_bits, _num_regs, regs) = analyze_ops(ops.iter().copied());

    // Sanity-check layout matches zenodo's program interface.
    // The 4 registers and their widths/types are the only contract `build`
    // must satisfy with the harness.
    assert!(
        regs.len() == 4,
        "expected 4 registers (target_x, target_y, offset_x, offset_y); got {}",
        regs.len()
    );
    for (i, r) in regs.iter().enumerate() {
        assert_eq!(
            r.len(),
            256,
            "register {i} should be 256 wide, got {}",
            r.len()
        );
    }
    for q in &regs[0] {
        assert!(
            matches!(q, QubitOrBit::Qubit(_)),
            "register 0 must be qubits"
        );
    }
    for q in &regs[1] {
        assert!(
            matches!(q, QubitOrBit::Qubit(_)),
            "register 1 must be qubits"
        );
    }
    for q in &regs[2] {
        assert!(matches!(q, QubitOrBit::Bit(_)), "register 2 must be bits");
    }
    for q in &regs[3] {
        assert!(matches!(q, QubitOrBit::Bit(_)), "register 3 must be bits");
    }

    println!("  total ops : {}", ops.len());
    println!("  qubits    : {}", total_qubits);
    println!("  bits      : {}", num_bits);

    println!("\n-- running correctness tests --");
    let (ok, avg_cliff, avg_tof, tot_tof, tot_cliff, n_shots, fail_reason) =
        run_tests(&ops, &regs, total_qubits, num_bits);
    if !ok {
        println!("\n!! correctness FAILED");
        if let Some(r) = &fail_reason {
            println!("  first failure          : {}", r);
        }
        let fail_note = match &fail_reason {
            Some(r) => format!("{note} | {r}"),
            None => note.clone(),
        };
        append_results_row(
            "FAIL",
            avg_tof,
            avg_cliff,
            total_qubits,
            ops.len(),
            &fail_note,
        );
        std::process::exit(1);
    }
    println!("  all {} shots OK", n_shots);

    println!("\n=== circuit metrics (secp256k1, n=256) ===");
    println!("  avg executed Toffoli  : {:.3}", avg_tof);
    println!("  avg executed Clifford : {:.3}", avg_cliff);
    println!(
        "  total Toffoli (sum)   : {} over {} shots",
        tot_tof, n_shots
    );
    println!("  total Clifford (sum)  : {}", tot_cliff);
    println!("  emitted ops           : {}", ops.len());
    println!("  qubits                : {}", total_qubits);

    append_results_row("OK", avg_tof, avg_cliff, total_qubits, ops.len(), &note);
    write_score(avg_tof, total_qubits);

    println!("\n=== experiment OK ===");
}
