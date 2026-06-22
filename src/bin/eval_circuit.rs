//! TRUSTED stage of the challenge harness.
//!
//! Reads the op stream produced by `build_circuit` from `ops.bin`,
//! re-simulates the circuit against the secp256k1 reference adder,
//! enforces the four validity checks (correctness, reversibility, phase,
//! ancilla cleanup), counts gates, writes `score.json`, and appends one
//! row to `results.tsv`.
//!
//! This binary deliberately does NOT import `quantum_ecc::point_add` —
//! contestant code never executes inside the trusted process. `ops.bin`
//! is treated as fully untrusted input and is bounds-checked before use.

use alloy_primitives::U256;
use quantum_ecc::circuit::{
    analyze_ops, BitId, Op, OperationType, QubitId, QubitOrBit, RegisterId,
};
use quantum_ecc::sim::Simulator;
use quantum_ecc::weierstrass_elliptic_curve::WeierstrassEllipticCurve;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};
use std::fs::{File, OpenOptions};
use std::io::{BufReader, Read, Write};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

const OPS_PATH: &str = "ops.bin";
// "Z" framing: 16-byte plaintext header (MAGIC + u64 count) then a zstd
// frame of the fixed-width records. The count is read before decompressing
// so we can bound memory, and we read exactly count*OP_BYTES bytes out of
// the decoder so a crafted frame cannot expand without bound.
const MAGIC: &[u8; 8] = b"QECCOPSZ";
// Cap the zstd window the decoder will accept (2^27 = 128 MiB). A forged
// ops.bin cannot force a huge decompression-window allocation.
const ZSTD_WINDOW_LOG_MAX: u32 = 27;
const FIELD_BYTES: usize = 8;
const OP_FIELDS: usize = 7;
const OP_BYTES: usize = OP_FIELDS * FIELD_BYTES;
// Resource caps. These are sanity limits to prevent a malicious ops.bin
// from OOM'ing the simulator (each qubit/bit costs 8 bytes in `Simulator::new`).
// Real circuits sit ~3k qubits and a few hundred million ops; caps are
// generous compared to that.
const MAX_OPS: u64 = 4_000_000_000;
const NUM_TESTS: usize = 9024;

// ─── Bounded ops.bin loader ────────────────────────────────────────────────
//
// Hand-rolled fixed-width LE framing. Per-op layout:
//   u32 kind (must be 0..=17 — rejects the rkyv jump-table-confusion attack
//             that bypassed the Toffoli counter in the Google challenge)
//   u32 _pad
//   u64 q_control2  (NO_QUBIT = u64::MAX = unused)
//   u64 q_control1
//   u64 q_target
//   u64 c_target    (NO_BIT = u64::MAX)
//   u64 c_condition
//   u64 r_target    (NO_REG = u64::MAX)
//
// After reassembly, each Op is fed to Op::validate() (upstream zkp_ecc
// post-incident hardening). validate() panics on:
//   - operand aliasing (CCX q q q etc. — would yield free non-reversible
//     resets, the ToB "strictly better exploit primitive")
//   - per-kind field-shape violations (e.g. R/Hmr with c_condition,
//     which would suppress phase randomization on dirty frees)
// We catch_unwind so a forged ops.bin produces an error, not a crash.

fn op_kind_from_u32(v: u32) -> Option<OperationType> {
    Some(match v {
        0 => OperationType::Neg,
        1 => OperationType::Register,
        2 => OperationType::AppendToRegister,
        3 => OperationType::BitInvert,
        4 => OperationType::BitStore0,
        5 => OperationType::BitStore1,
        6 => OperationType::X,
        7 => OperationType::Z,
        8 => OperationType::CX,
        9 => OperationType::CZ,
        10 => OperationType::Swap,
        11 => OperationType::R,
        12 => OperationType::Hmr,
        13 => OperationType::CCX,
        14 => OperationType::CCZ,
        15 => OperationType::PushCondition,
        16 => OperationType::PopCondition,
        17 => OperationType::DebugPrint,
        _ => return None,
    })
}

fn read_u64(bytes: &[u8], off: usize) -> u64 {
    u64::from_le_bytes(bytes[off..off + 8].try_into().unwrap())
}

fn load_ops(path: &str) -> Result<Vec<Op>, String> {
    let mut file = File::open(path).map_err(|e| format!("open {path}: {e}"))?;

    // Plaintext header: MAGIC + u64 op count. Read and validate before
    // decompressing so the op count bounds every allocation below.
    let mut header = [0u8; MAGIC.len() + 8];
    file.read_exact(&mut header)
        .map_err(|e| format!("{path}: too short to read header: {e}"))?;
    if &header[..MAGIC.len()] != MAGIC {
        return Err(format!("{path}: bad magic"));
    }
    let n = u64::from_le_bytes(header[MAGIC.len()..].try_into().unwrap());
    if n > MAX_OPS {
        return Err(format!("{path}: op count {n} exceeds cap {MAX_OPS}"));
    }
    let n = n as usize;

    // Stream-decompress the record body. We read exactly n * OP_BYTES bytes
    // out of the decoder, so a forged frame cannot expand without bound; the
    // window cap limits the decoder's own buffer allocation.
    let mut dec = zstd::stream::read::Decoder::new(BufReader::new(file))
        .map_err(|e| format!("{path}: zstd init: {e}"))?;
    dec.window_log_max(ZSTD_WINDOW_LOG_MAX)
        .map_err(|e| format!("{path}: zstd window cap: {e}"))?;

    let mut ops = Vec::with_capacity(n);
    let mut rec = [0u8; OP_BYTES];
    for i in 0..n {
        dec.read_exact(&mut rec)
            .map_err(|e| format!("op {i}: short read from compressed body: {e}"))?;
        let kind_raw = u32::from_le_bytes(rec[0..4].try_into().unwrap());
        let kind = op_kind_from_u32(kind_raw)
            .ok_or_else(|| format!("op {i}: unknown kind {kind_raw}"))?;
        // rec[4..8] are reserved padding for 8-byte alignment.
        let q_control2 = QubitId(read_u64(&rec, 8));
        let q_control1 = QubitId(read_u64(&rec, 16));
        let q_target = QubitId(read_u64(&rec, 24));
        let c_target = BitId(read_u64(&rec, 32));
        let c_condition = BitId(read_u64(&rec, 40));
        let r_target = RegisterId(read_u64(&rec, 48));

        let op = Op {
            kind,
            q_control2,
            q_control1,
            q_target,
            c_target,
            c_condition,
            r_target,
        };
        // Op::validate() panics on aliasing or per-kind field-shape errors.
        // Catch the unwind to convert into a clean rejection.
        let validated = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| op.validate()));
        if let Err(e) = validated {
            let msg = e
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| e.downcast_ref::<&'static str>().map(|s| s.to_string()))
                .unwrap_or_else(|| "validation panic".to_string());
            return Err(format!("op {i}: {msg}"));
        }
        ops.push(op);
    }

    // Reject trailing data: exactly n records must decompress, no more.
    let mut extra = [0u8; 1];
    match dec.read(&mut extra) {
        Ok(0) => {}
        Ok(_) => return Err(format!("{path}: trailing data after {n} ops")),
        Err(e) => return Err(format!("{path}: error checking for trailing data: {e}")),
    }
    Ok(ops)
}

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

// ─── Fiat-Shamir seed ──────────────────────────────────────────────────────
//
// SHAKE256 over the op stream. Determines test inputs, simulator RNG for
// R/Hmr phase randomization, etc.

fn fiat_shamir_seed(ops: &[Op]) -> sha3::Shake256Reader {
    let mut hasher = Shake256::default();
    hasher.update(b"quantum_ecc-fiat-shamir-v2");
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

// ─── Test runner ──────────────────────────────────────────────────────────

struct SeedReport {
    ok: bool,
    avg_cliff: f64,
    avg_tof: f64,
    tot_tof: u64,
    tot_cliff: u64,
    n_shots: usize,
    classical_failures: usize,
    phase_garbage_batches: usize,
    ancilla_garbage_batches: usize,
    fail_reason: Option<String>,
}

fn run_tests(
    ops: &[Op],
    layout_regs: &[Vec<QubitOrBit>],
    total_qubits: u64,
    num_bits: u64,
    mut xof: sha3::Shake256Reader,
    target_shots: usize,
) -> SeedReport {
    let curve = secp256k1();

    let mut targets = Vec::with_capacity(target_shots);
    let mut offsets = Vec::with_capacity(target_shots);
    let mut expected = Vec::with_capacity(target_shots);
    for _ in 0..target_shots {
        let mut rb = [[0u8; 32]; 2];
        // Disambiguate from std::io::Read (in scope for the zstd loader).
        XofReader::read(&mut xof, &mut rb[0]);
        XofReader::read(&mut xof, &mut rb[1]);
        let k1 = U256::from_le_bytes(rb[0]);
        let k2 = U256::from_le_bytes(rb[1]);
        let t = curve.mul(curve.gx, curve.gy, k1);
        let o = curve.mul(curve.gx, curve.gy, k2);
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
        let cond_mask: u64 = if bs == 64 { u64::MAX } else { (1u64 << bs) - 1 };

        sim.clear_for_shot();
        for shot in 0..bs {
            let i = batch * BATCH + shot;
            sim.set_register(&layout_regs[0], targets[i].0, shot);
            sim.set_register(&layout_regs[1], targets[i].1, shot);
            sim.set_register(&layout_regs[2], offsets[i].0, shot);
            sim.set_register(&layout_regs[3], offsets[i].1, shot);
        }

        sim.apply_iter(ops.iter());

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

        let phase = sim.phase & cond_mask;
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

        for register in layout_regs {
            for qb in register {
                if let QubitOrBit::Qubit(q) = *qb {
                    *sim.qubit_mut(q) = 0;
                }
            }
        }
        let mut garbage_q: Option<u64> = None;
        for q in 0..total_qubits {
            let v = sim.qubit(QubitId(q)) & cond_mask;
            if v != 0 {
                garbage_q = Some(q);
                break;
            }
        }
        if let Some(q) = garbage_q {
            ancilla_garbage_batches += 1;
            let v = sim.qubit(QubitId(q)) & cond_mask;
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

    let _ = num_bits;
    let denom = n.max(1) as f64;
    SeedReport {
        ok,
        avg_cliff: sim.stats.clifford_gates as f64 / denom,
        avg_tof: sim.stats.toffoli_gates as f64 / denom,
        tot_tof: sim.stats.toffoli_gates,
        tot_cliff: sim.stats.clifford_gates,
        n_shots: n,
        classical_failures,
        phase_garbage_batches,
        ancilla_garbage_batches,
        fail_reason,
    }
}

// ─── Output bookkeeping ────────────────────────────────────────────────────

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
    qubits: u64,
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

fn write_score(avg_tof: f64, qubits: u64) {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/score.json");
    let toffoli = avg_tof.round() as u64;
    let score = toffoli.saturating_mul(qubits);
    let body = format!(
        "{{\n  \"score\": {score},\n  \"metrics\": {{\n    \"toffoli\": {toffoli},\n    \"qubits\": {qubits}\n  }}\n}}\n"
    );
    if let Err(e) = std::fs::write(path, body) {
        eprintln!("warning: failed to write score.json: {e}");
    }
}

fn fail_and_exit(reason: &str, note: &str, ops_len: usize, total_qubits: u64) -> ! {
    eprintln!("\n!! eval FAILED: {reason}");
    let fail_note = if note.is_empty() {
        reason.to_string()
    } else {
        format!("{note} | {reason}")
    };
    append_results_row("FAIL", 0.0, 0.0, total_qubits, ops_len, &fail_note);
    std::process::exit(1);
}

fn main() {
    let note = parse_note();
    println!("=== quantum_ecc: eval_circuit (trusted stage) ===\n");

    let ops = match load_ops(OPS_PATH) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("!! could not load {OPS_PATH}: {e}");
            append_results_row("FAIL", 0.0, 0.0, 0, 0, &format!("{note} | load: {e}"));
            std::process::exit(1);
        }
    };
    println!("  loaded ops  : {}", ops.len());

    let (total_qubits, num_bits, _num_regs, regs) = analyze_ops(ops.iter());

    if regs.len() != 4 {
        fail_and_exit(
            &format!("expected 4 registers, got {}", regs.len()),
            &note,
            ops.len(),
            total_qubits,
        );
    }
    for (i, r) in regs.iter().enumerate() {
        if r.len() != 256 {
            fail_and_exit(
                &format!("register {i} should be 256 wide, got {}", r.len()),
                &note,
                ops.len(),
                total_qubits,
            );
        }
    }
    for q in &regs[0] {
        if !matches!(q, QubitOrBit::Qubit(_)) {
            fail_and_exit("register 0 must be qubits", &note, ops.len(), total_qubits);
        }
    }
    for q in &regs[1] {
        if !matches!(q, QubitOrBit::Qubit(_)) {
            fail_and_exit("register 1 must be qubits", &note, ops.len(), total_qubits);
        }
    }
    for q in &regs[2] {
        if !matches!(q, QubitOrBit::Bit(_)) {
            fail_and_exit("register 2 must be bits", &note, ops.len(), total_qubits);
        }
    }
    for q in &regs[3] {
        if !matches!(q, QubitOrBit::Bit(_)) {
            fail_and_exit("register 3 must be bits", &note, ops.len(), total_qubits);
        }
    }

    println!("  qubits      : {}", total_qubits);
    println!("  bits        : {}", num_bits);

    println!("\n-- correctness tests ({} shots) --", NUM_TESTS);
    let xof = fiat_shamir_seed(&ops);
    let r = run_tests(&ops, &regs, total_qubits, num_bits, xof, NUM_TESTS);
    println!("  tested shots            : {}", r.n_shots);
    println!("  classical mismatches    : {}", r.classical_failures);
    println!("  phase-garbage batches   : {}", r.phase_garbage_batches);
    println!("  ancilla-garbage batches : {}", r.ancilla_garbage_batches);
    if !r.ok {
        let reason = r.fail_reason.clone().unwrap_or_else(|| "(no detail)".into());
        println!("\n!! correctness FAILED: {reason}");
        let fail_note = format!("{note} | {reason}");
        append_results_row(
            "FAIL",
            r.avg_tof,
            r.avg_cliff,
            total_qubits,
            ops.len(),
            &fail_note,
        );
        std::process::exit(1);
    }
    println!("  all {} shots OK", r.n_shots);

    println!("\n=== circuit metrics (secp256k1, n=256) ===");
    println!("  avg executed Toffoli  : {:.3}", r.avg_tof);
    println!("  avg executed Clifford : {:.3}", r.avg_cliff);
    println!(
        "  total Toffoli (sum)   : {} over {} shots",
        r.tot_tof, r.n_shots
    );
    println!("  total Clifford (sum)  : {}", r.tot_cliff);
    println!("  emitted ops           : {}", ops.len());
    println!("  qubits                : {}", total_qubits);

    append_results_row("OK", r.avg_tof, r.avg_cliff, total_qubits, ops.len(), &note);
    write_score(r.avg_tof, total_qubits);

    println!("\n=== experiment OK ===");
}
