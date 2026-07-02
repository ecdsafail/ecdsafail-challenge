//! Analysis-only tool: reads `ops.bin` and reports circuit depth metrics
//! (toffoli/T-depth and total gate depth) via `circuit::analyze_depth`.
//!
//! Writes `depth.json` alongside `score.json`. It deliberately does NOT touch
//! `score.json`, run the simulator, or change the score — it only measures the
//! critical path of the already-built op stream so `analysis/cost_model.py` can
//! report a measured spacetime volume instead of a sequential upper bound.
//!
//! This is not the trusted scorer, so it uses a plain reader for our own
//! `ops.bin`; adversarial hardening lives in `eval_circuit`.
use quantum_ecc::circuit::{
    analyze_depth, analyze_ops, BitId, Op, OperationType, QubitId, RegisterId,
};
use std::io::{BufReader, Read};

const OPS_PATH: &str = "ops.bin";
const MAGIC: &[u8; 8] = b"QECCOPSZ";
const OP_BYTES: usize = 7 * 8;

fn read_u64(rec: &[u8], off: usize) -> u64 {
    u64::from_le_bytes(rec[off..off + 8].try_into().unwrap())
}

fn kind_from_u32(v: u32) -> Option<OperationType> {
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

fn load_ops(path: &str) -> Result<Vec<Op>, String> {
    let mut file = std::fs::File::open(path).map_err(|e| format!("open {path}: {e}"))?;
    let mut header = [0u8; MAGIC.len() + 8];
    file.read_exact(&mut header)
        .map_err(|e| format!("{path}: short header: {e}"))?;
    if &header[..MAGIC.len()] != MAGIC {
        return Err(format!("{path}: bad magic"));
    }
    let n = u64::from_le_bytes(header[MAGIC.len()..].try_into().unwrap()) as usize;
    let mut dec = zstd::stream::read::Decoder::new(BufReader::new(file))
        .map_err(|e| format!("{path}: zstd init: {e}"))?;
    let mut ops = Vec::with_capacity(n);
    let mut rec = [0u8; OP_BYTES];
    for i in 0..n {
        dec.read_exact(&mut rec)
            .map_err(|e| format!("op {i}: short read: {e}"))?;
        let kind = kind_from_u32(u32::from_le_bytes(rec[0..4].try_into().unwrap()))
            .ok_or_else(|| format!("op {i}: unknown kind"))?;
        ops.push(Op {
            kind,
            q_control2: QubitId(read_u64(&rec, 8)),
            q_control1: QubitId(read_u64(&rec, 16)),
            q_target: QubitId(read_u64(&rec, 24)),
            c_target: BitId(read_u64(&rec, 32)),
            c_condition: BitId(read_u64(&rec, 40)),
            r_target: RegisterId(read_u64(&rec, 48)),
        });
    }
    Ok(ops)
}

fn main() {
    let ops = match load_ops(OPS_PATH) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("!! could not load {OPS_PATH}: {e}");
            std::process::exit(1);
        }
    };
    let (num_qubits, num_bits, _regs, _r) = analyze_ops(ops.iter());
    let d = analyze_depth(ops.iter(), num_qubits as usize, num_bits as usize);

    println!("=== depth report ({} ops) ===", ops.len());
    println!("  qubits             : {num_qubits}");
    println!("  bits               : {num_bits}");
    println!("  toffoli_depth (T*)  : {}", d.toffoli_depth);
    println!("  gate_depth          : {}", d.gate_depth);
    if d.toffoli_depth > 0 {
        println!(
            "  gate parallelism    : {:.2}x  (emitted ops / gate_depth)",
            ops.len() as f64 / d.gate_depth.max(1) as f64
        );
    }

    let body = format!(
        "{{\n  \"toffoli_depth\": {},\n  \"gate_depth\": {},\n  \"qubits\": {}\n}}\n",
        d.toffoli_depth, d.gate_depth, num_qubits
    );
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/depth.json");
    if let Err(e) = std::fs::write(path, body) {
        eprintln!("warning: failed to write depth.json: {e}");
    } else {
        println!("  wrote depth.json");
    }
}
