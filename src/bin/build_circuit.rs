//! UNTRUSTED stage of the challenge harness.
//!
//! Invokes `point_add::build()` (the contestant's editable code) and
//! writes the resulting op stream to `ops.bin`. Nothing else runs here:
//! no simulation, no scoring, no `score.json`. The trusted `eval_circuit`
//! binary re-reads `ops.bin` from disk in a separate process so contestant
//! code cannot influence the score.

use quantum_ecc::circuit::Op;
use quantum_ecc::point_add;
use std::fs;
use std::path::Path;

const OPS_PATH: &str = "ops.bin";
const MAGIC: &[u8; 8] = b"QECCOPS1";

// Per-op layout (56 bytes, all little-endian):
//   u32  kind     | u32 _pad
//   u64  q_control2
//   u64  q_control1
//   u64  q_target
//   u64  c_target
//   u64  c_condition
//   u64  r_target
const OP_BYTES: usize = 56;

fn write_ops(ops: &[Op], path: &Path) -> std::io::Result<()> {
    let mut buf = Vec::with_capacity(MAGIC.len() + 8 + ops.len() * OP_BYTES);
    buf.extend_from_slice(MAGIC);
    buf.extend_from_slice(&(ops.len() as u64).to_le_bytes());
    for op in ops {
        buf.extend_from_slice(&(op.kind as u32).to_le_bytes());
        buf.extend_from_slice(&[0u8; 4]); // pad
        buf.extend_from_slice(&op.q_control2.0.to_le_bytes());
        buf.extend_from_slice(&op.q_control1.0.to_le_bytes());
        buf.extend_from_slice(&op.q_target.0.to_le_bytes());
        buf.extend_from_slice(&op.c_target.0.to_le_bytes());
        buf.extend_from_slice(&op.c_condition.0.to_le_bytes());
        buf.extend_from_slice(&op.r_target.0.to_le_bytes());
    }
    let tmp = path.with_extension("bin.tmp");
    fs::write(&tmp, &buf)?;
    fs::rename(&tmp, path)?;
    Ok(())
}

fn main() {
    println!("=== quantum_ecc: build_circuit (untrusted stage) ===\n");
    println!("-- building circuit --");
    let ops = point_add::build();
    println!("  emitted ops : {}", ops.len());

    let path = Path::new(OPS_PATH);
    if let Err(e) = write_ops(&ops, path) {
        eprintln!("error: failed to write {}: {}", OPS_PATH, e);
        std::process::exit(2);
    }
    println!(
        "  wrote       : {} ({} bytes)",
        OPS_PATH,
        ops.len() * OP_BYTES + MAGIC.len() + 8
    );
    println!("\n=== build_circuit OK ===");
}
