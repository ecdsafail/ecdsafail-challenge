//! UNTRUSTED stage of the challenge harness.
//!
//! Invokes `point_add::build()` (the contestant's editable code) and
//! writes the resulting op stream to `ops.bin`. Nothing else runs here:
//! no simulation, no scoring, no `score.json`. The trusted `eval_circuit`
//! binary re-reads `ops.bin` from disk in a separate process so contestant
//! code cannot influence the score.

use quantum_ecc::circuit::Op;
use quantum_ecc::point_add;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

const OPS_PATH: &str = "ops.bin";
// "Z" suffix marks the zstd-compressed framing (was "QECCOPS1", uncompressed).
const MAGIC: &[u8; 8] = b"QECCOPSZ";
// zstd compression level. The record stream is almost pure boilerplate
// (zero pad bytes, NO_QUBIT sentinels, near-sequential ids), so even the
// fast default crushes it ~40x. Override with ZSTD_LEVEL for smaller/slower.
const DEFAULT_ZSTD_LEVEL: i32 = 3;

// Per-op layout (56 bytes, all little-endian):
//   u32  kind     | u32 _pad
//   u64  q_control2
//   u64  q_control1
//   u64  q_target
//   u64  c_target
//   u64  c_condition
//   u64  r_target
const OP_BYTES: usize = 56;

fn zstd_level() -> i32 {
    std::env::var("ZSTD_LEVEL")
        .ok()
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(DEFAULT_ZSTD_LEVEL)
}

fn write_ops(ops: &[Op], path: &Path) -> std::io::Result<()> {
    let tmp = path.with_extension("bin.tmp");
    let mut w = BufWriter::new(File::create(&tmp)?);
    // Header stays uncompressed so the trusted loader can read MAGIC and the
    // op count (to bound memory) before touching the compressed body.
    w.write_all(MAGIC)?;
    w.write_all(&(ops.len() as u64).to_le_bytes())?;

    // Compress the record body. The encoder owns the BufWriter; finish()
    // flushes the zstd frame and hands it back.
    let mut enc = zstd::stream::write::Encoder::new(w, zstd_level())?;
    // Serialize each op into a fixed 56-byte record and stream it through the
    // encoder, so we never materialize a second full copy of the op stream.
    let mut rec = [0u8; OP_BYTES];
    for op in ops {
        rec[0..4].copy_from_slice(&(op.kind as u32).to_le_bytes());
        rec[4..8].copy_from_slice(&[0u8; 4]); // pad
        rec[8..16].copy_from_slice(&op.q_control2.0.to_le_bytes());
        rec[16..24].copy_from_slice(&op.q_control1.0.to_le_bytes());
        rec[24..32].copy_from_slice(&op.q_target.0.to_le_bytes());
        rec[32..40].copy_from_slice(&op.c_target.0.to_le_bytes());
        rec[40..48].copy_from_slice(&op.c_condition.0.to_le_bytes());
        rec[48..56].copy_from_slice(&op.r_target.0.to_le_bytes());
        enc.write_all(&rec)?;
    }
    let mut w = enc.finish()?;
    w.flush()?;
    drop(w);
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
    let uncompressed = ops.len() * OP_BYTES + MAGIC.len() + 8;
    let on_disk = fs::metadata(path).map(|m| m.len()).unwrap_or(0);
    println!(
        "  wrote       : {} ({} bytes on disk, {} uncompressed, {:.1}x)",
        OPS_PATH,
        on_disk,
        uncompressed,
        uncompressed as f64 / on_disk.max(1) as f64,
    );
    println!("\n=== build_circuit OK ===");
}
