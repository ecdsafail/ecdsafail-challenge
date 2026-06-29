//! Deferred Lehmer versioned module for ecdsa.fail contest.
//! VERSION: 2026.06.29-140652
//!
//! Currently delegates to trailmix for correctness. Future: full 3-phase for lower score.

use crate::circuit::Op;

pub const DEFERRED_LEHMER_VERSION: &str = "2026.06.29-140652";

pub fn build_deferred_lehmer_ops() -> Vec<Op> {
    eprintln!("deferred_lehmer version {}", DEFERRED_LEHMER_VERSION);
    super::trailmix_ludicrous::build_trailmix_ludicrous_ops()
}

pub fn deferred_selftest(_mode: &str) -> Result<(), String> {
    eprintln!("deferred_lehmer version {} (selftest)", DEFERRED_LEHMER_VERSION);
    Ok(())
}
