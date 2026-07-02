//! Library entry point shared by the two challenge binaries:
//!
//! * `build_circuit` (untrusted) calls `point_add::build` and serializes the
//!   resulting op stream to `ops.bin`. This is the only place contestant
//!   code executes.
//! * `eval_circuit` (trusted) reads `ops.bin`, re-simulates against the
//!   secp256k1 reference adder, validates correctness/reversibility/phase,
//!   counts gates, and writes `score.json`.
//!
//! Splitting the harness into two binaries keeps the trust boundary aligned
//! with the process boundary: contestant code cannot tamper with the
//! simulator, the test inputs, or the score because none of that runs in
//! `build_circuit`.

#[allow(dead_code)]
pub mod circuit;
pub mod point_add;
#[allow(dead_code)]
pub mod sim;
#[allow(dead_code)]
pub mod weierstrass_elliptic_curve;

// Bit-precise Kani proofs binding the arithmetic contract to the real Rust
// types. Compiled only under `cargo kani` (--cfg kani); invisible to the normal
// build and to benchmark.sh, so it cannot affect the circuit or the score.
#[cfg(kani)]
mod kani_proofs;
