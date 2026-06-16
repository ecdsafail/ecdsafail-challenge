//! Schrottenloher 2026 (arxiv 2606.02235) point-addition replication.
//!
//! This module implements the published recipe from
//! "Optimized Point Addition Circuits for Elliptic Curve Discrete
//! Logarithms" by A. Schrottenloher (2026). The paper reports reproducing
//! Babbush et al.'s secp256k1 point-addition cost (~1175 qubits /
//! 2^21.36 Toffoli) with a fully published logical-circuit recipe, at
//! slightly higher qubits and slightly fewer Toffoli than Babbush.
//!
//! Design choices here (matching the paper):
//!   - Approximate posture: top-MSB comparators, ~2^-13 tail failure
//!     rate on random inputs. This matches the paper's Section 4 and
//!     our existing tolerance (Shor's algorithm requires only a
//!     constant success probability per call).
//!   - Pseudo-Mersenne specialization (Alg 7, 10, 11) for
//!     secp256k1, q = 2^256 - f, f = 2^32 + 977. Generic-prime
//!     variants (Alg 6, 9) also provided.
//!   - In-place modular multiplication via the Khattar-Shutty-Gidney
//!     2025 "dialog" Extended Euclidean algorithm (Alg 2 + Alg 3),
//!     with the Fig.1 5-Toffoli 3-iter -> 5-bit garbage compression.
//!   - Top-level point-addition driver follows Haner-Roetteler 2020
//!     Alg 1 (same as Babbush, same as our existing
//!     `ec::deferred::ec_add_inplace_*`).
//!
//! The module lives alongside existing primitives; we do not modify
//! pre-existing `arith/mod_arith.rs`, `arith/rfold_mbu.rs`, or
//! `inversion/*` paths. Once validated, the
//! `pointadd::ec_add_inplace_schrottenloher` driver can be
//! benchmarked against `ec::deferred::ec_add_inplace_2eea`.

pub mod bezout_unpack;
pub mod gcd_compress;
pub mod gcd_compress5;
pub mod gcd_compress_jump;
pub mod gcd_jump;
pub mod gcd_pack;
pub mod jump_schedule;
pub mod mod_mul_eea;
pub mod msb_compare;
pub mod mul_pow2;
pub mod pm_prims;
pub mod pointadd;
pub mod wallace_square;

/// secp256k1 pseudo-Mersenne offset: q = 2^256 - `F_SECP256K1`.
/// Equal to our existing rfold R constant.
pub const F_SECP256K1: u64 = (1u64 << 32) + 977;

/// Width of the f-addition in pseudo-Mersenne add/double: enough to
/// hold f exactly. For secp256k1, f < 2^33.
pub const F_LSBS_SECP256K1: usize = 33;

/// Number of top MSBs to use in approximate comparators. The paper
/// uses 40-50 to drive a ~2^-13 tail failure rate; we pick 40 as a
/// conservative default.
pub const APPROX_CMP_MSBS_DEFAULT: usize = 40;
