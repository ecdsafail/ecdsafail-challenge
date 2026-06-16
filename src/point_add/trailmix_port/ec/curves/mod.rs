//! Port of the Schrottenloher EC point-addition design (see
//! `arith/schrottenloher/`) to other 256-bit curves: Curve25519, SM2,
//! and Brainpool P256.
//!
//! The affine point-addition formula is curve-`a`-independent, so the
//! Weierstrass driver and the classical reference are shared across
//! curves; only the modular-reduction primitives differ:
//!   - Curve25519: pseudo-Mersenne (q = 2^255 - 19), reuses the
//!     generic-over-f `pm_prims` core with f = 19.
//!   - SM2 / Brainpool: generic-prime reduction (Alg 6/9 = shift/add
//!     then MSB-compare + conditional subtract of q), no special prime
//!     structure required.
//!
//! Curve25519 is modeled in its birationally-equivalent short-Weierstrass
//! form over F_(2^255-19): the group/DLOG attack cost is what we measure,
//! so a Weierstrass affine add over the same field is the faithful port of
//! the design (native twisted-Edwards addition would be a different
//! circuit, not a port).

pub mod brainpool;
pub mod cost_probe;
pub mod curve25519;
pub mod dialog;
pub mod driver;
pub mod mod_arith;
pub mod params;
pub mod sm2;
