//! Safegcd hybrid entry point.
//!
//! Provides the SAFEGCD=1 env-var hook into the point-addition build.
//! Uses dialog GCD's proven Kaliski quotient internally.
//!
//! # Why hybrid and not pure safegcd
//!
//! Pure safegcd tracks the modular inverse through d/e coefficients using
//! the NOT+add+1 (Cuccaro) pattern: e = NOT(e) + d + 1 = d - e (mod 2^NX).
//!
//! The fundamental problem: Cuccaro arithmetic is mod 2^NX, but the
//! modular inverse tracking requires mod p arithmetic. When d < e, the
//! Cuccaro wraps: 2^NX + d - e. This 2^NX term does NOT reduce to 0 mod p
//! (2^NX mod p != 0), so the offset is irreducible. Over 741 iterations,
//! these offsets accumulate and corrupt the d/e coefficients.
//!
//! Multiple attempted fixes (all failed 9024/9024):
//! - NX enlargement (270 bits): 2^NX mod p still != 0, offset still accumulates
//! - Arithmetic shift: same issue, 2's complement != mod p
//! - Borrow capture via cuccaro_sub: borrow fixes wrap but corrected value
//!   still carries mod-2^NX artifacts, not mod-p purity
//! - Negate-before-dm (v2): fragile workaround, 6 extra Cuccaro ops per iteration
//!
//! The dialog GCD's Kaliski algorithm avoids this entirely - it uses ONLY
//! addition and bit shifts (no subtraction), so Cuccaro wrapping never occurs.

use super::*;
use crate::point_add::N;
use alloy_primitives::U256;

const DEFAULT_SAFEGCD_ITERS: usize = 741;

pub(crate) fn safegcd_iters() -> usize {
    std::env::var("SAFEGCD_ITERS").ok().and_then(|s| s.parse().ok()).unwrap_or(DEFAULT_SAFEGCD_ITERS)
}
pub(crate) fn safegcd_enabled() -> bool { std::env::var("SAFEGCD").ok().as_deref() == Some("1") }

pub(crate) fn emit_safegcd_raw_pa(
    b: &mut B, tx: &[QubitId], ty: &[QubitId], ox: &[BitId], oy: &[BitId], p: U256,
) {
    // HYBRID: Proven dialog GCD quotient (no Cuccaro wrapping issue)
    emit_dialog_gcd_raw_quotient(b, tx, ty, p);

    round84_emit_fused_square_xtail(b, tx, ty, ox, p);

    if dialog_fuse_c_form_enabled() {
        mod_add_triple_qb(b, tx, ox, p);
    } else {
        mod_sub_qb(b, tx, ox, p);
        mod_neg_inplace_fast(b, tx, p);
    }

    emit_dialog_gcd_raw_ipmul(b, tx, ty, p);

    mod_sub_qb(b, ty, oy, p);

    if dialog_fuse_x_restore_enabled() {
        mod_const_minus_reg_qb(b, tx, ox, p);
    } else {
        mod_neg_inplace_fast(b, tx, p);
        mod_add_qb(b, tx, ox, p);
    }
}
