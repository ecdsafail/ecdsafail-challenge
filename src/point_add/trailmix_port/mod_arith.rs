//! Exact modular arithmetic for secp256k1.
//!
//! Replaces the approximate rfold-based `poc_arith::{mod_add`, `mod_sub`,
//! `mod_mul`} with primitives that correctly reduce all results into
//! [0, p). Uses `compare_geq_const` + `controlled_sub_const` with a
//! caller-managed flag ancilla for reversibility.
//!
//! ## Design
//!
//! Each forward primitive takes a `flag` qubit (|0> on entry) that
//! records "did the reduction fire?" The caller holds onto the flag
//! until the reverse pass, which consumes it via self-inverse
//! `compare_geq_const`. This gives exact bidirectional reduction with
//! zero selfwire.
//!
//! Registers are 257 bits wide (bit 256 = overflow slot). Values are
//! maintained in [0, p) with a[256] = 0 after every primitive.
//!
//! All code is physical-only: no selfwire, no rfold, no R-on-non-zero.

use crate::point_add::trailmix_port::circuit::{Circuit, QReg};

/// secp256k1 prime p = 2^256 - 2^32 - 977, little-endian 32 bytes.
pub const SECP256K1_P_LE: [u8; 32] = [
    0x2F, 0xFC, 0xFF, 0xFF, 0xFE, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
];

/// Controlled `mod_add`: if ctrl=1, a += b mod p; else no-op.
///
/// Uses single-compare pattern: flag = (`a_post` >= p). When ctrl=0,
/// a is unchanged (still in [0, p)), so compare gives 0 → flag = 0,
/// and the sub is a no-op. When ctrl=1, a = `a_pre` + b ∈ [0, 2p),
/// compare gives (a >= p), flag records.
pub fn controlled_mod_add(
    circ: &mut Circuit,
    ctrl: &QReg,
    a: &[QReg],
    b: &[QReg],
    p_bytes: &[u8; 32],
    flag: &QReg,
) {
    crate::point_add::trailmix_port::arith::ripple_add::controlled_add(circ, ctrl, a, b);
    crate::point_add::trailmix_port::arith::compare::compare_geq_const(circ, a, p_bytes, flag);
    crate::point_add::trailmix_port::arith::const_add::controlled_sub_const(circ, flag, a, p_bytes);
}

pub fn controlled_mod_add_reverse(
    circ: &mut Circuit,
    ctrl: &QReg,
    a: &[QReg],
    b: &[QReg],
    p_bytes: &[u8; 32],
    flag: &QReg,
) {
    // Undo sub p: restores a to (a_pre + b) (post-add, pre-reduction).
    crate::point_add::trailmix_port::arith::const_add::controlled_add_const(circ, flag, a, p_bytes);
    // Self-inverse compare on unchanged a → flag back to 0.
    crate::point_add::trailmix_port::arith::compare::compare_geq_const(circ, a, p_bytes, flag);
    // Undo controlled integer add.
    crate::point_add::trailmix_port::arith::ripple_add::controlled_sub(circ, ctrl, a, b);
}

/// Controlled `mod_sub`: if ctrl=1, a -= b mod p; else no-op.
pub fn controlled_mod_sub(
    circ: &mut Circuit,
    ctrl: &QReg,
    a: &[QReg],
    b: &[QReg],
    p_bytes: &[u8; 32],
    flag: &QReg,
) {
    let n = a.len();
    crate::point_add::trailmix_port::arith::ripple_add::controlled_sub(circ, ctrl, a, b);
    // flag = ctrl AND (borrow bit).
    circ.ccx(ctrl, &a[n - 1], flag);
    // Add p if flag.
    crate::point_add::trailmix_port::arith::const_add::controlled_add_const(circ, flag, a, p_bytes);
}

/// Modular halving: a := a/2 mod p. Uses `parity_flag` as the
/// "was odd" indicator (caller-managed).
///
/// Pre: a in [0, p), a[256] = 0, `parity_flag` = |0>.
/// Post: a in [0, p), a[256] = 0, `parity_flag` = `a_pre`[0].
pub fn mod_halve(circ: &mut Circuit, a: &[QReg], p_bytes: &[u8; 32], parity_flag: &QReg) {
    // Record parity.
    circ.cx(&a[0], parity_flag);
    // If odd, add p (making a even).
    crate::point_add::trailmix_port::arith::const_add::controlled_add_const(circ, parity_flag, a, p_bytes);
    // Right shift. a is now (a + p*parity) / 2.
    crate::point_add::trailmix_port::arith::shift::right_shift(circ, a);
    // Post: a[256] = 0 (right_shift fills the top with 0).
}

/// Modular doubling: a := 2a mod p. Reuses `mod_add` with a itself.
/// Uses one flag ancilla.
pub fn mod_double(circ: &mut Circuit, a: &[QReg], p_bytes: &[u8; 32], flag: &QReg) {
    let p_val_pre = crate::point_add::trailmix_port::num_bigint::BigUint::from_bytes_le(p_bytes);
    {
        let a_for_capture: Vec<&QReg> = a.iter().collect();
        let p_val = p_val_pre.clone();
        circ.contract_capture(
            "mod_arith.mod_double",
            move |view, shot| -> Result<crate::point_add::trailmix_port::num_bigint::BigUint, String> {
                let mut v = crate::point_add::trailmix_port::num_bigint::BigUint::from(0u32);
                for (i, q) in a_for_capture.iter().enumerate() {
                    if view.contract_read_bit_shot(q, shot) {
                        v |= crate::point_add::trailmix_port::num_bigint::BigUint::from(1u32) << i;
                    }
                }
                if v >= p_val {
                    return Err(format!("a_pre = {v:#x} >= p"));
                }
                Ok(v)
            },
        );
    }
    // left_shift is exact doubling (bit 256 = old bit 255).
    // This could overshoot [0, p); reduce.
    crate::point_add::trailmix_port::arith::shift::left_shift(circ, a);
    // Now a holds 2*a_pre as a 257-bit value. a[256] = old bit 255.
    // Reduce: if a >= p, sub p.
    crate::point_add::trailmix_port::arith::compare::compare_geq_const(circ, a, p_bytes, flag);
    crate::point_add::trailmix_port::arith::const_add::controlled_sub_const(circ, flag, a, p_bytes);
    {
        let a_for_check: Vec<&QReg> = a.iter().collect();
        let p_val = p_val_pre;
        circ.contract_pop_and_check::<crate::point_add::trailmix_port::num_bigint::BigUint, _>(
            "mod_arith.mod_double",
            move |a_pre, view, shot| -> Result<(), String> {
                let mut a_post = crate::point_add::trailmix_port::num_bigint::BigUint::from(0u32);
                for (i, q) in a_for_check.iter().enumerate() {
                    if view.contract_read_bit_shot(q, shot) {
                        a_post |= crate::point_add::trailmix_port::num_bigint::BigUint::from(1u32) << i;
                    }
                }
                let expected = (a_pre * crate::point_add::trailmix_port::num_bigint::BigUint::from(2u32)) % &p_val;
                if a_post != expected {
                    return Err(format!(
                        "shot {shot}: a_post = {a_post:#x}, expected 2*a_pre mod p = {expected:#x}"
                    ));
                }
                Ok(())
            },
        );
    }
}

pub fn mod_double_reverse(circ: &mut Circuit, a: &[QReg], p_bytes: &[u8; 32], flag: &QReg) {
    crate::point_add::trailmix_port::arith::const_add::controlled_add_const(circ, flag, a, p_bytes);
    crate::point_add::trailmix_port::arith::compare::compare_geq_const(circ, a, p_bytes, flag);
    crate::point_add::trailmix_port::arith::shift::right_shift(circ, a);
}

// =====================================================================
// secp256k1-hardcoded mod_arith: exploits p = 2^256 - R structure.
//
// Key optimizations:
// - compare uses inline constant (7 zero bits in p → cheap carry chain)
// - controlled_sub_p = controlled_add(neg_p) where neg_p has 8 set bits
// - controlled_add_p = CX(ctrl, a[256]) + controlled_sub_R(ctrl, a[0..255])
//   where R has 7 set bits
// =====================================================================

// =====================================================================
// MBU mod_arith: no persistent flags, no reversal needed.
// Uses Lemma 4.1 from Luongo et al. (arXiv:2407.20167):
// The reduction flag 1[x+a >= p] equals 1[(x+a mod p) < a].
// The phase correction computes the EQUIVALENT comparison on
// the POST-reduction data, so no reversal is needed.
// =====================================================================

/// a += b mod p. MBU: flag is HMR'd immediately with phase
/// correction via `compare_less(a_reduced`, b). No persistent flag.
pub fn mod_add_mbu(circ: &mut Circuit, a: &[QReg], b: &[QReg], _p_bytes: &[u8; 32]) {
    // Step 1: integer add
    crate::point_add::trailmix_port::arith::ripple_add::add(circ, a, b);
    // Step 2: compare a >= p, store in flag
    let flag = circ.alloc_qreg("flag");
    crate::point_add::trailmix_port::arith::compare::compare_geq_p_secp256k1(circ, a, &flag);
    // Step 3: controlled sub p (flag as separate qubit, never modified)
    controlled_add_neg_p_secp256k1(circ, &flag, a);
    // a is now (a_old + b) mod p. flag = 1[a_old+b >= p] = 1[result < b]
    // (Lemma 4.1). MBU-verified compare_lt does HMR(flag) + declare,
    // and consumes/frees `flag` internally.
    crate::point_add::trailmix_port::arith::compare::compare_lt_phase_correction_mbu(circ, a, b, &flag);
}

/// a -= b mod p. MBU: flag HMR'd with phase correction.
/// Identity: 1[`a_old` < b] = 1[(a-b mod p) + b >= p].
/// Phase correction: temporarily add b to result, compare >= p.
pub fn mod_sub_mbu(circ: &mut Circuit, a: &[QReg], b: &[QReg], _p_bytes: &[u8; 32]) {
    let n = a.len();
    // Step 1: integer sub
    crate::point_add::trailmix_port::arith::ripple_add::sub(circ, a, b);
    // Step 2: flag = borrow = a[n-1]
    let flag = circ.alloc_qreg("flag");
    circ.cx(&a[n - 1], &flag);
    // Step 3: add p if borrow (correction)
    circ.cx(&flag, &a[n - 1]);
    let r = secp256k1_r_le();
    crate::point_add::trailmix_port::arith::const_add::controlled_sub_const(circ, &flag, &a[..n - 1], &r);
    // a = (a_old - b) mod p now. flag = 1[a_old < b] = 1[result + b >= p].
    // Temporarily add b, MBU compare >= p, undo. compare_geq_mbu handles
    // HMR(flag) + declare_identity internally.
    crate::point_add::trailmix_port::arith::ripple_add::add(circ, a, b); // a = result + b
                                               // compare_geq_p_secp256k1_phase_correction_mbu consumes and frees flag
    crate::point_add::trailmix_port::arith::compare::compare_geq_p_secp256k1_phase_correction_mbu(circ, a, flag);
    crate::point_add::trailmix_port::arith::ripple_add::sub(circ, a, b); // restore: a = result
}

#[must_use]
pub fn secp256k1_r_le() -> [u8; 32] {
    let mut r = [0u8; 32];
    r[0] = 0xD1;
    r[1] = 0x03;
    r[4] = 0x01;
    r
}

fn controlled_add_neg_p_secp256k1(circ: &mut Circuit, ctrl: &QReg, a: &[QReg]) {
    assert_eq!(a.len(), 257);
    // -p = R - 2^256. Add sparse R into the full 257-bit register so the
    // carry lands in a[256], then toggle the 2^256 term.
    let r = secp256k1_r_le();
    crate::point_add::trailmix_port::arith::const_add::controlled_add_const(circ, ctrl, a, &r);
    circ.cx(ctrl, &a[256]);
}

// ============================================================
// Polylog-ancilla classical-constant mod-p add/sub.
//
