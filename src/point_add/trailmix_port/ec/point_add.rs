//! EC point-add driver wired over the MBU/rfold primitives in
//! `rfold_mbu.rs`. Provides the out-of-place `ec_add_clean_out` entry
//! point and its `ec_addsub_clean_out_reverse` inverse. The in-place
//! wrapper is intentionally absent -- the y^2-linear single-inversion
//! design replaces any 2x / 4x Bennett wrap.
//!
//! The `*_deferred_w` family at the bottom of the file is the
//! working-register fallback used during the horner loop (the `_w`
//! suffix is historical — standing for "working copy preserved"; it
//! means the primitive leaves a[256] alive as an overflow bit and
//! the caller is responsible for finalizing it).

use crate::point_add::trailmix_port::circuit::{Cbit, Circuit, QReg};

fn low_pressure_creg_qload_enabled() -> bool {
    std::env::var("TRAILMIX_LOW_PRESSURE_CREG_QLOAD")
        .ok()
        .as_deref()
        != Some("0")
}

fn defer_y_materialization_enabled() -> bool {
    std::env::var("TRAILMIX_DEFER_Y_MATERIALIZE")
        .ok()
        .as_deref()
        != Some("0")
}

fn zero_dy_newdx_route_enabled() -> bool {
    std::env::var("TRAILMIX_ZERO_DY_NEWDX_ROUTE")
        .ok()
        .as_deref()
        == Some("1")
}

/// secp256k1 `R_const` = 2^32 + 977 as little-endian bytes.
#[must_use]
pub fn r_bytes() -> [u8; 32] {
    let mut r = [0u8; 32];
    r[0] = 0xD1;
    r[1] = 0x03; // 977 = 0x3D1
    r[4] = 0x01; // + 2^32
    r
}

/// Forward Horner loop: result += a * b (mod p, canonical [0, p)).
///
/// Internally uses rfold-approximate primitives whose output is in
/// [0, 2^256), then runs `reduce_once_secp256k1_from_rfold` to
/// canonicalize. The 1-qubit reduction flag is RETAINED and returned
/// to the caller — it must either be:
///   (a) consumed by the matching `horner_reverse` (which un-
///       canonicalizes first, then unwinds the rfold loop), or
///   (b) cleaned by `horner_canonical_flag_consume` if the result
///       won't be reversed (terminal output).
///
/// Returns the reduction flag (1 `QReg`).
pub fn horner_forward(circ: &mut Circuit, result: &[QReg], a: &[QReg], b: &[QReg]) -> QReg {
    let n = result.len();
    let prev = circ.current_section.clone();
    circ.set_section(&format!("{}/i={}:add", prev, n - 1));
    controlled_mod_add_deferred_w(circ, &b[n - 1], result, a);
    for i in (0..n - 1).rev() {
        circ.set_section(&format!("{prev}/i={i}:dbl"));
        mod_double_deferred_w(circ, result);
        circ.set_section(&format!("{prev}/i={i}:add"));
        controlled_mod_add_deferred_w(circ, &b[i], result, a);
    }
    // Canonicalize: result is in [0, 2^256) (rfold approximate);
    // bring it into [0, p) and retain the reduction bit.
    circ.set_section(&format!("{prev}/canon"));
    let flag = circ.alloc_qreg("horner_red_flag");
    crate::point_add::trailmix_port::rfold_mbu::reduce_once_secp256k1_from_rfold(circ, result, &flag);
    circ.set_section(&prev);
    flag
}

/// Reverse Horner loop: inverse of `horner_forward`. Consumes the
/// reduction flag emitted by the forward call, un-canonicalizes,
/// then unwinds the rfold-loop with `controlled_mod_sub_rfold_mbu`
/// and `mod_halve_rfold_mbu`. Frees `flag` at end.
pub fn horner_reverse(circ: &mut Circuit, result: &[QReg], a: &[QReg], b: &[QReg], flag: QReg) {
    let n = result.len();
    let prev = circ.current_section.clone();
    circ.set_section(&format!("{prev}/uncanon"));
    crate::point_add::trailmix_port::rfold_mbu::reduce_once_secp256k1_from_rfold_reverse(circ, result, &flag);
    circ.zero_and_free(flag);
    circ.set_section(&prev);
    for i in 0..n - 1 {
        crate::point_add::trailmix_port::rfold_mbu::controlled_mod_sub_rfold_mbu(circ, &b[i], result, a);
        crate::point_add::trailmix_port::rfold_mbu::mod_halve_rfold_mbu(circ, result);
    }
    crate::point_add::trailmix_port::rfold_mbu::controlled_mod_sub_rfold_mbu(circ, &b[n - 1], result, a);
}

/// `a += b * c (mod p)` — exact mod-p multiply-and-accumulate.
///
/// `horner_forward(result`, a, b) computes
///   result := 2^(n-1) * result + a*b   (mod p)
/// because of the Horner loop's structure. So we pre-multiply target
/// by `2^-(n-1) mod p` via raw `mod_halve_rfold_mbu` calls; the
/// subsequent `horner_forward`'s 2^(n-1) factor cancels.
///
/// Optimization vs `horner_reverse(target`, `0_REG`, `0_REG)`: `horner_reverse`
/// with all-zero inputs still emits Cuccaro gates for the no-op
/// `controlled_mod_subs` (~2.5K CCX wasted per iter × ~256 iters ≈ 640K
/// CCX wasted). Raw `mod_halve_rfold_mbu` calls skip that.
///
/// Note: `mod_halve_rfold_mbu`'s contract documents a "matching prior
/// `mod_double`" structural invariant; in practice the halve identity
/// fails with probability ~R/p ≈ 2^-224 for arbitrary canonical input
/// (a single bit mismatch on the rfold parity check). For 64 random
/// secp shots, this firing has probability ~64·2^-224 ≈ 0.
///
/// Sequence:
///   for _ in 0..n-1: `mod_halve_rfold_mbu(a)`     # a *= 2^-(n-1) mod p
///   let flag = `horner_forward(a`, b, c)            # a := `a_pre` + b*c (canonical)
///   `horner_canonical_flag_consume(a`, flag)        # free flag
pub fn mod_mac_inplace(circ: &mut Circuit, a: &[QReg], b: &[QReg], c: &[QReg]) {
    let n = a.len();
    assert_eq!(n, 257);
    assert_eq!(b.len(), n);
    assert_eq!(c.len(), n);
    let prev = circ.push_section("mod_mac");
    // Pre-multiply: a *= 2^-(n-1) via n-1 raw rfold halves.
    for _ in 0..(n - 1) {
        crate::point_add::trailmix_port::rfold_mbu::mod_halve_rfold_mbu(circ, a);
    }
    // a in rfold-approx form, value = a_pre * 2^-(n-1) mod p.
    let post_flag = horner_forward(circ, a, b, c);
    // a now canonical, value = a_pre + b*c mod p, retained flag.
    horner_canonical_flag_consume(circ, a, post_flag);
    circ.pop_section(&prev);
}

/// `a -= b * c (mod p)` — symmetric counterpart of `mod_mac_inplace`.
///
/// Sequence (mirrors `mod_mac`, swapping the order of the useful and
/// raw calls):
///   let flag = `horner_reverse(a`, b, c, `fresh_flag=0`)   # a := (a_pre-b*c)*2^-(n-1)
///   for _ in 0..n-1: `mod_double_rfold_mbu(a)`           # a *= 2^(n-1)
///   canonicalize a + free retained flag
pub fn mod_msc_inplace(circ: &mut Circuit, a: &[QReg], b: &[QReg], c: &[QReg]) {
    let n = a.len();
    assert_eq!(n, 257);
    assert_eq!(b.len(), n);
    assert_eq!(c.len(), n);
    let prev = circ.push_section("mod_msc");
    let pre_flag = circ.alloc_qreg("mod_msc_pre_flag");
    horner_reverse(circ, a, b, c, pre_flag);
    // a in rfold-approx form, value = (a_pre - b*c) * 2^-(n-1) mod p.
    // Multiply back by 2^(n-1) via raw mod_double_rfold_mbu calls,
    // then canonicalize.
    for _ in 0..(n - 1) {
        crate::point_add::trailmix_port::rfold_mbu::mod_double_rfold_mbu(circ, a);
    }
    // a in rfold-approx form, value = a_pre - b*c mod p.
    // Canonicalize with retained flag, then consume via
    // horner_canonical_flag_consume.
    let post_flag = circ.alloc_qreg("mod_msc_post_flag");
    crate::point_add::trailmix_port::rfold_mbu::reduce_once_secp256k1_from_rfold(circ, a, &post_flag);
    horner_canonical_flag_consume(circ, a, post_flag);
    circ.pop_section(&prev);
}

/// Cleans the canonicalization flag from a horner output that won't
/// be reversed (terminal output). The flag was set by
/// `reduce_once_secp256k1_from_rfold`'s `compare_geq_p` call; we
/// re-run the compare and consume via the existing
/// `compare_geq_p_secp256k1_consume` (which HMRs the flag with phase
/// correction).
///
/// Pre: `result` is canonical [0, p) (post-horner-canonicalize).
/// Post: flag is freed; `result` unchanged (still canonical).
pub fn horner_canonical_flag_consume(circ: &mut Circuit, result: &[QReg], flag: QReg) {
    // result is canonical, so a fresh compare_geq_p of result returns
    // 0. The retained `flag` value is thus equal to the output of a
    // fresh compare on result (both encode "did we need to reduce");
    // since result is now < p, the fresh compare returns 0, matching
    // flag's stale state. compare_geq_p_secp256k1_consume HMRs the
    // flag with the matching phase correction.
    crate::point_add::trailmix_port::arith::compare::compare_geq_p_secp256k1_consume(
        circ,
        &result[..result.len().min(257)],
        flag,
    );
}

/// Affine in-place point-add
/// (P=(tx,ty) -> P+Q, Q=(ox,oy) classical, preserved), where the slope inversion uses
/// the reversible shrunken-PZ divide (`shrunken_pz_state_machine::shrunken_pz_divide_forward` /
/// `shrunken_pz_divide_cancel`) -- no spooky pebbling, no `div_n/div_b` window. dx and dy stay
/// 257-bit through the divide (`shrunken_pz` needs the sign bit), so unlike the spooky path
/// there is no high-bit pop/re-push around the divides. Requires P.x != Q.x and
/// (ox - `new_x`) != 0 (generic addition; vertical/doubling excluded).
pub fn ec_add_inplace_shrunken_pz(
    circ: &mut Circuit,
    tx: &mut Vec<QReg>,
    ty: &mut Vec<QReg>,
    ox: &[Cbit],
    oy: &[Cbit],
) {
    use crate::point_add::trailmix_port::inversion::shrunken_pz_state_machine::{
        shrunken_pz_divide_cancel, shrunken_pz_divide_forward,
    };
    assert_eq!(tx.len(), 256, "tx is a 256-bit value register (P.x -> R.x)");
    assert_eq!(ty.len(), 256, "ty is a 256-bit value register (P.y -> R.y)");
    assert_eq!(ox.len(), 256);
    assert_eq!(oy.len(), 256);

    // Pad to 257-bit work registers (high overflow bit |0>) for the in-place mod
    // arithmetic; the shrunken-PZ divide needs the 257th sign bit. Both values are
    // canonical [0,p) on entry and exit, so the overflow bit is |0> and is freed
    // before return -- the public interface is 256-bit in/out.
    tx.push(circ.alloc_qreg("ec3.tx_ov"));
    ty.push(circ.alloc_qreg("ec3.ty_ov"));

    // Phase 1: ty := dy = oy - ty.
    circ.set_section("ec3.dy_build");
    if low_pressure_creg_qload_enabled() {
        mod_sub_from_creg_qload_w(circ, &ty[..], oy);
    } else {
        mod_sub_from_creg_w(circ, &ty[..], oy);
    }

    // Phase 2: tx := dx = ox - tx. (Keep 257 -- shrunken_pz divide needs the sign bit.)
    circ.set_section("ec3.dx_build");
    if low_pressure_creg_qload_enabled() {
        mod_sub_from_creg_qload_w(circ, &tx[..], ox);
    } else {
        mod_sub_from_creg_w(circ, &tx[..], ox);
    }

    // Phase 3: lambda = dy/dx (dx, dy preserved).
    circ.set_section("ec3.inv_fwd");
    let dx_inner = std::mem::take(tx);
    let dy_vec = std::mem::take(ty);
    let (dx_inner, dy_vec, lambda) = shrunken_pz_divide_forward(circ, dx_inner, dy_vec);
    *tx = dx_inner;
    *ty = dy_vec;

    if zero_dy_newdx_route_enabled() {
        // dy = lambda * dx. Zero it, then reuse ty as a qload scratch until
        // new_dy = lambda * new_dx is needed for the alt-witness cleanup.
        circ.set_section("ec3.dy_zero");
        crate::point_add::trailmix_port::rfold_mbu::mod_mul_rfold_mbu_undo(
            circ,
            &ty[..],
            &lambda,
            &tx[..],
        ); // undo the divide's dy := lambda * dx reconstruction

        circ.set_section("ec3.dx_clean");
        mod_sub_from_creg_scratch_qload_w(circ, &tx[..], ox, &ty[..]); // tx := ox - dx = tx_orig

        circ.set_section("ec3.new_x");
        mod_neg_inplace_w(circ, &tx[..]);
        mod_mac_inplace(circ, &tx[..], &lambda, &lambda);
        mod_sub_creg_scratch_qload_w(circ, &tx[..], ox, &ty[..]); // tx := new_x

        circ.set_section("ec3.alt.new_dx");
        mod_sub_from_creg_scratch_qload_w(circ, &tx[..], ox, &ty[..]); // tx := new_dx

        circ.set_section("ec3.alt.new_dy");
        mod_mac_inplace(circ, &ty[..], &lambda, &tx[..]); // ty := lambda * new_dx = new_dy
    } else {
        // Phase 4: tx := ox - dx = tx_orig.
        circ.set_section("ec3.dx_clean");
        mod_sub_from_creg_w(circ, &tx[..], ox);

        // Phase 5: tx := lambda^2 - tx_orig - ox = new_x.
        circ.set_section("ec3.new_x");
        mod_neg_inplace_w(circ, &tx[..]);
        mod_mac_inplace(circ, &tx[..], &lambda, &lambda);
        mod_sub_creg_w(circ, &tx[..], ox);

        // Phase 6: ty := dy + lambda*(tx_orig - new_x) - oy = new_y.
        // new_y = dy + lambda*(tx_orig - new_x) - oy. The intermediate
        // dx_diff = tx_orig - new_x = lambda^2 - ox - 2*new_x is computed IN PLACE in
        // tx (which holds new_x) -- NO separate 257-bit register. Its slot is exactly
        // what the qload temps reuse, so peak stays <=1050 AND the ox/oy adds become
        // O(n) (load/use/unload q-q) instead of the O(n^2) per-bit creg path.
        // tx: new_x -> dx_diff -> new_x, all exact (canonical [0,p) throughout).
        circ.set_section("ec3.new_y.dx_diff");
        mod_neg_inplace_w(circ, &tx[..]); // tx = -new_x
        mod_double_deferred_w(circ, &tx[..]); // tx = -2*new_x (rfold; mod_mac recanonicalizes)
        mod_mac_inplace(circ, &tx[..], &lambda, &lambda); // tx += lambda^2 (canonical out)
        mod_sub_creg_w(circ, &tx[..], ox); // tx -= ox => tx = dx_diff
        circ.set_section("ec3.new_y.build");
        mod_mac_inplace(circ, &ty[..], &lambda, &tx[..]); // ty += lambda*dx_diff
        if !defer_y_materialization_enabled() {
            mod_sub_creg_w(circ, &ty[..], oy); // ty = new_y
        }
        circ.set_section("ec3.new_y.dx_diff_clean");
        mod_add_creg_direct_w(circ, &tx[..], ox); // tx += ox (canonical), direct creg path
        mod_msc_inplace(circ, &tx[..], &lambda, &lambda); // tx -= lambda^2 => tx = -2*new_x (canonical)
        crate::point_add::trailmix_port::rfold_mbu::mod_halve_rfold_mbu(circ, &tx[..]); // tx = -new_x (halve of even -2*new_x)
        mod_neg_inplace_w(circ, &tx[..]); // tx = new_x

        // Phase 7: cancel lambda via the alt-witness lambda = new_dy/new_dx.
        circ.set_section("ec3.alt.new_dy");
        if !defer_y_materialization_enabled() {
            mod_add_creg_direct_w(circ, &ty[..], oy); // ty := new_y + oy = new_dy, direct creg path
        }
        circ.set_section("ec3.alt.new_dx");
        mod_sub_from_creg_w(circ, &tx[..], ox); // tx := ox - new_x = new_dx
    }
    circ.set_section("ec3.alt.cancel");
    let ndx_inner = std::mem::take(tx);
    let ndy_vec = std::mem::take(ty);
    let (ndx_inner, ndy_vec) = shrunken_pz_divide_cancel(circ, ndx_inner, ndy_vec, lambda);
    *tx = ndx_inner;
    *ty = ndy_vec;
    circ.set_section("ec3.alt.new_x_restore");
    if low_pressure_creg_qload_enabled() {
        mod_sub_from_creg_qload_w(circ, &tx[..], ox); // tx := ox - new_dx = new_x
    } else {
        mod_sub_from_creg_w(circ, &tx[..], ox); // tx := ox - new_dx = new_x
    }
    circ.set_section("ec3.alt.new_y_restore");
    if low_pressure_creg_qload_enabled() {
        mod_sub_creg_qload_w(circ, &ty[..], oy); // ty := new_dy - oy = new_y
    } else {
        mod_sub_creg_w(circ, &ty[..], oy); // ty := new_dy - oy = new_y
    }

    // Unpad: new_x/new_y are canonical [0,p), so the overflow bit is |0>. Drop it
    // to restore the 256-bit interface.
    circ.set_section("ec3.unpad");
    circ.zero_and_free(ty.pop().expect("ty padded to 257"));
    circ.zero_and_free(tx.pop().expect("tx padded to 257"));
    circ.set_section("ec3.done");
}

pub fn mod_double_deferred_w(circ: &mut Circuit, a: &[QReg]) {
    // Use rfold_mbu which HMRs the overflow bit so it doesn't
    // accumulate across Horner iterations. Tracker flags the
    // identity-based HMR (not yet verified).
    crate::point_add::trailmix_port::rfold_mbu::mod_double_rfold_mbu(circ, a);
}

/// `a -= c mod p` where `c` is a classical-bit register.
///
/// Mirrors `mod_sub_mbu`'s structure but with `add_creg`/`sub_creg` in
/// place of the QReg-QReg cuccaro add/sub. Avoids the 257-qubit
/// alloc + X-load that the caller would otherwise have to do.
/// Requires `a.len() == 257`.
pub fn mod_sub_creg_w(circ: &mut Circuit, a: &[QReg], c: &[crate::point_add::trailmix_port::circuit::Cbit]) {
    let n = a.len();
    assert_eq!(n, 257, "mod_sub_creg_w requires a.len() == 257");
    // Step 1: integer sub mod 2^n.
    crate::point_add::trailmix_port::arith::ripple_add::sub_creg(circ, a, c);
    // Step 2: flag = borrow = a[n-1].
    let flag = circ.alloc_qreg("creg.sub.flag");
    circ.cx(&a[n - 1], &flag);
    // Step 3: correction. CX flag→a[n-1] (clears the top bit when flag=1
    // since a[n-1] equals flag in the borrow case), then add p
    // (controlled_sub_const with -p = R fits in 256 bits).
    circ.cx(&flag, &a[n - 1]);
    let r = crate::point_add::trailmix_port::mod_arith::secp256k1_r_le();
    crate::point_add::trailmix_port::arith::const_add::controlled_sub_const(circ, &flag, &a[..n - 1], &r);
    // a = (a_old - c) mod p in [0, p). flag = 1[a_old < c] = 1[result + c >= p].
    // Step 4: add c back, phase-correction MBU HMRs flag against (a >= p),
    // sub c back.
    crate::point_add::trailmix_port::arith::ripple_add::add_creg(circ, a, c);
    crate::point_add::trailmix_port::arith::compare::compare_geq_p_secp256k1_phase_correction_mbu(circ, a, flag);
    crate::point_add::trailmix_port::arith::ripple_add::sub_creg(circ, a, c);
}

/// `a += c mod p` where `c` is a classical-bit register.
///
/// Exact mirror of `mod_sub_creg_w` with `add_creg`/`sub_creg` swapped
/// in the phase-correction bracket. Zero temp registers — only a single
/// 1-qubit flag, HMR-freed by the operand-free phase correction.
/// Requires `a.len() == 257`, `a` (and the value `c`) in `[0, p)`.
///
/// 1. `add_creg(a`, c)                    a = `a_old` + c (mod 2^257), in [0, 2p)
/// 2. flag = 1[a >= p]                  = 1[`a_old` + c >= p]
/// 3. if flag: a -= p                   a = (`a_old` + c) mod p, in [0, p)
/// 4. `sub_creg(a`, c)                    a = result - c (mod 2^257)
/// 5. phase-correction MBU: 1[a>=p]==flag, HMR-frees flag.
///    (flag=0 → a = `a_old` in [0,p), a<p;  flag=1 → a = `a_old` - p (mod
///     2^257) in [2^256+R, 2^257), a>=p — so the identity holds and
///     `compare_geq_p_secp256k1` (a correct general 257-bit comparator)
///     verifies it.)
/// 6. `add_creg(a`, c)                    a = result = (`a_old` + c) mod p
pub fn mod_add_creg_direct_w(circ: &mut Circuit, a: &[QReg], c: &[crate::point_add::trailmix_port::circuit::Cbit]) {
    let n = a.len();
    assert_eq!(n, 257, "mod_add_creg_direct_w requires a.len() == 257");
    // Step 1: integer add mod 2^n.
    crate::point_add::trailmix_port::arith::ripple_add::add_creg(circ, a, c);
    // Step 2: flag = (a >= p).
    let flag = circ.alloc_qreg("creg.add.flag");
    crate::point_add::trailmix_port::arith::compare::compare_geq_p_secp256k1(circ, a, &flag);
    // Step 3: if flag, a -= p (p as 33-byte LE constant, bit 256 = 0).
    let mut p_le33 = [0u8; 33];
    p_le33[..32].copy_from_slice(&crate::point_add::trailmix_port::mod_arith::SECP256K1_P_LE);
    crate::point_add::trailmix_port::arith::const_add::controlled_sub_const(circ, &flag, a, &p_le33);
    // a = (a_old + c) mod p in [0, p). flag = 1[a_old + c >= p].
    // Step 4-6: sub c, phase-correction MBU HMRs flag against (a >= p),
    // add c back.
    crate::point_add::trailmix_port::arith::ripple_add::sub_creg(circ, a, c);
    crate::point_add::trailmix_port::arith::compare::compare_geq_p_secp256k1_phase_correction_mbu(circ, a, flag);
    crate::point_add::trailmix_port::arith::ripple_add::add_creg(circ, a, c);
}

/// `a := c - a mod p` where `c` is a classical-bit register.
///
/// Same shape as `mod_sub_from_w`: negate then add. Uses the zero-temp direct
/// creg add to avoid allocating a 257-qubit qload temp in high-pressure phases.
pub fn mod_sub_from_creg_w(circ: &mut Circuit, a: &[QReg], c: &[crate::point_add::trailmix_port::circuit::Cbit]) {
    assert_eq!(a.len(), 257, "mod_sub_from_creg_w requires a.len() == 257");
    mod_neg_inplace_w(circ, a);
    mod_add_creg_direct_w(circ, a, c);
}

/// Overwrite `a` in place with `-a (mod p)`, output in `[0, p)`.
///
/// Unlike `mod_neg_clean` (which produces `p` when `a = 0` and so
/// can't be chained into primitives that assume inputs in
/// `[0, p)`), this uses the polylog identity p - x ≡ ~x + (p+1)
/// (mod 2^n) for n = 257 and x in [0, p). All ops are polylog-anc.
pub fn mod_neg_inplace_w(circ: &mut Circuit, a: &[QReg]) {
    assert_eq!(a.len(), 257);
    // Step 1: bit-flip all 257 bits → ~x.
    for q in a {
        circ.x(q);
    }
    // Step 2: add (p + 1) as a 33-byte LE constant. p + 1 = 0x30
    // followed by 30 bytes of 0xFF, byte 32 = 0. Highly structured
    // (one run + low-bit difference) — controlled_add_const's
    // runs-based path handles this in ~770 CCX (vs cuccaro's ~2.6K
    // for quantum-quantum).
    let mut p_plus_1: [u8; 33] = [0u8; 33];
    p_plus_1[..32].copy_from_slice(&crate::point_add::trailmix_port::mod_arith::SECP256K1_P_LE);
    p_plus_1[0] = p_plus_1[0].wrapping_add(1); // 0x2F -> 0x30
    crate::point_add::trailmix_port::arith::ripple_add::add_const(circ, a, &p_plus_1);
}

/// Load classical `c` (256 bits) into a fresh 257-qubit quantum temp (bit 256 =
/// |0>) via `x_if_bit` (0 Toffoli). Treats `c` (an EC-add input coordinate) as a
/// quantum operand -- NOT a constant. Caller unloads with `unload_creg_temp`.
fn load_creg_temp(circ: &mut Circuit, c: &[crate::point_add::trailmix_port::circuit::Cbit]) -> Vec<QReg> {
    let temp: Vec<QReg> = (0..257).map(|_| circ.alloc_qreg("creg.qload")).collect();
    for (i, b) in c.iter().enumerate().take(256) {
        circ.x_if_bit(&temp[i], *b);
    }
    temp
}

fn unload_creg_temp(circ: &mut Circuit, temp: Vec<QReg>, c: &[crate::point_add::trailmix_port::circuit::Cbit]) {
    for (i, b) in c.iter().enumerate().take(256) {
        circ.x_if_bit(&temp[i], *b);
    }
    for q in temp {
        circ.zero_and_free(q);
    }
}

/// `a += c (mod p)` via LOAD c into a quantum temp -> one O(n) q-q `mod_add_mbu`
/// -> UNLOAD. Replaces `mod_add_creg_w`'s O(n^2) per-bit `add_creg` (~n controlled
/// increments). +257 peak: use only where the section has >= ~257q headroom.
pub fn mod_add_creg_w(circ: &mut Circuit, a: &[QReg], c: &[crate::point_add::trailmix_port::circuit::Cbit]) {
    let temp = load_creg_temp(circ, c);
    crate::point_add::trailmix_port::mod_arith::mod_add_mbu(circ, a, &temp, &crate::point_add::trailmix_port::mod_arith::SECP256K1_P_LE);
    unload_creg_temp(circ, temp, c);
}

/// `a -= c (mod p)` -- q-q load/use/unload version of `mod_sub_creg_w`.
pub fn mod_sub_creg_qload_w(circ: &mut Circuit, a: &[QReg], c: &[crate::point_add::trailmix_port::circuit::Cbit]) {
    let temp = load_creg_temp(circ, c);
    crate::point_add::trailmix_port::mod_arith::mod_sub_mbu(circ, a, &temp, &crate::point_add::trailmix_port::mod_arith::SECP256K1_P_LE);
    unload_creg_temp(circ, temp, c);
}

/// `a := c - a (mod p)` -- q-q load/use/unload version of `mod_sub_from_creg_w`.
pub fn mod_sub_from_creg_qload_w(circ: &mut Circuit, a: &[QReg], c: &[crate::point_add::trailmix_port::circuit::Cbit]) {
    mod_neg_inplace_w(circ, a);
    mod_add_creg_w(circ, a, c);
}

fn load_creg_into_scratch(circ: &mut Circuit, scratch: &[QReg], c: &[crate::point_add::trailmix_port::circuit::Cbit]) {
    assert_eq!(scratch.len(), 257, "creg scratch must be 257 bits");
    for (i, b) in c.iter().enumerate().take(256) {
        circ.x_if_bit(&scratch[i], *b);
    }
}

fn unload_creg_from_scratch(circ: &mut Circuit, scratch: &[QReg], c: &[crate::point_add::trailmix_port::circuit::Cbit]) {
    for (i, b) in c.iter().enumerate().take(256) {
        circ.x_if_bit(&scratch[i], *b);
    }
}

pub fn mod_add_creg_scratch_qload_w(
    circ: &mut Circuit,
    a: &[QReg],
    c: &[crate::point_add::trailmix_port::circuit::Cbit],
    scratch: &[QReg],
) {
    load_creg_into_scratch(circ, scratch, c);
    crate::point_add::trailmix_port::mod_arith::mod_add_mbu(
        circ,
        a,
        scratch,
        &crate::point_add::trailmix_port::mod_arith::SECP256K1_P_LE,
    );
    unload_creg_from_scratch(circ, scratch, c);
}

pub fn mod_sub_creg_scratch_qload_w(
    circ: &mut Circuit,
    a: &[QReg],
    c: &[crate::point_add::trailmix_port::circuit::Cbit],
    scratch: &[QReg],
) {
    load_creg_into_scratch(circ, scratch, c);
    crate::point_add::trailmix_port::mod_arith::mod_sub_mbu(
        circ,
        a,
        scratch,
        &crate::point_add::trailmix_port::mod_arith::SECP256K1_P_LE,
    );
    unload_creg_from_scratch(circ, scratch, c);
}

pub fn mod_sub_from_creg_scratch_qload_w(
    circ: &mut Circuit,
    a: &[QReg],
    c: &[crate::point_add::trailmix_port::circuit::Cbit],
    scratch: &[QReg],
) {
    mod_neg_inplace_w(circ, a);
    mod_add_creg_scratch_qload_w(circ, a, c, scratch);
}

pub fn controlled_mod_add_deferred_w(circ: &mut Circuit, ctrl: &QReg, a: &[QReg], b: &[QReg]) {
    let n = a.len();
    assert_eq!(b.len(), n);
    crate::point_add::trailmix_port::rfold_mbu::controlled_mod_add_rfold_mbu(circ, ctrl, a, b);
}

#[cfg(test)]
mod tests {
    use super::{mod_mac_inplace, mod_msc_inplace};
    use crate::point_add::trailmix_port::circuit::{Circuit, QReg};
    use alloy_primitives::U256;
    use crate::point_add::trailmix_port::num_bigint::BigUint;
    use rand::RngCore;
    use zkp_ecc_lib::WeierstrassEllipticCurve;

    fn secp256k1() -> WeierstrassEllipticCurve {
        WeierstrassEllipticCurve {
            modulus: U256::from_str_radix(
                "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F",
                16,
            )
            .unwrap(),
            a: U256::from(0u64),
            b: U256::from(7u64),
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

    /// Verify mod_mac_inplace: a += b * c (mod p) on 16 random shots.
    /// All inputs canonical [0, p) on entry, output canonical on exit.
    #[test]
    #[ignore]
    fn mod_mac_inplace_random_secp() {
        use rand::RngCore;
        let p = BigUint::from_bytes_le(&crate::point_add::trailmix_port::mod_arith::SECP256K1_P_LE);
        let mut rng = rand::thread_rng();

        let mut circ = Circuit::new();
        let a = circ.alloc_qreg_bits("a", 257);
        let b = circ.alloc_qreg_bits("b", 257);
        let c = circ.alloc_qreg_bits("c", 257);

        let mut shots: Vec<(BigUint, BigUint, BigUint)> = Vec::with_capacity(16);
        for shot in 0..16 {
            let av = {
                let mut bs = [0u8; 32];
                rng.fill_bytes(&mut bs);
                BigUint::from_bytes_le(&bs) % &p
            };
            let bv = {
                let mut bs = [0u8; 32];
                rng.fill_bytes(&mut bs);
                BigUint::from_bytes_le(&bs) % &p
            };
            let cv = {
                let mut bs = [0u8; 32];
                rng.fill_bytes(&mut bs);
                BigUint::from_bytes_le(&bs) % &p
            };
            let to_bytes = |v: &BigUint| {
                let mut bs = v.to_bytes_le();
                bs.resize(32, 0);
                bs
            };
            circ.sim_load_reg_bytes_shot(&a[..256], &to_bytes(&av), shot);
            circ.sim_load_reg_bytes_shot(&b[..256], &to_bytes(&bv), shot);
            circ.sim_load_reg_bytes_shot(&c[..256], &to_bytes(&cv), shot);
            shots.push((av, bv, cv));
        }

        mod_mac_inplace(&mut circ, &a, &b, &c);

        let mut outs: Vec<QReg> = Vec::new();
        outs.extend(a);
        outs.extend(b);
        outs.extend(c);
        let (sim, det) = circ.destroy_sim(outs);
        let (a_d, rest) = det.split_at(257);
        let (b_d, c_d) = rest.split_at(257);
        for (shot, (av, bv, cv)) in shots.iter().enumerate() {
            let got_a = BigUint::from_bytes_le(&sim.read_bytes_shot(&a_d[..256], shot));
            let got_b = BigUint::from_bytes_le(&sim.read_bytes_shot(&b_d[..256], shot));
            let got_c = BigUint::from_bytes_le(&sim.read_bytes_shot(&c_d[..256], shot));
            let expected = (av + bv * cv) % &p;
            let a_bit256 = sim.read_bytes_shot(&a_d[256..257], shot)[0] & 1;
            assert_eq!(got_b, *bv, "shot {shot}: b mutated");
            assert_eq!(got_c, *cv, "shot {shot}: c mutated");
            assert_eq!(
                got_a, expected,
                "shot {shot}: a != a_pre + b*c mod p (expected {expected}, got {got_a})"
            );
            assert_eq!(a_bit256, 0, "shot {shot}: a bit 256 non-zero");
        }
    }

    /// Verify mod_msc_inplace: a -= b * c (mod p) on 16 random shots.
    #[test]
    #[ignore]
    fn mod_msc_inplace_random_secp() {
        use rand::RngCore;
        let p = BigUint::from_bytes_le(&crate::point_add::trailmix_port::mod_arith::SECP256K1_P_LE);
        let mut rng = rand::thread_rng();

        let mut circ = Circuit::new();
        let a = circ.alloc_qreg_bits("a", 257);
        let b = circ.alloc_qreg_bits("b", 257);
        let c = circ.alloc_qreg_bits("c", 257);

        let mut shots: Vec<(BigUint, BigUint, BigUint)> = Vec::with_capacity(16);
        for shot in 0..16 {
            let av = {
                let mut bs = [0u8; 32];
                rng.fill_bytes(&mut bs);
                BigUint::from_bytes_le(&bs) % &p
            };
            let bv = {
                let mut bs = [0u8; 32];
                rng.fill_bytes(&mut bs);
                BigUint::from_bytes_le(&bs) % &p
            };
            let cv = {
                let mut bs = [0u8; 32];
                rng.fill_bytes(&mut bs);
                BigUint::from_bytes_le(&bs) % &p
            };
            let to_bytes = |v: &BigUint| {
                let mut bs = v.to_bytes_le();
                bs.resize(32, 0);
                bs
            };
            circ.sim_load_reg_bytes_shot(&a[..256], &to_bytes(&av), shot);
            circ.sim_load_reg_bytes_shot(&b[..256], &to_bytes(&bv), shot);
            circ.sim_load_reg_bytes_shot(&c[..256], &to_bytes(&cv), shot);
            shots.push((av, bv, cv));
        }

        mod_msc_inplace(&mut circ, &a, &b, &c);

        let mut outs: Vec<QReg> = Vec::new();
        outs.extend(a);
        outs.extend(b);
        outs.extend(c);
        let (sim, det) = circ.destroy_sim(outs);
        let (a_d, rest) = det.split_at(257);
        let (_b_d, _c_d) = rest.split_at(257);
        for (shot, (av, bv, cv)) in shots.iter().enumerate() {
            let got_a = BigUint::from_bytes_le(&sim.read_bytes_shot(&a_d[..256], shot));
            let expected = (av + &p - (bv * cv) % &p) % &p;
            let a_bit256 = sim.read_bytes_shot(&a_d[256..257], shot)[0] & 1;
            assert_eq!(
                got_a, expected,
                "shot {shot}: a != a_pre - b*c mod p (expected {expected}, got {got_a})"
            );
            assert_eq!(a_bit256, 0, "shot {shot}: a bit 256 non-zero");
        }
    }

    /// Full in-place P+Q via the shrunken-PZ divide (`ec_add_inplace_shrunken_pz`).
    /// Exercises BOTH shrunken_pz divides together: the forward slope lambda=dy/dx and the
    /// alt-witness cancel lambda=new_dy/new_dx. 64 RANDOM secp points, NO schedule
    /// prefilter -- the schedule must handle random |dx|/|new_dx| (any miss is a
    /// schedule bug, not a skipped input). Checks new_x==R.x, new_y==R.y on all 64
    /// shots + phase-clean. 256-bit in/out.
    #[test]
    fn ec_add_inplace_shrunken_pz_random_64() {
        use crate::point_add::trailmix_port::circuit::Cbit;
        let curve = secp256k1();
        let mut rng = rand::thread_rng();
        let to_big = |u: U256| BigUint::from_bytes_le(&u.to_le_bytes::<32>());
        // RANDOM secp points -- NO schedule prefilter. If the shrunken-PZ schedule
        // can't handle a random |dx| / |new_dx| with high probability that is a
        // SCHEDULE BUG to fix, not a test input to skip. The only acceptable miss is
        // the ~2^-19 Shor tail.
        let mut cases: Vec<(U256, U256, U256, U256, U256, U256)> = Vec::with_capacity(64);
        while cases.len() < 64 {
            let draw = |rng: &mut rand::rngs::ThreadRng| -> U256 {
                U256::from(rng.next_u64())
                    ^ (U256::from(rng.next_u64()) << 64)
                    ^ (U256::from(rng.next_u64()) << 128)
                    ^ (U256::from(rng.next_u64()) << 192)
            };
            let (s_p, s_q) = (draw(&mut rng), draw(&mut rng));
            if s_p == U256::ZERO || s_q == U256::ZERO || s_p == s_q {
                continue;
            }
            let pp = curve.mul(curve.gx, curve.gy, s_p);
            let qq = curve.mul(curve.gx, curve.gy, s_q);
            if pp.0 == qq.0 {
                continue; // generic-add precondition P.x != Q.x (not a schedule filter)
            }
            let r = curve.add(pp.0, pp.1, qq.0, qq.1);
            cases.push((pp.0, pp.1, qq.0, qq.1, r.0, r.1));
        }

        let mut circ = Circuit::new();
        circ.set_max_qubit_peak(1300); // shrunken-PZ peak (re-measured after schedule fix)
        circ.set_section("ec3_test");
        let mut tx: Vec<QReg> = (0..256)
            .map(|i| circ.alloc_qreg(&format!("tx[{i}]")))
            .collect();
        let mut ty: Vec<QReg> = (0..256)
            .map(|i| circ.alloc_qreg(&format!("ty[{i}]")))
            .collect();
        let ox: Vec<Cbit> = (0..256).map(|_| circ.alloc_input_bit()).collect();
        let oy: Vec<Cbit> = (0..256).map(|_| circ.alloc_input_bit()).collect();
        let mut rs = Vec::with_capacity(64);
        for (shot, (px, py, qx, qy, rx, ry)) in cases.iter().enumerate() {
            circ.sim_load_reg_bytes_shot(&tx[..256], &px.to_le_bytes::<32>(), shot);
            circ.sim_load_reg_bytes_shot(&ty[..256], &py.to_le_bytes::<32>(), shot);
            circ.sim_load_bits_bytes_shot(&ox, &qx.to_le_bytes::<32>(), shot);
            circ.sim_load_bits_bytes_shot(&oy, &qy.to_le_bytes::<32>(), shot);
            rs.push((to_big(*rx), to_big(*ry)));
        }

        super::ec_add_inplace_shrunken_pz(&mut circ, &mut tx, &mut ty, &ox, &oy);
        let peak = circ.peak_qubits;

        {
            let (tx_r, ty_r, rsc) = (&tx, &ty, rs.clone());
            circ.contract_check("ec3_result", move |view, shot| {
                let rd = |reg: &[QReg]| -> BigUint {
                    let mut a = BigUint::from(0u32);
                    for j in 0..256 {
                        if view.contract_read_bit_shot(&reg[j], shot) {
                            a |= BigUint::from(1u32) << j;
                        }
                    }
                    a
                };
                let gx = rd(tx_r);
                let gy = rd(ty_r);
                if gx != rsc[shot].0 {
                    return Err(format!(
                        "shot {} new_x wrong: got {:x} want {:x}",
                        shot, gx, rsc[shot].0
                    ));
                }
                if gy != rsc[shot].1 {
                    return Err(format!(
                        "shot {} new_y wrong: got {:x} want {:x}",
                        shot, gy, rsc[shot].1
                    ));
                }
                Ok(())
            });
        }
        circ.assert_phase_clean();
        eprintln!(
            "ec_add_shrunken_pz: peak={} tof={} ops={}",
            peak,
            circ.executed_toffoli_shots / 64,
            circ.total_ops()
        );
        let mut outs = vec![];
        outs.extend(tx);
        outs.extend(ty);
        let _ = circ.destroy_sim(outs);
    }
}
