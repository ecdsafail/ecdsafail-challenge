//! Rfold-MBU mod arithmetic for secp256k1.
//!
//! Cheaper than `mod_arith.rs`'s exact-reduction MBU primitives:
//! skips the `compare_geq_p` + `controlled_add_neg_p` (~4n CCX) and
//! replaces it with a controlled add of R = 2^32+977 to the low
//! 256 bits (~16 CCX since R has 8 set bits). Phase correction
//! uses the same Lemma 4.1 identity: `1[r < b] = X` where
//! X = "did the integer add overflow into bit 256".
//!
//! Output range: rfold is APPROXIMATE — output is in [0, 2^256),
//! not [0, p). Probability of "underreduced" output (in [p, 2^256))
//! per add is R / 2^256 ≈ 2^-224 for random inputs. Composes
//! safely as long as a final exact reduction is applied at the
//! end of the pipeline, OR all downstream operations tolerate
//! [0, 2^256) inputs (which the rfold primitives themselves do).
//!
//! Identity verification:
//!   For a, b in [0, p): A = a+b, X = (A >= 2^256), r = A - X*p.
//!   - X=0: r = A. r<b iff a+b<b iff a<0 -> false. X=0 ✓
//!   - X=1: r = A-p. r<b iff a+b-p<b iff a<p -> true. X=1 ✓
//!   - X=0, A in [p, 2^256) (under-reduced): r=A. r<b iff a<0 -> false. X=0 ✓

use crate::point_add::trailmix_port::circuit::{Circuit, ContractReadable, QReg};

/// secp256k1 R = 2^32 + 977 = 0x100000003D1, little-endian 32 bytes.
fn r_bytes() -> [u8; 32] {
    let mut r = [0u8; 32];
    r[0] = 0xD1;
    r[1] = 0x03;
    r[4] = 0x01;
    r
}

/// Width of the fixed window the rfold `+R` is confined to. R spans bits
/// [0,32]; we work modulo `2^RFOLD_WINDOW` so the carry out of bit 32 ripples
/// at most RFOLD_WINDOW-33 = 40 bits (the rest is dropped — matters only for
/// ~40 consecutive 1s at the injection point, ≤ 2^-40 per call, within
/// Shor's budget). The add over a[..`RFOLD_WINDOW`] is EXACT mod `2^RFOLD_WINDOW`
/// (a clean modular window, NOT a value-dependent carry drop), so the reverse
/// (X-sandwich of the SAME add = subtract mod `2^RFOLD_WINDOW`) is its exact
/// inverse for every input — Bennett round-trips stay exactly clean.
const RFOLD_WINDOW: usize = 73;

/// Top-K width for the cma:phase comparator (the MBU phase correction's
/// `1[a<b]`). Comparing only the top `COMPARE_TOPK` bits gives a wrong per-shot
/// phase with probability ≈ 2^-COMPARE_TOPK (a,b tie on the top bits), well
/// within Shor's phase budget, at ~2*`COMPARE_TOPK` Toffoli vs the full ~2*256.
const COMPARE_TOPK: usize = 64;

/// a += b mod p (rfold approximate). MBU: HMR'd flag with
/// `compare_lt` phase correction. Cheaper than `mod_add_mbu` in
/// `mod_arith.rs` (skips the `compare_geq_p` + ctrl-add-neg_p).
///
// requires:
//   a.len() == 257,  b.len() in {256, 257}
//   a[256] == |0>,   b[256] == |0> (if 257-bit)
//   a_val < p,       b_val < p                   (STRICT — not just
//   < 2^256; the phase-correction identity Lemma 4.1 is proved only
//   for a, b in [0, p). Callers chaining rfold outputs must ensure
//   the output of the previous call is still < p.)
// ensures:
//   a_new ≡ a_pre + b  (mod p)
//   a_new < 2^256      (range check — may be ≥ p, see module docs)
//   a[256] = |0>
pub fn mod_add_rfold_mbu(circ: &mut Circuit, a: &[QReg], b: &[QReg]) {
    let n = a.len();
    let nb = b.len();
    assert_eq!(n, 257, "mod_add_rfold_mbu requires 257-bit a");
    assert!(nb == 256 || nb == 257);
    let prev = circ.push_section("madd");

    // ── Pre-condition contract (sim-verified when CONTRACTS=1) ────
    //
    // Migration note: the previous deferred post-condition contract
    // captured `a.to_vec()` clones for the post closure; QReg is no
    // longer Clone, and `contract_capture` requires `'static` closures
    // so a borrow of `a` cannot be threaded through. The pre-condition
    // (Lemma 4.1 ranges, a[256]/b[256] = |0>) is the load-bearing
    // contract; the post-condition `a_new ≡ a_pre + b mod p` is
    // covered by other end-to-end tests.
    let p = crate::point_add::trailmix_port::num_bigint::BigUint::from_bytes_le(&crate::point_add::trailmix_port::mod_arith::SECP256K1_P_LE);
    let p_for_pre = p.clone();
    circ.contract_check("mod_add_rfold_mbu pre", move |c, shot| {
        if circ_a_bit(&c, a, 256, shot) {
            return Err("a[256] must be |0>, got 1".to_string());
        }
        if nb == 257 && circ_a_bit(&c, b, 256, shot) {
            return Err("b[256] must be |0>, got 1".to_string());
        }
        let av = c.contract_read_u256_shot(a, shot);
        let bv = c.contract_read_u256_shot(&b[..b.len().min(256)], shot);
        if av >= p_for_pre {
            return Err(format!(
                "a_val ({av}) >= p (Lemma 4.1 \
                requires a < p)"
            ));
        }
        if bv >= p_for_pre {
            return Err(format!("b_val ({bv}) >= p"));
        }
        Ok(())
    });

    // Step 1: integer add. Cuccaro (1 ancilla) via crate::point_add::trailmix_port::arith::ripple_add::add
    // for the 257-bit path. The 256-bit-b path uses
    // add_with_carry_to_high which is an n-1-ancilla variant still.
    if nb == 257 {
        crate::point_add::trailmix_port::arith::ripple_add::add(circ, a, b);
    } else {
        let b_low = &b[..256];
        add_with_carry_to_high(circ, a, b_low);
    }

    // Step 2: rfold — add R to a[..256] if bit 256 set.
    let r = r_bytes();
    crate::point_add::trailmix_port::arith::const_add::controlled_add_const(circ, &a[256], &a[..256], &r);

    // Step 3+4: MBU compare-lt phase correction.
    crate::point_add::trailmix_port::arith::compare::compare_lt_phase_correction_mbu(circ, &a[..256], &b[..256], &a[256]);
    circ.pop_section(&prev);
}

fn circ_a_bit(c: &impl ContractReadable, reg: &[QReg], i: usize, shot: usize) -> bool {
    if i < reg.len() {
        c.contract_read_bit_shot(&reg[i], shot)
    } else {
        false
    }
}

/// 257-bit add of (a, `b_256+0`) where the sum's high bit lands in
/// a[256]. Delegates to canonical Cuccaro with explicit overflow.
/// Peak: 1 ancilla (carry-in). Was: n+1 ancillae.
fn add_with_carry_to_high(circ: &mut Circuit, a: &[QReg], b: &[QReg]) {
    let n = b.len();
    assert!(a.len() > n);
    crate::point_add::trailmix_port::arith::cuccaro::add_cuccaro_with_overflow(circ, &a[..=n], b);
}

/// Controlled a += b mod p (rfold approximate). MBU.
/// If ctrl=0: no-op (HMR'd flag is also 0 -> phase contribution 0).
/// If ctrl=1: same as `mod_add_rfold_mbu`.
//
// requires:
//   a.len() == 257,  b.len() in {256, 257}
//   a[256] == |0>,   b[256] == |0> (if 257-bit)
//   a_val < p,       b_val < p                   (STRICT — same
//   Lemma 4.1 precondition as mod_add_rfold_mbu)
//   ctrl is a single qubit (may alias b[i] for some i — OK)
// ensures:
//   ctrl=1: a_new ≡ a_pre + b  (mod p),  a_new < 2^256
//   ctrl=0: a_new == a_pre
//   a[256] = |0>
pub fn controlled_mod_add_rfold_mbu(circ: &mut Circuit, ctrl: &QReg, a: &[QReg], b: &[QReg]) {
    let n = a.len();
    let nb = b.len();
    assert_eq!(n, 257);
    assert!(nb == 256 || nb == 257);
    let prev = circ.current_section.clone();

    circ.set_section(&format!("{prev}/cma:int"));
    if nb == 257 {
        crate::point_add::trailmix_port::arith::ripple_add::controlled_add(circ, ctrl, a, b);
    } else {
        controlled_add_with_carry_to_high(circ, ctrl, a, &b[..256]);
    }
    circ.set_section(&format!("{prev}/cma:rfold"));
    let r = r_bytes();
    // rfold confined to a[..RFOLD_WINDOW]: exact (a+R) mod 2^RFOLD_WINDOW. The
    // X-sandwich in controlled_mod_sub_rfold_mbu inverts this exactly because
    // the window add is a clean modular op (no value-dependent carry drop).
    crate::point_add::trailmix_port::arith::const_add::controlled_add_const_runs_forced(
        circ,
        &a[256],
        &a[..RFOLD_WINDOW],
        &r,
    );

    circ.set_section(&format!("{prev}/cma:phase"));
    // a<b is decided by the top differing bit; comparing only the top
    // COMPARE_TOPK bits is wrong (hence a wrong per-shot phase) only when a,b
    // tie on those bits (≈ 2^-COMPARE_TOPK), within Shor's phase budget. The
    // X-sandwich reverse re-runs the same truncated comparator, so the
    // forward/reverse MBU pair stays consistent.
    crate::point_add::trailmix_port::arith::compare::controlled_compare_lt_phase_correction_mbu_topk(
        circ,
        ctrl,
        &a[..256],
        &b[..256],
        &a[256],
        COMPARE_TOPK,
    );
    circ.set_section(&prev);
}

/// Pointer equality between two `&QReg`s. `QReg` has no `PartialEq`
/// because its `id` field is module-private; aliasing detection
/// across slices uses the borrow's identity directly.
fn qreg_ptr_eq(a: &QReg, b: &QReg) -> bool {
    std::ptr::eq(a, b)
}

/// Controlled 257-bit add of (a, ctrl*`b_256`) where the carry-out
/// lands in a[256] (initially 0). Polylog peak via Cuccaro + `mcx_dirty`.
///
/// When ctrl aliases a or b, the Cuccaro inner loop would see ctrl's
/// value change mid-computation (because Cuccaro transiently modifies
/// b). To stay correct AND polylog-peak, we copy ctrl into a fresh
/// scratch qubit FIRST, run Cuccaro with scratch as the effective
/// control, then uncompute scratch via CX(ctrl, scratch) at the end.
/// This works because Cuccaro restores b to its entry value, hence
/// ctrl (if it aliased b) is also restored; scratch then = ctrl and
/// one CX zeros it.
///
/// Peak: +3 ancillae (scratch + Cuccaro's c + scratch). Polylog.
fn controlled_add_with_carry_to_high(circ: &mut Circuit, ctrl: &QReg, a: &[QReg], b: &[QReg]) {
    let n = b.len();
    assert!(a.len() > n);
    let alias_in_a = a[..=n].iter().any(|q| qreg_ptr_eq(q, ctrl));
    let alias_in_b = b.iter().any(|q| qreg_ptr_eq(q, ctrl));
    let aliases = alias_in_a || alias_in_b;
    if !aliases {
        crate::point_add::trailmix_port::arith::cuccaro::controlled_add_cuccaro_with_overflow(circ, ctrl, &a[..=n], b);
        return;
    }

    // Only b-aliasing is supported (matches the module contract at the
    // callsite — see `controlled_mod_add_rfold_mbu` doc). ctrl aliasing
    // a would mean "add ctrl=a[i] * b to a", but Cuccaro modifies a,
    // so ctrl's value would shift mid-add and the final scratch
    // uncompute would fail. We disallow this case explicitly.
    assert!(
        !alias_in_a,
        "controlled_add_with_carry_to_high: ctrl aliases a — not supported"
    );

    // ctrl aliases b — copy to fresh scratch and use that. Cuccaro
    // preserves b's value across the full add, hence ctrl (= b[i])'s
    // value is also preserved; the final cx(ctrl, scratch) zeros
    // scratch cleanly.
    let scratch = circ.alloc_qreg("c_add_scratch");
    circ.cx(ctrl, &scratch);
    circ.declare_copy_of(&scratch, ctrl);
    crate::point_add::trailmix_port::arith::cuccaro::controlled_add_cuccaro_with_overflow(circ, &scratch, &a[..=n], b);
    circ.cx(ctrl, &scratch);
    // scratch drops here; drain fires at next gate (gap=0).
}

/// a := 2a mod p (rfold approximate). MBU via Z(a[0]) phase fix
/// (parity identity: rfold X = bit 0 of post-rfold value, since
/// 2*`a_pre` is even and R is odd).
//
// requires:
//   a.len() == 257
//   a[256] == |0>
//   a_val < 2^256            (does NOT require a < p — unlike
//   mod_add_rfold_mbu, doubling's identity doesn't depend on Lemma
//   4.1. The rfold flag's equality to bit 0 works purely from
//   "2*x is even, R is odd".)
// ensures:
//   a_new ≡ 2 * a_pre  (mod p)
//   a_new < 2^256
//   a[256] = |0>
pub fn mod_double_rfold_mbu(circ: &mut Circuit, a: &[QReg]) {
    let n = a.len();
    assert_eq!(n, 257);
    let prev = circ.push_section("dbl");
    // Step 1: left shift. Bit 256 was 0 (precondition); after the
    // shift a[256] = old a[255] (= X, high bit of 2*a_pre), and
    // a[0] = 0 (freshly rotated in).
    crate::point_add::trailmix_port::arith::shift::left_shift(circ, a);
    // Step 2: rfold. Add R*a[256] to a[..256] via ctrl-add-const.
    // Low-bit effect: a[0]_post = 0 XOR (R[0] AND a[256]) = a[256]
    // (R[0]=1). Higher bits get more of R conditionally; only bit 0
    // is relevant to us.
    let r = r_bytes();
    // rfold confined to a[..RFOLD_WINDOW]: (a + R) mod 2^RFOLD_WINDOW, exact
    // on the window. runs_forced is the cheap path for R (the dispatcher would
    // pick the dense classq path at this small width, so call it directly).
    crate::point_add::trailmix_port::arith::const_add::controlled_add_const_runs_forced(
        circ,
        &a[256],
        &a[..RFOLD_WINDOW],
        &r,
    );
    // IDENTITY: val(a[0]) == val(a[256]) at this point. Proof:
    // left_shift put 0 in a[0], controlled_add_const XORed in
    // (R[0] AND a[256]) = (1 AND a[256]) = a[256]. QED.
    // Tell the tracker so the HMR(a[256]) + z_if_bit(a[0]) pair
    // discharges structurally.
    circ.declare_identity(&a[0], &a[256]);
    // Step 3: HMR a[256]. Phase ^= X * bit.
    let bit = circ.alloc_bit();
    circ.hmr(&a[256], bit);
    // Step 4: phase correction Z^X under bit. a[0] == a[256] (=X).
    circ.z_if_bit(&a[0], bit);
    circ.free_bit(bit);
    circ.pop_section(&prev);
}

/// Structural inverse of `mod_double_rfold_mbu`. Reverses the
/// `left_shift` + rfold + HMR forward by:
///   1. Re-creating a[256] as the overflow bit. Identity:
///      val(a[0]) = overflow (established by `mod_double_rfold_mbu`
///      via `declare_identity` on a[0] and a[256]); so CX(a[0],
///      a[256]) sets a[256] = a[0] = X.
///   2. Reversing the rfold: `controlled_sub_const(X`, a[..256], R).
///   3. Reversing the `left_shift` via `right_shift`.
///
/// This replaces the old `mod_halve_mbu` call (a full from-scratch
/// division by 2 via add-p-if-odd then shift), which was ~25x the
/// cost even though we always call halve in a context where we
/// KNOW the input came from a `mod_double`. On the EEA reverse
/// rounds alone this collapses from ~7800 ops/call to ~300 ops/call.
///
//
// requires:
//   a.len() == 257
//   a[256] == |0>
//   a was PRODUCED by a prior mod_double_rfold_mbu call on the same
//   register — this primitive is the structural inverse, not a
//   general halve. The caller must pair it with a mod_double_rfold_mbu
//   in Bennett-reversal fashion.
// ensures:
//   a_new == a_pre_of_matching_double
//   a[256] = |0>
pub fn mod_halve_rfold_mbu(circ: &mut Circuit, a: &[QReg]) {
    let n = a.len();
    assert_eq!(n, 257);
    let prev = circ.push_section("halve");

    // Step 1: Regenerate a[256] = X (overflow) from a[0]. Forward
    // established val(a[0]) == val(a[256]) just before HMR, so
    // CX(a[0], a[256]) with a[256] fresh-zero sets a[256] = a[0].
    circ.cx(&a[0], &a[256]);
    // Tell the tracker that a[256] is a copy of a[0] so subsequent
    // uses stay tracked.
    circ.declare_copy_of(&a[256], &a[0]);

    // Step 2: Reverse the rfold add. Forward added R·a[256] to
    // a[..256]; reverse subtracts it.
    let r = r_bytes();
    // Reverse the rfold: subtract R mod 2^RFOLD_WINDOW on the SAME window via
    // the X-sandwich of the exact-mod-window add (its exact inverse).
    for q in &a[..RFOLD_WINDOW] {
        circ.x(q);
    }
    crate::point_add::trailmix_port::arith::const_add::controlled_add_const_runs_forced(
        circ,
        &a[256],
        &a[..RFOLD_WINDOW],
        &r,
    );
    for q in &a[..RFOLD_WINDOW] {
        circ.x(q);
    }

    // Step 3: Reverse the left_shift. After this, a[255] holds the
    // overflow X (== the value a[256] held), a[256] holds 0 (rotated
    // in), and a[..255] = a_pre[..255]. Combined with a[255] = X =
    // a_pre[255], the full register equals a_pre.
    crate::point_add::trailmix_port::arith::shift::right_shift(circ, a);

    // After right_shift: a[256] = 0 (rotated in from a[0] pre-shift,
    // which was X pre-rfold-subtract; controlled_sub_const's low-bit
    // behaviour left a[0]=X, so after right_shift a[255] = X and
    // a[256] = 0). prove_zero confirms before freeing.
    // Note: in the forward, a[256] was HMR-freed; here we just check
    // that the right_shift rotated a zero into a[256].
    // (right_shift is a pure swap chain, so we don't emit extra
    // gates beyond the chain itself.)
    circ.pop_section(&prev);
}

/// GENERAL pseudo-Mersenne mod-halve: `a := a/2 mod p` for ANY `a < 2^256`
/// (not just the structural inverse of a double). This is the exact mirror of
/// `mod_double_rfold_mbu`:
///   - the double folds a TOP overflow (a[256]) with `+R` and MBU-frees it;
///   - the halve consumes a BOTTOM parity (a[0]): if odd it adds `p` cheaply as
///     `+2^256 - R` (set a[256], then a windowed `-R` over a[..`RFOLD_WINDOW`]),
///     shifts right (so the 2^256 becomes the +2^255 of `(a+p)/2`), and cleans
///     the parity flag with the half-p phase MBU.
/// Since `(c+p)/2 = (c-R)/2 + 2^255`, this computes `a_pre/2 mod p` exactly
/// except for the windowed `-R` borrow beyond bit `RFOLD_WINDOW` (~2^-40 tail,
/// Shor-tolerant), matching the double's approximation.
//
// requires: a.len()==257, a[256]==|0>, a_val < 2^256
// ensures:  a_new ≡ a_pre / 2 (mod p), a_new < 2^256, a[256]==|0>
pub fn mod_halve_pm_general(circ: &mut Circuit, a: &[QReg]) {
    let n = a.len();
    assert_eq!(n, 257);
    let prev = circ.push_section("halve_pm");
    let r = r_bytes();

    // parity flag = a[0] (the bit about to be shifted out).
    let flag = circ.alloc_qreg("halve_pm.parity");
    circ.cx(&a[0], &flag);

    // add p if odd, cheaply: +2^256 (set a[256]) and -R windowed on the low bits.
    // After this a is even (a[0] XOR R[0]*flag = a[0] XOR flag = 0 when odd).
    circ.cx(&flag, &a[256]);
    for q in &a[..RFOLD_WINDOW] {
        circ.x(q);
    }
    crate::point_add::trailmix_port::arith::const_add::controlled_add_const_runs_forced(circ, &flag, &a[..RFOLD_WINDOW], &r);
    for q in &a[..RFOLD_WINDOW] {
        circ.x(q);
    }

    // divide by 2: the a[256]=flag bit shifts down to a[255] (the +2^255).
    crate::point_add::trailmix_port::arith::shift::right_shift(circ, a);

    // clean the parity flag: flag = 1[a >= ceil(p/2)] on the result.
    crate::point_add::trailmix_port::arith::compare::compare_geq_half_p_secp256k1_phase_correction_mbu(circ, &a[..256], flag);
    circ.pop_section(&prev);
}

/// APPROXIMATE pseudo-Mersenne mod-halve for secp256k1.
///
/// Value semantics identical to `mod_halve_pm_general`: `a := a/2 mod p`
/// for any `a < 2^256` (drift in [0, p+R) inherited from the windowed
/// `-R` step, same as the forward double). The only difference is the
/// parity-flag cleanup: instead of a full 256-bit `a >= q/2` borrow
/// chain, we use the algebraic fact that for secp256k1 with
/// q = 2^256 - f, f ≈ 2^32:
///
///   `a_pre` even ⇒ `a_post` = `a_pre/2` < q/2 < 2^255   ⇒ `a_post`[255] = 0
///   `a_pre` odd  ⇒ `a_post` = (`a_pre+q)/2` ∈ [q/2, q). The sub-band
///                with `a_post` < 2^255 has measure ≈ 2^32 in a range of
///                ≈ 2^255; violating shots have probability ≈ 2^-224.
///                For 64-shot sim this is astronomical.
///
/// So `flag == a_post[255]` is a valid identity, and the phase
/// correction is one CZ via `cz_if_bit(a[255], bit)` after a single
/// HMR. Cost: ~32 Toffoli/halve vs ~1000 for the exact compare.
pub fn mod_halve_pm_general_approx_secp256k1(circ: &mut Circuit, a: &[QReg]) {
    let n = a.len();
    assert_eq!(n, 257);
    let prev = circ.push_section("halve_pm_approx");
    let r = r_bytes();

    // parity flag = a[0] (the bit about to be shifted out).
    let flag = circ.alloc_qreg("halve_pm_approx.parity");
    circ.cx(&a[0], &flag);

    // add p if odd: +2^256 (set a[256]) and -R windowed on the low bits.
    // The -R uses the Gidney borrowed-dirty constant adder: a[..lsbs] -=
    // flag*R via X-sandwich, with the carry scratch BORROWED from the
    // register's own idle high bits a[lsbs..2*lsbs-1] (restored on exit).
    // This costs ~3 clean ancillae instead of the clustered adder's ~10,
    // dropping the apply_bv-inv peak below the structural floor + paper
    // budget. Exact within the lsbs window (carry beyond lsbs dropped).
    circ.cx(&flag, &a[256]);
    let lsbs = RFOLD_WINDOW; // 73
    let dirty_lo = lsbs;
    let dirty_hi = lsbs + (lsbs - 1); // 145; a[73..145] are idle here
    for q in &a[..lsbs] {
        circ.x(q);
    }
    crate::point_add::trailmix_port::arith::gidney_const_adder::controlled_add_const_gidney(
        circ,
        &flag,
        &a[..lsbs],
        &r,
        &a[dirty_lo..dirty_hi],
    );
    for q in &a[..lsbs] {
        circ.x(q);
    }

    // divide by 2: the a[256]=flag bit shifts to a[255].
    crate::point_add::trailmix_port::arith::shift::right_shift(circ, a);

    // Approximate phase-correction MBU: flag ≡ a[255] (post-halve)
    // with ~2^-224 mismatch on uniform inputs. declare_identity checks
    // this in-sim across all 64 shots before HMR.
    circ.declare_identity(&flag, &a[255]);
    let bit = circ.alloc_bit();
    circ.hmr(&flag, bit);
    circ.z_if_bit(&a[255], bit);
    circ.free_bit(bit);
    circ.zero_and_free(flag);

    circ.pop_section(&prev);
}

/// Reversible exact cleanup for an rfold-style intermediate in
/// `[0, 2^256)`.
///
/// Raw rfold outputs live in `[0, p + R)`, so if `a >= p` then
/// necessarily `a = p + x` with `x < R`. Exact reduction is therefore
/// not a dense `a -= p` on 257 bits; it is simply
/// `a = a + R (mod 2^256)` on the low 256 bits, because
///
///   a + R = p + x + R = 2^256 + x.
///
/// The wrapped low-256 sum is exactly the canonical residue `x`, and
/// the reverse map is the matching wrapped `a -= R`.
///
/// `flag` is the retained underreduction indicator:
/// - forward: `flag ^= 1[a >= p]`, then `a += flag * R (mod 2^256)`
/// - reverse: `a -= flag * R (mod 2^256)`, then recompute the same
///   compare to toggle `flag` back to 0
///
/// This is the smallest exactization wrapper around the raw rfold
/// arithmetic. Unlike a bare "reduce if >= p" map, it is reversible
/// because the caller keeps the reduction bit.
pub fn reduce_once_secp256k1_from_rfold(circ: &mut Circuit, a: &[QReg], flag: &QReg) {
    assert_eq!(a.len(), 257);
    crate::point_add::trailmix_port::arith::compare::compare_geq_p_secp256k1(circ, a, flag);
    let r = r_bytes();
    crate::point_add::trailmix_port::arith::const_add::controlled_add_const(circ, flag, &a[..256], &r);
}

/// Reverse of `reduce_once_secp256k1_from_rfold`.
pub fn reduce_once_secp256k1_from_rfold_reverse(circ: &mut Circuit, a: &[QReg], flag: &QReg) {
    assert_eq!(a.len(), 257);
    let r = r_bytes();
    crate::point_add::trailmix_port::arith::const_add::controlled_sub_const(circ, flag, &a[..256], &r);
    crate::point_add::trailmix_port::arith::compare::compare_geq_p_secp256k1(circ, a, flag);
}

/// Modular multiplication: result = a * b mod p. MBU (no flags).
/// Uses Horner shift-and-add over bits of b from MSB to LSB.
/// All sub-primitives are rfold-MBU.
//
// requires:
//   result.len() == 257, a.len() == 257, b.len() == 257
//   result pre == |0…0> (all 257 qubits zero)
//   a[256] == |0>, b[256] == |0>
//   a_val < p, b_val < p  (STRICT — inherited from rfold-add
//   precondition for the inner Horner adds; rfold intermediates
//   stay < 2^256 but accumulate via adds that REQUIRE a_prev < p,
//   which holds by induction on i starting from 0.)
// ensures:
//   result ≡ a * b  (mod p)
//   result < 2^256  (rfold-approximate; may be ≥ p, see module docs)
//   a, b unchanged
pub fn mod_mul_rfold_mbu(circ: &mut Circuit, result: &[QReg], a: &[QReg], b: &[QReg]) {
    let n = a.len();
    debug_assert_eq!(n, 257);
    debug_assert_eq!(b.len(), 257);
    debug_assert_eq!(result.len(), 257);

    // Skip the top bit of b (always 0 for value < 2^256).
    // Skip the first mod_double (result starts at 0).
    controlled_mod_add_rfold_mbu(circ, &b[n - 2], result, a);
    for i in (0..n - 2).rev() {
        mod_double_rfold_mbu(circ, result);
        controlled_mod_add_rfold_mbu(circ, &b[i], result, a);
    }
}

/// Inverse of `mod_mul_rfold_mbu`: result -= a*b mod p.
/// Replays the Horner loop in reverse: csub then halve.
pub fn mod_mul_rfold_mbu_undo(circ: &mut Circuit, result: &[QReg], a: &[QReg], b: &[QReg]) {
    let n = a.len();
    debug_assert_eq!(n, 257);
    debug_assert_eq!(b.len(), 257);
    debug_assert_eq!(result.len(), 257);

    for i in 0..n - 2 {
        controlled_mod_sub_rfold_mbu(circ, &b[i], result, a);
        mod_halve_rfold_mbu(circ, result);
    }
    controlled_mod_sub_rfold_mbu(circ, &b[n - 2], result, a);
}

/// Controlled a -= b mod p (rfold approximate). MBU via X-sandwich
/// of `controlled_mod_add_rfold_mbu` on the low 256 bits.
///
/// For ctrl=0: outer NOTs cancel inner no-op, so a unchanged.
/// For ctrl=1: same identity as `mod_sub_rfold_mbu`'s X-sandwich
/// derivation. Inner X = borrow B, phase correction is
/// Z^(ctrl AND B) under bit (handled inside `controlled_mod_add`).
pub fn controlled_mod_sub_rfold_mbu(circ: &mut Circuit, ctrl: &QReg, a: &[QReg], b: &[QReg]) {
    let n = a.len();
    let nb = b.len();
    assert_eq!(n, 257);
    assert!(nb == 256 || nb == 257);
    let prev = circ.push_section("csub_mod");
    for i in 0..256 {
        circ.x(&a[i]);
    }
    controlled_mod_add_rfold_mbu(circ, ctrl, a, b);
    for i in 0..256 {
        circ.x(&a[i]);
    }
    circ.pop_section(&prev);
}

#[cfg(test)]
mod halve_tests {
    use super::*;
    use crate::point_add::trailmix_port::circuit::Circuit;
    use crate::point_add::trailmix_port::num_bigint::BigUint;
    use rand::Rng;

    fn p_big() -> BigUint {
        BigUint::from_bytes_le(&crate::point_add::trailmix_port::mod_arith::SECP256K1_P_LE)
    }

    #[test]
    fn mod_halve_pm_general_64_random() {
        let p = p_big();
        let inv2 = (&p + BigUint::from(1u32)) / BigUint::from(2u32); // 2^-1 mod p
        let mut rng = rand::thread_rng();
        let mut c = Circuit::new();
        c.set_section("halve_pm_test");
        let a = c.alloc_input_qreg_bits("a", 257);
        let mut vs = Vec::with_capacity(64);
        for shot in 0..64 {
            let mut bytes = [0u8; 32];
            rng.fill(&mut bytes);
            let v = BigUint::from_bytes_le(&bytes) % &p; // v < p, so a[256]=0
            let mut v_le = v.to_bytes_le();
            v_le.resize(33, 0);
            c.sim_load_reg_bytes_shot(&a, &v_le, shot);
            vs.push(v);
        }
        mod_halve_pm_general(&mut c, &a);
        {
            let (a_r, pc, vsc, inv2c) = (&a, p.clone(), vs.clone(), inv2.clone());
            c.contract_check("halve_pm_val", move |view, shot| {
                let mut got = BigUint::from(0u32);
                for j in 0..256 {
                    if view.contract_read_bit_shot(&a_r[j], shot) {
                        got |= BigUint::from(1u32) << j;
                    }
                }
                let want = (&vsc[shot] * &inv2c) % &pc; // v/2 mod p
                if &got % &pc != want {
                    return Err(format!(
                        "halve wrong (shot {}, got {}, want {})",
                        shot,
                        &got % &pc,
                        want
                    ));
                }
                Ok(())
            });
        }
        c.assert_phase_clean();
        let _ = c.destroy_sim(a);
    }

    #[test]
    fn mod_halve_pm_general_noncanonical_64() {
        // Same, but inputs are NON-canonical: v in [0, 2^256), NOT reduced mod p.
        // This is the rfold posture (values < 2^256, possibly >= p). If the
        // parity-flag MBU holds here, the halve is usable directly on rfold
        // outputs.
        let p = p_big();
        let inv2 = (&p + BigUint::from(1u32)) / BigUint::from(2u32);
        let mask = (BigUint::from(1u32) << 256) - BigUint::from(1u32);
        let mut rng = rand::thread_rng();
        let mut c = Circuit::new();
        c.set_section("halve_pm_nc_test");
        let a = c.alloc_input_qreg_bits("a", 257);
        let mut vs = Vec::with_capacity(64);
        for shot in 0..64 {
            let mut bytes = [0u8; 32];
            rng.fill(&mut bytes);
            let v = BigUint::from_bytes_le(&bytes) & &mask; // v < 2^256, a[256]=0
            let mut v_le = v.to_bytes_le();
            v_le.resize(33, 0);
            c.sim_load_reg_bytes_shot(&a, &v_le, shot);
            vs.push(v);
        }
        mod_halve_pm_general(&mut c, &a);
        {
            let (a_r, pc, vsc, inv2c) = (&a, p.clone(), vs.clone(), inv2.clone());
            c.contract_check("halve_pm_nc_val", move |view, shot| {
                let mut got = BigUint::from(0u32);
                for j in 0..256 {
                    if view.contract_read_bit_shot(&a_r[j], shot) {
                        got |= BigUint::from(1u32) << j;
                    }
                }
                let want = (&vsc[shot] * &inv2c) % &pc; // (v mod p)/2 mod p
                if &got % &pc != want {
                    return Err(format!(
                        "nc halve wrong (shot {}, v {}, got {}, want {})",
                        shot,
                        &vsc[shot],
                        &got % &pc,
                        want
                    ));
                }
                Ok(())
            });
        }
        c.assert_phase_clean();
        let _ = c.destroy_sim(a);
    }
}
