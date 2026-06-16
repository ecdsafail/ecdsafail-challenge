//! Pseudo-Mersenne approximate modular arithmetic per Schrottenloher
//! 2026 Sec. 4 (Alg 5-11).
//!
//! For q = 2^u - f with f << 2^u (secp256k1: u=256, f = 2^32+977),
//! mod-double and controlled mod-add can be implemented with one
//! n-bit shift / add and one ~bitlen(f)+padding-bit constant addition,
//! avoiding the full n-bit p-reduction.
use crate::point_add::trailmix_port::arith::schrottenloher::msb_compare;
use crate::point_add::trailmix_port::circuit::{Circuit, QReg};

/// Measurement-vent ancilla budget for the square's per-term controlled
/// mod-adds (Gidney hybrid register adder). The square runs after the
/// inversion frees its workspace, so it has ample headroom: `usize::MAX`
/// clamps to `n` (full Gidney 2n). Tuned down only if the square's vents would
/// lift the overall EC-add peak. (Matches the author resetting to the `GidneyAdder` backend in
/// `ControlledEfficientModularSquareAdd` "since we have lots of ancilla".)
const SQR_VENTS: usize = usize::MAX;

/// `reg[..lsbs] += ctrl * f` via the Gidney borrowed-dirty constant
/// adder, borrowing `reg[lsbs..2*lsbs-1]` (the register's own idle high
/// bits) as carry scratch — ~3 clean ancillae instead of the clustered
/// adder's ~10, so the apply-bitvector peak stays near the structural
/// floor. Exact within the `lsbs`-bit window (carry beyond dropped).
/// `ctrl` must lie outside `[..2*lsbs-1]` (e.g. the n=256 anc bit).
///
/// `f_vents > 0` selects the CLEAN-ancilla measurement-vented path: f is
/// materialized into a fresh `lsbs`-bit window and added with the Gidney
/// hybrid 2n adder (`3n-2-vents` Toffoli) instead of the borrowed-dirty 3n.
/// That drops the `XORCarries` restore pass (~n Toffoli/call) at a cost of
/// `+2*lsbs` peak qubits, so it is only used where headroom exists (the
/// square, which runs after the inversion frees its workspace). Value- and
/// phase-identical to the borrowed-dirty path (same mod-2^lsbs window).
fn gidney_cadd_f_window(
    circ: &mut Circuit,
    ctrl: &QReg,
    reg: &[QReg],
    lsbs: usize,
    f: &[u8],
    f_vents: usize,
) {
    if f_vents == 0 {
        let dh = lsbs + (lsbs - 1);
        debug_assert!(dh <= reg.len(), "register too short to borrow dirty bits");
        crate::point_add::trailmix_port::arith::gidney_const_adder::controlled_add_const_gidney(
            circ,
            ctrl,
            &reg[..lsbs],
            f,
            &reg[lsbs..dh],
        );
        return;
    }
    // Clean-vent path: materialize f, hybrid-add, unmaterialize, free.
    let fbit = |i: usize| -> bool { i / 8 < f.len() && (f[i / 8] >> (i % 8)) & 1 == 1 };
    let f_reg: Vec<QReg> = (0..lsbs).map(|_| circ.alloc_qreg("fmat")).collect();
    for i in 0..lsbs {
        if fbit(i) {
            circ.x(&f_reg[i]);
        }
    }
    crate::point_add::trailmix_port::arith::gidney_const_adder::controlled_hybrid_add(
        circ,
        ctrl,
        &reg[..lsbs],
        &f_reg,
        f_vents,
    );
    for i in 0..lsbs {
        if fbit(i) {
            circ.x(&f_reg[i]);
        }
    }
    for q in f_reg {
        circ.zero_and_free(q);
    }
}

/// `reg[..lsbs] -= ctrl * f` (X-sandwich of [`gidney_cadd_f_window`]).
fn gidney_csub_f_window(
    circ: &mut Circuit,
    ctrl: &QReg,
    reg: &[QReg],
    lsbs: usize,
    f: &[u8],
    f_vents: usize,
) {
    for q in &reg[..lsbs] {
        circ.x(q);
    }
    gidney_cadd_f_window(circ, ctrl, reg, lsbs, f, f_vents);
    for q in &reg[..lsbs] {
        circ.x(q);
    }
}

/// secp256k1 controlled mod-add venting ONLY the (n+1)-bit register add
/// (`f_vents=0` keeps the ~63-bit +f reduction on the peak-safe borrowed-dirty
/// path, avoiding its +63 constant materialization). This is the
/// qubit-efficient apply_bv-peak vent: +`vents` peak, -`vents` Toffoli/call.
pub fn controlled_mod_add_pm_secp256k1_regvents(
    circ: &mut Circuit,
    ctrl: &QReg,
    x: &[QReg],
    y: &[QReg],
    vents: usize,
) {
    let f_bytes = super::F_SECP256K1.to_le_bytes();
    let padding = 30usize;
    let f_bitlen = 33usize;
    controlled_mod_add_pm(
        circ,
        ctrl,
        x,
        y,
        &f_bytes,
        padding + f_bitlen,
        padding,
        vents,
        0,
    );
}

/// secp256k1 controlled mod-sub venting ONLY the register sub (see
/// [`controlled_mod_add_pm_secp256k1_regvents`]).
pub fn controlled_mod_sub_pm_secp256k1_regvents(
    circ: &mut Circuit,
    ctrl: &QReg,
    x: &[QReg],
    y: &[QReg],
    vents: usize,
) {
    let f_bytes = super::F_SECP256K1.to_le_bytes();
    let padding = 30usize;
    let f_bitlen = 33usize;
    controlled_mod_sub_pm_v(
        circ,
        ctrl,
        x,
        y,
        &f_bytes,
        padding + f_bitlen,
        padding,
        vents,
        0,
    );
}

/// In-place pseudo-Mersenne modular doubling: a := 2*a mod q.
///
/// Matches Schrottenloher 2026 Algorithm 7 / Schrottenloher's
/// `SpecialPrimeModularDouble`:
///
/// ```text
///     left_shift(a + anc1)     # a := 2*a, anc1 = carry-out
///     cadd(anc1, f, a[0..lsbs])
///     cx(a[0], anc1)           # clean anc1
/// ```
///
/// Pre: `a.len()` == n + 1 (one extra slot for the overflow / control
/// ancilla); a[0..n] holds x in [0, q); a[n] = |0>.
/// Post (success): a[0..n] holds 2*x mod q (with x < q/2 + 2^(n-1)/q
/// approximation, see below); a[n] = |0>.
///
/// `lsbs` is the width of the controlled f-addition slice. The paper's
/// default `padding = 30` gives `lsbs = bitlen(f) + 30`. Per-call
/// failure probability ~ 2^(-padding). For secp256k1, f has bitlen 33,
/// so `lsbs = 63` at padding=30.
///
/// "Failure" here means: the slice's carry-out at bit `lsbs-1` would
/// have been needed to produce the correct integer answer, but is
/// dropped. Whenever the slice + f does NOT exceed 2^lsbs, the output
/// is exact.
///
/// The reduction is also approximate in the "2x in [q, 2^n)" sense
/// inherent to Alg 7: when 2*x < 2^n, no f-add is applied even if
/// 2*x >= q. The output stays in [0, 2^n). On uniform x, the residual
/// representation drift is bounded by f/q ~ 2^-224 per call and is
/// dominated by the slice-truncation failure.
pub fn mod_double_pm(circ: &mut Circuit, a: &[QReg], f_bytes: &[u8], lsbs: usize, f_vents: usize) {
    let n = a.len() - 1;
    assert!(lsbs <= n, "lsbs must fit in the data slice");

    // 1) Multiply by 2 -- value bit i moves to slot i+1 in a; a[n]
    //    captures the old MSB (overflow flag).
    crate::point_add::trailmix_port::arith::shift::left_shift(circ, a);

    // 2) If anc (= a[n]) is set, add f to bottom lsbs bits of a. Carry
    //    out of bit lsbs-1 is dropped (paper's approximation). The
    //    borrowed-dirty Gidney constant adder (f_vents=0) borrows a's idle
    //    high bits; with f_vents>0 (square only) f is materialized and added
    //    with the cheaper measure-vented hybrid adder.
    gidney_cadd_f_window(circ, &a[n], a, lsbs, f_bytes, f_vents);

    // 3) Clean the ancilla. After step 1 we had a[0] = 0. After
    //    step 2 with f odd (bit 0 = 1), a[0] = a[n]. So CX clears a[n].
    circ.cx(&a[0], &a[n]);
}

/// Inverse of `mod_double_pm`. Net: a holds (a / 2) mod q where the
/// "/2" interprets a as a 2*x-representation. Equivalent to applying
/// the operations in reverse.
pub fn mod_double_pm_reverse(
    circ: &mut Circuit,
    a: &[QReg],
    f_bytes: &[u8],
    lsbs: usize,
    f_vents: usize,
) {
    let n = a.len() - 1;
    assert!(lsbs <= n, "lsbs must fit in the data slice");

    // Inverse of step 3: re-set the ancilla from a[0]. Since on the
    // forward pass a[0] = a[n] post-cleanup, this restores a[n] to
    // its value just before the cx.
    circ.cx(&a[0], &a[n]);

    // Inverse of step 2: controlled SUBTRACT f, matching the forward
    // window/vents so the Bennett round-trip is exact.
    gidney_csub_f_window(circ, &a[n], a, lsbs, f_bytes, f_vents);

    // Inverse of step 1: divide by 2 (right shift).
    crate::point_add::trailmix_port::arith::shift::right_shift(circ, a);
}

/// Convenience wrapper specialized to secp256k1 (q = 2^256 - `F_SECP256K1`).
/// Borrowed-dirty +f (no extra ancillae) -- for the `apply_bv` peak phase.
pub fn mod_double_pm_secp256k1(circ: &mut Circuit, a: &[QReg]) {
    mod_double_pm_secp256k1_vents(circ, a, 0);
}

/// Like [`mod_double_pm_secp256k1`] but vents the ~lsbs-bit +f reduction
/// (`f_vents` measurement-vent ancillae, each -1 Toffoli / +1 peak). For
/// the higher-qubit EC-add config's `apply_bv` peak phase.
pub fn mod_double_pm_secp256k1_vents(circ: &mut Circuit, a: &[QReg], f_vents: usize) {
    let f_bytes = super::F_SECP256K1.to_le_bytes();
    let padding = 30usize; // matches paper default
    let f_bitlen = 33usize;
    mod_double_pm(circ, a, &f_bytes, padding + f_bitlen, f_vents);
}

/// Convenience reverse for secp256k1 (borrowed-dirty +f).
pub fn mod_double_pm_secp256k1_reverse(circ: &mut Circuit, a: &[QReg]) {
    let f_bytes = super::F_SECP256K1.to_le_bytes();
    let padding = 30usize;
    let f_bitlen = 33usize;
    mod_double_pm_reverse(circ, a, &f_bytes, padding + f_bitlen, 0);
}

/// Measure-vented +f variants for the square (which has ancilla headroom):
/// the +f window uses the Gidney hybrid 2n adder instead of borrowed-dirty 3n.
pub fn mod_double_pm_secp256k1_fvented(circ: &mut Circuit, a: &[QReg]) {
    let f_bytes = super::F_SECP256K1.to_le_bytes();
    mod_double_pm(circ, a, &f_bytes, 30 + 33, usize::MAX);
}

/// Reverse of [`mod_double_pm_secp256k1_fvented`].
pub fn mod_double_pm_secp256k1_reverse_fvented(circ: &mut Circuit, a: &[QReg]) {
    let f_bytes = super::F_SECP256K1.to_le_bytes();
    mod_double_pm_reverse(circ, a, &f_bytes, 30 + 33, usize::MAX);
}

/// In-place controlled pseudo-Mersenne modular addition (Alg 10).
///
/// `y := y + ctrl * x  (mod q)`. Does NOT handle the `x + y = q` case
/// (which has probability ~1/q on random inputs).
///
/// Conventions:
///   - `x`, `y` are both (n+1)-bit; `x[n] = 0`, `y[n] = 0` pre-call;
///     `y[n]` is used as the carry slot (`anc_y`) and returned to 0.
///     `x[n]` stays 0 (`anc_x` is not touched).
///   - `f_bytes`: little-endian encoding of `f = 2^n - q`.
///   - `lsbs`: width of the conditional f-add (default `bitlen(f) +
///     padding`).
///   - `msbs`: width of the approximate top-MSB comparator used to
///     reset `y[n]`. Failure rate ~2^-msbs per call.
///
/// Cost shape: one (n+1)-bit Cuccaro controlled add (3*(n+1)
/// Toffoli) + one ~lsbs-bit controlled constant add + one 4*msbs+1
/// Toffoli geq-msbs comparator.
pub fn controlled_mod_add_pm(
    circ: &mut Circuit,
    ctrl: &QReg,
    x: &[QReg],
    y: &[QReg],
    f_bytes: &[u8],
    lsbs: usize,
    msbs: usize,
    vents: usize,
    f_vents: usize,
) {
    let n = x.len() - 1;
    assert_eq!(y.len(), n + 1, "x and y must have matching width n+1");
    assert!(lsbs <= n, "lsbs must fit in the data slice");
    assert!(msbs <= n, "msbs must fit in the data slice");

    // 1) y += ctrl * x  (n+1 bit add). y[n] (anc_y) captures overflow.
    //    x[n] (anc_x) is read as 0 and restored. `vents` measurement-vent
    //    ancillae turn carry-uncompute Toffolis into measurements (Gidney
    //    hybrid adder); 0 -> Cuccaro 3n, n -> Gidney 2n.
    crate::point_add::trailmix_port::arith::gidney_const_adder::controlled_hybrid_add(circ, ctrl, y, x, vents);

    // 2) If anc_y is set, add f to bottom lsbs of y. Carry out of bit
    //    lsbs-1 is dropped (paper's approximation). `f_vents>0` uses the
    //    measure-vented +f (square headroom); 0 = borrowed-dirty (peak).
    gidney_cadd_f_window(circ, &y[n], y, lsbs, f_bytes, f_vents);

    // 3) Approximate clt: if ctrl AND y_msbs < x_msbs, XOR anc_y.
    //    Erases anc_y in the overflow case (see Alg 10 case analysis).
    //
    //    Implementation: ctrl AND (y_top < x_top) = ctrl AND (x_top >= y_top)
    //    modulo the top-msbs equality tail (~2^-msbs).
    msb_compare::controlled_lt_msbs(circ, ctrl, &y[..n], &x[..n], msbs, &y[n]);
}

/// Inverse of `controlled_mod_add_pm`. Net: y := y - ctrl * x (mod q).
pub fn controlled_mod_add_pm_reverse(
    circ: &mut Circuit,
    ctrl: &QReg,
    x: &[QReg],
    y: &[QReg],
    f_bytes: &[u8],
    lsbs: usize,
    msbs: usize,
    vents: usize,
    f_vents: usize,
) {
    let n = x.len() - 1;
    assert_eq!(y.len(), n + 1);

    // Reverse of step 3: same operation (self-inverse on target via
    // double XOR — but here it RECREATES anc_y from a clean 0).
    //
    // On the forward pass, post-step-3 we had anc_y = 0. The reverse
    // recomputes the same predicate to restore anc_y to whatever value
    // it had pre-step-3 (= the overflow flag).
    msb_compare::controlled_lt_msbs(circ, ctrl, &y[..n], &x[..n], msbs, &y[n]);

    // Reverse of step 2: controlled SUBTRACT f, matching the forward
    // window/vents for an exact Bennett round-trip.
    gidney_csub_f_window(circ, &y[n], y, lsbs, f_bytes, f_vents);

    // Reverse of step 1: controlled SUBTRACT x from y. Our codebase has
    // `controlled_add` but no direct controlled_sub for q-q; emulate
    // via X-sandwich on y (twos-complement subtraction).
    controlled_sub_qq(circ, ctrl, y, x, vents);
}

/// Helper: `y -= ctrl * x` for (n+1)-bit q-q registers. Uses the
/// X-sandwich trick: `~y := ~y + ctrl*x` ⟹ `y := y - ctrl*x`. `vents`
/// sizes the hybrid adder's measurement-vent ancillae.
fn controlled_sub_qq(circ: &mut Circuit, ctrl: &QReg, y: &[QReg], x: &[QReg], vents: usize) {
    for q in y {
        circ.x(q);
    }
    crate::point_add::trailmix_port::arith::gidney_const_adder::controlled_hybrid_add(circ, ctrl, y, x, vents);
    for q in y {
        circ.x(q);
    }
}

/// secp256k1 convenience wrapper for Alg 10 (does NOT handle x+y=q).
/// Space-efficient path (vents=0 -> Cuccaro 3n register add), for the
/// `apply_bv` reconstruction phase where the peak leaves no free ancillae.
pub fn controlled_mod_add_pm_secp256k1(circ: &mut Circuit, ctrl: &QReg, x: &[QReg], y: &[QReg]) {
    controlled_mod_add_pm_secp256k1_vents(circ, ctrl, x, y, 0);
}

/// Like [`controlled_mod_add_pm_secp256k1`] but with a measurement-vent
/// ancilla budget for the register add (Gidney hybrid). Used by the square,
/// which runs after inversion with ample free ancillae.
pub fn controlled_mod_add_pm_secp256k1_vents(
    circ: &mut Circuit,
    ctrl: &QReg,
    x: &[QReg],
    y: &[QReg],
    vents: usize,
) {
    let f_bytes = super::F_SECP256K1.to_le_bytes();
    let padding = 30usize;
    let f_bitlen = 33usize;
    // Couple the +f vents to the register vents: the square (vents=MAX) gets
    // the measure-vented +f, apply_bv (vents=0) keeps borrowed-dirty.
    controlled_mod_add_pm(
        circ,
        ctrl,
        x,
        y,
        &f_bytes,
        padding + f_bitlen,
        padding,
        vents,
        vents,
    );
}

/// In-place controlled pseudo-Mersenne modular SUBTRACTION (Alg 11).
///
/// `y := y - ctrl * x  (mod q)`. Approximate (same drift bound as
/// `controlled_mod_add_pm`). Cost mirrors the forward mod-add: one
/// (n+1)-bit Cuccaro controlled-sub (X-sandwich of `controlled_add`),
/// one ~lsbs-bit controlled-sub-const, one msbs-comparator cleanup.
///
/// Conventions match `controlled_mod_add_pm`: y[n] is the working
/// borrow-flag slot (= 0 pre/post). x[n] stays 0.
pub fn controlled_mod_sub_pm(
    circ: &mut Circuit,
    ctrl: &QReg,
    x: &[QReg],
    y: &[QReg],
    f_bytes: &[u8],
    lsbs: usize,
    msbs: usize,
    f_vents: usize,
) {
    // Default: peak-safe register sub (vents = 0).
    controlled_mod_sub_pm_v(circ, ctrl, x, y, f_bytes, lsbs, msbs, 0, f_vents);
}

/// Like [`controlled_mod_sub_pm`] but with a separate measurement-vent
/// budget `vents` for the (n+1)-bit register sub (the `f_vents` budget
/// covers the ~lsbs-bit +f reduction).
#[allow(clippy::too_many_arguments)]
pub fn controlled_mod_sub_pm_v(
    circ: &mut Circuit,
    ctrl: &QReg,
    x: &[QReg],
    y: &[QReg],
    f_bytes: &[u8],
    lsbs: usize,
    msbs: usize,
    vents: usize,
    f_vents: usize,
) {
    let n = x.len() - 1;
    assert_eq!(y.len(), n + 1, "x and y must have matching width n+1");
    assert!(lsbs <= n);
    assert!(msbs <= n);

    // Step 1: y -= ctrl * x in n+1 bits via X-sandwich of controlled_add.
    // Pre: y[n] = 0. Post: y[n] = 1 iff borrow (y_pre < ctrl*x). Route the
    // register op through controlled_hybrid_add (vents measurement-vent
    // ancillae; 0 keeps the Cuccaro 3n peak-safe path) so it opens the same
    // `hybrid_cadd` section as the forward mod-add; otherwise this subtract
    // falls through to the enclosing `apply_bv_inv` leaf and makes the
    // forward/inverse profile look 4x asymmetric when the two directions are
    // actually symmetric.
    for q in y {
        circ.x(q);
    }
    crate::point_add::trailmix_port::arith::gidney_const_adder::controlled_hybrid_add(circ, ctrl, y, x, vents);
    for q in y {
        circ.x(q);
    }

    // Step 2: if borrow (y[n]=1), subtract f from y[..lsbs]. This is
    // the q-correction: adding q = 2^n - f when borrow means subtracting
    // 2^n (already captured by clearing the borrow flag in step 3) and
    // subtracting f from the low bits. Carry/borrow beyond lsbs is
    // dropped (paper's approximation).
    gidney_csub_f_window(circ, &y[n], y, lsbs, f_bytes, f_vents);

    // Step 3: clean y[n]. Identity: y[n] = ctrl AND 1[borrow].
    // borrow = 1[y_pre < ctrl*x] = 1[y_post + ctrl*x >= q] (since
    // y_post = y_pre - ctrl*x + q*borrow). Approximate via top-msbs
    // add-overflow: 1[y_post_top + x_top overflows top K bits].
    // This is the structural mirror of forward mod-add's
    // controlled_lt_msbs cleanup — top-K add-carry instead of
    // top-K sub-borrow.
    crate::point_add::trailmix_port::arith::ripple_add::controlled_add_overflow_msbs_phase_correction_mbu(
        circ,
        ctrl,
        &y[..n],
        &x[..n],
        &y[n],
        msbs,
    );
}

/// secp256k1 convenience wrapper for Alg 11.
pub fn controlled_mod_sub_pm_secp256k1(circ: &mut Circuit, ctrl: &QReg, x: &[QReg], y: &[QReg]) {
    controlled_mod_sub_pm_secp256k1_vents(circ, ctrl, x, y, 0);
}

/// Like [`controlled_mod_sub_pm_secp256k1`] but with a measurement-vent
/// budget for BOTH the (n+1)-bit register sub and the ~lsbs-bit +f
/// reduction. `vents` ancillae turn carry-uncompute Toffolis into
/// measurements (each -1 Toffoli, +1 peak). Used by the higher-qubit
/// EC-add config to spend a shared ancilla pool at the `apply_bv` peak.
pub fn controlled_mod_sub_pm_secp256k1_vents(
    circ: &mut Circuit,
    ctrl: &QReg,
    x: &[QReg],
    y: &[QReg],
    vents: usize,
) {
    let f_bytes = super::F_SECP256K1.to_le_bytes();
    let padding = 30usize;
    let f_bitlen = 33usize;
    // Couple the +f vents to the register vents (each clamps to its own
    // width-1 internally): vents=0 keeps the peak-safe borrowed-dirty path.
    controlled_mod_sub_pm_v(
        circ,
        ctrl,
        x,
        y,
        &f_bytes,
        padding + f_bitlen,
        padding,
        vents,
        vents,
    );
}

/// secp256k1 reverse for Alg 10 (space-efficient, vents=0).
pub fn controlled_mod_add_pm_secp256k1_reverse(
    circ: &mut Circuit,
    ctrl: &QReg,
    x: &[QReg],
    y: &[QReg],
) {
    controlled_mod_add_pm_secp256k1_reverse_vents(circ, ctrl, x, y, 0);
}

/// Like [`controlled_mod_add_pm_secp256k1_reverse`] but with a
/// measurement-vent ancilla budget (Gidney hybrid). Used by the square's
/// reverse direction.
pub fn controlled_mod_add_pm_secp256k1_reverse_vents(
    circ: &mut Circuit,
    ctrl: &QReg,
    x: &[QReg],
    y: &[QReg],
    vents: usize,
) {
    let f_bytes = super::F_SECP256K1.to_le_bytes();
    let padding = 30usize;
    let f_bitlen = 33usize;
    // Couple +f vents to register vents (square -> vented +f).
    controlled_mod_add_pm_reverse(
        circ,
        ctrl,
        x,
        y,
        &f_bytes,
        padding + f_bitlen,
        padding,
        vents,
        vents,
    );
}

/// Controlled mod-square-add for Alg 1 line 11's "in-place mul" style
/// use, OR controlled mod-square-sub for line 10. Implements Schrottenloher's
/// `ControlledSpecialPrimeModularSquareAdd`
/// (`special_mod_arithmetic.py:274`). Horner-MSB-first composing
/// `mod_double_pm` and `controlled_mod_add_pm`:
///
///   if `direction = Add`:    `output_reg` += ctrl * `y_reg^2` (mod q)
///   if `direction = Sub`:    `output_reg` -= ctrl * `y_reg^2` (mod q)
///
/// Conventions:
///   - `ctrl`: 1-qubit. Active-high.
///   - `y_reg`: 257-bit (n+1, bit n = anc = 0 pre/post).
///   - `output_reg`: 257-bit (n+1, bit n = anc = 0 pre/post). Modified
///     in place. Caller's pre-value gets the +/-y^2 contribution.
///
/// Cost shape: 1 ccx + 1 `controlled_mod_add_pm` + 1 ccx + 1
/// `mod_double_pm` per bit of y, plus n-1 setup doublings. Approximately
/// 2 * (n * mod_mul_pm-iter-cost) ≈ 2x `mod_mul` cost.
pub fn controlled_mod_square_add_pm_secp256k1(
    circ: &mut Circuit,
    ctrl: &QReg,
    y_reg: &[QReg],
    output_reg: &[QReg],
) {
    let n = 256usize;
    assert_eq!(y_reg.len(), n + 1);
    assert_eq!(output_reg.len(), n + 1);
    let prev = circ.push_section("sqr_add_pm");
    let c = circ.alloc_qreg("sqr.c");
    for _ in 0..(n - 1) {
        mod_double_pm_secp256k1_reverse_fvented(circ, output_reg);
    }
    for i in 0..n {
        circ.ccx(ctrl, &y_reg[n - 1 - i], &c);
        controlled_mod_add_pm_secp256k1_vents(circ, &c, y_reg, output_reg, SQR_VENTS);
        circ.ccx(ctrl, &y_reg[n - 1 - i], &c);
        if i < n - 1 {
            mod_double_pm_secp256k1_fvented(circ, output_reg);
        }
    }
    circ.zero_and_free(c);
    circ.pop_section(&prev);
}

/// In-place modular constant addition: `dst += val mod q`.
/// `dst` is 257 bits (with anc at bit 256, 0 pre/post). `val` is
/// classical, fits in 256 bits. Uses `controlled_mod_add_pm` with a
/// classically-loaded temp register; simple but not space-optimized.
pub fn mod_add_const_pm_secp256k1(circ: &mut Circuit, dst: &[QReg], val: &crate::point_add::trailmix_port::compat::num_bigint::BigUint) {
    let n = 256usize;
    assert_eq!(dst.len(), n + 1);
    let prev = circ.push_section("mod_add_const_pm");
    let temp = circ.alloc_qreg_bits("const_src", n + 1);
    let val_bytes = val.to_bytes_le();
    for i in 0..n {
        let byte_idx = i / 8;
        if byte_idx < val_bytes.len() && ((val_bytes[byte_idx] >> (i % 8)) & 1) == 1 {
            circ.x(&temp[i]);
        }
    }
    let unit = circ.alloc_qreg("const_unit");
    circ.x(&unit);
    controlled_mod_add_pm_secp256k1(circ, &unit, &temp, dst);
    circ.x(&unit);
    circ.zero_and_free(unit);
    for i in 0..n {
        let byte_idx = i / 8;
        if byte_idx < val_bytes.len() && ((val_bytes[byte_idx] >> (i % 8)) & 1) == 1 {
            circ.x(&temp[i]);
        }
    }
    for q in temp {
        circ.zero_and_free(q);
    }
    circ.pop_section(&prev);
}

/// In-place modular constant subtraction: `dst -= val mod q`.
/// Equivalent to `mod_add_const_pm_secp256k1(dst, q - val)`.
pub fn mod_sub_const_pm_secp256k1(circ: &mut Circuit, dst: &[QReg], val: &crate::point_add::trailmix_port::compat::num_bigint::BigUint) {
    let q: crate::point_add::trailmix_port::compat::num_bigint::BigUint =
        (crate::point_add::trailmix_port::compat::num_bigint::BigUint::from(1u32) << 256u32) - crate::point_add::trailmix_port::compat::num_bigint::BigUint::from(super::F_SECP256K1);
    let zero = crate::point_add::trailmix_port::compat::num_bigint::BigUint::from(0u32);
    let neg = if val == &zero {
        zero
    } else {
        (q.clone() - (val % &q)) % &q
    };
    mod_add_const_pm_secp256k1(circ, dst, &neg);
}

/// UNCONDITIONAL in-place pseudo-Mersenne modular add: `y += x mod q`.
///
/// The ctrl-free form of [`controlled_mod_add_pm_secp256k1`] (no `|1>`
/// control dummy). `x`, `y` are 257-bit (n+1); bit n is the overflow
/// slot (0 pre/post). Approximate, same drift as the controlled form.
/// Uses the unconditional Cuccaro register add (`crate::point_add::trailmix_port::arith::ripple_add::add`),
/// borrowed-dirty `+f` window, and uncontrolled top-msbs `lt` cleanup.
pub fn mod_add_pm_secp256k1(circ: &mut Circuit, x: &[QReg], y: &[QReg]) {
    let n = 256usize;
    assert_eq!(x.len(), n + 1);
    assert_eq!(y.len(), n + 1);
    let prev = circ.push_section("mod_add_pm_uncond");
    let f_bytes = super::F_SECP256K1.to_le_bytes();
    let padding = 30usize;
    let f_bitlen = 33usize;
    let lsbs = padding + f_bitlen; // 63
    let msbs = padding; // 30
                        // 1) y += x (n+1-bit Cuccaro). y[n] captures the overflow carry.
    crate::point_add::trailmix_port::arith::ripple_add::add(circ, y, x);
    // 2) if y[n]: add f to the bottom lsbs bits (carry beyond dropped).
    gidney_cadd_f_window(circ, &y[n], y, lsbs, &f_bytes, 0);
    // 3) clean y[n]: y[n] ^= (y_top < x_top).
    msb_compare::lt_msbs(circ, &y[..n], &x[..n], msbs, &y[n]);
    circ.pop_section(&prev);
}

/// UNCONDITIONAL in-place pseudo-Mersenne modular sub: `y -= x mod q`.
///
/// The ctrl-free form of [`controlled_mod_sub_pm_secp256k1`]. Structural
/// mirror of [`mod_add_pm_secp256k1`]: unconditional Cuccaro sub
/// (`crate::point_add::trailmix_port::arith::ripple_add::sub`), borrowed-dirty `-f` window, uncontrolled top-msbs
/// add-overflow cleanup.
pub fn mod_sub_pm_secp256k1(circ: &mut Circuit, x: &[QReg], y: &[QReg]) {
    let n = 256usize;
    assert_eq!(x.len(), n + 1);
    assert_eq!(y.len(), n + 1);
    let prev = circ.push_section("mod_sub_pm_uncond");
    let f_bytes = super::F_SECP256K1.to_le_bytes();
    let padding = 30usize;
    let f_bitlen = 33usize;
    let lsbs = padding + f_bitlen;
    let msbs = padding;
    // 1) y -= x (n+1-bit). y[n] captures the borrow.
    crate::point_add::trailmix_port::arith::ripple_add::sub(circ, y, x);
    // 2) if borrow (y[n]): sub f from the bottom lsbs bits.
    gidney_csub_f_window(circ, &y[n], y, lsbs, &f_bytes, 0);
    // 3) clean y[n] via the top-msbs add-overflow predicate.
    crate::point_add::trailmix_port::arith::ripple_add::add_overflow_msbs_phase_correction(
        circ,
        &y[..n],
        &x[..n],
        &y[n],
        msbs,
    );
    circ.pop_section(&prev);
}

/// Controlled modular negation: `if ctrl: dst := -dst mod q`.
/// In pseudo-Mersenne: -dst mod q = q - dst = 2^n - F - dst. Compute
/// via X-flip + cadd(F+1). Approximate (drift OK).
pub fn controlled_mod_neg_pm_secp256k1(circ: &mut Circuit, ctrl: &QReg, dst: &[QReg]) {
    let n = 256usize;
    assert_eq!(dst.len(), n + 1);
    let prev = circ.push_section("mod_neg_pm");
    // -dst mod q = q - dst = (2^n - F) - dst. So:
    //   dst := -dst + 2^n - F = (2^n - 1) - dst + (1 - F)
    //        = NOT(dst) + 1 - F
    //        = NOT(dst) + (q - F - 0 + 1)  ... eh.
    //
    // Simplest CORRECT recipe: ctrl-flip all bits + ctrl-add 1 to LSB
    // (= ctrl-twos-complement-negate), then ctrl-add (q - 2^n) = -F.
    // -F is NEGATIVE, can implement as ctrl-sub F.
    // Net: dst := (-dst) - F = -dst - F.
    //   = (2^n - dst) - F (interpreting twos-comp in n+1 bits)
    //   = q - dst.   ✓
    let f_bytes = super::F_SECP256K1.to_le_bytes();
    // ctrl-flip all n bits.
    for j in 0..n {
        circ.cx(ctrl, &dst[j]);
    }
    // ctrl-inc (= add 1 controlled).
    crate::point_add::trailmix_port::arith::khattar_gidney::cinc_khattar_gidney(circ, &dst[..=n], ctrl);
    // ctrl-sub F from low (F-bitlen-+padding) bits.
    let padding = 30usize;
    let f_bitlen = 33usize;
    let lsbs = f_bitlen + padding;
    crate::point_add::trailmix_port::arith::const_add::controlled_sub_const(circ, ctrl, &dst[..lsbs], &f_bytes);
    circ.pop_section(&prev);
}

/// EXACT modular halving wrapper: delegates to
/// `rfold_mbu::mod_halve_pm_general` (commit 3718821b — exact a/2 mod q
/// for any a < 2^256, peak < 600, MBU-clean).
pub fn mod_halve_pm_secp256k1_exact(circ: &mut Circuit, a: &[QReg]) {
    let n = 256usize;
    assert_eq!(a.len(), n + 1, "register must be 257 bits");
    crate::point_add::trailmix_port::arith::rfold_mbu::mod_halve_pm_general(circ, a);
}

/// APPROXIMATE halving wrapper: same value semantics as the exact form
/// (drift ≤ 1 q from the windowed -R, identical to `mod_double_pm`).
/// Recovers the parity flag via the `flag ≡ a[255]` identity instead of
/// a full 256-bit `a >= q/2` compare; cost ~32 Toffoli/call vs ~1000.
///
/// Used by the inverse-direction `apply_bitvector` in `mod_mul_eea_rev`.
pub fn mod_halve_pm_secp256k1(circ: &mut Circuit, a: &[QReg]) {
    let n = 256usize;
    assert_eq!(a.len(), n + 1, "register must be 257 bits");
    crate::point_add::trailmix_port::arith::rfold_mbu::mod_halve_pm_general_approx_secp256k1(circ, a);
}
pub fn controlled_mod_square_sub_pm_secp256k1(
    circ: &mut Circuit,
    ctrl: &QReg,
    y_reg: &[QReg],
    output_reg: &[QReg],
) {
    let n = 256usize;
    assert_eq!(y_reg.len(), n + 1);
    assert_eq!(output_reg.len(), n + 1);
    let prev = circ.push_section("sqr_sub_pm");
    let c = circ.alloc_qreg("sqr.c");
    // Gate-for-gate inverse of the add direction.
    for i in (0..n).rev() {
        if i < n - 1 {
            mod_double_pm_secp256k1_reverse_fvented(circ, output_reg);
        }
        circ.ccx(ctrl, &y_reg[n - 1 - i], &c);
        controlled_mod_add_pm_secp256k1_reverse_vents(circ, &c, y_reg, output_reg, SQR_VENTS);
        circ.ccx(ctrl, &y_reg[n - 1 - i], &c);
    }
    for _ in 0..(n - 1) {
        mod_double_pm_secp256k1_fvented(circ, output_reg);
    }
    circ.zero_and_free(c);
    circ.pop_section(&prev);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::point_add::trailmix_port::circuit::Circuit;
    use crate::point_add::trailmix_port::compat::num_bigint::BigUint;

    fn secp256k1_q() -> BigUint {
        (BigUint::from(1u32) << 256) - BigUint::from(super::super::F_SECP256K1)
    }

    fn write_x(circ: &mut Circuit, reg: &[QReg], value: &BigUint, shot: usize) {
        let mut buf = [0u8; 32];
        let bytes = value.to_bytes_le();
        for (i, b) in bytes.iter().take(32).enumerate() {
            buf[i] = *b;
        }
        circ.sim_load_reg_bytes_shot(&reg[..256], &buf, shot);
    }

    fn read_x(sim: &crate::point_add::trailmix_port::circuit::DestroyedSimState, reg: &[QReg], shot: usize) -> BigUint {
        let mut v = BigUint::from(0u32);
        for (i, q) in reg.iter().take(256).enumerate() {
            if sim.read_bit_shot(q, shot) == 1 {
                v |= BigUint::from(1u32) << i;
            }
        }
        v
    }

    /// Quick sanity test on mod_add_const_pm_secp256k1 for one shot
    /// per (val, dst) cross.
    #[test]
    fn mod_add_const_pm_secp256k1_sanity() {
        use rand::{thread_rng, Rng};

        let q = secp256k1_q();
        let mut circ = Circuit::new();
        let dst = circ.alloc_qreg_bits("dst", 257);
        let mut rng = thread_rng();

        let mut shot_data: Vec<(BigUint, BigUint, BigUint)> = Vec::new();
        for shot in 0..64 {
            let raw_d: [u8; 32] = rng.gen();
            let raw_v: [u8; 32] = rng.gen();
            let d = BigUint::from_bytes_le(&raw_d) % &q;
            let v = BigUint::from_bytes_le(&raw_v) % &q;
            write_x(&mut circ, &dst, &d, shot);
            shot_data.push((d.clone(), v.clone(), (&d + &v) % &q));
        }
        // We can only specialize on one val per call. Use shot 0's val
        // for all shots (since the call adds the SAME constant to all
        // shots' dst).
        let val = shot_data[0].1.clone();
        mod_add_const_pm_secp256k1(&mut circ, &dst, &val);
        let (sim, detached) = circ.destroy_sim(dst);
        for (shot, (d_in, _v, _w)) in shot_data.iter().enumerate() {
            let mut got = BigUint::from(0u32);
            for (i, qq) in detached[..256].iter().enumerate() {
                if sim.read_bit_shot(qq, shot) == 1 {
                    got |= BigUint::from(1u32) << i;
                }
            }
            let got_mod = &got % &q;
            let want = (d_in + &val) % &q;
            assert_eq!(got_mod, want, "shot {shot}: dst + val mod q mismatch");
        }
    }

    /// Controlled mod-square-add: 64 random shots, value-correct
    /// mod q. Demonstrates the primitive used in Alg 1 line 10/11.
    #[test]
    fn controlled_mod_square_add_pm_secp256k1_value_correct() {
        use rand::{thread_rng, Rng};

        let q = secp256k1_q();
        let mut circ = Circuit::new();
        circ.set_max_qubit_peak(780);
        let y_reg = circ.alloc_qreg_bits("yreg", 257);
        let out_reg = circ.alloc_qreg_bits("outreg", 257);
        let ctrl = circ.alloc_qreg("ctrl");

        let mut rng = thread_rng();
        let mut shot_data: Vec<(BigUint, BigUint, bool, BigUint)> = Vec::new();
        for shot in 0..64 {
            let raw_y: [u8; 32] = rng.gen();
            let raw_o: [u8; 32] = rng.gen();
            let y = BigUint::from_bytes_le(&raw_y) % &q;
            let o = BigUint::from_bytes_le(&raw_o) % &q;
            let c = (shot & 1) == 0;
            let want = if c {
                (&o + (&y * &y) % &q) % &q
            } else {
                o.clone()
            };
            write_x(&mut circ, &y_reg, &y, shot);
            write_x(&mut circ, &out_reg, &o, shot);
            if c {
                circ.sim_load_reg_bytes_shot(std::slice::from_ref(&ctrl), &[1u8], shot);
            }
            shot_data.push((y, o, c, want));
        }

        let ops0 = circ.ops.len();
        let ccx0 = circ.ccx_emitted;
        let ccz0 = circ.ccz_emitted;
        controlled_mod_square_add_pm_secp256k1(&mut circ, &ctrl, &y_reg, &out_reg);
        let sq_ops = circ.ops.len() - ops0;
        let sq_tof = (circ.ccx_emitted - ccx0) + (circ.ccz_emitted - ccz0);
        let peak = circ.peak_qubits;
        eprintln!(
            "  cost(controlled_mod_square_add_pm_secp256k1, n=257): ops={sq_ops} tof={sq_tof} peak+={peak}"
        );

        let mut outs: Vec<QReg> = Vec::new();
        outs.extend(y_reg);
        outs.extend(out_reg);
        outs.push(ctrl);
        let (sim, detached) = circ.destroy_sim(outs);

        for (shot, (_y, _o, _c, want)) in shot_data.iter().enumerate() {
            let mut got = BigUint::from(0u32);
            for (i, qq) in detached[257..257 + 256].iter().enumerate() {
                if sim.read_bit_shot(qq, shot) == 1 {
                    got |= BigUint::from(1u32) << i;
                }
            }
            let got_mod = &got % &q;
            assert_eq!(
                got_mod, *want,
                "shot {shot}: out_post mod q mismatch (got {got:#x})"
            );
        }
    }

    /// Sample 64 random pairs (x, y) in [0, q)^2 and verify that
    /// `controlled_mod_add_pm_secp256k1` computes `y := y + ctrl * x mod q`
    /// correctly on every shot. The ~2^-30 truncation tail and
    /// equality-tail failures are extremely unlikely to fire on 64
    /// random shots; if any do, we accept them by tolerance bound.
    #[test]
    fn controlled_mod_add_pm_secp256k1_random_inputs_match_reference() {
        use rand::rngs::StdRng;
        use rand::{Rng, SeedableRng};

        let q = secp256k1_q();
        let mut circ = Circuit::new();
        let x = circ.alloc_qreg_bits("x", 257);
        let y = circ.alloc_qreg_bits("y", 257);
        let ctrl = circ.alloc_qreg("ctrl");

        let mut rng = StdRng::seed_from_u64(0xa55a_1d4a_e1d0_a17e);
        let mut expected: Vec<BigUint> = Vec::with_capacity(64);
        let mut ctrl_v = [false; 64];
        for shot in 0..64 {
            let raw_x: [u8; 32] = rng.gen();
            let raw_y: [u8; 32] = rng.gen();
            let xv = BigUint::from_bytes_le(&raw_x) % &q;
            let yv = BigUint::from_bytes_le(&raw_y) % &q;
            write_x(&mut circ, &x, &xv, shot);
            write_x(&mut circ, &y, &yv, shot);
            ctrl_v[shot] = (shot & 1) == 0;
            if ctrl_v[shot] {
                circ.sim_load_reg_bytes_shot(std::slice::from_ref(&ctrl), &[1u8], shot);
            }
            expected.push(if ctrl_v[shot] {
                (&xv + &yv) % &q
            } else {
                yv.clone()
            });
        }

        let ops0 = circ.ops.len();
        let ccx0 = circ.ccx_emitted;
        let ccz0 = circ.ccz_emitted;
        controlled_mod_add_pm_secp256k1(&mut circ, &ctrl, &x, &y);
        let cma_ops = circ.ops.len() - ops0;
        let cma_tof = (circ.ccx_emitted - ccx0) + (circ.ccz_emitted - ccz0);
        let cma_peak = circ.peak_qubits;
        eprintln!(
            "  cost(controlled_mod_add_pm_secp256k1, n=257): ops={cma_ops} tof={cma_tof} peak+={cma_peak}"
        );

        let mut outs: Vec<QReg> = Vec::new();
        outs.extend(x);
        outs.extend(y);
        outs.push(ctrl);
        let (sim, detached) = circ.destroy_sim(outs);
        let y_d = &detached[257..514];

        let mut failures = 0usize;
        for (shot, want) in expected.iter().enumerate() {
            let mut got = BigUint::from(0u32);
            for (i, q_) in y_d.iter().take(256).enumerate() {
                if sim.read_bit_shot(q_, shot) == 1 {
                    got |= BigUint::from(1u32) << i;
                }
            }
            // anc_y must be cleaned to 0.
            let anc_y = sim.read_bit_shot(&y_d[256], shot);
            if &got != want || anc_y != 0 {
                eprintln!(
                    "shot {}: ctrl={} got={:#x} want={:#x} anc_y={}",
                    shot, ctrl_v[shot], got, want, anc_y
                );
                failures += 1;
            }
        }
        // 64 random shots; expected truncation+equality tail failures
        // are 64 * (2^-30 + 2^-30) ~ 0.
        assert_eq!(
            failures, 0,
            "controlled_mod_add_pm_secp256k1: {failures}/64 random shots failed"
        );
    }

    /// Unconditional `mod_add_pm_secp256k1`: 64 random inputs in [0, q),
    /// value-correct mod q, anc bit cleaned, source `x` preserved, and
    /// the phase tracker clean.
    #[test]
    fn mod_add_pm_secp256k1_uncond_random_inputs() {
        use rand::{thread_rng, Rng};

        let q = secp256k1_q();
        let mut circ = Circuit::new();
        let x = circ.alloc_qreg_bits("x", 257);
        let y = circ.alloc_qreg_bits("y", 257);

        let mut rng = thread_rng();
        let mut xs: Vec<BigUint> = Vec::with_capacity(64);
        let mut ys: Vec<BigUint> = Vec::with_capacity(64);
        for shot in 0..64 {
            let raw_x: [u8; 32] = rng.gen();
            let raw_y: [u8; 32] = rng.gen();
            let xv = BigUint::from_bytes_le(&raw_x) % &q;
            let yv = BigUint::from_bytes_le(&raw_y) % &q;
            write_x(&mut circ, &x, &xv, shot);
            write_x(&mut circ, &y, &yv, shot);
            xs.push(xv);
            ys.push(yv);
        }

        mod_add_pm_secp256k1(&mut circ, &x, &y);
        circ.assert_phase_clean();

        let mut outs: Vec<QReg> = Vec::new();
        outs.extend(x);
        outs.extend(y);
        let (sim, detached) = circ.destroy_sim(outs);
        let (x_d, y_d) = detached.split_at(257);

        for shot in 0..64 {
            let got_y = read_x(&sim, y_d, shot) % &q;
            let got_x = read_x(&sim, x_d, shot) % &q;
            let want_y = (&xs[shot] + &ys[shot]) % &q;
            assert_eq!(got_y, want_y, "shot {shot}: y += x mod q");
            assert_eq!(got_x, xs[shot], "shot {shot}: x preserved");
            assert_eq!(sim.read_bit_shot(&y_d[256], shot), 0, "shot {shot}: y anc");
            assert_eq!(sim.read_bit_shot(&x_d[256], shot), 0, "shot {shot}: x anc");
        }
    }

    /// Unconditional `mod_sub_pm_secp256k1`: 64 random inputs, value
    /// `y - x mod q`, anc cleaned, source preserved, phase clean.
    #[test]
    fn mod_sub_pm_secp256k1_uncond_random_inputs() {
        use rand::{thread_rng, Rng};

        let q = secp256k1_q();
        let mut circ = Circuit::new();
        let x = circ.alloc_qreg_bits("x", 257);
        let y = circ.alloc_qreg_bits("y", 257);

        let mut rng = thread_rng();
        let mut xs: Vec<BigUint> = Vec::with_capacity(64);
        let mut ys: Vec<BigUint> = Vec::with_capacity(64);
        for shot in 0..64 {
            let raw_x: [u8; 32] = rng.gen();
            let raw_y: [u8; 32] = rng.gen();
            let xv = BigUint::from_bytes_le(&raw_x) % &q;
            let yv = BigUint::from_bytes_le(&raw_y) % &q;
            write_x(&mut circ, &x, &xv, shot);
            write_x(&mut circ, &y, &yv, shot);
            xs.push(xv);
            ys.push(yv);
        }

        mod_sub_pm_secp256k1(&mut circ, &x, &y);
        circ.assert_phase_clean();

        let mut outs: Vec<QReg> = Vec::new();
        outs.extend(x);
        outs.extend(y);
        let (sim, detached) = circ.destroy_sim(outs);
        let (x_d, y_d) = detached.split_at(257);

        for shot in 0..64 {
            let got_y = read_x(&sim, y_d, shot) % &q;
            let got_x = read_x(&sim, x_d, shot) % &q;
            let want_y = (&ys[shot] + &q - &xs[shot]) % &q;
            assert_eq!(got_y, want_y, "shot {shot}: y -= x mod q");
            assert_eq!(got_x, xs[shot], "shot {shot}: x preserved");
            assert_eq!(sim.read_bit_shot(&y_d[256], shot), 0, "shot {shot}: y anc");
            assert_eq!(sim.read_bit_shot(&x_d[256], shot), 0, "shot {shot}: x anc");
        }
    }

    /// Sample 64 random inputs in [0, q) and verify that
    /// `mod_double_pm_secp256k1` produces `2*x mod q` on every shot.
    /// At padding=30 the per-shot truncation-failure rate is ~2^-30,
    /// so 64 shots should all succeed with overwhelming probability.
    /// The "2x in [q, 2^n)" representation-drift is bounded by
    /// f/q ~ 2^-224 per shot.
    #[test]
    fn mod_double_pm_secp256k1_random_inputs_match_reference() {
        use rand::rngs::StdRng;
        use rand::{Rng, SeedableRng};

        let q = secp256k1_q();
        let mut circ = Circuit::new();
        let a = circ.alloc_qreg_bits("a", 257);

        let mut rng = StdRng::seed_from_u64(0xb0d1_e5c2_a1f3_0011);
        let mut expected: Vec<BigUint> = Vec::with_capacity(64);
        for shot in 0..64 {
            let raw: [u8; 32] = rng.gen();
            let x = BigUint::from_bytes_le(&raw) % &q;
            write_x(&mut circ, &a, &x, shot);
            expected.push((&x * BigUint::from(2u32)) % &q);
        }

        let ops0 = circ.ops.len();
        let ccx0 = circ.ccx_emitted;
        let ccz0 = circ.ccz_emitted;
        mod_double_pm_secp256k1(&mut circ, &a);
        let pm_ops = circ.ops.len() - ops0;
        let pm_tof = (circ.ccx_emitted - ccx0) + (circ.ccz_emitted - ccz0);
        let pm_peak = circ.peak_qubits;
        eprintln!(
            "  cost(mod_double_pm_secp256k1, n=257): ops={pm_ops} tof={pm_tof} peak+={pm_peak}"
        );

        let (sim, detached) = circ.destroy_sim(a);
        for (shot, want) in expected.iter().enumerate() {
            let got = read_x(&sim, &detached, shot);
            assert_eq!(
                &got, want,
                "shot {}: got {:#x}, want {:#x}",
                shot, got, want
            );
            // Ancilla bit a[256] must be returned to zero.
            assert_eq!(
                sim.read_bit_shot(&detached[256], shot),
                0,
                "shot {}: a[256] not cleaned",
                shot
            );
        }
    }
}
