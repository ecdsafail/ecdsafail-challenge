//! Per-curve modular-arithmetic strategy for the EC-add port.
//!
//! The dialog (inverse/mul) and the Weierstrass driver are generic over a
//! `ModRed` implementation that supplies the field-modular primitives. Two
//! impls:
//!   - [`PseudoMersenne`] (q = 2^n - f, f small/odd): wraps the
//!     generic-over-f pseudo-Mersenne CORE in `arith::schrottenloher::pm_prims`
//!     (cheap windowed `+f`). Covers secp256k1 (f = 2^32+977) and Curve25519
//!     (f = 19).
//!   - [`GenericPrime`] (no special structure): compare-and-subtract-q via the
//!     approximate top-MSB comparators. Covers SM2 and Brainpool P256.
//!
//! All registers are `n+1` bits; bit `n` is the working overflow/borrow slot
//! (0 pre/post). `x[n]` stays 0 on the addend side. We never edit the
//! `schrottenloher`/`rfold_mbu` modules (a parallel session owns them) — this
//! module only consumes their stable `pub` primitives.

use crate::point_add::trailmix_port::arith::gidney_const_adder;
use crate::point_add::trailmix_port::arith::schrottenloher::{msb_compare, pm_prims};
use crate::point_add::trailmix_port::circuit::{Circuit, QReg};
use num_bigint::BigUint;
use num_traits::Zero;

/// Approximate-comparator top-bit count and `+f`-window padding. Per-call
/// reduction-drop / compare-tie probability ~ 2^-PADDING.
const PADDING: usize = 30;

/// The field-modular primitives the dialog + driver need. Every method
/// operates in place on `n+1`-bit registers (bit `n` = overflow slot).
pub trait ModRed {
    /// Field bit width (q < 2^n).
    fn n(&self) -> usize;
    /// Little-endian bytes of q.
    fn q_bytes(&self) -> &[u8];
    /// q as a `BigUint` (for classical-constant derivations).
    fn q(&self) -> &BigUint;

    /// `a := 2*a mod q` and its exact gate-inverse.
    fn mod_double(&self, circ: &mut Circuit, a: &[QReg]);
    fn mod_double_reverse(&self, circ: &mut Circuit, a: &[QReg]);

    /// `y := y + ctrl*x mod q` and its inverse (`y := y - ctrl*x`).
    fn controlled_mod_add(&self, circ: &mut Circuit, ctrl: &QReg, x: &[QReg], y: &[QReg]);
    fn controlled_mod_add_reverse(&self, circ: &mut Circuit, ctrl: &QReg, x: &[QReg], y: &[QReg]);
    /// `y := y - ctrl*x mod q`.
    fn controlled_mod_sub(&self, circ: &mut Circuit, ctrl: &QReg, x: &[QReg], y: &[QReg]);

    /// `a := a/2 mod q` (exact). Shared generic implementation by default.
    fn mod_halve(&self, circ: &mut Circuit, a: &[QReg]) {
        mod_halve_generic(circ, a, self.n(), self.q());
    }

    /// Unconditional q-q `y += x mod q` / `y -= x mod q` (driver coord steps).
    fn mod_add_uncond(&self, circ: &mut Circuit, x: &[QReg], y: &[QReg]);
    fn mod_sub_uncond(&self, circ: &mut Circuit, x: &[QReg], y: &[QReg]);

    /// `if ctrl: dst := -dst mod q`.
    fn controlled_mod_neg(&self, circ: &mut Circuit, ctrl: &QReg, dst: &[QReg]);

    /// `if ctrl: out := out - y^2 mod q` (driver step 10).
    fn controlled_mod_square_sub(&self, circ: &mut Circuit, ctrl: &QReg, y: &[QReg], out: &[QReg]);

    /// Clear the dialog mul's "zeroed" scratch register (value ≡ 0 mod q) to
    /// |0> so it can be freed. The DEFAULT (approximate-reduce impls like
    /// [`PseudoMersenne`]) leaves the 0-residue non-canonically as `q`, so it
    /// X-clears q's set bits (0 Toffoli). Exact-reduce impls
    /// ([`GenericPrime`]) override this to a no-op — their residue is already
    /// canonical |0>, and X-clearing q here would wrongly map 0 -> q and trip
    /// `prove_zero` on free.
    fn clear_zeroed(&self, circ: &mut Circuit, reg: &[QReg]) {
        let q_bytes = self.q_bytes();
        for i in 0..self.n() {
            let byte = i / 8;
            if byte < q_bytes.len() && (q_bytes[byte] >> (i % 8)) & 1 == 1 {
                circ.x(&reg[i]);
            }
        }
    }
}

// ===========================================================================
// Shared: exact generic modular halve (used by both impls).
// ===========================================================================

/// `a := a/2 mod q`, EXACT, for any odd prime q. `a` is `n+1` bits, value in
/// `[0, 2^n)`; bit `n` is a scratch slot (0 pre/post).
///
/// Recipe: `flag = parity(a)`; if odd add q (q odd => sum even), so the carry
/// lands in `a[n]`; `right_shift` divides by 2 (carry moves into `a[n-1]`,
/// the |0> low bit is dropped). The parity flag satisfies the exact identity
/// `flag == (a_post >= ceil(q/2))`, so it is uncomputed by recomputing that
/// compare-against-constant and XOR-cancelling — a clean reversible
/// compute/use/uncompute (no bare HMR).
pub fn mod_halve_generic(circ: &mut Circuit, a: &[QReg], n: usize, q: &BigUint) {
    assert_eq!(a.len(), n + 1, "halve register must be n+1 bits");
    let prev = circ.push_section("mod_halve_generic");
    let q_bytes = q.to_bytes_le();
    // ceil(q/2) = (q+1)/2 for odd q.
    let half_ceil = (q + BigUint::from(1u32)) >> 1u32;
    let half_bytes = half_ceil.to_bytes_le();

    let flag = circ.alloc_qreg("halve.parity");
    // flag = a[0] (parity). a[0] is consumed (set to 0) by the shift below.
    circ.cx(&a[0], &flag);
    // if odd: a += q over n+1 bits (carry into a[n]); sum is even => a[0]=0.
    crate::point_add::trailmix_port::arith::const_add::controlled_add_const(circ, &flag, a, &q_bytes);
    // a := a >> 1 : value/2 into a[0..n], a[n] -> 0.
    crate::point_add::trailmix_port::arith::shift::right_shift(circ, a);
    // clean flag: XOR (a_post >= ceil(q/2)); equals flag, cancels to 0.
    crate::point_add::trailmix_port::arith::compare::compare_geq_const(circ, &a[..n], &half_bytes, &flag);
    circ.zero_and_free(flag);
    circ.pop_section(&prev);
}

// ===========================================================================
// PseudoMersenne: q = 2^n - f, f small and ODD.
// ===========================================================================

/// Pseudo-Mersenne reduction wrapper. Wraps the generic-over-f CORE in
/// `pm_prims`. `f_vents = 0` (borrowed-dirty `+f`) everywhere — peak-cheapest;
/// the measure-vented variants are a later optimization.
pub struct PseudoMersenne {
    n: usize,
    f_bytes: Vec<u8>,
    lsbs: usize,
    msbs: usize,
    q: BigUint,
    q_bytes: Vec<u8>,
}

impl PseudoMersenne {
    /// `q = 2^n - f`. `lsbs = bitlen(f) + PADDING`, `msbs = PADDING`.
    #[must_use]
    pub fn new(n: usize, f: u64) -> Self {
        assert!(f & 1 == 1, "pseudo-Mersenne cleanup requires f odd");
        let f_big = BigUint::from(f);
        let f_bitlen = f_big.bits() as usize;
        let q = (BigUint::from(1u32) << (n as u32)) - &f_big;
        PseudoMersenne {
            n,
            f_bytes: f_big.to_bytes_le(),
            lsbs: f_bitlen + PADDING,
            msbs: PADDING,
            q_bytes: q.to_bytes_le(),
            q,
        }
    }
}

impl ModRed for PseudoMersenne {
    fn n(&self) -> usize {
        self.n
    }
    fn q_bytes(&self) -> &[u8] {
        &self.q_bytes
    }
    fn q(&self) -> &BigUint {
        &self.q
    }

    fn mod_double(&self, circ: &mut Circuit, a: &[QReg]) {
        pm_prims::mod_double_pm(circ, a, &self.f_bytes, self.lsbs, 0);
    }
    fn mod_double_reverse(&self, circ: &mut Circuit, a: &[QReg]) {
        pm_prims::mod_double_pm_reverse(circ, a, &self.f_bytes, self.lsbs, 0);
    }

    fn controlled_mod_add(&self, circ: &mut Circuit, ctrl: &QReg, x: &[QReg], y: &[QReg]) {
        pm_prims::controlled_mod_add_pm(
            circ,
            ctrl,
            x,
            y,
            &self.f_bytes,
            self.lsbs,
            self.msbs,
            0,
            0,
        );
    }
    fn controlled_mod_add_reverse(&self, circ: &mut Circuit, ctrl: &QReg, x: &[QReg], y: &[QReg]) {
        pm_prims::controlled_mod_add_pm_reverse(
            circ,
            ctrl,
            x,
            y,
            &self.f_bytes,
            self.lsbs,
            self.msbs,
            0,
            0,
        );
    }
    fn controlled_mod_sub(&self, circ: &mut Circuit, ctrl: &QReg, x: &[QReg], y: &[QReg]) {
        pm_prims::controlled_mod_sub_pm(circ, ctrl, x, y, &self.f_bytes, self.lsbs, self.msbs, 0);
    }

    /// Cheap approximate pseudo-Mersenne halve (mirror of the validated secp
    /// `rfold_mbu::mod_halve_pm_general_approx`, parameterized by n/f). Adds q
    /// when odd via `+2^n` (set a[n]) and a windowed borrowed-dirty `-f`, then
    /// `right_shift`, then the `flag == a[n-1]` MBU identity (~2^-(n-bitlen f)
    /// mismatch). No `compare_geq_const` — ~lsbs Toffoli vs the generic
    /// exact halve's Theta(n log n).
    fn mod_halve(&self, circ: &mut Circuit, a: &[QReg]) {
        let n = self.n;
        assert_eq!(a.len(), n + 1, "halve register must be n+1 bits");
        let lsbs = self.lsbs;
        // dirty scratch for the windowed adder: a[lsbs .. 2*lsbs-1] (borrowed,
        // restored). Requires 2*lsbs-1 <= n.
        assert!(2 * lsbs - 1 <= n, "halve dirty window overflows register");
        let prev = circ.push_section("mod_halve_pm_approx");
        let flag = circ.alloc_qreg("halve.parity");
        circ.cx(&a[0], &flag);
        // +q if odd = +2^n (set a[n]) and -f windowed on the low `lsbs` bits.
        circ.cx(&flag, &a[n]);
        for q in &a[..lsbs] {
            circ.x(q);
        }
        gidney_const_adder::controlled_add_const_gidney(
            circ,
            &flag,
            &a[..lsbs],
            &self.f_bytes,
            &a[lsbs..lsbs + (lsbs - 1)],
        );
        for q in &a[..lsbs] {
            circ.x(q);
        }
        // divide by 2: the a[n]=flag bit shifts into a[n-1].
        crate::point_add::trailmix_port::arith::shift::right_shift(circ, a);
        // approximate MBU cleanup: flag == a[n-1] post-halve.
        circ.declare_identity(&flag, &a[n - 1]);
        let bit = circ.alloc_bit();
        circ.hmr(&flag, bit);
        circ.z_if_bit(&a[n - 1], bit);
        circ.free_bit(bit);
        circ.zero_and_free(flag);
        circ.pop_section(&prev);
    }

    fn mod_add_uncond(&self, circ: &mut Circuit, x: &[QReg], y: &[QReg]) {
        // y += x via a |1> control (avoids depending on the private +f window
        // helper; correctness-equivalent to the uncond pm add).
        let one = circ.alloc_qreg("mod_add_uncond.one");
        circ.x(&one);
        self.controlled_mod_add(circ, &one, x, y);
        circ.x(&one);
        circ.zero_and_free(one);
    }
    fn mod_sub_uncond(&self, circ: &mut Circuit, x: &[QReg], y: &[QReg]) {
        let one = circ.alloc_qreg("mod_sub_uncond.one");
        circ.x(&one);
        self.controlled_mod_sub(circ, &one, x, y);
        circ.x(&one);
        circ.zero_and_free(one);
    }

    fn controlled_mod_neg(&self, circ: &mut Circuit, ctrl: &QReg, dst: &[QReg]) {
        // -dst mod q = q - dst = (2^n - f) - dst. Recipe: ctrl-flip n bits
        // (= 2^n-1-dst), ctrl-inc (= 2^n-dst), ctrl-sub f. Valid for any
        // pseudo-Mersenne (independent of f parity).
        let n = self.n;
        assert_eq!(dst.len(), n + 1);
        let prev = circ.push_section("mod_neg_pm");
        for j in 0..n {
            circ.cx(ctrl, &dst[j]);
        }
        crate::point_add::trailmix_port::arith::khattar_gidney::cinc_khattar_gidney(circ, &dst[..=n], ctrl);
        crate::point_add::trailmix_port::arith::const_add::controlled_sub_const(circ, ctrl, &dst[..self.lsbs], &self.f_bytes);
        circ.pop_section(&prev);
    }

    fn controlled_mod_square_sub(&self, circ: &mut Circuit, ctrl: &QReg, y: &[QReg], out: &[QReg]) {
        // Horner MSB-first square, gate-inverse of the add direction. Mirrors
        // controlled_mod_square_sub_pm_secp256k1 but parameterized over n.
        let n = self.n;
        assert_eq!(y.len(), n + 1);
        assert_eq!(out.len(), n + 1);
        let prev = circ.push_section("sqr_sub_pm");
        let c = circ.alloc_qreg("sqr.c");
        for i in (0..n).rev() {
            if i < n - 1 {
                self.mod_double_reverse(circ, out);
            }
            circ.ccx(ctrl, &y[n - 1 - i], &c);
            self.controlled_mod_add_reverse(circ, &c, y, out);
            circ.ccx(ctrl, &y[n - 1 - i], &c);
        }
        for _ in 0..(n - 1) {
            self.mod_double(circ, out);
        }
        circ.zero_and_free(c);
        circ.pop_section(&prev);
    }
}

// ===========================================================================
// GenericPrime: compare-and-subtract-q modular reduction (SM2, Brainpool),
// the general-prime counterpart to the PseudoMersenne fast-fold path.
// ===========================================================================

/// Generic-prime reduction: no special structure; reduce by comparing the
/// (n+1)-bit intermediate against q on the top `msbs` bits and conditionally
/// subtracting q.
pub struct GenericPrime {
    n: usize,
    q: BigUint,
    q_bytes: Vec<u8>,
    msbs: usize,
}

impl GenericPrime {
    #[must_use]
    pub fn new(n: usize, q: BigUint) -> Self {
        assert!(!q.is_zero());
        assert!(
            (&q & BigUint::from(1u32)) == BigUint::from(1u32),
            "GenericPrime requires q odd"
        );
        assert!(q < (BigUint::from(1u32) << (n as u32)), "q must be < 2^n");
        let q_bytes = q.to_bytes_le();
        GenericPrime {
            n,
            q_bytes,
            q,
            msbs: PADDING,
        }
    }
}

impl ModRed for GenericPrime {
    fn n(&self) -> usize {
        self.n
    }
    fn q_bytes(&self) -> &[u8] {
        &self.q_bytes
    }
    fn q(&self) -> &BigUint {
        &self.q
    }
    /// `a := 2*a mod q` (Alg 6), EXACT (output in [0, q)).
    ///
    /// `left_shift` then a Roetteler conditional-subtract reduction. The reduce
    /// flag (the borrow of `2a - q`) is uncomputed EXACTLY via parity (2a is
    /// even, q is odd, so the reduced result is even iff we did NOT subtract).
    /// Exactness keeps the value canonically in [0, q) -> no drift cascades to
    /// the dialog's downstream comparators.
    fn mod_double(&self, circ: &mut Circuit, a: &[QReg]) {
        let n = self.n;
        assert_eq!(a.len(), n + 1, "double register must be n+1 bits");
        let prev = circ.push_section("mod_double_generic");
        // 1. a := 2*a_old (in [0, 2q) for a_old in [0, q)); old MSB -> a[n].
        crate::point_add::trailmix_port::arith::shift::left_shift(circ, a);
        // 2. EXACT reduce by Roetteler conditional-subtract: subtract q
        //    unconditionally (a := 2a - q, two's complement), the borrow lands
        //    in a[n] = 1[2a < q]; capture it; add q back when borrow.
        crate::point_add::trailmix_port::arith::ripple_add::sub_const(circ, a, &self.q_bytes);
        let g = circ.alloc_qreg("gp.dbl_borrow");
        circ.cx(&a[n], &g); // g = borrow = 1[2a < q]
        crate::point_add::trailmix_port::arith::const_add::controlled_add_const(circ, &g, a, &self.q_bytes);
        // a is now in [0, q) EXACTLY; a[n] was cleared by the add-back carry
        // (when borrow) or was 0 (when no borrow).
        // 3. CLEANUP (EXACT): a_result is even iff borrow (2a even, q odd), so
        //    g == NOT a_result[0]. Clean g by XOR-ing (NOT a[0]).
        circ.x(&a[0]);
        circ.cx(&a[0], &g);
        circ.x(&a[0]);
        circ.zero_and_free(g);
        circ.pop_section(&prev);
    }

    /// Exact gate-inverse of `mod_double` (net: a := a/2 mod q interpreting a
    /// as a 2*x-representation).
    fn mod_double_reverse(&self, circ: &mut Circuit, a: &[QReg]) {
        let n = self.n;
        assert_eq!(a.len(), n + 1, "double register must be n+1 bits");
        let prev = circ.push_section("mod_double_generic_rev");
        // Inverse of cleanup: reconstruct g = NOT a_result[0].
        let g = circ.alloc_qreg("gp.dbl_borrow_rev");
        circ.x(&a[0]);
        circ.cx(&a[0], &g);
        circ.x(&a[0]);
        // Inverse of the add-back: subtract q when g.
        crate::point_add::trailmix_port::arith::const_add::controlled_sub_const(circ, &g, a, &self.q_bytes);
        // Inverse of the capture: a[n] == g (borrow) now; clear g.
        circ.cx(&a[n], &g);
        circ.zero_and_free(g);
        // Inverse of the unconditional subtract: add q back.
        crate::point_add::trailmix_port::arith::ripple_add::add_const(circ, a, &self.q_bytes);
        // Inverse of left_shift.
        crate::point_add::trailmix_port::arith::shift::right_shift(circ, a);
        circ.pop_section(&prev);
    }

    /// `y := y + ctrl*x mod q` (Alg 9), EXACT (output in [0, q)).
    ///
    /// Roetteler conditional-subtract reduction, drift-FREE:
    /// 1. `controlled_hybrid_add(ctrl, y, x, 0)` -> y = S = `y_old` + ctrl*x in
    ///    n+1 bits (S in [0, 2q) for reduced inputs).
    /// 2. `controlled_sub_const(ctrl, y, q)` -> y := S - q (two's complement);
    ///    the borrow lands in y[n] = ctrl AND 1[S < q].
    /// 3. capture g = y[n]; add q back when g (`controlled_add_const(g, y, q)`),
    ///    which restores y to S (the unreduced residue) and clears y[n] via the
    ///    carry. Net: `y_result` = (S < q ? S : S - q) in [0, q).
    /// 4. CLEANUP (EXACT): g = ctrl AND 1[S < q] == ctrl AND (`y_result` >= x)
    ///    (S < q <=> we kept `y_result` = S = `y_old` + x >= x; S >= q <=>
    ///    `y_result` = `y_old` + x - q < x since `y_old` < q). Clean via
    ///    `g ^= ctrl; g ^= ctrl AND (y_result < x)` = g ^= ctrl AND (y >= x).
    /// Output stays canonically in [0, q), so no drift cascades downstream.
    fn controlled_mod_add(&self, circ: &mut Circuit, ctrl: &QReg, x: &[QReg], y: &[QReg]) {
        let n = self.n;
        assert_eq!(x.len(), n + 1, "x must be n+1 bits");
        assert_eq!(y.len(), n + 1, "y must be n+1 bits");
        let prev = circ.push_section("cma_generic");
        gidney_const_adder::controlled_hybrid_add(circ, ctrl, y, x, 0);
        crate::point_add::trailmix_port::arith::const_add::controlled_sub_const(circ, ctrl, y, &self.q_bytes);
        let g = circ.alloc_qreg("gp.add_borrow");
        circ.cx(&y[n], &g); // g = ctrl AND 1[S < q]
        crate::point_add::trailmix_port::arith::const_add::controlled_add_const(circ, &g, y, &self.q_bytes);
        // Clean g: g == ctrl AND (y_result >= x).
        circ.cx(ctrl, &g);
        msb_compare::controlled_lt_msbs(circ, ctrl, &y[..n], &x[..n], self.msbs, &g);
        circ.zero_and_free(g);
        circ.pop_section(&prev);
    }

    /// Gate-inverse of `controlled_mod_add` (net: y := y - ctrl*x mod q).
    fn controlled_mod_add_reverse(&self, circ: &mut Circuit, ctrl: &QReg, x: &[QReg], y: &[QReg]) {
        let n = self.n;
        assert_eq!(x.len(), n + 1);
        assert_eq!(y.len(), n + 1);
        let prev = circ.push_section("cma_generic_rev");
        // Inverse of step 4 cleanup: reconstruct g.
        let g = circ.alloc_qreg("gp.add_borrow_rev");
        msb_compare::controlled_lt_msbs(circ, ctrl, &y[..n], &x[..n], self.msbs, &g);
        circ.cx(ctrl, &g);
        // Inverse of the add-back: subtract q when g.
        crate::point_add::trailmix_port::arith::const_add::controlled_sub_const(circ, &g, y, &self.q_bytes);
        // Inverse of the capture: y[n] == g (borrow) now; clear g.
        circ.cx(&y[n], &g);
        circ.zero_and_free(g);
        // Inverse of step 2: add q back when ctrl.
        crate::point_add::trailmix_port::arith::const_add::controlled_add_const(circ, ctrl, y, &self.q_bytes);
        // Inverse of step 1: y -= ctrl*x (X-sandwich of controlled add).
        for q in y {
            circ.x(q);
        }
        gidney_const_adder::controlled_hybrid_add(circ, ctrl, y, x, 0);
        for q in y {
            circ.x(q);
        }
        circ.pop_section(&prev);
    }

    /// `y := y - ctrl*x mod q` (Alg 11).
    ///
    /// For a generic prime, the borrow flag's clean post-state identity is a
    /// compare against q (NOT against 2^n), so the pseudo-Mersenne
    /// `controlled_add_overflow_msbs_phase_correction_mbu` cleanup (overflow
    /// past 2^n, valid only because PM has q ~ 2^n) does NOT close. Instead we
    /// compose validated, tracker-clean `GenericPrime` primitives:
    ///   `y - ctrl*x ≡ y + ctrl*(q - x)  (mod q)`.
    /// Negate x IN PLACE (q - x), controlled-mod-add into y, then negate x back
    /// (q - (q - x) = x). No extra n+1 temp — negation is its own inverse, so x
    /// is preserved, and the controlled-add preserves its (negated) addend.
    fn controlled_mod_sub(&self, circ: &mut Circuit, ctrl: &QReg, x: &[QReg], y: &[QReg]) {
        let n = self.n;
        assert_eq!(x.len(), n + 1);
        assert_eq!(y.len(), n + 1);
        let prev = circ.push_section("cms_generic");
        let one = circ.alloc_qreg("cms.one");
        circ.x(&one);
        // x := q - x.
        self.controlled_mod_neg(circ, &one, x);
        // y += ctrl*(q - x) mod q = y - ctrl*x mod q.
        self.controlled_mod_add(circ, ctrl, x, y);
        // x := q - (q - x) = x (restored; the add preserved x = q - x).
        self.controlled_mod_neg(circ, &one, x);
        circ.x(&one);
        circ.zero_and_free(one);
        circ.pop_section(&prev);
    }

    fn mod_add_uncond(&self, circ: &mut Circuit, x: &[QReg], y: &[QReg]) {
        let one = circ.alloc_qreg("mod_add_uncond.one");
        circ.x(&one);
        self.controlled_mod_add(circ, &one, x, y);
        circ.x(&one);
        circ.zero_and_free(one);
    }
    fn mod_sub_uncond(&self, circ: &mut Circuit, x: &[QReg], y: &[QReg]) {
        let one = circ.alloc_qreg("mod_sub_uncond.one");
        circ.x(&one);
        self.controlled_mod_sub(circ, &one, x, y);
        circ.x(&one);
        circ.zero_and_free(one);
    }

    /// `if ctrl: dst := q - dst mod q`. In-place recipe (mirror PM but with the
    /// full-width fold constant): ctrl-flip n bits (=2^n-1-dst), ctrl-inc
    /// (=2^n-dst), ctrl-sub (2^n-q) over the FULL n+1 bits (=> q-dst).
    ///
    /// The sub MUST span all n+1 bits: for dst=0 the ctrl-inc carries into
    /// dst[n] (2^n-dst = 2^n), and only an n+1-bit sub of (2^n-q) borrows that
    /// bit back to 0, giving the clean residue q (dst[n]=0). A sub over only
    /// the low n bits would leave dst[n]=1 (value 2^n+q) -- which corrupts the
    /// neg-add-neg `controlled_mod_sub` whenever the subtrahend is 0 (common in
    /// the dialog, where `x_reg` starts at 0).
    fn controlled_mod_neg(&self, circ: &mut Circuit, ctrl: &QReg, dst: &[QReg]) {
        let n = self.n;
        assert_eq!(dst.len(), n + 1);
        let prev = circ.push_section("mod_neg_generic");
        for j in 0..n {
            circ.cx(ctrl, &dst[j]);
        }
        crate::point_add::trailmix_port::arith::khattar_gidney::cinc_khattar_gidney(circ, &dst[..=n], ctrl);
        let two_n_minus_q = (BigUint::from(1u32) << (n as u32)) - &self.q;
        let fold_bytes = two_n_minus_q.to_bytes_le();
        crate::point_add::trailmix_port::arith::const_add::controlled_sub_const(circ, ctrl, &dst[..=n], &fold_bytes);
        circ.pop_section(&prev);
    }

    /// `if ctrl: out := out - y^2 mod q`. Copy of `PseudoMersenne`'s Horner
    /// structure, calling the `GenericPrime` double/add primitives.
    fn controlled_mod_square_sub(&self, circ: &mut Circuit, ctrl: &QReg, y: &[QReg], out: &[QReg]) {
        let n = self.n;
        assert_eq!(y.len(), n + 1);
        assert_eq!(out.len(), n + 1);
        let prev = circ.push_section("sqr_sub_generic");
        let c = circ.alloc_qreg("sqr.c");
        for i in (0..n).rev() {
            if i < n - 1 {
                self.mod_double_reverse(circ, out);
            }
            circ.ccx(ctrl, &y[n - 1 - i], &c);
            self.controlled_mod_add_reverse(circ, &c, y, out);
            circ.ccx(ctrl, &y[n - 1 - i], &c);
        }
        for _ in 0..(n - 1) {
            self.mod_double(circ, out);
        }
        circ.zero_and_free(c);
        circ.pop_section(&prev);
    }

    /// Exact Roetteler reduce leaves the dialog's 0-residue canonically as
    /// |0>; the register is already cleared, so this is a no-op. (X-clearing q
    /// would wrongly map 0 -> q and break `prove_zero` on the subsequent free.)
    fn clear_zeroed(&self, _circ: &mut Circuit, _reg: &[QReg]) {}
}

// ===========================================================================
// Solinas: structured (generalized-Mersenne) prime p = 2^n - f_s, where f_s
// is SPARSE (few set bits / few runs). The canonical case is SM2:
//   p = 2^256 - f_s,  f_s = 2^224 + 2^96 - 2^64 + 1
//             = 2^224 + (2^32 - 1)*2^64 + 1   (positive integer; ~34 set bits,
//             bit 224, the run [64,96), bit 0). f_s < 2^225.
// ===========================================================================

/// Solinas reduction. Same EXACT Roetteler conditional-subtract value/flag
/// semantics as [`GenericPrime`] (so the cleanup identities are identical and
/// the 0-residue is canonical |0>), but the dense `±p` constant arithmetic is
/// replaced by SPARSE folds of `c_subp = 2^n + f_s = 2^{n+1} - p`.
///
/// The two structural identities that make the substitution exact, both on the
/// FULL `n+1`-bit register (bit `n` = overflow slot):
///   * `sub_const(p)`  ==  `add_const(c_subp)`   (since `2^{n+1} - p = c_subp`,
///     and two's-complement `V - p = (V + c_subp) mod 2^{n+1}`). The borrow
///     `1[V < p]` lands in bit `n` of the result, exactly as for the dense
///     `sub_const(p)`. `c_subp` is sparse, so this routes to the cheap
///     run/sparse add backend instead of the Θ(n log n) dense path.
///   * `controlled_add_const(g, p)`  ==  `controlled_sub_const_sparse(g, c_subp)`
///     realized as the X-sandwich `~a; controlled_add_const(g, c_subp); ~a`
///     (since `+p ≡ -c_subp (mod 2^{n+1})`). The X-sandwich keeps the constant
///     SPARSE (`controlled_sub_const` would negate `c_subp` to the dense `p`).
///     For `g = 0` the inner add is a no-op and `~a;~a` is identity, so it is
///     correct for both flag values; the redundant-X auto-elide cancels the
///     boundary X-X pairs on bits the inner add never touches.
///
/// `controlled_mod_neg` similarly folds `f_s` (which is sparse) rather than the
/// dense `p`: `q - dst` = flip n bits + cinc + sparse `-f_s` over `n+1` bits.
pub struct Solinas {
    n: usize,
    p: BigUint,
    p_bytes: Vec<u8>,
    /// `f_s = 2^n - p` (sparse positive integer), little-endian bytes.
    fs_bytes: Vec<u8>,
    /// `c_subp = 2^n + f_s = 2^{n+1} - p` (sparse), little-endian bytes.
    c_subp_bytes: Vec<u8>,
    msbs: usize,
}

impl Solinas {
    /// `p = 2^n - f_s` with `f_s` sparse. `p` must be odd and `< 2^n`.
    #[must_use]
    pub fn new(n: usize, p: BigUint) -> Self {
        assert!(!p.is_zero());
        assert!(
            (&p & BigUint::from(1u32)) == BigUint::from(1u32),
            "Solinas requires p odd"
        );
        assert!(p < (BigUint::from(1u32) << (n as u32)), "p must be < 2^n");
        let two_n = BigUint::from(1u32) << (n as u32);
        let fs = &two_n - &p; // f_s = 2^n - p
        let c_subp = &two_n + &fs; // 2^n + f_s = 2^{n+1} - p
        Solinas {
            n,
            p_bytes: p.to_bytes_le(),
            fs_bytes: fs.to_bytes_le(),
            c_subp_bytes: c_subp.to_bytes_le(),
            p,
            msbs: PADDING,
        }
    }

    /// Add the SPARSE constant `c_subp` UNCONDITIONALLY via a |1> control
    /// routed through `controlled_add_const`. The unconditional
    /// `crate::point_add::trailmix_port::arith::ripple_add::add_const` would route a multi-bit constant to the DENSE
    /// Θ(n log^2 n) Vandaele path (it only special-cases popcount==1), throwing
    /// away `c_subp`'s sparsity; `controlled_add_const` instead dispatches the
    /// 4-run `c_subp` to the cheap `runs_forced` carry automaton. The extra
    /// |1> ancilla + 2 X gates are negligible.
    fn add_c_subp_uncond(&self, circ: &mut Circuit, a: &[QReg]) {
        let one = circ.alloc_qreg("sol.fold_one");
        circ.x(&one);
        crate::point_add::trailmix_port::arith::const_add::controlled_add_const(circ, &one, a, &self.c_subp_bytes);
        circ.x(&one);
        circ.zero_and_free(one);
    }

    /// `a := a - p` over the full `n+1`-bit register (two's complement), i.e.
    /// `a := (a + c_subp) mod 2^{n+1}`. The borrow `1[a_pre < p]` lands in
    /// `a[n]`. Sparse. Identical value+flag to `crate::point_add::trailmix_port::arith::ripple_add::sub_const(a, p)`.
    fn sub_p(&self, circ: &mut Circuit, a: &[QReg]) {
        self.add_c_subp_uncond(circ, a);
    }

    /// Inverse of [`sub_p`]: `a := a + p` over `n+1` bits = `(a - c_subp) mod
    /// 2^{n+1}`. X-sandwich keeps `c_subp` sparse.
    fn add_p(&self, circ: &mut Circuit, a: &[QReg]) {
        for q in a {
            circ.x(q);
        }
        self.add_c_subp_uncond(circ, a);
        for q in a {
            circ.x(q);
        }
    }

    /// `if g: a := a + p` over `n+1` bits. Realized as the X-sandwich
    /// `~a; controlled_add_const(g, c_subp); ~a` = `a -= g*c_subp` =
    /// `a += g*p (mod 2^{n+1})`. Sparse; correct for g=0 (no-op).
    fn controlled_add_p(&self, circ: &mut Circuit, g: &QReg, a: &[QReg]) {
        for q in a {
            circ.x(q);
        }
        crate::point_add::trailmix_port::arith::const_add::controlled_add_const(circ, g, a, &self.c_subp_bytes);
        for q in a {
            circ.x(q);
        }
    }

    /// `if g: a := a - p` over `n+1` bits (inverse of [`controlled_add_p`]) =
    /// `a += g*c_subp (mod 2^{n+1})`. Sparse.
    fn controlled_sub_p(&self, circ: &mut Circuit, g: &QReg, a: &[QReg]) {
        crate::point_add::trailmix_port::arith::const_add::controlled_add_const(circ, g, a, &self.c_subp_bytes);
    }

    /// `if ctrl: a := a - f_s` over the given slice (`f_s` sparse). X-sandwich
    /// of `controlled_add_const(ctrl, f_s)` keeps the constant sparse.
    fn controlled_sub_fs(&self, circ: &mut Circuit, ctrl: &QReg, a: &[QReg]) {
        for q in a {
            circ.x(q);
        }
        crate::point_add::trailmix_port::arith::const_add::controlled_add_const(circ, ctrl, a, &self.fs_bytes);
        for q in a {
            circ.x(q);
        }
    }
}

impl ModRed for Solinas {
    fn n(&self) -> usize {
        self.n
    }
    fn q_bytes(&self) -> &[u8] {
        &self.p_bytes
    }
    fn q(&self) -> &BigUint {
        &self.p
    }

    /// `a := 2*a mod p`, EXACT. Identical structure to `GenericPrime::mod_double`
    /// (Roetteler conditional-subtract; flag cleaned by parity), but the
    /// `sub_const(p)` / `controlled_add_const(g, p)` are the sparse `c_subp`
    /// folds (see [`Solinas::sub_p`] / [`Solinas::controlled_add_p`]).
    fn mod_double(&self, circ: &mut Circuit, a: &[QReg]) {
        let n = self.n;
        assert_eq!(a.len(), n + 1, "double register must be n+1 bits");
        let prev = circ.push_section("mod_double_solinas");
        crate::point_add::trailmix_port::arith::shift::left_shift(circ, a);
        self.sub_p(circ, a); // a := 2a - p; borrow -> a[n]
        let g = circ.alloc_qreg("sol.dbl_borrow");
        circ.cx(&a[n], &g); // g = borrow = 1[2a < p]
        self.controlled_add_p(circ, &g, a); // add p back when borrow
                                            // CLEANUP (EXACT): 2a is even, p odd => result even iff NOT subtracted,
                                            // so g == NOT a_result[0].
        circ.x(&a[0]);
        circ.cx(&a[0], &g);
        circ.x(&a[0]);
        circ.zero_and_free(g);
        circ.pop_section(&prev);
    }

    fn mod_double_reverse(&self, circ: &mut Circuit, a: &[QReg]) {
        let n = self.n;
        assert_eq!(a.len(), n + 1, "double register must be n+1 bits");
        let prev = circ.push_section("mod_double_solinas_rev");
        let g = circ.alloc_qreg("sol.dbl_borrow_rev");
        circ.x(&a[0]);
        circ.cx(&a[0], &g);
        circ.x(&a[0]);
        self.controlled_sub_p(circ, &g, a); // inverse of add-back
        circ.cx(&a[n], &g); // a[n] == g (borrow); clear g
        circ.zero_and_free(g);
        self.add_p(circ, a); // inverse of the unconditional sub
        crate::point_add::trailmix_port::arith::shift::right_shift(circ, a);
        circ.pop_section(&prev);
    }

    /// `y := y + ctrl*x mod p`, EXACT. Identical structure to
    /// `GenericPrime::controlled_mod_add` (Roetteler conditional-subtract;
    /// flag cleaned by `ctrl AND (y_result >= x)`), but the reduction
    /// `controlled_sub_const(ctrl, p)` / `controlled_add_const(g, p)` are the
    /// sparse `c_subp` folds.
    fn controlled_mod_add(&self, circ: &mut Circuit, ctrl: &QReg, x: &[QReg], y: &[QReg]) {
        let n = self.n;
        assert_eq!(x.len(), n + 1, "x must be n+1 bits");
        assert_eq!(y.len(), n + 1, "y must be n+1 bits");
        let prev = circ.push_section("cma_solinas");
        gidney_const_adder::controlled_hybrid_add(circ, ctrl, y, x, 0);
        // y := S - ctrl*p; borrow (when ctrl) lands in y[n].
        self.controlled_sub_p(circ, ctrl, y);
        let g = circ.alloc_qreg("sol.add_borrow");
        circ.cx(&y[n], &g); // g = ctrl AND 1[S < p]
        self.controlled_add_p(circ, &g, y); // add p back when borrow
                                            // Clean g: g == ctrl AND (y_result >= x).
        circ.cx(ctrl, &g);
        msb_compare::controlled_lt_msbs(circ, ctrl, &y[..n], &x[..n], self.msbs, &g);
        circ.zero_and_free(g);
        circ.pop_section(&prev);
    }

    fn controlled_mod_add_reverse(&self, circ: &mut Circuit, ctrl: &QReg, x: &[QReg], y: &[QReg]) {
        let n = self.n;
        assert_eq!(x.len(), n + 1);
        assert_eq!(y.len(), n + 1);
        let prev = circ.push_section("cma_solinas_rev");
        let g = circ.alloc_qreg("sol.add_borrow_rev");
        msb_compare::controlled_lt_msbs(circ, ctrl, &y[..n], &x[..n], self.msbs, &g);
        circ.cx(ctrl, &g);
        self.controlled_sub_p(circ, &g, y); // inverse of add-back
        circ.cx(&y[n], &g); // y[n] == g (borrow); clear g
        circ.zero_and_free(g);
        self.controlled_add_p(circ, ctrl, y); // inverse of the conditional sub-p
                                              // Inverse of step 1: y -= ctrl*x (X-sandwich of controlled add).
        for q in y {
            circ.x(q);
        }
        gidney_const_adder::controlled_hybrid_add(circ, ctrl, y, x, 0);
        for q in y {
            circ.x(q);
        }
        circ.pop_section(&prev);
    }

    /// `y := y - ctrl*x mod p` via `y + ctrl*(p - x)` (negate-add-negate),
    /// identical to `GenericPrime`; the sub-primitives are the sparse Solinas
    /// versions.
    fn controlled_mod_sub(&self, circ: &mut Circuit, ctrl: &QReg, x: &[QReg], y: &[QReg]) {
        let n = self.n;
        assert_eq!(x.len(), n + 1);
        assert_eq!(y.len(), n + 1);
        let prev = circ.push_section("cms_solinas");
        let one = circ.alloc_qreg("cms.one");
        circ.x(&one);
        self.controlled_mod_neg(circ, &one, x);
        self.controlled_mod_add(circ, ctrl, x, y);
        self.controlled_mod_neg(circ, &one, x);
        circ.x(&one);
        circ.zero_and_free(one);
        circ.pop_section(&prev);
    }

    fn mod_add_uncond(&self, circ: &mut Circuit, x: &[QReg], y: &[QReg]) {
        let one = circ.alloc_qreg("mod_add_uncond.one");
        circ.x(&one);
        self.controlled_mod_add(circ, &one, x, y);
        circ.x(&one);
        circ.zero_and_free(one);
    }
    fn mod_sub_uncond(&self, circ: &mut Circuit, x: &[QReg], y: &[QReg]) {
        let one = circ.alloc_qreg("mod_sub_uncond.one");
        circ.x(&one);
        self.controlled_mod_sub(circ, &one, x, y);
        circ.x(&one);
        circ.zero_and_free(one);
    }

    /// `if ctrl: dst := p - dst mod p`. Same recipe as `GenericPrime`
    /// (ctrl-flip n bits + cinc + sub `(2^n - p) = f_s` over `n+1` bits) but
    /// the `f_s` fold is sparse (see [`Solinas::controlled_sub_fs`]).
    fn controlled_mod_neg(&self, circ: &mut Circuit, ctrl: &QReg, dst: &[QReg]) {
        let n = self.n;
        assert_eq!(dst.len(), n + 1);
        let prev = circ.push_section("mod_neg_solinas");
        for j in 0..n {
            circ.cx(ctrl, &dst[j]);
        }
        crate::point_add::trailmix_port::arith::khattar_gidney::cinc_khattar_gidney(circ, &dst[..=n], ctrl);
        self.controlled_sub_fs(circ, ctrl, &dst[..=n]);
        circ.pop_section(&prev);
    }

    /// `if ctrl: out := out - y^2 mod p`. Horner, identical to `GenericPrime`
    /// (the Solinas double/add primitives are substituted automatically).
    fn controlled_mod_square_sub(&self, circ: &mut Circuit, ctrl: &QReg, y: &[QReg], out: &[QReg]) {
        let n = self.n;
        assert_eq!(y.len(), n + 1);
        assert_eq!(out.len(), n + 1);
        let prev = circ.push_section("sqr_sub_solinas");
        let c = circ.alloc_qreg("sqr.c");
        for i in (0..n).rev() {
            if i < n - 1 {
                self.mod_double_reverse(circ, out);
            }
            circ.ccx(ctrl, &y[n - 1 - i], &c);
            self.controlled_mod_add_reverse(circ, &c, y, out);
            circ.ccx(ctrl, &y[n - 1 - i], &c);
        }
        for _ in 0..(n - 1) {
            self.mod_double(circ, out);
        }
        circ.zero_and_free(c);
        circ.pop_section(&prev);
    }

    /// Exact reduce leaves the 0-residue canonically as |0> (same as
    /// `GenericPrime`), so clearing is a no-op.
    fn clear_zeroed(&self, _circ: &mut Circuit, _reg: &[QReg]) {}
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::point_add::trailmix_port::circuit::QReg;
    use num_bigint::BigUint;
    use rand::RngCore;

    fn rand_lt(rng: &mut dyn RngCore, q: &BigUint) -> BigUint {
        let mut b = [0u8; 32];
        loop {
            rng.fill_bytes(&mut b);
            let v = BigUint::from_bytes_le(&b);
            if &v < q {
                return v;
            }
        }
    }
    fn load(circ: &mut Circuit, reg: &[QReg], v: &BigUint, n: usize, shot: usize) {
        let mut buf = [0u8; 32];
        for (i, x) in v.to_bytes_le().iter().take(32).enumerate() {
            buf[i] = *x;
        }
        circ.sim_load_reg_bytes_shot(&reg[..n], &buf, shot);
    }
    fn read(
        sim: &crate::point_add::trailmix_port::circuit::DestroyedSimState,
        reg: &[QReg],
        n: usize,
        shot: usize,
    ) -> BigUint {
        let mut v = BigUint::zero();
        for i in 0..n {
            if sim.read_bit_shot(&reg[i], shot) == 1 {
                v |= BigUint::from(1u32) << i;
            }
        }
        v
    }

    /// Validate a `ModRed` impl's primitives + halve against the classical
    /// reference (value mod q) for the given curve, 64 random shots each.
    /// destroy_sim enforces clean teardown (ancillae back to |0>).
    fn check_curve<R: ModRed>(name: &str, mr: &R) {
        let n = mr.n();
        let q = mr.q().clone();
        let inv2 = (&q + BigUint::from(1u32)) >> 1u32; // (q+1)/2 = 2^-1 mod q
        let mut rng = rand::thread_rng();

        // mod_double: a -> 2a mod q.
        {
            let mut circ = Circuit::new();
            let a = circ.alloc_qreg_bits("a", n + 1);
            let mut want = Vec::new();
            for shot in 0..64 {
                let v = rand_lt(&mut rng, &q);
                load(&mut circ, &a, &v, n, shot);
                want.push((&v * 2u32) % &q);
            }
            mr.mod_double(&mut circ, &a);
            let (sim, d) = circ.destroy_sim(a);
            for shot in 0..64 {
                assert_eq!(
                    read(&sim, &d, n, shot) % &q,
                    want[shot],
                    "{name} mod_double shot {shot}"
                );
                assert_eq!(
                    sim.read_bit_shot(&d[n], shot),
                    0,
                    "{name} mod_double overflow slot"
                );
            }
        }
        // mod_halve: a -> a/2 mod q (exact).
        {
            let mut circ = Circuit::new();
            let a = circ.alloc_qreg_bits("a", n + 1);
            let mut want = Vec::new();
            for shot in 0..64 {
                let v = rand_lt(&mut rng, &q);
                load(&mut circ, &a, &v, n, shot);
                want.push((&v * &inv2) % &q);
            }
            mr.mod_halve(&mut circ, &a);
            let (sim, d) = circ.destroy_sim(a);
            for shot in 0..64 {
                assert_eq!(
                    read(&sim, &d, n, shot) % &q,
                    want[shot],
                    "{name} mod_halve shot {shot}"
                );
                assert_eq!(
                    sim.read_bit_shot(&d[n], shot),
                    0,
                    "{name} mod_halve overflow slot"
                );
            }
        }
        // controlled_mod_add (ctrl=1): y -> y + x mod q.
        {
            let mut circ = Circuit::new();
            let x = circ.alloc_qreg_bits("x", n + 1);
            let y = circ.alloc_qreg_bits("y", n + 1);
            let ctrl = circ.alloc_qreg("ctrl");
            for shot in 0..64 {
                circ.sim_load_reg_bytes_shot(std::slice::from_ref(&ctrl), &[1u8], shot);
            }
            let mut want = Vec::new();
            for shot in 0..64 {
                let xv = rand_lt(&mut rng, &q);
                let yv = rand_lt(&mut rng, &q);
                load(&mut circ, &x, &xv, n, shot);
                load(&mut circ, &y, &yv, n, shot);
                want.push((&xv + &yv) % &q);
            }
            mr.controlled_mod_add(&mut circ, &ctrl, &x, &y);
            let mut outs = x;
            outs.extend(y);
            outs.push(ctrl);
            let (sim, d) = circ.destroy_sim(outs);
            for shot in 0..64 {
                let yv = read(&sim, &d[n + 1..2 * (n + 1)], n, shot);
                assert_eq!(yv % &q, want[shot], "{name} controlled_mod_add shot {shot}");
            }
        }
        // controlled_mod_sub (ctrl=1): y -> y - x mod q.
        {
            let mut circ = Circuit::new();
            let x = circ.alloc_qreg_bits("x", n + 1);
            let y = circ.alloc_qreg_bits("y", n + 1);
            let ctrl = circ.alloc_qreg("ctrl");
            for shot in 0..64 {
                circ.sim_load_reg_bytes_shot(std::slice::from_ref(&ctrl), &[1u8], shot);
            }
            let mut want = Vec::new();
            for shot in 0..64 {
                let xv = rand_lt(&mut rng, &q);
                let yv = rand_lt(&mut rng, &q);
                load(&mut circ, &x, &xv, n, shot);
                load(&mut circ, &y, &yv, n, shot);
                want.push((&yv + &q - &xv) % &q);
            }
            mr.controlled_mod_sub(&mut circ, &ctrl, &x, &y);
            let mut outs = x;
            outs.extend(y);
            outs.push(ctrl);
            let (sim, d) = circ.destroy_sim(outs);
            for shot in 0..64 {
                let yv = read(&sim, &d[n + 1..2 * (n + 1)], n, shot);
                assert_eq!(yv % &q, want[shot], "{name} controlled_mod_sub shot {shot}");
            }
        }
        // controlled_mod_neg (ctrl=1): dst -> -dst mod q.
        {
            let mut circ = Circuit::new();
            let dst = circ.alloc_qreg_bits("dst", n + 1);
            let ctrl = circ.alloc_qreg("ctrl");
            for shot in 0..64 {
                circ.sim_load_reg_bytes_shot(std::slice::from_ref(&ctrl), &[1u8], shot);
            }
            let mut want = Vec::new();
            for shot in 0..64 {
                let v = rand_lt(&mut rng, &q);
                load(&mut circ, &dst, &v, n, shot);
                want.push((&q - &v) % &q);
            }
            mr.controlled_mod_neg(&mut circ, &ctrl, &dst);
            let mut outs = dst;
            outs.push(ctrl);
            let (sim, d) = circ.destroy_sim(outs);
            for shot in 0..64 {
                assert_eq!(
                    read(&sim, &d, n, shot) % &q,
                    want[shot],
                    "{name} controlled_mod_neg shot {shot}"
                );
            }
        }
        // controlled_mod_square_sub (ctrl=1): out -> out - y^2 mod q.
        {
            let mut circ = Circuit::new();
            let y = circ.alloc_qreg_bits("y", n + 1);
            let out = circ.alloc_qreg_bits("out", n + 1);
            let ctrl = circ.alloc_qreg("ctrl");
            for shot in 0..64 {
                circ.sim_load_reg_bytes_shot(std::slice::from_ref(&ctrl), &[1u8], shot);
            }
            let mut want = Vec::new();
            for shot in 0..64 {
                let yv = rand_lt(&mut rng, &q);
                let ov = rand_lt(&mut rng, &q);
                load(&mut circ, &y, &yv, n, shot);
                load(&mut circ, &out, &ov, n, shot);
                let y2 = (&yv * &yv) % &q;
                want.push((&ov + &q - &y2) % &q);
            }
            mr.controlled_mod_square_sub(&mut circ, &ctrl, &y, &out);
            let mut outs = y;
            outs.extend(out);
            outs.push(ctrl);
            let (sim, d) = circ.destroy_sim(outs);
            for shot in 0..64 {
                let ov = read(&sim, &d[n + 1..2 * (n + 1)], n, shot);
                assert_eq!(
                    ov % &q,
                    want[shot],
                    "{name} controlled_mod_square_sub shot {shot}"
                );
            }
        }
    }

    #[test]
    fn pseudo_mersenne_prims_secp256k1() {
        check_curve("secp256k1", &PseudoMersenne::new(256, (1u64 << 32) + 977));
    }

    #[test]
    fn pseudo_mersenne_prims_curve25519() {
        check_curve("curve25519", &PseudoMersenne::new(255, 19));
    }

    #[test]
    fn generic_prime_prims_sm2() {
        check_curve(
            "sm2",
            &GenericPrime::new(256, crate::point_add::trailmix_port::ec::curves::params::sm2().p),
        );
    }

    #[test]
    fn generic_prime_prims_brainpool() {
        check_curve(
            "brainpool",
            &GenericPrime::new(256, crate::point_add::trailmix_port::ec::curves::params::brainpoolp256r1().p),
        );
    }

    #[test]
    fn solinas_prims_sm2() {
        check_curve(
            "sm2-solinas",
            &Solinas::new(256, crate::point_add::trailmix_port::ec::curves::params::sm2().p),
        );
    }

    /// Regression: the dialog feeds the SUBTRAHEND x = 0 to controlled_mod_sub
    /// systematically (x_reg starts at 0). The neg-add-neg composition went
    /// through controlled_mod_neg(0); a low-n-bits-only fold there left dst[n]=1
    /// (value 2^n+q) and broke the e2e (gp.add_borrow prove_zero). The n+1-bit
    /// fold gives a clean q (= -0 mod q). Verify x=0 (and the y==x boundary)
    /// for both GenericPrime curves on EVERY shot. Fails before the fix.
    fn check_sub_neg_edges<R: ModRed>(name: &str, mr: &R) {
        let n = mr.n();
        let q = mr.q().clone();
        let mut rng = rand::thread_rng();

        // controlled_mod_sub with x = 0 (subtract 0 -> y unchanged).
        {
            let mut circ = Circuit::new();
            let x = circ.alloc_qreg_bits("x", n + 1);
            let y = circ.alloc_qreg_bits("y", n + 1);
            let ctrl = circ.alloc_qreg("ctrl");
            let mut want = Vec::new();
            for shot in 0..64 {
                circ.sim_load_reg_bytes_shot(std::slice::from_ref(&ctrl), &[1u8], shot);
                let yv = rand_lt(&mut rng, &q);
                load(&mut circ, &y, &yv, n, shot); // x stays 0
                want.push(yv);
            }
            mr.controlled_mod_sub(&mut circ, &ctrl, &x, &y);
            let mut outs = x;
            outs.extend(y);
            outs.push(ctrl);
            let (sim, d) = circ.destroy_sim(outs);
            for shot in 0..64 {
                let yv = read(&sim, &d[n + 1..2 * (n + 1)], n, shot);
                assert_eq!(yv % &q, want[shot], "{name} sub-x0 shot {shot}");
            }
        }
        // controlled_mod_neg with dst = 0 (-0 mod q == 0).
        {
            let mut circ = Circuit::new();
            let dst = circ.alloc_qreg_bits("dst", n + 1);
            let ctrl = circ.alloc_qreg("ctrl");
            for shot in 0..64 {
                circ.sim_load_reg_bytes_shot(std::slice::from_ref(&ctrl), &[1u8], shot);
            }
            mr.controlled_mod_neg(&mut circ, &ctrl, &dst);
            let mut outs = dst;
            outs.push(ctrl);
            let (sim, d) = circ.destroy_sim(outs);
            for shot in 0..64 {
                assert_eq!(
                    read(&sim, &d, n, shot) % &q,
                    BigUint::zero(),
                    "{name} neg-0 shot {shot}"
                );
                assert_eq!(
                    sim.read_bit_shot(&d[n], shot),
                    0,
                    "{name} neg-0 overflow slot shot {shot}"
                );
            }
        }
        // controlled_mod_sub with y == x (result 0).
        {
            let mut circ = Circuit::new();
            let x = circ.alloc_qreg_bits("x", n + 1);
            let y = circ.alloc_qreg_bits("y", n + 1);
            let ctrl = circ.alloc_qreg("ctrl");
            for shot in 0..64 {
                circ.sim_load_reg_bytes_shot(std::slice::from_ref(&ctrl), &[1u8], shot);
                let v = rand_lt(&mut rng, &q);
                load(&mut circ, &x, &v, n, shot);
                load(&mut circ, &y, &v, n, shot);
            }
            mr.controlled_mod_sub(&mut circ, &ctrl, &x, &y);
            let mut outs = x;
            outs.extend(y);
            outs.push(ctrl);
            let (sim, d) = circ.destroy_sim(outs);
            for shot in 0..64 {
                let yv = read(&sim, &d[n + 1..2 * (n + 1)], n, shot);
                assert_eq!(yv % &q, BigUint::zero(), "{name} sub-yx shot {shot}");
            }
        }
    }

    #[test]
    fn generic_prime_sub_neg_edges_sm2() {
        check_sub_neg_edges(
            "sm2",
            &GenericPrime::new(256, crate::point_add::trailmix_port::ec::curves::params::sm2().p),
        );
    }

    #[test]
    fn generic_prime_sub_neg_edges_brainpool() {
        check_sub_neg_edges(
            "brainpool",
            &GenericPrime::new(256, crate::point_add::trailmix_port::ec::curves::params::brainpoolp256r1().p),
        );
    }

    #[test]
    fn solinas_sub_neg_edges_sm2() {
        check_sub_neg_edges(
            "sm2-solinas",
            &Solinas::new(256, crate::point_add::trailmix_port::ec::curves::params::sm2().p),
        );
    }
}
