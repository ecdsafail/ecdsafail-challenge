//! Reversible unpacked PZ inversion as a bit-by-bit pipelined state machine
//! (design reference: `scripts/kaliski_test.py` `pz_big_step`). This supersedes
//! the full-division `shrunken_pz_primitives` module, whose coarser granularity
//! needed a fat quotient pad and did not handle large termination quotients.
//!
//! Per iteration (fixed count ~= sum of quotient bitlengths), gated on the state
//! flags so termination is intrinsic (no separate counter):
//!   DIVISION substep:  s = bitlen(A)-bitlen(B); align B<<s; if A>=B { A-=B;
//!                      `q_div` ^= 1<<s }; restore B>>s. A<B => `div_active=0`.
//!   MULTIPLY substep (pipelined): s = `ctz(q_mul)`; clear it; a += b<<s; restore.
//!                      `q_mul==0` => swap a,b; flip parity; `mul_active=0`.
//!   TRANSITION: q_div->q_mul; swap A,B; divide builds the NEXT quotient while
//!               the multiply drains the PREVIOUS. q pads are TINY (one quotient).
//! All shifts are `controlled_cyclic_rotate` (rotate-in-place, fixed width).
//! Up front: normalize x -> min(x, P-x) (sgn); final a corrected by parity ^ sgn.

#![allow(dead_code)]

use crate::point_add::trailmix_port::circuit::{Circuit, QReg};
use crate::point_add::trailmix_port::inversion::shrunken_pz_primitives::borrow_compare_refs;

fn env_usize(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(default)
}

fn trailmix_srot_width() -> usize {
    // The generated schedule's shift bounds need six bits on valid samples.
    // Keep an env override for experiments.
    env_usize("TRAILMIX_SROT_W", 6).max(1)
}

fn trailmix_counter_width() -> usize {
    if std::env::var("TRAILMIX_NO_COUNTER").ok().as_deref() == Some("1") {
        0
    } else {
        env_usize("TRAILMIX_COUNTER_W", 10)
    }
}

fn trailmix_q_width(wq: usize) -> usize {
    let w = wq.max(1);
    std::env::var("TRAILMIX_Q_CAP")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .map_or(w, |cap| w.min(cap.max(1)))
}

/// Per-step quotient width with SELECTIVE peak-targeting.
///
/// The global qubit peak at a `shrunken_pz` step is
///   2*max(wa,wb) + 2*max(wca,wcb) + q_width + FIXED.
/// A blunt global `TRAILMIX_Q_CAP` clamps q on ALL ~490 steps (most have
/// universal q in 23..38), but only the peak-binding step(s) need a smaller q
/// to lower the global peak. Clamping the rest just manufactures classical
/// misses (overflowed quotients) without helping the peak.
///
/// `TRAILMIX_Q_TARGET=T` instead gives each step a budget so that its working
/// width never exceeds T: `q <= T - 2*max(wa,wb) - 2*max(wca,wcb)`. Steps whose
/// other registers are small keep their full natural q (no miss); only the
/// wide-carry peak step(s) get q trimmed, and only by the minimum needed.
/// Falls back to `trailmix_q_width` (global cap) when `TRAILMIX_Q_TARGET` unset.
/// Cap the shared A/B register width (both A and B are resized to max(wa,wb)).
/// `TRAILMIX_AB_CAP` trims it on the steps where it would otherwise bind the peak.
fn trailmix_ab_width(wab: usize) -> usize {
    let w = wab.max(1);
    std::env::var("TRAILMIX_AB_CAP")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .map_or(w, |c| w.min(c.max(1)))
}

/// Cap the shared ca/cb cofactor register width (both resized to max(wca,wcb)).
/// `TRAILMIX_CACB_CAP` trims the dominant 2*245 carry pair at the peak step.
fn trailmix_cacb_width(wcacb: usize) -> usize {
    let w = wcacb.max(1);
    std::env::var("TRAILMIX_CACB_CAP")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .map_or(w, |c| w.min(c.max(1)))
}

fn trailmix_q_width_step(wq: usize, wa: usize, wb: usize, wca: usize, wcb: usize) -> usize {
    let natural = wq.max(1);
    let target = std::env::var("TRAILMIX_Q_TARGET")
        .ok()
        .and_then(|s| s.parse::<usize>().ok());
    let Some(target) = target else {
        return trailmix_q_width(wq);
    };
    // q budget is computed from the (possibly capped) A/B and ca/cb widths so the
    // working width 2*ab + 2*cacb + q meets `target` consistently with the resizes.
    let other = 2 * trailmix_ab_width(wa.max(wb)) + 2 * trailmix_cacb_width(wca.max(wcb));
    let budget = target.saturating_sub(other).max(1);
    // Still honor a global Q_CAP if both are set (take the tighter bound).
    let capped = natural.min(budget);
    std::env::var("TRAILMIX_Q_CAP")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .map_or(capped, |cap| capped.min(cap.max(1)))
        .max(1)
}

fn compute_active(c: &mut Circuit, counter: &[QReg]) -> QReg {
    let active = c.alloc_qreg("active");
    if counter.is_empty() {
        c.x(&active);
    } else {
        or_is_zero(c, counter, &active);
    }
    active
}

fn uncompute_active(c: &mut Circuit, counter: &[QReg], active: &QReg) {
    if counter.is_empty() {
        c.x(active);
    } else {
        or_is_zero(c, counter, active);
    }
}

/// `p + 1` (secp256k1 base field prime) as 33 LE bytes.
fn p_plus_1_bytes() -> Vec<u8> {
    vec![
        0x30, 0xfc, 0xff, 0xff, 0xfe, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x00,
    ]
}

/// Controlled field-negate `a := (p - a) mod p` IFF `g` (a in [0,p), 257-bit).
/// Self-inverse. `~a + (p+1) ≡ p - a (mod 2^257)`; canonical for a in [1,p).
/// (Relocated from `kaliski_spooky::unpacked` so `shrunken_pz` has no spooky-Kaliski dep.)
pub fn controlled_field_neg(c: &mut Circuit, g: &QReg, a: &[QReg]) {
    use crate::point_add::trailmix_port::arith::const_add::controlled_add_const;
    for q in a {
        c.cx(g, q);
    }
    controlled_add_const(c, g, a, &p_plus_1_bytes());
}

/// `s += bitlen(a) - bitlen(b)` (clz diff), bound by `bound`. After alignment in
/// the division substep, s is the shift to apply. Inverse: swap a,b.
/// LEAN `bit_length`: `s += bitlen(src)` (or `-=` if dec), via a reversible
/// prefix-AND ladder + gray-code deposit -- ~2n ccx (ladder build+unbuild) with
/// NO per-row position-equality. Supersedes the first-hit scan (~38 tof/row from
/// the per-row `toggle_on_cursor_eq_const` uncompute of `is_hit`).
///
/// Construction (MSB-first running flag `f_i` = "no 1 bit strictly above i"):
///   - prefix-AND ladder over ~src (X-bracketed) gives every `f_i` as a ladder
///     qubit, fully reversibly (fwd builds, rev unbuilds).
///   - deposit pos (init = n) ^= (i ^ (i+1)) gated on `f_i`, for i = n-1..0. The
///     gray differences telescope: pos collapses to the MSB index p (= bitlen-1).
///   - s += (pos + 1)  [bitlen]; then uncompute pos (re-run deposit) + ladder.
///
/// PRE: src nonzero (EEA gcd / nonzero quotient pad). For src==0 this returns
/// bitlen=1 (pos stays 0, +1); callers must not pass an all-zero src.
/// _middle core. Builds the prefix-AND ladder over ~src, deposits the MSB index
/// (= bitlen-1) into the caller's `pos` register (PRE: pos = |n>) in the FORWARD
/// sweep, runs `body` (which sees pos = MSB index), then unbuilds.
///
/// `body` returns whether the deposit should be UNDONE on the reverse sweep:
///   - `false` (DEFAULT, 3n): pos is KEPT at the MSB index -- the caller owns it
///     and must clear it later (e.g. via the SM's reverse). One consume = 3n.
///   - `true` (4n): the deposit is re-run on the reverse, returning pos to |n>.
///     Use when pos is a throwaway temp whose value was folded elsewhere in body.
///
/// The gray-code deposit is pure XOR (CX gated on a single flag materialized from
/// the prefix-AND with one ccx, then HMR-freed) -- so each consume is 1 toffoli
/// per position. Prefix build+unbuild = 2n; consume = n/sweep.
fn bit_length_lean_middle(
    circ: &mut Circuit,
    src: &[&QReg],
    pos: &[QReg],
    body: impl FnOnce(&mut Circuit) -> bool,
) {
    use crate::point_add::trailmix_port::arith::khattar_gidney::{kg_prefix_ancilla_count, KgPrefixAnd};
    let n = src.len();
    if n == 0 {
        body(circ);
        return;
    }
    // ~src (X-bracket); the prefix-AND reads the complemented bits.
    for q in src {
        circ.x(q);
    }
    // q = ~src MSB-first: q[j] = ~src[n-1-j]. The log*-ancilla KG streaming
    // prefix-AND gives, at layer i, AND(ctrls) = AND(q[0..i]) = "no 1 in top i
    // positions" = f_k ("no 1 strictly above k") for k = n-1-i. ctrls is 1-2 qubits
    // (KG conditionally-clean form), so the deposit consumes those controls directly.
    let qbits: Vec<&QReg> = src.iter().rev().copied().collect();
    let nanc = kg_prefix_ancilla_count(n);
    let anc_owned = circ.alloc_qreg_bits("bll.kganc", nanc);
    let anc: Vec<&QReg> = anc_owned.iter().collect();

    // Deposit at layer i (position k = n-1-i): gray-XOR (k ^ (k+1)) into pos gated
    // on f_k = AND(ctrls). For two controls this spends a CCX per gray bit and saves
    // the temporary flag qubit; that is the qubit-first tradeoff used in this branch.
    // pos starts at |n>; the gray differences telescope it to the MSB index p.
    // Self-inverse, so reverse undoes pos to |n>.
    fn prefix_controlled_x(circ: &mut Circuit, ctrls: &[&QReg], target: &QReg) {
        match ctrls {
            [] => circ.x(target),
            [c] => circ.cx(c, target),
            [a, b] => circ.ccx(a, b, target),
            _ => unreachable!("KG prefix ctrls is <=2 qubits"),
        }
    }

    fn deposit_step(
        circ: &mut Circuit,
        i: usize,
        ctrls: &[&QReg],
        pos: &[QReg],
        n: usize,
    ) {
        if i >= n {
            return; // i == n is the empty (k = -1) layer
        }
        let k = n - 1 - i;
        let gd = k ^ (k + 1);
        let bits: Vec<usize> = (0..pos.len()).filter(|&b| (gd >> b) & 1 == 1).collect();
        if bits.is_empty() {
            return;
        }
        for &b in &bits {
            prefix_controlled_x(circ, ctrls, &pos[b]);
        }
    }

    let kg = KgPrefixAnd::new(&qbits, &anc);
    let done = kg.forward(circ, |c, i, ctrls| deposit_step(c, i, ctrls, pos, n)); // pos -> p
    let clean = body(circ);
    if clean {
        // 4n: re-run the deposit on the reverse, returning pos to |n>.
        done.reverse(circ, |c, i, ctrls| deposit_step(c, i, ctrls, pos, n));
    } else {
        // 3n: unbuild the prefix only; pos stays at the MSB index (caller-owned).
        done.reverse(circ, |_, _, _| {});
    }
    drop(anc);
    for q in anc_owned {
        circ.zero_and_free(q);
    }
    for q in src {
        circ.x(q);
    }
}

/// Computes one bit of the modulo-`2^w` MSB-position value produced by
/// [`bit_length_lean_middle`], using one caller-owned target bit instead of a
/// full position register.
fn xor_bit_length_pos_bit(circ: &mut Circuit, src: &[&QReg], bit: usize, out: &QReg) {
    use crate::point_add::trailmix_port::arith::khattar_gidney::{
        kg_prefix_ancilla_count, KgPrefixAnd,
    };
    let n = src.len();
    if n == 0 {
        return;
    }
    for q in src {
        circ.x(q);
    }
    if ((n >> bit) & 1) == 1 {
        circ.x(out);
    }
    let qbits: Vec<&QReg> = src.iter().rev().copied().collect();
    let nanc = kg_prefix_ancilla_count(n);
    let anc_owned = circ.alloc_qreg_bits("bll.bit.kganc", nanc);
    let anc: Vec<&QReg> = anc_owned.iter().collect();

    fn deposit_step(
        circ: &mut Circuit,
        i: usize,
        ctrls: &[&QReg],
        bit: usize,
        out: &QReg,
        n: usize,
    ) {
        if i >= n {
            return;
        }
        let k = n - 1 - i;
        if ((k ^ (k + 1)) >> bit) & 1 == 0 {
            return;
        }
        match ctrls {
            [] => circ.x(out),
            [c] => circ.cx(c, out),
            [a, b] => circ.ccx(a, b, out),
            _ => unreachable!("KG prefix ctrls is <=2 qubits"),
        }
    }

    let kg = KgPrefixAnd::new(&qbits, &anc);
    let done = kg.forward(circ, |c, i, ctrls| {
        deposit_step(c, i, ctrls, bit, out, n)
    });
    done.reverse(circ, |_, _, _| {});
    drop(anc);
    for q in anc_owned {
        circ.zero_and_free(q);
    }
    for q in src {
        circ.x(q);
    }
}

fn fold_pos_bits_into_acc(circ: &mut Circuit, src: &[&QReg], acc: &[QReg], add: bool) {
    if acc.len() <= 5 && std::env::var("TRAILMIX_PREFIX_COUNT_FOLD").ok().as_deref() != Some("0") {
        fold_pos_by_prefix_count(circ, src, acc, add);
        return;
    }
    let bit = circ.alloc_qreg("clzm.posbit");
    for j in 0..acc.len() {
        xor_bit_length_pos_bit(circ, src, j, &bit);
        if add {
            ctrl_inc(circ, &bit, &acc[j..]);
        } else {
            ctrl_dec(circ, &bit, &acc[j..]);
        }
        xor_bit_length_pos_bit(circ, src, j, &bit);
    }
    circ.zero_and_free(bit);
}

fn add_const_mod_small(circ: &mut Circuit, acc: &[QReg], value: i64) {
    if acc.is_empty() {
        return;
    }
    use crate::point_add::trailmix_port::arith::ripple_add::add_const;
    let modulus = 1i128 << acc.len();
    let val = i128::from(value).rem_euclid(modulus) as u128;
    let bytes: Vec<u8> = (0..acc.len().div_ceil(8))
        .map(|i| (val >> (8 * i)) as u8)
        .collect();
    add_const(circ, acc, &bytes);
}

fn choose_dirty<'a>(ctrls: &[&QReg], target: &QReg, bank: &[&'a QReg]) -> Option<&'a QReg> {
    bank.iter()
        .copied()
        .find(|q| q.id() != target.id() && ctrls.iter().all(|c| c.id() != q.id()))
}

fn mcx_qubit_first(circ: &mut Circuit, ctrls: &[&QReg], target: &QReg, dirty_bank: &[&QReg]) {
    match ctrls {
        [] => circ.x(target),
        [c] => circ.cx(c, target),
        [a, b] => circ.ccx(a, b, target),
        _ => {
            if let Some(dirty) = choose_dirty(ctrls, target, dirty_bank) {
                crate::point_add::trailmix_port::arith::mcx::mcx_dirty(circ, ctrls, target, dirty);
            } else {
                crate::point_add::trailmix_port::arith::mcx::mcx_clean_k(circ, ctrls, target);
            }
        }
    }
}

fn inc_by_ctrls_no_clean(
    circ: &mut Circuit,
    prefix_ctrls: &[&QReg],
    acc: &[QReg],
    dirty_bank: &[&QReg],
) {
    for j in (0..acc.len()).rev() {
        let mut ctrls: Vec<&QReg> = Vec::with_capacity(prefix_ctrls.len() + j);
        ctrls.extend_from_slice(prefix_ctrls);
        ctrls.extend(acc[..j].iter());
        mcx_qubit_first(circ, &ctrls, &acc[j], dirty_bank);
    }
}

fn dec_by_ctrls_no_clean(
    circ: &mut Circuit,
    prefix_ctrls: &[&QReg],
    acc: &[QReg],
    dirty_bank: &[&QReg],
) {
    for q in acc {
        circ.x(q);
    }
    inc_by_ctrls_no_clean(circ, prefix_ctrls, acc, dirty_bank);
    for q in acc {
        circ.x(q);
    }
}

/// Qubit-first fold of the MSB position:
/// `pos = n - count_k(no source bit strictly above k)`.
///
/// This avoids the one-bit `clzm.posbit` scratch by applying the prefix
/// predicates directly to a small no-clean increment/decrement on `acc`.
fn fold_pos_by_prefix_count(circ: &mut Circuit, src: &[&QReg], acc: &[QReg], add: bool) {
    use crate::point_add::trailmix_port::arith::khattar_gidney::{
        kg_prefix_ancilla_count, KgPrefixAnd,
    };
    let n = src.len();
    if n == 0 || acc.is_empty() {
        return;
    }

    add_const_mod_small(circ, acc, if add { n as i64 } else { -(n as i64) });

    for q in src {
        circ.x(q);
    }
    let qbits: Vec<&QReg> = src.iter().rev().copied().collect();
    let nanc = kg_prefix_ancilla_count(n);
    let anc_owned = circ.alloc_qreg_bits("bll.count.kganc", nanc);
    let anc: Vec<&QReg> = anc_owned.iter().collect();
    let mut dirty_bank: Vec<&QReg> = Vec::with_capacity(src.len() + acc.len());
    dirty_bank.extend_from_slice(src);
    dirty_bank.extend(acc.iter());

    let kg = KgPrefixAnd::new(&qbits, &anc);
    let done = kg.forward(circ, |c, i, ctrls| {
        if i >= n {
            return;
        }
        if add {
            dec_by_ctrls_no_clean(c, ctrls, acc, &dirty_bank);
        } else {
            inc_by_ctrls_no_clean(c, ctrls, acc, &dirty_bank);
        }
    });
    done.reverse(circ, |_, _, _| {});

    drop(anc);
    for q in anc_owned {
        circ.zero_and_free(q);
    }
    for q in src {
        circ.x(q);
    }
}

/// `s += bitlen(src)` (or `-=` if dec). Qubit-first form: compute each bit of
/// the MSB index transiently and fold it into `s`, avoiding a full `bll.pos`
/// register at the bitlen peak.
fn bit_length_lean(circ: &mut Circuit, src: &[&QReg], s: &[QReg], dec: bool) {
    let n = src.len();
    if n == 0 {
        return;
    }
    let pbl = circ.push_section("p.bitlen");
    debug_assert!(
        (n as u64) <= (1u64 << (s.len() - 1)),
        "bit_length_lean: s width {} too small for n={n}",
        s.len()
    );
    if dec {
        for q in s {
            circ.x(q);
        }
    }
    let pfold = circ.push_section("p.countfold");
    fold_pos_bits_into_acc(circ, src, s, true); // s += MSB index
    circ.pop_section(&pfold);
    let pinc = circ.push_section("p.inc");
    crate::point_add::trailmix_port::arith::khattar_gidney::inc_khattar_gidney(circ, s);
    circ.pop_section(&pinc);
    if dec {
        for q in s {
            circ.x(q);
        }
    }
    circ.pop_section(&pbl);
}

fn lowq_clz_diff_const_fold_enabled() -> bool {
    if std::env::var("LOWQ_CLZ_DIFF_CONST_FOLD").ok().as_deref() != Some("1") {
        return false;
    }
    let target = std::env::var("TRAILMIX_Q_TARGET")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .expect("LOWQ_CLZ_DIFF_CONST_FOLD requires an integer TRAILMIX_Q_TARGET");
    assert!(
        matches!(target, 684 | 685),
        "LOWQ_CLZ_DIFF_CONST_FOLD is sealed to Q_TARGET 684/685"
    );
    true
}

fn lowq_hybrid_clz_enabled() -> bool {
    if std::env::var("LOWQ_HYBRID_CLZ").ok().as_deref() != Some("1") {
        return false;
    }
    assert_eq!(
        env_usize("TRAILMIX_SROT_W", 0),
        5,
        "LOWQ_HYBRID_CLZ requires the five-bit shift register"
    );
    assert_eq!(
        env_usize("TRAILMIX_THIN_CLZ_WINDOW", 0),
        78,
        "LOWQ_HYBRID_CLZ is sealed to the audited 78-bit windows"
    );
    assert!(
        matches!(env_usize("TRAILMIX_Q_TARGET", 0), 684 | 685),
        "LOWQ_HYBRID_CLZ is sealed to Q_TARGET 684/685"
    );
    true
}

fn lowq_exact_ctz_enabled() -> bool {
    if std::env::var("LOWQ_EXACT_CTZ").ok().as_deref() != Some("1") {
        return false;
    }
    assert_eq!(
        env_usize("TRAILMIX_SROT_W", 0),
        5,
        "LOWQ_EXACT_CTZ requires the five-bit shift register"
    );
    assert!(
        matches!(env_usize("TRAILMIX_Q_TARGET", 0), 684 | 685),
        "LOWQ_EXACT_CTZ is sealed to Q_TARGET 684/685"
    );
    true
}

/// One controlled fixed-distance shift layer. The forward direction is a
/// logical left shift on the promised branch because its top `distance` lanes
/// are zero. Reversing the pair order is the exact inverse for arbitrary data.
fn controlled_fixed_shift(
    circ: &mut Circuit,
    reg: &[QReg],
    control: &QReg,
    distance: usize,
    forward: bool,
) {
    if distance == 0 || distance >= reg.len() {
        return;
    }
    if forward {
        for hi in (distance..reg.len()).rev() {
            circ.cswap(control, &reg[hi], &reg[hi - distance]);
        }
    } else {
        for hi in distance..reg.len() {
            circ.cswap(control, &reg[hi], &reg[hi - distance]);
        }
    }
}

/// Toggle `out` iff the highest `prefix` lanes of `src` are all zero. The
/// peer register supplies restored dirty lenders and is unchanged.
fn toggle_zero_prefix_dirty(
    circ: &mut Circuit,
    src: &[QReg],
    prefix: usize,
    out: &QReg,
    peer: &[QReg],
) {
    use crate::point_add::trailmix_port::arith::mcx::mcx_dirty_ladder;

    assert!(prefix > 0 && prefix < src.len());
    let controls_owned = &src[src.len() - prefix..];
    for q in controls_owned {
        circ.x(q);
    }
    let controls: Vec<&QReg> = controls_owned.iter().collect();
    let dirty: Vec<&QReg> = peer.iter().take(prefix.saturating_sub(2)).collect();
    assert_eq!(
        dirty.len(),
        prefix.saturating_sub(2),
        "LOWQ_HYBRID_CLZ peer lender shortage"
    );
    mcx_dirty_ladder(circ, &controls, out, &dirty);
    for q in controls_owned.iter().rev() {
        circ.x(q);
    }
}

/// Toggle `out` iff `active` is set and the lowest `prefix` lanes of `src` are
/// all zero. Lenders may contain arbitrary data and are restored exactly.
fn toggle_active_zero_low_dirty(
    circ: &mut Circuit,
    src: &[QReg],
    prefix: usize,
    active: &QReg,
    out: &QReg,
    lenders: &[&QReg],
) {
    use crate::point_add::trailmix_port::arith::mcx::mcx_dirty_ladder;

    assert!(prefix > 0 && prefix < src.len());
    let controls_owned = &src[..prefix];
    for q in controls_owned {
        circ.x(q);
    }
    let mut controls: Vec<&QReg> = Vec::with_capacity(prefix + 1);
    controls.push(active);
    controls.extend(controls_owned.iter());
    let need = controls.len().saturating_sub(2);
    assert!(
        lenders.len() >= need,
        "LOWQ_EXACT_CTZ lender shortage: need={need} have={}",
        lenders.len()
    );
    mcx_dirty_ladder(circ, &controls, out, &lenders[..need]);
    for q in controls_owned.iter().rev() {
        circ.x(q);
    }
}

/// Compute `transcript = clz(src)` and normalize `src` to an MSB-one word.
/// Each branch bit controls one power-of-two shift and is retained until the
/// inverse restores `src`, so the map is bijective on the full basis space.
fn binary_clz_compute(
    circ: &mut Circuit,
    src: &[QReg],
    peer: &[QReg],
    transcript: &[QReg],
) {
    assert!(!src.is_empty() && src.len() <= (1usize << transcript.len()));
    for bit in (0..transcript.len()).rev() {
        let distance = 1usize << bit;
        if distance >= src.len() {
            continue;
        }
        toggle_zero_prefix_dirty(circ, src, distance, &transcript[bit], peer);
        controlled_fixed_shift(circ, src, &transcript[bit], distance, true);
    }
}

fn binary_clz_uncompute(
    circ: &mut Circuit,
    src: &[QReg],
    peer: &[QReg],
    transcript: &[QReg],
) {
    for bit in 0..transcript.len() {
        let distance = 1usize << bit;
        if distance >= src.len() {
            continue;
        }
        controlled_fixed_shift(circ, src, &transcript[bit], distance, false);
        toggle_zero_prefix_dirty(circ, src, distance, &transcript[bit], peer);
    }
}

/// PRE: `s=0`. Deposit `active*ctz(q)` directly into `s`, using `s` itself as
/// the branch transcript. The final left-shift sweep restores multi-hot q while
/// intentionally retaining s.
fn exact_multihot_ctz_deposit(
    circ: &mut Circuit,
    q: &[QReg],
    s: &[QReg],
    active: &QReg,
    lenders: &[&QReg],
) {
    assert_eq!(s.len(), 5, "LOWQ exact CTZ output width");
    let prev = circ.push_section("p.hctz.deposit");
    for bit in (0..s.len()).rev() {
        let distance = 1usize << bit;
        if distance >= q.len() {
            continue;
        }
        toggle_active_zero_low_dirty(circ, q, distance, active, &s[bit], lenders);
        controlled_fixed_shift(circ, q, &s[bit], distance, false);
    }
    for bit in 0..s.len() {
        let distance = 1usize << bit;
        if distance < q.len() {
            controlled_fixed_shift(circ, q, &s[bit], distance, true);
        }
    }
    circ.pop_section(&prev);
}

/// Exact gate inverse of `exact_multihot_ctz_deposit`.
/// PRE: `s=active*ctz(q)`. Restores q after the temporary normalization and
/// clears s to zero.
fn exact_multihot_ctz_erase(
    circ: &mut Circuit,
    q: &[QReg],
    s: &[QReg],
    active: &QReg,
    lenders: &[&QReg],
) {
    assert_eq!(s.len(), 5, "LOWQ exact CTZ output width");
    let prev = circ.push_section("p.hctz.erase");
    for bit in (0..s.len()).rev() {
        let distance = 1usize << bit;
        if distance < q.len() {
            controlled_fixed_shift(circ, q, &s[bit], distance, false);
        }
    }
    for bit in 0..s.len() {
        let distance = 1usize << bit;
        if distance >= q.len() {
            continue;
        }
        controlled_fixed_shift(circ, q, &s[bit], distance, true);
        toggle_active_zero_low_dirty(circ, q, distance, active, &s[bit], lenders);
    }
    circ.pop_section(&prev);
}

fn hybrid_transcript_width(max_window_len: usize) -> usize {
    let branch_bits = if max_window_len <= 1 {
        0
    } else {
        usize::BITS as usize - (max_window_len - 1).leading_zeros() as usize
    };
    branch_bits.max(5)
}

/// Deposit `active*(bitlen(a)-bitlen(b))` into the existing five-bit shift
/// register. Equal full register widths imply
/// `bitlen(a)-bitlen(b) = clz(b)-clz(a)` even when the audited low windows
/// differ. A single seven-bit transcript is reused sequentially.
fn hybrid_bitlen_diff_update(
    circ: &mut Circuit,
    a: &[QReg],
    b: &[QReg],
    lo_a: usize,
    lo_b: usize,
    target: &[QReg],
    active: &QReg,
    subtract_diff: bool,
) {
    use crate::point_add::trailmix_port::inversion::shrunken_pz_primitives::{
        ctrl_add, ctrl_sub,
    };

    assert_eq!(a.len(), b.len(), "LOWQ_HYBRID_CLZ requires equal full widths");
    assert_eq!(target.len(), 5, "LOWQ_HYBRID_CLZ target width");
    let prev = circ.push_section("p.hclz");
    let a_window = &a[lo_a.min(a.len() - 1)..];
    let b_window = &b[lo_b.min(b.len() - 1)..];
    let transcript = circ.alloc_qreg_bits(
        "hybrid.clz",
        hybrid_transcript_width(a_window.len().max(b_window.len())),
    );
    let target_refs: Vec<&QReg> = target.iter().collect();
    let low_refs: Vec<&QReg> = transcript[..target.len()].iter().collect();

    binary_clz_compute(circ, a_window, b, &transcript);
    if subtract_diff {
        ctrl_add(circ, active, &target_refs, &low_refs);
    } else {
        ctrl_sub(circ, active, &target_refs, &low_refs);
    }
    binary_clz_uncompute(circ, a_window, b, &transcript);

    binary_clz_compute(circ, b_window, a, &transcript);
    if subtract_diff {
        ctrl_sub(circ, active, &target_refs, &low_refs);
    } else {
        ctrl_add(circ, active, &target_refs, &low_refs);
    }
    binary_clz_uncompute(circ, b_window, a, &transcript);

    for q in transcript {
        circ.zero_and_free(q);
    }
    circ.pop_section(&prev);
}

fn hybrid_bitlen_diff_parity(
    circ: &mut Circuit,
    a: &[QReg],
    b: &[QReg],
    lo_a: usize,
    lo_b: usize,
    out: &QReg,
    active: &QReg,
) {
    assert_eq!(a.len(), b.len(), "LOWQ_HYBRID_CLZ requires equal full widths");
    let prev = circ.push_section("p.hclz.parity");
    let a_window = &a[lo_a.min(a.len() - 1)..];
    let b_window = &b[lo_b.min(b.len() - 1)..];
    let transcript = circ.alloc_qreg_bits(
        "hybrid.clz",
        hybrid_transcript_width(a_window.len().max(b_window.len())),
    );

    binary_clz_compute(circ, a_window, b, &transcript);
    circ.ccx(active, &transcript[0], out);
    binary_clz_uncompute(circ, a_window, b, &transcript);
    binary_clz_compute(circ, b_window, a, &transcript);
    circ.ccx(active, &transcript[0], out);
    binary_clz_uncompute(circ, b_window, a, &transcript);

    for q in transcript {
        circ.zero_and_free(q);
    }
    circ.pop_section(&prev);
}

/// `_middle` form of the clz-diff compute-USE-uncompute pattern: deposits the
/// `a` bitlen position into internal `pa`, folds the `b` position into it
/// transiently, FOLDS the diff
/// d = bitlen(a)-bitlen(b) (windowed) INTO `pa`, runs `body(circ, &pa)` with `pa`
/// holding the diff, then restores `pa` and un-deposits to |0>. No caller-supplied
/// diff register -- `pa` IS the diff, so nothing extra is live at the peak (this is
/// the `shrunken_pz_divide_forward` peak section). `w` sizes pa (must hold the window MSB
/// index and the signed diff). Scans un-nested (one KG ancilla set live at a time).
fn clz_diff_body_middle(
    circ: &mut Circuit,
    a: &[QReg],
    b: &[QReg],
    w: usize,
    lo_a: usize,
    lo_b: usize,
    body: impl FnOnce(&mut Circuit, &[QReg]),
) {
    use crate::point_add::trailmix_port::arith::ripple_add::add_const;
    let pbl = circ.push_section("p.bitlen");
    let aw: Vec<&QReg> = a[lo_a..a.len()].iter().collect();
    let bw: Vec<&QReg> = b[lo_b..b.len()].iter().collect();
    let pa = circ.alloc_qreg_bits("clzm.pa", w);
    let add_pa = |circ: &mut Circuit, pa: &[QReg], v: i64| {
        let val = i128::from(v).rem_euclid(1i128 << w) as u128;
        let bytes: Vec<u8> = (0..w.div_ceil(8)).map(|i| (val >> (8 * i)) as u8).collect();
        add_const(circ, pa, &bytes);
    };
    let na = aw.len();
    // Deposit pos_a. The b position is folded into pa below so no full pb register
    // is live at the bitlen peak.
    xor_const(circ, &pa, na);
    let pposa = circ.push_section("p.posa");
    bit_length_lean_middle(circ, &aw, &pa, |_| false); // pa = pos_a
    circ.pop_section(&pposa);

    let const_fold = lowq_clz_diff_const_fold_enabled();
    if const_fold {
        let pposb = circ.push_section("p.posb");
        fold_pos_bits_into_acc(circ, &bw, &pa, false); // pa -= pos_b
        circ.pop_section(&pposb);
        add_pa(circ, &pa, lo_a as i64 - lo_b as i64);
    } else {
        add_pa(circ, &pa, 1 + lo_a as i64);
        let pposb = circ.push_section("p.posb");
        fold_pos_bits_into_acc(circ, &bw, &pa, false); // pa -= pos_b
        circ.pop_section(&pposb);
        add_pa(circ, &pa, -(1 + lo_b as i64));
    }

    let pbody = circ.push_section("p.body");
    body(circ, &pa); // USE pa (= diff)
    circ.pop_section(&pbody);

    if const_fold {
        let pposb = circ.push_section("p.unposb");
        fold_pos_bits_into_acc(circ, &bw, &pa, true); // pa += pos_b
        circ.pop_section(&pposb);
        add_pa(circ, &pa, lo_b as i64 - lo_a as i64);
    } else {
        add_pa(circ, &pa, 1 + lo_b as i64);
        let pposb = circ.push_section("p.unposb");
        fold_pos_bits_into_acc(circ, &bw, &pa, true); // pa += pos_b
        circ.pop_section(&pposb);
        add_pa(circ, &pa, -(1 + lo_a as i64));
    }

    // un-deposit (self-inverse clean=false calls, reverse order).
    let pposa = circ.push_section("p.unposa");
    bit_length_lean_middle(circ, &aw, &pa, |_| false); // pa -> na
    circ.pop_section(&pposa);
    xor_const(circ, &pa, na); // pa -> 0
    for q in pa {
        circ.zero_and_free(q);
    }
    circ.pop_section(&pbl);
}

/// Rotate-LEFT `reg` in place by the quantum amount `s` (= reg << s, since the
/// aligned value's bitlen <= reg width so no nonzero bit wraps). Uses the ACYCLIC
/// `barrel_shift_inplace` (exactly `s.len()` layers, no wrap) rather than
/// `controlled_cyclic_rotate` (s.len()+1 full-width layers incl. a spurious
/// offset layer, + cyclic wrap churn): ~1.28x fewer cswaps. The no-wrap
/// precondition (top s bits of reg are |0>) is exactly the existing one.
/// forward=true is `<< s`; forward=false (restore) is `>> s`, Fredkin self-inverse.
fn rotate_left(circ: &mut Circuit, reg: &[QReg], s: &[QReg]) {
    crate::point_add::trailmix_port::arith::qshift_sub::barrel_shift_inplace(circ, reg, s, true);
}
fn rotate_right(circ: &mut Circuit, reg: &[QReg], s: &[QReg]) {
    crate::point_add::trailmix_port::arith::qshift_sub::barrel_shift_inplace(circ, reg, s, false);
}

/// `q[i] ^= active AND (s == i)` = `q ^= active·(1<<s)` -- the q-demux via KG
/// `unary_iterate_log_star` (~2 ccx/step) instead of a per-bit `eq_const_inplace` loop
/// (~58 tof/bit, ~30x more). active=0 => s masked to 0 => only i=0 gate fires,
/// `ANDed` with active=0 -> no-op. Self-inverse; `s` restored on exit.
fn set_bit_at_s_gated(circ: &mut Circuit, q_div: &[QReg], s: &[QReg], active: &QReg) {
    use crate::point_add::trailmix_port::arith::khattar_gidney::unary_iterate_log_star;
    let n_pad = q_div.len();
    if n_pad == 0 {
        return;
    }
    let prev = circ.push_section("p.demux");
    let sref: Vec<&QReg> = s.iter().collect();
    unary_iterate_log_star(circ, &sref, n_pad, |c, i, gate| {
        c.ccx(active, gate, &q_div[i]);
    });
    circ.pop_section(&prev);
}

fn add_noanc_refs(circ: &mut Circuit, a: &[&QReg], b: &[&QReg]) {
    crate::point_add::trailmix_port::arith::gidney_const_adder::hybrid_add_refs(circ, a, b, 0);
}

/// Unconditional `a -= b` (mod 2^len) via two's complement (X-bracket + add).
fn sub_refs(circ: &mut Circuit, a: &[&QReg], b: &[&QReg]) {
    for q in a {
        circ.x(q);
    }
    add_noanc_refs(circ, a, b);
    for q in a {
        circ.x(q);
    }
}

/// Controlled decrement `s -= 1` iff `g` (X-bracket + controlled increment).
fn ctrl_dec(circ: &mut Circuit, g: &QReg, s: &[QReg]) {
    use crate::point_add::trailmix_port::arith::khattar_gidney::cinc_khattar_gidney;
    for q in s {
        circ.x(q);
    }
    cinc_khattar_gidney(circ, s, g); // a=s, ctrl=g
    for q in s {
        circ.x(q);
    }
}

/// Controlled increment `s += 1` iff `g`.
fn ctrl_inc(circ: &mut Circuit, g: &QReg, s: &[QReg]) {
    use crate::point_add::trailmix_port::arith::khattar_gidney::cinc_khattar_gidney;
    cinc_khattar_gidney(circ, s, g);
}

/// Unconditional `a += b` (mod 2^len).
fn add_refs(circ: &mut Circuit, a: &[&QReg], b: &[&QReg]) {
    add_noanc_refs(circ, a, b);
}

/// Unpacked PZ state-machine registers. gcd pair (`a_gcd=A`, `b_gcd=B`) shrinks;
/// cofactor pair (ca=|a|, cb=|b|) grows. `q_div/q_mul` are the quotient pads
/// (~one quotient, ~26 bits each): `q_div` is built by the division (`q_div^=1`<<s),
/// swapped to `q_mul`, and DRAINED by the multiply (a += b<<`ctz(q_mul)`, clearing
/// it) -- the pipelined drain is what keeps the quotient record at one-quotient
/// size instead of a full ~256-bit tape. NOT removable (scripts/
/// `pz_fused_nopad_proto.py`: fusing gives the right inverse but s-recovery from
/// the cofactors mismatches ~30%, and an undrained pad accumulates a full tape).
pub struct PzSmRegs {
    pub a_gcd: Vec<QReg>,
    pub b_gcd: Vec<QReg>,
    pub ca: Vec<QReg>,
    pub cb: Vec<QReg>,
    pub q_div: Vec<QReg>,
    pub q_mul: Vec<QReg>,
}

/// Single-qubit state flags + sign. Invariant matches `pz_big_step`.
pub struct PzSmFlags {
    pub div_active: QReg,
    pub mul_active: QReg,
    pub offset: QReg,
    pub parity: QReg,
    pub sgn: QReg,
}

/// Load/unload the classical constant `c` into `reg` via X gates (self-inverse).
fn xor_const(circ: &mut Circuit, reg: &[QReg], c: usize) {
    for (j, q) in reg.iter().enumerate() {
        if (c >> j) & 1 == 1 {
            circ.x(q);
        }
    }
}

/// Magnitude compare `out ^= (a < b)` narrowed to the schedule window
/// `[lo, min(a.len, b.len))`. Used for the ALIGNED offset/o compares where a and
/// b share a bitlen (MSB guaranteed in [lo, hi) by the schedule), so the top bits
/// decide the order; a tie below `lo` (prob ~2^-(hi-lo) per the window width)
/// flips the result -- within the whole-pass tail tolerance. Forward and inverse
/// substeps call this with the same `lo`, so the (possibly-wrong) flag is
/// computed identically both ways and round-trips cleanly. Restores a,b.
/// NOT for the magnitude GATES (`g_mul/g_div)`: there A,B get arbitrarily close at
/// the div<->mul transition, so a deep tie is common, not a 2^-w tail.
fn narrow_lt(circ: &mut Circuit, a: &[QReg], b: &[QReg], out: &QReg, lo: usize) {
    let hi = a.len().min(b.len());
    let lo = lo.min(hi.saturating_sub(1));
    let ar: Vec<&QReg> = a[lo..hi].iter().collect();
    let br: Vec<&QReg> = b[lo..hi].iter().collect();
    borrow_compare_refs(circ, &ar, &br, out);
}

/// WINDOWED division substep: same as `division_substep_act` but the two clz
/// computations scan only the schedule's clz windows (`lo_a`/`lo_b` = window low
/// bounds for A/B) and the B<<s / restore rotates use `rot_bits` shift bits
/// (shift bound) instead of the full `s_rot` width. The offset-clean clz operates
/// on (A, `B_aligned`), both ~bitlen(A), so it reuses the A window (`lo_a`). For
/// in-schedule inputs this is gate-identical to `division_substep_act`.
#[allow(clippy::too_many_arguments)]
pub fn division_substep_windowed(
    circ: &mut Circuit,
    a: &[QReg],
    b: &[QReg],
    q_div: &[QReg],
    s_rot: &[QReg],
    offset: &QReg,
    gate_x: &[QReg],
    gate_y: &[QReg],
    counter: &[QReg],
    extra_lenders: &[&QReg],
    lo_a: usize,
    lo_b: usize,
    rot_bits: usize,
) {
    use crate::point_add::trailmix_port::inversion::shrunken_pz_primitives::ctrl_sub;
    let aref: Vec<&QReg> = a.iter().collect();
    let bref: Vec<&QReg> = b.iter().collect();
    let n_pad = q_div.len();
    let rb = rot_bits.min(s_rot.len());
    let w = s_rot.len();

    // diff = bitlen(A)-bitlen(B) (windowed _middle, folded into the clz's own pa);
    // mask s_rot = diff AND active.
    if lowq_hybrid_clz_enabled() {
        with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
            hybrid_bitlen_diff_update(circ, a, b, lo_a, lo_b, s_rot, active, false);
        });
    } else {
        clz_diff_body_middle(circ, a, b, w, lo_a, lo_b, |circ, diff| {
            with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
                for j in 0..w {
                    circ.ccx(active, &diff[j], &s_rot[j]);
                }
            });
        });
    }

    rotate_left(circ, b, &s_rot[0..rb]); // B <<= s if active (bounded rotator)

    // offset = active AND (A < B_aligned) -- narrowed (A,B_aligned share bitlen).
    {
        let or = circ.alloc_qreg("dg.offr");
        narrow_lt(circ, a, b, &or, lo_a);
        with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
            circ.ccx(active, &or, offset);
        });
        narrow_lt(circ, a, b, &or, lo_a);
        circ.zero_and_free(or);
    }
    rotate_right(circ, b, std::slice::from_ref(offset)); // B >>= 1 if offset
    ctrl_dec(circ, offset, s_rot); // s_rot -= 1 if offset => s_eff

    // clean offset via windowed _middle clz on (A, B_aligned) -> A window. The diff
    // lives in the clz's pa (this clz is the shrunken_pz_divide_forward peak section).
    if lowq_hybrid_clz_enabled() {
        with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
            hybrid_bitlen_diff_parity(circ, a, b, lo_a, lo_a, offset, active);
        });
    } else {
        clz_diff_body_middle(circ, a, b, w, lo_a, lo_a, |circ, diff| {
            with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
                circ.ccx(active, &diff[0], offset);
            });
        });
    }

    with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
        ctrl_sub(circ, active, &aref, &bref); // A -= B_aligned if active
    });

    with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
        set_bit_at_s_gated(circ, q_div, s_rot, active); // q_div ^= active·(1<<s_rot)
    });

    rotate_right(circ, b, &s_rot[0..rb]); // restore B >>= s_eff (bounded rotator)

    if lowq_exact_ctz_enabled() {
        let lenders: Vec<&QReg> = a
            .iter()
            .chain(b.iter())
            .chain(extra_lenders.iter().copied())
            .collect();
        with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
            exact_multihot_ctz_erase(circ, q_div, s_rot, active, &lenders);
        });
    } else {
        let t = circ.alloc_qreg_bits("dg.ctz", w);
        xor_const(circ, &t, n_pad);
        let rev: Vec<&QReg> = q_div.iter().rev().collect();
        bit_length_lean(circ, &rev, &t, true);
        let srr: Vec<&QReg> = s_rot.iter().collect();
        let tr: Vec<&QReg> = t.iter().collect();
        with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
            ctrl_sub(circ, active, &srr, &tr);
        });
        bit_length_lean(circ, &rev, &t, false);
        xor_const(circ, &t, n_pad);
        for lane in t {
            circ.zero_and_free(lane);
        }
    }
}

/// Gate-by-gate INVERSE of `division_substep_windowed` (for the backward pass).
/// Reverses the op sequence; the compute-use-uncompute blocks (clz-mask, offset,
/// offset-clean, q-demux) are self-inverse and run as-is; `rotate_left`<->right,
/// ctrl_sub->ctrl_add, ctrl_dec->ctrl_inc flip. Restores A += B<<`s_eff`, clears
/// the `q_div` bit, leaving `A/B/q_div/s/s_rot/offset` as before the forward step.
#[allow(clippy::too_many_arguments)]
pub fn division_substep_windowed_inv(
    circ: &mut Circuit,
    a: &[QReg],
    b: &[QReg],
    q_div: &[QReg],
    s_rot: &[QReg],
    offset: &QReg,
    gate_x: &[QReg],
    gate_y: &[QReg],
    counter: &[QReg],
    extra_lenders: &[&QReg],
    lo_a: usize,
    lo_b: usize,
    rot_bits: usize,
) {
    use crate::point_add::trailmix_port::inversion::shrunken_pz_primitives::ctrl_add;
    let aref: Vec<&QReg> = a.iter().collect();
    let bref: Vec<&QReg> = b.iter().collect();
    let n_pad = q_div.len();
    let rb = rot_bits.min(s_rot.len());
    let w = s_rot.len();

    // 12' reconstruct s_rot from the least-significant set quotient bit.
    if lowq_exact_ctz_enabled() {
        let lenders: Vec<&QReg> = a
            .iter()
            .chain(b.iter())
            .chain(extra_lenders.iter().copied())
            .collect();
        with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
            exact_multihot_ctz_deposit(circ, q_div, s_rot, active, &lenders);
        });
    } else {
        let t = circ.alloc_qreg_bits("dg.ctz", w);
        xor_const(circ, &t, n_pad);
        let rev: Vec<&QReg> = q_div.iter().rev().collect();
        bit_length_lean(circ, &rev, &t, true);
        let srr: Vec<&QReg> = s_rot.iter().collect();
        let tr: Vec<&QReg> = t.iter().collect();
        with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
            ctrl_add(circ, active, &srr, &tr);
        });
        bit_length_lean(circ, &rev, &t, false);
        xor_const(circ, &t, n_pad);
        for lane in t {
            circ.zero_and_free(lane);
        }
    }
    // 11' rotate_left (was rotate_right restore).
    rotate_left(circ, b, &s_rot[0..rb]);
    // 10' q_div demux (self-inverse XOR).
    with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
        set_bit_at_s_gated(circ, q_div, s_rot, active); // q_div ^= active·(1<<s_rot)
    });
    // 9' ctrl_sub -> ctrl_add (restore A += B_aligned).
    with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
        ctrl_add(circ, active, &aref, &bref);
    });
    // 8' offset clean (self-inverse, _middle); diff in the clz's pa.
    if lowq_hybrid_clz_enabled() {
        with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
            hybrid_bitlen_diff_parity(circ, a, b, lo_a, lo_a, offset, active);
        });
    } else {
        clz_diff_body_middle(circ, a, b, w, lo_a, lo_a, |circ, diff| {
            with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
                circ.ccx(active, &diff[0], offset);
            });
        });
    }
    // 7' ctrl_dec -> ctrl_inc.
    ctrl_inc(circ, offset, s_rot);
    // 6' rotate_left (was rotate_right by offset).
    rotate_left(circ, b, std::slice::from_ref(offset));
    // 5' offset compute (self-inverse) -- narrowed, same window as forward.
    {
        let or = circ.alloc_qreg("dg.offr");
        narrow_lt(circ, a, b, &or, lo_a);
        with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
            circ.ccx(active, &or, offset);
        });
        narrow_lt(circ, a, b, &or, lo_a);
        circ.zero_and_free(or);
    }
    // 4' rotate_right (was rotate_left B<<s).
    rotate_right(circ, b, &s_rot[0..rb]);
    // 3',2',1' clz-mask block (self-inverse, _middle) -- clears s_rot to |0>.
    if lowq_hybrid_clz_enabled() {
        with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
            hybrid_bitlen_diff_update(circ, a, b, lo_a, lo_b, s_rot, active, true);
        });
    } else {
        clz_diff_body_middle(circ, a, b, w, lo_a, lo_b, |circ, diff| {
            with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
                for j in 0..w {
                    circ.ccx(active, &diff[j], &s_rot[j]);
                }
            });
        });
    }
}

/// `out ^= (reg != 0)` (restores reg).
fn or_nonzero(circ: &mut Circuit, reg: &[QReg], out: &QReg) {
    use crate::point_add::trailmix_port::arith::mcx::mcx_clean_k;
    let prev = circ.push_section("p.ornz");
    for q in reg {
        circ.x(q);
    }
    let refs: Vec<&QReg> = reg.iter().collect();
    mcx_clean_k(circ, &refs, out); // out ^= (reg == 0)
    for q in reg {
        circ.x(q);
    }
    circ.x(out); // out ^= (reg != 0)
    circ.pop_section(&prev);
}

/// `out ^= (reg == 0)` via X-bracket + mcx (clean, self-inverse, restores reg).
fn or_is_zero(circ: &mut Circuit, reg: &[QReg], out: &QReg) {
    use crate::point_add::trailmix_port::arith::mcx::mcx_clean_k;
    let prev = circ.push_section("p.orz");
    for q in reg {
        circ.x(q);
    }
    let refs: Vec<&QReg> = reg.iter().collect();
    mcx_clean_k(circ, &refs, out); // out ^= (reg == 0)
    for q in reg {
        circ.x(q);
    }
    circ.pop_section(&prev);
}

fn or_is_zero_dirty(circ: &mut Circuit, reg: &[QReg], out: &QReg, dirty: &QReg) {
    use crate::point_add::trailmix_port::arith::mcx::mcx_dirty_any_k;
    debug_assert!(dirty.id() != out.id());
    debug_assert!(reg.iter().all(|q| q.id() != dirty.id()));
    let prev = circ.push_section("p.orz");
    for q in reg {
        circ.x(q);
    }
    let refs: Vec<&QReg> = reg.iter().collect();
    mcx_dirty_any_k(circ, &refs, out, dirty); // out ^= (reg == 0)
    for q in reg {
        circ.x(q);
    }
    circ.pop_section(&prev);
}

/// WINDOWED multiply substep: same as `multiply_substep_act` but the two clz
/// computations scan the schedule's cofactor clz windows. The `o` clz is on
/// (ca, cb<<s2), both ~bitlen(ca) -> ca window (`ca_window`). The s_rot-clean clz is
/// on (cb, ca) -> cb/ca windows. The cb<<s2 / restore rotates use `rot_bits`.
/// q (ctz) is small -> not windowed. Gate-identical for in-schedule inputs.
#[allow(clippy::too_many_arguments)]
pub fn multiply_substep_windowed(
    circ: &mut Circuit,
    a: &[QReg],
    b: &[QReg],
    q_mul: &[QReg],
    s_rot: &[QReg],
    off: &QReg,
    gate_x: &[QReg],
    gate_y: &[QReg],
    counter: &[QReg],
    extra_lenders: &[&QReg],
    ca_window: usize,
    cb_window: usize,
    rot_bits: usize,
) {
    use crate::point_add::trailmix_port::inversion::shrunken_pz_primitives::ctrl_add;
    let aref: Vec<&QReg> = a.iter().collect();
    let bref: Vec<&QReg> = b.iter().collect();
    let n_pad = q_mul.len();
    let rb = rot_bits.min(s_rot.len());
    let w = s_rot.len();

    if lowq_exact_ctz_enabled() {
        let lenders: Vec<&QReg> = a
            .iter()
            .chain(b.iter())
            .chain(extra_lenders.iter().copied())
            .collect();
        with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
            exact_multihot_ctz_deposit(circ, q_mul, s_rot, active, &lenders);
        });
    } else {
        let t = circ.alloc_qreg_bits("mg.ctz", w);
        let rev: Vec<&QReg> = q_mul.iter().rev().collect();
        xor_const(circ, &t, n_pad);
        bit_length_lean(circ, &rev, &t, true);
        with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
            for j in 0..w {
                circ.ccx(active, &t[j], &s_rot[j]);
            }
        });
        bit_length_lean(circ, &rev, &t, false);
        xor_const(circ, &t, n_pad);
        for lane in t {
            circ.zero_and_free(lane);
        }
    }

    with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
        set_bit_at_s_gated(circ, q_mul, s_rot, active); // q_mul ^= active·(1<<s_rot)
    });

    rotate_left(circ, b, &s_rot[0..rb]); // b <<= s if active (bounded rotator)
    with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
        ctrl_add(circ, active, &aref, &bref); // a += b<<s if active
    });

    // o = active AND (bitlen(ca) != bitlen(cb<<s2)) -- ca window, _middle; diff in
    // the clz's pa. This clz is the shrunken_pz_divide_forward peak section.
    if lowq_hybrid_clz_enabled() {
        with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
            hybrid_bitlen_diff_parity(circ, a, b, ca_window, ca_window, off, active);
        });
    } else {
        clz_diff_body_middle(circ, a, b, w, ca_window, ca_window, |circ, diff| {
            with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
                circ.ccx(active, &diff[0], off);
            });
        });
    }
    rotate_left(circ, b, std::slice::from_ref(off)); // b <<= 1 if o
    ctrl_inc(circ, off, s_rot);
    {
        let lt = circ.alloc_qreg("mg.cleanlt");
        narrow_lt(circ, a, b, &lt, ca_window);
        with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
            circ.ccx(active, &lt, off);
        });
        narrow_lt(circ, a, b, &lt, ca_window);
        circ.zero_and_free(lt);
    }
    rotate_right(circ, b, &s_rot[0..rb]); // restore b >>= s_eff (bounded rotator)

    // clean s_rot via _middle clz on (cb, ca): s_rot += (bitlen(cb)-bitlen(ca)).
    if lowq_hybrid_clz_enabled() {
        with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
            hybrid_bitlen_diff_update(circ, b, a, cb_window, ca_window, s_rot, active, false);
        });
    } else {
        clz_diff_body_middle(circ, b, a, w, cb_window, ca_window, |circ, diff| {
            let srr: Vec<&QReg> = s_rot.iter().collect();
            let ter: Vec<&QReg> = diff.iter().collect();
            with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
                ctrl_add(circ, active, &srr, &ter);
            });
        });
    }
}

/// Gate-by-gate INVERSE of `multiply_substep_windowed` (backward pass). Reverses
/// the sequence; clz/o/q-demux blocks are self-inverse; `rotate_left`<->right,
/// ctrl_add->ctrl_sub, ctrl_inc->ctrl_dec flip. Restores ca -= cb<<s2, re-sets
/// the `q_mul` bit.
#[allow(clippy::too_many_arguments)]
pub fn multiply_substep_windowed_inv(
    circ: &mut Circuit,
    a: &[QReg],
    b: &[QReg],
    q_mul: &[QReg],
    s_rot: &[QReg],
    off: &QReg,
    gate_x: &[QReg],
    gate_y: &[QReg],
    counter: &[QReg],
    extra_lenders: &[&QReg],
    ca_window: usize,
    cb_window: usize,
    rot_bits: usize,
) {
    use crate::point_add::trailmix_port::inversion::shrunken_pz_primitives::ctrl_sub;
    let aref: Vec<&QReg> = a.iter().collect();
    let bref: Vec<&QReg> = b.iter().collect();
    let n_pad = q_mul.len();
    let rb = rot_bits.min(s_rot.len());
    let w = s_rot.len();

    // 10' s_rot clean inverse: ctrl_add -> ctrl_sub (_middle); diff in the clz's pa.
    if lowq_hybrid_clz_enabled() {
        with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
            hybrid_bitlen_diff_update(circ, b, a, cb_window, ca_window, s_rot, active, true);
        });
    } else {
        clz_diff_body_middle(circ, b, a, w, cb_window, ca_window, |circ, diff| {
            let srr: Vec<&QReg> = s_rot.iter().collect();
            let ter: Vec<&QReg> = diff.iter().collect();
            with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
                ctrl_sub(circ, active, &srr, &ter);
            });
        });
    }
    // 9' rotate_left (was rotate_right restore).
    rotate_left(circ, b, &s_rot[0..rb]);
    // 8' clean-o block (self-inverse) -- narrowed, same window as forward.
    {
        let lt = circ.alloc_qreg("mg.cleanlt");
        narrow_lt(circ, a, b, &lt, ca_window);
        with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
            circ.ccx(active, &lt, off);
        });
        narrow_lt(circ, a, b, &lt, ca_window);
        circ.zero_and_free(lt);
    }
    // 7' ctrl_inc -> ctrl_dec.
    ctrl_dec(circ, off, s_rot);
    // 6' rotate_right (was rotate_left by o).
    rotate_right(circ, b, std::slice::from_ref(off));
    // 5' o clz block (self-inverse, _middle); diff in the clz's pa.
    if lowq_hybrid_clz_enabled() {
        with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
            hybrid_bitlen_diff_parity(circ, a, b, ca_window, ca_window, off, active);
        });
    } else {
        clz_diff_body_middle(circ, a, b, w, ca_window, ca_window, |circ, diff| {
            with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
                circ.ccx(active, &diff[0], off);
            });
        });
    }
    // 4' ctrl_add -> ctrl_sub (undo ca += cb<<s2).
    with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
        ctrl_sub(circ, active, &aref, &bref);
    });
    // 3' rotate_right (was rotate_left cb<<s2).
    rotate_right(circ, b, &s_rot[0..rb]);
    // 2' q_mul clear demux (self-inverse).
    with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
        set_bit_at_s_gated(circ, q_mul, s_rot, active); // q_mul ^= active·(1<<s_rot)
    });
    // 1' clear the least-significant-set-bit index from s_rot.
    if lowq_exact_ctz_enabled() {
        let lenders: Vec<&QReg> = a
            .iter()
            .chain(b.iter())
            .chain(extra_lenders.iter().copied())
            .collect();
        with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
            exact_multihot_ctz_erase(circ, q_mul, s_rot, active, &lenders);
        });
    } else {
        let t = circ.alloc_qreg_bits("mg.ctz", w);
        let rev: Vec<&QReg> = q_mul.iter().rev().collect();
        xor_const(circ, &t, n_pad);
        bit_length_lean(circ, &rev, &t, true);
        with_recomputed_counter_active_gate(circ, gate_x, gate_y, counter, |circ, active| {
            for j in 0..w {
                circ.ccx(active, &t[j], &s_rot[j]);
            }
        });
        bit_length_lean(circ, &rev, &t, false);
        xor_const(circ, &t, n_pad);
        for lane in t {
            circ.zero_and_free(lane);
        }
    }
}

// NEXT (reversible_pz_notes.md has the primitive mapping):
//   fn normalize_input(circ, x, sgn)               -- x -> min(x,P-x), set sgn
//   fn division_substep(circ, regs, flags, s, bound)
//   fn multiply_substep(circ, regs, flags, s, bound)
//   fn transition(circ, regs, flags)
//   fn iterate(circ, regs, flags, n_iters)         -- the fixed-count driver
//   fn recover_inverse(circ, regs, flags)          -- parity^sgn sign fix
//   test pz_sm_faithful  -- per-iter contract vs a Rust port of pz_big_step

// ===== shrunken_pz reversible inversion step driver (shared fwd/back, used by
// the round-trip test AND the EC-add) =====

// ---- shared forward/backward step helpers (used by the round-trip) ----

/// Like calling `gate_and_active` twice around `body`, but HOLDS the comparator
/// flag `lt=(x<y)` across the substep (which leaves x,y stationary) so the
/// full-width `borrow_compare` runs 2x not 4x. g = (x<y) AND active during body.
pub(crate) fn gate_hold(
    c: &mut Circuit,
    x: &[QReg],
    y: &[QReg],
    active: &QReg,
    g: &QReg,
    body: impl FnOnce(&mut Circuit, &QReg),
) {
    let lt = c.alloc_qreg("gh.lt");
    let xr: Vec<&QReg> = x.iter().collect();
    let yr: Vec<&QReg> = y.iter().collect();
    borrow_compare_refs(c, &xr, &yr, &lt); // lt = (x<y)
    c.ccx(&lt, active, g); // g = lt AND active
    body(c, g);
    c.ccx(&lt, active, g); // uncompute g
    borrow_compare_refs(c, &xr, &yr, &lt); // uncompute lt
    c.zero_and_free(lt);
}

/// Counter-gated variant of [`gate_hold`] that avoids keeping a separate
/// `active=(counter==0)` qubit live across the body. It computes
/// `g=(x<y)&(counter==0)`, runs `body`, then recomputes the same predicate to
/// clear `g`.
pub(crate) fn gate_hold_counter_active(
    c: &mut Circuit,
    x: &[QReg],
    y: &[QReg],
    counter: &[QReg],
    g: &QReg,
    body: impl FnOnce(&mut Circuit, &QReg),
) {
    let lt = c.alloc_qreg("gh.lt");
    let xr: Vec<&QReg> = x.iter().collect();
    let yr: Vec<&QReg> = y.iter().collect();
    borrow_compare_refs(c, &xr, &yr, &lt); // lt = (x<y)
    let toggle_gate = |c: &mut Circuit, lt: &QReg, g: &QReg| {
        if counter.is_empty() {
            c.cx(lt, g);
        } else {
            let act = c.alloc_qreg("gh.act");
            or_is_zero_dirty(c, counter, &act, lt);
            c.ccx(lt, &act, g);
            or_is_zero_dirty(c, counter, &act, lt);
            c.zero_and_free(act);
        }
    };
    toggle_gate(c, &lt, g);
    body(c, g);
    toggle_gate(c, &lt, g);
    borrow_compare_refs(c, &xr, &yr, &lt); // uncompute lt
    c.zero_and_free(lt);
}

fn with_counter_active_gate(
    c: &mut Circuit,
    lt: &QReg,
    counter: &[QReg],
    body: impl FnOnce(&mut Circuit, &QReg),
) {
    if counter.is_empty() {
        body(c, lt);
        return;
    }
    let g = c.alloc_qreg("gh.g");
    let toggle = |c: &mut Circuit, g: &QReg| {
        let act = c.alloc_qreg("gh.act");
        or_is_zero_dirty(c, counter, &act, lt);
        c.ccx(lt, &act, g);
        or_is_zero_dirty(c, counter, &act, lt);
        c.zero_and_free(act);
    };
    toggle(c, &g);
    body(c, &g);
    toggle(c, &g);
    c.zero_and_free(g);
}

fn with_recomputed_counter_active_gate(
    c: &mut Circuit,
    x: &[QReg],
    y: &[QReg],
    counter: &[QReg],
    body: impl FnOnce(&mut Circuit, &QReg),
) {
    let lt = c.alloc_qreg("gh.lt.re");
    let xr: Vec<&QReg> = x.iter().collect();
    let yr: Vec<&QReg> = y.iter().collect();
    borrow_compare_refs(c, &xr, &yr, &lt);
    with_counter_active_gate(c, &lt, counter, body);
    borrow_compare_refs(c, &xr, &yr, &lt);
    c.zero_and_free(lt);
}

/// Like [`gate_hold_counter_active`], but keeps only `lt=(x<y)` live across the
/// body. The body materializes `lt&(counter==0)` around the controlled segments
/// with [`with_counter_active_gate`], so the gate qubit is not live during
/// bit-length scans.
pub(crate) fn gate_hold_lazy_counter_active(
    c: &mut Circuit,
    x: &[QReg],
    y: &[QReg],
    counter: &[QReg],
    body: impl FnOnce(&mut Circuit, &QReg, &[QReg]),
) {
    let lt = c.alloc_qreg("gh.lt");
    let xr: Vec<&QReg> = x.iter().collect();
    let yr: Vec<&QReg> = y.iter().collect();
    borrow_compare_refs(c, &xr, &yr, &lt); // lt = (x<y)
    body(c, &lt, counter);
    borrow_compare_refs(c, &xr, &yr, &lt); // uncompute lt
    c.zero_and_free(lt);
}

/// done-counter (forward: counter += conv) / its inverse (counter -= conv),
/// conv = (A==0 & q==0). `done` is clean scratch (|0> at exit). User's recipe.
pub(crate) fn done_counter_fn(
    c: &mut Circuit,
    aa: &[QReg],
    qq: &[QReg],
    counter: &[QReg],
    inverse: bool,
) {
    if counter.is_empty() {
        return;
    }
    let done = c.alloc_qreg("done");
    let conv = |c: &mut Circuit, done: &QReg| {
        let az = c.alloc_qreg("d.az");
        let qz = c.alloc_qreg("d.qz");
        or_is_zero(c, aa, &az);
        or_is_zero(c, qq, &qz);
        c.ccx(&az, &qz, done); // done ^= (A==0 & q==0)
        or_is_zero(c, qq, &qz);
        or_is_zero(c, aa, &az);
        c.zero_and_free(qz);
        c.zero_and_free(az);
    };
    let cnz = |c: &mut Circuit, done: &QReg| {
        let z = c.alloc_qreg("d.cnz");
        or_nonzero(c, counter, &z);
        c.cx(&z, done); // done ^= (counter != 0)
        or_nonzero(c, counter, &z);
        c.zero_and_free(z);
    };
    if inverse {
        cnz(c, &done);
        ctrl_dec(c, &done, counter);
        conv(c, &done);
    } else {
        conv(c, &done);
        ctrl_inc(c, &done, counter);
        cnz(c, &done);
    }
    c.zero_and_free(done);
}

/// One forward (inverse=false) or backward (inverse=true) `shrunken_pz` step on the
/// dynamic-W registers at their current width. Resize is done by the caller.
#[allow(clippy::too_many_arguments)]
pub(crate) fn shrunken_pz_pass_step(
    c: &mut Circuit,
    aa: &[QReg],
    bb: &[QReg],
    cca: &[QReg],
    ccb: &[QReg],
    qq: &[QReg],
    counter: &[QReg],
    parity: &QReg,
    s_rot: &[QReg],
    off: &QReg,
    i: usize,
    inverse: bool,
) {
    use crate::point_add::trailmix_port::inversion::shrunken_pz_schedule::{reg_los, shift_bounds};
    fn rb(b: usize) -> usize {
        if b == 0 {
            1
        } else {
            64 - (b as u64).leading_zeros() as usize
        }
    }
    let (lo_a, lo_b, ca_window, cb_window, _) = reg_los(i);
    let (sdb, s2b) = shift_bounds(i);
    // Swap, gated g_swap=(q==0 & A!=0 & active). HOLD the (q==0)/(A!=0) flags
    // across the cswaps so or_nonzero(A)/or_is_zero(q) run 2x not 4x per step
    // (the swap preserves both predicates: q untouched, A_new=B_old!=0).
    let swap = |c: &mut Circuit, active: &QReg| {
        let qz = c.alloc_qreg("sw.qz");
        let anz = c.alloc_qreg("sw.anz");
        let t = c.alloc_qreg("sw.t");
        let g = c.alloc_qreg("g_swap");
        or_is_zero(c, qq, &qz);
        or_nonzero(c, aa, &anz);
        c.ccx(&qz, &anz, &t); // t = (q==0 & A!=0)
        c.ccx(&t, active, &g); // g_swap = t AND active
        for j in 0..aa.len() {
            c.cswap(&g, &aa[j], &bb[j]);
        }
        for j in 0..cca.len() {
            c.cswap(&g, &cca[j], &ccb[j]);
        }
        c.cx(&g, parity);
        c.ccx(&t, active, &g); // uncompute g (t,active preserved)
        c.ccx(&qz, &anz, &t); // uncompute t (qz held; anz=A_old!=0)
        or_nonzero(c, aa, &anz); // post-swap A=B_old!=0 -> clears anz
        or_is_zero(c, qq, &qz);
        c.zero_and_free(g);
        c.zero_and_free(t);
        c.zero_and_free(anz);
        c.zero_and_free(qz);
    };
    if inverse {
        done_counter_fn(c, aa, qq, counter, true);
        let active = compute_active(c, counter);
        swap(c, &active); // self-inverse
        uncompute_active(c, counter, &active);
        c.zero_and_free(active);
        {
            let lenders: Vec<&QReg> = cca.iter().chain(ccb.iter()).collect();
            division_substep_windowed_inv(
                c, aa, bb, qq, s_rot, off, cca, ccb, counter, &lenders, lo_a, lo_b, rb(sdb),
            );
        }
        {
            let lenders: Vec<&QReg> = aa.iter().chain(bb.iter()).collect();
            multiply_substep_windowed_inv(
                c,
                cca,
                ccb,
                qq,
                s_rot,
                off,
                aa,
                bb,
                counter,
                &lenders,
                ca_window,
                cb_window,
                rb(s2b),
            );
        }
    } else {
        {
            let lenders: Vec<&QReg> = aa.iter().chain(bb.iter()).collect();
            multiply_substep_windowed(
                c,
                cca,
                ccb,
                qq,
                s_rot,
                off,
                aa,
                bb,
                counter,
                &lenders,
                ca_window,
                cb_window,
                rb(s2b),
            );
        }
        {
            let lenders: Vec<&QReg> = cca.iter().chain(ccb.iter()).collect();
            division_substep_windowed(
                c, aa, bb, qq, s_rot, off, cca, ccb, counter, &lenders, lo_a, lo_b, rb(sdb),
            );
        }
        let active = compute_active(c, counter);
        swap(c, &active);
        uncompute_active(c, counter, &active);
        c.zero_and_free(active);
        done_counter_fn(c, aa, qq, counter, false);
    }
}

/// Resize a dynamic-W register to `target` bits: free high qubits (must be |0>)
/// or alloc fresh |0> ones, in place.
pub(crate) fn shrunken_pz_resize(c: &mut Circuit, reg: &mut Vec<QReg>, target: usize, name: &str) {
    while reg.len() > target {
        let q = reg.pop().unwrap();
        c.zero_and_free(q);
    }
    while reg.len() < target {
        let k = reg.len();
        reg.push(c.alloc_qreg(&format!("{name}[{k}]")));
    }
}

/// FORWARD `shrunken_pz` inversion driver. PRE: the registers hold the `S_0` state at width
/// `reg_widths(0)` -- A=p, B=|x| (sign-adjusted, < p/2), ca=0, cb=1, q=0,
/// counter=0, parity=1. Runs all `SHRUNKEN_PZ_NSTEPS` forward steps (resizing per step),
/// leaving the modular inverse of |x| in `ccb` (up to the `parity` bit: the true
/// value is `parity ? cb : p-cb`), with A=p, B=|x| at the EEA terminal. `s`,
/// `s_rot` (9 bits each), `off`, `parity`, `counter` (10 bits) are fixed-width.
#[allow(clippy::too_many_arguments)]
pub(crate) fn shrunken_pz_invert_forward(
    c: &mut Circuit,
    aa: &mut Vec<QReg>,
    bb: &mut Vec<QReg>,
    cca: &mut Vec<QReg>,
    ccb: &mut Vec<QReg>,
    qq: &mut Vec<QReg>,
    counter: &[QReg],
    parity: &QReg,
    s_rot: &[QReg],
    off: &QReg,
) {
    use crate::point_add::trailmix_port::inversion::shrunken_pz_schedule::{reg_widths, SHRUNKEN_PZ_NSTEPS};
    for i in 0..SHRUNKEN_PZ_NSTEPS {
        let (wa, wb, wca, wcb, wq) = reg_widths(i);
        let wab = trailmix_ab_width(wa.max(wb));
        let wcacb = trailmix_cacb_width(wca.max(wcb));
        shrunken_pz_resize(c, aa, wab, "A");
        shrunken_pz_resize(c, bb, wab, "B");
        shrunken_pz_resize(c, cca, wcacb, "ca");
        shrunken_pz_resize(c, ccb, wcacb, "cb");
        shrunken_pz_resize(c, qq, trailmix_q_width_step(wq, wa, wb, wca, wcb), "q");
        shrunken_pz_pass_step(
            c, aa, bb, cca, ccb, qq, counter, parity, s_rot, off, i, false,
        );
    }
}

/// BACKWARD `shrunken_pz` inversion driver (gate-for-gate inverse of `shrunken_pz_invert_forward`).
/// Restores the `S_0` state (A=p, B=|x|, ca=0, cb=1, q=0, counter=0, parity=1) and
/// uncomputes the inverse from `ccb`. Resizes back down per step.
#[allow(clippy::too_many_arguments)]
pub(crate) fn shrunken_pz_invert_backward(
    c: &mut Circuit,
    aa: &mut Vec<QReg>,
    bb: &mut Vec<QReg>,
    cca: &mut Vec<QReg>,
    ccb: &mut Vec<QReg>,
    qq: &mut Vec<QReg>,
    counter: &[QReg],
    parity: &QReg,
    s_rot: &[QReg],
    off: &QReg,
) {
    use crate::point_add::trailmix_port::inversion::shrunken_pz_schedule::{reg_widths, SHRUNKEN_PZ_NSTEPS};
    for i in (0..SHRUNKEN_PZ_NSTEPS).rev() {
        shrunken_pz_pass_step(
            c, aa, bb, cca, ccb, qq, counter, parity, s_rot, off, i, true,
        );
        if i > 0 {
            let (wa, wb, wca, wcb, wq) = reg_widths(i - 1);
            let wab = trailmix_ab_width(wa.max(wb));
            let wcacb = trailmix_cacb_width(wca.max(wcb));
            shrunken_pz_resize(c, aa, wab, "A");
            shrunken_pz_resize(c, bb, wab, "B");
            shrunken_pz_resize(c, cca, wcacb, "ca");
            shrunken_pz_resize(c, ccb, wcacb, "cb");
            shrunken_pz_resize(c, qq, trailmix_q_width_step(wq, wa, wb, wca, wcb), "q");
        }
    }
}

/// `lambda = dy / dx mod p`, with `dx` and `dy` PRESERVED. `dx`, `dy` are 257-bit
/// registers holding field elements in [0, p). Returns `(dx, dy, lambda)` -- dx
/// and dy unchanged (dy reconstructed via the HMR-ghost trick), lambda = dy·dx^-1
/// (257 bits, canonical). This is the shrunken_pz-native EC slope: the EEA consumes dx
/// (restored by the reverse), and dy is GHOSTED during the reverse so dy and
/// lambda are never both live across the inversion -> peak ~ EEA-peak + 256.
pub fn shrunken_pz_divide_forward(
    c: &mut Circuit,
    mut dx: Vec<QReg>,
    dy: Vec<QReg>,
) -> (Vec<QReg>, Vec<QReg>, Vec<QReg>) {
    use crate::point_add::trailmix_port::arith::compare::compare_geq_const;
    use crate::point_add::trailmix_port::arith::rfold_mbu::mod_mul_rfold_mbu;
    use crate::point_add::trailmix_port::inversion::shrunken_pz_schedule::reg_widths;
    use crate::point_add::trailmix_port::num_bigint::BigUint;
    assert_eq!(dx.len(), 257);
    assert_eq!(dy.len(), 257);
    // sgn = dx > p/2  <=>  dx >= (p+1)/2.
    let half_bytes = vec![
        0x18, 0xfe, 0xff, 0x7f, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f, 0x00,
    ];
    let p_bytes = crate::point_add::trailmix_port::mod_arith::SECP256K1_P_LE;

    // --- sign-adjust dx -> |dx| < p/2 (the schedule assumes |x| < p/2) ---
    let sgn = c.alloc_qreg("shpzdiv.sgn");
    compare_geq_const(c, &dx, &half_bytes, &sgn);
    controlled_field_neg(c, &sgn, &dx); // dx := (sgn ? p-dx : dx) = |dx|

    // --- set up the inversion S_0 state (B = |dx|, A = p, cb = 1, parity = 1) ---
    let (a0, b0, ca0, cb0, q0) = reg_widths(0);
    let (wg0, wc0) = (a0.max(b0), ca0.max(cb0));
    shrunken_pz_resize(c, &mut dx, wg0, "B"); // |dx| becomes the EEA B register
    let mut aa = c.alloc_qreg_bits("shpzdiv.A", wg0);
    let mut cca = c.alloc_qreg_bits("shpzdiv.ca", wc0);
    let mut ccb = c.alloc_qreg_bits("shpzdiv.cb", wc0);
    let mut qq = c.alloc_qreg_bits("shpzdiv.q", q0.max(1));
    let s_rot = c.alloc_qreg_bits("shpzdiv.srot", trailmix_srot_width());
    let off = c.alloc_qreg("shpzdiv.off");
    let parity = c.alloc_qreg("shpzdiv.par");
    let counter = c.alloc_qreg_bits("shpzdiv.ctr", trailmix_counter_width());
    let load_p = |c: &mut Circuit, reg: &[QReg]| {
        for (j, q) in reg.iter().enumerate() {
            if j < 256 && (p_bytes[j / 8] >> (j % 8)) & 1 == 1 {
                c.x(q);
            }
        }
    };
    load_p(c, &aa); // A = p
    c.x(&ccb[0]); // cb = 1
    c.x(&parity); // parity = 1

    // --- forward inversion: 1/|dx| in cb (up to the parity bit) ---
    shrunken_pz_invert_forward(
        c, &mut aa, &mut dx, &mut cca, &mut ccb, &mut qq, &counter, &parity, &s_rot, &off,
    );

    // --- TEAR DOWN the EEA pack before creating lambda. At convergence the PZ
    // state is A=0, B=1, ca=p, q=0 (all CONSTANTS) and cb=1/|dx| (the only data).
    // Free the constant registers (0-Toffoli uncompute) so only cb is live during
    // the multiply -- saves ~ca(258) qubits at the peak. Re-create them (cheap)
    // before the backward. ---
    let (ta, tb, tca, tq) = (aa.len(), dx.len(), cca.len(), qq.len());
    load_p(c, &cca); // ca: p -> 0
    c.x(&dx[0]); // B: 1 -> 0
    for q in std::mem::take(&mut aa) {
        c.zero_and_free(q); // A = 0
    }
    for q in std::mem::take(&mut dx) {
        c.zero_and_free(q); // B = 0
    }
    for q in std::mem::take(&mut cca) {
        c.zero_and_free(q); // ca = 0
    }
    for q in std::mem::take(&mut qq) {
        c.zero_and_free(q); // q = 0
    }

    // --- lambda = dy * (1/|dx|), parity/sign corrected (only cb live in the pack) ---
    let cb_w = ccb.len();
    shrunken_pz_resize(c, &mut ccb, 257, "cb"); // pad the inverse to 257 for mod_mul
    let lambda = c.alloc_qreg_bits("shpzdiv.lambda", 257);
    mod_mul_rfold_mbu(c, &lambda, &ccb[..257], &dy); // lambda_raw = dy * cb
    shrunken_pz_resize(c, &mut ccb, cb_w, "cb"); // restore width for the backward
                                                 // 1/dx = (-1)^{sgn + (1-parity)} * cb  ->  negate lambda when f = NOT(sgn^par).
    let f = c.alloc_qreg("shpzdiv.negf");
    c.cx(&sgn, &f);
    c.cx(&parity, &f);
    c.x(&f); // f = NOT(sgn XOR parity)
    controlled_field_neg(c, &f, &lambda);
    c.x(&f);
    c.cx(&parity, &f);
    c.cx(&sgn, &f); // uncompute f
    c.zero_and_free(f);

    // --- GHOST dy (HMR each bit, free 256q) so the reverse runs dy-free ---
    let mut ghosts = Vec::with_capacity(dy.len());
    for q in &dy {
        ghosts.push(c.hmr_ghost(q));
    }
    for q in dy {
        c.zero_and_free(q);
    }

    // --- RE-CREATE the constant pack (A=0, B=1, ca=p, q=0) for the backward ---
    aa = c.alloc_qreg_bits("shpzdiv.A", ta); // A = 0
    dx = c.alloc_qreg_bits("shpzdiv.B", tb);
    c.x(&dx[0]); // B = 1
    cca = c.alloc_qreg_bits("shpzdiv.ca", tca);
    load_p(c, &cca); // ca = p
    qq = c.alloc_qreg_bits("shpzdiv.q", tq); // q = 0

    // --- backward inversion: restore B = |dx|, uncompute cb/parity ---
    shrunken_pz_invert_backward(
        c, &mut aa, &mut dx, &mut cca, &mut ccb, &mut qq, &counter, &parity, &s_rot, &off,
    );

    // --- free the clean inversion ancillas (S_0: A=p, ca=0, cb=1, q=0, par=1) ---
    c.x(&parity);
    c.zero_and_free(parity);
    c.x(&ccb[0]); // cb: 1 -> 0
    load_p(c, &aa); // A: p -> 0
    for q in aa.into_iter().chain(cca).chain(ccb).chain(qq) {
        c.zero_and_free(q);
    }
    for q in s_rot.into_iter().chain(counter) {
        c.zero_and_free(q);
    }
    c.zero_and_free(off);

    // --- un-sign-adjust: |dx| -> dx, uncompute sgn ---
    shrunken_pz_resize(c, &mut dx, 257, "dx");
    controlled_field_neg(c, &sgn, &dx);
    compare_geq_const(c, &dx, &half_bytes, &sgn);
    c.zero_and_free(sgn);

    // --- reconstruct dy = lambda * dx and EXORCIZE the ghosts ---
    let dy_new = c.alloc_qreg_bits("shpzdiv.dy", 257);
    mod_mul_rfold_mbu(c, &dy_new, &lambda[..257], &dx);
    for (g, q) in ghosts.into_iter().zip(dy_new.iter()) {
        c.resolve_ghost(g, q);
    }

    (dx, dy_new, lambda)
}

/// CANCEL the `shrunken_pz` slope: given `lambda` = `new_dy` / `new_dx` (live, 257), drive it to
/// |0> and FREE it, with `new_dx` (dx) and `new_dy` (dy) PRESERVED. Returns
/// (`new_dx`, `new_dy`). By EC linearity `new_dy/new_dx` == lambda, so this is the
/// alt-witness cleanup that removes the slope ancilla after the point coordinates
/// are computed.
///
/// Mirror of `shrunken_pz_divide_forward`, but it GHOSTS lambda (not dy) up front so only
/// `new_dy` rides through the inversion as the passenger (peak = EEA-peak + 256, same
/// as forward). After inverting `new_dx` -> cb = `1/|new_dx`|, it recomputes
/// temp = `new_dy` * cb (parity/sign corrected) = `new_dy/new_dx` == lambda's original
/// value, resolves the lambda-ghost against temp (exorcizing it), uncomputes temp
/// via `mod_mul_rfold_mbu_undo`, then reverse-inverts to restore `new_dx`.
pub fn shrunken_pz_divide_cancel(
    c: &mut Circuit,
    mut dx: Vec<QReg>,
    dy: Vec<QReg>,
    lambda: Vec<QReg>,
) -> (Vec<QReg>, Vec<QReg>) {
    use crate::point_add::trailmix_port::arith::compare::compare_geq_const;
    use crate::point_add::trailmix_port::arith::rfold_mbu::{mod_mul_rfold_mbu, mod_mul_rfold_mbu_undo};
    use crate::point_add::trailmix_port::inversion::shrunken_pz_schedule::reg_widths;
    use crate::point_add::trailmix_port::num_bigint::BigUint;
    assert_eq!(dx.len(), 257);
    assert_eq!(dy.len(), 257);
    assert_eq!(lambda.len(), 257);
    let half_bytes = vec![
        0x18, 0xfe, 0xff, 0x7f, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f, 0x00,
    ];
    let p_bytes = crate::point_add::trailmix_port::mod_arith::SECP256K1_P_LE;

    // --- sign-adjust new_dx -> |new_dx| < p/2 ---
    let sgn = c.alloc_qreg("shpzcan.sgn");
    compare_geq_const(c, &dx, &half_bytes, &sgn);
    controlled_field_neg(c, &sgn, &dx);

    // --- GHOST lambda (HMR each bit, free 257q) so the inversion runs lambda-free;
    // new_dy is the sole 256-bit passenger (peak = EEA-peak + 256). ---
    let mut lam_ghosts = Vec::with_capacity(lambda.len());
    for q in &lambda {
        lam_ghosts.push(c.hmr_ghost(q));
    }
    for q in lambda {
        c.zero_and_free(q);
    }

    // --- set up the inversion S_0 (B = |new_dx|, A = p, cb = 1, parity = 1) ---
    let (a0, b0, ca0, cb0, q0) = reg_widths(0);
    let (wg0, wc0) = (a0.max(b0), ca0.max(cb0));
    shrunken_pz_resize(c, &mut dx, wg0, "B");
    let mut aa = c.alloc_qreg_bits("shpzcan.A", wg0);
    let mut cca = c.alloc_qreg_bits("shpzcan.ca", wc0);
    let mut ccb = c.alloc_qreg_bits("shpzcan.cb", wc0);
    let mut qq = c.alloc_qreg_bits("shpzcan.q", q0.max(1));
    let s_rot = c.alloc_qreg_bits("shpzcan.srot", trailmix_srot_width());
    let off = c.alloc_qreg("shpzcan.off");
    let parity = c.alloc_qreg("shpzcan.par");
    let counter = c.alloc_qreg_bits("shpzcan.ctr", trailmix_counter_width());
    let load_p = |c: &mut Circuit, reg: &[QReg]| {
        for (j, q) in reg.iter().enumerate() {
            if j < 256 && (p_bytes[j / 8] >> (j % 8)) & 1 == 1 {
                c.x(q);
            }
        }
    };
    load_p(c, &aa);
    c.x(&ccb[0]);
    c.x(&parity);

    // --- forward inversion: 1/|new_dx| in cb (passenger: new_dy) ---
    shrunken_pz_invert_forward(
        c, &mut aa, &mut dx, &mut cca, &mut ccb, &mut qq, &counter, &parity, &s_rot, &off,
    );

    // --- tear down the constant pack (A=0,B=1,ca=p,q=0); keep cb=1/|new_dx| ---
    let (ta, tb, tca, tq) = (aa.len(), dx.len(), cca.len(), qq.len());
    load_p(c, &cca);
    c.x(&dx[0]);
    for q in std::mem::take(&mut aa) {
        c.zero_and_free(q);
    }
    for q in std::mem::take(&mut dx) {
        c.zero_and_free(q);
    }
    for q in std::mem::take(&mut cca) {
        c.zero_and_free(q);
    }
    for q in std::mem::take(&mut qq) {
        c.zero_and_free(q);
    }

    // --- temp = new_dy * (1/|new_dx|), parity/sign corrected = new_dy/new_dx, the
    // original value of lambda. Resolve the lambda-ghost against it, then uncompute
    // temp. ---
    let cb_w = ccb.len();
    shrunken_pz_resize(c, &mut ccb, 257, "cb");
    let temp = c.alloc_qreg_bits("shpzcan.temp", 257);
    mod_mul_rfold_mbu(c, &temp, &ccb[..257], &dy); // temp_raw = dy * cb
    let f = c.alloc_qreg("shpzcan.negf");
    c.cx(&sgn, &f);
    c.cx(&parity, &f);
    c.x(&f); // f = NOT(sgn XOR parity)
    controlled_field_neg(c, &f, &temp); // temp = +/-(dy*cb) = new_dy/new_dx
    for (g, q) in lam_ghosts.into_iter().zip(temp.iter()) {
        c.resolve_ghost(g, q); // exorcize lambda (temp == lambda's value)
    }
    controlled_field_neg(c, &f, &temp); // un-correct: temp = dy*cb (raw)
    c.x(&f);
    c.cx(&parity, &f);
    c.cx(&sgn, &f); // uncompute f
    c.zero_and_free(f);
    mod_mul_rfold_mbu_undo(c, &temp, &ccb[..257], &dy); // temp -> 0
    for q in temp {
        c.zero_and_free(q);
    }
    shrunken_pz_resize(c, &mut ccb, cb_w, "cb");

    // --- re-create the pack, backward inversion (restore B=|new_dx|) ---
    aa = c.alloc_qreg_bits("shpzcan.A", ta);
    dx = c.alloc_qreg_bits("shpzcan.B", tb);
    c.x(&dx[0]);
    cca = c.alloc_qreg_bits("shpzcan.ca", tca);
    load_p(c, &cca);
    qq = c.alloc_qreg_bits("shpzcan.q", tq);
    shrunken_pz_invert_backward(
        c, &mut aa, &mut dx, &mut cca, &mut ccb, &mut qq, &counter, &parity, &s_rot, &off,
    );

    // --- free the clean inversion ancillas (S_0: A=p, ca=0, cb=1, q=0, par=1) ---
    c.x(&parity);
    c.zero_and_free(parity);
    c.x(&ccb[0]);
    load_p(c, &aa);
    for q in aa.into_iter().chain(cca).chain(ccb).chain(qq) {
        c.zero_and_free(q);
    }
    for q in s_rot.into_iter().chain(counter) {
        c.zero_and_free(q);
    }
    c.zero_and_free(off);

    // --- un-sign-adjust: |new_dx| -> new_dx, uncompute sgn ---
    shrunken_pz_resize(c, &mut dx, 257, "dx");
    controlled_field_neg(c, &sgn, &dx);
    compare_geq_const(c, &dx, &half_bytes, &sgn);
    c.zero_and_free(sgn);

    (dx, dy)
}

#[cfg(test)]
mod tests;
