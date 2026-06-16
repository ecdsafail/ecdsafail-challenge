//! Fast multiply-by-2^k mod q (secp256k1, q = 2^256 - f, f = 2^32+977), via the
//! reduction `tmp = carries*f; a[low] += tmp; clear tmp`:
//!
//!   2^k * x  ==  (x<<k mod 2^256)  +  carries * f   (mod q),  carries = x>>(256-k).
//!
//! One reduction (a single q-q add of carries*f) replaces the k single-doubling
//! f-folds. Wins for k>=3 (the q-q add ~ 2 const-folds, so it beats k folds once
//! k>~2). cx-clean of the carry slots needs f mod 2^k == 1 (true k<=4, f mod 16==1).
//! Windowed/approximate like `mod_double_pm` (carry beyond `lsbs` dropped).
//!
//! `tmp = carries*f` is a tiny carry-save multiply (partial products are just
//! copies of carries[i] at f's bit positions -- no Toffoli for products, only the
//! column compressors). It's planned as a gate tape; since cx/ccx are self-inverse
//! the uncompute is the tape replayed in reverse -- so the primitive is alloc-clean
//! with no hand-mirroring.

use crate::point_add::trailmix_port::arith::khattar_gidney::unary_iterate;
use crate::point_add::trailmix_port::arith::shift::{left_shift, right_shift};
use crate::point_add::trailmix_port::circuit::{Circuit, QReg};

/// As [`qrom_f_fold`] but adds `f · (addr interpreted as a k-bit SIGNED two's-
/// complement value)` — i.e. the top address bit subtracts `f·2^(k-1)`. Used for
/// the radix square's signed overflow fold (O ∈ {0,12..15} ≡ {0,-4..-1}).
pub fn qrom_f_fold_signed(circ: &mut Circuit, addr: &[QReg], x: &[QReg], lsbs: usize) {
    use crate::point_add::trailmix_port::arith::gidney_const_adder::hybrid_add_refs;
    let k = addr.len();
    assert!((2..=5).contains(&k));
    assert!(lsbs <= x.len());
    let f = super::F_SECP256K1 as i128;
    let modulus = 1i128 << lsbs;
    let scratch: Vec<QReg> = (0..lsbs).map(|_| circ.alloc_qreg("qffs.s")).collect();
    let addr_refs: Vec<&QReg> = addr.iter().collect();
    let xor_fa = |circ: &mut Circuit, scratch: &[QReg], a: usize, gate: &QReg| {
        let sa = if a >= (1 << (k - 1)) {
            a as i128 - (1 << k)
        } else {
            a as i128
        };
        let val = (f * sa).rem_euclid(modulus) as u128; // two's-complement lsbs-bit
        for (bit, sb) in scratch.iter().enumerate() {
            if (val >> bit) & 1 == 1 {
                circ.cx(gate, sb);
            }
        }
    };
    unary_iterate(circ, &addr_refs, 1 << k, |c, a, g| {
        xor_fa(c, &scratch, a, g)
    });
    let xr: Vec<&QReg> = x[0..lsbs].iter().collect();
    let sr: Vec<&QReg> = scratch.iter().collect();
    hybrid_add_refs(circ, &xr, &sr, lsbs - 1);
    unary_iterate(circ, &addr_refs, 1 << k, |c, a, g| {
        xor_fa(c, &scratch, a, g)
    });
    for q in scratch {
        circ.zero_and_free(q);
    }
}

/// `x[0..lsbs] += (subtract ? -1 : +1) · f · value(addr)` (mod 2^lsbs, carry
/// beyond `lsbs` dropped — windowed), via a 2^k-entry QROM lookup of the
/// precomputed constants `{f·a : a in 0..2^k}`, k = addr.len() (2..=5).
///
/// `f·O` is linear in O, but doing it as k controlled const-adds pays k separate
/// carry propagations through the full `lsbs` window. The QROM instead decodes
/// `addr` to a one-hot (`unary_iterate`, cheap), XORs the matching `f·a` into a
/// small scratch (free), and pays ONE vented add + one carry tail. The XOR body
/// is self-inverse, so re-running the iteration unloads the scratch — no HMR,
/// phase-clean by construction, alloc-clean.
pub fn qrom_f_fold(circ: &mut Circuit, addr: &[QReg], x: &[QReg], lsbs: usize, subtract: bool) {
    use crate::point_add::trailmix_port::arith::gidney_const_adder::hybrid_add_refs;
    let k = addr.len();
    assert!((2..=5).contains(&k), "qrom_f_fold: k must be 2..=5");
    assert!(lsbs <= x.len());
    let f = super::F_SECP256K1 as u128;
    // scratch spans the full window so the carry tail propagates; only the low
    // 33+k bits ever get written, the rest stay |0> (zero-padding the addend).
    let scratch: Vec<QReg> = (0..lsbs).map(|_| circ.alloc_qreg("qff.s")).collect();
    let s_hi = (33 + k).min(lsbs);
    let addr_refs: Vec<&QReg> = addr.iter().collect();
    let xor_fa = |circ: &mut Circuit, scratch: &[QReg], a: usize, gate: &QReg| {
        let fa = f * (a as u128);
        for (bit, sb) in scratch.iter().enumerate().take(s_hi) {
            if (fa >> bit) & 1 == 1 {
                circ.cx(gate, sb);
            }
        }
    };
    // load: scratch ^= f·addr (only the matching one-hot fires).
    unary_iterate(circ, &addr_refs, 1 << k, |c, a, g| {
        xor_fa(c, &scratch, a, g)
    });
    // x[0..lsbs] (±)= scratch.
    let xr: Vec<&QReg> = x[0..lsbs].iter().collect();
    let sr: Vec<&QReg> = scratch.iter().collect();
    if subtract {
        for q in &xr {
            circ.x(q);
        }
    }
    hybrid_add_refs(circ, &xr, &sr, lsbs - 1);
    if subtract {
        for q in &xr {
            circ.x(q);
        }
    }
    // unload (XOR body is self-inverse).
    unary_iterate(circ, &addr_refs, 1 << k, |c, a, g| {
        xor_fa(c, &scratch, a, g)
    });
    for q in scratch {
        circ.zero_and_free(q);
    }
}

const F_BITS: [usize; 7] = [0, 4, 6, 7, 8, 9, 32]; // set bits of f = 2^32 + 977

/// Operand of a planned gate: scratch[i], carries[i], or tmp[i].
#[derive(Clone, Copy)]
enum Opd {
    S(usize),
    C(usize),
    T(usize),
}

#[derive(Clone, Copy)]
enum Gate {
    Cx(Opd, Opd),
    Ccx(Opd, Opd, Opd),
}

/// Plan `tmp[0..tlen] = carries*f` as a carry-save multiply: alloc the scratch
/// qubits and record the cx/ccx tape. No gates emitted yet. Returns (tape, scratch).
fn plan_carries_times_f(circ: &mut Circuit, k: usize, tlen: usize) -> (Vec<Gate>, Vec<QReg>) {
    let mut scratch: Vec<QReg> = Vec::new();
    let mut tape: Vec<Gate> = Vec::new();
    // column contents as operands (refs), survivors written to tmp.
    let mut col: Vec<Vec<Opd>> = (0..tlen + 2).map(|_| Vec::new()).collect();
    // partial products: copy carries[i] into column i+j for each f-bit j.
    for &j in &F_BITS {
        for i in 0..k {
            let c = i + j;
            if c < tlen {
                let idx = scratch.len();
                scratch.push(circ.alloc_qreg("ctf_pp"));
                tape.push(Gate::Cx(Opd::C(i), Opd::S(idx))); // copy carries[i]
                col[c].push(Opd::S(idx));
            }
        }
    }
    for c in 0..tlen {
        let mut carries_up: Vec<Opd> = Vec::new();
        while col[c].len() > 1 {
            let cout = {
                let idx = scratch.len();
                scratch.push(circ.alloc_qreg("ctf_cy"));
                Opd::S(idx)
            };
            if col[c].len() >= 3 {
                let cin = col[c].pop().unwrap();
                let b = col[c].pop().unwrap();
                let a = col[c].pop().unwrap();
                // fa: cout=MAJ(a,b,cin); b=sum.  gates: ccx(a,b,cout) cx(a,b) ccx(b,cin,cout) cx(cin,b)
                tape.push(Gate::Ccx(a, b, cout));
                tape.push(Gate::Cx(a, b));
                tape.push(Gate::Ccx(b, cin, cout));
                tape.push(Gate::Cx(cin, b));
                carries_up.push(cout);
                col[c].push(b);
            } else {
                let b = col[c].pop().unwrap();
                let a = col[c].pop().unwrap();
                // ha: cout=a&b; b=a^b.  gates: ccx(a,b,cout) cx(a,b)
                tape.push(Gate::Ccx(a, b, cout));
                tape.push(Gate::Cx(a, b));
                carries_up.push(cout);
                col[c].push(b);
            }
        }
        if let Some(bit) = col[c].pop() {
            tape.push(Gate::Cx(bit, Opd::T(c))); // survivor -> tmp[c]
        }
        col[c + 1].extend(carries_up);
    }
    (tape, scratch)
}

fn emit_gate(circ: &mut Circuit, g: &Gate, carries: &[&QReg], tmp: &[QReg], scratch: &[QReg]) {
    let r = |o: &Opd| -> &QReg {
        match *o {
            Opd::S(i) => &scratch[i],
            Opd::C(i) => carries[i],
            Opd::T(i) => &tmp[i],
        }
    };
    match g {
        Gate::Cx(a, b) => circ.cx(r(a), r(b)),
        Gate::Ccx(a, b, c) => circ.ccx(r(a), r(b), r(c)),
    }
}

/// `a := (2^k * a) mod q` (secp256k1), windowed-approximate, ALLOC-CLEAN.
/// `a.len() == 256 + k`, `a[0..256]` = x < q, `a[256..256+k] = |0>`, k in 1..=4.
pub fn mod_mul_pow2_pm_secp256k1(circ: &mut Circuit, a: &[QReg], k: usize) {
    let n = 256usize;
    assert!((1..=4).contains(&k), "cx-clean needs f mod 2^k == 1 (k<=4)");
    assert_eq!(a.len(), n + k, "a must be n+k bits");
    let prev = circ.push_section("mul_pow2");

    // 1) shift left by k (Toffoli-free). carries = a[n..n+k] = top k bits of x.
    for _ in 0..k {
        left_shift(circ, a);
    }
    let carries: Vec<&QReg> = (0..k).map(|i| &a[n + i]).collect();

    // 2) tmp_ext[0..tlen] = carries*f; tmp_ext[tlen..] stay |0> (windowing pad).
    let tlen = 33 + k + 1;
    let lsbs = tlen + 30;
    let tmp: Vec<QReg> = (0..lsbs).map(|_| circ.alloc_qreg("mp2_tmp")).collect();
    let (tape, scratch) = plan_carries_times_f(circ, k, tlen);
    for g in &tape {
        emit_gate(circ, g, &carries, &tmp, &scratch);
    }

    // 3) a[0..lsbs] += tmp  (windowed q-q add; carry beyond lsbs dropped). Use the
    //    UNCONDITIONAL measurement-vented adder (~n Toffoli, full vents) -- the
    //    headroom supplies the clean vent ancillae; this add is the kept op (not
    //    reversed), so its internal HMR is self-contained.
    crate::point_add::trailmix_port::arith::gidney_const_adder::hybrid_add(circ, &a[0..lsbs], &tmp, lsbs - 1);

    // 4) uncompute tmp: replay the tape reversed (cx/ccx are self-inverse).
    for g in tape.iter().rev() {
        emit_gate(circ, g, &carries, &tmp, &scratch);
    }
    for q in scratch {
        circ.zero_and_free(q);
    }

    // 5) clear the k carry slots: a[0..k] == carries after the add (shifted low
    //    bits were 0; (carries*f) mod 2^k == carries since f mod 2^k == 1, k<=4).
    for i in 0..k {
        circ.cx(&a[i], &a[n + i]);
    }
    for q in tmp {
        circ.zero_and_free(q);
    }
    circ.pop_section(&prev);
}

/// Inverse of [`mod_mul_pow2_pm_secp256k1`]: `a := (a / 2^k) mod q` (maps
/// `2^k*x -> x`). Mirror of the forward: re-derive carries, rebuild tmp, `a -= tmp`
/// (X-sandwich hybrid_add), unbuild tmp, right-shift. Exactly inverts the forward.
pub fn mod_mul_pow2_pm_secp256k1_reverse(circ: &mut Circuit, a: &[QReg], k: usize) {
    let n = 256usize;
    assert!((1..=4).contains(&k));
    assert_eq!(a.len(), n + k);
    let prev = circ.push_section("mul_pow2_rev");
    let carries: Vec<&QReg> = (0..k).map(|i| &a[n + i]).collect();
    let tlen = 33 + k + 1;
    let lsbs = tlen + 30;

    // 1) re-derive carry slots (inverse of the forward's clear): a[n+i] ^= a[i].
    for i in 0..k {
        circ.cx(&a[i], &a[n + i]);
    }
    // 2) rebuild tmp = carries*f.
    let tmp: Vec<QReg> = (0..lsbs).map(|_| circ.alloc_qreg("mp2r_tmp")).collect();
    let (tape, scratch) = plan_carries_times_f(circ, k, tlen);
    for g in &tape {
        emit_gate(circ, g, &carries, &tmp, &scratch);
    }
    // 3) a[0..lsbs] -= tmp  (X-sandwich the vented add).
    for q in &a[0..lsbs] {
        circ.x(q);
    }
    crate::point_add::trailmix_port::arith::gidney_const_adder::hybrid_add(circ, &a[0..lsbs], &tmp, lsbs - 1);
    for q in &a[0..lsbs] {
        circ.x(q);
    }
    // 4) unbuild tmp.
    for g in tape.iter().rev() {
        emit_gate(circ, g, &carries, &tmp, &scratch);
    }
    for q in scratch {
        circ.zero_and_free(q);
    }
    for q in tmp {
        circ.zero_and_free(q);
    }
    // 5) right-shift by k (inverse of the forward's left-shift).
    for _ in 0..k {
        right_shift(circ, a);
    }
    circ.pop_section(&prev);
}

/// `x -= ctrl·(λ·2^j mod q)`, in place from λ's bits (no copy): the pre-wrapped
/// shift `λ·2^j mod q = λ[..n-j]<<j + f·λ[n-j..n]`. j in 0..=3. x is >= n+4 bits;
/// borrows run into the n+4 slack (bounded by assumption x±λ·2^j < 2^{n+1+slack}).
fn sub_lambda_pow2(circ: &mut Circuit, ctrl: &QReg, lam: &[QReg], x: &[QReg], j: usize) {
    use crate::point_add::trailmix_port::arith::gidney_const_adder::{controlled_add_const_gidney, hybrid_add_refs};
    let n = 256usize;
    let f_bytes = super::F_SECP256K1.to_le_bytes();
    let w = x.len(); // n+4
                     // sub1: x[j..w] -= ctrl·(λ[0..n-j])  (the shifted low part; borrow into slack).
                     //   X-sandwich the vented controlled add of λ[0..n-j] (zero-padded to w-j).
    let pad: Vec<QReg> = (0..(w - j) - (n - j))
        .map(|_| circ.alloc_qreg("slp_pad"))
        .collect();
    {
        let target: Vec<&QReg> = x[j..w].iter().collect();
        for q in &target {
            circ.x(q);
        }
        let addend: Vec<&QReg> = lam[0..n - j].iter().chain(pad.iter()).collect();
        // controlled (ctrl) vented add of addend into target, full vents.
        let s1 = circ.push_section("slp1");
        crate::point_add::trailmix_port::arith::gidney_const_adder::controlled_hybrid_add_refs(
            circ,
            ctrl,
            &target,
            &addend,
            target.len().saturating_sub(1),
        );
        circ.pop_section(&s1);
        for q in &target {
            circ.x(q);
        }
        let _ = hybrid_add_refs; // (kept import tidy)
    }
    for q in pad {
        circ.zero_and_free(q);
    }
    let s2 = circ.push_section("slp2");
    // sub2: x -= ctrl·f·λ[n-j..n] (the wrap of λ·2^j's top j bits). j==1 is one
    // const-sub; j>=2 is a tiny QROM of the gated j-bit slice f·(ctrl & λ[n-j..n]).
    let lsbs = 33 + j + 30;
    if j == 1 {
        let dh = lsbs + (lsbs - 1);
        let t = circ.alloc_qreg("slp_t");
        circ.ccx(ctrl, &lam[n - 1], &t); // t = ctrl & λ[255]
        for q in &x[0..lsbs] {
            circ.x(q);
        }
        controlled_add_const_gidney(circ, &t, &x[0..lsbs], &f_bytes, &x[lsbs..dh]);
        for q in &x[0..lsbs] {
            circ.x(q);
        }
        circ.ccx(ctrl, &lam[n - 1], &t);
        circ.zero_and_free(t);
    } else if j >= 2 {
        // temp[b] = ctrl & λ[n-j+b]; x -= f·temp via QROM (f·0 = 0 when ctrl=0).
        let temp: Vec<QReg> = (0..j).map(|_| circ.alloc_qreg("slp_qa")).collect();
        for b in 0..j {
            circ.ccx(ctrl, &lam[n - j + b], &temp[b]);
        }
        qrom_f_fold(circ, &temp, x, lsbs, true);
        for b in 0..j {
            circ.ccx(ctrl, &lam[n - j + b], &temp[b]);
        }
        for q in temp {
            circ.zero_and_free(q);
        }
    }
    circ.pop_section(&s2);
}

/// Per-window reduce (validated 0/200000 + 0/20000): after the 4 subs, x_low = U_low,
/// x[256..260] = O (= (-carry)&15, carry = #subs that overflowed bit 256, <=4). Restore
/// x to [~0,2^256) with x[256..260]=0:
///   (a) SIGNED fold using x[256..260] as controls: x += f·x256 + 2f·x257 + 4f·x258 − 8f·x259
///       (x259 = sign bit, subtracts 8f). After: x_low = value mod q (in [−4f, 2^256)).
///   (b) recompute carry from the now-stable post-fold x[208..256] via a top-slice add of
///       λ[208−j..256−j] (the top slice of λ·2^j mod q is just that λ-slice; f-corr is low),
///       then clear x[256..260] by a mod-16 add of carry (O + carry ≡ 0 mod 16). Uncompute.
fn reduce_window_radix4(
    circ: &mut Circuit,
    ctrl: &QReg,
    lam: &[QReg],
    x: &[QReg],
    w: usize,
    c: &QReg,
) {
    use crate::point_add::trailmix_port::arith::gidney_const_adder::{controlled_hybrid_add_refs, hybrid_add_refs};
    let lsf = 48usize; // fold low-window: holds 8f (~36b) + carry margin
                       // (a) signed fold via tiny QROM: x[0..lsf] += f·(signed O), O = x[256..260].
                       //     One 16-entry lookup + one add, vs the four margin-window const-adds.
    let sf = circ.push_section("redfold");
    qrom_f_fold_signed(circ, &x[256..260], x, lsf);
    circ.pop_section(&sf);
    let sr = circ.push_section("redrecomp");
    // (b) recompute carry into acc[slice..slice+4] from post-fold x[base..256] + Σ
    //     λ-slices (base = 256-slice; carry-in below `base` ignored, ~2^-slice/window).
    let slice = 28usize;
    let base = 256 - slice;
    let acc: Vec<QReg> = (0..slice + 4).map(|_| circ.alloc_qreg("sq4.acc")).collect();
    let pad: Vec<QReg> = (0..4).map(|_| circ.alloc_qreg("sq4.pad")).collect();
    for i in 0..slice {
        circ.cx(&x[base + i], &acc[i]);
    }
    for j in 0..4usize {
        circ.ccx(ctrl, &lam[4 * w + j], c);
        let addend: Vec<&QReg> = (0..slice)
            .map(|i| &lam[base - j + i])
            .chain(pad.iter())
            .collect();
        let target: Vec<&QReg> = acc.iter().collect();
        controlled_hybrid_add_refs(circ, c, &target, &addend, target.len() - 1);
        circ.ccx(ctrl, &lam[4 * w + j], c);
    }
    // clear x[256..260] : x[256..260] += carry (mod 16, carry in acc[48..52]).
    let xtop: Vec<&QReg> = x[256..260].iter().collect();
    let carry: Vec<&QReg> = acc[slice..slice + 4].iter().collect();
    hybrid_add_refs(circ, &xtop, &carry, xtop.len() - 1);
    // uncompute acc: reverse the Σ λ-slice adds, then un-copy.
    for j in (0..4usize).rev() {
        circ.ccx(ctrl, &lam[4 * w + j], c);
        let addend: Vec<&QReg> = (0..slice)
            .map(|i| &lam[base - j + i])
            .chain(pad.iter())
            .collect();
        let target: Vec<&QReg> = acc.iter().collect();
        for q in &target {
            circ.x(q);
        }
        controlled_hybrid_add_refs(circ, c, &target, &addend, target.len() - 1);
        for q in &target {
            circ.x(q);
        }
        circ.ccx(ctrl, &lam[4 * w + j], c);
    }
    for i in 0..slice {
        circ.cx(&x[base + i], &acc[i]);
    }
    for q in pad {
        circ.zero_and_free(q);
    }
    for q in acc {
        circ.zero_and_free(q);
    }
    circ.pop_section(&sr);
}

/// `x -= ctrl·λ² mod q` via the radix-2^4 reverse-Horner (validated structure,
/// scripts/radix4_square_proto.py). `x` is n+4 = 260 bits (overflow slack), λ is
/// n+1 = 257. LSB-first windows; per window /2^4 (skip first) + 4 pre-wrapped subs +
/// per-window reduce; x2^4 batch. Uses mod_mul_pow2[_reverse] for the clean /2^4, x2^4.
pub fn controlled_mod_square_sub_pm_secp256k1_radix4(
    circ: &mut Circuit,
    ctrl: &QReg,
    lam: &[QReg],
    x: &[QReg],
) {
    let n = 256usize;
    assert_eq!(lam.len(), n + 1);
    assert_eq!(x.len(), n + 4);
    let prev = circ.push_section("sqr_radix4");
    let c = circ.alloc_qreg("sqr4.c");
    let mut first = true;
    let prof = std::env::var("SQR4_PROF").is_ok();
    let (mut t_div, mut t_sub, mut t_red, mut t_mul) = (0u64, 0u64, 0u64, 0u64);
    let tof = |c: &Circuit| c.ccx_emitted + c.ccz_emitted;
    for w in 0..(n / 4) {
        if !first {
            let t = tof(circ);
            mod_mul_pow2_pm_secp256k1_reverse(circ, x, 4);
            t_div += tof(circ) - t;
        }
        first = false;
        let t = tof(circ);
        for j in 0..4usize {
            let bit = 4 * w + j;
            circ.ccx(ctrl, &lam[bit], &c); // c = ctrl & λ[bit]
            sub_lambda_pow2(circ, &c, lam, x, j);
            circ.ccx(ctrl, &lam[bit], &c);
        }
        t_sub += tof(circ) - t;
        let t = tof(circ);
        reduce_window_radix4(circ, ctrl, lam, x, w, &c);
        t_red += tof(circ) - t;
    }
    let t = tof(circ);
    for _ in 0..(n / 4 - 1) {
        mod_mul_pow2_pm_secp256k1(circ, x, 4);
    }
    t_mul += tof(circ) - t;
    if prof {
        eprintln!("SQR4_PROF: div={t_div} subs={t_sub} reduce={t_red} mul={t_mul}");
    }
    circ.zero_and_free(c);
    circ.pop_section(&prev);
}

/// `x -= ctrl·λ² mod q` via the schoolbook square with a running multiple
/// `m = λ·2^j mod q` (doubled by `mod_mul_pow2`, forward — no `÷2`, no Horner
/// doubling of x, no parity issue). Standard mod-subs of `m` keep x in [0,q).
/// Correct/simple baseline (uses a copy of λ for `m`); the in-place radix form is
/// the optimization. x and λ are n+1 = 257 bits.
pub fn controlled_mod_square_sub_pm_secp256k1_schoolbook(
    circ: &mut Circuit,
    ctrl: &QReg,
    lam: &[QReg],
    x: &[QReg],
) {
    let n = 256usize;
    assert_eq!(lam.len(), n + 1);
    assert_eq!(x.len(), n + 1);
    let prev = circ.push_section("sqr_schoolbook");
    // m = copy of λ (the running multiple).
    let m: Vec<QReg> = (0..n + 1).map(|_| circ.alloc_qreg("sqr_sb.m")).collect();
    for i in 0..n + 1 {
        circ.cx(&lam[i], &m[i]);
    }
    let c = circ.alloc_qreg("sqr_sb.c");
    for j in 0..n {
        circ.ccx(ctrl, &lam[j], &c); // c = ctrl & λ[j]
        super::pm_prims::controlled_mod_add_pm_secp256k1_reverse_vents(circ, &c, &m, x, usize::MAX); // x -= c·m  (mod q)
        circ.ccx(ctrl, &lam[j], &c);
        if j + 1 < n {
            mod_mul_pow2_pm_secp256k1(circ, &m, 1); // m *= 2 mod q
        }
    }
    // uncompute m: reverse the (n-1) doublings, then un-copy.
    for _ in 0..n - 1 {
        mod_mul_pow2_pm_secp256k1_reverse(circ, &m, 1);
    }
    for i in 0..n + 1 {
        circ.cx(&lam[i], &m[i]);
    }
    circ.zero_and_free(c);
    for q in m {
        circ.zero_and_free(q);
    }
    circ.pop_section(&prev);
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigUint;
    use rand::{thread_rng, RngCore};

    fn q() -> BigUint {
        (BigUint::from(1u32) << 256) - BigUint::from((1u64 << 32) + 977)
    }

    #[test]
    fn mul_pow2_value_clean_and_cost() {
        let qq = q();
        for k in [1usize, 2, 3, 4] {
            let n = 256;
            let mut circ = Circuit::new();
            circ.set_max_qubit_peak(4096);
            let a: Vec<QReg> = (0..n + k)
                .map(|i| circ.alloc_qreg(&format!("a{i}")))
                .collect();
            let mut rng = thread_rng();
            let mut xs: Vec<BigUint> = Vec::with_capacity(64);
            for shot in 0..64 {
                let mut bytes = [0u8; 32];
                rng.fill_bytes(&mut bytes);
                let x = BigUint::from_bytes_le(&bytes) % &qq;
                circ.sim_load_reg_bytes_shot(&a[..n], &x.to_bytes_le(), shot);
                xs.push(x);
            }
            let t0 = circ.ccx_emitted + circ.ccz_emitted;
            mod_mul_pow2_pm_secp256k1(&mut circ, &a, k);
            let tof = (circ.ccx_emitted + circ.ccz_emitted) - t0;

            // alloc-clean: a is the only thing that should be live (carries cleared).
            let (sim, det) = circ.destroy_sim(a.into_iter().collect());
            let mut bad = 0;
            for (shot, x) in xs.iter().enumerate() {
                let mut got = BigUint::from(0u32);
                for i in 0..n {
                    if sim.read_bit_shot(&det[i], shot) == 1 {
                        got |= BigUint::from(1u32) << i;
                    }
                }
                let want = (BigUint::from(1u32) << k) * x % &qq;
                if got % &qq != want {
                    bad += 1;
                }
            }
            eprintln!(
                "mul_pow2 k={k}: tof={tof} (vs k folds ~{}), value bad={bad}/64",
                k * 62
            );
            assert_eq!(bad, 0, "k={k}: {bad}/64 value mismatches");
        }
    }

    #[test]
    fn radix4_square_value() {
        let qq = q();
        let n = 256;
        let mut circ = Circuit::new();
        circ.set_max_qubit_peak(8192);
        let lam: Vec<QReg> = (0..n + 1)
            .map(|i| circ.alloc_qreg(&format!("l{i}")))
            .collect();
        let x: Vec<QReg> = (0..n + 4)
            .map(|i| circ.alloc_qreg(&format!("x{i}")))
            .collect();
        let ctrl = circ.alloc_qreg("ctrl");
        let mut rng = thread_rng();
        let mut ls = Vec::new();
        let mut xs = Vec::new();
        for shot in 0..64 {
            let mut lb = [0u8; 32];
            let mut xb = [0u8; 32];
            rng.fill_bytes(&mut lb);
            rng.fill_bytes(&mut xb);
            let l = BigUint::from_bytes_le(&lb) % &qq;
            let xv = BigUint::from_bytes_le(&xb) % &qq;
            circ.sim_load_reg_bytes_shot(&lam[..n], &l.to_bytes_le(), shot);
            circ.sim_load_reg_bytes_shot(&x[..n], &xv.to_bytes_le(), shot);
            circ.sim_load_reg_bytes_shot(std::slice::from_ref(&ctrl), &[1], shot);
            ls.push(l);
            xs.push(xv);
        }
        let t0 = circ.ccx_emitted + circ.ccz_emitted;
        controlled_mod_square_sub_pm_secp256k1_radix4(&mut circ, &ctrl, &lam, &x);
        let tof = (circ.ccx_emitted + circ.ccz_emitted) - t0;
        if std::env::var("SQR4_PROF").is_ok() {
            let mut rows: Vec<(String, u64)> = circ
                .executed_toffoli_by_section
                .iter()
                .filter(|(p, _)| p.contains("sqr_radix4"))
                .map(|(p, c)| (p.clone(), *c))
                .collect();
            rows.sort_by(|a, b| b.1.cmp(&a.1));
            for (p, c) in rows.iter().take(20) {
                eprintln!("  SEC {c:>9}  {p}");
            }
        }
        let mut outs: Vec<QReg> = Vec::new();
        outs.extend(x);
        outs.extend(lam);
        outs.push(ctrl);
        let (sim, det) = circ.destroy_sim(outs);
        let mut bad = 0;
        for shot in 0..64 {
            let mut got = BigUint::from(0u32);
            for i in 0..n {
                if sim.read_bit_shot(&det[i], shot) == 1 {
                    got |= BigUint::from(1u32) << i;
                }
            }
            let want = (&xs[shot] + &qq - (&ls[shot] * &ls[shot]) % &qq) % &qq;
            if got % &qq != want {
                bad += 1;
            }
        }
        eprintln!("radix4_square_value: tof={tof}, bad={bad}/64");
        assert_eq!(bad, 0, "{bad}/64 mismatches");
    }

    #[test]
    fn qrom_f_fold_value_and_cost() {
        let f = super::super::F_SECP256K1 as u128;
        let k = 4;
        let lsbs = 67;
        for sub in [false, true] {
            let mut circ = Circuit::new();
            circ.set_max_qubit_peak(4096);
            let addr: Vec<QReg> = (0..k).map(|i| circ.alloc_qreg(&format!("a{i}"))).collect();
            let x: Vec<QReg> = (0..lsbs)
                .map(|i| circ.alloc_qreg(&format!("x{i}")))
                .collect();
            let mut rng = thread_rng();
            let mut avs = Vec::new();
            let mut xvs = Vec::new();
            for shot in 0..64 {
                let a = (rng.next_u32() as usize) & ((1 << k) - 1);
                let xv = (rng.next_u64() as u128) & ((1u128 << lsbs) - 1);
                let ab: Vec<u8> = (0..((k + 7) / 8)).map(|b| (a >> (8 * b)) as u8).collect();
                let xb: Vec<u8> = (0..((lsbs + 7) / 8))
                    .map(|b| (xv >> (8 * b)) as u8)
                    .collect();
                circ.sim_load_reg_bytes_shot(&addr, &ab, shot);
                circ.sim_load_reg_bytes_shot(&x, &xb, shot);
                avs.push(a as u128);
                xvs.push(xv);
            }
            let t0 = circ.ccx_emitted + circ.ccz_emitted;
            qrom_f_fold(&mut circ, &addr, &x, lsbs, sub);
            let tof = (circ.ccx_emitted + circ.ccz_emitted) - t0;
            let mut outs: Vec<QReg> = Vec::new();
            outs.extend(x);
            outs.extend(addr);
            let (sim, det) = circ.destroy_sim(outs);
            let mask = (1u128 << lsbs) - 1;
            let mut bad = 0;
            for shot in 0..64 {
                let mut got = 0u128;
                for i in 0..lsbs {
                    if sim.read_bit_shot(&det[i], shot) == 1 {
                        got |= 1u128 << i;
                    }
                }
                let delta = (f * avs[shot]) & mask;
                let want = if sub {
                    (xvs[shot].wrapping_sub(delta)) & mask
                } else {
                    (xvs[shot] + delta) & mask
                };
                if got != want {
                    bad += 1;
                }
            }
            eprintln!("qrom_f_fold sub={sub}: tof={tof}, bad={bad}/64");
            assert_eq!(bad, 0, "sub={sub}: {bad}/64 mismatches");
        }
    }

    #[test]
    fn schoolbook_square_value_and_cost() {
        let qq = q();
        let n = 256;
        let mut circ = Circuit::new();
        circ.set_max_qubit_peak(8192);
        let lam: Vec<QReg> = (0..n + 1)
            .map(|i| circ.alloc_qreg(&format!("l{i}")))
            .collect();
        let x: Vec<QReg> = (0..n + 1)
            .map(|i| circ.alloc_qreg(&format!("x{i}")))
            .collect();
        let ctrl = circ.alloc_qreg("ctrl");
        let mut rng = thread_rng();
        let mut ls = Vec::new();
        let mut xs = Vec::new();
        for shot in 0..64 {
            let mut lb = [0u8; 32];
            let mut xb = [0u8; 32];
            rng.fill_bytes(&mut lb);
            rng.fill_bytes(&mut xb);
            let l = BigUint::from_bytes_le(&lb) % &qq;
            let xv = BigUint::from_bytes_le(&xb) % &qq;
            circ.sim_load_reg_bytes_shot(&lam[..n], &l.to_bytes_le(), shot);
            circ.sim_load_reg_bytes_shot(&x[..n], &xv.to_bytes_le(), shot);
            circ.sim_load_reg_bytes_shot(std::slice::from_ref(&ctrl), &[1], shot);
            ls.push(l);
            xs.push(xv);
        }
        let t0 = circ.ccx_emitted + circ.ccz_emitted;
        controlled_mod_square_sub_pm_secp256k1_schoolbook(&mut circ, &ctrl, &lam, &x);
        let tof = (circ.ccx_emitted + circ.ccz_emitted) - t0;
        let mut outs: Vec<QReg> = Vec::new();
        outs.extend(x);
        outs.extend(lam);
        outs.push(ctrl);
        let (sim, det) = circ.destroy_sim(outs);
        let mut bad = 0;
        for shot in 0..64 {
            let mut got = BigUint::from(0u32);
            for i in 0..n {
                if sim.read_bit_shot(&det[i], shot) == 1 {
                    got |= BigUint::from(1u32) << i;
                }
            }
            let want = (&xs[shot] + &qq - (&ls[shot] * &ls[shot]) % &qq) % &qq;
            if got % &qq != want {
                bad += 1;
            }
        }
        eprintln!("schoolbook_square: tof={tof}, value bad={bad}/64");
        assert_eq!(bad, 0, "{bad}/64 mismatches");
    }

    #[test]
    fn mul_pow2_roundtrip_clean() {
        // x2^k then /2^k must be identity (and alloc-clean).
        let qq = q();
        for k in [1usize, 2, 3, 4] {
            let n = 256;
            let mut circ = Circuit::new();
            circ.set_max_qubit_peak(4096);
            let a: Vec<QReg> = (0..n + k)
                .map(|i| circ.alloc_qreg(&format!("a{i}")))
                .collect();
            let mut rng = thread_rng();
            let mut xs: Vec<BigUint> = Vec::with_capacity(64);
            for shot in 0..64 {
                let mut bytes = [0u8; 32];
                rng.fill_bytes(&mut bytes);
                let x = BigUint::from_bytes_le(&bytes) % &qq;
                circ.sim_load_reg_bytes_shot(&a[..n], &x.to_bytes_le(), shot);
                xs.push(x);
            }
            mod_mul_pow2_pm_secp256k1(&mut circ, &a, k);
            mod_mul_pow2_pm_secp256k1_reverse(&mut circ, &a, k);
            let (sim, det) = circ.destroy_sim(a.into_iter().collect());
            let mut bad = 0;
            for (shot, x) in xs.iter().enumerate() {
                let mut got = BigUint::from(0u32);
                for i in 0..n {
                    if sim.read_bit_shot(&det[i], shot) == 1 {
                        got |= BigUint::from(1u32) << i;
                    }
                }
                if &got != x {
                    bad += 1;
                }
            }
            eprintln!("mul_pow2 roundtrip k={k}: bad={bad}/64");
            assert_eq!(bad, 0, "k={k}: roundtrip not identity, {bad}/64");
        }
    }
}
