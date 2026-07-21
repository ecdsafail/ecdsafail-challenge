//! Comparison primitives for secp256k1-sized registers: `>= const`,
//! `>= p` / `>= p/2`, physical and phase-corrected variants, built on the
//! Khattar-Gidney `compare_geq_theorem3` core. Extracted from `poc_arith`.

use crate::point_add::trailmix_port::circuit::{BorrowedQReg, Circuit, QReg};

/// Compare a >= val (classical constant), XOR result into flag.
/// Selfwire ripple-borrow: computes the carry-chain of
/// a + ~val + 1 using 2 transient ancillas.  `carry_out` = 1 iff
/// a >= val.  For each bit, the "b bit" is classical ~val[i].
pub fn compare_geq_const(circ: &mut Circuit, a: &[QReg], val: &[u8], flag: &QReg) {
    // Theorem 3 (Vandaele 2026): classical-quantum compare with 1 dirty
    // ancilla (polylog peak, log2(n)+2 for any constant). No n-qubit
    // temp register for the constant, so the peak stays logarithmic.
    let n = a.len();
    if n == 0 {
        circ.x(flag);
        return;
    }
    crate::point_add::trailmix_port::arith::khattar_gidney::compare_geq_theorem3(circ, a, val, flag);
}

/// Inline compare a >= `secp256k1_p` for 257-bit register.
/// Uses the exact identity for secp256k1:
///
///   p = 2^256 - R, where R = 2^32 + 977
///   x >= p  <=>  x + R overflows 256 bits
///
/// for `x = a[0..256)`. We realize the overflow predicate directly as:
///
///   a[256]
///   OR
///   (AND bits[33..255] AND
///      (a[32] OR (AND bits[10..31] AND (a[0..10) >= 47))))
///
/// where the low threshold comes from `2^10 - 977 = 47`.
///
/// The long ANDs use the Khattar-Gidney prefix decomposition rather
/// than the old `mcx_clean_k` recursion, which keeps the ancilla budget
/// small while making these all-ones checks linear-time.
pub fn compare_geq_p_secp256k1(circ: &mut Circuit, a: &[QReg], flag: &QReg) {
    compare_geq_p_secp256k1_inner(circ, a, BorrowedQReg::Borrowed(flag));
}

/// Consume variant: takes `flag` by value, frees it at last gate-touch
/// (before the uncompute pass allocates `kg_and_anc` ancillae, which
/// would advance `last_alloc_op_idx` past flag's last touch and trip the
/// strict-dealloc retention check).
pub fn compare_geq_p_secp256k1_consume(circ: &mut Circuit, a: &[QReg], flag: QReg) {
    compare_geq_p_secp256k1_inner(circ, a, BorrowedQReg::Owned(flag));
}

fn compare_geq_p_secp256k1_inner(circ: &mut Circuit, a: &[QReg], flag: BorrowedQReg<'_>) {
    assert!(a.len() == 257);
    use crate::point_add::trailmix_port::arith::khattar_gidney::xor_and_of_khattar_gidney;

    let low4_all_ones = circ.alloc_qreg("cmp_p_low4_all_ones");
    xor_and_of_khattar_gidney(circ, &a[..4], &low4_all_ones);

    let low_tail_or = circ.alloc_qreg("cmp_p_low_tail_or");
    circ.cx(&a[4], &low_tail_or);
    circ.cx(&low4_all_ones, &low_tail_or);
    circ.ccx(&a[4], &low4_all_ones, &low_tail_or);

    let low6_ge = circ.alloc_qreg("cmp_p_low6_ge");
    circ.ccx(&a[5], &low_tail_or, &low6_ge);

    let hi_or_67 = circ.alloc_qreg("cmp_p_hi_or_67");
    circ.cx(&a[6], &hi_or_67);
    circ.cx(&a[7], &hi_or_67);
    circ.ccx(&a[6], &a[7], &hi_or_67);

    let hi_or_89 = circ.alloc_qreg("cmp_p_hi_or_89");
    circ.cx(&a[8], &hi_or_89);
    circ.cx(&a[9], &hi_or_89);
    circ.ccx(&a[8], &a[9], &hi_or_89);

    let hi4_nonzero = circ.alloc_qreg("cmp_p_hi4_nonzero");
    circ.cx(&hi_or_67, &hi4_nonzero);
    circ.cx(&hi_or_89, &hi4_nonzero);
    circ.ccx(&hi_or_67, &hi_or_89, &hi4_nonzero);

    let low10_ge = circ.alloc_qreg("cmp_p_low10_ge");
    circ.cx(&low6_ge, &low10_ge);
    circ.cx(&hi4_nonzero, &low10_ge);
    circ.ccx(&low6_ge, &hi4_nonzero, &low10_ge);

    let mid_all_ones = circ.alloc_qreg("cmp_p_mid_all_ones");
    xor_and_of_khattar_gidney(circ, &a[10..32], &mid_all_ones);

    let mid_and_low = circ.alloc_qreg("cmp_p_mid_and_low");
    circ.ccx(&mid_all_ones, &low10_ge, &mid_and_low);

    let tail_or = circ.alloc_qreg("cmp_p_tail_or");
    circ.cx(&a[32], &tail_or);
    circ.cx(&mid_and_low, &tail_or);
    circ.ccx(&a[32], &mid_and_low, &tail_or);

    let high_all_ones = circ.alloc_qreg("cmp_p_high_all_ones");
    xor_and_of_khattar_gidney(circ, &a[33..256], &high_all_ones);

    let high_and_tail = circ.alloc_qreg("cmp_p_high_and_tail");
    circ.ccx(&high_all_ones, &tail_or, &high_and_tail);

    circ.cx(&a[256], &flag);
    circ.cx(&high_and_tail, &flag);
    circ.ccx(&a[256], &high_and_tail, &flag);
    // OWNED-flag consume path: free flag immediately after its last
    // gate-touch above. The uncompute pass below allocates kg_and_anc
    // inside xor_and_of_khattar_gidney; deferring the free until after
    // those allocs would trip the strict-dealloc retention check.
    if let BorrowedQReg::Owned(f) = flag {
        circ.zero_and_free(f);
    }

    // === MBU uncompute ===
    //
    // The 9 internal AND/OR-tree CCX pairs collapse to 1 CCX (forward
    // compute) + HMR + cz_if_bit (uncompute), saving 1 CCX per pair.
    //
    // For pure-AND targets (low6_ge, mid_and_low, high_and_tail), the
    // forward was a single CCX so the uncompute is straightforward:
    // declare_and_of(target, ctrl_a, ctrl_b); HMR(target); cz_if_bit.
    //
    // For OR-pattern targets (low_tail_or, hi_or_67, hi_or_89,
    // hi4_nonzero, low10_ge, tail_or), the forward was
    // `cx(p, t); cx(q, t); ccx(p, q, t)` giving t = p XOR q XOR (p AND q)
    // = p OR q. The uncompute peels off the linear part first, leaving
    // t = p AND q in the simulator (because (p OR q) XOR p XOR q = p AND q
    // in F2), then HMR + cz_if_bit discharges the AND obligation.
    mbu_uncompute_and(circ, high_and_tail, &high_all_ones, &tail_or);

    xor_and_of_khattar_gidney(circ, &a[33..256], &high_all_ones);
    drop(high_all_ones);

    mbu_uncompute_or(circ, tail_or, &a[32], &mid_and_low);

    mbu_uncompute_and(circ, mid_and_low, &mid_all_ones, &low10_ge);

    xor_and_of_khattar_gidney(circ, &a[10..32], &mid_all_ones);
    drop(mid_all_ones);

    mbu_uncompute_or(circ, low10_ge, &low6_ge, &hi4_nonzero);

    mbu_uncompute_or(circ, hi4_nonzero, &hi_or_67, &hi_or_89);

    mbu_uncompute_or(circ, hi_or_89, &a[8], &a[9]);

    mbu_uncompute_or(circ, hi_or_67, &a[6], &a[7]);

    mbu_uncompute_and(circ, low6_ge, &a[5], &low_tail_or);

    mbu_uncompute_or(circ, low_tail_or, &a[4], &low4_all_ones);

    xor_and_of_khattar_gidney(circ, &a[..4], &low4_all_ones);
    drop(low4_all_ones);
}

/// MBU uncompute of a pure-AND target: `target = p AND q` is replaced
/// by `HMR(target, bit); cz_if_bit(p, q, bit)` instead of the
/// reverse `ccx(p, q, target)`. Saves 1 CCX per call.
///
/// `target` enters with sim value `p AND q` (the forward CCX put it
/// there) and exits as |0> after HMR. `p` and `q` must NOT have
/// been re-versioned between the forward CCX and this call —
/// `declare_and_of` verifies the equality across all 64 sim shots.
fn mbu_uncompute_and(circ: &mut Circuit, target: QReg, p: &QReg, q: &QReg) {
    circ.declare_and_of(&target, p, q);
    let bit = circ.alloc_bit();
    circ.hmr(&target, bit);
    circ.cz_if_bit(p, q, bit);
    circ.free_bit(bit);
    drop(target);
}

/// MBU uncompute of an OR target whose forward was
/// `cx(p, target); cx(q, target); ccx(p, q, target)` (= `target = p OR q`).
///
/// Replaces the reverse `ccx(p, q, target); cx(q, target); cx(p, target)`
/// with `cx(q, target); cx(p, target); HMR(target); cz_if_bit(p, q, bit)`.
/// After the two CXs, `target = (p OR q) XOR q XOR p = p AND q` in F2,
/// matching the same MBU AND-discharge pattern. Saves 1 CCX per call.
fn mbu_uncompute_or(circ: &mut Circuit, target: QReg, p: &QReg, q: &QReg) {
    // Strip the linear part: target = p OR q -> target XOR q -> XOR p = p AND q.
    circ.cx(q, &target);
    circ.cx(p, &target);
    circ.declare_and_of(&target, p, q);
    let bit = circ.alloc_bit();
    circ.hmr(&target, bit);
    circ.cz_if_bit(p, q, bit);
    circ.free_bit(bit);
    drop(target);
}

/// Inline compare a >= ceil(p/2) for a 256-bit register.
///
/// With
///
///   ceil(p/2) = 2^255 - 2^31 - 488 = 2^255 - 2^31 - (2^9 - 24),
///
/// the predicate is:
///
///   a[255]
///   OR
///   (AND bits[32..254] AND
///      (a[31] OR (AND bits[9..30] AND (a[0..9) >= 24))))
pub fn compare_geq_half_p_secp256k1(circ: &mut Circuit, a: &[QReg], flag: &QReg) {
    compare_geq_half_p_secp256k1_inner(circ, a, BorrowedQReg::Borrowed(flag));
}

/// Consume variant: takes `flag` by value, frees it at last gate-touch
/// (before uncompute allocs, same reasoning as the consume version of
/// `compare_geq_p_secp256k1`).
pub fn compare_geq_half_p_secp256k1_consume(circ: &mut Circuit, a: &[QReg], flag: QReg) {
    compare_geq_half_p_secp256k1_inner(circ, a, BorrowedQReg::Owned(flag));
}

fn compare_geq_half_p_secp256k1_inner(circ: &mut Circuit, a: &[QReg], flag: BorrowedQReg<'_>) {
    assert!(a.len() == 256);
    use crate::point_add::trailmix_port::arith::khattar_gidney::xor_and_of_khattar_gidney;

    let hi_or_56 = circ.alloc_qreg("cmp_half_hi_or_56");
    circ.cx(&a[5], &hi_or_56);
    circ.cx(&a[6], &hi_or_56);
    circ.ccx(&a[5], &a[6], &hi_or_56);

    let hi_or_78 = circ.alloc_qreg("cmp_half_hi_or_78");
    circ.cx(&a[7], &hi_or_78);
    circ.cx(&a[8], &hi_or_78);
    circ.ccx(&a[7], &a[8], &hi_or_78);

    let hi4_nonzero = circ.alloc_qreg("cmp_half_hi4_nonzero");
    circ.cx(&hi_or_56, &hi4_nonzero);
    circ.cx(&hi_or_78, &hi4_nonzero);
    circ.ccx(&hi_or_56, &hi_or_78, &hi4_nonzero);

    let low5_ge24 = circ.alloc_qreg("cmp_half_low5_ge24");
    circ.ccx(&a[4], &a[3], &low5_ge24);

    let low9_ge24 = circ.alloc_qreg("cmp_half_low9_ge24");
    circ.cx(&hi4_nonzero, &low9_ge24);
    circ.cx(&low5_ge24, &low9_ge24);
    circ.ccx(&hi4_nonzero, &low5_ge24, &low9_ge24);

    let mid_all_ones = circ.alloc_qreg("cmp_half_mid_all_ones");
    xor_and_of_khattar_gidney(circ, &a[9..31], &mid_all_ones);

    let low_branch = circ.alloc_qreg("cmp_half_low_branch");
    circ.ccx(&mid_all_ones, &low9_ge24, &low_branch);

    let tail_or = circ.alloc_qreg("cmp_half_tail_or");
    circ.cx(&a[31], &tail_or);
    circ.cx(&low_branch, &tail_or);
    circ.ccx(&a[31], &low_branch, &tail_or);

    let high_all_ones = circ.alloc_qreg("cmp_half_high_all_ones");
    xor_and_of_khattar_gidney(circ, &a[32..255], &high_all_ones);

    let high_and_tail = circ.alloc_qreg("cmp_half_high_and_tail");
    circ.ccx(&high_all_ones, &tail_or, &high_and_tail);

    circ.cx(&a[255], &flag);
    circ.cx(&high_and_tail, &flag);
    circ.ccx(&a[255], &high_and_tail, &flag);
    // OWNED-flag consume path: free flag at its last touch, before
    // uncompute kg_and_anc allocs.
    if let BorrowedQReg::Owned(f) = flag {
        circ.zero_and_free(f);
    }

    // === MBU uncompute === (same pattern as compare_geq_p_secp256k1_inner;
    // see mbu_uncompute_and / mbu_uncompute_or for the algebra.)
    mbu_uncompute_and(circ, high_and_tail, &high_all_ones, &tail_or);

    xor_and_of_khattar_gidney(circ, &a[32..255], &high_all_ones);
    drop(high_all_ones);

    mbu_uncompute_or(circ, tail_or, &a[31], &low_branch);

    mbu_uncompute_and(circ, low_branch, &mid_all_ones, &low9_ge24);

    xor_and_of_khattar_gidney(circ, &a[9..31], &mid_all_ones);
    drop(mid_all_ones);

    mbu_uncompute_or(circ, low9_ge24, &hi4_nonzero, &low5_ge24);

    mbu_uncompute_and(circ, low5_ge24, &a[4], &a[3]);

    mbu_uncompute_or(circ, hi4_nonzero, &hi_or_56, &hi_or_78);

    mbu_uncompute_or(circ, hi_or_78, &a[7], &a[8]);

    mbu_uncompute_or(circ, hi_or_56, &a[5], &a[6]);
}

/// Compute-middle-uncompute variant: forward MAJ, set flag = (a>=b), invoke
/// `body` (which can use `flag` but must NOT touch `a` or `b` -- they are
/// scrambled during the middle), clear flag via the same cx(carry, flag),
/// then backward UMA. Cost ~2n CCX = HALF of `compute_compare_geq` + use +
/// `uncompute_compare_geq` (which would be ~4n).
///
/// Caller convention: on entry `flag` is |0>. Inside `body`, `flag` holds
/// (a >= b). On exit from `compare_geq_physical_middle`, `flag` is back to
/// |0> and `a`, `b` are restored exactly.
pub fn compare_geq_physical_middle<F: FnOnce(&mut Circuit, &QReg)>(
    circ: &mut Circuit,
    a: &[QReg],
    b: &[QReg],
    flag: &QReg,
    body: F,
) {
    let na = a.len();
    let nb = b.len();
    let n = na.max(nb);
    if n == 0 {
        circ.x(flag);
        body(circ, flag);
        circ.x(flag);
        return;
    }
    let prev = circ.push_section("cmp_middle");

    let carry = circ.alloc_qreg("carry");
    circ.x(&carry); // initial carry = 1

    let mut ext_a: Vec<QReg> = Vec::new();
    let mut ext_b: Vec<QReg> = Vec::new();

    // Forward MAJ pass.
    for i in 0..n {
        if i >= na {
            ext_a.push(circ.alloc_qreg("q"));
        }
        if i >= nb {
            ext_b.push(circ.alloc_qreg("q"));
        }
        let ai: &QReg = if i < na { &a[i] } else { &ext_a[i - na] };
        let bi: &QReg = if i < nb { &b[i] } else { &ext_b[i - nb] };
        circ.x(bi);
        circ.cx(&carry, bi);
        circ.cx(&carry, ai);
        circ.ccx(ai, bi, &carry);
    }

    circ.cx(&carry, flag); // flag = (a >= b)

    body(circ, flag); // callback uses flag

    circ.cx(&carry, flag); // XOR-clean flag back to |0>

    // Backward UMA pass.
    for i in (0..n).rev() {
        let ai: &QReg = if i < na { &a[i] } else { &ext_a[i - na] };
        let bi: &QReg = if i < nb { &b[i] } else { &ext_b[i - nb] };
        circ.ccx(ai, bi, &carry);
        circ.cx(&carry, ai);
        circ.cx(&carry, bi);
        circ.x(bi);
    }

    circ.x(&carry);
    circ.zero_and_free(carry);
    for q in ext_a {
        circ.zero_and_free(q);
    }
    for q in ext_b {
        circ.zero_and_free(q);
    }
    circ.pop_section(&prev);
}

/// Gidney measure-uncompute variant of [`compare_geq_physical_middle`].
///
/// Same contract: forward sets `flag = (a >= b)`, `body` may use `flag` but
/// must NOT touch `a`/`b` (scrambled during the middle), `flag` returns to
/// |0>, and `a`/`b` are restored exactly. The difference is the uncompute:
/// the `n` carries of the `a + ~b + 1` ripple are held in `n+1` ancillae and
/// each carry AND is erased by an X-basis measurement + a `CZ` on its two
/// (still-alive) AND inputs (Gidney 2018 measure-and-fixup, arXiv:1709.06648
/// Fig.3) rather than a Toffoli. Cost: `n` Toffoli vs `2n`; peak `+(n+1)`.
/// Use where the ancilla headroom exists (e.g. the GCD comparator).
pub fn compare_geq_gidney_middle<F: FnOnce(&mut Circuit, &QReg)>(
    circ: &mut Circuit,
    a: &[QReg],
    b: &[QReg],
    flag: &QReg,
    body: F,
) {
    let na = a.len();
    let nb = b.len();
    let n = na.max(nb);
    if n == 0 {
        circ.x(flag);
        body(circ, flag);
        circ.x(flag);
        return;
    }
    let prev = circ.push_section("cmp_gidney_middle");

    let mut ext_a: Vec<QReg> = Vec::new();
    let mut ext_b: Vec<QReg> = Vec::new();
    // Carry chain cy[0..=n]; cy[0] = carry-in = 1 (the +1 of a + ~b + 1).
    let mut cy: Vec<Option<QReg>> = Vec::with_capacity(n + 1);
    let c0 = circ.alloc_qreg("cmpg_cy");
    circ.x(&c0);
    cy.push(Some(c0));

    // Forward: compute each carry via a Gidney AND (held). Scrambles a[i],
    // b[i] into the AND inputs ta = a_i^c_i, tb = ~b_i^c_i.
    for i in 0..n {
        if i >= na {
            ext_a.push(circ.alloc_qreg("q"));
        }
        if i >= nb {
            ext_b.push(circ.alloc_qreg("q"));
        }
        let ai: &QReg = if i < na { &a[i] } else { &ext_a[i - na] };
        let bi: &QReg = if i < nb { &b[i] } else { &ext_b[i - nb] };
        let next = circ.alloc_qreg("cmpg_cy");
        let ci = cy[i].as_ref().unwrap();
        circ.x(bi); // bi = ~b_i
        circ.cx(ci, bi); // bi = ~b_i ^ c_i   (= tb)
        circ.cx(ci, ai); // ai = a_i  ^ c_i   (= ta)
        circ.ccx(ai, bi, &next); // next = ta & tb   (Gidney AND, 1 Toffoli)
        circ.cx(ci, &next); // next = c_i ^ (ta&tb) = c_{i+1}
        cy.push(Some(next));
    }

    circ.cx(cy[n].as_ref().unwrap(), flag); // flag = c_n = (a >= b)
    body(circ, flag);
    circ.cx(cy[n].as_ref().unwrap(), flag); // clean flag back to |0>

    // Reverse: measure-uncompute each carry AND, then restore a[i], b[i].
    for i in (0..n).rev() {
        let ai: &QReg = if i < na { &a[i] } else { &ext_a[i - na] };
        let bi: &QReg = if i < nb { &b[i] } else { &ext_b[i - nb] };
        let next = cy[i + 1].take().unwrap();
        // Undo the `c_i XOR`: next goes from c_{i+1} back to ta & tb.
        circ.cx(cy[i].as_ref().unwrap(), &next);
        // Measure-and-fixup AND erasure: HMR(next), then CZ(ta, tb).
        let mut g = circ.hmr_ghost(&next);
        circ.zero_and_free(next);
        circ.ghost_xor_cz(&mut g, ai, bi);
        circ.close_ghost(g);
        // Restore inputs: ta ^ c_i = a_i, tb ^ c_i = ~b_i, then ~b_i -> b_i.
        circ.cx(cy[i].as_ref().unwrap(), ai);
        circ.cx(cy[i].as_ref().unwrap(), bi);
        circ.x(bi);
    }

    let c0 = cy[0].take().unwrap();
    circ.x(&c0); // carry-in 1 -> 0
    circ.zero_and_free(c0);
    for q in ext_a {
        circ.zero_and_free(q);
    }
    for q in ext_b {
        circ.zero_and_free(q);
    }
    circ.pop_section(&prev);
}

/// MBU variant of `compare_lt_phase_correction`. Takes `q_to_hmr`
/// (the overflow bit whose phase-kick must be discharged), HMRs
/// it itself, and uses `declare_identity` so the tracker can follow
/// the obligation -> discharge match structurally.
///
/// Physically equivalent to: `caller.hmr(q_to_hmr`, bit);
/// `compare_lt_phase_correction(a`, b, bit). Same gate count
/// (+ 2 single-qubit X's on the compare-carry ancilla, and the
/// no-op `declare_identity`).
///
/// IDENTITY (proved at call site): `val(q_to_hmr)` = 1[a < b].
/// Caller must ensure this holds. For `rfold_mbu` callers with
/// a, b, `q_to_hmr` = overflow-of-add-pre-rfold: proved by the
/// case analysis in `rfold_mbu.rs`'s module header.
pub fn compare_lt_phase_correction_mbu(
    circ: &mut Circuit,
    a: &[QReg],
    b: &[QReg],
    q_to_hmr: &QReg,
) {
    let na = a.len();
    let nb = b.len();
    let n = na.max(nb);
    if n == 0 {
        // 0-width compare: 1[a<b] = 0. No phase. Still need to
        // HMR q_to_hmr -- but only if it's nonzero. If caller
        // insists on calling with n=0, they want trivial handling.
        let bit = circ.alloc_bit();
        circ.hmr(q_to_hmr, bit);
        circ.free_bit(bit);
        return;
    }

    let carry = circ.alloc_qreg("carry");
    circ.x(&carry);

    let mut ext_a: Vec<QReg> = Vec::new();
    let mut ext_b: Vec<QReg> = Vec::new();
    for i in 0..n {
        if i >= na {
            ext_a.push(circ.alloc_qreg("q"));
        }
        if i >= nb {
            ext_b.push(circ.alloc_qreg("q"));
        }
        let ai: &QReg = if i < na { &a[i] } else { &ext_a[i - na] };
        let bi: &QReg = if i < nb { &b[i] } else { &ext_b[i - nb] };
        circ.x(bi);
        circ.cx(&carry, bi);
        circ.cx(&carry, ai);
        circ.ccx(ai, bi, &carry);
    }
    // carry = 1[a >= b]. Flip to 1[a < b].
    circ.x(&carry);

    // Identity check is now inside Circuit::declare_identity.
    circ.declare_identity(q_to_hmr, &carry);

    let bit = circ.alloc_bit();
    circ.hmr(q_to_hmr, bit);
    circ.z_if_bit(&carry, bit);
    circ.free_bit(bit);

    // Restore carry to 1[a >= b] for backward UMA.
    circ.x(&carry);
    for i in (0..n).rev() {
        let ai: &QReg = if i < na { &a[i] } else { &ext_a[i - na] };
        let bi: &QReg = if i < nb { &b[i] } else { &ext_b[i - nb] };
        circ.ccx(ai, bi, &carry);
        circ.cx(&carry, ai);
        circ.cx(&carry, bi);
        circ.x(bi);
    }

    circ.x(&carry);
    // After X+MAJ+UMA+X, carry is physically |0>. Tell the tracker.
    drop(carry);
    drop(ext_a);
    drop(ext_b);
}

/// MBU variant of `controlled_compare_lt_phase_correction`.
/// Identity: `val(q_to_hmr)` = ctrl AND 1[a < b].
///
/// Needs an extra ancilla `match_q` to materialize the AND
/// (since `declare_identity` can only assert equality between
/// two qubits, not a logical expression).
pub fn controlled_compare_lt_phase_correction_mbu(
    circ: &mut Circuit,
    ctrl: &QReg,
    a: &[QReg],
    b: &[QReg],
    q_to_hmr: &QReg,
) {
    let na = a.len();
    let nb = b.len();
    let n = na.max(nb);
    if n == 0 {
        let bit = circ.alloc_bit();
        circ.hmr(q_to_hmr, bit);
        circ.free_bit(bit);
        return;
    }

    // COPY ctrl into a fresh ancilla BEFORE the forward MAJ runs,
    // because `ctrl` may alias a bit of `b` (e.g. when the outer
    // call passes the same register as addend and source-of-ctrl
    // -- see horner-squaring lsq = lambda*lambda). The MAJ's inner
    // loop does x(b[i]); cx(carry, b[i]); ... which would corrupt
    // the original ctrl qubit. We operate on ctrl_copy instead.
    let ctrl_copy = circ.alloc_qreg("ctrl_copy");
    circ.cx(ctrl, &ctrl_copy);

    let carry = circ.alloc_qreg("carry");
    circ.x(&carry);

    let mut ext_a: Vec<QReg> = Vec::new();
    let mut ext_b: Vec<QReg> = Vec::new();
    for i in 0..n {
        if i >= na {
            ext_a.push(circ.alloc_qreg("q"));
        }
        if i >= nb {
            ext_b.push(circ.alloc_qreg("q"));
        }
        let ai: &QReg = if i < na { &a[i] } else { &ext_a[i - na] };
        let bi: &QReg = if i < nb { &b[i] } else { &ext_b[i - nb] };
        circ.x(bi);
        circ.cx(&carry, bi);
        circ.cx(&carry, ai);
        circ.ccx(ai, bi, &carry);
    }
    // carry = 1[a >= b]. Flip to 1[a < b].
    circ.x(&carry);

    // Pin q_to_hmr's tracked value to ctrl_copy AND carry so the
    // upcoming cz_if_bit(ctrl_copy, carry, bit) structurally matches
    // the HMR obligation.
    circ.declare_and_of(q_to_hmr, &ctrl_copy, &carry);

    let bit = circ.alloc_bit();
    circ.hmr(q_to_hmr, bit);
    circ.cz_if_bit(&ctrl_copy, &carry, bit);
    circ.free_bit(bit);

    circ.x(&carry);
    for i in (0..n).rev() {
        let ai: &QReg = if i < na { &a[i] } else { &ext_a[i - na] };
        let bi: &QReg = if i < nb { &b[i] } else { &ext_b[i - nb] };
        circ.ccx(ai, bi, &carry);
        circ.cx(&carry, ai);
        circ.cx(&carry, bi);
        circ.x(bi);
    }

    circ.x(&carry);
    drop(carry);
    drop(ext_a);
    drop(ext_b);

    // Uncompute ctrl_copy now that ctrl's original value has been
    // restored by the backward MAJ.
    circ.cx(ctrl, &ctrl_copy);
    drop(ctrl_copy);
}

/// Top-K truncated variant of `controlled_compare_lt_phase_correction_mbu`.
/// Builds `1[a < b]` from only the top `k` bits (requires `a.len()==b.len()`):
/// the borrow chain starts at bit `n-k` with borrow-in 0, so the result equals
/// the true `1[a<b]` unless `a,b` agree on their top `k` bits (≈ 2^-k for
/// uniform inputs). That is a Shor-tolerant *phase* approximation; cost is
/// ~2k Toffoli instead of 2n. The reverse (X-sandwich in
/// `controlled_mod_sub_rfold_mbu`) re-runs this same truncated comparator, so
/// forward/reverse stay consistent. Only the per-shot phase is approximate;
/// the comparator's carry ancilla is computed-then-uncomputed exactly.
pub fn controlled_compare_lt_phase_correction_mbu_topk(
    circ: &mut Circuit,
    ctrl: &QReg,
    a: &[QReg],
    b: &[QReg],
    q_to_hmr: &QReg,
    k: usize,
) {
    let n = a.len();
    assert_eq!(n, b.len(), "topk comparator requires equal a/b lengths");
    if n == 0 || k == 0 {
        let bit = circ.alloc_bit();
        circ.hmr(q_to_hmr, bit);
        circ.free_bit(bit);
        return;
    }
    let k = k.min(n);
    let lo = n - k;

    let ctrl_copy = circ.alloc_qreg("ctrl_copy");
    circ.cx(ctrl, &ctrl_copy);

    let carry = circ.alloc_qreg("carry");
    circ.x(&carry);

    // Forward borrow chain over the top k bits only (assume no borrow from
    // below bit `lo`). carry := 1[a[lo..] >= b[lo..]].
    for i in lo..n {
        let ai = &a[i];
        let bi = &b[i];
        circ.x(bi);
        circ.cx(&carry, bi);
        circ.cx(&carry, ai);
        circ.ccx(ai, bi, &carry);
    }
    circ.x(&carry); // carry := 1[a[lo..] < b[lo..]] ≈ 1[a < b].

    circ.declare_and_of(q_to_hmr, &ctrl_copy, &carry);
    let bit = circ.alloc_bit();
    circ.hmr(q_to_hmr, bit);
    circ.cz_if_bit(&ctrl_copy, &carry, bit);
    circ.free_bit(bit);

    circ.x(&carry);
    for i in (lo..n).rev() {
        let ai = &a[i];
        let bi = &b[i];
        circ.ccx(ai, bi, &carry);
        circ.cx(&carry, ai);
        circ.cx(&carry, bi);
        circ.x(bi);
    }
    circ.x(&carry);
    drop(carry);

    circ.cx(ctrl, &ctrl_copy);
    drop(ctrl_copy);
}

/// MBU variant of `compare_geq_p_secp256k1_phase_correction`.
/// HMRs `q_to_hmr` after building `compare_carry` = 1[a >= p],
/// with `declare_identity` so the tracker sees the match.
///
/// IDENTITY (caller proves): `val(q_to_hmr)` = 1[a >= `p_secp256k1`].
/// Consumes and frees `q_to_hmr` internally.
pub fn compare_geq_p_secp256k1_phase_correction_mbu(
    circ: &mut Circuit,
    a: &[QReg],
    q_to_hmr: QReg,
) {
    assert!(a.len() == 257);
    // Use the specialized secp256k1 compare here, not the generic
    // theorem-3 builder. The specialized comparator already has a
    // proven 257-bit shape for p = 2^256 - 2^32 - 977 and is
    // self-inverse on its output qubit, which is exactly what this
    // MBU wrapper needs.
    let carry = circ.alloc_qreg("carry");
    compare_geq_p_secp256k1(circ, a, &carry);

    // IDENTITY: val(q_to_hmr) == val(carry) = 1[a >= p].
    circ.declare_identity(&q_to_hmr, &carry);
    let bit = circ.alloc_bit();
    circ.hmr(&q_to_hmr, bit);
    circ.z_if_bit(&carry, bit);
    circ.free_bit(bit);
    // Free q_to_hmr at its last gate-touch (hmr above) before the uncompute
    // section allocates kg_and_anc ancillae.
    drop(q_to_hmr);

    // Uncompute carry: after z_if_bit, carry = 1[a >= p] (value unchanged by
    // phase gate). We cannot pass carry directly to the second compare call
    // because that call allocates ancillae BEFORE it first touches carry,
    // which would advance last_alloc_op_idx past carry's last gate-touch
    // (z_if_bit) and trigger a wasteful-retention panic.
    //
    // Instead: build a fresh carry2 = 1[a >= p] via a second forward compare,
    // XOR carry2 into carry (zeroing carry since carry == carry2), free carry
    // at that last-touch, then uncompute carry2 via the consume variant.
    // carry is freed before carry2's compare inner allocs, carry2 is freed
    // at its last touch inside compare_geq_p_secp256k1_inner.
    let carry2 = circ.alloc_qreg("carry2");
    compare_geq_p_secp256k1(circ, a, &carry2);
    circ.cx(&carry2, &carry); // carry ^= carry2 = 0 (carry == carry2 = 1[a>=p])
    drop(carry); // free at last touch (cx above), before carry2 uncompute allocs
    compare_geq_p_secp256k1_consume(circ, a, carry2);
}

/// MBU variant of `compare_geq_half_p_secp256k1`.
/// HMRs `q_to_hmr` after building carry = 1[a >= ceil(p/2)].
/// Consumes and frees `q_to_hmr` internally.
pub fn compare_geq_half_p_secp256k1_phase_correction_mbu(
    circ: &mut Circuit,
    a: &[QReg],
    q_to_hmr: QReg,
) {
    let carry = circ.alloc_qreg("carry");
    compare_geq_half_p_secp256k1(circ, a, &carry);

    circ.declare_identity(&q_to_hmr, &carry);
    let bit = circ.alloc_bit();
    circ.hmr(&q_to_hmr, bit);
    circ.z_if_bit(&carry, bit);
    circ.free_bit(bit);
    // Free q_to_hmr at its last gate-touch (hmr above) before the uncompute
    // section allocates kg_and_anc ancillae.
    drop(q_to_hmr);

    // Same carry2 pattern as compare_geq_p_secp256k1_phase_correction_mbu:
    // build a fresh carry2 to zero carry before the uncompute allocs.
    let carry2 = circ.alloc_qreg("carry2");
    compare_geq_half_p_secp256k1(circ, a, &carry2);
    circ.cx(&carry2, &carry);
    drop(carry);
    compare_geq_half_p_secp256k1_consume(circ, a, carry2);
}
