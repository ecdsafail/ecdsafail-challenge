//! (refactor) Mechanically extracted from mod.rs. No logic changes.
use super::*;

// ═══════════════════════════════════════════════════════════════════════════
//  Cuccaro ripple-carry adder
// ═══════════════════════════════════════════════════════════════════════════
//
// Operates on two n-wide qubit registers `a` (addend, unchanged) and
// `acc` (accumulator, becomes a + acc mod 2^n). Also takes:
//   * c_in: one ancilla qubit, = 0 on entry, = 0 on exit (unchanged)
//   * z   : one ancilla qubit, = 0 on entry, = carry_out ⊕ z_in on exit
//           (i.e., the output carry is XORed into z; pass a fresh 0 bit
//           to receive the high bit)
//
// Based on Cuccaro et al. 2004 (arXiv:quant-ph/0410184), Figure 3.
//
// `MAJ(x, y, w)` triple:
//     CX(w, y)        # y ← y ⊕ w
//     CX(w, x)        # x ← x ⊕ w
//     CCX(x, y, w)    # w ← w ⊕ (x·y)        w becomes MAJ(w_old, y_old, x_old)
//
// `UMA(x, y, w)` triple (undoes MAJ, leaves sum bit in y):
//     CCX(x, y, w)
//     CX(w, x)
//     CX(x, y)

pub(crate) fn maj(b: &mut B, x: QubitId, y: QubitId, w: QubitId) {
    b.cx(w, y);
    b.cx(w, x);
    b.ccx(x, y, w);
}

pub(crate) fn uma(b: &mut B, x: QubitId, y: QubitId, w: QubitId) {
    b.ccx(x, y, w);
    b.cx(w, x);
    b.cx(x, y);
}

/// Fast Cuccaro add using carry ancillae + measurement-based UMA.
/// Same interface as `cuccaro_add` but uses n-1 carry ancillae so the
/// UMA sweep costs 0 Toffoli (measurement only). NOT emit_inverse-safe.
pub(crate) fn cuccaro_add_fast(b: &mut B, a: &[QubitId], acc: &[QubitId], c_in: QubitId) {
    let n = a.len();
    assert_eq!(n, acc.len());
    if n == 0 {
        return;
    }
    if n == 1 {
        b.cx(c_in, acc[0]);
        b.cx(a[0], acc[0]);
        return;
    }

    let carries = b.alloc_qubits(n - 1);

    // Forward MAJ sweep with carry ancillae.
    // Step 0: MAJ(c_in, acc[0], a[0]) → carry into carries[0]
    b.cx(a[0], acc[0]);
    b.cx(a[0], c_in);
    b.ccx(c_in, acc[0], carries[0]);
    b.cx(carries[0], a[0]);
    // Steps 1..n-2: MAJ(a[i-1], acc[i], a[i]) → carry into carries[i]
    for i in 1..n - 1 {
        b.cx(a[i], acc[i]);
        b.cx(a[i], a[i - 1]);
        b.ccx(a[i - 1], acc[i], carries[i]);
        b.cx(carries[i], a[i]);
    }

    // Final sum bit (same as original cuccaro_add)
    b.cx(a[n - 2], acc[n - 1]);
    b.cx(a[n - 1], acc[n - 1]);

    // Backward UMA sweep with measurement-based carry uncompute (0 Toffoli).
    for i in (1..n - 1).rev() {
        b.cx(carries[i], a[i]);
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        b.cz_if(a[i - 1], acc[i], m);
        b.cx(a[i], a[i - 1]);
        b.cx(a[i - 1], acc[i]);
    }
    // Step 0 UMA:
    b.cx(carries[0], a[0]);
    let m0 = b.alloc_bit();
    b.hmr(carries[0], m0);
    b.cz_if(c_in, acc[0], m0);
    b.cx(a[0], c_in);
    b.cx(c_in, acc[0]);

    b.free_vec(&carries);
}

/// Carry-BORROW twin of [`cuccaro_add_fast`]: identical gate sequence, but the
/// n-1 carry qubits are BORROWED from `carry_src` (which MUST be clean |0⟩ on
/// entry and is restored to |0⟩ on exit by the measurement-uncompute) instead
/// of freshly allocated. Flat Toffoli, zero new width at the peak — the carry
/// register is hosted on already-live but idle clean ancilla (e.g. the future
/// Kaliski m_hist transcript bits m_hist[iter+1..], guaranteed |0⟩ until their
/// own iteration writes them). `carry_src.len()` must be >= n-1.
pub(crate) fn cuccaro_add_fast_borrow(
    b: &mut B,
    a: &[QubitId],
    acc: &[QubitId],
    c_in: QubitId,
    carry_src: &[QubitId],
) {
    let n = a.len();
    assert_eq!(n, acc.len());
    if n == 0 {
        return;
    }
    if n == 1 {
        b.cx(c_in, acc[0]);
        b.cx(a[0], acc[0]);
        return;
    }
    assert!(carry_src.len() >= n - 1, "borrow carry_src too short");
    let carries = &carry_src[..n - 1];

    b.cx(a[0], acc[0]);
    b.cx(a[0], c_in);
    b.ccx(c_in, acc[0], carries[0]);
    b.cx(carries[0], a[0]);
    for i in 1..n - 1 {
        b.cx(a[i], acc[i]);
        b.cx(a[i], a[i - 1]);
        b.ccx(a[i - 1], acc[i], carries[i]);
        b.cx(carries[i], a[i]);
    }

    b.cx(a[n - 2], acc[n - 1]);
    b.cx(a[n - 1], acc[n - 1]);

    for i in (1..n - 1).rev() {
        b.cx(carries[i], a[i]);
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        b.cz_if(a[i - 1], acc[i], m);
        b.cx(a[i], a[i - 1]);
        b.cx(a[i - 1], acc[i]);
    }
    b.cx(carries[0], a[0]);
    let m0 = b.alloc_bit();
    b.hmr(carries[0], m0);
    b.cz_if(c_in, acc[0], m0);
    b.cx(a[0], c_in);
    b.cx(c_in, acc[0]);
    // carries are borrowed: restored to |0> by the measurement-uncompute,
    // NOT freed (they belong to the caller's m_hist register).
}

/// In-place addition `acc += a mod 2^n` on quantum n-bit registers.
/// * `c_in` is a fresh ancilla qubit at 0 on entry and returns to 0.
/// * `a` unchanged; `acc` becomes (a + acc) mod 2^n.
/// Pure mod-2^n: the high carry is discarded (no `z` ancilla). This is
/// honestly reversible because the last MAJ/UMA pair cancel out the
/// carry information on `a[n-1]`.
pub(crate) fn cuccaro_add(b: &mut B, a: &[QubitId], acc: &[QubitId], c_in: QubitId) {
    let n = a.len();
    assert_eq!(n, acc.len());
    if n == 0 {
        return;
    }
    if n == 1 {
        // acc[0] += a[0] + c_in  mod 2 ; c_in → 0
        b.cx(c_in, acc[0]);
        b.cx(a[0], acc[0]);
        return;
    }

    // Forward MAJ sweep.
    maj(b, c_in, acc[0], a[0]);
    for i in 1..n - 1 {
        maj(b, a[i - 1], acc[i], a[i]);
    }

    // Final sum bit: sum[n-1] = acc[n-1] XOR a[n-1] XOR carry_in_to_n-1,
    // where carry_in_to_n-1 is in a[n-2] after the MAJ sweep.
    b.cx(a[n - 2], acc[n - 1]);
    b.cx(a[n - 1], acc[n - 1]);

    // Reverse UMA sweep (skips the final MAJ since we didn't do it).
    for i in (1..n - 1).rev() {
        uma(b, a[i - 1], acc[i], a[i]);
    }
    uma(b, c_in, acc[0], a[0]);
}

/// Reverse of `cuccaro_add`: performs `acc -= a mod 2^n`.
/// Implemented as the exact inverse gate sequence of `cuccaro_add`.
pub(crate) fn cuccaro_sub(b: &mut B, a: &[QubitId], acc: &[QubitId], c_in: QubitId) {
    let n = a.len();
    assert_eq!(n, acc.len());
    if n == 0 {
        return;
    }
    if n == 1 {
        // Inverse of (cx c_in acc; cx a acc) is the same two gates in reverse.
        b.cx(a[0], acc[0]);
        b.cx(c_in, acc[0]);
        return;
    }

    // Inverse of `uma(c_in, acc[0], a[0])`, then the rest of UMA sweep
    // in reverse order.
    inv_uma(b, c_in, acc[0], a[0]);
    for i in 1..n - 1 {
        inv_uma(b, a[i - 1], acc[i], a[i]);
    }

    // Inverse of the final sum writes (both CX self-inverse; reverse order).
    b.cx(a[n - 1], acc[n - 1]);
    b.cx(a[n - 2], acc[n - 1]);

    // Inverse of the forward MAJ sweep.
    for i in (1..n - 1).rev() {
        inv_maj(b, a[i - 1], acc[i], a[i]);
    }
    inv_maj(b, c_in, acc[0], a[0]);
}

pub(crate) fn inv_maj(b: &mut B, x: QubitId, y: QubitId, w: QubitId) {
    // maj = CX(w,y); CX(w,x); CCX(x,y,w)
    // inv = CCX(x,y,w); CX(w,x); CX(w,y)
    b.ccx(x, y, w);
    b.cx(w, x);
    b.cx(w, y);
}

pub(crate) fn inv_uma(b: &mut B, x: QubitId, y: QubitId, w: QubitId) {
    // uma = CCX(x,y,w); CX(w,x); CX(x,y)
    // inv = CX(x,y); CX(w,x); CCX(x,y,w)
    b.cx(x, y);
    b.cx(w, x);
    b.ccx(x, y, w);
}

// ═══════════════════════════════════════════════════════════════════════════
//  Non-modular n-bit primitives
// ═══════════════════════════════════════════════════════════════════════════

/// Fast Cuccaro sub: `acc -= a mod 2^n` with measurement UMA (0 Toffoli
/// for UMA sweep). Exact gate-level inverse of `cuccaro_add_fast`.
pub(crate) fn cuccaro_sub_fast(b: &mut B, a: &[QubitId], acc: &[QubitId], c_in: QubitId) {
    let n = a.len();
    assert_eq!(n, acc.len());
    if n == 0 {
        return;
    }
    if n == 1 {
        b.cx(a[0], acc[0]);
        b.cx(c_in, acc[0]);
        return;
    }

    let carries = b.alloc_qubits(n - 1);

    // Forward inv_UMA sweep with carry ancillae (reversed UMA from cuccaro_sub).
    // Step 0:
    b.cx(c_in, acc[0]);
    b.cx(a[0], c_in);
    b.ccx(c_in, acc[0], carries[0]);
    b.cx(carries[0], a[0]);
    // Steps 1..n-2:
    for i in 1..n - 1 {
        b.cx(a[i - 1], acc[i]);
        b.cx(a[i], a[i - 1]);
        b.ccx(a[i - 1], acc[i], carries[i]);
        b.cx(carries[i], a[i]);
    }

    // Final sum bit (reversed from cuccaro_add)
    b.cx(a[n - 1], acc[n - 1]);
    b.cx(a[n - 2], acc[n - 1]);

    // Backward inv_MAJ sweep with measurement.
    for i in (1..n - 1).rev() {
        b.cx(carries[i], a[i]);
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        b.cz_if(a[i - 1], acc[i], m);
        b.cx(a[i], a[i - 1]);
        b.cx(a[i], acc[i]);
    }
    b.cx(carries[0], a[0]);
    let m0 = b.alloc_bit();
    b.hmr(carries[0], m0);
    b.cz_if(c_in, acc[0], m0);
    b.cx(a[0], c_in);
    b.cx(a[0], acc[0]);

    b.free_vec(&carries);
}

/// Carry-BORROW twin of [`cuccaro_sub_fast`]: see [`cuccaro_add_fast_borrow`].
/// Borrows `carry_src[..n-1]` (clean |0>, restored to |0>) as the carry block.
pub(crate) fn cuccaro_sub_fast_borrow(
    b: &mut B,
    a: &[QubitId],
    acc: &[QubitId],
    c_in: QubitId,
    carry_src: &[QubitId],
) {
    let n = a.len();
    assert_eq!(n, acc.len());
    if n == 0 {
        return;
    }
    if n == 1 {
        b.cx(a[0], acc[0]);
        b.cx(c_in, acc[0]);
        return;
    }
    assert!(carry_src.len() >= n - 1, "borrow carry_src too short");
    let carries = &carry_src[..n - 1];

    b.cx(c_in, acc[0]);
    b.cx(a[0], c_in);
    b.ccx(c_in, acc[0], carries[0]);
    b.cx(carries[0], a[0]);
    for i in 1..n - 1 {
        b.cx(a[i - 1], acc[i]);
        b.cx(a[i], a[i - 1]);
        b.ccx(a[i - 1], acc[i], carries[i]);
        b.cx(carries[i], a[i]);
    }

    b.cx(a[n - 1], acc[n - 1]);
    b.cx(a[n - 2], acc[n - 1]);

    for i in (1..n - 1).rev() {
        b.cx(carries[i], a[i]);
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        b.cz_if(a[i - 1], acc[i], m);
        b.cx(a[i], a[i - 1]);
        b.cx(a[i], acc[i]);
    }
    b.cx(carries[0], a[0]);
    let m0 = b.alloc_bit();
    b.hmr(carries[0], m0);
    b.cz_if(c_in, acc[0], m0);
    b.cx(a[0], c_in);
    b.cx(a[0], acc[0]);
    // carries borrowed: restored to |0>, NOT freed.
}

/// Borrow-carry `acc += a mod 2^n`: hosts the fast-Cuccaro carry register on
/// `carry_src` (clean |0>, restored). Flat Toffoli, no new peak width.
pub(crate) fn add_nbit_qq_fast_borrow(b: &mut B, a: &[QubitId], acc: &[QubitId], carry_src: &[QubitId]) {
    assert_eq!(a.len(), acc.len());
    let c_in = b.alloc_qubit();
    cuccaro_add_fast_borrow(b, a, acc, c_in, carry_src);
    b.free(c_in);
}

/// Borrow-carry `acc -= a mod 2^n`. Companion to [`add_nbit_qq_fast_borrow`].
pub(crate) fn sub_nbit_qq_fast_borrow(b: &mut B, a: &[QubitId], acc: &[QubitId], carry_src: &[QubitId]) {
    assert_eq!(a.len(), acc.len());
    let c_in = b.alloc_qubit();
    cuccaro_sub_fast_borrow(b, a, acc, c_in, carry_src);
    b.free(c_in);
}

/// Build a width-(n-1) clean-|0> carry register for a fast-Cuccaro add/sub of
/// width n, hosting as many bits as possible on `m_future` (clean |0> bits that
/// belong to the caller and are restored to |0> on exit), and freshly
/// allocating only the shortfall. Returns (full_carry_vec, fresh_count); the
/// caller must `free` the LAST `fresh_count` entries after the add/sub. Flat
/// Toffoli vs the all-fresh path; peak width drops by `min(n-1, m_future.len())`.
pub(crate) fn borrow_carry_register(
    b: &mut B,
    n: usize,
    m_future: &[QubitId],
) -> (Vec<QubitId>, usize) {
    let need = n.saturating_sub(1);
    let borrowed = need.min(m_future.len());
    let fresh_count = need - borrowed;
    let mut carries: Vec<QubitId> = Vec::with_capacity(need);
    carries.extend_from_slice(&m_future[..borrowed]);
    for _ in 0..fresh_count {
        carries.push(b.alloc_qubit());
    }
    (carries, fresh_count)
}

/// Max fresh carry qubits we will allocate on top of the m_future borrow
/// before falling back to the slow (carry-register-FREE, 1-ancilla in-place
/// Cuccaro). Tuned so the per-step peak stays <= the 2333 shift22 floor:
/// 2333 - (binder carrier 1175 + tmp 256 + init 512 + lam 256 + slack) leaves
/// headroom for ~120 fresh carries. When the m_future pool is too small to
/// cover (n-1) with <= this many fresh bits, we use slow Cuccaro for that one
/// call (flat WIDTH, +~n Toffoli) — only the few late iters with an exhausted
/// pool pay it, keeping the global Toffoli penalty small while still capping
/// the peak at the floor.
pub(crate) fn gz_max_fresh_carries() -> usize {
    // 131 = the max fresh-carry budget that keeps the per-step peak at the 2333
    // shift22 floor (132+ lets kal_bulk_step4 borrow-fast push to 2334+).
    // Swept empirically; minimizes the slow-Cuccaro fallback Toffoli at 2333.
    std::env::var("KAL_GZ_MAX_FRESH")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(131)
}

/// `acc += a mod 2^n`: borrow the carry register from `m_future` (clean |0>);
/// if the shortfall would exceed `gz_max_fresh_carries`, fall back to the slow
/// in-place Cuccaro (no carry register at all) so the peak stays at the floor.
pub(crate) fn add_nbit_qq_fast_mfut(b: &mut B, a: &[QubitId], acc: &[QubitId], m_future: &[QubitId]) {
    assert_eq!(a.len(), acc.len());
    let n = a.len();
    let need = n.saturating_sub(1);
    let borrowed = need.min(m_future.len());
    if need - borrowed > gz_max_fresh_carries() {
        // Pool too small: slow Cuccaro (1 ancilla, no carry register).
        let c_in = b.alloc_qubit();
        cuccaro_add(b, a, acc, c_in);
        b.free(c_in);
        return;
    }
    let (carries, fresh) = borrow_carry_register(b, n, m_future);
    let c_in = b.alloc_qubit();
    cuccaro_add_fast_borrow(b, a, acc, c_in, &carries);
    b.free(c_in);
    for &q in carries[carries.len() - fresh..].iter() {
        b.free(q);
    }
}

/// `acc -= a mod 2^n` with m_future borrow + slow-Cuccaro shortfall fallback.
pub(crate) fn sub_nbit_qq_fast_mfut(b: &mut B, a: &[QubitId], acc: &[QubitId], m_future: &[QubitId]) {
    assert_eq!(a.len(), acc.len());
    let n = a.len();
    let need = n.saturating_sub(1);
    let borrowed = need.min(m_future.len());
    if need - borrowed > gz_max_fresh_carries() {
        let c_in = b.alloc_qubit();
        cuccaro_sub(b, a, acc, c_in);
        b.free(c_in);
        return;
    }
    let (carries, fresh) = borrow_carry_register(b, n, m_future);
    let c_in = b.alloc_qubit();
    cuccaro_sub_fast_borrow(b, a, acc, c_in, &carries);
    b.free(c_in);
    for &q in carries[carries.len() - fresh..].iter() {
        b.free(q);
    }
}

/// Build a width-(n-1) clean-|0> carry register from TWO already-live clean
/// pools (`m_future` first, then `extra`), freshly allocating only the
/// shortfall. Both pools must be |0> on entry; the borrowed bits are restored
/// to |0> by the caller's measurement-uncompute. Returns (carry_vec,
/// fresh_count); the caller frees the LAST `fresh_count` entries. Adds ZERO
/// peak width for every bit drawn from the pools. See [`gz_late_recover`].
pub(crate) fn borrow_carry_register_pool(
    b: &mut B,
    n: usize,
    m_future: &[QubitId],
    extra: &[QubitId],
) -> (Vec<QubitId>, usize) {
    let need = n.saturating_sub(1);
    let mut carries: Vec<QubitId> = Vec::with_capacity(need);
    let take_mf = need.min(m_future.len());
    carries.extend_from_slice(&m_future[..take_mf]);
    if carries.len() < need {
        let take_ex = (need - carries.len()).min(extra.len());
        carries.extend_from_slice(&extra[..take_ex]);
    }
    let fresh_count = need - carries.len();
    for _ in 0..fresh_count {
        carries.push(b.alloc_qubit());
    }
    (carries, fresh_count)
}

/// `acc += a mod 2^n`: borrow carries from `m_future` THEN `extra` (both clean
/// |0>, restored on exit); slow-Cuccaro fallback only if the COMBINED pool is
/// still too small. Recovers the late-iter slow-fallback Toffoli at flat peak.
pub(crate) fn add_nbit_qq_fast_mfut_pool(
    b: &mut B,
    a: &[QubitId],
    acc: &[QubitId],
    m_future: &[QubitId],
    extra: &[QubitId],
) {
    assert_eq!(a.len(), acc.len());
    let n = a.len();
    let need = n.saturating_sub(1);
    let pool = m_future.len() + extra.len();
    let borrowed = need.min(pool);
    if need - borrowed > gz_max_fresh_carries() {
        let c_in = b.alloc_qubit();
        cuccaro_add(b, a, acc, c_in);
        b.free(c_in);
        return;
    }
    let (carries, fresh) = borrow_carry_register_pool(b, n, m_future, extra);
    let c_in = b.alloc_qubit();
    cuccaro_add_fast_borrow(b, a, acc, c_in, &carries);
    b.free(c_in);
    for &q in carries[carries.len() - fresh..].iter() {
        b.free(q);
    }
}

/// `acc -= a mod 2^n` twin of [`add_nbit_qq_fast_mfut_pool`].
pub(crate) fn sub_nbit_qq_fast_mfut_pool(
    b: &mut B,
    a: &[QubitId],
    acc: &[QubitId],
    m_future: &[QubitId],
    extra: &[QubitId],
) {
    assert_eq!(a.len(), acc.len());
    let n = a.len();
    let need = n.saturating_sub(1);
    let pool = m_future.len() + extra.len();
    let borrowed = need.min(pool);
    if need - borrowed > gz_max_fresh_carries() {
        let c_in = b.alloc_qubit();
        cuccaro_sub(b, a, acc, c_in);
        b.free(c_in);
        return;
    }
    let (carries, fresh) = borrow_carry_register_pool(b, n, m_future, extra);
    let c_in = b.alloc_qubit();
    cuccaro_sub_fast_borrow(b, a, acc, c_in, &carries);
    b.free(c_in);
    for &q in carries[carries.len() - fresh..].iter() {
        b.free(q);
    }
}

/// Fast `acc += a mod 2^n` using measurement-based Cuccaro.
pub(crate) fn add_nbit_qq_fast(b: &mut B, a: &[QubitId], acc: &[QubitId]) {
    assert_eq!(a.len(), acc.len());
    let c_in = b.alloc_qubit();
    cuccaro_add_fast(b, a, acc, c_in);
    b.free(c_in);
}

/// Fast `acc -= a mod 2^n` using measurement-based Cuccaro.
pub(crate) fn sub_nbit_qq_fast(b: &mut B, a: &[QubitId], acc: &[QubitId]) {
    assert_eq!(a.len(), acc.len());
    let c_in = b.alloc_qubit();
    cuccaro_sub_fast(b, a, acc, c_in);
    b.free(c_in);
}

/// `acc += a mod 2^n`. Caller must pre-extend both slices if they want the
/// top carry absorbed into the accumulator (i.e. pass n+1-bit slices with
/// top bits 0 to get a full n+1-bit add). The carry-out beyond the slice
/// is discarded via `R` on the `z` ancilla — safe when both inputs fit
/// in n-1 bits (as in our mod-p layer where both < 2p < 2^{n+1}).
pub(crate) fn add_nbit_qq(b: &mut B, a: &[QubitId], acc: &[QubitId]) {
    assert_eq!(a.len(), acc.len());
    let c_in = b.alloc_qubit();
    cuccaro_add(b, a, acc, c_in);
    b.free(c_in);
}

pub(crate) fn sub_nbit_qq(b: &mut B, a: &[QubitId], acc: &[QubitId]) {
    assert_eq!(a.len(), acc.len());
    let c_in = b.alloc_qubit();
    cuccaro_sub(b, a, acc, c_in);
    b.free(c_in);
}

pub(crate) fn centered_restoring_trial_subtract_clean(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    q_success: QubitId,
) {
    // Trial subtract for a centered-Euclid quotient bit. Compute the borrow,
    // copy out the success bit, then undo with the arithmetic inverse instead
    // of replaying the Cuccaro subtract wrapper through emit_inverse.
    assert_eq!(u.len(), v.len());
    let top_u = b.alloc_qubit();
    let top_v = b.alloc_qubit();
    let mut u_ext = u.to_vec();
    u_ext.push(top_u);
    let mut v_ext = v.to_vec();
    v_ext.push(top_v);
    sub_nbit_qq(b, &v_ext, &u_ext);
    b.cx(top_u, q_success);
    b.x(q_success);
    add_nbit_qq(b, &v_ext, &u_ext);
    b.free(top_v);
    b.free(top_u);
}

pub(crate) fn add_nbit_const(b: &mut B, acc: &[QubitId], c: U256) {
    let n = acc.len();
    let a = load_const(b, n, c);
    add_nbit_qq(b, &a, acc);
    unload_const(b, &a, c);
}

pub(crate) fn sub_nbit_const(b: &mut B, acc: &[QubitId], c: U256) {
    let n = acc.len();
    let a = load_const(b, n, c);
    sub_nbit_qq(b, &a, acc);
    unload_const(b, &a, c);
}

pub(crate) fn csub_nbit_const(b: &mut B, acc: &[QubitId], c: U256, ctrl: QubitId) {
    // acc -= (ctrl ? c : 0). Mirror of cadd_nbit_const.
    let n = acc.len();
    let a = b.alloc_qubits(n);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, a[i]);
        }
    }
    sub_nbit_qq(b, &a, acc);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, a[i]);
        }
    }
    b.free_vec(&a);
}

pub(crate) fn cadd_nbit_const(b: &mut B, acc: &[QubitId], c: U256, ctrl: QubitId) {
    // Conditional add of constant c, controlled by qubit ctrl.
    // Trick: load c into a qubit register via CX-from-ctrl gates
    // (so the loaded value is (ctrl ? c : 0)), then unconditional add,
    // then unload.
    let n = acc.len();
    let a = b.alloc_qubits(n);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, a[i]);
        }
    }
    add_nbit_qq(b, &a, acc);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, a[i]);
        }
    }
    b.free_vec(&a);
}

pub(crate) fn csub_nbit_const_fast(b: &mut B, acc: &[QubitId], c: U256, ctrl: QubitId) {
    let n = acc.len();
    let a = b.alloc_qubits(n);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, a[i]);
        }
    }
    sub_nbit_qq_fast(b, &a, acc);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, a[i]);
        }
    }
    b.free_vec(&a);
}

/// Controlled subtract of a classical constant without materializing the
/// `ctrl ? c : 0` addend.  This is the same measurement-uncomputed ripple idea
/// as [`sub_nbit_qq_fast`], but the carry/borrow recurrence is specialized to a
/// classical bit and the external control.  It saves the n-qubit loaded-constant
/// register at Kaliski halve peaks; for sparse secp256k1 `c=2^32+977` the CCX
/// count is essentially unchanged.
pub(crate) fn csub_nbit_const_direct_fast(b: &mut B, acc: &[QubitId], c: U256, ctrl: QubitId) {
    let n = acc.len();
    let cut = kal_carrytail_count(n, kal_carrytail_sub_enabled());
    csub_nbit_const_direct_fast_cut(b, acc, c, ctrl, cut);
}

/// `csub_nbit_const_direct_fast` with an explicit carry-tail `cut` (number of
/// borrow ancillae). Used by the shift22 reduction so its window is decoupled
/// from the global Kaliski carry-tail W. `cut` MUST be identical between a
/// forward call and its reversal (phase-parity law).
pub(crate) fn csub_nbit_const_direct_fast_cut(
    b: &mut B,
    acc: &[QubitId],
    c: U256,
    ctrl: QubitId,
    cut: usize,
) {
    let n = acc.len();
    if n == 0 {
        return;
    }
    if n == 1 {
        if bit(c, 0) {
            b.cx(ctrl, acc[0]);
        }
        return;
    }

    // CARRY-TAIL truncation: compute the borrow chain only through bit `cut-1`.
    // cut == n-1 (full chain) when truncation is disabled.  Phase-parity law:
    // forward sweep, difference XOR, and reverse uncompute all use this same
    // `cut`, so they are byte-identical width.
    let cut = cut.min(n - 1);
    let borrows = b.alloc_qubits(cut);

    // Forward borrow sweep. borrow_{i+1} = majority(!acc_i, k_i, borrow_i),
    // where k_i = ctrl when c_i=1 and 0 otherwise.
    let majfold = majfold_sub_enabled();
    for i in 0..cut {
        let target = borrows[i];
        let borrow_in = if i == 0 { None } else { Some(borrows[i - 1]) };
        if bit(c, i) {
            b.x(acc[i]);
            if let Some(bi) = borrow_in {
                // MAJ(!acc[i], ctrl, bi) -> fold to 1 CCX + free CX (bi pivot):
                // maj(a,b,d)=d^(a^d)&(b^d). Value identical -> backward Hmr unchanged.
                if majfold {
                    b.cx(bi, acc[i]);
                    b.cx(bi, ctrl);
                    b.ccx(acc[i], ctrl, target);
                    b.cx(bi, target);
                    b.cx(bi, ctrl);
                    b.cx(bi, acc[i]);
                } else {
                    b.ccx(acc[i], bi, target);
                    b.ccx(ctrl, acc[i], target);
                    b.ccx(ctrl, bi, target);
                }
            } else {
                b.ccx(acc[i], ctrl, target);
            }
            b.x(acc[i]);
        } else if let Some(bi) = borrow_in {
            b.x(acc[i]);
            b.ccx(acc[i], bi, target);
            b.x(acc[i]);
        }
    }

    // Difference bits: acc_i ^= k_i ^ borrow_i.  The const XOR is always exact;
    // the borrow XOR only applies for bits whose borrow-in was computed
    // (i-1 < cut, i.e. i <= cut).
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, acc[i]);
        }
        if i > 0 && i - 1 < cut {
            b.cx(borrows[i - 1], acc[i]);
        }
    }

    // Measurement-uncompute borrows in reverse.  For subtraction the post-sum
    // identity is borrow_{i+1} = majority(acc_i_final, k_i, borrow_i).
    for i in (0..cut).rev() {
        let m = b.alloc_bit();
        b.hmr(borrows[i], m);
        let borrow_in = if i == 0 { None } else { Some(borrows[i - 1]) };
        if bit(c, i) {
            if let Some(bi) = borrow_in {
                b.cz_if(acc[i], ctrl, m);
                b.cz_if(acc[i], bi, m);
                b.cz_if(ctrl, bi, m);
            } else {
                b.cz_if(acc[i], ctrl, m);
            }
        } else if let Some(bi) = borrow_in {
            b.cz_if(acc[i], bi, m);
        }
    }

    b.free_vec(&borrows);
}

pub(crate) fn cadd_nbit_const_fast(b: &mut B, acc: &[QubitId], c: U256, ctrl: QubitId) {
    let n = acc.len();
    let a = b.alloc_qubits(n);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, a[i]);
        }
    }
    add_nbit_qq_fast(b, &a, acc);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, a[i]);
        }
    }
    b.free_vec(&a);
}

/// Controlled add of a classical constant without a loaded addend register.
/// This is the carry analogue of [`csub_nbit_const_direct_fast`].
pub(crate) fn cadd_nbit_const_direct_fast(b: &mut B, acc: &[QubitId], c: U256, ctrl: QubitId) {
    let n = acc.len();
    let cut = kal_carrytail_count_c(n, c, kal_carrytail_add_enabled());
    cadd_nbit_const_direct_fast_cut(b, acc, c, ctrl, cut);
}

/// `cadd_nbit_const_direct_fast` with an explicit carry-tail `cut`. Used by the
/// shift22 reduction so its window is decoupled from the global Kaliski carry-tail
/// W. `cut` MUST be identical between a forward call and its reversal.
pub(crate) fn cadd_nbit_const_direct_fast_cut(
    b: &mut B,
    acc: &[QubitId],
    c: U256,
    ctrl: QubitId,
    cut: usize,
) {
    let n = acc.len();
    if n == 0 {
        return;
    }
    if n == 1 {
        if bit(c, 0) {
            b.cx(ctrl, acc[0]);
        }
        return;
    }

    // CARRY-TAIL truncation: compute the carry chain only through bit `cut-1`.
    // cut == n-1 (full chain) when truncation is disabled.  Phase-parity law:
    // forward sweep, sum XOR, and reverse uncompute all use this same `cut`.
    // Constant-aware window: dense constants (e.g. c=p+1 in mod_neg/mod_double,
    // top set bit 255) must NOT be truncated — a dropped carry corrupts the high
    // sum bits, poisons a downstream sign-bit comparison, leaves a flag dirty, and
    // its free()/R injects random phase. Anchoring k0 above c's top set bit gives
    // dense constants the full chain while sparse Solinas constants keep the tight
    // truncation. This is what makes the ADD path phase-clean.
    let cut = cut.min(n - 1);
    let carries = b.alloc_qubits(cut);

    // Forward carry sweep. carry_{i+1} = majority(acc_i, k_i, carry_i).
    let majfold = majfold_add_enabled();
    for i in 0..cut {
        let target = carries[i];
        let carry_in = if i == 0 { None } else { Some(carries[i - 1]) };
        if bit(c, i) {
            if let Some(ci) = carry_in {
                // MAJ(acc[i], ctrl, ci) -> fold to 1 CCX + free CX (ci pivot):
                // maj(a,b,d)=d^(a^d)&(b^d). Value identical -> backward Hmr unchanged.
                if majfold {
                    b.cx(ci, acc[i]);
                    b.cx(ci, ctrl);
                    b.ccx(acc[i], ctrl, target);
                    b.cx(ci, target);
                    b.cx(ci, ctrl);
                    b.cx(ci, acc[i]);
                } else {
                    b.ccx(acc[i], ci, target);
                    b.ccx(ctrl, acc[i], target);
                    b.ccx(ctrl, ci, target);
                }
            } else {
                b.ccx(acc[i], ctrl, target);
            }
        } else if let Some(ci) = carry_in {
            b.ccx(acc[i], ci, target);
        }
    }

    // Sum bits: acc_i ^= k_i ^ carry_i.  The const XOR is always exact; the
    // carry XOR only applies for bits whose carry-in was computed (i <= cut).
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, acc[i]);
        }
        if i > 0 && i - 1 < cut {
            b.cx(carries[i - 1], acc[i]);
        }
    }

    // Measurement-uncompute carries in reverse.  For addition the post-sum
    // identity is carry_{i+1} = majority(!acc_i_final, k_i, carry_i).
    for i in (0..cut).rev() {
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        let carry_in = if i == 0 { None } else { Some(carries[i - 1]) };
        if bit(c, i) {
            b.x(acc[i]);
            if let Some(ci) = carry_in {
                b.cz_if(acc[i], ctrl, m);
                b.cz_if(acc[i], ci, m);
                b.x(acc[i]);
                b.cz_if(ctrl, ci, m);
            } else {
                b.cz_if(acc[i], ctrl, m);
                b.x(acc[i]);
            }
        } else if let Some(ci) = carry_in {
            b.x(acc[i]);
            b.cz_if(acc[i], ci, m);
            b.x(acc[i]);
        }
    }

    b.free_vec(&carries);
}

pub(crate) fn add_nbit_const_fast(b: &mut B, acc: &[QubitId], c: U256) {
    let n = acc.len();
    let a = load_const(b, n, c);
    add_nbit_qq_fast(b, &a, acc);
    unload_const(b, &a, c);
}

pub(crate) fn sub_nbit_const_fast(b: &mut B, acc: &[QubitId], c: U256) {
    let n = acc.len();
    let a = load_const(b, n, c);
    sub_nbit_qq_fast(b, &a, acc);
    unload_const(b, &a, c);
}
