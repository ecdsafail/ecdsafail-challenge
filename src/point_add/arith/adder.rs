use super::*;

pub(crate) fn bit(c: U256, i: usize) -> bool {
    // alloy's U256::bit returns bool for index < 256.
    c.bit(i)
}

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

/// Same arithmetic as `cuccaro_add_fast`, but the carry lane is supplied by the
/// caller and must be clean on entry.  The HMR uncompute returns it to zero, so
/// Kaliski step4 can reuse clean high `tmp` lanes without increasing peak Q.
pub(crate) fn cuccaro_add_fast_borrowed_carries(
    b: &mut B,
    a: &[QubitId],
    acc: &[QubitId],
    c_in: QubitId,
    carries: &[QubitId],
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
    assert!(carries.len() >= n - 1);

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

/// Clean (X/CX/CCX only, emit_inverse-safe) Cuccaro add of an n-bit register
/// `a` into an (n+1)-bit accumulator `acc_ext`, capturing the carry-out into
/// `acc_ext[n]`. `acc_ext` may hold any (n+1)-bit value on entry; `c_in` is a
/// fresh ancilla at |0> that returns to |0>.
///
/// Unlike [`cuccaro_add`] (which discards the carry-out, omitting the top MAJ),
/// this runs the *full* n-step MAJ sweep so the carry-out is materialized in
/// `a[n-1]` after the sweep; we CX it into `acc_ext[n]`, then run the full UMA
/// sweep to write the sum bits and restore `a` and `c_in`. This is the
/// MAJ/UMA analogue of [`cuccaro_add_fast_low_to_ext`] (no measurement), so it
/// is safe inside `emit_inverse` blocks. `a` is preserved.
pub(crate) fn cuccaro_add_low_to_ext_clean(
    b: &mut B,
    a: &[QubitId],
    acc_ext: &[QubitId],
    c_in: QubitId,
) {
    let n = a.len();
    assert_eq!(acc_ext.len(), n + 1);
    if n == 0 {
        // acc_ext[0] += c_in.
        b.cx(c_in, acc_ext[0]);
        return;
    }

    // Full forward MAJ sweep (bits 0..=n-1). After this, a[n-1] holds the
    // carry-out of the whole addition.
    maj(b, c_in, acc_ext[0], a[0]);
    for i in 1..n {
        maj(b, a[i - 1], acc_ext[i], a[i]);
    }

    // Carry-out into the extension bit.
    b.cx(a[n - 1], acc_ext[n]);

    // Full reverse UMA sweep: writes sum bits into acc_ext[0..n], restores a
    // and c_in to their entry values.
    for i in (1..n).rev() {
        uma(b, a[i - 1], acc_ext[i], a[i]);
    }
    uma(b, c_in, acc_ext[0], a[0]);
}

/// Gate-level inverse of [`cuccaro_add_low_to_ext_clean`]: computes
/// `acc_ext := acc_ext - (a + c_in)` capturing the borrow-out into
/// `acc_ext[n]` (the same bit toggles, since add and subtract share the carry
/// identity under the running ext bit). `a` is preserved; `c_in` clean in/out.
pub(crate) fn cuccaro_sub_low_to_ext_clean(
    b: &mut B,
    a: &[QubitId],
    acc_ext: &[QubitId],
    c_in: QubitId,
) {
    let n = a.len();
    assert_eq!(acc_ext.len(), n + 1);
    if n == 0 {
        b.cx(c_in, acc_ext[0]);
        return;
    }

    // Inverse of the forward UMA sweep.
    inv_uma(b, c_in, acc_ext[0], a[0]);
    for i in 1..n {
        inv_uma(b, a[i - 1], acc_ext[i], a[i]);
    }

    // Inverse of the carry-out write (CX is self-inverse).
    b.cx(a[n - 1], acc_ext[n]);

    // Inverse of the forward MAJ sweep.
    for i in (1..n).rev() {
        inv_maj(b, a[i - 1], acc_ext[i], a[i]);
    }
    inv_maj(b, c_in, acc_ext[0], a[0]);
}


/// Hybrid Cuccaro add `acc += a mod 2^n` where the LOW `k` bits use the measured
/// fast adder (k carry lanes live => ~1 CCX/bit) and the HIGH `n-k` bits use the
/// coherent MAJ/UMA adder (0 carry ancilla => ~2 CCX/bit). The ripple carry out of
/// the low block is threaded into the high block via a single fresh `cout` qubit,
/// so peak ancilla = k (+1 transient cout) instead of n-1. `a` is preserved,
/// `c_in` clean in/out. Same arithmetic as [`cuccaro_add`]; the high carry-out is
/// discarded (mod 2^n). For k>=n this is just the all-fast adder; for k==0 the
/// all-coherent adder.
pub(crate) fn cuccaro_add_hybrid_lowfast(
    b: &mut B,
    a: &[QubitId],
    acc: &[QubitId],
    c_in: QubitId,
    k: usize,
) {
    let n = a.len();
    assert_eq!(n, acc.len());
    if n == 0 {
        return;
    }
    let k = k.min(n);
    if k == 0 {
        cuccaro_add(b, a, acc, c_in);
        return;
    }
    if k >= n {
        cuccaro_add_fast(b, a, acc, c_in);
        return;
    }
    // Low block: fast add of a[..k] into acc[..k], carry-out captured in `cout`.
    let cout = b.alloc_qubit();
    let mut acc_lo_ext: Vec<QubitId> = acc[..k].to_vec();
    acc_lo_ext.push(cout);
    cuccaro_add_fast_low_to_ext(b, &a[..k], &acc_lo_ext, c_in);
    // High block: coherent add of a[k..] into acc[k..] with carry-in = cout.
    // cuccaro_add consumes the carry-in (returns it to |0>) and discards its own
    // carry-out (mod 2^{n-k}) -- exactly the high half of a mod-2^n ripple add.
    cuccaro_add(b, &a[k..], &acc[k..], cout);
    b.free(cout);
}

/// Hybrid Cuccaro subtract `acc -= a mod 2^n`: LOW `k` bits measured-fast, HIGH
/// `n-k` bits coherent, borrow threaded through one fresh `cout`. Inverse-shaped
/// twin of [`cuccaro_add_hybrid_lowfast`]; same peak (k+1) and CCX profile.
pub(crate) fn cuccaro_sub_hybrid_lowfast(
    b: &mut B,
    a: &[QubitId],
    acc: &[QubitId],
    c_in: QubitId,
    k: usize,
) {
    let n = a.len();
    assert_eq!(n, acc.len());
    if n == 0 {
        return;
    }
    let k = k.min(n);
    if k == 0 {
        cuccaro_sub(b, a, acc, c_in);
        return;
    }
    if k >= n {
        cuccaro_sub_fast(b, a, acc, c_in);
        return;
    }
    let cout = b.alloc_qubit();
    let mut acc_lo_ext: Vec<QubitId> = acc[..k].to_vec();
    acc_lo_ext.push(cout);
    cuccaro_sub_fast_low_to_ext(b, &a[..k], &acc_lo_ext, c_in);
    cuccaro_sub(b, &a[k..], &acc[k..], cout);
    b.free(cout);
}

pub(crate) fn load_const(b: &mut B, n: usize, c: U256) -> Vec<QubitId> {
    let qs = b.alloc_qubits(n);
    for i in 0..n {
        if bit(c, i) {
            b.x(qs[i]);
        }
    }
    qs
}

pub(crate) fn unload_const(b: &mut B, qs: &[QubitId], c: U256) {
    for i in 0..qs.len() {
        if bit(c, i) {
            b.x(qs[i]);
        }
    }
    b.free_vec(qs);
}

pub(crate) fn load_bits(b: &mut B, bits: &[BitId]) -> Vec<QubitId> {
    let n = bits.len();
    let qs = b.alloc_qubits(n);
    for i in 0..n {
        // qs[i] ← bits[i] via conditional X
        b.x_if(qs[i], bits[i]);
    }
    qs
}

pub(crate) fn unload_bits(b: &mut B, qs: &[QubitId], bits: &[BitId]) {
    for i in 0..qs.len() {
        b.x_if(qs[i], bits[i]);
    }
    b.free_vec(qs);
}

/// Build an (n+1)-bit view by attaching a freshly-allocated 0 ancilla.
pub(crate) fn ext_reg(b: &mut B, reg: &[QubitId]) -> (Vec<QubitId>, QubitId) {
    let ovf = b.alloc_qubit();
    let mut r = reg.to_vec();
    r.push(ovf);
    (r, ovf)
}

/// Release the overflow ancilla (which must be 0 on exit).
pub(crate) fn unext_reg(b: &mut B, ovf: QubitId) {
    b.free(ovf);
}

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

/// Fast Cuccaro add into an extended accumulator where the source high bit is
/// known zero: `acc_ext += a + c_in (mod 2^(n+1))`.
pub(crate) fn cuccaro_add_fast_low_to_ext(b: &mut B, a: &[QubitId], acc_ext: &[QubitId], c_in: QubitId) {
    let n = a.len();
    assert_eq!(acc_ext.len(), n + 1);
    if n == 0 {
        b.cx(c_in, acc_ext[0]);
        return;
    }

    let carries = b.alloc_qubits(n);

    b.cx(a[0], acc_ext[0]);
    b.cx(a[0], c_in);
    b.ccx(c_in, acc_ext[0], carries[0]);
    b.cx(carries[0], a[0]);
    for i in 1..n {
        b.cx(a[i], acc_ext[i]);
        b.cx(a[i], a[i - 1]);
        b.ccx(a[i - 1], acc_ext[i], carries[i]);
        b.cx(carries[i], a[i]);
    }

    b.cx(a[n - 1], acc_ext[n]);

    for i in (1..n).rev() {
        b.cx(carries[i], a[i]);
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        b.cz_if(a[i - 1], acc_ext[i], m);
        b.cx(a[i], a[i - 1]);
        b.cx(a[i - 1], acc_ext[i]);
    }
    b.cx(carries[0], a[0]);
    let m0 = b.alloc_bit();
    b.hmr(carries[0], m0);
    b.cz_if(c_in, acc_ext[0], m0);
    b.cx(a[0], c_in);
    b.cx(c_in, acc_ext[0]);

    b.free_vec(&carries);
}

/// Fast Cuccaro subtract from an extended accumulator where the source high bit
/// is known zero: `acc_ext -= a + c_in (mod 2^(n+1))`.
pub(crate) fn cuccaro_sub_fast_low_to_ext(b: &mut B, a: &[QubitId], acc_ext: &[QubitId], c_in: QubitId) {
    let n = a.len();
    assert_eq!(acc_ext.len(), n + 1);
    if n == 0 {
        b.cx(c_in, acc_ext[0]);
        return;
    }

    let carries = b.alloc_qubits(n);

    b.cx(c_in, acc_ext[0]);
    b.cx(a[0], c_in);
    b.ccx(c_in, acc_ext[0], carries[0]);
    b.cx(carries[0], a[0]);
    for i in 1..n {
        b.cx(a[i - 1], acc_ext[i]);
        b.cx(a[i], a[i - 1]);
        b.ccx(a[i - 1], acc_ext[i], carries[i]);
        b.cx(carries[i], a[i]);
    }

    b.cx(a[n - 1], acc_ext[n]);

    for i in (1..n).rev() {
        b.cx(carries[i], a[i]);
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        b.cz_if(a[i - 1], acc_ext[i], m);
        b.cx(a[i], a[i - 1]);
        b.cx(a[i], acc_ext[i]);
    }
    b.cx(carries[0], a[0]);
    let m0 = b.alloc_bit();
    b.hmr(carries[0], m0);
    b.cz_if(c_in, acc_ext[0], m0);
    b.cx(a[0], c_in);
    b.cx(a[0], acc_ext[0]);

    b.free_vec(&carries);
}

/// Borrowed-carry form of [`cuccaro_add_fast_low_to_ext`].  The source has no
/// materialized high-zero pad lane: `acc_ext` is one bit wider than `a`, and
/// the caller supplies `a.len()` clean, pairwise-disjoint carry lanes.
pub(crate) fn cuccaro_add_fast_low_to_ext_borrowed_carries(
    b: &mut B,
    a: &[QubitId],
    acc_ext: &[QubitId],
    c_in: QubitId,
    carries: &[QubitId],
) {
    let n = a.len();
    assert_eq!(acc_ext.len(), n + 1);
    if n == 0 {
        b.cx(c_in, acc_ext[0]);
        return;
    }
    assert!(carries.len() >= n);

    b.cx(a[0], acc_ext[0]);
    b.cx(a[0], c_in);
    b.ccx(c_in, acc_ext[0], carries[0]);
    b.cx(carries[0], a[0]);
    for i in 1..n {
        b.cx(a[i], acc_ext[i]);
        b.cx(a[i], a[i - 1]);
        b.ccx(a[i - 1], acc_ext[i], carries[i]);
        b.cx(carries[i], a[i]);
    }

    b.cx(a[n - 1], acc_ext[n]);

    for i in (1..n).rev() {
        b.cx(carries[i], a[i]);
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        b.cz_if(a[i - 1], acc_ext[i], m);
        b.cx(a[i], a[i - 1]);
        b.cx(a[i - 1], acc_ext[i]);
    }
    b.cx(carries[0], a[0]);
    let m0 = b.alloc_bit();
    b.hmr(carries[0], m0);
    b.cz_if(c_in, acc_ext[0], m0);
    b.cx(a[0], c_in);
    b.cx(c_in, acc_ext[0]);
}

/// Borrowed-carry inverse of
/// [`cuccaro_add_fast_low_to_ext_borrowed_carries`].
pub(crate) fn cuccaro_sub_fast_low_to_ext_borrowed_carries(
    b: &mut B,
    a: &[QubitId],
    acc_ext: &[QubitId],
    c_in: QubitId,
    carries: &[QubitId],
) {
    let n = a.len();
    assert_eq!(acc_ext.len(), n + 1);
    if n == 0 {
        b.cx(c_in, acc_ext[0]);
        return;
    }
    assert!(carries.len() >= n);

    b.cx(c_in, acc_ext[0]);
    b.cx(a[0], c_in);
    b.ccx(c_in, acc_ext[0], carries[0]);
    b.cx(carries[0], a[0]);
    for i in 1..n {
        b.cx(a[i - 1], acc_ext[i]);
        b.cx(a[i], a[i - 1]);
        b.ccx(a[i - 1], acc_ext[i], carries[i]);
        b.cx(carries[i], a[i]);
    }

    b.cx(a[n - 1], acc_ext[n]);

    for i in (1..n).rev() {
        b.cx(carries[i], a[i]);
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        b.cz_if(a[i - 1], acc_ext[i], m);
        b.cx(a[i], a[i - 1]);
        b.cx(a[i], acc_ext[i]);
    }
    b.cx(carries[0], a[0]);
    let m0 = b.alloc_bit();
    b.hmr(carries[0], m0);
    b.cz_if(c_in, acc_ext[0], m0);
    b.cx(a[0], c_in);
    b.cx(a[0], acc_ext[0]);
}

/// Zero-carry-in specialization of
/// [`cuccaro_add_fast_low_to_ext_borrowed_carries`].  The omitted `c_in`
/// register is known zero: its only forward role is to preserve the original
/// low source bit until the measured carry clear.  After that clear `a[0]`
/// holds the same value, so it can control the phase correction directly.
pub(crate) fn cuccaro_add_fast_low_to_ext_borrowed_carries_no_cin(
    b: &mut B,
    a: &[QubitId],
    acc_ext: &[QubitId],
    carries: &[QubitId],
) {
    let n = a.len();
    assert_eq!(acc_ext.len(), n + 1);
    if n == 0 {
        return;
    }
    let gate_suffix = square_selfhost_gate_suffix_carries(n);
    let borrowed = n - gate_suffix;
    assert!(carries.len() >= borrowed);

    b.cx(a[0], acc_ext[0]);
    b.ccx(a[0], acc_ext[0], carries[0]);
    b.cx(carries[0], a[0]);
    for i in 1..borrowed {
        b.cx(a[i], acc_ext[i]);
        b.cx(a[i], a[i - 1]);
        b.ccx(a[i - 1], acc_ext[i], carries[i]);
        b.cx(carries[i], a[i]);
    }
    for i in borrowed..n {
        maj(b, a[i - 1], acc_ext[i], a[i]);
    }

    b.cx(a[n - 1], acc_ext[n]);

    for i in (borrowed..n).rev() {
        uma(b, a[i - 1], acc_ext[i], a[i]);
    }
    for i in (1..borrowed).rev() {
        b.cx(carries[i], a[i]);
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        b.cz_if(a[i - 1], acc_ext[i], m);
        b.cx(a[i], a[i - 1]);
        b.cx(a[i - 1], acc_ext[i]);
    }
    b.cx(carries[0], a[0]);
    let m0 = b.alloc_bit();
    b.hmr(carries[0], m0);
    b.cz_if(a[0], acc_ext[0], m0);
}

/// Zero-carry-in inverse of
/// [`cuccaro_add_fast_low_to_ext_borrowed_carries_no_cin`].
pub(crate) fn cuccaro_sub_fast_low_to_ext_borrowed_carries_no_cin(
    b: &mut B,
    a: &[QubitId],
    acc_ext: &[QubitId],
    carries: &[QubitId],
) {
    let n = a.len();
    assert_eq!(acc_ext.len(), n + 1);
    if n == 0 {
        return;
    }
    let gate_suffix = square_selfhost_gate_suffix_carries(n);
    let borrowed = n - gate_suffix;
    assert!(carries.len() >= borrowed);

    b.ccx(a[0], acc_ext[0], carries[0]);
    b.cx(carries[0], a[0]);
    for i in 1..borrowed {
        b.cx(a[i - 1], acc_ext[i]);
        b.cx(a[i], a[i - 1]);
        b.ccx(a[i - 1], acc_ext[i], carries[i]);
        b.cx(carries[i], a[i]);
    }
    for i in borrowed..n {
        inv_uma(b, a[i - 1], acc_ext[i], a[i]);
    }

    b.cx(a[n - 1], acc_ext[n]);

    for i in (borrowed..n).rev() {
        inv_maj(b, a[i - 1], acc_ext[i], a[i]);
    }
    for i in (1..borrowed).rev() {
        b.cx(carries[i], a[i]);
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        b.cz_if(a[i - 1], acc_ext[i], m);
        b.cx(a[i], a[i - 1]);
        b.cx(a[i], acc_ext[i]);
    }
    b.cx(carries[0], a[0]);
    let m0 = b.alloc_bit();
    b.hmr(carries[0], m0);
    b.cz_if(a[0], acc_ext[0], m0);
    b.cx(a[0], acc_ext[0]);
}

/// Add a materialized low prefix and an unmaterialized controlled high suffix.
/// The prefix's final carry is a valid controlled carry-in for the suffix.
pub(crate) fn cuccaro_add_fast_prefix_ctrl_suffix_no_cin(
    b: &mut B,
    prefix: &[QubitId],
    suffix: &[QubitId],
    acc: &[QubitId],
    ctrl: QubitId,
    carries: &[QubitId],
    scratch: QubitId,
) {
    let n = prefix.len();
    assert!(n > 0);
    assert!(!suffix.is_empty());
    assert_eq!(acc.len(), n + suffix.len());
    assert!(carries.len() >= n);

    b.cx(prefix[0], acc[0]);
    b.ccx(prefix[0], acc[0], carries[0]);
    b.cx(carries[0], prefix[0]);
    for i in 1..n {
        b.cx(prefix[i], acc[i]);
        b.cx(prefix[i], prefix[i - 1]);
        b.ccx(prefix[i - 1], acc[i], carries[i]);
        b.cx(carries[i], prefix[i]);
    }

    cuccaro_add_ctrl_lowq(b, suffix, &acc[n..], ctrl, prefix[n - 1], scratch);

    for i in (1..n).rev() {
        b.cx(carries[i], prefix[i]);
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        b.cz_if(prefix[i - 1], acc[i], m);
        b.cx(prefix[i], prefix[i - 1]);
        b.cx(prefix[i - 1], acc[i]);
    }
    b.cx(carries[0], prefix[0]);
    let m0 = b.alloc_bit();
    b.hmr(carries[0], m0);
    b.cz_if(prefix[0], acc[0], m0);
}

/// Inverse of [`cuccaro_add_fast_prefix_ctrl_suffix_no_cin`].
pub(crate) fn cuccaro_sub_fast_prefix_ctrl_suffix_no_cin(
    b: &mut B,
    prefix: &[QubitId],
    suffix: &[QubitId],
    acc: &[QubitId],
    ctrl: QubitId,
    carries: &[QubitId],
    scratch: QubitId,
) {
    let n = prefix.len();
    assert!(n > 0);
    assert!(!suffix.is_empty());
    assert_eq!(acc.len(), n + suffix.len());
    assert!(carries.len() >= n);

    b.ccx(prefix[0], acc[0], carries[0]);
    b.cx(carries[0], prefix[0]);
    for i in 1..n {
        b.cx(prefix[i - 1], acc[i]);
        b.cx(prefix[i], prefix[i - 1]);
        b.ccx(prefix[i - 1], acc[i], carries[i]);
        b.cx(carries[i], prefix[i]);
    }

    cuccaro_sub_ctrl_lowq(b, suffix, &acc[n..], ctrl, prefix[n - 1], scratch);

    for i in (1..n).rev() {
        b.cx(carries[i], prefix[i]);
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        b.cz_if(prefix[i - 1], acc[i], m);
        b.cx(prefix[i], prefix[i - 1]);
        b.cx(prefix[i], acc[i]);
    }
    b.cx(carries[0], prefix[0]);
    let m0 = b.alloc_bit();
    b.hmr(carries[0], m0);
    b.cz_if(prefix[0], acc[0], m0);
    b.cx(prefix[0], acc[0]);
}

/// Top-clean hybrid of [`cuccaro_add_fast_prefix_ctrl_suffix_no_cin`].  The top
/// `clean_top` bits of the materialized prefix ripple use the zero-ancilla
/// MAJ/UMA path, so the fast carry-lane region only spans `borrowed = n -
/// clean_top` low bits and needs `borrowed` carry lanes instead of `n`.  The
/// seam carry handed to the controlled suffix is still `prefix[n-1]`, exactly as
/// the baseline (the clean MAJ chain leaves the running carry in prefix[n-1]).
/// Value/phase identical to the baseline; only the low<->clean split moves.
pub(crate) fn cuccaro_add_fast_prefix_ctrl_suffix_no_cin_topclean(
    b: &mut B,
    prefix: &[QubitId],
    suffix: &[QubitId],
    acc: &[QubitId],
    ctrl: QubitId,
    carries: &[QubitId],
    scratch: QubitId,
    clean_top: usize,
) {
    let n = prefix.len();
    assert!(n > 0);
    assert!(!suffix.is_empty());
    assert_eq!(acc.len(), n + suffix.len());
    let clean_top = clean_top.min(n.saturating_sub(1));
    if clean_top == 0 {
        return cuccaro_add_fast_prefix_ctrl_suffix_no_cin(
            b, prefix, suffix, acc, ctrl, carries, scratch,
        );
    }
    let borrowed = n - clean_top;
    assert!(carries.len() >= borrowed);

    // Fast carry-lane MAJ for [0, borrowed).
    b.cx(prefix[0], acc[0]);
    b.ccx(prefix[0], acc[0], carries[0]);
    b.cx(carries[0], prefix[0]);
    for i in 1..borrowed {
        b.cx(prefix[i], acc[i]);
        b.cx(prefix[i], prefix[i - 1]);
        b.ccx(prefix[i - 1], acc[i], carries[i]);
        b.cx(carries[i], prefix[i]);
    }
    // Clean MAJ for [borrowed, n): carry rides in prefix[i].
    for i in borrowed..n {
        maj(b, prefix[i - 1], acc[i], prefix[i]);
    }

    cuccaro_add_ctrl_lowq(b, suffix, &acc[n..], ctrl, prefix[n - 1], scratch);

    // Clean UMA for [borrowed, n) in reverse.
    for i in (borrowed..n).rev() {
        uma(b, prefix[i - 1], acc[i], prefix[i]);
    }
    // Fast measured UMA for [0, borrowed) in reverse.
    for i in (1..borrowed).rev() {
        b.cx(carries[i], prefix[i]);
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        b.cz_if(prefix[i - 1], acc[i], m);
        b.cx(prefix[i], prefix[i - 1]);
        b.cx(prefix[i - 1], acc[i]);
    }
    b.cx(carries[0], prefix[0]);
    let m0 = b.alloc_bit();
    b.hmr(carries[0], m0);
    b.cz_if(prefix[0], acc[0], m0);
}

/// Inverse of [`cuccaro_add_fast_prefix_ctrl_suffix_no_cin_topclean`].
pub(crate) fn cuccaro_sub_fast_prefix_ctrl_suffix_no_cin_topclean(
    b: &mut B,
    prefix: &[QubitId],
    suffix: &[QubitId],
    acc: &[QubitId],
    ctrl: QubitId,
    carries: &[QubitId],
    scratch: QubitId,
    clean_top: usize,
) {
    let n = prefix.len();
    assert!(n > 0);
    assert!(!suffix.is_empty());
    assert_eq!(acc.len(), n + suffix.len());
    let clean_top = clean_top.min(n.saturating_sub(1));
    if clean_top == 0 {
        return cuccaro_sub_fast_prefix_ctrl_suffix_no_cin(
            b, prefix, suffix, acc, ctrl, carries, scratch,
        );
    }
    let borrowed = n - clean_top;
    assert!(carries.len() >= borrowed);

    // Fast carry-lane forward sweep for [0, borrowed).
    b.ccx(prefix[0], acc[0], carries[0]);
    b.cx(carries[0], prefix[0]);
    for i in 1..borrowed {
        b.cx(prefix[i - 1], acc[i]);
        b.cx(prefix[i], prefix[i - 1]);
        b.ccx(prefix[i - 1], acc[i], carries[i]);
        b.cx(carries[i], prefix[i]);
    }
    // Clean inv_uma for [borrowed, n) (forward, mirrors the body sub-topclean).
    for i in borrowed..n {
        inv_uma(b, prefix[i - 1], acc[i], prefix[i]);
    }

    cuccaro_sub_ctrl_lowq(b, suffix, &acc[n..], ctrl, prefix[n - 1], scratch);

    // Clean inv_maj for [borrowed, n) in reverse.
    for i in (borrowed..n).rev() {
        inv_maj(b, prefix[i - 1], acc[i], prefix[i]);
    }
    // Fast measured reverse sweep for [0, borrowed).
    for i in (1..borrowed).rev() {
        b.cx(carries[i], prefix[i]);
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        b.cz_if(prefix[i - 1], acc[i], m);
        b.cx(prefix[i], prefix[i - 1]);
        b.cx(prefix[i], acc[i]);
    }
    b.cx(carries[0], prefix[0]);
    let m0 = b.alloc_bit();
    b.hmr(carries[0], m0);
    b.cz_if(prefix[0], acc[0], m0);
    b.cx(prefix[0], acc[0]);
}


pub(crate) fn cuccaro_add_fast_windowed_low_to_ext(
    b: &mut B,
    a: &[QubitId],
    acc_ext: &[QubitId],
    c_in: QubitId,
    blocks: usize,
) {
    let n = a.len();
    assert_eq!(acc_ext.len(), n + 1);
    let ext_n = acc_ext.len();
    if ext_n == 0 {
        return;
    }
    let blocks = blocks.max(1).min(ext_n);
    if blocks == 1 {
        cuccaro_add_fast_low_to_ext(b, a, acc_ext, c_in);
        return;
    }

    let mut carry = c_in;
    let mut lo = 0usize;
    let mut couts: Vec<(QubitId, usize, QubitId)> = Vec::new();
    for blk in 0..blocks {
        let hi = ((blk + 1) * ext_n) / blocks;
        if hi <= lo {
            continue;
        }
        if blk == blocks - 1 || hi == ext_n {
            cuccaro_add_fast_low_to_ext(b, &a[lo..n], &acc_ext[lo..hi], carry);
            break;
        }
        let cout = b.alloc_qubit();
        let zero = b.alloc_qubit();
        let mut a_block: Vec<QubitId> = a[lo..hi].to_vec();
        a_block.push(zero);
        let mut acc_block: Vec<QubitId> = acc_ext[lo..hi].to_vec();
        acc_block.push(cout);
        let c_in = carry;
        cuccaro_add_fast(b, &a_block, &acc_block, carry);
        b.free(zero);
        couts.push((cout, hi, c_in));
        carry = cout;
        lo = hi;
    }

    for &(cout, p, c_in) in couts.iter().rev() {
        cmp_lt_into_fast_with_cin(b, &acc_ext[..p], &a[..p], c_in, cout);
        b.free(cout);
    }
}

/// Asymmetric 2-block fast windowed add `acc_ext(n+1) += a(n)` (carry into
/// acc_ext[n]). Both blocks measured-fast; the single inter-block carry is
/// uncomputed to |0> via one cmp_lt over the LOW `s` bits. Unlike the equal-split
/// [`cuccaro_add_fast_windowed_low_to_ext`], the boundary is placed at a chosen
/// `s`, so the cmp_lt uncompute cost (~s CCX) is decoupled from the (large) high
/// block: pick `s` small and the high block as wide as peak allows to capture most
/// of the fast saving for ~s overhead. Peak ancilla = max(s+1, n+1-s) carry lanes.
/// `a` preserved; `c_in` clean in/out. Phase-clean (carry returns to |0>).
pub(crate) fn cuccaro_add_fast_split2_low_to_ext(
    b: &mut B,
    a: &[QubitId],
    acc_ext: &[QubitId],
    c_in: QubitId,
    s: usize,
) {
    let n = a.len();
    assert_eq!(acc_ext.len(), n + 1);
    let ext_n = n + 1;
    // Degenerate: no useful boundary -> plain fast low_to_ext.
    if s == 0 || s >= ext_n - 1 {
        cuccaro_add_fast_low_to_ext(b, a, acc_ext, c_in);
        return;
    }
    // Low block [0,s): fast add capturing carry-out into `cout`.
    let cout = b.alloc_qubit();
    let zero = b.alloc_qubit();
    let mut a_block: Vec<QubitId> = a[0..s].to_vec();
    a_block.push(zero);
    let mut acc_block: Vec<QubitId> = acc_ext[0..s].to_vec();
    acc_block.push(cout);
    cuccaro_add_fast(b, &a_block, &acc_block, c_in);
    b.free(zero);
    // High block [s, n+1): fast add of a[s..n] into acc_ext[s..n+1], carry-in cout.
    cuccaro_add_fast_low_to_ext(b, &a[s..n], &acc_ext[s..ext_n], cout);
    // Uncompute cout to |0>: recompute the low-block carry (cmp_lt over s bits).
    cmp_lt_into_fast_with_cin(b, &acc_ext[..s], &a[..s], c_in, cout);
    b.free(cout);
}

/// Subtract twin of [`cuccaro_add_fast_split2_low_to_ext`]: `acc_ext -= a` with
/// borrow into acc_ext[n], borrow threaded across one boundary at `s` and
/// uncomputed via the borrow-form cmp_lt (X-conjugated). Phase-clean.
pub(crate) fn cuccaro_sub_fast_split2_low_to_ext(
    b: &mut B,
    a: &[QubitId],
    acc_ext: &[QubitId],
    c_in: QubitId,
    s: usize,
) {
    let n = a.len();
    assert_eq!(acc_ext.len(), n + 1);
    let ext_n = n + 1;
    if s == 0 || s >= ext_n - 1 {
        cuccaro_sub_fast_low_to_ext(b, a, acc_ext, c_in);
        return;
    }
    let bout = b.alloc_qubit();
    let zero = b.alloc_qubit();
    let mut a_block: Vec<QubitId> = a[0..s].to_vec();
    a_block.push(zero);
    let mut acc_block: Vec<QubitId> = acc_ext[0..s].to_vec();
    acc_block.push(bout);
    cuccaro_sub_fast(b, &a_block, &acc_block, c_in);
    b.free(zero);
    cuccaro_sub_fast_low_to_ext(b, &a[s..n], &acc_ext[s..ext_n], bout);
    // Borrow-form uncompute (mirrors the equal-split primitive's sub branch).
    for i in 0..s {
        b.x(a[i]);
    }
    cmp_lt_into_fast_with_cin(b, &a[..s], &acc_ext[..s], c_in, bout);
    for i in 0..s {
        b.x(a[i]);
    }
    b.free(bout);
}

pub(crate) fn cuccaro_sub_fast_windowed_low_to_ext(
    b: &mut B,
    a: &[QubitId],
    acc_ext: &[QubitId],
    c_in: QubitId,
    blocks: usize,
) {
    let n = a.len();
    assert_eq!(acc_ext.len(), n + 1);
    let ext_n = acc_ext.len();
    if ext_n == 0 {
        return;
    }
    let blocks = blocks.max(1).min(ext_n);
    if blocks == 1 {
        cuccaro_sub_fast_low_to_ext(b, a, acc_ext, c_in);
        return;
    }

    let mut borrow = c_in;
    let mut lo = 0usize;
    let mut bouts: Vec<(QubitId, usize, QubitId)> = Vec::new();
    for blk in 0..blocks {
        let hi = ((blk + 1) * ext_n) / blocks;
        if hi <= lo {
            continue;
        }
        if blk == blocks - 1 || hi == ext_n {
            cuccaro_sub_fast_low_to_ext(b, &a[lo..n], &acc_ext[lo..hi], borrow);
            break;
        }
        let bout = b.alloc_qubit();
        let zero = b.alloc_qubit();
        let mut a_block: Vec<QubitId> = a[lo..hi].to_vec();
        a_block.push(zero);
        let mut acc_block: Vec<QubitId> = acc_ext[lo..hi].to_vec();
        acc_block.push(bout);
        let b_in = borrow;
        cuccaro_sub_fast(b, &a_block, &acc_block, borrow);
        b.free(zero);
        bouts.push((bout, hi, b_in));
        borrow = bout;
        lo = hi;
    }

    for &(bout, p, b_in) in bouts.iter().rev() {
        for i in 0..p {
            b.x(a[i]);
        }
        cmp_lt_into_fast_with_cin(b, &a[..p], &acc_ext[..p], b_in, bout);
        for i in 0..p {
            b.x(a[i]);
        }
        b.free(bout);
    }
}


pub(crate) fn cuccaro_sub_fast_borrowed_carries(
    b: &mut B,
    a: &[QubitId],
    acc: &[QubitId],
    c_in: QubitId,
    carries: &[QubitId],
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
    assert!(carries.len() >= n - 1);

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
}

/// Zero-carry-in specialization of [`cuccaro_add_fast_borrowed_carries`]
/// (same-width, `acc += a mod 2^n`, no carry-out captured). The omitted `c_in`
/// register is *proven* |0> on entry: its only forward roles are (a) to seed the
/// MAJ chain at bit 0 with carry-in 0 and (b) to freeze the original `a[0]` until
/// the final measured UMA's phase correction. With c_in=0 the seed
/// `cx(c_in,acc[0]); cx(a[0],c_in); ccx(c_in,acc[0],c0)` collapses to
/// `ccx(a[0],acc[0],c0)`, and since c_in held `a[0]` (restored by the final
/// `cx(carries[0],a[0])` to its seed-time value) the final `cz_if(c_in,acc[0],m0)`
/// equals `cz_if(a[0],acc[0],m0)`. This is the same-width analogue of the proven
/// [`cuccaro_add_fast_low_to_ext_borrowed_carries_no_cin`]. Consumes NO `c_in`
/// qubit; `carries` must be clean on entry and is restored to |0>.
pub(crate) fn cuccaro_add_fast_borrowed_carries_no_cin(
    b: &mut B,
    a: &[QubitId],
    acc: &[QubitId],
    carries: &[QubitId],
) {
    let n = a.len();
    assert_eq!(n, acc.len());
    if n == 0 {
        return;
    }
    if n == 1 {
        // acc[0] += a[0] (c_in = 0); pure XOR, no carry lane needed.
        b.cx(a[0], acc[0]);
        return;
    }
    assert!(carries.len() >= n - 1);

    // Step 0 MAJ with c_in folded out (c_in == 0 == a[0]'s seed companion).
    b.cx(a[0], acc[0]);
    b.ccx(a[0], acc[0], carries[0]);
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
    // Step 0 UMA with c_in folded out. In the c_in form the tail is
    //   cz_if(c_in,acc[0],m0); cx(a[0],c_in); cx(c_in,acc[0])
    // where the pre-`cz_if` `cx(carries[0],a[0])` has restored a[0] to the
    // frozen c_in value, so `cz_if(c_in,..)` == `cz_if(a[0],..)`. The two
    // trailing CXs reset c_in (`cx(a[0],c_in)`) and then `cx(c_in,acc[0])`
    // with c_in already 0 — a no-op. Both drop out: NO trailing acc CX here.
    b.cx(carries[0], a[0]);
    let m0 = b.alloc_bit();
    b.hmr(carries[0], m0);
    b.cz_if(a[0], acc[0], m0);
}

/// Zero-carry-in inverse of [`cuccaro_add_fast_borrowed_carries_no_cin`]:
/// same-width `acc -= a mod 2^n`, derived from
/// [`cuccaro_sub_fast_borrowed_carries`] by folding out the proven-|0> `c_in`
/// exactly as in the add direction. Consumes NO `c_in` qubit; `carries` clean in
/// and restored to |0>.
pub(crate) fn cuccaro_sub_fast_borrowed_carries_no_cin(
    b: &mut B,
    a: &[QubitId],
    acc: &[QubitId],
    carries: &[QubitId],
) {
    let n = a.len();
    assert_eq!(n, acc.len());
    if n == 0 {
        return;
    }
    if n == 1 {
        // acc[0] -= a[0] (c_in = 0); pure XOR.
        b.cx(a[0], acc[0]);
        return;
    }
    assert!(carries.len() >= n - 1);

    // Step 0 with c_in folded out (the sub seed begins ccx(a[0],acc[0],c0)).
    b.ccx(a[0], acc[0], carries[0]);
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
    b.cz_if(a[0], acc[0], m0);
    b.cx(a[0], acc[0]);
}

/// Top-clean hybrid of [`cuccaro_add_fast_borrowed_carries_no_cin`].  Identical
/// arithmetic (`acc += a mod 2^n`, `a` preserved, `carries` restored to |0>),
/// but the top `clean_top` bit positions run the zero-ancilla MAJ/UMA
/// (CCX-uncomputed) while the low `borrowed = n - clean_top` positions keep
/// the fast measured-uncompute carry lanes.  This borrows only
/// `borrowed` carry lanes instead of `n - 1`, lowering transient peak by
/// `clean_top` qubits.  `clean_top` is clamped to `[0, n]`.  The caller
/// supplies `carries[..borrowed]` (same convention as the `_no_cin` baseline).
pub(crate) fn cuccaro_add_fast_borrowed_carries_no_cin_topclean(
    b: &mut B,
    a: &[QubitId],
    acc: &[QubitId],
    carries: &[QubitId],
    clean_top: usize,
) {
    let n = a.len();
    assert_eq!(n, acc.len());
    let clean_top = clean_top.min(n);
    let borrowed = n - clean_top;

    if borrowed == 0 {
        // Pure clean MAJ/UMA ripple (same as _no_cin, which already folds c_in=0).
        return cuccaro_add_fast_borrowed_carries_no_cin(b, a, acc, carries);
    }
    if clean_top == 0 {
        return cuccaro_add_fast_borrowed_carries_no_cin(b, a, acc, carries);
    }

    // Now 0 < borrowed < n and 0 < clean_top < n.
    // The fast region covers bits [0, borrowed) with borrowed MAJ gates
    // (step 0 + loop 1..borrowed), needing `borrowed` carry lanes.
    assert!(carries.len() >= borrowed);

    // ── Forward sweep ──────────────────────────────────────────────────────
    // Fast carry-lane MAJ for indices [0, borrowed).  After this sweep the
    // running carry lives in a[borrowed-1] (the seam).
    // Step 0 (c_in folded): cx(a0,acc0); ccx(a0,acc0,c0); cx(c0,a0)
    b.cx(a[0], acc[0]);
    b.ccx(a[0], acc[0], carries[0]);
    b.cx(carries[0], a[0]);
    // Steps 1..borrowed-1:
    for i in 1..borrowed {
        b.cx(a[i], acc[i]);
        b.cx(a[i], a[i - 1]);
        b.ccx(a[i - 1], acc[i], carries[i]);
        b.cx(carries[i], a[i]);
    }
    // Clean (zero-ancilla) MAJ for indices [borrowed, n).
    // At i=borrowed the carry-in x is a[borrowed-1], which holds the seam carry.
    for i in borrowed..n - 1 {
        maj(b, a[i - 1], acc[i], a[i]);
    }

    // Final sum bit (same as the _no_cin baseline):
    //   acc[n-1] ^= a[n-2]; acc[n-1] ^= a[n-1]
    b.cx(a[n - 2], acc[n - 1]);
    b.cx(a[n - 1], acc[n - 1]);

    // ── Reverse sweep ──────────────────────────────────────────────────────
    // Clean UMA for [borrowed, n): reverse of the clean MAJ above.
    for i in (borrowed..n - 1).rev() {
        uma(b, a[i - 1], acc[i], a[i]);
    }
    // Fast measured UMA for [0, borrowed): inverse of the fast MAJ above.
    for i in (1..borrowed).rev() {
        b.cx(carries[i], a[i]);
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        b.cz_if(a[i - 1], acc[i], m);
        b.cx(a[i], a[i - 1]);
        b.cx(a[i - 1], acc[i]);
    }
    // Step 0 UMA (c_in folded):
    b.cx(carries[0], a[0]);
    let m0 = b.alloc_bit();
    b.hmr(carries[0], m0);
    b.cz_if(a[0], acc[0], m0);
}

/// Top-clean hybrid of [`cuccaro_sub_fast_borrowed_carries_no_cin`].  Same
/// carry-lane savings as the add variant; this is its exact subtract analogue
/// (clean `inv_uma`/`inv_maj` on the top `clean_top` positions, fast measured
/// carry uncompute on the low region).
pub(crate) fn cuccaro_sub_fast_borrowed_carries_no_cin_topclean(
    b: &mut B,
    a: &[QubitId],
    acc: &[QubitId],
    carries: &[QubitId],
    clean_top: usize,
) {
    let n = a.len();
    assert_eq!(n, acc.len());
    let clean_top = clean_top.min(n);
    let borrowed = n - clean_top;

    if borrowed == 0 {
        return cuccaro_sub_fast_borrowed_carries_no_cin(b, a, acc, carries);
    }
    if clean_top == 0 {
        return cuccaro_sub_fast_borrowed_carries_no_cin(b, a, acc, carries);
    }

    // ── Forward sweep ──────────────────────────────────────────────────────
    // Fast carry-lane sub seed for index 0 (c_in folded):
    b.ccx(a[0], acc[0], carries[0]);
    b.cx(carries[0], a[0]);
    // Steps 1..borrowed-1:
    for i in 1..borrowed {
        b.cx(a[i - 1], acc[i]);
        b.cx(a[i], a[i - 1]);
        b.ccx(a[i - 1], acc[i], carries[i]);
        b.cx(carries[i], a[i]);
    }
    // Clean inv_uma for [borrowed, n):
    for i in borrowed..n - 1 {
        inv_uma(b, a[i - 1], acc[i], a[i]);
    }

    // Final sum bit (same as the _no_cin baseline):
    b.cx(a[n - 1], acc[n - 1]);
    b.cx(a[n - 2], acc[n - 1]);

    // ── Reverse sweep ──────────────────────────────────────────────────────
    // Clean inv_maj for [borrowed, n):
    for i in (borrowed..n - 1).rev() {
        inv_maj(b, a[i - 1], acc[i], a[i]);
    }
    // Fast measured inv_MAJ for [0, borrowed):
    for i in (1..borrowed).rev() {
        b.cx(carries[i], a[i]);
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        b.cz_if(a[i - 1], acc[i], m);
        b.cx(a[i], a[i - 1]);
        b.cx(a[i], acc[i]);
    }
    // Step 0 (c_in folded):
    b.cx(carries[0], a[0]);
    let m0 = b.alloc_bit();
    b.hmr(carries[0], m0);
    b.cz_if(a[0], acc[0], m0);
    b.cx(a[0], acc[0]);
}

/// Build-time differential self-test for the `_topclean` no-c_in ripple
/// variants.  For a spread of widths `n` and every `clean_top` in
/// `0..=min(n,8)`, checks (a) `acc == (acc_in (+|-) a_in) mod 2^n`, (b) `a`
/// restored, (c) every carry lane is |0>, (d) all other ancilla qubits are
/// |0>, (e) global phase is 0.  Returns `Err` with a precise diagnostic on
/// the first failure.  Invoked from `point_add::build()` under
/// `BODY_TOPCLEAN_SELFTEST=1`.
pub(crate) fn body_topclean_selftest() -> Result<(), String> {
    use crate::circuit::QubitOrBit;

    // Deterministic 64-bit LCG.
    fn lcg(s: &mut u64) -> u64 {
        *s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        *s
    }
    fn rand_u256(s: &mut u64, n: usize) -> U256 {
        let mut limbs = [0u64; 4];
        for limb in limbs.iter_mut() {
            *limb = lcg(s);
        }
        let mut v = U256::from_limbs(limbs);
        if n < 256 {
            let mask = (U256::from(1u64) << n) - U256::from(1u64);
            v &= mask;
        }
        v
    }

    let widths = [2usize, 3, 5, 8, 16, 33, 64, 65, 128, 255];
    for &n in widths.iter() {
        let max_clean_top = n.min(8);
        for clean_top in 0..=max_clean_top {
            for is_sub in [false, true] {
                let mut b = B::new();
                let a = b.alloc_qubits(n);
                let acc = b.alloc_qubits(n);
                let carries = b.alloc_qubits(n.saturating_sub(1));
                if is_sub {
                    cuccaro_sub_fast_borrowed_carries_no_cin_topclean(
                        &mut b, &a, &acc, &carries, clean_top,
                    );
                } else {
                    cuccaro_add_fast_borrowed_carries_no_cin_topclean(
                        &mut b, &a, &acc, &carries, clean_top,
                    );
                }
                let nq = b.next_qubit as usize;
                let nb = b.next_bit as usize;

                let a_qb: Vec<QubitOrBit> = a.iter().map(|&q| QubitOrBit::Qubit(q)).collect();
                let acc_qb: Vec<QubitOrBit> =
                    acc.iter().map(|&q| QubitOrBit::Qubit(q)).collect();

                let mut seed = Shake256::default();
                seed.update(b"body-topclean-selftest");
                seed.update(&(n as u64).to_le_bytes());
                seed.update(&(clean_top as u64).to_le_bytes());
                seed.update(&[is_sub as u8]);
                let mut xof = seed.finalize_xof();
                let mut sim = Simulator::new(nq.max(1), nb.max(1), &mut xof);

                let mut rs = 0x9e3779b97f4a7c15u64 ^ ((n as u64) << 32) ^ (clean_top as u64);
                let modulus = U256::from(1u64) << n;
                let mut a_in = [U256::ZERO; 64];
                let mut acc_in = [U256::ZERO; 64];
                for shot in 0..64usize {
                    let av = rand_u256(&mut rs, n);
                    let accv = rand_u256(&mut rs, n);
                    a_in[shot] = av;
                    acc_in[shot] = accv;
                    sim.set_register(&a_qb, av, shot);
                    sim.set_register(&acc_qb, accv, shot);
                }

                sim.apply_iter(b.ops.iter());

                if sim.phase != 0 {
                    return Err(format!(
                        "phase leak: n={n} clean_top={clean_top} sub={is_sub} phase={:#x}",
                        sim.phase
                    ));
                }
                for shot in 0..64usize {
                    let sum = if is_sub {
                        (acc_in[shot] + modulus - a_in[shot]) % modulus
                    } else {
                        (acc_in[shot] + a_in[shot]) % modulus
                    };
                    let got = sim.get_register(&acc_qb, shot);
                    if got != sum {
                        return Err(format!(
                            "value mismatch: n={n} clean_top={clean_top} sub={is_sub} shot={shot} got={got:#x} exp={sum:#x} a={:#x} acc={:#x}",
                            a_in[shot], acc_in[shot]
                        ));
                    }
                    let a_got = sim.get_register(&a_qb, shot);
                    if a_got != a_in[shot] {
                        return Err(format!(
                            "source not restored: n={n} clean_top={clean_top} sub={is_sub} shot={shot} a_got={a_got:#x} a_in={:#x}",
                            a_in[shot]
                        ));
                    }
                    // Every ancilla qubit (outside a/acc) must be |0>.
                    let used: std::collections::HashSet<u64> = a
                        .iter()
                        .chain(acc.iter())
                        .map(|q| q.0)
                        .collect();
                    for qid in 0..nq as u64 {
                        if used.contains(&qid) {
                            continue;
                        }
                        if (sim.qubit(QubitId(qid)) >> shot) & 1 != 0 {
                            return Err(format!(
                                "ancilla {qid} dirty: n={n} clean_top={clean_top} sub={is_sub} shot={shot}"
                            ));
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

/// Differential selftest for the prefix-ctrl-suffix top-clean adders against the
/// baseline (clean_top=0). Confirms value, source restoration, ancilla, and
/// phase for all clean_top, ctrl, both add/sub, across widths. Run via
/// `PREFIX_TOPCLEAN_SELFTEST=1`.
pub(crate) fn prefix_topclean_selftest() -> Result<(), String> {
    use crate::circuit::QubitOrBit;
    fn lcg(s: &mut u64) -> u64 {
        *s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        *s
    }
    fn rand_u256(s: &mut u64, n: usize) -> U256 {
        let mut limbs = [0u64; 4];
        for limb in limbs.iter_mut() {
            *limb = lcg(s);
        }
        let mut v = U256::from_limbs(limbs);
        if n < 256 {
            let mask = (U256::from(1u64) << n) - U256::from(1u64);
            v &= mask;
        }
        v
    }

    // (prefix_len, suffix_len). Total n = pl+sl must be <= 255 (U256 modulus).
    let shapes = [
        (2usize, 1usize), (3, 1), (3, 2), (5, 2), (8, 3), (16, 4),
        (33, 7), (64, 5), (128, 9), (180, 55), (240, 7), (244, 9),
    ];
    for &(pl, sl) in shapes.iter() {
        let n = pl + sl; // total a width
        let max_clean_top = pl.saturating_sub(1).min(8);
        for clean_top in 0..=max_clean_top {
            for is_sub in [false, true] {
                for ctrl_val in [0u8, 1u8] {
                    let mut b = B::new();
                    let prefix = b.alloc_qubits(pl);
                    let suffix = b.alloc_qubits(sl);
                    let acc = b.alloc_qubits(n);
                    let carries = b.alloc_qubits(pl); // >= borrowed
                    let scratch = b.alloc_qubit();
                    let ctrl = b.alloc_qubit();
                    if is_sub {
                        cuccaro_sub_fast_prefix_ctrl_suffix_no_cin_topclean(
                            &mut b, &prefix, &suffix, &acc, ctrl, &carries, scratch, clean_top,
                        );
                    } else {
                        cuccaro_add_fast_prefix_ctrl_suffix_no_cin_topclean(
                            &mut b, &prefix, &suffix, &acc, ctrl, &carries, scratch, clean_top,
                        );
                    }
                    let nq = b.next_qubit as usize;
                    let nb = b.next_bit as usize;

                    let prefix_qb: Vec<QubitOrBit> =
                        prefix.iter().map(|&q| QubitOrBit::Qubit(q)).collect();
                    let suffix_qb: Vec<QubitOrBit> =
                        suffix.iter().map(|&q| QubitOrBit::Qubit(q)).collect();
                    let acc_qb: Vec<QubitOrBit> =
                        acc.iter().map(|&q| QubitOrBit::Qubit(q)).collect();
                    let ctrl_qb = [QubitOrBit::Qubit(ctrl)];

                    let mut seed = Shake256::default();
                    seed.update(b"prefix-topclean-selftest");
                    seed.update(&(n as u64).to_le_bytes());
                    seed.update(&(clean_top as u64).to_le_bytes());
                    seed.update(&[is_sub as u8, ctrl_val]);
                    let mut xof = seed.finalize_xof();
                    let mut sim = Simulator::new(nq.max(1), nb.max(1), &mut xof);

                    let mut rs = 0x1234567u64 ^ ((n as u64) << 32) ^ ((clean_top as u64) << 8)
                        ^ ((ctrl_val as u64) << 16);
                    let modulus = U256::from(1u64) << n;
                    let mut a_in = [U256::ZERO; 64];
                    let mut acc_in = [U256::ZERO; 64];
                    for shot in 0..64usize {
                        let av = rand_u256(&mut rs, n);
                        let accv = rand_u256(&mut rs, n);
                        a_in[shot] = av;
                        acc_in[shot] = accv;
                        // In real GCD usage the prefix lanes hold `ctrl AND a_low`
                        // (loaded via ccx(ctrl, sub, gated)); the suffix is read
                        // through ctrl. So gate BOTH by ctrl in the fixture.
                        let pmask = (U256::from(1u64) << pl) - U256::from(1u64);
                        let prefix_load = if ctrl_val == 1 { av & pmask } else { U256::ZERO };
                        sim.set_register(&prefix_qb, prefix_load, shot);
                        sim.set_register(&suffix_qb, av >> pl, shot);
                        sim.set_register(&acc_qb, accv, shot);
                        if ctrl_val == 1 {
                            sim.set_register(&ctrl_qb, U256::from(1u64), shot);
                        }
                    }

                    sim.apply_iter(b.ops.iter());

                    if sim.phase != 0 {
                        return Err(format!(
                            "phase leak: pl={pl} sl={sl} clean_top={clean_top} sub={is_sub} ctrl={ctrl_val} phase={:#x}",
                            sim.phase
                        ));
                    }
                    for shot in 0..64usize {
                        // With both prefix and suffix gated by ctrl: delta = ctrl ? a : 0.
                        let delta = if ctrl_val == 1 { a_in[shot] } else { U256::ZERO };
                        let sum = if is_sub {
                            (acc_in[shot] + modulus - delta) % modulus
                        } else {
                            (acc_in[shot] + delta) % modulus
                        };
                        let got = sim.get_register(&acc_qb, shot);
                        if got != sum {
                            return Err(format!(
                                "value: pl={pl} sl={sl} ct={clean_top} sub={is_sub} ctrl={ctrl_val} shot={shot} got={got:#x} exp={sum:#x}"
                            ));
                        }
                        // prefix (gated by ctrl) and suffix must be restored.
                        let pmask = (U256::from(1u64) << pl) - U256::from(1u64);
                        let prefix_load = if ctrl_val == 1 { a_in[shot] & pmask } else { U256::ZERO };
                        let p_got = sim.get_register(&prefix_qb, shot);
                        let s_got = sim.get_register(&suffix_qb, shot);
                        if p_got != prefix_load || s_got != (a_in[shot] >> pl) {
                            return Err(format!(
                                "source: pl={pl} sl={sl} ct={clean_top} sub={is_sub} ctrl={ctrl_val} shot={shot} p_got={p_got:#x} s_got={s_got:#x}"
                            ));
                        }
                        let used: std::collections::HashSet<u64> = prefix
                            .iter()
                            .chain(suffix.iter())
                            .chain(acc.iter())
                            .chain(std::iter::once(&ctrl))
                            .map(|q| q.0)
                            .collect();
                        for qid in 0..nq as u64 {
                            if used.contains(&qid) {
                                continue;
                            }
                            if (sim.qubit(QubitId(qid)) >> shot) & 1 != 0 {
                                return Err(format!(
                                    "ancilla {qid} dirty: pl={pl} sl={sl} ct={clean_top} sub={is_sub} ctrl={ctrl_val} shot={shot}"
                                ));
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(())
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

/// Fredkin (controlled swap): swap (a, t) if ctrl. Decomposed as CX/CCX/CX.
pub(crate) fn cswap(b: &mut B, ctrl: QubitId, a: QubitId, t: QubitId) {
    if a == t {
        return;
    }
    assert!(
        ctrl != a && ctrl != t,
        "invalid CSWAP with control aliased to swapped wire"
    );
    b.cx(t, a);
    b.ccx(ctrl, a, t);
    b.cx(t, a);
}


/// flag ^= (u < v).  Non-destructive on u and v.
///
/// Uses a MAJ-only carry chain instead of the full sub+add pattern.
/// Identity: u < v iff carry-out of (~u + v) = 1, since
///   ~u + v = (2^n - 1 - u) + v = (v - u) + (2^n - 1)
/// which overflows 2^n iff v - u ≥ 1 iff v > u. We negate u in place,
/// run a forward MAJ sweep over (~u, v, c_in=0), capture u[n-1] (which
/// holds the high carry after the chain), then run the inverse MAJ
/// sweep + un-negate to restore u and v. Cost ≈ 2n CCX, half of the
/// previous sub+add (≈ 4n CCX).

// ═══════════════════════════════════════════════════════════════════════════
//  Primitives for the Kaliski port (qrisp-style)
// ═══════════════════════════════════════════════════════════════════════════

/// 3-controlled X with per-control polarity. Uses a borrowed scratch qubit
/// (must be supplied clean, returns clean).
// Gate: measured (Hmr) uncompute of the mcx3_polar scratch instead of a coherent
// CCX. Default OFF = byte-identical. Saves 1 CCX per controlled-ride MAJ/UMA leaf
// (the CAT4 lever tofprof found in materialized_add/sub_body on the live #1).
fn mcx3_measured_uncompute_enabled() -> bool {
    std::env::var("DIALOG_GCD_CTRL_LOWQ_MEASURED").ok().as_deref() == Some("1")
}

pub(crate) fn mcx3_polar(
    b: &mut B,
    c1: QubitId,
    p1: bool,
    c2: QubitId,
    p2: bool,
    c3: QubitId,
    p3: bool,
    target: QubitId,
    scratch: QubitId,
) {
    if !p1 {
        b.x(c1);
    }
    if !p2 {
        b.x(c2);
    }
    if !p3 {
        b.x(c3);
    }
    b.ccx(c1, c2, scratch);
    b.ccx(scratch, c3, target);
    if mcx3_measured_uncompute_enabled() {
        // Gidney measured AND-uncompute (cf. cuccaro_*_fast / ROUND84_FOLD): scratch
        // holds c1&c2 here; a free Hmr + classical CZ(c1,c2) correction uncomputes it
        // value/phase-exactly and resets scratch to |0>, saving the coherent CCX.
        let m = b.alloc_bit();
        b.hmr(scratch, m);
        b.cz_if(c1, c2, m);
    } else {
        b.ccx(c1, c2, scratch);
    }
    if !p3 {
        b.x(c3);
    }
    if !p2 {
        b.x(c2);
    }
    if !p1 {
        b.x(c1);
    }
}

pub(crate) fn ctrl_maj(b: &mut B, ctrl: QubitId, x: QubitId, y: QubitId, w: QubitId, scratch: QubitId) {
    b.ccx(ctrl, w, y);
    b.ccx(ctrl, w, x);
    mcx3_polar(b, ctrl, true, x, true, y, true, w, scratch);
}

pub(crate) fn ctrl_uma(b: &mut B, ctrl: QubitId, x: QubitId, y: QubitId, w: QubitId, scratch: QubitId) {
    mcx3_polar(b, ctrl, true, x, true, y, true, w, scratch);
    b.ccx(ctrl, w, x);
    b.ccx(ctrl, x, y);
}

pub(crate) fn ctrl_inv_maj(b: &mut B, ctrl: QubitId, x: QubitId, y: QubitId, w: QubitId, scratch: QubitId) {
    mcx3_polar(b, ctrl, true, x, true, y, true, w, scratch);
    b.ccx(ctrl, w, x);
    b.ccx(ctrl, w, y);
}

pub(crate) fn ctrl_inv_uma(b: &mut B, ctrl: QubitId, x: QubitId, y: QubitId, w: QubitId, scratch: QubitId) {
    b.ccx(ctrl, x, y);
    b.ccx(ctrl, w, x);
    mcx3_polar(b, ctrl, true, x, true, y, true, w, scratch);
}

pub(crate) fn cuccaro_add_ctrl_lowq(
    b: &mut B,
    a: &[QubitId],
    acc: &[QubitId],
    ctrl: QubitId,
    c_in: QubitId,
    scratch: QubitId,
) {
    let n = a.len();
    assert_eq!(n, acc.len());
    if n == 0 {
        return;
    }
    if n == 1 {
        b.ccx(ctrl, c_in, acc[0]);
        b.ccx(ctrl, a[0], acc[0]);
        return;
    }

    ctrl_maj(b, ctrl, c_in, acc[0], a[0], scratch);
    for i in 1..n - 1 {
        ctrl_maj(b, ctrl, a[i - 1], acc[i], a[i], scratch);
    }

    b.ccx(ctrl, a[n - 2], acc[n - 1]);
    b.ccx(ctrl, a[n - 1], acc[n - 1]);

    for i in (1..n - 1).rev() {
        ctrl_uma(b, ctrl, a[i - 1], acc[i], a[i], scratch);
    }
    ctrl_uma(b, ctrl, c_in, acc[0], a[0], scratch);
}

pub(crate) fn cuccaro_sub_ctrl_lowq(
    b: &mut B,
    a: &[QubitId],
    acc: &[QubitId],
    ctrl: QubitId,
    c_in: QubitId,
    scratch: QubitId,
) {
    let n = a.len();
    assert_eq!(n, acc.len());
    if n == 0 {
        return;
    }
    if n == 1 {
        b.ccx(ctrl, a[0], acc[0]);
        b.ccx(ctrl, c_in, acc[0]);
        return;
    }

    ctrl_inv_uma(b, ctrl, c_in, acc[0], a[0], scratch);
    for i in 1..n - 1 {
        ctrl_inv_uma(b, ctrl, a[i - 1], acc[i], a[i], scratch);
    }

    b.ccx(ctrl, a[n - 1], acc[n - 1]);
    b.ccx(ctrl, a[n - 2], acc[n - 1]);

    for i in (1..n - 1).rev() {
        ctrl_inv_maj(b, ctrl, a[i - 1], acc[i], a[i], scratch);
    }
    ctrl_inv_maj(b, ctrl, c_in, acc[0], a[0], scratch);
}

pub(crate) fn cucc_add_ctrl_lowq(b: &mut B, a: &[QubitId], acc: &[QubitId], ctrl: QubitId) {
    let c_in = b.alloc_qubit();
    let scratch = b.alloc_qubit();
    cuccaro_add_ctrl_lowq(b, a, acc, ctrl, c_in, scratch);
    b.free(scratch);
    b.free(c_in);
}

pub(crate) fn cucc_sub_ctrl_lowq(b: &mut B, a: &[QubitId], acc: &[QubitId], ctrl: QubitId) {
    let c_in = b.alloc_qubit();
    let scratch = b.alloc_qubit();
    cuccaro_sub_ctrl_lowq(b, a, acc, ctrl, c_in, scratch);
    b.free(scratch);
    b.free(c_in);
}


// ═══════════════════════════════════════════════════════════════════════════
//  Kaliski binary almost-inverse (qrisp-style, standard form)
// ═══════════════════════════════════════════════════════════════════════════
//
// Faithful port of `kaliski_mod_inv` from the qrisp reference at
// `quantum-elliptic-curve-logarithm/src/quantum/ec_arithmetic.py`.
//
// The function computes `v_in := v_in^{-1} mod p` in place, using a
// self-contained scratch region that is zeroed at function exit. Every
// per-iteration ancilla is uncomputed via the `conjugate` pattern or via
// classical invariants (e.g. `a ^= NOT s[0]` at the end of each iteration).
//
// Difference from qrisp: we work in STANDARD form, no Montgomery
// conversion. The final r register holds `-v_orig^{-1} * 2^{2n} mod p`
// instead of the Montgomery version. We compensate via a single in-place
// classical-constant multiplication by K = (2^{-2n}) mod p at function
// end, which gets us back to v_orig^{-1}.
//
// Assumption: v_in is a nonzero element of (Z/p)*. The test harness
// filters out the v_orig = 0 case before calling `build`, so we skip the

