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
    b.ccx(c1, c2, scratch);
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

