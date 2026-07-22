use super::*;

#[inline]
fn maj1_inputs_distinct(a: QubitId, k: QubitId, carry: QubitId, target: QubitId) -> bool {
    a != k && a != carry && a != target && k != carry && k != target && carry != target
}

#[inline]
fn fold_maj1_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_MAJ1").ok().as_deref() == Some("1")
}

fn emit_fold_maj1(b: &mut B, a: QubitId, k: QubitId, carry: QubitId, target: QubitId) {
    debug_assert!(maj1_inputs_distinct(a, k, carry, target));
    b.cx(carry, target);
    b.cx(carry, a);
    b.cx(carry, k);
    b.ccx(a, k, target);
    b.cx(carry, k);
    b.cx(carry, a);
}

fn emit_fold_majority(
    b: &mut B,
    a: QubitId,
    k: QubitId,
    carry: QubitId,
    target: QubitId,
    maj2: bool,
) {
    if fold_maj1_enabled() && maj1_inputs_distinct(a, k, carry, target) {
        emit_fold_maj1(b, a, k, carry, target);
    } else if maj2 {
        b.ccx(a, carry, target);
        b.cx(a, carry);
        b.ccx(k, carry, target);
        b.cx(a, carry);
    } else {
        b.ccx(a, carry, target);
        b.ccx(k, a, target);
        b.ccx(k, carry, target);
    }
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
    if n == 0 {
        return;
    }
    if n == 1 {
        if bit(c, 0) {
            b.cx(ctrl, acc[0]);
        }
        return;
    }

    let borrows = b.alloc_qubits(n - 1);

    // Forward borrow sweep. borrow_{i+1} = majority(!acc_i, k_i, borrow_i),
    // where k_i = ctrl when c_i=1 and 0 otherwise.
    for i in 0..n - 1 {
        let target = borrows[i];
        let borrow_in = if i == 0 { None } else { Some(borrows[i - 1]) };
        if bit(c, i) {
            b.x(acc[i]);
            if let Some(bi) = borrow_in {
                emit_fold_majority(b, acc[i], ctrl, bi, target, false);
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

    // Difference bits: acc_i ^= k_i ^ borrow_i.
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, acc[i]);
        }
        if i > 0 {
            b.cx(borrows[i - 1], acc[i]);
        }
    }

    // Measurement-uncompute borrows in reverse.  For subtraction the post-sum
    // identity is borrow_{i+1} = majority(acc_i_final, k_i, borrow_i).
    for i in (0..n - 1).rev() {
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
    if n == 0 {
        return;
    }
    if n == 1 {
        if bit(c, 0) {
            b.cx(ctrl, acc[0]);
        }
        return;
    }

    let carries = b.alloc_qubits(n - 1);

    // Forward carry sweep. carry_{i+1} = majority(acc_i, k_i, carry_i).
    for i in 0..n - 1 {
        let target = carries[i];
        let carry_in = if i == 0 { None } else { Some(carries[i - 1]) };
        if bit(c, i) {
            if let Some(ci) = carry_in {
                emit_fold_majority(b, acc[i], ctrl, ci, target, false);
            } else {
                b.ccx(acc[i], ctrl, target);
            }
        } else if let Some(ci) = carry_in {
            b.ccx(acc[i], ci, target);
        }
    }

    // Sum bits: acc_i ^= k_i ^ carry_i.
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, acc[i]);
        }
        if i > 0 {
            b.cx(carries[i - 1], acc[i]);
        }
    }

    // Measurement-uncompute carries in reverse.  For addition the post-sum
    // identity is carry_{i+1} = majority(!acc_i_final, k_i, carry_i).
    for i in (0..n - 1).rev() {
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

// ═══════════════════════════════════════════════════════════════════════════
//  Ancilla-light extended-carry constant adders (clean, emit_inverse-safe)
// ═══════════════════════════════════════════════════════════════════════════
//
// These add/subtract a classical constant `c` to an (n+1)-bit accumulator
// `acc_ext` (= n-bit register + a top extension bit), capturing the carry/borrow
// into `acc_ext[n]` — exactly like the load-a-full-(n+1)-register + Cuccaro
// pattern in `add_nbit_const`/`csub_nbit_const`, but the loaded constant register
// is only `n = acc_ext.len() - 1` qubits wide (not n+1). For the round84 Solinas
// constant c = 2^256 - p = 2^32 + 977, which has highest set bit 32 ≪ n, the
// low-n register trivially holds it, and the clean carry-capturing Cuccaro
// (`cuccaro_add/sub_low_to_ext_clean`, X/CX/CCX only) folds the overflow into
// `acc_ext[n]`. This drops the +1-qubit transient of the materialized 257-wide
// `load_const` at the mid-sub peak. All four are measurement-free, so they are
// safe to replay under `emit_inverse`.

/// `acc_ext := (acc_ext + c) mod 2^(n+1)` capturing carry into the top bit.
/// Drop-in value-replacement for `add_nbit_const` when the caller passes an
/// extended (n+1)-wide register and `c < 2^n`.
pub(crate) fn add_nbit_const_extcarry_clean(b: &mut B, acc_ext: &[QubitId], c: U256) {
    add_nbit_const_extcarry_clean_with_cin(b, acc_ext, c, None);
}

/// Same as [`add_nbit_const_extcarry_clean`] but optionally sources the Cuccaro
/// carry-in ancilla from a caller-supplied **clean (|0>) idle** qubit instead of
/// allocating a fresh one. When `borrow_cin = Some(q)`, `q` must be |0> on entry
/// and idle for the duration of this call; it is used as the carry-in slot and
/// returned to |0> (the clean MAJ/UMA sweep restores it). Sourcing the carry-in
/// from an existing live-but-idle lane removes the sole +1 fresh allocation that
/// pins the round84-lowq mid-sub peak at 1308 → 1307. Value-/phase-identical to
/// the fresh-ancilla path (the borrowed qubit plays the identical role).
pub(crate) fn add_nbit_const_extcarry_clean_with_cin(
    b: &mut B,
    acc_ext: &[QubitId],
    c: U256,
    borrow_cin: Option<QubitId>,
) {
    let ext = acc_ext.len();
    debug_assert!(ext >= 1);
    let n = ext - 1;
    let ca = load_const(b, n, c);
    let (c_in, fresh) = match borrow_cin {
        Some(q) => (q, false),
        None => (b.alloc_qubit(), true),
    };
    cuccaro_add_low_to_ext_clean(b, &ca, acc_ext, c_in);
    if fresh {
        b.free(c_in);
    }
    unload_const(b, &ca, c);
}

/// `acc_ext := (acc_ext - c) mod 2^(n+1)` capturing borrow into the top bit.
/// Drop-in value-replacement for `sub_nbit_const`.
pub(crate) fn sub_nbit_const_extcarry_clean(b: &mut B, acc_ext: &[QubitId], c: U256) {
    let ext = acc_ext.len();
    debug_assert!(ext >= 1);
    let n = ext - 1;
    let ca = load_const(b, n, c);
    let c_in = b.alloc_qubit();
    cuccaro_sub_low_to_ext_clean(b, &ca, acc_ext, c_in);
    b.free(c_in);
    unload_const(b, &ca, c);
}

/// Controlled `acc_ext += (ctrl ? c : 0)` (mod 2^(n+1)), carry into top bit.
/// The constant is loaded as `(ctrl ? c : 0)` via CX-from-ctrl, so the
/// unconditional clean adder realizes the controlled add. Drop-in for
/// `cadd_nbit_const`.
pub(crate) fn cadd_nbit_const_extcarry_clean(
    b: &mut B,
    acc_ext: &[QubitId],
    c: U256,
    ctrl: QubitId,
) {
    let ext = acc_ext.len();
    debug_assert!(ext >= 1);
    let n = ext - 1;
    let ca = b.alloc_qubits(n);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, ca[i]);
        }
    }
    let c_in = b.alloc_qubit();
    cuccaro_add_low_to_ext_clean(b, &ca, acc_ext, c_in);
    b.free(c_in);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, ca[i]);
        }
    }
    b.free_vec(&ca);
}

/// Controlled `acc_ext -= (ctrl ? c : 0)` (mod 2^(n+1)), borrow into top bit.
/// Drop-in for `csub_nbit_const`.
pub(crate) fn csub_nbit_const_extcarry_clean(
    b: &mut B,
    acc_ext: &[QubitId],
    c: U256,
    ctrl: QubitId,
) {
    csub_nbit_const_extcarry_clean_with_cin(b, acc_ext, c, ctrl, None);
}

/// Same as [`csub_nbit_const_extcarry_clean`] but optionally sources the Cuccaro
/// borrow-in ancilla from a caller-supplied clean (|0>) idle qubit. See
/// [`add_nbit_const_extcarry_clean_with_cin`] for the borrow contract. This is
/// the peak-binding call inside the round84-lowq mid-sub; borrowing its `c_in`
/// from the idle `a_ovf` lane drops the mid-sub peak 1308 → 1307.
pub(crate) fn csub_nbit_const_extcarry_clean_with_cin(
    b: &mut B,
    acc_ext: &[QubitId],
    c: U256,
    ctrl: QubitId,
    borrow_cin: Option<QubitId>,
) {
    let ext = acc_ext.len();
    debug_assert!(ext >= 1);
    let n = ext - 1;
    let ca = b.alloc_qubits(n);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, ca[i]);
        }
    }
    let (c_in, fresh) = match borrow_cin {
        Some(q) => (q, false),
        None => (b.alloc_qubit(), true),
    };
    cuccaro_sub_low_to_ext_clean(b, &ca, acc_ext, c_in);
    if fresh {
        b.free(c_in);
    }
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, ca[i]);
        }
    }
    b.free_vec(&ca);
}

pub(crate) fn add_nbit_const_direct_uncontrolled_fast(b: &mut B, acc: &[QubitId], c: U256) {
    let ctrl = b.alloc_qubit();
    b.x(ctrl);
    cadd_nbit_const_direct_fast(b, acc, c, ctrl);
    b.x(ctrl);
    b.free(ctrl);
}

pub(crate) fn sub_nbit_const_direct_uncontrolled_fast(b: &mut B, acc: &[QubitId], c: U256) {
    let ctrl = b.alloc_qubit();
    b.x(ctrl);
    csub_nbit_const_direct_fast(b, acc, c, ctrl);
    b.x(ctrl);
    b.free(ctrl);
}

pub(crate) fn add_nbit_const_fast(b: &mut B, acc: &[QubitId], c: U256) {
    if secp_direct_const_arith_enabled() {
        add_nbit_const_direct_uncontrolled_fast(b, acc, c);
        return;
    }
    let n = acc.len();
    let a = load_const(b, n, c);
    add_nbit_qq_fast(b, &a, acc);
    unload_const(b, &a, c);
}

pub(crate) fn sub_nbit_const_fast(b: &mut B, acc: &[QubitId], c: U256) {
    if secp_direct_const_arith_enabled() {
        sub_nbit_const_direct_uncontrolled_fast(b, acc, c);
        return;
    }
    let n = acc.len();
    let a = load_const(b, n, c);
    sub_nbit_qq_fast(b, &a, acc);
    unload_const(b, &a, c);
}

// ═══════════════════════════════════════════════════════════════════════════
//  Modular multiplication
// ═══════════════════════════════════════════════════════════════════════════
//
// Shift-and-add, MSB-to-LSB. `acc += x*y mod p`. Iteration:
//
//     for i from n-1 down to 0:
//         acc := 2*acc mod p
//         if y[i]:  acc := acc + x mod p
//
// For q*q mul, y[i] is a qubit; we implement the conditional add by
// CCX-copying x (gated on y[i]) into a temporary, adding, and
// uncopying. For q*b mul, y[i] is a classical bit and the copy is
// done with CX_if gates.

/// Fast `v := 2*v mod p` using measurement-based Cuccaro.
pub(crate) fn highest_set_bit(c: U256) -> usize {
    let mut hi = 0usize;
    for i in 0..256 {
        if bit(c, i) {
            hi = i;
        }
    }
    hi
}

pub(crate) fn double_carry_trunc_window() -> Option<usize> {
    std::env::var("KAL_DOUBLE_CARRY_TRUNC_W")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|&w| w > 0)
}

/// Carry/borrow-tail truncation window for the pseudomersenne overflow/underflow
/// FOLD adders (the controlled `acc[..LSBS] += c` / `-= c` correction after a
/// raw 256-bit add/sub in the materialized-special apply path). Default OFF.
/// Same idea as `double_carry_trunc_window`: the secp256k1 constant
/// c = 2^32+977 is 7-bit-sparse, so the fold's carry ripple can stop a small
/// window above bit 32. Forward (cadd) and inverse (csub) read the same window,
/// so the reverse apply exactly inverts the forward when no truncation triggers
/// (the regime selected by the co-tuned reroll).
pub(crate) fn fold_carry_trunc_window() -> Option<usize> {
    std::env::var("KAL_FOLD_CARRY_TRUNC_W")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|&w| w > 0)
}

/// Default-OFF lever: realize the per-position-controls majority carry/borrow
/// recurrence in 2 CCX instead of 3, with NO ancilla (and no measurement).
///
/// The 3-CCX block `target ^= maj(acc[i], cin, kc)` =
/// `acc·cin ⊕ kc·acc ⊕ kc·cin` is the genuine 3-distinct-input majority that
/// appears in [`cadd_per_position_controls_trunc`] /
/// [`csub_per_position_controls_trunc`] (the apply-phase fused double/halve
/// fold; per-position controls differ, so the single-`ctrl` De-Morgan AND-temp
/// does NOT apply here). It is exactly equal to
/// `acc·cin ⊕ kc·(acc ⊕ cin)`, which can be emitted as:
///   ccx(acc, cin, target);   // target ^= acc·cin
///   cx(acc, cin);            // cin' = acc ⊕ cin   (transient; FREE)
///   ccx(kc, cin, target);    // target ^= kc·(acc ⊕ cin)
///   cx(acc, cin);            // restore cin        (FREE)
/// = 2 CCX. `cin` is a borrow/carry ancilla that is read only at this position
/// (the next position reads `target`, not `cin`); it is restored before the
/// position completes, so the later sum-bit CX and the measurement-uncompute
/// (which read the *restored* `cin`) are untouched. Pure CCX/CX ⇒ no phase.
pub(crate) fn perpos_maj2_enabled() -> bool {
    std::env::var("DIALOG_GCD_PERPOS_MAJ2").ok().as_deref() == Some("1")
}

pub(crate) fn fold_maj2_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_MAJ2").ok().as_deref() == Some("1")
}

fn borrowed_const_fold_carries(
    b: &mut B,
    need: usize,
    borrowed: &[QubitId],
) -> (Vec<QubitId>, Vec<QubitId>) {
    let borrowed_len = borrowed.len().min(need);
    let owned = b.alloc_qubits(need - borrowed_len);
    let mut carries = Vec::with_capacity(need);
    carries.extend_from_slice(&borrowed[..borrowed_len]);
    carries.extend_from_slice(&owned);
    (carries, owned)
}

/// Carry-tail-truncated controlled add of a sparse classical constant.
///
/// Identical arithmetic to [`cadd_nbit_const_direct_fast`] except the forward
/// carry ripple (and the matching measurement-uncompute) is stopped `window`
/// bits above the constant's highest set bit `hi`. Carries `> hi + window`
/// are assumed 0; the corresponding high sum bits keep their input value.
/// This is exact unless a carry generated at/below `hi` propagates through an
/// unbroken run of `window + 1` ones in `acc` above `hi` — probability
/// ~2^-(window+1) per call for random `acc`. The carries `[0 ..= last]` follow
/// the exact same recurrence and post-sum identity as the full adder, so they
/// are returned cleanly to 0 (no phase / ancilla garbage); only the high sum
/// value is approximate.
pub(crate) fn cadd_nbit_const_direct_trunc_fast(
    b: &mut B,
    acc: &[QubitId],
    c: U256,
    ctrl: QubitId,
    window: usize,
) {
    cadd_nbit_const_direct_trunc_fast_borrowed_carries(b, acc, c, ctrl, window, &[]);
}

pub(crate) fn cadd_nbit_const_direct_trunc_fast_borrowed_carries(
    b: &mut B,
    acc: &[QubitId],
    c: U256,
    ctrl: QubitId,
    window: usize,
    borrowed_carries: &[QubitId],
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

    let hi = highest_set_bit(c);
    let last = core::cmp::min(n - 2, hi.saturating_add(window));
    let maj2 = fold_maj2_enabled();
    let (carries, owned_carries) = borrowed_const_fold_carries(b, last + 1, borrowed_carries);

    // Forward carry sweep, truncated at `last`. carry_{i+1} = maj(acc_i, k_i, carry_i).
    for i in 0..=last {
        let target = carries[i];
        let carry_in = if i == 0 { None } else { Some(carries[i - 1]) };
        if bit(c, i) {
            if let Some(ci) = carry_in {
                emit_fold_majority(b, acc[i], ctrl, ci, target, maj2);
            } else {
                b.ccx(acc[i], ctrl, target);
            }
        } else if let Some(ci) = carry_in {
            b.ccx(acc[i], ci, target);
        }
    }

    // Sum bits: acc_i ^= k_i ^ carry_{i-1}; carries above `last` are 0.
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, acc[i]);
        }
        if i > 0 && i - 1 <= last {
            b.cx(carries[i - 1], acc[i]);
        }
    }

    // Measurement-uncompute carries in reverse (same identity as the full adder).
    for i in (0..=last).rev() {
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

    b.free_vec(&owned_carries);
}

/// Carry-tail-truncated controlled subtract of a sparse classical constant.
/// Borrow analogue of [`cadd_nbit_const_direct_trunc_fast`]; the inverse used
/// by the apply-phase modular halve so that halve exactly inverts double when
/// neither truncation triggers (the regime selected by the co-tuned reroll).
pub(crate) fn csub_nbit_const_direct_trunc_fast(
    b: &mut B,
    acc: &[QubitId],
    c: U256,
    ctrl: QubitId,
    window: usize,
) {
    csub_nbit_const_direct_trunc_fast_borrowed_carries(b, acc, c, ctrl, window, &[]);
}

pub(crate) fn csub_nbit_const_direct_trunc_fast_borrowed_carries(
    b: &mut B,
    acc: &[QubitId],
    c: U256,
    ctrl: QubitId,
    window: usize,
    borrowed_carries: &[QubitId],
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

    let hi = highest_set_bit(c);
    let last = core::cmp::min(n - 2, hi.saturating_add(window));
    let maj2 = fold_maj2_enabled();
    let (borrows, owned_borrows) = borrowed_const_fold_carries(b, last + 1, borrowed_carries);

    // Forward borrow sweep, truncated at `last`.
    for i in 0..=last {
        let target = borrows[i];
        let borrow_in = if i == 0 { None } else { Some(borrows[i - 1]) };
        if bit(c, i) {
            b.x(acc[i]);
            if let Some(bi) = borrow_in {
                emit_fold_majority(b, acc[i], ctrl, bi, target, maj2);
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

    // Difference bits: acc_i ^= k_i ^ borrow_{i-1}; borrows above `last` are 0.
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, acc[i]);
        }
        if i > 0 && i - 1 <= last {
            b.cx(borrows[i - 1], acc[i]);
        }
    }

    // Measurement-uncompute borrows in reverse (same identity as the full sub).
    for i in (0..=last).rev() {
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

    b.free_vec(&owned_borrows);
}

fn special_fold_park_low_carries() -> usize {
    std::env::var("DIALOG_GCD_SPECIAL_FOLD_PARK_LOW_CARRIES")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0)
}

fn special_fold_park_low_carries_at_step(step: Option<usize>) -> usize {
    let mapped = step.and_then(|step| {
        let map =
            std::env::var("DIALOG_GCD_SPECIAL_FOLD_PARK_LOW_CARRIES_STEP_MAP").ok()?;
        map.split(',').rev().find_map(|entry| {
            let (raw_step, raw_value) = entry.trim().split_once(':')?;
            if raw_step.trim().parse::<usize>().ok()? != step {
                return None;
            }
            raw_value.trim().parse::<usize>().ok()
        })
    });
    mapped.unwrap_or_else(special_fold_park_low_carries)
}

fn cconst_nbit_direct_trunc_fast_parked(
    b: &mut B,
    acc: &[QubitId],
    c: U256,
    ctrl: QubitId,
    window: usize,
    park_low: usize,
    is_add: bool,
) {
    let n = acc.len();
    if n <= 1 {
        if n == 1 && bit(c, 0) {
            b.cx(ctrl, acc[0]);
        }
        return;
    }

    let hi = highest_set_bit(c);
    let last = core::cmp::min(n - 2, hi.saturating_add(window));
    let park_low = core::cmp::min(park_low, last.saturating_sub(hi));
    if park_low == 0 {
        if is_add {
            cadd_nbit_const_direct_trunc_fast(b, acc, c, ctrl, window);
        } else {
            csub_nbit_const_direct_trunc_fast(b, acc, c, ctrl, window);
        }
        return;
    }

    let split = last - park_low;
    let maj2 = fold_maj2_enabled();
    let prefix = b.alloc_qubits(split + 1);
    let kctrl = |i: usize| bit(c, i).then_some(ctrl);

    for i in 0..=split {
        let target = prefix[i];
        let carry_in = if i == 0 { None } else { Some(prefix[i - 1]) };
        if is_add {
            if let Some(kc) = kctrl(i) {
                if let Some(ci) = carry_in {
                    emit_fold_majority(b, acc[i], kc, ci, target, maj2);
                } else {
                    b.ccx(acc[i], kc, target);
                }
            } else if let Some(ci) = carry_in {
                b.ccx(acc[i], ci, target);
            }
        } else if let Some(kc) = kctrl(i) {
            b.x(acc[i]);
            if let Some(ci) = carry_in {
                emit_fold_majority(b, acc[i], kc, ci, target, maj2);
            } else {
                b.ccx(acc[i], kc, target);
            }
            b.x(acc[i]);
        } else if let Some(ci) = carry_in {
            b.x(acc[i]);
            b.ccx(acc[i], ci, target);
            b.x(acc[i]);
        }
    }

    for i in 0..=split {
        if let Some(kc) = kctrl(i) {
            b.cx(kc, acc[i]);
        }
        if i > 0 {
            b.cx(prefix[i - 1], acc[i]);
        }
    }

    for i in (0..park_low).rev() {
        let measured = b.alloc_bit();
        b.hmr(prefix[i], measured);
        let carry_in = if i == 0 { None } else { Some(prefix[i - 1]) };
        fold_postsum_carry_phase_uncompute(
            b,
            acc,
            kctrl(i),
            carry_in,
            measured,
            i,
            is_add,
        );
        b.free(prefix[i]);
    }

    let tail = b.alloc_qubits(park_low);
    let carry = |i: usize| {
        if i <= split {
            prefix[i]
        } else {
            tail[i - split - 1]
        }
    };
    for i in split + 1..=last {
        let target = carry(i);
        let carry_in = carry(i - 1);
        if is_add {
            if let Some(kc) = kctrl(i) {
                emit_fold_majority(b, acc[i], kc, carry_in, target, maj2);
            } else {
                b.ccx(acc[i], carry_in, target);
            }
        } else if let Some(kc) = kctrl(i) {
            b.x(acc[i]);
            emit_fold_majority(b, acc[i], kc, carry_in, target, maj2);
            b.x(acc[i]);
        } else {
            b.x(acc[i]);
            b.ccx(acc[i], carry_in, target);
            b.x(acc[i]);
        }
    }

    for i in split + 1..n {
        if let Some(kc) = kctrl(i) {
            b.cx(kc, acc[i]);
        }
        if i - 1 <= last {
            b.cx(carry(i - 1), acc[i]);
        }
    }

    for i in (split + 1..=last).rev() {
        let measured = b.alloc_bit();
        b.hmr(carry(i), measured);
        fold_postsum_carry_phase_uncompute(
            b,
            acc,
            kctrl(i),
            Some(carry(i - 1)),
            measured,
            i,
            is_add,
        );
        b.free(carry(i));
    }
    drop(tail);

    for i in 0..park_low {
        b.reacquire(prefix[i]);
        let carry_in = if i == 0 { None } else { Some(prefix[i - 1]) };
        fold_postsum_carry_compute(
            b,
            acc,
            kctrl(i),
            carry_in,
            prefix[i],
            i,
            is_add,
        );
    }

    for i in (0..=split).rev() {
        let measured = b.alloc_bit();
        b.hmr(prefix[i], measured);
        let carry_in = if i == 0 { None } else { Some(prefix[i - 1]) };
        fold_postsum_carry_phase_uncompute(
            b,
            acc,
            kctrl(i),
            carry_in,
            measured,
            i,
            is_add,
        );
        b.free(prefix[i]);
    }
}

pub(crate) fn cadd_nbit_const_direct_trunc_fast_releasing_scratch(
    b: &mut B,
    acc: &[QubitId],
    c: U256,
    ctrl: QubitId,
    window: usize,
    releasable_scratch: &[QubitId],
) {
    cadd_nbit_const_direct_trunc_fast_releasing_scratch_at_step(
        b,
        acc,
        c,
        ctrl,
        window,
        releasable_scratch,
        None,
    );
}

pub(crate) fn cadd_nbit_const_direct_trunc_fast_releasing_scratch_at_step(
    b: &mut B,
    acc: &[QubitId],
    c: U256,
    ctrl: QubitId,
    window: usize,
    releasable_scratch: &[QubitId],
    step: Option<usize>,
) {
    let park_low = special_fold_park_low_carries_at_step(step);
    if park_low == 0 || releasable_scratch.is_empty() {
        cadd_nbit_const_direct_trunc_fast_borrowed_carries(
            b,
            acc,
            c,
            ctrl,
            window,
            releasable_scratch,
        );
        return;
    }
    b.free_vec(releasable_scratch);
    cconst_nbit_direct_trunc_fast_parked(b, acc, c, ctrl, window, park_low, true);
    b.reacquire_vec(releasable_scratch);
}

pub(crate) fn csub_nbit_const_direct_trunc_fast_releasing_scratch(
    b: &mut B,
    acc: &[QubitId],
    c: U256,
    ctrl: QubitId,
    window: usize,
    releasable_scratch: &[QubitId],
) {
    csub_nbit_const_direct_trunc_fast_releasing_scratch_at_step(
        b,
        acc,
        c,
        ctrl,
        window,
        releasable_scratch,
        None,
    );
}

pub(crate) fn csub_nbit_const_direct_trunc_fast_releasing_scratch_at_step(
    b: &mut B,
    acc: &[QubitId],
    c: U256,
    ctrl: QubitId,
    window: usize,
    releasable_scratch: &[QubitId],
    step: Option<usize>,
) {
    let park_low = special_fold_park_low_carries_at_step(step);
    if park_low == 0 || releasable_scratch.is_empty() {
        csub_nbit_const_direct_trunc_fast_borrowed_carries(
            b,
            acc,
            c,
            ctrl,
            window,
            releasable_scratch,
        );
        return;
    }
    b.free_vec(releasable_scratch);
    cconst_nbit_direct_trunc_fast_parked(b, acc, c, ctrl, window, park_low, false);
    b.reacquire_vec(releasable_scratch);
}


pub(crate) fn cadd_per_position_controls_trunc(
    b: &mut B,
    acc: &[QubitId],
    controls: &[Option<QubitId>],
    last: usize,
) {
    let n = acc.len();
    debug_assert!(last < n);
    debug_assert!(controls.len() <= n);
    let kctrl = |i: usize| -> Option<QubitId> {
        if i < controls.len() {
            controls[i]
        } else {
            None
        }
    };
    let maj2 = perpos_maj2_enabled();
    let carries = b.alloc_qubits(last + 1);

    // Forward carry sweep, truncated at `last`. carry_i = maj(acc_i, k_i, carry_{i-1}).
    for i in 0..=last {
        let target = carries[i];
        let carry_in = if i == 0 { None } else { Some(carries[i - 1]) };
        if let Some(kc) = kctrl(i) {
            if let Some(ci) = carry_in {
                emit_fold_majority(b, acc[i], kc, ci, target, maj2);
            } else {
                b.ccx(acc[i], kc, target);
            }
        } else if let Some(ci) = carry_in {
            b.ccx(acc[i], ci, target);
        }
    }

    // Sum bits: acc_i ^= k_i ^ carry_{i-1}; carries above `last` are 0.
    for i in 0..n {
        if let Some(kc) = kctrl(i) {
            b.cx(kc, acc[i]);
        }
        if i > 0 && i - 1 <= last {
            b.cx(carries[i - 1], acc[i]);
        }
    }

    // Measurement-uncompute carries in reverse (free; same identity as the adder).
    for i in (0..=last).rev() {
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        let carry_in = if i == 0 { None } else { Some(carries[i - 1]) };
        if let Some(kc) = kctrl(i) {
            b.x(acc[i]);
            if let Some(ci) = carry_in {
                b.cz_if(acc[i], kc, m);
                b.cz_if(acc[i], ci, m);
                b.x(acc[i]);
                b.cz_if(kc, ci, m);
            } else {
                b.cz_if(acc[i], kc, m);
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

pub(crate) fn csub_per_position_controls_trunc(
    b: &mut B,
    acc: &[QubitId],
    controls: &[Option<QubitId>],
    last: usize,
) {
    let n = acc.len();
    debug_assert!(last < n);
    debug_assert!(controls.len() <= n);
    let kctrl = |i: usize| -> Option<QubitId> {
        if i < controls.len() {
            controls[i]
        } else {
            None
        }
    };
    let maj2 = perpos_maj2_enabled();
    let borrows = b.alloc_qubits(last + 1);

    // Forward borrow sweep, truncated at `last`.
    for i in 0..=last {
        let target = borrows[i];
        let borrow_in = if i == 0 { None } else { Some(borrows[i - 1]) };
        if let Some(kc) = kctrl(i) {
            b.x(acc[i]);
            if let Some(bi) = borrow_in {
                emit_fold_majority(b, acc[i], kc, bi, target, maj2);
            } else {
                b.ccx(acc[i], kc, target);
            }
            b.x(acc[i]);
        } else if let Some(bi) = borrow_in {
            b.x(acc[i]);
            b.ccx(acc[i], bi, target);
            b.x(acc[i]);
        }
    }

    // Difference bits: acc_i ^= k_i ^ borrow_{i-1}; borrows above `last` are 0.
    for i in 0..n {
        if let Some(kc) = kctrl(i) {
            b.cx(kc, acc[i]);
        }
        if i > 0 && i - 1 <= last {
            b.cx(borrows[i - 1], acc[i]);
        }
    }

    // Measurement-uncompute borrows in reverse (free; same identity as the sub).
    for i in (0..=last).rev() {
        let m = b.alloc_bit();
        b.hmr(borrows[i], m);
        let borrow_in = if i == 0 { None } else { Some(borrows[i - 1]) };
        if let Some(kc) = kctrl(i) {
            if let Some(bi) = borrow_in {
                b.cz_if(acc[i], kc, m);
                b.cz_if(acc[i], bi, m);
                b.cz_if(kc, bi, m);
            } else {
                b.cz_if(acc[i], kc, m);
            }
        } else if let Some(bi) = borrow_in {
            b.cz_if(acc[i], bi, m);
        }
    }

    b.free_vec(&borrows);
}

/// Default-OFF lever for the apply-phase fused double_y / halve_y fold ripple.
///
/// The fused fold `y ±= δ = c·e + 2c·d` (`c = 2^32+977`) has per-position
/// controls only at positions ≤ `hi_delta = 33`; positions `(33, last]` are a
/// pure carry/borrow PROPAGATION tail (constant bit 0). The baseline keeps all
/// eight fold ancillae (`e,d,h,xed,eord,n10` + the two overflow holders) live
/// for the WHOLE ripple — including across the wide high tail, which is the
/// double_y/halve_y high-water (`floor + 8 + 34 + W`, `W = KAL_DOUBLE_CARRY_TRUNC_W`).
///
/// This lever frees the FOUR purely-`e,d`-derived controls (`h,xed,eord,n10`)
/// after the active region `[0..=hi]` and before the high tail, then recomputes
/// them (cheap: free CX from `e,d`, plus one AND for `h`) just for the carry
/// uncompute pass. Net the high-tail high-water drops from `+8` to `+4` ancillae
/// (the two overflow holders `e,d` remain, plus the caller's two `ovf` qubits),
/// i.e. the fold floor falls by 4 qubits, value/phase-EXACT (identical arithmetic
/// and identical truncation `last`; only the ancilla lifetime is tightened).
pub(crate) fn fold_freed_tail_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_FREED_TAIL").ok().as_deref() == Some("1")
}

/// e,d-extension of the freed-tail lever (HYP-6 §4a). When ON (and the freed-tail
/// itself is ON), the fused-fold ripple ALSO releases the two base controls `e,d`
/// across the wide high tail — not just the four `e,d`-derived controls
/// (`h,xed,eord,n10`). `e,d` are dead as controls in the tail (all their fold
/// positions sit at ≤ `hi_delta = 33`), and both are recomputable from the live
/// overflow lanes via `d = ovf1 & s2`, `e = ovf1 ^ d ^ ovf2` (the SAME relation
/// holds in both the forward double_y fold and the inverse halve_y fold — see the
/// dispatch sites). Freeing `e,d` too drops the wide-tail high-water from `+4` to
/// `+2` ancillae, i.e. the fold floor falls a further 2 qubits (1220 → 1218 at
/// W=19), value/phase-EXACT (identical arithmetic; only ancilla lifetime tightens;
/// cost = a handful of CX + 1 CCX/call to re-derive `d`). Default OFF ⇒ the
/// freed-tail path is byte-identical to before this lever existed.
pub(crate) fn fold_freed_tail_ed_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_FREED_TAIL_ED")
        .ok()
        .as_deref()
        == Some("1")
}

/// Reuse four future-zero low-carry slots as the derived fold controls
/// (`h,xed,eord,n10`) while processing the sparse constant region. The original
/// control qubits are released before the 34-qubit low-carry lane is allocated
/// and restored only after carries 12..33 have been uncomputed. This is
/// value/phase exact and removes four qubits from the fused-fold high-water.
pub(crate) fn fold_host_derived_controls_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_HOST_DERIVED_CONTROLS")
        .ok()
        .as_deref()
        == Some("1")
}

fn fold_host_n10_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_HOST_N10")
        .ok()
        .as_deref()
        == Some("1")
}

fn fold_host_h_n10_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_HOST_H_N10")
        .ok()
        .as_deref()
        == Some("1")
}

fn fold_host_h_xed_n10_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_HOST_H_XED_N10")
        .ok()
        .as_deref()
        == Some("1")
}

fn fold_host_e_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_HOST_E")
        .ok()
        .as_deref()
        == Some("1")
}

fn fold_host_d_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_HOST_D")
        .ok()
        .as_deref()
        == Some("1")
}

/// Per-call carry window override for the FUSED FOLD only (`double_y`/`halve_y`),
/// decoupled from the GCD-walk's `KAL_DOUBLE_CARRY_TRUNC_W`. When set it caps the
/// fold ripple at `hi_delta + W_fold`; unset = inherit the GCD-walk window
/// (byte-identical base). Lowering it shrinks the fold high-water 1-for-1
/// (`floor + 42 + W_fold`) at the cost of a slightly higher fold-carry-escape
/// truncation rate (the same FS-island hazard class the shared window already
/// carries — see KAL_DOUBLE_CARRY_TRUNC_W).
pub(crate) fn fold_only_carry_trunc_window() -> Option<usize> {
    std::env::var("DIALOG_GCD_FOLD_CARRY_TRUNC_W")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|&w| w > 0)
}

fn fold_park_low_carries() -> usize {
    std::env::var("DIALOG_GCD_FOLD_PARK_LOW_CARRIES")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0)
}

pub(crate) fn fold_park_low_carries_at_step(step: Option<usize>) -> usize {
    let mapped = step.and_then(|step| {
        let map = std::env::var("DIALOG_GCD_FOLD_PARK_LOW_CARRIES_STEP_MAP").ok()?;
        map.split(',').rev().find_map(|entry| {
            let (raw_step, raw_value) = entry.trim().split_once(':')?;
            if raw_step.trim().parse::<usize>().ok()? != step {
                return None;
            }
            raw_value.trim().parse::<usize>().ok()
        })
    });
    mapped.unwrap_or_else(fold_park_low_carries)
}

pub(crate) fn fold_stream_controls_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_STREAM_CONTROLS")
        .ok()
        .as_deref()
        == Some("1")
        && fold_park_low_carries() >= 12
}

fn fold_host_streamed_control_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_HOST_STREAMED_CONTROL")
        .ok()
        .as_deref()
        == Some("1")
        && fold_park_low_carries() >= 13
}

fn fold_host_e_top_carry_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_HOST_E_TOP_CARRY")
        .ok()
        .as_deref()
        == Some("1")
}

fn fold_host_d_carry12_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_HOST_D_CARRY12")
        .ok()
        .as_deref()
        == Some("1")
        && fold_park_low_carries() >= 14
}

fn fold_host_ovf2_carry13_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_HOST_OVF2_CARRY13")
        .ok()
        .as_deref()
        == Some("1")
        && fold_park_low_carries() >= 15
}

fn fold_free_first_high_carry_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_FREE_FIRST_HIGH_CARRY")
        .ok()
        .as_deref()
        == Some("1")
}

fn fold_stream_profile_phase(b: &mut B, add_phase: &'static str, sub_phase: &'static str, is_add: bool) {
    if std::env::var("DIALOG_GCD_FOLD_PROFILE_PHASES").ok().as_deref() == Some("1") {
        b.set_phase(if is_add { add_phase } else { sub_phase });
    }
}

fn fold_postsum_carry_phase_uncompute(
    b: &mut B,
    acc: &[QubitId],
    kctrl: Option<QubitId>,
    carry_in: Option<QubitId>,
    measured: BitId,
    i: usize,
    is_add: bool,
) {
    if is_add {
        if let Some(kc) = kctrl {
            b.x(acc[i]);
            if let Some(ci) = carry_in {
                b.cz_if(acc[i], kc, measured);
                b.cz_if(acc[i], ci, measured);
                b.x(acc[i]);
                b.cz_if(kc, ci, measured);
            } else {
                b.cz_if(acc[i], kc, measured);
                b.x(acc[i]);
            }
        } else if let Some(ci) = carry_in {
            b.x(acc[i]);
            b.cz_if(acc[i], ci, measured);
            b.x(acc[i]);
        }
    } else if let Some(kc) = kctrl {
        if let Some(ci) = carry_in {
            b.cz_if(acc[i], kc, measured);
            b.cz_if(acc[i], ci, measured);
            b.cz_if(kc, ci, measured);
        } else {
            b.cz_if(acc[i], kc, measured);
        }
    } else if let Some(ci) = carry_in {
        b.cz_if(acc[i], ci, measured);
    }
}

fn fold_postsum_carry_compute(
    b: &mut B,
    acc: &[QubitId],
    kctrl: Option<QubitId>,
    carry_in: Option<QubitId>,
    target: QubitId,
    i: usize,
    is_add: bool,
) {
    if is_add {
        if let Some(kc) = kctrl {
            b.x(acc[i]);
            if let Some(ci) = carry_in {
                emit_fold_majority(
                    b,
                    acc[i],
                    kc,
                    ci,
                    target,
                    perpos_maj2_enabled(),
                );
                b.x(acc[i]);
            } else {
                b.ccx(acc[i], kc, target);
                b.x(acc[i]);
            }
        } else if let Some(ci) = carry_in {
            b.x(acc[i]);
            b.ccx(acc[i], ci, target);
            b.x(acc[i]);
        }
    } else if let Some(kc) = kctrl {
        if let Some(ci) = carry_in {
            emit_fold_majority(
                b,
                acc[i],
                kc,
                ci,
                target,
                perpos_maj2_enabled(),
            );
        } else {
            b.ccx(acc[i], kc, target);
        }
    } else if let Some(ci) = carry_in {
        b.ccx(acc[i], ci, target);
    }
}

fn fold_presum_carry_compute_and_sum(
    b: &mut B,
    acc: &[QubitId],
    kctrl: Option<QubitId>,
    carry_in: Option<QubitId>,
    target: QubitId,
    i: usize,
    is_add: bool,
    maj2: bool,
) {
    if is_add {
        if let Some(kc) = kctrl {
            if let Some(ci) = carry_in {
                emit_fold_majority(b, acc[i], kc, ci, target, maj2);
            } else {
                b.ccx(acc[i], kc, target);
            }
        } else if let Some(ci) = carry_in {
            b.ccx(acc[i], ci, target);
        }
    } else if let Some(kc) = kctrl {
        b.x(acc[i]);
        if let Some(ci) = carry_in {
            emit_fold_majority(b, acc[i], kc, ci, target, maj2);
        } else {
            b.ccx(acc[i], kc, target);
        }
        b.x(acc[i]);
    } else if let Some(ci) = carry_in {
        b.x(acc[i]);
        b.ccx(acc[i], ci, target);
        b.x(acc[i]);
    }
    if let Some(kc) = kctrl {
        b.cx(kc, acc[i]);
    }
    if let Some(ci) = carry_in {
        b.cx(ci, acc[i]);
    }
}

/// Build the secp256k1 fold per-position control vector `δ = c·e + 2c·d`
/// (`c = 2^32+977`) from the base controls `e,d` and the four derived controls
/// `h = e&d`, `xed = e^d`, `eord = e|d`, `n10 = ¬e&d`. Shared by the baseline
/// fused double_y/halve_y and the freed-tail lever so the arithmetic is identical.
pub(crate) fn secp_fold_controls(
    e: QubitId,
    d: QubitId,
    h: QubitId,
    xed: QubitId,
    eord: QubitId,
    n10: QubitId,
    hi_delta: usize,
    hi_c: usize,
) -> Vec<Option<QubitId>> {
    let mut controls: Vec<Option<QubitId>> = vec![None; hi_delta + 1];
    controls[0] = Some(e);
    controls[1] = Some(d);
    controls[4] = Some(e);
    controls[5] = Some(d);
    controls[6] = Some(e);
    controls[7] = Some(xed);
    controls[8] = Some(eord);
    controls[9] = Some(eord);
    controls[10] = Some(n10);
    controls[11] = Some(h);
    controls[hi_c] = Some(e); // bit 32
    controls[hi_delta] = Some(d); // bit 33
    controls
}

/// Freed-tail fold ripple (gated by [`fold_freed_tail_enabled`]). Value/phase
/// identical to `cadd_per_position_controls_trunc(acc, secp_fold_controls(...),
/// last)` but the four `e,d`-derived controls (`h,xed,eord,n10`) are released
/// before the wide high tail and recomputed only for the carry uncompute pass,
/// dropping the high-tail high-water by 4 ancillae. `is_add=false` runs the
/// borrow (subtract) variant for halve_y. `e`,`d` are read-only here; the caller
/// owns `h,xed,eord,n10` allocation/free — this routine consumes them via `free`
/// and the caller must NOT free them again (it re-derives `xed,eord,n10` from a
/// fresh alloc on return is NOT needed: this fn fully owns their lifetime).
pub(crate) fn fold_ripple_freed_tail(
    b: &mut B,
    acc: &[QubitId],
    e: QubitId,
    d: QubitId,
    h: QubitId,
    xed: QubitId,
    eord: QubitId,
    n10: QubitId,
    last: usize,
    is_add: bool,
) {
    // Without the e,d-extension `e,d` are held live across the whole ripple.
    fold_ripple_freed_tail_ed(
        b, acc, e, d, h, xed, eord, n10, None, None, last, is_add,
    );
}

/// Low-qubit fused-fold ripple that never materializes the four derived controls
/// simultaneously. A single ancilla walks through xed, eord, n10, and h in both
/// the forward and reverse low-carry sweeps.
pub(crate) fn fold_ripple_freed_tail_ed_streamed(
    b: &mut B,
    acc: &[QubitId],
    e: QubitId,
    d: QubitId,
    ed: Option<(QubitId, QubitId, QubitId)>,
    park_low: usize,
    last: usize,
    is_add: bool,
) {
    let free_ed = ed.is_some() && fold_freed_tail_ed_enabled();
    let n = acc.len();
    let hi_delta = 33usize;
    debug_assert!(last < n);
    debug_assert!(last > hi_delta, "freed-tail requires a nonempty high tail");
    let park_low = core::cmp::min(park_low, hi_delta);
    assert!(
        park_low >= 12,
        "streamed fold controls require at least 12 parked carries"
    );
    let host_streamed = fold_host_streamed_control_enabled();
    let host_e_top = free_ed && fold_host_e_top_carry_enabled();
    let host_d_carry12 =
        host_e_top && fold_host_d_carry12_enabled();
    let host_ovf2_carry13 =
        host_d_carry12 && fold_host_ovf2_carry13_enabled();
    let maj2 = perpos_maj2_enabled();
    let kctrl = |i: usize| match i {
        0 | 4 | 6 | 32 => Some(e),
        1 | 5 | 33 => Some(d),
        _ => None,
    };
    fold_stream_profile_phase(
        b,
        "dialog_gcd_streamed_double_active",
        "dialog_gcd_streamed_halve_active",
        is_add,
    );
    let low_chain_last = if host_e_top {
        hi_delta - 1
    } else {
        hi_delta
    };
    let mut low = b.alloc_qubits(
        low_chain_last + 1
            - usize::from(host_d_carry12)
            - usize::from(host_ovf2_carry13),
    );
    if host_d_carry12 {
        low.insert(12, d);
    }
    if host_ovf2_carry13 {
        let (_, ovf2, _) = ed.expect("host_ovf2_carry13 implies ed is Some");
        low.insert(13, ovf2);
        let (ovf1, _, _) = ed.expect("host_ovf2_carry13 implies ed is Some");
        b.cx(ovf1, ovf2);
        b.cx(d, ovf2);
        b.cx(e, ovf2);
    }
    let streamed_slot = low[park_low - 1];

    for i in 0..7 {
        let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
        fold_presum_carry_compute_and_sum(
            b,
            acc,
            kctrl(i),
            carry_in,
            low[i],
            i,
            is_add,
            maj2,
        );
    }
    let streamed = if host_streamed {
        streamed_slot
    } else {
        b.alloc_qubit()
    };
    b.cx(e, streamed);
    b.cx(d, streamed);
    fold_presum_carry_compute_and_sum(
        b,
        acc,
        Some(streamed),
        Some(low[6]),
        low[7],
        7,
        is_add,
        maj2,
    );
    b.ccx(e, d, streamed);
    for i in 8..10 {
        fold_presum_carry_compute_and_sum(
            b,
            acc,
            Some(streamed),
            Some(low[i - 1]),
            low[i],
            i,
            is_add,
            maj2,
        );
    }
    b.cx(e, streamed);
    fold_presum_carry_compute_and_sum(
        b,
        acc,
        Some(streamed),
        Some(low[9]),
        low[10],
        10,
        is_add,
        maj2,
    );
    b.cx(d, streamed);
    fold_presum_carry_compute_and_sum(
        b,
        acc,
        Some(streamed),
        Some(low[10]),
        low[11],
        11,
        is_add,
        maj2,
    );
    if host_streamed {
        b.ccx(e, d, streamed);
    }
    if host_d_carry12 {
        let (ovf1, _, s2) = ed.expect("host_d_carry12 implies ed is Some");
        b.ccx(ovf1, s2, d);
    }
    for i in 12..=low_chain_last {
        fold_presum_carry_compute_and_sum(
            b,
            acc,
            kctrl(i),
            Some(low[i - 1]),
            low[i],
            i,
            is_add,
            maj2,
        );
    }

    let free_first_high_carry = fold_free_first_high_carry_enabled()
        && park_low < low_chain_last
        && !(host_d_carry12 && park_low == 12)
        && !(host_ovf2_carry13 && park_low == 13);
    if free_first_high_carry {
        let m = b.alloc_bit();
        b.hmr(low[park_low], m);
        fold_postsum_carry_phase_uncompute(
            b,
            acc,
            kctrl(park_low),
            Some(low[park_low - 1]),
            m,
            park_low,
            is_add,
        );
        b.free(low[park_low]);
    }

    for i in (12..park_low).rev() {
        let m = b.alloc_bit();
        b.hmr(low[i], m);
        fold_postsum_carry_phase_uncompute(
            b,
            acc,
            kctrl(i),
            Some(low[i - 1]),
            m,
            i,
            is_add,
        );
        b.free(low[i]);
    }
    if host_d_carry12 {
        let (ovf1, _, s2) = ed.expect("host_d_carry12 implies ed is Some");
        b.reacquire(d);
        b.ccx(ovf1, s2, d);
    }
    if host_streamed {
        b.reacquire(streamed);
        b.ccx(e, d, streamed);
    }
    let m11 = b.alloc_bit();
    b.hmr(low[11], m11);
    fold_postsum_carry_phase_uncompute(
        b,
        acc,
        Some(streamed),
        Some(low[10]),
        m11,
        11,
        is_add,
    );
    b.free(low[11]);
    b.cx(d, streamed);
    let m10 = b.alloc_bit();
    b.hmr(low[10], m10);
    fold_postsum_carry_phase_uncompute(
        b,
        acc,
        Some(streamed),
        Some(low[9]),
        m10,
        10,
        is_add,
    );
    b.free(low[10]);
    b.cx(e, streamed);
    for i in (8..10).rev() {
        let m = b.alloc_bit();
        b.hmr(low[i], m);
        fold_postsum_carry_phase_uncompute(
            b,
            acc,
            Some(streamed),
            Some(low[i - 1]),
            m,
            i,
            is_add,
        );
        b.free(low[i]);
    }
    b.ccx(e, d, streamed);
    let m7 = b.alloc_bit();
    b.hmr(low[7], m7);
    fold_postsum_carry_phase_uncompute(
        b,
        acc,
        Some(streamed),
        Some(low[6]),
        m7,
        7,
        is_add,
    );
    b.free(low[7]);
    b.cx(d, streamed);
    b.cx(e, streamed);
    b.free(streamed);
    for i in (0..7).rev() {
        let m = b.alloc_bit();
        b.hmr(low[i], m);
        let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
        fold_postsum_carry_phase_uncompute(
            b,
            acc,
            kctrl(i),
            carry_in,
            m,
            i,
            is_add,
        );
        b.free(low[i]);
    }

    if host_ovf2_carry13 {
        let (ovf1, ovf2, _) = ed.expect("host_ovf2_carry13 implies ed is Some");
        b.reacquire(ovf2);
        b.cx(ovf1, ovf2);
        b.cx(d, ovf2);
        b.cx(e, ovf2);
    }

    if host_e_top {
        let (ovf1, ovf2, _) = ed.expect("host_e_top implies ed is Some");
        b.cx(ovf1, e);
        b.cx(d, e);
        b.cx(ovf2, e);
        fold_presum_carry_compute_and_sum(
            b,
            acc,
            Some(d),
            Some(low[hi_delta - 1]),
            e,
            hi_delta,
            is_add,
            maj2,
        );
    }

    if free_ed {
        let (ovf1, ovf2, s2) = ed.expect("free_ed implies ed is Some");
        if !host_e_top {
            b.cx(ovf1, e);
            b.cx(d, e);
            b.cx(ovf2, e);
            b.free(e);
        }
        let md = b.alloc_bit();
        b.hmr(d, md);
        b.cz_if(ovf1, s2, md);
        b.free(d);
    }

    fold_stream_profile_phase(
        b,
        "dialog_gcd_streamed_double_tail",
        "dialog_gcd_streamed_halve_tail",
        is_add,
    );
    let tail_len = last - hi_delta;
    let tail = b.alloc_qubits(tail_len);
    let cw = |i: usize| -> QubitId {
        if i < hi_delta {
            low[i]
        } else if i == hi_delta {
            if host_e_top {
                e
            } else {
                low[i]
            }
        } else {
            tail[i - hi_delta - 1]
        }
    };
    for i in hi_delta + 1..=last {
        if is_add {
            b.ccx(acc[i], cw(i - 1), cw(i));
        } else {
            b.x(acc[i]);
            b.ccx(acc[i], cw(i - 1), cw(i));
            b.x(acc[i]);
        }
    }
    for i in hi_delta + 1..n {
        if i - 1 <= last {
            b.cx(cw(i - 1), acc[i]);
        }
    }
    for i in (hi_delta + 1..=last).rev() {
        let m = b.alloc_bit();
        b.hmr(cw(i), m);
        let carry_in = cw(i - 1);
        if is_add {
            b.x(acc[i]);
            b.cz_if(acc[i], carry_in, m);
            b.x(acc[i]);
        } else {
            b.cz_if(acc[i], carry_in, m);
        }
        b.free(cw(i));
    }
    drop(tail);

    fold_stream_profile_phase(
        b,
        "dialog_gcd_streamed_double_reverse",
        "dialog_gcd_streamed_halve_reverse",
        is_add,
    );
    if free_ed {
        let (ovf1, ovf2, s2) = ed.expect("free_ed implies ed is Some");
        b.reacquire(d);
        b.ccx(ovf1, s2, d);
        if host_e_top {
            let m_top = b.alloc_bit();
            b.hmr(e, m_top);
            fold_postsum_carry_phase_uncompute(
                b,
                acc,
                Some(d),
                Some(low[hi_delta - 1]),
                m_top,
                hi_delta,
                is_add,
            );
            b.free(e);
        }
        b.reacquire(e);
        b.cx(ovf1, e);
        b.cx(d, e);
        b.cx(ovf2, e);
    }
    if host_ovf2_carry13 {
        let (ovf1, ovf2, _) = ed.expect("host_ovf2_carry13 implies ed is Some");
        b.cx(ovf1, ovf2);
        b.cx(d, ovf2);
        b.cx(e, ovf2);
    }

    for i in 0..park_low {
        if !(host_d_carry12 && i == 12)
            && !(host_ovf2_carry13 && i == 13)
        {
            b.reacquire(low[i]);
        }
    }
    for i in 0..7 {
        let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
        fold_postsum_carry_compute(
            b,
            acc,
            kctrl(i),
            carry_in,
            low[i],
            i,
            is_add,
        );
    }
    let streamed = if host_streamed {
        streamed_slot
    } else {
        b.alloc_qubit()
    };
    b.cx(e, streamed);
    b.cx(d, streamed);
    fold_postsum_carry_compute(
        b,
        acc,
        Some(streamed),
        Some(low[6]),
        low[7],
        7,
        is_add,
    );
    b.ccx(e, d, streamed);
    for i in 8..10 {
        fold_postsum_carry_compute(
            b,
            acc,
            Some(streamed),
            Some(low[i - 1]),
            low[i],
            i,
            is_add,
        );
    }
    b.cx(e, streamed);
    fold_postsum_carry_compute(
        b,
        acc,
        Some(streamed),
        Some(low[9]),
        low[10],
        10,
        is_add,
    );
    b.cx(d, streamed);
    fold_postsum_carry_compute(
        b,
        acc,
        Some(streamed),
        Some(low[10]),
        low[11],
        11,
        is_add,
    );
    if host_streamed {
        b.ccx(e, d, streamed);
    }
    if host_d_carry12 {
        let (ovf1, _, s2) = ed.expect("host_d_carry12 implies ed is Some");
        b.ccx(ovf1, s2, d);
    }
    for i in 12..park_low {
        fold_postsum_carry_compute(
            b,
            acc,
            kctrl(i),
            Some(low[i - 1]),
            low[i],
            i,
            is_add,
        );
    }
    if free_first_high_carry {
        b.reacquire(low[park_low]);
        fold_postsum_carry_compute(
            b,
            acc,
            kctrl(park_low),
            Some(low[park_low - 1]),
            low[park_low],
            park_low,
            is_add,
        );
    }

    for i in (12..=low_chain_last).rev() {
        let m = b.alloc_bit();
        b.hmr(low[i], m);
        fold_postsum_carry_phase_uncompute(
            b,
            acc,
            kctrl(i),
            Some(low[i - 1]),
            m,
            i,
            is_add,
        );
        b.free(low[i]);
    }
    if host_d_carry12 {
        let (ovf1, _, s2) = ed.expect("host_d_carry12 implies ed is Some");
        b.reacquire(d);
        b.ccx(ovf1, s2, d);
    }
    if host_ovf2_carry13 {
        let (ovf1, ovf2, _) = ed.expect("host_ovf2_carry13 implies ed is Some");
        b.reacquire(ovf2);
        b.cx(ovf1, ovf2);
        b.cx(d, ovf2);
        b.cx(e, ovf2);
    }
    if host_streamed {
        b.reacquire(streamed);
        b.ccx(e, d, streamed);
    }
    let m11 = b.alloc_bit();
    b.hmr(low[11], m11);
    fold_postsum_carry_phase_uncompute(
        b,
        acc,
        Some(streamed),
        Some(low[10]),
        m11,
        11,
        is_add,
    );
    b.free(low[11]);
    b.cx(d, streamed);
    let m10 = b.alloc_bit();
    b.hmr(low[10], m10);
    fold_postsum_carry_phase_uncompute(
        b,
        acc,
        Some(streamed),
        Some(low[9]),
        m10,
        10,
        is_add,
    );
    b.free(low[10]);
    b.cx(e, streamed);
    for i in (8..10).rev() {
        let m = b.alloc_bit();
        b.hmr(low[i], m);
        fold_postsum_carry_phase_uncompute(
            b,
            acc,
            Some(streamed),
            Some(low[i - 1]),
            m,
            i,
            is_add,
        );
        b.free(low[i]);
    }
    b.ccx(e, d, streamed);
    let m7 = b.alloc_bit();
    b.hmr(low[7], m7);
    fold_postsum_carry_phase_uncompute(
        b,
        acc,
        Some(streamed),
        Some(low[6]),
        m7,
        7,
        is_add,
    );
    b.free(low[7]);
    b.cx(d, streamed);
    b.cx(e, streamed);
    b.free(streamed);
    for i in (0..7).rev() {
        let m = b.alloc_bit();
        b.hmr(low[i], m);
        let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
        fold_postsum_carry_phase_uncompute(
            b,
            acc,
            kctrl(i),
            carry_in,
            m,
            i,
            is_add,
        );
        b.free(low[i]);
    }
    drop(low);
}

/// e,d-extension variant of [`fold_ripple_freed_tail`] (HYP-6 §4a). When
/// `ed = Some((ovf1, ovf2, s2))` AND [`fold_freed_tail_ed_enabled`], `e,d` are
/// additionally released across the wide high tail and recomputed from the live
/// overflow lanes (`d = ovf1 & s2`, `e = ovf1 ^ d ^ ovf2`) for the low uncompute
/// pass, dropping the tail high-water by 2 more ancillae. `ovf1, ovf2, s2` are
/// read-only and must be live & unchanged for the whole call. When `ed = None`
/// (or the knob is OFF) this is byte-identical to the plain freed-tail.
pub(crate) fn fold_ripple_freed_tail_ed(
    b: &mut B,
    acc: &[QubitId],
    e: QubitId,
    d: QubitId,
    h: QubitId,
    xed: QubitId,
    eord: QubitId,
    n10: QubitId,
    ed: Option<(QubitId, QubitId, QubitId)>,
    step: Option<usize>,
    last: usize,
    is_add: bool,
) {
    let configured_park_low = fold_park_low_carries_at_step(step);
    if fold_host_derived_controls_enabled() && configured_park_low <= 7 {
        fold_ripple_freed_tail_ed_hosted(
            b,
            acc,
            e,
            d,
            h,
            xed,
            eord,
            n10,
            ed,
            configured_park_low,
            last,
            is_add,
        );
        return;
    }
    if fold_stream_controls_enabled() && configured_park_low >= 12 {
        b.cx(h, n10);
        b.cx(d, n10);
        b.cx(h, eord);
        b.cx(xed, eord);
        b.cx(d, xed);
        b.cx(e, xed);
        b.free(n10);
        b.free(eord);
        b.free(xed);
        let mh = b.alloc_bit();
        b.hmr(h, mh);
        b.cz_if(e, d, mh);
        b.free(h);
        fold_ripple_freed_tail_ed_streamed(
            b,
            acc,
            e,
            d,
            ed,
            configured_park_low,
            last,
            is_add,
        );
        b.reacquire(h);
        b.ccx(e, d, h);
        b.reacquire(xed);
        b.cx(e, xed);
        b.cx(d, xed);
        b.reacquire(eord);
        b.cx(xed, eord);
        b.cx(h, eord);
        b.reacquire(n10);
        b.cx(d, n10);
        b.cx(h, n10);
        return;
    }

    let free_ed = ed.is_some() && fold_freed_tail_ed_enabled();
    let n = acc.len();
    let hi_delta = 33usize; // highest_set_bit(2^32+977)+1
    let hi_c = 32usize;
    debug_assert!(last < n);
    debug_assert!(last > hi_delta, "freed-tail requires a nonempty high tail");
    let controls = secp_fold_controls(e, d, h, xed, eord, n10, hi_delta, hi_c);
    let kctrl = |i: usize| controls.get(i).copied().flatten();
    let maj2 = perpos_maj2_enabled();
    let park_low = core::cmp::min(configured_park_low, hi_delta);
    let host_all_derived = fold_host_derived_controls_enabled() && park_low >= 15;
    let host_h_xed_n10 =
        (fold_host_h_xed_n10_enabled() || host_all_derived) && park_low >= 14;
    let host_h_n10 =
        (fold_host_h_n10_enabled() || host_h_xed_n10) && park_low >= 13;
    let host_xed = host_h_xed_n10;
    let host_eord = host_all_derived;
    let host_e = fold_host_e_enabled() && host_all_derived && free_ed && park_low >= 17;
    let host_d = fold_host_d_enabled() && host_e && park_low >= 18;
    let host_n10 = (fold_host_n10_enabled() || host_h_n10) && park_low >= 12;
    let stream_controls =
        fold_stream_controls_enabled() && park_low >= 12 && !host_n10;

    if stream_controls {
        b.cx(h, n10);
        b.cx(d, n10);
        b.cx(h, eord);
        b.cx(xed, eord);
        b.cx(d, xed);
        b.cx(e, xed);
        b.free(n10);
        b.free(eord);
        b.free(xed);
        let mh = b.alloc_bit();
        b.hmr(h, mh);
        b.cz_if(e, d, mh);
        b.free(h);
    }

    // Carry lane is split so the WIDE tail is allocated only AFTER the four
    // derived controls are freed (the peak instant). `low` = active-region
    // carries [0..=hi_delta]; `tail` = pure-propagation carries (hi_delta, last].
    // Index map: carry i -> low[i] (i<=hi_delta) else tail[i-hi_delta-1].
    let low = if host_h_n10 {
        b.cx(h, n10);
        b.cx(d, n10);
        b.free(n10);
        if host_eord {
            b.cx(h, eord);
            b.cx(xed, eord);
            b.free(eord);
        }
        if host_xed {
            b.cx(d, xed);
            b.cx(e, xed);
            b.free(xed);
        }
        b.ccx(e, d, h);
        b.free(h);
        if host_e {
            let (ovf1, ovf2, _) = ed.expect("host_e requires overflow controls");
            b.cx(ovf1, e);
            b.cx(d, e);
            b.cx(ovf2, e);
            b.free(e);
        }
        // When host_d is enabled, the live d qubit itself is carry slot 29.
        // It remains d through control bits 1 and 5, then is coherently cleared
        // immediately before carry 29 is generated into the same physical slot.
        let d_slot = host_d.then_some(d);
        let e_slot = if host_e {
            let slot = b.alloc_qubit();
            debug_assert_eq!(slot, e);
            Some(slot)
        } else {
            None
        };
        let h_slot = b.alloc_qubit();
        debug_assert_eq!(h_slot, h);
        let xed_slot = if host_xed {
            let slot = b.alloc_qubit();
            debug_assert_eq!(slot, xed);
            Some(slot)
        } else {
            None
        };
        let eord_slot = if host_eord {
            let slot = b.alloc_qubit();
            debug_assert_eq!(slot, eord);
            Some(slot)
        } else {
            None
        };
        let n10_slot = b.alloc_qubit();
        debug_assert_eq!(n10_slot, n10);
        let regular = b.alloc_qubits(
            hi_delta
                - 1
                - usize::from(host_xed)
                - usize::from(host_eord)
                - usize::from(host_e)
                - usize::from(host_d),
        );
        let mut low = Vec::with_capacity(hi_delta + 1);
        let mut next_regular = 0usize;
        for i in 0..=hi_delta {
            if i == 28 && host_d {
                low.push(n10_slot);
            } else if i == 29 && host_d {
                low.push(d_slot.expect("hosted d slot"));
            } else if i == 29 && host_e {
                low.push(n10_slot);
            } else if i == 30 {
                low.push(h_slot);
            } else if i == 31 && host_xed {
                low.push(xed_slot.expect("hosted xed slot"));
            } else if i == 32 && host_eord {
                low.push(eord_slot.expect("hosted eord slot"));
            } else if i == hi_delta {
                low.push(e_slot.unwrap_or(n10_slot));
            } else {
                low.push(regular[next_regular]);
                next_regular += 1;
            }
        }
        debug_assert_eq!(next_regular, regular.len());
        low
    } else if host_n10 {
        b.cx(h, n10);
        b.cx(d, n10);
        b.free(n10);
        let n10_slot = b.alloc_qubit();
        debug_assert_eq!(n10_slot, n10);
        let mut low = b.alloc_qubits(hi_delta);
        low.push(n10_slot);
        low
    } else {
        b.alloc_qubits(hi_delta + 1)
    };

    // ── 1. active region [0..=hi_delta]: parked carries (controls live) ──
    let mut tail_d = None;
    let mut streamed_forward = None;
    if host_h_n10 {
        // h and n10 are needed only at bits 11 and 10. Host them in future-zero
        // carry slots 30 and 33, then clear both before carry generation reaches
        // slot 30. Their original IDs are restored only after carries 30..33
        // have been uncomputed.
        let e_host = host_e.then_some(low[hi_delta]);
        let h_host = low[30];
        let xed_host = host_xed.then_some(low[31]);
        let eord_host = host_eord.then_some(low[32]);
        let d_host = host_d.then_some(low[29]);
        let n10_host = if host_d {
            low[28]
        } else if host_e {
            low[29]
        } else {
            low[hi_delta]
        };
        let e_ctrl = e_host.unwrap_or(e);
        let d_ctrl = d_host.unwrap_or(d);
        if let Some(e_host) = e_host {
            let (ovf1, ovf2, _) = ed.expect("host_e requires overflow controls");
            b.cx(ovf1, e_host);
            b.cx(d_ctrl, e_host);
            b.cx(ovf2, e_host);
        }
        b.ccx(e_ctrl, d_ctrl, h_host);
        if let Some(xed_host) = xed_host {
            b.cx(e_ctrl, xed_host);
            b.cx(d_ctrl, xed_host);
        }
        if let Some(eord_host) = eord_host {
            b.cx(xed_host.expect("hosted xed for eord"), eord_host);
            b.cx(h_host, eord_host);
        }
        b.cx(d_ctrl, n10_host);
        b.cx(h_host, n10_host);
        for i in 0..=hi_delta {
            if (host_d && i == 28)
                || (!host_d && host_e && i == 29)
                || (!host_e && i == 30)
            {
                b.cx(h_host, n10_host);
                b.cx(d_ctrl, n10_host);
                if let Some(eord_host) = eord_host {
                    b.cx(h_host, eord_host);
                    b.cx(xed_host.expect("hosted xed for eord"), eord_host);
                }
                if let Some(xed_host) = xed_host {
                    b.cx(d_ctrl, xed_host);
                    b.cx(e_ctrl, xed_host);
                }
                b.ccx(e_ctrl, d_ctrl, h_host);
            }
            if host_d && i == 29 {
                let (ovf1, _, s2) = ed.expect("host_d requires overflow controls");
                b.ccx(ovf1, s2, d_ctrl);
            }
            if host_e && i == hi_delta {
                let (ovf1, ovf2, _) = ed.expect("host_e requires overflow controls");
                b.cx(ovf2, e_ctrl);
                if host_d {
                    b.cx(tail_d.expect("tail d is live at bit 33"), e_ctrl);
                } else {
                    b.cx(d_ctrl, e_ctrl);
                }
                b.cx(ovf1, e_ctrl);
            }
            let target = low[i];
            let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
            let kc = match i {
                0 | 4 | 6 | 32 if host_e => Some(e_ctrl),
                1 | 5 if host_d => d_host,
                7 if host_xed => xed_host,
                8 | 9 if host_eord => eord_host,
                10 => Some(n10_host),
                11 => Some(h_host),
                33 if host_d => tail_d,
                _ => kctrl(i),
            };
            if is_add {
                if let Some(kc) = kc {
                    if let Some(ci) = carry_in {
                        emit_fold_majority(b, acc[i], kc, ci, target, maj2);
                    } else {
                        b.ccx(acc[i], kc, target);
                    }
                } else if let Some(ci) = carry_in {
                    b.ccx(acc[i], ci, target);
                }
            } else if let Some(kc) = kc {
                b.x(acc[i]);
                if let Some(ci) = carry_in {
                    emit_fold_majority(b, acc[i], kc, ci, target, maj2);
                } else {
                    b.ccx(acc[i], kc, target);
                }
                b.x(acc[i]);
            } else if let Some(ci) = carry_in {
                b.x(acc[i]);
                b.ccx(acc[i], ci, target);
                b.x(acc[i]);
            }
            if let Some(kc) = kc {
                b.cx(kc, acc[i]);
            }
            if i > 0 {
                b.cx(low[i - 1], acc[i]);
            }
            if host_d && i == hi_delta - 1 {
                // Carry 31 is dead after carry/sum bit 32. Park it early and
                // reuse its physical slot as d for bit 33 and the wide tail.
                // Its carry value is reconstructed transiently during cleanup.
                let measured = b.alloc_bit();
                b.hmr(low[31], measured);
                fold_postsum_carry_phase_uncompute(
                    b,
                    acc,
                    None,
                    Some(low[30]),
                    measured,
                    31,
                    is_add,
                );
                b.free(low[31]);
                let slot = b.alloc_qubit();
                debug_assert_eq!(slot, low[31]);
                let (ovf1, _, s2) = ed.expect("host_d requires overflow controls");
                b.ccx(ovf1, s2, slot);
                tail_d = Some(slot);
            }
        }
    } else if host_n10 {
        // n10 is needed only at bit 10. Host it on low[33], which remains |0>
        // until the carry sweep reaches that position, then clear the host and
        // continue the same ripple. The original n10 ID is restored after the
        // parked low carries have been released.
        let n10_host = low[hi_delta];
        b.cx(d, n10_host);
        b.cx(h, n10_host);
        for i in 0..=hi_delta {
            if i == hi_delta {
                b.cx(h, n10_host);
                b.cx(d, n10_host);
            }
            let target = low[i];
            let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
            let kc = if i == 10 { Some(n10_host) } else { kctrl(i) };
            if is_add {
                if let Some(kc) = kc {
                    if let Some(ci) = carry_in {
                        emit_fold_majority(b, acc[i], kc, ci, target, maj2);
                    } else {
                        b.ccx(acc[i], kc, target);
                    }
                } else if let Some(ci) = carry_in {
                    b.ccx(acc[i], ci, target);
                }
            } else if let Some(kc) = kc {
                b.x(acc[i]);
                if let Some(ci) = carry_in {
                    emit_fold_majority(b, acc[i], kc, ci, target, maj2);
                } else {
                    b.ccx(acc[i], kc, target);
                }
                b.x(acc[i]);
            } else if let Some(ci) = carry_in {
                b.x(acc[i]);
                b.ccx(acc[i], ci, target);
                b.x(acc[i]);
            }
            if let Some(kc) = kc {
                b.cx(kc, acc[i]);
            }
            if i > 0 {
                b.cx(low[i - 1], acc[i]);
            }
        }
    } else if stream_controls {
        for i in 0..7 {
            let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
            fold_presum_carry_compute_and_sum(
                b,
                acc,
                kctrl(i),
                carry_in,
                low[i],
                i,
                is_add,
                maj2,
            );
        }
        let streamed = b.alloc_qubit();
        b.cx(e, streamed);
        b.cx(d, streamed);
        fold_presum_carry_compute_and_sum(
            b,
            acc,
            Some(streamed),
            Some(low[6]),
            low[7],
            7,
            is_add,
            maj2,
        );
        b.ccx(e, d, streamed);
        for i in 8..10 {
            fold_presum_carry_compute_and_sum(
                b,
                acc,
                Some(streamed),
                Some(low[i - 1]),
                low[i],
                i,
                is_add,
                maj2,
            );
        }
        b.cx(e, streamed);
        fold_presum_carry_compute_and_sum(
            b,
            acc,
            Some(streamed),
            Some(low[9]),
            low[10],
            10,
            is_add,
            maj2,
        );
        b.cx(d, streamed);
        fold_presum_carry_compute_and_sum(
            b,
            acc,
            Some(streamed),
            Some(low[10]),
            low[11],
            11,
            is_add,
            maj2,
        );
        for i in 12..=hi_delta {
            fold_presum_carry_compute_and_sum(
                b,
                acc,
                kctrl(i),
                Some(low[i - 1]),
                low[i],
                i,
                is_add,
                maj2,
            );
        }
        streamed_forward = Some(streamed);
    } else {
        for i in 0..=hi_delta {
            let target = low[i];
            let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
            if is_add {
                if let Some(kc) = kctrl(i) {
                    if let Some(ci) = carry_in {
                        emit_fold_majority(b, acc[i], kc, ci, target, maj2);
                    } else {
                        b.ccx(acc[i], kc, target);
                    }
                } else if let Some(ci) = carry_in {
                    b.ccx(acc[i], ci, target);
                }
            } else if let Some(kc) = kctrl(i) {
                b.x(acc[i]);
                if let Some(ci) = carry_in {
                    emit_fold_majority(b, acc[i], kc, ci, target, maj2);
                } else {
                    b.ccx(acc[i], kc, target);
                }
                b.x(acc[i]);
            } else if let Some(ci) = carry_in {
                b.x(acc[i]);
                b.ccx(acc[i], ci, target);
                b.x(acc[i]);
            }
        }
        // ── 2. low sum bits [0..=hi_delta] while the controls are still live ──
        //   acc_i ^= k_i ^ carry_{i-1}. (The seam sum at hi_delta+1 is k=0 and is
        //   written in step 4b AFTER the tail carries are generated from original acc.)
        for i in 0..=hi_delta {
            if let Some(kc) = kctrl(i) {
                b.cx(kc, acc[i]);
            }
            if i > 0 {
                b.cx(low[i - 1], acc[i]);
            }
        }
    }

    // Optional peak lever: the tail only needs low[hi_delta] as its carry-in.
    // The lower parked carries are needed again later for low-carry cleanup, so
    // measurement-uncompute them now and recompute from the post-sum bits after
    // the tail has been freed. Parking just one carry is enough to drop the
    // fused double/halve high-water by one qubit.
    if park_low > 0 {
        if host_h_n10 {
            for i in (12..park_low).rev() {
                let m = b.alloc_bit();
                b.hmr(low[i], m);
                fold_postsum_carry_phase_uncompute(
                    b,
                    acc,
                    kctrl(i),
                    Some(low[i - 1]),
                    m,
                    i,
                    is_add,
                );
                b.free(low[i]);
            }

            let d_ctrl = if host_d {
                tail_d.expect("tail d is live during parked-carry cleanup")
            } else {
                d
            };
            let rev_e = if host_e {
                let (ovf1, ovf2, _) = ed.expect("host_e requires overflow controls");
                let rev_e = b.alloc_qubit();
                b.cx(ovf1, rev_e);
                b.cx(d_ctrl, rev_e);
                b.cx(ovf2, rev_e);
                Some(rev_e)
            } else {
                None
            };
            let e_ctrl = rev_e.unwrap_or(e);
            let rev_h = b.alloc_qubit();
            b.ccx(e_ctrl, d_ctrl, rev_h);
            let measured_h = b.alloc_bit();
            b.hmr(low[11], measured_h);
            fold_postsum_carry_phase_uncompute(
                b,
                acc,
                Some(rev_h),
                Some(low[10]),
                measured_h,
                11,
                is_add,
            );
            b.free(low[11]);

            let rev_n10 = b.alloc_qubit();
            debug_assert_eq!(rev_n10, low[11]);
            b.cx(d_ctrl, rev_n10);
            b.cx(rev_h, rev_n10);
            let measured_n10 = b.alloc_bit();
            b.hmr(low[10], measured_n10);
            fold_postsum_carry_phase_uncompute(
                b,
                acc,
                Some(rev_n10),
                Some(low[9]),
                measured_n10,
                10,
                is_add,
            );
            b.free(low[10]);
            b.cx(rev_h, rev_n10);
            b.cx(d_ctrl, rev_n10);
            b.free(rev_n10);

            let rev_xed = if host_xed {
                let rev_xed = b.alloc_qubit();
                b.cx(e_ctrl, rev_xed);
                b.cx(d_ctrl, rev_xed);
                Some(rev_xed)
            } else {
                None
            };
            let rev_eord = if host_eord {
                let rev_eord = b.alloc_qubit();
                b.cx(rev_xed.expect("hosted xed for eord"), rev_eord);
                b.cx(rev_h, rev_eord);
                Some(rev_eord)
            } else {
                None
            };
            if !host_xed {
                b.ccx(e, d, rev_h);
                b.free(rev_h);
            }
            for i in (0..10).rev() {
                let m = b.alloc_bit();
                b.hmr(low[i], m);
                let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
                let kc = match i {
                    0 | 4 | 6 if host_e => rev_e,
                    1 | 5 if host_d => Some(d_ctrl),
                    7 if host_xed => rev_xed,
                    8 | 9 if host_eord => rev_eord,
                    _ => kctrl(i),
                };
                fold_postsum_carry_phase_uncompute(
                    b,
                    acc,
                    kc,
                    carry_in,
                    m,
                    i,
                    is_add,
                );
                b.free(low[i]);
            }
            if let Some(rev_eord) = rev_eord {
                b.cx(rev_h, rev_eord);
                b.cx(rev_xed.expect("hosted xed for eord"), rev_eord);
                b.free(rev_eord);
            }
            if let Some(rev_xed) = rev_xed {
                b.cx(d_ctrl, rev_xed);
                b.cx(e_ctrl, rev_xed);
                b.free(rev_xed);
            }
            if host_xed {
                b.ccx(e_ctrl, d_ctrl, rev_h);
                b.free(rev_h);
            }
            if let Some(rev_e) = rev_e {
                let (ovf1, ovf2, _) = ed.expect("host_e requires overflow controls");
                b.cx(ovf2, rev_e);
                b.cx(d_ctrl, rev_e);
                b.cx(ovf1, rev_e);
                b.free(rev_e);
            }
        } else if host_n10 {
            for i in (11..park_low).rev() {
                let m = b.alloc_bit();
                b.hmr(low[i], m);
                fold_postsum_carry_phase_uncompute(
                    b,
                    acc,
                    kctrl(i),
                    Some(low[i - 1]),
                    m,
                    i,
                    is_add,
                );
                b.free(low[i]);
            }
            let rev_n10 = b.alloc_qubit();
            b.cx(d, rev_n10);
            b.cx(h, rev_n10);
            for i in (0..11).rev() {
                let m = b.alloc_bit();
                b.hmr(low[i], m);
                let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
                let kc = if i == 10 { Some(rev_n10) } else { kctrl(i) };
                fold_postsum_carry_phase_uncompute(
                    b, acc, kc, carry_in, m, i, is_add,
                );
                b.free(low[i]);
            }
            b.cx(h, rev_n10);
            b.cx(d, rev_n10);
            b.free(rev_n10);
        } else if stream_controls {
            let streamed = streamed_forward
                .take()
                .expect("streamed forward control is live");
            for i in (12..park_low).rev() {
                let m = b.alloc_bit();
                b.hmr(low[i], m);
                fold_postsum_carry_phase_uncompute(
                    b,
                    acc,
                    kctrl(i),
                    Some(low[i - 1]),
                    m,
                    i,
                    is_add,
                );
                b.free(low[i]);
            }
            let m11 = b.alloc_bit();
            b.hmr(low[11], m11);
            fold_postsum_carry_phase_uncompute(
                b,
                acc,
                Some(streamed),
                Some(low[10]),
                m11,
                11,
                is_add,
            );
            b.free(low[11]);
            b.cx(d, streamed);
            let m10 = b.alloc_bit();
            b.hmr(low[10], m10);
            fold_postsum_carry_phase_uncompute(
                b,
                acc,
                Some(streamed),
                Some(low[9]),
                m10,
                10,
                is_add,
            );
            b.free(low[10]);
            b.cx(e, streamed);
            for i in (8..10).rev() {
                let m = b.alloc_bit();
                b.hmr(low[i], m);
                fold_postsum_carry_phase_uncompute(
                    b,
                    acc,
                    Some(streamed),
                    Some(low[i - 1]),
                    m,
                    i,
                    is_add,
                );
                b.free(low[i]);
            }
            b.ccx(e, d, streamed);
            let m7 = b.alloc_bit();
            b.hmr(low[7], m7);
            fold_postsum_carry_phase_uncompute(
                b,
                acc,
                Some(streamed),
                Some(low[6]),
                m7,
                7,
                is_add,
            );
            b.free(low[7]);
            b.cx(d, streamed);
            b.cx(e, streamed);
            b.free(streamed);
            for i in (0..7).rev() {
                let m = b.alloc_bit();
                b.hmr(low[i], m);
                let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
                fold_postsum_carry_phase_uncompute(
                    b,
                    acc,
                    kctrl(i),
                    carry_in,
                    m,
                    i,
                    is_add,
                );
                b.free(low[i]);
            }
        } else {
            for i in (0..park_low).rev() {
                let m = b.alloc_bit();
                b.hmr(low[i], m);
                let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
                fold_postsum_carry_phase_uncompute(
                    b,
                    acc,
                    kctrl(i),
                    carry_in,
                    m,
                    i,
                    is_add,
                );
                b.free(low[i]);
            }
        }
    }
    debug_assert!(streamed_forward.is_none());

    // ── 3. free the four e,d-derived controls BEFORE allocating the wide tail ──
    //   (reset them to |0> via the free-CX uncompute + a measured AND-clear, then
    //   release the SAME qubits; reacquired+re-derived in step 6 so the caller
    //   sees them live & correct on return, exactly as the baseline ripple does.)
    if !stream_controls {
        if !host_n10 {
            b.cx(h, n10);
            b.cx(d, n10);
            b.free(n10);
        }
        if !host_eord {
            if host_h_n10 {
                let rev_h = b.alloc_qubit();
                b.ccx(e, d, rev_h);
                b.cx(rev_h, eord);
                b.ccx(e, d, rev_h);
                b.free(rev_h);
            } else {
                b.cx(h, eord);
            }
            if host_xed {
                let rev_xed = b.alloc_qubit();
                b.cx(e, rev_xed);
                b.cx(d, rev_xed);
                b.cx(rev_xed, eord);
                b.cx(d, rev_xed);
                b.cx(e, rev_xed);
                b.free(rev_xed);
            } else {
                b.cx(xed, eord);
            }
            b.free(eord);
        }
        if !host_xed {
            b.cx(d, xed);
            b.cx(e, xed);
            b.free(xed);
        }
        if !host_h_n10 {
            let mh = b.alloc_bit();
            b.hmr(h, mh);
            b.cz_if(e, d, mh);
            b.free(h);
        }
    }

    // ── 3b. e,d-extension: free e,d too (HYP-6 §4a) ──
    //   `e,d` are dead controls in the tail. Uncompute them to |0> (e first, since
    //   it is built from d), then release the SAME qubits; reacquired+re-derived
    //   in step 6 from the live overflow lanes (ovf1,ovf2,s2 are unchanged here).
    //   Uncompute mirrors the dispatch-site derivation: d = ovf1&s2 (CCX),
    //   e = ovf1 ^ d ^ ovf2 (3 CX). Reversing: e via the same 3 CX (d still live),
    //   then d via a measured AND-clear (ovf1,s2 unchanged ⇒ d == ovf1&s2 ⇒
    //   hmr+cz_if forces d→0, 0 Toffoli, phase-exact).
    if free_ed {
        let (ovf1, ovf2, s2) = ed.expect("free_ed implies ed is Some");
        if !host_e {
            b.cx(ovf1, e);
            b.cx(d, e);
            b.cx(ovf2, e);
            b.free(e);
        }
        if !host_d {
            let md = b.alloc_bit();
            b.hmr(d, md);
            b.cz_if(ovf1, s2, md);
            b.free(d);
        }
    }

    // Tail carries are allocated NOW (4 derived controls already freed ⇒ the
    // wide-lane high-water carries +4 ancillae instead of +8; +2 with the
    // e,d-extension since e,d are freed as well).
    let tail_len = last - hi_delta;
    let tail = b.alloc_qubits(tail_len);
    let cw = |i: usize| -> QubitId {
        if i <= hi_delta {
            low[i]
        } else {
            tail[i - hi_delta - 1]
        }
    };

    // ── 4a. high-tail carry generation (hi_delta, last]: pure propagation from
    //   ORIGINAL acc (acc[hi_delta+1..] untouched by step 2) ──
    for i in (hi_delta + 1)..=last {
        if is_add {
            b.ccx(acc[i], cw(i - 1), cw(i));
        } else {
            b.x(acc[i]);
            b.ccx(acc[i], cw(i - 1), cw(i));
            b.x(acc[i]);
        }
    }
    // ── 4b. high sum bits (hi_delta, last+1] (k=0, control-free): acc_i ^= carry_{i-1} ──
    for i in (hi_delta + 1)..n {
        if i - 1 <= last {
            b.cx(cw(i - 1), acc[i]);
        }
    }

    // ── 5. reverse uncompute the TAIL carries first (control-free), freeing
    //   them high→low so the wide lane shrinks before the derived controls return ──
    for i in (hi_delta + 1..=last).rev() {
        let m = b.alloc_bit();
        b.hmr(cw(i), m);
        let carry_in = cw(i - 1);
        if is_add {
            b.x(acc[i]);
            b.cz_if(acc[i], carry_in, m);
            b.x(acc[i]);
        } else {
            b.cz_if(acc[i], carry_in, m);
        }
        b.free(cw(i));
    }
    drop(tail);

    // ── 6. recompute the controls INTO THE SAME qubits (reacquire + re-derive)
    //   for the low uncompute pass; they stay live on return. ──
    //   e,d-extension: re-derive e,d FIRST (the four derived controls depend on
    //   them), from the live overflow lanes — exactly the dispatch-site formula:
    //   d = ovf1 & s2, e = ovf1 ^ d ^ ovf2.
    if free_ed && !host_d {
        let (ovf1, ovf2, s2) = ed.expect("free_ed implies ed is Some");
        b.reacquire(d);
        b.ccx(ovf1, s2, d);
        if !host_e {
            b.reacquire(e);
            b.cx(ovf1, e);
            b.cx(d, e);
            b.cx(ovf2, e);
        }
    }
    if host_h_n10 {
        if host_e {
            let measured = b.alloc_bit();
            b.hmr(low[hi_delta], measured);
            if host_d {
                fold_postsum_carry_phase_uncompute(
                    b,
                    acc,
                    tail_d,
                    Some(low[hi_delta - 1]),
                    measured,
                    hi_delta,
                    is_add,
                );
            } else {
                fold_postsum_carry_phase_uncompute(
                    b,
                    acc,
                    Some(d),
                    Some(low[hi_delta - 1]),
                    measured,
                    hi_delta,
                    is_add,
                );
            }
            b.free(low[hi_delta]);
            let (ovf1, ovf2, _) = ed.expect("host_e requires overflow controls");
            b.reacquire(e);
            b.cx(ovf1, e);
            if host_d {
                b.cx(tail_d.expect("tail d is live while restoring e"), e);
            } else {
                b.cx(d, e);
            }
            b.cx(ovf2, e);
        }
        if host_d {
            // low[31] carries d across the tail. Reconstruct carry 31 into a
            // temporary clean slot solely to phase-uncompute carry 32.
            let carry31 = b.alloc_qubit();
            fold_postsum_carry_compute(
                b,
                acc,
                None,
                Some(low[30]),
                carry31,
                31,
                is_add,
            );
            let measured32 = b.alloc_bit();
            b.hmr(low[32], measured32);
            fold_postsum_carry_phase_uncompute(
                b,
                acc,
                Some(e),
                Some(carry31),
                measured32,
                32,
                is_add,
            );
            b.free(low[32]);
            let measured31 = b.alloc_bit();
            b.hmr(carry31, measured31);
            fold_postsum_carry_phase_uncompute(
                b,
                acc,
                None,
                Some(low[30]),
                measured31,
                31,
                is_add,
            );
            b.free(carry31);

            for i in (28..=30).rev() {
                let measured = b.alloc_bit();
                b.hmr(low[i], measured);
                fold_postsum_carry_phase_uncompute(
                    b,
                    acc,
                    None,
                    Some(low[i - 1]),
                    measured,
                    i,
                    is_add,
                );
                b.free(low[i]);
            }

            // Move d from the borrowed carry-31 slot back to its original
            // carry-29 qubit using only Clifford gates.
            let d_tail = tail_d.expect("tail d is live during d restoration");
            b.reacquire(d);
            b.cx(d_tail, d);
            b.cx(d, d_tail);
            b.free(d_tail);
        } else {
            let high_start = if host_e { 29 } else { 30 };
            let high_end = if host_e { hi_delta - 1 } else { hi_delta };
            for i in (high_start..=high_end).rev() {
                let measured = b.alloc_bit();
                b.hmr(low[i], measured);
                let kc = match i {
                    32 => Some(e),
                    33 => Some(d),
                    _ => None,
                };
                fold_postsum_carry_phase_uncompute(
                    b,
                    acc,
                    kc,
                    Some(low[i - 1]),
                    measured,
                    i,
                    is_add,
                );
                b.free(low[i]);
            }
        }
    } else if host_n10 {
        let measured = b.alloc_bit();
        b.hmr(low[hi_delta], measured);
        fold_postsum_carry_phase_uncompute(
            b,
            acc,
            Some(d),
            Some(low[hi_delta - 1]),
            measured,
            hi_delta,
            is_add,
        );
        b.free(low[hi_delta]);
    }
    if stream_controls {
        for i in 0..park_low {
            b.reacquire(low[i]);
        }

        for i in 0..7 {
            let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
            fold_postsum_carry_compute(
                b,
                acc,
                kctrl(i),
                carry_in,
                low[i],
                i,
                is_add,
            );
        }

        // One temporary control walks through all four nonlinear predicates:
        // xed=e^d, eord=e|d=xed^(e&d), n10=eord^e, h=n10^d.
        let streamed = b.alloc_qubit();
        b.cx(e, streamed);
        b.cx(d, streamed);
        fold_postsum_carry_compute(
            b,
            acc,
            Some(streamed),
            Some(low[6]),
            low[7],
            7,
            is_add,
        );
        b.ccx(e, d, streamed);
        for i in 8..10 {
            fold_postsum_carry_compute(
                b,
                acc,
                Some(streamed),
                Some(low[i - 1]),
                low[i],
                i,
                is_add,
            );
        }
        b.cx(e, streamed);
        fold_postsum_carry_compute(
            b,
            acc,
            Some(streamed),
            Some(low[9]),
            low[10],
            10,
            is_add,
        );
        b.cx(d, streamed);
        fold_postsum_carry_compute(
            b,
            acc,
            Some(streamed),
            Some(low[10]),
            low[11],
            11,
            is_add,
        );
        for i in 12..park_low {
            fold_postsum_carry_compute(
                b,
                acc,
                kctrl(i),
                Some(low[i - 1]),
                low[i],
                i,
                is_add,
            );
        }

        // High active carries use only direct e/d controls. Keeping `streamed`
        // as h across this section avoids a second h derivation.
        for i in (12..=hi_delta).rev() {
            let m = b.alloc_bit();
            b.hmr(low[i], m);
            fold_postsum_carry_phase_uncompute(
                b,
                acc,
                kctrl(i),
                Some(low[i - 1]),
                m,
                i,
                is_add,
            );
            b.free(low[i]);
        }

        let m11 = b.alloc_bit();
        b.hmr(low[11], m11);
        fold_postsum_carry_phase_uncompute(
            b,
            acc,
            Some(streamed),
            Some(low[10]),
            m11,
            11,
            is_add,
        );
        b.free(low[11]);
        b.cx(d, streamed);

        let m10 = b.alloc_bit();
        b.hmr(low[10], m10);
        fold_postsum_carry_phase_uncompute(
            b,
            acc,
            Some(streamed),
            Some(low[9]),
            m10,
            10,
            is_add,
        );
        b.free(low[10]);
        b.cx(e, streamed);

        for i in (8..10).rev() {
            let m = b.alloc_bit();
            b.hmr(low[i], m);
            fold_postsum_carry_phase_uncompute(
                b,
                acc,
                Some(streamed),
                Some(low[i - 1]),
                m,
                i,
                is_add,
            );
            b.free(low[i]);
        }
        b.ccx(e, d, streamed);

        let m7 = b.alloc_bit();
        b.hmr(low[7], m7);
        fold_postsum_carry_phase_uncompute(
            b,
            acc,
            Some(streamed),
            Some(low[6]),
            m7,
            7,
            is_add,
        );
        b.free(low[7]);
        b.cx(d, streamed);
        b.cx(e, streamed);
        b.free(streamed);

        for i in (0..7).rev() {
            let m = b.alloc_bit();
            b.hmr(low[i], m);
            let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
            fold_postsum_carry_phase_uncompute(
                b,
                acc,
                kctrl(i),
                carry_in,
                m,
                i,
                is_add,
            );
            b.free(low[i]);
        }
        drop(low);

        // Restore the caller-visible controls only after the carry lane is gone.
        b.reacquire(h);
        b.ccx(e, d, h);
        b.reacquire(xed);
        b.cx(e, xed);
        b.cx(d, xed);
        b.reacquire(eord);
        b.cx(xed, eord);
        b.cx(h, eord);
        b.reacquire(n10);
        b.cx(d, n10);
        b.cx(h, n10);
    } else {
        b.reacquire(h);
        b.ccx(e, d, h);
        b.reacquire(xed);
        b.cx(e, xed);
        b.cx(d, xed);
        b.reacquire(eord);
        b.cx(xed, eord);
        b.cx(h, eord);
        b.reacquire(n10);
        b.cx(d, n10);
        b.cx(h, n10);

        if park_low > 0 {
            for i in 0..park_low {
                b.reacquire(low[i]);
                let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
                fold_postsum_carry_compute(
                    b,
                    acc,
                    kctrl(i),
                    carry_in,
                    low[i],
                    i,
                    is_add,
                );
            }
        }

        // ── 7. reverse uncompute the active carries [0..=hi_delta] ──
        let low_top = if host_d {
            27
        } else if host_e {
            28
        } else if host_h_n10 {
            29
        } else if host_n10 {
            hi_delta - 1
        } else {
            hi_delta
        };
        for i in (0..=low_top).rev() {
            let m = b.alloc_bit();
            b.hmr(low[i], m);
            let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
            fold_postsum_carry_phase_uncompute(
                b,
                acc,
                kctrl(i),
                carry_in,
                m,
                i,
                is_add,
            );
            b.free(low[i]);
        }
        drop(low);
    }
    // h, xed, eord, n10 are left LIVE and value-correct (= baseline post-ripple
    // state); the caller's normal derived-control uncompute block runs next.
}

fn fold_ripple_freed_tail_ed_hosted(
    b: &mut B,
    acc: &[QubitId],
    e: QubitId,
    d: QubitId,
    h: QubitId,
    xed: QubitId,
    eord: QubitId,
    n10: QubitId,
    ed: Option<(QubitId, QubitId, QubitId)>,
    park_low: usize,
    last: usize,
    is_add: bool,
) {
    let free_ed = ed.is_some() && fold_freed_tail_ed_enabled();
    let n = acc.len();
    let hi_delta = 33usize;
    let split = 11usize;
    debug_assert!(last < n);
    debug_assert!(last > hi_delta, "freed-tail requires a nonempty high tail");
    let maj2 = perpos_maj2_enabled();
    let park_low = core::cmp::min(park_low, split + 1);

    // Release the four derived controls before allocating the low-carry lane.
    b.cx(h, n10);
    b.cx(d, n10);
    b.free(n10);
    b.cx(h, eord);
    b.cx(xed, eord);
    b.cx(d, xed);
    b.cx(e, xed);
    b.free(eord);
    b.free(xed);
    let mh = b.alloc_bit();
    b.hmr(h, mh);
    b.cz_if(e, d, mh);
    b.free(h);

    let low = b.alloc_qubits(hi_delta + 1);
    // These slots stay zero until carry generation reaches bits 30..33.
    let h_host = low[30];
    let xed_host = low[31];
    let eord_host = low[32];
    let n10_host = low[33];
    b.ccx(e, d, h_host);
    b.cx(e, xed_host);
    b.cx(d, xed_host);
    b.cx(xed_host, eord_host);
    b.cx(h_host, eord_host);
    b.cx(d, n10_host);
    b.cx(h_host, n10_host);

    let hosted_kctrl = |i: usize| -> Option<QubitId> {
        match i {
            0 | 4 | 6 | 32 => Some(e),
            1 | 5 | 33 => Some(d),
            7 => Some(xed_host),
            8 | 9 => Some(eord_host),
            10 => Some(n10_host),
            11 => Some(h_host),
            _ => None,
        }
    };

    // Compute the carries that depend on the sparse derived controls.
    for i in 0..=split {
        let target = low[i];
        let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
        if is_add {
            if let Some(kc) = hosted_kctrl(i) {
                if let Some(ci) = carry_in {
                    emit_fold_majority(b, acc[i], kc, ci, target, maj2);
                } else {
                    b.ccx(acc[i], kc, target);
                }
            } else if let Some(ci) = carry_in {
                b.ccx(acc[i], ci, target);
            }
        } else if let Some(kc) = hosted_kctrl(i) {
            b.x(acc[i]);
            if let Some(ci) = carry_in {
                emit_fold_majority(b, acc[i], kc, ci, target, maj2);
            } else {
                b.ccx(acc[i], kc, target);
            }
            b.x(acc[i]);
        } else if let Some(ci) = carry_in {
            b.x(acc[i]);
            b.ccx(acc[i], ci, target);
            b.x(acc[i]);
        }
    }
    for i in 0..=split {
        if let Some(kc) = hosted_kctrl(i) {
            b.cx(kc, acc[i]);
        }
        if i > 0 {
            b.cx(low[i - 1], acc[i]);
        }
    }

    // Return the hosted controls to zero before their slots become carries.
    b.cx(h_host, n10_host);
    b.cx(d, n10_host);
    b.cx(h_host, eord_host);
    b.cx(xed_host, eord_host);
    b.cx(d, xed_host);
    b.cx(e, xed_host);
    let mh_host = b.alloc_bit();
    b.hmr(h_host, mh_host);
    b.cz_if(e, d, mh_host);

    // Continue through the control-free middle and the direct e/d high bits.
    for i in split + 1..=hi_delta {
        let target = low[i];
        let carry_in = Some(low[i - 1]);
        let kctrl = match i {
            32 => Some(e),
            33 => Some(d),
            _ => None,
        };
        if is_add {
            if let Some(kc) = kctrl {
                emit_fold_majority(
                    b,
                    acc[i],
                    kc,
                    carry_in.expect("high carry-in"),
                    target,
                    maj2,
                );
            } else {
                b.ccx(
                    acc[i],
                    carry_in.expect("high carry-in"),
                    target,
                );
            }
        } else if let Some(kc) = kctrl {
            b.x(acc[i]);
            emit_fold_majority(
                b,
                acc[i],
                kc,
                carry_in.expect("high carry-in"),
                target,
                maj2,
            );
            b.x(acc[i]);
        } else {
            b.x(acc[i]);
            b.ccx(
                acc[i],
                carry_in.expect("high carry-in"),
                target,
            );
            b.x(acc[i]);
        }
    }
    for i in split + 1..=hi_delta {
        let kctrl = match i {
            32 => Some(e),
            33 => Some(d),
            _ => None,
        };
        if let Some(kc) = kctrl {
            b.cx(kc, acc[i]);
        }
        b.cx(low[i - 1], acc[i]);
    }

    if park_low > 0 {
        let controls = secp_fold_controls(
            e, d, h_host, xed_host, eord_host, n10_host, hi_delta, 32,
        );
        for i in (0..park_low).rev() {
            let m = b.alloc_bit();
            b.hmr(low[i], m);
            let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
            fold_postsum_carry_phase_uncompute(
                b,
                acc,
                controls.get(i).copied().flatten(),
                carry_in,
                m,
                i,
                is_add,
            );
            b.free(low[i]);
        }
    }

    if free_ed {
        let (ovf1, ovf2, s2) = ed.expect("free_ed implies ed is Some");
        b.cx(ovf1, e);
        b.cx(d, e);
        b.cx(ovf2, e);
        b.free(e);
        let md = b.alloc_bit();
        b.hmr(d, md);
        b.cz_if(ovf1, s2, md);
        b.free(d);
    }

    let tail_len = last - hi_delta;
    let tail = b.alloc_qubits(tail_len);
    let cw = |i: usize| -> QubitId {
        if i <= hi_delta {
            low[i]
        } else {
            tail[i - hi_delta - 1]
        }
    };

    for i in (hi_delta + 1)..=last {
        if is_add {
            b.ccx(acc[i], cw(i - 1), cw(i));
        } else {
            b.x(acc[i]);
            b.ccx(acc[i], cw(i - 1), cw(i));
            b.x(acc[i]);
        }
    }
    for i in (hi_delta + 1)..n {
        if i - 1 <= last {
            b.cx(cw(i - 1), acc[i]);
        }
    }
    for i in (hi_delta + 1..=last).rev() {
        let m = b.alloc_bit();
        b.hmr(cw(i), m);
        let carry_in = cw(i - 1);
        if is_add {
            b.x(acc[i]);
            b.cz_if(acc[i], carry_in, m);
            b.x(acc[i]);
        } else {
            b.cz_if(acc[i], carry_in, m);
        }
        b.free(cw(i));
    }
    drop(tail);

    if free_ed {
        let (ovf1, ovf2, s2) = ed.expect("free_ed implies ed is Some");
        b.reacquire(d);
        b.ccx(ovf1, s2, d);
        b.reacquire(e);
        b.cx(ovf1, e);
        b.cx(d, e);
        b.cx(ovf2, e);
    }

    // Bits 12..33 need only direct e/d controls. Release them before restoring
    // the four derived-control qubits used by bits 0..11.
    for i in (split + 1..=hi_delta).rev() {
        let m = b.alloc_bit();
        b.hmr(low[i], m);
        let kctrl = match i {
            32 => Some(e),
            33 => Some(d),
            _ => None,
        };
        fold_postsum_carry_phase_uncompute(
            b,
            acc,
            kctrl,
            Some(low[i - 1]),
            m,
            i,
            is_add,
        );
        b.free(low[i]);
    }

    // The original control IDs were likely reused by low[0..3], so they cannot
    // be reacquired yet. Hold the reverse-pass controls in already-freed high
    // carry slots, then transfer them back after low[0..11] is released.
    let rev_h = b.alloc_qubit();
    let rev_xed = b.alloc_qubit();
    let rev_eord = b.alloc_qubit();
    let rev_n10 = b.alloc_qubit();
    b.ccx(e, d, rev_h);
    b.cx(e, rev_xed);
    b.cx(d, rev_xed);
    b.cx(rev_xed, rev_eord);
    b.cx(rev_h, rev_eord);
    b.cx(d, rev_n10);
    b.cx(rev_h, rev_n10);
    let controls =
        secp_fold_controls(e, d, rev_h, rev_xed, rev_eord, rev_n10, hi_delta, 32);

    if park_low > 0 {
        for i in 0..park_low {
            b.reacquire(low[i]);
            let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
            fold_postsum_carry_compute(
                b,
                acc,
                controls.get(i).copied().flatten(),
                carry_in,
                low[i],
                i,
                is_add,
            );
        }
    }

    for i in (0..=split).rev() {
        let m = b.alloc_bit();
        b.hmr(low[i], m);
        let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
        fold_postsum_carry_phase_uncompute(
            b,
            acc,
            controls.get(i).copied().flatten(),
            carry_in,
            m,
            i,
            is_add,
        );
        b.free(low[i]);
    }
    drop(low);

    b.reacquire(h);
    b.cx(rev_h, h);
    b.reacquire(xed);
    b.cx(rev_xed, xed);
    b.reacquire(eord);
    b.cx(rev_eord, eord);
    b.reacquire(n10);
    b.cx(rev_n10, n10);

    b.cx(rev_h, rev_n10);
    b.cx(d, rev_n10);
    b.cx(rev_h, rev_eord);
    b.cx(rev_xed, rev_eord);
    b.cx(d, rev_xed);
    b.cx(e, rev_xed);
    b.free(rev_n10);
    b.free(rev_eord);
    b.free(rev_xed);
    let mh_rev = b.alloc_bit();
    b.hmr(rev_h, mh_rev);
    b.cz_if(e, d, mh_rev);
    b.free(rev_h);
}
