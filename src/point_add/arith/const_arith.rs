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
                b.ccx(acc[i], kc, target);
                b.ccx(acc[i], ci, target);
                b.x(acc[i]);
                b.ccx(kc, ci, target);
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
            b.ccx(acc[i], kc, target);
            b.ccx(acc[i], ci, target);
            b.ccx(kc, ci, target);
        } else {
            b.ccx(acc[i], kc, target);
        }
    } else if let Some(ci) = carry_in {
        b.ccx(acc[i], ci, target);
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
    fold_ripple_freed_tail_ed(b, acc, e, d, h, xed, eord, n10, None, last, is_add);
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
    last: usize,
    is_add: bool,
) {
    let free_ed = ed.is_some() && fold_freed_tail_ed_enabled();
    let n = acc.len();
    let hi_delta = 33usize; // highest_set_bit(2^32+977)+1
    let hi_c = 32usize;
    debug_assert!(last < n);
    debug_assert!(last > hi_delta, "freed-tail requires a nonempty high tail");
    let controls = secp_fold_controls(e, d, h, xed, eord, n10, hi_delta, hi_c);
    let kctrl = |i: usize| controls.get(i).copied().flatten();
    let maj2 = perpos_maj2_enabled();
    let park_low = core::cmp::min(fold_park_low_carries(), hi_delta);

    // Carry lane is split so the WIDE tail is allocated only AFTER the four
    // derived controls are freed (the peak instant). `low` = active-region
    // carries [0..=hi_delta]; `tail` = pure-propagation carries (hi_delta, last].
    // Index map: carry i -> low[i] (i<=hi_delta) else tail[i-hi_delta-1].
    let low = b.alloc_qubits(hi_delta + 1);

    // ── 1. active region [0..=hi_delta]: parked carries (controls live) ──
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

    // Optional peak lever: the tail only needs low[hi_delta] as its carry-in.
    // The lower parked carries are needed again later for low-carry cleanup, so
    // measurement-uncompute them now and recompute from the post-sum bits after
    // the tail has been freed. Parking just one carry is enough to drop the
    // fused double/halve high-water by one qubit.
    if park_low > 0 {
        for i in (0..park_low).rev() {
            let m = b.alloc_bit();
            b.hmr(low[i], m);
            let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
            fold_postsum_carry_phase_uncompute(b, acc, kctrl(i), carry_in, m, i, is_add);
            b.free(low[i]);
        }
    }

    // ── 3. free the four e,d-derived controls BEFORE allocating the wide tail ──
    //   (reset them to |0> via the free-CX uncompute + a measured AND-clear, then
    //   release the SAME qubits; reacquired+re-derived in step 6 so the caller
    //   sees them live & correct on return, exactly as the baseline ripple does.)
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
        b.cx(ovf1, e);
        b.cx(d, e);
        b.cx(ovf2, e);
        b.free(e);
        let md = b.alloc_bit();
        b.hmr(d, md);
        b.cz_if(ovf1, s2, md);
        b.free(d);
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
    if free_ed {
        let (ovf1, ovf2, s2) = ed.expect("free_ed implies ed is Some");
        b.reacquire(d);
        b.ccx(ovf1, s2, d);
        b.reacquire(e);
        b.cx(ovf1, e);
        b.cx(d, e);
        b.cx(ovf2, e);
    }
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
            fold_postsum_carry_compute(b, acc, kctrl(i), carry_in, low[i], i, is_add);
        }
    }

    // ── 7. reverse uncompute the active carries [0..=hi_delta] ──
    for i in (0..=hi_delta).rev() {
        let m = b.alloc_bit();
        b.hmr(low[i], m);
        let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
        fold_postsum_carry_phase_uncompute(b, acc, kctrl(i), carry_in, m, i, is_add);
        b.free(low[i]);
    }
    drop(low);
    // h, xed, eord, n10 are left LIVE and value-correct (= baseline post-ripple
    // state); the caller's normal derived-control uncompute block runs next.
}
