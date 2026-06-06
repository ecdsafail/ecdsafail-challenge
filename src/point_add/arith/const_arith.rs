use super::*;

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
                b.ccx(acc[i], bi, target);
                b.ccx(ctrl, acc[i], target);
                b.ccx(ctrl, bi, target);
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
                b.ccx(acc[i], ci, target);
                b.ccx(ctrl, acc[i], target);
                b.ccx(ctrl, ci, target);
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
    let carries = b.alloc_qubits(last + 1);

    // Forward carry sweep, truncated at `last`. carry_{i+1} = maj(acc_i, k_i, carry_i).
    for i in 0..=last {
        let target = carries[i];
        let carry_in = if i == 0 { None } else { Some(carries[i - 1]) };
        if bit(c, i) {
            if let Some(ci) = carry_in {
                b.ccx(acc[i], ci, target);
                b.ccx(ctrl, acc[i], target);
                b.ccx(ctrl, ci, target);
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

    b.free_vec(&carries);
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
    let borrows = b.alloc_qubits(last + 1);

    // Forward borrow sweep, truncated at `last`.
    for i in 0..=last {
        let target = borrows[i];
        let borrow_in = if i == 0 { None } else { Some(borrows[i - 1]) };
        if bit(c, i) {
            b.x(acc[i]);
            if let Some(bi) = borrow_in {
                b.ccx(acc[i], bi, target);
                b.ccx(ctrl, acc[i], target);
                b.ccx(ctrl, bi, target);
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

    b.free_vec(&borrows);
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
    let carries = b.alloc_qubits(last + 1);

    // Forward carry sweep, truncated at `last`. carry_i = maj(acc_i, k_i, carry_{i-1}).
    for i in 0..=last {
        let target = carries[i];
        let carry_in = if i == 0 { None } else { Some(carries[i - 1]) };
        if let Some(kc) = kctrl(i) {
            if let Some(ci) = carry_in {
                b.ccx(acc[i], ci, target);
                b.ccx(kc, acc[i], target);
                b.ccx(kc, ci, target);
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
    let borrows = b.alloc_qubits(last + 1);

    // Forward borrow sweep, truncated at `last`.
    for i in 0..=last {
        let target = borrows[i];
        let borrow_in = if i == 0 { None } else { Some(borrows[i - 1]) };
        if let Some(kc) = kctrl(i) {
            b.x(acc[i]);
            if let Some(bi) = borrow_in {
                b.ccx(acc[i], bi, target);
                b.ccx(kc, acc[i], target);
                b.ccx(kc, bi, target);
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
