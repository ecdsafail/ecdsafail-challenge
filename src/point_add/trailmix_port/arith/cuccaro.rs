//! Cuccaro ripple-carry adders and subtractors (MAJ/UMA chains), including
//! the controlled, overflow-capturing, and 3n-Toffoli low-depth variants.
//! Extracted from the former `mbu_primitives` grab-bag.

use crate::point_add::trailmix_port::arith::mcx::{mcx_clean_k, mcx_dirty};
use crate::point_add::trailmix_port::circuit::{Circuit, QReg};

/// Cuccaro et al. (arXiv:quant-ph/0410184) in-place adder.
/// `a ← (a + b) mod 2^n` where a, b are both n-bit quantum registers.
/// `b` is preserved. 1 clean ancilla (carry-in) alloc'd internally.
///
/// Ancs: 1 (polylog ✓). Gates: 2n Toffoli + 4n CX.
/// Replaces `add_physical`'s n-1 AND ancillae (which violate HARD
/// RULE). MBU-based `add_physical` amortizes Toffolis via HMR+CZ in
/// UMA backward but pays with O(n) ancs; canonical Cuccaro does the
/// UMA with CCX, keeping ancs at O(1).
///
/// Structure:
///   1. MAJ cascade forward (n levels): each stage temporarily
///      stores `carry_i` in b[i]; a[i] becomes a XOR b XOR `carry_i`
///      (a partial sum).
///   2. UMA cascade reverse (n levels): restores b[i] to `b_i` and
///      finalizes a[i] = `sum_i`.
///
/// For overflow, use `add_cuccaro_with_overflow`; this variant
/// silently discards carry-out (mod 2^n arithmetic).
pub fn add_cuccaro(circ: &mut Circuit, a: &[QReg], b: &[QReg]) {
    let n = a.len();
    let nb = b.len();
    assert!(
        nb == n || nb == n - 1,
        "add_cuccaro: b must be same length as a, or 1 bit shorter \
         (treating b[n-1] as implicit 0); got a.len()={n}, b.len()={nb}"
    );

    // PRE: capture (a_pre, b_pre).
    if n > 0 {
        let a_for_capture: Vec<&QReg> = a.iter().collect();
        let b_for_capture: Vec<&QReg> = b.iter().collect();
        circ.contract_capture(
            "mbu.add_cuccaro.pre",
            move |view, shot| -> Result<(crate::point_add::trailmix_port::num_bigint::BigUint, crate::point_add::trailmix_port::num_bigint::BigUint), String> {
                let read = |regs: &[&QReg]| -> crate::point_add::trailmix_port::num_bigint::BigUint {
                    let mut v = crate::point_add::trailmix_port::num_bigint::BigUint::from(0u32);
                    for (i, q) in regs.iter().enumerate() {
                        if view.contract_read_bit_shot(q, shot) {
                            v |= crate::point_add::trailmix_port::num_bigint::BigUint::from(1u32) << i;
                        }
                    }
                    v
                };
                Ok((read(&a_for_capture), read(&b_for_capture)))
            },
        );
    }

    if n == 0 {
        return;
    }
    if n == 1 {
        if nb == 1 {
            circ.cx(&b[0], &a[0]);
        }
        // else: b is empty, b[0] is implicitly 0 → no-op.
        add_cuccaro_post_check(circ, a, b);
        return;
    }

    let c = circ.alloc_qreg("cuccaro_c");

    // MAJ_0.
    circ.cx(&b[0], &a[0]);
    circ.cx(&b[0], &c);
    circ.ccx(&c, &a[0], &b[0]);
    // MAJ_i for i in 1..n-1.
    for i in 1..n - 1 {
        circ.cx(&b[i], &a[i]);
        circ.cx(&b[i], &b[i - 1]);
        circ.ccx(&b[i - 1], &a[i], &b[i]);
    }
    // MAJ_{n-1} truncated. When b is 1 bit shorter, b[n-1] is
    // implicitly 0, so cx(b[n-1], a[n-1]) is a no-op and we skip it.
    if nb == n {
        circ.cx(&b[n - 1], &a[n - 1]);
    }

    // UMA_{n-1} truncated.
    circ.cx(&b[n - 2], &a[n - 1]);
    // UMA_i for i in (1..n-1).rev(): full UMA.
    for i in (1..n - 1).rev() {
        circ.ccx(&b[i - 1], &a[i], &b[i]);
        circ.cx(&b[i], &b[i - 1]);
        circ.cx(&b[i - 1], &a[i]);
    }
    circ.ccx(&c, &a[0], &b[0]);
    circ.cx(&b[0], &c);
    circ.cx(&c, &a[0]);

    // c drops here.

    add_cuccaro_post_check(circ, a, b);
}

/// Post-check for `add_cuccaro`: a == (`a_pre` + `b_pre`) mod 2^n; b unchanged.
fn add_cuccaro_post_check(circ: &mut Circuit, a: &[QReg], b: &[QReg]) {
    let n = a.len();
    let a_for_check: Vec<&QReg> = a.iter().collect();
    let b_for_check: Vec<&QReg> = b.iter().collect();
    circ.contract_pop_and_check::<(crate::point_add::trailmix_port::num_bigint::BigUint, crate::point_add::trailmix_port::num_bigint::BigUint), _>(
        "mbu.add_cuccaro.pre",
        move |cap, view, shot| -> Result<(), String> {
            let (a_pre, b_pre) = cap;
            let read = |regs: &[&QReg]| -> crate::point_add::trailmix_port::num_bigint::BigUint {
                let mut v = crate::point_add::trailmix_port::num_bigint::BigUint::from(0u32);
                for (i, q) in regs.iter().enumerate() {
                    if view.contract_read_bit_shot(q, shot) {
                        v |= crate::point_add::trailmix_port::num_bigint::BigUint::from(1u32) << i;
                    }
                }
                v
            };
            let a_post = read(&a_for_check);
            let b_post = read(&b_for_check);
            let modulus = crate::point_add::trailmix_port::num_bigint::BigUint::from(1u32) << n;
            let expected = (a_pre + b_pre) % &modulus;
            if a_post != expected {
                return Err(format!(
                    "add_cuccaro: a_post={a_post:#x}, expected (a_pre+b_pre) mod 2^{n} = {expected:#x} (a_pre={a_pre:#x}, b_pre={b_pre:#x})"
                ));
            }
            if &b_post != b_pre {
                return Err(format!(
                    "add_cuccaro: b changed {b_pre:#x}->{b_post:#x}"
                ));
            }
            Ok(())
        },
    );
}

/// Controlled Cuccaro add with a carry-window hook, taking register
/// slices. Lets callers pass slices that need explicit ordering (e.g. a
/// reversed slot view to operate on a BE-stored region as LE).
pub fn controlled_add_cuccaro_carry_window_refs(
    circ: &mut Circuit,
    ctrl: &QReg,
    a: &[&QReg],
    b: &[&QReg],
    window_hook: impl FnOnce(&mut Circuit, &QReg),
) {
    let n = b.len();
    assert_eq!(
        a.len(),
        n,
        "controlled_add_cuccaro_carry_window_refs: a/b length mismatch"
    );
    if n == 0 {
        return;
    }
    if n == 1 {
        // carry-out = ctrl AND a[0] AND b[0] (pre-add values).
        let cw = circ.alloc_qreg_bits("ccuc_cw1", 1);
        mcx_clean_k(circ, &[ctrl, a[0], b[0]], &cw[0]);
        window_hook(circ, &cw[0]);
        mcx_clean_k(circ, &[ctrl, a[0], b[0]], &cw[0]);
        drop(cw);
        circ.ccx(ctrl, b[0], a[0]);
        return;
    }

    let c = circ.alloc_qreg_bits("ccuc_cw_c", 1);
    let scratch = circ.alloc_qreg_bits("ccuc_cw_scratch", 1);
    let cccx = |circ: &mut Circuit, x: &QReg, y: &QReg, target: &QReg| {
        circ.ccx(ctrl, x, &scratch[0]);
        circ.ccx(&scratch[0], y, target);
        circ.ccx(ctrl, x, &scratch[0]);
    };

    // Forward MAJ cascade.
    circ.ccx(ctrl, b[0], a[0]);
    circ.ccx(ctrl, b[0], &c[0]);
    cccx(circ, &c[0], a[0], b[0]);
    for i in 1..n {
        circ.ccx(ctrl, b[i], a[i]);
        circ.ccx(ctrl, b[i], b[i - 1]);
        cccx(circ, b[i - 1], a[i], b[i]);
    }

    // Carry window: b[n-1] = ctrl AND (carry into bit n).
    window_hook(circ, b[n - 1]);

    // Backward UMA cascade.
    for i in (1..n).rev() {
        cccx(circ, b[i - 1], a[i], b[i]);
        circ.ccx(ctrl, b[i], b[i - 1]);
        circ.ccx(ctrl, b[i - 1], a[i]);
    }
    cccx(circ, &c[0], a[0], b[0]);
    drop(scratch);
    circ.ccx(ctrl, b[0], &c[0]);
    circ.ccx(ctrl, &c[0], a[0]);
}

/// Pure controlled add (3n CCX): if ctrl=1, a := a + b mod 2^n; else
/// a unchanged. Same semantics as [`controlled_add_cuccaro_mbu`] but
/// uses ~2.7x fewer Toffolis.
///
/// Semantics:
///   ctrl=1: a := (a + b) mod 2^n
///   ctrl=0: a unchanged
///   b, ctrl preserved in both cases.
///
/// Construction (same insight as [`crate::point_add::trailmix_port::arith::cuccaro_compare_act::
/// compare_and_sub_inplace_middle`], adapted to take an external
/// control instead of the captured compare-result):
///
///   FORWARD MAJ chain (1-qubit ripple, single carry register c=|0>):
///     per bit i:  CX(c, b[i]); CX(c, a[i]); CCX(a[i], b[i], c)
///     state post-bit i:
///       a[i] = `a_orig` XOR `c_in_i`
///       b[i] = `b_orig` XOR `c_in_i`
///       c    = `c_in_i` XOR `MAJ(c_in_i`, `a_orig`, `b_orig`) = `c_out_i` = `c_in`_{i+1}
///     After all n bits: c = carry-out of (a + b) >> n.
///     Cost: 1 CCX per bit, **n CCX total**.
///
///   REVERSE pass (gated on ctrl, descending i):
///     CCX(a[i], b[i], c)   ; restore c to `c_in_i`        (1 CCX)
///     CX(c, a[i])          ; a[i] := `a_orig`             (CX)
///     CCX(ctrl, b[i], a[i]); gated: a[i] XOR= `ctrl·b_i`   (1 CCX)
///                            ctrl=0: no-op  → a[i] stays `a_orig`
///                            ctrl=1: a[i] := `a_orig` XOR `b_orig` XOR `c_in_i`
///                                    = `sum_i`              ✓
///     CX(c, b[i])          ; b[i] := `b_orig`             (CX)
///     Cost: 2 CCX per bit, **2n CCX total**.
///
///   c is restored to |0> at the end (initial carry-in was 0, full ripple
///   unwinds to 0).
///
/// Total: **3n CCX** (vs 8n for `controlled_add_cuccaro_mbu`).
///
/// Why this works: the forward MAJ leaves b[i] = `b_orig` XOR `c_in_i`,
/// which is exactly the XOR-source needed to complete UMA via a single
/// `CCX(ctrl, b[i], a[i])`. The `CX(c, a)` in the reverse base case
/// unconditionally backs out the `c_in_i` contribution; the gated CCX
/// either adds in (b XOR `c_in`) when ctrl=1 (completing the sum) or
/// adds nothing when ctrl=0 (leaving `a_orig`).
///
/// Polylog peak: +1 ancilla (the single c qubit), no per-bit allocs.
/// `b` and `ctrl` are preserved.
///
/// Preconditions:
///   - `a.len()` == `b.len()` == n.
///   - ctrl NOT aliased with a or b (asserts otherwise).
pub fn controlled_add_cuccaro_3n(circ: &mut Circuit, ctrl: &QReg, a: &[QReg], b: &[QReg]) {
    let a_refs: Vec<&QReg> = a.iter().collect();
    let b_refs: Vec<&QReg> = b.iter().collect();
    controlled_add_cuccaro_3n_refs(circ, ctrl, &a_refs, &b_refs);
}

/// Reference-slice variant of [`controlled_add_cuccaro_3n`].
///
/// Same semantics (if ctrl=1, a := a+b mod 2^n; else unchanged), and same
/// **3n CCX** cost. Use when the caller already holds borrows (e.g. a
/// reverse-physical view onto an MSB-anchored packed register).
pub fn controlled_add_cuccaro_3n_refs(circ: &mut Circuit, ctrl: &QReg, a: &[&QReg], b: &[&QReg]) {
    let n = a.len();
    assert_eq!(
        b.len(),
        n,
        "controlled_add_cuccaro_3n_refs: a/b length mismatch"
    );

    let aliases_a = a.iter().any(|q| std::ptr::eq(*q, ctrl));
    let aliases_b = b.iter().any(|q| std::ptr::eq(*q, ctrl));
    assert!(
        !aliases_a,
        "controlled_add_cuccaro_3n_refs: ctrl aliases a -- unsupported"
    );
    assert!(
        !aliases_b,
        "controlled_add_cuccaro_3n_refs: ctrl aliases b -- unsupported"
    );

    // PRE: capture (a_pre, b_pre, ctrl_pre).
    if n > 0 {
        let a_for_capture: Vec<&QReg> = a.to_vec();
        let b_for_capture: Vec<&QReg> = b.to_vec();
        let ctrl_ref = ctrl;
        circ.contract_capture(
            "mbu.controlled_add_cuccaro_3n.pre",
            move |view, shot| -> Result<(crate::point_add::trailmix_port::num_bigint::BigUint, crate::point_add::trailmix_port::num_bigint::BigUint, bool), String> {
                let read = |regs: &[&QReg]| -> crate::point_add::trailmix_port::num_bigint::BigUint {
                    let mut v = crate::point_add::trailmix_port::num_bigint::BigUint::from(0u32);
                    for (i, q) in regs.iter().enumerate() {
                        if view.contract_read_bit_shot(q, shot) {
                            v |= crate::point_add::trailmix_port::num_bigint::BigUint::from(1u32) << i;
                        }
                    }
                    v
                };
                Ok((
                    read(&a_for_capture),
                    read(&b_for_capture),
                    view.contract_read_bit_shot(ctrl_ref, shot),
                ))
            },
        );
    }

    if n == 0 {
        return;
    }
    if n == 1 {
        // 1-bit case: a[0] ^= ctrl·b[0].
        circ.ccx(ctrl, b[0], a[0]);
        controlled_add_cuccaro_3n_post_check_refs(circ, ctrl, a, b);
        return;
    }

    let c = circ.alloc_qreg("ccuccaro3n_c");

    // Forward MAJ chain. Single carry register `c` ripples bit-by-bit.
    // Per bit i: CX(c, b); CX(c, a); CCX(a, b, c).
    for i in 0..n {
        circ.cx(&c, b[i]);
        circ.cx(&c, a[i]);
        circ.ccx(a[i], b[i], &c);
    }

    // Reverse pass, descending. Per bit i:
    //   CCX(a,b,c)       -- restore c to c_in_i
    //   CX(c, a)          -- a := a_orig (undo)
    //   CCX(ctrl, b, a)   -- gated: a XOR= ctrl·(b XOR c_in) = ctrl·b_orig XOR ctrl·c_in
    //   CX(c, b)          -- b := b_orig
    //
    // When ctrl=1, the chained CX(c,a) then CCX(ctrl,b,a) yields
    //   a_post = a_orig XOR (b_orig XOR c_in) = a_orig XOR b_orig XOR c_in = sum_i.
    // When ctrl=0, the CCX is a no-op, so a_post = a_orig.
    for i in (0..n).rev() {
        circ.ccx(a[i], b[i], &c);
        circ.cx(&c, a[i]);
        circ.ccx(ctrl, b[i], a[i]);
        circ.cx(&c, b[i]);
    }

    // c is back to |0> (initial carry-in was 0; ripple fully unwound).
    circ.zero_and_free(c);

    controlled_add_cuccaro_3n_post_check_refs(circ, ctrl, a, b);
}

/// LITERAL gate-by-gate inverse of `controlled_add_cuccaro_3n_refs`.
/// Emits the SAME gates in EXACT REVERSE order — not an algebraically
/// equivalent subtract circuit (like X-sandwich Cuccaro), but the
/// bit-for-bit inverted gate sequence.
///
/// This matters for drift cancellation in approximate-primitive
/// composition (e.g. Schrottenloher Alg 4's `apply_bitvector` inverse):
/// the forward primitive contributes drift that the X-sandwich form
/// cannot cancel, but the literal gate-inverse DOES cancel exactly.
///
/// Semantics on a state in the image of forward: takes (`a_post`, b, ctrl)
/// where `a_post` = `a_pre` + ctrl·b, returns (`a_pre`, b, ctrl). Cost
/// matches forward: 3n CCX.
pub fn controlled_add_cuccaro_3n_reverse_refs(
    circ: &mut Circuit,
    ctrl: &QReg,
    a: &[&QReg],
    b: &[&QReg],
) {
    let n = a.len();
    assert_eq!(
        b.len(),
        n,
        "controlled_add_cuccaro_3n_reverse_refs: a/b length mismatch"
    );

    let aliases_a = a.iter().any(|q| std::ptr::eq(*q, ctrl));
    let aliases_b = b.iter().any(|q| std::ptr::eq(*q, ctrl));
    assert!(!aliases_a, "ctrl aliases a -- unsupported");
    assert!(!aliases_b, "ctrl aliases b -- unsupported");

    if n == 0 {
        return;
    }
    if n == 1 {
        // 1-bit forward is just ccx(ctrl, b[0], a[0]); CCX is self-inverse.
        circ.ccx(ctrl, b[0], a[0]);
        return;
    }

    let c = circ.alloc_qreg("ccuccaro3n_c_rev");

    // Inverse of the forward's reverse pass (descending, with gate
    // order reversed within each bit):
    //   forward emitted, for i = n-1 down to 0:
    //     ccx(a,b,c); cx(c,a); ccx(ctrl,b,a); cx(c,b)
    //   inverse emits, for i = 0 up to n-1:
    //     cx(c,b); ccx(ctrl,b,a); cx(c,a); ccx(a,b,c)
    for i in 0..n {
        circ.cx(&c, b[i]);
        circ.ccx(ctrl, b[i], a[i]);
        circ.cx(&c, a[i]);
        circ.ccx(a[i], b[i], &c);
    }

    // Inverse of the forward MAJ chain (ascending, gates reversed
    // within bit):
    //   forward emitted, for i = 0 up to n-1:
    //     cx(c,b); cx(c,a); ccx(a,b,c)
    //   inverse emits, for i = n-1 down to 0:
    //     ccx(a,b,c); cx(c,a); cx(c,b)
    for i in (0..n).rev() {
        circ.ccx(a[i], b[i], &c);
        circ.cx(&c, a[i]);
        circ.cx(&c, b[i]);
    }

    circ.zero_and_free(c);
}

/// Convenience wrapper: literal gate-inverse of `controlled_add_cuccaro_3n`.
pub fn controlled_add_cuccaro_3n_reverse(circ: &mut Circuit, ctrl: &QReg, a: &[QReg], b: &[QReg]) {
    let a_refs: Vec<&QReg> = a.iter().collect();
    let b_refs: Vec<&QReg> = b.iter().collect();
    controlled_add_cuccaro_3n_reverse_refs(circ, ctrl, &a_refs, &b_refs);
}

/// Post-check for `controlled_add_cuccaro_3n_refs`.
fn controlled_add_cuccaro_3n_post_check_refs(
    circ: &mut Circuit,
    ctrl: &QReg,
    a: &[&QReg],
    b: &[&QReg],
) {
    let n = a.len();
    let a_for_check: Vec<&QReg> = a.to_vec();
    let b_for_check: Vec<&QReg> = b.to_vec();
    let ctrl_ref = ctrl;
    circ.contract_pop_and_check::<(crate::point_add::trailmix_port::num_bigint::BigUint, crate::point_add::trailmix_port::num_bigint::BigUint, bool), _>(
        "mbu.controlled_add_cuccaro_3n.pre",
        move |cap, view, shot| -> Result<(), String> {
            let (a_pre, b_pre, c_pre) = cap;
            let read = |regs: &[&QReg]| -> crate::point_add::trailmix_port::num_bigint::BigUint {
                let mut v = crate::point_add::trailmix_port::num_bigint::BigUint::from(0u32);
                for (i, q) in regs.iter().enumerate() {
                    if view.contract_read_bit_shot(q, shot) {
                        v |= crate::point_add::trailmix_port::num_bigint::BigUint::from(1u32) << i;
                    }
                }
                v
            };
            let a_post = read(&a_for_check);
            let b_post = read(&b_for_check);
            let c_post = view.contract_read_bit_shot(ctrl_ref, shot);
            let modulus = crate::point_add::trailmix_port::num_bigint::BigUint::from(1u32) << n;
            let addend = if *c_pre { b_pre.clone() } else { crate::point_add::trailmix_port::num_bigint::BigUint::from(0u32) };
            let expected = (a_pre + &addend) % &modulus;
            if a_post != expected {
                return Err(format!(
                    "ctrl_add_cuccaro_3n: a_post={:#x}, expected {:#x} (a_pre={:#x}, b_pre={:#x}, ctrl={})",
                    a_post, expected, a_pre, b_pre, u8::from(*c_pre),
                ));
            }
            if &b_post != b_pre {
                return Err(format!("ctrl_add_cuccaro_3n: b changed {b_pre:#x}->{b_post:#x}"));
            }
            if c_post != *c_pre {
                return Err(format!("ctrl_add_cuccaro_3n: ctrl changed {} -> {}", u8::from(*c_pre), u8::from(c_post)));
            }
            Ok(())
        },
    );
}

/// Variant of [`controlled_add_cuccaro`] that, in addition to
/// preserving b[..n-1] as the standard adder does, FREES `b[n-1]`
/// immediately after its last gate-touch inside Cuccaro UMA.
/// Caller asserts (via the free's sim mask check) that b[n-1]
/// was |0> on entry — UMA restores b[n-1] to its input value, so
/// this is the only valid case for consume.
///
/// Use when the caller's outer loop retires the top bit of the
/// b slice each iteration (e.g. `multi_sub`'s iter j retires b[L-1]
/// where L = slice length for that iter).
pub fn controlled_add_cuccaro_consume_top_b(
    circ: &mut Circuit,
    ctrl: &QReg,
    a: &[QReg],
    b: &[QReg],
) {
    let n = b.len();
    assert_eq!(
        a.len(),
        n,
        "controlled_add_cuccaro_consume_top_b: a/b length mismatch"
    );
    if n == 0 {
        return;
    }
    if n == 1 {
        // Only one bit — controlled_add of b[0] into a[0]; b[0] is
        // unchanged. Caller is expected to drop b[0] after this call.
        circ.ccx(ctrl, &b[0], &a[0]);
        return;
    }

    let c = circ.alloc_qreg_bits("ccuccaro_c", 1);

    // Streaming MBU-AND cccx (same pattern as
    // [`controlled_add_cuccaro_mbu`]).
    let mbu_cccx = |circ: &mut Circuit, x: &QReg, y: &QReg, target: &QReg| {
        let anc = circ.alloc_qreg_bits("ccuccaro_mbu_and", 1);
        circ.ccx(ctrl, x, &anc[0]);
        circ.ccx(&anc[0], y, target);
        let bit = circ.alloc_bit();
        circ.hmr(&anc[0], bit);
        circ.cz_if_bit(ctrl, x, bit);
        circ.free_bit(bit);
        drop(anc);
    };

    circ.ccx(ctrl, &b[0], &a[0]);
    circ.ccx(ctrl, &b[0], &c[0]);
    mbu_cccx(circ, &c[0], &a[0], &b[0]);

    for i in 1..n {
        circ.ccx(ctrl, &b[i], &a[i]);
        circ.ccx(ctrl, &b[i], &b[i - 1]);
        mbu_cccx(circ, &b[i - 1], &a[i], &b[i]);
    }

    // UMA cascade. Caller is expected to drop b[n-1] after this call.
    let top = n - 1;
    mbu_cccx(circ, &b[top - 1], &a[top], &b[top]);
    circ.ccx(ctrl, &b[top], &b[top - 1]);
    circ.ccx(ctrl, &b[top - 1], &a[top]);
    for i in (1..top).rev() {
        mbu_cccx(circ, &b[i - 1], &a[i], &b[i]);
        circ.ccx(ctrl, &b[i], &b[i - 1]);
        circ.ccx(ctrl, &b[i - 1], &a[i]);
    }
    mbu_cccx(circ, &c[0], &a[0], &b[0]);
    circ.ccx(ctrl, &b[0], &c[0]);
    circ.ccx(ctrl, &c[0], &a[0]);
}

/// Controlled Cuccaro adder with overflow. If ctrl=1:
/// `a_ext`[0..n] ← (a+b) mod 2^n, `a_ext`[n] ← (a+b) div 2^n.
/// If ctrl=0: unchanged.
///
/// Streaming MBU-AND form (same pattern as
/// [`controlled_add_cuccaro_mbu`]): each cccx (target ^= ctrl·x·y)
/// uses 2 CCX + HMR + `cz_if_bit` instead of 3-CCX clean-anc form.
/// Saves 1 CCX per cccx (2n+1 cccx invocations → ≈2n CCX saved per
/// adder). Polylog peak preserved (per-cccx anc allocated and freed).
pub fn controlled_add_cuccaro_with_overflow(
    circ: &mut Circuit,
    ctrl: &QReg,
    a_ext: &[QReg],
    b: &[QReg],
) {
    let n = b.len();
    assert!(a_ext.len() > n, "a_ext must have n+1 bits for overflow");
    if n == 0 {
        return;
    }
    let a = &a_ext[..n];
    let ovf = &a_ext[n];

    if n == 1 {
        // if ctrl: ovf ^= a[0]·b[0], a[0] ^= b[0]
        let dirty = circ.alloc_qreg_bits("c_ovf_scratch", 1);
        mcx_dirty(circ, &[ctrl, &a[0], &b[0]], ovf, &dirty[0]);
        circ.ccx(ctrl, &b[0], &a[0]);
        return;
    }

    let c = circ.alloc_qreg_bits("ccuccaro_c", 1);

    // Streaming MBU-AND cccx: target ^= ctrl AND x AND y. Allocates
    // a fresh anc per call, uncomputes via HMR + cz_if_bit. See
    // `controlled_add_cuccaro_mbu` for the identity-discharge proof.
    let mbu_cccx = |circ: &mut Circuit, x: &QReg, y: &QReg, target: &QReg| {
        let anc = circ.alloc_qreg_bits("ccuccaro_mbu_and", 1);
        circ.ccx(ctrl, x, &anc[0]); // anc = ctrl AND x
        circ.ccx(&anc[0], y, target); // target ^= ctrl·x·y
        let bit = circ.alloc_bit();
        circ.hmr(&anc[0], bit); // anc ← 0; obligation logged
        circ.cz_if_bit(ctrl, x, bit); // discharge AndOf(ctrl, x)
        circ.free_bit(bit);
        drop(anc);
    };

    // MAJ cascade.
    circ.ccx(ctrl, &b[0], &a[0]);
    circ.ccx(ctrl, &b[0], &c[0]);
    mbu_cccx(circ, &c[0], &a[0], &b[0]);
    for i in 1..n {
        circ.ccx(ctrl, &b[i], &a[i]);
        circ.ccx(ctrl, &b[i], &b[i - 1]);
        mbu_cccx(circ, &b[i - 1], &a[i], &b[i]);
    }

    // Capture overflow: b[n-1] currently holds carry_out (post-MAJ).
    circ.ccx(ctrl, &b[n - 1], ovf);

    // UMA cascade.
    for i in (1..n).rev() {
        mbu_cccx(circ, &b[i - 1], &a[i], &b[i]);
        circ.ccx(ctrl, &b[i], &b[i - 1]);
        circ.ccx(ctrl, &b[i - 1], &a[i]);
    }
    mbu_cccx(circ, &c[0], &a[0], &b[0]);
    circ.ccx(ctrl, &b[0], &c[0]);
    circ.ccx(ctrl, &c[0], &a[0]);
}

/// Cuccaro adder with explicit overflow bit. Receiver first:
///   `a_ext`: n+1 bits with high bit `|0⟩`; receives the sum in the low n bits
///            and the carry-out in `a_ext[n]`.
///   `b`: n bits, preserved addend.
pub fn add_cuccaro_with_overflow(circ: &mut Circuit, a_ext: &[QReg], b: &[QReg]) {
    let n = b.len();
    assert!(a_ext.len() > n, "a_ext must have n+1 bits for overflow");
    if n == 0 {
        return;
    }
    let a = &a_ext[..n];
    let ovf = &a_ext[n];

    if n == 1 {
        // 1-bit add with overflow: (a+b) mod 2 stored in a[0]; overflow = a AND b.
        // Actually: a[0] ← a[0] XOR b[0], ovf ← a·b.
        // Order matters: compute ovf first (using original values).
        circ.ccx(&a[0], &b[0], ovf);
        circ.cx(&b[0], &a[0]);
        return;
    }

    let c = circ.alloc_qreg("cuccaro_c");

    // MAJ cascade.
    circ.cx(&b[0], &a[0]);
    circ.cx(&b[0], &c);
    circ.ccx(&c, &a[0], &b[0]);
    for i in 1..n {
        circ.cx(&b[i], &a[i]);
        circ.cx(&b[i], &b[i - 1]);
        circ.ccx(&b[i - 1], &a[i], &b[i]);
    }

    // Capture overflow: after MAJ cascade, b[n-1] = carry_n.
    circ.cx(&b[n - 1], ovf);

    // UMA cascade.
    for i in (1..n).rev() {
        circ.ccx(&b[i - 1], &a[i], &b[i]);
        circ.cx(&b[i], &b[i - 1]);
        circ.cx(&b[i - 1], &a[i]);
    }
    circ.ccx(&c, &a[0], &b[0]);
    circ.cx(&b[0], &c);
    circ.cx(&c, &a[0]);

    // c drops here.
}
