//! Reversible UNPACKED PZ modular inversion. Separate registers A,B,|a|,|b| (no
//! cursor), reversible via the PZ cofactor ratio (no spooky pebbling). The whole
//! thing is built from ONE primitive -- restoring long division -- used forward
//! on the gcd pair and (its reverse) as the cofactor multiply.
//!
//! `long_division(A,B,q)`: A := A mod B, q := A // B (q starts |0>). Reversible.
//! `long_division_reverse(A,B,q)`: the inverse -- A := A + q*B, q := 0 (consumes
//! q). This is exactly the consuming multiply `|a| += q|b|` applied to (|a|,|b|).

use crate::point_add::trailmix_port::arith::cuccaro::{
    controlled_add_cuccaro_3n_refs, controlled_add_cuccaro_3n_refs_with_carry,
};
use crate::point_add::trailmix_port::arith::mcx::mcx_dirty_ladder;
use crate::point_add::trailmix_port::circuit::{Circuit, QReg};

/// The four unpacked PZ-EEA registers. gcd pair (`a_gcd=A`, `b_gcd=B`) shrinks;
/// cofactor pair (ca, cb) grows. Init A=P, B=dx, ca=0, cb=1; at the end
/// dx^{-1} == (cb - ca) mod p (one cofactor is 0, sign baked into which).
pub struct PzRegs {
    pub a_gcd: Vec<QReg>,
    pub b_gcd: Vec<QReg>,
    pub ca: Vec<QReg>,
    pub cb: Vec<QReg>,
}

fn borrow_compare_refs_inner(
    circ: &mut Circuit,
    v: &[&QReg],
    u: &[&QReg],
    active: Option<&QReg>,
    out: &QReg,
    borrowed_carry: Option<&QReg>,
) {
    let n = v.len();
    assert_eq!(u.len(), n);
    if n == 0 {
        return;
    }
    let pcmp = circ.push_section("p.cmp");
    for q in v {
        circ.x(q); // v -> ~v
    }
    let a = u; // accumulator
    let b = v; // = ~v
    let mut owned_carry = borrowed_carry.is_none().then(|| circ.alloc_qreg("bc.c"));
    let cc = borrowed_carry.unwrap_or_else(|| owned_carry.as_ref().unwrap());
    assert!(!std::ptr::eq(cc, out), "comparator carry aliases output");
    assert!(
        !v.iter().chain(u).any(|&q| std::ptr::eq(q, cc)),
        "comparator carry aliases an operand"
    );
    if borrowed_carry.is_some() {
        let active_ref = active;
        let out_ref = out;
        let carry_ref = cc;
        let v_for_capture = v.to_vec();
        let u_for_capture = u.to_vec();
        circ.contract_capture(
            "pz.borrow_compare.borrowed.pre",
            move |view, shot| -> Result<_, String> {
                let read = |regs: &[&QReg]| {
                    let mut value =
                        crate::point_add::trailmix_port::num_bigint::BigUint::from(0u32);
                    for (i, q) in regs.iter().enumerate() {
                        if view.contract_read_bit_shot(q, shot) {
                            value |= crate::point_add::trailmix_port::num_bigint::BigUint::from(
                                1u32,
                            ) << i;
                        }
                    }
                    value
                };
                Ok((
                    active_ref.is_none_or(|q| view.contract_read_bit_shot(q, shot)),
                    view.contract_read_bit_shot(carry_ref, shot),
                    view.contract_read_bit_shot(out_ref, shot),
                    read(&v_for_capture),
                    read(&u_for_capture),
                ))
            },
        );
    }
    circ.cx(b[0], a[0]);
    circ.cx(b[0], cc);
    circ.ccx(cc, a[0], b[0]);
    for i in 1..n {
        circ.cx(b[i], a[i]);
        circ.cx(b[i], b[i - 1]);
        circ.ccx(b[i - 1], a[i], b[i]);
    }
    if let Some(active) = active {
        circ.ccx(active, b[n - 1], out);
    } else {
        circ.cx(b[n - 1], out);
    }
    for i in (1..n).rev() {
        circ.ccx(b[i - 1], a[i], b[i]);
        circ.cx(b[i], b[i - 1]);
        circ.cx(b[i], a[i]);
    }
    circ.ccx(cc, a[0], b[0]);
    circ.cx(b[0], cc);
    circ.cx(b[0], a[0]);
    if borrowed_carry.is_none() {
        circ.zero_and_free(
            owned_carry
                .take()
                .expect("owned comparator carry disappeared"),
        );
    }
    for q in v {
        circ.x(q);
    }
    if let Some(borrowed_cc) = borrowed_carry {
        let active_ref = active;
        let out_ref = out;
        let carry_ref = borrowed_cc;
        let v_for_check = v.to_vec();
        let u_for_check = u.to_vec();
        circ.contract_pop_and_check::<
            (
                bool,
                bool,
                bool,
                crate::point_add::trailmix_port::num_bigint::BigUint,
                crate::point_add::trailmix_port::num_bigint::BigUint,
            ),
            _,
        >(
            "pz.borrow_compare.borrowed.pre",
            move |cap, view, shot| -> Result<(), String> {
                let (enabled_pre, carry_pre, out_pre, v_pre, u_pre) = cap;
                if *enabled_pre && *carry_pre {
                    return Err(
                        "borrowed comparator carry must be zero whenever active=1".to_string(),
                    );
                }
                let read = |regs: &[&QReg]| {
                    let mut value =
                        crate::point_add::trailmix_port::num_bigint::BigUint::from(0u32);
                    for (i, q) in regs.iter().enumerate() {
                        if view.contract_read_bit_shot(q, shot) {
                            value |= crate::point_add::trailmix_port::num_bigint::BigUint::from(
                                1u32,
                            ) << i;
                        }
                    }
                    value
                };
                let carry_post = view.contract_read_bit_shot(carry_ref, shot);
                let out_post = view.contract_read_bit_shot(out_ref, shot);
                let active_post = active_ref.is_none_or(|q| view.contract_read_bit_shot(q, shot));
                let expected_out = *out_pre ^ (*enabled_pre && v_pre < u_pre);
                if carry_post != *carry_pre
                    || active_post != *enabled_pre
                    || out_post != expected_out
                    || &read(&v_for_check) != v_pre
                    || &read(&u_for_check) != u_pre
                {
                    return Err("borrowed comparator contract failed".to_string());
                }
                Ok(())
            },
        );
    }
    circ.pop_section(&pcmp);
}

/// `out ^= (v < u)` (MAJ cascade + carry capture + un-MAJ), v,u restored. Refs
/// variant of `two_cursor::borrow_compare` (windows here are non-contiguous).
pub(crate) fn borrow_compare_refs(circ: &mut Circuit, v: &[&QReg], u: &[&QReg], out: &QReg) {
    borrow_compare_refs_inner(circ, v, u, None, out, None);
}

/// `out ^= (v < u)` using a caller-provided clean carry lane. The carry is
/// restored to zero before return.
pub(crate) fn borrow_compare_refs_with_carry(
    circ: &mut Circuit,
    v: &[&QReg],
    u: &[&QReg],
    out: &QReg,
    carry: &QReg,
) {
    assert!(!std::ptr::eq(out, carry), "comparator carry aliases output");
    assert!(
        !v.iter()
            .chain(u)
            .any(|&q| std::ptr::eq(q, out) || std::ptr::eq(q, carry)),
        "comparator output/carry aliases an operand"
    );
    borrow_compare_refs_inner(circ, v, u, None, out, Some(carry));
}

/// `out ^= active AND (v < u)`, with `v`, `u`, and `active` restored.
///
/// The active control is applied directly to the comparator's final carry. This
/// avoids retaining a separate `(v < u)` result qubit across the gated body.
pub(crate) fn borrow_compare_gated_refs(
    circ: &mut Circuit,
    v: &[&QReg],
    u: &[&QReg],
    active: &QReg,
    out: &QReg,
) {
    assert!(!std::ptr::eq(active, out), "gated compare control aliases output");
    assert!(
        !v.iter().chain(u).any(|&q| std::ptr::eq(q, active) || std::ptr::eq(q, out)),
        "gated compare control/output aliases an operand"
    );
    borrow_compare_refs_inner(circ, v, u, Some(active), out, None);
}

/// `out ^= active AND (v < u)` using a conditionally clean caller carry.
/// The carry may be arbitrary when `active=0`; it must be zero when `active=1`.
/// It is restored to its input value before return.
pub(crate) fn borrow_compare_gated_refs_with_carry(
    circ: &mut Circuit,
    v: &[&QReg],
    u: &[&QReg],
    active: &QReg,
    out: &QReg,
    carry: &QReg,
) {
    assert!(
        !std::ptr::eq(active, out),
        "gated compare control aliases output"
    );
    assert!(!std::ptr::eq(active, carry), "comparator carry aliases active");
    assert!(!std::ptr::eq(out, carry), "comparator carry aliases output");
    assert!(
        !v.iter()
            .chain(u)
            .any(|&q| std::ptr::eq(q, active) || std::ptr::eq(q, out) || std::ptr::eq(q, carry)),
        "gated compare control/output/carry aliases an operand"
    );
    borrow_compare_refs_inner(circ, v, u, Some(active), out, Some(carry));
}

/// Expose the final borrow bit only while the comparator is in its physical
/// middle. The caller-provided carry and both operands are restored exactly;
/// `body` must restore anything it borrows before returning.
///
/// This avoids materializing a persistent comparison-output lane when the
/// result is consumed immediately by a reversible multi-control toggle.
pub(crate) fn borrow_compare_middle_refs_with_carry<
    F: FnOnce(&mut Circuit, &QReg),
>(
    circ: &mut Circuit,
    v: &[&QReg],
    u: &[&QReg],
    carry: &QReg,
    body: F,
) {
    let n = v.len();
    assert_eq!(u.len(), n);
    assert!(n > 0, "borrow comparator middle requires a nonempty width");
    assert!(
        !v.iter()
            .chain(u)
            .any(|&q| std::ptr::eq(q, carry)),
        "comparator carry aliases an operand"
    );

    let pcmp = circ.push_section("p.cmp.middle");
    for q in v.iter().rev() {
        circ.x(q);
    }
    let a = u;
    let b = v;
    circ.cx(b[0], a[0]);
    circ.cx(b[0], carry);
    circ.ccx(carry, a[0], b[0]);
    for i in 1..n {
        circ.cx(b[i], a[i]);
        circ.cx(b[i], b[i - 1]);
        circ.ccx(b[i - 1], a[i], b[i]);
    }

    body(circ, b[n - 1]);

    for i in (1..n).rev() {
        circ.ccx(b[i - 1], a[i], b[i]);
        circ.cx(b[i], b[i - 1]);
        circ.cx(b[i], a[i]);
    }
    circ.ccx(carry, a[0], b[0]);
    circ.cx(b[0], carry);
    circ.cx(b[0], a[0]);
    for q in v {
        circ.x(q);
    }
    circ.pop_section(&pcmp);
}

/// a += b (mod 2^len) gated on `g`. Plain controlled Cuccaro (3n).
pub(crate) fn ctrl_add(c: &mut Circuit, g: &QReg, a: &[&QReg], b: &[&QReg]) {
    let prev = c.push_section("p.add");
    controlled_add_cuccaro_3n_refs(c, g, a, b);
    c.pop_section(&prev);
}

/// `a += g*b` with a caller carry that is clean only on the active branch.
pub(crate) fn ctrl_add_with_carry(
    c: &mut Circuit,
    g: &QReg,
    a: &[&QReg],
    b: &[&QReg],
    carry: &QReg,
) {
    let prev = c.push_section("p.add");
    controlled_add_cuccaro_3n_refs_with_carry(c, g, a, b, carry);
    c.pop_section(&prev);
}

/// a -= b (mod 2^len) gated on `g` (X-bracket + controlled add). PRE when g: a>=b.
pub(crate) fn ctrl_sub(c: &mut Circuit, g: &QReg, a: &[&QReg], b: &[&QReg]) {
    let prev = c.push_section("p.sub");
    for q in a {
        c.x(q);
    }
    controlled_add_cuccaro_3n_refs(c, g, a, b);
    for q in a {
        c.x(q);
    }
    c.pop_section(&prev);
}

/// `a -= g*b` with a caller carry that is clean only on the active branch.
pub(crate) fn ctrl_sub_with_carry(
    c: &mut Circuit,
    g: &QReg,
    a: &[&QReg],
    b: &[&QReg],
    carry: &QReg,
) {
    let prev = c.push_section("p.sub");
    for q in a {
        c.x(q);
    }
    controlled_add_cuccaro_3n_refs_with_carry(c, g, a, b, carry);
    for q in a {
        c.x(q);
    }
    c.pop_section(&prev);
}

/// `a += g*b (mod 2^n)` without allocating a carry qubit.
///
/// Each addend bit controls an increment of the corresponding suffix of `a`.
/// The increment is emitted from high to low so its controls see the
/// pre-increment lower suffix. Multi-controlled X gates borrow the other bits
/// of `b` as dirty lenders and restore them exactly. This route is intended for
/// the five-bit hybrid-CLZ transcript update, where saving Cuccaro's one clean
/// carry removes the global peak qubit at modest gate cost.
pub(crate) fn ctrl_add_dirty_lenders(
    c: &mut Circuit,
    g: &QReg,
    a: &[&QReg],
    b: &[&QReg],
) {
    let prev = c.push_section("p.add");
    assert_eq!(a.len(), b.len(), "ctrl_add_dirty_lenders width mismatch");
    let n = a.len();

    for i in 0..n {
        let dirty: Vec<&QReg> = b
            .iter()
            .enumerate()
            .filter_map(|(k, &q)| (k != i).then_some(q))
            .collect();
        for j in (i..n).rev() {
            let mut ctrls = Vec::with_capacity(j - i + 2);
            ctrls.push(g);
            ctrls.push(b[i]);
            ctrls.extend_from_slice(&a[i..j]);
            mcx_dirty_ladder(c, &ctrls, a[j], &dirty);
        }
    }
    c.pop_section(&prev);
}

/// `a -= g*b (mod 2^n)` using the allocation-free dirty-lender adder.
pub(crate) fn ctrl_sub_dirty_lenders(
    c: &mut Circuit,
    g: &QReg,
    a: &[&QReg],
    b: &[&QReg],
) {
    let prev = c.push_section("p.sub");
    for q in a {
        c.x(q);
    }
    ctrl_add_dirty_lenders(c, g, a, b);
    for q in a {
        c.x(q);
    }
    c.pop_section(&prev);
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct DirtyLenderExhaustiveReport {
    pub widths_checked: usize,
    pub input_states_checked: usize,
    pub gate_simulations: usize,
    pub width5_add_ops: usize,
    pub width5_add_toffoli: usize,
    pub width5_sub_ops: usize,
    pub width5_sub_toffoli: usize,
}

/// Interpret and exhaustively verify the emitted add/sub streams for widths
/// one through five. This covers every control, accumulator, and dirty-lender
/// state. Any allocation or non-classical gate fails closed.
#[doc(hidden)]
pub fn exhaustive_dirty_lender_check() -> DirtyLenderExhaustiveReport {
    use crate::circuit::{Op, OperationType};

    fn apply(ops: &[Op], mut state: u64) -> u64 {
        let bit = |state: u64, id: u64| ((state >> id) & 1) != 0;
        for op in ops {
            match op.kind {
                OperationType::X => state ^= 1u64 << op.q_target.0,
                OperationType::CX => {
                    if bit(state, op.q_control1.0) {
                        state ^= 1u64 << op.q_target.0;
                    }
                }
                OperationType::CCX => {
                    if bit(state, op.q_control1.0) && bit(state, op.q_control2.0) {
                        state ^= 1u64 << op.q_target.0;
                    }
                }
                other => panic!("dirty-lender primitive emitted unexpected gate {other:?}"),
            }
        }
        state
    }

    fn build(n: usize, subtract: bool) -> crate::point_add::B {
        let mut c = Circuit::new();
        let g = c.alloc_qreg("dirty-check.g");
        let a = c.alloc_qreg_bits("dirty-check.a", n);
        let b = c.alloc_qreg_bits("dirty-check.b", n);
        let a_refs: Vec<&QReg> = a.iter().collect();
        let b_refs: Vec<&QReg> = b.iter().collect();
        if subtract {
            ctrl_sub_dirty_lenders(&mut c, &g, &a_refs, &b_refs);
        } else {
            ctrl_add_dirty_lenders(&mut c, &g, &a_refs, &b_refs);
        }
        drop(a_refs);
        drop(b_refs);
        let builder = c.into_builder();
        assert_eq!(builder.next_qubit as usize, 2 * n + 1);
        assert_eq!(builder.peak_qubits as usize, 2 * n + 1);
        assert_eq!(builder.active_qubits as usize, 2 * n + 1);
        drop((g, a, b));
        builder
    }

    let mut input_states_checked = 0usize;
    let mut gate_simulations = 0usize;
    let mut width5_add_ops = 0usize;
    let mut width5_add_toffoli = 0usize;
    let mut width5_sub_ops = 0usize;
    let mut width5_sub_toffoli = 0usize;

    for n in 1..=5usize {
        let add = build(n, false);
        let sub = build(n, true);
        let states = 1usize << (2 * n + 1);
        input_states_checked += states;
        gate_simulations += 2 * states;
        if n == 5 {
            width5_add_ops = add.ops.len();
            width5_add_toffoli = add
                .ops
                .iter()
                .filter(|op| op.kind == OperationType::CCX)
                .count();
            width5_sub_ops = sub.ops.len();
            width5_sub_toffoli = sub
                .ops
                .iter()
                .filter(|op| op.kind == OperationType::CCX)
                .count();
        }

        let mask = (1u64 << n) - 1;
        for input in 0..states as u64 {
            let g = input & 1;
            let a = (input >> 1) & mask;
            let b = (input >> (n + 1)) & mask;
            for (subtract, builder) in [(false, &add), (true, &sub)] {
                let output = apply(&builder.ops, input);
                let got_g = output & 1;
                let got_a = (output >> 1) & mask;
                let got_b = (output >> (n + 1)) & mask;
                let want_a = if g == 0 {
                    a
                } else if subtract {
                    a.wrapping_sub(b) & mask
                } else {
                    a.wrapping_add(b) & mask
                };
                assert_eq!(got_g, g, "width={n} subtract={subtract}: control changed");
                assert_eq!(got_b, b, "width={n} subtract={subtract}: lender changed");
                assert_eq!(
                    got_a, want_a,
                    "width={n} subtract={subtract} g={g} a={a} b={b}"
                );
            }
        }
    }

    DirtyLenderExhaustiveReport {
        widths_checked: 5,
        input_states_checked,
        gate_simulations,
        width5_add_ops,
        width5_add_toffoli,
        width5_sub_ops,
        width5_sub_toffoli,
    }
}

/// Restoring long division. `a` (n qubits, value < 2^n), `b` (m qubits, 0<b<2^m),
/// `q` (n-m+1 qubits, |0>). After: a holds (a mod b) in [0,m), a[m..n)=0; q = a//b.
/// Per quotient position j (high to low): window w = a[j..j+m] ++ guard; set
/// q[j] = (w >= b); if q[j] subtract b from w. Reversible; reverse =
/// [`long_division_reverse`].
pub fn long_division(c: &mut Circuit, a: &[QReg], b: &[QReg], q: &[QReg]) {
    let n = a.len();
    let m = b.len();
    assert_eq!(q.len(), n - m + 1, "q width must be n-m+1");
    let bguard = c.alloc_qreg("ld.bguard"); // bext top bit (|0>)
    let wguard = c.alloc_qreg("ld.wguard"); // window top bit when j=n-m (|0>)
    let bext: Vec<&QReg> = b.iter().chain(std::iter::once(&bguard)).collect(); // m+1, top 0
    for j in (0..=n - m).rev() {
        let mut win: Vec<&QReg> = a[j..(j + m).min(n)].iter().collect();
        if j + m < n {
            win.push(&a[j + m]); // real high bit
        } else {
            win.push(&wguard); // top: separate alloc'd guard (disjoint from bext)
        }
        debug_assert_eq!(win.len(), m + 1);
        // q[j] = (win >= b): borrow_compare gives (win < b); X to flip.
        borrow_compare_refs(c, &win, &bext, &q[j]);
        c.x(&q[j]);
        // if q[j]: win -= b
        ctrl_sub(c, &q[j], &win, &bext);
    }
    c.zero_and_free(wguard);
    c.zero_and_free(bguard);
}

/// Inverse of [`long_division`]: a += q*b, q := 0. PRE: a = (orig a mod b),
/// q = orig a // b. This IS the consuming multiply `|a| += q|b|`.
pub fn long_division_reverse(c: &mut Circuit, a: &[QReg], b: &[QReg], q: &[QReg]) {
    let n = a.len();
    let m = b.len();
    assert_eq!(q.len(), n - m + 1, "q width must be n-m+1");
    let bguard = c.alloc_qreg("ld.bguard");
    let wguard = c.alloc_qreg("ld.wguard");
    let bext: Vec<&QReg> = b.iter().chain(std::iter::once(&bguard)).collect();
    for j in 0..=n - m {
        let mut win: Vec<&QReg> = a[j..(j + m).min(n)].iter().collect();
        if j + m < n {
            win.push(&a[j + m]);
        } else {
            win.push(&wguard);
        }
        // undo: if q[j], win += b ; then uncompute q[j] (X; re-compare).
        ctrl_add(c, &q[j], &win, &bext);
        c.x(&q[j]);
        borrow_compare_refs(c, &win, &bext, &q[j]); // q[j] -> 0
    }
    c.zero_and_free(wguard);
    c.zero_and_free(bguard);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::point_add::trailmix_port::num_bigint::BigUint;
    use rand::Rng;

    fn rd(view: &crate::point_add::trailmix_port::circuit::ContractSimView, reg: &[QReg], shot: usize) -> BigUint {
        let mut x = BigUint::from(0u32);
        for (j, qb) in reg.iter().enumerate() {
            if view.contract_read_bit_shot(qb, shot) {
                x |= BigUint::from(1u32) << j;
            }
        }
        x
    }

    /// long_division then long_division_reverse: a -> a mod b (q = a//b) -> a (q=0).
    #[test]
    fn long_division_roundtrip() {
        let n = 64usize;
        let m = 32usize;
        let mut rng = rand::thread_rng();
        let mut c = Circuit::new();
        c.set_max_qubit_peak(400);
        let a = c.alloc_qreg_bits("a", n);
        let b = c.alloc_qreg_bits("b", m);
        let q = c.alloc_qreg_bits("q", n - m + 1);
        // 64 shots: random a < 2^(n), b in [2^(m-1), 2^m) (nonzero high bit).
        let mut avs = Vec::new();
        let mut bvs = Vec::new();
        for shot in 0..64 {
            let av: BigUint = BigUint::from(rng.gen::<u64>() >> 1); // < 2^63
            let bv: BigUint = (BigUint::from(rng.gen::<u32>()) % (BigUint::from(1u32) << m as u32))
                | (BigUint::from(1u32) << (m as u32 - 1)); // normalized m-bit
            let mut al = av.to_bytes_le();
            al.resize(32, 0);
            c.sim_load_reg_bytes_shot(&a, &al, shot);
            let mut bl = bv.to_bytes_le();
            bl.resize(32, 0);
            c.sim_load_reg_bytes_shot(&b, &bl, shot);
            avs.push(av);
            bvs.push(bv);
        }
        long_division(&mut c, &a, &b, &q);
        {
            let (ar, qr, br, av2, bv2) = (&a, &q, &b, avs.clone(), bvs.clone());
            c.contract_check("ld_div", move |view, shot| {
                let rem = rd(&view, ar, shot);
                let quo = rd(&view, qr, shot);
                let bb = rd(&view, br, shot);
                let (av, bv) = (&av2[shot], &bv2[shot]);
                if bb != *bv {
                    return Err("b changed".into());
                }
                if rem != av % bv {
                    return Err(format!("rem wrong: {rem} != {}%{}", av, bv));
                }
                if quo != av / bv {
                    return Err(format!("quo wrong: {quo} != {}/{}", av, bv));
                }
                Ok(())
            });
        }
        long_division_reverse(&mut c, &a, &b, &q);
        {
            let (ar, qr, av2) = (&a, &q, avs.clone());
            c.contract_check("ld_rev", move |view, shot| {
                if rd(&view, ar, shot) != av2[shot] {
                    return Err("a not restored".into());
                }
                if rd(&view, qr, shot) != BigUint::from(0u32) {
                    return Err("q not cleared".into());
                }
                Ok(())
            });
        }
        c.assert_phase_clean();
        eprintln!(
            "LONG DIVISION roundtrip ok: peak {} q, {} tof",
            c.peak_qubits,
            c.executed_toffoli_shots / 64
        );
        let mut outs = vec![];
        outs.extend(a);
        outs.extend(b);
        outs.extend(q);
        let _ = c.destroy_sim(outs);
    }
}
