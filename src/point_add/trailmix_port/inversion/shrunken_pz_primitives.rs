//! Reversible UNPACKED PZ modular inversion. Separate registers A,B,|a|,|b| (no
//! cursor), reversible via the PZ cofactor ratio (no spooky pebbling). The whole
//! thing is built from ONE primitive -- restoring long division -- used forward
//! on the gcd pair and (its reverse) as the cofactor multiply.
//!
//! `long_division(A,B,q)`: A := A mod B, q := A // B (q starts |0>). Reversible.
//! `long_division_reverse(A,B,q)`: the inverse -- A := A + q*B, q := 0 (consumes
//! q). This is exactly the consuming multiply `|a| += q|b|` applied to (|a|,|b|).

use crate::point_add::trailmix_port::arith::cuccaro::controlled_add_cuccaro_3n_refs;
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
    not_gate: Option<&QReg>,
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
    let owned_carry = borrowed_carry.is_none().then(|| circ.alloc_qreg("bc.c"));
    let cc = borrowed_carry.unwrap_or_else(|| owned_carry.as_ref().unwrap());
    assert!(!std::ptr::eq(cc, out), "comparator carry aliases output");
    assert!(
        !v.iter().chain(u).any(|&q| std::ptr::eq(q, cc)),
        "comparator carry aliases an operand"
    );
    circ.cx(b[0], a[0]);
    circ.cx(b[0], cc);
    circ.ccx(cc, a[0], b[0]);
    for i in 1..n {
        circ.cx(b[i], a[i]);
        circ.cx(b[i], b[i - 1]);
        circ.ccx(b[i - 1], a[i], b[i]);
    }
    match (active, not_gate) {
        (Some(active), Some(not_gate)) => {
            assert!(!std::ptr::eq(active, not_gate), "comparator gates alias");
            assert!(
                !std::ptr::eq(not_gate, out),
                "comparator NOT gate aliases output"
            );
            assert!(
                !v.iter().chain(u).any(|&q| std::ptr::eq(q, not_gate)),
                "comparator NOT gate aliases a compared operand"
            );
            circ.x(not_gate);
            mcx_dirty_ladder(circ, &[active, not_gate, b[n - 1]], out, &[a[0]]);
            circ.x(not_gate);
        }
        (Some(active), None) => circ.ccx(active, b[n - 1], out),
        (None, None) => circ.cx(b[n - 1], out),
        (None, Some(_)) => panic!("comparator NOT gate requires an active gate"),
    }
    for i in (1..n).rev() {
        circ.ccx(b[i - 1], a[i], b[i]);
        circ.cx(b[i], b[i - 1]);
        circ.cx(b[i], a[i]);
    }
    circ.ccx(cc, a[0], b[0]);
    circ.cx(b[0], cc);
    circ.cx(b[0], a[0]);
    if let Some(cc) = owned_carry {
        circ.zero_and_free(cc);
    }
    for q in v {
        circ.x(q);
    }
    circ.pop_section(&pcmp);
}

/// `out ^= (v < u)` (MAJ cascade + carry capture + un-MAJ), v,u restored. Refs
/// variant of `two_cursor::borrow_compare` (windows here are non-contiguous).
pub(crate) fn borrow_compare_refs(circ: &mut Circuit, v: &[&QReg], u: &[&QReg], out: &QReg) {
    borrow_compare_refs_inner(circ, v, u, None, None, out, None);
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
    borrow_compare_refs_inner(circ, v, u, None, None, out, Some(carry));
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
    assert!(
        !std::ptr::eq(active, out),
        "gated compare control aliases output"
    );
    assert!(
        !v.iter()
            .chain(u)
            .any(|&q| std::ptr::eq(q, active) || std::ptr::eq(q, out)),
        "gated compare control/output aliases an operand"
    );
    borrow_compare_refs_inner(circ, v, u, Some(active), None, out, None);
}

/// `out ^= active AND (v < u)` using a caller-provided clean carry lane.
/// The carry is restored to zero before return.
pub(crate) fn borrow_compare_gated_refs_with_carry(
    circ: &mut Circuit,
    v: &[&QReg],
    u: &[&QReg],
    active: &QReg,
    out: &QReg,
    carry: &QReg,
) {
    assert!(
        !std::ptr::eq(active, carry),
        "comparator carry aliases active"
    );
    assert!(!std::ptr::eq(out, carry), "comparator carry aliases output");
    assert!(
        !v.iter()
            .chain(u)
            .any(|&q| std::ptr::eq(q, active) || std::ptr::eq(q, out) || std::ptr::eq(q, carry)),
        "gated compare control/output/carry aliases an operand"
    );
    borrow_compare_refs_inner(circ, v, u, Some(active), None, out, Some(carry));
}

/// `out ^= active AND NOT(not_gate) AND (v < u)` with a caller-provided clean
/// carry. `not_gate` must be outside the compared slices. All inputs and the
/// dirty ladder lender are restored exactly.
pub(crate) fn borrow_compare_gated_not_refs_with_carry(
    circ: &mut Circuit,
    v: &[&QReg],
    u: &[&QReg],
    active: &QReg,
    not_gate: &QReg,
    out: &QReg,
    carry: &QReg,
) {
    assert!(
        !v.is_empty(),
        "gated-NOT comparator requires a dirty lender"
    );
    assert!(!std::ptr::eq(active, not_gate), "comparator gates alias");
    assert!(
        !std::ptr::eq(active, carry),
        "comparator carry aliases active"
    );
    assert!(
        !std::ptr::eq(not_gate, carry),
        "comparator carry aliases NOT gate"
    );
    assert!(!std::ptr::eq(out, carry), "comparator carry aliases output");
    assert!(
        !v.iter().chain(u).any(|&q| {
            std::ptr::eq(q, active)
                || std::ptr::eq(q, not_gate)
                || std::ptr::eq(q, out)
                || std::ptr::eq(q, carry)
        }),
        "gated-NOT comparator control/output/carry aliases an operand"
    );
    borrow_compare_refs_inner(circ, v, u, Some(active), Some(not_gate), out, Some(carry));
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Q945Row364ComparatorReport {
    pub exhaustive_widths_checked: usize,
    pub unequal_states_checked: usize,
    pub equal_states_checked: usize,
    pub phase_cleanup_states_checked: usize,
    pub ancilla_cleanup_states_checked: usize,
    pub carry_restoration_states_checked: usize,
    pub exact_width: usize,
    pub exact_width_x_ops: usize,
    pub exact_width_cx_ops: usize,
    pub exact_width_ccx_ops: usize,
    pub exact_width_emitted_ops: usize,
    pub exact_width_external_qubits: usize,
    pub exact_width_extra_qubits: usize,
}

/// Exhaust the finite control/output truth table and all operand pairs for
/// widths one through four, then bind that parametric primitive to the exact
/// row-364 lower-80 emission. The 80-bit check is structural: no new qubit may
/// be allocated and only classical reversible gates are permitted.
#[doc(hidden)]
pub fn exhaustive_q945_row364_compare_check() -> Q945Row364ComparatorReport {
    use crate::circuit::{OperationType, QubitId};
    use crate::point_add::B;
    use crate::sim::Simulator;
    use sha3::{
        digest::{ExtendableOutput, Update},
        Shake128,
    };

    struct Harness {
        builder: B,
        a: Vec<u32>,
        b: Vec<u32>,
        active: u32,
        not_gate: u32,
        out: u32,
        carry: u32,
        external: Vec<u32>,
    }

    fn ids(reg: &[QReg]) -> Vec<u32> {
        reg.iter().map(QReg::id).collect()
    }

    fn build(width: usize) -> Harness {
        let mut c = Circuit::new();
        let a = c.alloc_qreg_bits("q945-row364.a", width);
        let b = c.alloc_qreg_bits("q945-row364.b", width);
        let active = c.alloc_qreg("q945-row364.active");
        let not_gate = c.alloc_qreg("q945-row364.a80");
        let out = c.alloc_qreg("q945-row364.out");
        let carry = c.alloc_qreg("q945-row364.b80");
        let ar: Vec<&QReg> = a.iter().collect();
        let br: Vec<&QReg> = b.iter().collect();
        borrow_compare_gated_not_refs_with_carry(
            &mut c, &ar, &br, &active, &not_gate, &out, &carry,
        );
        let a_ids = ids(&a);
        let b_ids = ids(&b);
        let external: Vec<u32> = a_ids
            .iter()
            .chain(&b_ids)
            .copied()
            .chain([active.id(), not_gate.id(), out.id(), carry.id()])
            .collect();
        Harness {
            builder: c.into_builder(),
            a: a_ids,
            b: b_ids,
            active: active.id(),
            not_gate: not_gate.id(),
            out: out.id(),
            carry: carry.id(),
            external,
        }
    }

    let mut unequal_states_checked = 0usize;
    let mut equal_states_checked = 0usize;
    let mut cleanup_states_checked = 0usize;
    for width in 1..=4usize {
        let harness = build(width);
        let limit = 1usize << width;
        let states: Vec<_> = (0..limit)
            .flat_map(|a| {
                (0..limit).flat_map(move |b| {
                    [false, true].into_iter().flat_map(move |active| {
                        [false, true].into_iter().flat_map(move |not_gate| {
                            [false, true]
                                .into_iter()
                                .map(move |out| (a, b, active, not_gate, out))
                        })
                    })
                })
            })
            .collect();
        for (batch, chunk) in states.chunks(64).enumerate() {
            let mut seed = Shake128::default();
            seed.update(b"q945-row364-exhaustive");
            seed.update(&(width as u64).to_le_bytes());
            seed.update(&(batch as u64).to_le_bytes());
            let mut xof = seed.finalize_xof();
            let mut sim = Simulator::new(
                harness.builder.next_qubit as usize,
                harness.builder.next_bit as usize,
                &mut xof,
            );
            let mut expected_a = vec![0u64; width];
            let mut expected_b = vec![0u64; width];
            let mut expected_active = 0u64;
            let mut expected_not_gate = 0u64;
            let mut expected_out = 0u64;
            for (shot, &(a, b, active, not_gate, out)) in chunk.iter().enumerate() {
                let mask = 1u64 << shot;
                for bit in 0..width {
                    if (a >> bit) & 1 == 1 {
                        *sim.qubit_mut(QubitId(u64::from(harness.a[bit]))) |= mask;
                        expected_a[bit] |= mask;
                    }
                    if (b >> bit) & 1 == 1 {
                        *sim.qubit_mut(QubitId(u64::from(harness.b[bit]))) |= mask;
                        expected_b[bit] |= mask;
                    }
                }
                if active {
                    *sim.qubit_mut(QubitId(u64::from(harness.active))) |= mask;
                    expected_active |= mask;
                }
                if not_gate {
                    *sim.qubit_mut(QubitId(u64::from(harness.not_gate))) |= mask;
                    expected_not_gate |= mask;
                }
                if out {
                    *sim.qubit_mut(QubitId(u64::from(harness.out))) |= mask;
                }
                if out ^ (active && !not_gate && a < b) {
                    expected_out |= mask;
                }
                if a == b {
                    equal_states_checked += 1;
                } else {
                    unequal_states_checked += 1;
                }
            }
            sim.apply_iter(harness.builder.ops.iter());
            for (bit, &id) in harness.a.iter().enumerate() {
                assert_eq!(sim.qubit(QubitId(u64::from(id))), expected_a[bit]);
            }
            for (bit, &id) in harness.b.iter().enumerate() {
                assert_eq!(sim.qubit(QubitId(u64::from(id))), expected_b[bit]);
            }
            assert_eq!(
                sim.qubit(QubitId(u64::from(harness.active))),
                expected_active
            );
            assert_eq!(
                sim.qubit(QubitId(u64::from(harness.not_gate))),
                expected_not_gate
            );
            assert_eq!(sim.qubit(QubitId(u64::from(harness.out))), expected_out);
            assert_eq!(sim.qubit(QubitId(u64::from(harness.carry))), 0);
            assert_eq!(sim.phase, 0, "Q945 row-364 comparator left phase garbage");
            for id in 0..harness.builder.next_qubit {
                if !harness.external.contains(&id) {
                    assert_eq!(sim.qubit(QubitId(u64::from(id))), 0);
                }
            }
            cleanup_states_checked += chunk.len();
        }
    }

    let exact_width = 80usize;
    let exact = build(exact_width);
    let exact_width_external_qubits = 2 * exact_width + 4;
    assert_eq!(exact.external.len(), exact_width_external_qubits);
    assert_eq!(
        exact.builder.next_qubit as usize,
        exact_width_external_qubits
    );
    assert_eq!(
        exact.builder.active_qubits as usize,
        exact_width_external_qubits
    );
    assert_eq!(
        exact.builder.peak_qubits as usize,
        exact_width_external_qubits
    );
    let exact_width_x_ops = exact
        .builder
        .ops
        .iter()
        .filter(|op| op.kind == OperationType::X)
        .count();
    let exact_width_cx_ops = exact
        .builder
        .ops
        .iter()
        .filter(|op| op.kind == OperationType::CX)
        .count();
    let exact_width_ccx_ops = exact
        .builder
        .ops
        .iter()
        .filter(|op| op.kind == OperationType::CCX)
        .count();
    let exact_width_emitted_ops = exact.builder.ops.len();
    assert_eq!(exact_width_x_ops, 2 * exact_width + 2);
    assert_eq!(exact_width_cx_ops, 4 * exact_width);
    assert_eq!(exact_width_ccx_ops, 2 * exact_width + 4);
    assert_eq!(
        exact_width_emitted_ops,
        exact_width_x_ops + exact_width_cx_ops + exact_width_ccx_ops
    );
    assert_eq!(exact_width_emitted_ops, 8 * exact_width + 6);
    assert_eq!(unequal_states_checked, 2_480);
    assert_eq!(equal_states_checked, 240);
    assert_eq!(cleanup_states_checked, 2_720);

    Q945Row364ComparatorReport {
        exhaustive_widths_checked: 4,
        unequal_states_checked,
        equal_states_checked,
        phase_cleanup_states_checked: cleanup_states_checked,
        ancilla_cleanup_states_checked: cleanup_states_checked,
        carry_restoration_states_checked: cleanup_states_checked,
        exact_width,
        exact_width_x_ops,
        exact_width_cx_ops,
        exact_width_ccx_ops,
        exact_width_emitted_ops,
        exact_width_external_qubits,
        exact_width_extra_qubits: exact.builder.peak_qubits as usize - exact_width_external_qubits,
    }
}

/// a += b (mod 2^len) gated on `g`. Plain controlled Cuccaro (3n).
pub(crate) fn ctrl_add(c: &mut Circuit, g: &QReg, a: &[&QReg], b: &[&QReg]) {
    let prev = c.push_section("p.add");
    controlled_add_cuccaro_3n_refs(c, g, a, b);
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

/// `a += g*b (mod 2^n)` without allocating a carry qubit.
///
/// Each addend bit controls an increment of the corresponding suffix of `a`.
/// The increment is emitted from high to low so its controls see the
/// pre-increment lower suffix. Multi-controlled X gates borrow the other bits
/// of `b` as dirty lenders and restore them exactly. This route is intended for
/// the five-bit hybrid-CLZ transcript update, where saving Cuccaro's one clean
/// carry removes the global peak qubit at modest gate cost.
pub(crate) fn ctrl_add_dirty_lenders(c: &mut Circuit, g: &QReg, a: &[&QReg], b: &[&QReg]) {
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
pub(crate) fn ctrl_sub_dirty_lenders(c: &mut Circuit, g: &QReg, a: &[&QReg], b: &[&QReg]) {
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

    fn rd(
        view: &crate::point_add::trailmix_port::circuit::ContractSimView,
        reg: &[QReg],
        shot: usize,
    ) -> BigUint {
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
