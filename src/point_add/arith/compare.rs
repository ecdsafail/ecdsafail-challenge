use super::*;

pub(crate) fn cmp_lt_into_fast(b: &mut B, u: &[QubitId], v: &[QubitId], flag: QubitId) {
    // The vented D1 core uses the slow (no-carries) comparator which
    // saves n peak qubits at cost of ~n CCX per call.
    if kal_vent_modadd_enabled() {
        cmp_lt_into(b, u, v, flag);
        return;
    }
    let n = u.len();
    assert_eq!(n, v.len());
    let c_in = b.alloc_qubit();
    let carries = b.alloc_qubits(n);
    for i in 0..n {
        b.x(u[i]);
    }

    // Forward MAJ sweep with carry ancillae
    b.cx(u[0], v[0]);
    b.cx(u[0], c_in);
    b.ccx(c_in, v[0], carries[0]);
    b.cx(carries[0], u[0]);
    for i in 1..n {
        b.cx(u[i], v[i]);
        b.cx(u[i], u[i - 1]);
        b.ccx(u[i - 1], v[i], carries[i]);
        b.cx(carries[i], u[i]);
    }

    b.cx(u[n - 1], flag);

    // Backward inv_MAJ with measurement
    for i in (1..n).rev() {
        b.cx(carries[i], u[i]);
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        b.cz_if(u[i - 1], v[i], m);
        b.cx(u[i], u[i - 1]);
        b.cx(u[i], v[i]);
    }
    b.cx(carries[0], u[0]);
    let m0 = b.alloc_bit();
    b.hmr(carries[0], m0);
    b.cz_if(c_in, v[0], m0);
    b.cx(u[0], c_in);
    b.cx(u[0], v[0]);

    for i in 0..n {
        b.x(u[i]);
    }
    b.free_vec(&carries);
    b.free(c_in);
}

#[allow(dead_code)] // retained reference/alternative impl; not on active build path
pub(crate) fn cmp_lt_into_fast_with_cin(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    c_in: QubitId,
    flag: QubitId,
) {
    let n = u.len();
    assert_eq!(n, v.len());
    assert!(!u.contains(&c_in));
    assert!(!v.contains(&c_in));
    assert_ne!(c_in, flag);
    assert!(!u.contains(&flag));
    assert!(!v.contains(&flag));
    let carries = b.alloc_qubits(n);
    for i in 0..n {
        b.x(u[i]);
    }

    b.cx(u[0], v[0]);
    b.cx(u[0], c_in);
    b.ccx(c_in, v[0], carries[0]);
    b.cx(carries[0], u[0]);
    for i in 1..n {
        b.cx(u[i], v[i]);
        b.cx(u[i], u[i - 1]);
        b.ccx(u[i - 1], v[i], carries[i]);
        b.cx(carries[i], u[i]);
    }

    b.cx(u[n - 1], flag);

    for i in (1..n).rev() {
        b.cx(carries[i], u[i]);
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        b.cz_if(u[i - 1], v[i], m);
        b.cx(u[i], u[i - 1]);
        b.cx(u[i], v[i]);
    }
    b.cx(carries[0], u[0]);
    let m0 = b.alloc_bit();
    b.hmr(carries[0], m0);
    b.cz_if(c_in, v[0], m0);
    b.cx(u[0], c_in);
    b.cx(u[0], v[0]);

    for i in 0..n {
        b.x(u[i]);
    }
    b.free_vec(&carries);
}

/// Like `cmp_lt_into_fast_with_cin` but the n-wide measured-uncompute carry lane
/// is supplied by the caller as borrowed clean (|0>) qubits (restored clean on
/// exit) instead of being allocated — so the comparator adds no peak qubits.
pub(crate) fn cmp_lt_into_fast_with_cin_borrowed_carries(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    c_in: QubitId,
    flag: QubitId,
    carries: &[QubitId],
) {
    let n = u.len();
    assert_eq!(n, v.len());
    assert!(carries.len() >= n);
    for i in 0..n {
        b.x(u[i]);
    }
    b.cx(u[0], v[0]);
    b.cx(u[0], c_in);
    b.ccx(c_in, v[0], carries[0]);
    b.cx(carries[0], u[0]);
    for i in 1..n {
        b.cx(u[i], v[i]);
        b.cx(u[i], u[i - 1]);
        b.ccx(u[i - 1], v[i], carries[i]);
        b.cx(carries[i], u[i]);
    }
    b.cx(u[n - 1], flag);
    for i in (1..n).rev() {
        b.cx(carries[i], u[i]);
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        b.cz_if(u[i - 1], v[i], m);
        b.cx(u[i], u[i - 1]);
        b.cx(u[i], v[i]);
    }
    b.cx(carries[0], u[0]);
    let m0 = b.alloc_bit();
    b.hmr(carries[0], m0);
    b.cz_if(c_in, v[0], m0);
    b.cx(u[0], c_in);
    b.cx(u[0], v[0]);
    for i in 0..n {
        b.x(u[i]);
    }
}

pub(crate) fn cmp_lt_fast_prefix_window_forward(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    c_in: QubitId,
    carries: &[QubitId],
    ctrl: QubitId,
    targets: &[(QubitId, usize)],
) {
    let n = u.len();
    assert_eq!(n, v.len());
    assert!(n > 0);
    assert!(carries.len() >= n);
    assert!(targets.iter().all(|&(_, p)| (1..=n).contains(&p)));
    assert!(targets.windows(2).all(|w| w[0].1 < w[1].1));

    b.cx(u[0], v[0]);
    b.cx(u[0], c_in);
    b.ccx(c_in, v[0], carries[0]);
    b.cx(carries[0], u[0]);
    let mut next_target = 0usize;
    while next_target < targets.len() && targets[next_target].1 == 1 {
        b.ccx(ctrl, u[0], targets[next_target].0);
        next_target += 1;
    }
    for i in 1..n {
        b.cx(u[i], v[i]);
        b.cx(u[i], u[i - 1]);
        b.ccx(u[i - 1], v[i], carries[i]);
        b.cx(carries[i], u[i]);
        while next_target < targets.len() && targets[next_target].1 == i + 1 {
            b.ccx(ctrl, u[i], targets[next_target].0);
            next_target += 1;
        }
    }
    assert_eq!(next_target, targets.len());
}

pub(crate) fn cmp_lt_fast_prefix_window_inverse(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    c_in: QubitId,
    carries: &[QubitId],
) {
    let n = u.len();
    assert_eq!(n, v.len());
    assert!(n > 0);
    assert!(carries.len() >= n);

    for i in (1..n).rev() {
        b.cx(carries[i], u[i]);
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        b.cz_if(u[i - 1], v[i], m);
        b.cx(u[i], u[i - 1]);
        b.cx(u[i], v[i]);
    }
    b.cx(carries[0], u[0]);
    let m0 = b.alloc_bit();
    b.hmr(carries[0], m0);
    b.cz_if(c_in, v[0], m0);
    b.cx(u[0], c_in);
    b.cx(u[0], v[0]);
}

/// Apply the HMR phase correction for `u < v + c_in` without an additional
/// quantum control. The nonlinear comparator replay executes only when the
/// classical HMR result is one.
pub(crate) fn cmp_lt_phase_conditioned_with_cin_borrowed_carries(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    c_in: QubitId,
    carries: &[QubitId],
    phase: BitId,
) {
    let n = u.len();
    assert_eq!(v.len(), n);
    assert!(n > 0);
    assert!(carries.len() >= n);

    b.push_condition(phase);
    for &q in u {
        b.x(q);
    }
    cmp_lt_fast_prefix_window_forward(b, u, v, c_in, carries, c_in, &[]);
    b.cz(u[n - 1], u[n - 1]);
    cmp_lt_fast_prefix_window_inverse(b, u, v, c_in, carries);
    for &q in u {
        b.x(q);
    }
    b.pop_condition();
}

pub(crate) fn cmp_lt_phase_conditioned(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    phase: BitId,
) {
    let n = u.len();
    assert_eq!(v.len(), n);
    assert!(n > 0);

    let c_in = b.alloc_qubit();
    b.push_condition(phase);
    for &q in u {
        b.x(q);
    }
    let carries = b.alloc_qubits(n);
    cmp_lt_fast_prefix_window_forward(b, u, v, c_in, &carries, c_in, &[]);
    b.cz(u[n - 1], u[n - 1]);
    cmp_lt_fast_prefix_window_inverse(b, u, v, c_in, &carries);
    b.free_vec(&carries);
    for &q in u {
        b.x(q);
    }
    b.pop_condition();
    b.free(c_in);
}


/// Slow (carry-array-free) `flag ^= (u < v + c_in)` comparator. Like
/// `cmp_lt_into` but threads a borrowed carry-IN qubit (left clean on exit)
/// through the bottom MAJ. Peak cost: 0 extra qubits beyond the supplied c_in
/// (the MAJ sweep works in place on `u`). Toffoli ~2n (no measured uncompute),
/// traded against the n-wide carry array the fast variant allocates.
pub(crate) fn cmp_lt_into_with_cin_slow(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    c_in: QubitId,
    flag: QubitId,
) {
    let n = u.len();
    assert_eq!(n, v.len());
    assert!(n > 0);
    for i in 0..n {
        b.x(u[i]);
    }
    maj(b, c_in, v[0], u[0]);
    for i in 1..n {
        maj(b, u[i - 1], v[i], u[i]);
    }
    b.cx(u[n - 1], flag);
    for i in (1..n).rev() {
        inv_maj(b, u[i - 1], v[i], u[i]);
    }
    inv_maj(b, c_in, v[0], u[0]);
    for i in 0..n {
        b.x(u[i]);
    }
}

pub(crate) fn cmp_lt_into(b: &mut B, u: &[QubitId], v: &[QubitId], flag: QubitId) {
    let n = u.len();
    assert_eq!(n, v.len());

    let c_in = b.alloc_qubit();

    // ~u in place (X is free in the metric).
    for i in 0..n {
        b.x(u[i]);
    }

    // Forward MAJ sweep — n MAJs (one more than cuccaro_add, which omits
    // the top one because it doesn't need the carry-out).
    maj(b, c_in, v[0], u[0]);
    for i in 1..n {
        maj(b, u[i - 1], v[i], u[i]);
    }
    // u[n-1] now holds the high carry = (u < v).
    b.cx(u[n - 1], flag);

    // Inverse sweep restores u and v to their (negated u) state.
    for i in (1..n).rev() {
        inv_maj(b, u[i - 1], v[i], u[i]);
    }
    inv_maj(b, c_in, v[0], u[0]);

    // Un-negate u.
    for i in 0..n {
        b.x(u[i]);
    }

    b.free(c_in);
}

