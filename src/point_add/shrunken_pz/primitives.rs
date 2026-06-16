use crate::circuit::QubitId;

use super::super::{cuccaro_add_ctrl_lowq, cuccaro_sub_ctrl_lowq, B};

pub(crate) fn borrow_compare_refs(b: &mut B, v: &[QubitId], u: &[QubitId], out: QubitId) {
    let n = v.len();
    assert_eq!(u.len(), n);
    if n == 0 {
        return;
    }
    b.set_phase("shrunken_pz_cmp");
    for &q in v {
        b.x(q);
    }
    let carry = b.alloc_qubit();
    b.cx(v[0], u[0]);
    b.cx(v[0], carry);
    b.ccx(carry, u[0], v[0]);
    for i in 1..n {
        b.cx(v[i], u[i]);
        b.cx(v[i], v[i - 1]);
        b.ccx(v[i - 1], u[i], v[i]);
    }
    b.cx(v[n - 1], out);
    for i in (1..n).rev() {
        b.ccx(v[i - 1], u[i], v[i]);
        b.cx(v[i], v[i - 1]);
        b.cx(v[i], u[i]);
    }
    b.ccx(carry, u[0], v[0]);
    b.cx(v[0], carry);
    b.cx(v[0], u[0]);
    b.release_zeroed(carry);
    for &q in v {
        b.x(q);
    }
}

pub(crate) fn ctrl_add(b: &mut B, g: QubitId, a: &[QubitId], addend: &[QubitId]) {
    assert_eq!(a.len(), addend.len());
    let c_in = b.alloc_qubit();
    let scratch = b.alloc_qubit();
    b.set_phase("shrunken_pz_add");
    cuccaro_add_ctrl_lowq(b, addend, a, g, c_in, scratch);
    b.release_zeroed(scratch);
    b.release_zeroed(c_in);
}

pub(crate) fn ctrl_sub(b: &mut B, g: QubitId, a: &[QubitId], subtrahend: &[QubitId]) {
    assert_eq!(a.len(), subtrahend.len());
    let c_in = b.alloc_qubit();
    let scratch = b.alloc_qubit();
    b.set_phase("shrunken_pz_sub");
    cuccaro_sub_ctrl_lowq(b, subtrahend, a, g, c_in, scratch);
    b.release_zeroed(scratch);
    b.release_zeroed(c_in);
}

pub(crate) fn long_division(b: &mut B, a: &[QubitId], divisor: &[QubitId], q: &[QubitId]) {
    let n = a.len();
    let m = divisor.len();
    assert_eq!(q.len(), n - m + 1, "q width must be n-m+1");
    let bguard = b.alloc_qubit();
    let wguard = b.alloc_qubit();
    let mut bext = divisor.to_vec();
    bext.push(bguard);
    for j in (0..=n - m).rev() {
        let mut win: Vec<QubitId> = a[j..(j + m).min(n)].to_vec();
        if j + m < n {
            win.push(a[j + m]);
        } else {
            win.push(wguard);
        }
        borrow_compare_refs(b, &win, &bext, q[j]);
        b.x(q[j]);
        ctrl_sub(b, q[j], &win, &bext);
    }
    b.release_zeroed(wguard);
    b.release_zeroed(bguard);
}

pub(crate) fn long_division_reverse(b: &mut B, a: &[QubitId], divisor: &[QubitId], q: &[QubitId]) {
    let n = a.len();
    let m = divisor.len();
    assert_eq!(q.len(), n - m + 1, "q width must be n-m+1");
    let bguard = b.alloc_qubit();
    let wguard = b.alloc_qubit();
    let mut bext = divisor.to_vec();
    bext.push(bguard);
    for j in 0..=n - m {
        let mut win: Vec<QubitId> = a[j..(j + m).min(n)].to_vec();
        if j + m < n {
            win.push(a[j + m]);
        } else {
            win.push(wguard);
        }
        ctrl_add(b, q[j], &win, &bext);
        b.x(q[j]);
        borrow_compare_refs(b, &win, &bext, q[j]);
    }
    b.release_zeroed(wguard);
    b.release_zeroed(bguard);
}
