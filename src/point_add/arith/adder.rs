use super::*;

pub(crate) fn bit(c: U256, i: usize) -> bool {

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

    b.free_vec(&carries);
}

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

pub(crate) fn cuccaro_add(b: &mut B, a: &[QubitId], acc: &[QubitId], c_in: QubitId) {
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

    maj(b, c_in, acc[0], a[0]);
    for i in 1..n - 1 {
        maj(b, a[i - 1], acc[i], a[i]);
    }

    b.cx(a[n - 2], acc[n - 1]);
    b.cx(a[n - 1], acc[n - 1]);

    for i in (1..n - 1).rev() {
        uma(b, a[i - 1], acc[i], a[i]);
    }
    uma(b, c_in, acc[0], a[0]);
}

pub(crate) fn cuccaro_sub(b: &mut B, a: &[QubitId], acc: &[QubitId], c_in: QubitId) {
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

    inv_uma(b, c_in, acc[0], a[0]);
    for i in 1..n - 1 {
        inv_uma(b, a[i - 1], acc[i], a[i]);
    }

    b.cx(a[n - 1], acc[n - 1]);
    b.cx(a[n - 2], acc[n - 1]);

    for i in (1..n - 1).rev() {
        inv_maj(b, a[i - 1], acc[i], a[i]);
    }
    inv_maj(b, c_in, acc[0], a[0]);
}

pub(crate) fn cuccaro_add_low_to_ext_clean(
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

    maj(b, c_in, acc_ext[0], a[0]);
    for i in 1..n {
        maj(b, a[i - 1], acc_ext[i], a[i]);
    }

    b.cx(a[n - 1], acc_ext[n]);

    for i in (1..n).rev() {
        uma(b, a[i - 1], acc_ext[i], a[i]);
    }
    uma(b, c_in, acc_ext[0], a[0]);
}

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

    inv_uma(b, c_in, acc_ext[0], a[0]);
    for i in 1..n {
        inv_uma(b, a[i - 1], acc_ext[i], a[i]);
    }

    b.cx(a[n - 1], acc_ext[n]);

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

pub(crate) fn ext_reg(b: &mut B, reg: &[QubitId]) -> (Vec<QubitId>, QubitId) {
    let ovf = b.alloc_qubit();
    let mut r = reg.to_vec();
    r.push(ovf);
    (r, ovf)
}

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

    b.free_vec(&carries);
}

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

pub(crate) fn cuccaro_add_fast_low_to_ext_topclean(
    b: &mut B,
    a: &[QubitId],
    acc_ext: &[QubitId],
    c_in: QubitId,
    clean_top: usize,
) {
    let n = a.len();
    assert_eq!(acc_ext.len(), n + 1);
    if n == 0 {
        b.cx(c_in, acc_ext[0]);
        return;
    }
    let clean_top = clean_top.min(n.saturating_sub(1));
    if clean_top == 0 {
        return cuccaro_add_fast_low_to_ext(b, a, acc_ext, c_in);
    }
    let borrowed = n - clean_top;
    let carries = b.alloc_qubits(borrowed);

    b.cx(a[0], acc_ext[0]);
    b.cx(a[0], c_in);
    b.ccx(c_in, acc_ext[0], carries[0]);
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
    b.cz_if(c_in, acc_ext[0], m0);
    b.cx(a[0], c_in);
    b.cx(c_in, acc_ext[0]);

    b.free_vec(&carries);
}

pub(crate) fn cuccaro_sub_fast_low_to_ext_topclean(
    b: &mut B,
    a: &[QubitId],
    acc_ext: &[QubitId],
    c_in: QubitId,
    clean_top: usize,
) {
    let n = a.len();
    assert_eq!(acc_ext.len(), n + 1);
    if n == 0 {
        b.cx(c_in, acc_ext[0]);
        return;
    }
    let clean_top = clean_top.min(n.saturating_sub(1));
    if clean_top == 0 {
        return cuccaro_sub_fast_low_to_ext(b, a, acc_ext, c_in);
    }
    let borrowed = n - clean_top;
    let carries = b.alloc_qubits(borrowed);

    b.cx(c_in, acc_ext[0]);
    b.cx(a[0], c_in);
    b.ccx(c_in, acc_ext[0], carries[0]);
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
    b.cz_if(c_in, acc_ext[0], m0);
    b.cx(a[0], c_in);
    b.cx(a[0], acc_ext[0]);

    b.free_vec(&carries);
}

pub(crate) fn cuccaro_add_fast_low_to_ext_borrowed_carries_topclean(
    b: &mut B,
    a: &[QubitId],
    acc_ext: &[QubitId],
    c_in: QubitId,
    carries: &[QubitId],
    clean_top: usize,
) {
    let n = a.len();
    assert_eq!(acc_ext.len(), n + 1);
    if n == 0 {
        b.cx(c_in, acc_ext[0]);
        return;
    }
    let clean_top = clean_top.min(n.saturating_sub(1));
    if clean_top == 0 {
        return cuccaro_add_fast_low_to_ext_borrowed_carries(b, a, acc_ext, c_in, carries);
    }
    let borrowed = n - clean_top;
    assert!(carries.len() >= borrowed);

    b.cx(a[0], acc_ext[0]);
    b.cx(a[0], c_in);
    b.ccx(c_in, acc_ext[0], carries[0]);
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
    b.cz_if(c_in, acc_ext[0], m0);
    b.cx(a[0], c_in);
    b.cx(c_in, acc_ext[0]);
}

pub(crate) fn cuccaro_sub_fast_low_to_ext_borrowed_carries_topclean(
    b: &mut B,
    a: &[QubitId],
    acc_ext: &[QubitId],
    c_in: QubitId,
    carries: &[QubitId],
    clean_top: usize,
) {
    let n = a.len();
    assert_eq!(acc_ext.len(), n + 1);
    if n == 0 {
        b.cx(c_in, acc_ext[0]);
        return;
    }
    let clean_top = clean_top.min(n.saturating_sub(1));
    if clean_top == 0 {
        return cuccaro_sub_fast_low_to_ext_borrowed_carries(b, a, acc_ext, c_in, carries);
    }
    let borrowed = n - clean_top;
    assert!(carries.len() >= borrowed);

    b.cx(c_in, acc_ext[0]);
    b.cx(a[0], c_in);
    b.ccx(c_in, acc_ext[0], carries[0]);
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
    b.cz_if(c_in, acc_ext[0], m0);
    b.cx(a[0], c_in);
    b.cx(a[0], acc_ext[0]);
}

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

        b.cx(a[0], acc[0]);
        return;
    }
    assert!(carries.len() >= n - 1);

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

    b.cx(carries[0], a[0]);
    let m0 = b.alloc_bit();
    b.hmr(carries[0], m0);
    b.cz_if(a[0], acc[0], m0);
}

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

        b.cx(a[0], acc[0]);
        return;
    }
    assert!(carries.len() >= n - 1);

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

    b.ccx(x, y, w);
    b.cx(w, x);
    b.cx(w, y);
}

pub(crate) fn inv_uma(b: &mut B, x: QubitId, y: QubitId, w: QubitId) {

    b.cx(x, y);
    b.cx(w, x);
    b.ccx(x, y, w);
}

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

pub(crate) fn cuccaro_add_ctrl_vented(
    b: &mut B, addend: &[QubitId], acc: &[QubitId], ctrl: QubitId, vent_pool: &[QubitId],
) {
    let n = addend.len();
    assert_eq!(n, acc.len());
    if n == 0 { return; }
    if n == 1 { b.ccx(ctrl, addend[0], acc[0]); return; }
    assert!(vent_pool.len() >= n - 1, "vented body needs n-1 borrowed vent lanes");
    for i in 1..n { b.cx(addend[i], acc[i]); }
    for i in (1..n-1).rev() { b.cx(addend[i], addend[i+1]); }
    for i in 0..n-1 {
        let anc = vent_pool[i];
        b.ccx(acc[i], addend[i], anc);
        b.cx(anc, addend[i+1]);
    }
    for i in (0..n-1).rev() {
        b.ccx(ctrl, addend[i+1], acc[i+1]);
        let anc = vent_pool[i];
        b.cx(anc, addend[i+1]);
        let m = b.alloc_bit();
        b.hmr(anc, m);
        b.cz_if(acc[i], addend[i], m);
    }
    for i in 1..n-1 { b.cx(addend[i], addend[i+1]); }
    b.ccx(ctrl, addend[0], acc[0]);
    for i in 1..n { b.cx(addend[i], acc[i]); }
}

pub(crate) fn cuccaro_sub_ctrl_vented(
    b: &mut B, subtrahend: &[QubitId], acc: &[QubitId], ctrl: QubitId, vent_pool: &[QubitId],
) {
    for &q in acc { b.x(q); }
    cuccaro_add_ctrl_vented(b, subtrahend, acc, ctrl, vent_pool);
    for &q in acc { b.x(q); }
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
