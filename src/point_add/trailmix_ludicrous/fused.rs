//! Fused double + controlled-double: `y := y * 2 * (1 + s2) mod q` with one
//! combined `(e+2d)*f` reduction instead of two separate `+f` const-adds.
//! AND-uncompute discharge is hmr+cz_if_bit; clear_and is inlined. This file
//! provides the full-clean fold (nv = L-1) plus the chunked and borrowed-dirty
//! variants (lower nv).

use super::arith::{F_SECP256K1, LSBS};
use super::{B, BExt};
use crate::circuit::{BitId, QubitId};

/// secp256k1 `(e+2d)*f` combined-fold addend control per low-bit position `p`
/// (encodes the bit pattern of `f` and `2f`, f = 2^32+977). 0 = None.
fn fold_ctl(p: usize) -> u8 {
    match p {
        0 | 4 | 6 | 32 => 1, // E
        1 | 5 | 33 => 2,     // D
        7 => 3,              // Xor (e^d)
        8 | 9 => 4,          // Or  (e|d)
        10 => 5,             // DnotE (d&~e)
        11 => 6,             // And (e&d)
        _ => 0,
    }
}

/// MBU AND-uncompute (HMR + conditional-CZ): `t` holds `a AND b` -> |0>, phase clean.
fn clear_and(circ: &mut B, t: &QubitId, a: &QubitId, b: &QubitId) {
    let bit = circ.alloc_bit();
    circ.hmr(*t, bit);
    circ.cz_if_bit(*a, *b, bit);
}

/// Carry-propagate `c` into the pure-propagation tail `y[..]` via a cascade of
/// prefix-controlled increments (`mcx_clean_k`, log* ancillae): the clean-tail
/// fold's tail [nv, L).
fn add_carry_into_tail_prefix(circ: &mut B, y: &[QubitId], c: &QubitId) {
    let t = y.len();
    for k in (1..t).rev() {
        let mut ctrls: Vec<&QubitId> = Vec::with_capacity(k + 1);
        ctrls.push(c);
        ctrls.extend(y[..k].iter());
        super::mcx::mcx_clean_k(circ, &ctrls, &y[k]);
    }
    circ.cx(*c, y[0]);
}

/// `tail_from = None` => full clean (top bit L-1 folded specially); `Some(nv)` =>
/// clean-tail: carry chain over [0, nv), then a prefix-controlled increment of the
/// pure-propagation tail [nv, L).
fn add_mf_fold_clean(circ: &mut B, e: &QubitId, d: &QubitId, y: &[QubitId]) {
    add_mf_fold_clean_tail(circ, e, d, y, None);
}

fn add_mf_fold_clean_tail(circ: &mut B, e: &QubitId, d: &QubitId, y: &[QubitId], tail_from: Option<usize>) {
    let l = y.len();
    assert!(l >= 2, "fold needs L >= 2");
    let loop_end = tail_from.unwrap_or(l - 1);
    const LAST_DERIVED: usize = 9;
    const LAST_AND: usize = 11;

    // Derived controls: cc = e&d, dne = d&~e, sxor = e^d, sor = e|d.
    let mut cc = Some(circ.alloc_qubit());
    circ.ccx(*e, *d, *cc.as_ref().unwrap());
    let mut dne = Some(circ.alloc_qubit());
    circ.x(*e);
    circ.ccx(*e, *d, *dne.as_ref().unwrap());
    circ.x(*e);
    let mut sxor = Some(circ.alloc_qubit());
    circ.cx(*e, *sxor.as_ref().unwrap());
    circ.cx(*d, *sxor.as_ref().unwrap());
    let mut sor = Some(circ.alloc_qubit());
    circ.cx(*sxor.as_ref().unwrap(), *sor.as_ref().unwrap());
    circ.cx(*cc.as_ref().unwrap(), *sor.as_ref().unwrap());

    // Resolve position -> addend control qubit (None when A[p]==0).
    fn fc<'a>(p: usize, e: &'a QubitId, d: &'a QubitId, cc: Option<&'a QubitId>, dne: Option<&'a QubitId>, sx: Option<&'a QubitId>, so: Option<&'a QubitId>) -> Option<&'a QubitId> {
        match fold_ctl(p) {
            1 => Some(e),
            2 => Some(d),
            3 => sx,
            4 => so,
            5 => dne,
            6 => cc,
            _ => None,
        }
    }

    // Forward Gidney-clean carry chain with inline addend sums.
    let mut cy: Vec<Option<QubitId>> = Vec::with_capacity(l - 1);
    let c1 = circ.alloc_qubit();
    if let Some(a0) = fc(0, e, d, cc.as_ref(), dne.as_ref(), sxor.as_ref(), sor.as_ref()) {
        circ.ccx(*a0, y[0], c1);
        circ.cx(*a0, y[0]);
    }
    cy.push(Some(c1));
    for i in 1..loop_end {
        let next = circ.alloc_qubit();
        {
            let ci = cy[i - 1].as_ref().unwrap();
            circ.cx(*ci, y[i]); // y[i] ^= carry_i
            if let Some(ai) = fc(i, e, d, cc.as_ref(), dne.as_ref(), sxor.as_ref(), sor.as_ref()) {
                circ.cx(*ai, *ci); // carry_i ^= addend_i
            }
            circ.ccx(y[i], *ci, next);
            if let Some(ai) = fc(i, e, d, cc.as_ref(), dne.as_ref(), sxor.as_ref(), sor.as_ref()) {
                circ.cx(*ai, *ci);
            }
            circ.cx(*ci, next); // carry_{i+1}
            if let Some(ai) = fc(i, e, d, cc.as_ref(), dne.as_ref(), sxor.as_ref(), sor.as_ref()) {
                circ.cx(*ai, y[i]); // inline sum
            }
        }
        cy.push(Some(next));
        if i == LAST_DERIVED {
            let so = sor.take().unwrap();
            circ.cx(*sxor.as_ref().unwrap(), so);
            circ.cx(*cc.as_ref().unwrap(), so);
            circ.zero_and_free(so);
            let sx = sxor.take().unwrap();
            circ.cx(*e, sx);
            circ.cx(*d, sx);
            circ.zero_and_free(sx);
        }
        if i == LAST_AND {
            let dn = dne.take().unwrap();
            circ.x(*e);
            clear_and(circ, &dn, e, d);
            circ.x(*e);
            circ.zero_and_free(dn);
            let c = cc.take().unwrap();
            clear_and(circ, &c, e, d);
            circ.zero_and_free(c);
        }
    }
    match tail_from {
        None => {
            // Top bit: y[L-1] += addend_{L-1} + cy_{L-1}.
            if let Some(at) = fc(l - 1, e, d, cc.as_ref(), dne.as_ref(), sxor.as_ref(), sor.as_ref()) {
                circ.cx(*at, y[l - 1]);
            }
            circ.cx(*cy[l - 2].as_ref().unwrap(), y[l - 1]);
        }
        Some(nv) => {
            // Pure-propagation tail [nv, L): y[nv..] += cy[nv-1]. cc/dne/sxor/sor are
            // all freed (LAST_AND/LAST_DERIVED < nv), so the mcx ancillae sit on top of
            // just the nv carries -> tail peak nv + log*(L-nv).
            add_carry_into_tail_prefix(circ, &y[nv..], cy[nv - 1].as_ref().unwrap());
        }
    }

    // Reverse: rebuild controls, AND-uncompute (hmr+cz) each carry.
    for i in (1..loop_end).rev() {
        if i == LAST_AND {
            let c = circ.alloc_qubit();
            circ.ccx(*e, *d, c);
            cc = Some(c);
            let dn = circ.alloc_qubit();
            circ.x(*e);
            circ.ccx(*e, *d, dn);
            circ.x(*e);
            dne = Some(dn);
        }
        if i == LAST_DERIVED {
            let sx = circ.alloc_qubit();
            circ.cx(*e, sx);
            circ.cx(*d, sx);
            let so = circ.alloc_qubit();
            circ.cx(sx, so);
            circ.cx(*cc.as_ref().unwrap(), so);
            sxor = Some(sx);
            sor = Some(so);
        }
        if let Some(ai) = fc(i, e, d, cc.as_ref(), dne.as_ref(), sxor.as_ref(), sor.as_ref()) {
            circ.cx(*ai, y[i]);
        }
        let next = cy[i].take().unwrap();
        let ci = cy[i - 1].take().unwrap();
        circ.cx(ci, next);
        if let Some(ai) = fc(i, e, d, cc.as_ref(), dne.as_ref(), sxor.as_ref(), sor.as_ref()) {
            circ.cx(*ai, ci);
        }
        // erase next: hmr + cz_if_bit(y[i], ci).
        let bit = circ.alloc_bit();
        circ.hmr(next, bit);
        circ.zero_and_free(next);
        circ.cz_if_bit(y[i], ci, bit);
        if let Some(ai) = fc(i, e, d, cc.as_ref(), dne.as_ref(), sxor.as_ref(), sor.as_ref()) {
            circ.cx(*ai, ci);
            circ.cx(*ai, y[i]);
        }
        cy[i - 1] = Some(ci);
    }
    // Reverse bit 0.
    let cy1 = cy[0].take().unwrap();
    if let Some(a0) = fc(0, e, d, cc.as_ref(), dne.as_ref(), sxor.as_ref(), sor.as_ref()) {
        circ.cx(*a0, y[0]);
        let bit = circ.alloc_bit();
        circ.hmr(cy1, bit);
        circ.zero_and_free(cy1);
        circ.cz_if_bit(y[0], *a0, bit);
        circ.cx(*a0, y[0]);
    } else {
        circ.zero_and_free(cy1);
    }

    // Uncompute the rebuilt derived controls.
    let sx = sxor.take().unwrap();
    let so = sor.take().unwrap();
    let cc = cc.take().unwrap();
    let dne = dne.take().unwrap();
    circ.x(*e);
    clear_and(circ, &dne, e, d);
    circ.x(*e);
    circ.zero_and_free(dne);
    circ.cx(sx, so);
    circ.cx(cc, so);
    circ.zero_and_free(so);
    circ.cx(*e, sx);
    circ.cx(*d, sx);
    circ.zero_and_free(sx);
    clear_and(circ, &cc, e, d);
    circ.zero_and_free(cc);
}

// ============================================================================
// CHUNKED fold: ceil(L/s_chunk) chunks, peak ~ s_chunk + L/s_chunk live carries
// (one chunk's internal carries plus the held boundary carries).
// ============================================================================

/// Build the 4 derived controls (e&d, e^d, e|d, d&~e).
fn build_fold_controls(circ: &mut B, e: &QubitId, d: &QubitId) -> (QubitId, QubitId, QubitId, QubitId) {
    let cc = circ.alloc_qubit();
    circ.ccx(*e, *d, cc);
    let sxor = circ.alloc_qubit();
    circ.cx(*e, sxor);
    circ.cx(*d, sxor);
    let sor = circ.alloc_qubit();
    circ.cx(sxor, sor);
    circ.cx(cc, sor);
    let dne = circ.alloc_qubit();
    circ.x(*e);
    circ.ccx(*e, *d, dne);
    circ.x(*e);
    (cc, sxor, sor, dne)
}

fn uncompute_fold_controls(circ: &mut B, e: &QubitId, d: &QubitId, cc: QubitId, sxor: QubitId, sor: QubitId, dne: QubitId) {
    circ.x(*e);
    clear_and(circ, &dne, e, d);
    circ.x(*e);
    circ.zero_and_free(dne);
    circ.cx(sxor, sor);
    circ.cx(cc, sor);
    circ.zero_and_free(sor);
    circ.cx(*e, sxor);
    circ.cx(*d, sxor);
    circ.zero_and_free(sxor);
    clear_and(circ, &cc, e, d);
    circ.zero_and_free(cc);
}

/// Position -> addend control qubit map.
fn fold_ctl_map(e: QubitId, d: QubitId, cc: QubitId, sxor: QubitId, sor: QubitId, dne: QubitId, l: usize) -> Vec<Option<QubitId>> {
    (0..l).map(|p| match fold_ctl(p) { 1 => Some(e), 2 => Some(d), 3 => Some(sxor), 4 => Some(sor), 5 => Some(dne), 6 => Some(cc), _ => None }).collect()
}

/// One chunk's clean add of the fold addend (ctl) into y, threaded cin/cout.
fn fold_chunk_clean(circ: &mut B, ctl: &[Option<QubitId>], y: &[QubitId], cin: Option<&QubitId>, cout: &QubitId) {
    let s = y.len();
    if s == 0 {
        if let Some(c) = cin { circ.cx(*c, *cout); }
        return;
    }
    let mut cy: Vec<Option<QubitId>> = (0..s - 1).map(|_| Some(circ.alloc_qubit())).collect();
    for i in 0..s {
        let on = ctl[i].as_ref();
        if i == 0 {
            let dst: QubitId = if s == 1 { *cout } else { *cy[0].as_ref().unwrap() };
            match cin {
                Some(c) => {
                    circ.cx(*c, y[0]);
                    if let Some(a) = on { circ.cx(*a, *c); }
                    circ.ccx(y[0], *c, dst);
                    if let Some(a) = on { circ.cx(*a, *c); }
                    circ.cx(*c, dst);
                }
                None => { if let Some(a) = on { circ.ccx(*a, y[0], dst); } }
            }
        } else {
            let ci: QubitId = *cy[i - 1].as_ref().unwrap();
            let dst: QubitId = if i == s - 1 { *cout } else { *cy[i].as_ref().unwrap() };
            circ.cx(ci, y[i]);
            if let Some(a) = on { circ.cx(*a, ci); }
            circ.ccx(y[i], ci, dst);
            if let Some(a) = on { circ.cx(*a, ci); }
            circ.cx(ci, dst);
        }
    }
    for i in 0..s {
        if let Some(a) = ctl[i].as_ref() { circ.cx(*a, y[i]); }
    }
    for i in (0..s - 1).rev() {
        let on = ctl[i].as_ref();
        if let Some(a) = on { circ.cx(*a, y[i]); }
        let next = cy[i].take().unwrap();
        if i == 0 {
            match cin {
                Some(c) => {
                    circ.cx(*c, next);
                    if let Some(a) = on { circ.cx(*a, *c); }
                    let bit = circ.alloc_bit();
                    circ.hmr(next, bit); circ.zero_and_free(next);
                    circ.cz_if_bit(y[0], *c, bit);
                    if let Some(a) = on { circ.cx(*a, *c); circ.cx(*a, y[0]); }
                }
                None => {
                    let bit = circ.alloc_bit();
                    circ.hmr(next, bit); circ.zero_and_free(next);
                    if let Some(a) = on { circ.cz_if_bit(y[0], *a, bit); }
                    if let Some(a) = on { circ.cx(*a, y[0]); }
                }
            }
        } else {
            let ci: QubitId = *cy[i - 1].as_ref().unwrap();
            circ.cx(ci, next);
            if let Some(a) = on { circ.cx(*a, ci); }
            let bit = circ.alloc_bit();
            circ.hmr(next, bit); circ.zero_and_free(next);
            circ.cz_if_bit(y[i], ci, bit);
            if let Some(a) = on { circ.cx(*a, ci); circ.cx(*a, y[i]); }
        }
    }
}

/// Gated-erase a boundary carry: materialize the addend into a temp, run the
/// uncontrolled gated erase on (y, temp), un-materialize.
fn fold_boundary_erase(circ: &mut B, ctl: &[Option<QubitId>], y: &[QubitId], cin: &QubitId, carry: QubitId) {
    let s = y.len();
    let temp: Vec<QubitId> = (0..s).map(|_| circ.alloc_qubit()).collect();
    for (i, c) in ctl.iter().enumerate() {
        if let Some(a) = c { circ.cx(*a, temp[i]); }
    }
    super::arith::erase_carry_gated_opt(circ, None, y, &temp, cin, &carry, None);
    circ.zero_and_free(carry);
    for (i, c) in ctl.iter().enumerate() {
        if let Some(a) = c { circ.cx(*a, temp[i]); }
    }
    for q in temp { circ.zero_and_free(q); }
}

fn add_mf_fold_chunked(circ: &mut B, e: &QubitId, d: &QubitId, y: &[QubitId], s_chunk: usize) {
    let l = y.len();
    let (cc, sxor, sor, dne) = build_fold_controls(circ, e, d);
    let ctl = fold_ctl_map(*e, *d, cc, sxor, sor, dne, l);
    let cin0 = circ.alloc_qubit();
    let nch = l.div_ceil(s_chunk);
    let mut boundary: Vec<QubitId> = Vec::with_capacity(nch);
    for j in 0..nch {
        let lo = j * s_chunk;
        let hi = ((j + 1) * s_chunk).min(l);
        let cout = circ.alloc_qubit();
        let cin = if j == 0 { cin0 } else { boundary[j - 1] };
        fold_chunk_clean(circ, &ctl[lo..hi], &y[lo..hi], Some(&cin), &cout);
        boundary.push(cout);
    }
    for j in (0..nch).rev() {
        let lo = j * s_chunk;
        let hi = ((j + 1) * s_chunk).min(l);
        let bnd = boundary.pop().expect("boundary present");
        let cin = if j == 0 { cin0 } else { boundary[j - 1] };
        fold_boundary_erase(circ, &ctl[lo..hi], &y[lo..hi], &cin, bnd);
    }
    circ.zero_and_free(cin0);
    uncompute_fold_controls(circ, e, d, cc, sxor, sor, dne);
}

// ============================================================================
// GRADUAL fold (borrowed-dirty tail).
// On-demand derived controls (one live at a time) + a clean window [0,nv) handing
// a carry to a dirty body over the high bits, whose carries are stored in borrowed
// bits and discharged by measure-and-correct: hmr to a bit + z_if_bit(t,bit)
// sandwiching the carry recompute = Z^(bit . carry).
// ============================================================================

enum OnCtl {
    None,
    E,
    D,
    Owned(QubitId),
}

/// Self-inverse derived-control build/clear (involution): q ^= ctl(e,d).
fn on_ctl_apply(circ: &mut B, e: &QubitId, d: &QubitId, k: u8, q: &QubitId) {
    match k {
        3 => {
            circ.cx(*e, *q);
            circ.cx(*d, *q);
        }
        4 => {
            circ.x(*e);
            circ.x(*d);
            circ.ccx(*e, *d, *q);
            circ.x(*q);
            circ.x(*e);
            circ.x(*d);
        }
        5 => {
            circ.x(*e);
            circ.ccx(*e, *d, *q);
            circ.x(*e);
        }
        6 => circ.ccx(*e, *d, *q),
        _ => {}
    }
}

fn on_ctl(circ: &mut B, e: &QubitId, d: &QubitId, p: usize) -> OnCtl {
    match fold_ctl(p) {
        1 => OnCtl::E,
        2 => OnCtl::D,
        k @ (3 | 4 | 5 | 6) => {
            let q = circ.alloc_qubit();
            on_ctl_apply(circ, e, d, k, &q);
            OnCtl::Owned(q)
        }
        _ => OnCtl::None,
    }
}

fn on_ctl_ref(c: &OnCtl, e: &QubitId, d: &QubitId) -> Option<QubitId> {
    match c {
        OnCtl::None => None,
        OnCtl::E => Some(*e),
        OnCtl::D => Some(*d),
        OnCtl::Owned(q) => Some(*q),
    }
}

fn on_ctl_free(circ: &mut B, e: &QubitId, d: &QubitId, p: usize, c: OnCtl) {
    if let OnCtl::Owned(q) = c {
        on_ctl_apply(circ, e, d, fold_ctl(p), &q);
        circ.zero_and_free(q);
    }
}

/// Recompute the L-1 carries of `y + A` and XOR them into `out` (the borrowed dirty
/// bits), restoring them. Per-position controls built on-demand.
fn xor_carries_perpos(circ: &mut B, e: &QubitId, d: &QubitId, base: usize, y: &[QubitId], out: &[QubitId], carry_in: Option<&QubitId>) {
    let n = y.len();
    fn ccx_cond(circ: &mut B, aq: Option<&QubitId>, c1: &QubitId, c2: &QubitId, t: &QubitId, g0: bool, g1: bool) {
        if let Some(a) = aq {
            if g0 {
                circ.cx(*a, *c1);
            }
            if g1 {
                circ.cx(*a, *c2);
            }
        }
        circ.ccx(*c1, *c2, *t);
        if let Some(a) = aq {
            if g0 {
                circ.cx(*a, *c1);
            }
            if g1 {
                circ.cx(*a, *c2);
            }
        }
    }
    for i in (1..n - 1).rev() {
        let c = on_ctl(circ, e, d, base + i);
        let aq = on_ctl_ref(&c, e, d);
        let g0 = aq.is_some();
        ccx_cond(circ, aq.as_ref(), &y[i], &out[i - 1], &out[i], g0, false);
        on_ctl_free(circ, e, d, base + i, c);
    }
    for i in 0..n - 1 {
        let c = on_ctl(circ, e, d, base + i);
        if let Some(a) = on_ctl_ref(&c, e, d) {
            circ.cx(a, out[i]);
        }
        on_ctl_free(circ, e, d, base + i, c);
    }
    {
        let c = on_ctl(circ, e, d, base);
        let aq = on_ctl_ref(&c, e, d);
        let g = aq.is_some();
        match carry_in {
            Some(cy) => ccx_cond(circ, aq.as_ref(), cy, &y[0], &out[0], g, g),
            None => {
                let cin = circ.alloc_qubit();
                ccx_cond(circ, aq.as_ref(), &cin, &y[0], &out[0], g, g);
                circ.zero_and_free(cin);
            }
        }
        on_ctl_free(circ, e, d, base, c);
    }
    for i in 1..n - 1 {
        let c = on_ctl(circ, e, d, base + i);
        let aq = on_ctl_ref(&c, e, d);
        let gi = aq.is_some();
        ccx_cond(circ, aq.as_ref(), &y[i], &out[i - 1], &out[i], gi, gi);
        on_ctl_free(circ, e, d, base + i, c);
    }
}

/// Borrowed-dirty carry-chain body. `carry_in` read-only (caller owns); `None` =>
/// carry-in 0. `dirty` (>= l-1 bits) are borrowed real-data bits used as transient
/// carry storage and restored.
fn dirty_body(circ: &mut B, e: &QubitId, d: &QubitId, base: usize, y: &[QubitId], dirty: &[QubitId], carry_in: Option<&QubitId>) {
    let l = y.len();
    assert!(l >= 2);
    assert!(dirty.len() >= l - 1, "need L-1 borrowed dirty bits");
    let mut cin_owned = if carry_in.is_none() { Some(circ.alloc_qubit()) } else { None };
    let mut bits: Vec<BitId> = Vec::with_capacity(l - 1); // bits[i] = measured cy_{i+1}
    let mut prev_new: Option<QubitId> = None;
    for i in 0..l - 1 {
        let new = circ.alloc_qubit();
        let anc = circ.alloc_qubit();
        let ctlh = on_ctl(circ, e, d, base + i);
        {
            let cyi: QubitId = if i == 0 {
                carry_in.copied().unwrap_or_else(|| *cin_owned.as_ref().unwrap())
            } else {
                *prev_new.as_ref().unwrap()
            };
            if let Some(ai) = on_ctl_ref(&ctlh, e, d) {
                circ.cx(ai, anc);
            }
            circ.cx(cyi, anc);
            circ.cx(cyi, y[i]);
            circ.ccx(y[i], anc, new);
            circ.cx(cyi, new); // new = carry_{i+1}
            circ.cx(new, dirty[i]); // store carry copy in borrowed bit
            circ.cx(cyi, anc);
            if let Some(ai) = on_ctl_ref(&ctlh, e, d) {
                circ.cx(ai, anc);
                circ.cx(ai, y[i]); // y[i] = sum_i
            }
        }
        on_ctl_free(circ, e, d, base + i, ctlh);
        circ.zero_and_free(anc);
        if i == 0 {
            if let Some(c) = cin_owned.take() {
                circ.zero_and_free(c);
            }
        } else {
            let b = circ.alloc_bit();
            circ.hmr(*prev_new.as_ref().unwrap(), b);
            circ.zero_and_free(prev_new.take().unwrap());
            bits.push(b);
        }
        prev_new = Some(new);
    }
    let cy_top = prev_new.take().unwrap(); // cy_{l-1}
    {
        let topc = on_ctl(circ, e, d, base + l - 1);
        if let Some(at) = on_ctl_ref(&topc, e, d) {
            circ.cx(at, y[l - 1]);
        }
        on_ctl_free(circ, e, d, base + l - 1, topc);
    }
    circ.cx(cy_top, y[l - 1]);
    let b = circ.alloc_bit();
    circ.hmr(cy_top, b);
    circ.zero_and_free(cy_top);
    bits.push(b);

    // discharge: dirty[i] currently = orig ^ cy_{i+1}; recompute restores it to
    // orig. z_if_bit(dirty[i],bit) before+after nets Z^(bit . cy_{i+1}).
    for i in 0..l - 1 {
        circ.z_if_bit(dirty[i], bits[i]);
    }
    for q in y {
        circ.x(*q);
    }
    xor_carries_perpos(circ, e, d, base, y, dirty, carry_in);
    for q in y {
        circ.x(*q);
    }
    for i in 0..l - 1 {
        circ.z_if_bit(dirty[i], bits[i]);
    }
}

/// Clean carry-chain forward over a window, holding all b carries (carries[b-1] =
/// carry handed into the next dirty window).
fn clean_window_fwd(circ: &mut B, e: &QubitId, d: &QubitId, base: usize, y: &[QubitId], carries: &[QubitId]) {
    let b = y.len();
    assert_eq!(carries.len(), b);
    {
        let c0 = on_ctl(circ, e, d, base);
        if let Some(a0) = on_ctl_ref(&c0, e, d) {
            circ.ccx(a0, y[0], carries[0]);
        }
        on_ctl_free(circ, e, d, base, c0);
    }
    for i in 1..b {
        let ci = on_ctl(circ, e, d, base + i);
        let ai = on_ctl_ref(&ci, e, d);
        circ.cx(carries[i - 1], y[i]);
        if let Some(a) = &ai {
            circ.cx(*a, carries[i - 1]);
        }
        circ.ccx(y[i], carries[i - 1], carries[i]);
        if let Some(a) = &ai {
            circ.cx(*a, carries[i - 1]);
        }
        circ.cx(carries[i - 1], carries[i]);
        on_ctl_free(circ, e, d, base + i, ci);
    }
    for i in 0..b {
        let ci = on_ctl(circ, e, d, base + i);
        if let Some(a) = on_ctl_ref(&ci, e, d) {
            circ.cx(a, y[i]);
        }
        on_ctl_free(circ, e, d, base + i, ci);
    }
}

/// Reverse of [`clean_window_fwd`]: erase all b held carries (single-term erase =
/// hmr + cz_if_bit).
fn clean_window_rev(circ: &mut B, e: &QubitId, d: &QubitId, base: usize, y: &[QubitId], carries: Vec<QubitId>) {
    let b = y.len();
    let mut cy: Vec<Option<QubitId>> = carries.into_iter().map(Some).collect();
    for i in (1..b).rev() {
        let ci_ctl = on_ctl(circ, e, d, base + i);
        let actl = on_ctl_ref(&ci_ctl, e, d);
        if let Some(ai) = &actl {
            circ.cx(*ai, y[i]);
        }
        let next = cy[i].take().unwrap();
        let ci = cy[i - 1].take().unwrap();
        circ.cx(ci, next);
        if let Some(ai) = &actl {
            circ.cx(*ai, ci);
        }
        let bit = circ.alloc_bit();
        circ.hmr(next, bit);
        circ.zero_and_free(next);
        circ.cz_if_bit(y[i], ci, bit);
        if let Some(ai) = &actl {
            circ.cx(*ai, ci);
            circ.cx(*ai, y[i]);
        }
        on_ctl_free(circ, e, d, base + i, ci_ctl);
        cy[i - 1] = Some(ci);
    }
    let cy0 = cy[0].take().unwrap();
    let c0 = on_ctl(circ, e, d, base);
    if let Some(a0) = on_ctl_ref(&c0, e, d) {
        circ.cx(a0, y[0]);
        let bit = circ.alloc_bit();
        circ.hmr(cy0, bit);
        circ.zero_and_free(cy0);
        circ.cz_if_bit(y[0], a0, bit);
        circ.cx(a0, y[0]);
    } else {
        circ.zero_and_free(cy0);
    }
    on_ctl_free(circ, e, d, base, c0);
}

/// Build the fused fold at exactly `nv` clean vents:
/// nv==L-1 => full-clean; nv>=prop_from => clean prefix-tail; else dirty gradual.
fn build_fold_at(circ: &mut B, e: &QubitId, d: &QubitId, y: &[QubitId], dirty: &[QubitId], nv: usize) {
    let l = y.len();
    if nv >= l - 1 {
        // nv == L-1 is full-clean; nv > L-1 only happens for the unset-schedule
        // fallback (i32::MAX) -- treat as full-clean too.
        add_mf_fold_clean(circ, e, d, y);
        return;
    }
    // prop_from = top set bit of the +f fold addend (33) + 1 = 34.
    const PROP_FROM: usize = 34;
    if nv >= 1 && nv >= PROP_FROM {
        add_mf_fold_clean_tail(circ, e, d, y, Some(nv));
        return;
    }
    if nv == 0 {
        dirty_body(circ, e, d, 0, y, dirty, None);
    } else {
        let carries: Vec<QubitId> = (0..nv).map(|_| circ.alloc_qubit()).collect();
        clean_window_fwd(circ, e, d, 0, &y[..nv], &carries);
        let cin = carries[nv - 1];
        dirty_body(circ, e, d, nv, &y[nv..], &dirty[nv..], Some(&cin));
        clean_window_rev(circ, e, d, 0, &y[..nv], carries);
    }
}

/// Dispatch the fused fold on the schedule code: -s = chunked; else nv = clean vents
/// (full-clean / clean-tail / dirty gradual via [`build_fold_at`]).
fn fused_fold(circ: &mut B, e: &QubitId, d: &QubitId, ylow: &[QubitId], dirty: &[QubitId]) {
    let code = super::next_fold();
    if code < 0 {
        add_mf_fold_chunked(circ, e, d, ylow, (-code) as usize);
    } else {
        build_fold_at(circ, e, d, ylow, dirty, code as usize);
    }
}

/// Fused double-then-controlled-double: `y := y * 2 * (1 + s2) mod q`. `y.len() == n
/// == 256`; the doubling uses two transient overflow bits (a 258-bit working view),
/// not a persistent reg slot.
pub fn fused_double_cdouble(circ: &mut B, s2: &QubitId, y: &[QubitId]) {
    let n = 256usize;
    assert_eq!(y.len(), n, "fused double expects 256-bit y (transient overflow)");
    let _ = F_SECP256K1;
    let hi = circ.alloc_qubit();
    let hi2 = circ.alloc_qubit();
    // 258-bit view: y[0..n] ++ hi (256) ++ hi2 (257).
    let mut w: Vec<QubitId> = y.to_vec();
    w.push(hi);
    w.push(hi2);
    // shift 1 (unconditional).
    for i in (1..w.len()).rev() {
        circ.swap(w[i], w[i - 1]);
    }
    // shift 2 (s2-controlled).
    for i in (1..w.len()).rev() {
        circ.cswap(*s2, w[i], w[i - 1]);
    }
    // combined fold: add (e+2d)*f into the low LSBS window. e=w[256], d=w[257].
    // dirty borrow = the coordinate's own bits just above the fold window.
    let borrow: Vec<QubitId> = y[LSBS..2 * LSBS - 1].to_vec();
    fused_fold(circ, &w[n], &w[n + 1], &y[..LSBS], &borrow);
    // clear carry bits via bit identity (post-fold w[0]==e, w[1]==d).
    circ.cx(y[0], w[n]); // clear hi (e)
    clear_and(circ, &w[n + 1], s2, &y[1]); // clear hi2 (d = s2 & y[1])
    circ.zero_and_free(hi);
    circ.zero_and_free(hi2);
}

/// Exact gate-inverse: `y := y / (2*(1+s2)) mod q`. Reverse of [`fused_double_cdouble`]:
/// compute the overflow bits, subtract m*f, shift right.
pub fn fused_double_cdouble_reverse(circ: &mut B, s2: &QubitId, y: &[QubitId]) {
    let n = 256usize;
    assert_eq!(y.len(), n, "fused halve expects 256-bit y (transient overflow)");
    let hi = circ.alloc_qubit();
    let hi2 = circ.alloc_qubit();
    let mut w: Vec<QubitId> = y.to_vec();
    w.push(hi);
    w.push(hi2);
    // reversed carry-clear: compute e=y[0]->hi, d=(s2&y[1])->hi2.
    circ.ccx(*s2, y[1], w[n + 1]); // d
    circ.cx(y[0], w[n]); // e
    // subtract m*f from the low window: X-sandwich the forward fold.
    let borrow: Vec<QubitId> = y[LSBS..2 * LSBS - 1].to_vec();
    for q in &y[..LSBS] {
        circ.x(*q);
    }
    fused_fold(circ, &w[n], &w[n + 1], &y[..LSBS], &borrow);
    for q in &y[..LSBS] {
        circ.x(*q);
    }
    // shift right (inverse of the two left shifts): s2-controlled then unconditional.
    for i in 1..w.len() {
        circ.cswap(*s2, w[i], w[i - 1]);
    }
    for i in 1..w.len() {
        circ.swap(w[i], w[i - 1]);
    }
    circ.zero_and_free(hi);
    circ.zero_and_free(hi2);
}
