//! Fused double + controlled-double: `y := y * 2 * (1 + s2) mod q` with one
//! combined `(e+2d)*f` reduction instead of two separate `+f` const-adds.
//! AND-uncompute discharge is hmr+cz_if_bit; clear_and is inlined. This file
//! provides the full-clean fold (nv = L-1) plus the chunked and borrowed-dirty
//! variants (lower nv).

use super::arith::{F_SECP256K1, LSBS};
use super::{B, BExt};
use crate::circuit::{BitId, QubitId};
use std::cell::Cell;

thread_local! {
    static FOLD_CALL_INDEX: Cell<usize> = const { Cell::new(0) };
}

pub(super) fn reset_fold_call_index() {
    FOLD_CALL_INDEX.with(|index| index.set(0));
}

fn next_fold_call_index() -> usize {
    FOLD_CALL_INDEX.with(|index| {
        let current = index.get();
        index.set(current + 1);
        current
    })
}

fn env_index_value(name: &str, index: usize) -> Option<usize> {
    std::env::var(name)
        .ok()
        .and_then(|value| {
            value
                .split(',')
                .filter_map(|item| item.trim().split_once(':'))
                .find_map(|(call, value)| {
                    (call.parse::<usize>().ok()? == index)
                        .then(|| value.parse::<usize>().ok())
                        .flatten()
                })
        })
}

fn fold_call_reserve(index: usize, default: usize) -> usize {
    let base = env_index_value("TLM_TARGET_FOLD_CALL_RESERVES", index).unwrap_or(default);
    env_index_value("TLM_TARGET_FOLD_CALL_RESERVE_OVERRIDES", index).unwrap_or(base)
}

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

/// Toggle `d AND NOT e` into `dne`, given the live intersection `cc = e AND d`.
/// The Boolean identity `d & !e = d ^ (e & d)` replaces one CCX with two CX.
/// This is an involution, so the same sequence clears `dne` after use.
fn toggle_dnot_e_from_intersection(
    circ: &mut B,
    d: &QubitId,
    cc: &QubitId,
    dne: &QubitId,
) {
    circ.cx(*d, *dne);
    circ.cx(*cc, *dne);
}

/// Carry-propagate `c` into the pure-propagation tail `y[..]` via a cascade of
/// prefix-controlled increments (`mcx_clean_k`, log* ancillae): the clean-tail
/// fold's tail [nv, L).
fn add_carry_into_tail_prefix(circ: &mut B, y: &[QubitId], c: &QubitId) {
    if std::env::var("TLM_FOLD_TAIL_CINC").ok().as_deref() == Some("1") {
        // Exact replacement for the old quadratic prefix cascade.
        let yv: Vec<QubitId> = y.to_vec();
        super::mcx::cinc_khattar_gidney(circ, &yv, c);
        return;
    }

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
    toggle_dnot_e_from_intersection(
        circ,
        d,
        cc.as_ref().unwrap(),
        dne.as_ref().unwrap(),
    );
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
            toggle_dnot_e_from_intersection(circ, d, cc.as_ref().unwrap(), &dn);
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
            toggle_dnot_e_from_intersection(circ, d, cc.as_ref().unwrap(), &dn);
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
    toggle_dnot_e_from_intersection(circ, d, &cc, &dne);
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
    toggle_dnot_e_from_intersection(circ, d, &cc, &dne);
    (cc, sxor, sor, dne)
}

fn uncompute_fold_controls(circ: &mut B, e: &QubitId, d: &QubitId, cc: QubitId, sxor: QubitId, sor: QubitId, dne: QubitId) {
    toggle_dnot_e_from_intersection(circ, d, &cc, &dne);
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
fn fold_boundary_erase(circ: &mut B, ctl: &[Option<QubitId>], y: &[QubitId], cin: Option<&QubitId>, carry: QubitId) {
    if std::env::var("TLM_FOLD_BOUNDARY_ZERO_DIRECT")
        .ok()
        .as_deref()
        == Some("1")
        && cin.is_some()
        && ctl.iter().all(Option::is_none)
    {
        fold_boundary_erase_zero_direct(circ, y, cin.expect("cin checked"), carry);
        return;
    }
    let s = y.len();
    let temp: Vec<QubitId> = (0..s).map(|_| circ.alloc_qubit()).collect();
    for (i, c) in ctl.iter().enumerate() {
        if let Some(a) = c { circ.cx(*a, temp[i]); }
    }
    match cin {
        Some(cin) => super::arith::erase_carry_gated_opt(circ, None, y, &temp, cin, &carry, None),
        None => {
            super::arith::erase_carry_gated_zero_cin_opt(circ, None, y, &temp, &carry, None);
            circ.zero_and_free(carry);
        }
    }
    for (i, c) in ctl.iter().enumerate() {
        if let Some(a) = c { circ.cx(*a, temp[i]); }
    }
    for q in temp { circ.zero_and_free(q); }
}

fn fold_boundary_erase_zero_direct(circ: &mut B, y: &[QubitId], cin: &QubitId, carry: QubitId) {
    let n = y.len();
    assert!(n >= 1, "zero boundary erase needs >= 1 bit");
    let bit = circ.alloc_bit();
    circ.hmr(carry, bit);
    circ.zero_and_free(carry);
    circ.push_condition(bit);

    let mut cy: Vec<Option<QubitId>> = Vec::with_capacity(n);
    let c0 = circ.alloc_qubit();
    circ.x(c0);
    circ.cx(*cin, c0);
    cy.push(Some(c0));
    for i in 0..n - 1 {
        let next = circ.alloc_qubit();
        let ci = cy[i].as_ref().unwrap();
        circ.cx(*ci, y[i]);
        circ.x(*ci);
        circ.ccx(y[i], *ci, next);
        circ.x(*ci);
        circ.cx(*ci, next);
        cy.push(Some(next));
    }
    {
        let i = n - 1;
        let ci = cy[i].as_ref().unwrap();
        circ.cx(*ci, y[i]);
        circ.neg();
        circ.x(*ci);
        circ.cz(y[i], *ci);
        circ.x(*ci);
        circ.z(*ci);
        circ.cx(*ci, y[i]);
    }
    for i in (0..n - 1).rev() {
        let next = cy[i + 1].take().unwrap();
        let ci = cy[i].as_ref().unwrap();
        circ.cx(*ci, next);
        let mbit = circ.alloc_bit();
        circ.hmr(next, mbit);
        circ.zero_and_free(next);
        circ.x(*ci);
        circ.cz_if_bit(y[i], *ci, mbit);
        circ.x(*ci);
        circ.cx(*ci, y[i]);
    }
    let c0 = cy[0].take().unwrap();
    circ.cx(*cin, c0);
    circ.x(c0);
    circ.zero_and_free(c0);
    circ.pop_condition();
}

fn add_mf_fold_chunked(circ: &mut B, e: &QubitId, d: &QubitId, y: &[QubitId], s_chunk: usize) {
    let l = y.len();
    let release_controls = std::env::var("TLM_FOLD_RELEASE_CONTROLS")
        .ok()
        .as_deref()
        == Some("1");
    let zero_cin = std::env::var("TLM_FOLD_CHUNK_ZERO_CIN")
        .ok()
        .as_deref()
        == Some("1");
    let mut controls = Some(build_fold_controls(circ, e, d));
    let (cc, sxor, sor, dne) = controls.expect("fold controls present");
    let mut ctl = fold_ctl_map(*e, *d, cc, sxor, sor, dne, l);
    let cin0 = (!zero_cin).then(|| circ.alloc_qubit());
    let nch = l.div_ceil(s_chunk);
    let last_control_chunk = 11usize.min(l - 1) / s_chunk;
    let mut boundary: Vec<QubitId> = Vec::with_capacity(nch);
    for j in 0..nch {
        let lo = j * s_chunk;
        let hi = ((j + 1) * s_chunk).min(l);
        let cout = circ.alloc_qubit();
        let cin = if j == 0 {
            cin0.as_ref()
        } else {
            Some(&boundary[j - 1])
        };
        fold_chunk_clean(circ, &ctl[lo..hi], &y[lo..hi], cin, &cout);
        boundary.push(cout);
        if release_controls && j == last_control_chunk && j + 1 < nch {
            let (cc, sxor, sor, dne) = controls.take().expect("fold controls present");
            uncompute_fold_controls(circ, e, d, cc, sxor, sor, dne);
        }
    }
    for j in (0..nch).rev() {
        if release_controls && j == last_control_chunk && controls.is_none() {
            let rebuilt = build_fold_controls(circ, e, d);
            ctl = fold_ctl_map(*e, *d, rebuilt.0, rebuilt.1, rebuilt.2, rebuilt.3, l);
            controls = Some(rebuilt);
        }
        let lo = j * s_chunk;
        let hi = ((j + 1) * s_chunk).min(l);
        let bnd = boundary.pop().expect("boundary present");
        let cin = if j == 0 {
            cin0.as_ref()
        } else {
            Some(&boundary[j - 1])
        };
        fold_boundary_erase(circ, &ctl[lo..hi], &y[lo..hi], cin, bnd);
    }
    if let Some(cin0) = cin0 {
        circ.zero_and_free(cin0);
    }
    let (cc, sxor, sor, dne) = controls.take().expect("fold controls restored");
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

/// Measurement-based clear for the nonlinear on-demand controls.
///
/// HMR contributes the phase `m * q`, where `m` is the measurement bit.  The
/// feedback below uses the algebraic-normal forms
///
///   e | d      = e + d + e*d
///   d & !e     = d + e*d
///   e & d      = e*d
///
/// over GF(2), so it contributes the same phase and cancels HMR exactly.
fn on_ctl_clear_nonlinear_hmr(
    circ: &mut B,
    e: &QubitId,
    d: &QubitId,
    k: u8,
    q: &QubitId,
) {
    let bit = circ.alloc_bit();
    circ.hmr(*q, bit);
    match k {
        4 => {
            circ.z_if_bit(*e, bit);
            circ.z_if_bit(*d, bit);
            circ.cz_if_bit(*e, *d, bit);
        }
        5 => {
            circ.z_if_bit(*d, bit);
            circ.cz_if_bit(*e, *d, bit);
        }
        6 => circ.cz_if_bit(*e, *d, bit),
        _ => unreachable!("HMR clear requires a nonlinear fold control"),
    }
}

fn on_ctl_free(circ: &mut B, e: &QubitId, d: &QubitId, p: usize, c: OnCtl) {
    if let OnCtl::Owned(q) = c {
        let k = fold_ctl(p);
        let hmr_disabled = std::env::var("TLM_FOLD_HMR_CONTROL_CLEANUP_DISABLE")
            .ok()
            .as_deref()
            == Some("1");
        if k == 3 || hmr_disabled {
            // XOR is cheaper reversibly; the flag preserves the exact old baseline.
            on_ctl_apply(circ, e, d, k, &q);
        } else {
            on_ctl_clear_nonlinear_hmr(circ, e, d, k, &q);
        }
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
    let call_index = next_fold_call_index();
    let timeline_start = circ.active_timeline.len();
    let entry_active = circ.active_qubits;
    let code = super::next_fold();
    let mut selected_nv = None;
    if code < 0 {
        let chunk = std::env::var("TLM_FOLD_CHUNK_FORCE")
            .ok()
            .and_then(|value| value.parse::<usize>().ok())
            .filter(|&value| value > 0)
            .unwrap_or((-code) as usize);
        add_mf_fold_chunked(circ, e, d, ylow, chunk);
    } else {
        let default_reserve = std::env::var("TLM_TARGET_FOLD_RESERVE")
            .ok()
            .and_then(|value| value.parse::<usize>().ok())
            .unwrap_or(4);
        let reserve = fold_call_reserve(call_index, default_reserve);
        let nv = super::target_qubit_headroom(circ)
            .map_or(code as usize, |headroom| {
                (code as usize).min(headroom.saturating_sub(reserve))
            });
        selected_nv = Some(nv);
        build_fold_at(circ, e, d, ylow, dirty, nv);
    }
    if std::env::var_os("TRACE_TLM_FOLD").is_some() {
        let local_peak = circ.active_timeline[timeline_start..]
            .iter()
            .map(|(_, active)| *active)
            .max()
            .unwrap_or(circ.active_qubits);
        eprintln!(
            "TLM_FOLD call={} phase={} code={} nv={} entry_active={} local_peak={} ops={}",
            call_index,
            circ.phase,
            code,
            selected_nv.map_or(-1, |value| value as i32),
            entry_active,
            local_peak,
            circ.current_ops_len(),
        );
    }
}

fn fused_fold_e_only(circ: &mut B, e: &QubitId, y: &[QubitId]) {
    let call_index = next_fold_call_index();
    let timeline_start = circ.active_timeline.len();
    let entry_active = circ.active_qubits;
    let code = super::next_fold();
    let default_reserve = std::env::var("TLM_TARGET_FOLD_RESERVE")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(4);
    let reserve = fold_call_reserve(call_index, default_reserve);
    let g = if code < 0 {
        0
    } else {
        super::target_qubit_headroom(circ)
            .map_or(code as usize, |headroom| {
                (code as usize).min(headroom.saturating_sub(reserve))
            })
            .min(LSBS - 1)
    };
    let f_bytes = F_SECP256K1.to_le_bytes();
    super::arith::add_f_window_pub(circ, e, y, LSBS, &f_bytes, Some(g));
    if std::env::var_os("TRACE_TLM_FOLD").is_some() {
        let local_peak = circ.active_timeline[timeline_start..]
            .iter()
            .map(|(_, active)| *active)
            .max()
            .unwrap_or(circ.active_qubits);
        eprintln!(
            "TLM_FOLD call={} phase={} code={} nv={} entry_active={} local_peak={} ops={}",
            call_index,
            circ.phase,
            code,
            g as i32,
            entry_active,
            local_peak,
            circ.current_ops_len(),
        );
    }
}

fn trace_fold_alloc(circ: &B, name: &str, stage: &str, i: usize) {
    if std::env::var_os("TRACE_TLM_FOLD_ALLOC").is_some() {
        eprintln!(
            "TLM_FOLD_ALLOC name={name} stage={stage} i={i} active={} ops={}",
            circ.active_qubits,
            circ.current_ops_len(),
        );
    }
}

/// Fused double-then-controlled-double: `y := y * 2 * (1 + s2) mod q`. `y.len() == n
/// == 256`; the doubling uses two transient overflow bits (a 258-bit working view),
/// not a persistent reg slot.
pub fn fused_double_cdouble(circ: &mut B, s2: &QubitId, y: &[QubitId]) {
    maybe_run_gradual_fold_nonlinear_control_hmr_selftest();
    let n = 256usize;
    assert_eq!(y.len(), n, "fused double expects 256-bit y (transient overflow)");
    let _ = F_SECP256K1;
    trace_fold_alloc(circ, "fwd_cdouble", "entry", usize::MAX);
    let hi = circ.alloc_qubit();
    trace_fold_alloc(circ, "fwd_cdouble", "after_hi", usize::MAX);
    let hi2 = circ.alloc_qubit();
    trace_fold_alloc(circ, "fwd_cdouble", "after_hi2", usize::MAX);
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

pub fn fused_double_only(circ: &mut B, y: &[QubitId]) {
    let n = 256usize;
    assert_eq!(y.len(), n, "fused double expects 256-bit y");
    trace_fold_alloc(circ, "fwd_only", "entry", usize::MAX);
    let hi = circ.alloc_qubit();
    trace_fold_alloc(circ, "fwd_only", "after_hi", usize::MAX);
    let mut w: Vec<QubitId> = y.to_vec();
    w.push(hi);
    for i in (1..w.len()).rev() {
        circ.swap(w[i], w[i - 1]);
    }
    fused_fold_e_only(circ, &w[n], y);
    circ.cx(y[0], w[n]);
    circ.zero_and_free(hi);
}

/// Exact gate-inverse: `y := y / (2*(1+s2)) mod q`. Reverse of [`fused_double_cdouble`]:
/// compute the overflow bits, subtract m*f, shift right.
pub fn fused_double_cdouble_reverse(circ: &mut B, s2: &QubitId, y: &[QubitId]) {
    maybe_run_gradual_fold_nonlinear_control_hmr_selftest();
    let n = 256usize;
    assert_eq!(y.len(), n, "fused halve expects 256-bit y (transient overflow)");
    trace_fold_alloc(circ, "rev_cdouble", "entry", usize::MAX);
    let hi = circ.alloc_qubit();
    trace_fold_alloc(circ, "rev_cdouble", "after_hi", usize::MAX);
    let hi2 = circ.alloc_qubit();
    trace_fold_alloc(circ, "rev_cdouble", "after_hi2", usize::MAX);
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

pub fn fused_double_only_reverse(circ: &mut B, y: &[QubitId]) {
    let n = 256usize;
    assert_eq!(y.len(), n, "fused halve expects 256-bit y");
    trace_fold_alloc(circ, "rev_only", "entry", usize::MAX);
    let hi = circ.alloc_qubit();
    trace_fold_alloc(circ, "rev_only", "after_hi", usize::MAX);
    let mut w: Vec<QubitId> = y.to_vec();
    w.push(hi);
    circ.cx(y[0], w[n]);
    for q in &y[..LSBS] {
        circ.x(*q);
    }
    fused_fold_e_only(circ, &w[n], y);
    for q in &y[..LSBS] {
        circ.x(*q);
    }
    for i in 1..w.len() {
        circ.swap(w[i], w[i - 1]);
    }
    circ.zero_and_free(hi);
}

/// Exhaust every `(e,d,HMR outcome)` basis case for each nonlinear derived
/// control.  This is available in normal builds because this worktree's wider
/// `cfg(test)` suite currently has unrelated compile failures.
fn gradual_fold_nonlinear_control_hmr_selftest() {
    use crate::circuit::OperationType;
    use crate::sim::Simulator;
    use sha3::{
        digest::{ExtendableOutput, Update},
        Shake128,
    };

    // Positions 8, 10, and 11 select OR, d&!e, and AND respectively.
    for &(position, kind) in &[(8usize, 4u8), (10, 5), (11, 6)] {
        assert_eq!(fold_ctl(position), kind);

        let mut circ = B::new();
        let e = circ.alloc_qubit();
        let d = circ.alloc_qubit();
        let q = circ.alloc_qubit();
        on_ctl_apply(&mut circ, &e, &d, kind, &q);
        on_ctl_free(&mut circ, &e, &d, position, OnCtl::Owned(q));

        assert_eq!(circ.active_qubits, 2, "owned control was not released");
        assert_eq!(circ.peak_qubits, 3, "cleanup increased peak width");
        assert_eq!(circ.next_bit, 1, "expected one HMR result bit");
        assert_eq!(
            circ.ops
                .iter()
                .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
                .count(),
            1,
            "cleanup must add no Toffoli-class gate",
        );
        assert_eq!(
            circ.ops
                .iter()
                .filter(|op| op.kind == OperationType::Hmr)
                .count(),
            1,
        );

        // Pack every (e,d) basis state sixteen times.  The fixed SHAKE seed
        // deterministically supplies both HMR outcomes for every state; the
        // coverage assertion below keeps that part of the proof explicit.
        let mut e_mask = 0u64;
        let mut d_mask = 0u64;
        for shot in 0..64usize {
            let state = shot & 3;
            e_mask |= ((state & 1) as u64) << shot;
            d_mask |= (((state >> 1) & 1) as u64) << shot;
        }

        let mut seed = Shake128::default();
        seed.update(b"gradual-fold-derived-control-hmr");
        seed.update(&[kind]);
        let mut xof = seed.finalize_xof();
        let mut sim = Simulator::new(
            circ.next_qubit as usize,
            circ.next_bit as usize,
            &mut xof,
        );
        *sim.qubit_mut(e) = e_mask;
        *sim.qubit_mut(d) = d_mask;
        sim.apply_iter(circ.ops.iter());

        assert_eq!(sim.qubit(e), e_mask, "e changed for control kind {kind}");
        assert_eq!(sim.qubit(d), d_mask, "d changed for control kind {kind}");
        assert_eq!(sim.qubit(q), 0, "owned control remained dirty for kind {kind}");
        assert_eq!(sim.phase, 0, "phase feedback failed for control kind {kind}");

        let measured = sim.bits[0];
        for state in 0..4usize {
            let mut outcomes = 0u8;
            for shot in (state..64usize).step_by(4) {
                outcomes |= 1 << ((measured >> shot) & 1);
            }
            assert_eq!(
                outcomes, 0b11,
                "HMR outcomes not exhaustive for kind {kind}, state {state}",
            );
        }
    }
}

fn maybe_run_gradual_fold_nonlinear_control_hmr_selftest() {
    if std::env::var_os("TLM_FOLD_HMR_CONTROL_SELFTEST").is_none() {
        return;
    }
    static SELFTEST: std::sync::Once = std::sync::Once::new();
    SELFTEST.call_once(|| {
        gradual_fold_nonlinear_control_hmr_selftest();
        eprintln!("TLM_FOLD_HMR_CONTROL_SELFTEST_OK");
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gradual_fold_nonlinear_control_hmr_cleanup_is_exact() {
        gradual_fold_nonlinear_control_hmr_selftest();
    }
}
