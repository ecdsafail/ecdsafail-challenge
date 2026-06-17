//! Reversible UNPACKED PZ modular inversion. Separate registers A,B,|a|,|b| (no
//! cursor), reversible via the PZ cofactor ratio (no spooky pebbling). The whole
//! thing is built from ONE primitive -- restoring long division -- used forward
//! on the gcd pair and (its reverse) as the cofactor multiply.
//!
//! `long_division(A,B,q)`: A := A mod B, q := A // B (q starts |0>). Reversible.
//! `long_division_reverse(A,B,q)`: the inverse -- A := A + q*B, q := 0 (consumes
//! q). This is exactly the consuming multiply `|a| += q|b|` applied to (|a|,|b|).

use crate::point_add::trailmix_port::arith::cuccaro::controlled_add_cuccaro_3n_refs;
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

/// `out ^= (v < u)` (MAJ cascade + carry capture + un-MAJ), v,u restored. Refs
/// variant of `two_cursor::borrow_compare` (windows here are non-contiguous).
pub(crate) fn borrow_compare_refs(circ: &mut Circuit, v: &[&QReg], u: &[&QReg], out: &QReg) {
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
    let cc = circ.alloc_qreg("bc.c");
    circ.cx(b[0], a[0]);
    circ.cx(b[0], &cc);
    circ.ccx(&cc, a[0], b[0]);
    for i in 1..n {
        circ.cx(b[i], a[i]);
        circ.cx(b[i], b[i - 1]);
        circ.ccx(b[i - 1], a[i], b[i]);
    }
    circ.cx(b[n - 1], out);
    for i in (1..n).rev() {
        circ.ccx(b[i - 1], a[i], b[i]);
        circ.cx(b[i], b[i - 1]);
        circ.cx(b[i], a[i]);
    }
    circ.ccx(&cc, a[0], b[0]);
    circ.cx(b[0], &cc);
    circ.cx(b[0], a[0]);
    circ.zero_and_free(cc);
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
