//! Streaming carry-save (partial Wallace) integer squarer: `r = y*y` where `r`
//! is 2n bits (|0> on entry). Exploits the squaring symmetry (only n(n-1)/2
//! cross-products + n free diagonals) and carry-save column compression (no
//! full-width carry ripple), so it costs ~n^2 Toffoli vs the Horner
//! full-multiply's ~2.7n^2 -- see `notes/wallace_square_finding.md`.
//!
//! This file: the value-correct carry-save core (garbage kept for a gate-by-gate
//! Bennett reverse). The streaming measurement-uncompute that bounds the live
//! scratch to ~n is layered on once the core is validated.

use crate::point_add::trailmix_port::circuit::{Circuit, QReg};

/// Reversible full adder (3:2 compressor). `b` becomes the sum `a^b^cin`;
/// `cout` (fresh |0>) becomes the carry `MAJ(a,b,cin)`; `a`,`cin` are unchanged
/// (carry-save garbage). 2 Toffoli.
pub(crate) fn fa(circ: &mut Circuit, a: &QReg, b: &QReg, cin: &QReg, cout: &QReg) {
    circ.ccx(a, b, cout); // cout = a & b_orig
    circ.cx(a, b); // b = a ^ b_orig
    circ.ccx(b, cin, cout); // cout = a&b_orig ^ (a^b_orig)&cin = MAJ
    circ.cx(cin, b); // b = a ^ b_orig ^ cin = sum
}

/// Reversible half adder. `b` becomes `a^b`; `cout` (fresh |0>) becomes `a&b`.
/// 1 Toffoli.
pub(crate) fn ha(circ: &mut Circuit, a: &QReg, b: &QReg, cout: &QReg) {
    circ.ccx(a, b, cout);
    circ.cx(a, b);
}

/// The (i,k) cross-product pairs (i<k) landing in result column `c=i+k+1`.
fn cross_pairs(n: usize) -> Vec<Vec<(usize, usize)>> {
    let mut cols = vec![Vec::new(); 2 * n + 2];
    for i in 0..n {
        for k in (i + 1)..n {
            cols[i + k + 1].push((i, k));
        }
    }
    cols
}

/// `r[0..2n] = y^2` (integer), `r` enters |0>. `garbage` collects every scratch
/// qubit (products + carry-save intermediates); the caller gate-reverses this
/// call to clean them. Value-correct reference; not yet scratch-bounded.
pub fn wallace_square_bennett(circ: &mut Circuit, y: &[QReg], r: &[QReg], garbage: &mut Vec<QReg>) {
    let n = y.len();
    assert_eq!(r.len(), 2 * n, "result must be 2n bits");
    let crosses = cross_pairs(n);
    let mut active: Vec<Vec<QReg>> = (0..2 * n + 2).map(|_| Vec::new()).collect();

    for c in 0..2 * n {
        // 1) materialize column c's fresh inputs.
        if c % 2 == 0 && c / 2 < n {
            let d = circ.alloc_qreg("wsq_diag");
            circ.cx(&y[c / 2], &d); // diagonal y_{c/2}
            active[c].push(d);
        }
        for &(i, k) in &crosses[c] {
            let p = circ.alloc_qreg("wsq_x");
            circ.ccx(&y[i], &y[k], &p); // cross y_i & y_k
            active[c].push(p);
        }
        // 2) carry-save compress -> one survivor bit + carries to c+1.
        let mut carries_up: Vec<QReg> = Vec::new();
        while active[c].len() > 1 {
            if active[c].len() >= 3 {
                let cin = active[c].pop().unwrap();
                let b = active[c].pop().unwrap();
                let a = active[c].pop().unwrap();
                let cout = circ.alloc_qreg("wsq_cy");
                fa(circ, &a, &b, &cin, &cout);
                carries_up.push(cout);
                active[c].push(b); // sum stays
                garbage.push(a);
                garbage.push(cin);
            } else {
                let b = active[c].pop().unwrap();
                let a = active[c].pop().unwrap();
                let cout = circ.alloc_qreg("wsq_cy");
                ha(circ, &a, &b, &cout);
                carries_up.push(cout);
                active[c].push(b);
                garbage.push(a);
            }
        }
        // 3) the survivor is r[c].
        if let Some(bit) = active[c].pop() {
            circ.cx(&bit, &r[c]);
            garbage.push(bit);
        }
        active[c + 1].extend(carries_up);
    }
    for col in &mut active {
        while let Some(q) = col.pop() {
            garbage.push(q);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{thread_rng, Rng};

    #[test]
    fn wallace_square_value_small_n() {
        // value-correctness of the carry-save squarer in the real sim, n=4..10.
        for n in [4usize, 5, 6, 8, 10] {
            let mut circ = Circuit::new();
            circ.set_max_qubit_peak(4096);
            let y: Vec<QReg> = (0..n).map(|i| circ.alloc_qreg(&format!("y{i}"))).collect();
            let r: Vec<QReg> = (0..2 * n)
                .map(|i| circ.alloc_qreg(&format!("r{i}")))
                .collect();
            // load 64 random y values, one per shot.
            let mut yvals = vec![0u64; 64];
            let mut rng = thread_rng();
            for (shot, yv) in yvals.iter_mut().enumerate() {
                *yv = rng.gen::<u64>() & ((1u64 << n) - 1);
                for (i, yq) in y.iter().enumerate() {
                    if ((*yv >> i) & 1) == 1 {
                        circ.sim_load_reg_bytes_shot(std::slice::from_ref(yq), &[1], shot);
                    }
                }
            }
            let mut garbage: Vec<QReg> = Vec::new();
            wallace_square_bennett(&mut circ, &y, &r, &mut garbage);

            // r first (read below); then y + garbage so destroy_sim sees every
            // live QReg (garbage is nonzero scratch -- not yet Bennett-reversed).
            let mut outs: Vec<QReg> = Vec::new();
            outs.extend(r);
            outs.extend(y);
            outs.extend(garbage);
            let (sim, det) = circ.destroy_sim(outs);
            for (shot, &yv) in yvals.iter().enumerate() {
                let mut got: u128 = 0;
                for (i, rq) in det.iter().take(2 * n).enumerate() {
                    if sim.read_bit_shot(rq, shot) == 1 {
                        got |= 1u128 << i;
                    }
                }
                let want = (yv as u128) * (yv as u128);
                assert_eq!(got, want, "n={n} shot={shot} y={yv}: got {got} want {want}");
            }
        }
    }
}
