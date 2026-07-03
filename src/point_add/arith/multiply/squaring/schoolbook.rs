use super::*;

// ─── 2-level Karatsuba variants (recursive on inner half-mults) ───
// Costs 2 extra z1_inner registers of ~2*(n/4+1) qubits each (~260 total for n=256).
// Higher peak qubits; use only at low-peak mul sites.

/// Symmetric schoolbook for squaring: x² = sum_i x[i]·2^(2i) + sum_{i<j} 2·x[i]·x[j]·2^(i+j).
/// Each cross-product is computed ONCE (instead of twice in full schoolbook),
/// halving the AND count + Cuccaro_add length. Saves ~130k CCX per squaring.
///
/// Row i layout (width n-i): bit 0 = diagonal x[i] at position 2i, bit 1 = 0
/// (gap), bit k+2 = cross-product (x[i] AND x[i+1+k]) at position i+(i+1+k)+1.
#[allow(dead_code)] // retained reference/alternative impl; not on active build path
pub(crate) fn schoolbook_square_symmetric(b: &mut B, x: &[QubitId], tmp_ext: &[QubitId]) {
    let n = x.len();
    debug_assert_eq!(tmp_ext.len(), 2 * n);
    for i in 0..n {
        // Width: bit 0 = diag at pos 2i, bit 1 = gap, bits 2..(n-i) = cross-
        // products at positions 2i+2..i+n. Last bit index = n-i, so width = n-i+1.
        // Edge case: i = n-1 has only the diagonal, width = 1.
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let num_cross = if i + 1 < n { n - i - 1 } else { 0 };
        // num_cross = number of cross-products in this row = width - 2 when width >= 2.
        let row = b.alloc_qubits(width);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            b.ccx(x[i], x[i + 1 + k], row[k + 2]);
        }
        let pad = b.alloc_qubit();
        let mut row_padded = row.clone();
        row_padded.push(pad);
        let slice: Vec<QubitId> = tmp_ext[2 * i..2 * i + width + 1].to_vec();
        let c_in = b.alloc_qubit();
        cuccaro_add_fast(b, &row_padded, &slice, c_in);
        b.free(c_in);
        b.free(pad);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            let m = b.alloc_bit();
            b.hmr(row[k + 2], m);
            b.cz_if(x[i], x[i + 1 + k], m);
        }
        b.free_vec(&row);
    }
}

#[allow(dead_code)] // retained reference/alternative impl; not on active build path
pub(crate) fn schoolbook_square_symmetric_lowq(b: &mut B, x: &[QubitId], tmp_ext: &[QubitId]) {
    let n = x.len();
    debug_assert_eq!(tmp_ext.len(), 2 * n);
    for i in 0..n {
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let num_cross = if i + 1 < n { n - i - 1 } else { 0 };
        let row = b.alloc_qubits(width);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            b.ccx(x[i], x[i + 1 + k], row[k + 2]);
        }
        let pad = b.alloc_qubit();
        let mut row_padded = row.clone();
        row_padded.push(pad);
        let slice: Vec<QubitId> = tmp_ext[2 * i..2 * i + width + 1].to_vec();
        let c_in = b.alloc_qubit();
        cuccaro_add(b, &row_padded, &slice, c_in);
        b.free(c_in);
        b.free(pad);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            let m = b.alloc_bit();
            b.hmr(row[k + 2], m);
            b.cz_if(x[i], x[i + 1 + k], m);
        }
        b.free_vec(&row);
    }
}

/// Like `schoolbook_square_symmetric` (fast, measurement UMA) but the per-row
/// Cuccaro carry lane is hosted on a caller-supplied clean register `host`
/// (returned clean) instead of a fresh allocation. Toffoli-identical to the
/// fast square, peak-identical to the lowq square — used for the z0 lobe of the
/// round84 Karatsuba square, where the not-yet-written z2 slice is clean scratch.
#[allow(dead_code)] // retained reference/alternative impl; not on active build path
pub(crate) fn schoolbook_square_symmetric_hosted(
    b: &mut B,
    x: &[QubitId],
    tmp_ext: &[QubitId],
    host: &[QubitId],
) {
    let n = x.len();
    debug_assert_eq!(tmp_ext.len(), 2 * n);
    if square_selfhost_safe_lane_reuse_enabled() {
        assert_qubit_slices_disjoint(&[x, tmp_ext, host]);
    }
    for i in 0..n {
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let num_cross = if i + 1 < n { n - i - 1 } else { 0 };
        let row = b.alloc_qubits(width);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            b.ccx(x[i], x[i + 1 + k], row[k + 2]);
        }
        let slice: Vec<QubitId> = tmp_ext[2 * i..2 * i + width + 1].to_vec();
        if square_selfhost_safe_lane_reuse_enabled() {
            // The z2 sibling host is clean and disjoint from x and z0.  It has
            // ample room for both the width carry lanes and one clean c_in.
            assert!(host.len() > width);
            cuccaro_add_fast_low_to_ext_borrowed_carries(
                b,
                &row,
                &slice,
                host[width],
                &host[..width],
            );
        } else {
            let pad = b.alloc_qubit();
            let mut row_padded = row.clone();
            row_padded.push(pad);
            let c_in = b.alloc_qubit();
            cuccaro_add_fast_borrowed_carries(
                b,
                &row_padded,
                &slice,
                c_in,
                &host[..row_padded.len() - 1],
            );
            b.free(c_in);
            b.free(pad);
        }
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            let m = b.alloc_bit();
            b.hmr(row[k + 2], m);
            b.cz_if(x[i], x[i + 1 + k], m);
        }
        b.free_vec(&row);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::Simulator;
    use sha3::{
        digest::{ExtendableOutput, Update, XofReader},
        Shake128,
    };

    fn get<R: XofReader>(sim: &Simulator<'_, R>, qs: &[QubitId], shot: usize) -> u64 {
        let mut v = 0u64;
        for (i, &q) in qs.iter().enumerate() {
            v |= ((sim.qubit(q) >> shot) & 1) << i;
        }
        v
    }

    /// `schoolbook_square_symmetric` writes `x²` (2n-bit) into a zero `tmp_ext`,
    /// leaving `x` unchanged and every internal ancilla back at |0>. Exhaustive
    /// over all n=4 inputs (x ∈ [0,16)).
    #[test]
    fn schoolbook_square_symmetric_computes_x_squared() {
        const N: usize = 4;
        let mut b = B::new();
        let x = b.alloc_qubits(N);
        let tmp = b.alloc_qubits(2 * N);
        schoolbook_square_symmetric(&mut b, &x, &tmp);
        let nq = b.next_qubit as usize;
        let nb = b.next_bit as usize;
        let inputs: std::collections::HashSet<u64> =
            x.iter().chain(tmp.iter()).map(|q| q.0).collect();

        let mut seed = Shake128::default();
        seed.update(b"sqsym-n4");
        let mut xof = seed.finalize_xof();
        let mut sim = Simulator::new(nq, nb, &mut xof);
        for shot in 0..(1 << N) {
            for (i, &q) in x.iter().enumerate() {
                if (shot >> i) & 1 != 0 {
                    *sim.qubit_mut(q) |= 1u64 << shot;
                }
            }
        }
        sim.apply_iter(b.ops.iter());
        assert_eq!(sim.phase, 0, "phase garbage");
        for shot in 0..(1 << N) {
            let xv = shot as u64;
            assert_eq!(get(&sim, &tmp, shot), xv * xv, "x={xv}: tmp != x^2");
            assert_eq!(get(&sim, &x, shot), xv, "x={xv} changed");
        }
        for q in 0..nq as u64 {
            if !inputs.contains(&q) {
                assert_eq!(sim.qubit(QubitId(q)), 0, "ancilla q{q} not clean");
            }
        }
    }
}
