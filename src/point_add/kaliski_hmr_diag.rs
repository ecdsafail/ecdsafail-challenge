//! Diagnostics for HMR sequences in generic vs specialized Kaliski prefix steps.
//!
//! This checks a concrete hypothesis: if the specialized prefix emits a
//! different HMR sequence from the generic prefix, then even classically
//! equivalent state updates could accumulate a phase mismatch under the
//! measurement-based uncompute scheme.

use crate::circuit::{OperationType, QubitId};

use super::{
    B, N, kaliski_iteration, kaliski_iteration_bulk_prefix3, SECP256K1_P,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HmrOp {
    pub q: QubitId,
}

fn extract_hmrs(ops: &[crate::circuit::Op]) -> Vec<HmrOp> {
    ops.iter()
        .filter(|op| matches!(op.kind, OperationType::Hmr))
        .map(|op| HmrOp { q: op.q_target })
        .collect()
}

fn build_generic_step0() -> Vec<HmrOp> {
    let mut b = B::new();
    let u = b.alloc_qubits(N);
    let v = b.alloc_qubits(N);
    let r = b.alloc_qubits(N);
    let s = b.alloc_qubits(N);
    let m = b.alloc_qubit();
    let f = b.alloc_qubit();
    // This is the actual generic iter-0 body; the f qubit starts at |0> in the
    // raw builder, but the emitted HMR pattern is what we want to compare.
    kaliski_iteration(&mut b, SECP256K1_P, &u, &v, &r, &s, m, f, 0);
    extract_hmrs(&b.ops)
}

fn build_special_step0() -> Vec<HmrOp> {
    let mut b = B::new();
    let u = b.alloc_qubits(N);
    let v = b.alloc_qubits(N);
    let r = b.alloc_qubits(N);
    let s = b.alloc_qubits(N);
    let m = b.alloc_qubit();
    kaliski_iteration_bulk_prefix3(&mut b, &u, &v, &r, &s, m, 0);
    extract_hmrs(&b.ops)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compare_hmr_sequences_step0() {
        let g = build_generic_step0();
        let s = build_special_step0();
        eprintln!("=== HMR diag: generic vs specialized step0 ===");
        eprintln!("generic HMR count     : {}", g.len());
        eprintln!("specialized HMR count : {}", s.len());
        let common = g.len().min(s.len());
        let mut first_diff = None;
        for i in 0..common {
            if g[i] != s[i] {
                first_diff = Some(i);
                break;
            }
        }
        eprintln!("first differing index : {:?}", first_diff);
        if let Some(i) = first_diff {
            eprintln!("generic[{}] = {:?}", i, g[i]);
            eprintln!("special[{}] = {:?}", i, s[i]);
        }
        eprintln!("============================================");
        assert!(g.len() >= s.len());
    }
}
