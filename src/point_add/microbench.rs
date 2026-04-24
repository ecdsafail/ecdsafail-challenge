//! Micro-benchmark module for exploring structural changes WITHOUT
//! touching the top-level point-add benchmark.
//!
//! Each function builds a sub-circuit in isolation, reports:
//!   - local peak qubits (upper bound on what the primitive adds on top of
//!     its input register footprint)
//!   - emitted CCX+CCZ ("Toffoli") count
//!
//! Run with:
//!   cargo test -p quantum_ecc microbench -- --nocapture
//!
//! All tests are opt-in and slow, so they are gated behind the
//! `MICROBENCH=1` environment variable.

use alloy_primitives::U256;

use super::{
    add_nbit_qq, cuccaro_add_fast, cuccaro_sub_fast, mod_mul_add_into_acc_karatsuba,
    mod_mul_add_into_acc_karatsuba2, mod_mul_add_into_acc_karatsuba_lowq,
    mod_mul_add_into_acc_schoolbook, mod_mul_write_into_zero_acc_karatsuba,
    mod_mul_write_into_zero_acc_karatsuba2, mod_mul_write_into_zero_acc_schoolbook,
    schoolbook_mul_into_addsub, schoolbook_mul_into_addsub_inverse,
    schoolbook_mul_into_addsub_lowq, schoolbook_mul_into_addsub_lowq_inverse, B, N, SECP256K1_P,
};
use crate::circuit::{Op, OperationType};

fn count_toffoli(ops: &[Op]) -> usize {
    ops.iter()
        .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
        .count()
}

struct Measured {
    total_qubits: u32,
    peak_qubits: u32,
    toffoli: usize,
    ops: usize,
}

fn measure<F: FnOnce(&mut B)>(f: F) -> Measured {
    let mut b = B::new();
    f(&mut b);
    Measured {
        total_qubits: b.next_qubit,
        peak_qubits: b.peak_qubits,
        toffoli: count_toffoli(&b.ops),
        ops: b.ops.len(),
    }
}

fn print_row(label: &str, m: &Measured) {
    println!(
        "microbench | {:<42} | toffoli={:>8} ops={:>9} peak_q={:>5} total_q={:>5}",
        label, m.toffoli, m.ops, m.peak_qubits, m.total_qubits
    );
}

fn enabled() -> bool {
    std::env::var("MICROBENCH")
        .ok()
        .map(|v| v != "0")
        .unwrap_or(false)
}

fn fill_x(b: &mut B, n: usize) -> Vec<super::QubitId> {
    // Allocate n qubits, pretend they are inputs (no initial X to keep it
    // simple; measurements here don't care about semantic value, only cost).
    b.alloc_qubits(n)
}

/// Build just `x + a` (non-modular, n-bit) via `cuccaro_add_fast` to
/// measure the Cuccaro-fast peak/Toffoli envelope.
fn bench_cuccaro_add_fast(n: usize) -> Measured {
    measure(|b| {
        let x = fill_x(b, n);
        let acc = fill_x(b, n);
        let c_in = b.alloc_qubit();
        cuccaro_add_fast(b, &x, &acc, c_in);
        b.free(c_in);
        b.free_vec(&acc);
        b.free_vec(&x);
    })
}

/// Same as above but using the non-fast `add_nbit_qq` (uses a single c_in
/// ancilla and no carry register). Trades Toffoli ↑ for peak ↓.
fn bench_cuccaro_add_slow(n: usize) -> Measured {
    measure(|b| {
        let x = fill_x(b, n);
        let acc = fill_x(b, n);
        add_nbit_qq(b, &x, &acc);
        b.free_vec(&acc);
        b.free_vec(&x);
    })
}

/// Reverse direction of cuccaro_add_fast (for symmetry sanity only).
fn bench_cuccaro_sub_fast(n: usize) -> Measured {
    measure(|b| {
        let x = fill_x(b, n);
        let acc = fill_x(b, n);
        let c_in = b.alloc_qubit();
        cuccaro_sub_fast(b, &x, &acc, c_in);
        b.free(c_in);
        b.free_vec(&acc);
        b.free_vec(&x);
    })
}

/// write_into_zero schoolbook: tmp_ext and inv done inside the primitive.
fn bench_mul_schoolbook_write() -> Measured {
    let p = SECP256K1_P;
    measure(|b| {
        let x = fill_x(b, N);
        let y = fill_x(b, N);
        let acc = fill_x(b, N); // assumed zero
        mod_mul_write_into_zero_acc_schoolbook(b, &acc, &x, &y, p);
        b.free_vec(&acc);
        b.free_vec(&y);
        b.free_vec(&x);
    })
}

fn bench_mul_add_schoolbook() -> Measured {
    let p = SECP256K1_P;
    measure(|b| {
        let x = fill_x(b, N);
        let y = fill_x(b, N);
        let acc = fill_x(b, N);
        mod_mul_add_into_acc_schoolbook(b, &acc, &x, &y, p);
        b.free_vec(&acc);
        b.free_vec(&y);
        b.free_vec(&x);
    })
}

fn bench_mul_karatsuba_write() -> Measured {
    let p = SECP256K1_P;
    measure(|b| {
        let x = fill_x(b, N);
        let y = fill_x(b, N);
        let acc = fill_x(b, N);
        mod_mul_write_into_zero_acc_karatsuba(b, &acc, &x, &y, p);
        b.free_vec(&acc);
        b.free_vec(&y);
        b.free_vec(&x);
    })
}

fn bench_mul_add_karatsuba() -> Measured {
    let p = SECP256K1_P;
    measure(|b| {
        let x = fill_x(b, N);
        let y = fill_x(b, N);
        let acc = fill_x(b, N);
        mod_mul_add_into_acc_karatsuba(b, &acc, &x, &y, p);
        b.free_vec(&acc);
        b.free_vec(&y);
        b.free_vec(&x);
    })
}

fn bench_mul_karatsuba2_write() -> Measured {
    let p = SECP256K1_P;
    measure(|b| {
        let x = fill_x(b, N);
        let y = fill_x(b, N);
        let acc = fill_x(b, N);
        mod_mul_write_into_zero_acc_karatsuba2(b, &acc, &x, &y, p);
        b.free_vec(&acc);
        b.free_vec(&y);
        b.free_vec(&x);
    })
}

/// Isolated schoolbook_mul_into_addsub forward+inverse pair.
fn bench_schoolbook_addsub_pair(lowq: bool) -> Measured {
    measure(|b| {
        let x = fill_x(b, N);
        let y = fill_x(b, N);
        let tmp_ext = b.alloc_qubits(2 * N);
        if lowq {
            schoolbook_mul_into_addsub_lowq(b, &x, &y, &tmp_ext);
            schoolbook_mul_into_addsub_lowq_inverse(b, &x, &y, &tmp_ext);
        } else {
            schoolbook_mul_into_addsub(b, &x, &y, &tmp_ext);
            schoolbook_mul_into_addsub_inverse(b, &x, &y, &tmp_ext);
        }
        b.free_vec(&tmp_ext);
        b.free_vec(&y);
        b.free_vec(&x);
    })
}

/// Forward-only schoolbook_mul_into_addsub (fast vs lowq).
fn bench_schoolbook_addsub_forward(lowq: bool) -> Measured {
    measure(|b| {
        let x = fill_x(b, N);
        let y = fill_x(b, N);
        let tmp_ext = b.alloc_qubits(2 * N);
        if lowq {
            schoolbook_mul_into_addsub_lowq(b, &x, &y, &tmp_ext);
        } else {
            schoolbook_mul_into_addsub(b, &x, &y, &tmp_ext);
        }
        b.free_vec(&tmp_ext);
        b.free_vec(&y);
        b.free_vec(&x);
    })
}

fn bench_mul_add_karatsuba2() -> Measured {
    let p = SECP256K1_P;
    measure(|b| {
        let x = fill_x(b, N);
        let y = fill_x(b, N);
        let acc = fill_x(b, N);
        mod_mul_add_into_acc_karatsuba2(b, &acc, &x, &y, p);
        b.free_vec(&acc);
        b.free_vec(&y);
        b.free_vec(&x);
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bench_adders() {
        if !enabled() {
            return;
        }
        for n in [64, 128, 256] {
            print_row(&format!("cuccaro_add_fast n={}", n), &bench_cuccaro_add_fast(n));
            print_row(&format!("cuccaro_add_slow n={}", n), &bench_cuccaro_add_slow(n));
            print_row(&format!("cuccaro_sub_fast n={}", n), &bench_cuccaro_sub_fast(n));
        }
    }

    #[test]
    fn bench_muls() {
        if !enabled() {
            return;
        }
        print_row("mul_schoolbook_write_zero", &bench_mul_schoolbook_write());
        print_row("mul_schoolbook_add", &bench_mul_add_schoolbook());
        print_row("mul_karatsuba1_write_zero", &bench_mul_karatsuba_write());
        print_row("mul_karatsuba1_add", &bench_mul_add_karatsuba());
        print_row("mul_karatsuba2_write_zero", &bench_mul_karatsuba2_write());
        print_row("mul_karatsuba2_add", &bench_mul_add_karatsuba2());
        print_row(
            "schoolbook_addsub_fast   (forward+inverse)",
            &bench_schoolbook_addsub_pair(false),
        );
        print_row(
            "schoolbook_addsub_lowq   (forward+inverse)",
            &bench_schoolbook_addsub_pair(true),
        );
        print_row(
            "schoolbook_addsub_fast   (forward only)",
            &bench_schoolbook_addsub_forward(false),
        );
        print_row(
            "schoolbook_addsub_lowq   (forward only)",
            &bench_schoolbook_addsub_forward(true),
        );

        // New lowq Karatsuba variant: should show +Toffoli, -peak_q vs karatsuba1.
        let lowq_add = measure(|b| {
            let p = SECP256K1_P;
            let x = fill_x(b, N);
            let y = fill_x(b, N);
            let acc = fill_x(b, N);
            mod_mul_add_into_acc_karatsuba_lowq(b, &acc, &x, &y, p);
            b.free_vec(&acc);
            b.free_vec(&y);
            b.free_vec(&x);
        });
        print_row("mul_karatsuba1_lowq_add", &lowq_add);
    }
}
