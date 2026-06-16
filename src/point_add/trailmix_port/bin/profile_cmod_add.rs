//! Measure single-call Tof cost of controlled_mod_add_pm_secp256k1
//! and mod_double_pm_secp256k1.

use trailmix::arith::schrottenloher::pm_prims::{
    controlled_mod_add_pm_secp256k1, mod_double_pm_secp256k1,
};
use trailmix::circuit::Circuit;

fn main() {
    let n = 256usize;
    let test = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "cmadd".to_string());

    let mut circ = Circuit::new();
    let ctrl = circ.alloc_qreg("ctrl");
    let x = circ.alloc_qreg_bits("x", n + 1);
    let y = circ.alloc_qreg_bits("y", n + 1);
    circ.x(&ctrl);
    for i in 0..n {
        circ.x(&x[i]);
    }
    for i in 0..n {
        circ.x(&y[i]);
    }

    let ops0 = circ.ops.len();
    let ccx0 = circ.ccx_emitted;
    let ccz0 = circ.ccz_emitted;

    match test.as_str() {
        "cmadd" => controlled_mod_add_pm_secp256k1(&mut circ, &ctrl, &x, &y),
        "mdbl" => mod_double_pm_secp256k1(&mut circ, &y),
        _ => panic!("test: cmadd | mdbl"),
    }

    let ops = circ.ops.len() - ops0;
    let tof = (circ.ccx_emitted - ccx0) + (circ.ccz_emitted - ccz0);
    let peak = circ.peak_qubits;
    println!("{:8} ops={:6} tof={:6} peak={}", test, ops, tof, peak);
    std::process::exit(0);
}
