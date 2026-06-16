//! Alg 3: Bezout reconstruction. Reads the garbage bit-vector from
//! `gcd_pack` in REVERSE, runs linear modular updates on (r, s) seeded
//! at (y, 0); ends at (0, y * x^-1 mod q).
//!
//! Per Schrottenloher 2026 (Section 3.2) and Schrottenloher's `ApplyBitVector`
//! (gcd.py:300) / `apply_bitvector` (`gcd_functions.py:149`).

use crate::point_add::trailmix_port::arith::schrottenloher::gcd_compress5::uncompress_classical_5;
use crate::point_add::trailmix_port::arith::schrottenloher::gcd_pack::{DIALOG_M, DIALOG_PACK};
use crate::point_add::trailmix_port::circuit::{Circuit, QReg};

/// Window-size- and vent-parameterized forward apply-bitvector (Horner
/// `tmp = x_orig * y mod q` driven by the dialog tape). `dialog_m = 5`
/// decompresses the base-3 packs per window (low qubit, +204k Toffoli of
/// compression); `dialog_m = 1` reads the RAW 2-bit pairs directly as
/// controls. `vents` is a measurement-vent ancilla budget spent at THIS
/// peak phase: each per-iter mod-double / controlled-add then turns up to
/// `vents` carry-uncompute Toffolis into measurements (-1 Toffoli, +1 peak;
/// the pool is reused across the sequential ops, so peak rises by `vents`,
/// not `vents` per op). `vents = 0` is the peak-safe Cuccaro path.
pub fn apply_bitvector_quantum_secp256k1_m(
    circ: &mut Circuit,
    garbage: &[QReg],
    x_reg: &[QReg],
    y_reg: &[QReg],
    dialog_m: usize,
    vents: usize,
) {
    use crate::point_add::trailmix_port::arith::schrottenloher::gcd_compress5::{
        compress_5iter_refs, compress_5iter_reverse_refs,
    };
    use crate::point_add::trailmix_port::arith::schrottenloher::gcd_pack::{dialog_pack, expected_iterations};
    use crate::point_add::trailmix_port::arith::schrottenloher::pm_prims::{
        controlled_mod_add_pm_secp256k1_vents, mod_double_pm_secp256k1_vents,
    };

    let n = 256usize;
    assert_eq!(x_reg.len(), n + 1, "x_reg must be n+1 = 257 bits");
    assert_eq!(y_reg.len(), n + 1, "y_reg must be n+1 = 257 bits");
    let pack = dialog_pack(dialog_m);
    let iters = expected_iterations(n);
    let garbage_len = iters / dialog_m * pack;
    assert!(
        garbage.len() >= garbage_len,
        "garbage must have at least {garbage_len} bits"
    );

    let prev = circ.push_section("apply_bv");

    // `vents` is the per-call measurement-vent budget spent at this peak
    // phase: it vents the 257-bit register add (threads x_reg) AND, coupled,
    // the ~63-bit +f reduction (materialized, +~125 peak). The two are
    // SEQUENTIAL within each mod-add, so peak rises by max(register_vents,
    // 63+f_vents), not the sum. With the small (m=3 SAT) dialog tape the
    // base peak is low enough that both fit under the budget.
    if dialog_m == 1 {
        for i in (0..iters).rev() {
            let off = 2 * i;
            mod_double_pm_secp256k1_vents(circ, y_reg, vents);
            controlled_mod_add_pm_secp256k1_vents(circ, &garbage[off], x_reg, y_reg, vents);
            for j in 0..=n {
                circ.cswap(&garbage[off + 1], &x_reg[j], &y_reg[j]);
            }
        }
        circ.pop_section(&prev);
        return;
    }

    let b0 = circ.alloc_qreg("apply.b0");
    let b0_and_b1 = circ.alloc_qreg("apply.b0_and_b1");

    // Per-window dialog: decompress one window's `pack`-bit pack into the
    // 2*dialog_m-bit (b0,b0&b1) pair layout ONCE per window (held in
    // num_anc = 2*dialog_m - pack extra ancilla), read each iter's pair via
    // cheap swaps, recompress at window exit. m=5 uses the base-3 arithmetic
    // compressor (10-bit window, 2 anc, borrows x_reg[0..8] dirty); m=3 uses
    // the cheap SAT compressor (6-bit window, 1 anc, no scratch).
    let two_m = 2 * dialog_m;
    let num_anc = two_m - pack;
    let mut win_anc: Vec<QReg> = Vec::new();

    for i in (0..iters).rev() {
        let pack_off = pack * (i / dialog_m);
        let slot = i % dialog_m;

        if slot == dialog_m - 1 {
            win_anc = (0..num_anc).map(|_| circ.alloc_qreg("apply_win")).collect();
            if dialog_m == 5 {
                let w: [&QReg; 10] = std::array::from_fn(|k| {
                    if k < pack {
                        &garbage[pack_off + k]
                    } else {
                        &win_anc[k - pack]
                    }
                });
                let dirty: [&QReg; 8] = std::array::from_fn(|j| &x_reg[j]);
                compress_5iter_reverse_refs(circ, &w, &dirty);
            } else {
                let w: [&QReg; 6] = std::array::from_fn(|k| {
                    if k < pack {
                        &garbage[pack_off + k]
                    } else {
                        &win_anc[k - pack]
                    }
                });
                crate::point_add::trailmix_port::arith::schrottenloher::gcd_compress::compress_3iter_reverse_refs(circ, &w);
            }
        }

        // Extract this iter's pair (b0, b0&b1) from the decompressed window.
        // Window bit k lives in the garbage pack (k < pack) or the freshly
        // decompressed ancilla (k >= pack).
        {
            let p0 = 2 * slot;
            let p1 = 2 * slot + 1;
            let r0: &QReg = if p0 < pack {
                &garbage[pack_off + p0]
            } else {
                &win_anc[p0 - pack]
            };
            circ.swap(&b0, r0);
            let r1: &QReg = if p1 < pack {
                &garbage[pack_off + p1]
            } else {
                &win_anc[p1 - pack]
            };
            circ.swap(&b0_and_b1, r1);
        }

        // y := 2y mod q ; y += b0 * x mod q ; cswap(b0_and_b1, x, y).
        mod_double_pm_secp256k1_vents(circ, y_reg, vents);
        controlled_mod_add_pm_secp256k1_vents(circ, &b0, x_reg, y_reg, vents);
        for j in 0..=n {
            circ.cswap(&b0_and_b1, &x_reg[j], &y_reg[j]);
        }

        // Put the pair back (b0, b0&b1 return to |0>).
        {
            let p0 = 2 * slot;
            let p1 = 2 * slot + 1;
            let r0: &QReg = if p0 < pack {
                &garbage[pack_off + p0]
            } else {
                &win_anc[p0 - pack]
            };
            circ.swap(&b0, r0);
            let r1: &QReg = if p1 < pack {
                &garbage[pack_off + p1]
            } else {
                &win_anc[p1 - pack]
            };
            circ.swap(&b0_and_b1, r1);
        }

        if slot == 0 {
            if dialog_m == 5 {
                let w: [&QReg; 10] = std::array::from_fn(|k| {
                    if k < pack {
                        &garbage[pack_off + k]
                    } else {
                        &win_anc[k - pack]
                    }
                });
                let dirty: [&QReg; 8] = std::array::from_fn(|j| &x_reg[j]);
                compress_5iter_refs(circ, &w, &dirty);
            } else {
                let w: [&QReg; 6] = std::array::from_fn(|k| {
                    if k < pack {
                        &garbage[pack_off + k]
                    } else {
                        &win_anc[k - pack]
                    }
                });
                crate::point_add::trailmix_port::arith::schrottenloher::gcd_compress::compress_3iter_refs(circ, &w);
            }
            for q in std::mem::take(&mut win_anc).into_iter().rev() {
                circ.zero_and_free(q);
            }
        }
    }
    debug_assert!(win_anc.is_empty(), "all apply_bv windows closed");

    // Strict-dealloc settle: the final window's recompress allocates internal
    // compressor ancilla after b0/b0_and_b1's last swap-touch, so the drops
    // below would not be alloc-clean. Both are |0> here, so cx into x_reg[0]
    // with distinct controls is a semantic no-op advancing their touch edge.
    circ.cx(&b0, &x_reg[0]);
    circ.cx(&b0_and_b1, &x_reg[0]);
    drop(b0);
    drop(b0_and_b1);
    circ.pop_section(&prev);
}

/// Window-size- and vent-parameterized inverse apply-bitvector (`tmp =
/// y * x^-1 mod q`). `dialog_m = 1` reads the RAW 2-bit pairs directly;
/// `dialog_m = 5` decompresses the base-3 packs. `vents` is the
/// measurement-vent budget for the per-iter controlled mod-sub at this
/// peak phase (mod-halve's reduction stays unvented). See
/// `apply_bitvector_quantum_secp256k1_m`.
pub fn apply_bitvector_quantum_secp256k1_inv_m(
    circ: &mut Circuit,
    garbage: &[QReg],
    x_reg: &[QReg],
    y_reg: &[QReg],
    dialog_m: usize,
    vents: usize,
) {
    use crate::point_add::trailmix_port::arith::schrottenloher::gcd_compress5::{
        compress_5iter_refs, compress_5iter_reverse_refs,
    };
    use crate::point_add::trailmix_port::arith::schrottenloher::gcd_pack::{dialog_pack, expected_iterations};
    use crate::point_add::trailmix_port::arith::schrottenloher::pm_prims::{
        controlled_mod_sub_pm_secp256k1_vents, mod_halve_pm_secp256k1,
    };
    use crate::point_add::trailmix_port::circuit::Capture;

    let n = 256usize;
    assert_eq!(x_reg.len(), n + 1, "x_reg must be n+1 = 257 bits");
    assert_eq!(y_reg.len(), n + 1, "y_reg must be n+1 = 257 bits");
    let pack = dialog_pack(dialog_m);
    let iters = expected_iterations(n);
    let garbage_len = iters / dialog_m * pack;
    assert!(garbage.len() >= garbage_len);

    let prev = circ.push_section("apply_bv_inv");

    // RAW dialog (dialog_m == 1): pairs are garbage[2i], garbage[2i+1],
    // read forward and used directly as controls (tape preserved).
    if dialog_m == 1 {
        for i in 0..iters {
            let off = 2 * i;
            for j in 0..=n {
                circ.cswap(&garbage[off + 1], &x_reg[j], &y_reg[j]);
            }
            controlled_mod_sub_pm_secp256k1_vents(circ, &garbage[off], x_reg, y_reg, vents);
            mod_halve_pm_secp256k1(circ, y_reg);
        }
        circ.pop_section(&prev);
        return;
    }

    let b0 = circ.alloc_qreg("apply_inv.b0");
    let b0_and_b1 = circ.alloc_qreg("apply_inv.b0_and_b1");

    // LOOP INVARIANT: (u, v) at the start of iter i equals
    // apply_bitvector_classical_reverse(u_pre, v_pre, garbage[..iter i], q)
    // applied to the entry (u_entry, v_entry). Capture (u, v) once
    // at entry and check post-iter against classical predicted state.
    let dbg = std::env::var("CIRC_APPLY_INV_DEBUG").is_ok();
    let q: crate::point_add::trailmix_port::compat::num_bigint::BigUint =
        (crate::point_add::trailmix_port::compat::num_bigint::BigUint::from(1u32) << 256u32) - crate::point_add::trailmix_port::compat::num_bigint::BigUint::from(super::F_SECP256K1);
    let entry_cap: Option<Capture<(crate::point_add::trailmix_port::compat::num_bigint::BigUint, crate::point_add::trailmix_port::compat::num_bigint::BigUint, Vec<u8>)>> = if dbg {
        let x_for: Vec<&QReg> = x_reg[..n].iter().collect();
        let y_for: Vec<&QReg> = y_reg[..n].iter().collect();
        let g_for: Vec<&QReg> = garbage[..garbage_len].iter().collect();
        Some(circ.contract_capture_handle::<(crate::point_add::trailmix_port::compat::num_bigint::BigUint, crate::point_add::trailmix_port::compat::num_bigint::BigUint, Vec<u8>), _>(
            "apply_inv_entry",
            move |view, shot| -> Result<(crate::point_add::trailmix_port::compat::num_bigint::BigUint, crate::point_add::trailmix_port::compat::num_bigint::BigUint, Vec<u8>), String> {
                let read = |regs: &[&QReg]| -> crate::point_add::trailmix_port::compat::num_bigint::BigUint {
                    let mut v = crate::point_add::trailmix_port::compat::num_bigint::BigUint::from(0u32);
                    for (i, q_) in regs.iter().enumerate() {
                        if view.contract_read_bit_shot(q_, shot) {
                            v |= crate::point_add::trailmix_port::compat::num_bigint::BigUint::from(1u32) << i;
                        }
                    }
                    v
                };
                let mut g = vec![0u8; g_for.len().div_ceil(8)];
                for (k, q_) in g_for.iter().enumerate() {
                    if view.contract_read_bit_shot(q_, shot) {
                        g[k / 8] |= 1 << (k % 8);
                    }
                }
                Ok((read(&x_for), read(&y_for), g))
            },
        ))
    } else {
        None
    };

    // Per-window dialog (forward garbage read): window opens at slot 0,
    // closes at slot M-1. See apply_bitvector_quantum_secp256k1.
    let two_m = 2 * dialog_m;
    let num_anc = two_m - pack;
    let mut win_anc: Vec<QReg> = Vec::new();

    for i in 0..iters {
        let pack_off = pack * (i / dialog_m);
        let slot = i % dialog_m;

        if slot == 0 {
            win_anc = (0..num_anc)
                .map(|_| circ.alloc_qreg("apply_inv_win"))
                .collect();
            if dialog_m == 5 {
                let w: [&QReg; 10] = std::array::from_fn(|k| {
                    if k < pack {
                        &garbage[pack_off + k]
                    } else {
                        &win_anc[k - pack]
                    }
                });
                let dirty: [&QReg; 8] = std::array::from_fn(|j| &x_reg[j]);
                compress_5iter_reverse_refs(circ, &w, &dirty);
            } else {
                let w: [&QReg; 6] = std::array::from_fn(|k| {
                    if k < pack {
                        &garbage[pack_off + k]
                    } else {
                        &win_anc[k - pack]
                    }
                });
                crate::point_add::trailmix_port::arith::schrottenloher::gcd_compress::compress_3iter_reverse_refs(circ, &w);
            }
        }

        // Extract this iter's pair from the decompressed window.
        {
            let p0 = 2 * slot;
            let p1 = 2 * slot + 1;
            let r0: &QReg = if p0 < pack {
                &garbage[pack_off + p0]
            } else {
                &win_anc[p0 - pack]
            };
            circ.swap(&b0, r0);
            let r1: &QReg = if p1 < pack {
                &garbage[pack_off + p1]
            } else {
                &win_anc[p1 - pack]
            };
            circ.swap(&b0_and_b1, r1);
        }

        // if b0_and_b1: cswap x, y.
        for j in 0..=n {
            circ.cswap(&b0_and_b1, &x_reg[j], &y_reg[j]);
        }

        // APPROXIMATE controlled mod-sub (Alg 11 mirror of Alg 10):
        // y -= b0 * x mod q. Borrow flag cleanup uses top-K add-overflow
        // (structural mirror of forward mod-add's top-K sub-borrow
        // cleanup), so the approximation behavior matches forward
        // mod-add exactly. ~900 Tof/call vs exact ~2800.
        controlled_mod_sub_pm_secp256k1_vents(circ, &b0, x_reg, y_reg, vents);

        // y *= 2^-1 mod q (approximate halve via flag==a[255] identity).
        mod_halve_pm_secp256k1(circ, y_reg);

        // Put the pair back.
        {
            let p0 = 2 * slot;
            let p1 = 2 * slot + 1;
            let r0: &QReg = if p0 < pack {
                &garbage[pack_off + p0]
            } else {
                &win_anc[p0 - pack]
            };
            circ.swap(&b0, r0);
            let r1: &QReg = if p1 < pack {
                &garbage[pack_off + p1]
            } else {
                &win_anc[p1 - pack]
            };
            circ.swap(&b0_and_b1, r1);
        }

        if slot == dialog_m - 1 {
            if dialog_m == 5 {
                let w: [&QReg; 10] = std::array::from_fn(|k| {
                    if k < pack {
                        &garbage[pack_off + k]
                    } else {
                        &win_anc[k - pack]
                    }
                });
                let dirty: [&QReg; 8] = std::array::from_fn(|j| &x_reg[j]);
                compress_5iter_refs(circ, &w, &dirty);
            } else {
                let w: [&QReg; 6] = std::array::from_fn(|k| {
                    if k < pack {
                        &garbage[pack_off + k]
                    } else {
                        &win_anc[k - pack]
                    }
                });
                crate::point_add::trailmix_port::arith::schrottenloher::gcd_compress::compress_3iter_refs(circ, &w);
            }
            for q in std::mem::take(&mut win_anc).into_iter().rev() {
                circ.zero_and_free(q);
            }
        }

        // INVARIANT check: after iter i, (u, v) should equal the
        // classical apply_reverse-style step on the captured entry
        // state, processing garbage iters 0..=i (forward read).
        if let Some(ref cap) = entry_cap {
            let x_for: Vec<&QReg> = x_reg[..n].iter().collect();
            let y_for: Vec<&QReg> = y_reg[..n].iter().collect();
            let q_c = q.clone();
            let iter_idx = i;
            cap.check(circ, move |captured, view, shot| -> Result<(), String> {
                let (u0, v0, g_at_entry) = captured;
                // Apply classical iter 0..=iter_idx (forward read).
                let mut u: crate::point_add::trailmix_port::compat::num_bigint::BigUint = u0.clone() % &q_c;
                let mut v: crate::point_add::trailmix_port::compat::num_bigint::BigUint = v0.clone() % &q_c;
                let two_inv = crate::point_add::trailmix_port::compat::num_bigint::BigUint::from(2u32)
                    .modpow(&(&q_c - crate::point_add::trailmix_port::compat::num_bigint::BigUint::from(2u32)), &q_c);
                for it in 0..=iter_idx {
                    let p_off = DIALOG_PACK * (it / DIALOG_M);
                    let mut pack = 0u8;
                    for k in 0..DIALOG_PACK {
                        let bit = (g_at_entry[(p_off + k) / 8] >> ((p_off + k) % 8)) & 1;
                        pack |= bit << k;
                    }
                    let decompressed = uncompress_classical_5(pack)
                        .ok_or_else(|| format!("bad pack at iter {it}"))?;
                    let slot = (it % DIALOG_M) * 2;
                    let b0_class = (decompressed >> slot) & 1;
                    let b0_and_b1_class = (decompressed >> (slot + 1)) & 1;
                    if b0_and_b1_class == 1 {
                        std::mem::swap(&mut u, &mut v);
                    }
                    if b0_class == 1 {
                        v = if v >= u { (v - &u) % &q_c } else { (v + &q_c - &u) % &q_c };
                    }
                    v = (&v * &two_inv) % &q_c;
                }
                let read = |regs: &[&QReg]| -> crate::point_add::trailmix_port::compat::num_bigint::BigUint {
                    let mut v = crate::point_add::trailmix_port::compat::num_bigint::BigUint::from(0u32);
                    for (i, q_) in regs.iter().enumerate() {
                        if view.contract_read_bit_shot(q_, shot) {
                            v |= crate::point_add::trailmix_port::compat::num_bigint::BigUint::from(1u32) << i;
                        }
                    }
                    v
                };
                let got_u = read(&x_for) % &q_c;
                let got_v = read(&y_for) % &q_c;
                if got_u != u || got_v != v {
                    return Err(format!(
                        "apply_inv iter {iter_idx} shot {shot}: want=({u:#x}, {v:#x}), got=({got_u:#x}, {got_v:#x})"
                    ));
                }
                Ok(())
            });
        }
    }

    // Strict-dealloc settle (see apply_bitvector_quantum_secp256k1): touch
    // b0/b0_and_b1 past the final swapper5's compressor allocs. |0> here.
    circ.cx(&b0, &x_reg[0]);
    circ.cx(&b0_and_b1, &x_reg[0]);
    drop(b0);
    drop(b0_and_b1);
    drop(entry_cap);
    circ.pop_section(&prev);
}
