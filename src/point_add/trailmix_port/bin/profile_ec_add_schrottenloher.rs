//! Per-section Toffoli breakdown of the Schrottenloher EC-add, with the
//! base-3 compression cost isolated via the `c5pack` markers inside
//! `compress_5iter_refs` / `_reverse_refs`. Tells us how much Toffoli the
//! base-3 packing costs (i.e. how much a raw 2-bit-per-pair dialog would
//! save) and at what qubit peak the current packed version sits.

use alloy_primitives::U256;
use rand::RngCore;
use trailmix::arith::schrottenloher::gcd_pack::u_padding;
use trailmix::arith::schrottenloher::pointadd::{
    ec_add_inplace_schrottenloher_jump_cfg, ec_add_inplace_schrottenloher_jump_lowtof_secp256k1,
    ec_add_inplace_schrottenloher_jump_secp256k1, ec_add_inplace_schrottenloher_secp256k1_m,
};
use trailmix::circuit::{Cbit, Circuit, QReg};
use trailmix::ec::point_add::ec_add_inplace_shrunken_pz;
use zkp_ecc_lib::WeierstrassEllipticCurve;

fn secp256k1() -> WeierstrassEllipticCurve {
    WeierstrassEllipticCurve {
        modulus: U256::from_str_radix(
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F",
            16,
        )
        .unwrap(),
        a: U256::from(0u64),
        b: U256::from(7u64),
        gx: U256::from_str_radix(
            "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798",
            16,
        )
        .unwrap(),
        gy: U256::from_str_radix(
            "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8",
            16,
        )
        .unwrap(),
        order: U256::from_str_radix(
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141",
            16,
        )
        .unwrap(),
    }
}

fn rand_case(curve: &WeierstrassEllipticCurve, rng: &mut impl RngCore) -> (U256, U256, U256, U256) {
    loop {
        let scalar = |rng: &mut dyn RngCore| {
            U256::from(rng.next_u64())
                ^ (U256::from(rng.next_u64()) << 64)
                ^ (U256::from(rng.next_u64()) << 128)
                ^ (U256::from(rng.next_u64()) << 192)
        };
        let (s_p, s_q) = (scalar(rng), scalar(rng));
        if s_p == U256::ZERO || s_q == U256::ZERO || s_p == s_q {
            continue;
        }
        let p = curve.mul(curve.gx, curve.gy, s_p);
        let q = curve.mul(curve.gx, curve.gy, s_q);
        if p.0 == q.0 {
            continue;
        }
        return (p.0, p.1, q.0, q.1);
    }
}

fn main() {
    const SHOTS: usize = 64;
    let curve = secp256k1();
    let mut rng = rand::thread_rng();
    let n = 256usize;
    let total = n + u_padding(n);

    let arg1 = std::env::args().nth(1).unwrap_or_else(|| "1".into());
    let arg2: usize = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(0);

    let mut circ = Circuit::new();
    circ.set_max_qubit_peak(1600);
    let mut outs: Vec<QReg> = Vec::new();

    if arg1 == "shrunkenpz" {
        eprintln!("[profile] shrunken-PZ variant");
        let mut tx: Vec<QReg> = (0..n)
            .map(|i| circ.alloc_qreg(&format!("tx[{i}]")))
            .collect();
        let mut ty: Vec<QReg> = (0..n)
            .map(|i| circ.alloc_qreg(&format!("ty[{i}]")))
            .collect();
        let ox: Vec<Cbit> = (0..n).map(|_| circ.alloc_input_bit()).collect();
        let oy: Vec<Cbit> = (0..n).map(|_| circ.alloc_input_bit()).collect();
        for shot in 0..SHOTS {
            let (px, py, qx, qy) = rand_case(&curve, &mut rng);
            circ.sim_load_reg_bytes_shot(&tx[..n], &px.to_le_bytes::<32>(), shot);
            circ.sim_load_reg_bytes_shot(&ty[..n], &py.to_le_bytes::<32>(), shot);
            circ.sim_load_bits_bytes_shot(&ox, &qx.to_le_bytes::<32>(), shot);
            circ.sim_load_bits_bytes_shot(&oy, &qy.to_le_bytes::<32>(), shot);
        }
        ec_add_inplace_shrunken_pz(&mut circ, &mut tx, &mut ty, &ox, &oy);
        outs.extend(std::mem::take(&mut tx));
        outs.extend(std::mem::take(&mut ty));
    } else {
        let mut x2: Vec<QReg> = (0..total)
            .map(|i| circ.alloc_qreg(&format!("x2[{i}]")))
            .collect();
        let y2: Vec<QReg> = {
            let mut v: Vec<QReg> = (0..n)
                .map(|i| circ.alloc_qreg(&format!("y2[{i}]")))
                .collect();
            v.push(circ.alloc_qreg("y2_anc"));
            v
        };
        let ox: Vec<Cbit> = (0..n).map(|_| circ.alloc_input_bit()).collect();
        let oy: Vec<Cbit> = (0..n).map(|_| circ.alloc_input_bit()).collect();

        for shot in 0..SHOTS {
            let (px, py, qx, qy) = rand_case(&curve, &mut rng);
            circ.sim_load_reg_bytes_shot(&x2[..n], &px.to_le_bytes::<32>(), shot);
            circ.sim_load_reg_bytes_shot(&y2[..n], &py.to_le_bytes::<32>(), shot);
            circ.sim_load_bits_bytes_shot(&ox, &qx.to_le_bytes::<32>(), shot);
            circ.sim_load_bits_bytes_shot(&oy, &qy.to_le_bytes::<32>(), shot);
        }

        // arg1 dialog_m: 1 = RAW (higher-qubit/lower-tof), 5 = packed (SP1-valid).
        // arg2 vents: apply_bv-peak measurement-vent budget (0 = peak-safe).
        if arg1 == "jump" {
            let jump = if arg2 == 0 { 2 } else { arg2 };
            eprintln!("[profile] JUMP low-qubit variant, jump = {jump}");
            ec_add_inplace_schrottenloher_jump_secp256k1(&mut circ, &mut x2, &y2, &ox, &oy, jump);
        } else if arg1 == "jumplowtof" {
            let jump = if arg2 == 0 { 2 } else { arg2 };
            eprintln!("[profile] JUMP low-tof variant (coupled venting), jump = {jump}");
            ec_add_inplace_schrottenloher_jump_lowtof_secp256k1(
                &mut circ, &mut x2, &y2, &ox, &oy, jump,
            );
        } else if arg1 == "jumpraw" {
            let jump = if arg2 == 0 { 2 } else { arg2 };
            eprintln!("[profile] JUMP raw (UNPACKED, no vent) variant, jump = {jump}");
            ec_add_inplace_schrottenloher_jump_cfg(
                &mut circ, &mut x2, &y2, &ox, &oy, jump, 0, 0, false, false,
            );
        } else {
            let dialog_m: usize = arg1.parse().unwrap_or(1);
            let vents = arg2;
            eprintln!("[profile] dialog_m = {dialog_m}, vents = {vents}");
            ec_add_inplace_schrottenloher_secp256k1_m(
                &mut circ, &mut x2, &y2, &ox, &oy, dialog_m, vents,
            );
        }
        outs.extend(std::mem::take(&mut x2));
        outs.extend(y2);
    }

    let total_tof = (circ.ccx_emitted + circ.ccz_emitted) as u64;
    let total_ops = circ.total_ops();
    let peak = circ.peak_qubits;

    // Sum tof for every section path containing the c5pack marker.
    let mut c5_tof: u64 = 0;
    let mut c5_ops: u64 = 0;
    for (path, count) in &circ.executed_toffoli_by_section {
        if path.contains("c5pack") {
            c5_tof += *count / (SHOTS as u64);
        }
    }
    for (path, count) in &circ.executed_ops_by_section {
        if path.contains("c5pack") {
            c5_ops += *count;
        }
    }

    println!("=== Schrottenloher EC-add (n=256) ===");
    println!("  total ops   : {total_ops}");
    println!("  total tof   : {total_tof}");
    println!("  peak qubits : {peak}");
    println!("  --- base-3 compression (c5pack) ---");
    println!("  c5pack ops  : {c5_ops}");
    println!(
        "  c5pack tof  : {c5_tof}   ({:.1}% of total)",
        100.0 * c5_tof as f64 / total_tof.max(1) as f64
    );
    println!("  --- projected raw-dialog (no c5pack) ---");
    println!(
        "  tof if raw  : ~{}   (saves {c5_tof})",
        total_tof.saturating_sub(c5_tof)
    );

    // Top-level roll-up (by 1st `/`-segment after ec_add_schr/...).
    let mut per_path_tof: std::collections::HashMap<String, u64> = std::collections::HashMap::new();
    for (path, count) in &circ.executed_toffoli_by_section {
        per_path_tof.insert(path.clone(), *count / (SHOTS as u64));
    }
    // Roll up by the LAST `/`-segment (the leaf op kind: csub, cswap,
    // mod_double, apply_bv, rev_cadd, c5pack, ...) — these are the
    // primitive-level buckets we want to compare.
    let mut leaf_tof: std::collections::BTreeMap<String, u64> = std::collections::BTreeMap::new();
    let strip_iter = |s: &str| -> String {
        match s.rfind('_') {
            Some(idx) if idx + 1 < s.len() && s[idx + 1..].chars().all(|c| c.is_ascii_digit()) => {
                s[..idx].to_string()
            }
            _ => s.to_string(),
        }
    };
    // Group key: leaf section, but controlled-add leaves keyed by
    // "<parent>/hybrid_cadd" so apply_bv vs GCD vs square stay separate.
    let key_of = |path: &str| -> String {
        let segs: Vec<&str> = path.split('/').collect();
        if segs.last() == Some(&"hybrid_cadd") && segs.len() >= 2 {
            format!("{}/hybrid_cadd", strip_iter(segs[segs.len() - 2]))
        } else {
            strip_iter(segs.last().unwrap_or(&path))
        }
    };
    for (path, tof) in &per_path_tof {
        *leaf_tof.entry(key_of(path)).or_insert(0) += tof;
    }
    // Per-group LOCAL peak occupancy = max section_peak over paths in the
    // group. Headroom = global_peak - local_peak = ancilla a leaf's ops
    // could use (e.g. an adder saves ~min(headroom, width) Toffoli/call)
    // WITHOUT raising the global peak.
    let mut leaf_localpeak: std::collections::HashMap<String, u32> =
        std::collections::HashMap::new();
    for (path, lp) in &circ.section_peak {
        let e = leaf_localpeak.entry(key_of(path)).or_insert(0);
        if *lp > *e {
            *e = *lp;
        }
    }
    let mut rows: Vec<(String, u64)> = leaf_tof.into_iter().collect();
    rows.sort_by(|a, b| b.1.cmp(&a.1));
    // Which leaves are ADDERS (saveable ~min(headroom,width) tof/call via
    // clean-ancilla / measure-uncompute) vs structural (cswap, etc.).
    let is_adder =
        |k: &str| k.contains("hybrid_cadd") || k.contains("gidney_cadd") || k.contains("cmp");
    println!("  --- expensive leaves: tof / local-peak / HEADROOM (global={peak}) ---");
    println!(
        "  {:<26} {:>10} {:>6} {:>10} {:>8} {}",
        "leaf", "tof", "tof%", "local_pk", "headroom", "kind"
    );
    for (k, tof) in rows.iter().take(20) {
        let lp = *leaf_localpeak.get(k).unwrap_or(&0);
        let head = (peak as i64) - (lp as i64);
        let kind = if is_adder(k) { "ADDER" } else { "struct" };
        println!(
            "  {:<26} {:>10} {:>5.1}% {:>10} {:>8} {}",
            k,
            tof,
            100.0 * *tof as f64 / total_tof.max(1) as f64,
            lp,
            head,
            kind
        );
    }

    let mut exact_peaks: Vec<(String, u32)> = circ
        .section_peak
        .iter()
        .map(|(section, live)| (section.clone(), *live))
        .collect();
    exact_peaks.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    println!("  --- exact section peaks ---");
    println!("  {:>10}  section", "local_pk");
    for (section, live) in exact_peaks.iter().take(25) {
        println!("  {:>10}  {section}", live);
    }

    let normalize_tag = |tag: &str| -> String {
        let mut out = String::new();
        let mut chars = tag.chars().peekable();
        while let Some(ch) = chars.next() {
            if ch == '[' {
                let mut saw_digit = false;
                while let Some(next) = chars.peek().copied() {
                    if next.is_ascii_digit() {
                        saw_digit = true;
                        chars.next();
                    } else {
                        break;
                    }
                }
                if saw_digit && chars.peek() == Some(&']') {
                    chars.next();
                    out.push_str("[]");
                } else {
                    out.push(ch);
                }
            } else {
                out.push(ch);
            }
        }
        out
    };
    let mut peak_tags = std::collections::BTreeMap::<String, usize>::new();
    for tag in &circ.peak_live_tags {
        *peak_tags.entry(normalize_tag(tag)).or_default() += 1;
    }
    let mut peak_tag_rows: Vec<(String, usize)> = peak_tags.into_iter().collect();
    peak_tag_rows.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    println!("  --- peak live tags ---");
    println!("  {:>10}  tag", "qubits");
    for (tag, count) in peak_tag_rows.iter().take(40) {
        println!("  {:>10}  {tag}", count);
    }

    let _ = circ.destroy_sim(outs);
}
