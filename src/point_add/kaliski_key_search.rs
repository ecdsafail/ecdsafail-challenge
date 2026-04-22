//! Search over small side-information keys for the hybrid Kaliski-jump moonshot.
//!
//! We already know the local 4-step joint transition family is tiny (125 joint
//! classes for `t=4`) and that `(u_low, v_low, cmp0, cmp1)` gives mean
//! ambiguity ≈ 1.74 and max ambiguity 4.
//!
//! This file brute-forces a small family of candidate side-information features
//! to see whether some equally-cheap key collapses the ambiguity even further.
//!
//! Result so far: on 300 random secp256k1 trajectories (~108k windows),
//! the triple of compare bits `(cmp0, cmp1, cmp2)` is by far the best cheap key:
//!
//!   key = (u_low, v_low, cmp0, cmp1, cmp2)
//!   mean ambiguity ≈ 1.035, max ambiguity = 2.
//!
//! This is the strongest classical signal yet that a practical hybrid batched
//! primitive might only need low bits plus three compare bits.

use std::collections::{BTreeMap, BTreeSet};

use alloy_primitives::U256;
use sha3::digest::{ExtendableOutput, Update, XofReader};

use super::SECP256K1_P;
use super::kaliski_jump::kaliski_step_uv;
use super::test_timeout::{check_deadline, two_min_deadline};

pub struct Sampler {
    reader: Box<dyn XofReader>,
    p: U256,
}

impl Sampler {
    pub fn new(seed: &[u8], p: U256) -> Self {
        let mut hasher = sha3::Shake128::default();
        hasher.update(seed);
        Self { reader: Box::new(hasher.finalize_xof()), p }
    }

    pub fn next(&mut self) -> U256 {
        loop {
            let mut buf = [0u8; 32];
            self.reader.read(&mut buf);
            let x = U256::from_le_slice(&buf);
            if x < self.p && !x.is_zero() {
                return x;
            }
        }
    }
}

#[derive(Clone, Debug)]
struct Row {
    low_u: u16,
    low_v: u16,
    cmp0: u8,
    cmp1: u8,
    cmp2: u8,
    u1_1: u8,
    v1_1: u8,
    u1_2: u8,
    v1_2: u8,
    u1_3: u8,
    v1_3: u8,
    u2_1: u8,
    v2_1: u8,
    u2_2: u8,
    v2_2: u8,
    seq_code: u16,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Feat {
    Cmp0, Cmp1, Cmp2,
    U1_1, V1_1,
    U1_2, V1_2,
    U1_3, V1_3,
    U2_1, V2_1,
    U2_2, V2_2,
}

fn feat_val(r: &Row, f: Feat) -> u8 {
    match f {
        Feat::Cmp0 => r.cmp0,
        Feat::Cmp1 => r.cmp1,
        Feat::Cmp2 => r.cmp2,
        Feat::U1_1 => r.u1_1,
        Feat::V1_1 => r.v1_1,
        Feat::U1_2 => r.u1_2,
        Feat::V1_2 => r.v1_2,
        Feat::U1_3 => r.u1_3,
        Feat::V1_3 => r.v1_3,
        Feat::U2_1 => r.u2_1,
        Feat::V2_1 => r.v2_1,
        Feat::U2_2 => r.u2_2,
        Feat::V2_2 => r.v2_2,
    }
}

fn build_rows(seed: &[u8], n_inputs: usize, w: usize, t: usize) -> Vec<Row> {
    let deadline = two_min_deadline();
    let mut sampler = Sampler::new(seed, SECP256K1_P);
    let mask = if w == 16 { U256::from(0xFFFFu64) } else { (U256::from(1u64) << w).wrapping_sub(U256::from(1u64)) };
    let mut rows = Vec::new();
    for input_idx in 0..n_inputs {
        if (input_idx & 31) == 0 { check_deadline(deadline, "kaliski_key_search::build_rows"); }
        let mut u = SECP256K1_P;
        let mut v = sampler.next();
        for _ in 0..742 {
            if v.is_zero() { break; }
            let mut us = vec![u];
            let mut vs = vec![v];
            let mut seq_code: u16 = 0;
            let mut uu = u;
            let mut vv = v;
            for i in 0..t {
                if vv.is_zero() { break; }
                let (nu, nv, kc) = kaliski_step_uv(uu, vv);
                seq_code |= (match kc {
                    super::kaliski_jump::KCase::UEven => 0u16,
                    super::kaliski_jump::KCase::VEven => 1u16,
                    super::kaliski_jump::KCase::UGtV  => 2u16,
                    super::kaliski_jump::KCase::VGtU  => 3u16,
                }) << (2 * i);
                uu = nu; vv = nv;
                us.push(uu); vs.push(vv);
            }
            let cmp = |i: usize| -> u8 {
                if i < us.len() && i < vs.len() { (us[i] > vs[i]) as u8 } else { 0 }
            };
            let ub = |i: usize, bits: usize| -> u8 {
                if i < us.len() { (us[i] & U256::from((1u64 << bits) - 1)).to::<u8>() } else { 0 }
            };
            let vb = |i: usize, bits: usize| -> u8 {
                if i < vs.len() { (vs[i] & U256::from((1u64 << bits) - 1)).to::<u8>() } else { 0 }
            };
            rows.push(Row {
                low_u: (u & mask).to::<u16>(),
                low_v: (v & mask).to::<u16>(),
                cmp0: cmp(0),
                cmp1: cmp(1),
                cmp2: cmp(2),
                u1_1: ub(1,1), v1_1: vb(1,1),
                u1_2: ub(1,2), v1_2: vb(1,2),
                u1_3: ub(1,3), v1_3: vb(1,3),
                u2_1: ub(2,1), v2_1: vb(2,1),
                u2_2: ub(2,2), v2_2: vb(2,2),
                seq_code,
            });
            let (u1, v1, _kc) = kaliski_step_uv(u, v);
            u = u1; v = v1;
        }
    }
    rows
}

#[derive(Debug, Clone)]
pub struct ComboStats {
    pub combo: Vec<Feat>,
    pub classes: usize,
    pub mean_seq_per_class: f64,
    pub max_seq_per_class: usize,
    pub singleton_classes: usize,
}

fn eval_combo(rows: &[Row], combo: &[Feat]) -> ComboStats {
    let mut map: BTreeMap<(u16, u16, Vec<u8>), BTreeSet<u16>> = BTreeMap::new();
    for r in rows {
        let extra: Vec<u8> = combo.iter().map(|&f| feat_val(r, f)).collect();
        map.entry((r.low_u, r.low_v, extra)).or_default().insert(r.seq_code);
    }
    let classes = map.len();
    let mut total = 0usize;
    let mut maxc = 0usize;
    let mut singles = 0usize;
    for seqs in map.values() {
        let c = seqs.len();
        total += c;
        if c > maxc { maxc = c; }
        if c == 1 { singles += 1; }
    }
    ComboStats {
        combo: combo.to_vec(),
        classes,
        mean_seq_per_class: total as f64 / classes as f64,
        max_seq_per_class: maxc,
        singleton_classes: singles,
    }
}

pub fn search_feature_combos(seed: &[u8], n_inputs: usize, w: usize, t: usize) -> Vec<ComboStats> {
    let deadline = two_min_deadline();
    let rows = build_rows(seed, n_inputs, w, t);
    let feats = [
        Feat::Cmp0, Feat::Cmp1, Feat::Cmp2,
        Feat::U1_1, Feat::V1_1,
        Feat::U1_2, Feat::V1_2,
        Feat::U1_3, Feat::V1_3,
        Feat::U2_1, Feat::V2_1,
        Feat::U2_2, Feat::V2_2,
    ];
    let mut out = Vec::new();

    fn rec(
        cur: &mut Vec<Feat>,
        idx: usize,
        left: usize,
        feats: &[Feat],
        rows: &[Row],
        out: &mut Vec<ComboStats>,
        deadline: std::time::Instant,
    ) {
        if left == 0 {
            out.push(eval_combo(rows, cur));
            return;
        }
        for i in idx..=feats.len() - left {
            if (out.len() & 255) == 0 { check_deadline(deadline, "kaliski_key_search::search_feature_combos"); }
            cur.push(feats[i]);
            rec(cur, i + 1, left - 1, feats, rows, out, deadline);
            cur.pop();
        }
    }

    for k in 1..=4 {
        let mut cur = Vec::new();
        rec(&mut cur, 0, k, &feats, &rows, &mut out, deadline);
    }

    out.sort_by(|a, b| {
        a.max_seq_per_class.cmp(&b.max_seq_per_class)
            .then_with(|| a.mean_seq_per_class.partial_cmp(&b.mean_seq_per_class).unwrap())
            .then_with(|| b.singleton_classes.cmp(&a.singleton_classes))
            .then_with(|| a.combo.len().cmp(&b.combo.len()))
    });
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn feature_search_top_results() {
        // Keep under 2 minutes by using 300 sampled trajectories. This still
        // gives >100k windows and was enough to discover the best combo.
        let res = search_feature_combos(b"kaliski-key-combo-seed-v1", 300, 8, 4);
        eprintln!("=== Kaliski feature-combo search (w=8,t=4, 300 inputs) ===");
        for item in res.iter().take(20) {
            eprintln!("combo={:?} classes={} mean={:.3} max={} singletons={}",
                item.combo, item.classes, item.mean_seq_per_class, item.max_seq_per_class, item.singleton_classes);
        }
        eprintln!("==========================================================");
        // The best known cheap key should show up at the top.
        assert!(!res.is_empty());
        assert!(res[0].max_seq_per_class <= 2);
    }
}
