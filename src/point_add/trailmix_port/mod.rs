pub mod circuit;
pub mod mod_arith;
pub mod rfold_mbu;

pub mod arith {
    pub mod compare;
    pub mod const_add;
    pub mod cuccaro;
    pub mod gidney_const_adder;
    pub mod khattar_gidney;
    pub mod mcx;
    pub mod qshift_sub;
    pub mod ripple_add;
    pub mod shift;

    pub mod rfold_mbu {
        pub use crate::point_add::trailmix_port::rfold_mbu::*;
    }
}

pub mod inversion {
    pub mod shrunken_pz_primitives;
    pub mod shrunken_pz_schedule;
    pub mod shrunken_pz_state_machine;
}

pub mod ec {
    pub mod point_add;
}

use alloy_primitives::U256;
use sha3::digest::{ExtendableOutput, XofReader};

use crate::circuit::{Op, OperationType, QubitId};
use crate::weierstrass_elliptic_curve::WeierstrassEllipticCurve;

const TRAILMIX_TAIL_NONCE_BITS: u32 = 48;
const TRAILMIX_NUM_TESTS: usize = 9024;

pub mod tracker {
    pub mod ghost {
        pub use crate::point_add::trailmix_port::circuit::Ghost;
    }
}

pub mod num_bigint {
    use std::fmt;
    use std::ops::{Add, BitAnd, BitOrAssign, Div, Mul, Rem, Shl, Shr, Sub};

    #[derive(Clone, Default, Debug, Eq, PartialEq, Ord, PartialOrd)]
    pub struct BigUint;

    impl BigUint {
        pub fn from_bytes_le(_bytes: &[u8]) -> Self {
            Self
        }

        pub fn to_bytes_le(&self) -> Vec<u8> {
            Vec::new()
        }
    }

    impl fmt::Display for BigUint {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            f.write_str("0")
        }
    }

    impl fmt::LowerHex for BigUint {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            if f.alternate() {
                f.write_str("0x0")
            } else {
                f.write_str("0")
            }
        }
    }

    impl From<u32> for BigUint {
        fn from(_value: u32) -> Self {
            Self
        }
    }

    impl From<u64> for BigUint {
        fn from(_value: u64) -> Self {
            Self
        }
    }

    impl Add for BigUint {
        type Output = BigUint;
        fn add(self, _rhs: BigUint) -> BigUint {
            BigUint
        }
    }

    impl Add<&BigUint> for BigUint {
        type Output = BigUint;
        fn add(self, _rhs: &BigUint) -> BigUint {
            BigUint
        }
    }

    impl Add<BigUint> for &BigUint {
        type Output = BigUint;
        fn add(self, _rhs: BigUint) -> BigUint {
            BigUint
        }
    }

    impl Add<&BigUint> for &BigUint {
        type Output = BigUint;
        fn add(self, _rhs: &BigUint) -> BigUint {
            BigUint
        }
    }

    impl Add<u32> for &BigUint {
        type Output = BigUint;
        fn add(self, _rhs: u32) -> BigUint {
            BigUint
        }
    }

    impl Add<u32> for BigUint {
        type Output = BigUint;
        fn add(self, _rhs: u32) -> BigUint {
            BigUint
        }
    }

    impl Sub for BigUint {
        type Output = BigUint;
        fn sub(self, _rhs: BigUint) -> BigUint {
            BigUint
        }
    }

    impl Sub<&BigUint> for BigUint {
        type Output = BigUint;
        fn sub(self, _rhs: &BigUint) -> BigUint {
            BigUint
        }
    }

    impl Sub<BigUint> for &BigUint {
        type Output = BigUint;
        fn sub(self, _rhs: BigUint) -> BigUint {
            BigUint
        }
    }

    impl Sub<&BigUint> for &BigUint {
        type Output = BigUint;
        fn sub(self, _rhs: &BigUint) -> BigUint {
            BigUint
        }
    }

    impl Mul for BigUint {
        type Output = BigUint;
        fn mul(self, _rhs: BigUint) -> BigUint {
            BigUint
        }
    }

    impl Mul<BigUint> for &BigUint {
        type Output = BigUint;
        fn mul(self, _rhs: BigUint) -> BigUint {
            BigUint
        }
    }

    impl Mul<&BigUint> for &BigUint {
        type Output = BigUint;
        fn mul(self, _rhs: &BigUint) -> BigUint {
            BigUint
        }
    }

    impl Rem<&BigUint> for BigUint {
        type Output = BigUint;
        fn rem(self, _rhs: &BigUint) -> BigUint {
            BigUint
        }
    }

    impl Rem<BigUint> for BigUint {
        type Output = BigUint;
        fn rem(self, _rhs: BigUint) -> BigUint {
            BigUint
        }
    }

    impl Rem<&BigUint> for &BigUint {
        type Output = BigUint;
        fn rem(self, _rhs: &BigUint) -> BigUint {
            BigUint
        }
    }

    impl Div for BigUint {
        type Output = BigUint;
        fn div(self, _rhs: BigUint) -> BigUint {
            BigUint
        }
    }

    impl BitAnd<&BigUint> for BigUint {
        type Output = BigUint;
        fn bitand(self, _rhs: &BigUint) -> BigUint {
            BigUint
        }
    }

    impl BitAnd<&BigUint> for &BigUint {
        type Output = BigUint;
        fn bitand(self, _rhs: &BigUint) -> BigUint {
            BigUint
        }
    }

    impl Shl<usize> for BigUint {
        type Output = BigUint;
        fn shl(self, _rhs: usize) -> BigUint {
            BigUint
        }
    }

    impl Shl<u32> for BigUint {
        type Output = BigUint;
        fn shl(self, _rhs: u32) -> BigUint {
            BigUint
        }
    }

    impl Shr<u32> for BigUint {
        type Output = BigUint;
        fn shr(self, _rhs: u32) -> BigUint {
            BigUint
        }
    }

    impl BitOrAssign<BigUint> for BigUint {
        fn bitor_assign(&mut self, _rhs: BigUint) {}
    }
}

fn set_default_env(name: &str, value: &str) {
    if std::env::var_os(name).is_none() {
        std::env::set_var(name, value);
    }
}

fn configure_sub1000_trailmix_route() {
    set_default_env("TRAILMIX_THIN_SCHEDULE", "1");
    set_default_env("TRAILMIX_THIN_SEED", "278");
    set_default_env("TRAILMIX_THIN_MARGIN", "0");
    set_default_env("TRAILMIX_THIN_VALIDATE", "500000");
    set_default_env("TRAILMIX_COUNTER_W", "8");
    set_default_env("TRAILMIX_Q_CAP", "19");
    set_default_env("TRAILMIX_Q_CAP_STEP_SET", "6:20,7:20,8:20,9:20,10:20,11:20,12:20,13:20,14:20,15:20,16:20,17:20,18:20,19:20,20:20,21:20,22:20,23:20,24:20,25:20,26:20,27:20,114:20,115:20,116:20,117:20,118:20,119:20,120:20,121:20,122:20,123:20,124:20,125:20,126:20,127:20,128:20,129:20,130:20,131:20,132:20,133:20,134:20,135:20,136:20,137:20,138:20,139:20,140:20,141:20,142:20,143:20,144:20,145:20,374:20,375:20,376:20,377:20,378:20,379:20,380:20,381:20,382:20,383:20,384:20,385:20,386:20,387:20,388:20,389:20,390:20,391:20,392:20,393:20,394:20,395:20,396:20,397:20,398:20,399:20,400:20,401:20,402:25,403:25,404:25,405:25,406:25,407:25,408:25,409:25,410:25,411:25,412:25,413:25,414:25,415:25,416:25,417:25,418:25,419:25,420:25,421:25,422:25,423:25");
    set_default_env("TRAILMIX_THIN_WIDTH_SET", "75:2:68,78:2:68,183:0:173,188:1:171,189:1:171,190:0:171,190:1:171,193:1:169,194:0:169,199:0:166,199:1:166,200:0:166,201:0:165,208:1:162,220:0:156,228:0:152,229:0:151,237:0:147,239:0:146,239:1:147,240:0:147,241:0:145,241:1:146,242:0:146,242:1:145,243:0:145,243:1:145,244:0:145,244:1:145,245:0:144,245:1:145,246:0:145,246:1:144,247:0:144,247:1:144,248:0:144,248:1:144,249:0:141,249:1:144,250:0:144,250:1:143,251:0:142,251:1:142,252:0:141,252:1:141,253:1:141,254:1:141,255:1:141,256:0:141,256:1:141,257:0:139,257:1:139,258:1:137,259:0:135,259:1:137,260:0:137,260:1:136,261:0:135,261:1:135,262:0:134,262:1:135,263:0:134,263:1:135,264:0:135,265:1:134,266:0:134,266:1:134,267:0:133,267:1:133,268:1:132,269:1:132,270:0:132,273:0:129,274:1:129,275:1:129,276:0:129,276:1:129,277:0:128,279:1:127,280:0:127,280:1:126,281:0:125,293:0:119,296:1:119,297:1:119,298:0:119,298:1:118,299:0:117,299:1:116,300:0:115,300:1:116,301:0:115,301:1:116,302:0:116,302:1:115,303:0:115,303:1:115,304:0:115,304:1:115,305:0:114,305:1:115,306:0:115,306:1:114,307:0:112,307:1:114,308:0:114,308:1:113,309:0:113,309:1:112,310:0:112,310:1:112,311:0:112,311:1:112,312:0:112,312:1:112,313:0:110,313:1:112,314:0:112,314:1:111,315:0:111,315:1:110,316:0:109,316:1:110,317:0:109,317:1:110,318:0:110,318:1:110,319:0:108,319:1:109,320:0:109,320:1:109,321:1:108,322:0:108,322:1:108,323:1:105,324:0:105,330:0:101,332:1:101,334:0:98,335:0:98,336:1:98,337:1:98,338:0:98,338:1:98,339:1:96,340:0:96,340:1:95,341:0:93,341:1:93,342:0:93,342:1:93,343:1:93,344:0:93,344:1:93,402:4:25,403:4:25,404:4:25,405:4:25,406:4:25,407:4:25,408:4:25,409:4:25,410:4:25,411:4:25,412:4:25,413:4:25,414:4:25,415:0:55,415:4:25,416:4:25,417:1:55,417:4:25,418:0:55,418:1:55,418:4:25,419:4:25,420:4:25,421:4:25,422:4:25,423:4:25");
    set_default_env("TRAILMIX_SROT_W", "5");
    set_default_env("TRAILMIX_TAIL_NONCE", "50");
}

#[derive(Clone, Debug, Default)]
struct TrailMixSupportReport {
    accepted_shots: usize,
    miss_factors: usize,
    repair_entries: usize,
    first_miss: Option<(usize, &'static str, usize)>,
    qcap_miss_factors: usize,
    qcap_entries: usize,
    first_qcap_miss: Option<(usize, &'static str, usize)>,
}

fn env_usize(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(default)
}

fn env_u64(name: &str, default: u64) -> u64 {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(default)
}

fn secp256k1() -> WeierstrassEllipticCurve {
    WeierstrassEllipticCurve {
        modulus: U256::from_str_radix(
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F",
            16,
        )
        .unwrap(),
        a: U256::from(0),
        b: U256::from(7),
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

fn sub_mod_p(a: U256, b: U256, p: U256) -> U256 {
    if a >= b {
        a - b
    } else {
        p - (b - a)
    }
}

fn support_report_for_xof(
    mut xof: sha3::Shake256Reader,
    target_draws: usize,
) -> TrailMixSupportReport {
    support_report_for_xof_limited(&mut xof, target_draws, None)
}

fn support_report_for_xof_limited(
    xof: &mut sha3::Shake256Reader,
    target_draws: usize,
    max_misses: Option<usize>,
) -> TrailMixSupportReport {
    let curve = secp256k1();
    let mut report = TrailMixSupportReport::default();
    let early_limit = max_misses.unwrap_or(usize::MAX);
    let qcap_check = std::env::var("TRAILMIX_SUPPORT_Q_CAP_CHECK")
        .ok()
        .and_then(|s| s.parse::<usize>().ok());
    let trace_qcap_steps = std::env::var("TRAILMIX_SUPPORT_Q_CAP_TRACE_STEPS")
        .ok()
        .as_deref()
        == Some("1");
    let trace_thin_sites = std::env::var("TRAILMIX_SUPPORT_THIN_TRACE_REPAIRS")
        .ok()
        .as_deref()
        == Some("1");
    for draw in 0..target_draws {
        let mut rb = [[0u8; 32]; 2];
        xof.read(&mut rb[0]);
        xof.read(&mut rb[1]);
        let k1 = U256::from_le_bytes(rb[0]);
        let k2 = U256::from_le_bytes(rb[1]);
        let t = curve.mul(curve.gx, curve.gy, k1);
        let o = curve.mul(curve.gx, curve.gy, k2);
        if t.0 == o.0 {
            continue;
        }
        if t.0.is_zero() && t.1.is_zero() {
            continue;
        }
        if o.0.is_zero() && o.1.is_zero() {
            continue;
        }
        let r = curve.add(t.0, t.1, o.0, o.1);
        report.accepted_shots += 1;

        let dx = sub_mod_p(t.0, o.0, curve.modulus);
        let c = sub_mod_p(o.0, r.0, curve.modulus);
        for (label, factor) in [("dx", dx), ("qx_minus_rx", c)] {
            let repairs =
                inversion::shrunken_pz_schedule::thin_factor_repairs_u256(factor);
            if repairs > 0 {
                if trace_thin_sites {
                    let sites = inversion::shrunken_pz_schedule::thin_factor_repair_sites_u256(
                        factor,
                        512,
                    );
                    eprintln!(
                        "TRAILMIX_THIN_REPAIR draw={} label={} repairs={} sites={:?}",
                        draw, label, repairs, sites
                    );
                }
                report.miss_factors += 1;
                report.repair_entries += repairs;
                if report.first_miss.is_none() {
                    report.first_miss = Some((draw, label, repairs));
                }
                if max_misses.is_some_and(|limit| report.miss_factors > limit) {
                    return report;
                }
            }
            if let Some(cap) = qcap_check {
                let entries = inversion::shrunken_pz_schedule::q_cap_overflow_entries_limited_u256(
                    factor,
                    cap,
                    early_limit,
                );
                if entries > 0 {
                    if trace_qcap_steps {
                        let steps = inversion::shrunken_pz_schedule::q_cap_overflow_steps_limited_u256(
                            factor,
                            cap,
                            32,
                        );
                        eprintln!(
                            "TRAILMIX_QCAP_OVERFLOW draw={} label={} cap={} entries={} steps={:?}",
                            draw, label, cap, entries, steps
                        );
                    }
                    report.qcap_miss_factors += 1;
                    report.qcap_entries += entries;
                    if report.first_qcap_miss.is_none() {
                        report.first_qcap_miss = Some((draw, label, entries));
                    }
                    if max_misses.is_some_and(|limit| report.qcap_miss_factors > limit) {
                        return report;
                    }
                }
            }
        }
    }
    report
}

fn tail_nonce_x_op(q: u32) -> Op {
    let mut op = Op::empty();
    op.kind = OperationType::X;
    op.q_target = QubitId(q.into());
    op
}

fn hash_tail_nonce(mut hasher: sha3::Shake256, nonce: u64, q0: u32, q1: u32) -> sha3::Shake256 {
    for i in 0..TRAILMIX_TAIL_NONCE_BITS {
        let q = if (nonce >> i) & 1 == 1 { q1 } else { q0 };
        let op = tail_nonce_x_op(q);
        crate::point_add::B::update_fiat_hash_op(&mut hasher, &op);
        crate::point_add::B::update_fiat_hash_op(&mut hasher, &op);
    }
    hasher
}

fn report_current_support(builder: &crate::point_add::B) {
    if std::env::var("TRAILMIX_SUPPORT_CHECK").ok().as_deref() != Some("1") {
        return;
    }
    let Some(hasher) = builder.clone_fiat_hash() else {
        eprintln!(
            "TRAILMIX_SUPPORT no hash stream; set POINT_ADD_HASH_OPS_LEN in count-only mode"
        );
        return;
    };
    let draws = env_usize("TRAILMIX_SUPPORT_SHOTS", TRAILMIX_NUM_TESTS);
    let report = support_report_for_xof(hasher.finalize_xof(), draws);
    eprintln!(
        "TRAILMIX_SUPPORT draws={} accepted={} miss_factors={} repair_entries={} first_miss={:?} qcap_miss_factors={} qcap_entries={} first_qcap_miss={:?}",
        draws,
        report.accepted_shots,
        report.miss_factors,
        report.repair_entries,
        report.first_miss,
        report.qcap_miss_factors,
        report.qcap_entries,
        report.first_qcap_miss
    );
}

fn search_tail_nonce(builder: &crate::point_add::B, q0: u32, q1: u32) {
    let limit = env_usize("TRAILMIX_TAIL_NONCE_SEARCH", 0);
    if limit == 0 {
        return;
    }
    let Some(base_hasher) = builder.clone_fiat_hash() else {
        eprintln!(
            "TRAILMIX_TAIL_SEARCH no hash stream; set POINT_ADD_HASH_OPS_LEN=base_ops+96 in count-only mode"
        );
        return;
    };
    let start = env_u64("TRAILMIX_TAIL_NONCE_START", 0);
    let draws = env_usize("TRAILMIX_TAIL_NONCE_SHOTS", TRAILMIX_NUM_TESTS);
    let trace = std::env::var("TRAILMIX_TAIL_NONCE_TRACE").is_ok();
    let trace_clean = std::env::var("TRAILMIX_TAIL_NONCE_TRACE_CLEAN")
        .ok()
        .as_deref()
        == Some("1");
    let continue_after_clean = std::env::var("TRAILMIX_TAIL_NONCE_CONTINUE")
        .ok()
        .as_deref()
        == Some("1");
    let early_miss = std::env::var("TRAILMIX_TAIL_NONCE_EARLY_MISS")
        .ok()
        .as_deref()
        == Some("1");
    let progress_every = env_usize("TRAILMIX_TAIL_NONCE_PROGRESS", 0);
    let default_threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    let threads = env_usize("TRAILMIX_TAIL_NONCE_THREADS", default_threads)
        .max(1)
        .min(limit.max(1));
    eprintln!(
        "TRAILMIX_TAIL_SEARCH_BEGIN start={} limit={} draws={} threads={} early_miss={} qcap_check={:?}",
        start,
        limit,
        draws,
        threads,
        early_miss,
        std::env::var("TRAILMIX_SUPPORT_Q_CAP_CHECK").ok()
    );

    let results: Vec<(Option<(u64, TrailMixSupportReport)>, Option<u64>)> =
        std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(threads);
            for tid in 0..threads {
                let base_hasher = base_hasher.clone();
                handles.push(scope.spawn(move || {
                    let mut best: Option<(u64, TrailMixSupportReport)> = None;
                    let mut clean: Option<u64> = None;
                    let mut off = tid;
                    while off < limit {
                        let nonce = start.wrapping_add(off as u64);
                        let hasher = hash_tail_nonce(base_hasher.clone(), nonce, q0, q1);
                        let mut xof = hasher.finalize_xof();
                        let report = support_report_for_xof_limited(
                            &mut xof,
                            draws,
                            early_miss.then_some(0),
                        );
                        if trace {
                            eprintln!(
                                "TRAILMIX_TAIL_SEARCH nonce={} miss_factors={} repair_entries={} first_miss={:?} qcap_miss_factors={} qcap_entries={} first_qcap_miss={:?}",
                                nonce, report.miss_factors, report.repair_entries, report.first_miss,
                                report.qcap_miss_factors, report.qcap_entries, report.first_qcap_miss
                            );
                        }
                        let better = best.as_ref().map_or(true, |(_, b)| {
                            (
                                report.miss_factors,
                                report.qcap_miss_factors,
                                report.qcap_entries,
                                report.repair_entries,
                            ) < (
                                b.miss_factors,
                                b.qcap_miss_factors,
                                b.qcap_entries,
                                b.repair_entries,
                            )
                        });
                        if better {
                            best = Some((nonce, report.clone()));
                        }
                        if progress_every != 0 && off % progress_every == 0 {
                            eprintln!(
                                "TRAILMIX_TAIL_SEARCH_PROGRESS searched_to={} nonce={} miss_factors={} qcap_miss_factors={} repair_entries={} qcap_entries={}",
                                off + 1,
                                nonce,
                                report.miss_factors,
                                report.qcap_miss_factors,
                                report.repair_entries,
                                report.qcap_entries
                            );
                        }
                        if report.miss_factors == 0 && report.qcap_miss_factors == 0 {
                            if trace_clean {
                                eprintln!("TRAILMIX_TAIL_SEARCH_CANDIDATE nonce={nonce}");
                            }
                            clean = Some(clean.map_or(nonce, |old| old.min(nonce)));
                            if !continue_after_clean {
                                break;
                            }
                        }
                        off += threads;
                    }
                    (best, clean)
                }));
            }
            handles
                .into_iter()
                .map(|h| h.join().expect("tail nonce search worker panicked"))
                .collect()
        });

    let mut best: Option<(u64, TrailMixSupportReport)> = None;
    let mut clean: Option<u64> = None;
    for (worker_best, worker_clean) in results {
        if let Some(nonce) = worker_clean {
            clean = Some(clean.map_or(nonce, |old| old.min(nonce)));
        }
        if let Some((nonce, report)) = worker_best {
            let better = best.as_ref().map_or(true, |(best_nonce, b)| {
                (
                    report.miss_factors,
                    report.qcap_miss_factors,
                    report.qcap_entries,
                    report.repair_entries,
                    nonce,
                ) < (
                    b.miss_factors,
                    b.qcap_miss_factors,
                    b.qcap_entries,
                    b.repair_entries,
                    *best_nonce,
                )
            });
            if better {
                best = Some((nonce, report));
            }
        }
    }
    if let Some((nonce, report)) = best {
        eprintln!(
            "TRAILMIX_TAIL_SEARCH_BEST nonce={} accepted={} miss_factors={} repair_entries={} first_miss={:?} qcap_miss_factors={} qcap_entries={} first_qcap_miss={:?} searched={} threads={}",
            nonce,
            report.accepted_shots,
            report.miss_factors,
            report.repair_entries,
            report.first_miss,
            report.qcap_miss_factors,
            report.qcap_entries,
            report.first_qcap_miss,
            limit,
            threads
        );
    }
    if let Some(nonce) = clean {
        eprintln!("TRAILMIX_TAIL_SEARCH_CLEAN nonce={nonce}");
    }
}

pub fn build_builder() -> crate::point_add::B {
    configure_sub1000_trailmix_route();

    let mut circ = circuit::Circuit::new();
    circ.set_section("trailmix_shrunken_pz");
    let mut tx = circ.alloc_qreg_bits("tx", 256);
    let mut ty = circ.alloc_qreg_bits("ty", 256);
    let ox: Vec<circuit::Cbit> = (0..256).map(|_| circ.alloc_input_bit()).collect();
    let oy: Vec<circuit::Cbit> = (0..256).map(|_| circ.alloc_input_bit()).collect();

    ec::point_add::ec_add_inplace_shrunken_pz(&mut circ, &mut tx, &mut ty, &ox, &oy);

    let mut out = std::mem::take(&mut tx);
    out.extend(std::mem::take(&mut ty));
    let out = circ.defragment(out);
    let tail_q0 = out[0].id();
    let tail_q1 = out[1].id();
    circ.declare_registers(&out[..256], &out[256..512], &ox, &oy);

    search_tail_nonce(&circ.b, tail_q0, tail_q1);

    if let Some(nonce) = std::env::var("TRAILMIX_TAIL_NONCE")
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
    {
        circ.set_section("trailmix_tail_nonce");
        for i in 0..TRAILMIX_TAIL_NONCE_BITS {
            let q = if (nonce >> i) & 1 == 1 {
                &out[1]
            } else {
                &out[0]
            };
            circ.x(q);
            circ.x(q);
        }
    }

    let _ = circ.destroy_sim(out);
    let mut builder = circ.into_builder();
    report_current_support(&builder);
    if std::env::var("TRACE_PHASE_OPS").is_ok() {
        use std::collections::BTreeMap;

        builder.close_counted_phase();
        let top_n = std::env::var("TRACE_PHASE_OPS_TOP")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(40);
        let mut rows = builder.counted_phase_rows.clone();
        rows.sort_by(|a, b| b.ops.cmp(&a.ops).then_with(|| a.phase.cmp(b.phase)));
        eprintln!("=== TrailMix count-only per-phase ops ===");
        eprintln!(
            "{:<56} {:>12} {:>12} {:>12} {:>12}",
            "phase", "ops", "toffoli", "hmr", "r"
        );
        for row in rows.into_iter().take(top_n) {
            eprintln!(
                "{:<56} {:>12} {:>12} {:>12} {:>12}",
                row.phase, row.ops, row.toffoli_ops, row.hmr_ops, row.r_ops
            );
        }
        let mut by_phase: BTreeMap<&'static str, crate::point_add::PhaseResource> =
            BTreeMap::new();
        for row in &builder.counted_phase_rows {
            let entry = by_phase
                .entry(row.phase)
                .or_insert(crate::point_add::PhaseResource {
                    phase: row.phase,
                    start: 0,
                    end: 0,
                    ops: 0,
                    toffoli_ops: 0,
                    ccx_ops: 0,
                    ccz_ops: 0,
                    hmr_ops: 0,
                    r_ops: 0,
                });
            entry.ops += row.ops;
            entry.toffoli_ops += row.toffoli_ops;
            entry.ccx_ops += row.ccx_ops;
            entry.ccz_ops += row.ccz_ops;
            entry.hmr_ops += row.hmr_ops;
            entry.r_ops += row.r_ops;
        }
        let mut agg: Vec<_> = by_phase.into_values().collect();
        agg.sort_by(|a, b| b.ops.cmp(&a.ops).then_with(|| a.phase.cmp(b.phase)));
        eprintln!("=== TrailMix aggregate per-phase ops ===");
        eprintln!(
            "{:<56} {:>12} {:>12} {:>12} {:>12}",
            "phase", "ops", "toffoli", "hmr", "r"
        );
        for row in agg.into_iter().take(top_n) {
            eprintln!(
                "{:<56} {:>12} {:>12} {:>12} {:>12}",
                row.phase, row.ops, row.toffoli_ops, row.hmr_ops, row.r_ops
            );
        }
    }
    if std::env::var("TRACE_PEAK").is_ok() || std::env::var("TRACE_PHASE_ACTIVE").is_ok() {
        builder.close_phase_active_region();
        eprintln!(
            "TRAILMIX_SHRUNKEN_PZ peak_qubits={} peak_phase='{}' ops={}",
            builder.peak_qubits,
            builder.peak_phase,
            builder.current_ops_len()
        );
        if std::env::var("TRACE_PHASE_ACTIVE").is_ok() {
            let top_n = std::env::var("TRACE_PHASE_ACTIVE_TOP")
                .ok()
                .and_then(|s| s.parse::<usize>().ok());
            let mut rows: Vec<_> = builder.phase_active_max.iter().collect();
            rows.sort_by(|a, b| b.1.cmp(a.1).then_with(|| a.0.cmp(b.0)));
            for (idx, (phase, active)) in rows.into_iter().enumerate() {
                if top_n.is_some_and(|limit| idx >= limit) {
                    break;
                }
                eprintln!("TRAILMIX_ACTIVE {:<48} {}", phase, active);
            }
        }
    }
    builder
}
