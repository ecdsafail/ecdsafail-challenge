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
    pub mod register_shared_eea;
    pub mod register_shared_eea_microkernels;
    pub mod register_shared_eea_reference;
    pub mod q944_dirty_catalytic_predicate;
    pub mod q944_gate_host_feasibility;
    pub mod q944_gate_host_lifecycle;
    pub mod q944_dirty_parity_microkernels;
    pub mod q944_full_structural;
    pub mod q944_quotient_witness;
    pub mod q945_local_hosts;
    pub mod q949_robust_envelope;
    pub mod shrunken_pz_primitives;
    pub mod shrunken_pz_schedule;
    pub mod shrunken_pz_state_machine;
}

pub mod ec {
    pub mod point_add;
}

use alloy_primitives::U256;
use sha3::digest::{ExtendableOutput, Update, XofReader};
use std::collections::BTreeMap;
use std::fmt::Write as _;

use crate::circuit::{Op, OperationType, QubitId};
use crate::weierstrass_elliptic_curve::WeierstrassEllipticCurve;

const TRAILMIX_TAIL_NONCE_BITS: u32 = 48;
const TRAILMIX_NUM_TESTS: usize = 9024;
const Q851_RESEARCH_SOURCE_ID: &str =
    "151230cd03cedcc6095b0eb10dcf74761e6982ec54c132086f3498a537ae8815";
pub const Q949_PROOF_DRAWS: usize = TRAILMIX_NUM_TESTS;
pub const Q949_PROOF_FACTORS: usize = 2 * Q949_PROOF_DRAWS;
pub const Q949_CENSUS_DIAGNOSTICS_SCHEMA: &str = "q949-census-diagnostics-v2";
pub const Q949_CENSUS_DIAGNOSTICS_OUT_ENV: &str = "Q949_CENSUS_DIAGNOSTICS_OUT";
pub const Q949_WIDTH_CONTEXTS_PER_FACTOR: usize =
    3 * inversion::shrunken_pz_schedule::SHRUNKEN_PZ_NSTEPS
        + 2 * (inversion::shrunken_pz_schedule::SHRUNKEN_PZ_NSTEPS - 1);

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
    set_default_env("TRAILMIX_THIN_CLZ_WINDOW", "78");
    set_default_env("TRAILMIX_THIN_MARGIN", "0");
    set_default_env("TRAILMIX_THIN_VALIDATE", "500000");
    set_default_env("TRAILMIX_COUNTER_W", "8");
    // Selective per-step peak target: clamp ONLY the peak-binding step's quotient
    // so the global peak drops 980 -> 979 while non-peak steps keep full q (vs a
    // blunt global Q_CAP=20 that clamps all ~490 steps and manufactures misses).
    // Q_CAP=99 neutralizes the old global clamp; TRAILMIX_Q_TARGET governs.
    // Q684 experiment: preserve the audited quotient widths, fuse sign with
    // parity, and remove the hybrid-CLZ carry allocation. Passenger-top reuse
    // is forbidden because lambda_raw is not guaranteed canonical.
    set_default_env("TRAILMIX_Q_CAP", "99");
    set_default_env("TRAILMIX_Q_TARGET", "683");
    set_default_env("TRAILMIX_SIGN_PARITY_Q_REUSE", "1");
    set_default_env("LOWQ_CLZ_DIFF_CONST_FOLD", "1");
    set_default_env("LOWQ_HYBRID_CLZ", "1");
    set_default_env("LOWQ_HYBRID_CLZ_KG_MCX", "1");
    set_default_env("LOWQ_HYBRID_CLZ_PREFIX_PARITY", "1");
    set_default_env("LOWQ_HYBRID_CLZ_NOALLOC_ADD", "1");
    set_default_env("LOWQ_EXACT_CTZ", "1");
    set_default_env("LOWQ_Q959_SELECTIVE_BORROW", "1");
    set_default_env("LOWQ_Q958_GATED_COMPARE", "1");
    set_default_env("LOWQ_Q957_TARGET683", "1");
    // Experimental Q956 counter-lane borrowing stays opt-in until the
    // support-preservation proof and authoritative profile are accepted.
    set_default_env("LOWQ_Q956_OFF_BORROW", "0");
    // Q954 composes Q955 with a schedule-certified counter[7]/s_rot alias and
    // canonical passenger-top release. It remains an explicit experiment.
    set_default_env("LOWQ_Q954_SROT_COUNTER7", "0");
    set_default_env("LOWQ_Q953_SROT_COUNTER67", "0");
    // Q949 replaces the persistent explicit counter with one done lane and an
    // affine terminal overlay in ca. Keep it opt-in until its WMI proof/profile
    // pair is complete.
    set_default_env("LOWQ_Q949_AFFINE_COUNTER", "0");
    // The two-stream robust row envelope changes every dynamic register width
    // and therefore the Fiat-Shamir operation hash. It remains independently
    // opt-in and requires a fresh support certificate before profiling.
    set_default_env("LOWQ_Q949_ROBUST_SYMMETRIC_SCHEDULE", "0");
    set_default_env("LOWQ_Q949_ROBUST_FRESH_SUPPORT_CERTIFIED", "0");
    set_default_env("LOWQ_PASSENGER_TOP_LIFETIME_EXPERIMENT", "0");
    set_default_env("LOWQ_BORROWED_TRANSCRIPT_EXPERIMENT", "0");
    set_default_env("LOWQ_BORROWED_TRANSCRIPT_FRESH_SUPPORT_CERTIFIED", "0");
    set_default_env("LOWQ_REVERSE_CA255_RELATIONAL_LOAN_EXPERIMENT", "0");
    set_default_env("LOWQ_Q948_DIRECT_HCLZ_PEAK_GUARD", "0");
    set_default_env("LOWQ_Q947_PASSENGER_DIRECT_HCLZ", "0");
    set_default_env("LOWQ_Q947_FRESH_SUPPORT_CERTIFIED", "0");
    set_default_env("LOWQ_Q946_SECOND_OWNERSHIP_RELEASE", "0");
    set_default_env("LOWQ_Q945_LOCAL_HOSTS", "0");
    set_default_env("LOWQ_Q945_DIRTY_PARITY_ARITHMETIC", "0");
    set_default_env("LOWQ_Q944_FULL_STRUCTURAL", "0");
    set_default_env("LOWQ_Q944_RESIDUAL_ONE_LANE_CUT", "0");
    set_default_env("TRAILMIX_SROT_W", "5");
    set_default_env("TRAILMIX_DEFER_Y_MATERIALIZE", "1");
    set_default_env("TRAILMIX_ZERO_DY_NEWDX_ROUTE", "1");
    // Q883 register-sharing candidate. Keep the route source-baked so local
    // trusted replay and the submission server build the same operation stream.
    set_default_env("TRAILMIX_REGISTER_SHARED_EEA", "1");
    // Q855 register-shared route. These exact arithmetic and lifetime
    // reductions were composed in the whole-circuit profile before baking.
    set_default_env("LOWQ_DIRECT_PREFIX_BITLEN", "1");
    set_default_env("LOWQ_REUSE_ZERO_CARRIES_FOR_PREFIX", "1");
    set_default_env("LOWQ_REUSE_ZERO_CARRIES_FOR_FULL_PREFIX_SCRATCH", "1");
    set_default_env("LOWQ_REUSE_ROTATED_BITLEN_SCRATCH", "1");
    set_default_env("LOWQ_REUSE_COEFFICIENT_COMPARATOR_SCRATCH", "1");
    set_default_env("LOWQ_Q839_PHASE_REMAINDER_SCRATCH", "1");
    set_default_env("LOWQ_FUSED_ZERO_PREFIX_BITLEN", "1");
    set_default_env("LOWQ_REUSE_COEFFICIENT_RAW_BITLEN_LOAN", "1");
    set_default_env("LOWQ_INPLACE_ROTATED_BITLEN_BOUNDARY", "1");
    set_default_env("LOWQ_LOAN_FUSED_PREFIX_SCRATCH", "1");
    // The Q847 package composes six independently proved lifetime and stream
    // choices. Explicit proof overrides remain available because defaults do
    // not overwrite caller-provided values.
    set_default_env("LOWQ_REGISTER_SHARED_REVERSE_DECREMENT_STREAM", "1");
    set_default_env("LOWQ_REUSE_LQ_AS_SWAP_OLD_R_LENGTH", "1");
    set_default_env("LOWQ_SPLIT_COEFFICIENT_ROTATION_LIFETIME", "1");
    set_default_env("LOWQ_REUSE_COEFFICIENT_LESS_THAN_LANES", "1");
    set_default_env("LOWQ_CALLER_SCRATCH_KG_REVERSE_DECREMENT", "1");
    set_default_env("LOWQ_REUSE_CLEAN_CHAIN_FOR_COEFFICIENT_ADD", "1");
    // Q845 source-bakes the proved covering cuts as one production route.
    set_default_env("LOWQ_REUSE_PRESERVED_DY_TOP_FOR_PREFIX", "1");
    set_default_env("LOWQ_MIXED_WIDTH_L_R_PRIME", "1");
    // Source-bake the independently proved Q851 stream cuts. Explicit values
    // still override these defaults in reduced proof binaries.
    set_default_env("LOWQ_PAIRED_BITLEN_SOURCE_COMPLEMENT", "1");
    set_default_env("LOWQ_COEFFICIENT_NONNEGATIVE_X_CANCEL", "1");
    set_default_env("LOWQ_Q845_LIFETIME_COEFFICIENT_FUSION", "1");
    set_default_env("LOWQ_FUSE_PROMISED_SWAP_SUPPORT_LIFETIME", "1");
    set_default_env("LOWQ_Q845_SWAP_ONLY_T_PRIME_LENGTH", "1");
    set_default_env("LOWQ_Q851_TRUNCATED_SWAP_ONLY_GUARD", "1");
    set_default_env("LOWQ_Q851_FIXED_SIGN_EVENT", "1");
    set_default_env("TRAILMIX_TAIL_NONCE", "4114534");
}

#[derive(Clone, Debug, Default)]
struct TrailMixSupportReport {
    accepted_shots: usize,
    miss_factors: usize,
    repair_entries: usize,
    first_miss: Option<(usize, &'static str, usize)>,
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
    let trace_coordinates = std::env::var("TRAILMIX_SUPPORT_COORDS")
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
            if trace_coordinates && repairs > 0 {
                let coordinates = inversion::shrunken_pz_schedule::
                    thin_factor_repair_coordinates_u256(factor);
                debug_assert_eq!(repairs, coordinates.len());
                for coordinate in coordinates {
                    eprintln!(
                        "TRAILMIX_SUPPORT_COORD draw={} factor={} step={} register={} observed={} available={} universal={}",
                        draw,
                        label,
                        coordinate.step,
                        coordinate.register,
                        coordinate.observed_width,
                        coordinate.available_width,
                        coordinate.universal_width,
                    );
                }
            }
            if repairs > 0 {
                report.miss_factors += 1;
                report.repair_entries += repairs;
                if report.first_miss.is_none() {
                    report.first_miss = Some((draw, label, repairs));
                }
                if max_misses.is_some_and(|limit| report.miss_factors > limit) {
                    return report;
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
        "TRAILMIX_SUPPORT draws={} accepted={} miss_factors={} repair_entries={} first_miss={:?}",
        draws,
        report.accepted_shots,
        report.miss_factors,
        report.repair_entries,
        report.first_miss
    );
}

fn q949_hex(bytes: &[u8]) -> String {
    let mut encoded = String::with_capacity(2 * bytes.len());
    for byte in bytes {
        write!(&mut encoded, "{byte:02x}").expect("write Q949 identity hex");
    }
    encoded
}

fn q949_finish_identity(hasher: sha3::Shake256) -> String {
    let mut output = [0u8; 32];
    hasher.finalize_xof().read(&mut output);
    q949_hex(&output)
}

fn q949_hash_component(hasher: &mut sha3::Shake256, name: &str, bytes: &[u8]) {
    hasher.update(&(name.len() as u64).to_le_bytes());
    hasher.update(name.as_bytes());
    hasher.update(&(bytes.len() as u64).to_le_bytes());
    hasher.update(bytes);
}

#[must_use]
pub fn q949_source_identity() -> String {
    Q851_RESEARCH_SOURCE_ID.to_owned()
}

#[must_use]
pub fn q949_schedule_identity() -> String {
    use inversion::shrunken_pz_schedule::{
        q949_effective_reg_los, q949_effective_reg_widths, shift_bounds,
        SHRUNKEN_PZ_NSTEPS,
    };

    let mut hasher = sha3::Shake256::default();
    hasher.update(b"q949-effective-schedule-identity-v2");
    hasher.update(&(SHRUNKEN_PZ_NSTEPS as u64).to_le_bytes());
    for row in 0..SHRUNKEN_PZ_NSTEPS {
        hasher.update(&(row as u64).to_le_bytes());
        for value in q949_effective_reg_widths(row) {
            hasher.update(&(value as u64).to_le_bytes());
        }
        for value in q949_effective_reg_los(row) {
            hasher.update(&(value as u64).to_le_bytes());
        }
        let (division, multiply) = shift_bounds(row);
        hasher.update(&(division as u64).to_le_bytes());
        hasher.update(&(multiply as u64).to_le_bytes());
    }
    q949_finish_identity(hasher)
}

#[must_use]
pub fn q949_route_identity() -> String {
    const ROUTE_ENV: &[&str] = &[
        "TRAILMIX_THIN_SCHEDULE",
        "TRAILMIX_THIN_TRAIN",
        "TRAILMIX_THIN_SEED",
        "TRAILMIX_THIN_CLZ_WINDOW",
        "TRAILMIX_THIN_MARGIN",
        "TRAILMIX_THIN_VALIDATE",
        "TRAILMIX_THIN_REPAIR_MARGIN",
        "TRAILMIX_Q_TARGET",
        "TRAILMIX_Q_CAP",
        "TRAILMIX_COUNTER_W",
        "TRAILMIX_SROT_W",
        "TRAILMIX_SIGN_PARITY_Q_REUSE",
        "TRAILMIX_PASSENGER_TOP_Q_REUSE",
        "TRAILMIX_Q_MODEL_GUARD",
        "TRAILMIX_AB_CAP",
        "TRAILMIX_CACB_CAP",
        "LOWQ_CLZ_DIFF_CONST_FOLD",
        "LOWQ_HYBRID_CLZ",
        "LOWQ_HYBRID_CLZ_KG_MCX",
        "LOWQ_HYBRID_CLZ_PREFIX_PARITY",
        "LOWQ_HYBRID_CLZ_NOALLOC_ADD",
        "LOWQ_EXACT_CTZ",
        "LOWQ_Q959_SELECTIVE_BORROW",
        "LOWQ_Q958_GATED_COMPARE",
        "LOWQ_Q957_TARGET683",
        "LOWQ_Q956_OFF_BORROW",
        "LOWQ_Q955_OFF_CANONICAL",
        "LOWQ_Q954_SROT_COUNTER7",
        "LOWQ_Q953_SROT_COUNTER67",
        "LOWQ_Q949_AFFINE_COUNTER",
        "LOWQ_Q949_ROBUST_SYMMETRIC_SCHEDULE",
        "LOWQ_PASSENGER_TOP_LIFETIME_EXPERIMENT",
        "LOWQ_BORROWED_TRANSCRIPT_EXPERIMENT",
        "LOWQ_REVERSE_CA255_RELATIONAL_LOAN_EXPERIMENT",
        "LOWQ_Q948_DIRECT_HCLZ_PEAK_GUARD",
        "LOWQ_Q947_PASSENGER_DIRECT_HCLZ",
        "LOWQ_Q946_SECOND_OWNERSHIP_RELEASE",
        "LOWQ_Q945_LOCAL_HOSTS",
        "LOWQ_Q945_DIRTY_PARITY_ARITHMETIC",
        "LOWQ_Q944_FULL_STRUCTURAL",
        "LOWQ_Q944_RESIDUAL_ONE_LANE_CUT",
        "TRAILMIX_DEFER_Y_MATERIALIZE",
        "TRAILMIX_ZERO_DY_NEWDX_ROUTE",
        "TRAILMIX_REGISTER_SHARED_EEA",
        "LOWQ_REGISTER_SHARED_REVERSE_DECREMENT_STREAM",
        "LOWQ_REUSE_LQ_AS_SWAP_OLD_R_LENGTH",
        "LOWQ_SPLIT_COEFFICIENT_ROTATION_LIFETIME",
        "LOWQ_REUSE_COEFFICIENT_LESS_THAN_LANES",
        "LOWQ_CALLER_SCRATCH_KG_REVERSE_DECREMENT",
        "LOWQ_REUSE_CLEAN_CHAIN_FOR_COEFFICIENT_ADD",
        "LOWQ_REUSE_PRESERVED_DY_TOP_FOR_PREFIX",
        "LOWQ_MIXED_WIDTH_L_R_PRIME",
        "LOWQ_SUB800_INPLACE_GUARD_ADDRESS",
        "TRAILMIX_TAIL_NONCE",
    ];

    let mut hasher = sha3::Shake256::default();
    hasher.update(b"q949-route-identity-v4");
    for name in ROUTE_ENV {
        let value = std::env::var(name).unwrap_or_else(|_| "<unset>".to_owned());
        q949_hash_component(&mut hasher, name, value.as_bytes());
    }
    q949_finish_identity(hasher)
}

fn q949_builder_fiat_hash(builder: &crate::point_add::B) -> sha3::Shake256 {
    builder.clone_fiat_hash().unwrap_or_else(|| {
        assert!(
            !builder.count_only,
            "Q949 proof cannot reconstruct a count-only hash"
        );
        let mut hasher = sha3::Shake256::default();
        hasher.update(b"quantum_ecc-fiat-shamir-v2");
        hasher.update(&(builder.ops.len() as u64).to_le_bytes());
        for op in &builder.ops {
            crate::point_add::B::update_fiat_hash_op(&mut hasher, op);
        }
        hasher
    })
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Q949ProofIdentity {
    pub source_id: String,
    pub schedule_id: String,
    pub route_id: String,
    pub tail_nonce: u64,
    pub tail_nonce_bits: u32,
    pub op_count: usize,
    pub op_stream_id: String,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Q925RegisterSharedProofIdentity {
    pub implementation_commit: String,
    pub point_add_tree: String,
    pub source_id: String,
    pub schedule_id: String,
    pub route_id: String,
    pub tail_nonce: u64,
    pub tail_nonce_bits: u32,
    pub op_count: usize,
    pub op_stream_id: String,
}

#[must_use]
pub fn q925_register_shared_schedule_identity() -> String {
    use inversion::register_shared_eea_reference::{
        reference_active_windows, REFERENCE_LENGTH_WIDTH, REFERENCE_R_LENGTH_WIDTH,
        REFERENCE_STEPS,
    };

    let mut hasher = sha3::Shake256::default();
    hasher.update(b"q925-register-shared-schedule-identity-v2");
    hasher.update(&(REFERENCE_STEPS as u64).to_le_bytes());
    hasher.update(&(REFERENCE_LENGTH_WIDTH as u64).to_le_bytes());
    hasher.update(&(REFERENCE_R_LENGTH_WIDTH as u64).to_le_bytes());
    for step in 1..=REFERENCE_STEPS {
        let windows = reference_active_windows(256, step);
        hasher.update(&(step as u64).to_le_bytes());
        for (start, end) in [
            windows.r_add_sub,
            windows.quotient_swap,
            windows.t_add_sub,
            windows.length_update_t,
            windows.length_update_r,
        ] {
            hasher.update(&(start as u64).to_le_bytes());
            hasher.update(&(end as u64).to_le_bytes());
        }
    }
    q949_finish_identity(hasher)
}

#[must_use]
pub fn q925_register_shared_proof_identity(
    builder: &crate::point_add::B,
    implementation_commit: &str,
    point_add_tree: &str,
) -> Q925RegisterSharedProofIdentity {
    assert_eq!(
        implementation_commit.len(),
        40,
        "Q925 implementation commit must be a full Git SHA-1"
    );
    assert!(
        implementation_commit
            .bytes()
            .all(|value| value.is_ascii_hexdigit()),
        "Q925 implementation commit must be hexadecimal"
    );
    assert_eq!(
        point_add_tree.len(),
        40,
        "Q925 point-add tree must be a full Git object ID"
    );
    assert!(
        point_add_tree
            .bytes()
            .all(|value| value.is_ascii_hexdigit()),
        "Q925 point-add tree must be hexadecimal"
    );
    let tail_nonce = std::env::var("TRAILMIX_TAIL_NONCE")
        .expect("Q925 identity requires TRAILMIX_TAIL_NONCE")
        .parse::<u64>()
        .expect("Q925 tail nonce must be an integer");
    assert!(
        tail_nonce < (1u64 << TRAILMIX_TAIL_NONCE_BITS),
        "Q925 tail nonce exceeds its encoded width"
    );
    if !builder.count_only {
        assert_eq!(builder.counted_ops, builder.ops.len());
    }

    let mut source_hasher = sha3::Shake256::default();
    source_hasher.update(b"q925-point-add-tree-source-binding-v2");
    source_hasher.update(point_add_tree.as_bytes());

    Q925RegisterSharedProofIdentity {
        implementation_commit: implementation_commit.to_owned(),
        point_add_tree: point_add_tree.to_owned(),
        source_id: q949_finish_identity(source_hasher),
        schedule_id: q925_register_shared_schedule_identity(),
        route_id: q949_route_identity(),
        tail_nonce,
        tail_nonce_bits: TRAILMIX_TAIL_NONCE_BITS,
        op_count: builder.counted_ops,
        op_stream_id: q949_finish_identity(q949_builder_fiat_hash(builder)),
    }
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum Q945SupportPhase {
    InvFwd,
    AltCancel,
}

impl Q945SupportPhase {
    pub const fn label(self) -> &'static str {
        match self {
            Self::InvFwd => "ec3.inv_fwd",
            Self::AltCancel => "ec3.alt.cancel",
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct Q945HclzSupportSite {
    pub phase: Q945SupportPhase,
    pub direction: inversion::shrunken_pz_schedule::Q949TraceDirection,
    pub row: usize,
    pub substep: inversion::q945_local_hosts::Q945Substep,
    pub form: inversion::q945_local_hosts::Q945HclzForm,
    pub host: inversion::q945_local_hosts::Q945Host,
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct Q945CarrySupportSite {
    pub phase: Q945SupportPhase,
    pub direction: inversion::shrunken_pz_schedule::Q949TraceDirection,
    pub row: usize,
    pub substep: inversion::q945_local_hosts::Q945Substep,
    pub host: inversion::q945_local_hosts::Q945Host,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Q945SupportSiteCount<S> {
    pub site: S,
    pub checks: usize,
    pub zero_entry_checks: usize,
    pub restoration_checks: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Q945SupportedHostMiss {
    pub draw: usize,
    pub factor_label: &'static str,
    pub factor: U256,
    pub phase: Q945SupportPhase,
    pub direction: inversion::shrunken_pz_schedule::Q949TraceDirection,
    pub row: usize,
    pub substep: inversion::q945_local_hosts::Q945Substep,
    pub form: Option<inversion::q945_local_hosts::Q945HclzForm>,
    pub host: inversion::q945_local_hosts::Q945Host,
    pub entry_value: bool,
    pub exit_value: bool,
    pub reason: &'static str,
    pub expected_lt: Option<bool>,
    pub full_lt: Option<bool>,
    pub route_lt: Option<bool>,
    pub boundary: inversion::shrunken_pz_schedule::Q945HostBoundaryState,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Q945SupportedNarrowCompareMiss {
    pub draw: usize,
    pub factor_label: &'static str,
    pub factor: U256,
    pub phase: Q945SupportPhase,
    pub coordinate: inversion::shrunken_pz_schedule::Q949NarrowCompareMiss,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Q925CommittedFactor {
    pub draw: usize,
    pub factor_label: &'static str,
    pub factor: U256,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Q925CommittedFactorCorpus {
    pub requested_draws: usize,
    pub accepted_draws: usize,
    pub rejected_draws: usize,
    pub factors: Vec<Q925CommittedFactor>,
}

/// Materialize the exact factor corpus committed by the builder's operation
/// stream. This analysis API preserves rejected-draw accounting and performs no
/// nonce search, so downstream width claims remain tied to one fixed circuit.
#[doc(hidden)]
pub fn q925_committed_factor_corpus(
    builder: &crate::point_add::B,
) -> Q925CommittedFactorCorpus {
    let requested_draws = Q949_PROOF_DRAWS;
    let mut xof = q949_builder_fiat_hash(builder).finalize_xof();
    let curve = secp256k1();
    let mut accepted_draws = 0usize;
    let mut rejected_draws = 0usize;
    let mut factors = Vec::with_capacity(Q949_PROOF_FACTORS);
    for draw in 0..requested_draws {
        let mut random = [[0u8; 32]; 2];
        xof.read(&mut random[0]);
        xof.read(&mut random[1]);
        let k1 = U256::from_le_bytes(random[0]);
        let k2 = U256::from_le_bytes(random[1]);
        let target = curve.mul(curve.gx, curve.gy, k1);
        let other = curve.mul(curve.gx, curve.gy, k2);
        if target.0 == other.0
            || (target.0.is_zero() && target.1.is_zero())
            || (other.0.is_zero() && other.1.is_zero())
        {
            rejected_draws += 1;
            continue;
        }
        accepted_draws += 1;
        let result = curve.add(target.0, target.1, other.0, other.1);
        factors.push(Q925CommittedFactor {
            draw,
            factor_label: "dx",
            factor: sub_mod_p(target.0, other.0, curve.modulus),
        });
        factors.push(Q925CommittedFactor {
            draw,
            factor_label: "qx_minus_rx",
            factor: sub_mod_p(other.0, result.0, curve.modulus),
        });
    }
    assert_eq!(accepted_draws + rejected_draws, requested_draws);
    assert_eq!(factors.len(), 2 * accepted_draws);
    Q925CommittedFactorCorpus {
        requested_draws,
        accepted_draws,
        rejected_draws,
        factors,
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Q945HostSupportReport {
    pub requested_draws: usize,
    pub accepted_draws: usize,
    pub rejected_draws: usize,
    pub factors_checked: usize,
    pub hclz_host_checks: usize,
    pub hclz_zero_entry_checks: usize,
    pub hclz_restoration_checks: usize,
    pub hclz_host_misses: usize,
    pub carry_host_checks: usize,
    pub carry_zero_entry_checks: usize,
    pub carry_restoration_checks: usize,
    pub carry_host_misses: usize,
    pub carry_semantic_checks: usize,
    pub carry_semantic_misses: usize,
    pub row364_checks: usize,
    pub row364_b80_zero_checks: usize,
    pub row364_identity_checks: usize,
    pub row364_misses: usize,
    pub row374_q24_checks: usize,
    pub row374_q24_zero_checks: usize,
    pub row374_q24_noncarry_touches: usize,
    pub row374_misses: usize,
    pub preterminal_counter_off_checks: usize,
    pub preterminal_counter_off_zero_checks: usize,
    pub preterminal_counter_off_misses: usize,
    pub row385_special_checks: usize,
    pub row385_special_zero_checks: usize,
    pub row385_special_misses: usize,
    pub width_misses: usize,
    pub clz_window_misses: usize,
    pub narrow_compare_checks: usize,
    pub narrow_compare_misses: usize,
    pub division_offset_compare_checks: usize,
    pub division_offset_compare_misses: usize,
    pub multiply_offset_cleanup_compare_checks: usize,
    pub multiply_offset_cleanup_compare_misses: usize,
    pub first_host_miss: Option<Q945SupportedHostMiss>,
    pub first_narrow_compare_miss: Option<Q945SupportedNarrowCompareMiss>,
    pub hclz_sites: Vec<Q945SupportSiteCount<Q945HclzSupportSite>>,
    pub carry_sites: Vec<Q945SupportSiteCount<Q945CarrySupportSite>>,
    pub earliest_terminal_row: usize,
    pub latest_terminal_row: usize,
    pub support_clean: bool,
}

#[must_use]
pub fn q949_proof_identity(builder: &crate::point_add::B) -> Q949ProofIdentity {
    let tail_nonce = std::env::var("TRAILMIX_TAIL_NONCE")
        .expect("Q949 identity requires TRAILMIX_TAIL_NONCE")
        .parse::<u64>()
        .expect("Q949 tail nonce must be an integer");
    assert!(
        tail_nonce < (1u64 << TRAILMIX_TAIL_NONCE_BITS),
        "Q949 tail nonce exceeds its encoded width"
    );
    if !builder.count_only {
        assert_eq!(builder.counted_ops, builder.ops.len());
    }
    let op_stream_id = q949_finish_identity(q949_builder_fiat_hash(builder));
    Q949ProofIdentity {
        source_id: q949_source_identity(),
        schedule_id: q949_schedule_identity(),
        route_id: q949_route_identity(),
        tail_nonce,
        tail_nonce_bits: TRAILMIX_TAIL_NONCE_BITS,
        op_count: builder.counted_ops,
        op_stream_id,
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Q949SupportedWidthMiss {
    pub draw: usize,
    pub factor_label: &'static str,
    pub factor: U256,
    pub coordinate: inversion::shrunken_pz_schedule::Q949WidthMiss,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Q949SupportedClzWindowMiss {
    pub draw: usize,
    pub factor_label: &'static str,
    pub factor: U256,
    pub coordinate: inversion::shrunken_pz_schedule::Q949ClzWindowMiss,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Q949WidthMissBucket {
    pub direction: inversion::shrunken_pz_schedule::Q949TraceDirection,
    pub phase: inversion::shrunken_pz_schedule::Q949WidthPhase,
    pub row: usize,
    pub register: &'static str,
    pub miss_count: usize,
    pub min_observed_width: usize,
    pub max_observed_width: usize,
    pub available_width: usize,
    pub max_excess: usize,
}

impl Q949WidthMissBucket {
    fn new(coordinate: inversion::shrunken_pz_schedule::Q949WidthMiss) -> Self {
        assert!(coordinate.observed_width > coordinate.available_width);
        Self {
            direction: coordinate.direction,
            phase: coordinate.phase,
            row: coordinate.row,
            register: coordinate.register,
            miss_count: 1,
            min_observed_width: coordinate.observed_width,
            max_observed_width: coordinate.observed_width,
            available_width: coordinate.available_width,
            max_excess: coordinate.observed_width - coordinate.available_width,
        }
    }

    fn record(&mut self, coordinate: inversion::shrunken_pz_schedule::Q949WidthMiss) {
        assert!(coordinate.observed_width > coordinate.available_width);
        assert_eq!(self.direction, coordinate.direction);
        assert_eq!(self.phase, coordinate.phase);
        assert_eq!(self.row, coordinate.row);
        assert_eq!(self.register, coordinate.register);
        assert_eq!(self.available_width, coordinate.available_width);
        self.miss_count += 1;
        self.min_observed_width = self.min_observed_width.min(coordinate.observed_width);
        self.max_observed_width = self.max_observed_width.max(coordinate.observed_width);
        self.max_excess = self
            .max_excess
            .max(coordinate.observed_width - coordinate.available_width);
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Q949ClzWindowMissBucket {
    pub direction: inversion::shrunken_pz_schedule::Q949TraceDirection,
    pub row: usize,
    pub register: &'static str,
    pub miss_count: usize,
    pub min_observed_width: usize,
    pub max_observed_width: usize,
    pub low: usize,
    pub available_width: usize,
    pub max_shortfall: usize,
}

impl Q949ClzWindowMissBucket {
    fn new(coordinate: inversion::shrunken_pz_schedule::Q949ClzWindowMiss) -> Self {
        assert!(coordinate.observed_width > 0);
        assert!(coordinate.observed_width <= coordinate.low);
        Self {
            direction: coordinate.direction,
            row: coordinate.row,
            register: coordinate.register,
            miss_count: 1,
            min_observed_width: coordinate.observed_width,
            max_observed_width: coordinate.observed_width,
            low: coordinate.low,
            available_width: coordinate.available_width,
            max_shortfall: coordinate.low + 1 - coordinate.observed_width,
        }
    }

    fn record(&mut self, coordinate: inversion::shrunken_pz_schedule::Q949ClzWindowMiss) {
        assert!(coordinate.observed_width > 0);
        assert!(coordinate.observed_width <= coordinate.low);
        assert_eq!(self.direction, coordinate.direction);
        assert_eq!(self.row, coordinate.row);
        assert_eq!(self.register, coordinate.register);
        assert_eq!(self.low, coordinate.low);
        assert_eq!(self.available_width, coordinate.available_width);
        self.miss_count += 1;
        self.min_observed_width = self.min_observed_width.min(coordinate.observed_width);
        self.max_observed_width = self.max_observed_width.max(coordinate.observed_width);
        self.max_shortfall = self
            .max_shortfall
            .max(coordinate.low + 1 - coordinate.observed_width);
    }
}

const Q949_WIDTH_REGISTERS: [&str; 5] = ["A", "B", "ca", "cb", "q"];
const Q949_CLZ_REGISTERS: [&str; 4] = ["A", "B", "ca", "cb"];

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Q949WidthDemandWitness {
    pub draw: usize,
    pub factor_label: &'static str,
    pub factor: U256,
    pub required_widths: [usize; 5],
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Q949JointWidthDemand {
    pub direction: inversion::shrunken_pz_schedule::Q949TraceDirection,
    pub phase: inversion::shrunken_pz_schedule::Q949WidthPhase,
    pub row: usize,
    pub observation_count: usize,
    pub required_widths: [usize; 5],
    pub available_widths: [usize; 5],
    pub minimum_fixed_capacity_sum: usize,
    pub slack_vs_target_683: i64,
    pub component_witnesses: [Q949WidthDemandWitness; 5],
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Q949ClzLowWitness {
    pub draw: usize,
    pub factor_label: &'static str,
    pub factor: U256,
    pub observed_width: usize,
    pub observed_widths: [usize; 4],
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Q949ClzLowBound {
    pub direction: inversion::shrunken_pz_schedule::Q949TraceDirection,
    pub row: usize,
    pub register: &'static str,
    pub observation_count: usize,
    pub nonzero_observation_count: usize,
    pub zero_observation_count: usize,
    pub minimum_nonzero_observed_width: Option<usize>,
    pub safe_low_upper_bound: Option<usize>,
    pub current_low: usize,
    pub available_width: usize,
    pub minimum_witness: Option<Q949ClzLowWitness>,
}

#[derive(Clone, Debug)]
struct Q949JointWidthDemandAccumulator {
    direction: inversion::shrunken_pz_schedule::Q949TraceDirection,
    phase: inversion::shrunken_pz_schedule::Q949WidthPhase,
    row: usize,
    observation_count: usize,
    required_widths: [usize; 5],
    available_widths: [usize; 5],
    component_witnesses: [Option<Q949WidthDemandWitness>; 5],
}

impl Q949JointWidthDemandAccumulator {
    fn new(
        observation: inversion::shrunken_pz_schedule::Q949WidthObservation,
        witness: Q949WidthDemandWitness,
    ) -> Self {
        let mut accumulator = Self {
            direction: observation.direction,
            phase: observation.phase,
            row: observation.row,
            observation_count: 0,
            required_widths: [0; 5],
            available_widths: observation.available_widths,
            component_witnesses: [None; 5],
        };
        accumulator.record(observation, witness);
        accumulator
    }

    fn record(
        &mut self,
        observation: inversion::shrunken_pz_schedule::Q949WidthObservation,
        witness: Q949WidthDemandWitness,
    ) {
        assert_eq!(self.direction, observation.direction);
        assert_eq!(self.phase, observation.phase);
        assert_eq!(self.row, observation.row);
        assert_eq!(self.available_widths, observation.available_widths);
        assert_eq!(witness.required_widths, observation.required_widths);
        self.observation_count += 1;
        for register in 0..5 {
            let required = observation.required_widths[register];
            let replace = required > self.required_widths[register]
                || (required == self.required_widths[register]
                    && self.component_witnesses[register].is_some_and(|current| {
                        q949_witness_rank(witness.draw, witness.factor_label)
                            < q949_witness_rank(current.draw, current.factor_label)
                    }));
            if replace || self.component_witnesses[register].is_none() {
                self.required_widths[register] = required;
                self.component_witnesses[register] = Some(witness);
            }
        }
    }

    fn finish(self) -> Q949JointWidthDemand {
        let minimum_fixed_capacity_sum = self.required_widths.iter().sum::<usize>();
        Q949JointWidthDemand {
            direction: self.direction,
            phase: self.phase,
            row: self.row,
            observation_count: self.observation_count,
            required_widths: self.required_widths,
            available_widths: self.available_widths,
            minimum_fixed_capacity_sum,
            slack_vs_target_683: inversion::shrunken_pz_schedule::Q949_TARGET_SUM as i64
                - minimum_fixed_capacity_sum as i64,
            component_witnesses: self.component_witnesses.map(|witness| {
                witness.expect("Q949 joint-width maximum is missing its witness")
            }),
        }
    }
}

#[derive(Clone, Debug)]
struct Q949ClzLowBoundAccumulator {
    direction: inversion::shrunken_pz_schedule::Q949TraceDirection,
    row: usize,
    register: &'static str,
    observation_count: usize,
    nonzero_observation_count: usize,
    zero_observation_count: usize,
    minimum_nonzero_observed_width: Option<usize>,
    current_low: usize,
    available_width: usize,
    minimum_witness: Option<Q949ClzLowWitness>,
}

impl Q949ClzLowBoundAccumulator {
    fn new(
        observation: inversion::shrunken_pz_schedule::Q949ClzWindowObservation,
        register: usize,
        witness: Q949ClzLowWitness,
    ) -> Self {
        let mut accumulator = Self {
            direction: observation.direction,
            row: observation.row,
            register: Q949_CLZ_REGISTERS[register],
            observation_count: 0,
            nonzero_observation_count: 0,
            zero_observation_count: 0,
            minimum_nonzero_observed_width: None,
            current_low: observation.lows[register],
            available_width: observation.available_widths[register],
            minimum_witness: None,
        };
        accumulator.record(observation, register, witness);
        accumulator
    }

    fn record(
        &mut self,
        observation: inversion::shrunken_pz_schedule::Q949ClzWindowObservation,
        register: usize,
        witness: Q949ClzLowWitness,
    ) {
        assert_eq!(self.direction, observation.direction);
        assert_eq!(self.row, observation.row);
        assert_eq!(self.register, Q949_CLZ_REGISTERS[register]);
        assert_eq!(self.current_low, observation.lows[register]);
        assert_eq!(self.available_width, observation.available_widths[register]);
        assert_eq!(witness.observed_width, observation.observed_widths[register]);
        assert_eq!(witness.observed_widths, observation.observed_widths);
        self.observation_count += 1;
        let observed_width = observation.observed_widths[register];
        if observed_width == 0 {
            self.zero_observation_count += 1;
            return;
        }
        self.nonzero_observation_count += 1;
        let replace = self
            .minimum_nonzero_observed_width
            .is_none_or(|minimum| observed_width < minimum)
            || (self.minimum_nonzero_observed_width == Some(observed_width)
                && self.minimum_witness.is_some_and(|current| {
                    q949_witness_rank(witness.draw, witness.factor_label)
                        < q949_witness_rank(current.draw, current.factor_label)
                }));
        if replace {
            self.minimum_nonzero_observed_width = Some(observed_width);
            self.minimum_witness = Some(witness);
        }
    }

    fn finish(self) -> Q949ClzLowBound {
        Q949ClzLowBound {
            direction: self.direction,
            row: self.row,
            register: self.register,
            observation_count: self.observation_count,
            nonzero_observation_count: self.nonzero_observation_count,
            zero_observation_count: self.zero_observation_count,
            minimum_nonzero_observed_width: self.minimum_nonzero_observed_width,
            safe_low_upper_bound: self
                .minimum_nonzero_observed_width
                .map(|width| width - 1),
            current_low: self.current_low,
            available_width: self.available_width,
            minimum_witness: self.minimum_witness,
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Q949SupportedTraceReport {
    pub requested_draws: usize,
    pub accepted_shots: usize,
    pub rejected_draws: usize,
    pub factors_checked: usize,
    pub forward_rows_checked: usize,
    pub backward_rows_checked: usize,
    pub row_bounds_checked: usize,
    pub entry_width_checks: usize,
    pub transient_width_checks: usize,
    pub post_swap_width_checks: usize,
    pub boundary_width_checks: usize,
    pub width_context_observations: usize,
    pub entry_width_misses: usize,
    pub transient_width_misses: usize,
    pub post_swap_width_misses: usize,
    pub boundary_width_misses: usize,
    pub width_misses: usize,
    pub clz_window_checks: usize,
    pub clz_window_misses: usize,
    pub clz_context_observations: usize,
    pub first_width_miss: Option<Q949SupportedWidthMiss>,
    pub first_clz_window_miss: Option<Q949SupportedClzWindowMiss>,
    pub width_miss_buckets: Vec<Q949WidthMissBucket>,
    pub clz_window_miss_buckets: Vec<Q949ClzWindowMissBucket>,
    pub width_excess_histogram: Vec<(usize, usize)>,
    pub clz_shortfall_histogram: Vec<(usize, usize)>,
    pub joint_width_demands: Vec<Q949JointWidthDemand>,
    pub clz_low_bounds: Vec<Q949ClzLowBound>,
    pub terminal_full_ca_checks: usize,
    pub reverse_row_380_relation_checks: usize,
    pub reverse_row_380_active_checks: usize,
    pub reverse_row_380_inactive_checks: usize,
    pub reverse_row_380_relation_failures: usize,
    pub earliest_terminal_row: usize,
    pub latest_terminal_row: usize,
    pub max_counter: usize,
}

fn q949_direction_label(
    direction: inversion::shrunken_pz_schedule::Q949TraceDirection,
) -> &'static str {
    use inversion::shrunken_pz_schedule::Q949TraceDirection;

    match direction {
        Q949TraceDirection::Forward => "forward",
        Q949TraceDirection::Reverse => "reverse",
    }
}

fn q949_width_phase_label(
    phase: inversion::shrunken_pz_schedule::Q949WidthPhase,
) -> &'static str {
    use inversion::shrunken_pz_schedule::Q949WidthPhase;

    match phase {
        Q949WidthPhase::Entry => "entry",
        Q949WidthPhase::Transient => "transient",
        Q949WidthPhase::PostSwap => "post_swap",
        Q949WidthPhase::Boundary => "boundary",
    }
}

fn q949_factor_label_index(factor_label: &str) -> usize {
    match factor_label {
        "dx" => 0,
        "qx_minus_rx" => 1,
        _ => panic!("unknown Q949 factor label: {factor_label}"),
    }
}

fn q949_witness_rank(draw: usize, factor_label: &str) -> (usize, usize) {
    (draw, q949_factor_label_index(factor_label))
}

fn q949_width_register_index(register: &str) -> usize {
    match register {
        "A" => 0,
        "B" => 1,
        "ca" => 2,
        "cb" => 3,
        "q" => 4,
        _ => panic!("unknown Q949 width register: {register}"),
    }
}

fn q949_joint_width_fits_current(demand: &Q949JointWidthDemand) -> bool {
    (0..5).all(|register| demand.required_widths[register] <= demand.available_widths[register])
}

fn q949_clz_low_is_safe(bound: &Q949ClzLowBound) -> bool {
    bound
        .safe_low_upper_bound
        .is_none_or(|safe_low| bound.current_low <= safe_low)
}

fn q949_clz_register_index(register: &str) -> usize {
    let index = q949_width_register_index(register);
    assert!(index < 4, "Q949 CLZ census cannot contain register q");
    index
}

fn q949_width_phase_miss_totals(report: &Q949SupportedTraceReport) -> [usize; 4] {
    use inversion::shrunken_pz_schedule::Q949WidthPhase;

    let mut totals = [0usize; 4];
    for bucket in &report.width_miss_buckets {
        let index = match bucket.phase {
            Q949WidthPhase::Entry => 0,
            Q949WidthPhase::Transient => 1,
            Q949WidthPhase::PostSwap => 2,
            Q949WidthPhase::Boundary => 3,
        };
        totals[index] += bucket.miss_count;
    }
    totals
}

fn q949_expected_width_context_keys() -> Vec<(
    inversion::shrunken_pz_schedule::Q949TraceDirection,
    inversion::shrunken_pz_schedule::Q949WidthPhase,
    usize,
)> {
    use inversion::shrunken_pz_schedule::{
        Q949TraceDirection, Q949WidthPhase, SHRUNKEN_PZ_NSTEPS,
    };

    let mut keys = Vec::with_capacity(Q949_WIDTH_CONTEXTS_PER_FACTOR);
    for phase in [
        Q949WidthPhase::Entry,
        Q949WidthPhase::Transient,
        Q949WidthPhase::PostSwap,
    ] {
        for row in 0..SHRUNKEN_PZ_NSTEPS {
            keys.push((Q949TraceDirection::Forward, phase, row));
        }
    }
    for row in 1..SHRUNKEN_PZ_NSTEPS {
        keys.push((Q949TraceDirection::Forward, Q949WidthPhase::Boundary, row));
    }
    for row in 0..SHRUNKEN_PZ_NSTEPS - 1 {
        keys.push((Q949TraceDirection::Reverse, Q949WidthPhase::Boundary, row));
    }
    assert_eq!(keys.len(), Q949_WIDTH_CONTEXTS_PER_FACTOR);
    keys
}

fn q949_expected_clz_low_keys(
    latest_terminal_row: usize,
) -> Vec<(
    inversion::shrunken_pz_schedule::Q949TraceDirection,
    usize,
    usize,
)> {
    use inversion::shrunken_pz_schedule::Q949TraceDirection;

    let mut keys = Vec::with_capacity(2 * (latest_terminal_row + 1) * 4);
    for direction in [Q949TraceDirection::Forward, Q949TraceDirection::Reverse] {
        for row in 0..=latest_terminal_row {
            for register in 0..4 {
                keys.push((direction, row, register));
            }
        }
    }
    keys
}

fn q949_validate_census_report(report: &Q949SupportedTraceReport) {
    use inversion::shrunken_pz_schedule::{
        q949_effective_reg_los, q949_effective_reg_widths, Q949TraceDirection,
        Q949WidthPhase, SHRUNKEN_PZ_NSTEPS,
    };

    assert_eq!(
        Q949_CENSUS_DIAGNOSTICS_SCHEMA,
        "q949-census-diagnostics-v2",
        "Q949 diagnostic schema changed without updating its validator"
    );
    assert_eq!(
        report.accepted_shots + report.rejected_draws,
        report.requested_draws,
        "Q949 draw accounting drift"
    );
    assert_eq!(
        report.factors_checked,
        2 * report.accepted_shots,
        "Q949 factor accounting drift"
    );
    let row_checks = report.factors_checked * SHRUNKEN_PZ_NSTEPS;
    assert_eq!(report.forward_rows_checked, row_checks);
    assert_eq!(report.backward_rows_checked, row_checks);
    assert_eq!(report.row_bounds_checked, row_checks);
    assert_eq!(report.entry_width_checks, row_checks * 5);
    assert_eq!(report.transient_width_checks, row_checks * 5);
    assert_eq!(report.post_swap_width_checks, row_checks * 5);
    assert_eq!(
        report.boundary_width_checks,
        report.factors_checked * 2 * (SHRUNKEN_PZ_NSTEPS - 1) * 5
    );
    let total_width_checks = report.entry_width_checks
        + report.transient_width_checks
        + report.post_swap_width_checks
        + report.boundary_width_checks;
    assert_eq!(
        report.width_context_observations * 5,
        total_width_checks,
        "Q949 width-context observation count drift"
    );
    assert_eq!(
        report.width_context_observations,
        report.factors_checked * Q949_WIDTH_CONTEXTS_PER_FACTOR,
        "Q949 width-context coverage drift"
    );
    assert_eq!(
        report.clz_context_observations * 4,
        report.clz_window_checks,
        "Q949 CLZ-context observation count drift"
    );
    assert_eq!(
        report.width_misses,
        report.entry_width_misses
            + report.transient_width_misses
            + report.post_swap_width_misses
            + report.boundary_width_misses,
        "Q949 width phase totals drift"
    );
    assert!(
        report.width_misses <= report.entry_width_checks * 3 + report.boundary_width_checks
    );
    assert!(report.clz_window_misses <= report.clz_window_checks);

    let expected_width_keys = q949_expected_width_context_keys();
    let actual_width_keys = report
        .joint_width_demands
        .iter()
        .map(|demand| (demand.direction, demand.phase, demand.row))
        .collect::<Vec<_>>();
    assert_eq!(
        actual_width_keys, expected_width_keys,
        "Q949 joint-width contexts are incomplete or not deterministically sorted"
    );
    assert_eq!(
        report.joint_width_demands.len(),
        Q949_WIDTH_CONTEXTS_PER_FACTOR
    );
    assert_eq!(
        report
            .joint_width_demands
            .iter()
            .map(|demand| demand.observation_count)
            .sum::<usize>(),
        report.width_context_observations,
        "Q949 joint-width observation total drift"
    );
    for demand in &report.joint_width_demands {
        assert_eq!(
            demand.observation_count, report.factors_checked,
            "Q949 joint-width context omitted a factor"
        );
        assert_eq!(
            demand.available_widths,
            q949_effective_reg_widths(demand.row)
        );
        assert_eq!(
            demand.minimum_fixed_capacity_sum,
            demand.required_widths.iter().sum::<usize>()
        );
        assert_eq!(
            demand.slack_vs_target_683,
            inversion::shrunken_pz_schedule::Q949_TARGET_SUM as i64
                - demand.minimum_fixed_capacity_sum as i64
        );
        for register in 0..5 {
            assert!(demand.required_widths[register] > 0);
            let witness = demand.component_witnesses[register];
            assert!(witness.draw < report.requested_draws);
            q949_factor_label_index(witness.factor_label);
            assert!(!witness.factor.is_zero());
            assert_eq!(
                witness.required_widths[register], demand.required_widths[register],
                "Q949 component maximum witness drift"
            );
        }
    }

    assert!(report.factors_checked > 0);
    assert!(report.earliest_terminal_row <= report.latest_terminal_row);
    assert!(report.latest_terminal_row < SHRUNKEN_PZ_NSTEPS);
    let expected_clz_keys = q949_expected_clz_low_keys(report.latest_terminal_row);
    let actual_clz_keys = report
        .clz_low_bounds
        .iter()
        .map(|bound| {
            (
                bound.direction,
                bound.row,
                q949_clz_register_index(bound.register),
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(
        actual_clz_keys, expected_clz_keys,
        "Q949 CLZ low-bound contexts are incomplete or not deterministically sorted"
    );
    assert_eq!(
        report
            .clz_low_bounds
            .iter()
            .map(|bound| bound.observation_count)
            .sum::<usize>(),
        report.clz_window_checks,
        "Q949 CLZ low-bound observation total drift"
    );
    for bounds in report.clz_low_bounds.chunks_exact(4) {
        assert!(bounds[0].observation_count > 0);
        for bound in &bounds[1..] {
            assert_eq!(bound.observation_count, bounds[0].observation_count);
        }
    }
    let clz_direction_span = (report.latest_terminal_row + 1) * 4;
    for index in 0..clz_direction_span {
        assert_eq!(
            report.clz_low_bounds[index].observation_count,
            report.clz_low_bounds[index + clz_direction_span].observation_count,
            "Q949 forward/reverse CLZ coverage drift"
        );
    }
    for bound in &report.clz_low_bounds {
        let register = q949_clz_register_index(bound.register);
        assert_eq!(
            bound.nonzero_observation_count + bound.zero_observation_count,
            bound.observation_count
        );
        assert_eq!(bound.current_low, q949_effective_reg_los(bound.row)[register]);
        assert_eq!(
            bound.available_width,
            q949_effective_reg_widths(bound.row)[register]
        );
        match (
            bound.minimum_nonzero_observed_width,
            bound.safe_low_upper_bound,
            bound.minimum_witness,
        ) {
            (Some(minimum), Some(safe_low), Some(witness)) => {
                assert!(bound.nonzero_observation_count > 0);
                assert!(minimum > 0);
                assert_eq!(safe_low, minimum - 1);
                assert_eq!(witness.observed_width, minimum);
                assert_eq!(witness.observed_widths[register], minimum);
                assert!(witness.draw < report.requested_draws);
                q949_factor_label_index(witness.factor_label);
                assert!(!witness.factor.is_zero());
            }
            (None, None, None) => assert_eq!(bound.nonzero_observation_count, 0),
            _ => panic!("Q949 CLZ minimum/witness presence drift"),
        }
    }

    let phase_totals = q949_width_phase_miss_totals(report);
    assert_eq!(
        phase_totals,
        [
            report.entry_width_misses,
            report.transient_width_misses,
            report.post_swap_width_misses,
            report.boundary_width_misses,
        ],
        "Q949 width landscape phase totals drift"
    );
    assert_eq!(
        report
            .width_miss_buckets
            .iter()
            .map(|bucket| bucket.miss_count)
            .sum::<usize>(),
        report.width_misses,
        "Q949 width landscape total drift"
    );
    assert_eq!(
        report.first_width_miss.is_some(),
        report.width_misses != 0,
        "Q949 first width miss presence drift"
    );

    for pair in report.width_miss_buckets.windows(2) {
        let lhs = (
            pair[0].direction,
            pair[0].phase,
            pair[0].row,
            q949_width_register_index(pair[0].register),
        );
        let rhs = (
            pair[1].direction,
            pair[1].phase,
            pair[1].row,
            q949_width_register_index(pair[1].register),
        );
        assert!(
            lhs < rhs,
            "Q949 width landscape keys are not unique and sorted"
        );
    }
    for bucket in &report.width_miss_buckets {
        let register = q949_width_register_index(bucket.register);
        assert!(bucket.row < SHRUNKEN_PZ_NSTEPS);
        assert!(bucket.miss_count > 0);
        assert!(bucket.min_observed_width <= bucket.max_observed_width);
        assert!(bucket.min_observed_width > bucket.available_width);
        assert_eq!(
            bucket.available_width,
            q949_effective_reg_widths(bucket.row)[register]
        );
        assert_eq!(
            bucket.max_excess,
            bucket.max_observed_width - bucket.available_width
        );
        if bucket.phase != Q949WidthPhase::Boundary {
            assert_eq!(bucket.direction, Q949TraceDirection::Forward);
        }
    }
    for demand in &report.joint_width_demands {
        for register in 0..5 {
            let has_miss_bucket = report.width_miss_buckets.iter().any(|bucket| {
                bucket.direction == demand.direction
                    && bucket.phase == demand.phase
                    && bucket.row == demand.row
                    && q949_width_register_index(bucket.register) == register
            });
            assert_eq!(
                demand.required_widths[register] > demand.available_widths[register],
                has_miss_bucket,
                "Q949 dense/sparse width landscape drift"
            );
        }
    }
    if let Some(first) = report.first_width_miss {
        assert!(first.draw < report.requested_draws);
        assert!(matches!(first.factor_label, "dx" | "qx_minus_rx"));
        assert!(!first.factor.is_zero());
        let coordinate = first.coordinate;
        assert!(report.width_miss_buckets.iter().any(|bucket| {
            bucket.direction == coordinate.direction
                && bucket.phase == coordinate.phase
                && bucket.row == coordinate.row
                && bucket.register == coordinate.register
                && (bucket.min_observed_width..=bucket.max_observed_width)
                    .contains(&coordinate.observed_width)
                && bucket.available_width == coordinate.available_width
        }));
    }

    assert_eq!(
        report
            .clz_window_miss_buckets
            .iter()
            .map(|bucket| bucket.miss_count)
            .sum::<usize>(),
        report.clz_window_misses,
        "Q949 CLZ landscape total drift"
    );
    assert_eq!(
        report.first_clz_window_miss.is_some(),
        report.clz_window_misses != 0,
        "Q949 first CLZ miss presence drift"
    );
    for pair in report.clz_window_miss_buckets.windows(2) {
        let lhs = (
            pair[0].direction,
            pair[0].row,
            q949_clz_register_index(pair[0].register),
        );
        let rhs = (
            pair[1].direction,
            pair[1].row,
            q949_clz_register_index(pair[1].register),
        );
        assert!(
            lhs < rhs,
            "Q949 CLZ landscape keys are not unique and sorted"
        );
    }
    for bucket in &report.clz_window_miss_buckets {
        let register = q949_clz_register_index(bucket.register);
        assert!(bucket.row < SHRUNKEN_PZ_NSTEPS);
        assert!(bucket.miss_count > 0);
        assert!(bucket.min_observed_width > 0);
        assert!(bucket.min_observed_width <= bucket.max_observed_width);
        assert!(bucket.max_observed_width <= bucket.low);
        assert_eq!(bucket.low, q949_effective_reg_los(bucket.row)[register]);
        assert_eq!(
            bucket.available_width,
            q949_effective_reg_widths(bucket.row)[register]
        );
        assert_eq!(
            bucket.max_shortfall,
            bucket.low + 1 - bucket.min_observed_width
        );
    }
    for bound in &report.clz_low_bounds {
        let has_miss_bucket = report.clz_window_miss_buckets.iter().any(|bucket| {
            bucket.direction == bound.direction
                && bucket.row == bound.row
                && bucket.register == bound.register
        });
        assert_eq!(
            !q949_clz_low_is_safe(bound),
            has_miss_bucket,
            "Q949 dense/sparse CLZ landscape drift"
        );
    }
    if let Some(first) = report.first_clz_window_miss {
        assert!(first.draw < report.requested_draws);
        assert!(matches!(first.factor_label, "dx" | "qx_minus_rx"));
        assert!(!first.factor.is_zero());
        let coordinate = first.coordinate;
        assert!(report.clz_window_miss_buckets.iter().any(|bucket| {
            bucket.direction == coordinate.direction
                && bucket.row == coordinate.row
                && bucket.register == coordinate.register
                && (bucket.min_observed_width..=bucket.max_observed_width)
                    .contains(&coordinate.observed_width)
                && bucket.low == coordinate.low
                && bucket.available_width == coordinate.available_width
        }));
    }
}

fn q949_census_status(report: &Q949SupportedTraceReport) -> &'static str {
    if report.requested_draws == Q949_PROOF_DRAWS
        && report.accepted_shots == Q949_PROOF_DRAWS
        && report.rejected_draws == 0
        && report.factors_checked == Q949_PROOF_FACTORS
        && report.width_misses == 0
        && report.clz_window_misses == 0
        && report.earliest_terminal_row <= report.latest_terminal_row
        && report.latest_terminal_row < inversion::shrunken_pz_schedule::SHRUNKEN_PZ_NSTEPS
        && report.max_counter
            == inversion::shrunken_pz_schedule::SHRUNKEN_PZ_NSTEPS
                - report.earliest_terminal_row
    {
        "pass"
    } else {
        "reject"
    }
}

fn q949_assert_diagnostic_identity(label: &str, identity: &str) {
    assert_eq!(identity.len(), 64, "Q949 diagnostic {label} length drift");
    assert!(
        identity
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte)),
        "Q949 diagnostic {label} is not lowercase hexadecimal"
    );
}

fn q949_census_diagnostics_json(
    identity: &Q949ProofIdentity,
    report: &Q949SupportedTraceReport,
) -> String {
    for (label, value) in [
        ("source_id", identity.source_id.as_str()),
        ("schedule_id", identity.schedule_id.as_str()),
        ("route_id", identity.route_id.as_str()),
        ("op_stream_id", identity.op_stream_id.as_str()),
    ] {
        q949_assert_diagnostic_identity(label, value);
    }

    let phase_totals = q949_width_phase_miss_totals(report);
    let status = q949_census_status(report);
    let coverage_complete = report.requested_draws == Q949_PROOF_DRAWS
        && report.accepted_shots == Q949_PROOF_DRAWS
        && report.rejected_draws == 0
        && report.factors_checked == Q949_PROOF_FACTORS;
    let passing_width_contexts = report
        .joint_width_demands
        .iter()
        .filter(|demand| q949_joint_width_fits_current(demand))
        .count();
    let safe_clz_low_contexts = report
        .clz_low_bounds
        .iter()
        .filter(|bound| q949_clz_low_is_safe(bound))
        .count();
    let mut output = String::with_capacity(
        2_048 + report.joint_width_demands.len() * 1_024
            + report.clz_low_bounds.len() * 384
            + report.width_miss_buckets.len() * 192
            + report.clz_window_miss_buckets.len() * 176,
    );
    write!(
        &mut output,
        concat!(
            "{{\"schema\":\"{}\",",
            "\"generator_model\":\"GPT-Codex\",\"source_bound\":true,",
            "\"identity\":{{\"source_id\":\"{}\",\"schedule_id\":\"{}\",",
            "\"route_id\":\"{}\",\"op_stream_id\":\"{}\",",
            "\"tail_nonce\":{},\"tail_nonce_bits\":{},\"op_count\":{}}},",
            "\"coverage\":{{\"requested_draws\":{},\"accepted_draws\":{},",
            "\"rejected_draws\":{},\"factors_checked\":{},",
            "\"expected_draws\":{},\"expected_factors\":{},",
            "\"complete\":{},",
            "\"forward_rows_checked\":{},\"backward_rows_checked\":{},",
            "\"row_bounds_checked\":{},\"entry_width_checks\":{},",
            "\"transient_width_checks\":{},\"post_swap_width_checks\":{},",
            "\"boundary_width_checks\":{},",
            "\"width_context_observations\":{},",
            "\"joint_width_contexts\":{},\"expected_joint_width_contexts\":{},",
            "\"clz_window_checks\":{},\"clz_context_observations\":{},",
            "\"clz_low_bound_contexts\":{}}},"
        ),
        Q949_CENSUS_DIAGNOSTICS_SCHEMA,
        identity.source_id,
        identity.schedule_id,
        identity.route_id,
        identity.op_stream_id,
        identity.tail_nonce,
        identity.tail_nonce_bits,
        identity.op_count,
        report.requested_draws,
        report.accepted_shots,
        report.rejected_draws,
        report.factors_checked,
        Q949_PROOF_DRAWS,
        Q949_PROOF_FACTORS,
        coverage_complete,
        report.forward_rows_checked,
        report.backward_rows_checked,
        report.row_bounds_checked,
        report.entry_width_checks,
        report.transient_width_checks,
        report.post_swap_width_checks,
        report.boundary_width_checks,
        report.width_context_observations,
        report.joint_width_demands.len(),
        Q949_WIDTH_CONTEXTS_PER_FACTOR,
        report.clz_window_checks,
        report.clz_context_observations,
        report.clz_low_bounds.len(),
    )
    .expect("write Q949 diagnostic coverage JSON");

    write!(
        &mut output,
        concat!(
            "\"joint_width_demand\":{{",
            "\"key\":[\"direction\",\"phase\",\"row\"],",
            "\"encoding\":\"dense_all_observed_contexts\",",
            "\"direction_order\":[\"forward\",\"reverse\"],",
            "\"phase_order\":[\"entry\",\"transient\",\"post_swap\",\"boundary\"],",
            "\"register_order\":[\"A\",\"B\",\"ca\",\"cb\",\"q\"],",
            "\"factor_order\":[\"dx\",\"qx_minus_rx\"],",
            "\"witness_tie_break\":\"lowest_draw_then_factor_order\",",
            "\"target_sum\":683,\"context_count\":{},",
            "\"expected_context_count\":{},\"observation_count\":{},",
            "\"passing_contexts\":{},\"failing_contexts\":{},\"contexts\":["
        ),
        report.joint_width_demands.len(),
        Q949_WIDTH_CONTEXTS_PER_FACTOR,
        report.width_context_observations,
        passing_width_contexts,
        report.joint_width_demands.len() - passing_width_contexts,
    )
    .expect("write Q949 joint-width diagnostic header");
    for (context_index, demand) in report.joint_width_demands.iter().enumerate() {
        if context_index != 0 {
            output.push(',');
        }
        write!(
            &mut output,
            concat!(
                "{{\"direction\":\"{}\",\"phase\":\"{}\",\"row\":{},",
                "\"observation_count\":{},",
                "\"required_widths\":[{},{},{},{},{}],",
                "\"available_widths\":[{},{},{},{},{}],",
                "\"minimum_fixed_capacity_sum\":{},",
                "\"slack_vs_target_683\":{},\"fits_current\":{},",
                "\"component_witnesses\":["
            ),
            q949_direction_label(demand.direction),
            q949_width_phase_label(demand.phase),
            demand.row,
            demand.observation_count,
            demand.required_widths[0],
            demand.required_widths[1],
            demand.required_widths[2],
            demand.required_widths[3],
            demand.required_widths[4],
            demand.available_widths[0],
            demand.available_widths[1],
            demand.available_widths[2],
            demand.available_widths[3],
            demand.available_widths[4],
            demand.minimum_fixed_capacity_sum,
            demand.slack_vs_target_683,
            q949_joint_width_fits_current(demand),
        )
        .expect("write Q949 joint-width diagnostic context");
        for (register, witness) in demand.component_witnesses.iter().enumerate() {
            if register != 0 {
                output.push(',');
            }
            write!(
                &mut output,
                concat!(
                    "{{\"register\":\"{}\",\"required_width\":{},",
                    "\"draw\":{},\"factor_label\":\"{}\",",
                    "\"factor\":\"0x{:064x}\",",
                    "\"required_widths\":[{},{},{},{},{}]}}"
                ),
                Q949_WIDTH_REGISTERS[register],
                demand.required_widths[register],
                witness.draw,
                witness.factor_label,
                witness.factor,
                witness.required_widths[0],
                witness.required_widths[1],
                witness.required_widths[2],
                witness.required_widths[3],
                witness.required_widths[4],
            )
            .expect("write Q949 joint-width maximum witness");
        }
        output.push_str("]}");
    }
    output.push_str("]},");

    write!(
        &mut output,
        concat!(
            "\"width_misses\":{{\"key\":[\"direction\",\"phase\",\"row\",\"register\"],",
            "\"encoding\":\"sparse_nonzero\",",
            "\"register_order\":[\"A\",\"B\",\"ca\",\"cb\",\"q\"],",
            "\"total\":{},\"bucket_count\":{},",
            "\"phase_totals\":{{\"entry\":{},\"transient\":{},",
            "\"post_swap\":{},\"boundary\":{}}},\"buckets\":["
        ),
        report.width_misses,
        report.width_miss_buckets.len(),
        phase_totals[0],
        phase_totals[1],
        phase_totals[2],
        phase_totals[3],
    )
    .expect("write Q949 width diagnostic header");
    for (index, bucket) in report.width_miss_buckets.iter().enumerate() {
        if index != 0 {
            output.push(',');
        }
        write!(
            &mut output,
            concat!(
                "{{\"direction\":\"{}\",\"phase\":\"{}\",\"row\":{},",
                "\"register\":\"{}\",\"count\":{},",
                "\"observed_width_min\":{},\"observed_width_max\":{},",
                "\"available_width\":{},\"max_excess\":{}}}"
            ),
            q949_direction_label(bucket.direction),
            q949_width_phase_label(bucket.phase),
            bucket.row,
            bucket.register,
            bucket.miss_count,
            bucket.min_observed_width,
            bucket.max_observed_width,
            bucket.available_width,
            bucket.max_excess,
        )
        .expect("write Q949 width diagnostic bucket");
    }
    output.push_str("],\"first_miss\":");
    if let Some(first) = report.first_width_miss {
        let coordinate = first.coordinate;
        write!(
            &mut output,
            concat!(
                "{{\"draw\":{},\"factor_label\":\"{}\",\"factor\":\"0x{:064x}\",",
                "\"direction\":\"{}\",\"phase\":\"{}\",\"row\":{},",
                "\"register\":\"{}\",\"observed_width\":{},",
                "\"available_width\":{},",
                "\"observed_widths\":[{},{},{},{},{}],",
                "\"available_widths\":[{},{},{},{},{}]}}"
            ),
            first.draw,
            first.factor_label,
            first.factor,
            q949_direction_label(coordinate.direction),
            q949_width_phase_label(coordinate.phase),
            coordinate.row,
            coordinate.register,
            coordinate.observed_width,
            coordinate.available_width,
            coordinate.observed_widths[0],
            coordinate.observed_widths[1],
            coordinate.observed_widths[2],
            coordinate.observed_widths[3],
            coordinate.observed_widths[4],
            coordinate.available_widths[0],
            coordinate.available_widths[1],
            coordinate.available_widths[2],
            coordinate.available_widths[3],
            coordinate.available_widths[4],
        )
        .expect("write Q949 first width miss");
    } else {
        output.push_str("null");
    }
    output.push_str("},");

    write!(
        &mut output,
        concat!(
            "\"clz_low_bounds\":{{\"key\":[\"direction\",\"row\",\"register\"],",
            "\"encoding\":\"dense_all_observed_contexts\",",
            "\"direction_order\":[\"forward\",\"reverse\"],",
            "\"register_order\":[\"A\",\"B\",\"ca\",\"cb\"],",
            "\"factor_order\":[\"dx\",\"qx_minus_rx\"],",
            "\"witness_tie_break\":\"lowest_draw_then_factor_order\",",
            "\"zero_width_semantics\":\"no_low_constraint\",",
            "\"context_count\":{},\"expected_context_count\":{},",
            "\"observation_count\":{},",
            "\"safe_contexts\":{},\"unsafe_contexts\":{},\"contexts\":["
        ),
        report.clz_low_bounds.len(),
        2 * (report.latest_terminal_row + 1) * 4,
        report.clz_window_checks,
        safe_clz_low_contexts,
        report.clz_low_bounds.len() - safe_clz_low_contexts,
    )
    .expect("write Q949 CLZ low-bound diagnostic header");
    for (context_index, bound) in report.clz_low_bounds.iter().enumerate() {
        if context_index != 0 {
            output.push(',');
        }
        write!(
            &mut output,
            concat!(
                "{{\"direction\":\"{}\",\"row\":{},\"register\":\"{}\",",
                "\"observation_count\":{},\"nonzero_observation_count\":{},",
                "\"zero_observation_count\":{},",
                "\"minimum_nonzero_observed_width\":"
            ),
            q949_direction_label(bound.direction),
            bound.row,
            bound.register,
            bound.observation_count,
            bound.nonzero_observation_count,
            bound.zero_observation_count,
        )
        .expect("write Q949 CLZ low-bound diagnostic context");
        if let Some(minimum) = bound.minimum_nonzero_observed_width {
            write!(&mut output, "{minimum}").expect("write Q949 CLZ minimum");
        } else {
            output.push_str("null");
        }
        output.push_str(",\"safe_low_upper_bound\":");
        if let Some(safe_low) = bound.safe_low_upper_bound {
            write!(&mut output, "{safe_low}").expect("write Q949 safe CLZ low");
        } else {
            output.push_str("null");
        }
        write!(
            &mut output,
            concat!(
                ",\"current_low\":{},\"available_width\":{},",
                "\"current_low_safe\":{},\"minimum_witness\":"
            ),
            bound.current_low,
            bound.available_width,
            q949_clz_low_is_safe(bound),
        )
        .expect("write Q949 CLZ low-bound values");
        if let Some(witness) = bound.minimum_witness {
            write!(
                &mut output,
                concat!(
                    "{{\"draw\":{},\"factor_label\":\"{}\",",
                    "\"factor\":\"0x{:064x}\",\"observed_width\":{},",
                    "\"observed_widths\":[{},{},{},{}]}}"
                ),
                witness.draw,
                witness.factor_label,
                witness.factor,
                witness.observed_width,
                witness.observed_widths[0],
                witness.observed_widths[1],
                witness.observed_widths[2],
                witness.observed_widths[3],
            )
            .expect("write Q949 CLZ minimum witness");
        } else {
            output.push_str("null");
        }
        output.push('}');
    }
    output.push_str("]},");

    write!(
        &mut output,
        concat!(
            "\"clz_window_misses\":{{\"key\":[\"direction\",\"row\",\"register\"],",
            "\"encoding\":\"sparse_nonzero\",",
            "\"register_order\":[\"A\",\"B\",\"ca\",\"cb\"],",
            "\"total\":{},\"bucket_count\":{},\"buckets\":["
        ),
        report.clz_window_misses,
        report.clz_window_miss_buckets.len(),
    )
    .expect("write Q949 CLZ diagnostic header");
    for (index, bucket) in report.clz_window_miss_buckets.iter().enumerate() {
        if index != 0 {
            output.push(',');
        }
        write!(
            &mut output,
            concat!(
                "{{\"direction\":\"{}\",\"row\":{},\"register\":\"{}\",",
                "\"count\":{},\"observed_width_min\":{},",
                "\"observed_width_max\":{},\"low\":{},",
                "\"available_width\":{},\"max_shortfall\":{}}}"
            ),
            q949_direction_label(bucket.direction),
            bucket.row,
            bucket.register,
            bucket.miss_count,
            bucket.min_observed_width,
            bucket.max_observed_width,
            bucket.low,
            bucket.available_width,
            bucket.max_shortfall,
        )
        .expect("write Q949 CLZ diagnostic bucket");
    }
    output.push_str("],\"first_miss\":");
    if let Some(first) = report.first_clz_window_miss {
        let coordinate = first.coordinate;
        write!(
            &mut output,
            concat!(
                "{{\"draw\":{},\"factor_label\":\"{}\",\"factor\":\"0x{:064x}\",",
                "\"direction\":\"{}\",\"row\":{},\"register\":\"{}\",",
                "\"observed_width\":{},\"low\":{},\"available_width\":{}}}"
            ),
            first.draw,
            first.factor_label,
            first.factor,
            q949_direction_label(coordinate.direction),
            coordinate.row,
            coordinate.register,
            coordinate.observed_width,
            coordinate.low,
            coordinate.available_width,
        )
        .expect("write Q949 first CLZ miss");
    } else {
        output.push_str("null");
    }
    write!(
        &mut output,
        concat!(
            "}},\"reverse_ca255_relation\":{{\"row\":380,",
            "\"checks\":{},\"active_checks\":{},\"inactive_checks\":{},",
            "\"failures\":{},",
            "\"identity\":\"ca[255]=NOT(active_and_ca_lt_cb)\"}},",
            "\"terminal\":{{\"full_ca_checks\":{},",
            "\"earliest_row\":{},\"latest_row\":{},\"max_counter\":{}}},",
            "\"proof_status\":\"{}\"}}\n"
        ),
        report.reverse_row_380_relation_checks,
        report.reverse_row_380_active_checks,
        report.reverse_row_380_inactive_checks,
        report.reverse_row_380_relation_failures,
        report.terminal_full_ca_checks,
        report.earliest_terminal_row,
        report.latest_terminal_row,
        report.max_counter,
        status,
    )
    .expect("write Q949 diagnostic trailer");

    assert!(
        output.starts_with("{\"schema\":\"q949-census-diagnostics-v2\","),
        "Q949 diagnostic schema prefix drift"
    );
    assert!(output.ends_with("}\n"), "Q949 diagnostic JSON trailer drift");
    output
}

fn q949_emit_census_diagnostics(
    path: &str,
    identity: &Q949ProofIdentity,
    report: &Q949SupportedTraceReport,
) {
    q949_validate_census_report(report);
    let phase_totals = q949_width_phase_miss_totals(report);
    let status = q949_census_status(report);
    eprintln!(
        concat!(
            "Q949_CENSUS_SUMMARY schema={} status={} draws={}/{} rejected={} ",
            "factors={}/{} width_contexts={} width_observations={} ",
            "width_misses={} width_buckets={} ",
            "width_by_phase=entry:{},transient:{},post_swap:{},boundary:{} ",
            "clz_contexts={} clz_observations={} clz_misses={} clz_buckets={}"
        ),
        Q949_CENSUS_DIAGNOSTICS_SCHEMA,
        status,
        report.accepted_shots,
        report.requested_draws,
        report.rejected_draws,
        report.factors_checked,
        Q949_PROOF_FACTORS,
        report.joint_width_demands.len(),
        report.width_context_observations,
        report.width_misses,
        report.width_miss_buckets.len(),
        phase_totals[0],
        phase_totals[1],
        phase_totals[2],
        phase_totals[3],
        report.clz_low_bounds.len(),
        report.clz_context_observations,
        report.clz_window_misses,
        report.clz_window_miss_buckets.len(),
    );
    if let Some(first) = report.first_width_miss {
        eprintln!(
            concat!(
                "Q949_CENSUS_FIRST_WIDTH draw={} factor={} factor_value=0x{:064x} ",
                "direction={} phase={} row={} register={} observed={} available={}"
            ),
            first.draw,
            first.factor_label,
            first.factor,
            q949_direction_label(first.coordinate.direction),
            q949_width_phase_label(first.coordinate.phase),
            first.coordinate.row,
            first.coordinate.register,
            first.coordinate.observed_width,
            first.coordinate.available_width,
        );
    }
    if let Some(first) = report.first_clz_window_miss {
        eprintln!(
            concat!(
                "Q949_CENSUS_FIRST_CLZ draw={} factor={} factor_value=0x{:064x} ",
                "direction={} row={} register={} observed={} low={} available={}"
            ),
            first.draw,
            first.factor_label,
            first.factor,
            q949_direction_label(first.coordinate.direction),
            first.coordinate.row,
            first.coordinate.register,
            first.coordinate.observed_width,
            first.coordinate.low,
            first.coordinate.available_width,
        );
    }

    let json = q949_census_diagnostics_json(identity, report);
    std::fs::write(path, &json)
        .unwrap_or_else(|error| panic!("write Q949 census diagnostics to {path}: {error}"));
    eprintln!(
        "Q949_CENSUS_DIAGNOSTICS schema={} path={} bytes={}",
        Q949_CENSUS_DIAGNOSTICS_SCHEMA,
        path,
        json.len()
    );
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct Q945SiteAccumulator {
    checks: usize,
    zero_entry_checks: usize,
    restoration_checks: usize,
}

fn q945_phase_for_factor(factor_label: &str) -> Q945SupportPhase {
    match factor_label {
        "dx" => Q945SupportPhase::InvFwd,
        "qx_minus_rx" => Q945SupportPhase::AltCancel,
        other => panic!("unknown Q945 support factor label: {other}"),
    }
}

/// Bind all 56 prospective persistent-gate host sites to the exact committed
/// Fiat-Shamir support. This is a feasibility census only: it does not enable
/// the Q944 route or make a semantic correctness claim.
#[doc(hidden)]
pub fn q944_gate_host_feasibility_check(
    builder: &crate::point_add::B,
) -> inversion::q944_gate_host_feasibility::Q944GateHostFeasibilityReport {
    use inversion::q944_gate_host_feasibility::Q944GateHostCensus;

    assert_eq!(
        std::env::var("LOWQ_Q945_LOCAL_HOSTS").ok().as_deref(),
        Some("1"),
        "Q944 gate-host census requires the Q945 local-host route"
    );
    assert_eq!(
        std::env::var("LOWQ_Q945_DIRTY_PARITY_ARITHMETIC")
            .ok()
            .as_deref(),
        Some("1"),
        "Q944 gate-host census requires the proved dirty-parity arithmetic"
    );
    assert_ne!(
        std::env::var("LOWQ_Q944_GATE_HOSTS").ok().as_deref(),
        Some("1"),
        "Q944 feasibility census must precede broad integration"
    );

    let requested_draws = env_usize("TRAILMIX_Q944_GATE_HOST_SHOTS", TRAILMIX_NUM_TESTS);
    assert_eq!(
        requested_draws, Q949_PROOF_DRAWS,
        "Q944 production gate-host census requires all {Q949_PROOF_DRAWS} draws"
    );
    let mut census = Q944GateHostCensus::new(requested_draws);
    let mut xof = q949_builder_fiat_hash(builder).finalize_xof();
    let curve = secp256k1();

    for draw in 0..requested_draws {
        let mut random = [[0u8; 32]; 2];
        xof.read(&mut random[0]);
        xof.read(&mut random[1]);
        let k1 = U256::from_le_bytes(random[0]);
        let k2 = U256::from_le_bytes(random[1]);
        let target = curve.mul(curve.gx, curve.gy, k1);
        let other = curve.mul(curve.gx, curve.gy, k2);
        if target.0 == other.0
            || (target.0.is_zero() && target.1.is_zero())
            || (other.0.is_zero() && other.1.is_zero())
        {
            census.record_rejected_draw();
            continue;
        }
        census.record_accepted_draw();
        let result = curve.add(target.0, target.1, other.0, other.1);
        for (factor_label, factor) in [
            ("dx", sub_mod_p(target.0, other.0, curve.modulus)),
            (
                "qx_minus_rx",
                sub_mod_p(other.0, result.0, curve.modulus),
            ),
        ] {
            let phase = q945_phase_for_factor(factor_label);
            let certificate = inversion::shrunken_pz_schedule::
                q949_affine_trace_certificate_u256(factor);
            // Preserve inherited semantic-route diagnostics, but do not let
            // them suppress the requested structural 56-site host census.
            census.record_inherited_diagnostics(
                certificate.width_misses,
                certificate.clz_window_misses,
                certificate.narrow_compare_misses,
            );
            census.record_factor(
                draw,
                factor_label,
                factor,
                phase,
                &certificate.q944_gate_call_observations,
            );
        }
    }

    census.finish()
}

/// Replay the exact committed Fiat-Shamir factors through the ideal PZ trace and
/// bind every Q945 local lender/carry to its value at the corresponding call
/// boundary. This returns misses instead of aborting so the first false host is
/// retained as a reproducible counterexample.
#[doc(hidden)]
pub fn q945_host_support_check(builder: &crate::point_add::B) -> Q945HostSupportReport {
    use inversion::q945_local_hosts::{Q945HclzForm, Q945StateRegister, Q945Substep};
    use inversion::shrunken_pz_schedule::Q949NarrowCompareSubstep;

    assert_eq!(
        std::env::var("LOWQ_Q945_LOCAL_HOSTS").ok().as_deref(),
        Some("1"),
        "Q945 host census requires the local-host route"
    );
    assert_ne!(
        std::env::var("LOWQ_Q946_FRESH_SUPPORT_CERTIFIED")
            .ok()
            .as_deref(),
        Some("1"),
        "Q945 host census cannot inherit a Q946 support certificate"
    );
    assert_ne!(
        std::env::var("LOWQ_Q947_FRESH_SUPPORT_CERTIFIED")
            .ok()
            .as_deref(),
        Some("1"),
        "Q945 host census cannot inherit a Q947 support certificate"
    );

    let requested_draws = env_usize("TRAILMIX_Q945_HOST_SUPPORT_SHOTS", TRAILMIX_NUM_TESTS);
    assert_eq!(
        requested_draws, Q949_PROOF_DRAWS,
        "Q945 production host census requires all {Q949_PROOF_DRAWS} draws"
    );
    let mut xof = q949_builder_fiat_hash(builder).finalize_xof();
    let curve = secp256k1();
    let mut accepted_draws = 0usize;
    let mut rejected_draws = 0usize;
    let mut factors_checked = 0usize;
    let mut hclz_host_checks = 0usize;
    let mut hclz_zero_entry_checks = 0usize;
    let mut hclz_restoration_checks = 0usize;
    let mut hclz_host_misses = 0usize;
    let mut carry_host_checks = 0usize;
    let mut carry_zero_entry_checks = 0usize;
    let mut carry_restoration_checks = 0usize;
    let mut carry_host_misses = 0usize;
    let mut carry_semantic_checks = 0usize;
    let mut carry_semantic_misses = 0usize;
    let mut row364_checks = 0usize;
    let mut row364_b80_zero_checks = 0usize;
    let mut row364_identity_checks = 0usize;
    let mut row364_misses = 0usize;
    let mut row374_q24_checks = 0usize;
    let mut row374_q24_zero_checks = 0usize;
    let mut row374_q24_noncarry_touches = 0usize;
    let mut row374_misses = 0usize;
    let mut preterminal_counter_off_checks = 0usize;
    let mut preterminal_counter_off_zero_checks = 0usize;
    let mut preterminal_counter_off_misses = 0usize;
    let mut row385_special_checks = 0usize;
    let mut row385_special_zero_checks = 0usize;
    let mut row385_special_misses = 0usize;
    let mut width_misses = 0usize;
    let mut clz_window_misses = 0usize;
    let mut narrow_compare_checks = 0usize;
    let mut narrow_compare_misses = 0usize;
    let mut division_offset_compare_checks = 0usize;
    let mut division_offset_compare_misses = 0usize;
    let mut multiply_offset_cleanup_compare_checks = 0usize;
    let mut multiply_offset_cleanup_compare_misses = 0usize;
    let mut first_host_miss = None;
    let mut first_narrow_compare_miss = None;
    let mut hclz_sites = BTreeMap::<Q945HclzSupportSite, Q945SiteAccumulator>::new();
    let mut carry_sites = BTreeMap::<Q945CarrySupportSite, Q945SiteAccumulator>::new();
    let mut earliest_terminal_row = usize::MAX;
    let mut latest_terminal_row = 0usize;

    for draw in 0..requested_draws {
        let mut random = [[0u8; 32]; 2];
        xof.read(&mut random[0]);
        xof.read(&mut random[1]);
        let k1 = U256::from_le_bytes(random[0]);
        let k2 = U256::from_le_bytes(random[1]);
        let target = curve.mul(curve.gx, curve.gy, k1);
        let other = curve.mul(curve.gx, curve.gy, k2);
        if target.0 == other.0
            || (target.0.is_zero() && target.1.is_zero())
            || (other.0.is_zero() && other.1.is_zero())
        {
            rejected_draws += 1;
            continue;
        }
        accepted_draws += 1;
        let result = curve.add(target.0, target.1, other.0, other.1);
        for (factor_label, factor) in [
            ("dx", sub_mod_p(target.0, other.0, curve.modulus)),
            (
                "qx_minus_rx",
                sub_mod_p(other.0, result.0, curve.modulus),
            ),
        ] {
            let phase = q945_phase_for_factor(factor_label);
            let certificate = inversion::shrunken_pz_schedule::
                q949_affine_trace_certificate_u256(factor);
            factors_checked += 1;
            width_misses += certificate.width_misses;
            clz_window_misses += certificate.clz_window_misses;
            narrow_compare_checks += certificate.narrow_compare_checks;
            narrow_compare_misses += certificate.narrow_compare_misses;
            division_offset_compare_checks += certificate.division_offset_compare_checks;
            division_offset_compare_misses += certificate.division_offset_compare_misses;
            multiply_offset_cleanup_compare_checks +=
                certificate.multiply_offset_cleanup_compare_checks;
            multiply_offset_cleanup_compare_misses +=
                certificate.multiply_offset_cleanup_compare_misses;
            earliest_terminal_row = earliest_terminal_row.min(certificate.first_terminal_row);
            latest_terminal_row = latest_terminal_row.max(certificate.first_terminal_row);
            if first_narrow_compare_miss.is_none() {
                first_narrow_compare_miss = certificate.first_narrow_compare_miss.map(
                    |coordinate| Q945SupportedNarrowCompareMiss {
                        draw,
                        factor_label,
                        factor,
                        phase,
                        coordinate,
                    },
                );
            }

            for observation in certificate.q945_hclz_host_observations {
                hclz_host_checks += 1;
                let zero = !observation.entry_value;
                let restored = observation.entry_value == observation.exit_value;
                hclz_zero_entry_checks += usize::from(zero);
                hclz_restoration_checks += usize::from(restored);
                let site = Q945HclzSupportSite {
                    phase,
                    direction: observation.direction,
                    row: observation.row,
                    substep: observation.substep,
                    form: observation.form,
                    host: observation.host,
                };
                let site_count = hclz_sites.entry(site).or_default();
                site_count.checks += 1;
                site_count.zero_entry_checks += usize::from(zero);
                site_count.restoration_checks += usize::from(restored);

                let failed = !zero || !restored;
                hclz_host_misses += usize::from(failed);
                if observation.host.register == Q945StateRegister::CounterOff {
                    preterminal_counter_off_checks += 1;
                    preterminal_counter_off_zero_checks += usize::from(zero);
                    preterminal_counter_off_misses += usize::from(failed);
                }
                if observation.row == 385 {
                    row385_special_checks += 1;
                    row385_special_zero_checks += usize::from(zero);
                    row385_special_misses += usize::from(failed);
                }
                if first_host_miss.is_none() && failed {
                    first_host_miss = Some(Q945SupportedHostMiss {
                        draw,
                        factor_label,
                        factor,
                        phase,
                        direction: observation.direction,
                        row: observation.row,
                        substep: observation.substep,
                        form: Some(observation.form),
                        host: observation.host,
                        entry_value: observation.entry_value,
                        exit_value: observation.exit_value,
                        reason: if !zero {
                            "hclz-host-nonzero-on-entry"
                        } else {
                            "hclz-host-not-restored"
                        },
                        expected_lt: None,
                        full_lt: None,
                        route_lt: None,
                        boundary: observation.boundary,
                    });
                }
            }

            for observation in certificate.q945_carry_host_observations {
                carry_host_checks += 1;
                carry_semantic_checks += 1;
                let zero = !observation.entry_value;
                let restored = observation.entry_value == observation.exit_value;
                let semantic = observation.boundary_reconstructed
                    && observation.expected_lt == observation.full_lt
                    && observation.expected_lt == observation.route_lt;
                carry_zero_entry_checks += usize::from(zero);
                carry_restoration_checks += usize::from(restored);
                carry_semantic_misses += usize::from(!semantic);
                let host_failed = !zero || !restored;
                carry_host_misses += usize::from(host_failed);
                let site = Q945CarrySupportSite {
                    phase,
                    direction: observation.direction,
                    row: observation.row,
                    substep: observation.substep,
                    host: observation.host,
                };
                let site_count = carry_sites.entry(site).or_default();
                site_count.checks += 1;
                site_count.zero_entry_checks += usize::from(zero);
                site_count.restoration_checks += usize::from(restored);

                let row364 = observation.row == 364
                    && observation.substep == Q945Substep::Division;
                if row364 {
                    row364_checks += 1;
                    row364_b80_zero_checks += usize::from(zero);
                    row364_identity_checks += usize::from(semantic);
                    row364_misses += usize::from(host_failed || !semantic);
                }
                let row374 = observation.row == 374
                    && observation.substep == Q945Substep::Division;
                if row374 {
                    assert_eq!(observation.host.register, Q945StateRegister::Q);
                    assert_eq!(observation.host.bit, 24);
                    row374_q24_checks += 1;
                    row374_q24_zero_checks += usize::from(zero);
                    row374_q24_noncarry_touches += observation.q24_noncarry_touches;
                    row374_misses += usize::from(
                        host_failed || !semantic || observation.q24_noncarry_touches != 0,
                    );
                }

                let failed = host_failed || !semantic || (row374 && observation.q24_noncarry_touches != 0);
                if first_host_miss.is_none() && failed {
                    let reason = if !observation.boundary_reconstructed {
                        "carry-boundary-not-reconstructed"
                    } else if !zero {
                        if row364 {
                            "row364-b80-nonzero"
                        } else if row374 {
                            "row374-q24-nonzero"
                        } else {
                            "carry-host-nonzero-on-entry"
                        }
                    } else if !restored {
                        "carry-host-not-restored"
                    } else if observation.expected_lt != observation.full_lt {
                        "full-comparison-mismatch"
                    } else if observation.expected_lt != observation.route_lt {
                        if row364 {
                            "row364-lower80-not-a80-identity-mismatch"
                        } else {
                            "narrow-comparison-mismatch"
                        }
                    } else {
                        "row374-q24-noncarry-touch"
                    };
                    first_host_miss = Some(Q945SupportedHostMiss {
                        draw,
                        factor_label,
                        factor,
                        phase,
                        direction: observation.direction,
                        row: observation.row,
                        substep: observation.substep,
                        form: None,
                        host: observation.host,
                        entry_value: observation.entry_value,
                        exit_value: observation.exit_value,
                        reason,
                        expected_lt: Some(observation.expected_lt),
                        full_lt: Some(observation.full_lt),
                        route_lt: Some(observation.route_lt),
                        boundary: observation.boundary,
                    });
                }
            }
        }
    }

    assert_eq!(accepted_draws + rejected_draws, requested_draws);
    assert_eq!(factors_checked, 2 * accepted_draws);
    assert_eq!(hclz_sites.len(), 208, "Q945 HCLZ site coverage drift");
    assert_eq!(carry_sites.len(), 56, "Q945 carry site coverage drift");
    assert!(hclz_sites.values().all(|site| site.checks == accepted_draws));
    assert!(carry_sites.values().all(|site| site.checks == accepted_draws));
    assert_eq!(hclz_host_checks, 208 * accepted_draws);
    assert_eq!(carry_host_checks, 56 * accepted_draws);
    assert_eq!(row364_checks, 4 * accepted_draws);
    assert_eq!(row374_q24_checks, 4 * accepted_draws);
    assert_eq!(preterminal_counter_off_checks, 104 * accepted_draws);
    assert_eq!(row385_special_checks, 16 * accepted_draws);
    assert_eq!(
        narrow_compare_checks,
        division_offset_compare_checks + multiply_offset_cleanup_compare_checks
    );
    assert_eq!(
        narrow_compare_misses,
        division_offset_compare_misses + multiply_offset_cleanup_compare_misses
    );
    assert_eq!(first_host_miss.is_some(), hclz_host_misses + carry_host_misses + carry_semantic_misses != 0 || row374_q24_noncarry_touches != 0);
    assert_eq!(first_narrow_compare_miss.is_some(), narrow_compare_misses != 0);

    let support_clean = rejected_draws == 0
        && width_misses == 0
        && clz_window_misses == 0
        && narrow_compare_misses == 0
        && hclz_host_misses == 0
        && carry_host_misses == 0
        && carry_semantic_misses == 0
        && row364_misses == 0
        && row374_misses == 0
        && preterminal_counter_off_misses == 0
        && row385_special_misses == 0;
    Q945HostSupportReport {
        requested_draws,
        accepted_draws,
        rejected_draws,
        factors_checked,
        hclz_host_checks,
        hclz_zero_entry_checks,
        hclz_restoration_checks,
        hclz_host_misses,
        carry_host_checks,
        carry_zero_entry_checks,
        carry_restoration_checks,
        carry_host_misses,
        carry_semantic_checks,
        carry_semantic_misses,
        row364_checks,
        row364_b80_zero_checks,
        row364_identity_checks,
        row364_misses,
        row374_q24_checks,
        row374_q24_zero_checks,
        row374_q24_noncarry_touches,
        row374_misses,
        preterminal_counter_off_checks,
        preterminal_counter_off_zero_checks,
        preterminal_counter_off_misses,
        row385_special_checks,
        row385_special_zero_checks,
        row385_special_misses,
        width_misses,
        clz_window_misses,
        narrow_compare_checks,
        narrow_compare_misses,
        division_offset_compare_checks,
        division_offset_compare_misses,
        multiply_offset_cleanup_compare_checks,
        multiply_offset_cleanup_compare_misses,
        first_host_miss,
        first_narrow_compare_miss,
        hclz_sites: hclz_sites
            .into_iter()
            .map(|(site, count)| Q945SupportSiteCount {
                site,
                checks: count.checks,
                zero_entry_checks: count.zero_entry_checks,
                restoration_checks: count.restoration_checks,
            })
            .collect(),
        carry_sites: carry_sites
            .into_iter()
            .map(|(site, count)| Q945SupportSiteCount {
                site,
                checks: count.checks,
                zero_entry_checks: count.zero_entry_checks,
                restoration_checks: count.restoration_checks,
            })
            .collect(),
        earliest_terminal_row,
        latest_terminal_row,
        support_clean,
    }
}

/// Replay every factor selected by the exact circuit Fiat-Shamir stream through
/// the ideal PZ transition and the explicit/affine terminal-counter differential.
/// This is deliberately separate from profiling: width and CLZ misses are
/// retained while the complete census runs, then any miss rejects the proof.
#[doc(hidden)]
pub fn q949_supported_trace_check(builder: &crate::point_add::B) -> Q949SupportedTraceReport {
    use inversion::shrunken_pz_schedule::{Q949TraceDirection, Q949WidthPhase};
    use std::collections::btree_map::Entry;

    assert_eq!(
        std::env::var("LOWQ_Q949_AFFINE_COUNTER").ok().as_deref(),
        Some("1"),
        "Q949 supported-trace proof requires the affine route"
    );
    let diagnostics_path = std::env::var(Q949_CENSUS_DIAGNOSTICS_OUT_ENV).unwrap_or_else(|_| {
        panic!("{Q949_CENSUS_DIAGNOSTICS_OUT_ENV} is required for the Q949 proof census")
    });
    assert!(
        !diagnostics_path.is_empty(),
        "{Q949_CENSUS_DIAGNOSTICS_OUT_ENV} cannot be empty"
    );
    let identity = q949_proof_identity(builder);
    let hasher = q949_builder_fiat_hash(builder);
    let requested_draws = env_usize("TRAILMIX_Q949_PROOF_SHOTS", TRAILMIX_NUM_TESTS);
    assert_eq!(
        requested_draws, Q949_PROOF_DRAWS,
        "Q949 production proof requires the full {Q949_PROOF_DRAWS}-draw census"
    );
    let curve = secp256k1();
    let mut xof = hasher.finalize_xof();
    let mut accepted_shots = 0usize;
    let mut rejected_draws = 0usize;
    let mut factors_checked = 0usize;
    let mut forward_rows_checked = 0usize;
    let mut backward_rows_checked = 0usize;
    let mut row_bounds_checked = 0usize;
    let mut entry_width_checks = 0usize;
    let mut transient_width_checks = 0usize;
    let mut post_swap_width_checks = 0usize;
    let mut boundary_width_checks = 0usize;
    let mut width_context_observations = 0usize;
    let mut entry_width_misses = 0usize;
    let mut transient_width_misses = 0usize;
    let mut post_swap_width_misses = 0usize;
    let mut boundary_width_misses = 0usize;
    let mut width_misses = 0usize;
    let mut clz_window_checks = 0usize;
    let mut clz_window_misses = 0usize;
    let mut clz_context_observations = 0usize;
    let mut first_width_miss = None;
    let mut first_clz_window_miss = None;
    let mut width_miss_buckets: BTreeMap<
        (Q949TraceDirection, Q949WidthPhase, usize, usize),
        Q949WidthMissBucket,
    > = BTreeMap::new();
    let mut clz_window_miss_buckets: BTreeMap<
        (Q949TraceDirection, usize, usize),
        Q949ClzWindowMissBucket,
    > = BTreeMap::new();
    let mut width_excess_histogram = BTreeMap::<usize, usize>::new();
    let mut clz_shortfall_histogram = BTreeMap::<usize, usize>::new();
    let mut joint_width_demands: BTreeMap<
        (Q949TraceDirection, Q949WidthPhase, usize),
        Q949JointWidthDemandAccumulator,
    > = BTreeMap::new();
    let mut clz_low_bounds: BTreeMap<
        (Q949TraceDirection, usize, usize),
        Q949ClzLowBoundAccumulator,
    > = BTreeMap::new();
    let mut terminal_full_ca_checks = 0usize;
    let mut reverse_row_380_relation_checks = 0usize;
    let mut reverse_row_380_active_checks = 0usize;
    let mut reverse_row_380_inactive_checks = 0usize;
    let mut reverse_row_380_relation_failures = 0usize;
    let mut earliest_terminal_row = usize::MAX;
    let mut latest_terminal_row = 0usize;
    let mut max_counter = 0usize;

    for draw in 0..requested_draws {
        let mut random = [[0u8; 32]; 2];
        xof.read(&mut random[0]);
        xof.read(&mut random[1]);
        let k1 = U256::from_le_bytes(random[0]);
        let k2 = U256::from_le_bytes(random[1]);
        let target = curve.mul(curve.gx, curve.gy, k1);
        let other = curve.mul(curve.gx, curve.gy, k2);
        if target.0 == other.0
            || (target.0.is_zero() && target.1.is_zero())
            || (other.0.is_zero() && other.1.is_zero())
        {
            rejected_draws += 1;
            continue;
        }
        let result = curve.add(target.0, target.1, other.0, other.1);
        accepted_shots += 1;
        for (factor_label, factor) in [
            ("dx", sub_mod_p(target.0, other.0, curve.modulus)),
            (
                "qx_minus_rx",
                sub_mod_p(other.0, result.0, curve.modulus),
            ),
        ] {
            let certificate = inversion::shrunken_pz_schedule::
                q949_affine_trace_certificate_u256(factor);
            factors_checked += 1;
            forward_rows_checked += certificate.rows_forward_checked;
            backward_rows_checked += certificate.rows_backward_checked;
            row_bounds_checked += certificate.row_bounds_checked;
            entry_width_checks += certificate.entry_width_checks;
            transient_width_checks += certificate.transient_width_checks;
            post_swap_width_checks += certificate.post_swap_width_checks;
            boundary_width_checks += certificate.boundary_width_checks;
            entry_width_misses += certificate.entry_width_misses;
            transient_width_misses += certificate.transient_width_misses;
            post_swap_width_misses += certificate.post_swap_width_misses;
            boundary_width_misses += certificate.boundary_width_misses;
            width_misses += certificate.width_misses;
            clz_window_checks += certificate.clz_window_checks;
            clz_window_misses += certificate.clz_window_misses;
            assert_eq!(
                certificate.width_observations.len(),
                Q949_WIDTH_CONTEXTS_PER_FACTOR,
                "Q949 factor width-context count drift"
            );
            width_context_observations += certificate.width_observations.len();
            for observation in &certificate.width_observations {
                let key = (observation.direction, observation.phase, observation.row);
                let witness = Q949WidthDemandWitness {
                    draw,
                    factor_label,
                    factor,
                    required_widths: observation.required_widths,
                };
                match joint_width_demands.entry(key) {
                    Entry::Vacant(entry) => {
                        entry.insert(Q949JointWidthDemandAccumulator::new(
                            *observation,
                            witness,
                        ));
                    }
                    Entry::Occupied(mut entry) => {
                        entry.get_mut().record(*observation, witness)
                    }
                }
            }
            assert_eq!(
                certificate.clz_window_observations.len() * 4,
                certificate.clz_window_checks,
                "Q949 factor CLZ-context count drift"
            );
            clz_context_observations += certificate.clz_window_observations.len();
            for observation in &certificate.clz_window_observations {
                for register in 0..4 {
                    let key = (observation.direction, observation.row, register);
                    let witness = Q949ClzLowWitness {
                        draw,
                        factor_label,
                        factor,
                        observed_width: observation.observed_widths[register],
                        observed_widths: observation.observed_widths,
                    };
                    match clz_low_bounds.entry(key) {
                        Entry::Vacant(entry) => {
                            entry.insert(Q949ClzLowBoundAccumulator::new(
                                *observation,
                                register,
                                witness,
                            ));
                        }
                        Entry::Occupied(mut entry) => {
                            entry.get_mut().record(*observation, register, witness)
                        }
                    }
                }
            }
            if first_width_miss.is_none() {
                first_width_miss = certificate.first_width_miss.map(|coordinate| {
                    Q949SupportedWidthMiss {
                        draw,
                        factor_label,
                        factor,
                        coordinate,
                    }
                });
            }
            if first_clz_window_miss.is_none() {
                first_clz_window_miss =
                    certificate
                        .first_clz_window_miss
                        .map(|coordinate| Q949SupportedClzWindowMiss {
                            draw,
                            factor_label,
                            factor,
                            coordinate,
                        });
            }
            assert_eq!(
                certificate.width_miss_coordinates.len(),
                certificate.width_misses,
                "Q949 factor width-miss coordinate count drift"
            );
            for coordinate in &certificate.width_miss_coordinates {
                let register = q949_width_register_index(coordinate.register);
                assert_eq!(
                    coordinate.observed_width,
                    coordinate.observed_widths[register]
                );
                assert_eq!(
                    coordinate.available_width,
                    coordinate.available_widths[register]
                );
                let excess = coordinate.observed_width - coordinate.available_width;
                assert!(excess > 0);
                *width_excess_histogram.entry(excess).or_default() += 1;
                let key = (
                    coordinate.direction,
                    coordinate.phase,
                    coordinate.row,
                    register,
                );
                match width_miss_buckets.entry(key) {
                    Entry::Vacant(entry) => {
                        entry.insert(Q949WidthMissBucket::new(*coordinate));
                    }
                    Entry::Occupied(mut entry) => entry.get_mut().record(*coordinate),
                }
            }
            assert_eq!(
                certificate.clz_window_miss_coordinates.len(),
                certificate.clz_window_misses,
                "Q949 factor CLZ-miss coordinate count drift"
            );
            for coordinate in &certificate.clz_window_miss_coordinates {
                let shortfall = coordinate.low + 1 - coordinate.observed_width;
                assert!(shortfall > 0);
                *clz_shortfall_histogram.entry(shortfall).or_default() += 1;
                let key = (
                    coordinate.direction,
                    coordinate.row,
                    q949_clz_register_index(coordinate.register),
                );
                match clz_window_miss_buckets.entry(key) {
                    Entry::Vacant(entry) => {
                        entry.insert(Q949ClzWindowMissBucket::new(*coordinate));
                    }
                    Entry::Occupied(mut entry) => entry.get_mut().record(*coordinate),
                }
            }
            terminal_full_ca_checks += certificate.terminal_full_ca_checks;
            reverse_row_380_relation_checks += certificate.reverse_row_380_relation_checks;
            reverse_row_380_active_checks += certificate.reverse_row_380_active_checks;
            reverse_row_380_inactive_checks += certificate.reverse_row_380_inactive_checks;
            reverse_row_380_relation_failures += certificate.reverse_row_380_relation_failures;
            earliest_terminal_row = earliest_terminal_row.min(certificate.first_terminal_row);
            latest_terminal_row = latest_terminal_row.max(certificate.first_terminal_row);
            max_counter = max_counter.max(certificate.max_counter);
        }
    }
    let report = Q949SupportedTraceReport {
        requested_draws,
        accepted_shots,
        rejected_draws,
        factors_checked,
        forward_rows_checked,
        backward_rows_checked,
        row_bounds_checked,
        entry_width_checks,
        transient_width_checks,
        post_swap_width_checks,
        boundary_width_checks,
        width_context_observations,
        entry_width_misses,
        transient_width_misses,
        post_swap_width_misses,
        boundary_width_misses,
        width_misses,
        clz_window_checks,
        clz_window_misses,
        clz_context_observations,
        first_width_miss,
        first_clz_window_miss,
        width_miss_buckets: width_miss_buckets.into_values().collect(),
        clz_window_miss_buckets: clz_window_miss_buckets.into_values().collect(),
        width_excess_histogram: width_excess_histogram.into_iter().collect(),
        clz_shortfall_histogram: clz_shortfall_histogram.into_iter().collect(),
        joint_width_demands: joint_width_demands
            .into_values()
            .map(Q949JointWidthDemandAccumulator::finish)
            .collect(),
        clz_low_bounds: clz_low_bounds
            .into_values()
            .map(Q949ClzLowBoundAccumulator::finish)
            .collect(),
        terminal_full_ca_checks,
        reverse_row_380_relation_checks,
        reverse_row_380_active_checks,
        reverse_row_380_inactive_checks,
        reverse_row_380_relation_failures,
        earliest_terminal_row,
        latest_terminal_row,
        max_counter,
    };

    q949_emit_census_diagnostics(&diagnostics_path, &identity, &report);
    assert_eq!(
        report.accepted_shots, Q949_PROOF_DRAWS,
        "Q949 census did not admit every draw; rejected={}",
        report.rejected_draws
    );
    assert_eq!(report.factors_checked, Q949_PROOF_FACTORS);
    assert_eq!(
        report.reverse_row_380_relation_checks,
        report.factors_checked,
        "Q949 reverse row-380 relation coverage drift"
    );
    assert_eq!(
        report.reverse_row_380_active_checks + report.reverse_row_380_inactive_checks,
        report.reverse_row_380_relation_checks,
        "Q949 reverse row-380 active/off partition drift"
    );
    assert!(
        report.reverse_row_380_relation_failures <= report.reverse_row_380_relation_checks,
        "Q949 reverse row-380 relation failure count drift"
    );
    assert_eq!(
        report.forward_rows_checked,
        Q949_PROOF_FACTORS * inversion::shrunken_pz_schedule::SHRUNKEN_PZ_NSTEPS
    );
    assert_eq!(report.forward_rows_checked, report.backward_rows_checked);
    assert_eq!(report.forward_rows_checked, report.row_bounds_checked);
    assert_eq!(
        report.entry_width_checks,
        Q949_PROOF_FACTORS * inversion::shrunken_pz_schedule::SHRUNKEN_PZ_NSTEPS * 5
    );
    assert_eq!(report.entry_width_checks, report.transient_width_checks);
    assert_eq!(report.entry_width_checks, report.post_swap_width_checks);
    assert_eq!(
        report.boundary_width_checks,
        Q949_PROOF_FACTORS
            * 2
            * (inversion::shrunken_pz_schedule::SHRUNKEN_PZ_NSTEPS - 1)
            * 5
    );
    assert_eq!(
        report
            .width_excess_histogram
            .iter()
            .map(|(_, count)| count)
            .sum::<usize>(),
        report.width_misses,
        "Q949 exact width-excess histogram does not cover every miss"
    );
    assert_eq!(
        report
            .clz_shortfall_histogram
            .iter()
            .map(|(_, count)| count)
            .sum::<usize>(),
        report.clz_window_misses,
        "Q949 exact CLZ-shortfall histogram does not cover every miss"
    );
    let diagnostic_rejection_report = std::env::var("Q949_DIAGNOSTIC_REJECTION_REPORT")
        .ok()
        .as_deref()
        == Some("1");
    if !diagnostic_rejection_report {
        assert_eq!(
            report.width_misses, 0,
            "Q949 support width census failed after all factors; entry={} transient={} post_swap={} boundary={} first={:?}",
            report.entry_width_misses,
            report.transient_width_misses,
            report.post_swap_width_misses,
            report.boundary_width_misses,
            report.first_width_miss,
        );
        assert_eq!(
            report.clz_window_misses, 0,
            "Q949 CLZ-window census failed after all factors; first={:?}",
            report.first_clz_window_miss,
        );
    }
    assert!(report.earliest_terminal_row <= report.latest_terminal_row);
    assert!(
        report.latest_terminal_row < inversion::shrunken_pz_schedule::SHRUNKEN_PZ_NSTEPS
    );
    assert_eq!(
        report.max_counter,
        inversion::shrunken_pz_schedule::SHRUNKEN_PZ_NSTEPS
            - report.earliest_terminal_row,
        "Q949 terminal count/earliest-row relation drift"
    );
    assert!(report.max_counter < 256, "Q949 affine terminal counter overflow");
    report
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
    let default_threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    let threads = env_usize("TRAILMIX_TAIL_NONCE_THREADS", default_threads)
        .max(1)
        .min(limit.max(1));

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
                                "TRAILMIX_TAIL_SEARCH nonce={} miss_factors={} repair_entries={} first_miss={:?}",
                                nonce, report.miss_factors, report.repair_entries, report.first_miss
                            );
                        }
                        let better = best.as_ref().map_or(true, |(_, b)| {
                            (report.miss_factors, report.repair_entries)
                                < (b.miss_factors, b.repair_entries)
                        });
                        if better {
                            best = Some((nonce, report.clone()));
                        }
                        if report.miss_factors == 0 {
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
                (report.miss_factors, report.repair_entries, nonce)
                    < (b.miss_factors, b.repair_entries, *best_nonce)
            });
            if better {
                best = Some((nonce, report));
            }
        }
    }
    if let Some((nonce, report)) = best {
        eprintln!(
            "TRAILMIX_TAIL_SEARCH_BEST nonce={} accepted={} miss_factors={} repair_entries={} first_miss={:?} searched={} threads={}",
            nonce,
            report.accepted_shots,
            report.miss_factors,
            report.repair_entries,
            report.first_miss,
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

    // Reserve the source-baked Q851 stream exactly. The fixed-sign event route
    // emits 838,845,965 records; avoiding the older 1,355,086,804-record
    // reservation is necessary for the 64 GiB submission-server envelope.
    // Vec can still grow for an explicitly overridden experimental route.
    let mut circ = circuit::Circuit::new_with_ops_capacity(838_845_965);
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
