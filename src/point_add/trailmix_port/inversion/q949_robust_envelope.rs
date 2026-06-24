use std::collections::BTreeSet;
use std::sync::OnceLock;

mod embedded_data {
    include!("q949_robust_envelope_data.rs");
}

pub const Q949_ROBUST_ENVELOPE_SCHEMA: &str = "q948-peak-safe-robust-row-envelope-v2";
pub const Q949_ROBUST_ENVELOPE_SHA256: &str =
    "ad72e9ef9d0be9b91f22fdc88fe7437fdedb3426ac032db63e6c28a6f6ee2e8e";
pub const Q949_ROBUST_PROJECTION_SHA256: &str =
    "b9f883b8e0ef831437c2ab2d1abf8378cd67fdab5a6b6d50b13e7fc8cc348465";
pub const Q949_ROBUST_SELECTION_RESULT_SHA256: &str =
    "c0f576c793e3c9dbcad53496cd7bd2b107a53ca8585e7a89a799567f49f9a697";
pub const Q949_ROBUST_SELECTION_JOB_ID: u64 = 71_581;
pub const Q949_ROBUST_ENVELOPE_GENERATION_JOB_ID: u64 = 71_581;
pub const Q949_ROBUST_ARTIFACT_BYTES: usize = 10_545_370;
pub const Q949_ROBUST_ENVELOPE_FRESH_VALIDITY_CLAIMED: bool = false;
pub const Q949_ROBUST_ROWS: usize = 530;
pub const Q949_ROBUST_TARGET_SUM: usize = 683;
pub const Q949_ROBUST_PEAK_SAFE_SUM: usize = 681;
pub const Q949_ROBUST_TRAINING_STREAMS: usize = 14;

pub const Q949_ROBUST_TRAINING_SOURCE_IDS: [&str; Q949_ROBUST_TRAINING_STREAMS] =
    embedded_data::TRAINING_SOURCE_IDS;
pub const Q949_ROBUST_TRAINING_SCHEDULE_IDS: [&str; Q949_ROBUST_TRAINING_STREAMS] =
    embedded_data::TRAINING_SCHEDULE_IDS;
pub const Q949_ROBUST_TRAINING_ROUTE_IDS: [&str; Q949_ROBUST_TRAINING_STREAMS] =
    embedded_data::TRAINING_ROUTE_IDS;
pub const Q949_ROBUST_TRAINING_OP_STREAM_IDS: [&str; Q949_ROBUST_TRAINING_STREAMS] =
    embedded_data::TRAINING_OP_STREAM_IDS;
pub const Q949_ROBUST_TRAINING_NONCES: [u64; Q949_ROBUST_TRAINING_STREAMS] =
    embedded_data::TRAINING_NONCES;
pub const Q949_ROBUST_TRAINING_JOB_IDS: [u64; Q949_ROBUST_TRAINING_STREAMS] =
    embedded_data::TRAINING_JOB_IDS;
pub const Q949_ROBUST_TRAINING_OP_COUNTS: [usize; Q949_ROBUST_TRAINING_STREAMS] =
    embedded_data::TRAINING_OP_COUNTS;
pub const Q949_ROBUST_TRAINING_DIAGNOSTICS_SHA256: [&str; Q949_ROBUST_TRAINING_STREAMS] =
    embedded_data::TRAINING_DIAGNOSTICS_SHA256;
pub const Q949_ROBUST_TRAINING_SOURCE_COMMITS: [&str; Q949_ROBUST_TRAINING_STREAMS] =
    embedded_data::TRAINING_SOURCE_COMMITS;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Q949RobustEnvelopeTrainingReport {
    pub artifact_bytes: usize,
    pub rows_checked: usize,
    pub width_component_witnesses_checked: usize,
    pub clz_contexts_checked: usize,
    pub clz_limiting_witnesses_checked: usize,
    pub minimum_pair_symmetric_slack: usize,
    pub maximum_pair_symmetric_sum: usize,
    pub training_streams_seen: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Q949RobustPhaseRequirement {
    ForwardBoundary,
    ForwardEntry,
    ForwardPostSwap,
    ForwardTransient,
    ReverseBoundary,
}

impl Q949RobustPhaseRequirement {
    fn index(self) -> usize {
        match self {
            Self::ForwardBoundary => 0,
            Self::ForwardEntry => 1,
            Self::ForwardPostSwap => 2,
            Self::ForwardTransient => 3,
            Self::ReverseBoundary => 4,
        }
    }
}

#[derive(Debug)]
struct Q949RobustEnvelope {
    requirements: [[usize; 5]; Q949_ROBUST_ROWS],
    phase_requirements: [[[usize; 5]; 5]; Q949_ROBUST_ROWS],
    widths: [[usize; 5]; Q949_ROBUST_ROWS],
    lows: [[usize; 5]; Q949_ROBUST_ROWS],
    clz_safe_low_upper_bounds: [[Option<usize>; 4]; Q949_ROBUST_ROWS],
    clz_bound_present: [[bool; 4]; Q949_ROBUST_ROWS],
    report: Q949RobustEnvelopeTrainingReport,
}

static ROBUST_ENVELOPE: OnceLock<Q949RobustEnvelope> = OnceLock::new();

fn assert_lower_hex(label: &str, value: &str, bytes: usize) {
    assert_eq!(value.len(), 2 * bytes, "{label} length drift");
    assert!(
        value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte)),
        "{label} is not lowercase hexadecimal"
    );
}

fn assert_compact_training_provenance() {
    assert_eq!(embedded_data::SELECTION_JOB_ID, Q949_ROBUST_SELECTION_JOB_ID);
    assert_eq!(
        embedded_data::SELECTION_RESULT_SHA256,
        Q949_ROBUST_SELECTION_RESULT_SHA256
    );
    assert_eq!(
        embedded_data::SELECTION_MAXIMUM_CARDINALITY,
        Q949_ROBUST_TRAINING_STREAMS
    );
    assert_eq!(embedded_data::PEAK_SAFE_PAIR_SYMMETRIC_CAP, 681);
    assert_eq!(
        embedded_data::TRAINING_STREAM_MASK.count_ones() as usize,
        Q949_ROBUST_TRAINING_STREAMS
    );
    let mut operation_streams = BTreeSet::new();
    for index in 0..Q949_ROBUST_TRAINING_STREAMS {
        for (label, value, bytes) in [
            ("training source ID", Q949_ROBUST_TRAINING_SOURCE_IDS[index], 32),
            (
                "training schedule ID",
                Q949_ROBUST_TRAINING_SCHEDULE_IDS[index],
                32,
            ),
            ("training route ID", Q949_ROBUST_TRAINING_ROUTE_IDS[index], 32),
            (
                "training operation-stream ID",
                Q949_ROBUST_TRAINING_OP_STREAM_IDS[index],
                32,
            ),
            (
                "training diagnostics SHA256",
                Q949_ROBUST_TRAINING_DIAGNOSTICS_SHA256[index],
                32,
            ),
            (
                "training source commit",
                Q949_ROBUST_TRAINING_SOURCE_COMMITS[index],
                20,
            ),
        ] {
            assert_lower_hex(label, value, bytes);
        }
        assert!(operation_streams.insert(Q949_ROBUST_TRAINING_OP_STREAM_IDS[index]));
        assert_eq!(embedded_data::TRAINING_NONCE_BITS[index], 48);
        assert!(Q949_ROBUST_TRAINING_NONCES[index] < (1u64 << 48));
        assert!(Q949_ROBUST_TRAINING_JOB_IDS[index] > 0);
        assert!(Q949_ROBUST_TRAINING_OP_COUNTS[index] > 0);
    }
    assert!(!Q949_ROBUST_ENVELOPE_FRESH_VALIDITY_CLAIMED);
}

fn parse_envelope() -> Q949RobustEnvelope {
    assert_eq!(embedded_data::ENVELOPE_SCHEMA, Q949_ROBUST_ENVELOPE_SCHEMA);
    assert_eq!(embedded_data::ENVELOPE_SHA256, Q949_ROBUST_ENVELOPE_SHA256);
    assert_eq!(embedded_data::PROJECTION_SHA256, Q949_ROBUST_PROJECTION_SHA256);
    assert_lower_hex(
        "Q949 robust projection SHA256",
        Q949_ROBUST_PROJECTION_SHA256,
        32,
    );
    assert_eq!(embedded_data::ARTIFACT_BYTES, Q949_ROBUST_ARTIFACT_BYTES);
    assert_compact_training_provenance();

    let requirements = embedded_data::ROW_REQUIREMENTS.map(|row| row.map(usize::from));
    let phase_requirements = embedded_data::PHASE_REQUIREMENTS
        .map(|row| row.map(|phase| phase.map(usize::from)));
    let mut widths = [[0usize; 5]; Q949_ROBUST_ROWS];
    let mut lows = [[0usize; 5]; Q949_ROBUST_ROWS];
    let mut clz_safe_low_upper_bounds = [[None; 4]; Q949_ROBUST_ROWS];
    let mut clz_bound_present = [[false; 4]; Q949_ROBUST_ROWS];
    let mut minimum_pair_symmetric_slack = usize::MAX;
    let mut maximum_pair_symmetric_sum = 0usize;
    let mut clz_contexts_checked = 0usize;
    let mut unconstrained_clz_contexts = 0usize;

    for row in 0..Q949_ROBUST_ROWS {
        let required = requirements[row];
        let ab = required[0].max(required[1]);
        let cacb = required[2].max(required[3]);
        let symmetric = [ab, ab, cacb, cacb, required[4]];
        assert!(symmetric.into_iter().all(|width| width > 0));
        assert!(symmetric[4] <= 99, "Q949 robust row {row} exceeds Q_CAP");
        let sum = symmetric.iter().sum::<usize>();
        assert!(
            sum <= Q949_ROBUST_PEAK_SAFE_SUM,
            "Q949 robust row {row} exceeds peak-safe capacity"
        );
        let slack = Q949_ROBUST_TARGET_SUM - sum;
        minimum_pair_symmetric_slack = minimum_pair_symmetric_slack.min(slack);
        maximum_pair_symmetric_sum = maximum_pair_symmetric_sum.max(sum);
        widths[row] = symmetric;

        for (phase_index, phase) in phase_requirements[row].iter().enumerate() {
            let absent = phase.iter().all(|&width| width == 0);
            assert!(
                !absent
                    || (row == 0 && phase_index == 0)
                    || (row == Q949_ROBUST_ROWS - 1 && phase_index == 4),
                "unexpected missing robust phase requirement"
            );
            assert!(
                phase
                    .iter()
                    .zip(symmetric)
                    .all(|(&phase_width, row_width)| phase_width == 0 || phase_width <= row_width),
                "robust phase requirement exceeds row allocation"
            );
        }

        for register in 0..4 {
            match embedded_data::CLZ_SAFE_LOW_UPPER_BOUNDS[row][register] {
                safe_low if safe_low >= 0 => {
                    let safe_low = safe_low as usize;
                    assert!(safe_low < symmetric[register]);
                    clz_bound_present[row][register] = true;
                    clz_safe_low_upper_bounds[row][register] = Some(safe_low);
                    lows[row][register] = safe_low;
                    clz_contexts_checked += 1;
                }
                -1 => {
                    clz_bound_present[row][register] = true;
                    clz_contexts_checked += 1;
                    unconstrained_clz_contexts += 1;
                }
                -2 => {}
                sentinel => panic!("unknown Q949 robust CLZ sentinel: {sentinel}"),
            }
        }
    }
    assert_eq!(minimum_pair_symmetric_slack, 2);
    assert_eq!(maximum_pair_symmetric_sum, 681);
    assert_eq!(
        minimum_pair_symmetric_slack,
        embedded_data::MINIMUM_PAIR_SYMMETRIC_SLACK
    );
    assert_eq!(
        maximum_pair_symmetric_sum,
        embedded_data::MAXIMUM_PAIR_SYMMETRIC_SUM
    );
    assert_eq!(clz_contexts_checked, embedded_data::CLZ_CONTEXTS);
    assert_eq!(unconstrained_clz_contexts, 2);
    assert_eq!(
        unconstrained_clz_contexts,
        embedded_data::UNCONSTRAINED_CLZ_CONTEXTS
    );

    Q949RobustEnvelope {
        requirements,
        phase_requirements,
        widths,
        lows,
        clz_safe_low_upper_bounds,
        clz_bound_present,
        report: Q949RobustEnvelopeTrainingReport {
            artifact_bytes: embedded_data::ARTIFACT_BYTES,
            rows_checked: Q949_ROBUST_ROWS,
            width_component_witnesses_checked: embedded_data::WIDTH_COMPONENT_WITNESSES,
            clz_contexts_checked,
            clz_limiting_witnesses_checked: embedded_data::CLZ_LIMITING_WITNESSES,
            minimum_pair_symmetric_slack,
            maximum_pair_symmetric_sum,
            training_streams_seen: embedded_data::TRAINING_STREAM_MASK.count_ones() as usize,
        },
    }
}

fn envelope() -> &'static Q949RobustEnvelope {
    ROBUST_ENVELOPE.get_or_init(parse_envelope)
}

#[must_use]
pub fn q949_robust_envelope_sha256() -> String {
    let _ = envelope();
    Q949_ROBUST_ENVELOPE_SHA256.to_owned()
}

#[must_use]
pub fn q949_robust_row_requirements(row: usize) -> [usize; 5] {
    envelope().requirements[row.min(Q949_ROBUST_ROWS - 1)]
}

pub fn q949_robust_phase_requirements(
    row: usize,
    phase: Q949RobustPhaseRequirement,
) -> Option<[usize; 5]> {
    let requirements = envelope().phase_requirements[row.min(Q949_ROBUST_ROWS - 1)][phase.index()];
    requirements.iter().any(|&width| width != 0).then_some(requirements)
}

#[must_use]
pub fn q949_robust_pair_symmetric_widths(row: usize) -> [usize; 5] {
    envelope().widths[row.min(Q949_ROBUST_ROWS - 1)]
}

#[must_use]
pub fn q949_robust_clz_lows(row: usize) -> [usize; 5] {
    envelope().lows[row.min(Q949_ROBUST_ROWS - 1)]
}

#[must_use]
pub fn q949_robust_clz_safe_low_upper_bounds(row: usize) -> [Option<usize>; 4] {
    envelope().clz_safe_low_upper_bounds[row.min(Q949_ROBUST_ROWS - 1)]
}

#[must_use]
pub fn q949_robust_clz_bound_present(row: usize) -> [bool; 4] {
    envelope().clz_bound_present[row.min(Q949_ROBUST_ROWS - 1)]
}

#[must_use]
pub fn q949_robust_envelope_training_check() -> Q949RobustEnvelopeTrainingReport {
    envelope().report
}
