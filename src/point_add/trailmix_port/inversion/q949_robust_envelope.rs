use std::collections::BTreeSet;
use std::sync::OnceLock;

// Keep the submitted point_add tree independent of research-only Cargo
// dependencies. This compact one-shot implementation hashes only the sealed
// provenance payload below; it does not affect the emitted circuit.
struct Sha256 {
    bytes: Vec<u8>,
}

impl Sha256 {
    fn new() -> Self {
        Self { bytes: Vec::new() }
    }

    fn update(&mut self, bytes: impl AsRef<[u8]>) {
        self.bytes.extend_from_slice(bytes.as_ref());
    }

    fn finalize(self) -> [u8; 32] {
        sha256(&self.bytes)
    }
}

fn sha256(input: &[u8]) -> [u8; 32] {
    const K: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1,
        0x923f82a4, 0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
        0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786,
        0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147,
        0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
        0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
        0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a,
        0x5b9cca4f, 0x682e6ff3, 0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
        0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
    ];

    let bit_len = u64::try_from(input.len())
        .expect("Q949 robust projection input exceeds u64")
        .checked_mul(8)
        .expect("Q949 robust projection bit length overflow");
    let mut padded = Vec::with_capacity(input.len() + 72);
    padded.extend_from_slice(input);
    padded.push(0x80);
    while padded.len() % 64 != 56 {
        padded.push(0);
    }
    padded.extend_from_slice(&bit_len.to_be_bytes());

    let mut state = [
        0x6a09e667u32,
        0xbb67ae85,
        0x3c6ef372,
        0xa54ff53a,
        0x510e527f,
        0x9b05688c,
        0x1f83d9ab,
        0x5be0cd19,
    ];
    for chunk in padded.chunks_exact(64) {
        let mut words = [0u32; 64];
        for (index, word) in words[..16].iter_mut().enumerate() {
            let offset = 4 * index;
            *word = u32::from_be_bytes(chunk[offset..offset + 4].try_into().unwrap());
        }
        for index in 16..64 {
            let s0 = words[index - 15].rotate_right(7)
                ^ words[index - 15].rotate_right(18)
                ^ (words[index - 15] >> 3);
            let s1 = words[index - 2].rotate_right(17)
                ^ words[index - 2].rotate_right(19)
                ^ (words[index - 2] >> 10);
            words[index] = words[index - 16]
                .wrapping_add(s0)
                .wrapping_add(words[index - 7])
                .wrapping_add(s1);
        }

        let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut h] = state;
        for index in 0..64 {
            let sum1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let choose = (e & f) ^ (!e & g);
            let temp1 = h
                .wrapping_add(sum1)
                .wrapping_add(choose)
                .wrapping_add(K[index])
                .wrapping_add(words[index]);
            let sum0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let majority = (a & b) ^ (a & c) ^ (b & c);
            let temp2 = sum0.wrapping_add(majority);

            h = g;
            g = f;
            f = e;
            e = d.wrapping_add(temp1);
            d = c;
            c = b;
            b = a;
            a = temp1.wrapping_add(temp2);
        }

        for (slot, value) in state.iter_mut().zip([a, b, c, d, e, f, g, h]) {
            *slot = slot.wrapping_add(value);
        }
    }

    let mut output = [0u8; 32];
    for (chunk, value) in output.chunks_exact_mut(4).zip(state) {
        chunk.copy_from_slice(&value.to_be_bytes());
    }
    output
}

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

fn hex(bytes: &[u8]) -> String {
    use std::fmt::Write as _;

    let mut output = String::with_capacity(2 * bytes.len());
    for byte in bytes {
        write!(&mut output, "{byte:02x}").expect("write Q949 robust projection SHA256");
    }
    output
}

fn assert_lower_hex(label: &str, value: &str, bytes: usize) {
    assert_eq!(value.len(), 2 * bytes, "{label} length drift");
    assert!(
        value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte)),
        "{label} is not lowercase hexadecimal"
    );
}

fn compact_projection_sha256() -> String {
    let mut hasher = Sha256::new();
    hasher.update(embedded_data::PROJECTION_DOMAIN);
    hasher.update(Q949_ROBUST_ENVELOPE_SHA256.as_bytes());
    hasher.update(Q949_ROBUST_SELECTION_RESULT_SHA256.as_bytes());
    hasher.update(embedded_data::SELECTION_JOB_ID.to_le_bytes());
    hasher.update((embedded_data::SELECTION_MAXIMUM_CARDINALITY as u64).to_le_bytes());
    hasher.update((embedded_data::PEAK_SAFE_PAIR_SYMMETRIC_CAP as u64).to_le_bytes());

    for index in 0..Q949_ROBUST_TRAINING_STREAMS {
        for value in [
            Q949_ROBUST_TRAINING_SOURCE_IDS[index],
            Q949_ROBUST_TRAINING_SCHEDULE_IDS[index],
            Q949_ROBUST_TRAINING_ROUTE_IDS[index],
            Q949_ROBUST_TRAINING_OP_STREAM_IDS[index],
        ] {
            hasher.update(value.as_bytes());
        }
        hasher.update(Q949_ROBUST_TRAINING_NONCES[index].to_le_bytes());
        hasher.update((embedded_data::TRAINING_NONCE_BITS[index] as u64).to_le_bytes());
        hasher.update((Q949_ROBUST_TRAINING_OP_COUNTS[index] as u64).to_le_bytes());
        hasher.update(Q949_ROBUST_TRAINING_DIAGNOSTICS_SHA256[index].as_bytes());
        hasher.update(Q949_ROBUST_TRAINING_JOB_IDS[index].to_le_bytes());
        hasher.update(Q949_ROBUST_TRAINING_SOURCE_COMMITS[index].as_bytes());
        hasher.update([1u8, 1, 0]);
    }
    for row in embedded_data::ROW_REQUIREMENTS {
        for value in row {
            hasher.update(value.to_le_bytes());
        }
    }
    for row in embedded_data::PHASE_REQUIREMENTS {
        for phase in row {
            for value in phase {
                hasher.update(value.to_le_bytes());
            }
        }
    }
    for row in embedded_data::CLZ_SAFE_LOW_UPPER_BOUNDS {
        for value in row {
            hasher.update(value.to_le_bytes());
        }
    }
    hex(&hasher.finalize())
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
    assert_eq!(compact_projection_sha256(), Q949_ROBUST_PROJECTION_SHA256);
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
