//! Coherent lowering of Luo et al.'s 835-qubit fixed-schedule EEA.
//!
//! The paper implementation uses measurement-based unary uncomputation and
//! exposes a resource-only placeholder for the inverse EEA.  This backend
//! instead embeds the fully decomposed X/CX/CCX step stream, emits it forward,
//! and emits the exact reversed stream for cleanup.  The surrounding divider
//! keeps the source live through multiplication and uses HMR only for product
//! transport, matching the executable register-shared lifecycle.

use crate::circuit::OperationType;
use crate::point_add::trailmix_port::circuit::{Circuit, QReg};
use std::io::Cursor;

const FIELD_WIDTH: usize = 257;
const VALUE_WIDTH: usize = 256;
const WORK_WIDTH: usize = 259;
const LT_WIDTH: usize = 8;
const LQ_WIDTH: usize = 9;
const SHIFT_WIDTH: usize = 9;
const LRP_WIDTH: usize = 8;
const AUX_WIDTH: usize = 12;
const LQ_ZERO_ENCODING: usize = (1 << LQ_WIDTH) - 1;
const LS_ZERO_ENCODING: usize = 258;
const LRP_ZERO_ENCODING: usize = (1 << LRP_WIDTH) - 1;
const CORE_WIDTH: usize = 568;
const DIRTY_REFERENCE_WIDTH: usize = 10;
const LOCAL_WIDTH: usize = CORE_WIDTH + DIRTY_REFERENCE_WIDTH;
const SCHEDULE_STEPS: usize = 1_616;
const STREAM_X_PER_TRAVERSAL: usize = 25_190_680;
const STREAM_CX_PER_TRAVERSAL: usize = 23_020_144;
// Includes the two emitted CCX gates for every clean-C3X MBU marker.
const STREAM_CCX_PER_TRAVERSAL: usize = 45_453_265;
const STREAM_HMR_PER_TRAVERSAL: usize = 5_378_204;
const STREAM_CZ_PER_TRAVERSAL: usize = 5_378_204;

const fn half_plus_one_le() -> [u8; 33] {
    let mut bytes = [0xff; 33];
    bytes[0] = 0x18;
    bytes[1] = 0xfe;
    bytes[3] = 0x7f;
    bytes[31] = 0x7f;
    bytes[32] = 0;
    bytes
}

const HALF_PLUS_ONE_LE: [u8; 33] = half_plus_one_le();

const STREAM_CHUNKS: [&[u8]; 36] = [
    include_bytes!("paper2607_exactwidth_data/chunk-0001-0045.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-0046-0090.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-0091-0135.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-0136-0180.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-0181-0225.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-0226-0270.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-0271-0315.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-0316-0360.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-0361-0405.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-0406-0450.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-0451-0495.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-0496-0540.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-0541-0585.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-0586-0630.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-0631-0675.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-0676-0720.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-0721-0765.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-0766-0810.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-0811-0855.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-0856-0900.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-0901-0945.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-0946-0990.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-0991-1035.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-1036-1080.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-1081-1125.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-1126-1170.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-1171-1215.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-1216-1260.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-1261-1305.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-1306-1350.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-1351-1395.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-1396-1440.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-1441-1485.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-1486-1530.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-1531-1575.zst"),
    include_bytes!("paper2607_exactwidth_data/chunk-1576-1616.zst"),
];

pub fn enabled() -> bool {
    std::env::var("PAPER2607_COHERENT_EEA").ok().as_deref() == Some("1")
}

fn read_u32(bytes: &[u8], offset: usize) -> u32 {
    u32::from_le_bytes(bytes[offset..offset + 4].try_into().expect("u32 record"))
}

fn decode_chunk(compressed: &[u8]) -> Vec<u8> {
    let decoded = zstd::stream::decode_all(Cursor::new(compressed))
        .expect("decode paper2607 primitive stream");
    assert!(decoded.len() >= 24, "truncated paper2607 stream header");
    assert_eq!(&decoded[..8], b"P26EEA2\0");
    assert_eq!(read_u32(&decoded, 8), VALUE_WIDTH as u32);
    assert_eq!(read_u32(&decoded, 12), LOCAL_WIDTH as u32);
    assert_eq!((decoded.len() - 24) % 8, 0, "partial paper2607 record");
    decoded
}

fn emit_record(circ: &mut Circuit, local: &[&QReg], word: u64) {
    let kind = (word & 0xf) as u8;
    let arity = ((word >> 4) & 0xf) as usize;
    let q0 = ((word >> 8) & 0x3ff) as usize;
    let q1 = ((word >> 18) & 0x3ff) as usize;
    let q2 = ((word >> 28) & 0x3ff) as usize;
    let q3 = ((word >> 38) & 0x3ff) as usize;
    let q4 = ((word >> 48) & 0x3ff) as usize;
    assert!(q0 < local.len());
    match (kind, arity) {
        (1, 1) => circ.x(local[q0]),
        (2, 2) => {
            assert!(q1 < local.len());
            circ.cx(local[q0], local[q1]);
        }
        (3, 3) => {
            assert!(q1 < local.len() && q2 < local.len());
            circ.ccx(local[q0], local[q1], local[q2]);
        }
        (4, 1) => circ.z(local[q0]),
        (5, 2) => {
            assert!(q1 < local.len());
            circ.cz(local[q0], local[q1]);
        }
        (6, 2) => {
            assert!(q1 < local.len());
            circ.swap(local[q0], local[q1]);
        }
        (7, 5) => {
            assert!(q1 < local.len() && q2 < local.len());
            assert!(q3 < local.len() && q4 < local.len());
            circ.ccx(local[q0], local[q1], local[q4]);
            circ.ccx(local[q2], local[q4], local[q3]);
            circ.clear_and(local[q4], local[q0], local[q1]);
        }
        _ => panic!("invalid paper2607 primitive kind={kind} arity={arity}"),
    }
}

struct Core {
    phase1: QReg,
    phase2: QReg,
    iteration: QReg,
    sign: QReg,
    work1: Vec<QReg>,
    work2: Vec<QReg>,
    l_t: Vec<QReg>,
    l_q: Vec<QReg>,
    l_s: Vec<QReg>,
    l_rp: Vec<QReg>,
    aux: Vec<QReg>,
}

struct Terminal {
    iteration: QReg,
    work2: Vec<QReg>,
    l_s: Vec<QReg>,
}

struct CanonicalTopLoan {
    restored: bool,
    context: &'static str,
}

impl Drop for CanonicalTopLoan {
    fn drop(&mut self) {
        assert!(
            self.restored || std::thread::panicking(),
            "{} canonical top loan dropped without restore",
            self.context
        );
    }
}

/// Lend a canonical field register's known-zero extension lane to the EEA.
/// The replacement lane need not retain physical identity because the 257th
/// lane is internal, canonical zero state rather than ABI-visible data.
fn loan_canonical_top(
    circ: &mut Circuit,
    register: &mut Vec<QReg>,
    context: &'static str,
) -> CanonicalTopLoan {
    assert_eq!(
        register.len(),
        FIELD_WIDTH,
        "{context} canonical register width"
    );
    let live_before = circ.b.active_qubits;
    let top = register.pop().expect("canonical top lane");
    circ.zero_and_free(top);
    assert_eq!(register.len(), FIELD_WIDTH - 1);
    assert_eq!(
        circ.b.active_qubits + 1,
        live_before,
        "{context} canonical top loan must free one qubit"
    );
    circ.lowq_passenger_top_releases += 1;
    CanonicalTopLoan {
        restored: false,
        context,
    }
}

fn restore_canonical_top(circ: &mut Circuit, register: &mut Vec<QReg>, mut loan: CanonicalTopLoan) {
    assert_eq!(
        register.len(),
        FIELD_WIDTH - 1,
        "{} shortened canonical register width",
        loan.context
    );
    let live_before = circ.b.active_qubits;
    register.push(circ.alloc_qreg(&format!("{}.restored", loan.context)));
    assert_eq!(register.len(), FIELD_WIDTH);
    assert_eq!(
        circ.b.active_qubits,
        live_before + 1,
        "{} canonical top restore must allocate one clean qubit",
        loan.context
    );
    assert!(
        circ.lowq_passenger_top_releases > 0,
        "passenger top loan state underflow"
    );
    circ.lowq_passenger_top_releases -= 1;
    loan.restored = true;
}

fn free_clean(circ: &mut Circuit, register: Vec<QReg>) {
    for lane in register {
        circ.zero_and_free(lane);
    }
}

fn toggle_constant(circ: &mut Circuit, register: &[QReg], value: usize) {
    for (index, lane) in register.iter().enumerate() {
        if (value >> index) & 1 != 0 {
            circ.x(lane);
        }
    }
}

fn toggle_initial_work1(circ: &mut Circuit, work1: &[QReg]) {
    use crate::point_add::trailmix_port::mod_arith::SECP256K1_P_LE;

    assert_eq!(work1.len(), WORK_WIDTH);
    circ.x(&work1[0]);
    for bit in 0..VALUE_WIDTH {
        if (SECP256K1_P_LE[bit / 8] >> (bit % 8)) & 1 != 0 {
            circ.x(&work1[WORK_WIDTH - 1 - bit]);
        }
    }
}

fn toggle_terminal_work1(circ: &mut Circuit, work1: &[QReg]) {
    use crate::point_add::trailmix_port::mod_arith::SECP256K1_P_LE;

    assert_eq!(work1.len(), WORK_WIDTH);
    for bit in 0..VALUE_WIDTH {
        if (SECP256K1_P_LE[bit / 8] >> (bit % 8)) & 1 != 0 {
            circ.x(&work1[bit]);
        }
    }
    circ.x(&work1[WORK_WIDTH - 1]);
}

fn local_wires<'a>(core: &'a Core, passenger: &'a [QReg]) -> Vec<&'a QReg> {
    assert!(
        passenger.len() >= DIRTY_REFERENCE_WIDTH,
        "paper2607 dirty-passenger lender shortage"
    );
    let mut wires = Vec::with_capacity(LOCAL_WIDTH);
    wires.extend([&core.phase1, &core.phase2, &core.iteration, &core.sign]);
    wires.extend(core.work1.iter());
    wires.extend(core.work2.iter());
    wires.extend(core.l_t.iter());
    wires.extend(core.l_q.iter());
    wires.extend(core.l_s.iter());
    wires.extend(core.l_rp.iter());
    wires.extend(core.aux.iter());
    wires.extend(passenger.iter().take(DIRTY_REFERENCE_WIDTH));
    assert_eq!(wires.len() - DIRTY_REFERENCE_WIDTH, CORE_WIDTH);
    assert_eq!(wires.len(), LOCAL_WIDTH);
    wires
}

fn count_stub_enabled(circ: &Circuit) -> bool {
    circ.b.count_only
        && std::env::var_os("POINT_ADD_HASH_OPS_LEN").is_none()
        && std::env::var("PAPER2607_COUNT_STUB").ok().as_deref() == Some("1")
}

fn emit_count_stub(circ: &mut Circuit) {
    circ.b
        .add_counted_kind(OperationType::X, STREAM_X_PER_TRAVERSAL);
    circ.b
        .add_counted_kind(OperationType::CX, STREAM_CX_PER_TRAVERSAL);
    circ.b
        .add_counted_kind(OperationType::CCX, STREAM_CCX_PER_TRAVERSAL);
    circ.b
        .add_counted_kind(OperationType::Hmr, STREAM_HMR_PER_TRAVERSAL);
    circ.b
        .add_counted_kind(OperationType::CZ, STREAM_CZ_PER_TRAVERSAL);
}

fn emit_forward(circ: &mut Circuit, core: &Core, passenger: &[QReg]) {
    if count_stub_enabled(circ) {
        emit_count_stub(circ);
        return;
    }
    let wires = local_wires(core, passenger);
    let mut expected_start = 1_u32;
    for compressed in STREAM_CHUNKS {
        let decoded = decode_chunk(compressed);
        let start = read_u32(&decoded, 16);
        let end = read_u32(&decoded, 20);
        assert_eq!(start, expected_start, "paper2607 chunk gap");
        assert!(end >= start && end <= SCHEDULE_STEPS as u32);
        for record in decoded[24..].chunks_exact(8) {
            emit_record(
                circ,
                &wires,
                u64::from_le_bytes(record.try_into().expect("primitive record")),
            );
        }
        expected_start = end + 1;
    }
    assert_eq!(expected_start, SCHEDULE_STEPS as u32 + 1);
}

fn emit_reverse(circ: &mut Circuit, core: &Core, passenger: &[QReg]) {
    if count_stub_enabled(circ) {
        emit_count_stub(circ);
        return;
    }
    let wires = local_wires(core, passenger);
    let mut expected_end = SCHEDULE_STEPS as u32;
    for compressed in STREAM_CHUNKS.iter().rev() {
        let decoded = decode_chunk(compressed);
        let start = read_u32(&decoded, 16);
        let end = read_u32(&decoded, 20);
        assert_eq!(end, expected_end, "paper2607 reverse chunk gap");
        assert!(start >= 1 && start <= end);
        for record in decoded[24..].chunks_exact(8).rev() {
            emit_record(
                circ,
                &wires,
                u64::from_le_bytes(record.try_into().expect("primitive record")),
            );
        }
        expected_end = start - 1;
    }
    assert_eq!(expected_end, 0);
}

fn rotation_swaps(width: usize, shift: usize) -> Vec<(usize, usize)> {
    let shift = shift % width;
    if shift == 0 {
        return Vec::new();
    }
    let mut seen = vec![false; width];
    let mut swaps = Vec::with_capacity(width - 1);
    for start in 0..width {
        if seen[start] {
            continue;
        }
        let mut cycle = Vec::new();
        let mut lane = start;
        while !seen[lane] {
            seen[lane] = true;
            cycle.push(lane);
            lane = (lane + shift) % width;
        }
        for &other in cycle.iter().skip(1) {
            swaps.push((cycle[0], other));
        }
    }
    swaps
}

fn canonicalize_terminal_work2(circ: &mut Circuit, terminal: &Terminal) {
    // l_s stores (shift - 1) mod 259.  Apply the missing unit rotation
    // directly, then use the encoded bits for the remaining rotation.
    for (left, right) in rotation_swaps(WORK_WIDTH, 1) {
        circ.swap(&terminal.work2[left], &terminal.work2[right]);
    }
    for (bit, control) in terminal.l_s.iter().enumerate() {
        for (left, right) in rotation_swaps(WORK_WIDTH, 1usize << bit) {
            circ.cswap(control, &terminal.work2[left], &terminal.work2[right]);
        }
    }
}

fn restore_terminal_work2_rotation(circ: &mut Circuit, terminal: &Terminal) {
    for (bit, control) in terminal.l_s.iter().enumerate().rev() {
        let swaps = rotation_swaps(WORK_WIDTH, 1usize << bit);
        for &(left, right) in swaps.iter().rev() {
            circ.cswap(control, &terminal.work2[left], &terminal.work2[right]);
        }
    }
    let unit = rotation_swaps(WORK_WIDTH, 1);
    for &(left, right) in unit.iter().rev() {
        circ.swap(&terminal.work2[left], &terminal.work2[right]);
    }
}

fn initialize(circ: &mut Circuit, mut dx: Vec<QReg>) -> Core {
    use super::register_shared_eea_microkernels::decrement_mod_2n;
    use super::shrunken_pz_state_machine::{bit_length_lean, controlled_field_neg};
    use crate::point_add::trailmix_port::arith::compare::compare_geq_const;

    assert_eq!(dx.len(), FIELD_WIDTH);
    let iteration = circ.alloc_qreg("paper2607.iteration");
    compare_geq_const(circ, &dx, &HALF_PLUS_ONE_LE, &iteration);
    controlled_field_neg(circ, &iteration, &dx);

    let mut l_rp = circ.alloc_qreg_bits("paper2607.l-rp", LRP_WIDTH);
    l_rp.push(circ.alloc_qreg("paper2607.l-rp.high-temporary"));
    let source: Vec<&QReg> = dx.iter().take(VALUE_WIDTH).collect();
    bit_length_lean(circ, &source, &l_rp, false);
    let length_scratch = circ.alloc_qreg_bits("paper2607.length-decrement", LRP_WIDTH);
    decrement_mod_2n(circ, &l_rp, &length_scratch);
    free_clean(circ, length_scratch);
    let l_rp_high = l_rp.pop().expect("paper2607 l_rp temporary high bit");
    circ.zero_and_free(l_rp_high);
    assert_eq!(l_rp.len(), LRP_WIDTH);

    dx.push(circ.alloc_qreg("paper2607.work2-pad0"));
    dx.push(circ.alloc_qreg("paper2607.work2-pad1"));
    dx.reverse();
    let work2 = dx;

    let work1 = circ.alloc_qreg_bits("paper2607.work1", WORK_WIDTH);
    toggle_initial_work1(circ, &work1);
    let phase1 = circ.alloc_qreg("paper2607.phase1");
    let phase2 = circ.alloc_qreg("paper2607.phase2");
    let sign = circ.alloc_qreg("paper2607.sign");
    let l_t = circ.alloc_qreg_bits("paper2607.l-t", LT_WIDTH);
    let l_q = circ.alloc_qreg_bits("paper2607.l-q", LQ_WIDTH);
    let l_s = circ.alloc_qreg_bits("paper2607.l-s", SHIFT_WIDTH);
    toggle_constant(circ, &l_q, LQ_ZERO_ENCODING);
    toggle_constant(circ, &l_s, LS_ZERO_ENCODING);
    let aux = circ.alloc_qreg_bits("paper2607.aux", AUX_WIDTH);

    Core {
        phase1,
        phase2,
        iteration,
        sign,
        work1,
        work2,
        l_t,
        l_q,
        l_s,
        l_rp,
        aux,
    }
}

fn release_terminal(circ: &mut Circuit, core: Core) -> Terminal {
    toggle_terminal_work1(circ, &core.work1);
    free_clean(circ, core.work1);
    toggle_constant(circ, &core.l_t, VALUE_WIDTH - 1);
    free_clean(circ, core.l_t);
    toggle_constant(circ, &core.l_q, LQ_ZERO_ENCODING);
    free_clean(circ, core.l_q);
    toggle_constant(circ, &core.l_rp, LRP_ZERO_ENCODING);
    free_clean(circ, core.l_rp);
    circ.zero_and_free(core.phase1);
    circ.zero_and_free(core.phase2);
    circ.zero_and_free(core.sign);
    free_clean(circ, core.aux);
    Terminal {
        iteration: core.iteration,
        work2: core.work2,
        l_s: core.l_s,
    }
}

fn rebuild_terminal(circ: &mut Circuit, terminal: Terminal) -> Core {
    let work1 = circ.alloc_qreg_bits("paper2607.work1.rebuilt", WORK_WIDTH);
    toggle_terminal_work1(circ, &work1);
    let l_t = circ.alloc_qreg_bits("paper2607.l-t.rebuilt", LT_WIDTH);
    toggle_constant(circ, &l_t, VALUE_WIDTH - 1);
    let l_q = circ.alloc_qreg_bits("paper2607.l-q.rebuilt", LQ_WIDTH);
    toggle_constant(circ, &l_q, LQ_ZERO_ENCODING);
    let l_rp = circ.alloc_qreg_bits("paper2607.l-rp.rebuilt", LRP_WIDTH);
    toggle_constant(circ, &l_rp, LRP_ZERO_ENCODING);
    Core {
        phase1: circ.alloc_qreg("paper2607.phase1.rebuilt"),
        phase2: circ.alloc_qreg("paper2607.phase2.rebuilt"),
        iteration: terminal.iteration,
        sign: circ.alloc_qreg("paper2607.sign.rebuilt"),
        work1,
        work2: terminal.work2,
        l_t,
        l_q,
        l_s: terminal.l_s,
        l_rp,
        aux: circ.alloc_qreg_bits("paper2607.aux.rebuilt", AUX_WIDTH),
    }
}

fn finish(circ: &mut Circuit, mut core: Core) -> Vec<QReg> {
    use super::register_shared_eea_microkernels::increment_mod_2n;
    use super::shrunken_pz_state_machine::{bit_length_lean, controlled_field_neg};
    use crate::point_add::trailmix_port::arith::compare::compare_geq_const;

    circ.zero_and_free(core.phase1);
    circ.zero_and_free(core.phase2);
    circ.zero_and_free(core.sign);
    toggle_initial_work1(circ, &core.work1);
    free_clean(circ, core.work1);
    free_clean(circ, core.l_t);
    toggle_constant(circ, &core.l_q, LQ_ZERO_ENCODING);
    free_clean(circ, core.l_q);
    toggle_constant(circ, &core.l_s, LS_ZERO_ENCODING);
    free_clean(circ, core.l_s);
    free_clean(circ, core.aux);

    core.work2.reverse();
    let pad1 = core.work2.pop().expect("paper2607 Work2 pad1");
    let pad0 = core.work2.pop().expect("paper2607 Work2 pad0");
    circ.zero_and_free(pad1);
    circ.zero_and_free(pad0);
    assert_eq!(core.work2.len(), FIELD_WIDTH);

    core.l_rp
        .push(circ.alloc_qreg("paper2607.l-rp.high-temporary.finish"));
    let length_scratch = circ.alloc_qreg_bits("paper2607.length-increment", LRP_WIDTH);
    increment_mod_2n(circ, &core.l_rp, &length_scratch);
    free_clean(circ, length_scratch);
    let source: Vec<&QReg> = core.work2.iter().take(VALUE_WIDTH).collect();
    bit_length_lean(circ, &source, &core.l_rp, true);
    free_clean(circ, core.l_rp);

    controlled_field_neg(circ, &core.iteration, &core.work2);
    compare_geq_const(circ, &core.work2, &HALF_PLUS_ONE_LE, &core.iteration);
    circ.zero_and_free(core.iteration);
    core.work2
}

fn toggle_inverse_sign(circ: &mut Circuit, terminal: &Terminal) {
    use super::shrunken_pz_state_machine::controlled_field_neg;

    // Canonical Work2 is t' || 000, so lane 256 is already the clean field
    // top required by the 257-bit modular arithmetic interface.
    assert_eq!(terminal.work2.len(), WORK_WIDTH);
    circ.x(&terminal.iteration);
    controlled_field_neg(circ, &terminal.iteration, &terminal.work2[..FIELD_WIDTH]);
    circ.x(&terminal.iteration);
}

pub fn divide_forward(
    circ: &mut Circuit,
    dx: Vec<QReg>,
    mut dy: Vec<QReg>,
) -> (Vec<QReg>, Vec<QReg>, Vec<QReg>) {
    use super::shrunken_pz_state_machine::{
        release_q955_canonical_lambda_top, restore_q955_canonical_lambda_top,
    };
    use crate::point_add::trailmix_port::arith::rfold_mbu::mod_mul_canonical_mbu;

    assert_eq!(dx.len(), FIELD_WIDTH);
    assert_eq!(dy.len(), FIELD_WIDTH);
    let released_dy_top = loan_canonical_top(circ, &mut dy, "paper2607 forward dy");
    let core = initialize(circ, dx);
    emit_forward(circ, &core, &dy);
    let mut terminal = release_terminal(circ, core);
    canonicalize_terminal_work2(circ, &terminal);
    toggle_inverse_sign(circ, &terminal);

    restore_canonical_top(circ, &mut dy, released_dy_top);
    let mut lambda = circ.alloc_qreg_bits("paper2607.lambda", FIELD_WIDTH);
    mod_mul_canonical_mbu(circ, &lambda, &terminal.work2[..FIELD_WIDTH], &dy);
    toggle_inverse_sign(circ, &terminal);
    restore_terminal_work2_rotation(circ, &terminal);
    release_q955_canonical_lambda_top(circ, &mut lambda);

    let dy_ghosts: Vec<_> = dy.iter().map(|lane| circ.hmr_ghost(lane)).collect();
    free_clean(circ, dy);
    let core = rebuild_terminal(circ, terminal);
    emit_reverse(circ, &core, &lambda);
    let dx = finish(circ, core);

    restore_q955_canonical_lambda_top(circ, &mut lambda);
    let dy = circ.alloc_qreg_bits("paper2607.dy-restored", FIELD_WIDTH);
    mod_mul_canonical_mbu(circ, &dy, &lambda, &dx);
    for (ghost, lane) in dy_ghosts.into_iter().zip(&dy) {
        circ.resolve_ghost(ghost, lane);
    }
    (dx, dy, lambda)
}

pub fn divide_cancel(
    circ: &mut Circuit,
    dx: Vec<QReg>,
    mut dy: Vec<QReg>,
    lambda: Vec<QReg>,
) -> (Vec<QReg>, Vec<QReg>) {
    use crate::point_add::trailmix_port::arith::rfold_mbu::{
        mod_mul_canonical_mbu, mod_mul_canonical_mbu_undo,
    };

    assert_eq!(dx.len(), FIELD_WIDTH);
    assert_eq!(dy.len(), FIELD_WIDTH);
    assert_eq!(lambda.len(), FIELD_WIDTH);
    let lambda_ghosts: Vec<_> = lambda.iter().map(|lane| circ.hmr_ghost(lane)).collect();
    free_clean(circ, lambda);

    let released_forward_dy_top = loan_canonical_top(circ, &mut dy, "paper2607 cancel-forward dy");
    let core = initialize(circ, dx);
    emit_forward(circ, &core, &dy);
    let mut terminal = release_terminal(circ, core);
    canonicalize_terminal_work2(circ, &terminal);
    toggle_inverse_sign(circ, &terminal);

    restore_canonical_top(circ, &mut dy, released_forward_dy_top);
    let quotient = circ.alloc_qreg_bits("paper2607.quotient-check", FIELD_WIDTH);
    mod_mul_canonical_mbu(circ, &quotient, &terminal.work2[..FIELD_WIDTH], &dy);
    for (ghost, lane) in lambda_ghosts.into_iter().zip(&quotient) {
        circ.resolve_ghost(ghost, lane);
    }
    mod_mul_canonical_mbu_undo(circ, &quotient, &terminal.work2[..FIELD_WIDTH], &dy);
    free_clean(circ, quotient);

    toggle_inverse_sign(circ, &terminal);
    restore_terminal_work2_rotation(circ, &terminal);
    let released_reverse_dy_top = loan_canonical_top(circ, &mut dy, "paper2607 cancel-reverse dy");
    let core = rebuild_terminal(circ, terminal);
    emit_reverse(circ, &core, &dy);
    let dx = finish(circ, core);
    restore_canonical_top(circ, &mut dy, released_reverse_dy_top);
    (dx, dy)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn apply_swaps<T>(values: &mut [T], swaps: &[(usize, usize)]) {
        for &(left, right) in swaps {
            values.swap(left, right);
        }
    }

    #[test]
    fn rotation_schedule_moves_each_lane_right() {
        for shift in [1, 2, 3, 17, 128, 256] {
            let mut lanes: Vec<_> = (0..WORK_WIDTH).collect();
            apply_swaps(&mut lanes, &rotation_swaps(WORK_WIDTH, shift));
            for source in 0..WORK_WIDTH {
                assert_eq!(lanes[(source + shift) % WORK_WIDTH], source);
            }
        }
    }

    #[test]
    fn rotation_schedule_reverses_exactly() {
        for shift in [1, 2, 3, 17, 128, 256] {
            let swaps = rotation_swaps(WORK_WIDTH, shift);
            let mut lanes: Vec<_> = (0..WORK_WIDTH).collect();
            apply_swaps(&mut lanes, &swaps);
            for &(left, right) in swaps.iter().rev() {
                lanes.swap(left, right);
            }
            assert_eq!(lanes, (0..WORK_WIDTH).collect::<Vec<_>>());
        }
    }

    #[test]
    fn embedded_stream_is_complete_and_primitive() {
        let mut expected_start = 1_u32;
        let mut records = 0_usize;
        for compressed in STREAM_CHUNKS {
            let decoded = decode_chunk(compressed);
            let start = read_u32(&decoded, 16);
            let end = read_u32(&decoded, 20);
            assert_eq!(start, expected_start);
            for record in decoded[24..].chunks_exact(8).step_by(10_003) {
                let word = u64::from_le_bytes(record.try_into().expect("primitive record"));
                let kind = word & 0xf;
                let arity = (word >> 4) & 0xf;
                assert!(matches!((kind, arity), (1, 1) | (2, 2) | (3, 3) | (7, 5)));
                assert!(((word >> 8) & 0x3ff) < LOCAL_WIDTH as u64);
            }
            records += (decoded.len() - 24) / 8;
            expected_start = end + 1;
        }
        assert_eq!(expected_start, SCHEDULE_STEPS as u32 + 1);
        assert!(records > 40_000_000);
    }
}
