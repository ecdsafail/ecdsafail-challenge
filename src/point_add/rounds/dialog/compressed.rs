//! Dialog-GCD compressed-sidecar path: the round763 block compressor, the
//! runway / composite scratch layout helpers, and the
//! `emit_dialog_gcd_compressed_sidecar_*` block-lifecycle emitters
//! (tobitvector / apply / ipmul / quotient). An alternate, lower-peak encoding
//! of the GCD transcript log; shares the raw-path config levers and comparators
//! from the parent `dialog` module.
use super::*;

pub(crate) fn round763_dedup_enabled() -> bool {
    // EXACT rewrite: the pair ccx(1,3->4) ... ccx(1,3->4) bracketing cx(1->0)
    // cancels (nothing between them touches 1/3/4), so it reduces to bare cx(1->0).
    // 2 CCX -> 0 per direction x ~1064 sites. Default OFF (op-stream reseed).
    std::env::var("DIALOG_GCD_ROUND763_DEDUP").ok().as_deref() == Some("1")
}

pub(crate) fn round763_compress_lever_enabled() -> bool {
    // Reachable-support rewrite of the round763 6->5 sidecar packer. Each raw
    // slot is (b0, b0_and_b1), with b0_and_b1 = b0 & (v<u), so state (0,1) is
    // unreachable on the verifier support. On that support, three CCX collapse
    // to CX and the compressor drops from 9 CCX to 4 CCX per direction.
    std::env::var("DIALOG_GCD_ROUND763_COMPRESS_LEVER")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn emit_dialog_gcd_round763_compressor(b: &mut B, block: &[QubitId]) {
    assert_eq!(block.len(), 6);
    if round763_compress_lever_enabled() {
        b.cx(block[5], block[3]);
        b.ccx(block[3], block[4], block[5]);
        b.cx(block[1], block[4]);
        b.cx(block[1], block[0]);
        b.ccx(block[4], block[5], block[1]);
        b.cx(block[0], block[2]);
        b.ccx(block[2], block[5], block[0]);
        b.ccx(block[0], block[1], block[5]);
        return;
    }
    b.ccx(block[4], block[5], block[3]);
    b.ccx(block[3], block[4], block[5]);
    b.ccx(block[1], block[2], block[4]);
    if round763_dedup_enabled() {
        b.cx(block[1], block[0]);
    } else {
        b.ccx(block[1], block[3], block[4]);
        b.cx(block[1], block[0]);
        b.ccx(block[1], block[3], block[4]);
    }
    b.ccx(block[4], block[5], block[1]);
    b.ccx(block[0], block[5], block[2]);
    b.ccx(block[2], block[5], block[0]);
    b.ccx(block[0], block[1], block[5]);
}

pub(crate) fn emit_dialog_gcd_round763_compressor_inverse(b: &mut B, block: &[QubitId]) {
    assert_eq!(block.len(), 6);
    if round763_compress_lever_enabled() {
        b.ccx(block[0], block[1], block[5]);
        b.ccx(block[2], block[5], block[0]);
        b.cx(block[0], block[2]);
        b.ccx(block[4], block[5], block[1]);
        b.cx(block[1], block[0]);
        b.cx(block[1], block[4]);
        b.ccx(block[3], block[4], block[5]);
        b.cx(block[5], block[3]);
        return;
    }
    b.ccx(block[0], block[1], block[5]);
    b.ccx(block[2], block[5], block[0]);
    b.ccx(block[0], block[5], block[2]);
    b.ccx(block[4], block[5], block[1]);
    if round763_dedup_enabled() {
        b.cx(block[1], block[0]);
    } else {
        b.ccx(block[1], block[3], block[4]);
        b.cx(block[1], block[0]);
        b.ccx(block[1], block[3], block[4]);
    }
    b.ccx(block[1], block[2], block[4]);
    b.ccx(block[3], block[4], block[5]);
    b.ccx(block[4], block[5], block[3]);
}

const DIALOG_GCD_K5_DATA_WIRES: [usize; 12] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12];
const DIALOG_GCD_K5_HEAD11_DATA_WIRES: [usize; 11] =
    [0, 1, 2, 4, 5, 6, 7, 8, 9, 11, 12];
#[derive(Clone, Copy)]
enum DialogGcdK5FableGate {
    X(usize),
    Cx(usize, usize),
    Ccx(usize, usize, usize),
}

const DIALOG_GCD_K5_FABLE_GATES: &[DialogGcdK5FableGate] = &[
    DialogGcdK5FableGate::Cx(7, 6),
    DialogGcdK5FableGate::Cx(6, 7),
    DialogGcdK5FableGate::Cx(7, 6),
    DialogGcdK5FableGate::Cx(8, 7),
    DialogGcdK5FableGate::Cx(7, 8),
    DialogGcdK5FableGate::Cx(9, 6),
    DialogGcdK5FableGate::X(10),
    DialogGcdK5FableGate::X(6),
    DialogGcdK5FableGate::Ccx(10, 6, 9),
    DialogGcdK5FableGate::X(6),
    DialogGcdK5FableGate::X(10),
    DialogGcdK5FableGate::Cx(9, 6),
    DialogGcdK5FableGate::Cx(11, 10),
    DialogGcdK5FableGate::X(6),
    DialogGcdK5FableGate::X(10),
    DialogGcdK5FableGate::Cx(5, 11),
    DialogGcdK5FableGate::Ccx(6, 10, 5),
    DialogGcdK5FableGate::Cx(5, 11),
    DialogGcdK5FableGate::X(10),
    DialogGcdK5FableGate::X(6),
    DialogGcdK5FableGate::Cx(11, 10),
    DialogGcdK5FableGate::Cx(10, 4),
    DialogGcdK5FableGate::X(6),
    DialogGcdK5FableGate::Cx(10, 11),
    DialogGcdK5FableGate::Ccx(6, 4, 10),
    DialogGcdK5FableGate::Cx(10, 11),
    DialogGcdK5FableGate::X(6),
    DialogGcdK5FableGate::Cx(10, 4),
    DialogGcdK5FableGate::Cx(7, 4),
    DialogGcdK5FableGate::X(10),
    DialogGcdK5FableGate::Cx(5, 7),
    DialogGcdK5FableGate::Cx(5, 9),
    DialogGcdK5FableGate::Ccx(10, 4, 5),
    DialogGcdK5FableGate::Cx(5, 9),
    DialogGcdK5FableGate::Cx(5, 7),
    DialogGcdK5FableGate::X(10),
    DialogGcdK5FableGate::Cx(7, 4),
    DialogGcdK5FableGate::Cx(10, 4),
    DialogGcdK5FableGate::Cx(7, 5),
    DialogGcdK5FableGate::X(4),
    DialogGcdK5FableGate::Cx(9, 11),
    DialogGcdK5FableGate::Ccx(4, 5, 9),
    DialogGcdK5FableGate::Cx(9, 11),
    DialogGcdK5FableGate::X(4),
    DialogGcdK5FableGate::Cx(7, 5),
    DialogGcdK5FableGate::Cx(10, 4),
    DialogGcdK5FableGate::Cx(9, 8),
    DialogGcdK5FableGate::X(10),
    DialogGcdK5FableGate::Ccx(10, 8, 9),
    DialogGcdK5FableGate::X(10),
    DialogGcdK5FableGate::Cx(9, 8),
    DialogGcdK5FableGate::Cx(4, 8),
    DialogGcdK5FableGate::Cx(11, 4),
    DialogGcdK5FableGate::X(8),
    DialogGcdK5FableGate::Cx(5, 11),
    DialogGcdK5FableGate::Ccx(8, 4, 5),
    DialogGcdK5FableGate::Cx(5, 11),
    DialogGcdK5FableGate::X(8),
    DialogGcdK5FableGate::Cx(11, 4),
    DialogGcdK5FableGate::Cx(4, 8),
    DialogGcdK5FableGate::Cx(7, 4),
    DialogGcdK5FableGate::X(4),
    DialogGcdK5FableGate::Ccx(5, 4, 7),
    DialogGcdK5FableGate::X(4),
    DialogGcdK5FableGate::Cx(7, 4),
    DialogGcdK5FableGate::X(8),
    DialogGcdK5FableGate::Ccx(4, 8, 6),
    DialogGcdK5FableGate::X(8),
    DialogGcdK5FableGate::X(6),
    DialogGcdK5FableGate::X(11),
    DialogGcdK5FableGate::Cx(4, 5),
    DialogGcdK5FableGate::Cx(4, 8),
    DialogGcdK5FableGate::Cx(4, 10),
    DialogGcdK5FableGate::Ccx(6, 11, 4),
    DialogGcdK5FableGate::Cx(4, 10),
    DialogGcdK5FableGate::Cx(4, 8),
    DialogGcdK5FableGate::Cx(4, 5),
    DialogGcdK5FableGate::X(11),
    DialogGcdK5FableGate::X(6),
    DialogGcdK5FableGate::Cx(9, 5),
    DialogGcdK5FableGate::X(10),
    DialogGcdK5FableGate::Ccx(10, 5, 9),
    DialogGcdK5FableGate::X(10),
    DialogGcdK5FableGate::Cx(9, 5),
    DialogGcdK5FableGate::X(8),
    DialogGcdK5FableGate::X(9),
    DialogGcdK5FableGate::Ccx(7, 8, 13),
    DialogGcdK5FableGate::Cx(4, 6),
    DialogGcdK5FableGate::Ccx(13, 9, 4),
    DialogGcdK5FableGate::Cx(4, 6),
    DialogGcdK5FableGate::Ccx(7, 8, 13),
    DialogGcdK5FableGate::X(9),
    DialogGcdK5FableGate::X(8),
    DialogGcdK5FableGate::X(4),
    DialogGcdK5FableGate::X(6),
    DialogGcdK5FableGate::X(11),
    DialogGcdK5FableGate::Ccx(4, 6, 13),
    DialogGcdK5FableGate::Ccx(13, 11, 10),
    DialogGcdK5FableGate::Ccx(4, 6, 13),
    DialogGcdK5FableGate::X(11),
    DialogGcdK5FableGate::X(6),
    DialogGcdK5FableGate::X(4),
    DialogGcdK5FableGate::X(10),
];

fn dialog_gcd_k5_fable_wire(data: &[QubitId; 13], ancilla: QubitId, wire: usize) -> QubitId {
    if wire == 13 {
        ancilla
    } else {
        debug_assert!(wire < data.len());
        data[wire]
    }
}

fn dialog_gcd_k5_emit_fable_gate(
    b: &mut B,
    data: &[QubitId; 13],
    ancilla: QubitId,
    gate: DialogGcdK5FableGate,
) {
    match gate {
        DialogGcdK5FableGate::X(a) => b.x(dialog_gcd_k5_fable_wire(data, ancilla, a)),
        DialogGcdK5FableGate::Cx(a, c) => b.cx(
            dialog_gcd_k5_fable_wire(data, ancilla, a),
            dialog_gcd_k5_fable_wire(data, ancilla, c),
        ),
        DialogGcdK5FableGate::Ccx(a, c, t) => b.ccx(
            dialog_gcd_k5_fable_wire(data, ancilla, a),
            dialog_gcd_k5_fable_wire(data, ancilla, c),
            dialog_gcd_k5_fable_wire(data, ancilla, t),
        ),
    }
}

fn dialog_gcd_k5_emit_fable_codec(
    b: &mut B,
    data: &[QubitId; 13],
    ancilla: QubitId,
    inverse: bool,
) {
    if inverse {
        for &gate in DIALOG_GCD_K5_FABLE_GATES.iter().rev() {
            dialog_gcd_k5_emit_fable_gate(b, data, ancilla, gate);
        }
    } else {
        for &gate in DIALOG_GCD_K5_FABLE_GATES {
            dialog_gcd_k5_emit_fable_gate(b, data, ancilla, gate);
        }
    }
}

fn emit_dialog_gcd_k5_clean_compressor(b: &mut B, data: &[QubitId; 13], ancilla: QubitId) {
    dialog_gcd_k5_emit_fable_codec(b, data, ancilla, false);
}

fn emit_dialog_gcd_k5_clean_compressor_inverse(
    b: &mut B,
    data: &[QubitId; 13],
    ancilla: QubitId,
) {
    dialog_gcd_k5_emit_fable_codec(b, data, ancilla, true);
}

fn dialog_gcd_k5_head11_enabled() -> bool {
    dialog_gcd_k5_clean_block_enabled()
        && dialog_gcd_active_iterations() >= 5
        && std::env::var("DIALOG_GCD_K5_HEAD11_CODEC")
            .ok()
            .as_deref()
            == Some("1")
}

fn dialog_gcd_k5_head11_apply_step_frame_enabled() -> bool {
    dialog_gcd_k5_head11_enabled()
        && std::env::var("DIALOG_GCD_K5_HEAD11_APPLY_STEP_FRAME")
            .ok()
            .as_deref()
            == Some("1")
}

fn emit_dialog_gcd_k5_head11_preconditioner(b: &mut B, data: &[QubitId; 13]) {
    b.x(data[0]);
    b.ccx(data[0], data[1], data[3]);
    b.ccx(data[2], data[3], data[0]);
    b.cx(data[0], data[3]);
}

fn emit_dialog_gcd_k5_head11_preconditioner_inverse(
    b: &mut B,
    data: &[QubitId; 13],
) {
    b.cx(data[0], data[3]);
    b.ccx(data[2], data[3], data[0]);
    b.ccx(data[0], data[1], data[3]);
    b.x(data[0]);
}

fn emit_dialog_gcd_k5_pair_encoder(b: &mut B, pair_raw: &[QubitId; 6]) {
    let core = [pair_raw[0], pair_raw[1], pair_raw[4], pair_raw[2], pair_raw[3]];
    b.cx(core[1], core[2]);
    b.cx(core[0], core[4]);
    b.x(core[3]);
    b.ccx(core[2], core[3], core[1]);
    b.cx(core[3], core[4]);
    b.ccx(core[3], core[4], core[0]);
    b.cx(core[2], core[4]);
    b.cx(core[0], core[3]);
    b.cx(core[3], core[2]);
    b.cx(core[3], core[4]);
    b.ccx(core[1], core[3], core[0]);
    b.cx(core[1], core[0]);
    b.cx(core[3], core[0]);
}

fn emit_dialog_gcd_k5_pair_encoder_inverse(b: &mut B, pair_raw: &[QubitId; 6]) {
    let core = [pair_raw[0], pair_raw[1], pair_raw[4], pair_raw[2], pair_raw[3]];
    b.cx(core[3], core[0]);
    b.cx(core[1], core[0]);
    b.ccx(core[1], core[3], core[0]);
    b.cx(core[3], core[4]);
    b.cx(core[3], core[2]);
    b.cx(core[0], core[3]);
    b.cx(core[2], core[4]);
    b.ccx(core[3], core[4], core[0]);
    b.cx(core[3], core[4]);
    b.ccx(core[2], core[3], core[1]);
    b.x(core[3]);
    b.cx(core[0], core[4]);
    b.cx(core[1], core[2]);
}

fn dialog_gcd_raw_s2(raw_block: &[QubitId], slot: usize) -> QubitId {
    raw_block[2 * dialog_gcd_sidecar_group_size() + slot]
}

fn dialog_gcd_block_raw_s2(
    raw_block: &[QubitId],
    block_steps: usize,
    slot: usize,
) -> QubitId {
    if dialog_gcd_k5_tail6_graph9_enabled() && block_steps == 6 {
        assert!(slot < DIALOG_GCD_K5_TAIL6_GRAPH9_STORED_STEPS);
        raw_block[2 * DIALOG_GCD_K5_TAIL6_GRAPH9_STORED_STEPS + slot]
    } else if dialog_gcd_k5_tail6_graph_enabled() && block_steps == 6 {
        assert!(slot < DIALOG_GCD_K5_TAIL6_GRAPH_STORED_STEPS);
        raw_block[2 * DIALOG_GCD_K5_TAIL6_GRAPH_STORED_STEPS + slot]
    } else if dialog_gcd_k5_tail7_enabled() && block_steps == 7 {
        assert!(slot < DIALOG_GCD_K5_TAIL7_STORED_STEPS);
        raw_block[2 * DIALOG_GCD_K5_TAIL7_STORED_STEPS + slot]
    } else {
        dialog_gcd_raw_s2(raw_block, slot)
    }
}

fn dialog_gcd_k5_pair01(raw_block: &[QubitId]) -> [QubitId; 6] {
    [
        raw_block[0],
        raw_block[1],
        raw_block[2],
        raw_block[3],
        dialog_gcd_raw_s2(raw_block, 0),
        dialog_gcd_raw_s2(raw_block, 1),
    ]
}

fn dialog_gcd_k5_pair23(raw_block: &[QubitId]) -> [QubitId; 6] {
    [
        raw_block[4],
        raw_block[5],
        raw_block[6],
        raw_block[7],
        dialog_gcd_raw_s2(raw_block, 2),
        dialog_gcd_raw_s2(raw_block, 3),
    ]
}

fn dialog_gcd_k5_data_from_raw(raw_block: &[QubitId]) -> [QubitId; 13] {
    [
        raw_block[1],
        dialog_gcd_raw_s2(raw_block, 0),
        raw_block[2],
        raw_block[3],
        dialog_gcd_raw_s2(raw_block, 1),
        raw_block[5],
        dialog_gcd_raw_s2(raw_block, 2),
        raw_block[6],
        raw_block[7],
        dialog_gcd_raw_s2(raw_block, 3),
        raw_block[8],
        raw_block[9],
        dialog_gcd_raw_s2(raw_block, 4),
    ]
}

fn dialog_gcd_k5_transfer_survivors(
    b: &mut B,
    compressed_block: &[QubitId],
    data: &[QubitId; 13],
    swap_host: bool,
) {
    assert_eq!(compressed_block.len(), 12);
    for (i, &wire) in DIALOG_GCD_K5_DATA_WIRES.iter().enumerate() {
        if swap_host {
            b.swap(compressed_block[i], data[wire]);
        } else {
            b.cx(compressed_block[i], data[wire]);
        }
    }
}

fn dialog_gcd_k5_head11_transfer_survivors(
    b: &mut B,
    compressed_block: &[QubitId],
    data: &[QubitId; 13],
    swap_host: bool,
) {
    assert_eq!(compressed_block.len(), DIALOG_GCD_K5_HEAD11_DATA_WIRES.len());
    for (i, &wire) in DIALOG_GCD_K5_HEAD11_DATA_WIRES.iter().enumerate() {
        if swap_host {
            b.swap(compressed_block[i], data[wire]);
        } else {
            b.cx(compressed_block[i], data[wire]);
        }
    }
}

fn dialog_gcd_k5_head11_compress_raw_to_block(
    b: &mut B,
    compressed_block: &[QubitId],
    raw_block: &[QubitId],
    swap_host: bool,
) {
    assert_eq!(compressed_block.len(), DIALOG_GCD_K5_HEAD11_DATA_WIRES.len());
    assert_eq!(raw_block.len(), 15);
    emit_dialog_gcd_k5_pair_encoder(b, &dialog_gcd_k5_pair01(raw_block));
    emit_dialog_gcd_k5_pair_encoder(b, &dialog_gcd_k5_pair23(raw_block));
    let data = dialog_gcd_k5_data_from_raw(raw_block);
    emit_dialog_gcd_k5_head11_preconditioner(b, &data);
    let ancilla = b.alloc_qubit();
    emit_dialog_gcd_k5_clean_compressor(b, &data, ancilla);
    b.free(ancilla);
    dialog_gcd_k5_head11_transfer_survivors(b, compressed_block, &data, swap_host);
}

fn dialog_gcd_k5_head11_decompress_block_to_raw(
    b: &mut B,
    compressed_block: &[QubitId],
    raw_block: &[QubitId],
    swap_host: bool,
) {
    assert_eq!(compressed_block.len(), DIALOG_GCD_K5_HEAD11_DATA_WIRES.len());
    assert_eq!(raw_block.len(), 15);
    let data = dialog_gcd_k5_data_from_raw(raw_block);
    dialog_gcd_k5_head11_transfer_survivors(b, compressed_block, &data, swap_host);
    let ancilla = b.alloc_qubit();
    emit_dialog_gcd_k5_clean_compressor_inverse(b, &data, ancilla);
    b.free(ancilla);
    emit_dialog_gcd_k5_head11_preconditioner_inverse(b, &data);
    emit_dialog_gcd_k5_pair_encoder_inverse(b, &dialog_gcd_k5_pair23(raw_block));
    emit_dialog_gcd_k5_pair_encoder_inverse(b, &dialog_gcd_k5_pair01(raw_block));
}

fn dialog_gcd_k5_head11_data_from_inplace_block(
    compressed_block: &[QubitId],
    data3: QubitId,
    data10: QubitId,
) -> [QubitId; 13] {
    assert_eq!(compressed_block.len(), DIALOG_GCD_K5_HEAD11_DATA_WIRES.len());
    let mut data = [data3; 13];
    for (i, &wire) in DIALOG_GCD_K5_HEAD11_DATA_WIRES.iter().enumerate() {
        data[wire] = compressed_block[i];
    }
    data[3] = data3;
    data[10] = data10;
    data
}

fn dialog_gcd_k5_copy_slot_to_frame(
    b: &mut B,
    frame: &[QubitId],
    b0: QubitId,
    b0_and_b1: QubitId,
    s2: QubitId,
) {
    assert_eq!(frame.len(), 3);
    b.cx(b0, frame[0]);
    b.cx(b0_and_b1, frame[1]);
    b.cx(s2, frame[2]);
}

fn dialog_gcd_k5_head11_toggle_step_frame(
    b: &mut B,
    compressed_block: &[QubitId],
    frame: &[QubitId],
    slot: usize,
) {
    assert_eq!(compressed_block.len(), DIALOG_GCD_K5_HEAD11_DATA_WIRES.len());
    assert_eq!(frame.len(), 3);
    assert!(slot < 5);
    let data3 = b.alloc_qubit();
    let data10 = b.alloc_qubit();
    let pair0 = b.alloc_qubit();
    let ancilla = b.alloc_qubit();
    let data = dialog_gcd_k5_head11_data_from_inplace_block(compressed_block, data3, data10);
    emit_dialog_gcd_k5_clean_compressor_inverse(b, &data, ancilla);
    emit_dialog_gcd_k5_head11_preconditioner_inverse(b, &data);
    match slot {
        0 | 1 => {
            let pair = [pair0, data[0], data[2], data[3], data[1], data[4]];
            emit_dialog_gcd_k5_pair_encoder_inverse(b, &pair);
            if slot == 0 {
                dialog_gcd_k5_copy_slot_to_frame(b, frame, pair[0], pair[0], pair[4]);
            } else {
                dialog_gcd_k5_copy_slot_to_frame(b, frame, pair[2], pair[3], pair[5]);
            }
            emit_dialog_gcd_k5_pair_encoder(b, &pair);
        }
        2 | 3 => {
            let pair = [pair0, data[5], data[7], data[8], data[6], data[9]];
            emit_dialog_gcd_k5_pair_encoder_inverse(b, &pair);
            if slot == 2 {
                dialog_gcd_k5_copy_slot_to_frame(b, frame, pair[0], pair[1], pair[4]);
            } else {
                dialog_gcd_k5_copy_slot_to_frame(b, frame, pair[2], pair[3], pair[5]);
            }
            emit_dialog_gcd_k5_pair_encoder(b, &pair);
        }
        4 => {
            dialog_gcd_k5_copy_slot_to_frame(b, frame, data[10], data[11], data[12]);
        }
        _ => unreachable!(),
    }
    emit_dialog_gcd_k5_head11_preconditioner(b, &data);
    emit_dialog_gcd_k5_clean_compressor(b, &data, ancilla);
    b.free(ancilla);
    b.free(pair0);
    b.free(data10);
    b.free(data3);
}

fn dialog_gcd_k5_head11_pair_encode_word(bits: &mut [bool; 15], slots: [usize; 2]) {
    let wire = [
        3 * slots[0],
        3 * slots[0] + 1,
        3 * slots[0] + 2,
        3 * slots[1],
        3 * slots[1] + 1,
    ];
    bits[wire[2]] ^= bits[wire[1]];
    bits[wire[4]] ^= bits[wire[0]];
    bits[wire[3]] ^= true;
    bits[wire[1]] ^= bits[wire[2]] && bits[wire[3]];
    bits[wire[4]] ^= bits[wire[3]];
    bits[wire[0]] ^= bits[wire[3]] && bits[wire[4]];
    bits[wire[4]] ^= bits[wire[2]];
    bits[wire[3]] ^= bits[wire[0]];
    bits[wire[2]] ^= bits[wire[3]];
    bits[wire[4]] ^= bits[wire[3]];
    bits[wire[0]] ^= bits[wire[1]] && bits[wire[3]];
    bits[wire[0]] ^= bits[wire[1]];
    bits[wire[0]] ^= bits[wire[3]];
}

fn dialog_gcd_k5_head11_code_word(pattern: u16) -> Option<u16> {
    let mut raw = std::array::from_fn::<_, 15, _>(|bit| (pattern >> bit) & 1 != 0);
    dialog_gcd_k5_head11_pair_encode_word(&mut raw, [0, 1]);
    dialog_gcd_k5_head11_pair_encode_word(&mut raw, [2, 3]);
    if raw[0] || raw[6] {
        return None;
    }

    const RAW_DATA_INDICES: [usize; 13] =
        [1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14];
    let mut wires = [false; 14];
    for (index, raw_index) in RAW_DATA_INDICES.into_iter().enumerate() {
        wires[index] = raw[raw_index];
    }
    wires[0] ^= true;
    wires[3] ^= wires[0] && wires[1];
    wires[0] ^= wires[2] && wires[3];
    wires[3] ^= wires[0];
    for &gate in DIALOG_GCD_K5_FABLE_GATES {
        match gate {
            DialogGcdK5FableGate::X(a) => wires[a] ^= true,
            DialogGcdK5FableGate::Cx(a, c) => wires[c] ^= wires[a],
            DialogGcdK5FableGate::Ccx(a, c, t) => wires[t] ^= wires[a] && wires[c],
        }
    }
    if wires[3] || wires[10] || wires[13] {
        return None;
    }
    Some(
        DIALOG_GCD_K5_HEAD11_DATA_WIRES
            .iter()
            .enumerate()
            .fold(0u16, |code, (index, &wire)| {
                code | (u16::from(wires[wire]) << index)
            }),
    )
}

pub(crate) fn dialog_gcd_k5_head11_supports(pattern: u16) -> bool {
    dialog_gcd_k5_head11_code_word(pattern).is_some()
}

pub(crate) fn dialog_gcd_k5_head11_codec_selftest() -> Result<(), String> {
    use sha3::digest::{ExtendableOutput, Update};

    let supported = (0u16..1 << 15)
        .filter(|&pattern| dialog_gcd_k5_head11_supports(pattern))
        .collect::<Vec<_>>();
    if supported.len() != 1 << DIALOG_GCD_K5_HEAD11_DATA_WIRES.len() {
        return Err(format!(
            "expected 2048 supported head words, got {}",
            supported.len()
        ));
    }
    let mut seen_codes = vec![false; 1 << DIALOG_GCD_K5_HEAD11_DATA_WIRES.len()];
    for &pattern in &supported {
        let code = dialog_gcd_k5_head11_code_word(pattern).expect("filtered support");
        if std::mem::replace(&mut seen_codes[code as usize], true) {
            return Err(format!("duplicate head code 0x{code:03x}"));
        }
    }

    let build_codec = |decompress: bool| {
        let mut b = B::new();
        let code = b.alloc_qubits(DIALOG_GCD_K5_HEAD11_DATA_WIRES.len());
        let raw = b.alloc_qubits(15);
        if decompress {
            dialog_gcd_k5_head11_decompress_block_to_raw(&mut b, &code, &raw, true);
        } else {
            dialog_gcd_k5_head11_compress_raw_to_block(&mut b, &code, &raw, true);
        }
        (b.ops, code, raw, b.next_qubit as usize, b.next_bit as usize)
    };
    let forward_codec = build_codec(false);
    let reverse_codec = build_codec(true);

    for batch_start in (0..supported.len()).step_by(64) {
        let patterns = &supported[batch_start..batch_start + 64];
        let mut raw_masks = [0u64; 15];
        let mut code_masks = [0u64; DIALOG_GCD_K5_HEAD11_DATA_WIRES.len()];
        for (shot, &pattern) in patterns.iter().enumerate() {
            let shot_bit = 1u64 << shot;
            for slot in 0..5 {
                if (pattern >> (3 * slot)) & 1 != 0 {
                    raw_masks[2 * slot] |= shot_bit;
                }
                if (pattern >> (3 * slot + 1)) & 1 != 0 {
                    raw_masks[2 * slot + 1] |= shot_bit;
                }
                if (pattern >> (3 * slot + 2)) & 1 != 0 {
                    raw_masks[10 + slot] |= shot_bit;
                }
            }
            let code = dialog_gcd_k5_head11_code_word(pattern).expect("supported pattern");
            for (index, mask) in code_masks.iter_mut().enumerate() {
                if (code >> index) & 1 != 0 {
                    *mask |= shot_bit;
                }
            }
        }

        let run = |decompress: bool| {
            let (ops, code, raw, num_qubits, num_bits) =
                if decompress { &reverse_codec } else { &forward_codec };
            let mut seed = sha3::Shake128::default();
            seed.update(b"dialog-gcd-k5-head11-codec-selftest");
            seed.update(&(batch_start as u64).to_le_bytes());
            seed.update(&[u8::from(decompress)]);
            let mut xof = seed.finalize_xof();
            let mut sim = Simulator::new(*num_qubits, *num_bits, &mut xof);
            sim.clear_for_shot();
            let source = if decompress {
                &code_masks[..]
            } else {
                &raw_masks[..]
            };
            let targets = if decompress { &code[..] } else { &raw[..] };
            for (&qubit, &mask) in targets.iter().zip(source.iter()) {
                *sim.qubit_mut(qubit) = mask;
            }
            sim.apply_iter(ops.iter());
            (
                code.iter().map(|&q| sim.qubit(q)).collect::<Vec<_>>(),
                raw.iter().map(|&q| sim.qubit(q)).collect::<Vec<_>>(),
                sim.phase,
            )
        };

        let (forward_code, forward_raw, forward_phase) = run(false);
        if forward_phase != 0 {
            return Err(format!(
                "forward phase garbage in batch {batch_start}: 0x{forward_phase:x}"
            ));
        }
        if forward_code != code_masks {
            return Err(format!(
                "forward code mismatch in batch {batch_start}: got {forward_code:x?}, want {code_masks:x?}"
            ));
        }
        if forward_raw.iter().any(|&mask| mask != 0) {
            return Err(format!(
                "forward raw garbage in batch {batch_start}: {forward_raw:x?}"
            ));
        }

        let (reverse_code, reverse_raw, reverse_phase) = run(true);
        if reverse_phase != 0 {
            return Err(format!(
                "reverse phase garbage in batch {batch_start}: 0x{reverse_phase:x}"
            ));
        }
        if reverse_code.iter().any(|&mask| mask != 0) {
            return Err(format!(
                "reverse code garbage in batch {batch_start}: {reverse_code:x?}"
            ));
        }
        if reverse_raw != raw_masks {
            return Err(format!(
                "reverse raw mismatch in batch {batch_start}: got {reverse_raw:x?}, want {raw_masks:x?}"
            ));
        }
    }
    Ok(())
}

fn dialog_gcd_k5_compress_raw_to_block(
    b: &mut B,
    compressed_block: &[QubitId],
    raw_block: &[QubitId],
    swap_host: bool,
) {
    assert_eq!(compressed_block.len(), 12);
    assert_eq!(raw_block.len(), 15);
    emit_dialog_gcd_k5_pair_encoder(b, &dialog_gcd_k5_pair01(raw_block));
    emit_dialog_gcd_k5_pair_encoder(b, &dialog_gcd_k5_pair23(raw_block));
    let data = dialog_gcd_k5_data_from_raw(raw_block);
    let ancilla = b.alloc_qubit();
    emit_dialog_gcd_k5_clean_compressor(b, &data, ancilla);
    b.free(ancilla);
    dialog_gcd_k5_transfer_survivors(b, compressed_block, &data, swap_host);
}

fn dialog_gcd_k5_decompress_block_to_raw(
    b: &mut B,
    compressed_block: &[QubitId],
    raw_block: &[QubitId],
    swap_host: bool,
) {
    assert_eq!(compressed_block.len(), 12);
    assert_eq!(raw_block.len(), 15);
    let data = dialog_gcd_k5_data_from_raw(raw_block);
    dialog_gcd_k5_transfer_survivors(b, compressed_block, &data, swap_host);
    let ancilla = b.alloc_qubit();
    emit_dialog_gcd_k5_clean_compressor_inverse(b, &data, ancilla);
    b.free(ancilla);
    emit_dialog_gcd_k5_pair_encoder_inverse(b, &dialog_gcd_k5_pair23(raw_block));
    emit_dialog_gcd_k5_pair_encoder_inverse(b, &dialog_gcd_k5_pair01(raw_block));
}

fn dialog_gcd_k5_compress_partial_raw_to_block(
    b: &mut B,
    compressed_block: &[QubitId],
    raw_block: &[QubitId],
    steps: usize,
    swap_host: bool,
) {
    assert!(steps <= DIALOG_GCD_HIGH_TAIL_ALIAS_GROUP_SIZE);
    assert_eq!(compressed_block.len(), 12);
    assert_eq!(raw_block.len(), 15);
    let base_bits = DIALOG_GCD_HIGH_TAIL_ALIAS_BLOCK_BITS;
    let raw_base = 2 * DIALOG_GCD_HIGH_TAIL_ALIAS_GROUP_SIZE;
    emit_dialog_gcd_round763_compressor(b, &raw_block[0..raw_base]);
    for i in 0..base_bits {
        if swap_host { b.swap(compressed_block[i], raw_block[i]); } else { b.cx(compressed_block[i], raw_block[i]); }
    }
    for slot in 0..steps {
        let s2 = dialog_gcd_raw_s2(raw_block, slot);
        if swap_host { b.swap(compressed_block[base_bits + slot], s2); } else { b.cx(compressed_block[base_bits + slot], s2); }
    }
}

fn dialog_gcd_k5_decompress_partial_block_to_raw(
    b: &mut B,
    compressed_block: &[QubitId],
    raw_block: &[QubitId],
    steps: usize,
    swap_host: bool,
) {
    assert!(steps <= DIALOG_GCD_HIGH_TAIL_ALIAS_GROUP_SIZE);
    assert_eq!(compressed_block.len(), 12);
    assert_eq!(raw_block.len(), 15);
    let base_bits = DIALOG_GCD_HIGH_TAIL_ALIAS_BLOCK_BITS;
    let raw_base = 2 * DIALOG_GCD_HIGH_TAIL_ALIAS_GROUP_SIZE;
    for i in 0..base_bits {
        if swap_host { b.swap(compressed_block[i], raw_block[i]); } else { b.cx(compressed_block[i], raw_block[i]); }
    }
    emit_dialog_gcd_round763_compressor_inverse(b, &raw_block[0..raw_base]);
    for slot in 0..steps {
        let s2 = dialog_gcd_raw_s2(raw_block, slot);
        if swap_host { b.swap(compressed_block[base_bits + slot], s2); } else { b.cx(compressed_block[base_bits + slot], s2); }
    }
}

fn dialog_gcd_k5_tail_pair1_enabled() -> bool {
    dialog_gcd_k5_clean_block_enabled()
        && !dialog_gcd_k5_tail7_enabled()
        && !dialog_gcd_k5_tail6_graph_enabled()
        && !dialog_gcd_k5_tail6_graph9_enabled()
        && dialog_gcd_active_iterations() % dialog_gcd_sidecar_group_size() == 2
        && std::env::var("DIALOG_GCD_K5_TAIL_PAIR1")
            .ok()
            .as_deref()
            == Some("1")
}

fn dialog_gcd_k5_tail6_graph9_enabled() -> bool {
    dialog_gcd_k5_clean_block_enabled()
        && dialog_gcd_active_iterations() >= 6
        && dialog_gcd_active_iterations() % dialog_gcd_sidecar_group_size() == 1
        && std::env::var("DIALOG_GCD_K5_TAIL6_GRAPH9_CODEC")
            .ok()
            .as_deref()
            == Some("1")
}

fn dialog_gcd_k5_release_decoded_block_bits() -> usize {
    if !dialog_gcd_k5_clean_block_enabled() || !dialog_gcd_apply_replay_swap_host_enabled() {
        return 0;
    }
    std::env::var("DIALOG_GCD_K5_RELEASE_DECODED_BLOCK_BITS")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(0)
}

fn dialog_gcd_k5_release_decoded_tail_bits() -> usize {
    std::env::var("DIALOG_GCD_K5_RELEASE_DECODED_TAIL_BITS")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or_else(dialog_gcd_k5_release_decoded_block_bits)
}

fn dialog_gcd_k5_release_scale_bits() -> usize {
    std::env::var("DIALOG_GCD_K5_RELEASE_SCALE_BITS")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(0)
}

fn dialog_gcd_k5_tail6_graph_enabled() -> bool {
    dialog_gcd_k5_clean_block_enabled()
        && dialog_gcd_active_iterations() >= 6
        && dialog_gcd_active_iterations() % dialog_gcd_sidecar_group_size() == 1
        && std::env::var("DIALOG_GCD_K5_TAIL6_GRAPH_CODEC")
            .ok()
            .as_deref()
            == Some("1")
}

fn dialog_gcd_k5_tail7_enabled() -> bool {
    dialog_gcd_k5_clean_block_enabled()
        && dialog_gcd_active_iterations() >= 7
        && dialog_gcd_active_iterations() % dialog_gcd_sidecar_group_size() == 2
        && std::env::var("DIALOG_GCD_K5_TAIL7_CODEC")
            .ok()
            .as_deref()
            == Some("1")
}

const DIALOG_GCD_K5_TAIL6_GRAPH_STORED_STEPS: usize = 3;
const DIALOG_GCD_K5_TAIL6_GRAPH_CODE_BITS: usize = 6;
const DIALOG_GCD_K5_TAIL6_GRAPH_RAW_CODE_MASKS: [u16; DIALOG_GCD_K5_TAIL6_GRAPH_CODE_BITS] =
    [0x0d4, 0x0d1, 0x040, 0x05f, 0x00d, 0x081];
const DIALOG_GCD_K5_TAIL6_GRAPH_CODE_CONSTANT: u8 = 0x26;
const DIALOG_GCD_K5_TAIL6_GRAPH_RAW_CONSTANT: u16 = 0x1dc;
const DIALOG_GCD_K5_TAIL6_GRAPH_RAW_CODE_DECODE_MASKS: [u8; 9] =
    [0x00, 0x3a, 0x03, 0x13, 0x26, 0x00, 0x04, 0x20, 0x00];
const DIALOG_GCD_K5_TAIL6_GRAPH_SELECTOR_RAW_MASK: u16 = 0x085;
const DIALOG_GCD_K5_TAIL6_GRAPH_SELECTOR_ANF: &[u16] = &[0x00, 0x02, 0x04, 0x05, 0x32];
pub(crate) const DIALOG_GCD_K5_TAIL6_GRAPH_SUPPORT: [u32; 32] = [
    0x24924, 0x24925, 0x24928, 0x24929, 0x2492b, 0x2492c, 0x2492d, 0x2492f,
    0x24944, 0x24945, 0x24947, 0x24948, 0x24949, 0x2494b, 0x2494d, 0x2494f,
    0x24958, 0x24959, 0x2495b, 0x2495c, 0x2495d, 0x2495f, 0x24965, 0x24967,
    0x24968, 0x24969, 0x2496b, 0x24978, 0x24979, 0x2497b, 0x2497d, 0x2497f,
];

const DIALOG_GCD_K5_TAIL6_GRAPH9_STORED_STEPS: usize = 4;
const DIALOG_GCD_K5_TAIL6_GRAPH9_CODE_BITS: usize = 9;
const DIALOG_GCD_K5_TAIL6_GRAPH9_RAW_CODE_MASKS: [u16; DIALOG_GCD_K5_TAIL6_GRAPH9_CODE_BITS] =
    [0x9dc, 0xb6a, 0x717, 0x404, 0xe92, 0x00c, 0xa17, 0x7af, 0xf44];
const DIALOG_GCD_K5_TAIL6_GRAPH9_RAW_CONSTANT: u16 = 0xc6c;
const DIALOG_GCD_K5_TAIL6_GRAPH9_RAW_CODE_DECODE_MASKS: [u16; 12] = [
    0x058, 0x000, 0x131, 0x111, 0x18e, 0x01b, 0x0d2, 0x000, 0x17d, 0x0a7,
    0x139, 0x000,
];
const DIALOG_GCD_K5_TAIL6_GRAPH9_SELECTOR_RAW_MASK: u16 = 0x71e;
const DIALOG_GCD_K5_TAIL6_GRAPH9_SELECTOR_PIVOT: usize = 1;
const DIALOG_GCD_K5_TAIL6_GRAPH9_SELECTOR_ANF: &[u16] = &[
    0x000, 0x002, 0x003, 0x006, 0x008, 0x00a, 0x011, 0x012, 0x020, 0x021,
    0x024, 0x040, 0x041, 0x042, 0x048, 0x060, 0x080, 0x081, 0x088, 0x0a0,
    0x100, 0x02c, 0x034, 0x064, 0x0a2, 0x0a4,
];
pub(crate) const DIALOG_GCD_K5_TAIL6_GRAPH9_SUPPORT: [u32; 75] = [
    0x24924, 0x24925, 0x24928, 0x24929, 0x2492b, 0x2492c, 0x2492d, 0x2492f,
    0x24944, 0x24945, 0x24947, 0x24948, 0x24949, 0x2494b, 0x2494d, 0x2494f,
    0x24958, 0x24959, 0x2495b, 0x2495c, 0x2495d, 0x2495f, 0x24965, 0x24967,
    0x24968, 0x24969, 0x2496b, 0x24978, 0x24979, 0x2497b, 0x2497d, 0x2497f,
    0x24a27, 0x24a29, 0x24a2b, 0x24a2d, 0x24a2f, 0x24a38, 0x24a3f, 0x24a45,
    0x24a47, 0x24a49, 0x24a4b, 0x24a4d, 0x24a58, 0x24a59, 0x24a5b, 0x24a5f,
    0x24a65, 0x24a68, 0x24a6b, 0x24a78, 0x24a7c, 0x24ac5, 0x24ac8, 0x24ac9,
    0x24acb, 0x24acd, 0x24add, 0x24ae9, 0x24af8, 0x24af9, 0x24b29, 0x24b3c,
    0x24b3d, 0x24b45, 0x24b49, 0x24b4b, 0x24b5f, 0x24b79, 0x24bc5, 0x24bc9,
    0x24bd8, 0x24be4, 0x24bf9,
];

const DIALOG_GCD_K5_TAIL7_STORED_STEPS: usize = 4;
const DIALOG_GCD_K5_TAIL7_CODE_BITS: usize = 5;
const DIALOG_GCD_K5_TAIL7_PACKED_CODE_MASKS: [u32; DIALOG_GCD_K5_TAIL7_CODE_BITS] =
    [0x8a0, 0x204, 0x80011, 0x38, 0x100402];
const DIALOG_GCD_K5_TAIL7_RAW_CODE_MASKS: [u16; DIALOG_GCD_K5_TAIL7_CODE_BITS] =
    [0x0a20, 0x0140, 0x0009, 0x020c, 0x0082];
const DIALOG_GCD_K5_TAIL7_CODE_CONSTANT: u8 = 1 << 4;
pub(crate) const DIALOG_GCD_K5_TAIL7_SUPPORT: [u32; 20] = [
    0x124924, 0x124925, 0x124929, 0x12492b, 0x124928, 0x12492d, 0x12492f,
    0x12494b, 0x124947, 0x124945, 0x12492c, 0x124958, 0x124949, 0x12495b,
    0x124959, 0x124967, 0x12495d, 0x124a4b, 0x12497f, 0x124979,
];
const DIALOG_GCD_K5_TAIL7_RAW_ANF: [&[u16]; 12] = [
    &[1, 4, 7, 10, 12, 24, 28],
    &[0, 16],
    &[0, 7, 8, 10, 12, 24, 28],
    &[1, 7, 10, 12, 24, 28],
    &[1, 8, 9, 26],
    &[],
    &[11],
    &[],
    &[2, 11],
    &[0, 1],
    &[0, 11],
    &[0],
];

fn dialog_gcd_toggle_mcx_with_dirty(
    b: &mut B,
    controls: &[QubitId],
    dirty: &[QubitId],
    target: QubitId,
) {
    assert!(!controls.contains(&target));
    assert!(controls
        .iter()
        .enumerate()
        .all(|(index, q)| !controls[..index].contains(q)));
    match controls.len() {
        0 => b.x(target),
        1 => b.cx(controls[0], target),
        2 => b.ccx(controls[0], controls[1], target),
        count => {
            assert!(dirty.len() >= count - 2);
            let bridge = dirty[0];
            assert_ne!(bridge, target);
            assert!(!controls.contains(&bridge));
            dialog_gcd_toggle_mcx_with_dirty(
                b,
                &controls[..count - 1],
                &dirty[1..],
                bridge,
            );
            b.ccx(bridge, controls[count - 1], target);
            dialog_gcd_toggle_mcx_with_dirty(
                b,
                &controls[..count - 1],
                &dirty[1..],
                bridge,
            );
            b.ccx(bridge, controls[count - 1], target);
        }
    }
}

fn dialog_gcd_toggle_anf_with_dirty(
    b: &mut B,
    code: &[QubitId],
    target: QubitId,
    dirty_pool: &[QubitId],
    terms: &[u16],
) {
    assert!(code.len() <= u16::BITS as usize);
    assert!(!code.contains(&target));
    for &mask in terms {
        let controls = code
            .iter()
            .enumerate()
            .filter_map(|(index, &q)| ((mask >> index) & 1 != 0).then_some(q))
            .collect::<Vec<_>>();
        let dirty = dirty_pool
            .iter()
            .copied()
            .filter(|q| *q != target && !controls.contains(q))
            .collect::<Vec<_>>();
        dialog_gcd_toggle_mcx_with_dirty(b, &controls, &dirty, target);
    }
}

fn dialog_gcd_k5_tail6_graph9_toggle_code_from_raw(
    b: &mut B,
    code: &[QubitId],
    raw_block: &[QubitId],
) {
    assert_eq!(code.len(), DIALOG_GCD_K5_TAIL6_GRAPH9_CODE_BITS);
    assert_eq!(raw_block.len(), 15);
    for (code_index, &mask) in DIALOG_GCD_K5_TAIL6_GRAPH9_RAW_CODE_MASKS
        .iter()
        .enumerate()
    {
        for raw_bit in 0..12 {
            if (mask >> raw_bit) & 1 != 0 {
                b.cx(raw_block[raw_bit], code[code_index]);
            }
        }
    }
}

fn dialog_gcd_k5_tail6_graph9_toggle_linear_raw_from_code(
    b: &mut B,
    code: &[QubitId],
    raw_block: &[QubitId],
) {
    assert_eq!(code.len(), DIALOG_GCD_K5_TAIL6_GRAPH9_CODE_BITS);
    assert_eq!(raw_block.len(), 15);
    for (raw_index, &mask) in DIALOG_GCD_K5_TAIL6_GRAPH9_RAW_CODE_DECODE_MASKS
        .iter()
        .enumerate()
    {
        if (DIALOG_GCD_K5_TAIL6_GRAPH9_RAW_CONSTANT >> raw_index) & 1 != 0 {
            b.x(raw_block[raw_index]);
        }
        for code_bit in 0..DIALOG_GCD_K5_TAIL6_GRAPH9_CODE_BITS {
            if (mask >> code_bit) & 1 != 0 {
                b.cx(code[code_bit], raw_block[raw_index]);
            }
        }
    }
}

fn dialog_gcd_k5_tail6_graph9_toggle_selector_fanout(
    b: &mut B,
    raw_block: &[QubitId],
) {
    assert_eq!(raw_block.len(), 15);
    assert_ne!(
        (DIALOG_GCD_K5_TAIL6_GRAPH9_SELECTOR_RAW_MASK
            >> DIALOG_GCD_K5_TAIL6_GRAPH9_SELECTOR_PIVOT)
            & 1,
        0
    );
    let pivot = raw_block[DIALOG_GCD_K5_TAIL6_GRAPH9_SELECTOR_PIVOT];
    for raw_index in 0..12 {
        if raw_index != DIALOG_GCD_K5_TAIL6_GRAPH9_SELECTOR_PIVOT
            && (DIALOG_GCD_K5_TAIL6_GRAPH9_SELECTOR_RAW_MASK >> raw_index) & 1 != 0
        {
            b.cx(pivot, raw_block[raw_index]);
        }
    }
}

fn dialog_gcd_k5_tail6_graph9_toggle_selector(
    b: &mut B,
    code: &[QubitId],
    raw_block: &[QubitId],
) {
    dialog_gcd_toggle_anf_with_dirty(
        b,
        code,
        raw_block[DIALOG_GCD_K5_TAIL6_GRAPH9_SELECTOR_PIVOT],
        raw_block,
        DIALOG_GCD_K5_TAIL6_GRAPH9_SELECTOR_ANF,
    );
}

fn dialog_gcd_k5_tail6_graph9_compress_raw_to_block(
    b: &mut B,
    code: &[QubitId],
    raw_block: &[QubitId],
) {
    dialog_gcd_k5_tail6_graph9_toggle_code_from_raw(b, code, raw_block);
    dialog_gcd_k5_tail6_graph9_toggle_linear_raw_from_code(b, code, raw_block);
    dialog_gcd_k5_tail6_graph9_toggle_selector_fanout(b, raw_block);
    dialog_gcd_k5_tail6_graph9_toggle_selector(b, code, raw_block);
}

fn dialog_gcd_k5_tail6_graph9_decompress_block_to_raw(
    b: &mut B,
    code: &[QubitId],
    raw_block: &[QubitId],
) {
    dialog_gcd_k5_tail6_graph9_toggle_selector(b, code, raw_block);
    dialog_gcd_k5_tail6_graph9_toggle_selector_fanout(b, raw_block);
    dialog_gcd_k5_tail6_graph9_toggle_linear_raw_from_code(b, code, raw_block);
    dialog_gcd_k5_tail6_graph9_toggle_code_from_raw(b, code, raw_block);
}

fn dialog_gcd_k5_tail6_graph9_raw_word(pattern: u32) -> u16 {
    let mut raw_word = 0u16;
    for slot in 0..DIALOG_GCD_K5_TAIL6_GRAPH9_STORED_STEPS {
        let digit = ((pattern >> (3 * slot)) & 7) as u16;
        raw_word |= (digit & 1) << (2 * slot);
        raw_word |= ((digit >> 1) & 1) << (2 * slot + 1);
        raw_word |= ((digit >> 2) & 1)
            << (2 * DIALOG_GCD_K5_TAIL6_GRAPH9_STORED_STEPS + slot);
    }
    raw_word
}

fn dialog_gcd_k5_tail6_graph9_code_word(raw_word: u16) -> u16 {
    DIALOG_GCD_K5_TAIL6_GRAPH9_RAW_CODE_MASKS
        .iter()
        .enumerate()
        .fold(0u16, |code, (index, &mask)| {
            code ^ ((((raw_word & mask).count_ones() & 1) as u16) << index)
        })
}

fn dialog_gcd_k5_tail6_graph9_selector_word(code: u16) -> u16 {
    DIALOG_GCD_K5_TAIL6_GRAPH9_SELECTOR_ANF
        .iter()
        .fold(0u16, |selector, &term| {
            selector ^ u16::from(code & term == term)
        })
}

fn dialog_gcd_k5_tail6_graph9_decode_word(code: u16) -> u16 {
    let selector = dialog_gcd_k5_tail6_graph9_selector_word(code);
    DIALOG_GCD_K5_TAIL6_GRAPH9_RAW_CODE_DECODE_MASKS
        .iter()
        .enumerate()
        .fold(DIALOG_GCD_K5_TAIL6_GRAPH9_RAW_CONSTANT, |raw, (index, &mask)| {
            let bit = ((code & mask).count_ones() & 1) as u16
                ^ (selector
                    & ((DIALOG_GCD_K5_TAIL6_GRAPH9_SELECTOR_RAW_MASK >> index) & 1));
            raw ^ (bit << index)
        })
}

pub(crate) fn dialog_gcd_k5_tail6_graph9_supports(pattern: u32) -> bool {
    if pattern >> 12 != 0x24 {
        return false;
    }
    let raw = dialog_gcd_k5_tail6_graph9_raw_word(pattern);
    let code = dialog_gcd_k5_tail6_graph9_code_word(raw);
    dialog_gcd_k5_tail6_graph9_decode_word(code) == raw
}

pub(crate) fn dialog_gcd_k5_tail6_graph9_codec_selftest() -> Result<(), String> {
    use sha3::digest::{ExtendableOutput, Update};

    for batch_start in (0..DIALOG_GCD_K5_TAIL6_GRAPH9_SUPPORT.len()).step_by(64) {
        let patterns = &DIALOG_GCD_K5_TAIL6_GRAPH9_SUPPORT
            [batch_start..(batch_start + 64).min(DIALOG_GCD_K5_TAIL6_GRAPH9_SUPPORT.len())];
        let mut raw_masks = [0u64; 15];
        let mut code_masks = [0u64; DIALOG_GCD_K5_TAIL6_GRAPH9_CODE_BITS];
        for (shot, &pattern) in patterns.iter().enumerate() {
            if !dialog_gcd_k5_tail6_graph9_supports(pattern) {
                return Err(format!("support pattern 0x{pattern:x} fails graph relation"));
            }
            let shot_bit = 1u64 << shot;
            let raw_word = dialog_gcd_k5_tail6_graph9_raw_word(pattern);
            for raw_bit in 0..12 {
                if (raw_word >> raw_bit) & 1 != 0 {
                    raw_masks[raw_bit] |= shot_bit;
                }
            }
            let code = dialog_gcd_k5_tail6_graph9_code_word(raw_word);
            for (index, mask) in code_masks.iter_mut().enumerate() {
                if (code >> index) & 1 != 0 {
                    *mask |= shot_bit;
                }
            }
        }

        let build_codec = |decompress: bool| {
            let mut b = B::new();
            let code = b.alloc_qubits(DIALOG_GCD_K5_TAIL6_GRAPH9_CODE_BITS);
            let raw = b.alloc_qubits(15);
            if decompress {
                dialog_gcd_k5_tail6_graph9_decompress_block_to_raw(&mut b, &code, &raw);
            } else {
                dialog_gcd_k5_tail6_graph9_compress_raw_to_block(&mut b, &code, &raw);
            }
            (b.ops, code, raw, b.next_qubit as usize, b.next_bit as usize)
        };

        let run = |decompress: bool| {
            let (ops, code, raw, num_qubits, num_bits) = build_codec(decompress);
            let mut seed = sha3::Shake128::default();
            seed.update(b"dialog-gcd-k5-tail6-graph9-codec-selftest");
            seed.update(&(batch_start as u64).to_le_bytes());
            let mut xof = seed.finalize_xof();
            let mut sim = Simulator::new(num_qubits, num_bits, &mut xof);
            sim.clear_for_shot();
            let source = if decompress { &code_masks[..] } else { &raw_masks[..] };
            let targets = if decompress { &code[..] } else { &raw[..] };
            for (&qubit, &mask) in targets.iter().zip(source.iter()) {
                *sim.qubit_mut(qubit) = mask;
            }
            sim.apply_iter(ops.iter());
            (
                code.iter().map(|&q| sim.qubit(q)).collect::<Vec<_>>(),
                raw.iter().map(|&q| sim.qubit(q)).collect::<Vec<_>>(),
                sim.phase,
            )
        };

        let active_mask = if patterns.len() == 64 {
            u64::MAX
        } else {
            (1u64 << patterns.len()) - 1
        };
        let (forward_code, forward_raw, forward_phase) = run(false);
        if forward_phase & active_mask != 0 {
            return Err(format!(
                "forward phase garbage in batch {batch_start}: 0x{:x}",
                forward_phase & active_mask
            ));
        }
        if forward_code
            .iter()
            .zip(code_masks.iter())
            .any(|(&got, &want)| (got ^ want) & active_mask != 0)
        {
            return Err(format!(
                "forward code mismatch in batch {batch_start}: got {forward_code:x?}, want {code_masks:x?}"
            ));
        }
        if forward_raw
            .iter()
            .any(|&mask| mask & active_mask != 0)
        {
            return Err(format!(
                "forward raw garbage in batch {batch_start}: {forward_raw:x?}"
            ));
        }

        let (reverse_code, reverse_raw, reverse_phase) = run(true);
        if reverse_phase & active_mask != 0 {
            return Err(format!(
                "reverse phase garbage in batch {batch_start}: 0x{:x}",
                reverse_phase & active_mask
            ));
        }
        if reverse_code
            .iter()
            .any(|&mask| mask & active_mask != 0)
        {
            return Err(format!(
                "reverse code garbage in batch {batch_start}: {reverse_code:x?}"
            ));
        }
        if reverse_raw
            .iter()
            .zip(raw_masks.iter())
            .any(|(&got, &want)| (got ^ want) & active_mask != 0)
        {
            return Err(format!(
                "reverse raw mismatch in batch {batch_start}: got {reverse_raw:x?}, want {raw_masks:x?}"
            ));
        }
    }
    Ok(())
}

fn dialog_gcd_k5_tail6_graph_toggle_code_from_raw(
    b: &mut B,
    code: &[QubitId],
    raw_block: &[QubitId],
) {
    assert_eq!(code.len(), DIALOG_GCD_K5_TAIL6_GRAPH_CODE_BITS);
    assert_eq!(raw_block.len(), 15);
    for (code_index, &mask) in DIALOG_GCD_K5_TAIL6_GRAPH_RAW_CODE_MASKS
        .iter()
        .enumerate()
    {
        for raw_bit in 0..9 {
            if (mask >> raw_bit) & 1 != 0 {
                b.cx(raw_block[raw_bit], code[code_index]);
            }
        }
        if (DIALOG_GCD_K5_TAIL6_GRAPH_CODE_CONSTANT >> code_index) & 1 != 0 {
            b.x(code[code_index]);
        }
    }
}

fn dialog_gcd_k5_tail6_graph_toggle_linear_raw_from_code(
    b: &mut B,
    code: &[QubitId],
    raw_block: &[QubitId],
) {
    assert_eq!(code.len(), DIALOG_GCD_K5_TAIL6_GRAPH_CODE_BITS);
    assert_eq!(raw_block.len(), 15);
    for (raw_index, &mask) in DIALOG_GCD_K5_TAIL6_GRAPH_RAW_CODE_DECODE_MASKS
        .iter()
        .enumerate()
    {
        if (DIALOG_GCD_K5_TAIL6_GRAPH_RAW_CONSTANT >> raw_index) & 1 != 0 {
            b.x(raw_block[raw_index]);
        }
        for code_bit in 0..DIALOG_GCD_K5_TAIL6_GRAPH_CODE_BITS {
            if (mask >> code_bit) & 1 != 0 {
                b.cx(code[code_bit], raw_block[raw_index]);
            }
        }
    }
}

fn dialog_gcd_k5_tail6_graph_toggle_selector_fanout(
    b: &mut B,
    raw_block: &[QubitId],
) {
    assert_eq!(raw_block.len(), 15);
    let pivot = raw_block[0];
    assert_eq!(DIALOG_GCD_K5_TAIL6_GRAPH_SELECTOR_RAW_MASK & 1, 1);
    for raw_index in 1..9 {
        if (DIALOG_GCD_K5_TAIL6_GRAPH_SELECTOR_RAW_MASK >> raw_index) & 1 != 0 {
            b.cx(pivot, raw_block[raw_index]);
        }
    }
}

fn dialog_gcd_k5_tail6_graph_toggle_selector(
    b: &mut B,
    code: &[QubitId],
    raw_block: &[QubitId],
) {
    dialog_gcd_toggle_anf_with_dirty(
        b,
        code,
        raw_block[0],
        raw_block,
        DIALOG_GCD_K5_TAIL6_GRAPH_SELECTOR_ANF,
    );
}

fn dialog_gcd_k5_tail6_graph_compress_raw_to_block(
    b: &mut B,
    code: &[QubitId],
    raw_block: &[QubitId],
) {
    dialog_gcd_k5_tail6_graph_toggle_code_from_raw(b, code, raw_block);
    dialog_gcd_k5_tail6_graph_toggle_linear_raw_from_code(b, code, raw_block);
    dialog_gcd_k5_tail6_graph_toggle_selector_fanout(b, raw_block);
    dialog_gcd_k5_tail6_graph_toggle_selector(b, code, raw_block);
}

fn dialog_gcd_k5_tail6_graph_decompress_block_to_raw(
    b: &mut B,
    code: &[QubitId],
    raw_block: &[QubitId],
) {
    dialog_gcd_k5_tail6_graph_toggle_selector(b, code, raw_block);
    dialog_gcd_k5_tail6_graph_toggle_selector_fanout(b, raw_block);
    dialog_gcd_k5_tail6_graph_toggle_linear_raw_from_code(b, code, raw_block);
    dialog_gcd_k5_tail6_graph_toggle_code_from_raw(b, code, raw_block);
}

pub(crate) fn dialog_gcd_k5_tail6_graph_codec_selftest() -> Result<(), String> {
    use sha3::digest::{ExtendableOutput, Update};

    let mut raw_masks = [0u64; 15];
    let mut code_masks = [0u64; DIALOG_GCD_K5_TAIL6_GRAPH_CODE_BITS];
    for shot in 0..64 {
        let pattern =
            DIALOG_GCD_K5_TAIL6_GRAPH_SUPPORT[shot % DIALOG_GCD_K5_TAIL6_GRAPH_SUPPORT.len()];
        let shot_bit = 1u64 << shot;
        let mut raw_word = 0u16;
        for slot in 0..DIALOG_GCD_K5_TAIL6_GRAPH_STORED_STEPS {
            if (pattern >> (3 * slot)) & 1 != 0 {
                raw_masks[2 * slot] |= shot_bit;
                raw_word |= 1 << (2 * slot);
            }
            if (pattern >> (3 * slot + 1)) & 1 != 0 {
                raw_masks[2 * slot + 1] |= shot_bit;
                raw_word |= 1 << (2 * slot + 1);
            }
            if (pattern >> (3 * slot + 2)) & 1 != 0 {
                let raw_index = 2 * DIALOG_GCD_K5_TAIL6_GRAPH_STORED_STEPS + slot;
                raw_masks[raw_index] |= shot_bit;
                raw_word |= 1 << raw_index;
            }
        }
        let code = DIALOG_GCD_K5_TAIL6_GRAPH_RAW_CODE_MASKS
            .iter()
            .enumerate()
            .fold(DIALOG_GCD_K5_TAIL6_GRAPH_CODE_CONSTANT, |packed, (index, &mask)| {
                packed ^ ((((raw_word & mask).count_ones() & 1) as u8) << index)
            });
        for (index, mask) in code_masks.iter_mut().enumerate() {
            if (code >> index) & 1 != 0 {
                *mask |= shot_bit;
            }
        }
    }

    let build_codec = |decompress: bool| {
        let mut b = B::new();
        let code = b.alloc_qubits(DIALOG_GCD_K5_TAIL6_GRAPH_CODE_BITS);
        let raw = b.alloc_qubits(15);
        if decompress {
            dialog_gcd_k5_tail6_graph_decompress_block_to_raw(&mut b, &code, &raw);
        } else {
            dialog_gcd_k5_tail6_graph_compress_raw_to_block(&mut b, &code, &raw);
        }
        (b.ops, code, raw, b.next_qubit as usize, b.next_bit as usize)
    };

    let run = |decompress: bool| {
        let (ops, code, raw, num_qubits, num_bits) = build_codec(decompress);
        let mut seed = sha3::Shake128::default();
        seed.update(b"dialog-gcd-k5-tail6-graph-codec-selftest");
        let mut xof = seed.finalize_xof();
        let mut sim = Simulator::new(num_qubits, num_bits, &mut xof);
        sim.clear_for_shot();
        let source = if decompress { &code_masks[..] } else { &raw_masks[..] };
        let targets = if decompress { &code[..] } else { &raw[..] };
        for (&qubit, &mask) in targets.iter().zip(source.iter()) {
            *sim.qubit_mut(qubit) = mask;
        }
        sim.apply_iter(ops.iter());
        (
            code.iter().map(|&q| sim.qubit(q)).collect::<Vec<_>>(),
            raw.iter().map(|&q| sim.qubit(q)).collect::<Vec<_>>(),
            sim.phase,
        )
    };

    let (forward_code, forward_raw, forward_phase) = run(false);
    if forward_phase != 0 {
        return Err(format!("forward phase garbage 0x{forward_phase:x}"));
    }
    if forward_code != code_masks {
        return Err(format!(
            "forward code mismatch: got {forward_code:x?}, want {code_masks:x?}"
        ));
    }
    if forward_raw.iter().any(|&mask| mask != 0) {
        return Err(format!("forward raw garbage: {forward_raw:x?}"));
    }

    let (reverse_code, reverse_raw, reverse_phase) = run(true);
    if reverse_phase != 0 {
        return Err(format!("reverse phase garbage 0x{reverse_phase:x}"));
    }
    if reverse_code.iter().any(|&mask| mask != 0) {
        return Err(format!("reverse code garbage: {reverse_code:x?}"));
    }
    if reverse_raw != raw_masks {
        return Err(format!(
            "reverse raw mismatch: got {reverse_raw:x?}, want {raw_masks:x?}"
        ));
    }
    Ok(())
}

fn dialog_gcd_k5_tail7_toggle_code_from_raw(
    b: &mut B,
    code: &[QubitId],
    raw_block: &[QubitId],
) {
    assert_eq!(code.len(), DIALOG_GCD_K5_TAIL7_CODE_BITS);
    assert_eq!(raw_block.len(), 15);
    for (code_index, &mask) in DIALOG_GCD_K5_TAIL7_RAW_CODE_MASKS.iter().enumerate() {
        for raw_bit in 0..12 {
            if (mask >> raw_bit) & 1 != 0 {
                b.cx(raw_block[raw_bit], code[code_index]);
            }
        }
        if (DIALOG_GCD_K5_TAIL7_CODE_CONSTANT >> code_index) & 1 != 0 {
            b.x(code[code_index]);
        }
    }
}

fn dialog_gcd_k5_tail7_toggle_raw_from_code(
    b: &mut B,
    code: &[QubitId],
    raw_block: &[QubitId],
) {
    assert_eq!(code.len(), DIALOG_GCD_K5_TAIL7_CODE_BITS);
    assert_eq!(raw_block.len(), 15);
    for (raw_index, terms) in DIALOG_GCD_K5_TAIL7_RAW_ANF.iter().enumerate() {
        dialog_gcd_toggle_anf_with_dirty(b, code, raw_block[raw_index], raw_block, terms);
    }
}

fn dialog_gcd_k5_tail7_compress_raw_to_block(
    b: &mut B,
    code: &[QubitId],
    raw_block: &[QubitId],
) {
    dialog_gcd_k5_tail7_toggle_code_from_raw(b, code, raw_block);
    dialog_gcd_k5_tail7_toggle_raw_from_code(b, code, raw_block);
}

fn dialog_gcd_k5_tail7_decompress_block_to_raw(
    b: &mut B,
    code: &[QubitId],
    raw_block: &[QubitId],
) {
    dialog_gcd_k5_tail7_toggle_raw_from_code(b, code, raw_block);
    dialog_gcd_k5_tail7_toggle_code_from_raw(b, code, raw_block);
}

pub(crate) fn dialog_gcd_k5_tail7_codec_selftest() -> Result<(), String> {
    use sha3::digest::{ExtendableOutput, Update};

    let mut raw_masks = [0u64; 15];
    let mut code_masks = [0u64; DIALOG_GCD_K5_TAIL7_CODE_BITS];
    for shot in 0..64 {
        let pattern = DIALOG_GCD_K5_TAIL7_SUPPORT[shot % DIALOG_GCD_K5_TAIL7_SUPPORT.len()];
        let shot_bit = 1u64 << shot;
        for slot in 0..DIALOG_GCD_K5_TAIL7_STORED_STEPS {
            if (pattern >> (3 * slot)) & 1 != 0 {
                raw_masks[2 * slot] |= shot_bit;
            }
            if (pattern >> (3 * slot + 1)) & 1 != 0 {
                raw_masks[2 * slot + 1] |= shot_bit;
            }
            if (pattern >> (3 * slot + 2)) & 1 != 0 {
                raw_masks[2 * DIALOG_GCD_K5_TAIL7_STORED_STEPS + slot] |= shot_bit;
            }
        }
        let code = DIALOG_GCD_K5_TAIL7_PACKED_CODE_MASKS
            .iter()
            .enumerate()
            .fold(0u8, |packed, (index, &mask)| {
                packed | ((((pattern & mask).count_ones() & 1) as u8) << index)
            });
        for (index, mask) in code_masks.iter_mut().enumerate() {
            if (code >> index) & 1 != 0 {
                *mask |= shot_bit;
            }
        }
    }

    let build_codec = |decompress: bool| {
        let mut b = B::new();
        let code = b.alloc_qubits(DIALOG_GCD_K5_TAIL7_CODE_BITS);
        let raw = b.alloc_qubits(15);
        if decompress {
            dialog_gcd_k5_tail7_decompress_block_to_raw(&mut b, &code, &raw);
        } else {
            dialog_gcd_k5_tail7_compress_raw_to_block(&mut b, &code, &raw);
        }
        (b.ops, code, raw, b.next_qubit as usize, b.next_bit as usize)
    };

    let run = |decompress: bool| {
        let (ops, code, raw, num_qubits, num_bits) = build_codec(decompress);
        let mut seed = sha3::Shake128::default();
        seed.update(b"dialog-gcd-k5-tail7-codec-selftest");
        let mut xof = seed.finalize_xof();
        let mut sim = Simulator::new(num_qubits, num_bits, &mut xof);
        sim.clear_for_shot();
        let source = if decompress { &code_masks[..] } else { &raw_masks[..] };
        let targets = if decompress { &code[..] } else { &raw[..] };
        for (&qubit, &mask) in targets.iter().zip(source.iter()) {
            *sim.qubit_mut(qubit) = mask;
        }
        sim.apply_iter(ops.iter());
        (
            code.iter().map(|&q| sim.qubit(q)).collect::<Vec<_>>(),
            raw.iter().map(|&q| sim.qubit(q)).collect::<Vec<_>>(),
            sim.phase,
        )
    };

    let (forward_code, forward_raw, forward_phase) = run(false);
    if forward_phase != 0 {
        return Err(format!("forward phase garbage 0x{forward_phase:x}"));
    }
    if forward_code != code_masks {
        return Err(format!(
            "forward code mismatch: got {forward_code:x?}, want {code_masks:x?}"
        ));
    }
    if forward_raw.iter().any(|&mask| mask != 0) {
        return Err(format!("forward raw garbage: {forward_raw:x?}"));
    }

    let (reverse_code, reverse_raw, reverse_phase) = run(true);
    if reverse_phase != 0 {
        return Err(format!("reverse phase garbage 0x{reverse_phase:x}"));
    }
    if reverse_code.iter().any(|&mask| mask != 0) {
        return Err(format!("reverse code garbage: {reverse_code:x?}"));
    }
    if reverse_raw != raw_masks {
        return Err(format!(
            "reverse raw mismatch: got {reverse_raw:x?}, want {raw_masks:x?}"
        ));
    }
    Ok(())
}

fn dialog_gcd_k5_tail_pair1_compress_raw_to_block(
    b: &mut B,
    compressed_block: &[QubitId],
    raw_block: &[QubitId],
    swap_host: bool,
) {
    assert_eq!(compressed_block.len(), 1);
    assert_eq!(raw_block.len(), 15);
    // Supported tail language:
    //   step 0 = (b0, b0_and_b1, s2) in {(0,0,1), (1,0,1)}
    //   step 1 = (0,0,1)
    // The sole code bit is step-0 b0.
    if swap_host {
        b.swap(compressed_block[0], raw_block[0]);
    } else {
        b.cx(compressed_block[0], raw_block[0]);
    }
    b.x(dialog_gcd_raw_s2(raw_block, 0));
    b.x(dialog_gcd_raw_s2(raw_block, 1));
}

fn dialog_gcd_k5_tail_pair1_decompress_block_to_raw(
    b: &mut B,
    compressed_block: &[QubitId],
    raw_block: &[QubitId],
    swap_host: bool,
) {
    assert_eq!(compressed_block.len(), 1);
    assert_eq!(raw_block.len(), 15);
    b.x(dialog_gcd_raw_s2(raw_block, 1));
    b.x(dialog_gcd_raw_s2(raw_block, 0));
    if swap_host {
        b.swap(compressed_block[0], raw_block[0]);
    } else {
        b.cx(compressed_block[0], raw_block[0]);
    }
}

pub(crate) fn emit_dialog_gcd_round763_compressed_block_swapper(
    b: &mut B,
    pair: &[QubitId],
    compressed_block: &[QubitId],
    scratch: QubitId,
    slot: usize,
) {
    assert_eq!(pair.len(), 2);
    assert_eq!(compressed_block.len(), 5);
    assert!(slot < 3);
    let mut block = compressed_block.to_vec();
    block.push(scratch);
    emit_dialog_gcd_round763_compressor_inverse(b, &block);
    b.swap(pair[0], block[2 * slot]);
    b.swap(pair[1], block[2 * slot + 1]);
    emit_dialog_gcd_round763_compressor(b, &block);
}

pub(crate) fn dialog_gcd_compressed_sidecar_blocks() -> usize {
    let group_size = dialog_gcd_sidecar_group_size();
    let blocks = (dialog_gcd_active_iterations() + group_size - 1) / group_size;
    if dialog_gcd_k5_tail7_enabled()
        || dialog_gcd_k5_tail6_graph_enabled()
        || dialog_gcd_k5_tail6_graph9_enabled()
    {
        blocks - 1
    } else {
        blocks
    }
}

fn dialog_gcd_compressed_sidecar_block_index(step: usize) -> usize {
    if dialog_gcd_k5_tail7_enabled()
        && step >= dialog_gcd_active_iterations() - 7
        || dialog_gcd_k5_tail6_graph_enabled()
            && step >= dialog_gcd_active_iterations() - 6
        || dialog_gcd_k5_tail6_graph9_enabled()
            && step >= dialog_gcd_active_iterations() - 6
    {
        dialog_gcd_compressed_sidecar_blocks() - 1
    } else {
        step / dialog_gcd_sidecar_group_size()
    }
}

fn dialog_gcd_compressed_sidecar_block_bits(block: usize) -> usize {
    if dialog_gcd_k5_head11_enabled() && block == 0 {
        DIALOG_GCD_K5_HEAD11_DATA_WIRES.len()
    } else if dialog_gcd_k5_tail6_graph9_enabled()
        && block + 1 == dialog_gcd_compressed_sidecar_blocks()
    {
        DIALOG_GCD_K5_TAIL6_GRAPH9_CODE_BITS
    } else if dialog_gcd_k5_tail6_graph_enabled()
        && block + 1 == dialog_gcd_compressed_sidecar_blocks()
    {
        DIALOG_GCD_K5_TAIL6_GRAPH_CODE_BITS
    } else if dialog_gcd_k5_tail7_enabled()
        && block + 1 == dialog_gcd_compressed_sidecar_blocks()
    {
        DIALOG_GCD_K5_TAIL7_CODE_BITS
    } else if dialog_gcd_k5_tail_pair1_enabled()
        && block + 1 == dialog_gcd_compressed_sidecar_blocks()
    {
        1
    } else {
        dialog_gcd_block_bits()
    }
}

fn dialog_gcd_compressed_sidecar_block_offset(block: usize) -> usize {
    (0..block)
        .map(dialog_gcd_compressed_sidecar_block_bits)
        .sum()
}

pub(crate) fn dialog_gcd_compressed_sidecar_bits() -> usize {
    (0..dialog_gcd_compressed_sidecar_blocks())
        .map(dialog_gcd_compressed_sidecar_block_bits)
        .sum()
}

pub(crate) fn dialog_gcd_compressed_sidecar_block(compressed_log: &[QubitId], step: usize) -> &[QubitId] {
    let block = dialog_gcd_compressed_sidecar_block_index(step);
    let start = dialog_gcd_compressed_sidecar_block_offset(block);
    let bits = dialog_gcd_compressed_sidecar_block_bits(block);
    &compressed_log[start..start + bits]
}

pub(crate) fn dialog_gcd_compressed_log_u_high_runway_enabled() -> bool {
    // Prototype, deliberately NOT enabled by configure_ecdsafail_submission_route.
    //
    // The wrapper used to allocate all of u and the complete compressed
    // transcript at once.  Instead, a late transcript suffix can use high u
    // lanes: those cells are not touched until forward replay has shrunk u below
    // their hosts, stay live across terminal-reuse apply, and are consumed by
    // reverse replay before u grows back into them.
    //
    // This is an experimental support-envelope optimization: it relies on the
    // same terminal convergence and width envelope as terminal reuse and
    // variable-width tobitvector.  Default OFF keeps the accepted route
    // byte-identical.
    // K=2: runway layout is now block_bits()-aware (8-bit stride), so it is safe
    // to host the wider K2 transcript blocks on u-high — this is the peak lever.
    std::env::var("DIALOG_GCD_COMPRESSED_LOG_U_HIGH_RUNWAY")
        .ok()
        .as_deref()
        == Some("1")
}

fn dialog_gcd_k5_constant_tail_stored_steps(block_steps: usize) -> Option<usize> {
    if dialog_gcd_k5_tail6_graph9_enabled() && block_steps == 6 {
        Some(DIALOG_GCD_K5_TAIL6_GRAPH9_STORED_STEPS)
    } else if dialog_gcd_k5_tail6_graph_enabled() && block_steps == 6 {
        Some(DIALOG_GCD_K5_TAIL6_GRAPH_STORED_STEPS)
    } else if dialog_gcd_k5_tail7_enabled() && block_steps == 7 {
        Some(DIALOG_GCD_K5_TAIL7_STORED_STEPS)
    } else {
        None
    }
}

fn dialog_gcd_k5_fixed_tail_apply_enabled() -> bool {
    (dialog_gcd_k5_tail7_enabled()
        || dialog_gcd_k5_tail6_graph_enabled()
        || dialog_gcd_k5_tail6_graph9_enabled())
        && (std::env::var("DIALOG_GCD_K5_FIXED_TAIL_APPLY")
            .ok()
            .as_deref()
            == Some("1")
            || std::env::var("DIALOG_GCD_K5_TAIL7_UNCONDITIONAL_APPLY")
                .ok()
                .as_deref()
                == Some("1"))
}

pub(crate) fn dialog_gcd_compressed_log_u_high_runway_blocks() -> usize {
    // Optional tuning cap for the prototype.  The uncapped layout parks the
    // longest suffix; lowering the cap is useful when balancing wrapper savings
    // against reverse-replay scratch pressure.  On the accepted a8d8d5a route,
    // 16 whole blocks is the largest prefix-independent tail runway before the
    // reverse add loses its cheap scratch host.  Keep larger schedules available
    // as an explicit experiment, but default the opt-in prototype to that safe
    // subset.
    std::env::var("DIALOG_GCD_COMPRESSED_LOG_U_HIGH_RUNWAY_BLOCKS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(16)
}

#[derive(Clone, Debug)]
pub(crate) struct DialogGcdCompressedLogUHighRunway {
    remapped_log: Vec<QubitId>,
    parked_u_indices: Vec<usize>,
}

pub(crate) fn dialog_gcd_slice_intersects(a: &[QubitId], b: &[QubitId]) -> bool {
    a.iter().any(|q| b.contains(q))
}

pub(crate) fn dialog_gcd_runway_layout() -> Vec<(usize, usize)> {
    // Leave the top six u lanes unparked.  The accepted a8d8d5a route hosts a
    // raw 3-step block there whenever the tail is wide enough; reserving those
    // lanes keeps that scratch host disjoint from parked transcript cells.
    let raw_block_bits = 2 * DIALOG_GCD_HIGH_TAIL_ALIAS_GROUP_SIZE;
    let Some(highest_host) = N.checked_sub(raw_block_bits + 1) else {
        return Vec::new();
    };
    let blocks = dialog_gcd_compressed_sidecar_blocks();

    // Find the longest whole-block suffix that fits.  Blocks are assigned in
    // forward order to descending u positions: the earliest parked block gets
    // the highest hosts because it is replayed last and therefore needs the
    // widest inactive-u threshold.
    let first_allowed = blocks.saturating_sub(dialog_gcd_compressed_log_u_high_runway_blocks());
    for first_block in first_allowed..blocks {
        let mut next_host = highest_host;
        let mut layout = Vec::with_capacity(
            (first_block..blocks)
                .map(dialog_gcd_compressed_sidecar_block_bits)
                .sum(),
        );
        let mut fits = true;
        for block in first_block..blocks {
            let (start, end) = dialog_gcd_compressed_sidecar_block_step_range(block);
            let active_threshold = (start..end)
                .map(dialog_gcd_tobitvector_active_width)
                .max()
                .unwrap_or(1);
            let block_offset = dialog_gcd_compressed_sidecar_block_offset(block);
            for slot in 0..dialog_gcd_compressed_sidecar_block_bits(block) {
                if next_host < active_threshold {
                    fits = false;
                    break;
                }
                layout.push((block_offset + slot, next_host));
                let Some(next) = next_host.checked_sub(1) else {
                    fits = false;
                    break;
                };
                next_host = next;
            }
            if !fits {
                break;
            }
        }
        if fits {
            return layout;
        }
    }
    Vec::new()
}

pub(crate) fn dialog_gcd_allocated_compressed_sidecar_bits() -> usize {
    if dialog_gcd_compressed_log_u_high_runway_enabled() {
        dialog_gcd_compressed_sidecar_bits() - dialog_gcd_runway_layout().len()
    } else {
        dialog_gcd_compressed_sidecar_bits()
    }
}

pub(crate) fn dialog_gcd_build_compressed_log_u_high_runway(
    u: &[QubitId],
    allocated_log: &[QubitId],
) -> Option<DialogGcdCompressedLogUHighRunway> {
    if !dialog_gcd_compressed_log_u_high_runway_enabled() {
        return None;
    }
    assert_eq!(u.len(), N);
    let layout = dialog_gcd_runway_layout();
    if layout.is_empty() {
        return None;
    }

    let expected_allocated = dialog_gcd_compressed_sidecar_bits() - layout.len();
    assert_eq!(allocated_log.len(), expected_allocated);
    let first_relocated = layout[0].0;
    assert_eq!(first_relocated, allocated_log.len());
    let mut remapped_log = allocated_log.to_vec();
    let mut parked_u_indices = Vec::with_capacity(layout.len());
    for (log_index, u_index) in layout {
        // These logical transcript cells are not needed until their late
        // forward blocks, when the width envelope guarantees that u[u_index] is
        // inactive and |0>.  Reverse consumes them before u grows back into the
        // same hosts.
        assert_eq!(log_index, remapped_log.len());
        remapped_log.push(u[u_index]);
        parked_u_indices.push(u_index);
    }
    assert_eq!(remapped_log.len(), dialog_gcd_compressed_sidecar_bits());
    Some(DialogGcdCompressedLogUHighRunway {
        remapped_log,
        parked_u_indices,
    })
}

pub(crate) fn dialog_gcd_release_terminal_u(
    b: &mut B,
    u: &[QubitId],
    runway: Option<&DialogGcdCompressedLogUHighRunway>,
) {
    for (index, &q) in u.iter().enumerate() {
        if runway.is_none_or(|r| !r.parked_u_indices.contains(&index)) {
            b.free(q);
        }
    }
}

pub(crate) fn dialog_gcd_reacquire_terminal_u(
    b: &mut B,
    u: &[QubitId],
    runway: Option<&DialogGcdCompressedLogUHighRunway>,
) {
    for (index, &q) in u.iter().enumerate() {
        if runway.is_none_or(|r| !r.parked_u_indices.contains(&index)) {
            b.reacquire(q);
        }
    }
}

pub(crate) fn dialog_gcd_runway_safe_future_prefix<'a>(
    future: Option<&'a [QubitId]>,
    u: &[QubitId],
    active_width: usize,
) -> Option<&'a [QubitId]> {
    let active_u = &u[..active_width];
    future
        .map(|slice| {
            let safe = slice
                .iter()
                .position(|q| active_u.contains(q))
                .unwrap_or(slice.len());
            &slice[..safe]
        })
        .filter(|slice| !slice.is_empty())
}

pub(crate) fn dialog_gcd_composite_scratch_enabled() -> bool {
    std::env::var("DIALOG_GCD_COMPOSITE_SCRATCH")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_borrow_current_block_enabled() -> bool {
    // The GCD-walk peak (compress_block / shift / reverse_add, all at the same
    // height) is pinned by the composite body-scratch DEFICIT: at the widest
    // (early) steps the materialized sub/add wants ~2*active_width-1 clean lanes
    // for gated+carries, but the only |0> borrow there is the unwritten
    // future-log (block k+1..), leaving a fresh-allocated deficit on top of the
    // resident tx+ty+u+log.
    //
    // Novel observation: the CURRENT block's own compressed cells are also |0>
    // for the entire duration of that block's steps -- forward they are written
    // only by compress_block AFTER every step, reverse they are decompressed
    // into raw_block BEFORE every step -- yet the future-carry slice deliberately
    // starts at block k+1 and never offers them. Folding block k's own cells into
    // the body-scratch borrow shrinks the deficit (a pure qubit relabel, 0 added
    // Toffoli) and is value-exact: the body's measured uncompute restores them to
    // |0> before compress_block/decompress consumes them.
    std::env::var("DIALOG_GCD_BORROW_CURRENT_BLOCK")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_borrow_current_s2_enabled() -> bool {
    // Successor lever to BORROW_CURRENT_BLOCK for the K2 path. The current step's
    // own shift2 (`s2`) cell is provably |0> across its sub/add body window
    // (forward: written only by the later shift phase; reverse: already
    // uncomputed by reverse_unshift) and is restored to |0> by the body's
    // measured uncompute before the shift/unshift consumer. Folding it into the
    // composite-scratch borrow removes one fresh-allocated deficit lane at the
    // width-clamped GCD-walk binder steps (where active_width is pinned at N and
    // the future-log borrow has already shrunk a block), dropping the three
    // compressed-block tobitvector near-binders one qubit. Pure relabel, 0 added
    // Toffoli, value-exact on the reachable GCD support. Default off keeps the
    // accepted op stream byte-identical.
    std::env::var("DIALOG_GCD_BORROW_CURRENT_S2")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_skip_zero_edge_cshift_enabled() -> bool {
    std::env::var("DIALOG_GCD_SKIP_ZERO_EDGE_CSHIFT")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_skip_zero_edge_tobit_cshift_enabled() -> bool {
    dialog_gcd_skip_zero_edge_cshift_enabled()
        || std::env::var("DIALOG_GCD_SKIP_ZERO_EDGE_TOBIT_CSHIFT")
            .ok()
            .as_deref()
            == Some("1")
}

pub(crate) fn dialog_gcd_skip_zero_edge_tobit_fwd_cshift_enabled() -> bool {
    dialog_gcd_skip_zero_edge_tobit_cshift_enabled()
        || std::env::var("DIALOG_GCD_SKIP_ZERO_EDGE_TOBIT_FWD_CSHIFT")
            .ok()
            .as_deref()
            == Some("1")
}

pub(crate) fn dialog_gcd_skip_zero_edge_tobit_rev_cshift_enabled() -> bool {
    dialog_gcd_skip_zero_edge_tobit_cshift_enabled()
        || std::env::var("DIALOG_GCD_SKIP_ZERO_EDGE_TOBIT_REV_CSHIFT")
            .ok()
            .as_deref()
            == Some("1")
}

pub(crate) fn dialog_gcd_skip_zero_edge_apply_cshift_enabled() -> bool {
    dialog_gcd_skip_zero_edge_cshift_enabled()
        || std::env::var("DIALOG_GCD_SKIP_ZERO_EDGE_APPLY_CSHIFT")
            .ok()
            .as_deref()
            == Some("1")
}

pub(crate) fn dialog_gcd_skip_zero_edge_apply_double_cshift_enabled() -> bool {
    dialog_gcd_skip_zero_edge_apply_cshift_enabled()
        || std::env::var("DIALOG_GCD_SKIP_ZERO_EDGE_APPLY_DOUBLE_CSHIFT")
            .ok()
            .as_deref()
            == Some("1")
}

pub(crate) fn dialog_gcd_skip_zero_edge_apply_halve_cshift_enabled() -> bool {
    dialog_gcd_skip_zero_edge_apply_cshift_enabled()
        || std::env::var("DIALOG_GCD_SKIP_ZERO_EDGE_APPLY_HALVE_CSHIFT")
            .ok()
            .as_deref()
            == Some("1")
}

pub(crate) fn dialog_gcd_borrow_zero_raw_future_enabled() -> bool {
    // During a block-lifecycle tobitvector body, not every raw transcript cell is
    // live yet. Forward pass: slots greater than the current slot are still |0>
    // until their later branch/shift phases. Reverse pass: those greater slots
    // have already been uncomputed back to |0> before this slot's reverse_add.
    // Borrowing those cells as composite scratch is a pure retiming of clean
    // storage: the measured add/sub body restores them before any future use.
    std::env::var("DIALOG_GCD_BORROW_ZERO_RAW_FUTURE")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) struct DialogGcdCompositeScratch {
    lanes: Vec<QubitId>,
    owned: Vec<QubitId>,
}

pub(crate) fn dialog_gcd_build_composite_scratch(
    b: &mut B,
    future: Option<&[QubitId]>,
    u: &[QubitId],
    v: &[QubitId],
    compressed_log: &[QubitId],
    raw_block: &[QubitId],
    active_width: usize,
    step: usize,
) -> DialogGcdCompositeScratch {
    // The selected add/sub body is the dominant consumer of this composite
    // scratch (gated host + borrowed carries). Under the no-physical-c_in body
    // it needs only 2*body_len-1 == 2*body_w-3 lanes (vs 2*active_width-1), and
    // for the untrimmed fastpath body_w == active_width, so the demand drops by
    // exactly 2 lanes — the -1 peak qubit after the gap lane is also reclaimed.
    let body_start = if dialog_gcd_odd_u_lowbit_fastpath_enabled() {
        1
    } else {
        0
    };
    let body_w = dialog_gcd_body_carry_trunc_width(active_width, step);
    let body_len = body_w.saturating_sub(body_start);
    let nocin = dialog_gcd_selected_body_nocin_enabled()
        && !dialog_gcd_selected_body_nocin_keep_pool()
        && body_start >= 1
        && body_len >= 1;
    let stream_suffix = dialog_gcd_selected_body_stream_suffix_bits(step, body_len);
    let want = if !dialog_gcd_raw_tobitvector_materialized_sub_enabled() {
        // Low-scratch CONTROLLED body (cucc_sub/add_ctrl_lowq): it allocates its
        // own c_in+scratch internally and IGNORES borrowed_carries entirely. The
        // only remaining consumer of this composite scratch is the branch-bits
        // comparator host (dialog_gcd_ccx_cmp_gt_truncated_into_width_hosted),
        // whose transient is c_in (1) + carries (compare_bits) = compare_bits+1
        // clean lanes. Sizing the scratch to that comparator need only (instead
        // of the materialized body's 2*active_width-1) collapses the `owned`
        // deficit that pins the GCD-walk peak. Never exceed the legacy ask, and
        // keep >= 1 so an empty borrow set still yields a valid (clean) slice.
        //
        // When the Gidney-vented controlled body is active, it ALSO consumes this
        // composite scratch: it vents its forward carry chain onto active_width-1
        // BORROWED |0> lanes (restored by the measured uncompute). So bump `want`
        // to cover both consumers — still <= the materialized 2*active_width-1, so
        // the peak stays at the baseline. This guarantees the vented body finds
        // enough borrow that it does NOT fresh-alloc (which would spike the peak).
        let compare_bits = dialog_gcd_compare_bits_for_step(step, active_width);
        let comparator_need = compare_bits + 1;
        let body_need = if dialog_gcd_ctrl_body_vented_enabled() {
            active_width.saturating_sub(1)
        } else {
            0
        };
        comparator_need.max(body_need).min(2 * active_width - 1).max(1)
    } else if nocin && stream_suffix >= 2 {
        2 * (body_len - stream_suffix) + 1
    } else if nocin && dialog_gcd_selected_body_stream_top_enabled(step, body_len) && body_len >= 2
    {
        2 * (body_len - 1)
    } else if nocin {
        // Match the body's exact host demand; never exceed the legacy ask.
        (2 * body_len - 1).min(2 * active_width - 1)
    } else {
        2 * active_width - 1
    };
    let mut lanes = Vec::with_capacity(want);
    let mut push = |q: QubitId| {
        if lanes.len() < want
            && !lanes.contains(&q)
            && !raw_block.contains(&q)
            && !u[..active_width].contains(&q)
            && !v[..active_width].contains(&q)
        {
            lanes.push(q);
        }
    };
    if let Some(future) = dialog_gcd_runway_safe_future_prefix(future, u, active_width) {
        for &q in future {
            push(q);
        }
    }
    if dialog_gcd_borrow_current_block_enabled() {
        // Current block's own compressed cells: |0> across this block's steps
        // (forward written only at compress_block, reverse decompressed before
        // steps). They sit just BELOW the future-carry slice's start (k+1) and
        // are otherwise idle scratch. Restored to |0> by the body's measured
        // uncompute. Skip any that the runway parked onto active u (excluded by
        // push's active-u guard anyway, but kept explicit for clarity).
        let block_cells = dialog_gcd_compressed_sidecar_block(compressed_log, step);
        for &q in block_cells {
            push(q);
        }
    }
    for &q in &v[active_width..] {
        push(q);
    }
    for &q in &u[active_width..] {
        if !compressed_log.contains(&q) {
            push(q);
        }
    }
    if dialog_gcd_borrow_current_s2_enabled() && !raw_block.is_empty() {
        // The CURRENT step's own K2 shift2 (`s2`) cell is |0> across this step's
        // body window: forward it is written only by the later SHIFT phase
        // (after the sub body), reverse it has just been uncomputed by
        // reverse_unshift (before the add body). It is restored to |0> by the
        // body's measured uncompute before either consumer runs. Folding it into
        // the body-scratch borrow shrinks the fresh deficit by one lane at the
        // width-clamped binder steps (the same retiming trick as the current-block
        // compressed cells; pure relabel, 0 added Toffoli). The `push` closure
        // excludes all raw_block cells, so add it explicitly with the same
        // operand/duplicate guards. Disjoint from b0/b0_and_b1 (different slot
        // offset) and from u/v (raw_block is its own register).
        let group_size = dialog_gcd_sidecar_group_size();
        let slot = step % group_size;
        let s2 = raw_block[2 * group_size + slot];
        if lanes.len() < want
            && !lanes.contains(&s2)
            && !u[..active_width].contains(&s2)
            && !v[..active_width].contains(&s2)
        {
            lanes.push(s2);
        }
        if dialog_gcd_trio_width_notch_enabled() && slot == 0 && group_size >= 2 {
            let sibling_s2 = raw_block[2 * group_size + 1];
            if lanes.len() < want
                && !lanes.contains(&sibling_s2)
                && !u[..active_width].contains(&sibling_s2)
                && !v[..active_width].contains(&sibling_s2)
            {
                lanes.push(sibling_s2);
            }
        }
    }
    if dialog_gcd_borrow_zero_raw_future_enabled() && !raw_block.is_empty() {
        let group_size = dialog_gcd_sidecar_group_size();
        let slot = step % group_size;
        let mut push_raw_zero = |q: QubitId| {
            if lanes.len() < want
                && !lanes.contains(&q)
                && !u[..active_width].contains(&q)
                && !v[..active_width].contains(&q)
            {
                lanes.push(q);
            }
        };
        for future_slot in (slot + 1)..group_size {
            push_raw_zero(raw_block[2 * future_slot]);
            push_raw_zero(raw_block[2 * future_slot + 1]);
            if dialog_gcd_k2_enabled() {
                push_raw_zero(raw_block[2 * group_size + future_slot]);
            }
        }
    }
    let owned = b.alloc_qubits(want - lanes.len());
    if std::env::var("PROBE_SCRATCH").is_ok() && active_width >= 254 {
        eprintln!(
            "SCRATCH step={} aw={} body_w={} body_len={} want={} borrowed={} owned={}",
            step,
            active_width,
            body_w,
            body_len,
            want,
            lanes.len(),
            owned.len()
        );
    }
    lanes.extend_from_slice(&owned);
    DialogGcdCompositeScratch { lanes, owned }
}

pub(crate) fn dialog_gcd_pick_runway_safe_borrow_slice<'a>(
    future: Option<&'a [QubitId]>,
    u: &'a [QubitId],
    compressed_log: &[QubitId],
    active_width: usize,
) -> Option<&'a [QubitId]> {
    if !dialog_gcd_compressed_log_u_high_runway_enabled() {
        return dialog_gcd_pick_borrow_slice(future, u, active_width);
    }

    let safe_future = dialog_gcd_runway_safe_future_prefix(future, u, active_width);
    if dialog_gcd_late_borrow_uv_high_enabled() && active_width >= 1 {
        let want = 2 * active_width - 1;
        let short = safe_future.map_or(true, |slice| slice.len() < want);
        if short && u.len() >= active_width + want {
            let candidate = &u[active_width..active_width + want];
            // Parked cells can still carry unread transcript data.  Be
            // conservative: only use an in-place high-u fallback when it is
            // disjoint from every logical transcript cell, including clean
            // parked cells already consumed by reverse replay.
            if !dialog_gcd_slice_intersects(candidate, compressed_log) {
                return Some(candidate);
            }
        }
    }
    safe_future
}

pub(crate) fn dialog_gcd_host_reverse_raw_block_enabled() -> bool {
    // K=2 originally disabled this because the non-pair raw block widened to 9
    // lanes while the host search assumed 6. The host search below is now
    // raw_block_len-aware, but keep K2 hosting behind a separate experiment knob.
    if dialog_gcd_k2_enabled()
        && std::env::var("DIALOG_GCD_K2_HOST_RAW_BLOCK")
            .ok()
            .as_deref()
            != Some("1")
    {
        return false;
    }
    std::env::var("DIALOG_GCD_HOST_REVERSE_RAW_BLOCK")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_k2_apply_inplace_raw_block_enabled() -> bool {
    dialog_gcd_k2_pair_compress_enabled()
        && std::env::var("DIALOG_GCD_K2_APPLY_INPLACE_RAW_BLOCK")
            .ok()
            .as_deref()
            == Some("1")
}

pub(crate) fn dialog_gcd_k5_free_clean_block_during_shift_enabled() -> bool {
    dialog_gcd_k5_clean_block_enabled()
        && std::env::var("DIALOG_GCD_K5_FREE_CLEAN_BLOCK_DURING_SHIFT")
            .ok()
            .as_deref()
            == Some("1")
}

pub(crate) fn dialog_gcd_reverse_raw_block_host<'a>(
    u: &'a [QubitId],
    compressed_log: &'a [QubitId],
    block: usize,
) -> Option<&'a [QubitId]> {
    if !dialog_gcd_host_reverse_raw_block_enabled() {
        return None;
    }
    let (start, _) = dialog_gcd_compressed_sidecar_block_step_range(block);
    let active_width = dialog_gcd_tobitvector_active_width(start);
    let want = 2 * active_width - 1;
    let raw_bits = dialog_gcd_raw_block_len();
    if u.len().saturating_sub(active_width) >= want + raw_bits {
        let candidate = &u[u.len() - raw_bits..];
        if !dialog_gcd_compressed_log_u_high_runway_enabled()
            || !dialog_gcd_slice_intersects(candidate, compressed_log)
        {
            return Some(candidate);
        }
    }
    let future_start = dialog_gcd_compressed_sidecar_block_offset(block + 1);
    let future = compressed_log.get(future_start..)?;
    if future.len() < want + raw_bits {
        return None;
    }
    if !dialog_gcd_compressed_log_u_high_runway_enabled() {
        return Some(&future[future.len() - raw_bits..]);
    }
    // Keep the raw host after the largest possible carry+gated prefix and away
    // from active u.  With remapped runway cells the old final-six shortcut can
    // alias the growing reverse u prefix.
    future[want..]
        .windows(raw_bits)
        .rev()
        .find(|candidate| !dialog_gcd_slice_intersects(candidate, &u[..active_width]))
}

pub(crate) fn dialog_gcd_forward_raw_block_host<'a>(
    u: &'a [QubitId],
    compressed_log: &'a [QubitId],
    block: usize,
) -> Option<&'a [QubitId]> {
    if !dialog_gcd_host_reverse_raw_block_enabled() {
        return None;
    }
    let (start, _) = dialog_gcd_compressed_sidecar_block_step_range(block);
    let active_width = dialog_gcd_tobitvector_active_width(start);
    let want = 2 * active_width - 1;
    let raw_bits = dialog_gcd_raw_block_len();
    let future_start = dialog_gcd_compressed_sidecar_block_offset(block + 1);
    if let Some(future) = compressed_log.get(future_start..) {
        if future.len() >= want + raw_bits {
            if !dialog_gcd_compressed_log_u_high_runway_enabled() {
                return Some(&future[future.len() - raw_bits..]);
            }
            if let Some(candidate) = future[want..]
                .windows(raw_bits)
                .rev()
                .find(|candidate| !dialog_gcd_slice_intersects(candidate, &u[..active_width]))
            {
                return Some(candidate);
            }
        }
    }
    if u.len().saturating_sub(active_width) >= want + raw_bits {
        let candidate = &u[u.len() - raw_bits..];
        if !dialog_gcd_compressed_log_u_high_runway_enabled()
            || !dialog_gcd_slice_intersects(candidate, compressed_log)
        {
            Some(candidate)
        } else {
            None
        }
    } else {
        None
    }
}

pub(crate) fn dialog_gcd_compressed_sidecar_future_carry_slice(
    compressed_log: &[QubitId],
    step: usize,
    active_width: usize,
) -> Option<&[QubitId]> {
    if !dialog_gcd_raw_tobitvector_borrow_future_log_carries_enabled() {
        return None;
    }
    let carry_need = active_width.saturating_sub(1);
    // When hosting the gated register too, request up to carry(n-1)+gated(n)=2n-1
    // clean slots; the consumer splits the returned slice. Graceful: never return
    // fewer than carry_need (so carry borrowing is preserved), never more than
    // what the future region holds.
    let want = if dialog_gcd_host_gated_enabled() {
        2 * active_width - 1
    } else {
        carry_need
    };
    let next_block = dialog_gcd_compressed_sidecar_block_index(step) + 1;
    let start = dialog_gcd_compressed_sidecar_block_offset(next_block);
    compressed_log
        .get(start..)
        .filter(|future| future.len() >= carry_need)
        .map(|future| &future[..future.len().min(want)])
}

pub(crate) fn dialog_gcd_compressed_sidecar_block_step_range(block: usize) -> (usize, usize) {
    if (dialog_gcd_k5_tail6_graph_enabled() || dialog_gcd_k5_tail6_graph9_enabled())
        && block + 1 == dialog_gcd_compressed_sidecar_blocks()
    {
        return (
            dialog_gcd_active_iterations() - 6,
            dialog_gcd_active_iterations(),
        );
    }
    if dialog_gcd_k5_tail7_enabled()
        && block + 1 == dialog_gcd_compressed_sidecar_blocks()
    {
        return (
            dialog_gcd_active_iterations() - 7,
            dialog_gcd_active_iterations(),
        );
    }
    let group_size = dialog_gcd_sidecar_group_size();
    let start = block * group_size;
    let end = (start + group_size).min(dialog_gcd_active_iterations());
    (start, end)
}

pub(crate) fn dialog_gcd_copy_compressed_block_to_raw(
    b: &mut B,
    compressed_block: &[QubitId],
    raw_block: &[QubitId],
    steps: usize,
) {
    if dialog_gcd_k5_head11_enabled()
        && steps == 5
        && compressed_block.len() == DIALOG_GCD_K5_HEAD11_DATA_WIRES.len()
    {
        dialog_gcd_k5_head11_decompress_block_to_raw(
            b,
            compressed_block,
            raw_block,
            dialog_gcd_apply_replay_swap_host_enabled(),
        );
        // At step 0, u=p and every nonzero field factor v satisfies u>v, so
        // b0_and_b1 == b0. Keep the duplicate lane zero during apply replay;
        // the caller aliases the control and can lend this cell as clean scratch.
        b.cx(raw_block[0], raw_block[1]);
        return;
    }
    if dialog_gcd_k5_tail6_graph9_enabled()
        && steps == 6
        && compressed_block.len() == DIALOG_GCD_K5_TAIL6_GRAPH9_CODE_BITS
    {
        dialog_gcd_k5_tail6_graph9_decompress_block_to_raw(b, compressed_block, raw_block);
        return;
    }
    if dialog_gcd_k5_tail6_graph_enabled()
        && steps == 6
        && compressed_block.len() == DIALOG_GCD_K5_TAIL6_GRAPH_CODE_BITS
    {
        dialog_gcd_k5_tail6_graph_decompress_block_to_raw(b, compressed_block, raw_block);
        return;
    }
    if dialog_gcd_k5_tail7_enabled()
        && steps == 7
        && compressed_block.len() == DIALOG_GCD_K5_TAIL7_CODE_BITS
    {
        dialog_gcd_k5_tail7_decompress_block_to_raw(b, compressed_block, raw_block);
        return;
    }
    if dialog_gcd_k5_tail_pair1_enabled() && steps == 2 && compressed_block.len() == 1 {
        dialog_gcd_k5_tail_pair1_decompress_block_to_raw(
            b,
            compressed_block,
            raw_block,
            dialog_gcd_apply_replay_swap_host_enabled(),
        );
        return;
    }
    if dialog_gcd_k5_clean_block_enabled() {
        if steps == 5 {
            dialog_gcd_k5_decompress_block_to_raw(
                b,
                compressed_block,
                raw_block,
                dialog_gcd_apply_replay_swap_host_enabled(),
            );
        } else {
            dialog_gcd_k5_decompress_partial_block_to_raw(
                b,
                compressed_block,
                raw_block,
                steps,
                dialog_gcd_apply_replay_swap_host_enabled(),
            );
        }
        return;
    }
    if dialog_gcd_k2_pair_compress_enabled() {
        dialog_gcd_k2_pair_copy_compressed_block_to_raw(b, compressed_block, raw_block, steps);
        return;
    }
    let base_bits = DIALOG_GCD_HIGH_TAIL_ALIAS_BLOCK_BITS; // 5
    let raw_base = 2 * dialog_gcd_sidecar_group_size(); // 6
    assert_eq!(compressed_block.len(), dialog_gcd_block_bits());
    assert_eq!(raw_block.len(), dialog_gcd_raw_block_len());
    let swap_host = dialog_gcd_apply_replay_swap_host_enabled();
    for i in 0..base_bits {
        if swap_host {
            b.swap(compressed_block[i], raw_block[i]);
        } else {
            b.cx(compressed_block[i], raw_block[i]);
        }
    }
    emit_dialog_gcd_round763_compressor_inverse(b, &raw_block[0..raw_base]);
    // K=2 shift2 tail: compressed[5..] -> raw[6..] (raw, no compression).
    for j in base_bits..dialog_gcd_block_bits() {
        let r = raw_base + (j - base_bits);
        if swap_host {
            b.swap(compressed_block[j], raw_block[r]);
        } else {
            b.cx(compressed_block[j], raw_block[r]);
        }
    }
}

pub(crate) fn dialog_gcd_clear_raw_block_copy(
    b: &mut B,
    compressed_block: &[QubitId],
    raw_block: &[QubitId],
    steps: usize,
) {
    if dialog_gcd_k5_head11_enabled()
        && steps == 5
        && compressed_block.len() == DIALOG_GCD_K5_HEAD11_DATA_WIRES.len()
    {
        // Reconstruct the duplicated step-0 branch bit before running the exact
        // inverse head codec.
        b.cx(raw_block[0], raw_block[1]);
        dialog_gcd_k5_head11_compress_raw_to_block(
            b,
            compressed_block,
            raw_block,
            dialog_gcd_apply_replay_swap_host_enabled(),
        );
        return;
    }
    if dialog_gcd_k5_tail6_graph9_enabled()
        && steps == 6
        && compressed_block.len() == DIALOG_GCD_K5_TAIL6_GRAPH9_CODE_BITS
    {
        dialog_gcd_k5_tail6_graph9_compress_raw_to_block(b, compressed_block, raw_block);
        return;
    }
    if dialog_gcd_k5_tail6_graph_enabled()
        && steps == 6
        && compressed_block.len() == DIALOG_GCD_K5_TAIL6_GRAPH_CODE_BITS
    {
        dialog_gcd_k5_tail6_graph_compress_raw_to_block(b, compressed_block, raw_block);
        return;
    }
    if dialog_gcd_k5_tail7_enabled()
        && steps == 7
        && compressed_block.len() == DIALOG_GCD_K5_TAIL7_CODE_BITS
    {
        dialog_gcd_k5_tail7_compress_raw_to_block(b, compressed_block, raw_block);
        return;
    }
    if dialog_gcd_k5_tail_pair1_enabled() && steps == 2 && compressed_block.len() == 1 {
        dialog_gcd_k5_tail_pair1_compress_raw_to_block(
            b,
            compressed_block,
            raw_block,
            dialog_gcd_apply_replay_swap_host_enabled(),
        );
        return;
    }
    if dialog_gcd_k5_clean_block_enabled() {
        if steps == 5 {
            dialog_gcd_k5_compress_raw_to_block(
                b,
                compressed_block,
                raw_block,
                dialog_gcd_apply_replay_swap_host_enabled(),
            );
        } else {
            dialog_gcd_k5_compress_partial_raw_to_block(
                b,
                compressed_block,
                raw_block,
                steps,
                dialog_gcd_apply_replay_swap_host_enabled(),
            );
        }
        return;
    }
    if dialog_gcd_k2_pair_compress_enabled() {
        dialog_gcd_k2_pair_clear_raw_block_copy(b, compressed_block, raw_block, steps);
        return;
    }
    let base_bits = DIALOG_GCD_HIGH_TAIL_ALIAS_BLOCK_BITS;
    let raw_base = 2 * dialog_gcd_sidecar_group_size();
    assert_eq!(compressed_block.len(), dialog_gcd_block_bits());
    assert_eq!(raw_block.len(), dialog_gcd_raw_block_len());
    let swap_host = dialog_gcd_apply_replay_swap_host_enabled();
    // Inverse of copy: clear the shift2 tail first, then recompress the base.
    for j in base_bits..dialog_gcd_block_bits() {
        let r = raw_base + (j - base_bits);
        if swap_host {
            b.swap(compressed_block[j], raw_block[r]);
        } else {
            b.cx(compressed_block[j], raw_block[r]);
        }
    }
    emit_dialog_gcd_round763_compressor(b, &raw_block[0..raw_base]);
    for i in 0..base_bits {
        if swap_host {
            b.swap(compressed_block[i], raw_block[i]);
        } else {
            b.cx(compressed_block[i], raw_block[i]);
        }
    }
}

pub(crate) fn dialog_gcd_k2_pair_inplace_raw_frame(
    compressed_block: &[QubitId],
    raw0: QubitId,
) -> [QubitId; 6] {
    assert_eq!(compressed_block.len(), DIALOG_GCD_HIGH_TAIL_ALIAS_BLOCK_BITS);
    [
        raw0,
        compressed_block[0],
        compressed_block[2],
        compressed_block[3],
        compressed_block[1],
        compressed_block[4],
    ]
}

pub(crate) fn dialog_gcd_k2_pair_inplace_decompress_block(
    b: &mut B,
    compressed_block: &[QubitId],
    raw0: QubitId,
    steps: usize,
) -> [QubitId; 6] {
    assert_eq!(steps, 2, "in-place K2 apply currently requires full pair blocks");
    let raw_frame = dialog_gcd_k2_pair_inplace_raw_frame(compressed_block, raw0);
    let core = dialog_gcd_k2_pair_core(&raw_frame);
    emit_dialog_gcd_k2_pair_core_encoder_inverse(b, &core);
    raw_frame
}

pub(crate) fn dialog_gcd_k2_pair_inplace_clear_block(
    b: &mut B,
    compressed_block: &[QubitId],
    raw0: QubitId,
    steps: usize,
) {
    assert_eq!(steps, 2, "in-place K2 apply currently requires full pair blocks");
    let raw_frame = dialog_gcd_k2_pair_inplace_raw_frame(compressed_block, raw0);
    let core = dialog_gcd_k2_pair_core(&raw_frame);
    emit_dialog_gcd_k2_pair_core_encoder(b, &core);
}

pub(crate) fn emit_dialog_gcd_compressed_sidecar_tobitvector_steps_block_lifecycle(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    compressed_log: &[QubitId],
    raw_block: &[QubitId],
) {
    assert_eq!(u.len(), N);
    assert_eq!(v.len(), N);
    assert!(raw_block.is_empty() || raw_block.len() == dialog_gcd_raw_block_len());
    assert!(compressed_log.len() >= dialog_gcd_compressed_sidecar_bits());

    for block in 0..dialog_gcd_compressed_sidecar_blocks() {
        let (start, end) = dialog_gcd_compressed_sidecar_block_step_range(block);
        let block_steps = end - start;
        let hosted_raw_block = dialog_gcd_forward_raw_block_host(u, compressed_log, block);
        let owned_raw_block =
            if dialog_gcd_host_reverse_raw_block_enabled() && hosted_raw_block.is_none() {
                b.alloc_qubits(dialog_gcd_raw_block_len())
            } else {
                Vec::new()
            };
        let raw_block = hosted_raw_block.unwrap_or_else(|| {
            if owned_raw_block.is_empty() {
                raw_block
            } else {
                &owned_raw_block
            }
        });
        for step in start..end {
            let slot = step - start;
            if dialog_gcd_k5_constant_tail_stored_steps(block_steps)
                .is_some_and(|stored_steps| slot >= stored_steps)
            {
                let active_width = dialog_gcd_tobitvector_active_width(step);
                let shift_width = dialog_gcd_tobitvector_shift_width(active_width, step);
                let v_shift = &v[..shift_width];
                b.set_phase("dialog_gcd_compressed_block_tobitvector_tail7_constant_shift");
                dialog_gcd_shift_right_assuming_even(b, v_shift);
                dialog_gcd_shift_right_assuming_even(b, v_shift);
                continue;
            }
            let b0 = raw_block[2 * slot];
            let b0_and_b1 = raw_block[2 * slot + 1];
            let active_width = dialog_gcd_tobitvector_active_width(step);
            let u_active = &u[..active_width];
            let v_active = &v[..active_width];
            let compare_bits = dialog_gcd_compare_bits_for_step(step, active_width);

            let future = dialog_gcd_compressed_sidecar_future_carry_slice(
                compressed_log,
                step,
                active_width,
            );
            let composite_scratch = dialog_gcd_composite_scratch_enabled().then(|| {
                dialog_gcd_build_composite_scratch(
                    b,
                    future,
                    u,
                    v,
                    compressed_log,
                    raw_block,
                    active_width,
                    step,
                )
            });
            let borrowed_carries = composite_scratch.as_ref().map_or_else(
                || {
                    dialog_gcd_pick_runway_safe_borrow_slice(
                        future,
                        u,
                        compressed_log,
                        active_width,
                    )
                },
                |scratch| Some(scratch.lanes.as_slice()),
            );

            b.set_phase("dialog_gcd_compressed_block_tobitvector_branch_bits");
            b.cx(v[0], b0);
            if dialog_gcd_fused_branch_bits_enabled() {
                // Fused path derives b0_and_b1 from the in-flight comparator carry
                // and never materializes a separate `cmp` ancilla. Allocating it
                // here would add a dead live-qubit at the branch_bits peak instant
                // (peak is measured by simultaneously-live count, not qubit-id reuse),
                // so it is allocated only on the non-fused branch below.
                if dialog_gcd_branch_bits_host_comparator_enabled() {
                    // Host the comparator's c_in+carries transient on the idle
                    // future-log slice (the same slice the subtract borrows below;
                    // it is unwritten at the comparator instant) so branch_bits no
                    // longer allocates its own peak qubit. Value-exact; the slice is
                    // returned clean by the measured uncompute sweep.
                    dialog_gcd_ccx_cmp_gt_truncated_into_width_hosted(
                        b,
                        u_active,
                        v_active,
                        b0,
                        b0_and_b1,
                        compare_bits,
                        borrowed_carries,
                    );
                } else {
                    dialog_gcd_ccx_cmp_gt_truncated_into_width(
                        b,
                        u_active,
                        v_active,
                        b0,
                        b0_and_b1,
                        compare_bits,
                    );
                }
            } else {
                let cmp = b.alloc_qubit();
                dialog_gcd_cmp_gt_truncated_into_width(b, u_active, v_active, cmp, compare_bits);
                b.ccx(b0, cmp, b0_and_b1);
                dialog_gcd_cmp_gt_truncated_into_width(b, u_active, v_active, cmp, compare_bits);
                b.free(cmp);
            }

            b.set_phase("dialog_gcd_compressed_block_tobitvector_cswap");
            let cswap_width = dialog_gcd_tobitvector_cswap_width(active_width, step);
            for (i, (&ui, &vi)) in u[..cswap_width]
                .iter()
                .zip(v[..cswap_width].iter())
                .enumerate()
            {
                if i == 0 && dialog_gcd_odd_u_lowbit_fastpath_enabled() {
                    continue;
                }
                cswap(b, b0_and_b1, ui, vi);
            }

            b.set_phase("dialog_gcd_compressed_block_tobitvector_subtract");
            dialog_gcd_controlled_sub_selected(b, u_active, v_active, b0, borrowed_carries, step);
            if std::env::var("DIALOG_GCD_FREE_SCRATCH_BEFORE_SHIFT")
                .ok()
                .as_deref()
                == Some("1")
            {
                if let Some(scratch) = composite_scratch.as_ref() {
                    b.free_vec(&scratch.owned);
                }
            }

            b.set_phase("dialog_gcd_compressed_block_tobitvector_shift");
            let shift_width = dialog_gcd_tobitvector_shift_width(active_width, step);
            let v_shift = &v[..shift_width];
            dialog_gcd_shift_right_assuming_even(b, v_shift);
            if dialog_gcd_k2_enabled() {
                // K=2: record shift2 = NOT v_active[0] (v still even after the
                // first shift) into the sidecar, then conditionally shift v_active
                // right once more. Free 1-bit shift is a relabel; this 2nd shift is
                // data-dependent (cswap cascade), ~aw CCX.
                let s2 = dialog_gcd_block_raw_s2(raw_block, block_steps, slot);
                let v0 = v_active[0];
                if std::env::var("DIALOG_GCD_K2_FORCE0").ok().as_deref() != Some("1") {
                    b.cx(v0, s2);
                    b.x(s2);
                }
                let pairs = v_shift.len().saturating_sub(1);
                for i in 0..pairs {
                    if dialog_gcd_skip_zero_edge_tobit_fwd_cshift_enabled() && i + 1 == pairs {
                        continue;
                    }
                    let (lo, hi) = (v_shift[i], v_shift[i + 1]);
                    cswap(b, s2, lo, hi);
                }
            }
            if std::env::var("DIALOG_GCD_FREE_SCRATCH_BEFORE_SHIFT")
                .ok()
                .as_deref()
                != Some("1")
            {
                if let Some(scratch) = composite_scratch.as_ref() {
                    b.free_vec(&scratch.owned);
                }
            }
        }

        b.set_phase("dialog_gcd_compressed_block_tobitvector_compress_block");
        let base_bits = DIALOG_GCD_HIGH_TAIL_ALIAS_BLOCK_BITS; // 5
        let compressed_block = dialog_gcd_compressed_sidecar_block(compressed_log, start);
        if dialog_gcd_compressed_log_u_high_runway_enabled() {
            // A parked forward block is first written only after its high-u
            // hosts have left the active prefix.
            assert!(
                !dialog_gcd_slice_intersects(
                    compressed_block,
                    &u[..dialog_gcd_tobitvector_active_width(start)]
                ),
                "compressed-log runway overlaps active forward u prefix at block {block}"
            );
        }
        if dialog_gcd_k5_head11_enabled()
            && start == 0
            && block_steps == 5
            && compressed_block.len() == DIALOG_GCD_K5_HEAD11_DATA_WIRES.len()
        {
            dialog_gcd_k5_head11_compress_raw_to_block(b, compressed_block, raw_block, true);
        } else if dialog_gcd_k5_tail6_graph9_enabled()
            && block_steps == 6
            && compressed_block.len() == DIALOG_GCD_K5_TAIL6_GRAPH9_CODE_BITS
        {
            dialog_gcd_k5_tail6_graph9_compress_raw_to_block(b, compressed_block, raw_block);
        } else if dialog_gcd_k5_tail6_graph_enabled()
            && block_steps == 6
            && compressed_block.len() == DIALOG_GCD_K5_TAIL6_GRAPH_CODE_BITS
        {
            dialog_gcd_k5_tail6_graph_compress_raw_to_block(b, compressed_block, raw_block);
        } else if dialog_gcd_k5_tail7_enabled()
            && block_steps == 7
            && compressed_block.len() == DIALOG_GCD_K5_TAIL7_CODE_BITS
        {
            dialog_gcd_k5_tail7_compress_raw_to_block(b, compressed_block, raw_block);
        } else if dialog_gcd_k5_tail_pair1_enabled()
            && end - start == 2
            && compressed_block.len() == 1
        {
            dialog_gcd_k5_tail_pair1_compress_raw_to_block(
                b,
                compressed_block,
                raw_block,
                true,
            );
        } else if dialog_gcd_k5_clean_block_enabled() {
            if end - start == 5 {
                dialog_gcd_k5_compress_raw_to_block(b, compressed_block, raw_block, true);
            } else {
                dialog_gcd_k5_compress_partial_raw_to_block(
                    b,
                    compressed_block,
                    raw_block,
                    end - start,
                    true,
                );
            }
        } else if dialog_gcd_k2_pair_compress_enabled() {
            dialog_gcd_k2_pair_clear_raw_block_copy(b, compressed_block, raw_block, end - start);
        } else {
            let raw_base = 2 * dialog_gcd_sidecar_group_size(); // 6
            emit_dialog_gcd_round763_compressor(b, &raw_block[0..raw_base]);
            for i in 0..base_bits {
                b.swap(raw_block[i], compressed_block[i]);
            }
            // K=2: stash the shift2 bits raw[raw_base..] into compressed_block[5..].
            for j in base_bits..dialog_gcd_block_bits() {
                b.swap(raw_block[raw_base + (j - base_bits)], compressed_block[j]);
            }
        }
        if !owned_raw_block.is_empty() {
            b.free_vec(&owned_raw_block);
        }
    }
}

pub(crate) fn emit_dialog_gcd_compressed_sidecar_tobitvector_steps_reverse_block_lifecycle(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    compressed_log: &[QubitId],
    raw_block: &[QubitId],
) {
    assert_eq!(u.len(), N);
    assert_eq!(v.len(), N);
    assert!(raw_block.is_empty() || raw_block.len() == dialog_gcd_raw_block_len());
    assert!(compressed_log.len() >= dialog_gcd_compressed_sidecar_bits());

    for block in (0..dialog_gcd_compressed_sidecar_blocks()).rev() {
        let (start, end) = dialog_gcd_compressed_sidecar_block_step_range(block);
        let block_steps = end - start;
        let compressed_block = dialog_gcd_compressed_sidecar_block(compressed_log, start);
        let hosted_raw_block = dialog_gcd_reverse_raw_block_host(u, compressed_log, block);
        let owned_raw_block =
            if dialog_gcd_host_reverse_raw_block_enabled() && hosted_raw_block.is_none() {
                b.alloc_qubits(dialog_gcd_raw_block_len())
            } else {
                Vec::new()
            };
        let raw_block = hosted_raw_block.unwrap_or_else(|| {
            if owned_raw_block.is_empty() {
                raw_block
            } else {
                &owned_raw_block
            }
        });

        b.set_phase("dialog_gcd_compressed_block_tobitvector_reverse_decompress_block");
        if dialog_gcd_compressed_log_u_high_runway_enabled() {
            // A parked block must be consumed while all of its high-u hosts are
            // outside this block's active prefix.
            assert!(
                !dialog_gcd_slice_intersects(
                    compressed_block,
                    &u[..dialog_gcd_tobitvector_active_width(start)]
                ),
                "compressed-log runway overlaps active reverse u prefix at block {block}"
            );
        }
        {
            let base_bits = DIALOG_GCD_HIGH_TAIL_ALIAS_BLOCK_BITS; // 5
            if dialog_gcd_k5_head11_enabled()
                && start == 0
                && block_steps == 5
                && compressed_block.len() == DIALOG_GCD_K5_HEAD11_DATA_WIRES.len()
            {
                dialog_gcd_k5_head11_decompress_block_to_raw(
                    b,
                    compressed_block,
                    raw_block,
                    true,
                );
            } else if dialog_gcd_k5_tail6_graph9_enabled()
                && block_steps == 6
                && compressed_block.len() == DIALOG_GCD_K5_TAIL6_GRAPH9_CODE_BITS
            {
                dialog_gcd_k5_tail6_graph9_decompress_block_to_raw(
                    b,
                    compressed_block,
                    raw_block,
                );
            } else if dialog_gcd_k5_tail6_graph_enabled()
                && block_steps == 6
                && compressed_block.len() == DIALOG_GCD_K5_TAIL6_GRAPH_CODE_BITS
            {
                dialog_gcd_k5_tail6_graph_decompress_block_to_raw(
                    b,
                    compressed_block,
                    raw_block,
                );
            } else if dialog_gcd_k5_tail7_enabled()
                && block_steps == 7
                && compressed_block.len() == DIALOG_GCD_K5_TAIL7_CODE_BITS
            {
                dialog_gcd_k5_tail7_decompress_block_to_raw(
                    b,
                    compressed_block,
                    raw_block,
                );
            } else if dialog_gcd_k5_tail_pair1_enabled()
                && end - start == 2
                && compressed_block.len() == 1
            {
                dialog_gcd_k5_tail_pair1_decompress_block_to_raw(
                    b,
                    compressed_block,
                    raw_block,
                    true,
                );
            } else if dialog_gcd_k5_clean_block_enabled() {
                if end - start == 5 {
                    dialog_gcd_k5_decompress_block_to_raw(b, compressed_block, raw_block, true);
                } else {
                    dialog_gcd_k5_decompress_partial_block_to_raw(
                        b,
                        compressed_block,
                        raw_block,
                        end - start,
                        true,
                    );
                }
            } else if dialog_gcd_k2_pair_compress_enabled() {
                dialog_gcd_k2_pair_copy_compressed_block_to_raw(
                    b,
                    compressed_block,
                    raw_block,
                    end - start,
                );
            } else {
                let raw_base = 2 * dialog_gcd_sidecar_group_size(); // 6
                for i in 0..base_bits {
                    b.swap(compressed_block[i], raw_block[i]);
                }
                emit_dialog_gcd_round763_compressor_inverse(b, &raw_block[0..raw_base]);
                // K=2: bring the shift2 bits compressed[5..] -> raw[raw_base..].
                for j in base_bits..dialog_gcd_block_bits() {
                    b.swap(compressed_block[j], raw_block[raw_base + (j - base_bits)]);
                }
            }
        }

        for step in (start..end).rev() {
            let slot = step - start;
            if dialog_gcd_k5_constant_tail_stored_steps(block_steps)
                .is_some_and(|stored_steps| slot >= stored_steps)
            {
                let active_width = dialog_gcd_tobitvector_active_width(step);
                let shift_width = dialog_gcd_tobitvector_shift_width(active_width, step);
                let v_shift = &v[..shift_width];
                b.set_phase(
                    "dialog_gcd_compressed_block_tobitvector_reverse_tail7_constant_unshift",
                );
                dialog_gcd_unshift_right_assuming_even(b, v_shift);
                dialog_gcd_unshift_right_assuming_even(b, v_shift);
                continue;
            }
            let b0 = raw_block[2 * slot];
            let b0_and_b1 = raw_block[2 * slot + 1];
            let active_width = dialog_gcd_tobitvector_active_width(step);
            let u_active = &u[..active_width];
            let v_active = &v[..active_width];
            let compare_bits = dialog_gcd_compare_bits_for_step(step, active_width);

            b.set_phase("dialog_gcd_compressed_block_tobitvector_reverse_unshift");
            let shift_width = dialog_gcd_tobitvector_shift_width(active_width, step);
            let v_shift = &v[..shift_width];
            if dialog_gcd_k2_enabled() {
                // mirror of forward K=2: conditional un-shift (reverse cswap order),
                // then uncompute s2 back to |0> (v_active[0] is restored after the
                // un-shift to the value s2 was derived from).
                let s2 = dialog_gcd_block_raw_s2(raw_block, block_steps, slot);
                let pairs = v_shift.len().saturating_sub(1);
                for i in (0..pairs).rev() {
                    if dialog_gcd_skip_zero_edge_tobit_rev_cshift_enabled() && i + 1 == pairs {
                        continue;
                    }
                    let (lo, hi) = (v_shift[i], v_shift[i + 1]);
                    cswap(b, s2, lo, hi);
                }
                let v0 = v_active[0];
                if std::env::var("DIALOG_GCD_K2_FORCE0").ok().as_deref() != Some("1") {
                    b.x(s2);
                    b.cx(v0, s2);
                }
            }
            dialog_gcd_unshift_right_assuming_even(b, v_shift);

            b.set_phase("dialog_gcd_compressed_block_tobitvector_reverse_add");
            let future = dialog_gcd_compressed_sidecar_future_carry_slice(
                compressed_log,
                step,
                active_width,
            );
            let composite_scratch = dialog_gcd_composite_scratch_enabled().then(|| {
                dialog_gcd_build_composite_scratch(
                    b,
                    future,
                    u,
                    v,
                    compressed_log,
                    raw_block,
                    active_width,
                    step,
                )
            });
            let borrowed_carries = composite_scratch.as_ref().map_or_else(
                || {
                    dialog_gcd_pick_runway_safe_borrow_slice(
                        future,
                        u,
                        compressed_log,
                        active_width,
                    )
                },
                |scratch| Some(scratch.lanes.as_slice()),
            );
            dialog_gcd_controlled_add_selected(b, u_active, v_active, b0, borrowed_carries, step);

            b.set_phase("dialog_gcd_compressed_block_tobitvector_reverse_cswap");
            let cswap_width = dialog_gcd_tobitvector_cswap_width(active_width, step);
            for (i, (&ui, &vi)) in u[..cswap_width]
                .iter()
                .zip(v[..cswap_width].iter())
                .enumerate()
            {
                if i == 0 && dialog_gcd_odd_u_lowbit_fastpath_enabled() {
                    continue;
                }
                cswap(b, b0_and_b1, ui, vi);
            }

            b.set_phase("dialog_gcd_compressed_block_tobitvector_reverse_branch_bits");
            if dialog_gcd_reverse_branch_conditional_replay_enabled() {
                let phase = b.alloc_bit();
                b.hmr(b0_and_b1, phase);
                dialog_gcd_cmp_gt_truncated_phase_conditioned_hosted(
                    b,
                    u_active,
                    v_active,
                    b0,
                    phase,
                    compare_bits,
                    borrowed_carries,
                );
            } else if dialog_gcd_fused_branch_bits_enabled() {
                // Fused path: no separate `cmp` ancilla (derives b0_and_b1 from the
                // comparator carry). Allocating it would add a dead live-qubit at the
                // reverse_branch_bits peak instant, so allocate only on the non-fused
                // branch below. See forward lifecycle for the rationale.
                if dialog_gcd_branch_bits_host_comparator_enabled() {
                    // Mirror of the forward path: host the comparator transient on
                    // the idle future-log slice (same slice the add borrowed above).
                    dialog_gcd_ccx_cmp_gt_truncated_into_width_hosted(
                        b,
                        u_active,
                        v_active,
                        b0,
                        b0_and_b1,
                        compare_bits,
                        borrowed_carries,
                    );
                } else {
                    dialog_gcd_ccx_cmp_gt_truncated_into_width(
                        b,
                        u_active,
                        v_active,
                        b0,
                        b0_and_b1,
                        compare_bits,
                    );
                }
            } else {
                let cmp = b.alloc_qubit();
                dialog_gcd_cmp_gt_truncated_into_width(b, u_active, v_active, cmp, compare_bits);
                b.ccx(b0, cmp, b0_and_b1);
                dialog_gcd_cmp_gt_truncated_into_width(b, u_active, v_active, cmp, compare_bits);
                b.free(cmp);
            }
            b.cx(v[0], b0);
            if let Some(scratch) = composite_scratch {
                b.free_vec(&scratch.owned);
            }
        }
        if !owned_raw_block.is_empty() {
            b.free_vec(&owned_raw_block);
        }
    }
}

pub(crate) fn emit_dialog_gcd_compressed_sidecar_apply_bitvector_block_lifecycle(
    b: &mut B,
    compressed_log: &[QubitId],
    x: &[QubitId],
    y: &[QubitId],
    p: U256,
    raw_block: &[QubitId],
) {
    assert_eq!(x.len(), N);
    assert_eq!(y.len(), N);
    let inplace_raw = dialog_gcd_k2_apply_inplace_raw_block_enabled();
    if inplace_raw {
        assert!(raw_block.is_empty());
    } else {
        assert_eq!(raw_block.len(), dialog_gcd_raw_block_len());
    }
    let inplace_raw0 = if inplace_raw {
        Some(b.alloc_qubit())
    } else {
        None
    };

    for block in (0..dialog_gcd_compressed_sidecar_blocks()).rev() {
        let (start, end) = dialog_gcd_compressed_sidecar_block_step_range(block);
        let block_steps = end - start;
        let compressed_block = dialog_gcd_compressed_sidecar_block(compressed_log, start);
        let head11_block = dialog_gcd_k5_head11_enabled()
            && start == 0
            && block_steps == 5
            && compressed_block.len() == DIALOG_GCD_K5_HEAD11_DATA_WIRES.len();
        let head11_step_frame =
            head11_block && dialog_gcd_k5_head11_apply_step_frame_enabled();
        let step_frame = head11_step_frame.then(|| b.alloc_qubits(3));
        if head11_step_frame {
            b.set_phase("dialog_gcd_compressed_block_apply_release_raw_host");
            b.free_vec(raw_block);
        }

        b.set_phase("dialog_gcd_compressed_block_apply_decompress_block");
        let raw_frame = if head11_step_frame {
            None
        } else {
            inplace_raw0.map(|raw0| {
            dialog_gcd_k2_pair_inplace_decompress_block(b, compressed_block, raw0, end - start)
            })
        };
        if raw_frame.is_none() && !head11_step_frame {
            dialog_gcd_copy_compressed_block_to_raw(b, compressed_block, raw_block, end - start);
        }
        if head11_block && !head11_step_frame {
            b.free(raw_block[1]);
        }
        let released_code_bits = if raw_frame.is_none()
            && !head11_step_frame
            && dialog_gcd_apply_replay_swap_host_enabled()
        {
            let requested = if block_steps == 6
                && compressed_block.len() == DIALOG_GCD_K5_TAIL6_GRAPH9_CODE_BITS
            {
                dialog_gcd_k5_release_decoded_tail_bits()
            } else {
                dialog_gcd_k5_release_decoded_block_bits()
            };
            requested.min(compressed_block.len())
        } else {
            0
        };
        let retained_code_bits = compressed_block.len() - released_code_bits;
        let released_code = &compressed_block[retained_code_bits..];
        b.free_vec(released_code);
        let scale_release_bits =
            dialog_gcd_k5_release_scale_bits().min(retained_code_bits);
        let scale_released_code =
            &compressed_block[retained_code_bits - scale_release_bits..retained_code_bits];
        let shift_clean_code =
            &compressed_block[..retained_code_bits - scale_release_bits];
        let raw = raw_frame.as_ref().map_or(raw_block, |frame| &frame[..]);
        let tail_clean_scratch = if dialog_gcd_k5_tail_pair1_enabled()
            && end - start == 2
            && compressed_block.len() == 1
        {
            raw.iter()
                .enumerate()
                .filter_map(|(index, &q)| {
                    (!matches!(index, 0 | 2 | 10 | 11)).then_some(q)
                })
                .chain(compressed_block.iter().copied())
                .collect::<Vec<_>>()
        } else {
            Vec::new()
        };
        let block_clean_scratch = if head11_step_frame {
            &[]
        } else if !tail_clean_scratch.is_empty() {
            tail_clean_scratch.as_slice()
        } else if !inplace_raw && dialog_gcd_apply_replay_swap_host_enabled() {
            &compressed_block[..retained_code_bits]
        } else {
            &[]
        };

        for step in (start..end).rev() {
            let slot = step - start;
            let constant_tail_stored_steps =
                dialog_gcd_k5_constant_tail_stored_steps(block_steps);
            if constant_tail_stored_steps.is_some_and(|stored_steps| slot >= stored_steps) {
                let stored_steps = constant_tail_stored_steps.expect("checked above");
                if !scale_released_code.is_empty() {
                    b.set_phase("dialog_gcd_compressed_block_apply_scale_release");
                    b.free_vec(scale_released_code);
                }
                b.set_phase("dialog_gcd_compressed_block_apply_tail7_constant_double_y");
                if dialog_gcd_k5_fixed_tail_apply_enabled() {
                    dialog_gcd_fixed_double_twice_y(b, y, p);
                    if !scale_released_code.is_empty() {
                        b.set_phase("dialog_gcd_compressed_block_apply_scale_reacquire");
                        b.reacquire_vec(scale_released_code);
                    }
                    continue;
                }
                let one = raw[12];
                if slot + 1 == block_steps {
                    b.x(one);
                }
                if dialog_gcd_apply_fused_fold_enabled() {
                    dialog_gcd_fused_double_y_at_step(b, y, p, one, Some(step));
                } else {
                    mod_double_inplace_fast(b, y, p);
                    cmod_double_inplace_lazy(b, y, p, one);
                }
                if slot == stored_steps {
                    b.x(one);
                }
                if !scale_released_code.is_empty() {
                    b.set_phase("dialog_gcd_compressed_block_apply_scale_reacquire");
                    b.reacquire_vec(scale_released_code);
                }
                continue;
            }
            if let Some(frame) = step_frame.as_ref() {
                b.set_phase("dialog_gcd_compressed_block_apply_load_step_frame");
                dialog_gcd_k5_head11_toggle_step_frame(b, compressed_block, frame, slot);
            }
            let (b0, b0_and_b1, s2) = if let Some(frame) = step_frame.as_ref() {
                (frame[0], frame[1], frame[2])
            } else {
                let b0 = raw[2 * slot];
                let b0_and_b1 = if head11_block && slot == 0 {
                    raw[0]
                } else {
                    raw[2 * slot + 1]
                };
                let s2 = dialog_gcd_block_raw_s2(raw, block_steps, slot);
                (b0, b0_and_b1, s2)
            };

            if !scale_released_code.is_empty() {
                b.set_phase("dialog_gcd_compressed_block_apply_scale_release");
                b.free_vec(scale_released_code);
            }
            b.set_phase("dialog_gcd_compressed_block_apply_double_y");
            let apply_k2 = dialog_gcd_k2_enabled()
                && std::env::var("DIALOG_GCD_K2_NO_APPLY").ok().as_deref() != Some("1");
            let free_clean_code = !inplace_raw
                && dialog_gcd_apply_replay_swap_host_enabled()
                && dialog_gcd_k5_free_clean_block_during_shift_enabled();
            if free_clean_code {
                b.free_vec(shift_clean_code);
            }
            if apply_k2 && dialog_gcd_apply_fused_fold_enabled() {
                // Fuse mod_double_inplace_fast + cmod_double_inplace_lazy into a
                // single shared carry chain (value-identical; see fn doc).
                dialog_gcd_fused_double_y_at_step(b, y, p, s2, Some(step));
            } else {
                mod_double_inplace_fast(b, y, p);
                if apply_k2 {
                    // mirror the forward K=2 second shift: conditional 2nd double of y.
                    // MUST use the lazy (Solinas, truncated) controlled double so it
                    // composes with the uncontrolled mod_double_inplace_fast above.
                    cmod_double_inplace_lazy(b, y, p, s2);
                }
            }
            if free_clean_code {
                b.reacquire_vec(shift_clean_code);
            }
            if !scale_released_code.is_empty() {
                b.set_phase("dialog_gcd_compressed_block_apply_scale_reacquire");
                b.reacquire_vec(scale_released_code);
            }

            b.set_phase("dialog_gcd_compressed_block_apply_cadd");
            if dialog_gcd_raw_apply_materialized_special_add_enabled() {
                let owned_clean_scratch = if inplace_raw || head11_step_frame {
                    b.alloc_qubits(dialog_gcd_block_bits())
                } else {
                    Vec::new()
                };
                let clean_scratch = if inplace_raw || head11_step_frame {
                    owned_clean_scratch.as_slice()
                } else {
                    block_clean_scratch
                };
                dialog_gcd_cmod_add_materialized_pseudomersenne_with_clean_scratch_at_step(
                    b,
                    y,
                    x,
                    b0,
                    p,
                    clean_scratch,
                    Some(step),
                );
                if inplace_raw || head11_step_frame {
                    b.free_vec(&owned_clean_scratch);
                }
            } else if dialog_gcd_raw_apply_direct_special_add_enabled() {
                dialog_gcd_cmod_add_pseudomersenne_lowq(b, y, x, b0, p);
            } else {
                cmod_add_qq_lowq(b, y, x, b0, p);
            }

            b.set_phase("dialog_gcd_compressed_block_apply_cswap");
            for (&xi, &yi) in x.iter().zip(y.iter()) {
                cswap(b, b0_and_b1, xi, yi);
            }
            if let Some(frame) = step_frame.as_ref() {
                b.set_phase("dialog_gcd_compressed_block_apply_clear_step_frame");
                dialog_gcd_k5_head11_toggle_step_frame(b, compressed_block, frame, slot);
            }
        }

        if !released_code.is_empty() {
            b.set_phase("dialog_gcd_compressed_block_apply_reacquire_block");
            b.reacquire_vec(released_code);
        }
        if head11_block && !head11_step_frame {
            b.reacquire(raw_block[1]);
        }
        b.set_phase("dialog_gcd_compressed_block_apply_clear_block_copy");
        if let Some(raw0) = inplace_raw0 {
            dialog_gcd_k2_pair_inplace_clear_block(b, compressed_block, raw0, end - start);
        } else if !head11_step_frame {
            dialog_gcd_clear_raw_block_copy(b, compressed_block, raw_block, end - start);
        }
        if let Some(frame) = step_frame {
            b.free_vec(&frame);
        }
        if head11_step_frame {
            b.set_phase("dialog_gcd_compressed_block_apply_reacquire_raw_host");
            b.reacquire_vec(raw_block);
        }
    }

    if let Some(raw0) = inplace_raw0 {
        b.free(raw0);
    }
}

pub(crate) fn emit_dialog_gcd_compressed_sidecar_apply_bitvector_reverse_exact_block_lifecycle(
    b: &mut B,
    compressed_log: &[QubitId],
    x: &[QubitId],
    y: &[QubitId],
    p: U256,
    raw_block: &[QubitId],
) {
    assert_eq!(x.len(), N);
    assert_eq!(y.len(), N);
    let inplace_raw = dialog_gcd_k2_apply_inplace_raw_block_enabled();
    if inplace_raw {
        assert!(raw_block.is_empty());
    } else {
        assert_eq!(raw_block.len(), dialog_gcd_raw_block_len());
    }
    let inplace_raw0 = if inplace_raw {
        Some(b.alloc_qubit())
    } else {
        None
    };

    for block in 0..dialog_gcd_compressed_sidecar_blocks() {
        let (start, end) = dialog_gcd_compressed_sidecar_block_step_range(block);
        let block_steps = end - start;
        let compressed_block = dialog_gcd_compressed_sidecar_block(compressed_log, start);
        let head11_block = dialog_gcd_k5_head11_enabled()
            && start == 0
            && block_steps == 5
            && compressed_block.len() == DIALOG_GCD_K5_HEAD11_DATA_WIRES.len();
        let head11_step_frame =
            head11_block && dialog_gcd_k5_head11_apply_step_frame_enabled();
        let step_frame = head11_step_frame.then(|| b.alloc_qubits(3));
        if head11_step_frame {
            b.set_phase("dialog_gcd_compressed_block_apply_reverse_release_raw_host");
            b.free_vec(raw_block);
        }

        b.set_phase("dialog_gcd_compressed_block_apply_reverse_decompress_block");
        let raw_frame = if head11_step_frame {
            None
        } else {
            inplace_raw0.map(|raw0| {
            dialog_gcd_k2_pair_inplace_decompress_block(b, compressed_block, raw0, end - start)
            })
        };
        if raw_frame.is_none() && !head11_step_frame {
            dialog_gcd_copy_compressed_block_to_raw(b, compressed_block, raw_block, end - start);
        }
        if head11_block && !head11_step_frame {
            b.free(raw_block[1]);
        }
        let released_code_bits = if raw_frame.is_none()
            && !head11_step_frame
            && dialog_gcd_apply_replay_swap_host_enabled()
        {
            let requested = if block_steps == 6
                && compressed_block.len() == DIALOG_GCD_K5_TAIL6_GRAPH9_CODE_BITS
            {
                dialog_gcd_k5_release_decoded_tail_bits()
            } else {
                dialog_gcd_k5_release_decoded_block_bits()
            };
            requested.min(compressed_block.len())
        } else {
            0
        };
        let retained_code_bits = compressed_block.len() - released_code_bits;
        let released_code = &compressed_block[retained_code_bits..];
        b.free_vec(released_code);
        let scale_release_bits =
            dialog_gcd_k5_release_scale_bits().min(retained_code_bits);
        let scale_released_code =
            &compressed_block[retained_code_bits - scale_release_bits..retained_code_bits];
        let shift_clean_code =
            &compressed_block[..retained_code_bits - scale_release_bits];
        let raw = raw_frame.as_ref().map_or(raw_block, |frame| &frame[..]);
        let tail_clean_scratch = if dialog_gcd_k5_tail_pair1_enabled()
            && end - start == 2
            && compressed_block.len() == 1
        {
            raw.iter()
                .enumerate()
                .filter_map(|(index, &q)| {
                    (!matches!(index, 0 | 2 | 10 | 11)).then_some(q)
                })
                .chain(compressed_block.iter().copied())
                .collect::<Vec<_>>()
        } else {
            Vec::new()
        };
        let block_clean_scratch = if head11_step_frame {
            &[]
        } else if !tail_clean_scratch.is_empty() {
            tail_clean_scratch.as_slice()
        } else if !inplace_raw && dialog_gcd_apply_replay_swap_host_enabled() {
            &compressed_block[..retained_code_bits]
        } else {
            &[]
        };

        for step in start..end {
            let slot = step - start;
            let constant_tail_stored_steps =
                dialog_gcd_k5_constant_tail_stored_steps(block_steps);
            if constant_tail_stored_steps.is_some_and(|stored_steps| slot >= stored_steps) {
                let stored_steps = constant_tail_stored_steps.expect("checked above");
                if !scale_released_code.is_empty() {
                    b.set_phase("dialog_gcd_compressed_block_apply_reverse_scale_release");
                    b.free_vec(scale_released_code);
                }
                b.set_phase(
                    "dialog_gcd_compressed_block_apply_reverse_tail7_constant_halve_y",
                );
                if dialog_gcd_k5_fixed_tail_apply_enabled() {
                    dialog_gcd_fixed_halve_twice_y(b, y, p);
                    if !scale_released_code.is_empty() {
                        b.set_phase(
                            "dialog_gcd_compressed_block_apply_reverse_scale_reacquire",
                        );
                        b.reacquire_vec(scale_released_code);
                    }
                    continue;
                }
                let one = raw[12];
                if slot == stored_steps {
                    b.x(one);
                }
                if dialog_gcd_apply_fused_fold_enabled()
                    && std::env::var("DIALOG_GCD_FUSE_HALVE_OFF").ok().as_deref() != Some("1")
                {
                    dialog_gcd_fused_halve_y_at_step(b, y, p, one, Some(step));
                } else {
                    mod_halve_inplace_fast(b, y, p);
                    cmod_halve_inplace_lazy(b, y, p, one);
                }
                if slot + 1 == block_steps {
                    b.x(one);
                }
                if !scale_released_code.is_empty() {
                    b.set_phase("dialog_gcd_compressed_block_apply_reverse_scale_reacquire");
                    b.reacquire_vec(scale_released_code);
                }
                continue;
            }
            if let Some(frame) = step_frame.as_ref() {
                b.set_phase("dialog_gcd_compressed_block_apply_reverse_load_step_frame");
                dialog_gcd_k5_head11_toggle_step_frame(b, compressed_block, frame, slot);
            }
            let (b0, b0_and_b1, s2) = if let Some(frame) = step_frame.as_ref() {
                (frame[0], frame[1], frame[2])
            } else {
                let b0 = raw[2 * slot];
                let b0_and_b1 = if head11_block && slot == 0 {
                    raw[0]
                } else {
                    raw[2 * slot + 1]
                };
                let s2 = dialog_gcd_block_raw_s2(raw, block_steps, slot);
                (b0, b0_and_b1, s2)
            };

            b.set_phase("dialog_gcd_compressed_block_apply_reverse_cswap");
            for (&xi, &yi) in x.iter().zip(y.iter()) {
                cswap(b, b0_and_b1, xi, yi);
            }

            b.set_phase("dialog_gcd_compressed_block_apply_reverse_csub");
            if dialog_gcd_raw_apply_reverse_materialized_special_sub_enabled() {
                let owned_clean_scratch = if inplace_raw || head11_step_frame {
                    b.alloc_qubits(dialog_gcd_block_bits())
                } else {
                    Vec::new()
                };
                let clean_scratch = if inplace_raw || head11_step_frame {
                    owned_clean_scratch.as_slice()
                } else {
                    block_clean_scratch
                };
                dialog_gcd_cmod_sub_materialized_pseudomersenne_with_clean_scratch_at_step(
                    b,
                    y,
                    x,
                    b0,
                    p,
                    clean_scratch,
                    Some(step),
                );
                if inplace_raw || head11_step_frame {
                    b.free_vec(&owned_clean_scratch);
                }
            } else if dialog_gcd_raw_apply_reverse_fast_sub_enabled() {
                cmod_sub_qq(b, y, x, b0, p);
            } else {
                cmod_sub_qq_lowq(b, y, x, b0, p);
            }

            if !scale_released_code.is_empty() {
                b.set_phase("dialog_gcd_compressed_block_apply_reverse_scale_release");
                b.free_vec(scale_released_code);
            }
            b.set_phase("dialog_gcd_compressed_block_apply_reverse_halve_y");
            let apply_k2 = dialog_gcd_k2_enabled()
                && std::env::var("DIALOG_GCD_K2_NO_APPLY").ok().as_deref() != Some("1");
            let free_clean_code = !inplace_raw
                && dialog_gcd_apply_replay_swap_host_enabled()
                && dialog_gcd_k5_free_clean_block_during_shift_enabled();
            if free_clean_code {
                b.free_vec(shift_clean_code);
            }
            if apply_k2
                && dialog_gcd_apply_fused_fold_enabled()
                && std::env::var("DIALOG_GCD_FUSE_HALVE_OFF").ok().as_deref() != Some("1")
            {
                // Fuse mod_halve_inplace_fast + cmod_halve_inplace_lazy into a
                // single shared borrow chain (exact inverse of the fused double;
                // see fn doc on dialog_gcd_fused_halve_y).
                dialog_gcd_fused_halve_y_at_step(b, y, p, s2, Some(step));
            } else {
                mod_halve_inplace_fast(b, y, p);
                if apply_k2 {
                    // mirror the forward K=2 second shift: conditional 2nd halve of y.
                    // MUST use the lazy (Solinas, truncated) controlled halve to match.
                    cmod_halve_inplace_lazy(b, y, p, s2);
                }
            }
            if free_clean_code {
                b.reacquire_vec(shift_clean_code);
            }
            if !scale_released_code.is_empty() {
                b.set_phase("dialog_gcd_compressed_block_apply_reverse_scale_reacquire");
                b.reacquire_vec(scale_released_code);
            }
            if let Some(frame) = step_frame.as_ref() {
                b.set_phase("dialog_gcd_compressed_block_apply_reverse_clear_step_frame");
                dialog_gcd_k5_head11_toggle_step_frame(b, compressed_block, frame, slot);
            }
        }

        if !released_code.is_empty() {
            b.set_phase("dialog_gcd_compressed_block_apply_reverse_reacquire_block");
            b.reacquire_vec(released_code);
        }
        if head11_block && !head11_step_frame {
            b.reacquire(raw_block[1]);
        }
        b.set_phase("dialog_gcd_compressed_block_apply_reverse_clear_block_copy");
        if let Some(raw0) = inplace_raw0 {
            dialog_gcd_k2_pair_inplace_clear_block(b, compressed_block, raw0, end - start);
        } else if !head11_step_frame {
            dialog_gcd_clear_raw_block_copy(b, compressed_block, raw_block, end - start);
        }
        if let Some(frame) = step_frame {
            b.free_vec(&frame);
        }
        if head11_step_frame {
            b.set_phase("dialog_gcd_compressed_block_apply_reverse_reacquire_raw_host");
            b.reacquire_vec(raw_block);
        }
    }

    if let Some(raw0) = inplace_raw0 {
        b.free(raw0);
    }
}

pub(crate) fn emit_dialog_gcd_compressed_sidecar_tobitvector_steps(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    compressed_log: &[QubitId],
    pair: &[QubitId],
    scratch: QubitId,
) {
    assert_eq!(u.len(), N);
    assert_eq!(v.len(), N);
    assert_eq!(pair.len(), 2);
    assert!(compressed_log.len() >= dialog_gcd_compressed_sidecar_bits());

    for step in 0..dialog_gcd_active_iterations() {
        let b0 = pair[0];
        let b0_and_b1 = pair[1];
        let cmp = b.alloc_qubit();
        let active_width = dialog_gcd_tobitvector_active_width(step);
        let u_active = &u[..active_width];
        let v_active = &v[..active_width];
        let compare_bits = dialog_gcd_compare_bits_for_step(step, active_width);

        b.set_phase("dialog_gcd_compressed_sidecar_tobitvector_branch_bits");
        b.cx(v[0], b0);
        if dialog_gcd_fused_branch_bits_enabled() {
            dialog_gcd_ccx_cmp_gt_truncated_into_width(
                b,
                u_active,
                v_active,
                b0,
                b0_and_b1,
                compare_bits,
            );
        } else {
            dialog_gcd_cmp_gt_truncated_into_width(b, u_active, v_active, cmp, compare_bits);
            b.ccx(b0, cmp, b0_and_b1);
            dialog_gcd_cmp_gt_truncated_into_width(b, u_active, v_active, cmp, compare_bits);
        }
        b.free(cmp);

        b.set_phase("dialog_gcd_compressed_sidecar_tobitvector_cswap");
        let cswap_width = dialog_gcd_tobitvector_cswap_width(active_width, step);
        for (i, (&ui, &vi)) in u[..cswap_width]
            .iter()
            .zip(v[..cswap_width].iter())
            .enumerate()
        {
            if i == 0 && dialog_gcd_odd_u_lowbit_fastpath_enabled() {
                continue;
            }
            cswap(b, b0_and_b1, ui, vi);
        }

        b.set_phase("dialog_gcd_compressed_sidecar_tobitvector_subtract");
        let borrowed_carries =
            dialog_gcd_compressed_sidecar_future_carry_slice(compressed_log, step, active_width);
        dialog_gcd_controlled_sub_selected(b, u_active, v_active, b0, borrowed_carries, step);

        b.set_phase("dialog_gcd_compressed_sidecar_tobitvector_shift");
        let shift_width = dialog_gcd_tobitvector_shift_width(active_width, step);
        dialog_gcd_shift_right_assuming_even(b, &v[..shift_width]);

        b.set_phase("dialog_gcd_compressed_sidecar_tobitvector_absorb_pair");
        let block = dialog_gcd_compressed_sidecar_block(compressed_log, step);
        emit_dialog_gcd_round763_compressed_block_swapper(
            b,
            pair,
            block,
            scratch,
            step % DIALOG_GCD_HIGH_TAIL_ALIAS_GROUP_SIZE,
        );
    }
}

pub(crate) fn emit_dialog_gcd_compressed_sidecar_tobitvector_steps_reverse(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    compressed_log: &[QubitId],
    pair: &[QubitId],
    scratch: QubitId,
) {
    assert_eq!(u.len(), N);
    assert_eq!(v.len(), N);
    assert_eq!(pair.len(), 2);
    assert!(compressed_log.len() >= dialog_gcd_compressed_sidecar_bits());

    for step in (0..dialog_gcd_active_iterations()).rev() {
        b.set_phase("dialog_gcd_compressed_sidecar_tobitvector_reverse_load_pair");
        let block = dialog_gcd_compressed_sidecar_block(compressed_log, step);
        emit_dialog_gcd_round763_compressed_block_swapper(
            b,
            pair,
            block,
            scratch,
            step % DIALOG_GCD_HIGH_TAIL_ALIAS_GROUP_SIZE,
        );

        let b0 = pair[0];
        let b0_and_b1 = pair[1];
        let cmp = b.alloc_qubit();
        let active_width = dialog_gcd_tobitvector_active_width(step);
        let u_active = &u[..active_width];
        let v_active = &v[..active_width];
        let compare_bits = dialog_gcd_compare_bits_for_step(step, active_width);

        b.set_phase("dialog_gcd_compressed_sidecar_tobitvector_reverse_unshift");
        let shift_width = dialog_gcd_tobitvector_shift_width(active_width, step);
        dialog_gcd_unshift_right_assuming_even(b, &v[..shift_width]);

        b.set_phase("dialog_gcd_compressed_sidecar_tobitvector_reverse_add");
        let borrowed_carries =
            dialog_gcd_compressed_sidecar_future_carry_slice(compressed_log, step, active_width);
        dialog_gcd_controlled_add_selected(b, u_active, v_active, b0, borrowed_carries, step);

        b.set_phase("dialog_gcd_compressed_sidecar_tobitvector_reverse_cswap");
        let cswap_width = dialog_gcd_tobitvector_cswap_width(active_width, step);
        for (i, (&ui, &vi)) in u[..cswap_width]
            .iter()
            .zip(v[..cswap_width].iter())
            .enumerate()
        {
            if i == 0 && dialog_gcd_odd_u_lowbit_fastpath_enabled() {
                continue;
            }
            cswap(b, b0_and_b1, ui, vi);
        }

        b.set_phase("dialog_gcd_compressed_sidecar_tobitvector_reverse_branch_bits");
        if dialog_gcd_fused_branch_bits_enabled() {
            dialog_gcd_ccx_cmp_gt_truncated_into_width(
                b,
                u_active,
                v_active,
                b0,
                b0_and_b1,
                compare_bits,
            );
        } else {
            dialog_gcd_cmp_gt_truncated_into_width(b, u_active, v_active, cmp, compare_bits);
            b.ccx(b0, cmp, b0_and_b1);
            dialog_gcd_cmp_gt_truncated_into_width(b, u_active, v_active, cmp, compare_bits);
        }
        b.free(cmp);
        b.cx(v[0], b0);
    }
}

pub(crate) fn emit_dialog_gcd_compressed_sidecar_apply_bitvector(
    b: &mut B,
    compressed_log: &[QubitId],
    x: &[QubitId],
    y: &[QubitId],
    p: U256,
    pair: &[QubitId],
    scratch: QubitId,
) {
    assert_eq!(x.len(), N);
    assert_eq!(y.len(), N);
    assert_eq!(pair.len(), 2);

    for step in (0..dialog_gcd_active_iterations()).rev() {
        b.set_phase("dialog_gcd_compressed_sidecar_apply_load_pair");
        let block = dialog_gcd_compressed_sidecar_block(compressed_log, step);
        emit_dialog_gcd_round763_compressed_block_swapper(
            b,
            pair,
            block,
            scratch,
            step % DIALOG_GCD_HIGH_TAIL_ALIAS_GROUP_SIZE,
        );

        let b0 = pair[0];
        let b0_and_b1 = pair[1];

        b.set_phase("dialog_gcd_compressed_sidecar_apply_double_y");
        mod_double_inplace_fast(b, y, p);

        b.set_phase("dialog_gcd_compressed_sidecar_apply_cadd");
        if dialog_gcd_raw_apply_materialized_special_add_enabled() {
            dialog_gcd_cmod_add_materialized_pseudomersenne_at_step(b, y, x, b0, p, Some(step));
        } else if dialog_gcd_raw_apply_direct_special_add_enabled() {
            dialog_gcd_cmod_add_pseudomersenne_lowq(b, y, x, b0, p);
        } else {
            cmod_add_qq_lowq(b, y, x, b0, p);
        }

        b.set_phase("dialog_gcd_compressed_sidecar_apply_cswap");
        for (&xi, &yi) in x.iter().zip(y.iter()) {
            cswap(b, b0_and_b1, xi, yi);
        }

        b.set_phase("dialog_gcd_compressed_sidecar_apply_unload_pair");
        let block = dialog_gcd_compressed_sidecar_block(compressed_log, step);
        emit_dialog_gcd_round763_compressed_block_swapper(
            b,
            pair,
            block,
            scratch,
            step % DIALOG_GCD_HIGH_TAIL_ALIAS_GROUP_SIZE,
        );
    }
}

pub(crate) fn emit_dialog_gcd_compressed_sidecar_apply_bitvector_reverse_exact(
    b: &mut B,
    compressed_log: &[QubitId],
    x: &[QubitId],
    y: &[QubitId],
    p: U256,
    pair: &[QubitId],
    scratch: QubitId,
) {
    assert_eq!(x.len(), N);
    assert_eq!(y.len(), N);
    assert_eq!(pair.len(), 2);

    for step in 0..dialog_gcd_active_iterations() {
        b.set_phase("dialog_gcd_compressed_sidecar_apply_reverse_load_pair");
        let block = dialog_gcd_compressed_sidecar_block(compressed_log, step);
        emit_dialog_gcd_round763_compressed_block_swapper(
            b,
            pair,
            block,
            scratch,
            step % DIALOG_GCD_HIGH_TAIL_ALIAS_GROUP_SIZE,
        );

        let b0 = pair[0];
        let b0_and_b1 = pair[1];

        b.set_phase("dialog_gcd_compressed_sidecar_apply_reverse_cswap");
        for (&xi, &yi) in x.iter().zip(y.iter()) {
            cswap(b, b0_and_b1, xi, yi);
        }

        b.set_phase("dialog_gcd_compressed_sidecar_apply_reverse_csub");
        if dialog_gcd_raw_apply_reverse_materialized_special_sub_enabled() {
            dialog_gcd_cmod_sub_materialized_pseudomersenne_at_step(b, y, x, b0, p, Some(step));
        } else if dialog_gcd_raw_apply_reverse_fast_sub_enabled() {
            cmod_sub_qq(b, y, x, b0, p);
        } else {
            cmod_sub_qq_lowq(b, y, x, b0, p);
        }

        b.set_phase("dialog_gcd_compressed_sidecar_apply_reverse_halve_y");
        mod_halve_inplace_fast(b, y, p);

        b.set_phase("dialog_gcd_compressed_sidecar_apply_reverse_unload_pair");
        let block = dialog_gcd_compressed_sidecar_block(compressed_log, step);
        emit_dialog_gcd_round763_compressed_block_swapper(
            b,
            pair,
            block,
            scratch,
            step % DIALOG_GCD_HIGH_TAIL_ALIAS_GROUP_SIZE,
        );
    }
}

pub(crate) fn emit_dialog_gcd_compressed_sidecar_ipmul_block_lifecycle(
    b: &mut B,
    factor: &[QubitId],
    target: &[QubitId],
    p: U256,
) {
    assert_eq!(factor.len(), N);
    assert_eq!(target.len(), N);

    let compressed_log = b.alloc_qubits(dialog_gcd_allocated_compressed_sidecar_bits());
    let raw_block = if dialog_gcd_host_reverse_raw_block_enabled() {
        Vec::new()
    } else {
        b.alloc_qubits(dialog_gcd_raw_block_len())
    };
    let u = b.alloc_qubits(N);
    let runway = dialog_gcd_build_compressed_log_u_high_runway(&u, &compressed_log);
    let replay_log = runway
        .as_ref()
        .map_or(compressed_log.as_slice(), |r| r.remapped_log.as_slice());
    b.set_phase("dialog_gcd_compressed_block_ipmul_load_p");
    for i in 0..N {
        if bit(p, i) {
            b.x(u[i]);
        }
    }

    b.set_phase("dialog_gcd_compressed_block_ipmul_tobitvector");
    emit_dialog_gcd_compressed_sidecar_tobitvector_steps_block_lifecycle(
        b, &u, factor, replay_log, &raw_block,
    );

    if dialog_gcd_raw_ipmul_terminal_reuse_enabled() {
        b.set_phase("dialog_gcd_compressed_block_ipmul_release_terminal_u");
        b.x(u[0]);
        dialog_gcd_release_terminal_u(b, &u, runway.as_ref());

        b.set_phase("dialog_gcd_compressed_block_ipmul_apply_bitvector_reuse_factor_zero");
        let inplace_apply_raw = dialog_gcd_k2_apply_inplace_raw_block_enabled();
        if inplace_apply_raw && !raw_block.is_empty() {
            b.free_vec(&raw_block);
        }
        let apply_raw_block = if !inplace_apply_raw && dialog_gcd_host_reverse_raw_block_enabled() {
            b.alloc_qubits(dialog_gcd_raw_block_len())
        } else {
            Vec::new()
        };
        emit_dialog_gcd_compressed_sidecar_apply_bitvector_block_lifecycle(
            b,
            replay_log,
            target,
            factor,
            p,
            if inplace_apply_raw {
                &[]
            } else if apply_raw_block.is_empty() {
                &raw_block
            } else {
                &apply_raw_block
            },
        );
        if !apply_raw_block.is_empty() {
            b.free_vec(&apply_raw_block);
        }

        if dialog_gcd_raw_ipmul_clear_p_residual_enabled() {
            b.set_phase("dialog_gcd_compressed_block_ipmul_clear_p_residual_source_lane");
            for i in 0..N {
                if bit(p, i) {
                    b.x(target[i]);
                }
            }
        }

        b.set_phase("dialog_gcd_compressed_block_ipmul_swap_product_into_target");
        for i in 0..N {
            b.swap(target[i], factor[i]);
        }

        if inplace_apply_raw && !raw_block.is_empty() {
            b.reacquire_vec(&raw_block);
        }

        b.set_phase("dialog_gcd_compressed_block_ipmul_reacquire_terminal_u");
        dialog_gcd_reacquire_terminal_u(b, &u, runway.as_ref());
        b.set_phase("dialog_gcd_compressed_block_ipmul_seed_terminal_u");
        b.x(u[0]);

        b.set_phase("dialog_gcd_compressed_block_ipmul_uncompute_tobitvector");
        emit_dialog_gcd_compressed_sidecar_tobitvector_steps_reverse_block_lifecycle(
            b, &u, factor, replay_log, &raw_block,
        );

        b.set_phase("dialog_gcd_compressed_block_ipmul_unload_p");
        for i in 0..N {
            if bit(p, i) {
                b.x(u[i]);
            }
        }
        if !b.k2_shift2_log.is_empty() {
            let log = std::mem::take(&mut b.k2_shift2_log);
            b.free_vec(&log);
        }
        b.free_vec(&u);
        if !raw_block.is_empty() {
            b.free_vec(&raw_block);
        }
        b.free_vec(&compressed_log);
        return;
    }

    let tmp = b.alloc_qubits(N);
    b.set_phase("dialog_gcd_compressed_block_ipmul_apply_bitvector");
    emit_dialog_gcd_compressed_sidecar_apply_bitvector_block_lifecycle(
        b, replay_log, target, &tmp, p, &raw_block,
    );

    b.set_phase("dialog_gcd_compressed_block_ipmul_swap_product_into_target");
    for i in 0..N {
        b.swap(target[i], tmp[i]);
    }

    b.set_phase("dialog_gcd_compressed_block_ipmul_free_zero_tmp");
    b.free_vec(&tmp);

    b.set_phase("dialog_gcd_compressed_block_ipmul_uncompute_tobitvector");
    emit_dialog_gcd_compressed_sidecar_tobitvector_steps_reverse_block_lifecycle(
        b, &u, factor, replay_log, &raw_block,
    );

    b.set_phase("dialog_gcd_compressed_block_ipmul_unload_p");
    for i in 0..N {
        if bit(p, i) {
            b.x(u[i]);
        }
    }
    if !b.k2_shift2_log.is_empty() {
        let log = std::mem::take(&mut b.k2_shift2_log);
        b.free_vec(&log);
    }
    b.free_vec(&u);
    b.free_vec(&raw_block);
    b.free_vec(&compressed_log);
}

pub(crate) fn emit_dialog_gcd_compressed_sidecar_ipmul(
    b: &mut B,
    factor: &[QubitId],
    target: &[QubitId],
    p: U256,
) {
    assert_eq!(factor.len(), N);
    assert_eq!(target.len(), N);

    if dialog_gcd_compressed_block_lifecycle_enabled() {
        emit_dialog_gcd_compressed_sidecar_ipmul_block_lifecycle(b, factor, target, p);
        return;
    }

    let compressed_log = b.alloc_qubits(dialog_gcd_compressed_sidecar_bits());
    let pair = b.alloc_qubits(2);
    let compressor_scratch = b.alloc_qubit();
    let u = b.alloc_qubits(N);
    b.set_phase("dialog_gcd_compressed_sidecar_ipmul_load_p");
    for i in 0..N {
        if bit(p, i) {
            b.x(u[i]);
        }
    }

    b.set_phase("dialog_gcd_compressed_sidecar_ipmul_tobitvector");
    emit_dialog_gcd_compressed_sidecar_tobitvector_steps(
        b,
        &u,
        factor,
        &compressed_log,
        &pair,
        compressor_scratch,
    );

    if dialog_gcd_raw_ipmul_terminal_reuse_enabled() {
        b.set_phase("dialog_gcd_compressed_sidecar_ipmul_release_terminal_u");
        b.x(u[0]);
        b.free_vec(&u);

        b.set_phase("dialog_gcd_compressed_sidecar_ipmul_apply_bitvector_reuse_factor_zero");
        emit_dialog_gcd_compressed_sidecar_apply_bitvector(
            b,
            &compressed_log,
            target,
            factor,
            p,
            &pair,
            compressor_scratch,
        );

        if dialog_gcd_raw_ipmul_clear_p_residual_enabled() {
            b.set_phase("dialog_gcd_compressed_sidecar_ipmul_clear_p_residual_source_lane");
            for i in 0..N {
                if bit(p, i) {
                    b.x(target[i]);
                }
            }
        }

        b.set_phase("dialog_gcd_compressed_sidecar_ipmul_swap_product_into_target");
        for i in 0..N {
            b.swap(target[i], factor[i]);
        }

        b.set_phase("dialog_gcd_compressed_sidecar_ipmul_reacquire_terminal_u");
        b.reacquire_vec(&u);
        b.set_phase("dialog_gcd_compressed_sidecar_ipmul_seed_terminal_u");
        b.x(u[0]);

        b.set_phase("dialog_gcd_compressed_sidecar_ipmul_uncompute_tobitvector");
        emit_dialog_gcd_compressed_sidecar_tobitvector_steps_reverse(
            b,
            &u,
            factor,
            &compressed_log,
            &pair,
            compressor_scratch,
        );

        b.set_phase("dialog_gcd_compressed_sidecar_ipmul_unload_p");
        for i in 0..N {
            if bit(p, i) {
                b.x(u[i]);
            }
        }
        b.free_vec(&u);
        b.free(compressor_scratch);
        b.free_vec(&pair);
        b.free_vec(&compressed_log);
        return;
    }

    let tmp = b.alloc_qubits(N);
    b.set_phase("dialog_gcd_compressed_sidecar_ipmul_apply_bitvector");
    emit_dialog_gcd_compressed_sidecar_apply_bitvector(
        b,
        &compressed_log,
        target,
        &tmp,
        p,
        &pair,
        compressor_scratch,
    );

    b.set_phase("dialog_gcd_compressed_sidecar_ipmul_swap_product_into_target");
    for i in 0..N {
        b.swap(target[i], tmp[i]);
    }

    b.set_phase("dialog_gcd_compressed_sidecar_ipmul_free_zero_tmp");
    b.free_vec(&tmp);

    b.set_phase("dialog_gcd_compressed_sidecar_ipmul_uncompute_tobitvector");
    emit_dialog_gcd_compressed_sidecar_tobitvector_steps_reverse(
        b,
        &u,
        factor,
        &compressed_log,
        &pair,
        compressor_scratch,
    );

    b.set_phase("dialog_gcd_compressed_sidecar_ipmul_unload_p");
    for i in 0..N {
        if bit(p, i) {
            b.x(u[i]);
        }
    }
    b.free_vec(&u);
    b.free(compressor_scratch);
    b.free_vec(&pair);
    b.free_vec(&compressed_log);
}

pub(crate) fn emit_dialog_gcd_compressed_sidecar_quotient_block_lifecycle(
    b: &mut B,
    factor: &[QubitId],
    target: &[QubitId],
    p: U256,
) {
    assert_eq!(factor.len(), N);
    assert_eq!(target.len(), N);

    let compressed_log = b.alloc_qubits(dialog_gcd_allocated_compressed_sidecar_bits());
    let raw_block = if dialog_gcd_host_reverse_raw_block_enabled() {
        Vec::new()
    } else {
        b.alloc_qubits(dialog_gcd_raw_block_len())
    };
    let u = b.alloc_qubits(N);
    let runway = dialog_gcd_build_compressed_log_u_high_runway(&u, &compressed_log);
    let replay_log = runway
        .as_ref()
        .map_or(compressed_log.as_slice(), |r| r.remapped_log.as_slice());
    b.set_phase("dialog_gcd_compressed_block_quotient_load_p");
    for i in 0..N {
        if bit(p, i) {
            b.x(u[i]);
        }
    }

    b.set_phase("dialog_gcd_compressed_block_quotient_tobitvector");
    emit_dialog_gcd_compressed_sidecar_tobitvector_steps_block_lifecycle(
        b, &u, factor, replay_log, &raw_block,
    );

    if dialog_gcd_raw_quotient_terminal_reuse_enabled() {
        b.set_phase("dialog_gcd_compressed_block_quotient_release_terminal_u");
        b.x(u[0]);
        dialog_gcd_release_terminal_u(b, &u, runway.as_ref());

        b.set_phase("dialog_gcd_compressed_block_quotient_apply_reverse_reuse_factor_zero");
        let inplace_apply_raw = dialog_gcd_k2_apply_inplace_raw_block_enabled();
        if inplace_apply_raw && !raw_block.is_empty() {
            b.free_vec(&raw_block);
        }
        let apply_raw_block = if !inplace_apply_raw && dialog_gcd_host_reverse_raw_block_enabled() {
            b.alloc_qubits(dialog_gcd_raw_block_len())
        } else {
            Vec::new()
        };
        emit_dialog_gcd_compressed_sidecar_apply_bitvector_reverse_exact_block_lifecycle(
            b,
            replay_log,
            factor,
            target,
            p,
            if inplace_apply_raw {
                &[]
            } else if apply_raw_block.is_empty() {
                &raw_block
            } else {
                &apply_raw_block
            },
        );
        if !apply_raw_block.is_empty() {
            b.free_vec(&apply_raw_block);
        }

        b.set_phase("dialog_gcd_compressed_block_quotient_swap_quotient_into_target");
        for i in 0..N {
            b.swap(target[i], factor[i]);
        }

        if inplace_apply_raw && !raw_block.is_empty() {
            b.reacquire_vec(&raw_block);
        }

        b.set_phase("dialog_gcd_compressed_block_quotient_reacquire_terminal_u");
        dialog_gcd_reacquire_terminal_u(b, &u, runway.as_ref());
        b.set_phase("dialog_gcd_compressed_block_quotient_seed_terminal_u");
        b.x(u[0]);

        b.set_phase("dialog_gcd_compressed_block_quotient_uncompute_tobitvector");
        emit_dialog_gcd_compressed_sidecar_tobitvector_steps_reverse_block_lifecycle(
            b, &u, factor, replay_log, &raw_block,
        );

        b.set_phase("dialog_gcd_compressed_block_quotient_unload_p");
        for i in 0..N {
            if bit(p, i) {
                b.x(u[i]);
            }
        }
        if !b.k2_shift2_log.is_empty() {
            let log = std::mem::take(&mut b.k2_shift2_log);
            b.free_vec(&log);
        }
        b.free_vec(&u);
        if !raw_block.is_empty() {
            b.free_vec(&raw_block);
        }
        b.free_vec(&compressed_log);
        return;
    }

    b.set_phase("dialog_gcd_compressed_block_quotient_apply_reverse");
    emit_dialog_gcd_compressed_sidecar_apply_bitvector_reverse_exact_block_lifecycle(
        b, replay_log, factor, target, p, &raw_block,
    );

    b.set_phase("dialog_gcd_compressed_block_quotient_uncompute_tobitvector");
    emit_dialog_gcd_compressed_sidecar_tobitvector_steps_reverse_block_lifecycle(
        b, &u, factor, replay_log, &raw_block,
    );

    b.set_phase("dialog_gcd_compressed_block_quotient_unload_p");
    for i in 0..N {
        if bit(p, i) {
            b.x(u[i]);
        }
    }
    if !b.k2_shift2_log.is_empty() {
        let log = std::mem::take(&mut b.k2_shift2_log);
        b.free_vec(&log);
    }
    b.free_vec(&u);
    b.free_vec(&raw_block);
    b.free_vec(&compressed_log);
}

pub(crate) fn emit_dialog_gcd_compressed_sidecar_quotient(
    b: &mut B,
    factor: &[QubitId],
    target: &[QubitId],
    p: U256,
) {
    assert_eq!(factor.len(), N);
    assert_eq!(target.len(), N);

    if dialog_gcd_compressed_block_lifecycle_enabled() {
        emit_dialog_gcd_compressed_sidecar_quotient_block_lifecycle(b, factor, target, p);
        return;
    }

    let compressed_log = b.alloc_qubits(dialog_gcd_compressed_sidecar_bits());
    let pair = b.alloc_qubits(2);
    let compressor_scratch = b.alloc_qubit();
    let u = b.alloc_qubits(N);
    b.set_phase("dialog_gcd_compressed_sidecar_quotient_load_p");
    for i in 0..N {
        if bit(p, i) {
            b.x(u[i]);
        }
    }

    b.set_phase("dialog_gcd_compressed_sidecar_quotient_tobitvector");
    emit_dialog_gcd_compressed_sidecar_tobitvector_steps(
        b,
        &u,
        factor,
        &compressed_log,
        &pair,
        compressor_scratch,
    );

    if dialog_gcd_raw_quotient_terminal_reuse_enabled() {
        b.set_phase("dialog_gcd_compressed_sidecar_quotient_release_terminal_u");
        b.x(u[0]);
        b.free_vec(&u);

        b.set_phase("dialog_gcd_compressed_sidecar_quotient_apply_reverse_reuse_factor_zero");
        emit_dialog_gcd_compressed_sidecar_apply_bitvector_reverse_exact(
            b,
            &compressed_log,
            factor,
            target,
            p,
            &pair,
            compressor_scratch,
        );

        b.set_phase("dialog_gcd_compressed_sidecar_quotient_swap_quotient_into_target");
        for i in 0..N {
            b.swap(target[i], factor[i]);
        }

        b.set_phase("dialog_gcd_compressed_sidecar_quotient_reacquire_terminal_u");
        b.reacquire_vec(&u);
        b.set_phase("dialog_gcd_compressed_sidecar_quotient_seed_terminal_u");
        b.x(u[0]);

        b.set_phase("dialog_gcd_compressed_sidecar_quotient_uncompute_tobitvector");
        emit_dialog_gcd_compressed_sidecar_tobitvector_steps_reverse(
            b,
            &u,
            factor,
            &compressed_log,
            &pair,
            compressor_scratch,
        );

        b.set_phase("dialog_gcd_compressed_sidecar_quotient_unload_p");
        for i in 0..N {
            if bit(p, i) {
                b.x(u[i]);
            }
        }
        b.free_vec(&u);
        b.free(compressor_scratch);
        b.free_vec(&pair);
        b.free_vec(&compressed_log);
        return;
    }

    b.set_phase("dialog_gcd_compressed_sidecar_quotient_apply_reverse");
    emit_dialog_gcd_compressed_sidecar_apply_bitvector_reverse_exact(
        b,
        &compressed_log,
        factor,
        target,
        p,
        &pair,
        compressor_scratch,
    );

    b.set_phase("dialog_gcd_compressed_sidecar_quotient_uncompute_tobitvector");
    emit_dialog_gcd_compressed_sidecar_tobitvector_steps_reverse(
        b,
        &u,
        factor,
        &compressed_log,
        &pair,
        compressor_scratch,
    );

    b.set_phase("dialog_gcd_compressed_sidecar_quotient_unload_p");
    for i in 0..N {
        if bit(p, i) {
            b.x(u[i]);
        }
    }
    b.free_vec(&u);
    b.free(compressor_scratch);
    b.free_vec(&pair);
    b.free_vec(&compressed_log);
}


pub(crate) fn emit_dialog_gcd_k2_pair_core_encoder(b: &mut B, core: &[QubitId]) {
    assert_eq!(core.len(), 5);
    // PROBE: 3-CCX reachable-support encoder (replaces 6-CCX). −3 scored CCX/call.
    b.cx(core[1], core[2]);
    b.cx(core[0], core[4]);
    b.x(core[3]);
    b.ccx(core[2], core[3], core[1]);
    b.cx(core[3], core[4]);
    b.ccx(core[3], core[4], core[0]);
    b.cx(core[2], core[4]);
    b.cx(core[0], core[3]);
    b.cx(core[3], core[2]);
    b.cx(core[3], core[4]);
    b.ccx(core[1], core[3], core[0]);
    b.cx(core[1], core[0]);
    b.cx(core[3], core[0]);
}

pub(crate) fn emit_dialog_gcd_k2_pair_core_encoder_inverse(b: &mut B, core: &[QubitId]) {
    assert_eq!(core.len(), 5);
    // Exact gate-reverse of the 3-CCX encoder (each op self-inverse).
    b.cx(core[3], core[0]);
    b.cx(core[1], core[0]);
    b.ccx(core[1], core[3], core[0]);
    b.cx(core[3], core[4]);
    b.cx(core[3], core[2]);
    b.cx(core[0], core[3]);
    b.cx(core[2], core[4]);
    b.ccx(core[3], core[4], core[0]);
    b.cx(core[3], core[4]);
    b.ccx(core[2], core[3], core[1]);
    b.x(core[3]);
    b.cx(core[0], core[4]);
    b.cx(core[1], core[2]);
}

pub(crate) fn dialog_gcd_k2_pair_core(raw_block: &[QubitId]) -> [QubitId; 5] {
    assert_eq!(raw_block.len(), 6);
    [
        raw_block[0], // first step b0
        raw_block[1], // first step b0_and_b1
        raw_block[4], // first step shift2
        raw_block[2], // second step b0
        raw_block[3], // second step b0_and_b1
    ]
}

pub(crate) fn dialog_gcd_k2_pair_copy_compressed_block_to_raw(
    b: &mut B,
    compressed_block: &[QubitId],
    raw_block: &[QubitId],
    steps: usize,
) {
    assert_eq!(compressed_block.len(), DIALOG_GCD_HIGH_TAIL_ALIAS_BLOCK_BITS);
    assert_eq!(raw_block.len(), 6);
    assert!((1..=2).contains(&steps));
    let swap_host = dialog_gcd_apply_replay_swap_host_enabled();
    if steps == 1 {
        let raw_encoded = [raw_block[0], raw_block[1], raw_block[4]];
        for (&c, &r) in compressed_block.iter().take(3).zip(raw_encoded.iter()) {
            if swap_host {
                b.swap(c, r);
            } else {
                b.cx(c, r);
            }
        }
        return;
    }
    let raw_encoded = [raw_block[1], raw_block[4], raw_block[2], raw_block[3], raw_block[5]];
    for (&c, &r) in compressed_block.iter().zip(raw_encoded.iter()) {
        if swap_host {
            b.swap(c, r);
        } else {
            b.cx(c, r);
        }
    }
    let core = dialog_gcd_k2_pair_core(raw_block);
    emit_dialog_gcd_k2_pair_core_encoder_inverse(b, &core);
}

pub(crate) fn dialog_gcd_k2_pair_clear_raw_block_copy(
    b: &mut B,
    compressed_block: &[QubitId],
    raw_block: &[QubitId],
    steps: usize,
) {
    assert_eq!(compressed_block.len(), DIALOG_GCD_HIGH_TAIL_ALIAS_BLOCK_BITS);
    assert_eq!(raw_block.len(), 6);
    assert!((1..=2).contains(&steps));
    let swap_host = dialog_gcd_apply_replay_swap_host_enabled();
    if steps == 1 {
        let raw_encoded = [raw_block[0], raw_block[1], raw_block[4]];
        for (&c, &r) in compressed_block.iter().take(3).zip(raw_encoded.iter()) {
            if swap_host {
                b.swap(c, r);
            } else {
                b.cx(c, r);
            }
        }
        return;
    }
    let core = dialog_gcd_k2_pair_core(raw_block);
    emit_dialog_gcd_k2_pair_core_encoder(b, &core);
    let raw_encoded = [raw_block[1], raw_block[4], raw_block[2], raw_block[3], raw_block[5]];
    for (&c, &r) in compressed_block.iter().zip(raw_encoded.iter()) {
        if swap_host {
            b.swap(c, r);
        } else {
            b.cx(c, r);
        }
    }
}

fn dialog_gcd_fixed_twice_fold(
    b: &mut B,
    y: &[QubitId],
    p: U256,
    e: QubitId,
    d: QubitId,
    is_add: bool,
) {
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));
    let h = b.alloc_qubit();
    b.ccx(e, d, h);
    let xed = b.alloc_qubit();
    b.cx(e, xed);
    b.cx(d, xed);
    let eord = b.alloc_qubit();
    b.cx(xed, eord);
    b.cx(h, eord);
    let n10 = b.alloc_qubit();
    b.cx(d, n10);
    b.cx(h, n10);

    let hi_c = highest_set_bit(c);
    let hi_delta = hi_c + 1;
    let controls = secp_fold_controls(e, d, h, xed, eord, n10, hi_delta, hi_c);
    let last = match fold_only_carry_trunc_window().or_else(double_carry_trunc_window) {
        Some(w) => core::cmp::min(y.len() - 2, hi_delta.saturating_add(w)),
        None => y.len() - 2,
    };
    if fold_freed_tail_enabled() && last > hi_delta {
        fold_ripple_freed_tail(b, y, e, d, h, xed, eord, n10, last, is_add);
    } else if is_add {
        cadd_per_position_controls_trunc(b, y, &controls, last);
    } else {
        csub_per_position_controls_trunc(b, y, &controls, last);
    }

    b.cx(h, n10);
    b.cx(d, n10);
    b.cx(h, eord);
    b.cx(xed, eord);
    b.cx(d, xed);
    b.cx(e, xed);
    b.free(n10);
    b.free(eord);
    b.free(xed);
    if dialog_gcd_fused_hclear_measured_enabled() {
        let measured = b.alloc_bit();
        b.hmr(h, measured);
        b.cz_if(e, d, measured);
    } else {
        b.ccx(e, d, h);
    }
    b.free(h);
}

fn dialog_gcd_fixed_double_twice_y(b: &mut B, y: &[QubitId], p: U256) {
    let n = y.len();
    debug_assert_eq!(n, 256);
    let ovf1 = b.alloc_qubit();
    b.swap(y[n - 1], ovf1);
    for i in (0..n - 1).rev() {
        b.swap(y[i], y[i + 1]);
    }
    let ovf2 = b.alloc_qubit();
    b.swap(y[n - 1], ovf2);
    for i in (0..n - 1).rev() {
        b.swap(y[i], y[i + 1]);
    }

    dialog_gcd_fixed_twice_fold(b, y, p, ovf2, ovf1, true);
    b.cx(y[0], ovf2);
    b.cx(y[1], ovf1);
    b.free(ovf2);
    b.free(ovf1);
}

fn dialog_gcd_fixed_halve_twice_y(b: &mut B, y: &[QubitId], p: U256) {
    let n = y.len();
    debug_assert_eq!(n, 256);
    let ovf2 = b.alloc_qubit();
    let ovf1 = b.alloc_qubit();
    b.cx(y[0], ovf2);
    b.cx(y[1], ovf1);

    dialog_gcd_fixed_twice_fold(b, y, p, ovf2, ovf1, false);
    for i in 0..n - 1 {
        b.swap(y[i], y[i + 1]);
    }
    b.swap(y[n - 1], ovf2);
    b.free(ovf2);
    for i in 0..n - 1 {
        b.swap(y[i], y[i + 1]);
    }
    b.swap(y[n - 1], ovf1);
    b.free(ovf1);
}

pub(crate) fn dialog_gcd_k5_tail7_fixed_apply_selftest() -> Result<(), String> {
    use sha3::digest::{ExtendableOutput, Update};

    let input_masks = (0..N)
        .map(|bit| {
            let x = (bit as u64)
                .wrapping_mul(0x9E37_79B9_7F4A_7C15)
                .wrapping_add(0xD1B5_4A32_D192_ED03);
            let x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            x ^ (x >> 27)
        })
        .collect::<Vec<_>>();

    for reverse in [false, true] {
        let build_one = |fixed: bool| {
            let mut b = B::new();
            let y = b.alloc_qubits(N);
            if fixed {
                if reverse {
                    dialog_gcd_fixed_halve_twice_y(&mut b, &y, SECP256K1_P);
                } else {
                    dialog_gcd_fixed_double_twice_y(&mut b, &y, SECP256K1_P);
                }
            } else {
                let one = b.alloc_qubit();
                b.x(one);
                if reverse {
                    dialog_gcd_fused_halve_y(&mut b, &y, SECP256K1_P, one);
                } else {
                    dialog_gcd_fused_double_y(&mut b, &y, SECP256K1_P, one);
                }
                b.x(one);
                b.free(one);
            }
            (b.ops, y, b.next_qubit as usize, b.next_bit as usize)
        };

        let run = |fixed: bool| {
            let (ops, y, num_qubits, num_bits) = build_one(fixed);
            let mut seed = sha3::Shake128::default();
            seed.update(b"dialog-gcd-k5-tail7-fixed-apply-selftest");
            seed.update(&[reverse as u8, fixed as u8]);
            let mut xof = seed.finalize_xof();
            let mut sim = Simulator::new(num_qubits, num_bits, &mut xof);
            sim.clear_for_shot();
            for (&q, &mask) in y.iter().zip(input_masks.iter()) {
                *sim.qubit_mut(q) = mask;
            }
            sim.apply_iter(ops.iter());
            let output = y.iter().map(|&q| sim.qubit(q)).collect::<Vec<_>>();
            let clean = (N..num_qubits).all(|q| sim.qubit(QubitId(q as u64)) == 0);
            (output, clean, sim.phase)
        };

        let (baseline, baseline_clean, baseline_phase) = run(false);
        let (fixed, fixed_clean, fixed_phase) = run(true);
        if !baseline_clean || baseline_phase != 0 {
            return Err(format!(
                "baseline dirty: reverse={reverse} clean={baseline_clean} phase=0x{baseline_phase:x}"
            ));
        }
        if !fixed_clean || fixed_phase != 0 {
            return Err(format!(
                "fixed dirty: reverse={reverse} clean={fixed_clean} phase=0x{fixed_phase:x}"
            ));
        }
        if baseline != fixed {
            let bit = baseline
                .iter()
                .zip(fixed.iter())
                .position(|(a, b)| a != b)
                .expect("different vectors have a differing bit");
            return Err(format!(
                "value mismatch: reverse={reverse} bit={bit} baseline=0x{:x} fixed=0x{:x}",
                baseline[bit], fixed[bit]
            ));
        }
    }
    Ok(())
}

pub(crate) fn dialog_gcd_fused_double_y(b: &mut B, y: &[QubitId], p: U256, s2: QubitId) {
    dialog_gcd_fused_double_y_at_step(b, y, p, s2, None);
}

pub(crate) fn dialog_gcd_fused_double_y_at_step(
    b: &mut B,
    y: &[QubitId],
    p: U256,
    s2: QubitId,
    step: Option<usize>,
) {
    let n = y.len();
    debug_assert_eq!(n, 256);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));

    // ── shift1 (unconditional left shift): ovf1 = old y[255]; y[0] = 0 ──
    let ovf1 = b.alloc_qubit();
    b.swap(y[n - 1], ovf1);
    for i in (0..n - 1).rev() {
        b.swap(y[i], y[i + 1]);
    }

    // ── cond-shift2 (left shift gated by s2) on the UNFOLDED register ──
    // ovf2 = s2 & top(Y0); y[0] = 0 (and y[1] = 0 iff s2, used by cleanup).
    let ovf2 = b.alloc_qubit();
    cswap(b, s2, y[n - 1], ovf2);
    for i in (0..n - 1).rev() {
        if dialog_gcd_skip_zero_edge_apply_double_cshift_enabled() && i == 0 {
            continue;
        }
        cswap(b, s2, y[i], y[i + 1]);
    }

    // ── derive the fold controls ──
    let e = b.alloc_qubit();
    let d = b.alloc_qubit();
    let hi_delta = highest_set_bit(c) + 1; // = 33 for secp256k1
    let last = match dialog_gcd_fused_fold_carry_trunc_window(step) {
        Some(w) => core::cmp::min(n - 2, hi_delta.saturating_add(w)),
        None => n - 2,
    };
    if fold_stream_controls_enabled() && fold_freed_tail_enabled() && last > hi_delta {
        if std::env::var("DIALOG_GCD_FOLD_PROFILE_PHASES").ok().as_deref() == Some("1") {
            b.set_phase("dialog_gcd_streamed_double_setup");
        }
        b.ccx(ovf1, s2, d);
        b.cx(ovf1, e);
        b.cx(d, e);
        b.cx(ovf2, e);
        fold_ripple_freed_tail_ed_streamed(
            b,
            y,
            e,
            d,
            Some((ovf1, ovf2, s2)),
            fold_park_low_carries_at_step(step),
            last,
            true,
        );
        if std::env::var("DIALOG_GCD_FOLD_PROFILE_PHASES").ok().as_deref() == Some("1") {
            b.set_phase("dialog_gcd_streamed_double_cleanup");
        }
    } else {
        let h = b.alloc_qubit();
        b.ccx(ovf1, s2, d);
        b.cx(ovf1, e);
        b.cx(d, e);
        b.cx(ovf2, e);
        b.ccx(ovf2, d, h);
        let xed = b.alloc_qubit();
        b.cx(e, xed);
        b.cx(d, xed);
        let eord = b.alloc_qubit();
        b.cx(xed, eord);
        b.cx(h, eord);
        let n10 = b.alloc_qubit();
        b.cx(d, n10);
        b.cx(h, n10);

        let mut controls: Vec<Option<QubitId>> = vec![None; hi_delta + 1];
        controls[0] = Some(e);
        controls[1] = Some(d);
        controls[4] = Some(e);
        controls[5] = Some(d);
        controls[6] = Some(e);
        controls[7] = Some(xed);
        controls[8] = Some(eord);
        controls[9] = Some(eord);
        controls[10] = Some(n10);
        controls[11] = Some(h);
        controls[highest_set_bit(c)] = Some(e);
        controls[hi_delta] = Some(d);
        if fold_freed_tail_enabled() && last > hi_delta {
            fold_ripple_freed_tail_ed(
                b,
                y,
                e,
                d,
                h,
                xed,
                eord,
                n10,
                Some((ovf1, ovf2, s2)),
                step,
                last,
                true,
            );
        } else {
            cadd_per_position_controls_trunc(b, y, &controls, last);
        }

        b.cx(h, n10);
        b.cx(d, n10);
        b.cx(h, eord);
        b.cx(xed, eord);
        b.cx(d, xed);
        b.cx(e, xed);
        b.free(n10);
        b.free(eord);
        b.free(xed);
        if dialog_gcd_fused_hclear_measured_enabled() {
            let m = b.alloc_bit();
            b.hmr(h, m);
            b.cz_if(ovf2, d, m);
        } else {
            b.ccx(ovf2, d, h);
        }
        b.free(h);
    }

    // ── cleanup: return the base and overflow controls to |0⟩ ──
    // Clear e via parity: y[0] == e.
    b.cx(y[0], e);
    // Clear d. Stock: ccx(s2, y[1], d) (d == s2 & y[1] post-fold). Measured
    // variant: d was set as `ovf1 & s2` and neither d, ovf1, nor s2 changed
    // since (ovf1 is an untouched overflow holder, s2 is the read-only gate, d
    // is used only as a control). So a Gidney measurement-uncompute on the
    // ORIGINAL set-controls is value-identical (forces d->0) and phase-exact
    // (d·rng cancels cz_if(ovf1, s2, ·)), at 0 Toffoli instead of 1.
    if dialog_gcd_fused_dclear_measured_enabled() {
        let m = b.alloc_bit();
        b.hmr(d, m);
        b.cz_if(ovf1, s2, m);
    } else {
        b.ccx(s2, y[1], d);
    }
    b.free(d);
    b.free(e);
    // Clear ovf1 == (s2 ? y[1] : y[0]).
    if dialog_gcd_fused_ovfclear_measured_enabled() {
        let m = b.alloc_bit();
        b.hmr(ovf1, m);
        b.cz_if(s2, y[1], m);
        b.x(s2);
        b.cz_if(s2, y[0], m);
        b.x(s2);
    } else {
        b.ccx(s2, y[1], ovf1);
        b.x(s2);
        b.ccx(s2, y[0], ovf1);
        b.x(s2);
    }
    b.free(ovf1);
    // Clear ovf2 == s2 & y[0].
    if dialog_gcd_fused_ovfclear_measured_enabled() {
        let m = b.alloc_bit();
        b.hmr(ovf2, m);
        b.cz_if(s2, y[0], m);
    } else {
        b.ccx(s2, y[0], ovf2);
    }
    b.free(ovf2);
}

pub(crate) fn dialog_gcd_fused_halve_y(b: &mut B, y: &[QubitId], p: U256, s2: QubitId) {
    dialog_gcd_fused_halve_y_at_step(b, y, p, s2, None);
}

pub(crate) fn dialog_gcd_fused_halve_y_at_step(
    b: &mut B,
    y: &[QubitId],
    p: U256,
    s2: QubitId,
    step: Option<usize>,
) {
    let n = y.len();
    debug_assert_eq!(n, 256);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));

    // ── recover the base fold controls directly from y_new ──
    let e = b.alloc_qubit();
    let d = b.alloc_qubit();
    let hi_delta = highest_set_bit(c) + 1; // = 33 for secp256k1
    let last = match dialog_gcd_fused_fold_carry_trunc_window(step) {
        Some(w) => core::cmp::min(n - 2, hi_delta.saturating_add(w)),
        None => n - 2,
    };
    let (ovf2, ovf1) =
        if fold_stream_controls_enabled() && fold_freed_tail_enabled() && last > hi_delta {
            if std::env::var("DIALOG_GCD_FOLD_PROFILE_PHASES").ok().as_deref() == Some("1") {
                b.set_phase("dialog_gcd_streamed_halve_setup");
            }
            b.cx(y[0], e);
            b.ccx(s2, y[1], d);
            let ovf2 = b.alloc_qubit();
            let ovf1 = b.alloc_qubit();
            b.ccx(e, s2, ovf2);
            b.cx(e, ovf1);
            let xed = b.alloc_qubit();
            b.cx(e, xed);
            b.cx(d, xed);
            b.ccx(s2, xed, ovf1);
            b.cx(d, xed);
            b.cx(e, xed);
            b.free(xed);
            fold_ripple_freed_tail_ed_streamed(
                b,
                y,
                e,
                d,
                Some((ovf1, ovf2, s2)),
                fold_park_low_carries_at_step(step),
                last,
                false,
            );
            if std::env::var("DIALOG_GCD_FOLD_PROFILE_PHASES").ok().as_deref() == Some("1") {
                b.set_phase("dialog_gcd_streamed_halve_cleanup");
            }
            (ovf2, ovf1)
        } else {
            let h = b.alloc_qubit();
            b.cx(y[0], e);
            b.ccx(s2, y[1], d);
            b.ccx(e, d, h);
            let xed = b.alloc_qubit();
            b.cx(e, xed);
            b.cx(d, xed);
            let eord = b.alloc_qubit();
            b.cx(xed, eord);
            b.cx(h, eord);
            let n10 = b.alloc_qubit();
            b.cx(d, n10);
            b.cx(h, n10);
            let ovf2 = b.alloc_qubit();
            let ovf1 = b.alloc_qubit();
            b.ccx(e, s2, ovf2);
            b.cx(e, ovf1);
            b.ccx(s2, xed, ovf1);

            let mut controls: Vec<Option<QubitId>> = vec![None; hi_delta + 1];
            controls[0] = Some(e);
            controls[1] = Some(d);
            controls[4] = Some(e);
            controls[5] = Some(d);
            controls[6] = Some(e);
            controls[7] = Some(xed);
            controls[8] = Some(eord);
            controls[9] = Some(eord);
            controls[10] = Some(n10);
            controls[11] = Some(h);
            controls[highest_set_bit(c)] = Some(e);
            controls[hi_delta] = Some(d);
            if fold_freed_tail_enabled() && last > hi_delta {
                fold_ripple_freed_tail_ed(
                    b,
                    y,
                    e,
                    d,
                    h,
                    xed,
                    eord,
                    n10,
                    Some((ovf1, ovf2, s2)),
                    step,
                    last,
                    false,
                );
            } else {
                csub_per_position_controls_trunc(b, y, &controls, last);
            }

            b.cx(h, n10);
            b.cx(d, n10);
            b.cx(h, eord);
            b.cx(xed, eord);
            b.cx(d, xed);
            b.cx(e, xed);
            b.free(n10);
            b.free(eord);
            b.free(xed);
            if dialog_gcd_fused_hclear_measured_enabled() {
                let m = b.alloc_bit();
                b.hmr(h, m);
                b.cz_if(e, d, m);
            } else {
                b.ccx(e, d, h);
            }
            b.free(h);
            (ovf2, ovf1)
        };

    // Clear e and d via the live overflow qubits (the register low bits are now
    // cleared by the csub, so we cannot read them off y any more):
    //   e == (s2 ? ovf2 : ovf1);   d == (s2 ? ovf1 : 0).
    if dialog_gcd_fused_halve_edclear_measured_enabled() {
        let me = b.alloc_bit();
        b.hmr(e, me);
        b.x(s2);
        b.cz_if(s2, ovf1, me);
        b.x(s2);
        b.cz_if(s2, ovf2, me);
        let md = b.alloc_bit();
        b.hmr(d, md);
        b.cz_if(s2, ovf1, md);
    } else {
        b.x(s2);
        b.ccx(s2, ovf1, e); // s2=0: e ^= ovf1
        b.x(s2);
        b.ccx(s2, ovf2, e); // s2=1: e ^= ovf2
        b.ccx(s2, ovf1, d); // s2=1: d ^= ovf1   (s2=0: d already 0)
    }
    b.free(e);
    b.free(d);

    // ── un-cond-shift2 (right shift gated by s2), re-inserting ovf2 at top ──
    for i in 0..n - 1 {
        if dialog_gcd_skip_zero_edge_apply_halve_cshift_enabled() && i == 0 {
            continue;
        }
        cswap(b, s2, y[i], y[i + 1]);
    }
    cswap(b, s2, y[n - 1], ovf2);
    // The boundary cswap already pulled the vacated top bit (0) into ovf2, so
    // ovf2 is |0> here. (A `ccx(s2, y[n-1], ovf2)` would WRONGLY re-set it to
    // s2&y[n-1] = e, dirtying the ancilla — the free's reset then masks the
    // value error but leaks global phase. So: no extra clear.)
    b.free(ovf2);

    // ── un-shift1 (unconditional right shift), re-inserting ovf1 at top ──
    for i in 0..n - 1 {
        b.swap(y[i], y[i + 1]);
    }
    b.swap(y[n - 1], ovf1);
    // The swap already pulled the vacated top bit (0) into ovf1, so ovf1 is |0>
    // here. (A `cx(y[n-1], ovf1)` would re-dirty it — see ovf2 note above.)
    b.free(ovf1);
}
