use super::{BitId, Op, OperationType, QubitId, RegisterId};
use std::fs::{remove_file, File};
use std::io::{BufReader, BufWriter, Read, Write};

type TmCircuit = super::trailmix_port::circuit::Circuit;
type TmOp = super::trailmix_port::circuit::Op;

pub fn run_product_column_ghost_canary() {
    fn bit(shot: usize, lane: usize, salt: u64) -> bool {
        let mut x = ((shot as u64) + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15);
        x ^= ((lane as u64) + 0x51).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        x ^= salt;
        x ^= x >> 30;
        x = x.wrapping_mul(0xbf58_476d_1ce4_e5b9);
        x ^= x >> 27;
        x = x.wrapping_mul(0x94d0_49bb_1331_11eb);
        ((x ^ (x >> 31)) & 1) == 1
    }

    let mut circ = TmCircuit::new();
    circ.set_max_qubit_peak(64);
    let carry = circ.alloc_input_qreg("carry_lsb");
    let xs: Vec<_> = (0..7)
        .map(|i| circ.alloc_input_qreg(&format!("x_{i}")))
        .collect();
    let ys: Vec<_> = (0..7)
        .map(|i| circ.alloc_input_qreg(&format!("y_{i}")))
        .collect();
    let product_bit = circ.alloc_qreg("product_column_bit");

    for shot in 0..64 {
        if bit(shot, 0, 0x1234) {
            circ.sim_load_reg_bytes_shot(std::slice::from_ref(&carry), &[1u8], shot);
        }
        for (lane, q) in xs.iter().chain(ys.iter()).enumerate() {
            if bit(shot, lane + 1, 0x9876) {
                circ.sim_load_reg_bytes_shot(std::slice::from_ref(q), &[1u8], shot);
            }
        }
    }

    let column_terms = [
        (0usize, 6usize),
        (1, 5),
        (2, 4),
        (3, 3),
        (4, 2),
        (5, 1),
        (6, 0),
    ];
    circ.cx(&carry, &product_bit);
    for (xi, yi) in column_terms {
        circ.ccx(&xs[xi], &ys[yi], &product_bit);
    }

    let mut ghost = circ.hmr_ghost(&product_bit);
    circ.zero_and_free(product_bit);
    circ.ghost_xor_z(&mut ghost, &carry);
    for (xi, yi) in column_terms {
        circ.ghost_xor_cz(&mut ghost, &xs[xi], &ys[yi]);
    }
    circ.close_ghost(ghost);
    circ.assert_phase_clean();

    let mut outs = vec![carry];
    outs.extend(xs);
    outs.extend(ys);
    let _ = circ.destroy_sim(outs);

    eprintln!(
        "trailmix_product_column_ghost_canary: PASS peak_qubits={} total_qubits={} total_ops={} ccx={}",
        circ.peak_qubits,
        circ.total_qubits(),
        circ.total_ops(),
        circ.ccx_count(),
    );
}

pub fn run_first_event_index_canary() {
    use super::trailmix_port::arith::mcx::mcx_clean_k;
    use super::trailmix_port::circuit::QReg;

    fn bit(shot: usize, lane: usize, salt: u64) -> bool {
        let mut x = ((shot as u64) + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15);
        x ^= ((lane as u64) + 0x73).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        x ^= salt;
        x ^= x >> 30;
        x = x.wrapping_mul(0xbf58_476d_1ce4_e5b9);
        x ^= x >> 27;
        x = x.wrapping_mul(0x94d0_49bb_1331_11eb);
        ((x ^ (x >> 31)) & 1) == 1
    }

    fn product_bit_for_shot(shot: usize, col: usize) -> bool {
        let mut v = bit(shot, 100 + col, 0xc001);
        for xi in 0..8 {
            if col >= xi {
                let yi = col - xi;
                if yi < 8 {
                    v ^= bit(shot, xi, 0xa11c) & bit(shot, 20 + yi, 0xb22d);
                }
            }
        }
        v
    }

    fn expected_first_event_masks(want_one: bool) -> (u64, [u64; 4]) {
        let mut found_mask = 0u64;
        let mut idx_masks = [0u64; 4];
        for shot in 0..64 {
            let mut hit_idx = None;
            for col in 0..16 {
                if product_bit_for_shot(shot, col) == want_one {
                    hit_idx = Some(col);
                    break;
                }
            }
            if let Some(col) = hit_idx {
                found_mask |= 1u64 << shot;
                for (bit_idx, mask) in idx_masks.iter_mut().enumerate() {
                    if ((col >> bit_idx) & 1) == 1 {
                        *mask |= 1u64 << shot;
                    }
                }
            }
        }
        (found_mask, idx_masks)
    }

    fn build_column(
        circ: &mut TmCircuit,
        carries: &[QReg],
        xs: &[QReg],
        ys: &[QReg],
        col: usize,
    ) -> (QReg, Vec<(usize, usize)>) {
        let product_bit = circ.alloc_qreg(&format!("product_col_{col}"));
        let mut terms = Vec::new();
        circ.cx(&carries[col], &product_bit);
        for xi in 0..xs.len() {
            if col >= xi {
                let yi = col - xi;
                if yi < ys.len() {
                    circ.ccx(&xs[xi], &ys[yi], &product_bit);
                    terms.push((xi, yi));
                }
            }
        }
        (product_bit, terms)
    }

    fn ghost_free_column(
        circ: &mut TmCircuit,
        product_bit: QReg,
        carry: &QReg,
        xs: &[QReg],
        ys: &[QReg],
        terms: &[(usize, usize)],
    ) {
        let mut ghost = circ.hmr_ghost(&product_bit);
        circ.zero_and_free(product_bit);
        circ.ghost_xor_z(&mut ghost, carry);
        for &(xi, yi) in terms {
            circ.ghost_xor_cz(&mut ghost, &xs[xi], &ys[yi]);
        }
        circ.close_ghost(ghost);
    }

    fn with_idx_const_controls<F: FnOnce(&mut TmCircuit, Vec<&QReg>)>(
        circ: &mut TmCircuit,
        idx: &[QReg],
        value: usize,
        event: &QReg,
        f: F,
    ) {
        for (bit_idx, q) in idx.iter().enumerate() {
            if ((value >> bit_idx) & 1) == 0 {
                circ.x(q);
            }
        }
        let mut ctrls = Vec::with_capacity(idx.len() + 1);
        ctrls.push(event);
        ctrls.extend(idx.iter());
        f(circ, ctrls);
        for (bit_idx, q) in idx.iter().enumerate().rev() {
            if ((value >> bit_idx) & 1) == 0 {
                circ.x(q);
            }
        }
    }

    fn toggle_found_if_event_and_index(
        circ: &mut TmCircuit,
        event: &QReg,
        idx: &[QReg],
        value: usize,
        found: &QReg,
    ) {
        with_idx_const_controls(circ, idx, value, event, |circ, ctrls| {
            mcx_clean_k(circ, &ctrls, found);
        });
    }

    fn toggle_index_if_event_and_not_found(
        circ: &mut TmCircuit,
        event: &QReg,
        found: &QReg,
        idx: &[QReg],
        value: usize,
    ) {
        for (bit_idx, q) in idx.iter().enumerate() {
            if ((value >> bit_idx) & 1) == 1 {
                circ.x(found);
                circ.ccx(event, found, q);
                circ.x(found);
            }
        }
    }

    fn first_event_forward(
        circ: &mut TmCircuit,
        event: &QReg,
        found: &QReg,
        idx: &[QReg],
        value: usize,
    ) {
        toggle_index_if_event_and_not_found(circ, event, found, idx, value);
        toggle_found_if_event_and_index(circ, event, idx, value, found);
    }

    fn first_event_reverse(
        circ: &mut TmCircuit,
        event: &QReg,
        found: &QReg,
        idx: &[QReg],
        value: usize,
    ) {
        toggle_found_if_event_and_index(circ, event, idx, value, found);
        toggle_index_if_event_and_not_found(circ, event, found, idx, value);
    }

    let mut circ = TmCircuit::new();
    circ.set_max_qubit_peak(96);
    let carries: Vec<_> = (0..16)
        .map(|i| circ.alloc_input_qreg(&format!("carry_{i}")))
        .collect();
    let xs: Vec<_> = (0..8)
        .map(|i| circ.alloc_input_qreg(&format!("x_{i}")))
        .collect();
    let ys: Vec<_> = (0..8)
        .map(|i| circ.alloc_input_qreg(&format!("y_{i}")))
        .collect();
    let first_one_found = circ.alloc_qreg("first_one_found");
    let first_one_idx = circ.alloc_qreg_bits("first_one_idx", 4);
    let first_zero_found = circ.alloc_qreg("first_zero_found");
    let first_zero_idx = circ.alloc_qreg_bits("first_zero_idx", 4);

    for shot in 0..64 {
        for (col, carry) in carries.iter().enumerate() {
            if bit(shot, 100 + col, 0xc001) {
                circ.sim_load_reg_bytes_shot(std::slice::from_ref(carry), &[1u8], shot);
            }
        }
        for (lane, q) in xs.iter().enumerate() {
            if bit(shot, lane, 0xa11c) {
                circ.sim_load_reg_bytes_shot(std::slice::from_ref(q), &[1u8], shot);
            }
        }
        for (lane, q) in ys.iter().enumerate() {
            if bit(shot, 20 + lane, 0xb22d) {
                circ.sim_load_reg_bytes_shot(std::slice::from_ref(q), &[1u8], shot);
            }
        }
    }

    for col in 0..16 {
        let (product_bit, terms) = build_column(&mut circ, &carries, &xs, &ys, col);
        first_event_forward(
            &mut circ,
            &product_bit,
            &first_one_found,
            &first_one_idx,
            col,
        );
        circ.x(&product_bit);
        first_event_forward(
            &mut circ,
            &product_bit,
            &first_zero_found,
            &first_zero_idx,
            col,
        );
        circ.x(&product_bit);
        ghost_free_column(&mut circ, product_bit, &carries[col], &xs, &ys, &terms);
    }

    let (one_found_expected, one_idx_expected) = expected_first_event_masks(true);
    assert_eq!(
        circ.sim_get_mask(first_one_found.id()),
        one_found_expected,
        "first-one found mask mismatch"
    );
    for (bit_idx, expected) in one_idx_expected.iter().enumerate() {
        assert_eq!(
            circ.sim_get_mask(first_one_idx[bit_idx].id()),
            *expected,
            "first-one index bit {bit_idx} mismatch"
        );
    }

    let (zero_found_expected, zero_idx_expected) = expected_first_event_masks(false);
    assert_eq!(
        circ.sim_get_mask(first_zero_found.id()),
        zero_found_expected,
        "first-zero found mask mismatch"
    );
    for (bit_idx, expected) in zero_idx_expected.iter().enumerate() {
        assert_eq!(
            circ.sim_get_mask(first_zero_idx[bit_idx].id()),
            *expected,
            "first-zero index bit {bit_idx} mismatch"
        );
    }

    for col in (0..16).rev() {
        let (product_bit, terms) = build_column(&mut circ, &carries, &xs, &ys, col);
        circ.x(&product_bit);
        first_event_reverse(
            &mut circ,
            &product_bit,
            &first_zero_found,
            &first_zero_idx,
            col,
        );
        circ.x(&product_bit);
        first_event_reverse(
            &mut circ,
            &product_bit,
            &first_one_found,
            &first_one_idx,
            col,
        );
        ghost_free_column(&mut circ, product_bit, &carries[col], &xs, &ys, &terms);
    }

    circ.zero_and_free(first_zero_found);
    for q in first_zero_idx {
        circ.zero_and_free(q);
    }
    circ.zero_and_free(first_one_found);
    for q in first_one_idx {
        circ.zero_and_free(q);
    }
    circ.assert_phase_clean();

    let mut outs = carries;
    outs.extend(xs);
    outs.extend(ys);
    let _ = circ.destroy_sim(outs);

    eprintln!(
        "trailmix_first_event_index_canary: PASS peak_qubits={} total_qubits={} total_ops={} ccx={} ccz={}",
        circ.peak_qubits,
        circ.total_qubits(),
        circ.total_ops(),
        circ.ccx_count(),
        circ.ccz_emitted,
    );
}

pub fn run_secp_shape_stream_compare_canary() {
    use super::trailmix_port::arith::mcx::mcx_clean_k;
    use super::trailmix_port::circuit::QReg;

    #[derive(Clone, Copy)]
    enum Event {
        Bit { col: usize, invert: bool },
        And { a: usize, b: usize },
    }

    fn bits_needed(len: usize) -> usize {
        let mut bits = 0usize;
        let mut cap = 1usize;
        while cap < len {
            bits += 1;
            cap <<= 1;
        }
        bits.max(1)
    }

    fn p_bit(col: usize) -> bool {
        match col {
            0..=9 => ((47usize >> col) & 1) == 1,
            10..=31 => true,
            32 => false,
            33..=255 => true,
            _ => false,
        }
    }

    fn z_bit_for_shot(shot: usize, col: usize) -> bool {
        match shot % 8 {
            0 => {
                if col <= 9 {
                    ((46usize >> col) & 1) == 1
                } else {
                    p_bit(col)
                }
            }
            1 => p_bit(col),
            2 => {
                if col <= 9 {
                    ((48usize >> col) & 1) == 1
                } else {
                    p_bit(col)
                }
            }
            3 => col == 256 || p_bit(col),
            4 => col != 200 && p_bit(col),
            5 => match col {
                32 => true,
                33..=255 => true,
                _ => false,
            },
            6 => match col {
                0..=9 => true,
                10..=31 => col != 20,
                32 => false,
                33..=255 => true,
                _ => false,
            },
            _ => {
                if col <= 9 {
                    ((48usize >> col) & 1) == 1
                } else {
                    matches!(col, 10..=31 | 33..=255)
                }
            }
        }
    }

    fn event_value(shot: usize, event: Event) -> bool {
        match event {
            Event::Bit { col, invert } => z_bit_for_shot(shot, col) ^ invert,
            Event::And { a, b } => z_bit_for_shot(shot, a) && z_bit_for_shot(shot, b),
        }
    }

    fn expected_event_masks(events: &[Event]) -> (u64, Vec<u64>) {
        let idx_bits = bits_needed(events.len());
        let mut found_mask = 0u64;
        let mut idx_masks = vec![0u64; idx_bits];
        for shot in 0..64 {
            let mut hit_idx = None;
            for (pos, event) in events.iter().copied().enumerate() {
                if event_value(shot, event) {
                    hit_idx = Some(pos);
                    break;
                }
            }
            if let Some(pos) = hit_idx {
                found_mask |= 1u64 << shot;
                for (bit_idx, mask) in idx_masks.iter_mut().enumerate() {
                    if ((pos >> bit_idx) & 1) == 1 {
                        *mask |= 1u64 << shot;
                    }
                }
            }
        }
        (found_mask, idx_masks)
    }

    fn with_idx_const_controls<F: FnOnce(&mut TmCircuit, Vec<&QReg>)>(
        circ: &mut TmCircuit,
        idx: &[QReg],
        value: usize,
        event: &QReg,
        f: F,
    ) {
        for (bit_idx, q) in idx.iter().enumerate() {
            if ((value >> bit_idx) & 1) == 0 {
                circ.x(q);
            }
        }
        let mut ctrls = Vec::with_capacity(idx.len() + 1);
        ctrls.push(event);
        ctrls.extend(idx.iter());
        f(circ, ctrls);
        for (bit_idx, q) in idx.iter().enumerate().rev() {
            if ((value >> bit_idx) & 1) == 0 {
                circ.x(q);
            }
        }
    }

    fn toggle_found_if_event_and_index(
        circ: &mut TmCircuit,
        event: &QReg,
        idx: &[QReg],
        value: usize,
        found: &QReg,
    ) {
        with_idx_const_controls(circ, idx, value, event, |circ, ctrls| {
            mcx_clean_k(circ, &ctrls, found);
        });
    }

    fn toggle_index_if_event_and_not_found(
        circ: &mut TmCircuit,
        event: &QReg,
        found: &QReg,
        idx: &[QReg],
        value: usize,
    ) {
        for (bit_idx, q) in idx.iter().enumerate() {
            if ((value >> bit_idx) & 1) == 1 {
                circ.x(found);
                circ.ccx(event, found, q);
                circ.x(found);
            }
        }
    }

    fn first_event_forward(
        circ: &mut TmCircuit,
        event: &QReg,
        found: &QReg,
        idx: &[QReg],
        value: usize,
    ) {
        toggle_index_if_event_and_not_found(circ, event, found, idx, value);
        toggle_found_if_event_and_index(circ, event, idx, value, found);
    }

    fn first_event_reverse(
        circ: &mut TmCircuit,
        event: &QReg,
        found: &QReg,
        idx: &[QReg],
        value: usize,
    ) {
        toggle_found_if_event_and_index(circ, event, idx, value, found);
        toggle_index_if_event_and_not_found(circ, event, found, idx, value);
    }

    fn ghost_free_and(circ: &mut TmCircuit, tmp: QReg, a: &QReg, b: &QReg) {
        let mut ghost = circ.hmr_ghost(&tmp);
        circ.zero_and_free(tmp);
        circ.ghost_xor_cz(&mut ghost, a, b);
        circ.close_ghost(ghost);
    }

    fn emit_event_forward(
        circ: &mut TmCircuit,
        z: &[QReg],
        event: Event,
        found: &QReg,
        idx: &[QReg],
        pos: usize,
    ) {
        match event {
            Event::Bit { col, invert } => {
                if invert {
                    circ.x(&z[col]);
                }
                first_event_forward(circ, &z[col], found, idx, pos);
                if invert {
                    circ.x(&z[col]);
                }
            }
            Event::And { a, b } => {
                let tmp = circ.alloc_qreg("event_and");
                circ.ccx(&z[a], &z[b], &tmp);
                first_event_forward(circ, &tmp, found, idx, pos);
                ghost_free_and(circ, tmp, &z[a], &z[b]);
            }
        }
    }

    fn emit_event_reverse(
        circ: &mut TmCircuit,
        z: &[QReg],
        event: Event,
        found: &QReg,
        idx: &[QReg],
        pos: usize,
    ) {
        match event {
            Event::Bit { col, invert } => {
                if invert {
                    circ.x(&z[col]);
                }
                first_event_reverse(circ, &z[col], found, idx, pos);
                if invert {
                    circ.x(&z[col]);
                }
            }
            Event::And { a, b } => {
                let tmp = circ.alloc_qreg("event_and");
                circ.ccx(&z[a], &z[b], &tmp);
                first_event_reverse(circ, &tmp, found, idx, pos);
                ghost_free_and(circ, tmp, &z[a], &z[b]);
            }
        }
    }

    fn assert_witness(circ: &TmCircuit, label: &str, found: &QReg, idx: &[QReg], events: &[Event]) {
        let (expected_found, expected_idx) = expected_event_masks(events);
        assert_eq!(
            circ.sim_get_mask(found.id()),
            expected_found,
            "{label} found mask mismatch"
        );
        for (bit_idx, expected) in expected_idx.iter().enumerate() {
            assert_eq!(
                circ.sim_get_mask(idx[bit_idx].id()),
                *expected,
                "{label} index bit {bit_idx} mismatch"
            );
        }
    }

    let high_nonzero_events: Vec<_> = (256..=276)
        .map(|col| Event::Bit { col, invert: false })
        .collect();
    let high_zero_events: Vec<_> = (33..=255)
        .rev()
        .map(|col| Event::Bit { col, invert: true })
        .collect();
    let mid_zero_events: Vec<_> = (10..=31)
        .rev()
        .map(|col| Event::Bit { col, invert: true })
        .collect();
    let low10_ge48_events = vec![
        Event::Bit {
            col: 6,
            invert: false,
        },
        Event::Bit {
            col: 7,
            invert: false,
        },
        Event::Bit {
            col: 8,
            invert: false,
        },
        Event::Bit {
            col: 9,
            invert: false,
        },
        Event::And { a: 5, b: 4 },
    ];

    let mut circ = TmCircuit::new();
    circ.set_max_qubit_peak(512);
    let z: Vec<_> = (0..=276)
        .map(|i| circ.alloc_input_qreg(&format!("z_{i}")))
        .collect();

    for shot in 0..64 {
        for (col, q) in z.iter().enumerate() {
            if z_bit_for_shot(shot, col) {
                circ.sim_load_reg_bytes_shot(std::slice::from_ref(q), &[1u8], shot);
            }
        }
    }

    let high_nonzero_found = circ.alloc_qreg("high_nonzero_found");
    let high_nonzero_idx =
        circ.alloc_qreg_bits("high_nonzero_idx", bits_needed(high_nonzero_events.len()));
    let high_zero_found = circ.alloc_qreg("high_zero_found");
    let high_zero_idx = circ.alloc_qreg_bits("high_zero_idx", bits_needed(high_zero_events.len()));
    let mid_zero_found = circ.alloc_qreg("mid_zero_found");
    let mid_zero_idx = circ.alloc_qreg_bits("mid_zero_idx", bits_needed(mid_zero_events.len()));
    let low10_ge48_found = circ.alloc_qreg("low10_ge48_found");
    let low10_ge48_idx =
        circ.alloc_qreg_bits("low10_ge48_idx", bits_needed(low10_ge48_events.len()));

    for (pos, event) in high_nonzero_events.iter().copied().enumerate() {
        emit_event_forward(
            &mut circ,
            &z,
            event,
            &high_nonzero_found,
            &high_nonzero_idx,
            pos,
        );
    }
    for (pos, event) in high_zero_events.iter().copied().enumerate() {
        emit_event_forward(&mut circ, &z, event, &high_zero_found, &high_zero_idx, pos);
    }
    for (pos, event) in mid_zero_events.iter().copied().enumerate() {
        emit_event_forward(&mut circ, &z, event, &mid_zero_found, &mid_zero_idx, pos);
    }
    for (pos, event) in low10_ge48_events.iter().copied().enumerate() {
        emit_event_forward(
            &mut circ,
            &z,
            event,
            &low10_ge48_found,
            &low10_ge48_idx,
            pos,
        );
    }

    assert_witness(
        &circ,
        "high_nonzero",
        &high_nonzero_found,
        &high_nonzero_idx,
        &high_nonzero_events,
    );
    assert_witness(
        &circ,
        "high_zero",
        &high_zero_found,
        &high_zero_idx,
        &high_zero_events,
    );
    assert_witness(
        &circ,
        "mid_zero",
        &mid_zero_found,
        &mid_zero_idx,
        &mid_zero_events,
    );
    assert_witness(
        &circ,
        "low10_ge48",
        &low10_ge48_found,
        &low10_ge48_idx,
        &low10_ge48_events,
    );

    for (pos, event) in low10_ge48_events.iter().copied().enumerate().rev() {
        emit_event_reverse(
            &mut circ,
            &z,
            event,
            &low10_ge48_found,
            &low10_ge48_idx,
            pos,
        );
    }
    for (pos, event) in mid_zero_events.iter().copied().enumerate().rev() {
        emit_event_reverse(&mut circ, &z, event, &mid_zero_found, &mid_zero_idx, pos);
    }
    for (pos, event) in high_zero_events.iter().copied().enumerate().rev() {
        emit_event_reverse(&mut circ, &z, event, &high_zero_found, &high_zero_idx, pos);
    }
    for (pos, event) in high_nonzero_events.iter().copied().enumerate().rev() {
        emit_event_reverse(
            &mut circ,
            &z,
            event,
            &high_nonzero_found,
            &high_nonzero_idx,
            pos,
        );
    }

    circ.zero_and_free(low10_ge48_found);
    for q in low10_ge48_idx {
        circ.zero_and_free(q);
    }
    circ.zero_and_free(mid_zero_found);
    for q in mid_zero_idx {
        circ.zero_and_free(q);
    }
    circ.zero_and_free(high_zero_found);
    for q in high_zero_idx {
        circ.zero_and_free(q);
    }
    circ.zero_and_free(high_nonzero_found);
    for q in high_nonzero_idx {
        circ.zero_and_free(q);
    }
    circ.assert_phase_clean();

    let _ = circ.destroy_sim(z);

    eprintln!(
        "trailmix_secp_shape_stream_compare_canary: PASS peak_qubits={} total_qubits={} total_ops={} ccx={} ccz={}",
        circ.peak_qubits,
        circ.total_qubits(),
        circ.total_ops(),
        circ.ccx_count(),
        circ.ccz_emitted,
    );
}

pub fn run_secp_shape_product_source_canary() {
    use super::trailmix_port::arith::mcx::mcx_clean_k;
    use super::trailmix_port::circuit::QReg;

    #[derive(Clone, Copy)]
    enum Event {
        Bit { col: usize, invert: bool },
        And { a: usize, b: usize },
    }

    fn bits_needed(len: usize) -> usize {
        let mut bits = 0usize;
        let mut cap = 1usize;
        while cap < len {
            bits += 1;
            cap <<= 1;
        }
        bits.max(1)
    }

    fn bit(shot: usize, lane: usize, salt: u64) -> bool {
        let mut x = ((shot as u64) + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15);
        x ^= ((lane as u64) + 0x97).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        x ^= salt;
        x ^= x >> 30;
        x = x.wrapping_mul(0xbf58_476d_1ce4_e5b9);
        x ^= x >> 27;
        x = x.wrapping_mul(0x94d0_49bb_1331_11eb);
        ((x ^ (x >> 31)) & 1) == 1
    }

    fn x_bit_for_shot(shot: usize, xi: usize) -> bool {
        bit(shot, xi, 0xa11c)
    }

    fn y_bit_for_shot(shot: usize, yi: usize) -> bool {
        bit(shot, 100 + yi, 0xb22d)
    }

    fn product_terms_value(shot: usize, col: usize, x_bits: usize, y_bits: usize) -> bool {
        let mut v = false;
        for xi in 0..x_bits {
            if col >= xi {
                let yi = col - xi;
                if yi < y_bits {
                    v ^= x_bit_for_shot(shot, xi) & y_bit_for_shot(shot, yi);
                }
            }
        }
        v
    }

    fn p_bit(col: usize) -> bool {
        match col {
            0..=9 => ((47usize >> col) & 1) == 1,
            10..=31 => true,
            32 => false,
            33..=255 => true,
            _ => false,
        }
    }

    fn z_bit_for_shot(shot: usize, col: usize) -> bool {
        match shot % 8 {
            0 => {
                if col <= 9 {
                    ((46usize >> col) & 1) == 1
                } else {
                    p_bit(col)
                }
            }
            1 => p_bit(col),
            2 => {
                if col <= 9 {
                    ((48usize >> col) & 1) == 1
                } else {
                    p_bit(col)
                }
            }
            3 => col == 256 || p_bit(col),
            4 => col != 200 && p_bit(col),
            5 => match col {
                32 => true,
                33..=255 => true,
                _ => false,
            },
            6 => match col {
                0..=9 => true,
                10..=31 => col != 20,
                32 => false,
                33..=255 => true,
                _ => false,
            },
            _ => {
                if col <= 9 {
                    ((48usize >> col) & 1) == 1
                } else {
                    matches!(col, 10..=31 | 33..=255)
                }
            }
        }
    }

    fn carry_bit_for_shot(shot: usize, col: usize, x_bits: usize, y_bits: usize) -> bool {
        z_bit_for_shot(shot, col) ^ product_terms_value(shot, col, x_bits, y_bits)
    }

    fn event_value(shot: usize, event: Event) -> bool {
        match event {
            Event::Bit { col, invert } => z_bit_for_shot(shot, col) ^ invert,
            Event::And { a, b } => z_bit_for_shot(shot, a) && z_bit_for_shot(shot, b),
        }
    }

    fn expected_event_masks(events: &[Event]) -> (u64, Vec<u64>) {
        let idx_bits = bits_needed(events.len());
        let mut found_mask = 0u64;
        let mut idx_masks = vec![0u64; idx_bits];
        for shot in 0..64 {
            let mut hit_idx = None;
            for (pos, event) in events.iter().copied().enumerate() {
                if event_value(shot, event) {
                    hit_idx = Some(pos);
                    break;
                }
            }
            if let Some(pos) = hit_idx {
                found_mask |= 1u64 << shot;
                for (bit_idx, mask) in idx_masks.iter_mut().enumerate() {
                    if ((pos >> bit_idx) & 1) == 1 {
                        *mask |= 1u64 << shot;
                    }
                }
            }
        }
        (found_mask, idx_masks)
    }

    fn with_idx_const_controls<F: FnOnce(&mut TmCircuit, Vec<&QReg>)>(
        circ: &mut TmCircuit,
        idx: &[QReg],
        value: usize,
        event: &QReg,
        f: F,
    ) {
        for (bit_idx, q) in idx.iter().enumerate() {
            if ((value >> bit_idx) & 1) == 0 {
                circ.x(q);
            }
        }
        let mut ctrls = Vec::with_capacity(idx.len() + 1);
        ctrls.push(event);
        ctrls.extend(idx.iter());
        f(circ, ctrls);
        for (bit_idx, q) in idx.iter().enumerate().rev() {
            if ((value >> bit_idx) & 1) == 0 {
                circ.x(q);
            }
        }
    }

    fn toggle_found_if_event_and_index(
        circ: &mut TmCircuit,
        event: &QReg,
        idx: &[QReg],
        value: usize,
        found: &QReg,
    ) {
        with_idx_const_controls(circ, idx, value, event, |circ, ctrls| {
            mcx_clean_k(circ, &ctrls, found);
        });
    }

    fn toggle_index_if_event_and_not_found(
        circ: &mut TmCircuit,
        event: &QReg,
        found: &QReg,
        idx: &[QReg],
        value: usize,
    ) {
        for (bit_idx, q) in idx.iter().enumerate() {
            if ((value >> bit_idx) & 1) == 1 {
                circ.x(found);
                circ.ccx(event, found, q);
                circ.x(found);
            }
        }
    }

    fn first_event_forward(
        circ: &mut TmCircuit,
        event: &QReg,
        found: &QReg,
        idx: &[QReg],
        value: usize,
    ) {
        toggle_index_if_event_and_not_found(circ, event, found, idx, value);
        toggle_found_if_event_and_index(circ, event, idx, value, found);
    }

    fn first_event_reverse(
        circ: &mut TmCircuit,
        event: &QReg,
        found: &QReg,
        idx: &[QReg],
        value: usize,
    ) {
        toggle_found_if_event_and_index(circ, event, idx, value, found);
        toggle_index_if_event_and_not_found(circ, event, found, idx, value);
    }

    fn build_product_column(
        circ: &mut TmCircuit,
        carries: &[QReg],
        xs: &[QReg],
        ys: &[QReg],
        col: usize,
    ) -> QReg {
        let product = circ.alloc_qreg("stream_product_bit");
        circ.cx(&carries[col], &product);
        for xi in 0..xs.len() {
            if col >= xi {
                let yi = col - xi;
                if yi < ys.len() {
                    circ.ccx(&xs[xi], &ys[yi], &product);
                }
            }
        }
        product
    }

    fn ghost_free_column(
        circ: &mut TmCircuit,
        product: QReg,
        carries: &[QReg],
        xs: &[QReg],
        ys: &[QReg],
        col: usize,
    ) {
        let mut ghost = circ.hmr_ghost(&product);
        circ.zero_and_free(product);
        circ.ghost_xor_z(&mut ghost, &carries[col]);
        for xi in 0..xs.len() {
            if col >= xi {
                let yi = col - xi;
                if yi < ys.len() {
                    circ.ghost_xor_cz(&mut ghost, &xs[xi], &ys[yi]);
                }
            }
        }
        circ.close_ghost(ghost);
    }

    fn ghost_free_and(circ: &mut TmCircuit, tmp: QReg, a: &QReg, b: &QReg) {
        let mut ghost = circ.hmr_ghost(&tmp);
        circ.zero_and_free(tmp);
        circ.ghost_xor_cz(&mut ghost, a, b);
        circ.close_ghost(ghost);
    }

    fn emit_event_forward(
        circ: &mut TmCircuit,
        carries: &[QReg],
        xs: &[QReg],
        ys: &[QReg],
        event: Event,
        found: &QReg,
        idx: &[QReg],
        pos: usize,
    ) {
        match event {
            Event::Bit { col, invert } => {
                let product = build_product_column(circ, carries, xs, ys, col);
                if invert {
                    circ.x(&product);
                }
                first_event_forward(circ, &product, found, idx, pos);
                if invert {
                    circ.x(&product);
                }
                ghost_free_column(circ, product, carries, xs, ys, col);
            }
            Event::And { a, b } => {
                let za = build_product_column(circ, carries, xs, ys, a);
                let zb = build_product_column(circ, carries, xs, ys, b);
                let tmp = circ.alloc_qreg("event_and");
                circ.ccx(&za, &zb, &tmp);
                first_event_forward(circ, &tmp, found, idx, pos);
                ghost_free_and(circ, tmp, &za, &zb);
                ghost_free_column(circ, zb, carries, xs, ys, b);
                ghost_free_column(circ, za, carries, xs, ys, a);
            }
        }
    }

    fn emit_event_reverse(
        circ: &mut TmCircuit,
        carries: &[QReg],
        xs: &[QReg],
        ys: &[QReg],
        event: Event,
        found: &QReg,
        idx: &[QReg],
        pos: usize,
    ) {
        match event {
            Event::Bit { col, invert } => {
                let product = build_product_column(circ, carries, xs, ys, col);
                if invert {
                    circ.x(&product);
                }
                first_event_reverse(circ, &product, found, idx, pos);
                if invert {
                    circ.x(&product);
                }
                ghost_free_column(circ, product, carries, xs, ys, col);
            }
            Event::And { a, b } => {
                let za = build_product_column(circ, carries, xs, ys, a);
                let zb = build_product_column(circ, carries, xs, ys, b);
                let tmp = circ.alloc_qreg("event_and");
                circ.ccx(&za, &zb, &tmp);
                first_event_reverse(circ, &tmp, found, idx, pos);
                ghost_free_and(circ, tmp, &za, &zb);
                ghost_free_column(circ, zb, carries, xs, ys, b);
                ghost_free_column(circ, za, carries, xs, ys, a);
            }
        }
    }

    fn assert_witness(circ: &TmCircuit, label: &str, found: &QReg, idx: &[QReg], events: &[Event]) {
        let (expected_found, expected_idx) = expected_event_masks(events);
        assert_eq!(
            circ.sim_get_mask(found.id()),
            expected_found,
            "{label} found mask mismatch"
        );
        for (bit_idx, expected) in expected_idx.iter().enumerate() {
            assert_eq!(
                circ.sim_get_mask(idx[bit_idx].id()),
                *expected,
                "{label} index bit {bit_idx} mismatch"
            );
        }
    }

    let high_nonzero_events: Vec<_> = (256..=276)
        .map(|col| Event::Bit { col, invert: false })
        .collect();
    let high_zero_events: Vec<_> = (33..=255)
        .rev()
        .map(|col| Event::Bit { col, invert: true })
        .collect();
    let mid_zero_events: Vec<_> = (10..=31)
        .rev()
        .map(|col| Event::Bit { col, invert: true })
        .collect();
    let low10_ge48_events = vec![
        Event::Bit {
            col: 6,
            invert: false,
        },
        Event::Bit {
            col: 7,
            invert: false,
        },
        Event::Bit {
            col: 8,
            invert: false,
        },
        Event::Bit {
            col: 9,
            invert: false,
        },
        Event::And { a: 5, b: 4 },
    ];

    let x_bits = 32usize;
    let y_bits = 32usize;

    let mut circ = TmCircuit::new();
    circ.set_max_qubit_peak(448);
    let carries: Vec<_> = (0..=276)
        .map(|i| circ.alloc_input_qreg(&format!("carry_lsb_{i}")))
        .collect();
    let xs: Vec<_> = (0..x_bits)
        .map(|i| circ.alloc_input_qreg(&format!("x_{i}")))
        .collect();
    let ys: Vec<_> = (0..y_bits)
        .map(|i| circ.alloc_input_qreg(&format!("y_{i}")))
        .collect();

    for shot in 0..64 {
        for (xi, q) in xs.iter().enumerate() {
            if x_bit_for_shot(shot, xi) {
                circ.sim_load_reg_bytes_shot(std::slice::from_ref(q), &[1u8], shot);
            }
        }
        for (yi, q) in ys.iter().enumerate() {
            if y_bit_for_shot(shot, yi) {
                circ.sim_load_reg_bytes_shot(std::slice::from_ref(q), &[1u8], shot);
            }
        }
        for (col, q) in carries.iter().enumerate() {
            if carry_bit_for_shot(shot, col, x_bits, y_bits) {
                circ.sim_load_reg_bytes_shot(std::slice::from_ref(q), &[1u8], shot);
            }
        }
    }

    let high_nonzero_found = circ.alloc_qreg("high_nonzero_found");
    let high_nonzero_idx =
        circ.alloc_qreg_bits("high_nonzero_idx", bits_needed(high_nonzero_events.len()));
    let high_zero_found = circ.alloc_qreg("high_zero_found");
    let high_zero_idx = circ.alloc_qreg_bits("high_zero_idx", bits_needed(high_zero_events.len()));
    let mid_zero_found = circ.alloc_qreg("mid_zero_found");
    let mid_zero_idx = circ.alloc_qreg_bits("mid_zero_idx", bits_needed(mid_zero_events.len()));
    let low10_ge48_found = circ.alloc_qreg("low10_ge48_found");
    let low10_ge48_idx =
        circ.alloc_qreg_bits("low10_ge48_idx", bits_needed(low10_ge48_events.len()));

    for (pos, event) in high_nonzero_events.iter().copied().enumerate() {
        emit_event_forward(
            &mut circ,
            &carries,
            &xs,
            &ys,
            event,
            &high_nonzero_found,
            &high_nonzero_idx,
            pos,
        );
    }
    for (pos, event) in high_zero_events.iter().copied().enumerate() {
        emit_event_forward(
            &mut circ,
            &carries,
            &xs,
            &ys,
            event,
            &high_zero_found,
            &high_zero_idx,
            pos,
        );
    }
    for (pos, event) in mid_zero_events.iter().copied().enumerate() {
        emit_event_forward(
            &mut circ,
            &carries,
            &xs,
            &ys,
            event,
            &mid_zero_found,
            &mid_zero_idx,
            pos,
        );
    }
    for (pos, event) in low10_ge48_events.iter().copied().enumerate() {
        emit_event_forward(
            &mut circ,
            &carries,
            &xs,
            &ys,
            event,
            &low10_ge48_found,
            &low10_ge48_idx,
            pos,
        );
    }

    assert_witness(
        &circ,
        "high_nonzero",
        &high_nonzero_found,
        &high_nonzero_idx,
        &high_nonzero_events,
    );
    assert_witness(
        &circ,
        "high_zero",
        &high_zero_found,
        &high_zero_idx,
        &high_zero_events,
    );
    assert_witness(
        &circ,
        "mid_zero",
        &mid_zero_found,
        &mid_zero_idx,
        &mid_zero_events,
    );
    assert_witness(
        &circ,
        "low10_ge48",
        &low10_ge48_found,
        &low10_ge48_idx,
        &low10_ge48_events,
    );

    for (pos, event) in low10_ge48_events.iter().copied().enumerate().rev() {
        emit_event_reverse(
            &mut circ,
            &carries,
            &xs,
            &ys,
            event,
            &low10_ge48_found,
            &low10_ge48_idx,
            pos,
        );
    }
    for (pos, event) in mid_zero_events.iter().copied().enumerate().rev() {
        emit_event_reverse(
            &mut circ,
            &carries,
            &xs,
            &ys,
            event,
            &mid_zero_found,
            &mid_zero_idx,
            pos,
        );
    }
    for (pos, event) in high_zero_events.iter().copied().enumerate().rev() {
        emit_event_reverse(
            &mut circ,
            &carries,
            &xs,
            &ys,
            event,
            &high_zero_found,
            &high_zero_idx,
            pos,
        );
    }
    for (pos, event) in high_nonzero_events.iter().copied().enumerate().rev() {
        emit_event_reverse(
            &mut circ,
            &carries,
            &xs,
            &ys,
            event,
            &high_nonzero_found,
            &high_nonzero_idx,
            pos,
        );
    }

    circ.zero_and_free(low10_ge48_found);
    for q in low10_ge48_idx {
        circ.zero_and_free(q);
    }
    circ.zero_and_free(mid_zero_found);
    for q in mid_zero_idx {
        circ.zero_and_free(q);
    }
    circ.zero_and_free(high_zero_found);
    for q in high_zero_idx {
        circ.zero_and_free(q);
    }
    circ.zero_and_free(high_nonzero_found);
    for q in high_nonzero_idx {
        circ.zero_and_free(q);
    }
    circ.assert_phase_clean();

    let mut outs = carries;
    outs.extend(xs);
    outs.extend(ys);
    let _ = circ.destroy_sim(outs);

    eprintln!(
        "trailmix_secp_shape_product_source_canary: PASS peak_qubits={} total_qubits={} total_ops={} ccx={} ccz={} x_bits={} y_bits={} carry_bits={}",
        circ.peak_qubits,
        circ.total_qubits(),
        circ.total_ops(),
        circ.ccx_count(),
        circ.ccz_emitted,
        x_bits,
        y_bits,
        277,
    );
}

pub fn run_cursor_pack_boundary_canary() {
    use super::trailmix_port::arith::cuccaro::{
        add_cuccaro, controlled_add_cuccaro_3n_refs, controlled_add_cuccaro_3n_reverse_refs,
    };
    use super::trailmix_port::arith::khattar_gidney::unary_iterate_log_star;
    use super::trailmix_port::circuit::{ContractSimView, QReg};

    fn load_bits(circ: &mut TmCircuit, regs: &[QReg], mut value: u64, shot: usize) {
        for q in regs {
            if (value & 1) == 1 {
                circ.sim_load_reg_bytes_shot(std::slice::from_ref(q), &[1u8], shot);
            }
            value >>= 1;
        }
    }

    fn read_bits(view: &ContractSimView<'_>, regs: &[&QReg], shot: usize) -> u64 {
        let mut value = 0u64;
        for (i, q) in regs.iter().enumerate() {
            if view.contract_read_bit_shot(q, shot) {
                value |= 1u64 << i;
            }
        }
        value
    }

    fn controlled_prefix_add_by_boundary(
        circ: &mut TmCircuit,
        boundary: &[QReg],
        a: &[QReg],
        b: &[QReg],
        reverse: bool,
    ) {
        let boundary_refs: Vec<&QReg> = boundary.iter().collect();
        unary_iterate_log_star(circ, &boundary_refs, a.len() + 1, |circ, prefix_len, gate| {
            if prefix_len == 0 {
                return;
            }
            let a_refs: Vec<&QReg> = a[..prefix_len].iter().collect();
            let b_refs: Vec<&QReg> = b[..prefix_len].iter().collect();
            if reverse {
                controlled_add_cuccaro_3n_reverse_refs(circ, gate, &a_refs, &b_refs);
            } else {
                controlled_add_cuccaro_3n_refs(circ, gate, &a_refs, &b_refs);
            }
        });
    }

    let mut fixed = TmCircuit::new();
    fixed.set_max_qubit_peak(96);
    let n = 16usize;
    let bnd = 9usize;
    let a: Vec<QReg> = (0..n)
        .map(|i| fixed.alloc_input_qreg(&format!("fixed_a{i}")))
        .collect();
    let b: Vec<QReg> = (0..n)
        .map(|i| fixed.alloc_input_qreg(&format!("fixed_b{i}")))
        .collect();
    for shot in 0..64 {
        let av = (shot as u64).wrapping_mul(0x9e37) ^ ((shot as u64) << 7);
        let bv = (shot as u64).wrapping_mul(0xbf59) ^ ((shot as u64) << 3);
        load_bits(&mut fixed, &a, av, shot);
        load_bits(&mut fixed, &b, bv, shot);
    }
    let a_hi_cap: Vec<&QReg> = a[bnd..].iter().collect();
    let b_all_cap: Vec<&QReg> = b.iter().collect();
    fixed.contract_capture(
        "cursor_pack.fixed.pre",
        move |view, shot| -> Result<(u64, u64), String> {
            Ok((
                read_bits(&view, &a_hi_cap, shot),
                read_bits(&view, &b_all_cap, shot),
            ))
        },
    );
    add_cuccaro(&mut fixed, &a[..bnd], &b[..bnd]);
    for q in &a[..bnd] {
        fixed.x(q);
    }
    add_cuccaro(&mut fixed, &a[..bnd], &b[..bnd]);
    for q in &a[..bnd] {
        fixed.x(q);
    }
    let a_hi_check: Vec<&QReg> = a[bnd..].iter().collect();
    let b_all_check: Vec<&QReg> = b.iter().collect();
    fixed.contract_pop_and_check::<(u64, u64), _>(
        "cursor_pack.fixed.pre",
        move |cap, view, shot| -> Result<(), String> {
            let (a_hi_pre, b_pre) = cap;
            let a_hi_post = read_bits(&view, &a_hi_check, shot);
            let b_post = read_bits(&view, &b_all_check, shot);
            if a_hi_post != *a_hi_pre {
                return Err(format!(
                    "fixed-boundary carry escaped on shot {shot}: {a_hi_pre:#x}->{a_hi_post:#x}"
                ));
            }
            if b_post != *b_pre {
                return Err(format!(
                    "fixed-boundary addend changed on shot {shot}: {b_pre:#x}->{b_post:#x}"
                ));
            }
            Ok(())
        },
    );
    fixed.assert_phase_clean();
    let mut fixed_outs = a;
    fixed_outs.extend(b);
    let _ = fixed.destroy_sim(fixed_outs);
    let fixed_peak = fixed.peak_qubits;
    let fixed_ops = fixed.total_ops();
    let fixed_ccx = fixed.ccx_count();

    let mut quantum = TmCircuit::new();
    quantum.set_max_qubit_peak(128);
    let n = 8usize;
    let a: Vec<QReg> = (0..n)
        .map(|i| quantum.alloc_input_qreg(&format!("q_a{i}")))
        .collect();
    let b: Vec<QReg> = (0..n)
        .map(|i| quantum.alloc_input_qreg(&format!("q_b{i}")))
        .collect();
    let boundary: Vec<QReg> = (0..4)
        .map(|i| quantum.alloc_input_qreg(&format!("boundary{i}")))
        .collect();
    for shot in 0..64 {
        let av = ((shot as u64).wrapping_mul(0x4d) ^ 0xa5) & ((1u64 << n) - 1);
        let bv = ((shot as u64).wrapping_mul(0x73) ^ 0x5a) & ((1u64 << n) - 1);
        let bnd = (shot % (n + 1)) as u8;
        load_bits(&mut quantum, &a, av, shot);
        load_bits(&mut quantum, &b, bv, shot);
        quantum.sim_load_reg_bytes_shot(&boundary, &[bnd], shot);
    }
    let a_roundtrip_pre: Vec<&QReg> = a.iter().collect();
    let b_roundtrip_pre: Vec<&QReg> = b.iter().collect();
    let boundary_roundtrip_pre: Vec<&QReg> = boundary.iter().collect();
    quantum.contract_capture(
        "cursor_pack.quantum.roundtrip_pre",
        move |view, shot| -> Result<(u64, u64, u64), String> {
            Ok((
                read_bits(&view, &a_roundtrip_pre, shot),
                read_bits(&view, &b_roundtrip_pre, shot),
                read_bits(&view, &boundary_roundtrip_pre, shot),
            ))
        },
    );
    let a_cap: Vec<&QReg> = a.iter().collect();
    let b_cap: Vec<&QReg> = b.iter().collect();
    let boundary_cap: Vec<&QReg> = boundary.iter().collect();
    quantum.contract_capture(
        "cursor_pack.quantum.pre",
        move |view, shot| -> Result<(u64, u64, u64), String> {
            Ok((
                read_bits(&view, &a_cap, shot),
                read_bits(&view, &b_cap, shot),
                read_bits(&view, &boundary_cap, shot),
            ))
        },
    );
    controlled_prefix_add_by_boundary(&mut quantum, &boundary, &a, &b, false);
    let a_after_add: Vec<&QReg> = a.iter().collect();
    let b_after_add: Vec<&QReg> = b.iter().collect();
    let boundary_after_add: Vec<&QReg> = boundary.iter().collect();
    quantum.contract_pop_and_check::<(u64, u64, u64), _>(
        "cursor_pack.quantum.pre",
        move |cap, view, shot| -> Result<(), String> {
            let (a_pre, b_pre, boundary_pre) = cap;
            let a_post = read_bits(&view, &a_after_add, shot);
            let b_post = read_bits(&view, &b_after_add, shot);
            let boundary_post = read_bits(&view, &boundary_after_add, shot);
            if boundary_post != *boundary_pre {
                return Err(format!(
                    "boundary drifted on shot {shot}: {boundary_pre}->{boundary_post}"
                ));
            }
            if b_post != *b_pre {
                return Err(format!("quantum addend changed on shot {shot}: {b_pre:#x}->{b_post:#x}"));
            }
            let prefix = (*boundary_pre as usize).min(n);
            let high_mask = !((1u64 << prefix) - 1) & ((1u64 << n) - 1);
            if (a_post & high_mask) != (*a_pre & high_mask) {
                return Err(format!(
                    "quantum-boundary carry escaped above prefix {prefix} on shot {shot}: a_pre={a_pre:#x} a_post={a_post:#x}"
                ));
            }
            Ok(())
        },
    );
    controlled_prefix_add_by_boundary(&mut quantum, &boundary, &a, &b, true);
    let a_final_check: Vec<&QReg> = a.iter().collect();
    let b_final_check: Vec<&QReg> = b.iter().collect();
    let boundary_final_check: Vec<&QReg> = boundary.iter().collect();
    quantum.contract_pop_and_check::<(u64, u64, u64), _>(
        "cursor_pack.quantum.roundtrip_pre",
        move |cap, view, shot| -> Result<(), String> {
            let (a_pre, b_pre, boundary_pre) = cap;
            let a_now = read_bits(&view, &a_final_check, shot);
            let b_now = read_bits(&view, &b_final_check, shot);
            let boundary_now = read_bits(&view, &boundary_final_check, shot);
            if a_now != *a_pre || b_now != *b_pre || boundary_now != *boundary_pre {
                return Err(format!(
                    "quantum-boundary reverse failed on shot {shot}: pre=({a_pre:#x},{b_pre:#x},{boundary_pre}) now=({a_now:#x},{b_now:#x},{boundary_now})"
                ));
            }
            Ok(())
        },
    );
    quantum.assert_phase_clean();
    let mut quantum_outs = a;
    quantum_outs.extend(b);
    quantum_outs.extend(boundary);
    let _ = quantum.destroy_sim(quantum_outs);

    eprintln!(
        "trailmix_cursor_pack_boundary_canary: PASS fixed_peak={} fixed_ops={} fixed_ccx={} quantum_peak={} quantum_ops={} quantum_ccx={}",
        fixed_peak,
        fixed_ops,
        fixed_ccx,
        quantum.peak_qubits,
        quantum.total_ops(),
        quantum.ccx_count(),
    );
}

pub fn build_trailmix_shrunken_pz_ops() -> Vec<Op> {
    let n = 256usize;
    let mut circ = TmCircuit::new();
    circ.disable_sim_contracts();
    circ.set_max_qubit_peak(1300);
    let ops_cap = std::env::var("TRAILMIX_BRIDGE_OPS_CAP")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(140_000_000);
    circ.set_ops_cap(ops_cap);

    let mut tx: Vec<_> = (0..n)
        .map(|i| circ.alloc_qreg(&format!("tx[{i}]")))
        .collect();
    let mut ty: Vec<_> = (0..n)
        .map(|i| circ.alloc_qreg(&format!("ty[{i}]")))
        .collect();
    let ox: Vec<_> = (0..n).map(|_| circ.alloc_input_bit()).collect();
    let oy: Vec<_> = (0..n).map(|_| circ.alloc_input_bit()).collect();

    super::trailmix_port::ec::point_add::ec_add_inplace_shrunken_pz(
        &mut circ, &mut tx, &mut ty, &ox, &oy,
    );
    circ.flush_pending_frees();

    assert_eq!(
        circ.ops_truncated,
        0,
        "TrailMix bridge op buffer truncated {} ops; raise TRAILMIX_BRIDGE_OPS_CAP above {}",
        circ.ops_truncated,
        circ.total_ops(),
    );
    eprintln!(
        "trailmix_bridge: total_ops={} peak_qubits={} total_qubits={} live_qubits={} total_bits={}",
        circ.total_ops(),
        circ.peak_qubits,
        circ.total_qubits(),
        circ.live_qubits(),
        circ.total_bits(),
    );
    if std::env::var("TRAILMIX_BRIDGE_PEAK_HIST").ok().as_deref() == Some("1") {
        print_peak_histogram(&circ);
    }
    if std::env::var("TRAILMIX_BRIDGE_PEAK_DETAIL").ok().as_deref() == Some("1") {
        print_peak_detail(&circ);
    }
    if std::env::var("TRAILMIX_BRIDGE_STATS_ONLY").ok().as_deref() == Some("1") {
        std::process::exit(0);
    }

    let mut out: Vec<_> = std::mem::take(&mut tx);
    out.extend(std::mem::take(&mut ty));
    let out = circ.defragment(out);

    circ.register(0);
    for q in &out[..n] {
        circ.append_qreg(q, 0);
    }
    circ.register(1);
    for q in &out[n..2 * n] {
        circ.append_qreg(q, 1);
    }
    circ.register(2);
    for b in &ox {
        circ.append_bit(*b, 2);
    }
    circ.register(3);
    for b in &oy {
        circ.append_bit(*b, 3);
    }

    let _public_output_ids: Vec<_> = out.into_iter().map(|q| q.detach()).collect();

    materialize_ops_low_peak(circ)
}

fn print_peak_histogram(circ: &TmCircuit) {
    let mut buckets = std::collections::BTreeMap::<String, usize>::new();
    for tag in &circ.peak_live_tags {
        let leaf = tag.rsplit('/').next().unwrap_or(tag);
        let class = leaf.split('[').next().unwrap_or(leaf).to_string();
        *buckets.entry(class).or_default() += 1;
    }
    let mut buckets = buckets.into_iter().collect::<Vec<_>>();
    buckets.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

    eprintln!("trailmix_bridge peak_hist top:");
    for (class, count) in buckets.iter().take(24) {
        eprintln!("  {class:32} {count}");
    }
}

fn print_peak_detail(circ: &TmCircuit) {
    let mut sections = std::collections::BTreeMap::<String, usize>::new();
    let mut classes = std::collections::BTreeMap::<String, usize>::new();
    for tag in &circ.peak_live_tags {
        let section = tag
            .rsplit_once('/')
            .map(|(section, _)| section)
            .unwrap_or("(root)");
        *sections.entry(section.to_string()).or_default() += 1;

        let class = tag
            .rsplit('/')
            .next()
            .unwrap_or(tag)
            .split('[')
            .next()
            .unwrap_or(tag)
            .to_string();
        *classes.entry(format!("{section}/{class}")).or_default() += 1;
    }

    let mut sections = sections.into_iter().collect::<Vec<_>>();
    sections.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    let mut classes = classes.into_iter().collect::<Vec<_>>();
    classes.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

    eprintln!("trailmix_bridge peak_detail section top:");
    for (section, count) in sections.iter().take(32) {
        eprintln!("  {count:4}  {section}");
    }
    eprintln!("trailmix_bridge peak_detail class top:");
    for (class, count) in classes.iter().take(64) {
        eprintln!("  {count:4}  {class}");
    }
}

fn materialize_ops_low_peak(circ: TmCircuit) -> Vec<Op> {
    let path = std::env::current_dir()
        .expect("current directory is available")
        .join(format!(".trailmix_bridge_ops_{}.tmp", std::process::id()));

    let mut count = 0usize;
    {
        let file = File::create(&path).expect("create TrailMix bridge temp op file");
        let mut writer = BufWriter::new(file);
        for op in circ.ops.iter().filter_map(|op| op.as_ref()) {
            write_op(&mut writer, &convert_op(op)).expect("write TrailMix bridge temp op");
            count += 1;
        }
        writer.flush().expect("flush TrailMix bridge temp op file");
    }
    drop(circ);

    let mut out = Vec::with_capacity(count);
    {
        let file = File::open(&path).expect("open TrailMix bridge temp op file");
        let mut reader = BufReader::new(file);
        for _ in 0..count {
            out.push(read_op(&mut reader).expect("read TrailMix bridge temp op"));
        }
    }
    let _ = remove_file(&path);
    out
}

fn write_op(mut writer: impl Write, op: &Op) -> std::io::Result<()> {
    writer.write_all(&(op.kind as u32).to_le_bytes())?;
    writer.write_all(&op.q_control2.0.to_le_bytes())?;
    writer.write_all(&op.q_control1.0.to_le_bytes())?;
    writer.write_all(&op.q_target.0.to_le_bytes())?;
    writer.write_all(&op.c_target.0.to_le_bytes())?;
    writer.write_all(&op.c_condition.0.to_le_bytes())?;
    writer.write_all(&op.r_target.0.to_le_bytes())?;
    Ok(())
}

fn read_op(mut reader: impl Read) -> std::io::Result<Op> {
    let kind = read_u32(&mut reader)?;
    Ok(Op {
        kind: operation_type_from_u32(kind),
        q_control2: QubitId(read_u64(&mut reader)?),
        q_control1: QubitId(read_u64(&mut reader)?),
        q_target: QubitId(read_u64(&mut reader)?),
        c_target: BitId(read_u64(&mut reader)?),
        c_condition: BitId(read_u64(&mut reader)?),
        r_target: RegisterId(read_u64(&mut reader)?),
    })
}

fn read_u32(mut reader: impl Read) -> std::io::Result<u32> {
    let mut bytes = [0u8; 4];
    reader.read_exact(&mut bytes)?;
    Ok(u32::from_le_bytes(bytes))
}

fn read_u64(mut reader: impl Read) -> std::io::Result<u64> {
    let mut bytes = [0u8; 8];
    reader.read_exact(&mut bytes)?;
    Ok(u64::from_le_bytes(bytes))
}

fn operation_type_from_u32(kind: u32) -> OperationType {
    match kind {
        0 => OperationType::Neg,
        1 => OperationType::Register,
        2 => OperationType::AppendToRegister,
        3 => OperationType::BitInvert,
        4 => OperationType::BitStore0,
        5 => OperationType::BitStore1,
        6 => OperationType::X,
        7 => OperationType::Z,
        8 => OperationType::CX,
        9 => OperationType::CZ,
        10 => OperationType::Swap,
        11 => OperationType::R,
        12 => OperationType::Hmr,
        13 => OperationType::CCX,
        14 => OperationType::CCZ,
        15 => OperationType::PushCondition,
        16 => OperationType::PopCondition,
        17 => OperationType::DebugPrint,
        _ => panic!("unknown operation kind {kind}"),
    }
}

fn convert_op(op: &TmOp) -> Op {
    let mut out = Op::empty();
    match *op {
        TmOp::Register(r) => {
            out.kind = OperationType::Register;
            out.r_target = RegisterId(r as u64);
        }
        TmOp::AppendQubit(q, r) => {
            out.kind = OperationType::AppendToRegister;
            out.q_target = QubitId(q as u64);
            out.r_target = RegisterId(r as u64);
        }
        TmOp::AppendBit(b, r) => {
            out.kind = OperationType::AppendToRegister;
            out.c_target = BitId(b as u64);
            out.r_target = RegisterId(r as u64);
        }
        TmOp::X(q) => {
            out.kind = OperationType::X;
            out.q_target = QubitId(q as u64);
        }
        TmOp::Z(q) => {
            out.kind = OperationType::Z;
            out.q_target = QubitId(q as u64);
        }
        TmOp::Cx(c, t) => {
            out.kind = OperationType::CX;
            out.q_control1 = QubitId(c as u64);
            out.q_target = QubitId(t as u64);
        }
        TmOp::Cz(a, b) => {
            out.kind = OperationType::CZ;
            out.q_control1 = QubitId(a as u64);
            out.q_target = QubitId(b as u64);
        }
        TmOp::Ccx(a, b, c) => {
            out.kind = OperationType::CCX;
            out.q_control2 = QubitId(a as u64);
            out.q_control1 = QubitId(b as u64);
            out.q_target = QubitId(c as u64);
        }
        TmOp::Ccz(a, b, c) => {
            out.kind = OperationType::CCZ;
            out.q_control2 = QubitId(a as u64);
            out.q_control1 = QubitId(b as u64);
            out.q_target = QubitId(c as u64);
        }
        TmOp::Swap(a, b) => {
            out.kind = OperationType::Swap;
            out.q_control1 = QubitId(a as u64);
            out.q_target = QubitId(b as u64);
        }
        TmOp::Hmr(q, b) => {
            out.kind = OperationType::Hmr;
            out.q_target = QubitId(q as u64);
            out.c_target = BitId(b as u64);
        }
        TmOp::R(q) => {
            out.kind = OperationType::R;
            out.q_target = QubitId(q as u64);
        }
        TmOp::Neg => {
            out.kind = OperationType::Neg;
        }
        TmOp::PushCondition(b) => {
            out.kind = OperationType::PushCondition;
            out.c_condition = BitId(b as u64);
        }
        TmOp::PopCondition => {
            out.kind = OperationType::PopCondition;
        }
        TmOp::BitInvert(b) => {
            out.kind = OperationType::BitInvert;
            out.c_target = BitId(b as u64);
        }
        TmOp::BitStore0(b) => {
            out.kind = OperationType::BitStore0;
            out.c_target = BitId(b as u64);
        }
        TmOp::BitStore1(b) => {
            out.kind = OperationType::BitStore1;
            out.c_target = BitId(b as u64);
        }
    }
    out
}
