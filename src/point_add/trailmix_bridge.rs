use super::{BitId, Op, OperationType, QubitId, RegisterId};

type TmCircuit = super::trailmix_port::circuit::Circuit;
type TmOp = super::trailmix_port::circuit::Op;

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

    circ.ops
        .iter()
        .filter_map(|op| op.as_ref())
        .map(convert_op)
        .collect()
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
