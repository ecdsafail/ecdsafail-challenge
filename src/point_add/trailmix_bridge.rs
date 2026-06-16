use super::{BitId, Op, OperationType, QubitId, RegisterId};
use std::fs::{remove_file, File};
use std::io::{BufReader, BufWriter, Read, Write};

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
    if std::env::var("TRAILMIX_BRIDGE_PEAK_HIST").ok().as_deref() == Some("1") {
        print_peak_histogram(&circ);
    }
    if std::env::var("TRAILMIX_BRIDGE_PEAK_DETAIL").ok().as_deref() == Some("1") {
        print_peak_detail(&circ);
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
        let section = tag.rsplit_once('/').map(|(section, _)| section).unwrap_or("(root)");
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
