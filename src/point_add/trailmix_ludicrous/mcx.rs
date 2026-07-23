
use super::{B, BExt};
use crate::circuit::{QubitId};
use std::sync::atomic::{AtomicU8, Ordering};

fn mbu_clear_and(circ: &mut B, t: &QubitId, c0: &QubitId, c1: &QubitId) {
    let bit = circ.alloc_bit();
    circ.hmr(*t, bit);
    circ.cz_if_bit(*c0, *c1, bit);
    circ.zero_and_free(*t);
}

// E284 (TLM_KG_INC_VENT=1): replace the reverse-pass AND-uncompute Toffolis in the KG
// increment with Gidney measurement-based uncomputation (mbu_clear_and, Clifford), bit-exact.
// Only ancillae dead after their uncompute are vented; live recursive temp_target toggles stay Ccx.
static KG_INC_VENT_FLAG: AtomicU8 = AtomicU8::new(2);
fn kg_inc_vent_enabled() -> bool {
    let c = KG_INC_VENT_FLAG.load(Ordering::Relaxed);
    if c != 2 {
        return c == 1;
    }
    let on = matches!(std::env::var("TLM_KG_INC_VENT").ok().as_deref(), Some("1"));
    KG_INC_VENT_FLAG.store(u8::from(on), Ordering::Relaxed);
    on
}

fn kg_get_layer_id(x: usize) -> usize {
    let mut layer_id = 0usize;
    let mut s = 0usize;
    while s <= x {
        s += (1usize << layer_id) + 1;
        layer_id += 1;
    }
    layer_id - 1
}

fn kg_start_layer(layer_id: usize) -> usize {
    let mut s = 0usize;
    for i in 0..layer_id {
        s += (1usize << i) + 1;
    }
    s
}

#[must_use]
pub fn kg_prefix_ancilla_count(n: usize) -> usize {
    if n <= 1 {
        return 0;
    }
    let targets_len = kg_get_layer_id(n - 1) + 1;
    if targets_len <= 2 {
        1
    } else {
        2 + kg_prefix_ancilla_count(targets_len)
    }
}

fn kg_apply_prefix_controlled_x(circ: &mut B, ctrls: &[&QubitId], target: &QubitId) {
    match ctrls {
        [] => circ.x(*target),
        [c] => circ.cx(**c, *target),
        [a, b] => circ.ccx(**a, **b, *target),
        _ => panic!("kg_apply_prefix_controlled_x: expected <=2 ctrls, got {}", ctrls.len()),
    }
}

fn kg_anc_index(len: usize, idx: isize) -> usize {
    if idx >= 0 {
        idx as usize
    } else {
        (len as isize + idx) as usize
    }
}

#[derive(Clone, Copy)]
enum KgPrefixOp<'a> {
    X(&'a QubitId),
    Ccx(&'a QubitId, &'a QubitId, &'a QubitId),
}

impl KgPrefixOp<'_> {
    fn emit(self, circ: &mut B) {
        match self {
            KgPrefixOp::X(q) => circ.x(*q),
            KgPrefixOp::Ccx(a, b, t) => circ.ccx(*a, *b, *t),
        }
    }
}

#[derive(Clone)]
struct KgPrefixLayer<'a> {
    ctrls: Vec<&'a QubitId>,
    ops: Vec<KgPrefixOp<'a>>,
}

fn kg_get_layers_for_prefix_and<'a>(
    q: &[&'a QubitId],
    inp_anc: &[&'a QubitId],
) -> Vec<KgPrefixLayer<'a>> {
    assert!(!q.is_empty(), "kg_get_layers_for_prefix_and: q must be non-empty");
    if q.len() == 1 {
        return vec![
            KgPrefixLayer { ctrls: Vec::new(), ops: Vec::new() },
            KgPrefixLayer { ctrls: vec![q[0]], ops: Vec::new() },
        ];
    }
    assert!(
        inp_anc.len() >= kg_prefix_ancilla_count(q.len()),
        "kg_get_layers_for_prefix_and: need {} ancillae for n={}, got {}",
        kg_prefix_ancilla_count(q.len()),
        q.len(),
        inp_anc.len(),
    );

    let n = q.len();
    let n_layers = kg_get_layer_id(q.len() - 1);
    let mut ret = vec![KgPrefixLayer { ctrls: Vec::new(), ops: Vec::new() }];
    let mut targets: Vec<&'a QubitId> = Vec::new();
    let mut anc: Vec<&'a QubitId> = vec![inp_anc[0]];

    for layer_id in 0..=n_layers {
        let st = kg_start_layer(layer_id);
        let en = n.min(kg_start_layer(layer_id + 1));

        let mut layer_ctrls = targets.clone();
        layer_ctrls.push(q[st]);
        ret.push(KgPrefixLayer { ctrls: layer_ctrls, ops: Vec::new() });

        for i in (st + 1)..en {
            let offset = i - st;
            let anc_len = anc.len();
            let q0 = q[i];
            let (q1, t) = if offset == 1 {
                (q[i - 1], anc[kg_anc_index(anc_len, -1)])
            } else {
                (
                    anc[kg_anc_index(anc_len, -(offset as isize - 1))],
                    anc[kg_anc_index(anc_len, -(offset as isize))],
                )
            };
            let mut ops = Vec::new();
            if std::ptr::eq(t, inp_anc[0]) {
                ops.push(KgPrefixOp::Ccx(q0, q1, t));
            } else {
                ops.push(KgPrefixOp::X(t));
                ops.push(KgPrefixOp::Ccx(q0, q1, t));
            }
            let mut ctrls = targets.clone();
            ctrls.push(t);
            ret.push(KgPrefixLayer { ctrls, ops });
        }

        let layer_len = en - st;
        let push_idx = kg_anc_index(anc.len(), 1 - layer_len as isize);
        targets.push(anc[push_idx]);

        let slice_start = kg_anc_index(anc.len(), 2 - layer_len as isize);
        let mut next_anc = anc[slice_start..].to_vec();
        next_anc.extend(q[st..en].iter());
        anc = next_anc;
    }

    if targets.len() <= 2 {
        return ret;
    }

    ret.push(KgPrefixLayer { ctrls: Vec::new(), ops: Vec::new() });
    let target_prefix_layers = kg_get_layers_for_prefix_and(&targets, &inp_anc[2..]);
    for layer_id in 1..=n_layers {
        let st = kg_start_layer(layer_id);
        let en = n.min(kg_start_layer(layer_id + 1));
        let target_prefix_targets = target_prefix_layers[layer_id].ctrls.clone();
        let ops_to_add = target_prefix_layers[layer_id].ops.clone();
        ret[st + 1].ops.extend_from_slice(&ops_to_add);

        let temp_target = if target_prefix_targets.len() == 1 {
            target_prefix_targets[0]
        } else {
            assert_eq!(target_prefix_targets.len(), 2);
            ret[st + 1].ops.push(KgPrefixOp::Ccx(
                target_prefix_targets[0],
                target_prefix_targets[1],
                inp_anc[1],
            ));
            inp_anc[1]
        };

        for i in st..en {
            let local = *ret[i + 1].ctrls.last().expect("empty local ctrl");
            ret[i + 1].ctrls = vec![temp_target, local];
        }

        if target_prefix_targets.len() == 2 {
            ret[en + 1].ops.push(KgPrefixOp::Ccx(
                target_prefix_targets[0],
                target_prefix_targets[1],
                temp_target,
            ));
        }
    }

    ret
}

fn xor_and_of_khattar_gidney_refs(circ: &mut B, bits: &[&QubitId], target: &QubitId) {
    match bits.len() {
        0 => {
            circ.x(*target);
            return;
        }
        1 => {
            circ.cx(*bits[0], *target);
            return;
        }
        2 => {
            circ.ccx(*bits[0], *bits[1], *target);
            return;
        }
        _ => {}
    }

    let anc_owned: Vec<QubitId> = (0..kg_prefix_ancilla_count(bits.len()))
        .map(|_| circ.alloc_qubit())
        .collect();
    let anc_refs: Vec<&QubitId> = anc_owned.iter().collect();
    let layers = kg_get_layers_for_prefix_and(bits, &anc_refs);

    for (i, layer) in layers.iter().enumerate() {
        if i > bits.len() {
            break;
        }
        for &op in &layer.ops {
            op.emit(circ);
        }
    }

    for (i, layer) in layers.iter().enumerate().rev() {
        if i > bits.len() {
            continue;
        }
        if i == bits.len() {
            kg_apply_prefix_controlled_x(circ, &layer.ctrls, target);
        }
        for &op in layer.ops.iter().rev() {
            op.emit(circ);
        }
    }
    drop(layers);
    drop(anc_refs);
    for q in anc_owned {
        circ.zero_and_free(q);
    }
}

pub fn mcx_clean_k(circ: &mut B, ctrls: &[&QubitId], target: &QubitId) {
    match ctrls.len() {
        0 => circ.x(*target),
        1 => circ.cx(*ctrls[0], *target),
        2 => circ.ccx(*ctrls[0], *ctrls[1], *target),
        3 => {
            let t = circ.alloc_qubit();
            circ.ccx(*ctrls[0], *ctrls[1], t);
            circ.ccx(t, *ctrls[2], *target);
            mbu_clear_and(circ, &t, ctrls[0], ctrls[1]);
        }
        4 => {
            let t01 = circ.alloc_qubit();
            let t23 = circ.alloc_qubit();
            circ.ccx(*ctrls[0], *ctrls[1], t01);
            circ.ccx(*ctrls[2], *ctrls[3], t23);
            circ.ccx(t01, t23, *target);
            mbu_clear_and(circ, &t23, ctrls[2], ctrls[3]);
            mbu_clear_and(circ, &t01, ctrls[0], ctrls[1]);
        }
        5 => {
            let t01 = circ.alloc_qubit();
            let t23 = circ.alloc_qubit();
            let t0123 = circ.alloc_qubit();
            circ.ccx(*ctrls[0], *ctrls[1], t01);
            circ.ccx(*ctrls[2], *ctrls[3], t23);
            circ.ccx(t01, t23, t0123);
            circ.ccx(t0123, *ctrls[4], *target);
            mbu_clear_and(circ, &t0123, &t01, &t23);
            mbu_clear_and(circ, &t23, ctrls[2], ctrls[3]);
            mbu_clear_and(circ, &t01, ctrls[0], ctrls[1]);
        }
        _ => {
            xor_and_of_khattar_gidney_refs(circ, ctrls, target);
        }
    }
}

pub fn inc_khattar_gidney(circ: &mut B, a: &[QubitId]) {
    let refs: Vec<&QubitId> = a.iter().collect();
    inc_khattar_gidney_refs_inner(circ, &refs, false);
}

pub fn cinc_khattar_gidney(circ: &mut B, a: &[QubitId], ctrl: &QubitId) {
    if a.is_empty() {
        return;
    }
    let mut combined: Vec<&QubitId> = Vec::with_capacity(a.len() + 1);
    combined.push(ctrl);
    combined.extend(a.iter());
    inc_khattar_gidney_refs_inner(circ, &combined, true);
}

fn inc_khattar_gidney_refs_inner(circ: &mut B, a: &[&QubitId], skip_lsb_x: bool) {
    let n = a.len();
    if n == 0 {
        return;
    }
    if n == 1 {
        if !skip_lsb_x {
            circ.x(*a[0]);
        }
        return;
    }

    let anc_owned: Vec<QubitId> = (0..kg_prefix_ancilla_count(n - 1))
        .map(|_| circ.alloc_qubit())
        .collect();
    let anc_refs: Vec<&QubitId> = anc_owned.iter().collect();
    let layers = kg_get_layers_for_prefix_and(&a[..n - 1], &anc_refs);

    for layer in &layers {
        for &op in &layer.ops {
            op.emit(circ);
        }
    }
    if !kg_inc_vent_enabled() {
        for (i, layer) in layers.iter().enumerate().rev() {
            if i < n && !(i == 0 && skip_lsb_x) {
                kg_apply_prefix_controlled_x(circ, &layer.ctrls, a[i]);
            }
            for &op in layer.ops.iter().rev() {
                op.emit(circ);
            }
        }
        drop(layers);
        drop(anc_refs);
        for q in anc_owned {
            circ.zero_and_free(q);
        }
        return;
    }

    // --- Vented reverse pass (E284, TLM_KG_INC_VENT=1) ---
    #[derive(Clone, Copy)]
    enum Step {
        X0(QubitId),
        X1(QubitId, QubitId),
        X2(QubitId, QubitId, QubitId),
        Xanc(QubitId),
        Uncmp(QubitId, QubitId, QubitId),
    }
    let mut plan: Vec<Step> = Vec::new();
    for (i, layer) in layers.iter().enumerate().rev() {
        if i < n && !(i == 0 && skip_lsb_x) {
            match layer.ctrls.as_slice() {
                [] => plan.push(Step::X0(*a[i])),
                [c] => plan.push(Step::X1(**c, *a[i])),
                [x, y] => plan.push(Step::X2(**x, **y, *a[i])),
                _ => panic!("inc_khattar_gidney vent: >2 prefix ctrls"),
            }
        }
        for &op in layer.ops.iter().rev() {
            match op {
                KgPrefixOp::X(t) => plan.push(Step::Xanc(*t)),
                KgPrefixOp::Ccx(x, y, t) => plan.push(Step::Uncmp(*x, *y, *t)),
            }
        }
    }
    drop(layers);
    drop(anc_refs);

    fn touches(s: &Step) -> [Option<u64>; 3] {
        match *s {
            Step::X0(t) => [Some(t.0), None, None],
            Step::X1(c, t) => [Some(c.0), Some(t.0), None],
            Step::X2(x, y, t) => [Some(x.0), Some(y.0), Some(t.0)],
            Step::Xanc(t) => [Some(t.0), None, None],
            Step::Uncmp(x, y, t) => [Some(x.0), Some(y.0), Some(t.0)],
        }
    }
    let mut occ: std::collections::HashMap<u64, Vec<usize>> = std::collections::HashMap::new();
    for (idx, s) in plan.iter().enumerate() {
        for q in touches(s).into_iter().flatten() {
            occ.entry(q).or_default().push(idx);
        }
    }
    let mut skip = vec![false; plan.len()];
    let mut vent_pure = vec![false; plan.len()];
    let mut vent_xc = vec![false; plan.len()];
    let mut vented_anc: std::collections::HashSet<u64> = std::collections::HashSet::new();
    for k in 0..plan.len() {
        if let Step::Uncmp(_, _, t) = plan[k] {
            let after: Vec<usize> = occ
                .get(&t.0)
                .map(|v| v.iter().copied().filter(|&j| j > k).collect())
                .unwrap_or_default();
            if after.is_empty() {
                vent_pure[k] = true;
                vented_anc.insert(t.0);
            } else if after.len() == 1
                && matches!(plan[after[0]], Step::Xanc(tt) if tt.0 == t.0)
            {
                vent_xc[k] = true;
                skip[after[0]] = true;
                vented_anc.insert(t.0);
            }
        }
    }
    for k in 0..plan.len() {
        if skip[k] {
            continue;
        }
        match plan[k] {
            Step::X0(t) => circ.x(t),
            Step::X1(c, t) => circ.cx(c, t),
            Step::X2(x, y, t) => circ.ccx(x, y, t),
            Step::Xanc(t) => circ.x(t),
            Step::Uncmp(x, y, t) => {
                if vent_pure[k] {
                    mbu_clear_and(circ, &t, &x, &y);
                } else if vent_xc[k] {
                    circ.x(t);
                    mbu_clear_and(circ, &t, &x, &y);
                } else {
                    circ.ccx(x, y, t);
                }
            }
        }
    }
    for q in anc_owned {
        if !vented_anc.contains(&q.0) {
            circ.zero_and_free(q);
        }
    }
}
