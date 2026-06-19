//! Multi-controlled X with log*(k) clean ancillae, built on the Khattar-Gidney
//! prefix-AND ladder. Used by the fold tail of the EC point-add circuit, which
//! propagates the fold carry into the high bits via a cascade of `mcx_clean_k`
//! calls.

use super::{B, BExt};
use crate::circuit::{QubitId};

/// Clear `t = c0 AND c1` back to |0> by measurement (HMR + cz_if_bit, 0
/// Toffoli). `c0`,`c1` must be alive/unmodified since the AND.
fn mbu_clear_and(circ: &mut B, t: &QubitId, c0: &QubitId, c1: &QubitId) {
    let bit = circ.alloc_bit();
    circ.hmr(*t, bit);
    circ.cz_if_bit(*c0, *c1, bit);
    circ.zero_and_free(*t);
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

/// Clean-ancilla count for the KG prefix layer decomposition.
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

/// `target ^= AND(bits)` via the KG prefix-AND ladder (2k-3 Toffoli, log*(k)
/// clean ancillae).
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

/// Clean `target ^= AND(ctrls)`; `ctrls` restored. log*(k) clean ancillae.
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
