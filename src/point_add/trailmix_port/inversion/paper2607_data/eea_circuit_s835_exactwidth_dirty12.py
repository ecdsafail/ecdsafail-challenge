import hashlib
import json
from functools import lru_cache
from pathlib import Path
from typing import Literal, Optional, Sequence

from qiskit import QuantumCircuit, QuantumRegister
from qiskit.circuit import Gate, Qubit

import eea_circuit_updated as _e

C_EEA = _e.C_EEA
N_CONFIG = _e.N_CONFIG
paper_len_width = _e.paper_len_width
paper_shift_width = _e.paper_shift_width
Nmax_steps = _e.Nmax_steps
active_windows = _e.active_windows
get_n_config = getattr(_e, "get_n_config")
set_measurement_uncompute = _e.set_measurement_uncompute
count_circuit_ops_recursive = getattr(_e, "count_circuit_ops_recursive", None)

_CERTIFIED_WINDOW_SHA256 = "3e1961f5550249604bf044edb65f1d1bc403ed75bd7178e283685ddb4f3cb880"
_CERTIFIED_WINDOW_PATH = Path(__file__).with_name("active_windows_1616.json")
_certified_window_bytes = _CERTIFIED_WINDOW_PATH.read_bytes()
if hashlib.sha256(_certified_window_bytes).hexdigest() != _CERTIFIED_WINDOW_SHA256:
    raise RuntimeError("secp256k1 active-window certificate hash mismatch")
_certified_window_table = json.loads(_certified_window_bytes)
if (
    _certified_window_table.get("schema") != "luo-secp256k1-active-windows-v2"
    or len(_certified_window_table.get("rows", ())) != 1616
):
    raise RuntimeError("invalid secp256k1 active-window certificate")
_CERTIFIED_WINDOW_ROWS = tuple(row["safe"] for row in _certified_window_table["rows"])

LT_WIDTH = 8
LQ_WIDTH = 9
LS_WIDTH = 9
LRP_WIDTH = 8
LS_MODULUS = 259
LS_ZERO = LS_MODULUS - 1
LRP_ZERO = (1 << LRP_WIDTH) - 1
CLEAN_AUX_SIZE = 12
DIRTY_PASSENGER_SIZE = 10


def __getattr__(name: str):
    return getattr(_e, name)


def _tight_unary_depth_for_labels(labels: Sequence[int]) -> int:
    labels = sorted(set(labels))
    if len(labels) <= 1:
        return 0
    bit = _e._split_bit(labels)
    z = [x for x in labels if ((x >> bit) & 1) == 0]
    o = [x for x in labels if ((x >> bit) & 1) == 1]
    return 1 + max(_tight_unary_depth_for_labels(z), _tight_unary_depth_for_labels(o))


def unary_iteration_tight(qc: QuantumCircuit, *, index_reg: Sequence[Qubit], labels: Sequence[int],
                          ctrl: Qubit, ancillas: Sequence[Qubit], leaf_fn, order: Literal["inc", "dec"] = "inc") -> None:
    labels = sorted(set(labels))
    if not labels:
        return
    need = _tight_unary_depth_for_labels(labels)
    if len(ancillas) < need:
        raise ValueError(f"tight unary iteration needs {need} ancillas, got {len(ancillas)}")
    def rec(sub_labels, g, depth):
        if len(sub_labels) == 1:
            leaf_fn(sub_labels[0], g); return
        b = _e._split_bit(sub_labels)
        z = [x for x in sub_labels if ((x >> b) & 1) == 0]
        o = [x for x in sub_labels if ((x >> b) & 1) == 1]
        h = ancillas[depth]
        _e._and_with_index_bit(qc, g, index_reg[b], h, 0)
        if order == "inc":
            rec(z, h, depth+1)
            qc.cx(g, h)
            rec(o, h, depth+1)
            qc.cx(g, h)
        else:
            qc.cx(g, h)
            rec(o, h, depth+1)
            qc.cx(g, h)
            rec(z, h, depth+1)
        _e._uncompute_and_with_index_bit(qc, g, index_reg[b], h, 0)
    rec(labels, ctrl, 0)


def dual_unary_iteration_tight(qc: QuantumCircuit, *, index_a: Sequence[Qubit], index_b: Sequence[Qubit], labels: Sequence[int],
                               ctrl_a: Qubit, ctrl_b: Qubit, ancillas_a: Sequence[Qubit], ancillas_b: Sequence[Qubit],
                               leaf_fn, order: Literal["inc", "dec"] = "inc") -> None:
    labels = sorted(set(labels))
    if not labels:
        return
    need = _tight_unary_depth_for_labels(labels)
    if len(ancillas_a) < need or len(ancillas_b) < need:
        raise ValueError(f"tight dual unary iteration needs {need} ancillas per endpoint")
    def rec(sub_labels, ga, gb, depth):
        if len(sub_labels) == 1:
            leaf_fn(sub_labels[0], ga, gb); return
        bit = _e._split_bit(sub_labels)
        z = [x for x in sub_labels if ((x >> bit) & 1) == 0]
        o = [x for x in sub_labels if ((x >> bit) & 1) == 1]
        ha = ancillas_a[depth]; hb = ancillas_b[depth]
        _e._and_with_index_bit(qc, ga, index_a[bit], ha, 0)
        _e._and_with_index_bit(qc, gb, index_b[bit], hb, 0)
        if order == "inc":
            rec(z, ha, hb, depth+1)
            qc.cx(ga, ha); qc.cx(gb, hb)
            rec(o, ha, hb, depth+1)
            qc.cx(gb, hb); qc.cx(ga, ha)
        else:
            qc.cx(ga, ha); qc.cx(gb, hb)
            rec(o, ha, hb, depth+1)
            qc.cx(gb, hb); qc.cx(ga, ha)
            rec(z, ha, hb, depth+1)
        _e._uncompute_and_with_index_bit(qc, gb, index_b[bit], hb, 0)
        _e._uncompute_and_with_index_bit(qc, ga, index_a[bit], ha, 0)
    rec(labels, ctrl_a, ctrl_b, 0)


def kg_prefix_ancilla_count(n: int) -> int:
    """Exact port of ``arith/khattar_gidney.rs::kg_prefix_ancilla_count``."""
    if n <= 1:
        return 0
    targets_len = _kg_get_layer_id(n - 1) + 1
    if targets_len <= 2:
        return 1
    return 2 + kg_prefix_ancilla_count(targets_len)


def _kg_get_layer_id(x: int) -> int:
    layer_id = 0
    start = 0
    while start <= x:
        start += (1 << layer_id) + 1
        layer_id += 1
    return layer_id - 1


def _kg_start_layer(layer_id: int) -> int:
    return sum((1 << i) + 1 for i in range(layer_id))


def _kg_get_layers_for_prefix_and(q: Sequence[Qubit], ancillas: Sequence[Qubit]):
    """Return the exact conditionally-clean KG layer schedule used by Rust."""
    q = list(q)
    ancillas = list(ancillas)
    if not q:
        raise ValueError("KG prefix input must be non-empty")
    if len(q) == 1:
        return [dict(ctrls=[], ops=[]), dict(ctrls=[q[0]], ops=[])]
    need = kg_prefix_ancilla_count(len(q))
    if len(ancillas) < need:
        raise ValueError(f"KG prefix needs {need} ancillas, got {len(ancillas)}")

    n = len(q)
    n_layers = _kg_get_layer_id(n - 1)
    layers = [dict(ctrls=[], ops=[])]
    targets: list[Qubit] = []
    anc = [ancillas[0]]

    for layer_id in range(n_layers + 1):
        start = _kg_start_layer(layer_id)
        end = min(n, _kg_start_layer(layer_id + 1))
        layers.append(dict(ctrls=targets + [q[start]], ops=[]))
        for i in range(start + 1, end):
            offset = i - start
            if offset == 1:
                q1, target = q[i - 1], anc[-1]
            else:
                q1, target = anc[-(offset - 1)], anc[-offset]
            ops = []
            if target is ancillas[0]:
                ops.append(("ccx", q[i], q1, target))
            else:
                ops.append(("x", target))
                ops.append(("ccx", q[i], q1, target))
            layers.append(dict(ctrls=targets + [target], ops=ops))

        layer_len = end - start
        targets.append(anc[1 - layer_len])
        anc = anc[2 - layer_len:] + q[start:end]

    if len(targets) <= 2:
        return layers

    layers.append(dict(ctrls=[], ops=[]))
    target_layers = _kg_get_layers_for_prefix_and(targets, ancillas[2:])
    for layer_id in range(1, n_layers + 1):
        start = _kg_start_layer(layer_id)
        end = min(n, _kg_start_layer(layer_id + 1))
        target_ctrls = list(target_layers[layer_id]["ctrls"])
        layers[start + 1]["ops"].extend(target_layers[layer_id]["ops"])
        if len(target_ctrls) == 1:
            temp_target = target_ctrls[0]
        elif len(target_ctrls) == 2:
            temp_target = ancillas[1]
            layers[start + 1]["ops"].append(
                ("ccx", target_ctrls[0], target_ctrls[1], temp_target)
            )
        else:
            raise AssertionError("KG recursive target prefix must expose one or two controls")
        for i in range(start, end):
            local = layers[i + 1]["ctrls"][-1]
            layers[i + 1]["ctrls"] = [temp_target, local]
        if len(target_ctrls) == 2:
            layers[end + 1]["ops"].append(
                ("ccx", target_ctrls[0], target_ctrls[1], temp_target)
            )
    return layers


def _kg_emit_op(qc: QuantumCircuit, op) -> None:
    if op[0] == "x":
        qc.x(op[1])
    elif op[0] == "ccx":
        qc.ccx(op[1], op[2], op[3])
    else:
        raise AssertionError(f"unknown KG op {op[0]}")


def _kg_emit_layers(qc: QuantumCircuit, layers, *, reverse: bool = False) -> None:
    layer_order = reversed(layers) if reverse else layers
    for layer in layer_order:
        op_order = reversed(layer["ops"]) if reverse else layer["ops"]
        for op in op_order:
            _kg_emit_op(qc, op)


def _kg_lowest_layer_touching(layers, changed: Sequence[Qubit]) -> Optional[int]:
    changed_ids = {id(q) for q in changed}
    for index, layer in enumerate(layers):
        for op in layer["ops"]:
            if any(id(q) in changed_ids for q in op[1:]):
                return index
    return None


def _kg_toggle_equality(qc: QuantumCircuit, *, base: Sequence[Qubit], c0: Qubit,
                        flag: Qubit, clean_temp: Qubit) -> None:
    controls = list(base) + [c0]
    if len(controls) == 1:
        qc.cx(controls[0], flag)
    elif len(controls) == 2:
        qc.ccx(controls[0], controls[1], flag)
    elif len(controls) == 3:
        qc.ccx(controls[0], controls[1], clean_temp)
        qc.ccx(clean_temp, controls[2], flag)
        qc.ccx(controls[0], controls[1], clean_temp)
    else:
        raise ValueError(f"KG equality expected at most three controls, got {len(controls)}")


def dual_unary_iteration_log_star(qc: QuantumCircuit, *,
                                  index_a: Sequence[Qubit], index_b: Sequence[Qubit],
                                  labels: Sequence[int], ancillas_a: Sequence[Qubit],
                                  ancillas_b: Sequence[Qubit], flag_a: Qubit,
                                  flag_b: Qubit, common_ctrl: Qubit, clean_temp: Qubit,
                                  leaf_fn, order: Literal["inc", "dec"] = "inc") -> None:
    """Dual exact KG unary iterator with synchronized Gray updates.

    Each callback sees cleanly materialized raw equality flags for both
    endpoints.  Prefix and equality ancillas, borrowed lanes, and endpoints
    are restored exactly on return.
    """
    labels = sorted(set(labels), reverse=(order == "dec"))
    if not labels:
        return
    if len(index_a) != len(index_b) or len(index_a) < 2:
        raise ValueError("dual KG iterator requires equal endpoint widths >= 2")
    n = len(index_a)
    # Fold the common control into each prefix input.  Keep it LAST so the
    # conditionally-clean KG schedule never borrows the shared Ctrl as a
    # target; both endpoint engines can then remain live simultaneously.
    # The prefix product is AND(c[n-1],...,c[1],Ctrl), while c[0] remains the
    # separate final control.
    need = kg_prefix_ancilla_count(n)
    if len(ancillas_a) < need or len(ancillas_b) < need:
        raise ValueError(f"dual KG iterator needs {need} ancillas per endpoint")

    def complement_for(index: Sequence[Qubit], value: int) -> None:
        for bit, lane in enumerate(index):
            if ((value >> bit) & 1) == 0:
                qc.x(lane)

    start = labels[0]
    complement_for(index_a, start)
    complement_for(index_b, start)
    bits_a = list(reversed(index_a))
    bits_b = list(reversed(index_b))
    prefix_a = bits_a[:-1] + [common_ctrl]
    prefix_b = bits_b[:-1] + [common_ctrl]
    layers_a = _kg_get_layers_for_prefix_and(prefix_a, ancillas_a[:need])
    layers_b = _kg_get_layers_for_prefix_and(prefix_b, ancillas_b[:need])
    for layers in (layers_a, layers_b):
        if any(op[-1] == common_ctrl for layer in layers for op in layer["ops"]):
            raise AssertionError("dual KG schedule must not target shared Ctrl")
    _kg_emit_layers(qc, layers_a)
    _kg_emit_layers(qc, layers_b)
    base_a = list(layers_a[len(prefix_a)]["ctrls"])
    base_b = list(layers_b[len(prefix_b)]["ctrls"])

    for position, label in enumerate(labels):
        _kg_toggle_equality(
            qc, base=base_a, c0=index_a[0], flag=flag_a, clean_temp=clean_temp,
        )
        _kg_toggle_equality(
            qc, base=base_b, c0=index_b[0], flag=flag_b, clean_temp=clean_temp,
        )
        leaf_fn(label, flag_a, flag_b)
        _kg_toggle_equality(
            qc, base=base_b, c0=index_b[0], flag=flag_b, clean_temp=clean_temp,
        )
        _kg_toggle_equality(
            qc, base=base_a, c0=index_a[0], flag=flag_a, clean_temp=clean_temp,
        )

        if position + 1 == len(labels):
            continue
        next_label = labels[position + 1]
        delta = label ^ next_label
        changed_a = [bits_a[n - 1 - bit] for bit in range(1, n) if (delta >> bit) & 1]
        changed_b = [bits_b[n - 1 - bit] for bit in range(1, n) if (delta >> bit) & 1]
        first_a = _kg_lowest_layer_touching(layers_a, changed_a)
        first_b = _kg_lowest_layer_touching(layers_b, changed_b)
        if first_b is not None:
            _kg_emit_layers(qc, layers_b[first_b:], reverse=True)
        if first_a is not None:
            _kg_emit_layers(qc, layers_a[first_a:], reverse=True)
        for bit in range(n):
            if (delta >> bit) & 1:
                qc.x(index_a[bit])
                qc.x(index_b[bit])
        if first_a is not None:
            _kg_emit_layers(qc, layers_a[first_a:])
        if first_b is not None:
            _kg_emit_layers(qc, layers_b[first_b:])

    _kg_emit_layers(qc, layers_b, reverse=True)
    _kg_emit_layers(qc, layers_a, reverse=True)
    complement_for(index_b, labels[-1])
    complement_for(index_a, labels[-1])


def _toggle_eq_const_under_ctrl_direct(qc: QuantumCircuit, *, endpoint: Sequence[Qubit], const: int, ctrl: Qubit, acc: Qubit, scratch: Sequence[Qubit]) -> None:
    # scratch supplies a temporary eq flag followed by mcx scratch.
    eq = scratch[0]
    pool = list(scratch[1:])
    _e.compute_eq_const(qc, endpoint, const, eq, pool)
    qc.ccx(ctrl, eq, acc)
    _e.compute_eq_const(qc, endpoint, const, eq, pool)


def _const_scratch(Scratch, width: int, carry: Qubit) -> list[Qubit]:
    # add_const_mod_2n expects width constant bits followed by one clean carry.
    return list(Scratch[:width]) + [carry]


def _controlled_adjacent_basis_swap(qc: QuantumCircuit, *, ctrl: Qubit,
                                    reg: Sequence[Qubit], a: int, b: int,
                                    scratch: Sequence[Qubit]) -> None:
    """Swap adjacent basis labels a/b under ctrl, restoring clean scratch."""
    diff = a ^ b
    if diff == 0 or diff & (diff - 1):
        raise ValueError("adjacent basis labels must differ in exactly one bit")
    target_bit = diff.bit_length() - 1
    controls = [ctrl]
    inverted: list[Qubit] = []
    for bit, qubit in enumerate(reg):
        if bit == target_bit:
            continue
        if ((a >> bit) & 1) == 0:
            qc.x(qubit)
            inverted.append(qubit)
        controls.append(qubit)
    _e.mcx_vchain(qc, controls, reg[target_bit], scratch)
    for qubit in reversed(inverted):
        qc.x(qubit)


def _controlled_basis_swap(qc: QuantumCircuit, *, ctrl: Qubit,
                           reg: Sequence[Qubit], a: int, b: int,
                           scratch: Sequence[Qubit]) -> None:
    """Exact controlled transposition of two computational-basis labels."""
    if a == b:
        return
    path = [a]
    current = a
    for bit in range(len(reg)):
        if ((a ^ b) >> bit) & 1:
            current ^= 1 << bit
            path.append(current)
    if path[-1] != b:
        raise AssertionError("basis-swap Gray path")
    edges = list(zip(path, path[1:]))
    for left, right in edges:
        _controlled_adjacent_basis_swap(
            qc, ctrl=ctrl, reg=reg, a=left, b=right, scratch=scratch,
        )
    for left, right in reversed(edges[:-1]):
        _controlled_adjacent_basis_swap(
            qc, ctrl=ctrl, reg=reg, a=left, b=right, scratch=scratch,
        )


def _controlled_zero_259_swap_linear(qc: QuantumCircuit, *, ctrl: Qubit,
                                     reg: Sequence[Qubit],
                                     scratch: Sequence[Qubit]) -> None:
    """Swap |0> and |259> with one high-control toggle, globally exactly.

    The difference word 259 has bits {0,1,8}.  Conjugating by
    x0 ^= x8; x1 ^= x8 maps it to the unit word 256, so the transposition
    needs one adjacent basis swap instead of a five-swap Gray palindrome.
    """
    if len(reg) != LS_WIDTH:
        raise ValueError("0/259 transposition requires a 9-bit register")
    qc.cx(reg[8], reg[0])
    qc.cx(reg[8], reg[1])
    _controlled_adjacent_basis_swap(
        qc, ctrl=ctrl, reg=reg, a=0, b=1 << 8, scratch=scratch,
    )
    qc.cx(reg[8], reg[1])
    qc.cx(reg[8], reg[0])


def inc_mod259_1ctrl(qc: QuantumCircuit, ctrl: Qubit,
                     reg: Sequence[Qubit], scratch: Sequence[Qubit]) -> None:
    """Controlled +1 on 0..258, extended to a permutation on all 9-bit words."""
    if len(reg) != LS_WIDTH:
        raise ValueError("mod-259 increment requires a 9-bit register")
    _e.inc_mod2n_1ctrl(qc, ctrl, list(reg), scratch[: LS_WIDTH - 1])
    _controlled_zero_259_swap_linear(qc, ctrl=ctrl, reg=reg, scratch=scratch)


def dec_mod259_1ctrl(qc: QuantumCircuit, ctrl: Qubit,
                     reg: Sequence[Qubit], scratch: Sequence[Qubit]) -> None:
    """Exact inverse of inc_mod259_1ctrl."""
    if len(reg) != LS_WIDTH:
        raise ValueError("mod-259 decrement requires a 9-bit register")
    _controlled_zero_259_swap_linear(qc, ctrl=ctrl, reg=reg, scratch=scratch)
    _e.dec_mod2n_1ctrl(qc, ctrl, list(reg), scratch[: LS_WIDTH - 1])


def _swap_zero_259_uncontrolled(qc: QuantumCircuit, reg: Sequence[Qubit],
                                one: Qubit, scratch: Sequence[Qubit]) -> None:
    """Swap basis labels 0 and 259, restoring a temporary constant-one bit."""
    qc.x(one)
    _controlled_zero_259_swap_linear(qc, ctrl=one, reg=reg, scratch=scratch)
    qc.x(one)


@lru_cache(maxsize=None)
def clean_c3x_mbu_gate() -> Gate:
    """Self-inverse C^3X with a clean temporary lowered by KMX HMR."""
    wires = QuantumRegister(5, "c3x")
    qc = QuantumCircuit(wires, name="CLEAN_C3X_MBU")
    qc.ccx(wires[0], wires[1], wires[4])
    qc.ccx(wires[2], wires[4], wires[3])
    qc.ccx(wires[0], wires[1], wires[4])
    return qc.to_gate()


def _dirty_c3x(qc: QuantumCircuit, a: Qubit, b: Qubit, c: Qubit, target: Qubit, dirty: Qubit) -> None:
    qc.append(clean_c3x_mbu_gate(), [a, b, c, target, dirty])


def _controlled_toffoli_dirty(qc: QuantumCircuit, ctrl: Qubit, a: Qubit, b: Qubit, target: Qubit, dirty: Qubit) -> None:
    _dirty_c3x(qc, ctrl, a, b, target, dirty)


def controlled_maj_dirty(qc: QuantumCircuit, ctrl: Qubit, a: Qubit, b: Qubit, c: Qubit, dirty: Qubit) -> None:
    qc.ccx(ctrl, a, b)
    qc.ccx(ctrl, a, c)
    _controlled_toffoli_dirty(qc, ctrl, c, b, a, dirty)


def controlled_uma_dirty(qc: QuantumCircuit, ctrl: Qubit, a: Qubit, b: Qubit, c: Qubit, dirty: Qubit) -> None:
    _controlled_toffoli_dirty(qc, ctrl, c, b, a, dirty)
    qc.ccx(ctrl, a, c)
    qc.ccx(ctrl, c, b)


def controlled_maj_inv_dirty(qc: QuantumCircuit, ctrl: Qubit, a: Qubit, b: Qubit, c: Qubit, dirty: Qubit) -> None:
    _controlled_toffoli_dirty(qc, ctrl, c, b, a, dirty)
    qc.ccx(ctrl, a, c)
    qc.ccx(ctrl, a, b)


def controlled_uma_inv_dirty(qc: QuantumCircuit, ctrl: Qubit, a: Qubit, b: Qubit, c: Qubit, dirty: Qubit) -> None:
    qc.ccx(ctrl, c, b)
    qc.ccx(ctrl, a, c)
    _controlled_toffoli_dirty(qc, ctrl, c, b, a, dirty)


def _apply_cell_dirty(qc: QuantumCircuit, mode: Literal["add", "sub"], pass_kind: Literal["first", "second"],
                      ctrl: Qubit, addend: Qubit, target: Qubit, carry: Qubit, dirty: Qubit) -> None:
    if mode == "add" and pass_kind == "first":
        controlled_maj_dirty(qc, ctrl, addend, target, carry, dirty)
    elif mode == "add" and pass_kind == "second":
        controlled_uma_dirty(qc, ctrl, addend, target, carry, dirty)
    elif mode == "sub" and pass_kind == "first":
        controlled_uma_inv_dirty(qc, ctrl, addend, target, carry, dirty)
    elif mode == "sub" and pass_kind == "second":
        controlled_maj_inv_dirty(qc, ctrl, addend, target, carry, dirty)
    else:
        raise ValueError("bad arithmetic cell mode/pass")


@lru_cache(maxsize=None)
def lc_swap_unary_gate(*, k: int, K: int, len_width: int, name: str = "LC_SWAP_S835_FAST") -> Gate:
    if k > K:
        raise ValueError("need k <= K")
    M = K - k + 1
    depth = _e.unary_depth(M)
    base = max(len_width, depth)
    scratch_size = base + 2
    Ctrl = QuantumRegister(1, "Ctrl")
    Direction = QuantumRegister(1, "Direction")
    Sign = QuantumRegister(1, "Sign")
    Work1 = QuantumRegister(M + 1, "Work1")
    l_t = QuantumRegister(len_width, "l_t")
    l_q = QuantumRegister(len_width, "l_q")
    Scratch = QuantumRegister(scratch_size, "Scratch")
    qc = _e._block_circuit(Ctrl, Direction, Sign, Work1, l_t, l_q, Scratch, name=name)
    carry = Scratch[base]
    direction_flag = Scratch[base + 1]
    cs = list(Scratch[:len_width]) + [carry]
    qc.append(_e.cuccaro_add_mod_2n_no_z_gate(len_width, name="ADD_lt_to_lq"), list(l_t) + list(l_q) + [carry])
    _e.add_const_mod_2n(qc, l_q, 3, cs)
    path = list(Scratch[:depth])
    def leaf(j: int, ej: Qubit) -> None:
        # Phase 2 inserts the next quotient bit at physical j.  Phase 3 removes
        # the current low quotient bit at physical j-1.  Direction (Phase1) is
        # retained by the caller, so this branch is exactly reversible.
        _e._and_with_index_bit(qc, ej, Direction[0], direction_flag, 0)
        _e.cswap_toffoli(qc, direction_flag, Sign[0], Work1[j - k + 1])
        qc.cx(ej, direction_flag)
        _e.cswap_toffoli(qc, direction_flag, Sign[0], Work1[j - k])
        qc.cx(ej, direction_flag)
        _e._uncompute_and_with_index_bit(qc, ej, Direction[0], direction_flag, 0)
    unary_iteration_tight(qc, index_reg=l_q, labels=list(range(k, K + 1)), ctrl=Ctrl[0], ancillas=path, leaf_fn=leaf, order="inc")
    _e.sub_const_mod_2n(qc, l_q, 3, cs)
    qc.append(_e.cuccaro_sub_mod_2n_no_z_gate(len_width, name="SUB_lt_from_lq"), list(l_t) + list(l_q) + [carry])
    return _e._finalize_block(qc)


@lru_cache(maxsize=None)
def lc_interval_addsub_unary_gate(*, n: int, k: int, K: int, len_width: int, shift_width: int,
                                  mode: Literal["add", "sub"], sign_update: bool,
                                  target: Literal["work1", "work2"], name: str) -> Gate:
    if k > K:
        raise ValueError("need k <= K")
    M = K - k + 1
    endpoint_width = max(len_width, shift_width)
    # Decode the complete interval.  Splitting a 2^d+1 interval into a 2^d
    # unary tree plus a special top label is unsound unless the tree is also
    # conditioned on the omitted high bit: the top endpoint otherwise aliases
    # label zero.  The full tree costs one additional path qubit per endpoint
    # and is injective over every in-range endpoint.
    labels_all_abs = list(range(k, K + 1))
    rel_count = len(labels_all_abs)
    labels_main = list(range(rel_count))
    top_special = False
    top_rel = rel_count - 1
    depth = _tight_unary_depth_for_labels(labels_main)
    # Layout note:
    #   anc_a/anc_b occupy the first 2*depth wires and are used only by
    #   the unary endpoint scans.  Endpoint affine transforms need
    #   endpoint_width scratch wires plus a carry.  For late steps the unary
    #   depth can be smaller than endpoint_width; placing carry immediately
    #   after the unary paths would then alias it with the constant-adder
    #   scratch.  We therefore place carry/acc/cell_pool after the larger of
    #   the unary-scratch region and the endpoint-transform scratch region.
    base = max(2 * depth, endpoint_width)
    scratch_size = base + 3
    Ctrl = QuantumRegister(1, "Ctrl")
    Sign = QuantumRegister(1, "Sign")
    Work1 = QuantumRegister(M, "Work1")
    Work2 = QuantumRegister(M, "Work2")
    l_t = QuantumRegister(len_width, "l_t")
    l_q = QuantumRegister(len_width, "l_q")
    l_s = QuantumRegister(shift_width, "l_s")
    Scratch = QuantumRegister(scratch_size, "Scratch")
    qc = _e._block_circuit(Ctrl, Sign, Work1, Work2, l_t, l_q, l_s, Scratch, name=name)
    anc_a = list(Scratch[:depth])
    anc_b = list(Scratch[depth:2*depth])
    carry = Scratch[base]
    acc = Scratch[base + 1]
    cell_pool = [Scratch[base + 2]]
    # Top-special equality controls reuse one clean unary-path wire as the
    # one-hot flag.  The remaining clean paths plus cell_pool form its MCX
    # scratch; this keeps the n=256 block within the 20-qubit shared pool.
    top_flag = Scratch[0]
    eq_scratch = [Scratch[base + 2]] + [q for q in Scratch[:base] if q != top_flag]
    cs = _const_scratch(Scratch, endpoint_width, carry)
    # Prepare L=(ell_t-1)+(ell_q-1)+4 and R=n+2-(ell_s-1).
    qc.append(_e.cuccaro_add_mod_2n_no_z_gate(len_width, name="ADD_lt_to_lq"), list(l_t) + list(l_q) + [carry])
    _e.add_const_mod_2n(qc, l_q, 4, cs[:len_width] + [carry])
    _e.const_minus_inplace(qc, l_s, n + 2, cs[:shift_width] + [carry])
    # Convert absolute endpoints to relative offsets in [0, K-k].
    _e.sub_const_mod_2n(qc, l_q, k, cs[:len_width] + [carry])
    _e.sub_const_mod_2n(qc, l_s, k, cs[:shift_width] + [carry])
    def qpair(j: int) -> tuple[Qubit, Qubit]:
        j_abs = k + j
        idx = j_abs - k
        if target == "work1":
            return Work2[idx], Work1[idx]
        if target == "work2":
            return Work1[idx], Work2[idx]
        raise ValueError("bad target")
    def leaf_first(j: int, rj: Qubit, lj: Qubit) -> None:
        addend, tgt = qpair(j)
        idx = j
        # Work1/Work2's r fields are big endian.  The low boundary R uses the
        # clean carry; cells toward L use the transformed lower addend bit as
        # the Cuccaro carry chain.
        if idx + 1 < rel_count:
            _apply_cell_dirty(
                qc, mode, "first", acc, addend, tgt, qpair(idx + 1)[0], cell_pool[0]
            )
        _apply_cell_dirty(qc, mode, "first", rj, addend, tgt, carry, cell_pool[0])
        if sign_update:
            qc.ccx(lj, addend, Sign[0])
        qc.cx(rj, acc)
        qc.cx(lj, acc)
    if top_special:
        addend, tgt = qpair(top_rel)
        _toggle_eq_const_under_ctrl_direct(qc, endpoint=l_s, const=top_rel, ctrl=Ctrl[0], acc=top_flag, scratch=eq_scratch)
        _apply_cell_dirty(qc, mode, "first", top_flag, addend, tgt, carry, cell_pool[0])
        qc.cx(top_flag, acc)
        _toggle_eq_const_under_ctrl_direct(qc, endpoint=l_s, const=top_rel, ctrl=Ctrl[0], acc=top_flag, scratch=eq_scratch)
        _toggle_eq_const_under_ctrl_direct(qc, endpoint=l_q, const=top_rel, ctrl=Ctrl[0], acc=top_flag, scratch=eq_scratch)
        if sign_update:
            qc.ccx(top_flag, addend, Sign[0])
        qc.cx(top_flag, acc)
        _toggle_eq_const_under_ctrl_direct(qc, endpoint=l_q, const=top_rel, ctrl=Ctrl[0], acc=top_flag, scratch=eq_scratch)
    dual_unary_iteration_tight(qc, index_a=l_s, index_b=l_q, labels=labels_main,
                            ctrl_a=Ctrl[0], ctrl_b=Ctrl[0], ancillas_a=anc_a,
                            ancillas_b=anc_b, leaf_fn=leaf_first, order="dec")
    def leaf_second(j: int, rj: Qubit, lj: Qubit) -> None:
        addend, tgt = qpair(j)
        idx = j
        qc.cx(lj, acc)
        qc.cx(rj, acc)
        if idx + 1 < rel_count:
            _apply_cell_dirty(
                qc, mode, "second", acc, addend, tgt, qpair(idx + 1)[0], cell_pool[0]
            )
        _apply_cell_dirty(qc, mode, "second", rj, addend, tgt, carry, cell_pool[0])
    dual_unary_iteration_tight(qc, index_a=l_s, index_b=l_q, labels=labels_main,
                            ctrl_a=Ctrl[0], ctrl_b=Ctrl[0], ancillas_a=anc_a,
                            ancillas_b=anc_b, leaf_fn=leaf_second, order="inc")
    if top_special:
        addend, tgt = qpair(top_rel)
        _toggle_eq_const_under_ctrl_direct(qc, endpoint=l_q, const=top_rel, ctrl=Ctrl[0], acc=top_flag, scratch=eq_scratch)
        qc.cx(top_flag, acc)
        _toggle_eq_const_under_ctrl_direct(qc, endpoint=l_q, const=top_rel, ctrl=Ctrl[0], acc=top_flag, scratch=eq_scratch)
        _toggle_eq_const_under_ctrl_direct(qc, endpoint=l_s, const=top_rel, ctrl=Ctrl[0], acc=top_flag, scratch=eq_scratch)
        qc.cx(top_flag, acc)
        _apply_cell_dirty(qc, mode, "second", top_flag, addend, tgt, carry, cell_pool[0])
        _toggle_eq_const_under_ctrl_direct(qc, endpoint=l_s, const=top_rel, ctrl=Ctrl[0], acc=top_flag, scratch=eq_scratch)
    _e.add_const_mod_2n(qc, l_s, k, cs[:shift_width] + [carry])
    _e.add_const_mod_2n(qc, l_q, k, cs[:len_width] + [carry])
    _e.const_minus_inplace(qc, l_s, n + 2, cs[:shift_width] + [carry])
    _e.sub_const_mod_2n(qc, l_q, 4, cs[:len_width] + [carry])
    qc.append(_e.cuccaro_sub_mod_2n_no_z_gate(len_width, name="SUB_lt_from_lq"), list(l_t) + list(l_q) + [carry])
    return _e._finalize_block(qc)


@lru_cache(maxsize=None)
def lc_prefix_addsub_unary_gate(*, k: int, K: int, len_width: int,
                                mode: Literal["add", "sub"], sign_update: bool,
                                target: Literal["work1", "work2"], name: str,
                                endpoint_offset: int = 2) -> Gate:
    if k > K:
        raise ValueError("need k <= K")
    M = K - k + 1
    depth = _e.unary_depth(M)
    base = max(depth, len_width)
    scratch_size = base + 3
    Ctrl = QuantumRegister(1, "Ctrl")
    Sign = QuantumRegister(1, "Sign")
    Work1 = QuantumRegister(M, "Work1")
    Work2 = QuantumRegister(M, "Work2")
    l_t = QuantumRegister(len_width, "l_t")
    Scratch = QuantumRegister(scratch_size, "Scratch")
    qc = _e._block_circuit(Ctrl, Sign, Work1, Work2, l_t, Scratch, name=name)
    path = list(Scratch[:depth])
    carry = Scratch[base]
    acc = Scratch[base + 1]
    cell_pool = [Scratch[base + 2]]
    cs = list(Scratch[:len_width]) + [carry]
    _e.add_const_mod_2n(qc, l_t, endpoint_offset, cs)
    def qpair(j: int) -> tuple[Qubit, Qubit]:
        idx = j - k
        if target == "work1":
            return Work2[idx], Work1[idx]
        if target == "work2":
            return Work1[idx], Work2[idx]
        raise ValueError("bad target")
    qc.cx(Ctrl[0], acc)
    def leaf_first(j: int, ej: Qubit) -> None:
        addend, tgt = qpair(j)
        if j == k:
            _apply_cell_dirty(qc, mode, "first", Ctrl[0], addend, tgt, carry, cell_pool[0])
        else:
            _apply_cell_dirty(qc, mode, "first", acc, addend, tgt, qpair(j - 1)[0], cell_pool[0])
        if sign_update:
            qc.ccx(ej, addend, Sign[0])
        qc.cx(ej, acc)
    unary_iteration_tight(qc, index_reg=l_t, labels=list(range(k, K + 1)), ctrl=Ctrl[0], ancillas=path, leaf_fn=leaf_first, order="inc")
    def leaf_second(j: int, ej: Qubit) -> None:
        addend, tgt = qpair(j)
        qc.cx(ej, acc)
        if j == k:
            _apply_cell_dirty(qc, mode, "second", Ctrl[0], addend, tgt, carry, cell_pool[0])
        else:
            _apply_cell_dirty(qc, mode, "second", acc, addend, tgt, qpair(j - 1)[0], cell_pool[0])
    unary_iteration_tight(qc, index_reg=l_t, labels=list(range(k, K + 1)), ctrl=Ctrl[0], ancillas=path, leaf_fn=leaf_second, order="dec")
    qc.cx(Ctrl[0], acc)
    _e.sub_const_mod_2n(qc, l_t, endpoint_offset, cs)
    return _e._finalize_block(qc)


def _upper_zero_map_controlled(qc: QuantumCircuit, *, ctrl: Qubit,
                               boundary_B: Sequence[Qubit], bits: Sequence[Qubit],
                               dirty: Sequence[Qubit], k: int, K: int,
                               scratch: Sequence[Qubit]) -> None:
    """Controlled upper-zero dirty map with one shared palindromic scan."""
    depth = _e.unary_depth(K - k + 1)
    if len(scratch) < depth + 2:
        raise ValueError("controlled upper-zero map scratch shortage")
    path = list(scratch[:depth])
    range_acc = scratch[depth]
    a_tmp = scratch[depth + 1]

    def compute_factor(bctrl: Qubit, bit: Qubit) -> None:
        # ctrl & !(bctrl & bit): out-of-range positions contribute the
        # multiplicative identity when active, while ctrl=0 is exact identity.
        qc.cx(ctrl, a_tmp)
        qc.ccx(bctrl, bit, a_tmp)

    def leaf_forward(j: int, bctrl: Qubit) -> None:
        idx = j - k
        if j == K:
            # At the pivot, a_K = ctrl xor ([K <= B] & bit_K).  Applying it
            # directly removes one compute/action/uncompute Toffoli.
            qc.cx(ctrl, dirty[idx])
            qc.ccx(bctrl, bits[idx], dirty[idx])
            return
        compute_factor(bctrl, bits[idx])
        qc.ccx(a_tmp, dirty[idx + 1], dirty[idx])
        compute_factor(bctrl, bits[idx])

    def leaf_reverse(j: int, bctrl: Qubit) -> None:
        idx = j - k
        compute_factor(bctrl, bits[idx])
        qc.ccx(a_tmp, dirty[idx + 1], dirty[idx])
        compute_factor(bctrl, bits[idx])

    labels = list(range(k, K + 1))

    def scan_forward(sub_labels: list[int], g: Qubit, level: int) -> None:
        if len(sub_labels) == 1:
            leaf_forward(sub_labels[0], range_acc)
            qc.cx(g, range_acc)
            return
        bit = _e._split_bit(sub_labels)
        zero = [j for j in sub_labels if ((j >> bit) & 1) == 0]
        one = [j for j in sub_labels if ((j >> bit) & 1) == 1]
        h = path[level]
        _e._and_with_index_bit(qc, g, boundary_B[bit], h, 0)
        scan_forward(zero, h, level + 1)
        qc.cx(g, h)
        scan_forward(one, h, level + 1)
        qc.cx(g, h)
        _e._uncompute_and_with_index_bit(qc, g, boundary_B[bit], h, 0)

    def scan_reverse(sub_labels: list[int], g: Qubit, level: int) -> None:
        if len(sub_labels) == 1:
            qc.cx(g, range_acc)
            leaf_reverse(sub_labels[0], range_acc)
            return
        bit = _e._split_bit(sub_labels)
        zero = [j for j in sub_labels if ((j >> bit) & 1) == 0]
        one = [j for j in sub_labels if ((j >> bit) & 1) == 1]
        h = path[level]
        _e._and_with_index_bit(qc, g, boundary_B[bit], h, 0)
        qc.cx(g, h)
        scan_reverse(one, h, level + 1)
        qc.cx(g, h)
        scan_reverse(zero, h, level + 1)
        _e._uncompute_and_with_index_bit(qc, g, boundary_B[bit], h, 0)

    def scan_palindrome(sub_labels: list[int], g: Qubit, level: int) -> None:
        if len(sub_labels) == 1:
            leaf_forward(sub_labels[0], range_acc)
            return
        bit = _e._split_bit(sub_labels)
        zero = [j for j in sub_labels if ((j >> bit) & 1) == 0]
        one = [j for j in sub_labels if ((j >> bit) & 1) == 1]
        h = path[level]
        _e._and_with_index_bit(qc, g, boundary_B[bit], h, 0)
        scan_forward(zero, h, level + 1)
        qc.cx(g, h)
        scan_palindrome(one, h, level + 1)
        qc.cx(g, h)
        scan_reverse(zero, h, level + 1)
        _e._uncompute_and_with_index_bit(qc, g, boundary_B[bit], h, 0)

    qc.cx(ctrl, range_acc)
    scan_palindrome(labels, ctrl, 0)
    qc.cx(ctrl, range_acc)


@lru_cache(maxsize=None)
def t_tail_zero_toggle_gate(*, n: int, len_width: int, shift_width: int,
                            name: str = "T_TAIL_ZERO_S835_FAST") -> Gate:
    """Toggle Tail iff Work2[A..=B] is zero for the dynamic t' tail."""
    work_size = n + 3
    labels = list(range(work_size))
    depth = _tight_unary_depth_for_labels(labels)
    map_need = _e.unary_depth(work_size) + 2

    def pivot_depth(sub_labels: list[int], pivot: int) -> int:
        if len(sub_labels) <= 1:
            return 0
        bit = _e._split_bit(sub_labels)
        branch = [j for j in sub_labels if ((j >> bit) & 1) == ((pivot >> bit) & 1)]
        return 1 + pivot_depth(branch, pivot)

    live_select_depth = pivot_depth(labels, labels[-1])

    Ctrl = QuantumRegister(1, "Ctrl")
    Tail = QuantumRegister(1, "Tail")
    Work1 = QuantumRegister(work_size, "Work1")
    Work2 = QuantumRegister(work_size, "Work2")
    l_t = QuantumRegister(len_width, "l_t")
    l_s = QuantumRegister(shift_width, "l_s")
    l_rp = QuantumRegister(len_width, "l_rp")
    map_offset = 0
    select_offset = map_need
    carry_offset = select_offset + live_select_depth
    Scratch = QuantumRegister(carry_offset + 1, "Scratch")
    qc = _e._block_circuit(Ctrl, Tail, Work1, Work2, l_t, l_s, l_rp, Scratch, name=name)
    length_carry = Scratch[carry_offset]

    def shift_lower_endpoint(forward: bool) -> None:
        # Adding two modulo 2^w is an increment of bits 1..w-1.
        if len_width <= 1:
            return
        upper = list(l_t[1:])
        ancillas = list(Scratch[:max(0, len(upper) - 1)])
        if forward:
            _e.inc_mod2n_uncontrolled(qc, upper, ancillas)
        else:
            _e.dec_mod2n_uncontrolled(qc, upper, ancillas)

    def reflect_upper_endpoint() -> None:
        # l_rp <- n-l_rp.  At n=256 the constant is the top bit of the
        # 9-bit endpoint, so its modular addition is a single X.
        for q in l_rp:
            qc.x(q)
        _e.inc_mod2n_uncontrolled(qc, l_rp, list(Scratch[:max(0, len_width - 1)]))
        if n == (1 << (len_width - 1)):
            qc.x(l_rp[len_width - 1])
        else:
            _e.add_const_mod_2n(
                qc, l_rp, n, list(Scratch[:len_width]) + [length_carry]
            )

    def transform_endpoints() -> None:
        # A=l_t+1 (after the appended zero lane) and
        # B=n+2-l_r'-l_s in zero-based physical coordinates.
        shift_lower_endpoint(True)
        qc.append(
            _e.cuccaro_add_mod_2n_no_z_gate(len_width, name="ADD_ls_to_lrp"),
            list(l_s[:len_width]) + list(l_rp) + [length_carry],
        )
        reflect_upper_endpoint()

    def restore_endpoints() -> None:
        reflect_upper_endpoint()
        qc.append(
            _e.cuccaro_sub_mod_2n_no_z_gate(len_width, name="SUB_ls_from_lrp"),
            list(l_s[:len_width]) + list(l_rp) + [length_carry],
        )
        shift_lower_endpoint(False)

    map_scratch = list(Scratch[map_offset:map_offset + map_need])
    # Only the path to the maximum label remains live across the central map.
    # Give those levels dedicated wires; all deeper selector levels are clean
    # before the map and can alias its scratch without widening the EEA step.
    select_path = (
        list(Scratch[select_offset:select_offset + live_select_depth])
        + map_scratch[:depth - live_select_depth]
    )

    def apply_upper_map() -> None:
        _upper_zero_map_controlled(
            qc, ctrl=Ctrl[0], boundary_B=l_rp, bits=Work2, dirty=Work1,
            k=0, K=work_size - 1, scratch=map_scratch,
        )

    def selected_leaf(j: int, ej: Qubit) -> None:
        qc.ccx(ej, Work1[j], Tail[0])

    def select_forward(sub_labels: list[int], g: Qubit, level: int) -> None:
        if len(sub_labels) == 1:
            selected_leaf(sub_labels[0], g)
            return
        bit = _e._split_bit(sub_labels)
        zero = [j for j in sub_labels if ((j >> bit) & 1) == 0]
        one = [j for j in sub_labels if ((j >> bit) & 1) == 1]
        h = select_path[level]
        _e._and_with_index_bit(qc, g, l_t[bit], h, 0)
        select_forward(zero, h, level + 1)
        qc.cx(g, h)
        select_forward(one, h, level + 1)
        qc.cx(g, h)
        _e._uncompute_and_with_index_bit(qc, g, l_t[bit], h, 0)

    def select_reverse(sub_labels: list[int], g: Qubit, level: int) -> None:
        if len(sub_labels) == 1:
            selected_leaf(sub_labels[0], g)
            return
        bit = _e._split_bit(sub_labels)
        zero = [j for j in sub_labels if ((j >> bit) & 1) == 0]
        one = [j for j in sub_labels if ((j >> bit) & 1) == 1]
        h = select_path[level]
        _e._and_with_index_bit(qc, g, l_t[bit], h, 0)
        qc.cx(g, h)
        select_reverse(one, h, level + 1)
        qc.cx(g, h)
        select_reverse(zero, h, level + 1)
        _e._uncompute_and_with_index_bit(qc, g, l_t[bit], h, 0)

    def select_map_palindrome(sub_labels: list[int], g: Qubit, level: int) -> None:
        if len(sub_labels) == 1:
            selected_leaf(sub_labels[0], g)
            apply_upper_map()
            selected_leaf(sub_labels[0], g)
            return
        bit = _e._split_bit(sub_labels)
        zero = [j for j in sub_labels if ((j >> bit) & 1) == 0]
        one = [j for j in sub_labels if ((j >> bit) & 1) == 1]
        h = select_path[level]
        _e._and_with_index_bit(qc, g, l_t[bit], h, 0)
        select_forward(zero, h, level + 1)
        qc.cx(g, h)
        select_map_palindrome(one, h, level + 1)
        qc.cx(g, h)
        select_reverse(zero, h, level + 1)
        _e._uncompute_and_with_index_bit(qc, g, l_t[bit], h, 0)

    transform_endpoints()
    select_map_palindrome(labels, Ctrl[0], 0)
    apply_upper_map()
    restore_endpoints()
    return _e._finalize_block(qc)


@lru_cache(maxsize=None)
def t_lower_borrow_toggle_gate(*, n: int, len_width: int,
                               name: str = "T_LOWER_BORROW_S835_FAST") -> Gate:
    """Toggle Neg by Tail times the exact borrow through the t prefix."""
    work_size = n + 3
    labels = list(range(1, work_size + 1))
    depth = _tight_unary_depth_for_labels(labels)
    base = max(depth, len_width)
    Ctrl = QuantumRegister(1, "Ctrl")
    Tail = QuantumRegister(1, "Tail")
    Neg = QuantumRegister(1, "Neg")
    Work1 = QuantumRegister(work_size, "Work1")
    Work2 = QuantumRegister(work_size, "Work2")
    l_t = QuantumRegister(len_width, "l_t")
    Scratch = QuantumRegister(base + 2, "Scratch")
    qc = _e._block_circuit(Ctrl, Tail, Neg, Work1, Work2, l_t, Scratch, name=name)
    carry = Scratch[base]
    active = Scratch[base + 1]

    # The first inverse-UMA pass of the controlled prefix subtractor stores
    # the borrow through position j in Work1[j].  Execute that pass without a
    # location control, use its intermediate value at the selected endpoint,
    # then reverse it.  The surrounding permutation cancels even when the
    # output control is inactive, so only the unary selector needs Ctrl&Tail.
    if len_width > 1:
        _e.inc_mod2n_uncontrolled(
            qc, l_t[1:], list(Scratch[:max(0, len_width - 2)])
        )
    qc.ccx(Ctrl[0], Tail[0], active)

    def first_pass_cell(idx: int) -> None:
        addend = Work1[idx]
        target = Work2[idx]
        carry_in = carry if idx == 0 else Work1[idx - 1]
        qc.cx(carry_in, target)
        qc.cx(addend, carry_in)
        qc.ccx(carry_in, target, addend)

    def leaf(j: int, ej: Qubit) -> None:
        idx = j - 1
        first_pass_cell(idx)
        qc.ccx(ej, Work1[idx], Neg[0])

    unary_iteration_tight(
        qc, index_reg=l_t, labels=labels, ctrl=active,
        ancillas=list(Scratch[:depth]), leaf_fn=leaf, order="inc",
    )

    for idx in range(work_size - 1, -1, -1):
        addend = Work1[idx]
        target = Work2[idx]
        carry_in = carry if idx == 0 else Work1[idx - 1]
        qc.ccx(carry_in, target, addend)
        qc.cx(addend, carry_in)
        qc.cx(carry_in, target)

    qc.ccx(Ctrl[0], Tail[0], active)
    if len_width > 1:
        _e.dec_mod2n_uncontrolled(
            qc, l_t[1:], list(Scratch[:max(0, len_width - 2)])
        )
    return _e._finalize_block(qc)

# Reuse the low-aux length update; it is already the paper dirty-work construction with live-range shared scratch.
import eea_circuit_s835_lowaux as _low
len_update_lt_unary_gate = _low.len_update_lt_unary_gate
len_update_lrp_unary_gate = _low.len_update_lrp_unary_gate


def _borrowed_c3x(qc: QuantumCircuit, a: Qubit, b: Qubit, c: Qubit,
                  target: Qubit, borrowed: Qubit) -> None:
    """Exact C3X using one unknown borrowed bit, restored with no phase."""
    qc.ccx(a, b, borrowed)
    qc.ccx(borrowed, c, target)
    qc.ccx(a, b, borrowed)
    qc.ccx(borrowed, c, target)


def _mcx_dirty_ladder(qc: QuantumCircuit, controls: Sequence[Qubit],
                      target: Qubit, dirty: Sequence[Qubit]) -> None:
    """Toggle ``target`` by all controls, restoring unknown dirty lenders.

    This is the exact ``4*k - 8``-CCX construction used by the Rust KMX
    lowerer in ``arith/mcx.rs``.  The first cascade includes the seed link;
    the second omits it, cancelling every dirty-seeded term while retaining
    the complete control product once.
    """
    k = len(controls)
    if k == 0:
        qc.x(target)
        return
    if k == 1:
        qc.cx(controls[0], target)
        return
    if k == 2:
        qc.ccx(controls[0], controls[1], target)
        return
    if len(dirty) < k - 2:
        raise ValueError(f"dirty MCX needs {k - 2} lenders, got {len(dirty)}")
    lenders = list(dirty[:k - 2])
    lanes = list(controls) + [target] + lenders
    if len({id(lane) for lane in lanes}) != len(lanes):
        raise ValueError("dirty MCX lanes must be distinct")

    def cascade(include_seed: bool) -> None:
        if include_seed:
            qc.ccx(controls[0], controls[1], lenders[0])
        for index in range(1, len(lenders)):
            qc.ccx(lenders[index - 1], controls[index + 1], lenders[index])
        qc.ccx(lenders[-1], controls[k - 1], target)
        for index in range(len(lenders) - 1, 0, -1):
            qc.ccx(lenders[index - 1], controls[index + 1], lenders[index])
        if include_seed:
            qc.ccx(controls[0], controls[1], lenders[0])

    cascade(True)
    cascade(False)


def _apply_cell_borrowed(qc: QuantumCircuit, mode: Literal["add", "sub"],
                         pass_kind: Literal["first", "second"], ctrl: Qubit,
                         addend: Qubit, target: Qubit, carry: Qubit,
                         borrowed: Qubit) -> None:
    def cmaj() -> None:
        qc.ccx(ctrl, addend, target)
        qc.ccx(ctrl, addend, carry)
        _borrowed_c3x(qc, ctrl, carry, target, addend, borrowed)

    def cuma() -> None:
        _borrowed_c3x(qc, ctrl, carry, target, addend, borrowed)
        qc.ccx(ctrl, addend, carry)
        qc.ccx(ctrl, carry, target)

    def cmaj_inv() -> None:
        _borrowed_c3x(qc, ctrl, carry, target, addend, borrowed)
        qc.ccx(ctrl, addend, carry)
        qc.ccx(ctrl, addend, target)

    def cuma_inv() -> None:
        qc.ccx(ctrl, carry, target)
        qc.ccx(ctrl, addend, carry)
        _borrowed_c3x(qc, ctrl, carry, target, addend, borrowed)

    table = {
        ("add", "first"): cmaj,
        ("add", "second"): cuma,
        ("sub", "first"): cuma_inv,
        ("sub", "second"): cmaj_inv,
    }
    try:
        table[(mode, pass_kind)]()
    except KeyError as exc:
        raise ValueError("bad borrowed arithmetic cell mode/pass") from exc


@lru_cache(maxsize=None)
def compact_lc_swap_gate(*, k: int, K: int,
                         name: str = "LC_SWAP_COMPACT") -> Gate:
    if k > K:
        raise ValueError("need k <= K")
    M = K - k + 1
    Ctrl = QuantumRegister(1, "Ctrl")
    Direction = QuantumRegister(1, "Direction")
    Sign = QuantumRegister(1, "Sign")
    Work1 = QuantumRegister(M + 1, "Work1")
    l_t = QuantumRegister(LT_WIDTH, "l_t")
    l_q = QuantumRegister(LQ_WIDTH, "l_q")
    depth = _tight_unary_depth_for_labels(list(range(k, K + 1)))
    base = max(LQ_WIDTH, depth)
    Scratch = QuantumRegister(base + 2, "Scratch")
    qc = _e._block_circuit(Ctrl, Direction, Sign, Work1, l_t, l_q, Scratch, name=name)
    path = list(Scratch[:depth])
    extension = Scratch[LQ_WIDTH - 1]
    carry = Scratch[base]
    direction_flag = Scratch[base + 1]
    qc.append(_e.cuccaro_add_mod_2n_no_z_gate(LQ_WIDTH, name="ADD_lt8_to_lq9"),
              list(l_t) + [extension] + list(l_q) + [carry])
    _e.add_const_mod_2n(qc, l_q, 3, list(Scratch[:LQ_WIDTH]) + [carry])

    def leaf(j: int, ej: Qubit) -> None:
        _e._and_with_index_bit(qc, ej, Direction[0], direction_flag, 0)
        _e.cswap_toffoli(qc, direction_flag, Sign[0], Work1[j - k + 1])
        qc.cx(ej, direction_flag)
        _e.cswap_toffoli(qc, direction_flag, Sign[0], Work1[j - k])
        qc.cx(ej, direction_flag)
        _e._uncompute_and_with_index_bit(qc, ej, Direction[0], direction_flag, 0)

    unary_iteration_tight(
        qc, index_reg=l_q, labels=list(range(k, K + 1)), ctrl=Ctrl[0],
        ancillas=path, leaf_fn=leaf, order="inc",
    )
    _e.sub_const_mod_2n(qc, l_q, 3, list(Scratch[:LQ_WIDTH]) + [carry])
    qc.append(_e.cuccaro_sub_mod_2n_no_z_gate(LQ_WIDTH, name="SUB_lt8_from_lq9"),
              list(l_t) + [extension] + list(l_q) + [carry])
    return _e._finalize_block(qc)


@lru_cache(maxsize=None)
def compact_interval_addsub_gate(*, n: int, k: int, K: int,
                                 mode: Literal["add", "sub"], sign_update: bool,
                                 target: Literal["work1", "work2"], name: str) -> Gate:
    if k > K:
        raise ValueError("need k <= K")
    M = K - k + 1
    Ctrl = QuantumRegister(1, "Ctrl")
    Sign = QuantumRegister(1, "Sign")
    Work1 = QuantumRegister(M, "Work1")
    Work2 = QuantumRegister(M, "Work2")
    l_t = QuantumRegister(LT_WIDTH, "l_t")
    l_q = QuantumRegister(LQ_WIDTH, "l_q")
    l_s = QuantumRegister(LS_WIDTH, "l_s")
    Dirty = QuantumRegister(DIRTY_PASSENGER_SIZE, "DirtyPassenger")
    Scratch = QuantumRegister(11, "Scratch")
    qc = _e._block_circuit(Ctrl, Sign, Work1, Work2, l_t, l_q, l_s,
                           Dirty, Scratch, name=name)
    kg_s = list(Scratch[0:3])
    kg_q = list(Scratch[3:6])
    eq_s = Scratch[6]
    eq_q = Scratch[7]
    carry = Scratch[8]
    acc = Scratch[9]
    extension = Scratch[10]
    cell_borrowed = Dirty[9]
    qc.append(_e.cuccaro_add_mod_2n_no_z_gate(LQ_WIDTH, name="ADD_lt8_to_lq9"),
              list(l_t) + [extension] + list(l_q) + [carry])
    affine_scratch = list(Scratch[:8]) + [extension, carry]
    _e.add_const_mod_2n(qc, l_q, 4, affine_scratch)
    _e.const_minus_inplace(qc, l_s, n + 2, affine_scratch)
    # In the modulo-259 encoding ell_s=0 is stored as integer 258.  The
    # affine endpoint reflection first maps that word to 0, whereas the
    # Aux22/v2 signed-sentinel endpoint is physical label 259.  This basis
    # transposition repairs exactly that case and is its own inverse.
    _swap_zero_259_uncontrolled(qc, l_s, extension, list(Scratch[:9]))

    def qpair(j: int) -> tuple[Qubit, Qubit]:
        idx = j - k
        if target == "work1":
            return Work2[idx], Work1[idx]
        if target == "work2":
            return Work1[idx], Work2[idx]
        raise ValueError("bad compact interval target")

    def leaf_first(j: int, sj: Qubit, qj: Qubit) -> None:
        addend, tgt = qpair(j)
        if j < K:
            next_addend, _ = qpair(j + 1)
            _apply_cell_borrowed(
                qc, mode, "first", acc, addend, tgt,
                next_addend, cell_borrowed,
            )
        _apply_cell_borrowed(
            qc, mode, "first", sj, addend, tgt, carry, cell_borrowed,
        )
        qc.cx(sj, acc)
        qc.cx(qj, acc)
        if sign_update:
            qc.ccx(qj, addend, Sign[0])

    dual_unary_iteration_log_star(
        qc, index_a=l_s, index_b=l_q, labels=list(range(k, K + 1)),
        ancillas_a=kg_s, ancillas_b=kg_q, flag_a=eq_s, flag_b=eq_q,
        common_ctrl=Ctrl[0], clean_temp=extension,
        leaf_fn=leaf_first, order="dec",
    )

    def leaf_second(j: int, sj: Qubit, qj: Qubit) -> None:
        addend, tgt = qpair(j)
        qc.cx(qj, acc)
        qc.cx(sj, acc)
        if j < K:
            next_addend, _ = qpair(j + 1)
            _apply_cell_borrowed(
                qc, mode, "second", acc, addend, tgt,
                next_addend, cell_borrowed,
            )
        _apply_cell_borrowed(
            qc, mode, "second", sj, addend, tgt, carry, cell_borrowed,
        )

    dual_unary_iteration_log_star(
        qc, index_a=l_s, index_b=l_q, labels=list(range(k, K + 1)),
        ancillas_a=kg_s, ancillas_b=kg_q, flag_a=eq_s, flag_b=eq_q,
        common_ctrl=Ctrl[0], clean_temp=extension,
        leaf_fn=leaf_second, order="inc",
    )
    _swap_zero_259_uncontrolled(qc, l_s, extension, list(Scratch[:9]))
    _e.const_minus_inplace(qc, l_s, n + 2, affine_scratch)
    _e.sub_const_mod_2n(qc, l_q, 4, affine_scratch)
    qc.append(_e.cuccaro_sub_mod_2n_no_z_gate(LQ_WIDTH, name="SUB_lt8_from_lq9"),
              list(l_t) + [extension] + list(l_q) + [carry])
    return _e._finalize_block(qc)


@lru_cache(maxsize=None)
def compact_prefix_addsub_gate(*, k: int, K: int,
                               mode: Literal["add", "sub"], sign_update: bool,
                               capture_borrow_sign: bool,
                               target: Literal["work1", "work2"], name: str) -> Gate:
    if k > K:
        raise ValueError("need k <= K")
    if k != 1 or K > 257:
        raise ValueError("compact T prefix is certified for physical labels 1..257")
    if sign_update:
        raise ValueError("compact T prefix sign update must use selected midpoint capture")
    M = K - k + 1
    Ctrl = QuantumRegister(1, "Ctrl")
    Sign = QuantumRegister(1, "Sign")
    Tail = QuantumRegister(1, "Tail")
    Work1 = QuantumRegister(M, "Work1")
    Work2 = QuantumRegister(M, "Work2")
    l_t = QuantumRegister(LT_WIDTH, "l_t")
    Borrowed = QuantumRegister(1, "Borrowed")
    # l_t is stored as truth-minus-one.  Keep it unmodified and decode
    # residues x=0..K-2 as physical cells j=x+2.  Physical cell 1 is the
    # unconditional lower boundary and is emitted explicitly.
    encoded_labels = list(range(0, K - 1))
    depth = _tight_unary_depth_for_labels(encoded_labels)
    base = max(depth, LT_WIDTH)
    Scratch = QuantumRegister(base + 2, "Scratch")
    qc = _e._block_circuit(Ctrl, Sign, Tail, Work1, Work2, l_t,
                           Borrowed, Scratch, name=name)
    path = list(Scratch[:depth])
    carry = Scratch[base]
    acc = Scratch[base + 1]

    def qpair(j: int) -> tuple[Qubit, Qubit]:
        idx = j - k
        if target == "work1":
            return Work2[idx], Work1[idx]
        if target == "work2":
            return Work1[idx], Work2[idx]
        raise ValueError("bad compact prefix target")

    def leaf_first(encoded: int, ej: Qubit) -> None:
        j = encoded + 2
        addend, tgt = qpair(j)
        previous_addend, _ = qpair(j - 1)
        _apply_cell_borrowed(
            qc, mode, "first", acc, addend, tgt,
            previous_addend, Borrowed[0],
        )
        if capture_borrow_sign:
            # After the first cell, addend stores the exact borrow through the
            # selected physical endpoint.  ej already contains Ctrl and the
            # endpoint equality, so Tail & ej & addend is the old retained
            # Neg predicate without a separate history bit or rescan.
            _borrowed_c3x(
                qc, Tail[0], ej, addend, Sign[0], Borrowed[0],
            )
        qc.cx(ej, acc)

    qc.cx(Ctrl[0], acc)
    addend1, tgt1 = qpair(1)
    # Scratch[0] is clean outside the unary tree, so the boundary cell uses
    # the clean MBU C3X lowering and returns it to zero before the tree starts.
    _apply_cell_dirty(
        qc, mode, "first", Ctrl[0], addend1, tgt1, carry, Scratch[0],
    )
    if encoded_labels:
        unary_iteration_tight(
            qc, index_reg=l_t, labels=encoded_labels, ctrl=Ctrl[0],
            ancillas=path, leaf_fn=leaf_first, order="inc",
        )

    def leaf_second(encoded: int, ej: Qubit) -> None:
        j = encoded + 2
        addend, tgt = qpair(j)
        qc.cx(ej, acc)
        previous_addend, _ = qpair(j - 1)
        _apply_cell_borrowed(
            qc, mode, "second", acc, addend, tgt,
            previous_addend, Borrowed[0],
        )

    if encoded_labels:
        unary_iteration_tight(
            qc, index_reg=l_t, labels=encoded_labels, ctrl=Ctrl[0],
            ancillas=path, leaf_fn=leaf_second, order="dec",
        )
    _apply_cell_dirty(
        qc, mode, "second", Ctrl[0], addend1, tgt1, carry, Scratch[0],
    )
    qc.cx(Ctrl[0], acc)
    return _e._finalize_block(qc)


def _apply_not_factor_with_borrowed(qc: QuantumCircuit, *, boundary_control: Qubit,
                                    data_bit: Qubit, neighbor: Optional[Qubit],
                                    target: Qubit, borrowed: Qubit) -> None:
    """Apply X or neighbor-controlled X under NOT(boundary_control & data_bit)."""
    if neighbor is None:
        qc.x(target)
        qc.cx(borrowed, target)
        qc.ccx(boundary_control, data_bit, borrowed)
        qc.cx(borrowed, target)
        qc.ccx(boundary_control, data_bit, borrowed)
    else:
        qc.cx(neighbor, target)
        qc.ccx(borrowed, neighbor, target)
        qc.ccx(boundary_control, data_bit, borrowed)
        qc.ccx(borrowed, neighbor, target)
        qc.ccx(boundary_control, data_bit, borrowed)


def _range_scan_tight(qc: QuantumCircuit, *, leq: bool,
                      boundary: Sequence[Qubit], k: int, K: int,
                      ctrl: Qubit, range_acc: Qubit,
                      path: Sequence[Qubit], leaf_fn,
                      order: Literal["inc", "dec"]) -> None:
    labels = list(range(k, K + 1))
    if leq and order == "inc":
        qc.cx(ctrl, range_acc)
        def wrapped(j: int, ej: Qubit) -> None:
            leaf_fn(j, range_acc)
            qc.cx(ej, range_acc)
        unary_iteration_tight(qc, index_reg=boundary, labels=labels, ctrl=ctrl,
                              ancillas=path, leaf_fn=wrapped, order=order)
    elif leq and order == "dec":
        def wrapped(j: int, ej: Qubit) -> None:
            qc.cx(ej, range_acc)
            leaf_fn(j, range_acc)
        unary_iteration_tight(qc, index_reg=boundary, labels=labels, ctrl=ctrl,
                              ancillas=path, leaf_fn=wrapped, order=order)
        qc.cx(ctrl, range_acc)
    elif not leq and order == "inc":
        def wrapped(j: int, ej: Qubit) -> None:
            qc.cx(ej, range_acc)
            leaf_fn(j, range_acc)
        unary_iteration_tight(qc, index_reg=boundary, labels=labels, ctrl=ctrl,
                              ancillas=path, leaf_fn=wrapped, order=order)
        qc.cx(ctrl, range_acc)
    elif not leq and order == "dec":
        qc.cx(ctrl, range_acc)
        def wrapped(j: int, ej: Qubit) -> None:
            leaf_fn(j, range_acc)
            qc.cx(ej, range_acc)
        unary_iteration_tight(qc, index_reg=boundary, labels=labels, ctrl=ctrl,
                              ancillas=path, leaf_fn=wrapped, order=order)
    else:
        raise ValueError("bad tight range-scan order")


def _upper_zero_map_borrowed(qc: QuantumCircuit, *, ctrl: Qubit,
                             boundary_B: Sequence[Qubit], bits: Sequence[Qubit],
                             dirty_map: Sequence[Qubit], borrowed: Qubit,
                             k: int, K: int, scratch: Sequence[Qubit]) -> None:
    depth = _tight_unary_depth_for_labels(list(range(k, K + 1)))
    if len(scratch) < depth + 1:
        raise ValueError("borrowed upper-zero map scratch shortage")
    path = list(scratch[:depth])
    range_acc = scratch[depth]

    def leaf_forward(j: int, bctrl: Qubit) -> None:
        idx = j - k
        _apply_not_factor_with_borrowed(
            qc, boundary_control=bctrl, data_bit=bits[idx],
            neighbor=None if j == K else dirty_map[idx + 1],
            target=dirty_map[idx], borrowed=borrowed,
        )

    def leaf_reverse(j: int, bctrl: Qubit) -> None:
        if j < K:
            idx = j - k
            _apply_not_factor_with_borrowed(
                qc, boundary_control=bctrl, data_bit=bits[idx],
                neighbor=dirty_map[idx + 1], target=dirty_map[idx],
                borrowed=borrowed,
            )

    _range_scan_tight(qc, leq=True, boundary=boundary_B, k=k, K=K, ctrl=ctrl,
                      range_acc=range_acc, path=path, leaf_fn=leaf_forward, order="inc")
    _range_scan_tight(qc, leq=True, boundary=boundary_B, k=k, K=K, ctrl=ctrl,
                      range_acc=range_acc, path=path, leaf_fn=leaf_reverse, order="dec")


def _lower_zero_map_borrowed(qc: QuantumCircuit, *, ctrl: Qubit,
                             boundary_A: Sequence[Qubit], bits: Sequence[Qubit],
                             dirty_map: Sequence[Qubit], borrowed: Qubit,
                             k: int, K: int, scratch: Sequence[Qubit]) -> None:
    depth = _tight_unary_depth_for_labels(list(range(k, K + 1)))
    if len(scratch) < depth + 1:
        raise ValueError("borrowed lower-zero map scratch shortage")
    path = list(scratch[:depth])
    range_acc = scratch[depth]

    def leaf_forward(j: int, bctrl: Qubit) -> None:
        idx = j - k
        _apply_not_factor_with_borrowed(
            qc, boundary_control=bctrl, data_bit=bits[idx],
            neighbor=None if j == k else dirty_map[idx - 1],
            target=dirty_map[idx], borrowed=borrowed,
        )

    def leaf_reverse(j: int, bctrl: Qubit) -> None:
        if j > k:
            idx = j - k
            _apply_not_factor_with_borrowed(
                qc, boundary_control=bctrl, data_bit=bits[idx],
                neighbor=dirty_map[idx - 1], target=dirty_map[idx],
                borrowed=borrowed,
            )

    _range_scan_tight(qc, leq=False, boundary=boundary_A, k=k, K=K, ctrl=ctrl,
                      range_acc=range_acc, path=path, leaf_fn=leaf_forward, order="dec")
    _range_scan_tight(qc, leq=False, boundary=boundary_A, k=k, K=K, ctrl=ctrl,
                      range_acc=range_acc, path=path, leaf_fn=leaf_reverse, order="inc")


def _highest_position_xor_write_borrowed(qc: QuantumCircuit, *, ctrl: Qubit,
                                         boundary_B: Sequence[Qubit], bits: Sequence[Qubit],
                                         dirty_map: Sequence[Qubit], target_len: Sequence[Qubit],
                                         borrowed: Qubit, k: int, K: int,
                                         scratch: Sequence[Qubit]) -> None:
    mask = (1 << len(target_len)) - 1

    def writes() -> None:
        for j in range(K, k, -1):
            _e.xor_const_into_reg_controls(
                qc, target_len, ((j - 1) ^ (j - 2)) & mask,
                ctrls=[ctrl, dirty_map[j - k]], scratch=scratch,
            )
        _e.xor_const_into_reg_controls(
            qc, target_len, ((k - 1) ^ mask) & mask,
            ctrls=[ctrl, dirty_map[0]], scratch=scratch,
        )

    _e.xor_const_into_reg_controls(qc, target_len, (K - 1) & mask,
                                   ctrls=[ctrl], scratch=scratch)
    writes()
    _upper_zero_map_borrowed(
        qc, ctrl=ctrl, boundary_B=boundary_B, bits=bits, dirty_map=dirty_map,
        borrowed=borrowed, k=k, K=K, scratch=scratch,
    )
    writes()
    _upper_zero_map_borrowed(
        qc, ctrl=ctrl, boundary_B=boundary_B, bits=bits, dirty_map=dirty_map,
        borrowed=borrowed, k=k, K=K, scratch=scratch,
    )


def _right_length_xor_write_borrowed(qc: QuantumCircuit, *, n: int, ctrl: Qubit,
                                     boundary_A: Sequence[Qubit], bits: Sequence[Qubit],
                                     dirty_map: Sequence[Qubit], target_len: Sequence[Qubit],
                                     borrowed: Qubit, k: int, K: int,
                                     scratch: Sequence[Qubit]) -> None:
    mask = (1 << len(target_len)) - 1

    def val(pos: int) -> int:
        return (n + 3 - pos) & mask

    def writes() -> None:
        for j in range(k, K):
            _e.xor_const_into_reg_controls(
                qc, target_len, val(j) ^ val(j + 1),
                ctrls=[ctrl, dirty_map[j - k]], scratch=scratch,
            )
        _e.xor_const_into_reg_controls(
            qc, target_len, val(K) ^ mask,
            ctrls=[ctrl, dirty_map[K - k]], scratch=scratch,
        )

    _e.xor_const_into_reg_controls(qc, target_len, val(k),
                                   ctrls=[ctrl], scratch=scratch)
    writes()
    _lower_zero_map_borrowed(
        qc, ctrl=ctrl, boundary_A=boundary_A, bits=bits, dirty_map=dirty_map,
        borrowed=borrowed, k=k, K=K, scratch=scratch,
    )
    writes()
    _lower_zero_map_borrowed(
        qc, ctrl=ctrl, boundary_A=boundary_A, bits=bits, dirty_map=dirty_map,
        borrowed=borrowed, k=k, K=K, scratch=scratch,
    )


@lru_cache(maxsize=None)
def compact_len_update_lt_gate(*, n: int, k: int, K: int,
                               name: str = "LEN_LT_COMPACT") -> Gate:
    M = K - k + 1
    Ctrl = QuantumRegister(1, "Ctrl")
    Work1 = QuantumRegister(M, "Work1")
    Work2 = QuantumRegister(M, "Work2")
    l_t = QuantumRegister(LT_WIDTH, "l_t")
    l_rp = QuantumRegister(LRP_WIDTH, "l_rp")
    Borrowed = QuantumRegister(1, "Borrowed")
    Scratch = QuantumRegister(11, "Scratch")
    qc = _e._block_circuit(Ctrl, Work1, Work2, l_t, l_rp, Borrowed, Scratch, name=name)
    extension = Scratch[10]
    boundary = list(l_rp) + [extension]
    map_scratch = list(Scratch[:10])
    _e.const_minus_inplace(qc, boundary, n + 2, map_scratch)
    _highest_position_xor_write_borrowed(
        qc, ctrl=Ctrl[0], boundary_B=boundary, bits=Work2, dirty_map=Work1,
        target_len=l_t, borrowed=Borrowed[0], k=k, K=K, scratch=map_scratch,
    )
    _highest_position_xor_write_borrowed(
        qc, ctrl=Ctrl[0], boundary_B=boundary, bits=Work1, dirty_map=Work2,
        target_len=l_t, borrowed=Borrowed[0], k=k, K=K, scratch=map_scratch,
    )
    _e.const_minus_inplace(qc, boundary, n + 2, map_scratch)
    return _e._finalize_block(qc)


@lru_cache(maxsize=None)
def compact_len_update_lrp_gate(*, n: int, k: int, K: int,
                                name: str = "LEN_LRP_COMPACT") -> Gate:
    M = K - k + 1
    Ctrl = QuantumRegister(1, "Ctrl")
    Work1 = QuantumRegister(M, "Work1")
    Work2 = QuantumRegister(M, "Work2")
    l_t = QuantumRegister(LT_WIDTH, "l_t")
    l_rp = QuantumRegister(LRP_WIDTH, "l_rp")
    Borrowed = QuantumRegister(1, "Borrowed")
    Scratch = QuantumRegister(11, "Scratch")
    qc = _e._block_circuit(Ctrl, Work1, Work2, l_t, l_rp, Borrowed, Scratch, name=name)
    extension = Scratch[10]
    boundary = list(l_t) + [extension]
    map_scratch = list(Scratch[:10])
    _e.add_const_mod_2n(qc, boundary, 3, map_scratch)
    _right_length_xor_write_borrowed(
        qc, n=n, ctrl=Ctrl[0], boundary_A=boundary, bits=Work1, dirty_map=Work2,
        target_len=l_rp, borrowed=Borrowed[0], k=k, K=K, scratch=map_scratch,
    )
    _right_length_xor_write_borrowed(
        qc, n=n, ctrl=Ctrl[0], boundary_A=boundary, bits=Work2, dirty_map=Work1,
        target_len=l_rp, borrowed=Borrowed[0], k=k, K=K, scratch=map_scratch,
    )
    _e.sub_const_mod_2n(qc, boundary, 3, map_scratch)
    return _e._finalize_block(qc)


@lru_cache(maxsize=None)
def compact_swap_work_and_len_gate(*, n: int, k4: int, K4: int,
                                   k5: int, K5: int,
                                   name: str = "SWAP_AND_LEN_COMPACT") -> Gate:
    work_size = n + 3
    Ctrl = QuantumRegister(1, "Ctrl")
    Work1 = QuantumRegister(work_size, "Work1")
    Work2 = QuantumRegister(work_size, "Work2")
    l_t = QuantumRegister(LT_WIDTH, "l_t")
    l_rp = QuantumRegister(LRP_WIDTH, "l_rp")
    Borrowed = QuantumRegister(1, "Borrowed")
    Scratch = QuantumRegister(11, "Scratch")
    qc = _e._block_circuit(Ctrl, Work1, Work2, l_t, l_rp, Borrowed, Scratch, name=name)
    for i in range(work_size):
        _e.cswap_toffoli(qc, Ctrl[0], Work1[i], Work2[i])
    gate_lt = compact_len_update_lt_gate(n=n, k=k4, K=K4)
    _e._append_with_optional_clbits(
        qc, gate_lt,
        [Ctrl[0]] + list(Work1[k4 - 1:K4]) + list(Work2[k4 - 1:K4])
        + list(l_t) + list(l_rp) + [Borrowed[0]] + list(Scratch),
    )
    gate_lrp = compact_len_update_lrp_gate(n=n, k=k5, K=K5)
    _e._append_with_optional_clbits(
        qc, gate_lrp,
        [Ctrl[0]] + list(Work1[k5 - 1:K5]) + list(Work2[k5 - 1:K5])
        + list(l_t) + list(l_rp) + [Borrowed[0]] + list(Scratch),
    )
    return _e._finalize_block(qc)


@lru_cache(maxsize=None)
def compact_tail_zero_gate(*, n: int,
                           name: str = "T_TAIL_ZERO_COMPACT") -> Gate:
    work_size = n + 3
    Ctrl = QuantumRegister(1, "Ctrl")
    Tail = QuantumRegister(1, "Tail")
    Work1 = QuantumRegister(work_size, "Work1")
    Work2 = QuantumRegister(work_size, "Work2")
    l_t = QuantumRegister(LT_WIDTH, "l_t")
    l_s = QuantumRegister(LS_WIDTH, "l_s")
    l_rp = QuantumRegister(LRP_WIDTH, "l_rp")
    Borrowed = QuantumRegister(1, "Borrowed")
    Scratch = QuantumRegister(10, "Scratch")
    qc = _e._block_circuit(Ctrl, Tail, Work1, Work2, l_t, l_s, l_rp,
                           Borrowed, Scratch, name=name)
    carry = Scratch[9]
    lrp_extended = list(l_rp) + [Borrowed[0]]
    affine_scratch = list(Scratch[:9]) + [carry]
    qc.append(_e.cuccaro_add_mod_2n_no_z_gate(LS_WIDTH, name="ADD_lrp8_to_ls9"),
              lrp_extended + list(l_s) + [carry])
    # The borrowed high addend contributes exactly 256 modulo 512.  Cancel it
    # without learning or changing the borrowed value.
    qc.cx(Borrowed[0], l_s[LS_WIDTH - 1])
    _e.const_minus_inplace(qc, l_s, n, affine_scratch)

    def selected_dirty_toggle() -> None:
        labels = list(range(0, work_size - 3))
        depth = _tight_unary_depth_for_labels(labels)

        def leaf(encoded_length: int, ej: Qubit) -> None:
            qc.ccx(ej, Work1[encoded_length + 2], Tail[0])

        unary_iteration_tight(
            qc, index_reg=l_t, labels=labels, ctrl=Ctrl[0],
            ancillas=list(Scratch[:depth]), leaf_fn=leaf, order="inc",
        )

    map_scratch = list(Scratch)
    selected_dirty_toggle()
    _upper_zero_map_borrowed(
        qc, ctrl=Ctrl[0], boundary_B=l_s, bits=Work2, dirty_map=Work1,
        borrowed=Borrowed[0], k=0, K=work_size - 1, scratch=map_scratch,
    )
    selected_dirty_toggle()
    _upper_zero_map_borrowed(
        qc, ctrl=Ctrl[0], boundary_B=l_s, bits=Work2, dirty_map=Work1,
        borrowed=Borrowed[0], k=0, K=work_size - 1, scratch=map_scratch,
    )

    _e.const_minus_inplace(qc, l_s, n, affine_scratch)
    qc.cx(Borrowed[0], l_s[LS_WIDTH - 1])
    qc.append(_e.cuccaro_sub_mod_2n_no_z_gate(LS_WIDTH, name="SUB_lrp8_from_ls9"),
              lrp_extended + list(l_s) + [carry])
    return _e._finalize_block(qc)


@lru_cache(maxsize=None)
def compact_lower_borrow_gate(*, n: int,
                              name: str = "T_LOWER_BORROW_COMPACT") -> Gate:
    work_size = n + 3
    Ctrl = QuantumRegister(1, "Ctrl")
    Tail = QuantumRegister(1, "Tail")
    Neg = QuantumRegister(1, "Neg")
    Work1 = QuantumRegister(work_size, "Work1")
    Work2 = QuantumRegister(work_size, "Work2")
    l_t = QuantumRegister(LT_WIDTH, "l_t")
    Borrowed = QuantumRegister(1, "Borrowed")
    Scratch = QuantumRegister(9, "Scratch")
    qc = _e._block_circuit(Ctrl, Tail, Neg, Work1, Work2, l_t,
                           Borrowed, Scratch, name=name)
    carry, active, eq = Scratch[:3]
    eq_pool = list(Scratch[3:])
    qc.ccx(Ctrl[0], Tail[0], active)

    def first_pass_cell(idx: int) -> None:
        addend = Work1[idx]
        target = Work2[idx]
        carry_in = carry if idx == 0 else Work1[idx - 1]
        qc.cx(carry_in, target)
        qc.cx(addend, carry_in)
        qc.ccx(carry_in, target, addend)

    for idx in range(work_size):
        first_pass_cell(idx)
        physical = idx + 1
        if 2 <= physical <= 257:
            _e.compute_eq_const(qc, l_t, physical - 2, eq, eq_pool)
            _borrowed_c3x(qc, active, eq, Work1[idx], Neg[0], Borrowed[0])
            _e.compute_eq_const(qc, l_t, physical - 2, eq, eq_pool)

    for idx in range(work_size - 1, -1, -1):
        addend = Work1[idx]
        target = Work2[idx]
        carry_in = carry if idx == 0 else Work1[idx - 1]
        qc.ccx(carry_in, target, addend)
        qc.cx(addend, carry_in)
        qc.cx(carry_in, target)
    qc.ccx(Ctrl[0], Tail[0], active)
    return _e._finalize_block(qc)

@lru_cache(maxsize=None)
def swap_work_and_len_unary_shared_gate(*, n: int, len_width: int, k4: int, K4: int,
                                        k5: int, K5: int, name: str = "SWAP_AND_LEN_S835_FAST") -> Gate:
    work_size = n + 3
    depth4 = _e.unary_depth(K4 - k4 + 1)
    depth5 = _e.unary_depth(K5 - k5 + 1)
    scratch4 = max(len_width + 1, depth4 + 2)
    scratch5 = max(len_width + 1, depth5 + 2)
    scratch_size = max(scratch4, scratch5)
    Ctrl = QuantumRegister(1, "Ctrl")
    Work1 = QuantumRegister(work_size, "Work1")
    Work2 = QuantumRegister(work_size, "Work2")
    l_t = QuantumRegister(len_width, "l_t")
    l_rp = QuantumRegister(len_width, "l_rp")
    Scratch = QuantumRegister(scratch_size, "Scratch")
    qc = _e._block_circuit(Ctrl, Work1, Work2, l_t, l_rp, Scratch, name=name)
    for i in range(work_size):
        _e.cswap_toffoli(qc, Ctrl[0], Work1[i], Work2[i])
    gate_lt = len_update_lt_unary_gate(n=n, k=k4, K=K4, len_width=len_width)
    _e._append_with_optional_clbits(qc, gate_lt, [Ctrl[0]] + list(Work1[k4 - 1:K4]) + list(Work2[k4 - 1:K4])
                                    + list(l_t) + list(l_rp) + list(Scratch[:scratch4]))
    gate_lrp = len_update_lrp_unary_gate(n=n, k=k5, K=K5, len_width=len_width)
    _e._append_with_optional_clbits(qc, gate_lrp, [Ctrl[0]] + list(Work1[k5 - 1:K5]) + list(Work2[k5 - 1:K5])
                                    + list(l_t) + list(l_rp) + list(Scratch[:scratch5]))
    return _e._finalize_block(qc)


def _fastdual_interval_scratch_size(n: int, k: int, K: int, len_width: int, shift_width: int) -> int:
    """Scratch size used by ``lc_interval_addsub_unary_gate``.

    This helper mirrors the scratch layout in ``lc_interval_addsub_unary_gate``.
    It is intentionally kept next to ``qiskit_paper_aux_size`` because the
    default Aux size used by the checkpointed counter must scale with this
    value.  For n=256 the worst case is 19 scratch qubits plus the temporary
    Ctrl bit, i.e. Aux=20.  For n=512 the unary path depth increases by one
    on each of the two endpoint scans, so the worst-case scratch is 21 and
    Aux must be 22.
    """
    if k > K:
        return 0
    endpoint_width = max(len_width, shift_width)
    rel_count = K - k + 1
    labels_main = list(range(rel_count))
    if rel_count > 1 and ((rel_count - 1) & (rel_count - 2)) == 0:
        # Same top-special split as lc_interval_addsub_unary_gate.
        labels_main = list(range(rel_count - 1))
    depth = _tight_unary_depth_for_labels(labels_main) if labels_main else 0
    base = max(2 * depth, endpoint_width)
    return base + 3


def _fastdual_prefix_scratch_size(k: int, K: int, len_width: int) -> int:
    if k > K:
        return 0
    depth = _e.unary_depth(K - k + 1)
    return max(depth, len_width) + 3


def _fastdual_interval_scratch_size(label_count: int, endpoint_width: int) -> int:
    """Scratch qubits used by lc_interval_addsub_unary_gate.

    The FASTDUAL interval Add/Sub block handles a one-more-than-a-power-of-two
    interval by pulling the top label out as a special endpoint.  Its two endpoint
    unary paths therefore have depth based on ``main_count`` rather than directly
    on ``label_count``.  The scratch layout in lc_interval_addsub_unary_gate is

        base = max(2*depth, endpoint_width)
        Scratch[base], Scratch[base+1], Scratch[base+2]

    so the number of scratch qubits needed by the block is ``base + 3``.
    This is 19 for n=256 but grows to 21 for n=384/512; the previous hard-coded
    lower bound of 19 caused the n=512 qubit-arity mismatch.
    """
    depth = _tight_unary_depth_for_labels(list(range(label_count))) if label_count > 1 else 0
    return max(2 * depth, endpoint_width) + 3


def fixed_schedule_shift_width(n: int, base_width: int, T_max: int) -> int:
    """Retain every post-terminal rotation without wrapping the pointer."""
    max_padding = max(1, T_max - 4 * n)
    return max(base_width, max_padding.bit_length())


def safe_active_windows(n: int, T: int) -> dict[str, tuple[int, int]]:
    """Return universally certified windows for secp256k1's fixed schedule."""
    if n == 256:
        if not 1 <= T <= len(_CERTIFIED_WINDOW_ROWS):
            raise ValueError(f"certified secp256k1 step out of range: {T}")
        row = _CERTIFIED_WINDOW_ROWS[T - 1]

        # A null certified window means the block control is unreachable on
        # every valid secp256k1 state at this step.  A singleton keeps the
        # generic controlled gate shape while adding no semantic assumption.
        def window(name: str) -> tuple[int, int]:
            value = row[name]
            return (1, 1) if value is None else (int(value[0]), int(value[1]))

        return {
            "r_addsub": window("r_addsub"),
            "swap": window("quotient_swap"),
            "t_addsub": window("t_addsub"),
            "len_update_lt": window("len_update_lt"),
            "len_update_lrp": window("len_update_lrp"),
        }
    try:
        return _e.active_windows(n, T)
    except ValueError:
        work_size = n + 3
        return {
            "r_addsub": (1, work_size),
            "swap": (1, work_size - 1),
            "t_addsub": (1, work_size),
            "len_update_lt": (1, work_size),
            "len_update_lrp": (1, work_size),
        }


@lru_cache(maxsize=None)
def compact_pre_shift_gate(*, work_size: int,
                           name: str = "PRE_SHIFT_MOD259") -> Gate:
    Phase1 = QuantumRegister(1, "Phase1")
    Phase2 = QuantumRegister(1, "Phase2")
    Work2 = QuantumRegister(work_size, "Work2")
    l_s = QuantumRegister(LS_WIDTH, "l_s")
    Scratch = QuantumRegister(10, "Scratch")
    qc = _e._block_circuit(Phase1, Phase2, Work2, l_s, Scratch, name=name)
    phase1_is0 = Scratch[0]
    both = Scratch[1]
    chain = list(Scratch[2:])

    qc.x(Phase1[0])
    qc.cx(Phase1[0], phase1_is0)
    qc.x(Phase1[0])
    for i in range(work_size - 1):
        _e.cswap_toffoli(qc, phase1_is0, Work2[i], Work2[i + 1])
    inc_mod259_1ctrl(qc, phase1_is0, l_s, chain)

    qc.ccx(phase1_is0, Phase2[0], both)
    _e.controlled_rotate_right_by_two(qc, both, list(Work2))
    dec_mod259_1ctrl(qc, both, l_s, chain)
    dec_mod259_1ctrl(qc, both, l_s, chain)
    qc.ccx(phase1_is0, Phase2[0], both)

    qc.x(Phase1[0])
    qc.cx(Phase1[0], phase1_is0)
    qc.x(Phase1[0])
    return _e._finalize_block(qc)


@lru_cache(maxsize=None)
def compact_post_shift_gate(*, work_size: int,
                            name: str = "POST_SHIFT_MOD259") -> Gate:
    Phase1 = QuantumRegister(1, "Phase1")
    Phase2 = QuantumRegister(1, "Phase2")
    Work2 = QuantumRegister(work_size, "Work2")
    l_s = QuantumRegister(LS_WIDTH, "l_s")
    Scratch = QuantumRegister(9, "Scratch")
    qc = _e._block_circuit(Phase1, Phase2, Work2, l_s, Scratch, name=name)
    both = Scratch[0]
    chain = list(Scratch[1:])

    for i in range(work_size - 1):
        _e.cswap_toffoli(qc, Phase1[0], Work2[i], Work2[i + 1])
    inc_mod259_1ctrl(qc, Phase1[0], l_s, chain)
    qc.ccx(Phase1[0], Phase2[0], both)
    _e.controlled_rotate_right_by_two(qc, both, list(Work2))
    dec_mod259_1ctrl(qc, both, l_s, chain)
    dec_mod259_1ctrl(qc, both, l_s, chain)
    qc.ccx(Phase1[0], Phase2[0], both)
    return _e._finalize_block(qc)


@lru_cache(maxsize=None)
def compact_phase_update_gate(name: str = "PHASE_UPDATE_COMPACT") -> Gate:
    Phase1 = QuantumRegister(1, "Phase1")
    Phase2 = QuantumRegister(1, "Phase2")
    Sign = QuantumRegister(1, "Sign")
    l_q = QuantumRegister(LQ_WIDTH, "l_q")
    l_rp = QuantumRegister(LRP_WIDTH, "l_rp")
    l_s = QuantumRegister(LS_WIDTH, "l_s")
    Scratch = QuantumRegister(11, "Scratch")
    qc = _e._block_circuit(Phase1, Phase2, Sign, l_q, l_rp, l_s, Scratch, name=name)
    z_lq, z_lrp, cond, tmp = Scratch[:4]
    pool = list(Scratch[4:])

    _e.compute_eq_const(qc, l_q, (1 << LQ_WIDTH) - 1, z_lq, pool)
    _e.compute_eq_const(qc, l_rp, LRP_ZERO, z_lrp, pool)
    qc.x(z_lrp)
    qc.ccx(z_lq, z_lrp, cond)
    qc.x(z_lrp)
    qc.cx(Sign[0], tmp)
    qc.cx(Phase1[0], tmp)
    qc.ccx(cond, tmp, Phase2[0])
    qc.cx(Phase1[0], tmp)
    qc.cx(Sign[0], tmp)
    qc.ccx(cond, Phase2[0], Sign[0])
    qc.x(z_lrp)
    qc.ccx(z_lq, z_lrp, cond)
    qc.x(z_lrp)
    _e.compute_eq_const(qc, l_rp, LRP_ZERO, z_lrp, pool)
    _e.compute_eq_const(qc, l_q, (1 << LQ_WIDTH) - 1, z_lq, pool)

    # Modulo-259 revisits the shift-zero sentinel during terminal padding.
    # Guard the phase transition with l_rp != 0 so padding remains frozen.
    _e.compute_eq_const(qc, l_s, LS_ZERO, z_lq, pool)
    _e.compute_eq_const(qc, l_rp, LRP_ZERO, z_lrp, pool)
    qc.x(z_lrp)
    qc.ccx(z_lq, z_lrp, cond)
    qc.x(z_lrp)
    qc.cx(cond, Phase1[0])
    qc.cx(cond, Phase2[0])
    qc.x(z_lrp)
    qc.ccx(z_lq, z_lrp, cond)
    qc.x(z_lrp)
    _e.compute_eq_const(qc, l_rp, LRP_ZERO, z_lrp, pool)
    _e.compute_eq_const(qc, l_s, LS_ZERO, z_lq, pool)
    return _e._finalize_block(qc)


def qiskit_paper_aux_size(n: int, len_width: int, shift_width: int, T_max: Optional[int] = None,
                          include_algorithm1: bool = False) -> int:
    if n != 256:
        raise ValueError("exact-width dirty12 route is certified only for secp256k1")
    return CLEAN_AUX_SIZE

def make_global_registers_noctrl(*, n: int, len_width: int, shift_width: int,
                                 T_max: Optional[int] = None, include_algorithm1: bool = False,
                                 aux_size: Optional[int] = None):
    work_size = n + 3
    Phase1 = QuantumRegister(1, "Phase1")
    Phase2 = QuantumRegister(1, "Phase2")
    Iter = QuantumRegister(1, "Iter")
    Sign = QuantumRegister(1, "Sign")
    Work1 = QuantumRegister(work_size, "Work1")
    Work2 = QuantumRegister(work_size, "Work2")
    l_t = QuantumRegister(LT_WIDTH, "l_t")
    l_q = QuantumRegister(LQ_WIDTH, "l_q")
    l_s = QuantumRegister(LS_WIDTH, "l_s")
    l_rp = QuantumRegister(LRP_WIDTH, "l_rp")
    if aux_size is None:
        aux_size = qiskit_paper_aux_size(n, len_width, shift_width, T_max, include_algorithm1)
    if aux_size != CLEAN_AUX_SIZE:
        raise ValueError(f"exact-width route requires Aux={CLEAN_AUX_SIZE}")
    Aux = QuantumRegister(aux_size, "Aux")
    Dirty = QuantumRegister(DIRTY_PASSENGER_SIZE, "DirtyPassenger")
    return Phase1, Phase2, Iter, Sign, Work1, Work2, l_t, l_q, l_s, l_rp, Aux, Dirty


def _make_condition(qc: QuantumCircuit, conditions, out: Qubit, scratch: Sequence[Qubit]) -> None:
    _e.compute_control(qc, conditions, out, scratch)


def _toggle_live_r_phase(qc: QuantumCircuit, *, phase1: Qubit,
                         l_rp: Sequence[Qubit], out: Qubit,
                         scratch: Sequence[Qubit]) -> None:
    """Toggle ``out`` by ``l_rp != 0 and phase1 == 0`` on valid EEA states.

    Length zero is encoded as all ones.  The Algorithm-3 terminal transition
    produces Phase1=Phase2=Sign=0, and padding preserves those controls.  Thus
    terminal and Phase1=1 are mutually exclusive on the block domain, making

        1 xor Phase1 xor [l_rp == 0]

    equal to ``[l_rp != 0] and not Phase1``.  Every operation targets ``out``,
    so a second invocation cleans it exactly.
    """
    qc.x(out)
    qc.cx(phase1, out)
    _e.compute_eq_const(qc, l_rp, (1 << len(l_rp)) - 1, out, scratch)


def append_one_step_T(qc: QuantumCircuit, *, T: int, n: int, len_width: int, shift_width: int,
                      Phase1, Phase2, Iter, Sign, Work1, Work2, l_t, l_q, l_s, l_rp,
                      Aux, Dirty) -> None:
    work_size = n + 3
    windows = safe_active_windows(n, T)
    k1, K1 = windows["r_addsub"]
    # The certified secp256k1 table already includes the live carry/sign lane.
    # Small-width fallback tests retain the historical one-lane repair.
    if n != 256:
        k1 = max(1, k1 - 1)
    k2, K2 = windows["swap"]
    k3, K3 = windows["t_addsub"]
    k4, K4 = windows["len_update_lt"]
    k5, K5 = windows["len_update_lrp"]
    ctrl = Aux[0]
    scratch = list(Aux[1:])
    pool = scratch
    # Pre-shift
    pre = compact_pre_shift_gate(work_size=work_size)
    _e._append_with_optional_clbits(qc, pre, [Phase1[0], Phase2[0]] + list(Work2)
                                    + list(l_s) + scratch[:10])
    # Terminal padding must only rotate Work2.  Fold l_rp!=0 and Phase1=0 into
    # the existing control and retain it across the complete R sequence.
    _toggle_live_r_phase(qc, phase1=Phase1[0], l_rp=l_rp, out=ctrl, scratch=scratch)
    rsub = compact_interval_addsub_gate(n=n, k=k1, K=K1,
                                        mode="sub", sign_update=True, target="work1", name="R_SUB_COMPACT")
    _e._append_with_optional_clbits(qc, rsub, [ctrl, Sign[0]] + list(Work1[k1-1:K1]) + list(Work2[k1-1:K1])
                                    + list(l_t) + list(l_q) + list(l_s) + list(Dirty) + scratch)
    # if the live R phase also has Phase2=1 then Sign ^= 1
    qc.ccx(ctrl, Phase2[0], Sign[0])
    # Convert ctrl from live-R to the R-add predicate by toggling it when
    # Phase1=0 and Phase2&Sign=1.  The clean C3X scratch is restored by the
    # primitive and remains available to the interval adder.
    qc.x(Phase1[0])
    _dirty_c3x(qc, Phase1[0], Phase2[0], Sign[0], ctrl, scratch[0])
    qc.x(Phase1[0])
    radd = compact_interval_addsub_gate(n=n, k=k1, K=K1,
                                        mode="add", sign_update=False, target="work1", name="R_ADD_COMPACT")
    _e._append_with_optional_clbits(qc, radd, [ctrl, Sign[0]] + list(Work1[k1-1:K1]) + list(Work2[k1-1:K1])
                                    + list(l_t) + list(l_q) + list(l_s) + list(Dirty) + scratch)
    qc.x(Phase1[0])
    _dirty_c3x(qc, Phase1[0], Phase2[0], Sign[0], ctrl, scratch[0])
    qc.x(Phase1[0])
    _toggle_live_r_phase(qc, phase1=Phase1[0], l_rp=l_rp, out=ctrl, scratch=scratch)
    # Swap: ctrl = Phase1 xor Phase2
    qc.cx(Phase1[0], ctrl); qc.cx(Phase2[0], ctrl)
    lcs = compact_lc_swap_gate(k=k2, K=K2)
    _e._append_with_optional_clbits(qc, lcs, [ctrl, Phase1[0], Sign[0]] + list(Work1[k2-1:K2+1]) + list(l_t) + list(l_q)
                                    + scratch[:lcs.num_qubits-(3+(K2-k2+2)+LT_WIDTH+LQ_WIDTH)])
    qc.cx(Phase2[0], ctrl); qc.cx(Phase1[0], ctrl)
    # l_q +/- updates.
    _make_condition(qc, [(Phase1[0], 1), (Phase2[0], 0)], ctrl, scratch)
    _e.dec_mod2n_1ctrl(qc, ctrl, list(l_q), scratch[:max(0,len_width-1)])
    _make_condition(qc, [(Phase1[0], 1), (Phase2[0], 0)], ctrl, scratch)
    _make_condition(qc, [(Phase1[0], 0), (Phase2[0], 1)], ctrl, scratch)
    _e.inc_mod2n_1ctrl(qc, ctrl, list(l_q), scratch[:max(0,len_width-1)])
    _make_condition(qc, [(Phase1[0], 0), (Phase2[0], 1)], ctrl, scratch)
    # Retain the tail predicate across the cancelling subtract/add pair.  The
    # restoring add captures its selected carry into Sign only after the
    # Sign-dependent subtract control has been uncomputed.
    tail_zero = scratch[-2]
    t_pool = [lane for lane in scratch if lane != tail_zero]
    tail_gate = compact_tail_zero_gate(n=n)
    tail_args = [Phase1[0], tail_zero] + list(Work1) + list(Work2) + list(l_t) + list(l_s) + list(l_rp) + [Dirty[1]] + [ctrl] + scratch[:-2]
    _e._append_with_optional_clbits(qc, tail_gate, tail_args)
    # T sub condition: Phase1=1 and (Phase2=1 or Sign=0)
    tmp = scratch[0]
    _make_condition(qc, [(Phase2[0], 0), (Sign[0], 1)], tmp, scratch[1:])
    _make_condition(qc, [(Phase1[0], 1), (tmp, 0)], ctrl, scratch[1:])
    _make_condition(qc, [(Phase2[0], 0), (Sign[0], 1)], tmp, scratch[1:])
    tsub = compact_prefix_addsub_gate(k=k3, K=K3,
                                      mode="sub", sign_update=False,
                                      capture_borrow_sign=False,
                                      target="work2", name="T_SUB_COMPACT")
    _e._append_with_optional_clbits(qc, tsub, [ctrl, Sign[0], tail_zero] + list(Work1[k3-1:K3]) + list(Work2[k3-1:K3])
                                    + list(l_t) + [Dirty[3]]
                                    + t_pool[:tsub.num_qubits-(3+2*(K3-k3+1)+LT_WIDTH+1)])
    _make_condition(qc, [(Phase2[0], 0), (Sign[0], 1)], tmp, scratch[1:])
    _make_condition(qc, [(Phase1[0], 1), (tmp, 0)], ctrl, scratch[1:])
    _make_condition(qc, [(Phase2[0], 0), (Sign[0], 1)], tmp, scratch[1:])
    qc.cx(Phase1[0], Sign[0])
    _make_condition(qc, [(Phase1[0], 1)], ctrl, scratch)
    tadd = compact_prefix_addsub_gate(k=k3, K=K3,
                                      mode="add", sign_update=False,
                                      capture_borrow_sign=True,
                                      target="work2", name="T_ADD_COMPACT")
    _e._append_with_optional_clbits(qc, tadd, [ctrl, Sign[0], tail_zero] + list(Work1[k3-1:K3]) + list(Work2[k3-1:K3])
                                    + list(l_t) + [Dirty[3]]
                                    + t_pool[:tadd.num_qubits-(3+2*(K3-k3+1)+LT_WIDTH+1)])
    _make_condition(qc, [(Phase1[0], 1)], ctrl, scratch)
    _e._append_with_optional_clbits(qc, tail_gate, tail_args)
    # Post-shift
    post = compact_post_shift_gate(work_size=work_size)
    _e._append_with_optional_clbits(qc, post, [Phase1[0], Phase2[0]] + list(Work2)
                                    + list(l_s) + scratch[:9])
    # Phase update
    pupdate = compact_phase_update_gate()
    _e._append_with_optional_clbits(qc, pupdate, [Phase1[0], Phase2[0], Sign[0]]
                                    + list(l_q) + list(l_rp) + list(l_s) + scratch)
    # End iteration every four steps.
    if T % 4 == 0:
        z_lq = scratch[0]; z_ls = scratch[1]; eq_pool = scratch[2:]
        _e.compute_eq_const(qc, l_q, (1 << LQ_WIDTH) - 1, z_lq, eq_pool)
        _e.compute_eq_const(qc, l_s, LS_ZERO, z_ls, eq_pool)
        # Termination is aligned to a four-step boundary.  During terminal
        # padding l_s returns to its modulo-259 zero sentinel only at offsets
        # 259 and 518; neither is divisible by four, and the certified horizon
        # is shorter than the 1036-step joint recurrence.  Therefore the
        # original two-flag end trigger remains exact without an l_rp guard.
        qc.ccx(z_lq, z_ls, ctrl)
        # The compact decoder needs all eleven Aux scratch lanes clean.  Keep
        # only the trigger bit live and serialize both equality tests away.
        _e.compute_eq_const(qc, l_s, LS_ZERO, z_ls, eq_pool)
        _e.compute_eq_const(qc, l_q, (1 << LQ_WIDTH) - 1, z_lq, eq_pool)
        # The original Section 4.5 bounds are unsafe.  These ranges come from
        # the pinned continuant certificate above; small-width tests still use
        # full scans because the certificate is specific to secp256k1.
        if n != 256:
            k4, K4, k5, K5 = 1, work_size, 1, work_size
        swlen = compact_swap_work_and_len_gate(
            n=n, k4=k4, K4=K4, k5=k5, K5=K5,
        )
        _e._append_with_optional_clbits(qc, swlen, [ctrl] + list(Work1) + list(Work2)
                                        + list(l_t) + list(l_rp) + [Dirty[4]] + scratch)
        qc.cx(ctrl, Iter[0])
        _e.compute_eq_const(qc, l_q, (1 << LQ_WIDTH) - 1, z_lq, eq_pool)
        _e.compute_eq_const(qc, l_s, LS_ZERO, z_ls, eq_pool)
        qc.ccx(z_lq, z_ls, ctrl)
        _e.compute_eq_const(qc, l_s, LS_ZERO, z_ls, eq_pool)
        _e.compute_eq_const(qc, l_q, (1 << LQ_WIDTH) - 1, z_lq, eq_pool)


def build_step_circuit(n:int, T:int, *, T_max:Optional[int]=None, aux_size:Optional[int]=None, measurement_uncompute:bool=True):
    cfg=get_n_config(n); lw=int(cfg['len_width']); T_max=int(T_max or cfg['T_max'])
    sw=LS_WIDTH
    if aux_size is None: aux_size=qiskit_paper_aux_size(n,lw,sw,T_max)
    set_measurement_uncompute(measurement_uncompute)
    regs=make_global_registers_noctrl(n=n,len_width=lw,shift_width=sw,T_max=T_max,aux_size=aux_size)
    qc=QuantumCircuit(*regs, name=f"S835_FASTDUAL_STEP_T{T}_{n}")
    Phase1,Phase2,Iter,Sign,Work1,Work2,l_t,l_q,l_s,l_rp,Aux,Dirty=regs
    append_one_step_T(qc,T=T,n=n,len_width=lw,shift_width=sw,Phase1=Phase1,Phase2=Phase2,Iter=Iter,Sign=Sign,Work1=Work1,Work2=Work2,l_t=l_t,l_q=l_q,l_s=l_s,l_rp=l_rp,Aux=Aux,Dirty=Dirty)
    return qc

if __name__ == '__main__':
    import argparse,json
    ap=argparse.ArgumentParser(); ap.add_argument('--n',type=int,default=256); ap.add_argument('--T',type=int,default=1); ap.add_argument('--count',action='store_true'); args=ap.parse_args()
    cfg=get_n_config(args.n); lw=int(cfg['len_width']); Tm=int(cfg['T_max'])
    sw=LS_WIDTH
    out={'n':args.n,'l_t_width':LT_WIDTH,'l_q_width':LQ_WIDTH,'l_s_width':LS_WIDTH,
         'l_rp_width':LRP_WIDTH,'T_max':Tm,'aux_size':qiskit_paper_aux_size(args.n,lw,sw,Tm),
         'dirty_passenger_size':DIRTY_PASSENGER_SIZE}
    qc=build_step_circuit(args.n,args.T,T_max=Tm)
    out['step_qubits']=qc.num_qubits; out['top_ops']={str(k):int(v) for k,v in qc.count_ops().items()}
    if args.count:
        out['ops']={str(k):int(v) for k,v in _e.count_circuit_ops_recursive(qc).items()}
    print(json.dumps(out,indent=2,sort_keys=True))
