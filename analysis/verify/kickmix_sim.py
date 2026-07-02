#!/usr/bin/env python3
"""A minimal, spec-faithful simulator for kickmix circuit files (.kmx).

Kickmix is the circuit format used by this challenge's source paper (Babbush et
al. 2026, arXiv:2603.28846v2; format spec in
original/zkp_ecc_zenodo_v2/docs/kickmix_{file_format,instruction_set}.md).
Because kickmix gates never create superpositions, a circuit can be simulated
along a single classical trajectory while tracking three things (per the spec):

  1. qubits : one bit per qubit id
  2. bits   : one bit per classical bit id
  3. phase  : +1 / -1 global sign (fuzz testing checks it stays +1)
  4. a condition stack for PUSH_CONDITION / POP_CONDITION

This mirrors exactly the semantics this repo's Rust simulator (src/sim.rs)
implements; re-deriving it independently in Python lets us validate the paper's
reference circuits (and reject its negative controls) without trusting either
side's code. Measurement (HMR) is simulated with a supplied RNG per the spec's
fuzz-testing rule: a random result, with phase kickback when the measured qubit
is ON and the result is ON.

Verified equivalent to the reference kickmix simulator
`original/zkp_ecc_zenodo_v2/lib/src/sim.rs` (Zenodo release for arXiv:2603.28846v2),
instruction for instruction: R = random phase-kickback then reset (`phase ^=
qubit & rng`, `qubit &= !cond`); HMR = random result, same phase kickback, reset;
CCX/CX/CZ/CCZ/Z/NEG/SWAP; the PushCondition/PopCondition stack; and all qubits
initialized to |0>. (The reference packs 64 shots into a u64 and tracks phase as
an XOR bit; here one trajectory tracks phase as +/-1 — equivalent per shot.)
"""
import re


class ParseError(Exception):
    pass


class Circuit:
    def __init__(self):
        self.instructions = []          # (name, targets:list[(kind,id)], cond:(k,id)|None)
        self.registers = {}             # reg_id -> list[(kind, id)] in increasing significance
        self.qubit_ids = set()
        self.bit_ids = set()

    @staticmethod
    def _target(tok):
        m = re.fullmatch(r"([qbr])(\d+)", tok)
        if not m:
            raise ParseError(f"bad target token: {tok!r}")
        return (m.group(1), int(m.group(2)))

    @classmethod
    def parse(cls, text):
        c = cls()
        for raw in text.splitlines():
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            parts = line.split()
            name = parts[0]
            if not re.fullmatch(r"[A-Z][A-Z0-9_]*", name):
                raise ParseError(f"bad instruction name: {name!r}")
            cond = None
            rest = parts[1:]
            if "if" in rest:
                i = rest.index("if")
                if i + 1 >= len(rest):
                    raise ParseError(f"'if' without condition bit in: {line!r}")
                cond = cls._target(rest[i + 1])
                if cond[0] != "b":
                    raise ParseError(f"condition must be a bit: {line!r}")
                rest = rest[:i]
            targets = [cls._target(t) for t in rest]
            for k, i in targets:
                if k == "q":
                    c.qubit_ids.add(i)
                elif k == "b":
                    c.bit_ids.add(i)
            if name == "APPEND_TO_REGISTER":
                (elem_k, elem_i), (rk, ri) = targets[0], targets[1]
                if rk != "r":
                    raise ParseError(f"APPEND_TO_REGISTER needs a register: {line!r}")
                c.registers.setdefault(ri, []).append((elem_k, elem_i))
            elif name == "REGISTER":
                (rk, ri) = targets[0]
                c.registers.setdefault(ri, [])
            if cond is not None:
                c.bit_ids.add(cond[1])
            c.instructions.append((name, targets, cond))
        return c

    # ---- register <-> integer helpers (little-endian, increasing significance) ----
    def reg_bits(self, reg_id):
        return self.registers[reg_id]

    def load_register(self, state, reg_id, value):
        for pos, (k, i) in enumerate(self.registers[reg_id]):
            b = (value >> pos) & 1
            (state.qubits if k == "q" else state.bits)[i] = b

    def read_register(self, state, reg_id):
        v = 0
        for pos, (k, i) in enumerate(self.registers[reg_id]):
            b = (state.qubits if k == "q" else state.bits).get(i, 0)
            v |= b << pos
        return v


class State:
    def __init__(self, n_qubits, n_bits):
        self.qubits = {i: 0 for i in range(n_qubits)}
        self.bits = {i: 0 for i in range(n_bits)}
        self.phase = 1
        self.cond = []


def simulate(circuit, state, rng):
    """Execute `circuit` in-place on `state`. `rng` must have random()->[0,1)."""
    q, b = state.qubits, state.bits

    def coin():
        return 1 if rng.random() < 0.5 else 0

    for name, targets, cond in circuit.instructions:
        # Metadata / control-flow ops ignore the condition stack (per spec).
        if name in ("APPEND_TO_REGISTER", "REGISTER", "DEBUG_PRINT"):
            continue
        if name == "PUSH_CONDITION":
            state.cond.append(bool(b.get(cond[1], 0)) if cond else True)
            continue
        if name == "POP_CONDITION":
            if state.cond:
                state.cond.pop()
            continue

        active = all(state.cond)
        if cond is not None:
            active = active and bool(b.get(cond[1], 0))
        if not active:
            continue

        ids = [i for (_k, i) in targets]
        if name == "X":
            q[ids[0]] ^= 1
        elif name == "CX":
            if q[ids[0]]:
                q[ids[1]] ^= 1
        elif name == "CCX":
            if q[ids[0]] and q[ids[1]]:
                q[ids[2]] ^= 1
        elif name == "SWAP":
            q[ids[0]], q[ids[1]] = q[ids[1]], q[ids[0]]
        elif name == "Z":
            if q[ids[0]]:
                state.phase = -state.phase
        elif name == "CZ":
            if q[ids[0]] and q[ids[1]]:
                state.phase = -state.phase
        elif name == "CCZ":
            if q[ids[0]] and q[ids[1]] and q[ids[2]]:
                state.phase = -state.phase
        elif name == "NEG":
            state.phase = -state.phase
        elif name == "R":
            # Reset. If the qubit was ON, the spec says the phase is randomized
            # (misuse); model that so bad circuits fail fuzzing.
            if q[ids[0]] == 1 and coin():
                state.phase = -state.phase
            q[ids[0]] = 0
        elif name == "HMR":
            # X-basis demolition measurement with fuzz-style random result and
            # phase kickback (spec: negate phase iff qubit ON and result ON).
            result = coin()
            if q[ids[0]] == 1 and result == 1:
                state.phase = -state.phase
            b[ids[1]] = result
            q[ids[0]] = 0
        elif name == "BIT_INVERT":
            b[ids[0]] ^= 1
        elif name == "BIT_STORE0":
            b[ids[0]] = 0
        elif name == "BIT_STORE1":
            b[ids[0]] = 1
        else:
            raise ParseError(f"unsupported instruction: {name}")
    return state
