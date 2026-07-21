#!/usr/bin/env python3
"""Flatten the paper's fixed-width EEA steps into a compact primitive stream."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib
import json
from pathlib import Path
import struct
import sys

import compression.zstd as zstd


KIND = {
    "x": 1,
    "cx": 2,
    "ccx": 3,
    "z": 4,
    "cz": 5,
    "swap": 6,
    "clean_c3x_mbu": 7,
}


def flatten(circuit, qmap=None):
    if qmap is None:
        qmap = {q: i for i, q in enumerate(circuit.qubits)}
    for item in circuit.data:
        op = item.operation
        qargs = [qmap[q] for q in item.qubits]
        name = op.name.lower()
        if name == "clean_c3x_mbu":
            yield name, qargs
            continue
        if item.clbits:
            raise RuntimeError(f"classical operands in unitary stream: {op.name}")
        if name in KIND:
            yield name, qargs
            continue
        definition = op.definition
        if definition is None:
            raise RuntimeError(f"opaque operation {op.name!r}")
        if definition.num_clbits:
            raise RuntimeError(f"dynamic definition in unitary stream: {op.name}")
        child_map = {q: qargs[i] for i, q in enumerate(definition.qubits)}
        yield from flatten(definition, child_map)


def pack_record(name: str, qargs: list[int]) -> bytes:
    if len(qargs) > 5:
        raise RuntimeError(f"primitive {name} has {len(qargs)} operands")
    q = qargs + [0, 0, 0, 0, 0]
    if any(x >= 1024 for x in qargs):
        raise RuntimeError(f"qubit index overflow in {name}: {qargs}")
    word = (KIND[name] | (len(qargs) << 4) | (q[0] << 8) | (q[1] << 18)
            | (q[2] << 28) | (q[3] << 38) | (q[4] << 48))
    return struct.pack("<Q", word)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--paper", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--start", type=int, default=1)
    ap.add_argument("--end", type=int, default=1481)
    ap.add_argument("--schedule-end", type=int, default=1481)
    ap.add_argument("--level", type=int, default=12)
    ap.add_argument("--module", default="eea_circuit_s835_fastdual")
    ap.add_argument("--aux-size", type=int, default=23)
    ap.add_argument("--expected-qubits", type=int, default=581)
    args = ap.parse_args()

    sys.path.insert(0, str(args.paper))
    eea = importlib.import_module(args.module)

    if args.start < 1 or args.end < args.start or args.end > args.schedule_end:
        raise SystemExit("invalid step range")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    counts = Counter()
    primitive_records = 0
    digest = hashlib.sha256()
    per_step = []
    with zstd.open(args.out, "wb", level=args.level) as stream:
        stream.write(b"P26EEA2\0")
        stream.write(struct.pack("<IIII", 256, args.expected_qubits, args.start, args.end))
        for step in range(args.start, args.end + 1):
            qc = eea.build_step_circuit(
                256,
                step,
                T_max=args.schedule_end,
                aux_size=args.aux_size,
                measurement_uncompute=False,
            )
            if qc.num_qubits != args.expected_qubits:
                raise RuntimeError(
                    f"step {step} width {qc.num_qubits}, expected {args.expected_qubits}"
                )
            step_count = 0
            step_counts = Counter()
            for name, qargs in flatten(qc):
                record = pack_record(name, qargs)
                stream.write(record)
                digest.update(record)
                counts[name] += 1
                step_counts[name] += 1
                primitive_records += 1
                step_count += 1
            per_step.append({"step": step, "records": step_count, "counts": dict(step_counts)})
            if step == args.start or step == args.end or step % 32 == 0:
                print(
                    f"step={step} records={step_count} total={primitive_records} "
                    f"ccx={counts['ccx']}",
                    flush=True,
                )

    report = {
        "schema": "paper2607-eea-primitive-stream-v3",
        "source_module": args.module,
        "n": 256,
        "qubits": args.expected_qubits,
        "aux_size": args.aux_size,
        "step_start": args.start,
        "step_end": args.end,
        "schedule_end": args.schedule_end,
        "measurement_uncompute": False,
        "record_bytes": 8,
        "records": primitive_records,
        "counts": dict(counts),
        "executed_toffoli": counts["ccx"] + 2 * counts["clean_c3x_mbu"],
        "raw_record_sha256": digest.hexdigest(),
        "compressed_bytes": args.out.stat().st_size,
        "per_step": per_step,
    }
    report_path = args.out.with_suffix(args.out.suffix + ".json")
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "per_step"}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
