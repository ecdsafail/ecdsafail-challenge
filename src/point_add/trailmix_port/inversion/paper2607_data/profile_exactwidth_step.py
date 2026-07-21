#!/usr/bin/env python3
"""Deterministic recursive profile for one exact-width paper2607 step."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import subprocess

import eea_circuit_s835_exactwidth_dirty12 as eea
import eea_circuit_updated as support


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def count_definition(circuit) -> Counter[str]:
    counts: Counter[str] = Counter()
    for item in circuit.data:
        operation = item.operation
        if operation.name in {"x", "cx", "ccx", "h", "cz", "measure", "reset", "u"}:
            counts[operation.name] += 1
        elif operation.definition is None:
            raise ValueError(f"unsupported primitive {operation.name!r}")
        else:
            counts.update(count_definition(operation.definition))
    return counts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--step", type=int, default=1470)
    args = parser.parse_args()

    circuit = eea.build_step_circuit(
        256,
        args.step,
        T_max=1616,
        aux_size=eea.CLEAN_AUX_SIZE,
        measurement_uncompute=False,
    )
    components: dict[str, Counter[str]] = defaultdict(Counter)
    for item in circuit.data:
        operation = item.operation
        if operation.definition is None:
            counts = Counter({operation.name: 1})
        else:
            counts = count_definition(operation.definition)
        components[operation.name].update(counts)

    support_path = Path(support.__file__).resolve()
    generator_path = Path(eea.__file__).resolve()
    upstream = support_path.parent
    commit = subprocess.check_output(
        ["git", "-C", str(upstream), "rev-parse", "HEAD"], text=True,
    ).strip()
    total = Counter(support.count_circuit_ops_recursive(circuit))
    report = {
        "schema": "paper2607-exactwidth-step-profile-v1",
        "step": args.step,
        "referenced_qubits": circuit.num_qubits,
        "clean_aux": eea.CLEAN_AUX_SIZE,
        "dirty_passenger": eea.DIRTY_PASSENGER_SIZE,
        "support": {
            "commit": commit,
            "path": str(support_path),
            "sha256": sha256(support_path),
        },
        "generator": {
            "path": str(generator_path),
            "sha256": sha256(generator_path),
        },
        "total": dict(sorted(total.items())),
        "components": {
            name: dict(sorted(counts.items()))
            for name, counts in sorted(components.items())
        },
    }
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
