#!/usr/bin/env python3
"""Verify and aggregate the certified paper2607 primitive shards."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import struct
import subprocess


MAGIC = b"P26EEA2\0"
FIELD_WIDTH = 256
LOCAL_WIDTH = 578
SCHEDULE_STEPS = 1616
SOURCE_MODULE = "eea_circuit_s835_exactwidth_dirty12"
AUX_SIZE = 12
NAME = re.compile(r"chunk-(\d{4})-(\d{4})\.zst$")


def read_exact(stream, size: int) -> bytes:
    value = stream.read(size)
    if len(value) != size:
        raise AssertionError(f"truncated stream: wanted {size}, got {len(value)}")
    return value


def verify_chunk(path: Path, expected_start: int) -> tuple[dict[str, object], int]:
    match = NAME.fullmatch(path.name)
    if match is None:
        raise AssertionError(f"unexpected chunk name: {path.name}")
    name_start, name_end = map(int, match.groups())
    report_path = path.with_suffix(path.suffix + ".json")
    report = json.loads(report_path.read_text(encoding="utf-8"))

    process = subprocess.Popen(
        ["zstd", "-q", "-dc", str(path)],
        stdout=subprocess.PIPE,
    )
    assert process.stdout is not None
    header = read_exact(process.stdout, 24)
    if header[:8] != MAGIC:
        raise AssertionError(f"{path.name}: wrong magic")
    field_width, local_width, start, end = struct.unpack("<IIII", header[8:])
    if (field_width, local_width) != (FIELD_WIDTH, LOCAL_WIDTH):
        raise AssertionError(f"{path.name}: wrong widths {(field_width, local_width)}")
    if (start, end) != (name_start, name_end):
        raise AssertionError(f"{path.name}: header/name range mismatch")
    if start != expected_start or not start <= end <= SCHEDULE_STEPS:
        raise AssertionError(f"{path.name}: noncontiguous range")

    digest = hashlib.sha256()
    payload_bytes = 0
    while True:
        block = process.stdout.read(8 * 1024 * 1024)
        if not block:
            break
        digest.update(block)
        payload_bytes += len(block)
    return_code = process.wait()
    if return_code != 0:
        raise AssertionError(f"{path.name}: zstd exited {return_code}")
    if payload_bytes % 8:
        raise AssertionError(f"{path.name}: partial primitive record")

    records = payload_bytes // 8
    checks = {
        "schema": "paper2607-eea-primitive-stream-v3",
        "n": FIELD_WIDTH,
        "qubits": LOCAL_WIDTH,
        "source_module": SOURCE_MODULE,
        "aux_size": AUX_SIZE,
        "step_start": start,
        "step_end": end,
        "schedule_end": SCHEDULE_STEPS,
        "measurement_uncompute": False,
        "records": records,
        "raw_record_sha256": digest.hexdigest(),
        "compressed_bytes": path.stat().st_size,
    }
    for key, expected in checks.items():
        if report.get(key) != expected:
            raise AssertionError(
                f"{path.name}: report {key}={report.get(key)!r}, expected {expected!r}"
            )

    per_step = report.get("per_step")
    if not isinstance(per_step, list) or len(per_step) != end - start + 1:
        raise AssertionError(f"{path.name}: malformed per-step report")
    if [int(row["step"]) for row in per_step] != list(range(start, end + 1)):
        raise AssertionError(f"{path.name}: noncontiguous per-step report")
    if sum(int(row["records"]) for row in per_step) != records:
        raise AssertionError(f"{path.name}: per-step record total mismatch")
    report_counts = Counter({key: int(value) for key, value in report["counts"].items()})
    per_step_counts: Counter[str] = Counter()
    for row in per_step:
        per_step_counts.update({key: int(value) for key, value in row["counts"].items()})
    if per_step_counts != report_counts:
        raise AssertionError(f"{path.name}: per-step primitive counts mismatch")
    if int(report["executed_toffoli"]) != (
        report_counts["ccx"] + 2 * report_counts["clean_c3x_mbu"]
    ):
        raise AssertionError(f"{path.name}: executed Toffoli total mismatch")

    report["file"] = path.name
    report["compressed_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    return report, end + 1


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--directory",
        type=Path,
        default=Path(__file__).resolve().parent,
    )
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    paths = sorted(args.directory.glob("chunk-*.zst"))
    expected_start = 1
    reports = []
    totals: Counter[str] = Counter()
    for path in paths:
        report, expected_start = verify_chunk(path, expected_start)
        reports.append(report)
        totals["records"] += int(report["records"])
        for kind, count in report["counts"].items():
            totals[kind] += int(count)

    if expected_start != SCHEDULE_STEPS + 1:
        raise AssertionError(
            f"incomplete schedule: next step {expected_start}, expected {SCHEDULE_STEPS + 1}"
        )
    if len(reports) != 36:
        raise AssertionError(f"wrong chunk count: {len(reports)}")

    kind7 = totals["clean_c3x_mbu"]
    aggregate = {
        "schema": "paper2607-eea-primitive-stream-aggregate-v1",
        "field_width": FIELD_WIDTH,
        "local_width": LOCAL_WIDTH,
        "source_module": SOURCE_MODULE,
        "aux_size": AUX_SIZE,
        "schedule_steps": SCHEDULE_STEPS,
        "chunk_count": len(reports),
        "records_per_traversal": totals["records"],
        "emitted_ops_per_traversal": totals["records"] + 3 * kind7,
        "executed_toffoli_per_traversal": totals["ccx"] + 2 * kind7,
        "four_traversal_emitted_ops": 4 * (totals["records"] + 3 * kind7),
        "four_traversal_executed_toffoli": 4 * (totals["ccx"] + 2 * kind7),
        "primitive_counts": {
            key: totals[key]
            for key in ("x", "cx", "ccx", "clean_c3x_mbu")
        },
        "chunks": [
            {
                key: report[key]
                for key in (
                    "file",
                    "step_start",
                    "step_end",
                    "records",
                    "raw_record_sha256",
                    "compressed_bytes",
                    "compressed_sha256",
                )
            }
            for report in reports
        ],
    }
    encoded = json.dumps(aggregate, indent=2, sort_keys=True) + "\n"
    if args.out is not None:
        args.out.write_text(encoded, encoding="utf-8")
    print(encoded, end="")


if __name__ == "__main__":
    main()
