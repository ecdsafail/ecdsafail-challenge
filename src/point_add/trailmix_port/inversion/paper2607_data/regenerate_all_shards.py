#!/usr/bin/env python3
"""Regenerate every paper2607 shard and publish them only after exact verification."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


SCHEDULE_END = 1616
CHUNK_SIZE = 45
EXPECTED = {
    "records_per_traversal": 138_231_827,
    "emitted_ops_per_traversal": 149_005_883,
    "executed_toffoli_per_traversal": 55_826_377,
    "four_traversal_emitted_ops": 596_023_532,
    "four_traversal_executed_toffoli": 223_305_508,
    "primitive_counts": {
        "x": 43_779_600,
        "cx": 42_217_202,
        "ccx": 48_643_673,
        "clean_c3x_mbu": 3_591_352,
    },
}


def ranges() -> list[tuple[int, int]]:
    return [
        (start, min(start + CHUNK_SIZE - 1, SCHEDULE_END))
        for start in range(1, SCHEDULE_END + 1, CHUNK_SIZE)
    ]


def generate_one(
    generator: Path,
    paper: Path,
    output_dir: Path,
    start: int,
    end: int,
) -> str:
    name = f"chunk-{start:04d}-{end:04d}.zst"
    output = output_dir / name
    command = [
        sys.executable,
        str(generator),
        "--paper",
        str(paper),
        "--out",
        str(output),
        "--start",
        str(start),
        "--end",
        str(end),
        "--schedule-end",
        str(SCHEDULE_END),
        "--module",
        "eea_circuit_s835_fastdual_aux22",
        "--aux-size",
        "22",
        "--expected-qubits",
        "581",
        "--level",
        "12",
    ]
    log = output.with_suffix(".log")
    with log.open("wb") as stream:
        result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
    if result.returncode:
        tail = log.read_text(errors="replace")[-4000:]
        raise RuntimeError(f"{name} failed with exit {result.returncode}:\n{tail}")
    report = json.loads(output.with_suffix(output.suffix + ".json").read_text())
    return f"{name}: records={report['records']} T={report['executed_toffoli']}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--paper", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=3)
    args = parser.parse_args()
    if args.jobs < 1:
        raise SystemExit("--jobs must be positive")

    directory = Path(__file__).resolve().parent
    generator = directory / "generate_eea_blob.py"
    verifier = directory / "verify_stream.py"
    with tempfile.TemporaryDirectory(
        prefix="paper2607-midpoint-", dir=directory.parent
    ) as temporary:
        output_dir = Path(temporary)
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = {
                pool.submit(generate_one, generator, args.paper, output_dir, start, end):
                (start, end)
                for start, end in ranges()
            }
            for future in as_completed(futures):
                print(future.result(), flush=True)

        manifest = output_dir / "aggregate_manifest.json"
        subprocess.run(
            [
                sys.executable,
                str(verifier),
                "--directory",
                str(output_dir),
                "--out",
                str(manifest),
            ],
            check=True,
            stdout=subprocess.DEVNULL,
        )
        aggregate = json.loads(manifest.read_text())
        for key, expected in EXPECTED.items():
            if aggregate.get(key) != expected:
                raise RuntimeError(
                    f"aggregate {key}={aggregate.get(key)!r}, expected {expected!r}"
                )

        # Every output has passed decompression, hash, range, count, and aggregate
        # checks. Only now replace the previously certified stream atomically.
        for start, end in ranges():
            name = f"chunk-{start:04d}-{end:04d}.zst"
            os.replace(output_dir / name, directory / name)
            os.replace(
                output_dir / f"{name}.json",
                directory / f"{name}.json",
            )
        os.replace(manifest, directory / "aggregate_manifest.json")
        shutil.rmtree(output_dir, ignore_errors=True)

    print(json.dumps(EXPECTED, sort_keys=True), flush=True)
    print("published 36 verified midpoint-fusion shards", flush=True)


if __name__ == "__main__":
    main()
