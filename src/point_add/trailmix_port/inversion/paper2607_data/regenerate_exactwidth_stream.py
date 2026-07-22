#!/usr/bin/env python3
"""Regenerate and certify the Q824 paper2607 primitive stream."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import os
from pathlib import Path
import shutil
import subprocess


SCHEDULE_END = 1616
CHUNK_STEPS = 45
PINNED_SUPPORT_COMMIT = "ac1ecffee14b5a977421b75669c52db6b4033646"
PINNED_SUPPORT_SHA256 = (
    "067d363deeabb6532b52f42eba884b0d184c5b74aa14d2c0d33e5579f668d277"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--python",
        type=Path,
        default=Path("/private/tmp/paper2607-venv/bin/python"),
    )
    parser.add_argument(
        "--support",
        type=Path,
        default=Path("/private/tmp/paper2607-upstream"),
    )
    parser.add_argument("--jobs", type=int, default=8)
    args = parser.parse_args()

    data = Path(__file__).resolve().parent
    generator = data / "generate_eea_blob.py"
    verifier = data.parent / "paper2607_exactwidth_data" / "verify_exactwidth_stream.py"
    support_source = args.support / "eea_circuit_updated.py"

    if args.out.exists():
        raise SystemExit(f"refusing to overwrite existing output: {args.out}")
    if args.jobs < 1:
        raise SystemExit("--jobs must be positive")
    if not args.python.is_file():
        raise SystemExit(f"missing Python interpreter: {args.python}")
    if not support_source.is_file():
        raise SystemExit(f"missing pinned support source: {support_source}")

    commit = subprocess.run(
        ["git", "-C", str(args.support), "rev-parse", "HEAD"],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    ).stdout.strip()
    if commit != PINNED_SUPPORT_COMMIT:
        raise SystemExit(f"support commit {commit}, expected {PINNED_SUPPORT_COMMIT}")
    support_hash = sha256(support_source)
    if support_hash != PINNED_SUPPORT_SHA256:
        raise SystemExit(
            f"support source hash {support_hash}, expected {PINNED_SUPPORT_SHA256}"
        )

    args.out.mkdir(parents=True)
    env = os.environ.copy()
    prior_pythonpath = env.get("PYTHONPATH")
    env["PYTHONPATH"] = os.pathsep.join(
        [str(args.support), str(data)]
        + ([prior_pythonpath] if prior_pythonpath else [])
    )

    ranges = [
        (start, min(start + CHUNK_STEPS - 1, SCHEDULE_END))
        for start in range(1, SCHEDULE_END + 1, CHUNK_STEPS)
    ]

    def generate(bounds: tuple[int, int]) -> tuple[int, int]:
        start, end = bounds
        output = args.out / f"chunk-{start:04d}-{end:04d}.zst"
        result = subprocess.run(
            [
                str(args.python),
                str(generator),
                "--paper",
                str(data),
                "--module",
                "eea_circuit_s835_exactwidth_dirty12",
                "--out",
                str(output),
                "--start",
                str(start),
                "--end",
                str(end),
                "--schedule-end",
                str(SCHEDULE_END),
                "--aux-size",
                "12",
                "--expected-qubits",
                "578",
            ],
            check=True,
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        output.with_suffix(output.suffix + ".log").write_text(
            result.stdout, encoding="utf-8"
        )
        return start, end

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(generate, bounds) for bounds in ranges]
        for future in as_completed(futures):
            start, end = future.result()
            print(f"PASS shard {start:04d}-{end:04d}", flush=True)

    aggregate = args.out / "aggregate.json"
    subprocess.run(
        [
            str(args.python),
            str(verifier),
            "--directory",
            str(args.out),
            "--out",
            str(aggregate),
        ],
        check=True,
        env=env,
    )
    shutil.copyfile(aggregate, args.out / "aggregate_manifest.json")
    print(f"PASS generated and certified {len(ranges)} shards in {args.out}")
    print(f"aggregate_sha256={sha256(aggregate)}")


if __name__ == "__main__":
    main()
