#!/usr/bin/env bash
# Benchmark a submission.
#
#   1. Wipe stale ops.bin / score.json so a contestant cannot pre-seed them.
#   2. Run build_circuit (UNTRUSTED — runs contestant code in src/point_add)
#      in its own process group; on exit, kill the whole group so a forked
#      child cannot survive and overwrite score.json after step 4.
#   3. Verify ops.bin exists. If build_circuit exited early or crashed,
#      the file is missing and we fail closed.
#   4. Run eval_circuit (TRUSTED — never imports contestant code) which
#      re-simulates the op stream, validates correctness/reversibility/
#      phase, and writes the canonical score.json.
#
# All command-line arguments are forwarded to eval_circuit (e.g. --note ...).
set -euo pipefail

# shellcheck disable=SC1091
. "$HOME/.cargo/env" 2>/dev/null || true

# 1. Clean slate.
rm -f ops.bin score.json

# Make sure both binaries are present (cheap rebuild — no-op if up to date).
RUSTFLAGS='-C linker=gcc' cargo build --release --bin build_circuit --bin eval_circuit

# 2. Run build_circuit in its own process group, then nuke the group.
#    `setsid` puts it in a fresh pgid; `kill -KILL -<pgid>` reaches every
#    child including double-forked daemons. We trap so we always reap.
#    On macOS `setsid` is missing — fall back to `set -m` + `kill %1`.
cleanup_pgid=""
cleanup() {
  if [[ -n "${cleanup_pgid}" ]]; then
    kill -KILL -"${cleanup_pgid}" 2>/dev/null || true
  fi
}
trap cleanup EXIT

if command -v setsid >/dev/null 2>&1; then
  setsid ./target/release/build_circuit &
  build_pid=$!
  cleanup_pgid="${build_pid}"
  set +e
  wait "${build_pid}"
  build_status=$?
  set -e
  kill -KILL -"${cleanup_pgid}" 2>/dev/null || true
  cleanup_pgid=""
else
  # Fallback: bash job control puts the background pipeline in its own pgid.
  set -m
  ./target/release/build_circuit &
  build_pid=$!
  cleanup_pgid="${build_pid}"
  set +e
  wait "${build_pid}"
  build_status=$?
  set -e
  kill -KILL -"${cleanup_pgid}" 2>/dev/null || true
  cleanup_pgid=""
  set +m
fi

if [[ "${build_status}" -ne 0 ]]; then
  echo "!! build_circuit exited with status ${build_status}" >&2
  exit "${build_status}"
fi

# 3. Verify ops.bin actually got produced.
if [[ ! -s ops.bin ]]; then
  echo "!! build_circuit did not produce ops.bin" >&2
  exit 1
fi

# 4. Trusted scoring stage.
./target/release/eval_circuit "$@"
