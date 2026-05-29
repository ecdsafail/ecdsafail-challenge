#!/usr/bin/env bash
# Ensure Rust + a working C linker are installed. Idempotent: no-op if
# `cargo` and `cc` are already on PATH.
set -euo pipefail

SUDO=""
if [[ ${EUID:-$(id -u)} -ne 0 ]] && command -v sudo >/dev/null 2>&1; then
  SUDO="sudo"
fi

# 1. System deps: a C compiler/linker, plus curl for the rustup bootstrap.
#    `cargo build` shells out to `cc` to link, so this is required even
#    after rustup itself succeeds.
need_cc=0
command -v cc >/dev/null 2>&1 || need_cc=1
command -v curl >/dev/null 2>&1 || need_cc=1

if [[ "${need_cc}" -eq 1 ]]; then
  if command -v apt-get >/dev/null 2>&1; then
    export DEBIAN_FRONTEND=noninteractive
    ${SUDO} apt-get update
    ${SUDO} apt-get install -y --no-install-recommends gcc libc6-dev curl ca-certificates
  elif command -v dnf >/dev/null 2>&1; then
    ${SUDO} dnf install -y gcc glibc-devel curl ca-certificates
  elif command -v yum >/dev/null 2>&1; then
    ${SUDO} yum install -y gcc glibc-devel curl ca-certificates
  elif command -v apk >/dev/null 2>&1; then
    ${SUDO} apk add --no-cache gcc musl-dev curl ca-certificates
  elif command -v pacman >/dev/null 2>&1; then
    ${SUDO} pacman -Sy --noconfirm gcc curl ca-certificates
  elif command -v zypper >/dev/null 2>&1; then
    ${SUDO} zypper --non-interactive install gcc glibc-devel curl ca-certificates
  elif command -v brew >/dev/null 2>&1; then
    : # macOS: cc comes from Xcode CLT, which `xcode-select --install` handles.
  else
    echo "setup.sh: no supported package manager found; install gcc manually" >&2
    exit 1
  fi
fi

# 2. Rust toolchain.
if ! command -v cargo >/dev/null 2>&1; then
  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
    | sh -s -- -y --profile minimal
fi

# shellcheck disable=SC1091
. "$HOME/.cargo/env" 2>/dev/null || true
