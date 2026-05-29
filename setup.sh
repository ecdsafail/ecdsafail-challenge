#!/usr/bin/env bash
# Ensure Rust + a working C linker are installed. Idempotent: no-op if
# `cargo` and a supported C compiler/linker are already on PATH.
set -euo pipefail

SUDO=""
if [[ ${EUID:-$(id -u)} -ne 0 ]] && command -v sudo >/dev/null 2>&1; then
  SUDO="sudo"
fi

find_c_compiler() {
  if [[ -n "${CC:-}" ]] && command -v "${CC}" >/dev/null 2>&1; then
    command -v "${CC}"
    return 0
  fi

  local candidate
  for candidate in gcc cc clang; do
    if command -v "${candidate}" >/dev/null 2>&1; then
      command -v "${candidate}"
      return 0
    fi
  done

  return 1
}

install_system_deps() {
  if command -v apt-get >/dev/null 2>&1; then
    export DEBIAN_FRONTEND=noninteractive
    ${SUDO} apt-get update
    ${SUDO} apt-get install -y --no-install-recommends gcc libc6-dev ca-certificates
  elif command -v dnf >/dev/null 2>&1; then
    ${SUDO} dnf install -y gcc glibc-devel ca-certificates
  elif command -v yum >/dev/null 2>&1; then
    ${SUDO} yum install -y gcc glibc-devel ca-certificates
  elif command -v apk >/dev/null 2>&1; then
    ${SUDO} apk add --no-cache gcc musl-dev ca-certificates
  elif command -v pacman >/dev/null 2>&1; then
    ${SUDO} pacman -Sy --noconfirm gcc ca-certificates
  elif command -v zypper >/dev/null 2>&1; then
    ${SUDO} zypper --non-interactive install gcc glibc-devel ca-certificates
  elif command -v brew >/dev/null 2>&1; then
    : # macOS: cc comes from Xcode CLT, which `xcode-select --install` handles.
  else
    return 1
  fi
}

# 1. System deps: a C compiler/linker.
#    Cargo needs a linker, and some crates may shell out through the `cc` crate.
if ! find_c_compiler >/dev/null 2>&1; then
  if ! install_system_deps; then
    cat >&2 <<'EOF'
setup.sh: failed to install system dependencies.

This environment needs a C compiler/linker before Rust can build this
repo. If package metadata cannot be downloaded, enable network/DNS for the
sandbox or use an image that already includes gcc (or clang), glibc-devel,
and ca-certificates.
EOF
    exit 1
  fi
fi

# Prefer a compiler that already exists. `benchmark.sh` repeats this detection
# because exported variables from setup.sh do not persist across separate runs.
compiler="$(find_c_compiler || true)"
if [[ -z "${compiler}" ]]; then
  echo "setup.sh: no C compiler found; install gcc or clang" >&2
  exit 1
fi
export CC="${compiler}"

# 2. Rust toolchain.
if ! command -v cargo >/dev/null 2>&1; then
  if ! curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
    | sh -s -- -y --profile minimal; then
    cat >&2 <<'EOF'
setup.sh: failed to install Rust with rustup.

If this sandbox has no outbound network/DNS, use an image that already includes
rustup/cargo and this repo's Rust toolchain.
EOF
    exit 1
  fi
fi

# shellcheck disable=SC1091
. "$HOME/.cargo/env" 2>/dev/null || true

if ! command -v cargo >/dev/null 2>&1; then
  echo "setup.sh: cargo is still not available after setup" >&2
  exit 1
fi
