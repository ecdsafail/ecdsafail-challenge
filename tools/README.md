# qmap — qubit-allocation profiler + viewer

A small, **opt-in** profiler and interactive HTML viewer for *where qubits go over
circuit time* in the point-add builder. Zero cost when disabled — the instrumentation
is gated entirely behind the `QMAP` env var.

It answers questions like: at any point in the circuit, how many qubits are live, and
which register role are they in (target_x, target_y, the GCD register `u`, the
transcript sidecar, scratch)?

## Contents

- `qmap_explorer.py` — interactive HTML viewer. Pure Python stdlib + a browser; no
  numpy/matplotlib/PIL required.
- `qmap_instrumentation.patch` — env-gated builder hooks. **No behavior change unless
  `QMAP` is set.**
- `README.md` — this file.

## 1. Install the instrumentation

```sh
git apply tools/qmap_instrumentation.patch
cargo build --release --bin build_circuit
```

The patch is against commit `94927be`. The hook sites (`alloc_qubit` / `free` /
`reacquire` / `push_op`, and the register allocations in the dialog GCD) are stable
across revisions, so if it doesn't apply cleanly to a newer tip it re-fits with minor
context adjustment — the changes are: a few struct fields, increment/decrement of a
running per-role counter at alloc/free, and a `set_role(...)` tag wrapped around each
register allocation.

## 2. Capture a profile

Whole circuit, anti-aliased (each snapshot records the **max** within its bucket, so a
coarse stride still captures true peaks), ~1 second:

```sh
QMAP=1 QMAP_STRIDE=500 QMAP_OUT=qmap.tsv ./target/release/build_circuit
```

True per-op resolution inside a span (no aliasing — snapshots every op in the window):

```sh
QMAP=1 QMAP_OP_START=1090000 QMAP_OP_END=1130000 QMAP_OUT=qmap_win.tsv ./target/release/build_circuit
```

Environment variables:

| var | meaning |
|---|---|
| `QMAP=1` | enable profiling (otherwise zero cost) |
| `QMAP_STRIDE=N` | snapshot every N ops (default 15000). Snapshots record the per-bucket max, so peaks survive coarse strides. |
| `QMAP_OP_START` / `QMAP_OP_END` | snapshot **every** op inside `[start, end)` — true 1-op resolution in a span |
| `QMAP_OUT=path` | output TSV path (default `/tmp/qmap.tsv`) |

Snapshots are O(1) (a running per-role count maintained at alloc/free), so even fine
capture is cheap.

## 3. View it

```sh
python3 tools/qmap_explorer.py qmap.tsv qmap.html "my run"
open qmap.html      # or xdg-open / just open the file in a browser
```

Explorer controls:

- **Stacked-by-role ↔ scratch-only** toggle
- **Per-role show/hide** chips (drop the idle I/O registers to magnify the rest, etc.)
- **Minimap** window-slider (drag the window to pan, drag empty space to select a span)
- **Wheel / trackpad zoom**, cursor-centered; **drag to pan**; double-click to reset
- **Y-axis auto-scales** to the currently visible window
- **Hover** for per-snapshot detail (op index, phase, per-role counts)

## TSV format

Tab-separated, one row per snapshot:

```
op_idx  phase  active  scr_res scr_live  tx_res tx_live  ty_res ty_live  u_res u_live  tr_res tr_live
```

`active` is total live qubits; the per-role pairs are `[reserved, live]` counts for
scratch / target_x / target_y / `u` / transcript. The viewer reads any file in this
format, so you can also generate it from your own tooling.
