//! Circuit builder: emits KMX text and manages qubit/bit allocation.
//!
//! Public APIs use typed `Qubit`/`Cbit` wrappers while the internal op
//! encoding stays compact `u32` ids. The builder tracks free pools for
//! recycling and emits REGISTER/APPEND declarations.

use std::fmt::Write;
use std::sync::Arc;

pub use crate::point_add::trailmix_port::tracker::phase_lattice::{AbsVal, PhaseTracker};

/// Build per-qubit 64-shot lanes that broadcast a SINGLE classical value
/// to all 64 shots. **You must call this explicitly when you want all
/// shots to see the same input** — `alloc_input_qreg_bits_with_lanes`
/// won't do it for you. This makes "I am only testing one input"
/// visible at the call site, instead of being a silent default.
///
/// `bytes` is little-endian: bit 0 of bytes[0] = qubit 0, etc.
/// Returns a `Vec<u64>` of length `n_qubits`, each `u64::MAX` (all-shot 1)
/// or 0 (all-shot 0) per the bit pattern.
#[must_use]
pub fn replicate_classical_to_lanes(bytes: &[u8], n_qubits: usize) -> Vec<u64> {
    (0..n_qubits)
        .map(|i| {
            let byte_idx = i / 8;
            let bit_idx = i % 8;
            let bit = if byte_idx < bytes.len() {
                (bytes[byte_idx] >> bit_idx) & 1
            } else {
                0
            };
            if bit == 1 {
                u64::MAX
            } else {
                0
            }
        })
        .collect()
}

/// Enter a scope, capturing source file and line. Returns a u32
/// `push_seq` that must be passed to `exit_scope`. Use `ScopeGuard`
/// for RAII-style pairing.
#[macro_export]
macro_rules! enter_scope {
    ($circ:expr, $name:expr) => {
        $circ.enter_scope_at($name, file!(), line!())
    };
}

/// Compact representation of a single kmx op. Each Op is 16 bytes
/// (aligned). Avoids the ~40–50 bytes per op that `Vec<String>`
/// would use, which dominates memory on 100M+ op circuits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Op {
    Register(u32),
    AppendQubit(u32, u32),
    AppendBit(u32, u32),
    X(u32),
    Z(u32),
    Cx(u32, u32),
    Cz(u32, u32),
    Ccx(u32, u32, u32),
    Ccz(u32, u32, u32),
    Swap(u32, u32),
    Hmr(u32, u32),
    R(u32),
    Neg,
    PushCondition(u32),
    PopCondition,
    BitInvert(u32),
    BitStore0(u32),
    BitStore1(u32),
}

impl Op {
    /// Render this op as one line of kmx text.
    #[must_use]
    pub fn kmx_string(&self) -> String {
        match *self {
            Op::Register(id) => format!("REGISTER r{id}"),
            Op::AppendQubit(q, r) => format!("APPEND_TO_REGISTER q{q} r{r}"),
            Op::AppendBit(b, r) => format!("APPEND_TO_REGISTER b{b} r{r}"),
            Op::X(q) => format!("X q{q}"),
            Op::Z(q) => format!("Z q{q}"),
            Op::Cx(c, t) => format!("CX q{c} q{t}"),
            Op::Cz(a, b) => format!("CZ q{a} q{b}"),
            Op::Ccx(a, b, c) => format!("CCX q{a} q{b} q{c}"),
            Op::Ccz(a, b, c) => format!("CCZ q{a} q{b} q{c}"),
            Op::Swap(a, b) => format!("SWAP q{a} q{b}"),
            Op::Hmr(q, b) => format!("HMR q{q} b{b}"),
            Op::R(q) => format!("R q{q}"),
            Op::Neg => "NEG".to_string(),
            Op::PushCondition(b) => format!("PUSH_CONDITION if b{b}"),
            Op::PopCondition => "POP_CONDITION".to_string(),
            Op::BitInvert(b) => format!("BIT_INVERT b{b}"),
            Op::BitStore0(b) => format!("BIT_STORE0 b{b}"),
            Op::BitStore1(b) => format!("BIT_STORE1 b{b}"),
        }
    }

    /// True if this op is a real gate (not a register declaration).
    #[must_use]
    pub fn is_gate(&self) -> bool {
        !matches!(
            self,
            Op::Register(_) | Op::AppendQubit(_, _) | Op::AppendBit(_, _)
        )
    }

    /// Qubits referenced (read or written) by this op. Empty for
    /// non-gate / non-qubit ops. Used by the strict-dealloc check.
    pub fn touched_qubits(&self) -> impl Iterator<Item = u32> {
        let mut buf: [u32; 3] = [0; 3];
        let n: usize = match *self {
            Op::X(q) | Op::Z(q) | Op::R(q) | Op::Hmr(q, _) => {
                buf[0] = q;
                1
            }
            Op::Cx(a, b) | Op::Cz(a, b) | Op::Swap(a, b) => {
                buf[0] = a;
                buf[1] = b;
                2
            }
            Op::Ccx(a, b, c) | Op::Ccz(a, b, c) => {
                buf[0] = a;
                buf[1] = b;
                buf[2] = c;
                3
            }
            _ => 0,
        };
        buf.into_iter().take(n)
    }
}

/// If `op` is one of the kickmix gates that supports an inline
/// `if b<bit>` suffix, return the collapsed kmx line; otherwise
/// return None. Used by `write_kmx` to fuse `PUSH_CONDITION` +
/// op + `POP_CONDITION` triples into a single kmx instruction.
///
/// HMR, R, BIT_*, REGISTER, `APPEND_TO_REGISTER`, and nested
/// PUSH/POP are NOT conditional-able; they fall through to the
/// 3-instruction form.
fn inline_conditional(op: &Op, bit: u32) -> Option<String> {
    match *op {
        Op::X(q) => Some(format!("X q{q} if b{bit}")),
        Op::Z(q) => Some(format!("Z q{q} if b{bit}")),
        Op::Cx(c, t) => Some(format!("CX q{c} q{t} if b{bit}")),
        Op::Cz(a, b) => Some(format!("CZ q{a} q{b} if b{bit}")),
        Op::Ccx(a, b, c) => Some(format!("CCX q{a} q{b} q{c} if b{bit}")),
        Op::Ccz(a, b, c) => Some(format!("CCZ q{a} q{b} q{c} if b{bit}")),
        Op::Swap(a, b) => Some(format!("SWAP q{a} q{b} if b{bit}")),
        Op::Neg => Some(format!("NEG if b{bit}")),
        _ => None,
    }
}

/// **Module-private** type-safe wrapper around a qubit id. Used ONLY
/// inside `circuit.rs` for the gate-emission machinery (Op enum,
/// sim/tracker internals). External modules cannot import or
/// construct this type — they must go through `QReg`.
#[repr(transparent)]
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash, PartialOrd, Ord)]
struct Qubit(u32);

impl Qubit {
    #[inline]
    fn raw(self) -> u32 {
        self.0
    }
}
impl From<u32> for Qubit {
    #[inline]
    fn from(n: u32) -> Self {
        Qubit(n)
    }
}
impl From<Qubit> for u32 {
    #[inline]
    fn from(q: Qubit) -> Self {
        q.0
    }
}

/// Type-safe wrapper around a classical bit id.
#[repr(transparent)]
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash, PartialOrd, Ord)]
pub struct Cbit(pub u32);

impl Cbit {
    #[inline]
    #[must_use]
    pub fn raw(self) -> u32 {
        self.0
    }
}
impl From<u32> for Cbit {
    #[inline]
    fn from(n: u32) -> Self {
        Cbit(n)
    }
}
impl From<Cbit> for u32 {
    #[inline]
    fn from(b: Cbit) -> Self {
        b.0
    }
}

pub type CReg = [Cbit];

/// Owning, drop-tracked **single-qubit** register. Linear: not `Copy`,
/// not `Clone`. Exactly one binding owns each qubit at any moment;
/// passing it around requires move semantics. When the `QReg` goes
/// out of scope, the qubit is queued onto the Circuit's `pending_frees`
/// and drain at the next gate emission (or via `flush_pending_frees`).
///
/// `QReg` is the ONLY way circuit-construction code can refer to a
/// qubit. The raw `u32` qubit-id stays internal to circuit.rs; no
/// other module can extract or construct it directly.
///
/// For multi-bit registers, use `Vec<QReg>` (via `alloc_qreg_bits`).
#[derive(Debug)]
pub struct QReg {
    /// The internal qubit id. **Module-private** to circuit.rs — no
    /// other module in the crate can read or set this field. The only
    /// way to reference a qubit is via `&QReg` passed to one of
    /// `Circuit`'s gate-emission methods.
    id: u32,
    /// Strong ref to the Circuit's pending-frees queue.
    pending: std::rc::Rc<std::cell::RefCell<Vec<u32>>>,
    /// When set, drop becomes a no-op. **Module-private**: only the
    /// post-`destroy_sim` path sets this on output qubits whose Circuit
    /// has already been consumed. There is no public `detach()` —
    /// every other path must free explicitly via `circ.zero_and_free`.
    detached: bool,
}

impl QReg {
    /// Crate-private accessor for the raw qubit id. Used by sibling
    /// modules that need to look up sim/tracker state for this `QReg`
    /// without going through a full Circuit gate method (e.g. the
    /// `tracker::ghost` machinery, which needs `sim_get_mask` + phase
    /// tracker promotion for a borrowed `QReg`). External callers
    /// have no path to a qubit id.
    #[inline]
    pub(crate) fn id(&self) -> u32 {
        self.id
    }

    /// Consume a public output qubit without scheduling a zero-free.
    ///
    /// The normal `Drop` path treats every still-owned `QReg` as scratch that
    /// must be returned to |0>. Bridge/export code registers output wires and
    /// then detaches them because those wires intentionally carry the public
    /// result into the external evaluator.
    #[inline]
    pub(crate) fn detach(mut self) -> u32 {
        self.detached = true;
        self.id
    }
}

impl Drop for QReg {
    fn drop(&mut self) {
        if self.detached {
            return;
        }
        self.pending.borrow_mut().push(self.id);
    }
}

/// Shared-ownership view of a `QReg`. Wraps `Rc<QReg>` so multiple
/// holders can keep the qubit alive without lifetime parameters in
/// the surrounding type. The underlying `QReg::drop` fires when the
/// LAST `SharedQReg` clone drops, queuing the free against the
/// strict-dealloc check just like a bare `QReg` would. Drop timing
/// is deterministic when all clones live inside one structure that
/// drops together (e.g. a state-machine storage tree).
///
/// Use `&*shared` or `shared.deref()` to get `&QReg` for
/// `Circuit`'s gate-emission methods.
#[derive(Clone, Debug)]
pub struct SharedQReg {
    inner: std::rc::Rc<QReg>,
}

impl SharedQReg {
    /// Wrap an owned `QReg`. After this call, the `QReg`'s drop is
    /// gated by Rc — it fires only when every clone of this
    /// `SharedQReg` has been dropped.
    #[must_use]
    pub fn new(q: QReg) -> Self {
        Self {
            inner: std::rc::Rc::new(q),
        }
    }

    /// Direct read of the underlying `QReg`. Same as `&*shared`,
    /// but explicit for readability.
    #[must_use]
    pub fn as_qreg(&self) -> &QReg {
        &self.inner
    }
}

impl std::ops::Deref for SharedQReg {
    type Target = QReg;
    fn deref(&self) -> &QReg {
        &self.inner
    }
}

/// Either an owned [`QReg`] (drops will queue its qubits for free) or
/// a borrow of one (drop is a no-op). Lets a single function signature
/// serve both consume and non-consume callers: pass `Owned(qreg)` when
/// you want the function to retire the register at its last gate-touch,
/// or `Borrowed(&qreg)` when the caller still owns it. Equivalent to
/// `Cow<'a, QReg>` semantically but with explicit drop behavior.
///
/// `Deref<Target=QReg>` routes all access through the underlying `QReg`
/// — `BorrowedQReg` deliberately does NOT add its own slice / Qubit
/// accessors, so any leakage of bare Qubit IDs goes through the `QReg`
/// API that the rest of the codebase already exposes (and which is a
/// known leak point that the linear-QReg discipline should eventually
/// tighten).
#[derive(Debug)]
pub enum BorrowedQReg<'a> {
    Owned(QReg),
    Borrowed(&'a QReg),
}

impl std::ops::Deref for BorrowedQReg<'_> {
    type Target = QReg;
    fn deref(&self) -> &QReg {
        match self {
            BorrowedQReg::Owned(q) => q,
            BorrowedQReg::Borrowed(q) => q,
        }
    }
}

impl From<QReg> for BorrowedQReg<'_> {
    fn from(q: QReg) -> Self {
        BorrowedQReg::Owned(q)
    }
}
impl<'a> From<&'a QReg> for BorrowedQReg<'a> {
    fn from(q: &'a QReg) -> Self {
        BorrowedQReg::Borrowed(q)
    }
}

/// Sim snapshot captured when streaming truncation clears the op
/// buffer. Lets the debugger replay the currently-retained tail of
/// ops from a correct starting point (the state at `op_idx` =
/// `ops_truncated_at_snapshot`).
#[derive(Clone, Debug)]
pub struct TruncationSnapshot {
    pub(crate) ops_truncated: u64,
    pub(crate) qubits: Vec<u64>,
    pub(crate) bits: Vec<u64>,
    pub(crate) phase: u64,
    pub(crate) cond_stack: Vec<u32>,
    pub(crate) hmr_counter: u64,
    pub(crate) r_on_nonzero_events: u64,
}

/// Read-only simulator state after the live simulator has been
/// destroyed. This is the only sanctioned way for non-debugger code
/// to inspect quantum register values.
#[derive(Clone, Debug)]
pub struct DestroyedSimState {
    qubits: Vec<u64>,
    bits: Vec<u64>,
    phase: u64,
    phase_errors: u32,
    phase_errors_by_tag: std::collections::HashMap<String, u32>,
}

impl std::ops::Index<usize> for DestroyedSimState {
    type Output = u64;

    fn index(&self, index: usize) -> &Self::Output {
        &self.qubits[index]
    }
}

impl DestroyedSimState {
    #[must_use]
    pub fn qubit_mask(&self, q: &QReg) -> u64 {
        self.qubits.get(q.id as usize).copied().unwrap_or(0)
    }

    #[must_use]
    pub fn read_bit(&self, q: &QReg) -> u8 {
        self.read_bit_shot(q, 0)
    }

    #[must_use]
    pub fn read_bit_shot(&self, q: &QReg, shot: usize) -> u8 {
        ((self.qubit_mask(q) >> shot) & 1) as u8
    }

    #[must_use]
    pub fn read_bytes(&self, qs: &[QReg]) -> Vec<u8> {
        self.read_bytes_shot(qs, 0)
    }

    #[must_use]
    pub fn read_bytes_shot(&self, qs: &[QReg], shot: usize) -> Vec<u8> {
        let nbytes = qs.len().div_ceil(8);
        let mut bytes = vec![0u8; nbytes];
        for (i, q) in qs.iter().enumerate() {
            if ((self.qubit_mask(q) >> shot) & 1) == 1 {
                bytes[i / 8] |= 1 << (i % 8);
            }
        }
        bytes
    }

    #[must_use]
    pub fn is_zero(&self, qs: &[QReg]) -> bool {
        qs.iter().all(|q| self.qubit_mask(q) == 0)
    }

    #[must_use]
    pub fn bit_len(&self, qs: &[QReg]) -> usize {
        for i in (0..qs.len()).rev() {
            if self.read_bit(&qs[i]) == 1 {
                return i + 1;
            }
        }
        0
    }

    #[must_use]
    pub fn bit_mask(&self, id: u32) -> u64 {
        self.bits.get(id as usize).copied().unwrap_or(0)
    }

    #[must_use]
    pub fn phase_mask(&self) -> u64 {
        self.phase
    }

    #[must_use]
    pub fn phase_error_count(&self) -> u32 {
        self.phase_errors
    }

    #[must_use]
    pub fn phase_error_counts_by_tag(&self) -> &std::collections::HashMap<String, u32> {
        &self.phase_errors_by_tag
    }
}

/// Restricted read-only simulator view used by contract checks.
/// This intentionally exposes only shot-local reads and does not
/// provide access to broader Circuit state.
#[derive(Clone, Copy)]
pub(crate) struct ContractSimView<'a> {
    qubits: &'a [u64],
    bits: &'a [u64],
    /// Pending (unresolved) spooky-pebble ghosts at the point the
    /// contract runs. Used by ghost-set contract assertions (count and
    /// per-ghost tape-bit masks). Only the test suite exercises these.
    #[cfg(test)]
    ghosts: &'a [crate::point_add::trailmix_port::tracker::ghost::GhostRecord],
}

impl ContractSimView<'_> {
    fn qubit_mask_raw(&self, id: u32) -> u64 {
        self.qubits.get(id as usize).copied().unwrap_or(0)
    }

    /// Number of unresolved ghosts at this point.
    #[cfg(test)]
    pub fn pending_ghost_count(&self) -> usize {
        self.ghosts.len()
    }

    /// The pending-ghost record for `id`, if present.
    #[cfg(test)]
    pub fn ghost(&self, id: u64) -> Option<&crate::point_add::trailmix_port::tracker::ghost::GhostRecord> {
        self.ghosts.iter().find(|g| g.id == id)
    }

    /// The tape bit ghost `id` stands for, on `shot` (its HMR'd
    /// qubit's value at creation). `None` if no such pending ghost.
    #[cfg(test)]
    pub fn ghost_value_shot(&self, id: u64, shot: usize) -> Option<bool> {
        self.ghost(id).map(|g| (g.mask_at_hmr >> shot) & 1 == 1)
    }

    #[allow(dead_code)]
    pub(crate) fn bit_mask(&self, id: u32) -> u64 {
        self.bits.get(id as usize).copied().unwrap_or(0)
    }

    pub fn read_u256_shot(&self, reg: &[QReg], shot: usize) -> num_bigint::BigUint {
        let mut v = num_bigint::BigUint::from(0u32);
        for (i, q) in reg.iter().enumerate() {
            if (self.qubit_mask_raw(q.id) >> shot) & 1 == 1 {
                v |= num_bigint::BigUint::from(1u32) << i;
            }
        }
        v
    }

    pub fn read_bit_shot(&self, q: &QReg, shot: usize) -> bool {
        (self.qubit_mask_raw(q.id) >> shot) & 1 == 1
    }

    pub fn contract_read_u256_shot(&self, reg: &[QReg], shot: usize) -> num_bigint::BigUint {
        self.read_u256_shot(reg, shot)
    }

    pub fn contract_read_bit_shot(&self, q: &QReg, shot: usize) -> bool {
        self.read_bit_shot(q, shot)
    }
}

/// Opaque handle to a captured contract. Methods `update` and `check`
/// take a `&mut Circuit` / `&Circuit` to re-run the closure against the
/// current sim state. Drop automatically removes the entry from the
/// deferred-contracts stack (via the shared Rc<`RefCell`<DeferredStack>>).
pub struct Capture<T: 'static> {
    id: u64,
    label: String,
    stack: std::rc::Rc<std::cell::RefCell<Vec<DeferredContract>>>,
    _marker: std::marker::PhantomData<T>,
}

impl<T: 'static> Capture<T> {
    /// Read-only verify: per-shot, f returns Err(msg) on violation.
    pub(crate) fn check<F>(&self, circ: &Circuit, f: F)
    where
        F: for<'a> FnMut(&T, ContractSimView<'a>, usize) -> Result<(), String>,
    {
        circ.capture_check_impl::<T, F>(self, f);
    }
}

impl<T: 'static> Drop for Capture<T> {
    fn drop(&mut self) {
        // Pop the entry from the shared stack. If contracts were disabled
        // at capture time, no entry exists — silently no-op.
        if let Ok(mut stack) = self.stack.try_borrow_mut() {
            if let Some(idx) = stack.iter().position(|d| d.id == self.id) {
                stack.swap_remove(idx);
            }
        }
    }
}

pub(crate) struct DeferredContract {
    id: u64,
    label: String,
    /// Type-erased Vec<T> where T is the per-shot capture type chosen
    /// by the `contract_capture` caller. Downcast at pop time when the
    /// caller's post closure tells us what T it expects.
    ///
    /// We only store DATA (Vec<T> for some 'static T), not closures.
    /// The pre and post closures run inline (synchronously, in the
    /// caller's stack frame) and only the captured T values cross the
    /// API boundary — so closure captures keep their natural caller
    /// lifetime without any `'static` lie.
    captured: Box<dyn std::any::Any>,
}

pub trait ContractReadable {
    fn contract_read_u256_shot(&self, reg: &[QReg], shot: usize) -> num_bigint::BigUint;
    fn contract_read_bit_shot(&self, q: &QReg, shot: usize) -> bool;
}

impl ContractReadable for ContractSimView<'_> {
    fn contract_read_u256_shot(&self, reg: &[QReg], shot: usize) -> num_bigint::BigUint {
        self.read_u256_shot(reg, shot)
    }

    fn contract_read_bit_shot(&self, q: &QReg, shot: usize) -> bool {
        self.read_bit_shot(q, shot)
    }
}

/// Tracks a logical register of qubits from alloc to free, so
/// downstream analysis (liveness, idle stretches, peak overlap)
/// can surface qubit-packing opportunities.
#[derive(Clone, Debug)]
pub struct TaggedRegion {
    pub name: String,
    pub qubits: Vec<u32>,
    pub start_op: usize,
    pub end_op: Option<usize>,
}

/// Captured per-emit state for an involutary gate, used to undo it
/// when the redundant-op auto-eliminator cancels it against an
/// inverse. Stored in `Circuit::elide_deltas` keyed by `op_idx`.
#[derive(Clone, Debug)]
struct ElideDelta {
    /// Per-operand `(qid, prev_last_touched_op, prev_last_op_on_qubit)`
    /// captured before this op's emission overwrote them. Up to 3
    /// operands per op (CCX/CCZ); slots 1-2 are unused for X/Z and
    /// slot 2 is unused for CX/CZ. Restored on elide so each operand's
    /// `last_touched_op` points back at whatever touched it before the
    /// elided pair, preserving strict-dealloc retention invariants.
    operand_prev: [Option<(u32, Option<u64>, Option<Op>)>; 3],
    /// Section that was current when the op was emitted; used to
    /// undo the increment to `executed_ops_by_section`.
    section: String,
}

pub struct Circuit {
    /// `None` slots are elided no-ops left behind by the redundant-op
    /// auto-eliminator (when an involutary gate cancels with its
    /// inverse). Slot indices are preserved on elide so absolute
    /// `op_idx` (= `ops_truncated` + slot) remains a stable reference for
    /// `last_touched_op`, `elide_deltas`, debugger replay, and the
    /// strict-dealloc retention check.
    ///
    /// Pre-allocated to `ops_cap` capacity at construction (and on
    /// `set_ops_cap`) so the streaming truncation kicks in at exactly
    /// `ops_cap` slots, before `Vec`'s growth strategy would double
    /// past it — which at 100M-op caps would require a 4GB transient
    /// allocation that overflows the 6GB process limit.
    pub ops: Vec<Option<Op>>,
    /// Streaming mode (set via `CIRC_OPS_CAP` env var): number of ops
    /// that have been truncated out of `self.ops` via periodic
    /// `clear()`. `total_ops() = ops_truncated + ops.len()`. Enables
    /// full-scale probes (e.g. EC add at secp) to fit in memory by
    /// dropping op history — downstream replay/debugger/kmx-emit is
    /// invalidated, but `peak_qubits` and op count stay correct.
    pub ops_truncated: u64,
    /// Bounded retained op-tail size for streaming/debugger replay.
    /// Defaults to `DEFAULT_OPS_CAP`, but can be overridden by the
    /// `CIRC_OPS_CAP` environment variable.
    ops_cap: usize,
    /// Persistent total CCX emission count; incremented on every
    /// `Op::Ccx(_)` push regardless of streaming truncation. Unlike
    /// `ccx_count()` which scans the (possibly-cleared) buffer, this
    /// counter is stable across full-scale probes.
    pub ccx_emitted: u64,
    /// Same for CCZ.
    pub ccz_emitted: u64,
    max_qubits_assert: Option<u32>,
    max_ops_assert: Option<u64>,
    /// Hard cap on peak live qubits. Unlike `max_qubits_assert`
    /// (which is purely opt-in via `CIRC_ASSERT_MAX_QUBITS`), this
    /// cap has a built-in default (`DEFAULT_MAX_QUBIT_PEAK = 2000`), a
    /// catastrophe backstop at the hard bound (per-test caps do the real
    /// regression detection). Override via the
    /// `CIRC_ASSERT_MAX_QUBIT_PEAK` env var; set it to 0 to disable.
    /// Fires as soon as `live` exceeds the cap (not at end-of-circuit).
    max_qubit_peak_assert: Option<u32>,
    /// `op_idx` of last gate that touched (read or wrote) each qubit.
    /// Used by `free_qubit`'s strict-dealloc check. `None` means the
    /// qubit was allocated but never touched by any gate — freeing
    /// such a qubit panics, since allocating a qubit and never
    /// writing it is always a bug (it represents wasted ancilla
    /// budget; remove the alloc instead).
    last_touched_op: Vec<Option<u64>>,
    /// Last gate Op that touched each qubit. Used by the redundant-op
    /// auto-eliminator in `push_gate_op` to detect involutary gates
    /// (X, Z, CX, CZ, CCX, CCZ) whose effect cancels — when all the
    /// operands of a new gate show the SAME prior gate as their
    /// `last_op`, the new emission and that prior emission are
    /// algebraic inverses and can be removed as a no-op pair.
    last_op_on_qubit: Vec<Option<Op>>,
    /// Per-emit undo info for involutary gates that are still candidates
    /// for cancellation (i.e. they are the `last_op_on_qubit` for every
    /// one of their operands). Indexed by `op_idx`. Pruned aggressively
    /// once any operand is touched by a different op (the entry can no
    /// longer be elided, since elide requires all operands of the new
    /// op to point at the same prior op).
    elide_deltas: std::collections::HashMap<u64, ElideDelta>,
    /// `op_idx` at the moment of the most recent `alloc_qubit`. Used as
    /// the cutoff in `free_qubit`'s wasteful-retention check: a qubit
    /// freed without having been touched since the last allocation
    /// could have been freed BEFORE that allocation, letting the
    /// allocator reuse its slot. We track the alloc edge rather
    /// than requiring strict gap=0 because a qubit freed N ops
    /// later with NO intervening allocation hasn't actually wasted
    /// anything — its slot wouldn't have been reused regardless.
    last_alloc_op_idx: u64,
    /// Queue of qubits whose owning [`QReg`] has been dropped but
    /// whose `free_qubit` call has not yet fired. Drained at the
    /// next gate emission (so the free's strict-dealloc check fires
    /// at the same `current_op_idx` as the gate that would have
    /// followed the register's last touch). Shared via Rc with all
    /// live `QRegs` (strong on the `QReg` side — outliving the Circuit
    /// is a hard error).
    pending_frees: std::rc::Rc<std::cell::RefCell<Vec<u32>>>,
    next_qubit: u32,
    next_bit: u32,
    // Lowest-id-first free pool: `alloc_qubit` hands out the smallest
    // free id, so freed slots are reused densely from the bottom. This
    // makes `defragment` able to compact live qubits into a contiguous
    // low-id block (the fuzzer's register convention needs input and
    // output to share qubit ids).
    free_qubits: std::collections::BTreeSet<u32>,
    free_bits: Vec<u32>,
    pub peak_qubits: u32,
    pub peak_bits: u32,
    pub peak_at_op: usize,
    pub peak_section: String,
    /// Diagnostic: per-section LOCAL peak occupancy = max total live qubits
    /// observed while that exact section path was current (captured at each
    /// alloc, where live rises). `global_peak - section_peak[s]` is the
    /// headroom available to that section's ops for measurement-venting
    /// without raising the global peak. Gate-neutral (profiling only).
    pub section_peak: std::collections::HashMap<String, u32>,
    /// Snapshot of live qubits at peak. When `peak_qubits` is updated
    /// (grows), we record the live-qubit set (all allocated minus all
    /// currently-freed). Post-construction tools use this to break
    /// down the peak by section tag.
    pub peak_live_qubits: Vec<u32>,
    /// Snapshot of the corresponding qubit tags at peak time. This is
    /// captured eagerly because qubits may be reused later with
    /// different tags.
    pub peak_live_tags: Vec<String>,
    pub live_series: Vec<(usize, u32)>,
    /// Classical simulation: 64 shots in parallel (matching the
    /// zenodo Simulator). Each qubit is a u64 bitmask — bit i is
    /// the qubit's value in shot i. HMR picks an independent u64 of
    /// random bits so all 64 scenarios get different random choices
    /// in one run. Phase cleanness requires `sim_phase == 0` (all
    /// 64 bits zero) at the end.
    pub(crate) sim: Option<Vec<u64>>,
    /// Classical bit simulation, 64 shots. bit i of `sim_bits`[b] =
    /// bit b's value in shot i.
    pub(crate) sim_bits: Option<Vec<u64>>,
    /// Condition stack (bit IDs). Effective condition mask at any
    /// point is AND of `sim_bits`[b] for each stacked b; empty stack
    /// = `u64::MAX` (unconditional). Used to implement `push_condition`
    /// / `pop_condition` as gate-level classical gating.
    pub(crate) sim_condition_stack: Vec<u32>,
    /// Count of R-on-nonzero events. Each R with qubit != 0 in some
    /// shot contributes phase kick in those shots (matching zenodo
    /// semantics). This counter is for reporting only.
    pub(crate) sim_phase_errors: u32,
    /// Per-section phase error counters.
    sim_phase_errors_by_tag: std::collections::HashMap<String, u32>,
    /// Current section label for phase-error attribution.
    pub current_section: String,
    /// Per-shot global phase as a u64 bitmask. Bit i of `sim_phase` is
    /// the phase of shot i (0 = +1, 1 = -1). `sim_phase` == 0 means
    /// all 64 shots are phase-clean.
    pub(crate) sim_phase: u64,
    /// HMR counter used to derive 64-bit random words for HMR and R ops.
    /// One u64 per op, shared across shots.
    pub(crate) sim_hmr_counter: u64,
    /// Seed for HMR PRNG. Drawn fresh per Circuit from system entropy so
    /// every run exercises a different HMR-bit pattern. The fuzz contract
    /// is that the circuit must be value-correct AND phase-clean for ALL
    /// HMR-bit patterns, so pinning the seed would mask real bugs.
    sim_hmr_seed: u64,
    /// Abstract-interpretation F2 polynomial tracker. Tracks every
    /// qubit/bit's symbolic value as an F2 polynomial over input
    /// atoms + fresh HMR random-bit atoms. At finalize, the global
    /// phase polynomial AND every unfreed HMR obligation must be
    /// the zero polynomial. See `phase_lattice.rs`.
    pub phase: PhaseTracker,
    /// (`op_idx`, `section_name`) pairs appended on each `set_section`
    /// call. Lets the debugger walk section boundaries and lets
    /// `goto @section_name` work.
    pub section_marks: Vec<(usize, String)>,
    /// Debugger dump hooks. Keyed by name; value is a closure that
    /// receives the Debugger and arg strings, returns text to print.
    /// Register via `Circuit::register_debug_dump(name, fn)`. Invoke
    /// from debugger REPL with `dump <name> [args]`.
    #[doc(hidden)]
    pub dump_hooks: std::collections::HashMap<String, crate::point_add::trailmix_port::tracker::debugger::DumpHook>,
    /// Snapshot of (`sim_qubits`, `sim_bits`) taken lazily at the first
    /// gate emission — captures state AFTER any `sim_load_reg_bytes`
    /// calls but BEFORE gate execution. This is the debugger's
    /// starting state for replay (cursor = 0).
    pub(crate) initial_sim_state: Option<(Vec<u64>, Vec<u64>)>,
    /// When streaming (`CIRC_OPS_CAP`) truncates the op buffer, we
    /// snapshot the full sim state just before the clear. The
    /// debugger uses this as its starting point so it can still
    /// replay the RETAINED tail of ops (the most recent N since the
    /// last truncation). Contains: (`ops_truncated_at_snapshot`,
    /// `sim_qubits`, `sim_bits`, `sim_phase`, `sim_condition_stack`,
    /// `sim_hmr_counter`, `sim_phase_errors`). If `None`, the debugger
    /// falls back to `initial_sim_state` (which only works when
    /// nothing was truncated).
    pub(crate) truncation_snapshot: Option<TruncationSnapshot>,
    /// Live debugger checkpoints, snapshotted by `push_gate_op` every
    /// `checkpoint_interval` ops. Each entry is the full sim state
    /// AFTER the first `op_idx` ops (where `op_idx` = `ops_truncated` +
    /// local index = absolute global op number). Lets `Debugger::attach`
    /// skip its O(N×M) replay on attach and reuse Circuit's already-run
    /// simulation. The same shape as `TruncationSnapshot` is used so the
    /// debugger consumes both uniformly.
    pub(crate) live_checkpoints: Vec<TruncationSnapshot>,
    /// Op-index spacing between live checkpoints. 10M chosen to balance
    /// memory cost (≈ 25 snapshots × 64 KB sim ≈ 1.6 MB at 250M ops)
    /// against `restore_snapshot_and_replay` worst-case forward steps
    /// (≤ 10M, ~1s on a release build).
    pub(crate) checkpoint_interval: usize,
    /// HMR seed at attach time — replicated here for the debugger
    /// since the field is private.
    pub initial_hmr_seed: u64,
    /// Opt-in register tags for liveness analysis. Populated via
    /// `tag_region` / `untag_region` at build time; consumed by
    /// debugger/profile tooling to report idle stretches and
    /// packing opportunities.
    pub regions: Vec<TaggedRegion>,
    /// Sum over all CCX/CCZ ops of `popcount(fire_mask)`. `fire_mask` is
    /// cond & c1 & c2 (for CCX) or cond & a & b & c (for CCZ). The
    /// per-shot avg executed Toffoli count is this / 64. Matches
    /// Google's "average executed Toffoli count across test cases"
    /// accounting (kickmix 2025 appendix).
    pub executed_toffoli_shots: u64,
    /// Per-section version of `executed_toffoli_shots`. Key = full
    /// section path at the time the op executed.
    pub executed_toffoli_by_section: std::collections::HashMap<String, u64>,
    /// Per-section total op count across the full streamed run.
    pub executed_ops_by_section: std::collections::HashMap<String, u64>,
    /// Explicit scope stack: (name, file, line) pushed by `enter_scope!`
    /// and popped by `exit_scope`. Each op records the depth of this
    /// stack at emission time so the debugger can map `op_idx` → source
    /// site without backtrace capture.
    pub scope_stack: Vec<ScopeFrame>,
    /// Start op index for each scope frame ever pushed. Indexed by the
    /// frame's `push_seq` field. Used by the debugger to print the
    /// scope chain for any `op_idx`.
    pub scope_frames_log: Vec<ScopeFrameLog>,
    /// For each op index, the `push_seq` of the innermost scope alive
    /// at that op. Use `scope_frames_log[op_scope[idx]]` to get the
    /// frame, then walk `parent` up to root. Length == `ops.len()`.
    pub op_scope: Vec<u32>,
    /// Per-qubit-ID name: auto-populated on alloc with
    /// `{current_section}/q{section_alloc_counter}`. Overridable via
    /// `Circuit::name_qubit(q, name)`. Recycled qubit IDs get a fresh
    /// tag on each alloc, so reading the tag shows the MOST RECENT
    /// purpose of that qubit slot.
    pub qubit_tags: Vec<Option<Arc<str>>>,
    /// History log of qubit (re)allocations, so the debugger can
    /// look up "what was q5388 at op 6294421?" — the tag when the
    /// qubit was allocated during the time window covering `op_idx`.
    pub qubit_alloc_log: Vec<QubitAllocEvent>,
    /// History log of classical-bit (re)allocations, parallel to
    /// `qubit_alloc_log`. Anonymous `alloc_bit` / `alloc_bits` calls
    /// record a `b{N}` tag so the debugger can still resolve raw
    /// `bN` references; the named variants
    /// `alloc_bit_named` / `alloc_bits_named` record a user-supplied
    /// label (with `name[i]` suffixing for the batch helper).
    pub bit_alloc_log: Vec<BitAllocEvent>,
    /// History log of spooky-pebble ghost create/resolve events, so the
    /// debugger can reconstruct the pending-ghost set at any cursor.
    /// Parallel to `qubit_alloc_log`. See `GhostEvent`.
    pub ghost_event_log: Vec<GhostEvent>,
    /// Interner for qubit tags. Maps each unique tag string to a single
    /// shared `Arc<str>`, so `qubit_tags`, `qubit_alloc_log`, etc. can
    /// share heap allocations instead of duplicating 30+ bytes per alloc
    /// event. At `ec_add` scale (~24M events, ~thousands of unique tags)
    /// this is the difference between fitting in 8G and OOM.
    tag_intern: std::collections::HashMap<String, Arc<str>>,
    /// Counter for `q{N}` suffixes within the current section. Resets
    /// in `set_section`. Used only for auto-generated tags.
    qubit_tag_counter: usize,
    deferred_contracts: std::rc::Rc<std::cell::RefCell<Vec<DeferredContract>>>,
    /// Monotonic id counter for `Capture<T>` handles.
    capture_next_id: u64,
    /// Interns `#[track_caller]`-captured gate-emission call sites into
    /// `scope_frames_log`. Key = (file, line, parent manual scope seq);
    /// value = the synthesized frame's `push_seq`. Bounded by the number
    /// of distinct gate call sites (a few thousand at `ec_add` scale) —
    /// NOT one entry per op, so `scope_frames_log` stays tiny while
    /// `op_scope` (one u32/op, already present) does the per-op work.
    /// This is what makes the debugger `src <op>` command resolve a
    /// <file:line> in `--release` without any manual `enter_scope!` tags.
    call_site_intern: std::collections::HashMap<(&'static str, u32, u32), u32>,
    /// Unresolved Ghost receipts from `hmr_ghost`. Each row mirrors
    /// the diagnostic fields of the live [`crate::point_add::trailmix_port::tracker::ghost::Ghost`]
    /// so the Circuit-level Drop and `assert_phase_clean` can produce
    /// a useful error even if the Ghost was leaked via
    /// `std::mem::forget`. Drained in `resolve_ghost`.
    pub(crate) pending_ghosts: Vec<crate::point_add::trailmix_port::tracker::ghost::GhostRecord>,
    /// Monotonic id counter for [`crate::point_add::trailmix_port::tracker::ghost::Ghost`].
    /// Each `hmr_ghost` allocates a fresh id; `resolve_ghost` matches
    /// on it when removing the row from `pending_ghosts`.
    pub(crate) next_ghost_id: u64,
}

#[derive(Clone, Debug)]
pub struct QubitAllocEvent {
    /// Op index at which this alloc fired (= `ops.len()` at the time).
    pub op_idx: usize,
    pub qubit: u32,
    /// Tag as of this alloc (scope-prefixed auto-name or caller override).
    ///
    /// Uses `Arc<str>` so repeated allocations of identically-tagged
    /// scratch qubits (e.g. "`cmp_ctrl_tmp`", "carry") share a single heap
    /// allocation. For ec_add-scale emissions with ~24M alloc events and
    /// only ~thousands of unique tag strings, this is a ~4x reduction
    /// in `qubit_alloc_log` memory (see `Circuit::intern_tag`).
    pub tag: Arc<str>,
}

#[derive(Clone, Debug)]
pub struct BitAllocEvent {
    /// Op index at which this alloc fired (= `ops.len()` at the time).
    pub op_idx: usize,
    pub bit: u32,
    /// Tag as of this alloc. Anonymous allocs record `b{N}`; named
    /// allocs record the user label (with `[i]` suffix for batches).
    /// Interned via the shared `tag_intern` map to share heap
    /// storage with the qubit alloc log.
    pub tag: Arc<str>,
}

/// One spooky-pebble ghost lifecycle event, logged for the time-travel
/// debugger so the pending-ghost set can be reconstructed at any cursor
/// position (replay events with `op_idx <= cursor`, create adds, resolve
/// removes). Mirrors the alloc-log pattern.
#[derive(Clone, Debug)]
pub struct GhostEvent {
    /// Global op index at which the event fired.
    pub op_idx: usize,
    /// true = `hmr_ghost` (create); false = `resolve_ghost` (resolve).
    pub create: bool,
    pub id: u64,
    pub anchor_id: u64,
    /// 64-shot tape-bit mask (create only; 0 for resolve).
    pub mask_at_hmr: u64,
    pub bit_raw: u32,
    /// Section path at the event.
    pub section: String,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct CircuitDebugStats {
    pub section_marks: usize,
    pub section_mark_bytes: usize,
    pub executed_ops_sections: usize,
    pub executed_tof_sections: usize,
    pub executed_section_key_bytes: usize,
    pub regions: usize,
    pub scope_frames_log: usize,
    pub op_scope: usize,
    pub qubit_tags_len: usize,
    pub qubit_tags_live: usize,
    pub qubit_alloc_log: usize,
    pub qubit_alloc_tag_bytes: usize,
}

#[derive(Clone, Debug)]
pub struct ScopeFrame {
    pub name: String,
    pub file: &'static str,
    pub line: u32,
    /// Sequence number of this push, matches index into `scope_frames_log`.
    pub push_seq: u32,
}

#[derive(Clone, Debug)]
pub struct ScopeFrameLog {
    pub name: String,
    pub file: &'static str,
    pub line: u32,
    pub parent: Option<u32>,
    pub start_op: usize,
    pub end_op: Option<usize>,
}

impl Default for Circuit {
    fn default() -> Self {
        Self::new()
    }
}

/// Install a process-level address-space cap (`RLIMIT_AS`) so the Rust
/// allocator returns failure (panic with backtrace at the allocation
/// site) before the script's 8 GiB systemd cap SIGKILLs us. Default
/// 6 GiB; override via `CIRC_MEM_RLIMIT_MB=<n>` (or `=0` to disable).
///
/// The crate's sole `unsafe` block: a libc `setrlimit` FFI call to set
/// `RLIMIT_AS` (no safe wrapper in std).
#[allow(unsafe_code)]
fn install_memory_rlimit_once() {
    use std::sync::Once;
    static INIT: Once = Once::new();
    INIT.call_once(|| {
        let limit_mb: u64 = std::env::var("CIRC_MEM_RLIMIT_MB")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(6000);
        if limit_mb == 0 {
            return;
        }
        let limit_bytes: u64 = limit_mb.saturating_mul(1024 * 1024);
        // Linux RLIMIT_AS = 9. rlim_t is u64 on 64-bit Linux.
        #[repr(C)]
        struct RLimit {
            rlim_cur: u64,
            rlim_max: u64,
        }
        const RLIMIT_AS: i32 = 9;
        unsafe extern "C" {
            fn setrlimit(resource: i32, rlim: *const RLimit) -> i32;
            fn getrlimit(resource: i32, rlim: *mut RLimit) -> i32;
        }
        // SAFETY: getrlimit/setrlimit are POSIX libc calls with stable
        // ABI on Linux (rlim_t == u64 on aarch64). Approved C FFI per
        // the no-unsafe-Rust feedback rule.
        unsafe {
            let mut existing = RLimit {
                rlim_cur: 0,
                rlim_max: 0,
            };
            let _ = getrlimit(RLIMIT_AS, &raw mut existing);
            // Only TIGHTEN the cap — never loosen an existing stricter limit.
            let new_cur = if existing.rlim_cur == 0 || limit_bytes < existing.rlim_cur {
                limit_bytes
            } else {
                existing.rlim_cur
            };
            let new_max = if existing.rlim_max == 0 || new_cur < existing.rlim_max {
                new_cur
            } else {
                existing.rlim_max
            };
            let rl = RLimit {
                rlim_cur: new_cur,
                rlim_max: new_max,
            };
            let _ = setrlimit(RLIMIT_AS, &raw const rl);
        }
    });
}

impl Circuit {
    const DEFAULT_OPS_CAP: usize = 100_000_000;
    const MIN_OPS_CAP: usize = 10_000;

    fn configured_ops_cap() -> usize {
        std::env::var("CIRC_OPS_CAP")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .map_or(Self::DEFAULT_OPS_CAP, |n| n.max(Self::MIN_OPS_CAP))
    }

    /// Override the streaming op-buffer cap for this Circuit instance.
    /// Caps below `MIN_OPS_CAP` are clamped up. Use in tests to keep
    /// the buffer small (the env-var route is per-process, this is
    /// per-Circuit and avoids the unsafe `std::env::set_var` API).
    pub fn set_ops_cap(&mut self, cap: usize) {
        self.ops_cap = cap.max(Self::MIN_OPS_CAP);
        // Resize the ops buffer's capacity to match — ensures the
        // streaming truncation kicks in at the configured cap before
        // Vec's doubling strategy would force a transient that
        // overflows the process memory limit.
        self.ops
            .reserve_exact(self.ops_cap.saturating_sub(self.ops.len()));
    }

    fn configured_max_qubits_assert() -> Option<u32> {
        std::env::var("CIRC_ASSERT_MAX_QUBITS")
            .ok()
            .and_then(|s| s.parse::<u32>().ok())
            .filter(|&n| n > 0)
    }

    fn configured_max_ops_assert() -> Option<u64> {
        std::env::var("CIRC_ASSERT_MAX_OPS")
            .ok()
            .and_then(|s| s.parse::<u64>().ok())
            .filter(|&n| n > 0)
    }

    /// Default cap for peak live qubits when `CIRC_ASSERT_MAX_QUBIT_PEAK`
    /// is unset. This is a CATASTROPHE BACKSTOP, not a cost target: it sits
    /// at the CLAUDE.md hard upper bound (2000) so a runaway alloc / infinite
    /// loop trips before it OOMs the 8G wrapper. Real regression detection is
    /// the job of the per-test `set_max_qubit_peak(...)` calls (e.g. EC-add
    /// ~1191, mod_mul_eea 1180, gcd_pack 870), which set tight, documented
    /// ceilings for cost-sensitive circuits. Set via env
    /// `CIRC_ASSERT_MAX_QUBIT_PEAK` to override, or 0 to disable.
    const DEFAULT_MAX_QUBIT_PEAK: u32 = 2000;

    fn configured_max_qubit_peak_assert() -> Option<u32> {
        match std::env::var("CIRC_ASSERT_MAX_QUBIT_PEAK") {
            Ok(s) => s.parse::<u32>().ok().filter(|&n| n > 0),
            Err(_) => Some(Self::DEFAULT_MAX_QUBIT_PEAK),
        }
    }

    /// Declare this circuit's expected peak-qubit ceiling explicitly (a per-test
    /// budget). Use in tests so each asserts its OWN known peak instead of relying
    /// on the global default — a regression that raises the peak fails that test.
    /// The `CIRC_ASSERT_MAX_QUBIT_PEAK` env var still takes precedence (so profile
    /// / diagnostic runs can relax or disable it); when the env is unset this value
    /// replaces the default.
    pub fn set_max_qubit_peak(&mut self, limit: u32) {
        if std::env::var("CIRC_ASSERT_MAX_QUBIT_PEAK").is_err() {
            self.max_qubit_peak_assert = Some(limit);
        }
    }

    #[must_use]
    pub fn new() -> Self {
        install_memory_rlimit_once();
        let ops_cap = Self::configured_ops_cap();
        // Fresh per-Circuit HMR seed from thread_rng. The fuzz contract
        // requires correctness for ALL HMR-bit patterns, so a fresh seed
        // every run is the right default.
        let hmr_seed: u64 = rand::Rng::gen(&mut rand::thread_rng());
        Circuit {
            // Pre-allocate the ops buffer to ops_cap so the streaming
            // truncation kicks in BEFORE Vec's doubling strategy
            // would force a transient that overflows the process
            // memory limit.
            ops: Vec::with_capacity(ops_cap),
            ops_truncated: 0,
            ops_cap,
            ccx_emitted: 0,
            ccz_emitted: 0,
            max_qubits_assert: Self::configured_max_qubits_assert(),
            max_ops_assert: Self::configured_max_ops_assert(),
            max_qubit_peak_assert: Self::configured_max_qubit_peak_assert(),
            last_touched_op: Vec::new(),
            last_op_on_qubit: Vec::new(),
            elide_deltas: std::collections::HashMap::new(),
            last_alloc_op_idx: 0,
            pending_frees: std::rc::Rc::new(std::cell::RefCell::new(Vec::new())),
            next_qubit: 0,
            next_bit: 0,
            free_qubits: std::collections::BTreeSet::new(),
            free_bits: Vec::new(),
            section_peak: std::collections::HashMap::new(),
            peak_qubits: 0,
            peak_bits: 0,
            peak_at_op: 0,
            peak_section: String::new(),
            peak_live_qubits: Vec::new(),
            peak_live_tags: Vec::new(),
            live_series: Vec::new(),
            sim: Some(vec![0u64; 4096]),
            sim_bits: Some(vec![0u64; 1024]),
            sim_condition_stack: Vec::new(),
            sim_phase_errors: 0,
            sim_phase_errors_by_tag: std::collections::HashMap::new(),
            current_section: String::from("(unset)"),
            sim_phase: 0,
            sim_hmr_counter: 0,
            sim_hmr_seed: hmr_seed,
            phase: PhaseTracker::new(),
            section_marks: Vec::new(),
            dump_hooks: std::collections::HashMap::new(),
            initial_sim_state: None,
            truncation_snapshot: None,
            live_checkpoints: Vec::new(),
            checkpoint_interval: 10_000_000,
            initial_hmr_seed: hmr_seed,
            regions: Vec::new(),
            executed_toffoli_shots: 0,
            executed_toffoli_by_section: std::collections::HashMap::new(),
            executed_ops_by_section: std::collections::HashMap::new(),
            scope_stack: Vec::new(),
            scope_frames_log: Vec::new(),
            op_scope: Vec::new(),
            qubit_tags: Vec::new(),
            tag_intern: std::collections::HashMap::new(),
            qubit_alloc_log: Vec::new(),
            bit_alloc_log: Vec::new(),
            ghost_event_log: Vec::new(),
            qubit_tag_counter: 0,
            deferred_contracts: std::rc::Rc::new(std::cell::RefCell::new(Vec::new())),
            capture_next_id: 0,
            call_site_intern: std::collections::HashMap::new(),
            pending_ghosts: Vec::new(),
            next_ghost_id: 0,
        }
    }

    /// Return (allocating once per distinct call site) the
    /// `scope_frames_log` `push_seq` for a `#[track_caller]`-captured
    /// gate-emission site `(file, line)`, nested under whatever manual
    /// `enter_scope!` frame is currently on the stack (or root). Reused
    /// across every op emitted from the same call site under the same
    /// manual scope, so this adds `O(#call_sites)` frames total, never
    /// O(#ops). Called from `push_gate_op` on the hot path; the common
    /// case is a single `HashMap` hit + return.
    fn intern_call_site(&mut self, file: &'static str, line: u32) -> u32 {
        let parent_seq = self.current_scope_seq();
        let parent = if parent_seq == u32::MAX {
            None
        } else {
            Some(parent_seq)
        };
        let key = (file, line, parent_seq);
        if let Some(&seq) = self.call_site_intern.get(&key) {
            return seq;
        }
        let push_seq = self.scope_frames_log.len() as u32;
        self.scope_frames_log.push(ScopeFrameLog {
            // Short, stable leaf name. The file:line is the load-bearing
            // part for `src`; the name just labels it as auto-captured.
            name: "gate@".to_string(),
            file,
            line,
            parent,
            start_op: self.ops.len(),
            end_op: None,
        });
        self.call_site_intern.insert(key, push_seq);
        push_seq
    }

    /// Push a scope frame with source location. Prefer the `enter_scope!`
    /// macro below, which captures file!() / line!() automatically.
    pub fn enter_scope_at(&mut self, name: &str, file: &'static str, line: u32) -> u32 {
        let push_seq = self.scope_frames_log.len() as u32;
        let parent = self.scope_stack.last().map(|f| f.push_seq);
        self.scope_frames_log.push(ScopeFrameLog {
            name: name.to_string(),
            file,
            line,
            parent,
            start_op: self.ops.len(),
            end_op: None,
        });
        self.scope_stack.push(ScopeFrame {
            name: name.to_string(),
            file,
            line,
            push_seq,
        });
        push_seq
    }

    /// Pop the innermost scope. `push_seq` must match the frame pushed.
    /// Mismatch panics — scope stack corruption is a bug, not a state
    /// we recover from.
    pub fn exit_scope(&mut self, push_seq: u32) {
        let top = self.scope_stack.pop().expect("exit_scope with empty stack");
        assert_eq!(
            top.push_seq, push_seq,
            "exit_scope push_seq mismatch (out-of-order)"
        );
        self.scope_frames_log[push_seq as usize].end_op = Some(self.ops.len());
    }

    /// Current innermost scope frame `push_seq`, or `u32::MAX` if no frame.
    fn current_scope_seq(&self) -> u32 {
        self.scope_stack.last().map_or(u32::MAX, |f| f.push_seq)
    }

    /// Render the scope chain (innermost first) for a given op index.
    #[must_use]
    pub fn op_source_trace(&self, op_idx: usize) -> Vec<(String, &'static str, u32)> {
        if op_idx >= self.op_scope.len() {
            return Vec::new();
        }
        let mut out = Vec::new();
        let mut cur = self.op_scope[op_idx];
        while cur != u32::MAX {
            let f = &self.scope_frames_log[cur as usize];
            out.push((f.name.clone(), f.file, f.line));
            cur = f.parent.unwrap_or(u32::MAX);
        }
        out
    }

    /// Tag a register for liveness analysis. Returns an index the
    /// caller passes to `untag_region` when the register is about
    /// to be freed. Cheap (one Vec push) — safe to sprinkle liberally.
    pub fn tag_region(&mut self, name: &str, qubits: &[QReg]) -> usize {
        let idx = self.regions.len();
        self.regions.push(TaggedRegion {
            name: name.to_string(),
            qubits: qubits.iter().map(|q| q.id).collect(),
            start_op: self.ops.len(),
            end_op: None,
        });
        idx
    }

    pub fn untag_region(&mut self, idx: usize) {
        self.regions[idx].end_op = Some(self.ops.len());
    }

    /// Lazily snapshot the sim state at the first gate emission.
    /// Captures whatever the caller loaded via `sim_load_reg_bytes` /
    /// `sim_load_bits_bytes`, so the debugger can replay from op 0
    /// with correct inputs.
    fn maybe_snapshot_initial(&mut self) {
        if self.initial_sim_state.is_some() {
            return;
        }
        if let (Some(q), Some(b)) = (&self.sim, &self.sim_bits) {
            self.initial_sim_state = Some((q.clone(), b.clone()));
        }
    }

    /// Total op count including any ops truncated in streaming mode.
    /// Use this instead of `self.ops.len()` when reporting op counts
    /// from long-running probes.
    #[must_use]
    pub fn total_ops(&self) -> u64 {
        self.ops_truncated + self.ops.len() as u64
    }

    #[must_use]
    pub fn elide_deltas_len(&self) -> usize {
        self.elide_deltas.len()
    }
    #[must_use]
    pub fn deferred_contracts_len(&self) -> usize {
        self.deferred_contracts.borrow().len()
    }

    /// Register a debugger dump hook. Invoke from debugger REPL with
    /// `dump <name> [args]`. Hook receives `&Debugger` (for qubit reads
    /// at current cursor + `section_marks` + scope chain) plus arg
    /// strings, returns the text to print.
    pub fn register_debug_dump<F>(&mut self, name: &str, hook: F)
    where
        F: Fn(&crate::point_add::trailmix_port::tracker::debugger::Debugger, &[&str]) -> String + Send + Sync + 'static,
    {
        self.dump_hooks.insert(name.to_string(), Box::new(hook));
    }

    #[must_use]
    pub fn debug_stats(&self) -> CircuitDebugStats {
        CircuitDebugStats {
            section_marks: self.section_marks.len(),
            section_mark_bytes: self.section_marks.iter().map(|(_, s)| s.len()).sum(),
            executed_ops_sections: self.executed_ops_by_section.len(),
            executed_tof_sections: self.executed_toffoli_by_section.len(),
            executed_section_key_bytes: self
                .executed_ops_by_section
                .keys()
                .map(std::string::String::len)
                .sum::<usize>()
                + self
                    .executed_toffoli_by_section
                    .keys()
                    .map(std::string::String::len)
                    .sum::<usize>(),
            regions: self.regions.len(),
            scope_frames_log: self.scope_frames_log.len(),
            op_scope: self.op_scope.len(),
            qubit_tags_len: self.qubit_tags.len(),
            qubit_tags_live: self.qubit_tags.iter().filter(|t| t.is_some()).count(),
            qubit_alloc_log: self.qubit_alloc_log.len(),
            qubit_alloc_tag_bytes: self.qubit_alloc_log.iter().map(|e| e.tag.len()).sum(),
        }
    }

    #[must_use]
    pub fn current_rss_kib() -> Option<u64> {
        let status = std::fs::read_to_string("/proc/self/status").ok()?;
        for line in status.lines() {
            if let Some(rest) = line.strip_prefix("VmRSS:") {
                return rest.split_whitespace().next().and_then(|s| s.parse().ok());
            }
        }
        None
    }

    fn truncate_debug_metadata(&mut self, boundary_op: usize) {
        self.section_marks.clear();
        self.section_marks
            .push((boundary_op, self.current_section.clone()));

        let live = self.next_qubit - self.free_qubits.len() as u32;
        self.live_series.clear();
        self.live_series.push((boundary_op, live));

        // Prune qubit_alloc_log entries that fall before the truncation
        // boundary AND are not the most-recent label for a currently-
        // live qubit. Persistent registers (e.g. pack0.cap_a, allocated
        // at op 0 and never freed) MUST keep their label so the debugger
        // can resolve `tag` / `p qreg` after a truncation at op 100M.
        // Without this carry-over the debugger renders "q770" instead of
        // "v4_inv.pack0.cap_a[0]", silently breaking post-truncation
        // inspection.
        let live: std::collections::HashSet<u32> = {
            let free: std::collections::HashSet<u32> = self.free_qubits.iter().copied().collect();
            (0..self.next_qubit).filter(|q| !free.contains(q)).collect()
        };
        let mut latest_pre: std::collections::HashMap<u32, usize> =
            std::collections::HashMap::new();
        for (i, ev) in self.qubit_alloc_log.iter().enumerate() {
            if ev.op_idx < boundary_op && live.contains(&ev.qubit) {
                latest_pre
                    .entry(ev.qubit)
                    .and_modify(|j| {
                        if self.qubit_alloc_log[*j].op_idx < ev.op_idx {
                            *j = i;
                        }
                    })
                    .or_insert(i);
            }
        }
        let keep_indices: std::collections::HashSet<usize> = latest_pre.values().copied().collect();
        let mut idx: usize = 0;
        self.qubit_alloc_log.retain(|ev| {
            let i = idx;
            idx += 1;
            ev.op_idx >= boundary_op || keep_indices.contains(&i)
        });

        // Parallel prune for bit_alloc_log: keep the most-recent label
        // for each currently-live cbit, drop earlier-than-boundary
        // entries otherwise. Without this, named cbits allocated early
        // would render as `b{N}` after a streaming truncation.
        let live_bits: std::collections::HashSet<u32> = {
            let free: std::collections::HashSet<u32> = self.free_bits.iter().copied().collect();
            (0..self.next_bit).filter(|b| !free.contains(b)).collect()
        };
        let mut latest_pre_bit: std::collections::HashMap<u32, usize> =
            std::collections::HashMap::new();
        for (i, ev) in self.bit_alloc_log.iter().enumerate() {
            if ev.op_idx < boundary_op && live_bits.contains(&ev.bit) {
                latest_pre_bit
                    .entry(ev.bit)
                    .and_modify(|j| {
                        if self.bit_alloc_log[*j].op_idx < ev.op_idx {
                            *j = i;
                        }
                    })
                    .or_insert(i);
            }
        }
        let keep_bit_indices: std::collections::HashSet<usize> =
            latest_pre_bit.values().copied().collect();
        let mut idx_b: usize = 0;
        self.bit_alloc_log.retain(|ev| {
            let i = idx_b;
            idx_b += 1;
            ev.op_idx >= boundary_op || keep_bit_indices.contains(&i)
        });
    }

    /// Push a gate Op and snapshot initial sim state if this is the
    /// first gate. All gate-emitting methods MUST call this (rather
    /// than self.ops.push directly) so the debugger's initial state
    /// captures pre-gate sim exactly. Register/append ops do not
    /// count as gates (they don't touch sim state).
    ///
    /// Returns `true` if the redundant-op auto-eliminator detected
    /// that this op cancels with a same-shape prior op (no operand
    /// changed since) and undid the prior emission. The caller MUST
    /// still apply the gate's sim/phase effect — for involutary gates
    /// the second application reverses the first, leaving sim and the
    /// phase tracker in their pre-prior state. The op itself is NOT
    /// pushed onto self.ops.
    ///
    /// `#[track_caller]` so `Location::caller()` resolves to the FIRST
    /// non-`#[track_caller]` frame above us — i.e. the primitive code
    /// that emitted the gate (`cursor_pack_v4.rs`, adders, …). The public
    /// gate methods (`x`/`cx`/`ccx`/…) and their `*_internal` helpers
    /// are also `#[track_caller]`, so the propagated location skips the
    /// Circuit wrapper layers and points at the real call site. This is
    /// what populates `scope_frames_log` in `--release` builds with no
    /// manual `enter_scope!` tags (the historical reason `src <op>` in
    /// the debugger always reported "no scope frame").
    #[track_caller]
    fn push_gate_op(&mut self, op: Op) -> bool {
        let caller = std::panic::Location::caller();
        let call_site_seq = self.intern_call_site(caller.file(), caller.line());
        // Drain any QReg-dropped qubits BEFORE we update current_op_idx
        // for the new gate. The drain calls free_qubit at the
        // current_op_idx of the LAST emitted gate, so a QReg dropped
        // right after its last touch gets gap=0.
        self.drain_pending_frees();
        self.maybe_snapshot_initial();
        // Live debugger checkpoint: at this point the PREVIOUS gate's
        // sim/phase apply has completed, so self.sim_* fields hold the
        // exact state AFTER `n_emitted` ops (= what cp.op_idx tracks
        // in the debugger). Snapshot every `checkpoint_interval` ops
        // so the debugger can attach in O(1) instead of replaying the
        // whole op stream.
        let n_emitted = self.ops_truncated as usize + self.ops.len();
        if let Some(sim) = self.sim.as_ref() {
            if self.checkpoint_interval > 0
                && n_emitted > 0
                && n_emitted.is_multiple_of(self.checkpoint_interval)
                && self
                    .live_checkpoints
                    .last()
                    .is_none_or(|c| (c.ops_truncated as usize) < n_emitted)
            {
                self.live_checkpoints.push(TruncationSnapshot {
                    ops_truncated: n_emitted as u64,
                    qubits: sim.clone(),
                    bits: self.sim_bits.clone().unwrap_or_default(),
                    phase: self.sim_phase,
                    cond_stack: self.sim_condition_stack.clone(),
                    hmr_counter: self.sim_hmr_counter,
                    r_on_nonzero_events: u64::from(self.sim_phase_errors),
                });
            }
        }
        self.phase.current_op_idx = self.ops.len() + self.ops_truncated as usize;
        let cur_op = self.phase.current_op_idx as u64;

        if let Some(prior_idx) = self.find_redundant_prior(&op) {
            self.elide_prior(prior_idx);
            return true;
        }

        // Ops that change the global condition mask or global phase
        // invalidate all pending elide candidates: a future "same"
        // gate on the same operands would actually fire under a
        // different cond and so wouldn't cancel its prior emission.
        if matches!(op, Op::PushCondition(_) | Op::PopCondition | Op::Neg) {
            self.elide_deltas.clear();
        }

        // Capture undo info before overwriting last_touched/last_op for
        // each operand. Only useful for involutary gates that could
        // become future elide candidates; for non-involutary ops we
        // skip the HashMap insertion entirely.
        let is_involutary = Self::normalized_for_redundancy(op).is_some();

        let mut delta_slots: [Option<(u32, Option<u64>, Option<Op>)>; 3] = [None, None, None];
        for (slot_idx, q) in op.touched_qubits().enumerate() {
            let qi = q as usize;
            assert!(
                q < self.next_qubit,
                "push_gate_op: op references qubit id q{} but next_qubit \
                 is only {}, peak_qubits={}, current_section='{}', \
                 op_idx={}, op={:?}",
                q,
                self.next_qubit,
                self.peak_qubits,
                self.current_section,
                cur_op,
                op,
            );
            if qi >= self.last_touched_op.len() {
                self.last_touched_op.resize(qi + 1, None);
            }
            if qi >= self.last_op_on_qubit.len() {
                self.last_op_on_qubit.resize(qi + 1, None);
            }
            let prev_lt = self.last_touched_op[qi];
            let prev_lop = self.last_op_on_qubit[qi];
            if is_involutary {
                delta_slots[slot_idx] = Some((q, prev_lt, prev_lop));
            }
            // Pruning: any op the OLD last_touched_op pointed at is no
            // longer a valid elide candidate (one of its operands just
            // moved on). Drop its delta to keep the HashMap bounded.
            if let Some(old_idx) = prev_lt {
                self.elide_deltas.remove(&old_idx);
            }
            self.last_touched_op[qi] = Some(cur_op);
            self.last_op_on_qubit[qi] = Some(op);
        }
        if is_involutary {
            self.elide_deltas.insert(
                cur_op,
                ElideDelta {
                    operand_prev: delta_slots,
                    section: self.current_section.clone(),
                },
            );
        }
        // Point this op at its auto-captured call-site leaf frame. That
        // frame's `parent` is the manual `enter_scope!` chain (if any),
        // so the debugger's `src` walk yields  call-site → manual scopes
        // → root, instead of the old behavior where this was the bare
        // manual `current_scope_seq()` (= u32::MAX whenever no caller had
        // sprinkled `enter_scope!`, which was always, in real circuits).
        self.op_scope.push(call_site_seq);
        *self
            .executed_ops_by_section
            .entry(self.current_section.clone())
            .or_insert(0) += 1;
        match op {
            Op::Ccx(_, _, _) => self.ccx_emitted += 1,
            Op::Ccz(_, _, _) => self.ccz_emitted += 1,
            _ => {}
        }

        // Streaming mode is unconditional: the circuit never keeps an
        // unbounded op buffer live. Each clear adds to ops_truncated
        // so total_ops() still reports the true count. The debugger
        // still works on the retained tail — we capture the sim state
        // AT truncation into `truncation_snapshot` so it can serve as
        // the debugger's starting state for replaying the tail.
        //
        // The truncate check runs BEFORE the push so that a push at
        // exactly `ops_cap` doesn't trigger Vec's doubling growth
        // (which would request a 2× transient allocation that
        // overflows the process memory limit at 100M-op caps).
        let cap = self.ops_cap;
        if self.ops.len() >= cap {
            let new_truncated = self.ops_truncated + self.ops.len() as u64;
            if let Some(sim) = &self.sim {
                self.truncation_snapshot = Some(TruncationSnapshot {
                    ops_truncated: new_truncated,
                    qubits: sim.clone(),
                    bits: self.sim_bits.clone().unwrap_or_default(),
                    phase: self.sim_phase,
                    cond_stack: self.sim_condition_stack.clone(),
                    hmr_counter: self.sim_hmr_counter,
                    r_on_nonzero_events: u64::from(self.sim_phase_errors),
                });
            }
            self.ops_truncated = new_truncated;
            // Drop live checkpoints that fall in or before the truncated
            // range. Their target ops are no longer in `self.ops`, so the
            // debugger can't replay forward from them. The
            // `truncation_snapshot` captured just above is the only
            // valid starting point at the new ops_start_idx.
            self.live_checkpoints
                .retain(|c| c.ops_truncated > new_truncated);
            self.ops.clear();
            self.op_scope.clear();
            // Drop all elide deltas: ops in the truncated range can no
            // longer be undone (they're not in self.ops anymore), and
            // their op_idx values are <= new_truncated.
            self.elide_deltas.clear();
            self.truncate_debug_metadata(new_truncated as usize);
        }
        self.ops.push(Some(op));

        // Memory watchdog: every ~1024 ops, read VmRSS and panic if over
        // threshold so we get a main-thread backtrace with section context
        // BEFORE systemd SIGKILLs us at the 8G unit cap (which produces no
        // useful diagnostic). Default threshold is 6 GiB; override via
        // CIRC_MEM_WATCHDOG_MB=<n>.
        let n = self.ops.len() + self.ops_truncated as usize;
        self.maybe_assert_max_ops(n as u64);
        if n.trailing_zeros() >= 10 {
            let limit_mb: u64 = std::env::var("CIRC_MEM_WATCHDOG_MB")
                .ok()
                .and_then(|s| s.parse().ok())
                .unwrap_or(6000);
            if let Ok(status) = std::fs::read_to_string("/proc/self/status") {
                for line in status.lines() {
                    if let Some(rest) = line.strip_prefix("VmRSS:") {
                        let kb: u64 = rest
                            .split_whitespace()
                            .next()
                            .and_then(|s| s.parse().ok())
                            .unwrap_or(0);
                        let mb = kb / 1024;
                        assert!(
                            mb <= limit_mb,
                            "CIRC_MEM_WATCHDOG: RSS {}MB > {}MB at op {}, section '{}', peak {}q",
                            mb,
                            limit_mb,
                            n,
                            self.current_section,
                            self.peak_qubits
                        );
                    }
                }
            }
        }
        false
    }

    /// Look up the previous Op that touched qubit `q`, if any.
    fn last_op_for(&self, q: u32) -> Option<Op> {
        self.last_op_on_qubit.get(q as usize).copied().flatten()
    }

    /// Normalize an involutary gate's operand ordering so that gates
    /// representing the same physical operation compare equal even
    /// when emitted with different operand orderings (CCX is symmetric
    /// in its two controls; CZ and CCZ are fully symmetric in their
    /// operands). Returns `None` for non-involutary ops.
    fn normalized_for_redundancy(op: Op) -> Option<Op> {
        match op {
            Op::X(_) | Op::Z(_) => Some(op),
            Op::Cx(c, t) => Some(Op::Cx(c, t)),
            Op::Cz(a, b) => {
                let (a, b) = if a <= b { (a, b) } else { (b, a) };
                Some(Op::Cz(a, b))
            }
            Op::Ccx(c1, c2, t) => {
                let (c1, c2) = if c1 <= c2 { (c1, c2) } else { (c2, c1) };
                Some(Op::Ccx(c1, c2, t))
            }
            Op::Ccz(a, b, c) => {
                let mut v = [a, b, c];
                v.sort_unstable();
                Some(Op::Ccz(v[0], v[1], v[2]))
            }
            _ => None,
        }
    }

    /// Flag involutary-gate redundancy. For each involutary gate,
    /// every qubit it touches must show that the SAME prior gate
    /// touched it, with no other operand having been modified in
    /// between. If all operands match the same normalized prior op,
    /// the new gate is the inverse of an uncancelled prior gate and
    /// the pair cancels to identity — a bug.
    /// If `op` is involutary and every one of its operands shows the
    /// SAME prior gate as its `last_op_on_qubit`, return that prior
    /// gate's `op_idx` — the pair cancels and can be auto-elided.
    fn find_redundant_prior(&self, op: &Op) -> Option<u64> {
        let norm_cur = Self::normalized_for_redundancy(*op)?;
        let mut prior_idx: Option<u64> = None;
        for q in op.touched_qubits() {
            let prev = self.last_op_for(q)?;
            let norm_prev = Self::normalized_for_redundancy(prev)?;
            if norm_prev != norm_cur {
                return None;
            }
            let prev_idx = self
                .last_touched_op
                .get(q as usize)
                .and_then(|x| x.as_ref())
                .copied()?;
            match prior_idx {
                None => prior_idx = Some(prev_idx),
                Some(p) if p != prev_idx => return None,
                _ => {}
            }
        }
        let p_idx = prior_idx?;
        // Need delta to undo the prior op. Delta is pruned aggressively
        // when ops overshadow each other; if it's gone, we can't safely
        // elide (the pre-prior values for restoring last_touched_op
        // are unknown). Skip the elide in that case.
        if !self.elide_deltas.contains_key(&p_idx) {
            return None;
        }
        // The prior op must still be in the retained op buffer
        // (post-truncation slot). Ops that have been streamed away are
        // not elidable — their counters have already escaped into
        // ops_truncated.
        if p_idx < self.ops_truncated {
            return None;
        }
        Some(p_idx)
    }

    /// Undo the involutary op at `prior_idx` so the about-to-arrive
    /// new op (which the caller will sim/phase-apply normally) ends up
    /// being the second half of an X-X / CCX-CCX cancellation pair.
    /// Sim/phase reversal is handled by the caller's normal apply
    /// (involutary gates: applying twice is identity).
    fn elide_prior(&mut self, prior_idx: u64) {
        let pos = (prior_idx - self.ops_truncated) as usize;
        let prior_op = self
            .ops
            .get(pos)
            .and_then(|s| *s)
            .expect("elide_prior: prior op missing from retained buffer");
        let delta = self
            .elide_deltas
            .remove(&prior_idx)
            .expect("elide_prior: delta missing for retained op");
        // Mark the prior op slot as a no-op (None). kmx emit and
        // sim/phase replay both skip None. The slot stays so op_idx
        // values used elsewhere (last_touched_op for other qubits,
        // op_scope, debugger references) remain stable.
        self.ops[pos] = None;
        // Restore each operand's last_touched_op and last_op_on_qubit
        // to its pre-prior value. Other qubits (not operands of the
        // elided pair) are unaffected.
        for slot in &delta.operand_prev {
            if let Some((q, prev_lt, prev_lop)) = *slot {
                let qi = q as usize;
                if qi < self.last_touched_op.len() {
                    self.last_touched_op[qi] = prev_lt;
                }
                if qi < self.last_op_on_qubit.len() {
                    self.last_op_on_qubit[qi] = prev_lop;
                }
                // If the pre-prior op idx points at an op whose delta
                // was pruned during emission of `prior_op`, that op is
                // re-exposed as the "candidate" for further elide
                // attempts but lacks a delta — so future elides
                // against it are silently skipped (no correctness
                // issue, just a missed optimization).
            }
        }
        // Counter undo for the prior emission.
        match prior_op {
            Op::Ccx(_, _, _) => self.ccx_emitted = self.ccx_emitted.saturating_sub(1),
            Op::Ccz(_, _, _) => self.ccz_emitted = self.ccz_emitted.saturating_sub(1),
            _ => {}
        }
        if let Some(c) = self.executed_ops_by_section.get_mut(&delta.section) {
            *c = c.saturating_sub(1);
        }
    }

    fn sim_ensure(&mut self, id: u32) {
        if let Some(ref mut v) = self.sim {
            let id = id as usize;
            if id >= v.len() {
                v.resize(id + 256, 0);
            }
        }
    }

    /// Read a qubit's 64-bit bitmask (bit i = shot i's value). 0 if
    /// sim is disabled.
    pub(crate) fn sim_get_mask(&self, id: u32) -> u64 {
        self.sim
            .as_ref()
            .map_or(0, |v| v.get(id as usize).copied().unwrap_or(0))
    }

    /// Snapshot the live simulator state into a `DestroyedSimState`
    /// for post-run inspection. `outputs` lists the `QRegs` the caller
    /// is keeping past circuit termination — each is detached to a
    /// plain `Vec<Qubit>` and returned alongside the snapshot. Any
    /// `QRegs` NOT in `outputs` must already be dropped; otherwise the
    /// call panics (qubit leak detector). For circuits that allocate
    /// only via `alloc_qubits` (raw `Vec<Qubit>`, no `QReg` ownership),
    /// pass `vec![]`.
    ///
    /// Returns `(snapshot, detached_outputs)` where `detached_outputs[i]`
    /// is the qubit list of `outputs[i]`.
    pub fn destroy_sim(&mut self, outputs: Vec<QReg>) -> (DestroyedSimState, Vec<QReg>) {
        // Phase obligations must be accounted for, exactly like QRegs: a circuit
        // torn down with destroy_sim must have NO pending ghosts (resolve_ghost
        // them, or hand them over via destroy_sim_ghosts). This is the empty-
        // ghosts case of destroy_sim_ghosts.
        assert!(
            self.pending_ghosts.is_empty(),
            "destroy_sim: {} ghost(s) still pending — resolve_ghost them, or use \
             destroy_sim_ghosts(outputs, ghosts) to hand over every ghost",
            self.pending_ghosts.len()
        );
        // Drain anything queued by QRegs already dropped before this
        // call. Each free goes through the strict-dealloc check; a
        // gap or sim-nonzero violation surfaces here as a panic.
        self.drain_pending_frees();
        // After drain, strong_count(&pending_frees) = 1 (Circuit) +
        // (one per remaining live QReg, since QReg is linear — no
        // Clone, no Copy). The caller must list every remaining QReg
        // as either an output (here) or have already dropped it.
        let strong = std::rc::Rc::strong_count(&self.pending_frees);
        let live = strong - 1;
        assert!(
            live == outputs.len(),
            "destroy_sim: {} QReg(s) live but {} listed as outputs. \
             Every remaining QReg must be passed as an output OR \
             dropped before destroy_sim.",
            live,
            outputs.len(),
        );
        // Mark each output detached: still a QReg, but with no auto-free
        // (the Circuit has been consumed, so there's no pending queue
        // for Drop to push to). This is the ONLY remaining call site
        // that flips `detached` — every other QReg must be explicitly
        // freed via `circ.zero_and_free`.
        let detached: Vec<QReg> = outputs
            .into_iter()
            .map(|mut q| {
                q.detached = true;
                q
            })
            .collect();
        let snapshot = DestroyedSimState {
            qubits: self.sim.clone().unwrap_or_default(),
            bits: self.sim_bits.clone().unwrap_or_default(),
            phase: self.sim_phase,
            phase_errors: self.sim_phase_errors,
            phase_errors_by_tag: self.sim_phase_errors_by_tag.clone(),
        };
        (snapshot, detached)
    }

    fn sim_bit_ensure(&mut self, id: u32) {
        if let Some(ref mut v) = self.sim_bits {
            let id = id as usize;
            if id >= v.len() {
                v.resize(id + 256, 0);
            }
        }
    }

    /// Read a classical bit's 64-bit bitmask.
    fn sim_bit_get_mask(&self, id: u32) -> u64 {
        self.sim_bits
            .as_ref()
            .map_or(0, |v| v.get(id as usize).copied().unwrap_or(0))
    }

    /// Read a bit's value in shot 0 (for legacy callers).
    /// Set a classical bit's value in shot 0 (broadcast to all shots
    /// if you want identical inputs across shots — use `sim_bit_set_mask`
    /// for per-shot control). Mostly used for unit tests.
    pub fn sim_bit_set(&mut self, id: u32, val: u8) {
        if let Some(ref mut v) = self.sim_bits {
            let id = id as usize;
            if id >= v.len() {
                v.resize(id + 256, 0);
            }
            // Broadcast: set all 64 shots to the same value.
            v[id] = if val & 1 == 1 { u64::MAX } else { 0 };
        }
    }

    /// Current condition mask (AND of all stacked bits). `u64::MAX`
    /// when stack is empty.
    fn sim_condition_mask(&self) -> u64 {
        let mut mask = u64::MAX;
        for &b in &self.sim_condition_stack {
            mask &= self.sim_bit_get_mask(b);
        }
        mask
    }

    /// Allocate a qubit with a descriptive name. The name is scoped
    /// under the current section (or `enter_scope!` frame). Qubit-id
    /// recycling from the free pool is fine — the NEW name takes
    /// over; the debugger's `qubit_name_at` walks the alloc log to
    /// show the right tag for a given `op_idx`.
    fn alloc_qubit(&mut self, name: &str) -> Qubit {
        // Drain QReg-dropped queue first: a QReg dropped between the
        // last gate emission and this alloc must be freed BEFORE we
        // advance last_alloc_op_idx, otherwise the deferred free's
        // strict-dealloc check will panic against last_alloc set by
        // *this* alloc rather than the prior one.
        self.drain_pending_frees();
        let q = self.free_qubits.pop_first().unwrap_or_else(|| {
            let q = self.next_qubit;
            self.next_qubit += 1;
            // Sanity bound: the circuit allocator should reuse freed slots
            // via the free_qubits pool, so next_qubit corresponds roughly
            // to peak_live_qubits over the run. Anything past 10M means a
            // QReg drop path is broken and IDs are not being reclaimed.
            // Override the cap via CIRC_NEXT_QUBIT_CAP=<n> if a legitimate
            // workload genuinely needs more.
            const DEFAULT_NEXT_QUBIT_CAP: u32 = 10_000_000;
            let cap: u32 = std::env::var("CIRC_NEXT_QUBIT_CAP")
                .ok()
                .and_then(|s| s.parse().ok())
                .unwrap_or(DEFAULT_NEXT_QUBIT_CAP);
            assert!(
                self.next_qubit <= cap,
                "alloc_qubit: next_qubit reached {} (cap {}). \
                 free_qubits pool size = {}, peak_qubits = {}, \
                 section = '{}'. The QReg drop->free reclaim path is \
                 broken: IDs are accumulating instead of reusing freed \
                 slots. Override via CIRC_NEXT_QUBIT_CAP=<n> if this is \
                 intentional.",
                self.next_qubit,
                cap,
                self.free_qubits.len(),
                self.peak_qubits,
                self.current_section,
            );
            q
        });
        let live = self.next_qubit - self.free_qubits.len() as u32;
        let new_peak = live > self.peak_qubits;
        {
            let cs = self.current_section.clone();
            let sp = self.section_peak.entry(cs).or_insert(0);
            if live > *sp {
                *sp = live;
            }
        }
        self.live_series
            .push((self.ops.len() + self.ops_truncated as usize, live));
        if let Some(ref mut v) = self.sim {
            let qi = q as usize;
            if qi >= v.len() {
                v.resize(qi + 256, 0);
            }
            v[qi] = 0;
        }
        let tag = self.scoped_tag(name);
        self.record_qubit_alloc(q, tag);
        // Reset last_touched_op: alloc by itself is NOT a touch.
        // free_qubit will panic on a qubit that's never been touched
        // by any gate (untouched alloc → wasted ancilla budget).
        let qi = q as usize;
        if qi >= self.last_touched_op.len() {
            self.last_touched_op.resize(qi + 1, None);
        } else {
            self.last_touched_op[qi] = None;
        }
        // Reset last_op_on_qubit so a fresh allocation can't trigger a
        // false-positive redundant-op flag against a stale Op left
        // behind by a previous owner of this slot.
        if qi >= self.last_op_on_qubit.len() {
            self.last_op_on_qubit.resize(qi + 1, None);
        } else {
            self.last_op_on_qubit[qi] = None;
        }
        // Record the alloc edge — free_qubit's wasteful-retention
        // check compares each freed qubit's last_touched_op against
        // this. If a qubit's last gate-touch precedes the most
        // recent alloc, it could (and should) have been freed BEFORE
        // that alloc, letting the allocator reuse its slot.
        self.last_alloc_op_idx = self.phase.current_op_idx as u64;
        if new_peak {
            self.peak_qubits = live;
            self.peak_at_op = self.ops.len() + self.ops_truncated as usize;
            self.peak_section = self.current_section.clone();
            self.snapshot_peak_live();
        }
        self.phase.on_alloc_qubit(q);
        self.maybe_assert_max_qubits(live);
        self.maybe_assert_max_qubit_peak(live);
        Qubit(q)
    }

    /// Prefix a local name with the innermost scope/section so
    /// debugger output unambiguously identifies the site.
    fn scoped_tag(&mut self, name: &str) -> String {
        let scope = if let Some(f) = self.scope_stack.last() {
            f.name.as_str()
        } else {
            self.current_section.as_str()
        };
        if scope.is_empty() {
            name.to_string()
        } else {
            format!("{scope}/{name}")
        }
    }

    /// Intern a tag string into a shared `Arc<str>`. Repeated allocations
    /// of identically-named scratch qubits share heap storage.
    fn intern_tag(&mut self, tag: String) -> Arc<str> {
        if let Some(existing) = self.tag_intern.get(&tag) {
            return existing.clone();
        }
        let arc: Arc<str> = Arc::from(tag.as_str());
        self.tag_intern.insert(tag, arc.clone());
        arc
    }

    fn record_qubit_alloc(&mut self, q: u32, tag: String) {
        let qi = q as usize;
        if qi >= self.qubit_tags.len() {
            self.qubit_tags.resize(qi + 1, None);
        }
        let interned = self.intern_tag(tag);
        self.qubit_tags[qi] = Some(interned.clone());
        self.qubit_alloc_log.push(QubitAllocEvent {
            op_idx: self.ops.len() + self.ops_truncated as usize,
            qubit: q,
            tag: interned,
        });
    }

    /// Record a classical-bit alloc event for the debugger's
    /// `bit_alloc_log`. Parallel to `record_qubit_alloc`. Tags are
    /// interned through the shared `tag_intern` map.
    fn record_bit_alloc(&mut self, b: u32, tag: String) {
        let interned = self.intern_tag(tag);
        self.bit_alloc_log.push(BitAllocEvent {
            op_idx: self.ops.len() + self.ops_truncated as usize,
            bit: b,
            tag: interned,
        });
    }

    /// Tag-resolution helper for cbits, mirroring `qubit_name_at_id`.
    /// Walks `bit_alloc_log` for the most recent alloc at or before
    /// `op_idx`. Falls back to `b{N}` if none exists.
    #[allow(dead_code)]
    #[must_use]
    pub fn cbit_name_at_id(&self, bit_id: u32, op_idx: usize) -> String {
        let mut best: Option<&str> = None;
        let mut best_op: usize = 0;
        for ev in &self.bit_alloc_log {
            if ev.bit == bit_id && ev.op_idx <= op_idx && (best.is_none() || ev.op_idx >= best_op) {
                best = Some(&*ev.tag);
                best_op = ev.op_idx;
            }
        }
        best.map_or_else(|| format!("b{bit_id}"), std::string::ToString::to_string)
    }

    /// Override the auto-generated tag for a single qubit with a
    /// semantic name. Records a new alloc-log entry so the debugger
    /// sees the name change at this op.
    #[allow(dead_code)]
    fn name_qubit(&mut self, q: Qubit, name: &str) {
        self.record_qubit_alloc(q.raw(), name.to_string());
    }

    /// Tag each qubit in `qs` as `name[0]`, `name[1]`, …
    #[allow(dead_code)]
    fn name_qubits(&mut self, qs: &[Qubit], name: &str) {
        for (i, &q) in qs.iter().enumerate() {
            self.record_qubit_alloc(q.raw(), format!("{name}[{i}]"));
        }
    }

    /// Return the most recent tag for qubit `q`, or a fallback.
    #[allow(dead_code)]
    fn qubit_name(&self, q: Qubit) -> String {
        self.qubit_tags
            .get(q.raw() as usize)
            .and_then(|t| t.as_ref())
            .map_or_else(|| format!("q{}", q.raw()), std::string::ToString::to_string)
    }

    /// Return the tag that was current for `q` at op `op_idx`.
    /// Walks the alloc log to find the most recent alloc event for
    /// this qubit at or before `op_idx`.
    #[allow(dead_code)]
    fn qubit_name_at(&self, q: Qubit, op_idx: usize) -> String {
        self.qubit_name_at_id(q.raw(), op_idx)
    }

    /// Public version of `qubit_name_at` that takes a raw qubit id
    /// (e.g. an entry of `peak_live_qubits`). Diagnostic only.
    #[must_use]
    pub fn qubit_name_at_id(&self, qubit_id: u32, op_idx: usize) -> String {
        let mut best: Option<&str> = None;
        let mut best_op: usize = 0;
        for ev in &self.qubit_alloc_log {
            if ev.qubit == qubit_id
                && ev.op_idx <= op_idx
                && (best.is_none() || ev.op_idx >= best_op)
            {
                best = Some(&*ev.tag);
                best_op = ev.op_idx;
            }
        }
        best.map_or_else(|| format!("q{qubit_id}"), std::string::ToString::to_string)
    }

    /// Allocate a qubit that has never been used before (skip
    /// the free pool). Useful when register qubits must not
    /// overlap with previously-freed internal state.
    #[allow(dead_code)]
    fn alloc_qubit_fresh(&mut self, name: &str) -> Qubit {
        // See alloc_qubit: drain queued QReg drops first so any
        // strict-dealloc panic surfaces against the prior alloc.
        self.drain_pending_frees();
        let q = self.next_qubit;
        self.next_qubit += 1;
        let live = self.next_qubit - self.free_qubits.len() as u32;
        let new_peak = live > self.peak_qubits;
        {
            let cs = self.current_section.clone();
            let sp = self.section_peak.entry(cs).or_insert(0);
            if live > *sp {
                *sp = live;
            }
        }
        self.live_series
            .push((self.ops.len() + self.ops_truncated as usize, live));
        let tag = self.scoped_tag(name);
        self.record_qubit_alloc(q, tag);
        let qi = q as usize;
        if qi >= self.last_touched_op.len() {
            self.last_touched_op.resize(qi + 1, None);
        } else {
            self.last_touched_op[qi] = None;
        }
        if qi >= self.last_op_on_qubit.len() {
            self.last_op_on_qubit.resize(qi + 1, None);
        } else {
            self.last_op_on_qubit[qi] = None;
        }
        self.last_alloc_op_idx = self.phase.current_op_idx as u64;
        if new_peak {
            self.peak_qubits = live;
            self.peak_at_op = self.ops.len() + self.ops_truncated as usize;
            self.peak_section = self.current_section.clone();
            self.snapshot_peak_live();
        }
        self.phase.on_alloc_qubit(q);
        self.maybe_assert_max_qubits(live);
        self.maybe_assert_max_qubit_peak(live);
        Qubit(q)
    }

    fn maybe_assert_max_ops(&mut self, total_ops: u64) {
        let Some(limit) = self.max_ops_assert else {
            return;
        };
        if total_ops <= limit {
            return;
        }

        let mut top_exact = self
            .executed_ops_by_section
            .iter()
            .map(|(section, ops)| (section.as_str(), *ops))
            .collect::<Vec<_>>();
        top_exact.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(b.0)));

        let mut prefix_ops = std::collections::BTreeMap::<String, u64>::new();
        for (section, ops) in &self.executed_ops_by_section {
            let head = section.split('/').next().unwrap_or(section).to_string();
            *prefix_ops.entry(head).or_default() += *ops;
        }
        let mut top_prefix = prefix_ops.into_iter().collect::<Vec<_>>();
        top_prefix.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

        let format_rows = |rows: &[(String, u64)]| -> String {
            rows.iter()
                .take(8)
                .fold(String::new(), |mut acc, (section, ops)| {
                    let _ = write!(acc, "\n  {ops:>12}  {section}");
                    acc
                })
        };
        let exact_rows = top_exact
            .iter()
            .take(8)
            .fold(String::new(), |mut acc, (section, ops)| {
                let _ = write!(acc, "\n  {ops:>12}  {section}");
                acc
            });
        let prefix_rows = format_rows(&top_prefix);

        let msg = format!(
            "MAX_OPS exceeded: total_ops={} limit={} section='{}' peak={} peak_section='{}'\
             \nTop section prefixes by emitted ops:{}\
             \nTop exact sections by emitted ops:{}",
            total_ops,
            limit,
            self.current_section,
            self.peak_qubits,
            self.peak_section,
            prefix_rows,
            exact_rows,
        );
        if std::env::var("DEBUG_ON_FAIL").is_ok() {
            eprintln!("{msg}");
            eprintln!("[debugger] DEBUG_ON_FAIL set — attaching at MAX_OPS panic");
            // Capture target BEFORE attach: attach moves self.ops into d,
            // so self.ops.len() afterwards is 0 and d.goto(0) would back
            // the cursor up to op #0. attach already sets cursor=ops.len()
            // = end-of-stream, which matches the panic point.
            let target = self.ops.len() + self.ops_truncated as usize;
            let mut d = crate::point_add::trailmix_port::debugger::Debugger::attach(self);
            d.goto(target);
            d.repl();
        }
        panic!("{}", msg);
    }

    /// Peak-qubit cap with a built-in default of 1175. See
    /// [`Circuit::max_qubit_peak_assert`] for semantics. Fires the
    /// instant `live` exceeds the cap, and attaches the debugger
    /// under `DEBUG_ON_FAIL` exactly like `maybe_assert_max_ops`.
    fn maybe_assert_max_qubit_peak(&mut self, live: u32) {
        let Some(limit) = self.max_qubit_peak_assert else {
            return;
        };
        if live <= limit {
            return;
        }

        let msg = format!(
            "MAX_QUBIT_PEAK exceeded: live={} limit={} section='{}' op={} peak={} peak_section='{}' \
             (override via CIRC_ASSERT_MAX_QUBIT_PEAK=<n>, set to 0 to disable; default {}). \
             NOTE: this aborts at the FIRST allocation that crosses the limit, so `live` is a \
             LOWER BOUND on the true peak (the circuit would go higher if allowed to continue) -- \
             do NOT read `live - limit` as 'how many qubits over'; profile the live-tag histogram \
             below to see what is alive here and reduce it.",
            live,
            limit,
            self.current_section,
            self.ops.len() + self.ops_truncated as usize,
            self.peak_qubits,
            self.peak_section,
            Self::DEFAULT_MAX_QUBIT_PEAK,
        );
        let mut tag_classes = std::collections::BTreeMap::<String, usize>::new();
        for tag in &self.peak_live_tags {
            let class = Self::normalized_tag_class(tag);
            *tag_classes.entry(class).or_default() += 1;
        }
        let mut tag_classes = tag_classes.into_iter().collect::<Vec<_>>();
        tag_classes.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
        if std::env::var("DEBUG_ON_FAIL").is_ok() {
            eprintln!("{msg}");
            eprintln!(
                "[debugger] peak_live_tag_class_histogram ({} classes):",
                tag_classes.len()
            );
            for (class, count) in &tag_classes {
                eprintln!("  {count:>5} {class}");
            }
            eprintln!("[debugger] DEBUG_ON_FAIL set — attaching at MAX_QUBIT_PEAK panic");
            // See MAX_OPS site: capture target BEFORE attach (attach moves ops).
            let target = self.ops.len() + self.ops_truncated as usize;
            let mut d = crate::point_add::trailmix_port::debugger::Debugger::attach(self);
            d.goto(target);
            d.repl();
        } else {
            eprintln!("{msg}");
            eprintln!(
                "[debugger] peak_live_tag_class_histogram ({} classes):",
                tag_classes.len()
            );
            for (class, count) in &tag_classes {
                eprintln!("  {count:>5} {class}");
            }
        }
        panic!("{msg}");
    }

    fn maybe_assert_max_qubits(&mut self, live: u32) {
        let Some(limit) = self.max_qubits_assert else {
            return;
        };
        if live <= limit {
            return;
        }

        let msg = format!(
            "MAX_QUBITS exceeded: live={} limit={} section='{}' op={} peak={} peak_section='{}'",
            live,
            limit,
            self.current_section,
            self.ops.len() + self.ops_truncated as usize,
            self.peak_qubits,
            self.peak_section,
        );
        let mut tag_classes = std::collections::BTreeMap::<String, usize>::new();
        for tag in &self.peak_live_tags {
            let class = Self::normalized_tag_class(tag);
            *tag_classes.entry(class).or_default() += 1;
        }
        let mut tag_classes = tag_classes.into_iter().collect::<Vec<_>>();
        tag_classes.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
        if std::env::var("DEBUG_ON_FAIL").is_ok() {
            eprintln!("{msg}");
            eprintln!("[debugger] peak_live_tag_class_histogram_top10:");
            for (class, count) in tag_classes.iter().take(10) {
                eprintln!("  {count:>5} {class}");
            }
            eprintln!("[debugger] DEBUG_ON_FAIL set — attaching");
            // See MAX_OPS site: capture target BEFORE attach (attach moves ops).
            let target = self.ops.len() + self.ops_truncated as usize;
            let mut d = crate::point_add::trailmix_port::debugger::Debugger::attach(self);
            d.goto(target);
            d.repl();
        } else {
            eprintln!("{msg}");
            eprintln!("[debugger] peak_live_tag_class_histogram_top10:");
            for (class, count) in tag_classes.iter().take(10) {
                eprintln!("  {count:>5} {class}");
            }
        }
        panic!("{msg}");
    }

    fn normalized_tag_class(tag: &str) -> String {
        let base = tag.rsplit('/').next().unwrap_or("(untagged)");
        let mut out = String::with_capacity(base.len());
        let mut chars = base.chars().peekable();
        while let Some(ch) = chars.next() {
            if ch == '[' {
                out.push('[');
                for next in chars.by_ref() {
                    if next == ']' {
                        out.push(']');
                        break;
                    }
                }
            } else {
                out.push(ch);
            }
        }
        out
    }

    /// Snapshot the live-qubit set at the current moment. Called by
    /// `alloc_qubit`/`alloc_qubit_fresh` whenever `peak_qubits` grows.
    /// Live = `{0 .. next_qubit} \ free_qubits` — the set of
    /// allocated-but-not-freed qubits. Used post-construction to
    /// break the peak down by section tag.
    fn snapshot_peak_live(&mut self) {
        let free: std::collections::HashSet<u32> = self.free_qubits.iter().copied().collect();
        self.peak_live_qubits.clear();
        self.peak_live_tags.clear();
        for q in 0..self.next_qubit {
            if free.contains(&q) {
                continue;
            }
            self.peak_live_qubits.push(q);
            self.peak_live_tags.push(
                self.qubit_tags
                    .get(q as usize)
                    .and_then(|t| t.as_ref())
                    .map_or_else(
                        || "(untagged)".to_string(),
                        std::string::ToString::to_string,
                    ),
            );
        }
    }

    /// Allocate `n` qubits tagged `name[0]`, `name[1]`, …, `name[n-1]`.
    #[allow(dead_code)]
    fn alloc_qubits(&mut self, name: &str, n: usize) -> Vec<Qubit> {
        (0..n)
            .map(|i| self.alloc_qubit(&format!("{name}[{i}]")))
            .collect()
    }

    /// Allocate a single-qubit `QReg`. The qubit is queued for free when
    /// the `QReg` drops (drain on next gate emission).
    pub fn alloc_qreg(&mut self, name: &str) -> QReg {
        let qubit = self.alloc_qubit(name);
        QReg {
            id: qubit.raw(),
            pending: std::rc::Rc::clone(&self.pending_frees),
            detached: false,
        }
    }

    /// Allocate `n` single-qubit `QRegs` as a `Vec<QReg>`. Each bit is
    /// independently owned and gets RAII drop. Use this for multi-bit
    /// registers — multi-bit `QReg` does not exist.
    pub fn alloc_qreg_bits(&mut self, name: &str, n: usize) -> Vec<QReg> {
        (0..n)
            .map(|i| {
                let qubit = self.alloc_qubit(&format!("{name}[{i}]"));
                QReg {
                    id: qubit.raw(),
                    pending: std::rc::Rc::clone(&self.pending_frees),
                    detached: false,
                }
            })
            .collect()
    }

    /// Re-tag an existing `Vec<QReg>` with a new semantic name. Each
    /// bit's tag becomes `{name}[{i}]` and a fresh alloc-log entry is
    /// recorded at the current op so the debugger and profiler see
    /// the rename at this point. The qubits' physical identity and
    /// allocator state are unchanged — this is the right tool when
    /// transferring ownership of a register from one logical role to
    /// another (e.g. forward EEA reusing `dx_ext` as `pair1_big` without
    /// a CX-move).
    pub fn relabel_qreg(&mut self, regs: &[QReg], name: &str) {
        for (i, reg) in regs.iter().enumerate() {
            self.record_qubit_alloc(reg.id, format!("{name}[{i}]"));
        }
    }

    /// Allocate a single qubit and wrap it in a `SharedQReg` (Rc-shared
    /// ownership). Use when the qubit needs to be shared into `mod_arith`
    /// row-builders or other state-machine storage trees that cannot
    /// hold lifetime-bound borrows.
    pub fn alloc_shared_qreg(&mut self, name: &str) -> SharedQReg {
        SharedQReg::new(self.alloc_qreg(name))
    }

    /// `n`-bit version of `alloc_shared_qreg`.
    pub fn alloc_shared_qreg_bits(&mut self, name: &str, n: usize) -> Vec<SharedQReg> {
        self.alloc_qreg_bits(name, n)
            .into_iter()
            .map(SharedQReg::new)
            .collect()
    }

    /// Allocate a single-qubit `QReg` whose bit is marked as an
    /// independent input variable for the phase lattice (a fresh F2
    /// atom). Use for input registers whose bits are independent F2
    /// variables.
    pub fn alloc_input_qreg(&mut self, name: &str) -> QReg {
        let qubit = self.alloc_qubit(name);
        self.phase.mark_input_qubit(qubit.raw());
        QReg {
            id: qubit.raw(),
            pending: std::rc::Rc::clone(&self.pending_frees),
            detached: false,
        }
    }

    /// Allocate `n` single-qubit `QRegs` as inputs (each bit is a fresh
    /// F2 atom in the phase lattice).
    pub fn alloc_input_qreg_bits(&mut self, name: &str, n: usize) -> Vec<QReg> {
        (0..n)
            .map(|i| {
                let qubit = self.alloc_qubit(&format!("{name}[{i}]"));
                self.phase.mark_input_qubit(qubit.raw());
                QReg {
                    id: qubit.raw(),
                    pending: std::rc::Rc::clone(&self.pending_frees),
                    detached: false,
                }
            })
            .collect()
    }

    /// Allocate `n` input `QRegs` AND seed the simulator with the given
    /// per-qubit 64-shot lanes. CLEAN init-time API; the only public
    /// way to set sim input state without going through the legacy
    /// `sim_load_reg_bytes_*` endpoints (which now assert init-time-only).
    ///
    /// `lanes[i]` is the 64-shot lane for input qubit i — bit `s` of
    /// `lanes[i]` is the value of qubit i in shot `s`.
    ///
    /// **The API DELIBERATELY does NOT broadcast a single classical
    /// value to all 64 shots.** If you want all shots to see the same
    /// input (e.g., for diagnostics), build the lanes explicitly via
    /// `replicate_classical_to_lanes(value)` — this makes the
    /// "I am only testing one input" intent visible at the call site.
    /// The default fuzz pattern is 64 different inputs, one per shot.
    pub fn alloc_input_qreg_bits_with_lanes(
        &mut self,
        name: &str,
        n: usize,
        lanes: &[u64],
    ) -> Vec<QReg> {
        assert!(
            self.ops.is_empty(),
            "alloc_input_qreg_bits_with_lanes: must be called before any \
             gate emission. self.ops.len() = {}",
            self.ops.len()
        );
        assert_eq!(
            lanes.len(),
            n,
            "alloc_input_qreg_bits_with_lanes: lanes.len()={} but n={}",
            lanes.len(),
            n
        );
        let regs = self.alloc_input_qreg_bits(name, n);
        // Direct lane copy via the internal sim-mutation helper (the
        // public sim_load endpoints are init-time-asserted; we are
        // init-time, but we have raw lanes not bytes).
        if let Some(sim) = self.sim.as_mut() {
            for (i, q) in regs.iter().enumerate() {
                let qr = q.id;
                if (qr as usize) >= sim.len() {
                    sim.resize(qr as usize + 256, 0);
                }
                sim[qr as usize] = lanes[i];
            }
        }
        regs
    }

    /// Drain the QReg-dropped queue, calling `free_qubit` for each
    /// pending qubit. Called automatically by `push_gate_op`; the
    /// public version lets the caller force a drain at logical
    /// boundary points (e.g. between sub-circuits).
    pub fn flush_pending_frees(&mut self) {
        self.drain_pending_frees();
    }

    fn drain_pending_frees(&mut self) {
        // `take` clears the queue first, so re-entrant frees during
        // release (e.g. tracker side-effects that drop another QReg)
        // don't loop indefinitely. free_qubit applies its full
        // strict-dealloc check (last-gate-touch gap=0 required) —
        // QReg drop doesn't relax that, it's still per-bit strict.
        let pending: Vec<u32> = self.pending_frees.borrow_mut().drain(..).collect();
        for qr in pending {
            self.free_qubit(Qubit(qr));
        }
    }
}

impl Drop for Circuit {
    fn drop(&mut self) {
        // Unresolved Ghost receipts are a phase leak: the HMR
        // obligation against the captured Cbit was never discharged.
        // Report BEFORE the QReg-leak check so the diagnostic points
        // at the HMR site (which the Ghost's own Drop also panics on,
        // but if the Ghost was `mem::forget`'d, this is the only
        // remaining signal).
        if !self.pending_ghosts.is_empty() && !std::thread::panicking() {
            let n = self.pending_ghosts.len();
            let first = &self.pending_ghosts[0];
            panic!(
                "Circuit dropped with {} unresolved Ghost(s); first \
                 is HMR at {} (op #{}, section '{}', anchor_id={}, \
                 bit b{}).",
                n,
                first.hmr_caller,
                first.hmr_op_idx,
                first.hmr_section,
                first.anchor_id,
                first.bit_raw,
            );
        }
        // QRegs hold strong Rcs to pending_frees. If any are still
        // alive when the Circuit is being dropped, that's a leak of
        // qubits whose free_qubit was never called. Use QReg::detach
        // to convert to plain Vec<Qubit> if you need to keep the
        // qubit IDs around past Circuit lifetime.
        let count = std::rc::Rc::strong_count(&self.pending_frees);
        if count > 1 {
            // count - 1 = number of live QRegs (Circuit holds 1).
            // Don't panic if we're already unwinding — that'd abort.
            assert!(
                std::thread::panicking(),
                "Circuit dropped while {} QReg(s) still alive — \
                 qubits leak (free_qubit was never called). \
                 Drop QRegs before the Circuit, or call \
                 QReg::detach to convert to plain Vec<Qubit>.",
                count - 1,
            );
        }
        // Drain any QRegs dropped during this scope (after the last
        // gate emission, before Circuit drop). The strict-dealloc
        // check still fires per-qubit; if a register's bits had
        // a gap, that panic surfaces here (not a silent leak).
        let pending_count = self.pending_frees.borrow().len();
        if pending_count > 0 && !std::thread::panicking() {
            let pending: Vec<u32> = self.pending_frees.borrow_mut().drain(..).collect();
            for qr in pending {
                self.free_qubit(Qubit(qr));
            }
        }
    }
}

impl Circuit {
    pub fn alloc_bit(&mut self) -> Cbit {
        let b = self.free_bits.pop().unwrap_or_else(|| {
            let b = self.next_bit;
            self.next_bit += 1;
            b
        });
        let live = self.next_bit - self.free_bits.len() as u32;
        if live > self.peak_bits {
            self.peak_bits = live;
        }
        self.phase.on_alloc_bit(b);
        // Anonymous allocs still log a `b{N}` tag so the debugger can
        // resolve raw `bN` references via cbit_id_from_tag and so
        // cbit_name_at_id has a non-fallback entry.
        let tag = format!("b{b}");
        self.record_bit_alloc(b, tag);
        Cbit(b)
    }

    pub fn alloc_bits(&mut self, n: u32) -> Vec<Cbit> {
        (0..n).map(|_| self.alloc_bit()).collect()
    }

    /// Allocate a classical bit with a descriptive name. Mirrors
    /// `alloc_qubit(name)` — records a `bit_alloc_log` entry under
    /// the current scope so the debugger can resolve the bit by
    /// label via `p creg <name>` or `p <name>`.
    pub fn alloc_bit_named(&mut self, name: &str) -> Cbit {
        let b = self.alloc_bit_raw_id();
        let tag = self.scoped_tag(name);
        self.record_bit_alloc(b, tag);
        Cbit(b)
    }

    /// Allocate `n` classical bits as a named batch. Each bit gets
    /// the tag `<name>[<i>]`, mirroring `alloc_qubits`.
    pub fn alloc_bits_named(&mut self, name: &str, n: u32) -> Vec<Cbit> {
        (0..n)
            .map(|i| {
                let b = self.alloc_bit_raw_id();
                let tag = self.scoped_tag(&format!("{name}[{i}]"));
                self.record_bit_alloc(b, tag);
                Cbit(b)
            })
            .collect()
    }

    /// Inner alloc step used by both `alloc_bit` and the named
    /// helpers. Pops a free id (or grows the pool), updates peak,
    /// and notifies the phase tracker. Caller is responsible for
    /// recording the tag via `record_bit_alloc`.
    fn alloc_bit_raw_id(&mut self) -> u32 {
        let b = self.free_bits.pop().unwrap_or_else(|| {
            let b = self.next_bit;
            self.next_bit += 1;
            b
        });
        let live = self.next_bit - self.free_bits.len() as u32;
        if live > self.peak_bits {
            self.peak_bits = live;
        }
        self.phase.on_alloc_bit(b);
        b
    }

    /// Allocate a qubit that represents an independent input variable
    /// for the phase lattice — a fresh F2 atom. Use for input
    /// registers whose bits are independent F2 variables.
    #[allow(dead_code)]
    fn alloc_input_qubit(&mut self, name: &str) -> Qubit {
        let q = self.alloc_qubit(name);
        self.phase.mark_input_qubit(q.raw());
        q
    }

    /// Allocate a classical bit with a fresh F2 atom.
    pub fn alloc_input_bit(&mut self) -> Cbit {
        let b = self.alloc_bit();
        self.phase.mark_input_bit(b.raw());
        b
    }

    #[allow(dead_code)]
    fn free_qubits_contains(&self, q: Qubit) -> bool {
        self.free_qubits.contains(&q.raw())
    }

    fn free_qubit(&mut self, q: Qubit) {
        let qr = q.raw();
        assert!(!self.free_qubits.contains(&qr), "double-free of qubit {qr}");
        // Wasteful-retention check: a qubit being freed must have
        // been touched by a gate AFTER the most recent allocation.
        // `None` means the qubit was allocated but never touched —
        // that's wasted ancilla, fix the alloc site. If
        // `last_touched_op[q] < last_alloc_op_idx`, the qubit was
        // idle across at least one allocation, meaning the allocator
        // had to grow when it could have reused this slot — fix the
        // code structure to free this qubit BEFORE that allocation.
        // Gaps with no intervening allocation are permitted (the
        // qubit's slot wouldn't have been reused regardless).
        {
            let tag = || {
                self.qubit_tags
                    .get(qr as usize)
                    .and_then(|t| t.as_ref())
                    .map_or_else(|| format!("q{qr}"), std::string::ToString::to_string)
            };
            let last_alloc = self.last_alloc_op_idx;
            let cur = self.phase.current_op_idx as u64;
            match self.last_touched_op.get(qr as usize).copied().flatten() {
                None => {
                    panic!(
                        "free_qubit q{qr} ({}): allocated but never \
                         touched by any gate. An untouched qubit is \
                         wasted ancilla — remove the alloc instead \
                         of freeing it.",
                        tag(),
                    );
                }
                Some(last) => {
                    if last < last_alloc {
                        // Identify the recent allocs whose op_idx is at or
                        // after `last_alloc` so the panic message points at
                        // the offending site, not just an op number.
                        let nearby = self
                            .qubit_alloc_log
                            .iter()
                            .rev()
                            .take_while(|ev| (ev.op_idx as u64) >= last)
                            .take(6)
                            .map(|ev| format!("@op{} q{}/{}", ev.op_idx, ev.qubit, ev.tag))
                            .collect::<Vec<_>>()
                            .join("; ");
                        let msg = format!(
                            "free_qubit q{qr} ({}): last gate-touch \
                             at op {last}, current op {cur}, but a \
                             later allocation happened at op \
                             {last_alloc}. The qubit was retained \
                             across that allocation — free it BEFORE \
                             the next alloc so the allocator can \
                             reuse its slot.\nRecent allocs after last touch: {}",
                            tag(),
                            nearby,
                        );
                        if std::env::var("DEBUG_ON_FAIL").is_ok() {
                            eprintln!("{msg}");
                            eprintln!(
                                "[debugger] DEBUG_ON_FAIL set — attaching at last touch op {last}"
                            );
                            let mut d = crate::point_add::trailmix_port::debugger::Debugger::attach(self);
                            d.goto(last as usize);
                            d.repl();
                        }
                        if std::env::var("TRAILMIX_BRIDGE_SKIP_ALLOC_AUDIT").is_err() {
                            panic!("{msg}");
                        }
                    }
                }
            }
        }
        // Sim-backed prove-zero, integrated. The simulator is always
        // live for the Circuit's lifetime, so prove_zero always runs:
        // assert q reads |0> on all 64 shots and inject Zero into the
        // tracker (DEBUG_ON_FAIL=1 drops into the REPL at the failing
        // op). The sim-disabled escape hatch was removed — it
        // weakened the strict-dealloc guarantee on a clean,
        // deterministic axis.
        self.prove_zero_raw(qr);
        self.phase.on_free_qubit(qr);
        self.free_qubits.insert(qr);
    }

    pub fn free_bit(&mut self, b: Cbit) {
        self.phase.on_free_bit(b.raw());
        self.free_bits.push(b.raw());
    }

    // Register declarations
    pub fn register(&mut self, id: u32) {
        *self
            .executed_ops_by_section
            .entry(self.current_section.clone())
            .or_insert(0) += 1;
        self.ops.push(Some(Op::Register(id)));
    }

    #[allow(dead_code)]
    fn append_qubit(&mut self, q: Qubit, reg: u32) {
        *self
            .executed_ops_by_section
            .entry(self.current_section.clone())
            .or_insert(0) += 1;
        self.ops.push(Some(Op::AppendQubit(q.raw(), reg)));
    }

    /// Append the single qubit underlying `q` to the output register
    /// `reg`. Use after `register(reg)` to declare which qubits hold
    /// the output bits of the program.
    pub fn append_qreg(&mut self, q: &QReg, reg: u32) {
        *self
            .executed_ops_by_section
            .entry(self.current_section.clone())
            .or_insert(0) += 1;
        self.ops.push(Some(Op::AppendQubit(q.id, reg)));
    }

    pub fn append_bit(&mut self, b: Cbit, reg: u32) {
        *self
            .executed_ops_by_section
            .entry(self.current_section.clone())
            .or_insert(0) += 1;
        self.ops.push(Some(Op::AppendBit(b.raw(), reg)));
    }

    // All gate sim semantics match the zenodo Simulator (lib/src/sim.rs):
    // all arithmetic is u64-bitwise, with a condition mask that gates
    // per-shot execution.

    // The gate entry points and their `*_internal` helpers are
    // `#[track_caller]` so the location captured in `push_gate_op`
    // skips these Circuit wrapper frames and resolves to the primitive
    // call site (cursor_pack_v4.rs, adders, …). `#[track_caller]` only
    // propagates through a frame if THAT frame is annotated, so every
    // link in  pub gate → *_internal → push_gate_op  must carry it.
    #[track_caller]
    pub fn x(&mut self, q: &QReg) {
        self.x_internal(Qubit(q.id));
    }
    #[track_caller]
    pub fn z(&mut self, q: &QReg) {
        self.z_internal(Qubit(q.id));
    }
    #[track_caller]
    pub fn cx(&mut self, ctrl: &QReg, tgt: &QReg) {
        self.cx_internal(Qubit(ctrl.id), Qubit(tgt.id));
    }
    #[track_caller]
    pub fn cz(&mut self, a: &QReg, b: &QReg) {
        self.cz_internal(Qubit(a.id), Qubit(b.id));
    }
    #[track_caller]
    pub fn ccx(&mut self, c1: &QReg, c2: &QReg, tgt: &QReg) {
        self.ccx_internal(Qubit(c1.id), Qubit(c2.id), Qubit(tgt.id));
    }
    #[track_caller]
    pub fn ccz(&mut self, a: &QReg, b: &QReg, c: &QReg) {
        self.ccz_internal(Qubit(a.id), Qubit(b.id), Qubit(c.id));
    }

    /// Marker for `emit_reverse_since`: current length of the retained
    /// op stream. Pass the value returned here to `emit_reverse_since`
    /// after running a forward block to emit its inverse.
    #[must_use]
    pub fn reverse_marker(&self) -> usize {
        self.ops.len()
    }

    /// Emit the inverse of a PURE-GATE block recorded since `start`
    /// (a prior `reverse_marker()`), uncomputing it.
    ///
    /// SUPPORTED: the self-inverse gates X/Z/Cx/Cz/Ccx/Ccz plus
    /// arbitrary internal alloc/dealloc. The block may freely allocate
    /// AND free ancillae (cinc, compare, rotate, bitlen, ...): a forward
    /// dealloc (R) is reversed to a fresh alloc of a |0> qubit (a
    /// dealloc is only legal on |0>, so its inverse is an alloc), and a
    /// forward alloc is reversed to a dealloc; an id-remap retargets the
    /// replayed gates so reused qubit ids resolve correctly. Inputs
    /// (allocated before the block) map identity and are untouched.
    ///
    /// So composite primitives can be generically reversed with NO
    /// hand-mirroring: run forward, then `emit_reverse_since(mark)`.
    ///
    /// NOT supported (panics): Hmr (irreversible measurement), Swap,
    /// Neg, Push/PopCondition, `BitStore`. Blocks using those need a
    /// hand-written inverse.
    #[track_caller]
    pub fn emit_reverse_since(&mut self, start: usize) {
        let block: Vec<Op> = self.ops[start..].iter().flatten().copied().collect();
        // Pre-scan: ids freed (R) within the block are "balanced" internal
        // ancillae — the reverse re-allocs them (at the R) and deallocs
        // them (at the AppendQubit). Ids ALLOCATED but NOT freed in the
        // block are "output" ancillae (returned live, e.g. NarrowDivAnc):
        // the reverse re-zeroes them via the inverse gates but leaves the
        // dealloc to the caller's QReg Drop — deallocing here would
        // double-free the still-live QReg wrapper.
        let freed_in_block: std::collections::HashSet<u32> = block
            .iter()
            .filter_map(|op| if let Op::R(q) = op { Some(*q) } else { None })
            .collect();
        // Pre-pair each PopCondition (by block index) with its matching
        // PushCondition's bit, so the reverse can re-emit the condition:
        // forward [Push(b); inner; Pop] reverses to [Push(b); rev-inner;
        // Pop] (forward Pop -> reverse Push(b), forward Push -> reverse
        // Pop). R/HMR never appear inside a condition (primitives use the
        // push_condition-safe fallback), so alloc-remap and conditions
        // don't interleave.
        let mut pop_bit: std::collections::HashMap<usize, u32> = std::collections::HashMap::new();
        {
            let mut pstack: Vec<u32> = Vec::new();
            for (i, op) in block.iter().enumerate() {
                match op {
                    Op::PushCondition(b) => pstack.push(*b),
                    Op::PopCondition => {
                        let b = pstack
                            .pop()
                            .expect("emit_reverse_since: unbalanced Push/PopCondition in block");
                        pop_bit.insert(i, b);
                    }
                    _ => {}
                }
            }
        }
        // id-remap: forward qubit id -> the reverse's current id for it.
        // A forward dealloc (R) reverses to a fresh alloc; a forward alloc
        // (AppendQubit) reverses to a dealloc. Gates retarget through the
        // map (identity for ids untouched by alloc/free in the block, e.g.
        // pre-allocated inputs).
        // id-remap as a per-id STACK: a forward id can be freed and re-allocated
        // many times in the block (e.g. a uniter's per-leaf ancilla reuses the
        // same id). Each forward R pushes a fresh reverse id; each forward alloc
        // pops one. A single-value map would OVERWRITE on the 2nd R and orphan
        // the prior fresh qubit (leaking ~one rev.anc per reuse).
        let mut remap: std::collections::HashMap<u32, Vec<u32>> = std::collections::HashMap::new();
        let map_id = |remap: &std::collections::HashMap<u32, Vec<u32>>, q: u32| -> u32 {
            remap.get(&q).and_then(|s| s.last()).copied().unwrap_or(q)
        };
        for (idx, op) in block.iter().enumerate().rev() {
            match *op {
                // forward Pop -> reverse Push(matched bit); forward Push
                // -> reverse Pop. Keeps the condition active around the
                // reversed inner ops.
                Op::PopCondition => {
                    let b = pop_bit[&idx];
                    self.push_gate_op(Op::PushCondition(b));
                    self.sim_condition_stack.push(b);
                    self.phase.on_push_condition(b);
                }
                Op::PushCondition(_) => {
                    self.push_gate_op(Op::PopCondition);
                    self.sim_condition_stack.pop();
                    self.phase.on_pop_condition();
                }
                Op::X(q) => self.x_internal(Qubit(map_id(&remap, q))),
                Op::Z(q) => self.z_internal(Qubit(map_id(&remap, q))),
                Op::Cx(c, t) => {
                    self.cx_internal(Qubit(map_id(&remap, c)), Qubit(map_id(&remap, t)));
                }
                Op::Cz(a, b) => {
                    self.cz_internal(Qubit(map_id(&remap, a)), Qubit(map_id(&remap, b)));
                }
                Op::Ccx(a, b, c) => self.ccx_internal(
                    Qubit(map_id(&remap, a)),
                    Qubit(map_id(&remap, b)),
                    Qubit(map_id(&remap, c)),
                ),
                Op::Ccz(a, b, c) => self.ccz_internal(
                    Qubit(map_id(&remap, a)),
                    Qubit(map_id(&remap, b)),
                    Qubit(map_id(&remap, c)),
                ),
                // forward dealloc -> reverse alloc of a fresh |0> qubit.
                // (A dealloc is only legal on |0>, so its inverse is an
                // alloc of a |0> qubit; remap retargets later gates.)
                Op::R(q) => {
                    let fresh = self.alloc_qubit("rev.anc").raw();
                    remap.entry(q).or_default().push(fresh);
                }
                // forward alloc -> reverse dealloc, but ONLY for balanced
                // internal ancillae (also freed in the block). Output
                // ancillae (alloc-only, returned live) are left for the
                // caller's QReg Drop — the inverse gates have already
                // re-zeroed them.
                Op::AppendQubit(q, _) => {
                    if freed_in_block.contains(&q) {
                        let cur = map_id(&remap, q);
                        self.prove_zero_raw(cur);
                        self.r_internal(Qubit(cur));
                        if let Some(s) = remap.get_mut(&q) {
                            s.pop();
                        }
                    }
                }
                Op::Register(_) | Op::AppendBit(_, _) => {}
                ref other => panic!(
                    "emit_reverse_since: unsupported op {other:?} in reverse block. \
                     Supported: X/Z/Cx/Cz/Ccx/Ccz (self-inverse), plus \
                     alloc/dealloc (handled via id-remap). NOT supported: \
                     Hmr (measurement), Swap, Neg, Push/PopCondition, \
                     BitStore — restructure to avoid these in a reversed block."
                ),
            }
        }
    }

    #[track_caller]
    fn x_internal(&mut self, q: Qubit) {
        let qr = q.raw();
        // Whether elided or not, we always run the sim/phase apply:
        // for elide it acts as the "second X" that undoes the prior
        // emission's first X (involutary self-inverse).
        let _elided = self.push_gate_op(Op::X(qr));
        if self.sim.is_some() {
            let cond = self.sim_condition_mask();
            self.sim_ensure(qr);
            self.sim.as_mut().unwrap()[qr as usize] ^= cond;
        }
        self.phase.on_x(qr);
    }

    #[track_caller]
    fn z_internal(&mut self, q: Qubit) {
        let qr = q.raw();
        let _elided = self.push_gate_op(Op::Z(qr));
        if self.sim.is_some() {
            let cond = self.sim_condition_mask();
            self.sim_ensure(qr);
            self.sim_phase ^= cond & self.sim.as_ref().unwrap()[qr as usize];
        }
        self.phase.on_z(qr);
    }

    #[track_caller]
    fn cx_internal(&mut self, ctrl: Qubit, tgt: Qubit) {
        let (cr, tr) = (ctrl.raw(), tgt.raw());
        let _elided = self.push_gate_op(Op::Cx(cr, tr));
        if self.sim.is_some() {
            let cond = self.sim_condition_mask();
            self.sim_ensure(cr);
            self.sim_ensure(tr);
            let c = self.sim.as_ref().unwrap()[cr as usize];
            self.sim.as_mut().unwrap()[tr as usize] ^= cond & c;
        }
        self.phase.on_cx(cr, tr);
    }

    #[track_caller]
    fn cz_internal(&mut self, a: Qubit, b: Qubit) {
        let (ar, br) = (a.raw(), b.raw());
        let _elided = self.push_gate_op(Op::Cz(ar, br));
        if self.sim.is_some() {
            let cond = self.sim_condition_mask();
            self.sim_ensure(ar);
            self.sim_ensure(br);
            let v = self.sim.as_ref().unwrap();
            self.sim_phase ^= cond & v[ar as usize] & v[br as usize];
        }
        self.phase.on_cz(ar, br);
    }

    #[track_caller]
    fn ccx_internal(&mut self, c1: Qubit, c2: Qubit, tgt: Qubit) {
        let (c1r, c2r, tr) = (c1.raw(), c2.raw(), tgt.raw());
        assert!(
            c1r != c2r && c1r != tr && c2r != tr,
            "non-physical CCX(q{c1r}, q{c2r}, q{tr}): all operands must be distinct"
        );
        self.push_gate_op(Op::Ccx(c1r, c2r, tr));
        if self.sim.is_some() {
            let cond = self.sim_condition_mask();
            self.sim_ensure(c1r);
            self.sim_ensure(c2r);
            self.sim_ensure(tr);
            let v = self.sim.as_ref().unwrap();
            let k = cond & v[c1r as usize] & v[c2r as usize];
            self.sim.as_mut().unwrap()[tr as usize] ^= k;
            let fired = u64::from(cond.count_ones());
            self.executed_toffoli_shots += fired;
            *self
                .executed_toffoli_by_section
                .entry(self.current_section.clone())
                .or_insert(0) += fired;
        }
        self.phase.on_ccx(c1r, c2r, tr);
    }

    #[track_caller]
    fn ccz_internal(&mut self, a: Qubit, b: Qubit, c: Qubit) {
        let (ar, br, cr) = (a.raw(), b.raw(), c.raw());
        self.push_gate_op(Op::Ccz(ar, br, cr));
        if self.sim.is_some() {
            let cond = self.sim_condition_mask();
            self.sim_ensure(ar);
            self.sim_ensure(br);
            self.sim_ensure(cr);
            let v = self.sim.as_ref().unwrap();
            let k = cond & v[ar as usize] & v[br as usize] & v[cr as usize];
            self.sim_phase ^= k;
            let fired = u64::from(cond.count_ones());
            self.executed_toffoli_shots += fired;
            *self
                .executed_toffoli_by_section
                .entry(self.current_section.clone())
                .or_insert(0) += fired;
        }
        self.phase.on_ccz(ar, br, cr);
    }

    /// Public &QReg-taking `declare_identity`. Forwards to the raw u32
    /// implementation in `phase_lattice.rs`.
    pub fn declare_identity(&mut self, q_a: &QReg, q_b: &QReg) {
        self.declare_identity_raw(q_a.id, q_b.id);
    }
    pub fn declare_copy_of(&mut self, q: &QReg, source: &QReg) {
        self.declare_copy_of_raw(q.id, source.id);
    }
    pub fn declare_and_of(&mut self, q: &QReg, a: &QReg, b: &QReg) {
        self.declare_and_of_raw(q.id, a.id, b.id);
    }
    pub fn declare_and3_of(&mut self, q: &QReg, a: &QReg, b: &QReg, c: &QReg) {
        self.declare_and3_of_raw(q.id, a.id, b.id, c.id);
    }
    /// `target := q1 XOR q2`. Sim-verified across 64 shots, then injected
    /// into the tracker as `XorOf(CopyOf(q1, v1), CopyOf(q2, v2))`. Used
    /// by snapshot-discharge primitives (e.g. `mod_mul_fused_v2`'s cma:phase
    /// replacement) where the HMR'd qubit's value is provably a 2-term XOR
    /// of two live qubits, and the discharge fires `z_if_bit(q1, bit);
    /// z_if_bit(q2, bit)` to deposit the matching `XorOf` leaves.
    pub fn declare_xor_of(&mut self, target: &QReg, q1: &QReg, q2: &QReg) {
        self.declare_xor_of_raw(target.id, q1.id, q2.id);
    }
    /// `target := q1 XOR q2 XOR q3`. Three-term XOR analogue of
    /// `declare_xor_of`. Sim-verified across 64 shots, then injected as
    /// `XorOf3(...)`. Use when the HMR'd qubit's value is a 3-term XOR
    /// of live qubits (e.g. snapshot-based discharge schemes where the
    /// algebraic identity involves a snapshot + an AND + a live data
    /// bit). Caller must ensure q1, q2, q3 are not modified between this
    /// declare and the corresponding `z_if_bit` discharges.
    pub fn declare_xor_of_three(&mut self, target: &QReg, q1: &QReg, q2: &QReg, q3: &QReg) {
        self.declare_xor_of_three_raw(target.id, q1.id, q2.id, q3.id);
    }

    /// Clear an AND-ancilla `t` (which currently holds `a AND b`, e.g.
    /// from a prior `ccx(a, b, t)`) back to |0>. Picks the uncompute
    /// strategy by inspecting the condition stack:
    ///
    ///   - no pushed condition: MBU clear — `declare_and_of` + HMR(t)
    ///     + `cz_if_bit(a, b)`. No Toffoli (the measurement replaces the
    ///     reverse CCX); this is the cheap path used everywhere outside
    ///     a `with_condition` block.
    ///   - inside a `push_condition` block: reversible clear
    ///     `ccx(a, b, t)`. HMR is forbidden under a condition, so pay
    ///     one Toffoli instead. `t` was `a&b`, so `t ^= a&b` → |0>.
    ///
    /// In both cases `t` is |0> on return. This does NOT free `t` — the
    /// caller still owns it; freeing only requires it already be |0>,
    /// which this guarantees.
    pub fn clear_and(&mut self, t: &QReg, a: &QReg, b: &QReg) {
        if self.sim_condition_stack.is_empty() {
            self.declare_and_of(t, a, b);
            let bit = self.alloc_bit();
            self.hmr(t, bit);
            self.cz_if_bit(a, b, bit);
            self.free_bit(bit);
        } else {
            self.ccx(a, b, t);
        }
    }
    /// `target := ctrl ? choice_a : choice_b`. Sim-verified across 64
    /// shots, then injected as `AbsVal::ChooseOf` into the tracker. Used
    /// by Luo-pack passthrough-shift primitives.
    pub fn declare_choose_of(
        &mut self,
        target: &QReg,
        ctrl: &QReg,
        choice_a: &QReg,
        choice_b: &QReg,
    ) {
        self.declare_choose_of_raw(target.id, ctrl.id, choice_a.id, choice_b.id);
    }

    /// Pure sim assertion: `q` reads classical `expected` (0 or 1)
    /// uniformly across all 64 shots. Does not modify the tracker or
    /// emit ops. `DEBUG_ON_FAIL=1` attaches the debugger on mismatch.
    pub fn assert_qubit_eq(&mut self, q: &QReg, expected: u8) {
        self.assert_qubit_eq_raw(q.id, expected);
    }

    /// Pure sim assertion: `qreg` (LSB at qreg[0]) encodes classical
    /// `expected` across all 64 shots.
    pub fn assert_reg_eq(&mut self, qreg: &[QReg], expected: u128) {
        for (i, q) in qreg.iter().enumerate() {
            let bit = ((expected >> i) & 1) as u8;
            self.assert_qubit_eq_raw(q.id, bit);
        }
    }

    /// Zero and free a qubit. Uses R (qubit must be |0⟩).
    // by-value IS the consumption: `q`'s Drop queues the free against the
    // pending-frees queue. Taking `&QReg` would skip the drop (qubit never freed).
    #[allow(clippy::needless_pass_by_value)]
    #[track_caller]
    pub fn zero_and_free(&mut self, q: QReg) {
        // Sim-backed: assert q is |0> across all 64 shots. Panics
        // with a clear message if the caller got here with a live
        // qubit; on success, tells the tracker q is Zero so the
        // subsequent R won't be flagged as R-on-nonzero. This is
        // the difference between a silent declare_zero (trust the
        // caller) and a runtime-verified assertion (prove via sim
        // evidence before injecting into the tracker).
        let qubit = Qubit(q.id);
        self.prove_zero_raw(qubit.raw());
        self.r_internal(qubit);
        // q drops here, queueing the free.
    }

    #[track_caller]
    pub fn swap(&mut self, a: &QReg, b: &QReg) {
        let (ar, br) = (a.id, b.id);
        self.push_gate_op(Op::Swap(ar, br));
        if self.sim.is_some() {
            // Per zenodo's 3-xor swap (conditional swap via cond mask):
            //   q_c1 ^= q_t
            //   q_t  ^= cond & q_c1
            //   q_c1 ^= q_t
            let cond = self.sim_condition_mask();
            self.sim_ensure(ar);
            self.sim_ensure(br);
            let s = self.sim.as_mut().unwrap();
            let mut qa = s[ar as usize];
            let mut qb = s[br as usize];
            qa ^= qb;
            qb ^= cond & qa;
            qa ^= qb;
            s[ar as usize] = qa;
            s[br as usize] = qb;
        }
        self.phase.on_swap(ar, br);
    }

    /// Compact + reorder live qubits into the contiguous id block
    /// `[0, n)` so that `slots[i]`'s VALUE ends up at qubit id `i`,
    /// emitting physical SWAPs. Returns fresh `QRegs` for ids `0..n-1`
    /// in the same order.
    ///
    /// `slots` MUST be the COMPLETE set of currently-live qubits
    /// (teardown-style); passing a partial set would let the swaps
    /// collide with un-passed live qubits and panics the all-live
    /// check below.
    ///
    /// Why: the zenodo fuzzer's register model writes the input AND
    /// reads the output from the SAME qubit ids. Internal free+realloc
    /// (e.g. `mod_mul_eea`'s GCD regrow) migrates a register's qubits
    /// to scattered high ids, so the emitted kmx's appended register
    /// (the END ids) no longer matches where the circuit read its
    /// input (the START ids 0..n). `defragment` restores the values to
    /// the canonical low ids, so input-ids == output-ids again.
    ///
    /// Two phases (relies on `alloc_qubit` handing out the lowest free
    /// id first — see `free_qubits: BTreeSet`):
    ///   1. Compaction: for each value at id `>= n`, alloc the lowest
    ///      free id (`< n`, since migrating freed a low slot), SWAP the
    ///      value in, free the old high id. After, the n values occupy
    ///      exactly `{0..n-1}`.
    ///   2. Sort: selection-swap so value `i` lands at id `i`.
    pub fn defragment(&mut self, mut slots: Vec<QReg>) -> Vec<QReg> {
        self.drain_pending_frees();
        let n = slots.len();
        if n == 0 {
            return slots;
        }
        let live = self.next_qubit as usize - self.free_qubits.len();
        assert_eq!(
            live, n,
            "defragment requires ALL live qubits (live={live}, passed={n})"
        );
        let n_u = n as u32;

        // Phase 1: move every value at id >= n into the lowest free id.
        let highs: Vec<usize> = (0..n).filter(|&v| slots[v].id >= n_u).collect();
        for v in highs {
            let lo = self.alloc_qreg("defrag");
            assert!(
                lo.id < n_u,
                "defragment phase 1: alloc gave id {} >= n {}",
                lo.id,
                n_u
            );
            self.swap(&lo, &slots[v]);
            let old = std::mem::replace(&mut slots[v], lo);
            self.zero_and_free(old);
        }

        // Phase 2: ids are now a permutation of 0..n-1; sort to identity.
        let mut pos: Vec<usize> = vec![0; n]; // pos[id] = value index at qubit `id`
        for v in 0..n {
            pos[slots[v].id as usize] = v;
        }
        for v in 0..n {
            let vu = v as u32;
            if slots[v].id == vu {
                continue;
            }
            let w = pos[v]; // value currently at id v
            let old_v_id = slots[v].id; // where value v currently lives (>= v)
            self.swap(&slots[v], &slots[w]);
            slots.swap(v, w); // keep slots[k] = the QReg now holding value k
            pos[v] = v;
            pos[old_v_id as usize] = w;
        }
        slots
    }

    /// Derive a 64-bit pseudo-random word from the HMR counter +
    /// seed (SplitMix-ish). Used for HMR and R ops so all 64 shots
    /// get independent random bits per op.
    fn sim_next_rng_u64(&mut self) -> u64 {
        self.sim_hmr_counter = self.sim_hmr_counter.wrapping_add(1);
        let mut z = self
            .sim_hmr_counter
            .wrapping_mul(0x9E37_79B9_7F4A_7C15)
            .wrapping_add(self.sim_hmr_seed);
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    #[track_caller]
    pub fn hmr(&mut self, q: &QReg, b: Cbit) {
        let (qr, br) = (q.id, b.raw());
        assert!(
            self.sim_condition_stack.is_empty(),
            "HMR q{qr} b{br} cannot be inside a push_condition block — \
             `HMR ... if cond` is forbidden.",
        );
        self.push_gate_op(Op::Hmr(qr, br));
        if self.sim.is_some() {
            let rng = self.sim_next_rng_u64();
            self.sim_ensure(qr);
            self.sim_bit_ensure(br);
            // Zenodo hmr: let r = rng & cond; bit = r; phase ^= qubit
            // & bit; qubit = 0. We're unconditional (asserted above),
            // so cond = u64::MAX.
            let qval = self.sim.as_ref().unwrap()[qr as usize];
            let r = rng;
            self.sim_bits.as_mut().unwrap()[br as usize] = r;
            self.sim_phase ^= qval & r;
            self.sim.as_mut().unwrap()[qr as usize] = 0;
        }
        // Symbolic F2 polynomial tracker: HMR introduces a fresh
        // atom for the random HMR bit, captures qubit-polynomial
        // XOR obligation on that bit, zeroes the qubit. Any bit that
        // gates a phase op (z_if_bit / cz_if_bit / ccz_if_bit) on
        // this HMR'd bit contributes back to its obligation; the
        // obligation must be the zero polynomial when the bit is
        // freed OR at assert_phase_clean.
        self.phase.on_hmr(qr, br);
    }

    /// Final build-time phase check. Panics if the global phase
    /// polynomial is non-zero OR if any HMR obligation has not been
    /// discharged. This is SOUND: F2-polynomial tracking over input
    /// atoms, cancellation via XOR, structural equality via `BTreeSet`
    /// of `BTreeSet` monomials. If the tracker returns zero, the phase
    /// is provably zero as an F2 polynomial in inputs. If a circuit
    /// uses a pattern the tracker can't reduce, the build fails with
    /// a list of uncancelled terms — restructure the circuit.
    pub fn assert_phase_clean(&self) {
        // Pre-check: pending Ghosts mean the caller skipped a
        // `resolve_ghost` somewhere. The tracker would also catch
        // the unmatched obligation, but the Ghost record carries
        // the HMR file:line which is the actually-useful diagnostic.
        if !self.pending_ghosts.is_empty() {
            let n = self.pending_ghosts.len();
            let first = &self.pending_ghosts[0];
            panic!(
                "assert_phase_clean: {} unresolved Ghost(s); first \
                 is HMR at {} (op #{}, section '{}', anchor_id={}, \
                 bit b{}). Resolve all Ghosts before asserting \
                 phase-clean.",
                n,
                first.hmr_caller,
                first.hmr_op_idx,
                first.hmr_section,
                first.anchor_id,
                first.bit_raw,
            );
        }
        self.phase.assert_clean();
    }

    /// Attach a time-travel debugger and enter its REPL over stdin.
    /// Blocks until the user quits. Call after the circuit is fully
    /// built. See `crate::point_add::trailmix_port::debugger::Debugger` for programmatic use.
    pub fn debug_repl(&mut self) {
        let mut d = crate::point_add::trailmix_port::debugger::Debugger::attach(self);
        d.repl();
    }

    // `prove_zero`, `declare_identity`, `declare_copy_of`, and
    // `declare_and_of` are defined as an `impl Circuit` block in
    // src/phase_lattice.rs so the tracker's `_unchecked` fact-
    // injection methods can stay module-private there. Each wrapper
    // requires an active sim, verifies the claimed identity across
    // all 64 shots, and panics on any mismatch — then forwards the
    // fact to the tracker. No path here to inject a fact blindly.

    // ── Contract / invariant checking ─────────────────────────────
    //
    // Primitives can annotate themselves with pre/post-conditions that
    // are verified at build time against the 64-shot sim.  Always on
    // when sim is active (no opt-out).  Any violated contract panics
    // immediately with (primitive label, shot index, current section,
    // op index, diagnostic).  The point: when a bug occurs, the FAILING
    // contract identifies the responsible primitive, not a downstream
    // symptom 300k ops later.
    //
    // Semantics:
    //  - `contracts_enabled()` is true iff sim is active.
    //  - `contract_read_u256_shot(reg, shot)` reads reg as little-endian
    //    unsigned integer for ONE shot.
    //  - `contract_read_signed_shot(reg, shot)` reads reg as two's-comp
    //    signed integer (top bit = sign).
    //  - `contract_check(label, F)` runs F over all 64 shots; F returns
    //    Result<(), String> where Err panics with shot + label.

    /// True iff contract checks are on. Cheap to call; callers can
    /// guard expensive per-shot reads behind it.
    pub(crate) fn contracts_enabled(&self) -> bool {
        self.sim.is_some()
    }

    fn contract_view(&self) -> ContractSimView<'_> {
        ContractSimView {
            qubits: self.sim.as_deref().unwrap_or(&[]),
            bits: self.sim_bits.as_deref().unwrap_or(&[]),
            #[cfg(test)]
            ghosts: &self.pending_ghosts,
        }
    }

    /// Run `check` over all 64 shots (`ops.len()` may be 0; ok). Each
    /// call returns Ok(()) or Err(detail). Panic loudly on first
    /// failure with label + shot + section + op index.
    pub(crate) fn contract_check<F>(&mut self, label: &str, mut check: F)
    where
        F: for<'a> FnMut(ContractSimView<'a>, usize) -> Result<(), String>,
    {
        if !self.contracts_enabled() {
            return;
        }
        let view = self.contract_view();
        for shot in 0..64 {
            if let Err(detail) = check(view, shot) {
                let msg = format!(
                    "CONTRACT [{}] shot {}: {} (section '{}', op #{})",
                    label,
                    shot,
                    detail,
                    self.current_section,
                    self.ops.len() + self.ops_truncated as usize
                );
                if std::env::var("DEBUG_ON_FAIL").is_ok() {
                    eprintln!("{msg}");
                    eprintln!("[debugger] DEBUG_ON_FAIL set — attaching");
                    // See MAX_OPS site: capture target BEFORE attach (attach moves ops).
                    let target = self.ops.len() + self.ops_truncated as usize;
                    let mut d = crate::point_add::trailmix_port::debugger::Debugger::attach(self);
                    d.goto(target);
                    d.repl();
                }
                panic!("{}", msg);
            }
        }
    }

    /// Capture per-shot pre-state. Pre runs inline (synchronously) and
    /// returns one T per shot; the Vec<T> is stored type-erased on the
    /// deferred-contract stack. Pop is `contract_pop_and_check::<T, _>`.
    ///
    /// Pre's closure may freely capture caller-scope refs (e.g. &[`QReg`]
    /// for the registers being read) — there is no `'static` bound,
    /// because the closure isn't stored, only its T outputs are.
    pub(crate) fn contract_capture<T, Pre>(&mut self, label: &str, mut pre: Pre)
    where
        T: 'static,
        Pre: for<'a> FnMut(ContractSimView<'a>, usize) -> Result<T, String>,
    {
        if !self.contracts_enabled() {
            return;
        }
        let view = self.contract_view();
        let mut captured: Vec<T> = Vec::with_capacity(64);
        for shot in 0..64 {
            match pre(view, shot) {
                Ok(value) => captured.push(value),
                Err(detail) => {
                    let msg = format!(
                        "CONTRACT [{} pre] shot {}: {} (section '{}', op #{})",
                        label,
                        shot,
                        detail,
                        self.current_section,
                        self.ops.len()
                    );
                    if std::env::var("DEBUG_ON_FAIL").is_ok() {
                        eprintln!("{msg}");
                        eprintln!("[debugger] DEBUG_ON_FAIL set — attaching");
                        // See MAX_OPS site: capture target BEFORE attach.
                        let target = self.ops.len();
                        let mut d = crate::point_add::trailmix_port::debugger::Debugger::attach(self);
                        d.goto(target);
                        d.repl();
                    }
                    panic!("{}", msg);
                }
            }
        }
        let id = self.capture_next_id;
        self.capture_next_id += 1;
        self.deferred_contracts.borrow_mut().push(DeferredContract {
            id,
            label: label.to_string(),
            captured: Box::new(captured),
        });
    }

    /// RAII-style capture handle: returns a Capture<T> token + pushes
    /// captured state onto the deferred contracts stack. Use
    /// `capture.update(circ, f)` and `capture.check(circ, f)` to re-check
    /// invariants at later points. Drop pops the entry automatically.
    pub(crate) fn contract_capture_handle<T, Pre>(&mut self, label: &str, pre: Pre) -> Capture<T>
    where
        T: 'static,
        Pre: for<'a> FnMut(ContractSimView<'a>, usize) -> Result<T, String>,
    {
        let id_pre = self.capture_next_id;
        self.contract_capture::<T, _>(label, pre);
        let id = if self.contracts_enabled() {
            self.capture_next_id - 1
        } else {
            id_pre
        };
        Capture {
            id,
            label: label.to_string(),
            stack: std::rc::Rc::clone(&self.deferred_contracts),
            _marker: std::marker::PhantomData,
        }
    }

    /// Read-only verify for a Capture<T> entry. Internal: called by `Capture::check`.
    pub(crate) fn capture_check_impl<T, F>(&self, cap: &Capture<T>, mut f: F)
    where
        T: 'static,
        F: for<'a> FnMut(&T, ContractSimView<'a>, usize) -> Result<(), String>,
    {
        if !self.contracts_enabled() {
            return;
        }
        let stack = self.deferred_contracts.borrow();
        let entry = stack
            .iter()
            .find(|d| d.id == cap.id)
            .unwrap_or_else(|| panic!("CONTRACT [{}]: capture id={} not found", cap.label, cap.id));
        let captured: &Vec<T> = entry.captured.downcast_ref::<Vec<T>>().unwrap_or_else(|| {
            panic!(
                "CONTRACT [{}]: type mismatch on capture_check_impl",
                entry.label
            )
        });
        let view = self.contract_view();
        for (shot, item) in captured.iter().enumerate() {
            if let Err(detail) = f(item, view, shot) {
                let msg = format!(
                    "CONTRACT [{} check] shot {}: {} (section '{}', op #{})",
                    entry.label,
                    shot,
                    detail,
                    self.current_section,
                    self.ops.len(),
                );
                panic!("{}", msg);
            }
        }
    }

    /// Pop the most recent deferred contract and run the caller-supplied
    /// post closure against the current sim state.
    ///
    /// The caller selects T (must match the T used in the matching
    /// `contract_capture` call) and passes a fresh closure that may
    /// freely borrow caller-scope `&[QReg]` refs — same lifetime story
    /// as `contract_capture`. The captured Vec<T> from the matching
    /// pre is downcast and threaded into the closure per shot.
    pub(crate) fn contract_pop_and_check<T, Post>(&mut self, label: &str, mut post: Post)
    where
        T: 'static,
        Post: for<'a> FnMut(&T, ContractSimView<'a>, usize) -> Result<(), String>,
    {
        if !self.contracts_enabled() {
            return;
        }
        let deferred = self
            .deferred_contracts
            .borrow_mut()
            .pop()
            .unwrap_or_else(|| panic!("CONTRACT [{label}]: empty deferred stack"));
        assert_eq!(
            deferred.label, label,
            "CONTRACT [{}]: deferred stack mismatch, found [{}]",
            label, deferred.label
        );
        let captured: Box<Vec<T>> = deferred.captured.downcast().unwrap_or_else(|_| {
            panic!(
                "CONTRACT [{label}]: captured-state type mismatch at pop. \
                 Either contract_capture and contract_pop_and_check used \
                 different T parameters, or stack frames are misaligned.",
            )
        });
        let view = self.contract_view();
        for (shot, item) in captured.iter().enumerate() {
            if let Err(detail) = post(item, view, shot) {
                let msg = format!(
                    "CONTRACT [{} post] shot {}: {} (section '{}', op #{})",
                    label,
                    shot,
                    detail,
                    self.current_section,
                    self.ops.len() + self.ops_truncated as usize
                );
                if std::env::var("DEBUG_ON_FAIL").is_ok() {
                    eprintln!("{msg}");
                    eprintln!("[debugger] DEBUG_ON_FAIL set — attaching");
                    // See MAX_OPS site: capture target BEFORE attach (attach moves ops).
                    let target = self.ops.len() + self.ops_truncated as usize;
                    let mut d = crate::point_add::trailmix_port::debugger::Debugger::attach(self);
                    d.goto(target);
                    d.repl();
                }
                panic!("{}", msg);
            }
        }
    }

    #[track_caller]
    pub fn neg(&mut self) {
        self.push_gate_op(Op::Neg);
        if self.sim.is_some() {
            let cond = self.sim_condition_mask();
            self.sim_phase ^= cond;
        }
        self.phase.on_neg();
    }

    #[track_caller]
    pub fn neg_if_bit(&mut self, bit: Cbit) {
        self.push_condition(bit);
        self.neg();
        self.pop_condition();
    }

    // Bit-conditioned gates. HMR bits are random, so gating
    // computational ops (X/CX/CCX) on them makes qubit VALUES
    // nondeterministic — UNLESS the gated block is Bennett
    // (compute + uncompute, net zero on all qubits), in which
    // case values are deterministic and only phase is affected.
    // Diagonal ops (Z/CZ/CCZ/NEG) are always safe to gate.
    // The table_lookup_3x3.kmx example uses CX-if-bit freely.
    #[track_caller]
    pub fn push_condition(&mut self, b: Cbit) {
        let br = b.raw();
        self.push_gate_op(Op::PushCondition(br));
        self.sim_condition_stack.push(br);
        self.phase.on_push_condition(br);
    }

    /// True iff at least one `push_condition` is active. Primitives
    /// that internally allocate + R-free ancillae (e.g. the
    /// Khattar–Gidney prefix-AND) consult this to fall back to a
    /// push_condition-safe construction, since `R` (and `HMR`) panic
    /// inside a conditional block.
    #[must_use]
    pub fn is_inside_push_condition(&self) -> bool {
        !self.sim_condition_stack.is_empty()
    }

    #[track_caller]
    pub fn pop_condition(&mut self) {
        self.push_gate_op(Op::PopCondition);
        self.sim_condition_stack.pop();
        self.phase.on_pop_condition();
    }

    /// Run a closure with the classical bit `bit` pushed onto the
    /// condition stack — every op emitted inside `f` fires iff
    /// `bit` is 1 at the push moment (kickmix's `PUSH_CONDITION`
    /// semantics are frozen-value, see the docs on line 530-ish of
    /// the instruction set spec).
    ///
    /// Preferred over manually bracketing `push_condition` / op /
    /// `pop_condition`: the closure scope makes the push/pop pairing
    /// syntactically impossible to miss, even across early returns
    /// or panics.
    pub fn with_condition<R>(&mut self, bit: Cbit, f: impl FnOnce(&mut Self) -> R) -> R {
        self.push_condition(bit);
        let result = f(self);
        self.pop_condition();
        result
    }

    /// Like `with_condition` but pushes a chain of classical bits
    /// (AND'd together as the effective condition). Pops them all on
    /// exit in reverse order to maintain push/pop nesting.
    pub fn with_conditions<R>(&mut self, bits: &[Cbit], f: impl FnOnce(&mut Self) -> R) -> R {
        for &b in bits {
            self.push_condition(b);
        }
        let result = f(self);
        for _ in bits {
            self.pop_condition();
        }
        result
    }

    // --- single-op classical-bit-conditioned shims ---
    // These are thin wrappers over `with_condition` for the
    // one-gate-per-condition case, which is common enough to deserve
    // the ergonomic name.

    // These are inlined (push_condition / gate / pop_condition) rather
    // than wrapped in `with_condition(|c| ...)` so `#[track_caller]`
    // propagates: a `#[track_caller]` annotation does NOT carry through
    // a closure invocation, so the closure form would make the debugger
    // `src` of a conditioned gate point at THIS file instead of the
    // primitive that called the shim. The inline form keeps the chain
    //   primitive → *_if_bit → *_internal → push_gate_op  all annotated.
    #[track_caller]
    pub fn x_if_bit(&mut self, q: &QReg, b: Cbit) {
        let qid = q.id;
        self.push_condition(b);
        self.x_internal(Qubit(qid));
        self.pop_condition();
    }
    #[track_caller]
    pub fn z_if_bit(&mut self, q: &QReg, bit: Cbit) {
        let qid = q.id;
        self.push_condition(bit);
        self.z_internal(Qubit(qid));
        self.pop_condition();
    }
    #[track_caller]
    pub fn cx_if_bit(&mut self, ctrl: &QReg, tgt: &QReg, bit: Cbit) {
        let (cid, tid) = (ctrl.id, tgt.id);
        self.push_condition(bit);
        self.cx_internal(Qubit(cid), Qubit(tid));
        self.pop_condition();
    }
    #[track_caller]
    pub fn cz_if_bit(&mut self, a: &QReg, b_qubit: &QReg, bit: Cbit) {
        let (aid, bid) = (a.id, b_qubit.id);
        self.push_condition(bit);
        self.cz_internal(Qubit(aid), Qubit(bid));
        self.pop_condition();
    }
    #[track_caller]
    pub fn ccx_if_bit(&mut self, c1: &QReg, c2: &QReg, tgt: &QReg, bit: Cbit) {
        let (c1id, c2id, tid) = (c1.id, c2.id, tgt.id);
        self.push_condition(bit);
        self.ccx_internal(Qubit(c1id), Qubit(c2id), Qubit(tid));
        self.pop_condition();
    }
    #[track_caller]
    pub fn ccz_if_bit(&mut self, a: &QReg, b: &QReg, c: &QReg, bit: Cbit) {
        let (aid, bid, cid) = (a.id, b.id, c.id);
        self.push_condition(bit);
        self.ccz_internal(Qubit(aid), Qubit(bid), Qubit(cid));
        self.pop_condition();
    }

    /// Write the bit pattern of a little-endian classical byte
    /// constant into a classical-bit register. Set bits get
    /// `BIT_STORE1`; unset bits stay |0⟩ from alloc. The register
    /// can then be used as a read-only "classical operand" in
    /// primitives like `add_cbits_mbu` without costing a 257-qubit
    /// scratch. Caller is responsible for clearing (via
    /// `clear_cbits_const`) and freeing before reuse.
    pub fn store_cbits_from_const(&mut self, bits: &[Cbit], val: &[u8]) {
        for (i, &b) in bits.iter().enumerate() {
            let bit = if i / 8 < val.len() {
                (val[i / 8] >> (i % 8)) & 1
            } else {
                0
            };
            if bit == 1 {
                self.bit_store1(b);
            }
        }
    }

    /// Inverse of `store_cbits_from_const`: clear set bits back to
    /// |0⟩ so the classical bits are ready for free/reuse.
    pub fn clear_cbits_const(&mut self, bits: &[Cbit], val: &[u8]) {
        for (i, &b) in bits.iter().enumerate() {
            let bit = if i / 8 < val.len() {
                (val[i / 8] >> (i % 8)) & 1
            } else {
                0
            };
            if bit == 1 {
                self.bit_store0(b);
            }
        }
    }

    /// Push a nested section `{current}/{sub}` and return the prior
    /// section name so the caller can restore it via `pop_section`.
    /// Used by primitives to get per-primitive op counts in profile
    /// output without the caller having to care.
    pub fn push_section(&mut self, sub: &str) -> String {
        let prev = self.current_section.clone();
        self.set_section(&format!("{prev}/{sub}"));
        prev
    }

    pub fn pop_section(&mut self, prev: &str) {
        self.set_section(prev);
    }

    pub fn set_section(&mut self, s: &str) {
        // If PHASE_TRACE is set, assert sim_phase is zero at every
        // section boundary. Reports which of the 64 shots broke
        // (bit-mask) so the caller can trace back the offending HMR
        // scenario.
        if std::env::var("PHASE_TRACE").is_ok() && self.sim.is_some() && self.sim_phase != 0 {
            let failing: Vec<usize> = (0..64).filter(|i| (self.sim_phase >> i) & 1 != 0).collect();
            eprintln!(
                "[PHASE_TRACE] sim_phase = {:#x} \
                (shots with bad phase: {:?}) entering section '{}' \
                (leaving '{}')",
                self.sim_phase, failing, s, self.current_section
            );
            panic!("sim_phase broken before section {s}");
        }
        // First set_section usually runs right before the first gate,
        // after all sim_load_* inputs are in place; snapshot for the
        // debugger's replay starting state.
        self.maybe_snapshot_initial();
        self.section_marks
            .push((self.ops.len() + self.ops_truncated as usize, s.to_string()));
        self.current_section = s.to_string();
        // Record entry-time live occupancy for this section so leaves that
        // never allocate (cswap, vents=0 adders) still get a local peak.
        let live = self.next_qubit - self.free_qubits.len() as u32;
        let sp = self.section_peak.entry(s.to_string()).or_insert(0);
        if live > *sp {
            *sp = live;
        }
        self.phase.current_section = s.to_string();
        // Reset per-section qubit-tag counter so the next alloc in
        // this section starts at q0 (scope-scoped naming).
        self.qubit_tag_counter = 0;
    }

    // `self` is read under cfg(test) (self.sim); not an associated fn.
    #[allow(dead_code, clippy::unused_self)]
    fn trace(&self, label: &str, ids: &[Qubit]) {
        #[cfg(test)]
        if self.sim.is_some() {
            eprintln!(
                "[TRACE] {:40} = <live-sim-readback-disabled> ({} bits)",
                label,
                ids.len()
            );
        }
        #[cfg(not(test))]
        let _ = (label, ids);
    }

    pub fn phase_summary(&self) {
        if self.sim_phase_errors > 0 {
            eprintln!(
                "  [PHASE] TOTAL: {} R-on-nonzero errors",
                self.sim_phase_errors
            );
            let mut tags: Vec<_> = self.sim_phase_errors_by_tag.iter().collect();
            tags.sort_by(|a, b| b.1.cmp(a.1));
            for (tag, count) in tags {
                eprintln!("    {tag:40} : {count}");
            }
        } else {
            eprintln!("  [PHASE] clean: 0 R-on-nonzero errors");
        }
    }

    /// Internal R primitive (resets a qubit by raw id). Used by drop
    /// machinery and `zero_and_free`; external callers should use
    /// `zero_and_free(q: QReg)` instead.
    #[track_caller]
    fn r_internal(&mut self, q: Qubit) {
        let qr = q.raw();
        assert!(
            self.sim_condition_stack.is_empty(),
            "R q{qr} cannot be inside a push_condition block — \
             `R ... if cond` is forbidden. Pop all conditions \
             before calling r().",
        );
        if self.sim.is_some() {
            self.sim_ensure(qr);
            let cond = self.sim_condition_mask();
            let rng = self.sim_next_rng_u64();
            let qval = self.sim.as_ref().unwrap()[qr as usize];
            // Zenodo R semantics: phase ^= qubit & rng & cond;
            // qubit = 0. On a truly-zero qubit this is a no-op; on a
            // non-zero qubit it contributes a random phase kick,
            // which is exactly what makes soundness testable — the
            // kick is retained in sim_phase so a bad R is detected.
            self.sim_phase ^= qval & rng & cond;
            if qval != 0 {
                self.sim_phase_errors += 1;
                let tag = self.current_section.clone();
                *self.sim_phase_errors_by_tag.entry(tag.clone()).or_insert(0) += 1;
            }
            self.sim.as_mut().unwrap()[qr as usize] = 0;
        }
        self.push_gate_op(Op::R(qr));
        self.phase.on_r(qr);
    }

    /// Assert no phase errors have occurred up to this point.
    /// Prints accumulated count and panics if any found.
    pub fn phase_assert(&self, label: &str) {
        assert!(
            self.sim_phase_errors == 0,
            "[PHASE ASSERT] {} phase errors at '{}' (op #{})",
            self.sim_phase_errors,
            label,
            self.ops.len() + self.ops_truncated as usize
        );
    }

    /// Assert the accumulated sim phase is 0 at this point.
    /// Panics if any bare HMR or mis-corrected MBUC has flipped
    /// the phase under the deterministic bit=qval model. This
    /// catches phase leaks at CIRCUIT GENERATION time, at the
    /// exact checkpoint they first break, rather than at the
    /// fuzz-test's global phase check at end-of-circuit.
    pub fn phase_assert_zero(&self, label: &str) {
        assert!(
            !(self.sim.is_some() && self.sim_phase != 0),
            "[PHASE ASSERT ZERO] sim_phase = {} at '{}' \
            (op #{}). A bare HMR or miscorrected MBUC flipped \
            the phase before this checkpoint.",
            self.sim_phase,
            label,
            self.ops.len() + self.ops_truncated as usize
        );
    }

    /// Sim-only: pre-load a classical byte value into consecutive
    /// qubits (LSB first). Used to initialise input registers with
    /// a specific test vector so that HMR/Z/CZ gate-tracking sees
    /// input-dependent phase kicks during circuit generation.
    /// Emits NO op — this just sets internal sim state. Must be
    /// called after qubits are allocated and before any ops that
    /// use them.
    /// Broadcast the same byte value into `reg`, replicating across
    /// all 64 shots. Use `sim_load_reg_bytes_shot` to load different
    /// values per shot.
    pub fn sim_load_reg_bytes(&mut self, reg: &[QReg], bytes: &[u8]) {
        // Hardening: enforce init-time discipline. Any gate op after
        // this would be a sign of mid-circuit sim mutation (a cheat
        // path used by previous agents to bake answers/witnesses).
        // Tests using `sim_load_reg_bytes_shot` for per-shot variation
        // call THAT directly with the same assertion.
        assert!(
            self.ops.is_empty(),
            "sim_load_reg_bytes called after {} gate ops were emitted. \
             This API is for INITIAL input loading only. Mid-circuit \
             sim mutation is forbidden — use the time-travel debugger \
             (DEBUG_ON_FAIL=1) for tracing instead.",
            self.ops.len()
        );
        for shot in 0..64 {
            self.sim_load_reg_bytes_shot_internal(reg, bytes, shot);
        }
    }

    /// Load bytes into `reg` for ONE shot only (leaves other shots
    /// untouched). Use to give each shot its own input.
    pub fn sim_load_reg_bytes_shot(&mut self, reg: &[QReg], bytes: &[u8], shot: usize) {
        // Same hardening as sim_load_reg_bytes — init-time only.
        assert!(
            self.ops.is_empty(),
            "sim_load_reg_bytes_shot called after {} gate ops were emitted. \
             This API is for INITIAL input loading only. Mid-circuit \
             sim mutation is forbidden — use the time-travel debugger \
             (DEBUG_ON_FAIL=1) for tracing instead.",
            self.ops.len()
        );
        self.sim_load_reg_bytes_shot_internal(reg, bytes, shot);
    }

    /// Internal helper: actual sim mutation, no init-time check.
    /// Used by both public entry points and by other init-time helpers.
    fn sim_load_reg_bytes_shot_internal(&mut self, reg: &[QReg], bytes: &[u8], shot: usize) {
        if self.sim.is_none() {
            return;
        }
        assert!(shot < 64, "shot must be < 64");
        let mask = 1u64 << shot;
        for (i, q) in reg.iter().enumerate() {
            let qr = q.id;
            self.sim_ensure(qr);
            let byte_idx = i / 8;
            let bit_idx = i % 8;
            let bit = if byte_idx < bytes.len() {
                (bytes[byte_idx] >> bit_idx) & 1
            } else {
                0
            };
            let cell = &mut self.sim.as_mut().unwrap()[qr as usize];
            if bit == 1 {
                *cell |= mask;
            } else {
                *cell &= !mask;
            }
        }
    }

    /// Broadcast bit values to all 64 shots.
    pub fn sim_load_bits_bytes(&mut self, bits: &[Cbit], bytes: &[u8]) {
        assert!(
            self.ops.is_empty(),
            "sim_load_bits_bytes called after {} gate ops were emitted. \
             Init-time only — no mid-circuit sim mutation.",
            self.ops.len()
        );
        for shot in 0..64 {
            self.sim_load_bits_bytes_shot_internal(bits, bytes, shot);
        }
    }

    /// Load bit values for ONE shot.
    pub fn sim_load_bits_bytes_shot(&mut self, bits: &[Cbit], bytes: &[u8], shot: usize) {
        assert!(
            self.ops.is_empty(),
            "sim_load_bits_bytes_shot called after {} gate ops were emitted. \
             Init-time only — no mid-circuit sim mutation.",
            self.ops.len()
        );
        self.sim_load_bits_bytes_shot_internal(bits, bytes, shot);
    }

    /// Internal helper, no init-time check.
    fn sim_load_bits_bytes_shot_internal(&mut self, bits: &[Cbit], bytes: &[u8], shot: usize) {
        if self.sim_bits.is_none() {
            return;
        }
        assert!(shot < 64, "shot must be < 64");
        let mask = 1u64 << shot;
        for (i, &b) in bits.iter().enumerate() {
            let br = b.raw();
            self.sim_bit_ensure(br);
            let byte_idx = i / 8;
            let bit_idx = i % 8;
            let bit = if byte_idx < bytes.len() {
                (bytes[byte_idx] >> bit_idx) & 1
            } else {
                0
            };
            let cell = &mut self.sim_bits.as_mut().unwrap()[br as usize];
            if bit == 1 {
                *cell |= mask;
            } else {
                *cell &= !mask;
            }
        }
    }

    // Classical bit operations (64-parallel).
    #[track_caller]
    pub fn bit_invert(&mut self, b: Cbit) {
        let br = b.raw();
        self.push_gate_op(Op::BitInvert(br));
        if self.sim_bits.is_some() {
            self.sim_bit_ensure(br);
            let cond = self.sim_condition_mask();
            self.sim_bits.as_mut().unwrap()[br as usize] ^= cond;
        }
        self.phase.on_bit_invert(br);
    }

    #[track_caller]
    pub fn bit_store0(&mut self, b: Cbit) {
        let br = b.raw();
        self.push_gate_op(Op::BitStore0(br));
        if self.sim_bits.is_some() {
            self.sim_bit_ensure(br);
            let cond = self.sim_condition_mask();
            self.sim_bits.as_mut().unwrap()[br as usize] &= !cond;
        }
        self.phase.on_bit_store0(br);
    }

    #[track_caller]
    pub fn bit_store1(&mut self, b: Cbit) {
        let br = b.raw();
        self.push_gate_op(Op::BitStore1(br));
        if self.sim_bits.is_some() {
            self.sim_bit_ensure(br);
            let cond = self.sim_condition_mask();
            self.sim_bits.as_mut().unwrap()[br as usize] |= cond;
        }
        self.phase.on_bit_store1(br);
    }

    // --- Compound classical-only ops ---
    //
    // Every op below is a short sequence of BIT_INVERT / BIT_STORE{0,1}
    // plus PUSH_CONDITION / POP_CONDITION nesting. They exist as
    // convenience wrappers so callers don't have to re-derive the
    // conditional gadget for common patterns (XOR, AND, OR, copy,
    // swap). Each is sim-verified; see circuit.rs tests.

    /// `dst ^= src` (classical bits). Equivalent to
    /// `BIT_INVERT dst if src`.
    pub fn bit_xor_into(&mut self, dst: Cbit, src: Cbit) {
        self.with_condition(src, |c| c.bit_invert(dst));
    }

    /// `dst = src` (classical bits), regardless of prior dst value.
    /// Implementation: `BIT_STORE0` dst; `BIT_STORE1` dst if src.
    pub fn bit_copy(&mut self, dst: Cbit, src: Cbit) {
        self.bit_store0(dst);
        self.with_condition(src, |c| c.bit_store1(dst));
    }

    /// `dst |= src`. Equivalent to `BIT_STORE1 dst if src`.
    pub fn bit_or_into(&mut self, dst: Cbit, src: Cbit) {
        self.with_condition(src, |c| c.bit_store1(dst));
    }

    /// `dst &= src`. Equivalent to `BIT_STORE0 dst if !src`, which
    /// expands to `BIT_INVERT src; BIT_STORE0 dst if src;
    /// BIT_INVERT src`. Requires `src` to be writable for the brief
    /// inversion window — callers using an input-only `src` should
    /// copy first.
    pub fn bit_and_into(&mut self, dst: Cbit, src: Cbit) {
        self.bit_invert(src);
        self.with_condition(src, |c| c.bit_store0(dst));
        self.bit_invert(src);
    }

    /// Swap two classical bits via the standard 3-XOR trick.
    pub fn bit_swap(&mut self, a: Cbit, b: Cbit) {
        self.bit_xor_into(a, b);
        self.bit_xor_into(b, a);
        self.bit_xor_into(a, b);
    }

    // Zero check: result = (all qubits in reg are 0).
    // Uses a multi-CX fan-in: OR all bits into result, then invert.
    // result must start at |0>. After: result = 1 if all bits are 0.
    // anc must be empty (n-2 qubits for n-qubit check).
    // All inputs (reg) are unchanged.
    pub fn zero_check(&mut self, reg: &[QReg], result: &QReg, anc: &[QReg]) {
        let n = reg.len();
        if n == 0 {
            self.x(result); // empty register is trivially zero
            return;
        }
        if n == 1 {
            self.cx(&reg[0], result);
            self.x(result);
            return;
        }
        // NOT all bits
        for q in reg {
            self.x(q);
        }
        if n == 2 {
            self.ccx(&reg[0], &reg[1], result);
        } else {
            assert!(anc.len() >= n - 2, "zero_check needs n-2 ancillae");
            self.ccx(&reg[0], &reg[1], &anc[0]);
            for i in 2..n - 1 {
                self.ccx(&anc[i - 2], &reg[i], &anc[i - 1]);
            }
            self.ccx(&anc[n - 3], &reg[n - 1], result);
            for i in (2..n - 1).rev() {
                self.ccx(&anc[i - 2], &reg[i], &anc[i - 1]);
            }
            self.ccx(&reg[0], &reg[1], &anc[0]);
        }
        for q in reg {
            self.x(q);
        }
    }

    // Controlled SWAP: if ctrl=1, swap a and b.
    // Uses CX+CCX decomposition (1 CCX + 2 CX).
    pub fn cswap(&mut self, ctrl: &QReg, a: &QReg, b: &QReg) {
        self.cx(b, a);
        self.ccx(ctrl, a, b);
        self.cx(b, a);
    }

    pub fn cond_right_shift(&mut self, ctrl: &QReg, reg: &[QReg]) {
        let n = reg.len();
        for i in 0..n - 1 {
            self.cswap(ctrl, &reg[i], &reg[i + 1]);
        }
    }

    pub fn cond_left_shift(&mut self, ctrl: &QReg, reg: &[QReg]) {
        let n = reg.len();
        for i in (1..n).rev() {
            self.cswap(ctrl, &reg[i], &reg[i - 1]);
        }
    }

    // Checkpoint for gate reversal
    #[must_use]
    pub fn checkpoint(&self) -> usize {
        self.ops.len()
    }

    #[must_use]
    pub fn op_count(&self) -> usize {
        self.ops
            .iter()
            .filter_map(|op| op.as_ref())
            .filter(|op| op.is_gate())
            .count()
    }

    #[must_use]
    pub fn ccx_count(&self) -> usize {
        self.ops
            .iter()
            .filter(|op| matches!(op, Some(Op::Ccx(_, _, _))))
            .count()
    }

    /// Count CCX patterns. Any non-real CCX (selfwire or overlap)
    /// is a bug — all CCXs must have three distinct qubits for
    /// physical realisability. Returns (selfwire, overlap, real).
    #[must_use]
    pub fn ccx_breakdown(&self) -> (usize, usize, usize) {
        let mut sw = 0usize;
        let mut overlap = 0usize;
        let mut real = 0usize;
        for op in self.ops.iter().filter_map(|op| op.as_ref()) {
            if let Op::Ccx(a, b, c) = *op {
                if a == b && b == c {
                    sw += 1;
                } else if b == c || a == c || a == b {
                    overlap += 1;
                } else {
                    real += 1;
                }
            }
        }
        (sw, overlap, real)
    }

    #[must_use]
    pub fn to_kmx(&self) -> String {
        // The op buffer is only the RETAINED tail: streaming truncation
        // (CIRC_OPS_CAP) `clear()`s earlier ops, counting them in
        // ops_truncated. Serializing only `self.ops` would silently emit
        // a kmx MISSING the truncated prefix — a wrong circuit. Fail loud
        // instead. (For >cap-op circuits, emit incrementally during build
        // before truncation clears the buffer.)
        assert_eq!(
            self.ops_truncated,
            0,
            "to_kmx: {} ops were streamed out of the retained buffer \
             (CIRC_OPS_CAP={}); the emitted kmx would be MISSING them. \
             Raise CIRC_OPS_CAP above total_ops()={} or stream the kmx \
             during build.",
            self.ops_truncated,
            self.ops_cap,
            self.total_ops(),
        );
        let mut out = String::new();
        for op in self.ops.iter().filter_map(|op| op.as_ref()) {
            writeln!(out, "{}", op.kmx_string()).unwrap();
        }
        out
    }

    /// Stream kmx lines into any `io::Write`. Preferred over `to_kmx`
    /// for large circuits — avoids materializing gigabytes of text
    /// in a single String.
    ///
    /// Collapses `PUSH_CONDITION if b<bit>` + single conditional-able op
    /// + `POP_CONDITION` into the kickmix inline form `<op> ... if b<bit>`,
    /// which is one kmx instruction instead of three. The in-memory op
    /// stream is unchanged; only the emitted kmx text is collapsed.
    ///
    /// # Errors
    /// Returns any `std::io::Error` produced while writing to `w`.
    pub fn write_kmx<W: std::io::Write>(&self, mut w: W) -> std::io::Result<()> {
        // See to_kmx: a truncated op buffer would emit a kmx missing the
        // streamed-out prefix. Fail loud rather than serialize a wrong
        // circuit.
        assert_eq!(
            self.ops_truncated,
            0,
            "write_kmx: {} ops were streamed out of the retained buffer \
             (CIRC_OPS_CAP={}); the emitted kmx would be MISSING them. \
             Raise CIRC_OPS_CAP above total_ops()={} or stream the kmx \
             during build.",
            self.ops_truncated,
            self.ops_cap,
            self.total_ops(),
        );
        let ops: Vec<&Op> = self.ops.iter().filter_map(|op| op.as_ref()).collect();
        let mut i = 0;
        while i < ops.len() {
            if i + 2 < ops.len() {
                if let (Op::PushCondition(b), inner, Op::PopCondition) =
                    (ops[i], ops[i + 1], ops[i + 2])
                {
                    if let Some(line) = inline_conditional(inner, *b) {
                        writeln!(w, "{line}")?;
                        i += 3;
                        continue;
                    }
                }
            }
            writeln!(w, "{}", ops[i].kmx_string())?;
            i += 1;
        }
        Ok(())
    }

    #[must_use]
    pub fn total_qubits(&self) -> u32 {
        self.next_qubit
    }

    #[must_use]
    pub fn live_qubits(&self) -> u32 {
        self.next_qubit - self.free_qubits.len() as u32
    }

    #[must_use]
    pub fn total_bits(&self) -> u32 {
        self.next_bit
    }
}

#[cfg(test)]
mod kmx_inline_conditional_tests {
    use super::*;

    fn write_kmx_to_string(c: &Circuit) -> String {
        let mut buf: Vec<u8> = Vec::new();
        c.write_kmx(&mut buf).unwrap();
        String::from_utf8(buf).unwrap()
    }

    #[test]
    fn cz_if_bit_collapses_to_inline() {
        let mut c = Circuit::new();
        let a = c.alloc_qreg("a");
        let b = c.alloc_qreg("b");
        let bit = c.alloc_bit();
        c.cz_if_bit(&a, &b, bit);

        let kmx = write_kmx_to_string(&c);
        assert!(
            kmx.lines()
                .any(|l| l.starts_with("CZ ") && l.contains("if b")),
            "expected an inline 'CZ ... if b...' line, got:\n{}",
            kmx
        );
        assert!(
            !kmx.lines()
                .any(|l| l.trim() == "PUSH_CONDITION if b0" || l.trim() == "POP_CONDITION"),
            "expected push/pop to be collapsed away, got:\n{}",
            kmx
        );
        c.free_bit(bit);
        let _ = c.destroy_sim(vec![a, b]);
    }

    #[test]
    fn unconditional_ops_unchanged() {
        let mut c = Circuit::new();
        let a = c.alloc_qreg("a");
        let b = c.alloc_qreg("b");
        c.cz(&a, &b);
        let kmx = write_kmx_to_string(&c);
        assert!(
            kmx.lines()
                .any(|l| l.starts_with("CZ ") && !l.contains("if b")),
            "expected a bare 'CZ q.. q..' line, got:\n{}",
            kmx
        );
        let _ = c.destroy_sim(vec![a, b]);
    }

    #[test]
    fn neg_if_bit_collapses() {
        let mut c = Circuit::new();
        let bit = c.alloc_bit();
        c.neg_if_bit(bit);
        let kmx = write_kmx_to_string(&c);
        assert!(
            kmx.lines().any(|l| l.starts_with("NEG if b")),
            "expected 'NEG if b..' line, got:\n{}",
            kmx
        );
        assert!(
            !kmx.lines().any(|l| l.trim() == "POP_CONDITION"),
            "POP_CONDITION should be collapsed:\n{}",
            kmx
        );
        c.free_bit(bit);
    }

    #[test]
    fn multi_op_pushcondition_block_not_collapsed() {
        // PUSH_CONDITION wrapping multiple distinct ops -> stays as 3+ instructions.
        let mut c = Circuit::new();
        let a = c.alloc_qreg("a");
        let b = c.alloc_qreg("b");
        let cc = c.alloc_qreg("cc");
        let bit = c.alloc_bit();
        c.with_condition(bit, |cir| {
            cir.cz(&a, &b);
            cir.cz(&a, &cc);
        });
        let kmx = write_kmx_to_string(&c);
        assert!(
            kmx.lines().any(|l| l.trim().starts_with("PUSH_CONDITION")),
            "multi-op block should preserve PUSH_CONDITION line:\n{}",
            kmx
        );
        assert!(
            kmx.lines().any(|l| l.trim() == "POP_CONDITION"),
            "multi-op block should preserve POP_CONDITION line:\n{}",
            kmx
        );
        c.free_bit(bit);
        let _ = c.destroy_sim(vec![a, b, cc]);
    }
}

#[cfg(test)]
mod declare_xor_of_tests {
    use super::*;

    /// End-to-end test for `declare_xor_of` via Circuit. Loads random
    /// values into q1, q2, computes target = q1 XOR q2 via two CXes,
    /// declares the identity (sim-verified across 64 shots), then HMRs
    /// target and discharges via two z_if_bit calls. `assert_phase_clean`
    /// must succeed.
    #[test]
    fn circuit_declare_xor_of_discharge_64lane() {
        let mut c = Circuit::new();
        // Use alloc_input_qreg so q1 / q2 are tracked as Top from the
        // start (sim_load_reg_bytes_shot can then place non-trivial
        // values; on_z later deposits CopyOf(q, ver) into the
        // obligation rather than short-circuiting on Zero).
        let q1 = c.alloc_input_qreg("q1");
        let q2 = c.alloc_input_qreg("q2");
        let t = c.alloc_qreg("t");
        // Random per-shot loads.
        use rand::Rng;
        let mut rng = rand::rngs::OsRng;
        for shot in 0..64 {
            if rng.gen::<bool>() {
                c.sim_load_reg_bytes_shot(std::slice::from_ref(&q1), &[1u8], shot);
            }
            if rng.gen::<bool>() {
                c.sim_load_reg_bytes_shot(std::slice::from_ref(&q2), &[1u8], shot);
            }
        }
        // t := q1 XOR q2. Two CXes; on_cx for the second flips the
        // tracker val into XorOf(...).
        c.cx(&q1, &t);
        c.cx(&q2, &t);
        // Sim-verify + tracker-inject. This OVERWRITES val(t) with the
        // XorOf(CopyOf(q1, v1), CopyOf(q2, v2)) flavor that we want for
        // the HMR obligation, regardless of whatever the CX chain
        // computed.
        c.declare_xor_of(&t, &q1, &q2);
        let bit = c.alloc_bit();
        c.hmr(&t, bit);
        c.z_if_bit(&q1, bit);
        c.z_if_bit(&q2, bit);
        c.free_bit(bit);
        c.assert_phase_clean();
        let _ = c.destroy_sim(vec![q1, q2, t]);
    }

    /// Sim-verification mismatch path: feeds a wrong XOR claim and
    /// expects a panic from the raw wrapper.
    #[test]
    #[should_panic(expected = "declare_xor_of")]
    fn circuit_declare_xor_of_sim_mismatch_panics() {
        let mut c = Circuit::new();
        let q1 = c.alloc_qreg("q1");
        let q2 = c.alloc_qreg("q2");
        let t = c.alloc_qreg("t");
        // Force shot 0 to t=1 while q1=q2=0 -> q1 XOR q2 = 0 != t.
        c.sim_load_reg_bytes_shot(std::slice::from_ref(&t), &[1u8], 0);
        c.declare_xor_of(&t, &q1, &q2);
    }

    /// End-to-end test for `declare_xor_of_three` via Circuit.
    #[test]
    fn circuit_declare_xor_of_three_discharge_64lane() {
        let mut c = Circuit::new();
        let q1 = c.alloc_input_qreg("q1");
        let q2 = c.alloc_input_qreg("q2");
        let q3 = c.alloc_input_qreg("q3");
        let t = c.alloc_qreg("t");
        use rand::Rng;
        let mut rng = rand::rngs::OsRng;
        for shot in 0..64 {
            if rng.gen::<bool>() {
                c.sim_load_reg_bytes_shot(std::slice::from_ref(&q1), &[1u8], shot);
            }
            if rng.gen::<bool>() {
                c.sim_load_reg_bytes_shot(std::slice::from_ref(&q2), &[1u8], shot);
            }
            if rng.gen::<bool>() {
                c.sim_load_reg_bytes_shot(std::slice::from_ref(&q3), &[1u8], shot);
            }
        }
        c.cx(&q1, &t);
        c.cx(&q2, &t);
        c.cx(&q3, &t);
        c.declare_xor_of_three(&t, &q1, &q2, &q3);
        let bit = c.alloc_bit();
        c.hmr(&t, bit);
        c.z_if_bit(&q1, bit);
        c.z_if_bit(&q2, bit);
        c.z_if_bit(&q3, bit);
        c.free_bit(bit);
        c.assert_phase_clean();
        let _ = c.destroy_sim(vec![q1, q2, q3, t]);
    }
}

#[cfg(test)]
mod cbit_naming_tests {
    use super::*;

    #[test]
    fn alloc_bit_anonymous_logs_default_tag() {
        // alloc_bit without a name still records a `b{N}` tag so the
        // debugger can resolve raw `bN` references via the alloc log.
        let mut c = Circuit::new();
        let b0 = c.alloc_bit();
        let b1 = c.alloc_bit();
        c.bit_store1(b0);
        c.bit_store0(b1);
        let op_idx = c.ops.len();
        assert_eq!(
            c.cbit_name_at_id(b0.raw(), op_idx),
            format!("b{}", b0.raw())
        );
        assert_eq!(
            c.cbit_name_at_id(b1.raw(), op_idx),
            format!("b{}", b1.raw())
        );
        c.free_bit(b0);
        c.free_bit(b1);
    }

    #[test]
    fn alloc_bit_named_records_tag() {
        let mut c = Circuit::new();
        let flag = c.alloc_bit_named("ready_flag");
        c.bit_store1(flag);
        let op_idx = c.ops.len();
        // Tag is scoped: with no active section/scope the tag is the
        // bare name; either way it must contain "ready_flag".
        let name = c.cbit_name_at_id(flag.raw(), op_idx);
        assert!(
            name == "ready_flag" || name.ends_with("/ready_flag"),
            "expected cbit name to be ready_flag (with optional scope), got {:?}",
            name,
        );
        c.free_bit(flag);
    }

    #[test]
    fn alloc_bits_named_indexes_each_bit() {
        let mut c = Circuit::new();
        let group = c.alloc_bits_named("ctrl", 4);
        assert_eq!(group.len(), 4);
        for (i, &b) in group.iter().enumerate() {
            if i % 2 == 0 {
                c.bit_store1(b);
            } else {
                c.bit_store0(b);
            }
        }
        let op_idx = c.ops.len();
        for (i, &b) in group.iter().enumerate() {
            let name = c.cbit_name_at_id(b.raw(), op_idx);
            let expected = format!("ctrl[{}]", i);
            assert!(
                name == expected || name.ends_with(&format!("/{}", expected)),
                "expected cbit[{}] name to be {:?}, got {:?}",
                i,
                expected,
                name,
            );
        }
        for b in group {
            c.free_bit(b);
        }
    }

    #[test]
    fn cbit_name_at_id_fallback_for_unknown() {
        let c = Circuit::new();
        // No allocs at all — the fallback is `b{N}`.
        assert_eq!(c.cbit_name_at_id(42, 0), "b42");
    }

    #[test]
    fn alloc_log_grows_with_named_and_anonymous() {
        let mut c = Circuit::new();
        let pre = c.bit_alloc_log.len();
        let anon = c.alloc_bit();
        let named = c.alloc_bit_named("trigger");
        let batch = c.alloc_bits_named("mask", 3);
        // 1 anon + 1 named + 3 batch = 5 new alloc log entries.
        assert_eq!(c.bit_alloc_log.len(), pre + 5);
        c.bit_store0(anon);
        c.bit_store0(named);
        for b in &batch {
            c.bit_store0(*b);
        }
        c.free_bit(anon);
        c.free_bit(named);
        for b in batch {
            c.free_bit(b);
        }
    }
}

#[cfg(test)]
mod redundancy_tests {
    use super::*;

    fn ccx_count_real(c: &Circuit) -> usize {
        c.ops
            .iter()
            .filter(|op| matches!(op, Some(Op::Ccx(_, _, _))))
            .count()
    }
    fn cx_count_real(c: &Circuit) -> usize {
        c.ops
            .iter()
            .filter(|op| matches!(op, Some(Op::Cx(_, _))))
            .count()
    }
    fn x_count_real(c: &Circuit) -> usize {
        c.ops
            .iter()
            .filter(|op| matches!(op, Some(Op::X(_))))
            .count()
    }

    #[test]
    fn x_after_x_elides() {
        let mut c = Circuit::new();
        let q = c.alloc_qreg("q");
        c.x(&q);
        c.x(&q);
        // Both X(q)'s elided to identity; sim should reflect q = |0>.
        assert_eq!(c.sim_get_mask(q.id), 0, "q should be |0> after X-X");
        assert_eq!(x_count_real(&c), 0, "no X ops should remain in self.ops");
        let _ = c.destroy_sim(vec![q]);
    }

    #[test]
    fn x_after_x_with_intervening_z_no_elide() {
        // Z and X are different ops, so X→Z→X doesn't trigger.
        let mut c = Circuit::new();
        let q = c.alloc_qreg("q");
        c.x(&q);
        c.z(&q);
        c.x(&q);
        // Both X's emit (Z between breaks the pair).
        assert_eq!(x_count_real(&c), 2);
        let _ = c.destroy_sim(vec![q]);
    }

    #[test]
    fn cx_after_cx_same_ctrl_elides() {
        let mut c = Circuit::new();
        let ctrl = c.alloc_qreg("ctrl");
        let tgt = c.alloc_qreg("tgt");
        c.cx(&ctrl, &tgt);
        c.cx(&ctrl, &tgt);
        assert_eq!(cx_count_real(&c), 0);
        let _ = c.destroy_sim(vec![ctrl, tgt]);
    }

    #[test]
    fn cx_after_cx_with_ctrl_modified_no_elide() {
        let mut c = Circuit::new();
        let ctrl = c.alloc_qreg("ctrl");
        let tgt = c.alloc_qreg("tgt");
        c.cx(&ctrl, &tgt);
        c.x(&ctrl); // modifies the control between the two CXs
        c.cx(&ctrl, &tgt);
        assert_eq!(cx_count_real(&c), 2);
        let _ = c.destroy_sim(vec![ctrl, tgt]);
    }

    #[test]
    fn cx_after_cx_with_target_modified_no_elide() {
        let mut c = Circuit::new();
        let ctrl = c.alloc_qreg("ctrl");
        let tgt = c.alloc_qreg("tgt");
        c.cx(&ctrl, &tgt);
        c.x(&tgt); // modifies the target between the two CXs
        c.cx(&ctrl, &tgt);
        assert_eq!(cx_count_real(&c), 2);
        let _ = c.destroy_sim(vec![ctrl, tgt]);
    }

    #[test]
    fn ccx_after_ccx_elides() {
        let mut c = Circuit::new();
        let a = c.alloc_qreg("a");
        let b = c.alloc_qreg("b");
        let t = c.alloc_qreg("t");
        c.ccx(&a, &b, &t);
        c.ccx(&a, &b, &t);
        assert_eq!(ccx_count_real(&c), 0);
        assert_eq!(c.ccx_emitted, 0);
        let _ = c.destroy_sim(vec![a, b, t]);
    }

    #[test]
    fn ccx_after_ccx_swapped_controls_elides() {
        // CCX is symmetric in its two controls; controls-swapped is
        // the same physical gate.
        let mut c = Circuit::new();
        let a = c.alloc_qreg("a");
        let b = c.alloc_qreg("b");
        let t = c.alloc_qreg("t");
        c.ccx(&a, &b, &t);
        c.ccx(&b, &a, &t);
        assert_eq!(ccx_count_real(&c), 0);
        let _ = c.destroy_sim(vec![a, b, t]);
    }

    #[test]
    fn ccx_after_ccx_with_one_ctrl_modified_no_elide() {
        let mut c = Circuit::new();
        let a = c.alloc_qreg("a");
        let b = c.alloc_qreg("b");
        let t = c.alloc_qreg("t");
        c.ccx(&a, &b, &t);
        c.x(&b);
        c.ccx(&a, &b, &t);
        assert_eq!(ccx_count_real(&c), 2);
        let _ = c.destroy_sim(vec![a, b, t]);
    }

    #[test]
    fn x_cx_x_does_not_elide() {
        // X(a); CX(a,b); X(a) — the CX touches a, breaking the chain.
        // Mathematically equivalent to negative-polarity CX, NOT identity.
        let mut c = Circuit::new();
        let a = c.alloc_qreg("a");
        let b = c.alloc_qreg("b");
        c.x(&a);
        c.cx(&a, &b);
        c.x(&a);
        assert_eq!(x_count_real(&c), 2, "X-CX-X must not elide");
        let _ = c.destroy_sim(vec![a, b]);
    }
}
