//! (refactor) Mechanically extracted from mod.rs. No logic changes.
use super::*;

pub(crate) struct B {
    pub ops: Vec<Op>,
    pub next_qubit: u64,
    pub next_bit: u64,
    pub next_register: u64,
    pub free_qubits: Vec<u64>,
    pub active_qubits: u32,
    pub peak_qubits: u32,
    pub peak_ops_idx: usize,
    pub peak_phase: &'static str,
    pub phase: &'static str,
    pub peak_log: Vec<(u32, &'static str, usize)>,
    pub phase_local_peaks: std::collections::BTreeMap<&'static str, (u32, usize)>,
    // (ops_len_at_transition, new_phase)
    pub phase_transitions: Vec<(usize, &'static str)>,
    // ── H201 diagnostic: TRACE_PEAK_OWNERS metadata-only owner tracking.
    // Default-off; populated only when env var TRACE_PEAK_OWNERS is set.
    // Each live qubit is associated with the phase that was active when it
    // was allocated (or with the explicit owner stack label if any).
    // Snapshots are recorded at every alloc that is within
    // TRACE_PEAK_OWNER_DELTA of the running peak; the final TRACE_PEAK block
    // filters them against the final peak and prints aggregates.
    pub owner_enabled: bool,
    pub owner_stack: Vec<&'static str>,
    pub owner_at_alloc: std::collections::BTreeMap<u64, &'static str>,
    // (active_count, phase_at_snapshot, ops_idx, owner_counts_grouped)
    pub owner_snapshots: Vec<(u32, &'static str, usize, std::collections::BTreeMap<&'static str, u32>)>,
}

impl B {
    pub(crate) fn new() -> Self {
        let owner_enabled = std::env::var("TRACE_PEAK_OWNERS").is_ok();
        Self {
            ops: Vec::new(),
            next_qubit: 0,
            next_bit: 0,
            next_register: 0,
            free_qubits: Vec::new(),
            active_qubits: 0,
            peak_qubits: 0,
            peak_ops_idx: 0,
            peak_phase: "",
            phase: "init",
            peak_log: Vec::new(),
            phase_local_peaks: std::collections::BTreeMap::new(),
            phase_transitions: Vec::new(),
            owner_enabled,
            owner_stack: Vec::new(),
            owner_at_alloc: std::collections::BTreeMap::new(),
            owner_snapshots: Vec::new(),
        }
    }
    /// Diagnostic helper: pushes a label onto the owner stack so subsequent
    /// allocations are attributed to that label (instead of the current
    /// phase name). Pops on Drop equivalent via paired call. METADATA-ONLY:
    /// has no effect when TRACE_PEAK_OWNERS is unset.
    #[allow(dead_code)]
    pub(crate) fn push_owner(&mut self, label: &'static str) {
        if self.owner_enabled {
            self.owner_stack.push(label);
        }
    }
    #[allow(dead_code)]
    pub(crate) fn pop_owner(&mut self) {
        if self.owner_enabled {
            self.owner_stack.pop();
        }
    }
    /// Scoped owner label: runs `f` with `label` active on the owner stack.
    /// METADATA-ONLY; no effect on emitted ops or qubit lifetimes.
    #[allow(dead_code)]
    pub(crate) fn with_owner<F: FnOnce(&mut B)>(&mut self, label: &'static str, f: F) {
        self.push_owner(label);
        f(self);
        self.pop_owner();
    }
    pub(crate) fn set_phase(&mut self, p: &'static str) {
        self.phase = p;
        self.phase_transitions.push((self.ops.len(), p));
    }
    pub(crate) fn alloc_qubit(&mut self) -> QubitId {
        self.active_qubits += 1;
        if self.active_qubits > self.peak_qubits {
            self.peak_qubits = self.active_qubits;
            self.peak_ops_idx = self.ops.len();
            self.peak_phase = self.phase;
            if std::env::var("TRACE_EACH_PEAK").is_ok() {
                eprintln!(
                    "PEAK active={} next_idx={} phase='{}' ops_idx={}",
                    self.active_qubits,
                    self.next_qubit,
                    self.phase,
                    self.ops.len()
                );
            }
        }
        if std::env::var("TRACE_PEAK").is_ok() && self.active_qubits + 10 >= self.peak_qubits {
            self.peak_log
                .push((self.active_qubits, self.phase, self.ops.len()));
        }
        if let Ok(prefix) = std::env::var("TRACE_PHASE_LOCAL_PEAK") {
            if !prefix.is_empty() && self.phase.starts_with(prefix.as_str()) {
                let entry = self
                    .phase_local_peaks
                    .entry(self.phase)
                    .or_insert((self.active_qubits, self.ops.len()));
                if self.active_qubits > entry.0 {
                    *entry = (self.active_qubits, self.ops.len());
                }
            }
        }
        let q = if let Some(q) = self.free_qubits.pop() {
            QubitId(q)
        } else {
            let q = self.next_qubit;
            self.next_qubit += 1;
            QubitId(q)
        };
        if self.owner_enabled {
            // Record this qubit's owner: top of owner_stack if present,
            // otherwise the current phase. Pure metadata.
            let owner: &'static str = self
                .owner_stack
                .last()
                .copied()
                .unwrap_or(self.phase);
            self.owner_at_alloc.insert(q.0, owner);
            // Take a near-peak snapshot at this allocation. The final
            // peak is unknown yet; we filter at print time using
            // TRACE_PEAK_OWNER_DELTA. We over-capture cheaply here:
            // snapshot every alloc within 64 of the running peak so we
            // never miss the final-peak band.
            if self.active_qubits + 64 >= self.peak_qubits {
                let mut counts: std::collections::BTreeMap<&'static str, u32> =
                    std::collections::BTreeMap::new();
                for (_qid, owner) in self.owner_at_alloc.iter() {
                    *counts.entry(*owner).or_insert(0) += 1;
                }
                self.owner_snapshots
                    .push((self.active_qubits, self.phase, self.ops.len(), counts));
            }
        }
        q
    }
    pub(crate) fn alloc_qubits(&mut self, n: usize) -> Vec<QubitId> {
        (0..n).map(|_| self.alloc_qubit()).collect()
    }
    pub(crate) fn alloc_bit(&mut self) -> BitId {
        let b = self.next_bit;
        self.next_bit += 1;
        BitId(b)
    }
    pub(crate) fn alloc_bits(&mut self, n: usize) -> Vec<BitId> {
        (0..n).map(|_| self.alloc_bit()).collect()
    }
    pub(crate) fn free(&mut self, q: QubitId) {
        self.r(q);
        self.free_qubits.push(q.0);
        if self.active_qubits > 0 {
            self.active_qubits -= 1;
        }
        if self.owner_enabled {
            self.owner_at_alloc.remove(&q.0);
        }
    }
    pub(crate) fn free_vec(&mut self, qs: &[QubitId]) {
        for &q in qs {
            self.free(q);
        }
    }
    pub(crate) fn declare_qubit_register(&mut self, qs: &[QubitId]) {
        let r = RegisterId(self.next_register);
        self.next_register += 1;
        for &q in qs {
            let mut op = Op::empty();
            op.kind = OperationType::AppendToRegister;
            op.q_target = q;
            op.r_target = r;
            self.ops.push(op);
        }
        let mut op = Op::empty();
        op.kind = OperationType::Register;
        op.r_target = r;
        self.ops.push(op);
    }
    pub(crate) fn declare_bit_register(&mut self, bs: &[BitId]) {
        let r = RegisterId(self.next_register);
        self.next_register += 1;
        for &b in bs {
            let mut op = Op::empty();
            op.kind = OperationType::AppendToRegister;
            op.c_target = b;
            op.r_target = r;
            self.ops.push(op);
        }
        let mut op = Op::empty();
        op.kind = OperationType::Register;
        op.r_target = r;
        self.ops.push(op);
    }
    pub(crate) fn x(&mut self, q: QubitId) {
        let mut op = Op::empty();
        op.kind = OperationType::X;
        op.q_target = q;
        self.ops.push(op);
    }
    pub(crate) fn cx(&mut self, ctrl: QubitId, tgt: QubitId) {
        let mut op = Op::empty();
        op.kind = OperationType::CX;
        op.q_control1 = ctrl;
        op.q_target = tgt;
        self.ops.push(op);
    }
    pub(crate) fn ccx(&mut self, c1: QubitId, c2: QubitId, tgt: QubitId) {
        let mut op = Op::empty();
        op.kind = OperationType::CCX;
        op.q_control2 = c1;
        op.q_control1 = c2;
        op.q_target = tgt;
        self.ops.push(op);
    }
    pub(crate) fn swap(&mut self, a: QubitId, b: QubitId) {
        let mut op = Op::empty();
        op.kind = OperationType::Swap;
        op.q_control1 = a;
        op.q_target = b;
        self.ops.push(op);
    }
    pub(crate) fn r(&mut self, q: QubitId) {
        let mut op = Op::empty();
        op.kind = OperationType::R;
        op.q_target = q;
        self.ops.push(op);
    }
    pub(crate) fn x_if(&mut self, q: QubitId, cond: BitId) {
        let mut op = Op::empty();
        op.kind = OperationType::X;
        op.q_target = q;
        op.c_condition = cond;
        self.ops.push(op);
    }
    // ── Measurement / phase / classical bit ops ──
    pub(crate) fn hmr(&mut self, q: QubitId, c: BitId) {
        let mut op = Op::empty();
        op.kind = OperationType::Hmr;
        op.q_target = q;
        op.c_target = c;
        self.ops.push(op);
    }
    // ── Classically-conditioned variants for all remaining gates ──
    pub(crate) fn cz_if(&mut self, a: QubitId, b: QubitId, cond: BitId) {
        let mut op = Op::empty();
        op.kind = OperationType::CZ;
        op.q_control1 = a;
        op.q_target = b;
        op.c_condition = cond;
        self.ops.push(op);
    }
    // Single-qubit classically-conditioned Z (sim: phase ^= cond & qubit(q)).
    // Equivalent to cz_if(c, q, cond) when c is a constant |1> ancilla.
    pub(crate) fn z_if(&mut self, q: QubitId, cond: BitId) {
        let mut op = Op::empty();
        op.kind = OperationType::Z;
        op.q_target = q;
        op.c_condition = cond;
        self.ops.push(op);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  emit_inverse: run a closure, pop the ops it emitted, and re-emit them
//  reversed.
//
//  The closure may contain `alloc_qubit` / `free` calls;
//  the R ops that `free` produces are SKIPPED during
//  reverse replay. This relies on the forward being "clean" — i.e. each
//  free lands on a qubit that the forward gates already drove to |0⟩
//  before the R. Under that invariant, the reverse gate sequence brings
//  the same qubit back to |0⟩ at the "alloc" point (pre-forward-allocation),
//  and the R we skipped is unnecessary.
//
//  The forward's internal alloc/free bookkeeping in the B's free
//  pool is NOT undone by the reverse — the pool state at reverse exit
//  equals the pool state at forward exit. Subsequent allocations in the
//  parent scope reuse those qubit IDs, seeing them at |0⟩ (as zeroed by
//  the reverse gate sequence).
// ═══════════════════════════════════════════════════════════════════════════
pub(crate) fn emit_inverse<F: FnOnce(&mut B)>(b: &mut B, f: F) {
    let start = b.ops.len();
    f(b);
    let end = b.ops.len();
    // Extract the forward slice and drop it from the builder.
    let fwd: Vec<_> = b.ops[start..end].to_vec();
    b.ops.truncate(start);
    for op in fwd.into_iter().rev() {
        match op.kind {
            OperationType::X
            | OperationType::Z
            | OperationType::CX
            | OperationType::CZ
            | OperationType::CCX
            | OperationType::CCZ
            | OperationType::Swap => b.ops.push(op),
            // R ops are the free markers. They're not directly reversible
            // as gates, but in a clean forward they're preceded by gates
            // that already zero the qubit. We skip them in reverse.
            OperationType::R => {}
            // Metadata ops (register declarations, debug prints) don't
            // affect state and shouldn't appear inside an emit_inverse
            // closure anyway, but skip them if they do.
            OperationType::Register
            | OperationType::AppendToRegister
            | OperationType::DebugPrint => {}
            _ => panic!(
                "emit_inverse: non-invertible op kind {:?} inside forward block",
                op.kind
            ),
        }
    }
}

pub const N: usize = 256;

/// secp256k1 prime:  p = 2^256 - 2^32 - 977.
pub const SECP256K1_P: U256 = U256::from_limbs([
    0xFFFFFFFEFFFFFC2F,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
]);
// ─── helpers: bit access on U256 ────────────────────────────────────────────

pub(crate) fn bit(c: U256, i: usize) -> bool {
    // alloy's U256::bit returns bool for index < 256.
    c.bit(i)
}

pub(crate) fn env_flag_enabled(name: &str, default: bool) -> bool {
    std::env::var(name)
        .ok()
        .map(|v| v != "0" && v.to_ascii_lowercase() != "false")
        .unwrap_or(default)
}

pub(crate) fn point_add_karatsuba_enabled() -> bool {
    env_flag_enabled("POINT_ADD_KARATSUBA", true)
}

/// Master switch for the "9n-floor" peak-drop construction (default OFF; when
/// OFF the emitted circuit is byte-identical to the f1d99ad C1 baseline).
///
/// The f1d99ad peak (2457) is a JOINT PIN across the pair2 Kaliski STEP-4
/// phases (kal_bulk_step4 / kal_step4 / bk_step4 / bk_bulk_step4 + their
/// step6_7_8 mod_double) AND the two binding multiplies (pair2_mul,
/// pair1_borrow_dx_mul1). Each STEP-4 holds a 256-wide `tmp` accumulator
/// co-resident with an ~255-wide measurement-Cuccaro carry register
/// (`sub_nbit_qq_fast` + `add_nbit_qq_fast`); each step6_7_8 holds a
/// `cadd_nbit_const_fast` const-register + carry register; each binding mul
/// holds the schoolbook Solinas-fold transient. Below the 2455-2457 cluster
/// the next floor is `shift22_step4` = 2333.
///
/// To move the GLOBAL peak below 2457 ALL THREE must drop simultaneously
/// (drop one alone -> another rebinds at 2457). When this flag is ON:
///   - STEP-4 sub/add: fast measurement-Cuccaro (n-1 carry register) ->
///     slow in-place Cuccaro (`add_nbit_qq`/`sub_nbit_qq`, 1 ancilla,
///     maj/uma carry rides in-place on tmp). Removes the ~255-carry register
///     from the binder at +~n CCX/op. Phase-clean + emit_inverse-safe (pure
///     CCX/CX, no measurement-vent).
///   - STEP 6/7/8 r-double: `mod_double_inplace_fast` (const reg + carries)
///     -> `mod_double_inplace_direct` (register-free direct const-add).
///   - pair2_mul + pair1_borrow_dx_mul1: schoolbook ->
///     `mod_mul_add_into_acc_schoolbook_lowscratch_fold` (the proven affine
///     y-mul peak-saver: lsx mul + lowscratch Solinas folds + direct doubles).
pub(crate) fn kal_gouzien_9n_enabled() -> bool {
    env_flag_enabled("KAL_GOUZIEN_9N", true)
}

/// Sub-knob: STEP-4 carry-register elimination (slow in-place Cuccaro). On by
/// default when the master flag is on. Diagnostic isolation: set to 0 to keep
/// fast Cuccaro at STEP-4 while still applying the mul/double drops.
pub(crate) fn gz_step4_slow() -> bool {
    kal_gouzien_9n_enabled() && env_flag_enabled("KAL_GZ_STEP4_SLOW", true)
}

/// Sub-knob: STEP-6/7/8 r-double register-free direct const-add.
pub(crate) fn gz_double_direct() -> bool {
    kal_gouzien_9n_enabled() && env_flag_enabled("KAL_GZ_DOUBLE_DIRECT", true)
}

/// Sub-knob: pair2_mul + pair1_borrow_dx_mul1 lowscratch Solinas fold.
pub(crate) fn gz_mul_lowscratch() -> bool {
    kal_gouzien_9n_enabled() && env_flag_enabled("KAL_GZ_MUL_LOWSCRATCH", true)
}

/// Sub-knob: late-iter Toffoli RECOVERY — widen the clean carry-borrow pool of
/// the STEP-4 q-q s-add/s-sub (the slow-Cuccaro fallback victim) by also
/// borrowing the GCD register `u`'s PROVABLY-|0> high bits.
///
/// Background: the 9n-floor (KAL_GOUZIEN_9N) hosts the fast-Cuccaro carry
/// register on clean future m_hist bits (`m_future`). The pool shrinks as the
/// walk advances, so once `(add_width-1) - |m_future|` exceeds KAL_GZ_MAX_FRESH
/// the s-add/s-sub falls back to a register-FREE in-place Cuccaro (+~n CCX).
/// That fallback fires for ALL iters with iter_idx > iters-1-KAL_GZ_MAX_FRESH
/// (~iter 277..iters) because s is full-width there while m_future is tiny —
/// it costs +~133k Toffoli (the +3.72% the 9n-floor paid for the 2333 peak).
///
/// Recovery: at the SAME late iters, the GCD walk has SHRUNK. The Kaliski
/// invariant guarantees `bitlen(u)+bitlen(v_w) <= 2n-iter_idx`, hence
/// individually `bitlen(u) <= 2n-iter_idx`, so for iter_idx>n the high bits
/// `u[2n-iter_idx..n)` (= iter_idx-n qubits) are PROVABLY |0>. They are already
/// allocated (part of the live `u` register) so borrowing them adds ZERO width.
/// They are not read between borrow and restore (the slow op is between `tmp`
/// and `s`; `u`'s value sits in [0..2n-iter_idx), untouched by the s-op and by
/// its own transform which runs strictly after (fwd) / before (bwd) the s-op).
/// The borrow is restored to |0> by the same measurement-uncompute as the
/// validated `cuccaro_*_fast_borrow`. The clean-bit boundary is a CLASSICAL
/// function of iter_idx (no data-dependent branch), so validity is structural.
///
/// Pool size becomes `|m_future| + (iter_idx-n) = iters-1-n` (~144-146) which is
/// constant across the late tail and keeps the fresh shortfall (255-pool ~110)
/// strictly below KAL_GZ_MAX_FRESH everywhere => NO slow Cuccaro at all, and
/// fewer fresh carries than before => per-step peak can only DROP, never rise.
pub(crate) fn gz_late_recover() -> bool {
    kal_gouzien_9n_enabled() && env_flag_enabled("KAL_GZ_LATE_RECOVER", true)
}

/// Sub-knob: shift22 (rshift22) Solinas reduction LOW-SCRATCH path. Replaces
/// the ~256-wide loaded-constant register of the shift22 STEP-3 const-add and
/// STEP-4 conditional const-sub (the 2333 binder) with Gidney venting
/// dirty-borrow const adders, and replaces the STEP-2 cuccaro spill-add's
/// ~257-wide `padded` register (the 2330 cuccaro_op_0 floor) with a narrow
/// k-bit spill add plus a venting dirty-borrow controlled-increment. Both cut
/// ~256 clean transient qubits at the affine y-mul Solinas-fold binder.
/// Gate KAL_GZ_SOLINAS_LOWSCRATCH (default on; set =0 to restore the byte-identical 1df6866 shift22).
pub(crate) fn gz_solinas_lowscratch() -> bool {
    kal_gouzien_9n_enabled() && env_flag_enabled("KAL_GZ_SOLINAS_LOWSCRATCH", true)
}

/// shift22 spill-fold COLLAPSE (default on): collapses the 5-op dirty shift22
/// spill fold into 2 ops by precomputing m=spill*977 and folding at pos 0 and
/// 32. Default ON; set SHIFT22_COLLAPSE=0 to restore the 5-op fold.
pub(crate) fn shift22_collapse() -> bool {
    env_flag_enabled("SHIFT22_COLLAPSE", true)
}

/// SOL_EXT_PRODUCT_POS32_FAST (default ON): use mod_add/sub_qq_fast for the
/// position-32 fold in mod_add/sub_solinas_ext_product. The Kaliski state is
/// not co-resident at the affine fold instant so the +512 transient is safe.
/// Set SOL_EXT_PRODUCT_POS32_FAST=0 to restore the slow fold.
pub(crate) fn sol_ext_product_pos32_fast() -> bool {
    gz_solinas_lowscratch() && env_flag_enabled("SOL_EXT_PRODUCT_POS32_FAST", true)
}

/// Provably-|0> high bits of the GCD register `u` at the late-iter STEP-4
/// instant: `u[2n-iter_idx .. n)` for iter_idx>n (else empty). Safe clean carry
/// donors for the s-add/s-sub borrow pool. See [`gz_late_recover`].
pub(crate) fn gz_u_clean_high(u: &[QubitId], iter_idx: usize) -> &[QubitId] {
    let n = u.len();
    if iter_idx <= n {
        return &u[..0];
    }
    let lo = (2 * n).saturating_sub(iter_idx).min(n);
    &u[lo..]
}

/// Sub-knob (lever C*): EARLY-iter carry-pool RELOCATION. The wide STEP-4
/// `v_w -= u` SUB (forward) / `v_w += u` ADD (backward) hosts its fast-Cuccaro
/// carry register on the GCD coefficient register `s`'s PROVABLY-|0> high bits
/// instead of the `m_future` (m_hist tail) pool, exactly mirroring
/// `gz_late_recover` (which hosts the LATE s-add/s-sub on u's |0> high bits).
///
/// Why this is the keystone: the dialog fold (`KAL_DIALOG_FOLD`) cannot relocate
/// the m_hist slots that the early SUB borrows as carries (they must stay clean
/// & dedicated). Moving those carries onto s-high frees ~n m_hist slots to fold,
/// dropping the inversion sweep peak at ~0 Toffoli (a pure carry-host relabel).
///
/// The s-high bits are clean by the SAME Kaliski invariant the late-recover and
/// the W-TRUNC layer already rely on: at iter i, the (r,s) step-3 cswap is width
/// `i+1` and the s +/-= tmp is width `i+2`, so `bitlen(s) <= i+1` at the v_w
/// SUB/ADD instant (forward: before the s-add; backward: after the s-sub undoes
/// it) ⇒ `s[i+1 .. n)` is |0> exactly when borrowed, and is restored to |0> by
/// the same measurement-uncompute as the validated `cuccaro_*_fast_borrow`.
/// The boundary is a CLASSICAL function of iter_idx (no data-dependent branch),
/// so validity is structural. BAKED DEFAULT ON (the validated C* construction);
/// `KAL_GZ_EARLY_RECOVER=0` restores the m_future-pool carry hosting.
pub(crate) fn gz_early_recover() -> bool {
    kal_gouzien_9n_enabled() && env_flag_enabled("KAL_GZ_EARLY_RECOVER", true)
}

/// Lowest index of `s` that is provably |0> at the STEP-4 v_w SUB/ADD instant of
/// iter `iter_idx`: `min(iter_idx+1, n)`. Single source of truth shared by the
/// hosting sites ([`gz_s_clean_high`]) and the dialog-fold slot map's carry
/// anchor (`dialog_fold_vw_slotmap`) so the two never disagree about how much of
/// the early SUB carry pool has been moved off `m_future`.
#[inline]
pub(crate) fn gz_s_clean_lo(iter_idx: usize, n: usize) -> usize {
    (iter_idx + 1).min(n)
}

/// Provably-|0> high bits of the coefficient register `s` at the STEP-4 v_w
/// SUB/ADD instant: `s[gz_s_clean_lo(iter_idx,n) .. n)`. Safe clean carry donors
/// for the wide early SUB/ADD borrow pool. See [`gz_early_recover`].
pub(crate) fn gz_s_clean_high(s: &[QubitId], iter_idx: usize) -> &[QubitId] {
    let n = s.len();
    let lo = gz_s_clean_lo(iter_idx, n);
    &s[lo..]
}

/// Provably-|0> high bits of the coefficient register `r` at the STEP-4 instant.
/// `bitlen(r) <= iter_idx+1` (symmetric to `s`: the step-3 cswap is width i+1 and
/// the s-add only writes r-via-transform over the low bits), so `r[i+1 .. n)` is
/// |0>. The SUB carry-borrow restores them to |0> before the transform reads r's
/// low bits and before the s-add reuses r-high. See [`gz_early_recover`].
pub(crate) fn gz_r_clean_high(r: &[QubitId], iter_idx: usize) -> &[QubitId] {
    let n = r.len();
    let lo = gz_s_clean_lo(iter_idx, n);
    &r[lo..]
}

/// STEP-4 load/SUB width: `(i<n ? n : 2n-i).min(wtrunc(i))`. Single source of
/// truth shared by the walk (kaliski_walk.rs) and the C* clean-pool sizing.
#[inline]
pub(crate) fn gz_load_width(iter_idx: usize, n: usize) -> usize {
    (if iter_idx < n { n } else { 2 * n - iter_idx }).min(kal_wtrunc_width(iter_idx, n))
}

/// Provably/empirically-|0> HIGH bits of the GCD register `u` available as clean
/// carry donors at the STEP-4 instant: `u[gz_load_width(i) .. n)`. The STEP-4
/// LOAD only reads `u[0..gz_load_width(i))` (it assumes the rest is 0 — the SAME
/// W-TRUNC envelope the circuit already validates 0/0/0 against), so `u`'s high
/// bits are |0> when borrowed, and the transform's cx/ccx/cx restores `u` exactly
/// so they are still |0> at the later s-add. Much wider than the provable
/// `gz_u_clean_high` for i<n (where that is empty), which is what lets the s-add
/// mid-iter borrow avoid pinning m_hist slots. Couples to the SAME island as the
/// existing load truncation (no new correctness risk beyond W-TRUNC).
pub(crate) fn gz_u_clean_high_wtrunc(u: &[QubitId], iter_idx: usize) -> &[QubitId] {
    let n = u.len();
    let lo = gz_load_width(iter_idx, n).min(n);
    &u[lo..]
}

/// Combined clean |0> carry-donor pool for the wide STEP-4 `v_w` SUB (fwd) / ADD
/// (bwd) under lever C*: s-high ++ r-high ++ u-high (W-TRUNC). `v_w`/`tmp` are the
/// operands (excluded); s, r are coefficients (high bits |0>) and u's high bits
/// are |0> by the W-TRUNC load envelope — none alias the SUB/ADD, and each is
/// restored to |0> by the same measurement-uncompute as the validated
/// late-recover. Drawn BEFORE `m_future`, so the m_hist tail keeps almost none of
/// the carry pool and folds away. Identical in fwd & bwd (same qubit ids) ⇒
/// reversible.
pub(crate) fn gz_vw_clean_pool(
    s: &[QubitId],
    r: &[QubitId],
    u: &[QubitId],
    iter_idx: usize,
) -> Vec<QubitId> {
    let mut v: Vec<QubitId> = Vec::new();
    v.extend_from_slice(gz_s_clean_high(s, iter_idx));
    v.extend_from_slice(gz_r_clean_high(r, iter_idx));
    v.extend_from_slice(gz_u_clean_high_wtrunc(u, iter_idx));
    v
}

/// Combined clean |0> carry-donor pool for the STEP-4 `s` add (fwd) / sub (bwd)
/// under lever C*: r-high ++ u-high (W-TRUNC) (`s` is the accumulator so s-high is
/// excluded). r-high is big early, u-high is big mid/late (anti-correlated), so
/// the pool covers the s-add need with ~0 `m_future` shortfall except in a narrow
/// √n boundary band near i≈n. See [`gz_early_recover`]; generalizes
/// [`gz_late_recover`]'s u-high hosting to the W-TRUNC envelope.
pub(crate) fn gz_s_clean_pool(r: &[QubitId], u: &[QubitId], iter_idx: usize) -> Vec<QubitId> {
    let mut v: Vec<QubitId> = Vec::new();
    v.extend_from_slice(gz_r_clean_high(r, iter_idx));
    v.extend_from_slice(gz_u_clean_high_wtrunc(u, iter_idx));
    v
}

/// Master switch for the stacked sub-2708-peak construction (default ON).
/// When ON (no env / != "0"):
///   - pair1 inverse borrows dx as in-place v_w  (drops pair1-backward carrier)
///   - pair2 inverse borrows tx as in-place v_w  (drops pair2-forward carrier)
///   - pair1_mul1 + pair2_mul use schoolbook (no Karatsuba z1_reg) so the mul
///     Solinas boundary frees ~258q  -> drops sol_sub6 / kara_z1_add below 2565
///   - Kaliski iters bumped (405/403) to clear the in-place-v + schoolbook +
///     direct-const-halve correctness/phase margin (9024-shot validated).
/// Net: global peak 2708 -> 2565, score 9.604e9 -> 9.361e9 (-11.34% vs baseline).
/// Set POINT_ADD_STACK_2565=0 to restore the byte-identical C1 circuit.
pub(crate) fn stack_2565_enabled() -> bool {
    env_flag_enabled("POINT_ADD_STACK_2565", true)
}

pub(crate) fn pair1_mul1_karatsuba_enabled(n: usize) -> bool {
    let min_n = std::env::var("POINT_ADD_KARATSUBA_MIN_N")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(256);
    // Stack default: pair1_mul1 is schoolbook (frees the Karatsuba z1_reg from
    // the binding sol_sub6 boundary). User can force Karatsuba via the knob.
    let karatsuba_default = !stack_2565_enabled();
    point_add_karatsuba_enabled()
        && n >= min_n
        && env_flag_enabled("KAL_PAIR1_MUL1_KARATSUBA", karatsuba_default)
}

pub(crate) fn direct_const_halve_enabled() -> bool {
    // The direct constant subtract halve is very slightly lower-peak by itself,
    // but older guarded Karatsuba attempts found that combining it with
    // pair1_mul1 Karatsuba can hit a phase-cleanliness cliff on alternate
    // seeds.  Prefer the revived Karatsuba win by default; both knobs remain
    // independently overrideable for diagnostics.
    env_flag_enabled("KAL_DIRECT_CONST_HALVE", !pair1_mul1_karatsuba_enabled(N))
}

pub(crate) fn pair1_mul2_karatsuba_enabled(n: usize) -> bool {
    let min_n = std::env::var("POINT_ADD_KARATSUBA_MIN_N")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(256);
    point_add_karatsuba_enabled()
        && n >= min_n
        && env_flag_enabled("KAL_PAIR1_MUL2_KARATSUBA", true)
}

pub(crate) fn pair2_mul_karatsuba_enabled(n: usize) -> bool {
    let min_n = std::env::var("POINT_ADD_KARATSUBA_MIN_N")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(256);
    // Stack default: pair2_mul is schoolbook (frees the Karatsuba z1_reg from
    // the binding kara_z1_add boundary). User can force Karatsuba via the knob.
    let karatsuba_default = !stack_2565_enabled();
    point_add_karatsuba_enabled()
        && n >= min_n
        && env_flag_enabled("KAL_PAIR2_MUL_KARATSUBA", karatsuba_default)
}
