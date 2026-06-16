//! Sound abstract-interpretation phase lattice with O(1) per gate.
//!
//! Tracks phase obligations from HMR via version-tagged qubit
//! references. Each qubit has a monotonic version counter
//! (incremented on every write). An MBUC obligation records the
//! CCX parents' (`qubit_id`, version) at HMR time; the matching
//! `cz_if_bit` checks those same (`qubit_id`, version) pairs are
//! still current.
//!
//! Lattice elements per qubit:
//!   Zero                        known |0⟩
//!   One                         known |1⟩
//!   CopyOf(q, v)                equals qubit q at version v
//!   AndOf(q1,v1,q2,v2)          equals q1@v1 AND q2@v2
//!   XorOf(a, b)                 depth-1 XOR of two sub-values
//!   Anchor(id)                  opaque equality class (`declare_identity`)
//!   Top                         unknown
//!
//! The tracker ALSO catches R-on-nonzero (zenodo R semantics
//! kicks `qubit & rng` into phase, which the tracker can't cancel).
//! Every R on a qubit whose abstract value isn't provably Zero is
//! reported at `assert_clean`.
//!
//! Trust assertions (`declare_identity` / `declare_copy_of` /
//! `declare_and_of`) let a caller PROVE facts the tracker cannot
//! derive structurally. Each call site carries a proof comment.
//! These are SYMBOLIC — they don't check at runtime. For a
//! stronger concrete-backed version see `Circuit::prove_zero`.

use std::collections::HashMap;

pub type QubitId = u32;
pub type Version = u64;
pub type Atom = u32;

#[derive(Clone, Copy, Debug, Default)]
pub struct PhaseTrackerStats {
    pub val_len: usize,
    pub obligations: usize,
    pub obligation_discharges: usize,
    pub r_on_nonzero: usize,
    pub hmr_atom_to_bit: usize,
    pub next_atom: Atom,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AbsVal {
    Zero,
    One,
    CopyOf(QubitId, Version),
    AndOf(QubitId, Version, QubitId, Version),
    /// 3-way AND of three qubits. Injected via `declare_and3_of`.
    /// Used for MBU of qubits computed as x AND y AND z (a common
    /// phase pattern for `ccz_if_bit` discharges).
    AndOf3(QubitId, Version, QubitId, Version, QubitId, Version),
    /// Depth-1 XOR of two sub-values (no nesting).
    XorOf(Box<AbsVal>, Box<AbsVal>),
    /// Depth-1 XOR of three sub-values (no nesting). Injected via
    /// `declare_xor_of_three`. The matcher in `assert_clean` accepts
    /// three unpaired discharges matching the three leaves in any
    /// order. Use when the HMR'd qubit's value is provably a 3-term
    /// XOR of values that can each be deposited via `z_if_bit` (or
    /// `cz_if_bit` / `ccz_if_bit` for AND-of leaves).
    XorOf3(Box<AbsVal>, Box<AbsVal>, Box<AbsVal>),
    /// Opaque equality class for `declare_identity`.
    Anchor(u64),
    /// Multiplexed value: `target = ctrl ? a : b` where `ctrl`, `a`,
    /// `b` are version-tagged qubit references. Injected via
    /// `declare_choose_of` -- used by the Luo-pack passthrough
    /// shift gadgets, where a controlled rotate gates a qubit's value
    /// on a quantum case-flag.
    ///
    /// Currently this variant carries the same "structural fact" role
    /// as CopyOf/AndOf: it tracks that the qubit is a known function
    /// of three live qubits, rather than Top. HMR discharge of a
    /// `ChooseOf` qubit is NOT supported -- callers must uncompute the
    /// choose-of structure (e.g. via the inverse passthrough shift)
    /// before HMR'ing this qubit.
    ChooseOf(
        /* ctrl */ QubitId,
        Version,
        /* a    */ QubitId,
        Version,
        /* b    */ QubitId,
        Version,
    ),
    Top,
}

impl AbsVal {
    fn is_zero(&self) -> bool {
        matches!(self, AbsVal::Zero)
    }
    fn is_top(&self) -> bool {
        matches!(self, AbsVal::Top)
    }
}

/// True iff the obligation's discharges structurally match its val.
/// Same logic as `assert_clean`'s check (pair-cancel + single/XorOf
/// match), but returns a bool instead of emitting an error.
fn obligation_is_clean(ob: &Obligation) -> bool {
    if ob.val.is_zero() {
        return true;
    }
    if ob.val.is_top() {
        return false;
    }
    let mut remaining: Vec<&(AbsVal, bool, usize)> = ob.discharges.iter().collect();
    remaining.sort_by(|a, b| format!("{:?}", a.0).cmp(&format!("{:?}", b.0)));
    let mut unpaired: Vec<&(AbsVal, bool, usize)> = Vec::new();
    let mut i = 0;
    while i < remaining.len() {
        if i + 1 < remaining.len()
            && remaining[i].0 == remaining[i + 1].0
            && remaining[i].1
            && remaining[i + 1].1
        {
            i += 2;
        } else {
            unpaired.push(remaining[i]);
            i += 1;
        }
    }
    if unpaired.len() == 1 && unpaired[0].0 == ob.val && unpaired[0].1 {
        return true;
    }
    if let AbsVal::XorOf(ref oa, ref ob_inner) = ob.val {
        if unpaired.len() == 2 {
            let (d0, d1) = (&unpaired[0].0, &unpaired[1].0);
            let (v0, v1) = (unpaired[0].1, unpaired[1].1);
            return (d0 == oa.as_ref() && d1 == ob_inner.as_ref() && v0 && v1)
                || (d0 == ob_inner.as_ref() && d1 == oa.as_ref() && v0 && v1);
        }
    }
    if let AbsVal::XorOf3(ref oa, ref ob_inner, ref oc) = ob.val {
        if unpaired.len() == 3 && unpaired.iter().all(|d| d.1) {
            let leaves = [oa.as_ref(), ob_inner.as_ref(), oc.as_ref()];
            let discharges = [&unpaired[0].0, &unpaired[1].0, &unpaired[2].0];
            return matches_as_multiset(&leaves, &discharges);
        }
    }
    false
}

/// True iff `discharges` and `leaves` (both of length 3) contain the
/// same multiset of `AbsVal`. Each leaf must be matched by exactly one
/// distinct discharge; used by the `XorOf3` obligation-match arm.
fn matches_as_multiset(leaves: &[&AbsVal; 3], discharges: &[&AbsVal; 3]) -> bool {
    // Brute-force over 3! permutations: with only 3 elements this is
    // 6 comparisons and stays simple/auditable.
    let perms = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    perms.iter().any(|p| {
        discharges[0] == leaves[p[0]]
            && discharges[1] == leaves[p[1]]
            && discharges[2] == leaves[p[2]]
    })
}

#[derive(Clone, Debug)]
struct Obligation {
    val: AbsVal,
    discharges: Vec<(AbsVal, bool, usize)>, // (val, versions_current, op_idx)
    section: String,
    hmr_q: u32,
    hmr_op_idx: usize,
    /// Multi-term ghost accounting (sim-mask based). `hmr_mask` is the
    /// 64-shot value of the vented qubit at HMR; `discharge_xor` is the
    /// running XOR of the discharge terms' 64-shot masks. When a
    /// `close_ghost`-style resolution fires, the tracker requires
    /// `discharge_xor == hmr_mask` (the discharges reproduce the vented
    /// value on every shot, so the HMR kickback cancels for any random
    /// bit) before removing the obligation. `None` for ordinary HMRs /
    /// single-term anchor resolves, which use structural matching.
    hmr_mask: Option<u64>,
    discharge_xor: u64,
}

pub struct PhaseTracker {
    val: Vec<AbsVal>,
    ver: Vec<Version>,
    cond_stack: Vec<u32>,
    obligations: HashMap<u32, Obligation>,
    direct_phase_nonzero: bool,
    /// R-on-nonzero samples (capped per section) and total counts.
    r_on_nonzero: Vec<(String, u32, AbsVal)>,
    r_on_nonzero_count: HashMap<String, u64>,
    next_anchor: u64,
    pub current_section: String,
    enabled: bool,
    pub next_atom: Atom,
    pub hmr_atom_to_bit: HashMap<Atom, u32>,
    /// Current op index, updated by `Circuit.push_gate_op`. Used for
    /// diagnostic op-idx tracking on obligations/discharges so the
    /// debugger REPL can `goto` the offending ops.
    pub current_op_idx: usize,
}

impl Default for PhaseTracker {
    fn default() -> Self {
        Self::new()
    }
}

impl PhaseTracker {
    #[must_use]
    pub fn new() -> Self {
        Self {
            val: Vec::new(),
            ver: Vec::new(),
            cond_stack: Vec::new(),
            obligations: HashMap::new(),
            direct_phase_nonzero: false,
            r_on_nonzero: Vec::new(),
            r_on_nonzero_count: HashMap::new(),
            next_anchor: 1,
            current_section: String::new(),
            enabled: true,
            next_atom: 1,
            hmr_atom_to_bit: HashMap::new(),
            current_op_idx: 0,
        }
    }

    #[must_use]
    pub fn interner_size(&self) -> usize {
        0
    }

    #[must_use]
    pub fn stats(&self) -> PhaseTrackerStats {
        PhaseTrackerStats {
            val_len: self.val.len(),
            obligations: self.obligations.len(),
            obligation_discharges: self
                .obligations
                .values()
                .map(|ob| ob.discharges.len())
                .sum(),
            r_on_nonzero: self.r_on_nonzero.len(),
            hmr_atom_to_bit: self.hmr_atom_to_bit.len(),
            next_atom: self.next_atom,
        }
    }

    fn ensure_qubit(&mut self, q: u32) {
        let q = q as usize;
        if self.val.len() <= q {
            self.val.resize(q + 1, AbsVal::Zero);
            self.ver.resize(q + 1, 0);
        }
    }

    pub fn fresh_atom(&mut self) -> Atom {
        let a = self.next_atom;
        self.next_atom += 1;
        a
    }

    pub fn mark_input_qubit(&mut self, q: u32) {
        if !self.enabled {
            return;
        }
        self.ensure_qubit(q);
        self.val[q as usize] = AbsVal::Top;
    }

    pub fn mark_input_bit(&mut self, _b: u32) {
        if !self.enabled {}
    }

    fn qval(&self, q: u32) -> &AbsVal {
        self.val.get(q as usize).unwrap_or(&AbsVal::Zero)
    }
    fn qver(&self, q: u32) -> Version {
        self.ver.get(q as usize).copied().unwrap_or(0)
    }

    fn write_qubit(&mut self, q: u32, v: AbsVal) {
        self.ensure_qubit(q);
        self.ver[q as usize] += 1;
        self.val[q as usize] = v;
    }

    // ── Fact injection (INTERNAL — use Circuit wrappers, not these)
    //
    // These `_unchecked` variants inject symbolic facts into the
    // tracker WITHOUT runtime verification. They are unsound on
    // their own. The ONLY supported callers are the Circuit
    // wrappers (`Circuit::declare_identity`, `declare_copy_of`,
    // `declare_and_of`), which precede each fact with a 64-shot
    // sim-time identity assertion so the injected symbolic fact
    // is always backed by concrete evidence on random inputs.
    // Direct calls to these `_unchecked` methods from outside
    // phase_lattice.rs / circuit.rs are a soundness bug.

    fn declare_identity_unchecked(&mut self, q_a: u32, q_b: u32) {
        if !self.enabled {
            return;
        }
        self.ensure_qubit(q_a);
        self.ensure_qubit(q_b);
        let id = self.next_anchor;
        self.next_anchor += 1;
        self.val[q_a as usize] = AbsVal::Anchor(id);
        self.val[q_b as usize] = AbsVal::Anchor(id);
    }

    fn declare_copy_of_unchecked(&mut self, q: u32, source: u32) {
        if !self.enabled {
            return;
        }
        self.ensure_qubit(q);
        self.ensure_qubit(source);
        self.val[q as usize] = self.copy_of_or_val(source);
    }

    fn declare_and_of_unchecked(&mut self, q: u32, a: u32, b: u32) {
        if !self.enabled {
            return;
        }
        self.ensure_qubit(q);
        self.ensure_qubit(a);
        self.ensure_qubit(b);
        let av = self.qval(a).clone();
        let bv = self.qval(b).clone();
        let new = match (av, bv) {
            (AbsVal::Zero, _) | (_, AbsVal::Zero) => AbsVal::Zero,
            (AbsVal::One, _) => self.copy_of_or_val(b),
            (_, AbsVal::One) => self.copy_of_or_val(a),
            _ => AbsVal::AndOf(a, self.qver(a), b, self.qver(b)),
        };
        self.val[q as usize] = new;
    }

    fn declare_and3_of_unchecked(&mut self, q: u32, a: u32, b: u32, c: u32) {
        if !self.enabled {
            return;
        }
        self.ensure_qubit(q);
        self.ensure_qubit(a);
        self.ensure_qubit(b);
        self.ensure_qubit(c);
        self.val[q as usize] = AbsVal::AndOf3(a, self.qver(a), b, self.qver(b), c, self.qver(c));
    }

    /// `target := ctrl ? choice_a : choice_b`. Records the multiplexed
    /// fact at the current versions of all three sources. See `AbsVal::`
    /// `ChooseOf` docs for usage constraints.
    /// `target := q1 XOR q2`. Records val(target) = XorOf(CopyOf(q1, v1),
    /// CopyOf(q2, v2)) where v1, v2 are the current versions of q1 and
    /// q2. This is the depth-1 XOR analogue of `declare_and_of` / `declare_identity`:
    /// it lets the tracker structurally match an HMR(target) obligation
    /// against two `z_if_bit(q1`, ...) + `z_if_bit(q2`, ...) discharges,
    /// via the `XorOf` branch of `assert_clean`.
    ///
    /// Constraints on the caller:
    ///  - Between this declare and the matching `z_if_bit` calls, neither
    ///    q1 nor q2 may be modified (any write bumps ver, making the
    ///    obligation's CopyOf(q, v) stale).
    ///  - target's val is overwritten — any prior tracker fact about
    ///    target is discarded.
    ///  - Sim verification (in the public wrapper `Circuit::declare_xor_of`)
    ///    asserts `sim_mask(target)` == `sim_mask(q1)` ^ `sim_mask(q2)` on all
    ///    64 shots before injection.
    ///
    /// Leaves are always recorded as `CopyOf(q, ver(q))`, never simplified
    /// to Zero/One/Anchor, so `on_z`'s fallthrough case (which deposits
    /// `CopyOf(q, ver(q))` for any non-Anchor/non-One val) lines up with
    /// the obligation's stored leaves. This means callers can ignore the
    /// current `AbsVal` of q1/q2 as long as their versions don't bump.
    fn declare_xor_of_unchecked(&mut self, target: u32, q1: u32, q2: u32) {
        if !self.enabled {
            return;
        }
        self.ensure_qubit(target);
        self.ensure_qubit(q1);
        self.ensure_qubit(q2);
        let leaf1 = AbsVal::CopyOf(q1, self.qver(q1));
        let leaf2 = AbsVal::CopyOf(q2, self.qver(q2));
        self.val[target as usize] = AbsVal::XorOf(Box::new(leaf1), Box::new(leaf2));
    }

    /// `target := q1 XOR q2 XOR q3` — three-term XOR analogue of
    /// `declare_xor_of`. Leaves are recorded as `CopyOf(q_i, ver(q_i))`
    /// to match `on_z`'s fallthrough discharge contribution; the
    /// XorOf3-branch of `assert_clean` accepts three unpaired
    /// discharges matching the leaves in any order.
    ///
    /// Same caller obligations as `declare_xor_of`: q1, q2, q3 must
    /// not be modified between this declare and the matching
    /// `z_if_bit` calls.
    fn declare_xor_of_three_unchecked(&mut self, target: u32, q1: u32, q2: u32, q3: u32) {
        if !self.enabled {
            return;
        }
        self.ensure_qubit(target);
        self.ensure_qubit(q1);
        self.ensure_qubit(q2);
        self.ensure_qubit(q3);
        let leaf1 = AbsVal::CopyOf(q1, self.qver(q1));
        let leaf2 = AbsVal::CopyOf(q2, self.qver(q2));
        let leaf3 = AbsVal::CopyOf(q3, self.qver(q3));
        self.val[target as usize] =
            AbsVal::XorOf3(Box::new(leaf1), Box::new(leaf2), Box::new(leaf3));
    }

    fn declare_choose_of_unchecked(&mut self, target: u32, ctrl: u32, a: u32, b: u32) {
        if !self.enabled {
            return;
        }
        self.ensure_qubit(target);
        self.ensure_qubit(ctrl);
        self.ensure_qubit(a);
        self.ensure_qubit(b);
        // Structural simplifications: if ctrl is known constant, fall
        // back to a simpler tracker fact.
        let ctrl_v = self.qval(ctrl).clone();
        let new = match ctrl_v {
            AbsVal::One => self.copy_of_or_val(a),
            AbsVal::Zero => self.copy_of_or_val(b),
            _ => AbsVal::ChooseOf(ctrl, self.qver(ctrl), a, self.qver(a), b, self.qver(b)),
        };
        self.val[target as usize] = new;
    }

    fn copy_of_or_val(&self, source: u32) -> AbsVal {
        match self.qval(source).clone() {
            AbsVal::Zero => AbsVal::Zero,
            AbsVal::One => AbsVal::One,
            AbsVal::Anchor(id) => AbsVal::Anchor(id),
            _ => AbsVal::CopyOf(source, self.qver(source)),
        }
    }

    /// Internal: inject `val(q) = Zero` into the tracker.
    /// Callers outside `phase_lattice.rs` / circuit.rs MUST go through
    /// `Circuit::prove_zero` (which panics unless sim is active and
    /// q reads |0> on all 64 shots). Direct use is a soundness hole.
    fn inject_zero_unchecked(&mut self, q: u32) {
        if !self.enabled {
            return;
        }
        self.ensure_qubit(q);
        self.val[q as usize] = AbsVal::Zero;
    }

    // ── Transfer functions ──────────────────────────────────────

    pub fn on_alloc_qubit(&mut self, q: u32) {
        if !self.enabled {
            return;
        }
        self.ensure_qubit(q);
        self.ver[q as usize] += 1;
        self.val[q as usize] = AbsVal::Zero;
    }
    pub fn on_alloc_bit(&mut self, _b: u32) {}
    pub fn on_free_qubit(&mut self, _q: u32) {}
    pub fn on_free_bit(&mut self, b: u32) {
        if !self.enabled {
            return;
        }
        // At free_bit time: if this bit has a CLEANLY MATCHED
        // obligation, remove it so the next alloc_bit returning this
        // ID doesn't accumulate spurious discharges against it.
        // If unmatched, leave it — assert_clean will report it.
        if let Some(ob) = self.obligations.get(&b) {
            if obligation_is_clean(ob) {
                self.obligations.remove(&b);
            }
        }
    }

    pub fn on_x(&mut self, q: u32) {
        if !self.enabled {
            return;
        }
        let old = self.qval(q).clone();
        let new = match old {
            AbsVal::Zero => AbsVal::One,
            AbsVal::One => AbsVal::Zero,
            _ => AbsVal::Top,
        };
        self.write_qubit(q, new);
    }

    pub fn on_cx(&mut self, c: u32, q: u32) {
        if !self.enabled {
            return;
        }
        self.ensure_qubit(c);
        self.ensure_qubit(q);
        let qv = self.qval(q).clone();
        let c_val = AbsVal::CopyOf(c, self.qver(c));
        let new = match qv {
            AbsVal::Zero => c_val,
            AbsVal::CopyOf(q2, v2) if q2 == c && v2 == self.qver(c) => AbsVal::Zero,
            AbsVal::CopyOf(..) | AbsVal::AndOf(..) => AbsVal::XorOf(Box::new(qv), Box::new(c_val)),
            AbsVal::XorOf(ref a, ref b) => {
                if **a == c_val {
                    (**b).clone()
                } else if **b == c_val {
                    (**a).clone()
                } else {
                    AbsVal::Top
                }
            }
            _ => AbsVal::Top,
        };
        self.write_qubit(q, new);
    }

    pub fn on_ccx(&mut self, a: u32, b: u32, q: u32) {
        if !self.enabled {
            return;
        }
        self.ensure_qubit(a);
        self.ensure_qubit(b);
        self.ensure_qubit(q);
        let qv = self.qval(q).clone();
        let and_val = AbsVal::AndOf(a, self.qver(a), b, self.qver(b));
        let new = match qv {
            AbsVal::Zero => and_val,
            _ if qv == and_val => AbsVal::Zero,
            AbsVal::CopyOf(..) | AbsVal::AndOf(..) => {
                AbsVal::XorOf(Box::new(qv), Box::new(and_val))
            }
            AbsVal::XorOf(ref va, ref vb) => {
                if **va == and_val {
                    (**vb).clone()
                } else if **vb == and_val {
                    (**va).clone()
                } else {
                    AbsVal::Top
                }
            }
            _ => AbsVal::Top,
        };
        self.write_qubit(q, new);
    }

    pub fn on_swap(&mut self, a: u32, b: u32) {
        if !self.enabled {
            return;
        }
        self.ensure_qubit(a);
        self.ensure_qubit(b);
        self.val.swap(a as usize, b as usize);
        self.ver.swap(a as usize, b as usize);
    }

    /// R: zenodo semantics `phase ^= qubit & rng; qubit = 0`. R on
    /// a tracker-Zero qubit is safe; R on anything else leaks a
    /// rng-dependent phase that can never cancel and is recorded
    /// here. Reported at `assert_clean`.
    pub fn on_r(&mut self, q: u32) {
        if !self.enabled {
            return;
        }
        self.ensure_qubit(q);
        let qv = &self.val[q as usize];
        if !qv.is_zero() {
            let sec = self.current_section.clone();
            let count = self.r_on_nonzero_count.entry(sec.clone()).or_insert(0);
            *count += 1;
            if *count <= 8 {
                self.r_on_nonzero.push((sec, q, qv.clone()));
            }
            if std::env::var("PL_TRACE").is_ok() && *count == 1 {
                eprintln!(
                    "[PL] on_r(q={}) val={:?} section='{}'",
                    q, qv, self.current_section
                );
            }
        }
        self.write_qubit(q, AbsVal::Zero);
    }

    // ── Phase transfers ─────────────────────────────────────────

    pub fn on_z(&mut self, q: u32) {
        if !self.enabled {
            return;
        }
        self.ensure_qubit(q);
        let qv = self.qval(q).clone();
        if qv.is_zero() {
            return;
        }
        let contrib = match qv {
            AbsVal::Anchor(_) => qv,
            AbsVal::One => AbsVal::One,
            _ => AbsVal::CopyOf(q, self.qver(q)),
        };
        self.discharge_to_cond(contrib);
    }

    pub fn on_cz(&mut self, a: u32, b: u32) {
        if !self.enabled {
            return;
        }
        let val = AbsVal::AndOf(a, self.qver(a), b, self.qver(b));
        self.discharge_to_cond(val);
    }

    pub fn on_ccz(&mut self, a: u32, b: u32, c: u32) {
        if !self.enabled {
            return;
        }
        let val = AbsVal::AndOf3(a, self.qver(a), b, self.qver(b), c, self.qver(c));
        self.discharge_to_cond(val);
    }

    pub fn on_neg(&mut self) {
        if !self.enabled {
            return;
        }
        self.discharge_to_cond(AbsVal::One);
    }

    fn discharge_to_cond(&mut self, val: AbsVal) {
        if val.is_zero() {
            return;
        }
        if self.cond_stack.is_empty() {
            self.direct_phase_nonzero = true;
        } else if self.cond_stack.len() == 1 {
            let bit = self.cond_stack[0];
            let valid = self.val_versions_current(&val);
            let idx = self.current_op_idx;
            if let Some(ob) = self.obligations.get_mut(&bit) {
                ob.discharges.push((val, valid, idx));
            }
        }
    }

    // ── HMR + condition stack ───────────────────────────────────

    pub fn on_hmr(&mut self, q: u32, b: u32) {
        if !self.enabled {
            return;
        }
        assert!(
            self.cond_stack.is_empty(),
            "HMR inside push_condition forbidden"
        );
        self.ensure_qubit(q);
        // Detect bit-ID reuse: if b already has an obligation, the
        // prior mbuc_free didn't clear it. This corrupts the next
        // HMR's accounting. Print loud diagnostic.
        if std::env::var("PL_TRACE").is_ok() {
            if let Some(old) = self.obligations.get(&b) {
                eprintln!(
                    "[PL] on_hmr(q={}, bit={}) at op={} \
                    OVERWRITES existing obligation from op={} \
                    (hmr_q={}, val={:?}, {} discharges)",
                    q,
                    b,
                    self.current_op_idx,
                    old.hmr_op_idx,
                    old.hmr_q,
                    old.val,
                    old.discharges.len()
                );
            }
        }
        let qv = self.qval(q).clone();
        self.obligations.insert(
            b,
            Obligation {
                val: qv,
                discharges: Vec::new(),
                section: self.current_section.clone(),
                hmr_q: q,
                hmr_op_idx: self.current_op_idx,
                hmr_mask: None,
                discharge_xor: 0,
            },
        );
        self.write_qubit(q, AbsVal::Zero);
    }

    /// Record the 64-shot sim mask of the vented qubit on its
    /// obligation, enabling sim-mask-based multi-term discharge
    /// (`accumulate_obligation_discharge` + `resolve_masked_obligation`).
    /// Called by `Circuit::hmr_ghost`.
    pub(crate) fn set_obligation_hmr_mask(&mut self, b: u32, mask: u64) {
        if !self.enabled {
            return;
        }
        if let Some(ob) = self.obligations.get_mut(&b) {
            ob.hmr_mask = Some(mask);
        }
    }

    /// XOR a discharge term's 64-shot sim mask into the obligation's
    /// running accumulator. Called by `Circuit::ghost_xor_{z,cz,ccz}`
    /// alongside the actual Z/CZ/CCZ phase gate.
    pub(crate) fn accumulate_obligation_discharge(&mut self, b: u32, mask: u64) {
        if !self.enabled {
            return;
        }
        if let Some(ob) = self.obligations.get_mut(&b) {
            ob.discharge_xor ^= mask;
        }
    }

    /// Resolve a multi-term ghost obligation: verify the accumulated
    /// discharge XOR reproduces the vented value on every shot, then
    /// remove the obligation. Returns `Err((hmr_mask, discharge_xor))`
    /// on mismatch (the caller panics with diagnostics + debugger hook).
    /// This is the sim-mask check that stands in for the structural
    /// matcher, which cannot reduce a multi-term XOR against an Anchor.
    pub(crate) fn resolve_masked_obligation(&mut self, b: u32) -> Result<(), (u64, u64)> {
        if !self.enabled {
            return Ok(());
        }
        match self.obligations.get(&b) {
            None => Ok(()), // already clean / sim disabled at hmr
            Some(ob) => {
                let hmr_mask = ob.hmr_mask.unwrap_or(0);
                if ob.discharge_xor == hmr_mask {
                    self.obligations.remove(&b);
                    Ok(())
                } else {
                    Err((hmr_mask, ob.discharge_xor))
                }
            }
        }
    }

    pub fn on_push_condition(&mut self, b: u32) {
        if !self.enabled {
            return;
        }
        self.cond_stack.push(b);
    }
    pub fn on_pop_condition(&mut self) {
        if !self.enabled {
            return;
        }
        self.cond_stack.pop();
    }

    pub fn on_bit_invert(&mut self, _b: u32) {}
    pub fn on_bit_store0(&mut self, _b: u32) {}
    pub fn on_bit_store1(&mut self, _b: u32) {}

    // ── Ghost / spooky-pebble accessors (PR-2) ──────────────────────
    //
    // See `tracker/ghost.rs`.
    // Used by `Circuit::hmr_ghost` / `Circuit::resolve_ghost` to
    // detach an HMR obligation from any specific (qubit, version)
    // pair and re-attach it at discharge time via a shared Anchor
    // id. `AbsVal::Anchor(id)` is always version-current
    // (`val_versions_current`'s `Anchor(_)` arm), so the obligation
    // survives arbitrarily many intervening ops between HMR and
    // discharge.

    /// Mint a fresh anchor id and overwrite `val[q] = Anchor(id)`.
    /// Caller (`Circuit::hmr_ghost`) immediately follows with
    /// `on_hmr(q, b)`, which captures `qval(q) = Anchor(id)` into
    /// the obligation. Returns the anchor id (also stored on the
    /// returned `Ghost` for the matching `re_anchor_for_resolve`).
    ///
    /// Does NOT bump `ver[q]` — the version bump happens inside
    /// `on_hmr`'s `write_qubit(q, Zero)` call. With the tracker
    /// disabled, returns 0 (no anchor minted; downstream
    /// `re_anchor_for_resolve` is also a no-op).
    pub(crate) fn anchor_qubit_and_get_id(&mut self, q: u32) -> u64 {
        if !self.enabled {
            return 0;
        }
        self.ensure_qubit(q);
        let id = self.next_anchor;
        self.next_anchor += 1;
        self.val[q as usize] = AbsVal::Anchor(id);
        id
    }

    /// Overwrite `val[r] = Anchor(id)` so the immediately-following
    /// `on_z(r)` deposits the matching `Anchor(id)` contribution
    /// into `obligation.discharges` (via `discharge_to_cond`).
    /// Caller (`Circuit::resolve_ghost`) is responsible for having
    /// sim-verified that `r`'s 64-shot mask matches the HMR-time
    /// mask of the anchored qubit BEFORE calling this — the anchor
    /// promotion is symbolic only.
    ///
    /// No-op when the tracker is disabled.
    pub(crate) fn re_anchor_for_resolve(&mut self, r: u32, id: u64) {
        if !self.enabled {
            return;
        }
        self.ensure_qubit(r);
        self.val[r as usize] = AbsVal::Anchor(id);
    }

    // ── Finalize ────────────────────────────────────────────────

    pub fn assert_clean(&self) {
        self.assert_clean_except(&std::collections::HashSet::new());
    }

    /// Like [`assert_clean`] but IGNORES obligations whose HMR cbit is in
    /// `skip_bits`. Used by `destroy_sim_ghosts`: if the lattice is clean
    /// except for the handed-over ghosts' obligations, then resolving those
    /// ghosts WOULD make it fully clean (so abandoning them is honest -- there
    /// is no extra, un-ghosted phase obligation hiding a real bug). The
    /// non-obligation checks (direct phase, R-on-nonzero) are NOT skipped.
    pub fn assert_clean_except(&self, skip_bits: &std::collections::HashSet<u32>) {
        use std::io::Write;
        if !self.enabled {
            return;
        }

        let mut errs: Vec<String> = Vec::new();

        if self.direct_phase_nonzero {
            errs.push(
                "direct phase (Z/CZ/NEG outside conditions) \
                is nonzero"
                    .into(),
            );
        }

        type RNonzeroBySec<'a> = HashMap<&'a String, (u64, Option<(u32, &'a AbsVal)>)>;
        if !self.r_on_nonzero_count.is_empty() {
            let mut by_sec: RNonzeroBySec = HashMap::new();
            for (sec, count) in &self.r_on_nonzero_count {
                by_sec.entry(sec).or_insert((0, None)).0 += *count;
            }
            for (sec, q, v) in &self.r_on_nonzero {
                let e = by_sec.entry(sec).or_insert((0, None));
                if e.1.is_none() {
                    e.1 = Some((*q, v));
                }
            }
            let mut ordered: Vec<_> = by_sec.into_iter().collect();
            ordered.sort_by_key(|(_, (c, _))| -(*c as i64));
            for (sec, (count, sample)) in ordered.iter().take(20) {
                let (sq, sv) = sample.unwrap_or((0, &AbsVal::Top));
                errs.push(format!(
                    "R on non-Zero qubit in [{sec}]: {count} event(s); \
                     sample q{sq} held val={sv:?}"
                ));
            }
        }

        for (bit, ob) in &self.obligations {
            if skip_bits.contains(bit) {
                continue; // this obligation belongs to a handed-over ghost
            }
            if ob.val.is_zero() {
                continue;
            }

            if ob.val.is_top() {
                errs.push(format!(
                    "b{}: [{}] q{}: HMR'd qubit had Top value",
                    bit, ob.section, ob.hmr_q
                ));
                continue;
            }

            // Step 1: pair-cancel identical valid discharges.
            let mut remaining: Vec<&(AbsVal, bool, usize)> = ob.discharges.iter().collect();
            remaining.sort_by(|a, b| format!("{:?}", a.0).cmp(&format!("{:?}", b.0)));
            let mut unpaired: Vec<&(AbsVal, bool, usize)> = Vec::new();
            let mut i = 0;
            while i < remaining.len() {
                if i + 1 < remaining.len()
                    && remaining[i].0 == remaining[i + 1].0
                    && remaining[i].1
                    && remaining[i + 1].1
                {
                    i += 2;
                } else {
                    unpaired.push(remaining[i]);
                    i += 1;
                }
            }

            if unpaired.len() == 1 && unpaired[0].0 == ob.val && unpaired[0].1 {
                continue;
            }

            if let AbsVal::XorOf(ref oa, ref ob_inner) = ob.val {
                if unpaired.len() == 2 {
                    let (d0, d1) = (&unpaired[0].0, &unpaired[1].0);
                    let (v0, v1) = (unpaired[0].1, unpaired[1].1);
                    let matched = (d0 == oa.as_ref() && d1 == ob_inner.as_ref() && v0 && v1)
                        || (d0 == ob_inner.as_ref() && d1 == oa.as_ref() && v0 && v1);
                    if matched {
                        continue;
                    }
                }
            }

            if let AbsVal::XorOf3(ref oa, ref ob_inner, ref oc) = ob.val {
                if unpaired.len() == 3 && unpaired.iter().all(|d| d.1) {
                    let leaves = [oa.as_ref(), ob_inner.as_ref(), oc.as_ref()];
                    let discharges = [&unpaired[0].0, &unpaired[1].0, &unpaired[2].0];
                    if matches_as_multiset(&leaves, &discharges) {
                        continue;
                    }
                }
            }

            let summary: Vec<String> = ob
                .discharges
                .iter()
                .map(|(v, ok, op)| {
                    format!("op{}:{:?}({})", op, v, if *ok { "ok" } else { "stale" })
                })
                .collect();
            errs.push(format!(
                "b{}: [{}] hmr_op={} q{}: obligation {:?}, {} \
                 discharge(s) [{}] — not structurally matched.",
                bit,
                ob.section,
                ob.hmr_op_idx,
                ob.hmr_q,
                ob.val,
                ob.discharges.len(),
                summary.join(", ")
            ));
        }

        if errs.is_empty() {
            return;
        }

        let stderr = std::io::stderr();
        let mut w = stderr.lock();
        writeln!(w, "[phase lattice] BUILD FAILED — {} error(s):", errs.len()).ok();
        for e in &errs {
            writeln!(w, "  {e}").ok();
        }
        w.flush().ok();
        drop(w);
        panic!("phase lattice: {} error(s)", errs.len());
    }

    fn val_versions_current(&self, v: &AbsVal) -> bool {
        match v {
            AbsVal::Zero | AbsVal::One | AbsVal::Top | AbsVal::Anchor(_) => true,
            AbsVal::CopyOf(q, ver) => self.qver(*q) == *ver,
            AbsVal::AndOf(q1, v1, q2, v2) => self.qver(*q1) == *v1 && self.qver(*q2) == *v2,
            AbsVal::AndOf3(q1, v1, q2, v2, q3, v3) => {
                self.qver(*q1) == *v1 && self.qver(*q2) == *v2 && self.qver(*q3) == *v3
            }
            AbsVal::ChooseOf(qc, vc, qa, va, qb, vb) => {
                self.qver(*qc) == *vc && self.qver(*qa) == *va && self.qver(*qb) == *vb
            }
            AbsVal::XorOf(a, b) => self.val_versions_current(a) && self.val_versions_current(b),
            AbsVal::XorOf3(a, b, c) => {
                self.val_versions_current(a)
                    && self.val_versions_current(b)
                    && self.val_versions_current(c)
            }
        }
    }
}

// ── Sound fact-injection API on Circuit ─────────────────────────────
// Lives here (inside phase_lattice.rs) so the `_unchecked` methods
// above can be private to this module. Each wrapper requires a live
// simulator, reads the claimed identity across all 64 sim shots, and
// panics on any mismatch — then (and only then) forwards the fact
// into the tracker. This is the only supported way to inject a
// symbolic fact from outside the tracker's own transfer functions.

impl crate::point_add::trailmix_port::circuit::Circuit {
    /// Assert `val(q_a)` == `val(q_b)` across all 64 sim shots, then
    /// anchor both to a fresh symbolic identity in the tracker.
    pub(crate) fn declare_identity_raw(&mut self, q_a: u32, q_b: u32) {
        let ma = self.sim_get_mask(q_a);
        let mb = self.sim_get_mask(q_b);
        if ma != mb {
            if std::env::var("DEBUG_ON_FAIL").is_ok() {
                let last_op = self.ops_truncated as usize + self.ops.len();
                let section = self.current_section.clone();
                eprintln!(
                    "[declare_identity] q{q_a} sim_mask={ma:#x} != q{q_b} sim_mask={mb:#x} (section '{section}', op #{last_op})"
                );
                eprintln!("[debugger] DEBUG_ON_FAIL set — attaching");
                let mut dbg = crate::point_add::trailmix_port::debugger::Debugger::attach(self);
                dbg.goto(last_op);
                dbg.repl();
            }
            panic!(
                "[declare_identity] q{} sim_mask={:#x} != \
                q{} sim_mask={:#x} (section '{}', op #{})",
                q_a,
                ma,
                q_b,
                mb,
                self.current_section,
                self.ops_truncated as usize + self.ops.len()
            );
        }
        self.phase.declare_identity_unchecked(q_a, q_b);
    }

    /// Assert val(q) == val(source) across all 64 sim shots, then
    /// inject CopyOf(source, ver) into the tracker.
    pub(crate) fn declare_copy_of_raw(&mut self, q: u32, source: u32) {
        let mq = self.sim_get_mask(q);
        let ms = self.sim_get_mask(source);
        assert!(
            mq == ms,
            "[declare_copy_of] q{} sim_mask={:#x} != \
            source q{} sim_mask={:#x} (section '{}', op #{})",
            q,
            mq,
            source,
            ms,
            self.current_section,
            self.ops_truncated as usize + self.ops.len()
        );
        self.phase.declare_copy_of_unchecked(q, source);
    }

    /// Assert val(q) == val(a) AND val(b) across all 64 sim shots,
    /// then inject AndOf(a, b) into the tracker.
    pub(crate) fn declare_and_of_raw(&mut self, q: u32, a: u32, b: u32) {
        let mq = self.sim_get_mask(q);
        let ma = self.sim_get_mask(a);
        let mb = self.sim_get_mask(b);
        let expected = ma & mb;
        assert!(
            mq == expected,
            "[declare_and_of] q{} sim_mask={:#x} != \
            (q{} AND q{}) sim_mask={:#x} \
            (section '{}', op #{})",
            q,
            mq,
            a,
            b,
            expected,
            self.current_section,
            self.ops_truncated as usize + self.ops.len()
        );
        self.phase.declare_and_of_unchecked(q, a, b);
    }

    /// Assert val(q) == val(a) AND val(b) AND val(c) across all 64
    /// sim shots, then inject AndOf3(a, b, c) into the tracker. Use
    /// this before HMR-freeing a qubit computed as a 3-way AND (e.g.
    /// for discharge via `ccz_if_bit(a, b, c, hmr_bit)`).
    pub(crate) fn declare_and3_of_raw(&mut self, q: u32, a: u32, b: u32, c: u32) {
        let mq = self.sim_get_mask(q);
        let ma = self.sim_get_mask(a);
        let mb = self.sim_get_mask(b);
        let mc = self.sim_get_mask(c);
        let expected = ma & mb & mc;
        assert!(
            mq == expected,
            "[declare_and3_of] q{} sim_mask={:#x} != \
            (q{} AND q{} AND q{}) sim_mask={:#x} \
            (section '{}', op #{})",
            q,
            mq,
            a,
            b,
            c,
            expected,
            self.current_section,
            self.ops_truncated as usize + self.ops.len()
        );
        self.phase.declare_and3_of_unchecked(q, a, b, c);
    }

    /// Assert val(target) == val(q1) XOR val(q2) across all 64 sim shots,
    /// then inject `XorOf(CopyOf(q1, v1), CopyOf(q2, v2))` into the tracker.
    /// See `declare_xor_of_unchecked` for the version-tracking semantics
    /// and caller obligations.
    pub(crate) fn declare_xor_of_raw(&mut self, target: u32, q1: u32, q2: u32) {
        let mt = self.sim_get_mask(target);
        let m1 = self.sim_get_mask(q1);
        let m2 = self.sim_get_mask(q2);
        let expected = m1 ^ m2;
        if mt != expected {
            if std::env::var("DEBUG_ON_FAIL").is_ok() {
                let last_op = self.ops_truncated as usize + self.ops.len();
                let section = self.current_section.clone();
                eprintln!(
                    "[declare_xor_of] q{target} sim_mask={mt:#x} != \
                    (q{q1} XOR q{q2}) sim_mask={expected:#x} (section '{section}', op #{last_op})"
                );
                eprintln!("[debugger] DEBUG_ON_FAIL set — attaching");
                let mut dbg = crate::point_add::trailmix_port::debugger::Debugger::attach(self);
                dbg.goto(last_op);
                dbg.repl();
            }
            panic!(
                "[declare_xor_of] q{} sim_mask={:#x} != \
                (q{} XOR q{}) sim_mask={:#x} (section '{}', op #{})",
                target,
                mt,
                q1,
                q2,
                expected,
                self.current_section,
                self.ops_truncated as usize + self.ops.len()
            );
        }
        self.phase.declare_xor_of_unchecked(target, q1, q2);
    }

    /// Assert val(target) == val(q1) XOR val(q2) XOR val(q3) across all
    /// 64 sim shots, then inject `XorOf3`(...) into the tracker.
    pub(crate) fn declare_xor_of_three_raw(&mut self, target: u32, q1: u32, q2: u32, q3: u32) {
        let mt = self.sim_get_mask(target);
        let m1 = self.sim_get_mask(q1);
        let m2 = self.sim_get_mask(q2);
        let m3 = self.sim_get_mask(q3);
        let expected = m1 ^ m2 ^ m3;
        if mt != expected {
            if std::env::var("DEBUG_ON_FAIL").is_ok() {
                let last_op = self.ops_truncated as usize + self.ops.len();
                let section = self.current_section.clone();
                eprintln!(
                    "[declare_xor_of_three] q{target} sim_mask={mt:#x} != \
                    (q{q1} XOR q{q2} XOR q{q3}) sim_mask={expected:#x} \
                    (section '{section}', op #{last_op})"
                );
                eprintln!("[debugger] DEBUG_ON_FAIL set — attaching");
                let mut dbg = crate::point_add::trailmix_port::debugger::Debugger::attach(self);
                dbg.goto(last_op);
                dbg.repl();
            }
            panic!(
                "[declare_xor_of_three] q{} sim_mask={:#x} != \
                (q{} XOR q{} XOR q{}) sim_mask={:#x} \
                (section '{}', op #{})",
                target,
                mt,
                q1,
                q2,
                q3,
                expected,
                self.current_section,
                self.ops_truncated as usize + self.ops.len()
            );
        }
        self.phase
            .declare_xor_of_three_unchecked(target, q1, q2, q3);
    }

    /// Assert val(target) == (ctrl ? a : b) across all 64 sim shots,
    /// then inject ChooseOf(ctrl, a, b) into the tracker. Used by the
    /// Luo-pack passthrough-shift gadgets.
    pub(crate) fn declare_choose_of_raw(&mut self, target: u32, ctrl: u32, a: u32, b: u32) {
        let mt = self.sim_get_mask(target);
        let mc = self.sim_get_mask(ctrl);
        let ma = self.sim_get_mask(a);
        let mb = self.sim_get_mask(b);
        let expected = (mc & ma) | ((!mc) & mb);
        assert!(
            mt == expected,
            "[declare_choose_of] q{} sim_mask={:#x} != \
            (q{} ? q{} : q{}) sim_mask={:#x} \
            (section '{}', op #{})",
            target,
            mt,
            ctrl,
            a,
            b,
            expected,
            self.current_section,
            self.ops_truncated as usize + self.ops.len()
        );
        self.phase.declare_choose_of_unchecked(target, ctrl, a, b);
    }

    /// Concrete-backed proof that q is |0> across all 64 sim shots.
    /// On success, inject val(q) = Zero into the tracker so downstream
    /// R (or free) doesn't flag a false positive. `DEBUG_ON_FAIL=1`
    /// drops into the debugger REPL at the failing op.
    /// Internal-only proof-of-zero by raw id. The `&QReg`-taking
    /// public wrapper lives on `Circuit` in `circuit.rs` where it can
    /// access the `QReg`'s id field directly.
    pub(crate) fn prove_zero_raw(&mut self, qr: u32) {
        let mask = self.sim_get_mask(qr);
        if mask != 0 {
            let alloc_tag = self
                .qubit_alloc_log
                .iter()
                .rev()
                .find(|ev| ev.qubit == qr)
                .map_or("<unknown>", |ev| &*ev.tag);
            let msg = format!(
                "[prove_zero] q{} expected |0> in all 64 shots \
                but sim_mask = {:#x} tag='{}' (section '{}', op #{})",
                qr,
                mask,
                alloc_tag,
                self.current_section,
                self.ops_truncated as usize + self.ops.len()
            );
            if std::env::var("DEBUG_ON_FAIL").is_ok() {
                eprintln!("{msg}");
                eprintln!("[debugger] DEBUG_ON_FAIL set — attaching");
                // Capture target BEFORE attach (attach moves self.ops).
                let target = self.ops_truncated as usize + self.ops.len();
                let mut d = crate::point_add::trailmix_port::debugger::Debugger::attach(self);
                d.goto(target);
                d.repl();
            }
            panic!("{}", msg);
        }
        self.phase.inject_zero_unchecked(qr);
    }

    /// Pure sim assertion that `qr` holds classical `expected` (0 or 1)
    /// uniformly across all 64 shots. Does NOT modify the tracker or
    /// emit any ops -- this is a check, not an injection. `DEBUG_ON_FAIL=1`
    /// drops into the debugger REPL at the failing op so callers can
    /// inspect the offending state.
    pub(crate) fn assert_qubit_eq_raw(&mut self, qr: u32, expected: u8) {
        let mask = self.sim_get_mask(qr);
        let want = if expected & 1 == 1 { !0u64 } else { 0u64 };
        if mask != want {
            let alloc_tag = self
                .qubit_alloc_log
                .iter()
                .rev()
                .find(|ev| ev.qubit == qr)
                .map_or("<unknown>", |ev| &*ev.tag);
            let msg = format!(
                "[assert_qubit_eq] q{} expected={} (want_mask={:#x}) but \
                 sim_mask={:#x} tag='{}' (section '{}', op #{})",
                qr,
                expected,
                want,
                mask,
                alloc_tag,
                self.current_section,
                self.ops_truncated as usize + self.ops.len(),
            );
            if std::env::var("DEBUG_ON_FAIL").is_ok() {
                eprintln!("{msg}");
                eprintln!("[debugger] DEBUG_ON_FAIL set — attaching");
                // Capture target BEFORE attach (attach moves self.ops).
                let target = self.ops_truncated as usize + self.ops.len();
                let mut d = crate::point_add::trailmix_port::debugger::Debugger::attach(self);
                d.goto(target);
                d.repl();
            }
            panic!("{msg}");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mbuc_pattern_clean() {
        let mut t = PhaseTracker::new();
        t.mark_input_qubit(1);
        t.mark_input_qubit(2);
        t.on_alloc_qubit(3);
        t.on_ccx(1, 2, 3);
        t.on_alloc_bit(10);
        t.on_hmr(3, 10);
        t.on_push_condition(10);
        t.on_cz(1, 2);
        t.on_pop_condition();
        t.assert_clean();
    }

    #[test]
    #[should_panic(expected = "error")]
    fn bare_hmr_leaks() {
        let mut t = PhaseTracker::new();
        t.mark_input_qubit(1);
        t.on_alloc_bit(5);
        t.on_hmr(1, 5);
        t.assert_clean();
    }

    #[test]
    fn copy_then_hmr_then_z_if_bit() {
        let mut t = PhaseTracker::new();
        t.mark_input_qubit(1);
        t.on_alloc_qubit(2);
        t.on_cx(1, 2);
        t.on_alloc_bit(10);
        t.on_hmr(2, 10);
        t.on_push_condition(10);
        t.on_z(1);
        t.on_pop_condition();
        t.assert_clean();
    }

    #[test]
    fn declare_identity_then_hmr_and_z() {
        let mut t = PhaseTracker::new();
        t.mark_input_qubit(0);
        t.mark_input_qubit(256);
        // Unit test of the tracker in isolation — no Circuit to run a
        // sim check, so use the _unchecked variant directly.
        t.declare_identity_unchecked(0, 256);
        t.on_alloc_bit(10);
        t.on_hmr(256, 10);
        t.on_push_condition(10);
        t.on_z(0);
        t.on_pop_condition();
        t.assert_clean();
    }

    #[test]
    fn declare_choose_of_records_versions() {
        // Build minimal tracker state: 4 qubits, all marked as inputs
        // so they are Top. After declare_choose_of, target's val is
        // ChooseOf with the version tags pinned at injection time.
        let mut t = PhaseTracker::new();
        t.mark_input_qubit(0); // ctrl
        t.mark_input_qubit(1); // a
        t.mark_input_qubit(2); // b
        t.on_alloc_qubit(3); // target
        t.declare_choose_of_unchecked(3, 0, 1, 2);
        match &t.val[3] {
            AbsVal::ChooseOf(qc, _, qa, _, qb, _) => {
                assert_eq!(*qc, 0);
                assert_eq!(*qa, 1);
                assert_eq!(*qb, 2);
            }
            other => panic!("expected ChooseOf, got {:?}", other),
        }
        // val_versions_current should be true initially.
        let v = t.val[3].clone();
        assert!(t.val_versions_current(&v));
        // Now bump version of `a` (qubit 1) -- choose-of should become
        // version-stale.
        t.write_qubit(1, AbsVal::Top);
        let v2 = t.val[3].clone();
        assert!(!t.val_versions_current(&v2));
    }

    #[test]
    fn declare_choose_of_constant_ctrl_simplifies() {
        // ctrl = One -> choose simplifies to copy_of_or_val(a).
        let mut t = PhaseTracker::new();
        t.on_alloc_qubit(0); // ctrl, starts Zero
        t.on_x(0); // -> One
        t.mark_input_qubit(1); // a (Top)
        t.mark_input_qubit(2); // b (Top)
        t.on_alloc_qubit(3);
        t.declare_choose_of_unchecked(3, 0, 1, 2);
        match &t.val[3] {
            AbsVal::CopyOf(q, _) => assert_eq!(*q, 1),
            other => panic!("expected CopyOf when ctrl=One, got {:?}", other),
        }
    }

    #[test]
    #[should_panic(expected = "phase lattice")]
    fn r_on_nonzero_is_caught() {
        let mut t = PhaseTracker::new();
        t.mark_input_qubit(1);
        t.on_r(1);
        t.assert_clean();
    }

    // ── Ghost-anchor accessor tests (PR-2) ─────────────────────────

    #[test]
    fn anchor_ids_are_distinct_per_qubit() {
        let mut t = PhaseTracker::new();
        t.mark_input_qubit(1);
        t.mark_input_qubit(2);
        let id_a = t.anchor_qubit_and_get_id(1);
        let id_b = t.anchor_qubit_and_get_id(2);
        assert_ne!(id_a, id_b);
        // And subsequent calls also mint fresh ids.
        let id_c = t.anchor_qubit_and_get_id(1);
        assert_ne!(id_a, id_c);
        assert_ne!(id_b, id_c);
        // val[1] now holds the most recent anchor.
        match &t.val[1] {
            AbsVal::Anchor(id) => assert_eq!(*id, id_c),
            other => panic!("expected Anchor, got {:?}", other),
        }
    }

    #[test]
    fn re_anchor_for_resolve_updates_val() {
        let mut t = PhaseTracker::new();
        t.on_alloc_qubit(7);
        assert!(matches!(t.val[7], AbsVal::Zero));
        t.re_anchor_for_resolve(7, 42);
        match &t.val[7] {
            AbsVal::Anchor(id) => assert_eq!(*id, 42),
            other => panic!("expected Anchor(42), got {:?}", other),
        }
    }

    #[test]
    fn declare_xor_of_then_hmr_and_two_z_discharge() {
        // target = q1 XOR q2; HMR target; z_if_bit(q1); z_if_bit(q2);
        // assert_clean. Matches the snapshot-discharge pattern used by
        // mod_mul_fused_v2's cma:phase replacement.
        let mut t = PhaseTracker::new();
        t.mark_input_qubit(1);
        t.mark_input_qubit(2);
        t.on_alloc_qubit(3); // target
                             // Unit tests use the _unchecked variant directly (no sim).
        t.declare_xor_of_unchecked(3, 1, 2);
        t.on_alloc_bit(10);
        t.on_hmr(3, 10);
        t.on_push_condition(10);
        t.on_z(1);
        t.on_z(2);
        t.on_pop_condition();
        t.on_free_bit(10);
        t.assert_clean();
    }

    #[test]
    #[should_panic(expected = "error")]
    fn declare_xor_of_stale_version_leaks() {
        // Tracker discovers a stale leaf: after declare_xor_of(target, q1,
        // q2), modify q1 (bumps version). z_if_bit(q1) deposits a leaf at
        // a NEW version that doesn't match the obligation's stored leaf.
        // assert_clean must report the unmatched obligation.
        let mut t = PhaseTracker::new();
        t.mark_input_qubit(1);
        t.mark_input_qubit(2);
        t.on_alloc_qubit(3);
        t.declare_xor_of_unchecked(3, 1, 2);
        t.on_alloc_bit(10);
        t.on_hmr(3, 10);
        // BUG: modify q1 between declare and discharge — version stale.
        t.on_x(1);
        t.on_push_condition(10);
        t.on_z(1);
        t.on_z(2);
        t.on_pop_condition();
        t.assert_clean();
    }

    #[test]
    #[should_panic(expected = "error")]
    fn declare_xor_of_missing_discharge_leaks() {
        // Only one z_if_bit fires; the XorOf obligation needs two.
        let mut t = PhaseTracker::new();
        t.mark_input_qubit(1);
        t.mark_input_qubit(2);
        t.on_alloc_qubit(3);
        t.declare_xor_of_unchecked(3, 1, 2);
        t.on_alloc_bit(10);
        t.on_hmr(3, 10);
        t.on_push_condition(10);
        t.on_z(1);
        t.on_pop_condition();
        t.assert_clean();
    }

    #[test]
    fn declare_xor_of_three_then_hmr_and_three_z_discharge() {
        // target = q1 XOR q2 XOR q3; HMR target; z_if_bit each;
        // assert_clean. End-to-end XorOf3 happy path.
        let mut t = PhaseTracker::new();
        t.mark_input_qubit(1);
        t.mark_input_qubit(2);
        t.mark_input_qubit(3);
        t.on_alloc_qubit(4);
        t.declare_xor_of_three_unchecked(4, 1, 2, 3);
        t.on_alloc_bit(10);
        t.on_hmr(4, 10);
        t.on_push_condition(10);
        t.on_z(1);
        t.on_z(2);
        t.on_z(3);
        t.on_pop_condition();
        t.on_free_bit(10);
        t.assert_clean();
    }

    #[test]
    fn declare_xor_of_three_accepts_permuted_discharges() {
        // The XorOf3 matcher must accept discharges in any order.
        let mut t = PhaseTracker::new();
        t.mark_input_qubit(1);
        t.mark_input_qubit(2);
        t.mark_input_qubit(3);
        t.on_alloc_qubit(4);
        t.declare_xor_of_three_unchecked(4, 1, 2, 3);
        t.on_alloc_bit(10);
        t.on_hmr(4, 10);
        t.on_push_condition(10);
        // Different order than the declare.
        t.on_z(3);
        t.on_z(1);
        t.on_z(2);
        t.on_pop_condition();
        t.on_free_bit(10);
        t.assert_clean();
    }

    #[test]
    #[should_panic(expected = "error")]
    fn declare_xor_of_three_partial_discharge_leaks() {
        // Two of three leaves discharged; obligation not matched.
        let mut t = PhaseTracker::new();
        t.mark_input_qubit(1);
        t.mark_input_qubit(2);
        t.mark_input_qubit(3);
        t.on_alloc_qubit(4);
        t.declare_xor_of_three_unchecked(4, 1, 2, 3);
        t.on_alloc_bit(10);
        t.on_hmr(4, 10);
        t.on_push_condition(10);
        t.on_z(1);
        t.on_z(2);
        t.on_pop_condition();
        t.assert_clean();
    }

    /// End-to-end tracker dance for the Ghost API: HMR with anchor
    /// promotion, then a deferred discharge on a re-anchored target,
    /// then assert_clean.
    #[test]
    fn ghost_anchor_discharge_matches_obligation() {
        let mut t = PhaseTracker::new();
        t.mark_input_qubit(1); // input "value-bearing" qubit q
        t.mark_input_qubit(2); // a logically-equal "witness" qubit r,
                               // computed later
                               // -- HMR side: anchor q, then HMR it on a fresh bit.
        let anchor_id = t.anchor_qubit_and_get_id(1);
        t.on_alloc_bit(10);
        t.on_hmr(1, 10);
        // Many intervening ops would happen here in real code.

        // -- Discharge side: re-anchor r to the same id, then
        //    z_if_bit(r, 10) := push_condition(10), on_z(r), pop.
        t.re_anchor_for_resolve(2, anchor_id);
        t.on_push_condition(10);
        t.on_z(2);
        t.on_pop_condition();
        t.on_free_bit(10);
        t.assert_clean();
    }

    #[test]
    #[should_panic(expected = "error")]
    fn ghost_anchor_without_discharge_leaks() {
        // Anchor + HMR but no discharge: assert_clean must report the
        // unmatched obligation. This is the "Ghost dropped without
        // resolve" failure mode at the tracker level.
        let mut t = PhaseTracker::new();
        t.mark_input_qubit(1);
        t.anchor_qubit_and_get_id(1);
        t.on_alloc_bit(10);
        t.on_hmr(1, 10);
        t.assert_clean();
    }
}
