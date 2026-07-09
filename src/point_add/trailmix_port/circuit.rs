use std::cell::RefCell;
use std::fmt;
use std::ops::Deref;
use std::rc::Rc;

use crate::circuit::{BitId, Op, OperationType, QubitId};
use crate::point_add::B;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash, PartialOrd, Ord)]
pub struct Cbit(pub u32);

impl Cbit {
    #[inline]
    pub fn raw(self) -> u32 {
        self.0
    }
}

#[derive(Debug)]
pub struct QReg {
    id: u32,
    pending: Rc<RefCell<Vec<u32>>>,
    detached: bool,
}

impl QReg {
    #[inline]
    pub(crate) fn id(&self) -> u32 {
        self.id
    }

    /// Create a non-owning register view with the same physical qubit id.
    /// Dropping the view must not return the underlying lane to the allocator.
    pub(crate) fn borrowed_alias(&self) -> Self {
        Self {
            id: self.id,
            pending: Rc::clone(&self.pending),
            detached: true,
        }
    }
}

impl Drop for QReg {
    fn drop(&mut self) {
        if !self.detached {
            self.pending.borrow_mut().push(self.id);
        }
    }
}

#[derive(Debug)]
pub enum BorrowedQReg<'a> {
    Owned(QReg),
    Borrowed(&'a QReg),
}

impl Deref for BorrowedQReg<'_> {
    type Target = QReg;

    fn deref(&self) -> &Self::Target {
        match self {
            BorrowedQReg::Owned(q) => q,
            BorrowedQReg::Borrowed(q) => q,
        }
    }
}

pub struct Ghost {
    bit: Cbit,
    consumed: bool,
}

impl Drop for Ghost {
    fn drop(&mut self) {
        if !self.consumed && !std::thread::panicking() {
            panic!("TrailMix ghost dropped without a matching close/resolve");
        }
    }
}

#[derive(Clone, Copy)]
pub struct ContractSimView<'a> {
    _phantom: std::marker::PhantomData<&'a ()>,
}

impl ContractSimView<'_> {
    pub fn read_bit_shot(&self, _q: &QReg, _shot: usize) -> bool {
        false
    }

    pub fn qubit_mask(&self, _q: &QReg) -> u64 {
        0
    }

    pub fn bit_mask(&self, _id: u32) -> u64 {
        0
    }

    pub fn phase_mask(&self) -> u64 {
        0
    }

    pub fn read_bytes_shot(&self, reg: &[QReg], _shot: usize) -> Vec<u8> {
        vec![0; reg.len().div_ceil(8)]
    }

    pub fn contract_read_u256_shot(
        &self,
        _reg: &[QReg],
        _shot: usize,
    ) -> crate::point_add::trailmix_port::num_bigint::BigUint {
        crate::point_add::trailmix_port::num_bigint::BigUint::default()
    }

    pub fn contract_read_bit_shot(&self, _q: &QReg, _shot: usize) -> bool {
        false
    }
}

pub trait ContractReadable {
    fn contract_read_u256_shot(
        &self,
        reg: &[QReg],
        shot: usize,
    ) -> crate::point_add::trailmix_port::num_bigint::BigUint;
    fn contract_read_bit_shot(&self, q: &QReg, shot: usize) -> bool;
}

impl ContractReadable for ContractSimView<'_> {
    fn contract_read_u256_shot(
        &self,
        reg: &[QReg],
        shot: usize,
    ) -> crate::point_add::trailmix_port::num_bigint::BigUint {
        self.contract_read_u256_shot(reg, shot)
    }

    fn contract_read_bit_shot(&self, q: &QReg, shot: usize) -> bool {
        self.contract_read_bit_shot(q, shot)
    }
}

pub struct DestroyedSimState;

impl DestroyedSimState {
    pub fn qubit_mask(&self, _q: &QReg) -> u64 {
        0
    }

    pub fn read_bit_shot(&self, _q: &QReg, _shot: usize) -> u8 {
        0
    }

    pub fn read_bytes_shot(&self, reg: &[QReg], _shot: usize) -> Vec<u8> {
        vec![0; reg.len().div_ceil(8)]
    }

    pub fn phase_mask(&self) -> u64 {
        0
    }

    pub fn bit_mask(&self, _id: u32) -> u64 {
        0
    }
}

pub struct Circuit {
    pub b: B,
    pub current_section: String,
    pub(crate) lowq_passenger_top_releases: u32,
    pub(crate) lowq_lambda_top_releases: u32,
    pub(crate) lowq_trace_register_widths: [usize; 7],
    pub(crate) lowq_trace_division_borrow: bool,
    pub(crate) lowq_trace_multiply_borrow: bool,
    pub(crate) lowq_q948_direct_hclz_peak_guard_active: bool,
    pending_frees: Rc<RefCell<Vec<u32>>>,
    named_peak_targets: Vec<u32>,
    live_qreg_names: std::collections::BTreeMap<u32, String>,
    next_ghost_id: u64,
}

impl Circuit {
    pub fn new() -> Self {
        Self::new_with_ops_capacity(0)
    }

    pub fn new_with_ops_capacity(ops_capacity: usize) -> Self {
        let b = if std::env::var("POINT_ADD_COUNT_ONLY").ok().as_deref() == Some("1") {
            B::new_count_only()
        } else if ops_capacity == 0 {
            B::new()
        } else {
            B::new_with_ops_capacity(ops_capacity)
        };
        Self {
            b,
            current_section: "trailmix".to_string(),
            lowq_passenger_top_releases: 0,
            lowq_lambda_top_releases: 0,
            lowq_trace_register_widths: [0; 7],
            lowq_trace_division_borrow: false,
            lowq_trace_multiply_borrow: false,
            lowq_q948_direct_hclz_peak_guard_active: false,
            pending_frees: Rc::new(RefCell::new(Vec::new())),
            named_peak_targets: Self::named_peak_targets_from_env(),
            live_qreg_names: std::collections::BTreeMap::new(),
            next_ghost_id: 0,
        }
    }

    pub fn into_builder(mut self) -> B {
        self.flush_pending_frees();
        self.b
    }

    pub fn total_ops(&self) -> u64 {
        self.b.current_ops_len() as u64
    }

    pub fn set_max_qubit_peak(&mut self, _peak: u32) {}

    pub fn contracts_enabled(&self) -> bool {
        false
    }

    pub fn contract_view(&self) -> ContractSimView<'_> {
        ContractSimView {
            _phantom: std::marker::PhantomData,
        }
    }

    pub fn contract_check<F>(&mut self, _label: &str, _check: F)
    where
        F: for<'a> FnMut(ContractSimView<'a>, usize) -> Result<(), String>,
    {
    }

    pub fn contract_capture<T, F>(&mut self, _label: &str, _pre: F)
    where
        F: for<'a> FnMut(ContractSimView<'a>, usize) -> Result<T, String>,
    {
    }

    pub fn contract_pop_and_check<T, F>(&mut self, _label: &str, _post: F)
    where
        F: for<'a> FnMut(&T, ContractSimView<'a>, usize) -> Result<(), String>,
    {
    }

    pub fn sim_load_reg_bytes_shot(&mut self, _reg: &[QReg], _bytes: &[u8], _shot: usize) {}

    pub fn sim_load_bits_bytes_shot(&mut self, _bits: &[Cbit], _bytes: &[u8], _shot: usize) {}

    pub fn assert_phase_clean(&self) {}

    pub fn destroy_sim(&mut self, outputs: Vec<QReg>) -> (DestroyedSimState, Vec<QReg>) {
        let outputs = outputs
            .into_iter()
            .map(|mut q| {
                q.detached = true;
                q
            })
            .collect();
        (DestroyedSimState, outputs)
    }

    pub fn alloc_qreg(&mut self, name: &str) -> QReg {
        self.flush_pending_frees();
        let previous_peak = self.b.peak_qubits;
        let q = self.b.alloc_qubit();
        if !self.named_peak_targets.is_empty() {
            assert!(
                self.live_qreg_names
                    .insert(q.0 as u32, name.to_owned())
                    .is_none(),
                "named peak trace reused a live physical qubit"
            );
            self.record_named_peak_plateau(name);
        }
        if std::env::var("TRACE_PEAK_NAMES").ok().as_deref() == Some("1")
            && self.b.peak_qubits > previous_peak
            && self.b.peak_qubits >= 800
        {
            eprintln!(
                "PEAK_NAME active={} phase={} ops_idx={} allocation={}",
                self.b.peak_qubits, self.b.peak_phase, self.b.peak_ops_idx, name
            );
        }
        QReg {
            id: q.0 as u32,
            pending: Rc::clone(&self.pending_frees),
            detached: false,
        }
    }

    pub fn alloc_qreg_bits(&mut self, name: &str, n: usize) -> Vec<QReg> {
        (0..n)
            .map(|i| self.alloc_qreg(&format!("{name}[{i}]")))
            .collect()
    }

    pub fn alloc_input_qreg_bits(&mut self, name: &str, n: usize) -> Vec<QReg> {
        self.alloc_qreg_bits(name, n)
    }

    pub fn alloc_bit(&mut self) -> Cbit {
        Cbit(self.b.alloc_bit().0 as u32)
    }

    pub fn alloc_input_bit(&mut self) -> Cbit {
        self.alloc_bit()
    }

    pub fn free_bit(&mut self, _b: Cbit) {}

    pub fn flush_pending_frees(&mut self) {
        let pending: Vec<u32> = self.pending_frees.borrow_mut().drain(..).collect();
        for q in pending {
            if !self.named_peak_targets.is_empty() {
                self.live_qreg_names
                    .remove(&q)
                    .expect("named peak trace freed an unnamed qubit");
            }
            self.b.free(QubitId(q.into()));
        }
    }

    pub fn zero_and_free(&mut self, mut q: QReg) {
        self.flush_pending_frees();
        if !self.named_peak_targets.is_empty() {
            self.live_qreg_names
                .remove(&q.id)
                .expect("named peak trace zero-freed an unnamed qubit");
        }
        self.b.free(QubitId(q.id.into()));
        q.detached = true;
    }

    fn named_component(name: &str) -> String {
        let Some(open) = name.rfind('[') else {
            return name.to_owned();
        };
        let Some(index) = name
            .strip_suffix(']')
            .and_then(|value| value.get(open + 1..))
        else {
            return name.to_owned();
        };
        if index.is_empty() || !index.bytes().all(|value| value.is_ascii_digit()) {
            return name.to_owned();
        }
        name[..open].to_owned()
    }

    fn named_peak_targets_from_env() -> Vec<u32> {
        let mut targets = Vec::new();
        if let Ok(value) = std::env::var("TRACE_NAMED_PEAK_TARGET") {
            targets.push(
                value
                    .parse::<u32>()
                    .expect("TRACE_NAMED_PEAK_TARGET integer"),
            );
        }
        if let Ok(value) = std::env::var("TRACE_NAMED_PEAK_TARGETS") {
            for item in value.split(',') {
                let item = item.trim();
                if item.is_empty() {
                    continue;
                }
                targets.push(
                    item.parse::<u32>()
                        .expect("TRACE_NAMED_PEAK_TARGETS comma-separated integers"),
                );
            }
        }
        if std::env::var("TRACE_LOWQ_OCCUPANCY_REPORT").ok().as_deref() == Some("1") {
            targets.extend([811, 823, 824]);
        }
        targets.sort_unstable();
        targets.dedup();
        targets
    }

    fn record_named_peak_plateau(&mut self, trigger_allocation: &str) {
        let target = self.b.active_qubits;
        if !self.named_peak_targets.contains(&target) {
            return;
        }
        assert_eq!(
            self.live_qreg_names.len(),
            self.b.active_qubits as usize,
            "named peak trace does not cover every live qubit"
        );
        let mut components = std::collections::BTreeMap::<String, usize>::new();
        for name in self.live_qreg_names.values() {
            *components.entry(Self::named_component(name)).or_default() += 1;
        }
        let live_components: Vec<_> = components.into_iter().collect();
        assert_eq!(
            live_components
                .iter()
                .map(|(_, lanes)| lanes)
                .sum::<usize>(),
            target as usize
        );
        let trigger_component = Self::named_component(trigger_allocation);
        let allocation_serial = self.b.allocation_serial;
        let ops_idx = self.b.current_ops_len();
        if let Some(plateau) = self.b.named_peak_plateaus.iter_mut().find(|plateau| {
            plateau.target == target
                && plateau.phase == self.b.phase
                && plateau.trigger_allocation == trigger_allocation
                && plateau.live_components == live_components
        }) {
            plateau.occurrences += 1;
            plateau.last_allocation_serial = allocation_serial;
            plateau.last_ops_idx = ops_idx;
            return;
        }
        self.b
            .named_peak_plateaus
            .push(crate::point_add::NamedPeakPlateau {
                target,
                phase: self.b.phase,
                trigger_allocation: trigger_allocation.to_owned(),
                trigger_component,
                occurrences: 1,
                first_allocation_serial: allocation_serial,
                last_allocation_serial: allocation_serial,
                first_ops_idx: ops_idx,
                last_ops_idx: ops_idx,
                live_components,
            });
    }

    pub fn register(&mut self, id: u32) {
        while self.b.next_register < id {
            self.b.next_register += 1;
        }
        let old = self.b.next_register;
        self.b.next_register = id;
        let mut op = Op::empty();
        op.kind = OperationType::Register;
        op.r_target = crate::circuit::RegisterId(id.into());
        self.b.push_op(op);
        self.b.next_register = old.max(id + 1);
    }

    pub fn append_qreg(&mut self, q: &QReg, reg: u32) {
        let mut op = Op::empty();
        op.kind = OperationType::AppendToRegister;
        op.q_target = QubitId(q.id.into());
        op.r_target = crate::circuit::RegisterId(reg.into());
        self.b.push_op(op);
    }

    pub fn append_bit(&mut self, bit: Cbit, reg: u32) {
        let mut op = Op::empty();
        op.kind = OperationType::AppendToRegister;
        op.c_target = BitId(bit.0.into());
        op.r_target = crate::circuit::RegisterId(reg.into());
        self.b.push_op(op);
    }

    pub fn declare_registers(&mut self, tx: &[QReg], ty: &[QReg], ox: &[Cbit], oy: &[Cbit]) {
        self.register(0);
        for q in tx {
            self.append_qreg(q, 0);
        }
        self.register(1);
        for q in ty {
            self.append_qreg(q, 1);
        }
        self.register(2);
        for &b in ox {
            self.append_bit(b, 2);
        }
        self.register(3);
        for &b in oy {
            self.append_bit(b, 3);
        }
    }

    pub fn defragment(&mut self, mut slots: Vec<QReg>) -> Vec<QReg> {
        self.flush_pending_frees();
        let n = slots.len();
        for v in 0..n {
            let want = v as u32;
            if slots[v].id == want {
                continue;
            }
            if let Some(w) = slots.iter().position(|q| q.id == want) {
                self.swap(&slots[v], &slots[w]);
                slots.swap(v, w);
            } else {
                self.b.reacquire(QubitId(want.into()));
                let lo = QReg {
                    id: want,
                    pending: Rc::clone(&self.pending_frees),
                    detached: false,
                };
                self.swap(&lo, &slots[v]);
                let old = std::mem::replace(&mut slots[v], lo);
                self.zero_and_free(old);
            }
        }
        slots
    }

    pub fn x(&mut self, q: &QReg) {
        self.flush_pending_frees();
        self.b.x(QubitId(q.id.into()));
    }

    pub fn z(&mut self, q: &QReg) {
        self.flush_pending_frees();
        let mut op = Op::empty();
        op.kind = OperationType::Z;
        op.q_target = QubitId(q.id.into());
        self.b.push_op(op);
    }

    pub fn cx(&mut self, ctrl: &QReg, tgt: &QReg) {
        self.flush_pending_frees();
        self.b.cx(QubitId(ctrl.id.into()), QubitId(tgt.id.into()));
    }

    pub fn cz(&mut self, a: &QReg, b: &QReg) {
        self.flush_pending_frees();
        self.b.cz(QubitId(a.id.into()), QubitId(b.id.into()));
    }

    pub fn ccx(&mut self, a: &QReg, b: &QReg, t: &QReg) {
        self.flush_pending_frees();
        self.b.ccx(
            QubitId(a.id.into()),
            QubitId(b.id.into()),
            QubitId(t.id.into()),
        );
    }

    pub fn ccz(&mut self, a: &QReg, b: &QReg, c: &QReg) {
        self.flush_pending_frees();
        let mut op = Op::empty();
        op.kind = OperationType::CCZ;
        op.q_control2 = QubitId(a.id.into());
        op.q_control1 = QubitId(b.id.into());
        op.q_target = QubitId(c.id.into());
        self.b.push_op(op);
    }

    pub fn swap(&mut self, a: &QReg, b: &QReg) {
        self.flush_pending_frees();
        self.b.swap(QubitId(a.id.into()), QubitId(b.id.into()));
    }

    pub fn cswap(&mut self, ctrl: &QReg, a: &QReg, b: &QReg) {
        self.cx(b, a);
        self.ccx(ctrl, a, b);
        self.cx(b, a);
    }

    pub fn hmr(&mut self, q: &QReg, bit: Cbit) {
        self.flush_pending_frees();
        self.b.hmr(QubitId(q.id.into()), BitId(bit.0.into()));
    }

    pub fn hmr_ghost(&mut self, q: &QReg) -> Ghost {
        let bit = self.alloc_bit();
        self.hmr(q, bit);
        self.next_ghost_id += 1;
        Ghost {
            bit,
            consumed: false,
        }
    }

    pub fn resolve_ghost(&mut self, mut g: Ghost, r: &QReg) {
        self.z_if_bit(r, g.bit);
        self.free_bit(g.bit);
        g.consumed = true;
    }

    pub fn ghost_xor_z(&mut self, g: &mut Ghost, r: &QReg) {
        self.z_if_bit(r, g.bit);
    }

    pub fn ghost_xor_cz(&mut self, g: &mut Ghost, a: &QReg, b: &QReg) {
        self.cz_if_bit(a, b, g.bit);
    }

    pub fn ghost_xor_ccz(&mut self, g: &mut Ghost, a: &QReg, b: &QReg, c: &QReg) {
        self.ccz_if_bit(a, b, c, g.bit);
    }

    pub fn close_ghost(&mut self, mut g: Ghost) {
        self.free_bit(g.bit);
        g.consumed = true;
    }

    pub fn x_if_bit(&mut self, q: &QReg, bit: Cbit) {
        self.flush_pending_frees();
        self.b.x_if(QubitId(q.id.into()), BitId(bit.0.into()));
    }

    pub fn z_if_bit(&mut self, q: &QReg, bit: Cbit) {
        self.flush_pending_frees();
        self.b.z_if(QubitId(q.id.into()), BitId(bit.0.into()));
    }

    pub fn cx_if_bit(&mut self, ctrl: &QReg, tgt: &QReg, bit: Cbit) {
        self.flush_pending_frees();
        let mut op = Op::empty();
        op.kind = OperationType::CX;
        op.q_control1 = QubitId(ctrl.id.into());
        op.q_target = QubitId(tgt.id.into());
        op.c_condition = BitId(bit.0.into());
        self.b.push_op(op);
    }

    pub fn cz_if_bit(&mut self, a: &QReg, b: &QReg, bit: Cbit) {
        self.flush_pending_frees();
        self.b.cz_if(
            QubitId(a.id.into()),
            QubitId(b.id.into()),
            BitId(bit.0.into()),
        );
    }

    pub fn ccx_if_bit(&mut self, a: &QReg, b: &QReg, t: &QReg, bit: Cbit) {
        self.flush_pending_frees();
        let mut op = Op::empty();
        op.kind = OperationType::CCX;
        op.q_control2 = QubitId(a.id.into());
        op.q_control1 = QubitId(b.id.into());
        op.q_target = QubitId(t.id.into());
        op.c_condition = BitId(bit.0.into());
        self.b.push_op(op);
    }

    pub fn ccz_if_bit(&mut self, a: &QReg, b: &QReg, c: &QReg, bit: Cbit) {
        self.flush_pending_frees();
        let mut op = Op::empty();
        op.kind = OperationType::CCZ;
        op.q_control2 = QubitId(a.id.into());
        op.q_control1 = QubitId(b.id.into());
        op.q_target = QubitId(c.id.into());
        op.c_condition = BitId(bit.0.into());
        self.b.push_op(op);
    }

    pub fn with_condition<R>(&mut self, bit: Cbit, f: impl FnOnce(&mut Self) -> R) -> R {
        self.flush_pending_frees();
        self.b.push_condition(BitId(bit.0.into()));
        let out = f(self);
        self.flush_pending_frees();
        self.b.pop_condition();
        out
    }

    pub fn with_conditions<R>(&mut self, bits: &[Cbit], f: impl FnOnce(&mut Self) -> R) -> R {
        for &bit in bits {
            self.flush_pending_frees();
            self.b.push_condition(BitId(bit.0.into()));
        }
        let out = f(self);
        for _ in bits {
            self.flush_pending_frees();
            self.b.pop_condition();
        }
        out
    }

    pub fn clear_and(&mut self, t: &QReg, a: &QReg, b: &QReg) {
        self.declare_and_of(t, a, b);
        let bit = self.alloc_bit();
        self.hmr(t, bit);
        self.cz_if_bit(a, b, bit);
        self.free_bit(bit);
    }

    pub fn declare_identity(&mut self, _q: &QReg, _source: &QReg) {}

    pub fn declare_copy_of(&mut self, _q: &QReg, _source: &QReg) {}

    pub fn declare_and_of(&mut self, _q: &QReg, _a: &QReg, _b: &QReg) {}

    pub fn declare_and3_of(&mut self, _q: &QReg, _a: &QReg, _b: &QReg, _c: &QReg) {}

    pub fn declare_xor_of(&mut self, _q: &QReg, _a: &QReg, _b: &QReg) {}

    pub fn declare_xor_of_three(&mut self, _q: &QReg, _a: &QReg, _b: &QReg, _c: &QReg) {}

    pub fn push_section(&mut self, sub: &str) -> String {
        let prev = self.current_section.clone();
        self.set_section(&format!("{prev}/{sub}"));
        prev
    }

    pub fn pop_section(&mut self, prev: &str) {
        self.set_section(prev);
    }

    pub fn set_section(&mut self, s: &str) {
        self.flush_pending_frees();
        self.current_section = s.to_string();
        let leaked: &'static str = Box::leak(s.to_string().into_boxed_str());
        self.b.set_phase(leaked);
    }

    pub fn to_kmx(&self) -> String {
        String::new()
    }
}

impl fmt::Debug for Circuit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Circuit")
            .field("current_section", &self.current_section)
            .field("peak_qubits", &self.b.peak_qubits)
            .finish()
    }
}
