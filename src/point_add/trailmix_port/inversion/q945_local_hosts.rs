//! Closed structural classification for the Q945 local-host route.
//!
//! This module binds the proposed lenders to the sealed Q949 allocation
//! schedule. It does not certify reachable-state zero claims. Q945 acceptance
//! remains blocked until the Q946 hardening prerequisites listed below are
//! integrated on top of this structural implementation.

use super::q949_robust_envelope::q949_robust_pair_symmetric_widths;

pub const Q945_TARGET: usize = 945;
pub const Q945_BASE_COMMIT: &str = "44649ea67d269d4457567a5525f716142886cff7";
pub const Q945_CONTEXTS_PER_CLASS: usize = 4;
pub const Q945_HCLZ_ROWS: [usize; 14] = [
    292, 293, 294, 301, 302, 303, 304, 336, 337, 338, 343, 344, 349, 385,
];
pub const Q945_NON_HCLZ_ROWS: [usize; 7] = [363, 364, 374, 375, 376, 379, 380];

pub const Q945_REQUIRED_Q946_INTEGRATIONS: [&str; 4] = [
    "q946-specific-fresh-certificate",
    "q946-exact-544-site-census",
    "q946-narrowed-comparison-oracle",
    "q946-composed-alias-proof",
];

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum Q945Substep {
    Division,
    Multiply,
}

impl Q945Substep {
    pub const ALL: [Self; 2] = [Self::Division, Self::Multiply];

    pub const fn label(self) -> &'static str {
        match self {
            Self::Division => "division",
            Self::Multiply => "multiply",
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum Q945HclzForm {
    Update,
    Parity,
}

impl Q945HclzForm {
    pub const ALL: [Self; 2] = [Self::Update, Self::Parity];

    pub const fn label(self) -> &'static str {
        match self {
            Self::Update => "update",
            Self::Parity => "parity",
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum Q945StateRegister {
    A,
    B,
    Ca,
    Cb,
    Q,
    CounterOff,
}

impl Q945StateRegister {
    pub const fn label(self) -> &'static str {
        match self {
            Self::A => "A",
            Self::B => "B",
            Self::Ca => "ca",
            Self::Cb => "cb",
            Self::Q => "q",
            Self::CounterOff => "counter[0]/off",
        }
    }

    const fn allocation_index(self) -> Option<usize> {
        match self {
            Self::A => Some(0),
            Self::B => Some(1),
            Self::Ca => Some(2),
            Self::Cb => Some(3),
            Self::Q => Some(4),
            Self::CounterOff => None,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct Q945Host {
    pub register: Q945StateRegister,
    pub bit: usize,
}

impl Q945Host {
    pub const fn new(register: Q945StateRegister, bit: usize) -> Self {
        Self { register, bit }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Q945HclzRoute {
    Borrow(Q945Host),
    Direct,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Q945CarryRoute {
    Borrow(Q945Host),
    Row364DivisionLower80 {
        carry: Q945Host,
        not_gate: Q945Host,
    },
}

pub const Q945_DIRECT_HCLZ_CLASSES: [(usize, Q945Substep, Q945HclzForm); 4] = [
    (294, Q945Substep::Division, Q945HclzForm::Parity),
    (337, Q945Substep::Multiply, Q945HclzForm::Parity),
    (343, Q945Substep::Multiply, Q945HclzForm::Parity),
    (344, Q945Substep::Multiply, Q945HclzForm::Parity),
];

pub const Q945_DIVISION_PARITY_HOSTS: [(usize, Q945Host); 13] = [
    (292, Q945Host::new(Q945StateRegister::Ca, 207)),
    (293, Q945Host::new(Q945StateRegister::Ca, 207)),
    (301, Q945Host::new(Q945StateRegister::Cb, 211)),
    (302, Q945Host::new(Q945StateRegister::Ca, 211)),
    (303, Q945Host::new(Q945StateRegister::Cb, 212)),
    (304, Q945Host::new(Q945StateRegister::Ca, 212)),
    (336, Q945Host::new(Q945StateRegister::Ca, 231)),
    (337, Q945Host::new(Q945StateRegister::Cb, 232)),
    (338, Q945Host::new(Q945StateRegister::Ca, 232)),
    (343, Q945Host::new(Q945StateRegister::Cb, 235)),
    (344, Q945Host::new(Q945StateRegister::Ca, 235)),
    (349, Q945Host::new(Q945StateRegister::Ca, 238)),
    (385, Q945Host::new(Q945StateRegister::Cb, 255)),
];

pub const Q945_MULTIPLY_PARITY_HOSTS: [(usize, Q945Host); 11] = [
    (292, Q945Host::new(Q945StateRegister::Q, 22)),
    (293, Q945Host::new(Q945StateRegister::A, 117)),
    (294, Q945Host::new(Q945StateRegister::B, 117)),
    (301, Q945Host::new(Q945StateRegister::A, 113)),
    (302, Q945Host::new(Q945StateRegister::B, 113)),
    (303, Q945Host::new(Q945StateRegister::A, 112)),
    (304, Q945Host::new(Q945StateRegister::B, 112)),
    (336, Q945Host::new(Q945StateRegister::B, 94)),
    (338, Q945Host::new(Q945StateRegister::B, 93)),
    (349, Q945Host::new(Q945StateRegister::A, 87)),
    (385, Q945Host::new(Q945StateRegister::B, 68)),
];

fn table_host(table: &[(usize, Q945Host)], row: usize) -> Option<Q945Host> {
    table
        .iter()
        .find_map(|(candidate, host)| (*candidate == row).then_some(*host))
}

pub fn q945_hclz_route(
    row: usize,
    substep: Q945Substep,
    form: Q945HclzForm,
) -> Q945HclzRoute {
    assert!(
        Q945_HCLZ_ROWS.contains(&row),
        "unclassified Q945 HCLZ row {row}"
    );
    if form == Q945HclzForm::Update {
        let host = if row == 385 {
            match substep {
                Q945Substep::Division => Q945Host::new(Q945StateRegister::Cb, 255),
                Q945Substep::Multiply => Q945Host::new(Q945StateRegister::B, 68),
            }
        } else {
            assert!(row <= 349, "Q945 off loan escaped the preterminal rows");
            Q945Host::new(Q945StateRegister::CounterOff, 0)
        };
        return Q945HclzRoute::Borrow(host);
    }

    let host = match substep {
        Q945Substep::Division => table_host(&Q945_DIVISION_PARITY_HOSTS, row),
        Q945Substep::Multiply => table_host(&Q945_MULTIPLY_PARITY_HOSTS, row),
    };
    match host {
        Some(host) => Q945HclzRoute::Borrow(host),
        None if Q945_DIRECT_HCLZ_CLASSES.contains(&(row, substep, form)) => {
            Q945HclzRoute::Direct
        }
        None => panic!(
            "unclassified Q945 HCLZ class row={row} substep={} form={}",
            substep.label(),
            form.label()
        ),
    }
}

pub fn q945_carry_route(row: usize, substep: Q945Substep) -> Q945CarryRoute {
    assert!(
        Q945_NON_HCLZ_ROWS.contains(&row),
        "unclassified Q945 borrowed-carry row {row}"
    );
    use Q945StateRegister::{A, B, Ca, Cb, Q};
    match (row, substep) {
        (363, Q945Substep::Multiply) => Q945CarryRoute::Borrow(Q945Host::new(A, 80)),
        (363, Q945Substep::Division) => Q945CarryRoute::Borrow(Q945Host::new(Cb, 247)),
        (364, Q945Substep::Multiply) => Q945CarryRoute::Borrow(Q945Host::new(B, 80)),
        (364, Q945Substep::Division) => Q945CarryRoute::Row364DivisionLower80 {
            carry: Q945Host::new(B, 80),
            not_gate: Q945Host::new(A, 80),
        },
        (374, Q945Substep::Multiply) => Q945CarryRoute::Borrow(Q945Host::new(B, 73)),
        (374, Q945Substep::Division) => Q945CarryRoute::Borrow(Q945Host::new(Q, 24)),
        (375, Q945Substep::Multiply) => Q945CarryRoute::Borrow(Q945Host::new(A, 72)),
        (375, Q945Substep::Division) => Q945CarryRoute::Borrow(Q945Host::new(Cb, 254)),
        (376, Q945Substep::Multiply) => Q945CarryRoute::Borrow(Q945Host::new(B, 72)),
        (376, Q945Substep::Division) => Q945CarryRoute::Borrow(Q945Host::new(Ca, 254)),
        (379, Q945Substep::Multiply) => Q945CarryRoute::Borrow(Q945Host::new(A, 71)),
        (379, Q945Substep::Division) => Q945CarryRoute::Borrow(Q945Host::new(Cb, 255)),
        (380, Q945Substep::Multiply) => Q945CarryRoute::Borrow(Q945Host::new(B, 71)),
        (380, Q945Substep::Division) => Q945CarryRoute::Borrow(Q945Host::new(Cb, 255)),
        _ => panic!(
            "unclassified Q945 borrowed-carry class row={row} substep={}",
            substep.label()
        ),
    }
}

fn assert_host_allocated(row: usize, host: Q945Host) {
    if let Some(index) = host.register.allocation_index() {
        let widths = q949_robust_pair_symmetric_widths(row);
        assert!(
            host.bit < widths[index],
            "Q945 host {}[{}] is outside row {row} allocation {}",
            host.register.label(),
            host.bit,
            widths[index]
        );
    } else {
        assert_eq!(host, Q945Host::new(Q945StateRegister::CounterOff, 0));
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Q945StaticHostReport {
    pub hclz_classes: usize,
    pub borrowed_hclz_classes: usize,
    pub direct_hclz_classes: usize,
    pub borrowed_hclz_sites: usize,
    pub direct_hclz_sites: usize,
    pub hclz_events: usize,
    pub borrowed_carry_classes: usize,
    pub borrowed_carry_calls: usize,
    pub non_hclz_events: usize,
    pub acceptance_prerequisites: usize,
}

pub fn assert_q945_static_host_table() -> Q945StaticHostReport {
    assert_eq!(Q945_BASE_COMMIT.len(), 40);
    assert!(Q945_HCLZ_ROWS.windows(2).all(|pair| pair[0] < pair[1]));
    assert!(Q945_NON_HCLZ_ROWS.windows(2).all(|pair| pair[0] < pair[1]));
    assert!(Q945_HCLZ_ROWS
        .iter()
        .all(|row| !Q945_NON_HCLZ_ROWS.contains(row)));

    let mut borrowed_hclz_classes = 0usize;
    let mut direct = Vec::new();
    for row in Q945_HCLZ_ROWS {
        for substep in Q945Substep::ALL {
            for form in Q945HclzForm::ALL {
                match q945_hclz_route(row, substep, form) {
                    Q945HclzRoute::Borrow(host) => {
                        assert_host_allocated(row, host);
                        borrowed_hclz_classes += 1;
                    }
                    Q945HclzRoute::Direct => direct.push((row, substep, form)),
                }
            }
        }
    }
    assert_eq!(borrowed_hclz_classes, 52);
    assert_eq!(direct, Q945_DIRECT_HCLZ_CLASSES);

    let mut borrowed_carry_classes = 0usize;
    for row in Q945_NON_HCLZ_ROWS {
        for substep in Q945Substep::ALL {
            match q945_carry_route(row, substep) {
                Q945CarryRoute::Borrow(host) => assert_host_allocated(row, host),
                Q945CarryRoute::Row364DivisionLower80 { carry, not_gate } => {
                    assert_eq!((row, substep), (364, Q945Substep::Division));
                    assert_eq!(carry, Q945Host::new(Q945StateRegister::B, 80));
                    assert_eq!(not_gate, Q945Host::new(Q945StateRegister::A, 80));
                    assert_host_allocated(row, carry);
                    assert_host_allocated(row, not_gate);
                }
            }
            borrowed_carry_classes += 1;
        }
    }
    assert_eq!(borrowed_carry_classes, 14);

    let hclz_classes = Q945_HCLZ_ROWS.len() * Q945Substep::ALL.len() * Q945HclzForm::ALL.len();
    let borrowed_hclz_sites = borrowed_hclz_classes * Q945_CONTEXTS_PER_CLASS;
    let direct_hclz_sites = direct.len() * Q945_CONTEXTS_PER_CLASS;
    let borrowed_carry_calls = borrowed_carry_classes * Q945_CONTEXTS_PER_CLASS;
    let non_hclz_events = borrowed_carry_calls * 2;
    assert_eq!(hclz_classes, 56);
    assert_eq!(borrowed_hclz_sites, 208);
    assert_eq!(direct_hclz_sites, 16);
    assert_eq!(borrowed_hclz_sites + direct_hclz_sites, 224);
    assert_eq!(borrowed_carry_calls, 56);
    assert_eq!(non_hclz_events, 112);

    Q945StaticHostReport {
        hclz_classes,
        borrowed_hclz_classes,
        direct_hclz_classes: direct.len(),
        borrowed_hclz_sites,
        direct_hclz_sites,
        hclz_events: borrowed_hclz_sites + direct_hclz_sites,
        borrowed_carry_classes,
        borrowed_carry_calls,
        non_hclz_events,
        acceptance_prerequisites: Q945_REQUIRED_Q946_INTEGRATIONS.len(),
    }
}
