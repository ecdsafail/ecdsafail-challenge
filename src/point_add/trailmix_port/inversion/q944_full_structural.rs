//! Closed production route table for the structural Q944 integration.
//!
//! Nine classes use paired support-zero lanes proved by WMICluster job 71935.
//! The five division classes without such a host use the quotient-witness
//! construction proved by job 71968. No class may silently fall back.

use super::q945_local_hosts::{Q945Host, Q945StateRegister, Q945Substep, Q945_NON_HCLZ_ROWS};

pub const Q944_GATE_HOST_CENSUS_COMMIT: &str =
    "9ee33567e2a4e20176300ede56ad01b3bf86fcab";
pub const Q944_GATE_HOST_CENSUS_TREE: &str =
    "46a78e5c5141fe3b1f169c739b699b18946c6015";
pub const Q944_GATE_HOST_CENSUS_JOB: usize = 71_935;
pub const Q944_QUOTIENT_WITNESS_COMMIT: &str =
    "1817f74f2fae27d622b46151fef45cf68acc0902";
pub const Q944_QUOTIENT_WITNESS_TREE: &str =
    "f5aad72a891e404006a929d853bdc45441aca866";
pub const Q944_QUOTIENT_WITNESS_JOB: usize = 71_968;
pub const Q944_QUOTIENT_WITNESS_BLOB: &str =
    "be5c4e15916eea2121fc120809dead8e0fbfd0a3";
pub const Q944_ORDINARY_CLASSES: usize = 9;
pub const Q944_QUOTIENT_CLASSES: usize = 5;
pub const Q944_ORDINARY_SITES: usize = 36;
pub const Q944_QUOTIENT_SITES: usize = 20;
pub const Q944_TOTAL_SITES: usize = Q944_ORDINARY_SITES + Q944_QUOTIENT_SITES;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Q944FullGateRoute {
    Ordinary { host: Q945Host, peer: Q945Host },
    QuotientWitness,
}

#[must_use]
pub fn q944_full_gate_route(row: usize, substep: Q945Substep) -> Q944FullGateRoute {
    use Q945StateRegister::{A, B, Ca, Cb};
    let pair = |host, host_bit, peer, peer_bit| Q944FullGateRoute::Ordinary {
        host: Q945Host::new(host, host_bit),
        peer: Q945Host::new(peer, peer_bit),
    };
    match (row, substep) {
        (363, Q945Substep::Division) => pair(Ca, 246, Cb, 246),
        (363, Q945Substep::Multiply) => pair(A, 80, B, 80),
        (364, Q945Substep::Division) => pair(Ca, 246, Cb, 246),
        (364, Q945Substep::Multiply) => pair(B, 80, A, 80),
        (374, Q945Substep::Division)
        | (375, Q945Substep::Division)
        | (376, Q945Substep::Division)
        | (379, Q945Substep::Division)
        | (380, Q945Substep::Division) => Q944FullGateRoute::QuotientWitness,
        (374, Q945Substep::Multiply) => pair(B, 73, A, 73),
        (375, Q945Substep::Multiply) => pair(A, 72, B, 72),
        (376, Q945Substep::Multiply) => pair(B, 72, A, 72),
        (379, Q945Substep::Multiply) => pair(B, 71, A, 71),
        (380, Q945Substep::Multiply) => pair(B, 71, A, 71),
        _ => panic!(
            "unclassified Q944 full structural class row={row} substep={}",
            substep.label()
        ),
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Q944FullStaticReport {
    pub rows: usize,
    pub classes: usize,
    pub ordinary_classes: usize,
    pub quotient_classes: usize,
    pub ordinary_sites: usize,
    pub quotient_sites: usize,
    pub total_sites: usize,
}

#[must_use]
pub fn assert_q944_full_static_route() -> Q944FullStaticReport {
    let mut ordinary = 0usize;
    let mut quotient = 0usize;
    for row in Q945_NON_HCLZ_ROWS {
        for substep in Q945Substep::ALL {
            match q944_full_gate_route(row, substep) {
                Q944FullGateRoute::Ordinary { host, peer } => {
                    assert_eq!(host.bit, peer.bit);
                    assert_ne!(host.register, peer.register);
                    ordinary += 1;
                }
                Q944FullGateRoute::QuotientWitness => {
                    assert_eq!(substep, Q945Substep::Division);
                    assert!([374, 375, 376, 379, 380].contains(&row));
                    quotient += 1;
                }
            }
        }
    }
    assert_eq!(ordinary, Q944_ORDINARY_CLASSES);
    assert_eq!(quotient, Q944_QUOTIENT_CLASSES);
    Q944FullStaticReport {
        rows: Q945_NON_HCLZ_ROWS.len(),
        classes: ordinary + quotient,
        ordinary_classes: ordinary,
        quotient_classes: quotient,
        ordinary_sites: 4 * ordinary,
        quotient_sites: 4 * quotient,
        total_sites: 4 * (ordinary + quotient),
    }
}
