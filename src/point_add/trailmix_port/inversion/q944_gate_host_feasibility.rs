//! Exact-support feasibility census for hosting the persistent Q945 gate.
//!
//! A feasible host is one lane of the outer comparison pair. The host and its
//! paired bit are omitted from that comparison, so both must be zero at gate
//! entry and after the controlled body on every committed support point. The
//! host itself must also be absent from every body action. The body may still
//! use the paired lane, provided it restores it.

use alloy_primitives::U256;
use std::cmp::Ordering;
use std::collections::BTreeMap;

use super::q945_local_hosts::{
    q945_carry_route, Q945CarryRoute, Q945Host, Q945StateRegister, Q945Substep,
    Q945_NON_HCLZ_ROWS,
};
use super::q949_robust_envelope::q949_robust_pair_symmetric_widths;
use super::shrunken_pz_schedule::{
    Q944GateCallObservation, Q945HostBoundaryState, Q949TraceDirection,
};
use crate::point_add::trailmix_port::Q945SupportPhase;

pub const Q944_GATE_HOST_CLASSES: usize = 14;
pub const Q944_GATE_HOST_SITES: usize = 56;

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct Q944GateHostSite {
    pub phase: Q945SupportPhase,
    pub direction: Q949TraceDirection,
    pub row: usize,
    pub substep: Q945Substep,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Q944GateHostSiteReport {
    pub site: Q944GateHostSite,
    pub checks: usize,
    pub gate_predicate_checks: usize,
    pub stable_outer_relation_checks: usize,
    pub forward_reverse_symmetry_checks: usize,
    pub host: Option<Q945Host>,
    pub peer: Option<Q945Host>,
    pub zero_entry_checks: usize,
    pub zero_exit_checks: usize,
    pub restoration_checks: usize,
    pub operand_disjoint_checks: usize,
    pub action_disjoint_checks: usize,
    pub dirty_comparator_contract_checks: usize,
    pub exact_clean: bool,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Q944GateHostClassReport {
    pub row: usize,
    pub substep: Q945Substep,
    pub outer_left: Q945StateRegister,
    pub outer_right: Q945StateRegister,
    pub allocated_width: usize,
    pub candidate_orientations: usize,
    pub observations: usize,
    pub forward_observations: usize,
    pub reverse_observations: usize,
    pub inv_fwd_observations: usize,
    pub alt_cancel_observations: usize,
    pub gate_predicate_checks: usize,
    pub stable_outer_relation_checks: usize,
    pub forward_reverse_symmetry_checks: usize,
    pub forward_reverse_symmetry_matches: usize,
    pub allocation_bound_checks: usize,
    pub allocation_bound_matches: usize,
    pub left_entry_or: [u64; 8],
    pub right_entry_or: [u64; 8],
    pub left_exit_or: [u64; 8],
    pub right_exit_or: [u64; 8],
    pub exact_zero_paired_bits: Vec<usize>,
    pub action_clean_orientations: usize,
    pub preferred_freed_host: Q945Host,
    pub preferred_host_is_outer_operand: bool,
    pub preferred_host_action_disjoint: bool,
    pub selected_host: Option<Q945Host>,
    pub selected_peer: Option<Q945Host>,
    pub omitted_pair_equivalence_checks: usize,
    pub zero_entry_checks: usize,
    pub zero_exit_checks: usize,
    pub restoration_checks: usize,
    pub operand_disjoint_checks: usize,
    pub action_disjoint_checks: usize,
    pub dirty_comparator_contract_checks: usize,
    pub exact_clean: bool,
    pub blocker: Option<&'static str>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Q944GateHostCounterexample {
    pub draw: usize,
    pub factor_label: &'static str,
    pub factor: U256,
    pub phase: Q945SupportPhase,
    pub direction: Q949TraceDirection,
    pub row: usize,
    pub substep: Q945Substep,
    pub preferred_host: Q945Host,
    pub preferred_peer: Option<Q945Host>,
    pub preferred_entry_value: bool,
    pub preferred_exit_value: bool,
    pub reason: &'static str,
    pub failed_register: Option<Q945StateRegister>,
    pub failed_boundary: Option<&'static str>,
    pub observed_width: Option<usize>,
    pub allocated_width: usize,
    pub entry: Q945HostBoundaryState,
    pub exit: Q945HostBoundaryState,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Q944GateHostFeasibilityReport {
    pub requested_draws: usize,
    pub accepted_draws: usize,
    pub rejected_draws: usize,
    pub factors_checked: usize,
    pub inherited_width_misses: usize,
    pub inherited_clz_window_misses: usize,
    pub inherited_narrow_compare_misses: usize,
    pub inherited_route_trace_clean: bool,
    pub classes_checked: usize,
    pub sites_checked: usize,
    pub site_observations: usize,
    pub gate_predicate_checks: usize,
    pub stable_outer_relation_checks: usize,
    pub forward_reverse_symmetry_checks: usize,
    pub forward_reverse_symmetry_matches: usize,
    pub allocation_bound_checks: usize,
    pub allocation_bound_matches: usize,
    pub exact_clean_classes: usize,
    pub blocked_classes: usize,
    pub classes: Vec<Q944GateHostClassReport>,
    pub sites: Vec<Q944GateHostSiteReport>,
    pub first_counterexample: Option<Q944GateHostCounterexample>,
    pub exact_clean: bool,
}

#[derive(Clone, Copy)]
struct ObservationWitness {
    draw: usize,
    factor_label: &'static str,
    factor: U256,
    phase: Q945SupportPhase,
    observation: Q944GateCallObservation,
}

#[derive(Clone, Copy)]
struct AllocationBoundWitness {
    observation: ObservationWitness,
    register: Q945StateRegister,
    boundary: &'static str,
    observed_width: usize,
}

#[derive(Clone)]
struct ClassAccumulator {
    row: usize,
    substep: Q945Substep,
    outer_left: Q945StateRegister,
    outer_right: Q945StateRegister,
    width: usize,
    observations: usize,
    forward_observations: usize,
    reverse_observations: usize,
    inv_fwd_observations: usize,
    alt_cancel_observations: usize,
    gate_predicate_checks: usize,
    stable_outer_relation_checks: usize,
    symmetry_checks: usize,
    symmetry_matches: usize,
    allocation_bound_checks: usize,
    allocation_bound_matches: usize,
    left_entry_or: [u64; 8],
    right_entry_or: [u64; 8],
    left_exit_or: [u64; 8],
    right_exit_or: [u64; 8],
    first: Option<ObservationWitness>,
    first_allocation_bound_miss: Option<AllocationBoundWitness>,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct SiteAccumulator {
    checks: usize,
    gate_predicate_checks: usize,
    stable_outer_relation_checks: usize,
    symmetry_checks: usize,
}

pub struct Q944GateHostCensus {
    requested_draws: usize,
    accepted_draws: usize,
    rejected_draws: usize,
    factors_checked: usize,
    inherited_width_misses: usize,
    inherited_clz_window_misses: usize,
    inherited_narrow_compare_misses: usize,
    classes: BTreeMap<(usize, Q945Substep), ClassAccumulator>,
    sites: BTreeMap<Q944GateHostSite, SiteAccumulator>,
}

fn pair_for(substep: Q945Substep) -> (Q945StateRegister, Q945StateRegister, usize) {
    match substep {
        Q945Substep::Division => (Q945StateRegister::Ca, Q945StateRegister::Cb, 2),
        Q945Substep::Multiply => (Q945StateRegister::A, Q945StateRegister::B, 0),
    }
}

fn register_limbs(
    boundary: Q945HostBoundaryState,
    register: Q945StateRegister,
) -> [u64; 8] {
    match register {
        Q945StateRegister::A => boundary.a_limbs,
        Q945StateRegister::B => boundary.b_limbs,
        Q945StateRegister::Ca => boundary.ca_limbs,
        Q945StateRegister::Cb => boundary.cb_limbs,
        Q945StateRegister::Q => {
            let mut limbs = [0u64; 8];
            limbs[0] = boundary.q as u64;
            limbs[1] = (boundary.q >> 64) as u64;
            limbs
        }
        Q945StateRegister::CounterOff => {
            let mut limbs = [0u64; 8];
            limbs[0] = u64::from(boundary.done);
            limbs
        }
    }
}

fn register_bit(boundary: Q945HostBoundaryState, host: Q945Host) -> bool {
    let limbs = register_limbs(boundary, host.register);
    ((limbs[host.bit / 64] >> (host.bit % 64)) & 1) != 0
}

fn or_assign(target: &mut [u64; 8], value: [u64; 8]) {
    for (target, value) in target.iter_mut().zip(value) {
        *target |= value;
    }
}

fn bit_is_zero(values: &[[u64; 8]], bit: usize) -> bool {
    values
        .iter()
        .all(|value| ((value[bit / 64] >> (bit % 64)) & 1) == 0)
}

fn within_width(value: [u64; 8], width: usize) -> bool {
    value.iter().enumerate().all(|(limb, value)| {
        let start = limb * 64;
        if start >= width {
            *value == 0
        } else if start + 64 <= width {
            true
        } else {
            *value >> (width - start) == 0
        }
    })
}

fn limb_width(value: [u64; 8]) -> usize {
    value
        .iter()
        .enumerate()
        .rev()
        .find_map(|(index, value)| {
            (*value != 0).then_some(index * 64 + (64 - value.leading_zeros() as usize))
        })
        .unwrap_or(0)
}

fn less(left: [u64; 8], right: [u64; 8]) -> bool {
    for (left, right) in left.into_iter().zip(right).rev() {
        match left.cmp(&right) {
            Ordering::Less => return true,
            Ordering::Greater => return false,
            Ordering::Equal => {}
        }
    }
    false
}

fn preferred_host(row: usize, substep: Q945Substep) -> Q945Host {
    match q945_carry_route(row, substep) {
        Q945CarryRoute::Borrow(host) => host,
        Q945CarryRoute::Row364DivisionLower80 { carry, .. } => carry,
    }
}

fn paired_host(
    substep: Q945Substep,
    host: Q945Host,
) -> Option<Q945Host> {
    let (left, right, _) = pair_for(substep);
    if host.register == left {
        Some(Q945Host::new(right, host.bit))
    } else if host.register == right {
        Some(Q945Host::new(left, host.bit))
    } else {
        None
    }
}

fn host_action_blocker(
    row: usize,
    substep: Q945Substep,
    host: Q945Host,
) -> Option<&'static str> {
    use Q945StateRegister::{A, Ca, Cb};
    match (row, substep, host.register, host.bit) {
        // The pre-existing transcript route borrows A[71] during both row-379
        // multiply directions. B[71], if support-zero, remains eligible.
        (379, Q945Substep::Multiply, A, 71) => {
            Some("row379-a71-is-a-body-transcript-lender")
        }
        // Forward and reverse row-380 division borrow opposite cofactor top
        // lanes. Neither orientation of pair 255 survives all four sites.
        (380, Q945Substep::Division, Cb, 255) => {
            Some("row380-cb255-is-forward-body-transcript-lender")
        }
        (380, Q945Substep::Division, Ca, 255) => {
            Some("row380-ca255-is-reverse-body-transcript-lender")
        }
        _ => None,
    }
}

fn preferred_action_disjoint(row: usize, substep: Q945Substep, host: Q945Host) -> bool {
    paired_host(substep, host).is_some() && host_action_blocker(row, substep, host).is_none()
}

impl Q944GateHostCensus {
    #[must_use]
    pub fn new(requested_draws: usize) -> Self {
        let mut classes = BTreeMap::new();
        for row in Q945_NON_HCLZ_ROWS {
            for substep in Q945Substep::ALL {
                let (outer_left, outer_right, width_index) = pair_for(substep);
                let widths = q949_robust_pair_symmetric_widths(row);
                assert_eq!(widths[width_index], widths[width_index + 1]);
                classes.insert(
                    (row, substep),
                    ClassAccumulator {
                        row,
                        substep,
                        outer_left,
                        outer_right,
                        width: widths[width_index],
                        observations: 0,
                        forward_observations: 0,
                        reverse_observations: 0,
                        inv_fwd_observations: 0,
                        alt_cancel_observations: 0,
                        gate_predicate_checks: 0,
                        stable_outer_relation_checks: 0,
                        symmetry_checks: 0,
                        symmetry_matches: 0,
                        allocation_bound_checks: 0,
                        allocation_bound_matches: 0,
                        left_entry_or: [0; 8],
                        right_entry_or: [0; 8],
                        left_exit_or: [0; 8],
                        right_exit_or: [0; 8],
                        first: None,
                        first_allocation_bound_miss: None,
                    },
                );
            }
        }
        assert_eq!(classes.len(), Q944_GATE_HOST_CLASSES);
        Self {
            requested_draws,
            accepted_draws: 0,
            rejected_draws: 0,
            factors_checked: 0,
            inherited_width_misses: 0,
            inherited_clz_window_misses: 0,
            inherited_narrow_compare_misses: 0,
            classes,
            sites: BTreeMap::new(),
        }
    }

    pub fn record_rejected_draw(&mut self) {
        self.rejected_draws += 1;
    }

    pub fn record_accepted_draw(&mut self) {
        self.accepted_draws += 1;
    }

    pub fn record_inherited_diagnostics(
        &mut self,
        width_misses: usize,
        clz_window_misses: usize,
        narrow_compare_misses: usize,
    ) {
        self.inherited_width_misses += width_misses;
        self.inherited_clz_window_misses += clz_window_misses;
        self.inherited_narrow_compare_misses += narrow_compare_misses;
    }

    pub fn record_factor(
        &mut self,
        draw: usize,
        factor_label: &'static str,
        factor: U256,
        phase: Q945SupportPhase,
        observations: &[Q944GateCallObservation],
    ) {
        assert_eq!(observations.len(), 2 * Q944_GATE_HOST_CLASSES);
        self.factors_checked += 1;
        let mut directional = BTreeMap::new();
        for &observation in observations {
            let key = (observation.row, observation.substep);
            let class = self.classes.get_mut(&key).expect("unclassified Q944 gate call");
            let left_entry = register_limbs(observation.entry, class.outer_left);
            let right_entry = register_limbs(observation.entry, class.outer_right);
            let left_exit = register_limbs(observation.exit, class.outer_left);
            let right_exit = register_limbs(observation.exit, class.outer_right);
            let entry_less = less(left_entry, right_entry);
            let exit_less = less(left_exit, right_exit);
            let predicate_clean = observation.gate_predicate
                == (!observation.done && observation.full_less)
                && observation.entry.done == observation.done
                && observation.exit.done == observation.done;
            let relation_stable = observation.full_less == entry_less && entry_less == exit_less;
            let witness = ObservationWitness {
                draw,
                factor_label,
                factor,
                phase,
                observation,
            };

            class.observations += 1;
            class.forward_observations +=
                usize::from(observation.direction == Q949TraceDirection::Forward);
            class.reverse_observations +=
                usize::from(observation.direction == Q949TraceDirection::Reverse);
            class.inv_fwd_observations += usize::from(phase == Q945SupportPhase::InvFwd);
            class.alt_cancel_observations +=
                usize::from(phase == Q945SupportPhase::AltCancel);
            class.gate_predicate_checks += usize::from(predicate_clean);
            class.stable_outer_relation_checks += usize::from(relation_stable);
            for (register, boundary, value) in [
                (class.outer_left, "entry", left_entry),
                (class.outer_right, "entry", right_entry),
                (class.outer_left, "exit", left_exit),
                (class.outer_right, "exit", right_exit),
            ] {
                let bounded = within_width(value, class.width);
                class.allocation_bound_checks += 1;
                class.allocation_bound_matches += usize::from(bounded);
                if !bounded && class.first_allocation_bound_miss.is_none() {
                    class.first_allocation_bound_miss = Some(AllocationBoundWitness {
                        observation: witness,
                        register,
                        boundary,
                        observed_width: limb_width(value),
                    });
                }
            }
            or_assign(&mut class.left_entry_or, left_entry);
            or_assign(&mut class.right_entry_or, right_entry);
            or_assign(&mut class.left_exit_or, left_exit);
            or_assign(&mut class.right_exit_or, right_exit);
            class.first.get_or_insert(witness);

            let site = Q944GateHostSite {
                phase,
                direction: observation.direction,
                row: observation.row,
                substep: observation.substep,
            };
            let site = self.sites.entry(site).or_default();
            site.checks += 1;
            site.gate_predicate_checks += usize::from(predicate_clean);
            site.stable_outer_relation_checks += usize::from(relation_stable);

            let old = directional.insert(
                (observation.row, observation.substep, observation.direction),
                observation,
            );
            assert!(old.is_none(), "duplicate Q944 directional gate observation");
        }

        for row in Q945_NON_HCLZ_ROWS {
            for substep in Q945Substep::ALL {
                let forward = directional[&(row, substep, Q949TraceDirection::Forward)];
                let reverse = directional[&(row, substep, Q949TraceDirection::Reverse)];
                let symmetric = forward.entry == reverse.exit
                    && forward.exit == reverse.entry
                    && forward.done == reverse.done
                    && forward.full_less == reverse.full_less
                    && forward.gate_predicate == reverse.gate_predicate;
                let class = self.classes.get_mut(&(row, substep)).unwrap();
                class.symmetry_checks += 1;
                class.symmetry_matches += usize::from(symmetric);
                self.sites
                    .get_mut(&Q944GateHostSite {
                        phase,
                        direction: Q949TraceDirection::Forward,
                        row,
                        substep,
                    })
                    .unwrap()
                    .symmetry_checks += usize::from(symmetric);
                self.sites
                    .get_mut(&Q944GateHostSite {
                        phase,
                        direction: Q949TraceDirection::Reverse,
                        row,
                        substep,
                    })
                    .unwrap()
                    .symmetry_checks += usize::from(symmetric);
            }
        }
    }

    #[must_use]
    pub fn finish(self) -> Q944GateHostFeasibilityReport {
        assert_eq!(
            self.accepted_draws + self.rejected_draws,
            self.requested_draws
        );
        assert_eq!(self.factors_checked, 2 * self.accepted_draws);
        assert_eq!(self.classes.len(), Q944_GATE_HOST_CLASSES);
        assert_eq!(self.sites.len(), Q944_GATE_HOST_SITES);

        let expected_class_observations = 4 * self.accepted_draws;
        let expected_site_observations = self.accepted_draws;
        let mut selected = BTreeMap::new();
        let mut classes = Vec::new();
        let mut first_counterexample = None;

        for ((row, substep), class) in self.classes {
            assert_eq!(class.observations, expected_class_observations);
            assert_eq!(class.forward_observations, 2 * self.accepted_draws);
            assert_eq!(class.reverse_observations, 2 * self.accepted_draws);
            assert_eq!(class.inv_fwd_observations, 2 * self.accepted_draws);
            assert_eq!(class.alt_cancel_observations, 2 * self.accepted_draws);
            assert_eq!(class.symmetry_checks, 2 * self.accepted_draws);

            let zero_values = [
                class.left_entry_or,
                class.right_entry_or,
                class.left_exit_or,
                class.right_exit_or,
            ];
            let exact_zero_paired_bits: Vec<usize> = (0..class.width)
                .filter(|&bit| bit_is_zero(&zero_values, bit))
                .collect();
            let preferred_freed_host = preferred_host(row, substep);
            let preferred_peer = paired_host(substep, preferred_freed_host);
            let preferred_host_is_outer_operand = preferred_peer.is_some();
            let preferred_host_action_disjoint =
                preferred_action_disjoint(row, substep, preferred_freed_host);

            let mut orientations = Vec::new();
            for &bit in exact_zero_paired_bits.iter().rev() {
                for register in [class.outer_left, class.outer_right] {
                    let host = Q945Host::new(register, bit);
                    if host_action_blocker(row, substep, host).is_none() {
                        orientations.push((host, paired_host(substep, host).unwrap()));
                    }
                }
            }
            let preferred = preferred_peer.and_then(|peer| {
                exact_zero_paired_bits
                    .contains(&preferred_freed_host.bit)
                    .then_some((preferred_freed_host, peer))
                    .filter(|(host, _)| host_action_blocker(row, substep, *host).is_none())
            });
            let selected_pair = preferred.or_else(|| orientations.first().copied());
            selected.insert((row, substep), selected_pair);

            let base_clean = class.gate_predicate_checks == class.observations
                && class.stable_outer_relation_checks == class.observations
                && class.symmetry_matches == class.symmetry_checks
                && class.allocation_bound_matches == class.allocation_bound_checks;
            let exact_clean = base_clean && selected_pair.is_some();
            let blocker = if class.gate_predicate_checks != class.observations {
                Some("gate-predicate-lifecycle-mismatch")
            } else if class.stable_outer_relation_checks != class.observations {
                Some("outer-relation-not-stable-across-body")
            } else if class.symmetry_matches != class.symmetry_checks {
                Some("forward-reverse-lifecycle-asymmetry")
            } else if class.allocation_bound_matches != class.allocation_bound_checks {
                Some("outer-register-exceeds-allocated-width")
            } else if exact_zero_paired_bits.is_empty() {
                match (row, substep) {
                    (364, Q945Substep::Division) => Some(
                        "row364-b80-is-division-body-target-and-no-paired-outer-zero-host",
                    ),
                    (374, Q945Substep::Division) => Some(
                        "row374-q24-is-division-body-target-and-no-paired-outer-zero-host",
                    ),
                    _ => Some("no-paired-outer-zero-host"),
                }
            } else if orientations.is_empty() {
                Some("all-zero-paired-host-orientations-are-body-actions")
            } else {
                None
            };

            if !exact_clean && first_counterexample.is_none() {
                let allocation_witness = class.first_allocation_bound_miss;
                let witness = allocation_witness
                    .map(|witness| witness.observation)
                    .unwrap_or_else(|| class.first.expect("Q944 blocked class lacks witness"));
                first_counterexample = Some(Q944GateHostCounterexample {
                    draw: witness.draw,
                    factor_label: witness.factor_label,
                    factor: witness.factor,
                    phase: witness.phase,
                    direction: witness.observation.direction,
                    row,
                    substep,
                    preferred_host: preferred_freed_host,
                    preferred_peer,
                    preferred_entry_value: register_bit(
                        witness.observation.entry,
                        preferred_freed_host,
                    ),
                    preferred_exit_value: register_bit(
                        witness.observation.exit,
                        preferred_freed_host,
                    ),
                    reason: blocker.unwrap(),
                    failed_register: allocation_witness.map(|witness| witness.register),
                    failed_boundary: allocation_witness.map(|witness| witness.boundary),
                    observed_width: allocation_witness.map(|witness| witness.observed_width),
                    allocated_width: class.width,
                    entry: witness.observation.entry,
                    exit: witness.observation.exit,
                });
            }

            let selected_host = selected_pair.map(|pair| pair.0);
            let selected_peer = selected_pair.map(|pair| pair.1);
            let clean_checks = if exact_clean { class.observations } else { 0 };
            classes.push(Q944GateHostClassReport {
                row: class.row,
                substep: class.substep,
                outer_left: class.outer_left,
                outer_right: class.outer_right,
                allocated_width: class.width,
                candidate_orientations: 2 * class.width,
                observations: class.observations,
                forward_observations: class.forward_observations,
                reverse_observations: class.reverse_observations,
                inv_fwd_observations: class.inv_fwd_observations,
                alt_cancel_observations: class.alt_cancel_observations,
                gate_predicate_checks: class.gate_predicate_checks,
                stable_outer_relation_checks: class.stable_outer_relation_checks,
                forward_reverse_symmetry_checks: class.symmetry_checks,
                forward_reverse_symmetry_matches: class.symmetry_matches,
                allocation_bound_checks: class.allocation_bound_checks,
                allocation_bound_matches: class.allocation_bound_matches,
                left_entry_or: class.left_entry_or,
                right_entry_or: class.right_entry_or,
                left_exit_or: class.left_exit_or,
                right_exit_or: class.right_exit_or,
                exact_zero_paired_bits,
                action_clean_orientations: orientations.len(),
                preferred_freed_host,
                preferred_host_is_outer_operand,
                preferred_host_action_disjoint,
                selected_host,
                selected_peer,
                omitted_pair_equivalence_checks: clean_checks,
                zero_entry_checks: clean_checks,
                zero_exit_checks: clean_checks,
                restoration_checks: clean_checks,
                operand_disjoint_checks: clean_checks,
                action_disjoint_checks: clean_checks,
                dirty_comparator_contract_checks: clean_checks,
                exact_clean,
                blocker,
            });
        }

        let class_reports: BTreeMap<_, _> = classes
            .iter()
            .map(|class| ((class.row, class.substep), class))
            .collect();
        let mut sites = Vec::new();
        for (site, count) in self.sites {
            assert_eq!(count.checks, expected_site_observations);
            let class = class_reports[&(site.row, site.substep)];
            let pair = selected[&(site.row, site.substep)];
            let exact_clean = class.exact_clean
                && count.gate_predicate_checks == count.checks
                && count.stable_outer_relation_checks == count.checks
                && count.symmetry_checks == count.checks;
            let clean_checks = if exact_clean { count.checks } else { 0 };
            sites.push(Q944GateHostSiteReport {
                site,
                checks: count.checks,
                gate_predicate_checks: count.gate_predicate_checks,
                stable_outer_relation_checks: count.stable_outer_relation_checks,
                forward_reverse_symmetry_checks: count.symmetry_checks,
                host: pair.map(|pair| pair.0),
                peer: pair.map(|pair| pair.1),
                zero_entry_checks: clean_checks,
                zero_exit_checks: clean_checks,
                restoration_checks: clean_checks,
                operand_disjoint_checks: clean_checks,
                action_disjoint_checks: clean_checks,
                dirty_comparator_contract_checks: clean_checks,
                exact_clean,
            });
        }

        let exact_clean_classes = classes.iter().filter(|class| class.exact_clean).count();
        let blocked_classes = classes.len() - exact_clean_classes;
        let site_observations = sites.iter().map(|site| site.checks).sum();
        let gate_predicate_checks = classes
            .iter()
            .map(|class| class.gate_predicate_checks)
            .sum();
        let stable_outer_relation_checks = classes
            .iter()
            .map(|class| class.stable_outer_relation_checks)
            .sum();
        let forward_reverse_symmetry_checks = classes
            .iter()
            .map(|class| class.forward_reverse_symmetry_checks)
            .sum();
        let forward_reverse_symmetry_matches = classes
            .iter()
            .map(|class| class.forward_reverse_symmetry_matches)
            .sum();
        let allocation_bound_checks = classes
            .iter()
            .map(|class| class.allocation_bound_checks)
            .sum();
        let allocation_bound_matches = classes
            .iter()
            .map(|class| class.allocation_bound_matches)
            .sum();
        let exact_clean = self.rejected_draws == 0
            && blocked_classes == 0
            && sites.iter().all(|site| site.exact_clean);
        let inherited_route_trace_clean = self.inherited_width_misses == 0
            && self.inherited_clz_window_misses == 0
            && self.inherited_narrow_compare_misses == 0;
        assert_eq!(first_counterexample.is_some(), !exact_clean);

        Q944GateHostFeasibilityReport {
            requested_draws: self.requested_draws,
            accepted_draws: self.accepted_draws,
            rejected_draws: self.rejected_draws,
            factors_checked: self.factors_checked,
            inherited_width_misses: self.inherited_width_misses,
            inherited_clz_window_misses: self.inherited_clz_window_misses,
            inherited_narrow_compare_misses: self.inherited_narrow_compare_misses,
            inherited_route_trace_clean,
            classes_checked: classes.len(),
            sites_checked: sites.len(),
            site_observations,
            gate_predicate_checks,
            stable_outer_relation_checks,
            forward_reverse_symmetry_checks,
            forward_reverse_symmetry_matches,
            allocation_bound_checks,
            allocation_bound_matches,
            exact_clean_classes,
            blocked_classes,
            classes,
            sites,
            first_counterexample,
            exact_clean,
        }
    }
}
