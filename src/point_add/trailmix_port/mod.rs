// HARD RULE: no `unsafe` Rust except for explicitly-approved C FFI.
// Every existing `unsafe` block must be either (a) part of an approved
// FFI call site with a per-site `#[allow(unsafe_code)]` and the call
// site recorded in feedback_no_unsafe.md, or (b) rewritten in safe
// Rust. New unsafe is forbidden without explicit user approval.
#![deny(unsafe_code)]

pub mod circuit;

pub mod arith;
pub mod ec;
pub mod inversion;
pub mod tracker;

// Backwards-compatible re-exports: callers still use
// `trailmix::poc_arith`, `crate::mbu_primitives`, etc.
pub use self::arith::{khattar_gidney, mod_arith, rfold_mbu, schrottenloher};
pub use self::ec::point_add;
pub use self::tracker::{debugger, phase_lattice};
