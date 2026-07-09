//! Round-level routines: the dialog-GCD inversion subsystem (raw and
//! compressed-sidecar variants, the per-step lever readers, and the fused
//! square+xtail helper) plus the top-level `emit_dialog_gcd_raw_pa` driver.
use super::*;

mod dialog;

pub(crate) use dialog::*;
