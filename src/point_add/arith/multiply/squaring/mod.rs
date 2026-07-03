//! Symmetric squaring: schoolbook symmetric square (+ low-peak and
//! self-hosted variants) and the windowed square-row machinery.
use super::*;

mod config;
mod row;
mod schoolbook;
mod selfhosted;

pub(crate) use config::*;
pub(crate) use row::*;
#[allow(unused_imports)] // schoolbook items are dead-code reference impls
pub(crate) use schoolbook::*;
pub(crate) use selfhosted::*;
