//! CK interpolation algorithms.
//!
//! Dispatches to the correct algorithm based on CK segment type.

pub mod evaluate;
pub mod slerp;

use crate::error::{Error, Result};
use crate::pointing::Pointing;
use muad_dib::kernel::ck_types::CkData;

/// Evaluate a CK segment's data at the given SCLK ticks.
pub fn evaluate_ck(data: &CkData, sclk: f64) -> Result<Pointing> {
    match data {
        CkData::Type1(d) => evaluate::evaluate_type1(d, sclk),
        CkData::Type3(d) => evaluate::evaluate_type3(d, sclk),
        CkData::Raw { ck_type, .. } => Err(Error::UnsupportedCkType {
            ck_type: *ck_type,
        }),
    }
}
