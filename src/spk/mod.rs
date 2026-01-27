//! SPK interpolation algorithms.
//!
//! Dispatches to the correct algorithm based on SPK segment type.

pub mod chain;
pub mod chebyshev;
pub mod hermite;
pub mod lagrange;
pub mod twobody;

use crate::error::{Error, Result};
use crate::state::State;
use muad_dib::kernel::spk_types::SpkData;

/// Evaluate an SPK segment's data at the given epoch.
///
/// Dispatches to the correct interpolation algorithm based on the data type.
pub fn evaluate_spk(data: &SpkData, epoch: f64) -> Result<State> {
    match data {
        SpkData::Type2(d) => chebyshev::evaluate_type2(d, epoch),
        SpkData::Type3(d) => chebyshev::evaluate_type3(d, epoch),
        SpkData::Type5(d) => twobody::evaluate_type5(d, epoch),
        SpkData::Type8(d) => lagrange::evaluate_type8(d, epoch),
        SpkData::Type9(d) => lagrange::evaluate_type9(d, epoch),
        SpkData::Type13(d) => hermite::evaluate_type13(d, epoch),
        SpkData::Raw { spk_type, .. } => Err(Error::UnsupportedSpkType {
            spk_type: *spk_type,
        }),
    }
}
