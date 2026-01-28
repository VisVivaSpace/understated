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

/// Select an interpolation window for unequally-spaced data.
///
/// Given a list of epochs (accessed via `epoch_at`), the evaluation epoch,
/// the polynomial degree, and total record count, returns `(start, end)` indices
/// for the interpolation window. Implements the CSPICE window selection algorithm.
#[allow(clippy::if_same_then_else)]
pub(crate) fn select_window_unequal(
    n: usize,
    degree: usize,
    epoch: f64,
    epoch_at: impl Fn(usize) -> f64,
) -> Result<(usize, usize)> {
    if n == 0 {
        return Err(Error::EpochOutOfRange {
            epoch,
            start: 0.0,
            end: 0.0,
        });
    }

    let start_epoch = epoch_at(0);
    let end_epoch = epoch_at(n - 1);

    if epoch < start_epoch || epoch > end_epoch {
        return Err(Error::EpochOutOfRange {
            epoch,
            start: start_epoch,
            end: end_epoch,
        });
    }

    // Binary search for the bracketing interval.
    let mut lower = 0;
    let mut upper = n - 1;
    while upper - lower > 1 {
        let mid = (lower + upper) / 2;
        if epoch_at(mid) <= epoch {
            lower = mid;
        } else {
            upper = mid;
        }
    }
    let high = lower + 1;

    let wndsiz = degree + 1;

    let first = if wndsiz % 2 == 1 {
        let near = if lower == 0 {
            lower
        } else if high >= n {
            lower
        } else if (epoch - epoch_at(lower)).abs() <= (epoch_at(high) - epoch).abs() {
            lower
        } else {
            high
        };
        let half = degree / 2;
        if near < half {
            0
        } else if near > n - 1 - (degree - half) {
            n - wndsiz
        } else {
            near - half
        }
    } else {
        let half = degree / 2;
        if lower < half {
            0
        } else if lower > n - 1 - (degree - half) {
            n - wndsiz
        } else {
            lower - half
        }
    };

    Ok((first, first + degree + 1))
}

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
