//! Leap second data and TDB/UTC conversion.
//!
//! Wraps muad-dib's [`LeapSecondData`](muad_dib::spice::lsk::LeapSecondData)
//! to provide UTC↔TDB conversion using the CSPICE-accurate `deltet()` algorithm.
//!
//! # NAIF Time System Relationships
//!
//! ```text
//! UTC -> TAI -> TT -> TDB
//!      +dAT  +32.184s  +periodic
//! ```

use crate::error::Result;
use crate::time::{format_calendar, format_iso8601, TimeFormat};
use crate::types::EpochTDB;

/// Leap second data extracted from an LSK file.
///
/// Wraps muad-dib's `LeapSecondData` and delegates all conversion methods.
#[derive(Debug, Clone)]
pub struct LeapSecondData {
    inner: muad_dib::spice::lsk::LeapSecondData,
}

impl LeapSecondData {
    /// Get the TAI-UTC offset (dAT) for a given TDB epoch.
    pub fn delta_at(&self, tdb: f64) -> f64 {
        self.inner.delta_at(tdb)
    }

    /// Compute ET - UTC (delta ET) for a given epoch.
    /// Matches CSPICE's `deltet_` exactly.
    pub fn deltet(&self, epoch: f64, epoch_type: muad_dib::spice::lsk::EpochType) -> f64 {
        self.inner.deltet(epoch, epoch_type)
    }

    /// Convert UTC seconds past J2000 to TDB.
    pub fn utc_to_tdb_seconds(&self, utc: f64) -> f64 {
        self.inner.utc_to_tdb_seconds(utc)
    }

    /// Convert TDB to UTC seconds past J2000.
    pub fn tdb_to_utc_seconds(&self, tdb: f64) -> f64 {
        self.inner.tdb_to_utc_seconds(tdb)
    }

    /// Convert a UTC time string to TDB epoch.
    pub fn utc_to_tdb(&self, utc_str: &str) -> Result<EpochTDB> {
        let utc_parsed = EpochTDB::parse(utc_str)?;
        let tdb = self.utc_to_tdb_seconds(utc_parsed.0);
        Ok(EpochTDB(tdb))
    }

    /// Convert a TDB epoch to a formatted UTC string.
    pub fn tdb_to_utc(&self, tdb: EpochTDB, format: TimeFormat) -> Result<String> {
        let utc_seconds = self.tdb_to_utc_seconds(tdb.0);

        let formatted = match format {
            TimeFormat::Iso8601 => format_iso8601(utc_seconds),
            TimeFormat::Calendar => format_calendar(utc_seconds),
            TimeFormat::JulianDate => {
                let jd = utc_seconds / 86400.0 + 2451545.0;
                format!("JD {:.6}", jd)
            }
        };

        Ok(formatted)
    }

    /// Access the DELTET/DELTA_T_A constant.
    pub fn delta_t_a(&self) -> f64 {
        self.inner.delta_t_a
    }

    /// Access the DELTET/K constant.
    pub fn k(&self) -> f64 {
        self.inner.k
    }

    /// Access the DELTET/EB constant.
    pub fn eb(&self) -> f64 {
        self.inner.eb
    }

    /// Access the DELTET/M constants.
    pub fn m(&self) -> [f64; 2] {
        self.inner.m
    }

    /// Access the leap second table.
    pub fn leap_seconds(&self) -> &[(f64, f64)] {
        &self.inner.leap_seconds
    }
}

impl From<muad_dib::spice::lsk::LeapSecondData> for LeapSecondData {
    fn from(inner: muad_dib::spice::lsk::LeapSecondData) -> Self {
        LeapSecondData { inner }
    }
}

/// Extract LeapSecondData from a muad-dib SpiceKernel.
pub fn extract_lsk(kernel: &muad_dib::kernel::SpiceKernel) -> Result<LeapSecondData> {
    use muad_dib::spice::lsk::LeapSecondExt;

    kernel
        .lsk_data()
        .map(LeapSecondData::from)
        .map_err(Into::into)
}

/// Try to extract LSK data, distinguishing "absent" from "malformed".
///
/// Returns `Ok(None)` if no LSK data is present in the kernel.
/// Returns `Err` if LSK data is present but malformed.
/// Returns `Ok(Some(...))` on success.
pub fn try_extract_lsk(kernel: &muad_dib::kernel::SpiceKernel) -> Result<Option<LeapSecondData>> {
    use muad_dib::spice::lsk::LeapSecondExt;

    if !kernel.has_lsk() {
        return Ok(None);
    }

    kernel
        .lsk_data()
        .map(|d| Some(LeapSecondData::from(d)))
        .map_err(Into::into)
}
