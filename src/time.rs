//! Time string parsing and formatting for TDB epochs.
//!
//! Delegates to muad-dib's time module for all parsing and formatting.
//! Re-exports `TimeFormat`, `format_iso8601`, `format_calendar`, and `tdb_to_calendar`.

use crate::error::{Error, Result};
use crate::types::EpochTDB;

pub use muad_dib::spice::time::{format_calendar, format_iso8601, tdb_to_calendar, TimeFormat};

impl EpochTDB {
    /// Parse a time string to TDB seconds past J2000.
    ///
    /// Auto-detects format. For UTC strings, use `LeapSecondData::utc_to_tdb_seconds()`
    /// to properly account for leap seconds.
    pub fn parse(time_str: &str) -> Result<EpochTDB> {
        muad_dib::types::EpochTDB::parse(time_str)
            .map(EpochTDB::from)
            .map_err(|_| Error::TimeParseError {
                input: time_str.to_string(),
            })
    }

    /// Parse with an explicit format hint.
    pub fn parse_with_format(time_str: &str, format: TimeFormat) -> Result<EpochTDB> {
        muad_dib::types::EpochTDB::parse_with_format(time_str, format)
            .map(EpochTDB::from)
            .map_err(|_| Error::TimeParseError {
                input: time_str.to_string(),
            })
    }
}
