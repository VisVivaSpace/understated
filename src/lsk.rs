//! Leap second data and TDB/UTC conversion.
//!
//! # NAIF Time System Relationships
//!
//! ```text
//! UTC -> TAI -> TT -> TDB
//!      +dAT  +32.184s  +periodic
//! ```

use crate::error::{Error, Result};
use crate::time::{format_calendar, format_iso8601, TimeFormat};
use crate::types::EpochTDB;

/// Leap second data extracted from an LSK file.
#[derive(Debug, Clone)]
pub struct LeapSecondData {
    /// DELTET/DELTA_T_A: typically 32.184 seconds.
    pub delta_t_a: f64,
    /// DELTET/K: TDB-TT relationship constant.
    pub k: f64,
    /// DELTET/EB: Earth's orbital eccentricity effect.
    pub eb: f64,
    /// DELTET/M: Mean anomaly constants [m0, m1].
    pub m: [f64; 2],
    /// Leap second entries: (TAI-UTC offset, epoch in TDB seconds past J2000).
    /// Sorted by epoch ascending.
    pub leap_seconds: Vec<(f64, f64)>,
}

impl LeapSecondData {
    /// Get the TAI-UTC offset (dAT) for a given TDB epoch.
    pub fn delta_at(&self, tdb: f64) -> f64 {
        let mut delta = 0.0;
        for (d, epoch) in &self.leap_seconds {
            if tdb >= *epoch {
                delta = *d;
            } else {
                break;
            }
        }
        delta
    }

    /// Convert TDB to TT.
    pub fn tdb_to_tt(&self, tdb: f64) -> f64 {
        let m = self.m[0] + self.m[1] * tdb;
        let e = m + self.eb * m.sin();
        let periodic = self.k * e.sin();
        tdb - periodic
    }

    /// Convert TT to TDB (iterative).
    pub fn tt_to_tdb(&self, tt: f64) -> f64 {
        let mut tdb = tt;
        for _ in 0..3 {
            let m = self.m[0] + self.m[1] * tdb;
            let e = m + self.eb * m.sin();
            let periodic = self.k * e.sin();
            tdb = tt + periodic;
        }
        tdb
    }

    /// Convert TDB to TAI.
    pub fn tdb_to_tai(&self, tdb: f64) -> f64 {
        let tt = self.tdb_to_tt(tdb);
        tt - self.delta_t_a
    }

    /// Convert TAI to TDB.
    pub fn tai_to_tdb(&self, tai: f64) -> f64 {
        let tt = tai + self.delta_t_a;
        self.tt_to_tdb(tt)
    }

    /// Convert TDB to UTC seconds past J2000.
    pub fn tdb_to_utc_seconds(&self, tdb: f64) -> f64 {
        let tai = self.tdb_to_tai(tdb);
        let delta_at = self.delta_at(tdb);
        tai - delta_at
    }

    /// Convert UTC seconds past J2000 to TDB.
    pub fn utc_to_tdb_seconds(&self, utc: f64) -> f64 {
        let tdb_approx = utc + self.delta_t_a + 10.0;
        let delta_at = self.delta_at(tdb_approx);
        let tai = utc + delta_at;
        self.tai_to_tdb(tai)
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
}

/// Extract LeapSecondData from a muad-dib SpiceKernel.
pub fn extract_lsk(kernel: &muad_dib::kernel::SpiceKernel) -> Result<LeapSecondData> {
    use muad_dib::spice::lsk::LeapSecondExt;

    kernel.lsk_data().ok_or(Error::MissingLskData)
        .map(|lsk| LeapSecondData {
            delta_t_a: lsk.delta_t_a,
            k: lsk.k,
            eb: lsk.eb,
            m: lsk.m,
            leap_seconds: lsk.leap_seconds,
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_lsk() -> LeapSecondData {
        // Approximate values from naif0012.tls
        let epoch_1972 = EpochTDB::parse("1972-JAN-01 00:00:00").unwrap().0;
        let epoch_1972_jul = EpochTDB::parse("1972-JUL-01 00:00:00").unwrap().0;
        let epoch_2017 = EpochTDB::parse("2017-JAN-01 00:00:00").unwrap().0;

        LeapSecondData {
            delta_t_a: 32.184,
            k: 1.657e-3,
            eb: 1.671e-2,
            m: [6.239996, 1.99096871e-7],
            leap_seconds: vec![
                (10.0, epoch_1972),
                (11.0, epoch_1972_jul),
                (37.0, epoch_2017),
            ],
        }
    }

    #[test]
    fn test_delta_at_lookup() {
        let lsk = make_test_lsk();
        assert_eq!(lsk.delta_at(-1e10), 0.0);
        assert_eq!(lsk.delta_at(1e9), 37.0);
    }

    #[test]
    fn test_tdb_tt_round_trip() {
        let lsk = make_test_lsk();
        let tdb = 0.0;
        let tt = lsk.tdb_to_tt(tdb);
        assert!((tdb - tt).abs() < 0.002);

        let tdb_back = lsk.tt_to_tdb(tt);
        assert!((tdb - tdb_back).abs() < 1e-9);
    }

    #[test]
    fn test_utc_tdb_round_trip() {
        let lsk = make_test_lsk();
        let utc_seconds = 0.0;
        let tdb = lsk.utc_to_tdb_seconds(utc_seconds);
        let utc_back = lsk.tdb_to_utc_seconds(tdb);
        assert!((utc_seconds - utc_back).abs() < 1.0);
    }

    #[test]
    fn test_utc_to_tdb_string() {
        let lsk = make_test_lsk();
        let tdb = lsk.utc_to_tdb("2000-01-01T12:00:00").unwrap();
        assert!(tdb.0 > 0.0);
        assert!(tdb.0 < 100.0);
    }

    #[test]
    fn test_tdb_to_utc_string() {
        let lsk = make_test_lsk();
        let utc = lsk.tdb_to_utc(EpochTDB(0.0), TimeFormat::Iso8601).unwrap();
        assert!(utc.starts_with("2000-01-01T"));
    }
}
