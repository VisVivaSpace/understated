//! Core newtypes for NAIF identifiers and time representations.

use std::fmt;

/// TDB seconds past J2000 epoch.
///
/// J2000 epoch is January 1, 2000, 12:00:00 TDB.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct EpochTDB(pub f64);

impl EpochTDB {
    /// J2000 epoch (TDB = 0).
    pub const J2000: EpochTDB = EpochTDB(0.0);

    /// Create from TDB seconds past J2000.
    #[inline]
    pub fn from_tdb_seconds(seconds: f64) -> Self {
        EpochTDB(seconds)
    }

    /// Get as TDB seconds past J2000.
    #[inline]
    pub fn as_tdb_seconds(self) -> f64 {
        self.0
    }
}

impl From<f64> for EpochTDB {
    fn from(value: f64) -> Self {
        EpochTDB(value)
    }
}

impl From<muad_dib::types::EpochTDB> for EpochTDB {
    fn from(e: muad_dib::types::EpochTDB) -> Self {
        EpochTDB(e.0)
    }
}

impl From<EpochTDB> for muad_dib::types::EpochTDB {
    fn from(e: EpochTDB) -> Self {
        muad_dib::types::EpochTDB(e.0)
    }
}

impl fmt::Display for EpochTDB {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} TDB", self.0)
    }
}

/// Spacecraft clock ticks.
///
/// SCLK times are instrument-specific tick counts, not convertible to TDB
/// without a SCLK kernel for the specific spacecraft.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Sclk(pub f64);

impl Sclk {
    /// Create from SCLK ticks.
    #[inline]
    pub fn from_ticks(ticks: f64) -> Self {
        Sclk(ticks)
    }

    /// Get as SCLK ticks.
    #[inline]
    pub fn as_ticks(self) -> f64 {
        self.0
    }
}

impl From<f64> for Sclk {
    fn from(value: f64) -> Self {
        Sclk(value)
    }
}

impl From<muad_dib::types::Sclk> for Sclk {
    fn from(s: muad_dib::types::Sclk) -> Self {
        Sclk(s.0)
    }
}

impl From<Sclk> for muad_dib::types::Sclk {
    fn from(s: Sclk) -> Self {
        muad_dib::types::Sclk(s.0)
    }
}

impl fmt::Display for Sclk {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} SCLK", self.0)
    }
}

/// NAIF body/frame identifier.
///
/// NAIF IDs follow conventions:
/// - Planets: x99 (e.g., 399 = Earth)
/// - Barycenters: x (e.g., 3 = Earth-Moon barycenter)
/// - Moons: x0y (e.g., 301 = Moon)
/// - Spacecraft: negative (e.g., -82 = Cassini)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct NaifId(pub i32);

impl NaifId {
    /// Sun
    pub const SUN: NaifId = NaifId(10);
    /// Solar System Barycenter
    pub const SSB: NaifId = NaifId(0);
    /// Earth-Moon Barycenter
    pub const EMB: NaifId = NaifId(3);
    /// Earth
    pub const EARTH: NaifId = NaifId(399);
    /// Moon
    pub const MOON: NaifId = NaifId(301);
    /// Mars Barycenter
    pub const MARS_BC: NaifId = NaifId(4);
    /// Mars
    pub const MARS: NaifId = NaifId(499);

    /// Check if this is a spacecraft (negative ID).
    #[inline]
    pub fn is_spacecraft(self) -> bool {
        self.0 < 0
    }

    /// Check if this is a planet (x99 pattern).
    #[inline]
    pub fn is_planet(self) -> bool {
        self.0 > 0 && self.0 % 100 == 99
    }

    /// Check if this is a barycenter (single digit 1-9 or 0 for SSB).
    #[inline]
    pub fn is_barycenter(self) -> bool {
        self.0 >= 0 && self.0 <= 9
    }
}

impl From<i32> for NaifId {
    fn from(value: i32) -> Self {
        NaifId(value)
    }
}

impl From<muad_dib::types::NaifId> for NaifId {
    fn from(id: muad_dib::types::NaifId) -> Self {
        NaifId(id.0)
    }
}

impl From<NaifId> for muad_dib::types::NaifId {
    fn from(id: NaifId) -> Self {
        muad_dib::types::NaifId(id.0)
    }
}

impl fmt::Display for NaifId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_epoch_tdb() {
        let epoch = EpochTDB::from_tdb_seconds(1e9);
        assert!((epoch.as_tdb_seconds() - 1e9).abs() < 1e-10);
        assert_eq!(format!("{}", EpochTDB(0.0)), "0 TDB");
    }

    #[test]
    fn test_sclk() {
        let sclk = Sclk::from_ticks(123456789.0);
        assert!((sclk.as_ticks() - 123456789.0).abs() < 1e-10);
        assert_eq!(format!("{}", sclk), "123456789 SCLK");
    }

    #[test]
    fn test_naif_id_classification() {
        assert!(NaifId(-82).is_spacecraft());
        assert!(!NaifId(-82).is_planet());

        assert!(NaifId(399).is_planet());
        assert!(!NaifId(399).is_spacecraft());

        assert!(NaifId(3).is_barycenter());
        assert!(NaifId(0).is_barycenter());
        assert!(!NaifId(399).is_barycenter());
    }

    #[test]
    fn test_naif_id_display() {
        assert_eq!(format!("{}", NaifId::EARTH), "399");
        assert_eq!(format!("{}", NaifId::SSB), "0");
    }
}
