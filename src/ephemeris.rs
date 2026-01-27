//! High-level facade for querying SPICE kernel data.
//!
//! The [`Ephemeris`] struct wraps `muad_dib::SpiceKernel` and provides the
//! primary user-facing API for state and attitude evaluation.
//!
//! # Example
//!
//! ```ignore
//! use understated::{Ephemeris, EpochTDB, NaifId, Sclk, TimeFormat};
//!
//! let eph = Ephemeris::load("de440s.bsp")?;
//!
//! let epoch = EpochTDB::parse("2020-01-01T12:00:00")?;
//! let state = eph.state_of(NaifId::EARTH, epoch, NaifId::SSB)?;
//! println!("Earth position: {:?} km", state.position);
//! ```

use std::path::Path;

use muad_dib::kernel::SpiceKernel;

use crate::ck::evaluate_ck;
use crate::error::{Error, Result};
use crate::lsk::{extract_lsk, LeapSecondData};
use crate::pointing::Pointing;
use crate::spk::chain;
use crate::state::State;
use crate::time::TimeFormat;
use crate::types::{EpochTDB, NaifId, Sclk};

/// High-level entry point for querying SPICE ephemeris and pointing data.
///
/// Wraps a `muad_dib::SpiceKernel` and provides clean methods for:
/// - State evaluation with automatic body chaining
/// - CK pointing evaluation
/// - UTC/TDB time conversion (when LSK data is loaded)
/// - Discovery of available bodies and instruments
pub struct Ephemeris {
    kernel: SpiceKernel,
    lsk: Option<LeapSecondData>,
}

impl Ephemeris {
    /// Load a single SPICE kernel file.
    ///
    /// Supports SPK (.bsp), CK (.bc), text PCK (.tpc), and LSK (.tls) files.
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let kernel = SpiceKernel::load(path)?;
        let lsk = extract_lsk(&kernel).ok();
        Ok(Ephemeris { kernel, lsk })
    }

    /// Load multiple SPICE kernel files.
    pub fn load_many<P: AsRef<Path>>(paths: &[P]) -> Result<Self> {
        let kernel = SpiceKernel::load_many(paths)?;
        let lsk = extract_lsk(&kernel).ok();
        Ok(Ephemeris { kernel, lsk })
    }

    /// Get the state of `target` relative to `center` at `epoch`.
    ///
    /// Automatically chains through intermediate center bodies when the SPK
    /// file doesn't directly provide the requested pair.
    ///
    /// # Arguments
    ///
    /// * `target` - NAIF ID of the target body
    /// * `epoch` - TDB seconds past J2000
    /// * `center` - NAIF ID of the center body (origin for state vector)
    pub fn state_of(&self, target: NaifId, epoch: EpochTDB, center: NaifId) -> Result<State> {
        chain::state_of(&self.kernel, target, epoch, center)
    }

    /// Get pointing for an instrument at the given SCLK time.
    ///
    /// Searches loaded CK data for a segment covering the instrument at the
    /// given spacecraft clock time.
    pub fn pointing_of(&self, instrument: NaifId, sclk: Sclk) -> Result<Pointing> {
        let md_instrument = muad_dib::types::NaifId(instrument.0);
        let sclk_ticks = sclk.as_ticks();

        let segment = self
            .kernel
            .ck_segments_for(md_instrument)
            .find(|seg| seg.initial_sclk <= sclk_ticks && sclk_ticks <= seg.final_sclk)
            .ok_or(Error::NoCoverage {
                body: instrument,
                epoch: sclk_ticks,
            })?;

        let view = self.kernel.ck_view(segment);
        let data = view.data();
        let mut pointing = evaluate_ck(data, sclk_ticks)?;
        pointing.frame = segment.frame_code;

        Ok(pointing)
    }

    /// Get all body IDs with SPK coverage.
    pub fn spk_bodies(&self) -> Vec<NaifId> {
        self.kernel
            .spk_bodies()
            .into_iter()
            .map(NaifId::from)
            .collect()
    }

    /// Get all instrument IDs with CK coverage.
    pub fn ck_instruments(&self) -> Vec<NaifId> {
        self.kernel
            .ck_instruments()
            .into_iter()
            .map(NaifId::from)
            .collect()
    }

    /// Convert a UTC time string to TDB epoch.
    ///
    /// Requires that an LSK (leap second kernel) file was loaded.
    pub fn utc_to_tdb(&self, utc: &str) -> Result<EpochTDB> {
        let lsk = self.lsk.as_ref().ok_or(Error::MissingLskData)?;
        lsk.utc_to_tdb(utc)
    }

    /// Convert a TDB epoch to a formatted UTC string.
    ///
    /// Requires that an LSK (leap second kernel) file was loaded.
    pub fn tdb_to_utc(&self, tdb: EpochTDB, format: TimeFormat) -> Result<String> {
        let lsk = self.lsk.as_ref().ok_or(Error::MissingLskData)?;
        lsk.tdb_to_utc(tdb, format)
    }

    /// Access the leap second data, if loaded.
    pub fn lsk_data(&self) -> Option<&LeapSecondData> {
        self.lsk.as_ref()
    }
}
