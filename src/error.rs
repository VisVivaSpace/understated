//! Error types for interpolation and time operations.

use crate::types::NaifId;

#[derive(thiserror::Error, Debug)]
pub enum Error {
    /// No coverage for body/instrument at the requested time.
    #[error("No coverage for body {body} at epoch {epoch}")]
    NoCoverage { body: NaifId, epoch: f64 },

    /// SPK segment type interpolation not implemented.
    #[error("SPK type {spk_type} interpolation not implemented")]
    UnsupportedSpkType { spk_type: i32 },

    /// CK segment type interpolation not implemented.
    #[error("CK type {ck_type} interpolation not implemented")]
    UnsupportedCkType { ck_type: i32 },

    /// Cannot parse time string.
    #[error("Cannot parse time string: '{input}'{}", reason.as_ref().map(|r| format!(": {r}")).unwrap_or_default())]
    TimeParseError {
        input: String,
        /// Upstream parse failure reason, if available.
        reason: Option<String>,
    },

    /// Leap second kernel (LSK) data required for TDB/UTC conversion.
    #[error("Leap second kernel (LSK) data required for TDB/UTC conversion")]
    MissingLskData,

    /// Epoch is outside the valid range for interpolation.
    #[error("Epoch {epoch} is outside segment coverage [{start}, {end}]")]
    EpochOutOfRange { epoch: f64, start: f64, end: f64 },

    /// I/O error from file operations.
    #[error(transparent)]
    Io(#[from] std::io::Error),

    /// Wrapped error from muad-dib.
    #[error("Kernel error: {0}")]
    Kernel(#[from] muad_dib::error::Error),
}

pub type Result<T> = std::result::Result<T, Error>;
