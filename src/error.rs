//! Error types for interpolation and time operations.

use crate::types::{FrameId, NaifId};

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

    /// States are in different reference frames.
    #[error("Cannot combine states in different frames: {frame_a} vs {frame_b}")]
    FrameMismatch { frame_a: FrameId, frame_b: FrameId },

    /// Invalid chain: self.target does not match other.center.
    #[error("Invalid chain: self.target ({self_target}) != other.center ({other_center})")]
    InvalidChain { self_target: NaifId, other_center: NaifId },

    /// States have different center bodies.
    #[error("Cannot subtract states with different centers: {center_a} vs {center_b}")]
    CenterMismatch { center_a: NaifId, center_b: NaifId },

    /// Center chain depth limit exceeded.
    #[error("Center chain depth limit ({limit}) exceeded for body {body}")]
    ChainDepthExceeded { body: NaifId, limit: usize },

    /// I/O error from file operations.
    #[error(transparent)]
    Io(#[from] std::io::Error),

    /// Wrapped error from muad-dib.
    #[error("Kernel error: {0}")]
    Kernel(#[from] muad_dib::error::Error),
}

pub type Result<T> = std::result::Result<T, Error>;
