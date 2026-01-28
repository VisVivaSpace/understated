//! **understated** — Spacecraft state and attitude interpolation from SPICE kernel data.
//!
//! This crate provides interpolation algorithms for evaluating SPK (state) and
//! CK (attitude) data from NAIF SPICE kernels. It wraps `muad_dib` for file I/O
//! and exposes a clean public API through the [`Ephemeris`] struct.
//!
//! # Supported Types
//!
//! **SPK (state):**
//! - Type 2: Chebyshev (position only, velocity via differentiation)
//! - Type 3: Chebyshev (position + velocity coefficients)
//! - Type 5: Two-body/Keplerian propagation
//! - Type 8: Lagrange (equal time steps)
//! - Type 9: Lagrange (unequal time steps)
//! - Type 13: Hermite (unequal time steps)
//!
//! **CK (attitude):**
//! - Type 1: Discrete pointing (nearest record)
//! - Type 3: Linear interpolation (rotation matrix, axis-angle; matches CSPICE CKE03)

pub mod error;
pub mod types;
pub mod state;
pub mod pointing;
pub mod coord;
pub mod time;
pub mod lsk;
pub mod rotation;
pub mod spk;
pub mod ck;
pub mod ephemeris;

pub use error::{Error, Result};
pub use types::{EpochTDB, FrameId, NaifId, Sclk};
pub use state::State;
pub use pointing::Pointing;
pub use coord::{Rectangular, Latitudinal, Spherical, Cylindrical};
pub use time::TimeFormat;
pub use lsk::LeapSecondData;
pub use ephemeris::Ephemeris;
