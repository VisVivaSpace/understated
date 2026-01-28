//! State vector: position and velocity with full relativity context.
//!
//! A state vector is only meaningful when you know:
//! - **target**: What body this state describes
//! - **center**: The origin of the coordinate system
//! - **frame**: The reference frame defining axis orientation
//!
//! # State Arithmetic
//!
//! - **Addition** (chain): `(SSB->Earth) + (Earth->Moon) = SSB->Moon`
//! - **Subtraction** (relative): `(Mars rel SSB) - (Earth rel SSB) = Mars rel Earth`
//! - **Negation** (reverse): `-(Earth->Moon) = Moon->Earth`

use std::ops::{Add, Neg, Sub};

use crate::error::{Error, Result};
use crate::types::{FrameId, NaifId};

/// State vector: position and velocity with full relativity context.
///
/// Position is in km, velocity is in km/s.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct State {
    /// Target body this state describes.
    pub target: NaifId,
    /// Center body (origin of coordinates).
    pub center: NaifId,
    /// Reference frame code (NAIF frame ID, e.g., J2000 = 1).
    pub frame: FrameId,
    /// Position vector [x, y, z] in km.
    pub position: [f64; 3],
    /// Velocity vector [vx, vy, vz] in km/s.
    pub velocity: [f64; 3],
}

impl State {
    /// Create a new state with full context.
    #[inline]
    pub fn new(
        target: NaifId,
        center: NaifId,
        frame: FrameId,
        position: [f64; 3],
        velocity: [f64; 3],
    ) -> Self {
        State {
            target,
            center,
            frame,
            position,
            velocity,
        }
    }

    /// Create a raw state without relativity context.
    ///
    /// For intermediate calculations only. Target, center, frame are set to
    /// placeholder values.
    #[inline]
    pub(crate) fn new_raw(position: [f64; 3], velocity: [f64; 3]) -> Self {
        State {
            target: NaifId(0),
            center: NaifId(0),
            frame: FrameId(0),
            position,
            velocity,
        }
    }

    /// Create a state with zero velocity.
    #[inline]
    pub fn from_position(target: NaifId, center: NaifId, frame: FrameId, position: [f64; 3]) -> Self {
        State {
            target,
            center,
            frame,
            position,
            velocity: [0.0, 0.0, 0.0],
        }
    }

    /// Position magnitude (distance from center).
    #[inline]
    pub fn distance(&self) -> f64 {
        let [x, y, z] = self.position;
        (x * x + y * y + z * z).sqrt()
    }

    /// Velocity magnitude (speed).
    #[inline]
    pub fn speed(&self) -> f64 {
        let [vx, vy, vz] = self.velocity;
        (vx * vx + vy * vy + vz * vz).sqrt()
    }

    /// Chain traversal returning `Result` instead of panicking.
    ///
    /// `(center->A) + (A->target) = (center->target)`
    ///
    /// Returns an error if frames don't match or the chain is invalid
    /// (i.e., `self.target != other.center`).
    pub fn try_add(&self, other: &State) -> Result<State> {
        if self.frame != other.frame {
            return Err(Error::FrameMismatch {
                frame_a: self.frame,
                frame_b: other.frame,
            });
        }
        if self.target != other.center {
            return Err(Error::InvalidChain {
                self_target: self.target,
                other_center: other.center,
            });
        }

        Ok(State {
            target: other.target,
            center: self.center,
            frame: self.frame,
            position: [
                self.position[0] + other.position[0],
                self.position[1] + other.position[1],
                self.position[2] + other.position[2],
            ],
            velocity: [
                self.velocity[0] + other.velocity[0],
                self.velocity[1] + other.velocity[1],
                self.velocity[2] + other.velocity[2],
            ],
        })
    }

    /// Relative motion returning `Result` instead of panicking.
    ///
    /// `(center->A) - (center->B) = (B->A)`
    ///
    /// Returns an error if frames or centers don't match.
    pub fn try_sub(&self, other: &State) -> Result<State> {
        if self.frame != other.frame {
            return Err(Error::FrameMismatch {
                frame_a: self.frame,
                frame_b: other.frame,
            });
        }
        if self.center != other.center {
            return Err(Error::CenterMismatch {
                center_a: self.center,
                center_b: other.center,
            });
        }

        Ok(State {
            target: self.target,
            center: other.target,
            frame: self.frame,
            position: [
                self.position[0] - other.position[0],
                self.position[1] - other.position[1],
                self.position[2] - other.position[2],
            ],
            velocity: [
                self.velocity[0] - other.velocity[0],
                self.velocity[1] - other.velocity[1],
                self.velocity[2] - other.velocity[2],
            ],
        })
    }
}

/// Chain traversal: `(center->A) + (A->target) = (center->target)`
///
/// # Panics
///
/// Panics if `self.frame != other.frame` or `self.target != other.center`.
/// Use [`State::try_add`] for a non-panicking alternative.
impl Add<&State> for State {
    type Output = State;

    fn add(self, other: &State) -> State {
        assert_eq!(
            self.frame, other.frame,
            "Cannot add states in different frames: {} vs {}",
            self.frame, other.frame
        );
        assert_eq!(
            self.target, other.center,
            "Invalid chain: self.target ({}) != other.center ({})",
            self.target, other.center
        );

        State {
            target: other.target,
            center: self.center,
            frame: self.frame,
            position: [
                self.position[0] + other.position[0],
                self.position[1] + other.position[1],
                self.position[2] + other.position[2],
            ],
            velocity: [
                self.velocity[0] + other.velocity[0],
                self.velocity[1] + other.velocity[1],
                self.velocity[2] + other.velocity[2],
            ],
        }
    }
}

impl Add for State {
    type Output = State;
    #[allow(clippy::op_ref)]
    fn add(self, other: State) -> State {
        self + &other
    }
}

/// Relative motion: `(center->A) - (center->B) = (B->A)`
///
/// # Panics
///
/// Panics if `self.frame != other.frame` or `self.center != other.center`.
/// Use [`State::try_sub`] for a non-panicking alternative.
impl Sub<&State> for State {
    type Output = State;

    fn sub(self, other: &State) -> State {
        assert_eq!(
            self.frame, other.frame,
            "Cannot subtract states in different frames: {} vs {}",
            self.frame, other.frame
        );
        assert_eq!(
            self.center, other.center,
            "Cannot subtract states with different centers: {} vs {}",
            self.center, other.center
        );

        State {
            target: self.target,
            center: other.target,
            frame: self.frame,
            position: [
                self.position[0] - other.position[0],
                self.position[1] - other.position[1],
                self.position[2] - other.position[2],
            ],
            velocity: [
                self.velocity[0] - other.velocity[0],
                self.velocity[1] - other.velocity[1],
                self.velocity[2] - other.velocity[2],
            ],
        }
    }
}

impl Sub for State {
    type Output = State;
    #[allow(clippy::op_ref)]
    fn sub(self, other: State) -> State {
        self - &other
    }
}

/// Reverse direction: `-(center->target)` becomes `(target->center)`.
impl Neg for State {
    type Output = State;

    fn neg(self) -> State {
        State {
            target: self.center,
            center: self.target,
            frame: self.frame,
            position: [-self.position[0], -self.position[1], -self.position[2]],
            velocity: [-self.velocity[0], -self.velocity[1], -self.velocity[2]],
        }
    }
}

/// Default state: SSB relative to SSB in frame 0 at the origin.
///
/// This is a zero-state sentinel — useful as an accumulator initial value.
impl Default for State {
    fn default() -> Self {
        State {
            target: NaifId(0),
            center: NaifId(0),
            frame: FrameId(0),
            position: [0.0, 0.0, 0.0],
            velocity: [0.0, 0.0, 0.0],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SSB: NaifId = NaifId(0);
    const EARTH: NaifId = NaifId(399);
    const MOON: NaifId = NaifId(301);
    const MARS: NaifId = NaifId(499);
    const J2000: FrameId = FrameId(1);

    #[test]
    fn test_distance_and_speed() {
        let state = State::new(EARTH, SSB, J2000, [3.0, 4.0, 0.0], [0.0, 0.0, 0.0]);
        assert!((state.distance() - 5.0).abs() < 1e-10);

        let state = State::new(EARTH, SSB, J2000, [0.0, 0.0, 0.0], [3.0, 4.0, 0.0]);
        assert!((state.speed() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_add_chain() {
        let ssb_to_earth = State::new(EARTH, SSB, J2000, [1.0, 2.0, 3.0], [0.1, 0.2, 0.3]);
        let earth_to_moon = State::new(MOON, EARTH, J2000, [4.0, 5.0, 6.0], [0.4, 0.5, 0.6]);
        let ssb_to_moon = ssb_to_earth + earth_to_moon;

        assert!((ssb_to_moon.position[0] - 5.0).abs() < 1e-10);
        assert!((ssb_to_moon.velocity[2] - 0.9).abs() < 1e-10);
        assert_eq!(ssb_to_moon.target, MOON);
        assert_eq!(ssb_to_moon.center, SSB);
    }

    #[test]
    fn test_sub_relative() {
        let ssb_to_mars = State::new(MARS, SSB, J2000, [5.0, 7.0, 9.0], [0.5, 0.7, 0.9]);
        let ssb_to_earth = State::new(EARTH, SSB, J2000, [1.0, 2.0, 3.0], [0.1, 0.2, 0.3]);
        let earth_to_mars = ssb_to_mars - ssb_to_earth;

        assert!((earth_to_mars.position[0] - 4.0).abs() < 1e-10);
        assert!((earth_to_mars.velocity[2] - 0.6).abs() < 1e-10);
        assert_eq!(earth_to_mars.target, MARS);
        assert_eq!(earth_to_mars.center, EARTH);
    }

    #[test]
    fn test_negate() {
        let ssb_to_earth =
            State::new(EARTH, SSB, J2000, [1.0, -2.0, 3.0], [-0.1, 0.2, -0.3]);
        let earth_to_ssb = -ssb_to_earth;

        assert!((earth_to_ssb.position[0] + 1.0).abs() < 1e-10);
        assert!((earth_to_ssb.position[1] - 2.0).abs() < 1e-10);
        assert_eq!(earth_to_ssb.target, SSB);
        assert_eq!(earth_to_ssb.center, EARTH);
    }

    #[test]
    fn test_add_ref() {
        let ssb_to_earth = State::new(EARTH, SSB, J2000, [1.0, 2.0, 3.0], [0.1, 0.2, 0.3]);
        let earth_to_moon = State::new(MOON, EARTH, J2000, [4.0, 5.0, 6.0], [0.4, 0.5, 0.6]);
        let sum = ssb_to_earth + &earth_to_moon;
        assert!((sum.position[0] - 5.0).abs() < 1e-10);
    }

    #[test]
    #[should_panic(expected = "Cannot add states in different frames")]
    fn test_add_frame_mismatch() {
        let s1 = State::new(EARTH, SSB, FrameId(1), [1.0, 2.0, 3.0], [0.1, 0.2, 0.3]);
        let s2 = State::new(MOON, EARTH, FrameId(2), [4.0, 5.0, 6.0], [0.4, 0.5, 0.6]);
        let _ = s1 + s2;
    }

    #[test]
    #[should_panic(expected = "Invalid chain")]
    fn test_add_invalid_chain() {
        let s1 = State::new(EARTH, SSB, J2000, [1.0, 2.0, 3.0], [0.1, 0.2, 0.3]);
        let s2 = State::new(MOON, SSB, J2000, [4.0, 5.0, 6.0], [0.4, 0.5, 0.6]);
        let _ = s1 + s2;
    }

    #[test]
    #[should_panic(expected = "Cannot subtract states with different centers")]
    fn test_sub_center_mismatch() {
        let s1 = State::new(MARS, SSB, J2000, [1.0, 2.0, 3.0], [0.1, 0.2, 0.3]);
        let s2 = State::new(MOON, EARTH, J2000, [4.0, 5.0, 6.0], [0.4, 0.5, 0.6]);
        let _ = s1 - s2;
    }

    #[test]
    fn test_default() {
        let state = State::default();
        assert_eq!(state.distance(), 0.0);
        assert_eq!(state.target, NaifId(0));
        assert_eq!(state.frame, FrameId(0));
    }

    #[test]
    fn test_try_add_ok() {
        let ssb_to_earth = State::new(EARTH, SSB, J2000, [1.0, 2.0, 3.0], [0.1, 0.2, 0.3]);
        let earth_to_moon = State::new(MOON, EARTH, J2000, [4.0, 5.0, 6.0], [0.4, 0.5, 0.6]);
        let result = ssb_to_earth.try_add(&earth_to_moon).unwrap();
        assert!((result.position[0] - 5.0).abs() < 1e-10);
        assert_eq!(result.target, MOON);
        assert_eq!(result.center, SSB);
    }

    #[test]
    fn test_try_add_frame_mismatch() {
        let s1 = State::new(EARTH, SSB, FrameId(1), [1.0, 2.0, 3.0], [0.1, 0.2, 0.3]);
        let s2 = State::new(MOON, EARTH, FrameId(2), [4.0, 5.0, 6.0], [0.4, 0.5, 0.6]);
        assert!(s1.try_add(&s2).is_err());
    }

    #[test]
    fn test_try_add_invalid_chain() {
        let s1 = State::new(EARTH, SSB, J2000, [1.0, 2.0, 3.0], [0.1, 0.2, 0.3]);
        let s2 = State::new(MOON, SSB, J2000, [4.0, 5.0, 6.0], [0.4, 0.5, 0.6]);
        assert!(s1.try_add(&s2).is_err());
    }

    #[test]
    fn test_try_sub_ok() {
        let ssb_to_mars = State::new(MARS, SSB, J2000, [5.0, 7.0, 9.0], [0.5, 0.7, 0.9]);
        let ssb_to_earth = State::new(EARTH, SSB, J2000, [1.0, 2.0, 3.0], [0.1, 0.2, 0.3]);
        let result = ssb_to_mars.try_sub(&ssb_to_earth).unwrap();
        assert!((result.position[0] - 4.0).abs() < 1e-10);
        assert_eq!(result.target, MARS);
        assert_eq!(result.center, EARTH);
    }

    #[test]
    fn test_try_sub_center_mismatch() {
        let s1 = State::new(MARS, SSB, J2000, [1.0, 2.0, 3.0], [0.1, 0.2, 0.3]);
        let s2 = State::new(MOON, EARTH, J2000, [4.0, 5.0, 6.0], [0.4, 0.5, 0.6]);
        assert!(s1.try_sub(&s2).is_err());
    }
}
