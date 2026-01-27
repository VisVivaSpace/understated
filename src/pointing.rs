//! Pointing data: quaternion and optional angular velocity with frame context.
//!
//! Quaternion uses SPICE convention: scalar-first [q0, q1, q2, q3].
//! Angular velocity is in radians per second.

/// Pointing data: quaternion and optional angular velocity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Pointing {
    /// Reference frame code (NAIF frame ID).
    pub frame: i32,
    /// Quaternion [q0, q1, q2, q3] (scalar-first).
    pub quaternion: [f64; 4],
    /// Angular velocity [wx, wy, wz] in rad/s, if available.
    pub angular_velocity: Option<[f64; 3]>,
}

impl Pointing {
    /// Create a new pointing with full context.
    #[inline]
    pub fn new(frame: i32, quaternion: [f64; 4], angular_velocity: Option<[f64; 3]>) -> Self {
        Pointing {
            frame,
            quaternion,
            angular_velocity,
        }
    }

    /// Create a raw pointing without frame context (frame = 0).
    #[inline]
    pub(crate) fn new_raw(quaternion: [f64; 4], angular_velocity: Option<[f64; 3]>) -> Self {
        Pointing {
            frame: 0,
            quaternion,
            angular_velocity,
        }
    }

    /// Create pointing from quaternion only.
    #[inline]
    pub fn from_quaternion(frame: i32, quaternion: [f64; 4]) -> Self {
        Pointing {
            frame,
            quaternion,
            angular_velocity: None,
        }
    }

    /// Scalar component of the quaternion.
    #[inline]
    pub fn scalar(&self) -> f64 {
        self.quaternion[0]
    }

    /// Vector component of the quaternion.
    #[inline]
    pub fn vector(&self) -> [f64; 3] {
        [self.quaternion[1], self.quaternion[2], self.quaternion[3]]
    }

    /// Check if the quaternion is normalized.
    pub fn is_normalized(&self) -> bool {
        let [q0, q1, q2, q3] = self.quaternion;
        let norm_sq = q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3;
        (norm_sq - 1.0).abs() < 1e-10
    }

    /// Return a copy with normalized quaternion.
    pub fn normalize(&self) -> Pointing {
        let [q0, q1, q2, q3] = self.quaternion;
        let norm = (q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3).sqrt();

        if norm < 1e-15 {
            return *self;
        }

        Pointing {
            frame: self.frame,
            quaternion: [q0 / norm, q1 / norm, q2 / norm, q3 / norm],
            angular_velocity: self.angular_velocity,
        }
    }
}

impl Default for Pointing {
    fn default() -> Self {
        Pointing {
            frame: 0,
            quaternion: [1.0, 0.0, 0.0, 0.0],
            angular_velocity: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const J2000: i32 = 1;

    #[test]
    fn test_normalized() {
        let p = Pointing::from_quaternion(J2000, [1.0, 0.0, 0.0, 0.0]);
        assert!(p.is_normalized());

        let unnorm = Pointing::from_quaternion(J2000, [2.0, 0.0, 0.0, 0.0]);
        assert!(!unnorm.is_normalized());

        let normalized = unnorm.normalize();
        assert!(normalized.is_normalized());
        assert!((normalized.quaternion[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_components() {
        let p = Pointing::new(J2000, [0.5, 0.5, 0.5, 0.5], Some([0.1, 0.2, 0.3]));
        assert!((p.scalar() - 0.5).abs() < 1e-10);
        assert_eq!(p.vector(), [0.5, 0.5, 0.5]);
        assert!(p.angular_velocity.is_some());
    }

    #[test]
    fn test_default() {
        let pointing = Pointing::default();
        assert!(pointing.is_normalized());
        assert_eq!(pointing.scalar(), 1.0);
        assert_eq!(pointing.frame, 0);
    }
}
