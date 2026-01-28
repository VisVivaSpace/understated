//! Coordinate system conversions.
//!
//! Provides type-safe conversions between coordinate systems:
//! - **Rectangular** (Cartesian): x, y, z
//! - **Latitudinal**: radius, longitude, latitude
//! - **Spherical**: radius, colatitude, longitude
//! - **Cylindrical**: r (radial), longitude, z
//!
//! All angles are in radians. Conversions follow NAIF CSPICE conventions.

/// Rectangular (Cartesian) coordinates.
///
/// Components are [x, y, z] in the same units as the source data (typically km).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rectangular(pub [f64; 3]);

impl Rectangular {
    #[inline]
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Rectangular([x, y, z])
    }

    #[inline]
    pub fn x(&self) -> f64 {
        self.0[0]
    }
    #[inline]
    pub fn y(&self) -> f64 {
        self.0[1]
    }
    #[inline]
    pub fn z(&self) -> f64 {
        self.0[2]
    }

    /// Euclidean norm.
    #[inline]
    pub fn magnitude(&self) -> f64 {
        (self.0[0].powi(2) + self.0[1].powi(2) + self.0[2].powi(2)).sqrt()
    }

    /// Convert to latitudinal coordinates.
    pub fn to_latitudinal(&self) -> Latitudinal {
        let [x, y, z] = self.0;
        let radius = self.magnitude();

        if radius == 0.0 {
            return Latitudinal {
                radius: 0.0,
                longitude: 0.0,
                latitude: 0.0,
            };
        }

        Latitudinal {
            radius,
            longitude: y.atan2(x),
            latitude: (z / radius).clamp(-1.0, 1.0).asin(),
        }
    }

    /// Convert to spherical coordinates.
    pub fn to_spherical(&self) -> Spherical {
        let [x, y, z] = self.0;
        let radius = self.magnitude();

        if radius == 0.0 {
            return Spherical {
                radius: 0.0,
                colatitude: 0.0,
                longitude: 0.0,
            };
        }

        Spherical {
            radius,
            colatitude: (z / radius).clamp(-1.0, 1.0).acos(),
            longitude: y.atan2(x),
        }
    }

    /// Convert to cylindrical coordinates.
    pub fn to_cylindrical(&self) -> Cylindrical {
        let [x, y, z] = self.0;
        Cylindrical {
            r: (x.powi(2) + y.powi(2)).sqrt(),
            longitude: y.atan2(x),
            z,
        }
    }
}

/// Latitudinal coordinates.
///
/// - `radius`: Distance from origin
/// - `longitude`: Angle in xy-plane from +X, range (-pi, pi]
/// - `latitude`: Angle from xy-plane toward +Z, range [-pi/2, pi/2]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Latitudinal {
    pub radius: f64,
    pub longitude: f64,
    pub latitude: f64,
}

impl Latitudinal {
    #[inline]
    pub fn new(radius: f64, longitude: f64, latitude: f64) -> Self {
        Latitudinal {
            radius,
            longitude,
            latitude,
        }
    }

    #[inline]
    pub fn longitude_deg(&self) -> f64 {
        self.longitude.to_degrees()
    }

    #[inline]
    pub fn latitude_deg(&self) -> f64 {
        self.latitude.to_degrees()
    }
}

impl From<Latitudinal> for Rectangular {
    fn from(lat: Latitudinal) -> Self {
        let cos_lat = lat.latitude.cos();
        Rectangular([
            lat.radius * cos_lat * lat.longitude.cos(),
            lat.radius * cos_lat * lat.longitude.sin(),
            lat.radius * lat.latitude.sin(),
        ])
    }
}

/// Spherical coordinates.
///
/// - `radius`: Distance from origin
/// - `colatitude`: Angle from +Z axis, range [0, pi]
/// - `longitude`: Angle in xy-plane from +X, range (-pi, pi]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Spherical {
    pub radius: f64,
    pub colatitude: f64,
    pub longitude: f64,
}

impl Spherical {
    #[inline]
    pub fn new(radius: f64, colatitude: f64, longitude: f64) -> Self {
        Spherical {
            radius,
            colatitude,
            longitude,
        }
    }

    #[inline]
    pub fn colatitude_deg(&self) -> f64 {
        self.colatitude.to_degrees()
    }

    #[inline]
    pub fn longitude_deg(&self) -> f64 {
        self.longitude.to_degrees()
    }
}

impl From<Spherical> for Rectangular {
    fn from(sph: Spherical) -> Self {
        let sin_co = sph.colatitude.sin();
        Rectangular([
            sph.radius * sin_co * sph.longitude.cos(),
            sph.radius * sin_co * sph.longitude.sin(),
            sph.radius * sph.colatitude.cos(),
        ])
    }
}

/// Cylindrical coordinates.
///
/// - `r`: Radial distance in xy-plane
/// - `longitude`: Angle in xy-plane from +X, range (-pi, pi]
/// - `z`: Height along z-axis
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cylindrical {
    pub r: f64,
    pub longitude: f64,
    pub z: f64,
}

impl Cylindrical {
    #[inline]
    pub fn new(r: f64, longitude: f64, z: f64) -> Self {
        Cylindrical { r, longitude, z }
    }

    #[inline]
    pub fn longitude_deg(&self) -> f64 {
        self.longitude.to_degrees()
    }
}

impl From<Cylindrical> for Rectangular {
    fn from(cyl: Cylindrical) -> Self {
        Rectangular([
            cyl.r * cyl.longitude.cos(),
            cyl.r * cyl.longitude.sin(),
            cyl.z,
        ])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI};

    const EPSILON: f64 = 1e-12;

    fn assert_near(a: f64, b: f64, msg: &str) {
        assert!(
            (a - b).abs() < EPSILON,
            "{}: {} != {} (diff={})",
            msg,
            a,
            b,
            (a - b).abs()
        );
    }

    #[test]
    fn test_rectangular_to_latitudinal() {
        let rect = Rectangular([6378.0, 0.0, 0.0]);
        let lat = rect.to_latitudinal();
        assert_near(lat.radius, 6378.0, "radius");
        assert_near(lat.longitude, 0.0, "longitude");
        assert_near(lat.latitude, 0.0, "latitude");

        let rect = Rectangular([0.0, 6378.0, 0.0]);
        let lat = rect.to_latitudinal();
        assert_near(lat.longitude, FRAC_PI_2, "longitude +Y");

        let rect = Rectangular([0.0, 0.0, 6378.0]);
        let lat = rect.to_latitudinal();
        assert_near(lat.latitude, FRAC_PI_2, "latitude +Z");
    }

    #[test]
    fn test_latitudinal_round_trip() {
        let original = Rectangular([1000.0, 2000.0, 3000.0]);
        let lat = original.to_latitudinal();
        let back: Rectangular = lat.into();
        assert_near(back.x(), original.x(), "x");
        assert_near(back.y(), original.y(), "y");
        assert_near(back.z(), original.z(), "z");
    }

    #[test]
    fn test_spherical_round_trip() {
        let original = Rectangular([1000.0, 2000.0, 3000.0]);
        let sph = original.to_spherical();
        let back: Rectangular = sph.into();
        assert_near(back.x(), original.x(), "x");
        assert_near(back.y(), original.y(), "y");
        assert_near(back.z(), original.z(), "z");
    }

    #[test]
    fn test_cylindrical_round_trip() {
        let original = Rectangular([1000.0, 2000.0, 3000.0]);
        let cyl = original.to_cylindrical();
        let back: Rectangular = cyl.into();
        assert_near(back.x(), original.x(), "x");
        assert_near(back.y(), original.y(), "y");
        assert_near(back.z(), original.z(), "z");
    }

    #[test]
    fn test_origin() {
        let origin = Rectangular([0.0, 0.0, 0.0]);
        assert_eq!(origin.to_latitudinal().radius, 0.0);
        assert_eq!(origin.to_spherical().radius, 0.0);
        assert_eq!(origin.to_cylindrical().r, 0.0);
    }

    #[test]
    fn test_near_pole_no_nan() {
        // When z ≈ radius, floating-point rounding could make z/radius slightly > 1.0
        // which would cause asin/acos to return NaN without clamping.
        let r = 6378.0_f64;
        // Construct a point where z/radius could exceed 1.0 due to rounding
        let rect = Rectangular([1e-15, 1e-15, r]);
        let lat = rect.to_latitudinal();
        assert!(!lat.latitude.is_nan(), "latitude should not be NaN near pole");
        assert!((lat.latitude - FRAC_PI_2).abs() < 1e-10, "should be near +90°");

        let sph = rect.to_spherical();
        assert!(!sph.colatitude.is_nan(), "colatitude should not be NaN near pole");
        assert!(sph.colatitude.abs() < 1e-10, "should be near 0° colatitude");
    }

    #[test]
    fn test_degree_conversions() {
        let lat = Latitudinal::new(1.0, PI / 4.0, PI / 6.0);
        assert_near(lat.longitude_deg(), 45.0, "lon deg");
        assert_near(lat.latitude_deg(), 30.0, "lat deg");

        let sph = Spherical::new(1.0, PI / 3.0, PI / 2.0);
        assert_near(sph.colatitude_deg(), 60.0, "colat deg");
        assert_near(sph.longitude_deg(), 90.0, "lon deg");

        let cyl = Cylindrical::new(1.0, -PI / 4.0, 0.0);
        assert_near(cyl.longitude_deg(), -45.0, "lon deg");
    }
}
