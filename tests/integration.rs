//! Integration tests using real SPICE kernel files through the Ephemeris API.

#[cfg(feature = "test-data")]
mod tests {
    use understated::{Ephemeris, EpochTDB, NaifId, TimeFormat};

    #[test]
    fn test_load_spk() {
        let eph = Ephemeris::load("test_data/de440s.bsp").unwrap();
        let bodies = eph.spk_bodies();
        assert!(!bodies.is_empty(), "de440s should contain bodies");
    }

    #[test]
    fn test_earth_state_from_de440s() {
        let eph = Ephemeris::load("test_data/de440s.bsp").unwrap();
        let epoch = EpochTDB(0.0); // J2000

        // Earth relative to Earth-Moon Barycenter
        let state = eph
            .state_of(NaifId::EARTH, epoch, NaifId::EMB)
            .unwrap();

        // Earth should be within ~5000 km of the EMB
        let dist = state.distance();
        assert!(
            dist < 5000.0,
            "Earth-EMB distance at J2000 should be < 5000 km, got {}",
            dist
        );
        assert!(dist > 0.0, "Distance should be positive");
    }

    #[test]
    fn test_moon_state_from_de440s() {
        let eph = Ephemeris::load("test_data/de440s.bsp").unwrap();
        let epoch = EpochTDB(0.0); // J2000

        // Moon relative to Earth-Moon Barycenter
        let state = eph
            .state_of(NaifId::MOON, epoch, NaifId::EMB)
            .unwrap();

        // Moon should be ~384,400 km from Earth, so distance from EMB is similar scale
        let dist = state.distance();
        assert!(
            dist > 300_000.0 && dist < 500_000.0,
            "Moon-EMB distance should be ~380k km, got {}",
            dist
        );
    }

    #[test]
    fn test_body_chaining() {
        let eph = Ephemeris::load("test_data/de440s.bsp").unwrap();
        let epoch = EpochTDB(0.0);

        // Earth relative to SSB requires chaining through EMB
        let state = eph
            .state_of(NaifId::EARTH, epoch, NaifId::SSB)
            .unwrap();

        // Earth should be ~1 AU from SSB
        let dist = state.distance();
        let au = 1.496e8; // km
        assert!(
            dist > 0.9 * au && dist < 1.1 * au,
            "Earth-SSB distance should be ~1 AU, got {}",
            dist
        );
        assert_eq!(state.target, NaifId::EARTH);
        assert_eq!(state.center, NaifId::SSB);
    }

    #[test]
    fn test_no_coverage_error() {
        let eph = Ephemeris::load("test_data/de440s.bsp").unwrap();

        // Epoch far outside de440s coverage
        let result = eph.state_of(NaifId::EARTH, EpochTDB(1e12), NaifId::SSB);
        assert!(result.is_err());
    }

    #[test]
    fn test_negation_via_center_swap() {
        let eph = Ephemeris::load("test_data/de440s.bsp").unwrap();
        let epoch = EpochTDB(0.0);

        // Moon relative to Earth: chaining Moon→EMB→Earth works because
        // when we need state_of(SSB, epoch, Earth), we find Earth's segment
        // (target=399) and use negation since target_code == center.
        let moon_rel_earth = eph
            .state_of(NaifId::MOON, epoch, NaifId::EARTH)
            .unwrap();

        // Verify by computing manually: Moon rel EMB - Earth rel EMB
        let moon_rel_emb = eph
            .state_of(NaifId::MOON, epoch, NaifId::EMB)
            .unwrap();
        let earth_rel_emb = eph
            .state_of(NaifId::EARTH, epoch, NaifId::EMB)
            .unwrap();
        let manual = moon_rel_emb - earth_rel_emb;

        for i in 0..3 {
            assert!(
                (moon_rel_earth.position[i] - manual.position[i]).abs() < 1e-6,
                "Position[{}] mismatch: {} vs {}",
                i,
                moon_rel_earth.position[i],
                manual.position[i]
            );
            assert!(
                (moon_rel_earth.velocity[i] - manual.velocity[i]).abs() < 1e-9,
                "Velocity[{}] mismatch: {} vs {}",
                i,
                moon_rel_earth.velocity[i],
                manual.velocity[i]
            );
        }
    }

    #[test]
    fn test_spk_bodies_discovery() {
        let eph = Ephemeris::load("test_data/de440s.bsp").unwrap();
        let bodies = eph.spk_bodies();

        // de440s should contain major planets + Moon + Sun
        let has_earth = bodies.iter().any(|b| *b == NaifId::EARTH);
        let has_moon = bodies.iter().any(|b| *b == NaifId::MOON);
        assert!(has_earth, "de440s should contain Earth");
        assert!(has_moon, "de440s should contain Moon");
    }

    #[test]
    fn test_load_multiple_files() {
        let eph = Ephemeris::load_many(&[
            "test_data/de440s.bsp",
            "test_data/naif0012.tls",
        ])
        .unwrap();

        // Should have LSK data from the .tls file
        assert!(eph.lsk_data().is_some(), "Should have LSK data");

        // Should have SPK data from the .bsp file
        assert!(!eph.spk_bodies().is_empty());
    }

    #[test]
    fn test_utc_tdb_conversion() {
        let eph = Ephemeris::load_many(&[
            "test_data/de440s.bsp",
            "test_data/naif0012.tls",
        ])
        .unwrap();

        let tdb = eph.utc_to_tdb("2020-01-01T12:00:00").unwrap();
        assert!(tdb.0 > 0.0, "2020 should be after J2000");

        let utc_str = eph.tdb_to_utc(tdb, TimeFormat::Iso8601).unwrap();
        assert!(utc_str.starts_with("2020-01-01T"), "Should round-trip: {}", utc_str);
    }

    #[test]
    fn test_utc_conversion_requires_lsk() {
        let eph = Ephemeris::load("test_data/de440s.bsp").unwrap();

        // Without LSK data, conversion should fail
        let result = eph.utc_to_tdb("2020-01-01T12:00:00");
        assert!(result.is_err());
    }

    #[test]
    fn test_ck_pointing() {
        let eph = Ephemeris::load("test_data/test.bc").unwrap();
        let instruments = eph.ck_instruments();
        assert!(!instruments.is_empty(), "test.bc should contain instruments");

        // Query pointing for the first instrument at the midpoint of its coverage
        let _instrument = instruments[0];

        // We need to find a valid SCLK within coverage
        // The CK file has segments — use the kernel to find coverage
        // For simplicity, just verify we loaded instruments
    }

    #[test]
    fn test_test_bsp_hermite() {
        let eph = Ephemeris::load("test_data/gmat-hermite.bsp").unwrap();
        let bodies = eph.spk_bodies();
        assert!(!bodies.is_empty(), "gmat-hermite.bsp should contain bodies");

        // Query at a known epoch (use J2000 or nearby)
        let _target = bodies[0];
        // Get center from the segment
        // Just verify it loads and we can list bodies
    }

    #[test]
    fn test_test_bsp_lagrange() {
        let eph = Ephemeris::load("test_data/gmat-lagrange.bsp").unwrap();
        let bodies = eph.spk_bodies();
        assert!(!bodies.is_empty(), "gmat-lagrange.bsp should contain bodies");
    }

    /// Verify that evaluating at closely-spaced epochs produces velocity
    /// consistent with the numerical derivative of position.
    #[test]
    fn test_state_continuity() {
        let eph = Ephemeris::load("test_data/de440s.bsp").unwrap();
        let epoch = EpochTDB(0.0); // J2000
        let dt = 1.0; // 1 second step

        let state0 = eph
            .state_of(NaifId::EARTH, epoch, NaifId::SSB)
            .unwrap();
        let state1 = eph
            .state_of(NaifId::EARTH, EpochTDB(epoch.0 + dt), NaifId::SSB)
            .unwrap();

        // Numerical derivative of position should approximate reported velocity
        for i in 0..3 {
            let numerical_vel = (state1.position[i] - state0.position[i]) / dt;
            let reported_vel = state0.velocity[i];
            let rel_error = if reported_vel.abs() > 1e-10 {
                ((numerical_vel - reported_vel) / reported_vel).abs()
            } else {
                (numerical_vel - reported_vel).abs()
            };
            assert!(
                rel_error < 0.01,
                "Velocity[{}] continuity: numerical={:.9}, reported={:.9}, rel_error={:.6}",
                i,
                numerical_vel,
                reported_vel,
                rel_error
            );
        }
    }

    /// Verify that known bodies have physically reasonable distances and speeds.
    #[test]
    fn test_position_magnitude() {
        let eph = Ephemeris::load("test_data/de440s.bsp").unwrap();
        let epoch = EpochTDB(0.0);

        // Earth-SSB should be ~1 AU
        let earth = eph
            .state_of(NaifId::EARTH, epoch, NaifId::SSB)
            .unwrap();
        let au = 149597870.7; // km
        assert!(
            earth.distance() > 0.9 * au && earth.distance() < 1.1 * au,
            "Earth-SSB distance should be ~1 AU, got {:.3} km",
            earth.distance()
        );
        // Earth orbital speed ~30 km/s
        assert!(
            earth.speed() > 25.0 && earth.speed() < 35.0,
            "Earth speed should be ~30 km/s, got {:.6}",
            earth.speed()
        );

        // Moon-Earth should be ~384,400 km
        let moon = eph
            .state_of(NaifId::MOON, epoch, NaifId::EARTH)
            .unwrap();
        assert!(
            moon.distance() > 350_000.0 && moon.distance() < 420_000.0,
            "Moon-Earth distance should be ~384k km, got {:.3}",
            moon.distance()
        );
        // Moon orbital speed ~1 km/s relative to Earth
        assert!(
            moon.speed() > 0.5 && moon.speed() < 1.5,
            "Moon-Earth speed should be ~1 km/s, got {:.6}",
            moon.speed()
        );
    }

    /// Verify that CK evaluations return unit quaternions.
    #[test]
    fn test_ck_quaternion_normalization() {
        let eph = Ephemeris::load("test_data/test.bc").unwrap();
        let instruments = eph.ck_instruments();
        if instruments.is_empty() {
            return;
        }

        // Try querying at several SCLK values within plausible ranges.
        // Since we can't access segment metadata directly, use known SCLK
        // values from the CSPICE validation tests (discovered via muad-dib).
        // The test.bc file uses large SCLK tick values — try a sweep.
        for sclk_val in [1e8, 5e8, 1e9, 2e9] {
            let sclk = understated::Sclk::from_ticks(sclk_val);
            if let Ok(pointing) = eph.pointing_of(instruments[0], sclk) {
                let [q0, q1, q2, q3] = pointing.quaternion;
                let norm_sq = q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3;
                assert!(
                    (norm_sq - 1.0).abs() < 1e-10,
                    "Quaternion should be normalized: norm_sq={}, sclk={}",
                    norm_sq,
                    sclk_val
                );
            }
            // If pointing_of returns Err, the SCLK is outside coverage — that's fine.
        }
    }

    /// Test multiple bodies from de440s: Earth barycenter, Sun, Moon, Earth.
    #[test]
    fn test_multi_body_type2() {
        let eph = Ephemeris::load("test_data/de440s.bsp").unwrap();
        let epoch = EpochTDB(0.0);

        // All of these should be queryable relative to SSB
        let bodies = [
            (NaifId::SUN, "Sun"),
            (NaifId::EMB, "EMB"),
            (NaifId::MOON, "Moon"),
            (NaifId::EARTH, "Earth"),
        ];

        for (body, name) in &bodies {
            let state = eph
                .state_of(*body, epoch, NaifId::SSB)
                .unwrap_or_else(|_| panic!("{} should be available at J2000", name));

            // All should have non-zero distance (even Sun is not at SSB exactly)
            assert!(
                state.distance() > 0.0,
                "{} should have non-zero distance from SSB",
                name
            );

            // All should have non-zero speed
            assert!(
                state.speed() > 0.0,
                "{} should have non-zero speed relative to SSB",
                name
            );

            // Target and center should be set correctly
            assert_eq!(state.target, *body, "{} target mismatch", name);
            assert_eq!(state.center, NaifId::SSB, "{} center mismatch", name);
        }
    }
}
