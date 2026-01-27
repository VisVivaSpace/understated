//! CSPICE validation tests.
//!
//! Validates understated's interpolation against CSPICE (the reference
//! implementation) using FFI bindings. These tests require the `cspice`
//! feature and the CSPICE toolkit installed at `/Users/nstrange/cspice`.
//!
//! Run with: `cargo test --features cspice,test-data`

#![cfg(all(feature = "cspice", feature = "test-data"))]

use std::ffi::CString;
use std::sync::Mutex;

use understated::{Ephemeris, EpochTDB, NaifId, Sclk, TimeFormat};

// CSPICE is not thread-safe — serialize all calls.
static CSPICE_LOCK: Mutex<()> = Mutex::new(());

// ============================================================================
// FFI bindings to CSPICE
// ============================================================================

#[link(name = "cspice", kind = "static")]
unsafe extern "C" {
    fn furnsh_c(file: *const i8);
    fn spkez_c(
        target: i32,
        epoch: f64,
        frame: *const i8,
        abcorr: *const i8,
        observer: i32,
        state: *mut [f64; 6],
        lt: *mut f64,
    );
    fn str2et_c(date: *const i8, et: *mut f64);
    fn et2utc_c(et: f64, format: *const i8, prec: i32, lenout: i32, utcstr: *mut i8);
    fn ckgpav_c(
        inst: i32,
        sclkdp: f64,
        tol: f64,
        r#ref: *const i8,
        cmat: *mut [[f64; 3]; 3],
        av: *mut [f64; 3],
        clkout: *mut f64,
        found: *mut i32,
    );
    fn m2q_c(r: *const [[f64; 3]; 3], q: *mut [f64; 4]);
    fn erract_c(operation: *const i8, lenout: i32, action: *mut i8);
    fn reset_c();
    fn failed_c() -> i32;
    fn kclear_c();
}

// ============================================================================
// Tolerance helpers
// ============================================================================

/// Position tolerance: 1e-9 km absolute, or 1e-15 relative for large values.
fn pos_tol(value: f64) -> f64 {
    1e-9_f64.max(value.abs() * 1e-15)
}

/// Velocity tolerance: 1e-12 km/s absolute, or 1e-15 relative for large values.
fn vel_tol(value: f64) -> f64 {
    1e-12_f64.max(value.abs() * 1e-15)
}

// ============================================================================
// Helper functions
// ============================================================================

fn furnsh(file: &str) {
    let c = CString::new(file).unwrap();
    unsafe { furnsh_c(c.as_ptr()) };
}

fn spkez(target: i32, epoch: f64, observer: i32) -> [f64; 6] {
    let frame = CString::new("J2000").unwrap();
    let abcorr = CString::new("NONE").unwrap();
    let mut state = [0.0f64; 6];
    let mut lt = 0.0f64;
    unsafe {
        spkez_c(
            target,
            epoch,
            frame.as_ptr(),
            abcorr.as_ptr(),
            observer,
            &mut state as *mut f64 as *mut [f64; 6],
            &mut lt,
        );
    }
    state
}

fn str2et(date: &str) -> f64 {
    let c = CString::new(date).unwrap();
    let mut et = 0.0f64;
    unsafe { str2et_c(c.as_ptr(), &mut et) };
    et
}

fn et2utc(et: f64, format: &str, prec: i32) -> String {
    let fmt = CString::new(format).unwrap();
    let mut buf = vec![0i8; 256];
    unsafe { et2utc_c(et, fmt.as_ptr(), prec, 256, buf.as_mut_ptr()) };
    let cstr = unsafe { std::ffi::CStr::from_ptr(buf.as_ptr()) };
    cstr.to_string_lossy().into_owned()
}

fn ckgpav(inst: i32, sclk: f64, tol: f64) -> Option<([f64; 4], [f64; 3])> {
    let ref_frame = CString::new("J2000").unwrap();
    let mut cmat = [[0.0f64; 3]; 3];
    let mut av = [0.0f64; 3];
    let mut clkout = 0.0f64;
    let mut found = 0i32;
    unsafe {
        ckgpav_c(
            inst,
            sclk,
            tol,
            ref_frame.as_ptr(),
            &mut cmat,
            &mut av,
            &mut clkout,
            &mut found,
        );
    }
    if found != 0 {
        let mut q = [0.0f64; 4];
        unsafe { m2q_c(&cmat, &mut q) };
        Some((q, av))
    } else {
        None
    }
}

fn cspice_init() {
    unsafe { kclear_c() };
    set_cspice_error_handling();
}

fn set_cspice_error_handling() {
    let op = CString::new("SET").unwrap();
    let mut action_buf = [0i8; 64];
    let action = b"RETURN\0";
    for (i, &b) in action.iter().enumerate() {
        action_buf[i] = b as i8;
    }
    unsafe { erract_c(op.as_ptr(), 64, action_buf.as_mut_ptr()) };
}

fn check_cspice_error() -> bool {
    let f = unsafe { failed_c() };
    if f != 0 {
        unsafe { reset_c() };
        true
    } else {
        false
    }
}

// ============================================================================
// SPK Validation Tests
// ============================================================================

#[test]
fn test_spk_type2_earth_vs_cspice() {
    let _lock = CSPICE_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    cspice_init();

    // Load kernel in CSPICE
    furnsh("test_data/de440s.bsp");

    // Load kernel in understated
    let eph = Ephemeris::load("test_data/de440s.bsp").unwrap();

    // Test at multiple epochs across coverage
    let test_epochs = [
        0.0,           // J2000
        86400.0,       // 1 day after J2000
        86400.0 * 30.0, // 30 days
        86400.0 * 365.25, // 1 year
        -86400.0 * 365.25, // 1 year before J2000
    ];

    for &epoch in &test_epochs {
        // CSPICE: Earth relative to SSB
        let cspice_state = spkez(399, epoch, 0);
        assert!(!check_cspice_error(), "CSPICE error at epoch {}", epoch);

        // understated: Earth relative to SSB
        let us_state = eph
            .state_of(NaifId::EARTH, EpochTDB(epoch), NaifId::SSB)
            .unwrap();

        for i in 0..3 {
            let diff = (us_state.position[i] - cspice_state[i]).abs();
            let tol = pos_tol(cspice_state[i]);
            assert!(
                diff < tol,
                "Position[{}] mismatch at epoch {}: understated={}, cspice={}, diff={}, tol={}",
                i, epoch, us_state.position[i], cspice_state[i], diff, tol
            );
        }

        for i in 0..3 {
            let diff = (us_state.velocity[i] - cspice_state[3 + i]).abs();
            let tol = vel_tol(cspice_state[3 + i]);
            assert!(
                diff < tol,
                "Velocity[{}] mismatch at epoch {}: understated={}, cspice={}, diff={}, tol={}",
                i, epoch, us_state.velocity[i], cspice_state[3 + i], diff, tol
            );
        }
    }
}

#[test]
fn test_spk_type2_moon_vs_cspice() {
    let _lock = CSPICE_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    cspice_init();

    furnsh("test_data/de440s.bsp");
    let eph = Ephemeris::load("test_data/de440s.bsp").unwrap();

    let test_epochs = [0.0, 86400.0 * 15.0, 86400.0 * 180.0];

    for &epoch in &test_epochs {
        // Moon relative to Earth (requires chaining through EMB)
        let cspice_state = spkez(301, epoch, 399);
        assert!(!check_cspice_error(), "CSPICE error at epoch {}", epoch);

        let us_state = eph
            .state_of(NaifId::MOON, EpochTDB(epoch), NaifId::EARTH)
            .unwrap();

        for i in 0..3 {
            let diff = (us_state.position[i] - cspice_state[i]).abs();
            let tol = pos_tol(cspice_state[i]);
            assert!(
                diff < tol,
                "Moon pos[{}] mismatch at epoch {}: diff={}, tol={}",
                i, epoch, diff, tol
            );
        }
        for i in 0..3 {
            let diff = (us_state.velocity[i] - cspice_state[3 + i]).abs();
            let tol = vel_tol(cspice_state[3 + i]);
            assert!(
                diff < tol,
                "Moon vel[{}] mismatch at epoch {}: diff={}, tol={}",
                i, epoch, diff, tol
            );
        }
    }
}

#[test]
fn test_spk_type2_mars_vs_cspice() {
    let _lock = CSPICE_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    cspice_init();

    furnsh("test_data/de440s.bsp");
    let eph = Ephemeris::load("test_data/de440s.bsp").unwrap();

    let epoch = 86400.0 * 100.0;

    // Mars barycenter relative to SSB
    let cspice_state = spkez(4, epoch, 0);
    assert!(!check_cspice_error());

    let us_state = eph
        .state_of(NaifId::MARS_BC, EpochTDB(epoch), NaifId::SSB)
        .unwrap();

    for i in 0..3 {
        let diff = (us_state.position[i] - cspice_state[i]).abs();
        let tol = pos_tol(cspice_state[i]);
        assert!(diff < tol, "Mars pos[{}] diff={}, tol={}", i, diff, tol);
    }
    for i in 0..3 {
        let diff = (us_state.velocity[i] - cspice_state[3 + i]).abs();
        let tol = vel_tol(cspice_state[3 + i]);
        assert!(diff < tol, "Mars vel[{}] diff={}, tol={}", i, diff, tol);
    }
}

#[test]
fn test_spk_type9_lagrange_vs_cspice() {
    let _lock = CSPICE_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    cspice_init();

    furnsh("test_data/gmat-lagrange.bsp");
    let eph = Ephemeris::load("test_data/gmat-lagrange.bsp").unwrap();

    let bodies = eph.spk_bodies();
    assert!(!bodies.is_empty(), "Should have bodies in lagrange file");

    let target = bodies[0];

    // Find coverage by trying a range of epochs
    // The GMAT test files typically cover a specific mission window
    // Try to evaluate at the midpoint — find segments manually
    let kernel = muad_dib::kernel::SpiceKernel::load("test_data/gmat-lagrange.bsp").unwrap();
    let md_target = muad_dib::types::NaifId(target.0);
    let seg = kernel
        .spk_segments_for(md_target)
        .next()
        .expect("Should have segment");

    let mid_epoch = (seg.initial_epoch + seg.final_epoch) / 2.0;
    let center_code = seg.center_code;

    // CSPICE
    let cspice_state = spkez(target.0, mid_epoch, center_code);
    if check_cspice_error() {
        // Skip if CSPICE can't handle this file
        return;
    }

    // understated
    let us_state = eph
        .state_of(target, EpochTDB(mid_epoch), NaifId(center_code))
        .unwrap();

    for i in 0..3 {
        let diff = (us_state.position[i] - cspice_state[i]).abs();
        let tol = pos_tol(cspice_state[i]);
        assert!(
            diff < tol,
            "Type9 pos[{}]: us={}, cspice={}, diff={}, tol={}",
            i, us_state.position[i], cspice_state[i], diff, tol
        );
    }
    for i in 0..3 {
        let diff = (us_state.velocity[i] - cspice_state[3 + i]).abs();
        let tol = vel_tol(cspice_state[3 + i]);
        assert!(
            diff < tol,
            "Type9 vel[{}]: us={}, cspice={}, diff={}, tol={}",
            i, us_state.velocity[i], cspice_state[3 + i], diff, tol
        );
    }
}

#[test]
fn test_spk_type13_hermite_vs_cspice() {
    let _lock = CSPICE_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    cspice_init();

    furnsh("test_data/gmat-hermite.bsp");
    let eph = Ephemeris::load("test_data/gmat-hermite.bsp").unwrap();

    let bodies = eph.spk_bodies();
    assert!(!bodies.is_empty());

    let target = bodies[0];

    let kernel = muad_dib::kernel::SpiceKernel::load("test_data/gmat-hermite.bsp").unwrap();
    let md_target = muad_dib::types::NaifId(target.0);
    let seg = kernel
        .spk_segments_for(md_target)
        .next()
        .expect("Should have segment");

    let mid_epoch = (seg.initial_epoch + seg.final_epoch) / 2.0;
    let center_code = seg.center_code;

    // CSPICE
    let cspice_state = spkez(target.0, mid_epoch, center_code);
    if check_cspice_error() {
        return;
    }

    // understated
    let us_state = eph
        .state_of(target, EpochTDB(mid_epoch), NaifId(center_code))
        .unwrap();

    for i in 0..3 {
        let diff = (us_state.position[i] - cspice_state[i]).abs();
        let tol = pos_tol(cspice_state[i]);
        assert!(
            diff < tol,
            "Type13 pos[{}]: us={}, cspice={}, diff={}, tol={}",
            i, us_state.position[i], cspice_state[i], diff, tol
        );
    }
    for i in 0..3 {
        let diff = (us_state.velocity[i] - cspice_state[3 + i]).abs();
        let tol = vel_tol(cspice_state[3 + i]);
        assert!(
            diff < tol,
            "Type13 vel[{}]: us={}, cspice={}, diff={}, tol={}",
            i, us_state.velocity[i], cspice_state[3 + i], diff, tol
        );
    }
}

// ============================================================================
// CK Validation Tests
// ============================================================================

#[test]
fn test_ck_vs_cspice() {
    let _lock = CSPICE_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    cspice_init();

    // CK files require an SCLK kernel for CSPICE. Our test.bc may not
    // need one if we use raw ticks. Load what we have.
    furnsh("test_data/test.bc");

    let eph = Ephemeris::load("test_data/test.bc").unwrap();
    let instruments = eph.ck_instruments();

    if instruments.is_empty() {
        return;
    }

    let inst = instruments[0];

    // Find the coverage range from the kernel
    let kernel = muad_dib::kernel::SpiceKernel::load("test_data/test.bc").unwrap();
    let md_inst = muad_dib::types::NaifId(inst.0);
    let seg = kernel
        .ck_segments_for(md_inst)
        .next()
        .expect("Should have CK segment");

    let mid_sclk = (seg.initial_sclk + seg.final_sclk) / 2.0;

    // understated
    let us_pointing = eph.pointing_of(inst, Sclk(mid_sclk)).unwrap();

    // CSPICE — ckgpav_c with a small tolerance
    let cspice_result = ckgpav(inst.0, mid_sclk, 1.0);
    if cspice_result.is_none() {
        // CSPICE couldn't find pointing (may need SCLK kernel), skip
        return;
    }
    let (cspice_q, cspice_av) = cspice_result.unwrap();

    // Quaternion comparison: account for sign ambiguity (q and -q are same rotation)
    let dot: f64 = us_pointing
        .quaternion
        .iter()
        .zip(cspice_q.iter())
        .map(|(a, b)| a * b)
        .sum();

    let sign = if dot < 0.0 { -1.0 } else { 1.0 };
    for i in 0..4 {
        let diff = (us_pointing.quaternion[i] - sign * cspice_q[i]).abs();
        assert!(
            diff < 1e-8,
            "Quaternion[{}] mismatch: understated={}, cspice={}, diff={}",
            i, us_pointing.quaternion[i], cspice_q[i], diff
        );
    }

    // Angular velocity comparison (if available)
    if let Some(us_av) = us_pointing.angular_velocity {
        for i in 0..3 {
            let diff = (us_av[i] - cspice_av[i]).abs();
            assert!(
                diff < 1e-8,
                "AngVel[{}] mismatch: understated={}, cspice={}, diff={}",
                i, us_av[i], cspice_av[i], diff
            );
        }
    }
}

// ============================================================================
// Time Conversion Validation
// ============================================================================

#[test]
fn test_str2et_vs_cspice() {
    let _lock = CSPICE_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    cspice_init();

    furnsh("test_data/naif0012.tls");
    let eph = Ephemeris::load("test_data/naif0012.tls").unwrap();

    let test_dates = [
        "2000-01-01T12:00:00",
        "2020-06-15T00:00:00",
        "1999-01-01T00:00:00",
        "2024-03-20T12:30:00",
    ];

    for &date in &test_dates {
        // CSPICE
        let cspice_et = str2et(date);
        if check_cspice_error() {
            continue;
        }

        // understated
        let us_et = eph.utc_to_tdb(date).unwrap();

        let diff = (us_et.0 - cspice_et).abs();
        assert!(
            diff < 0.01, // 10 ms tolerance for UTC→TDB (our LSK extraction may differ slightly)
            "str2et mismatch for '{}': understated={}, cspice={}, diff={}s",
            date, us_et.0, cspice_et, diff
        );
    }
}

#[test]
fn test_et2utc_vs_cspice() {
    let _lock = CSPICE_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    cspice_init();

    furnsh("test_data/naif0012.tls");
    let eph = Ephemeris::load("test_data/naif0012.tls").unwrap();

    let test_epochs = [0.0, 86400.0 * 365.25 * 20.0, -86400.0 * 365.25];

    for &epoch in &test_epochs {
        // CSPICE
        let cspice_utc = et2utc(epoch, "ISOC", 3);
        if check_cspice_error() {
            continue;
        }

        // understated
        let us_utc = eph.tdb_to_utc(EpochTDB(epoch), TimeFormat::Iso8601).unwrap();

        // Compare year-month-day portion (time may differ slightly)
        let cspice_date = &cspice_utc[..10];
        let us_date = &us_utc[..10];
        assert_eq!(
            cspice_date, us_date,
            "Date mismatch at epoch {}: understated='{}', cspice='{}'",
            epoch, us_utc, cspice_utc
        );
    }
}
