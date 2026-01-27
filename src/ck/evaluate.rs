//! CK Type 1 (discrete) and Type 3 (linear interpolation) evaluation.
//!
//! Type 3 evaluation matches CSPICE's CKE03 algorithm: rotation matrix
//! interpolation via axis-angle decomposition, not quaternion SLERP.

use crate::error::{Error, Result};
use crate::pointing::Pointing;
use crate::rotation::{axisar, m2q, mtxm, mxmt, q2m};
use muad_dib::kernel::ck_types::{Ck1Data, Ck3Data, PointingRecord};

/// Linear interpolation of angular velocity between two records.
///
/// Matches CSPICE CKE03: `AV = (1 - frac) * AV1 + frac * AV2`.
fn interpolate_angular_velocity(
    rec0: &PointingRecord,
    rec1: &PointingRecord,
    frac: f64,
) -> Option<[f64; 3]> {
    let av0 = rec0.angular_velocity?;
    let av1 = rec1.angular_velocity?;
    let w0 = 1.0 - frac;

    Some([
        w0 * av0[0] + frac * av1[0],
        w0 * av0[1] + frac * av1[1],
        w0 * av0[2] + frac * av1[2],
    ])
}

/// Evaluate CK Type 1: discrete pointing (most recent record at or before query).
pub fn evaluate_type1(data: &Ck1Data, sclk: f64) -> Result<Pointing> {
    if data.records.is_empty() {
        return Err(Error::EpochOutOfRange {
            epoch: sclk,
            start: 0.0,
            end: 0.0,
        });
    }

    let mut best_idx = None;
    for (i, record) in data.records.iter().enumerate() {
        if record.sclk <= sclk {
            best_idx = Some(i);
        } else {
            break;
        }
    }

    let idx = best_idx.ok_or(Error::EpochOutOfRange {
        epoch: sclk,
        start: data.records.first().unwrap().sclk,
        end: data.records.last().unwrap().sclk,
    })?;

    let record = &data.records[idx];
    Ok(Pointing::new_raw(
        record.quaternion(),
        record.angular_velocity,
    ))
}

/// Evaluate CK Type 3: linear interpolation between bracketing records.
///
/// Matches CSPICE's CKE03 algorithm:
/// 1. Convert bracketing quaternions to rotation matrices
/// 2. Compute relative rotation: `ROT = CMAT2^T * CMAT1`
/// 3. Extract rotation axis and angle
/// 4. Build partial rotation: `DELTA = axisar(axis, angle * frac)`
/// 5. Compose: `CMAT = CMAT1 * DELTA^T`
/// 6. Convert result back to quaternion
pub fn evaluate_type3(data: &Ck3Data, sclk: f64) -> Result<Pointing> {
    if data.records.is_empty() {
        return Err(Error::EpochOutOfRange {
            epoch: sclk,
            start: 0.0,
            end: 0.0,
        });
    }

    let start_sclk = data.records.first().unwrap().sclk;
    let end_sclk = data.records.last().unwrap().sclk;

    if sclk < start_sclk || sclk > end_sclk {
        return Err(Error::EpochOutOfRange {
            epoch: sclk,
            start: start_sclk,
            end: end_sclk,
        });
    }

    // Find lower bracketing record
    let mut lower_idx = 0;
    for (i, record) in data.records.iter().enumerate() {
        if record.sclk <= sclk {
            lower_idx = i;
        } else {
            break;
        }
    }

    let rec0 = &data.records[lower_idx];

    // CKE03: when t1 == t2 (same record or exact match), just convert and return.
    let upper_idx = (lower_idx + 1).min(data.records.len() - 1);
    let rec1 = &data.records[upper_idx];

    if rec0.sclk == rec1.sclk || upper_idx == lower_idx {
        // Single record — convert quaternion directly via q2m → m2q
        // to match the CSPICE path (which always goes through q2m).
        let cmat = q2m(&rec0.quaternion());
        let q = m2q(&cmat);
        return Ok(Pointing::new_raw(q, rec0.angular_velocity));
    }

    // Interpolation parameter
    let frac = (sclk - rec0.sclk) / (rec1.sclk - rec0.sclk);

    // Step 1: Convert quaternions to C-matrices
    let cmat1 = q2m(&rec0.quaternion());
    let cmat2 = q2m(&rec1.quaternion());

    // Step 2: Relative rotation ROT = CMAT2^T * CMAT1
    let rot = mtxm(&cmat2, &cmat1);

    // Step 3: Extract axis and angle
    let (axis, angle) = crate::rotation::raxisa(&rot);

    // Step 4: Partial rotation DELTA = axisar(axis, angle * frac)
    let delta = axisar(&axis, angle * frac);

    // Step 5: Compose CMAT = CMAT1 * DELTA^T
    let cmat = mxmt(&cmat1, &delta);

    // Step 6: Convert back to quaternion
    let q = m2q(&cmat);

    // Angular velocity: weighted average (same as CSPICE CKE03)
    let av = interpolate_angular_velocity(rec0, rec1, frac);

    Ok(Pointing::new_raw(q, av))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_1_SQRT_2;

    fn make_record(sclk: f64, q: [f64; 4], av: Option<[f64; 3]>) -> PointingRecord {
        PointingRecord {
            sclk,
            q0: q[0],
            q1: q[1],
            q2: q[2],
            q3: q[3],
            angular_velocity: av,
        }
    }

    #[test]
    fn test_evaluate_type1() {
        let data = Ck1Data {
            has_rates: true,
            records: vec![
                make_record(100.0, [1.0, 0.0, 0.0, 0.0], Some([0.0, 0.0, 0.1])),
                make_record(
                    200.0,
                    [FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0, 0.0],
                    Some([0.0, 0.0, 0.2]),
                ),
            ],
        };

        // At first record
        let p = evaluate_type1(&data, 100.0).unwrap();
        assert!((p.quaternion[0] - 1.0).abs() < 1e-10);

        // Between records - returns most recent
        let p = evaluate_type1(&data, 150.0).unwrap();
        assert!((p.quaternion[0] - 1.0).abs() < 1e-10);

        // At second record
        let p = evaluate_type1(&data, 200.0).unwrap();
        assert!((p.quaternion[0] - FRAC_1_SQRT_2).abs() < 1e-14);
    }

    #[test]
    fn test_evaluate_type3_exact_match() {
        let data = Ck3Data {
            has_rates: false,
            records: vec![
                make_record(100.0, [1.0, 0.0, 0.0, 0.0], None),
                make_record(200.0, [0.0, 1.0, 0.0, 0.0], None),
            ],
            interval_starts: vec![0],
        };

        // At first record: q2m → m2q roundtrip preserves identity quaternion
        let p = evaluate_type3(&data, 100.0).unwrap();
        assert!((p.quaternion[0] - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_evaluate_type3_midpoint() {
        // Identity and 180° rotation about x-axis
        let data = Ck3Data {
            has_rates: false,
            records: vec![
                make_record(100.0, [1.0, 0.0, 0.0, 0.0], None),
                make_record(200.0, [0.0, 1.0, 0.0, 0.0], None),
            ],
            interval_starts: vec![0],
        };

        // At midpoint: should be a 90° rotation about x-axis
        let p = evaluate_type3(&data, 150.0).unwrap();
        // Expected quaternion for 90° about x: (cos(45°), sin(45°), 0, 0)
        let expected_q0 = FRAC_1_SQRT_2;
        let expected_q1 = FRAC_1_SQRT_2;
        assert!(
            (p.quaternion[0] - expected_q0).abs() < 1e-14,
            "q0: {} vs {}",
            p.quaternion[0],
            expected_q0
        );
        assert!(
            (p.quaternion[1] - expected_q1).abs() < 1e-14,
            "q1: {} vs {}",
            p.quaternion[1],
            expected_q1
        );
    }

    #[test]
    fn test_evaluate_type3_out_of_range() {
        let data = Ck3Data {
            has_rates: false,
            records: vec![make_record(100.0, [1.0, 0.0, 0.0, 0.0], None)],
            interval_starts: vec![0],
        };
        assert!(matches!(
            evaluate_type3(&data, 50.0),
            Err(Error::EpochOutOfRange { .. })
        ));
    }
}
