//! CK Type 1 (discrete) and Type 3 (SLERP interpolated) evaluation.

use crate::ck::slerp::slerp;
use crate::error::{Error, Result};
use crate::pointing::Pointing;
use muad_dib::kernel::ck_types::{Ck1Data, Ck3Data, PointingRecord};

/// Linear interpolation of angular velocity between two records.
fn interpolate_angular_velocity(
    rec0: &PointingRecord,
    rec1: &PointingRecord,
    t: f64,
) -> Option<[f64; 3]> {
    let av0 = rec0.angular_velocity()?;
    let av1 = rec1.angular_velocity()?;

    Some([
        av0[0] + t * (av1[0] - av0[0]),
        av0[1] + t * (av1[1] - av0[1]),
        av0[2] + t * (av1[2] - av0[2]),
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
        record.angular_velocity(),
    ))
}

/// Evaluate CK Type 3: SLERP interpolation between bracketing records.
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

    // Exact match
    if (data.records[lower_idx].sclk - sclk).abs() < 1e-10 {
        let record = &data.records[lower_idx];
        return Ok(Pointing::new_raw(
            record.quaternion(),
            record.angular_velocity(),
        ));
    }

    let upper_idx = (lower_idx + 1).min(data.records.len() - 1);
    if upper_idx == lower_idx {
        let record = &data.records[lower_idx];
        return Ok(Pointing::new_raw(
            record.quaternion(),
            record.angular_velocity(),
        ));
    }

    let rec0 = &data.records[lower_idx];
    let rec1 = &data.records[upper_idx];

    let dt = rec1.sclk - rec0.sclk;
    let t = if dt.abs() < 1e-15 {
        0.0
    } else {
        (sclk - rec0.sclk) / dt
    };

    let q = slerp(&rec0.quaternion(), &rec1.quaternion(), t);
    let av = interpolate_angular_velocity(rec0, rec1, t);

    Ok(Pointing::new_raw(q, av))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_evaluate_type1() {
        let data = Ck1Data {
            has_rates: true,
            records: vec![
                PointingRecord {
                    sclk: 100.0,
                    q0: 1.0,
                    q1: 0.0,
                    q2: 0.0,
                    q3: 0.0,
                    av_x: Some(0.0),
                    av_y: Some(0.0),
                    av_z: Some(0.1),
                },
                PointingRecord {
                    sclk: 200.0,
                    q0: 0.707,
                    q1: 0.707,
                    q2: 0.0,
                    q3: 0.0,
                    av_x: Some(0.0),
                    av_y: Some(0.0),
                    av_z: Some(0.2),
                },
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
        assert!((p.quaternion[0] - 0.707).abs() < 1e-6);
    }

    #[test]
    fn test_evaluate_type3() {
        let data = Ck3Data {
            has_rates: false,
            records: vec![
                PointingRecord {
                    sclk: 100.0,
                    q0: 1.0,
                    q1: 0.0,
                    q2: 0.0,
                    q3: 0.0,
                    av_x: None,
                    av_y: None,
                    av_z: None,
                },
                PointingRecord {
                    sclk: 200.0,
                    q0: 0.0,
                    q1: 1.0,
                    q2: 0.0,
                    q3: 0.0,
                    av_x: None,
                    av_y: None,
                    av_z: None,
                },
            ],
            interval_starts: vec![0],
        };

        // At first record
        let p = evaluate_type3(&data, 100.0).unwrap();
        assert!((p.quaternion[0] - 1.0).abs() < 1e-10);

        // At midpoint - SLERP
        let p = evaluate_type3(&data, 150.0).unwrap();
        let expected = (PI / 4.0).cos();
        assert!((p.quaternion[0] - expected).abs() < 1e-6);
    }

    #[test]
    fn test_evaluate_type3_out_of_range() {
        let data = Ck3Data {
            has_rates: false,
            records: vec![PointingRecord {
                sclk: 100.0,
                q0: 1.0,
                q1: 0.0,
                q2: 0.0,
                q3: 0.0,
                av_x: None,
                av_y: None,
                av_z: None,
            }],
            interval_starts: vec![0],
        };
        assert!(matches!(evaluate_type3(&data, 50.0), Err(Error::EpochOutOfRange { .. })));
    }
}
