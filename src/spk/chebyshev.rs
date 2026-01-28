//! Chebyshev polynomial interpolation for SPK Types 2 and 3.
//!
//! Uses Clenshaw recurrence for numerical stability.

use crate::error::{Error, Result};
use crate::state::State;
use muad_dib::kernel::spk_types::{ChebyshevRecord, ChebyshevRecordWithVelocity, Spk2Data, Spk3Data};

/// Clenshaw recurrence for Chebyshev polynomial evaluation.
///
/// Given coefficients [c0, c1, ..., cn] and normalized argument s in [-1, 1],
/// computes c0*T0(s) + c1*T1(s) + ... + cn*Tn(s).
fn clenshaw(coeffs: &[f64], s: f64) -> f64 {
    let n = coeffs.len();
    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return coeffs[0];
    }

    let s2 = 2.0 * s;
    let mut b_k = 0.0;
    let mut b_k1 = 0.0;

    for i in (1..n).rev() {
        let b_k2 = b_k1;
        b_k1 = b_k;
        b_k = s2 * b_k1 - b_k2 + coeffs[i];
    }

    coeffs[0] + s * b_k - b_k1
}

/// Chebyshev derivative via T-polynomial coefficient method.
///
/// Computes derivative coefficients g[] via recurrence, then evaluates
/// with standard Clenshaw. Reference: Numerical Recipes Section 5.9.
fn clenshaw_derivative(coeffs: &[f64], s: f64, scale: f64) -> f64 {
    let n = coeffs.len();
    if n <= 1 {
        return 0.0;
    }

    let mut g = vec![0.0; n - 1];
    g[n - 2] = 2.0 * (n - 1) as f64 * coeffs[n - 1];

    for k in (0..n - 2).rev() {
        let g_k_plus_2 = if k + 2 < g.len() { g[k + 2] } else { 0.0 };
        g[k] = g_k_plus_2 + 2.0 * (k + 1) as f64 * coeffs[k + 1];
    }

    if !g.is_empty() {
        g[0] /= 2.0;
    }

    scale * clenshaw(&g, s)
}

/// Trait abstracting Chebyshev record access for Types 2 and 3.
trait ChebyshevRecordAccess {
    fn midpoint(&self) -> f64;
    fn radius(&self) -> f64;
}

impl ChebyshevRecordAccess for ChebyshevRecord {
    fn midpoint(&self) -> f64 { self.midpoint }
    fn radius(&self) -> f64 { self.radius }
}

impl ChebyshevRecordAccess for ChebyshevRecordWithVelocity {
    fn midpoint(&self) -> f64 { self.midpoint }
    fn radius(&self) -> f64 { self.radius }
}

/// Find the Chebyshev record containing `epoch` using binary search.
///
/// Records are sorted by midpoint. Returns the record index and
/// the normalized argument `s` in [-1, 1].
fn find_record_generic<R: ChebyshevRecordAccess>(
    records: &[R],
    init_epoch: f64,
    epoch: f64,
) -> Result<(usize, f64)> {
    if records.is_empty() {
        return Err(Error::EpochOutOfRange {
            epoch,
            start: init_epoch,
            end: init_epoch,
        });
    }

    // Binary search: find the record whose interval contains the epoch.
    // Records are sorted by midpoint and contiguous, so we search by midpoint.
    let n = records.len();
    let mut lo = 0usize;
    let mut hi = n;
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if records[mid].midpoint() + records[mid].radius() < epoch {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }

    // Check the candidate record (and its neighbor) for containment.
    let check_start = lo.saturating_sub(1);
    let check_end = n.min(lo + 2);
    for (idx, r) in records.iter().enumerate().take(check_end).skip(check_start) {
        let start = r.midpoint() - r.radius();
        let end = r.midpoint() + r.radius();
        if epoch >= start && epoch <= end {
            let s = (epoch - r.midpoint()) / r.radius();
            return Ok((idx, s));
        }
    }

    let first = &records[0];
    let last = &records[n - 1];
    Err(Error::EpochOutOfRange {
        epoch,
        start: first.midpoint() - first.radius(),
        end: last.midpoint() + last.radius(),
    })
}

/// Evaluate SPK Type 2 (Chebyshev position only, velocity via differentiation).
pub fn evaluate_type2(data: &Spk2Data, epoch: f64) -> Result<State> {
    let (idx, s) = find_record_generic(&data.records, data.init_epoch, epoch)?;
    let record = &data.records[idx];

    let x = clenshaw(&record.x_coeffs, s);
    let y = clenshaw(&record.y_coeffs, s);
    let z = clenshaw(&record.z_coeffs, s);

    let scale = 1.0 / record.radius;
    let vx = clenshaw_derivative(&record.x_coeffs, s, scale);
    let vy = clenshaw_derivative(&record.y_coeffs, s, scale);
    let vz = clenshaw_derivative(&record.z_coeffs, s, scale);

    Ok(State::new_raw([x, y, z], [vx, vy, vz]))
}

/// Evaluate SPK Type 3 (Chebyshev position + velocity coefficients).
pub fn evaluate_type3(data: &Spk3Data, epoch: f64) -> Result<State> {
    let (idx, s) = find_record_generic(&data.records, data.init_epoch, epoch)?;
    let record = &data.records[idx];

    let x = clenshaw(&record.x_coeffs, s);
    let y = clenshaw(&record.y_coeffs, s);
    let z = clenshaw(&record.z_coeffs, s);

    let vx = clenshaw(&record.vx_coeffs, s);
    let vy = clenshaw(&record.vy_coeffs, s);
    let vz = clenshaw(&record.vz_coeffs, s);

    Ok(State::new_raw([x, y, z], [vx, vy, vz]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use muad_dib::kernel::spk_types::{ChebyshevRecord, ChebyshevRecordWithVelocity};

    #[test]
    fn test_clenshaw_constant() {
        assert!((clenshaw(&[5.0], 0.0) - 5.0).abs() < 1e-10);
        assert!((clenshaw(&[5.0], 0.5) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_clenshaw_linear() {
        let coeffs = [3.0, 2.0]; // 3 + 2*s
        assert!((clenshaw(&coeffs, 0.0) - 3.0).abs() < 1e-10);
        assert!((clenshaw(&coeffs, 1.0) - 5.0).abs() < 1e-10);
        assert!((clenshaw(&coeffs, -1.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_clenshaw_quadratic() {
        // T_2(s) = 2s^2 - 1
        let coeffs = [0.0, 0.0, 1.0];
        assert!((clenshaw(&coeffs, 0.0) - (-1.0)).abs() < 1e-10);
        assert!((clenshaw(&coeffs, 1.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_derivative_linear() {
        // f(s) = 3 + 2*s, f'(s) = 2
        assert!((clenshaw_derivative(&[3.0, 2.0], 0.0, 1.0) - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_derivative_quadratic() {
        // P(s) = 4s^2 + 3s - 2, P'(s) = 8s + 3
        let coeffs = [0.0, 3.0, 2.0];
        assert!((clenshaw_derivative(&coeffs, 0.0, 1.0) - 3.0).abs() < 1e-10);
        assert!((clenshaw_derivative(&coeffs, 0.5, 1.0) - 7.0).abs() < 1e-10);
    }

    #[test]
    fn test_evaluate_type2() {
        let data = Spk2Data {
            init_epoch: 0.0,
            interval_length: 100.0,
            degree: 1,
            records: vec![ChebyshevRecord {
                midpoint: 50.0,
                radius: 50.0,
                x_coeffs: vec![1000.0, 10.0],
                y_coeffs: vec![2000.0, 20.0],
                z_coeffs: vec![3000.0, 30.0],
            }],
        };

        let state = evaluate_type2(&data, 50.0).unwrap();
        assert!((state.position[0] - 1000.0).abs() < 1e-14);

        let state = evaluate_type2(&data, 100.0).unwrap();
        assert!((state.position[0] - 1010.0).abs() < 1e-14);
        assert!((state.velocity[0] - 0.2).abs() < 1e-14); // 10 / 50
    }

    #[test]
    fn test_evaluate_type2_out_of_range() {
        let data = Spk2Data {
            init_epoch: 0.0,
            interval_length: 100.0,
            degree: 1,
            records: vec![ChebyshevRecord {
                midpoint: 50.0,
                radius: 50.0,
                x_coeffs: vec![1000.0],
                y_coeffs: vec![2000.0],
                z_coeffs: vec![3000.0],
            }],
        };
        assert!(matches!(evaluate_type2(&data, -10.0), Err(Error::EpochOutOfRange { .. })));
        assert!(matches!(evaluate_type2(&data, 110.0), Err(Error::EpochOutOfRange { .. })));
    }

    #[test]
    fn test_evaluate_type3() {
        let data = Spk3Data {
            init_epoch: 0.0,
            interval_length: 100.0,
            degree: 1,
            records: vec![ChebyshevRecordWithVelocity {
                midpoint: 50.0,
                radius: 50.0,
                x_coeffs: vec![1000.0, 10.0],
                y_coeffs: vec![2000.0, 20.0],
                z_coeffs: vec![3000.0, 30.0],
                vx_coeffs: vec![1.0, 0.1],
                vy_coeffs: vec![2.0, 0.2],
                vz_coeffs: vec![3.0, 0.3],
            }],
        };

        let state = evaluate_type3(&data, 50.0).unwrap();
        assert!((state.position[0] - 1000.0).abs() < 1e-14);
        assert!((state.velocity[0] - 1.0).abs() < 1e-14);
    }
}
