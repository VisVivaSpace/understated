//! Lagrange polynomial interpolation for SPK Types 8 and 9.

use crate::error::{Error, Result};
use crate::state::State;
use muad_dib::kernel::spk_types::{Spk8Data, Spk9Data, StateRecord};

/// Lagrange interpolation at a point (used in tests).
#[cfg(test)]
fn lagrange_interpolate(x_values: &[f64], y_values: &[f64], x: f64) -> f64 {
    let n = x_values.len();
    debug_assert_eq!(n, y_values.len());

    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return y_values[0];
    }

    let mut result = 0.0;
    for i in 0..n {
        let mut basis = 1.0;
        for j in 0..n {
            if i != j {
                let denom = x_values[i] - x_values[j];
                if denom.abs() > 1e-15 {
                    basis *= (x - x_values[j]) / denom;
                }
            }
        }
        result += y_values[i] * basis;
    }
    result
}

/// Compute Lagrange basis weights for all data points at evaluation point `x`.
///
/// Returns a vector of basis values `L_i(x)` such that the interpolated
/// value for any component is `sum(y_i * L_i(x))`.
fn lagrange_basis(x_values: &[f64], x: f64) -> Vec<f64> {
    let n = x_values.len();
    let mut basis = vec![1.0; n];
    for i in 0..n {
        for j in 0..n {
            if i != j {
                let denom = x_values[i] - x_values[j];
                if denom.abs() > 1e-15 {
                    basis[i] *= (x - x_values[j]) / denom;
                }
            }
        }
    }
    basis
}

/// Apply precomputed Lagrange basis to interpolate all 6 state components at once.
fn interpolate_state(states: &[StateRecord], basis: &[f64]) -> [f64; 6] {
    let mut result = [0.0f64; 6];
    for (i, s) in states.iter().enumerate() {
        let w = basis[i];
        result[0] += s.x * w;
        result[1] += s.y * w;
        result[2] += s.z * w;
        result[3] += s.vx * w;
        result[4] += s.vy * w;
        result[5] += s.vz * w;
    }
    result
}

/// Window selection for Type 8 (equally spaced).
///
/// Implements CSPICE algorithm: odd windows center on nearest epoch,
/// even windows use lower bracket.
#[allow(clippy::if_same_then_else)]
fn select_window_type8(data: &Spk8Data, epoch: f64) -> Result<(usize, usize)> {
    let n = data.states.len();
    if n == 0 {
        return Err(Error::EpochOutOfRange {
            epoch,
            start: data.start_epoch,
            end: data.start_epoch,
        });
    }

    let end_epoch = data.start_epoch + (n - 1) as f64 * data.step_size;
    if epoch < data.start_epoch || epoch > end_epoch {
        return Err(Error::EpochOutOfRange {
            epoch,
            start: data.start_epoch,
            end: end_epoch,
        });
    }

    let normalized = (epoch - data.start_epoch) / data.step_size;
    let lower = normalized.floor() as usize;
    let high = (lower + 1).min(n - 1);

    // CSPICE stores the polynomial degree, not the window size.
    // wndsiz (number of points) = degree + 1.
    let degree = data.window_size as usize;
    let wndsiz = degree + 1;

    let first = if wndsiz % 2 == 1 {
        let lower_epoch = data.start_epoch + lower as f64 * data.step_size;
        let high_epoch = data.start_epoch + high as f64 * data.step_size;
        let near = if lower == 0 {
            lower
        } else if high >= n {
            lower
        } else if (epoch - lower_epoch).abs() <= (high_epoch - epoch).abs() {
            lower
        } else {
            high
        };
        let half = degree / 2;
        if near < half {
            0
        } else if near > n - 1 - (degree - half) {
            n - wndsiz
        } else {
            near - half
        }
    } else {
        let half = degree / 2;
        if lower < half {
            0
        } else if lower > n - 1 - (degree - half) {
            n - wndsiz
        } else {
            lower - half
        }
    };

    Ok((first, first + degree + 1))
}

/// Window selection for Type 9 (unequally spaced).
fn select_window_type9(data: &Spk9Data, epoch: f64) -> Result<(usize, usize)> {
    let degree = data.window_size as usize;
    super::select_window_unequal(data.states.len(), degree, epoch, |i| data.states[i].epoch)
}

/// Evaluate SPK Type 8 (Lagrange, equally spaced).
pub fn evaluate_type8(data: &Spk8Data, epoch: f64) -> Result<State> {
    let (start_idx, end_idx) = select_window_type8(data, epoch)?;
    let window_states = &data.states[start_idx..end_idx];

    let epochs: Vec<f64> = (start_idx..end_idx)
        .map(|i| data.start_epoch + (i as f64) * data.step_size)
        .collect();

    let basis = lagrange_basis(&epochs, epoch);
    let c = interpolate_state(window_states, &basis);
    Ok(State::new_raw([c[0], c[1], c[2]], [c[3], c[4], c[5]]))
}

/// Evaluate SPK Type 9 (Lagrange, unequally spaced).
pub fn evaluate_type9(data: &Spk9Data, epoch: f64) -> Result<State> {
    let (start_idx, end_idx) = select_window_type9(data, epoch)?;
    let window_states = &data.states[start_idx..end_idx];

    let epochs: Vec<f64> = window_states.iter().map(|s| s.epoch).collect();

    let basis = lagrange_basis(&epochs, epoch);
    let c = interpolate_state(window_states, &basis);
    Ok(State::new_raw([c[0], c[1], c[2]], [c[3], c[4], c[5]]))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lagrange_linear() {
        let x = [0.0, 1.0];
        let y = [1.0, 3.0]; // y = 2x + 1
        assert!((lagrange_interpolate(&x, &y, 0.5) - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_lagrange_quadratic() {
        let x = [0.0, 1.0, 2.0];
        let y = [0.0, 1.0, 4.0]; // y = x^2
        assert!((lagrange_interpolate(&x, &y, 1.5) - 2.25).abs() < 1e-10);
    }

    #[test]
    fn test_evaluate_type8() {
        // window_size=3 means degree 3 (cubic), so 4 points needed
        let data = Spk8Data {
            start_epoch: 0.0,
            step_size: 10.0,
            window_size: 3,
            states: (0..4)
                .map(|i| StateRecord {
                    epoch: (i * 10) as f64,
                    x: (i * 10) as f64,
                    y: 0.0,
                    z: 0.0,
                    vx: 1.0,
                    vy: 0.0,
                    vz: 0.0,
                })
                .collect(),
        };

        let state = evaluate_type8(&data, 15.0).unwrap();
        assert!((state.position[0] - 15.0).abs() < 1e-12);
    }

    #[test]
    fn test_evaluate_type9() {
        // window_size=3 means degree 3 (cubic), so 4 points needed
        let data = Spk9Data {
            window_size: 3,
            states: vec![
                StateRecord { epoch: 0.0, x: 0.0, y: 0.0, z: 0.0, vx: 1.0, vy: 0.0, vz: 0.0 },
                StateRecord { epoch: 5.0, x: 5.0, y: 0.0, z: 0.0, vx: 1.0, vy: 0.0, vz: 0.0 },
                StateRecord { epoch: 15.0, x: 15.0, y: 0.0, z: 0.0, vx: 1.0, vy: 0.0, vz: 0.0 },
                StateRecord { epoch: 30.0, x: 30.0, y: 0.0, z: 0.0, vx: 1.0, vy: 0.0, vz: 0.0 },
            ],
        };

        let state = evaluate_type9(&data, 5.0).unwrap();
        assert!((state.position[0] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_out_of_range() {
        let data = Spk8Data {
            start_epoch: 0.0,
            step_size: 10.0,
            window_size: 1,
            states: vec![
                StateRecord { epoch: 0.0, x: 0.0, y: 0.0, z: 0.0, vx: 0.0, vy: 0.0, vz: 0.0 },
                StateRecord { epoch: 10.0, x: 10.0, y: 0.0, z: 0.0, vx: 0.0, vy: 0.0, vz: 0.0 },
            ],
        };
        assert!(matches!(evaluate_type8(&data, -5.0), Err(Error::EpochOutOfRange { .. })));
    }
}
