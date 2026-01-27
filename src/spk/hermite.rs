//! Hermite polynomial interpolation for SPK Type 13.
//!
//! Matches both position and velocity at each data point using divided
//! differences with duplicated epochs.

use crate::error::{Error, Result};
use crate::state::State;
use muad_dib::kernel::spk_types::Spk13Data;

/// Hermite interpolation for a single component.
fn hermite_interpolate(epochs: &[f64], values: &[f64], derivatives: &[f64], epoch: f64) -> f64 {
    let n = epochs.len();
    debug_assert_eq!(n, values.len());
    debug_assert_eq!(n, derivatives.len());

    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return values[0] + derivatives[0] * (epoch - epochs[0]);
    }

    let m = 2 * n;
    let mut z = vec![0.0; m];
    let mut q = vec![vec![0.0; m]; m];

    for i in 0..n {
        z[2 * i] = epochs[i];
        z[2 * i + 1] = epochs[i];
        q[2 * i][0] = values[i];
        q[2 * i + 1][0] = values[i];
        q[2 * i + 1][1] = derivatives[i];
        if i > 0 {
            q[2 * i][1] = (q[2 * i][0] - q[2 * i - 1][0]) / (z[2 * i] - z[2 * i - 1]);
        }
    }

    for j in 2..m {
        for i in j..m {
            let denom = z[i] - z[i - j];
            if denom.abs() > 1e-15 {
                q[i][j] = (q[i][j - 1] - q[i - 1][j - 1]) / denom;
            }
        }
    }

    let mut result = q[m - 1][m - 1];
    for i in (0..m - 1).rev() {
        result = result * (epoch - z[i]) + q[i][i];
    }
    result
}

/// Hermite derivative interpolation.
fn hermite_interpolate_derivative(
    epochs: &[f64],
    values: &[f64],
    derivatives: &[f64],
    epoch: f64,
) -> f64 {
    let n = epochs.len();
    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return derivatives[0];
    }

    let m = 2 * n;
    let mut z = vec![0.0; m];
    let mut q = vec![vec![0.0; m]; m];

    for i in 0..n {
        z[2 * i] = epochs[i];
        z[2 * i + 1] = epochs[i];
        q[2 * i][0] = values[i];
        q[2 * i + 1][0] = values[i];
        q[2 * i + 1][1] = derivatives[i];
        if i > 0 {
            q[2 * i][1] = (q[2 * i][0] - q[2 * i - 1][0]) / (z[2 * i] - z[2 * i - 1]);
        }
    }

    for j in 2..m {
        for i in j..m {
            let denom = z[i] - z[i - j];
            if denom.abs() > 1e-15 {
                q[i][j] = (q[i][j - 1] - q[i - 1][j - 1]) / denom;
            }
        }
    }

    // Derivative using product rule on Newton's form
    let mut result = 0.0;
    for i in 1..m {
        let mut d_prod = 0.0;
        for k in 0..i {
            let mut term = 1.0;
            for (j, &zj) in z[..i].iter().enumerate() {
                if j != k {
                    term *= epoch - zj;
                }
            }
            d_prod += term;
        }
        result += q[i][i] * d_prod;
    }
    result
}

/// Window selection for Type 13 (same algorithm as Type 9).
#[allow(clippy::if_same_then_else)]
fn select_window(data: &Spk13Data, epoch: f64) -> Result<(usize, usize)> {
    let n = data.states.len();
    if n == 0 {
        return Err(Error::EpochOutOfRange {
            epoch,
            start: 0.0,
            end: 0.0,
        });
    }

    let start_epoch = data.states.first().unwrap().epoch;
    let end_epoch = data.states.last().unwrap().epoch;

    if epoch < start_epoch || epoch > end_epoch {
        return Err(Error::EpochOutOfRange {
            epoch,
            start: start_epoch,
            end: end_epoch,
        });
    }

    let mut lower = 0;
    let mut upper = n - 1;
    while upper - lower > 1 {
        let mid = (lower + upper) / 2;
        if data.states[mid].epoch <= epoch {
            lower = mid;
        } else {
            upper = mid;
        }
    }
    let high = lower + 1;

    let wndsiz = data.window_size as usize;
    let degree = wndsiz - 1;

    let first = if wndsiz % 2 == 1 {
        let near = if lower == 0 {
            lower
        } else if high >= n {
            lower
        } else if (epoch - data.states[lower].epoch).abs()
            <= (data.states[high].epoch - epoch).abs()
        {
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

/// Evaluate SPK Type 13 (Hermite interpolation).
pub fn evaluate_type13(data: &Spk13Data, epoch: f64) -> Result<State> {
    let (start_idx, end_idx) = select_window(data, epoch)?;
    let window_states = &data.states[start_idx..end_idx];

    let epochs: Vec<f64> = window_states.iter().map(|s| s.epoch).collect();
    let x_vals: Vec<f64> = window_states.iter().map(|s| s.x).collect();
    let y_vals: Vec<f64> = window_states.iter().map(|s| s.y).collect();
    let z_vals: Vec<f64> = window_states.iter().map(|s| s.z).collect();
    let vx_vals: Vec<f64> = window_states.iter().map(|s| s.vx).collect();
    let vy_vals: Vec<f64> = window_states.iter().map(|s| s.vy).collect();
    let vz_vals: Vec<f64> = window_states.iter().map(|s| s.vz).collect();

    let x = hermite_interpolate(&epochs, &x_vals, &vx_vals, epoch);
    let y = hermite_interpolate(&epochs, &y_vals, &vy_vals, epoch);
    let z = hermite_interpolate(&epochs, &z_vals, &vz_vals, epoch);

    let vx = hermite_interpolate_derivative(&epochs, &x_vals, &vx_vals, epoch);
    let vy = hermite_interpolate_derivative(&epochs, &y_vals, &vy_vals, epoch);
    let vz = hermite_interpolate_derivative(&epochs, &z_vals, &vz_vals, epoch);

    Ok(State::new_raw([x, y, z], [vx, vy, vz]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use muad_dib::kernel::spk_types::StateRecord;

    #[test]
    fn test_hermite_linear() {
        // f(t) = t: f(0)=0, f'(0)=1, f(1)=1, f'(1)=1
        let epochs = [0.0, 1.0];
        let values = [0.0, 1.0];
        let derivatives = [1.0, 1.0];

        assert!((hermite_interpolate(&epochs, &values, &derivatives, 0.5) - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_evaluate_type13_at_data_points() {
        let data = Spk13Data {
            window_size: 2,
            states: vec![
                StateRecord { epoch: 0.0, x: 0.0, y: 0.0, z: 0.0, vx: 1.0, vy: 0.0, vz: 0.0 },
                StateRecord { epoch: 10.0, x: 10.0, y: 0.0, z: 0.0, vx: 1.0, vy: 0.0, vz: 0.0 },
            ],
        };

        let state = evaluate_type13(&data, 0.0).unwrap();
        assert!((state.position[0] - 0.0).abs() < 1e-6);

        let state = evaluate_type13(&data, 5.0).unwrap();
        assert!((state.position[0] - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_out_of_range() {
        let data = Spk13Data {
            window_size: 2,
            states: vec![
                StateRecord { epoch: 10.0, x: 0.0, y: 0.0, z: 0.0, vx: 0.0, vy: 0.0, vz: 0.0 },
                StateRecord { epoch: 20.0, x: 10.0, y: 0.0, z: 0.0, vx: 0.0, vy: 0.0, vz: 0.0 },
            ],
        };
        assert!(matches!(evaluate_type13(&data, 5.0), Err(Error::EpochOutOfRange { .. })));
    }
}
