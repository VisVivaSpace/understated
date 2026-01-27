//! Two-body (Keplerian) propagation for SPK Type 5.
//!
//! Uses universal variable formulation to handle elliptical, parabolic,
//! and hyperbolic orbits with Stumpff functions.

use crate::error::{Error, Result};
use crate::state::State;
use muad_dib::kernel::spk_types::Spk5Data;

const MAX_ITERATIONS: usize = 50;
const TOLERANCE: f64 = 1e-14;

/// Propagate a state using two-body dynamics.
pub fn propagate(r0: [f64; 3], v0: [f64; 3], gm: f64, dt: f64) -> State {
    if dt.abs() < 1e-12 {
        return State::new_raw(r0, v0);
    }

    let r_mag = (r0[0] * r0[0] + r0[1] * r0[1] + r0[2] * r0[2]).sqrt();
    let v_mag = (v0[0] * v0[0] + v0[1] * v0[1] + v0[2] * v0[2]).sqrt();

    if r_mag < 1e-10 {
        return State::new_raw(r0, v0);
    }

    let energy = v_mag * v_mag / 2.0 - gm / r_mag;
    let a = if energy.abs() < 1e-15 {
        f64::INFINITY
    } else {
        -gm / (2.0 * energy)
    };

    let r_dot = (r0[0] * v0[0] + r0[1] * v0[1] + r0[2] * v0[2]) / r_mag;

    let (r_new, v_new) = propagate_universal(r0, v0, gm, dt, a, r_mag, r_dot);
    State::new_raw(r_new, v_new)
}

fn propagate_universal(
    r0: [f64; 3],
    v0: [f64; 3],
    gm: f64,
    dt: f64,
    a: f64,
    r_mag: f64,
    r_dot: f64,
) -> ([f64; 3], [f64; 3]) {
    let sqrt_gm = gm.sqrt();
    let alpha = 1.0 / a;

    let mut chi = if a.is_finite() && a > 0.0 {
        sqrt_gm * dt.abs() * alpha.abs()
    } else {
        sqrt_gm * dt.abs() / r_mag
    };

    for _ in 0..MAX_ITERATIONS {
        let chi2 = chi * chi;
        let psi = chi2 * alpha;
        let (c2, c3) = stumpff(psi);

        let r = chi2 * c2 + r_dot / sqrt_gm * chi * (1.0 - psi * c3) + r_mag * (1.0 - psi * c2);

        let f_chi = r_mag * (1.0 - psi * c2) * chi
            + r_dot / sqrt_gm * chi2 * (1.0 - psi * c3)
            + chi2 * chi * c3
            - sqrt_gm * dt.abs();

        let f_prime = r;

        if f_prime.abs() < 1e-15 {
            break;
        }

        let delta = f_chi / f_prime;
        chi -= delta;

        if delta.abs() < TOLERANCE * chi.abs().max(1.0) {
            break;
        }
    }

    if dt < 0.0 {
        chi = -chi;
    }

    let chi2 = chi * chi;
    let psi = chi2 * alpha;
    let (c2, c3) = stumpff(psi);

    let f = 1.0 - chi2 * c2 / r_mag;
    let g = dt - chi2 * chi * c3 / sqrt_gm;

    let r_new = [
        f * r0[0] + g * v0[0],
        f * r0[1] + g * v0[1],
        f * r0[2] + g * v0[2],
    ];

    let r_new_mag = (r_new[0] * r_new[0] + r_new[1] * r_new[1] + r_new[2] * r_new[2]).sqrt();

    let f_dot = sqrt_gm * chi * (psi * c3 - 1.0) / (r_mag * r_new_mag);
    let g_dot = 1.0 - chi2 * c2 / r_new_mag;

    let v_new = [
        f_dot * r0[0] + g_dot * v0[0],
        f_dot * r0[1] + g_dot * v0[1],
        f_dot * r0[2] + g_dot * v0[2],
    ];

    (r_new, v_new)
}

/// Stumpff functions C2 and C3.
fn stumpff(psi: f64) -> (f64, f64) {
    if psi.abs() < 1e-10 {
        let c2 = 0.5 - psi / 24.0 + psi * psi / 720.0;
        let c3 = 1.0 / 6.0 - psi / 120.0 + psi * psi / 5040.0;
        (c2, c3)
    } else if psi > 0.0 {
        let sqrt_psi = psi.sqrt();
        let c2 = (1.0 - sqrt_psi.cos()) / psi;
        let c3 = (sqrt_psi - sqrt_psi.sin()) / (psi * sqrt_psi);
        (c2, c3)
    } else {
        let sqrt_neg_psi = (-psi).sqrt();
        let c2 = (1.0 - sqrt_neg_psi.cosh()) / psi;
        let c3 = (sqrt_neg_psi.sinh() - sqrt_neg_psi) / ((-psi) * sqrt_neg_psi);
        (c2, c3)
    }
}

fn find_nearest_state(data: &Spk5Data, epoch: f64) -> Result<(usize, f64)> {
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

    if epoch < start_epoch - 86400.0 * 365.0 || epoch > end_epoch + 86400.0 * 365.0 {
        return Err(Error::EpochOutOfRange {
            epoch,
            start: start_epoch,
            end: end_epoch,
        });
    }

    let mut best_idx = 0;
    let mut best_diff = (data.states[0].epoch - epoch).abs();

    for (i, state) in data.states.iter().enumerate() {
        let diff = (state.epoch - epoch).abs();
        if diff < best_diff {
            best_diff = diff;
            best_idx = i;
        }
    }

    let dt = epoch - data.states[best_idx].epoch;
    Ok((best_idx, dt))
}

/// Evaluate SPK Type 5 (two-body propagation).
pub fn evaluate_type5(data: &Spk5Data, epoch: f64) -> Result<State> {
    let (idx, dt) = find_nearest_state(data, epoch)?;
    let state = &data.states[idx];

    let r0 = [state.x, state.y, state.z];
    let v0 = [state.vx, state.vy, state.vz];

    Ok(propagate(r0, v0, data.gm, dt))
}

#[cfg(test)]
mod tests {
    use super::*;
    use muad_dib::kernel::spk_types::StateRecord;
    use std::f64::consts::PI;

    const EARTH_GM: f64 = 398600.435;

    #[test]
    fn test_stumpff_parabolic() {
        let (c2, c3) = stumpff(0.0);
        assert!((c2 - 0.5).abs() < 1e-10);
        assert!((c3 - 1.0 / 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_propagate_zero_time() {
        let r0 = [6678.0, 0.0, 0.0];
        let v0 = [0.0, 7.73, 0.0];
        let state = propagate(r0, v0, EARTH_GM, 0.0);
        assert!((state.position[0] - 6678.0).abs() < 1e-10);
    }

    #[test]
    fn test_circular_orbit() {
        let r = 6678.0;
        let v = (EARTH_GM / r).sqrt();
        let r0 = [r, 0.0, 0.0];
        let v0 = [0.0, v, 0.0];
        let period = 2.0 * PI * (r * r * r / EARTH_GM).sqrt();

        let state = propagate(r0, v0, EARTH_GM, period);
        assert!((state.position[0] - r).abs() < 1.0);
        assert!(state.position[1].abs() < 1.0);
    }

    #[test]
    fn test_evaluate_type5() {
        let data = Spk5Data {
            gm: EARTH_GM,
            states: vec![StateRecord {
                epoch: 0.0,
                x: 6678.0,
                y: 0.0,
                z: 0.0,
                vx: 0.0,
                vy: 7.73,
                vz: 0.0,
            }],
        };

        let state = evaluate_type5(&data, 0.0).unwrap();
        assert!((state.position[0] - 6678.0).abs() < 1e-6);
    }

    #[test]
    fn test_out_of_range() {
        let data = Spk5Data {
            gm: EARTH_GM,
            states: vec![StateRecord {
                epoch: 0.0,
                x: 6678.0,
                y: 0.0,
                z: 0.0,
                vx: 0.0,
                vy: 7.73,
                vz: 0.0,
            }],
        };
        assert!(matches!(evaluate_type5(&data, 1e10), Err(Error::EpochOutOfRange { .. })));
    }
}
