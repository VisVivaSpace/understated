//! Spherical Linear Interpolation (SLERP) for quaternions.
//!
//! Features antipodal handling and near-identity fallback to linear interpolation.

/// SLERP between two quaternions.
///
/// Interpolates along the shortest path on the 4D unit sphere.
///
/// # Arguments
///
/// * `q0` - Start quaternion [scalar, i, j, k]
/// * `q1` - End quaternion [scalar, i, j, k]
/// * `t` - Interpolation parameter [0, 1]
pub fn slerp(q0: &[f64; 4], q1: &[f64; 4], t: f64) -> [f64; 4] {
    fn normalize_quat(q: [f64; 4]) -> [f64; 4] {
        let norm = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
        if norm < 1e-15 {
            return q;
        }
        [q[0] / norm, q[1] / norm, q[2] / norm, q[3] / norm]
    }

    if t <= 0.0 {
        return normalize_quat(*q0);
    }
    if t >= 1.0 {
        return normalize_quat(*q1);
    }

    let mut dot = q0[0] * q1[0] + q0[1] * q1[1] + q0[2] * q1[2] + q0[3] * q1[3];

    // Antipodal handling: negate one quaternion to take shorter path
    let mut q1_adj = *q1;
    if dot < 0.0 {
        q1_adj = [-q1[0], -q1[1], -q1[2], -q1[3]];
        dot = -dot;
    }

    // Near-identity fallback to linear interpolation
    if dot > 0.9995 {
        let result = [
            q0[0] + t * (q1_adj[0] - q0[0]),
            q0[1] + t * (q1_adj[1] - q0[1]),
            q0[2] + t * (q1_adj[2] - q0[2]),
            q0[3] + t * (q1_adj[3] - q0[3]),
        ];
        return normalize_quat(result);
    }

    // Standard SLERP
    let theta_0 = dot.acos();
    let theta = theta_0 * t;

    let sin_theta = theta.sin();
    let sin_theta_0 = theta_0.sin();

    let s0 = theta.cos() - dot * sin_theta / sin_theta_0;
    let s1 = sin_theta / sin_theta_0;

    let result = [
        s0 * q0[0] + s1 * q1_adj[0],
        s0 * q0[1] + s1 * q1_adj[1],
        s0 * q0[2] + s1 * q1_adj[2],
        s0 * q0[3] + s1 * q1_adj[3],
    ];

    normalize_quat(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_endpoints() {
        let q0 = [1.0, 0.0, 0.0, 0.0];
        let half_angle = std::f64::consts::FRAC_PI_4;
        let q1 = [half_angle.cos(), half_angle.sin(), 0.0, 0.0];

        let r0 = slerp(&q0, &q1, 0.0);
        assert!((r0[0] - 1.0).abs() < 1e-10);

        let r1 = slerp(&q0, &q1, 1.0);
        assert!((r1[0] - q1[0]).abs() < 1e-10);
    }

    #[test]
    fn test_midpoint() {
        let q0 = [1.0, 0.0, 0.0, 0.0];
        let q1 = [0.0, 1.0, 0.0, 0.0]; // 180 deg around X

        let result = slerp(&q0, &q1, 0.5);
        let expected = (PI / 4.0).cos();
        assert!((result[0] - expected).abs() < 1e-6);
    }

    #[test]
    fn test_normalization() {
        let q0 = [1.0, 0.0, 0.0, 0.0];
        let q1 = [0.0, 1.0, 0.0, 0.0];

        for i in 0..=10 {
            let t = i as f64 / 10.0;
            let result = slerp(&q0, &q1, t);
            let norm = (result[0].powi(2) + result[1].powi(2) + result[2].powi(2) + result[3].powi(2)).sqrt();
            assert!((norm - 1.0).abs() < 1e-10, "Not normalized at t={}: {}", t, norm);
        }
    }

    #[test]
    fn test_antipodal() {
        let q0 = [1.0, 0.0, 0.0, 0.0];
        let q1 = [-1.0, 0.0, 0.0, 0.0]; // Same rotation

        let result = slerp(&q0, &q1, 0.5);
        let norm = (result[0].powi(2) + result[1].powi(2) + result[2].powi(2) + result[3].powi(2)).sqrt();
        assert!((norm - 1.0).abs() < 1e-10);
    }
}
