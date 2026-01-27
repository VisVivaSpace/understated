//! Rotation matrix utilities matching CSPICE algorithms.
//!
//! Provides quaternion↔matrix conversions, axis-angle extraction, and 3×3
//! matrix operations. All functions match the numerical paths of their CSPICE
//! counterparts to ensure agreement to machine precision.
//!
//! CSPICE equivalents:
//! - [`q2m`] → `q2m_c` — quaternion to rotation matrix
//! - [`m2q`] → `m2q_c` — rotation matrix to quaternion (Shepperd method)
//! - [`raxisa`] → `raxisa_c` — rotation axis and angle from matrix
//! - [`axisar`] → `axisar_c` — rotation matrix from axis and angle
//! - [`mtxm`] → `mtxm_c` — A^T * B
//! - [`mxmt`] → `mxmt_c` — A * B^T
//!
//! Matrices use row-major `[[f64; 3]; 3]` where `m[i][j]` is row i, column j.
//! Quaternions use SPICE convention: scalar-first `[q0, q1, q2, q3]`.

/// 3×3 rotation matrix type (row-major).
pub type Mat3 = [[f64; 3]; 3];

/// Convert a SPICE quaternion to a rotation matrix.
///
/// Matches CSPICE `q2m_c`. Handles non-unit quaternions by normalizing
/// internally (the "sharpening" step). Zero quaternion returns identity.
pub fn q2m(q: &[f64; 4]) -> Mat3 {
    let q01 = q[0] * q[1];
    let q02 = q[0] * q[2];
    let q03 = q[0] * q[3];
    let q12 = q[1] * q[2];
    let q13 = q[1] * q[3];
    let q23 = q[2] * q[3];
    let mut q1s = q[1] * q[1];
    let mut q2s = q[2] * q[2];
    let mut q3s = q[3] * q[3];

    // Sharpen: normalize if not already unit length
    let l2 = q[0] * q[0] + q1s + q2s + q3s;
    let (q01, q02, q03, q12, q13, q23) = if l2 != 1.0 && l2 != 0.0 {
        let s = 1.0 / l2;
        q1s *= s;
        q2s *= s;
        q3s *= s;
        (q01 * s, q02 * s, q03 * s, q12 * s, q13 * s, q23 * s)
    } else {
        (q01, q02, q03, q12, q13, q23)
    };

    [
        [
            1.0 - 2.0 * (q2s + q3s),
            2.0 * (q12 - q03),
            2.0 * (q13 + q02),
        ],
        [
            2.0 * (q12 + q03),
            1.0 - 2.0 * (q1s + q3s),
            2.0 * (q23 - q01),
        ],
        [
            2.0 * (q13 - q02),
            2.0 * (q23 + q01),
            1.0 - 2.0 * (q1s + q2s),
        ],
    ]
}

/// Convert a rotation matrix to a SPICE quaternion.
///
/// Matches CSPICE `m2q_c` (Shepperd method). Always returns a quaternion
/// with `q[0] >= 0`. Skips the `isrot` check since we only call this on
/// matrices we constructed ourselves.
pub fn m2q(m: &Mat3) -> [f64; 4] {
    let trace = m[0][0] + m[1][1] + m[2][2];
    let mtrace = 1.0 - trace;

    // Four candidates for the largest quaternion component squared * 4:
    let cc4 = trace + 1.0; // 4 * q0^2
    let s114 = mtrace + 2.0 * m[0][0]; // 4 * q1^2
    let s224 = mtrace + 2.0 * m[1][1]; // 4 * q2^2
    let s334 = mtrace + 2.0 * m[2][2]; // 4 * q3^2

    // Pick the branch with the largest denominator for numerical stability
    let (mut c, mut s1, mut s2, mut s3);

    if 1.0 <= cc4 {
        c = (cc4 * 0.25).sqrt();
        let factor = 1.0 / (c * 4.0);
        s1 = (m[2][1] - m[1][2]) * factor;
        s2 = (m[0][2] - m[2][0]) * factor;
        s3 = (m[1][0] - m[0][1]) * factor;
    } else if 1.0 <= s114 {
        s1 = (s114 * 0.25).sqrt();
        let factor = 1.0 / (s1 * 4.0);
        c = (m[2][1] - m[1][2]) * factor;
        s2 = (m[0][1] + m[1][0]) * factor;
        s3 = (m[0][2] + m[2][0]) * factor;
    } else if 1.0 <= s224 {
        s2 = (s224 * 0.25).sqrt();
        let factor = 1.0 / (s2 * 4.0);
        c = (m[0][2] - m[2][0]) * factor;
        s1 = (m[0][1] + m[1][0]) * factor;
        s3 = (m[1][2] + m[2][1]) * factor;
    } else {
        s3 = (s334 * 0.25).sqrt();
        let factor = 1.0 / (s3 * 4.0);
        c = (m[1][0] - m[0][1]) * factor;
        s1 = (m[0][2] + m[2][0]) * factor;
        s2 = (m[1][2] + m[2][1]) * factor;
    }

    // Polish: normalize if not already unit length
    let l2 = c * c + s1 * s1 + s2 * s2 + s3 * s3;
    if l2 != 1.0 {
        let polish = 1.0 / l2.sqrt();
        c *= polish;
        s1 *= polish;
        s2 *= polish;
        s3 *= polish;
    }

    // Ensure scalar part is non-negative
    if c > 0.0 {
        [c, s1, s2, s3]
    } else {
        [-c, -s1, -s2, -s3]
    }
}

/// Extract rotation axis and angle from a rotation matrix.
///
/// Matches CSPICE `raxisa_c`. Internally converts to quaternion via [`m2q`],
/// then extracts axis and angle from the quaternion representation.
///
/// Returns `(axis, angle)` where `axis` is a unit vector and `angle` is in
/// `[0, π]`. For the identity matrix, returns `([0, 0, 1], 0.0)`.
pub fn raxisa(m: &Mat3) -> ([f64; 3], f64) {
    let q = m2q(m);

    let vnorm = (q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();

    if vnorm == 0.0 {
        // Identity rotation
        ([0.0, 0.0, 1.0], 0.0)
    } else if q[0] == 0.0 {
        // Rotation by pi
        (
            [q[1], q[2], q[3]], // already unit from m2q
            std::f64::consts::PI,
        )
    } else {
        // General case: unitize the vector part
        let axis = [q[1] / vnorm, q[2] / vnorm, q[3] / vnorm];
        let angle = 2.0 * vnorm.atan2(q[0]);
        (axis, angle)
    }
}

/// Construct a rotation matrix from an axis and angle.
///
/// Matches CSPICE `axisar_c`. Rotates each basis vector by `angle` radians
/// about `axis` using the `vrotv` algorithm (Rodrigues' vector rotation
/// formula).
///
/// If `axis` is the zero vector, returns the identity matrix.
pub fn axisar(axis: &[f64; 3], angle: f64) -> Mat3 {
    let mut r = [[0.0f64; 3]; 3];

    // Identity basis vectors
    let basis: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

    for col in 0..3 {
        let rotated = vrotv(&basis[col], axis, angle);
        // Store as column: r[row][col]
        for row in 0..3 {
            r[row][col] = rotated[row];
        }
    }

    r
}

/// Rotate a vector about an axis by a given angle.
///
/// Matches CSPICE `vrotv_`. Decomposes the vector into components parallel
/// and perpendicular to the axis, then rotates the perpendicular component.
fn vrotv(v: &[f64; 3], axis: &[f64; 3], angle: f64) -> [f64; 3] {
    // Normalize axis
    let norm = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    if norm == 0.0 {
        return *v;
    }
    let x = [axis[0] / norm, axis[1] / norm, axis[2] / norm];

    // Project v onto axis: p = (v · x) * x
    let vdotx = v[0] * x[0] + v[1] * x[1] + v[2] * x[2];
    let p = [vdotx * x[0], vdotx * x[1], vdotx * x[2]];

    // Component orthogonal to axis: v1 = v - p
    let v1 = [v[0] - p[0], v[1] - p[1], v[2] - p[2]];

    // v2 = x × v1 (v1 rotated 90° about axis)
    let v2 = [
        x[1] * v1[2] - x[2] * v1[1],
        x[2] * v1[0] - x[0] * v1[2],
        x[0] * v1[1] - x[1] * v1[0],
    ];

    // result = p + cos(angle)*v1 + sin(angle)*v2
    let c = angle.cos();
    let s = angle.sin();
    [
        p[0] + c * v1[0] + s * v2[0],
        p[1] + c * v1[1] + s * v2[1],
        p[2] + c * v1[2] + s * v2[2],
    ]
}

/// Compute A^T * B for 3×3 matrices.
///
/// Matches CSPICE `mtxm_c`.
pub fn mtxm(a: &Mat3, b: &Mat3) -> Mat3 {
    let mut r = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = a[0][i] * b[0][j] + a[1][i] * b[1][j] + a[2][i] * b[2][j];
        }
    }
    r
}

/// Compute A * B^T for 3×3 matrices.
///
/// Matches CSPICE `mxmt_c`.
pub fn mxmt(a: &Mat3, b: &Mat3) -> Mat3 {
    let mut r = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = a[i][0] * b[j][0] + a[i][1] * b[j][1] + a[i][2] * b[j][2];
        }
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_1_SQRT_2, PI};

    fn assert_mat_eq(a: &Mat3, b: &Mat3, tol: f64) {
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (a[i][j] - b[i][j]).abs() < tol,
                    "m[{}][{}]: {} vs {} (diff {})",
                    i,
                    j,
                    a[i][j],
                    b[i][j],
                    (a[i][j] - b[i][j]).abs()
                );
            }
        }
    }

    fn assert_quat_eq(a: &[f64; 4], b: &[f64; 4], tol: f64) {
        // Quaternions q and -q represent the same rotation
        let dot: f64 = a.iter().zip(b.iter()).map(|(x, y)| x * y).sum();
        let sign = if dot < 0.0 { -1.0 } else { 1.0 };
        for i in 0..4 {
            assert!(
                (a[i] - sign * b[i]).abs() < tol,
                "q[{}]: {} vs {} (diff {})",
                i,
                a[i],
                sign * b[i],
                (a[i] - sign * b[i]).abs()
            );
        }
    }

    #[test]
    fn test_q2m_identity() {
        let q = [1.0, 0.0, 0.0, 0.0];
        let m = q2m(&q);
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        assert_mat_eq(&m, &id, 1e-15);
    }

    #[test]
    fn test_q2m_90deg_about_z() {
        // Q = (cos(pi/4), 0, 0, -sin(pi/4)) = (sqrt(2)/2, 0, 0, -sqrt(2)/2)
        // This is a rotation by pi/2 about the -z axis, giving [pi/2]_3.
        let q = [FRAC_1_SQRT_2, 0.0, 0.0, -FRAC_1_SQRT_2];
        let m = q2m(&q);
        // Expected: [[0, 1, 0], [-1, 0, 0], [0, 0, 1]]
        let expected = [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        assert_mat_eq(&m, &expected, 1e-15);
    }

    #[test]
    fn test_m2q_roundtrip() {
        let q_orig = [FRAC_1_SQRT_2, 0.0, 0.0, -FRAC_1_SQRT_2];
        let m = q2m(&q_orig);
        let q_back = m2q(&m);
        assert_quat_eq(&q_orig, &q_back, 1e-15);
    }

    #[test]
    fn test_raxisa_identity() {
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let (axis, angle) = raxisa(&id);
        assert_eq!(axis, [0.0, 0.0, 1.0]);
        assert!((angle).abs() < 1e-15);
    }

    #[test]
    fn test_axisar_roundtrip() {
        let axis = [0.0, 0.0, 1.0];
        let angle = PI / 4.0;
        let m = axisar(&axis, angle);
        let (axis_out, angle_out) = raxisa(&m);
        assert!((angle_out - angle).abs() < 1e-15);
        for i in 0..3 {
            assert!((axis_out[i] - axis[i]).abs() < 1e-15);
        }
    }

    #[test]
    fn test_mtxm() {
        // A^T * A should give identity for any rotation matrix
        let q = [0.5, 0.5, 0.5, 0.5];
        let a = q2m(&q);
        let ata = mtxm(&a, &a);
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        assert_mat_eq(&ata, &id, 1e-15);
    }

    #[test]
    fn test_mxmt() {
        // A * A^T should give identity for any rotation matrix
        let q = [0.5, 0.5, 0.5, 0.5];
        let a = q2m(&q);
        let aat = mxmt(&a, &a);
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        assert_mat_eq(&aat, &id, 1e-15);
    }
}
