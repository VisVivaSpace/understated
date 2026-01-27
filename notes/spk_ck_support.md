# SPK/CK Type Support Status

## SPK Types (Spacecraft/Planetary Kernel — State Vectors)

| Type | Name | Implemented | CSPICE-Validated | Test Data |
|------|------|:-----------:|:----------------:|-----------|
| 2 | Chebyshev (position only) | Yes | Yes | `de440s.bsp` |
| 3 | Chebyshev (position + velocity) | Yes | No | None — need a Type 3 BSP |
| 5 | Two-body / Keplerian | Yes | No | None — need a Type 5 BSP |
| 8 | Lagrange (equal time steps) | Yes | No | None — need a Type 8 BSP |
| 9 | Lagrange (unequal time steps) | Yes | Yes | `gmat-lagrange.bsp` |
| 13 | Hermite (unequal time steps) | Yes | Yes | `gmat-hermite.bsp` |

### Notes

- **Type 2** velocity is computed by differentiating the Chebyshev position polynomials (Clenshaw recurrence on derivatives). This is less precise than types that store velocity coefficients directly.
- **Type 3** stores separate Chebyshev coefficient sets for position and velocity. Not yet validated against CSPICE due to lack of test data.
- **Type 5** uses the universal variable formulation for two-body propagation, handling elliptical, parabolic, and hyperbolic orbits via Stumpff functions.
- **Types 8 and 9** use Lagrange interpolation. Type 8 assumes equal time steps; Type 9 allows unequal steps. Window selection follows the CSPICE algorithm (odd windows: nearest epoch; even windows: lower bracket).
- **Type 13** uses Hermite interpolation, matching both position and velocity at each data point. Window selection matches CSPICE.

## CK Types (C-Kernel — Attitude/Pointing)

| Type | Name | Implemented | CSPICE-Validated | Test Data |
|------|------|:-----------:|:----------------:|-----------|
| 1 | Discrete pointing (nearest) | Yes | No | None — need a Type 1 BC |
| 3 | SLERP quaternion interpolation | Yes | Yes | `test.bc` |

### Notes

- **Type 1** returns the nearest record without interpolation.
- **Type 3** uses SLERP (Spherical Linear Interpolation) between bracketing quaternion records, with antipodal handling and near-identity fallback to linear interpolation.
- Quaternions follow SPICE convention: scalar-first `[q0, q1, q2, q3]`.

## Test Data Gaps

The following types lack test data for CSPICE validation:

1. **SPK Type 3** — Need a BSP file containing Type 3 segments. Could be generated from planetary ephemeris tools or obtained from NAIF archives.
2. **SPK Type 5** — Need a BSP file with two-body propagation segments. These appear in some mission planning files.
3. **SPK Type 8** — Need a BSP file with equal-step Lagrange segments. Less common than Type 9.
4. **CK Type 1** — Need a BC file with discrete pointing segments. Some older mission files use Type 1.

## Validation Accuracy Targets

| Quantity | Target | Unit |
|----------|--------|------|
| Position | 1e-9 | km (1 micrometer) |
| Velocity | 1e-12 | km/s (1 nanometer/second) |
| Quaternion | ~1e-8 | (~0.00001 degrees) |
