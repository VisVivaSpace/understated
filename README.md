# understated

**Status: Experimental** — not all SPK/CK types are validated against CSPICE.

Spacecraft state and attitude interpolation from NAIF SPICE kernel data.

## What It Does

understated evaluates SPK (trajectory) and CK (attitude) data from SPICE kernel files. Given parsed segment data from [muad-dib](https://github.com/nstrange/muad-dib), it computes:

- **State vectors** (position + velocity) for planetary bodies and spacecraft at arbitrary TDB epochs
- **Pointing data** (quaternion + angular velocity) for instruments at arbitrary SCLK times
- **Coordinate transforms** between rectangular, latitudinal, spherical, and cylindrical systems
- **Time conversions** between TDB, UTC, calendar, ISO 8601, and Julian Date formats

## Architecture

```
understated (evaluation)  <--depends-->  muad-dib (file I/O)
```

- **muad-dib** parses SPICE kernel binary formats (DAF/SPK/CK/LSK) and provides raw segment data
- **understated** implements all interpolation algorithms and exposes a high-level `Ephemeris` API

## Supported Types

| SPK Type | Method | CSPICE-Validated |
|----------|--------|:----------------:|
| 2 | Chebyshev (position only) | Yes |
| 3 | Chebyshev (position + velocity) | No |
| 5 | Two-body / Keplerian | No |
| 8 | Lagrange (equal time steps) | No |
| 9 | Lagrange (unequal time steps) | Yes |
| 13 | Hermite (unequal time steps) | Yes |

| CK Type | Method | CSPICE-Validated |
|---------|--------|:----------------:|
| 1 | Discrete pointing | No |
| 3 | Linear interpolation (axis-angle, rotation matrix; matches CSPICE CKE03) | Yes |

## Usage

```rust
use understated::{Ephemeris, EpochTDB, NaifId};

let eph = Ephemeris::load("de440s.bsp")?;

let epoch = EpochTDB::parse("2020-01-01T12:00:00")?;
let state = eph.state_of(NaifId::EARTH, epoch, NaifId::SSB)?;

println!("Earth position: {:?} km", state.position);
println!("Earth velocity: {:?} km/s", state.velocity);
println!("Distance from SSB: {:.3} km", state.distance());
```

## Building

```bash
cargo build                    # Build the crate
cargo test                     # Run unit tests (no test data needed)
cargo test --features test-data  # Run integration tests with SPICE kernels
cargo clippy                   # Lint
```

## Examples

```bash
# Query planetary ephemeris
cargo run --example query_ephemeris --features test-data -- test_data/de440s.bsp

# Coordinate system conversions (no data files needed)
cargo run --example coordinate_transforms

# Time string parsing and formatting (no data files needed)
cargo run --example time_conversion
```

## Test Data

Binary test data files (`.bsp`, `.bc`, `.bpc`, `.hdf5`) are tracked with Git LFS. After cloning:

```bash
git lfs pull
```

Text kernel files (`.tls`, `.tpc`) and `test_data/README.md` are tracked as regular git objects.

## Validation

Where CSPICE validation is available, accuracy targets are:

| Quantity | Target |
|----------|--------|
| Position | 1 micrometer (1e-9 km) |
| Velocity | 1 nm/s (1e-12 km/s) |
| Quaternion | ~1e-8 (~0.00001 degrees) |

Run CSPICE validation tests (requires the [CSPICE toolkit](https://naif.jpl.nasa.gov/naif/toolkit_C.html)):

```bash
CSPICE_DIR=/path/to/cspice cargo test --features cspice,test-data
```

Set `CSPICE_DIR` to your CSPICE installation directory (the one containing `lib/`).
You can also export it in your shell profile to avoid repeating it:

```bash
export CSPICE_DIR=/path/to/cspice
cargo test --features cspice,test-data
```
