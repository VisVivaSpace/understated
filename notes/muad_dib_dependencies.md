# muad-dib API Dependencies

This document lists every muad-dib API that understated depends on. These APIs must remain stable when refactoring muad-dib into a pure I/O crate.

## Kernel Loading & Segment Access

### `muad_dib::kernel::SpiceKernel`

| Method | Used In | Purpose |
|--------|---------|---------|
| `load(path)` | `ephemeris.rs` | Load a single SPICE kernel file |
| `load_many(paths)` | `ephemeris.rs` | Load multiple kernel files |
| `spk_segments()` | `spk/chain.rs` | Iterate all SPK segments |
| `spk_segments_for(body)` | `spk/chain.rs` | Iterate SPK segments for a specific body |
| `spk_view(segment)` | `spk/chain.rs` | Get a view into segment data for interpolation |
| `ck_segments_for(instrument)` | `ephemeris.rs` | Find CK segments for an instrument |
| `ck_view(segment)` | `ephemeris.rs` | Get a view into CK segment data |
| `spk_bodies()` | `ephemeris.rs` | List all bodies with SPK coverage |
| `ck_instruments()` | `ephemeris.rs` | List all instruments with CK coverage |

## SPK Data Types

### `muad_dib::kernel::spk_types`

| Type | Used In | Purpose |
|------|---------|---------|
| `SpkData` | `spk/mod.rs` | Enum for dispatching to type-specific interpolation |
| `Spk2Data` | `spk/chebyshev.rs` | Type 2 Chebyshev data (position-only coefficients) |
| `Spk3Data` | `spk/chebyshev.rs` | Type 3 Chebyshev data (position + velocity coefficients) |
| `Spk5Data` | `spk/twobody.rs` | Type 5 two-body propagation parameters |
| `Spk8Data` | `spk/lagrange.rs` | Type 8 Lagrange data (equal time steps) |
| `Spk9Data` | `spk/lagrange.rs` | Type 9 Lagrange data (unequal time steps) |
| `Spk13Data` | `spk/hermite.rs` | Type 13 Hermite data (unequal time steps) |
| `StateRecord` | `spk/lagrange.rs`, `spk/hermite.rs` | Position + velocity record `{ epoch, state[6] }` |
| `ChebyshevRecord` | `spk/chebyshev.rs` | Chebyshev coefficient record (position) |
| `ChebyshevRecordWithVelocity` | `spk/chebyshev.rs` | Chebyshev coefficient record (position + velocity) |

### Fields accessed on SPK data structs

- `Spk2Data`: `records`, `record_size`, `initial_epoch`, `interval_length`
- `Spk3Data`: `records`, `record_size`, `initial_epoch`, `interval_length`
- `Spk5Data`: `gm`, `initial_epoch`, `state` (6-element array)
- `Spk8Data`: `records`, `step_size`, `initial_epoch`, `window_size`
- `Spk9Data`: `records`, `window_size`
- `Spk13Data`: `records`, `window_size`
- `StateRecord`: `epoch`, `state` (6-element array)
- `ChebyshevRecord`: `mid_epoch`, `radius`, `coefficients` (Vec of 3-element arrays)
- `ChebyshevRecordWithVelocity`: `mid_epoch`, `radius`, `position_coefficients`, `velocity_coefficients`

## CK Data Types

### `muad_dib::kernel::ck_types`

| Type | Used In | Purpose |
|------|---------|---------|
| `CkData` | `ck/mod.rs` | Enum for dispatching to type-specific evaluation |
| `Ck1Data` | `ck/evaluate.rs` | Type 1 discrete pointing data |
| `Ck3Data` | `ck/evaluate.rs` | Type 3 SLERP interpolation data |
| `PointingRecord` | `ck/evaluate.rs` | Quaternion + optional angular velocity record |

### Fields accessed on CK data structs

- `Ck1Data`: `records`
- `Ck3Data`: `records`, `interval_starts`
- `PointingRecord`: `sclk`, `quaternion` (4-element array), `angular_velocity` (Option of 3-element array)

## SPK Segment Metadata

### `SpkSegment` (from `spk_view` / `spk_segments`)

| Field | Type | Purpose |
|-------|------|---------|
| `target_code` | `i32` | NAIF ID of the target body |
| `center_code` | `i32` | NAIF ID of the center body |
| `frame_code` | `i32` | NAIF reference frame ID |
| `spk_type` | `i32` | SPK segment type (2, 3, 5, 8, 9, 13) |
| `initial_epoch` | `f64` | Start of segment coverage (TDB seconds past J2000) |
| `final_epoch` | `f64` | End of segment coverage (TDB seconds past J2000) |

## CK Segment Metadata

### `CkSegment` (from `ck_segments_for`)

| Field | Type | Purpose |
|-------|------|---------|
| `instrument_code` | `i32` | NAIF instrument/structure ID |
| `frame_code` | `i32` | NAIF reference frame ID |
| `initial_sclk` | `f64` | Start of segment coverage (SCLK ticks) |
| `final_sclk` | `f64` | End of segment coverage (SCLK ticks) |

## LSK / Time

### `muad_dib::spice::lsk::LeapSecondExt`

| Method | Used In | Purpose |
|--------|---------|---------|
| `lsk_data()` | `lsk.rs` | Extract leap second data from kernel pool |

Returns a struct with fields: `delta_t_a`, `k`, `eb`, `m` (2-element array), `leap_seconds` (Vec of tuples).

## Type Conversions

### `muad_dib::types::NaifId`

Used for conversion between understated's `NaifId` and muad-dib's `NaifId` via `From` implementations. The only field accessed is `.0: i32`.

## Summary

understated uses muad-dib exclusively for:

1. **File I/O** — Loading SPICE kernel files and parsing their binary formats
2. **Data extraction** — Providing parsed segment data (coefficients, records, parameters)
3. **Segment discovery** — Finding which segments cover a given body/instrument/epoch
4. **LSK extraction** — Reading leap second parameters from kernel pool

All interpolation, coordinate transforms, time formatting, and state chaining logic lives in understated.
