# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**understated** is a standalone Rust crate for interpolating spacecraft state (position + velocity) and attitude (quaternion + angular velocity) data from NAIF SPICE kernel segments. It extracts the interpolation logic from [muad-dib](/Users/nstrange/git/open_source/muad-dib) into a separate, publishable crate.

This crate does NOT read SPICE files directly — it relies on muad-dib for file parsing. If new parsing features are needed, create a markdown file explaining the requirement rather than adding file I/O here.

### Upstream Reference

- **muad-dib source:** `/Users/nstrange/git/open_source/muad-dib`
- **Interpolation code to extract:** `muad-dib/src/spice/interpolate/` (chebyshev, lagrange, hermite, twobody modules)
- **CK attitude code:** `muad-dib/src/spice/ck.rs` (SLERP, discrete pointing)
- **SPK state API:** `muad-dib/src/spice/spk.rs` (segment evaluation, body chaining)
- **Type support status:** `muad-dib/docs/SPK_CK_TYPE_SUPPORT.md`
- **Examples for feature reference:** `muad-dib/examples/` (query_ephemeris, coordinate_transforms, etc.)

### Scope

**SPK interpolation (state):**
- Type 2: Chebyshev (position only, velocity via differentiation)
- Type 3: Chebyshev (position + velocity coefficients)
- Type 5: Two-body/Keplerian propagation (universal variables)
- Type 8: Lagrange (equal time steps)
- Type 9: Lagrange (unequal time steps)
- Type 13: Hermite (unequal time steps, matches position + velocity)

**CK interpolation (attitude):**
- Type 1: Discrete pointing (nearest record)
- Type 3: SLERP quaternion interpolation between bracketing records

### Key Data Types (from muad-dib)

- `State` — target, center, frame, position `[f64; 3]`, velocity `[f64; 3]` with Add/Sub/Neg operators for body chaining
- `Pointing` — frame, quaternion `[f64; 4]` (scalar-first SPICE convention), optional angular velocity `[f64; 3]`

### Numerical Methods

- **Clenshaw recurrence** for Chebyshev evaluation (numerically stable)
- **Universal variable formulation** for two-body (handles elliptical/parabolic/hyperbolic)
- **Stumpff functions** for orbit type transitions
- **SLERP** with antipodal handling and near-identity fallback to linear interpolation
- **Window selection** must match CSPICE algorithm exactly (odd windows: nearest epoch; even windows: lower bracket)

## Build Commands

```bash
cargo build                    # Build the crate
cargo test                     # Run all tests
cargo test <test_name>         # Run a single test
cargo clippy                   # Lint
cargo doc --open               # Generate and view docs
```

Rust edition: 2024

## Testing

Use CSPICE at `/Users/nstrange/cspice` for validation testing only. Do NOT use CSPICE or any external library as a runtime dependency.

Validation accuracy targets (from muad-dib):
- Position: 1 micrometer (1e-9 km)
- Velocity: 1 nanometer/second (1e-12 km/s)
- Quaternion: ~1e-8 (~0.00001 degrees)

## Workflow Instructions

1. Write plan to `tasks/todo.md` with checkable items
2. For each step, list code that SHOULD NOT be modified
3. Check in before starting — allow discussion
4. Work through items, marking complete as you go
5. After major phases: add tests, commit, ask for review
6. Add review section to `tasks/todo.md` with summary
7. Find and fix root causes — no temporary fixes
8. Minimal, focused changes only
