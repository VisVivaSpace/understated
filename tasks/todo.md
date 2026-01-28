# understated - Implementation Progress

## Phase 0: Project Setup
- [x] Copy test_data/ from muad-dib
- [x] Add test_data/ to .gitignore
- [x] Add test-data feature flag

## Phase 1: Foundation
- [x] Update Cargo.toml — add muad-dib path dep, thiserror
- [x] Create src/error.rs
- [x] Create src/types.rs
- [x] Create src/state.rs
- [x] Create src/pointing.rs
- [x] Create src/lib.rs with module declarations and re-exports
- [x] Tests: State arithmetic, Pointing normalization

## Phase 2: Coordinates + Time
- [x] Create src/coord.rs
- [x] Create src/time.rs
- [x] Create src/lsk.rs
- [x] Tests: coordinate round-trips, epoch parsing, format round-trips

## Phase 3: SPK Interpolation
- [x] Create src/spk/mod.rs
- [x] Create src/spk/chebyshev.rs
- [x] Create src/spk/lagrange.rs
- [x] Create src/spk/hermite.rs
- [x] Create src/spk/twobody.rs
- [x] Tests: per-algorithm known-value tests

## Phase 4: CK Interpolation
- [x] Create src/ck/mod.rs
- [x] Create src/ck/evaluate.rs (CK type 1 discrete, type 3 axis-angle interpolation)
- [x] Tests: CK type 1/3

## Phase 5: Ephemeris Facade + Body Chaining
- [x] Create src/spk/chain.rs
- [x] Create src/ephemeris.rs
- [x] Update src/lib.rs with final exports
- [x] Integration tests with real SPK files

## Phase 6: CSPICE Validation
- [x] Add cspice feature flag
- [x] Create tests/cspice_validation.rs
- [x] Validate SPK/CK/time against CSPICE

## Phase 7: Documentation, Examples, and Housekeeping
- [x] Set up Git LFS for binary test data files
- [x] Port examples from muad-dib (query_ephemeris, coordinate_transforms, time_conversion)
- [x] Compile and run all examples
- [x] Add integration tests (state continuity, position magnitude, quaternion normalization, multi-body)
- [x] Create notes/spk_ck_support.md
- [x] Create notes/muad_dib_dependencies.md
- [x] Create README.md

## Phase 8: CK Algorithm Fix + Test Tolerance Audit

### Findings

**CK Type 3 — algorithm mismatch (root cause of test failure):**
- understated originally used quaternion SLERP; now uses CSPICE CKE03 axis-angle algorithm
- CSPICE CKE03 uses rotation-matrix axis-angle interpolation:
  q→R, relative rotation, extract axis-angle, partial rotation, compose
- Both are geodesic interpolation on SO(3), different floating-point paths
- Measured: max 6.274e-10 in rotation matrix elements (midpoint)
- At data records: exact match (0.0). Angular velocity: exact (2.7e-20)

**Time — tolerance too loose:**
- `test_str2et_vs_cspice` uses 1e-6 but actual max diff is 6.2e-11
- `test_et2utc_vs_cspice` — all strings match exactly

**SPK — already at machine precision (no changes needed):**
- Type 2: 7.5e-9 km pos (relative ~7.5e-17), 3.6e-15 km/s vel
- Type 9: 1.8e-12 km, 1.8e-15 km/s
- Type 13: 9.1e-13 km, 1.3e-14 km/s

### Tasks

- [x] Add rotation matrix utilities in `src/rotation.rs`:
      q2m, m2q (Shepperd), raxisa, axisar, mtxm, mxmt
- [x] Rewrite CK Type 3 `evaluate_type3()` to match CKE03 algorithm
- [x] Tighten `test_str2et_vs_cspice` tolerance from 1e-6 → 1e-9
- [x] Tighten `test_ck_vs_cspice` tolerance (expect machine precision)
- [x] Delete diagnostic `tests/tolerance_audit.rs`
- [x] Verify: `cargo test --features cspice,test-data` all pass
- [x] Verify: `cargo clippy` clean

**Phase 8 constraint (completed):** Only rotation.rs and ck/evaluate.rs were
modified. State, SPK, pointing, time, and lsk were left unchanged.

## Review

### Phase 7 Summary
- **Git LFS**: Configured to track `*.bsp`, `*.bc`, `*.bpc`, `*.hdf5` in `test_data/`. Text files (`.tls`, `.tpc`) remain regular git objects. Removed `/test_data` from `.gitignore`.
- **Examples**: 3 examples ported from muad-dib, adapted to understated's `Ephemeris` API. `kernel_pool` skipped (not in understated's scope).
- **Tests**: 4 new integration tests added — state continuity (velocity ≈ d(position)/dt), position magnitude (Earth ~1 AU, Moon ~384k km), CK quaternion normalization, and multi-body Type 2 queries. All 79 tests pass (54 unit + 8 CSPICE validation + 17 integration).
- **Notes**: SPK/CK support status table with validation gaps. Full muad-dib API dependency inventory.
- **README**: Experimental status, architecture diagram, supported types, build instructions, validation targets.
