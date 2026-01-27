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
- [x] Create src/ck/slerp.rs
- [x] Create src/ck/evaluate.rs
- [x] Tests: SLERP, CK type 1/3

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

## Review

### Phase 7 Summary
- **Git LFS**: Configured to track `*.bsp`, `*.bc`, `*.bpc`, `*.hdf5` in `test_data/`. Text files (`.tls`, `.tpc`) remain regular git objects. Removed `/test_data` from `.gitignore`.
- **Examples**: 3 examples ported from muad-dib, adapted to understated's `Ephemeris` API. `kernel_pool` skipped (not in understated's scope).
- **Tests**: 4 new integration tests added — state continuity (velocity ≈ d(position)/dt), position magnitude (Earth ~1 AU, Moon ~384k km), CK quaternion normalization, and multi-body Type 2 queries. All 82 tests pass (65 unit + 17 integration).
- **Notes**: SPK/CK support status table with validation gaps. Full muad-dib API dependency inventory.
- **README**: Experimental status, architecture diagram, supported types, build instructions, validation targets.
