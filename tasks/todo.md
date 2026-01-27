# understated - Implementation Progress

## Phase 0: Project Setup
- [ ] Copy test_data/ from muad-dib
- [ ] Add test_data/ to .gitignore
- [ ] Add test-data feature flag

## Phase 1: Foundation
- [ ] Update Cargo.toml — add muad-dib path dep, thiserror
- [ ] Create src/error.rs
- [ ] Create src/types.rs
- [ ] Create src/state.rs
- [ ] Create src/pointing.rs
- [ ] Create src/lib.rs with module declarations and re-exports
- [ ] Tests: State arithmetic, Pointing normalization

## Phase 2: Coordinates + Time
- [ ] Create src/coord.rs
- [ ] Create src/time.rs
- [ ] Create src/lsk.rs
- [ ] Tests: coordinate round-trips, epoch parsing, format round-trips

## Phase 3: SPK Interpolation
- [ ] Create src/spk/mod.rs
- [ ] Create src/spk/chebyshev.rs
- [ ] Create src/spk/lagrange.rs
- [ ] Create src/spk/hermite.rs
- [ ] Create src/spk/twobody.rs
- [ ] Tests: per-algorithm known-value tests

## Phase 4: CK Interpolation
- [ ] Create src/ck/mod.rs
- [ ] Create src/ck/slerp.rs
- [ ] Create src/ck/evaluate.rs
- [ ] Tests: SLERP, CK type 1/3

## Phase 5: Ephemeris Facade + Body Chaining
- [ ] Create src/spk/chain.rs
- [ ] Create src/ephemeris.rs
- [ ] Update src/lib.rs with final exports
- [ ] Integration tests with real SPK files

## Phase 6: CSPICE Validation
- [ ] Add cspice feature flag
- [ ] Create tests/cspice_validation.rs
- [ ] Validate SPK/CK/time against CSPICE
