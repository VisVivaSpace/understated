//! Recursive center-body traversal for SPK state evaluation.
//!
//! When requesting the state of body A relative to body B, the SPK file
//! may only contain A relative to some intermediate center C. This module
//! chains through the SPK segment hierarchy to compute arbitrary pairings.
//!
//! The algorithm computes both bodies relative to SSB (the root of the
//! SPK tree), then subtracts: `state(A, B) = state(A, SSB) - state(B, SSB)`.

use crate::error::{Error, Result};
use crate::spk::evaluate_spk;
use crate::state::State;
use crate::types::{EpochTDB, NaifId};
use muad_dib::kernel::SpiceKernel;

const SSB: i32 = 0;

/// Evaluate the state of `target` relative to `center` at `epoch`.
///
/// Searches loaded SPK segments for coverage and chains through
/// intermediate center bodies as needed.
pub fn state_of(
    kernel: &SpiceKernel,
    target: NaifId,
    epoch: EpochTDB,
    center: NaifId,
) -> Result<State> {
    if target == center {
        return Ok(State::new(target, center, 1, [0.0; 3], [0.0; 3]));
    }

    // Compute target relative to SSB
    let target_ssb = state_relative_to_ssb(kernel, target, epoch)?;

    // If center is SSB, we're done
    if center.0 == SSB {
        return Ok(target_ssb);
    }

    // Compute center relative to SSB
    let center_ssb = state_relative_to_ssb(kernel, center, epoch)?;

    // target rel center = (target rel SSB) - (center rel SSB)
    Ok(target_ssb - center_ssb)
}

/// Compute the state of a body relative to the Solar System Barycenter (SSB).
///
/// Follows the SPK segment chain upward: body → segment_center → ... → SSB.
fn state_relative_to_ssb(
    kernel: &SpiceKernel,
    body: NaifId,
    epoch: EpochTDB,
) -> Result<State> {
    if body.0 == SSB {
        return Ok(State::new(NaifId(SSB), NaifId(SSB), 1, [0.0; 3], [0.0; 3]));
    }

    let md_body = muad_dib::types::NaifId(body.0);
    let epoch_f = epoch.0;

    // Find a segment for this body
    let segment = kernel
        .spk_segments_for(md_body)
        .find(|seg| seg.initial_epoch <= epoch_f && epoch_f <= seg.final_epoch)
        .ok_or(Error::NoCoverage {
            body,
            epoch: epoch_f,
        })?;

    // Evaluate the segment
    let view = kernel.spk_view(segment);
    let data = view.data();
    let mut state = evaluate_spk(data, epoch_f)?;

    state.target = NaifId(segment.target_code);
    state.center = NaifId(segment.center_code);
    state.frame = segment.frame_code;

    // If segment center is SSB, done
    if segment.center_code == SSB {
        return Ok(state);
    }

    // Otherwise, chain upward: get segment_center relative to SSB
    let parent = state_relative_to_ssb(kernel, NaifId(segment.center_code), epoch)?;

    // parent = (seg_center rel SSB), state = (body rel seg_center)
    // result = (SSB → seg_center) + (seg_center → body) = (SSB → body)
    Ok(parent + state)
}
