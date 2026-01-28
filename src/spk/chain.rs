//! Recursive center-body traversal for SPK state evaluation.
//!
//! When requesting the state of body A relative to body B, the SPK file
//! may only contain A relative to some intermediate center C. This module
//! chains through the SPK segment hierarchy to compute arbitrary pairings.
//!
//! Strategy:
//! 1. If the target's segment center matches the requested center, return directly.
//! 2. Otherwise, find the nearest common ancestor in the SPK tree and
//!    chain both bodies to that ancestor, then subtract:
//!    `state(A, B) = state(A, ancestor) - state(B, ancestor)`.
//!    This minimizes floating-point error by avoiding unnecessary SSB round-trips.

use std::collections::HashSet;

use crate::error::{Error, Result};
use crate::spk::evaluate_spk;
use crate::state::State;
use crate::types::{EpochTDB, FrameId, NaifId};
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
        return Ok(State::new(target, center, FrameId::J2000, [0.0; 3], [0.0; 3]));
    }

    // First, try the direct path: evaluate target's segment and check
    // if its center matches the requested center.
    let state = evaluate_body(kernel, target, epoch)?;

    if state.center == center {
        return Ok(state);
    }

    // If the segment's center is the target we want (reversed query),
    // just negate.
    if state.target == center {
        return Ok(-state);
    }

    // Find common ancestor to minimize chaining depth and
    // accumulation of floating-point errors.
    let ancestor = find_common_ancestor(kernel, target, center, epoch);

    let target_anc = chain_to(kernel, target, epoch, ancestor)?;

    if center == ancestor {
        return Ok(target_anc);
    }

    let center_anc = chain_to(kernel, center, epoch, ancestor)?;

    // target rel center = (target rel ancestor) - (center rel ancestor)
    Ok(target_anc - center_anc)
}

/// Evaluate a single segment for the given body and apply context.
fn evaluate_body(kernel: &SpiceKernel, body: NaifId, epoch: EpochTDB) -> Result<State> {
    let md_body = muad_dib::types::NaifId(body.0);
    let epoch_f = epoch.0;

    let segment = kernel
        .spk_segments_for(md_body)
        .find(|seg| seg.initial_epoch <= epoch_f && epoch_f <= seg.final_epoch)
        .ok_or(Error::NoCoverage {
            body,
            epoch: epoch_f,
        })?;

    let view = kernel.spk_view(segment);
    let data = view.data();
    let mut state = evaluate_spk(data, epoch_f)?;

    state.target = NaifId::from(segment.target_code);
    state.center = NaifId::from(segment.center_code);
    state.frame = FrameId(segment.frame_code.0);

    Ok(state)
}

/// Walk the center chain for a body, returning the list of center body IDs.
fn center_chain(kernel: &SpiceKernel, body: NaifId, epoch: f64) -> Vec<NaifId> {
    let mut chain = vec![body];
    let mut current = body;
    for _ in 0..20 {
        // Safety limit to prevent infinite loops
        let md_body = muad_dib::types::NaifId(current.0);
        let seg = kernel
            .spk_segments_for(md_body)
            .find(|seg| seg.initial_epoch <= epoch && epoch <= seg.final_epoch);
        match seg {
            Some(s) => {
                let center = NaifId::from(s.center_code);
                chain.push(center);
                if center.0 == SSB {
                    break;
                }
                current = center;
            }
            None => break,
        }
    }
    chain
}

/// Find the nearest common ancestor of two bodies in the SPK segment tree.
/// Falls back to SSB if no closer ancestor is found.
fn find_common_ancestor(
    kernel: &SpiceKernel,
    target: NaifId,
    center: NaifId,
    epoch: EpochTDB,
) -> NaifId {
    let target_chain = center_chain(kernel, target, epoch.0);
    let center_set: HashSet<NaifId> = center_chain(kernel, center, epoch.0).into_iter().collect();

    // Find the first body in target's chain that also appears in center's chain.
    for &body in &target_chain {
        if center_set.contains(&body) {
            return body;
        }
    }

    NaifId(SSB)
}

/// Compute the state of a body relative to `ancestor` by following the segment chain.
fn chain_to(
    kernel: &SpiceKernel,
    body: NaifId,
    epoch: EpochTDB,
    ancestor: NaifId,
) -> Result<State> {
    if body == ancestor {
        return Ok(State::new(body, ancestor, FrameId::J2000, [0.0; 3], [0.0; 3]));
    }

    let state = evaluate_body(kernel, body, epoch)?;

    if state.center == ancestor {
        return Ok(state);
    }

    // Chain upward: get segment_center relative to ancestor
    let parent = chain_to(kernel, state.center, epoch, ancestor)?;

    // parent = (seg_center rel ancestor), state = (body rel seg_center)
    // result = (ancestor → seg_center) + (seg_center → body) = (ancestor → body)
    Ok(parent + state)
}
