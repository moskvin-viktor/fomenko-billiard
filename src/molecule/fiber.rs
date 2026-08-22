//! Layer 7 of the molecule engine: critical fibers.
//!
//! At a `genus_jump` or `split_merge` vertex the fiber is singular.  The
//! `CriticalFiber` record is the "few trajectories" object: the spine is the
//! finite separatrix graph (the banned corner orbits), the cylinders are the
//! bulk, and `(normalization_genus, n_nodes, arithmetic_genus)` pins the
//! homeomorphism type and its placement between the two regular sides.

use crate::table::Table;
use crate::torus::ConfocalParams;

use super::critical::{classify_transition, TransitionKind};

/// A critical (singular) fiber at a bifurcation value.
#[derive(Clone, Debug)]
pub struct CriticalFiber {
    pub lc: f32,
    /// `"genus_jump"` | `"split_merge"`.
    pub kind: &'static str,
    /// Genus of the lower side.
    pub normalization_genus: u32,
    /// Pinch points (glued point-pairs).
    pub n_nodes: u32,
    /// Must equal the higher side's genus.
    pub arithmetic_genus: u32,
    /// Human description of the separatrix graph.
    pub spine: String,
}

/// Build the critical fiber at a bifurcation value, `None` for non-bifurcations
/// (doc §9).
pub fn critical_fiber(table: &Table, cf: &ConfocalParams, lc: f32) -> Option<CriticalFiber> {
    let t = classify_transition(table, cf, lc, 1e-4);
    match t.kind {
        TransitionKind::GenusJump => {
            let g_lo = t.from.iter().copied().max().unwrap_or(0);
            let g_hi = t.to.iter().copied().max().unwrap_or(0);
            let (g_lo, g_hi) = if g_lo <= g_hi {
                (g_lo, g_hi)
            } else {
                (g_hi, g_lo)
            };
            Some(CriticalFiber {
                lc,
                kind: "genus_jump",
                normalization_genus: g_lo,
                n_nodes: g_hi - g_lo,
                arithmetic_genus: g_hi,
                spine: format!(
                    "figure-eight through the reflex-corner cone point (3-prong), \
                     1 collapsing cylinder; pinched genus-{g_hi} surface = \
                     genus-{g_lo} with {n} node(s)",
                    n = g_hi - g_lo
                ),
            })
        }
        TransitionKind::SplitMerge => Some(CriticalFiber {
            lc,
            kind: "split_merge",
            normalization_genus: 1,
            n_nodes: 1,
            arithmetic_genus: 1,
            spine: "atom-B saddle: two tori joined along a circle through a \
                    hyperbolic 2-orbit"
                .to_string(),
        }),
        _ => None,
    }
}
