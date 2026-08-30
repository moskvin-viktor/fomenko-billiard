//! Level → Manifold: the single table `λ → renderable manifold`.
//!
//! Implements `docs/level_to_manifold.md`.  The app asks "what is the manifold
//! at λ?" and gets back a flat object it can draw — a 2D torus point cloud, a
//! closed 1D curve (one circle per collapsed torus), a flat pseudo-integrable
//! surface, or nothing.  No cascade of if-else branches at the call site.
//!
//! The key rule (see `docs/degenerate_embedding.md`): the shape of a degenerate
//! layer is **read off the torus mapping**, never hardcoded.  Each sliding
//! start point is pushed through [`crate::torus::to_torus`]; at a degenerate
//! level one torus angle collapses to a constant and the other circulates, so
//! each torus degrades to a circle.  Circles are found by clustering the mapped
//! samples by `torus_index` (and, within a torus, by the collapsed angle — two
//! walls sliding on the same torus give two circles), then ordering each
//! cluster by its circulating angle.  At the separatrix `λ = b` the map returns
//! its sentinel; the sample is re-mapped a hair off the focal segment on both
//! sides, which lands on the two touching tori — the figure-eight falls out.

use crate::bifurcation::{classify, PhaseManifold};
use crate::confocal::ConfocalStructure;
use crate::domain::Domain;
use crate::phase3d::PhasePoint;
use crate::pseudo;
use macroquad::prelude::*;

use std::f32::consts::{PI, TAU};

/// Sampling densities for [`build_manifold`].
pub struct SampleConfig {
    /// Bounces per dense trajectory (2D torus / flat surface fill).
    pub max_steps: usize,
    /// Interior samples per trajectory segment.
    pub per_seg: usize,
    /// Caustic start points per connected component (dense fill).
    pub per_component: usize,
    /// Samples spread along a degenerate border piece (1D curve).
    pub curve_samples: usize,
}

impl Default for SampleConfig {
    fn default() -> Self {
        Self {
            max_steps: 200,
            per_seg: 8,
            per_component: 24,
            curve_samples: 256,
        }
    }
}

/// The phase manifold at a caustic level, as a directly renderable object.
pub enum Manifold {
    /// Generic integrable level: a 2D Liouville torus surface (1 or 2 tori),
    /// sampled as dense phase trajectories.  The `torus_index` baked into each
    /// point tells the embedder which lobe it belongs to.
    Torus2D {
        trajectories: Vec<Vec<PhasePoint>>,
    },
    /// Degenerate critical level: the manifold collapsed to closed 1D curves,
    /// one per collapsed torus (two touching circles at the separatrix — the
    /// figure-eight).  Each circle is an ordered loop of torus-angle points,
    /// embedded at draw time via `torus_embed` with the same scale/offset as
    /// the regular tori.
    Curve1D {
        circles: Vec<Vec<PhasePoint>>,
        /// The circles share a point (the hyperbolic pinch orbit): true exactly
        /// when the torus map hit its separatrix sentinel and each sliding
        /// sample re-mapped onto *both* touching tori.  Read off the mapping,
        /// never off the level value.
        pinched: bool,
    },
    /// Pseudo-integrable table (L/T/Z): a flat torus or genus-≥2 surface.
    FlatSurface {
        trajectories: Vec<Vec<pseudo::FlatPhasePoint>>,
        level: pseudo::Level,
        /// `Some` at a molecule critical value: the level is the one-sided
        /// limit and the renderer should pinch/collapse the embedding.
        singular: Option<crate::bifurcation::SingularInfo>,
    },
    /// Nothing reachable at this level.
    Empty,
}

impl Manifold {
    /// Number of connected components of a `Curve1D` (0 otherwise-shaped).
    pub fn n_components(&self) -> u32 {
        match self {
            Manifold::Curve1D { circles, .. } => circles.len() as u32,
            _ => 0,
        }
    }
}

/// The level → manifold table: one entry point for every confocal level.
///
/// Classification is delegated to [`crate::bifurcation::classify`]; this
/// function turns the classification into a *sampled, renderable* manifold:
///
/// - pseudo-integrable → [`Manifold::FlatSurface`] (flat-chart trajectories);
/// - degenerate critical layer → [`Manifold::Curve1D`] (circles read off the
///   torus mapping of the sliding starts);
/// - generic integrable → [`Manifold::Torus2D`] (dense caustic sampling);
/// - unreachable → [`Manifold::Empty`].
///
/// Non-confocal (polyline) domains have no confocal structure and no manifold
/// table; they return [`Manifold::Empty`].
pub fn build_manifold(domain: &Domain, lam: f32, config: &SampleConfig) -> Manifold {
    let Some(structure) = ConfocalStructure::of_domain(domain) else {
        return Manifold::Empty;
    };
    let cf = structure.cf;

    match classify(domain, &cf, lam) {
        PhaseManifold::Forbidden => Manifold::Empty,

        PhaseManifold::Flat { level, singular } => {
            let level = *level;
            match &level {
                pseudo::Level::Torus { .. } | pseudo::Level::GenusSurface { .. } => {
                    // At a singular value, sample at the same one-sided λ the
                    // level was classified at (mirrors the λ_hyp nudge below):
                    // the exact-critical caustic is tangentially degenerate.
                    let lam_s = singular.as_ref().map_or(lam, |s| s.side_lam);
                    let starts = crate::dense_caustic_starts(domain, lam_s, config.per_component);
                    let trajectories = starts
                        .iter()
                        .map(|&(p, v)| {
                            pseudo::sample_flat_trajectory(
                                domain,
                                p,
                                v,
                                config.max_steps,
                                config.per_seg,
                                &level,
                            )
                        })
                        .collect();
                    Manifold::FlatSurface {
                        trajectories,
                        level,
                        singular,
                    }
                }
                _ => Manifold::Empty,
            }
        }

        PhaseManifold::Degenerate { .. } => {
            let (circles, pinched) =
                degenerate_circles(&structure, domain, lam, config.curve_samples);
            if circles.is_empty() {
                // The border piece does not intersect the domain at all.
                Manifold::Empty
            } else {
                Manifold::Curve1D { circles, pinched }
            }
        }

        PhaseManifold::Tori { .. } => {
            // The hyperbola wall is a *regular* level, but sampling exactly
            // there is numerically degenerate: the caustic coincides with the
            // wall, so every start grazes it tangentially.  Sample a hair
            // inside the band instead — the torus is continuous in λ, so the
            // render is indistinguishable.
            let lam = match structure.lambda_hyp {
                Some(h) if (lam - h).abs() < 1e-3 => h + 1e-3,
                _ => lam,
            };
            let bounds = structure.torus_bounds();
            let starts = crate::dense_caustic_starts(domain, lam, config.per_component);
            let trajectories = starts
                .iter()
                .map(|&(p, v)| {
                    crate::phase3d::sample_trajectory_phase_dense(
                        domain,
                        p,
                        v,
                        config.max_steps,
                        config.per_seg,
                        bounds,
                    )
                })
                .collect();
            Manifold::Torus2D { trajectories }
        }
    }
}

// ----------------------------------------------------------------
// Degenerate layers: circles read off the torus mapping
// ----------------------------------------------------------------

/// Perpendicular nudge used to re-map separatrix samples (where `to_torus`
/// returns its sentinel): moves the sample off the focal segment onto each of
/// the two touching tori.  `lc` shifts from `b` to `b − NUDGE²`, safely past
/// the map's `sep_eps` while staying on the elliptic side.
const NUDGE: f32 = 1e-3;

/// Minimum samples for a cluster to count as a circle (smaller = noise).
const MIN_CIRCLE: usize = 8;

/// Circular gap (radians) that splits two curves sharing a torus (e.g. the
/// left and right hyperbola walls, pinned at collapsed θ₂ = 0 and π).
const SPLIT_GAP: f32 = 0.8;

/// Build the closed 1D curves of a degenerate caustic level.
///
/// Every sliding start point (border piece + tangent velocity, from
/// [`crate::confocal::critical_caustic_starts`]) is mapped through `to_torus`
/// in both directions of motion.  A collapsed torus angle comes back constant
/// (or non-finite from a zero-width libration — sanitised to a constant), the
/// surviving angle circulates.  Clustering by `torus_index`, then by the
/// collapsed angle, and ordering each cluster by the circulating angle yields
/// the circles — their number and connectivity fall out of the mapping.
///
/// Returns `(circles, pinched)`: `pinched` is true when the samples sat on the
/// separatrix itself (the map's sentinel fired and each sample landed on both
/// touching tori after the perpendicular nudge) — i.e. the circles meet at the
/// hyperbolic pinch orbit and form a figure-eight.
fn degenerate_circles(
    structure: &ConfocalStructure,
    domain: &Domain,
    lam: f32,
    n: usize,
) -> (Vec<Vec<PhasePoint>>, bool) {
    let starts = match crate::confocal::critical_caustic_starts(structure, domain, lam, n) {
        Some(s) => s,
        None => return (Vec::new(), false),
    };
    if starts.is_empty() {
        return (Vec::new(), false);
    }

    let cf = structure.cf;
    let bounds = structure.torus_bounds();
    let mut cache = crate::torus::TorusCache::default();

    // Map every start in both directions of motion (forward + return pass of
    // the sliding orbit — together they close the loop).
    let mut mapped: Vec<PhasePoint> = Vec::with_capacity(2 * starts.len());
    let mut n_samples = 0usize;
    let mut sentinel_hits = 0usize;
    for &(p, v) in &starts {
        for dir in [v, -v] {
            n_samples += 1;
            if map_degenerate_sample(p, dir, &cf, bounds, &mut cache, &mut mapped) {
                sentinel_hits += 1;
            }
        }
    }
    // Pinched when the sliding orbit genuinely lives on the separatrix: the
    // sentinel fired for (essentially) all samples, so every sample fanned out
    // onto both touching tori.
    let pinched = n_samples > 0 && 2 * sentinel_hits >= n_samples;

    // Cluster by torus index.
    let mut by_torus: Vec<(u32, Vec<PhasePoint>)> = Vec::new();
    for pt in mapped {
        match by_torus.iter_mut().find(|(idx, _)| *idx == pt.torus_index) {
            Some((_, v)) => v.push(pt),
            None => by_torus.push((pt.torus_index, vec![pt])),
        }
    }
    by_torus.sort_by_key(|(idx, _)| *idx);

    let mut circles = Vec::new();
    for (_, group) in by_torus {
        if group.len() < MIN_CIRCLE {
            continue;
        }

        // Which angle circulates?  The circulating angle covers the circle
        // roughly uniformly (small maximum circular gap); the collapsed angle
        // sits at one or a few constants (large maximum gap).
        let gap1 = max_circular_gap(group.iter().map(|p| p.theta1));
        let gap2 = max_circular_gap(group.iter().map(|p| p.theta2));
        let circ_is_theta1 = gap1 < gap2;

        // Two disjoint sliding orbits can share a torus (left + right wall)
        // at distinct collapsed-angle constants: split them apart.
        for mut sub in split_by_collapsed(group, !circ_is_theta1) {
            if sub.len() < MIN_CIRCLE {
                continue;
            }
            // Order the loop by the circulating angle.
            sub.sort_by(|a, b| {
                let (ka, kb) = if circ_is_theta1 {
                    (a.theta1, b.theta1)
                } else {
                    (a.theta2, b.theta2)
                };
                ka.total_cmp(&kb)
            });
            circles.push(sub);
        }
    }
    let pinched = pinched && circles.len() >= 2;
    (circles, pinched)
}

/// Map one sliding sample through the torus map, sanitising the collapsed
/// (non-finite) angle to a constant.  At the separatrix the map returns its
/// sentinel (`λc = b` exactly); re-map the sample nudged perpendicular to the
/// motion on *both* sides, which lands one sample on each touching torus.
///
/// Returns `true` iff the sentinel fired (the sample sat on the separatrix).
fn map_degenerate_sample(
    p: Vec2,
    v: Vec2,
    cf: &crate::torus::ConfocalParams,
    bounds: (f32, Option<f32>),
    cache: &mut crate::torus::TorusCache,
    out: &mut Vec<PhasePoint>,
) -> bool {
    if let Some(pt) = map_one(p, v, cf, bounds, cache) {
        out.push(pt);
        return false;
    }
    // Separatrix sentinel: nudge perpendicular to the motion, both sides.
    let perp = vec2(-v.y, v.x).normalize_or_zero() * NUDGE;
    for q in [p + perp, p - perp] {
        if let Some(pt) = map_one(q, v, cf, bounds, cache) {
            out.push(pt);
        }
    }
    true
}

/// One `to_torus` call → a sanitised [`PhasePoint`]; `None` on the sentinel.
fn map_one(
    p: Vec2,
    v: Vec2,
    cf: &crate::torus::ConfocalParams,
    bounds: (f32, Option<f32>),
    cache: &mut crate::torus::TorusCache,
) -> Option<PhasePoint> {
    let sample = crate::torus::PhaseSample::new(p.x, p.y, v.x, v.y);
    let mut params = crate::torus::TorusParams {
        confocal: cf,
        lam_wall: bounds.0,
        beta: bounds.1,
        sep_eps: 1e-9,
        cache,
    };
    let (_th1, th2, idx) = crate::torus::to_torus(&sample, &mut params);
    if idx == u32::MAX {
        return None;
    }
    // On a degenerate layer the collapsed degree of freedom is θ₁ by
    // convention (see `torus::map`).  Its libration has (near-)zero width, so
    // the raw phase is either non-finite or pure noise (π·w/w_full with a
    // vanishing w_full) — pin it to the constant 0 rather than trusting it.
    let th1 = 0.0;
    let th2 = if th2.is_finite() { th2 } else { 0.0 };
    Some(PhasePoint {
        x: p.x,
        y: p.y,
        theta: v.y.atan2(v.x) / PI,
        theta1: th1,
        theta2: th2,
        torus_index: idx,
    })
}

/// Largest circular gap between consecutive angle values in `[0, 2π)`.
/// Empty / single-value input counts as fully collapsed (gap = 2π).
fn max_circular_gap(vals: impl Iterator<Item = f32>) -> f32 {
    let mut v: Vec<f32> = vals.map(|x| x.rem_euclid(TAU)).collect();
    if v.len() < 2 {
        return TAU;
    }
    v.sort_by(f32::total_cmp);
    let mut gap = TAU - (v[v.len() - 1] - v[0]); // wrap-around gap
    for w in v.windows(2) {
        gap = gap.max(w[1] - w[0]);
    }
    gap
}

/// Split a torus cluster into sub-clusters by its *collapsed* angle: sort by
/// that angle and cut at circular gaps wider than [`SPLIT_GAP`].  A single
/// constant (possibly wrapped across 0/2π) stays one cluster; two constants
/// (e.g. θ₂ = 0 and π from the two hyperbola walls) split into two.
fn split_by_collapsed(mut pts: Vec<PhasePoint>, collapsed_is_theta1: bool) -> Vec<Vec<PhasePoint>> {
    let key = |p: &PhasePoint| {
        if collapsed_is_theta1 {
            p.theta1.rem_euclid(TAU)
        } else {
            p.theta2.rem_euclid(TAU)
        }
    };
    pts.sort_by(|a, b| key(a).total_cmp(&key(b)));
    let n = pts.len();

    // Group starts: index i such that the circular gap before pts[i] > SPLIT_GAP.
    let mut cuts: Vec<usize> = Vec::new();
    for i in 0..n {
        let prev = key(&pts[(i + n - 1) % n]);
        let cur = key(&pts[i]) + if i == 0 { TAU } else { 0.0 };
        if cur - prev > SPLIT_GAP {
            cuts.push(i);
        }
    }
    if cuts.is_empty() {
        return vec![pts];
    }

    let mut out = Vec::with_capacity(cuts.len());
    for (k, &start) in cuts.iter().enumerate() {
        let end = cuts[(k + 1) % cuts.len()];
        let mut group = Vec::new();
        let mut i = start;
        loop {
            group.push(pts[i]);
            i = (i + 1) % n;
            if i == end {
                break;
            }
        }
        out.push(group);
    }
    out
}
