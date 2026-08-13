//! The confocal structure of a billiard domain.
//!
//! One `Domain -> (all λ-relevant info)` analysis object.  Building it walks the
//! domain's quadric arcs once and derives everything downstream — the boundary
//! lambdas, the quadrilateral flag, the torus regime, the torus-mapping bounds,
//! and the inside-test — so the many functions that previously re-extracted
//! this per-domain information from a [`Domain`] now share a single source.

use crate::{domain, quadratic, torus::ConfocalParams};
use macroquad::prelude::*;

/// The confocal geometry of a domain, derived once from its quadric arcs.
#[derive(Clone, Copy, Debug)]
pub struct ConfocalStructure {
    /// The confocal family `(a, b)` the domain's quadrics live in.
    pub cf: ConfocalParams,
    /// λ_ell — the outer ellipse boundary (min `λ < b`).
    pub lambda_ell: f32,
    /// λ_hyp — the inner hyperbola boundary (max `b < λ < a`).
    /// `None` if the domain has no hyperbola wall (full-ellipse / L-shape).
    pub lambda_hyp: Option<f32>,
    /// Whether the boundary is exactly a confocal quadrilateral.
    pub is_quadrilateral: bool,
}

impl ConfocalStructure {
    /// Derive the structure from a domain's quadric boundary arcs.
    /// Returns `None` for domains with no confocal arcs (polyline).
    pub fn of_domain(dom: &domain::Domain) -> Option<Self> {
        let quadrics: Vec<&quadratic::ConfocalQuadric> = dom
            .segments
            .iter()
            .filter_map(|s| match s {
                domain::Segment::Quad { curve, .. } => Some(curve),
                _ => None,
            })
            .collect();
        if quadrics.is_empty() {
            return None;
        }

        let cf = ConfocalParams::new(quadrics[0].a_param, quadrics[0].b_param);
        let b = cf.b;
        let mut ell: Option<f32> = None;
        let mut hyp: Option<f32> = None;
        for &curve in &quadrics {
            let lam = curve.lambda;
            if lam < b {
                ell = Some(ell.map_or(lam, |e| e.min(lam)));
            } else if lam > b {
                hyp = Some(hyp.map_or(lam, |h| h.max(lam)));
            }
        }

        Some(Self {
            cf,
            lambda_ell: ell.unwrap_or(0.0),
            lambda_hyp: hyp,
            is_quadrilateral: quadrics.len() == 4 && ell.is_some() && hyp.is_some(),
        })
    }

    /// Torus-mapping bounds `(lam_wall, beta)`: `lam_wall` = λ_ell (or 0 for the
    /// full-ellipse fallback), `beta` = the hyperbola wall, `None` = full ellipse.
    /// Non-quadrilaterals fall back to the full-ellipse convention (`(0, None)`),
    /// matching how the torus map treats them.
    pub fn torus_bounds(&self) -> (f32, Option<f32>) {
        if self.is_quadrilateral {
            (self.lambda_ell, self.lambda_hyp)
        } else {
            (0.0, None)
        }
    }

    /// The torus regime for a caustic value `lam` (Λ for confocal).  Replaces the
    /// old `TorusRegime::from_value(domain, lam, beta)`.
    pub fn regime(&self, lam: f32) -> TorusRegime {
        if lam >= self.cf.b {
            // Hyperbolic caustic: single torus.
            TorusRegime::Single
        } else if self.is_quadrilateral {
            // Case C: ellipse caustic split by sign(y).
            TorusRegime::SplitByY
        } else {
            // Case A: ellipse caustic on full ellipse, split by sign(L).
            TorusRegime::SplitByL
        }
    }

    /// Membership test for a confocal quadrilateral.  A point is inside when it
    /// lies within the boundary ellipse and between the two hyperbola sheets.
    /// `slack` shrinks the region slightly inward.  Only valid for exact
    /// quadrilaterals; for non-quadrilaterals use [`domain::Domain::contains`].
    pub fn contains(&self, p: Vec2, slack: f32) -> bool {
        let Some(lambda_hyp) = self.lambda_hyp else {
            return false;
        };
        let a = self.cf.a;
        let b = self.cf.b;
        let lambda_ell = self.lambda_ell;
        // Ellipse: Q_λ_ell(p) < 0 means inside.
        let e_val = (b - lambda_ell) * p.x * p.x + (a - lambda_ell) * p.y * p.y
            - (a - lambda_ell) * (b - lambda_ell);
        if e_val >= -slack {
            return false;
        }
        let h_val = (b - lambda_hyp) * p.x * p.x + (a - lambda_hyp) * p.y * p.y
            - (a - lambda_hyp) * (b - lambda_hyp);
        h_val > slack
    }
}

/// How trajectories sit on the torus at a caustic value and domain shape.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum TorusRegime {
    /// Hyperbola caustic (Λ > B): Case B / C-hyperbola, one torus.
    Single,
    /// Ellipse caustic on a confocal quadrilateral (Case C): two tori, split by
    /// `sign(y)` into upper / lower regions.
    SplitByY,
    /// Ellipse caustic on a full ellipse (Case A): two tori, split by the sign
    /// of the angular momentum `x·vy − y·vx`.
    SplitByL,
}

impl TorusRegime {
    /// Number of disconnected tori for this regime.
    pub fn n_tori(self) -> u32 {
        match self {
            Self::Single => 1,
            Self::SplitByY | Self::SplitByL => 2,
        }
    }

    /// Classifier of the caustic: `"ellipse"` (Λ < B) or `"hyperbola"` (Λ > B).
    pub fn caustic_label(self) -> &'static str {
        match self {
            Self::SplitByY | Self::SplitByL => "ellipse",
            Self::Single => "hyperbola",
        }
    }

    /// Torus index for a phase-space point `(x, y, vx, vy)` under this regime.
    pub fn index_of(self, x: f32, y: f32, vx: f32, vy: f32) -> u32 {
        match self {
            Self::Single => 0,
            Self::SplitByY => {
                if y >= 0.0 {
                    0
                } else {
                    1
                }
            }
            Self::SplitByL => {
                if x * vy - y * vx > 0.0 {
                    0
                } else {
                    1
                }
            }
        }
    }
}

/// Sampling policy for [`caustic_starts`].
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum CausticSampling {
    /// One point per connected component of caustic ∩ domain.
    Sparse,
    /// `per_component` points spread along each component, both ±velocity
    /// directions (so the union fills the 2D torus).
    Dense,
}

impl CausticSampling {
    /// Justify `n_samples` and the run retirement for each policy.
    fn n_samples(self) -> usize {
        match self {
            Self::Sparse => 2048,
            Self::Dense => 4096,
        }
    }
}

/// Pick start points on the caustic inside the domain.
///
/// - [`CausticSampling::Sparse`] → one `(p, v)` per connected component.
/// - [`CausticSampling::Dense`] → spread `per_component` samples per component,
///   both directions of motion, so the union sweeps out the full 2D torus.
///
/// The inside-membership test is analytic when `structure.is_quadrilateral`
/// (via the structure), and ray-cast via `dom.contains` otherwise.
pub fn caustic_starts(
    structure: &ConfocalStructure,
    dom: &domain::Domain,
    lam: f32,
    sampling: CausticSampling,
    per_component: usize,
) -> Vec<(Vec2, Vec2)> {
    let quad = quadratic::confocal(structure.cf, lam);
    let n_samples = sampling.n_samples();
    let pts = quad.sample_boundary(n_samples);

    // Extra slack for component detection on non-quadrilaterals where the
    // containing ray-cast is unreliable near the boundary.
    let detect_slack = if structure.is_quadrilateral {
        0.03
    } else {
        0.06
    };
    let mut inside: Vec<bool> = Vec::with_capacity(n_samples);
    for &p in &pts {
        let in_domain = if structure.is_quadrilateral {
            structure.contains(p, detect_slack)
        } else {
            dom.contains(p)
        };
        inside.push(in_domain);
    }
    let n = pts.len();
    if !inside.iter().any(|&x| x) {
        return vec![];
    }

    let start = (0..n)
        .find(|&i| inside[i] && !inside[(i + n - 1) % n])
        .unwrap();

    let mut starts = Vec::new();
    let mut i = start;
    loop {
        while !inside[i] {
            i = (i + 1) % n;
            if i == start {
                return starts;
            }
        }
        let run_start = i;
        while inside[i] {
            i = (i + 1) % n;
            if i == start {
                break;
            }
        }
        let run_end = (i + n - 1) % n;
        let run_len = (run_end + n - run_start) % n + 1;

        match sampling {
            CausticSampling::Sparse => {
                // Ignore degenerate slivers (samples flickering on the boundary).
                if run_len >= 8 {
                    let mid = (run_start + run_len / 2) % n;
                    if let Some((p, vel)) = snap_velocity(&quad, pts[mid], structure, dom) {
                        starts.push((p, vel));
                    }
                }
            }
            CausticSampling::Dense => {
                // Spread `per_component` points, skipping boundary-adjacent samples.
                let margin = 8;
                let usable = run_len.saturating_sub(2 * margin);
                if usable > 0 {
                    let count = per_component.min(usable);
                    for k in 0..count {
                        let idx = (run_start
                            + margin
                            + (usable as f32 * k as f32 / count as f32) as usize)
                            % n;
                        if let Some((p, vel)) = snap_velocity(&quad, pts[idx], structure, dom) {
                            // Add BOTH directions of motion.
                            starts.push((p, vel));
                            starts.push((p, -vel));
                        }
                    }
                }
            }
        }
        if i == start {
            return starts;
        }
    }
}

/// Newton-project a point onto the quadric, then attach the caustic velocity.
/// Returns `None` if the projection lands outside the domain.
fn snap_velocity(
    quad: &quadratic::ConfocalQuadric,
    p: Vec2,
    structure: &ConfocalStructure,
    dom: &domain::Domain,
) -> Option<(Vec2, Vec2)> {
    // Snap onto the caustic curve along the gradient so the tangent velocity is
    // exact.
    let mut q = p;
    for _ in 0..8 {
        let v = quad.eval(q);
        if v.abs() < 1e-9 {
            break;
        }
        let g = quad.grad(q);
        let g2 = g.length_squared();
        if g2 < 1e-12 {
            break;
        }
        q -= (v / g2) * g;
    }
    // The snapped point must be genuinely inside (no slack): the projection can
    // push a sample onto the boundary.
    let genuinely_inside = if structure.is_quadrilateral {
        structure.contains(q, 0.0)
    } else {
        dom.contains(q)
    };
    if !genuinely_inside {
        return None;
    }
    let vel = quadratic::ConfocalQuadric::velocity_from_caustic(
        structure.cf,
        q,
        quad.lambda,
        vec2(0.0, 1.0),
    );
    Some((q, vel))
}
