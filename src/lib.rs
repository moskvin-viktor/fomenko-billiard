pub mod domain;
pub mod phase3d;
pub mod presets;
pub mod quadratic;
pub mod render;
pub mod torus;
pub mod torus_render;

pub const A: f32 = 4.0;
pub const B: f32 = 1.0;

/// How trajectories sit on the torus at the current Λ (caustic type) and
/// domain shape.  This is the single source of truth for "how many tori" and
/// "which observable splits them", so the renderer's `torus_index` rule and
/// the app's per-torus highlight picker agree by construction.
///
/// Frames the Jacobi–Moser cases of `docs/thorus_params.md`:
///
/// * [`TorusRegime::Single`] — a hyperbola caustic (Λ > B): Case B / C-hyperbola,
///   one connected torus, no split.
/// * [`TorusRegime::SplitByY`] — an ellipse caustic (Λ < B) on a confocal
///   quadrilateral (Case C): upper / lower region, split by `sign(y)`.
/// * [`TorusRegime::SplitByL`] — an ellipse caustic on a full ellipse (Case A,
///   e.g. the L-shape fallback with no hyperbola wall): split by the sign of
///   the angular momentum `x·vy − y·vx`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum TorusRegime {
    Single,
    SplitByY,
    SplitByL,
}

impl TorusRegime {
    /// Classify given `lam` (Λ for confocal, θ/π for polyline) for a domain.
    ///
    /// `beta` is the hyperbola-wall parameter from [`torus_bounds`]; `Some` means
    /// a confocal quadrilateral (Case C), `None` a full ellipse (Case A/B).
    pub fn from_value(is_confocal: bool, lam: f32, beta: Option<f32>) -> Self {
        if !is_confocal || lam >= B {
            // Polyline domains and hyperbola caustics are always one torus.
            Self::Single
        } else if beta.is_some() {
            // Case C: ellipse caustic on a confocal quadrilateral split by sign(y).
            Self::SplitByY
        } else {
            // Case A: ellipse caustic on a full ellipse split by sign(L).
            Self::SplitByL
        }
    }

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

    /// Torus index for a phase-space point `(x, y, vx, vy)` under this regime,
    /// mirroring the observable split used by a given embedding.
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

use macroquad::prelude::*;

/// Pick start points for a given domain at a given value of the second integral.
///
/// For confocal billiards (`is_confocal = true`), `lam` is the caustic parameter Λ.
/// We sample the caustic curve Q_Λ = 0 and pick one point in each connected
/// component of caustic ∩ domain.
///
/// For polyline billiards (`is_confocal = false`), `lam` is the normalised angle
/// θ/π ∈ [-1, 1], and `center` is a fixed interior point.  We return a single
/// trajectory start at `center` with velocity direction `θ = lam·π`.
pub fn get_start_points(
    a: f32,
    b: f32,
    lam: f32,
    domain: &domain::Domain,
    is_confocal: bool,
    center: Vec2,
) -> Vec<(Vec2, Vec2)> {
    if is_confocal {
        start_points_on_caustic(a, b, lam, domain)
    } else {
        let angle = lam * std::f32::consts::PI; // lam is θ/π ∈ [-1, 1]
        let v = vec2(angle.cos(), angle.sin());
        vec![(center, v)]
    }
}

/// Pick start points on the caustic inside the domain.
/// Returns one (position, velocity) for each connected component of
/// caustic ∩ domain. Ellipse caustics can have up to 4 components
/// (top-left, top-right, bottom-left, bottom-right arcs); hyperbola
/// caustics up to 2 (left and right branches).
pub fn start_points_on_caustic(
    a: f32,
    b: f32,
    lam: f32,
    domain: &domain::Domain,
) -> Vec<(Vec2, Vec2)> {
    let hint_dir = vec2(0.0, 1.0);
    let quad = quadratic::confocal(a, b, lam);
    let n_samples = 2048;
    let pts = quad.sample_boundary(n_samples);

    // Classify each sample. Prefer the exact confocal structure when the
    // domain is a confocal quadrilateral: it is robust at the boundary,
    // unlike ray casting. A small inward slack keeps only points that
    // are robustly interior, so thin lens-shaped slivers that merely
    // touch the boundary don't become spurious components.
    let bounds = confocal_bounds(domain);
    // Slack for component detection: keep only robustly-interior samples.
    // For non-quadrilaterals (L-shapes), use larger slack since ray casting
    // can be unreliable near the complex boundary.
    let slack = if bounds.is_some() { 0.03 } else { 0.06 };
    let inside: Vec<bool> = pts
        .iter()
        .map(|&p| match bounds {
            Some((le, lh)) => confocal_inside(a, b, p, le, lh, slack),
            None => domain.contains(p),
        })
        .collect();
    let n = pts.len();
    if !inside.iter().any(|&x| x) {
        return vec![];
    }

    // Rotate indices so we start inside a component.
    let start = (0..n)
        .find(|&i| inside[i] && !inside[(i + n - 1) % n])
        .unwrap();

    let mut starts = Vec::new();
    let mut i = start;
    loop {
        // Skip until the beginning of the next inside-run.
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
        let run_end = (i + n - 1) % n; // last inside index of this run
        let run_len = (run_end + n - run_start) % n + 1;

        // Ignore degenerate slivers: near intersections, `contains`
        // (ray casting) can flicker for a few samples and create
        // spurious micro-components right on the domain boundary.
        let min_run = 8;
        if run_len >= min_run {
            // Take the mid-sample: guaranteed interior by construction,
            // and away from the unreliable boundary-adjacent samples.
            let mid = (run_start + run_len / 2) % n;
            // Snap onto the caustic curve along the gradient so the
            // tangent velocity is exact.
            let p = project_to_curve(&quad, pts[mid]);
            // The snapped point must be genuinely inside (no slack): the
            // projection can push a mid-sample onto the boundary.
            let genuinely_inside = match bounds {
                Some((le, lh)) => confocal_inside(a, b, p, le, lh, 0.0),
                None => domain.contains(p),
            };
            if !genuinely_inside {
                continue;
            }
            let vel = quadratic::ConfocalQuadric::velocity_from_caustic(a, b, p, lam, hint_dir);
            starts.push((p, vel));
        }
        if i == start {
            return starts;
        }
    }
}

/// Densely sample start points along the caustic inside the domain.
///
/// Unlike [`start_points_on_caustic`] (which returns one point per connected
/// component), this returns *many* points spread along each component.  This is
/// used to fill out the 2D Liouville torus in the 3D phase-space view: each
/// start point traces a 1D curve, and the union of many such curves at the same
/// Λ sweeps out the full 2D invariant torus.
///
/// In action-angle coordinates the torus is the flat product S¹ × S¹; in the
/// (x, y, θ) embedding used here it appears as a warped 2D surface.
pub fn dense_caustic_starts(
    a: f32,
    b: f32,
    lam: f32,
    domain: &domain::Domain,
    per_component: usize,
) -> Vec<(Vec2, Vec2)> {
    let hint_dir = vec2(0.0, 1.0);
    let quad = quadratic::confocal(a, b, lam);
    let n_samples = 4096;
    let pts = quad.sample_boundary(n_samples);

    let bounds = confocal_bounds(domain);
    let slack = if bounds.is_some() { 0.03 } else { 0.06 };
    let inside: Vec<bool> = pts
        .iter()
        .map(|&p| match bounds {
            Some((le, lh)) => confocal_inside(a, b, p, le, lh, slack),
            None => domain.contains(p),
        })
        .collect();
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

        // Spread `per_component` points along this run, skipping the
        // unreliable boundary-adjacent samples at each end.
        let margin = 8;
        let usable = run_len.saturating_sub(2 * margin);
        if usable > 0 {
            let count = per_component.min(usable);
            for k in 0..count {
                let idx =
                    (run_start + margin + (usable as f32 * k as f32 / count as f32) as usize) % n;
                let p = project_to_curve(&quad, pts[idx]);
                let genuinely_inside = match bounds {
                    Some((le, lh)) => confocal_inside(a, b, p, le, lh, 0.0),
                    None => domain.contains(p),
                };
                if !genuinely_inside {
                    continue;
                }
                let vel = quadratic::ConfocalQuadric::velocity_from_caustic(a, b, p, lam, hint_dir);
                // Add BOTH directions of motion.  The caustic tangent has two
                // signs; a single sign fills only half the Liouville torus.
                starts.push((p, vel));
                starts.push((p, -vel));
            }
        }
        if i == start {
            return starts;
        }
    }
}

/// Extract `(λ_ellipse, λ_hyperbola)` of a confocal-quadrilateral domain.
/// Returns `None` for non-quadrilaterals (e.g. L-shapes with 6 arcs)
/// since the simple ellipse/hyperbola test is only valid for quadrilaterals.
fn confocal_bounds(domain: &domain::Domain) -> Option<(f32, f32)> {
    // Count quadric segments — must be exactly 4 for the simple test to be valid
    let mut ell = None;
    let mut hyp = None;
    let mut quad_count = 0;
    for seg in &domain.segments {
        if let domain::Segment::Quad { curve, .. } = seg {
            quad_count += 1;
            if curve.lambda < curve.b_param {
                ell = Some(curve.lambda);
            } else if curve.lambda < curve.a_param {
                hyp = Some(curve.lambda);
            }
        }
    }
    match (ell, hyp) {
        (Some(e), Some(h)) if quad_count == 4 => Some((e, h)),
        _ => None,
    }
}

/// Torus-mapping bounds `(lam_wall, beta)` for a domain.
///
/// * `lam_wall` — `λ₁` of the outer boundary (the ellipse wall).
/// * `beta` — `λ₂` of the hyperbola walls; `None` means the full ellipse.
///
/// For a confocal quadrilateral this is `(λ_ell, λ_hyp)`; for a full ellipse
/// (no hyperbola walls) it is `(0, None)`.  Non-quadrilaterals (L-shapes) fall
/// back to the full-ellipse convention.
pub fn torus_bounds(domain: &domain::Domain) -> (f32, Option<f32>) {
    match confocal_bounds(domain) {
        Some((ell, hyp)) => (ell, Some(hyp)),
        None => (0.0, None),
    }
}

/// Analytic inside-test for a confocal quadrilateral. A point is inside
/// when it lies within the boundary ellipse and between the two hyperbola
/// sheets. `slack` shrinks the region slightly inward, so borderline
/// slivers that merely touch the boundary are excluded.
fn confocal_inside(a: f32, b: f32, p: Vec2, lambda_ell: f32, lambda_hyp: f32, slack: f32) -> bool {
    // Ellipse: Q_λ_ell(p) < 0 means inside.  slack shrinks inward.
    let e_val = (b - lambda_ell) * p.x * p.x + (a - lambda_ell) * p.y * p.y
        - (a - lambda_ell) * (b - lambda_ell);
    let inside_ell = e_val < -slack;
    if !inside_ell {
        return false;
    }

    // Hyperbola: Q_λ_hyp(p) > 0 means between the two branches (inside).
    // For λ > B, the hyperbola opens left/right.  Q > 0 is the region
    // between the left and right sheets, which is the interior of the domain.
    let h_val = (b - lambda_hyp) * p.x * p.x + (a - lambda_hyp) * p.y * p.y
        - (a - lambda_hyp) * (b - lambda_hyp);
    let between_hyp = h_val > slack;
    between_hyp
}

/// Newton-project a point onto the quadric Q(x,y) = 0 along its gradient.
fn project_to_curve(quad: &quadratic::ConfocalQuadric, mut p: Vec2) -> Vec2 {
    for _ in 0..8 {
        let q = quad.eval(p);
        if q.abs() < 1e-9 {
            break;
        }
        let g = quad.grad(p);
        let g2 = g.length_squared();
        if g2 < 1e-12 {
            break;
        }
        p -= (q / g2) * g;
    }
    p
}
