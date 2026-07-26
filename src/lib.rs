pub mod domain;
pub mod phase3d;
pub mod presets;
pub mod quadratic;

pub const A: f32 = 4.0;
pub const B: f32 = 1.0;

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

/// Extract (λ_ellipse, λ_hyperbola) of a confocal-quadrilateral domain.
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
