use billiards::{domain, presets, torus::ConfocalParams};
use macroquad::prelude::*;

/// Inside test that handles both quadrilaterals and complex shapes.
/// Returns true if point is strictly inside the domain.
fn point_inside(p: Vec2, domain: &domain::Domain) -> bool {
    // Count quadric segments — only use analytic test for exactly 4 arcs.
    let mut ell_lambdas = Vec::new();
    let mut hyp_lambdas = Vec::new();
    let mut quad_count = 0;
    let mut aa = 0.0f32;
    let mut bb = 0.0f32;
    for seg in &domain.segments {
        if let domain::Segment::Quad { curve, .. } = seg {
            quad_count += 1;
            aa = curve.a_param;
            bb = curve.b_param;
            if curve.lambda < curve.b_param {
                ell_lambdas.push(curve.lambda);
            }
            if curve.lambda > curve.b_param {
                hyp_lambdas.push(curve.lambda);
            }
        }
    }

    if quad_count != 4 || ell_lambdas.is_empty() || hyp_lambdas.is_empty() {
        // Fall back to ray casting for non-quadrilaterals (e.g. L-shape).
        return domain.contains(p);
    }

    // Quadrilateral: the outer ellipse is the minimum λ, the inner hyperbola
    // is the minimum λ > B (the one that defines the right/left walls).
    // Quadrilateral: the outer ellipse is the minimum λ, the inner hyperbola
    // is the minimum λ > B (the one that defines the right/left walls).
    let lambda_ell = ell_lambdas.into_iter().reduce(f32::min).unwrap();
    let lambda_hyp = hyp_lambdas.into_iter().reduce(f32::min).unwrap();

    let inside_ell = (bb - lambda_ell) * p.x * p.x + (aa - lambda_ell) * p.y * p.y
        - (aa - lambda_ell) * (bb - lambda_ell)
        < 0.0;
    if !inside_ell {
        return false;
    }

    // Inside the hyperbola strip: Q_λ_hyp(p) > 0 (between the sheets).
    let h_val = (bb - lambda_hyp) * p.x * p.x + (aa - lambda_hyp) * p.y * p.y
        - (aa - lambda_hyp) * (bb - lambda_hyp);
    let between_hyp = h_val > 0.0;
    between_hyp
}

/// Generate a trajectory and check every bounce point is inside the domain.
fn check_trajectory_stays_inside(
    domain: &domain::Domain,
    p0: Vec2,
    v0: Vec2,
    max_steps: usize,
    label: &str,
    lam: f32,
    traj_idx: usize,
) {
    let segs = domain.trace(p0, v0, max_steps);
    if segs.is_empty() {
        // Zero-segment trajectories occur at degenerate/borderline start
        // points near the re-entrant corner.  Skip silently.
        return;
    }

    // Extract confocal bounds (if applicable) for analytic inside test.
    // Use the outer ellipse (min λ_ell) and the most restrictive hyperbola
    // (min λ_hyp > B, defining the right/left walls).
    let mut ell_lambdas = Vec::new();
    let mut hyp_lambdas = Vec::new();
    let mut aa = 0.0f32;
    let mut bb = 0.0f32;
    for seg in &domain.segments {
        if let domain::Segment::Quad { curve, .. } = seg {
            aa = curve.a_param;
            bb = curve.b_param;
            if curve.lambda < curve.b_param {
                ell_lambdas.push(curve.lambda);
            }
            if curve.lambda > curve.b_param {
                hyp_lambdas.push(curve.lambda);
            }
        }
    }

    // Use analytic test for quadrilaterals (4 quadric arcs only).
    // For complex shapes (L-shape with 6 arcs), the simple analytic test
    // using the outer hyperbola is WRONG because the left and right walls
    // are different. Fall back to ray casting.
    let quad_count = domain
        .segments
        .iter()
        .filter(|s| matches!(s, domain::Segment::Quad { .. }))
        .count();
    let use_analytic = !ell_lambdas.is_empty() && !hyp_lambdas.is_empty() && quad_count == 4;
    let lambda_ell = ell_lambdas.into_iter().reduce(f32::min);
    let lambda_hyp = hyp_lambdas.into_iter().reduce(f32::min);

    for (seg_i, &(p_start, hit)) in segs.iter().enumerate() {
        let mid = (p_start + hit) * 0.5;

        let inside = if use_analytic {
            let le = lambda_ell.unwrap();
            let lh = lambda_hyp.unwrap();
            let slack = 1e-2;
            let inside_ell = (bb - le) * mid.x * mid.x + (aa - le) * mid.y * mid.y
                - (aa - le) * (bb - le)
                < slack;
            if !inside_ell {
                false
            } else {
                let h_val =
                    (bb - lh) * mid.x * mid.x + (aa - lh) * mid.y * mid.y - (aa - lh) * (bb - lh);
                h_val > -slack
            }
        } else {
            domain.contains(mid)
        };

        assert!(
            inside,
            "{} | Λ={} | traj {} | seg {}: midpoint ({:.4}, {:.4}) outside domain! segment from ({:.4}, {:.4}) to ({:.4}, {:.4})",
            label, lam, traj_idx, seg_i, mid.x, mid.y, p_start.x, p_start.y, hit.x, hit.y,
        );
    }
}

/// Test trajectories for confocal billiards across the full Λ range.
fn test_confocal_trajectories(preset: &presets::Preset) {
    let domain = &preset.domain;
    let label = preset.label;

    // Extract boundary lambdas, using the same logic as the library.
    let cf = ConfocalParams::standard();
    let mut ell_lambdas = Vec::new();
    let mut hyp_lambdas = Vec::new();
    for seg in &domain.segments {
        if let domain::Segment::Quad { curve, .. } = seg {
            if curve.lambda < cf.b {
                ell_lambdas.push(curve.lambda);
            }
            if curve.lambda > cf.b {
                hyp_lambdas.push(curve.lambda);
            }
        }
    }
    assert!(
        !ell_lambdas.is_empty(),
        "{}: no ellipse boundary found",
        label
    );

    // For quadrilaterals, the outer ellipse has the minimum λ.
    // For L-shapes, take the most restrictive (minimum) ellipse λ.
    let lambda_ell = ell_lambdas.into_iter().reduce(f32::min).unwrap();
    let lambda_hyp = hyp_lambdas.into_iter().reduce(f32::min).unwrap();

    let e = 0.05;
    let ell_min = lambda_ell + e;
    let ell_max = cf.b - e;
    let hyp_min = lambda_hyp + e;
    let hyp_max = cf.a - e;

    let eps = 0.15;

    // Ellipse side
    let mut lam = ell_min;
    while lam < ell_max {
        let starts = billiards::get_start_points(lam, domain, true, vec2(0.0, 0.0));
        for (i, &(p, v)) in starts.iter().enumerate() {
            check_trajectory_stays_inside(domain, p, v, 300, label, lam, i);
        }
        lam += eps;
    }
    // Edge cases
    for &lam in &[ell_min, ell_max, (ell_min + ell_max) / 2.0] {
        let starts = billiards::get_start_points(lam, domain, true, vec2(0.0, 0.0));
        for (i, &(p, v)) in starts.iter().enumerate() {
            check_trajectory_stays_inside(domain, p, v, 300, label, lam, i);
        }
    }

    // Hyperbola side
    let mut lam = hyp_min;
    while lam < hyp_max {
        let starts = billiards::get_start_points(lam, domain, true, vec2(0.0, 0.0));
        for (i, &(p, v)) in starts.iter().enumerate() {
            check_trajectory_stays_inside(domain, p, v, 300, label, lam, i);
        }
        lam += eps;
    }
    for &lam in &[hyp_min, hyp_max, (hyp_min + hyp_max) / 2.0] {
        let starts = billiards::get_start_points(lam, domain, true, vec2(0.0, 0.0));
        for (i, &(p, v)) in starts.iter().enumerate() {
            check_trajectory_stays_inside(domain, p, v, 300, label, lam, i);
        }
    }
}

/// Test trajectories for polyline billiards across the full angle range.
fn test_polyline_trajectories(preset: &presets::Preset) {
    let domain = &preset.domain;
    let label = preset.label;
    let center = preset.start_center;

    // Verify starting center is inside
    assert!(
        domain.contains(center),
        "{}: start center ({}, {}) is outside domain!",
        label,
        center.x,
        center.y,
    );

    let eps = 0.15;
    let mut theta = -1.0;
    while theta <= 1.0 {
        let starts = billiards::get_start_points(theta, domain, false, center);
        for (i, &(p, v)) in starts.iter().enumerate() {
            check_trajectory_stays_inside(domain, p, v, 300, label, theta, i);
        }
        theta += eps;
    }
    for &theta in &[-1.0, -0.5, 0.0, 0.5, 1.0] {
        let starts = billiards::get_start_points(theta, domain, false, center);
        for (i, &(p, v)) in starts.iter().enumerate() {
            check_trajectory_stays_inside(domain, p, v, 300, label, theta, i);
        }
    }
}

// ---------------------------------------------------------------
// Individual tests — each preset gets its own test
// ---------------------------------------------------------------

#[test]
fn test_square_confocal_trajectories() {
    let preset = presets::all_presets()
        .into_iter()
        .find(|p| p.label.starts_with("Square: ellipse"))
        .expect("Square confocal preset not found");
    test_confocal_trajectories(&preset);
}

#[test]
fn test_thin_confocal_trajectories() {
    let preset = presets::all_presets()
        .into_iter()
        .find(|p| p.label.starts_with("Thin"))
        .expect("Thin confocal preset not found");
    test_confocal_trajectories(&preset);
}

#[test]
fn test_flat_confocal_trajectories() {
    let preset = presets::all_presets()
        .into_iter()
        .find(|p| p.label.starts_with("Flat"))
        .expect("Flat confocal preset not found");
    test_confocal_trajectories(&preset);
}

#[test]
fn test_lshape_confocal_trajectories() {
    let preset = presets::all_presets()
        .into_iter()
        .find(|p| p.label.contains("3π/2"))
        .expect("L-shape confocal preset not found");
    test_confocal_trajectories(&preset);
}

#[test]
fn test_square_polyline_trajectories() {
    let preset = presets::all_presets()
        .into_iter()
        .find(|p| p.label == "Square (polyline)")
        .expect("Square polyline preset not found");
    test_polyline_trajectories(&preset);
}

#[test]
fn test_lshape_polyline_trajectories() {
    let preset = presets::all_presets()
        .into_iter()
        .find(|p| p.label == "L-shape (polyline)")
        .expect("L-shape polyline preset not found");
    test_polyline_trajectories(&preset);
}

/// Verify that every start point produced by `start_points_on_caustic` is
/// inside the domain and lies on the correct caustic curve.
#[test]
fn test_all_lambda_start_points() {
    let configs = presets::all_presets();
    let eps: f32 = 0.1;

    for preset in &configs {
        let domain = &preset.domain;
        let has_quad = domain
            .segments
            .iter()
            .any(|seg| matches!(seg, domain::Segment::Quad { .. }));
        if !has_quad {
            continue;
        }

        let cf = ConfocalParams::standard();
        let mut ell_lambdas = Vec::new();
        let mut hyp_lambdas = Vec::new();
        for seg in &domain.segments {
            if let domain::Segment::Quad { curve, .. } = seg {
                if curve.lambda < cf.b {
                    ell_lambdas.push(curve.lambda);
                }
                if curve.lambda > cf.b {
                    hyp_lambdas.push(curve.lambda);
                }
            }
        }
        if ell_lambdas.is_empty() || hyp_lambdas.is_empty() {
            continue;
        }

        let lambda_ell = ell_lambdas.into_iter().reduce(f32::min).unwrap();
        let lambda_hyp = hyp_lambdas.into_iter().reduce(f32::min).unwrap();

        let e = 0.05f32;
        let ell_min = lambda_ell + e;
        let ell_max = cf.b - e;
        let mut lam = ell_min;
        while lam < ell_max {
            let starts = billiards::start_points_on_caustic(lam, domain);
            for (i, &(p, _)) in starts.iter().enumerate() {
                assert!(
                    point_inside(p, domain),
                    "{} | Λ={} | start {}: outside! pos=({}, {})",
                    preset.label,
                    lam,
                    i,
                    p.x,
                    p.y,
                );
                let q = billiards::quadratic::confocal(cf, lam).eval(p);
                assert!(
                    q.abs() < 1e-3,
                    "{} | Λ={} | start {}: off caustic! Q={} pos=({}, {})",
                    preset.label,
                    lam,
                    i,
                    q,
                    p.x,
                    p.y,
                );
            }
            lam += eps;
        }

        let hyp_min = lambda_hyp + e;
        let hyp_max = cf.a - e;
        let mut lam = hyp_min;
        while lam < hyp_max {
            let starts = billiards::start_points_on_caustic(lam, domain);
            for (i, &(p, _)) in starts.iter().enumerate() {
                assert!(
                    point_inside(p, domain),
                    "{} | Λ={} | start {}: outside! pos=({}, {})",
                    preset.label,
                    lam,
                    i,
                    p.x,
                    p.y,
                );
                let q = billiards::quadratic::confocal(cf, lam).eval(p);
                assert!(
                    q.abs() < 1e-3,
                    "{} | Λ={} | start {}: off caustic! Q={} pos=({}, {})",
                    preset.label,
                    lam,
                    i,
                    q,
                    p.x,
                    p.y,
                );
            }
            lam += eps;
        }
    }
}
