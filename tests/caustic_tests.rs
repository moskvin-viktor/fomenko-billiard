use billiards::domain;
use billiards::quadratic::ConfocalQuadric;
use macroquad::prelude::*;

const A: f32 = 4.0;
const B: f32 = 1.0;

/// Check if a point is inside a confocal quadrilateral analytically.
fn point_inside(p: Vec2, domain: &domain::Domain) -> bool {
    let (mut lambda_ell, mut lambda_hyp) = (f32::MAX, f32::MAX);
    let (mut aa, mut bb) = (0.0, 0.0);
    for seg in &domain.segments {
        if let domain::Segment::Quad { curve, .. } = seg {
            aa = curve.a_param;
            bb = curve.b_param;
            if curve.lambda < curve.b_param {
                lambda_ell = curve.lambda;
            }
            if curve.lambda > curve.b_param {
                lambda_hyp = curve.lambda;
            }
        }
    }
    if lambda_ell == f32::MAX || lambda_hyp == f32::MAX {
        return domain.contains(p);
    }
    let inside_ell = (bb - lambda_ell) * p.x * p.x + (aa - lambda_ell) * p.y * p.y
        - (aa - lambda_ell) * (bb - lambda_ell)
        < 0.0;
    let between_hyp = p.x.abs() < (aa - lambda_hyp).sqrt();
    inside_ell && between_hyp
}

fn test_start_point_inside(a: f32, b: f32, lam: f32, domain: &domain::Domain) {
    let starts = crate::start_points_on_caustic(a, b, lam, domain);
    for (i, &(p, _)) in starts.iter().enumerate() {
        assert!(
            point_inside(p, domain),
            "Start point {} for Λ={} is outside! pos=({}, {})",
            i,
            lam,
            p.x,
            p.y
        );
    }
}

#[test]
fn test_all_lambda_values() {
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

        let (mut lambda_ell, mut lambda_hyp) = (f32::MAX, f32::MAX);
        for seg in &domain.segments {
            if let domain::Segment::Quad { curve, .. } = seg {
                if curve.lambda < B {
                    lambda_ell = curve.lambda;
                }
                if curve.lambda > B {
                    lambda_hyp = curve.lambda;
                }
            }
        }

        let e = 0.05f32;
        let ell_min = lambda_ell + e;
        let ell_max = B - e;
        let mut lam = ell_min;
        while lam < ell_max {
            test_start_point_inside(A, B, lam, domain);
            lam += eps;
        }
        test_start_point_inside(A, B, ell_max, domain);

        let hyp_min = lambda_hyp + e;
        let hyp_max = A - e;
        let mut lam = hyp_min;
        while lam < hyp_max {
            test_start_point_inside(A, B, lam, domain);
            lam += eps;
        }
        test_start_point_inside(A, B, hyp_max, domain);
    }

    eprintln!("All tests passed!");
}
