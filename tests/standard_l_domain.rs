//! Focused tests for `Segment::intersect` on the standard L's quadric arcs.

use billiards::{domain, presets, torus::ConfocalParams};
use macroquad::prelude::*;

fn std_l() -> billiards::domain::Domain {
    presets::all_presets()
        .into_iter()
        .find(|p| p.label.contains("3π/2"))
        .expect("standard L preset")
        .domain
}

/// Map a (λ₁, λ₂) point in the first quadrant to physical (x, y).
fn lam_to_xy(cf: &ConfocalParams, lam1: f32, lam2: f32) -> Vec2 {
    let a = cf.a;
    let b = cf.b;
    let x = (((a - lam1) * (a - lam2)) / (a - b)).sqrt();
    let y = (((b - lam1) * (b - lam2)) / (b - a)).sqrt();
    vec2(x, y)
}

/// Each of the 6 arcs must be intersected by a ray aimed at its midpoint, and
/// the hit must lie on the correct branch (x>0 for the in-quadrant table).
#[test]
fn each_arc_intersects_at_midpoint() {
    let dom = std_l();
    assert_eq!(dom.segments.len(), 6, "standard L has 6 arcs");

    for (i, seg) in dom.segments.iter().enumerate() {
        let domain::Segment::Quad { curve, a, b } = seg else {
            panic!("seg {} should be a quad", i);
        };
        // Midpoint of the arc in λ-space, mapped to physical space.
        // For a hyperbola arc the endpoints share a λ₂; for an ellipse they
        // share a λ₁.  Use the geometric midpoint of a and b.
        let mid = (*a + *b) * 0.5;
        // Cast a ray from the origin toward the midpoint.
        let dir = mid.normalize();
        let hit = dom.intersect(vec2(0.0, 0.0), dir);
        assert!(
            hit.is_some(),
            "seg {} (λ={}) should be hit by ray toward midpoint {:?}",
            i,
            curve.lambda,
            mid
        );
        let (_t, _idx, hit_pt) = hit.unwrap();
        assert!(
            hit_pt.x > 0.0,
            "seg {} hit on wrong branch: {:?}",
            i,
            hit_pt
        );
    }
}

/// The base-bottom hyperbola arc (seg 1, β₁) must NOT accept hits on the far
/// (x<0) branch.  A ray from the origin pointing at the far branch should not
/// report a hit on this arc.
/// Rays from the origin into the first quadrant must intersect only the
/// correct (x>0) branch — a hit on a hyperbola's far sheet would be a sign
/// that branch-awareness is broken.
#[test]
fn hyperbola_rejects_far_branch() {
    let dom = std_l();
    // Sweep rays from the origin around the first quadrant.
    for i in 0..180 {
        let ang = std::f32::consts::FRAC_PI_2 * i as f32 / 180.0; // 0..90°
        let dir = vec2(ang.cos(), ang.sin());
        if let Some((_t, idx, hit)) = dom.intersect(vec2(0.0, 0.0), dir) {
            if let domain::Segment::Quad { curve, a, .. } = &dom.segments[idx] {
                assert!(
                    curve.branch_sign(hit) == curve.branch_sign(*a),
                    "ray θ={:.1}° hit wrong branch on seg {}",
                    ang * 180.0 / std::f32::consts::PI,
                    idx
                );
                assert!(
                    hit.x > 0.0,
                    "ray {}° hit x<0 (far branch): {:?}",
                    (ang * 180.0 / std::f32::consts::PI),
                    hit
                );
            }
        }
    }
}

/// Trajectories launched inside the standard L must bounce many times and
/// every segment midpoint must stay inside the domain (bounce points land on
/// the boundary, where point-in-domain is only tolerance-accurate).
#[test]
fn standard_l_traces() {
    let cf = ConfocalParams::standard();
    let dom = std_l();
    let p = lam_to_xy(&cf, 0.6, 1.7);
    let v = vec2(0.0, 1.0);
    let segs = dom.trace(p, v, 100);
    assert!(
        segs.len() >= 20,
        "trajectory should bounce many times, got {}",
        segs.len()
    );
    for &(a, b) in &segs {
        let mid = (a + b) * 0.5;
        assert!(dom.contains(mid), "segment midpoint outside: {:?}", mid);
    }
}

/// Occupancy: base and tall leg inside, notch outside.
#[test]
fn standard_l_occupancy() {
    let cf = ConfocalParams::standard();
    let dom = std_l();
    assert!(dom.contains(lam_to_xy(&cf, 0.6, 1.7)), "base point inside");
    assert!(
        dom.contains(lam_to_xy(&cf, 0.2, 2.3)),
        "tall leg point inside"
    );
    assert!(
        !dom.contains(lam_to_xy(&cf, 0.6, 2.3)),
        "notch point outside"
    );
    assert!(!dom.contains(lam_to_xy(&cf, -0.5, 2.0)), "outside ellipse");
}
