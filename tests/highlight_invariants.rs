//! Invariants the app's per-torus highlight behaviour must satisfy.
//!
//! These pin down the "one red trajectory per torus" logic that lives in the
//! app (`main.rs`): for any second-integral value, the trajectories drawn on
//! the billiard and the short phase-space highlights drawn on the torus must
//! pick out exactly the same set of *distinct* tori — and match the topology
//! (1 torus for a hyperbola caustic, 2 for an ellipse caustic in a confocal
//! quadrilateral).
//!
//! They are the behaviour the (upcoming) `rebuild()` refactor must preserve.

use billiards::{phase3d, presets, torus::ConfocalParams};
use macroquad::prelude::{vec2, Vec2};

/// Choose a valid Λ on the ellipse (Λ < B) side for a domain.
fn ellipse_lambda(domain: &billiards::domain::Domain) -> f32 {
    let cf = ConfocalParams::standard();
    let mut ell = f32::MAX;
    for seg in &domain.segments {
        if let billiards::domain::Segment::Quad { curve, .. } = seg {
            if curve.lambda < cf.b {
                ell = ell.min(curve.lambda);
            }
        }
    }
    (ell + 0.05 + cf.b - 0.05) / 2.0
}

/// Choose a valid Λ on the hyperbola (Λ > B) side for a domain.
fn hyperbola_lambda(domain: &billiards::domain::Domain) -> f32 {
    let cf = ConfocalParams::standard();
    let mut hyp = f32::MAX;
    for seg in &domain.segments {
        if let billiards::domain::Segment::Quad { curve, .. } = seg {
            if curve.lambda > cf.b {
                hyp = hyp.min(curve.lambda);
            }
        }
    }
    hyp + 0.2
}

/// Trace a 4-bounce highlight from a start, exactly as the app's
/// `build_highlights` does (skipping empty traces like the app's filter).
fn highlight(domain: &billiards::domain::Domain, p: Vec2, v: Vec2) -> Vec<phase3d::PhasePoint> {
    let bounds = billiards::torus_bounds(domain);
    phase3d::sample_trajectory_phase_dense(domain, p, v, 4, 8, bounds)
}

/// The app's `one_per_torus` torus key for a single start point (same rule as
/// `to_torus`): hyperbola/polyline → 0; quadrilateral + ellipse → sign(y);
/// full-ellipse fallback (L-shape) → sign of angular momentum.
fn torus_key(p: Vec2, v: Vec2, is_confocal: bool, lam: f32, is_quad: bool) -> u32 {
    let cf = ConfocalParams::standard();
    if !is_confocal || lam >= cf.b {
        0
    } else if is_quad {
        if p.y >= 0.0 {
            0
        } else {
            1
        }
    } else if p.x * v.y - p.y * v.x > 0.0 {
        0
    } else {
        1
    }
}

/// True if this preset is a confocal quadrilateral (its `torus_bounds` has a
/// hyperbola wall), as opposed to an L-shape (full-ellipse fallback).
fn is_quad_confocal(domain: &billiards::domain::Domain) -> bool {
    billiards::torus_bounds(domain).1.is_some()
}

/// True if the domain is a confocal quadrilateral boundary (exactly 4 arcs).
fn ellipse_boundary_is_quadrilateral(domain: &billiards::domain::Domain) -> bool {
    let quad_arcs = domain
        .segments
        .iter()
        .filter(|s| matches!(s, billiards::domain::Segment::Quad { .. }))
        .count();
    quad_arcs == 4
}

// ---------------------------------------------------------------------------
// Topology
// ---------------------------------------------------------------------------

/// Distinct `torus_index` values reached by 4-bounce highlights traced from
/// every start point, mirroring `build_highlights` (empty traces skipped).
fn distinct_tori_for_lambda(
    domain: &billiards::domain::Domain,
    lam: f32,
    is_confocal: bool,
    center: Vec2,
) -> Vec<u32> {
    let starts = billiards::get_start_points(lam, domain, is_confocal, center);
    assert!(!starts.is_empty(), "no start points for Λ={}", lam);

    let mut seen: Vec<u32> = Vec::new();
    for (p, v) in starts {
        let traj = highlight(domain, p, v);
        if let Some(pt) = traj.first() {
            if !seen.contains(&pt.torus_index) {
                seen.push(pt.torus_index);
            }
        }
    }
    seen
}

/// Hyperbola caustics (Λ > B) have a single connected torus, for every
/// confocal preset (quadrilateral or L-shape alike).
#[test]
fn hyperbola_caustic_is_single_torus() {
    for preset in presets::all_presets() {
        if !preset.is_confocal {
            continue;
        }
        let ids = distinct_tori_for_lambda(
            &preset.domain,
            hyperbola_lambda(&preset.domain),
            true,
            preset.start_center,
        );
        assert!(
            ids.len() == 1,
            "{} hyperbola: expected exactly 1 torus in highlights, got {:?}",
            preset.label,
            ids
        );
    }
}

/// An ellipse caustic in a confocal QUADRILATERAL splits the table into two
/// tori (upper / lower region).
#[test]
fn quadrilateral_ellipse_is_two_tori() {
    for preset in presets::all_presets() {
        if !preset.is_confocal {
            continue;
        }
        if !is_quad_confocal(&preset.domain) {
            // L-shape is not a quadrilateral; handled by the consistency test.
            continue;
        }
        let ids = distinct_tori_for_lambda(
            &preset.domain,
            ellipse_lambda(&preset.domain),
            true,
            preset.start_center,
        );
        assert!(
            ids.len() == 2,
            "{} ellipse: expected exactly 2 tori in highlights, got {:?}",
            preset.label,
            ids
        );
    }
}

/// Polyline domains have a single start trajectory → a single torus.
#[test]
fn polyline_is_single_torus() {
    for preset in presets::all_presets() {
        if preset.is_confocal {
            continue;
        }
        let ids = distinct_tori_for_lambda(&preset.domain, 0.2, false, preset.start_center);
        assert!(
            ids.len() == 1,
            "{}: expected exactly 1 torus in highlights, got {:?}",
            preset.label,
            ids
        );
    }
}

// ---------------------------------------------------------------------------
// 2D / 3D consistency (the refactor's core invariant)
// ---------------------------------------------------------------------------

/// The billiard draws one red trajectory per torus (picked by the app's
/// `one_per_torus`), while the 3D view draws short highlights. Both must
/// report the SAME number of tori and stay in sync: every torus the billiard
/// picks must have non-empty highlight points on the torus surface.
///
/// This is the invariant the user hit ("red on billiard but nothing on the
/// torus", "wrong torus"). It fails when a picked start point traces an empty
/// highlight (degenerate boundary start) so the 2D and 3D views diverge.
#[test]
fn billiard_picks_are_consistent_with_torus_highlights() {
    for preset in presets::all_presets() {
        let domain = &preset.domain;
        let center = preset.start_center;
        let lams: Vec<f32> = if preset.is_confocal {
            vec![ellipse_lambda(domain), hyperbola_lambda(domain)]
        } else {
            vec![0.2]
        };

        for lam in lams {
            let starts = billiards::get_start_points(lam, domain, preset.is_confocal, center);
            let is_quad = is_quad_confocal(domain);

            // Distinct keys the billiard would pick (one_per_torus, in order),
            // each mapped to the start that introduced it.
            let mut picked: Vec<(u32, usize)> = Vec::new(); // (key, start idx)
                                                            // Distinct tori that actually have points on the torus.
            let mut shown_tori: Vec<u32> = Vec::new();
            let mut empty_picks: Vec<(usize, Vec2)> = Vec::new();

            for (i, &(p, v)) in starts.iter().enumerate() {
                let key = torus_key(p, v, preset.is_confocal, lam, is_quad);
                if picked.iter().all(|&(k, _)| k != key) {
                    picked.push((key, i));
                }

                let traj = highlight(domain, p, v);
                if let Some(pt) = traj.first() {
                    if !shown_tori.contains(&pt.torus_index) {
                        shown_tori.push(pt.torus_index);
                    }
                } else {
                    empty_picks.push((i, p));
                }
            }

            assert!(
                picked.len() == shown_tori.len(),
                "{} Λ={}: 2D picks {} keys {:?} but 3D shows {} tori {:?}; \
                 empty-phase starts (idx → pos): {:?}",
                preset.label,
                lam,
                picked.len(),
                picked.iter().map(|&(k, _)| k).collect::<Vec<_>>(),
                shown_tori.len(),
                shown_tori,
                empty_picks
            );
        }
    }
}

/// The highlight torus indices must be a subset of the torus indices of the
/// dense phase-space fill (both must agree on which disconnected regions the
/// caustic visits).
#[test]
fn highlights_are_on_same_tori_as_dense_fill() {
    for preset in presets::all_presets() {
        if !preset.is_confocal {
            continue;
        }
        let domain = &preset.domain;
        let bounds = billiards::torus_bounds(domain);

        for lam in [ellipse_lambda(domain), hyperbola_lambda(domain)] {
            let dense = billiards::dense_caustic_starts(domain, lam, 24);
            assert!(
                !dense.is_empty(),
                "{} Λ={}: no dense starts",
                preset.label,
                lam
            );
            let dense_ids: std::collections::HashSet<u32> = dense
                .iter()
                .map(|&(p, v)| {
                    phase3d::sample_trajectory_phase_dense(domain, p, v, 4, 8, bounds)
                        .first()
                        .map(|q| q.torus_index)
                        .unwrap_or(u32::MAX)
                })
                .collect();

            let starts = billiards::get_start_points(lam, domain, true, preset.start_center);
            for (p, v) in starts {
                if let Some(pt) = highlight(domain, p, v).first() {
                    assert!(
                        dense_ids.contains(&pt.torus_index),
                        "{} Λ={}: highlight torus {} not among dense-fill tori {:?}",
                        preset.label,
                        lam,
                        pt.torus_index,
                        dense_ids
                    );
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Geometry invariants the refactor relies on
// ---------------------------------------------------------------------------

/// Bounce geometry: `reflect` must preserve the incoming vector's length
/// (unit speed, conserved H) for both line and quadric segments.
#[test]
fn reflect_is_an_isometry() {
    use billiards::{domain, quadratic};
    let cf = ConfocalParams::standard();
    let b = cf.b;
    let quad = quadratic::confocal(cf, 0.5);

    // A point on each wall where reflection happens.
    let line = domain::Segment::Line {
        a: vec2(-1.0, -1.0),
        b: vec2(1.0, -1.0),
    };
    let arc = domain::Segment::Quad {
        curve: quad,
        a: vec2(0.0, -1.0),
        b: vec2(0.0, 1.0),
    };
    let wall_points = [vec2(0.0, -1.0), vec2(0.0, b.sqrt())];
    let segs: [&domain::Segment; 2] = [&line, &arc];

    let n_ang = 360;
    for i in 0..n_ang {
        let th = std::f32::consts::PI * 2.0 * i as f32 / n_ang as f32;
        let dir = vec2(th.cos(), th.sin());
        let p = wall_points[i % 2];
        let r = segs[i % 2].reflect(p, dir);
        let len = r.length();
        assert!(
            (len - 1.0).abs() < 1e-4,
            "reflect on seg {}: length {} != 1 for dir (θ={})",
            i % 2,
            len,
            th
        );
    }
}

/// `torus_bounds` must report the correct (λ_wall, β) convention per domain
/// type — the hyperbola-wall quadrilateral vs the full-ellipse fallback for
/// L-shape and polyline domains.
#[test]
fn torus_bounds_conventions() {
    for preset in presets::all_presets() {
        let domain = &preset.domain;
        let (ell, beta) = billiards::torus_bounds(domain);

        // λ_wall is the outer ellipse boundary λ (the minimum λ_c in the quad
        // arcs), or 0.0 for the full-ellipse fallback of a non-quadrilateral.
        let mut ell_boundary: Option<f32> = None;
        for seg in &domain.segments {
            if let billiards::domain::Segment::Quad { curve, .. } = seg {
                if curve.lambda < ConfocalParams::standard().b {
                    ell_boundary = Some(ell_boundary.map_or(curve.lambda, |e| e.min(curve.lambda)));
                }
            }
        }
        let expect_hyp_wall = ellipse_boundary_is_quadrilateral(domain);
        if expect_hyp_wall {
            let boundary = ell_boundary.expect("quadrilateral has ellipse wall");
            assert!(
                (ell - boundary).abs() < 1e-4,
                "{}: expected λ_wall {} (ellipse boundary), got {}",
                preset.label,
                boundary,
                ell
            );
            assert!(
                beta.is_some(),
                "{}: quadrilateral should have a hyperbola wall (β), got None",
                preset.label
            );
        } else {
            assert!(
                (ell - 0.0).abs() < 1e-4 && beta.is_none(),
                "{}: expected full-ellipse fallback (0.0, None), got ({}, {:?})",
                preset.label,
                ell,
                beta
            );
        }
    }
}
