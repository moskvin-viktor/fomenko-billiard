//! Critical-layer behavior: trajectories must exist on critical (degenerate)
//! caustic layers, phase manifolds must exist on special layers, and the UI
//! must step through every layer (critical *and* regular) without skipping.
//!
//! A *degenerate* critical layer is one where the caustic collapses to a border
//! piece and the ball slides along it: the ellipse wall `λ_ell`, the focal
//! separatrix `b`, a hyperbola wall `λ_hyp`, or the focal axis `a`.  These are
//! the layers the old code returned zero trajectories for.  Molecule-critical
//! values (genus_jump, edge_swap, birth/death) are still *regular* levels and
//! are handled by the ordinary caustic sampler.

use billiards::presets::Preset;
use billiards::torus::ConfocalParams;
use macroquad::prelude::Vec2;

fn cf() -> ConfocalParams {
    ConfocalParams::standard()
}

fn confocal_presets() -> Vec<Preset> {
    billiards::presets::all_presets()
        .into_iter()
        .filter(|p| p.is_confocal)
        .collect()
}

/// The degenerate critical values of a domain: `b`, `a`, and any wall
/// (`λ_ell`, `λ_hyp`).  These are where the caustic collapses to a border piece.
fn degenerate_layers(preset: &Preset) -> Vec<f32> {
    let cf = cf();
    let s = billiards::confocal::ConfocalStructure::of_domain(&preset.domain).unwrap();
    let mut v = vec![cf.b, cf.a, s.lambda_ell];
    if let Some(h) = s.lambda_hyp {
        v.push(h);
    }
    v.retain(|&x| x > 0.0 && x <= cf.a + 1e-6);
    v.sort_by(|a, b| a.total_cmp(b));
    v.dedup();
    v
}

/// A point is inside (or on) the confocal domain, allowing boundary points so
/// wall-sliding starts pass.  The analytic structure test with a small slack
/// is exact for quadrilaterals/full ellipses; we fall back to ray-cast otherwise.
fn point_inside(preset: &Preset, p: Vec2) -> bool {
    let s = billiards::confocal::ConfocalStructure::of_domain(&preset.domain).unwrap();
    if s.is_quadrilateral || s.is_full_ellipse {
        // A small *negative* slack widens the region so exact boundary/wall
        // points are accepted (the strict test excludes them).
        s.contains(p, -5e-3)
    } else {
        preset.domain.contains(p)
    }
}

/// The interior degenerate segments — the focal separatrix `b` (x-axis segment)
/// and the focal axis `a` (vertical segment) — lie strictly inside the domain,
/// so the ball truly slides back and forth along them and `trace` produces
/// real trajectories.
fn is_interior_segment(lam: f32) -> bool {
    let cf = cf();
    (lam - cf.b).abs() < 1e-4 || (lam - cf.a).abs() < 1e-4
}

/// Critical/degenerate start points via the dedicated sampler (wall points lie
/// on the border, so this bypasses the strict-inside `get_start_points`).
fn critical_starts(preset: &Preset, lam: f32) -> Vec<(Vec2, Vec2)> {
    let s = billiards::confocal::ConfocalStructure::of_domain(&preset.domain).unwrap();
    billiards::confocal::critical_caustic_starts(&s, &preset.domain, lam, 64).unwrap_or_default()
}

// ---------------------------------------------------------------------------
// 1. Trajectories exist on critical (degenerate) layers
// ---------------------------------------------------------------------------

/// Every degenerate layer of every confocal preset yields finite trajectory
/// start points.  For the interior degenerate segments (`b`, `a`) the ball
/// slides back and forth and `trace` produces real trajectories; for the wall
/// layers the border points feed the collapsed phase sheet (rendered directly),
/// and the exact-corner slide may genuinely produce zero bounces.
#[test]
fn trajectories_exist_on_degenerate_layers() {
    for preset in confocal_presets() {
        let domain = &preset.domain;
        let layers = degenerate_layers(&preset);
        for &lam in &layers {
            let starts = critical_starts(&preset, lam);
            // Some degenerate pieces don't intersect the domain at all (e.g.
            // the focal segment for the L-shape) — skip those; the ones that do
            // must be finite and slide.
            if starts.is_empty() {
                continue;
            }

            for &(p, v) in starts.iter().take(16) {
                assert!(
                    p.x.is_finite() && p.y.is_finite() && v.x.is_finite() && v.y.is_finite(),
                    "{}: Λ={} start non-finite",
                    preset.label,
                    lam
                );
            }

            // Interior segments must actually slide (non-empty trace).
            if is_interior_segment(lam) {
                let traj = domain.trace(starts[0].0, starts[0].1, 50);
                assert!(
                    !traj.is_empty(),
                    "{}: Λ={} interior segment produced empty sliding trajectory",
                    preset.label,
                    lam
                );
            }
        }
    }
}

/// The degenerate start point lies *on* or *inside* the domain (the ball slides
/// along the border).  Uses the analytic membership test, which is exact for
/// wall points.
#[test]
fn degenerate_starts_are_inside_domain() {
    for preset in confocal_presets() {
        for &lam in &degenerate_layers(&preset) {
            let starts = critical_starts(&preset, lam);
            if starts.is_empty() {
                continue;
            }
            for &(p, _) in starts.iter().take(8) {
                assert!(
                    point_inside(&preset, p),
                    "{}: Λ={} degenerate start at ({:.3},{:.3}) outside domain",
                    preset.label,
                    lam,
                    p.x,
                    p.y
                );
            }
        }
    }
}

// ---------------------------------------------------------------------------
// 2. Phase manifolds exist on special layers
// ---------------------------------------------------------------------------

/// Every special layer is reachable and, for the degenerate ones, the critical
/// sheet (collapsed 2D phase manifold) exists and is non-empty.  Molecule-
/// critical regular levels are handled by the ordinary torus path.
#[test]
fn phase_manifold_exists_on_special_layers() {
    let cf = cf();
    for preset in confocal_presets() {
        let domain = &preset.domain;
        let layers = billiards::molecule::special_layers(domain, &cf);
        for &lam in &layers {
            // get_start_points must not panic.
            let _ = billiards::get_start_points(lam, domain, true, Vec2::ZERO);

            // Degenerate layers collapse to a 2D phase sheet.
            if degenerate_layers(&preset)
                .iter()
                .any(|&d| (d - lam).abs() < 1e-4)
            {
                let sheet = billiards::phase3d::critical_phase_sheet(domain, lam, 64);
                assert!(
                    sheet.is_some(),
                    "{}: Λ={} critical_phase_sheet returned None",
                    preset.label,
                    lam
                );
                // A degenerate piece may be empty when the border piece doesn't
                // intersect the domain (e.g. the focal segment for the L-shape),
                // but it must not be an error.
                let _ = sheet.unwrap();
            } else {
                // Regular layer: ordinary phase sampling stays finite.
                let bounds = billiards::confocal::ConfocalStructure::of_domain(domain)
                    .map(|s| s.torus_bounds())
                    .unwrap_or((0.0, None));
                let starts = billiards::dense_caustic_starts(domain, lam, 4);
                for (p, v) in starts.into_iter().take(2) {
                    let pts = billiards::phase3d::sample_trajectory_phase_dense(
                        domain, p, v, 20, 4, bounds,
                    );
                    for q in &pts {
                        assert!(
                            q.theta1.is_finite(),
                            "{}: Λ={} theta1 NaN",
                            preset.label,
                            lam
                        );
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// 3. UI steps through all layers (critical + regular)
// ---------------------------------------------------------------------------

/// The marker list that the UI steps through is strictly increasing, contains
/// every special layer, and includes regular layers in the reachable bands.
#[test]
fn markers_include_special_and_regular_layers() {
    let cf = cf();
    for preset in confocal_presets() {
        let markers = billiards::molecule_view::layer_markers(&preset);
        let specials = billiards::molecule::special_layers(&preset.domain, &cf);

        for &s in &specials {
            assert!(
                markers.iter().any(|m| (m.lam - s).abs() < 1e-4),
                "{}: special layer {} missing from markers",
                preset.label,
                s
            );
        }

        let n_regular = markers.len().saturating_sub(specials.len());
        assert!(
            n_regular >= 1,
            "{}: no regular layers between specials",
            preset.label
        );

        for w in markers.windows(2) {
            assert!(
                w[0].lam < w[1].lam,
                "{}: markers not strictly increasing: {} >= {}",
                preset.label,
                w[0].lam,
                w[1].lam
            );
        }
    }
}

/// The standard first domain (Square): special layers are exactly
/// `{b=1.0, hyperbola wall=2.5, a=4.0}`, each has trajectories, and the step
/// list visits the reachable elliptic + hyperbolic bands with regulars in
/// between (never the forbidden gap `(1.0, 2.5)`).
#[test]
fn square_all_three_degenerate_layers_have_trajectories_and_manifolds() {
    let cf = cf();
    let preset = confocal_presets()
        .into_iter()
        .find(|p| p.label.starts_with("Square: ellipse"))
        .expect("Square preset");
    let domain = &preset.domain;

    assert_eq!(
        billiards::molecule::special_layers(domain, &cf),
        vec![1.0, 2.5, 4.0],
        "Square special layers"
    );

    for &lam in &[1.0f32, 2.5, 4.0] {
        assert!(
            !critical_starts(&preset, lam).is_empty(),
            "Λ={lam} no critical starts"
        );
    }

    // The step list skips the forbidden gap and has regular layers on both
    // reachable sides of it.
    let layers = billiards::molecule::all_layers(domain, &cf);
    let ell_regular: Vec<f32> = layers
        .iter()
        .copied()
        .filter(|&l| l > 0.0 && l < 1.0)
        .collect();
    let hyp_regular: Vec<f32> = layers
        .iter()
        .copied()
        .filter(|&l| l > 2.5 && l < 4.0)
        .collect();
    assert!(
        !ell_regular.is_empty(),
        "Square: no elliptic regular layers"
    );
    assert!(
        ell_regular
            .iter()
            .all(|&l| !billiards::get_start_points(l, domain, true, Vec2::ZERO).is_empty()),
        "Square: an elliptic regular layer is unreachable"
    );
    assert!(
        !hyp_regular.is_empty(),
        "Square: no hyperbolic regular layers"
    );
    assert!(
        hyp_regular
            .iter()
            .all(|&l| !billiards::get_start_points(l, domain, true, Vec2::ZERO).is_empty()),
        "Square: a hyperbolic regular layer is unreachable"
    );
    // No layer lands in the forbidden gap (1, 2.5).
    assert!(
        !layers.iter().any(|&l| l > 1.0 + 1e-4 && l < 2.5 - 1e-4),
        "Square: a layer landed in the forbidden gap"
    );
}
