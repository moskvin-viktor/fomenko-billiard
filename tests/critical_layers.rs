//! Critical-layer behavior: trajectories must exist on critical (degenerate)
//! caustic layers, phase manifolds must exist on special layers, and the UI
//! must step through every layer (critical *and* regular) without skipping.
//!
//! A *degenerate* critical layer is one where the phase manifold collapses to a
//! 1D orbit: the ellipse wall `λ_ell`, the focal separatrix `b`, or the focal
//! axis `a`.  The hyperbola wall `λ_hyp` is *not* degenerate — the caustic
//! merely coincides with the wall, the accessible region keeps full area and
//! the level is an ordinary Liouville torus.  Molecule-critical values
//! (genus_jump, edge_swap, birth/death) are likewise *regular* levels handled
//! by the ordinary caustic sampler.

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

/// The degenerate critical values of a domain: `b`, `a`, and the ellipse wall
/// `λ_ell`.  These are where the phase manifold collapses to a 1D orbit.  The
/// hyperbola wall `λ_hyp` is deliberately *excluded*: it is a regular torus
/// level (the caustic coincides with the wall but the level keeps full area).
fn degenerate_layers(preset: &Preset) -> Vec<f32> {
    let cf = cf();
    let s = billiards::confocal::ConfocalStructure::of_domain(&preset.domain).unwrap();
    let mut v = vec![cf.b, cf.a, s.lambda_ell];
    v.retain(|&x| x >= 0.0 && x <= cf.a + 1e-6);
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

/// Every special layer is reachable and, for the degenerate ones, the manifold
/// table returns a collapsed 1D curve (or a flat/empty manifold when the piece
/// misses the domain or the table is pseudo-integrable) — never a full 2D
/// torus.  Molecule-critical regular levels are handled by the ordinary torus
/// path.
#[test]
fn phase_manifold_exists_on_special_layers() {
    use billiards::manifold::{build_manifold, Manifold, SampleConfig};
    let cf = cf();
    for preset in confocal_presets() {
        let domain = &preset.domain;
        let layers = billiards::molecule::special_layers(domain, &cf);
        for &lam in &layers {
            // get_start_points must not panic.
            let _ = billiards::get_start_points(lam, domain, true, Vec2::ZERO);

            // Degenerate layers collapse to closed 1D curves.
            if degenerate_layers(&preset)
                .iter()
                .any(|&d| (d - lam).abs() < 1e-4)
            {
                let m = build_manifold(domain, lam, &SampleConfig::default());
                assert!(
                    !matches!(m, Manifold::Torus2D { .. }),
                    "{}: Λ={} degenerate layer produced a full 2D torus",
                    preset.label,
                    lam
                );
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

/// Degenerate layers collapse to closed 1D circles with finite torus angles.
/// The circle count is read off the torus mapping (never hardcoded): the
/// separatrix `λ = b` is a figure-eight (2 circles), the focal axis `λ = a` a
/// single circle, and wall layers at least one circle per sliding wall.
#[test]
fn degenerate_layers_collapse_to_circles() {
    use billiards::manifold::{build_manifold, Manifold, SampleConfig};
    let cf = cf();
    for preset in confocal_presets() {
        let domain = &preset.domain;
        for &lam in &degenerate_layers(&preset) {
            let m = build_manifold(domain, lam, &SampleConfig::default());
            match m {
                Manifold::Curve1D { circles, pinched } => {
                    assert!(
                        !circles.is_empty(),
                        "{}: Λ={} Curve1D with zero circles",
                        preset.label,
                        lam
                    );
                    for (k, c) in circles.iter().enumerate() {
                        assert!(
                            c.len() >= 8,
                            "{}: Λ={} circle {} too sparse ({} points)",
                            preset.label,
                            lam,
                            k,
                            c.len()
                        );
                        for p in c {
                            assert!(
                                p.theta1.is_finite() && p.theta2.is_finite(),
                                "{}: Λ={} circle {} has non-finite torus angles",
                                preset.label,
                                lam,
                                k
                            );
                        }
                    }
                    if (lam - cf.b).abs() < 1e-4 {
                        assert_eq!(
                            circles.len(),
                            2,
                            "{}: Λ=b separatrix must be a figure-eight (2 circles)",
                            preset.label
                        );
                        assert!(
                            pinched,
                            "{}: Λ=b separatrix circles must be pinched",
                            preset.label
                        );
                    } else {
                        assert!(
                            !pinched,
                            "{}: Λ={} non-separatrix layer marked pinched",
                            preset.label,
                            lam
                        );
                    }
                }
                // The border piece may miss the domain entirely, and pseudo
                // tables classify as flat before the degenerate branch.
                Manifold::Empty | Manifold::FlatSurface { .. } => {}
                Manifold::Torus2D { .. } => panic!(
                    "{}: Λ={} degenerate layer produced a full 2D torus",
                    preset.label, lam
                ),
            }
        }
    }
}

/// Square preset: exact circle topology at each degenerate layer.
/// `λ = λ_ell = 0` → two wall circles (one per torus); `λ = b = 1` →
/// figure-eight (2 circles); `λ = a = 4` → 1 circle (focal axis).  The
/// hyperbola wall `λ_hyp = 2.5` is a *regular* level and must yield a full
/// 2D torus, not a collapsed curve.
#[test]
fn square_degenerate_circle_topology() {
    use billiards::manifold::{build_manifold, Manifold, SampleConfig};
    let preset = confocal_presets()
        .into_iter()
        .find(|p| p.label.starts_with("Square: ellipse"))
        .expect("Square preset");
    let domain = &preset.domain;

    let n_circles = |lam: f32| -> usize {
        match build_manifold(domain, lam, &SampleConfig::default()) {
            Manifold::Curve1D { circles, .. } => circles.len(),
            other => panic!(
                "Square Λ={lam}: expected Curve1D, got {}",
                match other {
                    Manifold::Torus2D { .. } => "Torus2D",
                    Manifold::FlatSurface { .. } => "FlatSurface",
                    Manifold::Empty => "Empty",
                    Manifold::Curve1D { .. } => unreachable!(),
                }
            ),
        }
    };

    assert_eq!(n_circles(0.0), 2, "Square Λ=λ_ell=0: two wall circles");
    assert_eq!(n_circles(1.0), 2, "Square Λ=b: figure-eight");
    assert_eq!(n_circles(4.0), 1, "Square Λ=a: single circle on focal axis");
    assert!(
        matches!(
            build_manifold(domain, 2.5, &SampleConfig::default()),
            Manifold::Torus2D { .. }
        ),
        "Square Λ_hyp=2.5: hyperbola wall is a regular level, expected Torus2D"
    );
}

// ---------------------------------------------------------------------------
// 2b. The torus geometry morphs smoothly onto the degenerate curves
// ---------------------------------------------------------------------------

/// `DomainScale::morph_to_level` reaches the exact degenerate limits and does
/// not jump across the separatrix:
/// - λ → λ_ell and λ → b: `r_minor → 0` (torus thins onto its equator ring);
/// - λ → b from the elliptic side: `gap → 2·r_major` (the two rings touch —
///   figure-eight);
/// - λ → a: `r_minor → 0` (torus thins onto its equator ring; the hole never
///   closes);
/// - mid-band levels are untouched;
/// - approaching b from both sides gives the same `r_minor` limit (0).
#[test]
fn torus_scale_morphs_continuously_to_degenerate_limits() {
    use billiards::torus_render::DomainScale;
    let cf = cf();
    let preset = confocal_presets()
        .into_iter()
        .find(|p| p.label.starts_with("Square: ellipse"))
        .expect("Square preset");
    let s = billiards::confocal::ConfocalStructure::of_domain(&preset.domain).unwrap();
    let base = DomainScale::of_extent(2.3);
    let morph = |lam: f32| DomainScale::of_extent(2.3).morph_to_level(&s, lam);

    // Exact limits.
    let at_b = morph(cf.b);
    assert!(at_b.r_minor.abs() < 1e-6, "Λ=b: r_minor must vanish");
    assert!(
        (at_b.gap - 2.0 * at_b.r_major).abs() < 1e-4,
        "Λ=b: gap must equal 2·r_major so the rings touch (gap={}, r_major={})",
        at_b.gap,
        at_b.r_major
    );
    let at_a = morph(cf.a);
    assert!(at_a.r_minor.abs() < 1e-6, "Λ=a: r_minor must vanish");
    assert!(
        at_a.r_major > 0.0,
        "Λ=a: r_major must survive (equator ring — the hole never closes)"
    );
    let at_ell = morph(s.lambda_ell);
    assert!(at_ell.r_minor.abs() < 1e-6, "Λ=λ_ell: r_minor must vanish");

    // Mid-band untouched.
    let mid = morph(0.5 * (s.lambda_ell + cf.b));
    assert!(
        (mid.r_minor - base.r_minor).abs() < 1e-6 && (mid.r_major - base.r_major).abs() < 1e-6,
        "mid elliptic band must be unmorphed"
    );

    // No jump across the separatrix: both one-sided limits agree with λ=b.
    let eps = 1e-4;
    let below = morph(cf.b - eps);
    let above = morph(cf.b + eps);
    assert!(
        below.r_minor < 1e-2 && above.r_minor < 1e-2,
        "r_minor must be near 0 on both sides of Λ=b (below={}, above={})",
        below.r_minor,
        above.r_minor
    );

    // Smooth ramps: no step bigger than the sweep would explain.  `gap` is
    // only meaningful on the two-torus (elliptic) side, so its continuity is
    // only checked there.
    let mut prev = morph(s.lambda_ell);
    let mut prev_lam = s.lambda_ell;
    let n = 200;
    for i in 1..=n {
        let lam = s.lambda_ell + (cf.a - s.lambda_ell) * i as f32 / n as f32;
        let cur = morph(lam);
        let mut pairs = vec![(prev.r_minor, cur.r_minor), (prev.r_major, cur.r_major)];
        if lam <= cf.b && prev_lam <= cf.b {
            pairs.push((prev.gap, cur.gap));
        }
        for (p, c) in pairs {
            assert!(
                (p - c).abs() < 0.2 * base.r_major,
                "scale jumped at Λ={lam}: {p} -> {c}"
            );
        }
        prev = cur;
        prev_lam = lam;
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

/// The standard first domain (Square): special layers are exactly the three
/// degenerate ones `{λ_ell=0.0, b=1.0, a=4.0}` (the hyperbola wall 2.5 is a
/// regular level with no marker), each has trajectories, and the step list
/// visits the reachable elliptic + hyperbolic bands with regulars in between
/// (never the forbidden gap `(1.0, 2.5)`).
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
        vec![0.0, 1.0, 4.0],
        "Square special layers"
    );

    for &lam in &[0.0f32, 1.0, 4.0] {
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
