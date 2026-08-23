//! Full-ellipse preset: the phase manifold is two tori below the focal
//! separatrix and one torus above (the A–B–A molecule).  This pins the
//! "2 tori → 1 torus" transition the app must show.

use billiards::torus::ConfocalParams;

fn ellipse_preset() -> billiards::presets::Preset {
    billiards::presets::all_presets()
        .into_iter()
        .find(|p| p.label.contains("Ellipse (confocal"))
        .expect("full-ellipse preset")
}

/// Below the separatrix (elliptic caustic) the full ellipse splits into two
/// tori (by sign of angular momentum L).
#[test]
fn ellipse_below_separatrix_is_two_tori() {
    let preset = ellipse_preset();
    let lam = 0.5; // 0 < lam < b
    let starts = billiards::dense_caustic_starts(&preset.domain, lam, 8);
    assert!(!starts.is_empty(), "no starts below separatrix");
    let bounds = billiards::confocal::ConfocalStructure::of_domain(&preset.domain)
        .map(|s| s.torus_bounds())
        .unwrap();
    // Full ellipse → beta is None → SplitByL regime → 2 tori.
    assert_eq!(bounds, (0.0, None));
    let regime = billiards::confocal::ConfocalStructure::of_domain(&preset.domain)
        .unwrap()
        .regime(lam);
    assert_eq!(regime, billiards::confocal::TorusRegime::SplitByL);
    assert_eq!(regime.n_tori(), 2);
}

/// Above the separatrix (hyperbolic caustic) the full ellipse is a single
/// torus.
#[test]
fn ellipse_above_separatrix_is_one_torus() {
    let preset = ellipse_preset();
    let lam = 2.5; // b < lam < a
    let regime = billiards::confocal::ConfocalStructure::of_domain(&preset.domain)
        .unwrap()
        .regime(lam);
    assert_eq!(regime, billiards::confocal::TorusRegime::Single);
    assert_eq!(regime.n_tori(), 1);
}

/// The phase manifold is finite (no panic/NaN) at the separatrix and both
/// sides — the "circle" (degenerate) and the 2→1 transition.
#[test]
fn ellipse_phase_is_finite_across_transition() {
    let preset = ellipse_preset();
    let cf = ConfocalParams::standard();
    let bounds = (0.0, None);
    for lam in [0.5f32, cf.b, 1.5, 2.5, 3.5] {
        let starts = billiards::dense_caustic_starts(&preset.domain, lam, 8);
        for (p, v) in starts.into_iter().take(3) {
            let pts = billiards::phase3d::sample_trajectory_phase_dense(
                &preset.domain,
                p,
                v,
                40,
                4,
                bounds,
            );
            for q in &pts {
                assert!(q.theta1.is_finite(), "theta1 NaN at Λ={lam}");
                assert!(q.theta2.is_finite(), "theta2 NaN at Λ={lam}");
            }
        }
    }
}

/// The *rendered* phase manifold must show the 2→1 transition: two distinct
/// torus lobes below the separatrix, one above.  This is what the 3D view
/// draws (`torus_ids`), so it pins requirement 3 ("2 tori → 1 torus in the
/// phase manifold").
#[test]
fn ellipse_rendered_tori_transition() {
    let preset = ellipse_preset();
    let cf = ConfocalParams::standard();
    let structure =
        billiards::confocal::ConfocalStructure::of_domain(&preset.domain).expect("structure");

    fn filled(domain: &billiards::domain::Domain, lam: f32) -> Vec<[f32; 4]> {
        let mut all = Vec::new();
        for (p, v) in billiards::dense_caustic_starts(domain, lam, 24) {
            all.extend(billiards::phase3d::sample_trajectory_phase_full(
                domain, p, v, 400, 12,
            ));
        }
        all
    }

    // Below the separatrix: two distinct rendered tori (ids 0 and 1).
    let below = filled(&preset.domain, 0.5);
    let regime = structure.regime(0.5);
    let ids = billiards::phase3d::torus_ids(&below, regime);
    let n0 = ids.iter().filter(|&&i| i == 0).count();
    let n1 = ids.iter().filter(|&&i| i == 1).count();
    assert!(
        n0 > 0 && n1 > 0 && n0 + n1 == below.len(),
        "full ellipse below b: expected TWO rendered tori, got ids {:?}",
        ids
    );

    // Above the separatrix: a single rendered torus (id 0).
    let above = filled(&preset.domain, 2.5);
    let regime = structure.regime(2.5);
    let ids = billiards::phase3d::torus_ids(&above, regime);
    let n0 = ids.iter().filter(|&&i| i == 0).count();
    assert!(
        n0 == above.len(),
        "full ellipse above b: expected ONE rendered torus, got ids {:?}",
        ids
    );

    // The degenerate separatrix itself is a circle (the two tori pinch onto the
    // focal segment); it must not panic and must be finite.
    let _ = billiards::phase3d::boundary_preimage_points(&preset.domain, cf.b, (0.0, None));
}
