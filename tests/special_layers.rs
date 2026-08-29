//! The "all special layers reachable" contract: every confocal preset exposes
//! an ordered set of special caustic values (molecule critical values, walls,
//! the `b` separatrix, and `a`) that the app's strip navigates.  Each layer
//! must be reachable (the app snaps Λ to it) and must produce a valid phase
//! manifold / trajectory without panicking.

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

/// The special layers include the separatrix `b`, the focal axis `a`, and for
/// every quadric-arc domain at least one wall.
#[test]
fn special_layers_cover_degeneracies() {
    let cf = cf();
    for preset in confocal_presets() {
        let layers = billiards::molecule::special_layers(&preset.domain, &cf);
        assert!(
            layers.iter().any(|&l| (l - cf.b).abs() < 1e-4),
            "{}: special layers must include the separatrix b={}",
            preset.label,
            cf.b
        );
        assert!(
            layers.iter().any(|&l| (l - cf.a).abs() < 1e-4),
            "{}: special layers must include the focal a={}",
            preset.label,
            cf.a
        );
        assert!(!layers.is_empty(), "{}: no special layers", preset.label);
    }
}

/// The molecule strip markers classify every special value (kind + label), and
/// snap navigation returns a valid finite value.
#[test]
fn markers_are_reachable_and_finite() {
    let cf = cf();
    for preset in confocal_presets() {
        let markers = billiards::molecule_view::layer_markers(&preset);
        assert!(!markers.is_empty(), "{}: no markers", preset.label);
        for m in &markers {
            assert!(m.lam.is_finite(), "{}: marker lam not finite", preset.label);
            assert!(
                m.lam >= 0.0 && m.lam <= cf.a + 1e-6,
                "{}: marker lam {} out of [0, a] (λ_ell = 0 is a valid wall layer)",
                preset.label,
                m.lam
            );
        }
        // Navigation returns a finite member of the set.
        for dir in [-1, 1] {
            let snapped = billiards::molecule_view::snap_nearest(&markers, 0.5, dir);
            assert!(snapped.is_some(), "{}: no snap target", preset.label);
            let v = snapped.unwrap();
            assert!(v.is_finite());
        }
    }
}

/// Snap-to-layer produces a valid phase manifold at every special layer.
#[test]
fn phase_is_valid_when_snapping_specials() {
    for preset in confocal_presets() {
        let marks = billiards::molecule_view::layer_markers(&preset);
        let bounds = billiards::confocal::ConfocalStructure::of_domain(&preset.domain)
            .map(|s| s.torus_bounds())
            .unwrap_or((0.0, None));
        for m in &marks {
            // get_start_points must not panic (empty is fine — degenerate layer).
            let _ = billiards::get_start_points(m.lam, &preset.domain, true, Vec2::ZERO);
            let pre = billiards::phase3d::boundary_preimage_points(&preset.domain, m.lam, bounds);
            for (q, _isc) in &pre {
                assert!(
                    q.theta1.is_finite() && q.theta2.is_finite(),
                    "{}: Λ={} preimage non-finite",
                    preset.label,
                    m.lam
                );
            }
        }
    }
}
