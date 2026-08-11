//! Integration tests for the second-integral range mapping.
//!
//! Exercises the public API of `billiards::second_integral` — that the slider's
//! normalized `[0, 1]` ↔ value mapping is a faithful parametrization of the
//! valid Λ / θ range, on both the ellipse and hyperbola sides.

use billiards::second_integral::{LambdaRange, SecondIntegralRange};

#[test]
fn polyline_is_angle_range() {
    let domain = billiards::domain::square();
    let r = SecondIntegralRange::for_domain(&domain, /* is_confocal */ false);
    assert!(
        matches!(r, SecondIntegralRange::Angle),
        "polyline domain should map to Angle, got something else"
    );
}

#[test]
fn angle_mapping_is_linear() {
    // Exercise the Angle mapping through a real polyline preset.
    let preset = billiards::presets::all_presets()
        .into_iter()
        .find(|p| !p.is_confocal)
        .expect("a polyline preset");
    let r = SecondIntegralRange::for_domain(&preset.domain, false);
    assert!((r.value_at_fraction(0.0) - (-1.0)).abs() < 1e-6);
    assert!((r.value_at_fraction(0.5) - 0.0).abs() < 1e-6);
    assert!((r.value_at_fraction(1.0) - 1.0).abs() < 1e-6);
    assert!((r.fraction_of_value(-1.0) - 0.0).abs() < 1e-6);
    assert!((r.fraction_of_value(0.0) - 0.5).abs() < 1e-6);
    assert!((r.fraction_of_value(1.0) - 1.0).abs() < 1e-6);
}

#[test]
fn confocal_fraction_round_trips() {
    let preset = billiards::presets::all_presets()
        .into_iter()
        .find(|p| p.is_confocal)
        .expect("a confocal preset");
    let range = SecondIntegralRange::for_domain(&preset.domain, true);

    // Ellipse-side midpoints round-trip within slider jitter.
    for t in [0.1f32, 0.25, 0.4] {
        let v = range.value_at_fraction(t);
        let t2 = range.fraction_of_value(v);
        assert!(
            (t - t2).abs() < 1e-4,
            "ellipse side: t={} -> v={} -> t2={}",
            t,
            v,
            t2
        );
    }
    // Hyperbola-side midpoints round-trip.
    for t in [0.6f32, 0.8, 0.95] {
        let v = range.value_at_fraction(t);
        let t2 = range.fraction_of_value(v);
        assert!(
            (t - t2).abs() < 1e-4,
            "hyperbola side: t={} -> v={} -> t2={}",
            t,
            v,
            t2
        );
    }
}

#[test]
fn clamp_jumps_gap_between_sides() {
    let preset = billiards::presets::all_presets()
        .into_iter()
        .find(|p| p.is_confocal)
        .expect("a confocal preset");
    let range = LambdaRange::of_domain(&preset.domain);
    let hyp_min = range.lambda_hyp + 0.05;

    // A value beyond the gap bottom clamps onto the hyperbola side rather than
    // snapping into the forbidden separatrix region.
    let from_above = range.clamp(hyp_min + 0.01, hyp_min);
    assert!(
        from_above >= hyp_min,
        "clamp upward should land on the hyperbola side, got {}",
        from_above
    );
}
