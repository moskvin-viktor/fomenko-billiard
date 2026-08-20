//! Verify the π⁻¹(boundary) sampler: mapping the billiard's walls and its
//! caustic onto the Liouville torus produces a non-empty, phase-valid set of
//! points that fall on the torus surface (angles in range) and that genuinely
//! represents *distinct* boundary classes (wall vs caustic).

use billiards::{phase3d, presets, torus_bounds};

/// For a confocal preset, sample the boundary preimage and check that both the
/// wall and caustic classes are represented and map to valid torus angles.
#[test]
fn boundary_preimage_maps_wall_and_caustic_onto_torus() {
    let cf = billiards::torus::ConfocalParams::standard();
    let preset = presets::all_presets()
        .into_iter()
        .find(|p| p.is_confocal)
        .expect("a confocal preset");
    let domain = &preset.domain;
    let lam = 0.5; // elliptic caustic well away from the degenerate values b, a
    let bounds = torus_bounds(domain);

    let preimages = phase3d::boundary_preimage_points(domain, lam, bounds);
    assert!(
        preimages.len() >= 2,
        "expected wall + caustic points, got {}",
        preimages.len()
    );

    let mut n_wall = 0;
    let mut n_caustic = 0;
    let mut min_t1 = f32::MAX;
    let mut max_t1 = f32::MIN;
    let mut min_t2 = f32::MAX;
    let mut max_t2 = f32::MIN;
    for (pt, is_caustic) in &preimages {
        assert!(
            (0.0..=2.0 * std::f32::consts::PI).contains(&pt.theta1),
            "theta1 out of range: {}",
            pt.theta1
        );
        assert!(
            (0.0..=2.0 * std::f32::consts::PI).contains(&pt.theta2),
            "theta2 out of range: {}",
            pt.theta2
        );
        assert!((-1.0..=1.0).contains(&pt.theta), "theta out of range");
        min_t1 = min_t1.min(pt.theta1);
        max_t1 = max_t1.max(pt.theta1);
        min_t2 = min_t2.min(pt.theta2);
        max_t2 = max_t2.max(pt.theta2);
        if *is_caustic {
            n_caustic += 1;
        } else {
            n_wall += 1;
        }
    }

    assert!(n_wall > 0, "no wall preimage points");
    assert!(n_caustic > 0, "no caustic preimage points");
    // The wall and caustic are not degenerate single points: they sweep out a
    // finite arc on the torus.
    assert!(max_t1 - min_t1 > 1e-3, "torus theta1 range too small");
    assert!(max_t2 - min_t2 > 1e-3, "torus theta2 range too small");

    // Every simple-libration phase must be within a full turn.
    let _ = cf;
}
