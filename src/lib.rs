pub mod app;
pub mod confocal;
pub mod domain;
pub mod molecule;
pub mod phase3d;
pub mod presets;
pub mod pseudo;
pub mod quadratic;
pub mod render;
pub mod second_integral;
pub mod table;
pub mod torus;
pub mod torus_render;
pub mod ui;

pub use confocal::TorusRegime;

use macroquad::prelude::*;

/// Pick start points for a given domain at a given value of the second integral.
///
/// For confocal billiards (`is_confocal = true`), `lam` is the caustic parameter Λ.
/// We sample the caustic curve Q_Λ = 0 and pick one point in each connected
/// component of caustic ∩ domain.
///
/// For polyline billiards (`is_confocal = false`), `lam` is the normalised angle
/// θ/π ∈ [-1, 1], and `center` is a fixed interior point.  We return a single
/// trajectory start at `center` with velocity direction `θ = lam·π`.
pub fn get_start_points(
    lam: f32,
    dom: &domain::Domain,
    is_confocal: bool,
    center: Vec2,
) -> Vec<(Vec2, Vec2)> {
    if is_confocal {
        let structure = confocal::ConfocalStructure::of_domain(dom)
            .expect("a confocal domain must have a confocal structure");
        confocal::caustic_starts(&structure, dom, lam, confocal::CausticSampling::Sparse, 0)
    } else {
        let angle = lam * std::f32::consts::PI; // lam is θ/π ∈ [-1, 1]
        let v = vec2(angle.cos(), angle.sin());
        vec![(center, v)]
    }
}

/// Dense per-component caustic start points, both ±velocity directions, so the
/// union sweeps out the 2D Liouville torus.
pub fn dense_caustic_starts(
    dom: &domain::Domain,
    lam: f32,
    per_component: usize,
) -> Vec<(Vec2, Vec2)> {
    let structure = confocal::ConfocalStructure::of_domain(dom)
        .expect("a confocal domain must have a confocal structure");
    confocal::caustic_starts(
        &structure,
        dom,
        lam,
        confocal::CausticSampling::Dense,
        per_component,
    )
}

/// One start point per connected component of the caustic ∩ domain.
pub fn start_points_on_caustic(lam: f32, dom: &domain::Domain) -> Vec<(Vec2, Vec2)> {
    let structure = confocal::ConfocalStructure::of_domain(dom)
        .expect("a confocal domain must have a confocal structure");
    confocal::caustic_starts(&structure, dom, lam, confocal::CausticSampling::Sparse, 0)
}

/// Torus-mapping bounds `(lam_wall, beta)` for a domain, in the torus-map
/// convention (full-ellipse fallback for non-quadrilaterals).  Thin wrapper
/// over [`confocal::ConfocalStructure::torus_bounds`].
pub fn torus_bounds(domain: &domain::Domain) -> (f32, Option<f32>) {
    confocal::ConfocalStructure::of_domain(domain)
        .map(|s| s.torus_bounds())
        .unwrap_or((0.0, None))
}
