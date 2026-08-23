use macroquad::prelude::vec2;

fn main() {
    let presets = billiards::presets::all_presets();
    let sq = presets
        .iter()
        .find(|p| p.label.starts_with("Square"))
        .unwrap();
    let dom = &sq.domain;
    let cf = billiards::torus::ConfocalParams::standard();

    let lam = 0.5f32;
    let (_t1, _t2, tbl) = billiards::torus::ConfocalStructure::of_domain(dom)
        .unwrap()
        .torus_bounds();
    // (1) wall loop: sample_boundary + velocities_for_lambda
    let mut wall_pts = Vec::new();
    for seg in &dom.segments {
        if let billiards::domain::Segment::Quad { curve, a, b } = seg {
            if (curve.lambda - 0.0).abs() < 1e-4 {
                // outer ellipse arcs only
                wall_pts.extend(billiards::domain::sample_wall_arc(curve, *a, *b, 40));
            }
        }
    }
    // also sample the whole boundary via sample_boundary (includes hyperbola+ellipse)
    let whole = dom.sample_boundary(24);
    println!("elliptic arcs from segments: {}", wall_pts.len());

    for p in &wall_pts {
        for vel in billiards::phase3d::velocities_for_lambda(*p, lam, cf.a, cf.b) {
            let s = billiards::torus::PhaseSample::new(p.x, p.y, vel.x, vel.y);
            let lc = billiards::torus::caustic(&s, &cf);
            let mut params = billiards::torus::TorusParams {
                confocal: &cf,
                lam_wall: bounds.0,
                beta: bounds.1,
                sep_eps: 1e-9,
                cache: &mut cache,
            };
            let (th1, _th2, tidx) = billiards::torus::to_torus(&s, &mut params);
            let tidx = if tidx == u32::MAX { 0 } else { tidx };
            let inside = (lc - lam).abs() < 1e-3;
            println!(
                "p=({:.2},{:.2}) vel=({:.2},{:.2}) lc={:.3} tidx={} inside={}",
                p.x, p.y, vel.x, vel.y, lc, tidx, inside
            );
        }
    }
}
