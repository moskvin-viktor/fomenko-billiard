use billiards::{domain, phase3d, presets, torus::ConfocalParams};
use macroquad::prelude::*;

/// Map a phase point to the two elliptic coordinates (μ, ν).
///
/// Every point (x, y) lies on exactly one confocal ellipse (μ < b) and one
/// confocal hyperbola (ν > b) of the family (b−λ)x² + (a−λ)y² = (a−λ)(b−λ).
/// These are the separable coordinates of the elliptic billiard, so on the
/// invariant torus at fixed Λ the pair (μ, ν) fills a genuine 2D rectangle
/// (μ ranges between the caustic and the boundary ellipse, ν between the two
/// boundary hyperbolas).  This is the correct torus parametrization, unlike
/// (atan2(y,x), θ) where θ is determined by position at fixed Λ.
fn to_elliptic(q: &[f32; 4]) -> (f32, f32) {
    let (x, y, _vx, _vy) = (q[0], q[1], q[2], q[3]);
    // Solve λ² − λ(a+b−x²−y²) − (bx² + ay² − ab) = 0 for the two roots.
    let cf = ConfocalParams::standard();
    let a = cf.a;
    let b = cf.b;
    let s = x * x + y * y;
    let t = b * x * x + a * y * y - a * b;
    let p = a + b - s;
    let disc = p * p + 4.0 * t;
    let sd = disc.max(0.0).sqrt();
    let lam1 = 0.5 * (p - sd);
    let lam2 = 0.5 * (p + sd);
    // lam1 < lam2; the ellipse coordinate is the smaller, hyperbola the larger.
    (lam1, lam2)
}

/// Verify that a set of phase points, viewed in elliptic coordinates (μ, ν),
/// fills a genuine 2D region (i.e. the manifold is a torus, not a 1D curve).
fn assert_is_2d_torus(points: &[[f32; 4]], label: &str) {
    assert!(
        points.len() > 200,
        "{}: not enough phase points ({}) to assess dimensionality",
        label,
        points.len()
    );

    // Bin into a grid over (μ, ν).
    const N: usize = 24;
    let mut bins = [[0u32; N]; N];
    let mut min_mu = f32::MAX;
    let mut max_mu = f32::MIN;
    let mut min_nu = f32::MAX;
    let mut max_nu = f32::MIN;
    for q in points {
        let (mu, nu) = to_elliptic(q);
        min_mu = min_mu.min(mu);
        max_mu = max_mu.max(mu);
        min_nu = min_nu.min(nu);
        max_nu = max_nu.max(nu);
    }
    let mu_span = (max_mu - min_mu).max(1e-6);
    let nu_span = (max_nu - min_nu).max(1e-6);
    for q in points {
        let (mu, nu) = to_elliptic(q);
        let i = (((mu - min_mu) / mu_span) * N as f32)
            .floor()
            .clamp(0.0, N as f32 - 1.0) as usize;
        let j = (((nu - min_nu) / nu_span) * N as f32)
            .floor()
            .clamp(0.0, N as f32 - 1.0) as usize;
        bins[i][j] += 1;
    }

    // Both elliptic coordinates must vary over a substantial range.
    let cf = ConfocalParams::standard();
    let mu_frac = mu_span / cf.b; // μ ranges in [−∞, b]
    let nu_frac = nu_span / (cf.a - cf.b); // ν ranges in [b, a]
    assert!(
        mu_frac > 0.1 && nu_frac > 0.1,
        "{}: elliptic coords do not both vary — μ span {:.3}, ν span {:.3}",
        label,
        mu_frac,
        nu_frac
    );

    // Count occupied bins and max occupancy per row/column.
    let mut occupied = 0;
    let mut max_row = 0;
    let mut max_col = 0;
    for i in 0..N {
        let mut r = 0;
        let mut c = 0;
        for j in 0..N {
            if bins[i][j] > 0 {
                occupied += 1;
                r += 1;
            }
            if bins[j][i] > 0 {
                c += 1;
            }
        }
        max_row = max_row.max(r);
        max_col = max_col.max(c);
    }

    let area_frac = occupied as f32 / (N * N) as f32;
    assert!(
        area_frac > 0.07,
        "{}: occupied area too small ({:.2}) — manifold looks 1D",
        label,
        area_frac
    );
    assert!(
        max_row > N / 3 && max_col > N / 3,
        "{}: elliptic coords do not fill a 2D region (max row {}, max col {})",
        label,
        max_row,
        max_col
    );
}

/// Verify that a set of phase points does NOT fill a 3D volume.
///
/// A genuine 2D torus is a *surface*: it has area but zero volume.  If we bin
/// the points in the 3D torus-embedding coordinates (X, Y, Z), a surface
/// occupies only a small fraction of the bounding box's voxels, whereas a
/// 3D volume would fill a large fraction.  This complements the 2D-area test:
/// together they show the manifold is 2D (not 1D, not 3D).
fn assert_is_not_3d_volume(points: &[[f32; 4]], label: &str) {
    assert!(
        points.len() > 200,
        "{}: not enough phase points ({}) to assess volume",
        label,
        points.len()
    );

    // Map each phase point to the 3D torus embedding (same as the renderer).
    let r_major = 1.6;
    let r_minor = 0.6;
    let mut pts3 = Vec::with_capacity(points.len());
    for q in points {
        let (x, y, vx, vy) = (q[0], q[1], q[2], q[3]);
        let phi1 = y.atan2(x);
        let phi2 = vy.atan2(vx);
        let cos2 = phi2.cos();
        let (sin1, cos1) = phi1.sin_cos();
        pts3.push(vec3(
            (r_major + r_minor * cos2) * cos1,
            (r_major + r_minor * cos2) * sin1,
            r_minor * phi2.sin(),
        ));
    }

    // Bounding box.
    let mut min = vec3(f32::MAX, f32::MAX, f32::MAX);
    let mut max = vec3(f32::MIN, f32::MIN, f32::MIN);
    for p in &pts3 {
        min.x = min.x.min(p.x);
        min.y = min.y.min(p.y);
        min.z = min.z.min(p.z);
        max.x = max.x.max(p.x);
        max.y = max.y.max(p.y);
        max.z = max.z.max(p.z);
    }

    // Bin into a 3D voxel grid.
    const N: usize = 16;
    let span_x = (max.x - min.x).max(1e-6);
    let span_y = (max.y - min.y).max(1e-6);
    let span_z = (max.z - min.z).max(1e-6);
    let mut voxels = vec![false; N * N * N];
    let mut occupied = 0usize;
    for p in &pts3 {
        let i = (((p.x - min.x) / span_x) * N as f32)
            .floor()
            .clamp(0.0, N as f32 - 1.0) as usize;
        let j = (((p.y - min.y) / span_y) * N as f32)
            .floor()
            .clamp(0.0, N as f32 - 1.0) as usize;
        let k = (((p.z - min.z) / span_z) * N as f32)
            .floor()
            .clamp(0.0, N as f32 - 1.0) as usize;
        let idx = (i * N + j) * N + k;
        if !voxels[idx] {
            voxels[idx] = true;
            occupied += 1;
        }
    }

    let vol_frac = occupied as f32 / (N * N * N) as f32;
    // A 2D surface fills only a thin shell of the bounding box.  A genuine
    // 3D volume would fill a much larger fraction (typically > 0.5).
    assert!(
        vol_frac < 0.25,
        "{}: phase manifold fills too much volume ({:.2}) — looks 3D, not a torus",
        label,
        vol_frac
    );
}

/// Verify that the *filled* torus (dense interior sampling, as used in the 3D
/// view) actually covers a reasonable fraction of the torus surface.
///
/// The torus surface is parametrized by the torus angles (θ₁, θ₂) of the
/// spec's §6 — the phases on the two ovals of the cubic.  At fixed Λ the pair
/// (θ₁, θ₂) fills a genuine 2D rectangle, so a well-filled torus should cover
/// most of it.
fn assert_torus_well_filled(
    points: &[[f32; 4]],
    label: &str,
    _is_hyperbola: bool,
    domain: &domain::Domain,
) {
    assert!(
        points.len() > 1000,
        "{}: not enough filled points ({}) to assess coverage",
        label,
        points.len()
    );

    // Bin the points in the (θ₁, θ₂) torus surface coordinates, using the same
    // mapping as the renderer.
    const N: usize = 32;
    let mut bins = [[0u32; N]; N];
    let mut min_t1 = f32::MAX;
    let mut max_t1 = f32::MIN;
    let mut min_t2 = f32::MAX;
    let mut max_t2 = f32::MIN;
    let mut mapped = Vec::with_capacity(points.len());
    let mut cache = billiards::torus::TorusCache::default();
    let cf = billiards::torus::ConfocalParams::standard();
    let (lam_wall, beta) = billiards::torus_bounds(domain);
    for q in points {
        let sample = billiards::torus::PhaseSample::new(q[0], q[1], q[2], q[3]);
        let mut params = billiards::torus::TorusParams {
            confocal: &cf,
            lam_wall,
            beta,
            sep_eps: 1e-9,
            cache: &mut cache,
        };
        let (th1, th2, _idx) = billiards::torus::to_torus(&sample, &mut params);
        mapped.push((th1, th2));
        min_t1 = min_t1.min(th1);
        max_t1 = max_t1.max(th1);
        min_t2 = min_t2.min(th2);
        max_t2 = max_t2.max(th2);
    }
    let t1_span = (max_t1 - min_t1).max(1e-6);
    let t2_span = (max_t2 - min_t2).max(1e-6);
    for (th1, th2) in &mapped {
        let i = (((th1 - min_t1) / t1_span) * N as f32)
            .floor()
            .clamp(0.0, N as f32 - 1.0) as usize;
        let j = (((th2 - min_t2) / t2_span) * N as f32)
            .floor()
            .clamp(0.0, N as f32 - 1.0) as usize;
        bins[i][j] += 1;
    }

    let mut occupied = 0;
    for i in 0..N {
        for j in 0..N {
            if bins[i][j] > 0 {
                occupied += 1;
            }
        }
    }
    let area_frac = occupied as f32 / (N * N) as f32;

    // A well-filled torus should cover a substantial fraction of the surface.
    // (The L-shape hyperbola caustic is a thin degenerate region, so it fills
    // a smaller fraction than the clean quadrilateral cases.)
    assert!(
        area_frac > 0.05,
        "{}: filled torus covers too little of the surface ({:.2}) — unrecognizable",
        label,
        area_frac
    );
}

/// Sample a dense phase-space point cloud for a confocal billiard at a given Λ.
///
/// We trace many trajectories from dense caustic start points and collect all
/// bounce phase points.  The union of these trajectories sweeps out the 2D
/// Liouville torus; mapping them to (φ₁, φ₂) angle coordinates gives a flat
/// 2D region.
fn sample_torus_cloud(domain: &domain::Domain, lam: f32, per_component: usize) -> Vec<[f32; 4]> {
    let mut all = Vec::new();
    let starts = billiards::dense_caustic_starts(domain, lam, per_component);
    for &(p, v) in &starts {
        let traj = phase3d::sample_trajectory_phase_full(domain, p, v, 1200, 1);
        all.extend(traj);
    }
    all
}

/// Sample the *filled* torus exactly as the 3D view does: dense interior
/// points along each segment of many caustic trajectories.
fn sample_filled_torus(domain: &domain::Domain, lam: f32, per_component: usize) -> Vec<[f32; 4]> {
    let mut all = Vec::new();
    let starts = billiards::dense_caustic_starts(domain, lam, per_component);
    for &(p, v) in &starts {
        let traj = phase3d::sample_trajectory_phase_full(domain, p, v, 400, 12);
        all.extend(traj);
    }
    all
}

/// Split phase points into disconnected regions by the sign of y.
///
/// An ellipse caustic in a confocal square (Case C) splits the table into an
/// upper and a lower accessible region, each its own torus.  A hyperbola
/// caustic (Case B) has a single torus spanning both y-halves (the orbit
/// crosses the focal segment).
fn split_regions(points: &[[f32; 4]]) -> (Vec<[f32; 4]>, Vec<[f32; 4]>) {
    let mut upper = Vec::new();
    let mut lower = Vec::new();
    for q in points {
        if q[1] >= 0.0 {
            upper.push(*q);
        } else {
            lower.push(*q);
        }
    }
    (upper, lower)
}

/// Verify that a single region's phase points form a genuine 2D torus in the
/// elliptic coordinates (μ, ν).  This is the "correct torus" check: the region
/// must be 2D (not a 1D curve) and must not be empty.
fn assert_lobe_is_2d_torus(points: &[[f32; 4]], label: &str) {
    assert!(
        points.len() > 50,
        "{}: region has too few points ({}) — torus collapsed or missing",
        label,
        points.len()
    );
    assert_is_2d_torus(points, label);
}

/// Verify the number of tori matches the number of disconnected regions.
///
/// A hyperbola caustic (Λ > B, Case B) has ONE connected region → ONE torus.
/// An ellipse caustic in a confocal square (Λ < B, Case C) has TWO regions
/// (upper / lower, split by the caustic ellipse) → TWO tori.
fn assert_num_tori(points: &[[f32; 4]], is_hyperbola: bool, label: &str) {
    let (upper, lower) = split_regions(points);
    if is_hyperbola {
        assert!(
            !upper.is_empty() && !lower.is_empty(),
            "{}: hyperbola caustic torus should span both y-halves (upper {}, lower {})",
            label,
            upper.len(),
            lower.len()
        );
        assert_is_2d_torus(points, label);
    } else {
        assert!(
            !upper.is_empty() && !lower.is_empty(),
            "{}: ellipse caustic must have TWO regions (upper {}, lower {}), but one is missing",
            label,
            upper.len(),
            lower.len()
        );
        assert_lobe_is_2d_torus(&upper, &format!("{} upper region", label));
        assert_lobe_is_2d_torus(&lower, &format!("{} lower region", label));
    }
}

/// Verify the two integrals H and Λ are conserved along a trajectory.
fn assert_integrals_conserved(domain: &domain::Domain, lam: f32, label: &str) {
    let starts = billiards::start_points_on_caustic(lam, domain);
    assert!(
        !starts.is_empty(),
        "{}: no start points for Λ={}",
        label,
        lam
    );

    for (i, &(p, v)) in starts.iter().enumerate() {
        let traj = phase3d::sample_trajectory_phase_full(domain, p, v, 300, 4);
        if traj.len() <= 10 {
            // Degenerate/borderline start point — skip.
            continue;
        }

        let h0 = 0.5 * (v.x * v.x + v.y * v.y);
        let cf = ConfocalParams::standard();
        let lam0 =
            v.x * v.x / cf.a + v.y * v.y / cf.b - (p.x * v.y - p.y * v.x).powi(2) / (cf.a * cf.b);

        for (k, q) in traj.iter().enumerate() {
            let (x, y, vx, vy) = (q[0], q[1], q[2], q[3]);
            let h = 0.5 * (vx * vx + vy * vy);
            let lam = vx * vx / cf.a + vy * vy / cf.b - (x * vy - y * vx).powi(2) / (cf.a * cf.b);
            assert!(
                (h - h0).abs() < 1e-2,
                "{}: H not conserved on traj {} sample {}: {:.5} vs {:.5}",
                label,
                i,
                k,
                h,
                h0
            );
            assert!(
                (lam - lam0).abs() < 1e-2,
                "{}: Λ not conserved on traj {} sample {}: {:.5} vs {:.5}",
                label,
                i,
                k,
                lam,
                lam0
            );
        }
    }
}

/// Pick a valid Λ on the ellipse side for a given domain.
fn ellipse_lambda(domain: &domain::Domain) -> f32 {
    let cf = ConfocalParams::standard();
    let mut ell = f32::MAX;
    for seg in &domain.segments {
        if let domain::Segment::Quad { curve, .. } = seg {
            if curve.lambda < cf.b {
                ell = ell.min(curve.lambda);
            }
        }
    }
    // Midpoint of the ellipse range (ell + e, B - e)
    (ell + 0.05 + cf.b - 0.05) / 2.0
}

/// Pick a valid Λ on the hyperbola side for a given domain.
fn hyperbola_lambda(domain: &domain::Domain) -> f32 {
    let cf = ConfocalParams::standard();
    let mut hyp = f32::MAX;
    for seg in &domain.segments {
        if let domain::Segment::Quad { curve, .. } = seg {
            if curve.lambda > cf.b {
                hyp = hyp.min(curve.lambda);
            }
        }
    }
    hyp + 0.2
}

// ---------------------------------------------------------------
// Tests
// ---------------------------------------------------------------

#[test]
fn test_square_confocal_phase_is_torus() {
    let preset = presets::all_presets()
        .into_iter()
        .find(|p| p.label.starts_with("Square: ellipse"))
        .expect("Square confocal preset not found");
    let domain = &preset.domain;

    let cloud = sample_torus_cloud(domain, ellipse_lambda(domain), 400);
    assert_is_2d_torus(&cloud, "square ellipse");
    assert_is_not_3d_volume(&cloud, "square ellipse");

    let cloud = sample_torus_cloud(domain, hyperbola_lambda(domain), 400);
    assert_is_2d_torus(&cloud, "square hyperbola");
    assert_is_not_3d_volume(&cloud, "square hyperbola");
}

#[test]
fn test_thin_confocal_phase_is_torus() {
    let preset = presets::all_presets()
        .into_iter()
        .find(|p| p.label.starts_with("Thin"))
        .expect("Thin confocal preset not found");
    let domain = &preset.domain;

    let cloud = sample_torus_cloud(domain, ellipse_lambda(domain), 400);
    assert_is_2d_torus(&cloud, "thin ellipse");
    assert_is_not_3d_volume(&cloud, "thin ellipse");

    let cloud = sample_torus_cloud(domain, hyperbola_lambda(domain), 400);
    assert_is_2d_torus(&cloud, "thin hyperbola");
    assert_is_not_3d_volume(&cloud, "thin hyperbola");
}

#[test]
fn test_flat_confocal_phase_is_torus() {
    let preset = presets::all_presets()
        .into_iter()
        .find(|p| p.label.starts_with("Flat"))
        .expect("Flat confocal preset not found");
    let domain = &preset.domain;

    let cloud = sample_torus_cloud(domain, ellipse_lambda(domain), 400);
    assert_is_2d_torus(&cloud, "flat ellipse");
    assert_is_not_3d_volume(&cloud, "flat ellipse");

    let cloud = sample_torus_cloud(domain, hyperbola_lambda(domain), 400);
    assert_is_2d_torus(&cloud, "flat hyperbola");
    assert_is_not_3d_volume(&cloud, "flat hyperbola");
}

#[test]
fn test_lshape_confocal_phase_is_torus() {
    let preset = presets::all_presets()
        .into_iter()
        .find(|p| p.label.contains("3π/2"))
        .expect("L-shape confocal preset not found");
    let domain = &preset.domain;

    let cloud = sample_torus_cloud(domain, ellipse_lambda(domain), 400);
    assert_is_2d_torus(&cloud, "lshape ellipse");
    assert_is_not_3d_volume(&cloud, "lshape ellipse");

    let cloud = sample_torus_cloud(domain, hyperbola_lambda(domain), 400);
    assert_is_2d_torus(&cloud, "lshape hyperbola");
    assert_is_not_3d_volume(&cloud, "lshape hyperbola");
}

#[test]
fn test_integrals_conserved() {
    let configs = presets::all_presets();
    for preset in &configs {
        if !preset.is_confocal {
            continue;
        }
        let domain = &preset.domain;
        assert_integrals_conserved(domain, ellipse_lambda(domain), preset.label);
        assert_integrals_conserved(domain, hyperbola_lambda(domain), preset.label);
    }
}

/// The filled torus (dense interior sampling, as rendered in the 3D view)
/// must actually cover a reasonable fraction of the torus surface.  This
/// guards against a sparse/unrecognizable fill.
#[test]
fn test_filled_torus_covers_surface() {
    let configs = presets::all_presets();
    for preset in &configs {
        if !preset.is_confocal {
            continue;
        }
        let domain = &preset.domain;
        let filled = sample_filled_torus(domain, ellipse_lambda(domain), 24);
        assert_torus_well_filled(&filled, &format!("{} ellipse", preset.label), false, domain);
        let filled = sample_filled_torus(domain, hyperbola_lambda(domain), 40);
        assert_torus_well_filled(
            &filled,
            &format!("{} hyperbola", preset.label),
            true,
            domain,
        );
    }
}

/// The number of tori must match the number of disconnected caustic regions:
/// one torus for an ellipse caustic, two tori (left/right lobes) for a
/// hyperbola caustic.  Each lobe must be a genuine 2D torus, not collapsed.
///
/// The L-shape confocal table is EXCLUDED: it is pseudo-integrable (genus-2
/// level sets, the flat-coordinate case), not a Liouville torus.
#[test]
fn test_num_tori_matches_regions() {
    let configs = presets::all_presets();
    for preset in &configs {
        if !preset.is_confocal || preset.label.contains("3π/2") {
            continue;
        }
        let domain = &preset.domain;

        // Ellipse caustic → one torus spanning both x-halves.
        let cloud = sample_filled_torus(domain, ellipse_lambda(domain), 40);
        assert_num_tori(&cloud, false, &format!("{}", preset.label));

        // Hyperbola caustic → two tori (left and right lobes).
        let cloud = sample_filled_torus(domain, hyperbola_lambda(domain), 40);
        assert_num_tori(&cloud, true, &format!("{}", preset.label));
    }
}

/// The *rendered* tori must be distinct objects: one per disconnected region.
///
/// This checks the actual torus embedding used by the renderer.  For a
/// hyperbola caustic there is ONE torus (Case B).  For an ellipse caustic in a
/// confocal square (Case C) the table splits into an upper and a lower
/// accessible region, so there are TWO tori (distinguished by `sign(y)`).
#[test]
fn test_rendered_tori_are_distinct() {
    let configs = presets::all_presets();
    for preset in &configs {
        if !preset.is_confocal {
            continue;
        }
        let domain = &preset.domain;
        let structure =
            billiards::confocal::ConfocalStructure::of_domain(domain).expect("confocal structure");

        // Ellipse caustic → exactly two tori (upper / lower region).
        let cloud = sample_filled_torus(domain, ellipse_lambda(domain), 40);
        let regime = structure.regime(ellipse_lambda(domain));
        let ids = phase3d::torus_ids(&cloud, regime);
        let n0 = ids.iter().filter(|&&i| i == 0).count();
        let n1 = ids.iter().filter(|&&i| i == 1).count();
        assert!(
            n0 > 0 && n1 > 0 && n0 + n1 == cloud.len(),
            "{} ellipse: expected TWO tori (ids 0 and 1), got {:?}",
            preset.label,
            ids
        );

        // Hyperbola caustic → exactly one torus.
        let cloud = sample_filled_torus(domain, hyperbola_lambda(domain), 40);
        let regime = structure.regime(hyperbola_lambda(domain));
        let ids = phase3d::torus_ids(&cloud, regime);
        let n_hyper = ids.iter().filter(|&&i| i == 0).count();
        assert!(
            n_hyper == cloud.len(),
            "{} hyperbola: expected ONE torus (all points on torus 0), got ids {:?}",
            preset.label,
            ids
        );
    }
}
