use macroquad::prelude::*;

/// A 3D phase space point (x, y, θ/π) where θ = atan2(vy, vx)
/// is the velocity direction angle, normalized to [-1, 1].
///
/// The torus angles `theta1`, `theta2` and `torus_index` are precomputed at
/// sampling time (once per trajectory, sharing a `TorusCache`) so the renderer
/// never rebuilds the expensive quadrature splines per point.
#[derive(Clone, Copy, Debug)]
pub struct PhasePoint {
    pub x: f32,
    pub y: f32,
    pub theta: f32,  // θ/π ∈ [-1, 1]
    pub theta1: f32, // torus phase on the first oval, ∈ [0, 2π)
    pub theta2: f32, // torus phase on the second oval, ∈ [0, 2π)
    pub torus_index: u32,
}

impl PhasePoint {
    /// Unit velocity direction `(cos θ, sin θ)` reconstructed from `theta`.
    /// Unit speed is enforced for free, matching the spec's convention.
    pub fn velocity(&self) -> Vec2 {
        let th = self.theta * std::f32::consts::PI;
        vec2(th.cos(), th.sin())
    }
}

// (The old `CriticalSheet` / `critical_phase_sheet` / `build_degenerate_orbit`
// hardcoded a single xz-plane circle for every degenerate layer.  Degenerate
// layers are now built by `crate::manifold::build_manifold` as
// `Manifold::Curve1D` — circles read off the torus mapping — and drawn by
// [`draw_curve1d`] below.)

/// The stride for drawing a cloud of `n` points under a `budget` cap, so the
/// drawn set `0, step, 2·step, …` has at most `budget` elements.
///
/// This bounds the per-frame `draw_circle` cost when orbiting (the profiled
/// ~150K-point wall): larger clouds are decimated, small clouds keep every
/// point (`step = 1`).
///
/// Returns the smallest `step ≥ 1` with `⌊n/step⌋ + 1 ≤ budget`.
pub fn decimation_step(n: usize, budget: usize) -> usize {
    assert!(budget >= 1, "budget must be >= 1");
    if n == 0 || n <= budget {
        return 1;
    }
    // Need ⌊n/step⌋ + 1 ≤ budget  ⇒  n/step < budget  ⇒  step > n/budget.
    // The smallest such integer is ⌊n/budget⌋ + 1.
    (n / budget).max(1) + 1
}

// ----------------------------------------------------------------
// 3D orbit camera
// ----------------------------------------------------------------
pub struct OrbitCamera3 {
    /// Rotation around vertical axis (yaw), in radians.
    azimuth: f32,
    /// Rotation around horizontal axis (pitch), in radians.
    elevation: f32,
    /// Distance from origin.
    distance: f32,
    /// Last mouse position for drag tracking.
    last_mouse: Option<Vec2>,
    /// Whether the user is currently dragging.
    dragging: bool,
}

impl OrbitCamera3 {
    pub fn new() -> Self {
        Self {
            azimuth: -0.6,
            elevation: 0.5,
            distance: 3.2,
            last_mouse: None,
            dragging: false,
        }
    }

    /// Handle mouse input for orbit control. Call once per frame.
    pub fn handle_input(&mut self) {
        let (mx, my) = (mouse_position().0, mouse_position().1);
        let mouse = vec2(mx, my);

        // Scroll to zoom
        let scroll = mouse_wheel();
        self.distance = (self.distance - scroll.1 * 0.3).max(1.5).min(20.0);

        // Right-click drag to orbit
        if is_mouse_button_pressed(MouseButton::Right) {
            self.dragging = true;
            self.last_mouse = Some(mouse);
        }
        if is_mouse_button_released(MouseButton::Right) {
            self.dragging = false;
            self.last_mouse = None;
        }

        if self.dragging {
            if let Some(last) = self.last_mouse {
                let dx = mouse.x - last.x;
                let dy = mouse.y - last.y;
                self.azimuth -= dx * 0.008;
                self.elevation = (self.elevation + dy * 0.008).clamp(-1.4, 1.4);
            }
            self.last_mouse = Some(mouse);
        }
    }

    /// Project a 3D point (world coordinates) to screen coordinates.
    /// Returns (screen_x, screen_y, depth) where depth is used for z-sorting.
    pub fn project(&self, p: Vec3, win_w: f32, win_h: f32) -> Vec3 {
        // Camera position on the sphere
        let cam_x = self.distance * self.elevation.cos() * self.azimuth.sin();
        let cam_y = self.distance * self.elevation.sin();
        let cam_z = self.distance * self.elevation.cos() * self.azimuth.cos();

        // Look-at vector (from camera to origin)
        let look = vec3(-cam_x, -cam_y, -cam_z).normalize();
        // Right vector (cross of look and world up)
        let world_up = vec3(0.0, 1.0, 0.0);
        let right = look.cross(world_up).normalize();
        // Up vector
        let up = right.cross(look).normalize();

        // Vector from camera to point
        let d = p - vec3(cam_x, cam_y, cam_z);

        // Project onto camera basis
        let sx = d.dot(right);
        let sy = d.dot(up);
        let sz = d.dot(look); // depth (negative = in front of camera)

        // Perspective projection
        let fov = 1.5;
        let near = 0.1;
        let scale = fov / (sz + near).max(0.1);

        let screen_x = win_w / 2.0 + sx * scale * win_h / 2.0;
        let screen_y = win_h / 2.0 - sy * scale * win_h / 2.0;

        vec3(screen_x, screen_y, sz)
    }

    /// A cheap fingerprint of the camera state, used to detect when the cached
    /// torus render must be regenerated (orbit/zoom changed).
    pub fn state_key(&self) -> (f32, f32, f32) {
        (self.azimuth, self.elevation, self.distance)
    }
}

// ----------------------------------------------------------------
// Phase space sampling
// ----------------------------------------------------------------

/// Sample a single trajectory densely, recording interior points along each
/// segment (not just bounce points), together with the precomputed torus
/// angles `(θ₁, θ₂, index)`.
///
/// This is what makes the torus surface look filled in the 3D view: sampling
/// `per_seg` points along every segment of many trajectories sweeps out the
/// whole 2D surface instead of leaving sparse dots at the bounces.
///
/// The torus mapping uses one shared `TorusCache` for the whole trajectory
/// (the spec's "precompute per torus, not per sample" rule) so the expensive
/// quadrature splines are built once, not per point.
pub fn sample_trajectory_phase_dense(
    domain: &crate::domain::Domain,
    p0: Vec2,
    v0: Vec2,
    max_steps: usize,
    per_seg: usize,
    bounds: (f32, Option<f32>),
) -> Vec<PhasePoint> {
    let mut pts = Vec::new();
    let mut p = p0;
    let mut v = v0;
    let mut cache = crate::torus::TorusCache::default();
    let cf = crate::torus::ConfocalParams::standard();
    for _ in 0..max_steps {
        let speed = v.length();
        if speed < 1e-12 {
            break;
        }
        let dir = v / speed;
        let (_t, idx, hit) = match domain.intersect(p, dir) {
            Some(r) => r,
            None => break,
        };
        for k in 0..per_seg {
            let f = k as f32 / per_seg as f32;
            let q = p + (hit - p) * f;
            let theta = v.y.atan2(v.x) / std::f32::consts::PI;
            let sample = crate::torus::PhaseSample::new(q.x, q.y, v.x, v.y);
            let mut params = crate::torus::TorusParams {
                confocal: &cf,
                lam_wall: bounds.0,
                beta: bounds.1,
                sep_eps: 1e-9,
                cache: &mut cache,
            };
            let (th1, th2, tidx) = crate::torus::to_torus(&sample, &mut params);
            let tidx = if tidx == u32::MAX { 0 } else { tidx };
            pts.push(PhasePoint {
                x: q.x,
                y: q.y,
                theta,
                theta1: th1,
                theta2: th2,
                torus_index: tidx,
            });
        }
        v = domain.reflect(hit, v, idx);
        p = hit + 1e-4 * v.normalize();
    }
    pts
}

/// Sample the full 4D phase space (x, y, vx, vy) along a trajectory,
/// including interior points between bounces.
///
/// This is used to verify the invariant manifold is a 2D torus: in the full
/// 4D phase space the two integrals H and Λ confine the motion to a 2D surface,
/// so a dense point cloud sampled here should have 2 large PCA eigenvalues and
/// 2 small ones.
pub fn sample_trajectory_phase_full(
    domain: &crate::domain::Domain,
    p0: Vec2,
    v0: Vec2,
    max_steps: usize,
    per_seg: usize,
) -> Vec<[f32; 4]> {
    let mut pts = Vec::new();
    let mut p = p0;
    let mut v = v0;
    for _ in 0..max_steps {
        let speed = v.length();
        if speed < 1e-12 {
            break;
        }
        let dir = v / speed;
        let (_t, idx, hit) = match domain.intersect(p, dir) {
            Some(r) => r,
            None => break,
        };
        // Sample interior points along this segment.
        for k in 0..per_seg {
            let f = k as f32 / per_seg as f32;
            let q = p + (hit - p) * f;
            pts.push([q.x, q.y, v.x, v.y]);
        }
        v = domain.reflect(hit, v, idx);
        p = hit + 1e-4 * v.normalize();
    }
    pts
}

/// Find the velocity direction(s) at position `pos` that give the target
/// second integral `lam` (with |v| = 1).
///
/// For a fixed Λ and position, there are exactly two unit velocity directions
/// (the two tangents to the confocal caustic through that point).  Sampling
/// both directions across many positions fills out the full 2D Liouville torus,
/// which is what the PCA dimensionality test needs.
pub fn velocities_for_lambda(pos: Vec2, lam: f32, a: f32, b: f32) -> Vec<Vec2> {
    let mut out = Vec::new();
    let n = 720;
    let mut prev = f_lam(0.0, pos, lam, a, b);
    for i in 1..=n {
        let th = 2.0 * std::f32::consts::PI * i as f32 / n as f32;
        let cur = f_lam(th, pos, lam, a, b);
        if prev * cur < 0.0 {
            // Bisect to refine the root.
            let mut lo = 2.0 * std::f32::consts::PI * (i - 1) as f32 / n as f32;
            let mut hi = th;
            for _ in 0..40 {
                let mid = 0.5 * (lo + hi);
                if f_lam(mid, pos, lam, a, b) * f_lam(lo, pos, lam, a, b) < 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            let th_root = 0.5 * (lo + hi);
            out.push(vec2(th_root.cos(), th_root.sin()));
        }
        prev = cur;
    }
    out
}

/// Λ(θ) − lam for a unit velocity at angle θ.
///
/// Must match `torus::confocal::caustic` exactly (the canonical conserved
/// quantity `Λ = b − L² + (a−b)vy²`, `L = x·vy − y·vx`) — an earlier version
/// used a differently-normalized `vx²/a + vy²/b − L²/(ab)` (equal to
/// `caustic()/a`, not `caustic()`), so it solved for `lam/a` instead of
/// `lam` and silently returned no roots whenever the true target `lam`
/// exceeded that scaled-down range.
fn f_lam(th: f32, pos: Vec2, lam: f32, a: f32, b: f32) -> f32 {
    let (vx, vy) = (th.cos(), th.sin());
    let l = pos.x * vy - pos.y * vx;
    b - l * l + (a - b) * vy * vy - lam
}

/// A length scale for how much of the billiard the trajectory is allowed to
/// reach at the current second integral `lam`.
///
/// The accessible region is the part of the domain outside the caustic
/// `Q_Λ = 0` (the trajectory never crosses the caustic).  As Λ changes, this
/// region swells or shrinks: near the boundary the caustic hugs the wall and
/// only a thin annulus is reachable (small torus); when the caustic collapses
/// inward the trajectory covers nearly the whole domain (big torus).  We
/// measure it by the mean distance of the caustic's in-domain points from the
/// domain centre, normalised so the torus radius tracks it.
pub fn accessible_region_scale(domain: &crate::domain::Domain, lam: f32) -> f32 {
    let cf = crate::torus::ConfocalParams::standard();
    // Degenerate caustics (Λ = b or Λ = a) collapse to a segment, where the
    // caustic sampler is ill-defined.  Instead of jumping to a fallback size
    // (which made the torus scale discontinuous when sweeping λ), evaluate at
    // a λ nudged just off the degeneracy — the caustic extent is continuous
    // there, so the returned scale varies smoothly through the critical level.
    const NUDGE: f32 = 2e-3;
    let lam = if (lam - cf.b).abs() < NUDGE {
        cf.b - NUDGE
    } else if (lam - cf.a).abs() < NUDGE {
        cf.a - NUDGE
    } else {
        lam
    };

    let quad = crate::quadratic::confocal(cf, lam);
    let pts = quad.sample_boundary(400);
    // The caustic's in-domain points: the accessible region is the annulus
    // between the caustic and the boundary, so the caustic's own extent tells
    // us how far the trajectory can reach inward.
    let mut sum = 0.0f32;
    let mut n = 0u32;
    for p in pts {
        if domain.contains(p) {
            sum += p.length();
            n += 1;
        }
    }
    if n == 0 {
        return crate::render::CachedDomain::new(domain).domain_extent;
    }
    let mean_r = sum / n as f32;
    // Normalise to the domain extent so the torus fills a comparable fraction
    // of the view; the caustic sits between the centre and the wall, so this
    // ranges from ~0 (tiny accessible region) to ~1 (whole domain).
    let extent = crate::render::CachedDomain::new(domain).domain_extent;
    (mean_r / extent).clamp(0.2, 1.0) * extent
}

// ----------------------------------------------------------------
// 3D drawing helpers
// ----------------------------------------------------------------

/// Standard torus surface (see §12 of the spec):
/// `((R + r cos θ₁) cos θ₂, (R + r cos θ₁) sin θ₂, r sin θ₁)`.
///
/// θ₂ is the azimuth about the z-axis (major / toroidal); θ₁ parametrizes the
/// tube cross-section (minor / poloidal).
fn embed_surface(theta1: f32, theta2: f32, r_major: f32, r_minor: f32) -> Vec3 {
    let (sin1, cos1) = theta1.sin_cos();
    let (sin2, cos2) = theta2.sin_cos();
    vec3(
        (r_major + r_minor * cos1) * cos2,
        (r_major + r_minor * cos1) * sin2,
        r_minor * sin1,
    )
}

/// Normal of the standard torus surface at `(θ₁, θ₂)` for shading.
fn surface_normal(theta1: f32, theta2: f32) -> Vec3 {
    let (sin1, cos1) = theta1.sin_cos();
    let (sin2, cos2) = theta2.sin_cos();
    vec3(cos1 * cos2, cos1 * sin2, sin1)
}

/// Embed a phase-space point onto a Liouville torus (using its precomputed
/// angles) and return both the 3D position and the outward surface normal.
///
/// `torus_index` was baked into the point at sampling time; each disconnected
/// region gets its own torus, offset in space by `scale.gap` so multiple tori
/// are distinct.  The scale carries the (possibly level-morphed, see
/// [`crate::torus_render::DomainScale::morph_to_level`]) radii and lobe gap,
/// so every drawer — tori, highlights, boundary preimage, degenerate curves —
/// embeds through the identical geometry.
pub fn torus_embed(
    pt: &PhasePoint,
    scale: &crate::torus_render::DomainScale,
    torus_index: u32,
) -> (Vec3, Vec3) {
    let pos = embed_surface(pt.theta1, pt.theta2, scale.r_major, scale.r_minor);
    let normal = surface_normal(pt.theta1, pt.theta2);
    // Offset each torus in space so multiple lobes are visually distinct.
    let offset = torus_index as f32 * scale.gap;
    (pos + vec3(offset, 0.0, 0.0), normal)
}

/// Assign each phase point to a torus ID, one per disconnected region.
///
/// The number of tori and the observable that splits them are given by the
/// [`TorusRegime`], so this always agrees with the app-side per-torus picker
/// (`one_per_torus`) and with the `torus_index` baked into sampled points.
pub fn torus_ids(points: &[[f32; 4]], regime: crate::confocal::TorusRegime) -> Vec<u32> {
    points
        .iter()
        .map(|q| regime.index_of(q[0], q[1], q[2], q[3]))
        .collect()
}

/// Draw short, bold red trajectories on the torus so one can see where each
/// torus's orbit actually lives on the 3D torus surface.
///
/// Each entry is the phase-space trace (already limited to a few bounces) of
/// one trajectory per disconnected torus, sampled from the *same* start point
/// as the trajectory drawn on the 2D billiard, so the two views correspond.
pub fn draw_torus_highlights(
    trajectories: &[Vec<PhasePoint>],
    cam: &OrbitCamera3,
    win_w: f32,
    win_h: f32,
) {
    draw_torus_highlights_scaled(
        trajectories,
        cam,
        win_w,
        win_h,
        &crate::torus_render::DomainScale::default(),
    );
}

/// [`draw_torus_highlights`] with an explicit torus scale (so the highlight
/// tracks the domain-derived torus size).
pub fn draw_torus_highlights_scaled(
    trajectories: &[Vec<PhasePoint>],
    cam: &OrbitCamera3,
    win_w: f32,
    win_h: f32,
    scale: &crate::torus_render::DomainScale,
) {
    // One representative per distinct torus index (trajectories never cross
    // between tori — each `torus_index` is a disconnected accessible region).
    let mut seen: Vec<u32> = Vec::new();
    for traj in trajectories {
        if traj.is_empty() {
            continue;
        }
        let idx = traj[0].torus_index;
        if seen.contains(&idx) {
            continue;
        }
        seen.push(idx);

        // Project each sample onto its torus; `torus_index` was baked in at
        // sampling time.
        let projected: Vec<Vec3> = traj
            .iter()
            .map(|pt| {
                let (p, _n) = torus_embed(pt, scale, pt.torus_index);
                cam.project(p, win_w, win_h)
            })
            .collect();

        // Bold red polyline through the consecutive samples.
        for win in projected.windows(2) {
            let a = win[0];
            let b = win[1];
            if a.z < -0.1 || b.z < -0.1 {
                continue;
            }
            draw_line(a.x, a.y, b.x, b.y, 4.0, RED);
        }

        // Bright marker dots so the short segment is easy to spot.
        for s in &projected {
            if s.z < -0.1 {
                continue;
            }
            draw_circle(s.x, s.y, 4.0, RED);
        }
    }
}

/// Draw the torus preimage π⁻¹ of the billiard boundary: every phase-space
/// point whose velocity is tangent to the caustic `Q_Λ = 0` and whose position
/// lies on the billiard boundary (the walls) or on the caustic itself (the
/// turning-point curves).
///
/// These are the curves on the torus that the *domain's* bounding curves map
/// to under the projection π, shown in a distinct color from the interior
/// point cloud so the size/shape of the billiard boundary is visible on the
/// abstract torus as it grows/shrinks with the domain.
pub fn boundary_preimage_points(
    domain: &crate::domain::Domain,
    lam: f32,
    bounds: (f32, Option<f32>),
) -> Vec<(PhasePoint, bool)> {
    let cf = crate::torus::ConfocalParams::standard();
    let mut cache = crate::torus::TorusCache::default();
    let mut out = Vec::new();

    // Map one physical point (position + one phase-space velocity, which is
    // always a caustic tangent) onto the torus.
    let mut map_velocity = |p: Vec2, vel: Vec2, is_caustic: bool| {
        let sample = crate::torus::PhaseSample::new(p.x, p.y, vel.x, vel.y);
        let mut params = crate::torus::TorusParams {
            confocal: &cf,
            lam_wall: bounds.0,
            beta: bounds.1,
            sep_eps: 1e-9,
            cache: &mut cache,
        };
        let (th1, th2, tidx) = crate::torus::to_torus(&sample, &mut params);
        let tidx = if tidx == u32::MAX { 0 } else { tidx };
        out.push((
            PhasePoint {
                x: p.x,
                y: p.y,
                theta: vel.y.atan2(vel.x) / std::f32::consts::PI,
                theta1: th1,
                theta2: th2,
                torus_index: tidx,
            },
            is_caustic,
        ));
    };

    // Elliptic border: on the elliptic side (λc < b) the outer ellipse is a
    // turning point of the λ₁-libration, where the ball glides tangent to the
    // ellipse.  This maps to a *full circle* on each of the two tori (upper arc
    // → torus 0, lower arc → torus 1).  The velocity is tangent to the outer
    // ellipse Q_λ_ell = 0, i.e. v ∥ ∇⊥Q_λ_ell (d₁ = 0), not a caustic tangent.
    let ell_wall = crate::confocal::ConfocalStructure::of_domain(domain)
        .map(|s| s.lambda_ell)
        .unwrap_or(0.0);
    if lam < cf.b - 1e-4 {
        let ell_quad = crate::quadratic::confocal(cf, ell_wall);
        // Sample the outer-ellipse arcs directly (denser and exact).
        let mut samples = Vec::new();
        for seg in &domain.segments {
            if let crate::domain::Segment::Quad { curve, a, b } = seg {
                if (curve.lambda - ell_wall).abs() < 1e-4 {
                    samples.extend(crate::domain::sample_wall_arc(curve, *a, *b, 40));
                }
            }
        }
        for p in samples {
            // Tangent to the outer ellipse (the λ₁-libration turning point).
            let g = ell_quad.grad(p);
            let tan = vec2(-g.y, g.x).normalize();
            map_velocity(p, tan, false);
            map_velocity(p, -tan, false);
        }
    }

    // Caustic arcs Q_Λ = 0 inside the domain (the turning-point curves).  At
    // a caustic point the velocity is exactly tangent to the caustic, given
    // analytically by rotating the gradient (the two ± signs are the two
    // phase-space sheets meeting at the turning point).
    if !((lam - cf.b).abs() < 0.03 || (lam - cf.a).abs() < 0.03) {
        let quad = crate::quadratic::confocal(cf, lam);
        for p in quad.sample_boundary(60) {
            if domain.contains(p) {
                let g = quad.grad(p);
                let tan = vec2(-g.y, g.x).normalize();
                map_velocity(p, tan, true);
                map_velocity(p, -tan, true);
            }
        }
    }

    out
}

/// Draw the boundary preimage points in distinct colours (cyan for the domain
/// wall, orange for the caustic), on top of the dimmed interior cloud.
pub fn draw_boundary_preimage(
    preimages: &[(PhasePoint, bool)],
    scale: &crate::torus_render::DomainScale,
    cam: &OrbitCamera3,
    win_w: f32,
    win_h: f32,
) {
    if preimages.is_empty() {
        return;
    }
    // Sort by camera depth so nearer boundary points draw on top.
    let mut pts: Vec<(Vec3, u32, bool)> = preimages
        .iter()
        .map(|(pt, is_ca)| {
            let (p, _n) = torus_embed(pt, scale, pt.torus_index);
            let s = cam.project(p, win_w, win_h);
            (s, pt.torus_index, *is_ca)
        })
        .collect();
    pts.sort_by(|a, b| {
        a.0.z
            .partial_cmp(&b.0.z)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    for (s, tidx, is_ca) in &pts {
        if s.z < -0.1 {
            continue;
        }
        // Keep boundary dots bright and opaque so they stand out against the
        // dimmed interior cloud.
        let depth = (-s.z).clamp(0.5, 5.0);
        let alpha = (0.65 + 0.35 * (1.0 - (depth - 0.5) / 4.5)).clamp(0.4, 1.0);
        // Separate the two torus lobes by alternating brightness.
        let dim = if *tidx % 2 == 0 { 1.0 } else { 0.6 };
        let (r, g, b) = if *is_ca {
            (255, 130, 20)
        } else {
            (70, 210, 255)
        };
        draw_circle(
            s.x,
            s.y,
            3.0,
            color_u8!(
                (r as f32 * dim) as u8,
                (g as f32 * dim) as u8,
                (b as f32 * dim) as u8,
                (alpha * 255.0) as u8
            ),
        );
    }
}

pub fn draw_phase_trajectories(
    trajectories: &[Vec<PhasePoint>],
    cam: &OrbitCamera3,
    win_w: f32,
    win_h: f32,
) {
    draw_phase_trajectories_scaled(
        trajectories,
        cam,
        win_w,
        win_h,
        &crate::torus_render::DomainScale::default(),
    );
}

/// [`draw_phase_trajectories`] with an explicit torus scale.
pub fn draw_phase_trajectories_scaled(
    trajectories: &[Vec<PhasePoint>],
    cam: &OrbitCamera3,
    win_w: f32,
    win_h: f32,
    scale: &crate::torus_render::DomainScale,
) {
    for traj in trajectories {
        if traj.len() < 2 {
            continue;
        }

        // Project all points onto the torus using the precomputed torus
        // index (baked in at sampling time).
        let projected: Vec<Vec3> = traj
            .iter()
            .map(|pt| {
                let (p, _n) = torus_embed(pt, scale, pt.torus_index);
                cam.project(p, win_w, win_h)
            })
            .collect();

        // Draw segments with depth-sorted hues
        let n = projected.len();
        for (i, chunk) in projected.windows(2).enumerate() {
            let a = chunk[0];
            let b = chunk[1];
            // Skip if behind camera
            if a.z < -0.1 || b.z < -0.1 {
                continue;
            }
            let t = i as f32 / n as f32;
            let hue = 30.0 + t * 200.0;
            let color = hsl_to_rgb(hue, 0.85, 0.55);
            // Depth-based alpha
            let depth = (-a.z.min(b.z)).clamp(0.5, 5.0);
            let alpha = (0.4 + 0.6 * (1.0 - (depth - 0.5) / 4.5)).min(1.0);
            let mut c = color;
            c.a = alpha;
            draw_line(a.x, a.y, b.x, b.y, 1.6, c);
        }
    }
}

/// Draw the closed 1D curves of a degenerate caustic level
/// (`Manifold::Curve1D`): one circle per collapsed torus, each an ordered loop
/// of torus-angle points.
///
/// Every point is embedded via [`torus_embed`] with the *same* (level-morphed)
/// scale and per-torus offset as the regular tori, so the circles are the
/// exact geometric limit of the collapsing tori:
///
/// - at λ_ell, λ = b and λ = a the morph takes `r_minor → 0`, so each circle
///   is the equator ring the thinning torus converges onto (θ₁ is always the
///   collapsing angle — see `torus::map`);
/// - at λ = b the morph also takes `gap → 2·r_major`, so the two rings touch
///   at one point — the figure-eight needs no extra positioning.
pub fn draw_curve1d(
    circles: &[Vec<PhasePoint>],
    scale: &crate::torus_render::DomainScale,
    cam: &OrbitCamera3,
    win_w: f32,
    win_h: f32,
) {
    const POINT_BUDGET: usize = 18_000;

    for circle in circles {
        if circle.len() < 2 {
            continue;
        }
        let step = decimation_step(circle.len(), POINT_BUDGET);

        // Project the loop points through the (morphed) torus embedding.
        let projected: Vec<Vec3> = circle
            .iter()
            .map(|pt| {
                let (p, _n) = torus_embed(pt, scale, pt.torus_index);
                cam.project(p, win_w, win_h)
            })
            .collect();

        // Closed polyline (last segment wraps back to the first point).
        let n = projected.len();
        for i in 0..n {
            let a = projected[i];
            let b = projected[(i + 1) % n];
            if a.z < -0.1 || b.z < -0.1 {
                continue;
            }
            let t = i as f32 / n as f32;
            let hue = 30.0 + t * 200.0;
            let color = crate::render::hsl_to_rgb(hue, 0.85, 0.55);
            let depth = (-a.z.min(b.z)).clamp(0.5, 5.0);
            let alpha = (0.4 + 0.6 * (1.0 - (depth - 0.5) / 4.5)).min(1.0);
            let mut c = color;
            c.a = alpha;
            draw_line(a.x, a.y, b.x, b.y, 2.0, c);
        }

        // Bright marker dots so the orbit is easy to spot.
        for (i, s) in projected.iter().enumerate() {
            if i % step != 0 || s.z < -0.1 {
                continue;
            }
            let depth = (-s.z).clamp(0.5, 5.0);
            let alpha = (0.5 + 0.3 * (1.0 - (depth - 0.5) / 4.5)).clamp(0.2, 0.9);
            let radius = 1.6 + 0.5 * (1.0 - (depth - 0.5) / 4.5);
            let (r, g, b) = (210, 120, 240); // violet: critical/singular
            draw_circle(s.x, s.y, radius, color_u8!(r, g, b, (alpha * 255.0) as u8));
        }
    }
}

/// Draw a dense set of phase-space points as small dots, filling out the
/// 2D Liouville torus surface.
///
/// Each trajectory contributes a 1D curve; the union of many trajectories at
/// the same Λ sweeps out the full 2D torus.  In action-angle coordinates the
/// torus is the flat product S¹ × S¹.  Here we embed it as a *standard donut*
/// in 3D: the two angle coordinates are the position angle φ₁ = atan2(y, x)
/// and the velocity angle φ₂ = θ = atan2(vy, vx).  Mapping (φ₁, φ₂) onto a
/// parametrized torus makes the invariant surface visibly a torus, instead of
/// the warped/flat-looking surface one gets from the raw (x, y, θ) embedding.
pub fn draw_phase_points(
    trajectories: &[Vec<PhasePoint>],
    cam: &OrbitCamera3,
    win_w: f32,
    win_h: f32,
) {
    draw_phase_points_scaled(
        trajectories,
        cam,
        win_w,
        win_h,
        &crate::torus_render::DomainScale::default(),
    );
}

/// [`draw_phase_points`] with an explicit torus scale.
pub fn draw_phase_points_scaled(
    trajectories: &[Vec<PhasePoint>],
    cam: &OrbitCamera3,
    win_w: f32,
    win_h: f32,
    scale: &crate::torus_render::DomainScale,
) {
    // Bound the per-raster point count so orbiting doesn't freeze the app.  A
    // large cloud is decimated with a uniform stride; the budget keeps the
    // profiled ~150K-draw_circle wall down to an interactive rate while
    // preserving spatial coverage.
    const POINT_BUDGET: usize = 18_000;
    let total: usize = trajectories.iter().map(|t| t.len()).sum();
    let step = decimation_step(total, POINT_BUDGET);

    for traj in trajectories {
        for (pt_idx, pt) in traj.iter().enumerate() {
            if pt_idx % step != 0 {
                continue;
            }
            let (p, normal) = torus_embed(pt, scale, pt.torus_index);
            let s = cam.project(p, win_w, win_h);
            if s.z < -0.1 {
                continue;
            }

            // Lambertian shading: brightness ∝ max(0, n · light).  This makes
            // the 3D curvature of the torus visible instead of a flat color.
            // The interior cloud is deliberately kept dim so the boundary
            // preimage (walls/caustic) drawn on top stays clearly visible.
            let light = vec3(0.4, 0.6, 0.7).normalize();
            let lambert = (normal.dot(light)).max(0.0);
            let shade = 0.16 + 0.38 * lambert;

            // Depth-based alpha and size (interior kept semi-transparent).
            let depth = (-s.z).clamp(0.5, 5.0);
            let alpha = (0.22 + 0.35 * (1.0 - (depth - 0.5) / 4.5)).clamp(0.05, 0.6);
            let radius = 1.4 + 0.6 * (1.0 - (depth - 0.5) / 4.5);

            // Warm base color, scaled by the lighting (dimmer than before).
            let r = (255.0 * shade) as u8;
            let g = (170.0 * shade) as u8;
            let b = (90.0 * shade) as u8;
            draw_circle(s.x, s.y, radius, color_u8!(r, g, b, (alpha * 255.0) as u8));
        }
    }
}

/// Draw the 3D axes and a bounding box to give spatial context.
pub fn draw_axes(cam: &OrbitCamera3, win_w: f32, win_h: f32, domain_extent: f32) {
    let scale = 1.0 / domain_extent.max(0.5);
    let axis_len = 1.3;

    // X axis (red)
    let x0 = cam.project(vec3(0.0, 0.0, 0.0), win_w, win_h);
    let x1 = cam.project(vec3(axis_len, 0.0, 0.0), win_w, win_h);
    if x0.z > -0.1 && x1.z > -0.1 {
        draw_line(x0.x, x0.y, x1.x, x1.y, 2.0, RED);
    }

    // Y axis (green)
    let y1 = cam.project(vec3(0.0, axis_len, 0.0), win_w, win_h);
    if x0.z > -0.1 && y1.z > -0.1 {
        draw_line(x0.x, x0.y, y1.x, y1.y, 2.0, GREEN);
    }

    // Z axis (blue) — velocity direction
    let z1 = cam.project(vec3(0.0, 0.0, axis_len), win_w, win_h);
    if x0.z > -0.1 && z1.z > -0.1 {
        draw_line(x0.x, x0.y, z1.x, z1.y, 2.0, BLUE);
    }

    // Axis labels
    draw_text("x", x1.x + 4.0, x1.y + 4.0, 14.0, RED);
    draw_text("y", y1.x + 4.0, y1.y + 4.0, 14.0, GREEN);
    draw_text("θ/π", z1.x + 4.0, z1.y + 4.0, 14.0, BLUE);

    // Draw a faint grid at θ = 0 (the billiard floor)
    let grid_n = 8;
    let grid_size = domain_extent * scale;
    for i in 0..=grid_n {
        let t = -grid_size + 2.0 * grid_size * i as f32 / grid_n as f32;
        // Line along x at fixed y, z=0
        let a = cam.project(vec3(t, -grid_size, 0.0), win_w, win_h);
        let b = cam.project(vec3(t, grid_size, 0.0), win_w, win_h);
        if a.z > -0.1 && b.z > -0.1 {
            draw_line(a.x, a.y, b.x, b.y, 0.5, color_u8!(80, 80, 120, 60));
        }
        // Line along y at fixed x, z=0
        let a = cam.project(vec3(-grid_size, t, 0.0), win_w, win_h);
        let b = cam.project(vec3(grid_size, t, 0.0), win_w, win_h);
        if a.z > -0.1 && b.z > -0.1 {
            draw_line(a.x, a.y, b.x, b.y, 0.5, color_u8!(80, 80, 120, 60));
        }
    }
}

fn hsl_to_rgb(hue: f32, s: f32, l: f32) -> Color {
    let h = hue / 360.0;
    let c = (1.0 - (2.0 * l - 1.0).abs()) * s;
    let x = c * (1.0 - ((h * 6.0) % 2.0 - 1.0).abs());
    let m = l - c / 2.0;
    let (r, g, b) = match (h * 6.0).floor() as i32 % 6 {
        0 => (c, x, 0.0),
        1 => (x, c, 0.0),
        2 => (0.0, c, x),
        3 => (0.0, x, c),
        4 => (x, 0.0, c),
        _ => (c, 0.0, x),
    };
    Color::new(r + m, g + m, b + m, 1.0)
}
