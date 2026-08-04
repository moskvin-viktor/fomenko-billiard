use crate::{A, B};
use macroquad::prelude::*;

/// A 3D phase space point (x, y, θ/π) where θ = atan2(vy, vx)
/// is the velocity direction angle, normalized to [-1, 1].
#[derive(Clone, Copy, Debug)]
pub struct PhasePoint {
    pub x: f32,
    pub y: f32,
    pub theta: f32, // θ/π ∈ [-1, 1]
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
}

// ----------------------------------------------------------------
// Phase space sampling
// ----------------------------------------------------------------

/// Sample phase space points from a billiard trajectory.
/// For each bounce point, records (x, y, θ) where θ = atan2(vy, vx)
/// is the velocity direction right after reflection.
pub fn sample_trajectory_phase(
    domain: &crate::domain::Domain,
    p0: Vec2,
    v0: Vec2,
    max_steps: usize,
) -> Vec<PhasePoint> {
    let mut pts = Vec::with_capacity(max_steps);
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
        // Record the bounce point: position is hit, velocity is v (just before reflection)
        let theta = v.y.atan2(v.x) / std::f32::consts::PI;
        pts.push(PhasePoint {
            x: hit.x,
            y: hit.y,
            theta,
        });
        // Reflect and continue
        v = domain.reflect(hit, v, idx);
        p = hit + 1e-4 * v.normalize();
    }
    pts
}

/// Sample phase space points from multiple trajectories, one per
/// start point on the caustic.
pub fn sample_all_trajectories_phase(
    domain: &crate::domain::Domain,
    starts: &[(Vec2, Vec2)],
    max_steps: usize,
) -> Vec<Vec<PhasePoint>> {
    starts
        .iter()
        .map(|&(p, v)| sample_trajectory_phase(domain, p, v, max_steps))
        .collect()
}

/// Sample a single trajectory densely, recording interior points along each
/// segment (not just bounce points).
///
/// This is what makes the torus surface look filled in the 3D view: sampling
/// `per_seg` points along every segment of many trajectories sweeps out the
/// whole 2D surface instead of leaving sparse dots at the bounces.
pub fn sample_trajectory_phase_dense(
    domain: &crate::domain::Domain,
    p0: Vec2,
    v0: Vec2,
    max_steps: usize,
    per_seg: usize,
) -> Vec<PhasePoint> {
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
        for k in 0..per_seg {
            let f = k as f32 / per_seg as f32;
            let q = p + (hit - p) * f;
            let theta = v.y.atan2(v.x) / std::f32::consts::PI;
            pts.push(PhasePoint {
                x: q.x,
                y: q.y,
                theta,
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
fn f_lam(th: f32, pos: Vec2, lam: f32, a: f32, b: f32) -> f32 {
    let (vx, vy) = (th.cos(), th.sin());
    let x = pos.x;
    let y = pos.y;
    vx * vx / a + vy * vy / b - (x * vy - y * vx).powi(2) / (a * b) - lam
}

// ----------------------------------------------------------------
// 3D drawing helpers
// ----------------------------------------------------------------

/// Map a phase-space point onto a standard torus embedding and return both the
/// 3D position and the outward surface normal.
///
/// The two angle coordinates are the **elliptic coordinates** (μ, ν) of the
/// point: every (x, y) lies on one confocal ellipse (μ < b) and one confocal
/// hyperbola (ν > b).  These are the separable coordinates of the elliptic
/// billiard, so at fixed Λ the pair (μ, ν) fills a genuine 2D rectangle — the
/// correct torus parametrization.  (Using (atan2(y,x), θ) instead collapses to
/// a 1D curve at fixed Λ, which is why the torus looked flat/unrecognizable.)
///
/// The outward normal of the standard torus at (φ₁, φ₂) is:
///   n = (cos φ₂ · cos φ₁, cos φ₂ · sin φ₁, sin φ₂)
/// which lets us shade the surface so its 3D curvature is visible.
fn torus_embed(pt: &PhasePoint, r_major: f32, r_minor: f32, is_hyperbola: bool) -> (Vec3, Vec3) {
    // Elliptic coordinates (μ, ν) from the confocal family.
    let (mu, nu) = elliptic_coords(pt.x, pt.y);
    // Map μ ∈ (−∞, b] onto the first torus angle.
    let phi1 = normalize_angle(mu, -2.0, B);
    // The hyperbola coordinate ν ∈ [b, a] only encodes the magnitude.  For a
    // hyperbola caustic (Λ > B) the two disconnected lobes (left/right
    // branches) have the SAME ν, so we use the sign of x to place them on
    // opposite halves of the torus.  For an ellipse caustic (Λ < B) the torus
    // is a single connected surface, so no sign is needed.
    let phi2 = if is_hyperbola && pt.x >= 0.0 {
        normalize_angle(nu, B, A) // right branch → one half
    } else if is_hyperbola {
        -normalize_angle(nu, B, A) // left branch → the other half
    } else {
        normalize_angle(nu, B, A)
    };
    let cos2 = phi2.cos();
    let (sin1, cos1) = phi1.sin_cos();
    let pos = vec3(
        (r_major + r_minor * cos2) * cos1,
        (r_major + r_minor * cos2) * sin1,
        r_minor * phi2.sin(),
    );
    let normal = vec3(cos2 * cos1, cos2 * sin1, phi2.sin());
    (pos, normal)
}

/// Solve the confocal-family equation for the two elliptic coordinates (μ, ν)
/// of a point (x, y):  (b−λ)x² + (a−λ)y² = (a−λ)(b−λ).
fn elliptic_coords(x: f32, y: f32) -> (f32, f32) {
    let a = A;
    let b = B;
    let s = x * x + y * y;
    let t = b * x * x + a * y * y - a * b;
    let p = a + b - s;
    let disc = p * p + 4.0 * t;
    let sd = disc.max(0.0).sqrt();
    let lam1 = 0.5 * (p - sd);
    let lam2 = 0.5 * (p + sd);
    (lam1, lam2) // μ = lam1 (ellipse), ν = lam2 (hyperbola)
}

/// Linearly map a value in [lo, hi] onto the angle range [−π, π].
fn normalize_angle(v: f32, lo: f32, hi: f32) -> f32 {
    let t = ((v - lo) / (hi - lo)).clamp(0.0, 1.0);
    -std::f32::consts::PI + 2.0 * std::f32::consts::PI * t
}

/// Draw a set of phase space trajectories as 3D curves on the torus.
pub fn draw_phase_trajectories(
    trajectories: &[Vec<PhasePoint>],
    cam: &OrbitCamera3,
    win_w: f32,
    win_h: f32,
    is_hyperbola: bool,
) {
    let r_major = 1.6;
    let r_minor = 0.6;

    for traj in trajectories {
        if traj.len() < 2 {
            continue;
        }

        // Project all points onto the torus
        let projected: Vec<Vec3> = traj
            .iter()
            .map(|pt| {
                let (p, _n) = torus_embed(pt, r_major, r_minor, is_hyperbola);
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
    is_hyperbola: bool,
) {
    // Torus radii — fixed size so it fills the view.
    let r_major = 1.6;
    let r_minor = 0.6;

    for traj in trajectories {
        for pt in traj {
            let (p, normal) = torus_embed(pt, r_major, r_minor, is_hyperbola);
            let s = cam.project(p, win_w, win_h);
            if s.z < -0.1 {
                continue;
            }

            // Lambertian shading: brightness ∝ max(0, n · light).  This makes
            // the 3D curvature of the torus visible instead of a flat color.
            let light = vec3(0.4, 0.6, 0.7).normalize();
            let lambert = (normal.dot(light)).max(0.0);
            let shade = 0.25 + 0.75 * lambert;

            // Depth-based alpha and size
            let depth = (-s.z).clamp(0.5, 5.0);
            let alpha = (0.3 + 0.6 * (1.0 - (depth - 0.5) / 4.5)).min(1.0);
            let radius = 1.8 + 0.8 * (1.0 - (depth - 0.5) / 4.5);

            // Warm base color, scaled by the lighting.
            let r = (255.0 * shade) as u8;
            let g = (200.0 * shade) as u8;
            let b = (110.0 * shade) as u8;
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
