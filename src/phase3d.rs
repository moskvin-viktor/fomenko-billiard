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
            distance: 5.0,
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

// ----------------------------------------------------------------
// 3D drawing helpers
// ----------------------------------------------------------------

/// Draw a set of phase space trajectories as 3D curves.
pub fn draw_phase_trajectories(
    trajectories: &[Vec<PhasePoint>],
    cam: &OrbitCamera3,
    win_w: f32,
    win_h: f32,
    domain_extent: f32,
) {
    // Scale factor so the billiard fits in the [-1, 1] box in x/y
    // and θ is in [-1, 1] as well.
    let scale = 1.0 / domain_extent.max(0.5);

    for traj in trajectories {
        if traj.len() < 2 {
            continue;
        }

        // Project all points
        let projected: Vec<Vec3> = traj
            .iter()
            .map(|pt| {
                let p = vec3(pt.x * scale, pt.y * scale, pt.theta);
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
