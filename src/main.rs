use billiards::{domain, phase3d, presets, quadratic, start_points_on_caustic, A, B};
use macroquad::prelude::*;
use std::cmp::Ordering;

// ----------------------------------------------------------------
// Camera
// ----------------------------------------------------------------
struct Camera {
    centre: Vec2,
    half_height: f32,
}

impl Camera {
    fn fit_domain(pts: &[Vec2]) -> Self {
        if pts.is_empty() {
            return Self {
                centre: vec2(0.0, 0.0),
                half_height: 2.0,
            };
        }
        let (mut min_x, mut max_x) = (f32::MAX, f32::MIN);
        let (mut min_y, mut max_y) = (f32::MAX, f32::MIN);
        for p in pts {
            if p.x < min_x {
                min_x = p.x;
            }
            if p.x > max_x {
                max_x = p.x;
            }
            if p.y < min_y {
                min_y = p.y;
            }
            if p.y > max_y {
                max_y = p.y;
            }
        }
        let cx = (min_x + max_x) / 2.0;
        let cy = (min_y + max_y) / 2.0;
        let half = ((max_x - min_x) / 2.0).max((max_y - min_y) / 2.0).max(0.5) * 1.2;
        Self {
            centre: vec2(cx, cy),
            half_height: half,
        }
    }

    #[inline(always)]
    fn world_to_screen(&self, p: Vec2, win_w: f32, win_h: f32) -> Vec2 {
        let aspect = win_w / win_h;
        let sx = ((p.x - self.centre.x) / (self.half_height * aspect) + 1.0) * 0.5 * win_w;
        let sy = (-(p.y - self.centre.y) / self.half_height + 1.0) * 0.5 * win_h;
        vec2(sx, sy)
    }
}

// ----------------------------------------------------------------
// Cached boundary
// ----------------------------------------------------------------
struct CachedDomain {
    boundary_pts: Vec<Vec2>,
    corners: Vec<Vec2>,
    camera: Camera,
    domain_extent: f32,
}

impl CachedDomain {
    fn new(domain: &domain::Domain) -> Self {
        let boundary_pts = domain.sample_boundary(30);
        let corners = domain.corners();
        let camera = Camera::fit_domain(&boundary_pts);

        // Compute domain extent for 3D view scaling
        let mut max_r = 0.0f32;
        for p in &boundary_pts {
            max_r = max_r.max(p.x.abs().max(p.y.abs()));
        }
        let domain_extent = max_r.max(0.5);

        Self {
            boundary_pts,
            corners,
            camera,
            domain_extent,
        }
    }
}

// ----------------------------------------------------------------
// Drawing
// ----------------------------------------------------------------
const BG: Color = color_u8!(15, 15, 35, 255);
const BOUNDARY: Color = color_u8!(180, 220, 255, 200);
const CORNER: Color = color_u8!(255, 200, 100, 200);

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

fn draw_boundary(cache: &CachedDomain, w: f32, h: f32) {
    let n = cache.boundary_pts.len();
    let mut sp: Vec<Vec2> = Vec::with_capacity(n);
    for p in &cache.boundary_pts {
        sp.push(cache.camera.world_to_screen(*p, w, h));
    }
    for i in 0..n {
        let j = (i + 1) % n;
        draw_line(sp[i].x, sp[i].y, sp[j].x, sp[j].y, 2.0, BOUNDARY);
    }
    for &c in &cache.corners {
        let s = cache.camera.world_to_screen(c, w, h);
        draw_circle(s.x, s.y, 4.0, CORNER);
    }
}

fn draw_trajectory(segs: &[(Vec2, Vec2)], cam: &Camera, w: f32, h: f32) {
    let n = segs.len();
    if n == 0 {
        return;
    }
    for (i, &(a, b)) in segs.iter().enumerate() {
        let t = i as f32 / n as f32;
        let hue = 30.0 + t * 200.0;
        let sa = cam.world_to_screen(a, w, h);
        let sb = cam.world_to_screen(b, w, h);
        draw_line(sa.x, sa.y, sb.x, sb.y, 1.8, hsl_to_rgb(hue, 0.85, 0.55));
    }
    for (i, &(_, b)) in segs.iter().enumerate().step_by(3) {
        let t = i as f32 / n as f32;
        let hue = 30.0 + t * 200.0;
        let sb = cam.world_to_screen(b, w, h);
        draw_circle(sb.x, sb.y, 2.5 + (1.0 - t) * 2.0, hsl_to_rgb(hue, 0.8, 0.7));
    }
}

fn draw_start_marker(pos: Vec2, cam: &Camera, w: f32, h: f32) {
    let s = cam.world_to_screen(pos, w, h);
    draw_circle_lines(s.x, s.y, 7.0, 2.5, YELLOW);
    draw_circle(s.x, s.y, 3.5, YELLOW);
}

/// Draw the degenerate caustic reference lines (focal segment at λ = B
/// and vertical segment at λ = A), plus the current caustic curve.
fn draw_caustic(a: f32, b: f32, lam: f32, domain: &domain::Domain, cam: &Camera, w: f32, h: f32) {
    // Draw degenerate caustic at λ = B: segment between foci (±c, 0)
    let c = (a - b).sqrt();
    let f1 = cam.world_to_screen(vec2(-c, 0.0), w, h);
    let f2 = cam.world_to_screen(vec2(c, 0.0), w, h);
    draw_line(f1.x, f1.y, f2.x, f2.y, 1.5, color_u8!(255, 200, 50, 80));
    draw_circle(f1.x, f1.y, 2.5, color_u8!(255, 200, 50, 100));
    draw_circle(f2.x, f2.y, 2.5, color_u8!(255, 200, 50, 100));

    // Draw degenerate caustic at λ = A: vertical segment x = 0, y ∈ [−√b, √b]
    let y = b.sqrt();
    let v1 = cam.world_to_screen(vec2(0.0, -y), w, h);
    let v2 = cam.world_to_screen(vec2(0.0, y), w, h);
    draw_line(v1.x, v1.y, v2.x, v2.y, 1.5, color_u8!(255, 200, 50, 80));
    draw_circle(v1.x, v1.y, 2.5, color_u8!(255, 200, 50, 100));
    draw_circle(v2.x, v2.y, 2.5, color_u8!(255, 200, 50, 100));

    // Draw the current caustic curve (skip if degenerate)
    if (lam - b).abs() < 0.03 || (lam - a).abs() < 0.03 {
        return;
    }
    let quad = quadratic::confocal(a, b, lam);
    let pts = quad.sample_boundary(80);
    for &p in &pts {
        if domain.contains(p) {
            let s = cam.world_to_screen(p, w, h);
            draw_circle(s.x, s.y, 1.8, color_u8!(255, 255, 100, 160));
        }
    }
}

/// Extract confocal boundary parameters from a domain.
fn get_boundary_lambdas(domain: &domain::Domain) -> (f32, f32) {
    let mut ell = -A + 1.0;
    let mut hyp = B + 1.0;
    for seg in &domain.segments {
        if let domain::Segment::Quad { curve, .. } = seg {
            if curve.lambda < B {
                ell = curve.lambda;
            }
            if curve.lambda > B {
                hyp = curve.lambda;
            }
        }
    }
    (ell, hyp)
}

/// Clamp Λ to the valid range, jumping over the gap between ellipse and
/// hyperbola ranges. Uses `prev` to determine direction when crossing.
///
/// Ellipse caustic: λ_ell < Λ < B  (inside boundary ellipse)
/// Hyperbola caustic: λ_hyp < Λ < A (inside boundary hyperbola)
fn clamp_lambda(lam: f32, prev: f32, domain: &domain::Domain) -> f32 {
    let (lambda_ell, lambda_hyp) = get_boundary_lambdas(domain);
    let e = 0.02;
    let ell_max = B - e;
    let hyp_min = lambda_hyp + e;
    let range_min = lambda_ell + e;
    let range_max = A - e;

    let clamped = lam.clamp(range_min, range_max);

    // If clamped lands in the gap between ellipse and hyperbola ranges,
    // jump to the appropriate side based on travel direction.
    if clamped > ell_max && clamped < hyp_min {
        match lam.partial_cmp(&prev).unwrap_or(Ordering::Equal) {
            Ordering::Greater => hyp_min,
            _ => ell_max,
        }
    } else {
        clamped
    }
}

/// Slider range info for a given domain.
struct SliderRange {
    ell_min: f32,
    ell_max: f32,
    hyp_min: f32,
    #[allow(dead_code)]
    hyp_max: f32,
    total_range: f32,
}

impl SliderRange {
    fn new(domain: &domain::Domain) -> Self {
        let (lambda_ell, lambda_hyp) = get_boundary_lambdas(domain);
        let e = 0.02;
        let ell_min = lambda_ell + e;
        let ell_max = B - e;
        let hyp_min = lambda_hyp + e;
        let _hyp_max = A - e;
        let ell_range = (ell_max - ell_min).max(0.0);
        let hyp_range = (_hyp_max - hyp_min).max(0.0);
        let total_range = ell_range + hyp_range;
        Self {
            ell_min,
            ell_max,
            hyp_min,
            hyp_max: _hyp_max,
            total_range,
        }
    }

    /// Map slider fraction t ∈ [0, 1] to lambda, skipping the gap.
    fn to_lambda(&self, t: f32) -> f32 {
        if self.total_range <= 0.0 {
            return 0.0;
        }
        let ell_range = self.ell_max - self.ell_min;
        let t = t.clamp(0.0, 1.0);
        let pos = t * self.total_range;
        if pos <= ell_range {
            self.ell_min + pos
        } else {
            self.hyp_min + (pos - ell_range)
        }
    }

    /// Map lambda to slider fraction t ∈ [0, 1].
    fn from_lambda(&self, lam: f32) -> f32 {
        if self.total_range <= 0.0 {
            return 0.0;
        }
        let ell_range = self.ell_max - self.ell_min;
        if lam <= self.ell_max {
            (lam - self.ell_min) / self.total_range
        } else if lam >= self.hyp_min {
            (ell_range + (lam - self.hyp_min)) / self.total_range
        } else {
            // In the gap — map to the nearest edge (ellipse end)
            ell_range / self.total_range
        }
    }
}

// ----------------------------------------------------------------
// Slider widget
// ----------------------------------------------------------------
struct Slider {
    x: f32,
    y: f32,
    width: f32,
    dragging: bool,
}

impl Slider {
    fn new() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            width: 0.0,
            dragging: false,
        }
    }

    /// Draw the slider and handle mouse interaction.
    /// Returns the new lambda value if the slider moved.
    fn update(&mut self, lambda: f32, range: &SliderRange) -> Option<f32> {
        let track_h = 4.0;
        let thumb_r = 8.0;
        let cy = self.y + track_h / 2.0;

        // Draw track (background)
        draw_rectangle(
            self.x,
            self.y,
            self.width,
            track_h,
            color_u8!(60, 60, 100, 180),
        );

        // Draw filled portion
        let t = range.from_lambda(lambda);
        let fill_w = t * self.width;
        if fill_w > 0.0 {
            draw_rectangle(
                self.x,
                self.y,
                fill_w,
                track_h,
                color_u8!(130, 130, 200, 220),
            );
        }

        // Draw thumb
        let thumb_x = self.x + t * self.width;
        draw_circle(thumb_x, cy, thumb_r, color_u8!(220, 220, 255, 255));
        draw_circle_lines(thumb_x, cy, thumb_r, 1.5, color_u8!(100, 100, 160, 200));

        // Draw gap indicator (a small break in the track)
        if range.ell_max < range.hyp_min {
            let gap_t = range.from_lambda(range.ell_max);
            let gap_x = self.x + gap_t * self.width;
            draw_line(gap_x, self.y - 2.0, gap_x, self.y + track_h + 2.0, 2.0, BG);
        }

        // Handle mouse interaction
        let mx = mouse_position().0;
        let my = mouse_position().1;

        if is_mouse_button_pressed(MouseButton::Left) {
            let dist = ((mx - thumb_x).powi(2) + (my - cy).powi(2)).sqrt();
            if dist <= thumb_r + 4.0
                || (mx >= self.x && mx <= self.x + self.width && (my - cy).abs() <= 12.0)
            {
                self.dragging = true;
            }
        }

        if self.dragging {
            if is_mouse_button_down(MouseButton::Left) {
                let t = ((mx - self.x) / self.width).clamp(0.0, 1.0);
                return Some(range.to_lambda(t));
            } else {
                self.dragging = false;
            }
        }

        None
    }
}

// ----------------------------------------------------------------
// Entry point
// ----------------------------------------------------------------
#[macroquad::main("Mathematical Billiards")]
async fn main() {
    let configs = presets::all_presets();
    let caches: Vec<CachedDomain> = configs
        .iter()
        .map(|p| CachedDomain::new(&p.domain))
        .collect();

    let mut idx = 0usize;
    let mut lambda: f32 = 0.2;
    let mut slider = Slider::new();
    let mut show_3d = false;
    let mut cam3d = phase3d::OrbitCamera3::new();
    let mut phase_trajectories: Vec<Vec<phase3d::PhasePoint>> = Vec::new();

    let trace_all = |domain: &domain::Domain, lam: f32| -> Vec<Vec<(Vec2, Vec2)>> {
        let starts = start_points_on_caustic(A, B, lam, domain);
        starts
            .iter()
            .map(|(p, v)| domain.trace(*p, *v, 300))
            .collect()
    };

    let mut trajectories = trace_all(&configs[idx].domain, lambda);
    let mut start_points = start_points_on_caustic(A, B, lambda, &configs[idx].domain);

    loop {
        let (w, h) = (screen_width(), screen_height());
        let cache = &caches[idx];
        let cam = &cache.camera;

        // Position the slider
        slider.x = w * 0.15;
        slider.y = h - 50.0;
        slider.width = w * 0.7;

        // Build slider range for the current domain
        let slider_range = SliderRange::new(&configs[idx].domain);

        // Toggle 3D phase space view
        if is_key_pressed(KeyCode::P) {
            show_3d = !show_3d;
            if show_3d {
                // Build phase space trajectories
                let starts = start_points_on_caustic(A, B, lambda, &configs[idx].domain);
                phase_trajectories =
                    phase3d::sample_all_trajectories_phase(&configs[idx].domain, &starts, 300);
            }
        }

        // Keyboard input
        if is_key_pressed(KeyCode::Tab) || is_key_pressed(KeyCode::Space) {
            idx = (idx + 1) % configs.len();
            lambda = 0.2;
            trajectories = trace_all(&configs[idx].domain, lambda);
            start_points = start_points_on_caustic(A, B, lambda, &configs[idx].domain);
            if show_3d {
                phase_trajectories = phase3d::sample_all_trajectories_phase(
                    &configs[idx].domain,
                    &start_points,
                    300,
                );
            }
        }

        if is_key_pressed(KeyCode::Up) || is_key_pressed(KeyCode::Right) {
            let next = clamp_lambda(lambda + 0.05, lambda, &configs[idx].domain);
            lambda = next;
            trajectories = trace_all(&configs[idx].domain, lambda);
            start_points = start_points_on_caustic(A, B, lambda, &configs[idx].domain);
            if show_3d {
                phase_trajectories = phase3d::sample_all_trajectories_phase(
                    &configs[idx].domain,
                    &start_points,
                    300,
                );
            }
        }
        if is_key_pressed(KeyCode::Down) || is_key_pressed(KeyCode::Left) {
            let next = clamp_lambda(lambda - 0.05, lambda, &configs[idx].domain);
            lambda = next;
            trajectories = trace_all(&configs[idx].domain, lambda);
            start_points = start_points_on_caustic(A, B, lambda, &configs[idx].domain);
            if show_3d {
                phase_trajectories = phase3d::sample_all_trajectories_phase(
                    &configs[idx].domain,
                    &start_points,
                    300,
                );
            }
        }

        // Slider input
        if let Some(new_lam) = slider.update(lambda, &slider_range) {
            if (new_lam - lambda).abs() > 0.0001 {
                lambda = new_lam;
                trajectories = trace_all(&configs[idx].domain, lambda);
                start_points = start_points_on_caustic(A, B, lambda, &configs[idx].domain);
                if show_3d {
                    phase_trajectories = phase3d::sample_all_trajectories_phase(
                        &configs[idx].domain,
                        &start_points,
                        300,
                    );
                }
            }
        }

        clear_background(BG);

        if show_3d {
            // ---- 3D phase space view ----
            cam3d.handle_input();

            // Draw axes and grid
            phase3d::draw_axes(&cam3d, w, h, cache.domain_extent);

            // Draw phase trajectories
            phase3d::draw_phase_trajectories(
                &phase_trajectories,
                &cam3d,
                w,
                h,
                cache.domain_extent,
            );

            // HUD
            let caustic_type = if lambda < B { "ellipse" } else { "hyperbola" };
            let info =
                format!(
                "{}  |  Λ = {:.3} ({})  |  {} trajs  |  [P] 2D  |  [Tab] next  |  right-drag orbit",
                configs[idx].label, lambda, caustic_type, start_points.len(),
            );
            draw_text(&info, 12.0, 28.0, 18.0, color_u8!(200, 200, 220, 220));

            draw_text(
                "3D phase space: (x, y, θ/π)  —  Liouville torus at fixed Λ",
                12.0,
                h - 12.0,
                14.0,
                color_u8!(150, 150, 170, 140),
            );
        } else {
            // ---- 2D billiard view ----
            draw_boundary(cache, w, h);
            draw_caustic(A, B, lambda, &configs[idx].domain, cam, w, h);

            for traj in &trajectories {
                draw_trajectory(traj, cam, w, h);
            }
            for &(p, _) in &start_points {
                draw_start_marker(p, cam, w, h);
            }

            let caustic_type = if lambda < B { "ellipse" } else { "hyperbola" };
            let total: usize = trajectories.iter().map(|t| t.len()).sum();
            let info = format!(
                "{}  |  Λ = {:.3} ({})  |  {} trajs, {} bounces  |  ↑↓ Λ  |  [P] 3D  |  Tab next",
                configs[idx].label,
                lambda,
                caustic_type,
                start_points.len(),
                total,
            );
            draw_text(&info, 12.0, 28.0, 18.0, color_u8!(200, 200, 220, 220));

            draw_text(
                "H = ½|v|²  |  Λ = vx²/a + vy²/b − (x·vy − y·vx)²/(ab)  |  v ⟂ ∇Q_Λ",
                12.0,
                h - 12.0,
                14.0,
                color_u8!(150, 150, 170, 140),
            );
        }

        next_frame().await;
    }
}
