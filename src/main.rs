use billiards::{domain, phase3d, presets, quadratic, A, B};
use macroquad::prelude::*;
use std::cmp::Ordering;

// ----------------------------------------------------------------
// Camera (2D)
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
// Drawing helpers
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

/// Draw caustic curves and degenerate reference lines (for confocal domains).
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

/// Draw a starting direction arrow for polyline domains.
fn draw_start_arrow(center: Vec2, angle: f32, cam: &Camera, w: f32, h: f32) {
    let dir = vec2(angle.cos(), angle.sin());
    let tip = center + dir * 0.3;
    let s1 = cam.world_to_screen(center, w, h);
    let s2 = cam.world_to_screen(tip, w, h);
    draw_line(s1.x, s1.y, s2.x, s2.y, 2.5, YELLOW);
    draw_circle(s1.x, s1.y, 3.5, YELLOW);
}

// ----------------------------------------------------------------
// Second-integral slider range
// ----------------------------------------------------------------

/// Describes the valid range for the second-integral slider.
enum SecondIntegralRange {
    /// Confocal caustic parameter Λ.
    Confocal {
        ell_min: f32,
        ell_max: f32,
        hyp_min: f32,
        #[allow(dead_code)]
        hyp_max: f32,
        total_range: f32,
    },
    /// Polyline angle θ/π ∈ [-1, 1].
    Angle,
}

impl SecondIntegralRange {
    fn for_domain(domain: &domain::Domain, is_confocal: bool) -> Self {
        if !is_confocal {
            return Self::Angle;
        }
        let e = 0.02;
        let (lambda_ell, lambda_hyp) = get_ell_hyp(domain);
        let ell_min = lambda_ell + e;
        let ell_max = B - e;
        let hyp_min = lambda_hyp + e;
        let hyp_max = A - e;
        let ell_range = (ell_max - ell_min).max(0.0);
        let hyp_range = (hyp_max - hyp_min).max(0.0);
        let total_range = ell_range + hyp_range;
        Self::Confocal {
            ell_min,
            ell_max,
            hyp_min,
            hyp_max,
            total_range,
        }
    }

    /// Map slider fraction t ∈ [0, 1] to the second integral value.
    fn to_value(&self, t: f32) -> f32 {
        match self {
            Self::Angle => -1.0 + 2.0 * t.clamp(0.0, 1.0),
            Self::Confocal {
                ell_min,
                ell_max,
                hyp_min,
                total_range,
                ..
            } => {
                if *total_range <= 0.0 {
                    return 0.0;
                }
                let ell_range = ell_max - ell_min;
                let t = t.clamp(0.0, 1.0);
                let pos = t * total_range;
                if pos <= ell_range {
                    ell_min + pos
                } else {
                    hyp_min + (pos - ell_range)
                }
            }
        }
    }

    /// Map a second integral value to slider fraction t ∈ [0, 1].
    fn from_value(&self, lam: f32) -> f32 {
        match self {
            Self::Angle => (lam.clamp(-1.0, 1.0) + 1.0) / 2.0,
            Self::Confocal {
                ell_min,
                ell_max,
                hyp_min,
                total_range,
                ..
            } => {
                if *total_range <= 0.0 {
                    return 0.0;
                }
                let ell_range = ell_max - ell_min;
                if lam <= *ell_max {
                    (lam - ell_min) / total_range
                } else if lam >= *hyp_min {
                    (ell_range + (lam - hyp_min)) / total_range
                } else {
                    ell_range / total_range
                }
            }
        }
    }
}

/// Extract (λ_ell, λ_hyp) boundary lambdas from a confocal domain.
fn get_ell_hyp(domain: &domain::Domain) -> (f32, f32) {
    let mut ell = -A + 1.0;
    let mut hyp = B + 1.0;
    for seg in &domain.segments {
        if let domain::Segment::Quad { curve, .. } = seg {
            if curve.lambda < B {
                ell = ell.min(curve.lambda);
            }
            if curve.lambda > B {
                hyp = hyp.max(curve.lambda);
            }
        }
    }
    (ell, hyp)
}

/// Clamp Λ to the valid range, jumping over the gap.
fn clamp_lambda(lam: f32, prev: f32, domain: &domain::Domain) -> f32 {
    let (lambda_ell, lambda_hyp) = get_ell_hyp(domain);
    let e = 0.02;
    let ell_max = B - e;
    let hyp_min = lambda_hyp + e;
    let range_min = lambda_ell + e;
    let range_max = A - e;

    let clamped = lam.clamp(range_min, range_max);
    if clamped > ell_max && clamped < hyp_min {
        match lam.partial_cmp(&prev).unwrap_or(Ordering::Equal) {
            Ordering::Greater => hyp_min,
            _ => ell_max,
        }
    } else {
        clamped
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

    fn update(&mut self, value: f32, range: &SecondIntegralRange, label: &str) -> Option<f32> {
        let track_h = 4.0;
        let thumb_r = 8.0;
        let cy = self.y + track_h / 2.0;

        // Track background
        draw_rectangle(
            self.x,
            self.y,
            self.width,
            track_h,
            color_u8!(60, 60, 100, 180),
        );

        // Filled portion
        let t = range.from_value(value);
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

        // Thumb
        let thumb_x = self.x + t * self.width;
        draw_circle(thumb_x, cy, thumb_r, color_u8!(220, 220, 255, 255));
        draw_circle_lines(thumb_x, cy, thumb_r, 1.5, color_u8!(100, 100, 160, 200));

        // Gap indicator (confocal-only)
        if let SecondIntegralRange::Confocal {
            ell_max, hyp_min, ..
        } = range
        {
            if ell_max < hyp_min {
                let gap_t = range.from_value(*ell_max);
                let gap_x = self.x + gap_t * self.width;
                draw_line(gap_x, self.y - 2.0, gap_x, self.y + track_h + 2.0, 2.0, BG);
            }
        }

        // Label and value text
        let value_text = match range {
            SecondIntegralRange::Angle => format!("{:.3}", value),
            SecondIntegralRange::Confocal { .. } => {
                let caustic_type = if value < B { "ellipse" } else { "hyperbola" };
                format!("{:.3} ({})", value, caustic_type)
            }
        };
        let info = format!("{} = {}", label, value_text);
        draw_text(
            &info,
            self.x - 60.0,
            self.y + track_h + 22.0,
            14.0,
            color_u8!(180, 180, 210, 200),
        );

        // Mouse interaction
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
                return Some(range.to_value(t));
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
    let mut second_int: f32 = 0.2; // Λ for confocal, θ/π for polyline
    let mut slider = Slider::new();
    let mut show_3d = false;
    let mut animate = false;
    let mut anim_dir: f32 = 1.0;
    let mut cam3d = phase3d::OrbitCamera3::new();
    let mut phase_trajectories: Vec<Vec<phase3d::PhasePoint>> = Vec::new();

    let build_trajectories = |domain: &domain::Domain, lam: f32, preset: &presets::Preset| {
        let starts =
            billiards::get_start_points(A, B, lam, domain, preset.is_confocal, preset.start_center);
        let trajs: Vec<Vec<(Vec2, Vec2)>> = starts
            .iter()
            .map(|(p, v)| domain.trace(*p, *v, 300))
            .collect();
        (starts, trajs)
    };

    // Build dense phase-space trajectories for the 3D torus view.
    // For confocal billiards we densely sample the caustic so the union of
    // trajectories fills the 2D Liouville torus.  For polyline billiards we
    // just use the single start point.
    let build_phase = |domain: &domain::Domain,
                       lam: f32,
                       preset: &presets::Preset|
     -> Vec<Vec<phase3d::PhasePoint>> {
        let starts = if preset.is_confocal {
            billiards::dense_caustic_starts(A, B, lam, domain, 24)
        } else {
            billiards::get_start_points(A, B, lam, domain, false, preset.start_center)
        };
        let bounds = billiards::torus_bounds(domain);
        // Sample interior points along each segment so the torus surface is
        // densely filled (not just sparse bounce dots).  The torus angles are
        // computed in here, once per trajectory, sharing a TorusCache.
        starts
            .iter()
            .map(|&(p, v)| phase3d::sample_trajectory_phase_dense(domain, p, v, 200, 8, bounds))
            .collect()
    };

    let (mut start_points, mut trajectories) =
        build_trajectories(&configs[idx].domain, second_int, &configs[idx]);

    loop {
        let (w, h) = (screen_width(), screen_height());
        let cache = &caches[idx];
        let cam = &cache.camera;
        let preset = &configs[idx];

        // Position slider
        slider.x = w * 0.15;
        slider.y = h - 50.0;
        slider.width = w * 0.7;

        let slider_range = SecondIntegralRange::for_domain(&preset.domain, preset.is_confocal);

        // Toggle animation
        if is_key_pressed(KeyCode::A) {
            animate = !animate;
        }

        // Toggle 3D
        if is_key_pressed(KeyCode::P) {
            show_3d = !show_3d;
            if show_3d {
                phase_trajectories = build_phase(&preset.domain, second_int, preset);
            }
        }

        // Tab / Space: next preset
        if is_key_pressed(KeyCode::Tab) || is_key_pressed(KeyCode::Space) {
            idx = (idx + 1) % configs.len();
            second_int = 0.2;
            let (s, t) = build_trajectories(&configs[idx].domain, second_int, &configs[idx]);
            start_points = s;
            trajectories = t;
            if show_3d {
                phase_trajectories = build_phase(&configs[idx].domain, second_int, &configs[idx]);
            }
        }

        // Keyboard: adjust second integral
        if configs[idx].is_confocal {
            if is_key_pressed(KeyCode::Up) || is_key_pressed(KeyCode::Right) {
                let next = clamp_lambda(second_int + 0.05, second_int, &configs[idx].domain);
                second_int = next;
                let (s, t) = build_trajectories(&configs[idx].domain, second_int, &configs[idx]);
                start_points = s;
                trajectories = t;
                if show_3d {
                    phase_trajectories =
                        build_phase(&configs[idx].domain, second_int, &configs[idx]);
                }
            }
            if is_key_pressed(KeyCode::Down) || is_key_pressed(KeyCode::Left) {
                let next = clamp_lambda(second_int - 0.05, second_int, &configs[idx].domain);
                second_int = next;
                let (s, t) = build_trajectories(&configs[idx].domain, second_int, &configs[idx]);
                start_points = s;
                trajectories = t;
                if show_3d {
                    phase_trajectories =
                        build_phase(&configs[idx].domain, second_int, &configs[idx]);
                }
            }
        } else {
            if is_key_pressed(KeyCode::Up) || is_key_pressed(KeyCode::Right) {
                second_int = (second_int + 0.05).clamp(-1.0, 1.0);
                let (s, t) = build_trajectories(&configs[idx].domain, second_int, &configs[idx]);
                start_points = s;
                trajectories = t;
                if show_3d {
                    phase_trajectories =
                        build_phase(&configs[idx].domain, second_int, &configs[idx]);
                }
            }
            if is_key_pressed(KeyCode::Down) || is_key_pressed(KeyCode::Left) {
                second_int = (second_int - 0.05).clamp(-1.0, 1.0);
                let (s, t) = build_trajectories(&configs[idx].domain, second_int, &configs[idx]);
                start_points = s;
                trajectories = t;
                if show_3d {
                    phase_trajectories =
                        build_phase(&configs[idx].domain, second_int, &configs[idx]);
                }
            }
        }

        // Slider input
        if let Some(new_val) =
            slider.update(second_int, &slider_range, preset.second_integral_label)
        {
            if (new_val - second_int).abs() > 0.0001 {
                second_int = new_val;
                let (s, t) = build_trajectories(&configs[idx].domain, second_int, &configs[idx]);
                start_points = s;
                trajectories = t;
                if show_3d {
                    phase_trajectories =
                        build_phase(&configs[idx].domain, second_int, &configs[idx]);
                }
            }
        }

        // Smooth animation: sweep the second integral back and forth.
        if animate {
            let speed = 0.4; // units per second
            let dt = get_frame_time();
            let step = speed * dt;

            if configs[idx].is_confocal {
                let (lambda_ell, lambda_hyp) = get_ell_hyp(&configs[idx].domain);
                let e = 0.05;
                let ell_min = lambda_ell + e;
                let ell_max = B - e;
                let hyp_min = lambda_hyp + e;
                let hyp_max = A - e;

                let mut next = second_int + anim_dir * step;
                // Bounce off the ends of the valid range.
                if next >= ell_max && next <= hyp_min {
                    // Crossing the gap: jump to the other side.
                    if anim_dir > 0.0 {
                        next = hyp_min;
                    } else {
                        next = ell_max;
                    }
                }
                if next >= hyp_max {
                    next = hyp_max;
                    anim_dir = -1.0;
                }
                if next <= ell_min {
                    next = ell_min;
                    anim_dir = 1.0;
                }
                second_int = next;
            } else {
                let mut next = second_int + anim_dir * step;
                if next >= 1.0 {
                    next = 1.0;
                    anim_dir = -1.0;
                }
                if next <= -1.0 {
                    next = -1.0;
                    anim_dir = 1.0;
                }
                second_int = next;
            }

            let (s, t) = build_trajectories(&configs[idx].domain, second_int, &configs[idx]);
            start_points = s;
            trajectories = t;
            if show_3d {
                phase_trajectories = build_phase(&configs[idx].domain, second_int, &configs[idx]);
            }
        }

        clear_background(BG);

        if show_3d {
            // ---- 3D phase space view ----
            cam3d.handle_input();
            phase3d::draw_axes(&cam3d, w, h, cache.domain_extent);
            // Dense points fill the 2D Liouville torus surface.
            phase3d::draw_phase_points(&phase_trajectories, &cam3d, w, h);
            phase3d::draw_phase_trajectories(&phase_trajectories, &cam3d, w, h);

            let info = format!(
                "{}  |  {} = {:.3}  |  {} trajs  |  [P] 2D  |  [A] anim {}  |  right-drag orbit",
                preset.label,
                preset.second_integral_label,
                second_int,
                phase_trajectories.len(),
                if animate { "on" } else { "off" },
            );
            draw_text(&info, 12.0, 28.0, 18.0, color_u8!(200, 200, 220, 220));
            draw_text(
                "3D phase space: (x, y, θ/π)  —  Liouville torus at fixed second integral",
                12.0,
                h - 12.0,
                14.0,
                color_u8!(150, 150, 170, 140),
            );
        } else {
            // ---- 2D billiard view ----
            draw_boundary(cache, w, h);

            if preset.is_confocal {
                draw_caustic(A, B, second_int, &preset.domain, cam, w, h);
            } else {
                // Draw the starting direction arrow
                let angle = second_int * std::f32::consts::PI;
                draw_start_arrow(preset.start_center, angle, cam, w, h);
            }

            for traj in &trajectories {
                draw_trajectory(traj, cam, w, h);
            }
            for &(p, _) in &start_points {
                draw_start_marker(p, cam, w, h);
            }

            let total: usize = trajectories.iter().map(|t| t.len()).sum();
            let info =
                format!(
                "{}  |  {} = {:.3}  |  {} trajs, {} bounces  |  ↑↓ adjust  |  [A] anim {}  |  [P] 3D  |  Tab next",
                preset.label, preset.second_integral_label, second_int, start_points.len(), total,
                if animate { "on" } else { "off" },
            );
            draw_text(&info, 12.0, 28.0, 18.0, color_u8!(200, 200, 220, 220));

            if preset.is_confocal {
                draw_text(
                    "H = ½|v|²  |  Λ = vx²/a + vy²/b − (x·vy − y·vx)²/(ab)  |  v ⟂ ∇Q_Λ",
                    12.0,
                    h - 12.0,
                    14.0,
                    color_u8!(150, 150, 170, 140),
                );
            } else {
                draw_text(
                    "H = ½|v|²  |  second integral = velocity angle θ = atan2(vy, vx)",
                    12.0,
                    h - 12.0,
                    14.0,
                    color_u8!(150, 150, 170, 140),
                );
            }
        }

        next_frame().await;
    }
}
