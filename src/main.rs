mod domain;
mod presets;
mod quadratic;

use macroquad::prelude::*;

const A: f32 = 4.0;
const B: f32 = 1.0;

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
}

impl CachedDomain {
    fn new(domain: &domain::Domain) -> Self {
        let boundary_pts = domain.sample_boundary(30);
        let corners = domain.corners();
        let camera = Camera::fit_domain(&boundary_pts);
        Self {
            boundary_pts,
            corners,
            camera,
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

fn draw_caustic(a: f32, b: f32, lam: f32, domain: &domain::Domain, cam: &Camera, w: f32, h: f32) {
    let quad = quadratic::confocal(a, b, lam);
    let pts = quad.sample_boundary(80);
    for &p in &pts {
        if domain.contains(p) {
            let s = cam.world_to_screen(p, w, h);
            draw_circle(s.x, s.y, 1.5, color_u8!(255, 255, 100, 120));
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

/// Clamp Λ to the valid range.
/// Ellipse caustic: λ_ell < Λ < b  (caustic inside boundary ellipse)
/// Hyperbola caustic: λ_hyp < Λ < a (caustic narrower than boundary hyperbola, inside)
fn clamp_lambda(lam: f32, domain: &domain::Domain) -> f32 {
    let (lambda_ell, lambda_hyp) = get_boundary_lambdas(domain);
    let e = 0.02;
    if lam >= B - e && lam <= B + e {
        if lam < B {
            B - 2.0 * e
        } else {
            B + 2.0 * e
        }
    } else if lam <= B {
        lam.clamp(lambda_ell + e, B - e)
    } else {
        lam.clamp(lambda_hyp + e, A - e)
    }
}

/// Pick a start point on the caustic that's inside the domain.
/// Returns (position, velocity).
fn start_on_caustic(a: f32, b: f32, lam: f32, domain: &domain::Domain) -> (Vec2, Vec2) {
    let (lambda_ell, lambda_hyp) = get_boundary_lambdas(domain);
    let hint_dir = vec2(0.0, 1.0);

    if lam < B {
        // ---- Ellipse caustic ----
        let rx = (a - lam).sqrt();
        let ry = (b - lam).sqrt();

        // Intersection of caustic ellipse (Λ) with boundary hyperbola (λ_hyp)
        let x2_r = ((a - lam) * (a - lambda_hyp)) / (a - b);
        let y2_r = ((b - lam) * (b - lambda_hyp)) / (b - a);

        if y2_r >= 0.0 && x2_r >= 0.0 {
            let y_int = y2_r.sqrt();
            let x_int = x2_r.sqrt();

            let angle_right = y_int.atan2(x_int);
            let angle_left = y_int.atan2(-x_int);

            let mid_upper = (angle_left + angle_right) / 2.0;
            let p_upper = vec2(rx * mid_upper.cos(), ry * mid_upper.sin());
            let vel =
                quadratic::ConfocalQuadric::velocity_from_caustic(a, b, p_upper, lam, hint_dir);
            return (p_upper, vel);
        }

        let p = vec2(rx, 0.0);
        let vel = quadratic::ConfocalQuadric::velocity_from_caustic(a, b, p, lam, hint_dir);
        (p, vel)
    } else {
        // ---- Hyperbola caustic ----
        // Λ > λ_hyp means the caustic is narrower than the boundary hyperbola.
        // The vertex at (√(a-Λ), 0) is inside the domain.
        let hyp_a = (a - lam).sqrt();
        let p = vec2(hyp_a, 0.0);
        let vel = quadratic::ConfocalQuadric::velocity_from_caustic(a, b, p, lam, hint_dir);
        (p, vel)
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

    let trace_current = |domain: &domain::Domain, lam: f32| {
        let (pos, vel) = start_on_caustic(A, B, lam, domain);
        domain.trace(pos, vel, 300)
    };

    let mut segments = trace_current(&configs[idx].domain, lambda);

    loop {
        let (w, h) = (screen_width(), screen_height());
        let cache = &caches[idx];
        let cam = &cache.camera;

        if is_key_pressed(KeyCode::Tab) || is_key_pressed(KeyCode::Space) {
            idx = (idx + 1) % configs.len();
            lambda = 0.2;
            segments = trace_current(&configs[idx].domain, lambda);
        }

        if is_key_pressed(KeyCode::Up) || is_key_pressed(KeyCode::Right) {
            lambda = clamp_lambda(lambda + 0.05, &configs[idx].domain);
            segments = trace_current(&configs[idx].domain, lambda);
        }
        if is_key_pressed(KeyCode::Down) || is_key_pressed(KeyCode::Left) {
            lambda = clamp_lambda(lambda - 0.05, &configs[idx].domain);
            segments = trace_current(&configs[idx].domain, lambda);
        }

        clear_background(BG);

        draw_boundary(cache, w, h);
        draw_caustic(A, B, lambda, &configs[idx].domain, cam, w, h);
        let (start_pos, _) = start_on_caustic(A, B, lambda, &configs[idx].domain);
        draw_trajectory(&segments, cam, w, h);
        draw_start_marker(start_pos, cam, w, h);

        let caustic_type = if lambda < B { "ellipse" } else { "hyperbola" };
        let info = format!(
            "{}  |  Λ = {:.3} ({})  |  {} bounces  |  ↑↓ Λ  |  Tab next",
            configs[idx].label,
            lambda,
            caustic_type,
            segments.len(),
        );
        draw_text(&info, 12.0, 28.0, 18.0, color_u8!(200, 200, 220, 220));

        draw_text(
            "H = ½|v|²  |  Λ = vx²/a + vy²/b − (x·vy − y·vx)²/(ab)  |  v ⟂ ∇Q_Λ",
            12.0,
            h - 12.0,
            14.0,
            color_u8!(150, 150, 170, 140),
        );

        next_frame().await;
    }
}
