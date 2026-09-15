//! 2D drawing for the billiard view.
//!
//! Owns the 2D camera, the cached boundary, and every `draw_*` primitive.  Kept
//! out of the binary's `main` so the app loop only wires input → state → draw
//! and the rendering stays independently readable (and testable against a
//! headless target if needed).

use crate::confocal::TorusRegime;
use crate::torus::ConfocalParams;
use crate::{domain, quadratic};
use macroquad::prelude::*;

// ----------------------------------------------------------------
// Camera (2D)
// ----------------------------------------------------------------
/// Fraction of the window height reserved at the bottom for the Λ slider and
/// (when shown) the molecule strip, kept clear of the billiard drawing.
const BOTTOM_UI_MARGIN: f32 = 0.2;

pub struct Camera {
    centre: Vec2,
    half_height: f32,
}

impl Camera {
    pub fn fit_domain(pts: &[Vec2]) -> Self {
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
    pub fn world_to_screen(&self, p: Vec2, win_w: f32, win_h: f32) -> Vec2 {
        // Reserve the bottom of the window for the slider / molecule strip so
        // a vertically-dominant domain (e.g. the "Thin" confocal ellipse)
        // never draws its own boundary into the UI band — the old mapping
        // filled the *whole* window height, so a tall enough shape's
        // boundary curve landed right on top of the Λ slider.
        let usable_h = win_h * (1.0 - BOTTOM_UI_MARGIN);
        let aspect = win_w / usable_h;
        let sx = ((p.x - self.centre.x) / (self.half_height * aspect) + 1.0) * 0.5 * win_w;
        let sy = (-(p.y - self.centre.y) / self.half_height + 1.0) * 0.5 * usable_h;
        vec2(sx, sy)
    }
}

// ----------------------------------------------------------------
// Cached boundary
// ----------------------------------------------------------------
pub struct CachedDomain {
    pub boundary_pts: Vec<Vec2>,
    pub corners: Vec<Vec2>,
    pub camera: Camera,
    pub domain_extent: f32,
}

impl CachedDomain {
    pub fn new(domain: &domain::Domain) -> Self {
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
pub const BG: Color = color_u8!(15, 15, 35, 255);
const BOUNDARY: Color = color_u8!(180, 220, 255, 200);
const CORNER: Color = color_u8!(255, 200, 100, 200);

pub(crate) fn hsl_to_rgb(hue: f32, s: f32, l: f32) -> Color {
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

pub fn draw_boundary(cache: &CachedDomain, w: f32, h: f32) {
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

pub fn draw_trajectory(segs: &[(Vec2, Vec2)], cam: &Camera, w: f32, h: f32) {
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

pub fn draw_start_marker(pos: Vec2, cam: &Camera, w: f32, h: f32) {
    let s = cam.world_to_screen(pos, w, h);
    draw_circle_lines(s.x, s.y, 7.0, 2.5, YELLOW);
    draw_circle(s.x, s.y, 3.5, YELLOW);
}

/// Draw caustic curves and degenerate reference lines (for confocal domains).
pub fn draw_caustic(
    cf: ConfocalParams,
    lam: f32,
    domain: &domain::Domain,
    cam: &Camera,
    w: f32,
    h: f32,
) {
    let a = cf.a;
    let b = cf.b;
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
    let quad = quadratic::confocal(cf, lam);
    let pts = quad.sample_boundary(80);
    for &p in &pts {
        if domain.contains(p) {
            let s = cam.world_to_screen(p, w, h);
            draw_circle(s.x, s.y, 1.8, color_u8!(255, 255, 100, 160));
        }
    }
}

/// Draw a starting direction arrow for polyline domains.
pub fn draw_start_arrow(center: Vec2, angle: f32, cam: &Camera, w: f32, h: f32) {
    let dir = vec2(angle.cos(), angle.sin());
    let tip = center + dir * 0.3;
    let s1 = cam.world_to_screen(center, w, h);
    let s2 = cam.world_to_screen(tip, w, h);
    draw_line(s1.x, s1.y, s2.x, s2.y, 2.5, YELLOW);
    draw_circle(s1.x, s1.y, 3.5, YELLOW);
}

/// Draw the first `n_bounces` segments of a trajectory bold red on the
/// billiard, so one can see exactly where each torus's orbit lives on the
/// domain itself.
pub fn draw_trajectory_red(segs: &[(Vec2, Vec2)], n_bounces: usize, cam: &Camera, w: f32, h: f32) {
    for &(a, b) in segs.iter().take(n_bounces) {
        let sa = cam.world_to_screen(a, w, h);
        let sb = cam.world_to_screen(b, w, h);
        draw_line(sa.x, sa.y, sb.x, sb.y, 3.0, RED);
        draw_circle(sa.x, sa.y, 3.5, RED);
    }
}

/// One representative trajectory per torus on the billiard, using the same
/// torus-index rule as the torus mapping (`TorusRegime`).
pub fn one_per_torus(start_points: &[(Vec2, Vec2)], regime: TorusRegime) -> Vec<usize> {
    let mut seen: Vec<u32> = Vec::new();
    let mut out = Vec::new();
    for (i, &(p, v)) in start_points.iter().enumerate() {
        let key = regime.index_of(p.x, p.y, v.x, v.y);
        if !seen.contains(&key) {
            seen.push(key);
            out.push(i);
        }
    }
    out
}
