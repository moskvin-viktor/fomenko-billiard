mod domain;
mod presets;

use domain::Domain;
use macroquad::prelude::*;
use presets::Preset;

// ----------------------------------------------------------------
// Coordinate helpers – auto-fit any domain in the window
// ----------------------------------------------------------------
struct Camera {
    centre: Vec2,
    half_height: f32,
}

impl Camera {
    fn fit_domain(domain: &Domain) -> Self {
        let pts = domain.sample_boundary(40);
        if pts.is_empty() {
            return Self {
                centre: vec2(0.0, 0.0),
                half_height: 2.0,
            };
        }
        let min_x = pts.iter().map(|p| p.x).reduce(f32::min).unwrap_or(-1.0);
        let max_x = pts.iter().map(|p| p.x).reduce(f32::max).unwrap_or(1.0);
        let min_y = pts.iter().map(|p| p.y).reduce(f32::min).unwrap_or(-1.0);
        let max_y = pts.iter().map(|p| p.y).reduce(f32::max).unwrap_or(1.0);

        let cx = (min_x + max_x) / 2.0;
        let cy = (min_y + max_y) / 2.0;
        let hw = (max_x - min_x) / 2.0;
        let hh = (max_y - min_y) / 2.0;

        // Add 20% padding
        let half = hh.max(hw) * 1.2;
        Self {
            centre: vec2(cx, cy),
            half_height: half.max(0.5),
        }
    }

    fn world_to_screen(&self, p: Vec2, win_w: f32, win_h: f32) -> Vec2 {
        let aspect = win_w / win_h;
        let sx = ((p.x - self.centre.x) / (self.half_height * aspect) + 1.0) * 0.5 * win_w;
        let sy = (-(p.y - self.centre.y) / self.half_height + 1.0) * 0.5 * win_h;
        vec2(sx, sy)
    }

    fn screen_to_world(&self, pos: Vec2, win_w: f32, win_h: f32) -> Vec2 {
        let aspect = win_w / win_h;
        let wx = (pos.x / win_w - 0.5) * 2.0 * self.half_height * aspect + self.centre.x;
        let wy = -(pos.y / win_h - 0.5) * 2.0 * self.half_height + self.centre.y;
        vec2(wx, wy)
    }
}

// ----------------------------------------------------------------
// Drawing helpers
// ----------------------------------------------------------------
const BG: Color = color_u8!(15, 15, 35, 255);
const BOUNDARY: Color = color_u8!(180, 220, 255, 200);
const CORNER: Color = color_u8!(255, 200, 100, 200);

fn draw_boundary(domain: &Domain, cam: &Camera, w: f32, h: f32) {
    let pts = domain.sample_boundary(30);
    let sp: Vec<Vec2> = pts.iter().map(|p| cam.world_to_screen(*p, w, h)).collect();

    if sp.is_empty() {
        return;
    }

    for i in 0..sp.len() {
        let j = (i + 1) % sp.len();
        draw_line(sp[i].x, sp[i].y, sp[j].x, sp[j].y, 2.0, BOUNDARY);
    }

    // Mark corners
    for c in domain.corners() {
        let s = cam.world_to_screen(c, w, h);
        draw_circle(s.x, s.y, 4.0, CORNER);
    }
}

fn draw_trajectory(segs: &[(Vec2, Vec2)], cam: &Camera, w: f32, h: f32) {
    let n = segs.len();
    for (i, (a, b)) in segs.iter().enumerate() {
        let t = i as f32 / n.max(1) as f32;
        let hue = 30.0 + t * 200.0;
        let color = hsl_to_rgb(hue, 0.85, 0.55);

        let sa = cam.world_to_screen(*a, w, h);
        let sb = cam.world_to_screen(*b, w, h);

        draw_line(sa.x, sa.y, sb.x, sb.y, 1.8, color);
        draw_circle(sb.x, sb.y, 2.5 + (1.0 - t) * 2.0, hsl_to_rgb(hue, 0.8, 0.7));
    }
}

fn draw_start_marker(pos: Vec2, cam: &Camera, w: f32, h: f32) {
    let s = cam.world_to_screen(pos, w, h);
    draw_circle_lines(s.x, s.y, 7.0, 2.5, YELLOW);
    draw_circle(s.x, s.y, 3.5, YELLOW);
}

fn draw_velocity_arrow(pos: Vec2, vel: Vec2, cam: &Camera, w: f32, h: f32) {
    let s = cam.world_to_screen(pos, w, h);
    let tip = cam.world_to_screen(pos + vel * 0.15, w, h);
    draw_line(s.x, s.y, tip.x, tip.y, 2.0, color_u8!(255, 255, 100, 180));
}

fn hsl_to_rgb(h: f32, s: f32, l: f32) -> Color {
    let h = h / 360.0;
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

// ----------------------------------------------------------------
// Random (simple LCG)
// ----------------------------------------------------------------
use std::cell::Cell;

thread_local! {
    static SEED: Cell<u64> = const { Cell::new(42) };
}

fn rand_f32() -> f32 {
    SEED.with(|seed| {
        let s = seed
            .get()
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        seed.set(s);
        ((s >> 33) as u32 as f32) / (u32::MAX as f32)
    })
}

// ----------------------------------------------------------------
// Entry point
// ----------------------------------------------------------------
#[macroquad::main("Mathematical Billiards")]
async fn main() {
    let configs = presets::all_presets();
    let mut idx = 0;
    let mut custom_pos: Option<Vec2> = None;
    let mut custom_vel: Option<Vec2> = None;

    let trace = |p: &Preset| -> (Vec<(Vec2, Vec2)>, bool) {
        let segs = p.domain.trace(p.start, p.vel, 300);
        let ok = !segs.is_empty();
        (segs, ok)
    };

    let (mut segments, mut valid) = trace(&configs[idx]);

    loop {
        let (w, h) = (screen_width(), screen_height());
        let cam = Camera::fit_domain(&configs[idx].domain);

        // ---- Input ----
        if is_key_pressed(KeyCode::Tab) || is_key_pressed(KeyCode::Space) {
            idx = (idx + 1) % configs.len();
            custom_pos = None;
            custom_vel = None;
            let r = trace(&configs[idx]);
            segments = r.0;
            valid = r.1;
        }

        if is_mouse_button_pressed(MouseButton::Left) {
            let mouse = mouse_position();
            let world = cam.screen_to_world(vec2(mouse.0, mouse.1), w, h);
            if configs[idx].domain.contains(world) {
                custom_pos = Some(world);
                let angle = rand_f32() * 2.0 * std::f32::consts::PI;
                let speed = 0.5 + rand_f32() * 0.8;
                custom_vel = Some(vec2(angle.cos(), angle.sin()) * speed);
                let segs = configs[idx].domain.trace(world, custom_vel.unwrap(), 300);
                segments = segs;
                valid = !segments.is_empty();
            }
        }

        // ---- Draw ----
        clear_background(BG);

        let start = custom_pos.unwrap_or(configs[idx].start);
        let vel = custom_vel.unwrap_or(configs[idx].vel);

        draw_boundary(&configs[idx].domain, &cam, w, h);

        if valid {
            draw_trajectory(&segments, &cam, w, h);
        }

        draw_start_marker(start, &cam, w, h);
        draw_velocity_arrow(start, vel, &cam, w, h);

        // ---- HUD ----
        let label = custom_pos
            .map(|_| "Custom (click inside to place)")
            .unwrap_or(configs[idx].label);

        let info = format!(
            "{}  |  1 bounce  |  [Tab] next  |  click to place start",
            label,
        );
        draw_text(&info, 12.0, 28.0, 18.0, color_u8!(200, 200, 220, 220));
        draw_text(
            "Single bounce  |  π/2 corner → return  |  3π/2 corner → deflect",
            12.0,
            h - 12.0,
            14.0,
            color_u8!(150, 150, 170, 140),
        );

        next_frame().await;
    }
}
