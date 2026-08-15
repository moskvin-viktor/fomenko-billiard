//! Sampling and 3D drawing for pseudo-integrable level sets.
//!
//! The torus path (`phase3d`) forces every level onto a donut via
//! `torus_embed`.  For the L that is wrong on both counts the user hit:
//! torus levels didn't show a torus, and genus-2 levels didn't look like
//! genus 2.  This module samples a trajectory through the **flat** chart
//! (`pseudo::to_flat`) and draws it with the correct embedding:
//!
//! - [`Level::Torus`] — the flat rectangle `(u1, u2)` with opposite edges
//!   identified → a donut (via `torus_angles`).
//! - [`Level::GenusSurface`] — the unfolded cross 12-gon with the four sheets
//!   stacked in z (via `cross_embed`).

use super::classify::Level;
use super::map::to_flat;
use super::render::{pretzel_embed, torus_angles};
use crate::torus::{ConfocalParams, PhaseSample};
use macroquad::prelude::*;

/// A phase-space point mapped into flat coordinates for rendering.
#[derive(Clone, Copy, Debug)]
pub struct FlatPhasePoint {
    /// The un-normalized length coordinate `u₁(λ₁)`.
    pub u1: f32,
    /// The un-normalized length coordinate `u₂(λ₂)`.
    pub u2: f32,
    /// The momentum sheet `(sign d₁, sign d₂)`.
    pub sheet: (i8, i8),
    /// The connected component (0 for a torus).
    pub component: u32,
}

/// Sample a trajectory densely through the flat chart, mapping every interior
/// point with `to_flat` for the given (precomputed) level.
///
/// The level is fixed along a trajectory (λc conserved), so it is built once
/// and reused for every sample.
pub fn sample_flat_trajectory(
    domain: &crate::domain::Domain,
    p0: Vec2,
    v0: Vec2,
    max_steps: usize,
    per_seg: usize,
    level: &Level,
) -> Vec<FlatPhasePoint> {
    let cf = ConfocalParams::standard();
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
            let sample = PhaseSample::new(q.x, q.y, v.x, v.y);
            if let Ok(flat) = to_flat(&sample, &cf, level) {
                pts.push(FlatPhasePoint {
                    u1: flat.u1,
                    u2: flat.u2,
                    sheet: flat.sheet,
                    component: flat.component,
                });
            }
        }
        v = domain.reflect(hit, v, idx);
        p = hit + 1e-4 * v.normalize();
    }
    pts
}

/// Draw a set of flat trajectories with the correct embedding for the level.
///
/// * `Level::Torus` — each point maps to torus angles `(θ1, θ2)` from the flat
///   rectangle, then to a donut (the doc's "revert to normalized angles").
/// * `Level::GenusSurface` — each point maps into the unfolded cross with its
///   sheet stacked in z.
pub fn draw_flat(
    trajectories: &[Vec<FlatPhasePoint>],
    level: &Level,
    cam: &crate::phase3d::OrbitCamera3,
    win_w: f32,
    win_h: f32,
) {
    match level {
        Level::Torus { shape, .. } => {
            let (w1, w2) = *shape;
            for traj in trajectories {
                for pt in traj {
                    let (t1, t2) = torus_angles(pt.u1, pt.u2, w1, w2);
                    let pos = donut(t1, t2);
                    let s = cam.project(pos, win_w, win_h);
                    if s.z < -0.1 {
                        continue;
                    }
                    draw_circle(s.x, s.y, 1.8, color_u8!(200, 120, 80, 200));
                }
            }
        }
        Level::GenusSurface { region, .. } => {
            // Cross moduli: A1 = u1(α₁), A2 = u1(α₂), B1 = u2(β₂), B2 = u2(β₃).
            let m = cross_moduli(region);
            for traj in trajectories {
                for pt in traj {
                    let pos = pretzel_embed(pt.u1, pt.u2, pt.sheet.0, pt.sheet.1, m);
                    let s = cam.project(pos, win_w, win_h);
                    if s.z < -0.1 {
                        continue;
                    }
                    // Color by sheet so the four layers are distinguishable.
                    let color = sheet_color(pt.sheet);
                    draw_circle(s.x, s.y, 1.8, color);
                }
            }
        }
        _ => {}
    }
}

/// Standard donut embedding (same as `phase3d::torus_embed` without the index
/// offset), used for torus levels.
fn donut(theta1: f32, theta2: f32) -> Vec3 {
    let (r_major, r_minor) = (1.6f32, 0.6f32);
    let (s1, c1) = theta1.sin_cos();
    let (s2, c2) = theta2.sin_cos();
    vec3(
        (r_major + r_minor * c1) * c2,
        (r_major + r_minor * c1) * s2,
        r_minor * s1,
    )
}

/// Cross moduli from the accessible region: `A1 = u1(α₁)`, `A2 = u1(α₂)`,
/// `B1 = u2(β₂)`, `B2 = u2(β₃)`.  These are the flat extents of the two legs.
fn cross_moduli(region: &super::classify::AccessibleRegion) -> super::render::CrossModuli {
    // The region's u1/u2 grid lines are monotone; the cross uses the full
    // extent of each leg.  For the standard L: A1 is the tall-leg width
    // (u1 at the reflex corner), A2 the base width, B1 the base height,
    // B2 the tall-leg height.
    let u1 = &region.u1;
    let u2 = &region.u2;
    let a1 = u1.get(1).copied().unwrap_or(0.0); // u1(α₁)
    let a2 = u1.last().copied().unwrap_or(a1); // u1(α₂)
    let b1 = u2.get(1).copied().unwrap_or(0.0); // u2(β₂)
    let b2 = u2.last().copied().unwrap_or(b1); // u2(β₃)
    super::render::CrossModuli { a1, a2, b1, b2 }
}

/// A color for a momentum sheet, so the four layers of the cross are distinct.
fn sheet_color(sheet: (i8, i8)) -> Color {
    let idx = ((sheet.0 > 0) as u8) * 2 + ((sheet.1 > 0) as u8);
    let hue = [0.0, 60.0, 180.0, 300.0][idx as usize];
    hsl(hue)
}

/// HSL → RGB helper (matches the renderer's palette).
fn hsl(hue: f32) -> Color {
    let s: f32 = 0.85;
    let l: f32 = 0.55;
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
