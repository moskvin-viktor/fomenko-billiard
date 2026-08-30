//! Sampling and 3D drawing for pseudo-integrable level sets.
//!
//! The torus path (`phase3d`) forces every level onto a donut via
//! `torus_embed`.  For the L that is wrong on both counts the user hit:
//! torus levels didn't show a torus, and genus-2 levels didn't look like
//! genus 2.  This module samples a trajectory through the **flat** chart
//! (`pseudo::to_flat`) and draws it with the **unified morph embedding**
//! ([`super::morph::unified_with_normal`]): one continuous surface family —
//! donut + shrinking handle — for torus and genus-2 levels alike, so the
//! genus bifurcations render as a smooth pinch instead of a surface swap.

use super::classify::Level;
use super::geometry::LevelGeometry;
use super::map::to_flat;
use super::morph::{pinch_point_3d, unified_with_normal, FlatMorph};
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

/// Sample the boundary preimage π⁻¹ of a pseudo-integrable level: the walls
/// and the caustic curves of the billiard, mapped through the flat chart.
/// Returns a `FlatPhasePoint` per sample, with `true` = caustic, `false` = wall.
pub fn sample_flat_boundary(
    domain: &crate::domain::Domain,
    lam: f32,
    level: &Level,
) -> Vec<(FlatPhasePoint, bool)> {
    let cf = ConfocalParams::standard();
    let mut out = Vec::new();

    let mut push = |p: Vec2, is_caustic: bool| {
        for vel in crate::phase3d::velocities_for_lambda(p, lam, cf.a, cf.b) {
            let sample = PhaseSample::new(p.x, p.y, vel.x, vel.y);
            if let Ok(flat) = to_flat(&sample, &cf, level) {
                out.push((
                    FlatPhasePoint {
                        u1: flat.u1,
                        u2: flat.u2,
                        sheet: flat.sheet,
                        component: flat.component,
                    },
                    is_caustic,
                ));
            }
        }
    };

    // Walls: the physical boundary of the billiard.
    for p in domain.sample_boundary(24) {
        push(p, false);
    }

    // Caustic arcs Q_Λ = 0 inside the domain (turning-point curves).  At a
    // caustic point the velocity is exactly tangent, given by rotating the
    // gradient; the two ± directions are the phase-space sheets meeting there.
    if !((lam - cf.b).abs() < 0.03 || (lam - cf.a).abs() < 0.03) {
        let quad = crate::quadratic::confocal(cf, lam);
        for p in quad.sample_boundary(60) {
            if domain.contains(p) {
                let g = quad.grad(p);
                let tan = vec2(-g.y, g.x).normalize();
                let sample_caustic = PhaseSample::new(p.x, p.y, tan.x, tan.y);
                if let Ok(flat) = to_flat(&sample_caustic, &cf, level) {
                    out.push((
                        FlatPhasePoint {
                            u1: flat.u1,
                            u2: flat.u2,
                            sheet: flat.sheet,
                            component: flat.component,
                        },
                        true,
                    ));
                }
                let sample_caustic = PhaseSample::new(p.x, p.y, -tan.x, -tan.y);
                if let Ok(flat) = to_flat(&sample_caustic, &cf, level) {
                    out.push((
                        FlatPhasePoint {
                            u1: flat.u1,
                            u2: flat.u2,
                            sheet: flat.sheet,
                            component: flat.component,
                        },
                        true,
                    ));
                }
            }
        }
    }

    out
}

/// Draw a set of flat trajectories on the unified morph embedding: a Lambert
/// depth-shaded point cloud plus hue-swept trajectory lines (matching the
/// smooth confocal tori), and — when the handle is nearly collapsed — a
/// bright marker at the pinch point so the singular level is legible.
pub fn draw_flat(
    trajectories: &[Vec<FlatPhasePoint>],
    geo: &LevelGeometry,
    morph: &FlatMorph,
    cam: &crate::phase3d::OrbitCamera3,
    win_w: f32,
    win_h: f32,
) {
    if geo.u1_extent <= 0.0 || geo.u2_extent <= 0.0 {
        return;
    }
    let light = vec3(0.4, 0.6, 0.7).normalize();

    // Dense shaded point cloud, decimated to a budget for interactivity.
    const POINT_BUDGET: usize = 18_000;
    let total: usize = trajectories.iter().map(|t| t.len()).sum();
    let step = crate::phase3d::decimation_step(total, POINT_BUDGET);
    for traj in trajectories {
        for (pt_idx, pt) in traj.iter().enumerate() {
            if pt_idx % step != 0 {
                continue;
            }
            let (pos, n) = unified_with_normal(pt.u1, pt.u2, pt.sheet.0, pt.sheet.1, geo, morph);
            let s = cam.project(pos, win_w, win_h);
            if s.z < -0.1 {
                continue;
            }
            let lambert = (n.dot(light)).max(0.0);
            let shade = 0.16 + 0.38 * lambert;
            let depth = (-s.z).clamp(0.5, 5.0);
            let alpha = (0.22 + 0.35 * (1.0 - (depth - 0.5) / 4.5)).clamp(0.05, 0.6);
            let radius = 1.4 + 0.6 * (1.0 - (depth - 0.5) / 4.5);
            let r = (255.0 * shade) as u8;
            let g = (170.0 * shade) as u8;
            let b = (90.0 * shade) as u8;
            draw_circle(s.x, s.y, radius, color_u8!(r, g, b, (alpha * 255.0) as u8));
        }
    }

    // Hue-swept trajectory lines so individual orbits are visible,
    // mirroring the smooth case's `draw_phase_trajectories`.
    for traj in trajectories {
        let n = traj.len();
        for (i, win) in traj.windows(2).enumerate() {
            let a = win[0];
            let b = win[1];
            let (pa, _na) = unified_with_normal(a.u1, a.u2, a.sheet.0, a.sheet.1, geo, morph);
            let (pb, _nb) = unified_with_normal(b.u1, b.u2, b.sheet.0, b.sheet.1, geo, morph);
            let sa = cam.project(pa, win_w, win_h);
            let sb = cam.project(pb, win_w, win_h);
            if sa.z < -0.1 || sb.z < -0.1 {
                continue;
            }
            let t = i as f32 / n.max(1) as f32;
            let hue = 30.0 + t * 200.0;
            let color = crate::render::hsl_to_rgb(hue, 0.85, 0.55);
            let depth = (-sa.z.min(sb.z)).clamp(0.5, 5.0);
            let alpha = (0.3 + 0.5 * (1.0 - (depth - 0.5) / 4.5)).clamp(0.1, 0.8);
            let mut c = color;
            c.a = alpha;
            draw_line(sa.x, sa.y, sb.x, sb.y, 1.5, c);
        }
    }

    // Singular pinch marker: when the handle is (nearly) collapsed, the level
    // is a pinched torus — flag the pinch point brightly.
    if morph.handle_t < 0.05 {
        if let Some(p) = pinch_point_3d(geo, morph) {
            let s = cam.project(p, win_w, win_h);
            if s.z > -0.1 {
                draw_circle(s.x, s.y, 7.0, color_u8!(255, 240, 120, 230));
                draw_circle(s.x, s.y, 3.5, color_u8!(255, 255, 255, 255));
            }
        }
    }
}

/// Draw bold red example trajectories on the flat surface (torus or genus-2),
/// mirroring the smooth case's `draw_torus_highlights`.  `highlights` are short
/// (few-bounce) traces, one per start point, embedded onto the surface.
pub fn draw_flat_highlights(
    highlights: &[Vec<FlatPhasePoint>],
    geo: &LevelGeometry,
    morph: &FlatMorph,
    cam: &crate::phase3d::OrbitCamera3,
    win_w: f32,
    win_h: f32,
) {
    for traj in highlights {
        if traj.len() < 2 {
            continue;
        }
        let projected: Vec<Vec3> = traj
            .iter()
            .map(|pt| {
                let p = unified_with_normal(pt.u1, pt.u2, pt.sheet.0, pt.sheet.1, geo, morph).0;
                cam.project(p, win_w, win_h)
            })
            .collect();
        for win in projected.windows(2) {
            let a = win[0];
            let b = win[1];
            if a.z < -0.1 || b.z < -0.1 {
                continue;
            }
            draw_line(a.x, a.y, b.x, b.y, 4.0, RED);
        }
        for s in &projected {
            if s.z < -0.1 {
                continue;
            }
            draw_circle(s.x, s.y, 4.0, RED);
        }
    }
}

/// Draw the boundary preimage π⁻¹ on the flat map: cyan for the domain walls,
/// orange for the caustic.  Same colour scheme as the smooth torus boundary,
/// so the two views agree.
pub fn draw_flat_boundary(
    boundary: &[(FlatPhasePoint, bool)],
    geo: &LevelGeometry,
    morph: &FlatMorph,
    cam: &crate::phase3d::OrbitCamera3,
    win_w: f32,
    win_h: f32,
) {
    if boundary.is_empty() {
        return;
    }
    let mut pts: Vec<(Vec3, bool)> = boundary
        .iter()
        .map(|(pt, is_ca)| {
            let p = unified_with_normal(pt.u1, pt.u2, pt.sheet.0, pt.sheet.1, geo, morph).0;
            (cam.project(p, win_w, win_h), *is_ca)
        })
        .collect();
    pts.sort_by(|a, b| {
        a.0.z
            .partial_cmp(&b.0.z)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    for (s, is_ca) in &pts {
        if s.z < -0.1 {
            continue;
        }
        let depth = (-s.z).clamp(0.5, 5.0);
        let alpha = (0.65 + 0.35 * (1.0 - (depth - 0.5) / 4.5)).clamp(0.4, 1.0);
        let (r, g, b) = if *is_ca {
            (255, 130, 20)
        } else {
            (70, 210, 255)
        };
        draw_circle(s.x, s.y, 3.0, color_u8!(r, g, b, (alpha * 255.0) as u8));
    }
}

