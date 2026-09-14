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
use super::map::{reflex_corners_in_flat, to_flat};
use super::morph::{pinch_point_3d, unified_with_normal, FlatMorph};
use crate::cached_render::CachedSurfaceRender;
use crate::torus::{ConfocalParams, PhaseSample};
use macroquad::prelude::*;

/// How close (in flat `(u1,u2)` units) a sample may come to a reflex-corner
/// image before the trajectory is terminated (see [`sample_flat_trajectory`]).
/// A reflex corner is a 3-pronged singularity (doc §8 of
/// `confocal_L_pseudo_integrable.md`): an orbit that reaches it has three
/// equally valid continuations and no canonical choice, so the billiard
/// reflection law silently picks one of them — the resulting direction can
/// differ by `O(1)` from what a trajectory passing just on the other side of
/// the corner would take. That is not a rendering artifact to smooth over;
/// the orbit itself is undefined past this point, so it must stop instead of
/// continuing through an arbitrary prong.
const CORNER_EPS: f32 = 0.03;

/// Above this ratio of (embedded 3D jump) / (raw flat-chart step), a
/// same-branch connecting line is a numerical artifact, not real motion —
/// see [`is_seam`].
///
/// The main lobe's angular speed is `π / lsplit` (tube) or `π / u2_extent`
/// (major), a per-*level* constant that is often well above 1 (e.g. `lsplit
/// ≈ 0.23` on one observed genus-2 level gives speed ≈ 13.6) — so *legitimate*
/// same-branch motion can already reach stretch ratios of ~10–20 on an
/// ordinarily-proportioned level; that is not a seam, just a steep chart.
/// Only right at a band end, where `u2_extent` (or `lsplit`) has shrunk
/// toward zero, does the ratio blow up to the hundreds for a physically tiny
/// step — measured ~270 at `u2_extent ≈ 0.013`. `100` sits with wide margin
/// above the legitimate case and below the pathological one.
const MAX_STRETCH_RATIO: f32 = 100.0;

/// Whether the straight line between two temporally-adjacent trajectory
/// samples' embedded positions would be a rendering artifact rather than a
/// real traced path: the embedded jump is wildly disproportionate to how far
/// the samples actually moved in the flat chart (see [`MAX_STRETCH_RATIO`]).
///
/// This used to also force a seam on every main-lobe/handle crossing, because
/// [`unified_with_normal`] collapsed the whole attach band onto a single
/// pinch point there, so a crossing step could produce an arbitrarily large,
/// spurious jump. That embedding bug is fixed (the handle now attaches along
/// the actual attach curve, matching the main lobe exactly at the seam), so a
/// legitimate crossing step is small like any other and the ratio test alone
/// catches genuine artifacts.
fn is_seam(a: &FlatPhasePoint, b: &FlatPhasePoint, pa: Vec3, pb: Vec3, _geo: &LevelGeometry) -> bool {
    let raw = ((a.u1 - b.u1).powi(2) + (a.u2 - b.u2).powi(2)).sqrt();
    (pa - pb).length() > MAX_STRETCH_RATIO * raw.max(1e-6)
}

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
    let corners = reflex_corners_in_flat(level);
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
                if near_reflex_corner(flat.u1, flat.u2, &corners) {
                    // Approaching the 3-pronged singularity: no continuation
                    // is canonical past this point (doc §8) — stop the whole
                    // trajectory here rather than let the reflection law pick
                    // an arbitrary prong and jump to an unrelated part of the
                    // surface.
                    return pts;
                }
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

/// Whether a flat sample sits within [`CORNER_EPS`] of any reflex-corner
/// image of the level.
fn near_reflex_corner(u1: f32, u2: f32, corners: &[(f32, f32)]) -> bool {
    corners
        .iter()
        .any(|&(cu1, cu2)| ((u1 - cu1).powi(2) + (u2 - cu2).powi(2)).sqrt() < CORNER_EPS)
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
            if is_seam(&a, &b, pa, pb, geo) {
                // The straight 3D line would cut through the surface's
                // interior instead of tracing along it; drop the segment,
                // not the points (the point cloud above already drew both
                // ends).
                continue;
            }
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
        let embedded: Vec<Vec3> = traj
            .iter()
            .map(|pt| unified_with_normal(pt.u1, pt.u2, pt.sheet.0, pt.sheet.1, geo, morph).0)
            .collect();
        for (win, pts) in embedded.windows(2).zip(traj.windows(2)) {
            let (pa, pb) = (win[0], win[1]);
            if is_seam(&pts[0], &pts[1], pa, pb, geo) {
                // Seam — see the comment in `draw_flat`.
                continue;
            }
            let sa = cam.project(pa, win_w, win_h);
            let sb = cam.project(pb, win_w, win_h);
            if sa.z < -0.1 || sb.z < -0.1 {
                continue;
            }
            draw_line(sa.x, sa.y, sb.x, sb.y, 4.0, RED);
        }
        for &p in &embedded {
            let s = cam.project(p, win_w, win_h);
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

/// Cached render-to-texture wrapper around [`draw_flat`]/[`draw_flat_boundary`],
/// on the same [`CachedSurfaceRender`] the smooth-torus path uses
/// (`torus_render::TorusRender`).  Before this, `draw_flat` was called
/// directly every frame regardless of camera motion; the smooth path's
/// "orbit freeze" fix (`docs/known_issues.md` #1) never covered the L-shape
/// view.  This gives the two paths identical caching: rasterize only when the
/// camera moved or the level's geometry/morph state changed, otherwise blit.
pub struct FlatRender {
    cache: CachedSurfaceRender,
}

/// Everything [`FlatRender::draw`] needs apart from the trajectories
/// themselves: the camera, window size, and the level's geometry/morph state.
/// Bundled for the same reason as `torus_render::DrawContext` — so the draw
/// call stays small despite the two paths needing the same shape of inputs.
pub struct FlatDrawContext<'a> {
    pub cam: &'a crate::phase3d::OrbitCamera3,
    pub win_w: f32,
    pub win_h: f32,
    pub geo: &'a LevelGeometry,
    pub morph: &'a FlatMorph,
}

impl FlatRender {
    pub fn new() -> Self {
        Self {
            cache: CachedSurfaceRender::new(),
        }
    }

    /// Access the cached offscreen texture (for tests / diagnostics).
    pub fn texture(&self) -> Option<&Texture2D> {
        self.cache.texture()
    }

    /// A cheap fingerprint of the flat geometry: point/line counts, a content
    /// sample, and the level's intrinsic geometry + morph state (both change
    /// with the caustic level `λc`, so a level sweep must re-rasterize even
    /// though the point cloud's raw `(u1,u2)` values don't move).
    fn geometry_key(
        trajectories: &[Vec<FlatPhasePoint>],
        boundary: &[(FlatPhasePoint, bool)],
        geo: &LevelGeometry,
        morph: &FlatMorph,
    ) -> (usize, usize, u32) {
        let n_trajs = trajectories.len();
        let n_pts: usize = trajectories.iter().map(|t| t.len()).sum::<usize>() + boundary.len();
        let mut sample = 0u32;
        for t in trajectories.iter().take(4) {
            if let Some(p) = t.first() {
                sample = sample
                    .wrapping_mul(31)
                    .wrapping_add(p.u1.to_bits() ^ p.u2.to_bits());
            }
        }
        for bits in [
            geo.origin.0.to_bits(),
            geo.origin.1.to_bits(),
            geo.split.to_bits(),
            geo.d1.to_bits(),
            geo.d2.to_bits(),
            geo.u1_extent.to_bits(),
            geo.u2_extent.to_bits(),
            morph.handle_t.to_bits(),
            morph.tube_collapse.to_bits(),
        ] {
            sample = sample.wrapping_mul(31).wrapping_add(bits);
        }
        (n_trajs, n_pts, sample)
    }

    /// Draw the flat surface, re-rasterizing into the offscreen target only
    /// when the camera or geometry/morph state changed; otherwise blit the
    /// cached texture.  Highlights are drawn uncached on top every frame
    /// (cheap — a handful of short traces), mirroring
    /// `TorusRender::draw_scaled`'s treatment of `draw_torus_highlights_scaled`.
    pub fn draw(
        &mut self,
        trajectories: &[Vec<FlatPhasePoint>],
        highlights: &[Vec<FlatPhasePoint>],
        boundary: &[(FlatPhasePoint, bool)],
        ctx: &FlatDrawContext,
    ) {
        let geom_key = Self::geometry_key(trajectories, boundary, ctx.geo, ctx.morph);
        self.cache
            .draw(ctx.cam.state_key(), ctx.win_w, ctx.win_h, geom_key, || {
                draw_flat(trajectories, ctx.geo, ctx.morph, ctx.cam, ctx.win_w, ctx.win_h);
                draw_flat_boundary(boundary, ctx.geo, ctx.morph, ctx.cam, ctx.win_w, ctx.win_h);
            });
        draw_flat_highlights(highlights, ctx.geo, ctx.morph, ctx.cam, ctx.win_w, ctx.win_h);
    }
}

impl Default for FlatRender {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pseudo::classify::classify_level;
    use crate::pseudo::geometry::level_geometry;
    use crate::pseudo::morph::flat_morph;
    use crate::table::Table;

    fn cf() -> ConfocalParams {
        ConfocalParams::new(4.0, 1.0)
    }

    fn standard_l_table() -> Table {
        Table::from_bounds(&[0.0, 0.4, 0.8], &[1.4, 2.0, 2.6], |l1, l2, _| {
            let tall = l1 <= 0.4 && l2 >= 2.0;
            let base = l1 <= 0.8 && l2 <= 2.0;
            tall || base
        })
    }

    fn pt(u1: f32, u2: f32, s1: i8, s2: i8) -> FlatPhasePoint {
        FlatPhasePoint {
            u1,
            u2,
            sheet: (s1, s2),
            component: 0,
        }
    }

    /// Regression for the "weird straight lines" bug: right next to a band
    /// end (`u2_extent` near zero), even a tiny same-sheet step in raw
    /// `(u1,u2)` must be flagged as a seam — proving the guard in
    /// `draw_flat`/`draw_flat_highlights` is load-bearing there, not just
    /// theoretical.
    #[test]
    fn test_near_death_step_is_a_seam() {
        let tab = standard_l_table();
        let lam = *tab.hyp.last().unwrap() - 1e-4;
        let level = classify_level(lam, &tab, &cf(), 1e-9, 1e-9);
        let geo = level_geometry(&level);
        let morph = flat_morph(&geo, &tab, &cf(), lam);

        let a = pt(0.1256, 0.0000, -1, 1);
        let b = pt(0.1239, 0.0017, -1, 1);
        let (pa, _) = unified_with_normal(a.u1, a.u2, a.sheet.0, a.sheet.1, &geo, &morph);
        let (pb, _) = unified_with_normal(b.u1, b.u2, b.sheet.0, b.sheet.1, &geo, &morph);
        assert!(
            is_seam(&a, &b, pa, pb, &geo),
            "a tiny same-sheet step near the death cap should still flag as a seam"
        );
    }

    /// A healthy mid-band level's same-sheet steps must never flag, so the
    /// guard never breaks a normal trajectory line into dotted fragments.
    #[test]
    fn test_mid_band_step_is_not_a_seam() {
        let tab = standard_l_table();
        let lam = 2.3f32;
        let level = classify_level(lam, &tab, &cf(), 1e-9, 1e-9);
        let geo = level_geometry(&level);
        let morph = flat_morph(&geo, &tab, &cf(), lam);

        let a = pt(0.10, 1.00, 1, 1);
        let b = pt(0.101, 1.005, 1, 1);
        let (pa, _) = unified_with_normal(a.u1, a.u2, a.sheet.0, a.sheet.1, &geo, &morph);
        let (pb, _) = unified_with_normal(b.u1, b.u2, b.sheet.0, b.sheet.1, &geo, &morph);
        assert!(
            !is_seam(&a, &b, pa, pb, &geo),
            "a tiny same-sheet step mid-band should not flag as a seam"
        );
    }

    /// Regression for the "red trajectory disconnected / torus boundary
    /// fractured" bug: on a genus-2 level with a narrow spine (`split` small
    /// relative to `u2_extent`), the main lobe's own angular speed is high,
    /// so normal same-branch per-sample motion can produce a large embedded
    /// jump — that must NOT be treated as a seam (a fixed distance cutoff
    /// wrongly flagged it; the ratio test doesn't, since the raw step is
    /// proportionally just as large).
    #[test]
    fn test_steep_same_branch_motion_is_not_a_seam() {
        let tab = standard_l_table();
        let lam = 1.2f32;
        let level = classify_level(lam, &tab, &cf(), 1e-9, 1e-9);
        let geo = level_geometry(&level);
        assert!(geo.split < 0.3, "fixture should have a narrow spine");
        let morph = flat_morph(&geo, &tab, &cf(), lam);

        // Mirrors a measured same-branch, same-sheet trajectory step at this
        // level: raw Δu ≈ 0.076 in each axis, comfortably inside the spine.
        let a = pt(0.0000, 1.0652, 1, 1);
        let b = pt(0.0756, 1.1408, 1, 1);
        assert!(a.u1 <= geo.split && b.u1 <= geo.split, "fixture should stay on the main lobe");
        let (pa, _) = unified_with_normal(a.u1, a.u2, a.sheet.0, a.sheet.1, &geo, &morph);
        let (pb, _) = unified_with_normal(b.u1, b.u2, b.sheet.0, b.sheet.1, &geo, &morph);
        assert!(
            !is_seam(&a, &b, pa, pb, &geo),
            "steep-but-legitimate same-branch motion should not flag as a seam"
        );
    }

    /// A tiny step that crosses the main-lobe/handle split must NOT flag as a
    /// seam, at any `u2` in the attach band — not just its midpoint. This is
    /// the regression for the fixed embedding bug: `unified_with_normal` now
    /// matches the main lobe's boundary value exactly at the split, for every
    /// `u2`, so a physically tiny step across it produces a physically tiny
    /// embedded jump, same as anywhere else.
    #[test]
    fn test_branch_crossing_is_no_longer_always_a_seam() {
        let tab = standard_l_table();
        let lam = 1.2f32;
        let level = classify_level(lam, &tab, &cf(), 1e-9, 1e-9);
        let geo = level_geometry(&level);
        let morph = flat_morph(&geo, &tab, &cf(), lam);

        let just_inside = geo.split - 1e-4;
        let just_outside = geo.split + 1e-4;

        for frac in [0.0f32, 0.5, 1.0] {
            let u2 = geo.attach.0 + frac * geo.d2;
            let a = pt(just_inside, u2, 1, 1);
            let b = pt(just_outside, u2, 1, 1);
            let (pa, _) = unified_with_normal(a.u1, a.u2, a.sheet.0, a.sheet.1, &geo, &morph);
            let (pb, _) = unified_with_normal(b.u1, b.u2, b.sheet.0, b.sheet.1, &geo, &morph);
            assert!(
                !is_seam(&a, &b, pa, pb, &geo),
                "a tiny step across the fixed split at u2={u2} should not flag as a seam"
            );
        }
    }

    /// Regression for the reported jumps at some λc: an orbit aimed straight
    /// at the reflex corner has no canonical continuation (doc §8, the
    /// 3-pronged singularity), so it must be terminated there instead of
    /// being silently reflected through one of the three prongs — which is
    /// what was producing the discontinuities.
    #[test]
    fn test_trajectory_terminates_near_reflex_corner() {
        let cf = ConfocalParams::standard();
        let domain = crate::domain::confocal_lshape_standard(cf, 0.4, 0.8, 1.4, 2.0, 2.6);
        // Segment 3 is the tall leg's inner ellipse arc, `a: reflex` — see
        // `confocal_lshape_standard`'s doc comment.
        let crate::domain::Segment::Quad { a: corner, .. } = domain.segments[3] else {
            panic!("segment 3 should be the tall-leg arc starting at the reflex corner");
        };
        let normal = domain.segments[3].inward_normal(corner);

        let tab = standard_l_table();
        let lam = 1.2f32; // mid-band genus-2: reflex corner accessible
        let level = classify_level(lam, &tab, &cf, 1e-9, 1e-9);

        let p0 = corner + normal * 0.15;
        let v0 = -normal; // aimed straight back at the corner

        let pts = sample_flat_trajectory(&domain, p0, v0, 30, 20, &level);

        let corners = reflex_corners_in_flat(&level);
        assert!(!corners.is_empty(), "fixture level should have a reflex corner");
        for pt in &pts {
            assert!(
                !near_reflex_corner(pt.u1, pt.u2, &corners),
                "a sample was emitted within CORNER_EPS of the reflex corner: ({}, {})",
                pt.u1,
                pt.u2
            );
        }
        assert!(
            pts.len() < 30 * 20,
            "trajectory aimed straight at the reflex corner should terminate before the step budget"
        );
        assert!(!pts.is_empty(), "should still record the approach up to the corner");
    }
}

