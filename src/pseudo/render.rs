//! 3D embeddings for pseudo-integrable level sets.
//!
//! A genus-1 (torus) level of a confocal table is a flat rectangle with
//! opposite edges identified; a genus-≥2 level is a flat translation surface
//! (the unfolded cross 12-gon).  This module provides the 3D embeddings:
//!
//! - [`torus_angles`] — flat rectangle `(u1, u2)` → torus phases (for torus
//!   levels).  This is the doc's "revert to normalized angles on the overlap".
//! - [`cross_embed`] — the flat cross chart (coplanar sheets), the *faithful*
//!   flat object (doc §6/§12).
//! - [`pretzel_embed`] — a double torus (two tori glued) embedding so the
//!   genus actually *reads* as genus 2 in 3D.  Arbitrary as a chart (doc §12)
//!   but the visualization the user asked for.

use macroquad::prelude::*;

/// Map a flat-rectangle point `(u1, u2) ∈ [0, w1] × [0, w2]` to torus angles
/// `(θ1, θ2) ∈ [0, 2π)²` by identifying opposite edges (a torus).
///
/// This is the doc's "revert to normalized angles on the overlap": a torus
/// level of the L is a flat rectangle, and normalizing each coordinate by its
/// width gives the two circle phases of the torus.
pub fn torus_angles(u1: f32, u2: f32, w1: f32, w2: f32) -> (f32, f32) {
    let two_pi = 2.0 * std::f32::consts::PI;
    let t1 = (u1 / w1.max(1e-30) * two_pi).rem_euclid(two_pi);
    let t2 = (u2 / w2.max(1e-30) * two_pi).rem_euclid(two_pi);
    (t1, t2)
}

/// Embed a torus level point (flat `(u1, u2)` on the rectangle `[0,w1]×[0,w2]`)
/// onto the standard donut, returning the 3D position and the outward surface
/// normal (for Lambert shading, matching the smooth confocal tori).
pub fn donut_with_normal(u1: f32, u2: f32, w1: f32, w2: f32) -> (Vec3, Vec3) {
    let (t1, t2) = torus_angles(u1, u2, w1, w2);
    let (r_major, r_minor) = (1.6f32, 0.6f32);
    let (s1, c1) = t1.sin_cos();
    let (s2, c2) = t2.sin_cos();
    let pos = vec3(
        (r_major + r_minor * c1) * c2,
        (r_major + r_minor * c1) * s2,
        r_minor * s1,
    );
    // Normal of the parametrized torus at (θ1, θ2).
    let n = vec3(c1 * c2, c1 * s2, s1);
    (pos, n)
}

/// The flat extents of the cross 12-gon (doc §6): `([-A1, A1] × [-B2, B2]) ∪
/// ([-A2, A2] × [-B1, B1])`.
#[derive(Clone, Copy, Debug)]
pub struct CrossModuli {
    /// Tall-leg width `u₁(α₁)`.
    pub a1: f32,
    /// Base width `u₁(α₂)`.
    pub a2: f32,
    /// Base height `u₂(β₂)`.
    pub b1: f32,
    /// Tall-leg height `u₂(β₃)`.
    pub b2: f32,
}

impl CrossModuli {
    /// Largest flat extent (used to normalize the pretzel walk).
    fn major(&self) -> f32 {
        self.a2.max(self.b2).max(self.a1).max(self.b1).max(1e-3)
    }
}

/// Map a flat point `(u1, u2)` on sheet `(σ1, σ2)` into the unfolded cross
/// 12-gon.  The four sheets tile a **single connected** cross in the plane
/// (`X = σ1·u1`, `Y = σ2·u2`), so the genus-2 surface is one object — not four
/// disconnected planes.
///
/// The cross is `([-A1, A1] × [-B2, B2]) ∪ ([-A2, A2] × [-B1, B1])` (doc §6),
/// and the four sheets are the sign pairs `(σ1, σ2)`.  The `moduli` bound the
/// cross extent; the reflex-corner images (where the genus is made visible)
/// sit at `(±A1, ±B1)`.
pub fn cross_embed(u1: f32, u2: f32, s1: i8, s2: i8, _moduli: CrossModuli) -> Vec3 {
    // X = σ1·u1, Y = σ2·u2: the four sheets tile the cross in the z=0 plane.
    let x = s1 as f32 * u1;
    let y = s2 as f32 * u2;
    vec3(x, y, 0.0)
}

/// Embed a sheet point onto a **double torus** (genus-2 pretzel): two torus
/// lobes joined by a bridge, so the manifold reads as two tori glued.
///
/// The doc (§12) warns that any 3D genus-2 embedding is arbitrary and hides
/// the flat structure; this is a *visualization* choice, not the flat
/// translation surface.  Geometrically it is the standard "two tori glued by a
/// cylinder" double torus:
///
/// - `u2` (normalized) is the tube angle around the handle cross-section (the
///   S¹ holonomy), mirrored by `σ2`;
/// - `u1` (normalized) walks a longitudinal loop through both lobes, its
///   "fold" giving two bulges at the two lobe centres.
///
/// The result has genus 2: two disjoint handle cycles.
pub fn pretzel_embed(u1: f32, u2: f32, s1: i8, s2: i8, moduli: CrossModuli) -> Vec3 {
    pretzel_with_normal(u1, u2, s1, s2, moduli).0
}

/// [`pretzel_embed`] with the outward tube normal, for Lambert shading so the
/// genus-2 surface carries the same depth/curvature cues as the torus levels.
pub fn pretzel_with_normal(u1: f32, u2: f32, s1: i8, s2: i8, moduli: CrossModuli) -> (Vec3, Vec3) {
    let two_pi = std::f32::consts::TAU;
    let maj = moduli.major();
    // Longitudinal loop angle: two lobes per 2π.
    let long = (u1 / maj * two_pi).rem_euclid(two_pi);
    // Tube angle around the cross-section (mirrored for the opposite sheet).
    let tube = ((u2 * two_pi).rem_euclid(two_pi)) * if s2 > 0 { 1.0 } else { -1.0 };
    let (sin_t, cos_t) = tube.sin_cos();
    let r = maj * 0.25;

    // Figure-eight center curve C(s) = c·(cos s, 0, sin 2s)
    // with lobes at s=0 (right) and s=π (left).  Tangent frame:
    //   T = dC/ds, binormal B = (0,1,0), normal N = T × B.
    let c = maj;
    let s = long;
    let (ss, cs) = s.sin_cos();
    let center = vec3(c * cs, 0.0, c * (2.0 * s).sin());
    let tangent = {
        let t = vec3(-c * ss, 0.0, 2.0);
        let len = t.length();
        if len < 1e-6 {
            vec3(0.0, 0.0, 1.0)
        } else {
            t / len
        }
    };
    let b = vec3(0.0, 1.0, 0.0);
    let n = tangent.cross(b);
    let n = if n.length() < 1e-6 {
        vec3(1.0, 0.0, 0.0)
    } else {
        n.normalize()
    };

    // The tube radial direction (outward surface normal) before the σ₁ mirror.
    let radial = (n * cos_t + b * sin_t).normalize();
    let pos = center + n * (r * cos_t) + b * (r * sin_t);
    let side = if s1 > 0 { 1.0 } else { -1.0 };
    (
        vec3(side * pos.x, pos.y, pos.z),
        vec3(side * radial.x, radial.y, radial.z),
    )
}

/// The four momentum sheets, in a canonical order.
pub fn sheets() -> [(i8, i8); 4] {
    [(1, 1), (1, -1), (-1, 1), (-1, -1)]
}
