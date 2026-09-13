//! Cached rendering of the Liouville torus.
//!
//! The profile showed the dominant cost is re-issuing ~150K `draw_circle` calls
//! every frame (≈49% in `buffer_update` + `draw_poly` + `QuadGl::geometry`),
//! even though the torus is static.  The fix is **render-to-texture**: rasterize
//! the torus into an offscreen `RenderTarget` only when the camera (Euler
//! angles) or the torus geometry changes, then blit that texture each frame.
//! The cache mechanics live in [`crate::cached_render::CachedSurfaceRender`]
//! (shared with the L-shape's `pseudo::view::FlatRender`); this module only
//! supplies the torus-specific geometry fingerprint and the rasterize closure.
//!
//! The camera is a separate concern (Euler angles in `OrbitCamera3`); here we
//! only decide *when* to re-rasterize and do the cheap blit.

use crate::cached_render::CachedSurfaceRender;
use crate::phase3d::{
    draw_phase_points_scaled, draw_phase_trajectories_scaled, OrbitCamera3, PhasePoint,
};
use macroquad::prelude::*;

/// A cheap fingerprint of the camera state, used to detect when a cached
/// render must be regenerated (orbit/zoom changed).
pub use crate::cached_render::camera_moved;

/// Radial extent assigned to the torus renderer, derived from the size of the
/// billiard domain.  Each torus lobe is centred on the origin and scaled to
/// fill (roughly) the projected domain, so when the domain grows or shrinks the
/// torus tracks it instead of staying a fixed absolute size (as it did before).
pub struct DomainScale {
    pub r_major: f32,
    pub r_minor: f32,
    /// Centre-to-centre distance between torus lobe `i` and lobe `i+1`.
    pub gap: f32,
}

impl DomainScale {
    /// From a domain that filled the 2D view with size `extent` (half-max
    /// coordinate magnitude), produce torus radii that fill a similar fraction
    /// of the 3D view.  The reference values match the torus drawn before this
    /// feature when `extent ≈ 2.3`, so existing renders change little.
    pub fn of_extent(extent: f32) -> Self {
        let r_major = 0.7 * extent;
        Self {
            r_major,
            r_minor: 0.26 * extent,
            gap: 3.0 * r_major,
        }
    }

    /// Smoothly deform the torus geometry toward the degenerate limit of the
    /// level λ, so sweeping λ makes the tori *continuously* collapse onto the
    /// 1D curves of the critical layers instead of jumping:
    ///
    /// The torus map's convention is that θ₁ is *always* the collapsing
    /// libration (see `torus::map::hyperbolic`), so every degeneration renders
    /// as the tube radius shrinking — the torus thins onto its equator ring,
    /// the hole never closes:
    ///
    /// - λ → λ_ell (ellipse wall): the wall–caustic libration collapses,
    ///   `r_minor → 0` — each torus thins onto its equator ring of radius
    ///   `r_major`.
    /// - λ → b (separatrix, from either side): `r_minor → 0` as above, and on
    ///   the two-tori side the lobes slide together, `gap → 2·r_major`, so at
    ///   λ = b the two equator rings touch at one point — the figure-eight is
    ///   the exact geometric limit.
    /// - λ → a (focal axis): the caustic–`a` libration collapses, `r_minor →
    ///   0` — the single torus thins onto its equator ring, which is the
    ///   focal-axis orbit.
    ///
    /// Each factor is a smoothstep of the distance to the critical value,
    /// normalised by a fraction of the containing band, so mid-band tori are
    /// untouched and every factor reaches 0 exactly at the critical level.
    pub fn morph_to_level(mut self, structure: &crate::confocal::ConfocalStructure, lam: f32) -> Self {
        /// Fraction of the band width over which the collapse happens.
        const WINDOW: f32 = 0.35;
        fn smoothstep(x: f32) -> f32 {
            let x = x.clamp(0.0, 1.0);
            x * x * (3.0 - 2.0 * x)
        }
        // Collapse factor for distance `d` to a critical value inside a band
        // of width `w`: 0 at the critical value, 1 past the morph window.
        let factor = |d: f32, w: f32| -> f32 {
            if w <= 0.0 {
                return 1.0;
            }
            smoothstep(d / (WINDOW * w))
        };

        let (a, b) = (structure.cf.a, structure.cf.b);
        let w_ell = b - structure.lambda_ell; // elliptic band width
        let w_hyp = a - b; // hyperbolic band width

        if lam <= b {
            // Elliptic side: collapse at the ellipse wall and at the separatrix.
            let s = factor(lam - structure.lambda_ell, w_ell) * factor(b - lam, w_ell);
            self.r_minor *= s;
            // Two lobes slide together approaching the separatrix.
            let s_sep = factor(b - lam, w_ell);
            self.gap = self.r_major * (2.0 + s_sep);
        } else {
            // Hyperbolic side: the θ₁ libration collapses both at the
            // separatrix and at the focal axis — the torus thins onto its
            // equator ring; the hole (r_major) is preserved.
            self.r_minor *= factor(lam - b, w_hyp) * factor(a - lam, w_hyp);
        }
        self
    }
}

/// Everything the renderer needs apart from its own cached state: the camera,
/// window size and the domain-derived torus scale.  Bundled so the draw call
/// stays small.
pub struct DrawContext<'a> {
    pub cam: &'a OrbitCamera3,
    pub win_w: f32,
    pub win_h: f32,
    pub scale: &'a DomainScale,
}

impl Default for DomainScale {
    fn default() -> Self {
        Self::of_extent(2.3)
    }
}

/// A scratch offscreen state structure for the torus point cloud, wrapping
/// the shared [`CachedSurfaceRender`] with torus-specific fingerprinting.
pub struct TorusRender {
    cache: CachedSurfaceRender,
}

impl TorusRender {
    pub fn new() -> Self {
        Self {
            cache: CachedSurfaceRender::new(),
        }
    }

    /// A cheap fingerprint of the torus geometry: number of trajectories, total
    /// point count, a sample angle, and the torus scale (radii/gap morph with
    /// λ, so a scale change must re-rasterize even when the point cloud is
    /// unchanged).  O(number of trajectories), not O(points).
    fn geometry_key(trajectories: &[Vec<PhasePoint>], scale: &DomainScale) -> (usize, usize, u32) {
        let n_trajs = trajectories.len();
        let n_pts: usize = trajectories.iter().map(|t| t.len()).sum();
        // Sample a few angles to catch same-count-but-different-geometry cases.
        let mut sample = 0u32;
        for t in trajectories.iter().take(4) {
            if let Some(p) = t.first() {
                sample = sample
                    .wrapping_mul(31)
                    .wrapping_add(p.theta1.to_bits() ^ p.theta2.to_bits());
            }
        }
        for bits in [
            scale.r_major.to_bits(),
            scale.r_minor.to_bits(),
            scale.gap.to_bits(),
        ] {
            sample = sample.wrapping_mul(31).wrapping_add(bits);
        }
        (n_trajs, n_pts, sample)
    }

    /// Access the cached offscreen texture (for tests / diagnostics).
    pub fn texture(&self) -> Option<&Texture2D> {
        self.cache.texture()
    }

    /// Draw the torus, re-rasterizing into the offscreen target only when the
    /// camera or geometry changed; otherwise blit the cached texture.
    pub fn draw(
        &mut self,
        trajectories: &[Vec<PhasePoint>],
        highlights: &[Vec<PhasePoint>],
        cam: &OrbitCamera3,
        win_w: f32,
        win_h: f32,
    ) {
        let no_boundary = Vec::<(PhasePoint, bool)>::new();
        let ctx = DrawContext {
            cam,
            win_w,
            win_h,
            scale: &DomainScale::default(),
        };
        self.draw_scaled(trajectories, highlights, &no_boundary, &ctx);
    }

    /// Draw the torus with a domain-derived scale (so it grows/shrinks with the
    /// billiard), and overlay the π⁻¹(boundary) curves on the surface.  See
    /// [`DomainScale`].
    pub fn draw_scaled(
        &mut self,
        trajectories: &[Vec<PhasePoint>],
        highlights: &[Vec<PhasePoint>],
        boundary: &[(PhasePoint, bool)],
        ctx: &DrawContext,
    ) {
        let geom_key = Self::geometry_key(trajectories, ctx.scale);
        self.cache.draw(
            ctx.cam.state_key(),
            ctx.win_w,
            ctx.win_h,
            geom_key,
            || {
                // The torus is drawn in screen space (projected by the orbit
                // camera), so we reuse the existing point/line drawing into
                // the target.
                draw_phase_points_scaled(trajectories, ctx.cam, ctx.win_w, ctx.win_h, ctx.scale);
                draw_phase_trajectories_scaled(
                    trajectories,
                    ctx.cam,
                    ctx.win_w,
                    ctx.win_h,
                    ctx.scale,
                );
                crate::phase3d::draw_boundary_preimage(
                    boundary, ctx.scale, ctx.cam, ctx.win_w, ctx.win_h,
                );
            },
        );

        // Bold red short trajectories on the torus surface, one per torus,
        // drawn on top of the cached texture each frame (cheap).
        crate::phase3d::draw_torus_highlights_scaled(
            highlights, ctx.cam, ctx.win_w, ctx.win_h, ctx.scale,
        );
    }
}

impl Default for TorusRender {
    fn default() -> Self {
        Self::new()
    }
}
