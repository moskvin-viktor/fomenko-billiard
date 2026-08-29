//! Cached rendering of the Liouville torus.
//!
//! The profile showed the dominant cost is re-issuing ~150K `draw_circle` calls
//! every frame (≈49% in `buffer_update` + `draw_poly` + `QuadGl::geometry`),
//! even though the torus is static.  The fix is **render-to-texture**: rasterize
//! the torus into an offscreen `RenderTarget` only when the camera (Euler
//! angles) or the torus geometry changes, then blit that texture each frame.
//!
//! The camera is a separate concern (Euler angles in `OrbitCamera3`); here we
//! only decide *when* to re-rasterize and do the cheap blit.

use crate::phase3d::{
    draw_phase_points_scaled, draw_phase_trajectories_scaled, OrbitCamera3, PhasePoint,
};
use macroquad::prelude::*;

/// A cheap fingerprint of the camera state, used to detect when the cached
/// torus render must be regenerated (orbit/zoom changed).
pub fn camera_moved(a: (f32, f32, f32), b: (f32, f32, f32), eps: f32) -> bool {
    (a.0 - b.0).abs() > eps || (a.1 - b.1).abs() > eps || (a.2 - b.2).abs() > eps
}

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

/// A scratch offscreen state structure for the torus point cloud.
pub struct TorusRender {
    /// Offscreen target we rasterize the torus into.
    target: Option<RenderTarget>,
    /// Camera state the current texture was rendered with.
    camera_key: (f32, f32, f32),
    /// Window size the texture was rendered at.
    size: (u32, u32),
    /// Fingerprint of the geometry the texture was rendered with.
    geom_key: (usize, usize, u32),
}

impl TorusRender {
    pub fn new() -> Self {
        Self {
            target: None,
            camera_key: (f32::NAN, f32::NAN, f32::NAN),
            size: (0, 0),
            geom_key: (usize::MAX, usize::MAX, u32::MAX),
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
        self.target.as_ref().map(|t| &t.texture)
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
        let size = (ctx.win_w as u32, ctx.win_h as u32);
        let cam_key = ctx.cam.state_key();
        let geom_key = Self::geometry_key(trajectories, ctx.scale);

        // Recreate the target if the window resized.
        if self.target.is_none() || self.size != size {
            self.target = Some(render_target(size.0, size.1));
            self.size = size;
            self.camera_key = (f32::NAN, f32::NAN, f32::NAN);
            self.geom_key = (usize::MAX, usize::MAX, u32::MAX);
        }

        // Re-rasterize only when the camera actually moved (past a jitter
        // threshold) or the geometry changed.  Sub-pixel drag would otherwise
        // trigger a full ~150K-point re-draw every frame, which is the
        // profiled orbit freeze.
        const CAM_EPS: f32 = 1e-4;
        let cam_changed = crate::torus_render::camera_moved(cam_key, self.camera_key, CAM_EPS);
        let needs_raster = cam_changed || geom_key != self.geom_key;

        if needs_raster {
            self.rasterize(trajectories, boundary, ctx);
            self.camera_key = cam_key;
            self.geom_key = geom_key;
        }

        // Blit the cached texture to the screen.
        if let Some(target) = &self.target {
            draw_texture_ex(
                &target.texture,
                0.0,
                0.0,
                WHITE,
                DrawTextureParams {
                    dest_size: Some(vec2(ctx.win_w, ctx.win_h)),
                    flip_y: true, // render targets are Y-flipped in macroquad
                    ..Default::default()
                },
            );
        }

        // Bold red short trajectories on the torus surface, one per torus,
        // drawn on top of the cached texture each frame (cheap).
        crate::phase3d::draw_torus_highlights_scaled(
            highlights, ctx.cam, ctx.win_w, ctx.win_h, ctx.scale,
        );
    }

    /// Rasterize the torus into the offscreen target.
    fn rasterize(
        &self,
        trajectories: &[Vec<PhasePoint>],
        boundary: &[(PhasePoint, bool)],
        ctx: &DrawContext,
    ) {
        let target = self.target.as_ref().expect("render target exists");
        let size = self.size;
        let (tw, th) = (size.0.max(1) as f32, size.1.max(1) as f32);

        // `cam.project` emits screen-pixel coordinates (origin top-left, y
        // down), i.e. `p.x ∈ [0, win_w]`, `p.y ∈ [0, win_h]`.  A plain
        // `Camera2D` on a render target uses the identity matrix (NDC), so
        // pixel coordinates would land off-target.  Instead map pixel space
        // onto the target: NDC_x = 2·x/tw − 1, NDC_y = −2·y/th + 1, which the
        // `Camera2D` zoom+offset reproduce exactly (the − on y restores y-down;
        // blitting with `flip_y` presents it upright).
        let cam2d = Camera2D {
            render_target: Some(target.clone()),
            zoom: vec2(2.0 / tw, -2.0 / th),
            offset: vec2(-1.0, 1.0),
            ..Default::default()
        };
        set_camera(&cam2d);

        clear_background(Color::new(0.0, 0.0, 0.0, 0.0));

        // The torus is drawn in screen space (projected by the orbit camera),
        // so we reuse the existing point/line drawing into the target.
        draw_phase_points_scaled(trajectories, ctx.cam, ctx.win_w, ctx.win_h, ctx.scale);
        draw_phase_trajectories_scaled(trajectories, ctx.cam, ctx.win_w, ctx.win_h, ctx.scale);
        crate::phase3d::draw_boundary_preimage(boundary, ctx.scale, ctx.cam, ctx.win_w, ctx.win_h);

        set_default_camera();
    }
}
