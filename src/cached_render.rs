//! Generic render-to-texture cache shared by the smooth-torus renderer
//! (`torus_render::TorusRender`) and the pseudo-integrable flat/genus-2
//! renderer (`pseudo::view::FlatRender`).
//!
//! Both draw a dense, decimated point cloud onto a 3D surface whose pixels
//! only need to change when the caustic level `λc` (or the camera) moves —
//! every other frame the same image is correct.  Re-rasterizing into an
//! offscreen `RenderTarget` only on camera-move or geometry-change, and
//! otherwise blitting the cached texture, was the fix for the ~150K-point
//! torus orbit-freeze (`docs/known_issues.md` #1).  Factoring the cache out
//! here means the L-shape path gets the identical fix instead of redrawing
//! its full decimated cloud unconditionally every frame.

use macroquad::prelude::*;

/// A cheap fingerprint of the camera state, used to detect when the cached
/// render must be regenerated (orbit/zoom changed).
pub fn camera_moved(a: (f32, f32, f32), b: (f32, f32, f32), eps: f32) -> bool {
    (a.0 - b.0).abs() > eps || (a.1 - b.1).abs() > eps || (a.2 - b.2).abs() > eps
}

/// Camera-move threshold below which a drag is jitter, not a real move.
const CAM_EPS: f32 = 1e-4;

/// Offscreen-texture cache: rasterize only when the camera or geometry
/// fingerprint changed, otherwise blit the cached texture.
pub struct CachedSurfaceRender {
    target: Option<RenderTarget>,
    camera_key: (f32, f32, f32),
    size: (u32, u32),
    geom_key: (usize, usize, u32),
}

impl CachedSurfaceRender {
    pub fn new() -> Self {
        Self {
            target: None,
            camera_key: (f32::NAN, f32::NAN, f32::NAN),
            size: (0, 0),
            geom_key: (usize::MAX, usize::MAX, u32::MAX),
        }
    }

    /// Access the cached offscreen texture (for tests / diagnostics).
    pub fn texture(&self) -> Option<&Texture2D> {
        self.target.as_ref().map(|t| &t.texture)
    }

    /// Draw via the cache.  `geom_key` fingerprints the current geometry
    /// (point/line counts, a content sample, and any scale/morph parameters
    /// that change the drawn shape) so a level change forces re-rasterization
    /// even with an unmoved camera.  `rasterize` draws exactly as it would to
    /// the screen — it runs inside a `Camera2D` already mapped so pixel
    /// coordinates from `OrbitCamera3::project` land correctly on the target.
    pub fn draw(
        &mut self,
        cam_key: (f32, f32, f32),
        win_w: f32,
        win_h: f32,
        geom_key: (usize, usize, u32),
        rasterize: impl FnOnce(),
    ) {
        let size = (win_w as u32, win_h as u32);

        // Recreate the target if the window resized.
        if self.target.is_none() || self.size != size {
            self.target = Some(render_target(size.0, size.1));
            self.size = size;
            self.camera_key = (f32::NAN, f32::NAN, f32::NAN);
            self.geom_key = (usize::MAX, usize::MAX, u32::MAX);
        }

        // Re-rasterize only when the camera actually moved (past a jitter
        // threshold) or the geometry changed.  Sub-pixel drag would otherwise
        // trigger a full re-draw of the decimated point cloud every frame,
        // which is the profiled orbit freeze.
        let cam_changed = camera_moved(cam_key, self.camera_key, CAM_EPS);
        let needs_raster = cam_changed || geom_key != self.geom_key;

        if needs_raster {
            let target = self.target.as_ref().expect("render target exists");
            let (tw, th) = (self.size.0.max(1) as f32, self.size.1.max(1) as f32);
            // `cam.project` emits screen-pixel coordinates (origin top-left, y
            // down), i.e. `p.x ∈ [0, win_w]`, `p.y ∈ [0, win_h]`.  A plain
            // `Camera2D` on a render target uses the identity matrix (NDC), so
            // pixel coordinates would land off-target.  Instead map pixel
            // space onto the target: NDC_x = 2·x/tw − 1, NDC_y = −2·y/th + 1,
            // which the `Camera2D` zoom+offset reproduce exactly (the − on y
            // restores y-down; blitting with `flip_y` presents it upright).
            let cam2d = Camera2D {
                render_target: Some(target.clone()),
                zoom: vec2(2.0 / tw, -2.0 / th),
                offset: vec2(-1.0, 1.0),
                ..Default::default()
            };
            set_camera(&cam2d);
            clear_background(Color::new(0.0, 0.0, 0.0, 0.0));
            rasterize();
            set_default_camera();

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
                    dest_size: Some(vec2(win_w, win_h)),
                    flip_y: true, // render targets are Y-flipped in macroquad
                    ..Default::default()
                },
            );
        }
    }
}

impl Default for CachedSurfaceRender {
    fn default() -> Self {
        Self::new()
    }
}
