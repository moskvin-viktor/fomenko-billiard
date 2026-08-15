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

use crate::phase3d::{draw_phase_points, draw_phase_trajectories, OrbitCamera3, PhasePoint};
use macroquad::prelude::*;

/// A cheap fingerprint of the camera state, used to detect when the cached
/// torus render must be regenerated (orbit/zoom changed).
pub fn camera_moved(a: (f32, f32, f32), b: (f32, f32, f32), eps: f32) -> bool {
    (a.0 - b.0).abs() > eps || (a.1 - b.1).abs() > eps || (a.2 - b.2).abs() > eps
}

/// Cached offscreen render of the torus point cloud.
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
    /// point count, and a sample angle.  O(number of trajectories), not O(points).
    fn geometry_key(trajectories: &[Vec<PhasePoint>]) -> (usize, usize, u32) {
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
        let size = (win_w as u32, win_h as u32);
        let cam_key = cam.state_key();
        let geom_key = Self::geometry_key(trajectories);

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
            self.rasterize(trajectories, cam, win_w, win_h);
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

        // Bold red short trajectories on the torus surface, one per torus,
        // drawn on top of the cached texture each frame (cheap).
        crate::phase3d::draw_torus_highlights(highlights, cam, win_w, win_h);
    }

    /// Rasterize the torus into the offscreen target.
    fn rasterize(
        &self,
        trajectories: &[Vec<PhasePoint>],
        cam: &OrbitCamera3,
        win_w: f32,
        win_h: f32,
    ) {
        let target = self.target.as_ref().expect("render target exists");
        let size = self.size;
        let (tw, th) = (size.0.max(1) as f32, size.1.max(1) as f32);

        // `cam.project` emits screen-pixel coordinates (origin top-left, y
        // down), i.e. `p.x ∈ [0, win_w]`, `p.y ∈ [0, win_h]`.  A plain
        // `Camera2D` on a render target uses the identity matrix (NDC), so
        // pixel coordinates would land off-target.  Instead map pixel space
        // onto the target: NDC_x = 2·x/tw − 1, NDC_y = −2·y/th + 1, which the
        // Camera2D zoom+offset reproduce exactly (the − on y restores y-down;
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
        draw_phase_points(trajectories, cam, win_w, win_h);
        draw_phase_trajectories(trajectories, cam, win_w, win_h);

        set_default_camera();
    }
}
