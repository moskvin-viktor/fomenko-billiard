//! Guarantee that the torus renderer actually draws something.
//!
//! Regression test: render-to-texture must produce non-background pixels, not
//! an empty frame.  Uses `#[macroquad::test]` which spawns a real window and
//! runs the async body inside a live macroquad context.
//!
//! Note: we read back the *render-target texture* directly (not the presented
//! screen via `get_screen_data`), which is deterministic under the test
//! harness regardless of when the frame is presented.

use billiards::{phase3d, presets, torus_render::TorusRender};
use macroquad::prelude::*;

/// Count non-background (non-transparent) pixels in an image.
fn count_drawn_pixels(img: &Image) -> u64 {
    let bytes = &img.bytes;
    let w = img.width as usize;
    let h = img.height as usize;
    let row = w * 4;
    let mut vis = 0u64;
    for y in 0..h {
        for x in 0..w {
            let i = y * row + x * 4;
            let (r, g, b) = (bytes[i] as i32, bytes[i + 1] as i32, bytes[i + 2] as i32);
            if r > 20 || g > 20 || b > 20 {
                vis += 1;
            }
        }
    }
    vis
}

/// Build the torus exactly as the app does, render it into the offscreen
/// target, and confirm the cached texture actually contains drawn pixels.
#[macroquad::test]
async fn torus_render_draws_something() {
    // Pick a confocal preset and a valid second integral.
    let preset = presets::all_presets()
        .into_iter()
        .find(|p| p.is_confocal)
        .expect("a confocal preset");
    let domain = &preset.domain;
    let lam = 0.5; // elliptic caustic, well inside (0, B)
    let bounds = billiards::torus_bounds(domain);
    let starts = billiards::dense_caustic_starts(domain, lam, 24);
    assert!(!starts.is_empty(), "should have caustic start points");

    let trajectories: Vec<Vec<phase3d::PhasePoint>> = starts
        .iter()
        .map(|&(p, v)| phase3d::sample_trajectory_phase_dense(domain, p, v, 200, 8, bounds))
        .collect();
    let n_pts: usize = trajectories.iter().map(|t| t.len()).sum();
    assert!(n_pts > 1000, "too few torus points ({})", n_pts);

    let (w, h) = (screen_width(), screen_height());
    let cam = phase3d::OrbitCamera3::new();

    let mut render = TorusRender::new();
    render.draw(&trajectories, &[], &cam, w, h);

    let texture = render.texture().expect("offscreen target texture exists");
    let img = texture.get_texture_data();
    let vis = count_drawn_pixels(&img);
    eprintln!("torus render: {} drawn pixels", vis);
    assert!(
        vis > 0,
        "no drawn pixels in the torus render target — nothing was rasterized"
    );
}
