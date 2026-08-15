//! Tests for the 3D orbit-performance fixes.
//!
//! The user reports "moving in 3D space is terribly slow — it freezes the
//! app."  The profiler (`torus_render.rs` header) shows the dominant cost is
//! re-issuing ~150K `draw_circle` calls whenever the orbit camera changes,
//! because the torus re-rasterizes every frame.  These tests pin the two pure,
//! headless-testable helpers that bound that cost:
//!
//! 1. [`phase3d::decimation_step`] — a stride so a rasterization never draws
//!    more than a fixed budget of points, while keeping coverage when the cloud
//!    is small.
//! 2. [`torus_render::camera_moved`] — a predicate so we re-rasterize only when
//!    the camera actually moved beyond a jitter threshold.

use billiards::{phase3d, torus_render};

// ---------------------------------------------------------------------------
// decimation_step
// ---------------------------------------------------------------------------

#[test]
fn decimation_keeps_everything_when_under_budget() {
    assert_eq!(phase3d::decimation_step(0, 100), 1);
    assert_eq!(phase3d::decimation_step(50, 100), 1);
    assert_eq!(phase3d::decimation_step(100, 100), 1);
}

#[test]
fn decimation_caps_point_draw_at_budget() {
    // For n=150_000, budget=20_000, the drawn count must be <= 20_000
    // (decimation_step picks a stride, then the renderer draws 0, step, 2·step…).
    let n = 150_000usize;
    let budget = 20_000usize;
    let step = phase3d::decimation_step(n, budget);
    assert!(step >= 2, "a cloud over budget must be decimated");
    let drawn = n / step + 1; // indices 0, step, 2·step, …
    assert!(
        drawn <= budget,
        "drawn {} should be <= budget {} (step {})",
        drawn,
        budget,
        step
    );
    // And it shouldn't decimate more than necessary.
    assert!(
        step <= n / budget + 1,
        "step {} too conservative for n/budget {}",
        step,
        n / budget
    );
}

#[test]
fn decimation_preserves_first_point_and_covers_range() {
    let n = 4096usize;
    let budget = 1024usize;
    let step = phase3d::decimation_step(n, budget);
    // Always starts at index 0.
    assert_eq!(0 % step, 0);
    // Coverage spread: the last kept index should not be the very first.
    let last = (n - 1) / step * step;
    assert!(last > 0, "should keep points beyond index 0");
}

// ---------------------------------------------------------------------------
// camera_moved
// ---------------------------------------------------------------------------

#[test]
fn camera_moved_detects_meaningful_motion() {
    use torus_render::camera_moved;
    // Identical -> not moved.
    assert!(!camera_moved((0.1, 0.2, 3.0), (0.1, 0.2, 3.0), 1e-4));
    // Sub-threshold jitter -> not moved.
    assert!(!camera_moved((0.1, 0.2, 3.0), (0.10004, 0.2, 3.0), 1e-3));
    // Real orbit -> moved.
    assert!(camera_moved((0.1, 0.2, 3.0), (0.3, 0.2, 3.0), 1e-3));
    assert!(camera_moved((0.1, 0.2, 3.0), (0.1, 0.2, 3.2), 1e-3));
}
