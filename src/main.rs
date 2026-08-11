use billiards::{domain, phase3d, presets, render, TorusRegime, A, B};
use macroquad::prelude::*;
use std::cmp::Ordering;

// The billiard's 2D drawing (Camera, CachedDomain, draw_*) lives in the
// `render` module.  This file is the app: input → state → draw loop.

// ----------------------------------------------------------------
// Second-integral slider range
// ----------------------------------------------------------------

/// Describes the valid range for the second-integral slider.
enum SecondIntegralRange {
    /// Confocal caustic parameter Λ.
    Confocal {
        ell_min: f32,
        ell_max: f32,
        hyp_min: f32,
        #[allow(dead_code)] // only used to size the range, not by the slider
        hyp_max: f32,
        total_range: f32,
    },
    /// Polyline angle θ/π ∈ [-1, 1].
    Angle,
}

impl SecondIntegralRange {
    fn for_domain(domain: &domain::Domain, is_confocal: bool) -> Self {
        if !is_confocal {
            return Self::Angle;
        }
        let (ell_min, ell_max, hyp_min, hyp_max) = LambdaRange::of_domain(domain).slider_bounds();
        let ell_range = (ell_max - ell_min).max(0.0);
        let hyp_range = (hyp_max - hyp_min).max(0.0);
        let total_range = ell_range + hyp_range;
        Self::Confocal {
            ell_min,
            ell_max,
            hyp_min,
            hyp_max,
            total_range,
        }
    }

    /// Map slider fraction t ∈ [0, 1] to the second integral value.
    /// Slider value for a fraction `t ∈ [0, 1]` of the track.
    fn value_at_fraction(&self, t: f32) -> f32 {
        match self {
            Self::Angle => -1.0 + 2.0 * t.clamp(0.0, 1.0),
            Self::Confocal {
                ell_min,
                ell_max,
                hyp_min,
                total_range,
                ..
            } => {
                if *total_range <= 0.0 {
                    return 0.0;
                }
                let ell_range = ell_max - ell_min;
                let t = t.clamp(0.0, 1.0);
                let pos = t * total_range;
                if pos <= ell_range {
                    ell_min + pos
                } else {
                    hyp_min + (pos - ell_range)
                }
            }
        }
    }

    /// Map a second integral value to slider fraction t ∈ [0, 1].
    /// Track fraction `∈ [0, 1]` for a slider value.
    fn fraction_of_value(&self, lam: f32) -> f32 {
        match self {
            Self::Angle => (lam.clamp(-1.0, 1.0) + 1.0) / 2.0,
            Self::Confocal {
                ell_min,
                ell_max,
                hyp_min,
                total_range,
                ..
            } => {
                if *total_range <= 0.0 {
                    return 0.0;
                }
                let ell_range = ell_max - ell_min;
                if lam <= *ell_max {
                    (lam - ell_min) / total_range
                } else if lam >= *hyp_min {
                    (ell_range + (lam - hyp_min)) / total_range
                } else {
                    ell_range / total_range
                }
            }
        }
    }
}

/// Inward margin used to stay clear of the degenerate boundaries at the ends
/// of the second-integral range (Λ = B separatrix, Λ close to A, and the
/// hyperbola walls).  Used by clamping and animation.
const SECOND_INT_EPS: f32 = 0.05;

/// Tighter margin used only by the slider, so the thumb can get closer to the
/// separatrix than the clamp/anim will allow.
const SLIDER_EPS: f32 = 0.02;

/// Λ change per arrow-key press.
const KEY_STEP: f32 = 0.05;

/// Λ swept per second during animation.
const ANIM_SPEED: f32 = 0.4;

/// Slider value change below which we ignore (avoid rebuilds on fp noise).
const SLIDER_JITTER: f32 = 0.0001;

/// The valid Λ range for a confocal domain.
///
/// Encapsulates the boundary-lambda extraction (λ_ell, λ_hyp) so the places
/// that reason about "where the slider / clamping / animation may go" share one
/// extraction instead of re-deriving it.  The inward margin `ε` is a *parameter*:
/// the slider and the clamp/anim legitimately use different values.
struct LambdaRange {
    /// λ_ell — the outer ellipse boundary (min `λ < B`).
    lambda_ell: f32,
    /// λ_hyp — the inner hyperbola boundary (max `λ > B`).
    lambda_hyp: f32,
}

impl LambdaRange {
    /// Extract (λ_ell, λ_hyp) from the domain's quadric boundary arcs.
    fn of_domain(domain: &domain::Domain) -> Self {
        let mut ell = f32::MAX;
        let mut hyp = f32::MIN;
        for seg in &domain.segments {
            if let domain::Segment::Quad { curve, .. } = seg {
                if curve.lambda < B {
                    ell = ell.min(curve.lambda);
                }
                if curve.lambda > B {
                    hyp = hyp.max(curve.lambda);
                }
            }
        }
        // Fall back to the raw confocal range if a boundary was missing.
        Self {
            lambda_ell: if ell == f32::MAX { 0.0 } else { ell },
            lambda_hyp: if hyp == f32::MIN { B + 1.0 } else { hyp },
        }
    }

    /// Values usable for the slider, with the slider's tight margin.
    fn slider_bounds(&self) -> (f32, f32, f32, f32) {
        let e = SLIDER_EPS;
        (
            self.lambda_ell + e, // ell_min
            B - e,               // ell_max
            self.lambda_hyp + e, // hyp_min
            A - e,               // hyp_max
        )
    }

    /// Clamp `lam` to the valid range, jumping over the forbidden gap between
    /// the ellipse and hyperbola sides (the separatrix near `Λ = B`).
    fn clamp(&self, lam: f32, prev: f32) -> f32 {
        let (ell_min, ell_max) = (self.lambda_ell + SECOND_INT_EPS, B - SECOND_INT_EPS);
        let (hyp_min, hyp_max) = (self.lambda_hyp + SECOND_INT_EPS, A - SECOND_INT_EPS);

        let clamped = lam.clamp(ell_min, hyp_max);
        if clamped > ell_max && clamped < hyp_min {
            match lam.partial_cmp(&prev).unwrap_or(Ordering::Equal) {
                Ordering::Greater => hyp_min,
                _ => ell_max,
            }
        } else {
            clamped
        }
    }
}

// ----------------------------------------------------------------
// Slider widget
// ----------------------------------------------------------------
struct Slider {
    x: f32,
    y: f32,
    width: f32,
    dragging: bool,
}

impl Slider {
    fn new() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            width: 0.0,
            dragging: false,
        }
    }

    fn update(&mut self, value: f32, range: &SecondIntegralRange, label: &str) -> Option<f32> {
        let track_h = 4.0;
        let thumb_r = 8.0;
        let cy = self.y + track_h / 2.0;

        // Track background
        draw_rectangle(
            self.x,
            self.y,
            self.width,
            track_h,
            color_u8!(60, 60, 100, 180),
        );

        // Filled portion
        let t = range.fraction_of_value(value);
        let fill_w = t * self.width;
        if fill_w > 0.0 {
            draw_rectangle(
                self.x,
                self.y,
                fill_w,
                track_h,
                color_u8!(130, 130, 200, 220),
            );
        }

        // Thumb
        let thumb_x = self.x + t * self.width;
        draw_circle(thumb_x, cy, thumb_r, color_u8!(220, 220, 255, 255));
        draw_circle_lines(thumb_x, cy, thumb_r, 1.5, color_u8!(100, 100, 160, 200));

        // Gap indicator (confocal-only)
        if let SecondIntegralRange::Confocal {
            ell_max, hyp_min, ..
        } = range
        {
            if ell_max < hyp_min {
                let gap_t = range.fraction_of_value(*ell_max);
                let gap_x = self.x + gap_t * self.width;
                draw_line(
                    gap_x,
                    self.y - 2.0,
                    gap_x,
                    self.y + track_h + 2.0,
                    2.0,
                    render::BG,
                );
            }
        }

        // Label and value text
        let value_text = match range {
            SecondIntegralRange::Angle => format!("{:.3}", value),
            SecondIntegralRange::Confocal { .. } => {
                let caustic_label = if value < B { "ellipse" } else { "hyperbola" };
                format!("{:.3} ({})", value, caustic_label)
            }
        };
        let info = format!("{} = {}", label, value_text);
        draw_text(
            &info,
            self.x - 60.0,
            self.y + track_h + 22.0,
            14.0,
            color_u8!(180, 180, 210, 200),
        );

        // Mouse interaction
        let mx = mouse_position().0;
        let my = mouse_position().1;

        if is_mouse_button_pressed(MouseButton::Left) {
            let dist = ((mx - thumb_x).powi(2) + (my - cy).powi(2)).sqrt();
            if dist <= thumb_r + 4.0
                || (mx >= self.x && mx <= self.x + self.width && (my - cy).abs() <= 12.0)
            {
                self.dragging = true;
            }
        }

        if self.dragging {
            if is_mouse_button_down(MouseButton::Left) {
                let t = ((mx - self.x) / self.width).clamp(0.0, 1.0);
                return Some(range.value_at_fraction(t));
            } else {
                self.dragging = false;
            }
        }

        None
    }
}

// ----------------------------------------------------------------
// Entry point
// ----------------------------------------------------------------
// ----------------------------------------------------------------
// Rebuildable derived-view state
// ----------------------------------------------------------------

/// Everything that must be recomputed together whenever the domain, the second
/// integral, or the 3D toggle changes. Grouping them into one type makes the
/// "rebuild on change" invariant explicit instead of a 3-line copy-paste
/// repeated at every input handler.
struct ViewState {
    /// Start points used to draw the 2D billiard trajectories.
    start_points: Vec<(Vec2, Vec2)>,
    /// 2D billiard trajectories: each is `(from, to)` per bounce.
    trajectories: Vec<Vec<(Vec2, Vec2)>>,
    /// Dense phase-space trajectories filling the 3D Liouville torus.
    /// Empty unless the 3D view is active (or `--burn` forced it).
    phase_trajectories: Vec<Vec<phase3d::PhasePoint>>,
    /// Short (few-bounce) phase-space highlights drawn red on the torus.
    torus_highlights: Vec<Vec<phase3d::PhasePoint>>,
}

impl ViewState {
    fn new(idx: usize, configs: &[presets::Preset], second_int: f32, show_3d: bool) -> Self {
        let mut s = Self {
            start_points: Vec::new(),
            trajectories: Vec::new(),
            phase_trajectories: Vec::new(),
            torus_highlights: Vec::new(),
        };
        s.rebuild(idx, configs, second_int, show_3d);
        s
    }

    /// Recompute the 2D trajectories *and* (if enabled) the 3D torus + red
    /// highlights for the given preset and second-integral value.
    fn rebuild(&mut self, idx: usize, configs: &[presets::Preset], second_int: f32, show_3d: bool) {
        let preset = &configs[idx];
        let domain = &preset.domain;

        // 2D billiard: one trajectory per start point on the caustic.
        let starts = billiards::get_start_points(
            A,
            B,
            second_int,
            domain,
            preset.is_confocal,
            preset.start_center,
        );
        self.start_points = starts.clone();
        self.trajectories = starts
            .iter()
            .map(|&(p, v)| domain.trace(p, v, 300))
            .collect();

        if show_3d {
            // Dense fill of the 3D Liouville torus.  Densely sample the caustic
            // so the union of trajectories sweeps out the full torus surface; a
            // polyline domain has only its single start point.
            let bounds = billiards::torus_bounds(domain);
            let dense_starts = if preset.is_confocal {
                billiards::dense_caustic_starts(A, B, second_int, domain, 24)
            } else {
                billiards::get_start_points(A, B, second_int, domain, false, preset.start_center)
            };
            self.phase_trajectories = dense_starts
                .iter()
                .map(|&(p, v)| phase3d::sample_trajectory_phase_dense(domain, p, v, 200, 8, bounds))
                .collect();

            // Red highlights on the torus: short (4-bounce) phase traces from
            // the *same* start points the 2D view drew, so both views agree on
            // where each torus lives.
            self.torus_highlights = starts
                .iter()
                .map(|&(p, v)| phase3d::sample_trajectory_phase_dense(domain, p, v, 4, 8, bounds))
                .filter(|t| !t.is_empty())
                .collect();
        } else {
            self.phase_trajectories.clear();
            self.torus_highlights.clear();
        }
    }
}

// ----------------------------------------------------------------
// Entry point
// ----------------------------------------------------------------
#[macroquad::main("Mathematical Billiards")]
async fn main() {
    let configs = presets::all_presets();
    let caches: Vec<render::CachedDomain> = configs
        .iter()
        .map(|p| render::CachedDomain::new(&p.domain))
        .collect();

    let mut idx = 0usize;
    let mut second_int: f32 = 0.2; // Λ for confocal, θ/π for polyline
    let mut slider = Slider::new();
    let mut show_3d = false;
    let mut animate = false;
    let mut anim_dir: f32 = 1.0;
    let mut cam3d = phase3d::OrbitCamera3::new();
    let mut view = ViewState::new(idx, &configs, second_int, show_3d);
    let mut torus_render = billiards::torus_render::TorusRender::new();

    // Profiling hook: `--burn N` forces 3D mode, builds the torus, runs N
    // frames, then exits — so `perf`/flamegraph get a bounded run.
    let mut burn_frames: Option<u32> = None;
    for arg in std::env::args().skip(1) {
        if let Some(n) = arg.strip_prefix("--burn=") {
            burn_frames = n.parse().ok();
        } else if arg == "--burn" {
            burn_frames = Some(300);
        }
    }
    if burn_frames.is_some() {
        show_3d = true;
        view.rebuild(idx, &configs, second_int, show_3d);
    }

    loop {
        let (w, h) = (screen_width(), screen_height());
        let cache = &caches[idx];
        let cam = &cache.camera;
        let preset = &configs[idx];

        // Position slider
        slider.x = w * 0.15;
        slider.y = h - 50.0;
        slider.width = w * 0.7;

        let slider_range = SecondIntegralRange::for_domain(&preset.domain, preset.is_confocal);

        // Toggle animation
        if is_key_pressed(KeyCode::A) {
            animate = !animate;
        }

        // Toggle 3D
        if is_key_pressed(KeyCode::P) {
            show_3d = !show_3d;
            view.rebuild(idx, &configs, second_int, show_3d);
        }

        // Tab / Space: next preset
        if is_key_pressed(KeyCode::Tab) || is_key_pressed(KeyCode::Space) {
            idx = (idx + 1) % configs.len();
            second_int = 0.2;
            view.rebuild(idx, &configs, second_int, show_3d);
        }

        // Keyboard: adjust second integral (↑/→ step, ↓/← step back).
        let step = if is_key_pressed(KeyCode::Up) || is_key_pressed(KeyCode::Right) {
            KEY_STEP
        } else if is_key_pressed(KeyCode::Down) || is_key_pressed(KeyCode::Left) {
            -KEY_STEP
        } else {
            0.0
        };
        if step != 0.0 {
            second_int = if preset.is_confocal {
                LambdaRange::of_domain(&preset.domain).clamp(second_int + step, second_int)
            } else {
                (second_int + step).clamp(-1.0, 1.0)
            };
            view.rebuild(idx, &configs, second_int, show_3d);
        }

        // Slider input
        if let Some(new_val) =
            slider.update(second_int, &slider_range, preset.second_integral_label)
        {
            if (new_val - second_int).abs() > SLIDER_JITTER {
                second_int = new_val;
                view.rebuild(idx, &configs, second_int, show_3d);
            }
        }

        // Smooth animation: sweep the second integral back and forth.
        if animate {
            let speed = ANIM_SPEED; // units per second
            let dt = get_frame_time();
            let step = speed * dt;

            if configs[idx].is_confocal {
                let le = SECOND_INT_EPS;
                let range = LambdaRange::of_domain(&configs[idx].domain);
                let ell_min = range.lambda_ell + le;
                let ell_max = B - le;
                let hyp_min = range.lambda_hyp + le;
                let hyp_max = A - le;

                let mut next = second_int + anim_dir * step;
                // Bounce off the ends of the valid range.
                if next >= ell_max && next <= hyp_min {
                    // Crossing the gap: jump to the other side.
                    if anim_dir > 0.0 {
                        next = hyp_min;
                    } else {
                        next = ell_max;
                    }
                }
                if next >= hyp_max {
                    next = hyp_max;
                    anim_dir = -1.0;
                }
                if next <= ell_min {
                    next = ell_min;
                    anim_dir = 1.0;
                }
                second_int = next;
            } else {
                let mut next = second_int + anim_dir * step;
                if next >= 1.0 {
                    next = 1.0;
                    anim_dir = -1.0;
                }
                if next <= -1.0 {
                    next = -1.0;
                    anim_dir = 1.0;
                }
                second_int = next;
            }

            view.rebuild(idx, &configs, second_int, show_3d);
        }

        clear_background(render::BG);

        if show_3d {
            // ---- 3D phase space view ----
            cam3d.handle_input();
            phase3d::draw_axes(&cam3d, w, h, cache.domain_extent);
            // Dense points fill the 2D Liouville torus surface.  Rasterized
            // into an offscreen target and blitted; re-rasterized only when the
            // camera or geometry changes.
            torus_render.draw(
                &view.phase_trajectories,
                &view.torus_highlights,
                &cam3d,
                w,
                h,
            );

            let info = format!(
                "{}  |  {} = {:.3}  |  {} trajs  |  [P] 2D  |  [A] anim {}  |  right-drag orbit",
                preset.label,
                preset.second_integral_label,
                second_int,
                view.phase_trajectories.len(),
                if animate { "on" } else { "off" },
            );
            draw_text(&info, 12.0, 28.0, 18.0, color_u8!(200, 200, 220, 220));
            draw_text(
                "3D phase space: (x, y, θ/π)  —  Liouville torus at fixed second integral",
                12.0,
                h - 12.0,
                14.0,
                color_u8!(150, 150, 170, 140),
            );
        } else {
            // ---- 2D billiard view ----
            render::draw_boundary(cache, w, h);

            if preset.is_confocal {
                render::draw_caustic(A, B, second_int, &preset.domain, cam, w, h);
            } else {
                // Draw the starting direction arrow
                let angle = second_int * std::f32::consts::PI;
                render::draw_start_arrow(preset.start_center, angle, cam, w, h);
            }

            for traj in &view.trajectories {
                render::draw_trajectory(traj, cam, w, h);
            }
            // Highlight one short (few-bounce) red trajectory per torus so one
            // can see exactly where each torus's orbit lives on the billiard.
            let regime = TorusRegime::from_value(
                preset.is_confocal,
                second_int,
                billiards::torus_bounds(&preset.domain).1,
            );
            let picks = render::one_per_torus(&view.start_points, regime);
            for &i in &picks {
                render::draw_trajectory_red(&view.trajectories[i], 4, cam, w, h);
            }
            for &(p, _) in &view.start_points {
                render::draw_start_marker(p, cam, w, h);
            }

            let total: usize = view.trajectories.iter().map(|t| t.len()).sum();
            let info =
                format!(
                "{}  |  {} = {:.3}  |  {} trajs, {} bounces  |  4 bounce/torus red  |  ↑↓  |  [A] anim {}  |  [P] 3D  |  Tab next",
                preset.label, preset.second_integral_label, second_int, view.start_points.len(), total,
                if animate { "on" } else { "off" },
            );
            draw_text(&info, 12.0, 28.0, 18.0, color_u8!(200, 200, 220, 220));

            if preset.is_confocal {
                draw_text(
                    "H = ½|v|²  |  Λ = vx²/a + vy²/b − (x·vy − y·vx)²/(ab)  |  v ⟂ ∇Q_Λ",
                    12.0,
                    h - 12.0,
                    14.0,
                    color_u8!(150, 150, 170, 140),
                );
            } else {
                draw_text(
                    "H = ½|v|²  |  second integral = velocity angle θ = atan2(vy, vx)",
                    12.0,
                    h - 12.0,
                    14.0,
                    color_u8!(150, 150, 170, 140),
                );
            }
        }

        next_frame().await;

        if let Some(n) = burn_frames {
            if n <= 1 {
                break;
            }
            burn_frames = Some(n - 1);
        }
    }
}
