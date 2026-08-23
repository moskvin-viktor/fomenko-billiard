//! The app: input → state → draw loop.
//!
//! Extracted from the old `main.rs` monolith.  `main.rs` is now a thin shell
//! that constructs an [`App`] and runs it; the loop state, the rebuild-on-change
//! [`ViewState`], and the input/animation/draw logic all live here.

use crate::second_integral::{LambdaRange, SecondIntegralRange};
use crate::torus::ConfocalParams;
use crate::{
    phase3d, presets, render,
    ui::{self, Slider},
};
use macroquad::prelude::*;

// Λ change per arrow-key press.
const KEY_STEP: f32 = 0.05;

// Λ swept per second during animation.
const ANIM_SPEED: f32 = 0.4;

// Rebuildable derived-view state
// ----------------------------------------------------------------
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
    /// π⁻¹(boundary): phase-space points on the domain walls / caustic,
    /// drawn in distinct colours on the torus so the billiard's own boundary is
    /// visible wrapped around the abstract torus.
    boundary_preimage: Vec<(phase3d::PhasePoint, bool)>,
    /// Flat-chart trajectories for a pseudo-integrable table (the L).  When
    /// non-empty, the 3D view renders these instead of the torus path.
    flat_trajectories: Vec<Vec<crate::pseudo::FlatPhasePoint>>,
    /// Short (few-bounce) example trajectories on the flat surface, drawn red
    /// (like the smooth case's `torus_highlights`).
    flat_highlights: Vec<Vec<crate::pseudo::FlatPhasePoint>>,
    /// π⁻¹(boundary) on the flat map: walls + caustic, drawn cyan/orange.
    flat_boundary: Vec<(crate::pseudo::FlatPhasePoint, bool)>,
    /// The classified level of the current flat view (torus vs genus-2).
    flat_level: Option<crate::pseudo::Level>,
    /// At a critical caustic level, the collapsed 2D phase manifold (a thin
    /// flat sheet) instead of a 3D torus.  None on regular levels.
    critical_sheet: Option<phase3d::CriticalSheet>,
}

impl ViewState {
    /// An empty view, before the first rebuild against real configs.
    fn empty() -> Self {
        Self {
            start_points: Vec::new(),
            trajectories: Vec::new(),
            phase_trajectories: Vec::new(),
            torus_highlights: Vec::new(),
            boundary_preimage: Vec::new(),
            flat_trajectories: Vec::new(),
            flat_highlights: Vec::new(),
            flat_boundary: Vec::new(),
            flat_level: None,
            critical_sheet: None,
        }
    }

    fn new(idx: usize, configs: &[presets::Preset], second_int: f32, show_3d: bool) -> Self {
        let mut s = Self::empty();
        s.rebuild(idx, configs, second_int, show_3d);
        s
    }

    /// Recompute the 2D trajectories *and* (if enabled) the 3D torus + red
    /// highlights for the given preset and second-integral value.
    fn rebuild(&mut self, idx: usize, configs: &[presets::Preset], second_int: f32, show_3d: bool) {
        let preset = &configs[idx];
        let domain = &preset.domain;

        // 2D billiard: one trajectory per start point on the caustic.  On a
        // degenerate critical layer the caustic collapses to a border piece and
        // the ball slides along it; use the dedicated critical sampler so the
        // 2D view still has start points (the phase sheet draws them).
        let start_points = if preset.is_confocal {
            if let Some(structure) = crate::confocal::ConfocalStructure::of_domain(domain) {
                if let Some(crit) =
                    crate::confocal::critical_caustic_starts(&structure, domain, second_int, 64)
                {
                    crit
                } else {
                    crate::get_start_points(
                        second_int,
                        domain,
                        preset.is_confocal,
                        preset.start_center,
                    )
                }
            } else {
                crate::get_start_points(second_int, domain, preset.is_confocal, preset.start_center)
            }
        } else {
            crate::get_start_points(second_int, domain, preset.is_confocal, preset.start_center)
        };
        self.start_points = start_points.clone();
        self.trajectories = start_points
            .iter()
            .map(|&(p, v)| domain.trace(p, v, 300))
            .collect();

        if show_3d {
            // On a critical caustic layer the phase manifold collapses to a
            // thin 2D sheet; show it instead of a 3D torus.
            let cf = crate::torus::ConfocalParams::standard();
            if let Some(sheet) = phase3d::critical_phase_sheet(domain, second_int, 64) {
                self.critical_sheet = Some(sheet);
                // Keep the 2D start points/trajectories; clear the torus path.
                self.phase_trajectories.clear();
                self.torus_highlights.clear();
                self.boundary_preimage.clear();
                self.flat_trajectories.clear();
                self.flat_highlights.clear();
                self.flat_boundary.clear();
                self.flat_level = None;
                return;
            }
            self.critical_sheet = None;

            // Pseudo-integrable table (the L): classify the level and sample
            // through the flat chart, so torus levels render as a torus and
            // genus-2 levels as the unfolded cross.
            let table = crate::table::Table::from_domain(domain, &cf);
            if let Some(tab) = table {
                let level = crate::pseudo::classify_level(second_int, &tab, &cf, 1e-9, 1e-9);
                match &level {
                    crate::pseudo::Level::Torus { .. }
                    | crate::pseudo::Level::GenusSurface { .. } => {
                        let dense_starts = crate::dense_caustic_starts(domain, second_int, 24);
                        self.flat_trajectories = dense_starts
                            .iter()
                            .map(|&(p, v)| {
                                crate::pseudo::sample_flat_trajectory(domain, p, v, 200, 8, &level)
                            })
                            .collect();
                        // Example trajectories: short (few-bounce) flat traces
                        // from the *same* 2D start points, drawn red on top.
                        self.flat_highlights = start_points
                            .iter()
                            .map(|&(p, v)| {
                                crate::pseudo::sample_flat_trajectory(domain, p, v, 4, 8, &level)
                            })
                            .filter(|t| !t.is_empty())
                            .collect();
                        // π⁻¹(boundary): walls + caustic on the flat map.
                        self.flat_boundary =
                            crate::pseudo::sample_flat_boundary(domain, second_int, &level);
                        self.flat_level = Some(level);
                        // Clear the torus path for this view.
                        self.phase_trajectories.clear();
                        self.torus_highlights.clear();
                        self.boundary_preimage.clear();
                        return;
                    }
                    _ => {}
                }
            }

            // Dense fill of the 3D Liouville torus.  Densely sample the caustic
            // so the union of trajectories sweeps out the full torus surface; a
            // polyline domain has only its single start point.
            let bounds = crate::confocal::ConfocalStructure::of_domain(domain)
                .map(|s| s.torus_bounds())
                .unwrap_or((0.0, None));
            let dense_starts = if preset.is_confocal {
                crate::dense_caustic_starts(domain, second_int, 24)
            } else {
                crate::get_start_points(second_int, domain, false, preset.start_center)
            };
            self.phase_trajectories = dense_starts
                .iter()
                .map(|&(p, v)| phase3d::sample_trajectory_phase_dense(domain, p, v, 200, 8, bounds))
                .collect();

            // Red highlights on the torus: short (4-bounce) phase traces from
            // the *same* start points the 2D view drew, so both views agree on
            // where each torus lives.
            self.torus_highlights = start_points
                .iter()
                .map(|&(p, v)| phase3d::sample_trajectory_phase_dense(domain, p, v, 4, 8, bounds))
                .filter(|t| !t.is_empty())
                .collect();

            // π⁻¹(boundary): the walls + caustic curves of the billiard,
            // mapped onto the torus in distinct colours.
            self.boundary_preimage = phase3d::boundary_preimage_points(domain, second_int, bounds);
        } else {
            self.phase_trajectories.clear();
            self.torus_highlights.clear();
            self.boundary_preimage.clear();
            self.flat_trajectories.clear();
            self.flat_highlights.clear();
            self.flat_boundary.clear();
            self.flat_level = None;
            self.critical_sheet = None;
        }
    }
}

// ----------------------------------------------------------------
// App
// ----------------------------------------------------------------

/// The whole application: owns the loop state and wires input → state → draw.
pub struct App {
    configs: Vec<presets::Preset>,
    caches: Vec<render::CachedDomain>,
    idx: usize,
    second_int: f32,
    slider: Slider,
    show_3d: bool,
    animate: bool,
    anim_dir: f32,
    cam3d: phase3d::OrbitCamera3,
    view: ViewState,
    torus_render: crate::torus_render::TorusRender,
    burn_frames: Option<u32>,
    show_molecule: bool,
    markers: Vec<crate::molecule_view::LayerMarker>,
}

impl App {
    pub fn new() -> Self {
        let configs = presets::all_presets();
        let caches: Vec<render::CachedDomain> = configs
            .iter()
            .map(|p| render::CachedDomain::new(&p.domain))
            .collect();

        let mut app = Self {
            configs,
            caches,
            idx: 0,
            second_int: 0.2, // Λ for confocal, θ/π for polyline
            slider: Slider::new(),
            show_3d: false,
            animate: false,
            anim_dir: 1.0,
            cam3d: phase3d::OrbitCamera3::new(),
            view: ViewState::empty(),
            torus_render: crate::torus_render::TorusRender::new(),
            burn_frames: None,
            show_molecule: true,
            markers: Vec::new(),
        };
        // Build the initial view against the real configs.
        app.view = ViewState::new(app.idx, &app.configs, app.second_int, app.show_3d);
        app.rebuild_markers();
        app
    }

    /// Parse the `--burn[=N]` profiling hook: force 3D mode, build the torus,
    /// run N frames, then exit — so `perf`/flamegraph get a bounded run.
    fn parse_burn(&mut self) {
        for arg in std::env::args().skip(1) {
            if let Some(n) = arg.strip_prefix("--burn=") {
                self.burn_frames = n.parse().ok();
            } else if arg == "--burn" {
                self.burn_frames = Some(300);
            }
        }
        if self.burn_frames.is_some() {
            self.show_3d = true;
            self.rebuild();
        }
    }

    fn rebuild(&mut self) {
        self.view
            .rebuild(self.idx, &self.configs, self.second_int, self.show_3d);
        self.rebuild_markers();
    }

    /// Recompute the special-layer markers for the current preset.
    fn rebuild_markers(&mut self) {
        self.markers = crate::molecule_view::layer_markers(&self.configs[self.idx]);
    }

    /// Snap the second integral to the previous / next special layer.
    fn snap_special(&mut self, dir: i32) {
        if self.markers.is_empty() {
            return;
        }
        let current = self.second_int;
        if let Some(next) = crate::molecule_view::snap_nearest(&self.markers, current, dir) {
            self.second_int = next;
            self.rebuild();
        }
    }

    fn handle_input(&mut self) {
        // Toggle animation
        if is_key_pressed(KeyCode::A) {
            self.animate = !self.animate;
        }

        // Toggle the molecule strip overlay
        if is_key_pressed(KeyCode::M) {
            self.show_molecule = !self.show_molecule;
        }

        // Snap to previous / next special layer
        if is_key_pressed(KeyCode::LeftBracket) {
            self.snap_special(-1);
        }
        if is_key_pressed(KeyCode::RightBracket) {
            self.snap_special(1);
        }

        // Toggle 3D
        if is_key_pressed(KeyCode::P) {
            self.show_3d = !self.show_3d;
            self.rebuild();
        }

        // Tab / Space: next preset
        if is_key_pressed(KeyCode::Tab) || is_key_pressed(KeyCode::Space) {
            self.idx = (self.idx + 1) % self.configs.len();
            self.second_int = 0.2;
            self.rebuild();
        }

        let preset = &self.configs[self.idx];
        let mut changed = false;

        // Keyboard: adjust second integral (↑/→ step, ↓/← step back).
        let step = if is_key_pressed(KeyCode::Up) || is_key_pressed(KeyCode::Right) {
            KEY_STEP
        } else if is_key_pressed(KeyCode::Down) || is_key_pressed(KeyCode::Left) {
            -KEY_STEP
        } else {
            0.0
        };
        if step != 0.0 {
            if preset.is_confocal {
                // Step through the special layers (walls, separatrix, focal
                // axis, molecule critical values) so none can be skipped, and
                // fall back to the fixed step when the preset has no layers.
                let dir = step.signum() as i32;
                if let Some(next) =
                    crate::molecule_view::snap_nearest(&self.markers, self.second_int, dir)
                {
                    if (next - self.second_int).abs() > 1e-6 {
                        self.second_int = next;
                        changed = true;
                    }
                } else {
                    self.second_int = LambdaRange::of_domain(&preset.domain)
                        .clamp(self.second_int + step, self.second_int);
                    changed = true;
                }
            } else {
                self.second_int = (self.second_int + step).clamp(-1.0, 1.0);
                changed = true;
            }
        }

        // Slider input
        let slider_range = SecondIntegralRange::for_domain(&preset.domain, preset.is_confocal);
        if let Some(new_val) =
            self.slider
                .update(self.second_int, &slider_range, preset.second_integral_label)
        {
            if (new_val - self.second_int).abs() > ui::SLIDER_JITTER {
                self.second_int = new_val;
                changed = true;
            }
        }

        if changed {
            self.rebuild();
        }
    }

    fn update_animation(&mut self) {
        if !self.animate {
            return;
        }
        let speed = ANIM_SPEED; // units per second
        let dt = get_frame_time();
        let step = speed * dt;

        if self.configs[self.idx].is_confocal {
            let range = LambdaRange::of_domain(&self.configs[self.idx].domain);
            let (ell_min, ell_max, hyp_min, hyp_max) = range.animation_bounds();

            let mut next = self.second_int + self.anim_dir * step;
            // Bounce off the ends of the valid range.
            if next >= ell_max && next <= hyp_min {
                // Crossing the gap: jump to the other side.
                if self.anim_dir > 0.0 {
                    next = hyp_min;
                } else {
                    next = ell_max;
                }
            }
            if next >= hyp_max {
                next = hyp_max;
                self.anim_dir = -1.0;
            }
            if next <= ell_min {
                next = ell_min;
                self.anim_dir = 1.0;
            }
            self.second_int = next;
        } else {
            let mut next = self.second_int + self.anim_dir * step;
            if next >= 1.0 {
                next = 1.0;
                self.anim_dir = -1.0;
            }
            if next <= -1.0 {
                next = -1.0;
                self.anim_dir = 1.0;
            }
            self.second_int = next;
        }

        self.rebuild();
    }

    fn draw(&mut self, w: f32, h: f32) {
        let cache = &self.caches[self.idx];
        let cam = &cache.camera;
        let preset = &self.configs[self.idx];

        clear_background(render::BG);

        if self.show_3d {
            // ---- 3D phase space view ----
            self.cam3d.handle_input();
            phase3d::draw_axes(&self.cam3d, w, h, cache.domain_extent);

            // Pseudo-integrable table (the L): draw the flat surface (torus for
            // torus levels, unfolded cross for genus-2).  Otherwise fall back to
            // the cached torus render.
            if let Some(sheet) = &self.view.critical_sheet {
                // Critical layer: the phase manifold collapsed to a thin 2D sheet.
                phase3d::draw_critical_sheet(sheet, &self.cam3d, w, h);
            } else if let Some(level) = &self.view.flat_level {
                crate::pseudo::draw_flat(&self.view.flat_trajectories, level, &self.cam3d, w, h);
                // Example trajectories (red) and boundary preimage (cyan/orange)
                // on the flat surface, matching the smooth torus view.
                crate::pseudo::draw_flat_highlights(
                    &self.view.flat_highlights,
                    level,
                    &self.cam3d,
                    w,
                    h,
                );
                crate::pseudo::draw_flat_boundary(
                    &self.view.flat_boundary,
                    level,
                    &self.cam3d,
                    w,
                    h,
                );
            } else {
                // Dense points fill the 2D Liouville torus surface.  Rasterized
                // into an offscreen target and blitted; re-rasterized only when
                // the camera or geometry changes.  The torus size is driven by
                // how much of the domain the trajectory is allowed to reach
                // at the current Λ (see `accessible_region_scale`), so as Λ is
                // animated the torus grows/shrinks with the accessible region.
                // The billiard's boundary (walls + caustic) is drawn on the
                // surface in distinct colours.
                let scale = crate::torus_render::DomainScale::of_extent(
                    phase3d::accessible_region_scale(&preset.domain, self.second_int),
                );
                let ctx = crate::torus_render::DrawContext {
                    cam: &self.cam3d,
                    win_w: w,
                    win_h: h,
                    scale: &scale,
                };
                self.torus_render.draw_scaled(
                    &self.view.phase_trajectories,
                    &self.view.torus_highlights,
                    &self.view.boundary_preimage,
                    &ctx,
                );
            }

            let info =
                format!(
                "{}  |  {} = {:.3}  |  {} trajs  |  [P] 2D  |  [A] anim {}  |  right-drag orbit",
                preset.label,
                preset.second_integral_label,
                self.second_int,
                self.view.flat_trajectories.len().max(self.view.phase_trajectories.len()),
                if self.animate { "on" } else { "off" },
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
                render::draw_caustic(
                    ConfocalParams::standard(),
                    self.second_int,
                    &preset.domain,
                    cam,
                    w,
                    h,
                );
            } else {
                // Draw the starting direction arrow
                let angle = self.second_int * std::f32::consts::PI;
                render::draw_start_arrow(preset.start_center, angle, cam, w, h);
            }

            for traj in &self.view.trajectories {
                render::draw_trajectory(traj, cam, w, h);
            }
            // Highlight one short (few-bounce) red trajectory per torus so one
            // can see exactly where each torus's orbit lives on the billiard.
            let regime = crate::confocal::ConfocalStructure::of_domain(&preset.domain)
                .map(|s| s.regime(self.second_int))
                .unwrap_or(crate::confocal::TorusRegime::Single);
            let picks = render::one_per_torus(&self.view.start_points, regime);
            for &i in &picks {
                render::draw_trajectory_red(&self.view.trajectories[i], 4, cam, w, h);
            }
            for &(p, _) in &self.view.start_points {
                render::draw_start_marker(p, cam, w, h);
            }

            let total: usize = self.view.trajectories.iter().map(|t| t.len()).sum();
            let info = format!(
                "{}  |  {} = {:.3}  |  {} trajs, {} bounces  |  4 bounce/torus red  |  ↑↓  |  [A] anim {}  |  [P] 3D  |  Tab next",
                preset.label, preset.second_integral_label, self.second_int, self.view.start_points.len(), total,
                if self.animate { "on" } else { "off" },
            );
            draw_text(&info, 12.0, 28.0, 18.0, color_u8!(200, 200, 220, 220));

            // Molecule strip: the Reeb graph over λ, with a marker per special
            // layer and a cursor at the current Λ.  `[` / `]` snap Λ to the
            // previous / next special layer, so every one is reachable.
            if self.show_molecule && preset.is_confocal {
                let n = crate::molecule_view::draw_strip(
                    &self.markers,
                    self.second_int,
                    w * 0.1,
                    h * 0.82,
                    w * 0.8,
                );
                draw_text(
                    &format!("{n} special layers  |  [ ] snap  |  [M] hide"),
                    w * 0.1,
                    h * 0.82 + 30.0,
                    14.0,
                    color_u8!(150, 150, 190, 200),
                );
            }

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
    }

    /// Run the main loop until the window closes (or `--burn` exhausts).
    pub async fn run(mut self) {
        self.parse_burn();

        loop {
            let (w, h) = (screen_width(), screen_height());

            // Position slider
            self.slider.x = w * 0.15;
            self.slider.y = h - 50.0;
            self.slider.width = w * 0.7;

            self.handle_input();
            self.update_animation();
            self.draw(w, h);

            next_frame().await;

            if let Some(n) = self.burn_frames {
                if n <= 1 {
                    break;
                }
                self.burn_frames = Some(n - 1);
            }
        }
    }
}
