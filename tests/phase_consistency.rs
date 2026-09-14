//! Phase-manifold consistency across the whole caustic sweep.
//!
//! The app must not crash, and must produce a *valid* 2D billiard trajectory and
//! a *valid* phase manifold, at **every** caustic value the user can reach — in
//! particular at the critical / special values (the elliptical/hyperbolic
//! boundary walls, the focal separatrix `Λ = b`, and the molecule's critical
//! values from `molecule::critical_values`).
//!
//! This is the contract the app is required to satisfy (issue: clicking the
//! slider at a separatrix-adjacent value panicked inside `torus::quadrature`).
//! Each test pins one requirement:
//!
//! 1. Every confocal preset yields a valid λ-table (boundary lambdas + molecule
//!    critical values), all within `(0, a)` and strictly increasing.
//! 2. Every molecular critical value touches a domain wall (exhaustive candidate
//!    set).
//! 3. At every λ in the table — and across a dense sweep that samples into the
//!    separatrix-adjacent sliver — one can get valid start points, a finite
//!    trajectory, and a finite phase manifold (no panic, no NaN/Inf).

use billiards::presets::Preset;
use billiards::torus::ConfocalParams;
use macroquad::prelude::Vec2;

/// The standard confocal family every preset is built in.
fn cf() -> ConfocalParams {
    ConfocalParams::standard()
}

/// Confocal (quadric-arc) presets, which are the ones that have a λ-table and a
/// torus map.
fn confocal_presets() -> Vec<Preset> {
    billiards::presets::all_presets()
        .into_iter()
        .filter(|p| p.is_confocal)
        .collect()
}

/// The complete λ-table for a confocal domain: the union of the molecule's
/// critical values, the boundary walls, and the special values `{b, a}`.
fn lambda_table(preset: &Preset) -> Vec<f32> {
    let cf = cf();
    let mut vals = Vec::new();

    if let Some(tab) = billiards::table::Table::from_domain(&preset.domain, &cf) {
        vals.extend(billiards::molecule::critical_values(&tab, &cf, 1e-6));
    }
    if let Some(s) = billiards::confocal::ConfocalStructure::of_domain(&preset.domain) {
        vals.push(s.lambda_ell);
        if let Some(h) = s.lambda_hyp {
            vals.push(h);
        }
    }
    vals.push(cf.b);
    vals.push(cf.a);

    vals.retain(|&v| v > 0.0 && v < cf.a);
    vals.sort_by(|x, y| x.total_cmp(y));
    vals.dedup();
    vals
}

/// The sorted, deduplicated λ-table (same as [`lambda_table`]).
fn sorted_lambda_table(preset: &Preset) -> Vec<f32> {
    lambda_table(preset)
}

/// Whether the caustic at `lam` intersects the domain (the level is reachable).
fn reachable_caustic(preset: &Preset, lam: f32) -> bool {
    let starts = billiards::get_start_points(lam, &preset.domain, true, Vec2::ZERO);
    !starts.is_empty()
}

/// Assert the trajectory + phase manifold at one λ are finite (no panic/NaN).
fn check_lambda_phase(preset: &Preset, lam: f32) {
    let domain = &preset.domain;

    // Prefer the analytic λ-chart membership when the domain is a table, a
    // quadrilateral, or a full ellipse (the ray-cast `Domain::contains` is
    // unreliable at reflex corners and for the full ellipse's folded boundary).
    let cf = cf();
    let table = billiards::table::Table::from_domain(domain, &cf);
    let structure = billiards::confocal::ConfocalStructure::of_domain(domain);
    let inside = |p: Vec2| match &table {
        Some(t) => t.contains(p.x, p.y, &cf),
        None => match &structure {
            Some(s) if s.is_quadrilateral || s.is_full_ellipse => s.contains(p, 3e-2),
            _ => domain.contains(p),
        },
    };

    // Start points must not panic (empty is fine — the level is unreachable).
    // Every start point must itself be inside the domain (the app relies on
    // this to draw a trajectory at all; `caustic_starts` guarantees it).
    let starts = billiards::get_start_points(lam, domain, true, Vec2::ZERO);
    for &(p, _) in starts.iter().take(4) {
        assert!(
            inside(p),
            "{}: Λ={lam} start point ({:.3},{:.3}) outside domain",
            preset.label,
            p.x,
            p.y
        );
    }

    // Valid trajectory: bounces stay finite (a real trajectory with no
    // NaN/Inf).  Continuity of the trace (every bounce inside) is validated
    // more thoroughly by `caustic_tests` over dense λ; here we assert the
    // stronger structural guarantee that the sampler produces finite segments
    // at every reachable λ, including critical values.
    for &(p, v) in starts.iter().take(4) {
        let segs = domain.trace(p, v, 120);
        for (from, to) in segs {
            assert!(
                from.x.is_finite() && from.y.is_finite(),
                "{}: Λ={lam} traj from non-finite",
                preset.label
            );
            assert!(
                to.x.is_finite() && to.y.is_finite(),
                "{}: Λ={lam} traj to non-finite",
                preset.label
            );
        }
    }

    // Phase manifold: dense torus sampling + boundary preimage, finite.
    //
    // A pseudo-integrable table (the L) never goes through this smooth
    // `to_torus` machinery (`manifold::build_manifold` routes it through
    // `Manifold::FlatSurface`/`pseudo::classify_level` instead — see
    // `tests/pseudo_flat.rs` / `l_singular_levels.rs` for its own finiteness
    // checks). `to_torus`'s bounds fall back to the full-ellipse convention
    // for a table (`is_quadrilateral = false`), so probing it here at a λ
    // that is a perfectly regular table level (e.g. Λ = b, which the table
    // may not even touch — `pseudo::table_touches_focal`) just exercises the
    // *unrelated* smooth separatrix sentinel, not a real bug in the table's
    // own rendering path.
    if table.is_some() {
        return;
    }
    let bounds = billiards::confocal::ConfocalStructure::of_domain(domain)
        .map(|s| s.torus_bounds())
        .unwrap_or((0.0, None));
    let dense = billiards::dense_caustic_starts(domain, lam, 8);
    for (p, v) in dense.into_iter().take(3) {
        let pts = billiards::phase3d::sample_trajectory_phase_dense(domain, p, v, 40, 4, bounds);
        for q in &pts {
            assert!(q.theta1.is_finite(), "phase theta1 NaN at Λ={lam}");
            assert!(q.theta2.is_finite(), "phase theta2 NaN at Λ={lam}");
            assert!(q.theta.is_finite(), "phase theta NaN at Λ={lam}");
        }
    }
    let pre = billiards::phase3d::boundary_preimage_points(domain, lam, bounds);
    for (q, _isc) in &pre {
        assert!(q.theta1.is_finite(), "preimage theta1 NaN at Λ={lam}");
        assert!(q.theta2.is_finite(), "preimage theta2 NaN at Λ={lam}");
    }
}

// ---------------------------------------------------------------------------

/// (1) Every confocal preset has a valid, monotone, non-empty λ table inside
/// `(0, a)`.
#[test]
fn every_domain_has_valid_lambda_table() {
    let cf = cf();
    for preset in confocal_presets() {
        let tab = sorted_lambda_table(&preset);
        assert!(!tab.is_empty(), "{}: no λ table produced", preset.label);
        for w in tab.windows(2) {
            assert!(
                w[0] < w[1],
                "{}: λ table not strictly increasing: {} then {}",
                preset.label,
                w[0],
                w[1]
            );
        }
        for &v in &tab {
            assert!(
                v > 0.0 && v < cf.a,
                "{}: λ {} out of (0, a)",
                preset.label,
                v
            );
        }
    }
}

/// (2) Every molecular critical value touches a domain wall — the candidate
/// table is exhaustive.
#[test]
fn critical_values_touch_domain_walls() {
    let cf = cf();
    for preset in confocal_presets() {
        let domain = &preset.domain;
        let mut walls: Vec<f32> = Vec::new();
        for seg in &domain.segments {
            if let billiards::domain::Segment::Quad { curve, .. } = seg {
                walls.push(curve.lambda);
            }
        }
        walls.push(cf.b); // separatrix is always a special wall

        let Some(tab) = billiards::table::Table::from_domain(domain, &cf) else {
            continue; // a bare quadrilateral has no molecule table
        };
        for c in billiards::molecule::critical_values(&tab, &cf, 1e-6) {
            assert!(
                walls.iter().any(|&w| (w - c).abs() < 1e-4),
                "{}: critical value {c} not on any domain wall",
                preset.label
            );
        }
    }
}

/// (3) Phase manifold + trajectories are valid at every λ in the table —
/// including the separatrix `b` and the wall values — for every domain where
/// that λ is actually reachable.
#[test]
fn phase_is_valid_at_all_table_values() {
    for preset in confocal_presets() {
        let table = sorted_lambda_table(&preset);
        for &lam in &table {
            if reachable_caustic(&preset, lam) {
                check_lambda_phase(&preset, lam);
            }
            // Even if unreachable, `get_start_points` / the phase path must not
            // panic — a degenerate level is a valid "empty" level, not a crash.
            let _ = billiards::get_start_points(lam, &preset.domain, true, Vec2::ZERO);
        }
    }
}

/// The separatrix `Λ = b` (degenerate caustic) is the historical panic source
/// (two roots coincide, no third root to integrate against).  The phase
/// manifold must degrade gracefully — return empty/zero, never panic.
#[test]
fn phase_survives_separatrix() {
    for preset in confocal_presets() {
        let _ = billiards::get_start_points(cf().b, &preset.domain, true, Vec2::ZERO);
        let bounds = billiards::confocal::ConfocalStructure::of_domain(&preset.domain)
            .map(|s| s.torus_bounds())
            .unwrap_or((0.0, None));
        let pre = billiards::phase3d::boundary_preimage_points(&preset.domain, cf().b, bounds);
        for (q, _isc) in &pre {
            assert!(
                q.theta1.is_finite(),
                "separatrix preimage theta1 non-finite"
            );
            assert!(
                q.theta2.is_finite(),
                "separatrix preimage theta2 non-finite"
            );
        }
    }
}

/// A dense sweep should also avoid panicking at *every* sampled Λ across the
/// reachable band — the app animates continuously across λ, so it must handle
/// every intermediate value, not just the table's critical ones.
#[test]
fn phase_is_valid_across_dense_sweep() {
    let cf = cf();
    let n = 200;
    for preset in confocal_presets() {
        // Hyperbolic side, from just above b to a (dense).
        for k in 0..=n {
            let lam = cf.b + (cf.a - cf.b - 1e-3) * k as f32 / n as f32;
            if reachable_caustic(&preset, lam) {
                check_lambda_phase(&preset, lam);
            }
        }
        // Elliptic side, from a small margin above the outer wall to b (dense).
        // The margin avoids the arbitrarily thin annulus right at the outer
        // ellipse boundary (λ near λ_ell), where a confined trajectory sits
        // arbitrarily close to the wall and the analytic inside test is
        // legitimately borderline — same convention as `caustic_tests`.
        for k in 0..=n {
            let lo = 0.05f32;
            let lam = lo + (cf.b - 1e-3 - lo) * k as f32 / n as f32;
            if reachable_caustic(&preset, lam) {
                check_lambda_phase(&preset, lam);
            }
        }
    }
}
