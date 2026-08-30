//! Bifurcations on the pseudo-integrable L-table: every critical layer must
//! produce a manifold (never Empty inside the reachable band), the genus jumps
//! must carry singular info and pinch data, and the unified morph embedding
//! must be continuous across every molecule transition.

use billiards::manifold::{build_manifold, Manifold, SampleConfig};
use billiards::molecule::TransitionKind;
use billiards::pseudo::{
    classify_level, flat_morph, level_geometry, pinch_point_3d, sheets, unified_embed, Level,
};
use billiards::table::Table;
use billiards::torus::ConfocalParams;

fn cf() -> ConfocalParams {
    ConfocalParams::standard()
}

/// The app's confocal L preset.
fn l_preset() -> billiards::presets::Preset {
    billiards::presets::all_presets()
        .into_iter()
        .find(|p| p.label.contains("L-shape (confocal"))
        .expect("confocal L preset exists")
}

fn l_table() -> Table {
    Table::from_domain(&l_preset().domain, &cf()).expect("L preset is a table domain")
}

// ---------------------------------------------------------------------------
// 1. Every marker yields a manifold; genus jumps carry singular pinch data.
// ---------------------------------------------------------------------------

#[test]
fn every_marker_yields_a_manifold() {
    let preset = l_preset();
    let cf = cf();
    let config = SampleConfig::default();
    let markers = billiards::bifurcation::critical_values(&preset.domain, &cf);
    assert!(!markers.is_empty());

    for &lam in &markers {
        let m = build_manifold(&preset.domain, lam, &config);
        assert!(
            !matches!(m, Manifold::Empty),
            "marker λ={lam} must not be Empty"
        );
    }
}

#[test]
fn genus_jumps_are_singular_pinched_levels() {
    let preset = l_preset();
    let config = SampleConfig::default();
    for lam in [0.4f32, 2.0] {
        let m = build_manifold(&preset.domain, lam, &config);
        let Manifold::FlatSurface {
            level, singular, ..
        } = m
        else {
            panic!("λ={lam} must be a FlatSurface, got something else");
        };
        let s = singular.expect("genus jump must carry SingularInfo");
        assert_eq!(s.kind, TransitionKind::GenusJump, "λ={lam}");
        assert!((s.lc - lam).abs() < 1e-5, "λ={lam}: lc={}", s.lc);
        // The carried level is the genus-side one-sided limit: a genus surface
        // with exactly one pinch (reflex image).
        let Level::GenusSurface { region, .. } = &level else {
            panic!("λ={lam}: singular level must be the genus-side limit");
        };
        let pinches: usize = region.reflex_by_component.iter().map(|v| v.len()).sum();
        assert_eq!(pinches, 1, "λ={lam}: exactly one pinch");
        // And the morph has (nearly) collapsed the handle there.
        let g = level_geometry(&level);
        let mo = flat_morph(&g, &l_table(), &cf(), s.side_lam);
        assert!(
            mo.handle_t < 0.05,
            "λ={lam}: handle_t={} must be collapsed",
            mo.handle_t
        );
    }
}

// ---------------------------------------------------------------------------
// 2. Marker list: death cap in, unreachable a out.
// ---------------------------------------------------------------------------

#[test]
fn markers_include_death_cap_and_exclude_a() {
    let preset = l_preset();
    let cf = cf();
    let markers = billiards::bifurcation::critical_values(&preset.domain, &cf);
    let top = *l_table().hyp.last().unwrap();
    assert!(
        markers.iter().any(|&v| (v - top).abs() < 1e-5),
        "death cap {top} must be a marker: {markers:?}"
    );
    assert!(
        !markers.iter().any(|&v| (v - cf.a).abs() < 1e-4),
        "a={} is unreachable and must not be a marker: {markers:?}",
        cf.a
    );
    for &v in &markers {
        assert!(v >= 0.0 && v <= top + 1e-5, "marker {v} outside [0, {top}]");
    }
}

// ---------------------------------------------------------------------------
// 3. Band ends: thin tori, collapsed on whichever axis actually vanishes
//    there (u1/tube at birth, u2/major at death).
// ---------------------------------------------------------------------------

#[test]
fn band_ends_are_thin_collapsed_tori() {
    let preset = l_preset();
    let tab = l_table();
    let cf = cf();
    let top = *tab.hyp.last().unwrap();
    let config = SampleConfig::default();

    for (lam, tube_should_vanish) in [(1e-4f32, true), (top - 1e-4, false)] {
        let m = build_manifold(&preset.domain, lam, &config);
        let Manifold::FlatSurface { level, .. } = m else {
            panic!("λ={lam}: band end must be a FlatSurface");
        };
        assert!(
            matches!(level, Level::Torus { .. }),
            "λ={lam}: band end is a (thin) torus level"
        );
        let g = level_geometry(&level);
        let mo = flat_morph(&g, &tab, &cf, lam);
        if tube_should_vanish {
            assert!(
                mo.tube_collapse < 0.05,
                "λ={lam}: tube_collapse={} must vanish at birth",
                mo.tube_collapse
            );
        } else {
            assert!(
                mo.major_collapse < 0.05,
                "λ={lam}: major_collapse={} must vanish at death",
                mo.major_collapse
            );
        }
    }
}

// ---------------------------------------------------------------------------
// 4. Morph continuity across the genus jumps.
// ---------------------------------------------------------------------------

/// Embed a fixed fractional probe grid of the *spine* (u1 ≤ split), which
/// exists on both sides of a genus jump.
fn spine_probes(lam: f32) -> Vec<macroquad::prelude::Vec3> {
    let level = classify_level(lam, &l_table(), &cf(), 1e-9, 1e-9);
    let g = level_geometry(&level);
    let m = flat_morph(&g, &l_table(), &cf(), lam);
    let mut out = Vec::new();
    for f1 in [0.0f32, 0.3, 0.7, 1.0] {
        for f2 in [0.1f32, 0.5, 0.9] {
            let u1 = g.origin.0 + f1 * (g.split - g.origin.0);
            let u2 = g.origin.1 + f2 * g.u2_extent;
            for &(s1, s2) in &sheets() {
                out.push(unified_embed(u1, u2, s1, s2, &g, &m));
            }
        }
    }
    out
}

#[test]
fn morph_is_continuous_across_genus_jumps() {
    for crit in [0.4f32, 2.0] {
        // (a) Spine probes barely move across the jump.
        let below = spine_probes(crit - 1e-3);
        let above = spine_probes(crit + 1e-3);
        for (i, (a, b)) in below.iter().zip(&above).enumerate() {
            assert!(
                (*a - *b).length() < 0.05,
                "crit {crit}: spine probe {i} jumped {} ({a:?} → {b:?})",
                (*a - *b).length()
            );
        }

        // (b) On the genus side, the handle sits within a shrinking ball
        // around the pinch point as λ → crit.
        let genus_side = |d: f32| if crit < 1.0 { crit + d } else { crit - d };
        let mut prev = f32::INFINITY;
        for d in [0.05f32, 0.01, 1e-3, 1e-4] {
            let lam = genus_side(d);
            let level = classify_level(lam, &l_table(), &cf(), 1e-9, 1e-9);
            let g = level_geometry(&level);
            let m = flat_morph(&g, &l_table(), &cf(), lam);
            let pp = pinch_point_3d(&g, &m).expect("genus side has a pinch");
            let mut maxd = 0.0f32;
            for p in [0.25f32, 0.75, 1.0] {
                for q in [0.0f32, 0.5, 1.0] {
                    let u1 = g.split + p * g.d1;
                    let u2 = g.attach.0 + q * g.d2;
                    for &(s1, s2) in &sheets() {
                        let pos = unified_embed(u1, u2, s1, s2, &g, &m);
                        maxd = maxd.max((pos - pp).length());
                    }
                }
            }
            assert!(
                maxd <= prev + 1e-6,
                "crit {crit}: handle extent must shrink monotonically ({maxd} > {prev})"
            );
            prev = maxd;
        }
        assert!(
            prev < 0.05,
            "crit {crit}: handle must collapse onto the pinch (extent {prev})"
        );
    }
}

// ---------------------------------------------------------------------------
// 5. Edge swaps: the level geometry is continuous across 0.8, 1.0 and 1.4.
// ---------------------------------------------------------------------------

#[test]
fn geometry_is_continuous_across_edge_swaps() {
    for crit in [0.8f32, 1.0, 1.4] {
        let geo = |lam: f32| {
            let level = classify_level(lam, &l_table(), &cf(), 1e-9, 1e-9);
            level_geometry(&level)
        };
        let a = geo(crit - 1e-3);
        let b = geo(crit + 1e-3);
        for (name, va, vb) in [
            ("u1_extent", a.u1_extent, b.u1_extent),
            ("u2_extent", a.u2_extent, b.u2_extent),
            ("d1", a.d1, b.d1),
            ("d2", a.d2, b.d2),
            ("split-origin", a.split - a.origin.0, b.split - b.origin.0),
        ] {
            assert!(
                (va - vb).abs() < 0.1,
                "crit {crit}: {name} jumped {va} → {vb}"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// 6. Singular info on every molecule critical value.
// ---------------------------------------------------------------------------

#[test]
fn every_molecule_critical_has_singular_info() {
    let tab = l_table();
    let cf = cf();
    let mut crits = billiards::molecule::critical_values(&tab, &cf, 1e-6);
    crits.push(*tab.hyp.last().unwrap());
    for lc in crits {
        let s = billiards::bifurcation::table_singular(&tab, &cf, lc)
            .unwrap_or_else(|| panic!("critical value {lc} must carry singular info"));
        assert!((s.lc - lc).abs() < 1e-5);
        assert!(
            (s.side_lam - lc).abs() < 1e-3,
            "side_lam {} too far from lc {lc}",
            s.side_lam
        );
        let level = classify_level(s.side_lam, &tab, &cf, 1e-9, 1e-9);
        assert!(
            !matches!(level, Level::Forbidden | Level::Separatrix),
            "one-sided limit at {lc} must be a real level"
        );
    }
}
