//! Flat-surface mapping for the genus-2 L: every phase point of an accessible
//! level must land on the flat translation surface (flat L / unfolded cross),
//! NOT on a donut `torus_embed`.
//!
//! The doc (`confocal_L_pseudo_integrable.md` §12) is explicit: a genus-2
//! surface has **no** flat torus picture — forcing it onto a donut is wrong and
//! renders as a garbled shape.  The correct object is the flat surface with the
//! `uᵢ = ∫ dλᵢ/√P` coordinates, where the flow has slope ±1 and the level set
//! unfolds to a translation surface.
//!
//! These tests pin the mapping (`to_flat`) over a dense trajectory on a
//! `GenusSurface` level: every point maps to a finite `(u1, u2)` with a sheet,
//! the sheets together form the flat L, and a segment between consecutive flat
//! points has slope ±1 (the doc's §13.2 "slope lock" sanity check, which is the
//! single strongest validation of the coordinates).

use billiards::pseudo::{classify_level, to_flat, Level};
use billiards::table::Table;
use billiards::torus::{ConfocalParams, PhaseSample};
use billiards::{domain, presets, pseudo::FlatError};

/// The standard-L preset and its confocal family.
fn standard_l() -> (domain::Domain, ConfocalParams) {
    let preset = presets::all_presets()
        .into_iter()
        .find(|p| p.label.contains("3π/2"))
        .expect("standard L preset");
    (preset.domain, ConfocalParams::standard())
}

fn table_of(dom: &domain::Domain) -> Table {
    Table::from_domain(dom, &ConfocalParams::standard()).expect("standard L table")
}

/// The accessible (genus-2) level for the L.
fn genus2_level(dom: &domain::Domain) -> Level {
    classify_level(0.7, &table_of(dom), &ConfocalParams::standard(), 1e-9, 1e-9)
}

/// A single straight bounce (segment) sampled densely, revealing the flat flow
/// direction.  `per_seg` points are taken along the bounce.
#[derive(Clone, Debug)]
struct FlatSegment {
    pts: Vec<FlatPoint>,
}

/// Map a trajectory's points through `to_flat`, grouped into straight segments
/// (each bounce).  Returns `Err` on the first point that fails.
fn flat_segments(
    dom: &domain::Domain,
    p0: macroquad::prelude::Vec2,
    v0: macroquad::prelude::Vec2,
    max_steps: usize,
    per_seg: usize,
    level: &Level,
) -> Result<Vec<FlatSegment>, FlatError> {
    let cf = ConfocalParams::standard();
    let mut segs = Vec::new();
    let mut p = p0;
    let mut v = v0;
    for _ in 0..max_steps {
        let speed = v.length();
        if speed < 1e-12 {
            break;
        }
        let dir = v / speed;
        let (_t, idx, hit) = match dom.intersect(p, dir) {
            Some(r) => r,
            None => break,
        };
        let mut pts = Vec::with_capacity(per_seg);
        for k in 0..per_seg {
            let f = k as f32 / per_seg as f32;
            let q = p + (hit - p) * f;
            let sample = PhaseSample::new(q.x, q.y, v.x, v.y);
            let flat = to_flat(&sample, &cf, level)?;
            pts.push(FlatPoint {
                u1: flat.u1,
                u2: flat.u2,
                sheet: flat.sheet,
                component: flat.component,
            });
        }
        segs.push(FlatSegment { pts });
        v = dom.reflect(hit, v, idx);
        p = hit + 1e-4 * v.normalize();
    }
    Ok(segs)
}

#[derive(Clone, Copy, Debug)]
struct FlatPoint {
    u1: f32,
    u2: f32,
    sheet: (i8, i8),
    component: u32,
}

/// A trajectory *on the λc=0.7 caustic*: launch from a start produced by
/// `dense_caustic_starts` at that level, so every phase point conserves
/// λc=0.7 and `to_flat` (built for λc=0.7) is the correct chart.
fn genus2_start() -> (
    domain::Domain,
    macroquad::prelude::Vec2,
    macroquad::prelude::Vec2,
) {
    let (dom, _) = standard_l();
    let starts = billiards::dense_caustic_starts(&dom, 0.7, 8);
    assert!(!starts.is_empty(), "no caustic starts at λc=0.7");
    let (p0, v0) = starts[0];
    (dom, p0, v0)
}

/// All flat points on a genus-2 accessible level must be finite and on the
/// same component (the standard L has a single connected accessible region).
#[test]
fn genus2_flat_mapping_is_finite_and_single_component() {
    let (dom, p0, v0) = genus2_start();
    let level = genus2_level(&dom);
    let Level::GenusSurface { region, .. } = &level else {
        panic!("λc=0.7 should be genus-2");
    };

    let segs = flat_segments(&dom, p0, v0, 60, 4, &level).expect("no flat errors");
    let n: usize = segs.iter().map(|s| s.pts.len()).sum();

    assert!(n > 50, "too few flat points: {}", n);
    for seg in &segs {
        for pt in &seg.pts {
            assert!(
                pt.u1.is_finite() && pt.u2.is_finite(),
                "flat coords must be finite"
            );
            assert!(pt.sheet.0 == 1 || pt.sheet.0 == -1);
            assert!(pt.sheet.1 == 1 || pt.sheet.1 == -1);
            assert_eq!(pt.component, 0, "L has one connected component");
        }
    }
    // The level really is genus 2 with one reflex corner.
    assert_eq!(region.genus, vec![2]);
    assert_eq!(region.n_components, 1);
}

/// Within a single straight bounce, consecutive flat points must have slope
/// ±1 (the doc's slope lock §13.2).  This validates the `u` coordinates, the
/// cubic, and the whole chart.
#[test]
fn flat_segments_have_slope_pm1() {
    let (dom, p0, v0) = genus2_start();
    let level = genus2_level(&dom);

    let segs = flat_segments(&dom, p0, v0, 4, 8, &level).expect("no flat errors");

    // Compare consecutive pairs *within each straight segment*: the flat image
    // of one straight bounce is a line with slope ±1 (the flow direction),
    // provided the segment does not cross a turning line mid-bounce.  Segments
    // that do cross a turning point are skipped (the doc's carry-and-flip note).
    let mut checked = 0;
    for seg in &segs {
        // Skip the segment if the sheet changes mid-segment (crosses a
        // turning point: a caustic-reflection at λc or a wall).
        let s0 = seg.pts[0].sheet;
        let crosses = seg.pts.iter().any(|p| p.sheet != s0);
        if crosses {
            continue;
        }
        for pair in seg.pts.windows(2) {
            let a = pair[0];
            let b = pair[1];
            let du1 = b.u1 - a.u1;
            let du2 = b.u2 - a.u2;
            if du1.abs() < 1e-6 && du2.abs() < 1e-6 {
                continue; // near-identical sample
            }
            let slope = if du1.abs() < 1e-9 {
                f32::INFINITY
            } else {
                du2 / du1
            };
            let is_pm1 = (slope - 1.0).abs() < 0.06 || (slope + 1.0).abs() < 0.06;
            assert!(
                is_pm1,
                "flat segment slope {:.3} not ±1 (du1={:.4}, du2={:.4})",
                slope, du1, du2
            );
            checked += 1;
        }
    }
    assert!(
        checked > 15,
        "not enough flat segments to validate slope ({checked})"
    );
}

/// A torus level (caustic shadows the corner) maps to `Level::Torus`, NOT a
/// genus surface — the ban of "angles → genus":
/// no reflex corner accessible → torus.
#[test]
fn no_corner_yields_torus_not_genus() {
    let (dom, _) = standard_l();
    // λc < α₁: strip, no reflex corner.
    let level = classify_level(
        0.3,
        &table_of(&dom),
        &ConfocalParams::standard(),
        1e-9,
        1e-9,
    );
    assert!(
        matches!(level, Level::Torus { .. }),
        "λc<α₁ must be a torus, not genus-2"
    );
    // β₂ < λc < β₃: strip, no reflex corner.
    let level = classify_level(
        2.3,
        &table_of(&dom),
        &ConfocalParams::standard(),
        1e-9,
        1e-9,
    );
    assert!(
        matches!(level, Level::Torus { .. }),
        "β₂<λc<β₃ must be a torus"
    );
}
