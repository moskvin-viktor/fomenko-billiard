//! Torus ↔ genus-2 distinction for the standard L, driven by the *reflex
//! corners present in the caustic-truncated accessible region*.
//!
//! The pseudo-integrable principle (Richards–Berry/Dragović–Radnović): a level
//! set is a torus iff the accessible region contains **no** reflex (3π/2)
//! corner; it is a genus-`1 + n` surface iff `n` reflex corners survive the
//! caustic cut.  For the standard L (one reflex corner at `(α₁, β₂)`) this
//! means:
//!
//! - when the caustic shadows the corner → **torus**;
//! - when the corner is accessible → **genus 2** (exactly one reflex corner).
//!
//! These tests pin that boundary against the *real* L preset (not a hand-built
//! `Table`), asserting the two `λc` regimes and the bifurcation values from
//! `docs/math/pseudo_integrable_genus.md` §7.

use billiards::pseudo::classify_level;
use billiards::table::Table;
use billiards::torus::ConfocalParams;
use billiards::{domain, presets, pseudo::Level};

/// The standard-L preset and its confocal family.
fn standard_l() -> (domain::Domain, ConfocalParams) {
    let preset = presets::all_presets()
        .into_iter()
        .find(|p| p.label.contains("3π/2"))
        .expect("standard L preset");
    (preset.domain, ConfocalParams::standard())
}

fn table_of(dom: &domain::Domain) -> Table {
    // `from_domain` reconstructs the standard L's 2×2 cell grid.
    Table::from_domain(dom, &ConfocalParams::standard()).expect("standard L table")
}

fn classify(dom: &domain::Domain, lc: f32) -> Level {
    classify_level(lc, &table_of(dom), &ConfocalParams::standard(), 1e-9, 1e-9)
}

// ---------------------------------------------------------------------------
// No reflex corner accessible → torus
// ---------------------------------------------------------------------------

/// λc < α₁: the elliptic caustic leaves only the vertical strip
/// `[0, λc] × [β₁, β₃]`, which contains no reflex corner → torus.
#[test]
fn caustic_below_corner_is_torus() {
    let (dom, _) = standard_l();
    // α₁ = 0.4.  A λc below it keeps only λ₁ ≤ 0.3 → rectangle → torus.
    let level = classify(&dom, 0.3);
    assert!(
        matches!(level, Level::Torus { .. }),
        "λc=0.3 (< α₁) should be a torus, got {:?}",
        level_kind(&level)
    );
}

/// β₂ < λc < β₃: the hyperbolic caustic forbids the whole base, leaving only
/// the tall-leg strip `[0, α₁] × [λc, β₃]`, which contains no reflex corner →
/// torus.
#[test]
fn caustic_above_corner_is_torus() {
    let (dom, _) = standard_l();
    // β₂ = 2.0, β₃ = 2.6.
    let level = classify(&dom, 2.3);
    assert!(
        matches!(level, Level::Torus { .. }),
        "λc=2.3 (β₂<λc<β₃) should be a torus, got {:?}",
        level_kind(&level)
    );
}

// ---------------------------------------------------------------------------
// Reflex corner accessible → genus 2 (exactly one corner for the L)
// ---------------------------------------------------------------------------

/// α₁ < λc < b: the caustic leaves the L (both legs) with the reflex corner at
/// `(α₁, β₂)` accessible → genus 2, exactly one reflex corner, one component.
#[test]
fn elliptical_caustic_with_corner_is_genus2() {
    let (dom, _) = standard_l();
    // α₁=0.4 < λc=0.7 < b=1.0.
    let level = classify(&dom, 0.7);
    match level {
        Level::GenusSurface { region, .. } => {
            assert_eq!(region.genus, vec![2], "genus should be 2");
            assert_eq!(region.n_components, 1, "one connected component");
            assert_eq!(
                region.reflex_by_component[0].len(),
                1,
                "exactly one reflex corner for the L"
            );
        }
        other => panic!("λc=0.7 should be genus-2, got {:?}", other),
    }
}

/// b < λc < β₁: the whole L is accessible, still with its one reflex corner →
/// genus 2.
#[test]
fn hyperbolic_caustic_with_corner_is_genus2() {
    let (dom, _) = standard_l();
    // b=1.0 < λc=1.2 < β₁=1.4.
    let level = classify(&dom, 1.2);
    match level {
        Level::GenusSurface { region, .. } => {
            assert_eq!(region.genus, vec![2], "genus should be 2");
            assert_eq!(region.n_components, 1);
            assert_eq!(region.reflex_by_component[0].len(), 1);
        }
        other => panic!("λc=1.2 should be genus-2, got {:?}", other),
    }
}

// ---------------------------------------------------------------------------
// Degenerate / forbidden
// ---------------------------------------------------------------------------

/// λc = b on the standard L is a *regular* genus-2 level, not a separatrix:
/// the table never touches the focal segment (ell max 0.8 < b < 1.4 hyp min),
/// so the caustic degeneracy happens outside the table and every u-interval
/// stays a positive distance from b.  `Level::Separatrix` is reserved for
/// focal-touching tables (see `table_touches_focal`).
#[test]
fn b_on_standard_l_is_regular_genus2() {
    let (dom, cf) = standard_l();
    assert!(!billiards::pseudo::table_touches_focal(&table_of(&dom), &cf));
    let level = classify(&dom, 1.0);
    match level {
        Level::GenusSurface { region, .. } => {
            assert_eq!(region.genus, vec![2], "λc=b: genus should be 2");
            assert_eq!(region.n_components, 1);
        }
        other => panic!("λc=b should be regular genus-2, got {:?}", level_kind(&other)),
    }
}

/// λc ≥ β₃ → nothing accessible.
#[test]
fn beyond_outer_wall_is_forbidden() {
    let (dom, _) = standard_l();
    let level = classify(&dom, 2.8);
    assert!(
        matches!(level, Level::Forbidden),
        "λc ≥ β₃ should be forbidden, got {:?}",
        level_kind(&level)
    );
}

// ---------------------------------------------------------------------------
// Bifurcation values: the topology switches exactly when the caustic passes
// the reflex corner (λc = α₁ and λc = β₂).
// ---------------------------------------------------------------------------

/// Just below α₁ → torus; just above α₁ → genus 2.
#[test]
fn bifurcation_at_alpha1() {
    let (dom, _) = standard_l();
    // α₁ = 0.4.  Below → torus, above → genus-2.
    let below = classify(&dom, 0.38);
    let above = classify(&dom, 0.42);
    assert!(
        matches!(below, Level::Torus { .. }),
        "below α₁ must be torus"
    );
    assert!(
        matches!(&above, Level::GenusSurface { .. }),
        "above α₁ must be genus-2"
    );
}

/// Just below β₂ → genus 2; just above β₂ → torus.
#[test]
fn bifurcation_at_beta2() {
    let (dom, _) = standard_l();
    // β₂ = 2.0.  Below → genus-2 (corner accessible), above → torus (base gone).
    let below = classify(&dom, 1.9);
    let above = classify(&dom, 2.1);
    assert!(
        matches!(&below, Level::GenusSurface { .. }),
        "below β₂ must be genus-2"
    );
    assert!(
        matches!(&above, Level::Torus { .. }),
        "above β₂ must be torus"
    );
}

fn level_kind(level: &Level) -> &'static str {
    match level {
        Level::Forbidden => "forbidden",
        Level::Separatrix => "separatrix",
        Level::Torus { .. } => "torus",
        Level::GenusSurface { .. } => "genus≥2",
    }
}
