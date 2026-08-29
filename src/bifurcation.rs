//! One source of truth for "what is the phase manifold at a caustic level λc".
//!
//! Previously the app, the molecule strip, and the tests each re-derived the
//! same classification with slightly different tolerances, and `ViewState::rebuild`
//! ran a fragile cascade (try critical sheet, else flat level, else torus).
//! This module centralises that decision:
//!
//! - [`classify`] — the single entry point: what is the manifold at this level?
//! - [`critical_values`] — every bifurcation point of a domain, in one place.
//! - [`degenerate_kind`] — which border piece a degenerate caustic collapsed onto.
//!
//! The low-level samplers (`caustic_starts`, `critical_caustic_starts`, the
//! pseudo `classify_level`) stay as point generators; this module orchestrates
//! them so callers never have to guess.

use crate::confocal::ConfocalStructure;
use crate::domain::Domain;
use crate::pseudo::Level;
use crate::table::Table;
use crate::torus::ConfocalParams;

/// The phase manifold at a caustic level `λc`.
#[derive(Clone, Debug)]
pub enum PhaseManifold {
    /// Generic integrable level: 1 or 2 Liouville tori.
    Tori {
        /// The torus regime (how many tori, and what splits them).
        regime: crate::confocal::TorusRegime,
    },
    /// Degenerate critical level: the manifold is a closed 1D orbit (a circle).
    Degenerate {
        /// Which border piece the caustic collapsed onto.
        kind: DegenerateKind,
    },
    /// Pseudo-integrable table (L/T/Z): a flat torus or a genus-≥2 surface.
    Flat {
        /// The classified flat level.
        level: Box<Level>,
    },
    /// Nothing reachable at this level.
    Forbidden,
}

/// Which border piece a degenerate caustic collapsed onto.
///
/// Note λ = λ_hyp (a hyperbola wall) is *not* degenerate: the caustic merely
/// coincides with the wall, the accessible region keeps its full area and the
/// level is an ordinary Liouville torus (trajectories graze the wall).  Only
/// the ellipse wall collapses the accessible annulus to zero width.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DegenerateKind {
    /// λ = b — the focal segment on the x-axis (degenerate ellipse).
    FocalSegment,
    /// λ = a — the vertical segment x = 0 (degenerate hyperbola).
    FocalAxis,
    /// λ = λ_ell — the outer ellipse wall.
    EllipseWall,
}

/// Classify the phase manifold at a caustic level.
///
/// Order of precedence:
/// 1. Pseudo-integrable table (L/T/Z) → `Flat` / `Forbidden`.
/// 2. Degenerate critical layer (wall, focal segment, focal axis) → `Degenerate`.
/// 3. Otherwise → `Tori`.
pub fn classify(domain: &Domain, cf: &ConfocalParams, lam: f32) -> PhaseManifold {
    // 1. Pseudo-integrable table.
    if let Some(tab) = Table::from_domain(domain, cf) {
        return match crate::pseudo::classify_level(lam, &tab, cf, 1e-9, 1e-9) {
            Level::Forbidden | Level::Separatrix => PhaseManifold::Forbidden,
            level => PhaseManifold::Flat {
                level: Box::new(level),
            },
        };
    }

    // 2. Degenerate critical layer.
    if let Some(kind) = degenerate_kind(domain, cf, lam) {
        return PhaseManifold::Degenerate { kind };
    }

    // 3. Generic integrable level.
    let structure = ConfocalStructure::of_domain(domain)
        .expect("a confocal domain must have a confocal structure");
    PhaseManifold::Tori {
        regime: structure.regime(lam),
    }
}

/// The degenerate border piece a caustic collapsed onto, `None` for regular levels.
///
/// A degenerate caustic is one where the phase manifold collapses to a 1D
/// orbit: the focal segment (`b`), the focal axis (`a`), or the ellipse wall
/// (`λ_ell`, where the accessible annulus between caustic and wall shrinks to
/// nothing).  The hyperbola wall `λ_hyp` is a *regular* level: the caustic
/// coincides with the wall but trajectories still fill a full 2D torus.
pub fn degenerate_kind(domain: &Domain, cf: &ConfocalParams, lam: f32) -> Option<DegenerateKind> {
    let structure = ConfocalStructure::of_domain(domain)?;
    let b = cf.b;
    let a = cf.a;

    if (lam - b).abs() < 1e-4 {
        return Some(DegenerateKind::FocalSegment);
    }
    if (lam - a).abs() < 1e-4 {
        return Some(DegenerateKind::FocalAxis);
    }
    if (lam - structure.lambda_ell).abs() < 1e-4 {
        return Some(DegenerateKind::EllipseWall);
    }
    None
}

/// Every bifurcation point of a domain, in one place: the molecule critical
/// values (from the rectilinear table when the domain is one), the ellipse
/// wall `λ_ell` (a degenerate layer, possibly at λ = 0), and the two
/// degeneracies `b` (focal separatrix) and `a` (focal axis).
///
/// The hyperbola wall `λ_hyp` is *not* listed: it is a regular torus level
/// (the caustic coincides with the wall but the level keeps full area), so it
/// gets no special marker or snapping.
///
/// This is the single source the app navigation, the molecule strip, and the
/// tests all use, so the tolerance/ordering can't drift between them.
pub fn critical_values(domain: &Domain, cf: &ConfocalParams) -> Vec<f32> {
    let mut vals = Vec::new();

    // Molecule critical values from the rectilinear λ-chart table (L, T, Z, …).
    if let Some(tab) = Table::from_domain(domain, cf) {
        vals.extend(crate::molecule::critical_values(&tab, cf, 1e-6));
    }

    // The ellipse wall (degenerate: the accessible annulus collapses there).
    if let Some(s) = ConfocalStructure::of_domain(domain) {
        vals.push(s.lambda_ell);
    }

    // The two degeneracies.
    vals.push(cf.b);
    vals.push(cf.a);

    vals.retain(|&v| v >= 0.0 && v <= cf.a + 1e-6);
    vals.sort_by(|x, y| x.total_cmp(y));
    vals.dedup();
    vals
}
