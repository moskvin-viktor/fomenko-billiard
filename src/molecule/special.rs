//! The special caustic layers of a domain, in ascending λ order.
//!
//! These are the layers the app must make reachable when navigating a domain:
//! the molecule's critical values (from the rectilinear λ-chart `Table` when
//! the domain is one), the boundary walls (ellipse / hyperbola), and the two
//! degeneracies `b` (focal separatrix) and `a` (vertical focal axis).  The same
//! list drives the molecule-side navigation and the phase-consistency tests, so
//! the app and the tests agree on what "all special layers" means.

use crate::table::Table;
use crate::torus::ConfocalParams;

/// The ordered, deduplicated special caustic values of a domain, all within
/// `(0, a)`.
pub fn special_layers(domain: &crate::domain::Domain, cf: &ConfocalParams) -> Vec<f32> {
    // Single source of truth: the bifurcation module owns every special value.
    crate::bifurcation::critical_values(domain, cf)
}

/// The complete ordered navigation list: every special layer *and* uniform
/// regular layers across the reachable sweep, deduplicated and strictly
/// increasing.
///
/// Regular layers are scattered uniformly across `(0, a)` and kept only where
/// the level is actually reachable (has a trajectory start point), so the UI
/// visits both the critical layers and enough regular layers in between, and
/// never steps onto an unreachable layer in a forbidden gap (e.g. the
/// separatrix gap `(b, λ_hyp)` of a quadrilateral).
pub fn all_layers(domain: &crate::domain::Domain, cf: &ConfocalParams) -> Vec<f32> {
    let specials = special_layers(domain, cf);
    let mut vals = specials.clone();

    // Uniform regular candidates across the whole caustic axis; keep those that
    // are reachable so forbidden gaps are skipped regardless of the preset's
    // reachable-band shape.
    const N: usize = 30;
    for k in 1..N {
        let lam = cf.a * k as f32 / N as f32;
        if crate::confocal::level_is_reachable(domain, cf, lam) {
            vals.push(lam);
        }
    }

    vals.retain(|&v| v > 0.0 && v <= cf.a + 1e-6);
    vals.sort_by(|x, y| x.total_cmp(y));
    vals.dedup();
    vals
}

/// The label of a special layer: which degeneracy / transition / wall it is.
/// Used by the 2D molecule strip to colour each marker.
///
/// For a value that is *not* one of the special layers, returns `"regular"` —
/// it is an ordinary reachable level in between the critical layers, not a
/// singularity.
pub fn layer_label(lam: f32, domain: &crate::domain::Domain, cf: &ConfocalParams) -> &'static str {
    if (lam - cf.b).abs() < 1e-4 {
        return "b (separatrix)";
    }
    if (lam - cf.a).abs() < 1e-4 {
        return "a (focal axis)";
    }
    if let Some(s) = crate::confocal::ConfocalStructure::of_domain(domain) {
        if (lam - s.lambda_ell).abs() < 1e-4 {
            return "ellipse wall";
        }
        if let Some(h) = s.lambda_hyp {
            if (lam - h).abs() < 1e-4 {
                return "hyperbola wall";
            }
        }
    }
    if let Some(tab) = Table::from_domain(domain, cf) {
        // Only classify as a transition if this is actually a molecule critical
        // value; otherwise it is a regular level.
        let crit = crate::molecule::critical_values(&tab, cf, 1e-6);
        if crit.iter().any(|&c| (c - lam).abs() < 1e-4) {
            let tr = crate::molecule::classify_transition(&tab, cf, lam, 1e-4);
            return tr.kind.label();
        }
    }
    "regular"
}
