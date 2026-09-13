//! Unified morphing embedding for pseudo-integrable levels.
//!
//! The old renderer picked a *different* surface per level kind — donut for
//! torus levels, pretzel for genus-2 — with incompatible axis conventions
//! (donut: u1→tube, u2→major; pretzel: u1→longitudinal, u2→tube), so the
//! embedding jumped discontinuously at every genus bifurcation.  This module
//! is one continuous family instead:
//!
//! - the **main lobe** is always a donut with `u1→tube`, `u2→major`
//!   (the repo-wide convention: θ₁/tube is the collapsing direction);
//! - the **handle** (the accessible cells beyond the reflex split line) is a
//!   small torus tangent to the main lobe at the pinch point, whose radii
//!   scale with [`FlatMorph::handle_t`].
//!
//! `handle_t` is driven by the *intrinsic* sliver sizes `d1`, `d2` of the
//! [`LevelGeometry`], so it is exactly 0 at both genus jumps (α₁: `d1→0`,
//! β₂: `d2→0`) and 1 mid-band: at a genus jump the handle has collapsed onto
//! the pinch point and the surface *is* the pinched torus, and the transition
//! is smooth from both sides.
//!
//! At the band ends the main lobe's *tube* thins (`tube_collapse`, a
//! smoothstep of the distance to either end, mirroring the smooth-table morph
//! `DomainScale::morph_to_level`'s `WINDOW = 0.35`).  The major radius is
//! **never** touched — this mirrors `DomainScale::morph_to_level` exactly,
//! which only ever shrinks `r_minor`/`gap` and leaves `r_major` fixed for
//! every degenerate limit (wall, separatrix, focal axis) regardless of which
//! one it is.
//!
//! An earlier version shrank the major radius instead at "death" (λ→hyp_max,
//! where `u2_extent → 0`, since `u2` plays the major-angle role here) on the
//! theory that dividing by a near-zero `u2_extent` would otherwise normalise
//! up to numerically-jittery noise. That noise never actually materialized —
//! `u1`/`u2` are computed from the same f32 source and stay well-conditioned
//! down to extents ~1e-2 — while shrinking the major radius did cause a real,
//! visible bug: once `r_major` dropped below `r` (tube), the ring torus
//! flipped into a self-intersecting spindle and the donut hole visibly
//! closed, well before the level was anywhere near singular (the window is a
//! coarse proxy on raw λ, not on the level's own shrinking extent). Fixing
//! the tube only, with the major radius structurally fixed at `R_MAJOR`,
//! makes that inversion impossible by construction.

use super::geometry::LevelGeometry;
use crate::table::Table;
use crate::torus::ConfocalParams;
use macroquad::prelude::*;

/// Main-lobe donut radii (match the old `donut_with_normal` so pure torus
/// levels look unchanged).
pub const R_MAJOR: f32 = 1.6;
pub const R_MINOR: f32 = 0.6;
/// Full-size handle loop / tube radii (reached mid-band, `handle_t = 1`).
const HANDLE_LOOP: f32 = 0.45;
const HANDLE_TUBE: f32 = 0.18;
/// Fraction of the λ band over which the tube collapses at birth/death.
const COLLAPSE_WINDOW: f32 = 0.35;
/// Fraction of the flat extent at which a handle sliver counts as full-size.
/// The sliver depths shrink like √(λ−λc) near a genus jump, so a generous
/// reference keeps the shrink visible over a reasonable λ range.
const SLIVER_REF: f32 = 0.5;

/// The morph state of the flat embedding at a caustic level.
#[derive(Clone, Copy, Debug)]
pub struct FlatMorph {
    /// Handle size ∈ [0, 1]: 0 at a genus jump (handle collapsed onto the
    /// pinch point), 1 mid-band (full genus-2 handle).
    pub handle_t: f32,
    /// Main-lobe tube (minor) radius factor ∈ [0, 1]: 0 at either band end
    /// (λ→0 or λ→hyp_max), 1 mid-band. The major radius is never scaled (see
    /// the module docs) so the hole can never close.
    pub tube_collapse: f32,
}

fn smoothstep(x: f32) -> f32 {
    let x = x.clamp(0.0, 1.0);
    x * x * (3.0 - 2.0 * x)
}

/// The main lobe's actual (tube, major) radii. `r_major` is always the fixed
/// `R_MAJOR` — see the module docs for why it must never be scaled — so
/// `r_major > R_MINOR >= r` unconditionally and the ring can never invert
/// into a self-intersecting spindle.
fn main_lobe_radii(m: &FlatMorph) -> (f32, f32) {
    (R_MINOR * m.tube_collapse, R_MAJOR)
}

/// The morph state at level `lam` of a table.
///
/// `handle_t` is intrinsic (from the level's own `d1`/`d2` slivers), so it
/// vanishes exactly when the handle does; `collapse` is extrinsic (distance
/// of λ to the reachable band ends `{0, hyp_max}`).
pub fn flat_morph(geo: &LevelGeometry, tab: &Table, _cf: &ConfocalParams, lam: f32) -> FlatMorph {
    let top = *tab.hyp.last().unwrap();

    let handle_t = if geo.d1 <= 0.0 || geo.d2 <= 0.0 {
        0.0
    } else {
        let r1 = geo.d1 / (SLIVER_REF * geo.u1_extent.max(1e-6));
        let r2 = geo.d2 / (SLIVER_REF * geo.u2_extent.max(1e-6));
        smoothstep(r1.min(r2))
    };

    // Thin the tube near *either* band end (birth at 0, death at hyp_max);
    // the major radius is fixed (see the module docs), so this only ever
    // makes the ring thinner, never inverts it.
    let birth = smoothstep(lam.max(0.0) / (COLLAPSE_WINDOW * top));
    let death = smoothstep((top - lam).max(0.0) / (COLLAPSE_WINDOW * top));
    let tube_collapse = birth * death;

    FlatMorph {
        handle_t,
        tube_collapse,
    }
}

/// Embed a flat point `(u1, u2)` on sheet `(σ1, σ2)` onto the morphing
/// surface, returning position and outward normal (for Lambert shading).
///
/// Main lobe (`u1 ≤ split`): a donut.  The four sheets tile it — `σ1` mirrors
/// the tube half-angle, `σ2` the major half-angle — so turning points
/// (`u = 0` / `u = extent`) are continuous seams.  The tube angle is offset by
/// π so the spine↔handle seam (`u1 = split`) sits on the *outer* equator.
///
/// Handle (`u1 > split`): a small torus whose center loop sticks radially out
/// of the main lobe at the pinch point, radii scaled by `m.handle_t`.  At
/// `handle_t = 0` every handle point coincides with the pinch point — the
/// pinched-torus singular level.
pub fn unified_with_normal(
    u1: f32,
    u2: f32,
    s1: i8,
    s2: i8,
    g: &LevelGeometry,
    m: &FlatMorph,
) -> (Vec3, Vec3) {
    use std::f32::consts::PI;

    let sg1 = s1 as f32;
    let sg2 = s2 as f32;
    let lu1 = (u1 - g.origin.0).clamp(0.0, g.u1_extent.max(0.0));
    let lu2 = (u2 - g.origin.1).clamp(0.0, g.u2_extent.max(0.0));
    let lsplit = (g.split - g.origin.0).max(1e-6);
    let h2 = g.u2_extent.max(1e-6);
    let (r, r_major) = main_lobe_radii(m);

    if lu1 <= lsplit + 1e-6 || g.d1 <= 0.0 || g.d2 <= 0.0 {
        // Main lobe: donut, u1→tube, u2→major.
        let tt = PI + sg1 * PI * (lu1.min(lsplit) / lsplit);
        let tm = sg2 * PI * (lu2 / h2);
        let (st, ct) = tt.sin_cos();
        let (sm, cm) = tm.sin_cos();
        let pos = vec3(
            (r_major + r * ct) * cm,
            (r_major + r * ct) * sm,
            r * st,
        );
        let n = vec3(ct * cm, ct * sm, st);
        return (pos, n);
    }

    // Handle: small torus tangent at the pinch point on the outer equator.
    let p = ((lu1 - lsplit) / g.d1).clamp(0.0, 1.0);
    let a_lo = g.attach.0 - g.origin.1;
    let q = ((lu2 - a_lo) / g.d2).clamp(0.0, 1.0);
    let psi = sg1 * PI * p; // loop angle (σ1 mirror closes the circle)
    let phi = sg2 * PI * q; // tube angle (σ2 mirror closes the circle)

    // Pinch point P: outer-equator main-lobe point at the attach-band middle.
    let tm = PI * (a_lo + 0.5 * g.d2) / h2;
    let e_r = vec3(tm.cos(), tm.sin(), 0.0);
    let e_t = vec3(-tm.sin(), tm.cos(), 0.0);
    let e_z = vec3(0.0, 0.0, 1.0);
    let pinch = (r_major + r) * e_r;

    let loop_r = HANDLE_LOOP * m.handle_t;
    let tube_r = HANDLE_TUBE * m.handle_t * m.tube_collapse.max(0.05);

    // Center loop through P: circle of radius loop_r around P + loop_r·e_r,
    // in the (e_r, e_z) plane; ψ = 0 at P.
    let center = pinch + loop_r * e_r;
    let radial = -psi.cos() * e_r + psi.sin() * e_z; // loop radial at ψ
    let n = phi.cos() * radial + phi.sin() * e_t;
    let pos = center + loop_r * radial + tube_r * n;
    (pos, n)
}

/// Position-only [`unified_with_normal`].
pub fn unified_embed(u1: f32, u2: f32, s1: i8, s2: i8, g: &LevelGeometry, m: &FlatMorph) -> Vec3 {
    unified_with_normal(u1, u2, s1, s2, g, m).0
}

/// The embedded pinch point (where the handle attaches / collapses), `None`
/// for levels without a pinch.  Used to draw a bright singular marker when
/// `handle_t` is small.
pub fn pinch_point_3d(g: &LevelGeometry, m: &FlatMorph) -> Option<Vec3> {
    if g.pinch.is_empty() {
        return None;
    }
    let a_lo = g.attach.0 - g.origin.1;
    let tm = std::f32::consts::PI * (a_lo + 0.5 * g.d2) / g.u2_extent.max(1e-6);
    let (r, r_major) = main_lobe_radii(m);
    Some(vec3(
        (r_major + r) * tm.cos(),
        (r_major + r) * tm.sin(),
        0.0,
    ))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::super::classify::classify_level;
    use super::super::geometry::level_geometry;
    use super::*;
    use crate::table::Table;
    use crate::torus::ConfocalParams;

    fn cf() -> ConfocalParams {
        ConfocalParams::new(4.0, 1.0)
    }

    fn standard_l_table() -> Table {
        Table::from_bounds(&[0.0, 0.4, 0.8], &[1.4, 2.0, 2.6], |l1, l2, _| {
            let tall = l1 <= 0.4 && l2 >= 2.0;
            let base = l1 <= 0.8 && l2 <= 2.0;
            tall || base
        })
    }

    fn morph_at(lam: f32) -> (LevelGeometry, FlatMorph) {
        let tab = standard_l_table();
        let level = classify_level(lam, &tab, &cf(), 1e-9, 1e-9);
        let g = level_geometry(&level);
        let m = flat_morph(&g, &tab, &cf(), lam);
        (g, m)
    }

    /// handle_t vanishes at both genus jumps and saturates mid-band.
    #[test]
    fn test_handle_t_profile() {
        let (_, m_a) = morph_at(0.4001);
        let (_, m_b) = morph_at(1.9999);
        let (_, m_mid) = morph_at(1.2);
        assert!(m_a.handle_t < 0.05, "α₁⁺: handle_t = {}", m_a.handle_t);
        assert!(m_b.handle_t < 0.05, "β₂⁻: handle_t = {}", m_b.handle_t);
        assert!(m_mid.handle_t > 0.9, "mid: handle_t = {}", m_mid.handle_t);
    }

    /// tube_collapse vanishes at *both* band ends (birth and death alike) and
    /// saturates mid-band. The major radius is never scaled — checked
    /// separately by `test_main_lobe_hole_never_closes`.
    #[test]
    fn test_collapse_profile() {
        let (_, m0) = morph_at(0.01);
        let (_, m1) = morph_at(2.59);
        let (_, mm) = morph_at(1.2);
        assert!(
            m0.tube_collapse < 0.05,
            "birth: tube_collapse = {}",
            m0.tube_collapse
        );
        assert!(
            m1.tube_collapse < 0.05,
            "death: tube_collapse = {}",
            m1.tube_collapse
        );
        assert!(mm.tube_collapse > 0.95, "mid: tube_collapse = {}", mm.tube_collapse);
    }

    /// At handle_t → 0 every handle point sits at the pinch point.
    #[test]
    fn test_handle_collapses_onto_pinch() {
        let (g, m) = morph_at(0.4001);
        let pp = pinch_point_3d(&g, &m).expect("genus level has a pinch");
        for p in [0.1f32, 0.5, 1.0] {
            for q in [0.0f32, 0.5, 1.0] {
                let u1 = g.split + p * g.d1;
                let u2 = g.attach.0 + q * g.d2;
                for &(s1, s2) in &super::super::render::sheets() {
                    let pos = unified_embed(u1, u2, s1, s2, &g, &m);
                    assert!(
                        (pos - pp).length() < 0.05,
                        "handle point {pos:?} far from pinch {pp:?}"
                    );
                }
            }
        }
    }

    /// Regression: the main lobe's major radius must never drop below its
    /// tube radius anywhere in the table's λ range — that inverts the ring
    /// torus into a self-intersecting spindle and the donut hole visibly
    /// closes (the reported bug). `main_lobe_radii` now makes this
    /// structurally impossible (`r_major` is the fixed `R_MAJOR`, never
    /// scaled), so this just guards the invariant going forward.
    #[test]
    fn test_main_lobe_hole_never_closes() {
        let tab = standard_l_table();
        let top = *tab.hyp.last().unwrap();
        let mut lam = 1e-3f32;
        while lam < top {
            let level = classify_level(lam, &tab, &cf(), 1e-9, 1e-9);
            let geo = level_geometry(&level);
            let m = flat_morph(&geo, &tab, &cf(), lam);
            let (r, r_major) = main_lobe_radii(&m);
            assert!(
                r_major >= r,
                "λ={lam}: r_major={r_major} < r={r} — donut hole closed"
            );
            lam += 0.01;
        }
    }

    /// The four sheets agree at the turning-point seams of the main lobe.
    #[test]
    fn test_sheet_seams_are_continuous() {
        let (g, m) = morph_at(1.2);
        // u1 = origin (tube seam) and u2 = origin / origin + extent (major seam).
        for &(u1, u2) in &[
            (g.origin.0, g.origin.1 + 0.3 * g.u2_extent),
            (g.origin.0 + 0.3 * g.split, g.origin.1),
            (g.origin.0 + 0.3 * g.split, g.origin.1 + g.u2_extent),
        ] {
            let base = unified_embed(u1, u2, 1, 1, &g, &m);
            for &(s1, s2) in &super::super::render::sheets() {
                // Seam points only: mirrored coordinates coincide there.
                let p = unified_embed(u1, u2, s1, s2, &g, &m);
                let same_u1_seam = (u1 - g.origin.0).abs() < 1e-6;
                let same_u2_seam = (u2 - g.origin.1).abs() < 1e-6
                    || (u2 - g.origin.1 - g.u2_extent).abs() < 1e-6;
                if (s1 == 1 || same_u1_seam) && (s2 == 1 || same_u2_seam) {
                    assert!(
                        (p - base).length() < 1e-4,
                        "seam mismatch at ({u1},{u2}) sheet ({s1},{s2})"
                    );
                }
            }
        }
    }
}
