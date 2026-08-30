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
//! At the band ends the *main lobe itself* thins — but unlike the smooth
//! confocal map (whose θ₁ is *always* the collapsing libration by
//! construction), the flat chart has no such guarantee: `u1` is the
//! ell-side length and `u2` the hyp-side length, and which one actually
//! shrinks depends on which wall is being approached.  At birth (λ→0, an
//! ell-type wall) `u1_extent → 0` — the existing tube convention is right.
//! At death (λ→hyp_max, the outer hyp-type wall) it is `u2_extent → 0`
//! instead — the *major* radius must collapse there, not the tube, or the
//! embedding normalises a near-zero extent up to a full angular sweep and
//! renders numerically-jittery noise instead of a clean degenerate limit.
//! `tube_collapse`/`major_collapse` are therefore two independent factors,
//! each a smoothstep of the distance to its own end, mirroring the
//! smooth-table morph (`DomainScale::morph_to_level`, WINDOW = 0.35 of the
//! band).

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
    /// Main-lobe tube (minor) radius factor ∈ [0, 1]: 0 at birth (λ→0, where
    /// `u1_extent` vanishes), 1 elsewhere.
    pub tube_collapse: f32,
    /// Main-lobe major radius factor ∈ [0, 1]: 0 at death (λ→hyp_max, where
    /// `u2_extent` vanishes), 1 elsewhere.
    pub major_collapse: f32,
}

fn smoothstep(x: f32) -> f32 {
    let x = x.clamp(0.0, 1.0);
    x * x * (3.0 - 2.0 * x)
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

    let tube_collapse = smoothstep(lam.max(0.0) / (COLLAPSE_WINDOW * top));
    let major_collapse = smoothstep((top - lam).max(0.0) / (COLLAPSE_WINDOW * top));

    FlatMorph {
        handle_t,
        tube_collapse,
        major_collapse,
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
    let r = R_MINOR * m.tube_collapse;
    let r_major = R_MAJOR * m.major_collapse;

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
    let tube_r = HANDLE_TUBE * m.handle_t * m.tube_collapse.min(m.major_collapse).max(0.05);

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
    let r = R_MINOR * m.tube_collapse;
    let r_major = R_MAJOR * m.major_collapse;
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

    /// tube_collapse vanishes at birth (u1 vanishes there), major_collapse
    /// vanishes at death (u2 vanishes there); both saturate mid-band.
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
            m1.major_collapse < 0.05,
            "death: major_collapse = {}",
            m1.major_collapse
        );
        assert!(mm.tube_collapse > 0.95 && mm.major_collapse > 0.95, "mid");
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
