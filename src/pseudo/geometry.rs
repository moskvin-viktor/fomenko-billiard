//! Intrinsic geometry of a pseudo-integrable level, keyed to its reflex
//! corner.
//!
//! The renderer used to read the flat moduli positionally from the region's
//! grid lines (`u1[1]`, `u1.last()`, …), which breaks whenever the caustic cut
//! inserts an extra grid line.  This module extracts the quantities the
//! embedding actually needs — where the surface splits into "spine + handle",
//! how deep/thick the handle is, and where the pinch (reflex image) sits —
//! from the *accessible cells* and the *reflex corner*, so they are stable
//! under grid refinement.
//!
//! Geometry of the standard L (one reflex corner):
//!
//! - the **spine** is the full-height part `u1 ∈ [origin.0, split]`;
//! - the **handle** is the extra leg `u1 ∈ (split, split + d1]`, occupying the
//!   `u2` band `attach` of thickness `d2`.
//!
//! The two genus jumps are visible directly: at `λc → α₁⁺` the handle depth
//! `d1 → 0`; at `λc → β₂⁻` its thickness `d2 → 0`.  Torus levels have
//! `d1 = d2 = 0` and `split = w1`.

use super::classify::Level;

/// The intrinsic flat geometry of a level, in raw `u`-coordinates.
#[derive(Clone, Debug)]
pub struct LevelGeometry {
    /// Lower-left corner `(u1_min, u2_min)` of the accessible region (a
    /// hyperbolic cut shifts `u2_min` above 0).
    pub origin: (f32, f32),
    /// Raw `u1` of the reflex grid line separating spine and handle; equals
    /// `origin.0 + u1_extent` (the full width) for torus levels.
    pub split: f32,
    /// Handle depth in `u1` (0 for torus levels).
    pub d1: f32,
    /// Handle thickness in `u2` (0 for torus levels).
    pub d2: f32,
    /// Total accessible extent in `u1`.
    pub u1_extent: f32,
    /// Total accessible extent in `u2` (spine height).
    pub u2_extent: f32,
    /// Raw `u2` band `[lo, hi]` the handle attaches over (`hi - lo = d2`).
    pub attach: (f32, f32),
    /// Flat images of the surviving reflex corners (the pinch points),
    /// raw `(u1, u2)`.
    pub pinch: Vec<(f32, f32)>,
}

/// Extract the [`LevelGeometry`] of a classified level.
///
/// Torus levels are a plain rectangle.  Genus surfaces are keyed to the reflex
/// corner: the spine is everything at `u1 ≤ split` (full `u2` extent for the
/// standard L), the handle everything beyond.  Levels with no surviving reflex
/// corner (degenerate slivers) degrade gracefully to the rectangle case.
pub fn level_geometry(level: &Level) -> LevelGeometry {
    match level {
        Level::Torus { shape, .. } => {
            let (w1, w2) = *shape;
            LevelGeometry {
                origin: (0.0, 0.0),
                split: w1,
                d1: 0.0,
                d2: 0.0,
                u1_extent: w1,
                u2_extent: w2,
                attach: (0.0, 0.0),
                pinch: Vec::new(),
            }
        }
        Level::GenusSurface { region, .. } => {
            // Bounding box of the accessible cells (grid lines outside the
            // accessible set — e.g. the dead band below a hyperbolic cut —
            // must not count).
            let mut u1_min = f32::INFINITY;
            let mut u1_max = f32::NEG_INFINITY;
            let mut u2_min = f32::INFINITY;
            let mut u2_max = f32::NEG_INFINITY;
            for (i, row) in region.accessible.iter().enumerate() {
                for (j, &acc) in row.iter().enumerate() {
                    if !acc {
                        continue;
                    }
                    u1_min = u1_min.min(region.u1[i]);
                    u1_max = u1_max.max(region.u1[i + 1]);
                    u2_min = u2_min.min(region.u2[j]);
                    u2_max = u2_max.max(region.u2[j + 1]);
                }
            }
            if !u1_min.is_finite() {
                // No accessible cells at all — empty rectangle.
                return LevelGeometry {
                    origin: (0.0, 0.0),
                    split: 0.0,
                    d1: 0.0,
                    d2: 0.0,
                    u1_extent: 0.0,
                    u2_extent: 0.0,
                    attach: (0.0, 0.0),
                    pinch: Vec::new(),
                };
            }

            let pinch: Vec<(f32, f32)> = region
                .reflex_by_component
                .iter()
                .flat_map(|v| v.iter().copied())
                .collect();

            // The reflex grid line splits spine from handle.  Without one
            // (shouldn't happen for a genus surface, but degrade gracefully)
            // the region is treated as a plain rectangle.
            let Some(&(split, _)) = pinch.first() else {
                return LevelGeometry {
                    origin: (u1_min, u2_min),
                    split: u1_max,
                    d1: 0.0,
                    d2: 0.0,
                    u1_extent: u1_max - u1_min,
                    u2_extent: u2_max - u2_min,
                    attach: (u2_min, u2_min),
                    pinch,
                };
            };

            // The handle: accessible cells strictly beyond the split line.
            // Its u2 range is the attach band.
            let mut at_lo = f32::INFINITY;
            let mut at_hi = f32::NEG_INFINITY;
            for (i, row) in region.accessible.iter().enumerate() {
                if region.u1[i] < split - 1e-6 {
                    continue;
                }
                for (j, &acc) in row.iter().enumerate() {
                    if !acc {
                        continue;
                    }
                    at_lo = at_lo.min(region.u2[j]);
                    at_hi = at_hi.max(region.u2[j + 1]);
                }
            }
            let (attach, d1, d2) = if at_lo.is_finite() {
                ((at_lo, at_hi), u1_max - split, at_hi - at_lo)
            } else {
                ((u2_min, u2_min), 0.0, 0.0)
            };

            LevelGeometry {
                origin: (u1_min, u2_min),
                split,
                d1,
                d2,
                u1_extent: u1_max - u1_min,
                u2_extent: u2_max - u2_min,
                attach,
                pinch,
            }
        }
        Level::Forbidden | Level::Separatrix => LevelGeometry {
            origin: (0.0, 0.0),
            split: 0.0,
            d1: 0.0,
            d2: 0.0,
            u1_extent: 0.0,
            u2_extent: 0.0,
            attach: (0.0, 0.0),
            pinch: Vec::new(),
        },
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::super::classify::classify_level;
    use super::*;
    use crate::table::Table;
    use crate::torus::ConfocalParams;

    fn cf() -> ConfocalParams {
        ConfocalParams::new(4.0, 1.0)
    }

    /// The app's standard L: ell = [0, 0.4, 0.8], hyp = [1.4, 2.0, 2.6].
    fn standard_l_table() -> Table {
        Table::from_bounds(&[0.0, 0.4, 0.8], &[1.4, 2.0, 2.6], |l1, l2, _| {
            let tall = l1 <= 0.4 && l2 >= 2.0;
            let base = l1 <= 0.8 && l2 <= 2.0;
            tall || base
        })
    }

    fn geo_at(lc: f32) -> LevelGeometry {
        let level = classify_level(lc, &standard_l_table(), &cf(), 1e-9, 1e-9);
        level_geometry(&level)
    }

    /// Mid-band genus-2 levels: positive handle depth and thickness, one pinch,
    /// finite everything.
    #[test]
    fn test_genus2_geometry_is_positive() {
        for lc in [0.6f32, 0.9, 1.0, 1.2, 1.7] {
            let g = geo_at(lc);
            assert_eq!(g.pinch.len(), 1, "λc={lc}: one pinch point");
            assert!(g.d1 > 1e-4, "λc={lc}: d1={} must be positive", g.d1);
            assert!(g.d2 > 1e-4, "λc={lc}: d2={} must be positive", g.d2);
            assert!(g.u1_extent > g.d1, "λc={lc}: spine has positive width");
            assert!(g.u2_extent > g.d2 - 1e-6, "λc={lc}: spine spans the handle");
            assert!(
                (g.attach.1 - g.attach.0 - g.d2).abs() < 1e-6,
                "λc={lc}: attach band thickness must equal d2"
            );
            for v in [g.split, g.d1, g.d2, g.u1_extent, g.u2_extent] {
                assert!(v.is_finite(), "λc={lc}: non-finite geometry");
            }
        }
    }

    /// Genus jump at α₁ = 0.4: handle depth d1 → 0 from above.
    #[test]
    fn test_d1_vanishes_at_alpha1() {
        let near = geo_at(0.401);
        let far = geo_at(0.7);
        assert!(
            near.d1 < 0.2 * far.d1,
            "d1 near α₁ ({}) must be much smaller than mid-band ({})",
            near.d1,
            far.d1
        );
    }

    /// Genus jump at β₂ = 2.0: handle thickness d2 → 0 from below.
    #[test]
    fn test_d2_vanishes_at_beta2() {
        let near = geo_at(1.999);
        let far = geo_at(1.7);
        assert!(
            near.d2 < 0.2 * far.d2,
            "d2 near β₂ ({}) must be much smaller than mid-band ({})",
            near.d2,
            far.d2
        );
    }

    /// Torus levels are plain rectangles: no handle, no pinch.
    #[test]
    fn test_torus_geometry() {
        for lc in [0.2f32, 2.3] {
            let g = geo_at(lc);
            assert!(g.pinch.is_empty(), "λc={lc}: torus has no pinch");
            assert_eq!(g.d1, 0.0);
            assert_eq!(g.d2, 0.0);
            assert!(g.u1_extent > 0.0 && g.u2_extent > 0.0);
            assert!((g.split - g.u1_extent).abs() < 1e-6);
        }
    }
}
