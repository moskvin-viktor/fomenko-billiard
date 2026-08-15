//! Per-sample mapping into flat coordinates for pseudo-integrable tables.
//!
//! This is the generalized `to_flat` of the doc (§11 of
//! `confocal_L_pseudo_integrable.md`): map a phase-space sample `(x, y, vx, vy)`
//! to `(u₁, u₂)` plus the sheet `(sign d₁, sign d₂)` and the fold quadrant
//! `(sign x, sign y)`.

use super::classify::Level;
use crate::table::Quadrant;
use crate::torus::{lam_dots, ConfocalParams, PhaseSample};

/// A sample mapped into flat coordinates.
#[derive(Clone, Copy, Debug)]
pub struct FlatSample {
    /// The un-normalized length coordinate `u₁(λ₁)`.
    pub u1: f32,
    /// The un-normalized length coordinate `u₂(λ₂)`.
    pub u2: f32,
    /// The momentum sheet `(sign d₁, sign d₂)` — which of the four diagonal
    /// directions the flow is moving in.
    pub sheet: (i8, i8),
    /// The fold quadrant `(sign x, sign y)` — which physical copy of the
    /// λ-chart the sample lives on.
    pub fold: Quadrant,
    /// Which connected component of the accessible region this sample is in.
    pub component: u32,
}

/// Errors that can occur during the flat mapping.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FlatError {
    /// The level is the separatrix (λc = b).
    Separatrix,
    /// The level is forbidden (no accessible points).
    Forbidden,
}

/// Map one phase-space sample to its flat coordinates.
///
/// Precondition: `level` was obtained from `classify_level(lc, …)` for the
/// sample's caustic value.  Callers should precompute the level once per
/// trajectory and pass it to every sample.
pub fn to_flat(
    s: &PhaseSample,
    cf: &ConfocalParams,
    level: &Level,
) -> Result<FlatSample, FlatError> {
    let (lam1, lam2) = crate::torus::confocal(s.x, s.y, cf);
    let (d1, d2) = lam_dots(s, lam1, lam2, cf);

    let (u1, u2) = match level {
        Level::Forbidden => return Err(FlatError::Forbidden),
        Level::Separatrix => return Err(FlatError::Separatrix),
        Level::Torus { u1, u2, .. } => (u1.u(lam1), u2.u(lam2)),
        Level::GenusSurface { u1, u2, .. } => (u1.u(lam1), u2.u(lam2)),
    };

    let sheet = (
        if d1 >= 0.0 { 1 } else { -1 },
        if d2 >= 0.0 { 1 } else { -1 },
    );
    let fold = Quadrant::from_signs(s.x.signum() as i8, s.y.signum() as i8);

    // Component: for a torus there's a single component (0).  For a genus
    // surface we look up the component of the cell containing (u1, u2).
    let component = match level {
        Level::Torus { .. } => 0,
        Level::GenusSurface { region, .. } => {
            let (i, j) = cell_of(&region.u1, &region.u2, u1, u2);
            region.component[i][j].unwrap_or(0)
        }
        _ => 0,
    };

    Ok(FlatSample {
        u1,
        u2,
        sheet,
        fold,
        component,
    })
}

/// Find the cell `(i, j)` containing the flat point `(u1, u2)` given the
/// monotone grid lines.
fn cell_of(u1_lines: &[f32], u2_lines: &[f32], u1: f32, u2: f32) -> (usize, usize) {
    let i = grid_index(u1_lines, u1);
    let j = grid_index(u2_lines, u2);
    (i, j)
}

/// Index of the cell containing `x` in a monotone grid `lines` (length n+1,
/// cells 0..n).  Clamped to the valid range.
fn grid_index(lines: &[f32], x: f32) -> usize {
    if lines.len() < 2 {
        return 0;
    }
    let n = lines.len() - 1;
    let mut lo = 0usize;
    let mut hi = n;
    while lo < hi {
        let mid = (lo + hi) / 2;
        if x >= lines[mid + 1] {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    lo.min(n - 1)
}

/// The reflex-corner positions in flat `(u₁, u₂)` space for a level.  A segment
/// passing within `eps` of any of these points should be terminated/flagged
/// (the 3-prong separatrix rule).
pub fn reflex_corners_in_flat(level: &Level) -> Vec<(f32, f32)> {
    match level {
        Level::GenusSurface { region, .. } => region
            .reflex_by_component
            .iter()
            .flatten()
            .copied()
            .collect(),
        _ => vec![],
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pseudo::classify::classify_level;
    use crate::table::Table;

    fn cf() -> ConfocalParams {
        ConfocalParams::new(4.0, 1.0)
    }

    fn l_table() -> Table {
        // Doc's standard L: ell = [0, α₁, α₂], hyp = [β₁, β₂, β₃]
        //   base = [0, α₂] × [β₁, β₂],  tall = [0, α₁] × [β₂, β₃]
        Table::from_bounds(&[0.0, 0.5, 1.5], &[2.0, 2.5, 3.0], |l1, l2, _| {
            let tall = l1 <= 0.5 && l2 >= 2.5;
            let base = l1 <= 1.5 && l2 <= 2.5;
            tall || base
        })
    }

    #[test]
    fn test_to_flat_basic() {
        let tab = l_table();
        let level = classify_level(0.7, &tab, &cf(), 1e-9, 1e-9);
        // A sample inside the base region of the L.
        // (x, y) with λ₁ ≈ 0.7, λ₂ ≈ 2.2 — pick a point in the base.
        // We'll just verify the mapping runs and returns finite values.
        let s = PhaseSample::new(1.0, 0.5, 0.0, 1.0);
        let fs = to_flat(&s, &cf(), &level).expect("mapping should succeed");
        assert!(fs.u1.is_finite() && fs.u2.is_finite());
        assert!(fs.sheet.0 == 1 || fs.sheet.0 == -1);
        assert!(fs.sheet.1 == 1 || fs.sheet.1 == -1);
    }

    #[test]
    fn test_to_flat_forbidden() {
        let tab = l_table();
        let level = classify_level(3.5, &tab, &cf(), 1e-9, 1e-9);
        let s = PhaseSample::new(1.0, 0.5, 0.0, 1.0);
        assert!(matches!(
            to_flat(&s, &cf(), &level),
            Err(FlatError::Forbidden)
        ));
    }

    #[test]
    fn test_reflex_corners() {
        let tab = l_table();
        let level = classify_level(0.7, &tab, &cf(), 1e-9, 1e-9);
        let corners = reflex_corners_in_flat(&level);
        assert_eq!(corners.len(), 1, "L elliptic genus-2 has one reflex corner");
    }
}
