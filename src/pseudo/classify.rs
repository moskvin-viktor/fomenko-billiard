//! Level classification for pseudo-integrable confocal tables.
//!
//! Given a caustic value `λc` and a [`Table`], determine the topology of the
//! accessible region and the flat moduli needed to map points into it.  This is
//! the generalized version of the doc's `classify_level` (§11 of
//! `docs/math/pseudo_integrable_genus.md`): instead of hardcoding the L-shape's
//! branching, it truncates the occupancy grid by the caustic, labels connected
//! components, and counts surviving reflex corners per component.

use super::ulength::{RootSide, ULength};
use crate::table::{Quadrant, Table};
use crate::torus::ConfocalParams;

/// The accessible flat region for a genus-≥2 level, after caustic truncation.
#[derive(Clone, Debug)]
pub struct AccessibleRegion {
    /// The `u₁` grid lines (monotone): `[u₁(ell[0]), u₁(ell[1]), …]`.
    pub u1: Vec<f32>,
    /// The `u₂` grid lines (monotone): `[u₂(hyp[0]), u₂(hyp[1]), …]`.
    pub u2: Vec<f32>,
    /// accessible[i][j] = cell (i, j) is reachable (occupied and not forbidden
    /// by the caustic cut).
    pub accessible: Vec<Vec<bool>>,
    /// Connected-component label per accessible cell (4-neighbour, respecting
    /// fold quadrants).  `None` for inaccessible cells.
    pub component: Vec<Vec<Option<u32>>>,
    /// Number of connected components.
    pub n_components: u32,
    /// Flat positions `(u₁, u₂)` of reflex corners that survive the cut,
    /// grouped by component index.
    pub reflex_by_component: Vec<Vec<(f32, f32)>>,
    /// Genus of each component = 1 + n_reflex_in_component.
    pub genus: Vec<u32>,
}

/// The topology of a level set.
#[derive(Clone, Debug)]
pub enum Level {
    /// No point of the table is accessible (λc beyond the outer hyperbola wall).
    Forbidden,
    /// λc = b: the caustic degenerates to the focal segment; periods diverge.
    Separatrix,
    /// The accessible region is a single rectangle → a torus.  Carries the two
    /// `ULength`s and the rectangle's flat size `(w1, w2)`.
    Torus {
        u1: ULength,
        u2: ULength,
        shape: (f32, f32),
    },
    /// The accessible region is not a rectangle → a genus-≥2 surface.
    GenusSurface {
        u1: ULength,
        u2: ULength,
        region: AccessibleRegion,
    },
}

/// Classify the level `λc` for a table.
///
/// * `sep_eps` — tolerance for `|λc − b|` (separatrix).
/// * `cut_eps` — tolerance for "a grid line equals λc" (bifurcation).
pub fn classify_level(
    lc: f32,
    tab: &Table,
    cf: &ConfocalParams,
    sep_eps: f32,
    cut_eps: f32,
) -> Level {
    let a = cf.a;
    let b = cf.b;

    // Separatrix — only a true singular level when the table actually touches
    // the focal segment.  Otherwise λc = b is a regular level of the flat
    // family (the caustic degenerates outside the table; every used interval
    // stays a positive distance from b) and must classify normally.
    if (lc - b).abs() < sep_eps && table_touches_focal(tab, cf) {
        return Level::Separatrix;
    }

    // Forbidden: λc beyond the outermost hyperbola wall.
    let hyp_max = *tab.hyp.last().unwrap();
    if lc >= hyp_max {
        return Level::Forbidden;
    }

    // Near λc = b the cubic P has a (near-)double root and the s-substitution
    // in `ULength` removes only one of the two coincident roots, so the raw
    // quadrature degenerates.  On a non-focal table (the only way we get here
    // with |λc − b| small) every *used* λ-interval stays a positive distance
    // from b, so build the ULengths at a λc nudged off b — error O(ε) — while
    // the cut / occupancy logic below keeps the true λc.
    let lc_u = if (lc - b).abs() < 1e-4 {
        if lc < b {
            b - 1e-4
        } else {
            b + 1e-4
        }
    } else {
        lc
    };

    // Build the two ULengths over the full accessible extent.
    //   elliptic (lc < b): λ₁ ∈ [0, lc] (turning at hi), λ₂ ∈ [hyp[0], hyp_max]
    //                      (both walls, root above = a).
    //   hyperbolic (lc > b): λ₁ ∈ [0, ell[ni]] (both walls), λ₂ ∈ [lc, hyp_max]
    //                      (turning at lo).
    let (u1, u2) = if lc < b {
        (
            ULength::new(0.0, lc_u, a, b, lc_u, RootSide::Hi, 512),
            ULength::new(tab.hyp[0], hyp_max, a, b, lc_u, RootSide::Hi, 512),
        )
    } else {
        (
            ULength::new(0.0, *tab.ell.last().unwrap(), a, b, lc_u, RootSide::Hi, 512),
            ULength::new(lc_u, hyp_max, a, b, lc_u, RootSide::Lo, 512),
        )
    };

    // Refined grid: insert λc as a cut line so partial cells become full cells.
    //   elliptic → cut in λ₁ (keep λ₁ ≤ lc); hyperbolic → cut in λ₂ (keep λ₂ ≥ lc).
    //
    // Only insert the cut when it lies strictly inside the table's extent:
    // outside it the cut cannot truncate any occupied cell, and inserting it
    // would create phantom grid cells beyond the table (whose midpoints used
    // to clamp onto boundary cells and read as accessible).
    let (ell_ref, hyp_ref) = if lc < b {
        let ell = if lc < *tab.ell.last().unwrap() {
            insert_sorted(&tab.ell, lc)
        } else {
            tab.ell.clone()
        };
        (ell, tab.hyp.clone())
    } else {
        let hyp = if lc > tab.hyp[0] {
            insert_sorted(&tab.hyp, lc)
        } else {
            tab.hyp.clone()
        };
        (tab.ell.clone(), hyp)
    };
    let ni = ell_ref.len() - 1;
    let nj = hyp_ref.len() - 1;

    // Occupancy + fold quadrant of each refined cell: it must lie inside an
    // occupied original cell and on the accessible side of the cut.
    let mut accessible = vec![vec![false; nj]; ni];
    let mut ref_quadrants = vec![vec![Quadrant::PP; nj]; ni];
    for i in 0..ni {
        let lam1_lo = ell_ref[i];
        let lam1_hi = ell_ref[i + 1];
        for j in 0..nj {
            let lam2_lo = hyp_ref[j];
            let lam2_hi = hyp_ref[j + 1];
            // The containing original cell; outside the table extent → unoccupied.
            let Some((oi, oj)) =
                containing_cell(tab, 0.5 * (lam1_lo + lam1_hi), 0.5 * (lam2_lo + lam2_hi))
            else {
                continue;
            };
            if !tab.occupied(oi, oj) {
                continue;
            }
            // Caustic bound (redundant after the cut, kept for safety).
            let keep = if lc < b {
                lam1_hi <= lc + cut_eps
            } else {
                lam2_lo >= lc - cut_eps
            };
            accessible[i][j] = keep;
            ref_quadrants[i][j] = tab.quadrant(oi, oj);
        }
    }

    // If nothing is accessible, it's forbidden (or a degenerate sliver).
    let any = accessible.iter().any(|row| row.iter().any(|&c| c));
    if !any {
        return Level::Forbidden;
    }

    // If the accessible cells form a single rectangle (contiguous, no holes,
    // single fold quadrant) → torus.
    if let Some((min_i, max_i, min_j, max_j)) = as_single_rectangle(&accessible, &ref_quadrants) {
        let w1 = u1.u(ell_ref[max_i + 1]) - u1.u(ell_ref[min_i]);
        let w2 = u2.u(hyp_ref[max_j + 1]) - u2.u(hyp_ref[min_j]);
        return Level::Torus {
            u1,
            u2,
            shape: (w1, w2),
        };
    }

    // Genus surface: label connected components (respecting fold quadrants),
    // then count surviving reflex corners per component.
    let (component, n_components) = label_components(&accessible, &ref_quadrants);

    // Flat grid lines (refined).
    let u1_lines: Vec<f32> = ell_ref.iter().map(|&e| u1.u(e)).collect();
    let u2_lines: Vec<f32> = hyp_ref.iter().map(|&h| u2.u(h)).collect();

    // Count reflex corners per component: a reflex vertex (i, j) survives if
    // its 4 adjacent cells are all accessible (the corner is interior to the
    // accessible region).  It belongs to the component of those cells.
    //
    // NOTE: reflex vertices are defined on the ORIGINAL grid; after refinement
    // the cut may have moved/split them.  We only count a reflex corner if the
    // vertex still has 3-of-4 occupied cells in the refined grid.
    let mut reflex_by_component: Vec<Vec<(f32, f32)>> = vec![Vec::new(); n_components as usize];
    for (i, j) in tab.reflex_vertices() {
        // The vertex (ell[i], hyp[j]) — locate it in the refined grid.
        let Some(ri) = find_in(&ell_ref, tab.ell[i]) else {
            continue;
        };
        let Some(rj) = find_in(&hyp_ref, tab.hyp[j]) else {
            continue;
        };
        // Its 4 adjacent refined cells: (ri-1..=ri, rj-1..=rj).  A reflex corner
        // of the accessible region has exactly 3 of the 4 occupied (accessible);
        // if the cut shadows it (≤2 accessible) it contributes nothing.
        if ri == 0 || rj == 0 || ri >= ni || rj >= nj {
            continue;
        }
        let cells = [(ri - 1, rj - 1), (ri - 1, rj), (ri, rj - 1), (ri, rj)];
        let n_accessible = cells.iter().filter(|&&(ci, cj)| accessible[ci][cj]).count();
        if n_accessible != 3 {
            continue;
        }
        let comp = component[ri - 1][rj - 1].unwrap();
        reflex_by_component[comp as usize].push((u1_lines[ri], u2_lines[rj]));
    }

    let genus: Vec<u32> = reflex_by_component
        .iter()
        .map(|r| 1 + r.len() as u32)
        .collect();

    let region = AccessibleRegion {
        u1: u1_lines,
        u2: u2_lines,
        accessible,
        component,
        n_components,
        reflex_by_component,
        genus,
    };

    Level::GenusSurface { u1, u2, region }
}

/// Insert `x` into a sorted `Vec<f32>`, keeping it sorted and deduplicated.
fn insert_sorted(v: &[f32], x: f32) -> Vec<f32> {
    let mut out = v.to_vec();
    if out.iter().any(|&e| (e - x).abs() < 1e-9) {
        return out;
    }
    let idx = out.partition_point(|&e| e < x);
    out.insert(idx, x);
    out
}

/// Does the table touch the focal segment (λ = b)?
///
/// Only then is λc = b a true separatrix (diverging periods, singular level).
/// A table whose grid stays strictly inside `(b, ·)` hyperbolically and
/// `(·, b)` elliptically never sees the degeneracy: λc = b is a regular level.
pub fn table_touches_focal(tab: &Table, cf: &ConfocalParams) -> bool {
    const EPS: f32 = 1e-6;
    tab.hyp[0] < cf.b + EPS || *tab.ell.last().unwrap() > cf.b - EPS
}

/// Find the original cell `(oi, oj)` containing the point `(λ₁, λ₂)`, `None`
/// if the point lies outside the table's grid extent.
fn containing_cell(tab: &Table, lam1: f32, lam2: f32) -> Option<(usize, usize)> {
    if lam1 < tab.ell[0]
        || lam1 > *tab.ell.last().unwrap()
        || lam2 < tab.hyp[0]
        || lam2 > *tab.hyp.last().unwrap()
    {
        return None;
    }
    let oi = tab
        .ell
        .partition_point(|&e| e <= lam1)
        .saturating_sub(1)
        .min(tab.ell.len() - 2);
    let oj = tab
        .hyp
        .partition_point(|&e| e <= lam2)
        .saturating_sub(1)
        .min(tab.hyp.len() - 2);
    Some((oi, oj))
}

/// Find the index of `x` in a sorted `Vec<f32>` (exact match within tolerance).
fn find_in(v: &[f32], x: f32) -> Option<usize> {
    v.iter().position(|&e| (e - x).abs() < 1e-9)
}

/// If the accessible cells form a single rectangle (contiguous block, no holes,
/// single fold quadrant), return its bounding cell indices
/// `(min_i, max_i, min_j, max_j)`.
fn as_single_rectangle(
    accessible: &[Vec<bool>],
    quadrants: &[Vec<Quadrant>],
) -> Option<(usize, usize, usize, usize)> {
    let (ni, nj) = (accessible.len(), accessible[0].len());
    let mut min_i = ni;
    let mut max_i = 0usize;
    let mut min_j = nj;
    let mut max_j = 0usize;
    let mut count = 0;
    let mut quad: Option<Quadrant> = None;
    for i in 0..ni {
        for j in 0..nj {
            if accessible[i][j] {
                count += 1;
                min_i = min_i.min(i);
                max_i = max_i.max(i);
                min_j = min_j.min(j);
                max_j = max_j.max(j);
                let q = quadrants[i][j];
                quad = Some(match quad {
                    None => q,
                    Some(prev) if prev == q => q,
                    Some(_) => return None, // multiple quadrants → not a rectangle
                });
            }
        }
    }
    if count == 0 {
        return None;
    }
    let (w, h) = (max_i - min_i + 1, max_j - min_j + 1);
    if count != w * h {
        return None; // has holes or is disconnected
    }
    Some((min_i, max_i, min_j, max_j))
}

/// Label connected components of the accessible cells (4-neighbour), respecting
/// fold quadrants (cells in different quadrants are never adjacent).
fn label_components(
    accessible: &[Vec<bool>],
    quadrants: &[Vec<Quadrant>],
) -> (Vec<Vec<Option<u32>>>, u32) {
    let (ni, nj) = (accessible.len(), accessible[0].len());
    let mut comp = vec![vec![None; nj]; ni];
    let mut next = 0u32;
    for i in 0..ni {
        for j in 0..nj {
            if !accessible[i][j] || comp[i][j].is_some() {
                continue;
            }
            let mut stack = vec![(i, j)];
            comp[i][j] = Some(next);
            while let Some((ci, cj)) = stack.pop() {
                let q = quadrants[ci][cj];
                for (di, dj) in [(1i32, 0i32), (-1, 0), (0, 1), (0, -1)] {
                    let ni2 = ci as i32 + di;
                    let nj2 = cj as i32 + dj;
                    if ni2 < 0 || nj2 < 0 || ni2 >= ni as i32 || nj2 >= nj as i32 {
                        continue;
                    }
                    let (ni2, nj2) = (ni2 as usize, nj2 as usize);
                    if accessible[ni2][nj2] && comp[ni2][nj2].is_none() && quadrants[ni2][nj2] == q
                    {
                        comp[ni2][nj2] = Some(next);
                        stack.push((ni2, nj2));
                    }
                }
            }
            next += 1;
        }
    }
    (comp, next)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

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
    fn test_l_elliptic_torus() {
        // λc < α₁ → thin vertical strip → torus.
        let tab = l_table();
        let level = classify_level(0.3, &tab, &cf(), 1e-9, 1e-9);
        assert!(
            matches!(level, Level::Torus { .. }),
            "λc<α₁ should be a torus"
        );
    }

    #[test]
    fn test_l_elliptic_genus2() {
        // α₁ < λc < b → L-shape → genus 2.
        let tab = l_table();
        let level = classify_level(0.7, &tab, &cf(), 1e-9, 1e-9);
        match level {
            Level::GenusSurface { region, .. } => {
                assert_eq!(region.genus, vec![2], "L elliptic genus should be 2");
                assert_eq!(region.n_components, 1);
                assert_eq!(region.reflex_by_component[0].len(), 1);
            }
            other => panic!("expected genus surface, got {:?}", other),
        }
    }

    #[test]
    fn test_l_hyperbolic_torus() {
        // β₂ < λc < β₃ → thin horizontal strip → torus.
        let tab = l_table();
        let level = classify_level(2.7, &tab, &cf(), 1e-9, 1e-9);
        assert!(
            matches!(level, Level::Torus { .. }),
            "β₂<λc<β₃ should be a torus"
        );
    }

    #[test]
    fn test_l_hyperbolic_genus2() {
        // b < λc < β₁ → whole L → genus 2.
        let tab = l_table();
        let level = classify_level(1.2, &tab, &cf(), 1e-9, 1e-9);
        match level {
            Level::GenusSurface { region, .. } => {
                assert_eq!(region.genus, vec![2], "L hyperbolic genus should be 2");
            }
            other => panic!("expected genus surface, got {:?}", other),
        }
    }

    #[test]
    fn test_l_forbidden() {
        let tab = l_table();
        let level = classify_level(3.5, &tab, &cf(), 1e-9, 1e-9);
        assert!(
            matches!(level, Level::Forbidden),
            "λc ≥ β₃ should be forbidden"
        );
    }

    #[test]
    fn test_l_separatrix() {
        // The doc L touches the focal segment (ell_last = 1.5 > b = 1), so
        // λc = b really is a separatrix there.
        let tab = l_table();
        assert!(table_touches_focal(&tab, &cf()));
        let level = classify_level(1.0, &tab, &cf(), 1e-9, 1e-9);
        assert!(
            matches!(level, Level::Separatrix),
            "λc = b should be separatrix on a focal-touching table"
        );
    }

    /// The app's standard L (ell up to 0.8 < b, hyp from 1.4 > b) never touches
    /// the focal segment: λc = b is a *regular* genus-2 level, not a separatrix.
    fn standard_l_table() -> Table {
        Table::from_bounds(&[0.0, 0.4, 0.8], &[1.4, 2.0, 2.6], |l1, l2, _| {
            let tall = l1 <= 0.4 && l2 >= 2.0;
            let base = l1 <= 0.8 && l2 <= 2.0;
            tall || base
        })
    }

    #[test]
    fn test_standard_l_b_is_regular_genus2() {
        let tab = standard_l_table();
        assert!(!table_touches_focal(&tab, &cf()));
        for lc in [1.0f32, 1.0 - 5e-5, 1.0 + 5e-5] {
            let level = classify_level(lc, &tab, &cf(), 1e-9, 1e-9);
            match level {
                Level::GenusSurface { region, u1, u2 } => {
                    assert_eq!(region.genus, vec![2], "λc={lc}: genus should be 2");
                    assert!(u1.total.is_finite() && u1.total > 0.0);
                    assert!(u2.total.is_finite() && u2.total > 0.0);
                    for &v in region.u1.iter().chain(region.u2.iter()) {
                        assert!(v.is_finite(), "λc={lc}: non-finite grid line");
                    }
                }
                other => panic!("λc={lc}: expected genus surface, got {other:?}"),
            }
        }
    }

    /// Phantom-cell regression: for λc ∈ (α₂, b) no cut is inserted beyond the
    /// table's elliptic extent, so the flat region's u₁ extent must stay
    /// bounded by u₁(α₂) as λc → b⁻ (previously fake cells [α₂, λc] appeared).
    #[test]
    fn test_no_phantom_cells_near_b() {
        let tab = standard_l_table();
        let cf = cf();
        for lc in [0.9f32, 0.99, 0.999] {
            let level = classify_level(lc, &tab, &cf, 1e-9, 1e-9);
            match level {
                Level::GenusSurface { region, u1, .. } => {
                    let u1_a2 = u1.u(0.8);
                    let last = *region.u1.last().unwrap();
                    assert!(
                        last <= u1_a2 + 1e-4,
                        "λc={lc}: u1 extent {last} exceeds u1(α₂)={u1_a2} (phantom cells)"
                    );
                }
                other => panic!("λc={lc}: expected genus surface, got {other:?}"),
            }
        }
    }

    /// Exact grid-edge classification stays finite and non-forbidden.
    #[test]
    fn test_standard_l_exact_edges() {
        let tab = standard_l_table();
        let cf = cf();
        for lc in [0.4f32, 0.8, 1.4, 2.0] {
            let level = classify_level(lc, &tab, &cf, 1e-9, 1e-9);
            assert!(
                !matches!(level, Level::Forbidden | Level::Separatrix),
                "λc={lc}: exact edge must classify as a real level, got {level:?}"
            );
        }
    }
}
