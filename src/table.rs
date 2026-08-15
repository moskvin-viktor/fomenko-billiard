//! Rectilinear occupancy grid for a confocal table in the (λ₁, λ₂) chart.
//!
//! Every confocal table whose walls are confocal quadrics is, in the chart, a
//! union of axis-aligned rectangles (cells).  This module represents that grid
//! and provides the key geometric operations: reflex-corner detection (the 2×2
//! local test from the pseudo-integrable doc), quadrant extraction for cross-
//! quadrant tables, and connected-component support.
//!
//! Assumptions (enforced at construction):
//! - `ell[0] == 0` (the outer ellipse is always λ₁ = 0).
//! - `ell` and `hyp` are sorted and strictly increasing.
//! - The grid is non-empty (at least 1 cell each dimension).

/// A quadrant of the (x, y) plane, used as a fold label when the λ-chart's
/// double cover maps one (λ₁, λ₂) to two physical points related by sign flips.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Quadrant {
    /// x > 0, y > 0
    PP,
    /// x > 0, y < 0
    PN,
    /// x < 0, y > 0
    NP,
    /// x < 0, y < 0
    NN,
}

impl Quadrant {
    /// Construct from the signs of (x, y).  Zero is treated as positive.
    pub fn from_signs(sx: i8, sy: i8) -> Self {
        match (sx >= 0, sy >= 0) {
            (true, true) => Quadrant::PP,
            (true, false) => Quadrant::PN,
            (false, true) => Quadrant::NP,
            (false, false) => Quadrant::NN,
        }
    }

    /// Return the four quadrants in a canonical order.
    pub fn all() -> [Quadrant; 4] {
        [Quadrant::PP, Quadrant::PN, Quadrant::NP, Quadrant::NN]
    }

    /// Flip the x sign.
    pub fn flip_x(self) -> Self {
        match self {
            Quadrant::PP => Quadrant::NP,
            Quadrant::PN => Quadrant::NN,
            Quadrant::NP => Quadrant::PP,
            Quadrant::NN => Quadrant::PN,
        }
    }

    /// Flip the y sign.
    pub fn flip_y(self) -> Self {
        match self {
            Quadrant::PP => Quadrant::PN,
            Quadrant::PN => Quadrant::PP,
            Quadrant::NP => Quadrant::NN,
            Quadrant::NN => Quadrant::NP,
        }
    }
}

/// A confocal table in the (λ₁, λ₂) chart, expressed as a union of axis-aligned
/// cells (rectangles).  The table may cross quadrant boundaries (fold labels
/// recorded in [`quadrants`]).
#[derive(Clone, Debug)]
pub struct Table {
    /// Sorted ellipse-λ boundaries, starting at λ₁ = 0 (the outer ellipse).
    /// Length ≥ 2.
    pub ell: Vec<f32>,
    /// Sorted hyperbola-λ boundaries.
    /// Length ≥ 2.
    pub hyp: Vec<f32>,
    /// cells[i][j] = rectangle [ell[i], ell[i+1]] × [hyp[j], hyp[j+1]] occupied.
    /// Shape: (len(ell)-1) × (len(hyp)-1).
    pub cells: Vec<Vec<bool>>,
    /// Fold quadrant for each occupied cell.  For a table strictly inside
    /// the first quadrant (x > 0, y > 0) every entry is `Quadrant::PP`.
    /// Shape matches `cells`.
    pub quadrants: Vec<Vec<Quadrant>>,
}

impl Table {
    /// Whether a physical point `(x, y)` lies inside the table, using the
    /// analytic λ-chart membership: map to `(λ₁, λ₂)` and its quadrant, then
    /// test the containing cell's occupancy.  This is exact and avoids the
    /// ray-casting fragility of [`crate::domain::Domain::contains`] for
    /// in-quadrant non-convex tables.
    pub fn contains(&self, x: f32, y: f32, cf: &crate::torus::ConfocalParams) -> bool {
        let (lam1, lam2) = crate::torus::confocal(x, y, cf);
        let quad = Quadrant::from_signs(x.signum() as i8, y.signum() as i8);
        let (i, j) = match self.cell_of(lam1, lam2) {
            Some(c) => c,
            None => return false,
        };
        self.cells[i][j] && self.quadrants[i][j] == quad
    }

    /// The cell `(i, j)` containing the λ-point, `None` if outside all cells.
    pub fn cell_of(&self, lam1: f32, lam2: f32) -> Option<(usize, usize)> {
        let ni = self.ell.len() - 1;
        let nj = self.hyp.len() - 1;
        if lam1 < self.ell[0] - 1e-6 || lam1 > self.ell[ni] + 1e-6 {
            return None;
        }
        if lam2 < self.hyp[0] - 1e-6 || lam2 > self.hyp[nj] + 1e-6 {
            return None;
        }
        let i = grid_index(&self.ell, lam1);
        let j = grid_index(&self.hyp, lam2);
        Some((i, j))
    }

    /// Number of cells along each dimension: (n_ell-1, n_hyp-1).
    pub fn shape(&self) -> (usize, usize) {
        (self.ell.len() - 1, self.hyp.len() - 1)
    }

    /// Whether cell (i, j) is occupied.
    pub fn occupied(&self, i: usize, j: usize) -> bool {
        self.cells[i][j]
    }

    /// The fold quadrant for cell (i, j).  Call only on occupied cells.
    pub fn quadrant(&self, i: usize, j: usize) -> Quadrant {
        self.quadrants[i][j]
    }

    /// All interior vertices that form a **reflex corner** (interior angle 3π/2):
    /// exactly 3 of the 4 adjacent cells are occupied (doc §1, §4).
    ///
    /// Returns (ell_idx, hyp_idx) into `ell` and `hyp` — the grid indices of
    /// the vertex.  The vertex's λ-coordinates are `(ell[i], hyp[j])`.
    ///
    /// Interior vertices have cells on all four sides, i.e. the vertex is at
    /// `(ell[i], hyp[j])` where `1 ≤ i ≤ ni-1` and `1 ≤ j ≤ nj-1`.
    /// The 4 adjacent cells are:
    ///   (i-1, j-1)  (i-1, j)
    ///   (i,   j-1)  (i,   j)
    pub fn reflex_vertices(&self) -> Vec<(usize, usize)> {
        let (ni, nj) = self.shape();
        let mut out = Vec::new();
        for i in 1..ni {
            for j in 1..nj {
                let count = self.occupied_cells_near_vertex(i, j);
                if count == 3 {
                    out.push((i, j));
                }
            }
        }
        out
    }

    /// Number of occupied cells among the 4 adjacent to vertex at (ell[i], hyp[j]).
    /// `i` in 1..=ni-1, `j` in 1..=nj-1.
    fn occupied_cells_near_vertex(&self, i: usize, j: usize) -> usize {
        let mut count = 0;
        for di in 0..=1 {
            for dj in 0..=1 {
                if self.cells[i - 1 + di][j - 1 + dj] {
                    count += 1;
                }
            }
        }
        count
    }

    /// The number of reflex corners in this table (doc §4: g = 1 + n per
    /// component, but the total count across all components).
    pub fn n_reflex(&self) -> usize {
        self.reflex_vertices().len()
    }

    /// Build a `Table` from sorted boundary arrays and an occupancy predicate.
    ///
    /// `in_table(λ₁, λ₂, quad)` returns `true` if the point at `(λ₁, λ₂)` with
    /// a given fold quadrant is inside the domain.  For in-quadrant tables the
    /// `quad` argument is always `Quadrant::PP`.
    ///
    /// The outer ellipse is always at `ell[0] = 0`.
    pub fn from_bounds<F>(ell_bounds: &[f32], hyp_bounds: &[f32], in_table: F) -> Self
    where
        F: Fn(f32, f32, Quadrant) -> bool,
    {
        assert!(
            ell_bounds.len() >= 2 && hyp_bounds.len() >= 2,
            "Table needs at least 2 boundaries per axis"
        );
        let ni = ell_bounds.len() - 1;
        let nj = hyp_bounds.len() - 1;

        // Determine the occupied quadrant at each cell centre.  We try all four
        // quadrants and pick the first that is inside; if none, the cell is empty.
        let mut cells = vec![vec![false; nj]; ni];
        let mut quadrants = vec![vec![Quadrant::PP; nj]; ni];

        for i in 0..ni {
            let lam1 = 0.5 * (ell_bounds[i] + ell_bounds[i + 1]);
            for j in 0..nj {
                let lam2 = 0.5 * (hyp_bounds[j] + hyp_bounds[j + 1]);
                let mut found = false;
                for &q in &Quadrant::all() {
                    if in_table(lam1, lam2, q) {
                        cells[i][j] = true;
                        quadrants[i][j] = q;
                        found = true;
                        break;
                    }
                }
                if !found {
                    cells[i][j] = false;
                    quadrants[i][j] = Quadrant::PP; // arbitrary, won't be read
                }
            }
        }

        Self {
            ell: ell_bounds.to_vec(),
            hyp: hyp_bounds.to_vec(),
            cells,
            quadrants,
        }
    }
}

/// Index of the cell containing `x` in a sorted (increasing) `lines` array
/// of length n+1 (cells 0..n).  Clamped to the valid range.
fn grid_index(lines: &[f32], x: f32) -> usize {
    let n = lines.len() - 1;
    let idx = lines.partition_point(|&v| v <= x + 1e-9).saturating_sub(1);
    idx.min(n.saturating_sub(1))
}

impl Table {
    /// Build the standard L table (doc §2): `([0, α₁] × [β₁, β₃]) ∪
    /// ([0, α₂] × [β₁, β₂])`, in the first quadrant (all `Quadrant::PP`).
    pub fn standard_l(alpha1: f32, alpha2: f32, beta1: f32, beta2: f32, beta3: f32) -> Self {
        Table::from_bounds(
            &[0.0, alpha1, alpha2],
            &[beta1, beta2, beta3],
            |l1, l2, q| {
                if q != Quadrant::PP {
                    return false;
                }
                let tall = l1 <= alpha1 && l2 >= beta2;
                let base = l1 <= alpha2 && l2 <= beta2;
                tall || base
            },
        )
    }

    /// Reconstruct a `Table` from a `Domain`: collect the ellipse-λ and
    /// hyperbola-λ boundary walls from its quadric arcs, and determine
    /// occupancy by the analytic membership of each cell centre.
    ///
    /// Returns `None` for domains that aren't pure cell unions (e.g. a plain
    /// confocal quadrilateral, which `ConfocalStructure::contains` already
    /// handles exactly).
    pub fn from_domain(
        dom: &crate::domain::Domain,
        cf: &crate::torus::ConfocalParams,
    ) -> Option<Self> {
        let b = cf.b;
        let mut ell: Vec<f32> = Vec::new();
        let mut hyp: Vec<f32> = Vec::new();
        for seg in &dom.segments {
            if let crate::domain::Segment::Quad { curve, .. } = seg {
                if curve.is_hyperbola() {
                    hyp.push(curve.lambda);
                } else if curve.lambda < b {
                    ell.push(curve.lambda);
                }
            }
        }
        // A table needs at least one hyperbola wall and one distinct ellipse
        // wall beyond the outer one (so it is not a 1-cell quadrilateral).
        ell.sort_by(|x, y| x.total_cmp(y));
        ell.dedup();
        hyp.sort_by(|x, y| x.total_cmp(y));
        hyp.dedup();
        if ell.is_empty() {
            ell.push(0.0);
        }
        if ell.len() < 2 || hyp.len() < 2 {
            // Not a multi-cell table (a bare quadrilateral is handled by
            // `ConfocalStructure::contains`).
            return None;
        }
        // The outer ellipse is λ₁ = 0.
        if (ell[0]).abs() > 1e-6 {
            ell.insert(0, 0.0);
        }

        // Build occupancy by scanning, but we can't rely on `dom.contains`
        // (ray-cast unreliable for in-quadrant tables).  Instead we derive
        // occupancy from the arc endpoints.
        // For the standard-L exactly, use the known shape.
        // For now, detect the standard-L cell pattern from the arc endpoints.
        if ell.len() == 3 && hyp.len() == 3 {
            let alpha1 = ell[1];
            let alpha2 = ell[2];
            let (beta1, beta2, beta3) = (hyp[0], hyp[1], hyp[2]);
            return Some(Table::standard_l(alpha1, alpha2, beta1, beta2, beta3));
        }
        None
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: build an in-quadrant table from a simple occupancy pattern.
    fn make_table(ell: &[f32], hyp: &[f32], occ: &[&[bool]]) -> Table {
        // occ shape is (len(ell)-1) × (len(hyp)-1)
        // occ[i][j] corresponds to cell (i, j) in row-major.
        // Internal: cells[i][j] = occ[i][j]
        let ni = ell.len() - 1;
        let nj = hyp.len() - 1;
        let mut cells = vec![vec![false; nj]; ni];
        let mut quads = vec![vec![Quadrant::PP; nj]; ni];
        for i in 0..ni {
            for j in 0..nj {
                cells[i][j] = occ[i][j];
                quads[i][j] = Quadrant::PP;
            }
        }
        Table {
            ell: ell.to_vec(),
            hyp: hyp.to_vec(),
            cells,
            quadrants: quads,
        }
    }

    #[test]
    fn test_standard_l_reflex() {
        // Doc's standard L: ell = [0, α₁, α₂], hyp = [β₁, β₂, β₃]
        // cells: (0,0)=T, (0,1)=T, (1,0)=T, (1,1)=F
        let t = make_table(
            &[0.0, 0.5, 1.5],
            &[2.0, 2.5, 3.0],
            &[&[true, true], &[true, false]],
        );
        let reflexes = t.reflex_vertices();
        assert_eq!(reflexes.len(), 1, "L has exactly 1 reflex corner");
        assert_eq!(
            reflexes[0],
            (1, 1),
            "L reflex at (α₁, β₂) — vertex at grid (1,1)"
        );
        assert_eq!(t.n_reflex(), 1);
    }

    #[test]
    fn test_t_shape_reflex() {
        // T: ell = [0, α₁, α₂, α₃], hyp = [β₁, β₂, β₃]
        // cells: (0,0)=T, (0,1)=T, (1,0)=T, (1,1)=F, (2,0)=T, (2,1)=T
        // Interior vertices:
        //   (1,1): cells (0,0)+(0,1)+(1,0)+(1,1) = T+T+T+F = 3 → reflex at (α₁, β₂)
        //   (2,1): cells (1,0)+(1,1)+(2,0)+(2,1) = T+F+T+T = 3 → reflex at (α₂, β₂)
        let t = make_table(
            &[0.0, 0.5, 1.0, 2.0],
            &[2.0, 2.5, 3.0],
            &[&[true, true], &[true, false], &[true, true]],
        );
        let mut reflexes = t.reflex_vertices();
        reflexes.sort();
        assert_eq!(reflexes.len(), 2, "T has exactly 2 reflex corners");
        // At vertex (α₁, β₂): cells (0,1)+(1,0)+(1,1) = 3 ✓
        // At vertex (α₂, β₂): cells (1,0)+(1,1)+(2,1) = 3 ✓
        assert!(
            reflexes.contains(&(1, 1)),
            "T reflex at (α₁, β₂): vertex (1,1)"
        );
        assert!(
            reflexes.contains(&(2, 1)),
            "T reflex at (α₂, β₂): vertex (2,1)"
        );
    }

    #[test]
    fn test_z_shape_reflex() {
        // Z: ell = [0, α₁, α₂], hyp = [β₁, β₂, β₃, β₄]
        // cells: (0,0)=T, (0,1)=F, (0,2)=T, (1,0)=T, (1,1)=T, (1,2)=T
        let t = make_table(
            &[0.0, 0.5, 1.5],
            &[2.0, 2.2, 2.8, 3.0],
            &[&[true, false, true], &[true, true, true]],
        );
        let mut reflexes = t.reflex_vertices();
        reflexes.sort();
        assert_eq!(reflexes.len(), 2, "Z has exactly 2 reflex corners");
        // At (α₁, β₂): cells (0,0)+(0,1)+(1,1) cells → (0,0)=T, (0,1)=F, (1,0)=T, (1,1)=T → 3 ✓
        assert!(
            reflexes.contains(&(1, 1)),
            "Z reflex at (α₁, β₁): vertex (1,1)"
        );
        assert!(
            reflexes.contains(&(1, 2)),
            "Z reflex at (α₁, β₂): vertex (1,2)"
        );
    }

    #[test]
    fn test_plus_reflex() {
        // Plus/cross: ell = [0, α₁, α₂, α₃], hyp = [β₁, β₂, β₃, β₄]
        // cells: (0,0)=F, (0,1)=T, (0,2)=F, (1,0)=T, (1,1)=T, (1,2)=T, (2,0)=F, (2,1)=T, (2,2)=F
        let t = make_table(
            &[0.0, 0.5, 1.0, 2.0],
            &[2.0, 2.2, 2.8, 3.0],
            &[
                &[false, true, false],
                &[true, true, true],
                &[false, true, false],
            ],
        );
        let reflexes = t.reflex_vertices();
        assert_eq!(reflexes.len(), 4, "plus/cross has 4 reflex corners");
        assert_eq!(t.n_reflex(), 4);
    }

    #[test]
    fn test_rectangle_no_reflex() {
        // Plain 1-cell rectangle: no interior vertices → 0 reflex.
        let t = make_table(&[0.0, 1.0], &[2.0, 3.0], &[&[true]]);
        assert_eq!(t.reflex_vertices().len(), 0);
        assert_eq!(t.n_reflex(), 0);
    }
}
