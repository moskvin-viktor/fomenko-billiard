//! Layers 2–3 of the molecule engine: the accessible-region mask and its
//! topology.
//!
//! At a caustic level `λc` the accessible region is the table intersected with
//! the caustic constraint (elliptic `λc < b` keeps `λ₁ ≤ λc`; hyperbolic
//! `λc > b` keeps `λ₂ ≥ λc`).  The cut line is inserted into the grid so no
//! cell straddles it, and every topological fact — components, corners, genus
//! `= 1 + n_reflex` — is read off the mask by inspecting the four cells around
//! each grid vertex.

use crate::table::Table;
use crate::torus::ConfocalParams;

/// A grid of cells in the `(λ₁, λ₂)` chart at a fixed caustic level, with the
/// caustic cut line inserted.
#[derive(Clone, Debug)]
pub struct Grid {
    /// Sorted `λ₁` edge coordinates (includes the caustic cut when elliptic).
    pub xs: Vec<f32>,
    /// Sorted `λ₂` edge coordinates (includes the caustic cut when hyperbolic).
    pub ys: Vec<f32>,
    /// `inside[i][j]` = cell `(i, j)` is accessible.
    /// Shape `(len(xs)-1) × (len(ys)-1)`.
    pub inside: Vec<Vec<bool>>,
}

impl Grid {
    /// Number of cells along each axis: `(len(xs)-1, len(ys)-1)`.
    pub fn shape(&self) -> (usize, usize) {
        (self.xs.len() - 1, self.ys.len() - 1)
    }
}

/// Build the accessible-region mask at caustic level `lc`.
///
/// Elliptic caustic (`lc < b`) keeps `λ₁ ≤ lc`; hyperbolic (`lc > b`) keeps
/// `λ₂ ≥ lc`.  The cut line is inserted into the grid, so the cell-centre test
/// is exact (no cell straddles the cut).
pub fn accessible_grid(table: &Table, cf: &ConfocalParams, lc: f32) -> Grid {
    let b = cf.b;
    let xs = if lc < b {
        insert_sorted(&table.ell, lc)
    } else {
        table.ell.clone()
    };
    let ys = if lc > b {
        insert_sorted(&table.hyp, lc)
    } else {
        table.hyp.clone()
    };
    let nx = xs.len() - 1;
    let ny = ys.len() - 1;
    let mut inside = vec![vec![false; ny]; nx];
    for i in 0..nx {
        let cx = 0.5 * (xs[i] + xs[i + 1]);
        for j in 0..ny {
            let cy = 0.5 * (ys[j] + ys[j + 1]);
            let Some((oi, oj)) = table.cell_of(cx, cy) else {
                continue;
            };
            if !table.occupied(oi, oj) {
                continue;
            }
            let keep = if lc < b { cx <= lc } else { cy >= lc };
            inside[i][j] = keep;
        }
    }
    Grid { xs, ys, inside }
}

/// The kind of a grid vertex, from the number/arrangement of its four
/// surrounding cells (doc §5).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VertexType {
    /// Exactly 1 of 4 cells inside: convex corner, interior angle `π/2`.
    Convex,
    /// Exactly 3 of 4 inside: reflex corner, angle `3π/2` — genus-carrying.
    Reflex,
    /// 2 diagonal cells inside: two parts kiss at a point (region splitting).
    Pinch,
    /// 2 adjacent cells inside: straight edge, no corner.
    Edge,
    /// All 4 cells inside: not a boundary vertex.
    Interior,
    /// No cells inside: not a boundary vertex.
    Exterior,
}

/// The four cells around grid vertex `(vi, vj)`, as `[SW, SE, NW, NE]`.
fn surrounding(grid: &Grid, vi: usize, vj: usize) -> [bool; 4] {
    let (nx, ny) = grid.shape();
    let mut out = [false; 4];
    for (k, (di, dj)) in [(-1i32, -1i32), (0, -1), (-1, 0), (0, 0)]
        .iter()
        .enumerate()
    {
        let ci = vi as i32 + di;
        let cj = vj as i32 + dj;
        if ci >= 0 && cj >= 0 && ci < nx as i32 && cj < ny as i32 {
            out[k] = grid.inside[ci as usize][cj as usize];
        }
    }
    out
}

/// Classify a grid vertex by its four surrounding cells (doc §5).
pub fn classify_vertex(grid: &Grid, vi: usize, vj: usize) -> VertexType {
    let s = surrounding(grid, vi, vj);
    let c = s.iter().filter(|&&b| b).count();
    match c {
        1 => VertexType::Convex,
        3 => VertexType::Reflex,
        2 => {
            let (sw, se, nw, ne) = (s[0], s[1], s[2], s[3]);
            if (sw && ne) || (se && nw) {
                VertexType::Pinch
            } else {
                VertexType::Edge
            }
        }
        4 => VertexType::Interior,
        _ => VertexType::Exterior,
    }
}

/// 4-connected flood fill of the accessible cells (a diagonal kiss does NOT
/// connect).  Returns a label per cell (`None` for inaccessible) and the number
/// of components.
pub fn components(grid: &Grid) -> (Vec<Vec<Option<u32>>>, u32) {
    let (nx, ny) = grid.shape();
    let mut label = vec![vec![None; ny]; nx];
    let mut next = 0u32;
    for i in 0..nx {
        for j in 0..ny {
            if !grid.inside[i][j] || label[i][j].is_some() {
                continue;
            }
            let mut stack = vec![(i, j)];
            label[i][j] = Some(next);
            while let Some((ci, cj)) = stack.pop() {
                for (di, dj) in [(1i32, 0i32), (-1, 0), (0, 1), (0, -1)] {
                    let ni2 = ci as i32 + di;
                    let nj2 = cj as i32 + dj;
                    if ni2 < 0 || nj2 < 0 || ni2 >= nx as i32 || nj2 >= ny as i32 {
                        continue;
                    }
                    let (ni2, nj2) = (ni2 as usize, nj2 as usize);
                    if grid.inside[ni2][nj2] && label[ni2][nj2].is_none() {
                        label[ni2][nj2] = Some(next);
                        stack.push((ni2, nj2));
                    }
                }
            }
            next += 1;
        }
    }
    (label, next)
}

/// The component of a reflex vertex: the label of one of its in-cells.
pub fn component_of_vertex(
    grid: &Grid,
    label: &[Vec<Option<u32>>],
    vi: usize,
    vj: usize,
) -> Option<u32> {
    let (nx, ny) = grid.shape();
    for (di, dj) in [(-1i32, -1i32), (0, -1), (-1, 0), (0, 0)] {
        let ci = vi as i32 + di;
        let cj = vj as i32 + dj;
        if ci >= 0 && cj >= 0 && ci < nx as i32 && cj < ny as i32 {
            if let Some(c) = label[ci as usize][cj as usize] {
                return Some(c);
            }
        }
    }
    None
}

/// The topology of an accessible grid: components, genus per component
/// (`= 1 + n_reflex`), and corner/pinch counts.
#[derive(Clone, Debug)]
pub struct Topology {
    /// Number of connected components.
    pub n_components: usize,
    /// Genus of each component, `1 + n_reflex` (doc §5).
    pub genus: Vec<u32>,
    /// Number of convex corners.
    pub convex: usize,
    /// Number of reflex corners.
    pub reflex: usize,
    /// Number of pinch vertices.
    pub pinch: usize,
    /// Whether nothing is accessible.
    pub empty: bool,
}

/// Read the topology off a mask (doc §5).
pub fn topology(grid: &Grid) -> Topology {
    let (nx, ny) = grid.shape();
    let (label, ncomp) = components(grid);
    let mut reflex_per_comp = vec![0usize; ncomp as usize];
    let mut convex = 0usize;
    let mut reflex = 0usize;
    let mut pinch = 0usize;
    for vi in 0..=nx {
        for vj in 0..=ny {
            match classify_vertex(grid, vi, vj) {
                VertexType::Convex => convex += 1,
                VertexType::Reflex => {
                    reflex += 1;
                    if let Some(c) = component_of_vertex(grid, &label, vi, vj) {
                        reflex_per_comp[c as usize] += 1;
                    }
                }
                VertexType::Pinch => pinch += 1,
                _ => {}
            }
        }
    }
    let genus: Vec<u32> = reflex_per_comp.iter().map(|&r| 1 + r as u32).collect();
    Topology {
        n_components: ncomp as usize,
        genus,
        convex,
        reflex,
        pinch,
        empty: ncomp == 0,
    }
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
