//! Layer 6 of the molecule engine: flat geometry and phase-manifold evolution.
//!
//! To *show* the evolution, emit at each sampled `λc` the accessible region as
//! a flat polyomino in real `(u₁, u₂)` coordinates — the actual geometry of the
//! phase manifold at that level — together with its topology.  Stepping `λc`
//! animates the fibers deforming and bifurcating.

use crate::pseudo::ulength::{RootSide, ULength};
use crate::table::Table;
use crate::torus::ConfocalParams;

use super::critical::critical_values;
use super::grid::{accessible_grid, topology};

/// A snapshot of the phase manifold at a caustic level.
#[derive(Clone, Debug)]
pub struct FiberFrame {
    pub lc: f32,
    pub n_components: usize,
    /// Sorted genus per component.
    pub genus: Vec<u32>,
    /// Accessible cells mapped to `(u₁, u₂)` boxes.
    pub flat_rects: Vec<(f32, f32, f32, f32)>,
    /// Whether this frame sits at a critical value (a bifurcation).
    pub is_critical: bool,
}

/// Map each accessible cell to a box in `(u₁, u₂)` (doc §8).
///
/// Turning-line edges (where a cell's far edge equals the caustic cut) are
/// singularity-safe via the `ULength` `s = √(r − λ)` substitution.
pub fn flat_rects(table: &Table, cf: &ConfocalParams, lc: f32) -> Vec<(f32, f32, f32, f32)> {
    let a = cf.a;
    let b = cf.b;
    let lam1_lo = table.ell[0];
    let lam1_hi = if lc < b {
        lc
    } else {
        *table.ell.last().unwrap()
    };
    let lam2_lo = if lc > b { lc } else { table.hyp[0] };
    let lam2_hi = *table.hyp.last().unwrap();

    let u1 = ULength::new(lam1_lo, lam1_hi, a, b, lc, RootSide::Hi, 512);
    let u2 = ULength::new(
        lam2_lo,
        lam2_hi,
        a,
        b,
        lc,
        if lc > b { RootSide::Lo } else { RootSide::Hi },
        512,
    );

    let mut boxes = Vec::new();
    for i in 0..table.ell.len() - 1 {
        let (l1, r1) = (table.ell[i], table.ell[i + 1]);
        for j in 0..table.hyp.len() - 1 {
            let (l2, r2) = (table.hyp[j], table.hyp[j + 1]);
            if !table.occupied(i, j) {
                continue;
            }
            // Clip the cell to the accessible band.
            let c1l = l1;
            let c1r = if lc < b { r1.min(lc) } else { r1 };
            let c2l = if lc > b { l2.max(lc) } else { l2 };
            let c2r = r2;
            if c1r <= c1l || c2r <= c2l {
                continue;
            }
            boxes.push((u1.u(c1l), u1.u(c1r), u2.u(c2l), u2.u(c2r)));
        }
    }
    boxes
}

/// Sample the evolution of the phase manifold across the sweep (doc §8).
pub fn evolution(table: &Table, cf: &ConfocalParams, n: usize) -> Vec<FiberFrame> {
    let top = *table.hyp.last().unwrap();
    let crit: Vec<f32> = critical_values(table, cf, 1e-6);
    let mut frames = Vec::new();
    for k in 0..n {
        let lc = 1e-3 + (top - 2e-3) * k as f32 / (n as f32 - 1.0);
        let t = topology(&accessible_grid(table, cf, lc));
        if t.empty {
            continue;
        }
        let mut genus = t.genus.clone();
        genus.sort_unstable();
        let is_critical = crit.iter().any(|&c| (lc - c).abs() < top / n as f32);
        frames.push(FiberFrame {
            lc,
            n_components: t.n_components,
            genus,
            flat_rects: flat_rects(table, cf, lc),
            is_critical,
        });
    }
    frames
}
