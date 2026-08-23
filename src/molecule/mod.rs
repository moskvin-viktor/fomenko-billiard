//! A unified molecule engine for confocal billiards.
//!
//! One implementation that, given any confocal table, produces the full
//! evolution of the phase manifolds as the caustic parameter `λc` sweeps —
//! serving both the Liouville-integrable tables (confocal rectangles, the
//! ellipse: genus-1 fibers, Fomenko atoms) and the pseudo-integrable reflex-
//! corner tables (L, T, Z, plus: genus-jumping fibers, stratum degenerations).
//! The output is the molecule (Reeb graph over the `λc`-axis) with every vertex
//! classified, plus per-level flat geometry you can animate.
//!
//! See `docs/confocal_molecule_engine.md` for the full design.
//!
//! - [`grid`] — Layers 2–3: the accessible-region mask and its topology.
//! - [`critical`] — Layer 4: critical values and transition classification.
//! - [`reeb`] — Layer 5: the Reeb-graph molecule.
//! - [`evolution`] — Layer 6: flat geometry and phase-manifold evolution.
//! - [`fiber`] — Layer 7: critical fibers.
//! - [`ellipse`] — the full-ellipse front-end.

pub mod critical;
pub mod ellipse;
pub mod evolution;
pub mod fiber;
pub mod grid;
pub mod reeb;
pub mod special;

pub use critical::{
    classify_transition, critical_values, topo_signature, Signature, Transition, TransitionKind,
};
pub use ellipse::ellipse_molecule;
pub use evolution::{evolution, flat_rects, FiberFrame};
pub use fiber::{critical_fiber, CriticalFiber};
pub use grid::{
    accessible_grid, classify_vertex, component_of_vertex, components, topology, Grid, Topology,
    VertexType,
};
pub use reeb::{build_molecule, Molecule, MoleculeEdge, MoleculeVertex};
pub use special::{all_layers, layer_label, special_layers};

use crate::table::Table;

/// The doc's `confocal_square`: one rectangle `[al, ar] × [bl, br]` → genus 1
/// (integrable).
pub fn confocal_square(_a: f32, _b: f32, al: f32, ar: f32, bl: f32, br: f32) -> Table {
    Table::from_bounds(&[al, ar], &[bl, br], |l1, l2, _q| {
        l1 >= al && l1 <= ar && l2 >= bl && l2 <= br
    })
}

/// The doc's `confocal_L`: a tall leg plus a base, with a reflex corner at
/// `(a₁, b₂)`.
#[allow(non_snake_case)]
pub fn confocal_L(_a: f32, _b: f32, a1: f32, a2: f32, b1: f32, b2: f32, b3: f32) -> Table {
    Table::from_bounds(&[0.0, a1, a2], &[b1, b2, b3], |l1, l2, _q| {
        let tall = l1 <= a1 && l2 >= b2;
        let base = l1 <= a2 && l2 <= b2;
        tall || base
    })
}
