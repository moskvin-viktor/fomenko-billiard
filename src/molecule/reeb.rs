//! Layer 5 of the molecule engine: assemble critical values (vertices) and the
//! open intervals between them (edges) into the Reeb graph over the `λc` axis.
//!
//! `build_molecule` is the single entry point for both worlds: a rectangle in →
//! an all-atom molecule (integrable); an L in → a molecule with two
//! `genus_jump` vertices (pseudo-integrable).  The builder does not know which
//! is which.

use crate::table::Table;
use crate::torus::ConfocalParams;

use super::critical::{classify_transition, critical_values, TransitionKind};
use super::grid::{accessible_grid, topology};

/// An open interval `(lo, hi)` of caustic values with constant topology.
#[derive(Clone, Debug)]
pub struct MoleculeEdge {
    pub lo: f32,
    pub hi: f32,
    pub n_components: usize,
    /// Sorted genus per component.
    pub genus: Vec<u32>,
    /// A representative caustic value in the interval.
    pub sample_lc: f32,
}

/// A critical value where the topology changes.
#[derive(Clone, Debug)]
pub struct MoleculeVertex {
    pub lc: f32,
    pub kind: TransitionKind,
    /// Genus multiset below.
    pub from: Vec<u32>,
    /// Genus multiset above.
    pub to: Vec<u32>,
}

/// The molecule: a Reeb graph over the `λc` axis.
#[derive(Clone, Debug)]
pub struct Molecule {
    pub vertices: Vec<MoleculeVertex>,
    pub edges: Vec<MoleculeEdge>,
}

impl Molecule {
    /// The doc's `summary()`: a compact string of vertices and edges.
    pub fn summary(&self) -> String {
        let mut parts = Vec::new();
        let mut vs = self.vertices.clone();
        vs.sort_by(|a, b| a.lc.total_cmp(&b.lc));
        let mut es = self.edges.clone();
        es.sort_by(|a, b| a.lo.total_cmp(&b.lo));
        let mut i = 0usize;
        for e in &es {
            while i < vs.len() && vs[i].lc <= e.lo + 1e-9 {
                parts.push(format!("[{}@{:.3}]", vs[i].kind.label(), vs[i].lc));
                i += 1;
            }
            let g = if e.genus.len() == 1 {
                format!("({},)", e.genus[0])
            } else {
                format!(
                    "({})",
                    e.genus
                        .iter()
                        .map(|g| g.to_string())
                        .collect::<Vec<_>>()
                        .join(", ")
                )
            };
            parts.push(format!("--(comp={},g={})--", e.n_components, g));
        }
        while i < vs.len() {
            parts.push(format!("[{}@{:.3}]", vs[i].kind.label(), vs[i].lc));
            i += 1;
        }
        parts.join(" ")
    }
}

/// Build the molecule for a table (doc §7).
pub fn build_molecule(table: &Table, cf: &ConfocalParams) -> Molecule {
    let crit = critical_values(table, cf, 1e-6);
    let top = *table.hyp.last().unwrap();
    let bounds: Vec<f32> = std::iter::once(0.0)
        .chain(crit.iter().copied())
        .chain(std::iter::once(top))
        .collect();

    let mut edges = Vec::new();
    for pair in bounds.windows(2) {
        let (lo, hi) = (pair[0], pair[1]);
        let mid = 0.5 * (lo + hi);
        let t = topology(&accessible_grid(table, cf, mid));
        if t.empty {
            continue;
        }
        let mut genus = t.genus.clone();
        genus.sort_unstable();
        edges.push(MoleculeEdge {
            lo,
            hi,
            n_components: t.n_components,
            genus,
            sample_lc: mid,
        });
    }

    let mut vertices = Vec::new();
    for &lc in &crit {
        let tr = classify_transition(table, cf, lc, 1e-4);
        if tr.kind == TransitionKind::EdgeSwap {
            continue;
        }
        // A birth/death at the very boundary coincides with the cap; keep the
        // cap and drop the coincident atom so the molecule reads cleanly.
        let at_boundary = (lc - edges.first().map_or(0.0, |e| e.lo)).abs() < 1e-6
            || (lc - edges.last().map_or(0.0, |e| e.hi)).abs() < 1e-6;
        if at_boundary && (tr.kind == TransitionKind::BirthA || tr.kind == TransitionKind::DeathA) {
            continue;
        }
        vertices.push(MoleculeVertex {
            lc,
            kind: tr.kind,
            from: tr.from,
            to: tr.to,
        });
    }

    // Endpoints: the outermost accessible edges cap with atom A.
    if let Some(first) = edges.first() {
        vertices.push(MoleculeVertex {
            lc: first.lo,
            kind: TransitionKind::AEnd,
            from: vec![],
            to: vec![],
        });
    }
    if let Some(last) = edges.last() {
        vertices.push(MoleculeVertex {
            lc: last.hi,
            kind: TransitionKind::AEnd,
            from: vec![],
            to: vec![],
        });
    }

    vertices.sort_by(|a, b| a.lc.total_cmp(&b.lc));
    Molecule { vertices, edges }
}
