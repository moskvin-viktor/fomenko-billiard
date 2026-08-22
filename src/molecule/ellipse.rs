//! The full-ellipse front-end (real-billiard molecules).
//!
//! The plain ellipse crosses the axes and the focal segment, so its `λ`-chart
//! is a fold cover, not a rectangle.  Rather than force it into a `Table`, wrap
//! it as a front-end that emits the same molecule data, reusing the shared
//! vertex kinds.  This is the point of unification: the ellipse's atom-B
//! molecule and the L's genus-jump molecule are the *same data structure*,
//! differing only in which vertex kinds appear.

use super::critical::TransitionKind;
use super::reeb::{Molecule, MoleculeEdge, MoleculeVertex};

/// The classical confocal-ellipse molecule, reproducing the known A – B – A
/// tree (Fokicheva/Fomenko).
pub fn ellipse_molecule(a: f32, b: f32) -> Molecule {
    let vertices = vec![
        MoleculeVertex {
            lc: 0.0,
            kind: TransitionKind::AEnd,
            from: vec![],
            to: vec![],
        },
        MoleculeVertex {
            lc: b,
            kind: TransitionKind::SplitMerge,
            from: vec![1, 1],
            to: vec![1],
        },
        MoleculeVertex {
            lc: a,
            kind: TransitionKind::AEnd,
            from: vec![],
            to: vec![],
        },
    ];
    let edges = vec![
        MoleculeEdge {
            lo: 0.0,
            hi: b,
            n_components: 2,
            genus: vec![1, 1],
            sample_lc: 0.5 * b,
        },
        MoleculeEdge {
            lo: b,
            hi: a,
            n_components: 1,
            genus: vec![1],
            sample_lc: 0.5 * (a + b),
        },
    ];
    Molecule { vertices, edges }
}
