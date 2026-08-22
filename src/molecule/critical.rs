//! Layer 4 of the molecule engine: the finite set of critical caustic values
//! and the classification of the transition at each.
//!
//! Combinatorial changes of the mask happen only when the cut line `λc` crosses
//! a table edge coordinate, or at the focal value `b`.  The candidate set is
//! therefore finite and explicit, and each candidate is classified by comparing
//! the topology just below and just above it.

use crate::table::Table;
use crate::torus::ConfocalParams;

use super::grid::{accessible_grid, topology};

/// A coarse topology signature, comparable across levels (doc §6).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Signature {
    /// Number of connected components.
    pub n_components: usize,
    /// Sorted genus multiset.
    pub genus: Vec<u32>,
    /// Number of pinch vertices.
    pub pinch: usize,
    /// Whether nothing is accessible.
    pub empty: bool,
}

/// The kind of a molecule vertex / transition (doc §6).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TransitionKind {
    /// The region appears (elliptic centre end) — a Fomenko atom A.
    BirthA,
    /// The region empties — a Fomenko atom A.
    DeathA,
    /// A reflex corner crossed the caustic cut: genus changed, a stratum
    /// degeneration (not an atom).
    GenusJump,
    /// Components merged/split through a pinch — an atom-B saddle.
    SplitMerge,
    /// A boundary edge swapped wall↔turning-line with no topology change.
    EdgeSwap,
    /// A molecule cap (outermost accessible edge).
    AEnd,
}

impl TransitionKind {
    /// The doc's string label.
    pub fn label(self) -> &'static str {
        match self {
            TransitionKind::BirthA => "birth_A",
            TransitionKind::DeathA => "death_A",
            TransitionKind::GenusJump => "genus_jump",
            TransitionKind::SplitMerge => "split_merge",
            TransitionKind::EdgeSwap => "edge_swap",
            TransitionKind::AEnd => "A_end",
        }
    }
}

/// A classified transition at a critical value (doc §6).
#[derive(Clone, Debug)]
pub struct Transition {
    pub kind: TransitionKind,
    /// Topology just below `lc`.
    pub below: Signature,
    /// Topology just above `lc`.
    pub above: Signature,
    /// Genus multiset below.
    pub from: Vec<u32>,
    /// Genus multiset above.
    pub to: Vec<u32>,
}

/// The finite set of critical caustic values: table edge coordinates plus the
/// focal value `b`, restricted to the swept range `(0, lam2_top - eps)`.
pub fn critical_values(table: &Table, cf: &ConfocalParams, eps: f32) -> Vec<f32> {
    let top = *table.hyp.last().unwrap();
    let mut cands: Vec<f32> = table
        .ell
        .iter()
        .chain(table.hyp.iter())
        .copied()
        .chain(std::iter::once(cf.b))
        .collect();
    cands.sort_by(|x, y| x.total_cmp(y));
    cands.dedup();
    cands
        .into_iter()
        .filter(|&c| c > 0.0 && c < top - eps)
        .collect()
}

/// The topology signature at a level.
pub fn topo_signature(table: &Table, cf: &ConfocalParams, lc: f32) -> Signature {
    let t = topology(&accessible_grid(table, cf, lc));
    let mut genus = t.genus.clone();
    genus.sort_unstable();
    Signature {
        n_components: t.n_components,
        genus,
        pinch: t.pinch,
        empty: t.empty,
    }
}

/// Classify the transition at a critical value by comparing the topology just
/// below and just above (doc §6).
///
/// Priority: birth/death, then component-count change or pinch → `split_merge`,
/// then genus-multiset change → `genus_jump`, else `edge_swap`.
pub fn classify_transition(table: &Table, cf: &ConfocalParams, lc: f32, eps: f32) -> Transition {
    let below = topo_signature(table, cf, lc - eps);
    let above = topo_signature(table, cf, lc + eps);

    let kind = if below.empty && !above.empty {
        TransitionKind::BirthA
    } else if above.empty && !below.empty {
        TransitionKind::DeathA
    } else if below.n_components != above.n_components || above.pinch > 0 || below.pinch > 0 {
        TransitionKind::SplitMerge
    } else if below.genus != above.genus {
        TransitionKind::GenusJump
    } else {
        TransitionKind::EdgeSwap
    };

    Transition {
        kind,
        from: below.genus.clone(),
        to: above.genus.clone(),
        below,
        above,
    }
}
