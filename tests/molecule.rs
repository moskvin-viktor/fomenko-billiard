//! Molecule-engine correctness tests.
//!
//! Drives the implementation of `src/molecule/` against the spec in
//! `docs/confocal_molecule_engine.md` §11 (worked molecules) and §12
//! (validation).  Each test targets one layer:
//!
//!   Layers 2–3  mask + topology             tests 1–7
//!   Layer  4    critical values/transitions tests 8–13
//!   Layer  5    molecule (Reeb graph)       tests 14–16
//!   Layer  6    flat geometry / evolution   tests 17–19
//!   Layer  7    critical fibers             tests 20–22
//!   ellipse     front-end                   test  23

use std::collections::HashSet;

use billiards::molecule::*;
use billiards::table::Table;
use billiards::torus::ConfocalParams;

/// The doc's worked family: a = 2.0, b = 1.0.
const A: f32 = 2.0;
const B: f32 = 1.0;

fn cf() -> ConfocalParams {
    ConfocalParams::new(A, B)
}

// ---------------------------------------------------------------------------
// Table library (doc §3)
// ---------------------------------------------------------------------------

/// confocal_square(a, b, al=0.05, ar=0.4, bl=1.2, br=1.8): one rectangle.
fn square_table() -> Table {
    confocal_square(A, B, 0.05, 0.4, 1.2, 1.8)
}

/// confocal_L(a, b, a1=0.25, a2=0.5, b1=1.2, b2=1.6, b3=1.9):
/// reflex corner at (a₁, b₂) = (0.25, 1.6).
fn l_table() -> Table {
    confocal_L(A, B, 0.25, 0.5, 1.2, 1.6, 1.9)
}

/// Two disjoint vertical strips with a λ₁-gap: above the gap wall the region
/// splits into two tori → exercises `split_merge` on a real Table (the doc's
/// in-quadrant L only ever produces `genus_jump`).
fn two_strip_table() -> Table {
    Table::from_bounds(&[0.05, 0.3, 0.5, 0.8], &[1.2, 1.9], |l1, _l2, _q| {
        l1 <= 0.3 || l1 >= 0.5
    })
}

// ---------------------------------------------------------------------------
// Layers 2 + 3 — mask and topology
// ---------------------------------------------------------------------------

/// λc < a₁: only the tall strip survives → single rectangle, genus 1.
#[test]
fn l_elliptic_below_corner_is_torus() {
    let g = accessible_grid(&l_table(), &cf(), 0.2);
    let t = topology(&g);
    assert_eq!(t.n_components, 1);
    assert_eq!(t.genus, vec![1]);
    assert_eq!(t.reflex, 0);
    assert!(!t.empty);
}

/// a₁ < λc < b: the L-shape survives with its reflex corner → genus 2.
#[test]
fn l_elliptic_with_corner_is_genus2() {
    let g = accessible_grid(&l_table(), &cf(), 0.7);
    let t = topology(&g);
    assert_eq!(t.n_components, 1);
    assert_eq!(t.genus, vec![2]);
    assert_eq!(t.reflex, 1);
}

/// b₂ < λc < b₃: only the tall strip survives → genus 1.
#[test]
fn l_hyperbolic_above_corner_is_torus() {
    let g = accessible_grid(&l_table(), &cf(), 1.7);
    let t = topology(&g);
    assert_eq!(t.n_components, 1);
    assert_eq!(t.genus, vec![1]);
    assert_eq!(t.reflex, 0);
}

/// λc ≥ lam2_top: nothing accessible.
#[test]
fn l_forbidden_beyond_top() {
    let g = accessible_grid(&l_table(), &cf(), 2.0);
    assert!(topology(&g).empty);
}

/// The square is a single genus-1 component at every accessible level.
#[test]
fn square_is_always_genus1() {
    for lc in [0.1f32, 0.3, 0.7, 1.1, 1.5] {
        let t = topology(&accessible_grid(&square_table(), &cf(), lc));
        assert_eq!(t.n_components, 1, "λc={lc}");
        assert_eq!(t.genus, vec![1], "λc={lc}");
        assert_eq!(t.reflex, 0, "λc={lc}");
    }
}

/// Doc §12.1 — genus two ways: g = 1 + n_reflex per component, where n_reflex
/// is counted independently via the 2×2 local vertex test.
#[test]
fn genus_equals_one_plus_reflex_per_component() {
    for lc in [0.2f32, 0.5, 0.7, 1.2, 1.7] {
        let g = accessible_grid(&l_table(), &cf(), lc);
        let t = topology(&g);
        let (label, ncomp) = components(&g);
        let (ni, nj) = (g.inside.len(), g.inside[0].len());
        let mut reflex_per_comp = vec![0usize; ncomp as usize];
        for vi in 0..=ni {
            for vj in 0..=nj {
                if classify_vertex(&g, vi, vj) == VertexType::Reflex {
                    let comp = component_of_vertex(&g, &label, vi, vj)
                        .expect("reflex vertex must be in a component");
                    reflex_per_comp[comp as usize] += 1;
                }
            }
        }
        assert_eq!(t.genus.len(), t.n_components);
        for (c, &g_comp) in t.genus.iter().enumerate() {
            assert_eq!(
                g_comp as usize,
                1 + reflex_per_comp[c],
                "genus of component {c} at λc={lc} must be 1 + n_reflex"
            );
        }
    }
}

/// Doc §12.2 — every component is a disk (V − E + F = 1).  A hole would break
/// the g = 1 + n_reflex hypothesis.
#[test]
fn every_component_is_a_disk() {
    for (tab, lcs) in [
        (&l_table(), &[0.2f32, 0.5, 0.7, 1.2, 1.7][..]),
        (&two_strip_table(), &[0.2f32, 0.7, 1.5][..]),
    ] {
        for &lc in lcs {
            let g = accessible_grid(tab, &cf(), lc);
            let (label, ncomp) = components(&g);
            let (ni, nj) = (g.inside.len(), g.inside[0].len());
            for c in 0..ncomp {
                let mut verts = HashSet::new();
                let mut edges = HashSet::new();
                let mut faces = 0usize;
                for i in 0..ni {
                    for j in 0..nj {
                        if label[i][j] != Some(c) {
                            continue;
                        }
                        faces += 1;
                        for &p in &[(i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1)] {
                            verts.insert(p);
                        }
                        for &(a, b) in &[
                            ((i, j), (i + 1, j)),
                            ((i, j + 1), (i + 1, j + 1)),
                            ((i, j), (i, j + 1)),
                            ((i + 1, j), (i + 1, j + 1)),
                        ] {
                            edges.insert((a, b));
                        }
                    }
                }
                let v = verts.len() as i64;
                let e = edges.len() as i64;
                let f = faces as i64;
                assert_eq!(
                    v - e + f,
                    1,
                    "component {c} at λc={lc} must be a disk, got V−E+F={}",
                    v - e + f
                );
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Layer 4 — critical values and transitions
// ---------------------------------------------------------------------------

#[test]
fn critical_values_of_l() {
    let crit = critical_values(&l_table(), &cf(), 1e-6);
    let expected = [0.25, 0.5, 1.0, 1.2, 1.6];
    assert_eq!(crit.len(), expected.len());
    for (c, e) in crit.iter().zip(expected.iter()) {
        assert!((c - e).abs() < 1e-4, "critical value {c} != {e}");
    }
}

#[test]
fn classify_l_transitions() {
    // Reflex corner crossed by the elliptic caustic: genus 1 → 2.
    let t = classify_transition(&l_table(), &cf(), 0.25, 1e-4);
    assert_eq!(t.kind, TransitionKind::GenusJump);
    assert_eq!(t.from, vec![1]);
    assert_eq!(t.to, vec![2]);

    // Reflex corner crossed by the hyperbolic caustic: genus 2 → 1.
    let t = classify_transition(&l_table(), &cf(), 1.6, 1e-4);
    assert_eq!(t.kind, TransitionKind::GenusJump);
    assert_eq!(t.from, vec![2]);
    assert_eq!(t.to, vec![1]);

    // Cosmetic wall crossings: no topology change.
    for lc in [0.5f32, 1.0, 1.2] {
        let t = classify_transition(&l_table(), &cf(), lc, 1e-4);
        assert_eq!(
            t.kind,
            TransitionKind::EdgeSwap,
            "λc={lc} should be edge_swap"
        );
    }
}

/// Doc §12.6 — an in-quadrant table not reaching the focal segment gives a
/// non-event (edge_swap) at λc = b, automatically.
#[test]
fn focal_value_is_edge_swap_for_in_quadrant_table() {
    let t = classify_transition(&l_table(), &cf(), B, 1e-4);
    assert_eq!(t.kind, TransitionKind::EdgeSwap);
}

/// Crossing the λ₁-gap wall merges/splits the two strips → split_merge (atom B).
#[test]
fn two_strip_table_splits_at_gap() {
    let t = classify_transition(&two_strip_table(), &cf(), 0.5, 1e-4);
    assert_eq!(t.kind, TransitionKind::SplitMerge);
    assert_eq!(t.from, vec![1]);
    assert_eq!(t.to, vec![1, 1]);
}

/// The square appears at its inner ellipse wall → birth_A.
#[test]
fn square_birth_at_inner_wall() {
    let t = classify_transition(&square_table(), &cf(), 0.05, 1e-4);
    assert_eq!(t.kind, TransitionKind::BirthA);
}

/// Doc §12.3 — every topology change point must coincide with a critical value.
#[test]
fn critical_values_are_exhaustive() {
    let crit = critical_values(&l_table(), &cf(), 1e-6);
    let mut prev: Option<Signature> = None;
    let mut prev_lc = 0.0f32;
    let n = 400;
    for k in 0..=n {
        let lc = 0.01 + 1.88 * k as f32 / n as f32;
        let sig = topo_signature(&l_table(), &cf(), lc);
        if let Some(p) = &prev {
            if sig != *p {
                let mid = 0.5 * (prev_lc + lc);
                assert!(
                    crit.iter().any(|&c| (c - mid).abs() < 0.02),
                    "topology changed near λc={mid:.3}, not in critical_values"
                );
            }
        }
        prev = Some(sig);
        prev_lc = lc;
    }
}

// ---------------------------------------------------------------------------
// Layer 5 — the molecule (Reeb graph)
// ---------------------------------------------------------------------------

/// Doc §11 — the L's molecule is A – genus_jump – genus_jump – A, with the
/// genus sequence 1 → 2 → 1 across the edges.
#[test]
fn l_molecule_structure() {
    let m = build_molecule(&l_table(), &cf());
    let kinds: Vec<TransitionKind> = m.vertices.iter().map(|v| v.kind).collect();
    assert_eq!(
        kinds,
        vec![
            TransitionKind::AEnd,
            TransitionKind::GenusJump,
            TransitionKind::GenusJump,
            TransitionKind::AEnd
        ]
    );
    assert!((m.vertices[1].lc - 0.25).abs() < 1e-4);
    assert!((m.vertices[2].lc - 1.6).abs() < 1e-4);

    // 6 edges (one per interval between consecutive critical values), genus
    // sequence 1,2,2,2,2,1.
    assert_eq!(m.edges.len(), 6);
    let genus_seq: Vec<u32> = m.edges.iter().map(|e| e.genus[0]).collect();
    assert_eq!(genus_seq, vec![1, 2, 2, 2, 2, 1]);

    // The summary reproduces the doc's A – jump – jump – A molecule.
    assert_eq!(
        m.summary(),
        "[A_end@0.000] --(comp=1,g=(1,))-- [genus_jump@0.250] \
         --(comp=1,g=(2,))-- --(comp=1,g=(2,))-- --(comp=1,g=(2,))-- \
         --(comp=1,g=(2,))-- [genus_jump@1.600] --(comp=1,g=(1,))-- [A_end@1.900]"
    );
}

/// Doc §11 — the square is fully integrable: only A_end caps, all edges genus 1.
#[test]
fn square_molecule_is_all_genus1() {
    let m = build_molecule(&square_table(), &cf());
    let kinds: Vec<TransitionKind> = m.vertices.iter().map(|v| v.kind).collect();
    assert_eq!(kinds, vec![TransitionKind::AEnd, TransitionKind::AEnd]);
    assert!(!m.edges.is_empty());
    for e in &m.edges {
        assert_eq!(e.genus, vec![1]);
        assert_eq!(e.n_components, 1);
    }
}

/// Doc §12.7 — every maximal genus-1 branch ends in an A_end.
#[test]
fn genus1_branches_end_in_a_end() {
    let m = build_molecule(&l_table(), &cf());
    assert_eq!(m.edges.first().unwrap().genus, vec![1]);
    assert_eq!(m.edges.last().unwrap().genus, vec![1]);
    assert_eq!(m.vertices.first().unwrap().kind, TransitionKind::AEnd);
    assert_eq!(m.vertices.last().unwrap().kind, TransitionKind::AEnd);
}

// ---------------------------------------------------------------------------
// Layer 6 — flat geometry and evolution
// ---------------------------------------------------------------------------

#[test]
fn flat_rects_are_finite_and_positive() {
    for lc in [0.2f32, 0.7, 1.2, 1.7] {
        let boxes = flat_rects(&l_table(), &cf(), lc);
        assert!(!boxes.is_empty(), "λc={lc} should be accessible");
        for &(u1l, u1r, u2l, u2r) in &boxes {
            assert!(u1l.is_finite() && u1r.is_finite() && u2l.is_finite() && u2r.is_finite());
            assert!(u1r > u1l, "positive u₁ width at λc={lc}");
            assert!(u2r > u2l, "positive u₂ width at λc={lc}");
        }
    }
    // Forbidden level → no boxes.
    assert!(flat_rects(&l_table(), &cf(), 2.0).is_empty());
}

/// One flat box per accessible cell — the flat polyomino tiles the region.
#[test]
fn flat_rects_match_accessible_cells() {
    for lc in [0.2f32, 0.7, 1.2, 1.7] {
        let g = accessible_grid(&l_table(), &cf(), lc);
        let n_accessible = g.inside.iter().flatten().filter(|&&b| b).count();
        let boxes = flat_rects(&l_table(), &cf(), lc);
        assert_eq!(
            boxes.len(),
            n_accessible,
            "one flat box per accessible cell at λc={lc}"
        );
    }
}

/// The evolution sweep passes through genus 1 → 2 → 1 with finite flat rects.
#[test]
fn evolution_tracks_genus() {
    let frames = evolution(&l_table(), &cf(), 200);
    assert!(!frames.is_empty());
    assert_eq!(frames.first().unwrap().genus[0], 1);
    assert_eq!(frames.last().unwrap().genus[0], 1);
    assert!(
        frames.iter().any(|f| f.genus[0] == 2),
        "must pass through genus 2"
    );
    for f in &frames {
        for &(a, b, c, d) in &f.flat_rects {
            assert!(a.is_finite() && b.is_finite() && c.is_finite() && d.is_finite());
        }
    }
}

// ---------------------------------------------------------------------------
// Layer 7 — critical fibers
// ---------------------------------------------------------------------------

#[test]
fn critical_fiber_genus_jump() {
    for lc in [0.25f32, 1.6] {
        let f = critical_fiber(&l_table(), &cf(), lc).expect("genus_jump fiber");
        assert_eq!(f.kind, "genus_jump");
        assert_eq!(f.normalization_genus, 1);
        assert_eq!(f.n_nodes, 1);
        assert_eq!(f.arithmetic_genus, 2);
        assert!(f.spine.contains("figure-eight"));
    }
}

#[test]
fn critical_fiber_none_for_edge_swap() {
    assert!(critical_fiber(&l_table(), &cf(), 0.5).is_none());
    assert!(critical_fiber(&l_table(), &cf(), 1.0).is_none());
}

/// Doc §12.5 — normalization_genus == max genus below, arithmetic_genus == max
/// genus above.
#[test]
fn arithmetic_genus_matches_topology() {
    let lc = 0.25;
    let f = critical_fiber(&l_table(), &cf(), lc).unwrap();
    let below = topo_signature(&l_table(), &cf(), lc - 1e-4);
    let above = topo_signature(&l_table(), &cf(), lc + 1e-4);
    assert_eq!(
        f.normalization_genus,
        below.genus.iter().max().copied().unwrap_or(0)
    );
    assert_eq!(
        f.arithmetic_genus,
        above.genus.iter().max().copied().unwrap_or(0)
    );
}

// ---------------------------------------------------------------------------
// Ellipse front-end (doc §10, §12.8)
// ---------------------------------------------------------------------------

/// The classical A – B – A tree, reproduced through the shared Molecule type.
#[test]
fn ellipse_molecule_is_aba() {
    let m = ellipse_molecule(A, B);
    let kinds: Vec<TransitionKind> = m.vertices.iter().map(|v| v.kind).collect();
    assert_eq!(
        kinds,
        vec![
            TransitionKind::AEnd,
            TransitionKind::SplitMerge,
            TransitionKind::AEnd
        ]
    );
    assert_eq!(m.edges.len(), 2);
    assert_eq!(m.edges[0].genus, vec![1, 1]); // two tori below the separatrix
    assert_eq!(m.edges[1].genus, vec![1]); // one torus above
    assert_eq!(m.edges[0].n_components, 2);
    assert_eq!(m.edges[1].n_components, 1);
    assert_eq!(
        m.summary(),
        "[A_end@0.000] --(comp=2,g=(1, 1))-- [split_merge@1.000] \
         --(comp=1,g=(1,))-- [A_end@2.000]"
    );
}
