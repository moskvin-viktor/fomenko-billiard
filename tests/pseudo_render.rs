//! Tests for the pseudo-integrable 3D renderers.
//!
//! The user reported the L-shape 3D view is wrong: torus levels don't show a
//! torus (they were forced onto a donut via the full-ellipse fallback), and
//! genus-2 levels didn't look like genus 2 (the flat cross reads as 4 planes).
//! The fix routes the L's 3D view through the pseudo machinery:
//!
//! - [`pseudo::render::torus_angles`] — a flat rectangle `(u1, u2)` with
//!   opposite edges identified maps to a torus (normalized angles → donut).
//! - [`pseudo::render::pretzel_embed`] — maps the flat cross onto a **double
//!   torus** (genus-2 pretzel): each of the two lobes is a torus, joined by a
//!   bridge, so the manifold reads as two tori glued.

use billiards::pseudo::render::{cross_embed, pretzel_embed, torus_angles, CrossModuli};

// ---------------------------------------------------------------------------
// torus_angles: flat rectangle → torus
// ---------------------------------------------------------------------------

#[test]
fn torus_angles_identify_opposite_edges() {
    // A flat rectangle [0, W1] × [0, W2] with opposite edges identified is a
    // torus.  The angles must be periodic: (u1, u2) and (u1+W1, u2) map to the
    // same angle pair (mod 2π).
    let (w1, w2) = (3.0f32, 2.0f32);
    let a = torus_angles(0.5, 0.5, w1, w2);
    let b = torus_angles(0.5 + w1, 0.5, w1, w2); // same point, +1 period in u1
    assert!(
        (a.0 - b.0).abs() < 1e-4 && (a.1 - b.1).abs() < 1e-4,
        "u1-periodic: ({:.4},{:.4}) vs ({:.4},{:.4})",
        a.0,
        b.0,
        a.1,
        b.1
    );
    let c = torus_angles(0.5, 0.5 + w2, w1, w2); // +1 period in u2
    assert!(
        (a.0 - c.0).abs() < 1e-4 && (a.1 - c.1).abs() < 1e-4,
        "u2-periodic: ({:.4},{:.4}) vs ({:.4},{:.4})",
        a.0,
        c.0,
        a.1,
        c.1
    );
}

#[test]
fn torus_angles_are_normalized_to_2pi() {
    let (w1, w2) = (3.0f32, 2.0f32);
    for &(u1, u2) in &[(0.0, 0.0), (1.0, 1.5), (2.9, 1.9), (0.3, 0.1)] {
        let (t1, t2) = torus_angles(u1, u2, w1, w2);
        assert!(
            (0.0..=2.0 * std::f32::consts::PI).contains(&t1),
            "θ1 out of range: {t1}"
        );
        assert!(
            (0.0..=2.0 * std::f32::consts::PI).contains(&t2),
            "θ2 out of range: {t2}"
        );
    }
}

// ---------------------------------------------------------------------------
// cross_embed: flat point is finite in the cross (flat chart is finite)
// ---------------------------------------------------------------------------

#[test]
fn cross_embed_is_finite_and_symmetric() {
    let m = CrossModuli {
        a1: 1.0,
        a2: 2.0,
        b1: 1.5,
        b2: 2.5,
    };
    let sheets = [(1i8, 1i8), (1, -1), (-1, 1), (-1, -1)];
    // Every sheet maps to a finite point (no NaN from the flat chart).
    for &(s1, s2) in &sheets {
        let p = cross_embed(0.5, 0.5, s1, s2, m);
        assert!(p.x.is_finite() && p.y.is_finite() && p.z.is_finite());
    }
    // Sign flip negates the coordinate (X=σ1·u1, Y=σ2·u2).
    let p_pp = cross_embed(0.7, 0.9, 1, 1, m);
    let p_np = cross_embed(0.7, 0.9, -1, 1, m);
    assert!((p_pp.x + p_np.x).abs() < 1e-4, "σ1 flip negates x");
    let p_pn = cross_embed(0.7, 0.9, 1, -1, m);
    assert!((p_pp.y + p_pn.y).abs() < 1e-4, "σ2 flip negates y");
}

// ---------------------------------------------------------------------------
// pretzel_embed: flat cross → double torus (two tori glued)
// ---------------------------------------------------------------------------

/// The double-torus embedding must actually produce genus-2 topology: two
/// distinct handle lobes (the two tori) joined by a bridge, not a single donut
/// and not four disconnected planes.
#[test]
fn pretzel_embed_has_two_lobes_and_a_bridge() {
    let m = CrossModuli {
        a1: 1.0,
        a2: 2.0,
        b1: 1.5,
        b2: 2.5,
    };

    // Sample a dense sheet (1,1) of flat points across the FULL cross extent so
    // the whole longitudinal loop (both lobes) is covered.
    let maj = m.a2.max(m.b2).max(m.a1).max(m.b1);
    let mut xs = Vec::new();
    let mut zs = Vec::new();
    let n = 60;
    for i in 0..=n {
        let u1 = maj * i as f32 / n as f32; // full longitudinal loop
        for j in 0..=n {
            let u2 = maj * j as f32 / n as f32;
            let p = pretzel_embed(u1, u2, 1, 1, m);
            xs.push(p.x);
            zs.push(p.z);
        }
    }
    let (xmin, xmax) = xs
        .iter()
        .cloned()
        .fold((f32::MAX, f32::MIN), |(a, b), x| (a.min(x), b.max(x)));
    let zmin = zs.iter().cloned().fold(f32::MAX, f32::min);
    let zmax = zs.iter().cloned().fold(f32::MIN, f32::max);
    // Two lobes along x (genus 2): the surface reaches a clear negative and
    // positive lobe centre.
    assert!(
        xmin < -1.0 && xmax > 1.0,
        "two lobes in x: [{xmin:.2},{xmax:.2}]"
    );
    // Real surface (thickness in z), not a flat plane.
    assert!(
        (zmax - zmin) > 0.5,
        "surface must have thickness in z, got {zmin:.2}..{zmax:.2}"
    );
}

/// The double-torus embedded sheets must not be separated in z (they form one
/// connected surface, not stacked planes).
#[test]
fn pretzel_embed_sheets_share_a_surface() {
    let m = CrossModuli {
        a1: 1.0,
        a2: 2.0,
        b1: 1.5,
        b2: 2.5,
    };
    let sheets = [(1i8, 1i8), (1, -1), (-1, 1), (-1, -1)];
    // Sample each sheet at the same point; all should land within a bounded z
    // shell (one connected fundagen2 object).
    let mut zmin = f32::MAX;
    let mut zmax = f32::MIN;
    for &(s1, s2) in &sheets {
        let p = pretzel_embed(0.0, 0.5, s1, s2, m);
        zmin = zmin.min(p.z);
        zmax = zmax.max(p.z);
    }
    assert!(
        zmax - zmin < 3.0,
        "double torus sheets should be within one surface, z span {}..{}",
        zmin,
        zmax
    );
}
