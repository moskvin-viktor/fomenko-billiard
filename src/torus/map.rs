//! The top-level phase-space → torus mapping.
//!
//! Maps one phase-space sample `(x, y, vx, vy)` to a point `(θ₁, θ₂,
//! torus_index)` on a Liouville torus.  Everything reduces to the single cubic
//!
//! ```text
//! w² = P(λ) = (a − λ)(b − λ)(λc − λ)
//! ```
//!
//! and the phases are built from its Abelian differentials `dλ/√P` via the
//! `Libration` quadrature tables.  This is a faithful port of the reference
//! implementation in `docs/thorus_params.md`.

use crate::torus::confocal::{caustic, confocal, lam_dots, nu_angle, ConfocalParams, PhaseSample};
use crate::torus::quadrature::Libration;

/// Per-trajectory cache of the precomputed `Libration`s.  `λc` is fixed along
/// a trajectory, so build the splines once and reuse them for every sample.
#[derive(Default)]
pub struct TorusCache {
    l1: Option<Libration>,
    l2: Option<Libration>,
}

impl TorusCache {
    /// Number of interpolation knots per `Libration`.  Higher = smoother phase
    /// but costlier to build (built once per trajectory, then O(log n) lookups).
    const KNOTS: usize = 512;

    /// The `λ₁` / `λ₂` `Libration` for this trajectory, building and caching it
    /// on first use.  Each oval is fixed along a trajectory (`λc` constant), so
    /// the spline is built at most once.
    fn l1(&mut self, lo: f32, hi: f32, roots: [f32; 3]) -> &Libration {
        self.l1
            .get_or_insert_with(|| Libration::new(lo, hi, roots, Self::KNOTS))
    }

    /// See [`Self::l1`].
    fn l2(&mut self, lo: f32, hi: f32, roots: [f32; 3]) -> &Libration {
        self.l2
            .get_or_insert_with(|| Libration::new(lo, hi, roots, Self::KNOTS))
    }
}

/// All the per-trajectory context `to_torus` needs, bundled so the hot mapping
/// takes two arguments instead of ten and the domain-shape parameters are never
/// re-derived or passed individually.
pub struct TorusParams<'a> {
    /// The confocal family `(a, b)`.
    pub confocal: &'a ConfocalParams,
    /// `λ₁` of the outer boundary (0 for the base ellipse).
    pub lam_wall: f32,
    /// `λ₂` of the hyperbola walls; `None` means the full ellipse.
    pub beta: Option<f32>,
    /// Separatrix tolerance for `|λc − b|`.
    pub sep_eps: f32,
    /// Reused quadrature cache (one per trajectory).
    pub cache: &'a mut TorusCache,
}

/// Values derived once per sample in `to_torus`, shared by the three case
/// functions so nothing is recomputed and the signatures stay small.
struct CaseParams {
    /// Confocal coordinates `(λ₁, λ₂)`.
    lam1: f32,
    lam2: f32,
    /// Caustic parameter `λc`.
    lc: f32,
    /// Time derivatives of `(λ₁, λ₂)`.
    d1: f32,
    d2: f32,
    /// The cubic roots `[a, b, λc]`.
    roots: [f32; 3],
}

/// Map one phase-space sample to `(θ₁, θ₂, torus_index)`.
///
/// `params` carries the confocal family, the domain-shape parameters and the
/// per-trajectory cache, so the signature stays small and the expensive
/// `Libration` splines are built once per trajectory rather than per sample.
pub fn to_torus(s: &PhaseSample, params: &mut TorusParams) -> (f32, f32, u32) {
    let cf = params.confocal;
    let (lam1, lam2) = confocal(s.x, s.y, cf);
    let lc = caustic(s, cf);
    let (d1, d2) = lam_dots(s, lam1, lam2, cf);
    let root = CaseParams {
        lam1,
        lam2,
        lc,
        d1,
        d2,
        roots: [cf.a, cf.b, lc],
    };

    // Separatrix: the torus degenerates here.  Return a sentinel index.
    if (lc - cf.b).abs() < params.sep_eps {
        return (0.0, 0.0, u32::MAX);
    }

    match params.beta {
        None if lc < cf.b => elliptic_full_ellipse(params, s, &root),
        Some(_) if lc < cf.b => elliptic_confocal_square(params, s, &root),
        _ => hyperbolic(params, s, &root),
    }
}

/// Case A — elliptic caustic on a full ellipse (no hyperbola walls): `β = None`.
/// `ν` circulates and the torus splits by the sign of the angular momentum.
fn elliptic_full_ellipse(
    params: &mut TorusParams,
    s: &PhaseSample,
    r: &CaseParams,
) -> (f32, f32, u32) {
    let cf = params.confocal;
    let l1 = params.cache.l1(params.lam_wall, r.lc, r.roots);
    let th1 = l1.theta(r.lam1, r.d1 > 0.0);

    let th2 = nu_angle(s.x, s.y, r.lam1, cf).rem_euclid(2.0 * std::f32::consts::PI);
    let idx = if s.x * s.vy - s.y * s.vx > 0.0 { 0 } else { 1 };
    (th1, th2, idx)
}

/// Case C — elliptic caustic on a confocal square (`β` wall): folded libration.
/// The accessible region splits into an upper and lower part, split by `sign(y)`.
fn elliptic_confocal_square(
    params: &mut TorusParams,
    s: &PhaseSample,
    r: &CaseParams,
) -> (f32, f32, u32) {
    let cf = params.confocal;
    // Compute `th1` fully before further borrowing `params.cache` for `l2`.
    let th1 = params
        .cache
        .l1(params.lam_wall, r.lc, r.roots)
        .theta(r.lam1, r.d1 > 0.0);

    let l2 = params.cache.l2(params.beta.unwrap(), cf.a, r.roots);
    let tau = if s.x >= 0.0 { 1.0 } else { -1.0 };
    let u2 = l2.w_full - tau * l2.tail(r.lam2);
    let frac = std::f32::consts::PI * u2 / (2.0 * l2.w_full);
    let th2 = if tau * r.d2 > 0.0 {
        frac
    } else {
        2.0 * std::f32::consts::PI - frac
    };
    let idx = if s.y >= 0.0 { 0 } else { 1 };
    (th1, th2, idx)
}

/// Cases B / C-hyperbola — hyperbolic caustic: both coordinates fold, and the
/// accessible region is a single torus.
///
/// Angle convention: θ₁ is always the angle whose libration *collapses* at the
/// degenerate levels of the band (matching the elliptic cases, where θ₁ is the
/// wall–caustic libration).  Here that is the `λ₂` libration between the
/// caustic `λc` and `a` — its width vanishes as `λc → a` (focal axis).  The
/// surviving `λ₁` cycle (wall → b → wall) is θ₂.  With this convention the
/// standard-donut embedder can always render a collapse as the *tube* radius
/// shrinking (torus thins onto its equator ring) instead of the hole closing.
fn hyperbolic(params: &mut TorusParams, s: &PhaseSample, r: &CaseParams) -> (f32, f32, u32) {
    let cf = params.confocal;

    // The collapsing `λ₂` libration (caustic ↔ a): θ₁.
    let th1 = {
        let h2 = params.cache.l2(r.lc, cf.a, r.roots);
        let tau = if s.x >= 0.0 { 1.0 } else { -1.0 };
        let u2 = h2.w_full - tau * h2.tail(r.lam2);
        let f2 = std::f32::consts::PI * u2 / (2.0 * h2.w_full);
        if tau * r.d2 > 0.0 {
            f2
        } else {
            2.0 * std::f32::consts::PI - f2
        }
    };

    // The surviving `λ₁` cycle (wall → b → wall): θ₂.
    let sig = if s.y >= 0.0 { 1.0 } else { -1.0 };
    let th2 = {
        let h1 = params.cache.l1(params.lam_wall, cf.b, r.roots);
        let u1 = h1.w_full + sig * h1.tail(r.lam1);
        let f1 = std::f32::consts::PI * u1 / (2.0 * h1.w_full);
        if -sig * r.d1 > 0.0 {
            f1
        } else {
            2.0 * std::f32::consts::PI - f1
        }
    };

    (th1, th2, 0)
}
