//! Liouville tori for confocal billiards, in the Jacobi normalization.
//!
//! Maps a phase-space sample `(x, y, vx, vy)` to a point `(θ₁, θ₂, torus_index)`
//! on a Liouville torus.  Everything reduces to the single cubic
//!
//! ```text
//! w² = P(λ) = (a − λ)(b − λ)(λc − λ)
//! ```
//!
//! and the phases are built from its Abelian differentials `dλ/√P`.  This is a
//! faithful port of the reference implementation in `docs/thorus_params.md`.

use macroquad::prelude::*;

// ---------------------------------------------------------------------------
// Coordinates
// ---------------------------------------------------------------------------

/// The confocal family `(a, b)` the whole module works in.  Every quadric,
/// phase-space function and torus map lives in one such family, so threading a
/// single value instead of two scalars cuts the noise on every signature.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConfocalParams {
    pub a: f32,
    pub b: f32,
}

impl ConfocalParams {
    pub fn new(a: f32, b: f32) -> Self {
        Self { a, b }
    }

    /// The confocal family in which every billiard in this crate lives
    /// (`crate::A`, `crate::B`).
    pub fn standard() -> Self {
        Self::new(crate::A, crate::B)
    }
}

/// A point in phase space: position `(x, y)` and velocity `(vx, vy)`.
/// Used so sampled points are not passed as a jumble of positional scalars.
#[derive(Clone, Copy, Debug)]
pub struct PhaseSample {
    pub x: f32,
    pub y: f32,
    pub vx: f32,
    pub vy: f32,
}

impl PhaseSample {
    pub fn new(x: f32, y: f32, vx: f32, vy: f32) -> Self {
        Self { x, y, vx, vy }
    }
}

/// Confocal coordinates `(λ₁, λ₂)` with `λ₁ ≤ b ≤ λ₂ ≤ a`, as roots of
/// `λ² − (a + b − x² − y²)λ + (ab − bx² − ay²) = 0`.
///
/// Sign-safe near the foci, where the discriminant `(λ₂ − λ₁)²` vanishes:
/// compute the larger-magnitude root by the standard formula and get the
/// other from the product `q / λ_big`.
pub fn confocal(x: f32, y: f32, cf: &ConfocalParams) -> (f32, f32) {
    let a = cf.a;
    let b = cf.b;
    let p = a + b - x * x - y * y; // λ₁ + λ₂
    let q = a * b - b * x * x - a * y * y; // λ₁ · λ₂
    let disc = (p * p - 4.0 * q).max(0.0); // (λ₂ − λ₁)²
    let big = if p >= 0.0 {
        0.5 * (p + disc.sqrt())
    } else {
        0.5 * (p - disc.sqrt())
    };
    let small = if big != 0.0 { q / big } else { 0.0 };
    if small <= big {
        (small, big)
    } else {
        (big, small)
    }
}

/// Caustic parameter `λc`.  Normalizes `v`, so any speed is fine.
pub fn caustic(s: &PhaseSample, cf: &ConfocalParams) -> f32 {
    let n = (s.vx * s.vx + s.vy * s.vy).sqrt();
    let (vx, vy) = (s.vx / n, s.vy / n);
    let l = s.x * vy - s.y * vx;
    cf.b - l * l + (cf.a - cf.b) * vy * vy
}

/// Unfolded angle around the ellipse `λ = λ₁`.  Exact:
/// `x = √(a − λ₁) cos ν`, `y = √(b − λ₁) sin ν`.
pub fn nu_angle(x: f32, y: f32, lam1: f32, cf: &ConfocalParams) -> f32 {
    let cx = (cf.a - lam1).max(1e-30).sqrt();
    let cy = (cf.b - lam1).max(1e-30).sqrt();
    (y / cy).atan2(x / cx)
}

/// Time derivatives of the confocal coordinates, by implicit differentiation
/// of the quadratic — no square roots, no branch ambiguity.
pub fn lam_dots(s: &PhaseSample, lam1: f32, lam2: f32, cf: &ConfocalParams) -> (f32, f32) {
    let a = cf.a;
    let b = cf.b;
    let d = lam1 - lam2;
    let d1 = 2.0 * (s.x * (b - lam1) * s.vx + s.y * (a - lam1) * s.vy) / d;
    let d2 = -2.0 * (s.x * (b - lam2) * s.vx + s.y * (a - lam2) * s.vy) / d;
    (d1, d2)
}

// ---------------------------------------------------------------------------
// The cubic and its integrals
// ---------------------------------------------------------------------------

/// `P(λ) / (roots[skip] − λ)`: the product of the other two factors.
fn p_without(lam: f32, roots: [f32; 3], skip: usize) -> f32 {
    let mut out = 1.0;
    for (j, r) in roots.iter().enumerate() {
        if j != skip {
            out *= r - lam;
        }
    }
    out
}

/// Monotone phase `w(λ)/W` on an oval of the cubic `P`.
///
/// Two cases, both removing the inverse-square-root endpoint singularities
/// exactly so the grid is uniform where the integrand is smooth:
///
/// * **Upper end is a root, lower end is a wall** (`λ = hi − s²`): the
///   substitution kills the singularity at `hi`; the wall at `lo` is not a
///   root of `P`, so `P > 0` there and needs no treatment.
/// * **Both ends are roots** (`λ = lo + (hi−lo) sin²φ`): the substitution kills
///   both endpoint singularities at once, leaving a smooth `2 dφ/√Q` integrand.
///
/// Precompute once per trajectory (`λc` is fixed), then every sample is an
/// `O(log n)` lookup instead of a quadrature call.
pub struct Libration {
    lo: f32,
    hi: f32,
    both_roots: bool,
    xs: Vec<f32>,
    gs: Vec<f32>,
    /// The full half-period `W = w(hi)`.
    pub w_full: f32,
}

impl Libration {
    pub fn new(lo: f32, hi: f32, roots: [f32; 3], knots: usize) -> Self {
        let lo_is_root = roots.iter().any(|&r| (r - lo).abs() < 1e-6);
        if lo_is_root {
            // Both ends are roots: λ = lo + (hi−lo) sin²φ, φ ∈ [0, π/2].
            // P(λ) = (λ−lo)(hi−λ)·Q(λ), dλ/dφ = 2(hi−lo) sinφ cosφ, so
            // dλ/√P = 2 dφ/√Q — smooth.  Q is the remaining linear factor.
            let r = roots
                .iter()
                .copied()
                .find(|&r| (r - lo).abs() > 1e-6 && (r - hi).abs() > 1e-6)
                .unwrap();
            // Sign σ such that P = σ·(λ−lo)(hi−λ)(r−λ).
            let lt = 0.5 * (lo + hi);
            let p_t = (roots[0] - lt) * (roots[1] - lt) * (roots[2] - lt);
            let sigma = p_t / ((lt - lo) * (hi - lt) * (r - lt));
            let q = |lam: f32| sigma * (r - lam);

            let mut xs = Vec::with_capacity(knots);
            let mut gs = Vec::with_capacity(knots);
            let dphi = (std::f32::consts::PI / 2.0) / (knots - 1) as f32;
            let mut acc = 0.0;
            let mut prev_f = 0.0;
            for i in 0..knots {
                let phi = dphi * i as f32;
                let lam = lo + (hi - lo) * phi.sin().powi(2);
                let f = 2.0 / q(lam).abs().sqrt();
                if i == 0 {
                    gs.push(0.0);
                } else {
                    acc += 0.5 * (f + prev_f) * dphi;
                    gs.push(acc);
                }
                xs.push(phi);
                prev_f = f;
            }
            let w_full = *gs.last().unwrap();
            Self {
                lo,
                hi,
                both_roots: true,
                xs,
                gs,
                w_full,
            }
        } else {
            // Upper end is a root, lower end is a wall: λ = hi − s².
            let k = argmin_abs(roots, hi);
            let smax = (hi - lo).max(0.0).sqrt();
            let mut xs = Vec::with_capacity(knots);
            let mut gs = Vec::with_capacity(knots);
            let mut prev_f = 0.0;
            let mut prev_s = 0.0;
            let mut acc = 0.0;
            for i in 0..knots {
                let s = smax * i as f32 / (knots - 1) as f32;
                let lam = hi - s * s;
                let f = 2.0 / p_without(lam, roots, k).abs().sqrt();
                if i == 0 {
                    gs.push(0.0);
                } else {
                    acc += 0.5 * (f + prev_f) * (s - prev_s);
                    gs.push(acc);
                }
                xs.push(s);
                prev_f = f;
                prev_s = s;
            }
            let w_full = *gs.last().unwrap();
            Self {
                lo,
                hi,
                both_roots: false,
                xs,
                gs,
                w_full,
            }
        }
    }

    /// `∫_lo^λ dλ'/√P`.
    pub fn w(&self, lam: f32) -> f32 {
        if self.both_roots {
            let t = ((lam - self.lo) / (self.hi - self.lo)).clamp(0.0, 1.0);
            let phi = t.sqrt().asin();
            interp(phi, &self.xs, &self.gs)
        } else {
            let s = (self.hi - lam).max(0.0).sqrt();
            let g = interp(s, &self.xs, &self.gs);
            self.w_full - g
        }
    }

    /// `∫_λ^hi dλ'/√P`.  This is what the fold formulas want.
    pub fn tail(&self, lam: f32) -> f32 {
        self.w_full - self.w(lam)
    }

    /// Phase `θ ∈ [0, 2π)` on the oval, choosing the increasing/decreasing branch.
    pub fn theta(&self, lam: f32, increasing: bool) -> f32 {
        let t = std::f32::consts::PI * self.w(lam) / self.w_full;
        if increasing {
            t
        } else {
            2.0 * std::f32::consts::PI - t
        }
    }
}

fn argmin_abs(roots: [f32; 3], target: f32) -> usize {
    let mut best = 0;
    let mut best_d = f32::MAX;
    for (j, r) in roots.iter().enumerate() {
        let d = (r - target).abs();
        if d < best_d {
            best_d = d;
            best = j;
        }
    }
    best
}

/// Linear interpolation with clamping to the grid range.
fn interp(x: f32, xs: &[f32], ys: &[f32]) -> f32 {
    let n = xs.len();
    if n == 0 {
        return 0.0;
    }
    if x <= xs[0] {
        return ys[0];
    }
    if x >= xs[n - 1] {
        return ys[n - 1];
    }
    // Binary search for the interval.
    let mut lo = 0usize;
    let mut hi = n - 1;
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if xs[mid] <= x {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let t = (x - xs[lo]) / (xs[hi] - xs[lo]).max(1e-30);
    ys[lo] + t * (ys[hi] - ys[lo])
}

// ---------------------------------------------------------------------------
// Top-level mapping
// ---------------------------------------------------------------------------

/// Per-trajectory cache of the precomputed `Libration`s.  `λc` is fixed along
/// a trajectory, so build the splines once and reuse them for every sample.
#[derive(Default)]
pub struct TorusCache {
    l1: Option<Libration>,
    l2: Option<Libration>,
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

    // Separatrix: the torus degenerates here.  Return a sentinel index.
    if (lc - cf.b).abs() < params.sep_eps {
        return (0.0, 0.0, u32::MAX);
    }

    match params.beta {
        None if lc < cf.b => elliptic_full_ellipse(params, s, lam1, lc, d1),
        Some(_) if lc < cf.b => elliptic_confocal_square(params, s, lam1, lam2, lc, d2),
        _ => hyperbolic(params, s, lam1, lam2, lc, d1, d2),
    }
}

/// Case A — elliptic caustic on a full ellipse (no hyperbola walls): `β = None`.
/// `ν` circulates and the torus splits by the sign of the angular momentum.
fn elliptic_full_ellipse(
    params: &mut TorusParams,
    s: &PhaseSample,
    lam1: f32,
    lc: f32,
    d1: f32,
) -> (f32, f32, u32) {
    let cf = params.confocal;
    let roots = [cf.a, cf.b, lc];
    if params.cache.l1.is_none() {
        params.cache.l1 = Some(Libration::new(params.lam_wall, lc, roots, 512));
    }
    let l1 = params.cache.l1.as_ref().unwrap();
    let th1 = l1.theta(lam1, d1 > 0.0);

    let th2 = nu_angle(s.x, s.y, lam1, cf).rem_euclid(2.0 * std::f32::consts::PI);
    let idx = if s.x * s.vy - s.y * s.vx > 0.0 { 0 } else { 1 };
    (th1, th2, idx)
}

/// Case C — elliptic caustic on a confocal square (`β` wall): folded libration.
/// The accessible region splits into an upper and lower part, split by `sign(y)`.
fn elliptic_confocal_square(
    params: &mut TorusParams,
    s: &PhaseSample,
    lam1: f32,
    lam2: f32,
    lc: f32,
    d2: f32,
) -> (f32, f32, u32) {
    let cf = params.confocal;
    let roots = [cf.a, cf.b, lc];
    let (d1, _bottom) = lam_dots(s, lam1, lam2, cf);
    if params.cache.l1.is_none() {
        params.cache.l1 = Some(Libration::new(params.lam_wall, lc, roots, 512));
    }
    if params.cache.l2.is_none() {
        let beta = params.beta.unwrap();
        params.cache.l2 = Some(Libration::new(beta, cf.a, roots, 512));
    }
    let l1 = params.cache.l1.as_ref().unwrap();
    let l2 = params.cache.l2.as_ref().unwrap();
    let th1 = l1.theta(lam1, d1 > 0.0);

    let tau = if s.x >= 0.0 { 1.0 } else { -1.0 };
    let u2 = l2.w_full - tau * l2.tail(lam2);
    let frac = std::f32::consts::PI * u2 / (2.0 * l2.w_full);
    let th2 = if tau * d2 > 0.0 {
        frac
    } else {
        2.0 * std::f32::consts::PI - frac
    };
    let idx = if s.y >= 0.0 { 0 } else { 1 };
    (th1, th2, idx)
}

/// Cases B / C-hyperbola — hyperbolic caustic: both coordinates fold, and the
/// accessible region is a single torus.
fn hyperbolic(
    params: &mut TorusParams,
    s: &PhaseSample,
    lam1: f32,
    lam2: f32,
    lc: f32,
    d1: f32,
    d2: f32,
) -> (f32, f32, u32) {
    let cf = params.confocal;
    let roots = [cf.a, cf.b, lc];
    if params.cache.l1.is_none() {
        params.cache.l1 = Some(Libration::new(params.lam_wall, cf.b, roots, 512));
    }
    if params.cache.l2.is_none() {
        params.cache.l2 = Some(Libration::new(lc, cf.a, roots, 512));
    }
    let h1 = params.cache.l1.as_ref().unwrap();
    let h2 = params.cache.l2.as_ref().unwrap();

    let sig = if s.y >= 0.0 { 1.0 } else { -1.0 };
    let u1 = h1.w_full + sig * h1.tail(lam1);
    let f1 = std::f32::consts::PI * u1 / (2.0 * h1.w_full);
    let th1 = if -sig * d1 > 0.0 {
        f1
    } else {
        2.0 * std::f32::consts::PI - f1
    };

    let tau = if s.x >= 0.0 { 1.0 } else { -1.0 };
    let u2 = h2.w_full - tau * h2.tail(lam2);
    let f2 = std::f32::consts::PI * u2 / (2.0 * h2.w_full);
    let th2 = if tau * d2 > 0.0 {
        f2
    } else {
        2.0 * std::f32::consts::PI - f2
    };

    (th1, th2, 0)
}

// ---------------------------------------------------------------------------
// Embedding
// ---------------------------------------------------------------------------

/// Standard torus embedding (see §12 of the spec):
/// `((R + r cos θ₁) cos θ₂, (R + r cos θ₁) sin θ₂, r sin θ₁)`.
///
/// `θ₂` is the azimuth about the z-axis (major / toroidal), `θ₁` parametrizes
/// the tube cross-section (minor / poloidal).
pub fn torus_embed(theta1: f32, theta2: f32, r_major: f32, r_minor: f32) -> Vec3 {
    let (sin1, cos1) = theta1.sin_cos();
    let (sin2, cos2) = theta2.sin_cos();
    vec3(
        (r_major + r_minor * cos1) * cos2,
        (r_major + r_minor * cos1) * sin2,
        r_minor * sin1,
    )
}

/// Normal of the standard torus at `(θ₁, θ₂)` for shading.
pub fn torus_normal(theta1: f32, theta2: f32) -> Vec3 {
    let (sin1, cos1) = theta1.sin_cos();
    let (sin2, cos2) = theta2.sin_cos();
    vec3(cos1 * cos2, cos1 * sin2, sin1)
}
