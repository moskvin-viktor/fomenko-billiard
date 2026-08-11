//! Confocal coordinates and phase-space helpers.
//!
//! The confocal family `(a, b)` and a phase-space sample `(x, y, vx, vy)` are
//! the two values threaded through the torus mapping; the functions here convert
//! between them and the confocal-coordinate view `(λ₁, λ₂)` used by
//! [`crate::torus::map`].

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
