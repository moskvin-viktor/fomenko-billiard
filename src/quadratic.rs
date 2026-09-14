use crate::torus::ConfocalParams;
use macroquad::prelude::*;

/// Confocal quadric family:
///
///   (b − λ)·x² + (a − λ)·y² = (a − λ)·(b − λ),   λ ≤ a
///
/// where ∞ > a > b > 0 are fixed constants.
/// The foci are at (±c, 0) with c² = a − b.
///
/// At λ = 0 we get the base ellipse x²/a + y²/b = 1.
/// At λ = b we get the degenerate union of:
///   - the segment [−c, +c] on the x-axis ("degenerate ellipse"), and
///   - the two horizontal rays (−∞, −c] ∪ [+c, ∞) ("degenerate hyperbola").
/// At λ = a we get the vertical segment x = 0, −√b ≤ y ≤ √b.
#[derive(Clone, Copy, Debug)]
pub struct ConfocalQuadric {
    pub a_param: f32, // a
    pub b_param: f32, // b
    pub lambda: f32,  // λ
}

impl ConfocalQuadric {
    /// Semi-focal distance c = √(a − b).
    #[allow(dead_code)]
    pub fn c(&self) -> f32 {
        (self.a_param - self.b_param).sqrt()
    }

    /// Evaluate Q(x, y).  The interior (inside the billiard) is Q(p) < 0.
    #[allow(dead_code)]
    pub fn eval(&self, p: Vec2) -> f32 {
        let a = self.a_param;
        let b = self.b_param;
        let lam = self.lambda;
        (b - lam) * p.x * p.x + (a - lam) * p.y * p.y - (a - lam) * (b - lam)
    }

    /// Gradient ∇Q.
    pub fn grad(&self, p: Vec2) -> Vec2 {
        let a = self.a_param;
        let b = self.b_param;
        let lam = self.lambda;
        vec2(2.0 * (b - lam) * p.x, 2.0 * (a - lam) * p.y)
    }

    pub fn inward_normal(&self, p: Vec2) -> Vec2 {
        -self.grad(p).normalize()
    }

    pub fn reflect(&self, p: Vec2, dir: Vec2) -> Vec2 {
        let n = self.inward_normal(p);
        dir - 2.0 * dir.dot(n) * n
    }

    /// Smallest positive t where ray p + t·dir hits Q = 0.
    pub fn intersect(&self, p: Vec2, dir: Vec2) -> Option<f32> {
        let (t1, t2) = self.intersect_roots(p, dir);
        t1.or(t2)
    }

    /// Both forward (`t > 0`) roots where ray `p + t·dir` hits `Q = 0`, in
    /// ascending order. A ray through a full closed quadric (e.g. the single
    /// ellipse `confocal_ellipse` splits into 4 quadrant arcs) crosses the
    /// *same* underlying curve twice; `intersect` only ever returns the
    /// nearer one, which is right for billiard tracing (the ball is always
    /// inside, so the near wall is the physical hit) but wrong for
    /// ray-casting containment tests, which need every crossing to get
    /// parity right.
    pub fn intersect_roots(&self, p: Vec2, dir: Vec2) -> (Option<f32>, Option<f32>) {
        let a = self.a_param;
        let b = self.b_param;
        let lam = self.lambda;
        let bx = b - lam;
        let ay = a - lam;

        let qa = bx * dir.x * dir.x + ay * dir.y * dir.y;
        let qb = 2.0 * (bx * p.x * dir.x + ay * p.y * dir.y);
        let qc = bx * p.x * p.x + ay * p.y * p.y - ay * bx;

        let eps = 1e-6;
        if qa.abs() < 1e-12 {
            if qb.abs() < 1e-12 {
                return (None, None);
            }
            let t = -qc / qb;
            return if t > eps { (Some(t), None) } else { (None, None) };
        }

        let disc = qb * qb - 4.0 * qa * qc;
        if disc < 0.0 {
            return (None, None);
        }
        let sd = disc.sqrt();
        let mut t1 = (-qb - sd) / (2.0 * qa);
        let mut t2 = (-qb + sd) / (2.0 * qa);
        if t1 > t2 {
            std::mem::swap(&mut t1, &mut t2);
        }
        let r1 = if t1 > eps { Some(t1) } else { None };
        let r2 = if t2 > eps { Some(t2) } else { None };
        (r1, r2)
    }

    /// Centre of the quadric (∇Q = 0).
    pub fn centre(&self) -> Vec2 {
        vec2(0.0, 0.0)
    }

    /// Whether this quadric is a hyperbola (`b < λ < a`), which has two
    /// branches (x > 0 and x < 0).  Ellipses (`λ < b`) are single-branched.
    pub fn is_hyperbola(&self) -> bool {
        self.lambda > self.b_param
    }

    /// The branch sign `sign(x)` of a point on this quadric.  For a hyperbola
    /// this selects which of the two sheets the point is on; for an ellipse it
    /// is always `+` (single branch).
    pub fn branch_sign(&self, p: Vec2) -> i8 {
        if self.is_hyperbola() {
            if p.x >= 0.0 {
                1
            } else {
                -1
            }
        } else {
            1
        }
    }

    /// Intersection points of two confocal quadrics with the same a, b.
    /// Returns 4 points [tr, br, bl, tl] when they exist.
    pub fn intersections(l1: &Self, l2: &Self) -> Option<[Vec2; 4]> {
        let a = l1.a_param;
        let b = l1.b_param;
        let lam1 = l1.lambda;
        let lam2 = l2.lambda;

        // Solve the linear system in x², y²:
        //   (b-λ₁)x² + (a-λ₁)y² = r₁,   rᵢ = (a-λᵢ)(b-λᵢ)
        //   (b-λ₂)x² + (a-λ₂)y² = r₂
        //
        // In matrix form:
        //   [ b-λ₁  a-λ₁ ] [x²]   [ r₁ ]
        //   [ b-λ₂  a-λ₂ ] [y²] = [ r₂ ]
        //
        // det = (b-λ₁)(a-λ₂) - (b-λ₂)(a-λ₁) = (a-b)(λ₂-λ₁)
        //
        // Cramer's rule:
        //   x² = (r₁·(a-λ₂) - r₂·(a-λ₁)) / det
        //      = (a-λ₁)(a-λ₂)(λ₂-λ₁) / ((a-b)(λ₂-λ₁))
        //      = (a-λ₁)(a-λ₂) / (a-b)
        //
        //   y² = ((b-λ₁)·r₂ - (b-λ₂)·r₁) / det
        //      = (b-λ₁)(b-λ₂)(λ₁-λ₂) / ((a-b)(λ₂-λ₁))
        //      = -(b-λ₁)(b-λ₂) / (a-b)
        //      = (b-λ₁)(b-λ₂) / (b-a)

        let x2 = ((a - lam1) * (a - lam2)) / (a - b);
        let y2 = ((b - lam1) * (b - lam2)) / (b - a);

        if x2 < 0.0 || y2 < 0.0 {
            return None;
        }

        let x = x2.sqrt();
        let y = y2.sqrt();

        Some([
            vec2(x, y),   // tr
            vec2(x, -y),  // br
            vec2(-x, -y), // bl
            vec2(-x, y),  // tl
        ])
    }

    /// Sample points on this quadric by casting rays from the origin.
    #[allow(dead_code)]
    pub fn sample_boundary(&self, n: usize) -> Vec<Vec2> {
        (0..n)
            .filter_map(|i| {
                let angle = 2.0 * std::f32::consts::PI * i as f32 / n as f32;
                let dir = vec2(angle.cos(), angle.sin());
                self.intersect(vec2(0.0, 0.0), dir).map(|t| dir * t)
            })
            .collect()
    }

    /// Velocity tangent to the caustic quadric Q_Λ(x,y) = 0 at point p.
    ///
    /// The tangent direction at p is perpendicular to ∇Q_Λ(p), so:
    ///   v = ± normalize( ∇Q_Λ(p) × (0,0,1) )  =  ± normalize( -∂Q/∂y, ∂Q/∂x )
    ///
    /// Gradient of the caustic: (2(b-Λ)x, 2(a-Λ)y)
    /// Tangent: (-(a-Λ)y, (b-Λ)x)  — rotates gradient 90° CW.
    ///
    /// We pick the sign closest to `hint_dir`.
    pub fn velocity_from_caustic(cf: ConfocalParams, p: Vec2, lambda: f32, hint_dir: Vec2) -> Vec2 {
        let a = cf.a;
        let b = cf.b;
        // ∇Q_Λ(p) = (2(b-Λ)x, 2(a-Λ)y)
        // Tangent (perpendicular, rotated 90° CW): (-(a-Λ)y, (b-Λ)x)
        let t = vec2(-(a - lambda) * p.y, (b - lambda) * p.x);
        let t = t.normalize();

        // Two possible directions: ±t
        // Pick the one closer to hint_dir
        let hint = hint_dir.normalize();
        if t.dot(hint) >= 0.0 {
            t
        } else {
            -t
        }
    }
}

// ---------------------------------------------------------------------------
// Convenience constructors
// ---------------------------------------------------------------------------

/// Create a confocal quadric: (b − λ)x² + (a − λ)y² = (a − λ)(b − λ).
pub fn confocal(cf: ConfocalParams, lambda: f32) -> ConfocalQuadric {
    ConfocalQuadric {
        a_param: cf.a,
        b_param: cf.b,
        lambda,
    }
}

/// Given an ellipse (λ₁) and a hyperbola (λ₂) from the same confocal family,
/// return 4 arcs `[(from, to, quadric)]` in CCW order:
///   0 = top (ellipse, right → left), 1 = left (hyperbola, top → bottom),
///   2 = bottom (ellipse, left → right), 3 = right (hyperbola, bottom → top).
pub fn quadrilateral_arcs(
    cf: ConfocalParams,
    lambda_ell: f32,
    lambda_hyp: f32,
) -> Vec<(Vec2, Vec2, ConfocalQuadric)> {
    let ell = confocal(cf, lambda_ell);
    let hyp = confocal(cf, lambda_hyp);

    let pts = ConfocalQuadric::intersections(&ell, &hyp)
        .expect("The two confocal quadrics must intersect");

    let [tr, br, bl, tl] = pts;

    vec![
        (tr, tl, ell), // top: ellipse right → left
        (tl, bl, hyp), // left: hyperbola top → bottom
        (bl, br, ell), // bottom: ellipse left → right
        (br, tr, hyp), // right: hyperbola bottom → top
    ]
}
