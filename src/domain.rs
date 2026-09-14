use crate::quadratic::ConfocalQuadric;
use crate::torus::ConfocalParams;
use macroquad::prelude::*;

/// A single smooth piece of the boundary.
#[derive(Clone, Debug)]
pub enum Segment {
    /// Straight wall from `a` to `b`.
    Line { a: Vec2, b: Vec2 },
    /// Confocal quadric arc from `a` to `b`.
    /// Only intersections that lie between the two endpoints are accepted.
    Quad {
        curve: ConfocalQuadric,
        a: Vec2,
        b: Vec2,
    },
}

impl Segment {
    pub fn inward_normal(&self, p: Vec2) -> Vec2 {
        match self {
            Segment::Line { a, b } => {
                let d = *b - *a;
                vec2(-d.y, d.x).normalize()
            }
            Segment::Quad { curve, .. } => curve.inward_normal(p),
        }
    }

    pub fn reflect(&self, p: Vec2, dir: Vec2) -> Vec2 {
        match self {
            Segment::Line { a, b } => {
                let d = *b - *a;
                let n = vec2(-d.y, d.x).normalize();
                dir - 2.0 * dir.dot(n) * n
            }
            Segment::Quad { curve, .. } => curve.reflect(p, dir),
        }
    }

    /// Count how many times the ray `p + t·dir` (`t > 0`) crosses this
    /// segment. Used by `Domain::contains`'s parity test, which — unlike
    /// billiard tracing — needs *every* crossing, not just the nearest: a
    /// ray through a domain wall built from a single closed quadric (e.g.
    /// `confocal_ellipse`'s 4 quadrant arcs sharing one curve) can cross that
    /// curve twice, and `Segment::intersect` only ever surfaces the nearer
    /// root.
    fn count_ray_crossings(&self, p: Vec2, dir: Vec2) -> usize {
        match self {
            Segment::Line { .. } => usize::from(
                self.intersect(p, dir)
                    .map(|(t, _, _)| t > 1e-8)
                    .unwrap_or(false),
            ),
            Segment::Quad { curve, a, b } => {
                let (r1, r2) = curve.intersect_roots(p, dir);
                [r1, r2]
                    .into_iter()
                    .flatten()
                    .filter(|&t| t > 1e-8 && on_arc(p + dir * t, *a, *b, curve))
                    .count()
            }
        }
    }

    /// Intersect ray `p + t.dir` with the segment.
    /// Returns `(t, hit, s)` where s is normalised position along the segment.
    pub fn intersect(&self, p: Vec2, dir: Vec2) -> Option<(f32, Vec2, f32)> {
        match self {
            Segment::Line { a, b } => {
                let ab = *b - *a;
                let denom = dir.x * ab.y - dir.y * ab.x;
                if denom.abs() < 1e-12 {
                    return None;
                }
                let ap = *a - p;
                let t = (ap.x * ab.y - ap.y * ab.x) / denom;
                let s = (ap.x * dir.y - ap.y * dir.x) / denom;
                if (0.0..=1.0).contains(&s) && t > 1e-8 {
                    Some((t, p + dir * t, s))
                } else {
                    None
                }
            }
            Segment::Quad { curve, a, b } => {
                let t = curve.intersect(p, dir)?;
                let hit = p + dir * t;
                if !on_arc(hit, *a, *b, curve) {
                    return None;
                }
                let s = arc_frac(hit, *a, *b, curve);
                Some((t, hit, s))
            }
        }
    }
}

/// True if `x` lies between `a` and `b` (inclusive), using the ordered pair.
fn between(x: f32, a: f32, b: f32) -> bool {
    let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
    x >= lo - 1e-4 && x <= hi + 1e-4
}

/// Shift `theta` by a multiple of `2π` so it lies within `π` of `reference`.
///
/// `atan2` only returns the principal range `(-π, π]`, so a conic-angle arc
/// that crosses the branch cut (e.g. a quadrant-III arc of a full ellipse,
/// which runs from `θ=π` down through `θ=-π` to `θ=-π/2`) reads as
/// non-monotonic if the raw `atan2` values are compared directly. Unwrapping
/// every angle to the branch nearest its arc's start point restores
/// monotonicity for any arc under half a turn (true of every quadric arc used
/// in this crate).
fn unwrap_near(theta: f32, reference: f32) -> f32 {
    let two_pi = std::f32::consts::TAU;
    let mut t = theta;
    while t - reference > std::f32::consts::PI {
        t -= two_pi;
    }
    while t - reference < -std::f32::consts::PI {
        t += two_pi;
    }
    t
}

/// Whether `hit` lies on the quadric arc from `a` to `b`, insensitive to the
/// polar angle from the origin.
///
/// For a hyperbola (`b < λ < a`) we restrict to the branch of the endpoints
/// (sign(x)) and use that `y` is monotone along a branch.  For an ellipse
/// (`λ < b`) we use the conic parameter angle `θ` with
/// `x = √(a−λ) cos θ, y = √(b−λ) sin θ` — the polar angle is *not* an affine
/// function of `θ` unless `a−λ = b−λ`, which is why the old origin-angle test
/// was wrong.  Both endpoints' coordinates are single-quadrant for the tables
/// in this crate; when they straddle quadrants the range test below still
/// holds because the arc spans monotonically.
fn on_arc(hit: Vec2, a: Vec2, b: Vec2, curve: &ConfocalQuadric) -> bool {
    if curve.is_hyperbola() {
        // Same branch, then y monotone along the branch.
        return curve.branch_sign(hit) == curve.branch_sign(a) && between(hit.y, a.y, b.y);
    }
    // Ellipse: conic angle θ.  `x = √(a−λ) cos θ, y = √(b−λ) sin θ`.  θ is
    // monotonic along any arc under half a turn, so a plain range test is
    // exact — but `atan2` only returns `(-π, π]`, so an arc crossing that
    // branch cut (e.g. a quadrant-III arc of a full ellipse) needs `b`/`hit`
    // unwrapped onto `a`'s branch first (`unwrap_near`) before comparing.
    // (The polar angle from the origin is NOT an affine function of θ, which
    // is why the old origin-angle test misread in-quadrant arcs.)
    let theta = |p: Vec2| {
        let ca = (curve.a_param - curve.lambda).max(1e-30).sqrt();
        let cb = (curve.b_param - curve.lambda).max(1e-30).sqrt();
        (p.y / cb).atan2(p.x / ca)
    };
    let ta = theta(a);
    let tb = unwrap_near(theta(b), ta);
    let th = unwrap_near(theta(hit), ta);
    between(th, ta, tb)
}

/// Normalised position `s ∈ [0, 1]` along the arc `a → b`.
fn arc_frac(hit: Vec2, a: Vec2, b: Vec2, curve: &ConfocalQuadric) -> f32 {
    if curve.is_hyperbola() {
        // y is monotone along a branch.
        if (b.y - a.y).abs() < 1e-9 {
            0.5
        } else {
            ((hit.y - a.y) / (b.y - a.y)).clamp(0.0, 1.0)
        }
    } else {
        // Ellipse: monotone in conic angle θ.
        let theta = |p: Vec2| {
            let ca = (curve.a_param - curve.lambda).max(1e-30).sqrt();
            let cb = (curve.b_param - curve.lambda).max(1e-30).sqrt();
            (p.y / cb).atan2(p.x / ca)
        };
        let ta = theta(a);
        let tb = unwrap_near(theta(b), ta);
        let th = unwrap_near(theta(hit), ta);
        if (ta - tb).abs() < 1e-9 {
            0.5
        } else {
            ((th - ta) / (tb - ta)).clamp(0.0, 1.0)
        }
    }
}

// ---------------------------------------------------------------------------
// Domain
// ---------------------------------------------------------------------------
#[derive(Clone)]
pub struct Domain {
    pub segments: Vec<Segment>,
}

impl Domain {
    pub fn new(segments: Vec<Segment>) -> Self {
        Self { segments }
    }

    pub fn corners(&self) -> Vec<Vec2> {
        self.segments
            .windows(2)
            .map(|w| match &w[0] {
                Segment::Line { b, .. } => *b,
                Segment::Quad { b, .. } => *b,
            })
            .collect()
    }

    pub fn contains(&self, point: Vec2) -> bool {
        // Try multiple ray directions for robustness in non-convex domains.
        // Single-direction ray casting can miss boundaries when the ray
        // grazes past a vertex (e.g. L-shape's top-right corner).
        for &ray_dir in &[
            vec2(1.0, 0.0),
            vec2(0.0, -1.0),
            vec2(0.5, -0.866),
            vec2(-0.5, -0.866),
        ] {
            let mut crossings = 0;
            for seg in &self.segments {
                crossings += seg.count_ray_crossings(point, ray_dir);
            }
            if crossings % 2 == 1 {
                return true;
            }
        }
        false
    }

    pub fn intersect(&self, p: Vec2, dir: Vec2) -> Option<(f32, usize, Vec2)> {
        let mut best: Option<(f32, usize, Vec2)> = None;
        for (i, seg) in self.segments.iter().enumerate() {
            if let Some((t, hit, _)) = seg.intersect(p, dir) {
                match &best {
                    Some((best_t, _, _)) if t >= *best_t => {}
                    _ => best = Some((t, i, hit)),
                }
            }
        }
        best
    }

    pub fn reflect(&self, p: Vec2, dir: Vec2, segment_idx: usize) -> Vec2 {
        // π/2 corner handling
        let corners: Vec<(Vec2, f32)> = self
            .segments
            .windows(2)
            .map(|w| {
                let corner = match &w[0] {
                    Segment::Line { b, .. } => *b,
                    Segment::Quad { b, .. } => *b,
                };
                let n1 = w[0].inward_normal(corner);
                let n2 = w[1].inward_normal(corner);
                let angle = n1.angle_between(n2).abs();
                (corner, angle)
            })
            .collect();

        let eps = 1e-4;
        for (corner, angle) in &corners {
            if p.distance(*corner) < eps {
                if (*angle - std::f32::consts::PI / 2.0).abs() < 0.1 {
                    return -dir;
                }
                break;
            }
        }
        self.segments[segment_idx].reflect(p, dir)
    }

    pub fn trace(&self, mut p: Vec2, mut v: Vec2, max_steps: usize) -> Vec<(Vec2, Vec2)> {
        let mut segs = Vec::with_capacity(max_steps);
        for _ in 0..max_steps {
            let speed = v.length();
            if speed < 1e-12 {
                break;
            }
            let dir = v / speed;
            let (_t, idx, hit) = match self.intersect(p, dir) {
                Some(r) => r,
                None => break,
            };
            segs.push((p, hit));
            v = self.reflect(hit, v, idx);
            p = hit + 1e-4 * v.normalize();
        }
        segs
    }

    /// Sample boundary points for drawing. Returns them in order.
    pub fn sample_boundary(&self, n: usize) -> Vec<Vec2> {
        let mut pts = Vec::new();
        for seg in &self.segments {
            match seg {
                Segment::Line { a, b } => {
                    pts.push(*a);
                    pts.push(*b);
                }
                Segment::Quad { curve, a, b } => {
                    sample_quadric_arc(curve, *a, *b, n)
                        .into_iter()
                        .for_each(|p| pts.push(p));
                }
            }
        }
        pts
    }
}

/// Sample points along a confocal-quadric arc from `a` to `b`, parameterized
/// by the **conic parameter** rather than the polar angle from the origin.
///
/// For an ellipse the conic angle `θ` (with `x = √(a−λ) cos θ`, `y = √(b−λ) sin
/// θ`) is monotonic along a single-quadrant arc; for a hyperbola the branch's
/// `y` is monotonic and `x` follows the branch.  This renders the arc exactly
/// where `intersect`/`on_arc` put it, instead of casting origin-rays (which hit
/// the wrong hyperbola branch or fill the ellipse interior).
/// Sample a confocal-quadric wall arc of the domain (`curve` between endpoints
/// `a` and `b`) into `n` points (inclusive of both ends).  Public so critical-
/// layer tools can slide a trajectory exactly along a boundary wall.
pub fn sample_wall_arc(curve: &ConfocalQuadric, a: Vec2, b: Vec2, n: usize) -> Vec<Vec2> {
    sample_quadric_arc(curve, a, b, n)
}

fn sample_quadric_arc(curve: &ConfocalQuadric, a: Vec2, b: Vec2, n: usize) -> Vec<Vec2> {
    let lam = curve.lambda;
    let ca = (curve.a_param - lam).max(1e-30);
    if curve.is_hyperbola() {
        // Branch x≥0 (or x<0): x²/(a−λ) − y²/(λ−b) = 1,
        // so x = ±√(ca)·√(1 + y²/(λ−b)), y monotone between endpoints.
        let d = (lam - curve.b_param).max(1e-30);
        let sign = if a.x >= 0.0 { 1.0 } else { -1.0 };
        let (y0, y1) = (a.y, b.y);
        let steps = n.max(3);
        let mut pts = Vec::with_capacity(steps + 1);
        for i in 0..=steps {
            let t = i as f32 / steps as f32;
            let y = y0 + (y1 - y0) * t;
            let x = sign * (ca * (1.0 + y * y / d)).max(0.0).sqrt();
            pts.push(vec2(x, y));
        }
        pts
    } else {
        // Ellipse: conic angle θ.
        let cb = (curve.b_param - lam).max(1e-30);
        let theta = |p: Vec2| (p.y / cb.sqrt()).atan2(p.x / ca.sqrt());
        let t0 = theta(a);
        let t1 = unwrap_near(theta(b), t0);
        let steps = n.max(3);
        let mut pts = Vec::with_capacity(steps + 1);
        for i in 0..=steps {
            let f = i as f32 / steps as f32;
            let th = t0 + (t1 - t0) * f;
            pts.push(vec2(ca.sqrt() * th.cos(), cb.sqrt() * th.sin()));
        }
        pts
    }
}

// ---------------------------------------------------------------------------
// Convenience constructors
// ---------------------------------------------------------------------------

/// Build a confocal quadrilateral domain from an ellipse (λ₁) and
/// a hyperbola (λ₂) of the same confocal family `cf`.
pub fn confocal_quad(cf: ConfocalParams, lambda_ell: f32, lambda_hyp: f32) -> Domain {
    let arcs = crate::quadratic::quadrilateral_arcs(cf, lambda_ell, lambda_hyp);
    Domain::new(
        arcs.iter()
            .map(|(from, to, curve)| Segment::Quad {
                curve: *curve,
                a: *from,
                b: *to,
            })
            .collect(),
    )
}

/// Build a full-ellipse domain: the single confocal ellipse `λ = lambda_ell`
/// (no hyperbola walls).  This is the classical confocal-ellipse billiard whose
/// phase manifold is two tori below the focal separatrix and one torus above
/// (the A–B–A molecule).
pub fn confocal_ellipse(cf: ConfocalParams, lambda_ell: f32) -> Domain {
    let ell = crate::quadratic::confocal(cf, lambda_ell);
    // The full ellipse is a single closed quadric arc; sample its 4 quadrant
    // arcs so `intersect`/`reflect` see a closed boundary.
    let a = cf.a;
    let b = cf.b;
    let ca = (a - lambda_ell).max(1e-30).sqrt();
    let cb = (b - lambda_ell).max(1e-30).sqrt();
    let pts = [
        vec2(ca, 0.0),  // right
        vec2(0.0, cb),  // top
        vec2(-ca, 0.0), // left
        vec2(0.0, -cb), // bottom
    ];
    let quad = |from: Vec2, to: Vec2| Segment::Quad {
        curve: ell,
        a: from,
        b: to,
    };
    Domain::new(vec![
        quad(pts[0], pts[1]),
        quad(pts[1], pts[2]),
        quad(pts[2], pts[3]),
        quad(pts[3], pts[0]),
    ])
}

/// L-shape built from confocal quadrics with a re-entrant 270° corner.
///
/// The shape is bounded by 6 confocal arcs in CCW order:
///   1. **Top** — outer ellipse from right-hyperbola to left-upper hyperbola
///   2. **Left-upper** — left-upper hyperbola from outer ellipse to step ellipse
///   3. **Step** — step ellipse from left-upper hyperbola to left-lower hyperbola
///   4. **Left-lower** — left-lower hyperbola from step ellipse to outer ellipse
///   5. **Bottom** — outer ellipse from left-lower hyperbola to right hyperbola
///   6. **Right** — right hyperbola from outer ellipse back to top
///
/// The notch (step) on the left side is controlled by `lambda_ell_step`
/// (an ellipse smaller than the outer one) and the two left-side hyperbolas,
/// creating a 270° re-entrant corner at the junction of arcs 3 and 4.
pub fn confocal_lshape(
    cf: ConfocalParams,
    lambda_ell_outer: f32,
    lambda_hyp_right: f32,
    lambda_hyp_left_upper: f32,
    lambda_hyp_left_lower: f32,
    lambda_ell_step: f32,
) -> Domain {
    let a = cf.a;
    let b = cf.b;
    let hyp_left_upper = lambda_hyp_left_upper.max(lambda_hyp_left_lower);
    let hyp_left_lower = if lambda_hyp_left_upper == lambda_hyp_left_lower {
        lambda_hyp_left_upper
    } else {
        lambda_hyp_left_lower
    };
    use crate::quadratic::ConfocalQuadric;

    // --- quadrics -----------------------------------------------------------
    let ell_outer = ConfocalQuadric {
        a_param: a,
        b_param: b,
        lambda: lambda_ell_outer,
    };
    let hyp_right = ConfocalQuadric {
        a_param: a,
        b_param: b,
        lambda: lambda_hyp_right,
    };
    let hyp_left_upper_q = ConfocalQuadric {
        a_param: a,
        b_param: b,
        lambda: hyp_left_upper,
    };
    let hyp_left_lower_q = ConfocalQuadric {
        a_param: a,
        b_param: b,
        lambda: hyp_left_lower,
    };
    let ell_step = ConfocalQuadric {
        a_param: a,
        b_param: b,
        lambda: lambda_ell_step,
    };

    // --- intersection points ------------------------------------------------
    // intersections returns [tr, br, bl, tl] where each is (±x, ±y).
    let [tr_or, br_or, _bl_or, _tl_or] = ConfocalQuadric::intersections(&ell_outer, &hyp_right)
        .expect("outer ellipse and right hyperbola must intersect");

    let [_tr_olu, _br_olu, _bl_olu, tl_olu] =
        ConfocalQuadric::intersections(&ell_outer, &hyp_left_upper_q)
            .expect("outer ellipse and left-upper hyperbola must intersect");

    let [_tr_slu, _br_slu, _bl_slu, tl_slu] =
        ConfocalQuadric::intersections(&ell_step, &hyp_left_upper_q)
            .expect("step ellipse and left-upper hyperbola must intersect");

    let [_tr_sll, _br_sll, bl_sll, _tl_sll] =
        ConfocalQuadric::intersections(&ell_step, &hyp_left_lower_q)
            .expect("step ellipse and left-lower hyperbola must intersect");

    let [_tr_oll, _br_oll, bl_oll, _tl_oll] =
        ConfocalQuadric::intersections(&ell_outer, &hyp_left_lower_q)
            .expect("outer ellipse and left-lower hyperbola must intersect");

    // --- 6 arcs in CCW order ------------------------------------------------
    Domain::new(vec![
        // 1. Top: outer ellipse, right-hyp → left-upper-hyp
        Segment::Quad {
            curve: ell_outer,
            a: tr_or,
            b: tl_olu,
        },
        // 2. Left-upper: left-upper hyperbola, outer-ellipse → step-ellipse
        Segment::Quad {
            curve: hyp_left_upper_q,
            a: tl_olu,
            b: tl_slu,
        },
        // 3. Step: step ellipse, left-upper-hyp → left-lower-hyp.
        //    Goes CCW from upper-left (tl_slu) to lower-left (bl_sll)
        //    through the left side of the step ellipse.
        Segment::Quad {
            curve: ell_step,
            a: tl_slu,
            b: bl_sll,
        },
        // 4. Left-lower: left-lower hyperbola, step-ellipse → outer-ellipse.
        //    Goes from the step-ellipse intersection (bl_sll) down to
        //    the outer-ellipse intersection (bl_oll), along the hyperbola.
        Segment::Quad {
            curve: hyp_left_lower_q,
            a: bl_sll,
            b: bl_oll,
        },
        // 5. Bottom: outer ellipse, left-lower-hyp → right-hyp
        Segment::Quad {
            curve: ell_outer,
            a: bl_oll,
            b: br_or,
        },
        // 6. Right: right hyperbola, outer-ellipse back to top
        Segment::Quad {
            curve: hyp_right,
            a: br_or,
            b: tr_or,
        },
    ])
}

/// Standard L-shape from the pseudo-integrable doc (§2): the region
/// `([0, α₁] × [β₁, β₃]) ∪ ([0, α₂] × [β₁, β₂])` in the (λ₁, λ₂) chart, with
/// `0 < α₁ < α₂ < b < β₁ < β₂ < β₃ < a`.
///
/// An **in-quadrant** table (strictly inside x>0, y>0) so the λ-chart is 1:1
/// (no fold labels needed).  Six corners: five convex, one reflex at
/// `(λ₁, λ₂) = (α₁, β₂)`.
///
/// Walls (CCW):
///   1. outer ellipse λ₁=0 (tall leg's left),
///   2. hyperbola β₁ (base bottom),
///   3. ellipse α₂ (base right),
///   4. hyperbola β₂ (step — top of the base),
///   5. ellipse α₁ (tall leg's right),
///   6. hyperbola β₃ (top).
pub fn confocal_lshape_standard(
    cf: ConfocalParams,
    alpha1: f32,
    alpha2: f32,
    beta1: f32,
    beta2: f32,
    beta3: f32,
) -> Domain {
    let a = cf.a;
    let b = cf.b;
    use crate::quadratic::ConfocalQuadric;

    let q = |lam| ConfocalQuadric {
        a_param: a,
        b_param: b,
        lambda: lam,
    };
    let ell_outer = q(0.0);
    let ell_alpha1 = q(alpha1);
    let ell_alpha2 = q(alpha2);
    let hyp_beta1 = q(beta1);
    let hyp_beta2 = q(beta2);
    let hyp_beta3 = q(beta3);

    // Corner at the intersection of an ellipse and a hyperbola: take the
    // top-right image (x>0, y>0), which is the one in the first quadrant.
    let corner = |ell: &ConfocalQuadric, hyp: &ConfocalQuadric| {
        ConfocalQuadric::intersections(ell, hyp).expect("confocal intersection")[0]
    };

    // Six corners, labelled by the (ellipse, hyperbola) walls meeting there.
    let tl0 = corner(&ell_outer, &hyp_beta1); // (λ₁, λ₂) = (0, β₁)
    let br_base = corner(&ell_alpha2, &hyp_beta1); // (α₂, β₁)
    let tr_base = corner(&ell_alpha2, &hyp_beta2); // (α₂, β₂)
    let reflex = corner(&ell_alpha1, &hyp_beta2); // (α₁, β₂) — reflex corner
    let tr_tall = corner(&ell_alpha1, &hyp_beta3); // (α₁, β₃)
    let tl_tall = corner(&ell_outer, &hyp_beta3); // (0, β₃)

    Domain::new(vec![
        // 1. Hyperbola β₁, outer ellipse → α₂-ellipse (base bottom)
        Segment::Quad {
            curve: hyp_beta1,
            a: tl0,
            b: br_base,
        },
        // 2. Ellipse α₂, from β₁ up to β₂ (base right)
        Segment::Quad {
            curve: ell_alpha2,
            a: br_base,
            b: tr_base,
        },
        // 3. Hyperbola β₂, from (α₂,β₂) back to (α₁,β₂) (step)
        Segment::Quad {
            curve: hyp_beta2,
            a: tr_base,
            b: reflex,
        },
        // 4. Ellipse α₁, from β₂ up to β₃ (tall leg right)
        Segment::Quad {
            curve: ell_alpha1,
            a: reflex,
            b: tr_tall,
        },
        // 5. Hyperbola β₃, from (α₁,β₃) to (0,β₃) (top)
        Segment::Quad {
            curve: hyp_beta3,
            a: tr_tall,
            b: tl_tall,
        },
        // 6. Outer ellipse, from (0,β₃) down to (0,β₁) (left wall, closes)
        Segment::Quad {
            curve: ell_outer,
            a: tl_tall,
            b: tl0,
        },
    ])
}

/// Simple axis-aligned square domain.
pub fn square() -> Domain {
    Domain::new(vec![
        Segment::Line {
            a: vec2(-1.0, -1.0),
            b: vec2(1.0, -1.0),
        },
        Segment::Line {
            a: vec2(1.0, -1.0),
            b: vec2(1.0, 1.0),
        },
        Segment::Line {
            a: vec2(1.0, 1.0),
            b: vec2(-1.0, 1.0),
        },
        Segment::Line {
            a: vec2(-1.0, 1.0),
            b: vec2(-1.0, -1.0),
        },
    ])
}

/// L-shape polyline.
pub fn lshape_poly() -> Domain {
    Domain::new(vec![
        Segment::Line {
            a: vec2(0.0, 0.0),
            b: vec2(2.0, 0.0),
        },
        Segment::Line {
            a: vec2(2.0, 0.0),
            b: vec2(2.0, 1.0),
        },
        Segment::Line {
            a: vec2(2.0, 1.0),
            b: vec2(1.0, 1.0),
        },
        Segment::Line {
            a: vec2(1.0, 1.0),
            b: vec2(1.0, 2.0),
        },
        Segment::Line {
            a: vec2(1.0, 2.0),
            b: vec2(0.0, 2.0),
        },
        Segment::Line {
            a: vec2(0.0, 2.0),
            b: vec2(0.0, 0.0),
        },
    ])
}
