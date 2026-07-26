use crate::quadratic::ConfocalQuadric;
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
                if s >= 0.0 && s <= 1.0 && t > 1e-8 {
                    Some((t, p + dir * t, s))
                } else {
                    None
                }
            }
            Segment::Quad { curve, a, b } => {
                let t = curve.intersect(p, dir)?;
                let hit = p + dir * t;
                if !between_angles(hit, *a, *b, curve.centre()) {
                    return None;
                }
                let centre = curve.centre();
                let ha = (hit - centre).y.atan2((hit - centre).x);
                let aa = (*a - centre).y.atan2((*a - centre).x);
                let ba = (*b - centre).y.atan2((*b - centre).x);
                Some((t, hit, angle_frac(ha, aa, ba)))
            }
        }
    }
}

fn norm_angle(mut ang: f32) -> f32 {
    while ang < 0.0 {
        ang += 2.0 * std::f32::consts::PI;
    }
    while ang >= 2.0 * std::f32::consts::PI {
        ang -= 2.0 * std::f32::consts::PI;
    }
    ang
}

/// Check whether `hit` lies on the shorter arc from `a` to `b` (CCW)
/// as seen from `centre`.
fn between_angles(hit: Vec2, a: Vec2, b: Vec2, centre: Vec2) -> bool {
    let ha = norm_angle((hit - centre).y.atan2((hit - centre).x));
    let aa = norm_angle((a - centre).y.atan2((a - centre).x));
    let ba = norm_angle((b - centre).y.atan2((b - centre).x));

    if ba > aa {
        ha >= aa - 0.05 && ha <= ba + 0.05
    } else {
        ha >= aa - 0.05 || ha <= ba + 0.05
    }
}

/// Fraction along the arc from aa to ba (CCW).
fn angle_frac(mut ha: f32, aa: f32, mut ba: f32) -> f32 {
    ha = norm_angle(ha);
    let aa = norm_angle(aa);
    ba = norm_angle(ba);

    if ba > aa {
        (ha - aa) / (ba - aa)
    } else {
        let ha_wrapped = if ha < aa {
            ha + 2.0 * std::f32::consts::PI
        } else {
            ha
        };
        (ha_wrapped - aa) / (ba + 2.0 * std::f32::consts::PI - aa)
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
            .filter_map(|w| match &w[0] {
                Segment::Line { b, .. } => Some(*b),
                Segment::Quad { b, .. } => Some(*b),
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
                if let Some((t, _, _)) = seg.intersect(point, ray_dir) {
                    if t > 1e-8 {
                        crossings += 1;
                    }
                }
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
            .filter_map(|w| {
                let corner = match &w[0] {
                    Segment::Line { b, .. } => *b,
                    Segment::Quad { b, .. } => *b,
                };
                let n1 = w[0].inward_normal(corner);
                let n2 = w[1].inward_normal(corner);
                let angle = n1.angle_between(n2).abs();
                Some((corner, angle))
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
                    let centre = curve.centre();
                    let aa = norm_angle((*a - centre).y.atan2((*a - centre).x));
                    let ba = norm_angle((*b - centre).y.atan2((*b - centre).x));
                    let mut delta = ba - aa;
                    if delta <= 0.0 {
                        delta += 2.0 * std::f32::consts::PI;
                    }
                    let steps = (delta / std::f32::consts::PI * n as f32).max(3.0) as usize;

                    for i in 0..=steps {
                        let frac = i as f32 / steps as f32;
                        let angle = aa + frac * delta;
                        let dir = vec2(angle.cos(), angle.sin());
                        if let Some(t) = curve.intersect(centre + 0.001 * dir, dir) {
                            pts.push(centre + dir * t);
                        }
                    }
                }
            }
        }
        pts
    }
}

// ---------------------------------------------------------------------------
// Convenience constructors
// ---------------------------------------------------------------------------

/// Build a confocal quadrilateral domain from an ellipse (λ₁) and
/// a hyperbola (λ₂) of the same confocal family (a, b).
pub fn confocal_quad(a: f32, b: f32, lambda_ell: f32, lambda_hyp: f32) -> Domain {
    let arcs = crate::quadratic::quadrilateral_arcs(a, b, lambda_ell, lambda_hyp);
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
    a: f32,
    b: f32,
    lambda_ell_outer: f32,
    lambda_hyp_right: f32,
    lambda_hyp_left_upper: f32,
    lambda_hyp_left_lower: f32,
    lambda_ell_step: f32,
) -> Domain {
    // left-upper & left-lower may be equal (single left hyperbola)
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
