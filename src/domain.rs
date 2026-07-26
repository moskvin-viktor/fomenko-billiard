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
        let ray_dir = vec2(1.0, 0.0);
        let mut crossings = 0;
        for seg in &self.segments {
            if let Some((t, _, _)) = seg.intersect(point, ray_dir) {
                if t > 1e-8 {
                    crossings += 1;
                }
            }
        }
        crossings % 2 == 1
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

/// L-shape built from confocal quadrics.
/// Uses 6 arcs: top/bottom ellipse (λ₁), right hyperbola (λ₂),
/// and splits the left hyperbola into two with λ₃ (upper) and λ₄ (lower),
/// creating a rectangular indentation.
#[allow(dead_code)]
pub fn confocal_lshape(
    a: f32,
    b: f32,
    lambda_ell: f32,
    lambda_hyp_right: f32,
    _lambda_hyp_left_upper: f32,
    _lambda_hyp_left_lower: f32,
) -> Domain {
    use crate::quadratic::ConfocalQuadric;

    let ell = ConfocalQuadric {
        a_param: a,
        b_param: b,
        lambda: lambda_ell,
    };
    let hyp_r = ConfocalQuadric {
        a_param: a,
        b_param: b,
        lambda: lambda_hyp_right,
    };

    let [tr, br, bl, tl] = ConfocalQuadric::intersections(&ell, &hyp_r)
        .expect("ellipse and right hyperbola must intersect");

    Domain::new(vec![
        Segment::Quad {
            curve: ell,
            a: tr,
            b: tl,
        },
        Segment::Quad {
            curve: hyp_r,
            a: tl,
            b: bl,
        },
        Segment::Quad {
            curve: ell,
            a: bl,
            b: br,
        },
        Segment::Quad {
            curve: hyp_r,
            a: br,
            b: tr,
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
