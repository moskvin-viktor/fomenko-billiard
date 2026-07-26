use macroquad::prelude::*;

/// A single smooth piece of the boundary.
#[derive(Clone, Debug)]
pub enum Segment {
    /// Straight wall from `a` to `b`.
    Line { a: Vec2, b: Vec2 },
    /// Circular arc with given centre, radius, and angular range (CCW, radians).
    Arc {
        centre: Vec2,
        radius: f32,
        start_angle: f32,
        end_angle: f32,
    },
}

impl Segment {
    /// Normal vector at a point on the segment, pointing inward.
    pub fn inward_normal(&self, p: Vec2) -> Vec2 {
        match self {
            Segment::Line { a, b } => {
                let d = *b - *a;
                // Rotate the tangent 90° counter-clockwise → inward for a
                // boundary traversed counter-clockwise.  For the stadium
                // segments are ordered CCW, so left-of-tangent is inward.
                vec2(-d.y, d.x).normalize()
            }
            Segment::Arc { centre, .. } => (*centre - p).normalize(),
        }
    }

    /// Reflect direction `dir` off this segment at point `p`.
    pub fn reflect(&self, p: Vec2, dir: Vec2) -> Vec2 {
        let n = self.inward_normal(p);
        dir - 2.0 * dir.dot(n) * n
    }

    /// Find the smallest `t > 0` such that the ray `p + t·dir` hits this segment,
    /// together with the hit point and the line parameter `s` (normalised along
    /// the segment, 0..1 for lines, angular fraction 0..1 for arcs).
    pub fn intersect(&self, p: Vec2, dir: Vec2) -> Option<(f32, Vec2, f32)> {
        match self {
            Segment::Line { a, b } => {
                let ab = *b - *a;
                let denom = dir.x * ab.y - dir.y * ab.x;
                if denom.abs() < 1e-12 {
                    return None; // parallel
                }
                let ap = *a - p;
                let t = (ap.x * ab.y - ap.y * ab.x) / denom;
                let s = (ap.x * dir.y - ap.y * dir.x) / denom;
                // s ∈ [0, 1] means we hit the segment, t > epsilon means forward
                if s >= 0.0 && s <= 1.0 && t > 1e-8 {
                    Some((t, p + dir * t, s))
                } else {
                    None
                }
            }
            Segment::Arc {
                centre,
                radius,
                start_angle,
                end_angle,
            } => {
                // Solve |p + t·dir - centre|² = radius²
                let oc = p - *centre;
                let qa = dir.dot(dir);
                let qb = 2.0 * oc.dot(dir);
                let qc = oc.dot(oc) - radius * radius;

                let disc = qb * qb - 4.0 * qa * qc;
                if disc < 0.0 {
                    return None;
                }
                let sd = disc.sqrt();
                let t1 = (-qb - sd) / (2.0 * qa);
                let t2 = (-qb + sd) / (2.0 * qa);

                for t in [t1, t2] {
                    if t <= 1e-8 {
                        continue;
                    }
                    let hit = p + dir * t;
                    let mut angle = (hit - *centre).y.atan2((hit - *centre).x);
                    // Normalise to [0, 2π)
                    if angle < 0.0 {
                        angle += 2.0 * std::f32::consts::PI;
                    }
                    let sa = if *start_angle < 0.0 {
                        *start_angle + 2.0 * std::f32::consts::PI
                    } else {
                        *start_angle
                    };
                    let ea = if *end_angle < 0.0 {
                        *end_angle + 2.0 * std::f32::consts::PI
                    } else {
                        *end_angle
                    };
                    // The arc runs from sa to ea in the CCW sense.
                    // For a CCW arc the shorter way, fold angle into [sa, ea].
                    // We check by seeing if angle is between sa and ea, handling wrap.
                    let delta = if ea > sa {
                        ea - sa
                    } else {
                        ea + 2.0 * std::f32::consts::PI - sa
                    };
                    let a_norm = if angle >= sa {
                        angle - sa
                    } else {
                        angle + 2.0 * std::f32::consts::PI - sa
                    };
                    if a_norm <= delta + 1e-6 {
                        let frac = a_norm / delta;
                        return Some((t, hit, frac));
                    }
                }
                None
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Domain – a closed piecewise-smooth boundary with corners.
//
// The boundary is traversed counter-clockwise. Consecutive segments meet
// at corners.  Interior angles at corners are π/2 or 3π/2.
// ---------------------------------------------------------------------------
#[derive(Clone, Debug)]
pub struct Domain {
    pub segments: Vec<Segment>,
}

impl Domain {
    pub fn new(segments: Vec<Segment>) -> Self {
        Self { segments }
    }

    /// The corner points between consecutive segments.
    pub fn corners(&self) -> Vec<Vec2> {
        self.segments
            .windows(2)
            .map(|w| match (&w[0], &w[1]) {
                (Segment::Line { b, .. }, _) => *b,
                (
                    Segment::Arc {
                        centre,
                        radius,
                        end_angle,
                        ..
                    },
                    _,
                ) => *centre + *radius * vec2(end_angle.cos(), end_angle.sin()),
            })
            .collect()
    }

    /// Check if a point is inside the domain (strictly).
    /// Uses the winding number / ray-casting approach: count crossings.
    pub fn contains(&self, point: Vec2) -> bool {
        // Shoot a ray to the right and count intersections with boundary segments
        let ray_dir = vec2(1.0, 0.0);
        let mut crossings = 0;
        for seg in &self.segments {
            if let Some((t, _hit, _)) = seg.intersect(point, ray_dir) {
                // Ignore hits at the exact start (t ≈ 0 isn't possible since we're inside)
                if t > 1e-8 {
                    crossings += 1;
                }
            }
        }
        crossings % 2 == 1
    }

    /// Find the nearest intersection with the boundary along the ray `p + t·dir`,
    /// returning `(t, segment_index, hit_point)`.
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

    /// Reflect direction off the segment at `segment_idx` at point `p`.
    /// Special handling for corners: if the hit is within ε of a corner
    /// with interior angle π/2, reflect back along the arrival direction
    /// (continuous extension).
    pub fn reflect(&self, p: Vec2, dir: Vec2, segment_idx: usize) -> Vec2 {
        // Check if we're near a corner
        let corners: Vec<(Vec2, f32)> = self
            .segments
            .windows(2)
            .map(|w| {
                let corner = match &w[0] {
                    Segment::Line { b, .. } => *b,
                    Segment::Arc {
                        centre,
                        radius,
                        end_angle,
                        ..
                    } => *centre + *radius * vec2(end_angle.cos(), end_angle.sin()),
                };
                // Estimate interior angle between the two segments at this corner
                let n1 = w[0].inward_normal(corner);
                let n2 = w[1].inward_normal(corner);
                let angle = n1.angle_between(n2).abs();
                (corner, angle)
            })
            .collect();

        let eps = 1e-4;
        for (corner, angle) in &corners {
            if p.distance(*corner) < eps {
                // π/2 corner: reflect back along arrival direction
                if (*angle - std::f32::consts::PI / 2.0).abs() < 0.1 {
                    return -dir;
                }
                // 3π/2 corner: standard reflection off the first segment is fine
                break;
            }
        }

        self.segments[segment_idx].reflect(p, dir)
    }

    /// Trace trajectory from start point `p` with initial velocity `v`.
    /// Returns segments as `(from, to)` world-space pairs.
    pub fn trace(&self, mut p: Vec2, mut v: Vec2, max_steps: usize) -> Vec<(Vec2, Vec2)> {
        let mut segs = Vec::with_capacity(max_steps);

        for _ in 0..max_steps {
            let speed = v.length();
            if speed < 1e-12 {
                break;
            }
            let dir = v / speed;

            let (_t, idx, hit) = match self.intersect(p, dir) {
                Some(result) => result,
                None => break,
            };

            segs.push((p, hit));

            v = self.reflect(hit, v, idx);
            // Nudge inward along the reflected direction so we don't
            // re-detect the same hit on the next iteration.
            let new_dir = v.normalize();
            p = hit + 1e-4 * new_dir;
        }

        segs
    }

    /// Sample boundary points for drawing (~n per segment).
    pub fn sample_boundary(&self, n: usize) -> Vec<Vec2> {
        let mut pts = Vec::new();
        for seg in &self.segments {
            match seg {
                Segment::Line { a, b } => {
                    pts.push(*a);
                    pts.push(*b);
                }
                Segment::Arc {
                    centre,
                    radius,
                    start_angle,
                    end_angle,
                } => {
                    let mut delta = end_angle - start_angle;
                    if delta <= 0.0 {
                        delta += 2.0 * std::f32::consts::PI;
                    }
                    let steps = (n as f32 * delta / (2.0 * std::f32::consts::PI)).max(3.0) as usize;
                    for i in 0..=steps {
                        let frac = i as f32 / steps as f32;
                        let angle = start_angle + frac * delta;
                        pts.push(*centre + *radius * vec2(angle.cos(), angle.sin()));
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

/// Bunimovich stadium: rectangle [-1, 1] × [-1, 1] capped with semicircular
/// arcs of radius 1 on the left and right.  Corners at π/2.
pub fn stadium() -> Domain {
    use std::f32::consts::FRAC_PI_2;

    let r = 1.0;
    Domain::new(vec![
        // Top straight: right → left
        Segment::Line {
            a: vec2(r, r),
            b: vec2(-r, r),
        },
        // Left semicircular arc (CCW from top to bottom)
        Segment::Arc {
            centre: vec2(-r, 0.0),
            radius: r,
            start_angle: FRAC_PI_2,
            end_angle: -FRAC_PI_2,
        },
        // Bottom straight: left → right
        Segment::Line {
            a: vec2(-r, -r),
            b: vec2(r, -r),
        },
        // Right semicircular arc (CCW from bottom to top)
        Segment::Arc {
            centre: vec2(r, 0.0),
            radius: r,
            start_angle: -FRAC_PI_2,
            end_angle: FRAC_PI_2,
        },
    ])
}

/// Rectangle [0, w] × [0, h].
pub fn rectangle(w: f32, h: f32) -> Domain {
    Domain::new(vec![
        Segment::Line {
            a: vec2(0.0, 0.0),
            b: vec2(w, 0.0),
        },
        Segment::Line {
            a: vec2(w, 0.0),
            b: vec2(w, h),
        },
        Segment::Line {
            a: vec2(w, h),
            b: vec2(0.0, h),
        },
        Segment::Line {
            a: vec2(0.0, h),
            b: vec2(0.0, 0.0),
        },
    ])
}

/// L-shaped domain: a 2×2 square with a 1×1 cutout in the top-right.
/// Corners at π/2 and 3π/2.
pub fn l_shape() -> Domain {
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
