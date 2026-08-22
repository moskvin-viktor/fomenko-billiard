//! Quadrature: tabulate the Abelian phase `w(λ) = ∫ dλ/√P` and interpolate it.
//!
//! The cubic `P` is fixed along a trajectory (`λc` constant), so each `Libration`
//! is built once — the expensive `sin`/`sqrt` quadrature — and then every sample
//! is an `O(log n)` lookup.  The two branches are the two endpoint-singularity
//! substitutions described on `Libration`.

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

        // A degenerate span: the interval collapses to a point (lo == hi), or
        // two roots coincide near the separatrix (e.g. lo == lc sitting within
        // 1e-6 of the focal root b), leaving no distinct third root to
        // integrate against.  Return a zero-length libration instead of
        // panicking — callers are expected to skip such a level via the
        // separatrix check in `to_torus`.
        if (hi - lo).abs() < 1e-9 || !has_third_root(roots, lo, hi) {
            return Self {
                lo,
                hi,
                both_roots: false,
                xs: vec![0.0],
                gs: vec![0.0],
                w_full: 0.0,
            };
        }

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

/// True if there is a third root distinct from both `lo` and `hi`.
fn has_third_root(roots: [f32; 3], lo: f32, hi: f32) -> bool {
    roots
        .iter()
        .any(|&r| (r - lo).abs() > 1e-6 && (r - hi).abs() > 1e-6)
}

/// Index of the root closest to `target` (ties → earliest).
fn argmin_abs(roots: [f32; 3], target: f32) -> usize {
    roots
        .iter()
        .enumerate()
        .min_by(|&(_, ra), &(_, rb)| (ra - target).abs().total_cmp(&(rb - target).abs()))
        .map(|(j, _)| j)
        .unwrap_or(0)
}

/// Linear interpolation of `(xs, ys)` at `x`, clamped to the grid range.
///
/// `xs` is strictly increasing, so `binary_search_by` finds either the exact
/// knot (`Ok(i)`) or the insertion point (`Err(i)` → `x` lies in `(xs[i-1],
/// xs[i])`).  Boundary clamping falls out of the same search.
fn interp(x: f32, xs: &[f32], ys: &[f32]) -> f32 {
    if xs.is_empty() {
        return 0.0;
    }

    // `i` is the index of the right endpoint of the interval containing `x`.
    let i = match xs.binary_search_by(|&v| v.total_cmp(&x)) {
        Ok(i) => i,
        Err(i) if i == xs.len() => xs.len() - 1, // x >= xs[n-1]
        Err(i) => i, // x <= xs[0] iff i == 0; else x in (xs[i-1], xs[i])
    };
    if i == 0 {
        return ys[0];
    }

    // Linear segment from xs[i-1] to xs[i], clamped so the upper edge returns
    // ys[n-1] exactly (unchanged from the original clamping behaviour).
    let (xa, xb) = (xs[i - 1], xs[i]);
    let t = ((x - xa) / (xb - xa).max(1e-30)).clamp(0.0, 1.0);
    ys[i - 1] + t * (ys[i] - ys[i - 1])
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A degenerate span (lo == hi, both coinciding with a root, e.g. the
    /// caustic sitting exactly on a wall) must not panic and must yield a
    /// zero-length libration.
    #[test]
    fn degenerate_span_does_not_panic() {
        let lo = 1.0;
        let hi = 1.0;
        let roots = [4.0, 1.0, 1.0]; // lo == hi == a root (b or lc)
        let lib = Libration::new(lo, hi, roots, 512);
        assert_eq!(lib.w_full, 0.0);
        assert!((lib.w(lo) - 0.0).abs() < 1e-9);
    }

    /// A normal both-roots span still integrates to a positive half-period.
    #[test]
    fn both_roots_span_is_positive() {
        let lo = 2.0;
        let hi = 4.0;
        let roots = [6.0, 1.0, 2.0]; // lo is a root (lc=2), hi=4 a wall, third=6
        let lib = Libration::new(lo, hi, roots, 512);
        assert!(lib.w_full > 0.0);
        assert!(lib.w_full.is_finite());
    }

    /// The reported panic: hyperbolic path with lo = lc sitting within 1e-6 of
    /// the focal root b, so the third candidate root (b) coincides with lo and
    /// no distinct third root exists.  Must not panic.
    #[test]
    fn hyperbolic_lo_near_focal_root_does_not_panic() {
        // a = 4, b = 1, lc = 1 + 5e-7 (within 1e-6 of b).  roots = [a, b, lc].
        let a = 4.0f32;
        let b = 1.0f32;
        let lc = b + 5e-7;
        let lo = lc;
        let hi = a;
        let roots = [a, b, lc];
        let lib = Libration::new(lo, hi, roots, 512);
        assert!(lib.w_full.is_finite(), "w_full must be finite, not panic");
    }
}
