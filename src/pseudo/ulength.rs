//! Un-normalized length coordinates `u(λ) = ∫ dλ/√P` for the pseudo-integrable
//! flat model.
//!
//! This is the doc's `ULength` (§11 of `confocal_L_pseudo_integrable.md`): the
//! `w`-integral of the torus doc *without* the `π/W` normalization, integrated
//! on a grid uniform in `s = √(r − λ)` where `r` is the root of `P` nearest the
//! interval.  That kills the inverse-square-root endpoint singularity when the
//! endpoint *is* a root (turning line), and is harmless when it isn't (wall).
//!
//! Unlike the torus `Libration`, this handles the **both-walls** case (no root
//! inside the interval, e.g. `[β₁, β₃]` with `r = a`) and the **root-at-lo**
//! case (hyperbolic turning at the lower end), which the L-table needs.

/// Which end of the interval the nearest root of `P` sits at.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RootSide {
    /// The root is at/above the upper end (`hi`): `λ = r − s²`.
    Hi,
    /// The root is at/below the lower end (`lo`): `λ = r + s²`.
    Lo,
}

/// `u(λ) = ∫_lo^λ dλ'/√P` on `[lo, hi]`, with `P = (a−λ)(b−λ)(lc−λ)`.
#[derive(Clone, Debug)]
pub struct ULength {
    lo: f32,
    hi: f32,
    /// The root of `P` used for the `s`-substitution.
    r: f32,
    /// Whether `u` is oriented from the `r`-side end (so `u(lo)=0`, `u(hi)=total`).
    flip: bool,
    /// The `s` grid (increasing).
    xs: Vec<f32>,
    /// Cumulative integral values on the grid, oriented from the `r`-side end.
    gs: Vec<f32>,
    /// The full integral `u(hi) = ∫_lo^hi dλ/√P`.
    pub total: f32,
}

impl ULength {
    /// Build the table for `[lo, hi]` in the family `(a, b, lc)`.
    ///
    /// * `RootSide::Hi` — nearest root is at/above `hi` (turning at the upper
    ///   end, or a wall with the root above).  Uses `λ = r − s²`.
    /// * `RootSide::Lo` — nearest root is at/below `lo` (turning at the lower
    ///   end).  Uses `λ = r + s²`.
    ///
    /// `knots` is the number of quadrature points (higher = smoother, costlier).
    pub fn new(lo: f32, hi: f32, a: f32, b: f32, lc: f32, side: RootSide, knots: usize) -> Self {
        let roots = [a, b, lc];
        // `s` runs from `s_start` (at the end where u=0) to `s_end` (where
        // u=total), so the grid is increasing and `u(lo)=0`, `u(hi)=total`.
        //
        //   RootSide::Hi:  λ = r − s², s = √(r−λ).  s_start at λ=hi, s_end at λ=lo.
        //                  g(s) = ∫_λ^hi (cumulative from the hi end); flip → u=total−g.
        //   RootSide::Lo:  λ = r + s², s = √(λ−r).  s_start at λ=lo, s_end at λ=hi.
        //                  g(s) = ∫_lo^λ (cumulative from the lo end); flip=false.
        let (r, s_start, s_end, flip) = match side {
            RootSide::Hi => {
                let r = roots
                    .iter()
                    .copied()
                    .filter(|&q| q >= hi - 1e-14)
                    .min_by(|x, y| x.total_cmp(y))
                    .unwrap_or(a);
                (r, (r - hi).max(0.0).sqrt(), (r - lo).max(0.0).sqrt(), true)
            }
            RootSide::Lo => {
                let r = roots
                    .iter()
                    .copied()
                    .filter(|&q| q <= lo + 1e-14)
                    .max_by(|x, y| x.total_cmp(y))
                    .unwrap_or(b);
                (r, (lo - r).max(0.0).sqrt(), (hi - r).max(0.0).sqrt(), false)
            }
        };

        // `other(λ) = |∏_{q ≠ r} (q − λ)|`.  With |dλ/ds| = 2s and
        // √|P| = s·√other, the integrand |dλ|/√|P| = 2 ds/√other — the s cancels,
        // so the inverse-square-root singularity at the root is removed.
        let other = |lam: f32| {
            let mut out = 1.0f32;
            for &q in &roots {
                if (q - r).abs() > 1e-14 {
                    out *= (q - lam).abs();
                }
            }
            out
        };
        let f = |lam: f32| 2.0 / other(lam).max(1e-300).sqrt();

        // λ as a function of s.
        let lam_of_s = |s: f32| if flip { r - s * s } else { r + s * s };

        let mut xs = Vec::with_capacity(knots);
        let mut gs = Vec::with_capacity(knots);
        let mut acc = 0.0;
        let mut prev_s = s_start;
        let mut prev_f = f(lam_of_s(s_start));
        gs.push(0.0);
        xs.push(s_start);
        for k in 1..knots {
            let s = s_start + (s_end - s_start) * k as f32 / (knots - 1) as f32;
            let cur_f = f(lam_of_s(s));
            acc += 0.5 * (cur_f + prev_f) * (s - prev_s);
            gs.push(acc);
            xs.push(s);
            prev_s = s;
            prev_f = cur_f;
        }
        let total = acc;

        Self {
            lo,
            hi,
            r,
            flip,
            xs,
            gs,
            total,
        }
    }

    /// `u(λ) = ∫_lo^λ dλ'/√P`, clamped to `[lo, hi]`.
    pub fn u(&self, lam: f32) -> f32 {
        let lam = lam.clamp(self.lo, self.hi);
        let s = (self.r - lam).abs().sqrt();
        let val = interp(s, &self.xs, &self.gs);
        if self.flip {
            self.total - val
        } else {
            val
        }
    }
}

/// Linear interpolation of `(xs, ys)` at `x`, clamped to the grid range.
fn interp(x: f32, xs: &[f32], ys: &[f32]) -> f32 {
    if xs.is_empty() {
        return 0.0;
    }
    let i = match xs.binary_search_by(|&v| v.total_cmp(&x)) {
        Ok(i) => i,
        Err(i) if i == xs.len() => xs.len() - 1,
        Err(i) => i,
    };
    if i == 0 {
        return ys[0];
    }
    let (xa, xb) = (xs[i - 1], xs[i]);
    let t = ((x - xa) / (xb - xa).max(1e-30)).clamp(0.0, 1.0);
    ys[i - 1] + t * (ys[i] - ys[i - 1])
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// u(λ) must be monotone increasing and equal to `total` at `hi`.
    #[test]
    fn test_monotone_and_endpoints() {
        let a = 4.0;
        let b = 1.0;
        let lc = 0.5;
        // Both-walls interval: [β₁, β₃] = [2.0, 3.0], nearest root is a=4.
        let ul = ULength::new(2.0, 3.0, a, b, lc, RootSide::Hi, 512);
        assert!((ul.u(2.0) - 0.0).abs() < 1e-6, "u(lo) should be 0");
        assert!(
            (ul.u(3.0) - ul.total).abs() < 1e-6,
            "u(hi) should equal total"
        );
        let mid = ul.u(2.5);
        assert!(mid > 0.0 && mid < ul.total, "u should be strictly inside");
        assert!(ul.total > 0.0, "total must be positive");
    }

    /// Root-at-lo case (hyperbolic turning at the lower end): [lc, β₃].
    #[test]
    fn test_root_at_lo() {
        let a = 4.0;
        let b = 1.0;
        let lc = 2.2;
        let ul = ULength::new(lc, 3.0, a, b, lc, RootSide::Lo, 512);
        assert!((ul.u(lc) - 0.0).abs() < 1e-6, "u(lo) should be 0");
        assert!(
            (ul.u(3.0) - ul.total).abs() < 1e-6,
            "u(hi) should equal total"
        );
        assert!(ul.total > 0.0);
    }

    /// Turning at the upper end (elliptic): [0, lc].
    #[test]
    fn test_root_at_hi() {
        let a = 4.0;
        let b = 1.0;
        let lc = 0.5;
        let ul = ULength::new(0.0, lc, a, b, lc, RootSide::Hi, 512);
        assert!((ul.u(0.0) - 0.0).abs() < 1e-6);
        assert!((ul.u(lc) - ul.total).abs() < 1e-6);
        assert!(ul.total > 0.0);
    }

    /// The integral must be finite and smooth even when the interval endpoint
    /// is exactly a root (turning line) — the whole point of the substitution.
    #[test]
    fn test_finite_at_turning_line() {
        let a = 4.0;
        let b = 1.0;
        let lc = 0.5;
        let ul = ULength::new(0.0, lc, a, b, lc, RootSide::Hi, 512);
        assert!(
            ul.total.is_finite(),
            "integral must be finite at a turning line"
        );
        assert!(ul.total > 0.0);
    }
}
