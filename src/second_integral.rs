//! The valid range of the second integral of motion, and its normalization.
//!
//! This is the "interface" the UI consumes: given a [`Domain`], what values of
//! the second integral are realizable, and how to interpolate a slider position
//! `t ∈ [0, 1]` over that range and back.  All the segmentation into the
//! ellipse / hyperbola sides and the degraded-boundary margins lives here, so
//! the binary only calls [`SecondIntegralRange::value_at_fraction`] and
//! [`SecondIntegralRange::fraction_of_value`].

use crate::{domain, torus::ConfocalParams};

/// Inward margin used to stay clear of the degenerate boundaries at the ends
/// of the second-integral range (Λ = B separatrix, Λ close to A, and the
/// hyperbola walls).  Used by clamping and animation.
const SECOND_INT_EPS: f32 = 0.05;

/// Tighter margin used only by the slider, so the thumb can get closer to the
/// separatrix than the clamp/anim will allow.
const SLIDER_EPS: f32 = 0.02;

/// The valid range of the second integral for a domain, plus the normalized
/// [`0,1]`↔value mapping a slider needs.
///
/// For confocal domains this is the caustic parameter Λ, segmented into an
/// elliptic side `(λ_ell, B)` and a hyperbolic side `(λ_hyp, A)` separated by
/// the forbidden gap at the separatrix Λ = B.  For polyline domains it is a
/// single normalized angle θ/π ∈ [-1, 1].
pub enum SecondIntegralRange {
    /// Confocal caustic parameter Λ, with the valid sub-ranges precomputed.
    Confocal {
        ell_min: f32,
        ell_max: f32,
        hyp_min: f32,
        hyp_max: f32,
        /// Combined `ell_range + hyp_range`, used to map `t ∈ [0,1]`.
        total_range: f32,
    },
    /// Polyline angle θ/π ∈ [-1, 1].
    Angle,
}

impl SecondIntegralRange {
    /// Build the range for a domain.  `is_confocal` selects the confocal Λ
    /// range versus the polyline angle.
    pub fn for_domain(domain: &domain::Domain, is_confocal: bool) -> Self {
        if !is_confocal {
            return Self::Angle;
        }
        let (ell_min, ell_max, hyp_min, hyp_max) = Self::confocal_bounds(domain);
        let ell_range = (ell_max - ell_min).max(0.0);
        let hyp_range = (hyp_max - hyp_min).max(0.0);
        let total_range = ell_range + hyp_range;
        Self::Confocal {
            ell_min,
            ell_max,
            hyp_min,
            hyp_max,
            total_range,
        }
    }

    /// The four usable Λ bounds `(ell_min, ell_max, hyp_min, hyp_max)`.
    ///
    /// Boundary lambdas come from the lib's shared segment-walk; the margins
    /// (`SLIDER_EPS`) leave room for the slider thumb near the separatrix.
    fn confocal_bounds(domain: &domain::Domain) -> (f32, f32, f32, f32) {
        let cf = ConfocalParams::standard();
        let (lambda_ell, lambda_hyp) = match crate::confocal::ConfocalStructure::of_domain(domain) {
            Some(s) => (s.lambda_ell, s.lambda_hyp),
            // No quadric arcs (polyline): no confocal constraint → full range.
            None => (0.0, Some(cf.b + 1.0)),
        };
        let e = SLIDER_EPS;
        // For a full ellipse (no hyperbola wall) the hyperbolic side starts just
        // above the focal separatrix `b`, not at an arbitrary `b + 1`.
        let hyp_start = lambda_hyp.unwrap_or(cf.b);
        (
            lambda_ell + e, // ell_min
            cf.b - e,       // ell_max
            hyp_start + e,  // hyp_min
            cf.a - e,       // hyp_max
        )
    }

    /// The value for a slider fraction `t ∈ [0, 1]` of the track.
    pub fn value_at_fraction(&self, t: f32) -> f32 {
        match self {
            Self::Angle => -1.0 + 2.0 * t.clamp(0.0, 1.0),
            Self::Confocal {
                ell_min,
                ell_max,
                hyp_min,
                total_range,
                ..
            } => {
                if *total_range <= 0.0 {
                    return 0.0;
                }
                let ell_range = ell_max - ell_min;
                let t = t.clamp(0.0, 1.0);
                let pos = t * total_range;
                if pos <= ell_range {
                    ell_min + pos
                } else {
                    hyp_min + (pos - ell_range)
                }
            }
        }
    }

    /// The track fraction `∈ [0, 1]` for a value.
    pub fn fraction_of_value(&self, lam: f32) -> f32 {
        match self {
            Self::Angle => (lam.clamp(-1.0, 1.0) + 1.0) / 2.0,
            Self::Confocal {
                ell_min,
                ell_max,
                hyp_min,
                total_range,
                ..
            } => {
                if *total_range <= 0.0 {
                    return 0.0;
                }
                let ell_range = ell_max - ell_min;
                if lam <= *ell_max {
                    (lam - ell_min) / total_range
                } else if lam >= *hyp_min {
                    (ell_range + (lam - hyp_min)) / total_range
                } else {
                    ell_range / total_range
                }
            }
        }
    }
}

/// The valid Λ range for a confocal domain, as a small value object.
///
/// Exposes the boundary lambdas so the app's clamping / animation logic can ask
/// "what values of Λ are reachable" without re-walking the domain.  This is the
/// same interface the slider uses.
pub struct LambdaRange {
    /// λ_ell — the outer ellipse boundary (min `λ < B`).
    pub lambda_ell: f32,
    /// λ_hyp — the inner hyperbola boundary (max `λ > B`).
    pub lambda_hyp: f32,
}

impl LambdaRange {
    /// Extract `(λ_ell, λ_hyp)` from the domain, or the full-range fallbacks
    /// (`0`, `B + 1`) when the domain has no confocal arcs.
    pub fn of_domain(domain: &domain::Domain) -> Self {
        let cf = ConfocalParams::standard();
        match crate::confocal::ConfocalStructure::of_domain(domain) {
            Some(s) => Self {
                lambda_ell: s.lambda_ell,
                // Full ellipse (no hyperbola wall): the hyperbolic side starts
                // just above the focal separatrix `b`.
                lambda_hyp: s.lambda_hyp.unwrap_or(cf.b),
            },
            None => Self {
                lambda_ell: 0.0,
                lambda_hyp: cf.b + 1.0,
            },
        }
    }

    /// The four usable Λ bounds `(ell_min, ell_max, hyp_min, hyp_max)` under
    /// the animation margin.  Used to bounce the animation off each side.
    pub fn animation_bounds(&self) -> (f32, f32, f32, f32) {
        let cf = ConfocalParams::standard();
        let e = SECOND_INT_EPS;
        (
            self.lambda_ell + e, // ell_min
            cf.b - e,            // ell_max
            self.lambda_hyp + e, // hyp_min
            cf.a - e,            // hyp_max
        )
    }

    /// Clamp `lam` to the valid range, jumping over the forbidden gap between
    /// the ellipse and hyperbola sides (the separatrix near Λ = B).
    pub fn clamp(&self, lam: f32, prev: f32) -> f32 {
        let cf = ConfocalParams::standard();
        let (ell_min, ell_max) = (self.lambda_ell + SECOND_INT_EPS, cf.b - SECOND_INT_EPS);
        let (hyp_min, hyp_max) = (self.lambda_hyp + SECOND_INT_EPS, cf.a - SECOND_INT_EPS);

        let clamped = lam.clamp(ell_min, hyp_max);
        if clamped > ell_max && clamped < hyp_min {
            match lam.partial_cmp(&prev).unwrap_or(std::cmp::Ordering::Equal) {
                std::cmp::Ordering::Greater => hyp_min,
                _ => ell_max,
            }
        } else {
            clamped
        }
    }
}
