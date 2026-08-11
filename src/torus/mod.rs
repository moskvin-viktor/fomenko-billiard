//! Liouville tori for confocal billiards, in the Jacobi normalization.
//!
//! Maps a phase-space sample `(x, y, vx, vy)` to a point `(θ₁, θ₂, torus_index)`
//! on a Liouville torus, via the phase `w(λ) = ∫ dλ/√P` on ovals of the cubic
//! `P(λ) = (a − λ)(b − λ)(λc − λ)`.
//!
//! - [`confocal`] — confocal coordinates and phase-space helpers.
//! - [`quadrature`] — tabulating and interpolating the Abelian phase.
//! - [`map`] — the top-level phase-space → torus mapping.

pub mod confocal;
pub mod map;
pub mod quadrature;

// Re-exports so the public surface is unchanged: callers use `torus::xyz`
// rather than reaching into the submodules.
pub use confocal::{caustic, confocal, lam_dots, nu_angle, ConfocalParams, PhaseSample};
pub use map::{to_torus, TorusCache, TorusParams};
pub use quadrature::Libration;
