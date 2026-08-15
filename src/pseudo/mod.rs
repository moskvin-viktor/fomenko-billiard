//! Pseudo-integrable confocal tables: flat-coordinate mapping for tables with
//! reflex (3π/2) corners.
//!
//! A table with `n` reflex corners has level sets of genus `1 + n` (per
//! connected component).  The correct coordinates are the un-normalized length
//! coordinates `uᵢ = ∫ dλᵢ/√P`, in which the flow has slope exactly ±1 and the
//! level set becomes a flat rectilinear billiard (a translation surface).
//!
//! - [`ulength`] — the `u(λ) = ∫ dλ/√P` quadrature table.
//! - [`classify`] — level classification (forbidden / torus / genus-≥2).
//! - [`map`] — per-sample mapping into flat coordinates.

pub mod classify;
pub mod map;
pub mod ulength;

pub use classify::{classify_level, AccessibleRegion, Level};
pub use map::{reflex_corners_in_flat, to_flat, FlatError, FlatSample};
pub use ulength::{RootSide, ULength};
