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
//! - [`render`] — 3D embeddings (torus angles, cross, double-torus pretzel).
//! - [`view`] — flat-chart sampling and drawing for the 3D view.
//!
//! See `docs/pseudo_integrable.md` for the end-to-end pipeline.

pub mod classify;
pub mod geometry;
pub mod map;
pub mod morph;
pub mod render;
pub mod ulength;
pub mod view;

pub use classify::{classify_level, table_touches_focal, AccessibleRegion, Level};
pub use geometry::{level_geometry, LevelGeometry};
pub use map::{reflex_corners_in_flat, to_flat, FlatError, FlatSample};
pub use morph::{flat_morph, pinch_point_3d, unified_embed, unified_with_normal, FlatMorph};
pub use render::{
    cross_embed, donut_with_normal, pretzel_embed, pretzel_with_normal, sheets, torus_angles,
    CrossModuli,
};
pub use ulength::{RootSide, ULength};
pub use view::{
    draw_flat, draw_flat_boundary, draw_flat_highlights, sample_flat_boundary,
    sample_flat_trajectory, FlatPhasePoint,
};
