use crate::domain;
use crate::domain::Domain;
use macroquad::prelude::*;

pub struct Preset {
    pub label: &'static str,
    pub domain: Domain,
    /// Label for the second integral slider: "Λ" for confocal, "θ/π" for polyline.
    pub second_integral_label: &'static str,
    /// Whether this is a confocal (quadric-based) billiard where we use caustic start points.
    pub is_confocal: bool,
    /// For polyline billiards, a fixed interior starting point.
    /// For confocal billiards, unused.
    pub start_center: Vec2,
}

/// Confocal family: (b − λ)x² + (a − λ)y² = (a − λ)(b − λ),  λ ≤ a.
///
///   λ <  b   →  ellipse
///   λ =  b   →  degenerate (segment between foci + horizontal rays)
///   b < λ < a →  hyperbola (opens left/right)
///   λ =  a   →  vertical segment
const A: f32 = 4.0;
const B: f32 = 1.0;

pub fn all_presets() -> Vec<Preset> {
    vec![
        Preset {
            label: "Square: ellipse λ=0, hyperbola λ=2.5",
            domain: domain::confocal_quad(A, B, 0.0, 2.5),
            second_integral_label: "Λ",
            is_confocal: true,
            start_center: vec2(0.0, 0.0),
        },
        Preset {
            label: "Thin: ellipse λ=-1, hyperbola λ=2.8",
            domain: domain::confocal_quad(A, B, -1.0, 2.8),
            second_integral_label: "Λ",
            is_confocal: true,
            start_center: vec2(0.0, 0.0),
        },
        Preset {
            label: "Flat: ellipse λ=0.5, hyperbola λ=2.2",
            domain: domain::confocal_quad(A, B, 0.5, 2.2),
            second_integral_label: "Λ",
            is_confocal: true,
            start_center: vec2(0.0, 0.0),
        },
        Preset {
            label: "L-shape (confocal, 3π/2 corner)",
            domain: domain::confocal_lshape(A, B, 0.0, 2.5, 3.5, 2.8, 0.7),
            second_integral_label: "Λ",
            is_confocal: true,
            start_center: vec2(0.0, 0.0),
        },
        Preset {
            label: "Square (polyline)",
            domain: domain::square(),
            second_integral_label: "θ/π",
            is_confocal: false,
            start_center: vec2(0.0, 0.0),
        },
        Preset {
            label: "L-shape (polyline)",
            domain: domain::lshape_poly(),
            second_integral_label: "θ/π",
            is_confocal: false,
            start_center: vec2(0.6, 0.6),
        },
    ]
}
