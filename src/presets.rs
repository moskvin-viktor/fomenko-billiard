use crate::domain;
use crate::domain::Domain;
use macroquad::prelude::*;

pub struct Preset {
    pub label: &'static str,
    pub domain: Domain,
    pub start: Vec2,
    pub vel: Vec2,
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
            start: vec2(0.0, 0.0),
            vel: vec2(0.4, 0.3),
        },
        Preset {
            label: "Thin: ellipse λ=-1, hyperbola λ=2.8",
            domain: domain::confocal_quad(A, B, -1.0, 2.8),
            start: vec2(0.0, 0.0),
            vel: vec2(0.3, 0.5),
        },
        Preset {
            label: "Flat: ellipse λ=0.5, hyperbola λ=2.2",
            domain: domain::confocal_quad(A, B, 0.5, 2.2),
            start: vec2(0.0, 0.0),
            vel: vec2(0.5, 0.2),
        },
        Preset {
            label: "L-shape (polyline)",
            domain: domain::lshape_poly(),
            start: vec2(0.3, 0.3),
            vel: vec2(0.5, 0.7),
        },
        Preset {
            label: "Square (polyline)",
            domain: domain::square(),
            start: vec2(0.0, 0.0),
            vel: vec2(0.4, 0.6),
        },
    ]
}
