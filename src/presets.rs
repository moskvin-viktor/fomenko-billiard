use crate::domain::{self, Domain};
use macroquad::prelude::*;

pub struct Preset {
    pub label: &'static str,
    pub domain: Domain,
    pub start: Vec2,
    pub vel: Vec2,
}

pub fn all_presets() -> Vec<Preset> {
    let stadium = domain::stadium();
    let rect = domain::rectangle(3.0, 2.0);
    let lshape = domain::l_shape();

    vec![
        Preset {
            label: "Bunimovich stadium",
            domain: stadium,
            start: vec2(-0.3, 0.2),
            vel: vec2(0.8, 0.4),
        },
        Preset {
            label: "Rectangle 3×2",
            domain: rect,
            start: vec2(0.5, 0.5),
            vel: vec2(0.6, 0.9),
        },
        Preset {
            label: "L-shape (π/2 & 3π/2 corners)",
            domain: lshape,
            start: vec2(0.3, 0.3),
            vel: vec2(0.5, 0.7),
        },
    ]
}

