//! UI widgets for the app.
//!
//! Extracted from the old `main.rs` monolith so the binary is a thin shell and
//! the drag-widget logic lives in a module of its own.

use crate::render;
use crate::second_integral::SecondIntegralRange;
use macroquad::prelude::*;

/// Slider value change below which we ignore (avoid rebuilds on fp noise).
pub const SLIDER_JITTER: f32 = 0.0001;

/// A horizontal drag slider over a [`SecondIntegralRange`].
pub struct Slider {
    pub x: f32,
    pub y: f32,
    pub width: f32,
    dragging: bool,
}

impl Slider {
    pub fn new() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            width: 0.0,
            dragging: false,
        }
    }

    pub fn update(&mut self, value: f32, range: &SecondIntegralRange, label: &str) -> Option<f32> {
        let track_h = 4.0;
        let thumb_r = 8.0;
        let cy = self.y + track_h / 2.0;

        // Track background
        draw_rectangle(
            self.x,
            self.y,
            self.width,
            track_h,
            color_u8!(60, 60, 100, 180),
        );

        // Filled portion
        let t = range.fraction_of_value(value);
        let fill_w = t * self.width;
        if fill_w > 0.0 {
            draw_rectangle(
                self.x,
                self.y,
                fill_w,
                track_h,
                color_u8!(130, 130, 200, 220),
            );
        }

        // Thumb
        let thumb_x = self.x + t * self.width;
        draw_circle(thumb_x, cy, thumb_r, color_u8!(220, 220, 255, 255));
        draw_circle_lines(thumb_x, cy, thumb_r, 1.5, color_u8!(100, 100, 160, 200));

        // Gap indicator (confocal-only)
        if let SecondIntegralRange::Confocal {
            ell_max, hyp_min, ..
        } = range
        {
            if ell_max < hyp_min {
                let gap_t = range.fraction_of_value(*ell_max);
                let gap_x = self.x + gap_t * self.width;
                draw_line(
                    gap_x,
                    self.y - 2.0,
                    gap_x,
                    self.y + track_h + 2.0,
                    2.0,
                    render::BG,
                );
            }
        }

        // Label and value text
        let value_text = match range {
            SecondIntegralRange::Angle => format!("{:.3}", value),
            SecondIntegralRange::Confocal { .. } => {
                let cf = crate::torus::ConfocalParams::standard();
                let caustic_label = if value < cf.b { "ellipse" } else { "hyperbola" };
                format!("{:.3} ({})", value, caustic_label)
            }
        };
        let info = format!("{} = {}", label, value_text);
        draw_text(
            &info,
            self.x - 60.0,
            self.y + track_h + 22.0,
            14.0,
            color_u8!(180, 180, 210, 200),
        );

        // Mouse interaction
        let mx = mouse_position().0;
        let my = mouse_position().1;

        if is_mouse_button_pressed(MouseButton::Left) {
            let dist = ((mx - thumb_x).powi(2) + (my - cy).powi(2)).sqrt();
            if dist <= thumb_r + 4.0
                || (mx >= self.x && mx <= self.x + self.width && (my - cy).abs() <= 12.0)
            {
                self.dragging = true;
            }
        }

        if self.dragging {
            if is_mouse_button_down(MouseButton::Left) {
                let t = ((mx - self.x) / self.width).clamp(0.0, 1.0);
                return Some(range.value_at_fraction(t));
            } else {
                self.dragging = false;
            }
        }

        None
    }
}
