//! 2D molecule strip: the Reeb graph over the λ-axis rendered as markers, and
//! the "all special layers" navigation.
//!
//! The strip shows the caustic axis `λ` horizontally, with a marker at every
//! special layer (walls, `b` separatrix, molecule critical values, `a` focal
//! axis), coloured by its kind, and a cursor at the current `Λ`.  This makes
//! every special layer visible and the current caustic's place in the molecule
//! explicit, and the navigation keys snap `Λ` to each special layer in turn.

use macroquad::prelude::*;

use crate::molecule::TransitionKind;
use crate::presets::Preset;
use crate::torus::ConfocalParams;

/// A marker on the strip.
pub struct LayerMarker {
    /// The caustic value.
    pub lam: f32,
    /// Human label.
    pub label: &'static str,
    /// Colour key (by kind).
    pub kind: TransitionKind,
}

/// The colour of a transition kind on the strip.
pub fn kind_color(kind: TransitionKind) -> Color {
    match kind {
        TransitionKind::BirthA | TransitionKind::DeathA | TransitionKind::AEnd => {
            color_u8!(120, 220, 120, 255) // green: atoms
        }
        TransitionKind::GenusJump => color_u8!(255, 160, 80, 255), // orange
        TransitionKind::SplitMerge => color_u8!(120, 160, 255, 255), // blue
        TransitionKind::EdgeSwap => color_u8!(160, 160, 180, 255), // grey
    }
}

/// Build the ordered layer markers for a preset.
///
/// Includes every special layer (walls, `b` separatrix, `a` focal axis,
/// molecule critical values) *and* uniform regular layers between them, so
/// stepping through the strip visits the whole sweep and can't skip a layer.
pub fn layer_markers(preset: &Preset) -> Vec<LayerMarker> {
    let cf = ConfocalParams::standard();
    let domain = &preset.domain;
    let layers = crate::molecule::all_layers(domain, &cf);

    let mut marks = Vec::with_capacity(layers.len());
    for lam in layers {
        let label = crate::molecule::layer_label(lam, domain, &cf);
        let kind = classify_kind(lam, domain, &cf);
        marks.push(LayerMarker { lam, label, kind });
    }
    marks
}

/// Classify a special value to a [`TransitionKind`] for colouring.  Degeneracies
/// and walls map to a representative atom; molecule critical values map to their
/// transition; ordinary (regular) reachable levels between criticals get an
/// edge-swap so they render neutrally rather than as a singularity.
fn classify_kind(lam: f32, domain: &crate::domain::Domain, cf: &ConfocalParams) -> TransitionKind {
    if (lam - cf.b).abs() < 1e-4 {
        return TransitionKind::SplitMerge;
    }
    if (lam - cf.a).abs() < 1e-4 {
        return TransitionKind::AEnd;
    }
    if let Some(tab) = crate::table::Table::from_domain(domain, cf) {
        // Only a molecule critical value is a real transition; other levels are
        // regular and should render neutrally.
        let crit = crate::molecule::critical_values(&tab, cf, 1e-6);
        if crit.iter().any(|&c| (c - lam).abs() < 1e-4) {
            return crate::molecule::classify_transition(&tab, cf, lam, 1e-4).kind;
        }
    }
    TransitionKind::EdgeSwap
}

/// Draw the strip across the view and label the current layer.  Returns the
/// number of special layers drawn.
pub fn draw_strip(markers: &[LayerMarker], current_lam: f32, x: f32, y: f32, width: f32) -> usize {
    if markers.is_empty() {
        return 0;
    }
    let lo = markers.first().unwrap().lam;
    let hi = markers.last().unwrap().lam;
    let span = (hi - lo).max(1e-3);

    // Axis line.
    draw_line(x, y, x + width, y, 2.0, color_u8!(120, 120, 150, 160));

    // Markers (all of them); only label those sufficiently spaced to avoid
    // overlapping text.  Alternate label rows to reduce collisions.
    let min_gap = 60.0; // px between labelled markers
    let mut prev_labelled_x = f32::MIN;
    for (i, m) in markers.iter().enumerate() {
        let t = ((m.lam - lo) / span).clamp(0.0, 1.0);
        let sx = x + t * width;
        draw_circle(sx, y, 5.0, kind_color(m.kind));
        if sx - prev_labelled_x >= min_gap {
            let row = (i % 2) as f32 * 14.0;
            draw_text(m.label, sx + 7.0, y - 6.0 - row, 11.0, kind_color(m.kind));
            prev_labelled_x = sx;
        }
    }

    // Current-position cursor + label of the nearest layer.
    let t = ((current_lam - lo) / span).clamp(0.0, 1.0);
    let sx = x + t * width;
    draw_line(sx, y - 18.0, sx, y + 18.0, 2.0, WHITE);
    draw_circle(sx, y, 6.0, WHITE);
    if let Some(m) = nearest(markers, current_lam) {
        draw_text(
            &format!("Λ={:.3}  {}", current_lam, m.label),
            sx + 8.0,
            y + 30.0,
            13.0,
            color_u8!(230, 230, 245, 255),
        );
    }

    markers.len()
}

/// The marker nearest to `current`.
pub fn nearest(markers: &[LayerMarker], current: f32) -> Option<&LayerMarker> {
    markers
        .iter()
        .min_by(|a, b| (a.lam - current).abs().total_cmp(&(b.lam - current).abs()))
}

/// The previous / next special layer value around `current`, wrapping.
pub fn snap_nearest(markers: &[LayerMarker], current: f32, dir: i32) -> Option<f32> {
    if markers.is_empty() {
        return None;
    }
    let n = markers.len() as i32;
    let mut best = 0;
    let mut best_d = f32::MAX;
    for (i, m) in markers.iter().enumerate() {
        let d = (m.lam - current).abs();
        if d < best_d {
            best_d = d;
            best = i as i32;
        }
    }
    let next = (best + dir + n) % n;
    Some(markers[next as usize].lam)
}
