// The billiard's 2D drawing (Camera, CachedDomain, draw_*) lives in the
// `render` module, the second-integral range math in the lib
// (`second_integral`), the drag widget in `ui`, and the app loop in `app`.
// This file is a thin shell: construct the app and run it.

use billiards::app::App;
use macroquad::prelude::*;

// The built-in default font (`ProggyClean.ttf`) is ASCII-only, so every π, λ,
// Λ, θ, ² etc in the HUD/labels rendered as a blank tofu box. Bundle a font
// with full Greek + math-symbol coverage and install it as the default so
// every existing `draw_text` call (there's no per-call font override
// anywhere in this codebase) picks it up automatically, on every platform —
// not just wherever a matching system font happens to be installed.
const FONT_BYTES: &[u8] = include_bytes!("../assets/fonts/NotoSansMono-Regular.ttf");

#[macroquad::main("Mathematical Billiards")]
async fn main() {
    let font = load_ttf_font_from_bytes(FONT_BYTES).expect("bundled font must parse");
    // Pre-rasterize every non-ASCII glyph this codebase actually draws (ASCII
    // is already cached by `load_ttf_font_from_bytes`), so the first frame
    // that uses them doesn't stall on cache misses.
    font.populate_font_cache(&"·—⁺⁻↑↓∇−≥⊥₁½²₂₃αβθλΛπ".chars().collect::<Vec<_>>(), 18);
    set_default_font(font);

    App::new().run().await;
}
