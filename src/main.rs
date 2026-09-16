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

// WebGL1 is macroquad's default on web, and miniquad 0.4.10's WebGL1 backend
// calls `gl.readBuffer` (a WebGL2-only method) when finishing a render pass
// on an offscreen render target — which is exactly how the 3D torus/molecule
// view is drawn (`cached_render::CachedSurfaceRender`). That throws an
// uncaught JS TypeError out of the wasm call, which never unwinds cleanly and
// leaves miniquad's internal event-handler lock permanently poisoned — so the
// 3D view silently renders nothing, and every input event afterward panics
// with "already borrowed". Forcing WebGL2 avoids the buggy WebGL1 code path
// entirely; it's supported by effectively all browsers macroquad targets.
fn window_conf() -> Conf {
    macroquad::miniquad::conf::Conf {
        window_title: "Mathematical Billiards".to_owned(),
        platform: macroquad::miniquad::conf::Platform {
            webgl_version: macroquad::miniquad::conf::WebGLVersion::WebGL2,
            ..Default::default()
        },
        ..Default::default()
    }
    .into()
}

#[macroquad::main(window_conf)]
async fn main() {
    let font = load_ttf_font_from_bytes(FONT_BYTES).expect("bundled font must parse");
    // Pre-rasterize every non-ASCII glyph this codebase actually draws (ASCII
    // is already cached by `load_ttf_font_from_bytes`), so the first frame
    // that uses them doesn't stall on cache misses.
    font.populate_font_cache(&"·—⁺⁻↑↓∇−≥⊥₁½²₂₃αβθλΛπ".chars().collect::<Vec<_>>(), 18);
    set_default_font(font);

    App::new().run().await;
}
