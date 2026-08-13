// The billiard's 2D drawing (Camera, CachedDomain, draw_*) lives in the
// `render` module, the second-integral range math in the lib
// (`second_integral`), the drag widget in `ui`, and the app loop in `app`.
// This file is a thin shell: construct the app and run it.

use billiards::app::App;

#[macroquad::main("Mathematical Billiards")]
async fn main() {
    App::new().run().await;
}
