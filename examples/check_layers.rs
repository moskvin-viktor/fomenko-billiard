fn main() {
    let cf = billiards::torus::ConfocalParams::standard();
    let presets = billiards::presets::all_presets();
    for p in presets.iter().filter(|p| p.is_confocal) {
        let layers = billiards::molecule::special_layers(&p.domain, &cf);
        println!("{}:", p.label);
        for l in &layers {
            let label = billiards::molecule::layer_label(*l, &p.domain, &cf);
            let starts = billiards::get_start_points(*l, &p.domain, true, macroquad::prelude::Vec2::ZERO);
            println!("  lam={:.4} label={} starts={}", l, label, starts.len());
        }
    }
}
