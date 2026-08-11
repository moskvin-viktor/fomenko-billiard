# Architecture

## Module Layout

```
src/
├── main.rs            — App entry, rendering, camera, input handling
├── lib.rs             — Constants, ConfocalStructure/TorusRegime re-exports, start-point selection
├── confocal.rs        — ConfocalStructure: domain→λ analysis, TorusRegime, caustic run-walk
├── domain.rs          — Domain & Segment types, ray intersection, reflection
├── quadratic.rs       — ConfocalQuadric type, intersection, reflection, geometry
├── presets.rs         — Predefined domain configurations
├── render.rs          — 2D billiard drawing: Camera, CachedDomain, draw_* primitives
├── second_integral.rs — Valid second-integral range + normalized slider mapping
├── phase3d.rs         — 3D phase-space sampling, torus embedding, drawing
├── torus_render.rs    — Cached offscreen rendering of the 3D torus
└── torus/
    ├── mod.rs         — Module facade + re-exports
    ├── confocal.rs    — ConfocalParams, PhaseSample, coordinate helpers
    ├── quadrature.rs  — Libration tables, interpolation (the Abelian phase)
    └── map.rs         — to_torus: phase-space → torus mapping (Cases A/B/C)
tests/
├── caustic_tests.rs        — Integration tests for caustic start-point validity
├── torus_tests.rs          — Torus-manifold topology tests
├── torus_render_test.rs    — Renderer smoke test (non-blank offscreen target)
└── highlight_invariants.rs — Per-torus highlight consistency invariants
```

## Module Dependencies

```
presets.rs  →  domain.rs  →  quadratic.rs
main.rs     →  domain.rs, quadratic.rs, presets.rs, lib.rs, phase3d.rs, render.rs
phase3d.rs  →  torus::{map, confocal}, torus_render
render.rs   →  domain.rs, quadratic.rs, lib.rs (TorusRegime)
lib.rs      →  domain.rs, quadratic.rs, torus

Inside torus/:  map.rs → confocal.rs, quadrature.rs
                mod.rs re-exports the three leaf modules.
```

No circular dependencies. `quadratic.rs` is leaf — no project imports of its own.