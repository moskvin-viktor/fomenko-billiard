# Architecture

## Module Layout

```
src/
├── main.rs            — App entry (thin shell constructing & running App)
├── lib.rs             — Module re-exports, caustic start-point selection
├── app.rs             — App loop: input → ViewState → draw
├── ui.rs              — Slider widget
├── confocal.rs        — ConfocalStructure: domain→λ analysis, TorusRegime, caustic starts
├── domain.rs          — Domain & Segment types, ray intersection, reflection,
│                          L-shape builders (confocal_lshape_standard)
├── quadratic.rs       — ConfocalQuadric type, intersection, reflection, geometry
├── table.rs           — Table: rectilinear λ-chart occupancy grid + Quadrant fold labels
├── presets.rs         — Predefined domain configurations
├── render.rs          — 2D billiard drawing: Camera, CachedDomain, draw_* primitives
├── second_integral.rs — Valid second-integral range + normalized slider mapping
├── phase3d.rs         — 3D phase-space sampling, torus embedding, orbit camera, drawing
├── torus_render.rs    — Cached offscreen point-cloud rendering. 3D torus
├── torus.rs           — Liouville-torus machinery (phase-space → torus mapping)
├── torus/
│   ├── mod.rs         — Module facade + re-exports
│   ├── confocal.rs    — ConfocalParams, PhaseSample, coordinate helpers
│   ├── quadrature.rs  — Libration tables, interpolation (the Abelian phase)
│   └── map.rs         — to_torus: phase-space → torus mapping (Cases A/B/C)
└── pseudo/            — Pseudo-integrable tables (reflex corners)
    ├── mod.rs         — Module facade + re-exports
    ├── ulength.rs     — ULength: u(λ)=∫dλ/√P quadrature table
    ├── classify.rs    — classify_level: Level (forbidden/separatrix/torus/genus)
    ├── map.rs         — to_flat: phase-space → flat coordinates + sheet
    ├── render.rs      — torus_angles, cross_embed, pretzel_embed (3D embeddings)
    └── view.rs        — sample_flat_trajectory, draw_flat (sampling + drawing)
tests/
├── caustic_tests.rs        — Integration tests for caustic start-point validity
├── torus_tests.rs          — Torus-manifold topology tests
├── torus_render_test.rs    — Renderer smoke test (non-blank offscreen target)
├── highlight_invariants.rs — Per-torus highlight consistency invariants
├── standard_l_domain.rs    — Standard L (in-quadrant) domain geometry + tracing
├── pseudo_topology.rs      — Torus ↔ genus-2 by reflex-corner presence; bifurcations
├── pseudo_flat.rs          — Flat-chart mapping (slope ±1, finite, single-component)
├── pseudo_render.rs        — 3D embedding helpers (torus angles, pretzel genus-2)
└── torus_perf.rs           — Orbit-performance helpers (decimation, camera throttle)
```

## Module Dependencies

```
presets.rs  →  domain.rs  →  quadratic.rs
app.rs      →  domain.rs, quadratic.rs, presets.rs, phase3d.rs, render.rs,
               ui.rs, second_integral.rs, torus.rs, torus_render.rs, pseudo, table
phase3d.rs  →  torus::{map, confocal}, torus_render
render.rs   →  domain.rs, quadratic.rs, confocal.rs (TorusRegime)
table.rs    →  torus (ConfocalParams)
pseudo/     →  table.rs, torus::{confocal}, quadratic
    classify.rs → table.rs, ulength.rs, torus
    map.rs      → classify.rs, table.rs, torus
    render.rs   → (geometry only)
    view.rs     → classify.rs, map.rs, render.rs, torus, phase3d
lib.rs      →  domain.rs, quadratic.rs, torus, pseudo, table, confocal

No circular dependencies.  quadratic.rs is leaf.  pseudo/ and table are the
pseudo-integrable extension; torus/ is the integrable (Liouville) part.
```

For the pseudo-integrable machinery and its two 3D renderers (torus vs.
genus-2), see [`pseudo_integrable.md`](pseudo_integrable.md).