# Architecture

## Module Layout

```
src/
├── main.rs        — App entry, rendering, camera, input handling
├── domain.rs      — Domain & Segment types, ray intersection, reflection
├── quadratic.rs   — ConfocalQuadric type, intersection, reflection, geometry
├── torus.rs       — Liouville-torus mapping (Jacobi normalization)
├── presets.rs     — Predefined domain configurations
├── phase3d.rs     — 3D phase-space sampling, torus embedding, drawing
tests/
├── caustic_tests.rs — Integration tests for caustic start-point validity
└── torus_tests.rs   — Torus-manifold tests
```

## Module Dependencies

```
presets.rs  →  domain.rs  →  quadratic.rs
main.rs     →  domain.rs, quadratic.rs, presets.rs, torus.rs, phase3d.rs
phase3d.rs  →  torus.rs
lib.rs      →  torus_bounds (torus mapping bounds for a domain)
tests       →  domain.rs, quadratic.rs, presets.rs, phase3d.rs, torus.rs
```

No circular dependencies. `quadratic.rs` is leaf — no project imports of its own.