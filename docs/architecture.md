# Architecture

## Module Layout

```
src/
├── main.rs        — App entry, rendering, camera, input handling
├── domain.rs      — Domain & Segment types, ray intersection, reflection
├── quadratic.rs   — ConfocalQuadric type, intersection, reflection, geometry
├── presets.rs     — Predefined domain configurations
tests/
└── caustic_tests.rs — Integration tests for caustic start-point validity
```

## Module Dependencies

```
presets.rs  →  domain.rs  →  quadratic.rs
main.rs     →  domain.rs, quadratic.rs, presets.rs
tests       →  domain.rs, quadratic.rs, presets.rs
```

No circular dependencies. `quadratic.rs` is leaf — no project imports of its own.