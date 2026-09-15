# Domain Reference

## `Segment` enum

Represents one smooth boundary piece.

### Variants

| Variant | Fields | Description |
|---------|--------|-------------|
| `Line` | `a: Vec2, b: Vec2` | Straight wall segment |
| `Quad` | `curve: ConfocalQuadric, a: Vec2, b: Vec2` | Confocal quadric arc |

### Methods

| Method | Signature | Description |
|--------|-----------|-------------|
| `inward_normal` | `(&self, p: Vec2) -> Vec2` | Unit normal pointing into the domain |
| `reflect` | `(&self, p: Vec2, dir: Vec2) -> Vec2` | Specular reflection at hit point |
| `intersect` | `(&self, p: Vec2, dir: Vec2) -> Option<(t, hit, s)>` | Ray-segment intersection. `t` = distance, `s` = parameter along segment (0..1) |

## `Domain` struct

| Field | Type | Description |
|-------|------|-------------|
| `segments` | `Vec<Segment>` | Ordered CCW boundary segments |

### Methods

| Method | Description |
|--------|-------------|
| `contains(p)` | Point-in-domain test via ray casting |
| `intersect(p, dir)` | Closest intersection with any segment → `(t, segment_idx, hit)` |
| `reflect(p, dir, idx)` | Reflect off segment `idx`, with π/2 corner special case |
| `trace(p, v, max_steps)` | Full trajectory: bounce `max_steps` times or until escape |
| `sample_boundary(n)` | Points for visual rendering, in boundary order |
| `corners()` | Junction points between consecutive segments |

## Constructor Functions (in `domain.rs`)

| Function | Description |
|----------|-------------|
| `confocal_quad(a, b, λ_ell, λ_hyp)` | Confocal quadrilateral from ellipse + hyperbola |
| `square()` | Unit square `[−1, 1]²` |
| `lshape_poly()` | L-shape with vertices `(0,0)→(2,0)→(2,1)→(1,1)→(1,2)→(0,2)` |