# Quadratic Reference

## `ConfocalQuadric` struct

| Field | Type | Description |
|-------|------|-------------|
| `a_param` | `f32` | Semi-major axis squared (a) |
| `b_param` | `f32` | Semi-minor axis squared (b) |
| `lambda` | `f32` | Confocal parameter (λ) |

### Methods

| Method | Description |
|--------|-------------|
| `c()` | Semi-focal distance `√(a − b)` |
| `eval(p)` | Evaluate `Q(p)`. Interior is `Q(p) < 0` |
| `grad(p)` | Gradient `∇Q(p)` |
| `inward_normal(p)` | Unit inward normal `−∇Q/|∇Q|` |
| `reflect(p, dir)` | Specular reflection off `Q = 0` |
| `intersect(p, dir)` | Smallest positive `t` where ray hits `Q = 0` |
| `centre()` | Always `(0, 0)` for this confocal family |
| `intersections(l1, l2)` | 4 intersection points `[tr, br, bl, tl]` of two confocal quadrics |
| `sample_boundary(n)` | Sample points by ray-casting from origin |
| `velocity_from_caustic(a, b, p, λ, hint)` | Tangent velocity on caustic, sign-picked toward `hint` |

## Free Functions (in `quadratic.rs`)

| Function | Description |
|----------|-------------|
| `confocal(a, b, λ)` | Construct a `ConfocalQuadric` |
| `quadrilateral_arcs(a, b, λ_ell, λ_hyp)` | 4 arcs `[(from, to, quadric)]` in CCW order for a confocal quadrilateral |