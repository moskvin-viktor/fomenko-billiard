# Mathematics

## Confocal Quadric Family

The billiard uses a confocal family of quadrics parameterized by `λ`:

```
(b − λ)·x² + (a − λ)·y² = (a − λ)·(b − λ),   λ ≤ a
```

Fixed constants: `a > b > 0`. Foci at `(±c, 0)` with `c² = a − b`.

In this project: `a = 4, b = 1`.

### Regime by λ

| λ range | Shape |
|---------|-------|
| `λ < b` | Ellipse |
| `λ = b` | Degenerate: segment `[−c, +c]` + horizontal rays |
| `b < λ < a` | Hyperbola (opens left/right) |
| `λ = a` | Vertical segment `x = 0, −√b ≤ y ≤ √b` |

## Integrals of Motion

For a billiard inside a confocal quadric, each trajectory conserves:

1. **Energy:** `H = ½|v|²` (speed remains constant after elastic reflection).
2. **Confocal constant:** `Λ = vx²/a + vy²/b − (x·vy − y·vx)² / (ab)`.

Every segment of the trajectory is tangent to the caustic quadric `Q_Λ(x,y) = 0`.

## Reflection Law

At each boundary point, the incoming velocity reflects specularly:

```
v' = v − 2(v · n̂) n̂
```

Where `n̂` is the **inward** normal at the hit point. For quadric arcs, the normal is derived from the gradient of `Q`:

```
n ∝ −∇Q = −( 2(b−λ)x, 2(a−λ)y )
```

## Caustic Tangent Velocity

Given a point `p` on the caustic `Q_Λ = 0`, the tangent direction (velocity) is perpendicular to `∇Q_Λ(p)`:

```
v ∝ ( −(a−Λ)·py,  (b−Λ)·px )
```

The sign is chosen to align with a `hint_dir` (default `(0, 1)`).

## Confocal Quadrilateral

A domain bounded by two ellipses can't intersect (they are nested), so instead the boundary uses:

- One ellipse arc (fixed `λ_ell < b`) — top and bottom.
- One hyperbola arc (fixed `λ_hyp > b`) — left and right.

These intersect at 4 points:

```
x² = (a − λ_ell)(a − λ_hyp) / (a − b)
y² = (b − λ_ell)(b − λ_hyp) / (b − a)
```

The 4 intersection points are: `tr (x, y)`, `br (x, −y)`, `bl (−x, −y)`, `tl (−x, y)`.

Boundary arcs (CCW order):
1. Top: ellipse `tr → tl`
2. Left: hyperbola `tl → bl`
3. Bottom: ellipse `bl → br`
4. Right: hyperbola `br → tr`