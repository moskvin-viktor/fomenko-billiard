# Degenerate Layer Phase-Space Embedding

## The Problem

At a degenerate (critical) caustic level — λ = b (focal separatrix), λ = a (focal
axis), λ = λ_ell (ellipse wall), λ = λ_hyp (hyperbola wall) — the phase manifold
collapses to a 1D closed orbit.  The current rendering (`CriticalSheet` /
`draw_critical_sheet`) hardcodes a single circle in the xz-plane for every
degenerate layer, which is topologically wrong:

- **Elliptical side of separatrix (λ < b)**: two tori → two disconnected circles.
- **Hyperbolic side (λ > b)**: one torus → one circle.
- **At the separatrix λ = b**: the two tori touch at a pinch point → a
  **figure-eight** (two circles meeting at a point).

The number of circles and their connectivity must be derived from the actual
phase-space mapping, not from a hardcoded shape.

## What a Correct Solution Needs

### 1. Count the circles from the torus regime

`ConfocalStructure::regime(lam)` already returns `TorusRegime`:

- `SplitByY` (elliptic, quadrilateral) → 2 tori → 2 circles.
- `SplitByL` (elliptic, full ellipse) → 2 tori → 2 circles.
- `Single` (hyperbolic) → 1 torus → 1 circle.

At the separatrix λ = b: the torus map returns `u32::MAX` (the sentinel in
`to_torus`).  There are 2 tori touching at a point → a figure-eight.

### 2. Embed each torus as a circle via the torus mapping

The torus mapping `to_torus` produces `(θ₁, θ₂, torus_index)` for any
phase-space sample.  At a degenerate level, the λ₁-libration (the elliptic
degree of freedom) collapses to a single point — only the λ₂ angle circulates.
So each torus degrades to a **circle** parametrised by the surviving angle:

```rust
// For a degenerate layer, the phase manifold projects to a circle:
//   θ₂(t) = t ∈ [0, 2π)   (the surviving angle)
//   θ₁    = constant       (the collapsed angle)
//   torus_index = 0 or 1   (which torus, from TorusRegime)
struct DegenerateCircle {
    theta2_samples: Vec<PhasePoint>,  // samples around the surviving angle
    torus_index: u32,
}
```

To build one circle for a torus:
1. Pick a start point on the caustic inside the domain (from
   `critical_caustic_starts` or `caustic_starts`).
2. Trace a short trajectory (a few bounces).
3. Record every sample's `theta2` via `to_torus` — it should sweep the full
   `[0, 2π)` range.
4. The torus_index tells you which circle it belongs to.

### 3. Embed the figure-eight at the separatrix

At λ = b, the two tori share a hyperbolic 2-orbit (the pinch point).  The
accessible region is the two tori touching at that point.  Build both circles
as above; they will have one point in common (the pinch).  Either draw both
separately and let the shared point visually connect them, or explicitly build
a single 8-shaped curve.

### 4. Embed each circle in 3D

Each circle lives in the 3D coordinate space of the phase manifold.  The
current `torus_embed` uses `(θ₁, θ₂, torus_index)` — for a degenerate level,
θ₁ is constant, so the embed is a curve on the torus surface.  Project that
curve into 3D via `torus_embed` with the same `r_major`/`r_minor` and
`torus_index` offset used for the regular tori, so the degenerate circles
appear in the same spatial position as the tori they collapsed from.

### 5. What to change

| File | What to change |
|------|----------------|
| `src/phase3d.rs` | Replace `CriticalSheet` / `build_degenerate_orbit` / `draw_critical_sheet` with a new builder that takes a `Domain`, `lam`, `TorusRegime`, produces `Vec<DegenerateCircle>`, and a new drawer that renders each circle as a polyline. |
| `src/bifurcation.rs` | `PhaseManifold::Degenerate` already carries `DegenerateKind`; the caller (`app.rs`) already matches on it. The kind is not needed for the embedding (the regime determines the number of circles), but may carry information about which wall the orbit slides along. |
| `src/app.rs` | `ViewState::rebuild` — construct the circles from `critical_caustic_starts` + `to_torus` + `TorusRegime` instead of calling `critical_phase_sheet`. |

### 6. What NOT to do

- Do NOT hardcode any shape (circle, figure-eight, line).  The shape must fall
  out of the torus mapping.
- Do NOT use the physical position `(x, y)` to determine the embedding — the
  embedding is in the torus angle space `(θ₁, θ₂)`, not in the billiard
  domain.
- Do NOT special-case the ellipse-wall and hyperbola-wall as different from
  the focal-axis / focal-segment.  All four are degenerate layers and should
  be handled by the same algorithm.

### 7. Verification

The degenerate-layer rendering must agree with the regular torus rendering at
nearby levels: as you sweep λ toward the separatrix, the two tori should
visually shrink onto two circles that meet at the pinch point.  The
`critical_layers` test suite (`tests/critical_layers.rs`) should be extended
to assert the number of circles and their finiteness.