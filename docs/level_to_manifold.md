# Level → Manifold: A Clean Architecture

## The Core Idea

Replace the current fragmented classification with a single table:

```
λ → a function that produces the manifold
```

The app asks "what is the manifold at λ?" and gets back a **closed curve** (1D)
or a **torus surface** (2D) that it can sample and render. No enums, no regime
classifiers, no cascade of if-else branches.

## The Flow

```
Domain + λ
    │
    ▼
LevelClassifier::classify(domain, λ)
    │
    ├── pseudo-integrable (L-table) → Level::Torus | GenusSurface | Forbidden
    │
    └── integrable (confocal) → what kind of torus?
              │
              ├── λ < b  (elliptic, 2 tori)  → 2 circles or 2 tori
              ├── λ = b  (separatrix)        → figure-eight
              ├── b < λ < a  (hyperbolic, 1) → 1 circle or 1 torus
              ├── λ = a  (focal axis)        → 1 circle
              └── λ = λ_ell / λ_hyp (walls)  → 1 circle (the wall itself)
```

The output is always a **flat object** that the renderer knows how to draw:

```rust
pub enum Manifold {
    /// A 2D torus surface, sampled as a dense point cloud.
    Torus2D {
        /// Points on the torus, each with torus angles and an index.
        points: Vec<PhasePoint>,
    },
    /// A closed 1D curve (circle or figure-eight) on the torus surface.
    Curve1D {
        /// The curve as a list of 3D points, embedded via torus_embed.
        /// Multiple disconnected curves are concatenated with a NaN sentinel.
        points: Vec<Vec3>,
        /// How many curves (1 for a circle, 2 for a figure-eight).
        n_components: u32,
    },
    /// A flat surface (pseudo-integrable table).
    FlatSurface {
        points: Vec<FlatPhasePoint>,
        level: pseudo::Level,
    },
    /// Nothing at this level.
    Empty,
}
```

## The Level → Manifold Table

The whole thing is one function:

```rust
fn build_manifold(domain: &Domain, lam: f32, config: &SampleConfig) -> Manifold
```

### For pseudo-integrable tables (L, T, Z)

Delegates to the existing `pseudo::classify_level` + `sample_flat_trajectory`.

### For integrable (confocal) levels

1. **Build the torus mapping** using `to_torus` with the domain's bounds.
2. **Sample the caustic** via `caustic_starts` (or `critical_caustic_starts` for
   degenerate levels).
3. **For each start point, trace a trajectory** and map every sample through
   `to_torus` → `PhasePoint { theta1, theta2, torus_index }`.

The manifold type falls out of the torus mapping results:

#### λ < b (elliptic, 2 tori)
- `caustic_starts` returns points on both components of the caustic.
- `to_torus` maps them to 2 distinct `torus_index` values (0 and 1).
- **Dense sampling** → `Manifold::Torus2D` with points on both tori.
- **At a degenerate level** (λ = λ_ell wall), the caustic is the wall itself:
  `critical_caustic_starts` + `to_torus` → θ₁ is constant, θ₂ circulates.
  The result is a **curve** on the torus surface → `Manifold::Curve1D` with
  `n_components=2` (one circle per torus, since there are 2 tori).

#### λ = b (separatrix)
- `to_torus` returns `u32::MAX` for every sample (the sentinel).
- The torus degrades to a **figure-eight**: two circles meeting at the pinch
  point.  Build by sampling short trajectories and mapping them through
  `to_torus` with `sep_eps` set to allow the mapping to work arbitrarily
  close to the separatrix, then clustering by `torus_index`.
- → `Manifold::Curve1D` with `n_components=2`, spatially touching.

#### b < λ < a (hyperbolic, 1 torus)
- `caustic_starts` returns points on the single hyperbola arc.
- `to_torus` maps all to `torus_index=0`.
- **Dense sampling** → `Manifold::Torus2D` with points on one torus.
- **At λ = a or λ = λ_hyp**: `critical_caustic_starts` + `to_torus` →
  θ₁ circulates, the torus is a single curve → `Manifold::Curve1D` with
  `n_components=1`.

## What Replaces What

| Current | New |
|---------|-----|
| `bifurcation::PhaseManifold` | `Manifold` (flat, no nested classifiers) |
| `bifurcation::DegenerateKind` | gone — the manifold type is determined by `to_torus` output, not by which wall |
| `bifurcation::classify()` | `build_manifold()` — returns a renderable object |
| `critical_phase_sheet` / `CriticalSheet` | `Manifold::Curve1D` — built by the same `build_manifold` function |
| `TorusRegime` | still useful for the torus-index split rule, but not as a manifold classifier |
| `ViewState::rebuild` cascade | single `match build_manifold(domain, lam) { ... }` |

## The Key Insight

The torus mapping `to_torus` already knows everything.  At a degenerate level
the torus angles collapse and the `torus_index` tells you how many circles
there are.  The manifold is **read off the torus mapping output**, not
classified by a separate decision tree.  This means:

- No special cases for walls vs focal segment vs focal axis.
- The number of circles is always `n_tori` from `TorusRegime::n_tori()`.
- The figure-eight at λ = b is just two circles that share a point (the pinch).
- The embedding uses `torus_embed` with the same scale/offset as the regular
  tori, so the circles appear at the same position as the collapsed tori.