# Pseudo-integrable tables in code: level classification and 3D rendering

This doc explains how the crate turns a confocal table with reflex (`3π/2`)
corners into the objects the app actually draws.  The physics is in
[`confocal_L_pseudo_integrable.md`](confocal_L_pseudo_integrable.md); here we
walk the code path end to end — from the `Domain` to the pixels — and cover the
two renderers (torus vs. genus-2) with the reasons each shape was chosen.

**TL;DR**
- A table is encoded as a **flat occupancy grid** in the `(λ₁, λ₂)` chart
  (`Table`).
- For a fixed caustic value `λc`, `classify_level` decides the *topology* of the
  accessible region: **forbidden / separatrix / torus / genus-≥2**, and returns
  the flat moduli.
- Each phase point maps to flat coordinates `(u₁, u₂)` plus a momentum sheet
  (`to_flat`).
- The 3D view draws torus levels as a **donut** (flat rectangle → normalized
  angles) and genus-2 levels as a **double torus / pretzel** (two tori glued),
  so the genus actually reads in 3D.

---

## 1. The data: `Table` (λ-chart occupancy grid)

`src/table.rs`.  Every confocal table (`λ`-walls are confocal quadrics) is, in
the `(λ₁, λ₂)` chart, a union of axis-aligned cells:

```
Table {
  ell: Vec<f32>,              // ellipse-λ boundaries, sorted, ell[0] = 0
  hyp: Vec<f32>,              // hyperbola-λ boundaries, sorted
  cells: Vec<Vec<bool>>,      // cells[i][j] = [ell[i],ell[i+1]]×[hyp[j],hyp[j+1]] occupied
  quadrants: Vec<Vec<Quadrant>>, // fold label (sign(x),sign(y)) per cell
}
```

Key facts used throughout:

- **Reflex corners fall out of a local 2×2 test.**  At an interior grid vertex
  `(ell[i], hyp[j])` the table occupies `m` of the 4 surrounding cells.  `m == 3`
  ⟺ a reflex (`3π/2`) corner.  `Table::reflex_vertices()` returns all of them.
  This is the doc §1/§4 "how many quadrants the table occupies" rule, and it is
  *shape-agnostic* — L, T, Z, plus all share the same code.
- **Membership is analytic.**  `Table::contains(x, y, cf)` maps the point to
  `(λ₁, λ₂)` + a quadrant and tests the containing cell.  This replaces the
  ray-cast `Domain::contains`, which is unreliable for in-quadrant non-convex
  tables (it was the source of the "far branch" false positives).
- The standard L preset is built by `domain::confocal_lshape_standard(α₁, α₂,
  β₁, β₂, β₃)` (a 6-arc table) and recovered from any such domain by
  `Table::from_domain`.

### Quadrants as fold labels

`Quadrant::{PP,PN,NP,NN}` = `sign(x), sign(y)`.  The λ-chart is a double cover
of the plane: one `(λ₁, λ₂)` maps to two physical points related by
`(x,y)→(−x,−y)`.  The quadrant disambiguates which copy.  For the in-quadrant
standard L every cell is `PP`; a cross-quadrant table carries the fold in
`quadrants`.  Connected-component labeling respects quadrants (cells in
different fold copies are never adjacent).

---

## 2. Level topology: `pseudo::classify_level`

`src/pseudo/classify.rs`.  Given a caustic `λc`, a `Table`, and the confocal
family, this returns the topology of the accessible set:

```rust
pub enum Level {
    Forbidden,          // λc past the outer hyperbola wall — nothing accessible
    Separatrix,         // λc = b — periods diverge
    Torus { u1, u2, shape },            // single rectangle → genus 1
    GenusSurface { u1, u2, region },    // not a rectangle → genus ≥ 2
}
```

Algorithm (general — no L-specific branching):

1. **Cut by the caustic.**  Insert `λc` as a grid line: elliptic caustic keeps
   cells with `λ₁ ≤ λc`; hyperbolic keeps `λ₂ ≥ λc`.  This is the doc §7
   "caustic truncates the region" rule.
2. **Rectangle test.**  If the surviving cells form a single rectangle (contiguous,
   no holes, one fold quadrant) → `Torus`.  This is the "no reflex corner
   accessible" case.
3. **Component labeling.**  Otherwise, 4-neighbour BFS (quadrant-aware) splits
   the accessible cells into components.
4. **Count surviving reflex corners per component.**  A reflex vertex survives
   if exactly 3 of its 4 adjacent refined cells are accessible.  Then
   `genus = 1 + n_reflex` per component (doc §4).

The returned `Torus::shape` is the flat rectangle's size `(w1, w2)`; the
`GenusSurface::region` (`AccessibleRegion`) carries the flat grid lines `u1`,
`u2`, the accessibility/component masks, and `reflex_by_component`
(flat positions of the cone points where the genus is made visible).

**Why this is right for the L:** the genus is a function of `λc`.  For the
standard L (`α₁=0.4, α₂=0.8, β₁=1.4, β₂=2.0, β₃=2.6`):

| λc | accessible set | level |
|---|---|---|
| `0 < λc < α₁` | vertical strip, no corner | torus |
| `α₁ < λc < b` | L, reflex corner accessible | genus 2 |
| `b < λc < β₁` | whole L, 1 reflex corner | genus 2 |
| `β₂ < λc < β₃` | horizontal strip (base gone), no corner | torus |
| `λc ≥ β₃` | — | forbidden |

The topology snaps exactly when the caustic passes a reflex corner
(bifurcation at `λc = α₁` and `λc = β₂`).  These cases are pinned in
`tests/pseudo_topology.rs`.

---

## 3. The flat chart: `u(λ) = ∫ dλ/√P`

`src/pseudo/ulength.rs`. The reason a genus-2 level is genuinely flat is that
the natural coordinates are the **un-normalized length integrals**

```
u₁(λ₁) = ∫₀^λ₁ dλ/√P ,   u₂(λ₂) = ∫_{β₁}^{λ₂} dλ/√P ,   P = (a−λ)(b−λ)(λc−λ)
```

`ULength` tabulates one such `uᵢ(λ) = ∫ dλ/√P` over a width-knot interval
`[lo, hi]`, integrating on a grid uniform in `s = √(r − λ)` where `r` is the
root of `P` nearest the interval.  That removes the inverse-square-root
endpoint singularity when the endpoint is a turning line and is harmless when it
isn't (a wall).  Two anchoring cases:
- `RootSide::Hi` — nearest root at/above `hi` (elliptic turning, or a wall with
  the root above);
- `RootSide::Lo` — nearest root at/below `lo` (hyperbolic turning).

`to_flat` (`src/pseudo/map.rs`) per phase sample:
1. confocal coords `(λ₁, λ₂)` and `λc`;
2. `u₁ = u1.u(λ₁)`, `u₂ = u2.u(λ₂)` using the level's tables;
3. momentum sheet `(sign d₁, sign d₂)` from the implicit λ̇ formulas;
4. component via the cell containing `(u₁, u₂)`.

It returns a `FlatSample { u1, u2, sheet, fold, component }`, or `FlatError`
(`Forbidden` / `Separatrix`) for non-physical levels.

---

## 4. Rendering: `src/pseudo/render.rs` and `src/pseudo/view.rs`

The torus path (`phase3d`/`torus_render`) forces every level onto a donut via
`torus_embed`.  For the L that is wrong on both counts the user reported: the
torus levels were fed through the full-ellipse fallback (Case A, split by
angular momentum) so no clean torus appeared, and genus-2 levels were smashed
onto a donut.  The pseudo renderers replace that with embeddings chosen by the
level.

### Torus levels → donut (`torus_angles`)

A torus level is a flat rectangle `[0, w₁] × [0, w₂]` with opposite edges
identified.  Normalize each coordinate by its width to get the two circle
phases of the torus:

```rust
torus_angles(u1, u2, w1, w2) -> (θ1, θ2) = (2π·u1/w1, 2π·u2/w2)   // mod 2π each
```

This is the doc's "revert to normalized angles on the overlap" (§12).  The
two phases embed onto the standard donut (same as `phase3d::torus_embed`),
so integral levels of the L show a true torus.

### Genus-2 levels: double torus / pretzel (`pretzel_embed`)

A flat **cross** embedded in 3D reads as "4 planes" (stacked or coplanar) — it
does not *look* like genus 2.  So for `Level::GenusSurface` we embed the flat
cross onto a **double torus**: two torus lobes joined by a bridge (a thickened
figure-eight), which has two independent handle cycles — it reads as genus 2
(two tori glued).

```rust
pretzel_embed(u1, u2, σ1, σ2, moduli) -> Vec3
```

- `u2` → tube angle around the cross-section (the S¹ holonomy), mirrored by `σ2`;
- `u1` → a longitudinal loop through both lobes (`moduli` scales it), `σ1`
  mirrors the figure;
- `moduli` is `CrossModuli { a1, a2, b1, b2 }` = `u(α₁), u(α₂), u(β₂), u(β₃)`,
  the flat extents of the two legs.

The four momentum sheets are color-coded.  The per-sheet coloring
distinguishes the layers, but all lie on the same connected transverse surface
(not stacked planes).

> **Honest caveat (from the doc §12).**  Any 3D embedding of genus 2 is a
> *visualization* choice — a pretzel with the handles is arbitrary and hides
> the flat translation structure.  The pretzel makes the genus *apparent* (two
> tori glued) at the cost of the clean flat chart.  The flat chart (`cross_embed`)
> is still available for the mathematically faithful picture of the level set.

### Sampling and drawing

`view::sample_flat_trajectory(domain, p0, v0, max_steps, per_seg, level)` traces
a physical billiard trajectory and maps every interior point with `to_flat`,
producing `Vec<FlatPhasePoint> { u1, u2, sheet, component }`.  `draw_flat`
takes the collected trajectories and the level and emits the correct 3D picture
(donut for torus, pretzel for genus-2), projected by the orbit camera.

---

## 5. App wiring

`src/app.rs`.  On a 3D rebuild, the app tries to build a `Table` from the
current domain:

```
Table::from_domain(domain, &cf)
    └─ classify_level(λc, &table, &cf) ─ Level
            ├─ Torus / GenusSurface → sample_flat_trajectory → draw_flat
            └─ (Forbidden / Separatrix / not a table) → fall back to torus path
```

Concretely in `ViewState::rebuild` (the `show_3d` branch):

1. `Table::from_domain(domain, &cf)` — if the domain is a multi-cell table
   (the L), continue; else the torus path.
2. `classify_level(second_int, &tab, &cf, sep_eps=1e-9, cut_eps=1e-9)`.
3. For `Torus` or `GenusSurface`, sample all dense caustic starts through the
   *flat* chart (`sample_flat_trajectory`) into `flat_trajectories`, record the
   `Level`, and clear the torus geometry.
4. In `App::draw`, if `flat_level` is present, `pseudo::draw_flat` renders the
   flat geometries; otherwise the cached torus renderer runs unchanged.

So the pseudo path is fully *additive*: quadrilateral (full torus) and polyline
domains keep the original torus pipeline; only tables with reflex corners switch
to the flat chart and its level-appropriate 3D shape.

---

## 6. Tests

- `tests/pseudo_topology.rs` — torus ↔ genus-2 by reflex-corner presence;
  bifurcation values at `λc = α₁` and `λc = β₂`.
- `tests/pseudo_flat.rs` — `to_flat` mapping is finite, single-component, and the
  flat segments have slope **±1** (doc §13.2 slope lock), validating the `u`
  chart, the cubic, and the coordinates together.
- `tests/pseudo_render.rs` — `torus_angles` identifies opposite edges (a torus);
  `pretzel_embed` produces two lobes and real thickness in z (genus 2), and the
  sheets lie on one connected surface.
- `tests/standard_l_domain.rs` — the L domain traces, occupancy is exact, and
  per-arc branch-aware intersection holds (`Segment::intersect`).
- `tests/torus_perf.rs` — the orbit-performance helper (`decimation_step`,
  `camera_moved`) that keeps the 3D view interactive.