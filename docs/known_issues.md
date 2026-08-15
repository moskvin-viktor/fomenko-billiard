# Known Issues & Remaining Work

Status of open problems and intentional follow-ups. Nothing here is a currently
failing test (all four suites are green); the L-shape 2D/3D divergence that used
to be item #1 was resolved by replacing the old all-quadrant `confocal_lshape`
preset with the in-quadrant standard L from `docs/confocal_L_pseudo_integrable.md`
(a real 6-arc table matching the doc's §2) and making quadric-arc intersection
branch-aware + conic-angle-parameterized.

## ~~Orbit freeze in 3D~~ (RESOLVED)
Orbiting/zooming the 3D phase-space view froze the app: the orbit camera re-
rasterized the full ~150K-point torus cloud into the offscreen target every
frame (~150K `draw_circle` calls).  Fixed by (1) a point-budget decimation
(`phase3d::decimation_step`, a uniform stride capping a rasterization at
`POINT_BUDGET` while preserving coverage) and (2) a camera-move threshold
(`torus_render::camera_moved`) so sub-jitter drag no longer triggers a re-
raster.

## ~~L 3D view forced onto a donut~~ (RESOLVED)
The L-shape's 3D phase-space view was routed entirely through the torus
machinery (`to_torus` full-ellipse fallback + `torus_embed` donut), so torus
levels didn't show a torus and genus-2 levels were forced onto a donut.  Now the
L's 3D view goes through the pseudo machinery (`pseudo/view.rs`):
- `Level::Torus` → flat rectangle `(u1,u2)` with opposite edges identified → a
  donut (`pseudo::render::torus_angles` + donut embed).
- `Level::GenusSurface` → **double torus (pretzel)** embedding
  (`pseudo::render::pretzel_embed`): two torus lobes joined by a bridge, so the
  manifold visibly reads as genus 2 (two tori glued).  The flat cross chart is
  still available via `cross_embed`, but the pretzel is what shows the genus in
  3D (per doc §12, any 3D genus-2 embedding is a visualization choice; a flat
  cross reads as planes).

## 2. ~~`main.rs` is still a thin-shell-plus-widgets monolith~~ (RESOLVED)
`main.rs` is now a thin shell: the drag widget moved to `src/ui.rs` (`Slider`),
the app loop + `ViewState` moved to `src/app.rs` (`App`), and `main` just
constructs and runs the app.

## 3. ~~`A` / `B` family constants — consumed but not the single clear home~~ (RESOLVED)
`ConfocalParams` is now the single first-class carrier of the confocal family
throughout the codebase. `ConfocalStructure` stores a `cf`, `quadratic`/`domain`/
`render` take `ConfocalParams` instead of `(a, b)` scalar pairs, and the dead
`_a`/`_b` params were dropped from `get_start_points`/`start_points_on_caustic`.
`ConfocalParams::standard()` is the single definition site (`a = 4.0, b = 1.0`);
the crate-root `A`/`B` constants are gone.

## 4. Remaining clippy hygiene (trivial, safe)
All pre-existing, none from recent refactors:
- `domain.rs:53` — manual `RangeInclusive::contains` → `(0.0..=1.0).contains(&s)`
- `domain.rs:133,183` — `.filter_map(..)` → `.map(..)`
- `phase3d.rs:45,62` — `OrbitCamera3` missing `Default`; min/max → `clamp`
- `quadratic.rs:14` — doc list indentation
- `torus_render.rs:28` — `TorusRender` missing `Default`
- `lib.rs` — unnecessary `let` in `confocal_inside` (should be gone now)
- Test-local warnings in `torus_tests.rs` (loop-indexed `bins`)

## 5. Duplicated boundary-lambda walk in `torus/map.rs` vs `confocal.rs`
`torus/map.rs` (`case.rs` helpers) and `lib.rs`/`confocal.rs` both walk the
`curve.b_param`/`a_param` thresholds in places. Mostly consolidated today
(`boundary_lambdas` → `ConfocalStructure`), but the seam still exists and is a
candidate for further unification.

## 6. Tests: repository convention is integration tests in `tests/`
Lib modules have no inline `#[cfg(test)]` — all behavioral tests live in
`tests/` (by design). Keep new tests there.

## 7. "One per torus" highlight selection is the thinnest-tested feature
The red-per-torus highlight (`one_per_torus`, `draw_torus_highlights`) has
invariant tests but no direct unit test that the *set* of highlighted `torus_index`
values matches `ConfocalStructure::regime`'s `n_tori()` for every preset/caustic.
Worth adding once #1 is resolved.