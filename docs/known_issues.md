# Known Issues & Remaining Work

Status of open problems and intentional follow-ups. The **first item is the only
known *bug* with a failing test**; the rest are structural/debt items.

## 1. L-shape 2D/3D divergence — the only failing test (BUG)
**Failing test:** `billiard_picks_are_consistent_with_torus_highlights`
(`tests/highlight_invariants.rs`).

At an **ellipse caustic** on the **L-shape** confocal preset, the 2D view and
the 3D torus disagree on the number of tori:

```
L-shape (confocal, 3π/2 corner) Λ=0.49999997:
2D picks 2 keys [1, 0] but 3D shows 1 tori [1];
empty-phase starts: [(1, (-1.60,-0.37)), (3, (1.83,-0.15))]
```

**Root cause(s):**
- `ViewState::rebuild` builds `torus_highlights` by tracing **only 4 bounces** and
  then silently `.filter(|t| !t.is_empty())`. For two of the L-shape's four
  caustic starts the 4-bounce phase trace comes back **empty** (start point too
  near a corner / degenerate), so those tori are dropped from the 3D view.
- Meanwhile the billiard's `one_per_torus` (in `render.rs`) counts *all* the
  starts (including the ones that trace empty highlights), so it still sees 2
  distinct tori. The two views diverge by construction.

**Likely fix direction:** make the 2D picker and the 3D highlights derive keys
from a *single* source of truth — pick one trajectory per torus **off the
filtered highlights** (the `torus_index` actually present) rather than counting
pre-filter starts. `ConfocalStructure` (added today) now centralizes the
`regime`/`contains`/`caustic_starts` logic, giving this fix one home.

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