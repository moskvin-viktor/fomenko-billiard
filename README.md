# Mathematical Billiards

An interactive demonstration of three pieces of classical mechanics and
dynamical systems theory, built as a real 2D/3D billiard simulator:

- **Dynamical billiards** — a point mass moves in a straight line inside a
  domain and reflects specularly off the boundary. Simple to state, rich in
  behavior: depending on the domain's shape the same rule produces motion
  that ranges from completely predictable to ergodic.
- **The Liouville–Arnold theorem** — for a billiard inside a confocal
  quadric (an ellipse or an ellipse-and-hyperbola quadrilateral), the motion
  conserves *two* independent integrals, not just one. Liouville's theorem
  says that's enough to confine every trajectory to a torus in phase space,
  foliating it completely. The app's 3D view renders that torus directly —
  press `P` and watch a 2D trajectory turn into a curve winding around a
  donut.
- **Fomenko's molecules** — push the same idea to a table with a *reflex*
  corner (an L-shape) and the tori break: the invariant surfaces become
  genus-2 (double tori), and the classical Liouville picture is no longer
  enough. Sweeping the conserved quantity Λ through its full range traces
  out a **Fomenko molecule** — a Reeb graph recording how the invariant
  surface's topology bifurcates — which the app draws as a strip along the
  bottom of the screen, live, as you move the slider.

The math behind all three is written up in [`docs/math/`](docs/math/README.md);
this README covers running the app.

## Running

```sh
cargo run --release
```

```sh
cargo test
```

## Controls

| Key / mouse | Action |
|---|---|
| `↑` / `→` | Increase the second integral (Λ, or the launch angle on a polyline table) |
| `↓` / `←` | Decrease it |
| drag the slider | Same, continuously |
| `Tab` / `Space` | Cycle to the next domain preset |
| `P` | Toggle the 3D phase-space (torus / molecule) view |
| right-drag (3D view) | Orbit the camera |
| `A` | Toggle animation (sweeps Λ automatically) |
| `M` | Toggle the molecule strip overlay |
| `[` / `]` | Snap to the previous / next special (critical) layer |

## Domain presets

| Label | Domain type |
|---|---|
| Square: ellipse λ=0, hyperbola λ=2.5 | Confocal quadrilateral |
| Thin: ellipse λ=-1, hyperbola λ=2.8 | Confocal quadrilateral |
| Flat: ellipse λ=0.5, hyperbola λ=2.2 | Confocal quadrilateral |
| L-shape (confocal, 3π/2 corner) | Pseudo-integrable L (reflex corner → genus 2) |
| Ellipse (confocal, full) | Single confocal ellipse |
| Square (polyline) | Plain axis-aligned square, for comparison |
| L-shape (polyline) | Plain L-shaped polygon, for comparison |

The last two (plain polygons, no confocal structure) are there as a control
group: their 3D view is a raw `(x, y, θ/π)` embedding rather than a Liouville
torus, so switching between them and the confocal presets makes it visible
just how much the two extra structures (confocality, the second integral)
buy you.

## Documentation

- [`docs/math/`](docs/math/README.md) — the theory: the confocal family, the
  two integrals of motion, the Liouville torus derivation, and the
  pseudo-integrable genus-2 / Fomenko-molecule generalization for tables
  with a reflex corner. See [`docs/math/references.md`](docs/math/references.md)
  for the literature it builds on.
- [`docs/implementation/`](docs/implementation/architecture.md) — the code:
  module layout and how the math above gets turned into what's on screen.

Build a browsable HTML copy of all of it with [mdBook](https://rust-lang.github.io/mdBook/):

```sh
cargo install mdbook mdbook-katex mdbook-mermaid
mdbook-mermaid install docs   # one-time: vendors the mermaid JS into docs/
mdbook build docs             # writes docs/book/index.html
mdbook serve docs             # or: live-reloading local server
```

## License

MIT — see [`LICENSE`](LICENSE).
