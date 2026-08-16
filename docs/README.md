# Mathematical Billiards

Interactive 2D billiard simulation inside confocal quadric domains, built with [macroquad](https://github.com/not-fl3/macroquad). Trajectories follow straight-line reflections off the boundary, and the system conserves two integrals of motion — energy `H` and the confocal constant `Λ`.

## Features

- **Confocal quadrilateral domains** — boundaries formed by two ellipse arcs and two hyperbola arcs from the same confocal family (`a = 4, b = 1`).
- **Caustic overlay** — visualize the confocal quadric `Q_Λ(x,y) = 0` that every trajectory is tangent to.
- **Real-time controls** — adjust the caustic parameter `Λ` with arrow keys, cycle presets with Tab/Space.
- **L-shaped and square polyline domains** — also supported for comparison.
- **Pseudo-integrable L (confocal, 3π/2 corner)** — the standard in-quadrant L
  from `confocal_L_pseudo_integrable.md`: a real 6-arc table whose level sets are
  genus-2 (flat-coordinate, translation-surface) rather than tori.
- **Colored trajectory rendering** — hue sweeps from warm to cool along each bounce sequence.

## Controls

| Key | Action |
|-----|--------|
| `↑` / `→` | Increase caustic parameter Λ |
| `↓` / `←` | Decrease caustic parameter Λ |
| `Tab` / `Space` | Cycle to next domain preset |

## Presets

| Label | Domain Type |
| Square: ellipse λ=0, hyperbola λ=2.5 | Confocal quadrilateral |
| Thin: ellipse λ=-1, hyperbola λ=2.8 | Confocal quadrilateral |
| Flat: ellipse λ=0.5, hyperbola λ=2.2 | Confocal quadrilateral |
| L-shape (confocal, 3π/2 corner) | Pseudo-integrable L (standard, in-quadrant) |
| Square (polyline) | Axis-aligned square |
| L-shape (polyline) | L-shaped polygon |

## 3D view

The `P` key toggles the 3D phase-space view.  For integrable (quadrilateral)
levels it shows the Liouville **torus**.  For the pseudo-integrable L it first
classifies the level and then shows either a **torus** (when the caustic shadows
the reflex corner) or a **double torus / pretzel** (when the reflex corner is
accessible → genus 2).  See [`pseudo_integrable.md`](pseudo_integrable.md) for
how this works.

## Running

```sh
cargo run --release
```

Tests:

```sh
cargo test
```

See [`architecture.md`](architecture.md) for the module layout,
[`pseudo_integrable.md`](pseudo_integrable.md) for the pseudo-integrable
renderers, and [`known_issues.md`](known_issues.md) for open bugs and follow-up
work.