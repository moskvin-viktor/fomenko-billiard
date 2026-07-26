# Mathematical Billiards

Interactive 2D billiard simulation inside confocal quadric domains, built with [macroquad](https://github.com/not-fl3/macroquad). Trajectories follow straight-line reflections off the boundary, and the system conserves two integrals of motion — energy `H` and the confocal constant `Λ`.

## Features

- **Confocal quadrilateral domains** — boundaries formed by two ellipse arcs and two hyperbola arcs from the same confocal family (`a = 4, b = 1`).
- **Caustic overlay** — visualize the confocal quadric `Q_Λ(x,y) = 0` that every trajectory is tangent to.
- **Real-time controls** — adjust the caustic parameter `Λ` with arrow keys, cycle presets with Tab/Space.
- **L-shaped and square polyline domains** — also supported for comparison.
- **Colored trajectory rendering** — hue sweeps from warm to cool along each bounce sequence.

## Controls

| Key | Action |
|-----|--------|
| `↑` / `→` | Increase caustic parameter Λ |
| `↓` / `←` | Decrease caustic parameter Λ |
| `Tab` / `Space` | Cycle to next domain preset |

## Presets

| Label | Domain Type |
|-------|-------------|
| Square: ellipse λ=0, hyperbola λ=2.5 | Confocal quadrilateral |
| Thin: ellipse λ=-1, hyperbola λ=2.8 | Confocal quadrilateral |
| Flat: ellipse λ=0.5, hyperbola λ=2.2 | Confocal quadrilateral |
| L-shape (polyline) | L-shaped polygon |
| Square (polyline) | Axis-aligned square |

## Running

```sh
cargo run --release
```

Tests:

```sh
cargo test
```