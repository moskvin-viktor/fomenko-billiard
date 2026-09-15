# Documentation

Two independent tracks:

- [`math/`](math/README.md) — the theory: confocal quadric billiards,
  Liouville tori, and (for tables with a reflex corner) the pseudo-integrable
  genus-2 generalization and its Fomenko molecule. Read this to understand
  *what* the simulation shows.
- [`implementation/`](implementation/architecture.md) — the code: module
  layout, dependencies, and how each piece of theory gets turned into pixels.
  Read this to understand *how* the crate is put together.

See the [project README](../README.md) for what the app does and how to run
it. This folder is also an [mdBook](https://rust-lang.github.io/mdBook/)
source (`book.toml` + `SUMMARY.md`) — run `mdbook build docs` from the repo
root to render it to a browsable static HTML site.
