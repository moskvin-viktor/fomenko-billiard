# Mathematics

The theory behind the simulation: confocal quadric billiards, their integrals
of motion, the Liouville tori this gives rise to, and — for tables with a
reflex corner — the pseudo-integrable, genus-2 generalization and its Fomenko
molecule.

Read in this order:

1. [`confocal_quadrics.md`](confocal_quadrics.md) — the confocal family, the
   two conserved integrals `H` and `Λ`, and the reflection law. Start here.
2. [`liouville_tori.md`](liouville_tori.md) — the full derivation of the
   Liouville torus for a confocal billiard: action-angle coordinates,
   the elliptic curve behind the second integral, and how many tori a given
   level has.
3. [`pseudo_integrable_genus.md`](pseudo_integrable_genus.md) — what happens
   once the table has a reflex (`3π/2`) corner: the invariant surfaces are no
   longer tori but genus-`g` surfaces, `g = 1 + n` in the number of
   (accessible) reflex corners.
4. [`bifurcation_atoms.md`](bifurcation_atoms.md) — the molecule itself: how
   the genus-`g` surfaces degenerate at critical caustic values, which
   degenerations are classical Fomenko atoms and which are genus-changing
   stratum-boundary events, and how the whole sweep assembles into a Reeb
   graph.

See [`references.md`](references.md) for the background literature these docs
build on, and [`../implementation/`](../implementation/architecture.md) for
how this theory maps onto the code.
