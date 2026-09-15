# References

The background literature the math docs in this folder build on.

## Liouville integrability & Fomenko's theory of atoms/molecules

- V. I. Arnold, *Mathematical Methods of Classical Mechanics*, 2nd ed.,
  Springer, 1989. — The Liouville–Arnold theorem: a compact invariant
  manifold of a completely integrable system is a torus, and the motion on
  it is quasi-periodic in action-angle coordinates. This is the theorem the
  whole 3D torus view in this project is a picture of.
- A. T. Fomenko and A. V. Bolsinov, *Integrable Hamiltonian Systems:
  Geometry, Topology, Classification*, CRC Press, 2004. — The source for
  "atoms" (the local structure of a singular fiber of the momentum map) and
  "molecules" (the Reeb graph of the whole foliation, decorated with atoms).
  `docs/math/bifurcation_atoms.md` classifies the confocal L's critical
  fibers in exactly this language.

## Confocal quadric billiards

- S. Tabachnikov, *Geometry and Billiards*, Student Mathematical Library
  vol. 30, AMS, 2005 (freely available from the author's website). —
  Chapter on billiards in conics: the confocal-quadric caustic family, the
  two integrals of motion, and the classical Jacobi/Chasles construction
  this project's `confocal_quadrics.md` and `liouville_tori.md` follow.
- V. V. Kozlov and D. V. Treshchëv, *Billiards: A Genetic Introduction to
  the Dynamics of Systems with Impacts*, Translations of Mathematical
  Monographs vol. 89, AMS, 1991.

## Pseudo-integrable billiards & flat surfaces

- A. N. Zemlyakov and A. B. Katok, "Topological transitivity of billiards
  in polygons," *Mathematical Notes* 18 (1975), 760–764. — The origin of
  polygonal billiards with irrational angles unfolding to a flat
  (translation) surface of genus > 1.
- P. J. Richens and M. V. Berry, "Pseudointegrable systems in classical and
  quantum mechanics," *Physica D* 2:3 (1981), 495–512. — Coined
  "pseudo-integrable": as many independent integrals as an integrable
  system, but level sets of higher genus instead of tori — exactly the
  regime `docs/math/pseudo_integrable_genus.md` and `bifurcation_atoms.md`
  work out for a confocal table with a reflex corner.
- A. Zorich, "Flat Surfaces," in *Frontiers in Number Theory, Physics, and
  Geometry I*, Springer, 2006, 439–586. — Survey of translation surfaces
  and the strata `H(2)`, `H(0)` referenced when the genus-2 level
  degenerates to a torus in `bifurcation_atoms.md`.

This list is deliberately short — the docs cite the specific results they use
inline; these are the books/papers to read for the surrounding theory.
