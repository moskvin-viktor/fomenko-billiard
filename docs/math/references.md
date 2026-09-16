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

## Billiards as integrable systems: Liouville foliations & molecules

- V. A. Moskvin, "Topology of Liouville Foliations for Integrable Billiards
  in Non-Convex Domains," *Moscow Univ. Math. Bull.* 73 (2018), 103–110.
  [doi:10.3103/S002713221803004X](https://doi.org/10.3103/S002713221803004X).
  — Liouville foliations for billiards on non-convex tables; the direct
  precursor to the reflex-corner / genus-2 case worked out for the confocal
  L in `pseudo_integrable_genus.md` and `bifurcation_atoms.md`. (By this
  project's author.)
- A. T. Fomenko and V. V. Vedyushkina, "Billiards and integrable systems,"
  *Russian Math. Surveys* 78:5 (2023), 881–954. — Survey of the
  billiard-as-integrable-system program this whole project is a worked
  example of: how a billiard's Liouville foliation and molecule are built,
  and where pseudo-integrable (reflex-corner) tables fit in.
- A. T. Fomenko and V. A. Kibkalo, "Topology of Liouville foliations of
  integrable billiards on table-complexes," *European Journal of
  Mathematics* 8:4 (2022), 1392. — Billiard "table-complexes" glued from
  several tables; the general framework the confocal L (as a single
  non-convex table) is a small case of.

## Confocal quadric billiards

- S. Tabachnikov, *Geometry and Billiards*, Student Mathematical Library
  vol. 30, AMS, 2005 (freely available from the author's website). —
  Chapter on billiards in conics: the confocal-quadric caustic family, the
  two integrals of motion, and the classical Jacobi/Chasles construction
  this project's `confocal_quadrics.md` and `liouville_tori.md` follow.
- V. V. Kozlov and D. V. Treshchëv, *Billiards: A Genetic Introduction to
  the Dynamics of Systems with Impacts*, Translations of Mathematical
  Monographs vol. 89, AMS, 1991.
