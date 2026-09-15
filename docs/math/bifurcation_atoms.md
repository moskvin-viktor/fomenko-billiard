# Bifurcation levels of the confocal L: critical fibers explicitly

**The idea, confirmed.** At a bifurcation value $\lambda_c = \lambda_*$ the invariant surface fails to be a smooth 2-manifold — but it fails *only along a finite graph of critical trajectories* (the orbits that run into the reflex corner). The rest of the critical fiber is regular. So the critical level is a **thin** object: a finite 1-dimensional spine with regular 2-dimensional pieces attached. This is exactly what lets you build its topology by hand, and this doc does that.

**The one twist worth stating up front.** In a genuinely Liouville-integrable system the regular fibers are all tori and the singular fibers are classified by *Fomenko atoms* (A, B, …). Here the regular fibers can be **genus 2** (previous doc), so the genus-changing bifurcations are **not atoms at all** — they are degenerations of the flat surface to the boundary of its stratum $\mathcal{H}(2)$ (a cylinder collapses). The genus-preserving events *are* ordinary atoms. Both are governed by a finite set of critical orbits, which is the unifying content.

Conventions as before: Jacobi family, $P(\lambda) = (a-\lambda)(b-\lambda)(\lambda_c-\lambda)$, unit speed, in-quadrant standard L with corners at $\lambda_1 = 0, \alpha_1, \alpha_2$ and $\lambda_2 = \beta_1, \beta_2, \beta_3$, reflex corner at $(\alpha_1, \beta_2)$, $0 < \alpha_1 < \alpha_2 < b < \beta_1 < \beta_2 < \beta_3 < a$.

---

## Contents

1. [The phase space and its foliation](#1-the-phase-space-and-its-foliation)
2. [What makes a level critical](#2-what-makes-a-level-critical)
3. [The bifurcation list for the standard L](#3-the-bifurcation-list-for-the-standard-l)
4. [Type I — genus jump (corner-crossing)](#4-type-i--genus-jump-corner-crossing)
5. [Type II — edge swap (caustic reaches a wall)](#5-type-ii--edge-swap-caustic-reaches-a-wall)
6. [Type III — endpoint limits](#6-type-iii--endpoint-limits)
7. [The critical fiber as a spine plus cylinders](#7-the-critical-fiber-as-a-spine-plus-cylinders)
8. [The molecule](#8-the-molecule)
9. [Why it is not a Fomenko atom](#9-why-it-is-not-a-fomenko-atom)
10. [Recipe: build any critical fiber](#10-recipe-build-any-critical-fiber)
11. [Sanity checks](#11-sanity-checks)

---

## 1. The phase space and its foliation

Fix unit speed. The phase space is the unit tangent bundle over the table — a compact **3-manifold** $M^3$ with coordinates $(x, y, \theta)$, $\theta$ the velocity angle, minus the reflex-corner orbits. The single integral $\lambda_c$ foliates it:

$$\pi : M^3 \longrightarrow [\,0,\ \beta_3\,], \qquad \pi(x,y,\theta) = \lambda_c.$$

Each fiber $F_{\lambda_c} = \pi^{-1}(\lambda_c)$ is a **2-dimensional** level set. For regular $\lambda_c$ it is the closed genus-$g$ surface of the previous doc (with $g = 1 + n_{\text{acc}}$, $n_{\text{acc}}$ = accessible reflex corners). For critical $\lambda_c$ it degenerates. The whole "special-levels topology" question is: describe the finitely many critical fibers, and how the regular fibers on either side limit onto them.

The right global summary is the **Reeb graph** (Fomenko's *molecule*): collapse each connected fiber to a point. Because there is one integral and the fibers are already 2-dimensional, the quotient is a 1-dimensional graph over the $\lambda_c$-axis — §8.

---

## 2. What makes a level critical

Two independent mechanisms, and it pays to keep them apart because they were conflated in my first pass at this.

**(a) Region degeneration — a mechanism of the *surface*.** As $\lambda_c$ moves, the accessible flat region $\mathcal{R}_u(\lambda_c) \subset \mathcal{L}_u$ changes shape (the caustic cut sweeps across). At special $\lambda_c$ its *combinatorial* type changes: a corner appears/disappears, or an edge switches from wall to turning line. This is a statement about the flat surface $S(\lambda_c)$ regardless of the direction of flow, and it is what moves the genus. In flat-geometry terms $S(\lambda_c)$ hits the **boundary of its stratum**: a cylinder's width goes to $0$, or a modulus goes to $\infty$.

**(b) Saddle connections — a mechanism of the *dynamics*.** On a *fixed* surface, the $45°$ direction is sometimes a saddle-connection direction (a $45°$ geodesic runs from cone point to cone point). This governs periodic-vs-minimal (the IET data of the previous doc) and varies with $\lambda_c$ on a complicated arithmetic set. It does **not** move the genus.

The $\lambda_c$-bifurcations of the *foliation* are the type-(a) events. Type (b) is a finer structure living *within* each regular edge of the molecule; it does not create molecule vertices. Everything below is type (a) unless flagged.

A clean geometric restatement of a type-(a) critical level: **the caustic conic $\lambda = \lambda_*$ passes through a corner of the table, or is tangent to (coincides with) a wall.** Both are visible directly in the $(\lambda_1,\lambda_2)$ rectangle picture as the cut line hitting a vertex or an edge of $\mathcal{L}$.

---

## 3. The bifurcation list for the standard L

Sweep $\lambda_c$ upward and read the accessible region off the L rectangle. The corner is at $(\alpha_1, \beta_2)$.

| $\lambda_c$ | accessible region | genus | event at this $\lambda_c$ |
|---|---|---|---|
| $\to 0^+$ | thin annulus (whispering gallery) | 1 | endpoint (Type III) |
| $0 < \lambda_c < \alpha_1$ | rectangle (cut left of corner) | 1 | — |
| $\lambda_c = \alpha_1$ | cut **through the reflex corner** | — | **genus jump 1→2 (Type I)** |
| $\alpha_1 < \lambda_c < \alpha_2$ | L, base foot capped by turning line | 2 | — |
| $\lambda_c = \alpha_2$ | caustic reaches inner base wall | 2 | edge swap (Type II) |
| $\alpha_2 < \lambda_c < b$ | full L | 2 | — |
| $\lambda_c = b$ | caustic type switches ell.→hyp. | 2 | **non-event** (see below) |
| $b < \lambda_c < \beta_1$ | full L | 2 | — |
| $\lambda_c = \beta_1$ | caustic reaches base bottom wall | 2 | edge swap (Type II) |
| $\beta_1 < \lambda_c < \beta_2$ | L, bottom strip removed | 2 | — |
| $\lambda_c = \beta_2$ | cut **through the reflex corner** | — | **genus jump 2→1 (Type I)** |
| $\beta_2 < \lambda_c < \beta_3$ | rectangle (tall leg only) | 1 | — |
| $\lambda_c = \beta_3$ | region empties | — | endpoint (Type III) |

**The $\lambda_c = b$ non-event.** One might expect the separatrix to bifurcate the fiber, as it does for the plain ellipse. It does **not** here, because the in-quadrant L never touches the focal segment: $\lambda_1 \le \alpha_2 < b$ and $\lambda_2 \ge \beta_1 > b$, so both integration ranges stay bounded away from the double root of $P$ at $\lambda = b$. The lengths $u_1, u_2$ remain finite, the surface deforms smoothly, and only the *name* of the caustic (ellipse vs. hyperbola) changes. **The separatrix is a bifurcation only for tables that actually reach the focal segment.** This is worth checking for your specific table before assuming $\lambda_c = b$ is special.

So for the standard L there are exactly **two topology-changing levels**, $\alpha_1$ and $\beta_2$, both corner-crossings, mirror images of each other.

---

## 4. Type I — genus jump (corner-crossing)

The main event. Take $\lambda_c = \alpha_1$ (the other is identical under ell.↔hyp.).

**Flat picture.** Below $\alpha_1$ the cut $u_1 = u_1(\lambda_c)$ sits left of the corner image $u_1 = A_1$: the region is the rectangle $[0, u_1(\lambda_c)] \times [0, B_2]$. Just above, the cut clears the corner and the base "foot" $[A_1, u_1(\lambda_c)] \times [0, B_1]$ appears, of width $u_1(\lambda_c) - A_1 \to 0^+$. On the unfolded 4-sheet surface this foot and its three mirror copies assemble into a **flat cylinder** whose width shrinks to $0$ as $\lambda_c \to \alpha_1^+$.

**A vanishing cylinder is a pinch.** A cylinder of circumference $\ell$ and width $h \to 0$, collapsed, identifies its two boundary curves: a handle closes up. Concretely $S(\lambda_c)$ for $\lambda_c > \alpha_1$ is genus 2; as $\lambda_c \to \alpha_1^+$ the core curve of this cylinder (a non-separating simple closed curve) shrinks to a point. The critical fiber is therefore

$$F_{\alpha_1} \;=\; \big(\text{genus-2 surface with that cycle pinched to a point}\big) \;=\; \text{a torus with two points identified (a node).}$$

Its normalization is the genus-1 surface that governs $\lambda_c < \alpha_1$; its arithmetic genus is 2, matching $\lambda_c > \alpha_1$. So the single nodal surface sits correctly between the two regular sides:

```
   λc < α₁                 λc = α₁                 λc > α₁
  ┌─────────┐            ┌─────────┐             ┌─────────┐
  │         │            │         │             │         │
  │ torus   │   ──►      │  torus  │•  ──►        │ genus 2 │
  │  (g=1)  │            │  with a │              │  (g=2)  │
  │         │            │  node   │             └────┐    │
  └─────────┘            └─────────┘•                 │foot│
                          two points                  └────┘
                          glued (•=•)             handle opened
```

**Where the "few trajectories" live.** The pinch point is the image of the reflex corner, which at $\lambda_c = \alpha_1$ lies *exactly on the caustic turning locus* $\lambda_1 = \lambda_c$. The orbits that pass through the node are precisely the $45°$ flow-lines into the $6\pi$ cone point — the **banned corner trajectories** of the previous doc. At a regular level none of them are on the fiber in a singular way; at the critical level they become the entire singular set. There are finitely many up to the surface's symmetry (the corner has three prongs, so three critical directions meet the node), and they are the whole obstruction to smoothness. This is the precise sense in which "a very small amount of trajectories" builds the special topology.

**Surgery description.** Passing $\lambda_c$ through $\alpha_1$ performs the surgery

$$T^2 \ \xrightarrow{\ \text{attach a 1-handle}\ }\ \Sigma_2,$$

with the node as the intermediate. In Morse-theory-of-the-momentum-map language this is a critical point of the map $\lambda_c$ whose fiber is non-Bott (the Hessian degenerates along the corner orbit), which is exactly why no atom label applies (§9).

---

## 5. Type II — edge swap (caustic reaches a wall)

At $\lambda_c = \alpha_2$ (and $\beta_1$) the caustic cut reaches a wall of the base leg. Below, the foot's far edge is a **turning line** (caustic tangency); above, it is a **wall**. Both are mirror boundaries in the unfolding, so the surface deforms **continuously and the genus does not change**. This is a genus-preserving critical level: mild, but still a molecule vertex, because the *cylinder decomposition* of the surface can change (a band that was bounded by a turning line becomes bounded by a wall, altering which saddle connections exist).

In flat-surface terms $S(\lambda_c)$ does not leave its stratum; it passes through a wall–turning-line tangency, a codimension-1 wall in the moduli of the *marked* surface but not a genus change. If you only care about the homeomorphism type of the fiber you may **merge these vertices into the adjacent edges**; if you care about the periodic-band structure (the IET combinatorics) keep them.

A quick test to decide whether a given wall-reaching level matters for *your* purpose: compare the genus and the number of cylinders on both sides. Same genus, same cylinder count → cosmetic, merge it. Same genus, different cylinder count → keep it as a Type II vertex.

---

## 6. Type III — endpoint limits

**$\lambda_c \to 0^+$ (whispering gallery).** The accessible annulus pinches to the boundary curve; the fiber degenerates to the two boundary-tracing orbits (clockwise / counterclockwise). In the molecule this is a terminal **atom A** (a torus shrinking onto a circle) — here genuinely a Liouville atom, since near the boundary the dynamics is a rotation and the fiber is an honest torus collapsing.

**$\lambda_c \to \beta_3^-$ (region empties).** The tall leg's caustic cut reaches its outer wall; the accessible region shrinks to nothing. Terminal atom A again, on the tall-leg side.

These two are the only genuinely *elliptic* (center-type) ends; they cap the molecule.

---

## 7. The critical fiber as a spine plus cylinders

The general construction, valid at any type-(a) critical level, and the concrete realization of your idea.

Fix $\lambda_c = \lambda_*$ and its flat surface $S_*$ (possibly already degenerate). The **separatrix graph** $\Gamma$ is the union of all $45°$ saddle connections — geodesics in the flow direction joining cone points (here the reflex-corner image to itself or to the vanishing-cylinder boundary). Two facts:

1. $\Gamma$ is a **finite** graph. Saddle connections in a fixed direction on a translation surface are isolated and, below any length bound, finite; the ones bounding the degenerating cylinder are $O(1)$ in number (typically 2–3 for the L).
2. $S_* \setminus \Gamma$ is a finite disjoint union of **open flat cylinders** $C_1, \dots, C_m$ — the maximal families of parallel closed $45°$ orbits.

Then the critical fiber is assembled as

$$F_{\lambda_*} \;=\; \Gamma \ \cup\ \overline{C_1} \ \cup \cdots \cup\ \overline{C_m}, \qquad \text{each } \overline{C_j} \text{ glued to } \Gamma \text{ along its two boundary circles.}$$

This is a **ribbon graph thickened by annuli** — a completely combinatorial object. You read its topology from the graph: Euler characteristic $\chi(F) = \chi(\Gamma)$ (cylinders contribute $0$), and orientability plus the boundary-gluing give the homeomorphism type. The "small amount of trajectories" is exactly $\Gamma$; the cylinders are the bulk, and they carry no bifurcation information (they persist smoothly to both neighboring regular levels).

For the Type I node at $\lambda_c = \alpha_1$: $\Gamma$ is a single figure-eight (the corner orbit's two prongs closing up through the node), $m = 1$ big cylinder (the surviving rectangle) whose two ends both attach to the figure-eight — giving the pinched torus. Draw $\Gamma$, attach the annulus, and you have $F_{\alpha_1}$ with no further choices.

---

## 8. The molecule

Collapsing each fiber to a point turns the foliation into a labeled graph over $[0, \beta_3]$. For the standard L, with Type II vertices merged (homeomorphism-type molecule):

```
        α₁ (genus jump 1→2)        β₂ (genus jump 2→1)
         │                          │
  A ─────●──────────────────────────●───── A
  ▲      genus 2 all along here      ▲
  │                                  │
  λc→0⁺                            λc→β₃⁻
  (whispering                     (tall-leg
   gallery)                        empties)
```

A path with two elliptic ends (atoms A) and two interior **genus-jump vertices** (not atoms). The middle edge carries genus-2 fibers; the two outer edges carry tori. Keeping the Type II vertices adds two valence-2 marks on the middle edge recording the cylinder-count changes at $\alpha_2$ and $\beta_1$.

For a table with more reflex corners the molecule is longer and can branch (when a caustic disconnects the region into components, the graph splits into parallel edges over that $\lambda_c$-interval, one per component, each with its own genus $1 + n_{\text{comp}}$). The rule is mechanical: **vertices = corner-crossings and region-splits; edge labels = (component, genus); the two ends of each maximal genus-1 branch cap with atom A.**

---

## 9. Why it is not a Fomenko atom

Fomenko–Zieschang theory classifies singular fibers of integrable systems that are **non-degenerate** in the Bott sense: the integral is a Morse–Bott function on each energy level, so critical submanifolds are non-degenerate and the local fiber is a product of an atom's surface with circles. Its regular fibers are **always tori**.

The corner-crossing violates the hypothesis in the sharpest possible way: the regular fibers on one side are **genus 2**, which no atom produces. The critical point of $\lambda_c$ sits *at the cone point*, where the level function is not Morse–Bott — the degeneracy is exactly the $6\pi$ excess. So the classical machinery does not apply, and the honest local model is the **flat-geometry one**: a family of translation surfaces in $\mathcal{H}(2)$ meeting the stratum boundary $\partial\mathcal{H}(2)$ (where a cylinder collapses and the surface drops to $\mathcal{H}(0) = $ the torus).

The two genus-preserving events *are* classical:

| level | classical? | model |
|---|---|---|
| $\lambda_c \to 0^+$, $\to \beta_3^-$ | yes | atom A (center / elliptic) |
| $\lambda_c = \alpha_2, \beta_1$ (edge swap) | yes (if kept) | atom B–like (orientable saddle) if a cylinder splits |
| $\lambda_c = \alpha_1, \beta_2$ (corner) | **no** | $\mathcal{H}(2) \to \partial\mathcal{H}(2)$ cylinder collapse; genus $2 \leftrightarrow 1$ |

The clean way to say it: **the confocal square is Liouville-integrable and its molecule is made of atoms; the confocal L is pseudo-integrable and its molecule additionally contains genus-jump vertices that are stratum degenerations.** The reflex corner is exactly the ingredient that pushes you out of the atom classification.

(This is the setting of Vedyushkina–Fomenko *billiard books* and Dragović–Radnović's *pseudo-integrable confocal billiards*: the former realize Fomenko atoms by gluing confocal pieces along arcs; the latter analyze precisely the genus-jump degenerations that appear once reflex corners are allowed.)

---

## 10. Recipe: build any critical fiber

Per candidate level $\lambda_* $:

1. **Is it type (a)?** Check whether the caustic conic $\lambda = \lambda_*$ passes through a corner (→ genus jump if a *reflex* corner, region-split if it pinches a neck) or is tangent to / coincides with a wall (→ edge swap). If neither, $\lambda_*$ is regular — no vertex. Confirm the focal value $\lambda_* = b$ against §3's non-event caveat.
2. **Accessible region on both sides.** Compute $\mathcal{R}_u(\lambda_* \mp \varepsilon)$ and their genera $1 + n_{\text{acc}}$ per component. The pair $(g_-, g_+)$ names the vertex.
3. **Locate the vanishing/created cylinder.** It is the region cell that has zero measure exactly at $\lambda_*$ (the shrinking foot, the pinched neck). Its core curve's homology class is the cycle being pinched.
4. **Separatrix graph $\Gamma$.** Shoot $45°$ rays from each reflex-corner image; keep those that terminate at a cone point. These finitely many saddle connections are $\Gamma$.
5. **Assemble** $F_{\lambda_*} = \Gamma \cup (\text{closures of the surviving cylinders})$, gluing each cylinder along its boundary circles. Read $\chi = \chi(\Gamma)$ and orientability → homeomorphism type.
6. **Place it in the molecule** as a vertex between the $(g_-)$ and $(g_+)$ edges, atom-labeled only if $g_- = g_+ = 1$.

For the standard L this yields, with no free choices: atoms A at the two ends, nodal (pinched-torus) critical fibers at $\alpha_1$ and $\beta_2$, genus-2 fibers between them.

---

## 11. Sanity checks

1. **Arithmetic/geometric genus match.** At a Type I node, arithmetic genus (normalization genus + number of glued point-pairs) must equal the *higher* neighboring genus, and the normalization genus must equal the *lower*. For $\alpha_1$: normalization $= 1$, one node, arithmetic $= 2$. ✓
2. **Euler characteristic across the vertex.** $\chi(F_{\lambda_*}) = \chi(\Gamma)$ must equal $2 - 2g_{\text{lower}} - (\text{\# nodes})$ computed from the normalization. Mismatch means a saddle connection was missed in $\Gamma$.
3. **Finite $u$-lengths at $b$.** For the in-quadrant table, verify $u_1, u_2$ stay finite as $\lambda_c \to b$ (integration ranges bounded away from the double root). If they diverge, your table touches the focal segment and $b$ *is* a separatrix vertex — recount.
4. **Cylinder persistence.** The bulk cylinders of a critical fiber must appear, deformed but present, on both neighboring regular fibers. A cylinder that exists on only one side means you mislabeled a Type II edge-swap as regular.
5. **Molecule endpoints are atoms A.** Every maximal genus-1 branch must terminate in a center (whispering-gallery or empty-region limit). A dangling genus-1 edge means a missing endpoint.
6. **Symmetry.** For a table symmetric under $y \to -y$, the molecule must be symmetric under the corresponding involution of $\lambda_c$; the two corner-crossings $\alpha_1, \beta_2$ are exchanged by the ell.↔hyp. duality and must carry identical vertex data.

---

## Appendix: one-paragraph summary

The integral $\lambda_c$ foliates the unit-tangent 3-manifold into 2-dimensional fibers; regular fibers are the genus-$1{+}n_{\text{acc}}$ surfaces of the previous doc, and the finitely many critical fibers are where the accessible flat region changes combinatorial type — precisely when the caustic conic passes through a corner or reaches a wall. A **reflex** corner-crossing changes the genus, and its critical fiber is a **pinched (nodal) surface**: the higher-genus surface with a cylinder collapsed to a point, the point being the image of the $6\pi$ cone point, which at the critical level lies on the caustic turning locus. The only trajectories that are singular there are the finitely many $45°$ orbits into the reflex corner — the banned corner orbits — so the critical fiber is genuinely thin: a finite separatrix graph with regular cylinders glued on, buildable by hand. Genus-preserving events (caustic reaching a wall, or a whispering-gallery/empty-region limit) are classical Fomenko atoms; the genus-jumps are not atoms but stratum-boundary degenerations $\mathcal{H}(2) \to \mathcal{H}(0)$. Assembling the vertices over the $\lambda_c$-axis gives the molecule: for the standard L, two elliptic ends and two genus-jump nodes.
