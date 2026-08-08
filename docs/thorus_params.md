# Liouville tori for confocal billiards

**What this is.** A complete recipe for taking a billiard trajectory sample `(x, y, vx, vy)` in a confocal domain and mapping it to a point `(θ₁, θ₂, torus_index)` on a Liouville torus, plus the geometry needed to know *how many* tori there are and where the count changes.

**Conventions.** Jacobi normalization of the confocal family, unit speed, and the cubic $P(\lambda)$ defined in §3. Everything reduces to that one cubic.

---

## Contents

1. [The confocal family](#1-the-confocal-family)
2. [Coordinates are roots of the family equation](#2-coordinates-are-roots-of-the-family-equation)
3. [The second integral and the elliptic curve](#3-the-second-integral-and-the-elliptic-curve)
4. [Turning points vs. folds — the thing that bites](#4-turning-points-vs-folds--the-thing-that-bites)
5. [Ovals, per regime](#5-ovals-per-regime)
6. [Building the angles](#6-building-the-angles)
7. [Counting tori](#7-counting-tori)
8. [Case A — full ellipse, elliptic caustic](#8-case-a--full-ellipse-elliptic-caustic)
9. [Case B — full ellipse, hyperbolic caustic](#9-case-b--full-ellipse-hyperbolic-caustic)
10. [Case C — confocal curvilinear square](#10-case-c--confocal-curvilinear-square)
11. [Bifurcations](#11-bifurcations)
12. [Drawing it](#12-drawing-it)
13. [Reference implementation](#13-reference-implementation)
14. [Numerical pitfalls](#14-numerical-pitfalls)
15. [Sanity checks](#15-sanity-checks)
16. [Summary tables](#16-summary-tables)

---

## 1. The confocal family

$$(b - \lambda)\,x^2 + (a - \lambda)\,y^2 = (a - \lambda)(b - \lambda), \qquad \lambda \le a,$$

with fixed $\infty > a > b > 0$. Foci at $(\pm c, 0)$, $c^2 = a - b$.

| $\lambda$ | curve |
|---|---|
| $\lambda < 0$ | ellipses larger than the base ellipse |
| $\lambda = 0$ | base ellipse $x^2/a + y^2/b = 1$ |
| $0 < \lambda < b$ | nested ellipses, shrinking as $\lambda$ grows |
| $\lambda = b$ | degenerate: the focal segment $[-c, c]$, plus the rays $\lvert x\rvert \ge c$ |
| $b < \lambda < a$ | hyperbolas, confocal with the ellipses |
| $\lambda = a$ | degenerate: the line $x = 0$ (inside the table, the chord $\lvert y\rvert \le \sqrt b$) |

**Why this normalization.** In the more familiar $x^2/A^2 + y^2/B^2 = 1$ convention you carry $\cosh\mu$ and $\cos\nu$ separately and end up with two different-looking separated equations. Here $A^2 \to a$, $B^2 \to b$, and

$$c^2\cosh^2\mu = a - \lambda_1, \qquad c^2\cos^2\nu = a - \lambda_2,$$

so the coordinates, the caustic parameter, and the branch points all live on the *same* $\lambda$-axis. One formula covers both degrees of freedom, and the underlying elliptic curve becomes visible. This is why the algebro-geometric literature (Moser; Veselov; Dragović–Radnović) writes the family this way.

---

## 2. Coordinates are roots of the family equation

Through every point of the plane pass exactly two members of the family — one ellipse, one hyperbola — and **their parameters are the coordinates**. Substitute $(x,y)$ into the family and solve for $\lambda$:

$$\boxed{\;\lambda^2 - (a + b - x^2 - y^2)\,\lambda + (ab - b x^2 - a y^2) = 0\;}$$

Roots $\lambda_1 \le b \le \lambda_2 \le a$. Inside the base ellipse, $\lambda_1 \in [0, b]$ and $\lambda_2 \in [b, a]$.

Handy identities (Vieta, then a two-line expansion):

$$\lambda_1 + \lambda_2 = a + b - x^2 - y^2, \qquad \lambda_1\lambda_2 = ab - bx^2 - ay^2,$$

$$x^2 = \frac{(a-\lambda_1)(a-\lambda_2)}{a-b}, \qquad y^2 = \frac{(b-\lambda_1)(\lambda_2-b)}{a-b}.$$

The second pair is the inverse map (up to the signs of $x$ and $y$) and makes an excellent round-trip test.

**Reading the degenerate quadrics.** The two degenerate members split neatly between the roots — which is exactly why the chart has folds:

- $\lambda_1 = b$ ⟺ on the focal segment ($y = 0$, $\lvert x\rvert \le c$)
- $\lambda_2 = b$ ⟺ on the outer rays ($y = 0$, $\lvert x\rvert \ge c$) — outside the table
- $\lambda_2 = a$ ⟺ on the $y$-axis ($x = 0$)
- $\lambda_1 = \lambda_2 = b$ ⟺ at a focus

**Angle reconstruction.** The chart is 2:1 or 4:1 over the plane, so you often need the unfolded angle:

$$\nu = \operatorname{atan2}\!\left(\frac{y}{\sqrt{b-\lambda_1}},\ \frac{x}{\sqrt{a-\lambda_1}}\right).$$

Exact, not an approximation: $x = \sqrt{a-\lambda_1}\cos\nu$, $y = \sqrt{b-\lambda_1}\sin\nu$ parametrizes the ellipse $\lambda = \lambda_1$.

---

## 3. The second integral and the elliptic curve

### The caustic parameter

For unit speed,

$$\boxed{\;\lambda_c = b - (x v_y - y v_x)^2 + (a-b)\,v_y^2\;}$$

is conserved along straight segments **and** across reflections off any member of the family. It is the classical Joachimsthal integral, equal to $b - K$ where $K$ is the product of angular momenta about the two foci.

For non-unit speed, divide the velocity-dependent part by $\lvert v\rvert^2$:

$$\lambda_c = b - \frac{(xv_y - yv_x)^2 - (a-b)v_y^2}{\lvert v\rvert^2}.$$

Simplest to normalize $v$ at the top of the mapping function and forget about it.

| $\lambda_c$ | caustic |
|---|---|
| $0 < \lambda_c < b$ | ellipse — orbit stays *outside* it |
| $\lambda_c = b$ | separatrix — orbit passes through a focus each bounce |
| $b < \lambda_c < a$ | hyperbola (both branches) — orbit stays *between* the branches |

### Separation

Hamilton–Jacobi separates, and both degrees of freedom obey the **same** equation:

$$p_{\lambda_i}^2 = \frac{\lambda_c - \lambda_i}{4\,(\lambda_i - a)(\lambda_i - b)}, \qquad i = 1, 2.$$

Sign check: for $\lambda_1 < b$ both numerator and denominator are positive when $\lambda_1 < \lambda_c$; for $\lambda_2 > b$ both are negative when $\lambda_2 > \lambda_c$. So the accessible set is $\lambda_1 \le \lambda_c$ (elliptic case) or $\lambda_2 \ge \lambda_c$ (hyperbolic case) — precisely "the orbit never crosses its caustic".

### The curve

Both ovals live on one cubic:

$$\boxed{\;w^2 = P(\lambda) := (a - \lambda)(b - \lambda)(\lambda_c - \lambda)\;}$$

$P \ge 0$ exactly on $(-\infty,\ \min(b,\lambda_c)] \cup [\max(b,\lambda_c),\ a]$ — one interval per degree of freedom. This is the Jacobi–Moser picture: the Liouville torus is the real part of the Jacobian of this curve, and the differentials $d\lambda/\sqrt{P}$ in §6 are its Abelian differentials.

### Velocities in the chart

Implicit differentiation of the quadratic — no square roots, no branch ambiguity:

$$\dot\lambda_1 = \frac{2\big[x(b - \lambda_1)v_x + y(a - \lambda_1)v_y\big]}{\lambda_1 - \lambda_2}, \qquad \dot\lambda_2 = \frac{2\big[x(b - \lambda_2)v_x + y(a - \lambda_2)v_y\big]}{\lambda_2 - \lambda_1}.$$

Equivalently $\dot\lambda_i = \pm\,2\sqrt{P(\lambda_i)}\,/\,\lvert\lambda_2 - \lambda_1\rvert$, and $\operatorname{sign}\dot\lambda_i = \operatorname{sign} p_{\lambda_i}$.

### Time is not the phase

$$dt = \frac{\lambda_2 - \lambda_1}{2\sqrt{P(\lambda_1)}}\,d\lambda_1.$$

The $(\lambda_2 - \lambda_1)$ factor couples the two degrees of freedom, so the phases built in §6 give a **wiggly** winding curve on the torus, not a straight line at constant speed. Same topology, same rotation number, same closed-vs-dense dichotomy — but not the linearized flow. See the end of §6 if you need the real thing.

---

## 4. Turning points vs. folds — the thing that bites

Every oval endpoint is a zero of $P$ or a wall, and $\dot\lambda \to 0$ at all of them. But they are **three different things**, and conflating them is the single most common way to get a torus that looks plausible and has the wrong rotation number.

| kind | where | what the orbit does | what the coordinate does |
|---|---|---|---|
| **Turning point** | $\lambda_i = \lambda_c$ | tangency to the caustic; genuinely reverses | $p_{\lambda_i} = 0$, smooth sign flip |
| **Wall** | $\lambda_i$ = boundary value | reflects; velocity reverses discontinuously | $p_{\lambda_i} \mapsto -p_{\lambda_i}$, and $P > 0$ there |
| **Fold** | $\lambda_1 = b$ or $\lambda_2 = a$ | *nothing* — sails through at full speed | chart is 2:1; must be unfolded |

The folds are the sneaky ones. At $\lambda_1 = b$,

$$b - \lambda_1 = \frac{(a-b)\,y^2}{\lambda_2 - b},$$

quadratic in $y$ — so $\lambda_1$ has a smooth *maximum* at $y = 0$ while the trajectory crosses the focal segment obliviously. Likewise $a - \lambda_2 = (a-b)x^2/(a-\lambda_1)$ at $x = 0$.

**Unfold with the geometric sign**, never with $\operatorname{sign}\dot\lambda$:

- fold at $\lambda_1 = b$ → label with $\operatorname{sign}(y)$
- fold at $\lambda_2 = a$ → label with $\operatorname{sign}(x)$

Skip the unfolding and you get a $\mathbb{Z}_2$ (or $\mathbb{Z}_2^2$) quotient of the true torus — still a torus, still a winding curve, but every rotation number is off by a factor of 2.

> **Useful fact.** In the hyperbolic-caustic regime, $\lambda_2 \ge \lambda_c > b$ means the orbit never touches $y = 0$ *outside* the foci. So every $y = 0$ crossing is a focal-segment crossing, and $\operatorname{sign}(y)$ is a clean label with no special cases.

---

## 5. Ovals, per regime

For the **full elliptic billiard** (boundary $\lambda_1 = 0$):

| regime | coord | interval | lower end | upper end |
|---|---|---|---|---|
| elliptic $0 < \lambda_c < b$ | $\lambda_1$ | $[0,\ \lambda_c]$ | wall | turning point |
| | $\lambda_2$ | $[b,\ a]$ | fold ($y=0$, outside foci) | fold ($x=0$) → use $\nu$ |
| hyperbolic $b < \lambda_c < a$ | $\lambda_1$ | $[0,\ b]$ | wall | fold ($\operatorname{sign} y$) |
| | $\lambda_2$ | $[\lambda_c,\ a]$ | turning point | fold ($\operatorname{sign} x$) |

**$\lambda_1$ never reaches $b$ in the elliptic regime.** The oval is $[0, \lambda_c]$, not $[0, b]$ — the value $b$ is simply not on this branch of the curve, since $p_{\lambda_1}^2 < 0$ beyond $\lambda_c$. This is what fixes $W_1$ in §6.

In the elliptic regime $\lambda_2$ traverses $b \to a \to b \to a \to b$ as $\nu$ runs once around $[0, 2\pi)$: it is a **double cover of the angle**, its endpoints are folds rather than turning points, and $p_{\lambda_2}$ *diverges* there while the motion is perfectly smooth. Use $\nu$ from §2; it is the honest circle.

For curvilinear polygons (§10) the intervals get cut by walls instead of running to the folds — but the classification of each endpoint into wall / turning point / fold is unchanged.

---

## 6. Building the angles

### Simple libration (a wall or turning point at each end, no fold)

With oval $[\lambda_{\min}, \lambda_{\max}]$,

$$w(\lambda) = \int_{\lambda_{\min}}^{\lambda} \frac{d\lambda'}{\sqrt{P(\lambda')}}, \qquad W = w(\lambda_{\max}),$$

$$\theta = \begin{cases}\pi\,w(\lambda)/W & \dot\lambda > 0\\[2pt] 2\pi - \pi\,w(\lambda)/W & \dot\lambda < 0\end{cases}$$

The increasing branch runs $\theta: 0 \to \pi$, the decreasing branch $\pi \to 2\pi \equiv 0$. Continuous, monotone along the flow, covers $[0,2\pi)$ exactly once per libration.

**Use $\pi w/W$, not $2\pi w/W$.** The factor $2\pi$ would run the increasing branch alone over the whole circle and then double-cover on the way back. The wrap at $\theta = \pi$ sits exactly at the turning point, where the branch label changes — that is the correct place for it, not a discontinuity.

### Libration through a fold

When the oval has a fold interior to the physical motion, build a monotone parameter $u$ covering the full unfolded range, then apply the same rule.

**Fold at $\lambda_1 = b$**, label $\sigma = \operatorname{sign}(y)$, wall at $\lambda_{\text{wall}}$:

$$W_1 = \int_{\lambda_{\text{wall}}}^{b}\frac{d\lambda}{\sqrt P}, \qquad u_1 = W_1 + \sigma\!\!\int_{\lambda_1}^{b}\!\frac{d\lambda}{\sqrt P} \in [0,\,2W_1],$$

$$\theta_1 = \frac{\pi u_1}{2W_1} \ \text{ if } \dot u_1 > 0, \quad \text{else } 2\pi - \frac{\pi u_1}{2W_1}, \qquad \operatorname{sign}\dot u_1 = -\sigma\operatorname{sign}\dot\lambda_1.$$

**Fold at $\lambda_2 = a$**, label $\tau = \operatorname{sign}(x)$, turning point or wall at $\lambda_{\text{end}}$:

$$W_2 = \int_{\lambda_{\text{end}}}^{a}\frac{d\lambda}{\sqrt P}, \qquad u_2 = W_2 - \tau\!\!\int_{\lambda_2}^{a}\!\frac{d\lambda}{\sqrt P} \in [0,\,2W_2],$$

$$\theta_2 = \frac{\pi u_2}{2W_2} \ \text{ if } \dot u_2 > 0, \quad \text{else } 2\pi - \frac{\pi u_2}{2W_2}, \qquad \operatorname{sign}\dot u_2 = +\tau\operatorname{sign}\dot\lambda_2.$$

These are the phases of $\mu \in [-\mu_0, \mu_0]$ and $\nu \in [\nu_c, \pi - \nu_c]$ rewritten in $\lambda$, with $\cos^2\nu_c = (a - \lambda_c)/(a-b)$.

### Circulation

Only $\lambda_2$ in the full-ellipse elliptic-caustic case circulates. Don't build a phase from $P$ — use $\theta_2 = \nu$ from §2 directly. The uniformized version is $\int d\nu/\sqrt{\lambda_c - b + (a-b)\cos^2\nu}$ normalized by its full period, but $\nu$ itself is already a perfectly good circle coordinate for drawing.

### Determining the branch

$\operatorname{sign}(\dot\lambda_i)$ from §3 is correct, but recomputing it per sample is fragile in two spots: at turning points it is genuinely $0$ and the sign is noise; near the foci $\lambda_2 - \lambda_1 \to 0$ makes the formula $0/0$.

**Robust approach:** compute the sign once at $t = 0$, carry it, and flip when $\lambda_i$ reaches an oval endpoint (detect by loss of monotonicity, or by $P(\lambda_i) < \epsilon$). With dense sampling along chords this is trivial. With sparse sampling, fall back to the formula and accept occasional glitches at the turns. For folded coordinates track $\operatorname{sign}(\dot u)$, not $\operatorname{sign}(\dot\lambda)$.

### If you need true action-angle variables

The phases above are a diffeomorphism of the torus, not the linearizing chart. For straight-line flow at constant speed you need

$$I_j = \frac{1}{2\pi}\oint p_{\lambda_j}\,d\lambda_j, \qquad \theta_j = \frac{\partial S}{\partial I_j},$$

which means inverting $\partial(I_1, I_2)/\partial(H, \lambda_c)$ numerically. Worth it for measuring frequencies or KAM-style perturbation work; unnecessary for drawing, since the rotation number is chart-independent.

---

## 7. Counting tori

> **Rule.** The number of Liouville tori over a given $\lambda_c$ is the number of connected components of the level set in phase space. Compute it as
>
> **(components of the accessible region in configuration space) × (sheets per component)**,
>
> with the sheets glued along walls, turning points, and folds.

Two facts make this mechanical:

**Fact 1 — the fibre is always four points.** Through any interior point of the accessible region pass exactly two tangent lines to the caustic conic, each carrying two orientations. The level set is a 4:1 cover of the accessible region, always. What differs between cases is only *how the four split*.

**Fact 2 — two gluing shapes cover everything.** An annulus with two boundary circles, doubled along both, is a torus:
$$S^1 \times [0,1] \ \cup_{S^1\times\partial[0,1]}\ S^1\times[0,1] \;=\; T^2.$$
A disk with four sides, taken in four copies and glued edge-to-edge, is also a torus — the standard unfolding of a rectangular billiard.

Each case below is one of these two shapes.

---

## 8. Case A — full ellipse, elliptic caustic

**Accessible region.** The annulus between caustic and boundary. Connected — the disk inside the caustic is never visited, but that does not disconnect anything. Only **one** region, so the second torus has to come from somewhere else.

**The split is in velocity, not position.** Sort the four directions at a point by the angular momentum about the centre,

$$L = x v_y - y v_x.$$

You get two and two — one orientation from *each* tangent line in each group. $L$ is a legitimate component label because it can only change sign by vanishing, and $L = 0$ forces

$$\lambda_c = b + (a-b)v_y^2 \ \ge\ b,$$

which is the hyperbolic regime. So for $\lambda_c < b$, $L$ never vanishes anywhere on the level set. It is not a new integral — just a locally constant function, i.e. exactly a component label. And $\operatorname{sign}\dot\nu = \operatorname{sign} L$, so it is the circulation sense.

**Result: 2 tori**, exchanged by time reversal $v \mapsto -v$ (equivalently $\nu \mapsto -\nu$), each a 2:1 cover of the annulus.

**The gluing.** The two sheets are $p_{\lambda_1} > 0$ (outward, caustic → wall) and $p_{\lambda_1} < 0$ (inward, wall → caustic). Glue:

- along the **caustic** circle, where $p_{\lambda_1} = 0$ and outward becomes inward smoothly;
- along the **boundary** circle, where reflection sends $p_{\lambda_1} \mapsto -p_{\lambda_1}$ — discontinuous in the plane, continuous in phase space.

Two annuli glued along *both* pairs of boundary circles is a torus. (One pair alone gives a sphere with two holes; the second gluing closes the handle.) Concretely the $[0,1]$ factor is $\lambda_1 \in [0,\lambda_c]$ doubling into the $\theta_1$ circle, and the $S^1$ factor is $\nu = \theta_2$, untouched.

**Coordinates.**
- $\theta_1$: simple libration on $[0, \lambda_c]$, $W_1 = \int_0^{\lambda_c} d\lambda/\sqrt P$. Wall at $0$, turning point at $\lambda_c$. Advances $2\pi$ per bounce.
- $\theta_2 = \nu$. Advances $2\pi\rho$ per bounce, $\rho$ = Poncelet rotation number.
- `torus_index` $= \operatorname{sign}(L) \in \{0, 1\}$.

---

## 9. Case B — full ellipse, hyperbolic caustic

**Accessible region.** The part of the table between the two branches of the caustic hyperbola — a curvilinear quadrilateral bounded by two hyperbola arcs and two arcs of the boundary ellipse. Connected, a topological disk.

**No circulation to split on.** The orbit crosses the focal segment, $L$ is not conserved, and all four sign-branches $(\pm p_{\lambda_1}, \pm p_{\lambda_2})$ are reachable from one another. The $y > 0$ and $y < 0$ halves are identified through the focal segment by the deck transformation $(\mu, \nu) \sim (-\mu, -\nu)$ at $\mu = 0$.

**Result: 1 torus**, a 4:1 cover of the quadrilateral — four copies glued edge-to-edge.

**Coordinates.** Both angles are folded librations (§6):
- $\theta_1$: fold at $\lambda_1 = b$ ($\operatorname{sign} y$), wall at $\lambda_1 = 0$. Advances $\pi$ per bounce.
- $\theta_2$: fold at $\lambda_2 = a$ ($\operatorname{sign} x$), turning point at $\lambda_2 = \lambda_c$.
- `torus_index` $= 0$ always.

---

## 10. Case C — confocal curvilinear square

**The table.** Bounded above and below by two arcs of the ellipse $\lambda_1 = 0$, left and right by the two branches of a confocal hyperbola $\lambda_2 = \beta$, $\beta \in (b, a)$. In the $(\lambda_1, \lambda_2)$ chart this table is literally a rectangle. All four corners are right angles, since confocal quadrics meet orthogonally — which is exactly what keeps the billiard integrable.

**An elliptic caustic splits it.** Take $\lambda_c \in (0, b)$. Tangency forces $\lambda_1 \le \lambda_c$, so the orbit stays *outside* the caustic ellipse. The caustic contributes two arcs to the table — an upper and a lower one, each running from the left wall to the right — cutting it into three pieces:

1. **upper region** — between the top boundary arc and the upper caustic arc: accessible
2. **middle lens** — inside the caustic: forbidden, no orbit with this $\lambda_c$ enters
3. **lower region** — between the lower caustic arc and the bottom boundary arc: accessible

The two accessible regions are disjoint, so $\operatorname{sign}(y)$ is locally constant on the level set.

**Each region carries one torus.** The upper region is a disk with four sides, and every side does the same job:

| side | value | effect |
|---|---|---|
| outer ellipse arc | $\lambda_1 = 0$ | wall: $p_{\lambda_1} \mapsto -p_{\lambda_1}$ |
| caustic arc | $\lambda_1 = \lambda_c$ | turning point: $p_{\lambda_1} = 0$ |
| left hyperbola arc | $\lambda_2 = \beta$, $x < 0$ | wall: $p_{\lambda_2} \mapsto -p_{\lambda_2}$ |
| right hyperbola arc | $\lambda_2 = \beta$, $x > 0$ | wall: $p_{\lambda_2} \mapsto -p_{\lambda_2}$ |

Four sides, four sign-branches, four copies glued — the rectangular unfolding again. Same for the lower region.

**Result: 2 tori**, four sheets each.

**Same count as Case A, entirely different reason.** Worth internalizing:

| | Case A (plain ellipse) | Case C (confocal square) |
|---|---|---|
| accessible set | one annulus | two disks |
| tori | 2 | 2 |
| sheets per torus | 2 | 4 |
| $\theta_2$ | circulates ($\nu$) | librates (wall to wall) |
| label | $\operatorname{sign}(L)$ | $\operatorname{sign}(y)$ |
| swaps the two tori | time reversal $v \to -v$ | reflection $y \to -y$ |
| preserves each torus | reflection | time reversal |

That last pair is the substantive difference. In Case C both coordinates librate, so $v \mapsto -v$ merely flips both sign-branches and lands on the *same* torus; the mirror symmetry is what pairs them. If the table is not symmetric about the $x$-axis (different hyperbola parameters above and below), the two tori are not even congruent — but there are still two.

**Coordinates.**
- $\theta_1$: simple libration on $[0, \lambda_c]$, $W_1 = \int_0^{\lambda_c} d\lambda/\sqrt P$.
- $\theta_2$: folded libration, fold at $\lambda_2 = a$ ($\operatorname{sign} x$), walls at $\lambda_2 = \beta$; so $W_2 = \int_\beta^a d\lambda/\sqrt P$ — the only change from Case B is the lower limit.
- `torus_index` $= \operatorname{sign}(y)$.

---

## 11. Bifurcations

Where the count or the structure changes rather than deforming. Branch on these explicitly in code.

**$\lambda_c \to b$ (separatrix).** The caustic ellipse flattens onto the focal segment. In Case A the two circulation tori merge; in Case C the forbidden lens collapses and the upper and lower regions touch along the focal segment. Either way **2 → 1**. The period diverges logarithmically ($W \to \infty$), so detect $\lvert\lambda_c - b\rvert < \epsilon$ and either skip or render the separatrix figure-eight explicitly.

**$\lambda_c \to 0^+$.** The caustic swells toward the boundary and the annulus (Case A) or the accessible strips (Case C) pinch. The tori degenerate to the whispering-gallery circles — the boundary run clockwise and counterclockwise. Two of them, visibly distinct: a good limiting check on the `torus_index` logic.

**Caustic stops fitting (Case C variant).** If the inner boundary is itself a confocal ellipse $\lambda_1 = \alpha$ rather than the focal region, then caustics with $\lambda_c > \alpha$ lie entirely inside the hole: no splitting, orbits bounce between the two ellipse arcs with $\lambda_1 \in [0, \alpha]$, **one torus**. Crossing $\lambda_c = \alpha$ downward is where the caustic emerges into the table and the region splits, **1 → 2**. This is a boundary-tangency bifurcation, not a separatrix — periods stay finite through it, unlike at $\lambda_c = b$. Handle it separately.

**Corners (curvilinear polygons only).** Orbits hitting a corner exactly are singular and belong to no torus. Measure zero, but numerically you should detect proximity to the corner points and drop those trajectories rather than let the reflection logic guess.

---

## 12. Drawing it

Standard embedding:

$$\big((R + r\cos\theta_1)\cos\theta_2,\ (R + r\cos\theta_1)\sin\theta_2,\ r\sin\theta_1\big).$$

Here $\theta_2$ is the azimuth about the $z$-axis — the **major / toroidal** angle — and $\theta_1$ parametrizes the tube cross-section (it modulates both distance from the axis and height) — the **minor / poloidal** angle. So assigning $\theta_1$ = librating per-chord phase and $\theta_2$ = circulating per-lap phase puts each where you want it: one chord = one trip around the tube, one lap of the ellipse = one trip around the donut. Keep the same assignment in Cases B and C so the picture varies continuously as $\lambda_c$ crosses $b$.

**The flat square is usually more informative.** Scatter $(\theta_1, \theta_2)$ on $[0,2\pi)^2$: the orbit fills a winding line whose slope is the rotation number $\rho$. Rational $\rho$ → closed Poncelet polygon → closed curve; irrational → densely fills the torus.

**Offset the two tori in space** when there are two — they occupy identical $(\theta_1,\theta_2)$ ranges and would otherwise overdraw.

**Sample along chords, not only at bounces.** A pipeline that records only reflection events gives a discrete point set — the billiard *map* orbit — not the flow curve. Sample each chord at $N$ equal arclength steps. (The map picture is still a useful cross-check: plot $(s, \sin\theta)$ at reflections and each torus collapses to an invariant curve.)

---

## 13. Reference implementation

Not production code — no caching strategy, fixed-order quadrature — but every formula above appears in it exactly once.

```python
import numpy as np
from numpy.polynomial.legendre import leggauss


# ---------- coordinates ----------

def confocal(x, y, a, b):
    """(lam1, lam2) with lam1 <= b <= lam2. Sign-safe near the foci."""
    p = a + b - x*x - y*y            # lam1 + lam2
    q = a*b - b*x*x - a*y*y          # lam1 * lam2
    disc = max(p*p - 4.0*q, 0.0)     # equals (lam2 - lam1)**2
    big = 0.5*(p + np.sqrt(disc)) if p >= 0 else 0.5*(p - np.sqrt(disc))
    small = q/big if big != 0.0 else 0.0
    return (small, big) if small <= big else (big, small)


def caustic(x, y, vx, vy, a, b):
    """lambda_c. Normalizes v, so any speed is fine."""
    n = np.hypot(vx, vy)
    vx, vy = vx/n, vy/n
    L = x*vy - y*vx
    return b - L*L + (a - b)*vy*vy


def nu_angle(x, y, lam1, a, b):
    """Unfolded angle around the ellipse lam = lam1. Exact."""
    return np.arctan2(y/np.sqrt(max(b - lam1, 1e-300)),
                      x/np.sqrt(max(a - lam1, 1e-300)))


def lam_dots(x, y, vx, vy, lam1, lam2, a, b):
    d = lam1 - lam2
    d1 = 2.0*(x*(b - lam1)*vx + y*(a - lam1)*vy)/d
    d2 = -2.0*(x*(b - lam2)*vx + y*(a - lam2)*vy)/d
    return d1, d2


# ---------- the cubic and its integrals ----------

def _P_without(lam, roots, skip):
    """P(lam) / (roots[skip] - lam): the product of the other two factors."""
    out = np.ones_like(np.asarray(lam, dtype=float))
    for j, r in enumerate(roots):
        if j != skip:
            out = out*(r - lam)
    return out


class Libration:
    """Monotone phase w(lam)/W on an oval whose UPPER end is a simple root of P.

    Substituting lam = hi - s**2 removes the inverse-square-root endpoint
    singularity exactly, so the grid is uniform in s -- which is where the
    integrand is smooth. Precompute once per trajectory (lambda_c is fixed),
    then every sample is an O(log n) lookup instead of a quadrature call.
    """

    def __init__(self, lo, hi, roots, knots=512):
        self.lo, self.hi, self.roots = lo, hi, roots
        self.k = int(np.argmin([abs(r - hi) for r in roots]))
        smax = np.sqrt(max(hi - lo, 0.0))
        s = np.linspace(0.0, smax, knots)
        lam = hi - s*s
        f = 2.0/np.sqrt(np.abs(_P_without(lam, roots, self.k)))
        g = np.concatenate([[0.0],
                            np.cumsum(0.5*(f[1:] + f[:-1])*np.diff(s))])
        self.s_grid, self.g_grid = s, g
        self.W = float(g[-1])          # the full half-period

    def w(self, lam):
        """int_lo^lam dlam'/sqrt(P)."""
        s = np.sqrt(max(self.hi - lam, 0.0))
        return self.W - float(np.interp(s, self.s_grid, self.g_grid))

    def tail(self, lam):
        """int_lam^hi dlam'/sqrt(P). This is what the fold formulas want."""
        return self.W - self.w(lam)

    def theta(self, lam, increasing):
        t = np.pi*self.w(lam)/self.W
        return t if increasing else 2.0*np.pi - t


def gauss_integral(lo, hi, roots, root_at_hi=True, n=256):
    """Standalone quadrature, if you'd rather not build a Libration."""
    k = int(np.argmin([abs(r - (hi if root_at_hi else lo)) for r in roots]))
    smax = np.sqrt(max(hi - lo, 0.0))
    s, w = leggauss(n)
    s = 0.5*smax*(s + 1.0)
    w = 0.5*smax*w
    lam = (hi - s*s) if root_at_hi else (lo + s*s)
    return float(np.sum(w*2.0/np.sqrt(np.abs(_P_without(lam, roots, k)))))


# ---------- top-level mapping ----------

def to_torus(x, y, vx, vy, a, b, cache=None, lam_wall=0.0, beta=None,
             sep_eps=1e-9):
    """Map one phase-space sample to (theta1, theta2, torus_index).

    lam_wall : lambda_1 of the outer boundary (0 for the base ellipse)
    beta     : lambda_2 of the hyperbola walls; None means the full ellipse
    cache    : dict reused across samples of the SAME trajectory
    """
    cache = {} if cache is None else cache
    n = np.hypot(vx, vy)
    vx, vy = vx/n, vy/n
    lam1, lam2 = confocal(x, y, a, b)
    lc = caustic(x, y, vx, vy, a, b)
    d1, d2 = lam_dots(x, y, vx, vy, lam1, lam2, a, b)
    roots = (a, b, lc)

    if abs(lc - b) < sep_eps:
        raise ValueError("separatrix: the torus degenerates here")

    # ---------------- elliptic caustic ----------------
    if lc < b:
        if "L1" not in cache:
            cache["L1"] = Libration(lam_wall, lc, roots)
        th1 = cache["L1"].theta(lam1, d1 > 0)

        if beta is None:                      # Case A -- nu circulates
            th2 = nu_angle(x, y, lam1, a, b) % (2.0*np.pi)
            idx = 0 if (x*vy - y*vx) > 0 else 1
        else:                                 # Case C -- folded libration
            if "L2" not in cache:
                cache["L2"] = Libration(beta, a, roots)
            L2 = cache["L2"]
            tau = 1.0 if x >= 0 else -1.0
            u2 = L2.W - tau*L2.tail(lam2)
            frac = np.pi*u2/(2.0*L2.W)
            th2 = frac if (tau*d2) > 0 else 2.0*np.pi - frac
            idx = 0 if y >= 0 else 1
        return th1, th2, idx

    # ---------------- hyperbolic caustic: both coords fold ----------------
    if "H1" not in cache:
        cache["H1"] = Libration(lam_wall, b, roots)
        cache["H2"] = Libration(lc, a, roots)
    H1, H2 = cache["H1"], cache["H2"]

    sig = 1.0 if y >= 0 else -1.0
    u1 = H1.W + sig*H1.tail(lam1)
    f1 = np.pi*u1/(2.0*H1.W)
    th1 = f1 if (-sig*d1) > 0 else 2.0*np.pi - f1

    tau = 1.0 if x >= 0 else -1.0
    u2 = H2.W - tau*H2.tail(lam2)
    f2 = np.pi*u2/(2.0*H2.W)
    th2 = f2 if (tau*d2) > 0 else 2.0*np.pi - f2

    return th1, th2, 0
```

**Interface note.** Keep `PhasePoint` as `(x, y, θ)` if you like and add a `velocity()` accessor returning $(\cos\theta, \sin\theta)$ — nothing is lost and unit speed is enforced for free. Have the mapping take `(x, y, vx, vy)` so it stays testable in isolation with synthetic inputs.

---

## 14. Numerical pitfalls

**Endpoint singularities.** Every $\int d\lambda/\sqrt P$ has an inverse-square-root endpoint. The substitution $\lambda = \lambda_* \mp s^2$ at a simple root $\lambda_*$ removes it *exactly*. For $\lambda_* = \lambda_c$:

$$\int \frac{d\lambda}{\sqrt{P(\lambda)}} = \int \frac{2\,ds}{\sqrt{(a - \lambda_c + s^2)(b - \lambda_c + s^2)}},$$

and analogously at $\lambda_* = b$ or $a$ with the remaining two factors. Plain Gauss–Legendre then converges fast. Adaptive quadrature on the raw integrand does not, and will quietly return garbage near the turning point.

**Precompute per torus, not per sample.** $\lambda_c$ is fixed along a trajectory, so build the monotone spline once and make each sample a lookup. Never call `quad` in the inner loop.

**Foci.** The quadratic's discriminant is $(\lambda_2 - \lambda_1)^2$, vanishing at $(\pm c, 0)$. Use the sign-safe form: compute the larger-magnitude root by the standard formula, then get the other from the product $q/\lambda_{\text{big}}$. Hyperbolic orbits with $\lambda_c$ near $b$ pass close to this.

**Separatrix.** $\lvert\lambda_c - b\rvert < \epsilon$ ⟹ $W \to \infty$. Detect and branch.

**Square-root guards.** $\sqrt{b - \lambda_1}$ and $\sqrt{a - \lambda_2}$ can see tiny negative arguments from round-off at the folds. Clamp at zero rather than letting `atan2` receive a NaN.

---

## 15. Sanity checks

Cheap, and each catches a distinct class of bug.

1. **Round-trip the coordinates.** Reconstruct $x^2, y^2$ from $(\lambda_1, \lambda_2)$ using §2 and compare. Catches root-ordering and sign errors.
2. **$\lambda_c$ is constant** to machine precision across a reflection. Catches a wrong reflection law or a non-confocal wall.
3. **$\Delta\theta_1$ per bounce** $= 2\pi$ (elliptic caustic) or $\pi$ (hyperbolic). Catches a wrong $W_1$ — in particular, using $[0,b]$ where $[0,\lambda_c]$ was needed.
4. **Rotation number.** The slope of the winding line on the flat square should equal the Poncelet rotation number. Measuring $2\rho$ means you used $\lambda_2$ where $\nu$ was needed, or skipped an unfolding.
5. **Poncelet closure.** Feed a caustic known to close after $n$ bounces and confirm the torus curve closes after exactly $n$.
6. **Component labels.** Verify $L$ never changes sign for $\lambda_c < b$ (Case A), and that $y$ never changes sign in Case C. A single flip means the caustic is misclassified or you crossed the separatrix.
7. **Whispering-gallery limit.** As $\lambda_c \to 0^+$ the two tori should collapse onto two distinct circles.

---

## 16. Summary tables

### Endpoint classification

| endpoint | root of $P$? | orbit does | coordinate does | handle by |
|---|---|---|---|---|
| turning point $\lambda_c$ | yes | reverses | $p = 0$, smooth flip | branch sign |
| wall | no | reflects | $p \mapsto -p$ | branch sign |
| fold $\lambda_1 = b$ | yes | passes through | 2:1 | $\operatorname{sign}(y)$ |
| fold $\lambda_2 = a$ | yes | passes through | 2:1 | $\operatorname{sign}(x)$ |

### Case summary

| case | region | tori | sheets/torus | $\theta_1$ | $\theta_2$ | index |
|---|---|---|---|---|---|---|
| A: ellipse, $\lambda_c < b$ | annulus | 2 | 2 | libration $[0,\lambda_c]$ | $\nu$ (circulates) | $\operatorname{sign} L$ |
| B: ellipse, $\lambda_c > b$ | quadrilateral | 1 | 4 | folded $[0,b]$ | folded $[\lambda_c,a]$ | $0$ |
| C: confocal square, $\lambda_c < b$ | two disks | 2 | 4 | libration $[0,\lambda_c]$ | folded $[\beta,a]$ | $\operatorname{sign} y$ |

### Everything in one place

$$c^2 = a - b, \qquad P(\lambda) = (a-\lambda)(b-\lambda)(\lambda_c-\lambda),$$

$$\lambda_c = b - (xv_y - yv_x)^2 + (a-b)v_y^2, \qquad p_{\lambda_i}^2 = \frac{\lambda_c - \lambda_i}{4(\lambda_i-a)(\lambda_i-b)},$$

$$\lambda^2 - (a + b - x^2 - y^2)\lambda + (ab - bx^2 - ay^2) = 0, \qquad w(\lambda) = \int \frac{d\lambda}{\sqrt{P(\lambda)}}.$$
