# L-shaped confocal billiards: pseudo-integrability and the genus count

**What this is.** A continuation of the confocal-billiard doc for tables with **reflex corners** (interior angle $3\pi/2$). The second integral survives, but the invariant surfaces are no longer tori: a table with $n$ reflex corners produces level sets of genus

$$\boxed{\;g = 1 + n\;}$$

(per connected component, and only counting reflex corners the caustic leaves accessible — both refinements below). The transformation recipe changes accordingly: no global angle pair, but a cleaner object instead — a **flat polygon on which the flow has slope exactly $\pm 1$**.

**Sanity check of the formula.** $n = 0$ (plain confocal rectangle) gives $g = 1$, the torus of the previous doc. So the naive guess $g = n$ cannot be right: it would assign the torus genus 0.

Conventions as before: Jacobi family, $P(\lambda) = (a-\lambda)(b-\lambda)(\lambda_c-\lambda)$, $\lambda_1 \le b \le \lambda_2 \le a$, unit speed.

---

## Contents

1. [Where a 3π/2 corner comes from](#1-where-a-3π2-corner-comes-from)
2. [The standard L table](#2-the-standard-l-table)
3. [What survives and what breaks](#3-what-survives-and-what-breaks)
4. [The genus count](#4-the-genus-count)
5. [The flat model: slope-±1 flow](#5-the-flat-model-slope-1-flow)
6. [The translation surface](#6-the-translation-surface)
7. [Level-set zoo: genus as a function of λc](#7-level-set-zoo-genus-as-a-function-of-λc)
8. [Banned trajectories are the 3-prong separatrices](#8-banned-trajectories-are-the-3-prong-separatrices)
9. [What replaces the rotation number](#9-what-replaces-the-rotation-number)
10. [The transformation recipe](#10-the-transformation-recipe)
11. [Reference implementation](#11-reference-implementation)
12. [Drawing it](#12-drawing-it)
13. [Sanity checks](#13-sanity-checks)
14. [Pointers to the literature](#14-pointers-to-the-literature)

---

## 1. Where a 3π/2 corner comes from

Every corner of a confocal table is a crossing of an ellipse $\lambda_1 = \text{const}$ with a hyperbola $\lambda_2 = \text{const}$, and confocal quadrics meet **orthogonally**. So locally every corner is a right-angled crossing dividing the neighbourhood into four quadrants, and the only question is how many quadrants the table occupies:

- **1 quadrant** → interior angle $\pi/2$ — convex corner (the only kind in the previous doc);
- **3 quadrants** → interior angle $3\pi/2$ — reflex corner, the new ingredient;
- 2 quadrants would mean the boundary passes straight through — not a corner.

No other angles are possible while the walls stay confocal, which is what keeps the analysis exact: the billiard remains separable in the $\lambda$-chart, wall by wall, and integrability fails only *globally*, through the topology.

---

## 2. The standard L table

To keep the chart 1:1, place the table strictly inside one quadrant of the plane ($x > 0$, $y > 0$), away from both axes and the focal segment. Then $(\lambda_1, \lambda_2)$ has no folds, and the table can be specified directly as a region in the $\lambda$-chart. Take

$$\mathcal{L} \;=\; \big([0,\alpha_1] \times [\beta_1,\beta_3]\big) \;\cup\; \big([0,\alpha_2] \times [\beta_1,\beta_2]\big),$$

with $0 < \alpha_1 < \alpha_2 < b < \beta_1 < \beta_2 < \beta_3 < a$. In the plane this is a curvilinear L: outer ellipse wall $\lambda_1 = 0$, two inner ellipse walls $\lambda_1 = \alpha_1, \alpha_2$, three hyperbola walls $\lambda_2 = \beta_1, \beta_2, \beta_3$.

```
 λ₂
 β₃ ┌───────┐
    │       │
    │  tall │
    │  leg  │
 β₂ │       └──────────────┐          ← reflex corner at (α₁, β₂)
    │                      │
    │        base          │
 β₁ └──────────────────────┘
    0       α₁             α₂   λ₁
```

Six corners: five convex, one reflex at $(\lambda_1, \lambda_2) = (\alpha_1, \beta_2)$, where the domain covers the three local quadrants $\{\lambda_1 < \alpha_1\} \cup \{\lambda_2 < \beta_2\}$.

Tables that cross the axes work too — you reinstate the fold labels $\operatorname{sign}(x), \operatorname{sign}(y)$ from the previous doc and the genus formula is unchanged — but the in-quadrant table is the clean setting and the rest of this doc assumes it.

---

## 3. What survives and what breaks

**Survives:**

- The integral. $\lambda_c = b - (xv_y - yv_x)^2 + (a-b)v_y^2$ is conserved along segments and across reflection off *any* confocal wall — the proof is local and never sees the corner. Every trajectory is still tangent to its caustic.
- Separation. $p_{\lambda_i}^2 = (\lambda_c - \lambda_i)/\big(4(\lambda_i - a)(\lambda_i - b)\big)$, ovals on the same cubic $P$, walls flip $p \mapsto -p$, turning points at $\lambda_c$. All local machinery is intact.
- The 4-point fibre. Two tangent lines to the caustic per point, two orientations each; the level set is a 4:1 cover of the accessible region.

**Breaks:**

- The gluing. Four copies of the region glued edge-to-edge is only a torus when every corner closes up smoothly. A reflex corner does not (§4), so the level set is a higher-genus surface.
- Arnold–Liouville. The theorem requires compact connected level sets *on which the flow has no equilibria and which admit as many commuting periodic flows as degrees of freedom*; a genus-$\ge 2$ surface admits no nonvanishing vector field at all ($\chi \ne 0$), so no action-angle chart exists, even locally-in-$\lambda_c$. The system is **pseudo-integrable** in the sense of Richens–Berry: a full set of integrals, non-toral level sets.
- The rotation number. Replaced by interval-exchange data (§9).

---

## 4. The genus count

Unfold: four copies of the accessible region, one per sign pair $(\sigma_1, \sigma_2) = (\operatorname{sign} p_{\lambda_1}, \operatorname{sign} p_{\lambda_2})$, glued along walls and turning lines. The surface is flat (the $u$-coordinates of §5 make this literal), so all curvature sits at the corner images, and Gauss–Bonnet does the counting.

**Convex corner.** Four copies of $\pi/2$ meet: cone angle $4 \cdot \tfrac{\pi}{2} = 2\pi$. Smooth point. Contributes nothing — which is why the previous doc never had to think about corners.

**Reflex corner.** Four copies of $3\pi/2$: cone angle $4 \cdot \tfrac{3\pi}{2} = 6\pi$. A genuine cone point with angle excess $6\pi - 2\pi = 4\pi$.

**Gauss–Bonnet** for a closed flat surface with cone points:

$$\sum_{\text{cone points}} (\theta_{\text{cone}} - 2\pi) = 2\pi\,(2g - 2).$$

With $n$ reflex corners: $n \cdot 4\pi = 2\pi(2g-2)$, hence

$$g = n + 1.$$

### Where that identity comes from

The continuous Gauss–Bonnet theorem for a closed oriented surface $S$ with a Riemannian metric is

$$\int_S K \, dA \;=\; 2\pi\,\chi(S) \;=\; 2\pi\,(2 - 2g),$$

where $K$ is Gaussian curvature and $\chi = 2 - 2g$ is the Euler characteristic. Our unfolded surface is **flat** — it is glued from Euclidean pieces (the four copies of the region, which the $u$-coordinates of §5 turn into genuine Euclidean polygons), so $K \equiv 0$ on the smooth part. All the curvature that the theorem demands has to be concentrated at the finitely many singular points, as *delta functions* of curvature. The clean way to make "curvature as a delta function" precise is the **angle defect**: at an isolated point where the total cone angle is $\theta_{\text{cone}}$, the concentrated curvature is

$$\int_{\{p\}} K\, dA \;=\; 2\pi - \theta_{\text{cone}} \;\equiv\; \delta_p,$$

the *defect*. A smooth flat point has $\theta_{\text{cone}} = 2\pi$ and $\delta_p = 0$, contributing nothing; a cone point has $\delta_p \ne 0$. (Sign convention: positive curvature = angle *deficit*, like a cone tip where less than $2\pi$ fits; negative curvature = angle *excess*, our case.) Substituting into Gauss–Bonnet turns the integral into a finite sum,

$$\sum_{p} \big(2\pi - \theta_{\text{cone}}(p)\big) \;=\; 2\pi\,(2 - 2g),$$

which, multiplying both sides by $-1$, is exactly the boxed identity above with $\theta_{\text{cone}} - 2\pi$ (the excess) on the left.

### Why one can localize curvature at a cone point

To see that a $\theta_{\text{cone}}$-cone really carries $2\pi - \theta_{\text{cone}}$ units of curvature, smooth the tip and apply the theorem to a small geodesic disk $D_\varepsilon$ around it. The metric version with boundary reads $\int_{D_\varepsilon} K\,dA + \oint_{\partial D_\varepsilon} k_g\, ds = 2\pi$, where $k_g$ is geodesic curvature of the boundary circle. For a genuine cone the sides are straight (geodesic) and the boundary arc of radius $r$ has length $\theta_{\text{cone}}\, r$, so $\oint k_g\,ds = \theta_{\text{cone}}$ (the total turning of the boundary is the cone angle, not $2\pi$). Hence the interior integral is $2\pi - \theta_{\text{cone}}$ regardless of $\varepsilon$ — the curvature is genuinely concentrated at the tip, independent of how you smooth it.

### The count for confocal corners

Our cone angles are forced by the confocal geometry. Every corner is an orthogonal crossing, so it locally splits the plane into four right-angle ($\pi/2$) quadrants; the table occupies $m$ of them, giving interior angle $m\pi/2$. The 4-sheet unfolding places one copy of the corner in each sheet, so the **total** cone angle assembled around the corner image is

$$\theta_{\text{cone}} \;=\; 4 \times (\text{interior angle}) \;=\; 4 \cdot \frac{m\pi}{2} \;=\; 2m\pi,$$

and the defect is $\delta = 2\pi - 2m\pi = 2\pi(1-m)$.

| corner | quadrants $m$ | interior angle | $\theta_{\text{cone}}$ | defect $\delta$ |
|---|---|---|---|---|
| convex | 1 | $\pi/2$ | $2\pi$ | $0$ |
| reflex | 3 | $3\pi/2$ | $6\pi$ | $-4\pi$ |

Only reflex corners contribute, each with defect $-4\pi$ (excess $+4\pi$). Summing $n$ of them and inserting into $\sum_p \delta_p = 2\pi(2-2g)$:

$$n \cdot (-4\pi) = 2\pi(2 - 2g) \quad\Longleftrightarrow\quad -2n = 2 - 2g \quad\Longleftrightarrow\quad g = n + 1.$$

### The translation-surface refinement

There is a sharper statement worth knowing, because it tells you *what kind* of singularity each reflex corner is, not just how much curvature it carries. On a translation surface the flat structure comes from a holomorphic 1-form $\omega$ (the abelian differential $dz$ in local flat coordinates), and a cone point of angle $2\pi(k+1)$ is exactly a **zero of order $k$** of $\omega$. The Riemann–Roch / degree count for a genus-$g$ surface says the zeros of any abelian differential total

$$\sum_{\text{zeros}} k_i \;=\; 2g - 2.$$

Our reflex corners have $\theta_{\text{cone}} = 6\pi = 2\pi(2+1)$, hence $k = 2$ each — double zeros, the *3-pronged* singularities of §8. With $n$ of them, $\sum k_i = 2n = 2g-2$, recovering $g = n+1$ once more and identifying the stratum: the L-table ($n=1$) lives in $\mathcal{H}(2)$, a single double zero; a $T$ or $Z$ ($n=2$) lives in $\mathcal{H}(2,2)$; and so on. The three prongs at each corner are the $k+1 = 3$ outgoing horizontal (here, $45°$) separatrices, which is the analytic reason the banned corner orbit has no unique continuation.

| table | $n$ | $g$ |
|---|---|---|
| confocal rectangle | 0 | 1 (torus) |
| L | 1 | 2 |
| T or Z | 2 | 3 |
| plus/cross | 4 | 5 |

**Two refinements**, both inherited from the counting rule of the previous doc:

1. **Per component.** If the caustic disconnects the accessible region, each component unfolds separately: $g = 1 + n_{\text{comp}}$ with $n_{\text{comp}}$ the reflex corners inside *that* component. Total level set = disjoint union.
2. **Only accessible corners count.** A reflex corner shadowed by the caustic (inside the forbidden zone) contributes nothing. §7 works this out for the L.

One hypothesis behind the clean formula: the 4-sheet unfolding of each component must be connected. For any table with hyperbola walls both coordinates librate (walls or turning points at both ends), so all four sign pairs communicate and this holds automatically. The circulation splitting of Case A (plain ellipse) cannot occur here.

---

## 5. The flat model: slope-±1 flow

The right chart is the **un-normalized length coordinates** — the $w$-integrals of the previous doc *without* the $\pi/W$ normalization:

$$u_1(\lambda_1) = \int_0^{\lambda_1} \frac{d\lambda}{\sqrt{P(\lambda)}}, \qquad u_2(\lambda_2) = \int_{\beta_1}^{\lambda_2} \frac{d\lambda}{\sqrt{P(\lambda)}}.$$

Two facts make these the canonical choice.

**Rectangles map to rectangles.** Each $u_i$ depends only on $\lambda_i$, so any product region in $\lambda$ maps to a product region in $u$. The curvilinear L maps to a genuinely **flat L**:

$$\mathcal{L}_u = \big([0,A_1] \times [0,B_2]\big) \cup \big([0,A_2] \times [0,B_1]\big),$$

$$A_i = u_1(\alpha_i), \qquad B_1 = u_2(\beta_2), \qquad B_2 = u_2(\beta_3).$$

The moduli $A_1, A_2, B_1, B_2$ depend on $\lambda_c$ through $P$ — each level set gets its own flat L.

**The flow has slope exactly $\pm 1$.** From $\dot\lambda_i = \sigma_i\, 2\sqrt{P(\lambda_i)}/(\lambda_2 - \lambda_1)$,

$$\frac{du_1}{dt} = \frac{2\sigma_1}{\lambda_2 - \lambda_1}, \qquad \frac{du_2}{dt} = \frac{2\sigma_2}{\lambda_2 - \lambda_1} \qquad\Longrightarrow\qquad \frac{du_2}{du_1} = \frac{\sigma_2}{\sigma_1} = \pm 1.$$

The awkward $(\lambda_2 - \lambda_1)$ coupling that made the torus phases wiggly cancels in the *ratio*. After the time change $ds = 2\,dt/(\lambda_2 - \lambda_1)$ the motion is unit-speed straight lines at $45°$.

**So the entire billiard, at fixed $\lambda_c$, is isomorphic to the flat billiard in $\mathcal{L}_u$ with direction locked to the diagonals.** Walls and turning lines both act as mirrors ($u_i$ reverses); the four sign pairs $(\sigma_1,\sigma_2)$ are the four diagonal directions; corners map to corners with the same angles, so the reflex corner of the table is the reflex corner of $\mathcal{L}_u$. Everything hard about the confocal geometry has been absorbed into four numbers.

---

## 6. The translation surface

Unfold the flat L by the two reflections $u_1 \mapsto -u_1$ and $u_2 \mapsto -u_2$, one copy per sign pair, drawn as

$$X = \sigma_1 u_1, \qquad Y = \sigma_2 u_2.$$

The four copies tile a **cross-shaped 12-gon**

$$\big([-A_1,A_1] \times [-B_2,B_2]\big) \cup \big([-A_2,A_2] \times [-B_1,B_1]\big),$$

```
            ┌───────┐
            │  −,+  │  +,+           Y
   ┌────────┤       ├────────┐
   │        │       │        │
   │  −,+   │       │   +,+  │
   ├────────┼───────┼────────┤ ─ X
   │  −,−   │       │   +,−  │
   │        │       │        │
   └────────┤       ├────────┘
            │  −,−  │  +,−
            └───────┘
```

with **opposite parallel edges identified by translation** (right outer edge $X = A_2$ glues to $X = -A_2$, top $Y = B_2$ to $Y = -B_2$, and each step edge $X = \pm A_1$, $Y = \pm B_1$ to its mirror image). Reflecting twice is translating — that is why the identifications are translations and why the reflected flow becomes a single **translation flow** in the fixed direction $45°$.

This object is a **translation surface**: genus 2, one singular point (the four images of the reflex corner glue into one $6\pi$ cone point), i.e. the stratum $\mathcal{H}(2)$ — the same family as the L-shaped tables of Veech/McMullen fame. Every fact below is imported wholesale from that theory.

Consistency check on the genus without Gauss–Bonnet: Euler characteristic of the identified 12-gon comes out $\chi = -2$, so $g = 2$. Either computation, same answer.

---

## 7. Level-set zoo: genus as a function of λc

The caustic truncates $\mathcal{L}_u$ exactly as it truncated the rectangle: an elliptic caustic imposes $\lambda_1 \le \lambda_c$, i.e. $u_1 \le u_1(\lambda_c)$, a vertical cut; a hyperbolic caustic imposes $\lambda_2 \ge \lambda_c$, a horizontal cut from below. The cut edge is a turning line, which glues exactly like a wall — so only the *shape* of the truncated region matters, and the question is always: **does the reflex corner survive the cut?**

For the standard L (recall the reflex corner sits at $(\alpha_1, \beta_2)$, i.e. at $(A_1, B_1)$ in flat coordinates):

| $\lambda_c$ | accessible region | topology |
|---|---|---|
| $0 < \lambda_c < \alpha_1$ | vertical strip $[0, u_1(\lambda_c)] \times [0, B_2]$ — rectangle | **torus** |
| $\lambda_c = \alpha_1$ | caustic ellipse passes through the reflex corner | bifurcation |
| $\alpha_1 < \lambda_c < b$ | L-shape, right wall of the base possibly replaced by the caustic line (for $\lambda_c < \alpha_2$) | **genus 2** |
| $\lambda_c = b$ | separatrix, $u$-lengths diverge | degenerate |
| $b < \lambda_c \le \beta_1$ | whole L | **genus 2** |
| $\beta_1 < \lambda_c < \beta_2$ | L cut from below by the caustic hyperbola; reflex corner still inside | **genus 2** |
| $\lambda_c = \beta_2$ | caustic hyperbola passes through the reflex corner | bifurcation |
| $\beta_2 < \lambda_c < \beta_3$ | strip $[0, A_1] \times [u_2(\lambda_c), B_2]$ — rectangle (the base is entirely forbidden) | **torus** |
| $\lambda_c \ge \beta_3$ | empty | — |

So the genus is a *function of the level*, jumping $1 \to 2 \to 1$ as $\lambda_c$ sweeps through, with jumps exactly at the two values where the caustic conic passes through the reflex corner — a new kind of bifurcation, distinct from the separatrix ($\lambda_c = b$, where periods diverge) and from the boundary-tangency of the previous doc. At the corner-crossing values the periods stay finite; only the topology snaps.

For a table with several reflex corners, run the same analysis corner by corner: at each level, $g = 1 + \#\{\text{reflex corners in the accessible region}\}$ per component. Components arise the same way as before — a caustic arc can pinch the region in two (e.g. a hyperbolic caustic threading a U-shaped table), and then each side gets its own surface with its own count.

---

## 8. Banned trajectories are the 3-prong separatrices

Why the reflex corner forces a ban while convex corners never did:

**Convex corner** ($\pi/2$): the unfolded surface is *smooth* there (cone angle $2\pi$), so the corner orbit has a unique continuation — unfold, continue the straight line, fold back. In the table this is the familiar rule "reflect in both walls", i.e. return antiparallel. Nothing needs banning.

**Reflex corner** ($3\pi/2$): the cone angle is $6\pi$, and a cone point of angle $6\pi$ has, in each fixed direction, **three** incoming and **three** outgoing straight rays (a zero of order 2 of the abelian differential — a 3-pronged singularity). An orbit arriving at the corner has three equally valid continuations; no canonical choice exists, so the orbit is declared singular. These are exactly your "invalid trajectories".

They are few but structurally important:

- **Measure zero**: countably many saddle connections per surface; a random initial condition never hits them.
- **They organize everything else**: two orbits passing the corner on opposite sides separate at rate $O(1)$ — not exponentially (the flow has zero entropy), but enough to shred any would-be invariant circle. This splitting *is* the mechanism by which the torus fails.
- **Numerically**: detect proximity of the straight segment to the corner point and terminate/flag the trajectory, exactly as the corner-handling note of the previous doc — but here it is a matter of correctness, not merely hygiene, since the reflection logic would otherwise silently pick one of the three prongs.

---

## 9. What replaces the rotation number

On the torus the dynamics was a rigid rotation; here the first-return map to a transversal (say the segment $X = 0^+$ in the unfolded picture) is an **interval exchange transformation** (IET) — a piecewise translation of $[0, \text{length}]$ with finitely many pieces (four to five for $\mathcal{H}(2)$). The data replacing the single number $\rho$:

- the **moduli** $(A_1, A_2, B_1, B_2)$, functions of $\lambda_c$ — they fix the IET's interval lengths;
- the fixed direction ($45°$ — the direction never varies; only the surface does).

Consequences, imported from translation-surface theory:

- **Periodic ⇔ closed.** If the flow in the $45°$ direction is periodic (all orbits close), the level set decomposes into cylinders of parallel closed orbits — the analogue of a rational rotation number, but now the surface splits into *finitely many bands with different periods*, separated by saddle connections. Poncelet-style closure still exists, level by level.
- **Minimal but not uniquely ergodic is possible.** Unlike an irrational rotation, an IET can be minimal (every orbit dense) while carrying several ergodic measures. A single trajectory may then equidistribute with respect to none of them. Practically: two long trajectories on the *same* level set can paint visibly different densities. This is not a bug in your code.
- **Zero entropy, weak mixing at most.** The flow is never chaotic in the exponential sense; the complexity is combinatorial.
- **Arithmetic sensitivity.** Whether a given $\lambda_c$ gives periodic, minimal-uniquely-ergodic, or minimal-non-uniquely-ergodic behaviour depends on the arithmetic of the moduli, which vary transcendentally with $\lambda_c$. Expect all behaviours to occur as $\lambda_c$ sweeps its range, interleaved in a complicated set. (When the moduli hit the special proportions of Veech/McMullen surfaces, the dichotomy is clean: every direction is either periodic or uniquely ergodic; generic $\lambda_c$ has no reason to be special.)

A practical substitute for the rotation-number check of the previous doc: for each $\lambda_c$, record the itinerary of wall labels, or measure the fraction of time spent in the tall leg vs. the base — for periodic levels these lock to rationals, and their variation with $\lambda_c$ makes a nice devil's-staircase-like plot.

---

## 10. The transformation recipe

Per sample $(x, y, v_x, v_y)$, at fixed table $(\alpha_1, \alpha_2, \beta_1, \beta_2, \beta_3)$:

1. **Coordinates.** $(\lambda_1, \lambda_2)$ from the confocal quadratic (sign-safe form). In-quadrant table ⇒ no fold labels needed.
2. **Level.** $\lambda_c$ from the velocity formula (normalize $v$ first).
3. **Classify the level** using §7: forbidden / torus / genus-2, and the accessible flat region $\mathcal{R}_u \subseteq \mathcal{L}_u$ with its moduli. If torus: fall back to the previous doc's normalized angles, done.
4. **Sheet.** $(\sigma_1, \sigma_2) = (\operatorname{sign}\dot\lambda_1, \operatorname{sign}\dot\lambda_2)$ from the implicit-differentiation formulas; carry-and-flip between samples for robustness, as before.
5. **Flat position.** $u_1(\lambda_1)$, $u_2(\lambda_2)$ by the precomputed cumulative integrals (substitution $\lambda = \lambda_* \mp s^2$ at any endpoint that is a root of $P$; endpoints that are walls are regular).
6. **Output** either
   - **table picture**: $(u_1, u_2)$ in $\mathcal{R}_u$ plus the direction $(\sigma_1, \sigma_2)$ (draw as color), or
   - **surface picture**: $(X, Y) = (\sigma_1 u_1, \sigma_2 u_2)$ in the cross-shaped polygon, edges identified by translation.

There is deliberately **no step normalizing to $[0, 2\pi)$**: on a genus-2 surface no such global pair of circles exists, and the honest coordinates are the flat ones. The `torus_index` of the previous API generalizes to `(component, sheet)`.

---

## 11. Reference implementation

Builds on the `confocal`, `caustic`, `lam_dots` and the $s$-substitution idea from the previous doc; restated here so this file is self-contained.

```python
import numpy as np


# ---------- basics (as in the torus doc) ----------

def confocal(x, y, a, b):
    p = a + b - x*x - y*y
    q = a*b - b*x*x - a*y*y
    disc = max(p*p - 4.0*q, 0.0)
    big = 0.5*(p + np.sqrt(disc)) if p >= 0 else 0.5*(p - np.sqrt(disc))
    small = q/big if big != 0.0 else 0.0
    return (min(small, big), max(small, big))


def caustic_level(x, y, vx, vy, a, b):
    n = np.hypot(vx, vy); vx, vy = vx/n, vy/n
    L = x*vy - y*vx
    return b - L*L + (a - b)*vy*vy


def lam_dots(x, y, vx, vy, lam1, lam2, a, b):
    d = lam1 - lam2
    return (2.0*(x*(b - lam1)*vx + y*(a - lam1)*vy)/d,
            -2.0*(x*(b - lam2)*vx + y*(a - lam2)*vy)/d)


# ---------- cumulative length u(lambda) ----------

class ULength:
    """u(lam) = int_lo^lam dl/sqrt(|P|) on [lo, hi], P = (a-l)(b-l)(lc-l).

    Integrates on a grid uniform in s = sqrt(r - l), where r is the root of P
    nearest the interval (above hi or at hi). That kills the inverse-sqrt
    endpoint singularity when hi IS a root (turning line), and is harmless
    when it isn't (wall). If instead the root sits at/below lo (hyperbolic
    turning at the lower end), pass root_side='lo'.
    """

    def __init__(self, lo, hi, a, b, lc, root_side='hi', knots=800):
        self.lo, self.hi = lo, hi
        roots = np.array([a, b, lc])
        if root_side == 'hi':
            r = roots[roots >= hi - 1e-14].min()
            s_lo, s_hi = np.sqrt(r - hi), np.sqrt(r - lo)
            s = np.linspace(s_lo, s_hi, knots)
            lam = r - s*s
        else:
            r = roots[roots <= lo + 1e-14].max()
            s_lo, s_hi = np.sqrt(hi - r), np.sqrt(lo - r)
            s = np.linspace(s_lo, s_hi, knots)
            lam = r + s*s
        other = np.abs(np.prod([q - lam for q in roots if q != r], axis=0))
        f = 2.0/np.sqrt(np.maximum(other, 1e-300))
        g = np.concatenate([[0.0],
                            np.cumsum(0.5*(f[1:] + f[:-1])*np.abs(np.diff(s)))])
        # g is cumulative from the r-side end; orient so u(lo)=0, u(hi)=total
        if root_side == 'hi':
            self._s, self._g, self._flip = s, g, True
        else:
            self._s, self._g, self._flip = s, g, False
        self.total = float(g[-1])
        self._r, self._side = r, root_side

    def u(self, lam):
        lam = min(max(lam, self.lo), self.hi)
        s = np.sqrt(abs(self._r - lam))
        val = float(np.interp(s, self._s, self._g))
        return (self.total - val) if self._flip else val


# ---------- level classification for the standard L ----------

def classify_level(lc, a, b, tab, sep_eps=1e-9):
    """tab = (alpha1, alpha2, beta1, beta2, beta3).  Returns a dict with
    kind in {'forbidden','torus','genus2'} and the flat moduli."""
    a1, a2, b1, b2, b3 = tab
    if abs(lc - b) < sep_eps:
        return {'kind': 'separatrix'}

    if lc < b:                      # elliptic caustic: lam1 <= lc
        if lc <= a1:
            U1 = ULength(0.0, lc, a, b, lc)            # turning at hi
            U2 = ULength(b1, b3, a, b, lc, 'hi')       # walls both ends
            return {'kind': 'torus', 'U1': U1, 'U2': U2,
                    'shape': [(U1.total, U2.total)]}
        hi1 = min(lc, a2)
        U1 = ULength(0.0, hi1, a, b, lc)
        U2 = ULength(b1, b3, a, b, lc, 'hi')
        A1 = ULength(0.0, a1, a, b, lc).total
        B1 = ULength(b1, b2, a, b, lc).total
        return {'kind': 'genus2', 'U1': U1, 'U2': U2,
                'A1': A1, 'A2': U1.total, 'B1': B1, 'B2': U2.total}

    # hyperbolic caustic: lam2 >= lc
    if lc >= b3:
        return {'kind': 'forbidden'}
    if lc >= b2:                    # base forbidden -> rectangle
        U1 = ULength(0.0, a1, a, b, lc)
        U2 = ULength(lc, b3, a, b, lc, 'lo')           # turning at lo
        return {'kind': 'torus', 'U1': U1, 'U2': U2,
                'shape': [(U1.total, U2.total)]}
    lo2 = max(lc, b1)
    U1 = ULength(0.0, a2, a, b, lc)
    U2 = ULength(lo2, b3, a, b, lc, 'lo' if lc > b1 else 'hi')
    A1 = ULength(0.0, a1, a, b, lc).total
    B1 = ULength(lo2, b2, a, b, lc, 'lo' if lc > b1 else 'hi').total
    return {'kind': 'genus2', 'U1': U1, 'U2': U2,
            'A1': A1, 'A2': U1.total, 'B1': B1, 'B2': U2.total}


# ---------- per-sample mapping ----------

def to_flat(x, y, vx, vy, a, b, tab, level=None):
    """-> (u1, u2, sheet, level_dict).  Draw (u1,u2) colored by sheet, or
    (sheet[0]*u1, sheet[1]*u2) in the unfolded cross polygon."""
    n = np.hypot(vx, vy); vx, vy = vx/n, vy/n
    lam1, lam2 = confocal(x, y, a, b)
    lc = caustic_level(x, y, vx, vy, a, b)
    if level is None:
        level = classify_level(lc, a, b, tab)
    if level['kind'] in ('forbidden', 'separatrix'):
        raise ValueError(level['kind'])
    d1, d2 = lam_dots(x, y, vx, vy, lam1, lam2, a, b)
    sheet = (1 if d1 >= 0 else -1, 1 if d2 >= 0 else -1)
    return level['U1'].u(lam1), level['U2'].u(lam2), sheet, level
```

Same caveats as before: precompute `classify_level` once per trajectory; carry the sheet signs between samples instead of trusting `d1, d2` near turning lines; guard $|\lambda_c - b|$; and terminate trajectories that approach the reflex corner — in flat coordinates that is proximity of the segment to the point $(A_1, B_1)$, which is a cleaner test than anything in the $(x,y)$ plane.

---

## 12. Drawing it

**No donut.** A genus-2 surface has no flat embedding in $\mathbb{R}^3$ and no canonical round picture; a pretzel embedding is possible but arbitrary and hides the flat structure that carries all the dynamics. Use the flat pictures:

1. **The flat L, four colors.** Plot $(u_1, u_2)$ inside $\mathcal{R}_u$, colored by sheet. Trajectories are $45°$ segments bouncing off the sides. This is the direct generalization of the flat square of the torus doc — indeed for $n = 0$ it *is* that picture, up to the $\pi/W$ rescaling of each axis.
2. **The unfolded cross.** Plot $(\sigma_1 u_1, \sigma_2 u_2)$ in the 12-gon of §6 and mark the edge identifications (matching arrows on opposite parallel edges). Trajectories are parallel $45°$ segments that exit one edge and re-enter at the identified point of the partner edge. Mark the four images of the reflex corner — visually, orbits shear around these points, which is the genus made visible.
3. **Level stack.** The most informative global picture: a column of flat-L thumbnails indexed by $\lambda_c$, showing the region deform, the corner get swallowed at $\lambda_c = \alpha_1$ and $\beta_2$, and the topology jump torus → genus 2 → torus.

For torus levels, revert to the previous doc's normalized $(\theta_1, \theta_2)$ so the two docs' pictures agree on the overlap.

---

## 13. Sanity checks

1. **λc constancy** across reflections off *all five* wall types — catches a wall that is not actually confocal.
2. **Slope lock.** In $(u_1, u_2)$, every sampled segment must have slope $\pm 1$ to quadrature accuracy. This single check validates the coordinates, the cubic, and the integrals at once; it is the analogue (and strengthening) of the old rotation-number test.
3. **Cone angle by simulation.** Launch a tight pencil of parallel orbits aimed just left and just right of the reflex corner; they must land $O(1)$ apart afterwards (3-prong splitting). Around any convex corner the pencil must reassemble.
4. **Genus-1 regression.** For $\lambda_c < \alpha_1$ or $\lambda_c \in (\beta_2, \beta_3)$, the output must reproduce the torus doc's picture exactly after rescaling axes by $\pi/W_i$.
5. **Bifurcation values.** The topology must switch precisely at $\lambda_c = \alpha_1$ and $\lambda_c = \beta_2$ — the caustic through the corner. Off-by-one errors in the classification table show up here.
6. **Periodic levels.** Hunt (by bisection in $\lambda_c$) for a level where the itinerary repeats; on it, verify the surface decomposes into parallel period bands — cylinders — rather than a single closed curve. A closed *single* curve on a genus-2 level means the unfolding is wrong.

---

## 14. Pointers to the literature

Informal, for orientation:

- **Richens & Berry (1981)** — coined *pseudo-integrable* for polygonal billiards with rational angles whose invariant surfaces have genus $\ge 2$; the $g = n+1$ phenomenon for L-shaped flat tables is already here.
- **Zemlyakov & Katok (1975)** — the unfolding construction (reflections → translation surface) used in §6.
- **Dragović & Radnović (2014–)** — *pseudo-integrable billiards within confocal conics*: precisely this setting (confocal tables with $3\pi/2$ corners), including the genus computations, the reduction to interval exchanges, and the arithmetic phenomena of §9.
- **Veech (1989), McMullen (2003)** — L-shaped translation surfaces in $\mathcal{H}(2)$, the optimal dynamical dichotomies when the moduli are special. Your surfaces live in the same stratum with $\lambda_c$-dependent moduli.

---

## Appendix: the one-paragraph summary

A reflex corner leaves every *local* structure of the confocal billiard intact — integral, caustics, separation — and breaks exactly one *global* thing: the four-sheet unfolding stops closing up smoothly, acquiring a $6\pi$ cone point per reflex corner. Gauss–Bonnet then forces genus $1 + n$, killing action-angle variables and rotation numbers ("pseudo-integrable"). The correct replacement for the torus chart is the pair of unnormalized lengths $u_i = \int d\lambda_i/\sqrt{P}$, in which each level set becomes a flat L-shaped billiard with the direction frozen at $45°$ — a translation surface in $\mathcal{H}(2)$ whose moduli vary with $\lambda_c$ — and the dynamics becomes an interval exchange. The banned corner trajectories are the 3-prong separatrices of the cone point, and the genus itself is a function of the level: the caustic can shadow the corner and restore the torus, with bifurcations exactly when the caustic conic passes through the reflex corner.
