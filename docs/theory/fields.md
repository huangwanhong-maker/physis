# Fields: Their Bases, Representations and Laws

*Status: draft (2026-10-29), [decision 0012](../../development/decisions/0012-fields-as-sections.md).* The theory of physis's `fields` module, of the geometric complexes and cochain fields of `discretisation`, and of the first balance law of `continuum`. Every statement that the code relies on is either proved here or cited, and every one that can be checked numerically is held by a check, named beside it.

## 1. Fields as Sections

### 1.1 The Four Axes of a Field

**Definition 1.1 (field).** A *field* is a section $\varphi\in\Gamma(E)$ of a bundle $E\to B$ over a *base* $B$, whose fibre carries a physical dimension $\delta$ (its values lie in $E_b\otimes\mathbb R_\delta$). It is held for computation by a *representation* and comes from evidence of a *provenance*.

The four are kept apart because each changes while the others stay:

| Axis | What it says | Examples |
|---|---|---|
| base | where the field is defined | a region of space; a surface; a curve; a worldline; time alone; phase space; a complex that is itself the world |
| fibre | what values it takes, and in what unit | a scalar; a vector on named axes; a symmetric tensor; a $k$-form; a label in a finite set; a measure; a connection |
| representation | how it is held | a formula; a spectral expansion; samples with a stated reconstruction; a cochain on a cellulation; a composite of parts; a uniform value |
| provenance | what kind of evidence its values are | physical; archive; model; assimilated; synthetic |

*Rationale.* Scalar, vector, tensor, form and label fields, on regions, surfaces, complexes or phase space, are unified by the base and the fibre, not by a class per kind. The representation is the axis along which computation differs: a spectral field samples a sphere by rows, and a cochain pairs with chains and has no value at a point without a reconstruction. Keeping it apart from the base is what lets a Regge geometry (edge lengths on a complex: a metric, not a cochain) and a 1-cochain (a circulation) share a base and differ in fibre and representation.

### 1.2 Evaluations and Their Qualities

**Definition 1.2 (evaluation).** A query is a pair $(t,x)$ of an instant and a point of the base. Its *evaluation* is either a value with quality *valid*, possibly with a *note* that qualifies it, or one of three refusals, each with its reason:

- *outside*: $x\notin B$, so the field has no value there (the World Magnetic Model below 1 km under the ellipsoid; salinity in the air);
- *missing*: $x\in B$, but the data that would give the value are absent;
- *invalid*: $x\in B$, but the value is ill-defined there (a direction where the vector vanishes).

A refusal is never a zero. A field that is physically defined only on its material region (formalism §3.2) has no value elsewhere. A zero would be a claim about the world.

### 1.3 Capabilities

What a representation can do is a set of *capabilities*. Each is a C++20 concept on the representation's type, not a method that every field must have:

- evaluation at a point, and in batches;
- an exact spatial derivative;
- pairing with chains, a coboundary, a Hodge star;
- sampling a sphere faster than point by point;
- transfer to another cellulation;
- committed versions in time.

`capabilities_of<F>()` reads them off the type, a type may state more (`extra_capabilities`), and `described(f)` reports them at run time. A composition keeps a capability only when it provides it. A sum reads a sphere by rows when both its parts do, and it has no derivative of its own.

## 2. Bases and Their Domains

**Definition 2.1 (domain).** A *domain* $D$ of points $P$ is a map $q\mapsto\mathrm{Cover}(q)$ to "covered" or (quality, reason). Its *support* is $\{x: \mathrm{Cover}(t,x)\text{ covered}\}$ at each $t$.

The domains written are:

- `Everywhere`;
- `Shell` (inner $\le\lvert x-c\rvert\le$ outer, each side with its reason);
- `Box`;
- `Interval` (on a line);
- `During` (an interval of time);
- `AllOf`, the intersection, whose refusal is that of the first part that refuses.

Reasons are views of static text, so that evaluating costs no allocation. A domain that names a body is the body's domain's to write: the WMM's region above WGS 84, the Moon's lowest ground.

## 3. Representations

### 3.1 Formulas, and Their Exact Derivatives

An analytic field is a formula $f(t,x)$ on a domain, written once as a template over the number type of $x$.

**Proposition 3.1.** Evaluated on dual numbers $x+\epsilon_i e_i$ ($\epsilon_i\epsilon_j=0$), a formula composed of arithmetic and differentiable elementary functions gives $f(x)+\sum_i\epsilon_i\,\partial_i f(x)$. Its derivatives are exact up to the rounding of their evaluation.

*Proof.* Forward-mode automatic differentiation (Griewank and Walther 2008, ch. 3): the dual numbers form the ring $\mathbb R[\epsilon_1,\ldots,\epsilon_n]/(\epsilon_i\epsilon_j)$, each elementary $g$ is extended by $g(a+\sum b_i\epsilon_i)=g(a)+g'(a)\sum b_i\epsilon_i$, and the chain rule makes composition a ring homomorphism on first-order jets. $\square$

The Jacobian $J$ of a vector field gives $\nabla\cdot v=\operatorname{tr}J$ and $(\nabla\times v)_i=\varepsilon_{ijk}J_{kj}$. For a field without an exact derivative, central differences of step $h$ along each axis give $\partial_i f$ with truncation error $\tfrac{h^2}{6}\partial_i^3f$ and rounding error of order $\varepsilon\lvert f\rvert/h$. Every check that uses them derives its limit from the first, as below (§3.2), and keeps the step large enough for the second to be negligible.

### 3.2 Potential Fields Outside a Sphere

The spectral representation is physis `analysis`'s internal expansion (decision 0010): $B=-\nabla V$ with
$$
V(r,\psi,\lambda)=a\sum_{n=1}^{N}\Big(\frac ar\Big)^{n+1}\sum_{m=0}^{n}(g_n^m\cos m\lambda+h_n^m\sin m\lambda)\,\breve P_n^m(\sin\psi),
$$
on the shell $r\ge a$ outside its sources. `analysis` gives the components north $X$, east $Y$ and down $Z$ at planetocentric latitude $\psi$ and longitude $\lambda$. On the frame's axes,
$$
B=X\hat n+Y\hat e-Z\hat u,\qquad \hat u=(\cos\psi\cos\lambda,\cos\psi\sin\lambda,\sin\psi),\quad \hat e=(-\sin\lambda,\cos\lambda,0),\quad \hat n=\hat u\times\hat e=(-\sin\psi\cos\lambda,-\sin\psi\sin\lambda,\cos\psi).
$$
At a pole, the meridian is $\lambda=\operatorname{atan2}(y,x)$ of the point given, as `analysis`'s is.

**Proposition 3.2.** Outside its sources a potential field has $\nabla\cdot B=0$ and $\nabla\times B=0$.

*Proof.* $\nabla\times\nabla V=0$ for any $C^2$ potential, and each term $r^{-(n+1)}Y_n^m$ is harmonic, so $\nabla\cdot B=-\nabla^2V=0$ (Kellogg 1929, ch. VII). $\square$

The checks hold both by central differences. For a field that the dipole dominates, $\lvert\partial^3B\rvert\approx3\cdot4\cdot5\,\lvert B\rvert/r^3$, so a step $h=10^{-4}r$ gives a divergence below about $30(h/r)^2\lvert B\rvert/r=3\times10^{-7}\lvert B\rvert/r$:

- WMM2025 measures $2.3\times10^{-7}$ (geophysics, `check_wmm`);
- an expansion of degree 20 is held below $10^{-5}$ at $h=2\times10^{-4}a$ (`check_fields`).

On a sphere about the centre, the field is read row by row (`analysis` `field_grid`): the Legendre functions are computed once per latitude, in $O(LN^2+LMN)$ for $L$ rows and $M$ columns, against $O(LMN^2)$ point by point. The two agree to $10^{-13}$ (degree 20) and $7\times10^{-15}$ (T2015_449 at 30 km).

### 3.3 Uniform Fields

A uniform field has one value on its domain. Examples are the interplanetary field over the Moon, and the gravity of a small world.

## 4. Reductions and Composition

### 4.1 Reductions

**Definition 4.1 (reduction).** A *reduction* is a map $\rho\colon F\to F'$ of fibres at each point: $\lvert v\rvert$; the component $v\cdot e_i$ on the frame's axes; $v\cdot d$; the radial part $v\cdot\hat r$ about a centre $c$; the horizontal part $\lvert v-(v\cdot\hat r)\hat r\rvert$; the components on a local basis (the east, north and up of a sphere), or on given axes $Rv$.

**Proposition 4.2 (invariance).** Let a rigid motion $x'=Qx+b$ act on a vector field by $v'(x')=Qv(x)$. Then $\lvert v\rvert$ is invariant. The radial and horizontal parts about $c$ are invariant when $c$ is carried with the points ($c'=Qc+b$). A component on the frame's axes is not invariant.

*Proof.* $\lvert Qv\rvert=\lvert v\rvert$ for orthogonal $Q$, and $\hat r'=Q\hat r$, so that $v'\cdot\hat r'=v\cdot\hat r$. But $Qv\cdot e_i\ne v\cdot e_i$ in general. $\square$

Each reduction therefore declares whether it is invariant, and a drawing or a reading of a component must name its frame. This is formalism §3.6's remark on invariance, and the scene model's rule (the infrastructure's ADR 0053).

### 4.2 Sums, Pullbacks, Frames and Restrictions

The ways of making fields from fields are these:

- **Sum.** $a+b$ is defined on the intersection of the bases, and refuses where either part does, with that part's reason. Its parts must have one dimension. Its parts keep their provenance, and its own is the most derived of theirs, in the order physical, archive, model, assimilated, synthetic. That order runs from the record of a measurement to a value with no law behind it. It is an order of kinds of evidence, not of accuracy.
- **Pullback** along a map of bases $m\colon P\to Q$: $(m^*f)(t,p)=f(t,m(t,p))$, refused where $m$ is. The standard atmosphere, a field on geometric height, is carried into space by the height of each point above the geoid.
- **Change of frame** (passive, formalism §6.1). For the motion $x'=Rx+b$ from the field's frame to another, the same field is described by
$$
\varphi'(t,x')=\varphi(t,R^{\mathsf T}(x'-b)),\qquad v'(t,x')=R\,v(t,R^{\mathsf T}(x'-b)),
$$
and the derivatives turn as $\nabla'\varphi'=R\nabla\varphi$ and $J'=RJR^{\mathsf T}$ (chain rule, $\partial x/\partial x'=R^{\mathsf T}$).
- **Restriction** to a domain $D$: the field where $D$ covers, refused elsewhere with $D$'s reason.

## 5. Samplings and Sampled Fields

### 5.1 Samplings

A *sampling* is the restriction of formalism Definition 3.10 made finite: an ordered set of points of the base with a layout.

| Sampling | Points | Index |
|---|---|---|
| `SphereGrid` | cells of equal angle about a centre: row $i$ at latitude $90^\circ-(i+\tfrac12)\,180^\circ/n_\text{lat}$, column $j$ at longitude $-180^\circ+(j+\tfrac12)\,360^\circ/n_\text{lon}$ | $i\,n_\text{lon}+j$ |
| `PlaneGrid` | $o+\frac{i+1/2}{n_u}u+\frac{j+1/2}{n_v}v$ | $i+n_uj$ |
| `VolumeGrid` | $c+A^{\mathsf T}\big(((i+\tfrac12)/n_0-\tfrac12)e_0,\ldots\big)$, the rows of $A$ orthonormal axes | $i+n_0(j+n_1k)$ |
| `PointList` | given | given |

`sample(f, g, t)` reads $f$ at every point, by the representation's own method where it has one (a sphere by rows). It counts the answers by quality, and each reason or note once with the number of points it qualifies.

### 5.2 Sampled Fields

Samples are not a field until a reconstruction says what the field is between them. `Sampled` reconstructs a volume grid trilinearly between the cells' centres.

**Proposition 5.1.** Trilinear interpolation reproduces every affine field $f(x)=a+g\cdot x$ exactly.

*Proof.* The interpolant is $\sum_{\alpha\in\{0,1\}^3}w_\alpha f(x_\alpha)$ with $w_\alpha=\prod_k s_k^{\alpha_k}(1-s_k)^{1-\alpha_k}$. The weights sum to 1, and $\sum_\alpha w_\alpha x_\alpha=x$, since each is a tensor product of the one-dimensional identities $(1-s)x_0+sx_1=x$. Both are linear in $f$. $\square$

In the outer half of a boundary cell, between the last centres and the box's faces, the value is extended from the nearest centres, with a note. Beyond the box there is none. A point whose eight samples are not all valid has no value, with the reason of the first that is not.

## 6. Integral Curves

**Definition 6.1.** A *field line* of $v$ through a seed is the curve with $dx/ds=v/\lvert v\rvert$, parametrised by arc length; a *streamline* solves $dx/d\tau=v$.

Both are integrated by the classical fourth-order Runge–Kutta method (Hairer, Nørsett and Wanner 1993, §II.1), with a fixed step, or with a step that is a fraction of the distance from a centre, so that lines about a body take steps in proportion to their scale. A curve ends, in each direction, at one of four stops, which it reports:

- **the edge of the base.** A stage of the step has no value. The step is halved, up to `refine` times, so that the curve ends within $2^{-\text{refine}}$ of a step of the edge, and the end carries the edge's quality and reason;
- **a null of the field**, where $\lvert v\rvert\le$ zero and the direction is undefined;
- **the steps**, `max_steps` taken;
- **the length**, `max_length` reached.

**Proposition 6.2.** For the rigid rotation $v=\omega\times x$ ($\omega=\Omega\hat z$) and a step $h$ in $\tau$, one RK4 step multiplies the horizontal position, as a complex number, by
$$
R(i\theta)=1+i\theta-\frac{\theta^2}{2}-\frac{i\theta^3}{6}+\frac{\theta^4}{24},\qquad \theta=\Omega h,\qquad \lvert R(i\theta)\rvert^2=1-\frac{\theta^6}{72}+\frac{\theta^8}{576}.
$$
After $N$ steps the radius is $\lvert R\rvert^N$ times its first, and the angle $N\arg R=N(\theta-\theta^5/120+O(\theta^7))$.

*Proof.* On $\dot z=i\Omega z$, RK4 is multiplication by its stability function $R(w)=\sum_{k\le4}w^k/k!$ at $w=i\theta$ (Hairer and Wanner 1996, §IV.2). Expanding $\lvert R\rvert^2=(1-\theta^2/2+\theta^4/24)^2+(\theta-\theta^3/6)^2$ and $\arg R$ gives the rest. $\square$

The check holds a streamline's radius and angle to these values to $10^{-12}$ (`check_fields`). A field line of the same field (unit speed) turns by $h/r$ a step, and the proposition holds with $\theta=h/r$.

## 7. Versions in Time

A field that a law advances is published as a sequence of *versions*, each a field committed at the instant it describes and never changed afterwards. A query at $t$ reads the latest version committed at or before $t$; before the first, it is *missing*. A reader holds the version it read, a shared pointer to a constant field, so that a version committed while it reads does not change what it sees. The observers released at one instant then read one state (formalism Proposition 8.6).

## 8. Fields Behind One Type

`AnyField<V, P>` holds a field of value type $V$ on points $P$ whatever its representation. It keeps:

- its descriptor, with the capabilities read off its type;
- its batches;
- its exact derivative, when it has one;
- its own way of sampling a sphere.

The call through the erased type is made once for a batch or a sampling, not once a point, so that a registry of a world's fields (Cosmopoiesis's `fields.hpp`) costs no more than the fields. Numerical kernels keep the static types.

## 9. Geometric Complexes, and Cochains on Boxes

### 9.1 Topology and Geometry

The topology of a cellulation is morphe's CW complex. Its geometry is additional data, a position for each vertex in a frame. The topology fixes $d$, the geometry the Hodge stars, and the two are checked against each other.

**Definition 9.1.** The coboundary $d_k\colon C^k\to C^{k+1}$ is the transpose of the cellular boundary $\partial_{k+1}$, held as a sparse matrix with the integer incidences as values.

Then $d_{k+1}d_k=(\partial_{k+1}\partial_{k+2})^{\mathsf T}=0$, which morphe validates on every complex (Desbrun, Hirani, Leok and Marsden 2005, §3).

The measures of the cells are these:

- an edge's length;
- a 2-cell's vector area $\tfrac12\sum_i p_i\times p_{i+1}$ round its attaching word, whose magnitude is its area and whose direction is its unit normal (the word's orientation by the right-hand rule), with the centroid of its fan of triangles;
- a 3-cell's signed volume $\tfrac13\sum_f[c:f]\,(c_f\cdot A_f)$ by the divergence theorem applied to $x/3$, with its centroid from the pyramids on its faces. A pyramid's centroid lies a quarter of the way from its base to its apex.

These are exact for cells with planar faces, as every cell of a box is.

**Proposition 9.2 (orientation).** Let a complex of dimension $n$ be embedded in $\mathbb R^n$. Let each top cell $i$ be regular: a closed ball whose boundary sphere is the union of its facets, each entering $\partial i$ once, with incidence $\pm1$. Give each $(n-1)$-cell $f$ the conormal $\nu_f$ of its own orientation: a 2-cell's normal by the right-hand rule round its attaching word; an edge's direction turned a quarter in the plane; the line's direction at a vertex. Give each top cell the orientation $o_i=\operatorname{sign}\mathrm{vol}(i)$ of its signed measure. Then $d(i,f)=o_i\,\sigma(i,f)$, where $\sigma(i,f)=+1$ when $\nu_f$ points out of $i$ and $-1$ when it points in.

*Proof.* For $n=1$ directly: $\partial e=\text{head}-\text{tail}$, $o_e=\operatorname{sign}(x_\text{head}-x_\text{tail})$, and the conormal $+x$ points out of $e$ at its right end and into it at its left.

For $n\ge2$, $\partial_ni=\sum_f[i:f]\,f$ is an $(n-1)$-cycle of the cellulation of the boundary sphere $S^{n-1}$ by the facets. That cellulation has no $n$-cells, so its $(n-1)$-cycles are its homology $H_{n-1}(S^{n-1};\mathbb Z)\cong\mathbb Z$. This is generated by the fundamental class, in which each facet appears once, with coefficient $\sigma(i,f)$ relative to its own orientation. A cycle whose coefficients are $\pm1$ is therefore $\pm$ that class: $[i:f]=\varepsilon_i\,\sigma(i,f)$ for one sign $\varepsilon_i$.

Let $A_f$ be $\nu_f$ times the facet's measure. The divergence theorem applied to $x/n$ gives $\lvert\mathrm{vol}(i)\rvert=\tfrac1n\sum_f\sigma(i,f)\,(c_f\cdot A_f)$. The signed measure $\tfrac1n\sum_f[i:f]\,(c_f\cdot A_f)$ is then $\varepsilon_i\lvert\mathrm{vol}(i)\rvert$, so $\varepsilon_i=o_i$. Finally $d(i,f)=[i:f]$, since $d=\partial^{\mathsf T}$. $\square$

`measure()` lets a caller check this identity, and the checks do (`check_geometric`).

### 9.2 Boxes and Their Cochains

A box is the product of intervals cut at given coordinates (morphe `product`). Each cell knows the axes it spans and the indices of its least vertex, so that a point is located by a binary search along each axis. A box refined by $r$ cuts each interval into $r$ parts, and its cells are nested in the coarse box's.

A $k$-cochain on a box pairs with $k$-chains exactly. A value at a point needs a reconstruction, which a box gives in two degrees:

- for $k=0$, the multilinear interpolation of the cell's vertices: Whitney's 0-forms on a cube, the tensor products of the hat functions;
- for $k=n$, a density, the cell's value over its measure.

Other degrees are paired, not evaluated, and a query at a point refuses, saying why.

**Proposition 9.3 (transfer).** On a box refined from another, the de Rham map of the Whitney interpolant, $T=\mathcal R'\circ\mathcal W$, commutes with the coboundary: $Td=dT$ in every degree. For top cochains, $T$ keeps the total.

*Proof.* The box's cochain complex is the graded tensor product of its intervals' complexes, and $T$ is the tensor product of one-dimensional transfers:

- along an axis a cell spans, a fine edge takes the coarse edge's value times its share of the length, the integral of a constant density;
- along an axis it does not span, a fine vertex takes the linear interpolation of the coarse vertices about it.

On an interval, both commute with $d$: the sum of the shares of the fine edges in a coarse edge is the difference of the interpolated end values. A tensor product of chain maps is a chain map. The shares of a coarse top cell's fine cells sum to 1, which keeps the total. $\square$

This is formalism Proposition 6.2, which the checks hold in every degree (`check_geometric`).

## 10. A Balance Law: Heat Conduction

### 10.1 The Law and Its Discretisation

The energy density $u$, an $n$-form, balances its flux $J$, an $(n-1)$-form, and its sources $s$:
$$
\partial_tu+dJ=s,\qquad J=-k\star dT,\qquad u=\rho c\,T
$$
(Fourier 1822). On a cellulation the state is the energy $u_i$ of each top cell. The flux $J_f$ is the heat through each facet along its conormal. The closure is the two-point flux of discrete exterior calculus with a diagonal Hodge star: across an interior facet between the cell $a$ behind it and $b$ ahead,
$$
J_f=G_f(T_a-T_b),\qquad G_f=\frac{A_f}{d_a/k_a+d_b/k_b},
$$
with $A_f$ the facet's measure and $d_a,d_b$ the distances from the centroids to it along its conormal. These are the two half-cells' conductances in series. The cell's temperature is $T_i=u_i/C_i$, with heat capacity $C_i=\rho c\,V_i$ of its material.

A boundary facet is one of:

- insulated;
- held at a temperature, with $G=k_aA_f/d_a$;
- given a flux density;
- exchanging with an outside temperature through a coefficient $h$, with $G=A_f/(1/h+d_a/k_a)$.

**Proposition 10.1 (exact affine states).** On a box whose cells each hold one material, and whose interfaces between materials are planes of cuts across the axis of conduction, a steady temperature that is affine in each material along that axis, with its flux continuous across the interfaces, is reproduced exactly.

*Proof.* On a box the line between two centroids is along the facet's conormal. Within a material the two-point flux of an affine $T$ is therefore exactly $-k\,(\partial T/\partial n)\,A_f$. Across an interface the series conductance is the exact conductance of the two half-layers. So the discrete fluxes are the continuous ones, which balance in each cell. $\square$

The checks hold a linear profile and a two-material profile to $10^{-13}$ (`check_heat`).

The semi-discrete system is $C\dot T=-LT+b$, with $L=d\,G\,d^{\mathsf T}$ symmetric positive semidefinite. It is definite once any facet is held or exchanging. $b$ is the forcing of the boundary and the sources. Duplicate contributions to an entry of $L$ are summed at assembly, so that each row's diagonal is its own.

### 10.2 The θ-Step and Conservation

A θ-step of calx `integrators` (calx decision 0011) gives the new temperatures $T^*$:
$$
(C+\theta hL)\,T^*=(C-(1-\theta)hL)\,T+hb.
$$
The energy is then updated by the fluxes of the weighted temperatures $T_\theta=\theta T^*+(1-\theta)T$:
$$
u_i'=u_i-h\sum_f\sigma(i,f)\,J_f(T_\theta)+hs_i,\qquad T'=u'/C.
$$

**Proposition 10.2 (conservation whatever the solver).** In floating point, the total energy changes over a step by $h$ times the boundary's net inflow and the sources, up to the rounding of the sums. Between two materials, what leaves one enters the other. Both hold whatever the tolerance of the linear solve for $T^*$.

*Proof.* Each interior facet's $J_f$ enters the update of the cell behind it with $\sigma=+1$ and of the cell ahead with $\sigma=-1$. The same number is used both times, so in $\sum_iu_i'$ every interior term cancels exactly, whatever $T_\theta$ is. Only the boundary terms and the sources remain. The same cancellation, restricted to the cells of one material, leaves its boundary terms and its interface terms, which the other material has with opposite signs. $\square$

In exact arithmetic $T'=T^*$, which is the θ-method itself. The solve's tolerance bounds only how far $T'$ is from it. The check makes the solve loose ($10^{-6}$) and finds the budget closed to $10^{-17}$ of the stored energy (`check_heat`).

### 10.3 Stability and Dissipation

**Proposition 10.3.** For an insulated body without sources, $E=\tfrac12\sum_iC_iT_i^2$ does not increase over a step when $\theta\ge\tfrac12$, for every $h$. When $\theta=0$, it does not increase if $h\lambda_\text{max}\le2$, where $\lambda_\text{max}$ is the largest eigenvalue of $C^{-1}L$.

*Proof.* On an eigenvector of $C^{-1}L$, with eigenvalue $\lambda\ge0$ and $z=h\lambda$, the step multiplies by $R(z)=(1-(1-\theta)z)/(1+\theta z)$. The eigenvectors are $C$-orthogonal, so $E$ is a sum of their squared amplitudes, and $\lvert R(z)\rvert\le1$ for all $z\ge0$ exactly when $\theta\ge\tfrac12$, or when $\theta<\tfrac12$ and $z\le2/(1-2\theta)$ (Hairer and Wanner 1996, §IV.3). $\square$

**Proposition 10.4 (an explicit step that is always stable).** $h\le\min_iC_i/L_{ii}$ gives $h\lambda_\text{max}\le2$.

*Proof.* By Gershgorin's theorem every eigenvalue of $C^{-1}L$ lies in a disc about $L_{ii}/C_i$ of radius $\sum_{j\ne i}\lvert L_{ij}\rvert/C_i\le L_{ii}/C_i$ (diagonal dominance of a Laplacian with nonnegative boundary terms). So $\lambda_\text{max}\le\max_i2L_{ii}/C_i$. $\square$

The checks hold the following (`check_heat`):

- the dissipation, at each step;
- the bound;
- the explicit scheme's instability just above $2/\lambda_\text{max}$;
- convergence to the decaying sine mode of the heat equation: second order in space on uniform boxes (the two-point flux's supraconvergence) and in time at $\theta=\tfrac12$, first order in time at $\theta=1$.

## 11. References

- Desbrun, M., Hirani, A. N., Leok, M. and Marsden, J. E. (2005), Discrete exterior calculus, arXiv:math/0508341.
- Fourier, J. (1822), *Théorie analytique de la chaleur*, Firmin Didot, Paris.
- Griewank, A. and Walther, A. (2008), *Evaluating Derivatives*, 2nd ed., SIAM.
- Hairer, E., Nørsett, S. P. and Wanner, G. (1993), *Solving Ordinary Differential Equations I*, 2nd ed., Springer.
- Hairer, E. and Wanner, G. (1996), *Solving Ordinary Differential Equations II*, 2nd ed., Springer.
- Kellogg, O. D. (1929), *Foundations of Potential Theory*, Springer.
- Whitney, H. (1957), *Geometric Integration Theory*, Princeton University Press.
- Arnold, D. N., Falk, R. S. and Winther, R. (2006), Finite element exterior calculus, homological techniques, and applications, *Acta Numerica* 15, 1–155.
