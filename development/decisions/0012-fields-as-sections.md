# 0012. Fields as Sections, Geometric Complexes, and the First Balance Law

- **Status:** Accepted, provisionally (boxes as the first geometric complexes; heat conduction as the first balance law; capabilities as concepts, not yet a solver interface)
- **Date:** 2026-10-29
- **Decided by:** the user asked for "a unified way, model various fields (tensor, scalar, even manifold or complexes), and their formalism. Then the visualizations for them", with the evolution of fields "in structured way, like PDE … but separated"; a reviewer's design was read critically, and its separation of base from representation, its capabilities and its conserved state were kept, its general solver and its virtual interfaces not. Claude wrote this module and this record. Their use is Cosmopoiesis's [decision 0007](https://github.com/huangwanhong-maker/Cosmopoiesis-Lab/blob/main/development/decisions/0007-the-worlds-fields-and-their-drawing.md) and the infrastructure's [ADR 0053](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0053-fields-and-their-visual-expression.md).

## 1. Context

The infrastructure's domains each had their models as functions of their own: WMM2025 at a geodetic place, the standard atmosphere at a height, the Moon's crust at a radius and an angle, normal gravity at a place. Each application that drew or sampled one wrote its own loop, grid, refusal and colouring, with its own conventions. There was no common answer to three questions:

- what a field is, whatever the model;
- where it has a value and why not;
- how it is held for computation.

The cochains of decision 0008 had a topology and no geometry, and there was no law that advanced a field in time.

## 2. Decision

### 2.1 The Module `fields`

A new module, `fields` ([fields/README.md](../../fields/README.md); [docs/theory/fields.md](../../docs/theory/fields.md)), holds what a field is independently of how it is stored, sampled, solved or drawn:

- **A field is a section of a bundle over a base.** Its base, fibre, representation and provenance are kept apart in its descriptor (theory §1).
- **An evaluation is a value or one of three refusals.** The refusals are outside, missing and invalid, each with its reason, and never a zero. A valid value may carry a note (theory §1.2).
- **Capabilities are C++20 concepts on the representation**, not virtual methods every field must have. They are point evaluation, batches, an exact derivative, pairing with chains, a coboundary, a Hodge star, fast sampling of a sphere, transfer and versions. They are read off the type and reported in the descriptor (theory §1.3).
- **Representations.** They are formulas, with exact derivatives by calx's dual numbers; potential fields over `analysis`'s expansions, read by rows on spheres; uniform fields; samples on a volume grid with trilinear reconstruction; and versions committed in time (theory §3, §5, §7).
- **Reductions and compositions.** Reductions declare their invariance. Sums give a composite of the most derived provenance of its parts. There are pullbacks along a map of bases, changes of frame that turn vectors and Jacobians, and restrictions (theory §4).
- **Samplings.** These are spheres, planes, volumes and point lists, each with its order of points, and sampled answers count their reasons and notes (theory §5).
- **Integral curves.** Field lines and streamlines are integrated by RK4, each end saying why it ended. The step is halved at the edge of the base (theory §6).
- **Type erasure for registries.** `AnyField` makes one call per batch (theory §8).

### 2.2 Geometric Complexes and Cochain Fields in `discretisation`

`discretisation` gains `geometric.hpp` (theory §9):

- coboundaries in every degree from morphe's chain complex;
- the measures of cells: lengths, vector areas, signed volumes and centroids;
- facets with their conormals, and the orientation identity $d(i,f)=o_i\sigma(i,f)$;
- boxes as morphe products, with location of points and refinement.

It also gains `cochain_field.hpp`. A cochain on a box is a field: it pairs with chains, and is evaluated at points in degrees 0 and $n$ by Whitney's reconstruction. Its transfer to a refinement commutes with $d$.

### 2.3 The First Balance Law in `continuum`

`continuum` is written, its first content heat conduction as a balance law in conservation form ([continuum/README.md](../../continuum/README.md); theory §10):

- the two-point flux with the series conductance of half-cells;
- four kinds of boundary;
- the θ-step of calx (calx decision 0011);
- the energy then updated by the fluxes, so that it is conserved to rounding whatever the linear solver's tolerance.

The law is separate from the field and from the solver. The state is the conserved quantity on the cells. The law gives the flux and the rate. The time integrator is calx's. The result is published as committed versions of a cochain field.

### 2.4 Dependencies

`fields` depends on `units`, `geometry` and `analysis`, and on calx `differentiation`. `discretisation` now depends on `fields`, and `continuum` on `discretisation`, `fields` and calx `integrators` and `linalg` (architecture §2).

## 3. Alternatives Considered

- **A class for each kind of field** (scalar, vector, tensor, form, complex). Rejected: the kind is two axes, base and fibre, and the class hierarchy would multiply by the representations.
- **A vector-calculus language with one general `solve()`.** Rejected: a solver depends on the law's structure (conservation form, symmetry, stiffness), which a general interface would hide. Each law states its discretisation and uses calx's integrators for its time.
- **Virtual interfaces for capabilities.** Rejected: they would put a virtual call in every evaluation of a kernel. Concepts keep the static types for the kernels, and `AnyField` erases them once per batch where a registry needs one type.
- **Regge edge lengths as a 1-cochain.** Rejected: a metric is not a cochain. It does not pair with chains linearly or transform by $d$. It is a field of another fibre on the same base.
- **Numerics in physis beside the laws.** Rejected: the θ-method, conjugate gradients and the sparse rows are calx's, with no physics in them.

## 4. Consequences

- physis has 31 check programs. The three new ones are:
  - `check_physis_fields` (49 checks): a dipole by dual numbers; central differences against the derived limit; potential fields against their degree-1 closed form, and on spheres by rows; reductions; compositions; samplings; curves against RK4's stability function; versions; samples with their reconstruction; type erasure;
  - `check_physis_geometric` (18 checks): coboundaries, measures, facets and the orientation identity on boxes of one, two and three dimensions; the transfer commuting with $d$;
  - `check_physis_heat` (17 checks): affine and two-material states exact; the budget closed to $10^{-17}$ of the stored energy with a loose solve; per-material budgets; dissipation; the Gershgorin step; instability above it; convergence of second order in space and in time.
- The domains present their models as fields: geodesy's normal gravity, geophysics's WMM2025 (on its base, and continued outward), meteorology's standard atmosphere on its line of heights, and planetary science's lunar crust.
- A bug was found and fixed while the heat law was written: duplicate entries of the conductance Laplacian were stored apart, so that the diagonal read for the explicit bound was a part of it.
- A sum's capabilities are the intersection of its parts' fast samplings, so that a composite of two spectral fields still reads a sphere by rows.
