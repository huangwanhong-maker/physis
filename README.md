# physis

*Geometry and theoretical physics.* The name comes from the Greek *φύσις*, "nature".

## 1. Purpose

physis implements general physical laws, and the geometry they are written in. It does not depend on any particular body.

It is where the project's two kinds of mathematics meet physical law:

- the structures of [morphe](https://github.com/huangwanhong-maker/morphe), such as groups and cell complexes;
- the numerical methods of [calx](https://github.com/huangwanhong-maker/calx).

physis is one of the foundation libraries of [The Science Infrastructures](https://github.com/huangwanhong-maker/The-Science-Infrastructures). It stands on its own.

## 2. Planned Modules

| Module | Content | Status |
|---|---|---|
| [`units`](units/) | Physical dimensions as a free abelian group over the SI base units; quantities and conversions; CODATA constants, versioned | draft: dimensions and quantities |
| [`symmetry`](symmetry/) | Lie groups and algebras as physics uses them: U(1), SU(N), SO(N), the Lorentz group and SL(2, ℂ); Haar measures, characters, the Weyl integration formula; representations and the recoupling (3nj) symbols, classical and quantum | draft: U(1), SU(2), SU(N), Lorentz, exact 3j/6j/9j, quantum 6j; point groups of labelled points (decision 0011) |
| [`geometry`](geometry/) | Euclidean, affine and projective geometry; spherical and ellipsoidal geometry; rotations (SO(3), SE(3), quaternions); manifolds, metrics, connections, geodesics, curvature, differential forms; information geometry; map projections | draft: metrics, curvature, geodesics; rotations, moving frames and the frame graph with its holonomy (decision 0009) |
| [`spacetime`](spacetime/) | Minkowski space and causal structure; exact solutions of general relativity in charts; discrete spacetimes: Regge geometry, causal sets, causal dynamical triangulations | draft: Minkowski, exact solutions, Regge geometry, causal sets, emergent dimensions |
| [`computational_geometry`](computational_geometry/) | Convex hulls, triangulations, Voronoi diagrams, meshes; spatial indices; robust predicates | draft: points on the sphere, their spherical Delaunay triangulation, angles and areas; the faces of a plane graph from its rotation system, polylines (decision 0018) |
| [`analysis`](analysis/) | Real and complex analysis as computation needs it; function spaces, operators, spectral theory; Fourier and harmonic analysis on groups and spheres; integral transforms and equations; calculus of variations | draft: potential fields outside a sphere (internal spherical-harmonic expansions at points and on grids, to high degree at every latitude), their spectra, random fields with a given spectrum and the part an estimate leaves unresolved (decision 0010); the covariance of such a field between any two points, of its potential and of its field, exact or tabulated at one radius (decision 0013); stationary random fields of a squared-exponential covariance read at any point by the spectral representation (decision 0015) |
| [`dynamics`](dynamics/) | State spaces, vector fields and flows; linear systems and their exact propagators; hybrid systems (flow, guard, jump); discrete-time, delay and stochastic dynamics; stability, bifurcation, Lyapunov exponents, invariant measures | draft: linear propagators, hybrid systems; piecewise-deterministic processes with forced jumps, exact linear elimination and its memory kernel (decision 0018) |
| `control` | Systems with inputs and outputs: controllability, observability, feedback, optimal and model-predictive control, system identification | planned |
| [`quantum`](quantum/) | Hilbert spaces, states and operators; spins and qubits; unitary and open (GKSL) dynamics, as Liouvillians or generalised Bloch vectors; quantum trajectories; measurement; mean-field closures | draft: the open two-level system and its trajectories |
| [`statistical_mechanics`](statistical_mechanics/) | Ensembles and partition functions; lattice and spin models (Ising, Potts, Heisenberg, transverse-field); Glauber and Metropolis dynamics; mean-field theory and phase transitions; Markov jump processes and their irreversibility | draft: master equations (stationary laws by GTH, relaxation, entropy production and cycle affinities), first passage and phase-type laws, the entropy of mixing (decision 0016); ensembles and spin models planned; open systems composed at relations, lumpability with its discrepancy and coupling bound, anchored jump processes (decision 0018) |
| [`fields`](fields/) | What a field is, independently of its representation: a section over a base, with its fibre, representation and provenance; evaluations that refuse with a reason; capabilities as concepts; reductions, sums, pullbacks, changes of frame and restrictions; samplings; integral curves; versions in time | draft: formulas with exact derivatives, potential fields, uniform and sampled fields, curves, versions, type erasure for registries (decision 0012) |
| [`discretisation`](discretisation/) | Finite-difference, finite-volume, finite-element and spectral methods. Structure-preserving discretisation: discrete exterior calculus on morphe's complexes. Geometric and Lie-group integrators. | draft: cochains, potentials and the Hodge decomposition (dense and by conjugate gradients); geometric complexes (coboundaries in every degree, measures, facets and their orientation), boxes, and cochains on them as fields with their transfer to a refinement (decision 0012); a JavaScript stand-in for the spectral transform on the sphere and longitude-latitude grids is in [`js/`](js/); Hodge Laplacians of every degree with diagonal stars (decision 0018) |
| [`field_theory`](field_theory/) | Classical fields and Noether's theorem; lattice field theory on complexes; Monte Carlo and its statistics; Feynman graphs and perturbation theory; free fields; renormalisation; conformal and topological field theories | draft: scalar fields on complexes, statistics, Feynman graphs, invariant theories of every rank, large-N solutions |
| [`gauge`](gauge/) | Connections and curvature; lattice gauge theory on any 2-complex; Monte Carlo of gauge fields; Wilson and Polyakov loops; character expansions; exact two-dimensional Yang–Mills; BF theory; Chern–Simons | draft: gauge fields on CW complexes with couplings by 2-cell and frozen links, updates (the heat bath of finite groups among them), observables, Migdal–Witten, flat connections, the exact contraction of ℤ₂ theories with derivatives |
| [`quantum_gravity`](quantum_gravity/) | Spin networks; spin foams (Ponzano–Regge, Turaev–Viro, Barrett–Crane, EPRL); group field theory and tensor models; matrix models; causal dynamical triangulations; causal-set dynamics | draft: Turaev–Viro, BF amplitudes, group field theories, spin foams, random geometries, condensates, tensor and matrix models; stranded theories on real fields over finite groups (decision 0011) |
| `electromagnetism` | Maxwell's equations, potentials and gauges, radiation, fields in media; on complexes by discrete exterior calculus | planned |
| [`mechanics`](mechanics/) | Lagrangian and Hamiltonian mechanics; rigid-body rotation; two-body and n-body problems; perturbation theory | draft: the two-body problem: elements, Kepler's equation, universal propagation, the frames of an orbit (decision 0009); motion from a discrete action (variational integrators, discrete Noether, forces, the Euclidean action) and rigid bodies (mass properties, the Lie group variational integrator, wrenches from a potential) (decision 0013); planar trees of rigid links on a fixed or floating base, by their Lagrangian and by Newton–Euler, with explicit integrators, and compliant contact with the ground (decision 0014); the Euler–Lagrange equations of any Lagrangian solved for the acceleration by jets, and a tree in a fluid by strip theory, its added mass in the Lagrangian and its drag a force (decision 0017) |
| `gravitation` | Potential theory; spherical-harmonic gravity fields; tidal potentials. General relativity: metrics, post-Newtonian dynamics, light propagation. | planned |
| `thermodynamics` | Equilibrium thermodynamics, equations of state, phase changes. The common base for seawater and moist air. | planned |
| [`continuum`](continuum/) | Balance laws of continuous media in conservation form. Navier–Stokes equations. Geophysical fluid dynamics: rotating and stratified flow, shallow-water and primitive equations, Boussinesq and anelastic approximations. Waves. Transport. | draft: heat conduction as a balance law, conservative whatever the solver (decision 0012); kinematic waves on networks of links by Godunov's scheme, and advection and diffusion on graphs (decision 0016) |
| [`chemistry`](chemistry/) | Reaction kinetics; chemical thermodynamics; coupling of reactions to transport | draft: reaction networks, their conservation laws and balance, mass action, and their Markov chains solved exactly (decision 0011); reaction–diffusion on graphs, Turing's instability and the pinning of fronts on a chain (decision 0016); Doi–Peliti operators equal to the master equation (decision 0018) |
| `radiation` | Radiative transfer: absorption, emission and scattering | planned |

## 3. Principles

The principles are set out in full in [docs/architecture.md](docs/architecture.md).

- **General laws, not particular bodies.** physis knows potentials, fluids and radiation. It does not know about the Earth, the Moon or a particular ocean.
- **Explicit units.** Units are explicit at every interface.
- **Intrinsic geometry where possible.** Coordinates are explicit charts on top of it.
- **Cited derivations.** Every implemented law cites its derivation.

## 4. Relation to Sibling Libraries

physis depends on calx and on morphe.

- **Inside the infrastructure**, it uses the copies of calx and morphe already present under `foundations/`.
- **Built on its own**, it fetches them through the build system.

No library of the project depends on the infrastructure's domains or applications.

## 5. Status

**Draft (2026-10-24).** The first C++ modules, `units`, `dynamics` and `quantum` (the open two-level system), are written and checked, as the neuroscience domain's spiking networks need them ([ADR 0004](development/decisions/0004-the-continuum-units-and-dynamics.md)); with them, the theoretical physics of [decision 0005](development/decisions/0005-symmetry-spacetime-fields-and-gauge.md): `symmetry`, `geometry` (metrics, curvature and geodesics), `spacetime`, `field_theory`, `gauge` and `quantum_gravity`, on morphe's complexes, with their conventions in the [glossary](docs/glossary.md); and the generation of spacetime by field theories of [decision 0006](development/decisions/0006-the-generation-of-spacetime-by-field-theories.md), from an action through its Feynman graphs, spin foams and condensates to the continuum; and the exact contraction of ℤ₂ gauge theories with a coupling for each 2-cell, with the heat bath of finite groups ([decision 0007](development/decisions/0007-exact-z2-contraction-and-finite-heat-bath.md)); discrete exterior calculus and the Hodge decomposition in `discretisation`, and points on the sphere in `computational_geometry` ([decision 0008](development/decisions/0008-cochains-hodge-and-the-sphere.md)); rotations, moving frames and the frame graph with its holonomy in `geometry`, and the two-body problem in `mechanics` ([decision 0009](development/decisions/0009-rotations-frames-and-the-two-body-problem.md)); fields as sections over their bases, geometric complexes and cochains on boxes, and heat conduction as the first balance law in `continuum` ([decision 0012](development/decisions/0012-fields-as-sections.md), with its [theory](docs/theory/fields.md)); motion from a discrete action and rigid bodies in `mechanics`, and the covariance of a random potential field in `analysis` ([decision 0013](development/decisions/0013-discrete-mechanics-and-the-law-of-a-field-known-in-part.md), with its [theory](docs/theory/discrete-mechanics.md)); articulated bodies in the plane and compliant contact in `mechanics` ([decision 0014](development/decisions/0014-articulated-bodies-in-the-plane-and-compliant-contact.md)); stationary random fields read at any point in `analysis` ([decision 0015](development/decisions/0015-stationary-random-fields.md)); bodies in a fluid and the Euler–Lagrange equations solved for the acceleration in `mechanics` ([decision 0017](development/decisions/0017-bodies-in-a-fluid-and-lagrangian-dynamics.md)); kinematic waves on networks and transport on graphs in `continuum`, master equations, first passage and the entropy of mixing in `statistical_mechanics`, and reaction–diffusion on graphs in `chemistry` ([decision 0016](development/decisions/0016-kinematic-waves-master-equations-and-reaction-diffusion.md)); a provisional JavaScript stand-in for `discretisation` is in [`js/`](js/) ([ADR 0003](development/decisions/0003-provisional-javascript-stand-in.md)). The rest is designed and planned; the place of each module in the whole of the project's mathematics is the infrastructure's [map of the mathematics](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/docs/mathematics.md).

Inside The Science Infrastructures, build and check with CMake 3.20 or later (calx and morphe are taken from `../calx` and `../morphe`):

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## 6. License

physis is released under the MIT License. See [`LICENSE`](LICENSE).
