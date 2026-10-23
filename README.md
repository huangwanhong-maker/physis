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
| [`symmetry`](symmetry/) | Lie groups and algebras as physics uses them: U(1), SU(N), SO(N), the Lorentz group and SL(2, ℂ); Haar measures, characters, the Weyl integration formula; representations and the recoupling (3nj) symbols, classical and quantum | draft: U(1), SU(2), SU(N), Lorentz, exact 3j/6j/9j, quantum 6j |
| [`geometry`](geometry/) | Euclidean, affine and projective geometry; spherical and ellipsoidal geometry; rotations (SO(3), SE(3), quaternions); manifolds, metrics, connections, geodesics, curvature, differential forms; information geometry; map projections | draft: metrics, curvature, geodesics |
| [`spacetime`](spacetime/) | Minkowski space and causal structure; exact solutions of general relativity in charts; discrete spacetimes: Regge geometry, causal sets, causal dynamical triangulations | draft: Minkowski, exact solutions, Regge geometry, causal sets |
| `computational_geometry` | Convex hulls, triangulations, Voronoi diagrams, meshes; spatial indices; robust predicates | planned |
| `analysis` | Real and complex analysis as computation needs it; function spaces, operators, spectral theory; Fourier and harmonic analysis on groups and spheres; integral transforms and equations; calculus of variations | planned |
| [`dynamics`](dynamics/) | State spaces, vector fields and flows; linear systems and their exact propagators; hybrid systems (flow, guard, jump); discrete-time, delay and stochastic dynamics; stability, bifurcation, Lyapunov exponents, invariant measures | draft: linear propagators, hybrid systems |
| `control` | Systems with inputs and outputs: controllability, observability, feedback, optimal and model-predictive control, system identification | planned |
| [`quantum`](quantum/) | Hilbert spaces, states and operators; spins and qubits; unitary and open (GKSL) dynamics, as Liouvillians or generalised Bloch vectors; quantum trajectories; measurement; mean-field closures | draft: the open two-level system and its trajectories |
| `statistical_mechanics` | Ensembles and partition functions; lattice and spin models (Ising, Potts, Heisenberg, transverse-field); Glauber and Metropolis dynamics; mean-field theory and phase transitions | planned |
| `discretisation` | Finite-difference, finite-volume, finite-element and spectral methods. Structure-preserving discretisation: discrete exterior calculus on morphe's complexes. Geometric and Lie-group integrators. | planned; a JavaScript stand-in for the spectral transform on the sphere and longitude-latitude grids is in [`js/`](js/) |
| [`field_theory`](field_theory/) | Classical fields and Noether's theorem; lattice field theory on complexes; Monte Carlo and its statistics; Feynman graphs and perturbation theory; free fields; renormalisation; conformal and topological field theories | draft: scalar fields on complexes, Markov-chain statistics, Feynman graphs, zero-dimensional theory |
| [`gauge`](gauge/) | Connections and curvature; lattice gauge theory on any 2-complex; Monte Carlo of gauge fields; Wilson and Polyakov loops; character expansions; exact two-dimensional Yang–Mills; BF theory; Chern–Simons | draft: gauge fields on CW complexes, updates, observables, Migdal–Witten, flat connections |
| [`quantum_gravity`](quantum_gravity/) | Spin networks; spin foams (Ponzano–Regge, Turaev–Viro, Barrett–Crane, EPRL); group field theory and tensor models; matrix models; causal dynamical triangulations; causal-set dynamics | draft: Turaev–Viro, BF amplitudes, coloured tensor models, matrix-model gluings |
| `electromagnetism` | Maxwell's equations, potentials and gauges, radiation, fields in media; on complexes by discrete exterior calculus | planned |
| `mechanics` | Lagrangian and Hamiltonian mechanics; rigid-body rotation; two-body and n-body problems; perturbation theory | planned |
| `gravitation` | Potential theory; spherical-harmonic gravity fields; tidal potentials. General relativity: metrics, post-Newtonian dynamics, light propagation. | planned |
| `thermodynamics` | Equilibrium thermodynamics, equations of state, phase changes. The common base for seawater and moist air. | planned |
| `continuum` | Navier–Stokes equations. Geophysical fluid dynamics: rotating and stratified flow, shallow-water and primitive equations, Boussinesq and anelastic approximations. Waves. Transport. | planned |
| `chemistry` | Reaction kinetics; chemical thermodynamics; coupling of reactions to transport | planned |
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

**Draft (2026-10-23).** The first C++ modules, `units`, `dynamics` and `quantum` (the open two-level system), are written and checked, as the neuroscience domain's spiking networks need them ([ADR 0004](development/decisions/0004-the-continuum-units-and-dynamics.md)); with them, the theoretical physics of [decision 0005](development/decisions/0005-symmetry-spacetime-fields-and-gauge.md): `symmetry`, `geometry` (metrics, curvature and geodesics), `spacetime`, `field_theory`, `gauge` and `quantum_gravity`, on morphe's complexes, with their conventions in the [glossary](docs/glossary.md); a provisional JavaScript stand-in for `discretisation` is in [`js/`](js/) ([ADR 0003](development/decisions/0003-provisional-javascript-stand-in.md)). The rest is designed and planned; the place of each module in the whole of the project's mathematics is the infrastructure's [map of the mathematics](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/docs/mathematics.md).

Inside The Science Infrastructures, build and check with CMake 3.20 or later (calx and morphe are taken from `../calx` and `../morphe`):

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## 6. License

physis is released under the MIT License. See [`LICENSE`](LICENSE).
