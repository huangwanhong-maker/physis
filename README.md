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
| `units` | Physical dimensions as a free abelian group over the SI base units; quantities and conversions; CODATA constants, versioned | planned |
| `geometry` | Euclidean, spherical and ellipsoidal geometry. Rotations: SO(3), SE(3), quaternions. Differential and Riemannian geometry: manifolds, metrics, connections, geodesics, curvature, differential forms. Map projections. | planned |
| `discretisation` | Finite-difference, finite-volume, finite-element and spectral methods. Structure-preserving discretisation: discrete exterior calculus on morphe's complexes. Geometric and Lie-group integrators. | planned; a JavaScript stand-in for the spectral transform on the sphere and longitude-latitude grids is in [`js/`](js/) |
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

**Planned (2026-10-21), with a first stand-in (2026-10-22).** Design documents exist, and a provisional JavaScript stand-in for `discretisation`, tested, in [`js/`](js/) ([ADR 0003](development/decisions/0003-provisional-javascript-stand-in.md)). There is no C++ yet.

## 6. License

physis is released under the MIT License. See [`LICENSE`](LICENSE).
