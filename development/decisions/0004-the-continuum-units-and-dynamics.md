# 0004. The Continuum: New Modules, and the First C++ in Units and Dynamics

- **Status:** Accepted, provisionally (the representation of units in §2.3 is a first answer to open decision 1)
- **Date:** 2026-10-23
- **Decided by:** the user asked that mathematics common to the sciences be designed as a whole in the foundations and written as the neuroscience domain's spiking networks need it, saying that they would themselves make the dynamics abstract (states abstract, dynamics such as PDEs implemented by various models). Claude designed physis's part and wrote this record. The whole map is The Science Infrastructures' [ADR 0035](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0035-the-mathematical-foundations.md) and [map of the mathematics](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/docs/mathematics.md).

## 1. Context

physis was planned as geometry and the laws of physics. The map of the mathematics gives it the continuum where structures meet numbers: geometry, analysis, dynamics, control and discretisation, beside the physics written in them. A neuron, a box model of the climate and a discretised ocean are all dynamical systems; their flows, events and invariants deserve one home, which calx (knowing only arrays and functions) is not, and which no domain should own.

## 2. Decision

1. **Six new modules in the plan:** `computational_geometry` (hulls, triangulations, meshes, spatial indices, robust predicates), `analysis` (function spaces and operators, harmonic analysis on groups and spheres, integral transforms and equations, calculus of variations), `dynamics` (state spaces, flows, linear propagators, hybrid, discrete-time, delay and stochastic systems, stability and bifurcation), `control` (systems with inputs and outputs), `quantum` (states and operators, unitary and open GKSL dynamics, quantum trajectories, measurement), which the user asked for so that quantum neurons, Ising-like and open-system, can be simulated, and `statistical_mechanics` (ensembles, spin models, Glauber dynamics, mean field). The module order: `computational_geometry` on `geometry`; `analysis` and `dynamics` on `units`; `control` on `dynamics`; `quantum` on `geometry` and `dynamics`; `statistical_mechanics` on `thermodynamics`, `dynamics` and `quantum`; `discretisation` on `geometry`, `analysis` and `dynamics`; `mechanics` on `dynamics` too.
2. **Written now, header-only, status `draft`:** `units` (dimensions over ℤ⁷, quantities read with their units) and `dynamics` (exact propagators of linear systems by Van Loan's augmented exponential, scaled by the components' physical sizes; the concepts of linear, nonlinear, gated and hybrid systems; one step by calx's integrators or by Rush and Larsen's method; one trajectory with its events).
3. **Units at run time:** a dimension is a vector of seven integer exponents, checked where quantities are read (from a specification, a file, a user); inside, values are SI doubles. Compile-time types stay open.
4. **Layout, names and build** as calx's decision 0004: `<module>/include/physis/<module>/`, the namespace `physis::<module>`, the target `physis::<module>`, CMake 3.20 and C++20, checks on calx's `verification`; inside the infrastructure, physis adds the sibling `../calx` when no larger project has.

## 3. Alternatives Considered

- **Dynamics in calx:** calx's integrators know a function from vectors to vectors; a system with a guard, a jump and a dwell, and a propagator built from a model's matrices, are structure that calx does not hold.
- **Dynamics in the neuroscience domain:** the abstraction serves every domain that simulates, and the user asked for it to be separate.
- **Units as compile-time types now:** they would type every quantity of every model; run-time dimensions at the boundaries give the safety that matters (a conductance never read as a current) at no cost inside the loops, and leave the choice open.

## 4. Consequences

- physis has its first C++ beside its JavaScript stand-in: 2 check programs, 14 checks.
- The neuroscience domain's neuron models are physis hybrid systems; other domains state their systems the same way.
- physis's README, architecture and CLAUDE.md list the new modules and the build.
