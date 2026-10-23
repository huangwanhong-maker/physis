# 0005. Symmetry, Spacetime, Fields, Gauge and Quantum Gravity

- **Status:** Accepted, provisionally (the conventions in §2.4 answer open decision 2 in part; §2.5 answers open decision 4)
- **Date:** 2026-10-23
- **Decided by:** the user asked for theoretical physics to be written in parallel with morphe's algebra and topology, its structure planned and only what is needed now implemented: spacetime models, the infrastructure of field theory (quantum field theory and group field theory especially), and tools for research on symmetry and gauge, on the 2-complexes and CW complexes of morphe. Claude designed physis's part and wrote this record; morphe's part is its decision 0004, calx's its decision 0005, and the composition The Science Infrastructures' [ADR 0039](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0039-complexes-symmetry-and-fields.md).

## 1. Context

physis planned geometry and the classical laws, and since decision 0004 the continuum and quantum theory. It had no home for the theories the user now studies: symmetry as Lie groups and their representations; spacetime as an object (continuous or discrete); fields on spacetimes and on complexes, quantised by path integrals; gauge fields; and the approaches to quantum gravity in which these meet, spin foams and group field theory. The infrastructure's open decision 10 left to physis which modules of theoretical physics it adds.

## 2. Decision

### 2.1 Six New Modules

| Module | Content (plan) | Written now (`draft`) |
|---|---|---|
| `symmetry` | Lie groups and algebras as physics uses them: U(1), SU(N), SO(N), Sp(N), the Lorentz and Poincaré groups, SL(2, ℂ); exponential maps, Haar measures, characters and the Weyl integration formula; representations, Wigner D-functions, Clebsch–Gordan series and the recoupling (3nj) symbols, and their quantum deformations at roots of unity; symmetry breaking and Noether currents | U(1); SU(2) as unit quaternions, with its algebra, exponential, characters and Haar measure; SU(N) as matrices, with reunitarisation and Haar sampling; the Lorentz group and the covering SL(2, ℂ) → SO⁺(1, 3); the 3j, 6j and 9j symbols, exact as signed square roots of rationals (on morphe's integers); the quantum 6j symbol at q = e^{iπ/r} |
| `spacetime` | Spacetimes as objects: Minkowski space and its causal structure; exact solutions as metric fields in charts; causal structure and horizons; discrete spacetimes (Regge geometry, causal sets, causal dynamical triangulations) | Minkowski space, causal relations and proper time; a catalogue of exact solutions (Minkowski, Schwarzschild in two charts, Reissner–Nordström, Kerr, FLRW for any curvature and scale factor, de Sitter, anti-de Sitter); Regge geometry of simplices and complexes (volumes by Cayley–Menger, dihedral angles, deficit angles, the Regge action); causal sets sprinkled into Minkowski diamonds, with their relations, links, chains and dimension estimators |
| `field_theory` | Classical fields (Lagrangians, Euler–Lagrange equations, Noether's theorem); lattice field theory on complexes and lattices; Monte Carlo and its statistics; perturbation theory (Feynman graphs, symmetry factors, integrals); free fields and propagators; renormalisation group; conformal and topological field theories; fields on curved spacetimes | scalar fields on the vertices of a complex (free and φ⁴), their exact free propagator, Metropolis and heat-bath updates; Markov-chain statistics (binning, autocorrelation time); vacuum Feynman graphs of φ^k theories with their symmetry factors, and zero-dimensional field theory as the check of the expansion |
| `gauge` | Connections and curvature; lattice gauge theory on any 2-complex; Monte Carlo of gauge fields; observables; character and strong-coupling expansions; exact two-dimensional Yang–Mills; BF theory and flat connections; Chern–Simons (later) | gauge fields on a CW complex's edges with values in a finite group, U(1), SU(2) or SU(N); holonomies along attaching words; gauge transformations; the Wilson action; Metropolis, heat-bath and overrelaxation updates; Wilson and Polyakov loops; the topological charge in two dimensions; the character expansion by the Weyl integration formula; Migdal and Witten's exact partition function on closed surfaces; flat connections counted |
| `quantum_gravity` | Spin networks; spin foams (Ponzano–Regge, Turaev–Viro, Barrett–Crane, EPRL); group field theory and tensor models; matrix models and random surfaces; causal dynamical triangulations; causal-set dynamics | the Turaev–Viro state sum on a triangulated 3-manifold (Ponzano–Regge regularised at a root of unity); BF amplitudes of 2-complexes over finite groups; the Feynman graphs of coloured tensor models and group field theories (coloured graphs), their dual complexes, jackets and Gurau degree, melonic graphs; the gluings of polygons of a matrix model and the Harer–Zagier numbers |
| `electromagnetism` | Maxwell's equations, potentials and gauges, radiation from sources, fields in media; on complexes by discrete exterior calculus (Yee's scheme as its instance) | planned |

`geometry` is written in part: metric fields in charts, Christoffel symbols, the Riemann, Ricci, Einstein and Weyl tensors and the Kretschmann scalar, computed exactly to rounding by forward differentiation of the metric (calx decision 0005), and geodesics as dynamical systems. Rotations and the rest of its plan wait.

### 2.2 The Order of the Modules

```
units
  └─ symmetry                 groups and representations (calx linalg, random; morphe arithmetic)
       └─ geometry            metrics, curvature, geodesics (calx differentiation, integrators)
            └─ spacetime      exact solutions, causal structure; Regge geometry, causal sets (morphe topology)
                 └─ field_theory       fields on complexes; Monte Carlo; Feynman graphs (morphe combinatorics)
                      └─ gauge         gauge fields on 2-complexes (morphe topology, algebra)
                           └─ quantum_gravity    spin foams, GFT, tensor and matrix models
electromagnetism: geometry, discretisation, gauge
gravitation: geometry, mechanics, spacetime
```

*Rationale.* A gauge field is a field, so `gauge` stands on `field_theory`'s lattices and Monte Carlo; group field theory is a field theory on group manifolds whose Feynman amplitudes are lattice gauge theories (spin foams) on its Feynman 2-complexes, so it stands above `gauge`, beside the spin foams it generates. Lie groups move from `geometry` to `symmetry`, beneath it, since physics uses them as symmetries before it uses them as spaces; `geometry` takes its rotations from `symmetry`.

### 2.3 Where the Theories Meet the Complexes

- A scalar field is a 0-cochain on a complex; its kinetic term lives on the 1-cells.
- A lattice gauge field is a function from the 1-cells to a group; the holonomy of a 2-cell is the product along its attaching word, and its trace enters the action. Every 2-complex of morphe, not only hypercubic lattices, carries a gauge theory.
- A discrete spacetime is a simplicial complex with edge lengths (Regge) or a partial order (causal set).
- The Feynman graphs of a group field theory are coloured graphs whose bicoloured cycles bound the 2-cells of a 2-complex, the dual 2-skeleton of the pseudomanifold they encode; the amplitude is a BF or spin-foam state sum on that 2-complex.

### 2.4 Conventions

- **Metric signature** (−, +, +, +), with the timelike coordinate x⁰ = ct first; the Riemann tensor R^a_{bcd} = ∂_c Γ^a_{db} − ∂_d Γ^a_{cb} + Γ^a_{ce} Γ^e_{db} − Γ^a_{de} Γ^e_{cb}, Ricci R_{bd} = R^a_{bad} (Misner, Thorne and Wheeler 1973, the "+ + +" convention).
- **Units of spacetime models.** Coordinates are lengths in metres (time as ct); a mass enters as its gravitational length GM/c² and a charge as √(G/(4πε₀)) Q/c², in metres; the conversions are `units`' constants. These are SI quantities, so invariant 3 holds.
- **Dimensionless theories.** Lattice field theories, gauge groups, recoupling symbols and state sums are dimensionless (lattice units, a = 1; ħ = 1); physical scale enters only where a caller sets a lattice spacing or a mass, and is converted there.
- **Spins** are passed as twice their value (two_j), so that half-integers are exact integers.
- **Gauge transformations** act on a link from tail to head as U_e ↦ g_tail U_e g_head⁻¹, and a holonomy multiplies the link variables in the order of the attaching word.

These are recorded in `docs/glossary.md`, created with this decision.

### 2.5 The Scope of General Relativity

Open decision 4 is answered: physis keeps a general metric toolkit (`geometry`, `spacetime`), not only post-Newtonian corrections; post-Newtonian dynamics stays planned in `gravitation`, and Einstein's equations as an evolution problem (numerical relativity) with it.

## 3. Alternatives Considered

- **One module for all of theoretical physics beyond the classical** (fields, gauge and gravity together): it would grow without inner order; the six modules each have one kind of object.
- **Group field theory inside `field_theory`:** its amplitudes are gauge theories on 2-complexes, so it would need `gauge` above it; placing it with the spin foams keeps the order acyclic.
- **Lie groups inside `geometry`,** as the map of the mathematics first had them: the recoupling theory of SU(2) and the characters of compact groups are representation theory, not geometry, and `gauge` and `quantum` need them without metrics.
- **Spacetime models inside `gravitation`:** exact solutions, causal sets and Regge geometry are arenas for every field theory, and the dynamics of gravity is a separate concern.
- **Physical units inside lattice theories:** every lattice computation would carry a spacing it does not use; dimensionless at the core and converted at the boundary is the practice of the field (Montvay and Münster 1994, ch. 1).

## 4. Consequences

- physis's plan has twenty-one modules; its README, architecture and CLAUDE.md list them, their order and what is written. The infrastructure's map of the mathematics moves MSC 22 (Lie groups) to `symmetry` and names the new modules under MSC 70, 81 and 83.
- The new modules are header-only and checked on calx's `verification`, against closed forms (Kretschmann scalars, free propagators, Migdal–Witten partition functions, Harer–Zagier numbers, the Turaev–Viro invariant of the 3-sphere) and identities (orthogonality of 6j symbols, Biedenharn–Elliott, Gauss–Bonnet, Schläfli, gauge invariance).
- The Complex Atlas, an application of the infrastructure, shows complexes, their invariants and the fields on them (ADR 0039).
