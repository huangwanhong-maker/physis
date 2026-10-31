# 0018. Open Systems Composed, What Coarse Descriptions Forget, Hybrid Processes, Memory, and Fields on Cells

- **Status:** Accepted, provisionally. The coupling bound is for the aggregation by the stationary law. Elimination is for linear systems. Hodge stars are diagonal.
- **Date:** 2026-10-31
- **Decided by:** the user, who asked for the formalism of their manuscript on urban landscapes to be implemented ("some module may resuable, like CW-complex, QFT, etc."); the infrastructure's [ADR 0062](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0062-landscapes-as-cellular-configurations.md). Claude wrote these modules and this record.

## 1. Context

A landscape whose stores, rules and passages change is a stochastic process assembled from parts. To be exact about it, physis needed eight things it lacked:
- the composition of parts at their interfaces;
- the test of whether a coarser description is a process of its own, and its error when it is not;
- a record of jumps that can be replayed;
- flows between jumps with jumps forced at a boundary;
- the memory left by eliminating variables;
- the identity between counting particles and following their master equation;
- the Laplacians of every degree on a complex with physical stars;
- the faces of a plane network.

## 2. Decision

### 2.1 Composition (`statistical_mechanics/composition.hpp`)

An open system has states, ports (functions of its state) and moves grouped into mechanisms. Two systems compose at a relation between their ports. The composite's states are the fibre product, the pairs whose ports are related. Local moves are lifted where they stay in it and refused at the interface where they would not. Joint mechanisms, declared on the components they read, are added once, at the first composition that holds all of them. With no ports the composite is the free product. Composition is associative with this rule. It differs from Baez, Fong and Pollard's gluing of boundary states (2016, *J. Math. Phys.* 57, 033301): here the parts keep their own state spaces, and the interface is a relation.

### 2.2 Lumpability (`statistical_mechanics/lumping.hpp`)

- **Tests.** Strong lumpability of a partition (Kemeny and Snell 1960, §6.3) and eventwise lumpability, mechanism by mechanism, each with the states where it fails and by how much.
- **The lumped generator.** The intertwining defect, and the path density of a coarse record.
- **Dwell laws.** The time a block lasts, a phase-type law.
- **When lumpability fails:**
  - `aggregated`, the aggregation by the stationary law, which keeps the blocks' stationary law;
  - `outward_discrepancy` δ, which bounds the total-variation distance between the projected and aggregated path laws over [0, T] by 1 − exp(−δT) (a maximal coupling; Lindvall 1992).

### 2.3 Jump Processes (`statistical_mechanics/jump_process.hpp`)

Gillespie's algorithm (1977) on occurrences, each with its mechanism, match, rate and successor. The record is anchored, so it holds the initial state and every event, and the state at any time can be read from it.

### 2.4 Piecewise-Deterministic Processes (`dynamics/piecewise_deterministic.hpp`)

Davis's processes (1984, *J. R. Stat. Soc. B* 46, 353). A flow, rates of spontaneous jumps, and a domain whose boundary forces a jump with a reset. The flow is integrated by RK4 augmented with the integrated rate Λ. Events are located by bisection at Λ = −ln U and at the boundary. A forced jump is a boundary condition of the generator, not a rate.

### 2.5 Elimination (`dynamics/elimination.hpp`)

For ẋ = Ax + By and ẏ = Cx + Dy, the exact reduced equation of x: its memory kernel K(τ) = B e^{Dτ} C, the noise term B e^{Dt} y(0), the full matrix, and the Markovian limit A − BD⁻¹C (Zwanzig 1961; Mori 1965; for linear systems, variation of constants).

### 2.6 Doi–Peliti Operators (`chemistry/doi_peliti.hpp`)

Creation and annihilation operators on a truncated Fock space. The operator of a reaction network, with κ = c_r/Πα!, equals the generator of its master equation on the counts reachable within the truncation (Doi 1976; Peliti 1985). Nothing quantum is in it.

### 2.7 Hodge Laplacians (`discretisation/hodge_laplacian.hpp`)

Δ_p on the cochains of any morphe chain complex, with diagonal stars. Its spectrum, eigenvectors orthonormal for the weighted inner product, and harmonic space, whose dimension is b_p (Eckmann 1944; Lim 2020). Constitutive laws are Hodge stars: heat in rooms is −Δ₃ with stars 1/G on faces and 1/C on rooms; walkers on a network are −Δ₀ with stars C on junctions and G on segments.

### 2.8 Plane Faces (`computational_geometry/planar.hpp`)

The faces of a plane graph traced from its rotation system, with V − E + F = 2c over the face orbits. Which faces are bounded by a simple circuit, their areas, and the polyline helpers: length, point at a distance, split, nearest point, crossing of segments.

## 3. Alternatives Considered

- **Composition as Baez–Fong–Pollard's cospans.** It identifies boundary states, which cannot say that a lane stands only under a permit; a relation on component states can.
- **Exact lumpability alone.** It answers yes or no. The coupling bound says how long an inexact description can be trusted, and the dwell law what it remembers.
- **Elimination by simulation of the full system.** It cannot separate memory from noise; the reduced equation does, exactly.

## 4. Consequences

- **Checks:**
  - `check_lumping` (17), `check_composition` (11), `check_jump_process` (6);
  - `check_piecewise_deterministic` (9), `check_elimination` (5);
  - `check_doi_peliti` (11), `check_hodge_laplacian` (10), `check_planar` (17).
- The urban-science domain's `landscape` module stands on these (the infrastructure's ADR 0062). Its memory of the shops is computed by calx's θ-method on the detailed-balance form of the same system, and checked against this elimination.
