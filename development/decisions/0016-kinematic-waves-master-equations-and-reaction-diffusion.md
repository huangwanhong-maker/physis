# 0016. Kinematic Waves, Master Equations, and Reaction–Diffusion on Graphs

- **Status:** Accepted, provisionally (first-order kinematic waves with concave fluxes and entries that refuse rather than queue; dense master-equation routines for a few thousand states; a run-time Runge–Kutta method kept in `chemistry` until calx has one)
- **Date:** 2026-10-31
- **Decided by:** the user asked for streets to be analysed as an oceanographer analyses the sea: crowds as waves on street networks, movement among places as a process whose irreversibility is circulation, the mixing of uses along a street, and the coupling between holdings as reaction and diffusion, in the infrastructure's urban-science domain. Under the rule that general methods are written in the foundations, Claude wrote these methods and this record; the domain names the streets, the places and the rates.

## 1. Context

The street domain needs four general methods that name no object in nature:

- a conservation law whose flux saturates, solved on networks of links joined at nodes, for crowds that queue and spill back;
- Markov jump processes on places, with their stationary law, their entropy production and the currents round their cycles, their first passage and absorption;
- the entropy of mixing of kinds over places, and its decomposition over groups of places;
- reaction–diffusion on a graph: when diffusion destabilises a homogeneous state, and when a front between two states is held by the discreteness of the sites.

physis had heat conduction as its only balance law (decision 0012), mass action and the exact Markov chains of molecule numbers in `chemistry` (decision 0011), and `statistical_mechanics` only planned.

## 2. Decision

### 2.1 `continuum`: Kinematic Waves and Transport on Graphs

| Header | Content |
|---|---|
| `kinematic_waves.hpp` | `Flux` (Greenshields's and the triangular law), `godunov_flux` (min(demand, supply)), `riemann` (the entropy solution in closed form); `Link`, `KinematicNetwork`: cells on links, entries, exits, junctions with turning fractions under a non-FIFO rule (crowds) or FIFO (Daganzo 1995), laws scheduled in time, a step refused beyond the CFL bound, and a budget closing to rounding |
| `transport.hpp` | `GraphTransport`: a scalar on the vertices of a graph, exchanged by conductances and upwinded advective fluxes in flux form; Heun's step under the explicit bound of positivity; the steady state by a dense solve |

*Rationale.* Godunov's flux for a concave law is min(D(ρ_l), S(ρ_r)), the cell transmission model of traffic, so that one formula serves the interior of a link and its nodes: every rule at a node is a sharing of supplies among demands. The non-FIFO rule constrains each movement only by its own outgoing supply, as people who overtake do; FIFO blocks a whole incoming stream behind its most constrained branch, as lanes do. An entry refuses what its link cannot take and counts it, rather than queuing it, so that a state is its cells' densities alone.

### 2.2 `statistical_mechanics`: Its First Content

| Header | Content |
|---|---|
| `master_equation.hpp` | `generator`, `irreducible`, `stationary` (the GTH algorithm), `evolve` (calx's `expm`), `relative_entropy`; `irreversibility`: currents, entropy production (Schnakenberg eq. 7.6), a fundamental cycle basis from morphe's spanning forest with each cycle's affinity and current (eq. 7.11), Kolmogorov's criterion |
| `first_passage.hpp` | `PhaseType` (density, survival, moments, quantiles, the atom at zero), `phase_type`, `absorption_probabilities`, `mean_first_passage` |
| `mixing.hpp` | `mixing` (H(Q), H(Q\|X), the mutual information M, normalised and by place), `mixing_decomposition` (between and within groups, exactly) |

*Rationale.* GTH eliminates states without forming differences, so each π_i has a small relative error even when π spans tens of orders of magnitude, where replacing an equation of πQ = 0 by Σπ = 1 loses the small entries; the checks hold a birth–death chain over forty-six orders of magnitude to 2.4 × 10⁻¹⁴. Affinities are taken over a fundamental basis because they are additive over the cycle space: the basis decides detailed balance and decomposes σ, and has |E| − |V| + 1 members where enumerating every cycle would be exponential.

The module's planned ensembles, spin models and Monte Carlo stay planned; this content uses no other module of physis.

### 2.3 `chemistry`: Reaction–Diffusion on Graphs

| Header | Content |
|---|---|
| `reaction_diffusion.hpp` | `turing_modes` (the growth of each mode of a graph's Laplacian, by calx's symmetric eigenproblem and general eigenvalues); `turing_two_species` (the conditions, band, critical ratio and fastest k² in closed form); `path_laplacian`; `mckean_pinning`, `mckean_front`, `mckean_rate`; `runge_kutta_4` with run-time sizes; `ReactionDiffusion` (mass action of a `Network` at each vertex with diffusion on the graph) |

The pinning interval of McKean's caricature on a chain is derived in the header: the stationary front is a decaying solution of a linear recurrence on each side, the two meeting at the interface, which gives a_± = ½(1 ± 1/√(1 + 4D)). Keener (1987) and Fáth (1998) are cited for the phenomenon; the checks read the closed form off a direct integration that does not know it.

The fastest-growing k² of two species is found where dλ/dk² = 0, a quadratic in k² derived in the header, rather than at the minimum of h(k²), which is the onset wavenumber and differs from it away from onset.

## 3. Alternatives Considered

- **A general finite-volume solver for conservation laws** (any flux, any reconstruction). Not yet: a concave flux makes Godunov's flux a minimum of two monotone functions, and the street domain needs no other; a non-concave flux or a second-order reconstruction can join beside it.
- **Queues at entries**, a point store before each link. Rejected for now: refused demand is counted in the budget, and a queue would add a state beyond the cells. It may come when a domain needs to keep what waits.
- **Enumerating every cycle** for affinities. Rejected: exponential, and redundant, as the fundamental cycles span them.
- **Sparse master equations** (Krylov exponentials, iterative stationary solvers). Not yet: the places of a street are hundreds, within the dense routines' reach.
- **calx's fixed-size integrators** for reaction–diffusion. Not usable: a graph's size is known at run time. A run-time RK4 is written in `chemistry`, provisionally, to move down to calx with its integrators.

## 4. Consequences

- `check_physis_kinematic_waves` (24), `check_physis_transport` (15), `check_physis_master_equation` (25), `check_physis_first_passage` (15), `check_physis_mixing` (8) and `check_physis_reaction_diffusion` (21), against closed forms, brute force and Monte Carlo; no fixtures were needed.
- `chemistry` now links calx `linalg`; `statistical_mechanics` links calx `linalg` and morphe `combinatorics`, and is added to the root CMakeLists after `quantum`.
- The infrastructure's urban-science domain builds its streets on these methods: crowds as kinematic waves, movement among places as a master equation whose circulation is its irreversibility, the mixing of uses, and the persistence of regimes along a street as pinned fronts.
