# 0022. The Laws beneath Physical Computing: Chaos and Bifurcations, Waves and Rotating Flow, Networks of Conductances, Fluctuations, and Kinetics of Any Form

- **Status:** Accepted, provisionally. (0021 is left to the fluids of the vessels, which the infrastructure's ADR 0070 named.)
- **Date:** 2026-11-06
- **Decided by:**
  - **The request.** The user asked the infrastructure to "begin implement formalism and simulators" of its physical-computing domain ([ADR 0076](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0076-physical-computing.md)), the fluids, chemical and molecular families first, with "extenable, revisable, resuable" "the most important".
  - **What Claude did.** By the rule that a general method is written once in the foundation where it belongs, everything those simulators need that names no object in nature was written here first, with its checks, and this record.

## 1. Context

The domain's simulators model substrates: a tank, a convection cell, a rotating annulus, a BZ medium, a DNA circuit. What they stand on is general: the Lyapunov spectrum of a flow, the Floquet multipliers of a periodic system, the onset of convection, two-layer quasi-geostrophic flow, Kirchhoff's laws on a graph, the Potts model, the statistics of currents, reaction–diffusion with kinetics that are not mass action, and the compilation of polynomial ODEs into reaction networks. None of it was in physis, and every domain of the Earth and of life will need some of it.

## 2. Decision

| Module | Header | Content | Checks |
|---|---|---|---|
| `dynamics` | `flows.hpp` | Flows of run-time dimension (RK4, `integrate`), Jacobians by central differences, Lyapunov spectra of flows and maps by Benettin's method, Kaplan–Yorke's dimension; test systems with known invariants: Lorenz 63 and 96, Rössler, Hénon, the logistic map, Rosenzweig and MacArthur's predator–prey with its closed-form Hopf point | `check_chaos`: Lorenz's sum −41/3, its λ₁ 0.9056, Hénon's ln 0.3, the logistic map's ln 2 |
| `dynamics` | `floquet.hpp` | Monodromy matrices and Floquet multipliers; the damped parametric oscillator | `check_chaos`: Liouville's formula, the Mathieu tongue to O(ε²), the damped tip 4γ/ω₀ |
| `dynamics` | `continuation.hpp` | Pseudo-arclength continuation of equilibria (Keller), with folds and Hopf points located by bisection | `check_bifurcation`: x' = μ − x², the cusp's ±2/(3√3), the Brusselator's b = 1 + a² |
| `dynamics` | `kuramoto.hpp` | Kuramoto's networks, all to all or on a graph; the order parameter; Lorentzian frequencies by quantiles | `check_bifurcation`: r = √(1 − 2γ/K) |
| `continuum` | `waves.hpp` | Gravity–capillary dispersion, viscous damping, a basin of the damped wave equation with mirror walls | `check_waves`: exact grid modes, the leapfrog's conserved energy |
| `continuum` | `convection.hpp` | The marginal Rayleigh number by shooting with the exact propagator; critical values; Lorenz's truncation | `check_convection`: 1707.762, 1100.65, 27π⁴/4 |
| `continuum` | `two_layer.hpp` | Two-layer quasi-geostrophic flow, pseudo-spectral, with a channel by images | `check_two_layer`: Phillips's growth rates, energy and enstrophy kept |
| `continuum` | `capillarity.hpp` | Young–Laplace, Hagen–Poiseuille, Washburn | `check_capillarity` |
| `discretisation` | `kirchhoff.hpp` | Networks of conductances on graphs with fixed potentials and injections, dense or by conjugate gradients | `check_kirchhoff`: series, parallel, Wheatstone, Thomson |
| `statistical_mechanics` | `potts.hpp` | The q-state Potts model on a torus, Metropolis with uniform or neighbour proposals, exact enumeration, domains | `check_potts`: Metropolis against enumeration |
| `statistical_mechanics` | `fluctuations.hpp` | Current statistics by the tilted generator, the thermodynamic uncertainty relation, Landauer's bound | `check_fluctuations`: the biased ring, Gillespie, 200 random chains |
| `units` | `constants.hpp` | The seven exact defining constants of the SI | in `check_fluctuations` |
| `chemistry` | `field_kinetics.hpp` | Reaction–diffusion with kinetics of any form, vertex by vertex, on a graph, a torus or a grid with no-flux walls (a dish); a Strang-split step implicit in the kinetics (Alexander's SDIRK, damped Newton, recursive halving, non-negative species); `integrate_stiff` for well-mixed stiff networks | `check_field_kinetics`: diffusion's decay, Fisher–KPP with Bramson's correction, the split step's order and stability, the dish's walls keeping the total |
| `chemistry` | `dual_rail.hpp` | Polynomial ODEs compiled to networks of mass action by dual rails, exactly | `check_dual_rail`: d(p − q)/dt = P(p − q), Lorenz carried |

Three choices made on the way, each found by a check:

1. **Kinetics with poles need constraints.** The Oregonator's term (u − q)/(u + q) has a pole at u = −q. An implicit stage solved by plain Newton converged to a root beyond it, which the chemistry forbids. The split step's Newton is therefore damped, its failures halve the step, and species declared concentrations keep their iterates non-negative.
2. **A channel by images.** A periodic box admits modes uniform across a channel, which walls forbid and which can grow fastest. The two-layer model keeps the state odd in y on request, which is exactly a channel between free-slip walls.
3. **Pulled fronts converge slowly.** A Fisher front's speed approaches 2√(Dk) only logarithmically. The checks compare it with Bramson's correction, not with the limit.
4. **A torus can mimic a spiral.** On a periodic grid a plane wave can circulate for ever, so a sustained wave there shows nothing about rotation. Media that must not wrap are given no-flux walls (`SparseLaplacian::grid`), as a dish has.

## 3. Alternatives Considered

- **These methods in the domain.** Rejected: they name no object in nature.
- **Calx's integrators of fixed size.** Kept for small fixed systems. Run-time sizes need `flows.hpp`; chemistry's `runge_kutta_4` stays provisional, as written.
- **A general stiff integrator in calx.** Deferred: the split step serves kinetics on graphs, which is where the stiffness has appeared. A general Rosenbrock or BDF method belongs to calx when a second user appears.

## 4. Consequences

- physis gains 12 headers (one of them `units`'s constants) and 11 check programs, 124 checks in all; its whole suite passes, 1176 checks.
- The planned items of `dynamics` (Lyapunov exponents, bifurcation) and `statistical_mechanics` (the Potts model) are written in part.
- Open: the first Lyapunov coefficient of a Hopf point; branch switching at branch points; adaptive steps for the flows.
