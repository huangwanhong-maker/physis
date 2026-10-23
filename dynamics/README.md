# dynamics: Dynamical Systems

*Status: draft (2026-10-23): linear propagators and hybrid systems.* Discrete-time maps, delay and stochastic dynamics, equilibria, stability, bifurcation, Lyapunov exponents and invariant measures are planned (the infrastructure's [map of the mathematics](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/docs/mathematics.md), §5).

## 1. Purpose

A dynamical system as an object: its state space, its flow, its events and what it keeps invariant. The methods that integrate a vector field knowing nothing of it are calx's (`integrators`); a field equation in space becomes a dynamical system through `discretisation`. Domains state their systems here: a neuron, a box model of the climate, an epidemic.

## 2. Concepts

- **Linear systems** x' = A x + B u + c are integrated exactly over a step by the exponential of the augmented matrix [[A, B, c], [0, 0, 0]] h (Van Loan 1978), whose blocks give e^{Ah} and the integrals of e^{As} against the input and the constant (Rotter and Diesmann 1999). Components may be held still. In SI units the components of one state can differ by ten orders of magnitude, so the matrix is made dimensionless by the components' typical sizes, rounded to powers of two, before it is exponentiated.
- **Hybrid systems** flow between events, fire when the state enters a guard set (tested on the grid), jump, and dwell for a while with some components held: the integrate-and-fire neuron and its refractory period are the type. A system may also declare a silent jump, `settle`, applied after each step's flow and emitting no event: a second jump set (Goebel, Sanfelice and Teel 2012), with which a guard is given hysteresis, so that an upward crossing is counted once.
- **Gated systems** declare, for each gate, its steady state and time constant, and are integrated by Rush and Larsen's exponential method.

## 3. Content

| Header | Content |
|---|---|
| `linear.hpp` | `Propagator<N, M>`, `propagator(A, B, c, h, held, sizes)` |
| `hybrid.hpp` | The concepts `System`, `LinearFlow`, `NonlinearFlow`, `Gated`, `Hybrid`, `Settling`; `linear_propagator<S>`; `step<S>` by Euler, Heun, RK4, Dormand–Prince or Rush–Larsen with held components; `simulate<S>`, one trajectory with its events |

## 4. Checks

`tests/check_dynamics.cpp`: the propagator of a membrane and its synaptic current against the closed form at equal and unequal time constants, and with held components; scaling by the components' sizes against none; the events of an integrate-and-fire system on the grid against its closed-form first passage; Rush and Larsen's step, exact for frozen coefficients and of order one otherwise.

## 5. Dependencies

physis `units` (conceptually); calx `linalg` and `integrators`.

## 6. References

- Van Loan, C. F. (1978), Computing integrals involving the matrix exponential, *IEEE Trans. Autom. Control* 23, 395–404.
- Rotter, S. and Diesmann, M. (1999), Exact digital simulation of time-invariant linear systems with applications to neuronal modeling, *Biol. Cybern.* 81, 381–402.
- Goebel, R., Sanfelice, R. G. and Teel, A. R. (2012), *Hybrid Dynamical Systems*, Princeton University Press.
- Rush, S. and Larsen, H. (1978), A practical algorithm for solving dynamic membrane equations, *IEEE Trans. Biomed. Eng.* 25, 389–392.

## 7. Open Questions

1. Events located between grid points (by interpolation or root finding) for systems whose event times matter more finely than the step.
2. Stochastic dynamics: the system here, the statistics in the inference library, the methods in calx.
