# 0017. Bodies in a Fluid, and the Euler–Lagrange Equations Solved for the Acceleration

- **Status:** Accepted, provisionally (strip theory: the reactive force of each strip's added mass, and quadratic drag; no wake)
- **Date:** 2026-10-31
- **Decided by:** the user, who asked for three robots, "fish for ocean" among them; the infrastructure's [ADR 0059](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0059-robots-in-their-environments.md). Claude wrote these modules and this record.

## 1. Context

A swimming robot is a planar tree (decision 0014) in water. The water acts on it in two ways. The reactive force comes from the kinetic energy of the water it sets moving, its added mass. The resistive force is drag. The added mass belongs in the Lagrangian, as the fluid's kinetic energy, so that the variational integrator moves the body and the water together. But then the explicit methods cannot use the tree's Newton–Euler forward dynamics, which know only the body's mass matrix.

## 2. Decision

### 2.1 The Euler–Lagrange Equations Solved for the Acceleration (`mechanics/lagrangian.hpp`)

For any Lagrangian $L(q, v)$ written over its number type and a generalized force $Q$,
$$\frac{\partial^2L}{\partial v\,\partial v}\,a = Q + \frac{\partial L}{\partial q} - \frac{\partial^2L}{\partial v\,\partial q}\,v.$$
One evaluation on calx's jets of order 2 in the $2N$ arguments gives the gradient and the Hessian's blocks, and calx's small dense solve gives the acceleration. Semi-implicit Euler and Runge–Kutta are written on it. For a tree alone it is Newton–Euler's acceleration, to $6\times10^{-16}$.

### 2.2 A Tree in a Fluid (`mechanics/fluid.hpp`)

Strip theory ([theory](../../docs/theory/discrete-mechanics.md) §8). Each link carries a slender part of a length and a height, cut into strips. Each strip moves relative to the fluid, a uniform current $u$, with components $w_n$ across the part and $w_t$ along it:
- **the reactive force:** the fluid's kinetic energy $T_f = \sum \tfrac12(\rho\pi h^2/4)\,ds\,w_n^2$, the added mass of a plate per unit length (Newman 1977, §4.13; Lighthill 1960), a term of the Lagrangian (`TreeInFluid`);
- **the resistive force:** $f = -\tfrac12\rho C_n h\,ds\,|w_n|w_n\,n - \tfrac12\rho C_f(2h)\,ds\,|w_t|w_t\,e$ (Morison et al. 1950; Taylor 1952), its generalized force $J^{\mathsf T}f$ at each strip, $|w|$ regularised for jets.

## 3. Alternatives Considered

- **Added mass as a force** (the acceleration times $m_a$, moved to the left). Rejected: it puts the acceleration in a force, which the variational integrator cannot take. As the fluid's kinetic energy it is a term of the action, and its forces follow.
- **Newton–Euler with an added-mass matrix.** Rejected: the matrix's derivatives, the forces of its change, would be written by hand. Jets give them from the one Lagrangian.
- **Lighthill's large-amplitude elongated-body theory, or a vortex method.** Not yet: strip theory is the first model a swimming robot is designed with. The trailing edge's momentum flux, which gives a real fish most of its thrust, is open (theory §8).

## 4. Consequences

- `check_physis_fluid` (9) checks:
  - the Lagrangian's acceleration against Newton–Euler's;
  - a plate's added mass across and its absence along;
  - quadratic drag's closed forms across and along;
  - a plate carried toward a current's speed as quadratic drag allows;
  - the energy of the body and the water together kept without drag;
  - the variational integrator against Runge–Kutta within its order-2 error.
- The robotics domain's fish is built on them (ADR 0059).
