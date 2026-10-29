# continuum

*Status: draft (2026-10-29), [decision 0012](../development/decisions/0012-fields-as-sections.md).* Balance laws of continuous media on cellulations, in conservation form. The first is heat conduction. The Navier–Stokes equations, geophysical fluid dynamics (rotating and stratified flow, shallow-water and primitive equations), waves and transport are planned.

## 1. Concepts

A balance law says that a conserved quantity $u$, an $n$-form, changes by the divergence of its flux $J$, an $(n-1)$-form, and by its sources: $\partial_tu+dJ=s$. A constitutive relation closes it, such as Fourier's law $J=-k\star dT$ with $u=\rho cT$. On a cellulation the state is the quantity on each top cell and the flux is a number on each facet. Each interior flux enters two cells with opposite signs, so the discrete law conserves exactly what the continuous one does.

The law, the time integrator and the field are separate:

- the law gives the fluxes and the rate;
- the integrator, calx's θ-method, advances the state;
- the result is published as committed versions of a cochain field (physis `fields`, `discretisation`).

## 2. Contents

| Header | Content |
|---|---|
| `heat.hpp` | `Material` (conductivity, density, specific heat), `Boundary` (insulated, temperature, flux, exchange), `Budget`; `HeatConduction` on a box: the two-point flux with half-cells' conductances in series, the conductance Laplacian, `flux`, `rate`, `outflow`, `step` (a θ-step, then the energy updated by the fluxes), `steady`, `explicit_step_bound` (Gershgorin's), `quadratic` (½ Σ C T²) |

## 3. Checks

`tests/check_heat.cpp` (17 checks):

- **Exact states.** Two cells against their closed form. A linear profile, and a two-material profile with continuous flux, steady and exact to $10^{-13}$.
- **Budgets.** With a solve of tolerance $10^{-6}$, the budget closes to $10^{-17}$ of the stored energy. Each material's budget, with the heat through their interface, closes too. An insulated body conserves its energy.
- **Stability.** $\tfrac12\sum CT^2$ does not increase at θ ≥ ½. Gershgorin's step is stable, and the explicit scheme just above $2/\lambda_\text{max}$ is not.
- **Convergence** to the decaying sine mode: second order in space, second in time at θ = ½ and first at θ = 1.
- **As a field.** The temperature is a cochain field with committed versions.

## 4. Dependencies

`discretisation` (geometric complexes, boxes, cochain fields), `fields` (versions); calx `integrators` (the θ-method) and `linalg` (sparse rows, conjugate gradients).

## 5. References

- Fourier, J. (1822), *Théorie analytique de la chaleur*, Firmin Didot, Paris.
- Desbrun, M., Hirani, A. N., Leok, M. and Marsden, J. E. (2005), Discrete exterior calculus, arXiv:math/0508341.
- Eymard, R., Gallouët, T. and Herbin, R. (2000), Finite volume methods, in *Handbook of Numerical Analysis* VII, 713–1020 (the two-point flux and its convergence).
- Hairer, E. and Wanner, G. (1996), *Solving Ordinary Differential Equations II*, 2nd ed., Springer.

## 6. Open Questions

1. **Cellulations other than boxes**, where the line between centroids is not along the conormal and the two-point flux is not consistent: a multipoint flux, or a circumcentric dual.
2. **Nonlinear laws** (a conductivity that depends on temperature; advection): a Newton iteration in the step.
3. **The coupling of laws across interfaces** (formalism Proposition 3.4), where a Robin condition joins two media's own laws.
