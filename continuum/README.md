# continuum

*Status: draft (2026-10-29), [decision 0012](../development/decisions/0012-fields-as-sections.md); kinematic waves and transport on graphs (2026-10-31), [decision 0016](../development/decisions/0016-kinematic-waves-master-equations-and-reaction-diffusion.md).* Balance laws of continuous media on cellulations, in conservation form. The first is heat conduction; then kinematic waves on networks of links, and the advection and diffusion of a scalar on a graph. The Navier–Stokes equations, geophysical fluid dynamics (rotating and stratified flow, shallow-water and primitive equations) and other waves are planned.

## 1. Concepts

A balance law says that a conserved quantity $u$, an $n$-form, changes by the divergence of its flux $J$, an $(n-1)$-form, and by its sources: $\partial_tu+dJ=s$. A constitutive relation closes it, such as Fourier's law $J=-k\star dT$ with $u=\rho cT$. On a cellulation the state is the quantity on each top cell and the flux is a number on each facet. Each interior flux enters two cells with opposite signs, so the discrete law conserves exactly what the continuous one does.

The law, the time integrator and the field are separate:

- the law gives the fluxes and the rate;
- the integrator, calx's θ-method, advances the state;
- the result is published as committed versions of a cochain field (physis `fields`, `discretisation`).

The same form holds in one dimension on a network: a density $\rho$ per metre of link carried at a flow $q(\rho)$, $\partial_t\rho+\partial_xq(\rho)=0$ (Lighthill and Whitham 1955; Richards 1956). With a concave $q$, waves travel at $q'(\rho)$ and shocks at $[q]/[\rho]$; Godunov's flux between two cells is $\min(D(\rho_l),S(\rho_r))$, the demand of the cell behind and the supply of the cell ahead (Daganzo 1994; Lebacque 1996). At the nodes, entries admit what the first cell can take, exits take what their supply allows, and junctions share each outgoing link's supply among the demands sent to it by turning fractions, with or without first-in-first-out (Daganzo 1995; Tampère et al. 2011). On a graph whose vertices are cells, a scalar is exchanged across edges by diffusion through conductances and by upwinded advection, again in flux form, so that its content is conserved to rounding.

## 2. Contents

| Header | Content |
|---|---|
| `heat.hpp` | `Material` (conductivity, density, specific heat), `Boundary` (insulated, temperature, flux, exchange), `Budget`; `HeatConduction` on a box: the two-point flux with half-cells' conductances in series, the conductance Laplacian, `flux`, `rate`, `outflow`, `step` (a θ-step, then the energy updated by the fluxes), `steady`, `explicit_step_bound` (Gershgorin's), `quadratic` (½ Σ C T²) |
| `kinematic_waves.hpp` | `Flux` (Greenshields's and the triangular law: `q`, `dq`, `critical`, `capacity`, `demand`, `supply`, `max_speed`), `godunov_flux`, `riemann` (the entropy solution in closed form); `Link`, `KinematicNetwork` (entries, exits, junctions with turning fractions under non-FIFO or FIFO rules, scheduled laws; `start`, `step` with its `Budget`, `cfl`, `link_flows`, `content`) |
| `transport.hpp` | `GraphTransport`: volumes, edges, conductances, advective fluxes; `flux`, `rate`, `stable_step` (the explicit bound of positivity), `step` (Heun, strong-stability preserving), `content`, `matrix`, `steady` (dense, with its content) |

## 3. Checks

`tests/check_heat.cpp` (17 checks):

- **Exact states.** Two cells against their closed form. A linear profile, and a two-material profile with continuous flux, steady and exact to $10^{-13}$.
- **Budgets.** With a solve of tolerance $10^{-6}$, the budget closes to $10^{-17}$ of the stored energy. Each material's budget, with the heat through their interface, closes too. An insulated body conserves its energy.
- **Stability.** $\tfrac12\sum CT^2$ does not increase at θ ≥ ½. Gershgorin's step is stable, and the explicit scheme just above $2/\lambda_\text{max}$ is not.
- **Convergence** to the decaying sine mode: second order in space, second in time at θ = ½ and first at θ = 1.
- **As a field.** The temperature is a cochain field with committed versions.

`tests/check_kinematic_waves.cpp` (24 checks):

- **The Riemann problem.** Godunov's flux equals q of the Riemann solution at ξ = 0 for 40000 random pairs of both laws; Greenshields's shocks at $v(1-(\rho_l+\rho_r)/\rho_j)$, Lax's entropy condition, the fan's density of wave speed ξ, the triangular fan at $\rho_c$; the integral form $\int\rho\,d\xi=\Xi(\rho_l+\rho_r)+q(\rho_l)-q(\rho_r)$.
- **Convergence** on one link: order 1 in L1 at a shock (time-averaged, as a discrete shock's profile cycles), and $\Delta x\log(\ell/\Delta x)$ at a centred rarefaction (Harabetian 1988).
- **Networks.** A merge, a FIFO diverge and a choked exit: the budget closes to $10^{-14}$ of what passed, densities stay in $[0,\rho_j]$, and demand is refused once the queue reaches an entry. FIFO keeps turning fractions exactly; non-FIFO lets a free movement pass; merges share in proportion to demand; random junctions never exceed a supply or a demand.
- **Spillback** behind a bottleneck at the Rankine–Hugoniot speed (relative error $5\times10^{-5}$); a scheduled opening discharging the queue at the new capacity; a scheduled faster law tightening the CFL bound.

`tests/check_transport.cpp` (15 checks): content kept to rounding under any advection; a uniform state steady exactly under a circulation; the discrete heat kernel $e^{-2gt}I_n(2gt)$ reached at order 2 in time, with the variance $2gt$ exact; a pulse's mean moving at $u/V$; positivity at the explicit bound and its loss beyond; the maximum principle under a circulation; the steady state against the linear system and a long integration.

## 4. Dependencies

`discretisation` (geometric complexes, boxes, cochain fields), `fields` (versions); calx `integrators` (the θ-method) and `linalg` (sparse rows, conjugate gradients; dense solves for the steady state on a graph).

## 5. References

- Fourier, J. (1822), *Théorie analytique de la chaleur*, Firmin Didot, Paris.
- Desbrun, M., Hirani, A. N., Leok, M. and Marsden, J. E. (2005), Discrete exterior calculus, arXiv:math/0508341.
- Eymard, R., Gallouët, T. and Herbin, R. (2000), Finite volume methods, in *Handbook of Numerical Analysis* VII, 713–1020 (the two-point flux and its convergence).
- Hairer, E. and Wanner, G. (1996), *Solving Ordinary Differential Equations II*, 2nd ed., Springer.
- Daganzo, C. F. (1994), The cell transmission model, *Transp. Res. B* 28, 269–287; (1995), part II: network traffic, *Transp. Res. B* 29, 79–93.
- Greenshields, B. D. (1935), A study of traffic capacity, *Proc. Highway Res. Board* 14, 448–477.
- Harabetian, E. (1988), Rarefactions and large time behavior for parabolic equations and monotone schemes, *Commun. Math. Phys.* 114, 527–536.
- Jin, W.-L. and Zhang, H. M. (2003), On the distribution schemes for determining flows through a merge, *Transp. Res. B* 37, 521–540.
- Lebacque, J.-P. (1996), The Godunov scheme and what it means for first order traffic flow models, *Proc. 13th Int. Symp. Transportation and Traffic Theory*, 647–677.
- LeVeque, R. J. (2002), *Finite Volume Methods for Hyperbolic Problems*, Cambridge University Press.
- Lighthill, M. J. and Whitham, G. B. (1955), On kinematic waves II: a theory of traffic flow on long crowded roads, *Proc. R. Soc. A* 229, 317–345.
- Newell, G. F. (1993), A simplified theory of kinematic waves in highway traffic, *Transp. Res. B* 27, 281–313.
- Richards, P. I. (1956), Shock waves on the highway, *Oper. Res.* 4, 42–51.
- Shu, C.-W. and Osher, S. (1988), Efficient implementation of essentially non-oscillatory shock-capturing schemes, *J. Comput. Phys.* 77, 439–471.
- Tampère, C. M. J., Corthout, R., Cattrysse, D. and Immers, L. H. (2011), A generic class of first order node models for dynamic macroscopic simulation of traffic flows, *Transp. Res. B* 45, 289–309.

## 6. Open Questions

1. **Cellulations other than boxes**, where the line between centroids is not along the conormal and the two-point flux is not consistent: a multipoint flux, or a circumcentric dual.
2. **Nonlinear laws** (a conductivity that depends on temperature; advection): a Newton iteration in the step.
3. **The coupling of laws across interfaces** (formalism Proposition 3.4), where a Robin condition joins two media's own laws.
4. **Queues at entries.** An entry refuses what its link cannot take and counts it; a point queue (a vertical store before the link) would keep it, at the cost of a state beyond the cells.
5. **Higher order and other fluxes** for kinematic waves: a MUSCL reconstruction, and fluxes that are not concave (crowds at very high density), where Godunov's flux is no longer min(demand, supply).
