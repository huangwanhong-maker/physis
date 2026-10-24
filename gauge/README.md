# gauge: Gauge Fields

*Status: draft (2026-10-24): lattice gauge theory on 2-complexes, with a coupling for each 2-cell and links fixed at will; its Monte Carlo, observables and exact results; the exact contraction of ℤ₂ theories with their derivatives ([decision 0007](../development/decisions/0007-exact-z2-contraction-and-finite-heat-bath.md)).* Connections and curvature in the continuum, strong-coupling expansions beyond two dimensions, fermions, BF and Chern–Simons theories as field theories, and gauge theories in four dimensions at scale are planned ([decision 0005](../development/decisions/0005-symmetry-spacetime-fields-and-gauge.md)).

## 1. Purpose

Tools for research on gauge symmetry: a gauge field on any CW complex of morphe (lattice tori, surfaces given by words, dual 2-complexes of triangulations, the 2-complexes of coloured graphs), with a finite group, ℤ_N, U(1), SU(2) or SU(N); its Monte Carlo; its gauge-invariant observables; and the exact results against which all of it is checked.

## 2. Concepts

- **A lattice gauge field** assigns to each 1-cell a group element, the transport from its tail to its head. The holonomy of a 2-cell is the transport along its attaching word from its base (in the word's order), its lattice curvature. Gauge transformations act by U_e ↦ g_tail U_e g_head⁻¹ and conjugate holonomies, so the Wilson action β Σ (1 − Re tr U_f/N) and the traces of closed transports (Wilson and Polyakov loops) are invariant.
- **Linearity of the local action.** When an edge occurs at most once in each word, the action is −β/N Re tr(U_e Σ_e) with the staple Σ_e; heat bath (von Mises for U(1); Creutz and Kennedy–Pendleton for SU(2); Cabibbo–Marinari for SU(N)) and overrelaxation apply. When it occurs twice in a word (a one-vertex surface), only Metropolis, which recomputes holonomies, does; the field knows which edges are which.
- **The character expansion** writes a plaquette's weight as Σ_R c_R χ_R, the coefficients by the Weyl integration formula. On a closed orientable surface, integrating the links gives Migdal and Witten's Z = Σ_R d_R^χ (c_R/d_R)^F; on a disc the plaquettes are independent and Wilson loops obey an exact area law.
- **Couplings and fixed links.** The action β Σ_f w_f (1 − Re tr U_f / N) may weigh each 2-cell by its own w_f, and links may be frozen: a gauge fixing, a boundary condition, an observation. Updates leave frozen links alone.
- **The heat bath of a finite group** draws a link among the |G| elements in proportion to their weights, each computed by setting it: exact on any complex and for any couplings, and aperiodic where ℤ₂'s Metropolis, whose one proposal is the flip, is not.
- **The exact contraction of a ℤ₂ theory** sums the free links along an order of the 2-cells, with states in the character basis of the live links (Drouffe and Zuber 1983). With no link fixed, a state lives on the ℤ/2 1-cycles of the cut graph: a superposition of ℤ₂ spin networks, of dimension 2^{b₁}, so the cost is that of the widest cut, as for any transfer matrix (Baxter 1982, ch. 7). Pairing forward and backward states gives every mean, two-point functions, cut laws and the weights of the spin networks along which the character expansion's surfaces cross a cut. The same pass over calx jets gives log Z with its gradient and Hessian in grouped couplings, exactly.
- **The U(1) topological charge** Q = (1/2π) Σ_f o_f arg U_f on a closed orientable surface (with the orientation morphe's surface classification gives) is an integer.
- **Flat connections** of a finite group are counted by homomorphisms of the fundamental group: Π_c |G|^{V_c − 1} |Hom(π_1, G)|.

## 3. Content

| Header | Content |
|---|---|
| `group.hpp` | The concepts `GaugeGroup`, `StapleGroup` and `FiniteGaugeGroup`; `ZN{n}`, `U1Group`, `SU2Group`, `SUNGroup{n}`, `FiniteGroupAdapter{&group, weight}` (identity, mul, inv, re_trace, random, near, class_angle, name; staples; order) |
| `field.hpp` | `GaugeField<G>(group, X)`: `links()`, `set`, `randomise(stream)`, `linear(e)`, `transport(path)`, `holonomy(f)`, `plaquette(f)`, `mean_plaquette()`, `action(β)`, `wilson_loop(path)`, `transform(g)`, `local_plaquettes(e)`, `staple(e)`; `set_couplings(w)`, `coupling(f)`, `link_coupling(e)`; `freeze(e, u)`, `thaw(e)`, `frozen(e)` |
| `updates.hpp` | `metropolis_link`, `metropolis_sweep(U, β, stream, eps, hits)`, `heat_bath_link` (U(1), SU(2), SU(N)), `heat_bath_sweep(U, β, stream)` (staple groups, and every finite group), `overrelaxation_sweep(U)` (U(1), SU(2)), `SweepStatistics` |
| `contraction.hpp` | `Z2Theory {links, plaquettes, beta, fixed}` and `Z2Theory::of(X, β, fixed)`; `contraction_schedule`; `contract(T, order, two_points, cuts, networks, limits)` → `Contraction {log_z, link, plaquette, two_point, cuts (products, the law of negatives, spin networks), width, largest}`; `log_z_derivatives<M, K>(T, order, of, parameters)` (calx jets); `log_z(T, order)`; `log_cosh`; `ContractionLimits`, `Intractable` |
| `observables.hpp` | `mean_plaquette`, `polyakov_loop`, `plaquettes`, `plaquette_angles`, `link_angles`, `topological_charge(U, orientation)` |
| `exact.hpp` | `expansion_u1(β)`, `expansion_su2(β)`, `expansion_zn(N, β)`; `migdal_witten(expansion, χ, F)`, `disc(expansion, F)`, `disc_wilson_loop(expansion, A)`, `free_plaquette_sun(N, β)`; `count_flat_connections(X, G)`, `count_flat_connections_by_enumeration(X, G)` |
| `run.hpp` | `GroupSpec {kind: "Z" \| "U1" \| "SU2" \| "SU", n}`, `RunOptions`, `RunResult`; `run(X, spec, options)` (a group chosen at run time), `exact_mean_plaquette(X, spec, β)` |

## 4. Checks

`tests/check_gauge.cpp` (22 checks): gauge invariance of the action, plaquettes and Wilson loops for ℤ_3, S_3, U(1), SU(2) and SU(3) on a lattice torus, the one-vertex genus-2 surface and the dual 2-complex of ∂Δ⁴; staples against the local action; the heat-bath samplers against their exact moments on both branches; Monte Carlo against Migdal–Witten for U(1) and SU(2) (heat bath, also at strong coupling), SU(2) on one-vertex surfaces (Metropolis, non-linear links) and ℤ_3, all within 4 errors, and the heat bath of finite groups for ℤ_3 and ℤ_2 against the character formula; the disc's plaquette and Wilson loops of areas 1, 3, 4 against the area law for U(1) and SU(2), and SU(3) (Cabibbo–Marinari) against the Weyl integral; ℤ_2 and ℤ_3 enumerated exactly against the character formula; the integer topological charge, and |Q| = 1 for uniform flux; flat S_3 and ℤ_2 connections through π_1 and by enumeration (486, 3888, 16, 4).

`tests/check_contraction.cpp` (22 checks): on a 3 × 3 torus with a coupling for each plaquette and three links fixed, the contraction against brute force (log Z, every ⟨h_f⟩ and ⟨u_e⟩, 170 two-point functions, the law at a cut); log Z against Migdal and Witten on the 4 × 4 torus and the surface of genus 2; the 8192 flat connections at strong coupling; the spin networks at every cut are ℤ/2 1-cycles; the jets' gradient as Σ ⟨h_f⟩ and their Hessian against differences; Metropolis and the heat bath with couplings and frozen links against the contraction; and a circle in no 2-cell, whose Wilson loop the heat bath finds to vanish while ℤ₂'s Metropolis never moves it.

## 5. Dependencies

physis `field_theory` (statistics) and `symmetry` (the groups, Haar measures, Weyl integration); morphe `topology` (complexes, surfaces, fundamental groups) and `algebra` (finite groups, homomorphisms); calx `random`, `quadrature` (through `symmetry`) and `differentiation` (the jets of the contraction).

## 6. References

- Wilson, K. G. (1974), Confinement of quarks, *Phys. Rev. D* 10, 2445–2459.
- Creutz, M. (1983), *Quarks, Gluons and Lattices*, Cambridge University Press; Creutz, M. (1980), Monte Carlo study of quantized SU(2) gauge theory, *Phys. Rev. D* 21, 2308–2315; Creutz, M. (1987), Overrelaxation and Monte Carlo simulation, *Phys. Rev. D* 36, 515–519.
- Montvay, I. and Münster, G. (1994), *Quantum Fields on a Lattice*, Cambridge University Press, ch. 3–4.
- Kennedy, A. D. and Pendleton, B. J. (1985), Improved heatbath method for Monte Carlo calculations in lattice gauge theories, *Phys. Lett. B* 156, 393–399.
- Cabibbo, N. and Marinari, E. (1982), A new method for updating SU(N) matrices in computer simulations of gauge theories, *Phys. Lett. B* 119, 387–390.
- Best, D. J. and Fisher, N. I. (1979), Efficient simulation of the von Mises distribution, *Appl. Statist.* 28, 152–157.
- Drouffe, J.-M. and Zuber, J.-B. (1983), Strong coupling and mean field methods in lattice gauge theories, *Phys. Rep.* 102, 1–119.
- Baxter, R. J. (1982), *Exactly Solved Models in Statistical Mechanics*, Academic Press, ch. 7 (transfer matrices).
- Savit, R. (1980), Duality in field theory and statistical systems, *Rev. Mod. Phys.* 52, 453–487 (ℤ₂ gauge theory, its character expansion and its duality).
- Migdal, A. A. (1975), Recursion equations in gauge field theories, *Sov. Phys. JETP* 42, 413–418; Witten, E. (1991), On quantum gauge theories in two dimensions, *Commun. Math. Phys.* 141, 153–209.
- Lüscher, M. (1982), Topology of lattice gauge fields, *Commun. Math. Phys.* 85, 39–48.
- Dijkgraaf, R. and Witten, E. (1990), Topological gauge theories and group cohomology, *Commun. Math. Phys.* 129, 393–429.

## 7. Open Questions

1. **Four-dimensional lattices at scale**: the update loops are general (any complex) and not yet fast; a hypercubic special case, parallel checkerboard sweeps on calx's scheduler, and the computation stack's devices.
2. **SU(N) beyond the disc** in the exact checks: SU(3) characters for Migdal–Witten on closed surfaces.
3. **Fermions and the continuum**: Wilson and staggered fermions, and connections on principal bundles in physis `geometry`.
4. **Contraction beyond ℤ₂ and beyond a given order**: other abelian groups in the character basis, non-abelian groups through intertwiners, and orders found by a tree decomposition of the dual graph instead of the caller's.
