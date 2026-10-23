# gauge: Gauge Fields

*Status: draft (2026-10-23): lattice gauge theory on 2-complexes, its Monte Carlo, observables and exact results.* Connections and curvature in the continuum, strong-coupling expansions beyond two dimensions, fermions, BF and Chern–Simons theories as field theories, and gauge theories in four dimensions at scale are planned ([decision 0005](../development/decisions/0005-symmetry-spacetime-fields-and-gauge.md)).

## 1. Purpose

Tools for research on gauge symmetry: a gauge field on any CW complex of morphe (lattice tori, surfaces given by words, dual 2-complexes of triangulations, the 2-complexes of coloured graphs), with a finite group, ℤ_N, U(1), SU(2) or SU(N); its Monte Carlo; its gauge-invariant observables; and the exact results against which all of it is checked.

## 2. Concepts

- **A lattice gauge field** assigns to each 1-cell a group element, the transport from its tail to its head. The holonomy of a 2-cell is the transport along its attaching word from its base (in the word's order), its lattice curvature. Gauge transformations act by U_e ↦ g_tail U_e g_head⁻¹ and conjugate holonomies, so the Wilson action β Σ (1 − Re tr U_f/N) and the traces of closed transports (Wilson and Polyakov loops) are invariant.
- **Linearity of the local action.** When an edge occurs at most once in each word, the action is −β/N Re tr(U_e Σ_e) with the staple Σ_e; heat bath (von Mises for U(1); Creutz and Kennedy–Pendleton for SU(2); Cabibbo–Marinari for SU(N)) and overrelaxation apply. When it occurs twice in a word (a one-vertex surface), only Metropolis, which recomputes holonomies, does; the field knows which edges are which.
- **The character expansion** writes a plaquette's weight as Σ_R c_R χ_R, the coefficients by the Weyl integration formula. On a closed orientable surface, integrating the links gives Migdal and Witten's Z = Σ_R d_R^χ (c_R/d_R)^F; on a disc the plaquettes are independent and Wilson loops obey an exact area law.
- **The U(1) topological charge** Q = (1/2π) Σ_f o_f arg U_f on a closed orientable surface (with the orientation morphe's surface classification gives) is an integer.
- **Flat connections** of a finite group are counted by homomorphisms of the fundamental group: Π_c |G|^{V_c − 1} |Hom(π_1, G)|.

## 3. Content

| Header | Content |
|---|---|
| `group.hpp` | The concepts `GaugeGroup` and `StapleGroup`; `ZN{n}`, `U1Group`, `SU2Group`, `SUNGroup{n}`, `FiniteGroupAdapter{&group, weight}` (identity, mul, inv, re_trace, random, near, class_angle, name; staples) |
| `field.hpp` | `GaugeField<G>(group, X)`: `links()`, `set`, `randomise(stream)`, `linear(e)`, `transport(path)`, `holonomy(f)`, `plaquette(f)`, `mean_plaquette()`, `action(β)`, `wilson_loop(path)`, `transform(g)`, `local_plaquettes(e)`, `staple(e)` |
| `updates.hpp` | `metropolis_link`, `metropolis_sweep(U, β, stream, eps, hits)`, `heat_bath_link` (U(1), SU(2), SU(N)), `heat_bath_sweep(U, β, stream)`, `overrelaxation_sweep(U)` (U(1), SU(2)), `SweepStatistics` |
| `observables.hpp` | `mean_plaquette`, `polyakov_loop`, `plaquettes`, `plaquette_angles`, `link_angles`, `topological_charge(U, orientation)` |
| `exact.hpp` | `expansion_u1(β)`, `expansion_su2(β)`, `expansion_zn(N, β)`; `migdal_witten(expansion, χ, F)`, `disc(expansion, F)`, `disc_wilson_loop(expansion, A)`, `free_plaquette_sun(N, β)`; `count_flat_connections(X, G)`, `count_flat_connections_by_enumeration(X, G)` |
| `run.hpp` | `GroupSpec {kind: "Z" \| "U1" \| "SU2" \| "SU", n}`, `RunOptions`, `RunResult`; `run(X, spec, options)` (a group chosen at run time), `exact_mean_plaquette(X, spec, β)` |

## 4. Checks

`tests/check_gauge.cpp` (20 checks): gauge invariance of the action, plaquettes and Wilson loops for ℤ_3, S_3, U(1), SU(2) and SU(3) on a lattice torus, the one-vertex genus-2 surface and the dual 2-complex of ∂Δ⁴; staples against the local action; the heat-bath samplers against their exact moments on both branches; Monte Carlo against Migdal–Witten for U(1) and SU(2) (heat bath, also at strong coupling), SU(2) on one-vertex surfaces (Metropolis, non-linear links) and ℤ_3, all within 4 errors; the disc's plaquette and Wilson loops of areas 1, 3, 4 against the area law for U(1) and SU(2), and SU(3) (Cabibbo–Marinari) against the Weyl integral; ℤ_2 and ℤ_3 enumerated exactly against the character formula; the integer topological charge, and |Q| = 1 for uniform flux; flat S_3 and ℤ_2 connections through π_1 and by enumeration (486, 3888, 16, 4).

## 5. Dependencies

physis `field_theory` (statistics) and `symmetry` (the groups, Haar measures, Weyl integration); morphe `topology` (complexes, surfaces, fundamental groups) and `algebra` (finite groups, homomorphisms); calx `random`, `quadrature` (through `symmetry`).

## 6. References

- Wilson, K. G. (1974), Confinement of quarks, *Phys. Rev. D* 10, 2445–2459.
- Creutz, M. (1983), *Quarks, Gluons and Lattices*, Cambridge University Press; Creutz, M. (1980), Monte Carlo study of quantized SU(2) gauge theory, *Phys. Rev. D* 21, 2308–2315; Creutz, M. (1987), Overrelaxation and Monte Carlo simulation, *Phys. Rev. D* 36, 515–519.
- Montvay, I. and Münster, G. (1994), *Quantum Fields on a Lattice*, Cambridge University Press, ch. 3–4.
- Kennedy, A. D. and Pendleton, B. J. (1985), Improved heatbath method for Monte Carlo calculations in lattice gauge theories, *Phys. Lett. B* 156, 393–399.
- Cabibbo, N. and Marinari, E. (1982), A new method for updating SU(N) matrices in computer simulations of gauge theories, *Phys. Lett. B* 119, 387–390.
- Best, D. J. and Fisher, N. I. (1979), Efficient simulation of the von Mises distribution, *Appl. Statist.* 28, 152–157.
- Drouffe, J.-M. and Zuber, J.-B. (1983), Strong coupling and mean field methods in lattice gauge theories, *Phys. Rep.* 102, 1–119.
- Migdal, A. A. (1975), Recursion equations in gauge field theories, *Sov. Phys. JETP* 42, 413–418; Witten, E. (1991), On quantum gauge theories in two dimensions, *Commun. Math. Phys.* 141, 153–209.
- Lüscher, M. (1982), Topology of lattice gauge fields, *Commun. Math. Phys.* 85, 39–48.
- Dijkgraaf, R. and Witten, E. (1990), Topological gauge theories and group cohomology, *Commun. Math. Phys.* 129, 393–429.

## 7. Open Questions

1. **Four-dimensional lattices at scale**: the update loops are general (any complex) and not yet fast; a hypercubic special case, parallel checkerboard sweeps on calx's scheduler, and the computation stack's devices.
2. **SU(N) beyond the disc** in the exact checks: SU(3) characters for Migdal–Witten on closed surfaces.
3. **Fermions and the continuum**: Wilson and staggered fermions, and connections on principal bundles in physis `geometry`.
