# chemistry: Reaction Kinetics and Reaction–Diffusion

*Status: draft (2026-10-29): reaction networks, mass action and their Markov chains (decision 0011); reaction–diffusion on graphs, Turing's instability and the pinning of fronts (2026-10-31, decision 0016).* Chemical thermodynamics, and the coupling of reactions to transport in continua, are planned (README §2).

## 1. Purpose

The laws by which reactions proceed, stated so that a domain names only its species and its constants: which molecules exist, their formulas and rates belong to a chemistry domain; how any network of reactions moves its amounts, deterministically or molecule by molecule, belongs here.

## 2. Concepts

- **Networks and their stoichiometry.** Reactions Σ αᵢ Xᵢ → Σ βᵢ Xᵢ, the stoichiometric matrix N = β − α, and its left kernel: the conserved combinations, a lattice whose basis over ℤ comes from the Smith normal form (Horn and Jackson 1972; Feinberg 2019).
- **Balance.** Given each species' composition A and charge z, a network balances its elements and charge exactly when AN = 0 and zᵀN = 0.
- **Mass action.** dc/dt = N v(c), v_r = k_r Π cᵢ^{αᵢᵣ} (Guldberg and Waage 1864), integrated by calx.
- **The Markov chain of molecule numbers.** Reaction r fires at a_r(n) = c_r Π C(nᵢ, αᵢᵣ), with c_r = k_r Π αᵢᵣ! / Ω^{|α_r|−1}: occurrences, not ordered tuples (Gillespie 1976, 1977; Anderson and Kurtz 2015). On finite reachable states, the generator is written exactly, and the law of the absorbing state is solved exactly. Paths are sampled by Gillespie's direct method.
- **Symmetry.** A permutation of species that carries the reactions onto themselves commutes with the generator. A law invariant under it stays invariant, though each path need not be.
- **Reaction–diffusion on graphs.** Species at the vertices of a graph react by mass action and diffuse along its edges, dc_v/dt = f(c_v) − Σ_w L_vw D c_w. About a homogeneous steady state, the mode of the Laplacian with eigenvalue Λ grows at the largest real part of eig(J − ΛD); diffusion destabilises a stable state when the species diffuse unequally enough (Turing 1952; on networks Othmer and Scriven 1971, Nakao and Mikhailov 2010). For two species the conditions, the unstable band, the critical ratio of diffusivities and the fastest mode are in closed form (Murray 2003, §2.3).
- **Fronts and their pinning.** Between the two stable states of a bistable medium on a chain of sites, a front moves unless the sites' discreteness holds it (Keener 1987; Fáth 1998). For McKean's piecewise-linear caricature the stationary front is μⁿ/(1 + μ) on one side and 1 − μ^{1−n}/(1 + μ) on the other, μ the decaying root of D μ² − (2D + 1) μ + D = 0, so that it is pinned exactly for thresholds a in (½(1 − 1/√(1 + 4D)), ½(1 + 1/√(1 + 4D))) (derived in `reaction_diffusion.hpp`).

## 3. Content

| Header | Content |
|---|---|
| `kinetics.hpp` | `Reaction`, `Network`; `stoichiometric_matrix`, `conservation_laws`, `Balance`, `balance`; `mass_action`; `stochastic_constant`, `propensity`, `fire`, `drift`, `StateSpace`, `reachable`, `generator`, `absorption`, `gillespie`; `permute`, `image`, `symmetry_of` |
| `reaction_diffusion.hpp` | `TuringModes`, `turing_growth`, `turing_modes` (on any graph's Laplacian); `TwoSpeciesTuring`, `turing_two_species`, `turing_two_species_growth` (closed forms); `path_laplacian`; `Pinning`, `mckean_pinning`, `mckean_front`, `mckean_rate`; `runge_kutta_4` (run-time sizes, provisional); `ReactionDiffusion` (mass action of a `Network` with diffusion on a graph: `rate`, `step`, `diffusive_step`) |

## 4. Checks

`tests/check_kinetics.cpp` (14 checks), on abstract networks:
- conservation laws and balance;
- mass action against closed forms;
- the neutral network's exact fixation probability i/N;
- a generator commuting with a symmetry of species;
- Dynkin's drift at every reachable state;
- the Poisson stationary law of a birth–death network sampled by Gillespie's method.

`tests/check_reaction_diffusion.cpp` (21 checks):
- the path Laplacian's eigenvalues in closed form;
- the two-species formulas against the Schnakenberg kinetics' own Turing condition in a, b and d, the exchange of species, the band's ends as roots, the fastest mode against a scan;
- each mode's growth on a chain equal to the dispersion relation at its eigenvalue, the unstable modes exactly those in the band;
- Schnakenberg's kinetics as a network on a chain of 54 sites settling into the predicted pattern (mode 5, three peaks);
- `rate` against `mass_action` and a dense L D c; the totals of a conserving network kept;
- McKean's front, integrated from a step, at rest on the closed form with its interface sites at a₊ and a₋; no flip just inside the interval; invasion just outside.

## 5. Dependencies

morphe `arithmetic` and `algebra` (rationals, the Smith form); calx `random`; calx `linalg` (the symmetric eigenproblem and general eigenvalues); calx `integrators` in the checks.

## 6. References

- Anderson, D. F. and Kurtz, T. G. (2015). *Stochastic Analysis of Biochemical Systems*. Springer.
- Feinberg, M. (2019). *Foundations of Chemical Reaction Network Theory*. Springer.
- Gillespie, D. T. (1976). A general method for numerically simulating the stochastic time evolution of coupled chemical reactions. *J. Comput. Phys.* 22, 403–434.
- Gillespie, D. T. (1977). Exact stochastic simulation of coupled chemical reactions. *J. Phys. Chem.* 81, 2340–2361.
- Horn, F. and Jackson, R. (1972). General mass action kinetics. *Arch. Rational Mech. Anal.* 47, 81–116.
- Fáth, G. (1998). Propagation failure of traveling waves in a discrete bistable medium. *Physica D* 116, 176–190.
- Keener, J. P. (1987). Propagation and its failure in coupled systems of discrete excitable cells. *SIAM J. Appl. Math.* 47, 556–572.
- McKean, H. P. (1970). Nagumo's equation. *Adv. Math.* 4, 209–223.
- Murray, J. D. (2003). *Mathematical Biology II: Spatial Models and Biomedical Applications*, 3rd ed. Springer.
- Nakao, H. and Mikhailov, A. S. (2010). Turing patterns in network-organized activator–inhibitor systems. *Nat. Phys.* 6, 544–550.
- Othmer, H. G. and Scriven, L. E. (1971). Instability and dynamic pattern in cellular networks. *J. Theor. Biol.* 32, 507–537.
- Schnakenberg, J. (1979). Simple chemical reaction systems with limit cycle behaviour. *J. Theor. Biol.* 81, 389–400.
- Turing, A. M. (1952). The chemical basis of morphogenesis. *Phil. Trans. R. Soc. B* 237, 37–72.

## 7. Open Questions

1. **Run-time integrators.** `runge_kutta_4` belongs in calx once its integrators take run-time sizes; implicit–explicit schemes would lift the diffusive step bound on fine graphs.
2. **Sparse Laplacians.** `turing_modes` is dense (O(n³)): large graphs need a few extreme eigenpairs by Lanczos.
3. **Smooth bistable kinetics.** The pinning interval is closed-form only for the piecewise-linear caricature; for cubic kinetics it is computed (Keener 1987), not yet here.
