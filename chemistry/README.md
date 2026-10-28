# chemistry: Reaction Kinetics

*Status: draft (2026-10-29): reaction networks, mass action and their Markov chains (decision 0011).* Chemical thermodynamics and the coupling of reactions to transport are planned (README §2).

## 1. Purpose

The laws by which reactions proceed, stated so that a domain names only its species and its constants: which molecules exist, their formulas and rates belong to a chemistry domain; how any network of reactions moves its amounts, deterministically or molecule by molecule, belongs here.

## 2. Concepts

- **Networks and their stoichiometry.** Reactions Σ αᵢ Xᵢ → Σ βᵢ Xᵢ, the stoichiometric matrix N = β − α, and its left kernel: the conserved combinations, a lattice whose basis over ℤ comes from the Smith normal form (Horn and Jackson 1972; Feinberg 2019).
- **Balance.** Given each species' composition A and charge z, a network balances its elements and charge exactly when AN = 0 and zᵀN = 0.
- **Mass action.** dc/dt = N v(c), v_r = k_r Π cᵢ^{αᵢᵣ} (Guldberg and Waage 1864), integrated by calx.
- **The Markov chain of molecule numbers.** Reaction r fires at a_r(n) = c_r Π C(nᵢ, αᵢᵣ), with c_r = k_r Π αᵢᵣ! / Ω^{|α_r|−1}: occurrences, not ordered tuples (Gillespie 1976, 1977; Anderson and Kurtz 2015). On finite reachable states, the generator is written exactly, and the law of the absorbing state is solved exactly. Paths are sampled by Gillespie's direct method.
- **Symmetry.** A permutation of species that carries the reactions onto themselves commutes with the generator. A law invariant under it stays invariant, though each path need not be.

## 3. Content

| Header | Content |
|---|---|
| `kinetics.hpp` | `Reaction`, `Network`; `stoichiometric_matrix`, `conservation_laws`, `Balance`, `balance`; `mass_action`; `stochastic_constant`, `propensity`, `fire`, `drift`, `StateSpace`, `reachable`, `generator`, `absorption`, `gillespie`; `permute`, `image`, `symmetry_of` |

## 4. Checks

`tests/check_kinetics.cpp` (14 checks), on abstract networks:
- conservation laws and balance;
- mass action against closed forms;
- the neutral network's exact fixation probability i/N;
- a generator commuting with a symmetry of species;
- Dynkin's drift at every reachable state;
- the Poisson stationary law of a birth–death network sampled by Gillespie's method.

## 5. Dependencies

morphe `arithmetic` and `algebra` (rationals, the Smith form); calx `random`; calx `integrators` in the checks.

## 6. References

- Anderson, D. F. and Kurtz, T. G. (2015). *Stochastic Analysis of Biochemical Systems*. Springer.
- Feinberg, M. (2019). *Foundations of Chemical Reaction Network Theory*. Springer.
- Gillespie, D. T. (1976). A general method for numerically simulating the stochastic time evolution of coupled chemical reactions. *J. Comput. Phys.* 22, 403–434.
- Gillespie, D. T. (1977). Exact stochastic simulation of coupled chemical reactions. *J. Phys. Chem.* 81, 2340–2361.
- Horn, F. and Jackson, R. (1972). General mass action kinetics. *Arch. Rational Mech. Anal.* 47, 81–116.
