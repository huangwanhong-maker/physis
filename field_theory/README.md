# field_theory: Fields and Their Quantisation

*Status: draft (2026-10-23): scalar fields on complexes, Markov-chain statistics, vacuum Feynman graphs, zero-dimensional field theory; invariant theories of every rank and their large-N solutions ([decision 0006](../development/decisions/0006-the-generation-of-spacetime-by-field-theories.md)).* Classical field theory (Lagrangians, Euler–Lagrange equations, Noether's theorem), free fields and propagators in the continuum, Feynman graphs with external legs and their integrals, the renormalisation group, conformal and topological field theories, and fields on curved spacetimes are planned ([decision 0005](../development/decisions/0005-symmetry-spacetime-fields-and-gauge.md)).

## 1. Purpose

The infrastructure of field theory on which `gauge` and `quantum_gravity` stand: fields as cochains on the cell complexes of morphe, their Euclidean path integrals sampled by Markov chains, the statistics that give those samples their errors, and the combinatorics of perturbation theory. Everything is in lattice units (a = 1, ħ = 1); physical scale enters where a caller sets a lattice spacing.

## 2. Concepts

- **A scalar field is a 0-cochain** on a CW complex. Its kinetic term ½ |dφ|² lives on the 1-cells, so the free action is ½ φᵀ(Δ + m²)φ with Δ the graph Laplacian of the 1-skeleton, on any complex, not only hypercubic lattices; φ⁴ adds (λ/4!) Σ φ⁴ on the vertices.
- **The free propagator** is (Δ + m²)⁻¹: by dense elimination on a small complex, and by the momentum sum on a periodic lattice, which must agree.
- **Markov-chain Monte Carlo** samples e^{−S}: Metropolis for any action, the exact Gaussian heat bath for the free field. Successive samples are correlated, so a mean's error is √(2 τ_int σ²/n), with τ_int estimated by Sokal's automatic window, or by binning; the jackknife carries errors through nonlinear functions (Binder's cumulant).
- **Feynman graphs** are the isomorphism classes of Wick pairings; each graph Γ arises (k!)^n n!/|Aut Γ| times, so Z/Z₀ = Σ_Γ (−λ)^n/|Aut Γ| and log Z sums the connected graphs. Pairings are grouped by labelled multigraph, with a closed-form multiplicity, and the classes found by canonical relabelling.
- **Zero-dimensional φ⁴** is the integral ∫ dφ e^{−φ²/2 − λφ⁴/4!}: its perturbation series has zero radius of convergence, and since Z is completely monotone in λ the error of each truncation is below the first omitted term, e^{−3/(2λ)} at the optimum. The exact series coefficients (morphe's rationals) and those of log Z check the graphs.

- **Invariant theories generate spacetime** ([decision 0006](../development/decisions/0006-the-generation-of-spacetime-by-field-theories.md)). A complex field of rank r (scalar, vector, matrix, tensor) with U(N)^r-invariant interactions, its bubbles: its Feynman expansion runs over all propagator matchings σ₀, each graph an (r+1)-coloured graph weighted N^{F₀}, with isomorphism classes, |Aut Γ| and exact series of Z and log Z in N (the matrix graphs are surfaces weighted N^{2−2g}, the tensor graphs pseudomanifolds weighted by Gurau's degree). The dominant graphs resum at large N: the vector model's gap equation; the eigenvalue density and planar free energy of the Hermitian quartic matrix model (Brézin, Itzykson, Parisi and Zuber 1978), whose critical point g_c = −1/48 (with V = M²/2 + gM⁴) is continuous two-dimensional quantum gravity, γ = −1/2; the melonic two-point function G = 1 + x G^m with its Fuss–Catalan series and γ = 1/2, branched polymers.

## 3. Content

| Header | Content |
|---|---|
| `scalar.hpp` | `ScalarTheory {m2, lambda}`; `neighbours(X)`, `laplacian(X)`, `free_propagator(X, m2)`, `lattice_coordinates`, `lattice_propagator(sides, m2, r)`; `action`, `local_action`; `metropolis_sweep(T, N, phi, stream, step)`, `heat_bath_sweep(T, N, phi, stream)`; `magnetisation`, `mean_square`, `binder_cumulant` |
| `statistics.hpp` | `Estimate {mean, error, tau, n}`; `mean`, `variance`, `autocorrelation_time(x, c = 6)`, `estimate(x)`, `binned(x, bins)`, `jackknife(series, f, bins)` |
| `feynman.hpp` | `FeynmanGraph {order, adjacency, pairings, automorphisms, connected, edges, loops, weight()}`; `labelled_vacuum_graphs(k, n, f)`, `vacuum_graphs(k, n)`, `wick_pairings(k, n, f)`, `graph_of` |
| `zero_dimensional.hpp` | `zero_dimensional_partition(λ, m2)`, `zero_dimensional_moment(2k, λ, m2)`, `partition_series(N)` (exact), `log_series(a)` (exact), `partial_sum` |
| `invariant.hpp` | `Interaction{name, bubble, coupling, scaling}`, `InvariantTheory{name, rank, interactions}`; `complex_scalar`, `vector_model`, `matrix_model`, `quartic_melonic_model`, `wheel_model`, `necklace_model`, `theory_by_name`; `expand(T, counts, classes)` → `Expansion{all, connected, classes: FeynmanClass{sigma0, graph, labelled, automorphisms, faces, scaled, components, connected, degree}}`; `series(T, order)` → `Series{z, free_energy, connected}`; `leading`, `evaluate`, `leading_power` |
| `large_n.hpp` | `vector_two_point`, `vector_critical`, `vector_free_energy`; `bipz_critical`, `bipz_a2`, `bipz_density`, `bipz_free_energy`, `bipz_series`, `bipz_moment`; `melonic_critical`, `melonic_two_point`, `melonic_coefficients`; `exponent_from_series` |

## 4. Checks

`tests/check_field_theory.cpp` (20 checks): the dense propagator against the momentum sum on 4×4, 3×4×5 and 2×6 tori, and the defining equation on the subdivided 7-vertex torus (42 vertices, degrees 4 to 12); the heat bath's two-point function and susceptibility against the propagator, and Binder's cumulant of the Gaussian field, within 4 errors; Metropolis against the propagator; φ⁴ on one point against the integral; the autocorrelation time of an AR(1) series against (1 + ρ)/(2(1 − ρ)); Wick pairings grouped with the closed-form multiplicity (φ³, φ⁴, φ⁶); the symmetry factors 1/8 and 1/16, 1/48, 1/128; Σ 1/|Aut| = (kn − 1)!!/((k!)^n n!); the connected φ⁴ graphs (1, 2, 4, 10, 28 of them) against log Z to order 5 (−1/8, 1/12, −11/96, 17/72, −619/960); every truncation of the series within its next term, and the optimal one within e^{−3/(2λ)}.

`tests/check_invariant.cpp` and `tests/check_large_n.cpp` ([decision 0006](../development/decisions/0006-the-generation-of-spacetime-by-field-theories.md)): Gaussian moments in every rank (m!, the rising factorial N(N+1)⋯, the permutation sums of complex matrices with Catalan leading terms); log Z equal to the sum over connected graphs; N^{2−2g} for matrices and Gurau's formula for tensors; the melonic leading free energies 2^{n−1}Cat_n/n and D(2D)^{n−1}Cat_n/n; ∫ρ = 1 and the moments of the eigenvalue density; the BIPZ series against its closed form and against planar ribbon-graph counts; the critical exponents −1/2 and 1/2 from 400 coefficients.

## 5. Dependencies

morphe `topology` (the complexes) and `arithmetic` (exact series); calx `linalg`, `random`, `quadrature`.

## 6. References

- Montvay, I. and Münster, G. (1994), *Quantum Fields on a Lattice*, Cambridge University Press, ch. 1–2.
- Desbrun, M., Hirani, A. N., Leok, M. and Marsden, J. E. (2005), Discrete exterior calculus, arXiv:math/0508341.
- Metropolis, N., Rosenbluth, A. W., Rosenbluth, M. N., Teller, A. H. and Teller, E. (1953), Equation of state calculations by fast computing machines, *J. Chem. Phys.* 21, 1087–1092.
- Madras, N. and Sokal, A. D. (1988), The pivot algorithm, *J. Stat. Phys.* 50, 109–186; Sokal, A. D. (1997), Monte Carlo methods in statistical mechanics, in *Functional Integration*, Plenum, 131–192.
- Flyvbjerg, H. and Petersen, H. G. (1989), Error estimates on averages of correlated data, *J. Chem. Phys.* 91, 461–466.
- Efron, B. (1982), *The Jackknife, the Bootstrap and Other Resampling Plans*, SIAM.
- Binder, K. (1981), Finite size scaling analysis of Ising model block distribution functions, *Z. Phys. B* 43, 119–140.
- Itzykson, C. and Zuber, J.-B. (1980), *Quantum Field Theory*, McGraw-Hill, §6-1.
- Zinn-Justin, J. (2002), *Quantum Field Theory and Critical Phenomena*, 4th ed., Oxford, §9.1.
- Dyson, F. J. (1952), Divergence of perturbation theory in quantum electrodynamics, *Phys. Rev.* 85, 631–632.

## 7. Open Questions

1. **Feynman graphs with external legs** and their integrals (momentum space, Schwinger parameters), and graph isomorphism by partition refinement (morphe `combinatorics`) beyond the brute-force orders.
2. **Hybrid Monte Carlo** for theories whose local updates are slow (critical φ⁴, fermions), with calx's integrators.
3. **Where Markov-chain statistics belong** once the inference foundation exists (its `processes` and `bayes`): they are here because field theory needs them first.
