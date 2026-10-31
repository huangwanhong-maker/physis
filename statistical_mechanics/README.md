# statistical_mechanics

*Status: draft (2026-10-31), [decision 0016](../development/decisions/0016-kinematic-waves-master-equations-and-reaction-diffusion.md).* Markov jump processes on finitely many states: their master equation and stationary law, their irreversibility as currents round cycles, first passage and absorption, and the entropy of mixing. Ensembles and partition functions, lattice and spin models, Glauber and Metropolis dynamics, and mean-field theory are planned (README §2).

## 1. Purpose

The statistical laws of systems that jump among states, stated so that a domain names only its states and its rates: which places a walker moves among, and how often, belong to a domain; how any such process relaxes, how far it is from equilibrium, how long it takes to arrive and where it ends, belong here. The entropy of mixing measures how far a population of kinds spread over places is from being mixed.

## 2. Concepts

- **The master equation.** Rates $k_{ij}\ge 0$ of jumping from $i$ to $j$ give the generator $Q$ ($Q_{ij}=k_{ij}$, $Q_{ii}=-\sum_{j\ne i}k_{ij}$), and a distribution evolves as $p(t)=p(0)e^{Qt}$ (van Kampen 2007, ch. V). An irreducible chain has one stationary law $\pi Q=0$, found by the GTH algorithm, which eliminates states one by one and forms no difference of positive numbers, so that every $\pi_i$ keeps a small relative error however small it is (Grassmann, Taksar and Heyman 1985).
- **Irreversibility.** The stationary current $J_{ij}=\pi_ik_{ij}-\pi_jk_{ji}$ is divergence-free, so it is a sum of currents round the fundamental cycles of a spanning tree (morphe's spanning forest): one for each chord. The entropy production $\sigma=\tfrac12\sum J_{ij}\ln(\pi_ik_{ij}/\pi_jk_{ji})\ge 0$ equals $\sum_C J(C)A(C)$, $A(C)=\ln\prod_C k_{ab}/k_{ba}$ the cycle's affinity, which no $\pi$ enters (Schnakenberg 1976, eqs. 7.6 and 7.11). Detailed balance holds exactly when every affinity vanishes: Kolmogorov's criterion. Movement whose rates follow the gradient of a potential cannot circulate.
- **The H-theorem.** $D(p(t)\,\|\,\pi)$ never increases.
- **First passage.** The time to leave a set of transient states has a phase-type law, density $\alpha e^{T\tau}t$, survival $\alpha e^{T\tau}\mathbf 1$, moments $k!\,\alpha(-T)^{-k}\mathbf 1$ (Neuts 1981); where it leaves to is row $a$ of $(-T)^{-1}R$; the mean time to a target solves $\sum_jk_{ij}(m_j-m_i)=-1$.
- **Mixing.** For counts of kinds at places, $M=H(Q)-H(Q|X)=\sum_xp_xD(p(\cdot|x)\,\|\,p(\cdot))$, the mutual information between place and kind: zero exactly when every place holds the whole's proportions (Theil and Finizza 1971). Over groups of places, $M=M_\text{between}+\sum_gp_gM_g$ exactly (the chain rule), and an exchange of compositions among places never raises $M$ (data processing).

## 3. Content

| Header | Content |
|---|---|
| `master_equation.hpp` | `generator`, `irreducible`, `stationary` (GTH), `evolve` (calx's `expm`), `relative_entropy`; `Cycle`, `Irreversibility`, `irreversibility` (currents, entropy production, a fundamental cycle basis with affinities and currents, Kolmogorov's criterion) |
| `first_passage.hpp` | `PhaseType` (`density`, `survival`, `moment`, `quantile`, `atom`), `phase_type`, `absorption_probabilities`, `mean_first_passage` |
| `mixing.hpp` | `Mixing`, `mixing` (H(Q), H(Q\|X), M, its normalised value and local terms); `MixingDecomposition`, `mixing_decomposition` (between and within groups) |

## 4. Checks

- `tests/check_master_equation.cpp` (25 checks): birth–death stationary laws against detailed balance entry by entry over forty-six orders of magnitude; π Q = 0 on random chains; a reducible chain refused; driven rings of 3 and 7 states against $A=n\ln(k_+/k_-)$, $J=(k_+-k_-)/n$; Schnakenberg's identity on 60 random graphs, every edge's term non-negative, one cycle per chord; rates from energies giving the Boltzmann law with σ = 0, and one rate changed breaking Kolmogorov's criterion; the H-theorem; `evolve` against explicit Euler (order 1) and against the spectral decomposition of a reversible generator; the relative entropy's closed form and support.
- `tests/check_first_passage.cpp` (15 checks): the exponential and Erlang laws (density, survival, moments, quantiles); an atom at zero; 40000 sampled paths against the mean, the median and the law of the absorbing state; gambler's ruin; mean first passage on a symmetric ring, $m_i=i(n-i)/2k$, Kac's identity $\sum_jk_{ij}m_j=1/\pi_i-1$ for every target, and the phase-type mean; a trapped transient state refused.
- `tests/check_mixing.cpp` (8 checks): a hand-computed table; $M=0$ for proportional counts and $M=H(Q)$ under separation; the identities on random tables; the chain rule over random partitions; invariance under permuting places and kinds; the data-processing inequality under random doubly stochastic exchanges.

## 5. Dependencies

calx `linalg` (dense matrices, `solve`, `inverse`, `expm`, `Rectangular`; `symmetric_eigen` in the checks); morphe `combinatorics` (the spanning forest of the graph of rates). No other module of physis: the planned ensembles will add `thermodynamics`, `dynamics` and `quantum` (architecture §2).

## 6. References

- Cover, T. M. and Thomas, J. A. (2006). *Elements of Information Theory*, 2nd ed. Wiley.
- Grassmann, W. K., Taksar, M. I. and Heyman, D. P. (1985). Regenerative analysis and steady state distributions for Markov chains. *Oper. Res.* 33, 1107–1116.
- Latouche, G. and Ramaswami, V. (1999). *Introduction to Matrix Analytic Methods in Stochastic Modeling*. SIAM.
- Mora, R. and Ruiz-Castillo, J. (2011). Entropy-based segregation indices. *Sociol. Methodol.* 41, 159–194.
- Neuts, M. F. (1981). *Matrix-Geometric Solutions in Stochastic Models*. Johns Hopkins University Press.
- Schnakenberg, J. (1976). Network theory of microscopic and macroscopic behavior of master equation systems. *Rev. Mod. Phys.* 48, 571–585.
- Theil, H. and Finizza, A. J. (1971). A note on the measurement of racial integration of schools by means of informational concepts. *J. Math. Sociol.* 1, 187–193.
- van Kampen, N. G. (2007). *Stochastic Processes in Physics and Chemistry*, 3rd ed. Elsevier.

## 7. Open Questions

1. **Large chains.** Every routine is dense, O(n³): a few thousand states. Sparse generators would need Krylov exponentials (expmv) and iterative stationary solvers; uniformisation is an alternative for evolve.
2. **Fluctuations of currents.** The large deviations of time-integrated currents (the tilted generator's leading eigenvalue) and the fluctuation theorem's symmetry are the natural next step beyond the mean entropy production.
3. **The ensembles** of the module's planned scope: partition functions, spin models and their Monte Carlo, beside these processes.
