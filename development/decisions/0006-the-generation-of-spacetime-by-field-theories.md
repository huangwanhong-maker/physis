# 0006. The Generation of Spacetime by Field Theories

- **Status:** Accepted, provisionally (the framework is a first design; the large-N and condensate methods are the simplest of their kind)
- **Date:** 2026-10-23
- **Decided by:** the user set out a formalism of emergent spacetime, an action S[φ] whose perturbative expansion generates 2-complexes (spin foams, discrete spacetimes) and whose condensed phase gives a continuous spacetime, asked whether the libraries can investigate it and show it (spin foams above all), and asked that it hold for group field theory, quantum field theory and other theories alike. Claude designed the framework and wrote this record; the composition in the infrastructure is its [ADR 0040](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0040-the-generation-of-spacetime.md).

## 1. Context

Decision 0005 gave physis the pieces: fields and Feynman graphs (`field_theory`), gauge fields on 2-complexes (`gauge`), state sums and coloured tensor models (`quantum_gravity`), and spacetimes discrete and continuous (`spacetime`). It did not give the passage from one to the next: an action as an object, from which Feynman rules follow; the spin foam as an object, with spins and amplitudes; the phase in which the quanta condense; and the continuum that results, with numbers to compare.

The passage is the same in a family of theories, ordered by the rank of the field (Di Francesco, Ginsparg and Zinn-Justin 1995; Gurau 2016, *Random Tensors*; Oriti 2009, 2017):

| Rank | Theory | Feynman graphs | Discrete spacetime | Large N or condensate | Continuum |
|---|---|---|---|---|---|
| 0 | complex scalar, \|φ\|^{2k} | graphs | the graph | mean field | — |
| 1 | vector, U(N) | 2-coloured graphs | cacti | the gap equation | — |
| 2 | matrix | ribbon graphs | surfaces (quadrangulations) | eigenvalue density | 2D quantum gravity (γ = −1/2) |
| d ≥ 3 | tensor | (d+1)-coloured graphs | pseudomanifolds | melonic resummation | branched polymers (γ = 1/2) |
| d | group field theory | (d+1)-coloured graphs | spin foams | condensates | Friedmann dynamics with a bounce |

## 2. Decision

### 2.1 Five Stages for Every Theory

1. **The action.** A complex field of rank r with indices in a set of size N (U(N)^r-invariant theories) or in a group G (group field theories, with the gauge projection); interactions given by **bubbles**, bipartite r-coloured graphs whose vertices are the fields φ and φ̄ and whose colour-i edges contract the i-th indices; couplings λ_B.
2. **The perturbative expansion.** With the Gaussian measure ⟨φ_a φ̄_b⟩ = δ_ab, a Wick contraction of n_B copies of each bubble is a matching σ₀ of black to white vertices, the propagators, colour 0. The Feynman graph is the (r+1)-coloured graph so formed; its amplitude is N^{F₀} for U(N) theories, F₀ the number of faces of colours (0, i), and a BF or spin-foam state sum on its 2-complex for group field theories. Series in λ and N follow exactly, as integers and rationals.
3. **The discrete spacetime.** The Feynman graph's 2-complex (its vertices, edges and bicoloured cycles) and its pseudomanifold (morphe `topology`); for matrices a surface of genus g, weighted N^{2−2g}; for tensors a pseudomanifold of degree ω; for group field theories a **spin foam**, the 2-complex with spins on its faces and amplitudes on its vertices, whose vertices are quantum simplices.
4. **The large-N or condensed phase.** The saddle point or Schwinger–Dyson equation of the dominant graphs: the gap equation of vector models, the eigenvalue density of matrix models (Brézin, Itzykson, Parisi and Zuber 1978), the melonic two-point function of tensor models (Bonzom, Gurau, Riello and Rivasseau 2011); for group field theories the condensate, a mean field σ(g_I) whose relational dynamics is cosmology (Gielen, Oriti and Sindoni 2013; Oriti, Sindoni and Wilson-Ewing 2016).
5. **The continuum.** Critical couplings and exponents (the string susceptibility), the effective equations of the condensate (Friedmann's with a bounce at a critical density), and the emergent dimensions of the discrete spacetimes: spectral (by diffusion) and Hausdorff (by the growth of balls).

### 2.2 Placement

| Module | New content |
|---|---|
| morphe `topology` | bubbles (bipartite D-coloured graphs: quartic melonic, pillow, necklace, the cycles of matrix models, the pairs of vector models) and Feynman graphs glued from bubbles by a propagator matching |
| physis `field_theory` | `invariant.hpp`: theories of rank r, their Feynman expansion and its exact series; `large_n.hpp`: gap equation, eigenvalue density and planar free energy, melonic resummation, critical points and exponents |
| physis `quantum_gravity` | `gft.hpp`: group field theories (Boulatov, Ooguri, quartic melonic) as invariant theories over a group, their Feynman amplitudes; `spin_foam.hpp`: spin foams on 2-complexes, vertex amplitudes from boundary spin networks, the geometry of quantum tetrahedra and its large-spin (Regge) limit, sampling of spins; `random_geometry.hpp`: random planar quadrangulations (2D quantum gravity) and random melonic graphs (branched polymers); `condensate.hpp`: condensate cosmology and the effective Friedmann equation |
| physis `spacetime` | `dimensions.hpp`: spectral and Hausdorff dimensions of discrete spaces |

The real scalar φ^k theory, whose Wick contractions are not bipartite, keeps its graphs in `field_theory/feynman.hpp`.

### 2.3 Conventions

Fields are complex; the Gaussian measure is normalised to ⟨φ_a φ̄_b⟩ = δ_ab; an interaction is the trace invariant I_B(φ, φ̄) of its bubble, entering as e^{−λ_B I_B}; the large-N scaling λ_B → λ_B N^{−s_B} is stated with each model. Condensate cosmology is in Planck units (ħ = c = 1, G a parameter), its time the relational clock of a massless scalar field. They are added to the glossary.

### 2.4 Planned

Four-dimensional spin foams beyond BF (Barrett–Crane, EPRL, on SL(2, ℂ)); causal dynamical triangulations; causal-set dynamics; the functional renormalisation group of tensor models and group field theories; double-scaling limits as flows to the continuum.

## 3. Alternatives Considered

- **A pipeline for group field theory alone.** The user asked for quantum field theory as well; the stages are the same, and matrix and tensor models are where the continuum limit is known exactly, so they check the machinery that group field theory then uses.
- **A new module for emergence.** The stages already have homes (fields, gauge, quantum gravity, spacetime); a module of their composition would own nothing. The composition is shown in an application (ADR 0040).
- **Real fields for every rank.** Wick contractions of real fields are not bipartite and their graphs are not coloured; complex fields give the coloured graphs whose combinatorics morphe holds, as is standard for tensor models (Gurau 2016).

## 4. Consequences

- `field_theory`, `quantum_gravity` and `spacetime` grow by the files of §2.2, each with its checks against exact results (moments of Gaussian vectors and matrices, Catalan and Fuss–Catalan numbers, the BIPZ series and its critical exponent, the identities of the condensate's Friedmann equation, return probabilities of lattice random walks).
- The Complex Atlas gains a mode that walks the five stages for a chosen theory, and a view of spin foams (ADR 0040).
