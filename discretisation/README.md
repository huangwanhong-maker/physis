# discretisation: Structure-Preserving Discretisation

*Status: draft (2026-10-24): discrete exterior calculus on morphe's CW complexes, with real cochains, their coboundaries, potentials and the Hodge decomposition ([decision 0008](../development/decisions/0008-cochains-hodge-and-the-sphere.md)).* Finite differences, finite volumes, finite elements, spectral methods in C++, and geometric and Lie-group integrators are planned. A JavaScript stand-in for the spectral transform on the sphere and for longitude–latitude grids is in [`js/`](../js/) ([decision 0003](../development/decisions/0003-provisional-javascript-stand-in.md)).

## 1. Purpose

Discretisations that keep the structure of the continuous problem: d∘d = 0, Stokes's theorem, the orthogonal splitting of forms. The first is discrete exterior calculus on morphe's complexes, which the economics domain uses for the circulation of value and electromagnetism on complexes will use for its fields.

## 2. Concepts

- **Cochains.** A k-cochain is a real number on each k-cell; the coboundary is the transpose of the cellular boundary. (d₀φ)(e) = φ(head) − φ(tail), and (d₁ω)(f) is the integral of ω round the attaching word of f, so that d₁d₀ = 0 and Stokes's theorem hold by construction (Desbrun, Hirani, Leok and Marsden 2005).
- **Potentials and the tree gauge.** Integrating a 1-cochain along a spanning forest gives a potential whose coboundary agrees with it on the forest. The values left on the chords are the integrals round the fundamental cycles: complete invariants modulo exact cochains, all zero exactly when the cochain is exact.
- **The Hodge decomposition.** With diagonal Hodge stars (weights on edges and 2-cells), a 1-cochain splits orthogonally into d₀α (exact: a potential's gradient), δ₁γ (co-exact: curl round the 2-cells) and h (harmonic: closed and co-closed, carrying the integrals round the cycles no 2-cell fills, of dimension b₁) (Eckmann 1944; Lim 2020).
  - The dense method solves the two normal equations by calx's solve, made definite by bases of their kernels: the components, and the 2-cycles of the subcomplex of the 2-cells and their edges (morphe).
  - The sparse method, beyond 1500 vertices or 2-cells, uses calx's conjugate gradients on the singular, consistent systems.
- **Flows with conductances.** A flow is a 1-chain. With a conductance c_e > 0 on each edge (Ohm's law), its drop F / c is a cochain, decomposed with W₁ = c; multiplied back by c, the three parts are flows, orthogonal in the energy Σ F² / c. The driven flow c d₀α carries the whole boundary and the other two none (Kirchhoff's current law, Kirchhoff 1847); the local flow is the same amount round every edge of a 2-cell. The pairing Σ F ω of a flow with a cochain splits into the pairings of like parts.

## 3. Content

| Header | Content |
|---|---|
| `cochains.hpp` | `Cochain`; `d0`, `d1`, `d0_transpose`, `d1_transpose`; `integral(ω, path)`; `potential(X, ω)` → `Potential{phi, chords, defect, in_tree}`; `components(X)` |
| `hodge.hpp` | `Weights{edge, face}`, `HodgeOptions{method, dense_limit, tolerance}`; `hodge_decomposition(X, ω, W, options)` and `hodge_decompositions` (many at once) → `HodgeDecomposition{exact, coexact, harmonic, potential, curl, norms, sparse, iterations, residual}`; `harmonic_basis(X, W)`; `inner1`; `flow_decomposition(X, F, conductance, face_weights, options)` → `FlowDecomposition{driven, local, global, potential, circulation, energies, drop}` |

## 4. Checks

`tests/check_cochains.cpp` (107 checks):

- **On eight complexes** (the 4 × 3 lattice torus, the surface of genus 2, the projective plane, the Klein bottle, the dunce hat, the octahedral sphere, a grid disc and a graph):
  - d₁d₀ = 0 and Stokes;
  - an exact cochain's potential recovered;
  - with random weights, the parts adding up to ω, orthogonal, the harmonic part closed and co-closed, the co-exact part carrying the whole curvature, and the harmonic part decomposing into itself;
  - the dense and sparse methods agreeing to 10⁻⁸;
  - the harmonic dimension equal to b₁(X; ℝ);
  - a flow with random conductances: its driven part carrying its whole boundary and the circulating parts none, the energies adding up, its pairing with a cochain splitting into like parts.
- **A 60 × 60 lattice torus** decomposed by the sparse method.
- **Closed forms:** a circle's chord carrying its integral; the weighted cycle's harmonic part c/w_e; a triangle's integral ln 0.6, harmonic while unfilled and curvature once filled; Ohm's law on two parallel conductances.

## 5. Dependencies

morphe `topology` (CW complexes, homology with generators), `combinatorics` (spanning forests), `algebra`; calx `linalg` (dense solve, CSR, conjugate gradients). In the checks, calx `random` and `verification`.

## 6. References

- Desbrun, M., Hirani, A. N., Leok, M. and Marsden, J. E. (2005), Discrete exterior calculus, arXiv:math/0508341.
- Eckmann, B. (1944), Harmonische Funktionen und Randwertaufgaben in einem Komplex, *Comment. Math. Helv.* 17, 240–255; Dodziuk, J. (1976), Finite-difference approach to the Hodge theory of harmonic forms, *Amer. J. Math.* 98, 79–104.
- Lim, L.-H. (2020), Hodge Laplacians on graphs, *SIAM Review* 62, 685–715.
- Jiang, X., Lim, L.-H., Yao, Y. and Ye, Y. (2011), Statistical ranking and combinatorial Hodge theory, *Math. Program.* 127, 203–244.
- Kirchhoff, G. (1847), Ueber die Auflösung der Gleichungen, auf welche man bei der Untersuchung der linearen Vertheilung galvanischer Ströme geführt wird, *Ann. Phys. Chem.* 72, 497–508.

## 7. Open Questions

1. **Non-diagonal Hodge stars** (circumcentric duals of simplicial meshes), for fields on curved surfaces.
2. **Higher degrees**: 2-cochains on 3-complexes, for electromagnetism in three dimensions.
3. **A sparse factorisation** beside conjugate gradients, for many right-hand sides on one large complex.
