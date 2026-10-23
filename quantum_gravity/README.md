# quantum_gravity: Spin Foams, Group Field Theory and Random Geometry

*Status: draft (2026-10-23): the Turaev–Viro state sum on triangulated 3-manifolds; BF theory with finite groups on 2-complexes; the Feynman expansion of coloured tensor models and group field theories; matrix models and the Harer–Zagier numbers.* Spin networks evaluated in general, the Barrett–Crane and EPRL spin foams, causal dynamical triangulations, causal-set dynamics and the condensates of group field theory are planned ([decision 0005](../development/decisions/0005-symmetry-spacetime-fields-and-gauge.md)).

## 1. Purpose

The approaches to quantum gravity in which geometry is summed over as combinatorial structure: spin foams, which are lattice gauge theories on 2-complexes; group field theories, whose Feynman graphs generate those 2-complexes; and the matrix and tensor models of random geometry, their simplest instances. The complexes are morphe's (`topology`: simplicial complexes, coloured graphs, CW complexes); the groups and their recoupling theory are physis `symmetry`'s.

## 2. Concepts

- **The Ponzano–Regge model and its regularisation.** A sum over spins on the edges of a triangulation of products of 6j symbols, one per tetrahedron, is a path integral of three-dimensional Euclidean gravity, since the 6j symbol is asymptotically cos(S_Regge + π/4)/√(12πV) (Ponzano and Regge 1968). With the quantum group U_q(sl₂) at q = e^{iπ/r} the spins are bounded by (r − 2)/2 and the sum, normalised by w^{−2} per vertex, is invariant under the Pachner moves: the Turaev–Viro invariant, equal to |τ_r(M)|² for the Witten–Reshetikhin–Turaev invariant τ_r (Turaev and Viro 1992; Turaev 1994; Roberts 1995).
- **The sign convention** used is TV_r(M) = w^{−2V} Σ Π_e (−1)^{2j}[2j+1] Π_triangles (−1)^{j₁+j₂+j₃} Π_tetrahedra {6j}_q with Racah–Wigner-normalised quantum 6j symbols. In the 2–3 move the signs of the new edge and triangles combine to (−1)^{x+S}, the sign of the Biedenharn–Elliott identity; it is the convention of Burton, Maria and Spreer (2018) with the Δ factors moved into the symbols. The checks confirm it on five triangulations of S³, two of S² × S¹ and on lens spaces.
- **Triangulations in the generalised sense** (`Triangulation3`): tetrahedra with their six edges and four triangles, built from a simplicial complex or from a 4-coloured graph, whose tetrahedra are its vertices, triangles its edges and edges its bicoloured cycles. On a coloured graph the Turaev–Viro sum is the amplitude of that Feynman graph in the Boulatov group field theory with the quantum group.
- **BF theory with a finite group** on a 2-complex: the δ-function lattice gauge theory, Z_BF = |G|^{F−E} #flat = |G|^{χ−1} |Hom(π1, G)| per component. On a surface it is Σ_ρ (dim ρ)^χ (Mednykh's formula; Witten's two-dimensional gauge theory at zero area); on the dual 2-skeleton of a triangulated 3-manifold with v vertices it is |G|^{v−1} |Hom(π1 M, G)|, the Ponzano–Regge divergence made finite. Its topological part is the Dijkgraaf–Witten invariant |Hom(π1, G)|/|G|.
- **Coloured tensor models and group field theories.** The Feynman graphs at order (λλ̄)^p are the (D+1)-coloured graphs on 2p vertices; the free energy's coefficient is (1/p!) Σ over connected tuples (σ_0 = id, σ_1, …, σ_D) of the amplitude, which is N^F in the tensor model. With λ scaled by N^{−D(D−1)/4}, a graph carries N^{D − 2ω/(D−1)!}, ω Gurau's degree, so the 1/N expansion is led by the melonic graphs, counted with a root by the Fuss–Catalan numbers (Gurau 2011; Bonzom, Gurau, Riello and Rivasseau 2011). With a finite group the amplitude is BF theory on the graph's 2-complex (the Boulatov model, 1992, for D = 3).
- **Matrix models.** The GUE moments ⟨tr M^{2n}⟩ = Σ_g ε_g(n) N^{n+1−2g} count the orientable gluings of a 2n-gon by genus; the vertices of a gluing π are the cycles of γπ. The numbers ε_g(n) obey the recursion of Harer and Zagier (1986), and ε_0(n) is the Catalan number ('t Hooft 1974; Di Francesco, Ginsparg and Zinn-Justin 1995).

## 3. Content

| Header | Content |
|---|---|
| `turaev_viro.hpp` | `Triangulation3` (vertices, edges, triangles, tetrahedra); `triangulation(SimplicialComplex)`, `triangulation(ColouredGraph)`; `StateSum turaev_viro(const Triangulation3&, int r)` and its overloads for a `SimplicialComplex` and a 4-coloured `ColouredGraph` (value, colourings visited); `lens_space_tv(p, r)`, |τ_r(L(p, 1))|² |
| `bf.hpp` | `holonomy`; `count_flat` (through π1) and `count_flat_brute` (enumeration); `bf_partition_function`; `dijkgraaf_witten` |
| `tensor_models.hpp` | `for_each_coloured_graph(D, p, f)` (σ_0 = id); `amplitude_exponent`; `free_energy(D, p)` by degree; `fuss_catalan(m, p)`; `boulatov_amplitude(graph, finite group)` |
| `matrix_models.hpp` | `for_each_pairing`, `gluing_vertices`, `gluing_genus`, `gluing_complex`; `gluings_by_genus(n)`; `harer_zagier(g, n)`; `gue_moment(n, N)`; `gue_sample`, `trace_power` |

## 4. Checks

- `tests/check_turaev_viro.cpp`: TV_r(S³) = w^{−2} for r = 3..6 on ∂Δ⁴, the suspension of ∂Δ³, the join C₃ * C₃ and the elementary melon, and on a non-melonic coloured graph of S³ for r = 3..7; TV_r(S² × S¹) = 1 on a coloured graph found by enumeration (spherical bubbles, π1 = ℤ) for r = 3..7, and at r = 3 on a staircase triangulation of ∂Δ³ × C₃; the lens spaces L(p, 1), p = 2..5, built as quotients of joins of cycles by a free rotation (two different triangulations of L(3, 1)), against |τ_r(L(p, 1))|² for r = 3..7; TV₃(ℝP³) = 0.
- `tests/check_quantum_gravity.cpp`: flat connections through π1 against enumeration on eight complexes; Z_BF of surfaces of genus 0..3 against Mednykh's formula for S3 and Q8; Boulatov amplitudes of two S³ graphs; the Dijkgraaf–Witten invariant of a lens space; Gurau's scaling on 14 782 connected coloured graphs (D = 2, 3, 4); the melonic free energy 1, 2, 22/3, 35 for D = 3, p = 1..4; ε_g(n) against the Harer–Zagier recursion for n ≤ 7 and ε_0(n) against the Catalan numbers; the genus of each gluing against morphe's classification of surfaces; GUE moments of order 2, 4, 6 by Monte Carlo.

## 5. Dependencies

physis `symmetry` (quantum 6j symbols and dimensions); morphe `topology` (simplicial complexes, coloured graphs, CW complexes, π1, surfaces), `algebra` (finite groups and homomorphisms) and `arithmetic`; calx `random` (GUE samples), and in the checks `verification`.

## 6. References

- Ponzano, G. and Regge, T. (1968), Semiclassical limit of Racah coefficients, in *Spectroscopic and Group Theoretical Methods in Physics*, North-Holland, 1–58.
- Turaev, V. G. and Viro, O. Ya. (1992), State sum invariants of 3-manifolds and quantum 6j-symbols, *Topology* 31, 865–902.
- Kirillov, A. N. and Reshetikhin, N. Yu. (1989), Representations of the algebra U_q(sl(2)), q-orthogonal polynomials and invariants of links, in *Infinite-Dimensional Lie Algebras and Groups*, World Scientific, 285–339.
- Turaev, V. G. (1994), *Quantum Invariants of Knots and 3-Manifolds*, de Gruyter.
- Roberts, J. (1995), Skein theory and Turaev–Viro invariants, *Topology* 34, 771–787.
- Jeffrey, L. C. (1992), Chern–Simons–Witten invariants of lens spaces and torus bundles, and the semiclassical approximation, *Commun. Math. Phys.* 147, 563–604.
- Burton, B. A., Maria, C. and Spreer, J. (2018), Algorithms and complexity for Turaev–Viro invariants, *Discrete Comput. Geom.* 59, 365–399.
- Boulatov, D. V. (1992), A model of three-dimensional lattice gravity, *Mod. Phys. Lett. A* 7, 1629–1646.
- Ooguri, H. (1992), Topological lattice models in four dimensions, *Mod. Phys. Lett. A* 7, 2799–2810.
- Dijkgraaf, R. and Witten, E. (1990), Topological gauge theories and group cohomology, *Commun. Math. Phys.* 129, 393–429.
- Mednykh, A. D. (1978), Determination of the number of nonequivalent coverings over a compact Riemann surface, *Soviet Math. Dokl.* 19, 318–320.
- Freidel, L. and Louapre, D. (2004), Ponzano–Regge model revisited I: gauge fixing, observables and interacting spinning particles, *Class. Quantum Grav.* 21, 5685–5726.
- Oriti, D. (2009), The group field theory approach to quantum gravity, in *Approaches to Quantum Gravity*, ed. D. Oriti, Cambridge University Press, 310–331.
- Gurau, R. (2011), The 1/N expansion of colored tensor models, *Ann. Henri Poincaré* 12, 829–847.
- Bonzom, V., Gurau, R., Riello, A. and Rivasseau, V. (2011), Critical behavior of colored tensor models in the large N limit, *Nucl. Phys. B* 853, 174–195.
- 't Hooft, G. (1974), A planar diagram theory for strong interactions, *Nucl. Phys. B* 72, 461–473.
- Harer, J. and Zagier, D. (1986), The Euler characteristic of the moduli space of curves, *Invent. Math.* 85, 457–485.
- Di Francesco, P., Ginsparg, P. and Zinn-Justin, J. (1995), 2D gravity and random matrices, *Phys. Rep.* 254, 1–133.

## 7. Open Questions

1. **The EPRL and Barrett–Crane spin foams** of four-dimensional gravity need the unitary representations of SL(2, ℂ) and its Spin(4) analogues (physis `symmetry`, planned), and their vertex amplitudes are integrals that will need calx's quadrature on non-compact groups.
2. **Faster state sums.** Enumeration is exponential in the number of edges; tree decompositions of the dual graph (Burton, Maria and Spreer 2018) make the Turaev–Viro sum fixed-parameter tractable.
3. **Causal dynamical triangulations** (Ambjørn, Jurkiewicz and Loll): Monte Carlo over foliated triangulations with Pachner moves, on physis `field_theory`'s Markov chains and morphe's complexes.
4. **Causal-set dynamics**: sequential growth models and the Benincasa–Dowker action, on physis `spacetime`'s causal sets.
5. **Group field theory condensates** and their cosmology: the mean-field (Gross–Pitaevskii) equations of a GFT on SU(2)^4.
6. **Non-orientable pseudomanifolds** from non-bipartite coloured graphs, and the Boulatov model with compact Lie groups beyond finite ones, with its divergences regularised.
