# 0011. Point Groups, Stranded Field Theories on Finite Groups, and Reaction Kinetics

- **Status:** Accepted, provisionally (point groups are found within a relative tolerance, for finite sets off a line; the stranded theories are real fields with the identity propagator and no gauge projection in their Wick expansion; the kinetics' exact absorption is dense beyond its acyclic case)
- **Date:** 2026-10-29
- **Decided by:** the user asked for a chemistry workbench whose first plug-in executes their manuscript *An Inquiry into Chemical Structure Generation through Symmetry Breaking* (v1.4) and verifies its formalism, under the rule that general methods are written in the foundations. Claude wrote the methods here and this record; the workbench is the infrastructure's [ADR 0051](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0051-the-chemistry-workbench.md).

## 1. Context

The manuscript has two formalisms. The first, a group-equivariant grammar of chemical graphs, needs:
- the symmetry groups of molecules, with the responses they allow;
- reaction networks, as mass action and as Markov chains.

The second, a group-field-inspired generator of 2-complexes, needs:
- Boulatov's tetrahedral interaction on a *real, uncoloured* field over a finite group;
- its Feynman graphs with their 2-complexes.

physis had group field theories in their coloured and tensorial forms only (decision 0006), no point groups, and a `chemistry` module only planned.

## 2. Decision

### 2.1 `symmetry`: Point Groups

| Header | Content |
|---|---|
| `point_group.hpp` | `LabelledPoints`; `congruences(A, B)`, the isometries carrying one labelled set onto another (found from the images of an independent triple, or of a pair and the normal for a planar set); `point_group(A)` → `PointGroup`: its elements (identity first), multiplication table as a morphe `FiniteGroup`, Schoenflies symbol (C∞v and D∞h on a line), symmetry elements for drawing, `chiral()`, `fixed_vectors(axial)` (the mean of tr Q or det Q · tr Q: Neumann's principle); `generate`, `rotation`, `reflection` |

*Rationale.* The group is built from matrices, not from the points' permutation, which a planar set does not render faithfully: the reflection in its plane moves no point. Handing the table to morphe gives conjugacy classes and exact character tables (morphe decision 0006) without a second implementation.

### 2.2 `quantum_gravity`: Stranded Field Theories on Finite Groups

| Header | Content |
|---|---|
| `stranded.hpp` | `StrandedVertex::from_words` (one word of index names per factor, each used twice); `tetrahedral_vertex()` (Boulatov's), `melonic_vertex()`; `interaction(v, G, ψ)`, exact over any number type; `slot_classes`, `translate`, `gauge_project`; `for_each_pairing`, `stranded_graph` (faces, 2-complex as morphe's `CWComplex`, connectedness); `vertex_automorphisms`, `wreath_maps`, `feynman_classes` (classes under Aut(v) ≀ Sₙ, with their pairings and automorphisms); `gaussian_moment`, ⟨Vⁿ⟩ by brute force |

*Rationale.* With independent standard normal values, Wick's theorem gives ⟨Vⁿ⟩ = Σ over pairings of |G|^F, F the closed strands. That identity, checked against the brute-force moment, is what makes "the field generates labelled 2-complexes" exact. Classes are named by their least image under the wreath product, whose order |Aut v|ⁿ n! is small. The canonical forms of labelled graphs would not do: disjoint copies of one vertex leave refinement cells too large to order.

The amplitudes with the gauge projection come from `bf.hpp`, as for the coloured models.

### 2.3 `chemistry`: Reaction Kinetics (the Module's First Content)

| Header | Content |
|---|---|
| `kinetics.hpp` | `Reaction`, `Network` (species, reactions with rate constants as rationals); `stoichiometric_matrix`, `conservation_laws` (a basis over ℤ of the left kernel, by morphe's Smith form), `balance` (AN and zᵀN); `mass_action` (for calx's integrators); the Markov chain: `stochastic_constant` (Gillespie's c = k α!/Ω^{\|α\|−1}), `propensity`, `drift`, `reachable`, `generator`, `absorption` (exact in rationals or in floating point: along a topological order when the states are acyclic, from the occupation measure otherwise), `gillespie` (calx's streams); symmetry: `permute`, `image`, `symmetry_of` |

The module names no species: its checks use abstract networks. Chemical thermodynamics and the coupling of reactions to transport stay planned.

## 3. Checks

- **`symmetry/tests/check_point_groups.cpp` (61 checks).**
  - For each of 26 groups, C1 to Ih, the orbit of four generic labelled points recovers the group's order, Schoenflies symbol and number of conjugacy classes (Cotton 1990).
  - It also recovers the polar groups, the groups fixing an axial vector, and the chiral groups.
  - A regular hexagon is D6h. Points on a line give D∞h or C∞v.
  - Mirror-image tetrahedra with four labels are congruent only by one improper isometry; a rotated copy by one proper isometry.
- **`quantum_gravity/tests/check_stranded.cpp` (22 checks).**
  - V_tet[1] = 64 and V_tet[(−1)^a] = 0 over Z₂, and V_tet along the eight characters (64 exactly when s₁ = s₃).
  - Invariance under (g, h, g) for Z₂, Z₃ and S₃, and its failure elsewhere.
  - The gauge projection annihilating the odd functional over Z₂.
  - Wick sums equal to brute-force moments for n ≤ 2 over 1, Z₂, Z₃ (n = 1 over S₃), n = 3 over 1 and Z₂, and the melonic vertex.
  - The vertex's 4 automorphisms; class sizes |Aut v|ⁿ n!/|Aut Γ| for n ≤ 3.
  - The 24 two-vertex gluings: 4 with six faces (Betti numbers (1, 0, 3), BF amplitude 8 over Z₂ and 216 over S₃), 4 with four, 16 with three.
- **`chemistry/tests/check_kinetics.cpp` (14 checks).**
  - Conservation laws in the left kernel; the law ±(1, 1, 1, 2) of an autocatalytic network; balance caught when broken.
  - Mass action against closed forms for X ⇌ Y and 2X → Y; the stochastic constant.
  - The neutral network X + Y → 2X, X + Y → 2Y fixing X with probability i/N exactly; its generator commuting with the exchange.
  - Dynkin's drift at every reachable state.
  - Gillespie's paths against the Poisson stationary law of a birth–death network.

## 4. Consequences

- The Chemistry Workbench's first plug-in computes its point groups, its field theory and its kinetics here (ADR 0051).
- A chemistry domain, when it is founded, will find its general methods already in place.
- `symmetry` now links morphe's `algebra`.
