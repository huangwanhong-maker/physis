# 0008. Cochains and the Hodge Decomposition; Points on the Sphere

- **Status:** Accepted, provisionally (diagonal Hodge stars only; the sparse method's tolerance is fixed by hand)
- **Date:** 2026-10-24
- **Decided by:** the user asked for an economics domain that renders and investigates the circulation of value on a map and in spacetime, and set the rule that a domain writes no general method of its own. Claude wrote the general methods here and this record; the domain's use of them is the infrastructure's [ADR 0044](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0044-value-in-space-and-time.md).

## 1. Context

The economics domain carries value on networks: a transport on each leg, whose logarithm is a 1-cochain, and a flow on each leg, a 1-chain. Its questions are those of discrete exterior calculus:

- the integral of a cochain round a loop;
- whether a cochain is exact, and its potential;
- the orthogonal split of a cochain into exact, co-exact and harmonic parts, on a complex whose faces are given.

Its faces come from geography: the triangles of a set of places on the sphere. Discrete exterior calculus on morphe's complexes was planned in `discretisation`, where only a JavaScript stand-in for spectral transforms was written (decision 0003). Computational geometry was planned in `computational_geometry`, with nothing written.

## 2. Decision

### 2.1 `discretisation`: Cochains

`cochains.hpp` provides real k-cochains on a morphe `CWComplex`:

- the coboundaries d₀ and d₁ (d₁ the integral round each 2-cell's word) and their transposes, so that d₁d₀ = 0 and Stokes's theorem hold by construction;
- the integral along an edge path;
- the potential of a 1-cochain along a spanning forest of the 1-skeleton (morphe's), with its defect and its chord values: the tree gauge, and the potential exactly when the defect vanishes.

### 2.2 `discretisation`: the Hodge Decomposition

`hodge.hpp` splits a 1-cochain orthogonally, for diagonal Hodge stars (weights on edges and 2-cells), into d₀α + δ₁γ + h, with δ₁ = W₁⁻¹d₁ᵀW₂. The two normal equations, d₀ᵀW₁d₀ α = d₀ᵀW₁ω and d₁W₁⁻¹d₁ᵀ γ′ = d₁ω, are singular; their kernels are the locally constant functions and the 2-cycles.

- **Dense method** (up to 1500 vertices and 2-cells): calx's solve, made definite by the kernels' bases. The bases are morphe's: the components of a spanning forest, and the free generators of H₂ of the subcomplex of the 2-cells and their edges, which is all a 2-cycle can use.
- **Sparse method** (beyond): calx's conjugate gradients on the singular, consistent systems, with the potential projected off the components. The exact and co-exact parts are unique whatever solution is reached.

`harmonic_basis` gives an orthonormal basis of the harmonic space, of dimension b₁(X; ℝ). Many cochains are decomposed with one assembly.

**Flows with conductances.** A flow is a 1-chain, not a cochain, and its decomposition needs a metric. `flow_decomposition` takes a conductance c_e > 0 on each edge (Ohm's law), decomposes the drop F / c with W₁ = c and multiplies back by c. The three flows are orthogonal in the energy Σ F² / c. The driven one, c d₀α, carries the whole boundary of F: the other two are divergence-free, which is Kirchhoff's current law. The local one is the same amount round every edge of each 2-cell. The pairing of a flow with any 1-cochain, which needs no metric, splits into the pairings of like parts. With unit conductances, every edge would carry the same driven flow for the same difference of potential, whatever its size; on networks whose links differ by orders of magnitude, such as the world's trade, the driven and circulating parts would then each exceed the flow they make up.

### 2.3 `computational_geometry`: Points on the Sphere

`sphere.hpp` provides:

- unit vectors from latitude and longitude;
- the great-circle angle, accurate at every separation;
- the signed area of a spherical triangle (Van Oosterom and Strackee);
- the spherical Delaunay triangulation, as the convex hull of the points in ℝ³ (Brown 1979), built incrementally with a deterministic perturbation of 10⁻⁹. Degenerate configurations, such as the cube's faces with four points on each circle, are triangulated in one valid way and lose no point. Coincident points are merged, and whether the hull surrounds the centre is reported.

## 3. Alternatives Considered

- **The domain's own least squares and triangulation.** Rejected by the user's rule: electromagnetism and fluids on complexes will need the same cochains, and every map the same triangulation.
- **A borrowed geometry library (CGAL, qhull).** CGAL's spherical triangulations are GPL; qhull would sit behind our interface well but is a large dependency for a hull of a few thousand points. The incremental hull is short, exact in its tests and checked against closed forms.
- **The sparse method everywhere.** The dense method is exact to rounding and needs no tolerance, so it is kept where it is cheap; the two are checked against each other.

## 4. Consequences

- Two modules gain C++: `discretisation` (cochains, the Hodge decomposition, flows with conductances) and `computational_geometry` (points on the sphere). Their checks:
  - **cochains and Hodge:** d₁d₀ = 0 and Stokes on eight complexes; potentials recovered; the three parts adding up, orthogonal, closed and co-closed, and idempotent; the dense and sparse methods agreeing; harmonic dimensions equal to b₁(X; ℝ) on the torus, genus 2, the projective plane, the Klein bottle, the dunce hat, the octahedral sphere, a grid disc and a graph; a 60 × 60 lattice torus decomposed by the sparse method; the weighted cycle's closed form; the triangle's integral, harmonic while unfilled and curvature once filled;
  - **flows:** on the eight complexes, with random conductances, Kirchhoff's current law, the energies adding up and the pairing splitting into like parts; Ohm's law on two parallel conductances;
  - **points on the sphere:** the octahedron and the cube; 500 random points with F = 2V − 4, E = 3V − 6, every triangle anticlockwise and Delaunay, areas adding up to 4π; a hemisphere that does not cover the sphere; duplicates merged.
- `discretisation` now holds C++ beside its JavaScript stand-in, and its README says which is which.
