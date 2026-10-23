# symmetry: Lie Groups and Their Representations

*Status: draft (2026-10-23): U(1), SU(2), SU(N), the Lorentz group and SL(2, ℂ); the exact recoupling theory of SU(2); quantum 6j symbols at roots of unity; Weyl integration and character expansions.* Sp(N), the Poincaré group, the representations of SL(2, ℂ), Wigner D-functions, symmetry breaking and Noether currents are planned ([decision 0005](../development/decisions/0005-symmetry-spacetime-fields-and-gauge.md)).

## 1. Purpose

The groups that physics uses as symmetries, with their algebras, invariant measures and representations: the gauge groups of lattice gauge theory (physis `gauge`), the recoupling theory of spin networks and spin foams (physis `quantum_gravity`), and the Lorentz group of spacetime (physis `spacetime`). Finite and finitely presented groups are morphe's (`algebra`); a group here is a Lie group, computed in floating point except where its representation theory is exact.

## 2. Concepts

- **SU(2) is the 3-sphere of unit quaternions.** U = a0 I + i a·σ; the product is the quaternion product, the inverse the conjugate, and every element is conjugate to diag(e^{iθ/2}, e^{−iθ/2}) with a0 = cos(θ/2), θ ∈ [0, 2π] its class angle. The generators T_k = −(i/2)σ_k exponentiate to rotations, and R(U)_ij = ½ tr(σ_i U σ_j U†) is the double cover SU(2) → SO(3), active and right-handed. The irreducible representations are the spins j, with characters χ_j(θ) = sin((2j+1)θ/2)/sin(θ/2), evaluated as Chebyshev polynomials U_{2j}(cos(θ/2)). The Haar measure is the uniform measure on S³: a normalised Gaussian 4-vector.
- **SU(N)** as matrices of run-time order, with the generalised Gell-Mann generators t_a (tr t_a t_b = δ_ab/2), the exponential exp(−i x^a t_a) through the real form, reunitarisation by Gram–Schmidt and the phase of the determinant, and Haar sampling by the Gram–Schmidt orthonormalisation of a Ginibre matrix (Mezzadri's phase correction is built into Gram–Schmidt's positive diagonal).
- **The Lorentz group** SO⁺(1, 3) in signature (−, +, +, +) with x⁰ = ct, covered by SL(2, ℂ) through X = x^μ σ_μ ↦ A X A†: Λ(A)^μ_ν = ½ tr(σ_μ A σ_ν A†), with kernel ±I. Two boosts in different directions compose to a boost and a rotation, the Wigner rotation behind Thomas precession, whose angle has a closed form.
- **Recoupling.** Spins are passed as two_j. Clebsch–Gordan coefficients, 3j, 6j and 9j symbols are signed square roots of rationals, computed exactly on morphe's integers (Racah's formulas) and returned as `SignedSqrt`. The 9j symbol is exact too: in its sum over three 6j symbols the triads containing the summation spin appear in pairs, so only the Δ factors of its six rows and columns remain under the root. Sums of commensurable surds are added exactly (`sum`), which turns orthogonality and the Racah and Biedenharn–Elliott identities into exact checks.
- **Quantum 6j symbols** at q = e^{iπ/r}: Racah's formula with quantum factorials [n]! (Kirillov–Reshetikhin), for the spins j ≤ (r − 2)/2 of admissible triads (j₁ + j₂ + j₃ ≤ r − 2); the total quantum dimension w² = Σ [2j+1]² = r / (2 sin²(π/r)). These are the weights of the Turaev–Viro state sum, the finite regularisation of the Ponzano–Regge model of three-dimensional quantum gravity.
- **Weyl integration.** The Haar integral of a class function reduces to the maximal torus with the density |Δ|²; with it, the character expansion f = Σ c_R χ_R of Wilson's Boltzmann weight gives c_n = I_n(β) for U(1) and c_j = 2(2j+1) I_{2j+1}(β)/β for SU(2), the inputs of exact two-dimensional Yang–Mills.

## 3. Content

| Header | Content |
|---|---|
| `u1.hpp` | `U1` (angle in (−π, π]), `wrap_angle`, `character(n, u)`, `haar_u1` |
| `su2.hpp` | `SU2` (a0, a1, a2, a3) with product, inverse, `trace`, `re_trace`, `class_angle`, `matrix`, quaternion sums and `reprojected`; `pauli`, `from_matrix`, `exp_su2`, `rotation_su2`, `log_su2`, `adjoint` (→ SO(3)), `character`, `character_su2(two_j, θ)`, `haar_su2` |
| `sun.hpp` | `SUN` (complex N×N: product, sums, `adjoint`, `trace`, `re_trace`, `determinant`, `reunitarised`, `unitarity_defect`); `generators(N)`, `exp_matrix`, `exp_sun(N, x)`, `haar_sun(N, stream)` |
| `lorentz.hpp` | `Matrix4`, `minkowski_metric`, `boost(φ, n)`, `rotation(n, θ)`, `multiply`, `lorentz_inverse`, `lorentz_defect`; `cover(A)` SL(2, ℂ) → SO⁺(1, 3), `boost_sl2c`; `decompose` (Λ = B R), `rotation_angle`, `wigner_angle` |
| `recoupling.hpp` | `SignedSqrt`, `surd`, `sum` (exact sums of commensurable surds), `isqrt`, `rational_sqrt`; `factorial` (per-thread table), `triangle`, `delta_squared`; `clebsch_gordan`, `wigner_3j`, `wigner_6j`, `wigner_9j` |
| `quantum.hpp` | `quantum_integer`, `quantum_factorial`, `quantum_dimension`, `quantum_dimension_signed`, `total_quantum_dimension_squared`, `admissible`, `quantum_delta`, `quantum_6j(r, …)` |
| `weyl.hpp` | `integrate_u1`, `integrate_su2` (class angle), `integrate_su3` (two eigenphases), `coefficient_u1`, `coefficient_su2` |

## 4. Checks

- `tests/check_groups.cpp` (36 checks): the SU(2) product against its matrices, associativity, inverses, the double cover as a homomorphism against ½ tr(σ_i U σ_j U†), exp and log; characters against their closed form, at ±I, and their Clebsch–Gordan series; Haar moments of SU(2), U(1) and SU(3) by Monte Carlo, within five standard errors; su(N) generators, exp_sun (and its agreement with exp_su2), reunitarisation; the Lorentz group: η preserved, rapidities added, the inverse, the cover of boosts and rotations, the cover as a homomorphism with kernel ±I, the Wigner rotation against its closed form; Weyl integration: orthonormal characters of U(1) and SU(2), character coefficients against modified Bessel functions, the SU(3) moments E|tr U|² = 1 and E|tr U|⁴ = 2.
- `tests/check_recoupling.cpp` (21 checks): known values; exact Clebsch–Gordan and 6j orthogonality; the symmetries of 3j, 6j and 9j symbols; Racah's sum rule and the Biedenharn–Elliott identity, exactly; the 9j symbol against its sum in floating point and its reduction when a spin vanishes; the Ponzano–Regge formula at spins 20–27 and its error halving when the spins double; the total quantum dimension, the orthogonality, the Biedenharn–Elliott identity and the tetrahedral symmetry of quantum 6j symbols, and their classical limit.

## 5. Dependencies

calx `linalg` (the exponential of SU(N)), `random` (Haar sampling), `quadrature` (Weyl integration); morphe `arithmetic` (exact recoupling). In the checks, calx `verification`.

## 6. References

- Racah, G. (1942), Theory of complex spectra II, *Phys. Rev.* 62, 438–462.
- Edmonds, A. R. (1957), *Angular Momentum in Quantum Mechanics*, Princeton University Press.
- Varshalovich, D. A., Moskalev, A. N. and Khersonskii, V. K. (1988), *Quantum Theory of Angular Momentum*, World Scientific.
- Ponzano, G. and Regge, T. (1968), Semiclassical limit of Racah coefficients, in *Spectroscopic and Group Theoretical Methods in Physics*, North-Holland, 1–58.
- Roberts, J. (1999), Classical 6j-symbols and the tetrahedron, *Geom. Topol.* 3, 21–66.
- Kirillov, A. N. and Reshetikhin, N. Yu. (1989), Representations of the algebra U_q(sl(2)), q-orthogonal polynomials and invariants of links, in *Infinite-Dimensional Lie Algebras and Groups*, World Scientific, 285–339.
- Turaev, V. G. and Viro, O. Ya. (1992), State sum invariants of 3-manifolds and quantum 6j-symbols, *Topology* 31, 865–902.
- Mezzadri, F. (2007), How to generate random matrices from the classical compact groups, *Notices AMS* 54, 592–604.
- Weyl, H. (1926), Theorie der Darstellung kontinuierlicher halb-einfacher Gruppen durch lineare Transformationen III, *Math. Z.* 24, 377–395; Bröcker, T. and tom Dieck, T. (1985), *Representations of Compact Lie Groups*, Springer.
- Drouffe, J.-M. and Zuber, J.-B. (1983), Strong coupling and mean field methods in lattice gauge theories, *Phys. Rep.* 102, 1–119.
- Georgi, H. (1999), *Lie Algebras in Particle Physics*, 2nd ed., Westview; Sexl, R. U. and Urbantke, H. K. (2001), *Relativity, Groups, Particles*, Springer.

## 7. Open Questions

1. **Characters of SU(N) for N ≥ 3** by the Weyl character formula, for the character expansion of SU(3) gauge theory; the integral over SU(3) is here, the characters not yet.
2. **Wigner D-functions** and the harmonic analysis on SU(2) and the sphere, shared with physis `analysis`.
3. **Representations of SL(2, ℂ)** (the principal series), for the EPRL spin foam.
4. **The sign conventions of state sums** with quantum 6j symbols (Turaev–Viro): fixed and checked in `quantum_gravity`, by invariance under change of triangulation.
5. **Speed.** The exact 6j symbol costs integers of O(j log j) bits; a cache, or Schulten–Gordon recursion in floating point, may be needed for large spins in state sums.
