# Glossary

*Status: draft, 2026-10-23.* The sign, orientation and normalisation conventions physis fixes once (CLAUDE.md §6), and the notation of its modules with the names in the code. Created with [decision 0005](../development/decisions/0005-symmetry-spacetime-fields-and-gauge.md); conventions of modules written earlier are added as they are fixed.

## 1. Conventions Fixed Once

| Convention | Choice | Where it shows |
|---|---|---|
| Units | SI internally, radians for angles; quantities read with their units at the boundaries | `units` |
| Metric signature | (−, +, +, +), the timelike coordinate x⁰ = ct first | `geometry`, `spacetime`, `symmetry` (Lorentz) |
| Units of spacetime models | coordinates in metres; a mass as its gravitational length GM/c², a charge as √(G/(4πε₀)) Q/c², in metres | `spacetime` |
| Christoffel symbols | Γ^a_{bc} = ½ g^{ad}(∂_b g_{dc} + ∂_c g_{db} − ∂_d g_{bc}) | `geometry` |
| Riemann tensor | R^a_{bcd} = ∂_c Γ^a_{db} − ∂_d Γ^a_{cb} + Γ^a_{ce} Γ^e_{db} − Γ^a_{de} Γ^e_{cb}; Ricci R_{bd} = R^a_{bad} (Misner, Thorne and Wheeler 1973, "+ + +") | `geometry` |
| Rotations | active, right-handed | `symmetry`, `geometry` |
| SU(2) | U = a₀ I + i a·σ, a unit quaternion; su(2) spanned by T_k = −(i/2) σ_k, exp(θ n·T) the active rotation by θ about n; class angle θ with a₀ = cos(θ/2); the cover SU(2) → SO(3) by R(U)_ij = ½ tr(σ_i U σ_j U†) | `symmetry` |
| Spins | passed as twice their value (`two_j`), so that half-integers are exact integers | `symmetry`, `quantum_gravity` |
| Recoupling symbols | Racah's formulas with the Condon–Shortley phase; 6j symbols in the Racah–Wigner normalisation | `symmetry` |
| Quantum integers | [n] = sin(nπ/r)/sin(π/r) at q = e^{iπ/r} | `symmetry`, `quantum_gravity` |
| Lattice theories | dimensionless (lattice units, a = 1; ħ = 1); physical scale enters where a caller sets a lattice spacing or a mass | `field_theory`, `gauge` |
| Gauge transformations | U_e ↦ g_tail U_e g_head⁻¹ | `gauge` |
| Holonomy of a 2-cell | the product of the link variables along its attaching word, from its base: U_f = Π U_e^{±1} | `gauge` |
| Wilson action | S = β Σ_f (1 − Re tr U_f / N) | `gauge` |
| Fourier transform, spherical harmonics | to be fixed with `analysis` and `discretisation` | — |

## 2. Notation

| Symbol | Meaning | Units | Code |
|---|---|---|---|
| g_{ab}, g^{ab} | metric and its inverse | dimensionless in charts of lengths | `MetricField`, `inverse` |
| Γ^a_{bc}, R^a_{bcd}, R_{ab}, R, G_{ab}, C_{abcd} | Christoffel symbols, Riemann, Ricci, scalar, Einstein, Weyl | m⁻¹, m⁻² | `christoffel`, `curvature` |
| K = R_{abcd} R^{abcd} | Kretschmann scalar | m⁻⁴ | `kretschmann` |
| M, Q, a, ℓ | gravitational length, charge length, spin length, curvature radius | m | `spacetime` models |
| ε_h, θ | deficit angle at a hinge, dihedral angle | rad | `regge.hpp` |
| j, two_j, m | spin, twice the spin, magnetic number | — | `recoupling.hpp` |
| {j₁ j₂ j₃; j₄ j₅ j₆} | 6j symbol | — | `wigner_6j`, `quantum_6j` |
| χ_j(θ) | character of the spin-j representation of SU(2) | — | `character_su2` |
| c_R(β) | character-expansion coefficient of the Boltzmann weight | — | `gauge` exact |
| β | inverse gauge coupling | — | `gauge` |
| U_e, U_f | link variable, plaquette (2-cell) holonomy | — | `gauge` |
| w² = Σ_j [2j+1]² | total quantum dimension | — | `total_quantum_dimension_squared` |
| ω | Gurau degree of a coloured graph | — | morphe `ColouredGraph::degree` |
| ε_g(n) | Harer–Zagier numbers: gluings of a 2n-gon into genus g | — | `quantum_gravity` |
