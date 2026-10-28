# Glossary

*Status: draft, 2026-10-28.* The sign, orientation and normalisation conventions physis fixes once (CLAUDE.md §6), and the notation of its modules with the names in the code. Created with [decision 0005](../development/decisions/0005-symmetry-spacetime-fields-and-gauge.md); conventions of modules written earlier are added as they are fixed.

## 1. Conventions Fixed Once

| Convention | Choice | Where it shows |
|---|---|---|
| Units | SI internally, radians for angles; quantities read with their units at the boundaries | `units` |
| Metric signature | (−, +, +, +), the timelike coordinate x⁰ = ct first | `geometry`, `spacetime`, `symmetry` (Lorentz) |
| Units of spacetime models | coordinates in metres; a mass as its gravitational length GM/c², a charge as √(G/(4πε₀)) Q/c², in metres | `spacetime` |
| Christoffel symbols | Γ^a_{bc} = ½ g^{ad}(∂_b g_{dc} + ∂_c g_{db} − ∂_d g_{bc}) | `geometry` |
| Riemann tensor | R^a_{bcd} = ∂_c Γ^a_{db} − ∂_d Γ^a_{cb} + Γ^a_{ce} Γ^e_{db} − Γ^a_{de} Γ^e_{cb}; Ricci R_{bd} = R^a_{bad} (Misner, Thorne and Wheeler 1973, "+ + +") | `geometry` |
| Rotations | active, right-handed; `rotation_matrix(axis, θ)` the active rotation, `frame_rotation(axis, θ)` the passive one (SOFA's Rᵢ(θ), which turns the axes by θ and so a fixed vector's coordinates by −θ) | `symmetry`, `geometry` |
| Quaternions | Hamilton's, stored (w, x, y, z), q = cos(θ/2) + sin(θ/2) n acting as the active rotation by θ about n; q and −q are one rotation; `Rotation::su2()` gives the SU(2) element (w, −x, −y, −z) of `symmetry`'s convention | `geometry` |
| Euler sequences | named by their axes in order of application ("ZYX"); intrinsic: each turn about the body's new axis, R = R₁(α) R₂(β) R₃(γ); extrinsic: about the fixed axes, the reverse product; the middle angle in [−π/2, π/2] (Tait–Bryan) or [0, π] (proper) | `geometry` |
| Motions between frames | x_target = R x_source + b, with rates v_target = R v + Ṙ x + ḃ; a `StateTransform` holds (R, Ṙ, b, ḃ); composition (A ∘ B) applies B first; the angular velocity ω = vee(Ṙ Rᵀ), in the target's axes | `geometry` |
| Frame graph | a tree; each frame's edge is the motion from it to its parent; a relation's edge the motion from its first frame to its second by another route; the holonomy of a relation is the motion once round its loop, from its first frame, through the relation and back along the tree | `geometry` |
| Classical orbital elements | p, e, i, Ω, ω, ν (Vallado 2013); the undefined angles of circular or equatorial orbits replaced by u = ω + ν, ϖ = Ω + ω or λ = Ω + ω + ν, which `Elements::kind` names | `mechanics` |
| SU(2) | U = a₀ I + i a·σ, a unit quaternion; su(2) spanned by T_k = −(i/2) σ_k, exp(θ n·T) the active rotation by θ about n; class angle θ with a₀ = cos(θ/2); the cover SU(2) → SO(3) by R(U)_ij = ½ tr(σ_i U σ_j U†) | `symmetry` |
| Spins | passed as twice their value (`two_j`), so that half-integers are exact integers | `symmetry`, `quantum_gravity` |
| Recoupling symbols | Racah's formulas with the Condon–Shortley phase; 6j symbols in the Racah–Wigner normalisation | `symmetry` |
| Quantum integers | [n] = sin(nπ/r)/sin(π/r) at q = e^{iπ/r} | `symmetry`, `quantum_gravity` |
| Lattice theories | dimensionless (lattice units, a = 1; ħ = 1); physical scale enters where a caller sets a lattice spacing or a mass | `field_theory`, `gauge` |
| Gauge transformations | U_e ↦ g_tail U_e g_head⁻¹ | `gauge` |
| Holonomy of a 2-cell | the product of the link variables along its attaching word, from its base: U_f = Π U_e^{±1} | `gauge` |
| Wilson action | S = β Σ_f (1 − Re tr U_f / N) | `gauge` |
| Invariant field theories (decision 0006) | complex fields φ of rank r, Gaussian measure with ⟨φ_a φ̄_b⟩ = δ_ab; an interaction is the trace invariant I_B of its bubble, entering as e^{−λ_B I_B}; a Feynman graph's N-power is its number of faces of colours (0, i), colour 0 the propagators | `field_theory`, `quantum_gravity` |
| Bubbles | bipartite r-coloured graphs, black vertices φ, white vertices φ̄, colour i contracting the i-th indices | morphe `bubbles.hpp` |
| Condensate cosmology | Planck units, ħ = c = 1 with G a parameter; time is relational, the value of a massless scalar field φ | `quantum_gravity` |
| Spherical harmonics (decision 0010) | Schmidt's semi-normalised associated Legendre functions without the Condon–Shortley phase, P̆ₙᵐ (calx `special_functions`); fully normalised P̄ₙᵐ = √(2n + 1) P̆ₙᵐ for gravity; ψ the planetocentric latitude, λ the east longitude | `analysis` |
| Internal potential | V = a Σₙ (a/r)ⁿ⁺¹ Σₘ (gₙᵐ cos mλ + hₙᵐ sin mλ) P̆ₙᵐ(sin ψ), B = −∇V; components north X, east Y, down Z; gₙᵐ, hₙᵐ in the field's unit | `analysis` |
| Fourier transform | to be fixed with the transforms of `analysis` and `discretisation` | — |

## 2. Notation

| Symbol | Meaning | Units | Code |
|---|---|---|---|
| g_{ab}, g^{ab} | metric and its inverse | dimensionless in charts of lengths | `MetricField`, `inverse` |
| Γ^a_{bc}, R^a_{bcd}, R_{ab}, R, G_{ab}, C_{abcd} | Christoffel symbols, Riemann, Ricci, scalar, Einstein, Weyl | m⁻¹, m⁻² | `christoffel`, `curvature` |
| K = R_{abcd} R^{abcd} | Kretschmann scalar | m⁻⁴ | `kretschmann` |
| M, Q, a, ℓ | gravitational length, charge length, spin length, curvature radius | m | `spacetime` models |
| gₙᵐ, hₙᵐ, a | Gauss coefficients of an internal potential, its reference radius | the field's unit; m | `HarmonicPotential` |
| X, Y, Z | north, east, down components of a potential field | the field's unit | `FieldValue::north`, `east`, `down` |
| Rₙ(r) | Lowes–Mauersberger spectrum, the power of degree n at radius r | the field's unit squared | `lowes_spectrum` |
| A nᵝ ρ²ⁿ, (1 − ρ) a | a spectral law and its effective depth | the field's unit squared; m | `SpectralLaw`, `SpectralLaw::depth` |
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
| R, Ṙ, b, ḃ | a motion's rotation, its rate, its translation and its rate | —, s⁻¹, m, m/s | `StateTransform::r, r_dot, b, b_dot` |
| θ, d | the angle of a motion's screw and its slide along the axis (Chasles) | rad, m | `Deviation::angle, slide` |
| p, e, i, Ω, ω, ν | semi-latus rectum, eccentricity, inclination, right ascension of the ascending node, argument of periapsis, true anomaly | m, —, rad | `Elements::p, e, i, raan, argp, nu` |
| E, H, M | eccentric, hyperbolic and mean anomalies | rad | `eccentric_anomaly`, `hyperbolic_anomaly`, `mean_anomaly` |
| χ, ψ, c₂, c₃ | universal anomaly, χ²/a, Stumpff functions | m^½, —, — | `propagate`, `stumpff` |
| μ | the attracting body's gravitational parameter GM | m³ s⁻² | passed by the caller |
