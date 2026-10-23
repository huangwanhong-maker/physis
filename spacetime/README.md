# spacetime: Spacetimes as Objects, Continuous and Discrete

*Status: draft (2026-10-23): Minkowski space, exact solutions, Regge geometry, causal sets, emergent dimensions* ([decision 0005](../development/decisions/0005-symmetry-spacetime-fields-and-gauge.md)). Lorentzian Regge calculus, causal dynamical triangulations, the Benincasa–Dowker action, sprinkling into curved spacetimes, Penrose diagrams and horizons as objects are planned.

## 1. Purpose

The arenas of physics: spacetimes as smooth manifolds with exact metrics, and as discrete structures (simplicial complexes with edge lengths, partial orders) that approaches to quantum gravity put in their place. Field theories are written on them (`field_theory`, `gauge`); the dynamics of gravity itself is `gravitation`'s.

## 2. Concepts

- **Minkowski space** ℝ^{1,d−1} with η = diag(−1, 1, …, 1), x⁰ = ct a length; its causal relations (timelike, null, spacelike; future and past), proper time along worldlines, and the volume ζ_d τ^d of an Alexandrov interval. Lorentz transformations are `symmetry`'s.
- **Exact solutions** are metric fields of `geometry`, with parameters as lengths (a mass as GM/c², a charge as √(G/4πε₀) Q/c², a spin as J/Mc): their curvature and geodesics come from `geometry`, and closed forms (horizons, innermost stable orbits, Kretschmann scalars) are given beside them.
- **Regge geometry** (Regge 1961): a simplicial complex whose simplices are flat and whose metric is fixed by its squared edge lengths; curvature sits on the hinges (codimension 2) as deficit angles ε_h = 2π − Σθ, and the action is S = Σ A_h ε_h, with ∂S/∂s_e = Σ (∂A_h/∂s_e) ε_h by Schläfli's identity. Volumes come from the Gram matrix (or the Cayley–Menger determinant), dihedral angles from an embedding of each simplex by its Cholesky factor.
- **Causal sets** (Bombelli, Lee, Meyer and Sorkin 1987): a locally finite partial order, here sprinkled into an Alexandrov interval of Minkowski space by a Poisson process, so that number measures volume and order gives causal structure. The order is held as bitsets of futures and pasts; the ordering fraction gives the Myrheim–Meyer dimension, the longest chain the proper time, and links and interval abundances the input of causal-set actions.

- **Emergent dimensions** (`dimensions.hpp`, [decision 0006](../development/decisions/0006-the-generation-of-spacetime-by-field-theories.md)) of discrete spaces (1-skeletons, dual graphs of top cells, coloured graphs): the spectral dimension from the exact return probability of the lazy random walk, the Hausdorff dimension from the growth of balls. They are the continuum's diagnostics in causal dynamical triangulations, random maps and branched polymers.

## 3. Content

| Header | Content |
|---|---|
| `minkowski.hpp` | `Event<d>`, `Minkowski<d>` (a metric field); `interval`, `Relation`, `relation`, `chronologically_precedes`, `causally_precedes`, `proper_time`; `alexandrov_constant`, `alexandrov_volume` |
| `models.hpp` | `gravitational_length`; `MinkowskiCartesian`, `MinkowskiSpherical`; `Schwarzschild`, `SchwarzschildEddingtonFinkelstein`, `ReissnerNordstrom`, `Kerr` (horizons, `isco`, `kretschmann`); `FLRW<ScaleFactor>` with `PowerLaw` and `Exponential`; `DeSitterStatic`, `DeSitterFlat`, `AntiDeSitterGlobal` |
| `regge.hpp` | `SquaredLengths`, `gram`, `simplex_volume`, `cayley_menger_volume`, `embed`, `dihedral_angle`, `volume_gradient`; `EdgeLengths`; `ReggeComplex` (hinges, interior, `angle_sums`, `deficits`, `hinge_volumes`, `action` with an optional boundary term, `action_gradient`) |
| `causal_set.hpp` | `sprinkle_diamond<d>`, `expected_ordering_fraction`, `myrheim_meyer_dimension`; `CausalSet<d>` (`precedes`, `relations`, `ordering_fraction`, `interval_size`, `links`, `interval_abundances`, `longest_chain`) |
| `dimensions.hpp` | `Graph` (compressed rows), `Graph::from_edges`, `Graph::from`; `skeleton_graph`, `dual_graph`, `coloured_graph`; `return_probability`, `sample_vertices`; `spectral_dimension(g, SpectralOptions)` → `SpectralDimension{t, P, ds}`; `hausdorff_dimension(g, HausdorffOptions)` → `HausdorffDimension{r, volume, dh}` |

## 4. Checks

- `tests/check_models.cpp`: causal relations and the boost-invariance of the interval; the twins' proper times; ζ₂ = 1/2 and ζ₄ = π/24; Kretschmann scalars of Schwarzschild (in both charts, inside the horizon too), Reissner–Nordström and Kerr against their closed forms; Ricci-flatness of the vacuum solutions; a trace-free Weyl tensor; constant curvature R_abcd = ±(g_ac g_bd − g_ad g_bc)/ℓ² of de Sitter and anti-de Sitter; the Friedmann equations from FLRW's Einstein tensor; conformal flatness of FLRW and de Sitter; Kerr's innermost stable orbits at a = 0 and a = M; a circular orbit's period; the advance of the periapsis against the exact quadrature ∫ √(p/(p − 6 − 2e cos χ)) dχ and its weak-field limit 6πM/p; the photon sphere; conservation of energy, angular momentum and norm.
- `tests/check_discrete.cpp`: volumes of regular and random simplices by Gram and Cayley–Menger; dihedral angles of regular simplices and the sine law; ∂V/∂s and Schläfli's identity against central differences; Gauss–Bonnet on the octahedron, the 6-vertex projective plane and the 7-vertex torus; a flat triangulated cube (zero deficits, the boundary term 6πn); the boundary of the 4-simplex and the gradient of its Regge action; the Poisson number of a sprinkling; ordering fractions in two and four dimensions and the Myrheim–Meyer dimension; links and interval abundances against the order; the longest chain in two dimensions against Baik, Deift and Johansson's 2√N − 1.7711 N^{1/6}.

`tests/check_dimensions.cpp`: return probabilities on lattice tori in one to four dimensions against the Fourier sum; d_s and d_H near d; ball counts against exact lattice counts; the self-dual square and cubic lattices; random trees near 4/3 and 2 (finite-size estimates with documented tolerances).

## 5. Dependencies

physis `geometry`; calx `random` (and, in the checks, `quadrature`).

## 6. References

- Misner, C. W., Thorne, K. S. and Wheeler, J. A. (1973), *Gravitation*, Freeman, ch. 27, 31, 33, 42.
- Wald, R. M. (1984), *General Relativity*, University of Chicago Press.
- Carroll, S. M. (2004), *Spacetime and Geometry*, Addison-Wesley, ch. 5–8.
- Griffiths, J. B. and Podolský, J. (2009), *Exact Space-Times in Einstein's General Relativity*, Cambridge University Press.
- Henry, R. C. (2000), Kretschmann scalar for a Kerr–Newman black hole, *ApJ* 535, 350–353.
- Bardeen, J. M., Press, W. H. and Teukolsky, S. A. (1972), Rotating black holes, *ApJ* 178, 347–369.
- Cutler, C., Kennefick, D. and Poisson, E. (1994), Gravitational radiation reaction for bound motion around a Schwarzschild black hole, *Phys. Rev. D* 50, 3816–3835.
- Regge, T. (1961), General relativity without coordinates, *Nuovo Cimento* 19, 558–571.
- Hartle, J. B. and Sorkin, R. (1981), Boundary terms in the action for the Regge calculus, *Gen. Relativ. Gravit.* 13, 541–549.
- Milnor, J. (1994), The Schläfli differential equality, in *Collected Papers* 1, Publish or Perish.
- Bombelli, L., Lee, J., Meyer, D. and Sorkin, R. D. (1987), Space-time as a causal set, *Phys. Rev. Lett.* 59, 521–524.
- Myrheim, J. (1978), Statistical geometry, CERN preprint TH-2538.
- Meyer, D. A. (1988), *The Dimension of Causal Sets*, PhD thesis, MIT.
- Brightwell, G. and Gregory, R. (1991), Structure of random discrete spacetime, *Phys. Rev. Lett.* 66, 260–263.
- Baik, J., Deift, P. and Johansson, K. (1999), On the distribution of the length of the longest increasing subsequence of random permutations, *J. Amer. Math. Soc.* 12, 1119–1178.
- Surya, S. (2019), The causal set approach to quantum gravity, *Living Rev. Relativ.* 22, 5.

## 7. Open Questions

1. Lorentzian Regge calculus (timelike hinges, complex dihedral angles) and its use by causal dynamical triangulations and spin foams.
2. `ReggeComplex` takes plain vertex lists now; it will take morphe's simplicial complexes once `topology` is written.
3. Sprinkling into curved spacetimes (conformally flat ones first: de Sitter, FLRW) and the Benincasa–Dowker action from the interval abundances.
4. Constants G and c here until physis `units` has its versioned CODATA constants.
