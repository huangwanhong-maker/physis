# fields

*Status: draft (2026-10-29), [decision 0012](../development/decisions/0012-fields-as-sections.md).* What a field is, independently of how it is stored, sampled, solved or drawn: a section of a bundle over a base, with its fibre, representation and provenance, and the capabilities of its representation. The theory, with the propositions the checks hold, is [docs/theory/fields.md](../docs/theory/fields.md).

## 1. Concepts

A **field** is a section $\varphi\in\Gamma(E)$ of a bundle over a **base**, with values in a **fibre** that carries a physical dimension. It is held by a **representation** and comes from evidence of a **provenance**. The four are kept apart in its `Descriptor`, because each changes while the others stay: the same base carries a metric and a cochain, the same field is a formula or its samples, and the same values are a model's or an assimilation's.

A query (an instant and a point) is answered by an `Evaluation`: a value, possibly with a note, or a refusal with its reason. The refusal is *outside* (not in the base), *missing* (the data are absent) or *invalid* (the value is ill-defined there), and never a zero.

What a representation can do is a set of **capabilities**, each a C++20 concept on its type: a point evaluation, batches, an exact derivative, pairing with chains, a coboundary, a Hodge star, a fast sampling of a sphere, transfer, versions. They are read off the type, so a kernel keeps its static types. `AnyField` erases them once a batch, for a registry.

## 2. Contents

| Header | Content |
|---|---|
| `field.hpp` | `Query`, `Evaluation`, `Quality`; `Descriptor` with `Base`, `Fibre`, `Representation`, `Provenance`, `Part`; `Capability` and the concepts `Field`, `PointEvaluable`, `BatchEvaluable`, `SpatiallyDifferentiable`, `Versioned`; `capabilities_of`, `described`; `evaluate_many`; `divergence`, `curl` |
| `domains.hpp` | `Cover` and the concept `DomainOf`; `Everywhere`, `Shell`, `Box`, `Interval`, `During`, `AllOf` |
| `analytic.hpp` | `Analytic` (a formula, its exact derivative by calx's dual numbers), `analytic<V, P>()`, `Uniform`, `central_differences` |
| `potential.hpp` | `PotentialField`: B = −∇V of `analysis`'s expansions on the frame's axes, refused inside its reference sphere, read on spheres by rows |
| `reductions.hpp` | the concept `Reduction`; `Magnitude`, `Component`, `Along`, `Radial`, `Horizontal`, `LocalOfSphere`, `OnAxes`, each saying whether it is invariant |
| `compose.hpp` | `reduce`, `sum` (a composite of the most derived provenance), `pullback`, `in_frame` (vectors and Jacobians turned), `restrict` |
| `sampling.hpp` | `SphereGrid`, `PlaneGrid`, `VolumeGrid`, `PointList`; `Samples` (values, qualities, reasons and notes counted); `sample(f, g, t)` |
| `curves.hpp` | `integral_curve`: field lines and streamlines by RK4, each end's stop and reason |
| `history.hpp` | `History`: committed versions in time, a reader holding the version it read |
| `sampled.hpp` | `Sampled`: samples on a volume grid with trilinear reconstruction |
| `any.hpp` | `AnyField<V, P>`, `ScalarField`, `VectorField` |

## 3. Checks

`tests/check_fields.cpp` (49 checks):

- **Formulas.** A dipole's divergence and curl by dual numbers vanish to $10^{-14}$, and central differences agree with them within their derived limit.
- **Potential fields.** An expansion of degree 1 is the dipole's formula, a sphere read by rows equals the expansion point by point to $10^{-13}$, and the divergence and curl of degree 20 vanish by differences. Below the reference sphere the field is refused with its reason, and the potential's gradient is the field.
- **Reductions and compositions.** Their invariance is as declared. Sums give their parts' provenance, and read a sphere by rows when both parts do. A pullback carries a profile on height into space. A change of frame turns vectors and Jacobians.
- **Samplings.** Their orders, and their counts of reasons and notes.
- **Curves.** A dipole's field line lies on $r=L\cos^2\psi$ and ends at the expansion's sphere within $2^{-12}$ of a step. A streamline of a rigid rotation follows RK4's stability function, in radius and in phase, to $10^{-12}$.
- **Versions and samples.** A reader keeps its version while another is committed. Trilinear samples reproduce affine fields exactly.
- **Type erasure.** It keeps the descriptor, the capabilities and the fast path.

## 4. Dependencies

`units`, `geometry` (vectors, matrices, rotations), `analysis` (potential fields), calx `differentiation` (dual numbers). In the checks, calx `random` and `verification`.

## 5. References

- Griewank, A. and Walther, A. (2008), *Evaluating Derivatives*, 2nd ed., SIAM.
- Hairer, E., Nørsett, S. P. and Wanner, G. (1993), *Solving Ordinary Differential Equations I*, 2nd ed., Springer.
- Mac Lane, S. and Moerdijk, I. (1992), *Sheaves in Geometry and Logic*, Springer (sections of bundles, restriction).

## 6. Open Questions

1. **Fields on phase space and on the radiation bundle**, whose bases are not regions of space: their samplings and reductions (moments, integrals over directions).
2. **Batches on devices**: an evaluation of a batch by calx's compute stack rather than a loop.
3. **Uncertainty as a fibre**: a field whose values are distributions, beside its provenance.
