# geometry: Metrics, Curvature and Geodesics

*Status: draft (2026-10-28): metrics, curvature, geodesics* ([decision 0005](../development/decisions/0005-symmetry-spacetime-fields-and-gauge.md)); *rotations, moving frames and the frame graph with its holonomy* ([decision 0009](../development/decisions/0009-rotations-frames-and-the-two-body-problem.md)). Euclidean, affine, projective, spherical and ellipsoidal geometry, connections other than Levi-Civita's, differential forms, information geometry and map projections are planned (the ellipsoid's own geodesics and projections are geodesy's, which names the Earth's figure).

## 1. Purpose

The pseudo-Riemannian geometry of manifolds in charts, as the laws of physics are written in it: a metric stated once as code, and its connection, curvature and geodesics computed from it. `spacetime` builds its models on it; `continuum`, `gravitation` and `discretisation` will.

## 2. Concepts

- **A metric field** is a symmetric nondegenerate form g_ab(x) in a chart. Here it is an object with `n` and a templated `operator()`, evaluated on doubles for its components and on calx's jets for their derivatives (calx decision 0005), so that every derivative is exact to rounding and no metric needs its derivatives written by hand. The chart belongs to the type; changes of chart are explicit.
- **Curvature** follows algebraically from g, ∂g and ∂∂g: the Christoffel symbols, the Riemann tensor R^a_bcd = ∂_c Γ^a_db − ∂_d Γ^a_cb + Γ^a_ce Γ^e_db − Γ^a_de Γ^e_cb (the sign convention of decision 0005, §2.4, after Misner, Thorne and Wheeler), the Ricci tensor R_bd = R^a_bad, the scalar, the Einstein tensor, the Weyl tensor and the Kretschmann scalar.
- **Rotations** are unit quaternions (Hamilton's, (w, x, y, z), acting as active rotations), with their matrices, axis and angle, the exponential map from rotation vectors, and Euler angles of the twelve sequences in both conventions (intrinsic, the axes turning with the body; extrinsic, fixed); `frame_rotation` is the passive rotation of SOFA's R₁, R₂, R₃. The double cover by SU(2) is `symmetry`'s, and the conversions agree with it.
- **Moving frames.** A frame's coordinates are carried to another's by a rigid motion and its rate, x_target = R x_source + b and v_target = R v + Ṙ x + ḃ (to second order, with R̈ and b̈). Its size, for comparing motions, is the pair of invariants of Chasles's theorem: the angle of its screw and its slide along the axis.
- **The frame graph.** A world's frames, each joined to its parent by a motion that depends on time, form a tree; the motion between any two is composed along the path through their lowest common ancestor. A relation is a second route between two frames; the holonomy of the loop it closes is the identity when the two routes agree, and its angle and slide say how far they do not.
- **Geodesics** are the solutions of ẍ^a + Γ^a_bc ẋ^b ẋ^c = 0, integrated as a first-order system in (x, u) by calx's one-step methods; g(u, u) and the momenta p_a = g_ab u^b of coordinates the metric does not depend on are conserved, and measure the integration's error.

## 3. Content

| Header | Content |
|---|---|
| `metric.hpp` | `Matrix`, `Point`, `Tensor2/3/4`; the concept `MetricField`; `determinant`, `inverse`, `contract`, `negative_directions` (the inertia of a form); `Euclidean<N>`, `PolarPlane`, `RoundSphere` |
| `curvature.hpp` | `Curvature<n>` (g, g⁻¹, ∂g, Γ, Riemann, Ricci, scalar, Einstein); `christoffel` (first order, by dual numbers), `curvature` (second order, by jets); `lowered_riemann`, `raise_all`, `full_contraction`, `kretschmann`, `weyl`, `einstein_mixed` |
| `geodesic.hpp` | `Phase<n>` (x, u), `position`, `velocity`, `phase`; `geodesic_field`, `norm`, `momentum`, `complete_velocity` (u⁰ from g(u, u) = κ); `geodesic_step` (RK4 or Dormand–Prince), `integrate_geodesic` with an observer |
| `rotation.hpp` | `Vec3`, `Mat3` and their operators (inline namespace `operators`); `rotation_matrix`, `frame_rotation`, `skew`, `vee`; `EulerSequence`, `EulerConvention`; `Rotation` (from and to matrices by Shepperd's method, quaternions, SU(2), Euler angles; `exp`, `log`, `slerp`, `distance`); `orthonormalised` |
| `kinematics.hpp` | `Transform`, `StateTransform`, `KinematicTransform`; `Deviation` and `deviation` (Chasles's angle and slide); `differentiated` |
| `frame_graph.hpp` | `FrameGraph<Time>`: `root`, `add`, `relate`, `find`, `require`, `path`, `transform`, `holonomy` |

## 4. Checks

`tests/check_frames.cpp` (36 checks): rotations' conversions and round trips, Shepperd's method near 0 and π, the double cover against `symmetry`, the Euler sequences at their singular angles, slerp, motions' velocities and accelerations against finite differences, the frame graph's paths and transforms against composition by hand, holonomies of a consistent loop and of a planted disagreement.

`tests/check_geometry.cpp`: inverse, determinant and inertia of small forms; the plane in polar coordinates is flat, with Γ^r_φφ = −r and Γ^φ_rφ = 1/r; the round sphere of radius a has R = 2/a², R^θ_φθφ = sin²θ, a vanishing Einstein tensor and K = 4/a⁴; a great circle closes after one circumference, keeping g(u, u) and p_φ. The curvature of the spacetime models and their geodesics are checked in `spacetime`.

## 5. Dependencies

calx `differentiation` (jets) and `integrators`; physis `symmetry` (SU(2)).

## 6. References

- Misner, C. W., Thorne, K. S. and Wheeler, J. A. (1973), *Gravitation*, Freeman, ch. 8, 11, 13, 14.
- Wald, R. M. (1984), *General Relativity*, University of Chicago Press, ch. 3.
- Carroll, S. M. (2004), *Spacetime and Geometry*, Addison-Wesley, ch. 3.
- O'Neill, B. (1983), *Semi-Riemannian Geometry with Applications to Relativity*, Academic Press.
- Golub, G. H. and Van Loan, C. F. (2013), *Matrix Computations*, 4th ed., Johns Hopkins, §3.4.

## 7. Open Questions

1. Frames (tetrads) and Cartan's structure equations beside coordinate components, for spinors and for numerical relativity.
2. Charts as objects with explicit transition maps, and tensors that carry their chart.
3. Geodesic events (horizons, turning points) located by root finding rather than by the caller.
