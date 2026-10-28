# 0009. Rotations, Moving Frames and Their Holonomy; the Two-Body Problem

- **Status:** Accepted, provisionally (the frame graph's time is a template parameter, its edges functions of it; a graph of frames that change their parents is not provided)
- **Date:** 2026-10-28
- **Decided by:** the user asked for the spacetime of Cosmopoiesis, "our spacetime model … and the coordinate systems and their conversions, like aircraft and spacecraft navigation systems", and set the rule that a domain writes no general method of its own. Claude wrote the general methods here and this record; the domains' use of them is the infrastructure's [ADR 0049](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0049-the-spacetime-of-the-earth-and-the-moon.md).

## 1. Context

The domains of the Earth and the sky were provisional JavaScript. Their C++ needed, before anything that names a body:

- rotations in their several forms (matrices, unit quaternions, Euler sequences of every convention, the exponential map) with the conversions between them made exact and safe at their singularities;
- rigid motions of space and their rates, so that a velocity carried between rotating frames picks up its ω × r;
- a structure in which a world's frames are joined, the transformation between any two composed along the path between them, and the disagreement of two routes measured;
- the two-body problem, on which the navigation of spacecraft rests: elements and states, Kepler's equation, propagation on any conic, and the frames an orbit carries.

`geometry` held metrics, curvature and geodesics; `mechanics` was planned.

## 2. Decision

### 2.1 `geometry`: Rotations, Motions and the Frame Graph

| Header | Content |
|---|---|
| `rotation.hpp` | `Vec3`, `Mat3` and their operators (in the inline namespace `operators`); `rotation_matrix` (active) and `frame_rotation` (passive, SOFA's Rᵢ); `EulerSequence` (the twelve sequences, intrinsic and extrinsic); `Rotation`, a unit quaternion with `from_matrix` (Shepperd's method), `exp` and `log`, `from_euler` and `euler`, `from_su2` and `su2` (the double cover, physis's `symmetry`), `slerp`, `distance`; `orthonormalised` |
| `kinematics.hpp` | `Transform`, `StateTransform` (a rotation, a translation and their rates: x_target = R x_source + b, and v_target = R v + Ṙ x + ḃ), `KinematicTransform` (to second order); `deviation`, the invariants of a motion by Chasles's theorem (the angle and the slide along its screw axis) |
| `frame_graph.hpp` | `FrameGraph<Time>`: a tree of frames, each joined to its parent by a motion that depends on the time; `transform(from, to, t)` composed along the unique path through the lowest common ancestor (as NAIF's SPICE composes its frames); `relate`, a second route between two frames, closing a loop; `holonomy(relation, t)`, the motion once round it |

*Rationale.* A world's frames form a groupoid: between any two there is one transformation, and composition along any path must agree. A tree makes the composition unique; a relation that closes a loop makes the disagreement of two routes measurable, as the holonomy of the loop, the identity when the routes agree. Its size is stated by the motion's invariants under conjugation (its angle and its slide), not by its translation, which depends on where the loop is based. The navigation bridge uses this to show how far two standards for the same frame disagree: the equinox and CIO routes from the GCRS to the ITRS (well under a milliarcsecond, the celestial pole offsets), the IAU's model of the Moon against DE440's libration (some arcseconds), and GMST of 1982 against IAU 2006's in TEME (tens of milliarcseconds).

### 2.2 `mechanics`: the Two-Body Problem

| Header | Content |
|---|---|
| `two_body.hpp` | `Elements` (p, e, i, Ω, ω, ν, and which of them stand in for undefined ones); `elements(r, v, μ)` and `state(elements, μ)` (Vallado's RV2COE and COE2RV); Kepler's equation on ellipses (Halley's step from Danby's start) and hyperbolas, Barker's equation for parabolas, the anomalies; `propagate(r, v, Δt, μ)` by universal variables and the Stumpff functions; `orbit_frame` (RTN, LVLH, VNC, perifocal) with its rate under the central force or a given acceleration |

The constants of bodies (the Earth's GM, the Moon's) are not here: the domains pass μ.

## 3. Checks

- `geometry/tests/check_frames.cpp` (36 checks): rotations' conversions and their round trips, Shepperd's method at and near θ = 0 and π, the double cover against `symmetry`; the twelve Euler sequences in both conventions, at their singular middle angle; slerp's midpoint on the geodesic; velocities and accelerations in moving frames against finite differences; composition and inversion of motions with their rates; the frame graph's paths and transforms against composition by hand; a consistent loop's holonomy vanishing, and a planted disagreement recovered as its angle and slide.
- `mechanics/tests/check_two_body.cpp` (30 checks): elements against Vallado's own RV2COE (through python-sgp4's copy) at his Example 2-5 and 119 random orbits, to 10⁻¹²; round trips through every kind of orbit; Kepler's equation to rounding; propagation against calx's Dormand–Prince integration of the equations of motion for a low orbit, a Molniya orbit and a hyperbolic flyby; energy and angular momentum kept; the frames' rates against central differences.

## 4. Consequences

- The metrology, geodesy, astronomy and navigation domains state their frames as edges of one `FrameGraph<metrology::time::Instant>`; Cosmopoiesis composes them.
- `mechanics` begins with the two-body problem; Lagrangian and Hamiltonian mechanics, rigid bodies and the n-body problem remain planned, and will sit beside it.
- The JavaScript of the visualization layer interpolates rotations sampled from this C++ (slerp of quaternions) rather than recomputing them.
