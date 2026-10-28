# mechanics: Classical Mechanics

*Status: draft (2026-10-28): the two-body problem.* Lagrangian and Hamiltonian mechanics, rigid-body rotation, the n-body problem and perturbation theory are planned (README §2).

## 1. Purpose

The motion of bodies under forces, stated so that a domain names only its bodies and its constants: the Earth's GM and its orbits belong to the navigation bridge and to astronomy, the conic a body follows about any centre of attraction to physis.

## 2. Concepts

- **The two-body problem.** A point moving about a centre of attraction μ/r² follows a conic, fixed by six classical elements: the semi-latus rectum p, the eccentricity e, the inclination i, the right ascension of the ascending node Ω, the argument of periapsis ω and the true anomaly ν. The elements of a state follow from its angular momentum, its node and its eccentricity vector; where an angle is undefined (a circular orbit has no periapsis, an equatorial one no node) the element taken in its place is named (Vallado 2013, Algorithms 9–10).
- **Kepler's equation** E − e sin E = M on an ellipse and e sinh H − H = M on a hyperbola, solved by Halley's and Newton's steps from Danby's start; a parabola's by Barker's equation in closed form.
- **Propagation by universal variables.** A state is carried over a time on any conic by the universal anomaly χ, the Stumpff functions c₂ and c₃ and the Lagrange coefficients (Vallado 2013, Algorithm 8; Battin 1999, §4.5): one method for ellipses, parabolas and hyperbolas, with no singularity at e = 1.
- **Frames of an orbit.** Right-handed triads moving with the body, as motions to the inertial frame with their rates: RTN (radial, transverse, normal), LVLH (local vertical, local horizontal, as NASA flies the ISS: z to the centre, y against the orbit's normal), VNC (velocity, normal, co-normal) and the perifocal PQW. Their rates follow from differentiating the vectors that define them, under the central force or any acceleration given.

## 3. Content

| Header | Content |
|---|---|
| `two_body.hpp` | `Elements`, `OrbitKind`; `elements(r, v, μ)`, `state(elements, μ)`, `perifocal_to_inertial`; `eccentric_anomaly`, `hyperbolic_anomaly`, `mean_anomaly`, `true_anomaly` and the conversions between anomalies; `stumpff`, `propagate(r, v, Δt, μ)`; `OrbitFrame`, `orbit_frame(kind, r, v, a, μ)` |

## 4. Checks

`tests/check_two_body.cpp` (30 checks): elements against Vallado's own RV2COE at his Example 2-5 and 119 random orbits (to 10⁻¹²); round trips through every kind of orbit; Kepler's equation to rounding for e up to 0.999 and on hyperbolas, Barker's equation; propagation against calx's Dormand–Prince integration of the equations of motion (a low orbit, a Molniya orbit, a hyperbolic flyby, to 10⁻⁹ of the radius over a day), energy and angular momentum kept, a period returning the body to its place; the frames' orthonormality and their rates against central differences along the orbit. The reference values are made by `tests/make_fixtures.py` with python-sgp4's copy of Vallado's code.

## 5. Dependencies

physis `geometry` (rotations, motions and their rates); calx `integrators` and `random` in the checks.

## 6. References

- Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*, 4th ed. Microcosm Press.
- Battin, R. H. (1999). *An Introduction to the Mathematics and Methods of Astrodynamics*, revised ed. AIAA.
- Danby, J. M. A. (1988). *Fundamentals of Celestial Mechanics*, 2nd ed. Willmann-Bell.
- CCSDS 500.0-G-4 (2019). *Navigation Data — Definitions and Conventions* (local orbital frames).
