# mechanics: Classical Mechanics

*Status: draft (2026-10-29): the two-body problem (decision 0009); discrete mechanics and rigid bodies ([decision 0013](../development/decisions/0013-discrete-mechanics-and-the-law-of-a-field-known-in-part.md), [theory](../docs/theory/discrete-mechanics.md)).* Constrained discrete Lagrangians, contact, multibody chains, the n-body problem and perturbation theory are planned (README §2).

## 1. Purpose

The motion of bodies under forces, stated so that a domain names only its bodies and its constants: the Earth's GM and its orbits belong to the navigation bridge and to astronomy, the conic a body follows about any centre of attraction to physis.

## 2. Concepts

- **Motion from a discrete action.** A system is its Lagrangian, written once over its number type. Its motion is the stationary path of a discrete action $\sum_k L_d(q_k, q_{k+1})$, by the midpoint or the trapezoidal rule, stepped by the discrete Euler–Lagrange equations with the exact Jacobian of calx's jets. The map is symplectic, and every symmetry's momentum is conserved exactly (discrete Noether). Forces not from a potential enter by the discrete Lagrange–d'Alembert principle, and coupling is a term of the action ([theory §1–3](../docs/theory/discrete-mechanics.md)).
- **The path integral is the same action.** Its time imaginary, the discrete action of a closed path is `field_theory`'s lattice action on the cycle of the time steps, and the lattice's $\langle x^2\rangle$ tends to the quantum oscillator's ([theory §5](../docs/theory/discrete-mechanics.md)).
- **Rigid bodies.** Mass properties are composed from parts by Steiner's theorem. Motion on $\mathbb R^3 \times SO(3)$ is by the Lie group variational integrator, its rotation solved in Cayley's parameter. The free body keeps its energy exactly (Moser and Veselov's discrete top), and a wrench comes from a potential $U(x, R)$ by jets ([theory §4](../docs/theory/discrete-mechanics.md)).
- **The two-body problem.** A point moving about a centre of attraction μ/r² follows a conic, fixed by six classical elements: the semi-latus rectum p, the eccentricity e, the inclination i, the right ascension of the ascending node Ω, the argument of periapsis ω and the true anomaly ν. The elements of a state follow from its angular momentum, its node and its eccentricity vector; where an angle is undefined (a circular orbit has no periapsis, an equatorial one no node) the element taken in its place is named (Vallado 2013, Algorithms 9–10).
- **Kepler's equation** E − e sin E = M on an ellipse and e sinh H − H = M on a hyperbola, solved by Halley's and Newton's steps from Danby's start; a parabola's by Barker's equation in closed form.
- **Propagation by universal variables.** A state is carried over a time on any conic by the universal anomaly χ, the Stumpff functions c₂ and c₃ and the Lagrange coefficients (Vallado 2013, Algorithm 8; Battin 1999, §4.5): one method for ellipses, parabolas and hyperbolas, with no singularity at e = 1.
- **Frames of an orbit.** Right-handed triads moving with the body, as motions to the inertial frame with their rates: RTN (radial, transverse, normal), LVLH (local vertical, local horizontal, as NASA flies the ISS: z to the centre, y against the orbit's normal), VNC (velocity, normal, co-normal) and the perifocal PQW. Their rates follow from differentiating the vectors that define them, under the central force or any acceleration given.

## 3. Content

| Header | Content |
|---|---|
| `variational.hpp` | `DiscreteLagrangian<N, L>` (midpoint, trapezoid; its value, $D_1$, $D_2$, $D_{12}$ by jets), `VariationalIntegrator` (Newton on the forced discrete Euler–Lagrange equations), `momentum`, `velocity` (the Legendre transform and its inverse), `energy`, `momentum_map`, `euclidean_action` |
| `rigid_body.hpp` | `MassProperties`, `solid_box`, `point_mass`, `combine`, `principal_axes`; `RigidBodyState`, `Wrench`, `potential_wrench`, `potential_energy`, `place`, `angular_velocity`, `kinetic_energy`, `angular_momentum`; `RigidBodyIntegrator` (the Lie group variational integrator) |
| `two_body.hpp` | `Elements`, `OrbitKind`; `elements(r, v, μ)`, `state(elements, μ)`, `perifocal_to_inertial`; `eccentric_anomaly`, `hyperbolic_anomaly`, `mean_anomaly`, `true_anomaly` and the conversions between anomalies; `stumpff`, `propagate(r, v, Δt, μ)`; `OrbitFrame`, `orbit_frame(kind, r, v, a, μ)` |

## 4. Checks

`tests/check_variational.cpp` (25 checks):
- the trapezoidal integrator is Störmer–Verlet to rounding;
- the midpoint rule on the oscillator conserves $H$ exactly and turns by $2\arctan(h\omega/2)$;
- the angular momentum of an orbit is conserved to $4\times10^{-13}$ at a step that misses its energy by $10^{-3}$;
- two bodies coupled by an interaction conserve their momentum to rounding;
- two oscillators and their spring follow their discrete normal modes exactly;
- a charge in a magnetic field moves on an exact discrete circle, and a change of gauge shifts its momenta by $e\nabla\chi$ and leaves its motion;
- the damped oscillator and hyperbolic motion converge at order 2;
- the relativistic Legendre transform and its inverse;
- symplecticity and symmetry;
- the Euclidean action is the lattice action, and the lattice's $\langle x^2\rangle$ tends at order $a^2$ to $(\hbar/2m\omega)\coth(\beta\hbar\omega/2)$.

`tests/check_rigid_body.cpp` (18 checks):
- composed boxes are the box, and a dumbbell's inertia is $\mu(|d|^2I - dd^{\mathsf T})$;
- principal axes;
- a potential's wrench against differences;
- the free body conserves $R\Pi$, $|\Pi|$ and its energy to rounding, stays on $SO(3)$, and $\Pi(T)$ converges at order 2 to Euler's equations;
- the symmetric top converges at order 2 to its closed form;
- a small box librates at $\sqrt{3(\mu/r^3)(I_2 - I_1)/I_3}$;
- orbit and attitude conserve their total angular momentum to rounding while exchanging it.

`tests/check_two_body.cpp` (30 checks): elements against Vallado's own RV2COE at his Example 2-5 and 119 random orbits (to 10⁻¹²); round trips through every kind of orbit; Kepler's equation to rounding for e up to 0.999 and on hyperbolas, Barker's equation; propagation against calx's Dormand–Prince integration of the equations of motion (a low orbit, a Molniya orbit, a hyperbolic flyby, to 10⁻⁹ of the radius over a day), energy and angular momentum kept, a period returning the body to its place; the frames' orthonormality and their rates against central differences along the orbit. The reference values are made by `tests/make_fixtures.py` with python-sgp4's copy of Vallado's code.

## 5. Dependencies

physis `geometry` (rotations, motions and their rates); calx `differentiation` (jets) and `linalg` (small dense systems, the symmetric eigenproblem); calx `integrators` and `random`, and physis `field_theory`, in the checks.

## 6. References

- Lee, T., Leok, M. and McClamroch, N. H. (2007). Lie group variational integrators for the full body problem. *Comput. Methods Appl. Mech. Engrg.* 196, 2907–2924.
- Marsden, J. E. and West, M. (2001). Discrete mechanics and variational integrators. *Acta Numerica* 10, 357–514.
- Moser, J. and Veselov, A. P. (1991). Discrete versions of some classical integrable systems and factorization of matrix polynomials. *Commun. Math. Phys.* 139, 217–243.
- Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*, 4th ed. Microcosm Press.
- Battin, R. H. (1999). *An Introduction to the Mathematics and Methods of Astrodynamics*, revised ed. AIAA.
- Danby, J. M. A. (1988). *Fundamentals of Celestial Mechanics*, 2nd ed. Willmann-Bell.
- CCSDS 500.0-G-4 (2019). *Navigation Data — Definitions and Conventions* (local orbital frames).
