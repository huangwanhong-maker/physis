# 0013. Discrete Mechanics, Rigid Bodies, and the Law of a Field Known in Part

- **Status:** Accepted, provisionally (two quadrature rules; rigid bodies without constraints or contact; the covariance of internal potential fields only)
- **Date:** 2026-10-29
- **Decided by:** the user asked, after the unified fields, for physical agents "not distinguished from other physical objects", their dynamics and their coupling in one world, with the question "what about action S based dynamics … the field theory one … path integral even needed?", and "different dynamics system, and even can translate between them"; a long discussion with another assistant was read critically, and its embodied agents, its triple relation of world, bodies and knowledge and its multi-fidelity dynamics were kept, its symbolic variation of actions and its solver framework not. Claude wrote these modules and this record. Their use is the infrastructure's [ADR 0054](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0054-bodies-dynamics-and-knowledge.md): the robotics domain's vehicles and Cosmopoiesis's rover.

## 1. Context

physis's `mechanics` held the two-body problem alone (decision 0009). Bodies in a world need three things from the foundation:
- a way to state their dynamics once and integrate it without losing what it conserves;
- the rotation of a body, and its coupling to the field it moves in;
- the forces that do not come from a potential.

The user's question was whether dynamics should come from an action, as a field theory's does, and whether a path integral is needed. A mission that learns a field from a body's readings also needs the law of the field it does not know: the covariance of a random potential field between any two points.

## 2. Decision

### 2.1 Discrete Mechanics (`mechanics/variational.hpp`)

Dynamics is stated as a Lagrangian, written once over its number type, and integrated as the stationary path of a discrete action ([theory](../../docs/theory/discrete-mechanics.md) §1–3):
- the midpoint and trapezoidal discrete Lagrangians;
- the discrete Euler–Lagrange equations in position–momentum form, solved by Newton's method on the exact $D_{12}$ from calx's jets;
- forces by the discrete Lagrange–d'Alembert principle;
- the continuous Legendre transform and its inverse, the energy;
- the momentum map of a symmetry, conserved exactly by discrete Noether.

Coupling is a term of the action: two bodies and their interaction, or a charge and its vector potential.

The answer to the user's question: dynamics from an action, yes. A path integral is not needed for classical motion. It is the same discrete action used a second way: `euclidean_action` gives the Euclidean action of a closed path, which for a quadratic Lagrangian is `field_theory`'s free lattice action on the cycle of the time steps (theory §5).

### 2.2 Rigid Bodies (`mechanics/rigid_body.hpp`)

- mass properties composed from parts (boxes, points) by Steiner's theorem, with principal axes from calx's eigensolver;
- the Lie group variational integrator of Lee, Leok and McClamroch, its rotation solved in Cayley's parameter (theory §4, Lemma 4.1);
- wrenches from a potential $U(x, R)$, its force and moment by jets.

The free body is Moser and Veselov's discrete top, which keeps its energy exactly (theory, Theorem 4.2). MacCullagh's potential couples orbit and attitude, the total angular momentum conserved to rounding.

### 2.3 The Covariance of a Random Potential Field (`analysis/covariance.hpp`)

The law of decision 0010's random fields, at any two points:
- $\operatorname{Cov}(V(x), V(y)) = a^2\sum_n\sigma_n^2\rho^{n+1}P_n(\hat x\cdot\hat y)$, by Schmidt's addition theorem;
- $\operatorname{Cov}(B(x), B(y))$ as its mixed Hessian, by the chain rule over five sums of Legendre functions and their derivatives, which hold at coincident points;
- `TabulatedCovariance`, the sums at one radius tabulated in the angle, for a magnetometer's readings and a map at its height.

It is checked against the sum over the coefficients to rounding, and the table against the sums.

### 2.4 Dependencies

`mechanics` now uses calx `differentiation` and `linalg` (the small dense solve of calx decision 0012), besides `geometry`. Its checks also use `field_theory`, `integrators` and `random`. `analysis` is unchanged in its dependencies.

## 3. Alternatives Considered

- **Symbolic variation of an action** (an `ActionFunctional::vary` returning equations, as the discussion proposed). Rejected: jets give the variation's values and Jacobian exactly at the points needed, with no algebra system to own.
- **Explicit integrators of the equations of motion** (Runge–Kutta on Hamilton's equations). Kept in calx for what has no structure. A mechanical system is integrated from its action, so that its momenta and its energy behave.
- **Quaternion or exponential coordinates for the rigid body's step.** Rejected for Cayley's parameter, in which the step's equation is rational and its Jacobian exact. The attitude is still kept as physis's quaternion rotation.
- **A general multi-fidelity solver framework.** Not written. The fidelities and the translation between them are the robotics domain's, for its vehicles (ADR 0054). A framework waits until a second science needs one.
- **Kriging in physis.** Rejected: the conditioning of a Gaussian field on observations is statistics. It is tyche's `spatial`; physis gives the covariance.

## 4. Consequences

- `mechanics`'s checks: discrete mechanics (25) and rigid bodies (18), each a closed form or an identity of the discrete theory. `analysis` gains `check_physis_covariance` (6).
- The rounding of a discrete velocity carries the positions' over $h$: a body far from its coordinates' origin is integrated about a nearby origin, as robotics's vehicles do.
- Holonomic constraints in the variational integrator, contact, higher orders and multibody chains are open (theory §6).
