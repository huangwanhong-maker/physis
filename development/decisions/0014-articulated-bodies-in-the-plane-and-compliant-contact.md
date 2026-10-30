# 0014. Articulated Bodies in the Plane, and Compliant Contact

- **Status:** Accepted, provisionally (planar trees of revolute joints; contact compliant, with a ground line)
- **Date:** 2026-10-30
- **Decided by:** the user asked for robots that could be built in the world, "has a theoretical physics base", with "multiple the dynamics coupling, as we said we may has many methods, like action S based", and chose a rabbit for the Moon; the representation was discussed (a kinematic graph of links and joints, the state in SE(2) × Tᵏ first, the geometry attached to the links). Claude wrote these modules and this record. Their use is the infrastructure's [ADR 0058](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0058-a-body-with-joints.md): the robotics domain's `legged` and Cosmopoiesis's ROBOT page.

## 1. Context

Decision 0013 gave `mechanics` motion from a discrete action and the free rigid body. A robot with legs is several rigid bodies joined, which touch the ground and leave it. The foundation lacked three things:
- the kinematic tree: its coordinates, its kinetic energy written once over any number type, and its equations of motion computed in more than one way, so that each method of integration can be held to the others;
- explicit integrators of the same equations, for the cheap end of the fidelities, beside the variational integrator;
- contact with the ground that a variational integrator can take as a force.

The theory note listed "multibody chains" and "contact and impacts" as not yet here (its §6, now §9).

## 2. Decision

### 2.1 Planar Trees (`mechanics/multibody.hpp`)

A tree of rigid links joined by revolute joints, on a fixed base or a floating one ([theory](../../docs/theory/discrete-mechanics.md) §6):
- coordinates $q = (x, z, \theta_0, q_1, \dots, q_{n-1})$ on a floating base, a chart of $SE(2)\times T^{n-1}$;
- the kinematics of the links and of any point on them, and its Jacobian;
- the Lagrangian $L = T - V$, with each joint's armature, written once over its number type, so that the variational integrator of decision 0013 moves the tree;
- Newton–Euler: the recursive algorithm for the generalized force an acceleration needs, the mass matrix from the links' Jacobians, the forward dynamics by calx's small dense solve;
- the centre of mass, the linear and angular momenta, the energy;
- semi-implicit Euler and Runge–Kutta (calx's) on Newton–Euler's equations.

The two formulations are one motion: Newton–Euler's $M\ddot q + h$ equals the Euler–Lagrange expression of the Lagrangian by jets, to rounding.

### 2.2 Compliant Contact (`mechanics/contact.hpp`)

A ground line $z = z_0 + s\,x$ and a law of contact at a point (theory §7):
- the normal force of Hunt and Crossley, $f_n = k\,\delta^n(1 + c\,\dot\delta) \ge 0$;
- Coulomb's friction regularised by $\tanh(v_t/v_s)$.

Both are written over any number type, so the variational integrator takes them as forces by Lagrange–d'Alembert, and the explicit integrators on doubles.

### 2.3 Dependencies

`mechanics` now uses calx `integrators` (the classical Runge–Kutta method), besides `differentiation` and `linalg`.

## 3. Alternatives Considered

- **Three dimensions at once**, $SE(3)\times T^k$ with spatial vectors (Featherstone's). Deferred: a hopping robot in its plane shows every part of the problem (a floating base, joints, contact, a controller) with a third of the coordinates, and its checks are closed forms. The spatial tree follows the planar one when a robot needs it, the fish first.
- **Only the variational integrator.** Rejected: a robot's controller is designed and run on cheap explicit models, and the comparison of methods on one motion is part of what is learnt. The variational integrator stays the structured end.
- **Rigid contact** (complementarity, Moreau's time-stepping, or nonsmooth variational collision integrators). Deferred: compliance is what a real foot on regolith has, its parameters are measurable, and it is smooth enough for Newton's method. Rigid contact is a later choice of the fidelities.
- **The articulated-body algorithm** (Featherstone's $O(n)$ forward dynamics). Not needed at six links: the mass matrix and a dense solve cost less than the controller.
- **A general robot description format** (URDF). Not adopted as the representation: the tree is a list of links with parents, joints and inertias, which a URDF reader can fill later.

## 4. Consequences

- `mechanics`'s checks gain `check_physis_multibody` (16): the double pendulum's mass matrix against its closed form; the two formulations equal to $1.9\times10^{-16}$; the forward dynamics' inverse; the Jacobian against jets; discrete Noether on a free tree (the momenta to rounding over 4000 steps); the fall under gravity by $Mgh$ a step; the orders 1, 2 and 4 of the three integrators; a box at rest pressed in by $(mg/2k)^{1/n}$, an impact's restitution against Hunt and Crossley's, a slide stopping within Coulomb's distance.
- The robotics domain's rabbit is built on them (ADR 0058).
- Open (theory §9): the spatial tree; rigid contact and impacts; holonomic loops (a closed chain); a ground that is a field of heights, not a line.
