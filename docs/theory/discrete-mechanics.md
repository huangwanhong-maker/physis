# Discrete Mechanics: Motion from a Discrete Action

*Status: draft (2026-10-30). Serves learning and use.* The theory of `mechanics/variational.hpp` and `mechanics/rigid_body.hpp` ([decision 0013](../../development/decisions/0013-discrete-mechanics-and-the-law-of-a-field-known-in-part.md)), and of `mechanics/multibody.hpp` and `mechanics/contact.hpp` ([decision 0014](../../development/decisions/0014-articulated-bodies-in-the-plane-and-compliant-contact.md)), and of `mechanics/lagrangian.hpp` and `mechanics/fluid.hpp` ([decision 0017](../../development/decisions/0017-bodies-in-a-fluid-and-lagrangian-dynamics.md)). Equation numbers in the code refer to the sections here and to the sources cited.

## 1. The Discrete Action and Its Stationary Paths

**Definition 1.1.** Let $Q = \mathbb R^N$ be a configuration space and $L(t, q, v)$ a Lagrangian. A *discrete Lagrangian* $L_d(t_k, q_k, q_{k+1})$ approximates the action $\int_{t_k}^{t_k+h} L\,dt$ along the path between $q_k$ and $q_{k+1}$. Two rules are used, with $v = (q_{k+1} - q_k)/h$:

- the midpoint rule, $L_d = h\,L\big(t_k + \tfrac h2,\ \tfrac{q_k + q_{k+1}}{2},\ v\big)$;
- the trapezoidal rule, $L_d = \tfrac h2\big[L(t_k, q_k, v) + L(t_k + h, q_{k+1}, v)\big]$.

The *discrete action* of a sequence $(q_0, \dots, q_n)$ is $S_d = \sum_k L_d(t_k, q_k, q_{k+1})$.

**Proposition 1.2 (discrete Euler–Lagrange equations).** A sequence with fixed ends makes $S_d$ stationary exactly when, for $0 < k < n$,
$$D_2 L_d(q_{k-1}, q_k) + D_1 L_d(q_k, q_{k+1}) = 0.$$
With the discrete Legendre transforms $p_k = -D_1 L_d(q_k, q_{k+1})$ and $p_{k+1} = D_2 L_d(q_k, q_{k+1})$, these equations are a map $(q_k, p_k) \mapsto (q_{k+1}, p_{k+1})$. One first solves $p_k + D_1 L_d(q_k, q_{k+1}) = 0$ for $q_{k+1}$, then evaluates the second transform (Marsden and West 2001, §1.3–1.6).

*Proof.* $q_k$ appears in two terms of $S_d$. Its variation is the left side, and the transforms match the two terms. $\square$

**Computation.** $L$ is written once over its number type. Evaluated on calx's jets of order 2 in the $2N$ arguments $(q_k, q_{k+1})$, $L_d$ gives $D_1$, $D_2$ and the mixed block $D_{12} = \partial^2 L_d/\partial q_k\partial q_{k+1}$, exact to rounding. $D_{12}$ is the Jacobian of Newton's method for $q_{k+1}$. The variation of the action is computed, never derived by hand.

**Proposition 1.3.** For $L = \tfrac12 v^{\mathsf T} M v - V(q)$:
- the trapezoidal rule's map is Störmer–Verlet;
- the midpoint rule's map is the implicit midpoint rule, which on a quadratic $V$ conserves $H = \tfrac12 p^{\mathsf T}M^{-1}p + V$ exactly and turns each normal mode by $2\arctan(h\omega/2)$.

*Proof.* Trapezoid: $D_1 L_d = -Mv - \tfrac h2\nabla V(q_k)$, so $v = M^{-1}(p_k - \tfrac h2 \nabla V(q_k)) = M^{-1}p_{k+1/2}$. Then $p_{k+1} = Mv - \tfrac h2 \nabla V(q_{k+1})$, which is velocity Verlet.

Midpoint: $p_k = Mv + \tfrac h2\nabla V(\bar q)$ and $p_{k+1} = Mv - \tfrac h2\nabla V(\bar q)$, with $\bar q$ the midpoint. So $q_{k+1} = q_k + hM^{-1}\tfrac{p_k + p_{k+1}}2$ and $p_{k+1} = p_k - h\nabla V(\bar q)$, the implicit midpoint rule. It conserves every quadratic invariant (Hairer, Lubich and Wanner 2006, Theorem IV.2.1). On a linear system it is the Cayley transform, which rotates by $\theta$ with $\tan(\theta/2) = h\omega/2$. $\square$

## 2. Noether's Theorem, Discretely

**Theorem 2.1 (Marsden and West 2001, §1.3).** Suppose $L_d$ is invariant under a one-parameter group acting on $Q$ with generator $\xi_Q$. Then the momentum map $J(q_k, p_k) = \langle p_k, \xi_Q(q_k)\rangle$ is conserved exactly by every step, whatever $h$.

The checks hold it in three cases:
- the angular momentum of an orbit of $e = 0.5$, to $4\times10^{-13}$ over $10^5$ steps of a hundredth of its period, at a step so coarse that the energy is wrong by a part in a thousand;
- the total momentum of two bodies coupled by an anharmonic interaction;
- the uniform motion of their centre of mass.

*Remark (rounding).* The discrete velocity $(q_{k+1} - q_k)/h$ carries the rounding of the positions over $h$. A conserved momentum therefore changes by a few ulps of $m|q|/h$ a step, and the checks measure it so. A body far from the origin of its coordinates is integrated about a nearby origin; robotics's vehicles and the rigid body's centre of mass do so.

**Coupling is in the action.** Two systems and their interaction are one Lagrangian, $L = L_a + L_b + L_{ab}$. Its discrete Euler–Lagrange equations are the coupled equations. A field couples to a particle the same way: a charge in a magnetic field has $L = \tfrac12 m|v|^2 + e\,A(q)\cdot v$.

**Proposition 2.2 (gauge).** Let $A' = A + \nabla\chi$ with $\chi$ quadratic, as between the symmetric gauge $A = \tfrac12 B(-y, x)$ and Landau's $B(-y, 0)$, where $\chi = -\tfrac12 Bxy$. Then the midpoint rule's $L_d$ changes by exactly $e[\chi(q_{k+1}) - \chi(q_k)]$. So the motion is the same, and the momenta differ by $e\nabla\chi(q_k)$.

*Proof.* $\nabla\chi$ is linear, so the midpoint rule integrates $\nabla\chi\cdot\dot q$ along the segment exactly. A discrete total difference has $D_1 + D_2$ equal to the Legendre shifts and changes no equation of motion. $\square$

In a uniform field the discrete orbit is an exact circle, turned by $-2\arctan(h\Omega/2)$ each step with $\Omega = eB/m$. The characteristic equation of the discrete Euler–Lagrange equations has roots of modulus one.

## 3. Forces Not from a Potential

**The discrete Lagrange–d'Alembert principle** (Marsden and West 2001, §3, on forcing) is
$$\delta S_d + \sum_k \big(f_d^-\cdot\delta q_k + f_d^+\cdot\delta q_{k+1}\big) = 0,$$
with the midpoint forces $f_d^\pm = \tfrac h2 F(t_k + \tfrac h2, \bar q, v)$. The forced transforms are $p_k + D_1 L_d + f_d^- = 0$ and $p_{k+1} = D_2 L_d + f_d^+$.

The force's dependence on $q_{k+1}$ enters Newton's Jacobian by dual numbers. The scheme is of order 2: the damped oscillator holds to it.

## 4. Rigid Bodies

**Mass properties.** A body is a union of parts. Masses add, and the centre is the mass-weighted mean. Each part's inertia about its own centre is carried to the body's by Steiner's theorem, $J = \sum_i J_i + m_i(|d_i|^2 I - d_id_i^{\mathsf T})$. A box of sides $a, b, c$ has $J = \tfrac m{12}\operatorname{diag}(b^2+c^2, a^2+c^2, a^2+b^2)$.

**The Lie group variational integrator.** The configuration is $(x, R) \in \mathbb R^3 \times SO(3)$, with $\dot R = R\,\Omega^\wedge$ and $\Pi = J\Omega$. Lee, Leok and McClamroch (2007, §3) take the discrete Lagrangian
$$L_d = \frac m{2h}|x_{k+1} - x_k|^2 + \frac1h\operatorname{tr}\big[(I - F_k)J_d\big] - \frac h2(U_k + U_{k+1}), \qquad F_k = R_k^{\mathsf T}R_{k+1},\quad J_d = \tfrac12\operatorname{tr}(J)\,I - J.$$
Its equations are
$$h\,\Pi_k^\wedge + \tfrac{h^2}2 M_k^\wedge = F_kJ_d - J_dF_k^{\mathsf T}, \qquad R_{k+1} = R_kF_k, \qquad \Pi_{k+1} = F_k^{\mathsf T}\big(\Pi_k + \tfrac h2 M_k\big) + \tfrac h2 M_{k+1},$$
with Störmer–Verlet for $(x, p)$. $R$ stays on $SO(3)$ because it is only ever multiplied by rotations.

**Lemma 4.1.** For symmetric $A$: $a^\wedge A + Aa^\wedge = ((\operatorname{tr}A)\,I - A)a)^\wedge$ and $[a^\wedge, b^\wedge] = (a\times b)^\wedge$. With Cayley's parametrisation
$$F = (I - f^\wedge)^{-1}(I + f^\wedge) = I + c\,(f^\wedge + f^\wedge f^\wedge), \qquad c = \frac2{1+|f|^2},$$
it follows that
$$\operatorname{vee}(FJ_d - J_dF^{\mathsf T}) = c\,(Jf + f\times Jf).$$

*Proof.* $(f^\wedge)^{\mathsf T} = -f^\wedge$ and $(f^{\wedge2})^{\mathsf T} = f^{\wedge2}$, so
$$FJ_d - J_dF^{\mathsf T} = c\,[f^\wedge J_d + J_d f^\wedge + f^{\wedge2}J_d - J_d f^{\wedge2}].$$
The first identity with $A = J_d$, whose $(\operatorname{tr}J_d)I - J_d = J$, gives $(Jf)^\wedge$. Write $f^{\wedge2}J_d - J_df^{\wedge2} = f^\wedge(f^\wedge J_d + J_df^\wedge) - (f^\wedge J_d + J_df^\wedge)f^\wedge$. This is $[f^\wedge, (Jf)^\wedge] = (f\times Jf)^\wedge$. $\square$

The step's rotation therefore solves three smooth rational equations, $c\,(Jf + f\times Jf) = h\,\Pi'$. Newton's method on its exact Jacobian, from $f = \tfrac h2 J^{-1}\Pi'$, meets them in three iterations. The rotation turns through $2\arctan|f|$, and its quaternion is $(1, f)$ normalised.

**Theorem 4.2 (the free body keeps its energy).** Without a moment the step is Moser and Veselov's discrete rigid body (1991). It conserves:
- the spatial angular momentum $R\Pi$ exactly;
- the Casimir $|\Pi|$ exactly;
- the energy $E = \tfrac12\Pi\cdot J^{-1}\Pi$ exactly.

*Proof.* $R_{k+1}\Pi_{k+1} = R_kF_kF_k^{\mathsf T}\Pi_k = R_k\Pi_k$, and $F_k$ is orthogonal. For the energy: Moser and Veselov show that their system keeps the spectral invariants of the Lax matrix $M + \lambda J_d^2$, $M = \Pi^\wedge$, as Manakov's continuous top does. The coefficient of $\lambda$ in $\operatorname{tr}(M + \lambda J_d^2)^3$ gives $\operatorname{tr}(M^2J_d^2) = \Pi^{\mathsf T}J_d^2\Pi - |\Pi|^2\operatorname{tr}J_d^2$.

In three dimensions $J_d$ and $J^{-1} = ((\operatorname{tr}J_d)I - J_d)^{-1}$ are diagonal in one basis, with entries $d_i$ and $1/(s - d_i)$, where $s = \sum d_i$. The condition $d_i^2 = \alpha + \beta/(s - d_i)$ for all $i$ says that the $d_i$ are the roots of $x^3 - sx^2 - \alpha x + (\alpha s + \beta)$. Its second coefficient is $-s = -\sum d_i$ already, and $\alpha$ and $\beta$ fix the other two. So $J_d^2 = \alpha I + \beta J^{-1}$, and $\Pi^{\mathsf T}J_d^2\Pi = \alpha|\Pi|^2 + 2\beta E$. $\square$

The checks see the energy conserved to $4\times10^{-14}$ over $10^5$ steps. The time is the only thing the discrete body gets wrong, at order 2.

**A wrench from a potential.** The force $f = -\partial U/\partial x$ and the moment $M = -\partial U/\partial\eta$, for the variation $R(\eta) = R(I + \eta^\wedge)$, are one function's derivatives, both by jets. When $U(Qx, QR) = U(x, R)$, the total angular momentum $x\times p + R\Pi$ is conserved exactly (Lee et al. 2007, §3; Theorem 2.1 with the rotations of space).

MacCullagh's potential of a body in a central field is
$$U = -\frac{\mu m}r - \frac\mu{2r^3}\big(\operatorname{tr}J - 3\,n^{\mathsf T}Jn\big), \qquad n = R^{\mathsf T}x/r.$$
It couples orbit and attitude:
- a small body librates in pitch at $\sqrt{3(\mu/r^3)(I_2 - I_1)/I_3}$;
- a tumbling one exchanges angular momentum with its orbit while the sum is conserved to rounding.

*Remark.* The libration formula is the limit of a small body. A pitched body feels a tangential force $\propto \theta(I_2 - I_1)/r^4$, which drives its orbit's epicycle near resonance ($\omega_p \approx 1.07\,n$ for the checks' box), and the epicycle drives the pitch. At a size of $3\times10^{-2}r$ this shifts the libration by $2\times10^{-3}$; at $3\times10^{-3}r$, by $10^{-5}$.

## 5. The Euclidean Action and the Path Integral

The same discrete action with imaginary time weights the Euclidean path integral, $Z = \int e^{-S_E/\hbar}\,\mathcal Dx$, over the closed paths of the time cycle $\beta\hbar = na$ (Feynman and Hibbs 1965, §10; Creutz and Freedman 1981).

**Proposition 5.1.** For the oscillator, $L_E = \tfrac12 m\dot x^2 + \tfrac12 m\omega^2x^2$, the trapezoidal $S_E$ of a closed path is $\hbar$ times physis `field_theory`'s free lattice action on the cycle of the time steps. The field is $\varphi = x/\sqrt{a\hbar/m}$ and the lattice mass is $m^2_{\text{lat}} = (a\omega)^2$. So
$$\langle x^2\rangle = \frac{a\hbar}m\,G(0) \longrightarrow \frac\hbar{2m\omega}\coth\frac{\beta\hbar\omega}2 \quad (a \to 0),$$
at order $a^2$.

*Proof.* $S_E = \sum_k\big[m(x_{k+1} - x_k)^2/2a + \tfrac a2 m\omega^2 x_k^2\big]$. Substituting $x = \varphi\sqrt{a\hbar/m}$ gives $\hbar[\tfrac12\sum(\Delta\varphi)^2 + \tfrac12(a\omega)^2\sum\varphi^2]$. $\square$

The classical motion is the stationary point of the one action, and the quantum statistics the integral over it. No other formalism is needed: a path integral is the discrete action used a second way.

## 6. Articulated Bodies in the Plane

**Definition 6.1.** A *planar tree* is a finite set of rigid links $0, \dots, n-1$, each $i > 0$ joined to a parent $p(i) < i$ by a revolute joint. Link $i$ has a mass $m_i$, a centre $c_i$ and an inertia $I_i$ about it in its own frame, and its joint has an *armature* $a_i \ge 0$. The joint of link $i$ sits at the point $r_i$ of its parent's frame. On a *floating base* the coordinates are
$$q = (x, z, \theta_0, q_1, \dots, q_{n-1}),$$
a chart of $SE(2)\times T^{n-1}$: the base's reference point and angle, then the joints' angles. On a *fixed base* the base is held, and $q = (q_1, \dots, q_{n-1})$.

The plane is $(x, z)$ with $z$ up, and an angle turns $x$ toward $z$. A point $X$ of link $i$'s frame is at $o_i + R(\theta_i)X$, where
$$\theta_i = \theta_{p(i)} + q_i + \alpha_i, \qquad o_i = o_{p(i)} + R(\theta_{p(i)})\,r_i,$$
and $\alpha_i$ is the joint's offset. A point's velocity is $\dot o_i + \omega_i\,{\perp}(R(\theta_i)X)$, with ${\perp}(a_x, a_z) = (-a_z, a_x)$ and $\omega_i = \dot\theta_i$.

**The Lagrangian.**
$$L(q, v) = \sum_i\Big(\tfrac12 m_i|\dot c_i|^2 + \tfrac12 I_i\,\omega_i^2 - m_i\,g\,c_{i,z}\Big) + \sum_{i>0}\tfrac12 a_i\,v_i^2 .$$
It is written once over its number type, so that the variational integrator of §1–3 moves the tree.

*Remark (armature).* A joint driven through a gear of ratio $G$ by a rotor of inertia $J_r$ adds the kinetic energy $\tfrac12 J_r(G\dot q_i + \omega_{p(i)})^2$. The armature $a_i = G^2J_r$ keeps its leading term and drops the coupling $J_r\,G\,\dot q_i\,\omega_{p(i)}$, as is usual and small for large $G$ (Featherstone 2008, §9.6).

**Proposition 6.2 (two formulations of one motion).** Let $J_i$ and $J_{\omega i}$ be the Jacobians of $c_i$ and $\omega_i$ in $v$. Then
$$M(q) = \sum_i\big(m_iJ_i^{\mathsf T}J_i + I_iJ_{\omega i}^{\mathsf T}J_{\omega i}\big) + \operatorname{diag}(0, 0, 0, a_1, \dots, a_{n-1}),$$
and the Euler–Lagrange expression $\frac{d}{dt}\partial_vL - \partial_qL$ equals $M(q)\ddot q + h(q, v)$, where $h$ is the generalized force that the recursive Newton–Euler algorithm needs for the acceleration $\ddot q = 0$ (Luh, Walker and Paul 1980; Featherstone 2008, §5.3 and §6.1).

*Proof.* $T = \tfrac12 v^{\mathsf T}Mv$ by the definitions of the Jacobians. The Euler–Lagrange equations of $T - V$ are then $M\ddot q + \dot Mv - \partial_q T + \partial_q V$, which is Newton–Euler's balance of each link's force and moment projected on the coordinates by virtual work. $\square$

The checks hold the two to $1.9\times10^{-16}$ relative on a floating tree of six links. They take the Euler–Lagrange expression from the Lagrangian on calx's jets, so neither is derived by hand. Forward dynamics solves $M\ddot q = Q - h$ by calx's small dense solve.

**Forces.** A force $f$ at a point $P$ of link $i$ has the generalized force $Q = J_P(q)^{\mathsf T}f$ (virtual work), and a joint's torque $\tau$ is $Q_i = \tau$.

**Proposition 6.3 (the free tree's momenta).** Without gravity, the Lagrangian of a floating tree is invariant under translations and rotations of the plane. The translations $x \mapsto x + s$, $z \mapsto z + s$ and the rotation about the origin, generated by $\xi_Q = (-z, x, 1, 0, \dots, 0)$, have the momentum maps $p_x$, $p_z$ and
$$J = -z\,p_x + x\,p_z + p_{\theta_0}.$$
By Theorem 2.1 the variational integrator keeps all three exactly. Under uniform gravity $p_x$ is still kept, and $p_z$ falls by exactly $M_{\text{tot}}\,g\,h$ each step.

*Proof.* A rigid motion of the plane moves every link alike and leaves $|\dot c_i|$, $\omega_i$ and the joints' angles unchanged. The midpoint rule inherits the invariance, since the motion is linear in the coordinates it acts on. The gravity term $-g\sum m_ic_{i,z}$ is linear in $z$, so its midpoint discretisation contributes exactly $-M_{\text{tot}}\,g\,h$ to each step's change of $p_z$. $\square$

The checks see the momenta change by a few ulps of $M|q|/h$ a step (the rounding of §2's remark) and drift by $10^{-10}$ of the largest over 4000 steps. The energy stays within the order-2 error of the step.

**The explicit methods.** On Newton–Euler's equations, $\dot q = v$, $\dot v = M^{-1}(Q - h)$:
- semi-implicit Euler, $v_{k+1} = v_k + h\,M^{-1}(Q - h)(q_k, v_k)$, $q_{k+1} = q_k + h\,v_{k+1}$, of order 1;
- calx's classical Runge–Kutta method, of order 4.

The checks measure the orders on the double pendulum: 1, 2 for the midpoint variational integrator, and 4.

## 7. Contact with a Compliant Ground

**Definition 7.1.** A *ground line* is $z = z_0 + s\,x$. Its outward normal is $n = (-s, 1)/\sqrt{1 + s^2}$ and its tangent $t = (1, s)/\sqrt{1 + s^2}$. A point $p$ moving at $\dot p$ is pressed into it by $\delta = -n\cdot(p - (0, z_0))$ when $\delta > 0$.

**The normal force** (Hunt and Crossley 1975):
$$f_n = k\,\delta^n\,(1 + c\,\dot\delta), \qquad f_n \ge 0 .$$
The exponent $n = 3/2$ is Hertz's for a sphere on an elastic half-space. The damping is proportional to the force, so the force rises from zero at first touch, with no jump as a spring and dashpot would give.

**Proposition 7.2 (restitution).** For small $c\,v_0$, an impact at speed $v_0$ rebounds at $e\,v_0$, with
$$e \approx 1 - \tfrac23\,c\,v_0 .$$

*Proof sketch.* The energy lost in the cycle is $\oint k\,c\,\delta^n\dot\delta\,d\delta$. To first order in $c$ it is evaluated on the elastic cycle, where $\dot\delta^2 = v_0^2 - \tfrac{2k}{m(n+1)}\delta^{n+1}$. With $w = \delta^{n+1}$ each half of the cycle gives $\frac{k}{n+1}\int_0^{w_{\max}}\sqrt{v_0^2 - Aw}\,dw = \tfrac13 m v_0^3$, $A = \frac{2k}{m(n+1)}$, so the loss is $\tfrac23\,c\,m\,v_0^3$ for every $n$ (Hunt and Crossley 1975; Marhefka and Orin 1999). It is also $(1 - e^2)\tfrac12 mv_0^2 \approx (1 - e)\,mv_0^2$, and the result follows. $\square$

So a landing's loss grows with its speed, as a real foot's does. The checks drop a box at 1 m/s and see its rebound within $5\times10^{-3}$ of the estimate.

**Friction.** Coulomb's law $f_t = -\mu f_n\operatorname{sign}(v_t)$, with $v_t = t\cdot\dot p$, is regularised as
$$f_t = -\mu\,f_n\tanh(v_t/v_s).$$
Below the slip speed $v_s$ the ground acts like a stiff viscous damper. A point the friction holds therefore creeps, at a speed of order $v_s$ times the fraction of the friction cone it uses. $v_s$ is the model's error and is stated with it (Stewart 2000, on compliant and rigid models of contact). The checks slide a box at 1 m/s on $\mu = 0.5$ and see it stop within $10^{-2}$ of Coulomb's distance $v_0^2/2\mu g$.

**Proposition 7.3 (rest).** A body of mass $m$ resting on two equal contacts is pressed in by $\delta = (mg/2k)^{1/n}$ at each.

The checks see it to $10^{-6}$ after the impact's oscillation has died.

The force at a point becomes generalized forces by $J_P^{\mathsf T}(f_n n + f_t t)$. Both laws are smooth enough for calx's jets, so the variational integrator takes them as forces (§3) and Newton's method meets each step. Where the step's equation is nonsmooth (a contact's first touch, a $\tanh$ at the slip speed), the robotics domain falls back for one step to semi-implicit Euler and counts it.

*Caveat.* The ground does not move or deform beyond $\delta$: there is no sinkage or bulldozing of a soil (terramechanics), no rolling resistance, and no stiction of a body at rest.

## 8. Bodies in a Fluid

**8.1 The Euler–Lagrange equations solved for the acceleration.** For an autonomous $L(q, v)$, $\frac d{dt}\partial_vL = (\partial_v\partial_vL)\,a + (\partial_v\partial_qL)\,v$. So the equations $\frac d{dt}\partial_vL - \partial_qL = Q$ give
$$(\partial_v\partial_vL)\,a = Q + \partial_qL - (\partial_v\partial_qL)\,v.$$
On jets of order 2 in $(q, v)$ one evaluation of $L$ gives every term exactly (`lagrangian_acceleration`). For a tree's Lagrangian it is Newton–Euler's acceleration: Proposition 6.2, by another road, which the checks hold to $6\times10^{-16}$.

**8.2 The reactive force.** A plate of span $h$ moving across itself at $w_n$ through still fluid sets the fluid moving with the kinetic energy $\tfrac12\,(\rho\pi h^2/4)\,w_n^2$ per unit length. This is the two-dimensional added mass of a plate (Newman 1977, §4.13). Along itself it moves no fluid. For a slender body cut into strips,
$$T_f = \sum_{\text{strips}} \tfrac12\,\frac{\rho\pi h^2}{4}\,ds\;w_n^2, \qquad w = \dot p - u,$$
is a term of the Lagrangian $L = T + T_f - V$. This is Lighthill's (1960) reactive model, strip by strip. Its Euler–Lagrange equations give the added-mass matrix and the forces of its change, and the variational integrator moves the body and the water together. In a uniform, steady current $u$ the term depends on $\dot p - u$ only, so a body at rest in the current feels no reactive force.

**Proposition 8.1.** Without drag, $E = T + T_f$ is conserved.

*Proof.* $L = T + T_f$ is quadratic in $v$ and autonomous, so its energy $v\cdot\partial_vL - L = L$ is the Noether quantity of time translation. $\square$

The checks see it kept to $2\times10^{-13}$ over 2 s by Runge–Kutta on §8.1's equations, for a plate turning and moving.

**8.3 The resistive force.** Each strip feels quadratic drag across itself and skin friction along it:
$$f = -\tfrac12\rho C_n h\,ds\,|w_n|w_n\,n - \tfrac12\rho C_f(2h)\,ds\,|w_t|w_t\,e,$$
with $|w|$ regularised as $\sqrt{w^2 + \varepsilon^2}$, $\varepsilon = 10^{-4}$ m/s. For a plate moving across itself, $(m + m_a)\dot w = -k|w|w$ with $k = \tfrac12\rho C_nh\ell$. Hence $w(t) = w_0/(1 + k|w_0|t/(m + m_a))$, and along itself the same with $m$ alone. The checks hold both to $5\times10^{-8}$. In a current a plate at rest is pushed by $ku^2/(m + m_a)$ and carried toward $u$ hyperbolically, as the same law says.

*Caveat.* Each strip is treated as if alone in a two-dimensional flow. The flow the body has disturbed is not resolved: its wake, the vortices a tail sheds, and the momentum a trailing edge sheds into the wake (Lighthill 1971). A real fish swims mostly by that last, reactive thrust. A body driven by strip theory swims by drag and needs a larger beat for the same speed; the robotics domain's fish measures how much larger.

## 9. What Is Not Yet Here

- Holonomic constraints in the variational integrator. calx's RATTLE holds them for separable systems; a constrained discrete Lagrangian (Marsden and West, §3.4) would hold them for any, closed chains among them.
- Rigid contact and impacts: nonsmooth discrete mechanics (Fetecau, Marsden, Ortiz and West 2003), or complementarity and Moreau's time-stepping (Stewart 2000).
- Higher order: Galerkin discrete Lagrangians and the composition of symmetric steps.
- The spatial tree, $SE(3)\times T^k$, on the rigid body's integrator and Featherstone's spatial vectors, for robotics's arms, suspensions and swimmers.
- A ground that is a field of heights, with its slopes, instead of a line.
- A fluid's wake: Lighthill's large-amplitude elongated-body theory, its trailing-edge momentum flux, or a vortex method.

## References

- Creutz, M. and Freedman, B. (1981), A statistical approach to quantum mechanics, *Ann. Phys.* 132, 427–462.
- Feynman, R. P. and Hibbs, A. R. (1965), *Quantum Mechanics and Path Integrals*, McGraw-Hill.
- Featherstone, R. (2008), *Rigid Body Dynamics Algorithms*, Springer.
- Lighthill, M. J. (1960), Note on the swimming of slender fish, *J. Fluid Mech.* 9, 305–317; (1971), Large-amplitude elongated-body theory of fish locomotion, *Proc. R. Soc. B* 179, 125–138.
- Morison, J. R., O'Brien, M. P., Johnson, J. W. and Schaaf, S. A. (1950), The force exerted by surface waves on piles, *Petroleum Transactions* 189, 149–154.
- Newman, J. N. (1977), *Marine Hydrodynamics*, MIT Press.
- Taylor, G. I. (1952), Analysis of the swimming of long and narrow animals, *Proc. R. Soc. A* 214, 158–183.
- Fetecau, R. C., Marsden, J. E., Ortiz, M. and West, M. (2003), Nonsmooth Lagrangian mechanics and variational collision integrators, *SIAM J. Appl. Dyn. Syst.* 2, 381–416.
- Hairer, E., Lubich, C. and Wanner, G. (2006), *Geometric Numerical Integration*, 2nd ed., Springer.
- Hunt, K. H. and Crossley, F. R. E. (1975), Coefficient of restitution interpreted as damping in vibroimpact, *J. Appl. Mech.* 42, 440–445.
- Lee, T., Leok, M. and McClamroch, N. H. (2007), Lie group variational integrators for the full body problem, *Comput. Methods Appl. Mech. Engrg.* 196, 2907–2924.
- Luh, J. Y. S., Walker, M. W. and Paul, R. P. C. (1980), On-line computational scheme for mechanical manipulators, *J. Dyn. Syst. Meas. Control* 102, 69–76.
- Marhefka, D. W. and Orin, D. E. (1999), A compliant contact model with nonlinear damping for simulation of robotic systems, *IEEE Trans. Syst. Man Cybern. A* 29, 566–572.
- Marsden, J. E. and West, M. (2001), Discrete mechanics and variational integrators, *Acta Numerica* 10, 357–514.
- Moser, J. and Veselov, A. P. (1991), Discrete versions of some classical integrable systems and factorization of matrix polynomials, *Commun. Math. Phys.* 139, 217–243.
- Stewart, D. E. (2000), Rigid-body dynamics with friction and impact, *SIAM Rev.* 42, 3–39.
