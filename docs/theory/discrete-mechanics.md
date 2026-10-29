# Discrete Mechanics: Motion from a Discrete Action

*Status: draft (2026-10-29). Serves learning and use.* The theory of `mechanics/variational.hpp` and `mechanics/rigid_body.hpp` ([decision 0013](../../development/decisions/0013-discrete-mechanics-and-the-law-of-a-field-known-in-part.md)). Equation numbers in the code refer to the sections here and to the sources cited.

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

## 6. What Is Not Yet Here

- Holonomic constraints in the variational integrator. calx's RATTLE holds them for separable systems; a constrained discrete Lagrangian (Marsden and West, §3.4) would hold them for any.
- Contact and impacts: nonsmooth discrete mechanics (Fetecau, Marsden, Ortiz and West 2003).
- Higher order: Galerkin discrete Lagrangians and the composition of symmetric steps.
- Multibody chains on the rigid body's integrator, for robotics's arms and suspensions.

## References

- Creutz, M. and Freedman, B. (1981), A statistical approach to quantum mechanics, *Ann. Phys.* 132, 427–462.
- Feynman, R. P. and Hibbs, A. R. (1965), *Quantum Mechanics and Path Integrals*, McGraw-Hill.
- Fetecau, R. C., Marsden, J. E., Ortiz, M. and West, M. (2003), Nonsmooth Lagrangian mechanics and variational collision integrators, *SIAM J. Appl. Dyn. Syst.* 2, 381–416.
- Hairer, E., Lubich, C. and Wanner, G. (2006), *Geometric Numerical Integration*, 2nd ed., Springer.
- Lee, T., Leok, M. and McClamroch, N. H. (2007), Lie group variational integrators for the full body problem, *Comput. Methods Appl. Mech. Engrg.* 196, 2907–2924.
- Marsden, J. E. and West, M. (2001), Discrete mechanics and variational integrators, *Acta Numerica* 10, 357–514.
- Moser, J. and Veselov, A. P. (1991), Discrete versions of some classical integrable systems and factorization of matrix polynomials, *Commun. Math. Phys.* 139, 217–243.
