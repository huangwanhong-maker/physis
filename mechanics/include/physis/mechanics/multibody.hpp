// Articulated bodies in a plane: a kinematic tree of rigid links joined by revolute joints, on a fixed or a
// floating base (decision 0014; docs/theory/discrete-mechanics.md §6).
//
// Conventions. The plane is (x, z), z up. An angle θ is measured from +x toward +z, so that a positive
// rotation carries x toward z and a moment r × f = r_x f_z − r_z f_x is positive when it increases θ. A
// link's frame has its origin at its joint (the base's at its reference point) and is turned by θ; a point
// X of the link's frame is at o + R(θ) X, R(θ) = [[cos θ, −sin θ], [sin θ, cos θ]]. The perpendicular
// ⊥a = (−a_z, a_x) is the velocity of the point a when the plane turns at unit rate about the origin.
//
// Coordinates. q = (x, z, θ₀, q₁, …, q_{n−1}) on a floating base: the base's reference point and angle,
// then the joints' angles in the order of the links; on a fixed base q = (q₁, …, q_{n−1}) and the base is
// held at its joint's place. A link's parent comes before it; its angle is θ_i = θ_{p(i)} + q_i + α_i
// (α_i its offset), and its joint is at the point r_i of its parent's frame. The configuration space is
// SE(2) × Tⁿ⁻¹ on a floating base, and the coordinates its chart (θ₀ and the q_i taken on ℝ).
//
// Two formulations of one motion, checked against each other to rounding (tests/check_multibody.cpp):
//  - the Lagrangian L(q, v) = Σ ½ m_i |ċ_i|² + ½ I_i ω_i² + Σ ½ a_j v_j² − Σ m_i g c_{i,z}, written once
//    over any number type (doubles, calx jets) so that the discrete action of variational.hpp moves it
//    and its Euler–Lagrange equations follow by differentiation;
//  - Newton–Euler: the recursive Newton–Euler algorithm (Luh, Walker and Paul 1980, J. Dyn. Syst. Meas.
//    Control 102, 69–76; Featherstone 2008, Rigid Body Dynamics Algorithms, §5.3) gives the generalized
//    force Q = M(q) a + h(q, v) that an acceleration a needs: with a = 0 the bias h (Coriolis, centrifugal
//    and gravity), and the mass matrix M(q) = Σ m_i J_iᵀ J_i + I_i J_ωiᵀ J_ωi + diag(a) from the links'
//    Jacobians (Featherstone, §6.1).
// The point Jacobian J_P maps v to a point's velocity; Jᵀ f is the generalized force of a force f there,
// by virtual work. Forward dynamics solves M a = Q − h.
//
// A joint's armature a_j is the inertia its actuator adds on its own coordinate, a gear's ratio squared
// times its rotor's inertia (G² J_r): the rotor's coupling to its parent's rotation, J_r ω_p, is dropped,
// as is usual and small (Featherstone, §9.6).
//
// Explicit integrators of (q, v) for the cheap end of the fidelities: semi-implicit (symplectic) Euler,
// order 1, and calx's classical Runge–Kutta, order 4, on the same Newton–Euler equations. The variational
// integrator of the Lagrangian (variational.hpp) is the structured end: symplectic, of order 2, keeping a
// symmetry's momentum exactly.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <calx/integrators/runge_kutta.hpp>
#include <calx/linalg/small.hpp>
#include <physis/mechanics/variational.hpp>

namespace physis::mechanics {

enum class Base { fixed, floating };

// A rigid link of a planar tree.
struct PlanarLink {
  std::string name;
  int parent = -1;                  // −1 for the base, otherwise an earlier link
  std::array<double, 2> joint{};    // the joint's place in the parent's frame (m); the base's place if fixed
  double offset = 0;                // the link's angle relative to its parent's when its joint is at zero (rad)
  double mass = 0;                  // kg
  std::array<double, 2> centre{};   // the centre of mass in the link's frame (m)
  double inertia = 0;               // about the centre, about the plane's normal (kg m²)
  double armature = 0;              // the actuator's inertia on the joint's coordinate (kg m²)
};

// A link's place and motion: its frame's origin and angle, their rates, its centre of mass and its velocity.
template <class S>
struct PlanarLinkState {
  std::array<S, 2> o{}, odot{};
  S theta{}, omega{};
  std::array<S, 2> c{}, cdot{};
};

template <class S>
inline std::array<S, 2> rotate2(const S& theta, const std::array<double, 2>& x) {
  using std::cos, std::sin;
  const S c = cos(theta), s = sin(theta);
  return {c * x[0] - s * x[1], s * x[0] + c * x[1]};
}
template <class S>
inline std::array<S, 2> perp2(const std::array<S, 2>& a) { return {-a[1], a[0]}; }
inline double cross2(const std::array<double, 2>& a, const std::array<double, 2>& b) { return a[0] * b[1] - a[1] * b[0]; }

template <std::size_t N>
class PlanarTree {
 public:
  PlanarTree(Base base, std::vector<PlanarLink> links, double gravity) : base_(base), links_(std::move(links)), g_(gravity) {
    if (links_.empty()) throw std::invalid_argument("physis::mechanics::PlanarTree: at least the base link");
    if (links_[0].parent != -1) throw std::invalid_argument("physis::mechanics::PlanarTree: the first link is the base");
    for (std::size_t i = 1; i < links_.size(); i++)
      if (links_[i].parent < 0 || std::size_t(links_[i].parent) >= i)
        throw std::invalid_argument("physis::mechanics::PlanarTree: a link's parent comes before it");
    for (const auto& l : links_)
      if (!(l.mass >= 0) || !(l.inertia >= 0) || !(l.armature >= 0)) throw std::invalid_argument("physis::mechanics::PlanarTree: masses and inertias are not negative");
    if (base_dofs() + links_.size() - 1 != N) throw std::invalid_argument("physis::mechanics::PlanarTree: N is the base's coordinates and one a joint");
    ancestors_.assign(links_.size(), {});
    for (std::size_t i = 0; i < links_.size(); i++)
      for (int k = int(i); k > 0; k = links_[std::size_t(k)].parent) ancestors_[i].push_back(std::size_t(k));
  }

  Base base() const { return base_; }
  bool floating() const { return base_ == Base::floating; }
  std::size_t base_dofs() const { return floating() ? 3 : 0; }
  const std::vector<PlanarLink>& links() const { return links_; }
  double gravity() const { return g_; }
  void set_gravity(double g) { g_ = g; }
  // The coordinate of link i's joint (i ≥ 1).
  std::size_t coordinate(std::size_t i) const { return base_dofs() + i - 1; }

  double total_mass() const {
    double m = 0;
    for (const auto& l : links_) m += l.mass;
    return m;
  }

  // ---- kinematics, on any number type ------------------------------------------------------------------
  template <class S>
  std::vector<PlanarLinkState<S>> kinematics(const std::array<S, N>& q, const std::array<S, N>& v) const {
    std::vector<PlanarLinkState<S>> st(links_.size());
    auto& b = st[0];
    if (floating()) {
      if constexpr (N >= 3) {
        b.o = {q[0], q[1]};
        b.odot = {v[0], v[1]};
        b.theta = q[2] + links_[0].offset;
        b.omega = v[2];
      }
    } else {
      b.o = {S(links_[0].joint[0]), S(links_[0].joint[1])};
      b.odot = {S(0.0), S(0.0)};
      b.theta = S(links_[0].offset);
      b.omega = S(0.0);
    }
    centre(b, links_[0]);
    for (std::size_t i = 1; i < links_.size(); i++) {
      const auto& l = links_[i];
      const auto& p = st[std::size_t(l.parent)];
      auto& s = st[i];
      const std::array<S, 2> r = rotate2(p.theta, l.joint);
      const std::array<S, 2> pr = perp2(r);
      s.o = {p.o[0] + r[0], p.o[1] + r[1]};
      s.odot = {p.odot[0] + p.omega * pr[0], p.odot[1] + p.omega * pr[1]};
      const std::size_t k = coordinate(i);
      s.theta = p.theta + q[k] + l.offset;
      s.omega = p.omega + v[k];
      centre(s, l);
    }
    return st;
  }
  template <class S>
  std::vector<PlanarLinkState<S>> kinematics(const std::array<S, N>& q) const {
    std::array<S, N> v;
    v.fill(S(0.0));
    return kinematics(q, v);
  }

  // A point X of link i's frame, its place and its velocity.
  template <class S>
  std::array<S, 2> point(const std::vector<PlanarLinkState<S>>& st, std::size_t i, const std::array<double, 2>& X) const {
    const auto d = rotate2(st[i].theta, X);
    return {st[i].o[0] + d[0], st[i].o[1] + d[1]};
  }
  template <class S>
  std::array<S, 2> point_velocity(const std::vector<PlanarLinkState<S>>& st, std::size_t i, const std::array<double, 2>& X) const {
    const auto d = perp2(rotate2(st[i].theta, X));
    return {st[i].odot[0] + st[i].omega * d[0], st[i].odot[1] + st[i].omega * d[1]};
  }

  // ---- the Lagrangian ---------------------------------------------------------------------------------------
  template <class S>
  S kinetic(const std::array<S, N>& q, const std::array<S, N>& v) const {
    const auto st = kinematics(q, v);
    S T = S(0.0);
    for (std::size_t i = 0; i < links_.size(); i++) {
      const auto& l = links_[i];
      T += (st[i].cdot[0] * st[i].cdot[0] + st[i].cdot[1] * st[i].cdot[1]) * (0.5 * l.mass) + st[i].omega * st[i].omega * (0.5 * l.inertia);
      if (i > 0 && l.armature > 0) T += v[coordinate(i)] * v[coordinate(i)] * (0.5 * l.armature);
    }
    return T;
  }
  template <class S>
  S potential(const std::array<S, N>& q) const {
    const auto st = kinematics(q);
    S V = S(0.0);
    for (std::size_t i = 0; i < links_.size(); i++) V += st[i].c[1] * (links_[i].mass * g_);
    return V;
  }
  template <class S>
  S lagrangian(const std::array<S, N>& q, const std::array<S, N>& v) const { return kinetic(q, v) - potential(q); }

  // ---- Jacobians and generalized forces ----------------------------------------------------------------------
  // The columns ∂p/∂q_k of a point p of link i (its place given): the base's translations, and for the base's
  // angle and each joint on the way from the base to i, the turn ⊥(p − o_k) about that joint.
  template <class S>
  std::array<std::array<S, N>, 2> jacobian(const std::vector<PlanarLinkState<S>>& st, std::size_t i, const std::array<S, 2>& p) const {
    std::array<std::array<S, N>, 2> J;
    for (auto& row : J) row.fill(S(0.0));
    if (floating()) {
      if constexpr (N >= 3) {
        J[0][0] = S(1.0), J[1][1] = S(1.0);
        const auto d = perp2(std::array<S, 2>{p[0] - st[0].o[0], p[1] - st[0].o[1]});
        J[0][2] = d[0], J[1][2] = d[1];
      }
    }
    for (std::size_t k : ancestors_[i]) {
      const auto d = perp2(std::array<S, 2>{p[0] - st[k].o[0], p[1] - st[k].o[1]});
      J[0][coordinate(k)] = d[0], J[1][coordinate(k)] = d[1];
    }
    return J;
  }
  // The generalized force Jᵀ f of a force f at the point p of link i.
  template <class S>
  std::array<S, N> generalized_force(const std::vector<PlanarLinkState<S>>& st, std::size_t i, const std::array<S, 2>& p, const std::array<S, 2>& f) const {
    const auto J = jacobian(st, i, p);
    std::array<S, N> Q;
    for (std::size_t k = 0; k < N; k++) Q[k] = J[0][k] * f[0] + J[1][k] * f[1];
    return Q;
  }

  // ---- Newton–Euler ------------------------------------------------------------------------------------------
  // M(q) = Σ m_i J_iᵀ J_i + I_i J_ωiᵀ J_ωi + diag(armature).
  Block<N> mass_matrix(const Coordinates<N>& q) const {
    const auto st = kinematics(q);
    Block<N> M{};
    for (std::size_t i = 0; i < links_.size(); i++) {
      const auto& l = links_[i];
      const auto J = jacobian(st, i, st[i].c);
      std::array<double, N> Jw{};
      if constexpr (N >= 3)
        if (floating()) Jw[2] = 1;
      for (std::size_t k : ancestors_[i]) Jw[coordinate(k)] = 1;
      for (std::size_t a = 0; a < N; a++)
        for (std::size_t b = 0; b < N; b++) M[a][b] += l.mass * (J[0][a] * J[0][b] + J[1][a] * J[1][b]) + l.inertia * Jw[a] * Jw[b];
    }
    for (std::size_t i = 1; i < links_.size(); i++) M[coordinate(i)][coordinate(i)] += links_[i].armature;
    return M;
  }

  // Q = M(q) a + h(q, v): the recursive Newton–Euler algorithm, gravity included if asked.
  Coordinates<N> inverse_dynamics(const Coordinates<N>& q, const Coordinates<N>& v, const Coordinates<N>& a, bool with_gravity = true) const {
    const std::size_t n = links_.size();
    const auto st = kinematics(q, v);
    std::vector<std::array<double, 2>> oddot(n), cddot(n);
    std::vector<double> alpha(n);
    oddot[0] = {0, 0}, alpha[0] = 0;
    if constexpr (N >= 3)
      if (floating()) oddot[0] = {a[0], a[1]}, alpha[0] = a[2];
    for (std::size_t i = 0; i < n; i++) {
      if (i > 0) {
        const std::size_t p = std::size_t(links_[i].parent);
        const std::array<double, 2> r = {st[i].o[0] - st[p].o[0], st[i].o[1] - st[p].o[1]};
        const auto pr = perp2(r);
        const double w2 = st[p].omega * st[p].omega;
        oddot[i] = {oddot[p][0] + alpha[p] * pr[0] - w2 * r[0], oddot[p][1] + alpha[p] * pr[1] - w2 * r[1]};
        alpha[i] = alpha[p] + a[coordinate(i)];
      }
      const std::array<double, 2> d = {st[i].c[0] - st[i].o[0], st[i].c[1] - st[i].o[1]};
      const auto pd = perp2(d);
      const double w2 = st[i].omega * st[i].omega;
      cddot[i] = {oddot[i][0] + alpha[i] * pd[0] - w2 * d[0], oddot[i][1] + alpha[i] * pd[1] - w2 * d[1]};
    }
    const double g = with_gravity ? g_ : 0.0;
    // Backward: the force F_i and the couple τ_i that link i's parent exerts on it at its joint.
    std::vector<std::array<double, 2>> F(n), accF(n, {0, 0});
    std::vector<double> tau(n), accM(n, 0.0);
    for (std::size_t j = n; j-- > 0;) {
      const auto& l = links_[j];
      F[j] = {l.mass * cddot[j][0] + accF[j][0], l.mass * (cddot[j][1] + g) + accF[j][1]};
      const std::array<double, 2> oc = {st[j].o[0] - st[j].c[0], st[j].o[1] - st[j].c[1]};
      tau[j] = l.inertia * alpha[j] - cross2(oc, F[j]) + accM[j];
      if (j > 0) {
        const std::size_t p = std::size_t(l.parent);
        accF[p][0] += F[j][0], accF[p][1] += F[j][1];
        const std::array<double, 2> rc = {st[j].o[0] - st[p].c[0], st[j].o[1] - st[p].c[1]};
        accM[p] += tau[j] + cross2(rc, F[j]);
      }
    }
    Coordinates<N> Q{};
    if constexpr (N >= 3)
      if (floating()) Q[0] = F[0][0], Q[1] = F[0][1], Q[2] = tau[0];
    for (std::size_t i = 1; i < n; i++) Q[coordinate(i)] = tau[i] + links_[i].armature * a[coordinate(i)];
    return Q;
  }
  Coordinates<N> bias(const Coordinates<N>& q, const Coordinates<N>& v) const { return inverse_dynamics(q, v, Coordinates<N>{}); }

  // a = M⁻¹ (Q − h).
  Coordinates<N> forward_dynamics(const Coordinates<N>& q, const Coordinates<N>& v, const Coordinates<N>& Q) const {
    const auto h = bias(q, v);
    Coordinates<N> r;
    for (std::size_t k = 0; k < N; k++) r[k] = Q[k] - h[k];
    return calx::linalg::solve<N>(mass_matrix(q), r);
  }

  // ---- what the motion keeps ----------------------------------------------------------------------------------
  std::array<double, 2> centre_of_mass(const Coordinates<N>& q) const {
    const auto st = kinematics(q);
    std::array<double, 2> c{0, 0};
    for (std::size_t i = 0; i < links_.size(); i++) c[0] += links_[i].mass * st[i].c[0], c[1] += links_[i].mass * st[i].c[1];
    const double m = total_mass();
    return {c[0] / m, c[1] / m};
  }
  std::array<double, 2> linear_momentum(const Coordinates<N>& q, const Coordinates<N>& v) const {
    const auto st = kinematics(q, v);
    std::array<double, 2> p{0, 0};
    for (std::size_t i = 0; i < links_.size(); i++) p[0] += links_[i].mass * st[i].cdot[0], p[1] += links_[i].mass * st[i].cdot[1];
    return p;
  }
  // About a point a: Σ (c_i − a) × m_i ċ_i + I_i ω_i (the actuators' rotors are not counted: their armature
  // is not a body that turns with the links).
  double angular_momentum(const Coordinates<N>& q, const Coordinates<N>& v, const std::array<double, 2>& about) const {
    const auto st = kinematics(q, v);
    double L = 0;
    for (std::size_t i = 0; i < links_.size(); i++) {
      const std::array<double, 2> r = {st[i].c[0] - about[0], st[i].c[1] - about[1]};
      L += links_[i].mass * cross2(r, st[i].cdot) + links_[i].inertia * st[i].omega;
    }
    return L;
  }
  double energy(const Coordinates<N>& q, const Coordinates<N>& v) const { return kinetic(q, v) + potential(q); }

 private:
  template <class S>
  static void centre(PlanarLinkState<S>& s, const PlanarLink& l) {
    const auto d = rotate2(s.theta, l.centre);
    const auto pd = perp2(d);
    s.c = {s.o[0] + d[0], s.o[1] + d[1]};
    s.cdot = {s.odot[0] + s.omega * pd[0], s.odot[1] + s.omega * pd[1]};
  }

  Base base_;
  std::vector<PlanarLink> links_;
  double g_;
  std::vector<std::vector<std::size_t>> ancestors_;  // the links whose joints move link i: i, its parent, …, the base's child
};

// The tree's Lagrangian as variational.hpp takes it, L(t, q, v), with any further potential V(q) of the
// system (springs, limits) given as a function over any number type.
template <std::size_t N, class ExtraPotential>
struct TreeLagrangian {
  const PlanarTree<N>* tree;
  ExtraPotential extra;
  template <class S>
  S operator()(double, const std::array<S, N>& q, const std::array<S, N>& v) const { return tree->lagrangian(q, v) - extra(q); }
};
struct NoPotential {
  template <class S, std::size_t N>
  S operator()(const std::array<S, N>&) const { return S(0.0); }
};

// ---- explicit integrators of (q, v) ------------------------------------------------------------------------
// Q(t, q, v) is the generalized force of everything but gravity and the tree's inertia (actuators, contact,
// springs written as forces), on doubles.

// Semi-implicit Euler: v₁ = v + h a(q, v), q₁ = q + h v₁. Order 1; symplectic for a separable Lagrangian.
template <std::size_t N, class Force>
void semi_implicit_euler(const PlanarTree<N>& tree, const Force& Q, double& t, Coordinates<N>& q, Coordinates<N>& v, double h) {
  const auto a = tree.forward_dynamics(q, v, Q(t, q, v));
  for (std::size_t k = 0; k < N; k++) v[k] += h * a[k];
  for (std::size_t k = 0; k < N; k++) q[k] += h * v[k];
  t += h;
}

// The classical Runge–Kutta method (calx integrators) on the state (q, v), the force held at the step's
// time as a controller's is (zero-order hold) and evaluated at each stage's state otherwise. Order 4.
template <std::size_t N, class Force>
void runge_kutta4(const PlanarTree<N>& tree, const Force& Q, double& t, Coordinates<N>& q, Coordinates<N>& v, double h) {
  calx::integrators::Vec<2 * N> x;
  for (std::size_t k = 0; k < N; k++) x[k] = q[k], x[N + k] = v[k];
  const double t0 = t;
  calx::integrators::rk4<2 * N>(
      [&](const calx::integrators::Vec<2 * N>& y, calx::integrators::Vec<2 * N>& dy) {
        Coordinates<N> qq, vv;
        for (std::size_t k = 0; k < N; k++) qq[k] = y[k], vv[k] = y[N + k];
        const auto a = tree.forward_dynamics(qq, vv, Q(t0, qq, vv));
        for (std::size_t k = 0; k < N; k++) dy[k] = vv[k], dy[N + k] = a[k];
      },
      x, h);
  for (std::size_t k = 0; k < N; k++) q[k] = x[k], v[k] = x[N + k];
  t += h;
}

}  // namespace physis::mechanics
