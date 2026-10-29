// Discrete mechanics: dynamics from a discrete action (decision 0013; docs/theory/discrete-mechanics.md).
//
// A mechanical system on coordinates q ∈ ℝᴺ is given by its Lagrangian L(t, q, v), written once as a
// template over its number type, so that it is evaluated on doubles and on calx's jets. Its action over
// a path is approximated step by step: the discrete action S_d = Σ_k L_d(t_k, q_k, q_{k+1}), with a
// discrete Lagrangian L_d ≈ ∫ L dt over a step h by a quadrature rule,
//
//   midpoint    L_d = h L(t_k + h/2, (q_k + q_{k+1})/2, (q_{k+1} − q_k)/h),
//   trapezoid   L_d = (h/2) [L(t_k, q_k, v) + L(t_k + h, q_{k+1}, v)],   v = (q_{k+1} − q_k)/h.
//
// The paths that make S_d stationary, their ends held, satisfy the discrete Euler–Lagrange equations
// D₂L_d(q_{k−1}, q_k) + D₁L_d(q_k, q_{k+1}) = 0. With the discrete Legendre transforms
// p_k = −D₁L_d(q_k, q_{k+1}) and p_{k+1} = D₂L_d(q_k, q_{k+1}) they become a map (q_k, p_k) ↦ (q_{k+1},
// p_{k+1}): solve the first for q_{k+1}, evaluate the second (Marsden and West 2001, Discrete mechanics and
// variational integrators, Acta Numerica 10, 357–514, §1.3–1.6). L_d evaluated on calx jets of order 2 in
// its 2N arguments gives D₁, D₂ and the mixed block D₁₂ that Newton's method needs, exact to rounding:
// the variation of the action is computed, not derived by hand.
//
// What the scheme keeps, from its being the stationary point of a discrete action (Marsden and West,
// §1.3 and §2.6):
//   - the map is symplectic;
//   - Noether's theorem holds discretely: if L_d is invariant under a group acting on q, the momentum
//     map J = ⟨p_k, ξ_Q(q_k)⟩ is conserved exactly by every step, whatever h;
//   - both rules are of order 2, and the energy's error stays bounded rather than drifting.
// In floating point the discrete velocity (q_{k+1} − q_k)/h carries the rounding of the positions over h,
// so a conserved momentum changes by a few ulps of m|q|/h a step: a body far from the origin of its
// coordinates (a rover 1.7 × 10⁶ m from the Moon's centre) is integrated in coordinates about a nearby
// origin, and its place in the world is that origin's frame composed with them.
// The trapezoidal rule applied to L = ½ vᵀ M v − V(q) is Störmer–Verlet (calx integrators), which the
// checks hold to rounding.
//
// Forces that do not come from a potential (a motor, friction, drag) enter by the discrete
// Lagrange–d'Alembert principle: δS_d + Σ (f_d⁻ · δq_k + f_d⁺ · δq_{k+1}) = 0, with the midpoint forces
// f_d^± = (h/2) F(t_k + h/2, (q_k + q_{k+1})/2, (q_{k+1} − q_k)/h), so that
//   p_k + D₁L_d(q_k, q_{k+1}) + f_d⁻ = 0,   p_{k+1} = D₂L_d(q_k, q_{k+1}) + f_d⁺   (Marsden and West, §3, on forcing).
//
// Coupling is in the action. Two systems and their interaction are one Lagrangian, L = L_a + L_b + L_ab,
// whose discrete Euler–Lagrange equations are the coupled equations of motion: the interaction is a term
// of the action, not a force added after each system has moved. A field coupled to a particle (the
// vector potential's q A·v) is such a term.
//
// The same discrete action, its time imaginary, weights the Euclidean path integral: euclidean_action
// sums the trapezoidal L_d of a Euclidean Lagrangian (the kinetic energy plus the potential) round a
// closed path, the exponent of e^{−S_E/ħ} (Feynman and Hibbs 1965, §10; Creutz and Freedman 1981, Ann.
// Phys. 132, 427). For a quadratic Lagrangian it is physis field_theory's free scalar field on the cycle of
// the time steps; the classical motion is the stationary point of the one and the quantum statistics the
// integral over the other.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <calx/differentiation/jet.hpp>
#include <calx/linalg/small.hpp>

namespace physis::mechanics {

namespace cd = calx::differentiation;

template <std::size_t N>
using Coordinates = std::array<double, N>;
template <std::size_t N>
using Block = std::array<std::array<double, N>, N>;

enum class Quadrature { midpoint, trapezoid };

// No force: the default of an integrator.
struct NoForce {
  template <class S, std::size_t N>
  std::array<S, N> operator()(double, const std::array<S, N>&, const std::array<S, N>&) const { return {}; }
};

// A discrete Lagrangian from a Lagrangian, a step and a quadrature rule.
template <std::size_t N, class Lagrangian>
class DiscreteLagrangian {
 public:
  DiscreteLagrangian(Lagrangian L, double h, Quadrature rule = Quadrature::midpoint) : L_(std::move(L)), h_(h), rule_(rule) {
    if (!(h > 0) || !std::isfinite(h)) throw std::invalid_argument("physis::mechanics::DiscreteLagrangian: a positive step");
  }

  double step() const { return h_; }
  Quadrature rule() const { return rule_; }
  const Lagrangian& lagrangian() const { return L_; }

  // L_d(t, q₀, q₁) on any number type.
  template <class S>
  S operator()(double t, const std::array<S, N>& q0, const std::array<S, N>& q1) const {
    std::array<S, N> v;
    for (std::size_t i = 0; i < N; i++) v[i] = (q1[i] - q0[i]) / h_;
    if (rule_ == Quadrature::midpoint) {
      std::array<S, N> mid;
      for (std::size_t i = 0; i < N; i++) mid[i] = (q0[i] + q1[i]) * 0.5;
      return L_(t + 0.5 * h_, mid, v) * h_;
    }
    return (L_(t, q0, v) + L_(t + h_, q1, v)) * (0.5 * h_);
  }

  // Its value and derivatives at (q₀, q₁): D₁, D₂ and D₁₂[i][j] = ∂²L_d/∂q₀ᵢ∂q₁ⱼ.
  struct Derivatives {
    double value = 0;
    Coordinates<N> d1{}, d2{};
    Block<N> d12{};
  };
  Derivatives derivatives(double t, const Coordinates<N>& q0, const Coordinates<N>& q1) const {
    std::array<double, 2 * N> x;
    for (std::size_t i = 0; i < N; i++) x[i] = q0[i], x[N + i] = q1[i];
    const auto z = cd::variables<2 * N, 2>(x);
    std::array<cd::Jet<2 * N>, N> a, b;
    for (std::size_t i = 0; i < N; i++) a[i] = z[i], b[i] = z[N + i];
    const cd::Jet<2 * N> L = (*this)(t, a, b);
    Derivatives d;
    d.value = L.v;
    for (std::size_t i = 0; i < N; i++) {
      d.d1[i] = L.d[i];
      d.d2[i] = L.d[N + i];
      for (std::size_t j = 0; j < N; j++) d.d12[i][j] = L.h[i][N + j];
    }
    return d;
  }

 private:
  Lagrangian L_;
  double h_;
  Quadrature rule_;
};

// A state of the discrete flow: the time, the coordinates and the discrete momenta. Each step adds h to
// t; a clock that must not gather the rounding of a million additions sets t = t₀ + k h itself.
template <std::size_t N>
struct PhaseState {
  double t = 0;
  Coordinates<N> q{}, p{};
};

// How a step's implicit equation was solved.
struct Convergence {
  std::size_t iterations = 0;
  double change = 0;  // the last Newton step, |Δq₁|∞
};

// ---- the continuous Legendre transform ------------------------------------------------------------------

// p = ∂L/∂v at (t, q, v).
template <std::size_t N, class Lagrangian>
Coordinates<N> momentum(const Lagrangian& L, double t, const Coordinates<N>& q, const Coordinates<N>& v) {
  std::array<cd::Dual<N>, N> qs;
  for (std::size_t i = 0; i < N; i++) qs[i] = cd::Dual<N>(q[i]);
  const auto vs = cd::variables<N, 1>(v);
  const auto l = L(t, qs, vs);
  Coordinates<N> p;
  for (std::size_t i = 0; i < N; i++) p[i] = l.d[i];
  return p;
}

// The velocity v with ∂L/∂v(t, q, v) = p, by Newton's method on the Hessian ∂²L/∂v² (a regular
// Lagrangian's, invertible); from v = guess. It stops when the residual is within the tolerance of |p|
// or the correction within it of |v|, whichever comes first, so that neither p = 0 nor v = 0 stalls it.
template <std::size_t N, class Lagrangian>
Coordinates<N> velocity(const Lagrangian& L, double t, const Coordinates<N>& q, const Coordinates<N>& p, Coordinates<N> guess = {},
                        double tolerance = 1e-14, std::size_t max_iterations = 50) {
  std::array<cd::Jet<N>, N> qs;
  for (std::size_t i = 0; i < N; i++) qs[i] = cd::Jet<N>(q[i]);
  Coordinates<N> v = guess;
  for (std::size_t it = 0; it < max_iterations; it++) {
    const auto l = L(t, qs, cd::variables<N, 2>(v));
    Coordinates<N> r;
    double residual = 0, scale = 0;
    for (std::size_t i = 0; i < N; i++) r[i] = p[i] - l.d[i], residual = std::max(residual, std::fabs(r[i])), scale = std::max(scale, std::fabs(p[i]));
    if (residual <= tolerance * scale) return v;
    const Coordinates<N> dv = calx::linalg::solve<N>(l.h, r);
    double change = 0, size = 0;
    for (std::size_t i = 0; i < N; i++) v[i] += dv[i], change = std::max(change, std::fabs(dv[i])), size = std::max(size, std::fabs(v[i]));
    if (change <= tolerance * size) return v;
  }
  throw std::runtime_error("physis::mechanics::velocity: the Legendre transform did not converge (is the Lagrangian regular?)");
}

// The energy E = p · v − L at (t, q, v): the Hamiltonian, conserved when L does not depend on t.
template <std::size_t N, class Lagrangian>
double energy(const Lagrangian& L, double t, const Coordinates<N>& q, const Coordinates<N>& v) {
  const Coordinates<N> p = momentum<N>(L, t, q, v);
  double e = -L(t, q, v);
  for (std::size_t i = 0; i < N; i++) e += p[i] * v[i];
  return e;
}

// The momentum map of a one-parameter group with infinitesimal generator ξ_Q(q): J = ⟨p, ξ_Q(q)⟩.
template <std::size_t N>
double momentum_map(const Coordinates<N>& p, const Coordinates<N>& generator) {
  double j = 0;
  for (std::size_t i = 0; i < N; i++) j += p[i] * generator[i];
  return j;
}

// ---- the variational integrator --------------------------------------------------------------------------

template <std::size_t N, class Lagrangian, class Force = NoForce>
class VariationalIntegrator {
 public:
  VariationalIntegrator(DiscreteLagrangian<N, Lagrangian> Ld, Force F = {}, double tolerance = 1e-13, std::size_t max_iterations = 50)
      : Ld_(std::move(Ld)), F_(std::move(F)), tol_(tolerance), max_(max_iterations) {}

  const DiscreteLagrangian<N, Lagrangian>& discrete_lagrangian() const { return Ld_; }
  double step_size() const { return Ld_.step(); }

  // The state at (t, q) moving with velocity v: its discrete momentum is the continuous one there, the
  // usual start of a variational integrator (Marsden and West, §2.6).
  PhaseState<N> start(double t, const Coordinates<N>& q, const Coordinates<N>& v) {
    last_v_ = v;
    have_v_ = true;
    return {t, q, momentum<N>(Ld_.lagrangian(), t, q, v)};
  }

  // One step: q₁ from p + D₁L_d(q, q₁) + f_d⁻(q, q₁) = 0 by Newton's method from q + h v (v the last
  // step's velocity), then p₁ = D₂L_d(q, q₁) + f_d⁺. Newton stops when its correction is within the
  // tolerance of the larger of |q₁| and the step's displacement |q₁ − q|: the rounding of q₁ is relative to
  // the one, the accuracy of the motion to the other.
  Convergence step(PhaseState<N>& s) {
    const double h = Ld_.step();
    if (!have_v_) last_v_ = velocity<N>(Ld_.lagrangian(), s.t, s.q, s.p), have_v_ = true;
    Coordinates<N> q1;
    for (std::size_t i = 0; i < N; i++) q1[i] = s.q[i] + h * last_v_[i];
    Convergence c;
    for (;;) {
      if (c.iterations >= max_) throw std::runtime_error("physis::mechanics::VariationalIntegrator: the step's implicit equation did not converge");
      const auto d = Ld_.derivatives(s.t, s.q, q1);
      Block<N> J = d.d12;
      Coordinates<N> r;
      for (std::size_t i = 0; i < N; i++) r[i] = -(s.p[i] + d.d1[i]);
      add_force(s.t, s.q, q1, r, J);
      const Coordinates<N> dq = calx::linalg::solve<N>(J, r);
      c.iterations++;
      double change = 0, size = 0;
      for (std::size_t i = 0; i < N; i++) {
        q1[i] += dq[i];
        change = std::max(change, std::fabs(dq[i]));
        size = std::max({size, std::fabs(q1[i]), std::fabs(q1[i] - s.q[i])});
      }
      c.change = change;
      if (change <= tol_ * size) break;
    }
    const auto d = Ld_.derivatives(s.t, s.q, q1);
    Coordinates<N> p1 = d.d2;
    if constexpr (!std::is_same_v<Force, NoForce>) {
      const Coordinates<N> f = force_values(s.t, s.q, q1);
      for (std::size_t i = 0; i < N; i++) p1[i] += f[i];
    }
    for (std::size_t i = 0; i < N; i++) last_v_[i] = (q1[i] - s.q[i]) / h;
    s.q = q1;
    s.p = p1;
    s.t += h;
    return c;
  }

  // The velocity of the last step, (q_{k+1} − q_k)/h.
  const Coordinates<N>& last_velocity() const { return last_v_; }

 private:
  // f_d^± = (h/2) F(t + h/2, (q₀ + q₁)/2, (q₁ − q₀)/h), equal for both ends.
  Coordinates<N> force_values(double t, const Coordinates<N>& q0, const Coordinates<N>& q1) const {
    const double h = Ld_.step();
    Coordinates<N> mid, v;
    for (std::size_t i = 0; i < N; i++) mid[i] = 0.5 * (q0[i] + q1[i]), v[i] = (q1[i] - q0[i]) / h;
    const auto F = F_(t + 0.5 * h, mid, v);
    Coordinates<N> f;
    for (std::size_t i = 0; i < N; i++) f[i] = 0.5 * h * F[i];
    return f;
  }
  // Adds −f_d⁻ to the residual and −∂f_d⁻/∂q₁ to the Jacobian (the residual is −(p + D₁L_d + f_d⁻)).
  void add_force(double t, const Coordinates<N>& q0, const Coordinates<N>& q1, Coordinates<N>& r, Block<N>& J) const {
    if constexpr (!std::is_same_v<Force, NoForce>) {
      const double h = Ld_.step();
      const auto z = cd::variables<N, 1>(q1);
      std::array<cd::Dual<N>, N> mid, v;
      for (std::size_t i = 0; i < N; i++) mid[i] = (z[i] + q0[i]) * 0.5, v[i] = (z[i] - q0[i]) / h;
      const auto F = F_(t + 0.5 * h, mid, v);
      for (std::size_t i = 0; i < N; i++) {
        r[i] -= 0.5 * h * F[i].v;
        for (std::size_t j = 0; j < N; j++) J[i][j] += 0.5 * h * F[i].d[j];
      }
    }
  }

  DiscreteLagrangian<N, Lagrangian> Ld_;
  Force F_;
  double tol_;
  std::size_t max_;
  Coordinates<N> last_v_{};
  bool have_v_ = false;
};

// ---- the Euclidean action ----------------------------------------------------------------------------------

// The Euclidean action of a closed path x₀, …, x_{n−1} (x_n = x₀) with time step a: the trapezoidal
// discrete Lagrangian of the Euclidean Lagrangian L_E(q, v) summed round it,
// S_E = Σ_k (a/2)[L_E(x_k, v_k) + L_E(x_{k+1}, v_k)], v_k = (x_{k+1} − x_k)/a.
template <std::size_t N, class EuclideanLagrangian>
double euclidean_action(const EuclideanLagrangian& LE, double a, const std::vector<Coordinates<N>>& path) {
  const std::size_t n = path.size();
  if (n < 2) throw std::invalid_argument("physis::mechanics::euclidean_action: a closed path of two points or more");
  double s = 0;
  for (std::size_t k = 0; k < n; k++) {
    const auto& x0 = path[k];
    const auto& x1 = path[(k + 1) % n];
    Coordinates<N> v;
    for (std::size_t i = 0; i < N; i++) v[i] = (x1[i] - x0[i]) / a;
    s += 0.5 * a * (LE(x0, v) + LE(x1, v));
  }
  return s;
}

}  // namespace physis::mechanics
