// Checks of articulated bodies in a plane and of compliant contact (mechanics/multibody.hpp,
// mechanics/contact.hpp; decision 0014, docs/theory/discrete-mechanics.md §6–7).
//
// - the double pendulum's mass matrix against its closed form, and Newton–Euler's columns of M;
// - the two formulations of one motion: on a floating tree of six links with inertias and armatures, the
//   generalized force M a + h of the recursive Newton–Euler algorithm equals the Euler–Lagrange expression
//   ∂²L/∂v² a + ∂²L/∂v∂q v − ∂L/∂q of the Lagrangian differentiated by calx jets, to rounding;
// - a point's Jacobian against the jets' derivative of its place;
// - discrete Noether on the free tree: the variational integrator keeps the momentum maps of translation
//   and rotation to rounding over thousands of steps, and under gravity the vertical momentum falls by
//   exactly M g h a step;
// - the orders of the three integrators on the double pendulum: semi-implicit Euler 1, the midpoint
//   variational integrator 2, Runge–Kutta 4;
// - compliant contact: a box comes to rest pressed in by (m g / 2k)^{1/n} at each corner; an impact's
//   restitution is Hunt and Crossley's 1 − (2/3) c v₀; a box sliding on friction μ stops in v₀²/2μg, less
//   the regularisation's creep.

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <string>
#include <vector>

#include <calx/differentiation/jet.hpp>
#include <calx/verification/check.hpp>
#include <physis/mechanics/contact.hpp>
#include <physis/mechanics/multibody.hpp>
#include <physis/mechanics/variational.hpp>

using calx::verification::Suite;
namespace pm = physis::mechanics;
namespace cd = calx::differentiation;

namespace {

constexpr double PI = std::numbers::pi;

// A double pendulum of point masses m1, m2 on massless rods l1, l2, hung from a fixed pivot.
pm::PlanarTree<2> double_pendulum(double m1, double m2, double l1, double l2, double g) {
  std::vector<pm::PlanarLink> L(3);
  L[0].name = "pivot";
  L[1].name = "upper", L[1].parent = 0, L[1].mass = m1, L[1].centre = {l1, 0};
  L[2].name = "lower", L[2].parent = 1, L[2].joint = {l1, 0}, L[2].mass = m2, L[2].centre = {l2, 0};
  return pm::PlanarTree<2>(pm::Base::fixed, L, g);
}

// A floating tree of six links like a small quadruped's side: a body, a three-link hind leg and a
// two-link foreleg, with inertias and armatures.
pm::PlanarTree<8> creature(double g) {
  std::vector<pm::PlanarLink> L(6);
  L[0] = {"body", -1, {0, 0}, 0.0, 1.3, {0.01, 0.0}, 0.0108, 0};
  L[1] = {"thigh", 0, {-0.12, -0.02}, -0.5, 0.16, {0.04, 0.0}, 2.1e-4, 8e-4};
  L[2] = {"shank", 1, {0.10, 0.0}, -2.0, 0.08, {0.05, 0.0}, 1.0e-4, 8e-4};
  L[3] = {"foot", 2, {0.12, 0.0}, 2.6, 0.05, {0.04, 0.0}, 4.0e-5, 3e-4};
  L[4] = {"arm", 0, {0.11, -0.03}, -2.0, 0.06, {0.03, 0.0}, 3.0e-5, 2e-4};
  L[5] = {"forearm", 4, {0.07, 0.0}, 1.0, 0.04, {0.04, 0.0}, 2.5e-5, 2e-4};
  return pm::PlanarTree<8>(pm::Base::floating, L, g);
}

template <std::size_t N>
double max_abs(const pm::Coordinates<N>& a) {
  double m = 0;
  for (double x : a) m = std::max(m, std::fabs(x));
  return m;
}

// The Euler–Lagrange expression Q = ∂²L/∂v² a + ∂²L/∂v∂q v − ∂L/∂q, by jets of order 2 in (q, v).
template <std::size_t N>
pm::Coordinates<N> euler_lagrange(const pm::PlanarTree<N>& tree, const pm::Coordinates<N>& q, const pm::Coordinates<N>& v, const pm::Coordinates<N>& a) {
  std::array<double, 2 * N> x;
  for (std::size_t i = 0; i < N; i++) x[i] = q[i], x[N + i] = v[i];
  const auto z = cd::variables<2 * N, 2>(x);
  std::array<cd::Jet<2 * N>, N> qq, vv;
  for (std::size_t i = 0; i < N; i++) qq[i] = z[i], vv[i] = z[N + i];
  const auto L = tree.lagrangian(qq, vv);
  pm::Coordinates<N> Q{};
  for (std::size_t i = 0; i < N; i++) {
    double s = -L.d[i];
    for (std::size_t j = 0; j < N; j++) s += L.h[N + i][N + j] * a[j] + L.h[N + i][j] * v[j];
    Q[i] = s;
  }
  return Q;
}

struct Free {
  template <class S, std::size_t N>
  std::array<S, N> operator()(double, const std::array<S, N>&, const std::array<S, N>&) const {
    std::array<S, N> f;
    f.fill(S(0.0));
    return f;
  }
};

// A box (one floating link) on the ground, touching at its two lower corners.
struct BoxOnGround {
  const pm::PlanarTree<3>* tree;
  pm::GroundLine ground;
  pm::CompliantContact law;
  std::array<std::array<double, 2>, 2> corners;
  template <class S>
  std::array<S, 3> operator()(double, const std::array<S, 3>& q, const std::array<S, 3>& v) const {
    const auto st = tree->kinematics(q, v);
    std::array<S, 3> Q;
    Q.fill(S(0.0));
    for (const auto& X : corners) {
      const auto p = tree->point(st, 0, X);
      const auto c = pm::contact_force(ground, law, p, tree->point_velocity(st, 0, X));
      if (!c.touching) continue;
      const auto f = tree->generalized_force(st, 0, p, c.force);
      for (std::size_t k = 0; k < 3; k++) Q[k] += f[k];
    }
    return Q;
  }
};

}  // namespace

int main() {
  Suite s;

  // ---- the double pendulum's mass matrix --------------------------------------------------------------------
  {
    const double m1 = 1.3, m2 = 0.7, l1 = 0.9, l2 = 0.6, g = 9.81;
    const auto tree = double_pendulum(m1, m2, l1, l2, g);
    double worst = 0, worst_columns = 0;
    for (double q2 : {-2.0, -0.4, 0.3, 1.7}) {
      const pm::Coordinates<2> q = {0.6, q2};
      const auto M = tree.mass_matrix(q);
      const double c = std::cos(q2);
      const double closed[2][2] = {{(m1 + m2) * l1 * l1 + m2 * l2 * l2 + 2 * m2 * l1 * l2 * c, m2 * l2 * l2 + m2 * l1 * l2 * c},
                                   {m2 * l2 * l2 + m2 * l1 * l2 * c, m2 * l2 * l2}};
      for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++) worst = std::max(worst, std::fabs(M[i][j] - closed[i][j]));
      for (std::size_t k = 0; k < 2; k++) {
        pm::Coordinates<2> e{};
        e[k] = 1;
        const auto col = tree.inverse_dynamics(q, {0, 0}, e, false);
        for (std::size_t i = 0; i < 2; i++) worst_columns = std::max(worst_columns, std::fabs(col[i] - M[i][k]));
      }
    }
    s.check("double pendulum: the mass matrix against its closed form (kg m²)", worst, 1e-14);
    s.check("double pendulum: Newton–Euler with a unit acceleration and no velocity or gravity gives M's columns (kg m²)", worst_columns, 1e-14);
  }

  // ---- two formulations: Newton–Euler and the Lagrangian's Euler–Lagrange equations ----------------------
  {
    const auto tree = creature(1.62);
    double worst = 0, scale = 0;
    const pm::Coordinates<8> qs[3] = {{0.3, 0.25, 0.1, -0.4, -2.1, 2.5, -2.2, 1.1},
                                      {-1.2, 0.4, -0.35, 0.2, -1.4, 2.9, -1.6, 0.6},
                                      {0.0, 0.15, 0.8, -1.0, -2.5, 2.0, -2.6, 1.6}};
    const pm::Coordinates<8> vs[3] = {{0.5, -0.3, 1.2, -3.0, 2.0, 4.0, -1.0, 2.5},
                                      {-1.5, 2.0, -0.4, 6.0, -5.0, 1.0, 3.0, -2.0},
                                      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
    const pm::Coordinates<8> as[3] = {{1.0, -2.0, 0.5, 30.0, -20.0, 10.0, 5.0, -8.0},
                                      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                      {-3.0, 4.0, 2.0, -1.0, 7.0, -6.0, 2.0, 1.0}};
    for (int k = 0; k < 3; k++) {
      const auto ne = tree.inverse_dynamics(qs[k], vs[k], as[k]);
      const auto el = euler_lagrange(tree, qs[k], vs[k], as[k]);
      for (std::size_t i = 0; i < 8; i++) worst = std::max(worst, std::fabs(ne[i] - el[i])), scale = std::max(scale, std::fabs(el[i]));
    }
    s.check("two formulations: on a floating tree of six links, Newton–Euler's M a + h equals the Lagrangian's Euler–Lagrange expression by jets (relative)",
            worst / scale, 1e-13);

    // Forward dynamics inverts it.
    const auto a = tree.forward_dynamics(qs[0], vs[0], tree.inverse_dynamics(qs[0], vs[0], as[0]));
    double back = 0;
    for (std::size_t i = 0; i < 8; i++) back = std::max(back, std::fabs(a[i] - as[0][i]));
    s.check("forward dynamics: M⁻¹ (Q − h) returns the acceleration Q was made from (rad/s², m/s²)", back, 1e-10);

    // The point Jacobian against the jets' derivative of the point's place.
    const auto st = tree.kinematics(qs[0]);
    const std::array<double, 2> toe = {0.09, 0.0};
    const auto J = tree.jacobian(st, 3, tree.point(st, 3, toe));
    const auto z = cd::variables<8, 1>(qs[0]);
    const auto stj = tree.kinematics(z);
    const auto pj = tree.point(stj, 3, toe);
    double dj = 0;
    for (std::size_t k = 0; k < 8; k++) dj = std::max({dj, std::fabs(J[0][k] - pj[0].d[k]), std::fabs(J[1][k] - pj[1].d[k])});
    s.check("jacobian: a point's columns ⊥(p − o_k) against the jets' derivative of its place (m)", dj, 1e-15);
  }

  // ---- discrete Noether on a free tree ------------------------------------------------------------------------
  {
    const auto tree = creature(0.0);
    using L = pm::TreeLagrangian<8, pm::NoPotential>;
    const double h = 2e-3;
    pm::VariationalIntegrator<8, L> vi(pm::DiscreteLagrangian<8, L>(L{&tree, {}}, h));
    pm::Coordinates<8> q0 = {0.1, 0.3, 0.2, -0.4, -2.1, 2.5, -2.2, 1.1}, v0 = {0.4, -0.2, 2.0, -6.0, 5.0, 8.0, -4.0, 3.0};
    auto st = vi.start(0, q0, v0);
    const auto rot = [](const pm::PhaseState<8>& x) { return -x.q[1] * x.p[0] + x.q[0] * x.p[1] + x.p[2]; };
    // A momentum map changes by rounding alone: the discrete velocity (q₁ − q₀)/h carries the rounding of
    // the coordinates over h, so a step's change is a few ulps of M|q|/h (variational.hpp), and over
    // thousands of steps it wanders like their sum.
    double step_p = 0, step_J = 0, drift_p = 0, drift_J = 0, E0 = tree.energy(q0, v0), dE = 0;
    const double px0 = st.p[0], pz0 = st.p[1], J0 = rot(st);
    for (int k = 0; k < 4000; k++) {
      const auto before = st;
      vi.step(st);
      double qmax = 0;
      for (std::size_t i = 0; i < 8; i++) qmax = std::max(qmax, std::fabs(st.q[i]));
      const double ulp_p = 2.2e-16 * tree.total_mass() * qmax / h, ulp_J = ulp_p * (std::fabs(st.q[0]) + std::fabs(st.q[1]) + 1);
      step_p = std::max({step_p, std::fabs(st.p[0] - before.p[0]) / ulp_p, std::fabs(st.p[1] - before.p[1]) / ulp_p});
      step_J = std::max(step_J, std::fabs(rot(st) - rot(before)) / ulp_J);
      drift_p = std::max({drift_p, std::fabs(st.p[0] - px0), std::fabs(st.p[1] - pz0)});
      drift_J = std::max(drift_J, std::fabs(rot(st) - J0));
      dE = std::max(dE, std::fabs(tree.energy(st.q, vi.last_velocity()) - E0) / E0);
    }
    double psize = 0;
    for (double x : st.p) psize = std::max(psize, std::fabs(x));
    s.check("noether: the free tree's discrete linear momentum changes a step by rounding alone (in ulps of M|q|/h)", step_p, 8);
    s.check("noether: its discrete angular momentum about the origin, −z pₓ + x p_z + p_θ, likewise (in ulps of M|q|(1 + |x| + |z|)/h)", step_J, 8);
    s.check("noether: over 4000 steps the momenta drift by rounding alone (relative to the largest momentum)", std::max(drift_p, drift_J) / psize, 1e-10);
    s.check("noether: while its energy, not a momentum map of the step, stays within its order-2 error (relative)", dE, 2e-3);

    // Under gravity the vertical momentum falls by exactly M g h a step.
    const auto heavy = creature(1.62);
    pm::VariationalIntegrator<8, L> vg(pm::DiscreteLagrangian<8, L>(L{&heavy, {}}, h));
    auto sg = vg.start(0, q0, v0);
    double fall = 0;
    for (int k = 0; k < 500; k++) {
      const double before = sg.p[1];
      vg.step(sg);
      fall = std::max(fall, std::fabs((before - sg.p[1]) - heavy.total_mass() * 1.62 * h));
    }
    s.check("noether: under gravity the discrete vertical momentum falls by M g h each step (relative to M g h)", fall / (heavy.total_mass() * 1.62 * h), 1e-10);
  }

  // ---- the orders of the three integrators -------------------------------------------------------------------
  {
    const auto tree = double_pendulum(1.0, 0.8, 0.7, 0.5, 9.81);
    const pm::Coordinates<2> q0 = {1.0, -0.5}, v0 = {0.0, 0.0};
    const double T = 0.5;
    const auto none = [](double, const pm::Coordinates<2>&, const pm::Coordinates<2>&) { return pm::Coordinates<2>{}; };
    const auto run = [&](int method, double h) {
      double t = 0;
      pm::Coordinates<2> q = q0, v = v0;
      const int n = int(std::lround(T / h));
      if (method == 2) {
        using L = pm::TreeLagrangian<2, pm::NoPotential>;
        pm::VariationalIntegrator<2, L> vi(pm::DiscreteLagrangian<2, L>(L{&tree, {}}, h));
        auto st = vi.start(0, q0, v0);
        for (int k = 0; k < n; k++) vi.step(st);
        return st.q;
      }
      for (int k = 0; k < n; k++) {
        if (method == 1) pm::semi_implicit_euler(tree, none, t, q, v, h);
        else pm::runge_kutta4(tree, none, t, q, v, h);
      }
      return q;
    };
    const auto ref = run(4, T / 8192);
    const auto order = [&](int method, double h) {
      const auto e1 = run(method, h), e2 = run(method, h / 2);
      const double a = std::max(std::fabs(e1[0] - ref[0]), std::fabs(e1[1] - ref[1]));
      const double b = std::max(std::fabs(e2[0] - ref[0]), std::fabs(e2[1] - ref[1]));
      return std::log2(a / b);
    };
    s.check("orders: semi-implicit Euler on the double pendulum is of order 1 (the measured order − 1)", order(1, T / 256) - 1, 0.1);
    s.check("orders: the midpoint variational integrator, of order 2", order(2, T / 64) - 2, 0.1);
    s.check("orders: Runge–Kutta, of order 4", order(4, T / 32) - 4, 0.15);
  }

  // ---- compliant contact ----------------------------------------------------------------------------------------
  {
    // A box of 2 kg, 0.2 m by 0.1 m, comes to rest on two corners.
    const double m = 2.0, g = 9.81;
    std::vector<pm::PlanarLink> L(1);
    L[0] = {"box", -1, {0, 0}, 0.0, m, {0, 0}, m * (0.04 + 0.01) / 12, 0};
    const pm::PlanarTree<3> box(pm::Base::floating, L, g);
    pm::CompliantContact law;
    law.stiffness = 2e6, law.exponent = 1.5, law.damping = 3.0, law.friction = 0.5, law.slip_speed = 0.005;
    const BoxOnGround Q{&box, {0.0, 0.0}, law, {{{-0.1, -0.05}, {0.1, -0.05}}}};
    double t = 0;
    pm::Coordinates<3> q = {0, 0.051, 0}, v = {0, 0, 0};
    for (int k = 0; k < 40000; k++) pm::semi_implicit_euler(box, Q, t, q, v, 5e-5);
    const double expected = std::pow(m * g / 2 / law.stiffness, 1 / law.exponent);
    s.check("contact: a box at rest is pressed in by (m g / 2k)^{1/n} at each corner (relative)", (0.05 - q[1]) / expected - 1, 1e-6);

    // An impact without gravity at v₀ = 1 m/s: its restitution against Hunt and Crossley's 1 − (2/3) c v₀.
    const pm::PlanarTree<3> drifting(pm::Base::floating, L, 0.0);
    pm::CompliantContact soft = law;
    soft.damping = 0.06, soft.friction = 0.0;
    const BoxOnGround Qi{&drifting, {0.0, 0.0}, soft, {{{-0.1, -0.05}, {0.1, -0.05}}}};
    t = 0, q = {0, 0.0505, 0}, v = {0, -1.0, 0};
    for (int k = 0; k < 20000 && !(q[1] > 0.0505 && v[1] > 0); k++) pm::runge_kutta4(drifting, Qi, t, q, v, 1e-6);
    s.check("contact: an impact at 1 m/s rebounds with Hunt and Crossley's restitution 1 − (2/3) c v₀", v[1] - (1 - 2.0 / 3 * soft.damping * 1.0), 5e-3);

    // A box sliding at 1 m/s on μ = 0.5 stops in v₀²/2μg (0.102 m), the regularisation's creep apart.
    t = 0, q = {0, 0.05 - expected, 0}, v = {1.0, 0, 0};
    double x_stop = 0;
    for (int k = 0; k < 200000; k++) {
      pm::semi_implicit_euler(box, Q, t, q, v, 2e-5);
      if (v[0] < 0.02) { x_stop = q[0]; break; }
    }
    const double coulomb = 1.0 / (2 * law.friction * g) * (1 - 0.02 * 0.02);
    s.check("contact: a box sliding at 1 m/s on μ = 0.5 slows to 0.02 m/s within Coulomb's distance (relative)", x_stop / coulomb - 1, 0.01);
  }

  return s.done();
}
