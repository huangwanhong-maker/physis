// Checks of rigid bodies (mechanics/rigid_body.hpp; decision 0013, docs/theory/discrete-mechanics.md §4).
//
// - Mass properties: a box cut into 24 boxes and recombined by the parallel axis theorem is the box; a
//   dumbbell's inertia is μ (|d|² I − d dᵀ); the principal axes of a rotated tensor are the rotation.
// - potential_wrench: its force and moment against differences of the potential.
// - The free body, its inertia not diagonal: the spatial angular momentum R Π, the Casimir |Π| and the
//   energy are conserved to rounding, the attitude stays on SO(3), and Π(T) converges at order 2 to
//   Dormand–Prince on Euler's equations. The energy is exact because the free flow is Moser and Veselov's
//   discrete rigid body (1991, Commun. Math. Phys. 139, 217–243), which keeps tr M² and tr M² J_d², M = Π^:
//   that is |Π|² and Πᵀ J_d² Π, and in three dimensions J_d² = α I + β J⁻¹ (the d_i are the roots of a cubic
//   whose first coefficient tr J_d already fixes; Vieta), so Πᵀ J_d² Π = α |Π|² + 2β E.
// - The symmetric top against its closed form R(t) = exp(t L^/I₁) R₀ exp(t (1/I₃ − 1/I₁) Π₃ e₃^), at order 2.
// - Orbit and attitude coupled by MacCullagh's potential: a small body's pitch librates at
//   √(3 (μ/r³)(I₂ − I₁)/I₃) on the circular orbit of the whole potential,
//   and a tumbling body on an eccentric orbit exchanges angular momentum with it while the total is
//   conserved to rounding.

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <string>
#include <vector>

#include <calx/integrators/runge_kutta.hpp>
#include <calx/verification/check.hpp>
#include <physis/mechanics/rigid_body.hpp>

using calx::verification::Suite;
namespace pm = physis::mechanics;
namespace pg = physis::geometry;
using namespace physis::geometry::operators;
using pg::Mat3;
using pg::Rotation;
using pg::Vec3;

namespace {

double largest(const Mat3& a) {
  double m = 0;
  for (const auto& r : a)
    for (double x : r) m = std::max(m, std::fabs(x));
  return m;
}
double largest(const Vec3& a) { return std::max({std::fabs(a[0]), std::fabs(a[1]), std::fabs(a[2])}); }

Mat3 conjugate(const Mat3& Q, const Vec3& d) {
  const Mat3 D{{{d[0], 0, 0}, {0, d[1], 0}, {0, 0, d[2]}}};
  return Q * D * pg::transpose(Q);
}

// MacCullagh's potential of a body of mass m and inertia J (body axes) at x from a point mass μ/G:
// U = −μ m/r − (μ/2r³)(tr J − 3 nᵀ J n), n = Rᵀ x / r, over any number type.
struct MacCullagh {
  double mu, m;
  Mat3 J;
  template <class S>
  S operator()(const std::array<S, 3>& x, const std::array<std::array<S, 3>, 3>& R) const {
    using std::sqrt;
    const S r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
    const S r = sqrt(r2);
    std::array<S, 3> n;
    for (int i = 0; i < 3; i++) n[i] = (R[0][i] * x[0] + R[1][i] * x[1] + R[2][i] * x[2]) / r;
    S nJn = 0.0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) nJn += n[i] * J[i][j] * n[j];
    const double trace = J[0][0] + J[1][1] + J[2][2];
    return -mu * m / r - mu / (2.0 * r * r2) * (trace - 3.0 * nJn);
  }
};

}  // namespace

int main() {
  Suite s;
  const double pi = std::numbers::pi;
  const auto none = [](double, const Vec3&, const Rotation&) { return pm::Wrench{}; };

  // ---- mass properties ---------------------------------------------------------------------------------------
  {
    const Vec3 lo{-0.2, -0.1, 0.0}, hi{0.4, 0.3, 0.15};
    const double rho = 2700;
    const auto whole = pm::solid_box(lo, hi, rho);
    std::vector<pm::MassProperties> parts;
    const std::array<int, 3> n{4, 3, 2};
    for (int i = 0; i < n[0]; i++)
      for (int j = 0; j < n[1]; j++)
        for (int k = 0; k < n[2]; k++) {
          const auto at = [&](int a, int c) { return lo[a] + (hi[a] - lo[a]) * c / n[a]; };
          parts.push_back(pm::solid_box({at(0, i), at(1, j), at(2, k)}, {at(0, i + 1), at(1, j + 1), at(2, k + 1)}, rho));
        }
    const auto joined = pm::combine(parts);
    s.check("mass: a box cut into 24 boxes and joined by the parallel axis theorem is the box (mass, centre and inertia, relative)",
            std::max({std::fabs(joined.mass / whole.mass - 1), largest(joined.centre - whole.centre) / 0.4,
                      largest(joined.inertia - whole.inertia) / largest(whole.inertia)}),
            1e-14);
    // A dumbbell: J = μ (|d|² I − d dᵀ) about its centre, μ = m₁m₂/(m₁ + m₂), d = r₁ − r₂.
    const double m1 = 2, m2 = 3, mu = m1 * m2 / (m1 + m2);
    const Vec3 r1{1, 0, 0.5}, r2{-1, 1, 0};
    const auto bell = pm::combine({pm::point_mass(r1, m1), pm::point_mass(r2, m2)});
    const Vec3 d = r1 - r2;
    Mat3 expected{};
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) expected[i][j] = mu * ((i == j ? pg::dot(d, d) : 0) - d[i] * d[j]);
    s.check("mass: a dumbbell's inertia is μ (|d|² I − d dᵀ)", largest(bell.inertia - expected), 1e-14);
    // Principal axes of Q diag(1, 2, 3) Qᵀ.
    const Rotation Q = Rotation::about(Vec3{1, 2, -0.5}, 0.8);
    const auto P = pm::principal_axes(conjugate(Q.matrix(), {1, 2, 3}));
    const Mat3 A = P.axes.matrix(), q = Q.matrix();
    double axes = 0;
    for (int k = 0; k < 3; k++) axes = std::max(axes, std::fabs(std::fabs(A[0][k] * q[0][k] + A[1][k] * q[1][k] + A[2][k] * q[2][k]) - 1));
    s.check("mass: the principal moments of Q diag(1, 2, 3) Qᵀ", largest(P.moments - Vec3{1, 2, 3}), 1e-14);
    s.check("mass: and its principal axes are Q's columns (|cos| − 1)", axes, 1e-14);
  }

  // ---- the wrench of a potential, against its differences ------------------------------------------------------------
  {
    const MacCullagh U{1.0, 2.0, conjugate(Rotation::about(Vec3{0, 1, 1}, 0.4).matrix(), {0.1, 0.2, 0.25})};
    const Vec3 x{0.9, -0.4, 0.3};
    const Rotation R = Rotation::about(Vec3{1, -1, 2}, 1.1);
    const auto w = pm::potential_wrench(U, x, R);
    const double e = 1e-6;
    double force = 0, moment = 0;
    for (int i = 0; i < 3; i++) {
      Vec3 a = x, b = x;
      a[i] += e, b[i] -= e;
      force = std::max(force, std::fabs(w.force[i] + (pm::potential_energy(U, a, R) - pm::potential_energy(U, b, R)) / (2 * e)) / largest(w.force));
      Vec3 t{0, 0, 0};
      t[i] = e;
      const double up = pm::potential_energy(U, x, R * Rotation::exp(t)), down = pm::potential_energy(U, x, R * Rotation::exp(-t));
      moment = std::max(moment, std::fabs(w.moment[i] + (up - down) / (2 * e)) / largest(w.moment));
    }
    s.check("wrench: the force −∂U/∂x of MacCullagh's potential, by jets, against differences (relative)", force, 1e-8);
    s.check("wrench: and its moment −∂U/∂η, R(η) = R exp(η^), against differences (relative)", moment, 1e-8);
  }

  // ---- the free body ---------------------------------------------------------------------------------------------
  {
    pm::MassProperties body;
    body.mass = 1;
    body.inertia = conjugate(Rotation::about(Vec3{0.3, -1, 0.7}, 1.3).matrix(), {1, 2, 3});
    const Vec3 Pi0{1.0, 0.1, 0.5};
    const double h = 0.01;
    const pm::RigidBodyIntegrator lg(body, h);
    pm::RigidBodyState st;
    st.R = Rotation::about(Vec3{1, 1, 0}, 0.3);
    st.Pi = Pi0;
    const Vec3 L0 = st.R.apply(st.Pi);
    const double C0 = pg::norm(st.Pi), E0 = pm::kinetic_energy(body, st);
    double dL = 0, dC = 0, dE = 0;
    std::size_t iterations = 0;
    const int steps = 100000;
    for (int k = 0; k < steps; k++) {
      iterations = std::max(iterations, lg.step(st, none).iterations);
      dL = std::max(dL, pg::norm(st.R.apply(st.Pi) - L0) / pg::norm(L0));
      dC = std::max(dC, std::fabs(pg::norm(st.Pi) / C0 - 1));
      dE = std::max(dE, std::fabs(pm::kinetic_energy(body, st) / E0 - 1));
    }
    s.check("free body: the spatial angular momentum R Π is conserved to rounding (10⁵ steps, relative)", dL, 1e-11);
    s.check("free body: the Casimir |Π| to rounding (relative)", dC, 1e-11);
    s.check("free body: the attitude stays on SO(3) (|R Rᵀ − I|)", pg::orthogonality_defect(st.R.matrix()), 1e-14);
    s.check("free body: and the energy, exactly: the flow is Moser and Veselov's, and Πᵀ J_d² Π = α|Π|² + 2βE (relative)", dE, 1e-12);
    s.require("free body: Newton's method on Cayley's parameter takes four iterations at most", iterations <= 4, std::to_string(iterations) + " iterations");
    // Π(T) against Dormand–Prince on Euler's equations Π̇ = Π × J⁻¹ Π.
    const double T = 4;
    std::array<double, 3> y{Pi0[0], Pi0[1], Pi0[2]};
    const auto f = [&](const std::array<double, 3>& Pi, std::array<double, 3>& dPi) {
      const Vec3 w = pm::angular_velocity(body, {Pi[0], Pi[1], Pi[2]});
      const Vec3 c = pg::cross({Pi[0], Pi[1], Pi[2]}, w);
      dPi = {c[0], c[1], c[2]};
    };
    double hs = 1e-3;
    for (int i = 0; i < 40; i++) calx::integrators::dopri5<3>(f, y, T / 40, {1.0, 1.0, 1.0}, calx::integrators::Tolerance{1e-13, 1e-13}, hs);
    const auto error = [&](int n) {
      const pm::RigidBodyIntegrator g(body, T / n);
      pm::RigidBodyState a;
      a.Pi = Pi0;
      for (int i = 0; i < n; i++) g.step(a, none);
      return pg::norm(a.Pi - Vec3{y[0], y[1], y[2]});
    };
    const double e1 = error(100), e2 = error(200), e3 = error(400);
    s.check("free body: Π(T) converges to Euler's equations at order 2 (log₂ of the fall − 2)",
            std::max(std::fabs(std::log2(e1 / e2) - 2), std::fabs(std::log2(e2 / e3) - 2)), 0.05);
  }

  // ---- the symmetric top --------------------------------------------------------------------------------------------
  {
    const double I1 = 1.0, I3 = 2.5, T = 5;
    pm::MassProperties body;
    body.mass = 1;
    body.inertia = {{{I1, 0, 0}, {0, I1, 0}, {0, 0, I3}}};
    const Rotation R0 = Rotation::about(Vec3{0.2, 1, -0.4}, 0.9);
    const Vec3 Pi0{0.6, -0.3, 1.1};
    const Vec3 L = R0.apply(Pi0);
    const Rotation exact = Rotation::exp((T / I1) * L) * R0 * Rotation::exp(Vec3{0, 0, T * (1 / I3 - 1 / I1) * Pi0[2]});
    const auto error = [&](int n) {
      const pm::RigidBodyIntegrator g(body, T / n);
      pm::RigidBodyState a;
      a.R = R0, a.Pi = Pi0;
      for (int i = 0; i < n; i++) g.step(a, none);
      return (a.R.inverse() * exact).angle();
    };
    const double e1 = error(200), e2 = error(400), e3 = error(800);
    s.check("top: the attitude of the symmetric top converges to exp(t L^/I₁) R₀ exp(t (1/I₃ − 1/I₁) Π₃ e₃^) at order 2 (log₂ of the fall − 2)",
            std::max(std::fabs(std::log2(e1 / e2) - 2), std::fabs(std::log2(e2 / e3) - 2)), 0.05);
    s.check("top: at 800 steps it is within 10⁻⁴ rad of it", e3, 1e-4, "rad");
  }

  // ---- orbit and attitude: the gravity gradient ---------------------------------------------------------------------
  {
    // A small box on a circular orbit of the whole potential (μ = 1, r = 1), its least axis radial, its pitch
    // disturbed. The formula is the limit of a small body: a pitched body feels a tangential force
    // ∝ θ (I₂ − I₁)/r⁴, which drives the orbit's epicycle near resonance (ω_p ≈ 1.07 n), and the epicycle the
    // pitch; at a size of 3 × 10⁻² r that shifts the libration by 2 × 10⁻³, at 3 × 10⁻³ r (here) by 10⁻⁵.
    const auto body = pm::solid_box({-0.0015, -0.001, -0.0005}, {0.0015, 0.001, 0.0005}, 1 / (0.003 * 0.002 * 0.001));
    const MacCullagh U{1.0, body.mass, body.inertia};
    const auto wrench = [&](double, const Vec3& x, const Rotation& R) { return pm::potential_wrench(U, x, R); };
    const double I1 = body.inertia[0][0], I2 = body.inertia[1][1], I3 = body.inertia[2][2];
    const double expected = 2 * pi / std::sqrt(3 * (I2 - I1) / I3);  // 3 μ/r³ = 3
    const pm::RigidBodyIntegrator lg(body, 2 * pi / 4000);
    pm::RigidBodyState st;
    st.x = {1, 0, 0};
    const double v = std::sqrt(-wrench(0, st.x, Rotation{}).force[0] / body.mass);  // v²/r = the inward force over m
    st.p = {0, body.mass * v, 0};
    st.R = Rotation::about(pg::Axis::z, 0.01);
    st.Pi = body.inertia * Vec3{0, 0, v};
    std::vector<double> crossings;
    double before = 0.01, t_before = 0;
    while (crossings.size() < 21) {
      lg.step(st, wrench);
      const Vec3 radial = pg::unit(st.x), b1 = st.R.apply({1, 0, 0});
      const double pitch = std::atan2(pg::cross(radial, b1)[2], pg::dot(radial, b1));
      if ((before > 0) != (pitch > 0)) crossings.push_back(t_before + (st.t - t_before) * before / (before - pitch));
      before = pitch, t_before = st.t;
    }
    const double period = 2 * (crossings.back() - crossings.front()) / static_cast<double>(crossings.size() - 1);
    s.check("gravity gradient: a small box's pitch librates at √(3 (μ/r³)(I₂ − I₁)/I₃) (ten periods, relative)", std::fabs(period / expected - 1), 1e-4);

    // A tumbling body on an orbit of e = 0.2: orbit and attitude exchange angular momentum, the total conserved.
    const auto big = pm::solid_box({-0.15, -0.1, -0.05}, {0.15, 0.1, 0.05}, 1 / (0.3 * 0.2 * 0.1));
    const MacCullagh V{1.0, big.mass, big.inertia};
    const auto wrench2 = [&](double, const Vec3& x, const Rotation& R) { return pm::potential_wrench(V, x, R); };
    const pm::RigidBodyIntegrator g(big, 2 * pi / 2000);
    pm::RigidBodyState b;
    b.x = {0.8, 0, 0};
    b.p = {0, std::sqrt(1.2 / 0.8), 0};
    b.R = Rotation::about(Vec3{1, 2, 3}, 0.7);
    b.Pi = big.inertia * Vec3{0.7, -0.4, 1.3};
    const Vec3 L0 = pm::angular_momentum(b), orbit0 = pg::cross(b.x, b.p);
    const double E0 = pm::kinetic_energy(big, b) + pm::potential_energy(V, b.x, b.R);
    double total = 0, exchange = 0, energy = 0;
    for (int k = 0; k < 40000; k++) {
      g.step(b, wrench2);
      total = std::max(total, pg::norm(pm::angular_momentum(b) - L0) / pg::norm(L0));
      exchange = std::max(exchange, pg::norm(pg::cross(b.x, b.p) - orbit0) / pg::norm(L0));
      energy = std::max(energy, std::fabs((pm::kinetic_energy(big, b) + pm::potential_energy(V, b.x, b.R)) / E0 - 1));
    }
    s.check("gravity gradient: orbit and attitude together conserve the total angular momentum x × p + R Π to rounding (20 orbits, relative)", total, 1e-11);
    s.require("gravity gradient: while they exchange it (the orbit's changes by more than 10⁻⁶ of the total)", exchange > 1e-6,
              "the orbit's change " + std::to_string(exchange));
    s.check("gravity gradient: the energy's error stays small (relative, h = T/2000)", energy, 1e-5);
  }

  return s.done();
}
