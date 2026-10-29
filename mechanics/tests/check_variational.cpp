// Checks of discrete mechanics (mechanics/variational.hpp; decision 0013, docs/theory/discrete-mechanics.md).
//
// Each check is a closed form or an identity of the discrete theory, not a comparison with another
// integrator's numbers:
// - the trapezoidal rule's integrator of ½|v|² − V is Störmer–Verlet, to rounding (Marsden and West 2001,
//   §2.6);
// - the midpoint rule's, on a quadratic Lagrangian, is the implicit midpoint rule: the oscillator's energy
//   is conserved exactly and its phase turns by 2 atan(hω/2) each step;
// - Noether: the angular momentum of an orbit is conserved to rounding at a step so coarse that the
//   energy is wrong by a part in a thousand; the total momentum of two bodies coupled by an interaction
//   term, and the uniform motion of their centre of mass;
// - coupling in the action: two oscillators and the spring between them, whose normal modes the
//   discrete flow turns exactly;
// - a charge in a magnetic field, by the term e A·v: its discrete orbit is an exact circle turned by
//   −2 atan(hΩ/2) each step, and a change of gauge changes the momenta by e ∇χ and the motion not at all;
// - the Lagrange–d'Alembert principle: the damped oscillator at order 2;
// - a Lagrangian not quadratic in v, the free relativistic particle in a uniform field: the Legendre
//   transform and its inverse, the energy γmc² + V, and hyperbolic motion with a transverse momentum
//   at order 2;
// - the map is symplectic (the pendulum, by differences) and symmetric;
// - the Euclidean action of the oscillator is ħ times physis field_theory's free lattice action on the
//   cycle of its time steps, and the lattice's ⟨x²⟩ tends at order a² to the quantum oscillator's
//   (ħ/2mω) coth(βħω/2).

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <string>
#include <vector>

#include <calx/integrators/symplectic.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <morphe/topology/constructions.hpp>
#include <physis/field_theory/scalar.hpp>
#include <physis/mechanics/variational.hpp>

using calx::verification::Suite;
namespace pm = physis::mechanics;
namespace ci = calx::integrators;

namespace {

template <std::size_t N, class L, class F = pm::NoForce>
auto integrator(L lagrangian, double h, pm::Quadrature rule = pm::Quadrature::midpoint, F force = {}) {
  return pm::VariationalIntegrator<N, L, F>(pm::DiscreteLagrangian<N, L>(std::move(lagrangian), h, rule), std::move(force));
}

template <std::size_t N>
double distance(const pm::Coordinates<N>& a, const pm::Coordinates<N>& b) {
  double m = 0;
  for (std::size_t i = 0; i < N; i++) m = std::max(m, std::fabs(a[i] - b[i]));
  return m;
}

// The centre of the circle through three points of the plane.
pm::Coordinates<2> circumcentre(const pm::Coordinates<2>& a, const pm::Coordinates<2>& b, const pm::Coordinates<2>& c) {
  const double d = 2 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]));
  const double A = a[0] * a[0] + a[1] * a[1], B = b[0] * b[0] + b[1] * b[1], C = c[0] * c[0] + c[1] * c[1];
  return {(A * (b[1] - c[1]) + B * (c[1] - a[1]) + C * (a[1] - b[1])) / d, (A * (c[0] - b[0]) + B * (a[0] - c[0]) + C * (b[0] - a[0])) / d};
}

}  // namespace

int main() {
  Suite s;
  const double pi = std::numbers::pi;

  // ---- the trapezoidal rule is Störmer–Verlet -------------------------------------------------------------------
  {
    const auto L = [](double, const auto& q, const auto& v) {
      using std::sqrt;
      return 0.5 * (v[0] * v[0] + v[1] * v[1]) + 1.0 / sqrt(q[0] * q[0] + q[1] * q[1]);
    };
    const double e = 0.6, h = 2 * pi / 400;
    auto vi = integrator<2>(L, h, pm::Quadrature::trapezoid);
    auto st = vi.start(0, {1 - e, 0}, {0, std::sqrt((1 + e) / (1 - e))});
    ci::Vec<2> q{1 - e, 0}, p{0, std::sqrt((1 + e) / (1 - e))};
    const ci::Vec<2> inv{1, 1};
    const auto grad = [](const ci::Vec<2>& x, ci::Vec<2>& g) {
      const double r3 = std::pow(x[0] * x[0] + x[1] * x[1], 1.5);
      g = {x[0] / r3, x[1] / r3};
    };
    double worst = 0;
    std::size_t iterations = 0;
    for (int k = 0; k < 4000; k++) {
      iterations = std::max(iterations, vi.step(st).iterations);
      ci::stormer_verlet<2>(inv, grad, q, p, h);
      worst = std::max({worst, distance<2>(st.q, q), distance<2>(st.p, p)});
    }
    s.check("trapezoid: the integrator of ½|v|² − V is Störmer–Verlet (ten orbits of e = 0.6; max |Δq|, |Δp|)", worst, 1e-10);
    s.require("trapezoid: its step's equation is linear in q₁, and Newton's method meets it in two iterations at most", iterations <= 2,
              std::to_string(iterations) + " iterations");
  }

  // ---- the midpoint rule on the oscillator ----------------------------------------------------------------------
  {
    const double m = 2.0, k = 8.0, w = std::sqrt(k / m), h = 0.05;
    const auto L = [=](double, const auto& q, const auto& v) { return 0.5 * m * v[0] * v[0] - 0.5 * k * q[0] * q[0]; };
    auto vi = integrator<1>(L, h);
    auto st = vi.start(0, {1.0}, {0.0});
    const auto H = [&] { return st.p[0] * st.p[0] / (2 * m) + 0.5 * k * st.q[0] * st.q[0]; };
    const double H0 = H();
    double worst = 0;
    const int steps = 100000;
    for (int i = 0; i < steps; i++) vi.step(st), worst = std::max(worst, std::fabs(H() - H0) / H0);
    s.check("midpoint: on the oscillator it is the implicit midpoint rule, and conserves p²/2m + kq²/2 exactly (10⁵ steps, relative)", worst, 1e-12);
    const double theta = 2 * std::atan(h * w / 2);
    const double phase = std::atan2(-st.p[0] / (m * w), st.q[0]);
    s.check("midpoint: its phase turns by 2 atan(hω/2) each step (after 10⁵ steps, rad)", std::fabs(std::remainder(phase - steps * theta, 2 * pi)), 1e-8);
  }

  // ---- Noether: an orbit's angular momentum, whatever the step --------------------------------------------------
  {
    const auto L = [](double, const auto& q, const auto& v) {
      using std::sqrt;
      return 0.5 * (v[0] * v[0] + v[1] * v[1]) + 1.0 / sqrt(q[0] * q[0] + q[1] * q[1]);
    };
    const double e = 0.5, h = 2 * pi / 100;  // a hundred steps an orbit
    auto vi = integrator<2>(L, h);
    auto st = vi.start(0, {1 - e, 0}, {0, std::sqrt((1 + e) / (1 - e))});
    const auto J = [&] { return pm::momentum_map<2>(st.p, {-st.q[1], st.q[0]}); };  // rotations: ξ_Q(q) = (−y, x)
    const auto E = [&] { return 0.5 * (st.p[0] * st.p[0] + st.p[1] * st.p[1]) - 1 / std::hypot(st.q[0], st.q[1]); };
    const double J0 = J(), E0 = E();
    double dJ = 0, first = 0, last = 0;
    for (int orbit = 0; orbit < 1000; orbit++)
      for (int i = 0; i < 100; i++) {
        vi.step(st);
        dJ = std::max(dJ, std::fabs(J() / J0 - 1));
        const double dE = std::fabs(E() / E0 - 1);
        if (orbit < 100) first = std::max(first, dE);
        if (orbit >= 900) last = std::max(last, dE);
      }
    s.check("noether: the angular momentum ⟨p, ξ_Q(q)⟩ of an orbit of e = 0.5 is conserved to rounding (10⁵ steps of T/100, relative)", dJ, 1e-11);
    s.check("noether: at that step the energy's error is not small (relative, at least 10⁻³), but bounded: the last hundred orbits' over the first's − 1",
            first >= 1e-3 ? std::fabs(last / first - 1) : 1.0, 0.1);
  }

  // ---- coupling in the action: two bodies and their interaction --------------------------------------------------
  {
    // Two bodies on a line, and an anharmonic interaction between them: L = L_a + L_b + L_ab.
    const double ma = 1.0, mb = 2.5, kappa = 3.0, mu = 0.7;
    const auto La = [=](double, const auto&, const auto& v) { return 0.5 * ma * v[0] * v[0]; };
    const auto Lb = [=](double, const auto&, const auto& v) { return 0.5 * mb * v[1] * v[1]; };
    const auto Lab = [=](double, const auto& q, const auto&) {
      const auto d = q[0] - q[1];
      return -(0.5 * kappa * d * d + 0.25 * mu * d * d * d * d);
    };
    const auto L = [=](double t, const auto& q, const auto& v) { return La(t, q, v) + Lb(t, q, v) + Lab(t, q, v); };
    auto vi = integrator<2>(L, 0.05);
    auto st = vi.start(0, {0.0, 1.3}, {0.8, -0.1});
    // Exact in exact arithmetic. In floating point the discrete velocity (q₁ − q₀)/h carries the rounding of
    // the positions over h, so the momentum changes each step by a few ulps of m|q|/h (the centre of mass
    // drifts to |q| ≈ 160 m here), with no part proportional to a power of h.
    const double P0 = st.p[0] + st.p[1], X0 = (ma * st.q[0] + mb * st.q[1]) / (ma + mb);
    double ulps = 0, drift = 0, dX = 0, P = P0;
    const int steps = 20000;
    for (int k = 1; k <= steps; k++) {
      vi.step(st);
      const double scale = std::fabs(st.p[0]) + std::fabs(st.p[1]);
      const double rounding = (ma + mb) * std::max(std::fabs(st.q[0]), std::fabs(st.q[1])) / 0.05 * 0x1p-52;
      ulps = std::max(ulps, std::fabs(st.p[0] + st.p[1] - P) / rounding);
      P = st.p[0] + st.p[1];
      drift = std::max(drift, std::fabs(P - P0) / scale);
      const double X = (ma * st.q[0] + mb * st.q[1]) / (ma + mb), travelled = k * 0.05 * P0 / (ma + mb);
      dX = std::max(dX, std::fabs(X - X0 - travelled) / std::fabs(travelled));
    }
    s.check("coupling: two bodies and the interaction between them conserve the total momentum to rounding (its change a step, in ulps of (m_a + m_b)|q|/h)", ulps, 4);
    s.check("coupling: over 20 000 steps it drifts by rounding alone (relative to |p_a| + |p_b|)", drift, 1e-10);
    s.check("coupling: and their centre of mass moves uniformly, X_k = X₀ + k h P/M (relative to the distance travelled)", dX, 1e-10);

    // Two oscillators and the spring between them: the normal modes (q₁ ± q₂)/2 turn independently, each by
    // 2 atan(hω/2), ω² = k/m and (k + 2κ)/m. Started from q = (1, 0) at rest, the energy beats between them.
    const double m = 1.0, k = 1.0, c = 0.05, h = 0.1;
    const auto Lo = [=](double, const auto& q, const auto& v) {
      const auto d = q[0] - q[1];
      return 0.5 * m * (v[0] * v[0] + v[1] * v[1]) - 0.5 * k * (q[0] * q[0] + q[1] * q[1]) - 0.5 * c * d * d;
    };
    auto vo = integrator<2>(Lo, h);
    auto so = vo.start(0, {1.0, 0.0}, {0.0, 0.0});
    const double ts = 2 * std::atan(h * std::sqrt(k / m) / 2), ta = 2 * std::atan(h * std::sqrt((k + 2 * c) / m) / 2);
    double worst = 0, lowest = 1;
    for (int n = 1; n <= 10000; n++) {
      vo.step(so);
      const double sym = 0.5 * std::cos(n * ts), anti = 0.5 * std::cos(n * ta);
      worst = std::max({worst, std::fabs(so.q[0] - (sym + anti)), std::fabs(so.q[1] - (sym - anti))});
      lowest = std::min(lowest, (so.p[0] * so.p[0] / m + k * so.q[0] * so.q[0]) / k);  // the first oscillator's energy, over its start
    }
    s.check("coupling: two oscillators and their spring follow their discrete normal modes exactly (10⁴ steps, m)", worst, 1e-10, "m");
    s.check("coupling: the first gives its energy to the second at the beat (its least energy, over its first)", lowest, 1e-3);
  }

  // ---- a charge in a magnetic field: the coupling e A·v, and its gauge -------------------------------------------
  {
    const double m = 1.2, e = 0.8, B = 1.5, Omega = e * B / m, h = 0.1;
    // The symmetric gauge A = ½ B (−y, x) and Landau's A = B (−y, 0) differ by ∇χ, χ = −½ B x y.
    const auto Ls = [=](double, const auto& q, const auto& v) { return 0.5 * m * (v[0] * v[0] + v[1] * v[1]) + e * 0.5 * B * (-q[1] * v[0] + q[0] * v[1]); };
    const auto Ll = [=](double, const auto& q, const auto& v) { return 0.5 * m * (v[0] * v[0] + v[1] * v[1]) - e * B * q[1] * v[0]; };
    auto vs = integrator<2>(Ls, h);
    auto vl = integrator<2>(Ll, h);
    auto ss = vs.start(0, {0.3, -0.2}, {0.9, 0.4});
    auto sl = vl.start(0, {0.3, -0.2}, {0.9, 0.4});
    std::vector<pm::Coordinates<2>> path{ss.q};
    double same = 0, momenta = 0;
    for (int k = 0; k < 2000; k++) {
      vs.step(ss), vl.step(sl);
      path.push_back(ss.q);
      same = std::max(same, distance<2>(ss.q, sl.q));
      const pm::Coordinates<2> shift{ss.p[0] + e * (-0.5 * B * ss.q[1]), ss.p[1] + e * (-0.5 * B * ss.q[0])};
      momenta = std::max(momenta, distance<2>(shift, sl.p));
    }
    const auto c = circumcentre(path[0], path[1], path[2]);
    const double rho = std::hypot(path[0][0] - c[0], path[0][1] - c[1]);
    double radius = 0, angle = 0;
    const double theta = -2 * std::atan(h * Omega / 2);  // a positive charge in B ẑ turns clockwise
    for (std::size_t k = 0; k + 1 < path.size(); k++) {
      const double ax = path[k][0] - c[0], ay = path[k][1] - c[1], bx = path[k + 1][0] - c[0], by = path[k + 1][1] - c[1];
      radius = std::max(radius, std::fabs(std::hypot(bx, by) / rho - 1));
      angle = std::max(angle, std::fabs(std::atan2(ax * by - ay * bx, ax * bx + ay * by) - theta));
    }
    s.check("charge: its discrete orbit in a uniform B is an exact circle (2000 steps, relative radius)", radius, 1e-12);
    s.check("charge: turned by −2 atan(hΩ/2) each step, Ω = eB/m (rad)", angle, 1e-12);
    s.check("gauge: the motion in the symmetric and in Landau's gauge is the same (m)", same, 1e-12, "m");
    s.check("gauge: their canonical momenta differ by e ∇χ, χ = −½ B x y", momenta, 1e-12);
  }

  // ---- Lagrange–d'Alembert: the damped oscillator ----------------------------------------------------------------
  {
    const double m = 1.0, k = 4.0, c = 0.4, g = c / (2 * m), wd = std::sqrt(k / m - g * g), T = 5;
    const auto L = [=](double, const auto& q, const auto& v) { return 0.5 * m * v[0] * v[0] - 0.5 * k * q[0] * q[0]; };
    const auto F = [=](double, const auto&, const auto& v) { return std::array{-c * v[0]}; };
    const double exact = std::exp(-g * T) * (std::cos(wd * T) + g / wd * std::sin(wd * T));
    const auto error = [&](int n) {
      auto vi = integrator<1>(L, T / n, pm::Quadrature::midpoint, F);
      auto st = vi.start(0, {1.0}, {0.0});
      for (int i = 0; i < n; i++) vi.step(st);
      return std::fabs(st.q[0] - exact);
    };
    const double e1 = error(100), e2 = error(200), e3 = error(400);
    s.check("forces: the damped oscillator by the discrete Lagrange–d'Alembert principle, order 2 (log₂ of the fall − 2)",
            std::max(std::fabs(std::log2(e1 / e2) - 2), std::fabs(std::log2(e2 / e3) - 2)), 0.05);
  }

  // ---- a Lagrangian not quadratic in v: the relativistic particle in a uniform field -------------------------------
  {
    // c = 1, m = 1; the field's force m a along x: L = −m√(1 − |v|²) + m a x.
    const double mass = 1.0, a = 1.0;
    const auto L = [=](double, const auto& q, const auto& v) {
      using std::sqrt;
      return -mass * sqrt(1 - (v[0] * v[0] + v[1] * v[1])) + mass * a * q[0];
    };
    const pm::Coordinates<2> q0{0, 0}, v0{0.6, 0.3};
    const double gamma = 1 / std::sqrt(1 - 0.45);
    const auto p = pm::momentum<2>(L, 0, q0, v0);
    s.check("legendre: p = ∂L/∂v = γ m v for the relativistic particle", std::max(std::fabs(p[0] - gamma * mass * 0.6), std::fabs(p[1] - gamma * mass * 0.3)), 1e-15);
    s.check("legendre: its inverse returns the velocity, by Newton's method on ∂²L/∂v²", distance<2>(pm::velocity<2>(L, 0, q0, p), v0), 1e-14);
    s.check("legendre: the energy p·v − L is γ m c² − m a x", std::fabs(pm::energy<2>(L, 0, {0.4, 0}, v0) - (gamma * mass - mass * a * 0.4)), 1e-14);
    // Hyperbolic motion with a transverse momentum P: p_x = m a t, p_y = P, ε = √(m² + P²),
    // x = (√(ε² + (a t)²) − ε)/a, y = (P/a) asinh(a t/ε) (m = 1), from the instant t₀ at which p_x = 0.
    const double P = p[1], t0 = p[0] / (mass * a);
    const double eps = std::sqrt(mass * mass + P * P), T = 2;
    const auto X = [&](double t) { return (std::sqrt(eps * eps + (a * t) * (a * t)) - eps) / a; };
    const auto Y = [&](double t) { return P / a * std::asinh(a * t / eps); };
    const auto error = [&](int n) {
      auto vi = integrator<2>(L, T / n);
      auto st = vi.start(0, {0, 0}, v0);
      for (int i = 0; i < n; i++) vi.step(st);
      return std::hypot(st.q[0] - (X(t0 + T) - X(t0)), st.q[1] - (Y(t0 + T) - Y(t0)));
    };
    const double e1 = error(100), e2 = error(200), e3 = error(400);
    s.check("relativity: hyperbolic motion with a transverse momentum, at order 2 (log₂ of the fall − 2)",
            std::max(std::fabs(std::log2(e1 / e2) - 2), std::fabs(std::log2(e2 / e3) - 2)), 0.05);
  }

  // ---- the map is symplectic and symmetric: the pendulum ------------------------------------------------------------
  {
    const auto L = [](double, const auto& q, const auto& v) {
      using std::cos;
      return 0.5 * v[0] * v[0] + cos(q[0]);
    };
    const auto step = [&](double q, double p, double h) {
      auto vi = integrator<1>(L, h);
      pm::PhaseState<1> st{0, {q}, {p}};
      vi.step(st);
      return std::array<double, 2>{st.q[0], st.p[0]};
    };
    const double d = 1e-6, h = 0.3;
    const auto a = step(1.0 + d, 0.4, h), b = step(1.0 - d, 0.4, h), c = step(1.0, 0.4 + d, h), e = step(1.0, 0.4 - d, h);
    const double J = ((a[0] - b[0]) * (c[1] - e[1]) - (a[1] - b[1]) * (c[0] - e[0])) / (4 * d * d);
    s.check("symplectic: one midpoint step of the pendulum keeps area (the Jacobian's determinant − 1, by differences)", std::fabs(J - 1), 1e-8);
    const auto f = step(1.0, 0.4, h);
    const auto r = step(f[0], -f[1], h);
    s.check("symmetric: a step and its reverse return to the start", std::max(std::fabs(r[0] - 1.0), std::fabs(-r[1] - 0.4)), 1e-13);
  }

  // ---- the Euclidean action and the quantum oscillator -------------------------------------------------------------
  {
    const double m = 1.3, w = 0.7, hbar = 1.05;
    const auto LE = [=](const auto& x, const auto& v) { return 0.5 * m * v[0] * v[0] + 0.5 * m * w * w * x[0] * x[0]; };
    // A closed path of n points, time step a, against the lattice action of φ = x / √(aħ/m) at m² = (aω)².
    const std::size_t n = 24;
    const double a = 0.15, unit = std::sqrt(a * hbar / m);
    calx::random::Stream z(20261029, 0x70617468u, 0, 0);  // "path"
    std::vector<pm::Coordinates<1>> path(n);
    std::vector<double> phi(n);
    for (std::size_t k = 0; k < n; k++) path[k] = {z.normal()}, phi[k] = path[k][0] / unit;
    const auto X = morphe::topology::torus_lattice({n});
    const double SE = pm::euclidean_action<1>(LE, a, path);
    const double S = physis::field_theory::action({(a * w) * (a * w), 0.0}, X, phi);
    s.check("euclidean: S_E/ħ of a closed path is field_theory's free lattice action on the cycle of its steps, m² = (aω)² (relative)", std::fabs(SE / hbar / S - 1), 1e-13);
    // ⟨x²⟩ = (aħ/m) G(0) on the cycle of n steps, β ħ = n a, against (ħ/2mω) coth(βħω/2).
    const double beta_hbar = 4.0, quantum = hbar / (2 * m * w) / std::tanh(beta_hbar * w / 2);
    const auto error = [&](std::size_t steps) {
      const double as = beta_hbar / static_cast<double>(steps);
      return std::fabs(as * hbar / m * physis::field_theory::lattice_propagator({steps}, (as * w) * (as * w), {0}) - quantum);
    };
    const double e1 = error(64), e2 = error(128), e3 = error(256);
    s.check("path integral: the lattice's ⟨x²⟩ tends to the quantum oscillator's (ħ/2mω) coth(βħω/2) at order a² (log₂ of the fall − 2)",
            std::max(std::fabs(std::log2(e1 / e2) - 2), std::fabs(std::log2(e2 / e3) - 2)), 0.02);
    s.check("path integral: at 1024 steps it is within 10⁻⁵ of it (relative)", error(1024) / quantum, 1e-5);
  }

  return s.done();
}
