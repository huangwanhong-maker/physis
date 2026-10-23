// Checks of physis/quantum: the exponential of a 2x2 complex matrix against calx's on its real
// form; the propagator between detections against its closed form on resonance; the
// stationary state of the master equation against its exact flow run long, and, without
// dephasing, against the textbook population; the Bloch ball kept; and the bond between the
// two descriptions: the mean waiting time of the unravelling equals the inverse of the master
// equation's emission rate.

#include <physis/quantum/two_level.hpp>
#include <calx/linalg/expm.hpp>
#include <calx/quadrature/gauss_kronrod.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <string>

namespace pq = physis::quantum;
namespace pd = physis::dynamics;
using pq::Complex;

int main() {
  calx::verification::Suite s;

  // ---- exp2 against expm of the real form [[X, -Y], [Y, X]] of M = X + iY ------------------------------------------
  {
    double worst = 0;
    calx::random::Stream z(7, 0, 0, 0);
    auto check = [&](const pq::Matrix2& M) {
      calx::linalg::Matrix R(4);
      for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++) {
          R(i, j) = R(i + 2, j + 2) = M[i][j].real();
          R(i + 2, j) = M[i][j].imag();
          R(i, j + 2) = -M[i][j].imag();
        }
      const auto E = calx::linalg::expm(R);
      const auto e = pq::exp2(M);
      double size = 0, diff = 0;
      for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++) {
          size = std::max(size, std::abs(e[i][j]));
          diff = std::max(diff, std::abs(e[i][j] - Complex(E(i, j), E(i + 2, j))));
        }
      worst = std::max(worst, diff / size);
    };
    for (int k = 0; k < 200; k++) {
      const double scale = std::pow(10.0, -3 + 4 * z.uniform());
      pq::Matrix2 M;
      for (auto& row : M)
        for (auto& x : row) x = scale * Complex(2 * z.uniform() - 1, 2 * z.uniform() - 1);
      check(M);
    }
    check({{{Complex(0.3, 0.1), Complex(1, 0)}, {Complex(0, 0), Complex(0.3, 0.1)}}});  // a Jordan block: q = 0
    s.check("quantum: exp2 against calx's expm of the real form, 200 random matrices and a Jordan block (relative)", worst, 1e-13);
  }
  // ---- the propagator between detections on resonance, against its closed form ----------------------------------------
  {
    // From |g>, with Omega' = sqrt(Omega^2 - gamma^2 / 4):
    //   c_g = e^{-gamma t / 4} (cos(Omega' t / 2) + gamma / (2 Omega') sin(Omega' t / 2)),
    //   c_e = -i (Omega / Omega') e^{-gamma t / 4} sin(Omega' t / 2).
    const pq::TwoLevel q{400, 0, 100, 0};
    const double W = std::sqrt(q.Omega * q.Omega - q.gamma * q.gamma / 4);
    double worst = 0;
    for (int k = 1; k <= 100; k++) {
      const double t = k * 1e-3;
      const auto c = pq::apply(pq::no_jump_propagator(q, t), {1.0, 0.0});
      const double e = std::exp(-q.gamma * t / 4);
      const Complex cg(e * (std::cos(W * t / 2) + q.gamma / (2 * W) * std::sin(W * t / 2)), 0);
      const Complex ce(0, -q.Omega / W * e * std::sin(W * t / 2));
      worst = std::max({worst, std::abs(c[0] - cg), std::abs(c[1] - ce)});
    }
    s.check("quantum: the propagator between detections on resonance against its closed form", worst, 1e-13);
  }
  // ---- the stationary state of the master equation ---------------------------------------------------------------------
  {
    const pq::TwoLevel q{300, 120, 100, 40};
    const auto P = pd::linear_propagator<pq::Bloch>(q, 1e-3);
    pd::State<3> r{0, 0, -1};
    for (int k = 0; k < 20000; k++) P.apply(r, 0.0);
    const auto st = pq::stationary(q);
    s.check("quantum: the stationary Bloch vector against the master equation's exact flow after 20 s",
            std::max({std::fabs(r[0] - st[0]), std::fabs(r[1] - st[1]), std::fabs(r[2] - st[2])}), 1e-12);
    const pq::TwoLevel c{300, 120, 100, 0};
    const double textbook = (c.Omega * c.Omega / 4) / (c.Delta * c.Delta + c.gamma * c.gamma / 4 + c.Omega * c.Omega / 2);
    s.check("quantum: without dephasing, rho_ee = (Omega^2 / 4) / (Delta^2 + gamma^2 / 4 + Omega^2 / 2) (relative)",
            pq::excited(pq::stationary(c)) / textbook - 1, 1e-14);
    // The flow keeps the Bloch ball: |r| <= 1 from every pure state.
    calx::random::Stream z(11, 0, 0, 0);
    double most = 0;
    for (int k = 0; k < 50; k++) {
      const double a = std::acos(2 * z.uniform() - 1), b = 2 * std::numbers::pi * z.uniform();
      pd::State<3> x{std::sin(a) * std::cos(b), std::sin(a) * std::sin(b), std::cos(a)};
      for (int j = 0; j < 400; j++) {
        P.apply(x, 0.0);
        most = std::max(most, std::hypot(x[0], x[1], x[2]));
      }
    }
    s.require("quantum: the master equation keeps every state in the Bloch ball", most <= 1 + 1e-12, "largest |r| " + calx::verification::shortest(most));
  }
  // ---- the unravelling and the master equation: the mean waiting time is the inverse rate ------------------------------
  {
    double worst = 0;
    for (const pq::TwoLevel& q : {pq::TwoLevel{400, 0, 100, 0}, pq::TwoLevel{50, 0, 100, 0}, pq::TwoLevel{300, 200, 100, 0}, pq::TwoLevel{1000, -50, 30, 0}}) {
      const double R = pq::emission_rate(q);
      auto survival = [&](double t) {
        const auto c = pq::apply(pq::no_jump_propagator(q, t), {1.0, 0.0});
        return std::norm(c[0]) + std::norm(c[1]);
      };
      const auto mean = calx::quadrature::integrate(survival, 0, 80 / R, 1e-15, 1e-13, 4000);
      worst = std::max(worst, std::fabs(mean.value * R - 1));
    }
    s.check("quantum: the mean waiting time of the unravelling times the master equation's emission rate, less 1, for four drives", worst, 1e-11);
  }
  return s.done();
}
