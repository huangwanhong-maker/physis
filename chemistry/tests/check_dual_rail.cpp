// Checks of dual-rail compilation (dual_rail.hpp; decision 0022).
//
// - Exactness: for Lorenz's system the network's mass action gives d(pᵢ − qᵢ)/dt = Pᵢ(p − q) at random positive
//   states (to rounding), with every rail's concentration non-negative.
// - The network carries Lorenz's trajectory: integrated from p = x₀⁺, q = x₀⁻, its differences follow Lorenz's to
//   t = 1 (within 1e−8 relative to the attractor's size).
// - A harmonic oscillator x' = y, y' = −x compiles to 4 species and turns at angular frequency 1 (to 1e−8 after
//   one turn).
// - Every reaction's rate is positive and its catalysts are unchanged.

#include <physis/chemistry/dual_rail.hpp>
#include <physis/dynamics/flows.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <numbers>

using namespace physis::chemistry;
using morphe::arithmetic::Integer;
using morphe::arithmetic::Rational;

int main() {
  calx::verification::Suite s;
  const auto P = lorenz_polynomials(10, 28, Rational(Integer(8), Integer(3)));
  const auto D = dual_rail(P);
  s.require("Lorenz: six rails", D.network.species.size() == 6);
  bool positive = true, catalysts = true;
  for (const auto& r : D.network.reactions) {
    positive = positive && Rational(0) < r.rate;
    // Each reaction changes exactly one species by +1, or (annihilation) two by −1.
    long long changed = 0;
    for (std::size_t i = 0; i < D.network.species.size(); i++) changed += std::llabs(D.network.change(i, &r - D.network.reactions.data()));
    catalysts = catalysts && (changed == 1 || r.name.rfind("annihilation", 0) == 0);
  }
  s.require("every rate positive", positive);
  s.require("each reaction feeds one rail, its other species catalysts (or annihilates a pair)", catalysts);

  calx::random::Stream rng(9, 0, 0, 0);
  double worst = 0;
  for (int t = 0; t < 50; t++) {
    std::vector<double> c(6), dc(6), x(3);
    for (double& v : c) v = rng.uniform(0, 20);
    for (int i = 0; i < 3; i++) x[i] = c[D.plus[i]] - c[D.minus[i]];
    mass_action(D.network, c, dc);
    for (int i = 0; i < 3; i++) {
      const double target = P.evaluate(i, x);
      worst = std::max(worst, std::fabs(dc[D.plus[i]] - dc[D.minus[i]] - target) / (1 + std::fabs(target)));
    }
  }
  s.check("exactness: d(p − q)/dt = P(p − q) at random states (relative)", worst, 1e-12);

  {
    physis::dynamics::State c{1, 0, 1, 0, 20, 0};  // x = 1, y = 1, z = 20
    physis::dynamics::State L{1, 1, 20};
    auto crn = [&](double, const physis::dynamics::State& y, physis::dynamics::State& dy) { mass_action(D.network, y, dy); };
    double least = 1e300;
    physis::dynamics::integrate(crn, 0.0, 1.0, c, 1e-4, [&](double, const physis::dynamics::State& y) { for (double v : y) least = std::min(least, v); });
    physis::dynamics::integrate(physis::dynamics::Lorenz{}, 0.0, 1.0, L, 1e-4);
    double err = 0;
    for (int i = 0; i < 3; i++) err = std::max(err, std::fabs(c[D.plus[i]] - c[D.minus[i]] - L[i]));
    s.check("the network carries Lorenz's trajectory to t = 1 (relative to the attractor's size 20)", err / 20, 1e-8);
    s.require("its concentrations stay non-negative", least >= 0);
  }
  {
    PolynomialSystem H;
    H.variables = {"x", "y"};
    H.rates = {{{1, {0, 1}}}, {{-1, {1, 0}}}};
    const auto R = dual_rail(H, 5);
    s.require("oscillator: four rails", R.network.species.size() == 4);
    physis::dynamics::State c{1, 0, 0, 0};
    auto crn = [&](double, const physis::dynamics::State& y, physis::dynamics::State& dy) { mass_action(R.network, y, dy); };
    physis::dynamics::integrate(crn, 0.0, 2 * std::numbers::pi, c, 1e-4);
    s.check("oscillator: after one turn x = 1, y = 0", std::hypot(c[0] - c[1] - 1, c[2] - c[3]), 1e-8);
  }
  return s.done();
}
