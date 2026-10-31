// Checks of piecewise-deterministic Markov processes (piecewise_deterministic.hpp).
//
// The sawtooth with resets: on (0, 1) the state rises at unit speed, jumps to 0 at the rate λ = 2, and is
// forced back to 0 on reaching 1. Its stationary density is λ e^{−λx}/(1 − e^{−λ}) on [0, 1) (the balance
// p′ = −λp, with the flux p(1) at the boundary returning to 0), so:
// - the time average of x is 1/λ − e^{−λ}/(1 − e^{−λ}), and the forced jumps are the fraction e^{−λ} of all;
// - Dynkin's formula for f = sin 2πx, which meets the boundary condition f(1) = f(0): the time average of
//   A f = f′ − λ f vanishes;
// - for f = x², which does not (f(1) = 1 ≠ f(0) = 0), the time average of A f is the boundary flux
//   p(1) = λ e^{−λ}/(1 − e^{−λ}) — the forced jumps per unit time: the forced jump is a boundary condition
//   on the domain of A, and dropping it leaves a term.
// With no rate the resets come every unit of time exactly; with no jumps the flow x′ = −x is integrated
// to the accuracy of RK4.

#include <cmath>
#include <numbers>
#include <vector>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <physis/dynamics/piecewise_deterministic.hpp>

using calx::verification::Suite;
namespace pd = physis::dynamics;

int main() {
  Suite s;
  calx::random::Stream rng(23, 4, 5, 6);
  const double lambda = 2.0;
  pd::PiecewiseDeterministic P;
  P.flow = [](const std::vector<double>&, std::vector<double>& dx) { dx[0] = 1; };
  P.rate = [&](const std::vector<double>&) { return lambda; };
  P.jump = [](const std::vector<double>&, calx::random::Stream&) { return std::vector<double>{0.0}; };
  P.inside = [](const std::vector<double>& x) { return 1 - x[0]; };
  P.reset = [](const std::vector<double>&, calx::random::Stream&) { return std::vector<double>{0.0}; };
  const double T = 20000, h = 5e-3;
  const auto path = pd::simulate(P, {0.0}, T, h, rng, 1u << 30, 10000000);
  // Time averages, exact on each segment between jumps, where x(t) = x₀ + t.
  const double tau = 2 * std::numbers::pi;
  double mx = 0, af_in = 0, af_out = 0, t = 0, x0 = 0, drift = 0;
  const auto segment = [&](double x, double d) {
    const double y = x + d;
    mx += x * d + d * d / 2;
    af_in += (std::sin(tau * y) - std::sin(tau * x)) - lambda * (std::cos(tau * x) - std::cos(tau * y)) / tau;
    af_out += (y * y - x * x) - lambda * (y * y * y - x * x * x) / 3;
  };
  for (const auto& e : path.events) {
    segment(x0, e.time - t);
    drift = std::max(drift, std::fabs(e.before[0] - (x0 + e.time - t)));
    t = e.time, x0 = e.after[0];
  }
  segment(x0, path.end - t);
  mx /= path.end, af_in /= path.end, af_out /= path.end;
  s.check("the flow between jumps is followed exactly", drift, 1e-9);
  const double q = std::exp(-lambda), mean = 1 / lambda - q / (1 - q), flux = lambda * q / (1 - q);
  std::size_t forced = 0;
  for (const auto& e : path.events) forced += e.forced;
  s.check("the time average of x is 1/λ − e^{−λ}/(1 − e^{−λ})", mx - mean, 0.006);
  s.check("the forced jumps are the fraction e^{−λ} of all", static_cast<double>(forced) / path.events.size() - q, 0.01);
  s.check("Dynkin: f = sin 2πx meets the boundary condition and E A f = 0", af_in, 0.04);
  s.check("f = x² does not: E A f is the boundary flux λ e^{−λ}/(1 − e^{−λ})", af_out - flux, 0.01);
  s.check("which is the forced jumps per unit time", forced / path.end - flux, 0.015);

  pd::PiecewiseDeterministic D = P;
  D.rate = nullptr;
  D.jump = nullptr;
  const auto saw = pd::simulate(D, {0.0}, 10.5, 0.01, rng);
  double worst = 0;
  for (std::size_t k = 0; k < saw.events.size(); k++) worst = std::max(worst, std::fabs(saw.events[k].time - (k + 1.0)));
  s.require("with no rate, ten resets in 10.5", saw.events.size() == 10);
  s.check("each exactly one unit after the last", worst, 1e-8);

  pd::PiecewiseDeterministic E;
  E.flow = [](const std::vector<double>& x, std::vector<double>& dx) { dx[0] = -x[0]; };
  const auto decay = pd::simulate(E, {1.0}, 1.0, 0.01, rng);
  s.check("with no jumps, x′ = −x integrated to e^{−1}", decay.states.back()[0] - std::exp(-1.0), 1e-9);
  return s.done();
}
