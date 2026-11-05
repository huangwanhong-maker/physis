// Checks of continuation and of Kuramoto's networks (continuation.hpp, kuramoto.hpp; decision 0022).
//
// - The fold of x' = μ − x²: followed from μ = 1 down the branch x = √μ, the curve turns at μ = 0 (to 1e−8), and
//   comes back up the branch x = −√μ, which is unstable.
// - x' = μ + x − x³ has folds at μ = ±2/(3√3), x = ∓1/√3 (to 1e−8): the hysteresis loop of a bistable switch.
// - The Brusselator, x' = a − (b + 1)x + x²y, y' = bx − x²y, loses its equilibrium (a, b/a) at the Hopf point
//   b = 1 + a² with frequency a (to 1e−7), and has no fold.
// - Kuramoto: for Lorentzian frequencies of half-width γ, the order parameter settles at √(1 − 2γ/K) at K = 2K_c
//   (within 0.03 for N = 2000) and stays incoherent at K = K_c/2 (r below 0.1).

#include <physis/dynamics/continuation.hpp>
#include <physis/dynamics/kuramoto.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <complex>

using namespace physis::dynamics;

int main() {
  calx::verification::Suite s;
  {
    auto f = [](const State& x, double mu, State& out) { out[0] = mu - x[0] * x[0]; };
    ContinuationOptions o;
    o.direction = -1, o.step = 0.02, o.mu_max = 1.2, o.max_points = 400;
    const auto B = continue_equilibria(f, State{1.0}, 1.0, o);
    bool found = false;
    double mu = 1;
    for (const auto& b : B.bifurcations)
      if (b.kind == Bifurcation::Kind::fold) found = true, mu = b.mu;
    s.require("fold: x' = μ − x² turns at a fold", found);
    s.check("fold: at μ = 0", mu, 1e-8);
    bool unstable_lower = false;
    for (const auto& p : B.points)
      if (p.x[0] < -0.1) unstable_lower = !p.stable;
    s.require("fold: the lower branch x = −√μ is unstable", unstable_lower);
  }
  {
    auto f = [](const State& x, double mu, State& out) { out[0] = mu + x[0] - x[0] * x[0] * x[0]; };
    ContinuationOptions o;
    o.direction = +1, o.step = 0.02, o.mu_min = -1.0, o.mu_max = 1.0, o.max_points = 2000;
    const auto B = continue_equilibria(f, State{-1.5}, -1.0, o);
    std::vector<double> folds;
    for (const auto& b : B.bifurcations)
      if (b.kind == Bifurcation::Kind::fold) folds.push_back(b.mu);
    s.require("cusp: two folds on the S-shaped branch", folds.size() == 2);
    const double m = 2 / (3 * std::sqrt(3.0));
    if (folds.size() == 2) {
      s.check("cusp: the first fold at μ = +2/(3√3)", folds[0] - m, 1e-8);
      s.check("cusp: the second at μ = −2/(3√3)", folds[1] + m, 1e-8);
    }
  }
  {
    const double a = 1.3;
    auto f = [&](const State& x, double b, State& out) {
      out[0] = a - (b + 1) * x[0] + x[0] * x[0] * x[1];
      out[1] = b * x[0] - x[0] * x[0] * x[1];
    };
    ContinuationOptions o;
    o.direction = +1, o.step = 0.01, o.mu_max = 4.0, o.max_points = 2000;
    const auto B = continue_equilibria(f, State{a, 1.0 / a}, 1.0, o);
    const Bifurcation* h = nullptr;
    std::size_t folds = 0;
    for (const auto& b : B.bifurcations) {
      if (b.kind == Bifurcation::Kind::hopf && !h) h = &b;
      if (b.kind == Bifurcation::Kind::fold) folds++;
    }
    s.require("Brusselator: a Hopf point and no fold", h != nullptr && folds == 0);
    if (h) {
      s.check("Brusselator: the Hopf point b = 1 + a²", h->mu - (1 + a * a), 1e-7);
      s.check("Brusselator: its frequency a", h->frequency - a, 1e-7);
    }
  }
  {
    const std::size_t N = 2000;
    const double gamma = 0.5;
    Kuramoto k{lorentzian_frequencies(N, 0.0, gamma), 0.0, {}};
    State theta(N);
    for (std::size_t i = 0; i < N; i++) theta[i] = 2 * std::numbers::pi * std::fmod(0.61803398875 * i, 1.0);
    auto settle = [&](double K) {
      k.coupling = K;
      State th = theta;
      integrate(k, 0.0, 200.0, th, 0.01);
      double r = 0;
      const int m = 200;
      for (int q = 0; q < m; q++) {
        integrate(k, 0.0, 0.25, th, 0.01);
        r += std::abs(order_parameter(th));
      }
      return r / m;
    };
    const double Kc = kuramoto_critical_coupling(1 / (std::numbers::pi * gamma));
    s.check("Kuramoto: K_c = 2γ for Lorentzian frequencies", Kc - 2 * gamma, 1e-15);
    s.check("Kuramoto: r = √(1 − 2γ/K) at K = 2K_c", settle(2 * Kc) - kuramoto_lorentzian_order(2 * Kc, gamma), 0.03);
    s.require("Kuramoto: incoherent at K = K_c/2 (r < 0.1)", settle(0.5 * Kc) < 0.1);
  }
  return s.done();
}
