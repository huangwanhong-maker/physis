// Checks of current fluctuations, the thermodynamic uncertainty relation and Landauer's bound
// (fluctuations.hpp; decision 0022).
//
// - A biased walk on a ring of 5 states, rates a forward and b backward, counting +1 forward: J̄ = a − b and
//   D₂ = a + b exactly (a difference of two Poisson processes), σ = (a − b) ln(a/b), and the TUR ratio
//   (a + b) ln(a/b) / (2(a − b)) ≥ 1.
// - Gillespie's simulation of a random 4-state chain: the variance of the current over long runs, divided by
//   the time, agrees with D₂ (within 4 standard errors), and its mean with J̄.
// - The relation holds for 200 random chains (the least ratio ≥ 1); near equilibrium the entropy current,
//   increments ln(k_ij / k_ji), approaches equality (within 1e−2 at a drive of 1e−3).
// - Landauer's bound at 300 K: k_B T ln 2 = 2.8710e−21 J.

#include <physis/statistical_mechanics/fluctuations.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <vector>

using namespace physis::statistical_mechanics;

int main() {
  calx::verification::Suite s;
  {
    const std::size_t n = 5;
    const double a = 2.0, b = 0.5;
    calx::linalg::Matrix k(n), d(n);
    for (std::size_t i = 0; i < n; i++) {
      const std::size_t f = (i + 1) % n, g = (i + n - 1) % n;
      k(i, f) = a, k(i, g) = b;
      d(i, f) = 1, d(i, g) = -1;
    }
    const auto c = current_statistics(k, d);
    s.check("ring: J̄ = a − b", c.mean_rate - (a - b), 1e-13);
    s.check("ring: D₂ = a + b", c.variance_rate - (a + b), 1e-12);
    s.check("ring: σ = (a − b) ln(a/b)", c.entropy_production - (a - b) * std::log(a / b), 1e-12);
    s.check("ring: the TUR ratio (a + b) ln(a/b) / (2(a − b))", c.tur_ratio - (a + b) * std::log(a / b) / (2 * (a - b)), 1e-12);
  }
  calx::random::Stream rng(2026, 11, 6, 0);
  auto random_chain = [&](std::size_t n, double bias) {
    calx::linalg::Matrix k(n), d(n);
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = i + 1; j < n; j++) {
        const double base = rng.uniform(0.2, 2.0), drive = bias * rng.uniform(-1, 1);
        k(i, j) = base * std::exp(drive / 2), k(j, i) = base * std::exp(-drive / 2);
        const double inc = rng.uniform(-1, 1);
        d(i, j) = inc, d(j, i) = -inc;
      }
    return std::pair{k, d};
  };
  {
    const auto [k, d] = random_chain(4, 2.0);
    const auto c = current_statistics(k, d);
    const double T = 200;
    const int runs = 4000;
    double sum = 0, sum2 = 0;
    std::size_t state = 0;
    for (int r = 0; r < runs; r++) {
      double t = 0, J = 0;
      while (true) {
        double total = 0;
        for (std::size_t j = 0; j < 4; j++) total += j == state ? 0 : k(state, j);
        t += -std::log(rng.uniform()) / total;
        if (t > T) break;
        double u = rng.uniform() * total;
        std::size_t next = 0;
        for (std::size_t j = 0; j < 4; j++) {
          if (j == state) continue;
          next = j;
          if ((u -= k(state, j)) <= 0) break;
        }
        J += d(state, next), state = next;
      }
      sum += J, sum2 += J * J;
    }
    const double mean = sum / runs, var = sum2 / runs - mean * mean;
    // The variance of a sample variance of near-normal data is about 2σ⁴/(runs − 1).
    const double se_var = var * std::sqrt(2.0 / (runs - 1));
    s.check("Gillespie: the variance of the current over t equals D₂ (in standard errors)", (var - c.variance_rate * T) / se_var, 4);
    s.check("Gillespie: and its mean J̄ t (in standard errors)", (mean - c.mean_rate * T) / std::sqrt(var / runs), 4);
  }
  {
    double least = INFINITY;
    for (int r = 0; r < 200; r++) {
      const auto [k, d] = random_chain(3 + r % 4, 3.0);
      least = std::min(least, current_statistics(k, d).tur_ratio);
    }
    s.require("TUR: D₂ σ / (2 J̄²) ≥ 1 for 200 random chains", least >= 1 - 1e-12, "least ratio " + std::to_string(least));
    // Equality is approached, in linear response, by the entropy current itself: increments ln(k_ij / k_ji).
    auto [k, d] = random_chain(4, 1e-3);
    for (std::size_t i = 0; i < 4; i++)
      for (std::size_t j = 0; j < 4; j++)
        if (i != j) d(i, j) = std::log(k(i, j) / k(j, i));
    s.check("TUR: near equilibrium the entropy current saturates it (ratio − 1, at a drive of 1e−3)", current_statistics(k, d).tur_ratio - 1, 1e-2);
  }
  s.check("Landauer: k_B T ln 2 at 300 K (J)", landauer_bound(300) - 2.8710e-21, 5e-25, "J");
  s.check("the precision floor: 2 k_B T / ε² at 300 K, ε = 1% (J)", precision_cost_bound(300, 0.01) / (2 * 1.380649e-23 * 300 / 1e-4) - 1, 1e-15);
  return s.done();
}
