// Networks of phase oscillators: Kuramoto's model and its order parameter (decision 0022).
//
// N oscillators with natural frequencies ωᵢ, coupled through their phase differences (Kuramoto 1975; 1984,
// Chemical Oscillations, Waves, and Turbulence, ch. 5):
//
//   θᵢ' = ωᵢ + (K/N) Σⱼ sin(θⱼ − θᵢ)              (all to all), or
//   θᵢ' = ωᵢ + K Σⱼ Aᵢⱼ sin(θⱼ − θᵢ)               (on a graph of weights A).
//
// With the order parameter r e^{iψ} = (1/N) Σ e^{iθⱼ}, the all-to-all sum is K r sin(ψ − θᵢ), O(N) per step.
// For a unimodal symmetric law g of the frequencies, partial synchrony appears at K_c = 2/(π g(0)); for the
// Lorentzian of half-width γ, K_c = 2γ and the stationary order parameter is r = √(1 − 2γ/K) above it, exactly
// in the limit N → ∞ (Kuramoto 1984; on its invariant manifold, Ott and Antonsen 2008, Chaos 18, 037113).
// lorentzian_frequencies draws them as deterministic quantiles γ tan(π(i + ½)/N − π/2), so that a finite
// network's law is as close to the Lorentzian as N allows. Checks: tests/check_bifurcation.cpp.
#pragma once

#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <vector>

#include <calx/linalg/dense.hpp>

#include "flows.hpp"

namespace physis::dynamics {

inline std::complex<double> order_parameter(const State& theta) {
  std::complex<double> z(0, 0);
  for (double t : theta) z += std::polar(1.0, t);
  return z / static_cast<double>(theta.size());
}

struct Kuramoto {
  std::vector<double> omega;
  double coupling = 0;            // K
  calx::linalg::Matrix weights;   // A; empty (order 0) for all to all
  void operator()(double, const State& theta, State& d) const {
    const std::size_t n = theta.size();
    if (weights.n == 0) {
      const auto z = order_parameter(theta);
      const double r = std::abs(z), psi = std::arg(z);
      for (std::size_t i = 0; i < n; i++) d[i] = omega[i] + coupling * r * std::sin(psi - theta[i]);
      return;
    }
    for (std::size_t i = 0; i < n; i++) {
      double s = 0;
      for (std::size_t j = 0; j < n; j++)
        if (weights(i, j) != 0) s += weights(i, j) * std::sin(theta[j] - theta[i]);
      d[i] = omega[i] + coupling * s;
    }
  }
};

inline double kuramoto_critical_coupling(double density_at_zero) { return 2 / (std::numbers::pi * density_at_zero); }

// r of the Lorentzian of half-width γ at coupling K: √(1 − 2γ/K) above K_c = 2γ, 0 below.
inline double kuramoto_lorentzian_order(double K, double gamma) { return K > 2 * gamma ? std::sqrt(1 - 2 * gamma / K) : 0.0; }

inline std::vector<double> lorentzian_frequencies(std::size_t n, double centre, double gamma) {
  std::vector<double> w(n);
  for (std::size_t i = 0; i < n; i++) w[i] = centre + gamma * std::tan(std::numbers::pi * (static_cast<double>(i) + 0.5) / static_cast<double>(n) - std::numbers::pi / 2);
  return w;
}

}  // namespace physis::dynamics
