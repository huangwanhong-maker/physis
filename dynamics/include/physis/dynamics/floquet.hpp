// Linear systems with periodic coefficients: the monodromy matrix and its Floquet multipliers
// (decision 0022).
//
// For x' = A(t) x with A(t + T) = A(t), the fundamental matrix Φ(t) (Φ(0) = I) satisfies
// Φ(t + T) = Φ(t) M, M = Φ(T) the monodromy matrix (Floquet 1883; Coddington and Levinson 1955, Theory
// of Ordinary Differential Equations, ch. 3, §5). Its eigenvalues μ are the Floquet multipliers: the
// zero solution is asymptotically stable when every |μ| < 1 and unstable when one |μ| > 1, and a
// multiplier −1 marks a subharmonic response, of period 2T. Liouville's formula gives
// det M = exp ∫₀ᵀ tr A(t) dt, a check that holds whatever A's eigenvalues.
//
// The Mathieu equation x'' + (δ + ε cos t) x = 0, written as a first-order system, is the type: its
// first tongue of instability leaves δ = 1/4 between δ = 1/4 ± ε/2 − ε²/8 + O(ε³) (from the
// boundaries a = 1 ± q − q²/8 of the standard form y'' + (a − 2q cos 2s) y = 0; Abramowitz and Stegun
// 1964, §20.2.25). A damping term 2γx' lifts the tongue: in x'' + 2γx' + ω₀²(1 − ε cos 2ω₀t) x = 0
// the first tongue's tip is at ε = 4γ/ω₀ to first order in γ (Landau and Lifshitz, Mechanics, §27).
// Those are the Faraday instability's equations, mode by mode (Benjamin and Ursell 1954).
//
// The fundamental matrix is integrated by the classical Runge–Kutta method in a given number of steps
// over the period; the error in M is O((T/steps)⁴). Checks: tests/check_chaos.cpp.
#pragma once

#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

#include <calx/linalg/dense.hpp>
#include <calx/linalg/eigenvalues.hpp>

#include "flows.hpp"

namespace physis::dynamics {

struct Floquet {
  calx::linalg::Matrix monodromy;
  std::vector<std::complex<double>> multipliers;
  double spectral_radius = 0;  // max |μ|
  double log_growth = 0;       // ln(max |μ|) / T: the Floquet exponent with the largest real part
};

// A(t, Matrix&) writes A at time t; n the dimension, T the period.
template <class A>
Floquet floquet(A&& coefficients, std::size_t n, double period, std::size_t steps = 2000) {
  calx::linalg::Matrix a(n);
  auto field = [&](double t, const State& y, State& dy) {
    coefficients(t, a);
    // y holds Φ column by column: dΦ/dt = A Φ.
    for (std::size_t c = 0; c < n; c++)
      for (std::size_t i = 0; i < n; i++) {
        double s = 0;
        for (std::size_t j = 0; j < n; j++) s += a(i, j) * y[c * n + j];
        dy[c * n + i] = s;
      }
  };
  State y(n * n, 0.0);
  for (std::size_t c = 0; c < n; c++) y[c * n + c] = 1;
  const double h = period / static_cast<double>(steps);
  for (std::size_t s = 0; s < steps; s++) rk4_step(field, h * static_cast<double>(s), y, h);
  Floquet out;
  out.monodromy = calx::linalg::Matrix(n);
  for (std::size_t c = 0; c < n; c++)
    for (std::size_t i = 0; i < n; i++) out.monodromy(i, c) = y[c * n + i];
  out.multipliers = calx::linalg::eigenvalues(out.monodromy);
  for (const auto& m : out.multipliers) out.spectral_radius = std::max(out.spectral_radius, std::abs(m));
  out.log_growth = std::log(out.spectral_radius) / period;
  return out;
}

// x'' + 2γ x' + (w0² − w1² cos Ωt) x = 0 as a first-order system: a damped, parametrically driven
// oscillator (the Mathieu equation when γ = 0).
inline Floquet parametric_oscillator(double gamma, double w0_squared, double w1_squared, double Omega, std::size_t steps = 2000) {
  return floquet([&](double t, calx::linalg::Matrix& a) {
    a(0, 0) = 0, a(0, 1) = 1;
    a(1, 0) = -(w0_squared - w1_squared * std::cos(Omega * t)), a(1, 1) = -2 * gamma;
  }, 2, 2 * 3.141592653589793238462643383279502884 / Omega, steps);
}

}  // namespace physis::dynamics
