// The fluctuations of currents in Markov jump processes, the thermodynamic uncertainty relation, and Landauer's
// bound: the physical floors of the cost of a computation (decision 0022).
//
// A current. In a chain with rates k_ij (i → j, as in master_equation.hpp), a current counts each jump i → j
// with an increment d_ij = −d_ji: J_t = Σ over the jumps up to t of d. In the stationary state its mean grows
// at the rate J̄ = Σ_{i≠j} π_i k_ij d_ij and its variance at a rate D₂ = lim Var J_t / t. Both are derivatives at
// s = 0 of the scaled cumulant generating function λ(s), the largest eigenvalue of the tilted generator
// (L_s f)(i) = Σ_j k_ij e^{s d_ij} f(j) − k_i f(i) (Lebowitz and Spohn 1999, J. Stat. Phys. 95, 333; Touchette
// 2009, Phys. Rep. 478, 1, §3.4). Expanding L_s = L + s L₁ + (s²/2) L₂ about its right eigenvector 1 and left
// eigenvector π:
//
//   λ′(0) = π L₁ 1 = J̄,   λ″(0) = π L₂ 1 + 2 π L₁ r₁,   L r₁ = J̄ 1 − L₁ 1,  π·r₁ = 0,
//
// with (L₁)_ij = k_ij d_ij, (L₂)_ij = k_ij d_ij². D₂ = λ″(0). The linear system is solved densely with one
// equation replaced by π·r₁ = 0.
//
// The thermodynamic uncertainty relation (Barato and Seifert 2015, Phys. Rev. Lett. 114, 158101; proved for
// every current of every chain in its stationary state by Gingrich, Horowitz, Perunov and England 2016, Phys.
// Rev. Lett. 116, 120601):
//
//   D₂ σ / (2 J̄²) ≥ 1,   σ the entropy production per unit time in units of k_B,
//
// so that measuring a current to a relative precision ε (ε² = Var J_t / J̄² t² after a time t) costs a
// dissipation of at least 2 k_B T / ε². tur_ratio returns the left side; its equality is approached near
// equilibrium.
//
// Landauer's bound (Landauer 1961, IBM J. Res. Dev. 5, 183; Bérut et al. 2012, Nature 483, 187): erasing one
// bit at temperature T dissipates at least k_B T ln 2.
//
// Dense: O(n³), for chains of a few thousand states at most. Checks: tests/check_fluctuations.cpp.
#pragma once

#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <vector>

#include <calx/linalg/dense.hpp>

#include <physis/units/constants.hpp>

#include "master_equation.hpp"

namespace physis::statistical_mechanics {

struct CurrentStatistics {
  double mean_rate = 0;        // J̄
  double variance_rate = 0;    // D₂ = lim Var J_t / t
  double entropy_production = 0;  // σ, in units of k_B per unit time
  double tur_ratio = 0;        // D₂ σ / (2 J̄²) ≥ 1; infinite when J̄ = 0
};

// k: rates i → j (k_ij > 0 ⇔ k_ji > 0, irreducible); d: the increments, antisymmetric.
inline CurrentStatistics current_statistics(const calx::linalg::Matrix& k, const calx::linalg::Matrix& d) {
  const std::size_t n = k.n;
  if (d.n != n) throw std::invalid_argument("physis::statistical_mechanics::current_statistics: one increment per pair of states");
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++)
      if (std::fabs(d(i, j) + d(j, i)) > 1e-12 * (1 + std::fabs(d(i, j)))) throw std::invalid_argument("physis::statistical_mechanics::current_statistics: increments antisymmetric, d_ij = −d_ji");
  const auto irr = irreversibility(k);
  const auto& pi = irr.pi;
  CurrentStatistics out;
  out.entropy_production = irr.entropy_production;
  std::vector<double> L1one(n, 0.0);
  double J = 0, pL2 = 0;
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++)
      if (i != j && k(i, j) > 0) {
        L1one[i] += k(i, j) * d(i, j);
        J += pi[i] * k(i, j) * d(i, j);
        pL2 += pi[i] * k(i, j) * d(i, j) * d(i, j);
      }
  out.mean_rate = J;
  // L r₁ = J 1 − L₁ 1 with π·r₁ = 0: the last equation replaced by the gauge.
  calx::linalg::Matrix A(n), b(n);
  for (std::size_t i = 0; i < n; i++) {
    double out_rate = 0;
    for (std::size_t j = 0; j < n; j++)
      if (j != i) A(i, j) = k(i, j), out_rate += k(i, j);
    A(i, i) = -out_rate;
    b(i, 0) = J - L1one[i];
  }
  for (std::size_t j = 0; j < n; j++) A(n - 1, j) = pi[j];
  b(n - 1, 0) = 0;
  const auto r = calx::linalg::solve(A, b);
  double pL1r = 0;
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++)
      if (i != j && k(i, j) > 0) pL1r += pi[i] * k(i, j) * d(i, j) * r(j, 0);
  out.variance_rate = pL2 + 2 * pL1r;
  out.tur_ratio = J == 0 ? INFINITY : out.variance_rate * out.entropy_production / (2 * J * J);
  return out;
}

// k_B T ln 2 per bit erased, J.
inline double landauer_bound(double temperature) { return physis::units::si::boltzmann * temperature * std::numbers::ln2; }

// The least dissipation (J) of measuring a current to relative precision ε at temperature T: 2 k_B T / ε².
inline double precision_cost_bound(double temperature, double relative_precision) {
  return 2 * physis::units::si::boltzmann * temperature / (relative_precision * relative_precision);
}

}  // namespace physis::statistical_mechanics
