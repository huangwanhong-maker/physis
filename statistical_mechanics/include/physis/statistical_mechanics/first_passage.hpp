// First passage and absorption in Markov jump processes: phase-type laws, absorption probabilities
// and mean first-passage times (decision 0016).
//
// The law (Neuts 1981, Matrix-Geometric Solutions in Stochastic Models, ch. 2; Latouche and Ramaswami
// 1999, Introduction to Matrix Analytic Methods in Stochastic Modeling, ch. 2). Rates k_ij as in
// master_equation.hpp. Given a set of transient states (in a chosen order), let T be the generator
// restricted to them: T_ab = k_{t_a t_b} (a ≠ b), T_aa = −Σ_{j≠t_a} k_{t_a j}, the total rate out of
// t_a in the whole chain; t = −T 1 is then the rate from each transient state to the states outside.
// Started from α on the transient states (Σ α ≤ 1, the defect an atom at time 0), the time τ of the
// first exit from them has the phase-type law
//
//   density α exp(Tτ) t,   survival P(τ > s) = α exp(Ts) 1,   E τᵏ = k! α (−T)^{−k} 1
//
// (Neuts 1981, ch. 2; Latouche and Ramaswami 1999, ch. 2), finite exactly when every transient state
// reaches the outside, i.e. T is non-singular. The law of the state reached first outside the transient set, from transient state a,
// is row a of (−T)⁻¹ R, R_ab = k_{t_a, b} the rates into the states outside (the embedded chain's
// absorption probabilities written with rates). The mean time to reach a target, m_i, solves
// Σ_j k_ij (m_j − m_i) = −1 for i ≠ target, m_target = 0: the phase-type mean with every other state
// transient.
//
// Accuracy. Dense: O(n³) per expm and per solve; quantile() bisects the survival function to a
// relative 10⁻¹³ (some fifty exponentials). Checks: tests/check_first_passage.cpp.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include <calx/linalg/dense.hpp>
#include <calx/linalg/expm.hpp>
#include <calx/linalg/rectangular.hpp>
#include <physis/statistical_mechanics/master_equation.hpp>

namespace physis::statistical_mechanics {

struct PhaseType {
  calx::linalg::Matrix T;
  std::vector<double> alpha, exit;  // the initial law on the transient states; t = −T 1

  // α exp(Tτ) t (0 for τ < 0); the atom at 0, 1 − Σ α, is not part of it.
  double density(double tau) const {
    if (tau < 0) return 0;
    const auto v = row_exp(tau);
    double d = 0;
    for (std::size_t a = 0; a < v.size(); a++) d += v[a] * exit[a];
    return d;
  }
  // α exp(Tτ) 1 = P(time > τ) (1 for τ < 0).
  double survival(double tau) const {
    if (tau < 0) return 1;
    const auto v = row_exp(tau);
    double s = 0;
    for (double x : v) s += x;
    return s;
  }
  // E τᵏ = k! α (−T)^{−k} 1 (the atom at 0 counted; E τ⁰ = 1).
  double moment(int k) const {
    if (k < 0) throw std::invalid_argument("physis::statistical_mechanics::PhaseType::moment: k ≥ 0");
    if (k == 0) return 1;
    calx::linalg::Matrix minus = -1.0 * T;
    const calx::linalg::Matrix inv = calx::linalg::inverse(minus);
    std::vector<double> x(T.n, 1.0);
    double factorial = 1;
    for (int r = 1; r <= k; r++) {
      x = calx::linalg::multiply(inv, x);
      factorial *= r;
    }
    double m = 0;
    for (std::size_t a = 0; a < T.n; a++) m += alpha[a] * x[a];
    return factorial * m;
  }
  // The least τ with P(time ≤ τ) ≥ p, p ∈ [0, 1), by bisection on the survival function.
  double quantile(double p) const {
    if (!(p >= 0 && p < 1)) throw std::invalid_argument("physis::statistical_mechanics::PhaseType::quantile: p in [0, 1)");
    if (p <= atom()) return 0;
    const auto cdf = [&](double s) { return 1 - survival(s); };
    double hi = std::max(moment(1), 1e-300), lo = 0;
    for (int i = 0; cdf(hi) < p; i++) {
      if (i > 2000) throw std::runtime_error("physis::statistical_mechanics::PhaseType::quantile: no bracket");
      lo = hi, hi *= 2;
    }
    while (hi - lo > 1e-13 * hi) {
      const double mid = 0.5 * (lo + hi);
      (cdf(mid) < p ? lo : hi) = mid;
    }
    return hi;
  }
  // 1 − Σ α: the probability of starting outside the transient states.
  double atom() const {
    double s = 0;
    for (double a : alpha) s += a;
    return std::max(0.0, 1 - s);
  }

 private:
  std::vector<double> row_exp(double tau) const {
    const calx::linalg::Matrix E = calx::linalg::expm(tau * T);
    std::vector<double> v(T.n, 0.0);
    for (std::size_t a = 0; a < T.n; a++) {
      if (alpha[a] == 0) continue;
      for (std::size_t b = 0; b < T.n; b++) v[b] += alpha[a] * E(a, b);
    }
    return v;
  }
};

namespace detail {

// The positions of `transient` among the states (n where absent), checking that they are distinct
// and in range, and that every one of them reaches a state outside.
inline std::vector<std::size_t> transient_slots(const calx::linalg::Matrix& k, std::span<const std::size_t> transient, const char* where) {
  check_rates(k, where);
  std::vector<std::size_t> slot(k.n, k.n);
  for (std::size_t a = 0; a < transient.size(); a++) {
    if (transient[a] >= k.n || slot[transient[a]] != k.n) throw std::invalid_argument(std::string(where) + ": transient states distinct and in range");
    slot[transient[a]] = a;
  }
  std::vector<char> outside(k.n, 0);
  for (std::size_t i = 0; i < k.n; i++) outside[i] = slot[i] == k.n;
  const auto reach = reaching(k, outside);
  for (std::size_t s : transient)
    if (!reach[s]) throw std::invalid_argument(std::string(where) + ": a transient state never leaves the transient states");
  return slot;
}

inline calx::linalg::Matrix restricted(const calx::linalg::Matrix& k, std::span<const std::size_t> transient) {
  const calx::linalg::Matrix Q = generator(k);
  calx::linalg::Matrix T(transient.size());
  for (std::size_t a = 0; a < transient.size(); a++)
    for (std::size_t b = 0; b < transient.size(); b++) T(a, b) = Q(transient[a], transient[b]);
  return T;
}

}  // namespace detail

// The phase-type law of the time to leave `transient`, started from alpha on them.
inline PhaseType phase_type(const calx::linalg::Matrix& k, std::span<const std::size_t> transient, std::span<const double> alpha) {
  const auto slot = detail::transient_slots(k, transient, "physis::statistical_mechanics::phase_type");
  if (alpha.size() != transient.size()) throw std::invalid_argument("physis::statistical_mechanics::phase_type: one initial probability for each transient state");
  double total = 0;
  for (double a : alpha) {
    if (!(a >= 0)) throw std::invalid_argument("physis::statistical_mechanics::phase_type: initial probabilities ≥ 0");
    total += a;
  }
  if (total > 1 + 1e-12) throw std::invalid_argument("physis::statistical_mechanics::phase_type: initial probabilities summing to at most 1");
  PhaseType P;
  P.T = detail::restricted(k, transient);
  P.alpha.assign(alpha.begin(), alpha.end());
  P.exit.assign(transient.size(), 0.0);
  for (std::size_t a = 0; a < transient.size(); a++)
    for (std::size_t j = 0; j < k.n; j++)
      if (slot[j] == k.n) P.exit[a] += k(transient[a], j);
  return P;
}

// Row a: the law of the first state outside `transient` reached from transient[a], over the states
// `absorbing` (columns in their order): (−T)⁻¹ R. A row sums to less than 1 when the chain can leave
// the transient states into a state not listed.
inline calx::linalg::Rectangular absorption_probabilities(const calx::linalg::Matrix& k, std::span<const std::size_t> transient, std::span<const std::size_t> absorbing) {
  const auto slot = detail::transient_slots(k, transient, "physis::statistical_mechanics::absorption_probabilities");
  std::vector<char> seen(k.n, 0);
  for (std::size_t b : absorbing) {
    if (b >= k.n || slot[b] != k.n || seen[b]++) throw std::invalid_argument("physis::statistical_mechanics::absorption_probabilities: absorbing states distinct, in range and not transient");
  }
  const calx::linalg::Matrix T = detail::restricted(k, transient);
  calx::linalg::Rectangular R(transient.size(), absorbing.size());
  for (std::size_t a = 0; a < transient.size(); a++)
    for (std::size_t b = 0; b < absorbing.size(); b++) R(a, b) = k(transient[a], absorbing[b]);
  return calx::linalg::Rectangular::from_square(calx::linalg::inverse(-1.0 * T)) * R;
}

// m_i, the mean time to reach `target` from i (m_target = 0); throws unless every state reaches it.
inline std::vector<double> mean_first_passage(const calx::linalg::Matrix& k, std::size_t target) {
  detail::check_rates(k, "physis::statistical_mechanics::mean_first_passage");
  if (target >= k.n) throw std::invalid_argument("physis::statistical_mechanics::mean_first_passage: no such target");
  std::vector<std::size_t> others;
  for (std::size_t i = 0; i < k.n; i++)
    if (i != target) others.push_back(i);
  std::vector<double> m(k.n, 0.0);
  if (others.empty()) return m;
  detail::transient_slots(k, others, "physis::statistical_mechanics::mean_first_passage");
  const calx::linalg::Matrix inv = calx::linalg::inverse(-1.0 * detail::restricted(k, others));
  const std::vector<double> x = calx::linalg::multiply(inv, std::vector<double>(others.size(), 1.0));
  for (std::size_t a = 0; a < others.size(); a++) m[others[a]] = x[a];
  return m;
}

}  // namespace physis::statistical_mechanics
