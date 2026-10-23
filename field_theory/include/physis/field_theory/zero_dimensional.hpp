// Zero-dimensional field theory: the φ⁴ integral and its perturbation series.
//
// The "field theory in zero dimensions" (Hurst 1952; Zinn-Justin 2002, Quantum Field Theory
// and Critical Phenomena, 4th ed., §9.1; Kleinert and Schulte-Frohlinde 2001, Critical
// Properties of φ⁴-Theories, ch. 16) is the integral
//
//   Z(λ) = (1/√(2π)) ∫ dφ exp(-½ m² φ² - (λ/4!) φ⁴),     m² = 1 by default,
//
// the partition function of φ⁴ theory on one point. Expanded in λ by Wick's theorem,
// ⟨φ^{4n}⟩₀ = (4n - 1)!!, so
//
//   Z(λ) = Σ_n z_n λ^n,  z_n = (-1)^n (4n - 1)!! / ((4!)^n n!),
//
// a series with zero radius of convergence (Dyson 1952): its terms decrease until
// n ≈ 3/(2λ) and then grow factorially. Z is completely monotone in λ > 0 (it is the Laplace
// transform of a positive measure), so the error of the truncated series is bounded by the
// first omitted term, and optimal truncation leaves an error of order e^{-3/(2λ)}. The
// coefficients of log Z = W(λ), the sum over connected vacuum graphs, follow exactly from
// those of Z by the logarithm of a power series (feynman.hpp computes the same numbers from
// the graphs).

#pragma once

#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <vector>

#include <calx/quadrature/gauss_kronrod.hpp>
#include <morphe/arithmetic/rational.hpp>

namespace physis::field_theory {

using morphe::arithmetic::Integer;
using morphe::arithmetic::Rational;

namespace detail {

// ∫ f over |φ| <= 40/√m² (beyond which e^{-m² φ²/2} < 1e-300), on pieces of width 1/(2√m²)
// so that the adaptive rule sees the peak at 0 whatever the integrand (a rule on the whole
// line places one node near 0, where φ^{2k} e^{-S} vanishes, and is fooled).
template <class F>
double integrate_line(F&& f, double m2) {
  const double h = 0.5 / std::sqrt(m2);
  double s = 0;
  for (int k = -80; k < 80; k++) s += calx::quadrature::integrate(f, k * h, (k + 1) * h, 1e-300, 1e-14, 200).value;
  return s;
}

}  // namespace detail

// Z(λ) by adaptive quadrature.
inline double zero_dimensional_partition(double lambda, double m2 = 1) {
  if (lambda < 0 || !(m2 > 0)) throw std::invalid_argument("physis::field_theory::zero_dimensional_partition: λ >= 0, m² > 0");
  auto f = [&](double p) { return std::exp(-0.5 * m2 * p * p - lambda / 24 * p * p * p * p); };
  return detail::integrate_line(f, m2) / std::sqrt(2 * std::numbers::pi);
}

// ⟨φ^{2k}⟩ in the interacting zero-dimensional theory, by quadrature.
inline double zero_dimensional_moment(int two_k, double lambda, double m2 = 1) {
  if (lambda < 0 || !(m2 > 0)) throw std::invalid_argument("physis::field_theory::zero_dimensional_moment: λ >= 0, m² > 0");
  auto w = [&](double p) { return std::exp(-0.5 * m2 * p * p - lambda / 24 * p * p * p * p); };
  const double num = detail::integrate_line([&](double p) { return std::pow(p, two_k) * w(p); }, m2);
  return num / detail::integrate_line(w, m2);
}

// The coefficients z_0..z_N of Z(λ) at m² = 1, exactly.
inline std::vector<Rational> partition_series(std::size_t N) {
  std::vector<Rational> z;
  for (std::size_t n = 0; n <= N; n++) {
    const Integer num = morphe::arithmetic::double_factorial(static_cast<long long>(4 * n) - 1);
    const Integer den = morphe::arithmetic::pow(Integer(24), n) * morphe::arithmetic::factorial(n);
    z.push_back(Rational(n % 2 ? -num : num, den));
  }
  return z;
}

// The logarithm of a power series a(λ) with a_0 = 1: b = log a, by n b_n = n a_n - Σ_{k<n} k b_k a_{n-k}.
inline std::vector<Rational> log_series(const std::vector<Rational>& a) {
  if (a.empty() || a[0] != Rational(1)) throw std::invalid_argument("physis::field_theory::log_series: a_0 = 1");
  std::vector<Rational> b(a.size(), Rational(0));
  for (std::size_t n = 1; n < a.size(); n++) {
    Rational s = Rational(static_cast<long long>(n)) * a[n];
    for (std::size_t k = 1; k < n; k++) s -= Rational(static_cast<long long>(k)) * b[k] * a[n - k];
    b[n] = s / Rational(static_cast<long long>(n));
  }
  return b;
}

// The partial sum Σ_{n <= N} c_n λ^n of a series.
inline double partial_sum(const std::vector<Rational>& c, std::size_t N, double lambda) {
  double s = 0, p = 1;
  for (std::size_t n = 0; n <= N && n < c.size(); n++, p *= lambda) s += c[n].to_double() * p;
  return s;
}

}  // namespace physis::field_theory
