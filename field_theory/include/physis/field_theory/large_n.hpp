// The large-N phases of invariant field theories: the resummed dominant graphs, their critical
// points and exponents (physis decision 0006, §2.1, stage 4).
//
// - The vector model (rank 1), S = φ̄·φ + (λ/N)(φ̄·φ)²: at large N the dominant graphs are the
//   cacti, and the two-point function G = ⟨φ̄·φ⟩/N solves the gap equation
//   G (μ + 2λG) = 1 (Stanley 1968; Moshe and Zinn-Justin 2003, Phys. Rep. 385, 69, §2), with
//   ∂_λ (log Z)/N = −G² at leading order, so that the free energy's coefficients of (−λ)^n are
//   2^{n−1} Cat_n / n.
// - The Hermitian matrix model (rank 2), Z = ∫ dM e^{−N tr(M²/2 + g M⁴)}: at large N the
//   planar graphs, whose sum is fixed by the eigenvalue density on [−2a, 2a] (Brézin, Itzykson,
//   Parisi and Zuber 1978, Commun. Math. Phys. 59, 35, eqs. 16–23):
//     ρ(λ) = (1/π)(1/2 + 4g a² + 2g λ²) √(4a² − λ²),   12 g a⁴ + a² − 1 = 0,
//     E₀(g) = −lim N⁻² log Z/Z₀ = (a² − 1)(9 − a²)/24 − (1/2) log a²
//           = −Σ_{k>=1} (−12 g)^k (2k − 1)! / (k! (k + 2)!).
//   Its critical point g_c = −1/48 (the coupling of g M⁴/4 is then −1/12) is the continuum limit
//   of random surfaces: E₀ ~ (g − g_c)^{2−γ} with the string susceptibility γ = −1/2, pure
//   two-dimensional quantum gravity (Kazakov 1985; David 1985; Di Francesco, Ginsparg and
//   Zinn-Justin 1995, §2.5).
// - Tensor models (rank D >= 3): the melonic graphs dominate, and the two-point function solves
//   G = 1 + x G^m: m = 2, x = −2λ for the quartic melonic model with one bubble (either black
//   vertex of the bubble may be the external one), −2Dλ with one of each colour; m = D + 1,
//   x = λλ̄ for the coloured simplicial model (Bonzom, Gurau, Riello and
//   Rivasseau 2011, Nucl. Phys. B 853, 174). Its coefficients are the Fuss–Catalan numbers
//   (1/(mp + 1)) C(mp + 1, p); its critical point x_c = (m − 1)^{m−1}/m^m, G_c = m/(m − 1),
//   with G − G_c ~ (x_c − x)^{1/2}: branched polymers, γ = 1/2 (Gurau and Ryan 2014).
//
// Exponents from series: if c_k ~ C x_c^{−k} k^{−α}, then α is the limit of
// k (1 − (c_k / c_{k−1}) x_c) corrected by Richardson extrapolation; γ = 3 − α for a free
// energy and γ = 2 − α for a two-point function.

#pragma once

#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <vector>

#include <calx/quadrature/gauss_kronrod.hpp>
#include <calx/solvers/roots.hpp>
#include <morphe/arithmetic/rational.hpp>

namespace physis::field_theory {

using morphe::arithmetic::Integer;
using morphe::arithmetic::Rational;

// ---- the vector model ----------------------------------------------------------------------------------------------

// The large-N two-point function: the root of G (μ + 2λG) = 1 that tends to 1/μ as λ → 0.
inline double vector_two_point(double lambda, double mu = 1) {
  if (!(mu > 0)) throw std::invalid_argument("physis::field_theory::vector_two_point: μ > 0");
  if (lambda == 0) return 1 / mu;
  const double d = mu * mu + 8 * lambda;
  if (d < 0) throw std::invalid_argument("physis::field_theory::vector_two_point: beyond the critical coupling λ_c = −μ²/8");
  return 2 / (mu + std::sqrt(d));  // = (−μ + √d)/(4λ), without cancellation
}
// The vector model's critical coupling, where the two roots meet.
inline double vector_critical(double mu = 1) { return -mu * mu / 8; }

// The leading free energy's coefficients of (−λ)^n, divided by N: 2^{n−1} Cat_n / n, n >= 1.
inline std::vector<Rational> vector_free_energy(std::size_t order) {
  std::vector<Rational> w(order + 1, Rational(0));
  for (std::size_t n = 1; n <= order; n++) {
    const Integer cat = morphe::arithmetic::exact_div(morphe::arithmetic::binomial(2 * static_cast<long long>(n), static_cast<long long>(n)), Integer(static_cast<long long>(n) + 1));
    w[n] = Rational(morphe::arithmetic::pow(Integer(2), n - 1) * cat, Integer(static_cast<long long>(n)));
  }
  return w;
}

// ---- the Hermitian quartic matrix model (BIPZ) ------------------------------------------------------------------------

inline constexpr double bipz_critical = -1.0 / 48;

// a²(g), the root of 12 g a⁴ + a² − 1 = 0 continuous from a² = 1 at g = 0, for g >= g_c.
inline double bipz_a2(double g) {
  if (g < bipz_critical) throw std::invalid_argument("physis::field_theory::bipz_a2: g below the critical coupling −1/48");
  if (g == 0) return 1;
  return 2 / (1 + std::sqrt(1 + 48 * g));  // = (−1 + √(1 + 48g))/(24g)
}
// The planar eigenvalue density at λ (zero outside [−2a, 2a]).
inline double bipz_density(double lambda, double g) {
  const double a2 = bipz_a2(g), d = 4 * a2 - lambda * lambda;
  if (d <= 0) return 0;
  return (0.5 + 4 * g * a2 + 2 * g * lambda * lambda) * std::sqrt(d) / std::numbers::pi;
}
// The planar free energy E₀(g) − E₀(0).
inline double bipz_free_energy(double g) {
  const double a2 = bipz_a2(g);
  return (a2 - 1) * (9 - a2) / 24 - 0.5 * std::log(a2);
}
// Its Taylor coefficients e_1 .. e_order, exactly: e_k = −(−12)^k (2k − 1)! / (k! (k + 2)!).
inline std::vector<Rational> bipz_series(std::size_t order) {
  std::vector<Rational> e(order + 1, Rational(0));
  for (std::size_t k = 1; k <= order; k++) {
    const Integer num = morphe::arithmetic::pow(Integer(-12), k) * morphe::arithmetic::factorial(2 * k - 1);
    e[k] = -Rational(num, morphe::arithmetic::factorial(k) * morphe::arithmetic::factorial(k + 2));
  }
  return e;
}
// The planar moment ∫ λ^{2j} ρ(λ) dλ, by quadrature.
inline double bipz_moment(int two_j, double g) {
  const double a = std::sqrt(bipz_a2(g));
  return calx::quadrature::integrate([&](double l) { return std::pow(l, two_j) * bipz_density(l, g); }, -2 * a, 2 * a, 1e-13, 1e-12).value;
}

// ---- melonic resummation ---------------------------------------------------------------------------------------------

// x_c = (m − 1)^{m−1} / m^m, the critical point of G = 1 + x G^m.
inline double melonic_critical(int m) {
  if (m < 2) throw std::invalid_argument("physis::field_theory::melonic_critical: m >= 2");
  return std::pow(m - 1.0, m - 1) / std::pow(static_cast<double>(m), m);
}
// The branch of G = 1 + x G^m analytic at x = 0, for x <= x_c (Brent's method on a bracket).
inline double melonic_two_point(double x, int m) {
  const double xc = melonic_critical(m);
  if (x > xc * (1 + 1e-14)) throw std::invalid_argument("physis::field_theory::melonic_two_point: x beyond the critical point");
  if (x == 0) return 1;
  auto f = [&](double G) { return 1 + x * std::pow(G, m) - G; };
  if (x < 0) return calx::solvers::brent(f, 0.0, 1.0).x;
  const double Gc = m / (m - 1.0);
  if (f(Gc) > 0) return Gc;  // at x_c, to rounding
  return calx::solvers::brent(f, 1.0, Gc).x;
}
// Its Taylor coefficients, the Fuss–Catalan numbers (1/(mp + 1)) C(mp + 1, p), p = 0..order.
inline std::vector<Integer> melonic_coefficients(int m, std::size_t order) {
  std::vector<Integer> c;
  for (std::size_t p = 0; p <= order; p++) {
    const long long mp = static_cast<long long>(m) * static_cast<long long>(p);
    c.push_back(morphe::arithmetic::exact_div(morphe::arithmetic::binomial(mp + 1, static_cast<long long>(p)), Integer(mp + 1)));
  }
  return c;
}

// ---- exponents from series --------------------------------------------------------------------------------------------

// α in c_k ~ C x_c^{−k} k^{−α}, from the terms k <= K of a series whose critical point is known:
// a_k = k (1 − x_c c_k / c_{k−1}) → α + O(1/k), extrapolated by Richardson in 1/k to second order.
inline double exponent_from_series(const std::vector<double>& c, double xc) {
  const std::size_t K = c.size() - 1;
  if (K < 8) throw std::invalid_argument("physis::field_theory::exponent_from_series: at least 8 terms");
  auto a = [&](std::size_t k) { return static_cast<double>(k) * (1 - xc * c[k] / c[k - 1]); };
  // a_k = α + b/k + c/k²: eliminate b and c with k, k − 1, k − 2.
  const double k = static_cast<double>(K);
  const double a0 = a(K), a1 = a(K - 1), a2 = a(K - 2);
  return (k * k * a0 - 2 * (k - 1) * (k - 1) * a1 + (k - 2) * (k - 2) * a2) / 2;
}

}  // namespace physis::field_theory
