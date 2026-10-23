// The recoupling theory of SU(2), exactly: Clebsch–Gordan coefficients and the Wigner 3j, 6j
// and 9j symbols.
//
// Spins and their projections are passed as twice their values (two_j, two_m; physis decision
// 0005, §2.4), so that half-integers are exact integers. Every symbol is a signed square root
// of a rational number, and is returned exactly as one (`SignedSqrt`), computed on morphe's
// integers; value() gives the double.
//
// - **Clebsch–Gordan coefficients** ⟨j₁ m₁ j₂ m₂ | J M⟩ by Racah's single-sum formula
//   (Racah 1942, Phys. Rev. 62, 438, eq. 16'; Varshalovich, Moskalev and Khersonskii 1988,
//   Quantum Theory of Angular Momentum, §8.2.1, eq. 3), with the Condon–Shortley phase.
// - **3j symbols** (j₁ j₂ j₃; m₁ m₂ m₃) = (−1)^{j₁−j₂−m₃} (2j₃ + 1)^{−1/2} ⟨j₁ m₁ j₂ m₂ | j₃ −m₃⟩
//   (Edmonds 1957, Angular Momentum in Quantum Mechanics, eq. 3.7.3).
// - **6j symbols** {j₁ j₂ j₃; j₄ j₅ j₆} by Racah's formula
//   Δ(j₁j₂j₃) Δ(j₁j₅j₆) Δ(j₄j₂j₆) Δ(j₄j₅j₃) Σ_z (−1)^z (z+1)! / [Π_i (z − α_i)! Π_k (β_k − z)!],
//   α = (j₁+j₂+j₃, j₁+j₅+j₆, j₄+j₂+j₆, j₄+j₅+j₃), β = (j₁+j₂+j₄+j₅, j₂+j₃+j₅+j₆, j₃+j₁+j₆+j₄),
//   Δ(abc)² = (a+b−c)! (a−b+c)! (−a+b+c)! / (a+b+c+1)! (Racah 1942, eq. 36; Edmonds eq. 6.3.7).
//   It vanishes unless its four triads (j₁j₂j₃), (j₁j₅j₆), (j₄j₂j₆), (j₄j₅j₃) satisfy the
//   triangle conditions; it is the tetrahedral symbol of the tetrahedron whose opposite edges
//   are (j₁, j₄), (j₂, j₅), (j₃, j₆).
// - **9j symbols** by the sum over three 6j symbols (Edmonds eq. 6.4.3),
//   Σ_x (−1)^{2x} (2x+1) {j₁ j₄ j₇; j₈ j₉ x} {j₂ j₅ j₈; j₄ x j₆} {j₃ j₆ j₉; x j₁ j₂}.
//   The triads containing x appear twice in each product, so their Δ factors square to
//   rationals and the 9j symbol is the square root of the product of the Δ² of its six rows
//   and columns times a rational: it is exact as well.
//
// Sums of these symbols (orthogonality, the Racah and Biedenharn–Elliott identities) are sums
// of square roots of rationals whose ratios are squares; `sum` adds such surds exactly, and
// refuses incommensurable ones.
//
// Cost: the 6j symbol is O(min β − max α) terms of integers of O(j log j) bits. The table of
// factorials is kept per thread.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <deque>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <morphe/arithmetic/integer.hpp>
#include <morphe/arithmetic/rational.hpp>

namespace physis::symmetry {

using morphe::arithmetic::Integer;
using morphe::arithmetic::Rational;

// ⌊√n⌋ for n ≥ 0, by Newton's iteration on integers.
inline Integer isqrt(const Integer& n) {
  if (n.sign() < 0) throw std::domain_error("physis::symmetry::isqrt: negative argument");
  if (n.is_zero()) return 0;
  Integer x = morphe::arithmetic::pow(Integer(2), (n.bit_length() + 1) / 2 + 1);
  for (;;) {
    Integer y = (x + n / x) / 2;
    if (y >= x) return x;
    x = std::move(y);
  }
}
// √q when q is the square of a rational.
inline std::optional<Rational> rational_sqrt(const Rational& q) {
  if (q.sign() < 0) return std::nullopt;
  const Integer a = isqrt(q.numerator()), b = isqrt(q.denominator());
  if (a * a != q.numerator() || b * b != q.denominator()) return std::nullopt;
  return Rational(a, b);
}

// sign · √square.
struct SignedSqrt {
  int sign = 0;     // −1, 0 or +1
  Rational square;  // ≥ 0; zero exactly when sign is 0

  bool is_zero() const { return sign == 0; }
  double value() const { return sign * std::sqrt(square.to_double()); }
  std::string str() const {
    if (sign == 0) return "0";
    const auto r = rational_sqrt(square);
    const std::string body = r ? r->str() : "sqrt(" + square.str() + ")";
    return sign < 0 ? "-" + body : body;
  }
  friend bool operator==(const SignedSqrt& a, const SignedSqrt& b) { return a.sign == b.sign && a.square == b.square; }
  friend SignedSqrt operator*(const SignedSqrt& a, const SignedSqrt& b) {
    if (a.sign == 0 || b.sign == 0) return {};
    return {a.sign * b.sign, a.square * b.square};
  }
  friend SignedSqrt operator*(const Rational& r, const SignedSqrt& a) {
    if (r.is_zero() || a.sign == 0) return {};
    return {r.sign() * a.sign, r * r * a.square};
  }
  SignedSqrt operator-() const { return {-sign, square}; }
};

// r as a surd.
inline SignedSqrt surd(const Rational& r) { return r.is_zero() ? SignedSqrt{} : SignedSqrt{r.sign(), r * r}; }

// The exact sum of surds whose squares have rational-square ratios; throws std::domain_error
// for incommensurable ones (such as √2 + √3).
inline SignedSqrt sum(const std::vector<SignedSqrt>& terms) {
  const SignedSqrt* base = nullptr;
  for (const auto& t : terms)
    if (t.sign) { base = &t; break; }
  if (!base) return {};
  Rational total = 0;
  for (const auto& t : terms) {
    if (!t.sign) continue;
    const auto r = rational_sqrt(t.square / base->square);
    if (!r) throw std::domain_error("physis::symmetry::sum: incommensurable surds " + t.str() + " and " + base->str());
    total += t.sign > 0 ? *r : -*r;
  }
  if (total.is_zero()) return {};
  return {total.sign(), total * total * base->square};
}

// n!, from a table kept per thread. A deque, so that a reference returned earlier stays valid
// while the table grows within one expression.
inline const Integer& factorial(int n) {
  if (n < 0) throw std::domain_error("physis::symmetry::factorial: negative argument");
  thread_local std::deque<Integer> table{Integer(1)};
  while (static_cast<int>(table.size()) <= n) table.push_back(table.back() * Integer(static_cast<long long>(table.size())));
  return table[n];
}

// Whether spins (as two_j) a, b, c satisfy the triangle conditions: |a − b| ≤ c ≤ a + b and
// a + b + c even (so that j₁ + j₂ + j₃ is an integer).
inline bool triangle(int a, int b, int c) {
  return a >= 0 && b >= 0 && c >= 0 && (a + b + c) % 2 == 0 && c >= std::abs(a - b) && c <= a + b;
}

// Δ(abc)² for a triad satisfying the triangle conditions.
inline Rational delta_squared(int a, int b, int c) {
  return Rational(factorial((a + b - c) / 2) * factorial((a - b + c) / 2) * factorial((-a + b + c) / 2), factorial((a + b + c) / 2 + 1));
}

namespace detail {

inline bool parity_ok(int two_j, int two_m) { return std::abs(two_m) <= two_j && (two_j - two_m) % 2 == 0; }

// The sum in Racah's formula for {a b c; d e f} (arguments as two_j), assuming its four triads
// are admissible.
inline Rational racah_sum(int a, int b, int c, int d, int e, int f) {
  const int al[4] = {(a + b + c) / 2, (a + e + f) / 2, (d + b + f) / 2, (d + e + c) / 2};
  const int be[3] = {(a + b + d + e) / 2, (b + c + e + f) / 2, (c + a + f + d) / 2};
  const int lo = *std::max_element(al, al + 4), hi = *std::min_element(be, be + 3);
  Rational s = 0;
  for (int z = lo; z <= hi; z++) {
    Integer den = 1;
    for (int i = 0; i < 4; i++) den *= factorial(z - al[i]);
    for (int k = 0; k < 3; k++) den *= factorial(be[k] - z);
    const Rational term(factorial(z + 1), den);
    s += (z % 2) ? -term : term;
  }
  return s;
}

}  // namespace detail

// {a b c; d e f}, arguments as two_j.
inline SignedSqrt wigner_6j(int a, int b, int c, int d, int e, int f) {
  if (!triangle(a, b, c) || !triangle(a, e, f) || !triangle(d, b, f) || !triangle(d, e, c)) return {};
  const Rational s = detail::racah_sum(a, b, c, d, e, f);
  if (s.is_zero()) return {};
  return {s.sign(), s * s * delta_squared(a, b, c) * delta_squared(a, e, f) * delta_squared(d, b, f) * delta_squared(d, e, c)};
}

// ⟨j1 m1 j2 m2 | J M⟩, arguments as twice the values.
inline SignedSqrt clebsch_gordan(int j1, int m1, int j2, int m2, int J, int M) {
  using detail::parity_ok;
  if (M != m1 + m2 || !triangle(j1, j2, J) || !parity_ok(j1, m1) || !parity_ok(j2, m2) || !parity_ok(J, M)) return {};
  // All in units of one: the arguments of the factorials.
  const int a = (j1 + j2 - J) / 2, b = (j1 - m1) / 2, c = (j2 + m2) / 2, d = (J - j2 + m1) / 2, e = (J - j1 - m2) / 2;
  const int lo = std::max({0, -d, -e}), hi = std::min({a, b, c});
  Rational s = 0;
  for (int k = lo; k <= hi; k++) {
    const Rational term(1, factorial(k) * factorial(a - k) * factorial(b - k) * factorial(c - k) * factorial(d + k) * factorial(e + k));
    s += (k % 2) ? -term : term;
  }
  if (s.is_zero()) return {};
  const Rational pref = Rational(Integer(J + 1) * factorial((J + j1 - j2) / 2) * factorial((J - j1 + j2) / 2) * factorial((j1 + j2 - J) / 2),
                                 factorial((j1 + j2 + J) / 2 + 1)) *
                        Rational(factorial((J + M) / 2) * factorial((J - M) / 2) * factorial((j1 - m1) / 2) * factorial((j1 + m1) / 2) *
                                 factorial((j2 - m2) / 2) * factorial((j2 + m2) / 2));
  return {s.sign(), s * s * pref};
}

// (j1 j2 j3; m1 m2 m3), arguments as twice the values.
inline SignedSqrt wigner_3j(int j1, int j2, int j3, int m1, int m2, int m3) {
  if (m1 + m2 + m3 != 0) return {};
  const SignedSqrt cg = clebsch_gordan(j1, m1, j2, m2, j3, -m3);
  if (cg.is_zero()) return {};
  const int phase = (j1 - j2 - m3) / 2;
  return {(phase % 2 == 0 ? 1 : -1) * cg.sign, cg.square / Rational(j3 + 1)};
}

// {j1 j2 j3; j4 j5 j6; j7 j8 j9}, arguments as two_j.
inline SignedSqrt wigner_9j(int j1, int j2, int j3, int j4, int j5, int j6, int j7, int j8, int j9) {
  if (!triangle(j1, j2, j3) || !triangle(j4, j5, j6) || !triangle(j7, j8, j9) || !triangle(j1, j4, j7) || !triangle(j2, j5, j8) ||
      !triangle(j3, j6, j9))
    return {};
  const int lo = std::max({std::abs(j1 - j9), std::abs(j4 - j8), std::abs(j2 - j6)});
  const int hi = std::min({j1 + j9, j4 + j8, j2 + j6});
  Rational t = 0;
  for (int x = lo; x <= hi; x++) {
    if (!triangle(j1, j9, x) || !triangle(j8, j4, x) || !triangle(j2, j6, x)) continue;
    const Rational term = Rational(x + 1) * detail::racah_sum(j1, j4, j7, j8, j9, x) * detail::racah_sum(j2, j5, j8, j4, x, j6) *
                          detail::racah_sum(j3, j6, j9, x, j1, j2) * delta_squared(j1, j9, x) * delta_squared(j8, j4, x) *
                          delta_squared(j2, j6, x);
    t += (x % 2) ? -term : term;
  }
  if (t.is_zero()) return {};
  const Rational D = delta_squared(j1, j2, j3) * delta_squared(j4, j5, j6) * delta_squared(j7, j8, j9) * delta_squared(j1, j4, j7) *
                     delta_squared(j2, j5, j8) * delta_squared(j3, j6, j9);
  return {t.sign(), t * t * D};
}

}  // namespace physis::symmetry
