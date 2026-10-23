// The circle group U(1) = {e^{iθ}}: the gauge group of electromagnetism and the simplest
// compact Lie group.
//
// An element is held by its angle θ, kept in (−π, π]; the product adds angles. Its
// irreducible representations are one-dimensional, labelled by n ∈ ℤ, with characters
// χ_n(θ) = e^{inθ}; they are orthonormal for the Haar measure dθ / 2π (Peter–Weyl, Bröcker
// and tom Dieck 1985, Representations of Compact Lie Groups, ch. II–III). The trace of the
// defining representation is e^{iθ}, whose real part cos θ enters Wilson's action.

#pragma once

#include <cmath>
#include <complex>
#include <numbers>

#include <calx/random/stream.hpp>

namespace physis::symmetry {

// θ reduced to (−π, π].
inline double wrap_angle(double theta) {
  constexpr double two_pi = 2 * std::numbers::pi;
  double t = std::remainder(theta, two_pi);
  if (t <= -std::numbers::pi) t += two_pi;
  return t;
}

struct U1 {
  double theta = 0;  // the angle, rad, in (−π, π]

  static U1 identity() { return {}; }
  static U1 from_angle(double t) { return {wrap_angle(t)}; }
  friend U1 operator*(const U1& a, const U1& b) { return {wrap_angle(a.theta + b.theta)}; }
  U1 inverse() const { return {wrap_angle(-theta)}; }
  std::complex<double> value() const { return std::polar(1.0, theta); }
  // Re tr U / N for the defining representation (N = 1).
  double re_trace() const { return std::cos(theta); }
};

// χ_n(U) = e^{inθ}.
inline std::complex<double> character(int n, const U1& u) { return std::polar(1.0, n * u.theta); }

// A draw from the Haar measure: θ uniform on (−π, π].
inline U1 haar_u1(calx::random::Stream& z) { return U1::from_angle(std::numbers::pi * (2 * z.uniform() - 1)); }

}  // namespace physis::symmetry
