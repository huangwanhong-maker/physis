// The group SU(2) of unit quaternions: spin, isospin, and the gauge group of the weak force
// and of three-dimensional gravity.
//
// An element is U = a0 I + i (a1 σ1 + a2 σ2 + a3 σ3) with a0² + |a|² = 1, the Pauli matrices
// σ_k; so SU(2) is the 3-sphere S³ ⊂ ℝ⁴ and the product is that of quaternions:
//   (a0, a)(b0, b) = (a0 b0 − a·b, a0 b + b0 a − a × b).
// Real combinations of elements (the staples of lattice gauge theory) are quaternions, real
// multiples of elements, so the struct also carries the arithmetic of ℍ.
//
// - **The algebra** su(2) has the basis T_k = −(i/2) σ_k, with [T_j, T_k] = ε_jkl T_l; the
//   exponential of x·T is the rotation by the angle |x| about x / |x|,
//   exp(x·T) = cos(|x|/2) I − i sin(|x|/2) (x̂·σ).
// - **The double cover** SU(2) → SO(3): R(U)_ij = ½ tr(σ_i U σ_j U†). With the convention
//   above, U = exp(θ n·T) maps to the active (right-handed) rotation by θ about n, and U and
//   −U map to the same rotation (Sakurai and Napolitano 2017, Modern Quantum Mechanics,
//   §3.3).
// - **Classes and characters.** Every element is conjugate to diag(e^{iθ/2}, e^{−iθ/2}),
//   θ ∈ [0, 2π], with a0 = cos(θ/2): its class angle, the angle of the rotation it covers.
//   The irreducible representations are the spins j ∈ ½ℕ, of dimension 2j + 1, with
//   characters χ_j(θ) = sin((2j + 1)θ/2) / sin(θ/2) = U_{2j}(cos(θ/2)), the Chebyshev
//   polynomial of the second kind, evaluated here by its recurrence so that θ = 0 and 2π
//   need no limit (Bröcker and tom Dieck 1985, II.5). Spins are passed as two_j = 2j.
// - **Haar measure.** The normalised volume of S³: a uniform point is a normalised Gaussian
//   4-vector (Muller 1959, Comm. ACM 2(4), 19–20).

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <numbers>

#include <calx/random/stream.hpp>

namespace physis::symmetry {

using Complex = std::complex<double>;
// A 2×2 complex matrix, row-major: m[i][j].
using Matrix2c = std::array<std::array<Complex, 2>, 2>;
// A 3×3 real matrix, row-major.
using Matrix3 = std::array<std::array<double, 3>, 3>;
using Vector3 = std::array<double, 3>;

// The Pauli matrices σ_1, σ_2, σ_3.
inline const std::array<Matrix2c, 3>& pauli() {
  static const std::array<Matrix2c, 3> s{{
      {{{Complex(0, 0), Complex(1, 0)}, {Complex(1, 0), Complex(0, 0)}}},
      {{{Complex(0, 0), Complex(0, -1)}, {Complex(0, 1), Complex(0, 0)}}},
      {{{Complex(1, 0), Complex(0, 0)}, {Complex(0, 0), Complex(-1, 0)}}},
  }};
  return s;
}

struct SU2 {
  double a0 = 1, a1 = 0, a2 = 0, a3 = 0;  // U = a0 I + i a·σ

  static SU2 identity() { return {}; }
  friend SU2 operator*(const SU2& x, const SU2& y) {
    return {x.a0 * y.a0 - x.a1 * y.a1 - x.a2 * y.a2 - x.a3 * y.a3,
            x.a0 * y.a1 + y.a0 * x.a1 - (x.a2 * y.a3 - x.a3 * y.a2),
            x.a0 * y.a2 + y.a0 * x.a2 - (x.a3 * y.a1 - x.a1 * y.a3),
            x.a0 * y.a3 + y.a0 * x.a3 - (x.a1 * y.a2 - x.a2 * y.a1)};
  }
  // U⁻¹ = U†, the conjugate quaternion (for a unit quaternion).
  SU2 inverse() const { return {a0, -a1, -a2, -a3}; }
  SU2 adjoint() const { return inverse(); }
  double trace() const { return 2 * a0; }
  // Re tr U / N for the defining representation (N = 2).
  double re_trace() const { return a0; }
  // The class angle θ ∈ [0, 2π], a0 = cos(θ/2).
  double class_angle() const { return 2 * std::acos(std::clamp(a0 / norm(), -1.0, 1.0)); }

  // The quaternion algebra ℍ.
  friend SU2 operator+(const SU2& x, const SU2& y) { return {x.a0 + y.a0, x.a1 + y.a1, x.a2 + y.a2, x.a3 + y.a3}; }
  friend SU2 operator-(const SU2& x, const SU2& y) { return {x.a0 - y.a0, x.a1 - y.a1, x.a2 - y.a2, x.a3 - y.a3}; }
  friend SU2 operator*(double s, const SU2& x) { return {s * x.a0, s * x.a1, s * x.a2, s * x.a3}; }
  SU2& operator+=(const SU2& y) { return *this = *this + y; }
  // The quaternion norm; 1 on SU(2), and |V| = √det V for a real multiple V of an element.
  double norm() const { return std::sqrt(a0 * a0 + a1 * a1 + a2 * a2 + a3 * a3); }
  // The nearest element of SU(2): the quaternion divided by its norm.
  SU2 reprojected() const {
    const double n = norm();
    return {a0 / n, a1 / n, a2 / n, a3 / n};
  }

  Matrix2c matrix() const {
    return {{{Complex(a0, a3), Complex(a2, a1)}, {Complex(-a2, a1), Complex(a0, -a3)}}};
  }
};

// The element of a 2×2 matrix of the form a0 I + i a·σ (a0 = Re tr / 2, a_k = Im tr(σ_k M) / 2).
inline SU2 from_matrix(const Matrix2c& m) {
  return {0.5 * (m[0][0] + m[1][1]).real(), 0.5 * (m[0][1] + m[1][0]).imag(), 0.5 * (m[0][1] - m[1][0]).real(),
          0.5 * (m[0][0] - m[1][1]).imag()};
}

// exp(x·T), T_k = −(i/2) σ_k: the rotation by |x| about x / |x|.
inline SU2 exp_su2(const Vector3& x) {
  const double t = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
  // sin(t/2)/t, by its series where t is small.
  const double s = t < 1e-4 ? 0.5 - t * t / 48 : std::sin(t / 2) / t;
  return {std::cos(t / 2), -s * x[0], -s * x[1], -s * x[2]};
}
// The rotation by angle θ about the unit axis n.
inline SU2 rotation_su2(const Vector3& n, double theta) { return exp_su2({theta * n[0], theta * n[1], theta * n[2]}); }

// The principal logarithm: x with exp_su2(x) = U and |x| = θ ∈ [0, 2π]. At U = −I (θ = 2π)
// the axis is undetermined and (2π, 0, 0) is returned.
inline Vector3 log_su2(const SU2& u) {
  const double v = std::sqrt(u.a1 * u.a1 + u.a2 * u.a2 + u.a3 * u.a3);
  const double t = 2 * std::atan2(v, u.a0);
  if (v == 0) return {u.a0 >= 0 ? 0.0 : 2 * std::numbers::pi, 0.0, 0.0};
  const double f = -t / v;
  return {f * u.a1, f * u.a2, f * u.a3};
}

// The rotation R(U)_ij = ½ tr(σ_i U σ_j U†) of ℝ³ covered by U (Rodrigues' form of it).
inline Matrix3 adjoint(const SU2& u) {
  const double w = u.a0, x = -u.a1, y = -u.a2, z = -u.a3;  // the rotation quaternion w + x i + y j + z k
  return {{{1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)},
           {2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)},
           {2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)}}};
}

// χ_j at the element whose a0 = cos(θ/2) is given: U_{2j}(a0), Chebyshev's recurrence
// U_{n+1} = 2 a0 U_n − U_{n−1}.
inline double character_at(int two_j, double a0) {
  double u0 = 1, u1 = 2 * a0;
  if (two_j == 0) return u0;
  for (int n = 1; n < two_j; n++) {
    const double u2 = 2 * a0 * u1 - u0;
    u0 = u1;
    u1 = u2;
  }
  return u1;
}
// χ_j(U) for U ∈ SU(2).
inline double character(int two_j, const SU2& u) { return character_at(two_j, u.a0); }
// χ_j at the class angle θ.
inline double character_su2(int two_j, double theta) { return character_at(two_j, std::cos(theta / 2)); }

// A draw from the Haar measure: a normalised Gaussian 4-vector.
inline SU2 haar_su2(calx::random::Stream& z) {
  for (;;) {
    const SU2 q{z.normal(), z.normal(), z.normal(), z.normal()};
    if (q.norm() > 1e-12) return q.reprojected();
  }
}

}  // namespace physis::symmetry
