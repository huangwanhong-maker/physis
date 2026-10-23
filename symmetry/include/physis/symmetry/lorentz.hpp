// The Lorentz group: the symmetries of Minkowski space that fix the origin.
//
// Coordinates are x = (x⁰, x¹, x², x³) with x⁰ = ct, and the metric is η = diag(−1, 1, 1, 1)
// (physis decision 0005, §2.4). A Lorentz transformation is a real 4×4 matrix Λ with
// Λᵀ η Λ = η, acting on contravariant components as x'^μ = Λ^μ_ν x^ν; the proper
// orthochronous ones (det Λ = 1, Λ⁰₀ ≥ 1) form SO⁺(1, 3), and Λ⁻¹ = η Λᵀ η.
//
// - **Boosts** of rapidity φ along the unit vector n: Λ⁰₀ = cosh φ, Λ⁰_i = Λ^i₀ = sinh φ n_i,
//   Λ^i_j = δ_ij + (cosh φ − 1) n_i n_j, which maps the rest frame to one moving with
//   velocity c tanh φ n. Rapidities of collinear boosts add.
// - **The universal cover** SL(2, ℂ) → SO⁺(1, 3). A point x is the Hermitian matrix
//   X = x^μ σ_μ (σ_0 = I, σ_k Pauli), with det X = −η(x, x); A ∈ SL(2, ℂ) acts by
//   X ↦ A X A†, which preserves det X, so that Λ(A)^μ_ν = ½ tr(σ_μ A σ_ν A†). The kernel is
//   {±I}. The boost above is covered by A = exp(φ n·σ / 2) = cosh(φ/2) I + sinh(φ/2) n·σ, and
//   the rotation by θ about n by the SU(2) element exp(−iθ n·σ / 2) (su2.hpp) (Penrose and
//   Rindler 1984, Spinors and Space-Time 1, §1.2; Sexl and Urbantke 2001, Relativity, Groups,
//   Particles, ch. 9).
// - **The Wigner rotation.** Two boosts in different directions compose to a boost and a
//   rotation, Λ(u₂) Λ(u₁) = Λ(u₃) R, the rotation behind Thomas precession. With
//   γ = γ₁ γ₂ (1 + u₁·u₂ / c²) the Lorentz factor of the composition, the angle ε of R obeys
//   cos ε = (1 + γ + γ₁ + γ₂)² / ((1 + γ)(1 + γ₁)(1 + γ₂)) − 1
//   (Sexl and Urbantke 2001, §4.3; O'Donnell and Visser 2011, Eur. J. Phys. 32, 1033).

#pragma once

#include <array>
#include <cmath>
#include <complex>

#include "su2.hpp"

namespace physis::symmetry {

using Matrix4 = std::array<std::array<double, 4>, 4>;

// η = diag(−1, 1, 1, 1).
inline Matrix4 minkowski_metric() {
  Matrix4 e{};
  e[0][0] = -1;
  e[1][1] = e[2][2] = e[3][3] = 1;
  return e;
}

inline Matrix4 multiply(const Matrix4& x, const Matrix4& y) {
  Matrix4 z{};
  for (int i = 0; i < 4; i++)
    for (int k = 0; k < 4; k++)
      for (int j = 0; j < 4; j++) z[i][j] += x[i][k] * y[k][j];
  return z;
}

// Λ⁻¹ = η Λᵀ η.
inline Matrix4 lorentz_inverse(const Matrix4& L) {
  Matrix4 m{};
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) m[i][j] = ((i == 0) == (j == 0) ? 1 : -1) * L[j][i];
  return m;
}

// max |(Λᵀ η Λ − η)_ij|: zero to rounding for a Lorentz transformation.
inline double lorentz_defect(const Matrix4& L) {
  const Matrix4 e = minkowski_metric();
  double d = 0;
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) {
      double s = 0;
      for (int k = 0; k < 4; k++) s += L[k][i] * e[k][k] * L[k][j];
      d = std::max(d, std::fabs(s - e[i][j]));
    }
  return d;
}

// The boost of rapidity phi along the unit vector n.
inline Matrix4 boost(double phi, const Vector3& n) {
  const double c = std::cosh(phi), s = std::sinh(phi);
  Matrix4 L{};
  L[0][0] = c;
  for (int i = 0; i < 3; i++) {
    L[0][i + 1] = L[i + 1][0] = s * n[i];
    for (int j = 0; j < 3; j++) L[i + 1][j + 1] = (i == j ? 1 : 0) + (c - 1) * n[i] * n[j];
  }
  return L;
}

// The rotation by theta about the unit vector n, acting on the spatial components.
inline Matrix4 rotation(const Vector3& n, double theta) {
  const Matrix3 R = adjoint(rotation_su2(n, theta));
  Matrix4 L{};
  L[0][0] = 1;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) L[i + 1][j + 1] = R[i][j];
  return L;
}

inline Matrix2c multiply(const Matrix2c& x, const Matrix2c& y) {
  Matrix2c z{};
  for (int i = 0; i < 2; i++)
    for (int k = 0; k < 2; k++)
      for (int j = 0; j < 2; j++) z[i][j] += x[i][k] * y[k][j];
  return z;
}
inline Matrix2c dagger(const Matrix2c& x) { return {{{std::conj(x[0][0]), std::conj(x[1][0])}, {std::conj(x[0][1]), std::conj(x[1][1])}}}; }

// Λ(A)^μ_ν = ½ tr(σ_μ A σ_ν A†), σ_0 = I.
inline Matrix4 cover(const Matrix2c& A) {
  std::array<Matrix2c, 4> s;
  s[0] = {{{Complex(1, 0), Complex(0, 0)}, {Complex(0, 0), Complex(1, 0)}}};
  for (int k = 0; k < 3; k++) s[k + 1] = pauli()[k];
  const Matrix2c Ad = dagger(A);
  Matrix4 L{};
  for (int nu = 0; nu < 4; nu++) {
    const Matrix2c X = multiply(multiply(A, s[nu]), Ad);
    for (int mu = 0; mu < 4; mu++) {
      const Matrix2c P = multiply(s[mu], X);
      L[mu][nu] = 0.5 * (P[0][0] + P[1][1]).real();
    }
  }
  return L;
}

// exp(φ n·σ / 2), which covers boost(φ, n).
inline Matrix2c boost_sl2c(double phi, const Vector3& n) {
  const double c = std::cosh(phi / 2), s = std::sinh(phi / 2);
  return {{{Complex(c + s * n[2], 0), Complex(s * n[0], -s * n[1])}, {Complex(s * n[0], s * n[1]), Complex(c - s * n[2], 0)}}};
}

// Λ = B R, a boost after a rotation (the polar decomposition of a proper orthochronous Λ):
// B carries the time axis where Λ does.
struct BoostAndRotation {
  double rapidity = 0;
  Vector3 direction{1, 0, 0};
  Matrix4 boost, rotation;
};
inline BoostAndRotation decompose(const Matrix4& L) {
  BoostAndRotation d;
  d.rapidity = std::acosh(std::max(1.0, L[0][0]));
  const double v = std::sqrt(L[1][0] * L[1][0] + L[2][0] * L[2][0] + L[3][0] * L[3][0]);
  if (v > 0) d.direction = {L[1][0] / v, L[2][0] / v, L[3][0] / v};
  d.boost = boost(d.rapidity, d.direction);
  d.rotation = multiply(lorentz_inverse(d.boost), L);
  return d;
}

// The angle of the spatial rotation in a Lorentz transformation of the form diag(1, R).
inline double rotation_angle(const Matrix4& R) {
  return std::acos(std::clamp((R[1][1] + R[2][2] + R[3][3] - 1) / 2, -1.0, 1.0));
}

// The Wigner angle of boost(phi2, n2) · boost(phi1, n1), by its closed form.
inline double wigner_angle(double phi1, const Vector3& n1, double phi2, const Vector3& n2) {
  const double g1 = std::cosh(phi1), g2 = std::cosh(phi2);
  const double dot = std::tanh(phi1) * std::tanh(phi2) * (n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2]);
  const double g = g1 * g2 * (1 + dot);
  const double c = (1 + g + g1 + g2) * (1 + g + g1 + g2) / ((1 + g) * (1 + g1) * (1 + g2)) - 1;
  return std::acos(std::clamp(c, -1.0, 1.0));
}

}  // namespace physis::symmetry
