// Integrals of class functions over compact groups, by the Weyl integration formula, and the
// character expansion of a class function.
//
// A class function f (f(g U g⁻¹) = f(U)) depends only on the eigenvalues, and its Haar
// integral reduces to one over the maximal torus with the Weyl density |Δ|², Δ the Vandermonde
// determinant of the eigenvalues (Weyl 1926; Bröcker and tom Dieck 1985, IV.1; Bump 2013, Lie
// Groups, ch. 17):
//
// - U(1): ∫ f dU = (1/2π) ∫_{−π}^{π} f(θ) dθ;
// - SU(2), f a function of the class angle θ ∈ [0, 2π] (eigenvalues e^{±iθ/2}):
//   ∫ f dU = (1/π) ∫_0^{2π} sin²(θ/2) f(θ) dθ;
// - SU(3), f a function of two eigenphases (θ₁, θ₂), the third being −θ₁ − θ₂:
//   ∫ f dU = (1/(3! (2π)²)) ∫∫ Π_{i<k} |e^{iθ_i} − e^{iθ_k}|² f dθ₁ dθ₂.
//
// The character expansion of a class function, f = Σ_R c_R χ_R, has coefficients
// c_R = ∫ f χ_R* dU, by the orthonormality of characters (Peter–Weyl). The Boltzmann weight of
// Wilson's lattice action is expanded so: for U(1), e^{β cos θ} = Σ_n I_n(β) e^{inθ}; for SU(2),
// e^{(β/2) tr U} = (2/β) Σ_j (2j + 1) I_{2j+1}(β) χ_j(U), with the modified Bessel functions I_n
// (Drouffe and Zuber 1983, Phys. Rep. 102, 1–119, §2.1). These are the inputs of exact
// two-dimensional Yang–Mills (physis gauge).
//
// The integrals are computed by calx's adaptive Gauss–Kronrod rule to the requested
// tolerances (the integrands are smooth and periodic).

#pragma once

#include <cmath>
#include <complex>
#include <numbers>

#include <calx/quadrature/gauss_kronrod.hpp>

#include "su2.hpp"

namespace physis::symmetry {

struct Tolerance {
  double absolute = 1e-14, relative = 1e-13;
};

// ∫_{U(1)} f dU for f(θ).
template <class F>
double integrate_u1(F&& f, Tolerance tol = {}) {
  constexpr double pi = std::numbers::pi;
  return calx::quadrature::integrate(f, -pi, pi, tol.absolute, tol.relative, 4000).value / (2 * pi);
}

// ∫_{SU(2)} f dU for f(θ) of the class angle.
template <class F>
double integrate_su2(F&& f, Tolerance tol = {}) {
  constexpr double pi = std::numbers::pi;
  auto g = [&](double t) {
    const double s = std::sin(t / 2);
    return s * s * f(t);
  };
  return calx::quadrature::integrate(g, 0.0, 2 * pi, tol.absolute, tol.relative, 4000).value / pi;
}

// ∫_{SU(3)} f dU for f(θ₁, θ₂) of two eigenphases.
template <class F>
double integrate_su3(F&& f, Tolerance tol = {1e-12, 1e-11}) {
  constexpr double pi = std::numbers::pi;
  auto inner = [&](double t1) {
    auto g = [&](double t2) {
      const double t3 = -t1 - t2;
      const double d12 = 2 * std::sin((t1 - t2) / 2), d13 = 2 * std::sin((t1 - t3) / 2), d23 = 2 * std::sin((t2 - t3) / 2);
      return d12 * d12 * d13 * d13 * d23 * d23 * f(t1, t2);
    };
    return calx::quadrature::integrate(g, -pi, pi, tol.absolute, tol.relative, 4000).value;
  };
  return calx::quadrature::integrate(inner, -pi, pi, tol.absolute, tol.relative, 4000).value / (6 * 4 * pi * pi);
}

// c_n = ∫ f χ_n* dU on U(1), for f(θ).
template <class F>
std::complex<double> coefficient_u1(int n, F&& f, Tolerance tol = {}) {
  const double re = integrate_u1([&](double t) { return f(t) * std::cos(n * t); }, tol);
  const double im = integrate_u1([&](double t) { return -f(t) * std::sin(n * t); }, tol);
  return {re, im};
}

// c_j = ∫ f χ_j dU on SU(2) (its characters are real), for f(θ) of the class angle.
template <class F>
double coefficient_su2(int two_j, F&& f, Tolerance tol = {}) {
  return integrate_su2([&](double t) { return f(t) * character_su2(two_j, t); }, tol);
}

}  // namespace physis::symmetry
