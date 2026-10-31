// The memory left by eliminating variables from a linear system: the exact reduced equation of the
// variables kept, with its memory kernel and the term from the eliminated variables' initial state, and its
// limit when the eliminated variables are fast (decision 0018).
//
// The identity (Zwanzig 1961, Phys. Rev. 124, 983; Mori 1965, Prog. Theor. Phys. 33, 423; for a linear system
// it is variation of constants: Chorin, Hald and Kupferman 2000, PNAS 97, 2968). With x kept and y eliminated,
//
//   ẋ = A x + B y,  ẏ = C x + D y   ⟹   ẋ(t) = A x(t) + B e^{Dt} y(0) + ∫₀ᵗ K(t − s) x(s) ds,  K(τ) = B e^{Dτ} C,
//
// since y(t) = e^{Dt} y(0) + ∫₀ᵗ e^{D(t−s)} C x(s) ds. The reduced equation is exact and not Markovian: x's
// future depends on its past through K, and on y(0) through the second term, which in a statistical reading is
// the noise. When D is stable and fast the kernel concentrates at τ = 0 with weight ∫₀^∞ K = −B D⁻¹ C, and the
// Markovian limit ẋ ≈ (A − B D⁻¹ C) x is adiabatic elimination; the error is of the order of the ratio of the
// time scales.
//
// Dense blocks: A (m × m), B (m × k), C (k × m), D (k × k). Checks: tests/check_elimination.cpp.

#pragma once

#include <cstddef>
#include <stdexcept>
#include <vector>

#include <calx/linalg/dense.hpp>
#include <calx/linalg/expm.hpp>
#include <calx/linalg/rectangular.hpp>

namespace physis::dynamics {

struct LinearElimination {
  calx::linalg::Matrix A, D;
  calx::linalg::Rectangular B, C;

  void validate() const {
    if (B.rows != A.n || B.cols != D.n || C.rows != D.n || C.cols != A.n)
      throw std::invalid_argument("physis::dynamics::LinearElimination: blocks A (m×m), B (m×k), C (k×m), D (k×k)");
  }
  // K(τ) = B e^{Dτ} C, an m × m matrix.
  calx::linalg::Matrix kernel(double tau) const {
    validate();
    const auto E = calx::linalg::expm(tau * D);
    calx::linalg::Matrix K(A.n);
    for (std::size_t i = 0; i < A.n; i++)
      for (std::size_t a = 0; a < D.n; a++) {
        double be = 0;
        for (std::size_t b = 0; b < D.n; b++) be += B(i, b) * E(b, a);
        if (be == 0) continue;
        for (std::size_t j = 0; j < A.n; j++) K(i, j) += be * C(a, j);
      }
    return K;
  }
  // B e^{Dt} y0.
  std::vector<double> noise(double t, const std::vector<double>& y0) const {
    validate();
    const auto E = calx::linalg::expm(t * D);
    std::vector<double> ey(D.n, 0.0), out(A.n, 0.0);
    for (std::size_t a = 0; a < D.n; a++)
      for (std::size_t b = 0; b < D.n; b++) ey[a] += E(a, b) * y0[b];
    for (std::size_t i = 0; i < A.n; i++)
      for (std::size_t a = 0; a < D.n; a++) out[i] += B(i, a) * ey[a];
    return out;
  }
  // The full system's matrix [[A, B], [C, D]].
  calx::linalg::Matrix full() const {
    validate();
    const std::size_t m = A.n, k = D.n;
    calx::linalg::Matrix M(m + k);
    for (std::size_t i = 0; i < m; i++)
      for (std::size_t j = 0; j < m; j++) M(i, j) = A(i, j);
    for (std::size_t i = 0; i < m; i++)
      for (std::size_t a = 0; a < k; a++) M(i, m + a) = B(i, a);
    for (std::size_t a = 0; a < k; a++)
      for (std::size_t j = 0; j < m; j++) M(m + a, j) = C(a, j);
    for (std::size_t a = 0; a < k; a++)
      for (std::size_t b = 0; b < k; b++) M(m + a, m + b) = D(a, b);
    return M;
  }
  // The Markovian (adiabatic) limit A − B D⁻¹ C; D must be invertible.
  calx::linalg::Matrix markovian() const {
    validate();
    // Z = D⁻¹ C column by column, each column of C solved as the first column of a square right-hand side.
    calx::linalg::Matrix out = A;
    for (std::size_t j = 0; j < A.n; j++) {
      calx::linalg::Matrix rhs(D.n);
      for (std::size_t a = 0; a < D.n; a++) rhs(a, 0) = C(a, j);
      const auto Z = calx::linalg::solve(D, rhs);
      for (std::size_t i = 0; i < A.n; i++)
        for (std::size_t a = 0; a < D.n; a++) out(i, j) -= B(i, a) * Z(a, 0);
    }
    return out;
  }
};

}  // namespace physis::dynamics
