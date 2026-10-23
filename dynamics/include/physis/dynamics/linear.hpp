// Linear time-invariant systems x' = A x + B u + c, and their exact propagators.
//
// Over a step h with the input u held constant, the solution is
//
//   x(t + h) = e^{A h} x(t) + (int_0^h e^{A s} ds) (B u + c),
//
// and both terms come from one exponential of the augmented matrix
//
//   M = [[A, B, c], [0, 0, 0], [0, 0, 0]] h,   e^M = [[P, Q, r], [0, I, 0], [0, 0, 1]],
//
// whose upper blocks are P = e^{A h}, Q = int_0^h e^{A s} ds B and r = int_0^h e^{A s} ds c
// (Van Loan 1978, "Computing integrals involving the matrix exponential", IEEE Trans.
// Autom. Control 23, 395-404, Theorem 1). A need not be invertible. The propagator is exact
// up to the accuracy of the exponential (calx/linalg, about 1e-15 relative); applied step
// by step it reproduces the solution at the grid points without error growing with the
// number of steps beyond rounding (Rotter and Diesmann 1999, "Exact digital simulation of
// time-invariant linear systems with applications to neuronal modeling", Biol. Cybern. 81,
// 381-402).
//
// Held components: a component may be held still (its row of A, B and c set to zero), as
// a neuron's potential is clamped while it is refractory and the other components flow on.
//
// Scales: in SI units the components of one state can differ by ten orders of magnitude (a
// potential of 1e-2 V driven by a current of 1e-10 A through 1e8 ohm), and the augmented
// matrix then has entries of 1e5 beside entries of 1e-2: its norm forces many squarings,
// each of which loses accuracy in the small entries. The caller may give each component's
// and each input's typical size; the matrix is made dimensionless by the diagonal
// similarity D^-1 M D with those sizes rounded to powers of two (so the scaling itself is
// exact), exponentiated, and scaled back: exp(M) = D exp(D^-1 M D) D^-1.

#pragma once

#include <array>
#include <cmath>
#include <cstddef>

#include <calx/linalg/dense.hpp>
#include <calx/linalg/expm.hpp>

namespace physis::dynamics {

using calx::linalg::Matrix;

// x <- P x + Q u + r, for a state of N components and M inputs.
template <std::size_t N, std::size_t M = 1>
struct Propagator {
  std::array<double, N * N> P{};
  std::array<double, N * M> Q{};
  std::array<double, N> r{};

  void apply(std::array<double, N>& x, const std::array<double, M>& u) const {
    std::array<double, N> y{};
    for (std::size_t i = 0; i < N; i++) {
      double s = r[i];
      for (std::size_t j = 0; j < N; j++) s += P[i * N + j] * x[j];
      for (std::size_t k = 0; k < M; k++) s += Q[i * M + k] * u[k];
      y[i] = s;
    }
    x = y;
  }
  // The single-input case.
  void apply(std::array<double, N>& x, double u) const requires(M == 1) { apply(x, std::array<double, 1>{u}); }
};

// The nearest power of two to a positive size (1 for anything else).
inline double power_of_two(double size) {
  return size > 0 && std::isfinite(size) ? std::ldexp(1.0, static_cast<int>(std::lround(std::log2(size)))) : 1.0;
}

// A (N x N), B (N x M, row-major), c (N), the step h, the components held still, and the
// typical sizes of the components and of the inputs (all 1 if not given).
template <std::size_t N, std::size_t M = 1>
Propagator<N, M> propagator(const Matrix& A, const std::array<double, N * M>& B, const std::array<double, N>& c, double h,
                            const std::array<bool, N>& held = {}, const std::array<double, N>& xsize = {},
                            const std::array<double, M>& usize = {}) {
  constexpr std::size_t K = N + M + 1;
  std::array<double, K> d{};
  for (std::size_t i = 0; i < N; i++) d[i] = power_of_two(xsize[i]);
  for (std::size_t k = 0; k < M; k++) d[N + k] = power_of_two(usize[k]);
  d[N + M] = 1;
  // X = D^-1 M D: X(i, j) = M(i, j) d_j / d_i.
  Matrix X(K);
  for (std::size_t i = 0; i < N; i++) {
    if (held[i]) continue;
    for (std::size_t j = 0; j < N; j++) X(i, j) = A(i, j) * h * d[j] / d[i];
    for (std::size_t k = 0; k < M; k++) X(i, N + k) = B[i * M + k] * h * d[N + k] / d[i];
    X(i, N + M) = c[i] * h / d[i];
  }
  const Matrix E = calx::linalg::expm(X);
  Propagator<N, M> p;
  for (std::size_t i = 0; i < N; i++) {
    for (std::size_t j = 0; j < N; j++) p.P[i * N + j] = E(i, j) * d[i] / d[j];
    for (std::size_t k = 0; k < M; k++) p.Q[i * M + k] = E(i, N + k) * d[i] / d[N + k];
    p.r[i] = E(i, N + M) * d[i];
  }
  return p;
}

}  // namespace physis::dynamics
