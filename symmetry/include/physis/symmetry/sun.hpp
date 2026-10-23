// The special unitary groups SU(N): the gauge group of the strong force (N = 3), and of
// lattice gauge theory in general.
//
// An element is an N×N complex matrix U with U U† = I and det U = 1; the order N is chosen
// at run time. The same type holds any complex matrix of order N, since sums of elements
// (staples) and the generators are needed beside the elements themselves.
//
// - **Generators.** The generalised Gell-Mann matrices λ_a, a = 1 … N² − 1: symmetric and
//   antisymmetric off-diagonal ones and N − 1 diagonal ones (Georgi 1999, Lie Algebras in
//   Particle Physics, §13; Bertlmann and Krammer 2008, J. Phys. A 41, 235303). The Hermitian
//   generators t_a = λ_a / 2 satisfy tr(t_a t_b) = δ_ab / 2, and exp(−i x^a t_a) ∈ SU(N); for
//   N = 2 this is physis::symmetry::exp_su2.
// - **The exponential** of an anti-Hermitian matrix X = −i x^a t_a is computed as that of its
//   real form [[Re X, −Im X], [Im X, Re X]] of order 2N by calx's Padé approximant, accurate
//   to about 10⁻¹⁴.
// - **Reunitarisation.** Rounding makes a long product drift off the group; the rows are
//   orthonormalised by modified Gram–Schmidt, and the phase of the determinant, e^{iφ}, is
//   removed by the factor e^{−iφ/N}.
// - **Haar measure.** The Gram–Schmidt orthonormalisation of a Ginibre matrix (independent
//   standard complex Gaussian entries) is Haar-distributed on U(N): Mezzadri's phase
//   correction (Mezzadri 2007, Notices AMS 54, 592–604) is built in, since Gram–Schmidt gives
//   the triangular factor a positive diagonal. Removing the phase of the determinant as above
//   then gives the Haar measure on SU(N), because U(N) = (SU(N) × U(1)) / ℤ_N and the centre
//   preserves the measure.
// - **Determinant** by Gaussian elimination with partial pivoting, O(N³).

#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include <calx/linalg/expm.hpp>
#include <calx/random/stream.hpp>

namespace physis::symmetry {

struct SUN {
  using Complex = std::complex<double>;

  std::size_t n = 0;
  std::vector<Complex> a;  // row-major

  SUN() = default;
  explicit SUN(std::size_t order) : n(order), a(order * order) {}
  static SUN identity(std::size_t order) {
    SUN m(order);
    for (std::size_t i = 0; i < order; i++) m(i, i) = 1;
    return m;
  }
  Complex& operator()(std::size_t i, std::size_t j) { return a[i * n + j]; }
  const Complex& operator()(std::size_t i, std::size_t j) const { return a[i * n + j]; }

  friend SUN operator*(const SUN& x, const SUN& y) {
    SUN z(x.n);
    for (std::size_t i = 0; i < x.n; i++)
      for (std::size_t k = 0; k < x.n; k++) {
        const Complex xik = x(i, k);
        for (std::size_t j = 0; j < x.n; j++) z(i, j) += xik * y(k, j);
      }
    return z;
  }
  friend SUN operator+(SUN x, const SUN& y) {
    for (std::size_t i = 0; i < x.a.size(); i++) x.a[i] += y.a[i];
    return x;
  }
  friend SUN operator-(SUN x, const SUN& y) {
    for (std::size_t i = 0; i < x.a.size(); i++) x.a[i] -= y.a[i];
    return x;
  }
  friend SUN operator*(Complex s, SUN x) {
    for (auto& v : x.a) v *= s;
    return x;
  }
  SUN& operator+=(const SUN& y) { return *this = *this + y; }

  SUN adjoint() const {
    SUN m(n);
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++) m(i, j) = std::conj((*this)(j, i));
    return m;
  }
  // U⁻¹ = U† on the group.
  SUN inverse() const { return adjoint(); }
  Complex trace() const {
    Complex t = 0;
    for (std::size_t i = 0; i < n; i++) t += (*this)(i, i);
    return t;
  }
  // Re tr U / N.
  double re_trace() const { return trace().real() / static_cast<double>(n); }

  Complex determinant() const {
    std::vector<Complex> m = a;
    Complex det = 1;
    for (std::size_t k = 0; k < n; k++) {
      std::size_t p = k;
      for (std::size_t i = k + 1; i < n; i++)
        if (std::abs(m[i * n + k]) > std::abs(m[p * n + k])) p = i;
      if (m[p * n + k] == Complex(0)) return 0;
      if (p != k) {
        for (std::size_t j = 0; j < n; j++) std::swap(m[k * n + j], m[p * n + j]);
        det = -det;
      }
      det *= m[k * n + k];
      for (std::size_t i = k + 1; i < n; i++) {
        const Complex f = m[i * n + k] / m[k * n + k];
        for (std::size_t j = k; j < n; j++) m[i * n + j] -= f * m[k * n + j];
      }
    }
    return det;
  }

  // The nearest-by-construction element of SU(N): rows orthonormalised, determinant's phase removed.
  SUN reunitarised() const {
    SUN u = *this;
    for (std::size_t i = 0; i < n; i++) {
      for (std::size_t k = 0; k < i; k++) {
        Complex dot = 0;  // <row k, row i>
        for (std::size_t j = 0; j < n; j++) dot += std::conj(u(k, j)) * u(i, j);
        for (std::size_t j = 0; j < n; j++) u(i, j) -= dot * u(k, j);
      }
      double norm = 0;
      for (std::size_t j = 0; j < n; j++) norm += std::norm(u(i, j));
      norm = std::sqrt(norm);
      if (norm == 0) throw std::domain_error("physis::symmetry::SUN::reunitarised: rows are linearly dependent");
      for (std::size_t j = 0; j < n; j++) u(i, j) /= norm;
    }
    const Complex phase = std::polar(1.0, -std::arg(u.determinant()) / static_cast<double>(n));
    return phase * u;
  }
  // max |(U U†)_ij − δ_ij|.
  double unitarity_defect() const {
    const SUN p = (*this) * adjoint();
    double d = 0;
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++) d = std::max(d, std::abs(p(i, j) - Complex(i == j ? 1.0 : 0.0)));
    return d;
  }
};

// The Hermitian generators t_a = λ_a / 2 of su(N), a = 0 … N² − 2: first the symmetric and
// antisymmetric pairs (j < k) in the order of (j, k), then the diagonal ones.
inline std::vector<SUN> generators(std::size_t n) {
  using C = SUN::Complex;
  std::vector<SUN> t;
  for (std::size_t j = 0; j < n; j++)
    for (std::size_t k = j + 1; k < n; k++) {
      SUN s(n), a(n);
      s(j, k) = s(k, j) = 0.5;
      a(j, k) = C(0, -0.5);
      a(k, j) = C(0, 0.5);
      t.push_back(s);
      t.push_back(a);
    }
  for (std::size_t l = 1; l < n; l++) {
    SUN d(n);
    const double c = std::sqrt(2.0 / (static_cast<double>(l) * (l + 1))) / 2;
    for (std::size_t j = 0; j < l; j++) d(j, j) = c;
    d(l, l) = -c * static_cast<double>(l);
    t.push_back(d);
  }
  return t;
}

// exp(X) of any complex matrix, through its real form.
inline SUN exp_matrix(const SUN& X) {
  const std::size_t n = X.n;
  calx::linalg::Matrix R(2 * n);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) {
      R(i, j) = R(i + n, j + n) = X(i, j).real();
      R(i + n, j) = X(i, j).imag();
      R(i, j + n) = -X(i, j).imag();
    }
  const calx::linalg::Matrix E = calx::linalg::expm(R);
  SUN U(n);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) U(i, j) = SUN::Complex(E(i, j), E(i + n, j));
  return U;
}

// exp(−i x^a t_a) for coordinates x of length N² − 1.
inline SUN exp_sun(std::size_t n, const std::vector<double>& x) {
  const auto t = generators(n);
  if (x.size() != t.size()) throw std::invalid_argument("physis::symmetry::exp_sun: N^2 - 1 coordinates are needed");
  SUN X(n);
  for (std::size_t a = 0; a < t.size(); a++) X += SUN::Complex(0, -x[a]) * t[a];
  return exp_matrix(X);
}

// A draw from the Haar measure on SU(N).
inline SUN haar_sun(std::size_t n, calx::random::Stream& z) {
  SUN g(n);
  for (auto& v : g.a) v = SUN::Complex(z.normal(), z.normal()) / std::sqrt(2.0);
  return g.reunitarised();
}

}  // namespace physis::symmetry
