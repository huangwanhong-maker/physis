// Metric fields in a chart: the pseudo-Riemannian structure of a manifold, as code.
//
// A metric field on an n-dimensional manifold is a smooth assignment of a symmetric,
// nondegenerate bilinear form g_x to each tangent space (O'Neill 1983, Semi-Riemannian
// Geometry, ch. 3; Wald 1984, General Relativity, §2.3). In a chart with coordinates x^a it
// is the symmetric matrix of components g_ab(x) = g(∂_a, ∂_b). Its signature (the number of
// negative directions) is constant on a connected manifold: 0 for a Riemannian metric, 1
// for a Lorentzian one (physis decision 0005, §2.4: signature (-, +, +, +), the timelike
// coordinate first).
//
// Here a metric field is an object with
//   static constexpr std::size_t n;                                     // the dimension
//   template <class T> std::array<std::array<T, n>, n> operator()(const std::array<T, n>& x) const;
// its parameters (a mass, a radius, a scale factor) being members. Written once as a
// template over the number type T, it is evaluated on doubles for its components and on
// calx's jets for their first and second derivatives (curvature.hpp). The code of a metric
// calls the elementary functions unqualified, after `using std::sin;` and the like.
//
// The chart is part of the metric's type: the same manifold in two charts is two metric
// fields, related by an explicit change of coordinates (physis invariant 4).

#pragma once

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <stdexcept>
#include <utility>

namespace physis::geometry {

template <class T, std::size_t n>
using Matrix = std::array<std::array<T, n>, n>;
template <std::size_t n>
using Point = std::array<double, n>;
template <std::size_t n>
using Tensor2 = std::array<std::array<double, n>, n>;
template <std::size_t n>
using Tensor3 = std::array<Tensor2<n>, n>;
template <std::size_t n>
using Tensor4 = std::array<Tensor3<n>, n>;

template <class M>
concept MetricField = requires(const M& g, const Point<M::n>& x) {
  { M::n } -> std::convertible_to<std::size_t>;
  { g(x) } -> std::convertible_to<Tensor2<M::n>>;
};

// The determinant and inverse of a small matrix, by Gauss-Jordan elimination with partial
// pivoting (Golub and Van Loan 2013, §3.4); throws when the matrix is singular.
template <std::size_t n>
inline double determinant(Tensor2<n> a) {
  double det = 1;
  for (std::size_t k = 0; k < n; k++) {
    std::size_t p = k;
    for (std::size_t i = k + 1; i < n; i++)
      if (std::fabs(a[i][k]) > std::fabs(a[p][k])) p = i;
    if (a[p][k] == 0) return 0;
    if (p != k) { std::swap(a[p], a[k]); det = -det; }
    det *= a[k][k];
    for (std::size_t i = k + 1; i < n; i++) {
      const double f = a[i][k] / a[k][k];
      for (std::size_t j = k; j < n; j++) a[i][j] -= f * a[k][j];
    }
  }
  return det;
}

template <std::size_t n>
inline Tensor2<n> inverse(Tensor2<n> a) {
  Tensor2<n> b{};
  for (std::size_t i = 0; i < n; i++) b[i][i] = 1;
  for (std::size_t k = 0; k < n; k++) {
    std::size_t p = k;
    for (std::size_t i = k + 1; i < n; i++)
      if (std::fabs(a[i][k]) > std::fabs(a[p][k])) p = i;
    if (a[p][k] == 0) throw std::domain_error("physis::geometry::inverse: the metric is degenerate here");
    std::swap(a[p], a[k]);
    std::swap(b[p], b[k]);
    const double d = a[k][k];
    for (std::size_t j = 0; j < n; j++) { a[k][j] /= d; b[k][j] /= d; }
    for (std::size_t i = 0; i < n; i++) {
      if (i == k || a[i][k] == 0) continue;
      const double f = a[i][k];
      for (std::size_t j = 0; j < n; j++) { a[i][j] -= f * a[k][j]; b[i][j] -= f * b[k][j]; }
    }
  }
  return b;
}

// g(u, v) = g_ab u^a v^b.
template <std::size_t n>
inline double contract(const Tensor2<n>& g, const Point<n>& u, const Point<n>& v) {
  double s = 0;
  for (std::size_t a = 0; a < n; a++)
    for (std::size_t b = 0; b < n; b++) s += g[a][b] * u[a] * v[b];
  return s;
}

// The number of negative directions of a symmetric form, by Sylvester's law of inertia: the
// signs of the pivots of symmetric elimination (with diagonal pivoting, and a symmetric
// pair rotation when every remaining diagonal entry vanishes).
template <std::size_t n>
inline std::size_t negative_directions(Tensor2<n> a) {
  std::size_t negative = 0;
  for (std::size_t k = 0; k < n; k++) {
    std::size_t p = k;
    for (std::size_t i = k + 1; i < n; i++)
      if (std::fabs(a[i][i]) > std::fabs(a[p][p])) p = i;
    if (std::fabs(a[p][p]) < 1e-300) {
      // Every remaining diagonal entry is zero: add row and column j to k, for a j with a_kj != 0.
      std::size_t j = n;
      for (std::size_t i = k + 1; i < n && j == n; i++)
        if (a[k][i] != 0) j = i;
      if (j == n) continue;  // a zero row: a degenerate direction, counted neither way
      for (std::size_t i = 0; i < n; i++) a[k][i] += a[j][i];
      for (std::size_t i = 0; i < n; i++) a[i][k] += a[i][j];
      p = k;
    }
    if (p != k) {
      std::swap(a[p], a[k]);
      for (auto& row : a) std::swap(row[p], row[k]);
    }
    if (a[k][k] < 0) negative++;
    for (std::size_t i = k + 1; i < n; i++) {
      const double f = a[i][k] / a[k][k];
      for (std::size_t j = k; j < n; j++) a[i][j] -= f * a[k][j];
    }
    for (std::size_t i = k + 1; i < n; i++) a[k][i] = a[i][k] = 0;
  }
  return negative;
}

// ---- standard Riemannian metrics ----------------------------------------------------------------------------------------

// Euclidean space in Cartesian coordinates.
template <std::size_t N>
struct Euclidean {
  static constexpr std::size_t n = N;
  template <class T>
  Matrix<T, n> operator()(const std::array<T, n>&) const {
    Matrix<T, n> g{};
    for (std::size_t i = 0; i < n; i++) g[i][i] = T(1.0);
    return g;
  }
};

// The Euclidean plane in polar coordinates (r, φ): ds² = dr² + r² dφ².
struct PolarPlane {
  static constexpr std::size_t n = 2;
  template <class T>
  Matrix<T, n> operator()(const std::array<T, n>& x) const {
    Matrix<T, n> g{};
    g[0][0] = T(1.0);
    g[1][1] = x[0] * x[0];
    return g;
  }
};

// The round 2-sphere of radius a in colatitude and longitude (θ, φ): ds² = a² (dθ² + sin²θ dφ²).
struct RoundSphere {
  static constexpr std::size_t n = 2;
  double radius = 1;  // m
  template <class T>
  Matrix<T, n> operator()(const std::array<T, n>& x) const {
    using std::sin;
    Matrix<T, n> g{};
    const T s = sin(x[0]);
    g[0][0] = T(radius * radius);
    g[1][1] = radius * radius * s * s;
    return g;
  }
};

}  // namespace physis::geometry
