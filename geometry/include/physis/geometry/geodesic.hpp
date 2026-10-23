// Geodesics of a metric field: the curves that parallel-transport their own tangent.
//
// In a chart a geodesic x^a(λ) with affine parameter λ obeys (MTW 1973, eq. 13.37; Carroll
// 2004, eq. 3.44)
//   d²x^a/dλ² + Γ^a_bc (dx^b/dλ)(dx^c/dλ) = 0,
// which, with u^a = dx^a/dλ, is the first-order system on the phase space (x, u) of
// dimension 2n
//   dx^a/dλ = u^a,   du^a/dλ = -Γ^a_bc u^b u^c,
// a dynamical system in the sense of physis `dynamics`, integrated here by calx's one-step
// methods: the classical Runge-Kutta method with a fixed step, or Dormand and Prince's
// embedded pair with error control.
//
// Along an exact geodesic g(u, u) is constant (-1 for a timelike geodesic parametrised by
// proper time, 0 for a null one; MTW eq. 16.4 for c = 1, here with x^0 = ct so that λ is a
// length), and for every Killing field ξ the quantity g(ξ, u) is conserved. When the metric
// does not depend on a coordinate x^a, ∂_a is a Killing field and the covariant component
// p_a = g_ab u^b is conserved (the energy -p_t and angular momentum p_φ of a static,
// axisymmetric spacetime). Their drift measures the error of the integration.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include <calx/integrators/runge_kutta.hpp>

#include "curvature.hpp"

namespace physis::geometry {

// A point of phase space: (x^0 .. x^{n-1}, u^0 .. u^{n-1}).
template <std::size_t n>
using Phase = std::array<double, 2 * n>;

template <std::size_t n>
inline Point<n> position(const Phase<n>& y) {
  Point<n> x;
  for (std::size_t a = 0; a < n; a++) x[a] = y[a];
  return x;
}
template <std::size_t n>
inline Point<n> velocity(const Phase<n>& y) {
  Point<n> u;
  for (std::size_t a = 0; a < n; a++) u[a] = y[n + a];
  return u;
}
template <std::size_t n>
inline Phase<n> phase(const Point<n>& x, const Point<n>& u) {
  Phase<n> y;
  for (std::size_t a = 0; a < n; a++) { y[a] = x[a]; y[n + a] = u[a]; }
  return y;
}

// The vector field of the geodesic equation at y.
template <MetricField M>
inline void geodesic_field(const M& metric, const Phase<M::n>& y, Phase<M::n>& dy) {
  constexpr std::size_t n = M::n;
  const auto Gamma = christoffel(metric, position<n>(y));
  for (std::size_t a = 0; a < n; a++) {
    dy[a] = y[n + a];
    double s = 0;
    for (std::size_t b = 0; b < n; b++)
      for (std::size_t c = 0; c < n; c++) s += Gamma[a][b][c] * y[n + b] * y[n + c];
    dy[n + a] = -s;
  }
}

// g(u, u) at the point of y.
template <MetricField M>
inline double norm(const M& metric, const Phase<M::n>& y) {
  constexpr std::size_t n = M::n;
  return contract(metric(position<n>(y)), velocity<n>(y), velocity<n>(y));
}

// The covariant component p_a = g_ab u^b, conserved when the metric does not depend on x^a.
template <MetricField M>
inline double momentum(const M& metric, const Phase<M::n>& y, std::size_t a) {
  constexpr std::size_t n = M::n;
  const auto g = metric(position<n>(y));
  double s = 0;
  for (std::size_t b = 0; b < n; b++) s += g[a][b] * y[n + b];
  return s;
}

// Completes u by its component u^0 so that g(u, u) = kappa (-1: timelike, 0: null), taking
// the root with u^0 > 0, which is future-directed in a chart whose x^0 increases to the future.
template <MetricField M>
inline Point<M::n> complete_velocity(const M& metric, const Point<M::n>& x, Point<M::n> u, double kappa) {
  constexpr std::size_t n = M::n;
  const auto g = metric(x);
  // g_00 (u^0)² + 2 b u^0 + c = kappa.
  double b = 0, c = -kappa;
  for (std::size_t i = 1; i < n; i++) {
    b += g[0][i] * u[i];
    for (std::size_t j = 1; j < n; j++) c += g[i][j] * u[i] * u[j];
  }
  const double a = g[0][0];
  if (std::fabs(a) < 1e-300) {
    if (b == 0) throw std::domain_error("physis::geometry::complete_velocity: no solution for u^0");
    u[0] = -c / (2 * b);
    return u;
  }
  const double disc = b * b - a * c;
  if (disc < 0) throw std::domain_error("physis::geometry::complete_velocity: no real u^0 (the spatial part is not allowed for this kind)");
  const double r1 = (-b + std::sqrt(disc)) / a, r2 = (-b - std::sqrt(disc)) / a;
  u[0] = std::max(r1, r2);
  return u;
}

enum class Stepper { rk4, dopri5 };

// One step of affine parameter h from y; with dopri5, tol sets the error control and
// `sub` carries the substep to try next.
template <MetricField M>
inline void geodesic_step(const M& metric, Phase<M::n>& y, double h, Stepper method = Stepper::rk4, calx::integrators::Tolerance tol = {1e-12, 1e-12}) {
  constexpr std::size_t N = 2 * M::n;
  auto f = [&](const Phase<M::n>& z, Phase<M::n>& dz) { geodesic_field(metric, z, dz); };
  if (method == Stepper::rk4) {
    calx::integrators::rk4<N>(f, y, h);
  } else {
    std::array<double, N> scale;
    for (std::size_t i = 0; i < N; i++) scale[i] = std::max(1.0, std::fabs(y[i]));
    double sub = 0;
    calx::integrators::dopri5<N>(f, y, h, scale, tol, sub);
  }
}

// Integrates `steps` steps of h, calling observe(lambda, y) at the start and after each step;
// observe may return false to stop.
template <MetricField M, class Observer>
inline void integrate_geodesic(const M& metric, Phase<M::n>& y, double h, std::size_t steps, Observer&& observe,
                               Stepper method = Stepper::rk4, calx::integrators::Tolerance tol = {1e-12, 1e-12}) {
  double lambda = 0;
  if (!observe(lambda, static_cast<const Phase<M::n>&>(y))) return;
  for (std::size_t k = 0; k < steps; k++) {
    geodesic_step(metric, y, h, method, tol);
    lambda += h;
    if (!observe(lambda, static_cast<const Phase<M::n>&>(y))) return;
  }
}

}  // namespace physis::geometry
