// The Euler–Lagrange equations of any Lagrangian, solved for the acceleration (decision 0017;
// docs/theory/discrete-mechanics.md §8).
//
// For L(q, v) and a generalized force Q, the equations d/dt ∂L/∂v − ∂L/∂q = Q read
//
//   (∂²L/∂v∂v) a = Q + ∂L/∂q − (∂²L/∂v∂q) v,
//
// since d/dt ∂L/∂v = (∂²L/∂v∂v) a + (∂²L/∂v∂q) v for an autonomous L. The Lagrangian is written once over
// its number type; on calx's jets of order 2 in the 2N arguments (q, v) one evaluation gives ∂L/∂q and the
// blocks of the Hessian, exact to rounding, and calx's small dense solve the acceleration. So any Lagrangian
// can be moved by the explicit methods too: a tree with the kinetic energy of the fluid it moves (fluid.hpp),
// whose mass matrix Newton–Euler does not know.
//
// For a tree alone it is the tree's forward dynamics (the checks hold them equal to rounding).

#pragma once

#include <array>
#include <cstddef>

#include <calx/differentiation/jet.hpp>
#include <calx/integrators/runge_kutta.hpp>
#include <calx/linalg/small.hpp>

namespace physis::mechanics {

template <std::size_t N, class Lagrangian>
std::array<double, N> lagrangian_acceleration(const Lagrangian& L, double t, const std::array<double, N>& q, const std::array<double, N>& v,
                                              const std::array<double, N>& Q) {
  std::array<double, 2 * N> z;
  for (std::size_t k = 0; k < N; k++) z[k] = q[k], z[N + k] = v[k];
  const auto x = calx::differentiation::variables<2 * N, 2>(z);
  using J = calx::differentiation::Taylor<2 * N, 2>;
  std::array<J, N> qj, vj;
  for (std::size_t k = 0; k < N; k++) qj[k] = x[k], vj[k] = x[N + k];
  const J l = L(t, qj, vj);
  std::array<std::array<double, N>, N> M;
  std::array<double, N> rhs;
  for (std::size_t i = 0; i < N; i++) {
    rhs[i] = Q[i] + l.d[i];
    for (std::size_t j = 0; j < N; j++) {
      M[i][j] = l.h[N + i][N + j];
      rhs[i] -= l.h[N + i][j] * v[j];
    }
  }
  return calx::linalg::solve<N>(M, rhs);
}

// Semi-implicit Euler and the classical Runge–Kutta method on a Lagrangian's equations, Q(t, q, v) the
// generalized force of everything not in L, held at the step's time in Runge–Kutta's stages' time.
template <std::size_t N, class Lagrangian, class Force>
void semi_implicit_euler_lagrangian(const Lagrangian& L, const Force& Q, double& t, std::array<double, N>& q, std::array<double, N>& v, double h) {
  const auto a = lagrangian_acceleration<N>(L, t, q, v, Q(t, q, v));
  for (std::size_t k = 0; k < N; k++) v[k] += h * a[k];
  for (std::size_t k = 0; k < N; k++) q[k] += h * v[k];
  t += h;
}

template <std::size_t N, class Lagrangian, class Force>
void runge_kutta4_lagrangian(const Lagrangian& L, const Force& Q, double& t, std::array<double, N>& q, std::array<double, N>& v, double h) {
  calx::integrators::Vec<2 * N> x;
  for (std::size_t k = 0; k < N; k++) x[k] = q[k], x[N + k] = v[k];
  const double t0 = t;
  calx::integrators::rk4<2 * N>(
      [&](const calx::integrators::Vec<2 * N>& y, calx::integrators::Vec<2 * N>& dy) {
        std::array<double, N> qq, vv;
        for (std::size_t k = 0; k < N; k++) qq[k] = y[k], vv[k] = y[N + k];
        const auto a = lagrangian_acceleration<N>(L, t0, qq, vv, Q(t0, qq, vv));
        for (std::size_t k = 0; k < N; k++) dy[k] = vv[k], dy[N + k] = a[k];
      },
      x, h);
  for (std::size_t k = 0; k < N; k++) q[k] = x[k], v[k] = x[N + k];
  t += h;
}

}  // namespace physis::mechanics
