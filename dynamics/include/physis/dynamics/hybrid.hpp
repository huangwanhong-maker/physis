// Hybrid dynamical systems: a flow between events, a guard that triggers an event, a jump
// that the event applies, and a dwell after it.
//
// A hybrid system on R^n (Goebel, Sanfelice and Teel 2012, Hybrid Dynamical Systems,
// Princeton, ch. 1-2) here has
//
// - a flow x' = f(x; p, u) driven by one input u held constant over each step: *linear*,
//   x' = A x + B u + c, integrated exactly (linear.hpp), or *nonlinear*, integrated by a
//   one-step method (calx/integrators), or by the exponential method of Rush and Larsen
//   when the system declares gates;
// - a guard G = {x : guard(p, x)}, tested at the grid points t_k = k h: the event is placed
//   at the first grid point at which the state lies in G (grid-constrained events; the
//   error in its time is below h);
// - a jump x <- R(x) applied at the event;
// - a dwell of duration dwell(p) after the event, during which the components marked
//   `held` are held still and the guard is not tested (the refractory period of a neuron).
//
// A system is a type with static members, checked by the concepts below; its parameters
// are a struct, so that one type serves every parameter set and the compiler generates the
// arithmetic of each system. Integration of one trajectory is `simulate`; a library that
// integrates millions of copies at once (neuroscience's spiking networks) uses the same
// pieces, step and propagator, in its own loops.

#pragma once

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <vector>

#include <calx/integrators/runge_kutta.hpp>
#include <calx/linalg/dense.hpp>

#include "linear.hpp"

namespace physis::dynamics {

template <std::size_t N>
using State = std::array<double, N>;

enum class Method { exact, euler, heun, rk4, dopri5, rush_larsen };

template <class S>
concept System = requires {
  { S::n } -> std::convertible_to<std::size_t>;
  typename S::Params;
};

template <class S>
concept LinearFlow = System<S> && requires(const typename S::Params& p, Matrix& A, State<S::n>& B, State<S::n>& c) {
  { S::system(p, A, B, c) };
};

template <class S>
concept NonlinearFlow = System<S> && requires(const typename S::Params& p, const State<S::n>& x, double u, State<S::n>& dx) {
  { S::derivative(p, x, u, dx) };
};

// Gates for Rush and Larsen's method: for each component marked in S::gate, its steady state
// and time constant at the present state; the other components flow by derivative().
template <class S>
concept Gated = NonlinearFlow<S> && requires(const typename S::Params& p, const State<S::n>& x, State<S::n>& inf, State<S::n>& tau) {
  { S::gates(p, x, inf, tau) };
  { S::gate } -> std::convertible_to<std::array<bool, S::n>>;
};

template <class S>
concept Hybrid = (LinearFlow<S> || NonlinearFlow<S>) && requires(const typename S::Params& p, State<S::n>& x, const State<S::n>& cx) {
  { S::guard(p, cx) } -> std::convertible_to<bool>;
  { S::jump(p, x) };
  { S::dwell(p) } -> std::convertible_to<double>;
  { S::held } -> std::convertible_to<std::array<bool, S::n>>;
};

// The exact propagator of a linear system over h, with `held` components still. A system
// may declare the typical sizes of its components (S::scale) and of its input
// (S::input_scale), which keep the exponential accurate in SI units (linear.hpp).
template <LinearFlow S>
Propagator<S::n> linear_propagator(const typename S::Params& p, double h, const std::array<bool, S::n>& held = {}) {
  Matrix A(S::n);
  State<S::n> B{}, c{};
  S::system(p, A, B, c);
  std::array<double, S::n> xsize{};
  std::array<double, 1> usize{};
  if constexpr (requires { S::scale; }) for (std::size_t i = 0; i < S::n; i++) xsize[i] = S::scale[i];
  if constexpr (requires { S::input_scale; }) usize[0] = S::input_scale;
  return propagator<S::n, 1>(A, B, c, h, held, xsize, usize);
}

// The settings of a one-step method.
struct Stepping {
  Method method = Method::rk4;
  int substeps = 1;                              // fixed-step methods take h / substeps
  calx::integrators::Tolerance tolerance{1e-6, 1e-6};  // dopri5, relative to each component's scale
};

// One step of h of a nonlinear system, the `held` components still. `sub` carries
// Dormand-Prince's substep from step to step; `scale` is each component's typical size.
template <NonlinearFlow S>
void step(const typename S::Params& p, State<S::n>& x, double u, double h, const Stepping& how, const std::array<bool, S::n>& held,
          double& sub, const State<S::n>& scale) {
  auto f = [&](const State<S::n>& y, State<S::n>& dy) {
    S::derivative(p, y, u, dy);
    for (std::size_t i = 0; i < S::n; i++)
      if (held[i]) dy[i] = 0;
  };
  const int k = how.substeps > 0 ? how.substeps : 1;
  const double hs = h / k;
  switch (how.method) {
    case Method::euler: for (int i = 0; i < k; i++) calx::integrators::euler<S::n>(f, x, hs); break;
    case Method::heun: for (int i = 0; i < k; i++) calx::integrators::heun<S::n>(f, x, hs); break;
    case Method::rk4: for (int i = 0; i < k; i++) calx::integrators::rk4<S::n>(f, x, hs); break;
    case Method::dopri5: calx::integrators::dopri5<S::n>(f, x, h, scale, how.tolerance, sub); break;
    case Method::rush_larsen:
      if constexpr (Gated<S>) {
        for (int i = 0; i < k; i++) {
          // Gates move exactly towards their steady states with the coefficients frozen at
          // the start of the substep; the rest moves by forward Euler (Rush and Larsen 1978,
          // IEEE Trans. Biomed. Eng. 25, 389-392). Order 1.
          State<S::n> inf{}, tau{}, dx{};
          S::gates(p, x, inf, tau);
          f(x, dx);
          for (std::size_t c = 0; c < S::n; c++) {
            if (held[c]) continue;
            if (S::gate[c]) x[c] = inf[c] + (x[c] - inf[c]) * std::exp(-hs / tau[c]);
            else x[c] += hs * dx[c];
          }
        }
        break;
      } else {
        throw std::invalid_argument("physis::dynamics::step: Rush and Larsen's method needs gates");
      }
    case Method::exact: throw std::invalid_argument("physis::dynamics::step: a nonlinear flow has no exact propagator here");
  }
}

// One trajectory on the grid t_k = k h.
template <std::size_t N>
struct Trajectory {
  double h = 0;
  std::vector<double> t;              // sampled times
  std::vector<State<N>> x;            // sampled states
  std::vector<double> events;         // the times of the events
};

// Simulates from x0 over [0, T] with input u(t) held over each step (evaluated at its
// start), recording every `every` steps. Linear systems are integrated exactly whatever
// `how` says; nonlinear ones by `how`.
template <Hybrid S>
Trajectory<S::n> simulate(const typename S::Params& p, State<S::n> x, const std::function<double(double)>& u, double T, double h,
                          const Stepping& how = {}, int every = 1, const State<S::n>& scale = [] { State<S::n> s; s.fill(1); return s; }()) {
  Trajectory<S::n> out;
  out.h = h;
  const auto steps = static_cast<long long>(std::llround(T / h));
  const auto dwell_steps = static_cast<long long>(std::llround(S::dwell(p) / h));
  long long dwelling = 0;
  double sub = 0;
  const std::array<bool, S::n> none{};
  std::array<Propagator<S::n>, 2> P;
  if constexpr (LinearFlow<S>) P = {linear_propagator<S>(p, h), linear_propagator<S>(p, h, S::held)};
  for (long long k = 0; k <= steps; k++) {
    if (every > 0 && k % every == 0) {
      out.t.push_back(k * h);
      out.x.push_back(x);
    }
    if (k == steps) break;
    const double uk = u(k * h);
    const bool hold = dwelling > 0;
    if constexpr (LinearFlow<S>) P[hold ? 1 : 0].apply(x, uk);
    else step<S>(p, x, uk, h, how, hold ? S::held : none, sub, scale);
    if (hold) dwelling--;
    else if (S::guard(p, x)) {
      out.events.push_back((k + 1) * h);
      S::jump(p, x);
      dwelling = dwell_steps;
    }
  }
  return out;
}

}  // namespace physis::dynamics
