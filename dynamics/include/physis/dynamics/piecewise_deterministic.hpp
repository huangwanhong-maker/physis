// Piecewise-deterministic Markov processes: a flow between jumps, jumps at a rate along it, and forced jumps
// when the flow reaches the boundary of its domain (decision 0018).
//
// The process (Davis 1984, J. R. Statist. Soc. B 46, 353–388; Davis 1993, Markov Models and Optimization,
// ch. 2). On an open domain E = {x : g(x) > 0} of ℝⁿ: between jumps x follows x′ = F(x); jumps come at the
// rate λ(x) along the flow and redraw x from a kernel R(x, ·); when the flow reaches the boundary ∂E
// (g = 0) the state is redrawn at once from a boundary kernel R_∂(x, ·). The extended generator is
//
//   (A f)(x) = ∇f(x)·F(x) + λ(x) ∫ [f(y) − f(x)] R(x, dy)    for x ∈ E,
//
// on functions f that meet the boundary condition f(x) = ∫ f(y) R_∂(x, dy) at x ∈ ∂E (Davis 1993, Thm 26.14):
// the forced jump is a condition on the domain of A, not a term of A. Dynkin's formula d/dt E f(x_t) =
// E (A f)(x_t) holds for f in the domain; for an f that breaks the boundary condition the flux through ∂E
// adds a term, which is why the condition cannot be dropped (tests/check_piecewise_deterministic.cpp).
//
// The construction (Davis 1984, §2). From x, draw ε ~ Exp(1); follow the flow while integrating
// Λ(t) = ∫₀ᵗ λ(x_s) ds; jump by R when Λ reaches ε, by R_∂ when g first reaches 0, whichever comes first.
// Here the flow and Λ are integrated together by classical Runge–Kutta steps of a fixed h; within a step the
// time at which Λ = ε, or g = 0, is located by bisection on the step length, re-integrating from the step's
// start, to a relative 10⁻¹⁰ of h. The error is that of the integrator, O(h⁴) per unit time for smooth F and
// λ; the times of jumps are exact up to it.
//
// The state is a std::vector<double>; F, λ, g and the kernels are given as functions. Checks:
// tests/check_piecewise_deterministic.cpp.

#pragma once

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <calx/random/stream.hpp>

namespace physis::dynamics {

struct PiecewiseDeterministic {
  using Vec = std::vector<double>;
  std::function<void(const Vec& x, Vec& dx)> flow;                    // F(x)
  std::function<double(const Vec& x)> rate;                           // λ(x) ≥ 0; empty for none
  std::function<Vec(const Vec& x, calx::random::Stream& rng)> jump;    // a draw from R(x, ·)
  std::function<double(const Vec& x)> inside;                         // g(x) > 0 in the domain; empty for ℝⁿ
  std::function<Vec(const Vec& x, calx::random::Stream& rng)> reset;   // a draw from R_∂(x, ·)
};

struct PdmpEvent {
  double time = 0;
  bool forced = false;  // at the boundary, not at the rate
  std::vector<double> before, after;
};

struct PdmpPath {
  std::vector<PdmpEvent> events;
  std::vector<double> times;              // sample times (every `every` steps, and at events)
  std::vector<std::vector<double>> states;
  double end = 0;
  bool stopped = false;  // too many events (a possible accumulation)
};

namespace detail {

// One RK4 step of the flow augmented with Λ: y = (x, Λ).
inline void pdmp_step(const PiecewiseDeterministic& P, std::vector<double>& x, double& L, double h) {
  const std::size_t n = x.size();
  std::vector<double> k1(n), k2(n), k3(n), k4(n), y(n);
  const auto lam = [&](const std::vector<double>& z) { return P.rate ? P.rate(z) : 0.0; };
  P.flow(x, k1);
  const double l1 = lam(x);
  for (std::size_t i = 0; i < n; i++) y[i] = x[i] + 0.5 * h * k1[i];
  P.flow(y, k2);
  const double l2 = lam(y);
  for (std::size_t i = 0; i < n; i++) y[i] = x[i] + 0.5 * h * k2[i];
  P.flow(y, k3);
  const double l3 = lam(y);
  for (std::size_t i = 0; i < n; i++) y[i] = x[i] + h * k3[i];
  P.flow(y, k4);
  const double l4 = lam(y);
  for (std::size_t i = 0; i < n; i++) x[i] += h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
  L += h / 6 * (l1 + 2 * l2 + 2 * l3 + l4);
}

}  // namespace detail

// A path from x0 over [0, t_end] with steps of h, sampled every `every` steps; at most max_events jumps.
inline PdmpPath simulate(const PiecewiseDeterministic& P, std::vector<double> x0, double t_end, double h, calx::random::Stream& rng,
                         std::size_t every = 1, std::size_t max_events = 1000000) {
  if (!P.flow) throw std::invalid_argument("physis::dynamics::simulate: a flow");
  if (!(h > 0) || !(t_end >= 0)) throw std::invalid_argument("physis::dynamics::simulate: a step h > 0 and an end ≥ 0");
  if (P.rate && !P.jump) throw std::invalid_argument("physis::dynamics::simulate: a rate needs a jump kernel");
  if (P.inside && !P.reset) throw std::invalid_argument("physis::dynamics::simulate: a boundary needs a reset kernel");
  if (P.inside && !(P.inside(x0) > 0)) throw std::invalid_argument("physis::dynamics::simulate: the start lies outside the domain");
  PdmpPath out;
  std::vector<double> x = std::move(x0);
  double t = 0, L = 0, eps = rng.exponential();
  std::size_t steps = 0;
  out.times.push_back(0), out.states.push_back(x);
  const auto outside = [&](const std::vector<double>& z) { return P.inside && !(P.inside(z) > 0); };
  while (t < t_end) {
    const double dt = std::min(h, t_end - t);
    std::vector<double> y = x;
    double M = L;
    detail::pdmp_step(P, y, M, dt);
    const bool fires = M >= eps, leaves = outside(y);
    if (!fires && !leaves) {
      x = std::move(y), L = M, t += dt;
      if (++steps % every == 0 || t >= t_end) out.times.push_back(t), out.states.push_back(x);
      continue;
    }
    // Locate the first of the two within the step: the largest s whose state is inside with Λ < ε.
    double lo = 0, hi = dt;
    while (hi - lo > 1e-10 * dt) {
      const double mid = 0.5 * (lo + hi);
      std::vector<double> z = x;
      double K = L;
      detail::pdmp_step(P, z, K, mid);
      if (K >= eps || outside(z)) hi = mid;
      else lo = mid;
    }
    std::vector<double> z = x;
    double K = L;
    detail::pdmp_step(P, z, K, hi);
    // At the boundary the jump is forced (the rate cannot act outside the domain); otherwise Λ reached ε.
    PdmpEvent e;
    e.time = t + hi;
    e.forced = outside(z);
    e.before = z;
    e.after = e.forced ? P.reset(z, rng) : P.jump(z, rng);
    if (P.inside && !(P.inside(e.after) > 0)) throw std::runtime_error("physis::dynamics::simulate: a jump lands outside the domain");
    x = e.after;
    t = e.time;
    L = 0, eps = rng.exponential();
    out.events.push_back(std::move(e));
    out.times.push_back(t), out.states.push_back(x);
    if (out.events.size() >= max_events) { out.stopped = true; break; }
  }
  out.end = t;
  return out;
}

}  // namespace physis::dynamics
