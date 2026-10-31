// Advection and diffusion of a scalar on a graph, by finite volumes (decision 0016).
//
// The law. A concentration c (an amount per unit volume) on the vertices of a graph, each vertex a
// cell of volume V_i, exchanges through the edges. Across an edge e = (i, j) with conductance g_e
// (m³/s: a diffusive exchange) and advective volume flux u_e (m³/s, positive from i to j), the amount
// carried from i to j per second is, upwinded (LeVeque 2002, Finite Volume Methods for Hyperbolic
// Problems, ch. 4 and 9; Patankar 1980, Numerical Heat Transfer and Fluid Flow, §5.2):
//
//   F_e = g_e (c_i − c_j) + max(u_e, 0) c_i + min(u_e, 0) c_j,
//
// and V_i dc_i/dt = −Σ_{e out of i} F_e + Σ_{e into i} F_e + s_i, s_i the amount a source adds per
// second. The scheme is in flux form, so the content Σ V c is kept to rounding whatever the advection,
// changing by Σ s alone; a uniform concentration is steady exactly when the advection is divergence-
// free at every vertex (Σ_in u = Σ_out u), as a circulation is.
//
// Time. step() is Heun's method (the explicit trapezoidal rule, RK2), strong-stability preserving with
// the explicit Euler step's bound (Shu and Osher 1988, J. Comput. Phys. 77, 439; Gottlieb, Shu and
// Tadmor 2001, SIAM Rev. 43, 89): for h ≤ min_i V_i / (Σ_{e∋i} g_e + Σ_{e leaving i downwind} |u_e|),
// I + hA has no negative entry, so Euler, and therefore Heun, keeps concentrations non-negative and,
// for divergence-free advection, within their initial range (a discrete maximum principle).
// stable_step() returns that bound; a larger step is not refused, and may oscillate.
//
// The steady state with given sources (Σ s = 0) and content is the solution of the linear system
// A c = −V⁻¹ s with one equation replaced by Σ V c = content, solved densely by calx (small graphs:
// O(n³)); it is unique when the graph of exchanges is connected.
//
// Accuracy: first order in space (upwind), second in time. Checks: tests/check_transport.cpp.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include <calx/linalg/dense.hpp>

namespace physis::continuum {

struct GraphTransport {
  std::vector<double> volume;                     // V_i, m³ (or any measure of a vertex)
  std::vector<std::array<std::size_t, 2>> edges;  // (i, j)
  std::vector<double> conductance, advection;     // per edge: g_e ≥ 0 and u_e (m³/s); empty means zero

  // dc/dt for concentrations c and sources s (amount per second at each vertex; empty: none).
  std::vector<double> rate(std::span<const double> c, std::span<const double> source = {}) const {
    validate(c, source, "physis::continuum::GraphTransport::rate", true);
    const std::vector<double> F = flux(c);
    std::vector<double> r(volume.size(), 0.0);
    for (std::size_t e = 0; e < edges.size(); e++) {
      r[edges[e][0]] -= F[e];
      r[edges[e][1]] += F[e];
    }
    for (std::size_t i = 0; i < r.size(); i++) r[i] = (r[i] + (source.empty() ? 0.0 : source[i])) / volume[i];
    return r;
  }

  // F_e, the amount carried along each edge from i to j per second.
  std::vector<double> flux(std::span<const double> c) const {
    validate(c, {}, "physis::continuum::GraphTransport::flux", true);
    std::vector<double> F(edges.size());
    for (std::size_t e = 0; e < edges.size(); e++) {
      const double ci = c[edges[e][0]], cj = c[edges[e][1]];
      const double g = conductance.empty() ? 0.0 : conductance[e], u = advection.empty() ? 0.0 : advection[e];
      F[e] = g * (ci - cj) + std::max(u, 0.0) * ci + std::min(u, 0.0) * cj;
    }
    return F;
  }

  // The explicit Euler bound min_i V_i / (outflow rate of i), s (infinite when nothing flows).
  double stable_step() const {
    validate({}, {}, "physis::continuum::GraphTransport::stable_step");
    std::vector<double> out(volume.size(), 0.0);
    for (std::size_t e = 0; e < edges.size(); e++) {
      const double g = conductance.empty() ? 0.0 : conductance[e], u = advection.empty() ? 0.0 : advection[e];
      out[edges[e][0]] += g + std::max(u, 0.0);
      out[edges[e][1]] += g + std::max(-u, 0.0);
    }
    double h = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < volume.size(); i++)
      if (out[i] > 0) h = std::min(h, volume[i] / out[i]);
    return h;
  }

  // One step of Heun's method (RK2), the sources held over it.
  std::vector<double> step(std::span<const double> c, double h, std::span<const double> source = {}) const {
    if (!(h > 0 && std::isfinite(h))) throw std::invalid_argument("physis::continuum::GraphTransport::step: a step h > 0");
    const std::vector<double> k1 = rate(c, source);
    std::vector<double> y(c.size());
    for (std::size_t i = 0; i < c.size(); i++) y[i] = c[i] + h * k1[i];
    const std::vector<double> k2 = rate(y, source);
    for (std::size_t i = 0; i < c.size(); i++) y[i] = c[i] + 0.5 * h * (k1[i] + k2[i]);
    return y;
  }

  // Σ V c.
  double content(std::span<const double> c) const {
    if (c.size() != volume.size()) throw std::invalid_argument("physis::continuum::GraphTransport::content: one concentration for each vertex");
    double s = 0;
    for (std::size_t i = 0; i < c.size(); i++) s += volume[i] * c[i];
    return s;
  }

  // A, with dc/dt = A c + V⁻¹ s (dense; small graphs).
  calx::linalg::Matrix matrix() const {
    validate({}, {}, "physis::continuum::GraphTransport::matrix");
    const std::size_t n = volume.size();
    calx::linalg::Matrix A(n);
    for (std::size_t e = 0; e < edges.size(); e++) {
      const std::size_t i = edges[e][0], j = edges[e][1];
      const double g = conductance.empty() ? 0.0 : conductance[e], u = advection.empty() ? 0.0 : advection[e];
      const double a = g + std::max(u, 0.0), b = g - std::min(u, 0.0);  // F_e = a c_i − b c_j
      A(i, i) -= a, A(i, j) += b;
      A(j, i) += a, A(j, j) -= b;
    }
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++) A(i, j) /= volume[i];
    return A;
  }

  // The steady state with sources s (Σ s = 0) and the given content Σ V c.
  std::vector<double> steady(std::span<const double> source, double content) const {
    validate({}, source, "physis::continuum::GraphTransport::steady");
    const std::size_t n = volume.size();
    if (n == 0) return {};
    double total = 0, scale = 0;
    for (std::size_t i = 0; i < source.size(); i++) total += source[i], scale += std::fabs(source[i]);
    if (std::fabs(total) > 1e-12 * std::max(scale, 1e-300)) throw std::invalid_argument("physis::continuum::GraphTransport::steady: the sources sum to zero (else no state is steady)");
    calx::linalg::Matrix M = matrix(), B(n);
    for (std::size_t i = 0; i < n; i++) B(i, 0) = source.empty() ? 0.0 : -source[i] / volume[i];
    for (std::size_t j = 0; j < n; j++) M(n - 1, j) = volume[j];
    B(n - 1, 0) = content;
    calx::linalg::Matrix X;
    try {
      X = calx::linalg::solve(M, B);
    } catch (const std::runtime_error&) {
      throw std::invalid_argument("physis::continuum::GraphTransport::steady: no unique steady state (the graph of exchanges is not connected)");
    }
    std::vector<double> c(n);
    for (std::size_t i = 0; i < n; i++) c[i] = X(i, 0);
    return c;
  }

  void validate(std::span<const double> c, std::span<const double> source, const char* where, bool concentrations = false) const {
    const std::size_t n = volume.size();
    for (double v : volume)
      if (!(v > 0 && std::isfinite(v))) throw std::invalid_argument(std::string(where) + ": volumes > 0");
    if (!conductance.empty() && conductance.size() != edges.size()) throw std::invalid_argument(std::string(where) + ": one conductance for each edge");
    if (!advection.empty() && advection.size() != edges.size()) throw std::invalid_argument(std::string(where) + ": one advective flux for each edge");
    for (std::size_t e = 0; e < edges.size(); e++) {
      if (edges[e][0] >= n || edges[e][1] >= n || edges[e][0] == edges[e][1]) throw std::invalid_argument(std::string(where) + ": an edge between two distinct vertices");
      if (!conductance.empty() && !(conductance[e] >= 0)) throw std::invalid_argument(std::string(where) + ": conductances ≥ 0");
    }
    if ((concentrations || !c.empty()) && c.size() != n) throw std::invalid_argument(std::string(where) + ": one concentration for each vertex");
    if (!source.empty() && source.size() != n) throw std::invalid_argument(std::string(where) + ": one source for each vertex");
  }
};

}  // namespace physis::continuum
