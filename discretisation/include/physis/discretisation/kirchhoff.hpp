// Networks of conductances on a graph: Kirchhoff's laws with given potentials and injections
// (decision 0022).
//
// A network is a graph whose edges e = (a, b) carry conductances g_e > 0: an electrical circuit of
// resistors, a network of tubes carrying a slow flow (Hagen–Poiseuille conductances), a porous medium.
// The flow along e is F_e = g_e (φ_a − φ_b) (Ohm's law), and at every vertex whose potential is not fixed
// the flows balance the injection s_v (Kirchhoff's current law; Kirchhoff 1847):
//
//   Σ_{e ∋ v} ± F_e = s_v,  that is  (L φ)_v = s_v,  L = Bᵀ G B the weighted graph Laplacian.
//
// With some potentials fixed (the terminals), the others solve the reduced system L_ff φ_f = s_f − L_fc φ_c,
// which is positive definite exactly when every component of the free vertices touches a terminal; a
// network without one is refused with that reason. The solution minimises the dissipation Σ g_e (Δφ_e)²
// among potentials with the terminals' values (Thomson's principle), and the effective conductance between
// two terminals held one unit apart is the total flow between them.
//
// This is the graph form of the flows with conductances on a 2-complex in hodge.hpp, kept apart because a
// graph has no faces; a network whose conductances change by a local rule (Physarum, coupled learning) is
// solved by it at every step. Small networks are solved densely, larger ones by calx's conjugate gradients
// to a relative residual of 1e−13. Checks: tests/check_kirchhoff.cpp.
#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include <calx/linalg/conjugate_gradient.hpp>
#include <calx/linalg/dense.hpp>
#include <calx/linalg/sparse.hpp>

namespace physis::discretisation {

struct ConductanceEdge {
  std::size_t a = 0, b = 0;
  double conductance = 1;
};

struct ConductanceNetwork {
  std::size_t vertices = 0;
  std::vector<ConductanceEdge> edges;
  std::size_t add(std::size_t a, std::size_t b, double g) {
    edges.push_back({a, b, g});
    return edges.size() - 1;
  }
};

struct KirchhoffSolution {
  std::vector<double> potential;  // at every vertex
  std::vector<double> flow;       // along every edge, positive from a to b
  double dissipation = 0;         // Σ g (Δφ)²
};

// fixed: (vertex, potential) pairs; injection: s_v at each vertex (ignored at fixed ones), or empty for none.
inline KirchhoffSolution solve_kirchhoff(const ConductanceNetwork& net, const std::vector<std::pair<std::size_t, double>>& fixed,
                                         const std::vector<double>& injection = {}) {
  const std::size_t n = net.vertices;
  if (!injection.empty() && injection.size() != n) throw std::invalid_argument("physis::discretisation::solve_kirchhoff: one injection per vertex");
  std::vector<char> is_fixed(n, 0);
  KirchhoffSolution out;
  out.potential.assign(n, 0.0);
  for (const auto& [v, phi] : fixed) {
    if (v >= n) throw std::out_of_range("physis::discretisation::solve_kirchhoff: a fixed vertex out of range");
    is_fixed[v] = 1, out.potential[v] = phi;
  }
  for (const auto& e : net.edges)
    if (e.a >= n || e.b >= n || !(e.conductance > 0)) throw std::invalid_argument("physis::discretisation::solve_kirchhoff: edges between vertices, with positive conductances");
  // Number the free vertices.
  std::vector<std::size_t> index(n, n);
  std::size_t m = 0;
  for (std::size_t v = 0; v < n; v++)
    if (!is_fixed[v]) index[v] = m++;
  // Each free component must touch a terminal: a search from the terminals.
  {
    std::vector<std::vector<std::size_t>> adj(n);
    for (const auto& e : net.edges) adj[e.a].push_back(e.b), adj[e.b].push_back(e.a);
    std::vector<char> seen(is_fixed);
    std::vector<std::size_t> stack;
    for (std::size_t v = 0; v < n; v++)
      if (is_fixed[v]) stack.push_back(v);
    while (!stack.empty()) {
      const std::size_t v = stack.back();
      stack.pop_back();
      for (std::size_t w : adj[v])
        if (!seen[w]) seen[w] = 1, stack.push_back(w);
    }
    for (std::size_t v = 0; v < n; v++)
      if (!seen[v]) throw std::invalid_argument("physis::discretisation::solve_kirchhoff: a part of the network touches no fixed potential, so its potentials are undetermined");
  }
  std::vector<double> rhs(m, 0.0);
  for (std::size_t v = 0; v < n; v++)
    if (!is_fixed[v] && !injection.empty()) rhs[index[v]] = injection[v];
  if (m > 0) {
    std::vector<std::uint32_t> rows, cols;
    std::vector<double> values;
    std::vector<double> diag(m, 0.0);
    for (const auto& e : net.edges) {
      const bool fa = is_fixed[e.a], fb = is_fixed[e.b];
      if (!fa) diag[index[e.a]] += e.conductance;
      if (!fb) diag[index[e.b]] += e.conductance;
      if (!fa && !fb) {
        rows.push_back(static_cast<std::uint32_t>(index[e.a])), cols.push_back(static_cast<std::uint32_t>(index[e.b])), values.push_back(-e.conductance);
        rows.push_back(static_cast<std::uint32_t>(index[e.b])), cols.push_back(static_cast<std::uint32_t>(index[e.a])), values.push_back(-e.conductance);
      } else if (!fa) {
        rhs[index[e.a]] += e.conductance * out.potential[e.b];
      } else if (!fb) {
        rhs[index[e.b]] += e.conductance * out.potential[e.a];
      }
    }
    std::vector<double> x;
    if (m <= 256) {
      calx::linalg::Matrix A(m);
      for (std::size_t i = 0; i < m; i++) A(i, i) = diag[i];
      for (std::size_t k = 0; k < rows.size(); k++) A(rows[k], cols[k]) += values[k];
      calx::linalg::Matrix b(m);
      // solve(A, B) takes a square right-hand side; one column used.
      for (std::size_t i = 0; i < m; i++) b(i, 0) = rhs[i];
      const auto X = calx::linalg::solve(A, b);
      x.resize(m);
      for (std::size_t i = 0; i < m; i++) x[i] = X(i, 0);
    } else {
      for (std::size_t i = 0; i < m; i++) rows.push_back(static_cast<std::uint32_t>(i)), cols.push_back(static_cast<std::uint32_t>(i)), values.push_back(diag[i]);
      std::vector<std::uint64_t> order;
      const auto A = calx::linalg::compress<std::uint32_t>(m, m, rows, cols, &order);
      std::vector<double> a(values.size());
      for (std::size_t k = 0; k < a.size(); k++) a[k] = values[order[k]];
      const auto r = calx::linalg::conjugate_gradient(A, std::span<const double>(a), std::span<const double>(rhs), 1e-13);
      if (!r.converged) throw std::runtime_error("physis::discretisation::solve_kirchhoff: conjugate gradients did not converge");
      x = r.x;
    }
    for (std::size_t v = 0; v < n; v++)
      if (!is_fixed[v]) out.potential[v] = x[index[v]];
  }
  out.flow.resize(net.edges.size());
  for (std::size_t k = 0; k < net.edges.size(); k++) {
    const auto& e = net.edges[k];
    const double d = out.potential[e.a] - out.potential[e.b];
    out.flow[k] = e.conductance * d;
    out.dissipation += e.conductance * d * d;
  }
  return out;
}

// The net flow leaving each vertex through the edges: Σ outgoing − Σ incoming.
inline std::vector<double> divergence(const ConductanceNetwork& net, const std::vector<double>& flow) {
  std::vector<double> d(net.vertices, 0.0);
  for (std::size_t k = 0; k < net.edges.size(); k++) d[net.edges[k].a] += flow[k], d[net.edges[k].b] -= flow[k];
  return d;
}

// The effective conductance between s and t: the flow from s with φ_s = 1, φ_t = 0.
inline double effective_conductance(const ConductanceNetwork& net, std::size_t s, std::size_t t) {
  const auto sol = solve_kirchhoff(net, {{s, 1.0}, {t, 0.0}});
  return divergence(net, sol.flow)[s];
}

}  // namespace physis::discretisation
