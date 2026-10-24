// Real cochains on a CW complex: the discrete differential forms of discrete exterior calculus.
//
// A k-cochain on a CW complex X (morphe topology) is a real number on each k-cell: a 0-cochain a
// function on the vertices, a 1-cochain a value on each oriented edge (a discrete 1-form), a
// 2-cochain a value on each 2-cell. The coboundary is the transpose of the cellular boundary
// (Hatcher 2002, §3.1; Desbrun, Hirani, Leok and Marsden 2005, §3):
//
//   (d₀ φ)(e) = φ(head e) − φ(tail e),     (d₁ ω)(f) = Σ_{steps of f's word} ±ω(e),
//
// the second the integral of ω round the attaching word of f, so that d₁ d₀ = 0 (a word is a
// closed path) and Stokes's theorem ∫_{∂f} ω = (d₁ ω)(f) holds by definition. The integral of a
// 1-cochain along any edge path is the signed sum of its values.
//
// A 1-cochain is closed when d₁ ω = 0 and exact when ω = d₀ φ. An exact cochain is closed, and its
// integral round every closed path vanishes; conversely a 1-cochain whose integral round every
// closed path vanishes is exact, with φ its integral from a root along a spanning forest, unique up
// to a constant on each component. So the integrals round the fundamental cycles of a spanning
// forest (one for each edge outside it) are complete invariants of ω modulo exact cochains: the
// tree gauge of a lattice gauge theory with group (ℝ, +), and coordinates on H¹(X; ℝ) when ω is
// closed. Multiplicative transports U > 0 enter through ω = log U.

#pragma once

#include <morphe/combinatorics/graph.hpp>
#include <morphe/topology/cw_complex.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace physis::discretisation {

using morphe::topology::CWComplex;
using morphe::topology::Path;
using Cochain = std::vector<double>;

// d₀: functions on the vertices to 1-cochains.
inline Cochain d0(const CWComplex& X, const Cochain& phi) {
  if (phi.size() != X.vertices) throw std::invalid_argument("physis::discretisation::d0: one value per vertex");
  Cochain w(X.edges.size());
  for (std::size_t e = 0; e < X.edges.size(); e++) w[e] = phi[X.edges[e][1]] - phi[X.edges[e][0]];
  return w;
}

// d₁: 1-cochains to 2-cochains, the integral round each 2-cell's word.
inline Cochain d1(const CWComplex& X, const Cochain& omega) {
  if (omega.size() != X.edges.size()) throw std::invalid_argument("physis::discretisation::d1: one value per edge");
  Cochain c(X.faces.size(), 0.0);
  for (std::size_t f = 0; f < X.faces.size(); f++)
    for (const auto& s : X.faces[f].boundary) c[f] += s.sign * omega[s.edge];
  return c;
}

// The transposes: d₀ᵀ takes a value on each edge to its net inflow at each vertex (head minus
// tail), the boundary of a 1-chain; d₁ᵀ takes a value on each 2-cell to its edges, with the signs
// of the words, the boundary of a 2-chain.
inline Cochain d0_transpose(const CWComplex& X, const Cochain& w) {
  if (w.size() != X.edges.size()) throw std::invalid_argument("physis::discretisation::d0_transpose: one value per edge");
  Cochain v(X.vertices, 0.0);
  for (std::size_t e = 0; e < X.edges.size(); e++) {
    v[X.edges[e][1]] += w[e];
    v[X.edges[e][0]] -= w[e];
  }
  return v;
}

inline Cochain d1_transpose(const CWComplex& X, const Cochain& c) {
  if (c.size() != X.faces.size()) throw std::invalid_argument("physis::discretisation::d1_transpose: one value per 2-cell");
  Cochain w(X.edges.size(), 0.0);
  for (std::size_t f = 0; f < X.faces.size(); f++)
    for (const auto& s : X.faces[f].boundary) w[s.edge] += s.sign * c[f];
  return w;
}

// The integral of a 1-cochain along an edge path.
inline double integral(const Cochain& omega, const Path& p) {
  double s = 0;
  for (const auto& step : p) s += step.sign * omega.at(step.edge);
  return s;
}

// The integral of ω from the root of each component of a spanning forest of X's 1-skeleton: a
// potential φ with d₀ φ = ω on the forest's edges (φ = 0 at each root). Its defect,
// max_e |ω(e) − (d₀ φ)(e)|, is zero exactly when ω is exact; the forest is morphe's.
struct Potential {
  Cochain phi;
  Cochain chords;          // ω − d₀ φ: zero on the forest, the fundamental-cycle integrals on the chords
  double defect = 0;
  std::vector<char> in_tree;
};

inline Potential potential(const CWComplex& X, const Cochain& omega) {
  if (omega.size() != X.edges.size()) throw std::invalid_argument("physis::discretisation::potential: one value per edge");
  const auto G = X.one_skeleton();
  const auto F = morphe::combinatorics::spanning_forest(G);
  Potential P;
  P.phi.assign(X.vertices, 0.0);
  P.in_tree = F.in_tree;
  // Vertices in order of depth, so that each parent is set before its children.
  std::vector<std::size_t> order(X.vertices);
  for (std::size_t v = 0; v < X.vertices; v++) order[v] = v;
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return F.depth[a] < F.depth[b]; });
  for (std::size_t v : order) {
    const std::size_t e = F.parent_edge[v];
    if (e == morphe::combinatorics::SpanningForest::none) continue;
    // e joins v to its parent: φ(v) = φ(parent) ± ω(e).
    P.phi[v] = X.edges[e][1] == v ? P.phi[X.edges[e][0]] + omega[e] : P.phi[X.edges[e][1]] - omega[e];
  }
  P.chords = omega;
  const Cochain dphi = d0(X, P.phi);
  for (std::size_t e = 0; e < X.edges.size(); e++) {
    P.chords[e] -= dphi[e];
    P.defect = std::max(P.defect, std::fabs(P.chords[e]));
  }
  return P;
}

// The number of components of X's 1-skeleton, and the component of each vertex.
inline std::vector<std::size_t> components(const CWComplex& X) {
  return morphe::combinatorics::spanning_forest(X.one_skeleton()).component;
}

}  // namespace physis::discretisation
