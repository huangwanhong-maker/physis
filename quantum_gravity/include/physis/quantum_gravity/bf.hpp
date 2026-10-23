// BF theory with a finite gauge group on a 2-complex: the spin foam of flat connections.
//
// BF theory has the action ∫ tr(B ∧ F); integrating out B forces the curvature F to vanish,
// so the path integral is a sum over flat connections modulo gauge (Horowitz 1989; Baez 1996,
// Spin foam models, Class. Quantum Grav. 15, 1827). On a 2-complex X with a group G it is the
// lattice gauge theory whose plaquette weight is the δ function at the identity:
//
//   Z_BF(X; G) = ∫ Π_e dg_e Π_f δ(g_f),   g_f the holonomy along the attaching word of f,
//
// with the normalised Haar measure dg = 1/|G| Σ_g and δ(g) = |G| [g = 1], so that
// Z_BF = |G|^{F − E} · #flat. A flat connection is determined, up to the gauge freedom at the
// vertices, by the homomorphism π1(X) → G of its holonomies (Hatcher 2002, §1.A; Lyndon and
// Schupp 1977), and per connected component
//
//   #flat = |G|^{V − 1} |Hom(π1(X), G)|,   so   Z_BF = |G|^{χ(X) − 1} |Hom(π1(X), G)|.
//
// The ratio |Hom(π1 M, G)| / |G| is the untwisted Dijkgraaf–Witten invariant of a manifold
// (Dijkgraaf and Witten 1990; Freed and Quinn 1993). On a closed surface of genus g,
// Mednykh's formula gives Z_BF = Σ_ρ (dim ρ)^{χ}, Witten's two-dimensional gauge theory at
// zero area. On the dual 2-skeleton of a triangulated 3-manifold with v vertices, χ of the
// 2-complex is v, and Z_BF = |G|^{v − 1} |Hom(π1 M, G)|: the divergence of the Ponzano–Regge
// model, one factor per vertex, here finite (Freidel and Louapre 2004).
//
// count_flat_brute enumerates all |G|^E assignments (tiny complexes, for the checks);
// count_flat goes through π1 and morphe's homomorphism count, exponential only in the number
// of generators of π1 after simplification.

#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <morphe/algebra/finite_group.hpp>
#include <morphe/arithmetic/rational.hpp>
#include <morphe/topology/cw_complex.hpp>
#include <morphe/topology/invariants.hpp>

namespace physis::quantum_gravity {

using morphe::algebra::FiniteGroup;
using morphe::arithmetic::Integer;
using morphe::arithmetic::Rational;
using morphe::topology::CWComplex;

// The holonomy of a 2-cell under an assignment of group elements to the edges.
inline std::uint32_t holonomy(const CWComplex& X, const FiniteGroup& G, const std::vector<std::uint32_t>& link, std::size_t face) {
  std::uint32_t g = FiniteGroup::identity();
  for (const auto& s : X.faces[face].boundary) g = G.mul(g, s.sign > 0 ? link[s.edge] : G.inv(link[s.edge]));
  return g;
}

// The number of flat connections, by enumeration of all |G|^E assignments.
inline Integer count_flat_brute(const CWComplex& X, const FiniteGroup& G) {
  const std::size_t E = X.edges.size();
  std::vector<std::uint32_t> link(E, 0);
  // Faces checked as soon as their last edge is assigned.
  std::vector<std::vector<std::size_t>> due(E + 1);
  for (std::size_t f = 0; f < X.faces.size(); f++) {
    std::size_t last = 0;
    for (const auto& s : X.faces[f].boundary) last = std::max(last, s.edge + 1);
    due[last].push_back(f);
  }
  for (std::size_t f : due[0])
    if (holonomy(X, G, link, f) != 0) return 0;
  Integer count = 0;
  // Iterative backtracking over edges 0..E-1: next[k] is the next value to try on edge k.
  std::size_t k = 0;
  std::vector<std::uint32_t> next(E + 1, 0);
  for (;;) {
    if (k == E) {
      count += 1;
      if (k == 0) break;
      k--;
      continue;
    }
    if (next[k] >= G.order()) {
      next[k] = 0;
      if (k == 0) break;
      k--;
      continue;
    }
    link[k] = next[k]++;
    bool ok = true;
    for (std::size_t f : due[k + 1])
      if (holonomy(X, G, link, f) != 0) { ok = false; break; }
    if (ok) k++;
  }
  return count;
}

// The number of flat connections: the product over connected components of
// |G|^{V_i − 1} |Hom(π1(X_i), G)|.
inline Integer count_flat(const CWComplex& X, const FiniteGroup& G) {
  const auto F = morphe::combinatorics::spanning_forest(X.one_skeleton());
  Integer total = 1;
  const Integer order(static_cast<long long>(G.order()));
  for (std::size_t c = 0; c < F.components(); c++) {
    std::size_t vertices = 0;
    for (std::size_t v = 0; v < X.vertices; v++) vertices += F.component[v] == c;
    const auto pi1 = morphe::topology::fundamental_group(X, F.roots[c]);
    const auto homs = morphe::algebra::count_homomorphisms(pi1.simplified.presentation, G);
    total *= morphe::arithmetic::pow(order, vertices - 1) * Integer(static_cast<unsigned long long>(homs));
  }
  return total;
}

// Z_BF = |G|^{F − E} #flat, the δ-function state sum with normalised Haar measure.
inline Rational bf_partition_function(const CWComplex& X, const FiniteGroup& G) {
  const long long exponent = static_cast<long long>(X.faces.size()) - static_cast<long long>(X.edges.size());
  return Rational(count_flat(X, G)) * morphe::arithmetic::pow(Rational(Integer(static_cast<long long>(G.order()))), exponent);
}

// The untwisted Dijkgraaf–Witten invariant |Hom(π1, G)| / |G| of a connected complex.
inline Rational dijkgraaf_witten(const CWComplex& X, const FiniteGroup& G) {
  const auto pi1 = morphe::topology::fundamental_group(X);
  return Rational(Integer(static_cast<unsigned long long>(morphe::algebra::count_homomorphisms(pi1.simplified.presentation, G))),
                  Integer(static_cast<long long>(G.order())));
}

}  // namespace physis::quantum_gravity
