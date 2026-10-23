// Exact results of lattice gauge theory on 2-complexes.
//
// 1. The character expansion (Drouffe and Zuber 1983, Phys. Rep. 102, 1, §2; Montvay and
//    Münster 1994, §3.4). The Boltzmann weight of a plaquette is a class function and so a sum
//    of characters, e^{β Re tr U / N} = Σ_R c_R(β) χ_R(U), with c_R = ∫ dU e^{β Re tr U/N} χ_R(U)*
//    (Haar measure of mass 1), computed here by the Weyl integration formula (physis
//    symmetry). For U(1), c_n = I_n(β); for SU(2), c_j = 2(2j+1) I_{2j+1}(β)/β; for ℤ_N,
//    c_n = Σ_{m ≡ n (N)} I_m(β).
//
// 2. Two-dimensional Yang–Mills on a closed orientable surface (Migdal 1975, Sov. Phys. JETP
//    42, 413; Witten 1991, Commun. Math. Phys. 141, 153; Rusakov 1990). Integrating the links
//    of a cellulation with V vertices, E edges and F 2-cells, by ∫ dU χ_R(A U) χ_S(U⁻¹ B) =
//    δ_RS χ_R(A B)/d_R and ∫ dU χ_R(U A U⁻¹ B) = χ_R(A) χ_R(B)/d_R,
//      Z = ∫ Π_e dU_e Π_f e^{β Re tr U_f/N} = Σ_R d_R^{χ} (c_R/d_R)^F,  χ = V - E + F,
//    and the mean plaquette is (1/F) ∂_β log Z = Σ_R w_R (c_R'/c_R) / Σ_R w_R with
//    w_R = d_R^χ (c_R/d_R)^F. The result depends on the surface only through χ and F:
//    two-dimensional Yang–Mills is invariant under subdivision of the area.
//
// 3. On a disc (a planar complex with one boundary circle) the holonomies of the 2-cells are
//    independent Haar variables after fixing a gauge on a maximal tree, so the mean plaquette
//    is c_0'/c_0 and a Wilson loop in the fundamental representation around A 2-cells is
//    (c_f / (d_f c_0))^A: the area law with string tension -log(c_f/(d_f c_0)) (Wilson 1974).
//
// 4. Flat connections of a finite group G: the assignments of elements to edges under which
//    every 2-cell's holonomy is 1. Gauge-fixing a spanning tree of each component, they
//    correspond to homomorphisms of its fundamental group into G, so their number is
//    Π_components |G|^{V_c - 1} |Hom(π_1(X_c), G)|; at β → ∞ the ℤ_N or finite-group theory
//    counts them, and they are the states of BF theory (Dijkgraaf and Witten 1990).

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

#include <morphe/algebra/finite_group.hpp>
#include <morphe/topology/invariants.hpp>
#include <physis/symmetry/weyl.hpp>

namespace physis::gauge {

using morphe::topology::CWComplex;

struct Representation {
  int label = 0;          // n for U(1) and ℤ_N, two_j for SU(2)
  double dimension = 1;   // d_R
  double c = 0, dc = 0;   // c_R(β) and dc_R/dβ
};

struct CharacterExpansion {
  std::string group;
  double beta = 0;
  std::vector<Representation> reps;
  std::size_t trivial = 0, fundamental = 0;  // indices into reps
};

// U(1): n = -nmax..nmax, nmax = ⌈β⌉ + 40 by default (c_n = I_n(β) falls faster than (β/2)^n/n!).
inline CharacterExpansion expansion_u1(double beta, int nmax = -1) {
  if (nmax < 0) nmax = static_cast<int>(std::ceil(beta)) + 40;
  CharacterExpansion X{"U(1)", beta, {}, 0, 0};
  for (int n = -nmax; n <= nmax; n++) {
    const double c = symmetry::integrate_u1([&](double t) { return std::exp(beta * std::cos(t)) * std::cos(n * t); });
    const double dc = symmetry::integrate_u1([&](double t) { return std::cos(t) * std::exp(beta * std::cos(t)) * std::cos(n * t); });
    if (n == 0) X.trivial = X.reps.size();
    if (n == 1) X.fundamental = X.reps.size();
    X.reps.push_back({n, 1, c, dc});
  }
  return X;
}

// SU(2): two_j = 0..jmax, weight e^{β Re tr U / 2} = e^{β cos(θ/2)}.
inline CharacterExpansion expansion_su2(double beta, int two_jmax = -1) {
  if (two_jmax < 0) two_jmax = static_cast<int>(std::ceil(beta)) + 40;
  CharacterExpansion X{"SU(2)", beta, {}, 0, 1};
  for (int tj = 0; tj <= two_jmax; tj++) {
    const double c = symmetry::coefficient_su2(tj, [&](double t) { return std::exp(beta * std::cos(t / 2)); });
    const double dc = symmetry::coefficient_su2(tj, [&](double t) { return std::cos(t / 2) * std::exp(beta * std::cos(t / 2)); });
    X.reps.push_back({tj, static_cast<double>(tj + 1), c, dc});
  }
  return X;
}

// ℤ_N: n = 0..N-1, characters e^{2πink/N}.
inline CharacterExpansion expansion_zn(std::uint32_t N, double beta) {
  CharacterExpansion X{"Z" + std::to_string(N), beta, {}, 0, N > 1 ? 1u : 0u};
  for (std::uint32_t n = 0; n < N; n++) {
    double c = 0, dc = 0;
    for (std::uint32_t k = 0; k < N; k++) {
      const double t = 2 * std::numbers::pi * k / N, w = std::exp(beta * std::cos(t));
      c += w * std::cos(n * t);
      dc += std::cos(t) * w * std::cos(n * t);
    }
    X.reps.push_back({static_cast<int>(n), 1, c / N, dc / N});
  }
  return X;
}

struct ExactResult {
  double log_z = 0;           // log Z, the links' Haar measures of mass 1
  double mean_plaquette = 0;  // ⟨Re tr U_f / N⟩
};

// Migdal and Witten's partition function on a closed orientable surface of Euler
// characteristic chi cellulated by F 2-cells.
inline ExactResult migdal_witten(const CharacterExpansion& X, long long chi, std::size_t F) {
  std::vector<double> logw;
  for (const auto& r : X.reps) {
    if (!(r.c > 0)) { logw.push_back(-std::numeric_limits<double>::infinity()); continue; }
    logw.push_back(static_cast<double>(chi) * std::log(r.dimension) + static_cast<double>(F) * (std::log(r.c) - std::log(r.dimension)));
  }
  const double top = *std::max_element(logw.begin(), logw.end());
  double sum = 0, num = 0;
  for (std::size_t i = 0; i < logw.size(); i++) {
    const double w = std::exp(logw[i] - top);
    sum += w;
    if (w > 0) num += w * X.reps[i].dc / X.reps[i].c;
  }
  return {top + std::log(sum), num / sum};
}

// The disc: independent 2-cells.
inline ExactResult disc(const CharacterExpansion& X, std::size_t F) {
  const auto& r = X.reps[X.trivial];
  return {static_cast<double>(F) * std::log(r.c), r.dc / r.c};
}
// ⟨Re tr W / N⟩ of a fundamental Wilson loop around A 2-cells of a disc.
inline double disc_wilson_loop(const CharacterExpansion& X, std::size_t area) {
  const auto& f = X.reps[X.fundamental];
  return std::pow(f.c / (f.dimension * X.reps[X.trivial].c), static_cast<double>(area));
}

// ⟨Re tr U / N⟩ of one free SU(N) plaquette, N = 2 or 3, by the Weyl integration formula.
inline double free_plaquette_sun(std::size_t N, double beta) {
  if (N == 2) {
    const auto X = expansion_su2(beta, 0);
    return X.reps[0].dc / X.reps[0].c;
  }
  if (N != 3) throw std::invalid_argument("physis::gauge::free_plaquette_sun: N = 2 or 3");
  auto rt = [](double t1, double t2) { return (std::cos(t1) + std::cos(t2) + std::cos(t1 + t2)) / 3; };
  const double num = symmetry::integrate_su3([&](double t1, double t2) { return rt(t1, t2) * std::exp(beta * rt(t1, t2)); });
  const double den = symmetry::integrate_su3([&](double t1, double t2) { return std::exp(beta * rt(t1, t2)); });
  return num / den;
}

// The number of flat connections of a finite group: Π_c |G|^{V_c - 1} |Hom(π_1(X_c), G)|.
inline morphe::arithmetic::Integer count_flat_connections(const CWComplex& X, const morphe::algebra::FiniteGroup& G) {
  using morphe::arithmetic::Integer;
  const auto F = morphe::combinatorics::spanning_forest(X.one_skeleton());
  std::vector<std::size_t> size(F.components(), 0);
  for (std::size_t v = 0; v < X.vertices; v++) size[F.component[v]]++;
  Integer total = 1;
  for (std::size_t c = 0; c < F.components(); c++) {
    const auto pi1 = morphe::topology::fundamental_group(X, F.roots[c]);
    total *= morphe::arithmetic::pow(Integer(static_cast<long long>(G.order())), size[c] - 1) *
             Integer(static_cast<unsigned long long>(morphe::algebra::count_homomorphisms(pi1.presentation, G)));
  }
  return total;
}

// The same by enumerating all |G|^E assignments (small complexes only).
inline std::uint64_t count_flat_connections_by_enumeration(const CWComplex& X, const morphe::algebra::FiniteGroup& G) {
  const std::size_t E = X.edges.size();
  std::vector<std::uint32_t> g(E, 0);
  std::uint64_t count = 0;
  for (;;) {
    bool flat = true;
    for (const auto& f : X.faces) {
      std::uint32_t h = 0;
      for (const auto& s : f.boundary) h = G.mul(h, s.sign > 0 ? g[s.edge] : G.inv(g[s.edge]));
      if (h != 0) { flat = false; break; }
    }
    count += flat;
    std::size_t e = 0;
    while (e < E && ++g[e] == G.order()) g[e++] = 0;
    if (e == E) break;
  }
  return count;
}

}  // namespace physis::gauge
