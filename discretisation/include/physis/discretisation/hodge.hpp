// The Hodge decomposition of 1-cochains on a CW complex of dimension at most 2.
//
// With inner products ⟨a, b⟩₁ = Σ_e w_e a_e b_e on 1-cochains and ⟨a, b⟩₂ = Σ_f w_f a_f b_f on
// 2-cochains (diagonal Hodge stars, weights > 0; all 1 by default), every 1-cochain splits
// orthogonally and uniquely (Eckmann 1944; Dodziuk 1976; Lim 2020, SIAM Rev. 62, 685, Thm 5.3):
//
//   ω = d₀ α + δ₁ γ + h,     δ₁ = W₁⁻¹ d₁ᵀ W₂ the adjoint of d₁,
//
//   exact part d₀ α:     the gradient of a potential; its integral round every closed path is 0;
//   co-exact part δ₁ γ:  curl, circulation round the 2-cells; d₁ (δ₁ γ) is the whole curvature d₁ ω;
//   harmonic part h:     d₁ h = 0 and d₀ᵀ W₁ h = 0; the harmonic space has dimension b₁(X; ℝ), and
//                        h carries the integrals of ω round the cycles no 2-cell fills.
//
// Read on flows (a 1-chain identified with a 1-cochain), the three parts are the potential-driven
// flow, the circulation round faces and the circulation round holes (HodgeRank: Jiang, Lim, Yao
// and Ye 2011, Math. Program. 127, 203; games: Candogan, Menache, Ozdaglar and Parrilo 2011).
//
// The method. α solves d₀ᵀ W₁ d₀ α = d₀ᵀ W₁ ω and, with γ′ = W₂ γ, γ′ solves
// d₁ W₁⁻¹ d₁ᵀ γ′ = d₁ ω: the normal equations of the two orthogonal projections. Both matrices
// are singular, with kernels the locally constant functions (one per component) and the real
// 2-cycles (ker d₁ᵀ, of dimension b₂); the right-hand sides are orthogonal to them.
//
//   dense (small complexes): adding JJᵀ for a basis J of each kernel makes the matrices positive
//     definite without changing the solution in the complement; the bases are morphe's (the
//     components of a spanning forest, and the free generators of H₂(X; ℤ), which span the real
//     2-cycles of a 2-complex), the solves calx's, O(V³ + F³), exact to rounding.
//   sparse (large complexes): calx's conjugate gradients on the singular systems, consistent, from
//     0, with the potential projected off the components. The exact and co-exact parts are unique
//     whatever solution is reached (d₀ kills the constants, d₁ᵀ the 2-cycles), so no basis of the
//     2-cycles is needed; the curl γ is then a solution, not the one orthogonal to the 2-cycles.
//     O(iterations · cells), to a relative residual of 10⁻¹².
//
// The method is chosen by size (dense up to 1500 vertices and 1500 2-cells) unless one is asked
// for. Many cochains are decomposed with one assembly.

#pragma once

#include <physis/discretisation/cochains.hpp>
#include <calx/linalg/conjugate_gradient.hpp>
#include <calx/linalg/dense.hpp>
#include <calx/linalg/sparse.hpp>
#include <morphe/topology/invariants.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace physis::discretisation {

struct Weights {
  Cochain edge;  // w_e > 0, empty for all 1
  Cochain face;  // w_f > 0, empty for all 1
  double e(std::size_t i) const { return edge.empty() ? 1.0 : edge[i]; }
  double f(std::size_t i) const { return face.empty() ? 1.0 : face[i]; }
};

struct HodgeOptions {
  enum Method { automatic, dense, sparse };
  Method method = automatic;
  std::size_t dense_limit = 1500;  // the largest V and F solved densely by `automatic`
  double tolerance = 1e-12;        // the relative residual of the sparse method
};

struct HodgeDecomposition {
  Cochain exact, coexact, harmonic;  // ω = exact + coexact + harmonic
  Cochain potential;                 // α: exact = d₀ α, summing to zero on each component
  Cochain curl;                      // γ: coexact = δ₁ γ (orthogonal to the 2-cycles by the dense method)
  double exact_norm2 = 0, coexact_norm2 = 0, harmonic_norm2 = 0;  // in ⟨·,·⟩₁
  bool sparse = false;               // the method used
  std::size_t iterations = 0;        // of the sparse method, both solves
  double residual = 0;               // the larger relative residual of the sparse method's two solves
};

inline double inner1(const Cochain& a, const Cochain& b, const Weights& W = {}) {
  double s = 0;
  for (std::size_t i = 0; i < a.size(); i++) s += W.e(i) * a[i] * b[i];
  return s;
}

namespace hodge_detail {

// Solutions of A x = b for each b, in batches of at most n right-hand sides (calx's square solve).
inline std::vector<std::vector<double>> solve_many(const calx::linalg::Matrix& A, const std::vector<std::vector<double>>& bs) {
  const std::size_t n = A.n;
  std::vector<std::vector<double>> xs(bs.size(), std::vector<double>(n, 0.0));
  if (n == 0) return xs;
  for (std::size_t start = 0; start < bs.size(); start += n) {
    const std::size_t m = std::min(n, bs.size() - start);
    calx::linalg::Matrix B(n);
    for (std::size_t j = 0; j < m; j++)
      for (std::size_t i = 0; i < n; i++) B(i, j) = bs[start + j][i];
    const auto X = calx::linalg::solve(A, B);
    for (std::size_t j = 0; j < m; j++)
      for (std::size_t i = 0; i < n; i++) xs[start + j][i] = X(i, j);
  }
  return xs;
}

// A basis of the real 2-cycles, ker d₁ᵀ: the free generators of H₂ (no 3-cells: H₂ = Z₂). A 2-cycle
// lives on the 2-cells and the edges of their words, so the homology is taken of that subcomplex
// alone: the rest of the 1-skeleton, which may be far larger (every trade link of the world beside
// the few hundred sides of a map's triangles), changes nothing in dimension 2.
inline std::vector<Cochain> two_cycles(const CWComplex& X) {
  std::vector<Cochain> Z;
  if (X.faces.empty()) return Z;
  constexpr std::size_t none = static_cast<std::size_t>(-1);
  std::vector<std::size_t> edge_index(X.edges.size(), none), vertex_index(X.vertices, none);
  CWComplex Y;
  auto vertex = [&](std::size_t v) {
    if (vertex_index[v] == none) vertex_index[v] = Y.vertices++;
    return vertex_index[v];
  };
  for (const auto& f : X.faces)
    for (const auto& s : f.boundary)
      if (edge_index[s.edge] == none) {
        edge_index[s.edge] = Y.edges.size();
        const std::size_t t = vertex(X.edges[s.edge][0]), h = vertex(X.edges[s.edge][1]);
        Y.edges.push_back({t, h});
      }
  for (const auto& f : X.faces) {
    morphe::topology::Face g;
    g.base = vertex(f.base);
    for (const auto& s : f.boundary) g.boundary.push_back({edge_index[s.edge], s.sign});
    Y.faces.push_back(std::move(g));
  }
  const auto H = morphe::topology::homology(Y, true);
  if (H.size() < 3) return Z;
  for (const auto& g : H[2].free_generators) {
    Cochain z(X.faces.size(), 0.0);
    for (const auto& entry : g) z[entry.index] = entry.value.to_double();
    Z.push_back(std::move(z));
  }
  return Z;
}

// The entries of the two Laplacians: (row, column, value), duplicates adding up.
struct Triplets {
  std::vector<std::uint32_t> r, c;
  std::vector<double> v;
  void add(std::size_t i, std::size_t j, double x) {
    r.push_back(static_cast<std::uint32_t>(i));
    c.push_back(static_cast<std::uint32_t>(j));
    v.push_back(x);
  }
};

inline Triplets vertex_laplacian(const CWComplex& X, const Weights& W) {
  Triplets T;
  for (std::size_t e = 0; e < X.edges.size(); e++) {
    const std::size_t t = X.edges[e][0], h = X.edges[e][1];
    if (t == h) continue;
    const double w = W.e(e);
    T.add(t, t, w), T.add(h, h, w), T.add(t, h, -w), T.add(h, t, -w);
  }
  return T;
}

// d₁ W₁⁻¹ d₁ᵀ: Σ_e (1/w_e) ε_{fe} ε_{ge}, with ε the exponent sums of the words.
inline Triplets face_laplacian(const CWComplex& X, const Weights& W) {
  const std::size_t E = X.edges.size(), F = X.faces.size();
  std::vector<std::vector<std::pair<std::size_t, double>>> faces_of(E);  // edge -> (2-cell, exponent sum)
  for (std::size_t f = 0; f < F; f++) {
    std::vector<std::pair<std::size_t, double>> sums;
    for (const auto& s : X.faces[f].boundary) {
      auto it = std::find_if(sums.begin(), sums.end(), [&](const auto& p) { return p.first == s.edge; });
      if (it == sums.end()) sums.push_back({s.edge, static_cast<double>(s.sign)});
      else it->second += s.sign;
    }
    for (const auto& [e, c] : sums)
      if (c != 0) faces_of[e].push_back({f, c});
  }
  Triplets T;
  for (std::size_t e = 0; e < E; e++)
    for (const auto& [f, a] : faces_of[e])
      for (const auto& [g, b] : faces_of[e]) T.add(f, g, a * b / W.e(e));
  return T;
}

inline calx::linalg::Matrix dense_of(const Triplets& T, std::size_t n) {
  calx::linalg::Matrix A(n);
  for (std::size_t k = 0; k < T.v.size(); k++) A(T.r[k], T.c[k]) += T.v[k];
  return A;
}

struct SparseMatrix {
  calx::linalg::Csr<> A;
  std::vector<double> a;
};

inline SparseMatrix sparse_of(const Triplets& T, std::size_t n) {
  SparseMatrix S;
  std::vector<std::uint64_t> order;
  S.A = calx::linalg::compress<std::uint32_t>(n, n, T.r, T.c, &order);
  S.a.resize(order.size());
  for (std::size_t k = 0; k < order.size(); k++) S.a[k] = T.v[order[k]];
  return S;
}

}  // namespace hodge_detail

// The decompositions of several 1-cochains, with one assembly of the two systems.
inline std::vector<HodgeDecomposition> hodge_decompositions(const CWComplex& X, const std::vector<Cochain>& omegas, const Weights& W = {},
                                                            const HodgeOptions& options = {}) {
  using namespace hodge_detail;
  if (X.dimension() > 2) throw std::invalid_argument("physis::discretisation::hodge_decomposition: complexes of dimension at most 2");
  const std::size_t V = X.vertices, E = X.edges.size(), F = X.faces.size();
  for (const auto& w : omegas)
    if (w.size() != E) throw std::invalid_argument("physis::discretisation::hodge_decomposition: one value per edge");
  const bool sparse = options.method == HodgeOptions::sparse ||
                      (options.method == HodgeOptions::automatic && (V > options.dense_limit || F > options.dense_limit));
  std::vector<HodgeDecomposition> out(omegas.size());
  const auto comp = components(X);
  // The exact parts: d₀ᵀ W₁ d₀ α = d₀ᵀ W₁ ω.
  std::vector<std::vector<double>> rhs;
  for (const auto& omega : omegas) {
    Cochain weighted(E);
    for (std::size_t e = 0; e < E; e++) weighted[e] = W.e(e) * omega[e];
    rhs.push_back(d0_transpose(X, weighted));
  }
  const Triplets L0 = vertex_laplacian(X, W);
  std::vector<std::vector<double>> alphas;
  if (!sparse) {
    calx::linalg::Matrix A = dense_of(L0, V);
    for (std::size_t u = 0; u < V; u++)
      for (std::size_t v = 0; v < V; v++)
        if (comp[u] == comp[v]) A(u, v) += 1;
    alphas = solve_many(A, rhs);
  } else {
    // The kernel: the normalised indicators of the components.
    std::size_t count = 0;
    for (std::size_t c : comp) count = std::max(count, c + 1);
    std::vector<std::vector<double>> kernel(count, std::vector<double>(V, 0.0));
    std::vector<double> size(count, 0.0);
    for (std::size_t v = 0; v < V; v++) size[comp[v]] += 1;
    for (std::size_t v = 0; v < V; v++) kernel[comp[v]][v] = 1 / std::sqrt(size[comp[v]]);
    const SparseMatrix S = sparse_of(L0, V);
    for (std::size_t k = 0; k < omegas.size(); k++) {
      const auto R = calx::linalg::conjugate_gradient(S.A, std::span<const double>(S.a), std::span<const double>(rhs[k]), options.tolerance, 0, kernel);
      out[k].iterations += R.iterations;
      out[k].residual = std::max(out[k].residual, R.relative_residual);
      alphas.push_back(R.x);
    }
  }
  for (std::size_t k = 0; k < omegas.size(); k++) {
    out[k].potential = alphas[k];
    out[k].exact = d0(X, alphas[k]);
    out[k].coexact.assign(E, 0.0);
    out[k].curl.assign(F, 0.0);
    out[k].sparse = sparse;
  }
  // The co-exact parts: d₁ W₁⁻¹ d₁ᵀ γ′ = d₁ ω; γ = W₂⁻¹ γ′.
  if (F) {
    std::vector<std::vector<double>> curls;
    for (const auto& omega : omegas) curls.push_back(d1(X, omega));
    const Triplets L2 = face_laplacian(X, W);
    std::vector<std::vector<double>> gps;
    if (!sparse) {
      calx::linalg::Matrix B = dense_of(L2, F);
      for (const auto& z : two_cycles(X))
        for (std::size_t f = 0; f < F; f++)
          if (z[f] != 0)
            for (std::size_t g = 0; g < F; g++) B(f, g) += z[f] * z[g];
      gps = solve_many(B, curls);
    } else {
      const SparseMatrix S = sparse_of(L2, F);
      for (std::size_t k = 0; k < omegas.size(); k++) {
        const auto R = calx::linalg::conjugate_gradient(S.A, std::span<const double>(S.a), std::span<const double>(curls[k]), options.tolerance);
        out[k].iterations += R.iterations;
        out[k].residual = std::max(out[k].residual, R.relative_residual);
        gps.push_back(R.x);
      }
    }
    for (std::size_t k = 0; k < omegas.size(); k++) {
      const Cochain lifted = d1_transpose(X, gps[k]);
      for (std::size_t e = 0; e < E; e++) out[k].coexact[e] = lifted[e] / W.e(e);
      for (std::size_t f = 0; f < F; f++) out[k].curl[f] = gps[k][f] / W.f(f);
    }
  }
  for (std::size_t k = 0; k < omegas.size(); k++) {
    auto& H = out[k];
    H.harmonic.resize(E);
    for (std::size_t e = 0; e < E; e++) H.harmonic[e] = omegas[k][e] - H.exact[e] - H.coexact[e];
    H.exact_norm2 = inner1(H.exact, H.exact, W);
    H.coexact_norm2 = inner1(H.coexact, H.coexact, W);
    H.harmonic_norm2 = inner1(H.harmonic, H.harmonic, W);
  }
  return out;
}

inline HodgeDecomposition hodge_decomposition(const CWComplex& X, const Cochain& omega, const Weights& W = {}, const HodgeOptions& options = {}) {
  return hodge_decompositions(X, {omega}, W, options).front();
}

// An orthonormal basis (in ⟨·,·⟩₁) of the harmonic 1-cochains: the harmonic parts of the
// indicators of the chords of a spanning forest, which span the harmonic space, orthonormalised
// by Gram–Schmidt with a relative tolerance. Its size is b₁(X; ℝ).
inline std::vector<Cochain> harmonic_basis(const CWComplex& X, const Weights& W = {}, double tolerance = 1e-9, const HodgeOptions& options = {}) {
  const std::size_t E = X.edges.size();
  const auto P = potential(X, Cochain(E, 0.0));
  std::vector<Cochain> chords;
  for (std::size_t e = 0; e < E; e++)
    if (!P.in_tree[e]) {
      Cochain c(E, 0.0);
      c[e] = 1;
      chords.push_back(std::move(c));
    }
  std::vector<Cochain> basis;
  for (auto& D : hodge_decompositions(X, chords, W, options)) {
    Cochain h = std::move(D.harmonic);
    for (int pass = 0; pass < 2; pass++)  // twice, for orthogonality to rounding
      for (const auto& b : basis) {
        const double p = inner1(h, b, W);
        for (std::size_t i = 0; i < E; i++) h[i] -= p * b[i];
      }
    const double n = std::sqrt(inner1(h, h, W));
    if (n <= tolerance) continue;
    for (double& x : h) x /= n;
    basis.push_back(std::move(h));
  }
  return basis;
}

// Flows with conductances: Ohm's and Kirchhoff's laws on a 2-complex (Kirchhoff 1847; Lim 2020, §6).
//
// A flow F is a 1-chain, the amount carried along each edge per unit of time. With conductances
// c_e > 0 its drop is the 1-cochain E = F / c and its energy Σ F_e² / c_e = ⟨E, E⟩₁ with W₁ = c.
// Decomposing E in that inner product and multiplying back by c gives three flows,
//
//   F = c d₀ α + c δ₁ γ + c h,
//
// orthogonal in Σ F G / c, their energies adding up:
//
//   driven flow c d₀ α:  Ohm's law for the potential α; it carries the whole boundary d₀ᵀ F, since
//                        the other two are divergence-free (d₀ᵀ W₁ δ₁ γ = d₀ᵀ W₁ h = 0: Kirchhoff's
//                        current law); it is the flow the boundary alone would drive through the
//                        conductances;
//   local flow c δ₁ γ = d₁ᵀ W₂ γ:  a circulation round each 2-cell, the same amount on every edge of
//                        it (with face weights 1, γ itself);
//   global flow c h:     a circulation round the holes.
//
// The pairing of a flow with a 1-cochain ω needs no metric: Σ F_e ω_e. With ω decomposed in the same
// inner product it splits into three pairings of like parts, Σ driven·ω_exact + Σ local·ω_coexact +
// Σ global·ω_harmonic, the cross terms vanishing by orthogonality. With every c = 1 this is the plain
// decomposition read on flows.
struct FlowDecomposition {
  Cochain driven, local, global;  // three flows adding up to F
  Cochain potential;              // α: driven = c d₀ α, summing to zero on each component
  Cochain circulation;            // W₂ γ: the local flow round each 2-cell
  double energy = 0, driven_energy = 0, local_energy = 0, global_energy = 0;  // Σ F² / c and its parts
  HodgeDecomposition drop;        // the decomposition of E = F / c, with W₁ = c
};

inline FlowDecomposition flow_decomposition(const CWComplex& X, const Cochain& flow, const Cochain& conductance = {}, const Cochain& face_weights = {},
                                            const HodgeOptions& options = {}) {
  const std::size_t E = X.edges.size();
  if (flow.size() != E || (!conductance.empty() && conductance.size() != E))
    throw std::invalid_argument("physis::discretisation::flow_decomposition: a flow and a conductance for each edge");
  Weights W{conductance, face_weights};
  for (std::size_t e = 0; e < W.edge.size(); e++)
    if (!(W.edge[e] > 0)) throw std::invalid_argument("physis::discretisation::flow_decomposition: conductances are positive");
  Cochain drop(E);
  for (std::size_t e = 0; e < E; e++) drop[e] = flow[e] / W.e(e);
  FlowDecomposition D;
  D.drop = hodge_decomposition(X, drop, W, options);
  D.driven.resize(E), D.local.resize(E), D.global.resize(E);
  for (std::size_t e = 0; e < E; e++) {
    D.driven[e] = W.e(e) * D.drop.exact[e];
    D.local[e] = W.e(e) * D.drop.coexact[e];
    D.global[e] = W.e(e) * D.drop.harmonic[e];
  }
  D.potential = D.drop.potential;
  D.circulation.resize(D.drop.curl.size());
  for (std::size_t f = 0; f < D.circulation.size(); f++) D.circulation[f] = W.f(f) * D.drop.curl[f];
  D.energy = inner1(drop, drop, W);
  D.driven_energy = D.drop.exact_norm2;
  D.local_energy = D.drop.coexact_norm2;
  D.global_energy = D.drop.harmonic_norm2;
  return D;
}

}  // namespace physis::discretisation
