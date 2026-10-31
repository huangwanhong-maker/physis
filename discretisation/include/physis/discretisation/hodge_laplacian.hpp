// The Hodge Laplacians of every degree on a finite cell complex, with diagonal Hodge stars: their matrices,
// spectra and harmonic spaces; and the reading of a constitutive law as a Hodge star, which makes diffusion
// on cells a Laplacian of the complex (decision 0018).
//
// The operators (Eckmann 1944; Dodziuk 1976, Amer. J. Math. 98, 79; Lim 2020, SIAM Rev. 62, 685). On the
// real cochains C^p of a complex with cellular coboundary d_p = ∂_{p+1}ᵀ : C^p → C^{p+1} (morphe's chain
// complex of any dimension: a CWComplex, a RegularComplex, a simplicial complex) and inner products
// ⟨a, b⟩_p = Σ_σ w_σ a_σ b_σ with weights w > 0 (diagonal Hodge stars ⋆_p = W_p), the adjoint is
// d_p* = W_p⁻¹ d_pᵀ W_{p+1} and
//
//   Δ_p = d_{p−1} d_{p−1}* + d_p* d_p,   ⟨φ, Δ_p φ⟩_p = ‖d_{p−1}* φ‖² + ‖d_p φ‖²  ≥ 0,
//
// self-adjoint for ⟨·,·⟩_p. Its spectrum is that of the symmetric W_p^{1/2} Δ_p W_p^{−1/2}; the eigenvectors
// returned, u_k = W_p^{−1/2} v_k, are orthonormal for ⟨·,·⟩_p. ker Δ_p (the harmonic p-cochains) has dimension
// the Betti number b_p(X; ℝ) whatever the weights: Hodge's theorem for a finite complex. An eigenvalue is a
// frequency of variation over the complex relative to the chosen weights — a wavelength only when the
// weights carry the geometry.
//
// Constitutive laws as Hodge stars (Tonti 1975; Bossavit 1998, Computational Electromagnetism, ch. 5;
// Desbrun, Hirani, Leok and Marsden 2005, §2). Heat in the top cells of a complex — a room's heat content
// Q_c, an extensive p-cochain on the p-cells — conducted through their shared facets of conductance G_f, with
// capacities C_c and the outside at reference temperature, obeys dQ/dt = −d G dᵀ (C⁻¹ Q) with d = d_{p−1}:
// exactly dQ/dt = −Δ_p Q for the stars w_f = 1/G_f on the facets and w_c = 1/C_c on the cells (d_p = 0 on
// top cells, so Δ_p = d_{p−1} d_{p−1}*). A facet bounding one cell only (an outer wall) conducts to the
// outside; its row of dᵀ has a single entry. Walkers or tracers on the vertices of a network diffusing
// through its edges are the case p = 0 with the roles of the stars exchanged: dq/dt = −Δ₀ q for the
// potential q = C⁻¹ n, w_v = C_v, w_e = G_e. Conservation: ⟨1, Δ_p Q⟩ over the cells is the flux through the
// facets on the boundary alone, interior facets cancelling (each carries [c:f] = −[c′:f]).
//
// Dense: O(n_p³) for the spectrum; a few thousand cells. Checks: tests/check_hodge_laplacian.cpp.

#pragma once

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <calx/linalg/dense.hpp>
#include <calx/linalg/eigen.hpp>
#include <calx/linalg/rectangular.hpp>
#include <morphe/algebra/chain_complex.hpp>

namespace physis::discretisation {

// d_p as an n_{p+1} × n_p matrix: (d_p φ)(τ) = Σ_σ [τ : σ] φ(σ).
inline calx::linalg::Rectangular coboundary(const morphe::algebra::ChainComplex& C, std::size_t p) {
  const std::size_t np = p < C.n.size() ? C.n[p] : 0, nq = p + 1 < C.n.size() ? C.n[p + 1] : 0;
  calx::linalg::Rectangular d(nq, np);
  if (nq == 0) return d;
  const auto& B = C.boundary[p + 1];
  for (std::size_t t = 0; t < nq; t++)
    for (const auto& e : B.columns[t]) d(t, e.index) = static_cast<double>(e.value.to_int64());
  return d;
}

struct HodgeLaplacian {
  std::size_t degree = 0;
  calx::linalg::Matrix laplacian;  // Δ_p acting on p-cochains (column vectors)
  std::vector<double> values;      // ascending
  calx::linalg::Matrix vectors;    // column k: the eigenvector of values[k], ⟨u_k, u_l⟩_p = δ_kl
  std::vector<double> weight;      // the star w_p used
  // The number of eigenvalues below tol times the largest: the harmonic dimension.
  std::size_t harmonic(double tol = 1e-9) const {
    const double top = values.empty() ? 0 : std::fabs(values.back());
    std::size_t k = 0;
    for (double v : values) k += std::fabs(v) <= tol * std::max(top, 1.0);
    return k;
  }
};

// Δ_p with stars weights[q] (empty for all 1) on q-cochains, q = p − 1, p, p + 1.
inline HodgeLaplacian hodge_laplacian(const morphe::algebra::ChainComplex& C, std::size_t p, const std::vector<std::vector<double>>& weights = {}) {
  if (p >= C.n.size()) throw std::invalid_argument("physis::discretisation::hodge_laplacian: no cochains of that degree");
  const auto star = [&](std::size_t q) {
    const std::size_t n = q < C.n.size() ? C.n[q] : 0;
    std::vector<double> w(n, 1.0);
    if (q < weights.size() && !weights[q].empty()) {
      if (weights[q].size() != n) throw std::invalid_argument("physis::discretisation::hodge_laplacian: one weight for each cell");
      for (std::size_t i = 0; i < n; i++) {
        if (!(weights[q][i] > 0)) throw std::invalid_argument("physis::discretisation::hodge_laplacian: weights > 0");
        w[i] = weights[q][i];
      }
    }
    return w;
  };
  const std::size_t n = C.n[p];
  const auto wp = star(p);
  HodgeLaplacian H;
  H.degree = p;
  H.weight = wp;
  H.laplacian = calx::linalg::Matrix(n);
  // Up part: d_p* d_p = W_p⁻¹ d_pᵀ W_{p+1} d_p.
  if (p + 1 < C.n.size()) {
    const auto d = coboundary(C, p);
    const auto w1 = star(p + 1);
    for (std::size_t t = 0; t < d.rows; t++)
      for (std::size_t i = 0; i < n; i++) {
        const double a = d(t, i);
        if (a == 0) continue;
        for (std::size_t j = 0; j < n; j++)
          if (d(t, j) != 0) H.laplacian(i, j) += a * w1[t] * d(t, j) / wp[i];
      }
  }
  // Down part: d_{p−1} d_{p−1}* = d_{p−1} W_{p−1}⁻¹ d_{p−1}ᵀ W_p.
  if (p > 0) {
    const auto d = coboundary(C, p - 1);  // n_p × n_{p−1}
    const auto w0 = star(p - 1);
    for (std::size_t s = 0; s < d.cols; s++)
      for (std::size_t i = 0; i < n; i++) {
        const double a = d(i, s);
        if (a == 0) continue;
        for (std::size_t j = 0; j < n; j++)
          if (d(j, s) != 0) H.laplacian(i, j) += a * d(j, s) * wp[j] / w0[s];
      }
  }
  // The spectrum of S = W^{1/2} Δ W^{−1/2}, symmetric.
  calx::linalg::Matrix S(n);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) S(i, j) = std::sqrt(wp[i]) * H.laplacian(i, j) / std::sqrt(wp[j]);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = i + 1; j < n; j++) S(i, j) = S(j, i) = 0.5 * (S(i, j) + S(j, i));
  const auto E = calx::linalg::symmetric_eigen(S);
  H.values = E.values;
  H.vectors = calx::linalg::Matrix(n);
  for (std::size_t k = 0; k < n; k++)
    for (std::size_t i = 0; i < n; i++) H.vectors(i, k) = E.vectors(i, k) / std::sqrt(wp[i]);
  return H;
}

}  // namespace physis::discretisation
