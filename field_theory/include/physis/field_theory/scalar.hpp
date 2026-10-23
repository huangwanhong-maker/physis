// A real scalar field on the vertices of a cell complex: free and φ⁴ theory.
//
// The field is a 0-cochain φ: vertices → ℝ of a CW complex X, in lattice units (a = 1,
// ħ = 1; physis decision 0005, §2.4). Its Euclidean action is
//
//   S[φ] = ½ Σ_e (φ_head(e) - φ_tail(e))² + Σ_x (½ m² φ_x² + (λ/4!) φ_x⁴)
//        = ½ φᵀ (Δ + m²) φ + (λ/4!) Σ_x φ_x⁴,
//
// the kinetic term on the 1-cells (a discrete exterior derivative squared: ½ |dφ|², Desbrun,
// Hirani, Leok and Marsden 2005, arXiv:math/0508341, §5) and Δ = dᵀd the graph Laplacian of
// the 1-skeleton (loops contribute nothing, parallel edges each once). On the periodic
// hypercubic lattice this is the standard lattice action of Montvay and Münster 1994,
// Quantum Fields on a Lattice, eqs. (1.120) and (2.1)–(2.4), up to the normalisation of φ.
//
// The free field (λ = 0, m² > 0) is Gaussian, with two-point function
// ⟨φ_x φ_y⟩ = G(x, y) = ((Δ + m²)⁻¹)_{xy}, which on the periodic lattice ℤ_{L_1} × ... × ℤ_{L_d}
// is the momentum sum (Montvay and Münster, eq. 1.124)
//
//   G(x, y) = (1/V) Σ_p cos(p·(x - y)) / (m² + Σ_μ 4 sin²(p_μ/2)),  p_μ = 2π n_μ / L_μ.
//
// Updates of one site x at a time leave e^{-S} invariant:
// - Metropolis (Metropolis et al. 1953): φ_x → φ_x + δ, δ uniform in [-step, step], accepted
//   with probability min(1, e^{-ΔS});
// - the heat bath of the free field: given its neighbours, φ_x is Gaussian with mean
//   Σ_{y~x} φ_y / (deg x + m²) and variance 1 / (deg x + m²), and is drawn exactly.
// Observables per configuration: M = Σ_x φ_x, ⟨φ²⟩ = (1/V) Σ φ_x², and from a series of them
// the susceptibility χ = ⟨M²⟩/V (the symmetric phase, ⟨M⟩ = 0) and Binder's cumulant
// U = 1 - ⟨M⁴⟩/(3⟨M²⟩²) (Binder 1981, Z. Phys. B 43, 119), 0 for a Gaussian.

#pragma once

#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <vector>

#include <calx/linalg/dense.hpp>
#include <calx/random/stream.hpp>
#include <morphe/topology/cw_complex.hpp>

namespace physis::field_theory {

using morphe::topology::CWComplex;

struct ScalarTheory {
  double m2 = 1;      // the squared mass, lattice units
  double lambda = 0;  // the quartic coupling (λ/4!) φ⁴
};

// The neighbours of each vertex along the 1-cells, with multiplicity, loops omitted.
struct Neighbours {
  std::vector<std::vector<std::size_t>> of;
  std::size_t size() const { return of.size(); }
};

inline Neighbours neighbours(const CWComplex& X) {
  Neighbours N;
  N.of.resize(X.vertices);
  for (const auto& e : X.edges)
    if (e[0] != e[1]) N.of[e[0]].push_back(e[1]), N.of[e[1]].push_back(e[0]);
  return N;
}

// The graph Laplacian Δ = dᵀd of the 1-skeleton, dense.
inline calx::linalg::Matrix laplacian(const CWComplex& X) {
  calx::linalg::Matrix L(X.vertices);
  for (const auto& e : X.edges)
    if (e[0] != e[1]) {
      L(e[0], e[0]) += 1, L(e[1], e[1]) += 1;
      L(e[0], e[1]) -= 1, L(e[1], e[0]) -= 1;
    }
  return L;
}

// The free propagator (Δ + m²)⁻¹ by dense Gaussian elimination: O(V³), for small complexes.
inline calx::linalg::Matrix free_propagator(const CWComplex& X, double m2) {
  if (!(m2 > 0)) throw std::invalid_argument("physis::field_theory::free_propagator: m² > 0");
  calx::linalg::Matrix A = laplacian(X);
  for (std::size_t i = 0; i < A.n; i++) A(i, i) += m2;
  return calx::linalg::solve(A, calx::linalg::Matrix::identity(A.n));
}

// The coordinates of vertex v of morphe's torus_lattice(sides): row-major, the last side fastest.
inline std::vector<std::size_t> lattice_coordinates(std::size_t v, const std::vector<std::size_t>& sides) {
  std::vector<std::size_t> x(sides.size());
  for (std::size_t mu = sides.size(); mu-- > 0;) x[mu] = v % sides[mu], v /= sides[mu];
  return x;
}

// G(x, y) on the periodic lattice by the momentum sum, for the displacement r = x - y.
inline double lattice_propagator(const std::vector<std::size_t>& sides, double m2, const std::vector<long long>& r) {
  std::size_t V = 1;
  for (auto L : sides) V *= L;
  double G = 0;
  std::vector<std::size_t> n(sides.size(), 0);
  for (std::size_t k = 0; k < V; k++) {
    double phase = 0, denom = m2;
    for (std::size_t mu = 0; mu < sides.size(); mu++) {
      const double p = 2 * std::numbers::pi * static_cast<double>(n[mu]) / static_cast<double>(sides[mu]);
      phase += p * static_cast<double>(r[mu]);
      const double s = std::sin(p / 2);
      denom += 4 * s * s;
    }
    G += std::cos(phase) / denom;
    for (std::size_t mu = sides.size(); mu-- > 0;) {
      if (++n[mu] < sides[mu]) break;
      n[mu] = 0;
    }
  }
  return G / static_cast<double>(V);
}

inline double action(const ScalarTheory& T, const CWComplex& X, const std::vector<double>& phi) {
  double s = 0;
  for (const auto& e : X.edges) {
    const double d = phi[e[1]] - phi[e[0]];
    s += 0.5 * d * d;
  }
  for (double f : phi) s += 0.5 * T.m2 * f * f + T.lambda / 24 * f * f * f * f;
  return s;
}

// The part of the action that depends on φ_x, at the value f.
inline double local_action(const ScalarTheory& T, const Neighbours& N, const std::vector<double>& phi, std::size_t x, double f) {
  double s = 0.5 * T.m2 * f * f + T.lambda / 24 * f * f * f * f;
  for (std::size_t y : N.of[x]) s += 0.5 * (f - phi[y]) * (f - phi[y]);
  return s;
}

// One Metropolis sweep over the vertices in order; returns the number of accepted moves.
inline std::size_t metropolis_sweep(const ScalarTheory& T, const Neighbours& N, std::vector<double>& phi, calx::random::Stream& z,
                                    double step = 1.0) {
  std::size_t accepted = 0;
  for (std::size_t x = 0; x < N.size(); x++) {
    const double old = phi[x], proposal = old + step * (2 * z.uniform() - 1);
    const double dS = local_action(T, N, phi, x, proposal) - local_action(T, N, phi, x, old);
    if (dS <= 0 || z.uniform() < std::exp(-dS)) phi[x] = proposal, accepted++;
  }
  return accepted;
}

// One exact heat-bath sweep of the free field (λ = 0).
inline void heat_bath_sweep(const ScalarTheory& T, const Neighbours& N, std::vector<double>& phi, calx::random::Stream& z) {
  if (T.lambda != 0) throw std::invalid_argument("physis::field_theory::heat_bath_sweep: the free field only (λ = 0)");
  for (std::size_t x = 0; x < N.size(); x++) {
    double s = 0;
    for (std::size_t y : N.of[x]) s += phi[y];
    const double a = static_cast<double>(N.of[x].size()) + T.m2;
    phi[x] = s / a + z.normal() / std::sqrt(a);
  }
}

inline double magnetisation(const std::vector<double>& phi) {
  double m = 0;
  for (double f : phi) m += f;
  return m;
}
inline double mean_square(const std::vector<double>& phi) {
  double s = 0;
  for (double f : phi) s += f * f;
  return s / static_cast<double>(phi.size());
}

// Binder's cumulant from the means of M² and M⁴.
inline double binder_cumulant(double m2, double m4) { return 1 - m4 / (3 * m2 * m2); }

}  // namespace physis::field_theory
