// Matrix models and random surfaces: the gluings of a polygon and the Harer–Zagier numbers.
//
// 't Hooft (1974, Nucl. Phys. B 72, 461) saw that the Feynman graphs of an N × N matrix field
// are ribbon graphs, drawn on surfaces, a graph of genus g weighted N^{χ}: the 1/N expansion
// is a sum over topologies, and in the double-scaling limit the sum over triangulated
// surfaces of two-dimensional quantum gravity (Di Francesco, Ginsparg and Zinn-Justin 1995,
// Phys. Rep. 254, 1). It is the group field theory of dimension 2 with trivial group, the
// first rung of the tensor models (tensor_models.hpp).
//
// The simplest instance: in the Gaussian unitary ensemble (Hermitian M, density ∝ e^{−tr M²/2},
// so ⟨M_ij M_kl⟩ = δ_il δ_jk), Wick's theorem expands ⟨tr M^{2n}⟩ as a sum over the pairings
// of the 2n sides of a 2n-gon, each glued orientably into a closed surface whose vertices are
// the index sums: ⟨tr M^{2n}⟩ = Σ_g ε_g(n) N^{n + 1 − 2g}, ε_g(n) the number of gluings of
// genus g. A pairing π glues side i to side π(i) reversed, so corner i is identified with
// corner π(i) + 1: the vertices are the cycles of γ π, γ the rotation i ↦ i + 1, and
// χ = V − n + 1 = 2 − 2g. Harer and Zagier (1986, Invent. Math. 85, 457–485) found the
// recursion
//   (n + 1) ε_g(n) = 2(2n − 1) ε_g(n − 1) + (n − 1)(2n − 1)(2n − 3) ε_{g−1}(n − 2),
// with ε_0(n) the Catalan number; it computes the virtual Euler characteristics of the moduli
// spaces of curves.

#pragma once

#include <complex>
#include <cstddef>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include <calx/random/stream.hpp>
#include <morphe/arithmetic/integer.hpp>
#include <morphe/topology/constructions.hpp>

namespace physis::quantum_gravity {

using Pairing = std::vector<std::size_t>;  // a fixed-point-free involution of {0..2n-1}

// Calls f(π) for every pairing of 2n sides: (2n − 1)!! of them.
template <class F>
void for_each_pairing(std::size_t n, F&& f) {
  Pairing pi(2 * n, SIZE_MAX);
  std::function<void()> go = [&]() {
    std::size_t i = 0;
    while (i < 2 * n && pi[i] != SIZE_MAX) i++;
    if (i == 2 * n) { f(static_cast<const Pairing&>(pi)); return; }
    for (std::size_t j = i + 1; j < 2 * n; j++)
      if (pi[j] == SIZE_MAX) {
        pi[i] = j, pi[j] = i;
        go();
        pi[i] = pi[j] = SIZE_MAX;
      }
  };
  go();
}

// The number of vertices of the glued surface: the cycles of γ π.
inline std::size_t gluing_vertices(const Pairing& pi) {
  const std::size_t m = pi.size();
  std::vector<char> seen(m, 0);
  std::size_t cycles = 0;
  for (std::size_t s = 0; s < m; s++) {
    if (seen[s]) continue;
    cycles++;
    for (std::size_t c = s; !seen[c]; c = (pi[c] + 1) % m) seen[c] = 1;
  }
  return cycles;
}
inline std::size_t gluing_genus(const Pairing& pi) {
  const std::size_t n = pi.size() / 2;
  return (n + 1 - gluing_vertices(pi)) / 2;
}

// The glued surface as a CW complex (morphe's polygon complex of the word a b ... A ...).
inline morphe::topology::CWComplex gluing_complex(const Pairing& pi) {
  morphe::topology::Polygon P;
  std::map<std::size_t, std::string> label;
  for (std::size_t i = 0; i < pi.size(); i++) {
    const std::size_t first = std::min(i, pi[i]);
    if (!label.count(first)) label[first] = "x" + std::to_string(label.size() + 1);
    P.push_back({label[first], i == first ? 1 : -1});
  }
  return morphe::topology::polygon_complex(std::vector<morphe::topology::Polygon>{P}, "gluing");
}

// ε_g(n) for all g by enumerating the pairings.
inline std::vector<morphe::arithmetic::Integer> gluings_by_genus(std::size_t n) {
  std::vector<morphe::arithmetic::Integer> eps(n / 2 + 1, 0);
  for_each_pairing(n, [&](const Pairing& pi) { eps[gluing_genus(pi)] += 1; });
  return eps;
}

// ε_g(n) by the Harer–Zagier recursion.
inline morphe::arithmetic::Integer harer_zagier(long long g, long long n) {
  using morphe::arithmetic::Integer;
  if (g < 0 || n < 0 || 2 * g > n) return 0;
  if (n == 0) return g == 0 ? 1 : 0;
  static thread_local std::map<std::pair<long long, long long>, Integer> memo;
  if (auto it = memo.find({g, n}); it != memo.end()) return it->second;
  const Integer r = morphe::arithmetic::exact_div(
      Integer(2 * (2 * n - 1)) * harer_zagier(g, n - 1) + Integer((n - 1) * (2 * n - 1) * (2 * n - 3)) * harer_zagier(g - 1, n - 2), Integer(n + 1));
  memo[{g, n}] = r;
  return r;
}

// ⟨tr M^{2n}⟩ in the GUE of size N: Σ_g ε_g(n) N^{n + 1 − 2g}.
inline morphe::arithmetic::Integer gue_moment(long long n, long long N) {
  morphe::arithmetic::Integer m = 0;
  for (long long g = 0; 2 * g <= n; g++) m += harer_zagier(g, n) * morphe::arithmetic::pow(morphe::arithmetic::Integer(N), static_cast<unsigned long long>(n + 1 - 2 * g));
  return m;
}

// A sample of the GUE of size N: M_ii ~ N(0, 1), M_ij = (x + i y)/√2 above the diagonal.
inline std::vector<std::complex<double>> gue_sample(std::size_t N, calx::random::Stream& z) {
  std::vector<std::complex<double>> M(N * N);
  for (std::size_t i = 0; i < N; i++) {
    M[i * N + i] = z.normal();
    for (std::size_t j = i + 1; j < N; j++) {
      const std::complex<double> x(z.normal() / std::sqrt(2.0), z.normal() / std::sqrt(2.0));
      M[i * N + j] = x;
      M[j * N + i] = std::conj(x);
    }
  }
  return M;
}

// tr M^k.
inline double trace_power(const std::vector<std::complex<double>>& M, std::size_t N, std::size_t k) {
  std::vector<std::complex<double>> P(N * N, 0.0), Q(N * N);
  for (std::size_t i = 0; i < N; i++) P[i * N + i] = 1;
  for (std::size_t s = 0; s < k; s++) {
    for (std::size_t i = 0; i < N; i++)
      for (std::size_t j = 0; j < N; j++) {
        std::complex<double> v = 0;
        for (std::size_t l = 0; l < N; l++) v += P[i * N + l] * M[l * N + j];
        Q[i * N + j] = v;
      }
    P.swap(Q);
  }
  double t = 0;
  for (std::size_t i = 0; i < N; i++) t += P[i * N + i].real();
  return t;
}

}  // namespace physis::quantum_gravity
