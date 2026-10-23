// Regge geometry: a piecewise-flat spacetime as a simplicial complex with edge lengths.
//
// Regge (1961, Nuovo Cimento 19, 558) replaced a smooth geometry by a simplicial complex
// whose simplices are flat and whose metric is fixed by the lengths of its edges. Curvature
// is concentrated on the hinges, the simplices of codimension 2: around a hinge h in the
// interior, the dihedral angles θ_h of the k-simplices that contain it fall short of 2π by
// the deficit angle
//   ε_h = 2π - Σ_{σ ∋ h} θ_h(σ),
// and the Einstein-Hilbert action ½ ∫ R √g becomes the Regge action
//   S = Σ_h A_h ε_h,
// A_h the (k-2)-volume of the hinge (Regge 1961, eq. 3.4; Misner, Thorne and Wheeler 1973,
// ch. 42). In two dimensions the hinges are vertices (A_h = 1) and Σ_v ε_v = 2π χ is the
// Gauss-Bonnet theorem; S is then topological. The variation of S is, by Schläfli's
// identity Σ_h A_h dθ_h = 0 for each simplex (Schläfli 1858; Milnor 1994), simply
//   ∂S/∂s_e = Σ_h (∂A_h/∂s_e) ε_h,
// so Regge's equations are ∂S/∂s_e = 0. At a boundary hinge the angle π - Σθ (the
// extrinsic curvature of the boundary) enters the boundary term of Hartle and Sorkin (1981).
//
// A simplex's geometry is given by its squared edge lengths s_ij. With vertex 0 as base, its
// Gram matrix G_ab = (s_0a + s_0b - s_ab) / 2 (a, b = 1..k) gives its volume
//   V = sqrt(det G) / k!,
// equal to the Cayley-Menger determinant's (-1)^{k+1} det CM / (2^k (k!)²) = V² (Blumenthal
// 1953, §40), and an embedding in ℝ^k by the Cholesky factor G = L Lᵀ, from which the
// dihedral angle at the hinge opposite the vertices i and j is the angle between the
// components of p_i and p_j orthogonal to the hinge. The derivative of the volume follows
// from ∂ det G = det G tr(G⁻¹ ∂G):
//   ∂V/∂s_ab = -(V/2) (G⁻¹)_ab (a ≠ b, both ≠ 0),   ∂V/∂s_0a = (V/2) Σ_b (G⁻¹)_ab.
//
// Signature: Euclidean here; every simplex must be nondegenerate (G positive definite), or
// the functions throw. Lorentzian Regge calculus (complex angles, timelike and spacelike
// hinges; Sorkin 1975, 2019) is planned. Complexity: a simplex of dimension k costs O(k³);
// the complex, O(F k³) for F facets.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <numbers>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace physis::spacetime {

using SquaredLengths = std::vector<std::vector<double>>;  // s[i][j], symmetric, zero diagonal
using Vertices = std::vector<std::size_t>;                 // the labels of a simplex's vertices

namespace detail {

// Determinant of a small square matrix by elimination with partial pivoting.
inline double det(std::vector<std::vector<double>> a) {
  const std::size_t n = a.size();
  double d = 1;
  for (std::size_t k = 0; k < n; k++) {
    std::size_t p = k;
    for (std::size_t i = k + 1; i < n; i++)
      if (std::fabs(a[i][k]) > std::fabs(a[p][k])) p = i;
    if (a[p][k] == 0) return 0;
    if (p != k) { std::swap(a[p], a[k]); d = -d; }
    d *= a[k][k];
    for (std::size_t i = k + 1; i < n; i++) {
      const double f = a[i][k] / a[k][k];
      for (std::size_t j = k; j < n; j++) a[i][j] -= f * a[k][j];
    }
  }
  return d;
}

inline std::vector<std::vector<double>> inverse(std::vector<std::vector<double>> a) {
  const std::size_t n = a.size();
  std::vector<std::vector<double>> b(n, std::vector<double>(n, 0.0));
  for (std::size_t i = 0; i < n; i++) b[i][i] = 1;
  for (std::size_t k = 0; k < n; k++) {
    std::size_t p = k;
    for (std::size_t i = k + 1; i < n; i++)
      if (std::fabs(a[i][k]) > std::fabs(a[p][k])) p = i;
    if (a[p][k] == 0) throw std::domain_error("physis::spacetime::regge: a degenerate simplex");
    std::swap(a[p], a[k]);
    std::swap(b[p], b[k]);
    const double dk = a[k][k];
    for (std::size_t j = 0; j < n; j++) { a[k][j] /= dk; b[k][j] /= dk; }
    for (std::size_t i = 0; i < n; i++) {
      if (i == k || a[i][k] == 0) continue;
      const double f = a[i][k];
      for (std::size_t j = 0; j < n; j++) { a[i][j] -= f * a[k][j]; b[i][j] -= f * b[k][j]; }
    }
  }
  return b;
}

inline double factorial(std::size_t k) {
  double f = 1;
  for (std::size_t i = 2; i <= k; i++) f *= static_cast<double>(i);
  return f;
}

}  // namespace detail

// The Gram matrix of a simplex with vertex 0 as base, from its squared edge lengths.
inline std::vector<std::vector<double>> gram(const SquaredLengths& s) {
  const std::size_t k = s.size() - 1;
  std::vector<std::vector<double>> G(k, std::vector<double>(k));
  for (std::size_t a = 1; a <= k; a++)
    for (std::size_t b = 1; b <= k; b++) G[a - 1][b - 1] = 0.5 * (s[0][a] + s[0][b] - s[a][b]);
  return G;
}

// The k-volume of a k-simplex (1 for a point).
inline double simplex_volume(const SquaredLengths& s) {
  const std::size_t k = s.size() - 1;
  if (k == 0) return 1;
  const double d = detail::det(gram(s));
  if (!(d > 0)) throw std::domain_error("physis::spacetime::regge: the edge lengths violate the triangle inequalities (or the simplex is flat)");
  return std::sqrt(d) / detail::factorial(k);
}

// The same volume by the Cayley-Menger determinant.
inline double cayley_menger_volume(const SquaredLengths& s) {
  const std::size_t k = s.size() - 1;
  std::vector<std::vector<double>> cm(k + 2, std::vector<double>(k + 2, 1.0));
  cm[0][0] = 0;
  for (std::size_t i = 0; i <= k; i++)
    for (std::size_t j = 0; j <= k; j++) cm[i + 1][j + 1] = s[i][j];
  const double f = detail::factorial(k);
  const double v2 = ((k + 1) % 2 ? -1.0 : 1.0) * detail::det(cm) / (std::pow(2.0, static_cast<double>(k)) * f * f);
  if (!(v2 > 0)) throw std::domain_error("physis::spacetime::regge: the Cayley-Menger determinant has the wrong sign");
  return std::sqrt(v2);
}

// Coordinates of the k+1 vertices in ℝ^k: vertex 0 at the origin, the others from the
// Cholesky factor of the Gram matrix.
inline std::vector<std::vector<double>> embed(const SquaredLengths& s) {
  const auto G = gram(s);
  const std::size_t k = G.size();
  std::vector<std::vector<double>> L(k, std::vector<double>(k, 0.0));
  for (std::size_t i = 0; i < k; i++)
    for (std::size_t j = 0; j <= i; j++) {
      double sum = G[i][j];
      for (std::size_t m = 0; m < j; m++) sum -= L[i][m] * L[j][m];
      if (i == j) {
        if (!(sum > 0)) throw std::domain_error("physis::spacetime::regge: the Gram matrix is not positive definite");
        L[i][i] = std::sqrt(sum);
      } else {
        L[i][j] = sum / L[j][j];
      }
    }
  std::vector<std::vector<double>> p(k + 1, std::vector<double>(k, 0.0));
  for (std::size_t a = 1; a <= k; a++) p[a] = L[a - 1];
  return p;
}

// The interior dihedral angle of a k-simplex (k >= 2) at its hinge opposite the vertices i
// and j: the angle between its facets opposite i and opposite j.
inline double dihedral_angle(const SquaredLengths& s, std::size_t i, std::size_t j) {
  const std::size_t k = s.size() - 1;
  if (k < 2 || i == j || i > k || j > k) throw std::invalid_argument("physis::spacetime::regge: a dihedral angle needs k >= 2 and two distinct vertices");
  const auto p = embed(s);
  std::vector<std::size_t> hinge;
  for (std::size_t v = 0; v <= k; v++)
    if (v != i && v != j) hinge.push_back(v);
  const auto& b = p[hinge[0]];
  // An orthonormal basis of the hinge's directions (modified Gram-Schmidt).
  std::vector<std::vector<double>> basis;
  for (std::size_t m = 1; m < hinge.size(); m++) {
    std::vector<double> w(k);
    for (std::size_t c = 0; c < k; c++) w[c] = p[hinge[m]][c] - b[c];
    for (const auto& e : basis) {
      double dot = 0;
      for (std::size_t c = 0; c < k; c++) dot += w[c] * e[c];
      for (std::size_t c = 0; c < k; c++) w[c] -= dot * e[c];
    }
    double norm = 0;
    for (double x : w) norm += x * x;
    norm = std::sqrt(norm);
    for (double& x : w) x /= norm;
    basis.push_back(w);
  }
  auto normal_part = [&](std::size_t v) {
    std::vector<double> w(k);
    for (std::size_t c = 0; c < k; c++) w[c] = p[v][c] - b[c];
    for (const auto& e : basis) {
      double dot = 0;
      for (std::size_t c = 0; c < k; c++) dot += w[c] * e[c];
      for (std::size_t c = 0; c < k; c++) w[c] -= dot * e[c];
    }
    return w;
  };
  const auto u = normal_part(i), w = normal_part(j);
  double uu = 0, ww = 0, uw = 0;
  for (std::size_t c = 0; c < k; c++) { uu += u[c] * u[c]; ww += w[c] * w[c]; uw += u[c] * w[c]; }
  return std::atan2(std::sqrt(std::max(0.0, uu * ww - uw * uw)), uw);
}

// ∂V/∂s_ab for every pair a < b, as a symmetric matrix with zero diagonal.
inline std::vector<std::vector<double>> volume_gradient(const SquaredLengths& s) {
  const std::size_t k = s.size() - 1;
  std::vector<std::vector<double>> dV(k + 1, std::vector<double>(k + 1, 0.0));
  if (k == 0) return dV;
  const double V = simplex_volume(s);
  const auto Gi = detail::inverse(gram(s));
  for (std::size_t a = 1; a <= k; a++) {
    double row = 0;
    for (std::size_t b = 1; b <= k; b++) row += Gi[a - 1][b - 1];
    dV[0][a] = dV[a][0] = 0.5 * V * row;
    for (std::size_t b = a + 1; b <= k; b++) dV[a][b] = dV[b][a] = -0.5 * V * Gi[a - 1][b - 1];
  }
  return dV;
}

// Squared edge lengths of the vertices of a simplex from a function s(u, v) of vertex labels.
template <class Lengths>
inline SquaredLengths squared_lengths(const Vertices& simplex, Lengths&& s) {
  const std::size_t m = simplex.size();
  SquaredLengths q(m, std::vector<double>(m, 0.0));
  for (std::size_t i = 0; i < m; i++)
    for (std::size_t j = i + 1; j < m; j++) q[i][j] = q[j][i] = s(simplex[i], simplex[j]);
  return q;
}

// A table of squared edge lengths keyed by unordered pairs of vertex labels.
class EdgeLengths {
 public:
  void set(std::size_t a, std::size_t b, double squared) { table_[key(a, b)] = squared; }
  double operator()(std::size_t a, std::size_t b) const {
    const auto it = table_.find(key(a, b));
    if (it == table_.end()) throw std::out_of_range("physis::spacetime::EdgeLengths: no length for an edge");
    return it->second;
  }
  const std::map<std::pair<std::size_t, std::size_t>, double>& table() const { return table_; }

 private:
  static std::pair<std::size_t, std::size_t> key(std::size_t a, std::size_t b) { return a < b ? std::pair{a, b} : std::pair{b, a}; }
  std::map<std::pair<std::size_t, std::size_t>, double> table_;
};

// A pure simplicial complex of dimension k >= 2, its hinges, and its Regge geometry under
// any assignment of squared edge lengths.
class ReggeComplex {
 public:
  explicit ReggeComplex(std::vector<Vertices> facets) : facets_(std::move(facets)) {
    if (facets_.empty()) throw std::invalid_argument("physis::spacetime::ReggeComplex: no facets");
    k_ = facets_[0].size() - 1;
    if (k_ < 2) throw std::invalid_argument("physis::spacetime::ReggeComplex: the dimension must be at least 2");
    std::map<Vertices, std::size_t> faces;  // codimension-1 faces and their numbers of cofacets
    std::map<Vertices, std::size_t> hinge_index;
    for (std::size_t f = 0; f < facets_.size(); f++) {
      auto& v = facets_[f];
      if (v.size() != k_ + 1) throw std::invalid_argument("physis::spacetime::ReggeComplex: facets of different dimensions");
      std::sort(v.begin(), v.end());
      for (std::size_t i = 0; i <= k_; i++) {
        Vertices face;
        for (std::size_t m = 0; m <= k_; m++)
          if (m != i) face.push_back(v[m]);
        faces[face]++;
      }
      for (std::size_t i = 0; i <= k_; i++)
        for (std::size_t j = i + 1; j <= k_; j++) {
          Vertices h;
          for (std::size_t m = 0; m <= k_; m++)
            if (m != i && m != j) h.push_back(v[m]);
          auto [it, inserted] = hinge_index.try_emplace(h, hinges_.size());
          if (inserted) { hinges_.push_back(h); star_.emplace_back(); }
          star_[it->second].push_back({f, i, j});
        }
      for (std::size_t i = 0; i <= k_; i++)
        for (std::size_t j = i + 1; j <= k_; j++) edges_.emplace(v[i], v[j]);
    }
    interior_.assign(hinges_.size(), true);
    for (std::size_t h = 0; h < hinges_.size(); h++)
      for (const auto& inc : star_[h]) {
        for (std::size_t drop : {inc.i, inc.j}) {
          Vertices face;
          for (std::size_t m = 0; m <= k_; m++)
            if (m != drop) face.push_back(facets_[inc.facet][m]);
          if (faces[face] != 2) interior_[h] = false;
        }
      }
  }

  std::size_t dimension() const { return k_; }
  const std::vector<Vertices>& facets() const { return facets_; }
  const std::vector<Vertices>& hinges() const { return hinges_; }
  bool interior(std::size_t h) const { return interior_[h]; }
  // Every edge (a, b), a < b.
  std::vector<std::pair<std::size_t, std::size_t>> edges() const { return {edges_.begin(), edges_.end()}; }

  // Σ θ over the k-simplices around each hinge.
  template <class Lengths>
  std::vector<double> angle_sums(Lengths&& s) const {
    std::vector<double> sum(hinges_.size(), 0.0);
    for (std::size_t h = 0; h < hinges_.size(); h++)
      for (const auto& inc : star_[h]) sum[h] += dihedral_angle(squared_lengths(facets_[inc.facet], s), inc.i, inc.j);
    return sum;
  }
  // The deficit angle 2π - Σθ at each interior hinge, and π - Σθ at each boundary hinge.
  template <class Lengths>
  std::vector<double> deficits(Lengths&& s) const {
    auto e = angle_sums(s);
    for (std::size_t h = 0; h < e.size(); h++) e[h] = (interior_[h] ? 2 : 1) * std::numbers::pi - e[h];
    return e;
  }
  // The (k-2)-volume of each hinge.
  template <class Lengths>
  std::vector<double> hinge_volumes(Lengths&& s) const {
    std::vector<double> A(hinges_.size());
    for (std::size_t h = 0; h < hinges_.size(); h++) A[h] = simplex_volume(squared_lengths(hinges_[h], s));
    return A;
  }
  // S = Σ_interior A_h ε_h, and with the boundary term Σ_boundary A_h (π - Σθ) when asked.
  template <class Lengths>
  double action(Lengths&& s, bool boundary_term = false) const {
    const auto e = deficits(s);
    const auto A = hinge_volumes(s);
    double S = 0;
    for (std::size_t h = 0; h < e.size(); h++)
      if (interior_[h] || boundary_term) S += A[h] * e[h];
    return S;
  }
  // ∂S/∂s_e for every edge (a, b), a < b, by Schläfli's identity: Σ_h (∂A_h/∂s_e) ε_h.
  template <class Lengths>
  std::map<std::pair<std::size_t, std::size_t>, double> action_gradient(Lengths&& s, bool boundary_term = false) const {
    std::map<std::pair<std::size_t, std::size_t>, double> g;
    for (const auto& e : edges_) g[e] = 0;
    const auto e = deficits(s);
    for (std::size_t h = 0; h < hinges_.size(); h++) {
      if (!interior_[h] && !boundary_term) continue;
      const auto& v = hinges_[h];
      if (v.size() < 2) continue;  // a point: its volume is 1, and does not vary
      const auto dA = volume_gradient(squared_lengths(v, s));
      for (std::size_t i = 0; i < v.size(); i++)
        for (std::size_t j = i + 1; j < v.size(); j++) g[{v[i], v[j]}] += dA[i][j] * e[h];
    }
    return g;
  }

 private:
  struct Incidence {
    std::size_t facet, i, j;  // the facet, and the positions of its two vertices not in the hinge
  };
  std::size_t k_ = 0;
  std::vector<Vertices> facets_, hinges_;
  std::vector<std::vector<Incidence>> star_;
  std::vector<bool> interior_;
  std::set<std::pair<std::size_t, std::size_t>> edges_;
};

}  // namespace physis::spacetime
