// The Turaev–Viro state sum: the Ponzano–Regge model of three-dimensional Euclidean quantum
// gravity, regularised at a root of unity, which is a topological invariant of closed
// 3-manifolds.
//
// Ponzano and Regge (1968) observed that the asymptotics of the 6j symbol of SU(2) is
// cos(S_Regge + π/4) / √(12πV): a sum over spins on the edges of a triangulation of products
// of 6j symbols is a path integral of Regge's discrete gravity. Their sum diverges; Turaev
// and Viro (1992, Topology 31, 865–902) replaced SU(2) by the quantum group U_q(sl₂) at
// q = e^{iπ/r}, which leaves finitely many spins j = 0, ½, …, (r − 2)/2, and proved the sum
// invariant under the Pachner moves, hence a topological invariant TV_r(M). Turaev and Walker
// proved TV_r(M) = |τ_r(M)|², τ_r the Witten–Reshetikhin–Turaev invariant (Turaev 1994,
// Quantum Invariants of Knots and 3-Manifolds, ch. VII; Roberts 1995, Topology 34, 771–787).
//
// The state sum. A colouring assigns to each edge e of a triangulation a spin, here two_j_e
// in 0..r−2, admissible on every triangle (triangle conditions, integral sum, j₁ + j₂ + j₃ ≤
// r − 2). With the quantum 6j symbol in the Racah–Wigner normalisation (physis symmetry,
// quantum.hpp; Kirillov and Reshetikhin 1989), the convention used here is
//
//   TV_r(M) = w^{−2V} Σ_colourings Π_edges (−1)^{2j_e} [2j_e + 1]
//                                  Π_triangles (−1)^{j₁+j₂+j₃} Π_tetrahedra {6j}_q,
//
// w² = Σ_j [2j + 1]² = r / (2 sin²(π/r)), V the number of vertices. Settling the signs: in
// the 2–3 Pachner move the three new triangles about the new edge x and its weight contribute
// (−1)^{2x} (−1)^{3x + Σ T + Σ B}, the one triangle removed (−1)^{Σ E}, and their ratio is
// (−1)^{x + S}, S the sum of the nine outer spins: exactly the sign of the Biedenharn–Elliott
// identity Σ_x (−1)^{S+x} [2x+1] {6j}{6j}{6j} = {6j}{6j}; orthogonality of the 6j symbols,
// with the same edge and triangle signs, gives the 1–4 move with the factor w^{−2}. This is
// the convention of Burton, Maria and Spreer (2018, Discrete Comput. Geom. 59, 365–399),
// whose triangle weight (−1)^{j₁+j₂+j₃} Δ² carries the Δ factors that the Racah–Wigner
// symbols hold here. The checks confirm it by independence of the triangulation and by the
// values on S³ (w^{−2}), S² × S¹ (1) and the lens spaces (|τ_r|²).
//
// Triangulations are taken in the generalised sense of Turaev and Viro (tetrahedra glued
// along faces, edges and vertices possibly identified within a tetrahedron), as a list of
// tetrahedra with their six edges and four triangles, built from a simplicial complex or from
// a 4-coloured graph (the Feynman graph of the Boulatov group field theory, whose amplitude
// with the quantum group is this state sum on the graph's dual triangulation).
//
// Cost: the colourings are enumerated by backtracking, edges ordered tetrahedron by
// tetrahedron so that triangles close early and inadmissible branches are cut; exponential
// in the number of edges in the worst case. The 6j symbols are tabulated once per r.

#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <morphe/topology/coloured_graph.hpp>
#include <morphe/topology/simplicial.hpp>
#include <physis/symmetry/quantum.hpp>

namespace physis::quantum_gravity {

// A triangulation in the generalised sense: tetrahedra with their edges and triangles.
// Edge slots of a tetrahedron with vertices 0..3: 01, 02, 03, 12, 13, 23; triangle slots:
// the face opposite vertex 0, 1, 2, 3.
struct Triangulation3 {
  std::size_t vertices = 0, edges = 0;
  std::vector<std::array<std::size_t, 3>> triangles;  // the edges of each triangle
  std::vector<std::array<std::size_t, 6>> tetrahedra;  // the edges of each tetrahedron
};

// From a pure 3-dimensional simplicial complex.
inline Triangulation3 triangulation(const morphe::topology::SimplicialComplex& K) {
  if (K.dimension() != 3 || !K.is_pure()) throw std::invalid_argument("physis::quantum_gravity::triangulation: a pure 3-dimensional complex");
  Triangulation3 T;
  T.vertices = K.count(0);
  T.edges = K.count(1);
  for (const auto& t : K.simplices[2]) T.triangles.push_back({K.find({t[0], t[1]}), K.find({t[0], t[2]}), K.find({t[1], t[2]})});
  for (const auto& s : K.simplices[3]) {
    std::array<std::size_t, 6> e{};
    int k = 0;
    for (int i = 0; i < 4; i++)
      for (int j = i + 1; j < 4; j++) e[static_cast<std::size_t>(k++)] = K.find({s[static_cast<std::size_t>(i)], s[static_cast<std::size_t>(j)]});
    T.tetrahedra.push_back(e);
  }
  return T;
}

// From a 4-coloured graph (D = 3): a tetrahedron for each vertex, its vertex of colour i
// opposite its facet of colour i; the edge between the tetrahedron's vertices i and j is the
// bicoloured cycle of the other two colours through the graph vertex, the triangle opposite
// vertex l is the edge of colour l, the vertices of the triangulation are the 3-bubbles.
inline Triangulation3 triangulation(const morphe::topology::ColouredGraph& G) {
  G.validate();
  if (G.colours != 4) throw std::invalid_argument("physis::quantum_gravity::triangulation: a 4-coloured graph");
  const std::size_t p = G.p;
  // face[k][l][v]: the index of the (k, l)-cycle through vertex v (k < l).
  std::array<std::array<std::vector<std::size_t>, 4>, 4> face;
  Triangulation3 T;
  for (std::size_t k = 0; k < 4; k++)
    for (std::size_t l = k + 1; l < 4; l++) {
      face[k][l].assign(2 * p, SIZE_MAX);
      for (const auto& cyc : G.faces(k, l)) {
        for (std::size_t b : cyc) {
          face[k][l][b] = T.edges;
          face[k][l][p + G.sigma[k][b]] = T.edges;
        }
        T.edges++;
      }
    }
  auto edge_of = [&](std::size_t k, std::size_t l, std::size_t v) { return k < l ? face[k][l][v] : face[l][k][v]; };
  T.vertices = G.bubbles();
  // Triangles: the graph edges (colour l, black b); edges: the (l, m)-cycles through b.
  for (std::size_t l = 0; l < 4; l++)
    for (std::size_t b = 0; b < p; b++) {
      std::array<std::size_t, 3> tri{};
      std::size_t n = 0;
      for (std::size_t m = 0; m < 4; m++)
        if (m != l) tri[n++] = edge_of(l, m, b);
      T.triangles.push_back(tri);
    }
  for (std::size_t v = 0; v < 2 * p; v++) {
    std::array<std::size_t, 6> e{};
    std::size_t n = 0;
    for (std::size_t i = 0; i < 4; i++)
      for (std::size_t j = i + 1; j < 4; j++) {
        std::size_t k = 0, l = 0, m = 0;
        for (std::size_t c = 0; c < 4; c++)
          if (c != i && c != j) (m++ == 0 ? k : l) = c;
        e[n++] = edge_of(k, l, v);  // the edge {i, j} is the cycle of the complementary colours
      }
    T.tetrahedra.push_back(e);
  }
  return T;
}

struct StateSum {
  double value = 0;            // TV_r(M)
  std::uint64_t colourings = 0;  // admissible colourings visited
  int r = 0;
};

inline StateSum turaev_viro(const Triangulation3& T, int r) {
  if (r < 3) throw std::invalid_argument("physis::quantum_gravity::turaev_viro: r >= 3");
  const int n = r - 1;  // two_j in 0..r-2
  // The 6j table, indexed by (a, b, c, d, e, f) in base n, with the tetrahedron's edge slots
  // 01, 02, 03, 12, 13, 23 mapped to {a b c; d e f} as a = 01, b = 02, c = 12, d = 23, e = 13,
  // f = 03, so that the triads (abc), (aef), (dbf), (dec) are the four triangles.
  std::vector<double> six(static_cast<std::size_t>(n) * n * n * n * n * n);
  auto key = [n](int a, int b, int c, int d, int e, int f) { return ((((static_cast<std::size_t>(a) * n + b) * n + c) * n + d) * n + e) * n + f; };
  for (int a = 0; a < n; a++)
    for (int b = 0; b < n; b++)
      for (int c = 0; c < n; c++)
        for (int d = 0; d < n; d++)
          for (int e = 0; e < n; e++)
            for (int f = 0; f < n; f++) six[key(a, b, c, d, e, f)] = symmetry::quantum_6j(r, a, b, c, d, e, f);
  std::vector<double> edge_weight(static_cast<std::size_t>(n));
  for (int j = 0; j < n; j++) edge_weight[static_cast<std::size_t>(j)] = symmetry::quantum_dimension_signed(j, r);

  // Order the edges tetrahedron by tetrahedron (breadth first through shared edges).
  const std::size_t E = T.edges;
  std::vector<std::size_t> order, position(E, SIZE_MAX);
  std::vector<std::vector<std::size_t>> tets_of_edge(E);
  for (std::size_t t = 0; t < T.tetrahedra.size(); t++)
    for (auto e : T.tetrahedra[t]) tets_of_edge[e].push_back(t);
  std::vector<char> tet_seen(T.tetrahedra.size(), 0);
  for (std::size_t start = 0; start < T.tetrahedra.size(); start++) {
    if (tet_seen[start]) continue;
    std::vector<std::size_t> queue{start};
    tet_seen[start] = 1;
    for (std::size_t qi = 0; qi < queue.size(); qi++) {
      for (auto e : T.tetrahedra[queue[qi]]) {
        if (position[e] != SIZE_MAX) continue;
        position[e] = order.size();
        order.push_back(e);
        for (auto t2 : tets_of_edge[e])
          if (!tet_seen[t2]) tet_seen[t2] = 1, queue.push_back(t2);
      }
    }
  }
  for (std::size_t e = 0; e < E; e++)
    if (position[e] == SIZE_MAX) position[e] = order.size(), order.push_back(e);
  // The triangles and tetrahedra completed when each position is assigned.
  std::vector<std::vector<std::size_t>> tri_done(E), tet_done(E);
  for (std::size_t f = 0; f < T.triangles.size(); f++) {
    std::size_t last = 0;
    for (auto e : T.triangles[f]) last = std::max(last, position[e]);
    tri_done[last].push_back(f);
  }
  for (std::size_t t = 0; t < T.tetrahedra.size(); t++) {
    std::size_t last = 0;
    for (auto e : T.tetrahedra[t]) last = std::max(last, position[e]);
    tet_done[last].push_back(t);
  }
  std::vector<int> colour(E, 0);
  StateSum S;
  S.r = r;
  std::function<void(std::size_t, double)> go = [&](std::size_t k, double weight) {
    if (k == E) {
      S.value += weight;
      S.colourings++;
      return;
    }
    const std::size_t e = order[k];
    for (int j = 0; j < n; j++) {
      colour[e] = j;
      double w = weight * edge_weight[static_cast<std::size_t>(j)];
      bool ok = true;
      for (auto f : tri_done[k]) {
        const auto& tr = T.triangles[f];
        const int a = colour[tr[0]], b = colour[tr[1]], c = colour[tr[2]];
        if (!symmetry::admissible(a, b, c, r)) { ok = false; break; }
        if (((a + b + c) / 2) % 2) w = -w;
      }
      if (!ok) continue;
      for (auto t : tet_done[k]) {
        const auto& x = T.tetrahedra[t];
        w *= six[key(colour[x[0]], colour[x[1]], colour[x[3]], colour[x[5]], colour[x[4]], colour[x[2]])];
      }
      if (w == 0) continue;
      go(k + 1, w);
    }
  };
  go(0, 1.0);
  S.value *= std::pow(symmetry::total_quantum_dimension_squared(r), -static_cast<double>(T.vertices));
  return S;
}

inline StateSum turaev_viro(const morphe::topology::SimplicialComplex& K, int r) { return turaev_viro(triangulation(K), r); }
inline StateSum turaev_viro(const morphe::topology::ColouredGraph& G, int r) { return turaev_viro(triangulation(G), r); }

// |τ_r(L(p, 1))|², the square of the modulus of the Witten–Reshetikhin–Turaev invariant of the
// lens space L(p, 1), the manifold glued from two solid tori by S T^p S: with
// S_{0j} = √(2/r) sin((2j+1)π/r) and the twists θ_j = e^{2πi j(j+1)/r},
// τ_r(L(p, 1)) = Σ_j S_{0j}² θ_j^p up to a phase (Jeffrey 1992, Commun. Math. Phys. 147,
// 563–604; Turaev 1994, ch. VII). By the theorem of Turaev and Walker, TV_r(L(p, 1)) equals it.
inline double lens_space_tv(int p, int r) {
  double re = 0, im = 0;
  for (int two_j = 0; two_j <= r - 2; two_j++) {
    const double s = std::sqrt(2.0 / r) * std::sin((two_j + 1) * std::numbers::pi / r);
    const double h = two_j * (two_j + 2) / 4.0 / r;  // j (j + 1) / r
    re += s * s * std::cos(2 * std::numbers::pi * p * h);
    im += s * s * std::sin(2 * std::numbers::pi * p * h);
  }
  return re * re + im * im;
}

}  // namespace physis::quantum_gravity
