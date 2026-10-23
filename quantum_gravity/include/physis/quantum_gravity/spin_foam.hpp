// Spin foams: three-dimensional SU(2) BF theory (Ponzano–Regge) on 2-complexes, and its
// quantum deformation at a root of unity (Turaev–Viro).
//
// A spin foam is a 2-complex whose faces carry spins and whose vertices carry amplitudes
// (Baez 1998, Class. Quantum Grav. 15, 1827–1858; Perez 2013, Living Rev. Relativ. 16, 3,
// §3). Here the 2-complex is dual to a closed triangulated 3-manifold: a spin-foam vertex for
// each tetrahedron, a spin-foam edge for each triangle, a spin-foam face for each edge of the
// triangulation. It is built either from a simplicial complex, as morphe's dual 2-complex, or
// from a 4-coloured graph, the Feynman graph of the Boulatov group field theory, as the
// graph's 2-complex of bicoloured cycles (Boulatov 1992, Mod. Phys. Lett. A 7, 1629–1646;
// Freidel 2005, Int. J. Theor. Phys. 44, 1769–1783). Its spins j_f are those of the edges of
// the triangulation, passed as two_j.
//
// The amplitude of a spin assignment is that of the Turaev–Viro state sum, term by term (see
// turaev_viro.hpp for the signs):
//
//   A(j) = w^{−2V} Π_faces (−1)^{2j_f} [2j_f + 1]_q Π_edges (−1)^{j₁+j₂+j₃} Π_vertices {6j}_q,
//
// a spin-foam edge being a triangle with its three spins. The vertex amplitude is the
// evaluation of the vertex's boundary spin network, the tetrahedral graph whose nodes are the
// four triangles and whose links are the six edges with their spins: the 6j symbol. The sum
// of A over admissible spins is TV_r. At r = 0 the weights are those of Ponzano and Regge
// (1968): [2j + 1] → 2j + 1, the classical 6j symbol, no normalisation (the sum diverges and
// is not taken).
//
// Geometry. In the Ponzano–Regge reading each spin is a length, ℓ = j + ½, and each vertex a
// quantum tetrahedron. When the six lengths bound a Euclidean tetrahedron (positive
// Cayley–Menger determinant), its volume V and exterior dihedral angles θ_e give the
// asymptotics (Ponzano and Regge 1968; Roberts 1999, Geom. Topol. 3, 21–66)
//   {6j} ≈ cos(Σ_e (j_e + ½) θ_e + π/4) / √(12πV),
// the cosine of Regge's action: the semiclassical limit of the spin foam is discrete general
// relativity. The Regge geometry of a whole spin assignment (deficit angles at the edges of
// the triangulation, Regge action Σ ℓ_e ε_e) follows when every tetrahedron is Euclidean.
//
// Sampling. The amplitudes alternate in sign, so Monte Carlo samples spins with weight |A|
// (Metropolis, single-face changes and changes of the spins on all edges at a vertex of the
// triangulation, which together connect every admissible assignment of the same parity class
// in H¹(M; ℤ/2): changes of one spin by one unit break the parity of its triangles), and
// reports ⟨sign A⟩_{|A|} = Σ A / Σ |A|, the factor that reweights to the true state sum. The
// chain is ergodic within the cohomology class of its start (the zero assignment) and exact
// in distribution for |A| there.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <vector>

#include <calx/random/stream.hpp>
#include <morphe/topology/coloured_graph.hpp>
#include <morphe/topology/dual.hpp>
#include <morphe/topology/simplicial.hpp>
#include <physis/spacetime/regge.hpp>
#include <physis/symmetry/quantum.hpp>
#include <physis/symmetry/recoupling.hpp>

#include "turaev_viro.hpp"

namespace physis::quantum_gravity {

// The edge slots of a tetrahedron with vertices 0..3: 01, 02, 03, 12, 13, 23.
inline constexpr std::array<std::array<int, 2>, 6> tetrahedron_slots{{{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}};

struct SpinFoam {
  morphe::topology::CWComplex complex;  // the 2-complex: vertices ↔ tetrahedra, faces ↔ edges of the triangulation
  Triangulation3 triangulation;         // the dual triangulation
  std::vector<std::size_t> face_edge;   // face f ↔ edge face_edge[f] of the triangulation
  std::vector<std::size_t> edge_face;   // its inverse
  std::vector<std::size_t> vertex_tetrahedron;              // spin-foam vertex v ↔ tetrahedron
  std::vector<std::array<std::size_t, 6>> vertex_faces;     // the faces at v in slot order 01, 02, 03, 12, 13, 23
  std::vector<std::array<std::size_t, 4>> tetrahedron_vertices;  // the triangulation's vertices of each tetrahedron, 0..V−1

  std::size_t faces() const { return face_edge.size(); }
  std::size_t vertices() const { return vertex_tetrahedron.size(); }
  // Spins on the faces as colours of the triangulation's edges, and back.
  std::vector<int> edge_colours(const std::vector<int>& two_j) const {
    if (two_j.size() != faces()) throw std::invalid_argument("physis::quantum_gravity::SpinFoam: one spin per face");
    std::vector<int> c(triangulation.edges);
    for (std::size_t f = 0; f < faces(); f++) c[face_edge[f]] = two_j[f];
    return c;
  }
  // The six spins of vertex v in slot order.
  std::array<int, 6> vertex_spins(const std::vector<int>& two_j, std::size_t v) const {
    std::array<int, 6> s{};
    for (int k = 0; k < 6; k++) s[static_cast<std::size_t>(k)] = two_j.at(vertex_faces[v][static_cast<std::size_t>(k)]);
    return s;
  }
};

namespace detail {

inline void finish(SpinFoam& F) {
  const auto& T = F.triangulation;
  F.edge_face.assign(T.edges, SIZE_MAX);
  for (std::size_t f = 0; f < F.face_edge.size(); f++) {
    if (F.edge_face[F.face_edge[f]] != SIZE_MAX) throw std::invalid_argument("physis::quantum_gravity::spin_foam: two faces around one edge");
    F.edge_face[F.face_edge[f]] = f;
  }
  for (std::size_t e = 0; e < T.edges; e++)
    if (F.edge_face[e] == SIZE_MAX) throw std::invalid_argument("physis::quantum_gravity::spin_foam: an edge of the triangulation has no face (the manifold has a boundary)");
  F.vertex_faces.clear();
  for (std::size_t v = 0; v < F.vertex_tetrahedron.size(); v++) {
    std::array<std::size_t, 6> fs{};
    for (int k = 0; k < 6; k++) fs[static_cast<std::size_t>(k)] = F.edge_face[T.tetrahedra[F.vertex_tetrahedron[v]][static_cast<std::size_t>(k)]];
    F.vertex_faces.push_back(fs);
  }
}

}  // namespace detail

// From a closed triangulated 3-manifold: morphe's dual 2-complex (dual vertices in the order
// of the tetrahedra, dual faces one per edge).
inline SpinFoam spin_foam(const morphe::topology::SimplicialComplex& K) {
  SpinFoam F;
  F.triangulation = triangulation(K);
  const auto D = morphe::topology::dual_two_complex(K);
  F.complex = D.complex;
  F.complex.name = "spin foam on " + K.name;
  for (const auto& h : D.face_simplex) F.face_edge.push_back(K.find(h));
  for (std::size_t v = 0; v < D.complex.vertices; v++) F.vertex_tetrahedron.push_back(v);
  // The triangulation's vertices, renumbered 0..V−1 in the order of K's 0-simplices.
  std::vector<std::size_t> number(K.vertices, SIZE_MAX);
  for (std::size_t i = 0; i < K.count(0); i++) number[K.simplices[0][i][0]] = i;
  for (const auto& t : K.simplices[3]) F.tetrahedron_vertices.push_back({number[t[0]], number[t[1]], number[t[2]], number[t[3]]});
  detail::finish(F);
  return F;
}

// From a 4-coloured graph: its 2-complex of bicoloured cycles, whose faces are listed in the
// order of the triangulation's edges and whose vertices are the tetrahedra.
inline SpinFoam spin_foam(const morphe::topology::ColouredGraph& G) {
  if (G.colours != 4) throw std::invalid_argument("physis::quantum_gravity::spin_foam: a 4-coloured graph");
  SpinFoam F;
  F.triangulation = triangulation(G);
  F.complex = G.two_complex();
  F.complex.name = "spin foam of a 4-coloured graph";
  for (std::size_t f = 0; f < F.complex.faces.size(); f++) F.face_edge.push_back(f);
  for (std::size_t v = 0; v < 2 * G.p; v++) F.vertex_tetrahedron.push_back(v);
  // The vertex of colour i of tetrahedron v is the 3-bubble without colour i through v.
  const std::uint32_t all = 15u;
  std::array<std::vector<std::size_t>, 4> comp;
  std::array<std::size_t, 4> offset{};
  std::size_t total = 0;
  for (std::size_t i = 0; i < 4; i++) {
    std::size_t n = 0;
    comp[i] = G.components(all & ~(1u << i), &n);
    offset[i] = total;
    total += n;
  }
  for (std::size_t v = 0; v < 2 * G.p; v++) F.tetrahedron_vertices.push_back({offset[0] + comp[0][v], offset[1] + comp[1][v], offset[2] + comp[2][v], offset[3] + comp[3][v]});
  detail::finish(F);
  return F;
}

// ---- weights ----------------------------------------------------------------------------------------------------------

// The face weight (−1)^{2j}[2j + 1]_q, or (−1)^{2j}(2j + 1) at r = 0.
inline double face_weight(int two_j, int r) { return r == 0 ? (two_j % 2 ? -1.0 : 1.0) * (two_j + 1) : symmetry::quantum_dimension_signed(two_j, r); }
// Whether a triad is admissible: the triangle conditions, and j₁ + j₂ + j₃ ≤ r − 2 when r > 0.
inline bool triad(int a, int b, int c, int r) { return r == 0 ? symmetry::triangle(a, b, c) : symmetry::admissible(a, b, c, r); }
// The vertex amplitude {6j} of six spins in slot order 01, 02, 03, 12, 13, 23: the 6j symbol
// {01 02 12; 23 13 03}, whose triads (01 02 12), (01 13 03), (23 02 03), (23 13 12) are the
// four triangles.
inline double vertex_amplitude(const std::array<int, 6>& s, int r) {
  if (r == 0) return symmetry::wigner_6j(s[0], s[1], s[3], s[5], s[4], s[2]).value();
  return symmetry::quantum_6j(r, s[0], s[1], s[3], s[5], s[4], s[2]);
}

struct Amplitude {
  double value = 0;          // A(j), normalised by w^{−2V} when r > 0
  bool admissible = false;   // every triangle admissible
  double faces = 1, edges = 1, vertices = 1, normalisation = 1;  // the factors
  std::vector<double> vertex;  // the amplitude of each spin-foam vertex
};

// The amplitude of a spin assignment (two_j per face) at level r (r = 0: Ponzano–Regge).
inline Amplitude amplitude(const SpinFoam& F, const std::vector<int>& two_j, int r) {
  if (r != 0 && r < 3) throw std::invalid_argument("physis::quantum_gravity::amplitude: r = 0 or r >= 3");
  const auto c = F.edge_colours(two_j);
  Amplitude A;
  for (int x : two_j)
    if (x < 0 || (r > 0 && x > r - 2)) return A;
  for (const auto& tr : F.triangulation.triangles) {
    const int a = c[tr[0]], b = c[tr[1]], d = c[tr[2]];
    if (!triad(a, b, d, r)) return A;
    if (((a + b + d) / 2) % 2) A.edges = -A.edges;
  }
  A.admissible = true;
  for (int x : two_j) A.faces *= face_weight(x, r);
  for (std::size_t v = 0; v < F.vertices(); v++) {
    A.vertex.push_back(vertex_amplitude(F.vertex_spins(two_j, v), r));
    A.vertices *= A.vertex.back();
  }
  if (r > 0) A.normalisation = std::pow(symmetry::total_quantum_dimension_squared(r), -static_cast<double>(F.triangulation.vertices));
  A.value = A.normalisation * A.faces * A.edges * A.vertices;
  return A;
}

// Calls f(two_j, A) for every admissible assignment at level r > 0 (backtracking over the
// faces in an order that closes triangles early); returns their number.
inline std::uint64_t for_each_state(const SpinFoam& F, int r, const std::function<void(const std::vector<int>&, double)>& f) {
  if (r < 3) throw std::invalid_argument("physis::quantum_gravity::for_each_state: r >= 3");
  const auto& T = F.triangulation;
  const std::size_t E = T.edges;
  // Edges in the order of the tetrahedra, as turaev_viro() orders them.
  std::vector<std::size_t> order, position(E, SIZE_MAX);
  for (const auto& t : T.tetrahedra)
    for (auto e : t)
      if (position[e] == SIZE_MAX) position[e] = order.size(), order.push_back(e);
  for (std::size_t e = 0; e < E; e++)
    if (position[e] == SIZE_MAX) position[e] = order.size(), order.push_back(e);
  std::vector<std::vector<std::size_t>> tri_done(E);
  for (std::size_t t = 0; t < T.triangles.size(); t++) {
    std::size_t last = 0;
    for (auto e : T.triangles[t]) last = std::max(last, position[e]);
    tri_done[last].push_back(t);
  }
  std::vector<int> colour(E, 0), two_j(F.faces(), 0);
  std::uint64_t count = 0;
  std::function<void(std::size_t)> go = [&](std::size_t k) {
    if (k == E) {
      for (std::size_t g = 0; g < F.faces(); g++) two_j[g] = colour[F.face_edge[g]];
      count++;
      f(two_j, amplitude(F, two_j, r).value);
      return;
    }
    for (int j = 0; j <= r - 2; j++) {
      colour[order[k]] = j;
      bool ok = true;
      for (auto t : tri_done[k]) {
        const auto& tr = T.triangles[t];
        if (!symmetry::admissible(colour[tr[0]], colour[tr[1]], colour[tr[2]], r)) { ok = false; break; }
      }
      if (ok) go(k + 1);
    }
  };
  go(0);
  return count;
}

// ---- quantum tetrahedra -------------------------------------------------------------------------------------------------

// A spin-foam vertex as a quantum tetrahedron, and its boundary spin network.
struct QuantumTetrahedron {
  std::array<int, 6> two_j{};                         // spins in slot order 01, 02, 03, 12, 13, 23
  std::array<std::array<int, 2>, 6> links{};          // the boundary spin network: link k joins the two triangles (nodes 0..3, the triangle opposite vertex i) that share edge k
  std::array<double, 6> lengths{};                    // ℓ = j + ½
  bool euclidean = false;                             // whether these lengths bound a Euclidean tetrahedron
  double volume = 0;                                  // when Euclidean
  std::array<double, 6> dihedral{};                   // interior dihedral angles at the edges, when Euclidean
  double six_j = 0;                                   // the vertex amplitude at level r (r = 0: classical)
  double ponzano_regge = std::numeric_limits<double>::quiet_NaN();  // cos(Σ(j+½)θ_ext + π/4)/√(12πV), when Euclidean
  double envelope = std::numeric_limits<double>::quiet_NaN();       // 1/√(12πV)
};

// The geometry alone (lengths, links, volume, dihedral angles, the Ponzano–Regge value),
// without the 6j symbol.
inline QuantumTetrahedron tetrahedron_geometry(const std::array<int, 6>& two_j) {
  QuantumTetrahedron Q;
  Q.two_j = two_j;
  for (int k = 0; k < 6; k++) {
    const auto [a, b] = tetrahedron_slots[static_cast<std::size_t>(k)];
    int other[2], n = 0;
    for (int i = 0; i < 4; i++)
      if (i != a && i != b) other[n++] = i;
    Q.links[static_cast<std::size_t>(k)] = {other[0], other[1]};
    Q.lengths[static_cast<std::size_t>(k)] = two_j[static_cast<std::size_t>(k)] / 2.0 + 0.5;
  }
  spacetime::SquaredLengths s(4, std::vector<double>(4, 0.0));
  for (int k = 0; k < 6; k++) {
    const auto [a, b] = tetrahedron_slots[static_cast<std::size_t>(k)];
    s[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] = s[static_cast<std::size_t>(b)][static_cast<std::size_t>(a)] = Q.lengths[static_cast<std::size_t>(k)] * Q.lengths[static_cast<std::size_t>(k)];
  }
  // Euclidean exactly when the Gram matrix is positive definite (its leading minors positive).
  const auto G = spacetime::gram(s);
  const double m1 = G[0][0], m2 = G[0][0] * G[1][1] - G[0][1] * G[1][0];
  Q.euclidean = m1 > 0 && m2 > 0 && spacetime::detail::det(G) > 0;
  if (!Q.euclidean) return Q;
  Q.volume = spacetime::simplex_volume(s);
  double phase = std::numbers::pi / 4;
  for (int k = 0; k < 6; k++) {
    const auto& o = Q.links[static_cast<std::size_t>(k)];
    Q.dihedral[static_cast<std::size_t>(k)] = spacetime::dihedral_angle(s, static_cast<std::size_t>(o[0]), static_cast<std::size_t>(o[1]));
    phase += Q.lengths[static_cast<std::size_t>(k)] * (std::numbers::pi - Q.dihedral[static_cast<std::size_t>(k)]);
  }
  Q.envelope = 1 / std::sqrt(12 * std::numbers::pi * Q.volume);
  Q.ponzano_regge = std::cos(phase) * Q.envelope;
  return Q;
}

// The quantum tetrahedron with its vertex amplitude at level r (r = 0: the classical 6j).
inline QuantumTetrahedron quantum_tetrahedron(const std::array<int, 6>& two_j, int r) {
  QuantumTetrahedron Q = tetrahedron_geometry(two_j);
  Q.six_j = vertex_amplitude(two_j, r);
  return Q;
}

inline QuantumTetrahedron vertex_geometry(const SpinFoam& F, const std::vector<int>& two_j, std::size_t v, int r) {
  return quantum_tetrahedron(F.vertex_spins(two_j, v), r);
}

// ---- Regge geometry of a spin assignment --------------------------------------------------------------------------------

struct ReggeGeometry {
  std::size_t euclidean = 0;          // the number of tetrahedra that are Euclidean
  bool complete = false;              // all of them
  std::vector<double> deficit;        // per edge of the triangulation, 2π − Σ θ, when complete
  std::vector<double> length;         // per edge, j + ½
  double action = std::numeric_limits<double>::quiet_NaN();  // Σ ℓ_e ε_e, when complete
};

inline ReggeGeometry regge_geometry(const SpinFoam& F, const std::vector<int>& two_j) {
  const auto& T = F.triangulation;
  const auto c = F.edge_colours(two_j);
  ReggeGeometry R;
  R.length.resize(T.edges);
  for (std::size_t e = 0; e < T.edges; e++) R.length[e] = c[e] / 2.0 + 0.5;
  std::vector<double> sum(T.edges, 0.0);
  for (std::size_t v = 0; v < F.vertices(); v++) {
    const auto Q = tetrahedron_geometry(F.vertex_spins(two_j, v));
    if (!Q.euclidean) continue;
    R.euclidean++;
    const auto& tet = T.tetrahedra[F.vertex_tetrahedron[v]];
    for (int k = 0; k < 6; k++) sum[tet[static_cast<std::size_t>(k)]] += Q.dihedral[static_cast<std::size_t>(k)];
  }
  R.complete = R.euclidean == F.vertices();
  if (!R.complete) return R;
  R.deficit.resize(T.edges);
  R.action = 0;
  for (std::size_t e = 0; e < T.edges; e++) {
    R.deficit[e] = 2 * std::numbers::pi - sum[e];
    R.action += R.length[e] * R.deficit[e];
  }
  return R;
}

// ---- sampling -----------------------------------------------------------------------------------------------------------

struct SampleOptions {
  int r = 4;                      // level, r >= 3
  std::size_t thermalise = 200;   // sweeps before measuring
  std::size_t sweeps = 2000;      // sweeps measured
  std::uint64_t seed = 1;
};

struct Sample {
  std::vector<int> last;                          // the final assignment (two_j per face)
  std::vector<double> mean_spin;                  // ⟨j_f⟩ under |A|, per face
  std::vector<std::vector<double>> histogram;     // per face, the frequency of each two_j in 0..r−2
  double mean_sign = 0;                           // ⟨sign A⟩ under |A|: Σ A / Σ |A|
  double positive = 0;                            // the fraction of measured states with A > 0
  double acceptance = 0;
  std::vector<double> history;                    // log |A| after each measured sweep
};

// Metropolis with weight |A| at level r: each sweep proposes a new spin on every face in turn
// and a change by ±1 (independently on each edge) of the spins on all edges at each vertex of
// the triangulation; loops at a vertex are left unchanged so that parity is kept.
inline Sample sample(const SpinFoam& F, const SampleOptions& o) {
  if (o.r < 3) throw std::invalid_argument("physis::quantum_gravity::sample: r >= 3");
  const int r = o.r, top = r - 2;
  const auto& T = F.triangulation;
  const std::size_t nf = F.faces();
  // Faces of each triangle, triangles and spin-foam vertices of each face.
  std::vector<std::vector<std::size_t>> tris_of(nf), verts_of(nf);
  for (std::size_t t = 0; t < T.triangles.size(); t++)
    for (auto e : T.triangles[t]) {
      auto& v = tris_of[F.edge_face[e]];
      if (std::find(v.begin(), v.end(), t) == v.end()) v.push_back(t);
    }
  for (std::size_t v = 0; v < F.vertices(); v++)
    for (auto f : F.vertex_faces[v]) {
      auto& w = verts_of[f];
      if (std::find(w.begin(), w.end(), v) == w.end()) w.push_back(v);
    }
  // The faces at each vertex of the triangulation (edges with one end there).
  std::vector<std::vector<std::size_t>> star(T.vertices);
  {
    std::vector<std::array<std::size_t, 2>> ends(T.edges, {SIZE_MAX, SIZE_MAX});
    for (std::size_t v = 0; v < F.vertices(); v++) {
      const auto& tet = T.tetrahedra[F.vertex_tetrahedron[v]];
      for (int k = 0; k < 6; k++) {
        const auto [a, b] = tetrahedron_slots[static_cast<std::size_t>(k)];
        ends[tet[static_cast<std::size_t>(k)]] = {F.tetrahedron_vertices[v][static_cast<std::size_t>(a)], F.tetrahedron_vertices[v][static_cast<std::size_t>(b)]};
      }
    }
    for (std::size_t e = 0; e < T.edges; e++) {
      const auto [a, b] = ends[e];
      if (a == SIZE_MAX || a == b) continue;  // loops keep their spin under vertex moves
      star[a].push_back(F.edge_face[e]);
      star[b].push_back(F.edge_face[e]);
    }
  }
  std::vector<int> j(nf, 0);
  // The local factors that change with the faces in `changed`: |face weights|, triangle
  // admissibility and signs, vertex amplitudes; returns the log of |weight| and its sign
  // (−inf when inadmissible).
  auto local = [&](const std::vector<std::size_t>& changed, double& log_abs, int& sign) {
    std::vector<std::size_t> tris, verts;
    for (auto f : changed) {
      for (auto t : tris_of[f]) tris.push_back(t);
      for (auto v : verts_of[f]) verts.push_back(v);
    }
    std::sort(tris.begin(), tris.end());
    tris.erase(std::unique(tris.begin(), tris.end()), tris.end());
    std::sort(verts.begin(), verts.end());
    verts.erase(std::unique(verts.begin(), verts.end()), verts.end());
    log_abs = 0;
    sign = 1;
    for (auto f : changed) {
      if (j[f] < 0 || j[f] > top) { log_abs = -INFINITY; return; }
      const double w = face_weight(j[f], r);
      log_abs += std::log(std::fabs(w));
      if (w < 0) sign = -sign;
    }
    for (auto t : tris) {
      const auto& tr = T.triangles[t];
      const int a = j[F.edge_face[tr[0]]], b = j[F.edge_face[tr[1]]], c = j[F.edge_face[tr[2]]];
      if (!symmetry::admissible(a, b, c, r)) { log_abs = -INFINITY; return; }
      if (((a + b + c) / 2) % 2) sign = -sign;
    }
    for (auto v : verts) {
      const double x = vertex_amplitude(F.vertex_spins(j, v), r);
      if (x == 0) { log_abs = -INFINITY; return; }
      log_abs += std::log(std::fabs(x));
      if (x < 0) sign = -sign;
    }
  };
  calx::random::Stream z(o.seed, 0x666f616du, 0, 0);  // "foam"
  Sample S;
  S.histogram.assign(nf, std::vector<double>(static_cast<std::size_t>(top + 1), 0.0));
  S.mean_spin.assign(nf, 0.0);
  std::size_t proposed = 0, accepted = 0, measured = 0;
  double sum_sign = 0;
  const auto full = amplitude(F, j, r);
  double log_total = std::log(std::fabs(full.value));
  int total_sign = full.value < 0 ? -1 : 1;
  auto propose = [&](const std::vector<std::size_t>& faces, const std::vector<int>& next) {
    double lo, ln;
    int so, sn;
    local(faces, lo, so);
    std::vector<int> old(faces.size());
    for (std::size_t k = 0; k < faces.size(); k++) old[k] = j[faces[k]], j[faces[k]] = next[k];
    local(faces, ln, sn);
    proposed++;
    if (std::isfinite(ln) && (ln >= lo || z.uniform() < std::exp(ln - lo))) {
      accepted++;
      log_total += ln - lo;
      total_sign *= so * sn;
      return;
    }
    for (std::size_t k = 0; k < faces.size(); k++) j[faces[k]] = old[k];
  };
  const std::size_t total_sweeps = o.thermalise + o.sweeps;
  for (std::size_t sweep = 0; sweep < total_sweeps; sweep++) {
    for (std::size_t f = 0; f < nf; f++) {
      if (top == 0) break;
      int x = static_cast<int>(z.below(static_cast<std::uint64_t>(top)));  // a new value, different from the present one
      if (x >= j[f]) x++;
      propose({f}, {x});
    }
    for (std::size_t v = 0; v < T.vertices; v++) {
      if (star[v].empty()) continue;
      std::vector<std::size_t> faces;
      std::vector<int> next;
      for (auto f : star[v])
        if (std::find(faces.begin(), faces.end(), f) == faces.end()) faces.push_back(f);
      for (auto f : faces) next.push_back(j[f] + (z.uniform() < 0.5 ? 1 : -1));
      propose(faces, next);
    }
    if (sweep < o.thermalise) continue;
    measured++;
    for (std::size_t f = 0; f < nf; f++) {
      S.histogram[f][static_cast<std::size_t>(j[f])] += 1;
      S.mean_spin[f] += j[f] / 2.0;
    }
    sum_sign += total_sign;
    if (total_sign > 0) S.positive += 1;
    S.history.push_back(log_total);
  }
  for (std::size_t f = 0; f < nf; f++) {
    S.mean_spin[f] /= static_cast<double>(measured);
    for (auto& h : S.histogram[f]) h /= static_cast<double>(measured);
  }
  S.mean_sign = sum_sign / static_cast<double>(measured);
  S.positive /= static_cast<double>(measured);
  S.acceptance = proposed ? static_cast<double>(accepted) / static_cast<double>(proposed) : 0;
  S.last = j;
  return S;
}

}  // namespace physis::quantum_gravity
