// Cell complexes realised in space: their coboundaries in every degree, their measures, their facets,
// and boxes (decision 0012; docs/theory/fields.md §9).
//
// The topology of a cellulation is morphe's CW complex; its geometry is additional data, a position for
// each vertex in a frame. The two are kept apart (the topology fixes d, the geometry the Hodge stars),
// and checked against each other.
//
// Coboundaries. d_k: Cᵏ → Cᵏ⁺¹ is the transpose of the cellular boundary ∂_{k+1} (morphe's chain
// complex), held as a sparse matrix with the integer incidences as values; d_{k+1} d_k = 0 because
// ∂_k ∂_{k+1} = 0, which morphe validates (Desbrun, Hirani, Leok and Marsden 2005, §3).
//
// Measures. An edge's length; a 2-cell's vector area ½ Σᵢ pᵢ × pᵢ₊₁ round its attaching word, whose
// magnitude is its area and whose direction (the word's orientation by the right-hand rule) its unit
// normal, its centroid that of the fan of triangles from its first vertex; a 3-cell's signed volume
// ⅓ Σ_f [c : f] (c_f · A_f) by the divergence theorem with the field x/3, over its faces with their
// incidences, and its centroid from the pyramids on its faces with apex at the mean of their centroids
// (a pyramid's centroid lies a quarter of the way from its base to its apex). These are exact for cells
// with planar faces, as every cell of a box is.
//
// Facets. In a complex of dimension n embedded in space of dimension n (a box of three dimensions, a
// planar complex, an interval on a line), each (n − 1)-cell has a conormal: a unit vector across it (a
// 2-cell's normal; an edge's direction turned a quarter in the plane; a vertex's, along the line). It
// points out of one top cell and into the other, or out of the one a boundary facet bounds; an
// interior facet thus has a cell behind it and one ahead. The orientation o of a top cell (the sign of
// its signed volume) and the incidences of d agree with this geometry: d(i, f) = o_i σ(i, f), σ = +1
// when the conormal points out of i, which measure() lets a caller check.
//
// Boxes. A box is the product of intervals cut at given coordinates (morphe product), each cell knowing
// the axes it spans and the indices of its least vertex, so that a point is located by a binary search
// along each axis. A box refined by r cuts each interval into r parts; its cells are nested in the
// coarse box's.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <calx/linalg/sparse.hpp>
#include <morphe/topology/constructions.hpp>
#include <physis/discretisation/cochains.hpp>
#include <physis/geometry/rotation.hpp>

namespace physis::discretisation {

using geometry::Vec3;
inline constexpr std::size_t no_cell = std::numeric_limits<std::size_t>::max();

// ---- coboundaries ---------------------------------------------------------------------------------------------

struct Coboundary {
  std::size_t k = 0;  // d_k: Cᵏ → Cᵏ⁺¹
  calx::linalg::Csr<std::uint32_t> pattern;
  std::vector<double> value;

  Cochain apply(const Cochain& c) const {
    if (c.size() != pattern.cols) throw std::invalid_argument("physis::discretisation::Coboundary::apply: one value per k-cell");
    Cochain out(pattern.rows, 0.0);
    for (std::size_t r = 0; r < pattern.rows; r++)
      for (std::uint64_t e = pattern.begin(r); e < pattern.end(r); e++) out[r] += value[e] * c[pattern.col[e]];
    return out;
  }
  Cochain apply_transpose(const Cochain& c) const {
    if (c.size() != pattern.rows) throw std::invalid_argument("physis::discretisation::Coboundary::apply_transpose: one value per (k + 1)-cell");
    Cochain out(pattern.cols, 0.0);
    for (std::size_t r = 0; r < pattern.rows; r++)
      for (std::uint64_t e = pattern.begin(r); e < pattern.end(r); e++) out[pattern.col[e]] += value[e] * c[r];
    return out;
  }
};

inline Coboundary coboundary(const CWComplex& X, std::size_t k) {
  if (k + 1 > X.dimension()) throw std::invalid_argument("physis::discretisation::coboundary: d_k needs cells of dimension k + 1");
  const auto C = X.chain_complex();
  const auto& D = C.boundary[k + 1];  // ∂_{k+1}: columns the (k + 1)-cells
  std::vector<std::uint32_t> r, c;
  std::vector<double> v;
  for (std::size_t i = 0; i < D.cols(); i++)
    for (const auto& e : D.columns[i]) {
      r.push_back(static_cast<std::uint32_t>(i));
      c.push_back(static_cast<std::uint32_t>(e.index));
      v.push_back(static_cast<double>(e.value.to_int64()));
    }
  Coboundary d;
  d.k = k;
  std::vector<std::uint64_t> order;
  d.pattern = calx::linalg::compress<std::uint32_t>(D.cols(), C.n[k], r, c, &order);
  d.value.resize(v.size());
  for (std::size_t e = 0; e < v.size(); e++) d.value[e] = v[order[e]];
  return d;
}

// ---- geometry -------------------------------------------------------------------------------------------------

struct GeometricComplex {
  CWComplex topology;
  std::vector<Vec3> position;  // of each vertex, m, on the frame's axes
  std::string frame = "box";

  std::size_t dimension() const { return topology.dimension(); }
};

struct Measures {
  std::size_t dimension = 0;
  std::vector<double> length;                 // of each edge
  std::vector<double> area;                   // of each 2-cell
  std::vector<Vec3> normal;                   // of each 2-cell: unit, by its word's orientation
  std::vector<double> volume;                 // of each 3-cell
  std::vector<int> orientation;               // of each top cell: the sign of its signed measure
  std::array<std::vector<Vec3>, 4> centroid;  // of the cells of each dimension

  // The k-cells' measures: 1 for vertices, lengths, areas, volumes.
  std::vector<double> of(std::size_t k) const {
    switch (k) {
      case 0: return std::vector<double>(centroid[0].size(), 1.0);
      case 1: return length;
      case 2: return area;
      case 3: return volume;
    }
    throw std::invalid_argument("physis::discretisation::Measures: cells of dimension at most 3");
  }
};

namespace geometric_detail {

// The vertices of a 2-cell's word, in order, from its base.
inline std::vector<std::size_t> word_vertices(const CWComplex& X, const morphe::topology::Face& f) {
  std::vector<std::size_t> v{f.base};
  std::size_t at = f.base;
  for (const auto& s : f.boundary) {
    at = s.sign > 0 ? X.edges[s.edge][1] : X.edges[s.edge][0];
    v.push_back(at);
  }
  v.pop_back();  // the word is closed: its last vertex is its base
  return v;
}

}  // namespace geometric_detail

inline Measures measure(const GeometricComplex& K) {
  using namespace geometry::operators;
  const CWComplex& X = K.topology;
  if (K.position.size() != X.vertices) throw std::invalid_argument("physis::discretisation::measure: a position for each vertex");
  Measures M;
  M.dimension = X.dimension();
  if (M.dimension > 3) throw std::invalid_argument("physis::discretisation::measure: complexes of dimension at most 3");
  M.centroid[0] = K.position;
  for (const auto& e : X.edges) {
    const Vec3 a = K.position[e[0]], b = K.position[e[1]];
    M.length.push_back(geometry::norm(b - a));
    M.centroid[1].push_back(0.5 * (a + b));
  }
  std::vector<Vec3> vector_area;
  for (const auto& f : X.faces) {
    const auto v = geometric_detail::word_vertices(X, f);
    Vec3 A{}, c{};
    double total = 0;
    for (std::size_t i = 0; i < v.size(); i++) A = A + 0.5 * geometry::cross(K.position[v[i]], K.position[v[(i + 1) % v.size()]]);
    const double area = geometry::norm(A);
    const Vec3 n = area > 0 ? (1 / area) * A : Vec3{};
    for (std::size_t i = 1; i + 1 < v.size(); i++) {
      const Vec3 p0 = K.position[v[0]], p1 = K.position[v[i]], p2 = K.position[v[i + 1]];
      const double t = 0.5 * geometry::dot(geometry::cross(p1 - p0, p2 - p0), n);
      c = c + t * ((1.0 / 3) * (p0 + p1 + p2));
      total += t;
    }
    if (!(total != 0)) {  // a degenerate cell: its vertices' mean
      c = {};
      for (std::size_t i : v) c = c + K.position[i];
      c = (1.0 / static_cast<double>(v.size())) * c;
    } else {
      c = (1 / total) * c;
    }
    vector_area.push_back(A);
    M.area.push_back(area);
    M.normal.push_back(n);
    M.centroid[2].push_back(c);
  }
  if (M.dimension == 3) {
    for (const auto& cell : X.higher[0]) {
      Vec3 o{};
      for (const auto& e : cell) o = o + M.centroid[2][e.index];
      o = (1.0 / static_cast<double>(cell.size())) * o;
      double signed_volume = 0;
      Vec3 c{};
      for (const auto& e : cell) {
        const double inc = static_cast<double>(e.value.to_int64());
        const double v = inc * geometry::dot(M.centroid[2][e.index] - o, vector_area[e.index]) / 3;
        signed_volume += v;
        c = c + v * (o + 0.75 * (M.centroid[2][e.index] - o));
      }
      M.volume.push_back(std::fabs(signed_volume));
      M.orientation.push_back(signed_volume >= 0 ? 1 : -1);
      M.centroid[3].push_back(signed_volume != 0 ? (1 / signed_volume) * c : o);
    }
  }
  // The orientations of top cells below dimension 3: against the plane's normal (2) or the line (1).
  if (M.dimension == 2) {
    const Vec3 plane = M.normal.empty() ? Vec3{0, 0, 1} : M.normal[0];
    for (const Vec3& n : M.normal) M.orientation.push_back(geometry::dot(n, plane) >= 0 ? 1 : -1);
  } else if (M.dimension == 1) {
    M.orientation.assign(X.edges.size(), 1);
  } else if (M.dimension == 0) {
    M.orientation.assign(X.vertices, 1);
  }
  return M;
}

// The facets of a complex of dimension n in space of dimension n, with their conormals and their cells.
struct Facets {
  std::size_t n = 0;                              // the top dimension; the facets have dimension n − 1
  std::vector<double> measure;                    // area, length, or 1
  std::vector<Vec3> conormal;                     // unit, across the facet
  std::vector<Vec3> centroid;
  std::vector<std::array<std::size_t, 2>> cells;  // [behind, ahead]: the conormal points out of the first and into the second; no_cell where none

  bool boundary(std::size_t f) const { return cells[f][0] == no_cell || cells[f][1] == no_cell; }
  // The one cell of a boundary facet, and +1 if the conormal points out of it.
  std::pair<std::size_t, int> bounded(std::size_t f) const { return cells[f][0] != no_cell ? std::pair{cells[f][0], 1} : std::pair{cells[f][1], -1}; }
};

inline Facets facets(const GeometricComplex& K, const Measures& M) {
  using namespace geometry::operators;
  const CWComplex& X = K.topology;
  const std::size_t n = M.dimension;
  if (n == 0) throw std::invalid_argument("physis::discretisation::facets: a complex of dimension at least 1");
  Facets F;
  F.n = n;
  const std::size_t count = X.cells(n - 1);
  F.cells.assign(count, {no_cell, no_cell});
  if (n == 3) {
    F.measure = M.area;
    F.conormal = M.normal;
  } else if (n == 2) {
    const Vec3 plane = M.normal[0];
    for (std::size_t e = 0; e < X.edges.size(); e++) {
      const Vec3 t = K.position[X.edges[e][1]] - K.position[X.edges[e][0]];
      F.measure.push_back(M.length[e]);
      F.conormal.push_back(geometry::unit(geometry::cross(t, plane)));
    }
  } else {
    const Vec3 line = geometry::unit(K.position[X.edges[0][1]] - K.position[X.edges[0][0]]);
    F.measure.assign(count, 1.0);
    F.conormal.assign(count, line);
  }
  F.centroid = M.centroid[n - 1];
  const Coboundary d = coboundary(X, n - 1);
  for (std::size_t i = 0; i < d.pattern.rows; i++)
    for (std::uint64_t e = d.pattern.begin(i); e < d.pattern.end(i); e++) {
      const std::size_t f = d.pattern.col[e];
      const bool out = geometry::dot(F.centroid[f] - M.centroid[n][i], F.conormal[f]) > 0;
      std::size_t& slot = F.cells[f][out ? 0 : 1];
      if (slot != no_cell) throw std::invalid_argument("physis::discretisation::facets: a facet with two cells on one side");
      slot = i;
    }
  return F;
}

// ---- boxes ----------------------------------------------------------------------------------------------------

struct BoxComplex {
  struct Cell {
    std::uint8_t axes = 0;         // bit a set when the cell spans axis a
    std::array<int, 3> corner{};   // the indices of its least vertex along each axis
  };

  GeometricComplex complex;
  std::array<std::vector<double>, 3> coordinates;  // each axis's cuts, increasing
  std::array<int, 3> n{};                          // cells along each axis
  std::array<std::vector<Cell>, 4> cells;          // the cells of each dimension, as morphe numbers them

  std::size_t dimension() const { return complex.dimension(); }

  std::size_t find(std::size_t k, std::uint8_t axes, const std::array<int, 3>& corner) const {
    const auto& index = lookup_[k];
    const std::uint64_t key = key_of(axes, corner);
    const auto it = std::lower_bound(index.begin(), index.end(), std::pair<std::uint64_t, std::size_t>{key, 0});
    return it != index.end() && it->first == key ? it->second : no_cell;
  }

  // The top cell holding a point, and the point's coordinates in it in [0, 1] along each axis it spans.
  struct Location {
    std::size_t cell = no_cell;
    std::array<int, 3> corner{};
    std::array<double, 3> local{};
  };
  std::optional<Location> locate(const Vec3& x) const {
    Location L;
    for (int a = 0; a < 3; a++) {
      const auto& c = coordinates[static_cast<std::size_t>(a)];
      if (n[static_cast<std::size_t>(a)] == 0) continue;
      if (!(x[static_cast<std::size_t>(a)] >= c.front() && x[static_cast<std::size_t>(a)] <= c.back())) return std::nullopt;
      const auto it = std::upper_bound(c.begin(), c.end(), x[static_cast<std::size_t>(a)]);
      const int i = std::clamp(static_cast<int>(it - c.begin()) - 1, 0, n[static_cast<std::size_t>(a)] - 1);
      L.corner[static_cast<std::size_t>(a)] = i;
      L.local[static_cast<std::size_t>(a)] = (x[static_cast<std::size_t>(a)] - c[static_cast<std::size_t>(i)]) / (c[static_cast<std::size_t>(i) + 1] - c[static_cast<std::size_t>(i)]);
    }
    L.cell = find(dimension(), full_axes(), L.corner);
    return L;
  }

  std::uint8_t full_axes() const {
    std::uint8_t m = 0;
    for (int a = 0; a < 3; a++)
      if (n[static_cast<std::size_t>(a)] > 0) m = static_cast<std::uint8_t>(m | (1u << a));
    return m;
  }

  // Indexes the cells by (axes, corner); called by box_complex.
  void index() {
    for (std::size_t k = 0; k < 4; k++) {
      lookup_[k].clear();
      for (std::size_t i = 0; i < cells[k].size(); i++) lookup_[k].emplace_back(key_of(cells[k][i].axes, cells[k][i].corner), i);
      std::sort(lookup_[k].begin(), lookup_[k].end());
    }
  }

 private:
  static std::uint64_t key_of(std::uint8_t axes, const std::array<int, 3>& c) {
    return (static_cast<std::uint64_t>(axes) << 60) | (static_cast<std::uint64_t>(c[0]) << 40) | (static_cast<std::uint64_t>(c[1]) << 20) | static_cast<std::uint64_t>(c[2]);
  }
  std::array<std::vector<std::pair<std::uint64_t, std::size_t>>, 4> lookup_;
};

// The box with the given cuts along each axis (an axis with one cut is a point of the box: no cells
// along it), in a frame.
inline BoxComplex box_complex(std::array<std::vector<double>, 3> coordinates, std::string frame = "box") {
  BoxComplex B;
  for (int a = 0; a < 3; a++) {
    const auto& c = coordinates[static_cast<std::size_t>(a)];
    if (c.empty()) throw std::invalid_argument("physis::discretisation::box_complex: at least one cut along each axis");
    for (std::size_t i = 1; i < c.size(); i++)
      if (!(c[i] > c[i - 1])) throw std::invalid_argument("physis::discretisation::box_complex: cuts increase along each axis");
    if (c.size() > (1u << 20)) throw std::invalid_argument("physis::discretisation::box_complex: at most 2²⁰ cells along an axis");
    B.n[static_cast<std::size_t>(a)] = static_cast<int>(c.size()) - 1;
  }
  B.coordinates = std::move(coordinates);
  const auto I = [&](int a) { return morphe::topology::interval(static_cast<std::size_t>(B.n[static_cast<std::size_t>(a)])); };
  CWComplex X = morphe::topology::product(morphe::topology::product(I(0), I(1)), I(2));
  X.name = "box " + std::to_string(B.n[0]) + "×" + std::to_string(B.n[1]) + "×" + std::to_string(B.n[2]);
  const std::size_t ny = static_cast<std::size_t>(B.n[1]) + 1, nz = static_cast<std::size_t>(B.n[2]) + 1;
  // Vertex v of (I₀ × I₁) × I₂ is (i |I₁₀| + j) |I₂₀| + k.
  std::vector<std::array<int, 3>> ijk(X.vertices);
  B.complex.position.resize(X.vertices);
  for (std::size_t v = 0; v < X.vertices; v++) {
    const std::size_t k = v % nz, j = (v / nz) % ny, i = v / (nz * ny);
    ijk[v] = {static_cast<int>(i), static_cast<int>(j), static_cast<int>(k)};
    B.complex.position[v] = {B.coordinates[0][i], B.coordinates[1][j], B.coordinates[2][k]};
  }
  // Each cell's axes and least vertex, from its vertices.
  auto cell_of = [&](const std::vector<std::size_t>& vs) {
    std::array<int, 3> lo{1 << 30, 1 << 30, 1 << 30}, hi{-1, -1, -1};
    for (std::size_t v : vs)
      for (int a = 0; a < 3; a++) {
        lo[static_cast<std::size_t>(a)] = std::min(lo[static_cast<std::size_t>(a)], ijk[v][static_cast<std::size_t>(a)]);
        hi[static_cast<std::size_t>(a)] = std::max(hi[static_cast<std::size_t>(a)], ijk[v][static_cast<std::size_t>(a)]);
      }
    BoxComplex::Cell c;
    c.corner = lo;
    for (int a = 0; a < 3; a++)
      if (hi[static_cast<std::size_t>(a)] > lo[static_cast<std::size_t>(a)]) c.axes = static_cast<std::uint8_t>(c.axes | (1u << a));
    return c;
  };
  for (std::size_t v = 0; v < X.vertices; v++) B.cells[0].push_back(cell_of({v}));
  for (const auto& e : X.edges) B.cells[1].push_back(cell_of({e[0], e[1]}));
  std::vector<std::vector<std::size_t>> face_vertices;
  for (const auto& f : X.faces) {
    face_vertices.push_back(geometric_detail::word_vertices(X, f));
    B.cells[2].push_back(cell_of(face_vertices.back()));
  }
  if (!X.higher.empty())
    for (const auto& cell : X.higher[0]) {
      std::vector<std::size_t> vs;
      for (const auto& e : cell) vs.insert(vs.end(), face_vertices[e.index].begin(), face_vertices[e.index].end());
      B.cells[3].push_back(cell_of(vs));
    }
  B.complex.topology = std::move(X);
  B.complex.frame = std::move(frame);
  B.index();
  return B;
}

// The box [lower, upper] cut into n equal cells along each axis.
inline BoxComplex uniform_box(const Vec3& lower, const Vec3& upper, const std::array<int, 3>& n, std::string frame = "box") {
  std::array<std::vector<double>, 3> c;
  for (int a = 0; a < 3; a++) {
    const int m = n[static_cast<std::size_t>(a)];
    if (m < 0) throw std::invalid_argument("physis::discretisation::uniform_box: no negative counts");
    for (int i = 0; i <= m; i++)
      c[static_cast<std::size_t>(a)].push_back(m == 0 ? lower[static_cast<std::size_t>(a)]
                                                      : (i == m ? upper[static_cast<std::size_t>(a)] : lower[static_cast<std::size_t>(a)] + (upper[static_cast<std::size_t>(a)] - lower[static_cast<std::size_t>(a)]) * i / m));
  }
  return box_complex(std::move(c), std::move(frame));
}

// Each interval of the box cut into r equal parts.
inline BoxComplex refine(const BoxComplex& b, int r) {
  if (r < 1) throw std::invalid_argument("physis::discretisation::refine: r ≥ 1");
  std::array<std::vector<double>, 3> c;
  for (int a = 0; a < 3; a++) {
    const auto& x = b.coordinates[static_cast<std::size_t>(a)];
    c[static_cast<std::size_t>(a)].push_back(x[0]);
    for (std::size_t i = 1; i < x.size(); i++)
      for (int s = 1; s <= r; s++) c[static_cast<std::size_t>(a)].push_back(s == r ? x[i] : x[i - 1] + (x[i] - x[i - 1]) * s / r);
  }
  return box_complex(std::move(c), b.complex.frame);
}

}  // namespace physis::discretisation
