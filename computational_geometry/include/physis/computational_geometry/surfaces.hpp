// Triangle surfaces in space: their points welded, their area and bounds, their sections by a plane as
// oriented loops, and the labelled points that come within a distance of one another (decision 0019).
//
// A triangle surface is a list of points and of triangles, three indices each, oriented by the order of their
// corners: the normal n = (b − a) × (c − a) points to the side from which a, b, c turn counter-clockwise. A
// closed surface oriented outward bounds a solid; its mass properties are physis::mechanics::solid_polyhedron,
// its topology morphe's (topology/surfaces.hpp, on the simplicial complex of its triangles).
//
// Welding. Points closer than a tolerance ε are made one: the points are hashed to cells of side ε, each is
// compared with those already kept in its own and the 26 neighbouring cells, and joins the first within ε.
// The relation is not transitive, so the result depends on the order of the points when they lie in chains
// spaced by less than ε; for the seams of a mesh, whose duplicates coincide to rounding, it does not.
// Triangles left with a repeated corner are dropped.
//
// Sections. The plane P = {x : N·(x − o) = 0}, N a unit normal, puts each corner on one side by the sign of
// s = N·(x − o), a corner with s = 0 counted on the positive side (a symbolic perturbation: the plane is
// moved by an infinitesimal toward −N). Then a triangle meets the plane in no segment or in exactly one,
// joining the crossing points of two of its edges, and a crossing point is named by its edge, so that
// segments are joined into chains exactly, without a tolerance: at an edge shared by two triangles, one
// segment ends where the next begins. The direction of a segment is read from the signs of its corners and
// their order, not from the crossing points, which coincide where a corner lies on the plane. On a closed 2-manifold every chain closes into a loop (Minetto et al.
// 2017, An optimal algorithm for 3D triangle mesh slicing, Computer-Aided Design 92, 1, §3). Each segment is
// directed along N × n. For a closed surface oriented outward, N × (N × n) = −n_⊥ is the inward direction in
// the plane, so the solid lies on the left of every segment seen from +N: the loops go counter-clockwise
// round the pieces of the section and clockwise round its holes, and the area of the section is the sum of
// the loops' shoelace areas in any basis (u, v) of the plane with u × v = N (Green's theorem). The section's
// area integrates over the planes' offsets to the solid's volume, ∫ A(h) dh = V (Cavalieri's principle);
// the checks hold the two together. On a surface that is not closed, or not a manifold where it is cut, a
// chain may stop at an edge with one triangle; it is returned open.
//
// Close pairs. Of points carrying labels, the pairs of labels some of whose points lie within a distance d,
// with the least distance between them: the points are sorted into cells of side d, and each cell compared
// with itself and its 13 forward neighbours, so that every pair of points within d is examined once. The
// cost is linear in the points times the points of other labels in a neighbourhood. Checks:
// tests/check_surfaces.cpp.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include <physis/computational_geometry/planar.hpp>

namespace physis::computational_geometry {

using Point3 = std::array<double, 3>;
using Triangle = std::array<std::size_t, 3>;

struct TriangleSurface {
  std::vector<Point3> points;
  std::vector<Triangle> triangles;
};

struct Box3 {
  Point3 lower{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
  Point3 upper{-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
  bool empty() const { return !(lower[0] <= upper[0]); }
};

namespace detail {
inline Point3 sub3(const Point3& a, const Point3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
inline Point3 cross3(const Point3& a, const Point3& b) { return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}; }
inline double dot3(const Point3& a, const Point3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
inline double norm3(const Point3& a) { return std::sqrt(dot3(a, a)); }

// The cell of a grid of side h holding a point, and a hash of cells for unordered maps.
using Cell = std::array<std::int64_t, 3>;
inline Cell cell_of(const Point3& p, double h) {
  return {static_cast<std::int64_t>(std::floor(p[0] / h)), static_cast<std::int64_t>(std::floor(p[1] / h)), static_cast<std::int64_t>(std::floor(p[2] / h))};
}
struct CellHash {
  std::size_t operator()(const Cell& c) const {
    std::uint64_t h = 1469598103934665603ull;
    for (auto x : c) h = (h ^ static_cast<std::uint64_t>(x)) * 1099511628211ull;
    return static_cast<std::size_t>(h ^ (h >> 29));
  }
};
}  // namespace detail

// ---- measures ------------------------------------------------------------------------------------------------

inline Box3 bounds(const std::vector<Point3>& points) {
  Box3 b;
  for (const auto& p : points)
    for (int i = 0; i < 3; i++) b.lower[i] = std::min(b.lower[i], p[i]), b.upper[i] = std::max(b.upper[i], p[i]);
  return b;
}

// The normal (b − a) × (c − a) of a triangle, of length twice its area.
inline Point3 triangle_normal(const TriangleSurface& S, const Triangle& t) {
  return detail::cross3(detail::sub3(S.points[t[1]], S.points[t[0]]), detail::sub3(S.points[t[2]], S.points[t[0]]));
}

inline double surface_area(const TriangleSurface& S) {
  double a = 0;
  for (const auto& t : S.triangles) a += detail::norm3(triangle_normal(S, t));
  return a / 2;
}

// The points of S within ε of one another made one (see above); `map`, when given, receives each old point's new index.
inline TriangleSurface welded(const TriangleSurface& S, double epsilon, std::vector<std::size_t>* map = nullptr) {
  if (!(epsilon > 0)) throw std::invalid_argument("physis::computational_geometry::welded: a positive tolerance");
  TriangleSurface W;
  std::vector<std::size_t> to(S.points.size());
  std::unordered_map<detail::Cell, std::vector<std::size_t>, detail::CellHash> cells;
  for (std::size_t p = 0; p < S.points.size(); p++) {
    const auto c = detail::cell_of(S.points[p], epsilon);
    std::size_t found = SIZE_MAX;
    for (int di = -1; di <= 1 && found == SIZE_MAX; di++)
      for (int dj = -1; dj <= 1 && found == SIZE_MAX; dj++)
        for (int dk = -1; dk <= 1 && found == SIZE_MAX; dk++) {
          const auto it = cells.find(detail::Cell{c[0] + di, c[1] + dj, c[2] + dk});
          if (it == cells.end()) continue;
          for (std::size_t q : it->second)
            if (detail::norm3(detail::sub3(W.points[q], S.points[p])) <= epsilon) { found = q; break; }
        }
    if (found == SIZE_MAX) {
      found = W.points.size();
      W.points.push_back(S.points[p]);
      cells[c].push_back(found);
    }
    to[p] = found;
  }
  for (const auto& t : S.triangles) {
    const Triangle u{to.at(t[0]), to.at(t[1]), to.at(t[2])};
    if (u[0] != u[1] && u[1] != u[2] && u[0] != u[2]) W.triangles.push_back(u);
  }
  if (map) *map = std::move(to);
  return W;
}

// ---- sections by a plane -------------------------------------------------------------------------------------

struct Plane3 {
  Point3 origin{};
  Point3 normal{0, 0, 1};  // of unit length
};

struct SectionChain {
  std::vector<Point3> points;  // in order; a closed loop does not repeat its first point
  bool closed = false;
};

// The section of S by the plane, as chains of crossing points directed along N × n (see above).
inline std::vector<SectionChain> plane_section(const TriangleSurface& S, const Plane3& P) {
  const double nn = detail::norm3(P.normal);
  if (!(std::fabs(nn - 1) < 1e-9)) throw std::invalid_argument("physis::computational_geometry::plane_section: a unit normal");
  std::vector<double> s(S.points.size());
  for (std::size_t p = 0; p < S.points.size(); p++) s[p] = detail::dot3(P.normal, detail::sub3(S.points[p], P.origin));
  const auto side = [&](std::size_t p) { return s[p] >= 0; };
  using Edge = std::pair<std::size_t, std::size_t>;
  const auto edge = [](std::size_t a, std::size_t b) { return a < b ? Edge{a, b} : Edge{b, a}; };
  const auto crossing = [&](const Edge& e) {
    const double t = s[e.first] / (s[e.first] - s[e.second]);
    const auto& a = S.points[e.first];
    const auto& b = S.points[e.second];
    return Point3{a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2])};
  };
  // Each segment, from the crossing of one edge to the crossing of another.
  std::map<Edge, std::vector<Edge>> next;  // from an edge, the edges its segments go to
  std::map<Edge, int> incoming;
  for (const auto& t : S.triangles) {
    const bool a = side(t[0]), b = side(t[1]), c = side(t[2]);
    if (a == b && b == c) continue;
    // The corner alone on its side, i; the segment joins the edges (i, i+1) and (i−1, i). Its direction along
    // N × n is fixed by the pattern of signs and the triangle's order alone, which keeps it exact where the
    // crossing points coincide (a corner on the plane): from (i, i+1) to (i−1, i) when i is on the positive
    // side, the other way when it is on the negative side.
    const int i = (a != b && a != c) ? 0 : (b != a && b != c) ? 1 : 2;
    const std::size_t ti = t[i], tn = t[(i + 1) % 3], tp = t[(i + 2) % 3];
    Edge from = edge(ti, tn), to = edge(tp, ti);
    if (!side(ti)) std::swap(from, to);
    next[from].push_back(to);
    incoming[to]++;
    incoming.try_emplace(from, 0);
  }
  std::vector<SectionChain> chains;
  std::map<Edge, std::size_t> used;  // segments taken from each edge
  auto follow = [&](Edge start) {
    SectionChain c;
    Edge e = start;
    while (true) {
      c.points.push_back(crossing(e));
      auto it = next.find(e);
      if (it == next.end() || used[e] >= it->second.size()) break;
      const Edge f = it->second[used[e]++];
      if (f == start) { c.closed = true; break; }
      e = f;
    }
    return c;
  };
  // Open chains first, from edges where more segments leave than arrive; then the loops.
  for (const auto& [e, out] : next)
    for (std::size_t k = incoming[e]; k < out.size(); k++)
      if (used[e] < out.size()) chains.push_back(follow(e));
  for (const auto& [e, out] : next)
    while (used[e] < out.size()) chains.push_back(follow(e));
  return chains;
}

// The points of a chain in the plane's basis (u, v): (u·(x − o), v·(x − o)).
inline std::vector<Point2> in_plane(const std::vector<Point3>& points, const Point3& origin, const Point3& u, const Point3& v) {
  std::vector<Point2> out;
  out.reserve(points.size());
  for (const auto& p : points) {
    const Point3 d = detail::sub3(p, origin);
    out.push_back({detail::dot3(u, d), detail::dot3(v, d)});
  }
  return out;
}

// The area of the section: the sum of its closed loops' shoelace areas in a basis with u × v = N.
inline double section_area(const std::vector<SectionChain>& chains, const Plane3& P) {
  Point3 u = std::fabs(P.normal[0]) < 0.9 ? Point3{1, 0, 0} : Point3{0, 1, 0};
  u = detail::cross3(u, P.normal);
  const double l = detail::norm3(u);
  u = {u[0] / l, u[1] / l, u[2] / l};
  const Point3 v = detail::cross3(P.normal, u);
  double a = 0;
  for (const auto& c : chains)
    if (c.closed) a += ring_area(in_plane(c.points, P.origin, u, v));
  return a;
}

// ---- close pairs ---------------------------------------------------------------------------------------------

struct ClosePair {
  std::size_t a = 0, b = 0;  // labels, a < b
  double distance = 0;       // the least distance between their points, ≤ d
};

inline std::vector<ClosePair> close_pairs(const std::vector<Point3>& points, const std::vector<std::size_t>& labels, double d) {
  if (points.size() != labels.size()) throw std::invalid_argument("physis::computational_geometry::close_pairs: a label for each point");
  if (!(d > 0)) throw std::invalid_argument("physis::computational_geometry::close_pairs: a positive distance");
  const std::size_t n = points.size();
  std::vector<detail::Cell> cell(n);
  for (std::size_t p = 0; p < n; p++) cell[p] = detail::cell_of(points[p], d);
  std::vector<std::size_t> order(n);
  for (std::size_t p = 0; p < n; p++) order[p] = p;
  std::sort(order.begin(), order.end(), [&](std::size_t x, std::size_t y) { return cell[x] != cell[y] ? cell[x] < cell[y] : x < y; });
  std::unordered_map<detail::Cell, std::pair<std::size_t, std::size_t>, detail::CellHash> range;  // cell → [begin, end) in order
  for (std::size_t i = 0; i < n;) {
    std::size_t j = i;
    while (j < n && cell[order[j]] == cell[order[i]]) j++;
    range[cell[order[i]]] = {i, j};
    i = j;
  }
  std::map<std::pair<std::size_t, std::size_t>, double> best;
  const double d2 = d * d;
  auto compare = [&](std::size_t p, std::size_t q) {
    if (labels[p] == labels[q]) return;
    const Point3 r = detail::sub3(points[p], points[q]);
    const double r2 = detail::dot3(r, r);
    if (r2 > d2) return;
    const auto k = labels[p] < labels[q] ? std::make_pair(labels[p], labels[q]) : std::make_pair(labels[q], labels[p]);
    auto [it, fresh] = best.try_emplace(k, r2);
    if (!fresh) it->second = std::min(it->second, r2);
  };
  // The 13 forward neighbours: offsets (di, dj, dk) lexicographically after (0, 0, 0).
  std::vector<std::array<int, 3>> forward;
  for (int di = -1; di <= 1; di++)
    for (int dj = -1; dj <= 1; dj++)
      for (int dk = -1; dk <= 1; dk++)
        if (std::array<int, 3>{di, dj, dk} > std::array<int, 3>{0, 0, 0}) forward.push_back({di, dj, dk});
  for (const auto& [k, r] : range) {
    const auto& c = cell[order[r.first]];
    for (std::size_t x = r.first; x < r.second; x++)
      for (std::size_t y = x + 1; y < r.second; y++) compare(order[x], order[y]);
    for (const auto& o : forward) {
      const auto it = range.find(detail::Cell{c[0] + o[0], c[1] + o[1], c[2] + o[2]});
      if (it == range.end()) continue;
      for (std::size_t x = r.first; x < r.second; x++)
        for (std::size_t y = it->second.first; y < it->second.second; y++) compare(order[x], order[y]);
    }
  }
  std::vector<ClosePair> out;
  out.reserve(best.size());
  for (const auto& [k, r2] : best) out.push_back({k.first, k.second, std::sqrt(r2)});
  return out;
}

}  // namespace physis::computational_geometry
