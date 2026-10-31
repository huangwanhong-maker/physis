// The faces of a plane graph whose edges are polylines: the rotation system read from the directions in
// which edges leave each vertex, the faces as the orbits of the face permutation, their signed areas, and
// which of them are bounded by a simple circuit — the 2-cells of a regular cell structure (decision 0018).
//
// The construction (Edmonds 1960, Notices AMS 7, 646; Mohar and Thomassen 2001, Graphs on Surfaces, §3.2). Each
// edge e = (u, v) gives two darts, u → v and v → u. At each vertex the darts leaving it are ordered
// counter-clockwise by the direction of their first segment (the rotation σ). The face on the left of a dart
// u → v continues with the dart leaving v that comes just before v → u in that order (the next one clockwise
// from it): φ(d) = σ⁻¹(α(d)), α reversing a dart. The faces are the orbits of φ. For a plane embedding — edges
// meeting only at their ends, which is assumed, not checked — each component is a map on the sphere, so
// V − E + F = 2c with F the orbits, one of them each component's outer boundary walk (Euler's formula
// V − E + F = 1 + c counts a single unbounded face instead); each bounded face is traversed
// counter-clockwise, with positive shoelace area over its polylines, each outer walk clockwise, negative.
//
// A face is simple when its boundary walk visits no vertex twice. A simple bounded face is a disc whose
// boundary is a circuit, so it can be a 2-cell of a regular complex (morphe regular_complex.hpp); a face with
// a dangling edge or a bridge inside it is not, and stays a hole. With them, the polylines' plane geometry: length,
// the point at an arc length, a cut, the nearest point, a ring's area, crossings and containment. Cost O(E log E). Checks:
// tests/check_planar.cpp.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace physis::computational_geometry {

using Point2 = std::array<double, 2>;

struct PlaneEdge {
  std::size_t u = 0, v = 0;
  std::vector<Point2> line;  // from u's point to v's point, both included; empty for a straight segment
};

struct PlaneFace {
  std::vector<std::size_t> vertices;               // in order round the face
  std::vector<std::pair<std::size_t, int>> edges;  // (edge, +1 when walked from u to v, −1 otherwise)
  double area = 0;                                 // signed: > 0 for bounded faces
  bool simple = false;                             // no vertex visited twice
};

// ---- polylines -------------------------------------------------------------------------------------------------

inline double polyline_length(const std::vector<Point2>& l) {
  double s = 0;
  for (std::size_t k = 0; k + 1 < l.size(); k++) s += std::hypot(l[k + 1][0] - l[k][0], l[k + 1][1] - l[k][1]);
  return s;
}
// The point at arc length s from the first point (clamped to the ends).
inline Point2 point_at(const std::vector<Point2>& l, double s) {
  if (l.empty()) throw std::invalid_argument("physis::computational_geometry::point_at: an empty polyline");
  for (std::size_t k = 0; k + 1 < l.size(); k++) {
    const double d = std::hypot(l[k + 1][0] - l[k][0], l[k + 1][1] - l[k][1]);
    if (s <= d && d > 0) return {l[k][0] + (l[k + 1][0] - l[k][0]) * s / d, l[k][1] + (l[k + 1][1] - l[k][1]) * s / d};
    s -= d;
  }
  return l.back();
}
// The two pieces of a polyline cut at arc length s, the cut point ending the first and starting the second.
inline std::pair<std::vector<Point2>, std::vector<Point2>> split_polyline(const std::vector<Point2>& l, double s) {
  const double total = polyline_length(l);
  if (!(s > 0 && s < total)) throw std::invalid_argument("physis::computational_geometry::split_polyline: a cut strictly inside the polyline");
  std::vector<Point2> a{l[0]}, b;
  double run = 0;
  std::size_t k = 0;
  for (; k + 1 < l.size(); k++) {
    const double d = std::hypot(l[k + 1][0] - l[k][0], l[k + 1][1] - l[k][1]);
    if (run + d >= s && d > 0) break;
    run += d;
    a.push_back(l[k + 1]);
  }
  const Point2 p = point_at(l, s);
  a.push_back(p);
  b.push_back(p);
  for (std::size_t q = k + 1; q < l.size(); q++) b.push_back(l[q]);
  return {a, b};
}
// The arc length of the point of the polyline nearest to p, and the distance to it.
inline std::pair<double, double> nearest_on_polyline(const std::vector<Point2>& l, Point2 p) {
  double best = std::numeric_limits<double>::infinity(), at = 0, run = 0;
  for (std::size_t k = 0; k + 1 < l.size(); k++) {
    const double dx = l[k + 1][0] - l[k][0], dy = l[k + 1][1] - l[k][1], d2 = dx * dx + dy * dy;
    double u = d2 > 0 ? ((p[0] - l[k][0]) * dx + (p[1] - l[k][1]) * dy) / d2 : 0;
    u = std::clamp(u, 0.0, 1.0);
    const double dist = std::hypot(l[k][0] + u * dx - p[0], l[k][1] + u * dy - p[1]);
    if (dist < best) best = dist, at = run + u * std::sqrt(d2);
    run += std::sqrt(d2);
  }
  return {at, best};
}
// The signed area of a closed ring (counter-clockwise positive), by the shoelace formula.
inline double ring_area(const std::vector<Point2>& r) {
  double twice = 0;
  for (std::size_t k = 0; k < r.size(); k++) {
    const auto& a = r[k];
    const auto& b = r[(k + 1) % r.size()];
    twice += a[0] * b[1] - b[0] * a[1];
  }
  return twice / 2;
}
// Whether the open segments ab and cd cross at a point interior to both (touching at ends does not count).
inline bool segments_cross(Point2 a, Point2 b, Point2 c, Point2 d) {
  const auto orient = [](Point2 p, Point2 q, Point2 r) { return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]); };
  const double o1 = orient(a, b, c), o2 = orient(a, b, d), o3 = orient(c, d, a), o4 = orient(c, d, b);
  return ((o1 > 0 && o2 < 0) || (o1 < 0 && o2 > 0)) && ((o3 > 0 && o4 < 0) || (o3 < 0 && o4 > 0));
}
// Whether a point lies strictly inside a ring (even–odd rule).
inline bool inside_ring(const std::vector<Point2>& r, Point2 p) {
  bool in = false;
  for (std::size_t i = 0, j = r.size() - 1; i < r.size(); j = i++)
    if ((r[i][1] > p[1]) != (r[j][1] > p[1]) && p[0] < (r[j][0] - r[i][0]) * (p[1] - r[i][1]) / (r[j][1] - r[i][1]) + r[i][0]) in = !in;
  return in;
}

namespace detail {

inline std::vector<Point2> polyline(const std::vector<Point2>& points, const PlaneEdge& e) {
  if (!e.line.empty()) return e.line;
  return {points.at(e.u), points.at(e.v)};
}

}  // namespace detail

// The faces of the plane graph; isolated vertices make no face.
inline std::vector<PlaneFace> plane_faces(const std::vector<Point2>& points, const std::vector<PlaneEdge>& edges) {
  const std::size_t E = edges.size();
  // Dart 2e is u → v, dart 2e + 1 is v → u.
  std::vector<std::vector<Point2>> geom(2 * E);
  std::vector<std::size_t> tail(2 * E), head(2 * E);
  for (std::size_t e = 0; e < E; e++) {
    const auto& ed = edges[e];
    if (ed.u >= points.size() || ed.v >= points.size()) throw std::invalid_argument("physis::computational_geometry::plane_faces: an edge leaves the points");
    if (ed.u == ed.v) throw std::invalid_argument("physis::computational_geometry::plane_faces: a loop");
    auto line = detail::polyline(points, ed);
    if (line.size() < 2) throw std::invalid_argument("physis::computational_geometry::plane_faces: a polyline of two points at least");
    geom[2 * e] = line;
    std::reverse(line.begin(), line.end());
    geom[2 * e + 1] = std::move(line);
    tail[2 * e] = ed.u, head[2 * e] = ed.v;
    tail[2 * e + 1] = ed.v, head[2 * e + 1] = ed.u;
  }
  // The rotation at each vertex: darts leaving it by the angle of their first segment.
  const auto angle = [&](std::size_t d) {
    const auto& g = geom[d];
    std::size_t k = 1;
    while (k + 1 < g.size() && g[k][0] == g[0][0] && g[k][1] == g[0][1]) k++;
    return std::atan2(g[k][1] - g[0][1], g[k][0] - g[0][0]);
  };
  std::vector<std::vector<std::size_t>> around(points.size());
  for (std::size_t d = 0; d < 2 * E; d++) around[tail[d]].push_back(d);
  std::vector<std::size_t> position(2 * E);
  for (auto& list : around) {
    std::sort(list.begin(), list.end(), [&](std::size_t a, std::size_t b) { return angle(a) < angle(b); });
    for (std::size_t k = 0; k < list.size(); k++) position[list[k]] = k;
  }
  // φ(d): at v = head(d), the dart before the reverse of d in counter-clockwise order.
  const auto next = [&](std::size_t d) {
    const std::size_t r = d ^ 1, v = tail[r];
    const auto& list = around[v];
    const std::size_t k = position[r];
    return list[(k + list.size() - 1) % list.size()];
  };
  std::vector<char> seen(2 * E, 0);
  std::vector<PlaneFace> faces;
  for (std::size_t d0 = 0; d0 < 2 * E; d0++) {
    if (seen[d0]) continue;
    PlaneFace f;
    double twice = 0;
    std::size_t d = d0;
    do {
      seen[d] = 1;
      f.vertices.push_back(tail[d]);
      f.edges.push_back({d / 2, d % 2 == 0 ? 1 : -1});
      const auto& g = geom[d];
      for (std::size_t k = 0; k + 1 < g.size(); k++) twice += g[k][0] * g[k + 1][1] - g[k + 1][0] * g[k][1];
      d = next(d);
    } while (d != d0);
    f.area = twice / 2;
    std::vector<std::size_t> v = f.vertices;
    std::sort(v.begin(), v.end());
    f.simple = std::adjacent_find(v.begin(), v.end()) == v.end();
    faces.push_back(std::move(f));
  }
  return faces;
}

}  // namespace physis::computational_geometry
