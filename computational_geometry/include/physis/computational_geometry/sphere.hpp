// Points on the unit sphere: their spherical Delaunay triangulation, distances and areas.
//
// The spherical Delaunay triangulation of points on S² is the boundary of their convex hull in
// ℝ³: a triangle of the points is Delaunay (its circumcircle on the sphere holds no other point)
// exactly when its plane has every other point on the side of the centre (Brown 1979, Inf.
// Process. Lett. 9, 223; Renka 1997, ACM TOMS 23, 416). When the points do not all lie in one
// hemisphere the hull contains the centre, and its triangles tile the sphere: a 2-sphere with
// F = 2V − 4 triangles and E = 3V − 6 edges, whose spherical areas add up to 4π.
//
// The hull is built incrementally (Clarkson and Shor 1989; de Berg et al. 2008, ch. 11): a first
// tetrahedron, then each point in turn removes the triangles it sees and joins itself to the
// horizon. O(n²) in the worst case, plenty for the few thousand places a map of the world needs.
// Orientation tests are in long double on the points moved by a deterministic perturbation of
// 10⁻⁹ (simulation of simplicity in its simplest form, Edelsbrunner and Mücke 1990): four or more
// points on one circle, where the Delaunay triangulation is not unique, are triangulated in one
// of the valid ways, and no point is lost to a zero test. Points within 10⁻⁹ rad of one already
// taken are merged with it.
//
// Positions are geocentric latitude and longitude on a sphere, in radians: a map's places, not a
// geodetic survey (the geodesy domain holds the ellipsoid).

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <numbers>
#include <stdexcept>
#include <utility>
#include <vector>

namespace physis::computational_geometry {

using Vec3 = std::array<double, 3>;

inline Vec3 unit_vector(double latitude, double longitude) {
  return {std::cos(latitude) * std::cos(longitude), std::cos(latitude) * std::sin(longitude), std::sin(latitude)};
}

inline double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
inline Vec3 cross(const Vec3& a, const Vec3& b) { return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}; }

// The great-circle angle between two unit vectors, accurate at every separation.
inline double angle(const Vec3& a, const Vec3& b) { return std::atan2(std::sqrt(dot(cross(a, b), cross(a, b))), dot(a, b)); }

// The area of the spherical triangle (a, b, c) on the unit sphere, signed by orientation:
// tan(E/2) = a · (b × c) / (1 + a·b + b·c + c·a) (Van Oosterom and Strackee 1983, IEEE Trans.
// Biomed. Eng. 30, 125).
inline double spherical_area(const Vec3& a, const Vec3& b, const Vec3& c) {
  return 2 * std::atan2(dot(a, cross(b, c)), 1 + dot(a, b) + dot(b, c) + dot(c, a));
}

struct SphericalTriangulation {
  std::vector<Vec3> points;                          // the points taken (merged ones once)
  std::vector<std::size_t> point_of;                 // each input point -> its point
  std::vector<std::array<std::size_t, 3>> triangles; // anticlockwise seen from outside
  std::vector<std::array<std::size_t, 2>> edges;     // u < v
  bool covers_sphere = false;                        // the hull contains the centre
};

namespace sphere_detail {

inline long double orient(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
  // The signed volume of (b − a, c − a, d − a): positive when d is on the side the anticlockwise
  // triangle (a, b, c) faces.
  const long double bx = b[0] - a[0], by = b[1] - a[1], bz = b[2] - a[2];
  const long double cx = c[0] - a[0], cy = c[1] - a[1], cz = c[2] - a[2];
  const long double dx = d[0] - a[0], dy = d[1] - a[1], dz = d[2] - a[2];
  return bx * (cy * dz - cz * dy) - by * (cx * dz - cz * dx) + bz * (cx * dy - cy * dx);
}

}  // namespace sphere_detail

inline SphericalTriangulation spherical_delaunay(const std::vector<Vec3>& input) {
  using sphere_detail::orient;
  const double tolerance = 0;
  SphericalTriangulation T;
  for (const auto& p : input) {
    std::size_t same = T.points.size();
    for (std::size_t q = 0; q < T.points.size(); q++)
      if (angle(p, T.points[q]) < 1e-9) { same = q; break; }
    if (same == T.points.size()) T.points.push_back(p);
    T.point_of.push_back(same);
  }
  // The perturbed points on which every test is made.
  std::vector<Vec3> P = T.points;
  for (std::size_t i = 0; i < P.size(); i++) {
    std::uint64_t h = 0x9e3779b97f4a7c15ull * (i + 1);
    for (int k = 0; k < 3; k++) {
      h ^= h >> 31;
      h *= 0xbf58476d1ce4e5b9ull;
      h ^= h >> 29;
      P[i][k] += 1e-9 * (static_cast<double>(h >> 11) / 9007199254740992.0 - 0.5);
    }
    const double r = std::sqrt(dot(P[i], P[i]));
    for (double& x : P[i]) x /= r;
  }
  const std::size_t n = P.size();
  if (n < 4) throw std::invalid_argument("physis::computational_geometry::spherical_delaunay: at least four distinct points");
  // A first tetrahedron, well conditioned: the point farthest from the line of the first two, then
  // the point farthest from their plane.
  std::size_t i0 = 0, i1 = 1, i2 = n, i3 = n;
  double best = 0;
  for (std::size_t k = 2; k < n; k++) {
    const Vec3 u{P[i1][0] - P[i0][0], P[i1][1] - P[i0][1], P[i1][2] - P[i0][2]}, w{P[k][0] - P[i0][0], P[k][1] - P[i0][1], P[k][2] - P[i0][2]};
    const Vec3 c = cross(u, w);
    if (dot(c, c) > best) { best = dot(c, c); i2 = k; }
  }
  long double most = 0;
  for (std::size_t k = 2; k < n && i2 != n; k++) {
    const long double o = orient(P[i0], P[i1], P[i2], P[k]);
    if (k != i2 && (o < 0 ? -o : o) > most) { most = o < 0 ? -o : o; i3 = k; }
  }
  if (i2 == n || i3 == n || most < 1e-15) throw std::invalid_argument("physis::computational_geometry::spherical_delaunay: the points lie in one plane");
  std::vector<std::array<std::size_t, 3>> faces;
  std::vector<char> alive;
  auto add = [&](std::size_t a, std::size_t b, std::size_t c) {
    faces.push_back({a, b, c});
    alive.push_back(1);
  };
  if (orient(P[i0], P[i1], P[i2], P[i3]) > 0) std::swap(i1, i2);  // (i0, i1, i2) faces away from i3
  add(i0, i1, i2);
  add(i0, i3, i1);
  add(i1, i3, i2);
  add(i2, i3, i0);
  for (std::size_t p = 0; p < n; p++) {
    if (p == i0 || p == i1 || p == i2 || p == i3) continue;
    std::vector<std::size_t> visible;
    for (std::size_t f = 0; f < faces.size(); f++)
      if (alive[f] && orient(P[faces[f][0]], P[faces[f][1]], P[faces[f][2]], P[p]) > tolerance) visible.push_back(f);
    if (visible.empty()) continue;  // inside the hull (not on the sphere): no triangle sees it
    // The horizon: directed edges of visible faces whose reverse is not on a visible face.
    std::map<std::pair<std::size_t, std::size_t>, int> count;
    for (std::size_t f : visible) {
      const auto& t = faces[f];
      for (int k = 0; k < 3; k++) count[{t[k], t[(k + 1) % 3]}]++;
      alive[f] = 0;
    }
    for (const auto& [e, c] : count)
      if (!count.count({e.second, e.first})) add(e.first, e.second, p);
  }
  for (std::size_t f = 0; f < faces.size(); f++)
    if (alive[f]) T.triangles.push_back(faces[f]);
  std::map<std::pair<std::size_t, std::size_t>, int> edges;
  for (const auto& t : T.triangles)
    for (int k = 0; k < 3; k++) edges[std::minmax(t[k], t[(k + 1) % 3])]++;
  for (const auto& [e, c] : edges) T.edges.push_back({e.first, e.second});
  // The hull contains the centre when the origin is inside every face's plane.
  T.covers_sphere = true;
  const Vec3 origin{0, 0, 0};
  for (const auto& t : T.triangles) T.covers_sphere = T.covers_sphere && orient(P[t[0]], P[t[1]], P[t[2]], origin) < -1e-12;
  return T;
}

}  // namespace physis::computational_geometry
