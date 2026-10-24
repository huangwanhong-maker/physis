// Checks of points on the sphere (sphere.hpp).
//
// The six vertices of the octahedron give its 8 faces; the eight of the cube, whose faces put
// four points on each of six circles, give 12 triangles, a valid choice among the degenerate
// ones. On 500 random points: a 2-sphere (F = 2V − 4, E = 3V − 6, every edge on two triangles),
// every triangle anticlockwise from outside and empty of other points beyond its plane (the
// Delaunay property), the areas adding up to 4π. Points within one hemisphere do not cover the
// sphere; duplicate points are merged. Distances and areas against closed forms: a quarter of
// the equator is π/2, the octant triangle has area π/2.

#include <physis/computational_geometry/sphere.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <map>
#include <numbers>
#include <string>
#include <vector>

namespace pc = physis::computational_geometry;

namespace {

bool closed_sphere(const pc::SphericalTriangulation& T) {
  std::map<std::pair<std::size_t, std::size_t>, int> directed;
  for (const auto& t : T.triangles)
    for (int k = 0; k < 3; k++) directed[{t[k], t[(k + 1) % 3]}]++;
  for (const auto& [e, c] : directed)
    if (c != 1 || directed.count({e.second, e.first}) == 0) return false;
  return true;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  constexpr double pi = std::numbers::pi;
  const pc::Vec3 x{1, 0, 0}, y{0, 1, 0}, zz{0, 0, 1};
  s.check("a quarter of the equator is π/2", pc::angle(x, y) - pi / 2, 1e-15);
  s.check("the octant triangle has area π/2", pc::spherical_area(x, y, zz) - pi / 2, 1e-14);

  // The octahedron and the cube.
  {
    const std::vector<pc::Vec3> oct{{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
    const auto T = pc::spherical_delaunay(oct);
    s.require("the octahedron's six vertices: its 8 faces, a closed sphere", T.triangles.size() == 8 && T.edges.size() == 12 && closed_sphere(T) && T.covers_sphere);
    std::vector<pc::Vec3> cube;
    const double c = 1 / std::sqrt(3.0);
    for (int i = 0; i < 8; i++) cube.push_back({(i & 1 ? c : -c), (i & 2 ? c : -c), (i & 4 ? c : -c)});
    const auto C = pc::spherical_delaunay(cube);
    s.require("the cube's eight vertices (four on each face circle): 12 triangles, a closed sphere", C.triangles.size() == 12 && C.edges.size() == 18 && closed_sphere(C));
  }
  // Random points.
  {
    calx::random::Stream z(20261024, 0x73706865u, 0, 0);  // "sphe"
    std::vector<pc::Vec3> pts;
    for (int i = 0; i < 500; i++) {
      const double u = 2 * z.uniform() - 1, phi = 2 * pi * z.uniform();
      pts.push_back(pc::unit_vector(std::asin(u), phi));
    }
    const auto T = pc::spherical_delaunay(pts);
    const std::size_t V = T.points.size();
    s.require("500 random points: F = 2V − 4 and E = 3V − 6", T.triangles.size() == 2 * V - 4 && T.edges.size() == 3 * V - 6,
              std::to_string(T.triangles.size()) + " triangles, " + std::to_string(T.edges.size()) + " edges");
    s.require("... every edge on two triangles, oppositely oriented", closed_sphere(T));
    double area = 0, worst = 0;
    bool anticlockwise = true;
    for (const auto& t : T.triangles) {
      const double a = pc::spherical_area(T.points[t[0]], T.points[t[1]], T.points[t[2]]);
      anticlockwise = anticlockwise && a > 0;
      area += a;
      const pc::Vec3 n = pc::cross({T.points[t[1]][0] - T.points[t[0]][0], T.points[t[1]][1] - T.points[t[0]][1], T.points[t[1]][2] - T.points[t[0]][2]},
                                   {T.points[t[2]][0] - T.points[t[0]][0], T.points[t[2]][1] - T.points[t[0]][1], T.points[t[2]][2] - T.points[t[0]][2]});
      for (const auto& p : T.points) {
        const double side = pc::dot(n, {p[0] - T.points[t[0]][0], p[1] - T.points[t[0]][1], p[2] - T.points[t[0]][2]});
        worst = std::max(worst, side);
      }
    }
    s.require("... every triangle anticlockwise seen from outside", anticlockwise);
    s.check("... no point beyond any triangle's plane (the Delaunay property)", worst, 1e-8);
    s.check("... the areas add up to 4π", area - 4 * pi, 1e-9);
    // One hemisphere, and duplicates.
    std::vector<pc::Vec3> north;
    for (const auto& p : pts)
      if (p[2] > 0.1) north.push_back(p);
    s.require("points in one hemisphere do not cover the sphere", !pc::spherical_delaunay(north).covers_sphere);
    auto twice = pts;
    twice.push_back(pts[7]);
    const auto D = pc::spherical_delaunay(twice);
    s.require("a duplicate point is merged with its twin", D.points.size() == V && D.point_of.back() == D.point_of[7]);
  }
  return s.done();
}
