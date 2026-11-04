// Checks of triangle surfaces in space (surfaces.hpp; decision 0019).
//
// - Welding: a cube given as a soup of 36 corners, jittered by 10⁻¹⁰, welds to 8 points and 12 triangles;
//   its area is 6.
// - Sections of the cube: a level plane cuts one counter-clockwise loop of area 1, through the faces' diagonals too; the plane through the top
//   face, moved infinitesimally down by the perturbation, still cuts area 1, and the plane through the bottom
//   face cuts nothing; the plane through the centre normal to (1, 1, 1) cuts the regular hexagon of side
//   √2/2, area 3√3/4.
// - Cavalieri's principle: the sections of a subdivided icosahedron integrate, by the midpoint rule over 4000
//   planes, to the volume of physis::mechanics::solid_polyhedron, along an axis and along an oblique
//   direction; and so do a torus's, whose sections through its hole are a counter-clockwise outer loop and a
//   clockwise inner one, two loops of opposite signs.
// - A surface that is not closed: a square of two triangles cut across gives one open chain of three points.
// - Close pairs: three unit cubes on a line, gaps 0.5 and 2.5; and 2000 random points of seven labels against
//   every pair of points compared directly.

#include <array>
#include <cmath>
#include <map>
#include <numbers>
#include <vector>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <physis/computational_geometry/surfaces.hpp>
#include <physis/mechanics/rigid_body.hpp>

using calx::verification::Suite;
namespace cg = physis::computational_geometry;
namespace pm = physis::mechanics;

namespace {

cg::TriangleSurface cube(double x0 = 0) {
  cg::TriangleSurface S;
  for (int k = 0; k < 8; k++) S.points.push_back({x0 + (k & 1 ? 1.0 : 0.0), k & 2 ? 1.0 : 0.0, k & 4 ? 1.0 : 0.0});
  const int quads[6][4] = {{0, 2, 3, 1}, {4, 5, 7, 6}, {0, 1, 5, 4}, {2, 6, 7, 3}, {0, 4, 6, 2}, {1, 3, 7, 5}};
  for (const auto& q : quads) {
    S.triangles.push_back({std::size_t(q[0]), std::size_t(q[1]), std::size_t(q[2])});
    S.triangles.push_back({std::size_t(q[0]), std::size_t(q[2]), std::size_t(q[3])});
  }
  return S;
}

// The icosahedron subdivided `levels` times, its points pushed onto the unit sphere; oriented outward.
cg::TriangleSurface icosphere(int levels) {
  const double t = (1 + std::sqrt(5.0)) / 2;
  cg::TriangleSurface S;
  S.points = {{-1, t, 0}, {1, t, 0}, {-1, -t, 0}, {1, -t, 0}, {0, -1, t}, {0, 1, t}, {0, -1, -t}, {0, 1, -t}, {t, 0, -1}, {t, 0, 1}, {-t, 0, -1}, {-t, 0, 1}};
  S.triangles = {{0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11}, {1, 5, 9}, {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
                 {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9}, {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1}};
  auto unit = [](cg::Point3 p) { const double l = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]); return cg::Point3{p[0] / l, p[1] / l, p[2] / l}; };
  for (auto& p : S.points) p = unit(p);
  for (int l = 0; l < levels; l++) {
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> mid;
    auto middle = [&](std::size_t a, std::size_t b) {
      const auto k = a < b ? std::make_pair(a, b) : std::make_pair(b, a);
      auto it = mid.find(k);
      if (it != mid.end()) return it->second;
      const auto& p = S.points[a];
      const auto& q = S.points[b];
      S.points.push_back(unit({(p[0] + q[0]) / 2, (p[1] + q[1]) / 2, (p[2] + q[2]) / 2}));
      return mid[k] = S.points.size() - 1;
    };
    std::vector<cg::Triangle> next;
    for (const auto& tr : S.triangles) {
      const std::size_t a = middle(tr[0], tr[1]), b = middle(tr[1], tr[2]), c = middle(tr[2], tr[0]);
      next.push_back({tr[0], a, c});
      next.push_back({tr[1], b, a});
      next.push_back({tr[2], c, b});
      next.push_back({a, b, c});
    }
    S.triangles = std::move(next);
  }
  return S;
}

// A torus about the z axis, radii R and r, m × n quads, the tube's angle offset by half a step so that no
// vertex lies on z = 0; oriented outward.
cg::TriangleSurface torus(double R, double r, int m, int n) {
  cg::TriangleSurface S;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++) {
      const double u = 2 * std::numbers::pi * i / m, v = 2 * std::numbers::pi * (j + 0.5) / n;
      S.points.push_back({(R + r * std::cos(v)) * std::cos(u), (R + r * std::cos(v)) * std::sin(u), r * std::sin(v)});
    }
  auto at = [&](int i, int j) { return std::size_t(((i + m) % m) * n + (j + n) % n); };
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++) {
      S.triangles.push_back({at(i, j), at(i + 1, j), at(i + 1, j + 1)});
      S.triangles.push_back({at(i, j), at(i + 1, j + 1), at(i, j + 1)});
    }
  return S;
}

pm::MassProperties solid(const cg::TriangleSurface& S) {
  std::vector<physis::geometry::Vec3> P(S.points.begin(), S.points.end());
  std::vector<std::array<std::size_t, 3>> T(S.triangles.begin(), S.triangles.end());
  return pm::solid_polyhedron(P, T, 1);
}

// ∫ A(h) dh over the extent of S along the unit direction N, by the midpoint rule on k planes.
double cavalieri(const cg::TriangleSurface& S, const cg::Point3& N, int k) {
  double lo = 1e300, hi = -1e300;
  for (const auto& p : S.points) {
    const double h = N[0] * p[0] + N[1] * p[1] + N[2] * p[2];
    lo = std::min(lo, h), hi = std::max(hi, h);
  }
  const double dh = (hi - lo) / k;
  double v = 0;
  for (int i = 0; i < k; i++) {
    const double h = lo + (i + 0.5) * dh;
    const cg::Plane3 P{{h * N[0], h * N[1], h * N[2]}, N};
    v += cg::section_area(cg::plane_section(S, P), P) * dh;
  }
  return v;
}

}  // namespace

int main() {
  Suite s;
  {
    const auto C = cube();
    cg::TriangleSurface soup;
    calx::random::Stream z(19, 1, 0, 0);
    for (const auto& t : C.triangles) {
      for (int k = 0; k < 3; k++) {
        auto p = C.points[t[k]];
        for (auto& x : p) x += 1e-10 * (z.uniform() - 0.5);
        soup.points.push_back(p);
      }
      const std::size_t b = soup.points.size() - 3;
      soup.triangles.push_back({b, b + 1, b + 2});
    }
    std::vector<std::size_t> map;
    const auto W = cg::welded(soup, 1e-8, &map);
    s.require("a cube as a soup of 36 corners welds to 8 points and 12 triangles", W.points.size() == 8 && W.triangles.size() == 12 && map.size() == 36);
    s.check("its area is 6", cg::surface_area(W) - 6, 1e-8);
    const auto B = cg::bounds(C.points);
    s.require("its bounds are the unit cube", B.lower == cg::Point3{0, 0, 0} && B.upper == cg::Point3{1, 1, 1});
  }
  {
    const auto C = cube();
    const cg::Plane3 mid{{0, 0, 0.3}, {0, 0, 1}};
    const auto L = cg::plane_section(C, mid);
    s.require("the cube cut at z = 0.3: one closed loop, through the four sides and the diagonals that cut them", L.size() == 1 && L[0].closed && L[0].points.size() == 8);
    s.check("of area +1 (counter-clockwise seen from +z)", cg::section_area(L, mid) - 1, 1e-15);
    const cg::Plane3 top{{0, 0, 1}, {0, 0, 1}}, bottom{{0, 0, 0}, {0, 0, 1}};
    s.check("the plane through the top face, moved down infinitesimally: area 1", cg::section_area(cg::plane_section(C, top), top) - 1, 1e-15);
    s.require("the plane through the bottom face: nothing", cg::plane_section(C, bottom).empty());
    const double k = 1 / std::sqrt(3.0);
    const cg::Plane3 oblique{{0.5, 0.5, 0.5}, {k, k, k}};
    const auto H = cg::plane_section(C, oblique);
    s.require("the plane through the centre normal to (1, 1, 1): one closed loop", H.size() == 1 && H[0].closed);
    s.check("the regular hexagon of side √2/2, area 3√3/4", cg::section_area(H, oblique) - 3 * std::sqrt(3.0) / 4, 1e-14);
  }
  {
    const auto S = icosphere(4);
    const double V = solid(S).mass;
    s.check("Cavalieri: a sphere of 5120 triangles, its sections along z integrate to its volume, relative", cavalieri(S, {0, 0, 1}, 4000) / V - 1, 2e-6);
    const double k = 1 / std::sqrt(14.0);
    s.check("and along (1, 2, 3)/√14, relative", cavalieri(S, {k, 2 * k, 3 * k}, 4000) / V - 1, 2e-6);
  }
  {
    const auto T = torus(2, 0.6, 64, 32);
    const cg::Plane3 eq{{0, 0, 0}, {0, 0, 1}};
    const auto L = cg::plane_section(T, eq);
    int positive = 0, negative = 0;
    double total = 0;
    for (const auto& c : L) {
      const double a = cg::ring_area(cg::in_plane(c.points, {0, 0, 0}, {1, 0, 0}, {0, 1, 0}));
      (a > 0 ? positive : negative)++;
      total += a;
    }
    s.require("a torus cut through its hole: two closed loops, one counter-clockwise and one clockwise", L.size() == 2 && L[0].closed && L[1].closed && positive == 1 && negative == 1);
    s.check("their areas add to the annulus π((R + r)² − (R − r)²) within the polygon's error, relative", total / (std::numbers::pi * (2.6 * 2.6 - 1.4 * 1.4)) - 1, 1e-2);
    const double V = solid(T).mass;
    s.check("Cavalieri on the torus along z, relative", cavalieri(T, {0, 0, 1}, 4000) / V - 1, 2e-6);
    s.check("and along x, through the hole, relative", cavalieri(T, {1, 0, 0}, 4000) / V - 1, 2e-6);
  }
  {
    cg::TriangleSurface Q;
    Q.points = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
    Q.triangles = {{0, 1, 2}, {0, 2, 3}};
    const auto L = cg::plane_section(Q, {{0.5, 0, 0}, {1, 0, 0}});
    s.require("a square of two triangles cut across: one open chain of three points", L.size() == 1 && !L[0].closed && L[0].points.size() == 3);
  }
  {
    std::vector<cg::Point3> P;
    std::vector<std::size_t> label;
    for (int c = 0; c < 3; c++) {
      const auto C = cube(c == 0 ? 0.0 : c == 1 ? 1.5 : 5.0);
      for (const auto& p : C.points) P.push_back(p), label.push_back(std::size_t(c));
    }
    const auto a = cg::close_pairs(P, label, 0.6), b = cg::close_pairs(P, label, 0.4), c = cg::close_pairs(P, label, 3.6);
    s.require("three cubes, gaps 0.5 and 2.5: within 0.6 only the first two, at 0.5", a.size() == 1 && a[0].a == 0 && a[0].b == 1 && std::fabs(a[0].distance - 0.5) < 1e-15);
    s.require("within 0.4 none", b.empty());
    s.require("within 3.6 the first two and the last two, not the first and the last (4.0 apart)",
              c.size() == 2 && c[0].a == 0 && c[0].b == 1 && c[1].a == 1 && c[1].b == 2 && std::fabs(c[1].distance - 2.5) < 1e-15);

    calx::random::Stream z(19, 2, 0, 0);
    P.clear();
    label.clear();
    for (int i = 0; i < 2000; i++) {
      P.push_back({z.uniform(-5, 5), z.uniform(-5, 5), z.uniform(-1, 1)});
      label.push_back(std::size_t(z.uniform() * 7));
    }
    const double d = 0.35;
    std::map<std::pair<std::size_t, std::size_t>, double> brute;
    for (std::size_t i = 0; i < P.size(); i++)
      for (std::size_t j = i + 1; j < P.size(); j++) {
        if (label[i] == label[j]) continue;
        const double r = std::hypot(P[i][0] - P[j][0], P[i][1] - P[j][1], P[i][2] - P[j][2]);
        if (r > d) continue;
        const auto k = std::minmax(label[i], label[j]);
        auto [it, fresh] = brute.try_emplace({k.first, k.second}, r);
        if (!fresh) it->second = std::min(it->second, r);
      }
    const auto found = cg::close_pairs(P, label, d);
    double worst = found.size() == brute.size() ? 0 : 1;
    for (const auto& p : found) {
      const auto it = brute.find({p.a, p.b});
      worst = std::max(worst, it == brute.end() ? 1.0 : std::fabs(it->second - p.distance));
    }
    s.check("2000 random points of seven labels: the close pairs and their least distances, as every pair compared", worst, 1e-15);
  }
  return s.done();
}
