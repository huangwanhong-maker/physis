// Checks of the faces of a plane graph (planar.hpp).
//
// - A 3 × 2 grid of unit squares: six bounded faces of area 1, simple, traversed counter-clockwise, and one
//   outer face of area −6; V − E + F = 12 − 17 + 7 = 2.
// - A dead end inside a square: the square's face visits the corner twice (not simple), its area still 1.
// - Edges as polylines: a square whose sides bow outwards by a polyline point each has the area of the
//   polygon through all the points; a digon of two bowed edges is a simple face of two vertices.
// - Two components: V − E + F = 2 × 2 with F the face orbits, an outer walk each.
// - Polylines: length, the point at an arc length, a cut into two whose lengths add, the nearest point, a ring's
//   area, crossing segments (and touching ones that do not cross), a point inside a ring.

#include <cmath>
#include <vector>

#include <calx/verification/check.hpp>
#include <physis/computational_geometry/planar.hpp>

using calx::verification::Suite;
namespace cg = physis::computational_geometry;

int main() {
  Suite s;
  {
    std::vector<cg::Point2> P;
    for (int y = 0; y <= 2; y++)
      for (int x = 0; x <= 3; x++) P.push_back({double(x), double(y)});
    std::vector<cg::PlaneEdge> E;
    for (int y = 0; y <= 2; y++)
      for (int x = 0; x < 3; x++) E.push_back({std::size_t(y * 4 + x), std::size_t(y * 4 + x + 1), {}});
    for (int y = 0; y < 2; y++)
      for (int x = 0; x <= 3; x++) E.push_back({std::size_t(y * 4 + x), std::size_t((y + 1) * 4 + x), {}});
    const auto F = cg::plane_faces(P, E);
    std::size_t bounded = 0, simple = 0;
    double worst = 0, outer = 0;
    for (const auto& f : F) {
      if (f.area > 0) bounded++, simple += f.simple, worst = std::max(worst, std::fabs(f.area - 1));
      else outer = f.area;
    }
    s.require("a 3 × 2 grid: six bounded faces, all simple, and one outer", bounded == 6 && simple == 6 && F.size() == 7);
    s.check("each of area 1", worst, 1e-15);
    s.check("the outer face of area −6", outer + 6, 1e-15);
    s.require("Euler: V − E + F = 2", static_cast<long>(P.size()) - static_cast<long>(E.size()) + static_cast<long>(F.size()) == 2);
  }
  {
    std::vector<cg::Point2> P{{0, 0}, {1, 0}, {1, 1}, {0, 1}, {0.5, 0.5}};
    std::vector<cg::PlaneEdge> E{{0, 1, {}}, {1, 2, {}}, {2, 3, {}}, {3, 0, {}}, {0, 4, {}}};
    const auto F = cg::plane_faces(P, E);
    const cg::PlaneFace* inner = nullptr;
    for (const auto& f : F)
      if (f.area > 0) inner = &f;
    s.require("a dead end inside a square: the face visits a corner twice", inner && !inner->simple && inner->vertices.size() == 6);
    s.check("and its area is still 1", inner ? inner->area - 1 : 1, 1e-15);
  }
  {
    std::vector<cg::Point2> P{{0, 0}, {2, 0}, {2, 2}, {0, 2}};
    std::vector<cg::PlaneEdge> E{{0, 1, {{0, 0}, {1, -0.5}, {2, 0}}}, {1, 2, {{2, 0}, {2.5, 1}, {2, 2}}}, {2, 3, {{2, 2}, {1, 2.5}, {0, 2}}}, {3, 0, {{0, 2}, {-0.5, 1}, {0, 0}}}};
    const auto F = cg::plane_faces(P, E);
    double area = 0;
    for (const auto& f : F)
      if (f.area > 0) area = f.area;
    s.check("a square with bowed sides: the area of the polygon through all its points, 4 + 4 × 0.5", area - 6, 1e-14);
    std::vector<cg::PlaneEdge> D{{0, 1, {{0, 0}, {1, -0.5}, {2, 0}}}, {1, 0, {{2, 0}, {1, 0.5}, {0, 0}}}};
    const auto G = cg::plane_faces(P, D);
    bool digon = false;
    for (const auto& f : G)
      if (f.area > 0) digon = f.simple && f.vertices.size() == 2 && std::fabs(f.area - 1) < 1e-14;
    s.require("a digon of two bowed edges is a simple face of two vertices, of area 1", digon);
  }
  {
    std::vector<cg::Point2> P{{0, 0}, {1, 0}, {0, 1}, {5, 5}, {6, 5}, {5, 6}};
    std::vector<cg::PlaneEdge> E{{0, 1, {}}, {1, 2, {}}, {2, 0, {}}, {3, 4, {}}, {4, 5, {}}, {5, 3, {}}};
    const auto F = cg::plane_faces(P, E);
    std::size_t outers = 0;
    for (const auto& f : F) outers += f.area < 0;
    s.require("two components: V − E + F = 2 × 2 with F the face orbits, an outer walk each", 6 - 6 + static_cast<long>(F.size()) == 4 && outers == 2);
  }
  {
    const std::vector<cg::Point2> l{{0, 0}, {3, 0}, {3, 4}};
    s.check("a polyline's length", cg::polyline_length(l) - 7, 1e-15);
    const auto p = cg::point_at(l, 5);
    s.check("the point at arc length 5", std::hypot(p[0] - 3, p[1] - 2), 1e-15);
    const auto [a, b] = cg::split_polyline(l, 5);
    s.check("cut at 5, its pieces' lengths add up", cg::polyline_length(a) - 5 + cg::polyline_length(b) - 2, 1e-15);
    s.require("and meet at the cut", a.back() == b.front() && a.size() == 3 && b.size() == 2);
    const auto [at, dist] = cg::nearest_on_polyline(l, {4, 1});
    s.check("the nearest point to (4, 1) is at arc length 4, distance 1", std::fabs(at - 4) + std::fabs(dist - 1), 1e-15);
    s.check("a ring's area", cg::ring_area({{0, 0}, {2, 0}, {2, 3}, {0, 3}}) - 6, 1e-15);
    s.require("crossing segments cross; touching ones do not", cg::segments_cross({0, 0}, {2, 2}, {0, 2}, {2, 0}) && !cg::segments_cross({0, 0}, {1, 1}, {1, 1}, {2, 0}));
    s.require("a point inside a ring, and one outside", cg::inside_ring({{0, 0}, {2, 0}, {2, 3}, {0, 3}}, {1, 1}) && !cg::inside_ring({{0, 0}, {2, 0}, {2, 3}, {0, 3}}, {3, 1}));
  }
  return s.done();
}
