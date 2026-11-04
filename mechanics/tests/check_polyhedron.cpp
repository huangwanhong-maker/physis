// Checks of the solid bounded by a closed triangle surface (mechanics/rigid_body.hpp, solid_polyhedron;
// decision 0019).
//
// - A box as twelve triangles is solid_box: its mass, centre and inertia to rounding; the same box far from
//   the origin (10⁴ away), its diagonals cut the other way, is still the box.
// - A box turned by a rotation R has the inertia R J Rᵀ.
// - A regular tetrahedron of edge a has an isotropic inertia, m a²/20 about every axis through its centre,
//   and the volume a³/(6√2).
// - A surface oriented inward encloses no positive volume, and is refused.

#include <array>
#include <cmath>
#include <numbers>
#include <vector>

#include <calx/verification/check.hpp>
#include <physis/mechanics/rigid_body.hpp>

using calx::verification::Suite;
namespace pm = physis::mechanics;
namespace pg = physis::geometry;
using pg::Mat3;
using pg::Vec3;
using Tri = std::array<std::size_t, 3>;

namespace {

// The box [lo, hi] as twelve triangles oriented outward; `flip` cuts each face along its other diagonal.
void box(const Vec3& lo, const Vec3& hi, bool flip, std::vector<Vec3>& P, std::vector<Tri>& T) {
  P.clear();
  T.clear();
  for (int k = 0; k < 8; k++) P.push_back({k & 1 ? hi[0] : lo[0], k & 2 ? hi[1] : lo[1], k & 4 ? hi[2] : lo[2]});
  // Faces as quads (a, b, c, d) counter-clockwise seen from outside.
  const int quads[6][4] = {{0, 2, 3, 1}, {4, 5, 7, 6}, {0, 1, 5, 4}, {2, 6, 7, 3}, {0, 4, 6, 2}, {1, 3, 7, 5}};
  for (const auto& q : quads) {
    if (!flip) {
      T.push_back({std::size_t(q[0]), std::size_t(q[1]), std::size_t(q[2])});
      T.push_back({std::size_t(q[0]), std::size_t(q[2]), std::size_t(q[3])});
    } else {
      T.push_back({std::size_t(q[0]), std::size_t(q[1]), std::size_t(q[3])});
      T.push_back({std::size_t(q[1]), std::size_t(q[2]), std::size_t(q[3])});
    }
  }
}

double largest(const Mat3& a) {
  double m = 0;
  for (const auto& r : a)
    for (double x : r) m = std::max(m, std::fabs(x));
  return m;
}
double difference(const Mat3& a, const Mat3& b) {
  double m = 0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) m = std::max(m, std::fabs(a[i][j] - b[i][j]));
  return m;
}

}  // namespace

int main() {
  Suite s;
  {
    std::vector<Vec3> P;
    std::vector<Tri> T;
    const Vec3 lo{1, -1, 2}, hi{3, 0.5, 2.25};
    const double rho = 7.5;
    box(lo, hi, false, P, T);
    const auto A = pm::solid_polyhedron(P, T, rho);
    const auto B = pm::solid_box(lo, hi, rho);
    s.check("a box of twelve triangles: its mass is the box's", A.mass / B.mass - 1, 1e-14);
    s.check("its centre", std::hypot(A.centre[0] - B.centre[0], A.centre[1] - B.centre[1], A.centre[2] - B.centre[2]), 1e-14);
    s.check("its inertia, relative", difference(A.inertia, B.inertia) / largest(B.inertia), 1e-14);

    const Vec3 far{1e4, -2e4, 3e4};
    box({lo[0] + far[0], lo[1] + far[1], lo[2] + far[2]}, {hi[0] + far[0], hi[1] + far[1], hi[2] + far[2]}, true, P, T);
    const auto C = pm::solid_polyhedron(P, T, rho);
    s.check("the box 10⁴ from the origin, cut along the other diagonals: its mass", C.mass / B.mass - 1, 1e-12);
    s.check("its inertia, relative", difference(C.inertia, B.inertia) / largest(B.inertia), 1e-10);

    const auto R = pg::Rotation::about(Vec3{1 / std::sqrt(5.25), 2 / std::sqrt(5.25), -0.5 / std::sqrt(5.25)}, 0.7).matrix();
    box(lo, hi, false, P, T);
    for (auto& p : P) p = {R[0][0] * p[0] + R[0][1] * p[1] + R[0][2] * p[2], R[1][0] * p[0] + R[1][1] * p[1] + R[1][2] * p[2], R[2][0] * p[0] + R[2][1] * p[1] + R[2][2] * p[2]};
    const auto D = pm::solid_polyhedron(P, T, rho);
    Mat3 expected{};
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        for (int k = 0; k < 3; k++)
          for (int l = 0; l < 3; l++) expected[i][j] += R[i][k] * B.inertia[k][l] * R[j][l];
    s.check("a box turned by R: its inertia is R J Rᵀ, relative", difference(D.inertia, expected) / largest(expected), 1e-13);
  }
  {
    // A regular tetrahedron of edge 2√2, its faces oriented outward from its centre.
    std::vector<Vec3> P{{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}};
    std::vector<Tri> T{{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
    for (auto& t : T) {
      const Vec3 a = P[t[0]], b = P[t[1]], c = P[t[2]];
      const Vec3 n = pg::cross({b[0] - a[0], b[1] - a[1], b[2] - a[2]}, {c[0] - a[0], c[1] - a[1], c[2] - a[2]});
      if (pg::dot(n, a) < 0) std::swap(t[1], t[2]);
    }
    const double a = 2 * std::numbers::sqrt2, rho = 3;
    const auto M = pm::solid_polyhedron(P, T, rho);
    const double V = a * a * a / (6 * std::numbers::sqrt2);
    s.check("a regular tetrahedron: its volume a³/(6√2)", M.mass / (rho * V) - 1, 1e-15);
    const double I = M.mass * a * a / 20;
    Mat3 iso{};
    for (int i = 0; i < 3; i++) iso[i][i] = I;
    s.check("its inertia m a²/20 about every axis through its centre, relative", difference(M.inertia, iso) / I, 1e-14);

    std::vector<Tri> inward = T;
    for (auto& t : inward) std::swap(t[1], t[2]);
    bool refused = false;
    try {
      (void)pm::solid_polyhedron(P, inward, rho);
    } catch (const std::invalid_argument&) {
      refused = true;
    }
    s.require("the same surface oriented inward is refused", refused);
  }
  return s.done();
}
