// Checks of physis/symmetry/point_group.hpp. For each of 26 point groups, from C1 to Ih, the group
// is generated from its generators, and the orbit of two generic points with different labels
// under it, a set whose symmetry group is exactly that group, is handed to point_group(): it must
// return the group's order and Schoenflies symbol, a multiplication table that satisfies the group
// axioms, the number of conjugacy classes the group has (Cotton 1990, Appendix IIA: the character
// tables), the polar groups as those that fix a vector (C1, Cs, Cn, Cnv), the groups that fix an
// axial vector (C1, Ci, Cs, Cn, Cnh, S2n), and the chiral groups as the proper ones (C1, Cn, Dn,
// T, O, I). Then congruences: two configurations of a tetrahedron with four distinct labels,
// mirror images, are congruent only by improper isometries, one each way.

#include <physis/symmetry/point_group.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace ps = physis::symmetry;
namespace ma = morphe::algebra;

namespace {

const double pi = std::acos(-1.0);
const ps::Vec3 X{1, 0, 0}, Y{0, 1, 0}, Z{0, 0, 1};

ps::Mat3 inv() { return {{{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}}}; }
ps::Mat3 C(int n, ps::Vec3 axis = Z) { return ps::rotation(axis, 2 * pi / n); }
ps::Mat3 S(int n) { return ps::pg_detail::mul(ps::reflection(Z), C(n)); }  // rotoreflection about z
ps::Vec3 unit(ps::Vec3 v) { return ps::pg_detail::scale(1 / ps::pg_detail::norm(v), v); }

struct Case {
  std::string name;
  std::vector<ps::Mat3> generators;
  std::size_t order, classes;
  bool polar, axial, chiral;
};

std::vector<Case> cases() {
  const double phi = (1 + std::sqrt(5.0)) / 2;
  const ps::Vec3 diag = unit({1, 1, 1}), five = unit({0, 1, phi});
  const auto sigma_v = ps::reflection(X), sigma_h = ps::reflection(Z);
  const auto C2x = C(2, X);
  std::vector<Case> c = {
      {"C1", {}, 1, 1, true, true, true},
      {"Cs", {sigma_h}, 2, 2, true, true, false},
      {"Ci", {inv()}, 2, 2, false, true, false},
      {"C2", {C(2)}, 2, 2, true, true, true},
      {"C3", {C(3)}, 3, 3, true, true, true},
      {"C4", {C(4)}, 4, 4, true, true, true},
      {"C6", {C(6)}, 6, 6, true, true, true},
      {"C2v", {C(2), sigma_v}, 4, 4, true, false, false},
      {"C3v", {C(3), sigma_v}, 6, 3, true, false, false},
      {"C6v", {C(6), sigma_v}, 12, 6, true, false, false},
      {"C2h", {C(2), sigma_h}, 4, 4, false, true, false},
      {"C3h", {C(3), sigma_h}, 6, 6, false, true, false},
      {"S4", {S(4)}, 4, 4, false, true, false},
      {"S6", {S(6)}, 6, 6, false, true, false},
      {"D2", {C(2), C2x}, 4, 4, false, false, true},
      {"D3", {C(3), C2x}, 6, 3, false, false, true},
      {"D2h", {C(2), C2x, sigma_h}, 8, 8, false, false, false},
      {"D3h", {C(3), C2x, sigma_h}, 12, 6, false, false, false},
      {"D6h", {C(6), C2x, sigma_h}, 24, 12, false, false, false},
      {"D2d", {C(2), C2x, S(4)}, 8, 5, false, false, false},
      {"D3d", {C(3), C2x, S(6)}, 12, 6, false, false, false},
      {"T", {C(2), C(3, diag)}, 12, 4, false, false, true},
      {"Td", {C(2), C(3, diag), ps::reflection(unit({1, -1, 0}))}, 24, 5, false, false, false},
      {"Th", {C(2), C(3, diag), inv()}, 24, 8, false, false, false},
      {"O", {C(4), C(3, diag)}, 24, 5, false, false, true},
      {"Oh", {C(4), C(3, diag), inv()}, 48, 10, false, false, false},
      {"I", {C(5, five), C(3, diag)}, 60, 5, false, false, true},
      {"Ih", {C(5, five), C(3, diag), inv()}, 120, 10, false, false, false},
  };
  return c;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  // Four generic points with distinct labels: fixed by no element, and not coplanar, so that even
  // the orbit under the trivial group spans space (two or three points would lie in a line or a
  // plane, whose own reflections would enlarge the group).
  const std::vector<ps::Vec3> generic{{0.731, 0.247, 0.418}, {-0.193, 0.862, -0.351}, {0.402, -0.517, 0.611}, {-0.655, -0.141, -0.298}};
  for (const auto& c : cases()) {
    const auto elements = ps::generate(c.generators);
    ps::LabelledPoints A;
    for (const auto& g : elements)
      for (std::size_t k = 0; k < generic.size(); k++) A.add(std::string(1, static_cast<char>('a' + k)), ps::pg_detail::apply(g, generic[k]));
    // Shift the whole set, so that the centroid is found rather than assumed.
    for (auto& x : A.x) x = {x[0] + 3.1, x[1] - 1.7, x[2] + 0.4};
    ps::PointGroup G;
    try {
      G = ps::point_group(A);
    } catch (const std::exception& e) {
      s.require(c.name + ": point group", false, e.what());
      continue;
    }
    const std::size_t classes = ma::conjugacy_classes(G.group).size();
    const bool ok = elements.size() == c.order && G.order() == c.order && G.schoenflies == c.name && classes == c.classes;
    s.require(c.name + ": order, symbol and classes", ok,
              G.schoenflies + ", |G| = " + std::to_string(G.order()) + ", " + std::to_string(classes) + " classes");
    s.require(c.name + ": polar, axial, chiral",
              (G.fixed_vectors(false) > 0) == c.polar && (G.fixed_vectors(true) > 0) == c.axial && G.chiral() == c.chiral,
              "fixed vectors " + std::to_string(G.fixed_vectors(false)) + ", axial " + std::to_string(G.fixed_vectors(true)));
  }

  // The orbit of one point of a planar set: the reflection in the plane is found though it moves
  // no point (a regular hexagon of equal labels: D6h, 24 elements, 6 points).
  {
    ps::LabelledPoints hexagon;
    for (int k = 0; k < 6; k++) hexagon.add("x", {std::cos(k * pi / 3), std::sin(k * pi / 3), 0});
    const auto G = ps::point_group(hexagon);
    s.require("regular hexagon: D6h", G.schoenflies == "D6h" && G.order() == 24, G.schoenflies);
    ps::LabelledPoints line;
    line.add("x", {0, 0, -1}), line.add("y", {0, 0, 0}), line.add("x", {0, 0, 1});
    ps::LabelledPoints polar_line;
    polar_line.add("x", {0, 0, -1}), polar_line.add("y", {0, 0, 0}), polar_line.add("z", {0, 0, 1.2});
    s.require("points on a line: D∞h and C∞v", ps::point_group(line).schoenflies == "D∞h" && ps::point_group(polar_line).schoenflies == "C∞v");
  }

  // Congruences of mirror-image tetrahedral configurations with four distinct labels.
  {
    const double h = 1 / std::sqrt(3.0);
    ps::LabelledPoints plus, minus;
    const std::vector<std::string> labels{"w", "x", "y", "z"};
    const std::vector<ps::Vec3> v{{h, h, h}, {h, -h, -h}, {-h, h, -h}, {-h, -h, h}};
    plus.add("c", {0, 0, 0}), minus.add("c", {0, 0, 0});
    for (int i = 0; i < 4; i++) {
      plus.add(labels[i], v[i]);
      minus.add(labels[i], {-v[i][0], v[i][1], v[i][2]});  // reflected in the plane x = 0
    }
    const auto self = ps::point_group(plus);
    const auto between = ps::congruences(plus, minus);
    const auto back = ps::congruences(minus, plus);
    s.require("four distinct labels on a tetrahedron: C1, chiral", self.schoenflies == "C1" && self.chiral());
    s.require("its mirror image is congruent to it by one isometry, improper", between.size() == 1 && between[0].det() == -1 && back.size() == 1 && back[0].det() == -1);
    // A proper rotation of the configuration is congruent by exactly one proper isometry.
    ps::LabelledPoints turned = plus;
    const auto R = ps::rotation(unit({0.3, -0.5, 0.8}), 1.234);
    for (auto& x : turned.x) x = ps::pg_detail::apply(R, x);
    const auto rot = ps::congruences(plus, turned);
    s.require("a rotated copy is congruent by one proper isometry", rot.size() == 1 && rot[0].det() == 1 && ps::pg_detail::max_abs_diff(rot[0].Q, R) < 1e-9);
  }
  return s.done();
}
