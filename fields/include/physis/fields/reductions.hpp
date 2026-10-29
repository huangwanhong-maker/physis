// Reductions of a fibre: what is shown or measured of a field's value (decision 0012;
// docs/theory/fields.md §4; formalism Definition 3.10).
//
// A reduction ρ: F → F′ maps a field's value at a point to a value of another fibre: the magnitude of
// a vector, its component along an axis, its radial and horizontal parts about a centre, its
// components on a local basis. Each says whether it is *invariant*: a magnitude, or a radial
// component about a given centre, is the same in every frame; a component on a frame's axes names that
// frame, and a drawing of it must name it too (formalism §3.6, Remark on invariance).
//
// A reduction is a value with operator()(query, value), the fibre it produces from the one it is given,
// a name and the flag `invariant`.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <string_view>

#include <physis/fields/field.hpp>

namespace physis::fields {

template <class R>
concept Reduction = requires(const R& r, const Query<typename R::point_type>& q, const typename R::input& v, const Fibre& f) {
  typename R::input;
  typename R::output;
  typename R::point_type;
  { r(q, v) } -> std::convertible_to<typename R::output>;
  { r.fibre(f) } -> std::same_as<Fibre>;
  { r.name() } -> std::convertible_to<std::string_view>;
  { R::invariant } -> std::convertible_to<bool>;
};

namespace reduction_detail {
inline Fibre scalar_of(const Fibre& f) {
  Fibre s;
  s.kind = FibreKind::scalar;
  s.dimension = f.dimension;
  return s;
}
}  // namespace reduction_detail

// |v|.
struct Magnitude {
  using input = Vec3;
  using output = double;
  using point_type = Vec3;
  static constexpr bool invariant = true;
  double operator()(const Query<Vec3>&, const Vec3& v) const { return geometry::norm(v); }
  Fibre fibre(const Fibre& f) const { return reduction_detail::scalar_of(f); }
  std::string_view name() const { return "magnitude"; }
};

// v's component along one of the frame's axes.
struct Component {
  using input = Vec3;
  using output = double;
  using point_type = Vec3;
  static constexpr bool invariant = false;
  int axis = 0;
  double operator()(const Query<Vec3>&, const Vec3& v) const { return v[static_cast<std::size_t>(axis)]; }
  Fibre fibre(const Fibre& f) const {
    Fibre s = reduction_detail::scalar_of(f);
    s.basis = f.basis;
    return s;
  }
  std::string_view name() const {
    constexpr std::array<std::string_view, 3> n{"x component", "y component", "z component"};
    return n[static_cast<std::size_t>(axis)];
  }
};

// v · d for a unit direction d on the frame's axes.
struct Along {
  using input = Vec3;
  using output = double;
  using point_type = Vec3;
  static constexpr bool invariant = false;
  Vec3 direction{0, 0, 1};
  std::string_view label = "along a direction";
  double operator()(const Query<Vec3>&, const Vec3& v) const { return geometry::dot(v, direction); }
  Fibre fibre(const Fibre& f) const { return reduction_detail::scalar_of(f); }
  std::string_view name() const { return label; }
};

// v · r̂, r̂ = (x − c)/|x − c|: outward positive.
struct Radial {
  using input = Vec3;
  using output = double;
  using point_type = Vec3;
  static constexpr bool invariant = true;
  Vec3 centre{};
  double operator()(const Query<Vec3>& q, const Vec3& v) const {
    const Vec3 r{q.x[0] - centre[0], q.x[1] - centre[1], q.x[2] - centre[2]};
    const double n = geometry::norm(r);
    return n > 0 ? geometry::dot(v, r) / n : 0.0;
  }
  Fibre fibre(const Fibre& f) const { return reduction_detail::scalar_of(f); }
  std::string_view name() const { return "radial"; }
};

// |v − (v · r̂) r̂|: the part across the radius.
struct Horizontal {
  using input = Vec3;
  using output = double;
  using point_type = Vec3;
  static constexpr bool invariant = true;
  Vec3 centre{};
  double operator()(const Query<Vec3>& q, const Vec3& v) const {
    const Vec3 r{q.x[0] - centre[0], q.x[1] - centre[1], q.x[2] - centre[2]};
    const double n = geometry::norm(r);
    if (!(n > 0)) return geometry::norm(v);
    const double vr = geometry::dot(v, r) / n;
    return std::sqrt(std::max(0.0, geometry::dot(v, v) - vr * vr));
  }
  Fibre fibre(const Fibre& f) const { return reduction_detail::scalar_of(f); }
  std::string_view name() const { return "horizontal"; }
};

// v on the local east, north and up of a sphere about a centre, at the point (at a pole, the meridian
// of atan2(y, x)).
struct LocalOfSphere {
  using input = Vec3;
  using output = Vec3;
  using point_type = Vec3;
  static constexpr bool invariant = false;
  Vec3 centre{};
  Vec3 operator()(const Query<Vec3>& q, const Vec3& v) const {
    const Vec3 r{q.x[0] - centre[0], q.x[1] - centre[1], q.x[2] - centre[2]};
    const double psi = std::atan2(r[2], std::hypot(r[0], r[1])), lambda = std::atan2(r[1], r[0]);
    const double sp = std::sin(psi), cp = std::cos(psi), sl = std::sin(lambda), cl = std::cos(lambda);
    const Vec3 up{cp * cl, cp * sl, sp}, east{-sl, cl, 0}, north{-sp * cl, -sp * sl, cp};
    return {geometry::dot(v, east), geometry::dot(v, north), geometry::dot(v, up)};
  }
  Fibre fibre(const Fibre& f) const {
    Fibre g = f;
    g.basis = "east, north, up of the sphere";
    return g;
  }
  std::string_view name() const { return "local east, north, up"; }
};

// R v: v on the axes whose directions are R's rows (a cube's axes, a vehicle's body axes).
struct OnAxes {
  using input = Vec3;
  using output = Vec3;
  using point_type = Vec3;
  static constexpr bool invariant = false;
  Mat3 rows = geometry::identity3();
  std::string_view basis = "the given axes";
  Vec3 operator()(const Query<Vec3>&, const Vec3& v) const {
    using namespace geometry::operators;
    return rows * v;
  }
  Fibre fibre(const Fibre& f) const {
    Fibre g = f;
    g.basis = std::string(basis);
    return g;
  }
  std::string_view name() const { return basis; }
};

}  // namespace physis::fields
