// The bases of fields as regions with a test of membership (decision 0012; docs/theory/fields.md §2).
//
// A field is defined on its base and nowhere else: outside it a query has no value, which is not the
// value zero (formalism §3.2). A domain is anything that answers, for a query, whether it is covered,
// and if not why: covers(q) → Cover. The few here name no body; a domain of the World Magnetic Model
// (heights above WGS 84, the model's lifetime) is the geophysics domain's, built from these or written
// beside them.
//
// Reasons are views of static text: a domain made with its own words keeps them for as long as the
// program runs.

#pragma once

#include <cmath>
#include <limits>
#include <string_view>
#include <tuple>

#include <physis/fields/field.hpp>

namespace physis::fields {

struct Cover {
  Quality quality = Quality::valid;
  std::string_view reason{};
  constexpr bool ok() const noexcept { return quality == Quality::valid; }
};

template <class D, class P>
concept DomainOf = requires(const D& d, const Query<P>& q) {
  { d.covers(q) } -> std::same_as<Cover>;
};

// Every point at every instant.
struct Everywhere {
  template <class P>
  constexpr Cover covers(const Query<P>&) const noexcept { return {}; }
};

// The points between two spheres about a centre, inner ≤ |x − c| ≤ outer.
struct Shell {
  Vec3 centre{};
  double inner = 0;
  double outer = std::numeric_limits<double>::infinity();
  std::string_view below = "below the inner sphere of the field's shell";
  std::string_view above = "above the outer sphere of the field's shell";

  Cover covers(const Query<Vec3>& q) const noexcept {
    const double r = std::hypot(q.x[0] - centre[0], q.x[1] - centre[1], q.x[2] - centre[2]);
    if (!(r >= inner)) return {Quality::outside, below};
    if (!(r <= outer)) return {Quality::outside, above};
    return {};
  }
};

// The points of an axis-aligned box of the field's frame, lower ≤ x ≤ upper in each coordinate.
struct Box {
  Vec3 lower{};
  Vec3 upper{};
  std::string_view beyond = "outside the field's box";

  Cover covers(const Query<Vec3>& q) const noexcept {
    for (int i = 0; i < 3; i++)
      if (!(q.x[i] >= lower[i] && q.x[i] <= upper[i])) return {Quality::outside, beyond};
    return {};
  }
};

// An interval of a line, lower ≤ x ≤ upper.
struct Interval {
  double lower = -std::numeric_limits<double>::infinity();
  double upper = std::numeric_limits<double>::infinity();
  std::string_view below = "below the field's interval";
  std::string_view above = "above the field's interval";

  Cover covers(const Query<double>& q) const noexcept {
    if (!(q.x >= lower)) return {Quality::outside, below};
    if (!(q.x <= upper)) return {Quality::outside, above};
    return {};
  }
};

// The instants of a span of time, from ≤ t ≤ to, at every point.
struct During {
  double from = -std::numeric_limits<double>::infinity();
  double to = std::numeric_limits<double>::infinity();
  std::string_view before = "before the field's span of time";
  std::string_view after = "after the field's span of time";

  template <class P>
  constexpr Cover covers(const Query<P>& q) const noexcept {
    if (!(q.t >= from)) return {Quality::outside, before};
    if (!(q.t <= to)) return {Quality::outside, after};
    return {};
  }
};

// The intersection of domains: the first that refuses gives the reason.
template <class... D>
struct AllOf {
  std::tuple<D...> parts;

  template <class P>
  constexpr Cover covers(const Query<P>& q) const {
    Cover c;
    std::apply([&](const auto&... d) { ((c.ok() ? (void)(c = d.covers(q)) : (void)0), ...); }, parts);
    return c;
  }
};
template <class... D>
AllOf<D...> all_of(D... d) { return {std::tuple<D...>(std::move(d)...)}; }

}  // namespace physis::fields
