// Integral curves of vector fields: field lines and streamlines (decision 0012; docs/theory/fields.md §6).
//
// A field line of v through a seed is the curve whose tangent is v everywhere, parametrised by arc
// length: dx/ds = v/|v|. A streamline (the path of a particle in a steady flow) is dx/dτ = v. Both are
// integrated by the classical fourth-order Runge–Kutta method (Hairer, Nørsett and Wanner 1993,
// §II.1), with a step either fixed or a fraction of the distance from a centre (so that lines about a
// body take steps in proportion to their scale).
//
// A curve ends, in each direction, for one of four reasons, which it reports: it reaches the edge of
// the field's base (a stage of the step has no value: the curve's last step is halved until it fits,
// up to `refine` times, so that the curve ends within 2^−refine of a step of the edge); it reaches a
// null of the field (|v| ≤ zero, where the direction is undefined); it has taken max_steps; it has
// reached max_length.

#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string_view>
#include <type_traits>
#include <vector>

#include <physis/fields/field.hpp>

namespace physis::fields {

struct CurveOptions {
  enum class Direction : std::uint8_t { forward, backward, both };
  Direction direction = Direction::both;
  bool along_unit = true;       // a field line (dx/ds = v/|v|); false: a streamline (dx/dτ = v)
  double step = 0;              // a fixed step (m, or s for a streamline), or 0 for step_fraction
  double step_fraction = 0.02;  // the step as a fraction of the distance from `about`
  Vec3 about{};
  std::size_t max_steps = 4000;
  double max_length = std::numeric_limits<double>::infinity();
  double zero = 0;  // |v| at or below it is a null of the field
  int refine = 12;  // halvings of a step that leaves the base
};

enum class CurveStop : std::uint8_t { edge, null, steps, length };

constexpr std::string_view name(CurveStop s) {
  switch (s) {
    case CurveStop::edge: return "edge";
    case CurveStop::null: return "null";
    case CurveStop::steps: return "steps";
    case CurveStop::length: return "length";
  }
  return "?";
}

struct CurveEnd {
  CurveStop stop = CurveStop::steps;
  Quality quality = Quality::valid;  // at the edge: why the field has no value beyond it
  std::string_view reason{};
};

// The points from the backward end to the forward end, the seed among them, with the arc length from
// the backward end.
struct Curve {
  std::vector<Vec3> points;
  std::vector<double> arc;
  CurveEnd backward, forward;
  std::size_t seed_index = 0;
  bool seeded = true;  // false when the seed itself has no value: the curve is empty
  Quality seed_quality = Quality::valid;
  std::string_view seed_reason{};
};

namespace curves_detail {

struct Stage {
  bool ok = false;
  bool null = false;
  Vec3 k{};
  Quality quality = Quality::valid;
  std::string_view reason{};
};

template <class F>
Stage stage(const F& f, double t, const Vec3& x, double sign, const CurveOptions& o) {
  const auto e = f.at({t, x});
  Stage s;
  if (!e.ok()) {
    s.quality = e.quality;
    s.reason = e.reason;
    return s;
  }
  const double n = geometry::norm(e.value);
  if (!(n > o.zero)) {
    s.null = true;
    return s;
  }
  const double scale = o.along_unit ? sign / n : sign;
  s.ok = true;
  s.k = {scale * e.value[0], scale * e.value[1], scale * e.value[2]};
  return s;
}

// One direction from the seed: the points after the seed, and how it ended.
template <class F>
void half(const F& f, const Vec3& seed, double t, double sign, const CurveOptions& o, std::vector<Vec3>& out, CurveEnd& end) {
  using namespace geometry::operators;
  Vec3 x = seed;
  double length = 0;
  for (std::size_t n = 0; n < o.max_steps; n++) {
    const Stage k1 = stage(f, t, x, sign, o);
    if (k1.null) {
      end = {CurveStop::null, Quality::valid, {}};
      return;
    }
    if (!k1.ok) {  // cannot happen after the seed: every point kept had a value
      end = {CurveStop::edge, k1.quality, k1.reason};
      return;
    }
    const double h0 = o.step > 0 ? o.step : o.step_fraction * geometry::norm(x - o.about);
    double h = h0;
    for (int halving = 0;; halving++) {
      const Stage k2 = stage(f, t, x + (h / 2) * k1.k, sign, o);
      const Stage k3 = k2.ok ? stage(f, t, x + (h / 2) * k2.k, sign, o) : k2;
      const Stage k4 = k3.ok ? stage(f, t, x + h * k3.k, sign, o) : k3;
      const Vec3 next = x + (h / 6) * (k1.k + 2.0 * k2.k + 2.0 * k3.k + k4.k);
      const Stage there = k4.ok ? stage(f, t, next, sign, o) : k4;
      if (there.ok || there.null) {
        length += geometry::norm(next - x);
        out.push_back(next);
        x = next;
        break;
      }
      if (halving >= o.refine) {
        end = {CurveStop::edge, there.quality, there.reason};
        return;
      }
      h /= 2;
    }
    if (length >= o.max_length) {
      end = {CurveStop::length, Quality::valid, {}};
      return;
    }
  }
  end = {CurveStop::steps, Quality::valid, {}};
}

}  // namespace curves_detail

template <PointEvaluable F>
  requires std::is_same_v<typename F::value_type, Vec3> && std::is_same_v<typename F::point_type, Vec3>
Curve integral_curve(const F& f, const Vec3& seed, double t, const CurveOptions& o = {}) {
  Curve c;
  const auto at_seed = f.at({t, seed});
  if (!at_seed.ok()) {
    c.seeded = false;
    c.seed_quality = at_seed.quality;
    c.seed_reason = at_seed.reason;
    return c;
  }
  std::vector<Vec3> back, ahead;
  if (o.direction != CurveOptions::Direction::forward) curves_detail::half(f, seed, t, -1.0, o, back, c.backward);
  if (o.direction != CurveOptions::Direction::backward) curves_detail::half(f, seed, t, 1.0, o, ahead, c.forward);
  c.points.reserve(back.size() + 1 + ahead.size());
  for (auto it = back.rbegin(); it != back.rend(); ++it) c.points.push_back(*it);
  c.seed_index = c.points.size();
  c.points.push_back(seed);
  for (const Vec3& x : ahead) c.points.push_back(x);
  c.arc.assign(c.points.size(), 0.0);
  for (std::size_t i = 1; i < c.points.size(); i++) {
    const Vec3& a = c.points[i - 1];
    const Vec3& b = c.points[i];
    c.arc[i] = c.arc[i - 1] + std::hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
  }
  return c;
}

}  // namespace physis::fields
