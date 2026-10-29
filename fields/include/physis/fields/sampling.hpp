// Samplings: the finite sets of points on which a field is read (decision 0012; docs/theory/fields.md §5).
//
// A sampling is the restriction r of a visual expression (formalism Definition 3.10) made finite, and
// the support of every map, slice and cube: an ordered set of points of the base with a layout.
//
//   SphereGrid   cells of equal angle on a sphere about a centre, nlat rows from the north and nlon
//                columns from 180° W, the field read at the cells' centres: row i at latitude
//                90° − (i + ½)·180°/nlat, column j at longitude −180° + (j + ½)·360°/nlon, index i·nlon + j
//                (an image's order, north up, west left);
//   PlaneGrid    nu × nv cells of a parallelogram, origin + s u + r v for s, r in [0, 1], read at the
//                cells' centres, index i + nu·j;
//   VolumeGrid   n₀ × n₁ × n₂ cells of a box about a centre along three orthonormal axes, read at the
//                cells' centres, index i + n₀ (j + n₁ k), the first axis fastest (Fortran's order);
//   PointList    points in a given order.
//
// sample(f, g, t) reads f at every point of g at instant t and counts the answers by quality and by
// reason. A representation that reads a sampling faster than point by point (a spectral expansion on a
// sphere: the Legendre functions once per row) provides f.sample(g, t), which sample() then uses.

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <span>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

#include <physis/fields/field.hpp>

namespace physis::fields {

using namespace physis::geometry::operators;

struct SphereGrid {
  Vec3 centre{};
  double radius = 0;
  int nlat = 0, nlon = 0;

  std::size_t size() const { return static_cast<std::size_t>(nlat) * static_cast<std::size_t>(nlon); }
  double latitude(int i) const { return std::numbers::pi / 2 - (i + 0.5) * std::numbers::pi / nlat; }
  double longitude(int j) const { return -std::numbers::pi + (j + 0.5) * 2 * std::numbers::pi / nlon; }
  Vec3 point(std::size_t index) const {
    const int i = static_cast<int>(index / static_cast<std::size_t>(nlon)), j = static_cast<int>(index % static_cast<std::size_t>(nlon));
    const double p = latitude(i), l = longitude(j);
    return centre + radius * Vec3{std::cos(p) * std::cos(l), std::cos(p) * std::sin(l), std::sin(p)};
  }
};

struct PlaneGrid {
  Vec3 origin{}, u{}, v{};
  int nu = 0, nv = 0;

  std::size_t size() const { return static_cast<std::size_t>(nu) * static_cast<std::size_t>(nv); }
  Vec3 point(std::size_t index) const {
    const double s = (static_cast<double>(index % static_cast<std::size_t>(nu)) + 0.5) / nu;
    const double r = (static_cast<double>(index / static_cast<std::size_t>(nu)) + 0.5) / nv;
    return origin + s * u + r * v;
  }
};

struct VolumeGrid {
  Vec3 centre{};
  Mat3 axes = geometry::identity3();  // rows: the grid's three axes, unit vectors on the frame's axes
  Vec3 edges{};                       // the box's edge along each axis (m)
  std::array<int, 3> n{};             // cells along each axis

  std::size_t size() const { return static_cast<std::size_t>(n[0]) * static_cast<std::size_t>(n[1]) * static_cast<std::size_t>(n[2]); }
  std::size_t index(int i, int j, int k) const { return static_cast<std::size_t>(i) + static_cast<std::size_t>(n[0]) * (static_cast<std::size_t>(j) + static_cast<std::size_t>(n[1]) * static_cast<std::size_t>(k)); }
  // The offset of a cell's centre from the box's centre, along the grid's axes.
  Vec3 offset(int i, int j, int k) const {
    return {((i + 0.5) / n[0] - 0.5) * edges[0], ((j + 0.5) / n[1] - 0.5) * edges[1], ((k + 0.5) / n[2] - 0.5) * edges[2]};
  }
  // A point given along the grid's axes, on the frame's axes.
  Vec3 from_axes(const Vec3& o) const { return centre + geometry::transpose(axes) * o; }
  // A point of the frame, along the grid's axes.
  Vec3 to_axes(const Vec3& x) const { return axes * (x - centre); }
  Vec3 point(std::size_t index) const {
    const int i = static_cast<int>(index % static_cast<std::size_t>(n[0]));
    const std::size_t rest = index / static_cast<std::size_t>(n[0]);
    const int j = static_cast<int>(rest % static_cast<std::size_t>(n[1])), k = static_cast<int>(rest / static_cast<std::size_t>(n[1]));
    return from_axes(offset(i, j, k));
  }
};

struct PointList {
  std::vector<Vec3> points;
  std::size_t size() const { return points.size(); }
  Vec3 point(std::size_t index) const { return points[index]; }
};

template <class G>
concept Sampling = requires(const G& g, std::size_t i) {
  { g.size() } -> std::convertible_to<std::size_t>;
  { g.point(i) } -> std::convertible_to<Vec3>;
};

// The answers on a sampling, with their counts by quality and by reason.
template <class V>
struct Samples {
  std::vector<V> values;
  std::vector<Quality> quality;
  std::array<std::size_t, 4> counts{};                             // by Quality
  std::vector<std::pair<std::string_view, std::size_t>> reasons;  // each reason with the number of points it refused, in order of appearance
  std::vector<std::pair<std::string_view, std::size_t>> notes;    // each note with the number of valid points it qualified
  std::vector<std::uint16_t> why;                                 // per point: 0, or 1 + the index of its reason (refused) or note (valid)

  std::size_t size() const { return values.size(); }
  std::size_t valid() const { return counts[0]; }
  std::string_view reason(std::size_t i) const { return quality[i] != Quality::valid && why[i] ? reasons[why[i] - 1u].first : std::string_view{}; }
  std::string_view note(std::size_t i) const { return quality[i] == Quality::valid && why[i] ? notes[why[i] - 1u].first : std::string_view{}; }

  // Every point valid, without notes.
  static Samples all_valid(std::vector<V> v) {
    Samples s;
    s.quality.assign(v.size(), Quality::valid);
    s.why.assign(v.size(), 0);
    s.counts[0] = v.size();
    s.values = std::move(v);
    return s;
  }
};

namespace sampling_detail {
// Counts `what` in the list and returns 1 + its index, or 0 for no text.
inline std::uint16_t tally(std::vector<std::pair<std::string_view, std::size_t>>& list, std::string_view what) {
  if (what.empty()) return 0;
  for (std::size_t k = 0; k < list.size(); k++)
    if (list[k].first == what) {
      list[k].second++;
      return static_cast<std::uint16_t>(k + 1);
    }
  if (list.size() >= 65534) throw std::length_error("physis::fields: more than 65534 distinct reasons in one sampling");
  list.emplace_back(what, 1);
  return static_cast<std::uint16_t>(list.size());
}
}  // namespace sampling_detail

// Samples from evaluations, in order.
template <class V>
Samples<V> samples_of(std::span<const Evaluation<V>> answers) {
  Samples<V> s;
  s.values.resize(answers.size());
  s.quality.resize(answers.size());
  s.why.resize(answers.size());
  for (std::size_t i = 0; i < answers.size(); i++) {
    const Evaluation<V>& a = answers[i];
    s.values[i] = a.ok() ? a.value : V{};
    s.quality[i] = a.quality;
    s.counts[static_cast<std::size_t>(a.quality)]++;
    s.why[i] = a.ok() ? sampling_detail::tally(s.notes, a.note) : sampling_detail::tally(s.reasons, a.reason);
  }
  return s;
}

// The queries of a sampling at an instant.
template <Sampling G>
std::vector<Query<Vec3>> queries_of(const G& g, double t) {
  std::vector<Query<Vec3>> q(g.size());
  for (std::size_t i = 0; i < q.size(); i++) q[i] = {t, g.point(i)};
  return q;
}

// A field read on a sampling: by the representation's own method when it has one for this sampling.
template <PointEvaluable F, Sampling G>
  requires std::is_same_v<typename F::point_type, Vec3>
Samples<typename F::value_type> sample(const F& f, const G& g, double t) {
  if constexpr (requires { { f.sample(g, t) } -> std::same_as<Samples<typename F::value_type>>; }) {
    return f.sample(g, t);
  } else {
    const auto q = queries_of(g, t);
    const auto answers = evaluate_many(f, std::span<const Query<Vec3>>(q));
    return samples_of(std::span<const Evaluation<typename F::value_type>>(answers));
  }
}

}  // namespace physis::fields
