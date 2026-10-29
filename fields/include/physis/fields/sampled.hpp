// Fields held as samples on a volume grid, with a stated reconstruction (decision 0012;
// docs/theory/fields.md §5).
//
// Samples are not a field until a reconstruction says what the field is between them. Here it is
// trilinear interpolation between the cells' centres, which reproduces every affine field exactly; in
// the outer half of the boundary cells, between the last centres and the box's faces, the value is
// extended from the nearest centres along each axis, and says so; beyond the box there is none. A point
// whose eight surrounding samples are not all valid has no value, with the reason of the first that is
// not.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>

#include <physis/fields/field.hpp>
#include <physis/fields/sampling.hpp>

namespace physis::fields {

template <class V>
class Sampled {
 public:
  using value_type = V;
  using point_type = Vec3;

  Sampled(Descriptor d, VolumeGrid grid, Samples<V> samples) : d_(std::move(d)), g_(grid), s_(std::move(samples)) {
    if (s_.size() != g_.size()) throw std::invalid_argument("physis::fields::Sampled: one sample for each cell of the grid");
    d_.representation.kind = RepresentationKind::sampled;
    if (d_.representation.about.empty()) d_.representation.about = "samples at cells' centres, trilinear between them";
  }

  const Descriptor& descriptor() const { return d_; }
  const VolumeGrid& grid() const { return g_; }
  const Samples<V>& samples() const { return s_; }

  Evaluation<V> at(const Query<Vec3>& q) const {
    const Vec3 o = g_.to_axes(q.x);
    std::array<int, 3> i0{};
    std::array<double, 3> w{};
    bool held = false;
    for (int a = 0; a < 3; a++) {
      const int n = g_.n[static_cast<std::size_t>(a)];
      // The coordinate in cells, the centres at 0 … n − 1.
      const double u = o[static_cast<std::size_t>(a)] / g_.edges[static_cast<std::size_t>(a)] * n + n / 2.0 - 0.5;
      if (!(u >= -0.5 && u <= n - 0.5)) return refused<V>(Quality::outside, "outside the sampled box");
      double c = u;
      if (c < 0 || c > n - 1) {
        c = std::clamp(c, 0.0, static_cast<double>(n - 1));
        held = true;
      }
      const int lo = n == 1 ? 0 : std::min(static_cast<int>(std::floor(c)), n - 2);
      i0[static_cast<std::size_t>(a)] = lo;
      w[static_cast<std::size_t>(a)] = n == 1 ? 0.0 : c - lo;
    }
    V value{};
    for (int corner = 0; corner < 8; corner++) {
      std::array<int, 3> idx{};
      double weight = 1;
      for (int a = 0; a < 3; a++) {
        const int up = (corner >> a) & 1;
        const int n = g_.n[static_cast<std::size_t>(a)];
        idx[static_cast<std::size_t>(a)] = std::min(i0[static_cast<std::size_t>(a)] + up, n - 1);
        weight *= up ? w[static_cast<std::size_t>(a)] : 1 - w[static_cast<std::size_t>(a)];
      }
      if (weight == 0) continue;
      const std::size_t k = g_.index(idx[0], idx[1], idx[2]);
      if (s_.quality[k] != Quality::valid) {
        const std::string_view why = s_.reason(k);
        return refused<V>(s_.quality[k], why.empty() ? std::string_view("a sample about the point has no value") : why);
      }
      if constexpr (std::is_same_v<V, double>) value += weight * s_.values[k];
      else
        for (int c = 0; c < 3; c++) value[static_cast<std::size_t>(c)] += weight * s_.values[k][static_cast<std::size_t>(c)];
    }
    return valid(value, held ? std::string_view("held from the nearest samples at the box's edge") : std::string_view{});
  }

 private:
  Descriptor d_;
  VolumeGrid g_;
  Samples<V> s_;
};

}  // namespace physis::fields
