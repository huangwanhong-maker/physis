// Cochains on boxes as fields, and their transfer to a refinement (decision 0012; docs/theory/fields.md
// §9; formalism Propositions 3.6 and 6.2).
//
// A k-cochain is a number on each k-cell: an integral of a k-form over it. It pairs with k-chains
// exactly, ⟨c, Σ aᵢ σᵢ⟩ = Σ aᵢ c(σᵢ); a value at a point is not among its data and needs a reconstruction,
// which a box gives for two degrees:
//   k = 0   the multilinear interpolation of the vertices of the cell holding the point (Whitney's 0-forms on
//           a cube, tensor products of the hat functions);
//   k = n   a density, the cell's value over its measure, constant in the cell.
// Other degrees are paired, not evaluated; a query at a point refuses, saying why.
//
// Transfer. On a box refined from another (every coarse cut among the fine ones), the de Rham map of the
// Whitney interpolant, T = R′ ∘ W, is the tensor product of its one-dimensional parts: along an axis a
// cell spans, a fine edge takes its share (its length over the coarse edge's) of the coarse edge holding
// it, the integral of a constant density; along an axis it does not span, a fine vertex takes the linear
// interpolation of the coarse vertices about it. Each part commutes with the one-dimensional d, so the
// product commutes with the box's d, which is their graded tensor product: T d = d T in every degree, and
// for top cochains T keeps the total (Proposition 6.2).

#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

#include <physis/discretisation/geometric.hpp>
#include <physis/fields/field.hpp>

namespace physis::discretisation {

namespace fields = physis::fields;

class CochainField {
 public:
  using value_type = double;
  using point_type = Vec3;
  static constexpr std::uint32_t extra_capabilities = static_cast<std::uint32_t>(fields::Capability::chain) |
                                                      static_cast<std::uint32_t>(fields::Capability::coboundary) |
                                                      static_cast<std::uint32_t>(fields::Capability::transfer);

  CochainField(fields::Descriptor d, std::shared_ptr<const BoxComplex> box, std::size_t degree, Cochain values)
      : d_(std::move(d)), box_(std::move(box)), k_(degree), c_(std::move(values)) {
    if (k_ > box_->dimension()) throw std::invalid_argument("physis::discretisation::CochainField: a degree beyond the box's dimension");
    if (c_.size() != box_->cells[k_].size()) throw std::invalid_argument("physis::discretisation::CochainField: one value per k-cell");
    d_.representation.kind = fields::RepresentationKind::cochain;
    d_.fibre.kind = fields::FibreKind::form;
    d_.fibre.degree = static_cast<int>(k_);
    d_.base.kind = fields::BaseKind::complex;
    d_.base.dimension = static_cast<int>(box_->dimension());
    if (box_->dimension() > 0) measures_ = measure(box_->complex);
  }

  const fields::Descriptor& descriptor() const { return d_; }
  const BoxComplex& box() const { return *box_; }
  std::shared_ptr<const BoxComplex> shared_box() const { return box_; }
  std::size_t degree() const { return k_; }
  const Cochain& values() const { return c_; }

  // ⟨c, chain⟩ for a k-chain given as (cell, coefficient) pairs.
  double pair(const std::vector<std::pair<std::size_t, double>>& chain) const {
    double s = 0;
    for (const auto& [cell, a] : chain) s += a * c_.at(cell);
    return s;
  }

  // dc, a (k + 1)-cochain field on the same box.
  CochainField coboundary() const {
    fields::Descriptor d = d_;
    d.name = "d(" + d_.name + ")";
    return CochainField(std::move(d), box_, k_ + 1, discretisation::coboundary(box_->complex.topology, k_).apply(c_));
  }

  fields::Evaluation<double> at(const fields::Query<Vec3>& q) const {
    const std::size_t n = box_->dimension();
    if (k_ != 0 && k_ != n) return fields::refused<double>(fields::Quality::invalid, "a cochain of this degree is paired with chains, not read at points");
    const auto L = box_->locate(q.x);
    if (!L) return fields::refused<double>(fields::Quality::outside, "outside the box");
    if (k_ == n) return fields::valid(c_[L->cell] / measures_.of(n)[L->cell], "a density, constant in its cell");
    // Multilinear in the cell: the vertices at the corner plus each subset of the spanned axes.
    double value = 0;
    const std::uint8_t axes = box_->full_axes();
    for (int corner = 0; corner < 8; corner++) {
      if ((corner & ~static_cast<int>(axes)) != 0) continue;
      std::array<int, 3> at = L->corner;
      double w = 1;
      for (int a = 0; a < 3; a++) {
        if (!(axes & (1u << a))) continue;
        const int up = (corner >> a) & 1;
        at[static_cast<std::size_t>(a)] += up;
        w *= up ? L->local[static_cast<std::size_t>(a)] : 1 - L->local[static_cast<std::size_t>(a)];
      }
      value += w * c_[box_->find(0, 0, at)];
    }
    return fields::valid(value);
  }

 private:
  fields::Descriptor d_;
  std::shared_ptr<const BoxComplex> box_;
  std::size_t k_;
  Cochain c_;
  Measures measures_;
};

// T: Cᵏ(coarse) → Cᵏ(fine), for a fine box whose cuts include the coarse box's.
inline Cochain transfer(const BoxComplex& coarse, const BoxComplex& fine, std::size_t k, const Cochain& c) {
  if (c.size() != coarse.cells[k].size()) throw std::invalid_argument("physis::discretisation::transfer: one value per coarse k-cell");
  // Along each axis, for each fine cut: the coarse interval holding it and the position in it.
  std::array<std::vector<std::pair<int, double>>, 3> at;
  for (int a = 0; a < 3; a++) {
    const auto& C = coarse.coordinates[static_cast<std::size_t>(a)];
    const auto& F = fine.coordinates[static_cast<std::size_t>(a)];
    if (C.front() != F.front() || C.back() != F.back()) throw std::invalid_argument("physis::discretisation::transfer: the boxes differ in extent");
    for (double x : F) {
      if (C.size() == 1) {
        at[static_cast<std::size_t>(a)].push_back({0, 0.0});
        continue;
      }
      const int i = std::clamp(static_cast<int>(std::upper_bound(C.begin(), C.end(), x) - C.begin()) - 1, 0, static_cast<int>(C.size()) - 2);
      at[static_cast<std::size_t>(a)].push_back({i, (x - C[static_cast<std::size_t>(i)]) / (C[static_cast<std::size_t>(i) + 1] - C[static_cast<std::size_t>(i)])});
    }
    for (double x : C)
      if (!std::binary_search(F.begin(), F.end(), x)) throw std::invalid_argument("physis::discretisation::transfer: a coarse cut not among the fine");
  }
  Cochain out(fine.cells[k].size(), 0.0);
  for (std::size_t f = 0; f < out.size(); f++) {
    const auto& cell = fine.cells[k][f];
    // Per axis: the coarse choices with their weights.
    std::array<std::array<std::pair<int, double>, 2>, 3> choices{};
    std::array<int, 3> count{};
    for (int a = 0; a < 3; a++) {
      const auto& A = at[static_cast<std::size_t>(a)];
      const int i = cell.corner[static_cast<std::size_t>(a)];
      if (cell.axes & (1u << a)) {
        // A fine edge [i, i + 1] inside the coarse edge holding its midpoint: its share of it.
        const auto& C = coarse.coordinates[static_cast<std::size_t>(a)];
        const auto& F = fine.coordinates[static_cast<std::size_t>(a)];
        const double mid = 0.5 * (F[static_cast<std::size_t>(i)] + F[static_cast<std::size_t>(i) + 1]);
        const int e = std::clamp(static_cast<int>(std::upper_bound(C.begin(), C.end(), mid) - C.begin()) - 1, 0, static_cast<int>(C.size()) - 2);
        choices[static_cast<std::size_t>(a)][0] = {e, (F[static_cast<std::size_t>(i) + 1] - F[static_cast<std::size_t>(i)]) / (C[static_cast<std::size_t>(e) + 1] - C[static_cast<std::size_t>(e)])};
        count[static_cast<std::size_t>(a)] = 1;
      } else {
        // A fine vertex: the linear interpolation of the coarse vertices about it.
        const auto [lo, s] = A[static_cast<std::size_t>(i)];
        if (coarse.coordinates[static_cast<std::size_t>(a)].size() == 1) {
          choices[static_cast<std::size_t>(a)][0] = {0, 1.0};
          count[static_cast<std::size_t>(a)] = 1;
        } else {
          choices[static_cast<std::size_t>(a)][0] = {lo, 1 - s};
          choices[static_cast<std::size_t>(a)][1] = {lo + 1, s};
          count[static_cast<std::size_t>(a)] = 2;
        }
      }
    }
    double v = 0;
    for (int p = 0; p < count[0]; p++)
      for (int q = 0; q < count[1]; q++)
        for (int r = 0; r < count[2]; r++) {
          const double w = choices[0][static_cast<std::size_t>(p)].second * choices[1][static_cast<std::size_t>(q)].second * choices[2][static_cast<std::size_t>(r)].second;
          if (w == 0) continue;
          const std::array<int, 3> corner{choices[0][static_cast<std::size_t>(p)].first, choices[1][static_cast<std::size_t>(q)].first, choices[2][static_cast<std::size_t>(r)].first};
          const std::size_t source = coarse.find(k, cell.axes, corner);
          if (source == no_cell) throw std::logic_error("physis::discretisation::transfer: no coarse cell under a fine one");
          v += w * c[source];
        }
    out[f] = v;
  }
  return out;
}

}  // namespace physis::discretisation
