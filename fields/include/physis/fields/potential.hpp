// Potential fields outside a sphere as fields (decision 0012; docs/theory/fields.md §3): the spectral
// representation of analysis's internal expansions.
//
// A PotentialField is B = −∇V of a HarmonicPotential (physis analysis), on the shell r ≥ inner about
// the expansion's centre (by default r ≥ a, its reference sphere, outside which it holds), with its
// components on the axes of the body-fixed frame the expansion is written in. analysis gives north,
// east and down at a point of planetocentric latitude ψ and east longitude λ; on the frame's axes
//
//   B = X n̂ + Y ê − Z û,   û = (cos ψ cos λ, cos ψ sin λ, sin ψ),
//                          ê = (−sin λ, cos λ, 0),
//                          n̂ = (−sin ψ cos λ, −sin ψ sin λ, cos ψ),
//
// at the poles with the meridian λ = atan2(y, x) the point gives (analysis's convention there). On a
// sphere about its centre it is read row by row (analysis field_grid: the Legendre functions once per
// latitude), O(L N²) + O(L M N) for L rows and M columns instead of O(L M N²).

#pragma once

#include <cmath>
#include <memory>
#include <numbers>
#include <span>
#include <utility>

#include <physis/analysis/potential_field.hpp>
#include <physis/fields/domains.hpp>
#include <physis/fields/field.hpp>
#include <physis/fields/sampling.hpp>

namespace physis::fields {

class PotentialField {
 public:
  using value_type = Vec3;
  using point_type = Vec3;
  static constexpr std::uint32_t extra_capabilities = static_cast<std::uint32_t>(Capability::sphere_grid);

  // Below the reference sphere a query is refused: the expansion does not hold inside its sources. The
  // expansion is shared, not copied, among the copies of the field.
  PotentialField(Descriptor d, std::shared_ptr<const analysis::HarmonicPotential> p, int to_degree = -1)
      : d_(std::move(d)), p_(std::move(p)), degree_(to_degree < 0 ? p_->degree : std::min(to_degree, p_->degree)) {
    d_.representation.kind = RepresentationKind::spectral;
    d_.fibre.kind = FibreKind::vector;
    if (d_.fibre.basis.empty()) d_.fibre.basis = "frame";
    domain_.inner = p_->radius;
    domain_.below = "below the expansion's reference sphere, inside its sources";
  }
  PotentialField(Descriptor d, analysis::HarmonicPotential p, int to_degree = -1)
      : PotentialField(std::move(d), std::make_shared<const analysis::HarmonicPotential>(std::move(p)), to_degree) {}

  const Descriptor& descriptor() const { return d_; }
  const analysis::HarmonicPotential& expansion() const { return *p_; }
  int degree() const { return degree_; }

  Evaluation<Vec3> at(const Query<Vec3>& q) const {
    const Cover c = domain_.covers(q);
    if (!c.ok()) return refused<Vec3>(c.quality, c.reason);
    const double r = geometry::norm(q.x);
    const double psi = std::atan2(q.x[2], std::hypot(q.x[0], q.x[1])), lambda = std::atan2(q.x[1], q.x[0]);
    return valid(cartesian(analysis::field_at(*p_, r, psi, lambda, degree_), psi, lambda));
  }

  // The potential V at a point (the field's unit times metres).
  Evaluation<double> potential(const Query<Vec3>& q) const {
    const Cover c = domain_.covers(q);
    if (!c.ok()) return refused<double>(c.quality, c.reason);
    const double r = geometry::norm(q.x);
    return valid(analysis::field_at(*p_, r, std::atan2(q.x[2], std::hypot(q.x[0], q.x[1])), std::atan2(q.x[1], q.x[0]), degree_).potential);
  }

  // A sphere about the expansion's centre, row by row.
  Samples<Vec3> sample(const SphereGrid& g, double t) const {
    if (g.centre != Vec3{} || !(g.radius >= p_->radius)) return generic(g, t);
    std::vector<double> lats(static_cast<std::size_t>(g.nlat)), lons(static_cast<std::size_t>(g.nlon));
    for (int i = 0; i < g.nlat; i++) lats[static_cast<std::size_t>(i)] = g.latitude(i);
    for (int j = 0; j < g.nlon; j++) lons[static_cast<std::size_t>(j)] = g.longitude(j);
    const auto grid = analysis::field_grid(*p_, g.radius, lats, lons, degree_);
    std::vector<Vec3> v(grid.size());
    for (std::size_t k = 0; k < grid.size(); k++) v[k] = cartesian(grid[k], lats[k / lons.size()], lons[k % lons.size()]);
    return Samples<Vec3>::all_valid(std::move(v));
  }
  template <Sampling G>
  Samples<Vec3> sample(const G& g, double t) const { return generic(g, t); }

  // North, east and down at (ψ, λ) on the frame's axes.
  static Vec3 cartesian(const analysis::FieldValue& f, double psi, double lambda) {
    const double sp = std::sin(psi), cp = std::cos(psi), sl = std::sin(lambda), cl = std::cos(lambda);
    const Vec3 up{cp * cl, cp * sl, sp}, east{-sl, cl, 0}, north{-sp * cl, -sp * sl, cp};
    return f.north * north + f.east * east - f.down * up;
  }

 private:
  template <Sampling G>
  Samples<Vec3> generic(const G& g, double t) const {
    const auto q = queries_of(g, t);
    std::vector<Evaluation<Vec3>> a(q.size());
    for (std::size_t i = 0; i < q.size(); i++) a[i] = at(q[i]);
    return samples_of(std::span<const Evaluation<Vec3>>(a));
  }

  Descriptor d_;
  std::shared_ptr<const analysis::HarmonicPotential> p_;
  int degree_;
  Shell domain_;
};

}  // namespace physis::fields
