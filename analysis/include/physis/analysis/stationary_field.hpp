// A stationary random field of a given covariance, read at any point of space (decision 0015).
//
// The law. A field g on ℝᴰ with mean 0 and the squared-exponential covariance
//
//   C(r) = σ² exp(−|r|² / 2ℓ²),
//
// ℓ its correlation length: the field changes on that scale and is smooth below it. Its spectral density is
// the normal density of covariance I/ℓ² (Bochner's theorem: C is σ² times the characteristic function of
// that normal law at r).
//
// The spectral representation (Shinozuka and Deodatis 1991, Appl. Mech. Rev. 44, 191–204): a realisation is
//
//   g(x) = σ √(2/M) Σₘ cos(kₘ·(x − x₀) + φₘ),
//
// with M wave vectors kₘ drawn from the spectral density and phases φₘ uniform on [0, 2π), all from one
// counter-based stream, so that a seed gives the same field on every machine. For every M its covariance
// over the realisations is C exactly: E_φ[cos(k·x + φ) cos(k·y + φ)] = ½ cos(k·(x − y)), and E_k of that is
// ½ C(x − y)/σ². Its law tends to the Gaussian as M grows (the central limit theorem); M = 256 by default.
//
// It is read at any point, with no grid: a field over a whole body's surface at a metre's scale costs M
// cosines a point. The origin x₀ (the stationary law does not depend on it) keeps the phases accurate far
// from the coordinates' origin: on the Moon's surface, |x| = 1.7 × 10⁶ m.
//
// Caveat. A body's surface is a sphere, not a plane; the field is the restriction of a field of space, so
// its covariance along the surface is C of the chord, which differs from C of the arc by a part in (ℓ/R)²,
// nothing at a metre's scale on a world.

#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <stdexcept>
#include <vector>

#include <calx/random/stream.hpp>

namespace physis::analysis {

template <std::size_t D>
class StationaryField {
 public:
  using Point = std::array<double, D>;

  // σ the standard deviation, ℓ the correlation length, `seed` and `stream` the realisation (a field of a
  // world takes the world's seed and a stream of its own), M the number of waves, x₀ the origin.
  StationaryField(double sigma, double length, std::uint64_t seed, std::uint32_t stream = 0, int waves = 256, Point origin = {})
      : sigma_(sigma), length_(length), origin_(origin) {
    if (!(sigma >= 0) || !(length > 0) || waves < 1) throw std::invalid_argument("physis::analysis::StationaryField: σ ≥ 0, ℓ > 0 and at least one wave");
    calx::random::Stream s(seed, 0x57A7u, stream, 0u);
    k_.resize(std::size_t(waves));
    phase_.resize(std::size_t(waves));
    for (int m = 0; m < waves; m++) {
      for (std::size_t i = 0; i < D; i++) k_[std::size_t(m)][i] = s.normal() / length;
      phase_[std::size_t(m)] = s.uniform(0.0, 2 * std::numbers::pi);
    }
    amplitude_ = sigma * std::sqrt(2.0 / waves);
  }

  double sigma() const { return sigma_; }
  double length() const { return length_; }
  int waves() const { return int(k_.size()); }
  const Point& origin() const { return origin_; }

  // The law's covariance at a separation r.
  double covariance(const Point& r) const {
    double r2 = 0;
    for (double c : r) r2 += c * c;
    return sigma_ * sigma_ * std::exp(-r2 / (2 * length_ * length_));
  }

  // The realisation at x.
  double operator()(const Point& x) const {
    Point d;
    for (std::size_t i = 0; i < D; i++) d[i] = x[i] - origin_[i];
    double g = 0;
    for (std::size_t m = 0; m < k_.size(); m++) {
      double a = phase_[m];
      for (std::size_t i = 0; i < D; i++) a += k_[m][i] * d[i];
      g += std::cos(a);
    }
    return amplitude_ * g;
  }

 private:
  double sigma_, length_, amplitude_ = 0;
  Point origin_;
  std::vector<Point> k_;
  std::vector<double> phase_;
};

}  // namespace physis::analysis
