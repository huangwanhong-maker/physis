// Minkowski space ℝ^{1,d-1}: its metric, its causal structure and proper time.
//
// Minkowski space of dimension d is ℝ^d with the flat metric η = diag(-1, 1, ..., 1) in
// inertial coordinates (x^0, x^1, ..., x^{d-1}), x^0 = ct a length (physis decision 0005,
// §2.4). The squared interval between events p and q is
//   σ²(p, q) = η_ab (q - p)^a (q - p)^b = -(Δx^0)² + Σ_i (Δx^i)²,
// and q lies in the chronological future I⁺(p) when σ² < 0 and Δx^0 > 0, on the future
// light cone when σ² = 0 and Δx^0 > 0, and is spacelike to p when σ² > 0 (Naber 2012, The
// Geometry of Minkowski Spacetime, 2nd ed., §1.3-1.4; Wald 1984, §8.1). The causal future
// J⁺(p) is the union of I⁺(p), its light cone and p. These relations are invariant under
// the orthochronous Lorentz group (physis `symmetry`), whose transformations are not
// written here.
//
// The proper time along a timelike worldline is the length ∫ sqrt(-η(ẋ, ẋ)) dλ; along a
// polyline of events it is the sum over segments of sqrt(-σ²), and since x^0 = ct it is the
// length cτ, in metres.
//
// The Alexandrov interval (causal diamond) I(p, q) = I⁺(p) ∩ I⁻(q) of two events at proper
// time τ has volume ζ_d τ^d, ζ_d = π^{(d-1)/2} / (2^{d-1} d Γ((d+1)/2)) (Myrheim 1978;
// Meyer 1988): τ²/2 in two dimensions, π τ⁴ / 24 in four.

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace physis::spacetime {

// An event in inertial coordinates (ct, x^1, ..., x^{d-1}), metres.
template <std::size_t d>
using Event = std::array<double, d>;

// The flat metric of ℝ^{1,d-1} in inertial coordinates, as a metric field of physis
// `geometry`.
template <std::size_t d>
struct Minkowski {
  static constexpr std::size_t n = d;
  template <class T>
  std::array<std::array<T, n>, n> operator()(const std::array<T, n>&) const {
    std::array<std::array<T, n>, n> g{};
    g[0][0] = T(-1.0);
    for (std::size_t i = 1; i < n; i++) g[i][i] = T(1.0);
    return g;
  }
};

// σ²(p, q) = -(Δx^0)² + Σ (Δx^i)², m².
template <std::size_t d>
inline double interval(const Event<d>& p, const Event<d>& q) {
  double s = -(q[0] - p[0]) * (q[0] - p[0]);
  for (std::size_t i = 1; i < d; i++) s += (q[i] - p[i]) * (q[i] - p[i]);
  return s;
}

enum class Relation { equal, timelike_future, timelike_past, null_future, null_past, spacelike };

// Where q lies relative to p. Squared intervals within `tolerance` (m²) of zero count as null.
template <std::size_t d>
inline Relation relation(const Event<d>& p, const Event<d>& q, double tolerance = 0) {
  const double s = interval(p, q);
  const double dt = q[0] - p[0];
  if (s > tolerance) return Relation::spacelike;
  if (std::fabs(s) <= tolerance) {
    if (dt == 0) return Relation::equal;  // σ² = 0 with Δx^0 = 0 means q = p (within tolerance)
    return dt > 0 ? Relation::null_future : Relation::null_past;
  }
  return dt > 0 ? Relation::timelike_future : Relation::timelike_past;
}

// q ∈ I⁺(p): q is in the chronological future of p.
template <std::size_t d>
inline bool chronologically_precedes(const Event<d>& p, const Event<d>& q) {
  return q[0] > p[0] && interval(p, q) < 0;
}

// q ∈ J⁺(p) \ {p}: q is in the causal future of p and distinct from it. The comparison
// (Δx^0)² >= |Δx|² is made without square roots, so it is exact for representable inputs
// up to the rounding of the squares.
template <std::size_t d>
inline bool causally_precedes(const Event<d>& p, const Event<d>& q) {
  const double dt = q[0] - p[0];
  if (dt <= 0) return false;
  double r2 = 0;
  for (std::size_t i = 1; i < d; i++) r2 += (q[i] - p[i]) * (q[i] - p[i]);
  return r2 <= dt * dt;
}

// The proper time cτ (m) along a polyline of events, each future-directed and causal from the
// one before (a null segment contributes nothing); throws otherwise.
template <std::size_t d>
inline double proper_time(const std::vector<Event<d>>& worldline) {
  double tau = 0;
  for (std::size_t k = 1; k < worldline.size(); k++) {
    const double s = interval(worldline[k - 1], worldline[k]);
    if (worldline[k][0] < worldline[k - 1][0] || s > 0)
      throw std::domain_error("physis::spacetime::proper_time: a segment of the worldline is not future-directed and causal");
    tau += std::sqrt(-s);
  }
  return tau;
}

// ζ_d, the volume of the Alexandrov interval of unit proper time in ℝ^{1,d-1}.
inline double alexandrov_constant(std::size_t d) {
  const double dd = static_cast<double>(d);
  return std::pow(std::numbers::pi, (dd - 1) / 2) / (std::pow(2.0, dd - 1) * dd * std::tgamma((dd + 1) / 2));
}

// The volume of the Alexandrov interval of proper time tau (m^d).
inline double alexandrov_volume(std::size_t d, double tau) { return alexandrov_constant(d) * std::pow(tau, static_cast<double>(d)); }

}  // namespace physis::spacetime
