// Monte Carlo updates of lattice gauge fields, each leaving e^{-S} dU invariant.
//
// - Metropolis (Metropolis et al. 1953; Creutz, Jacobs and Rebbi 1979, Phys. Rev. D 20,
//   1915): U_e → U' drawn symmetrically near U_e (a random element R near 1 times U_e, or
//   another element of a finite group), accepted with probability min(1, e^{-ΔS}); ΔS is
//   recomputed from the holonomies of the 2-cells containing e, so it is exact for any
//   complex, also when e occurs twice in one word.
// - Heat bath, where the local action is linear, -β/N Re tr(U Σ): U is drawn from
//   e^{β/N Re tr(U Σ)} dU directly.
//   · U(1): with Σ = r e^{iψ}, φ = θ + ψ has the von Mises density e^{κ cos φ}, κ = β r,
//     sampled by Best and Fisher's rejection (1979, Appl. Statist. 28, 152–157).
//   · SU(2): with Σ = k V (a quaternion of norm k, V ∈ SU(2)), X = U V has density
//     √(1 - x0²) e^{a x0} in x0, a = β k, and is uniform in the other directions; x0 is
//     drawn by Kennedy and Pendleton's method (1985, Phys. Lett. B 156, 393) for a >= 2 and by
//     Creutz's (1980, Phys. Rev. D 21, 2308) below, both exact; then U = X V⁻¹.
//   · SU(N): Cabibbo and Marinari's method (1982, Phys. Lett. B 119, 387): U ← R U for an
//     SU(2) element R acting on rows (i, j), for every pair i < j, each R drawn by the SU(2)
//     heat bath with the quaternion projection q of the (i, j) block of U Σ, since
//     Re tr(R U Σ) = 2 (r q)_0 + const and a = 2 β |q| / N.
//   Edges where the action is not linear are updated by Metropolis.
// - Heat bath of a finite group (ℤ_N, or any finite group by its table): U_e drawn among the |G|
//   elements with probability ∝ e^{β Σ_{f ∋ e} w_f Re tr U_f / N}, each weight computed by setting
//   the element in turn; exact for any complex (an edge met twice in a word included) and any
//   couplings. It is also aperiodic, which Metropolis on ℤ₂ is not: its one proposal is the flip,
//   so a link whose weight does not depend on it (in no 2-cell, or where its occurrences cancel)
//   flips at every sweep, and a loop of n such links is multiplied by (−1)ⁿ at every sweep: the
//   chain is periodic, and a loop of even length never changes.
// - Overrelaxation (Adler 1981; Creutz 1987, Phys. Rev. D 36, 515), microcanonical: U ↦ the
//   reflection that keeps Re tr(U Σ): for U(1), θ ↦ -θ - 2ψ; for SU(2), U ↦ V⁻¹ U⁻¹ V⁻¹. It
//   is not ergodic alone and is mixed with heat-bath sweeps.

#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <numbers>
#include <type_traits>
#include <vector>

#include "field.hpp"

namespace physis::gauge {

struct SweepStatistics {
  std::size_t proposed = 0, accepted = 0;
  double acceptance() const { return proposed ? static_cast<double>(accepted) / static_cast<double>(proposed) : 1.0; }
  SweepStatistics& operator+=(const SweepStatistics& o) {
    proposed += o.proposed, accepted += o.accepted;
    return *this;
  }
};

// One Metropolis update of link e, with `hits` proposals.
template <GaugeGroup G>
SweepStatistics metropolis_link(GaugeField<G>& U, std::size_t e, double beta, Stream& z, double eps, int hits = 1) {
  SweepStatistics st;
  const auto& g = U.group();
  for (int h = 0; h < hits; h++) {
    const auto old = U.link(e);
    const double before = U.local_plaquettes(e);
    U.set(e, g.near(old, eps, z));
    const double dS = -beta * (U.local_plaquettes(e) - before);
    st.proposed++;
    if (dS <= 0 || z.uniform() < std::exp(-dS)) st.accepted++;
    else U.set(e, old);
  }
  return st;
}

template <GaugeGroup G>
SweepStatistics metropolis_sweep(GaugeField<G>& U, double beta, Stream& z, double eps = 0.5, int hits = 1) {
  SweepStatistics st;
  for (std::size_t e = 0; e < U.links().size(); e++)
    if (!U.frozen(e)) st += metropolis_link(U, e, beta, z, eps, hits);
  return st;
}

namespace detail {

// x with density ∝ e^{κ cos x} on (-π, π] (Best and Fisher 1979).
inline double von_mises(double kappa, Stream& z) {
  constexpr double pi = std::numbers::pi;
  if (kappa < 1e-8) return pi * (2 * z.uniform() - 1);
  const double tau = 1 + std::sqrt(1 + 4 * kappa * kappa);
  const double rho = (tau - std::sqrt(2 * tau)) / (2 * kappa);
  const double r = (1 + rho * rho) / (2 * rho);
  for (;;) {
    const double u1 = z.uniform(), u2 = z.uniform(), u3 = z.uniform();
    const double c0 = std::cos(pi * u1);
    const double f = (1 + r * c0) / (r + c0);
    const double c = kappa * (r - f);
    if (c * (2 - c) - u2 > 0 || std::log(c / u2) + 1 - c >= 0) return (u3 > 0.5 ? 1 : -1) * std::acos(std::clamp(f, -1.0, 1.0));
  }
}

// x0 in [-1, 1] with density ∝ √(1 - x0²) e^{a x0}.
inline double su2_x0(double a, Stream& z) {
  if (a < 2) {
    // Creutz: x0 from e^{a x0} by inversion, accepted with probability √(1 - x0²).
    for (;;) {
      const double y = z.uniform();
      const double x0 = a < 1e-10 ? 2 * y - 1 : 1 + std::log(y + (1 - y) * std::exp(-2 * a)) / a;
      const double u = z.uniform();
      if (u * u <= 1 - x0 * x0) return x0;
    }
  }
  // Kennedy and Pendleton.
  for (;;) {
    const double r1 = z.uniform(), r2 = z.uniform(), r3 = z.uniform(), r4 = z.uniform();
    const double c = std::cos(2 * std::numbers::pi * r2);
    const double l2 = -(std::log(r1) + c * c * std::log(r3)) / (2 * a);
    if (r4 * r4 <= 1 - l2) return 1 - 2 * l2;
  }
}

// X ∈ SU(2) with density ∝ e^{a X_0} dX.
inline symmetry::SU2 su2_heat_bath(double a, Stream& z) {
  const double x0 = su2_x0(a, z);
  const double r = std::sqrt(std::max(0.0, 1 - x0 * x0));
  // A uniform direction on the sphere.
  const double c = 2 * z.uniform() - 1, phi = 2 * std::numbers::pi * z.uniform();
  const double s = std::sqrt(std::max(0.0, 1 - c * c));
  return {x0, r * s * std::cos(phi), r * s * std::sin(phi), r * c};
}

}  // namespace detail

// Heat-bath draws of one link from its staple.
inline symmetry::U1 heat_bath_link(const U1Group&, const std::complex<double>& sigma, double beta, Stream& z) {
  const double psi = std::arg(sigma);
  return symmetry::U1::from_angle(detail::von_mises(beta * std::abs(sigma), z) - psi);
}
inline symmetry::SU2 heat_bath_link(const SU2Group& g, const symmetry::SU2& sigma, double beta, Stream& z) {
  const double k = sigma.norm();
  if (k < 1e-300) return g.random(z);
  const symmetry::SU2 V = (1 / k) * sigma;
  return (detail::su2_heat_bath(beta * k, z) * V.inverse()).reprojected();
}
inline symmetry::SUN heat_bath_link(const SUNGroup& g, symmetry::SUN U, const symmetry::SUN& sigma, double beta, Stream& z) {
  const std::size_t n = g.n;
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = i + 1; j < n; j++) {
      const symmetry::SUN W = U * sigma;
      const symmetry::Matrix2c w{{{W(i, i), W(i, j)}, {W(j, i), W(j, j)}}};
      const symmetry::SU2 q = symmetry::from_matrix(w);
      const double k = q.norm();
      symmetry::SU2 r;
      if (k < 1e-300) r = symmetry::haar_su2(z);
      else r = detail::su2_heat_bath(2 * beta * k / static_cast<double>(n), z) * ((1 / k) * q).inverse();
      const auto m = r.matrix();
      for (std::size_t c = 0; c < n; c++) {
        const auto ui = U(i, c), uj = U(j, c);
        U(i, c) = m[0][0] * ui + m[0][1] * uj;
        U(j, c) = m[1][0] * ui + m[1][1] * uj;
      }
    }
  return U.reunitarised();
}

// One heat-bath sweep; edges where the action is not linear, or whose 2-cells carry different
// couplings, are updated by Metropolis; frozen links are left alone.
template <StapleGroup G>
SweepStatistics heat_bath_sweep(GaugeField<G>& U, double beta, Stream& z, double eps = 0.5) {
  SweepStatistics st;
  for (std::size_t e = 0; e < U.links().size(); e++) {
    if (U.frozen(e)) continue;
    const double w = U.link_coupling(e);
    if (!U.linear(e) || std::isnan(w)) {
      st += metropolis_link(U, e, beta, z, eps, 4);
      continue;
    }
    if (U.occurrences(e).empty()) {
      U.set(e, U.group().random(z));  // a link in no 2-cell is free: Haar
      continue;
    }
    const auto sigma = U.staple(e);
    if constexpr (std::is_same_v<G, SUNGroup>) U.set(e, heat_bath_link(U.group(), U.link(e), sigma, beta * w, z));
    else U.set(e, heat_bath_link(U.group(), sigma, beta * w, z));
    st.proposed++, st.accepted++;
  }
  return st;
}

// One heat-bath sweep of a finite group; frozen links are left alone.
template <FiniteGaugeGroup G>
SweepStatistics heat_bath_sweep(GaugeField<G>& U, double beta, Stream& z, double = 0.5) {
  SweepStatistics st;
  const std::size_t n = U.group().order();
  std::vector<double> w(n);
  for (std::size_t e = 0; e < U.links().size(); e++) {
    if (U.frozen(e)) continue;
    double top = -std::numeric_limits<double>::infinity();
    for (std::size_t g = 0; g < n; g++) {
      U.set(e, static_cast<typename G::Element>(g));
      w[g] = beta * U.local_plaquettes(e);
      top = std::max(top, w[g]);
    }
    double total = 0;
    for (double& x : w) total += x = std::exp(x - top);
    std::size_t g = 0;
    for (double y = z.uniform() * total; g + 1 < n && y >= w[g]; g++) y -= w[g];
    U.set(e, static_cast<typename G::Element>(g));
    st.proposed++, st.accepted++;
  }
  return st;
}

// One overrelaxation sweep (U(1) and SU(2)), over the linear edges.
inline void overrelaxation_link(GaugeField<U1Group>& U, std::size_t e) {
  const auto sigma = U.staple(e);
  if (std::abs(sigma) < 1e-300) return;
  U.set(e, symmetry::U1::from_angle(-U.link(e).theta - 2 * std::arg(sigma)));
}
inline void overrelaxation_link(GaugeField<SU2Group>& U, std::size_t e) {
  const auto sigma = U.staple(e);
  const double k = sigma.norm();
  if (k < 1e-300) return;
  const symmetry::SU2 Vi = ((1 / k) * sigma).inverse();
  U.set(e, (Vi * U.link(e).inverse() * Vi).reprojected());
}
template <class G>
void overrelaxation_sweep(GaugeField<G>& U)
  requires std::is_same_v<G, U1Group> || std::is_same_v<G, SU2Group>
{
  for (std::size_t e = 0; e < U.links().size(); e++)
    if (U.linear(e) && !U.frozen(e) && !U.occurrences(e).empty() && !std::isnan(U.link_coupling(e))) overrelaxation_link(U, e);
}

}  // namespace physis::gauge
