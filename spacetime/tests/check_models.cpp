// Checks of physis/spacetime's continuous models: Minkowski space's causal relations, the
// invariance of the interval under a boost, and proper time; the curvature of every exact
// solution against its closed forms (Kretschmann scalars of Schwarzschild, Reissner-Nordström
// and Kerr; Ricci-flatness of the vacuum solutions, inside the horizon too in Eddington-
// Finkelstein coordinates; the constant curvature of de Sitter and anti-de Sitter; the
// Friedmann equations from FLRW's Einstein tensor; the vanishing Weyl tensor of conformally
// flat spacetimes); and Schwarzschild's geodesics: a circular orbit, the relativistic advance
// of the periapsis against its exact quadrature and its weak-field limit 6πM/p, the photon
// sphere, and the conservation of energy, angular momentum and norm.

#include <physis/geometry/curvature.hpp>
#include <physis/geometry/geodesic.hpp>
#include <physis/spacetime/minkowski.hpp>
#include <physis/spacetime/models.hpp>
#include <calx/quadrature/gauss_kronrod.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <string>
#include <vector>

namespace pg = physis::geometry;
namespace ps = physis::spacetime;

namespace {

const double pi = std::numbers::pi;

template <std::size_t n>
double max_abs(const pg::Tensor2<n>& t) {
  double m = 0;
  for (const auto& r : t)
    for (double v : r) m = std::max(m, std::fabs(v));
  return m;
}
template <std::size_t n>
double max_abs(const pg::Tensor4<n>& t) {
  double m = 0;
  for (const auto& a : t)
    for (const auto& b : a)
      for (const auto& c : b)
        for (double v : c) m = std::max(m, std::fabs(v));
  return m;
}

// The largest |R_abcd - s (g_ac g_bd - g_ad g_bc) / ℓ²|, relative to the largest such component.
template <std::size_t n>
double constant_curvature_deviation(const pg::Curvature<n>& k, double sign, double ell) {
  const auto R = pg::lowered_riemann(k);
  double m = 0, size = 0;
  for (std::size_t a = 0; a < n; a++)
    for (std::size_t b = 0; b < n; b++)
      for (std::size_t c = 0; c < n; c++)
        for (std::size_t d = 0; d < n; d++) {
          const double want = sign * (k.g[a][c] * k.g[b][d] - k.g[a][d] * k.g[b][c]) / (ell * ell);
          m = std::max(m, std::fabs(R[a][b][c][d] - want));
          size = std::max(size, std::fabs(want));
        }
  return m / size;
}

// Ricci-flatness measured against cancellation: max |R_bd| / max Σ_a |R^a_bad|, the size of the
// terms whose sum R_bd is.
template <std::size_t n>
double ricci_cancellation(const pg::Curvature<n>& k) {
  double m = 0, size = 0;
  for (std::size_t b = 0; b < n; b++)
    for (std::size_t d = 0; d < n; d++) {
      double t = 0;
      for (std::size_t a = 0; a < n; a++) t += std::fabs(k.Riemann[a][b][a][d]);
      m = std::max(m, std::fabs(k.Ricci[b][d]));
      size = std::max(size, t);
    }
  return m / size;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- Minkowski space ---------------------------------------------------------------------------------------------------
  {
    using E = ps::Event<4>;
    const E o{0, 0, 0, 0};
    s.require("spacetime: relations to the origin: timelike, null and spacelike, future and past",
              ps::relation(o, E{2, 1, 0, 0}) == ps::Relation::timelike_future && ps::relation(o, E{-2, 0, 1, 1}) == ps::Relation::timelike_past &&
                  ps::relation(o, E{5, 3, 4, 0}) == ps::Relation::null_future && ps::relation(o, E{-1, 0, 0, 1}) == ps::Relation::null_past &&
                  ps::relation(o, E{1, 2, 0, 0}) == ps::Relation::spacelike && ps::relation(o, o) == ps::Relation::equal);
    // A boost of rapidity 0.7 along an oblique direction preserves every interval.
    const double eta = 0.7, ch = std::cosh(eta), sh = std::sinh(eta);
    const double nx = 0.6, ny = 0.0, nz = 0.8;
    auto boost = [&](const E& x) {
      const double par = nx * x[1] + ny * x[2] + nz * x[3];
      const double t = ch * x[0] - sh * par, p = -sh * x[0] + ch * par;
      return E{t, x[1] + (p - par) * nx, x[2] + (p - par) * ny, x[3] + (p - par) * nz};
    };
    const E p{0.3, -1.2, 0.4, 2.0}, q{2.9, 0.1, -0.7, 1.1};
    s.check("spacetime: the interval is invariant under a boost (relative)", ps::interval(boost(p), boost(q)) / ps::interval(p, q) - 1, 1e-14);
    // The twins: out at speed 3/5 for 5 m of ct and back, against staying at home.
    const std::vector<ps::Event<2>> travelling{{0, 0}, {5, 3}, {10, 0}}, home{{0, 0}, {10, 0}};
    s.check("spacetime: the travelling twin ages cτ = 8 m against 10 m at home", ps::proper_time(travelling) - 8 + (ps::proper_time(home) - 10), 1e-14);
    s.check("spacetime: Alexandrov interval of unit proper time: area 1/2 in 2 dimensions", ps::alexandrov_constant(2) - 0.5, 1e-15);
    s.check("spacetime: and volume π/24 in 4 (relative)", ps::alexandrov_constant(4) / (pi / 24) - 1, 1e-15);
  }
  // ---- Schwarzschild, in both charts -------------------------------------------------------------------------------------
  {
    const ps::Schwarzschild bh{1.5};
    const ps::SchwarzschildEddingtonFinkelstein ef{1.5};
    double kdev = 0, ricci = 0, kdev_ef = 0, ricci_ef = 0;
    for (double r : {3.2, 5.0, 11.0, 40.0})
      for (double th : {0.4, 1.2, 2.5}) {
        const auto k = pg::curvature(bh, pg::Point<4>{0.0, r, th, 0.3});
        const double K = bh.kretschmann(r);
        kdev = std::max(kdev, std::fabs(pg::kretschmann(k) / K - 1));
        ricci = std::max(ricci, ricci_cancellation(k));
      }
    // Eddington-Finkelstein coordinates reach inside the horizon r = 2M = 3.
    for (double r : {0.8, 2.0, 3.0, 7.5}) {
      const auto k = pg::curvature(ef, pg::Point<4>{1.0, r, 1.0, 0.0});
      const double K = bh.kretschmann(r);
      kdev_ef = std::max(kdev_ef, std::fabs(pg::kretschmann(k) / K - 1));
      ricci_ef = std::max(ricci_ef, ricci_cancellation(k));
    }
    s.check("spacetime: Schwarzschild's Kretschmann scalar is 48 M^2 / r^6 (relative, 12 points)", kdev, 1e-12);
    s.check("spacetime: Schwarzschild is Ricci-flat (max |R_bd| / max Σ_a |R^a_bad|)", ricci, 1e-13);
    s.check("spacetime: in ingoing Eddington-Finkelstein coordinates too, inside the horizon (relative)", kdev_ef, 1e-12);
    s.check("spacetime: and Ricci-flat there (max |R_bd| / max Σ_a |R^a_bad|)", ricci_ef, 1e-13);
    const auto k = pg::curvature(bh, pg::Point<4>{0.0, 6.0, 1.0, 0.0});
    const auto C = pg::weyl(k);
    const auto R = pg::lowered_riemann(k);
    double dev = 0;
    for (std::size_t a = 0; a < 4; a++)
      for (std::size_t b = 0; b < 4; b++)
        for (std::size_t c = 0; c < 4; c++)
          for (std::size_t d = 0; d < 4; d++) dev = std::max(dev, std::fabs(C[a][b][c][d] - R[a][b][c][d]));
    s.check("spacetime: in vacuum the Weyl tensor is the Riemann tensor (max difference / max |R_abcd|)", dev / max_abs(R), 1e-14);
    s.require("spacetime: Schwarzschild is Lorentzian outside and, in Eddington-Finkelstein coordinates, inside the horizon",
              pg::negative_directions(bh(pg::Point<4>{0, 5, 1, 0})) == 1 && pg::negative_directions(ef(pg::Point<4>{0, 1, 1, 0})) == 1);
  }
  // ---- Reissner-Nordström and Kerr ---------------------------------------------------------------------------------------
  {
    const ps::ReissnerNordstrom rn{1.0, 0.6};
    double kdev = 0, scalar = 0, trace = 0;
    for (double r : {2.5, 4.0, 9.0}) {
      const auto k = pg::curvature(rn, pg::Point<4>{0.0, r, 0.8, 0.0});
      kdev = std::max(kdev, std::fabs(pg::kretschmann(k) / rn.kretschmann(r) - 1));
      scalar = std::max(scalar, std::fabs(k.scalar) * r * r);
      const auto C = pg::weyl(k);
      for (std::size_t b = 0; b < 4; b++)
        for (std::size_t d = 0; d < 4; d++) {
          double t = 0;
          for (std::size_t a = 0; a < 4; a++)
            for (std::size_t c = 0; c < 4; c++) t += k.ginv[a][c] * C[a][b][c][d];
          trace = std::max(trace, std::fabs(t) * r * r);
        }
    }
    s.check("spacetime: Reissner-Nordström's Kretschmann scalar 8 (6M²r² - 12MQ²r + 7Q⁴) / r⁸ (relative)", kdev, 1e-12);
    s.check("spacetime: its Ricci scalar vanishes (the field's stress is trace-free; |R| r^2)", scalar, 1e-13);
    s.check("spacetime: its Weyl tensor is trace-free (max |g^ac C_abcd| r^2)", trace, 1e-13);
    s.check("spacetime: its horizons r± = M ± sqrt(M² - Q²) make g_tt, a sum of terms of order 1, vanish", std::fabs(rn(pg::Point<4>{0, rn.outer_horizon(), 1, 0})[0][0]) + std::fabs(rn(pg::Point<4>{0, rn.inner_horizon(), 1, 0})[0][0]), 1e-14);

    const ps::Kerr kerr{1.0, 0.7};
    double kk = 0, kr = 0;
    for (double r : {2.2, 3.5, 8.0})
      for (double th : {0.3, 1.1, pi / 2, 2.6}) {
        const auto k = pg::curvature(kerr, pg::Point<4>{0.0, r, th, 0.0});
        const double K = kerr.kretschmann(r, th);
        kk = std::max(kk, std::fabs(pg::kretschmann(k) - K) / std::max(std::fabs(K), 1e-3 / std::pow(r, 6)));
        kr = std::max(kr, ricci_cancellation(k));
      }
    s.check("spacetime: Kerr's Kretschmann scalar against Henry's closed form (relative, 12 points)", kk, 1e-11);
    s.check("spacetime: Kerr is Ricci-flat (max |R_bd| / max Σ_a |R^a_bad|)", kr, 1e-13);
    s.check("spacetime: Kerr's outer horizon makes g^rr vanish (Δ = 0)", kerr.outer_horizon() * kerr.outer_horizon() - 2 * kerr.outer_horizon() + 0.49, 1e-15);
    s.require("spacetime: Kerr's innermost stable orbits: 6M without spin, M and 9M prograde and retrograde at a = M",
              std::fabs(ps::Kerr{1, 0}.isco() - 6) < 1e-12 && std::fabs(ps::Kerr{1, 1}.isco(true) - 1) < 1e-12 && std::fabs(ps::Kerr{1, 1}.isco(false) - 9) < 1e-12);
  }
  // ---- de Sitter and anti-de Sitter --------------------------------------------------------------------------------------
  {
    const double ell = 3.0;
    const auto ks = pg::curvature(ps::DeSitterStatic{ell}, pg::Point<4>{0.0, 1.7, 1.0, 0.2});
    const auto kf = pg::curvature(ps::DeSitterFlat{ell}, pg::Point<4>{0.8, 0.3, -1.0, 2.0});
    const auto ka = pg::curvature(ps::AntiDeSitterGlobal{ell}, pg::Point<4>{0.0, 4.2, 2.0, 0.0});
    s.check("spacetime: de Sitter, static patch: R_abcd = (g_ac g_bd - g_ad g_bc) / ℓ² (max deviation, relative)", constant_curvature_deviation(ks, 1, ell), 1e-14);
    s.check("spacetime: de Sitter, flat slicing: the same", constant_curvature_deviation(kf, 1, ell), 1e-14);
    s.check("spacetime: anti-de Sitter: R_abcd = -(g_ac g_bd - g_ad g_bc) / ℓ²", constant_curvature_deviation(ka, -1, ell), 1e-14);
    s.check("spacetime: scalar curvatures 12/ℓ² and -12/ℓ² (relative)", std::fabs(ks.scalar * ell * ell / 12 - 1) + std::fabs(kf.scalar * ell * ell / 12 - 1) + std::fabs(ka.scalar * ell * ell / 12 + 1), 1e-13);
    s.check("spacetime: de Sitter is conformally flat (max |C_abcd| ℓ²)", max_abs(pg::weyl(ks)) * ell * ell, 1e-14);
  }
  // ---- FLRW: the Friedmann equations ---------------------------------------------------------------------------------------
  {
    // a = (ct / t0)^p with curvature K: ȧ = p a / x, ä = p (p - 1) a / x², x = ct.
    const double t0 = 2.0, p = 2.0 / 3, K = 0.01;
    const ps::FLRW<ps::PowerLaw> u{K, ps::PowerLaw{t0, p}};
    double dev = 0, wdev = 0;
    for (double x : {0.7, 2.0, 5.0})
      for (double r : {0.5, 3.0}) {
        const auto k = pg::curvature(u, pg::Point<4>{x, r, 1.0, 0.4});
        const double a = std::pow(x / t0, p), H = p / x, acc = p * (p - 1) / (x * x);
        const double G00 = 3 * (H * H + K / (a * a)), Gii = -(2 * acc + H * H + K / (a * a));
        const auto Gm = pg::einstein_mixed(k);
        dev = std::max({dev, std::fabs(k.Einstein[0][0] / G00 - 1), std::fabs(Gm[1][1] / Gii - 1), std::fabs(Gm[2][2] / Gii - 1),
                        std::fabs(Gm[3][3] / Gii - 1), std::fabs(k.scalar / (6 * (acc + H * H + K / (a * a))) - 1)});
        wdev = std::max(wdev, max_abs(pg::weyl(k)) / std::fabs(k.scalar));
      }
    s.check("spacetime: FLRW (power law, K = 0.01/m²): G_00 = 3 (H² + K/a²), G^i_i = -(2ä/a + H² + K/a²), R = 6 (ä/a + H² + K/a²) (relative)", dev, 1e-12);
    s.check("spacetime: FLRW is conformally flat (max |C_abcd| / |R|)", wdev, 1e-13);
    const ps::FLRW<ps::Exponential> inflation{0.0, ps::Exponential{0.5}};
    const auto k = pg::curvature(inflation, pg::Point<4>{1.3, 2.0, 0.7, 0.0});
    s.check("spacetime: exponential expansion is de Sitter: R = 12 H² (relative)", k.scalar / (12 * 0.25) - 1, 1e-13);
  }
  // ---- Schwarzschild geodesics -------------------------------------------------------------------------------------------
  {
    const double M = 1.0;
    const ps::Schwarzschild bh{M};
    // A circular orbit at r = 10M: dφ/d(ct) = sqrt(M / r³), u^t = 1 / sqrt(1 - 3M/r).
    {
      const double r = 10, Omega = std::sqrt(M / (r * r * r)), ut = 1 / std::sqrt(1 - 3 * M / r);
      pg::Phase<4> y = pg::phase<4>({0, r, pi / 2, 0}, {ut, 0, 0, Omega * ut});
      const double tau = 2 * pi / Omega / ut;  // proper time of one revolution
      const std::size_t steps = 4000;
      double rdev = 0;
      pg::integrate_geodesic(bh, y, tau / steps, steps, [&](double, const pg::Phase<4>& z) { rdev = std::max(rdev, std::fabs(z[1] - r)); return true; });
      s.check("spacetime: a circular orbit at r = 10M keeps its radius (max |Δr| / M, RK4)", rdev, 1e-9);
      s.check("spacetime: and revolves once in ct = 2π / sqrt(M/r³) (relative)", std::fabs(y[3] - 2 * pi) + std::fabs(y[0] * Omega / (2 * pi) - 1), 1e-9);
    }
    // A precessing orbit of semi-latus rectum p = 30 M and eccentricity 0.4, from periapsis
    // (Cutler, Kennefick and Poisson 1994, Phys. Rev. D 50, 3816, eqs. 2.4-2.10).
    {
      const double P = 30, e = 0.4;
      const double E = std::sqrt(((P - 2) * (P - 2) - 4 * e * e) / (P * (P - 3 - e * e)));
      const double L = P * M / std::sqrt(P - 3 - e * e);
      const double rp = P * M / (1 + e);
      pg::Phase<4> y = pg::phase<4>({0, rp, pi / 2, 0}, {E / (1 - 2 * M / rp), 0, 0, L / (rp * rp)});
      s.check("spacetime: the orbit's initial velocity is timelike and normalised, |g(u,u) + 1|", pg::norm(bh, y) + 1, 1e-14);
      const double E0 = -pg::momentum(bh, y, 0), L0 = pg::momentum(bh, y, 3);
      double edrift = 0, ldrift = 0, ndrift = 0;
      // Integrate until u^r turns from negative to positive (the next periapsis), then locate
      // the turning point by the secant method on the step length.
      const double h = 2.0;
      bool left = false;
      pg::Phase<4> prev = y;
      for (int k = 0; k < 100000; k++) {
        prev = y;
        pg::geodesic_step(bh, y, h, pg::Stepper::dopri5, {1e-13, 1e-13});
        edrift = std::max(edrift, std::fabs(-pg::momentum(bh, y, 0) - E0));
        ldrift = std::max(ldrift, std::fabs(pg::momentum(bh, y, 3) - L0));
        ndrift = std::max(ndrift, std::fabs(pg::norm(bh, y) + 1));
        if (y[5] < 0) left = true;
        if (left && y[5] >= 0) break;
      }
      auto ur_after = [&](double dl) {
        pg::Phase<4> z = prev;
        pg::geodesic_step(bh, z, dl, pg::Stepper::dopri5, {1e-14, 1e-14});
        return z;
      };
      double a = 0, b = h, fa = prev[5], fb = y[5];
      for (int it = 0; it < 60 && std::fabs(b - a) > 1e-13; it++) {
        const double c = b - fb * (b - a) / (fb - fa);
        const double fc = ur_after(c)[5];
        a = b, fa = fb, b = c, fb = fc;
      }
      const double dphi = ur_after(b)[3];
      // The exact angle between periapses: ∫_0^{2π} sqrt(p / (p - 6 - 2e cos χ)) dχ.
      const auto exact = calx::quadrature::integrate([&](double chi) { return std::sqrt(P / (P - 6 - 2 * e * std::cos(chi))); }, 0, 2 * pi, 1e-15, 1e-15);
      s.check("spacetime: the advance of the periapsis (p = 30M, e = 0.4) against the exact quadrature (relative)", (dphi - 2 * pi) / (exact.value - 2 * pi) - 1, 1e-8);
      s.check("spacetime: energy -p_t conserved along the orbit", edrift, 1e-11);
      s.check("spacetime: angular momentum p_φ conserved", ldrift / L0, 1e-11);
      s.check("spacetime: g(u, u) = -1 kept", ndrift, 1e-11);
      // The weak field: (Δφ - 2π) / (6πM/p) = 1 + (18 + e²) / (4p) + O(1/p²), from the expansion of the
      // integrand: Δφ - 2π = 6π/p + (3π/2)(18 + e²)/p² + O(1/p³).
      const double Pw = 1e4;
      const auto weak = calx::quadrature::integrate([&](double chi) { return std::sqrt(Pw / (Pw - 6 - 2 * e * std::cos(chi))); }, 0, 2 * pi, 1e-16, 1e-16);
      s.check("spacetime: at p = 10^4 M the advance is 6πM/p to first order (relative, less the 1/p term)", (weak.value - 2 * pi) / (6 * pi / Pw) - 1 - (18 + e * e) / (4 * Pw), 1e-6);
    }
    // The photon sphere: a null circular orbit at r = 3M has no radial acceleration.
    {
      auto radial_acceleration = [&](double r) {
        const double f = 1 - 2 * M / r;
        const pg::Phase<4> y = pg::phase<4>({0, r, pi / 2, 0}, {1 / f, 0, 0, 1 / (r * std::sqrt(f))});  // null: f (u^t)² = r² (u^φ)²
        pg::Phase<4> dy;
        pg::geodesic_field(bh, y, dy);
        return std::pair{dy[5], pg::norm(bh, y)};
      };
      const auto [a3, n3] = radial_acceleration(bh.photon_sphere());
      const auto [a4, n4] = radial_acceleration(4 * M);
      s.check("spacetime: the photon sphere r = 3M: a null circular orbit, |d²r/dλ²|", a3, 1e-15);
      s.require("spacetime: and at r = 4M a null orbit is not circular", std::fabs(a4) > 1e-3 && std::fabs(n3) < 1e-14 && std::fabs(n4) < 1e-14, "d²r/dλ² = " + std::to_string(a4));
    }
  }
  return s.done();
}
