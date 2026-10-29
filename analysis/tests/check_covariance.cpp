// Checks of the covariance of random potential fields (analysis/covariance.hpp).
//
// - The covariance of V and of B at pairs of points (on and above the sphere, near and far, at a pole,
//   at one point) against the sum over the coefficients, Σₖ σₖ² bₖ(x) bₖ(y)ᵀ, each bₖ the field of one unit
//   coefficient by field_at: the addition theorem and the chain rule, to rounding.
// - At degree 720: the variance at the reference sphere is Σₙ E[Rₙ] (Lowes), and the field's
//   covariance is the mixed second difference of the potential's.
// - The table at one radius against the exact sums, for pairs from a metre to 15 km apart; and pairs at
//   another radius, or beyond the table, summed exactly.

#include <physis/analysis/covariance.hpp>
#include <physis/analysis/potential_field.hpp>

#include <calx/verification/check.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <vector>

using namespace physis::analysis;
using Point = PotentialCovariance::Point;
using Block = PotentialCovariance::Block;

namespace {

Point cartesian(double r, double lat, double lon) { return {r * std::cos(lat) * std::cos(lon), r * std::cos(lat) * std::sin(lon), r * std::sin(lat)}; }

// The field of field_at in Cartesian components: north, east and down turned into the body frame.
Point field_cartesian(const HarmonicPotential& p, double r, double lat, double lon) {
  const FieldValue f = field_at(p, r, lat, lon);
  const double cp = std::cos(lat), sp = std::sin(lat), cl = std::cos(lon), sl = std::sin(lon);
  const Point north{-sp * cl, -sp * sl, cp}, east{-sl, cl, 0}, down{-cp * cl, -cp * sl, -sp};
  Point b;
  for (int i = 0; i < 3; i++) b[i] = f.north * north[i] + f.east * east[i] + f.down * down[i];
  return b;
}

struct Place {
  double r, lat, lon;
};

}  // namespace

int main() {
  calx::verification::Suite s;
  const double pi = std::numbers::pi, a = 1737.4e3;

  // ---- against the sum over the coefficients -------------------------------------------------------------------
  {
    const int N = 16;
    std::vector<double> spectrum(N + 1, 0.0);
    for (int n = 1; n <= N; n++) spectrum[n] = 1e4 / n * std::pow(0.9, 2.0 * n);  // nT²
    const PotentialCovariance K(spectrum, 2, N, a);
    const std::vector<std::array<Place, 2>> pairs{
        {{{a, 0.3, 1.1}, {1.05 * a, 0.25, 1.3}}},        // near, at two radii
        {{{a, 0.3, 1.1}, {a, 0.3, 1.1}}},                // one point: the variance
        {{{1.2 * a, -0.7, 2.0}, {a, 0.6, -1.4}}},        // far apart
        {{{a, pi / 2, 0.0}, {1.1 * a, 1.3, 0.4}}},       // a pole
        {{{a, 0.2, 0.5}, {1.01 * a, -0.2, 0.5 + pi}}}};  // nearly antipodal
    double field = 0, potential = 0;
    for (const auto& pr : pairs) {
      const auto& [x, y] = pr;
      Block expected{};
      double v = 0;
      for (int n = 2; n <= N; n++) {
        const double var = spectrum[n] / ((n + 1.0) * (2.0 * n + 1.0));
        for (int m = 0; m <= n; m++)
          for (int part = 0; part < (m ? 2 : 1); part++) {
            HarmonicPotential p = HarmonicPotential::zero(N, a);
            (part ? p.h : p.g)[HarmonicPotential::index(n, m)] = 1;
            const Point bx = field_cartesian(p, x.r, x.lat, x.lon), by = field_cartesian(p, y.r, y.lat, y.lon);
            for (int i = 0; i < 3; i++)
              for (int j = 0; j < 3; j++) expected[i][j] += var * bx[i] * by[j];
            v += var * field_at(p, x.r, x.lat, x.lon).potential * field_at(p, y.r, y.lat, y.lon).potential;
          }
      }
      const Block H = K.field(cartesian(x.r, x.lat, x.lon), cartesian(y.r, y.lat, y.lon));
      double scale = 0, worst = 0;
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) scale = std::max(scale, std::fabs(expected[i][j]));
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) worst = std::max(worst, std::fabs(H[i][j] - expected[i][j]) / scale);
      field = std::max(field, worst);
      potential = std::max(potential, std::fabs(K.potential(cartesian(x.r, x.lat, x.lon), cartesian(y.r, y.lat, y.lon)) / v - 1));
    }
    s.check("covariance: Cov(V(x), V(y)) = a² Σ σₙ² ρⁿ⁺¹ Pₙ(x̂·ŷ) is the sum over the coefficients (five pairs, relative)", potential, 1e-12);
    s.check("covariance: Cov(B(x), B(y)), its mixed Hessian, is Σₖ σₖ² bₖ(x) bₖ(y)ᵀ (five pairs: near, one point, far, a pole, antipodes; relative)", field, 1e-12);
  }

  // ---- at the degrees of a crustal field ---------------------------------------------------------------------------
  {
    const int N = 720;
    std::vector<double> spectrum(N + 1, 0.0);
    for (int n = 1; n <= N; n++) spectrum[n] = 30.0 * n * std::pow(0.995, 2.0 * n);  // nT², a crust-like law
    const PotentialCovariance K(spectrum, 1, N, a);
    double sum = 0;
    for (int n = 1; n <= N; n++) sum += spectrum[n];
    const Point x = cartesian(a, 0.4, -0.8);
    s.check("covariance: at degree 720 the variance E|B|² at the reference sphere is Σ E[Rₙ] (Lowes; relative)", std::fabs(K.total_variance(x) / sum - 1), 1e-12);
    // The field's covariance against the mixed second difference of the potential's, 3 km apart, 2 km up.
    const Point y = cartesian(a + 2e3, 0.4 + 3e3 / a, -0.8);
    const Block H = K.field(x, y);
    const double h = 5.0;  // m: the truncation, of order (h n/a)² ≈ 4 × 10⁻⁶ at n = 720, above the rounding (10⁻⁹)
    double worst = 0, scale = 0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) scale = std::max(scale, std::fabs(H[i][j]));
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        Point xp = x, xm = x, yp = y, ym = y;
        xp[i] += h, xm[i] -= h, yp[j] += h, ym[j] -= h;
        const double d = (K.potential(xp, yp) - K.potential(xp, ym) - K.potential(xm, yp) + K.potential(xm, ym)) / (4 * h * h);
        worst = std::max(worst, std::fabs(H[i][j] - d) / scale);
      }
    s.check("covariance: at degree 720 Cov(B(x), B(y)) is the mixed second difference of Cov(V(x), V(y)) (relative)", worst, 1e-5);

    // The table at a mast's height, 0.01 rad (17 km) in 4000 steps: N δ = 1.8 × 10⁻³.
    const double r = a + 1.5;
    const TabulatedCovariance T(K, r, 0.01, 4000);
    const double variance = K.total_variance(cartesian(r, 0, 0));
    double table = 0;
    for (int k = 0; k <= 200; k++) {
      const double d = k == 0 ? 0.0 : std::pow(10.0, 4.2 * k / 200.0);  // 1 m to 15.8 km, and the point itself
      const Point p = cartesian(r, 0.3, 1.1), q = cartesian(r, 0.3 + d / r * std::cos(0.7 * k), 1.1 + d / r * std::sin(0.7 * k) / std::cos(0.3));
      const Block A = T.field(p, q), B = K.field(p, q);
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) table = std::max(table, std::fabs(A[i][j] - B[i][j]) / variance);
    }
    s.check("table: at one radius, interpolated in the angle, the covariance is the exact sums' (1 m to 16 km; over the variance)", table, 1e-10);
    const Point p = cartesian(r, 0.3, 1.1);
    const Block far = T.field(p, cartesian(r, 0.33, 1.1)), other = T.field(p, cartesian(r + 50, 0.3001, 1.1));
    const Block far_exact = K.field(p, cartesian(r, 0.33, 1.1)), other_exact = K.field(p, cartesian(r + 50, 0.3001, 1.1));
    bool same = true;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) same = same && far[i][j] == far_exact[i][j] && other[i][j] == other_exact[i][j];
    s.require("table: a pair beyond it, or at another radius, is summed exactly", same);
  }

  return s.done();
}
