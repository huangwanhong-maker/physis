// Checks of the stationary random field read at any point (analysis/stationary_field.hpp).
//
// Over 4000 realisations (seeds) of a field on ℝ³ with σ = 2 and ℓ = 1.5: its mean is 0 and its variance
// σ²; its covariance at separations of ℓ and 3ℓ, along different axes, is the law's (isotropy); its
// kurtosis is near the normal law's 3. Each within four of its standard errors. A seed fixes a
// realisation and another stream gives another; read on the Moon's sphere, 1.7 × 10⁶ m from the origin of
// the coordinates, it keeps its law.

#include <physis/analysis/stationary_field.hpp>

#include <calx/verification/check.hpp>

#include <array>
#include <cmath>
#include <string>

using physis::analysis::StationaryField;
using P = StationaryField<3>::Point;

namespace {

struct Moments {
  double mean = 0, variance = 0, kurtosis = 0, cov_l = 0, cov_3l = 0;
};

// The moments over n realisations at x, and the covariances with x + ℓ e₁ and x + 3ℓ e₂.
Moments moments(int n, P x, P origin) {
  const double l = 1.5;
  double s1 = 0, s2 = 0, s4 = 0, c1 = 0, c3 = 0;
  for (int i = 1; i <= n; i++) {
    const StationaryField<3> g(2.0, l, std::uint64_t(i), 0, 256, origin);
    const double a = g(x), b = g({x[0] + l, x[1], x[2]}), c = g({x[0], x[1] + 3 * l, x[2]});
    s1 += a, s2 += a * a, s4 += a * a * a * a, c1 += a * b, c3 += a * c;
  }
  Moments m;
  m.mean = s1 / n;
  m.variance = s2 / n;
  m.kurtosis = (s4 / n) / (m.variance * m.variance);
  m.cov_l = c1 / n, m.cov_3l = c3 / n;
  return m;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  const int n = 4000;
  const double sigma2 = 4.0, se = 1 / std::sqrt(double(n));
  const Moments m = moments(n, {0.3, -0.7, 1.1}, {0, 0, 0});
  s.check("stationary: the mean over 4000 realisations is 0 (in σ/√n)", m.mean / (2.0 * se), 4);
  s.check("stationary: the variance is σ² (relative, in its standard error √(2/n))", (m.variance / sigma2 - 1) / (std::sqrt(2.0) * se), 4);
  const double r1 = std::exp(-0.5), r3 = std::exp(-4.5);
  s.check("stationary: the covariance at a separation ℓ along x is σ² e^{−1/2} (in its standard error σ²√((1 + ρ²)/n))",
          (m.cov_l / sigma2 - r1) / (std::sqrt(1 + r1 * r1) * se), 4);
  s.check("stationary: at 3ℓ along y, σ² e^{−9/2}: the law is isotropic and falls as it should (in its standard error)",
          (m.cov_3l / sigma2 - r3) / (std::sqrt(1 + r3 * r3) * se), 4);
  s.check("stationary: its kurtosis is near the normal law's 3 (256 waves; the excess of a sum of M is −3/(2M))", m.kurtosis - 3, 0.3);
  {
    const StationaryField<3> a(1.0, 1.0, 42, 1), b(1.0, 1.0, 42, 1), c(1.0, 1.0, 42, 2);
    const P x{0.25, 0.5, -1.0};
    s.require("stationary: a seed and a stream fix the realisation, the same on every machine", a(x) == b(x));
    s.require("stationary: another stream is another realisation", a(x) != c(x));
  }
  {
    const P far{1737.4e3, 0.0, 0.0};
    const Moments f = moments(n, far, far);
    s.check("stationary: read on the Moon's sphere about an origin there, its variance is still σ² (in its standard error)",
            (f.variance / sigma2 - 1) / (std::sqrt(2.0) * se), 4);
    s.check("stationary: and its covariance at ℓ (in its standard error)", (f.cov_l / sigma2 - r1) / (std::sqrt(1 + r1 * r1) * se), 4);
  }
  return s.done();
}
