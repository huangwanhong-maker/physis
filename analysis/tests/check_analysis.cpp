// Checks of analysis: potential fields from their coefficients, their spectra, and random fields.
//
// - The field of a dipole against its closed form.
// - An expansion to degree 449 (a closed formula, make_fixtures.py) against pyshtools at points on and
//   above the sphere; within a hundredth of a degree of the poles, where pyshtools's horizontal
//   components lose 10⁻⁹, against the exact field in 40-digit arithmetic (degree 30) and pyshtools's
//   radial component; its Lowes–Mauersberger spectrum against pyshtools's.
// - The field against the gradient of the potential, by central differences; the grid against the
//   points; SHTOOLS's format read back.
// - Random fields: their spectrum against the law they were drawn from; one seed, one field.
// - The spectral law: recovered exactly from itself, and closely from a realisation.
// - The unresolved part: with an estimate that damps the high degrees, estimate and unresolved part
//   together carry the law's spectrum, and no power is added where the estimate already has it.

#include <physis/analysis/potential_field.hpp>
#include <physis/analysis/random_field.hpp>

#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

using namespace physis::analysis;

namespace {

constexpr double pi = std::numbers::pi, deg = pi / 180, a = 1737400.0;

HarmonicPotential formula(int degree) {
  HarmonicPotential p = HarmonicPotential::zero(degree, a);
  for (int n = 1; n <= degree; n++)
    for (int m = 0; m <= n; m++) {
      const std::size_t k = HarmonicPotential::index(n, m);
      p.g[k] = 1e-9 * std::sin(1.3 * n + 0.7 * m + 0.1) / std::pow(n + 1.0, 1.2);
      p.h[k] = m ? 1e-9 * std::cos(0.9 * n - 1.7 * m) / std::pow(n + 1.0, 1.2) : 0.0;
    }
  return p;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- a dipole -------------------------------------------------------------------------------------------
  {
    HarmonicPotential d = HarmonicPotential::zero(1, a);
    const double g10 = -3e-5;
    d.g[HarmonicPotential::index(1, 0)] = g10;
    double worst = 0;
    for (const double lat : {-1.5, -0.3, 0.0, 0.8, 1.57})
      for (const double r : {a, 2.5 * a}) {
        const FieldValue f = field_at(d, r, lat, 0.7);
        const double q3 = std::pow(a / r, 3);
        worst = std::max({worst, std::fabs(f.north + q3 * g10 * std::cos(lat)), std::fabs(f.east), std::fabs(f.down + 2 * q3 * g10 * std::sin(lat)),
                          std::fabs(f.potential - a * (a / r) * (a / r) * g10 * std::sin(lat)) / a});
      }
    s.check("dipole: X = −(a/r)³ g₁₀ cos ψ, Y = 0, Z = −2 (a/r)³ g₁₀ sin ψ, V = a (a/r)² g₁₀ sin ψ", worst / std::fabs(g10), 1e-15);
  }

  // ---- against pyshtools ------------------------------------------------------------------------------------
  {
    std::ifstream in(PHYSIS_ANALYSIS_FIXTURES "/potential-reference.txt");
    s.require("pyshtools: the reference values are there", bool(in));
    const HarmonicPotential low = formula(30), high = formula(449);
    std::string line, kind;
    double worst_low = 0, worst_high = 0, worst_polar = 0, worst_exact = 0, worst_spectrum = 0;
    int count = 0, exact = 0;
    const std::vector<double> R = lowes_spectrum(high);
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream f(line);
      f >> kind;
      if (kind == "field" || kind == "exact") {
        int n = 0;
        double lat = 0, lon = 0, r = 0, br = 0, bt = 0, bp = 0;
        f >> n >> lat >> lon >> r >> br >> bt >> bp;
        const FieldValue v = field_at(n == 30 ? low : high, r, lat * deg, lon * deg);
        const double scale = std::sqrt(br * br + bt * bt + bp * bp);
        const double radial = std::fabs(v.radial() - br) / scale;
        const double e = std::max({radial, std::fabs(-v.north - bt) / scale, std::fabs(v.east - bp) / scale});
        if (kind == "exact") {
          worst_exact = std::max(worst_exact, e);
          exact++;
        } else if (std::fabs(lat) > 89.9) {
          worst_polar = std::max(worst_polar, radial);  // pyshtools's horizontal components lose 10⁻⁹ there
          count++;
        } else {
          (n == 30 ? worst_low : worst_high) = std::max(n == 30 ? worst_low : worst_high, e);
          count++;
        }
      } else if (kind == "spectrum") {
        int n = 0;
        double value = 0;
        f >> n >> value;
        worst_spectrum = std::max(worst_spectrum, std::fabs(R[n] / value - 1));
      }
    }
    s.require("pyshtools: 48 points and 6 exact ones read", count == 48 && exact == 6, std::to_string(count) + ", " + std::to_string(exact));
    s.check("pyshtools: the field of degree 30 at 18 points on and above the sphere", worst_low, 1e-12);
    s.check("pyshtools: the field of degree 449 at 18 points (rounding grows as n²ε: 5×10⁻¹¹ at 449)", worst_high, 1e-10);
    s.check("pyshtools: B_r within a hundredth of a degree of the poles, at degrees 30 and 449", worst_polar, 1e-10);
    s.check("exact: the field of degree 30 within a hundredth of a degree of the poles, against 40-digit arithmetic", worst_exact, 1e-13);
    s.check("pyshtools: the Lowes–Mauersberger spectrum", worst_spectrum, 1e-12);
  }

  // ---- the field is the gradient of the potential ------------------------------------------------------------
  {
    const HarmonicPotential p = formula(30);
    double worst = 0;
    for (const double lat : {-1.2, 0.3, 1.1}) {
      const double lon = 2.0, r = a + 50e3, hr = 1.0, ha = 1e-6;
      const FieldValue f = field_at(p, r, lat, lon);
      const double dVdr = (field_at(p, r + hr, lat, lon).potential - field_at(p, r - hr, lat, lon).potential) / (2 * hr);
      const double dVdp = (field_at(p, r, lat + ha, lon).potential - field_at(p, r, lat - ha, lon).potential) / (2 * ha);
      const double dVdl = (field_at(p, r, lat, lon + ha).potential - field_at(p, r, lat, lon - ha).potential) / (2 * ha);
      const double scale = f.total();
      worst = std::max({worst, std::fabs(f.down - dVdr) / scale, std::fabs(f.north + dVdp / r) / scale,
                        std::fabs(f.east + dVdl / (r * std::cos(lat))) / scale});
    }
    s.check("gradient: (X, Y, Z) = −(∂V/r∂ψ, ∂V/r cos ψ ∂λ, −∂V/∂r) by central differences", worst, 1e-6);
  }

  // ---- the grid and the points ----------------------------------------------------------------------------
  {
    const HarmonicPotential p = formula(449);
    const std::vector<double> lats{-89.99 * deg, -40 * deg, 0.0, 33 * deg, 89.995 * deg}, lons{0.0, 1.0, 2.5, 4.0, 6.0};
    const std::vector<FieldValue> grid = field_grid(p, a + 10e3, lats, lons);
    double worst = 0;
    for (std::size_t i = 0; i < lats.size(); i++)
      for (std::size_t j = 0; j < lons.size(); j++) {
        const FieldValue f = field_at(p, a + 10e3, lats[i], lons[j]), g = grid[i * lons.size() + j];
        const double scale = f.total();
        worst = std::max({worst, std::fabs(f.north - g.north) / scale, std::fabs(f.east - g.east) / scale, std::fabs(f.down - g.down) / scale,
                          std::fabs(f.potential - g.potential) / ((a + 10e3) * scale)});
      }
    s.check("grid: the field by rows equals the field at the points, at degree 449", worst, 1e-10);
  }

  // ---- SHTOOLS's format --------------------------------------------------------------------------------------
  {
    const HarmonicPotential p = formula(3);
    const auto path = std::filesystem::temp_directory_path() / "physis-analysis-check.sh";
    {
      std::ofstream out(path);
      out.precision(17);
      out << "1737.4, 3\n";
      for (int n = 0; n <= 3; n++)
        for (int m = 0; m <= n; m++) out << n << " " << m << " " << p.g[HarmonicPotential::index(n, m)] * 1e9 << " " << p.h[HarmonicPotential::index(n, m)] * 1e9 << "\n";
    }
    const HarmonicPotential q = read_shtools(path.string(), 1e-9, true);
    double worst = 0;
    for (std::size_t k = 0; k < p.g.size(); k++) worst = std::max({worst, std::fabs(p.g[k] - q.g[k]), std::fabs(p.h[k] - q.h[k])});
    s.require("SHTOOLS format: the radius in kilometres and the degree read from the header", q.radius == a && q.degree == 3);
    s.check("SHTOOLS format: the coefficients read back (in nanotesla, scaled to tesla)", worst / 1e-9, 1e-15);
    std::filesystem::remove(path);
  }

  // ---- random fields -----------------------------------------------------------------------------------------
  {
    SpectralLaw law;
    law.amplitude = 2e-18, law.exponent = 1.0, law.ratio = 0.995;
    std::vector<double> target(401);
    for (int n = 1; n <= 400; n++) target[n] = law(n);
    const HarmonicPotential f = random_potential(target, 50, 400, a, 20261028, 1);
    const std::vector<double> R = lowes_spectrum(f);
    double mean = 0, var = 0;
    int K = 0;
    for (int n = 50; n <= 400; n++) {
      mean += R[n] / target[n];
      var += 2.0 / (2 * n + 1);
      K++;
    }
    mean /= K;
    const double sd = std::sqrt(var) / K;
    s.check("random field: its spectrum averages the law's over degrees 50–400 (within four standard deviations)", std::fabs(mean - 1) / (4 * sd), 1.0);
    s.require("random field: nothing below the first degree asked for", R[49] == 0 && R[1] == 0);
    const HarmonicPotential again = random_potential(target, 50, 400, a, 20261028, 1), other = random_potential(target, 50, 400, a, 20261028, 2);
    s.require("random field: one seed and key, one field; another key, another", again.g == f.g && again.h == f.h && other.g != f.g);

    // The law from itself, and from the realisation.
    const SpectralLaw exact = fit_spectral_law(target, 20, 380);
    s.check("spectral law: recovered from its own spectrum (A, β, ρ)",
            std::max({std::fabs(exact.amplitude / law.amplitude - 1), std::fabs(exact.exponent - law.exponent), std::fabs(exact.ratio - law.ratio)}), 1e-8);
    // From one realisation, ρ to within four of its standard errors: s² (XᵀX)⁻¹ for the columns
    // (1, log n, 2n), s the fit's residual.
    const SpectralLaw fitted = fit_spectral_law(R, 50, 400);
    calx::linalg::Matrix XtX(3);
    for (int n = 50; n <= 400; n++) {
      const double x[3] = {1.0, std::log(double(n)), 2.0 * n};
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) XtX(i, j) += x[i] * x[j];
    }
    const calx::linalg::Matrix inv = calx::linalg::solve(XtX, calx::linalg::Matrix::identity(3));
    const double se = fitted.rms_log_residual * std::sqrt(inv(2, 2));
    s.check("spectral law: log ρ from one realisation, within four standard errors", std::fabs(std::log(fitted.ratio) - std::log(law.ratio)) / (4 * se), 1.0);

    // The unresolved part of an estimate that damps the high degrees.
    HarmonicPotential estimate = f;
    for (int n = 1; n <= 400; n++) {
      const double w = 1.0 / (1.0 + std::pow(n / 200.0, 4));
      for (int m = 0; m <= n; m++) estimate.g[HarmonicPotential::index(n, m)] *= w, estimate.h[HarmonicPotential::index(n, m)] *= w;
    }
    std::vector<double> law_spectrum(601);
    for (int n = 1; n <= 600; n++) law_spectrum[n] = law(n);
    const std::vector<double> u = unresolved_spectrum(estimate, law, 600), Re = lowes_spectrum(estimate);
    double worst_sum = 0;
    for (int n = 50; n <= 600; n++) worst_sum = std::max(worst_sum, std::fabs((n <= 400 ? Re[n] : 0.0) + u[n] - std::max(law(n), n <= 400 ? Re[n] : 0.0)) / law(n));
    s.check("unresolved part: the estimate's power and the unresolved power make the law's, degree by degree", worst_sum, 1e-12);
    bool none_added = true;
    for (int n = 50; n <= 400; n++)
      if (Re[n] >= law(n) && u[n] != 0) none_added = false;
    s.require("unresolved part: no power added where the estimate already has the law's", none_added);
    const HarmonicPotential both = estimate + unresolved_part(estimate, law, 600, 7, 3, 50);
    const std::vector<double> Rb = lowes_spectrum(both);
    double m2 = 0;
    int K2 = 0;
    for (int n = 300; n <= 600; n++) m2 += Rb[n] / law(n), K2++;
    s.check("unresolved part: estimate and realisation together average the law's spectrum over degrees 300–600", std::fabs(m2 / K2 - 1), 0.03);
  }

  // ---- truncation and sums ---------------------------------------------------------------------------------
  {
    const HarmonicPotential p = formula(40), t = p.truncated(10);
    const FieldValue a1 = field_at(p, a, 0.4, 1.0, 10), a2 = field_at(t, a, 0.4, 1.0);
    s.require("truncation: an expansion truncated to degree 10 is the first ten degrees of it", a1.north == a2.north && a1.down == a2.down);
    const HarmonicPotential sum = t + p;
    const FieldValue f = field_at(sum, a, 0.4, 1.0), g = field_at(p, a, 0.4, 1.0);
    s.check("sum: the field of a sum is the sum of the fields", std::fabs(f.down - (g.down + a2.down)) / std::fabs(g.down), 1e-12);
  }
  return s.done();
}
