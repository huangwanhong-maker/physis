// Random potential fields with a given spectrum, and the part of a field an estimate leaves
// unresolved (decision 0010).
//
// An isotropic Gaussian potential field on the sphere has independent coefficients gₙᵐ, hₙᵐ, each
// normal with mean 0 and a variance σₙ² that depends on the degree alone; its expected
// Lowes–Mauersberger spectrum at the reference sphere is then
//
//   E[Rₙ] = (n + 1)(2n + 1) σₙ²,
//
// since degree n has 2n + 1 coefficients (gₙ⁰ and the pairs gₙᵐ, hₙᵐ). A spectrum fixes the field's
// law; a seed fixes a realisation, the same on every machine (calx's counter-based streams).
//
// The law of a crustal field. Magnetised sources at a radius ρa below the reference sphere, random
// and uncorrelated, give a spectrum close to Rₙ = A nᵝ ρ²ⁿ: the power law of the source statistics
// times the attenuation of degree n over the depth (1 − ρ)a, which the field of every source at that
// radius suffers (the forms of Jackson 1994 and of Voorhies, Sabaka and Purucker 2002 are of this
// kind; β depends on the statistics of the sources). fit_spectral_law fits log A, β and log ρ by
// least squares on log Rₙ over a band of degrees, a linear problem.
//
// The unresolved part. An estimate of a field from data (an assimilation, a regularised inversion)
// is, under Gaussian assumptions, its conditional mean; the truth is the estimate plus an error
// independent of it, and by the law of total variance the error's power at degree n is the prior's
// less the estimate's: max(0, R_law(n) − R_estimate(n)), and the law's own beyond the estimate's
// degree. unresolved_part draws such an error: added to the estimate, it gives a field with the
// law's spectrum that agrees with the estimate wherever the estimate is certain, a realisation of
// the scales the data could not reach. It is synthetic, and says so.
//
// Caveat. The law is fitted where the estimate resolves the field; a band chosen too high takes the
// estimate's smoothing for the crust's spectrum. The fitted ρ is an effective depth of the sources
// and of the estimate's regularisation together, not a measured depth.

#pragma once

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <calx/linalg/dense.hpp>
#include <calx/random/stream.hpp>

#include "potential_field.hpp"

namespace physis::analysis {

// Rₙ = A nᵝ ρ²ⁿ at the reference sphere.
struct SpectralLaw {
  double amplitude = 0;  // A
  double exponent = 0;   // β
  double ratio = 1;      // ρ, the sources' radius over the reference sphere's
  double rms_log_residual = 0;  // of the fit, in natural log units
  int n_from = 0, n_to = 0;     // the band fitted
  double operator()(int n) const { return amplitude * std::pow(double(n), exponent) * std::pow(ratio, 2.0 * n); }
  // The depth of the sources below a sphere of radius a: (1 − ρ) a.
  double depth(double a) const { return (1 - ratio) * a; }
};

// The law fitted to spectrum[n] over n_from ≤ n ≤ n_to (positive values only).
inline SpectralLaw fit_spectral_law(const std::vector<double>& spectrum, int n_from, int n_to) {
  if (n_from < 1 || n_to >= int(spectrum.size()) || n_to - n_from < 3) throw std::invalid_argument("physis::analysis: a band of at least four degrees within the spectrum");
  calx::linalg::Matrix A(3), B(3);
  int used = 0;
  for (int n = n_from; n <= n_to; n++) {
    if (!(spectrum[n] > 0)) continue;
    const double x[3] = {1.0, std::log(double(n)), 2.0 * n};
    const double y = std::log(spectrum[n]);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) A(i, j) += x[i] * x[j];
      B(i, 0) += x[i] * y;
    }
    used++;
  }
  if (used < 4) throw std::invalid_argument("physis::analysis: too few positive degrees to fit");
  const calx::linalg::Matrix X = calx::linalg::solve(A, B);
  SpectralLaw law;
  law.amplitude = std::exp(X(0, 0));
  law.exponent = X(1, 0);
  law.ratio = std::exp(X(2, 0));
  law.n_from = n_from, law.n_to = n_to;
  double ss = 0;
  for (int n = n_from; n <= n_to; n++)
    if (spectrum[n] > 0) {
      const double r = std::log(spectrum[n]) - std::log(law(n));
      ss += r * r;
    }
  law.rms_log_residual = std::sqrt(ss / used);
  return law;
}

// A random expansion to degree n_to with E[Rₙ] = variance_power[n] at the reference sphere for
// n_from ≤ n ≤ n_to (and nothing below n_from), from the stream (seed, key).
inline HarmonicPotential random_potential(const std::vector<double>& spectrum, int n_from, int n_to, double radius,
                                          std::uint64_t seed, std::uint32_t key = 0) {
  if (n_from < 1 || n_to < n_from || n_to >= int(spectrum.size())) throw std::invalid_argument("physis::analysis: degrees outside the spectrum");
  HarmonicPotential p = HarmonicPotential::zero(n_to, radius);
  calx::random::Stream rng(seed, key, 0x5eed, 0);
  for (int n = n_from; n <= n_to; n++) {
    const double sd = std::sqrt(std::max(0.0, spectrum[n]) / ((n + 1.0) * (2.0 * n + 1.0)));
    for (int m = 0; m <= n; m++) {
      const std::size_t k = HarmonicPotential::index(n, m);
      p.g[k] = sd * rng.normal();
      p.h[k] = m == 0 ? 0.0 : sd * rng.normal();
    }
  }
  return p;
}

// The spectrum of the part an estimate leaves unresolved, to degree n_to: the law's power less the
// estimate's, never below zero, and the law's own beyond the estimate's degree.
inline std::vector<double> unresolved_spectrum(const HarmonicPotential& estimate, const SpectralLaw& law, int n_to) {
  const std::vector<double> R = lowes_spectrum(estimate);
  std::vector<double> u(n_to + 1, 0.0);
  for (int n = 1; n <= n_to; n++) u[n] = std::max(0.0, law(n) - (n <= estimate.degree ? R[n] : 0.0));
  return u;
}

// A realisation of that part, from n_from (the first degree the estimate does not fully resolve,
// by default 1) to n_to.
inline HarmonicPotential unresolved_part(const HarmonicPotential& estimate, const SpectralLaw& law, int n_to, std::uint64_t seed,
                                         std::uint32_t key = 0, int n_from = 1) {
  return random_potential(unresolved_spectrum(estimate, law, n_to), n_from, n_to, estimate.radius, seed, key);
}

}  // namespace physis::analysis
