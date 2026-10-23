// Estimates from Markov chains: means with their errors, autocorrelation times, jackknife.
//
// Successive configurations of a Markov chain Monte Carlo are correlated, so the error of a
// mean is not σ/√n. For a series x_1..x_n of an observable with variance σ² and normalised
// autocorrelation ρ(t), the variance of the mean is (2 τ_int σ²)/n, with the integrated
// autocorrelation time τ_int = 1/2 + Σ_{t>=1} ρ(t) (Madras and Sokal 1988, J. Stat. Phys.
// 50, 109–186, eqs. 2.15–2.19; Sokal 1997, Monte Carlo Methods in Statistical Mechanics,
// eqs. 3.10–3.19). The sum is cut at the self-consistent window W, the smallest W with
// W >= c τ_int(W), c = 6 by default (Sokal's automatic windowing), which keeps the bias
// O(e^{-c}) and the variance of the estimate of τ_int about 2(2W+1)/n τ_int².
//
// Binning (Flyvbjerg and Petersen 1989, J. Chem. Phys. 91, 461–466) gives the same error
// without a model of ρ: the series is cut into B bins whose means are nearly independent
// once a bin is long compared with τ_int, and the error is the standard error of the bin
// means. The jackknife (Efron 1982, ch. 3; Young 2015, arXiv:1210.3781, §3) carries the
// error to a nonlinear function of several means, such as the Binder cumulant: over B bins,
// f is evaluated on the B leave-one-bin-out means f_i, and the error is
// √((B-1)/B Σ (f_i - f̄)²).

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <span>
#include <stdexcept>
#include <vector>

namespace physis::field_theory {

struct Estimate {
  double mean = 0;
  double error = 0;  // one standard deviation of the mean
  double tau = 0.5;  // integrated autocorrelation time, in steps of the series
  std::size_t n = 0;
};

inline double mean(std::span<const double> x) {
  if (x.empty()) return 0;
  double s = 0;
  for (double v : x) s += v;
  return s / static_cast<double>(x.size());
}

inline double variance(std::span<const double> x) {
  if (x.size() < 2) return 0;
  const double m = mean(x);
  double s = 0;
  for (double v : x) s += (v - m) * (v - m);
  return s / static_cast<double>(x.size() - 1);
}

// The integrated autocorrelation time by Sokal's automatic window (c = 6).
inline double autocorrelation_time(std::span<const double> x, double c = 6) {
  const std::size_t n = x.size();
  if (n < 4) return 0.5;
  const double m = mean(x);
  double c0 = 0;
  for (double v : x) c0 += (v - m) * (v - m);
  c0 /= static_cast<double>(n);
  if (c0 <= 0) return 0.5;
  double tau = 0.5;
  for (std::size_t t = 1; t < n / 2; t++) {
    double ct = 0;
    for (std::size_t i = 0; i + t < n; i++) ct += (x[i] - m) * (x[i + t] - m);
    ct /= static_cast<double>(n - t);
    tau += ct / c0;
    if (static_cast<double>(t) >= c * tau) break;
  }
  return std::max(tau, 0.5);
}

// The mean with its error from the autocorrelation time: σ² 2 τ_int / n.
inline Estimate estimate(std::span<const double> x, double c = 6) {
  Estimate e;
  e.n = x.size();
  e.mean = mean(x);
  if (x.size() < 2) return e;
  e.tau = autocorrelation_time(x, c);
  double v = 0;
  for (double y : x) v += (y - e.mean) * (y - e.mean);
  v /= static_cast<double>(x.size());
  e.error = std::sqrt(2 * e.tau * v / static_cast<double>(x.size()));
  return e;
}

// The mean with the standard error of `bins` bin means (the last partial bin dropped).
inline Estimate binned(std::span<const double> x, std::size_t bins = 32) {
  Estimate e;
  e.n = x.size();
  e.mean = mean(x);
  bins = std::min(bins, x.size());
  if (bins < 2) return e;
  const std::size_t width = x.size() / bins;
  std::vector<double> means(bins);
  for (std::size_t b = 0; b < bins; b++) means[b] = mean(x.subspan(b * width, width));
  e.error = std::sqrt(variance(means) / static_cast<double>(bins));
  // The autocorrelation time implied by the binned variance: σ_bin² ≈ 2 τ σ² / width.
  const double v = variance(x);
  e.tau = v > 0 ? 0.5 * e.error * e.error * static_cast<double>(x.size()) / v : 0.5;
  return e;
}

// The jackknife of f(means of several series of equal length) over `bins` bins.
inline Estimate jackknife(const std::vector<std::vector<double>>& series, const std::function<double(const std::vector<double>&)>& f,
                          std::size_t bins = 32) {
  if (series.empty()) throw std::invalid_argument("physis::field_theory::jackknife: no series");
  const std::size_t n = series[0].size(), k = series.size();
  for (const auto& s : series)
    if (s.size() != n) throw std::invalid_argument("physis::field_theory::jackknife: series of unequal length");
  bins = std::min(bins, n);
  Estimate e;
  e.n = n;
  std::vector<double> total(k, 0);
  for (std::size_t j = 0; j < k; j++)
    for (double v : series[j]) total[j] += v;
  std::vector<double> all(k);
  for (std::size_t j = 0; j < k; j++) all[j] = total[j] / static_cast<double>(n);
  e.mean = f(all);
  if (bins < 2) return e;
  const std::size_t width = n / bins, used = width * bins;
  std::vector<double> fi(bins);
  for (std::size_t b = 0; b < bins; b++) {
    std::vector<double> m(k);
    for (std::size_t j = 0; j < k; j++) {
      double s = 0;
      for (std::size_t i = 0; i < used; i++)
        if (i / width != b) s += series[j][i];
      m[j] = s / static_cast<double>(used - width);
    }
    fi[b] = f(m);
  }
  const double fbar = mean(fi);
  double s = 0;
  for (double v : fi) s += (v - fbar) * (v - fbar);
  e.error = std::sqrt(static_cast<double>(bins - 1) / static_cast<double>(bins) * s);
  return e;
}

}  // namespace physis::field_theory
