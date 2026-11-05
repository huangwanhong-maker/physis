// Continuation of equilibria in a parameter, with their folds and Hopf points (decision 0022).
//
// The branch. Equilibria of x' = f(x, μ), x ∈ ℝⁿ, form curves in (x, μ). Pseudo-arclength continuation (Keller
// 1977, in Applications of Bifurcation Theory, Academic Press, 359–384; Kuznetsov 2004, Elements of Applied
// Bifurcation Theory, 3rd ed., §10.2) follows such a curve through its folds: from a point y = (x, μ) with unit
// tangent t, the predictor y + s t is corrected by Newton's method on the bordered system
//
//   f(y) = 0,   t·(y − y_predicted) = 0,
//
// whose Jacobian [f_x f_μ; tᵀ] is regular at regular points and at simple folds. The new tangent solves
// [f_x f_μ; t_oldᵀ] t = (0, 1), which keeps the orientation; the step grows after easy corrections and halves
// after hard ones. Derivatives are by central differences (flows.hpp's jacobian), or the caller's.
//
// Bifurcations along it. A fold (saddle-node) is where the tangent's μ-component changes sign. A Hopf point is
// where a pair of complex eigenvalues of f_x crosses the imaginary axis: the largest real part among the
// eigenvalues with imaginary part above a tolerance changes sign. Both are located by bisection on the arclength
// between the two points that bracket the change, each trial point corrected onto the curve, to the tolerance on
// the test function. A Hopf point records its frequency, the crossing pair's |Im λ|. Whether a Hopf bifurcation is
// super- or subcritical needs the first Lyapunov coefficient, which this module does not compute.
//
// Checks: tests/check_bifurcation.cpp — the fold of x' = μ − x² at μ = 0, the two folds of x' = μ + x − x³ at
// μ = ±2/(3√3), and the Brusselator's Hopf point b = 1 + a² with frequency a.
#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <calx/linalg/dense.hpp>
#include <calx/linalg/eigenvalues.hpp>

#include "flows.hpp"

namespace physis::dynamics {

struct ContinuationOptions {
  double step = 0.01;          // the initial arclength step
  double min_step = 1e-8, max_step = 0.1;
  double mu_min = -1e300, mu_max = 1e300;
  std::size_t max_points = 5000;
  int direction = +1;          // the initial sense in μ
  double tolerance = 1e-10;    // Newton's, on |f|
  double imaginary_tolerance = 1e-8;  // eigenvalues with |Im λ| above it are a complex pair
};

struct BranchPoint {
  State x;
  double mu = 0;
  bool stable = false;       // every eigenvalue of f_x has a negative real part
  double leading = 0;        // the largest real part of f_x's eigenvalues
};

struct Bifurcation {
  enum class Kind { fold, hopf } kind = Kind::fold;
  State x;
  double mu = 0;
  double frequency = 0;      // Hopf: |Im λ| of the crossing pair
};

struct Branch {
  std::vector<BranchPoint> points;
  std::vector<Bifurcation> bifurcations;
  std::string stopped;       // why the continuation ended
};

namespace continuation_detail {

// f as a field of y = (x, μ) ∈ ℝⁿ⁺¹ with values in ℝⁿ.
template <class F>
struct Augmented {
  F& f;
  std::size_t n;
  void eval(const State& y, State& out) const {
    const State x(y.begin(), y.begin() + static_cast<std::ptrdiff_t>(n));
    out.assign(n, 0.0);
    f(x, y[n], out);
  }
  // [f_x f_μ] by central differences, n × (n + 1), row-major.
  std::vector<double> jacobian(const State& y, double delta = 1e-7) const {
    std::vector<double> J(n * (n + 1));
    State yp(y), ym(y), fp, fm;
    for (std::size_t j = 0; j <= n; j++) {
      const double h = delta * (1 + std::fabs(y[j]));
      yp[j] = y[j] + h, ym[j] = y[j] - h;
      eval(yp, fp), eval(ym, fm);
      for (std::size_t i = 0; i < n; i++) J[i * (n + 1) + j] = (fp[i] - fm[i]) / (2 * h);
      yp[j] = ym[j] = y[j];
    }
    return J;
  }
};

// Solves the bordered (n + 1) × (n + 1) system [J; bᵀ] z = r.
inline State bordered_solve(const std::vector<double>& J, const State& b, const State& r, std::size_t n) {
  calx::linalg::Matrix A(n + 1), R(n + 1);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j <= n; j++) A(i, j) = J[i * (n + 1) + j];
  for (std::size_t j = 0; j <= n; j++) A(n, j) = b[j];
  for (std::size_t i = 0; i <= n; i++) R(i, 0) = r[i];
  const auto Z = calx::linalg::solve(A, R);
  State z(n + 1);
  for (std::size_t i = 0; i <= n; i++) z[i] = Z(i, 0);
  return z;
}

inline void normalise(State& t) {
  double s = 0;
  for (double v : t) s += v * v;
  s = std::sqrt(s);
  for (double& v : t) v /= s;
}

}  // namespace continuation_detail

// f(x, μ, out) writes f(x, μ). Starts from an approximate equilibrium x0 at μ0, corrected first at fixed μ0.
template <class F>
Branch continue_equilibria(F&& f, State x0, double mu0, const ContinuationOptions& o = {}) {
  using namespace continuation_detail;
  const std::size_t n = x0.size();
  Augmented<F> A{f, n};
  Branch out;
  // Newton at fixed μ0.
  {
    State y(x0);
    y.push_back(mu0);
    for (int it = 0; it < 50; it++) {
      State fy;
      A.eval(y, fy);
      double norm = 0;
      for (double v : fy) norm = std::max(norm, std::fabs(v));
      if (norm < o.tolerance) break;
      State b(n + 1, 0.0), r(n + 1, 0.0);
      b[n] = 1;
      for (std::size_t i = 0; i < n; i++) r[i] = -fy[i];
      const State d = bordered_solve(A.jacobian(y), b, r, n);
      for (std::size_t i = 0; i < n; i++) y[i] += d[i];
      if (it == 49) throw std::runtime_error("physis::dynamics::continue_equilibria: no equilibrium near the start");
    }
    x0.assign(y.begin(), y.begin() + static_cast<std::ptrdiff_t>(n));
  }
  State y(x0);
  y.push_back(mu0);
  // The initial tangent: [f_x f_μ; e_μᵀ] t = (0, 1), oriented by `direction`.
  State t;
  {
    State b(n + 1, 0.0), r(n + 1, 0.0);
    b[n] = 1, r[n] = 1;
    t = bordered_solve(A.jacobian(y), b, r, n);
    normalise(t);
    if ((t[n] > 0 ? 1 : -1) != o.direction) for (double& v : t) v = -v;
  }
  auto describe = [&](const State& yy) {
    BranchPoint p;
    p.x.assign(yy.begin(), yy.begin() + static_cast<std::ptrdiff_t>(n));
    p.mu = yy[n];
    const auto J = A.jacobian(yy);
    calx::linalg::Matrix fx(n);
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++) fx(i, j) = J[i * (n + 1) + j];
    double lead = -1e300, hopf = -1e300, freq = 0;
    for (const auto& e : calx::linalg::eigenvalues(fx)) {
      lead = std::max(lead, e.real());
      if (std::fabs(e.imag()) > o.imaginary_tolerance && e.real() > hopf) hopf = e.real(), freq = std::fabs(e.imag());
    }
    p.leading = lead;
    p.stable = lead < 0;
    return std::tuple{p, hopf, freq};
  };
  // Corrects a prediction onto the curve; returns nothing if Newton fails.
  auto correct = [&](const State& predicted, const State& tangent) -> std::optional<State> {
    State z(predicted);
    for (int it = 0; it < 15; it++) {
      State fz;
      A.eval(z, fz);
      double norm = 0;
      for (double v : fz) norm = std::max(norm, std::fabs(v));
      double along = 0;
      for (std::size_t j = 0; j <= n; j++) along += tangent[j] * (z[j] - predicted[j]);
      if (norm < o.tolerance && std::fabs(along) < o.tolerance) return z;
      State r(n + 1);
      for (std::size_t i = 0; i < n; i++) r[i] = -fz[i];
      r[n] = -along;
      const State d = bordered_solve(A.jacobian(z), tangent, r, n);
      for (std::size_t j = 0; j <= n; j++) z[j] += d[j];
      for (double v : z)
        if (!std::isfinite(v)) return std::nullopt;
    }
    return std::nullopt;
  };
  auto tangent_at = [&](const State& yy, const State& old) {
    State r(n + 1, 0.0);
    r[n] = 1;
    State tt = bordered_solve(A.jacobian(yy), old, r, n);
    normalise(tt);
    return tt;
  };
  auto [p0, h0, f0] = describe(y);
  out.points.push_back(p0);
  double hopf_prev = h0;
  double s = o.step;
  while (out.points.size() < o.max_points) {
    State pred(n + 1);
    for (std::size_t j = 0; j <= n; j++) pred[j] = y[j] + s * t[j];
    auto next = correct(pred, t);
    if (!next) {
      s /= 2;
      if (s < o.min_step) { out.stopped = "the corrector failed at the least step"; break; }
      continue;
    }
    const State tn = tangent_at(*next, t);
    auto [p, hopf, freq] = describe(*next);
    // A fold: the tangent's μ-component changes sign. Bisect on the arclength from y.
    auto bisect = [&](auto test, double lo, double hi) {
      State best = *next;
      for (int it = 0; it < 60 && hi - lo > 1e-14; it++) {
        const double mid = 0.5 * (lo + hi);
        State pm(n + 1);
        for (std::size_t j = 0; j <= n; j++) pm[j] = y[j] + mid * t[j];
        auto z = correct(pm, t);
        if (!z) break;
        best = *z;
        if (test(*z) < 0) hi = mid;  // the test's sign flips between y (≥ 0) and next (< 0), normalised by the caller
        else lo = mid;
      }
      return best;
    };
    if ((tn[n] > 0) != (t[n] > 0)) {
      const double sign0 = t[n] > 0 ? 1 : -1;
      const State at = bisect([&](const State& z) { return sign0 * tangent_at(z, t)[n]; }, 0, s);
      Bifurcation b;
      b.kind = Bifurcation::Kind::fold;
      b.x.assign(at.begin(), at.begin() + static_cast<std::ptrdiff_t>(n));
      b.mu = at[n];
      out.bifurcations.push_back(b);
    }
    if (hopf_prev > -1e299 && hopf > -1e299 && (hopf_prev < 0) != (hopf < 0)) {
      const double sign0 = hopf_prev < 0 ? -1 : 1;
      const State at = bisect([&](const State& z) { return sign0 * std::get<1>(describe(z)); }, 0, s);
      Bifurcation b;
      b.kind = Bifurcation::Kind::hopf;
      b.x.assign(at.begin(), at.begin() + static_cast<std::ptrdiff_t>(n));
      b.mu = at[n];
      b.frequency = std::get<2>(describe(at));
      out.bifurcations.push_back(b);
    }
    hopf_prev = hopf;
    y = *next, t = tn;
    out.points.push_back(p);
    if (y[n] < o.mu_min || y[n] > o.mu_max) { out.stopped = "the parameter left its range"; break; }
    s = std::min(o.max_step, s * 1.3);
  }
  if (out.stopped.empty()) out.stopped = "the most points";
  return out;
}

}  // namespace physis::dynamics
