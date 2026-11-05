// Flows of vector fields of any dimension, their linearisation, and the Lyapunov spectra of flows and
// maps (decision 0022).
//
// A vector field is a callable f(t, x, dx) writing dx = f(t, x) for x ∈ ℝⁿ, n known only at run time
// (calx's integrators are for a dimension fixed at compile time). rk4_step is the classical fourth-order
// Runge–Kutta step; integrate() steps from t0 to t1 with a step h, shortening the last so as to land on
// t1, and hands each state to an observer.
//
// Linearisation. jacobian_vector(f, t, x, v) is the directional derivative J(t, x) v by central
// differences, [f(x + εv) − f(x − εv)] / 2ε with ε = δ (1 + |x|) / |v|: its error is O(ε²) from the
// truncation and O(u / ε) from rounding (u the unit round-off), about 1e−10 relative at δ = 1e−6.
// A caller with an exact derivative passes its own instead (lyapunov_spectrum's second form).
//
// Lyapunov spectra (Benettin, Galgani, Giorgilli and Strelcyn 1980, Meccanica 15, 9–20 and 21–30; the
// "standard method", Parker and Chua 1989, Practical Numerical Algorithms for Chaotic Systems, §6). k
// tangent vectors v_j are carried with the state, v' = J(x) v for a flow, v ↦ DF(x) v for a map, and
// every m steps orthonormalised by the modified Gram–Schmidt process, which writes the propagated frame
// as Q R; the logarithms of R's diagonal accumulate, and their means over the time are the exponents
// λ_1 ≥ … ≥ λ_k. With k = n their sum is the mean of the divergence tr J along the path, which for a
// field of constant divergence (Lorenz's, −(σ + 1 + β)) is exact whatever the path: the check of the
// method. The finite-time exponents converge as 1/√T for chaotic paths; the spectrum depends on the
// initial state only through the transient.
//
// Kaplan and Yorke's dimension (1979, Lecture Notes in Mathematics 730, 204–227): with j the largest
// index such that λ_1 + … + λ_j ≥ 0, D_KY = j + (λ_1 + … + λ_j) / |λ_{j+1}|.
//
// Checks: tests/check_chaos.cpp.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <calx/linalg/dense.hpp>

namespace physis::dynamics {

using State = std::vector<double>;

namespace detail {
inline double norm(const State& v) {
  double s = 0;
  for (double x : v) s += x * x;
  return std::sqrt(s);
}
}  // namespace detail

// One step of the classical Runge–Kutta method for x' = f(t, x).
template <class F>
void rk4_step(F&& f, double t, State& x, double h) {
  const std::size_t n = x.size();
  State k1(n), k2(n), k3(n), k4(n), y(n);
  f(t, x, k1);
  for (std::size_t i = 0; i < n; i++) y[i] = x[i] + 0.5 * h * k1[i];
  f(t + 0.5 * h, y, k2);
  for (std::size_t i = 0; i < n; i++) y[i] = x[i] + 0.5 * h * k2[i];
  f(t + 0.5 * h, y, k3);
  for (std::size_t i = 0; i < n; i++) y[i] = x[i] + h * k3[i];
  f(t + h, y, k4);
  for (std::size_t i = 0; i < n; i++) x[i] += h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
}

// Steps from t0 to t1 (h > 0), the last step shortened to land on t1; observe(t, x) after each step.
// Returns the number of steps.
template <class F, class Observer>
std::size_t integrate(F&& f, double t0, double t1, State& x, double h, Observer&& observe) {
  if (!(h > 0)) throw std::invalid_argument("physis::dynamics::integrate: a positive step");
  std::size_t steps = 0;
  double t = t0;
  while (t < t1) {
    const double s = std::min(h, t1 - t);
    if (s <= 1e-15 * std::max(1.0, std::fabs(t1))) break;
    rk4_step(f, t, x, s);
    t = (t1 - t - s <= 1e-15 * std::max(1.0, std::fabs(t1))) ? t1 : t + s;
    steps++;
    observe(t, static_cast<const State&>(x));
  }
  return steps;
}

template <class F>
std::size_t integrate(F&& f, double t0, double t1, State& x, double h) {
  return integrate(f, t0, t1, x, h, [](double, const State&) {});
}

// J(t, x) v by central differences.
template <class F>
State jacobian_vector(F&& f, double t, const State& x, const State& v, double delta = 1e-6) {
  const std::size_t n = x.size();
  const double nv = detail::norm(v);
  State out(n, 0.0);
  if (nv == 0) return out;
  const double eps = delta * (1 + detail::norm(x)) / nv;
  State xp(x), xm(x), fp(n), fm(n);
  for (std::size_t i = 0; i < n; i++) xp[i] += eps * v[i], xm[i] -= eps * v[i];
  f(t, xp, fp);
  f(t, xm, fm);
  for (std::size_t i = 0; i < n; i++) out[i] = (fp[i] - fm[i]) / (2 * eps);
  return out;
}

// The Jacobian matrix, column by column.
template <class F>
calx::linalg::Matrix jacobian(F&& f, double t, const State& x, double delta = 1e-6) {
  const std::size_t n = x.size();
  calx::linalg::Matrix J(n);
  State e(n, 0.0);
  for (std::size_t j = 0; j < n; j++) {
    e[j] = 1;
    const State c = jacobian_vector(f, t, x, e, delta);
    for (std::size_t i = 0; i < n; i++) J(i, j) = c[i];
    e[j] = 0;
  }
  return J;
}

struct LyapunovOptions {
  std::size_t count = 0;              // exponents wanted, k ≤ n; 0 for all n
  double step = 0.01;                 // h of the flow; ignored for maps
  double transient = 0;               // time (flows) or iterations (maps) discarded first
  double time = 1000;                 // time (flows) or iterations (maps) averaged over
  std::size_t renormalise_every = 10; // steps between orthonormalisations
  std::size_t record_every = 0;       // renormalisations between recorded running estimates; 0: none
};

struct LyapunovSpectrum {
  std::vector<double> exponents;             // λ_1 ≥ … ≥ λ_k, per unit time (flows) or per iteration (maps)
  double time = 0;                           // the time averaged over
  std::vector<double> times;                 // where running estimates were recorded
  std::vector<std::vector<double>> running;  // the running estimates there
  State final_state;
};

namespace detail {
// Modified Gram–Schmidt on the columns v[0..k); returns the norms (the diagonal of R).
inline std::vector<double> orthonormalise(std::vector<State>& v) {
  std::vector<double> r(v.size());
  for (std::size_t j = 0; j < v.size(); j++) {
    for (std::size_t i = 0; i < j; i++) {
      double d = 0;
      for (std::size_t c = 0; c < v[j].size(); c++) d += v[i][c] * v[j][c];
      for (std::size_t c = 0; c < v[j].size(); c++) v[j][c] -= d * v[i][c];
    }
    r[j] = norm(v[j]);
    if (!(r[j] > 0)) throw std::runtime_error("physis::dynamics::lyapunov_spectrum: a tangent vector collapsed");
    for (double& c : v[j]) c /= r[j];
  }
  return r;
}

inline std::size_t count_of(const LyapunovOptions& o, std::size_t n) {
  const std::size_t k = o.count == 0 ? n : o.count;
  if (k > n) throw std::invalid_argument("physis::dynamics::lyapunov_spectrum: at most n exponents");
  return k;
}

inline std::vector<State> initial_frame(std::size_t n, std::size_t k) {
  std::vector<State> v(k, State(n, 0.0));
  for (std::size_t j = 0; j < k; j++) v[j][j] = 1;
  return v;
}
}  // namespace detail

// The spectrum of the flow of f, with the tangent dynamics given by jv(t, x, v) = J(t, x) v.
template <class F, class JV>
LyapunovSpectrum lyapunov_spectrum(F&& f, JV&& jv, State x, const LyapunovOptions& o) {
  const std::size_t n = x.size(), k = detail::count_of(o, n);
  const double h = o.step;
  if (!(h > 0) || !(o.time > 0) || o.renormalise_every == 0) throw std::invalid_argument("physis::dynamics::lyapunov_spectrum: positive step, time and interval");
  if (o.transient > 0) integrate(f, 0.0, o.transient, x, h);
  // The augmented system: the state, then the k tangent vectors.
  auto augmented = [&](double t, const State& y, State& dy) {
    const State xs(y.begin(), y.begin() + static_cast<std::ptrdiff_t>(n));
    State fx(n);
    f(t, xs, fx);
    std::copy(fx.begin(), fx.end(), dy.begin());
    for (std::size_t j = 0; j < k; j++) {
      const State vj(y.begin() + static_cast<std::ptrdiff_t>(n * (j + 1)), y.begin() + static_cast<std::ptrdiff_t>(n * (j + 2)));
      const State w = jv(t, xs, vj);
      std::copy(w.begin(), w.end(), dy.begin() + static_cast<std::ptrdiff_t>(n * (j + 1)));
    }
  };
  std::vector<State> v = detail::initial_frame(n, k);
  State y(n * (k + 1));
  std::copy(x.begin(), x.end(), y.begin());
  std::vector<double> sums(k, 0.0);
  LyapunovSpectrum out;
  const std::size_t steps = static_cast<std::size_t>(std::llround(o.time / h));
  double t = o.transient;
  std::size_t renormalisations = 0;
  for (std::size_t s = 1; s <= steps; s++) {
    for (std::size_t j = 0; j < k; j++) std::copy(v[j].begin(), v[j].end(), y.begin() + static_cast<std::ptrdiff_t>(n * (j + 1)));
    rk4_step(augmented, t, y, h);
    t += h;
    for (std::size_t j = 0; j < k; j++) std::copy(y.begin() + static_cast<std::ptrdiff_t>(n * (j + 1)), y.begin() + static_cast<std::ptrdiff_t>(n * (j + 2)), v[j].begin());
    if (s % o.renormalise_every == 0 || s == steps) {
      const auto r = detail::orthonormalise(v);
      for (std::size_t j = 0; j < k; j++) sums[j] += std::log(r[j]);
      renormalisations++;
      if (o.record_every && renormalisations % o.record_every == 0) {
        const double elapsed = static_cast<double>(s) * h;
        std::vector<double> est(k);
        for (std::size_t j = 0; j < k; j++) est[j] = sums[j] / elapsed;
        out.times.push_back(elapsed);
        out.running.push_back(std::move(est));
      }
    }
  }
  out.time = static_cast<double>(steps) * h;
  out.exponents.resize(k);
  for (std::size_t j = 0; j < k; j++) out.exponents[j] = sums[j] / out.time;
  out.final_state.assign(y.begin(), y.begin() + static_cast<std::ptrdiff_t>(n));
  return out;
}

// The spectrum of the flow of f, its tangent dynamics by central differences.
template <class F>
LyapunovSpectrum lyapunov_spectrum(F&& f, State x, const LyapunovOptions& o) {
  return lyapunov_spectrum(f, [&](double t, const State& xs, const State& v) { return jacobian_vector(f, t, xs, v); }, std::move(x), o);
}

// The spectrum of the map x ↦ F(x), F(x, out); o.transient and o.time count iterations.
template <class Map>
LyapunovSpectrum lyapunov_spectrum_map(Map&& F, State x, const LyapunovOptions& o) {
  const std::size_t n = x.size(), k = detail::count_of(o, n);
  auto field = [&](double, const State& xs, State& out) { F(xs, out); };
  State next(n);
  for (std::size_t s = 0; s < static_cast<std::size_t>(o.transient); s++) F(x, next), x.swap(next);
  std::vector<State> v = detail::initial_frame(n, k);
  std::vector<double> sums(k, 0.0);
  const std::size_t steps = static_cast<std::size_t>(o.time);
  for (std::size_t s = 1; s <= steps; s++) {
    for (std::size_t j = 0; j < k; j++) v[j] = jacobian_vector(field, 0.0, x, v[j]);
    F(x, next);
    x.swap(next);
    if (s % o.renormalise_every == 0 || s == steps) {
      const auto r = detail::orthonormalise(v);
      for (std::size_t j = 0; j < k; j++) sums[j] += std::log(r[j]);
    }
  }
  LyapunovSpectrum out;
  out.time = static_cast<double>(steps);
  out.exponents.resize(k);
  for (std::size_t j = 0; j < k; j++) out.exponents[j] = sums[j] / out.time;
  out.final_state = x;
  return out;
}

// Kaplan and Yorke's dimension of a spectrum sorted in decreasing order.
inline double kaplan_yorke(const std::vector<double>& lambda) {
  double sum = 0;
  for (std::size_t j = 0; j < lambda.size(); j++) {
    if (sum + lambda[j] < 0) return j == 0 ? 0.0 : static_cast<double>(j) + sum / std::fabs(lambda[j]);
    sum += lambda[j];
  }
  return static_cast<double>(lambda.size());
}

// ---- test systems with known invariants ----------------------------------------------------------------------------------

// Lorenz (1963): x' = σ(y − x), y' = x(r − z) − y, z' = xy − bz; divergence −(σ + 1 + b).
struct Lorenz {
  double sigma = 10, r = 28, b = 8.0 / 3;
  void operator()(double, const State& x, State& dx) const {
    dx[0] = sigma * (x[1] - x[0]);
    dx[1] = x[0] * (r - x[2]) - x[1];
    dx[2] = x[0] * x[1] - b * x[2];
  }
  State jacobian_vector(const State& x, const State& v) const {
    return {sigma * (v[1] - v[0]), (r - x[2]) * v[0] - v[1] - x[0] * v[2], x[1] * v[0] + x[0] * v[1] - b * v[2]};
  }
};

// Lorenz (1996): x_i' = (x_{i+1} − x_{i−2}) x_{i−1} − x_i + F, cyclic; divergence −n.
struct Lorenz96 {
  double forcing = 8;
  void operator()(double, const State& x, State& dx) const {
    const std::size_t n = x.size();
    for (std::size_t i = 0; i < n; i++) dx[i] = (x[(i + 1) % n] - x[(i + n - 2) % n]) * x[(i + n - 1) % n] - x[i] + forcing;
  }
};

// Rössler (1976): x' = −y − z, y' = x + ay, z' = b + z(x − c).
struct Rossler {
  double a = 0.2, b = 0.2, c = 5.7;
  void operator()(double, const State& x, State& dx) const {
    dx[0] = -x[1] - x[2];
    dx[1] = x[0] + a * x[1];
    dx[2] = b + x[2] * (x[0] - c);
  }
};

// Hénon (1976): (x, y) ↦ (1 − a x² + y, b x); its Jacobian's determinant is −b everywhere.
struct Henon {
  double a = 1.4, b = 0.3;
  void operator()(const State& x, State& out) const {
    out[0] = 1 - a * x[0] * x[0] + x[1];
    out[1] = b * x[0];
  }
};

// Rosenzweig and MacArthur's populations (1963; Rosenzweig 1971, Science 171, 385): a prey of logistic growth and a
// predator of saturating (Holling II) appetite,
//   x' = r x (1 − x/K) − a x y / (1 + a h x),   y' = e a x y / (1 + a h x) − m y.
// The coexistence equilibrium x* = m / (a (e − m h)) loses stability at a Hopf point where the prey's isocline
// peaks, x* = (K − 1/(a h))/2, so at K_H = 2x* + 1/(a h): enriching the prey destabilises the system (the paradox
// of enrichment). The form is general: interacting molecules, cells and species alike.
struct RosenzweigMacArthur {
  double r = 1, K = 2, a = 1, h = 1, e = 0.5, m = 0.2;
  void operator()(double, const State& x, State& dx) const {
    const double eat = a * x[0] * x[1] / (1 + a * h * x[0]);
    dx[0] = r * x[0] * (1 - x[0] / K) - eat;
    dx[1] = e * eat - m * x[1];
  }
  double prey_at_coexistence() const { return m / (a * (e - m * h)); }
  double predator_at_coexistence() const {
    const double x = prey_at_coexistence();
    return r * (1 - x / K) * (1 + a * h * x) / a;
  }
  double hopf_capacity() const { return 2 * prey_at_coexistence() + 1 / (a * h); }
};

// The logistic map x ↦ r x (1 − x); at r = 4 its exponent is ln 2.
struct Logistic {
  double r = 4;
  void operator()(const State& x, State& out) const { out[0] = r * x[0] * (1 - x[0]); }
};

}  // namespace physis::dynamics
