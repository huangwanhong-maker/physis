// Chains on a classified set (decision 0020; the infrastructure's ADR 0067 and its design, §3.2): the states are
// classified at nested levels (a classification of divisions, groups, classes), with a law π on them, and at rate
// λ_ℓ a state is drawn again from π within its block at level ℓ. What the chain remembers of a state, level by level,
// then decays at rates that are sums of the coarser levels' rates.
//
// The classification. Levels ℓ = 0, …, L: block[ℓ][x] is the block of state x at level ℓ; level 0 is one block (the
// whole set), level L the states themselves, and each level refines the one before (two states in one block at level
// ℓ + 1 share their block at ℓ). The tree may be as uneven as a real classification is.
//
// The chain. Let Π_ℓ be the conditional expectation, under π, onto the functions of the level-ℓ block:
// Π_ℓ(x, y) = [B_ℓ(x) = B_ℓ(y)] π(y) / π(B_ℓ(x)), a Markov kernel that draws y from π within x's block. The
// generator is L = Σ_{ℓ<L} λ_ℓ (Π_ℓ − I). Each Π_ℓ keeps π and is self-adjoint in L²(π), so the chain is reversible
// with law π.
//
// Theorem (its modes). Let V_0 be the constants and, for m = 1, …, L, V_m the functions of the level-m block whose
// mean vanishes on every level-(m − 1) block. These are orthogonal in L²(π), they sum to all functions, and
//     L f = −μ_m f on V_m,   μ_m = λ_0 + … + λ_{m−1},
// so P_τ = e^{τL} = Σ_m e^{−μ_m τ} Q_m with Q_m = Π_m − Π_{m−1} (Q_0 = Π_0) the projection onto V_m; equivalently
//     P_τ(x, y) = Σ_m (e^{−μ_m τ} − e^{−μ_{m+1} τ}) Π_m(x, y),   e^{−μ_{L+1} τ} = 0,
// a mixture over the deepest level the chain still holds. Proof: for f ∈ V_m, Π_ℓ f = f when ℓ ≥ m (f is a function of
// the level-m block, so of every finer one) and Π_ℓ f = 0 when ℓ < m (f has mean zero on every level-(m−1) block, so on
// every coarser one); hence L f = −Σ_{ℓ<m} λ_ℓ f. ∎ So a level is remembered no longer than any coarser one
// (μ_1 ≤ … ≤ μ_L): a premises keeps its division at least as long as its class. Balanced trees with uniform π are the
// case of the design's Proposition 3.2, up to its parameters.
//
// What follows from the closed form: the kernel at any τ; the log-likelihood of observed transitions and its gradient
// in the rates; the rates' maximum likelihood (Newton's method in log λ, Laplace's errors); the weights with which each
// mode carries a transition (the design's foam, Theorem 4.3), e^{−μ_m τ} Q_m(x, y) / P_τ(x, y), signed in general; and
// the expected number of unseen changes between two observations, by the endpoint-conditioned expectation of Hobolth
// and Jensen (2005, Statistical Applications in Genetics and Molecular Biology 4, 18), whose integral the modes make
// closed: E[N | x, y, τ] = Σ_{a,b} I_ab(τ) [−μ_b δ_ab Q_b(x, y) − Σ_i Q_a(x, i) L_ii Q_b(i, y)] / P_τ(x, y), with
// I_ab(τ) = ∫_0^τ e^{−μ_a s − μ_b (τ − s)} ds. A redraw that returns the same state is not a change.
//
// Costs: the kernel O(L) a pair, the expected changes O(L² N) a pair for N states. Checks: tests/check_nested_chain.cpp
// (the closed form against calx's matrix exponential, the modes as eigenvectors, the expected changes against Van
// Loan's integral, the rates recovered from simulated transitions).

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include <calx/linalg/dense.hpp>

namespace physis::statistical_mechanics {

struct Classification {
  std::vector<std::vector<std::size_t>> block;  // block[ℓ][x], ℓ = 0 … L
  std::vector<double> law;                      // π, positive, summing to 1
  std::size_t states() const { return law.size(); }
  std::size_t levels() const { return block.empty() ? 0 : block.size() - 1; }  // L
};

// The classification checked: level 0 one block, level L the states, each level refining the one before, π a law.
inline void validate(const Classification& C) {
  const std::size_t N = C.states();
  if (C.block.size() < 2) throw std::invalid_argument("physis::statistical_mechanics::Classification: at least two levels (the whole set and the states)");
  for (const auto& b : C.block)
    if (b.size() != N) throw std::invalid_argument("physis::statistical_mechanics::Classification: a level does not classify every state");
  for (std::size_t x = 0; x < N; x++)
    if (C.block[0][x] != C.block[0][0]) throw std::invalid_argument("physis::statistical_mechanics::Classification: level 0 is not one block");
  std::vector<std::size_t> seen(N, 0);
  for (std::size_t x = 0; x < N; x++) {
    if (C.block.back()[x] >= N || seen[C.block.back()[x]]++) throw std::invalid_argument("physis::statistical_mechanics::Classification: the last level is not the states themselves");
  }
  for (std::size_t l = 0; l + 1 < C.block.size(); l++)
    for (std::size_t x = 0; x < N; x++)
      for (std::size_t y = x + 1; y < N; y++)
        if (C.block[l + 1][x] == C.block[l + 1][y] && C.block[l][x] != C.block[l][y])
          throw std::invalid_argument("physis::statistical_mechanics::Classification: level " + std::to_string(l + 1) + " does not refine level " + std::to_string(l));
  double s = 0;
  for (double p : C.law) {
    if (!(p > 0)) throw std::invalid_argument("physis::statistical_mechanics::Classification: the law is not positive");
    s += p;
  }
  if (std::fabs(s - 1) > 1e-9) throw std::invalid_argument("physis::statistical_mechanics::Classification: the law does not sum to 1");
}

// The chain of a classification with its rates λ_0 … λ_{L−1} (per unit of time): its blocks' masses and its memories.
struct NestedChain {
  Classification C;
  std::vector<double> lambda;                // λ_ℓ, ℓ < L
  std::vector<double> mu;                    // μ_m, m = 0 … L
  std::vector<std::vector<double>> mass;     // mass[ℓ][b] = π(block b at level ℓ)

  NestedChain(Classification c, std::vector<double> rates) : C(std::move(c)), lambda(std::move(rates)) {
    validate(C);
    const std::size_t L = C.levels();
    if (lambda.size() != L) throw std::invalid_argument("physis::statistical_mechanics::NestedChain: one rate for each level but the last");
    for (double r : lambda)
      if (!(r >= 0) || !std::isfinite(r)) throw std::invalid_argument("physis::statistical_mechanics::NestedChain: a rate is negative or not finite");
    mu.assign(L + 1, 0.0);
    for (std::size_t m = 1; m <= L; m++) mu[m] = mu[m - 1] + lambda[m - 1];
    mass.resize(L + 1);
    for (std::size_t l = 0; l <= L; l++) {
      std::size_t B = 0;
      for (auto b : C.block[l]) B = std::max(B, b + 1);
      mass[l].assign(B, 0.0);
      for (std::size_t x = 0; x < C.states(); x++) mass[l][C.block[l][x]] += C.law[x];
    }
  }
  std::size_t levels() const { return C.levels(); }
  // Π_m(x, y): y drawn from π within x's level-m block.
  double redraw(std::size_t m, std::size_t x, std::size_t y) const {
    return C.block[m][x] == C.block[m][y] ? C.law[y] / mass[m][C.block[m][x]] : 0.0;
  }
  // Q_m(x, y), the projection onto the m-th mode.
  double mode(std::size_t m, std::size_t x, std::size_t y) const { return m == 0 ? redraw(0, x, y) : redraw(m, x, y) - redraw(m - 1, x, y); }
  // P_τ(x, y).
  double kernel(std::size_t x, std::size_t y, double tau) const {
    const std::size_t L = levels();
    double p = 0;
    for (std::size_t m = 0; m <= L; m++) {
      const double here = std::exp(-mu[m] * tau), next = m == L ? 0.0 : std::exp(-mu[m + 1] * tau);
      p += (here - next) * redraw(m, x, y);
    }
    return p;
  }
  calx::linalg::Matrix kernel_matrix(double tau) const {
    calx::linalg::Matrix P(C.states());
    for (std::size_t x = 0; x < C.states(); x++)
      for (std::size_t y = 0; y < C.states(); y++) P(x, y) = kernel(x, y, tau);
    return P;
  }
  calx::linalg::Matrix generator() const {
    const std::size_t N = C.states();
    calx::linalg::Matrix G(N);
    for (std::size_t l = 0; l < levels(); l++)
      for (std::size_t x = 0; x < N; x++)
        for (std::size_t y = 0; y < N; y++) G(x, y) += lambda[l] * (redraw(l, x, y) - (x == y ? 1.0 : 0.0));
    return G;
  }
  // The weight with which each mode carries the transition x → y over τ: e^{−μ_m τ} Q_m(x, y) / P_τ(x, y), summing to 1.
  std::vector<double> mode_weights(std::size_t x, std::size_t y, double tau) const {
    const double p = kernel(x, y, tau);
    std::vector<double> w(levels() + 1);
    for (std::size_t m = 0; m <= levels(); m++) w[m] = std::exp(-mu[m] * tau) * mode(m, x, y) / p;
    return w;
  }
  // The expected number of changes (redraws to another state) between observing x and, τ later, y.
  double expected_changes(std::size_t x, std::size_t y, double tau) const {
    const std::size_t L = levels(), N = C.states();
    std::vector<double> diag(N, 0.0);  // L_ii
    for (std::size_t i = 0; i < N; i++)
      for (std::size_t l = 0; l < L; l++) diag[i] += lambda[l] * (redraw(l, i, i) - 1.0);
    const auto I = [&](std::size_t a, std::size_t b) {
      const double ma = mu[a], mb = mu[b];
      if (std::fabs(ma - mb) < 1e-12 * std::max(1.0, std::fabs(ma))) return tau * std::exp(-ma * tau);
      return (std::exp(-mb * tau) - std::exp(-ma * tau)) / (ma - mb);
    };
    double s = 0;
    for (std::size_t a = 0; a <= L; a++)
      for (std::size_t b = 0; b <= L; b++) {
        double c = a == b ? -mu[b] * mode(b, x, y) : 0.0;
        for (std::size_t i = 0; i < N; i++) {
          const double qa = mode(a, x, i);
          if (qa != 0) c -= qa * diag[i] * mode(b, i, y);
        }
        s += I(a, b) * c;
      }
    return s / kernel(x, y, tau);
  }
};

// An observed transition: state x, then y a time τ later, seen `count` times.
struct Transition {
  std::size_t x = 0, y = 0;
  double tau = 1;
  double count = 1;
};

// The log-likelihood of transitions under the chain, and its gradient in the rates.
inline double log_likelihood(const NestedChain& X, const std::vector<Transition>& T, std::vector<double>* gradient = nullptr) {
  const std::size_t L = X.levels();
  if (gradient) gradient->assign(L, 0.0);
  double s = 0;
  for (const auto& t : T) {
    const double p = X.kernel(t.x, t.y, t.tau);
    if (!(p > 0)) return -std::numeric_limits<double>::infinity();
    s += t.count * std::log(p);
    if (gradient) {
      // dP/dλ_ℓ = Σ_{m>ℓ} (−τ) e^{−μ_m τ} Q_m(x, y).
      std::vector<double> term(L + 1);
      for (std::size_t m = 0; m <= L; m++) term[m] = -t.tau * std::exp(-X.mu[m] * t.tau) * X.mode(m, t.x, t.y);
      double tail = 0;
      for (std::size_t l = L; l-- > 0;) {
        tail += term[l + 1];
        (*gradient)[l] += t.count * tail / p;
      }
    }
  }
  return s;
}

// The rates' maximum likelihood, with Laplace's errors and the memories' (by the delta method).
struct NestedFit {
  std::vector<double> lambda, lambda_error;  // per unit of time
  std::vector<double> mu, mu_error;          // μ_m, m = 0 … L
  std::vector<bool> at_bound;                // a rate the data drive to zero (its error then not given)
  std::vector<double> gradient;              // of the log-likelihood in log λ, at the maximum
  double log_likelihood = 0;
  int iterations = 0;
  bool converged = false;
};

namespace detail {

// The inverse of a small symmetric positive definite matrix, by Gauss–Jordan; empty if singular.
inline std::vector<std::vector<double>> inverse(std::vector<std::vector<double>> A) {
  const std::size_t n = A.size();
  std::vector<std::vector<double>> B(n, std::vector<double>(n, 0.0));
  for (std::size_t i = 0; i < n; i++) B[i][i] = 1;
  for (std::size_t c = 0; c < n; c++) {
    std::size_t p = c;
    for (std::size_t r = c + 1; r < n; r++)
      if (std::fabs(A[r][c]) > std::fabs(A[p][c])) p = r;
    if (std::fabs(A[p][c]) < 1e-300) return {};
    std::swap(A[p], A[c]), std::swap(B[p], B[c]);
    const double d = A[c][c];
    for (std::size_t k = 0; k < n; k++) A[c][k] /= d, B[c][k] /= d;
    for (std::size_t r = 0; r < n; r++)
      if (r != c) {
        const double f = A[r][c];
        for (std::size_t k = 0; k < n; k++) A[r][k] -= f * A[c][k], B[r][k] -= f * B[c][k];
      }
  }
  return B;
}

}  // namespace detail

// Newton's method in θ = log λ from `start` (positive rates), its Hessian by central differences of the exact
// gradient, steps halved until the likelihood rises; a rate that the data push to zero stops at `floor`, and one within a
// hundred times the floor that the likelihood does not push up is reported at its bound, with no error (the memories'
// errors then from the free rates alone).
inline NestedFit fit_rates(const Classification& C, const std::vector<Transition>& T, std::vector<double> start, int max_iterations = 100, double floor = 1e-8) {
  const std::size_t L = C.levels();
  if (start.size() != L) throw std::invalid_argument("physis::statistical_mechanics::fit_rates: one starting rate for each level but the last");
  std::vector<double> theta(L);
  for (std::size_t l = 0; l < L; l++) theta[l] = std::log(std::max(start[l], floor));
  const auto value = [&](const std::vector<double>& th, std::vector<double>* g) {
    std::vector<double> lam(L);
    for (std::size_t l = 0; l < L; l++) lam[l] = std::exp(th[l]);
    std::vector<double> gl;
    const double v = log_likelihood(NestedChain(C, lam), T, g ? &gl : nullptr);
    if (g) {
      g->assign(L, 0.0);
      for (std::size_t l = 0; l < L; l++) (*g)[l] = gl[l] * lam[l];
    }
    return v;
  };
  const auto hessian = [&](const std::vector<double>& th) {
    std::vector<std::vector<double>> H(L, std::vector<double>(L, 0.0));
    const double h = 1e-5;
    for (std::size_t j = 0; j < L; j++) {
      auto up = th, dn = th;
      up[j] += h, dn[j] -= h;
      std::vector<double> gu, gd;
      value(up, &gu), value(dn, &gd);
      for (std::size_t i = 0; i < L; i++) H[i][j] = (gu[i] - gd[i]) / (2 * h);
    }
    for (std::size_t i = 0; i < L; i++)
      for (std::size_t j = i + 1; j < L; j++) H[i][j] = H[j][i] = (H[i][j] + H[j][i]) / 2;
    return H;
  };
  // Converged when Newton's step over the free rates would promise to gain less than 1e-9 in the log-likelihood,
  // ½ g_F·(−H_F)⁻¹ g_F: a test in the likelihood's own units, whatever the rates' scale. A rate within a hundred times the
  // floor whose gradient does not point above it is held there (the bound's Karush–Kuhn–Tucker condition).
  const auto at_floor = [&](const std::vector<double>& th, const std::vector<double>& gr, std::size_t i) { return th[i] <= std::log(100 * floor) && gr[i] < 1e-5; };
  const auto promised = [&](const std::vector<double>& th, const std::vector<double>& gr, const std::vector<std::vector<double>>& H) {
    std::vector<std::size_t> free;
    for (std::size_t i = 0; i < L; i++)
      if (!at_floor(th, gr, i)) free.push_back(i);
    if (free.empty()) return 0.0;
    std::vector<std::vector<double>> A(free.size(), std::vector<double>(free.size()));
    for (std::size_t a = 0; a < free.size(); a++)
      for (std::size_t b = 0; b < free.size(); b++) A[a][b] = -H[free[a]][free[b]];
    const auto inv = detail::inverse(A);
    if (inv.empty()) return std::numeric_limits<double>::infinity();
    double d = 0;
    for (std::size_t a = 0; a < free.size(); a++)
      for (std::size_t b = 0; b < free.size(); b++) d += gr[free[a]] * inv[a][b] * gr[free[b]];
    return d < 0 ? std::numeric_limits<double>::infinity() : d / 2;
  };
  NestedFit F;
  std::vector<double> g;
  double v = value(theta, &g);
  for (F.iterations = 0; F.iterations < max_iterations; F.iterations++) {
    auto H = hessian(theta);
    if (promised(theta, g, H) < 1e-9) { F.converged = true; break; }
    std::vector<std::vector<double>> negH(L, std::vector<double>(L));
    for (std::size_t i = 0; i < L; i++)
      for (std::size_t j = 0; j < L; j++) negH[i][j] = -H[i][j];
    auto inv = detail::inverse(negH);
    std::vector<double> step(L, 0.0);
    bool newton = !inv.empty();
    if (newton) {
      for (std::size_t i = 0; i < L; i++)
        for (std::size_t j = 0; j < L; j++) step[i] += inv[i][j] * g[j];
      double dir = 0;
      for (std::size_t i = 0; i < L; i++) dir += step[i] * g[i];
      newton = dir > 0;  // an ascent direction
    }
    if (!newton) step = g;  // steepest ascent where the Hessian is not negative definite
    double scale = 1, best = v;
    std::vector<double> next = theta;
    for (int k = 0; k < 40; k++, scale /= 2) {
      std::vector<double> trial(L);
      for (std::size_t i = 0; i < L; i++) trial[i] = std::max(theta[i] + scale * std::clamp(step[i], -5.0, 5.0), std::log(floor));
      const double tv = value(trial, nullptr);
      if (tv > best) { best = tv, next = trial; break; }
    }
    double moved = 0;
    for (std::size_t i = 0; i < L; i++) moved = std::max(moved, std::fabs(next[i] - theta[i]));
    theta = next;
    const double gain = best - v;
    v = value(theta, &g);
    if (moved < 1e-10 && gain <= 1e-12) break;  // Newton stalls: the coordinates are polished below
  }
  // Where Newton stalls (a rate near its floor while another is free), each rate in turn is maximised alone by golden
  // section in θ between the floor and eight e-folds above its value, until a sweep gains nothing.
  for (int sweep = 0; !F.converged && sweep < 60; sweep++) {
    const double before = v;
    for (std::size_t l = 0; l < L; l++) {
      double a = std::log(floor), b = std::max(theta[l], std::log(floor)) + 8;
      const double phi = (std::sqrt(5.0) - 1) / 2;
      const auto at = [&](double t) { auto th = theta; th[l] = t; return value(th, nullptr); };
      double c = b - phi * (b - a), d = a + phi * (b - a), fc = at(c), fd = at(d);
      for (int k = 0; k < 80 && b - a > 1e-10; k++) {
        if (fc > fd) b = d, d = c, fd = fc, c = b - phi * (b - a), fc = at(c);
        else a = c, c = d, fc = fd, d = a + phi * (b - a), fd = at(d);
      }
      const double t = (a + b) / 2, ft = at(t);
      if (ft > v) theta[l] = t, v = ft;
    }
    v = value(theta, &g);
    F.iterations++;
    if (promised(theta, g, hessian(theta)) < 1e-9) F.converged = true;
    else if (v - before < 1e-12) break;
  }
  F.log_likelihood = v;
  F.gradient = g;
  F.lambda.resize(L), F.at_bound.assign(L, false);
  for (std::size_t l = 0; l < L; l++) F.lambda[l] = std::exp(theta[l]), F.at_bound[l] = at_floor(theta, g, l);
  // Laplace: the covariance of θ the inverse of −H; of λ, diag(λ) Cov(θ) diag(λ).
  auto H = hessian(theta);
  std::vector<std::vector<double>> negH(L, std::vector<double>(L));
  for (std::size_t i = 0; i < L; i++)
    for (std::size_t j = 0; j < L; j++) negH[i][j] = -H[i][j];
  const auto cov = detail::inverse(negH);
  F.lambda_error.assign(L, std::nan(""));
  F.mu.assign(L + 1, 0.0), F.mu_error.assign(L + 1, 0.0);
  for (std::size_t m = 1; m <= L; m++) F.mu[m] = F.mu[m - 1] + F.lambda[m - 1];
  if (!cov.empty()) {
    for (std::size_t l = 0; l < L; l++) F.lambda_error[l] = F.at_bound[l] ? std::nan("") : F.lambda[l] * std::sqrt(std::max(0.0, cov[l][l]));
    for (std::size_t m = 1; m <= L; m++) {
      double var = 0;
      for (std::size_t i = 0; i < m; i++)
        for (std::size_t j = 0; j < m; j++)
          if (!F.at_bound[i] && !F.at_bound[j]) var += F.lambda[i] * F.lambda[j] * cov[i][j];
      F.mu_error[m] = std::sqrt(std::max(0.0, var));
    }
  } else {
    for (std::size_t m = 1; m <= L; m++) F.mu_error[m] = std::nan("");
  }
  return F;
}

}  // namespace physis::statistical_mechanics
