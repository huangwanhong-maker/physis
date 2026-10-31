// Checks of first passage (first_passage.hpp) against closed forms and Monte Carlo.
//
// - One transient state left at rate λ: the exponential law's density, survival, moments k!/λᵏ and
//   quantiles −ln(1 − p)/λ.
// - A chain of m transient states passed at rate λ: Erlang's density λᵐτ^{m−1}e^{−λτ}/(m−1)!, its survival
//   Σ_{j<m} e^{−λτ}(λτ)ʲ/j!, its moments (m+k−1)!/((m−1)! λᵏ), and quantiles inverting the survival.
// - An atom at zero when the initial law sums to less than 1.
// - A random chain of four transient and two absorbing states, against 40000 paths sampled with calx's
//   streams: the mean time, the probability of being absorbed by the median, and the law of the
//   absorbing state, each within four standard errors.
// - Gambler's ruin with rates λ up and μ down: absorption at N from i with probability
//   (1 − rⁱ)/(1 − rᴺ), r = μ/λ.
// - Mean first passage on a symmetric ring of n states at rate k each way: m_i = i(n − i)/(2k); Kac's
//   identity Σ_j k_ij m_j = 1/π_i − 1 (the mean return time 1/(π_i q_i)) on random chains for every
//   target; and the phase-type mean equal to the mean first passage.
// - A transient state that never leaves is refused.

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <physis/statistical_mechanics/first_passage.hpp>

using calx::linalg::Matrix;
using calx::verification::Suite;
namespace sm = physis::statistical_mechanics;

namespace {

double factorial(int n) { return std::tgamma(n + 1.0); }

}  // namespace

int main() {
  Suite s;
  calx::random::Stream rng(20261030, 19, 0, 0);

  // ---- the exponential and Erlang laws -------------------------------------------------------------------------------------------
  {
    const double lambda = 2.5;
    Matrix k(2);
    k(0, 1) = lambda;
    const std::vector<std::size_t> transient{0};
    const std::vector<double> alpha{1.0};
    const auto P = sm::phase_type(k, transient, alpha);
    double worst = 0;
    for (double tau : {0.0, 0.1, 0.5, 1.3, 4.0}) {
      worst = std::max(worst, std::fabs(P.density(tau) / (lambda * std::exp(-lambda * tau)) - 1));
      worst = std::max(worst, std::fabs(P.survival(tau) / std::exp(-lambda * tau) - 1));
    }
    for (int m = 1; m <= 4; m++) worst = std::max(worst, std::fabs(P.moment(m) / (factorial(m) / std::pow(lambda, m)) - 1));
    for (double p : {0.01, 0.5, 0.9, 0.999}) worst = std::max(worst, std::fabs(P.quantile(p) / (-std::log(1 - p) / lambda) - 1));
    s.check("exponential: density, survival, moments k!/λᵏ and quantiles −ln(1 − p)/λ (relative)", worst, 1e-12);

    const int m = 4;
    const double l = 1.5;
    Matrix e(m + 1);
    for (int i = 0; i < m; i++) e(i, i + 1) = l;
    const std::vector<std::size_t> tr{0, 1, 2, 3};
    const std::vector<double> a{1, 0, 0, 0};
    const auto Er = sm::phase_type(e, tr, a);
    double dens = 0, surv = 0, mom = 0, quant = 0;
    for (double tau : {0.2, 1.0, 2.7, 6.0}) {
      dens = std::max(dens, std::fabs(Er.density(tau) / (std::pow(l, m) * std::pow(tau, m - 1) * std::exp(-l * tau) / factorial(m - 1)) - 1));
      double sv = 0;
      for (int j = 0; j < m; j++) sv += std::exp(-l * tau) * std::pow(l * tau, j) / factorial(j);
      surv = std::max(surv, std::fabs(Er.survival(tau) / sv - 1));
    }
    for (int r = 1; r <= 3; r++) mom = std::max(mom, std::fabs(Er.moment(r) / (factorial(m + r - 1) / (factorial(m - 1) * std::pow(l, r))) - 1));
    for (double p : {0.05, 0.5, 0.95}) quant = std::max(quant, std::fabs(Er.survival(Er.quantile(p)) - (1 - p)));
    s.check("Erlang(4, λ): the density λᵐτ^{m−1}e^{−λτ}/(m−1)! (relative)", dens, 1e-12);
    s.check("Erlang(4, λ): the survival Σ_{j<m} e^{−λτ}(λτ)ʲ/j! (relative)", surv, 1e-12);
    s.check("Erlang(4, λ): the moments (m+k−1)!/((m−1)! λᵏ) (relative)", mom, 1e-12);
    s.check("Erlang(4, λ): S(quantile(p)) = 1 − p", quant, 1e-12);

    const std::vector<double> part{0.8};
    const auto At = sm::phase_type(k, transient, part);
    s.require("an initial law of mass 0.8: survival(0) = 0.8, an atom 0.2 at zero, quantile(0.1) = 0, E τ⁰ = 1",
              std::fabs(At.survival(0) - 0.8) < 1e-15 && std::fabs(At.atom() - 0.2) < 1e-15 && At.quantile(0.1) == 0 && At.moment(0) == 1 &&
                  std::fabs(At.moment(1) - 0.8 / lambda) < 1e-15);
  }

  // ---- Monte Carlo ----------------------------------------------------------------------------------------------------------------------
  {
    const std::size_t n = 6;
    Matrix k(n);
    for (std::size_t i = 0; i < 4; i++)
      for (std::size_t j = 0; j < n; j++)
        if (i != j && rng.uniform() < 0.7) k(i, j) = rng.uniform(0.2, 2.0);
    for (std::size_t i = 0; i < 3; i++) k(i, i + 1) = std::max(k(i, i + 1), 0.3);
    k(3, 4) = std::max(k(3, 4), 0.5);  // a path from every transient state to an exit
    const std::vector<std::size_t> transient{0, 1, 2, 3}, absorbing{4, 5};
    const std::vector<double> alpha{0.4, 0.3, 0.2, 0.1};
    const auto P = sm::phase_type(k, transient, alpha);
    const auto B = sm::absorption_probabilities(k, transient, absorbing);
    const double mean = P.moment(1), sd = std::sqrt(P.moment(2) - mean * mean), median = P.quantile(0.5);
    const int paths = 40000;
    double sum = 0;
    int by_median = 0, at4 = 0;
    for (int path = 0; path < paths; path++) {
      double u = rng.uniform(), t = 0;
      std::size_t x = 0;
      while (x < 3 && u >= alpha[x]) u -= alpha[x++];
      while (x < 4) {
        double q = 0;
        for (std::size_t j = 0; j < n; j++)
          if (j != x) q += k(x, j);
        t += rng.exponential() / q;
        double v = rng.uniform() * q;
        std::size_t y = 0;
        for (; y < n; y++) {
          if (y == x) continue;
          if (v < k(x, y)) break;
          v -= k(x, y);
        }
        x = std::min(y, n - 1);
      }
      sum += t;
      by_median += t <= median;
      at4 += x == 4;
    }
    double expected4 = 0;
    for (std::size_t a = 0; a < 4; a++) expected4 += alpha[a] * B(a, 0);
    const double z_mean = (sum / paths - mean) / (sd / std::sqrt(paths));
    const double z_median = (static_cast<double>(by_median) / paths - 0.5) / std::sqrt(0.25 / paths);
    const double z_law = (static_cast<double>(at4) / paths - expected4) / std::sqrt(expected4 * (1 - expected4) / paths);
    s.check("Monte Carlo, 40000 paths: the mean time α(−T)⁻¹1 (standard errors)", z_mean, 4);
    s.check("Monte Carlo: the probability of absorption by the median, ½ (standard errors)", z_median, 4);
    s.check("Monte Carlo: the law of the absorbing state, α(−T)⁻¹R (standard errors)", z_law, 4);
    double rows = 0;
    for (std::size_t a = 0; a < 4; a++) rows = std::max(rows, std::fabs(B(a, 0) + B(a, 1) - 1));
    s.check("absorption probabilities: each row sums to 1 when every exit is listed", rows, 1e-14);
  }

  // ---- gambler's ruin --------------------------------------------------------------------------------------------------------------------
  {
    const std::size_t N = 10;
    const double up = 1.2, down = 0.8, r = down / up;
    Matrix k(N + 1);
    for (std::size_t i = 1; i < N; i++) k(i, i + 1) = up, k(i, i - 1) = down;
    std::vector<std::size_t> transient;
    for (std::size_t i = 1; i < N; i++) transient.push_back(i);
    const std::vector<std::size_t> absorbing{0, N};
    const auto B = sm::absorption_probabilities(k, transient, absorbing);
    double worst = 0;
    for (std::size_t a = 0; a < transient.size(); a++) {
      const double i = static_cast<double>(transient[a]);
      worst = std::max(worst, std::fabs(B(a, 1) - (1 - std::pow(r, i)) / (1 - std::pow(r, static_cast<double>(N)))));
      worst = std::max(worst, std::fabs(B(a, 0) + B(a, 1) - 1));
    }
    s.check("gambler's ruin: P(reach N | i) = (1 − rⁱ)/(1 − rᴺ), r = μ/λ", worst, 1e-14);
  }

  // ---- mean first passage --------------------------------------------------------------------------------------------------------------------
  {
    const std::size_t n = 9;
    const double rate = 0.7;
    Matrix k(n);
    for (std::size_t i = 0; i < n; i++) k(i, (i + 1) % n) = rate, k((i + 1) % n, i) = rate;
    const auto m = sm::mean_first_passage(k, 0);
    double worst = 0;
    for (std::size_t i = 0; i < n; i++) worst = std::max(worst, std::fabs(m[i] - static_cast<double>(i * (n - i)) / (2 * rate)));
    s.check("a symmetric ring: m_i = i(n − i)/(2k)", worst, 1e-12, "s");

    double kac = 0, ph = 0;
    for (int trial = 0; trial < 10; trial++) {
      const std::size_t size = 4 + trial % 5;
      Matrix kr(size);
      for (std::size_t i = 0; i < size; i++)
        for (std::size_t j = 0; j < size; j++)
          if (i != j && (rng.uniform() < 0.5 || j == (i + 1) % size)) kr(i, j) = std::exp(rng.uniform(-1.5, 1.5));
      const auto pi = sm::stationary(kr);
      for (std::size_t target = 0; target < size; target++) {
        const auto mt = sm::mean_first_passage(kr, target);
        double sum = 0;
        for (std::size_t j = 0; j < size; j++)
          if (j != target) sum += kr(target, j) * mt[j];
        kac = std::max(kac, std::fabs(sum / (1 / pi[target] - 1) - 1));
      }
      std::vector<std::size_t> others;
      for (std::size_t i = 1; i < size; i++) others.push_back(i);
      std::vector<double> start(others.size(), 0.0);
      start.back() = 1;
      ph = std::max(ph, std::fabs(sm::phase_type(kr, others, start).moment(1) - sm::mean_first_passage(kr, 0)[size - 1]));
    }
    s.check("Kac's identity Σ_j k_ij m_j = 1/π_i − 1 on random chains, every target (relative)", kac, 1e-12);
    s.check("the phase-type mean equals the mean first passage", ph, 1e-12, "s");

    Matrix trapped(3);
    trapped(0, 1) = trapped(1, 0) = 1;
    trapped(2, 0) = 1;
    bool refused = false;
    try {
      const std::vector<std::size_t> tr{0, 1};
      const std::vector<double> al{1, 0};
      sm::phase_type(trapped, tr, al);
    } catch (const std::invalid_argument&) {
      refused = true;
    }
    s.require("a transient state that never leaves is refused", refused);
  }

  return s.done();
}
