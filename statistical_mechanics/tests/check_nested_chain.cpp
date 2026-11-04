// Checks of chains on a classified set (nested_chain.hpp).
//
// On an uneven classification of 9 states in 3 families, 5 groups (a family of one group, one of one class) with an
// uneven law:
// - the closed form of P_τ equals calx's matrix exponential of the generator, at three times;
// - the chain is reversible with its law, and its law stationary;
// - each mode's projection is an eigenvector of the generator with eigenvalue −μ_m, μ_m the sum of the coarser rates,
//   and the projections sum to the identity and are orthogonal (Q_a Q_b = δ_ab Q_b);
// - the mode weights of a transition sum to 1;
// - the expected number of changes between two observations equals Van Loan's integral ∫ e^{sL} L_off e^{(τ−s)L} ds,
//   from the matrix exponential of the augmented generator, at every pair;
// - the gradient of the log-likelihood equals central differences;
// - rates recovered from 60 000 simulated transitions within three standard errors, the memories in increasing order;
// - faulty classifications refused (a level that does not refine, a law that is not one).

#include <cmath>
#include <string>
#include <vector>

#include <calx/linalg/dense.hpp>
#include <calx/linalg/expm.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <physis/statistical_mechanics/nested_chain.hpp>

using calx::verification::Suite;
namespace sm = physis::statistical_mechanics;

namespace {

// Families {0,1,2,3}, {4,5,6,7}, {8}; groups {0,1}, {2,3}, {4,5,6,7}, {8}; the states themselves.
sm::Classification uneven() {
  sm::Classification C;
  C.block = {{0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 1, 1, 1, 2}, {0, 0, 1, 1, 2, 2, 2, 2, 3}, {0, 1, 2, 3, 4, 5, 6, 7, 8}};
  C.law = {0.05, 0.15, 0.1, 0.1, 0.2, 0.05, 0.05, 0.1, 0.2};
  return C;
}

double max_abs_diff(const calx::linalg::Matrix& A, const calx::linalg::Matrix& B) {
  double d = 0;
  for (std::size_t i = 0; i < A.a.size(); i++) d = std::max(d, std::fabs(A.a[i] - B.a[i]));
  return d;
}

calx::linalg::Matrix scaled(const calx::linalg::Matrix& A, double s) {
  auto B = A;
  for (double& v : B.a) v *= s;
  return B;
}

bool refused(const sm::Classification& C) {
  try {
    sm::validate(C);
  } catch (const std::invalid_argument&) {
    return true;
  }
  return false;
}

}  // namespace

int main() {
  Suite s;
  const auto C = uneven();
  const sm::NestedChain X(C, {0.07, 0.2, 0.9});
  const std::size_t N = C.states();
  const auto G = X.generator();
  // The closed form against the matrix exponential.
  double worst = 0;
  for (double tau : {0.3, 1.0, 4.0}) worst = std::max(worst, max_abs_diff(X.kernel_matrix(tau), calx::linalg::expm(scaled(G, tau))));
  s.require("P_τ in closed form equals exp(τL) at τ = 0.3, 1, 4", worst < 1e-12, std::to_string(worst));
  // Reversible, its law stationary.
  double balance = 0, stationary = 0;
  for (std::size_t x = 0; x < N; x++) {
    double col = 0;
    for (std::size_t y = 0; y < N; y++) balance = std::max(balance, std::fabs(C.law[x] * G(x, y) - C.law[y] * G(y, x))), col += C.law[y] * G(y, x);
    stationary = std::max(stationary, std::fabs(col));
  }
  s.require("detailed balance with the law, and the law stationary", balance < 1e-15 && stationary < 1e-15);
  // The modes.
  bool eigen = true, sum = true, orthogonal = true;
  for (std::size_t m = 0; m <= C.levels(); m++) {
    calx::linalg::Matrix Q(N);
    for (std::size_t x = 0; x < N; x++)
      for (std::size_t y = 0; y < N; y++) Q(x, y) = X.mode(m, x, y);
    const auto GQ = G * Q;
    eigen = eigen && max_abs_diff(GQ, scaled(Q, -X.mu[m])) < 1e-14;
    for (std::size_t b = 0; b <= C.levels(); b++) {
      calx::linalg::Matrix R(N);
      for (std::size_t x = 0; x < N; x++)
        for (std::size_t y = 0; y < N; y++) R(x, y) = X.mode(b, x, y);
      const auto QR = Q * R;
      orthogonal = orthogonal && max_abs_diff(QR, b == m ? Q : calx::linalg::Matrix(N)) < 1e-14;
    }
  }
  for (std::size_t x = 0; x < N; x++)
    for (std::size_t y = 0; y < N; y++) {
      double t = 0;
      for (std::size_t m = 0; m <= C.levels(); m++) t += X.mode(m, x, y);
      sum = sum && std::fabs(t - (x == y ? 1.0 : 0.0)) < 1e-15;
    }
  s.require("each mode an eigenspace of the generator, eigenvalue −(λ_0 + … + λ_{m−1}): μ = 0, 0.07, 0.27, 1.17",
            eigen && std::fabs(X.mu[1] - 0.07) < 1e-15 && std::fabs(X.mu[2] - 0.27) < 1e-15 && std::fabs(X.mu[3] - 1.17) < 1e-15);
  s.require("the modes' projections sum to the identity and are orthogonal", sum && orthogonal);
  double wsum = 0;
  for (std::size_t x = 0; x < N; x++)
    for (std::size_t y = 0; y < N; y++) {
      double t = 0;
      for (double w : X.mode_weights(x, y, 1.7)) t += w;
      wsum = std::max(wsum, std::fabs(t - 1));
    }
  s.require("the mode weights of every transition sum to 1", wsum < 1e-12);
  // The expected number of changes against Van Loan's integral.
  {
    const double tau = 2.5;
    calx::linalg::Matrix A(2 * N);
    for (std::size_t i = 0; i < N; i++)
      for (std::size_t j = 0; j < N; j++) {
        A(i, j) = tau * G(i, j), A(N + i, N + j) = tau * G(i, j);
        if (i != j) A(i, N + j) = tau * G(i, j);
      }
    const auto E = calx::linalg::expm(A);
    const auto P = X.kernel_matrix(tau);
    double d = 0;
    for (std::size_t x = 0; x < N; x++)
      for (std::size_t y = 0; y < N; y++) d = std::max(d, std::fabs(X.expected_changes(x, y, tau) - E(x, N + y) / P(x, y)));
    s.require("the expected changes between two observations equal Van Loan's integral, at every pair", d < 1e-9, std::to_string(d));
  }
  // The gradient of the log-likelihood against central differences.
  std::vector<sm::Transition> T = {{0, 0, 1, 5}, {0, 1, 1, 2}, {1, 3, 2, 1}, {4, 8, 1, 1}, {8, 8, 3, 4}, {5, 6, 0.5, 2}};
  {
    std::vector<double> g;
    const double v = sm::log_likelihood(X, T, &g);
    double d = 0;
    for (std::size_t l = 0; l < C.levels(); l++) {
      auto up = X.lambda, dn = X.lambda;
      up[l] += 1e-6, dn[l] -= 1e-6;
      const double fd = (sm::log_likelihood(sm::NestedChain(C, up), T) - sm::log_likelihood(sm::NestedChain(C, dn), T)) / 2e-6;
      d = std::max(d, std::fabs(fd - g[l]));
    }
    s.require("the log-likelihood's gradient in the rates equals central differences", std::isfinite(v) && d < 1e-6, std::to_string(d));
  }
  // Rates recovered from simulated transitions: x from the law, the deepest level kept drawn from its probability,
  // y drawn from the law within x's block there.
  {
    calx::random::Stream rng(20261104, 20, 0, 0);
    const std::vector<double> truth = {0.05, 0.3, 1.2};
    const sm::NestedChain Y(C, truth);
    std::vector<std::vector<std::size_t>> cell(N * N);
    std::vector<double> counts(N * N, 0.0);
    const auto draw = [&](const std::vector<double>& w) {
      double u = rng.uniform(), t = 0;
      for (std::size_t i = 0; i < w.size(); i++) if ((t += w[i]) >= u) return i;
      return w.size() - 1;
    };
    const double tau = 1;
    for (int k = 0; k < 60000; k++) {
      const std::size_t x = draw(C.law);
      std::vector<double> keep(C.levels() + 1);
      for (std::size_t m = 0; m <= C.levels(); m++) keep[m] = std::exp(-Y.mu[m] * tau) - (m == C.levels() ? 0.0 : std::exp(-Y.mu[m + 1] * tau));
      const std::size_t m = draw(keep);
      std::vector<double> w(N);
      for (std::size_t y = 0; y < N; y++) w[y] = Y.redraw(m, x, y);
      counts[x * N + draw(w)] += 1;
    }
    std::vector<sm::Transition> sim;
    for (std::size_t x = 0; x < N; x++)
      for (std::size_t y = 0; y < N; y++)
        if (counts[x * N + y] > 0) sim.push_back({x, y, tau, counts[x * N + y]});
    const auto F = sm::fit_rates(C, sim, {0.1, 0.1, 0.1});
    bool within = F.converged;
    std::string detail;
    for (std::size_t l = 0; l < C.levels(); l++) {
      within = within && std::fabs(F.lambda[l] - truth[l]) < 3 * F.lambda_error[l];
      detail += std::to_string(F.lambda[l]) + "±" + std::to_string(F.lambda_error[l]) + " ";
    }
    s.require("rates recovered from 60 000 simulated transitions within three standard errors", within, detail);
    s.require("... the memories in increasing order: a coarser level is kept at least as long", F.mu[1] <= F.mu[2] && F.mu[2] <= F.mu[3]);
    // A level with no memory of its own: transitions drawn with λ_1 = 0, its rate driven to the bound.
    const sm::NestedChain Z(C, {0.05, 0.0, 1.2});
    std::vector<double> c2(N * N, 0.0);
    for (int k = 0; k < 60000; k++) {
      const std::size_t x = draw(C.law);
      std::vector<double> keep(C.levels() + 1);
      for (std::size_t m = 0; m <= C.levels(); m++) keep[m] = std::exp(-Z.mu[m] * tau) - (m == C.levels() ? 0.0 : std::exp(-Z.mu[m + 1] * tau));
      const std::size_t m = draw(keep);
      std::vector<double> w(N);
      for (std::size_t y = 0; y < N; y++) w[y] = Z.redraw(m, x, y);
      c2[x * N + draw(w)] += 1;
    }
    std::vector<sm::Transition> sim2;
    for (std::size_t x = 0; x < N; x++)
      for (std::size_t y = 0; y < N; y++)
        if (c2[x * N + y] > 0) sim2.push_back({x, y, tau, c2[x * N + y]});
    const auto F2 = sm::fit_rates(C, sim2, {0.1, 0.1, 0.1});
    s.require("a level with no memory of its own: its rate found at most a few hundredths, or at its bound with no error, the fit converged",
              F2.converged && (F2.at_bound[1] ? std::isnan(F2.lambda_error[1]) : F2.lambda[1] < 3 * F2.lambda_error[1] + 1e-3),
              std::to_string(F2.lambda[1]) + (F2.at_bound[1] ? " at its bound" : ""));
  }
  // Faulty classifications.
  {
    auto bad = uneven();
    bad.block[2][3] = 3;  // state 3 in group 3 with state 8 though in another family
    bad.block[2][8] = 3;
    auto law = uneven();
    law.law[0] = 0.5;
    s.require("a level that does not refine the one before, and a law that is not one, are refused", refused(bad) && refused(law));
  }
  return s.done();
}
