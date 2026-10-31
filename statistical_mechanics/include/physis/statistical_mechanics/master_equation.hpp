// Markov jump processes on finitely many states: the master equation, its stationary law, its
// relaxation, and its irreversibility as currents round cycles (decision 0016).
//
// The law (van Kampen 2007, Stochastic Processes in Physics and Chemistry, 3rd ed., ch. V). States
// 0 … n−1, a rate k_ij ≥ 0 of jumping from i to j (i ≠ j; the diagonal of the matrix is ignored). A
// distribution p (a row vector) obeys the master equation
//
//   dp_j/dt = Σ_i (p_i k_ij − p_j k_ji),  i.e.  dp/dt = p Q,  Q_ij = k_ij (i ≠ j),  Q_ii = −Σ_{j≠i} k_ij,
//
// so p(t) = p(0) exp(Qt) (calx's expm). When the graph of positive rates is strongly connected
// (irreducible) there is one stationary law π Q = 0, Σ π = 1, and π > 0.
//
// The stationary law is found by the GTH algorithm (Grassmann, Taksar and Heyman 1985, Oper. Res. 33,
// 1107): states are eliminated one by one, the rates of the censored chain on the states that remain
// being k_ij + k_im k_mj / Σ_{l<m} k_ml, and π recovered by back-substitution. It forms no difference
// of positive numbers, so each π_i has a small relative error even when π spans hundreds of orders
// of magnitude, where replacing an equation of π Q = 0 by Σ π = 1 would lose the small entries. O(n³).
//
// Irreversibility (Schnakenberg 1976, Rev. Mod. Phys. 48, 571). The stationary current along an edge
// is J_ij = π_i k_ij − π_j k_ji, and the entropy production (in units of k_B per unit time)
//
//   σ = ½ Σ_{i≠j} J_ij ln(π_i k_ij / π_j k_ji) ≥ 0    (Schnakenberg eq. 7.6),
//
// each term non-negative, zero exactly under detailed balance. As J is divergence-free, it is a sum of
// cycle currents: choosing a spanning tree of the graph of positive rates (morphe's spanning forest),
// each chord (i, j), i < j, closes one fundamental cycle i → j → (the tree's path back to i), J(C) the
// chord's current and A(C) = ln Π_{a→b ∈ C} k_ab / k_ba its affinity, which no π enters. Then
//
//   σ = Σ_C J(C) A(C)    (Schnakenberg eq. 7.11),
//
// and the chain is in detailed balance, π_i k_ij = π_j k_ji for every edge, exactly when every A(C)
// vanishes: Kolmogorov's criterion, which the fundamental cycles suffice to test, affinities being
// additive over the cycle space. The relative entropy D(p(t) ‖ π) never increases (the H-theorem of
// the master equation; van Kampen 2007, §V.5).
//
// Limits. Dense: O(n³) time and O(n²) memory in stationary, evolve and irreversibility; a few thousand
// states at most. Checks: tests/check_master_equation.cpp.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <deque>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include <calx/linalg/dense.hpp>
#include <calx/linalg/expm.hpp>
#include <morphe/combinatorics/graph.hpp>

namespace physis::statistical_mechanics {

namespace detail {

// Rates k_ij ≥ 0, finite, off the diagonal of a square matrix of order ≥ 1.
inline void check_rates(const calx::linalg::Matrix& k, const char* where) {
  if (k.n == 0 || k.a.size() != k.n * k.n) throw std::invalid_argument(std::string(where) + ": a square matrix of rates of order ≥ 1");
  for (std::size_t i = 0; i < k.n; i++)
    for (std::size_t j = 0; j < k.n; j++)
      if (i != j && !(k(i, j) >= 0 && std::isfinite(k(i, j)))) throw std::invalid_argument(std::string(where) + ": rates k_ij finite and ≥ 0");
}

// The states that reach `targets` along positive rates (backward search), targets included.
inline std::vector<char> reaching(const calx::linalg::Matrix& k, const std::vector<char>& targets) {
  std::vector<char> seen = targets;
  std::deque<std::size_t> q;
  for (std::size_t i = 0; i < k.n; i++)
    if (seen[i]) q.push_back(i);
  while (!q.empty()) {
    const std::size_t j = q.front();
    q.pop_front();
    for (std::size_t i = 0; i < k.n; i++)
      if (!seen[i] && i != j && k(i, j) > 0) seen[i] = 1, q.push_back(i);
  }
  return seen;
}

}  // namespace detail

// Whether the graph of positive rates is strongly connected.
inline bool irreducible(const calx::linalg::Matrix& k) {
  detail::check_rates(k, "physis::statistical_mechanics::irreducible");
  std::vector<char> start(k.n, 0);
  start[0] = 1;
  const auto back = detail::reaching(k, start);
  std::vector<char> fwd(k.n, 0);
  fwd[0] = 1;
  std::deque<std::size_t> q{0};
  while (!q.empty()) {
    const std::size_t i = q.front();
    q.pop_front();
    for (std::size_t j = 0; j < k.n; j++)
      if (!fwd[j] && j != i && k(i, j) > 0) fwd[j] = 1, q.push_back(j);
  }
  for (std::size_t i = 0; i < k.n; i++)
    if (!back[i] || !fwd[i]) return false;
  return true;
}

// Q: Q_ij = k_ij, Q_ii = −Σ_{j≠i} k_ij.
inline calx::linalg::Matrix generator(const calx::linalg::Matrix& k) {
  detail::check_rates(k, "physis::statistical_mechanics::generator");
  calx::linalg::Matrix Q(k.n);
  for (std::size_t i = 0; i < k.n; i++) {
    double out = 0;
    for (std::size_t j = 0; j < k.n; j++)
      if (j != i) Q(i, j) = k(i, j), out += k(i, j);
    Q(i, i) = -out;
  }
  return Q;
}

// π Q = 0, Σ π = 1, by the GTH algorithm; throws unless irreducible.
inline std::vector<double> stationary(const calx::linalg::Matrix& k) {
  if (!irreducible(k)) throw std::invalid_argument("physis::statistical_mechanics::stationary: the chain is not irreducible");
  const std::size_t n = k.n;
  calx::linalg::Matrix A = k;
  for (std::size_t i = 0; i < n; i++) A(i, i) = 0;
  for (std::size_t m = n; m-- > 1;) {
    double s = 0;
    for (std::size_t j = 0; j < m; j++) s += A(m, j);
    if (!(s > 0)) throw std::runtime_error("physis::statistical_mechanics::stationary: a state with no exit in the censored chain");
    for (std::size_t i = 0; i < m; i++) A(i, m) /= s;
    for (std::size_t i = 0; i < m; i++) {
      const double aim = A(i, m);
      if (aim == 0) continue;
      for (std::size_t j = 0; j < m; j++)
        if (j != i) A(i, j) += aim * A(m, j);
    }
  }
  std::vector<double> pi(n, 0.0);
  pi[0] = 1;
  double total = 1;
  for (std::size_t m = 1; m < n; m++) {
    double v = 0;
    for (std::size_t i = 0; i < m; i++) v += pi[i] * A(i, m);
    pi[m] = v;
    total += v;
  }
  for (double& x : pi) x /= total;
  return pi;
}

// p0 exp(Q t), t ≥ 0.
inline std::vector<double> evolve(const calx::linalg::Matrix& k, std::span<const double> p0, double t) {
  if (p0.size() != k.n) throw std::invalid_argument("physis::statistical_mechanics::evolve: one probability for each state");
  if (!(t >= 0 && std::isfinite(t))) throw std::invalid_argument("physis::statistical_mechanics::evolve: a time t ≥ 0");
  const calx::linalg::Matrix E = calx::linalg::expm(t * generator(k));
  std::vector<double> p(k.n, 0.0);
  for (std::size_t i = 0; i < k.n; i++) {
    if (p0[i] == 0) continue;
    for (std::size_t j = 0; j < k.n; j++) p[j] += p0[i] * E(i, j);
  }
  return p;
}

// D(p ‖ q) = Σ p_i ln(p_i/q_i) in nats; terms with p_i = 0 vanish; +∞ when some p_i > 0 = q_i.
inline double relative_entropy(std::span<const double> p, std::span<const double> q) {
  if (p.size() != q.size()) throw std::invalid_argument("physis::statistical_mechanics::relative_entropy: p and q of one length");
  double d = 0;
  for (std::size_t i = 0; i < p.size(); i++) {
    if (p[i] < 0 || q[i] < 0) throw std::invalid_argument("physis::statistical_mechanics::relative_entropy: non-negative p and q");
    if (p[i] == 0) continue;
    if (q[i] == 0) return std::numeric_limits<double>::infinity();
    d += p[i] * std::log(p[i] / q[i]);
  }
  return d;
}

struct Cycle {
  std::vector<std::size_t> states;  // i0 → i1 → … → i_{L−1} → i0 (i0 not repeated); i0 → i1 is the chord
  double affinity = 0;              // A(C) = Σ ln(k_ab / k_ba) round the cycle
  double current = 0;               // J(C) = J_{i0 i1}, the chord's stationary current
};

struct Irreversibility {
  std::vector<double> pi;           // the stationary distribution
  calx::linalg::Matrix current;     // J_ij = π_i k_ij − π_j k_ji
  double entropy_production = 0;    // σ = ½ Σ J_ij ln(π_i k_ij / π_j k_ji) ≥ 0 (Schnakenberg eq. 7.6), per unit time
  std::vector<Cycle> cycles;        // a fundamental cycle basis of the graph of positive rates; σ = Σ_C J(C) A(C) (eq. 7.11)
  bool detailed_balance = false;    // Kolmogorov's criterion: every |A(C)| ≤ tolerance
};

// Requires k_ij > 0 ⇔ k_ji > 0, and irreducibility.
inline Irreversibility irreversibility(const calx::linalg::Matrix& k, double tolerance = 1e-10) {
  detail::check_rates(k, "physis::statistical_mechanics::irreversibility");
  const std::size_t n = k.n;
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++)
      if (i != j && (k(i, j) > 0) != (k(j, i) > 0)) throw std::invalid_argument("physis::statistical_mechanics::irreversibility: every positive rate needs its reverse (k_ij > 0 ⇔ k_ji > 0)");
  Irreversibility out;
  out.pi = stationary(k);
  out.current = calx::linalg::Matrix(n);
  morphe::combinatorics::Multigraph G;
  G.vertices = n;
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = i + 1; j < n; j++) {
      if (!(k(i, j) > 0)) continue;
      const double fwd = out.pi[i] * k(i, j), bwd = out.pi[j] * k(j, i);
      out.current(i, j) = fwd - bwd;
      out.current(j, i) = bwd - fwd;
      out.entropy_production += (fwd - bwd) * (std::log(fwd) - std::log(bwd));
      G.edges.push_back({i, j});
    }
  const auto F = morphe::combinatorics::spanning_forest(G);
  out.detailed_balance = true;
  for (std::size_t e = 0; e < G.edges.size(); e++) {
    if (F.in_tree[e]) continue;
    const std::size_t i = G.edges[e][0], j = G.edges[e][1];
    Cycle C;
    C.states = {i, j};
    std::size_t at = j;
    for (const auto& step : F.path(j, i, G)) {
      at = G.edges[step.edge][step.sign > 0 ? 1 : 0];
      if (at != i) C.states.push_back(at);
    }
    for (std::size_t s = 0; s < C.states.size(); s++) {
      const std::size_t a = C.states[s], b = C.states[(s + 1) % C.states.size()];
      C.affinity += std::log(k(a, b)) - std::log(k(b, a));
    }
    C.current = out.current(i, j);
    if (std::fabs(C.affinity) > tolerance) out.detailed_balance = false;
    out.cycles.push_back(std::move(C));
  }
  return out;
}

}  // namespace physis::statistical_mechanics
