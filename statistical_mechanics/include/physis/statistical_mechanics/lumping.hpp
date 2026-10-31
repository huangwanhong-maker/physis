// Lumping a Markov jump process: when a partition of its states makes the process of blocks Markov
// (strong lumpability), the test whole and mechanism by mechanism, the lumped rates and the intertwining of
// the generators; and, when it fails, the exact law of the blocks observed — the aggregated process, whose
// path density is the effective action of the coarse description and whose dwell in a block is phase-type
// (decision 0018).
//
// Strong lumpability (Kemeny and Snell 1960, Finite Markov Chains, §6.3; in continuous time Buchholz 1994,
// J. Appl. Probab. 31, 59–75). States 0…n−1, rates k_ij, a partition C into blocks. The process C(Z_t) is a
// Markov chain for every initial law when, for every two states i, i′ of one block A and every block
// B ≠ A, Σ_{j∈B} k_ij = Σ_{j∈B} k_i′j; its rates k̄_AB are those common sums. Equivalently the generators
// intertwine, Q (f∘C) = (Q̄ f)∘C for every f on the blocks, and then exp(Qt)(f∘C) = (exp(Q̄t) f)∘C: the
// coarse observables evolve by the coarse generator alone.
//
// Mechanisms. When the rates are a sum Σ_m k^(m) of declared mechanisms (event kinds), lumpability of each
// suffices — the lumped rates add — and is not necessary: the failures of two mechanisms can cancel in
// their sum (two states of a block, one leaving by mechanism 1, the other by mechanism 2, at one rate).
//
// The aggregated process (Colquhoun and Hawkes 1981, Proc. R. Soc. Lond. B 211, 205–235; Fredkin and Rice
// 1986, J. Appl. Probab. 23, 208–214). The blocks are observed whether or not they lump. A coarse path that
// starts in block b₀, stays τ₀, jumps to b₁, …, and is still in b_m when observation ends after τ_m has, with
// respect to Lebesgue measure on the jump times and counting measure on the blocks, the density
//
//   p = p₀|_{b₀} e^{Q_{b₀b₀} τ₀} Q_{b₀b₁} e^{Q_{b₁b₁} τ₁} ⋯ Q_{b_{m−1}b_m} e^{Q_{b_m b_m} τ_m} 1,
//
// Q_{AB} the block of the generator from A to B (Q_{AA} keeps the diagonal −Σ_j k_ij). Its negative logarithm
// is the effective action of the coarse description relative to that reference: −log of the conditional
// expectation of the fine path density given the coarse path, made explicit. Under strong lumpability it is
// the lumped chain's −log p̄ for every initial law; otherwise it depends on the whole past, and the dwell in
// a block entered with law α is phase-type with sub-generator Q_{AA} (first_passage.hpp), not exponential:
// eliminating the states within a block leaves memory.
//
// Limits: dense, O(n²) for the tests and O(m n³) for a path of m jumps. Checks: tests/check_lumping.cpp.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

#include <calx/linalg/dense.hpp>
#include <calx/linalg/expm.hpp>

#include "first_passage.hpp"
#include "master_equation.hpp"

namespace physis::statistical_mechanics {

// A partition: the block of each state, blocks numbered 0…b−1, each nonempty.
inline std::size_t block_count(std::span<const std::size_t> block, std::size_t n) {
  if (block.size() != n) throw std::invalid_argument("physis::statistical_mechanics: one block for each state");
  std::size_t b = 0;
  for (std::size_t x : block) b = std::max(b, x + 1);
  std::vector<char> used(b, 0);
  for (std::size_t x : block) used[x] = 1;
  for (char u : used)
    if (!u) throw std::invalid_argument("physis::statistical_mechanics: the blocks must be numbered 0…b−1, each nonempty");
  return b;
}

struct Lumpability {
  bool lumpable = true;
  double discrepancy = 0;   // the largest |Σ_{j∈B} k_ij − Σ_{j∈B} k_i′j| over i, i′ in one block A, B ≠ A
  std::size_t from = 0, other = 0, to_block = 0;  // where it is largest: states i, i′ and the block B
};

// Rates into each block from each state: R(i, B) = Σ_{j∈B, j≠i} k_ij.
inline std::vector<std::vector<double>> rates_into_blocks(const calx::linalg::Matrix& k, std::span<const std::size_t> block) {
  const std::size_t b = block_count(block, k.n);
  std::vector<std::vector<double>> R(k.n, std::vector<double>(b, 0.0));
  for (std::size_t i = 0; i < k.n; i++)
    for (std::size_t j = 0; j < k.n; j++)
      if (j != i) R[i][block[j]] += k(i, j);
  return R;
}

// Strong lumpability of the rates k under the partition, to a tolerance relative to the largest rate.
inline Lumpability strongly_lumpable(const calx::linalg::Matrix& k, std::span<const std::size_t> block, double tolerance = 1e-12) {
  detail::check_rates(k, "physis::statistical_mechanics::strongly_lumpable");
  const auto R = rates_into_blocks(k, block);
  double scale = 0;
  for (std::size_t i = 0; i < k.n; i++)
    for (std::size_t j = 0; j < k.n; j++)
      if (i != j) scale = std::max(scale, k(i, j));
  Lumpability L;
  const std::size_t b = R.empty() ? 0 : R[0].size();
  std::vector<std::size_t> first(b, k.n);
  for (std::size_t i = 0; i < k.n; i++)
    if (first[block[i]] == k.n) first[block[i]] = i;
  for (std::size_t i = 0; i < k.n; i++) {
    const std::size_t r = first[block[i]];
    for (std::size_t B = 0; B < b; B++) {
      if (B == block[i]) continue;
      const double d = std::fabs(R[i][B] - R[r][B]);
      if (d > L.discrepancy) L.discrepancy = d, L.from = r, L.other = i, L.to_block = B;
    }
  }
  L.lumpable = L.discrepancy <= tolerance * std::max(scale, 1e-300);
  return L;
}

// Lumpability of each mechanism (k^(m), summed to the process's rates): lumpable when every one is; the
// report is that of the worst.
inline Lumpability eventwise_lumpable(const std::vector<calx::linalg::Matrix>& mechanisms, std::span<const std::size_t> block, double tolerance = 1e-12) {
  Lumpability worst;
  for (const auto& k : mechanisms) {
    const auto L = strongly_lumpable(k, block, tolerance);
    if (!L.lumpable) worst.lumpable = false;
    if (L.discrepancy > worst.discrepancy) {
      const bool keep = worst.lumpable;
      worst = L;
      worst.lumpable = keep && L.lumpable;
    }
  }
  return worst;
}

// The lumped rates k̄_AB, read from the first state of each block (meaningful when strongly lumpable).
inline calx::linalg::Matrix lumped(const calx::linalg::Matrix& k, std::span<const std::size_t> block) {
  const auto R = rates_into_blocks(k, block);
  const std::size_t b = R.empty() ? 0 : R[0].size();
  calx::linalg::Matrix out(b);
  std::vector<char> done(b, 0);
  for (std::size_t i = 0; i < k.n; i++) {
    if (done[block[i]]) continue;
    done[block[i]] = 1;
    for (std::size_t B = 0; B < b; B++)
      if (B != block[i]) out(block[i], B) = R[i][B];
  }
  return out;
}

// Approximate lumping. When the blocks do not lump, a coarse chain can still be chosen — the aggregation by
// weights w (the stationary law, say): k̄_AB = Σ_{i∈A} w_i R(i, B) / Σ_{i∈A} w_i, which is exact when the
// chain lumps and, when w is the stationary law, keeps the blocks' stationary law (π̄_A k̄_AB is the stationary
// flux from A to B, and the fluxes balance on every block) — and its error bounded. With the outward discrepancy
//
//   δ = max_i Σ_{B ≠ C(i)} |R(i, B) − k̄_{C(i) B}|   (per unit time),
//
// a coupling that lets the two chains jump together at the smaller rate toward each block, and the fine chain
// move freely inside a block, keeps their block paths equal until a first mismatch whose hazard is at most δ;
// so over [0, T] the laws of the whole block path of the fine chain and of the coarse chain's path differ in
// total variation by at most 1 − e^{−δT} ≤ δT (the coupling inequality; Lindvall 1992, Lectures on the
// Coupling Method, ch. I). The bound is conservative and linear in T; it is zero exactly when the chain lumps
// to k̄.
inline calx::linalg::Matrix aggregated(const calx::linalg::Matrix& k, std::span<const std::size_t> block, std::span<const double> weight) {
  if (weight.size() != k.n) throw std::invalid_argument("physis::statistical_mechanics::aggregated: one weight for each state");
  const auto R = rates_into_blocks(k, block);
  const std::size_t b = R.empty() ? 0 : R[0].size();
  calx::linalg::Matrix out(b);
  std::vector<double> mass(b, 0.0);
  for (std::size_t i = 0; i < k.n; i++) {
    if (!(weight[i] >= 0)) throw std::invalid_argument("physis::statistical_mechanics::aggregated: weights ≥ 0");
    mass[block[i]] += weight[i];
    for (std::size_t B = 0; B < b; B++)
      if (B != block[i]) out(block[i], B) += weight[i] * R[i][B];
  }
  for (std::size_t A = 0; A < b; A++) {
    if (!(mass[A] > 0)) throw std::invalid_argument("physis::statistical_mechanics::aggregated: every block needs positive weight");
    for (std::size_t B = 0; B < b; B++) out(A, B) /= mass[A];
  }
  return out;
}

// δ = max_i Σ_{B ≠ C(i)} |R(i, B) − k̄_{C(i) B}|, the outward discrepancy of a coarse chain k̄ from the fine one.
inline double outward_discrepancy(const calx::linalg::Matrix& k, std::span<const std::size_t> block, const calx::linalg::Matrix& coarse) {
  const auto R = rates_into_blocks(k, block);
  double worst = 0;
  for (std::size_t i = 0; i < k.n; i++) {
    double s = 0;
    for (std::size_t B = 0; B < coarse.n; B++)
      if (B != block[i]) s += std::fabs(R[i][B] - coarse(block[i], B));
    worst = std::max(worst, s);
  }
  return worst;
}

// max over states i and blocks B of |(Q 1_B)(i) − (Q̄ 1_B)(C(i))|, the indicators 1_B spanning the functions
// on the blocks: zero exactly when the generators intertwine.
inline double intertwining_defect(const calx::linalg::Matrix& k, std::span<const std::size_t> block, const calx::linalg::Matrix& coarse) {
  const auto Q = generator(k), Qb = generator(coarse);
  double worst = 0;
  for (std::size_t i = 0; i < k.n; i++)
    for (std::size_t B = 0; B < coarse.n; B++) {
      double fine = 0;
      for (std::size_t j = 0; j < k.n; j++)
        if (block[j] == B) fine += Q(i, j);
      worst = std::max(worst, std::fabs(fine - Qb(block[i], B)));
    }
  return worst;
}

// A coarse path: the blocks visited and the time in each; the last time is censored (still there).
struct CoarsePath {
  std::vector<std::size_t> blocks;
  std::vector<double> dwell;
};

// log p of a coarse path under the fine chain started from p0 (the aggregated process). −log p is the
// effective action. −∞ when the path is impossible.
inline double coarse_log_density(const calx::linalg::Matrix& k, std::span<const std::size_t> block, std::span<const double> p0, const CoarsePath& path) {
  if (p0.size() != k.n) throw std::invalid_argument("physis::statistical_mechanics::coarse_log_density: one probability for each state");
  if (path.blocks.empty() || path.blocks.size() != path.dwell.size()) throw std::invalid_argument("physis::statistical_mechanics::coarse_log_density: a block and a dwell for each step");
  block_count(block, k.n);
  const auto Q = generator(k);
  // The states of each block, and the sub-generator exponentials as needed.
  const auto states_of = [&](std::size_t B) {
    std::vector<std::size_t> s;
    for (std::size_t i = 0; i < k.n; i++)
      if (block[i] == B) s.push_back(i);
    return s;
  };
  double logp = 0;
  std::vector<std::size_t> here = states_of(path.blocks[0]);
  std::vector<double> v(here.size());
  for (std::size_t a = 0; a < here.size(); a++) v[a] = p0[here[a]];
  for (std::size_t step = 0; step < path.blocks.size(); step++) {
    const double tau = path.dwell[step];
    if (!(tau >= 0 && std::isfinite(tau))) throw std::invalid_argument("physis::statistical_mechanics::coarse_log_density: dwell times ≥ 0");
    calx::linalg::Matrix S(here.size());
    for (std::size_t a = 0; a < here.size(); a++)
      for (std::size_t c = 0; c < here.size(); c++) S(a, c) = Q(here[a], here[c]) * tau;
    const auto E = calx::linalg::expm(S);
    std::vector<double> w(here.size(), 0.0);
    for (std::size_t a = 0; a < here.size(); a++)
      if (v[a] != 0)
        for (std::size_t c = 0; c < here.size(); c++) w[c] += v[a] * E(a, c);
    if (step + 1 == path.blocks.size()) {
      double total = 0;
      for (double x : w) total += x;
      return total > 0 ? logp + std::log(total) : -std::numeric_limits<double>::infinity();
    }
    const auto there = states_of(path.blocks[step + 1]);
    std::vector<double> u(there.size(), 0.0);
    for (std::size_t a = 0; a < here.size(); a++)
      if (w[a] != 0)
        for (std::size_t c = 0; c < there.size(); c++) u[c] += w[a] * Q(here[a], there[c]);
    double norm = 0;
    for (double x : u) norm += x;
    if (!(norm > 0)) return -std::numeric_limits<double>::infinity();
    for (double& x : u) x /= norm;
    logp += std::log(norm);
    here = there;
    v = std::move(u);
  }
  return logp;
}

// The law of the dwell in block B entered with law alpha on its states: phase-type with sub-generator Q_BB.
inline PhaseType dwell_law(const calx::linalg::Matrix& k, std::span<const std::size_t> block, std::size_t B, std::span<const double> alpha) {
  std::vector<std::size_t> inside;
  for (std::size_t i = 0; i < k.n; i++)
    if (block[i] == B) inside.push_back(i);
  return phase_type(k, inside, alpha);
}

}  // namespace physis::statistical_mechanics
