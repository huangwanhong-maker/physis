// The entropy of mixing, and the departure from the fully mixed state (decision 0016).
//
// The law. Counts n_xq of kinds q at places x (a table, places × kinds), N = Σ n, p_xq = n_xq/N, with
// marginals p_x and p_q. The entropy of the kinds and its conditional on the place,
//
//   H(Q) = −Σ_q p_q ln p_q,   H(Q|X) = Σ_x p_x H(p(·|x)),
//
// differ by the mutual information between place and kind,
//
//   M = H(Q) − H(Q|X) = Σ_x p_x D(p(·|x) ‖ p(·)) = Σ_{x,q} p_xq ln(p_xq / p_x p_q) ≥ 0,
//
// zero exactly when every place holds the kinds in the whole's proportions, and H(Q) exactly when each
// place holds one kind (Theil and Finizza 1971, J. Math. Sociol. 1, 187: the information index of
// segregation; Mora and Ruiz-Castillo 2011, Sociol. Methodol. 41, 159; Cover and Thomas 2006, Elements
// of Information Theory, 2nd ed., §2.4). Over a partition of the places into groups, the chain rule
// of mutual information (Cover and Thomas 2006, Theorem 2.5.2), the group being a function of the
// place, gives exactly
//
//   M = M_between + Σ_g p_g M_g,
//
// M_between the M of the table summed over each group, M_g that of group g's places among themselves.
// M is invariant under permuting places or kinds, and does not increase when places exchange their
// compositions by a stochastic map (the data-processing inequality; Cover and Thomas 2006, Theorem 2.8.1).
//
// Computation. M is summed as its non-negative local terms p_x D(p(·|x) ‖ p(·)), each ratio formed as
// n_xq N / (n_x n_q), exact for integer counts below 2⁵³ in their products: M is zero to rounding for
// proportional counts. All in nats. Checks: tests/check_mixing.cpp.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

#include <calx/linalg/rectangular.hpp>

namespace physis::statistical_mechanics {

struct Mixing {
  double entropy = 0;         // H(Q), nats
  double conditional = 0;     // H(Q|X)
  double departure = 0;       // M
  double normalised = 0;      // M / H(Q) ∈ [0, 1] (0 if H(Q) = 0)
  std::vector<double> local;  // p_x D(p(·|x) ‖ p(·)) per place; sums to M
};

inline Mixing mixing(const calx::linalg::Rectangular& counts) {
  if (counts.rows == 0 || counts.cols == 0) throw std::invalid_argument("physis::statistical_mechanics::mixing: a table of places × kinds");
  std::vector<double> row(counts.rows, 0.0), col(counts.cols, 0.0);
  double N = 0;
  for (std::size_t x = 0; x < counts.rows; x++)
    for (std::size_t q = 0; q < counts.cols; q++) {
      const double n = counts(x, q);
      if (!(n >= 0 && std::isfinite(n))) throw std::invalid_argument("physis::statistical_mechanics::mixing: counts finite and ≥ 0");
      row[x] += n, col[q] += n, N += n;
    }
  if (!(N > 0)) throw std::invalid_argument("physis::statistical_mechanics::mixing: a table with something in it");
  Mixing m;
  for (double c : col)
    if (c > 0) m.entropy -= (c / N) * std::log(c / N);
  m.local.assign(counts.rows, 0.0);
  for (std::size_t x = 0; x < counts.rows; x++) {
    if (row[x] == 0) continue;
    double h = 0, d = 0;
    for (std::size_t q = 0; q < counts.cols; q++) {
      const double n = counts(x, q);
      if (n == 0) continue;
      h -= (n / row[x]) * std::log(n / row[x]);
      d += (n / N) * std::log((n * N) / (row[x] * col[q]));
    }
    m.conditional += (row[x] / N) * h;
    m.local[x] = std::max(d, 0.0);  // a KL divergence; a negative is rounding
    m.departure += m.local[x];
  }
  m.normalised = m.entropy > 0 ? std::min(1.0, m.departure / m.entropy) : 0.0;
  return m;
}

struct MixingDecomposition {
  double total = 0, between = 0, within = 0;  // total = between + within; within = Σ_g p_g M_g
  std::vector<double> within_group;           // M_g, the departure among group g's places (0 for an empty group)
};

// group[x] ∈ 0 … G−1 for each place; G is the largest label plus one.
inline MixingDecomposition mixing_decomposition(const calx::linalg::Rectangular& counts, std::span<const std::size_t> group) {
  if (group.size() != counts.rows) throw std::invalid_argument("physis::statistical_mechanics::mixing_decomposition: one group for each place");
  MixingDecomposition d;
  d.total = mixing(counts).departure;
  std::size_t G = 0;
  for (std::size_t g : group) G = std::max(G, g + 1);
  calx::linalg::Rectangular by_group(G, counts.cols);
  std::vector<double> mass(G, 0.0);
  double N = 0;
  for (std::size_t x = 0; x < counts.rows; x++)
    for (std::size_t q = 0; q < counts.cols; q++) by_group(group[x], q) += counts(x, q), mass[group[x]] += counts(x, q), N += counts(x, q);
  d.between = mixing(by_group).departure;
  d.within_group.assign(G, 0.0);
  for (std::size_t g = 0; g < G; g++) {
    if (mass[g] == 0) continue;
    std::vector<std::size_t> places;
    for (std::size_t x = 0; x < counts.rows; x++)
      if (group[x] == g) places.push_back(x);
    calx::linalg::Rectangular sub(places.size(), counts.cols);
    for (std::size_t r = 0; r < places.size(); r++)
      for (std::size_t q = 0; q < counts.cols; q++) sub(r, q) = counts(places[r], q);
    d.within_group[g] = mixing(sub).departure;
    d.within += (mass[g] / N) * d.within_group[g];
  }
  return d;
}

}  // namespace physis::statistical_mechanics
