// Causal sets: spacetime as a locally finite partial order, sprinkled into Minkowski space.
//
// A causal set is a set C with a partial order ≺ that is locally finite: every interval
// I(x, y) = {z : x ≺ z ≺ y} is finite (Bombelli, Lee, Meyer and Sorkin 1987, Phys. Rev.
// Lett. 59, 521; Surya 2019, Living Rev. Relativ. 22, 5). The causal-set hypothesis is that
// order plus number is geometry: the order gives the causal structure, and the number of
// elements the volume, at a density ρ. A causal set faithful to a Lorentzian manifold is
// obtained by sprinkling: a Poisson process of density ρ in the manifold, ordered by its
// causal relation, so that the number of elements in a region of volume V is Poisson with
// mean ρV, which no lattice can achieve without breaking Lorentz invariance.
//
// Here elements are sprinkled into an Alexandrov interval of ℝ^{1,d-1} (minkowski.hpp) of
// proper time τ, with tips at ct = ∓τ/2 on the axis, and ordered by causal precedence. The
// order is kept as bitsets: the future and the past of each element, so that a relation is
// one bit, the interval between two elements the popcount of an intersection, and a link
// (x ≺ y with nothing between) an empty intersection. Memory 2N²/8 bytes; links and
// interval abundances cost O(R N / 64) for R relations, the longest chain O(N²).
//
// Estimators of geometry from order alone:
// - the ordering fraction r = R / C(N, 2), whose expectation in a d-dimensional interval is
//     f(d) = Γ(d + 1) Γ(d/2) / (2 Γ(3d/2))
//   (Myrheim 1978, CERN preprint TH-2538; Meyer 1988, PhD thesis, MIT): 1 in one dimension,
//   1/2 in two, 1/10 in four; inverting f gives the Myrheim-Meyer dimension;
// - the longest chain, whose length L grows as m_d (ρV)^{1/d} with m_2 = 2 (Brightwell and
//   Gregory 1991, Phys. Rev. Lett. 66, 260; Bollobás and Winkler 1988): in two dimensions the
//   order of a sprinkled interval is that of a random permutation, and L is the longest
//   increasing subsequence, whose mean is 2√N - 1.7711 N^{1/6} + o(N^{1/6}) (Baik, Deift
//   and Johansson 1999, J. Amer. Math. Soc. 12, 1119);
// - the abundances N_k of intervals with k elements inside (N_0 = the links), the input of
//   the Benincasa-Dowker action (planned).

#pragma once

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

#include <calx/random/stream.hpp>

#include "minkowski.hpp"

namespace physis::spacetime {

// Elements sprinkled at density rho (per m^d) into the Alexandrov interval of proper time
// tau (m) in ℝ^{1,d-1}, its tips at ct = -tau/2 and +tau/2 on the axis; their number is
// Poisson with mean rho ζ_d tau^d.
template <std::size_t d>
inline std::vector<Event<d>> sprinkle_diamond(double tau, double rho, calx::random::Stream& z) {
  const auto count = z.poisson(rho * alexandrov_volume(d, tau));
  std::vector<Event<d>> events;
  events.reserve(static_cast<std::size_t>(count));
  const double h = tau / 2;
  while (static_cast<std::int64_t>(events.size()) < count) {
    Event<d> e;
    for (std::size_t a = 0; a < d; a++) e[a] = z.uniform(-h, h);
    double r2 = 0;
    for (std::size_t a = 1; a < d; a++) r2 += e[a] * e[a];
    const double room = h - std::fabs(e[0]);
    if (r2 < room * room) events.push_back(e);
  }
  return events;
}

// f(d), the expected ordering fraction of an interval of ℝ^{1,d-1}, for real d >= 1.
inline double expected_ordering_fraction(double d) {
  return 0.5 * std::exp(std::lgamma(d + 1) + std::lgamma(d / 2) - std::lgamma(1.5 * d));
}

// The Myrheim-Meyer dimension: the d in [1, 20] with f(d) = r, by bisection (f decreases).
inline double myrheim_meyer_dimension(double r) {
  double lo = 1, hi = 20;
  if (r >= expected_ordering_fraction(lo)) return lo;
  if (r <= expected_ordering_fraction(hi)) return hi;
  for (int it = 0; it < 200 && hi - lo > 1e-13; it++) {
    const double mid = 0.5 * (lo + hi);
    (expected_ordering_fraction(mid) > r ? lo : hi) = mid;
  }
  return 0.5 * (lo + hi);
}

template <std::size_t d>
class CausalSet {
 public:
  // The causal set of the given events, ordered by causal precedence in ℝ^{1,d-1}. The
  // elements are renumbered in order of x^0, so that i ≺ j implies i < j.
  explicit CausalSet(std::vector<Event<d>> events) : events_(std::move(events)) {
    std::sort(events_.begin(), events_.end(), [](const Event<d>& a, const Event<d>& b) { return a[0] < b[0]; });
    n_ = events_.size();
    words_ = (n_ + 63) / 64;
    future_.assign(n_ * words_, 0);
    past_.assign(n_ * words_, 0);
    for (std::size_t i = 0; i < n_; i++)
      for (std::size_t j = i + 1; j < n_; j++)
        if (causally_precedes(events_[i], events_[j])) {
          future_[i * words_ + j / 64] |= std::uint64_t(1) << (j % 64);
          past_[j * words_ + i / 64] |= std::uint64_t(1) << (i % 64);
          relations_++;
        }
  }

  std::size_t size() const { return n_; }
  const std::vector<Event<d>>& events() const { return events_; }
  // i ≺ j.
  bool precedes(std::size_t i, std::size_t j) const { return (future_[i * words_ + j / 64] >> (j % 64)) & 1; }
  // The number of related pairs R.
  std::size_t relations() const { return relations_; }
  // R / C(N, 2).
  double ordering_fraction() const { return n_ < 2 ? 0.0 : static_cast<double>(relations_) / (0.5 * static_cast<double>(n_) * static_cast<double>(n_ - 1)); }
  // |I(i, j)|, the number of elements strictly between i and j (i ≺ j).
  std::size_t interval_size(std::size_t i, std::size_t j) const {
    std::size_t c = 0;
    for (std::size_t w = 0; w < words_; w++) c += std::popcount(future_[i * words_ + w] & past_[j * words_ + w]);
    return c;
  }
  // The links: related pairs with nothing between them (the Hasse diagram).
  std::vector<std::pair<std::size_t, std::size_t>> links() const {
    std::vector<std::pair<std::size_t, std::size_t>> l;
    for (std::size_t i = 0; i < n_; i++)
      for (std::size_t j = i + 1; j < n_; j++)
        if (precedes(i, j) && empty_between(i, j)) l.emplace_back(i, j);
    return l;
  }
  // N_k for k = 0 .. kmax: the number of related pairs with exactly k elements between them.
  std::vector<std::size_t> interval_abundances(std::size_t kmax) const {
    std::vector<std::size_t> N(kmax + 1, 0);
    for (std::size_t i = 0; i < n_; i++)
      for (std::size_t j = i + 1; j < n_; j++)
        if (precedes(i, j)) {
          const std::size_t k = interval_size(i, j);
          if (k <= kmax) N[k]++;
        }
    return N;
  }
  // The number of elements of the longest chain x_1 ≺ x_2 ≺ ... ≺ x_L.
  std::size_t longest_chain() const {
    std::vector<std::size_t> L(n_, 1);
    std::size_t best = n_ ? 1 : 0;
    for (std::size_t j = 0; j < n_; j++) {
      for (std::size_t i = 0; i < j; i++)
        if (precedes(i, j)) L[j] = std::max(L[j], L[i] + 1);
      best = std::max(best, L[j]);
    }
    return best;
  }

 private:
  bool empty_between(std::size_t i, std::size_t j) const {
    for (std::size_t w = 0; w < words_; w++)
      if (future_[i * words_ + w] & past_[j * words_ + w]) return false;
    return true;
  }

  std::vector<Event<d>> events_;
  std::size_t n_ = 0, words_ = 0, relations_ = 0;
  std::vector<std::uint64_t> future_, past_;
};

}  // namespace physis::spacetime
