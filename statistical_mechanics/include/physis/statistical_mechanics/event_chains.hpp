// The probability of a record of events: chains of states on cells, joined at events (decision 0020; the
// infrastructure's ADR 0067 and its design, Definition 4.1 and Theorem 4.2).
//
// The model. A record has cells and pieces.
//   A cell (a segment of a thing's worldline, a lens of a relation) has a finite state space and a propagator P, the
//   Markov kernel over its life from its first state to its last: e^{τ(L − H)} for its generator L, its time τ and the
//   hazards H of what would have ended it (sub-stochastic where the record's events of a modelled kind did not happen).
//   A piece (an event, or one independent part of one) has the cells that end at it (in), the cells that begin at it
//   (out) and a kernel K(x_in, x_out): a transition, a law (a birth: no in), a sum (a death or the record's end: no out,
//   K = 1), a likelihood (an observation: in and out the same thing, K = diag ℓ), or a rule that joins several cells.
//   Pieces are given in the order of time; every cell begins at one piece and ends at a later one.
// Under the design's axiom A5, cells change independently between events, so the joint law of every cell's first and
// last states is the product of the kernels and propagators, and the probability of the record is its contraction:
//     Z = Σ_x Π_pieces K(x_in, x_out) Π_cells P_c(x_c^first, x_c^last).
//
// The contraction. Along the order of time, the live cells are held in factors, each a table over the last states of
// the cells it joins; cells come to share a factor only where a piece joins them. At a piece, the factors holding its
// in-cells are multiplied, the kernel summed against their in-states, and each out-cell carried to its last state by
// its propagator. A factor left with no live cell is a number of Z. The cost is the size of the largest factor, which
// is set by how many live cells events have joined (the design's §11); a factor beyond `limit` entries is refused.
// Tables are rescaled as they go, their scales summed into log Z, so long records do not underflow.
//
// Checks: tests/check_event_chains.cpp (Z against the enumeration of every state of every cell; a cell cut by an
// identity piece, its propagator split as P(τ₁)P(τ₂), the same Z (the design's Theorem 4.5); a record with stochastic
// kernels and no observation, Z = 1 (Proposition 4.6); factors joined only where pieces join their cells).

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <calx/linalg/dense.hpp>

namespace physis::statistical_mechanics {

struct ChainCell {
  std::size_t states = 0;
  calx::linalg::Matrix propagator;  // (first state, last state)
};

struct ChainPiece {
  std::vector<std::size_t> in, out;  // the cells ending here and beginning here
  std::vector<double> kernel;        // K over (in states…, out states…), row-major, the last index fastest
};

struct EventChains {
  std::vector<ChainCell> cells;
  std::vector<ChainPiece> pieces;  // in the order of time
};

struct ChainContraction {
  double log_z = 0;
  std::size_t largest = 1;  // entries of the largest factor met
};

// Every cell begun once and ended once, later; every kernel and propagator of its cells' sizes.
inline void validate(const EventChains& M) {
  const std::size_t C = M.cells.size();
  std::vector<long> begun(C, -1), ended(C, -1);
  for (std::size_t c = 0; c < C; c++)
    if (M.cells[c].states == 0 || M.cells[c].propagator.n != M.cells[c].states)
      throw std::invalid_argument("physis::statistical_mechanics::EventChains: cell " + std::to_string(c) + " has no states or a propagator of another size");
  for (std::size_t p = 0; p < M.pieces.size(); p++) {
    const auto& q = M.pieces[p];
    if (q.in.empty() && q.out.empty()) throw std::invalid_argument("physis::statistical_mechanics::EventChains: piece " + std::to_string(p) + " joins no cell");
    std::size_t size = 1;
    for (auto c : q.in) {
      if (c >= C) throw std::invalid_argument("physis::statistical_mechanics::EventChains: a cell out of range");
      if (ended[c] >= 0) throw std::invalid_argument("physis::statistical_mechanics::EventChains: cell " + std::to_string(c) + " ends twice");
      if (begun[c] < 0) throw std::invalid_argument("physis::statistical_mechanics::EventChains: cell " + std::to_string(c) + " ends before it begins");
      ended[c] = static_cast<long>(p), size *= M.cells[c].states;
    }
    for (auto c : q.out) {
      if (c >= C) throw std::invalid_argument("physis::statistical_mechanics::EventChains: a cell out of range");
      if (begun[c] >= 0) throw std::invalid_argument("physis::statistical_mechanics::EventChains: cell " + std::to_string(c) + " begins twice");
      begun[c] = static_cast<long>(p), size *= M.cells[c].states;
    }
    if (q.kernel.size() != size) throw std::invalid_argument("physis::statistical_mechanics::EventChains: piece " + std::to_string(p) + "'s kernel is not of its cells' sizes");
  }
  for (std::size_t c = 0; c < C; c++)
    if (begun[c] < 0 || ended[c] < 0) throw std::invalid_argument("physis::statistical_mechanics::EventChains: cell " + std::to_string(c) + " is not both begun and ended");
}

namespace detail {

struct Factor {
  std::vector<std::size_t> cells, dims;
  std::vector<double> table;  // row-major over cells, the last fastest
  std::size_t size() const { return table.size(); }
};

inline std::vector<std::size_t> strides(const std::vector<std::size_t>& dims) {
  std::vector<std::size_t> s(dims.size(), 1);
  for (std::size_t k = dims.size(); k-- > 1;) s[k - 1] = s[k] * dims[k];
  return s;
}

// The product of two factors, over the cells of the first and then those of the second it lacks.
inline Factor multiply(const Factor& A, const Factor& B, std::size_t limit) {
  Factor F;
  F.cells = A.cells, F.dims = A.dims;
  std::vector<long> inA(B.cells.size(), -1);
  for (std::size_t j = 0; j < B.cells.size(); j++) {
    const auto it = std::find(A.cells.begin(), A.cells.end(), B.cells[j]);
    if (it != A.cells.end()) inA[j] = it - A.cells.begin();
    else F.cells.push_back(B.cells[j]), F.dims.push_back(B.dims[j]);
  }
  std::size_t n = 1;
  for (auto d : F.dims) {
    n *= d;
    if (n > limit) throw std::length_error("physis::statistical_mechanics::contract: a factor beyond the limit of entries");
  }
  F.table.assign(n, 0.0);
  const auto sA = strides(A.dims), sB = strides(B.dims);
  std::vector<std::size_t> where(B.cells.size());  // the position of each of B's cells in F
  for (std::size_t j = 0, extra = A.cells.size(); j < B.cells.size(); j++) where[j] = inA[j] >= 0 ? static_cast<std::size_t>(inA[j]) : extra++;
  std::vector<std::size_t> idx(F.cells.size(), 0);
  for (std::size_t k = 0; k < n; k++) {
    std::size_t a = 0, b = 0;
    for (std::size_t i = 0; i < A.cells.size(); i++) a += idx[i] * sA[i];
    for (std::size_t j = 0; j < B.cells.size(); j++) b += idx[where[j]] * sB[j];
    F.table[k] = A.table[a] * B.table[b];
    for (std::size_t i = idx.size(); i-- > 0;) {
      if (++idx[i] < F.dims[i]) break;
      idx[i] = 0;
    }
  }
  return F;
}

}  // namespace detail

inline ChainContraction contract(const EventChains& M, std::size_t limit = std::size_t{1} << 22) {
  validate(M);
  using detail::Factor;
  ChainContraction R;
  std::vector<std::optional<Factor>> factors;
  std::vector<long> owner(M.cells.size(), -1);
  for (const auto& q : M.pieces) {
    // The factors holding the in-cells, multiplied.
    std::set<long> held;
    for (auto c : q.in) held.insert(owner[c]);
    Factor F;
    F.table = {1.0};
    for (long f : held) {
      F = detail::multiply(F, *factors[static_cast<std::size_t>(f)], limit);
      factors[static_cast<std::size_t>(f)].reset();
    }
    // The kernel summed against the in-states: the result over F's other cells, then the out-cells.
    std::vector<long> pos(q.in.size());
    for (std::size_t i = 0; i < q.in.size(); i++) pos[i] = std::find(F.cells.begin(), F.cells.end(), q.in[i]) - F.cells.begin();
    Factor G;
    std::vector<std::size_t> rest;
    for (std::size_t i = 0; i < F.cells.size(); i++)
      if (std::find(q.in.begin(), q.in.end(), F.cells[i]) == q.in.end()) rest.push_back(i), G.cells.push_back(F.cells[i]), G.dims.push_back(F.dims[i]);
    std::size_t outs = 1;
    for (auto c : q.out) G.cells.push_back(c), G.dims.push_back(M.cells[c].states), outs *= M.cells[c].states;
    std::size_t n = 1;
    for (auto d : G.dims) {
      n *= d;
      if (n > limit) throw std::length_error("physis::statistical_mechanics::contract: a factor beyond the limit of entries");
    }
    G.table.assign(n, 0.0);
    std::vector<std::size_t> idx(F.cells.size(), 0);
    for (std::size_t k = 0; k < F.table.size(); k++) {
      const double v = F.table[k];
      if (v != 0) {
        std::size_t r = 0, a = 0;
        for (std::size_t j = 0; j < rest.size(); j++) r = r * F.dims[rest[j]] + idx[rest[j]];
        for (std::size_t i = 0; i < q.in.size(); i++) a = a * M.cells[q.in[i]].states + idx[static_cast<std::size_t>(pos[i])];
        const double* K = q.kernel.data() + a * outs;
        double* out = G.table.data() + r * outs;
        for (std::size_t o = 0; o < outs; o++) out[o] += v * K[o];
      }
      for (std::size_t i = idx.size(); i-- > 0;) {
        if (++idx[i] < F.dims[i]) break;
        idx[i] = 0;
      }
    }
    // Each out-cell carried to its last state.
    for (std::size_t j = rest.size(); j < G.cells.size(); j++) {
      const auto& P = M.cells[G.cells[j]].propagator;
      const auto s = detail::strides(G.dims);
      const std::size_t d = G.dims[j], stride = s[j], outer = n / (d * stride);
      std::vector<double> col(d), res(d);
      for (std::size_t o = 0; o < outer; o++)
        for (std::size_t in = 0; in < stride; in++) {
          const std::size_t base = o * d * stride + in;
          for (std::size_t a = 0; a < d; a++) col[a] = G.table[base + a * stride];
          for (std::size_t e = 0; e < d; e++) {
            double t = 0;
            for (std::size_t a = 0; a < d; a++) t += col[a] * P(a, e);
            res[e] = t;
          }
          for (std::size_t e = 0; e < d; e++) G.table[base + e * stride] = res[e];
        }
    }
    R.largest = std::max(R.largest, G.size());
    // Rescaled, its scale into log Z.
    double big = 0;
    for (double v : G.table) big = std::max(big, std::fabs(v));
    if (!(big > 0)) {
      R.log_z = -std::numeric_limits<double>::infinity();
      return R;
    }
    for (double& v : G.table) v /= big;
    R.log_z += std::log(big);
    for (auto c : q.in) owner[c] = -1;
    if (G.cells.empty()) {
      R.log_z += std::log(G.table[0]);
      continue;
    }
    for (auto c : G.cells) owner[c] = static_cast<long>(factors.size());
    factors.push_back(std::move(G));
  }
  return R;
}

}  // namespace physis::statistical_mechanics
