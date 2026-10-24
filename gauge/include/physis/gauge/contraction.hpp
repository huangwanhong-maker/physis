// The ℤ₂ lattice gauge theory of a 2-complex, computed exactly by contraction.
//
// The theory. On a 2-complex X a ℤ₂ gauge field is a sign u_e = ±1 on each edge, and the holonomy
// of a 2-cell is the product of the signs along its attaching word; over ℤ₂ an edge met twice
// cancels, so h_f(u) = Π u_e over the edges occurring an odd number of times in the word of f.
// The weight of a field is
//
//   W(u) = Π_f e^{β_f h_f(u)},
//
// the Wilson weight of field.hpp with a coupling for each 2-cell (for ℤ₂, Re tr U = ±1, and
// e^{−S} = e^{−Σ β_f} W). Some links may be fixed: a gauge fixing (a maximal tree, or the temporal
// gauge of a history), a boundary condition, or an observation. Z = Σ_u W(u) sums the free links
// with the counting measure; with Haar measures of mass 1 (exact.hpp) log Z is less by
// E_free log 2.
//
// The method. Taking the 2-cells in a given order, a link enters the contraction before its first
// 2-cell and leaves after its last; at each point of the order the partial sum of W over the links
// that have left is a function of the live links, the forward state F, and the partial sum over the
// 2-cells still to come is the backward state B. States are kept in the character basis of
// ℤ₂^live, f = Σ_S f̂_S χ_S, χ_S(u) = Π_{e ∈ S} u_e (Drouffe and Zuber 1983, Phys. Rep. 102, 1,
// §2; for ℤ₂, c₀ = cosh β and c₁ = sinh β, expansion_zn), where every step is elementary:
//
//   a 2-cell with weight e^{β c χ_V} (V its free links, c the product of its fixed ones) multiplies
//     by cosh β + c sinh β χ_V: f̂_S ← cosh β f̂_S + c sinh β f̂_{S ⊕ V};
//   summing a link x out keeps the S without x, doubled; a new link changes nothing.
//
// With no link fixed, every state is gauge invariant, and its characters are the ℤ/2 1-cycles of
// the cut: the graph whose edges are the live links and whose vertices are the vertices of X
// they join (gauge invariance at each vertex of the cut forces an even degree there). The state is a
// superposition of ℤ₂ spin networks on the cut, of dimension 2^{b₁ of the cut}; a forward state
// grows as the order crosses more independent cycles, so the cost is set by the widest cut, as in
// any transfer matrix (Baxter 1982, Exactly Solved Models, ch. 7; the tree-width of the dual
// order). Fixed links are constants, and conditioning on them is free.
//
// The readings. Every expectation is a pairing at one point of the order,
// Σ_u F(u) B(u) = 2^L Σ_S F̂_S B̂_S: a 2-cell's mean holonomy where it is applied, a link's mean where
// it leaves, a two-point function by carrying the insertion χ_x from the first link to the second,
// and at a cut the joint law of the live links (dense, when few) or the weight with which the
// character expansion's surfaces cross the cut along each spin network S, w_S = F̂_S B̂_S / Σ F̂ B̂.
// Forward states are kept at checkpoints and recomputed segment by segment for the backward pass.
//
// The derivatives. log Z and its gradient and Hessian in the couplings, grouped into parameters,
// come from one forward pass whose coefficients are calx jets (Taylor<N, 2>): the gradient is
// Σ ⟨h_f⟩ over a parameter's 2-cells and the Hessian the covariance of those sums, exact.
//
// A state larger than the limits throws Intractable: the caller samples instead (updates.hpp).

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <calx/differentiation/jet.hpp>
#include <morphe/topology/cw_complex.hpp>

namespace physis::gauge {

struct Intractable : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// A ℤ₂ gauge theory: its links, the links each 2-cell's holonomy multiplies, a coupling for each
// 2-cell, and the fixed links (0 free, ±1 fixed; empty: all free).
struct Z2Theory {
  std::size_t links = 0;
  std::vector<std::vector<std::size_t>> plaquettes;
  std::vector<double> beta;
  std::vector<int> fixed;

  int fixed_value(std::size_t e) const { return fixed.empty() ? 0 : fixed[e]; }

  // From a CW complex: a 2-cell multiplies the edges of its word that occur an odd number of times.
  static Z2Theory of(const morphe::topology::CWComplex& X, std::vector<double> beta, std::vector<int> fixed = {}) {
    if (beta.size() != X.faces.size()) throw std::invalid_argument("physis::gauge::Z2Theory::of: one coupling per 2-cell");
    Z2Theory T;
    T.links = X.edges.size();
    for (const auto& f : X.faces) {
      std::map<std::size_t, int> parity;
      for (const auto& s : f.boundary) parity[s.edge] ^= 1;
      std::vector<std::size_t> odd;
      for (const auto& [e, p] : parity)
        if (p) odd.push_back(e);
      T.plaquettes.push_back(std::move(odd));
    }
    T.beta = std::move(beta);
    T.fixed = std::move(fixed);
    return T;
  }
};

struct ContractionLimits {
  std::size_t support = std::size_t{1} << 20;  // characters in one state
  std::size_t segment = 32;                    // operations between checkpoints of the forward pass
  std::size_t dense = 22;                      // live links at a cut whose joint law is enumerated
};

// A cut after the first `after` 2-cells of the order. The links in `keep` stay live there; each set
// in `products` is read there, its mean product and the law of how many are negative.
struct ContractionCut {
  std::size_t after = 0;
  std::vector<std::size_t> keep;
  std::vector<std::vector<std::size_t>> products;
};

// The order of the contraction, with each link's slot (its bit in a character) while live.
struct ContractionSchedule {
  enum Kind : unsigned char { introduce, apply, eliminate, observe };
  struct Op {
    Kind kind;
    std::size_t index;  // a link (introduce, eliminate), a 2-cell (apply), a cut (observe)
  };
  std::vector<Op> ops;
  std::vector<int> slot;                          // link -> its bit while live; -1 for a fixed link
  std::vector<std::vector<std::size_t>> live_at;  // cut -> the links live there
  std::size_t width = 0;                          // the most links live at once
};

inline ContractionSchedule contraction_schedule(const Z2Theory& T, const std::vector<std::size_t>& order, const std::vector<ContractionCut>& cuts = {}) {
  ContractionSchedule S;
  const std::size_t n = T.links, F = order.size();
  if (F != T.plaquettes.size()) throw std::invalid_argument("physis::gauge::contraction_schedule: the order must take every 2-cell once");
  S.slot.assign(n, -1);
  std::vector<std::size_t> remaining(n, 0);
  for (std::size_t f : order)
    for (std::size_t e : T.plaquettes.at(f)) remaining[e]++;
  // Each link kept for a cut stays live until that cut; a kept link not yet used by then enters there.
  std::vector<long> keep_until(n, -1);
  std::vector<std::vector<std::size_t>> cuts_at(F + 1), deferred(F + 1), enter_at(F + 1);
  for (std::size_t q = 0; q < cuts.size(); q++) {
    if (cuts[q].after > F) throw std::invalid_argument("physis::gauge::contraction_schedule: a cut after more 2-cells than there are");
    cuts_at[cuts[q].after].push_back(q);
    for (std::size_t e : cuts[q].keep) keep_until[e] = std::max(keep_until[e], static_cast<long>(cuts[q].after));
  }
  std::uint64_t used = 0;
  std::vector<std::size_t> occupant(64, n);
  std::vector<char> entered(n, 0);
  auto enter = [&](std::size_t e) {
    if (entered[e] || T.fixed_value(e)) return;
    entered[e] = 1;
    // 63 at most: the character with every bit set marks an empty slot of the states' tables.
    if (__builtin_popcountll(used) >= 63) throw Intractable("more than 63 links live at one cut");
    const int b = __builtin_ctzll(~used);
    used |= std::uint64_t{1} << b;
    S.slot[e] = b;
    occupant[static_cast<std::size_t>(b)] = e;
    S.ops.push_back({ContractionSchedule::introduce, e});
    S.width = std::max<std::size_t>(S.width, static_cast<std::size_t>(__builtin_popcountll(used)));
  };
  auto leave = [&](std::size_t e) {
    S.ops.push_back({ContractionSchedule::eliminate, e});
    used &= ~(std::uint64_t{1} << S.slot[e]);
  };
  auto done_with = [&](std::size_t e, std::size_t k) {  // after k 2-cells
    if (T.fixed_value(e)) return;
    if (keep_until[e] >= static_cast<long>(k)) deferred[static_cast<std::size_t>(keep_until[e])].push_back(e);
    else leave(e);
  };
  auto observe = [&](std::size_t k) {
    for (std::size_t q : cuts_at[k]) {
      for (std::size_t e : cuts[q].keep)
        if (!entered[e]) {
          enter(e);
          if (remaining[e] == 0 && !T.fixed_value(e)) deferred[k].push_back(e);
        }
      S.ops.push_back({ContractionSchedule::observe, q});
      for (std::size_t b = 0; b < 64; b++)
        if ((used >> b) & 1) S.live_at[q].push_back(occupant[b]);
    }
    for (std::size_t e : deferred[k]) leave(e);
  };
  S.live_at.assign(cuts.size(), {});
  // Links in no 2-cell and kept by no cut are summed out at once (a factor 2 each).
  for (std::size_t e = 0; e < n; e++)
    if (remaining[e] == 0 && keep_until[e] < 0 && !T.fixed_value(e)) {
      enter(e);
      leave(e);
    }
  observe(0);
  for (std::size_t k = 0; k < F; k++) {
    const std::size_t f = order[k];
    for (std::size_t e : T.plaquettes[f]) enter(e);
    S.ops.push_back({ContractionSchedule::apply, f});
    for (std::size_t e : T.plaquettes[f])
      if (--remaining[e] == 0) done_with(e, k + 1);
    observe(k + 1);
  }
  return S;
}

// log cosh β, without overflow: log c₀(β), the free energy of one 2-cell of the ℤ₂ theory.
inline double log_cosh(double b) {
  const double a = std::fabs(b);
  return a + std::log1p(std::exp(-2 * a)) - std::log(2.0);
}

namespace contraction_detail {

// A mean of signs lies in [−1, 1]; a ratio of pairings may leave it by a rounding error.
inline double sign_mean(double x) { return std::clamp(x, -1.0, 1.0); }

// The coefficients of a state: an open-addressing table from character (a mask of live slots) to
// coefficient, of any value type (doubles, or calx jets for the derivatives). The mask with every
// bit set marks an empty slot.
template <class Value>
class Table {
 public:
  static constexpr std::uint64_t none = ~std::uint64_t{0};

  explicit Table(std::size_t expected = 4) { reserve(expected); }

  std::size_t size() const { return used_; }

  void add(std::uint64_t key, const Value& value) {
    if (2 * (used_ + 1) > keys_.size()) grow();
    const std::size_t m = keys_.size() - 1;
    std::size_t i = mix(key) & m;
    while (keys_[i] != none && keys_[i] != key) i = (i + 1) & m;
    if (keys_[i] == none) {
      keys_[i] = key;
      used_++;
    }
    values_[i] += value;
  }
  Value get(std::uint64_t key) const {
    const std::size_t m = keys_.size() - 1;
    for (std::size_t i = mix(key) & m; keys_[i] != none; i = (i + 1) & m)
      if (keys_[i] == key) return values_[i];
    return Value(0.0);
  }
  template <class Fn>
  void each(Fn fn) const {
    for (std::size_t i = 0; i < keys_.size(); i++)
      if (keys_[i] != none) fn(keys_[i], values_[i]);
  }
  template <class Fn>
  void each_value(Fn fn) {
    for (std::size_t i = 0; i < keys_.size(); i++)
      if (keys_[i] != none) fn(values_[i]);
  }

 private:
  std::vector<std::uint64_t> keys_;
  std::vector<Value> values_;
  std::size_t used_ = 0;

  static std::size_t mix(std::uint64_t x) {
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return static_cast<std::size_t>(x);
  }
  void reserve(std::size_t n) {
    std::size_t capacity = 8;
    while (capacity < 2 * n + 2) capacity <<= 1;
    keys_.assign(capacity, none);
    values_.assign(capacity, Value(0.0));
    used_ = 0;
  }
  void grow() {
    Table g(2 * used_ + 2);
    each([&](std::uint64_t k, const Value& v) { g.add(k, v); });
    *this = std::move(g);
  }
};

inline double magnitude(double v) { return std::fabs(v); }
template <std::size_t N, unsigned K>
inline double magnitude(const calx::differentiation::Taylor<N, K>& v) { return std::fabs(v.v); }

// f ← f · (a + b χ_V)
template <class Value>
inline Table<Value> times(const Table<Value>& f, const Value& a, const Value& b, std::uint64_t V) {
  Table<Value> g(V ? 2 * f.size() : f.size());
  if (V == 0) {
    const Value ab = a + b;
    f.each([&](std::uint64_t S, const Value& x) { g.add(S, ab * x); });
  } else {
    f.each([&](std::uint64_t S, const Value& x) {
      g.add(S, a * x);
      g.add(S ^ V, b * x);
    });
  }
  return g;
}

template <class Value>
inline Table<Value> summed_out(const Table<Value>& f, std::uint64_t x) {
  Table<Value> g(f.size());
  f.each([&](std::uint64_t S, const Value& a) {
    if (!(S & x)) g.add(S, a);
  });
  return g;
}

// A state with its scale: the function is e^{scale} Σ_S f̂_S χ_S.
template <class Value>
struct State {
  Table<Value> f;
  double scale = 0;

  State() { f.add(0, Value(1.0)); }

  void normalise() {
    double m = 0;
    f.each([&](std::uint64_t, const Value& a) { m = std::max(m, magnitude(a)); });
    if (m == 0) return;
    f.each_value([&](Value& a) { a *= 1 / m; });
    scale += std::log(m);
  }
  void sum_out(std::uint64_t x) {
    f = summed_out(f, x);
    scale += std::log(2.0);
  }
};

}  // namespace contraction_detail

// What the contraction yields.
struct Contraction {
  double log_z = 0;                    // log Σ_u W(u) over the free links
  std::vector<double> link;            // ⟨u_e⟩ (a fixed link: its value)
  std::vector<double> plaquette;       // ⟨h_f⟩
  std::vector<double> two_point;       // ⟨u_x u_y⟩ for each pair asked
  struct CutLaw {
    std::size_t live = 0, support = 0;
    std::vector<double> products;      // the mean product of each set read
    std::vector<double> negative;      // the law of how many products are negative, if dense
    double mean_negative = 0;          // its mean, always
    std::vector<std::pair<std::vector<std::size_t>, double>> networks;  // spin networks (sets of live links), weights
  };
  std::vector<CutLaw> cuts;
  std::size_t width = 0, largest = 0;  // the most live links, the largest state
};

// Contracts the theory along an order of its 2-cells. Two-points are pairs of free links, the
// first entering no later than the second leaves; `networks` asks for that many of the heaviest
// spin networks at each cut.
inline Contraction contract(const Z2Theory& T, const std::vector<std::size_t>& order,
                            const std::vector<std::array<std::size_t, 2>>& two_points = {},
                            const std::vector<ContractionCut>& cuts = {}, std::size_t networks = 0,
                            const ContractionLimits& limits = {}) {
  using namespace contraction_detail;
  using St = State<double>;
  using Tb = Table<double>;
  const ContractionSchedule P = contraction_schedule(T, order, cuts);
  const std::size_t n = T.links, N = P.ops.size(), SEG = std::max<std::size_t>(1, limits.segment);
  auto bit = [&](std::size_t e) { return P.slot[e] < 0 ? std::uint64_t{0} : std::uint64_t{1} << P.slot[e]; };
  auto cell_of = [&](std::size_t f, std::uint64_t& V, double& c) {
    V = 0;
    c = 1;
    for (std::size_t e : T.plaquettes[f]) {
      if (const int x = T.fixed_value(e)) c *= x;
      else V ^= bit(e);
    }
  };
  auto multiply = [&](St& s, std::size_t f) {
    std::uint64_t V;
    double c;
    cell_of(f, V, c);
    const double b = T.beta[f];
    s.f = times(s.f, 1.0, c * std::tanh(b), V);
    s.scale += log_cosh(b);
    s.normalise();
  };
  auto step = [&](St& s, std::size_t i) {
    const auto& op = P.ops[i];
    if (op.kind == ContractionSchedule::apply) multiply(s, op.index);
    else if (op.kind == ContractionSchedule::eliminate) s.sum_out(bit(op.index));
  };
  Contraction E;
  E.width = P.width;
  E.link.assign(n, 0);
  for (std::size_t e = 0; e < n; e++) E.link[e] = T.fixed_value(e);
  E.plaquette.assign(T.plaquettes.size(), 0);
  E.cuts.assign(cuts.size(), {});

  std::vector<std::size_t> enter_op(n, N), leave_op(n, N);
  for (std::size_t i = 0; i < N; i++) {
    if (P.ops[i].kind == ContractionSchedule::introduce) enter_op[P.ops[i].index] = i;
    if (P.ops[i].kind == ContractionSchedule::eliminate) leave_op[P.ops[i].index] = i;
  }
  std::set<std::size_t> starts, ends;
  for (const auto& [x, y] : two_points)
    if (x != y && P.slot[x] >= 0 && P.slot[y] >= 0) {
      starts.insert(enter_op[x] + 1);
      ends.insert(leave_op[y]);
    }

  // Forward, with a checkpoint every SEG operations and the states two-points start from.
  std::vector<St> checkpoint((N + SEG - 1) / SEG + 1);
  std::map<std::size_t, St> start_state;
  St s;
  for (std::size_t i = 0; i < N; i++) {
    if (i % SEG == 0) checkpoint[i / SEG] = s;
    if (starts.count(i)) start_state[i] = s;
    step(s, i);
    E.largest = std::max(E.largest, s.f.size());
    if (s.f.size() > limits.support) throw Intractable("a state of more than " + std::to_string(limits.support) + " characters");
  }
  E.log_z = std::log(s.f.get(0)) + s.scale;

  // Backward, segment by segment, the forward states of each segment recomputed from its checkpoint.
  std::map<std::size_t, std::pair<St, St>> end_state;  // op -> (F before it, B after it)
  St b;
  for (std::size_t seg = (N + SEG - 1) / SEG; seg-- > 0;) {
    const std::size_t lo = seg * SEG, hi = std::min(N, lo + SEG);
    std::vector<St> local;
    local.reserve(hi - lo);
    St x = checkpoint[seg];
    for (std::size_t i = lo; i < hi; i++) {
      local.push_back(x);
      step(x, i);
    }
    for (std::size_t i = hi; i-- > lo;) {
      const auto& op = P.ops[i];
      const Tb& f = local[i - lo].f;
      if (op.kind == ContractionSchedule::apply) {
        std::uint64_t V;
        double c;
        cell_of(op.index, V, c);
        const double t = std::tanh(T.beta[op.index]);
        double num = 0, den = 0;
        b.f.each([&](std::uint64_t S, double bv) {
          const double fs = f.get(S), fv = f.get(S ^ V);
          num += bv * (t * fs + c * fv);
          den += bv * (fs + c * t * fv);
        });
        E.plaquette[op.index] = sign_mean(num / den);
        multiply(b, op.index);
      } else if (op.kind == ContractionSchedule::eliminate) {
        const std::uint64_t xb = bit(op.index);
        if (ends.count(i)) end_state[i] = {local[i - lo], b};
        double num = 0, den = 0;
        b.f.each([&](std::uint64_t S, double bv) {
          num += bv * f.get(S ^ xb);
          den += bv * f.get(S);
        });
        E.link[op.index] = sign_mean(num / den);
      } else if (op.kind == ContractionSchedule::introduce) {
        b.sum_out(bit(op.index));
      } else {  // observe
        auto& law = E.cuts[op.index];
        const auto& live = P.live_at[op.index];
        const auto& sets = cuts[op.index].products;
        law.live = live.size();
        law.support = f.size();
        double den = 0;
        std::vector<std::pair<std::uint64_t, double>> w;
        b.f.each([&](std::uint64_t S, double bv) {
          const double fs = f.get(S);
          den += bv * fs;
          if (networks) w.push_back({S, bv * fs});
        });
        if (networks) {
          for (auto& [S, v] : w) v /= den;
          std::sort(w.begin(), w.end(), [](const auto& p, const auto& q) { return std::fabs(p.second) > std::fabs(q.second); });
          for (std::size_t k = 0; k < std::min(networks, w.size()); k++) {
            std::vector<std::size_t> set;
            for (std::size_t e : live)
              if (w[k].first & bit(e)) set.push_back(e);
            law.networks.push_back({std::move(set), w[k].second});
          }
        }
        for (const auto& set : sets) {
          std::uint64_t V = 0;
          double c = 1;
          for (std::size_t e : set) {
            if (const int v = T.fixed_value(e)) c *= v;
            else V ^= bit(e);
          }
          double num = 0;
          b.f.each([&](std::uint64_t S, double bv) { num += bv * f.get(S ^ V); });
          law.products.push_back(sign_mean(c * num / den));
          law.mean_negative += (1 - law.products.back()) / 2;
        }
        // The whole law, when the live links can be enumerated: from characters to values by the
        // Walsh–Hadamard transform (bit k set: the k-th live link is −1).
        if (!sets.empty() && live.size() <= limits.dense) {
          const std::size_t L = live.size(), size = std::size_t{1} << L;
          std::vector<double> fv(size, 0), bv(size, 0);
          auto compress = [&](std::uint64_t S) {
            std::size_t m = 0;
            for (std::size_t k = 0; k < L; k++)
              if (S & bit(live[k])) m |= std::size_t{1} << k;
            return m;
          };
          f.each([&](std::uint64_t S, double a) { fv[compress(S)] = a; });
          b.f.each([&](std::uint64_t S, double a) { bv[compress(S)] = a; });
          auto wht = [&](std::vector<double>& v) {
            for (std::size_t h = 1; h < size; h <<= 1)
              for (std::size_t i0 = 0; i0 < size; i0 += h << 1)
                for (std::size_t j = i0; j < i0 + h; j++) {
                  const double u = v[j], u2 = v[j + h];
                  v[j] = u + u2;
                  v[j + h] = u - u2;
                }
          };
          wht(fv);
          wht(bv);
          std::vector<long> position(n, -1);
          for (std::size_t k = 0; k < L; k++) position[live[k]] = static_cast<long>(k);
          law.negative.assign(sets.size() + 1, 0);
          double total = 0;
          for (std::size_t xv = 0; xv < size; xv++) {
            const double p = fv[xv] * bv[xv];
            if (p <= 0) continue;
            std::size_t count = 0;
            for (const auto& set : sets) {
              int product = 1;
              for (std::size_t e : set) {
                const int v = T.fixed_value(e);
                product *= v ? v : ((xv >> position[e]) & 1 ? -1 : 1);
              }
              count += product < 0;
            }
            law.negative[count] += p;
            total += p;
          }
          for (double& c : law.negative) c /= total;
        }
      }
    }
  }

  // Two-point functions by insertion: χ_x carried from just after x enters to where y leaves.
  for (const auto& [x, y] : two_points) {
    if (x == y) { E.two_point.push_back(1); continue; }
    if (P.slot[x] < 0 || P.slot[y] < 0) { E.two_point.push_back(E.link[x] * E.link[y]); continue; }
    const std::size_t i1 = enter_op[x] + 1, i2 = leave_op[y];
    if (i1 > i2) throw std::invalid_argument("physis::gauge::contract: a two-point whose second link leaves before the first enters");
    St g = start_state.at(i1);
    Tb moved(g.f.size());
    g.f.each([&](std::uint64_t S, double a) { moved.add(S ^ bit(x), a); });
    g.f = std::move(moved);
    for (std::size_t i = i1; i < i2; i++) step(g, i);
    const auto& [F2, B2] = end_state.at(i2);
    const std::uint64_t y2 = bit(y);
    double num = 0, den = 0;
    B2.f.each([&](std::uint64_t S, double bv) {
      num += bv * g.f.get(S ^ y2);
      den += bv * F2.f.get(S);
    });
    E.two_point.push_back(sign_mean(std::exp(g.scale - F2.scale) * num / den));
  }
  return E;
}

// ---- log Z and its derivatives in the couplings --------------------------------------------------

// log Z with its gradient and Hessian in M parameters, the coupling of 2-cell f being
// parameters[of[f]] (T.beta is not read), by one forward pass over calx jets. With order 1 the
// Hessian is left zero.
template <std::size_t M>
struct LogZDerivatives {
  double log_z = 0;
  std::array<double, M> gradient{};
  std::array<std::array<double, M>, M> hessian{};
  std::size_t largest = 0;
};

template <std::size_t M, unsigned K = 2>
inline LogZDerivatives<M> log_z_derivatives(const Z2Theory& T, const std::vector<std::size_t>& order, const std::vector<std::size_t>& of,
                                            const std::array<double, M>& parameters, const ContractionLimits& limits = {}) {
  using namespace contraction_detail;
  using J = calx::differentiation::Taylor<M, K>;
  const ContractionSchedule P = contraction_schedule(T, order);
  auto bit = [&](std::size_t e) { return P.slot[e] < 0 ? std::uint64_t{0} : std::uint64_t{1} << P.slot[e]; };
  const auto beta = calx::differentiation::variables<M, K>(parameters);
  std::array<J, M> ch, sh;
  for (std::size_t p = 0; p < M; p++) {
    ch[p] = cosh(beta[p]);
    sh[p] = sinh(beta[p]);
  }
  State<J> s;
  LogZDerivatives<M> out;
  for (const auto& op : P.ops) {
    if (op.kind == ContractionSchedule::apply) {
      std::uint64_t V = 0;
      double c = 1;
      for (std::size_t e : T.plaquettes[op.index]) {
        if (const int x = T.fixed_value(e)) c *= x;
        else V ^= bit(e);
      }
      const std::size_t p = of.at(op.index);
      s.f = times(s.f, ch[p], sh[p] * c, V);
      s.normalise();
    } else if (op.kind == ContractionSchedule::eliminate) {
      s.sum_out(bit(op.index));
    }
    out.largest = std::max(out.largest, s.f.size());
    if (s.f.size() > limits.support) throw Intractable("a state of more than " + std::to_string(limits.support) + " characters");
  }
  const J lz = log(s.f.get(0));
  out.log_z = lz.v + s.scale;
  for (std::size_t p = 0; p < M; p++) out.gradient[p] = lz.d[p];
  if constexpr (K == 2)
    for (std::size_t p = 0; p < M; p++)
      for (std::size_t q = 0; q < M; q++) out.hessian[p][q] = lz.h[p][q];
  return out;
}

// log Z alone, with the couplings of T, by one forward pass; `largest` receives the largest state.
inline double log_z(const Z2Theory& T, const std::vector<std::size_t>& order, const ContractionLimits& limits, std::size_t& largest) {
  using namespace contraction_detail;
  const ContractionSchedule P = contraction_schedule(T, order);
  auto bit = [&](std::size_t e) { return P.slot[e] < 0 ? std::uint64_t{0} : std::uint64_t{1} << P.slot[e]; };
  State<double> s;
  for (const auto& op : P.ops) {
    if (op.kind == ContractionSchedule::apply) {
      std::uint64_t V = 0;
      double c = 1;
      for (std::size_t e : T.plaquettes[op.index]) {
        if (const int x = T.fixed_value(e)) c *= x;
        else V ^= bit(e);
      }
      const double b = T.beta[op.index];
      s.f = times(s.f, 1.0, c * std::tanh(b), V);
      s.scale += log_cosh(b);
      s.normalise();
    } else if (op.kind == ContractionSchedule::eliminate) {
      s.sum_out(bit(op.index));
    }
    largest = std::max(largest, s.f.size());
    if (s.f.size() > limits.support) throw Intractable("a state of more than " + std::to_string(limits.support) + " characters");
  }
  return std::log(s.f.get(0)) + s.scale;
}

inline double log_z(const Z2Theory& T, const std::vector<std::size_t>& order, const ContractionLimits& limits = {}) {
  std::size_t largest = 0;
  return log_z(T, order, limits, largest);
}

}  // namespace physis::gauge
