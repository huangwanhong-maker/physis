// Reaction networks and their kinetics: stoichiometry, conservation and balance, mass action, and
// the continuous-time Markov chain of the numbers of molecules.
//
// A reaction network has species 1..s and reactions r: Σ_i α_ir X_i → Σ_i β_ir X_i, each with a
// rate constant k_r (Horn and Jackson 1972, Arch. Rational Mech. Anal. 47, 81; Feinberg 2019,
// Foundations of Chemical Reaction Network Theory, ch. 2). Its stoichiometric matrix is
// N = β − α ∈ ℤ^{s × m}, column r the change reaction r makes.
//
// Conservation. A vector y with yᵀN = 0 is conserved by every reaction: yᵀn is constant along
// every trajectory, deterministic or stochastic. The conserved y form the left kernel of N, a
// saturated lattice; conservation_laws() returns a basis of it over ℤ from the Smith normal form of
// Nᵀ (morphe). Given the composition of the species (a matrix A of atoms of each element in each
// species) and their charges z, the network balances elements and charge exactly when AN = 0 and
// zᵀN = 0: the rows of A and z must lie in that lattice.
//
// Mass action. In a well-mixed volume the concentrations obey dc/dt = N v(c) with
// v_r(c) = k_r Π_i c_i^{α_ir} (Guldberg and Waage 1864).
//
// The Markov chain (Gillespie 1976, J. Comput. Phys. 22, 403; 1977, J. Phys. Chem. 81, 2340;
// Anderson and Kurtz 2015, Stochastic Analysis of Biochemical Systems, ch. 1). With n_i molecules
// of X_i in a volume Ω, reaction r occurs at the rate
//
//   a_r(n) = c_r h_r(n),  h_r(n) = Π_i C(n_i, α_ir),  c_r = k_r Π_i α_ir! / Ω^{|α_r| − 1},
//
// h_r the number of distinct sets of reactant molecules, and takes n to n + N_r. The generator is
// (Lf)(n) = Σ_r a_r(n) [f(n + N_r) − f(n)]; with f the coordinate n_i it gives Dynkin's drift
// d/dt E[n] = N E[a(n)], which equals N v(E[n]/Ω)Ω only for networks of first order. As m injective
// matches of a pattern with |Aut| symmetries are m/|Aut| occurrences, h_r counts occurrences, not
// ordered tuples: the convention matters, and is the one of Gillespie's.
//
// The chain is solved exactly where its reachable states are finite (a positive conservation law
// bounds them): generator() writes it, and absorption() gives the law of the absorbing state from a
// start, by the embedded jump chain, exactly in rationals or in floating point, when the reachable
// states admit no cycle (each reaction strictly lowers some potential, as in irreversible
// networks). gillespie() samples paths with calx's counter-based streams.
//
// Symmetry. A permutation σ of the species is a symmetry of the network when it carries the
// reactions, with their constants, onto the reactions. Then a_{σr}(σn) = a_r(n) and the generator
// commutes with the permutation of states; a law invariant under σ stays invariant, though each
// path need not.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <functional>
#include <map>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <calx/random/stream.hpp>
#include <morphe/algebra/smith.hpp>
#include <morphe/arithmetic/rational.hpp>

namespace physis::chemistry {

using morphe::arithmetic::Integer;
using morphe::arithmetic::Rational;

struct Reaction {
  std::string name;
  std::vector<std::pair<std::size_t, unsigned>> reactants, products;  // (species, coefficient)
  Rational rate = 1;                                                    // k_r
  unsigned order() const {
    unsigned o = 0;
    for (const auto& [i, a] : reactants) o += a;
    return o;
  }
};

struct Network {
  std::vector<std::string> species;
  std::vector<Reaction> reactions;

  std::size_t find(const std::string& name) const {
    for (std::size_t i = 0; i < species.size(); i++)
      if (species[i] == name) return i;
    throw std::invalid_argument("physis::chemistry::Network: no species " + name);
  }
  // Adds a reaction written with species names: {{"A", 1}, {"L", 1}} → {{"L", 2}}.
  Reaction& add(std::string name, const std::vector<std::pair<std::string, unsigned>>& lhs,
                const std::vector<std::pair<std::string, unsigned>>& rhs, Rational k) {
    Reaction r;
    r.name = std::move(name);
    for (const auto& [s, a] : lhs) r.reactants.push_back({find(s), a});
    for (const auto& [s, b] : rhs) r.products.push_back({find(s), b});
    r.rate = std::move(k);
    reactions.push_back(std::move(r));
    return reactions.back();
  }
  // α_ir and β_ir.
  long long alpha(std::size_t i, std::size_t r) const {
    long long a = 0;
    for (const auto& [j, c] : reactions[r].reactants) a += j == i ? c : 0;
    return a;
  }
  long long beta(std::size_t i, std::size_t r) const {
    long long b = 0;
    for (const auto& [j, c] : reactions[r].products) b += j == i ? c : 0;
    return b;
  }
  long long change(std::size_t i, std::size_t r) const { return beta(i, r) - alpha(i, r); }
};

// N = β − α, species by reactions.
inline morphe::algebra::Matrix stoichiometric_matrix(const Network& net) {
  morphe::algebra::Matrix N(net.species.size(), net.reactions.size());
  for (std::size_t i = 0; i < net.species.size(); i++)
    for (std::size_t r = 0; r < net.reactions.size(); r++) N(i, r) = Integer(net.change(i, r));
  return N;
}

// A basis over ℤ of {y : yᵀN = 0}.
inline std::vector<std::vector<Integer>> conservation_laws(const Network& net) {
  const auto N = stoichiometric_matrix(net);
  const auto S = morphe::algebra::smith(N.transpose(), {.V = true});
  std::vector<std::vector<Integer>> out;
  for (std::size_t j = S.rank(); j < net.species.size(); j++) {
    std::vector<Integer> y(net.species.size());
    for (std::size_t i = 0; i < net.species.size(); i++) y[i] = S.V(i, j);
    out.push_back(std::move(y));
  }
  return out;
}

// A N and zᵀ N: the imbalance of each element and of charge in each reaction (zero when balanced).
struct Balance {
  morphe::algebra::Matrix elements;  // elements × reactions
  std::vector<Integer> charge;       // per reaction
  bool balanced() const { return elements.is_zero() && std::all_of(charge.begin(), charge.end(), [](const Integer& c) { return c.is_zero(); }); }
};
inline Balance balance(const Network& net, const morphe::algebra::Matrix& composition, const std::vector<Integer>& charges) {
  if (composition.cols != net.species.size() || charges.size() != net.species.size())
    throw std::invalid_argument("physis::chemistry::balance: composition or charges of the wrong size");
  const auto N = stoichiometric_matrix(net);
  Balance b;
  b.elements = composition * N;
  b.charge.assign(net.reactions.size(), Integer(0));
  for (std::size_t r = 0; r < net.reactions.size(); r++)
    for (std::size_t i = 0; i < net.species.size(); i++) b.charge[r] += charges[i] * N(i, r);
  return b;
}

// dc/dt = N v(c) by mass action, with the rate constants as doubles.
template <class Vector>
void mass_action(const Network& net, const Vector& c, Vector& dc) {
  for (std::size_t i = 0; i < net.species.size(); i++) dc[i] = 0;
  for (std::size_t r = 0; r < net.reactions.size(); r++) {
    double v = net.reactions[r].rate.to_double();
    for (const auto& [i, a] : net.reactions[r].reactants) v *= std::pow(c[i], static_cast<double>(a));
    for (const auto& [i, a] : net.reactions[r].reactants) dc[i] -= a * v;
    for (const auto& [i, b] : net.reactions[r].products) dc[i] += b * v;
  }
}

// ---- the Markov chain of molecule numbers ------------------------------------------------------------------------------------

using State = std::vector<long long>;

// c_r = k_r Π α! / Ω^{|α| − 1}, exactly.
inline Rational stochastic_constant(const Reaction& r, long long volume) {
  Rational c = r.rate;
  for (const auto& [i, a] : r.reactants)
    for (unsigned k = 2; k <= a; k++) c = c * Rational(static_cast<long long>(k));
  const unsigned order = r.order();
  if (order == 0) return c * Rational(volume);
  for (unsigned k = 1; k < order; k++) c = c / Rational(volume);
  return c;
}

// a_r(n) = c_r Π C(n_i, α_ir), in the number type T (Rational or double).
template <class T = Rational>
T propensity(const Network& net, std::size_t r, const State& n, long long volume = 1) {
  T h = T(1);
  for (const auto& [i, a] : net.reactions[r].reactants) {
    if (n[i] < static_cast<long long>(a)) return T(0);
    // C(n, a) = n (n−1) … (n−a+1) / a!
    long long num = 1, den = 1;
    for (unsigned k = 0; k < a; k++) num *= n[i] - k, den *= k + 1;
    h = h * T(num / den);
  }
  if constexpr (std::is_same_v<T, Rational>) return stochastic_constant(net.reactions[r], volume) * h;
  else return stochastic_constant(net.reactions[r], volume).to_double() * h;
}

inline State fire(const Network& net, std::size_t r, State n) {
  for (std::size_t i = 0; i < n.size(); i++) n[i] += net.change(i, r);
  return n;
}

// The drift Σ_r a_r(n) N_r of the generator at n.
template <class T = Rational>
std::vector<T> drift(const Network& net, const State& n, long long volume = 1) {
  std::vector<T> d(net.species.size(), T(0));
  for (std::size_t r = 0; r < net.reactions.size(); r++) {
    const T a = propensity<T>(net, r, n, volume);
    for (std::size_t i = 0; i < d.size(); i++) d[i] = d[i] + a * T(net.change(i, r));
  }
  return d;
}

struct StateSpace {
  std::vector<State> states;           // the start first
  std::map<State, std::size_t> index;
  std::size_t find(const State& n) const {
    const auto it = index.find(n);
    return it == index.end() ? states.size() : it->second;
  }
};

// The states reachable from n0, breadth first (throws past `limit`).
inline StateSpace reachable(const Network& net, const State& n0, std::size_t limit = 1000000) {
  StateSpace S;
  S.states.push_back(n0);
  S.index[n0] = 0;
  for (std::size_t k = 0; k < S.states.size(); k++) {
    const State n = S.states[k];
    for (std::size_t r = 0; r < net.reactions.size(); r++) {
      if (net.reactions[r].rate == Rational(0) || propensity<Rational>(net, r, n) == Rational(0)) continue;
      const State m = fire(net, r, n);
      if (S.index.count(m)) continue;
      if (S.states.size() >= limit) throw std::runtime_error("physis::chemistry::reachable: more states than the limit");
      S.index[m] = S.states.size();
      S.states.push_back(m);
    }
  }
  return S;
}

// The generator on the reachable states: Q[x] = {(y, rate)} for y ≠ x; the diagonal is minus the
// row's sum.
template <class T = Rational>
std::vector<std::map<std::size_t, T>> generator(const Network& net, const StateSpace& S, long long volume = 1) {
  std::vector<std::map<std::size_t, T>> Q(S.states.size());
  for (std::size_t x = 0; x < S.states.size(); x++)
    for (std::size_t r = 0; r < net.reactions.size(); r++) {
      const T a = propensity<T>(net, r, S.states[x], volume);
      if (a == T(0)) continue;
      const std::size_t y = S.find(fire(net, r, S.states[x]));
      if (y == S.states.size()) throw std::invalid_argument("physis::chemistry::generator: a transition leaves the state space");
      if (y != x) Q[x][y] = Q[x][y] + a;
    }
  return Q;
}

// The law of the absorbing state reached from the start (states where no reaction can occur), by
// the embedded jump chain: along a topological order when the reachable states admit no cycle,
// otherwise from the occupation measure g = e₀ᵀ(I − P_TT)⁻¹ of the transient states, solved by
// Gaussian elimination (dense: for a few thousand states at most). Every reachable state must lead
// to absorption.
template <class T = Rational>
std::map<State, T> absorption(const Network& net, const State& n0, long long volume = 1, std::size_t limit = 1000000) {
  const StateSpace S = reachable(net, n0, limit);
  const auto Q = generator<T>(net, S, volume);
  const std::size_t n = S.states.size();
  const auto absorbing = [&](std::size_t x) { return Q[x].empty(); };
  // Kahn's topological order of the transition graph.
  std::vector<std::size_t> indeg(n, 0), order;
  for (std::size_t x = 0; x < n; x++)
    for (const auto& [y, a] : Q[x]) indeg[y]++;
  std::deque<std::size_t> ready;
  for (std::size_t x = 0; x < n; x++)
    if (!indeg[x]) ready.push_back(x);
  while (!ready.empty()) {
    const std::size_t x = ready.front();
    ready.pop_front();
    order.push_back(x);
    for (const auto& [y, a] : Q[x])
      if (!--indeg[y]) ready.push_back(y);
  }
  if (order.size() != n) {
    // Transient states, and the system (I − P_TT)ᵀ g = e₀ for the expected visits g.
    std::vector<std::size_t> transient, slot(n, n);
    for (std::size_t x = 0; x < n; x++)
      if (!absorbing(x)) slot[x] = transient.size(), transient.push_back(x);
    const std::size_t m = transient.size();
    std::vector<std::vector<T>> M(m, std::vector<T>(m + 1, T(0)));
    for (std::size_t a = 0; a < m; a++) M[a][a] = T(1);
    std::vector<T> total(n, T(0));
    for (std::size_t x = 0; x < n; x++)
      for (const auto& [y, a] : Q[x]) total[x] = total[x] + a;
    for (std::size_t a = 0; a < m; a++)
      for (const auto& [y, rate] : Q[transient[a]])
        if (slot[y] < n) M[slot[y]][a] = M[slot[y]][a] - rate / total[transient[a]];
    if (slot[0] < n) M[slot[0]][m] = T(1);
    for (std::size_t col = 0; col < m; col++) {
      std::size_t pivot = col;
      if constexpr (std::is_floating_point_v<T>) {
        for (std::size_t r = col + 1; r < m; r++)
          if (std::fabs(M[r][col]) > std::fabs(M[pivot][col])) pivot = r;
      } else {
        while (pivot < m && M[pivot][col] == T(0)) pivot++;
      }
      if (pivot == m || M[pivot][col] == T(0)) throw std::invalid_argument("physis::chemistry::absorption: some states never reach absorption");
      std::swap(M[col], M[pivot]);
      const T inv = T(1) / M[col][col];
      for (std::size_t k = col; k <= m; k++) M[col][k] = M[col][k] * inv;
      for (std::size_t r = 0; r < m; r++) {
        if (r == col || M[r][col] == T(0)) continue;
        const T f = M[r][col];
        for (std::size_t k = col; k <= m; k++) M[r][k] = M[r][k] - f * M[col][k];
      }
    }
    std::map<State, T> out;
    if (slot[0] == n) {
      out[S.states[0]] = T(1);
      return out;
    }
    for (std::size_t a = 0; a < m; a++) {
      const T g = M[a][m];
      if (g == T(0)) continue;
      for (const auto& [y, rate] : Q[transient[a]])
        if (slot[y] == n) out[S.states[y]] = out[S.states[y]] + g * rate / total[transient[a]];
    }
    return out;
  }
  std::vector<T> p(n, T(0));
  p[0] = T(1);
  std::map<State, T> out;
  for (auto x : order) {
    if (p[x] == T(0)) continue;
    T total = T(0);
    for (const auto& [y, a] : Q[x]) total = total + a;
    if (total == T(0)) {
      out[S.states[x]] = p[x];
      continue;
    }
    for (const auto& [y, a] : Q[x]) p[y] = p[y] + p[x] * a / total;
  }
  return out;
}

// A path of the chain from n0 until t_end or absorption, by Gillespie's direct method; at most
// `max_events` events. The times and states are appended to the vectors.
inline void gillespie(const Network& net, State n, double t_end, calx::random::Stream& rng, std::vector<double>& times,
                      std::vector<State>& states, long long volume = 1, std::size_t max_events = 10000000) {
  double t = 0;
  times.push_back(t);
  states.push_back(n);
  std::vector<double> a(net.reactions.size());
  for (std::size_t event = 0; event < max_events; event++) {
    double total = 0;
    for (std::size_t r = 0; r < a.size(); r++) total += a[r] = propensity<double>(net, r, n, volume);
    if (total <= 0) return;
    t += rng.exponential() / total;
    if (t > t_end) return;
    double u = rng.uniform() * total;
    std::size_t r = 0;
    while (r + 1 < a.size() && u >= a[r]) u -= a[r++];
    n = fire(net, r, std::move(n));
    times.push_back(t);
    states.push_back(n);
  }
}

// ---- symmetry -----------------------------------------------------------------------------------------------------------------

using SpeciesPermutation = std::vector<std::size_t>;  // species i ↦ sigma[i]

inline State permute(const SpeciesPermutation& sigma, const State& n) {
  State m(n.size());
  for (std::size_t i = 0; i < n.size(); i++) m[sigma[i]] = n[i];
  return m;
}

// The reaction σ carries reaction r to, or reactions.size() if its image is not a reaction with the
// same constant.
inline std::size_t image(const Network& net, const SpeciesPermutation& sigma, std::size_t r) {
  const auto sorted = [&](const std::vector<std::pair<std::size_t, unsigned>>& side) {
    std::vector<std::pair<std::size_t, unsigned>> s;
    for (const auto& [i, a] : side) s.push_back({sigma[i], a});
    std::sort(s.begin(), s.end());
    return s;
  };
  const auto plain = [](std::vector<std::pair<std::size_t, unsigned>> s) {
    std::sort(s.begin(), s.end());
    return s;
  };
  const auto lhs = sorted(net.reactions[r].reactants), rhs = sorted(net.reactions[r].products);
  for (std::size_t q = 0; q < net.reactions.size(); q++)
    if (plain(net.reactions[q].reactants) == lhs && plain(net.reactions[q].products) == rhs && net.reactions[q].rate == net.reactions[r].rate) return q;
  return net.reactions.size();
}

// Whether σ is a symmetry of the network: it permutes the reactions, constants kept.
inline bool symmetry_of(const Network& net, const SpeciesPermutation& sigma) {
  std::vector<char> hit(net.reactions.size(), 0);
  for (std::size_t r = 0; r < net.reactions.size(); r++) {
    const std::size_t q = image(net, sigma, r);
    if (q == net.reactions.size() || hit[q]++) return false;
  }
  return true;
}

}  // namespace physis::chemistry
