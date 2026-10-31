// Composing open Markov jump processes at an interface: the joint states compatible there, each part's
// events lifted to them, joint events declared once on the parts they involve, and the composite rates; the
// order in which three systems are composed leaves the generator unchanged when every event is allocated
// once (decision 0018).
//
// An open system (in the spirit of Willems's interconnection of behaviours, 2007, IEEE Control Syst. Mag.
// 27(6), 46–99; in stochastic process algebra, Hillston 1996, A Compositional Approach to Performance
// Modelling; the stochastic automata networks of Plateau 1985, ACM SIGMETRICS, whose generators are Kronecker
// sums) has finitely many states, ports — each a map from its states to the values it shows there — and
// events: moves s → s′ at rates, each labelled with the mechanism it belongs to. Two systems with ports p and
// q and a compatibility relation R between their values compose on
//
//   Z = {(s, t) : (p(s), q(t)) ∈ R},   the fibre product Z₁ ×_B Z₂ when R is equality on B,
//
// a local event s → s′ lifting to (s, t) → (s′, t) when the successor is still compatible, and otherwise
// refused and reported as an illegal lift: one side alone cannot break the interface. A joint event (both
// sides moving at once, or one side's rate depending on the other's state) is a mechanism on the components
// it involves, declared once: it is added at the first composition whose result contains all of them. The
// composite rate between two joint states is the sum of the events that make that move; with no constraint
// and no joint events the generator is the Kronecker sum Q₁ ⊗ I + I ⊗ Q₂.
//
// This is not the composition of open Markov processes of Baez, Fong and Pollard (2016, J. Math. Phys. 57,
// 033301), where two processes are glued along shared states (a cospan of finite sets, the composite's state
// set a pushout). Here the state set is a product constrained at the interface: two parts that each have a
// state at every moment, like two stores and the street between them.
//
// Coherence. Composed in either order, three systems give the same joint states (each a state of every
// component) and, every mechanism being allocated once, the same rate for each mechanism between each two
// joint states: the generator is the same, composition associative up to the canonical bijection. Checked
// rather than assumed: rates_by_components() lists the rates keyed by the components' states, so the two
// orders can be compared. Limits: the joint states are enumerated, |Z₁| |Z₂| at most. Checks:
// tests/check_composition.cpp.

#pragma once

#include <cstddef>
#include <functional>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <calx/linalg/dense.hpp>

namespace physis::statistical_mechanics {

struct Move {
  std::string mechanism;
  std::size_t from = 0, to = 0;
  double rate = 0;
};

// A state of a composite, as the state of each component, by name.
using ComponentState = std::map<std::string, std::size_t>;

// A joint mechanism: on the components it names, from a joint state the successors (of those components; the
// others unchanged) and their rates.
struct JointMechanism {
  std::string name;
  std::set<std::string> components;
  std::function<std::vector<std::pair<ComponentState, double>>(const ComponentState&)> moves;
};

struct OpenSystem {
  std::string name;
  std::vector<ComponentState> states;            // each state as the states of the components
  std::map<std::string, std::vector<int>> ports;  // port ↦ the value it shows in each state
  std::vector<Move> moves;
  std::set<std::string> components() const {
    std::set<std::string> c;
    if (!states.empty())
      for (const auto& [k, v] : states[0]) c.insert(k);
    return c;
  }
  // The rates k_ij (summed over mechanisms) as a matrix.
  calx::linalg::Matrix rates() const {
    calx::linalg::Matrix k(states.size());
    for (const auto& m : moves)
      if (m.from != m.to) k(m.from, m.to) += m.rate;
    return k;
  }
  // The rates of one mechanism.
  calx::linalg::Matrix rates(const std::string& mechanism) const {
    calx::linalg::Matrix k(states.size());
    for (const auto& m : moves)
      if (m.mechanism == mechanism && m.from != m.to) k(m.from, m.to) += m.rate;
    return k;
  }
  std::size_t find(const ComponentState& s) const {
    for (std::size_t i = 0; i < states.size(); i++)
      if (states[i] == s) return i;
    return states.size();
  }
  // (mechanism, from, to) ↦ rate, with states written as the components' states: independent of how the
  // system was composed.
  std::map<std::tuple<std::string, ComponentState, ComponentState>, double> rates_by_components() const {
    std::map<std::tuple<std::string, ComponentState, ComponentState>, double> out;
    for (const auto& m : moves) out[{m.mechanism, states[m.from], states[m.to]}] += m.rate;
    return out;
  }
};

// An atomic system of n states, its own component.
inline OpenSystem atomic(const std::string& name, std::size_t n, std::map<std::string, std::vector<int>> ports, std::vector<Move> moves) {
  OpenSystem S;
  S.name = name;
  for (std::size_t i = 0; i < n; i++) S.states.push_back({{name, i}});
  for (const auto& [p, v] : ports)
    if (v.size() != n) throw std::invalid_argument("physis::statistical_mechanics::atomic: port " + p + " needs a value for each state");
  S.ports = std::move(ports);
  for (auto& m : moves) {
    if (m.from >= n || m.to >= n || !(m.rate >= 0)) throw std::invalid_argument("physis::statistical_mechanics::atomic: a move between states at a rate ≥ 0");
    if (m.mechanism.find('/') == std::string::npos) m.mechanism = name + "/" + m.mechanism;
  }
  S.moves = std::move(moves);
  return S;
}

struct Composition {
  OpenSystem system;
  std::vector<std::string> illegal;  // local moves refused because their successor breaks the interface
};

// a ⊙ b along port pa of a and port pb of b, joint states those with compatible(pa value, pb value); the ports
// used are consumed, the others kept. With both port names empty the two are juxtaposed with no constraint (every
// pair of states joint) and keep all their ports. Joint mechanisms are added when this composition is the first
// to contain all their components.
inline Composition compose(const OpenSystem& a, const OpenSystem& b, const std::string& pa, const std::string& pb,
                           const std::function<bool(int, int)>& compatible, const std::vector<JointMechanism>& joint = {}) {
  const bool free = pa.empty() && pb.empty();
  if (!free && (pa.empty() || pb.empty())) throw std::invalid_argument("physis::statistical_mechanics::compose: a port on both sides, or on neither");
  const std::vector<int> none_a(a.states.size(), 0), none_b(b.states.size(), 0);
  const auto ia = free ? a.ports.end() : a.ports.find(pa), ib = free ? b.ports.end() : b.ports.find(pb);
  if (!free && (ia == a.ports.end() || ib == b.ports.end())) throw std::invalid_argument("physis::statistical_mechanics::compose: no such port");
  const std::vector<int>& va = free ? none_a : ia->second;
  const std::vector<int>& vb = free ? none_b : ib->second;
  const auto ca = a.components(), cb = b.components();
  for (const auto& c : ca)
    if (cb.count(c)) throw std::invalid_argument("physis::statistical_mechanics::compose: the systems share the component " + c);
  Composition out;
  OpenSystem& S = out.system;
  S.name = "(" + a.name + " ⊙ " + b.name + ")";
  std::map<std::pair<std::size_t, std::size_t>, std::size_t> index;
  std::vector<std::pair<std::size_t, std::size_t>> pairs;
  for (std::size_t s = 0; s < a.states.size(); s++)
    for (std::size_t t = 0; t < b.states.size(); t++)
      if (free || compatible(va[s], vb[t])) {
        index[{s, t}] = pairs.size();
        pairs.push_back({s, t});
        ComponentState joint_state = a.states[s];
        joint_state.insert(b.states[t].begin(), b.states[t].end());
        S.states.push_back(std::move(joint_state));
      }
  for (const auto& [name, values] : a.ports)
    if (name != pa)
      for (const auto& [s, t] : pairs) S.ports[name].push_back(values[s]);
  for (const auto& [name, values] : b.ports)
    if (name != pb) {
      if (S.ports.count(name)) throw std::invalid_argument("physis::statistical_mechanics::compose: both systems have a port " + name);
      for (const auto& [s, t] : pairs) S.ports[name].push_back(values[t]);
    }
  // Lifts of the local moves.
  for (const auto& m : a.moves)
    for (std::size_t t = 0; t < b.states.size(); t++) {
      const auto f = index.find({m.from, t});
      if (f == index.end()) continue;
      const auto g = index.find({m.to, t});
      if (g == index.end()) { out.illegal.push_back(m.mechanism); continue; }
      S.moves.push_back({m.mechanism, f->second, g->second, m.rate});
    }
  for (const auto& m : b.moves)
    for (std::size_t s = 0; s < a.states.size(); s++) {
      const auto f = index.find({s, m.from});
      if (f == index.end()) continue;
      const auto g = index.find({s, m.to});
      if (g == index.end()) { out.illegal.push_back(m.mechanism); continue; }
      S.moves.push_back({m.mechanism, f->second, g->second, m.rate});
    }
  // Joint mechanisms whose components are all here for the first time.
  std::set<std::string> all = ca;
  all.insert(cb.begin(), cb.end());
  const auto within = [](const std::set<std::string>& need, const std::set<std::string>& have) {
    for (const auto& c : need)
      if (!have.count(c)) return false;
    return true;
  };
  for (const auto& J : joint) {
    if (!within(J.components, all) || within(J.components, ca) || within(J.components, cb)) continue;
    for (std::size_t i = 0; i < S.states.size(); i++) {
      ComponentState sub;
      for (const auto& c : J.components) sub[c] = S.states[i].at(c);
      for (const auto& [succ, rate] : J.moves(sub)) {
        if (!(rate >= 0)) throw std::invalid_argument("physis::statistical_mechanics::compose: a joint rate ≥ 0");
        if (rate == 0) continue;
        ComponentState next = S.states[i];
        for (const auto& [c, v] : succ) {
          if (!J.components.count(c)) throw std::invalid_argument("physis::statistical_mechanics::compose: joint mechanism " + J.name + " moves a component it does not name");
          next[c] = v;
        }
        const std::size_t j = S.find(next);
        if (j == S.states.size()) { out.illegal.push_back(J.name); continue; }
        if (j != i) S.moves.push_back({J.name, i, j, rate});
      }
    }
  }
  return out;
}

}  // namespace physis::statistical_mechanics
