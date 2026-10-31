// Markov jump processes on any state space, given by the occurrences of their events at each state: the
// stochastic simulation of a path and its record, anchored at the initial state, from which the path is
// reconstructed (decision 0018).
//
// The process. At a state z each mechanism m (an event kind) has its occurrences — its matches — at z, each
// with a rate λ_{m,o}(z) ≥ 0 and a law of successors. The backward generator is
//
//   (A f)(z) = Σ_m Σ_{o ∈ Match_m(z)} λ_{m,o}(z) ∫ [f(z′) − f(z)] P_{m,o}(z, dz′),
//
// so the total rate at z is the sum over every mechanism and every match: a rule that applies at several
// places fires at each. The direct method (Gillespie 1977, J. Phys. Chem. 81, 2340–2361) samples the path:
// wait an exponential time of rate Λ(z) = Σ λ, choose an occurrence with probability λ/Λ, draw its successor.
// When Λ(z) = 0 the state is absorbing and the path stays.
//
// The record. It holds the initial state and, for each event, its time, its mechanism, its match and the
// successor drawn: a list of event names alone would not determine the path when a rule has several matches
// or a random successor. replay() rebuilds the state at any time from the record alone.
//
// Explosion. The construction makes sense up to the first accumulation of jumps; with rates growing fast
// enough (a birth rate n², whose mean holding times sum to Σ 1/n² < ∞) that happens in finite time with
// probability one, and simulate() stops at a declared number of events and says so rather than pretending to
// reach the end. A Foster–Lyapunov bound A V ≤ c V with V norm-like excludes it (Meyn and Tweedie 1993,
// Adv. Appl. Probab. 25, 518–548, §2).
//
// Checks: tests/check_jump_process.cpp.

#pragma once

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <calx/random/stream.hpp>

namespace physis::statistical_mechanics {

template <class State>
struct Occurrence {
  std::string mechanism;  // the event kind
  std::string match;      // where it applies (a description of the match)
  double rate = 0;
  std::function<State(const State&, calx::random::Stream&)> successor;
};

template <class State>
struct Record {
  State initial;
  struct Event {
    double time = 0;
    std::string mechanism, match;
    State successor;
  };
  std::vector<Event> events;
  double end = 0;          // the time the record reaches
  bool stopped = false;    // the event limit was reached before `end` (a possible explosion)
  bool absorbed = false;   // a state with no occurrence was reached

  // The state at time t ≤ end, from the record alone.
  const State& at(double t) const {
    const State* z = &initial;
    for (const auto& e : events) {
      if (e.time > t) break;
      z = &e.successor;
    }
    return *z;
  }
};

// Samples a path from z0 until t_end, or until max_events events. occurrences(z) lists the occurrences at z.
template <class State, class Occurrences>
Record<State> simulate(State z0, Occurrences&& occurrences, double t_end, calx::random::Stream& rng, std::size_t max_events = 1000000) {
  if (!(t_end >= 0)) throw std::invalid_argument("physis::statistical_mechanics::simulate: an end time ≥ 0");
  Record<State> R;
  R.initial = z0;
  State z = std::move(z0);
  double t = 0;
  while (true) {
    const std::vector<Occurrence<State>> occ = occurrences(z);
    double total = 0;
    for (const auto& o : occ) {
      if (!(o.rate >= 0) || !std::isfinite(o.rate)) throw std::invalid_argument("physis::statistical_mechanics::simulate: an occurrence's rate is not finite and ≥ 0");
      total += o.rate;
    }
    if (total == 0) { R.absorbed = true; R.end = t_end; return R; }
    if (R.events.size() >= max_events) { R.stopped = true; R.end = t; return R; }
    const double dt = rng.exponential() / total;
    if (t + dt > t_end) { R.end = t_end; return R; }
    t += dt;
    double u = rng.uniform() * total;
    std::size_t k = 0;
    for (; k + 1 < occ.size(); k++) {
      if (u < occ[k].rate) break;
      u -= occ[k].rate;
    }
    while (occ[k].rate == 0 && k > 0) k--;
    z = occ[k].successor(z, rng);
    R.events.push_back({t, occ[k].mechanism, occ[k].match, z});
  }
}

}  // namespace physis::statistical_mechanics
