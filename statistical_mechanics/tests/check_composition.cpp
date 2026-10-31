// Checks of the composition of open Markov jump processes (composition.hpp).
//
// - Two systems with no constraint and no joint events: the composite rates are the Kronecker sum
//   Q₁ ⊗ I + I ⊗ Q₂, entry by entry.
// - A fibre product: a switch and a counter that must agree on a bit; the joint states those that agree; the
//   switch flipping alone refused as an illegal lift, a joint event flipping both admitted.
// - Coherence: three systems (A, B, C; A–B and B–C interfaces) with local events and joint mechanisms on {A, B},
//   {B, C} and {A, B, C}, composed as (A ⊙ B) ⊙ C and as A ⊙ (B ⊙ C): the same joint states and the same rate
//   for every mechanism between every two joint states, so the same generator; each joint mechanism added
//   once (its total rate the same as when declared on the full product directly).

#include <cmath>
#include <map>
#include <string>
#include <vector>

#include <calx/verification/check.hpp>
#include <physis/statistical_mechanics/composition.hpp>
#include <physis/statistical_mechanics/master_equation.hpp>

using calx::verification::Suite;
namespace sm = physis::statistical_mechanics;

int main() {
  Suite s;
  const auto always = [](int, int) { return true; };
  const auto equal = [](int a, int b) { return a == b; };

  {
    const auto A = sm::atomic("A", 2, {{"p", {0, 0}}}, {{"flip", 0, 1, 1.5}, {"flip", 1, 0, 0.5}});
    const auto B = sm::atomic("B", 3, {{"q", {0, 0, 0}}}, {{"up", 0, 1, 2.0}, {"up", 1, 2, 3.0}, {"down", 2, 0, 0.25}});
    const auto C = sm::compose(A, B, "p", "q", always);
    const auto F = sm::compose(A, B, "", "", always);
    const auto k = C.system.rates(), ka = A.rates(), kb = B.rates();
    double worst = 0;
    for (std::size_t i = 0; i < C.system.states.size(); i++)
      for (std::size_t j = 0; j < C.system.states.size(); j++) {
        if (i == j) continue;
        const auto& x = C.system.states[i];
        const auto& y = C.system.states[j];
        const std::size_t a = x.at("A"), b = x.at("B"), a2 = y.at("A"), b2 = y.at("B");
        const double want = (b == b2 ? ka(a, a2) : 0) + (a == a2 ? kb(b, b2) : 0);
        worst = std::max(worst, std::fabs(k(i, j) - want));
      }
    s.require("no constraint: every pair of states is joint", C.system.states.size() == 6 && C.illegal.empty());
    s.check("the composite rates are the Kronecker sum", worst, 1e-15);
    s.require("juxtaposed with no ports named: the same rates, and every port kept", F.system.rates_by_components() == C.system.rates_by_components() && F.system.ports.size() == 2 && C.system.ports.empty());
  }
  {
    const auto Sw = sm::atomic("switch", 2, {{"bit", {0, 1}}}, {{"flip", 0, 1, 1.0}, {"flip", 1, 0, 1.0}});
    const auto Ct = sm::atomic("counter", 4, {{"bit", {0, 1, 0, 1}}}, {{"tick", 0, 2, 1.0}, {"tick", 2, 0, 1.0}, {"tick", 1, 3, 1.0}, {"tick", 3, 1, 1.0}});
    sm::JointMechanism both{"both", {"switch", "counter"}, [](const sm::ComponentState& z) {
      std::vector<std::pair<sm::ComponentState, double>> out;
      out.push_back({{{"switch", 1 - z.at("switch")}, {"counter", z.at("counter") ^ 1}}, 0.5});
      return out;
    }};
    const auto C = sm::compose(Sw, Ct, "bit", "bit", equal, {both});
    s.require("a fibre product: the joint states are those that agree on the bit", C.system.states.size() == 4);
    bool refused = !C.illegal.empty();
    for (const auto& m : C.illegal) refused = refused && m == "switch/flip";
    s.require("the switch flipping alone breaks the interface: refused", refused);
    std::size_t joint = 0, ticks = 0;
    for (const auto& m : C.system.moves) joint += m.mechanism == "both", ticks += m.mechanism == "counter/tick";
    s.require("the joint flip is admitted, and the counter's ticks, which keep the bit", joint == 4 && ticks == 4);
    s.require("the composite is irreducible", sm::irreducible(C.system.rates()));
  }
  {
    // A: 2 states, port to B its state. B: 3 states, port to A b % 2, port to C b ≥ 1. C: 2 states, port its state.
    const auto A = sm::atomic("A", 2, {{"toB", {0, 1}}}, {{"rest", 1, 0, 0.3}});
    const auto B = sm::atomic("B", 3, {{"toA", {0, 1, 0}}, {"toC", {0, 1, 1}}}, {{"climb", 1, 2, 0.8}, {"fall", 2, 0, 0.6}});
    const auto C = sm::atomic("C", 2, {{"toB", {0, 1}}}, {});
    const auto ab_ok = [](int a, int b) { return a <= b; };  // A may be 1 only where B's bit is 1
    const auto bc_ok = [](int b, int c) { return c <= b; };
    const sm::JointMechanism jab{"ab", {"A", "B"}, [](const sm::ComponentState& z) {
      std::vector<std::pair<sm::ComponentState, double>> out;
      if (z.at("A") == 0 && z.at("B") == 0) out.push_back({{{"A", 1}, {"B", 1}}, 1.1});
      return out;
    }};
    const sm::JointMechanism jbc{"bc", {"B", "C"}, [](const sm::ComponentState& z) {
      std::vector<std::pair<sm::ComponentState, double>> out;
      if (z.at("B") >= 1 && z.at("C") == 0) out.push_back({{{"C", 1}}, 0.9});
      if (z.at("C") == 1) out.push_back({{{"C", 0}}, 0.4});
      return out;
    }};
    const sm::JointMechanism jabc{"abc", {"A", "B", "C"}, [](const sm::ComponentState& z) {
      std::vector<std::pair<sm::ComponentState, double>> out;
      if (z.at("A") == 1 && z.at("C") == 1) out.push_back({{{"A", 0}, {"B", 0}, {"C", 0}}, 0.7});
      return out;
    }};
    const std::vector<sm::JointMechanism> J{jab, jbc, jabc};
    const auto left = sm::compose(sm::compose(A, B, "toB", "toA", ab_ok, J).system, C, "toC", "toB", bc_ok, J);
    const auto right = sm::compose(A, sm::compose(B, C, "toC", "toB", bc_ok, J).system, "toB", "toA", ab_ok, J);
    std::set<sm::ComponentState> sl(left.system.states.begin(), left.system.states.end()), sr(right.system.states.begin(), right.system.states.end());
    s.require("both orders have the same joint states", sl == sr && !sl.empty(), std::to_string(sl.size()) + " and " + std::to_string(sr.size()));
    const auto rl = left.system.rates_by_components(), rr = right.system.rates_by_components();
    s.require("and the same rate for every mechanism between every two joint states", rl == rr, std::to_string(rl.size()) + " and " + std::to_string(rr.size()) + " entries");
    std::map<std::string, double> total;
    for (const auto& [key, r] : rl) total[std::get<0>(key)] += r;
    s.require("each joint mechanism allocated once", total.count("ab") && total.count("bc") && total.count("abc") && total.count("A/rest") && total.count("B/climb"));
    // Declared directly on the product of the three, the joint mechanisms have the same rates.
    double ab_direct = 0;
    for (const auto& z : sl)
      if (z.at("A") == 0 && z.at("B") == 0) {
        auto n = z;
        n["A"] = 1, n["B"] = 1;
        if (sl.count(n)) ab_direct += 1.1;
      }
    s.check("the A–B mechanism's total rate is what it is on the joint states directly", total["ab"] - ab_direct, 1e-12);
  }
  return s.done();
}
