// Checks of the simulation of jump processes and their records (jump_process.hpp).
//
// - A birth–death chain on 0…6: the time spent in each state over a long path within 0.01 of the stationary
//   law (master_equation.hpp's GTH), every event's successor a neighbour.
// - Matches: k switches each flipping at rate r — one mechanism with k matches — make k r T events on average.
// - The record determines the path: at(t) replays the state at every event time from the initial state alone.
// - Explosion: births at rate n² from n = 1 accumulate before Σ 1/n² = π²/6 on average; the simulation stops
//   at its event limit and says so, the time reached averaging π²/6 over 200 paths.

#include <cmath>
#include <numbers>
#include <string>
#include <vector>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <physis/statistical_mechanics/jump_process.hpp>
#include <physis/statistical_mechanics/master_equation.hpp>

using calx::verification::Suite;
namespace sm = physis::statistical_mechanics;

int main() {
  Suite s;
  calx::random::Stream rng(5, 0, 0, 1);
  {
    const int N = 6;
    const auto up = [](int n) { return 1.0 + 0.3 * n; };
    const auto down = [](int n) { return 0.8 * n; };
    calx::linalg::Matrix k(N + 1);
    for (int n = 0; n <= N; n++) {
      if (n < N) k(n, n + 1) = up(n);
      if (n > 0) k(n, n - 1) = down(n);
    }
    const auto pi = sm::stationary(k);
    const auto occ = [&](const int& n) {
      std::vector<sm::Occurrence<int>> out;
      if (n < N) out.push_back({"birth", "", up(n), [](const int& m, calx::random::Stream&) { return m + 1; }});
      if (n > 0) out.push_back({"death", "", down(n), [](const int& m, calx::random::Stream&) { return m - 1; }});
      return out;
    };
    const double T = 20000;
    const auto R = sm::simulate<int>(0, occ, T, rng, 10000000);
    std::vector<double> time(N + 1, 0.0);
    double t = 0;
    int z = R.initial;
    bool neighbours = true;
    for (const auto& e : R.events) {
      time[z] += e.time - t;
      neighbours = neighbours && std::abs(e.successor - z) == 1;
      t = e.time, z = e.successor;
    }
    time[z] += T - t;
    double worst = 0;
    for (int n = 0; n <= N; n++) worst = std::max(worst, std::fabs(time[n] / T - pi[n]));
    s.check("a birth–death chain: the time in each state is its stationary probability", worst, 0.01);
    s.require("every successor a neighbour", neighbours && !R.stopped);
    bool replay = true;
    int w = R.initial;
    for (const auto& e : R.events) {
      replay = replay && R.at(std::nextafter(e.time, 0.0)) == w && R.at(e.time) == e.successor;
      w = e.successor;
    }
    s.require("the record replays the path from the initial state alone", replay);
  }
  {
    const std::size_t k = 7;
    const double r = 0.6, T = 5000;
    const auto occ = [&](const std::vector<int>& z) {
      std::vector<sm::Occurrence<std::vector<int>>> out;
      for (std::size_t i = 0; i < z.size(); i++)
        out.push_back({"flip", "switch " + std::to_string(i), r, [i](const std::vector<int>& x, calx::random::Stream&) { auto y = x; y[i] ^= 1; return y; }});
      return out;
    };
    const auto R = sm::simulate<std::vector<int>>(std::vector<int>(k, 0), occ, T, rng, 10000000);
    const double n = static_cast<double>(R.events.size()), want = k * r * T;
    s.check("one mechanism with 7 matches fires 7 r T times on average (relative)", (n - want) / want, 4 * std::sqrt(want) / want);
  }
  {
    const auto occ = [](const long long& n) {
      return std::vector<sm::Occurrence<long long>>{{"birth", "", static_cast<double>(n) * static_cast<double>(n), [](const long long& m, calx::random::Stream&) { return m + 1; }}};
    };
    double mean = 0;
    bool stopped = true;
    for (int trial = 0; trial < 200; trial++) {
      const auto R = sm::simulate<long long>(1, occ, 40.0, rng, 100000);
      stopped = stopped && R.stopped && R.end < 40;
      mean += R.end / 200;
    }
    s.require("births at rate n² accumulate: every path stops at the event limit before t = 40", stopped);
    s.check("the time reached averages π²/6", mean - std::numbers::pi * std::numbers::pi / 6, 0.25);
  }
  return s.done();
}
