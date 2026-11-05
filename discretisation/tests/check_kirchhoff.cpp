// Checks of networks of conductances (kirchhoff.hpp; decision 0022).
//
// - Series and parallel: conductances g₁, g₂ give g₁g₂/(g₁ + g₂) and g₁ + g₂.
// - A balanced Wheatstone bridge carries nothing across its bridge.
// - Kirchhoff's current law at every free vertex, with injections, to rounding.
// - Thomson's principle: moving a free potential raises the dissipation.
// - The conjugate gradients' branch (a 24 × 24 grid, 574 free vertices) agrees with a dense solve of the same
//   reduced system (1e−10 relative), and a 2 × 2 grid's corner-to-corner conductance is 1 (by symmetry).
// - A part of the network touching no fixed potential is refused.

#include <physis/discretisation/kirchhoff.hpp>

#include <calx/linalg/dense.hpp>
#include <calx/verification/check.hpp>

#include <cmath>

using namespace physis::discretisation;

int main() {
  calx::verification::Suite s;
  {
    ConductanceNetwork n;
    n.vertices = 3;
    n.add(0, 1, 2.0), n.add(1, 2, 3.0);
    s.check("series: g₁g₂/(g₁ + g₂)", effective_conductance(n, 0, 2) - 1.2, 1e-14);
    ConductanceNetwork p;
    p.vertices = 2;
    p.add(0, 1, 2.0), p.add(0, 1, 3.0);
    s.check("parallel: g₁ + g₂", effective_conductance(p, 0, 1) - 5.0, 1e-14);
  }
  {
    // A bridge: 0 → 1 → 3 and 0 → 2 → 3 with ratios equal; bridge 1–2.
    ConductanceNetwork n;
    n.vertices = 4;
    n.add(0, 1, 1.0), n.add(1, 3, 2.0), n.add(0, 2, 3.0), n.add(2, 3, 6.0);
    const std::size_t bridge = n.add(1, 2, 5.0);
    const auto sol = solve_kirchhoff(n, {{0, 1.0}, {3, 0.0}});
    s.check("Wheatstone: the balanced bridge carries nothing", sol.flow[bridge], 1e-14);
  }
  {
    ConductanceNetwork n;
    n.vertices = 6;
    n.add(0, 1, 1.0), n.add(1, 2, 0.5), n.add(2, 3, 2.0), n.add(1, 4, 1.5), n.add(4, 3, 0.7), n.add(4, 5, 1.1), n.add(5, 2, 0.3);
    const std::vector<double> inj{0, 0.2, -0.1, 0, 0.05, 0.3};
    const auto sol = solve_kirchhoff(n, {{0, 1.0}, {3, -0.5}}, inj);
    const auto div = divergence(n, sol.flow);
    double worst = 0;
    for (std::size_t v : {1u, 2u, 4u, 5u}) worst = std::max(worst, std::fabs(div[v] - inj[v]));
    s.check("Kirchhoff's current law at the free vertices", worst, 1e-14);
    auto dissipation = [&](std::vector<double> phi) {
      double d = 0;
      for (const auto& e : n.edges) d += e.conductance * std::pow(phi[e.a] - phi[e.b], 2);
      return d;
    };
    const auto sol0 = solve_kirchhoff(n, {{0, 1.0}, {3, -0.5}});
    auto moved = sol0.potential;
    moved[4] += 0.01;
    s.require("Thomson's principle: moving a free potential raises the dissipation", dissipation(moved) > sol0.dissipation);
  }
  {
    const std::size_t w = 24;
    ConductanceNetwork n;
    n.vertices = w * w;
    for (std::size_t y = 0; y < w; y++)
      for (std::size_t x = 0; x < w; x++) {
        const double g = 1 + 0.5 * std::sin(0.3 * x + 0.7 * y);
        if (x + 1 < w) n.add(y * w + x, y * w + x + 1, g);
        if (y + 1 < w) n.add(y * w + x, (y + 1) * w + x, g * 0.8);
      }
    const auto sol = solve_kirchhoff(n, {{0, 1.0}, {w * w - 1, 0.0}});
    // The same reduced system, densely.
    const std::size_t m = w * w - 2;
    calx::linalg::Matrix A(m), b(m);
    auto idx = [&](std::size_t v) { return v - 1; };
    std::vector<double> rhs(m, 0.0);
    for (const auto& e : n.edges) {
      const bool fa = e.a == 0 || e.a == w * w - 1, fb = e.b == 0 || e.b == w * w - 1;
      if (!fa) A(idx(e.a), idx(e.a)) += e.conductance;
      if (!fb) A(idx(e.b), idx(e.b)) += e.conductance;
      if (!fa && !fb) A(idx(e.a), idx(e.b)) -= e.conductance, A(idx(e.b), idx(e.a)) -= e.conductance;
      else if (!fa) rhs[idx(e.a)] += e.conductance * (e.b == 0 ? 1.0 : 0.0);
      else if (!fb) rhs[idx(e.b)] += e.conductance * (e.a == 0 ? 1.0 : 0.0);
    }
    for (std::size_t i = 0; i < m; i++) b(i, 0) = rhs[i];
    const auto X = calx::linalg::solve(A, b);
    double err = 0;
    for (std::size_t v = 1; v + 1 < w * w; v++) err = std::max(err, std::fabs(sol.potential[v] - X(idx(v), 0)));
    s.check("conjugate gradients agree with a dense solve (574 free vertices)", err, 1e-10);
    ConductanceNetwork sq;
    sq.vertices = 4;
    sq.add(0, 1, 1), sq.add(1, 3, 1), sq.add(0, 2, 1), sq.add(2, 3, 1);
    s.check("a square's corner-to-corner conductance is 1", effective_conductance(sq, 0, 3) - 1, 1e-14);
  }
  {
    ConductanceNetwork n;
    n.vertices = 4;
    n.add(0, 1, 1.0), n.add(2, 3, 1.0);
    bool refused = false;
    try {
      solve_kirchhoff(n, {{0, 1.0}, {1, 0.0}});
    } catch (const std::invalid_argument&) {
      refused = true;
    }
    s.require("a part touching no fixed potential is refused", refused);
  }
  return s.done();
}
