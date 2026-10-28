// Checks of physis/chemistry/kinetics.hpp on abstract networks (species X, Y, …; elements e₁, e₂):
// stoichiometric matrices and their conservation laws, a basis of the left kernel over ℤ; element
// and charge balance, AN = 0 and zᵀN = 0, and its failure on an unbalanced reaction; mass action
// against closed forms (first-order decay, the approach to the equilibrium of X ⇌ Y, and second-
// order dimerisation 2X → Y, 1/x = 1/x₀ + 2kt); the stochastic constants c = k α!/Ω^{|α|−1}; the
// generator's drift Σ_y Q(x, y)(y − x) equal to N a(x) state by state; the exact absorption law of
// the neutral network X + Y → 2X, X + Y → 2Y, in which X fixes with probability i/N from i of N
// (the Moran model's martingale), with a symmetric law from a symmetric start; the generator
// commuting with the exchange X ↔ Y, a symmetry of the network; and Gillespie's paths of the
// birth–death network ∅ → X, X → ∅ against its Poisson stationary law.

#include <physis/chemistry/kinetics.hpp>
#include <calx/integrators/runge_kutta.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

namespace pc = physis::chemistry;
using morphe::arithmetic::Integer;
using morphe::arithmetic::Rational;

int main() {
  calx::verification::Suite s;

  // ---- stoichiometry and conservation -----------------------------------------------------------------------------------------
  {
    pc::Network net;
    net.species = {"X", "Y", "Z"};
    net.add("bind", {{"X", 1}, {"Y", 1}}, {{"Z", 1}}, 2);
    net.add("unbind", {{"Z", 1}}, {{"X", 1}, {"Y", 1}}, 1);
    const auto laws = pc::conservation_laws(net);
    bool ok = laws.size() == 2;
    const auto N = pc::stoichiometric_matrix(net);
    for (const auto& y : laws)
      for (std::size_t r = 0; r < 2 && ok; r++) {
        Integer sum = 0;
        for (std::size_t i = 0; i < 3; i++) sum += y[i] * N(i, r);
        ok = ok && sum.is_zero();
      }
    s.require("X + Y ⇌ Z: two conservation laws, in the left kernel of N", ok);

    // A network whose left kernel is spanned by (1, 1, 1, 2): the autocatalytic network of the
    // chemistry checks, with species named abstractly.
    pc::Network f;
    f.species = {"P", "X", "Y", "W"};
    f.add("x", {{"P", 1}, {"X", 1}}, {{"X", 2}}, 1);
    f.add("y", {{"P", 1}, {"Y", 1}}, {{"Y", 2}}, 1);
    f.add("w", {{"X", 1}, {"Y", 1}}, {{"W", 1}}, 1);
    const auto fl = pc::conservation_laws(f);
    const bool primitive = fl.size() == 1 && (fl[0] == std::vector<Integer>{1, 1, 1, 2} || fl[0] == std::vector<Integer>{-1, -1, -1, -2});
    s.require("P + X → 2X, P + Y → 2Y, X + Y → W: the one law ±(1, 1, 1, 2)", primitive);

    // Balance: species built from elements e₁, e₂ with charges.
    morphe::algebra::Matrix A(2, 4);  // e₁, e₂ in P, X, Y, W
    for (std::size_t j = 0; j < 3; j++) A(0, j) = 2, A(1, j) = 1;
    A(0, 3) = 4, A(1, 3) = 2;
    s.require("element and charge balance: AN = 0, zᵀN = 0", pc::balance(f, A, {0, 0, 0, 0}).balanced());
    A(1, 3) = 1;
    s.require("an unbalanced composition is caught", !pc::balance(f, A, {0, 0, 0, 0}).balanced());
    A(1, 3) = 2;
    s.require("an unbalanced charge is caught", !pc::balance(f, A, {0, 1, 0, 0}).balanced());
  }

  // ---- mass action -----------------------------------------------------------------------------------------------------------
  {
    pc::Network net;
    net.species = {"X", "Y"};
    net.add("forward", {{"X", 1}}, {{"Y", 1}}, 3);
    net.add("back", {{"Y", 1}}, {{"X", 1}}, 1);
    calx::integrators::Vec<2> c{1, 0};
    const auto f = [&](const calx::integrators::Vec<2>& x, calx::integrators::Vec<2>& dx) { pc::mass_action(net, x, dx); };
    for (int i = 0; i < 2000; i++) calx::integrators::rk4<2>(f, c, 1e-3);
    // x(t) = 1/4 + 3/4 e^{−4t}
    s.check("X ⇌ Y (k₊ = 3, k₋ = 1) at t = 2 against 1/4 + 3/4 e^{−4t}", c[0] - (0.25 + 0.75 * std::exp(-8.0)), 1e-12);

    pc::Network dim;
    dim.species = {"X", "Y"};
    dim.add("dimerise", {{"X", 2}}, {{"Y", 1}}, Rational(1, 2));
    calx::integrators::Vec<2> d{2, 0};
    const auto g = [&](const calx::integrators::Vec<2>& x, calx::integrators::Vec<2>& dx) { pc::mass_action(dim, x, dx); };
    for (int i = 0; i < 1000; i++) calx::integrators::rk4<2>(g, d, 1e-3);
    // dx/dt = −2k x², 1/x = 1/x₀ + 2kt: x(1) = 1/(1/2 + 1) = 2/3
    s.check("2X → Y (k = 1/2) at t = 1 against 1/x = 1/x₀ + 2kt", d[0] - 2.0 / 3, 1e-12);
    s.require("stochastic constants: c = k·2!/Ω for 2X → Y", pc::stochastic_constant(dim.reactions[0], 10) == Rational(1, 10));
  }

  // ---- the Markov chain ------------------------------------------------------------------------------------------------------
  {
    pc::Network moran;
    moran.species = {"X", "Y"};
    moran.add("X wins", {{"X", 1}, {"Y", 1}}, {{"X", 2}}, 1);
    moran.add("Y wins", {{"X", 1}, {"Y", 1}}, {{"Y", 2}}, 1);
    const long long N = 12;
    bool fix = true;
    for (long long i = 1; i < N; i++) {
      const auto law = pc::absorption(moran, {i, N - i});
      const auto it = law.find({N, 0});
      fix = fix && it != law.end() && it->second == Rational(Integer(i), Integer(N));
    }
    s.require("X + Y → 2X, X + Y → 2Y: X fixes with probability i/N exactly (N = 12)", fix);
    s.require("the exchange X ↔ Y is a symmetry of the network", pc::symmetry_of(moran, {1, 0}));

    // The generator commutes with the exchange on the states reached from (6, 6).
    const auto S = pc::reachable(moran, {6, 6});
    const auto Q = pc::generator(moran, S);
    bool commute = true;
    for (std::size_t x = 0; x < S.states.size(); x++)
      for (const auto& [y, a] : Q[x]) {
        const std::size_t sx = S.find(pc::permute({1, 0}, S.states[x])), sy = S.find(pc::permute({1, 0}, S.states[y]));
        const auto it = Q[sx].find(sy);
        commute = commute && it != Q[sx].end() && it->second == a;
      }
    s.require("its generator commutes with the exchange of states", commute);

    // Dynkin's drift, state by state, on a network of second order.
    pc::Network net;
    net.species = {"X", "Y", "Z"};
    net.add("a", {{"X", 2}}, {{"Y", 1}}, 3);
    net.add("b", {{"X", 1}, {"Y", 1}}, {{"Z", 1}}, 2);
    net.add("c", {{"Z", 1}}, {{"X", 3}}, 1);
    const auto R = pc::reachable(net, {5, 1, 0}, 100000);
    const auto G = pc::generator(net, R, 4);
    bool drift = true;
    for (std::size_t x = 0; x < R.states.size(); x++) {
      std::vector<Rational> d(3, Rational(0));
      for (const auto& [y, a] : G[x])
        for (std::size_t i = 0; i < 3; i++) d[i] = d[i] + a * Rational(R.states[y][i] - R.states[x][i]);
      drift = drift && d == pc::drift(net, R.states[x], 4);
    }
    s.require("Σ_y Q(x, y)(y − x) = N a(x) at every reachable state (" + std::to_string(R.states.size()) + ")", drift);

    // Gillespie: ∅ → X at λ, X → ∅ at μ per molecule; the stationary law is Poisson(λ/μ).
    pc::Network bd;
    bd.species = {"X"};
    bd.add("birth", {}, {{"X", 1}}, 20);
    bd.add("death", {{"X", 1}}, {}, 2);
    calx::random::Stream rng(2026, 10, 28, 0);
    std::vector<double> t;
    std::vector<pc::State> path;
    pc::gillespie(bd, {0}, 5000, rng, t, path);
    double area = 0, square = 0;
    for (std::size_t k = 0; k + 1 < t.size(); k++) {
      if (t[k] < 50) continue;  // after the transient
      const double dt = t[k + 1] - std::max(t[k], 50.0);
      area += dt * path[k][0];
      square += dt * path[k][0] * path[k][0];
    }
    const double T = t.back() - 50, mean = area / T, var = square / T - mean * mean;
    s.check("birth–death: time-averaged mean against λ/μ = 10", mean - 10, 0.25);
    s.check("and variance against λ/μ = 10 (Poisson)", var - 10, 1.0);
  }
  return s.done();
}
