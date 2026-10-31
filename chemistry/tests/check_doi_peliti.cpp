// Checks of the Doi–Peliti representation (doi_peliti.hpp).
//
// - The ladder operators on the box {0…6}²: [a_i, a_i†] = 1 below the top occupation (broken only there),
//   [a_i, a_j†] = 0 for i ≠ j, a†a the number of each mode, ⟨1| a† = ⟨1| below the top.
// - The network 2A → B, B → 2A, A + B → 3A (A + 2B conserved, second order with identical reactants):
//   the operator Σ_r L_r equals the transpose of kinetics' generator on the reachable states from (6, 0),
//   entry by entry, in volumes 1 and 3 — Gillespie's c_r h_r(n) and the Doi–Peliti κ n^{(α)} are one
//   chain — every column of the operator sums to zero (⟨1| L = 0), and ⟨1| a_A† a_A L |n⟩ is Dynkin's drift
//   Σ_r a_r(n) N_Ar.

#include <cmath>
#include <vector>

#include <calx/verification/check.hpp>
#include <physis/chemistry/doi_peliti.hpp>

using calx::linalg::Matrix;
using calx::verification::Suite;
namespace pc = physis::chemistry;

int main() {
  Suite s;
  const pc::Fock F{2, 6};
  const std::size_t n = F.size();
  {
    double comm = 0, cross = 0, number = 0, projection = 0;
    for (std::size_t i = 0; i < 2; i++) {
      const auto a = pc::annihilation(F, i), ad = pc::creation(F, i);
      const auto C = a * ad - ad * a;
      for (std::size_t k = 0; k < n; k++) {
        const auto occ = F.occupation(k);
        if (occ[i] < 6)
          for (std::size_t m = 0; m < n; m++) comm = std::max(comm, std::fabs(C(m, k) - (m == k ? 1.0 : 0.0)));
        const auto N = ad * a;
        number = std::max(number, std::fabs(N(k, k) - static_cast<double>(occ[i])));
        // ⟨1| a† |n⟩ = 1 unless n_i is at the top.
        double col = 0;
        for (std::size_t m = 0; m < n; m++) col += ad(m, k);
        if (occ[i] < 6) projection = std::max(projection, std::fabs(col - 1));
      }
      const auto b = pc::annihilation(F, 1 - i), bd = pc::creation(F, 1 - i);
      const auto D = a * bd - bd * a;
      for (double v : D.a) cross = std::max(cross, std::fabs(v));
    }
    s.check("[a_i, a_i†] = 1 below the top occupation", comm, 1e-15);
    s.check("[a_i, a_j†] = 0 for i ≠ j", cross, 1e-15);
    s.check("a†a counts each mode", number, 1e-15);
    s.check("⟨1| a† = ⟨1| below the top", projection, 1e-15);
  }
  pc::Network net;
  net.species = {"A", "B"};
  net.add("dimerise", {{"A", 2}}, {{"B", 1}}, pc::Rational(3, 2));
  net.add("split", {{"B", 1}}, {{"A", 2}}, pc::Rational(7, 10));
  net.add("convert", {{"A", 1}, {"B", 1}}, {{"A", 3}}, pc::Rational(2, 5));
  const pc::State n0{6, 0};
  const auto S = pc::reachable(net, n0);
  for (long long volume : {1LL, 3LL}) {
    const auto L = pc::network_operator(F, net, volume);
    const auto Q = pc::generator<double>(net, S, volume);
    double worst = 0, columns = 0;
    for (std::size_t x = 0; x < S.states.size(); x++) {
      const std::size_t kx = F.index({S.states[x][0], S.states[x][1]});
      double out = 0;
      for (std::size_t y = 0; y < S.states.size(); y++) {
        const std::size_t ky = F.index({S.states[y][0], S.states[y][1]});
        const auto it = Q[x].find(y);
        const double q = it == Q[x].end() ? 0.0 : it->second;
        if (y != x) worst = std::max(worst, std::fabs(L(ky, kx) - q)), out += q;
      }
      worst = std::max(worst, std::fabs(L(kx, kx) + out));
      double col = 0;
      for (std::size_t m = 0; m < n; m++) col += L(m, kx);
      columns = std::max(columns, std::fabs(col));
    }
    s.require("the reachable states from (6, 0): (6,0) (4,1) (2,2) (0,3)", S.states.size() == 4);
    s.check(volume == 1 ? "Σ L_r is the transpose of the master equation's generator (Ω = 1)" : "and in a volume Ω = 3", worst, 1e-12);
    s.check(volume == 1 ? "⟨1| L = 0: probability is kept (Ω = 1)" : "⟨1| L = 0 (Ω = 3)", columns, 1e-12);
  }
  {
    const auto L = pc::network_operator(F, net, 1);
    const auto NA = pc::creation(F, 0) * pc::annihilation(F, 0);
    const auto NL = NA * L;
    double worst = 0;
    for (const auto& st : S.states) {
      const std::size_t k = F.index({st[0], st[1]});
      double lhs = 0;
      for (std::size_t m = 0; m < n; m++) lhs += NL(m, k);
      double drift = 0;
      for (std::size_t r = 0; r < net.reactions.size(); r++) drift += pc::propensity<double>(net, r, st, 1) * static_cast<double>(net.change(0, r));
      worst = std::max(worst, std::fabs(lhs - drift));
    }
    s.check("⟨1| a_A† a_A L |n⟩ is Dynkin's drift of n_A", worst, 1e-12);
  }
  return s.done();
}
