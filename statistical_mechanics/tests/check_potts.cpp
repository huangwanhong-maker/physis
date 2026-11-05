// Checks of the Potts model (potts.hpp; decision 0022).
//
// - Energy: all sites in one state give −J × (pairs); on a 4 × 4 lattice in two vertical stripes of two
//   columns each, −J × 24 (von Neumann: 32 pairs, 8 across the stripes' two boundaries).
// - Exact enumeration: the mean energy of a 3 × 3 lattice, q = 3, von Neumann, at T = 1.5, and q = 2, Moore,
//   at T = 2.5, by enumeration; Metropolis with uniform proposals agrees (within 4 standard errors by batch
//   means).
// - Coarsening: with neighbour proposals at T = 0 the energy never increases.
// - Domains: the stripes are two cells of area 8, each touching one other; a checkerboard of 3 × 3 blocks on a
//   6 × 6 torus is four cells, each with two neighbours.

#include <physis/statistical_mechanics/potts.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <vector>

using namespace physis::statistical_mechanics;

int main() {
  calx::verification::Suite s;
  {
    PottsLattice L(4, 4, 3, 1.0);
    s.check("energy: one state everywhere is −J × 32 pairs", L.energy() + 32, 1e-12);
    for (std::size_t y = 0; y < 4; y++)
      for (std::size_t x = 0; x < 4; x++) L.set_state(L.site(x, y), x < 2 ? 0 : 1);
    s.check("energy: two stripes, −J × 24 (8 of the 32 pairs cross a boundary)", L.energy() + 24, 1e-12);
    const auto d = L.domains();
    s.require("domains: the stripes are two cells of area 8", d.area.size() == 2 && d.area[0] == 8 && d.area[1] == 8);
    s.require("domains: each touching the other alone", d.sides[0] == 1 && d.sides[1] == 1);
  }
  {
    PottsLattice L(6, 6, 4, 1.0);
    for (std::size_t y = 0; y < 6; y++)
      for (std::size_t x = 0; x < 6; x++) L.set_state(L.site(x, y), (x / 3) + 2 * (y / 3));
    const auto d = L.domains();
    bool ok = d.area.size() == 4;
    for (std::size_t c = 0; ok && c < 4; c++) ok = d.area[c] == 9 && d.sides[c] == 2;
    s.require("domains: four blocks on a torus, each with two neighbours", ok);
  }
  auto metropolis = [&](std::size_t q, Neighbourhood nb, double T, std::uint64_t seed) {
    PottsLattice L(3, 3, q, 1.0, nb);
    calx::random::Stream rng(seed, 1, 2, 3);
    for (int k = 0; k < 2000; k++) L.sweep(T, PottsProposal::uniform, rng);
    const int batches = 50, per = 4000;
    std::vector<double> means;
    for (int b = 0; b < batches; b++) {
      double e = 0;
      for (int k = 0; k < per; k++) L.sweep(T, PottsProposal::uniform, rng), e += L.energy();
      means.push_back(e / per);
    }
    double m = 0, v = 0;
    for (double x : means) m += x;
    m /= batches;
    for (double x : means) v += (x - m) * (x - m);
    return std::pair<double, double>{m, std::sqrt(v / (batches - 1) / batches)};
  };
  {
    const double exact = potts_exact_mean_energy(3, 3, 3, 1.0, Neighbourhood::von_neumann, 1.5);
    const auto [m, se] = metropolis(3, Neighbourhood::von_neumann, 1.5, 17);
    s.check("Metropolis against enumeration, q = 3, von Neumann, T = 1.5 (in standard errors)", (m - exact) / se, 4);
    const double exact2 = potts_exact_mean_energy(3, 3, 2, 1.0, Neighbourhood::moore, 2.5);
    const auto [m2, se2] = metropolis(2, Neighbourhood::moore, 2.5, 19);
    s.check("Metropolis against enumeration, q = 2, Moore, T = 2.5 (in standard errors)", (m2 - exact2) / se2, 4);
  }
  {
    PottsLattice L(40, 40, 1600, 1.0, Neighbourhood::moore);
    for (std::size_t i = 0; i < 1600; i++) L.set_state(i, i);
    calx::random::Stream rng(5, 0, 0, 0);
    double e = L.energy();
    bool never_up = true;
    for (int k = 0; k < 60; k++) {
      L.sweep(0, PottsProposal::neighbour, rng);
      const double f = L.energy();
      if (f > e + 1e-12) never_up = false;
      e = f;
    }
    s.require("coarsening at T = 0: the energy never increases", never_up);
    s.require("coarsening: cells merge (fewer than a third remain after 60 sweeps)", L.domains().area.size() < 1600 / 3);
  }
  return s.done();
}
