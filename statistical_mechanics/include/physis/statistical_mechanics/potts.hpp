// The q-state Potts model on a periodic square lattice: Metropolis dynamics, exact enumeration of small
// lattices, and the domains of a configuration (decision 0022).
//
// Each site carries a state s ∈ {0, …, q − 1} and the energy is H = −J Σ_{⟨ij⟩} δ(s_i, s_j) over the pairs
// of a neighbourhood: the four nearest (von Neumann) or the eight nearest and next nearest (Moore)
// (Potts 1952; Wu 1982, Rev. Mod. Phys. 54, 235). With q = 2 it is the Ising model at coupling J/2.
//
// Dynamics. A sweep visits L² sites drawn at random; at each a new state is proposed and accepted with
// probability min(1, e^{−ΔH/T}) (Metropolis et al. 1953), or when ΔH ≤ 0 at T = 0. Two proposals:
// - uniform, a state drawn uniformly from the other q − 1: symmetric, so the chain is reversible with
//   respect to the Boltzmann law e^{−H/T}/Z, the check of the dynamics against exact enumeration;
// - neighbour, the state of a randomly chosen neighbour: the proposal of the Potts model of grain growth
//   and of foams (Anderson, Srolovitz, Grest and Sahni 1984; Glazier, Anderson and Grest 1990), with
//   large q, in which each state is a cell, and T small. It is not symmetric, and models coarsening, not
//   equilibrium: at T = 0 the energy never increases.
// The Moore neighbourhood lessens the lattice's anisotropy of boundary energy, which pins boundaries
// along its axes at T = 0 (Glazier et al. 1990).
//
// Domains. The cells of a configuration are the connected sets of equal state (four-connected); each has
// its area (sites) and its sides, the number of distinct cells it touches across an edge. On a torus the
// mean number of sides of a cellular pattern with three cells at each vertex is 6 (Euler's formula); the
// lattice admits four at some vertices, so the mean is near 6, not exactly 6.
//
// Draws are from calx's keyed streams: the same seed gives the same run. Checks: tests/check_potts.cpp.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <calx/random/stream.hpp>

namespace physis::statistical_mechanics {

enum class Neighbourhood { von_neumann, moore };
enum class PottsProposal { uniform, neighbour };

class PottsLattice {
 public:
  PottsLattice(std::size_t width, std::size_t height, std::size_t q, double coupling = 1, Neighbourhood nb = Neighbourhood::von_neumann)
      : w_(width), h_(height), q_(q), J_(coupling), nb_(nb), s_(width * height, 0) {
    if (width < 3 || height < 3 || q < 2) throw std::invalid_argument("physis::statistical_mechanics::PottsLattice: at least 3 × 3 sites and 2 states");
  }

  std::size_t width() const { return w_; }
  std::size_t height() const { return h_; }
  std::size_t states() const { return q_; }
  std::size_t site(std::size_t x, std::size_t y) const { return (y % h_) * w_ + (x % w_); }
  std::size_t state(std::size_t i) const { return s_[i]; }
  void set_state(std::size_t i, std::size_t v) {
    if (v >= q_) throw std::out_of_range("physis::statistical_mechanics::PottsLattice: a state below q");
    s_[i] = v;
  }
  const std::vector<std::size_t>& states_of_sites() const { return s_; }

  // The neighbours of site i in the lattice's neighbourhood.
  std::vector<std::size_t> neighbours(std::size_t i) const {
    const std::size_t x = i % w_, y = i / w_;
    std::vector<std::size_t> n{site(x + 1, y), site(x + w_ - 1, y), site(x, y + 1), site(x, y + h_ - 1)};
    if (nb_ == Neighbourhood::moore)
      for (std::size_t dy : {std::size_t(1), h_ - 1})
        for (std::size_t dx : {std::size_t(1), w_ - 1}) n.push_back(site(x + dx, y + dy));
    return n;
  }

  double energy() const {
    double e = 0;
    for (std::size_t i = 0; i < s_.size(); i++)
      for (std::size_t j : neighbours(i))
        if (s_[i] == s_[j]) e -= J_;
    return 0.5 * e;  // each pair twice
  }

  // The change of energy if site i took state v.
  double delta(std::size_t i, std::size_t v) const {
    double d = 0;
    for (std::size_t j : neighbours(i)) d += J_ * ((s_[j] == s_[i] ? 1.0 : 0.0) - (s_[j] == v ? 1.0 : 0.0));
    return d;
  }

  // One sweep of L² Metropolis steps at temperature T; returns the number accepted.
  std::size_t sweep(double temperature, PottsProposal proposal, calx::random::Stream& rng) {
    std::size_t accepted = 0;
    for (std::size_t k = 0; k < s_.size(); k++) {
      const std::size_t i = static_cast<std::size_t>(rng.below(s_.size()));
      std::size_t v;
      if (proposal == PottsProposal::uniform) {
        v = static_cast<std::size_t>(rng.below(q_ - 1));
        if (v >= s_[i]) v++;
      } else {
        const auto n = neighbours(i);
        v = s_[n[static_cast<std::size_t>(rng.below(n.size()))]];
        if (v == s_[i]) continue;
      }
      const double d = delta(i, v);
      const bool accept = d <= 0 || (temperature > 0 && rng.uniform() < std::exp(-d / temperature));
      if (accept) s_[i] = v, accepted++;
    }
    return accepted;
  }

  struct Domains {
    std::vector<std::size_t> label;              // each site's cell
    std::vector<std::size_t> area;               // each cell's sites
    std::vector<std::size_t> sides;              // each cell's distinct neighbouring cells
    std::vector<std::size_t> state;              // each cell's state
  };

  // The four-connected cells of the configuration (periodic).
  Domains domains() const {
    const std::size_t n = s_.size(), none = static_cast<std::size_t>(-1);
    Domains d;
    d.label.assign(n, none);
    std::vector<std::size_t> stack;
    for (std::size_t i = 0; i < n; i++) {
      if (d.label[i] != none) continue;
      const std::size_t c = d.area.size();
      d.area.push_back(0), d.state.push_back(s_[i]);
      d.label[i] = c, stack.push_back(i);
      while (!stack.empty()) {
        const std::size_t v = stack.back();
        stack.pop_back();
        d.area[c]++;
        const std::size_t x = v % w_, y = v / w_;
        for (std::size_t u : {site(x + 1, y), site(x + w_ - 1, y), site(x, y + 1), site(x, y + h_ - 1)})
          if (d.label[u] == none && s_[u] == s_[v]) d.label[u] = c, stack.push_back(u);
      }
    }
    std::vector<std::vector<std::size_t>> touching(d.area.size());
    for (std::size_t i = 0; i < n; i++) {
      const std::size_t x = i % w_, y = i / w_;
      for (std::size_t u : {site(x + 1, y), site(x, y + 1)})
        if (d.label[u] != d.label[i]) touching[d.label[i]].push_back(d.label[u]), touching[d.label[u]].push_back(d.label[i]);
    }
    d.sides.resize(d.area.size());
    for (std::size_t c = 0; c < touching.size(); c++) {
      auto& t = touching[c];
      std::sort(t.begin(), t.end());
      d.sides[c] = static_cast<std::size_t>(std::unique(t.begin(), t.end()) - t.begin());
    }
    return d;
  }

 private:
  std::size_t w_, h_, q_;
  double J_;
  Neighbourhood nb_;
  std::vector<std::size_t> s_;
};

// The exact mean energy ⟨H⟩ of a lattice of at most 20 sites at temperature T > 0, by enumeration of all q^N
// configurations.
inline double potts_exact_mean_energy(std::size_t width, std::size_t height, std::size_t q, double coupling, Neighbourhood nb, double temperature) {
  PottsLattice L(width, height, q, coupling, nb);
  const std::size_t n = width * height;
  double total = 1;
  for (std::size_t i = 0; i < n; i++) total *= static_cast<double>(q);
  if (total > 4e6) throw std::invalid_argument("physis::statistical_mechanics::potts_exact_mean_energy: at most about four million configurations");
  std::vector<std::size_t> digits(n, 0);
  // The energies are small integers times J; subtract the ground state for stability.
  double z = 0, ez = 0;
  const double e0 = -coupling * static_cast<double>(n * (nb == Neighbourhood::moore ? 4 : 2));
  for (std::uint64_t c = 0; c < static_cast<std::uint64_t>(total); c++) {
    for (std::size_t i = 0; i < n; i++) L.set_state(i, digits[i]);
    const double e = L.energy();
    const double w = std::exp(-(e - e0) / temperature);
    z += w, ez += w * e;
    for (std::size_t i = 0; i < n; i++) {
      if (++digits[i] < q) break;
      digits[i] = 0;
    }
  }
  return ez / z;
}

}  // namespace physis::statistical_mechanics
