// The Doi–Peliti representation of a reaction network: creation and annihilation operators on the
// occupation numbers, the reaction operators they compose, and the identity of the result with the master
// equation of the network's Markov chain (kinetics.hpp) — the field-theoretic notation of counted entities,
// with no quantum mechanics in it (decision 0018).
//
// The algebra (Doi 1976, J. Phys. A 9, 1465 and 1479; Peliti 1985, J. Physique 46, 1469; Täuber 2014,
// Critical Dynamics, ch. 9). For s species, the occupation states |n⟩, n ∈ ℕˢ, and on them
//
//   a_i† |n⟩ = |n + e_i⟩,   a_i |n⟩ = n_i |n − e_i⟩,   [a_i, a_j†] = δ_ij,
//
// unnormalised (the Doi convention): a_i† a_i |n⟩ = n_i |n⟩. A law P is the vector |P⟩ = Σ P(n)|n⟩; the
// projection ⟨1| = Σ ⟨n| gives total probability, ⟨1| a_i† = ⟨1|, and expectations are ⟨1| g(a†a) |P⟩ —
// a sum of probabilities, not the Born rule. A reaction Σ α_i X_i → Σ β_i X_i of constant κ acts by
//
//   L_r = κ [Π (a_i†)^{β_i} − Π (a_i†)^{α_i}] Π a_i^{α_i},
//
// since Π a^α |n⟩ = Π n_i^{(α_i)} |n − α⟩ with n^{(α)} the falling factorial: L_r |n⟩ = κ n^{(α)} (|n + β − α⟩
// − |n⟩), gain minus loss, and d|P⟩/dt = Σ_r L_r |P⟩ is the master equation. Its propensity κ n^{(α)} counts
// ordered tuples of reactants; Gillespie's c_r h_r(n) counts sets, c_r = k_r Π α_i! / Ω^{|α|−1} (kinetics.hpp),
// so κ = k_r / Ω^{|α|−1}: the two conventions describe the same chain, and network_operator() uses kinetics'
// constants so that the operator equals the transpose of kinetics' generator exactly.
//
// Truncation. The operators are matrices on the box {0, …, N}ˢ, a_i† sending n_i = N to 0. The truncated
// algebra keeps [a_i, a_i†] = 1 below the top and breaks it at n_i = N; a network whose reachable states stay
// below the top (a conservation law bounding them) is represented exactly. Dense: (N + 1)ˢ states, a few
// thousand at most. Checks: tests/check_doi_peliti.cpp.

#pragma once

#include <cstddef>
#include <stdexcept>
#include <vector>

#include <calx/linalg/dense.hpp>

#include "kinetics.hpp"

namespace physis::chemistry {

// The box {0, …, cap}^species, states indexed in mixed radix (species 0 fastest).
struct Fock {
  std::size_t species = 1, cap = 1;
  std::size_t size() const {
    std::size_t m = 1;
    for (std::size_t i = 0; i < species; i++) m *= cap + 1;
    return m;
  }
  std::size_t index(const std::vector<long long>& n) const {
    std::size_t k = 0, w = 1;
    for (std::size_t i = 0; i < species; i++) {
      if (n[i] < 0 || static_cast<std::size_t>(n[i]) > cap) return size();
      k += static_cast<std::size_t>(n[i]) * w;
      w *= cap + 1;
    }
    return k;
  }
  std::vector<long long> occupation(std::size_t k) const {
    std::vector<long long> n(species);
    for (std::size_t i = 0; i < species; i++) n[i] = static_cast<long long>(k % (cap + 1)), k /= cap + 1;
    return n;
  }
};

// Matrices act on column vectors |P⟩: entry (m, n) = ⟨m| O |n⟩.
inline calx::linalg::Matrix creation(const Fock& F, std::size_t i) {
  calx::linalg::Matrix A(F.size());
  for (std::size_t k = 0; k < F.size(); k++) {
    auto n = F.occupation(k);
    n[i] += 1;
    const std::size_t m = F.index(n);
    if (m < F.size()) A(m, k) = 1;
  }
  return A;
}
inline calx::linalg::Matrix annihilation(const Fock& F, std::size_t i) {
  calx::linalg::Matrix A(F.size());
  for (std::size_t k = 0; k < F.size(); k++) {
    auto n = F.occupation(k);
    if (n[i] == 0) continue;
    const double c = static_cast<double>(n[i]);
    n[i] -= 1;
    A(F.index(n), k) = c;
  }
  return A;
}
inline calx::linalg::Matrix power(const calx::linalg::Matrix& A, unsigned k) {
  calx::linalg::Matrix P = calx::linalg::Matrix::identity(A.n);
  for (unsigned q = 0; q < k; q++) P = P * A;
  return P;
}

// L_r = κ [Π (a†)^β − Π (a†)^α] Π a^α, composed from the ladder operators.
inline calx::linalg::Matrix reaction_operator(const Fock& F, const Network& net, std::size_t r, double kappa) {
  calx::linalg::Matrix down = calx::linalg::Matrix::identity(F.size()), gain = down, loss = down;
  for (std::size_t i = 0; i < F.species; i++) {
    const auto a = static_cast<unsigned>(net.alpha(i, r)), b = static_cast<unsigned>(net.beta(i, r));
    if (a) down = power(annihilation(F, i), a) * down;
    if (b) gain = power(creation(F, i), b) * gain;
    if (a) loss = power(creation(F, i), a) * loss;
  }
  return kappa * (gain * down - loss * down);
}

// Σ_r L_r with κ_r chosen so that κ_r n^{(α)} = c_r h_r(n), Gillespie's propensity in a volume Ω.
inline calx::linalg::Matrix network_operator(const Fock& F, const Network& net, long long volume = 1) {
  if (net.species.size() != F.species) throw std::invalid_argument("physis::chemistry::network_operator: one mode for each species");
  calx::linalg::Matrix L(F.size());
  for (std::size_t r = 0; r < net.reactions.size(); r++) {
    // c_r h_r(n) = c_r Π n^{(α)}/α! , so κ = c_r / Π α!.
    double factorials = 1;
    for (std::size_t i = 0; i < F.species; i++)
      for (long long q = 2; q <= net.alpha(i, r); q++) factorials *= static_cast<double>(q);
    const double kappa = stochastic_constant(net.reactions[r], volume).to_double() / factorials;
    L = L + reaction_operator(F, net, r, kappa);
  }
  return L;
}

}  // namespace physis::chemistry
