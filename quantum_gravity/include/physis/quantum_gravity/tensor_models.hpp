// Coloured tensor models and group field theories: their Feynman expansion.
//
// A coloured group field theory in dimension D (Gurau 2011, Ann. Henri Poincaré 12, 829;
// Oriti 2009, in Approaches to Quantum Gravity, ch. 16) has D + 1 complex fields φ_c on G^D,
// a Gaussian measure (the propagator, with the gauge projection for a group G) and the two
// interactions λ Π_c φ_c and λ̄ Π_c φ̄_c, each a D-simplex whose facets are the fields. Its
// Feynman graphs are the (D+1)-coloured graphs of morphe's ColouredGraph: a black vertex for
// each λ, a white one for each λ̄, an edge of colour c for each propagator of φ_c. With G
// trivial and indices in {1..N} it is the coloured tensor model, invariant under U(N)^{D+1}.
//
// Expansion. At order (λλ̄)^p the Wick contractions of colour c are the bijections from the
// black to the white vertices, so the labelled graphs are the tuples (σ_0, ..., σ_D) of
// permutations, (p!)^{D+1} of them, each weighted 1/(p!)² by the exponential series; fixing
// σ_0 = id by relabelling the whites, the free energy log Z has the coefficient
//   F_p = (1/p!) Σ_{connected tuples with σ_0 = id} A(Γ),
// which is Σ_Γ A(Γ)/|Aut Γ| over unlabelled connected graphs, and (rooting at one of the p
// black vertices) 1/p times the sum over rooted graphs.
//
// Amplitudes. In the tensor model every face (bicoloured cycle) is a free index sum, so
// A(Γ) = N^{F(Γ)}, and with the coupling scaled as λ N^{−D(D−1)/4} the order (λλ̄)^p carries
// N^{F − D(D−1)p/2} = N^{D − 2ω/(D−1)!}, ω Gurau's degree (morphe topology): the 1/N
// expansion, led by the melonic graphs (ω = 0), counted with a root by the Fuss–Catalan
// numbers (1/((D+1)p + 1)) C((D+1)p + 1, p) (Bonzom, Gurau, Riello and Rivasseau 2011, Nucl.
// Phys. B 853, 174). With a finite group G (the Boulatov model of three-dimensional gravity
// when D = 3, Boulatov 1992; Ooguri 1992 for D = 4, both for compact Lie groups) the
// amplitude is the BF state sum on the graph's 2-complex, bf.hpp, = |G|^{v − 1} |Hom(π1, G)|
// for D = 3, v the number of vertices of the dual triangulation (3-bubbles).
//
// The enumeration visits (p!)^D tuples: p ≤ 5 for D = 3.

#pragma once

#include <cstddef>
#include <map>
#include <numeric>
#include <vector>

#include <morphe/algebra/finite_group.hpp>
#include <morphe/arithmetic/rational.hpp>
#include <morphe/topology/coloured_graph.hpp>

#include "bf.hpp"

namespace physis::quantum_gravity {

using morphe::topology::ColouredGraph;

// Calls f(Γ) for every (D+1)-coloured graph with p black vertices and σ_0 the identity.
template <class F>
void for_each_coloured_graph(std::size_t D, std::size_t p, F&& f) {
  std::vector<std::size_t> id(p);
  std::iota(id.begin(), id.end(), 0);
  ColouredGraph G{D + 1, p, std::vector<std::vector<std::size_t>>(D + 1, id)};
  // An odometer over σ_1, ..., σ_D, each running through all permutations.
  for (;;) {
    f(static_cast<const ColouredGraph&>(G));
    std::size_t c = 1;
    while (c <= D && !std::next_permutation(G.sigma[c].begin(), G.sigma[c].end())) c++;  // wraps to the identity
    if (c > D) return;
  }
}

// The power of N carried by a graph in the tensor model with the scaled coupling:
// F − D(D−1)p/2 = D − 2ω/(D−1)!.
inline long long amplitude_exponent(const ColouredGraph& G) {
  const long long D = static_cast<long long>(G.dimension());
  return static_cast<long long>(G.face_count()) - D * (D - 1) * static_cast<long long>(G.p) / 2;
}

// The coefficients of the free energy at order (λλ̄)^p, by degree ω: F_p = Σ_ω c_ω N^{D − 2ω/(D−1)!}.
inline std::map<std::size_t, Rational> free_energy(std::size_t D, std::size_t p) {
  std::map<std::size_t, Integer> count;
  for_each_coloured_graph(D, p, [&](const ColouredGraph& G) {
    if (G.connected()) count[G.degree()] += 1;
  });
  std::map<std::size_t, Rational> out;
  const Integer pf = morphe::arithmetic::factorial(p);
  for (const auto& [omega, n] : count) out[omega] = Rational(n, pf);
  return out;
}

// The Fuss–Catalan number (1/(m p + 1)) C(m p + 1, p): rooted melonic (m)-coloured graphs
// with 2p vertices, m = D + 1.
inline Integer fuss_catalan(long long m, long long p) {
  return morphe::arithmetic::exact_div(morphe::arithmetic::binomial(m * p + 1, p), Integer(m * p + 1));
}

// The amplitude of a graph in the coloured group field theory with a finite group: BF theory
// on its 2-complex (normalised Haar measure, δ propagator projections).
inline Rational boulatov_amplitude(const ColouredGraph& G, const FiniteGroup& group) { return bf_partition_function(G.two_complex(), group); }

}  // namespace physis::quantum_gravity
