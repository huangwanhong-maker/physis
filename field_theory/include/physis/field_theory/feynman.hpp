// Vacuum Feynman graphs of φ^k theory, with their symmetry factors.
//
// With the interaction (λ/k!) φ^k, the n-th order of Z/Z₀ is (-λ/k!)^n/n! ⟨(φ^k)^n⟩₀, and by
// Wick's theorem ⟨(φ^k)^n⟩₀ is the number of perfect matchings (Wick pairings) of the kn
// half-edges of n vertices with k legs each, (kn - 1)!!. A pairing is a Feynman graph with
// labelled vertices and legs. The group S_k ≀ S_n, of order (k!)^n n!, permutes the legs at
// each vertex and the vertices; its orbits are the isomorphism classes of graphs and the
// stabiliser of a pairing is the automorphism group of its graph, so a graph Γ arises from
// (k!)^n n! / |Aut Γ| pairings and
//
//   Z/Z₀ = Σ_Γ (-λ)^n / |Aut Γ|,    log(Z/Z₀) = Σ_{Γ connected} (-λ)^n / |Aut Γ|
//
// (Itzykson and Zuber 1980, Quantum Field Theory, §6-1-3; Zinn-Justin 2002, §9.1; the
// symmetry factor of Γ is |Aut Γ|).
//
// Pairings are grouped by the labelled multigraph they produce, its adjacency matrix A
// (a_vv the number of self-loops at v, each using two legs): a pairing produces A in
// Π_v k! / (Π_v 2^{a_vv} a_vv! Π_{v<w} a_vw!) ways. Labelled multigraphs are enumerated
// directly, and grouped into isomorphism classes by their least relabelling (n! of them:
// brute force, fit for n <= 6). wick_pairings() enumerates the pairings themselves for the
// smallest orders, to check the grouping.

#pragma once

#include <algorithm>
#include <cstddef>
#include <functional>
#include <map>
#include <numeric>
#include <vector>

#include <morphe/arithmetic/rational.hpp>

namespace physis::field_theory {

using Adjacency = std::vector<std::vector<int>>;

struct FeynmanGraph {
  std::size_t order = 0;      // n, the number of vertices
  Adjacency adjacency;        // canonical; adjacency[v][v] counts self-loops at v
  morphe::arithmetic::Integer pairings;      // the Wick pairings that produce the graph
  morphe::arithmetic::Integer automorphisms; // |Aut Γ| = (k!)^n n! / pairings
  bool connected = false;
  std::size_t edges = 0, loops = 0;  // propagators, and independent cycles (E - n + components)
  morphe::arithmetic::Rational weight() const { return morphe::arithmetic::Rational(morphe::arithmetic::Integer(1), automorphisms); }
};

namespace detail {

inline std::size_t components(const Adjacency& a) {
  const std::size_t n = a.size();
  std::vector<std::size_t> p(n);
  std::iota(p.begin(), p.end(), 0);
  std::function<std::size_t(std::size_t)> root = [&](std::size_t x) { return p[x] == x ? x : p[x] = root(p[x]); };
  for (std::size_t v = 0; v < n; v++)
    for (std::size_t w = v + 1; w < n; w++)
      if (a[v][w]) p[root(v)] = root(w);
  std::size_t c = 0;
  for (std::size_t v = 0; v < n; v++) c += root(v) == v;
  return c;
}

inline Adjacency canonical(const Adjacency& a) {
  const std::size_t n = a.size();
  std::vector<std::size_t> perm(n);
  std::iota(perm.begin(), perm.end(), 0);
  Adjacency best;
  do {
    Adjacency b(n, std::vector<int>(n));
    for (std::size_t v = 0; v < n; v++)
      for (std::size_t w = 0; w < n; w++) b[v][w] = a[perm[v]][perm[w]];
    if (best.empty() || b < best) best = std::move(b);
  } while (std::next_permutation(perm.begin(), perm.end()));
  return best;
}

// The Wick pairings that produce the labelled multigraph a.
inline morphe::arithmetic::Integer pairings_of(const Adjacency& a, std::size_t k) {
  using morphe::arithmetic::Integer;
  const std::size_t n = a.size();
  Integer num = morphe::arithmetic::pow(morphe::arithmetic::factorial(k), n), den = 1;
  for (std::size_t v = 0; v < n; v++) {
    den *= morphe::arithmetic::pow(Integer(2), static_cast<unsigned long long>(a[v][v])) * morphe::arithmetic::factorial(static_cast<unsigned long long>(a[v][v]));
    for (std::size_t w = v + 1; w < n; w++) den *= morphe::arithmetic::factorial(static_cast<unsigned long long>(a[v][w]));
  }
  return morphe::arithmetic::exact_div(num, den);
}

}  // namespace detail

// Calls f(a) for every labelled vacuum multigraph of φ^k with n vertices.
inline void labelled_vacuum_graphs(std::size_t k, std::size_t n, const std::function<void(const Adjacency&)>& f) {
  Adjacency a(n, std::vector<int>(n, 0));
  std::vector<int> left(n, static_cast<int>(k));
  // Fill the upper triangle pair by pair: (0,0), (0,1), ..., (0,n-1), (1,1), ...
  std::function<void(std::size_t, std::size_t)> go = [&](std::size_t v, std::size_t w) {
    if (v == n) { f(a); return; }
    if (w == n) {
      if (left[v] == 0) go(v + 1, v + 1);
      return;
    }
    if (v == w) {
      for (int loops = 0; 2 * loops <= left[v]; loops++) {
        a[v][v] = loops, left[v] -= 2 * loops;
        go(v, w + 1);
        left[v] += 2 * loops, a[v][v] = 0;
      }
    } else {
      for (int m = 0; m <= std::min(left[v], left[w]); m++) {
        a[v][w] = a[w][v] = m, left[v] -= m, left[w] -= m;
        go(v, w + 1);
        left[v] += m, left[w] += m, a[v][w] = a[w][v] = 0;
      }
    }
  };
  if (n == 0) { f(a); return; }
  go(0, 0);
}

// The vacuum graphs of φ^k at order n, one per isomorphism class.
inline std::vector<FeynmanGraph> vacuum_graphs(std::size_t k, std::size_t n) {
  using morphe::arithmetic::Integer;
  std::map<Adjacency, Integer> classes;
  labelled_vacuum_graphs(k, n, [&](const Adjacency& a) { classes[detail::canonical(a)] += detail::pairings_of(a, k); });
  const Integer group = morphe::arithmetic::pow(morphe::arithmetic::factorial(k), n) * morphe::arithmetic::factorial(n);
  std::vector<FeynmanGraph> out;
  for (const auto& [a, count] : classes) {
    FeynmanGraph g;
    g.order = n;
    g.adjacency = a;
    g.pairings = count;
    g.automorphisms = morphe::arithmetic::exact_div(group, count);
    const std::size_t c = detail::components(a);
    g.connected = n > 0 && c == 1;
    g.edges = k * n / 2;
    g.loops = g.edges - n + c;
    out.push_back(std::move(g));
  }
  return out;
}

// Calls f(partner) for every perfect matching of the kn legs (leg v k + i of vertex v);
// (kn - 1)!! of them, so only for the smallest orders.
inline void wick_pairings(std::size_t k, std::size_t n, const std::function<void(const std::vector<std::size_t>&)>& f) {
  const std::size_t m = k * n;
  std::vector<std::size_t> partner(m, SIZE_MAX);
  std::function<void()> go = [&]() {
    std::size_t i = 0;
    while (i < m && partner[i] != SIZE_MAX) i++;
    if (i == m) { f(partner); return; }
    for (std::size_t j = i + 1; j < m; j++)
      if (partner[j] == SIZE_MAX) {
        partner[i] = j, partner[j] = i;
        go();
        partner[i] = partner[j] = SIZE_MAX;
      }
  };
  if (m % 2 == 0) go();
}

// The labelled multigraph of a pairing.
inline Adjacency graph_of(const std::vector<std::size_t>& partner, std::size_t k, std::size_t n) {
  Adjacency a(n, std::vector<int>(n, 0));
  for (std::size_t i = 0; i < partner.size(); i++)
    if (i < partner[i]) {
      const std::size_t v = i / k, w = partner[i] / k;
      if (v == w) a[v][v]++;
      else a[v][w]++, a[w][v]++;
    }
  return a;
}

}  // namespace physis::field_theory
