// Field theories of stranded vertices on a real field over a finite group: the interaction as a
// polynomial, its symmetries, and its Feynman graphs with their 2-complexes.
//
// The field. A real field of rank d on a finite group G is a function ψ: G^d → ℝ, |G|^d numbers.
// Its Gaussian measure exp(−½ Σ_x ψ(x)²) makes the values independent standard normals, so the
// propagator is ⟨ψ(x) ψ(y)⟩ = δ_{xy}: it identifies the d slots of two field factors one by one.
// This is the original, uncoloured form of the simplicial group field theories (Boulatov 1992,
// Mod. Phys. Lett. A 7, 1629, without the gauge projection; physis decision 0006 holds their
// coloured and tensorial forms, gft.hpp).
//
// The vertex. A stranded vertex is a monomial of k field factors whose kd slots are joined in
// pairs, every index summed once: V[ψ] = Σ_{indices ∈ G} Π_f ψ(slots of f). It is given by words,
// one per factor, each letter an index used exactly twice. Boulatov's tetrahedral vertex is
//
//   V_tet[ψ] = Σ_{a,b,c,d,e,f ∈ G} ψ(a,b,c) ψ(c,d,e) ψ(e,b,f) ψ(f,d,a):
//
// four triangles, any two sharing one edge, glued into a tetrahedron.
//
// Symmetries. Left translation of slot s by g_s, (Uψ)(x) = ψ(g_1⁻¹x_1, …, g_d⁻¹x_d), leaves V invariant
// when every index sits in slots translated alike, for then summing over g⁻¹·index is a bijection;
// for every G the invariance is exactly that when |G| > 1, by evaluating on a field supported on one
// point. So the slots fall into classes, those joined by some strand, and V is invariant under
// G^{#classes}. For the tetrahedral vertex the strands join slot 1 to slot 3 and slot 2 to slot 2:
// the classes are {1, 3} and {2}, and the group is G × G, of which the simultaneous translation of
// all three slots is the diagonal (the same holds for right translations). The gauge projection of
// group field theories averages over the simultaneous right translation,
// (Pψ)(x) = |G|⁻¹ Σ_h ψ(x_1 h, …, x_d h); for abelian G it coincides with the average over the
// simultaneous left translation, so it annihilates every functional odd under that translation.
//
// Feynman graphs. By Wick's theorem (Isserlis 1918), ⟨V^n⟩ is a sum over the (kn − 1)!! perfect
// matchings of the kn field factors of n copies of the vertex, the Wick pairings. Each pairing
// identifies slots; the strands then close into cycles, the faces, and every face carries a free
// sum over G:
//
//   ⟨V^n⟩ = Σ_{pairings} |G|^{F(pairing)},
//
// the identity check_stranded.cpp verifies against the moment computed by brute force. A pairing's
// 2-complex has a vertex for each copy of the interaction, an edge for each propagator and a 2-cell
// for each face, attached along the propagators it runs through; morphe's topology gives its
// homology, physis's bf.hpp its amplitude with the gauge projection. The pairings of n copies fall
// into classes under the group Aut(v) ≀ S_n, which permutes the copies and maps each copy by an
// automorphism of the vertex (a permutation of the factors carrying strands to strands, slot by
// slot); a class Γ holds |Aut v|^n n! / |Aut Γ| pairings, so that
//
//   Z/Z₀ = Σ_n ((−λ)^n / n!) ⟨V^n⟩ = Σ_Γ (−λ |Aut v|)^{n(Γ)} |G|^{F(Γ)} / |Aut Γ|
//
// as a formal series. A class is named by the least image of its pairings under the group, which is
// small (|Aut v|^n n!), and |Aut Γ| is the number of its elements that fix a pairing.

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <morphe/topology/invariants.hpp>

#include "bf.hpp"
#include "gft.hpp"

namespace physis::quantum_gravity {

struct StrandedVertex {
  std::string name;
  std::size_t factors = 0, rank = 0;
  std::vector<std::size_t> strand;   // strand[f·rank + s]: the slot (f'·rank + s') it is joined to
  std::vector<std::size_t> index;    // index[f·rank + s]: the number of the summed index there
  std::vector<std::string> words;    // the factors as words of index names

  std::size_t slots() const { return factors * rank; }
  std::size_t indices() const { return slots() / 2; }

  // From one word per factor, all of one length, each letter used exactly twice.
  static StrandedVertex from_words(std::string name, const std::vector<std::string>& words) {
    StrandedVertex v;
    v.name = std::move(name);
    v.words = words;
    v.factors = words.size();
    if (words.empty()) throw std::invalid_argument("physis::quantum_gravity::StrandedVertex: no factors");
    v.rank = words[0].size();
    std::map<char, std::vector<std::size_t>> where;
    for (std::size_t f = 0; f < words.size(); f++) {
      if (words[f].size() != v.rank) throw std::invalid_argument("physis::quantum_gravity::StrandedVertex: factors of different ranks");
      for (std::size_t s = 0; s < v.rank; s++) where[words[f][s]].push_back(f * v.rank + s);
    }
    v.strand.assign(v.slots(), 0);
    v.index.assign(v.slots(), 0);
    std::size_t k = 0;
    for (const auto& [letter, at] : where) {
      if (at.size() != 2) throw std::invalid_argument(std::string("physis::quantum_gravity::StrandedVertex: index ") + letter + " is not used exactly twice");
      v.strand[at[0]] = at[1];
      v.strand[at[1]] = at[0];
      v.index[at[0]] = v.index[at[1]] = k++;
    }
    return v;
  }
};

// Boulatov's tetrahedral vertex ψ(a,b,c)ψ(c,d,e)ψ(e,b,f)ψ(f,d,a).
inline StrandedVertex tetrahedral_vertex() { return StrandedVertex::from_words("tetrahedral", {"abc", "cde", "ebf", "fda"}); }
// The quartic melonic vertex of colour 1 in rank 3, ψ(a,b,c)ψ(d,b,c)ψ(d,e,f)ψ(a,e,f).
inline StrandedVertex melonic_vertex() { return StrandedVertex::from_words("quartic melonic", {"abc", "dbc", "def", "aef"}); }

// ---- fields on G^d ---------------------------------------------------------------------------------------------------------

// A point of G^d as a number, x_1 |G|^{d−1} + … + x_d.
inline std::size_t point(const std::vector<std::uint32_t>& x, std::size_t order) {
  std::size_t p = 0;
  for (auto g : x) p = p * order + g;
  return p;
}
inline std::vector<std::uint32_t> coordinates(std::size_t p, std::size_t order, std::size_t rank) {
  std::vector<std::uint32_t> x(rank);
  for (std::size_t s = rank; s-- > 0;) x[s] = static_cast<std::uint32_t>(p % order), p /= order;
  return x;
}
inline std::size_t field_size(std::size_t order, std::size_t rank) {
  std::size_t n = 1;
  for (std::size_t s = 0; s < rank; s++) n *= order;
  return n;
}

// V[ψ], summing every index over G; T is the field's value type (Rational, double).
template <class T>
T interaction(const StrandedVertex& v, const FiniteGroup& G, const std::vector<T>& psi) {
  const std::size_t q = G.order(), m = v.indices();
  if (psi.size() != field_size(q, v.rank)) throw std::invalid_argument("physis::quantum_gravity::interaction: a field of the wrong size");
  std::vector<std::uint32_t> value(m, 0), x(v.rank);
  T total{};
  for (;;) {
    T term = T(1);
    for (std::size_t f = 0; f < v.factors && term != T(0); f++) {
      for (std::size_t s = 0; s < v.rank; s++) x[s] = value[v.index[f * v.rank + s]];
      term = term * psi[point(x, q)];
    }
    total = total + term;
    std::size_t i = 0;
    while (i < m && ++value[i] == q) value[i++] = 0;
    if (i == m) break;
  }
  return total;
}

// The classes of slots joined by strands: V is invariant under independent translations of the
// classes.
inline std::vector<std::vector<std::size_t>> slot_classes(const StrandedVertex& v) {
  std::vector<std::size_t> p(v.rank);
  std::iota(p.begin(), p.end(), 0);
  std::function<std::size_t(std::size_t)> root = [&](std::size_t x) { return p[x] == x ? x : p[x] = root(p[x]); };
  for (std::size_t i = 0; i < v.slots(); i++) p[root(i % v.rank)] = root(v.strand[i] % v.rank);
  std::map<std::size_t, std::vector<std::size_t>> classes;
  for (std::size_t s = 0; s < v.rank; s++) classes[root(s)].push_back(s);
  std::vector<std::vector<std::size_t>> out;
  for (auto& [r, c] : classes) out.push_back(c);
  std::sort(out.begin(), out.end());  // by least slot
  return out;
}

// (Uψ)(x) = ψ(g_1⁻¹ x_1, …, g_d⁻¹ x_d): left translation of each slot.
template <class T>
std::vector<T> translate(const std::vector<T>& psi, const FiniteGroup& G, const std::vector<std::uint32_t>& g) {
  const std::size_t q = G.order(), d = g.size();
  std::vector<T> out(psi.size());
  for (std::size_t p = 0; p < psi.size(); p++) {
    auto x = coordinates(p, q, d);
    for (std::size_t s = 0; s < d; s++) x[s] = G.mul(G.inv(g[s]), x[s]);
    out[p] = psi[point(x, q)];
  }
  return out;
}

// The gauge projection (Pψ)(x) = |G|⁻¹ Σ_h ψ(x_1 h, …, x_d h).
template <class T>
std::vector<T> gauge_project(const std::vector<T>& psi, const FiniteGroup& G, std::size_t rank) {
  const std::size_t q = G.order();
  std::vector<T> out(psi.size());
  for (std::size_t p = 0; p < psi.size(); p++) {
    T s{};
    for (std::uint32_t h = 0; h < q; h++) {
      auto x = coordinates(p, q, rank);
      for (auto& g : x) g = G.mul(g, h);
      s = s + psi[point(x, q)];
    }
    out[p] = s / T(static_cast<long long>(q));
  }
  return out;
}

// ---- Wick pairings and their graphs ----------------------------------------------------------------------------------------

using Pairing = std::vector<std::size_t>;  // partner[i] for each field factor i = copy·k + factor

// Every perfect matching of m items, (m − 1)!! of them.
inline void for_each_pairing(std::size_t m, const std::function<void(const Pairing&)>& f) {
  if (m % 2) return;
  Pairing partner(m, m);
  std::function<void()> go = [&]() {
    std::size_t a = 0;
    while (a < m && partner[a] != m) a++;
    if (a == m) { f(partner); return; }
    for (std::size_t b = a + 1; b < m; b++) {
      if (partner[b] != m) continue;
      partner[a] = b, partner[b] = a;
      go();
      partner[a] = partner[b] = m;
    }
  };
  go();
}

struct StrandedGraph {
  std::size_t copies = 0;
  Pairing partner;
  std::vector<std::array<std::size_t, 2>> propagators;  // (factor, factor), the first the lesser
  CWComplex complex;                                    // its 2-complex
  bool connected = false;
  std::size_t faces() const { return complex.faces.size(); }
};

// The faces of a pairing: each closed strand, the propagators it crosses in order.
inline StrandedGraph stranded_graph(const StrandedVertex& v, std::size_t n, const Pairing& partner) {
  const std::size_t k = v.factors, d = v.rank;
  if (partner.size() != n * k) throw std::invalid_argument("physis::quantum_gravity::stranded_graph: a pairing of the wrong size");
  StrandedGraph S;
  S.copies = n;
  S.partner = partner;
  std::vector<std::size_t> edge_of(n * k);
  CWComplex& X = S.complex;
  X.name = v.name + " Feynman 2-complex";
  X.vertices = n;
  for (std::size_t i = 0; i < n * k; i++)
    if (i < partner[i]) {
      edge_of[i] = edge_of[partner[i]] = S.propagators.size();
      S.propagators.push_back({i, partner[i]});
      X.edges.push_back({i / k, partner[i] / k});
    }
  // Half-strands (factor i, slot s) as i·d + s; inside a copy the vertex's strands join them.
  std::vector<char> seen(n * k * d, 0);
  for (std::size_t start = 0; start < n * k * d; start++) {
    if (seen[start]) continue;
    morphe::topology::Face face;
    face.base = (start / d) / k;
    std::size_t x = start;
    do {
      seen[x] = 1;
      const std::size_t i = x / d, s = x % d, copy = i / k;
      const std::size_t local = v.strand[(i % k) * d + s];
      const std::size_t y = (copy * k + local / d) * d + local % d;  // across the vertex
      seen[y] = 1;
      const std::size_t j = y / d, t = y % d;
      const std::size_t e = edge_of[j];
      face.boundary.push_back({e, S.propagators[e][0] == j ? 1 : -1});  // along the propagator from j
      x = partner[j] * d + t;
    } while (x != start);
    X.faces.push_back(std::move(face));
  }
  X.validate();
  // Connected: the copies joined by propagators.
  std::vector<std::size_t> p(n);
  std::iota(p.begin(), p.end(), 0);
  std::function<std::size_t(std::size_t)> root = [&](std::size_t a) { return p[a] == a ? a : p[a] = root(p[a]); };
  for (const auto& e : X.edges) p[root(e[0])] = root(e[1]);
  std::size_t c = 0;
  for (std::size_t a = 0; a < n; a++) c += root(a) == a;
  S.connected = c <= 1;
  return S;
}

// The automorphisms of the vertex: permutations of its factors carrying strands to strands, slot
// by slot.
inline std::vector<std::vector<std::size_t>> vertex_automorphisms(const StrandedVertex& v) {
  std::vector<std::vector<std::size_t>> out;
  std::vector<std::size_t> sigma(v.factors);
  std::iota(sigma.begin(), sigma.end(), 0);
  do {
    bool ok = true;
    for (std::size_t i = 0; i < v.slots() && ok; i++) {
      const std::size_t j = v.strand[i];
      const std::size_t fi = sigma[i / v.rank] * v.rank + i % v.rank, fj = sigma[j / v.rank] * v.rank + j % v.rank;
      ok = v.strand[fi] == fj;
    }
    if (ok) out.push_back(sigma);
  } while (std::next_permutation(sigma.begin(), sigma.end()));
  return out;
}

// The elements of Aut(v) ≀ S_n as maps of the field factors: copy c's factor f goes to copy σ(c)'s
// factor α_c(f). Their number is |Aut v|^n n!, a few hundred for the orders computed here.
inline std::vector<std::vector<std::size_t>> wreath_maps(const StrandedVertex& v, std::size_t n) {
  const auto aut = vertex_automorphisms(v);
  const std::size_t k = v.factors;
  std::vector<std::vector<std::size_t>> out;
  std::vector<std::size_t> sigma(n);
  std::iota(sigma.begin(), sigma.end(), 0);
  do {
    std::vector<std::size_t> choice(n, 0);
    for (;;) {
      std::vector<std::size_t> map(n * k);
      for (std::size_t c = 0; c < n; c++)
        for (std::size_t f = 0; f < k; f++) map[c * k + f] = sigma[c] * k + aut[choice[c]][f];
      out.push_back(std::move(map));
      std::size_t i = 0;
      while (i < n && ++choice[i] == aut.size()) choice[i++] = 0;
      if (i == n) break;
    }
  } while (std::next_permutation(sigma.begin(), sigma.end()));
  return out;
}

struct FeynmanClass {
  Pairing representative;         // the least pairing of the class
  StrandedGraph graph;
  std::size_t pairings = 0;       // the Wick pairings in the class
  std::size_t automorphisms = 0;  // |Aut Γ|: the elements of Aut(v) ≀ S_n that fix the pairing
};

// The classes of the Wick pairings of n copies of the vertex, with their counts: each pairing's class
// is named by its least image under Aut(v) ≀ S_n, and its automorphisms are counted among them.
inline std::vector<FeynmanClass> feynman_classes(const StrandedVertex& v, std::size_t n) {
  const auto W = wreath_maps(v, n);
  std::map<Pairing, std::size_t> index;
  std::vector<FeynmanClass> out;
  Pairing image(n * v.factors);
  for_each_pairing(n * v.factors, [&](const Pairing& p) {
    Pairing least = p;
    std::size_t fixing = 0;
    for (const auto& map : W) {
      for (std::size_t i = 0; i < p.size(); i++) image[map[i]] = map[p[i]];
      if (image == p) fixing++;
      if (image < least) least = image;
    }
    const auto it = index.find(least);
    if (it != index.end()) { out[it->second].pairings++; return; }
    index.emplace(least, out.size());
    out.push_back({least, stranded_graph(v, n, least), 1, fixing});
  });
  return out;
}

// ⟨V^n⟩ for ψ(x) independent standard normals, by brute force: every assignment of all indices of
// the n copies, each contributing E[Π ψ] = Π_x (m_x − 1)!! over the points x met m_x times (zero if
// some m_x is odd).
inline Integer gaussian_moment(const StrandedVertex& v, const FiniteGroup& G, std::size_t n) {
  const std::size_t q = G.order(), m = v.indices() * n, size = field_size(q, v.rank);
  std::vector<std::uint32_t> value(m, 0), x(v.rank);
  std::vector<unsigned> mult(size);
  std::vector<Integer> dfact{1, 1};  // (j − 1)!! for j = 0, 1 (unused odd), 2, …
  Integer total = 0;
  for (;;) {
    std::fill(mult.begin(), mult.end(), 0u);
    for (std::size_t c = 0; c < n; c++)
      for (std::size_t f = 0; f < v.factors; f++) {
        for (std::size_t s = 0; s < v.rank; s++) x[s] = value[c * v.indices() + v.index[f * v.rank + s]];
        mult[point(x, q)]++;
      }
    Integer term = 1;
    for (auto k : mult) {
      if (k % 2) { term = 0; break; }
      while (dfact.size() <= k) dfact.push_back(dfact.size() % 2 ? dfact[dfact.size() - 1] : dfact[dfact.size() - 2] * Integer(static_cast<long long>(dfact.size() - 1)));
      term *= dfact[k];
    }
    total += term;
    std::size_t i = 0;
    while (i < m && ++value[i] == q) value[i++] = 0;
    if (i == m) break;
  }
  return total;
}

}  // namespace physis::quantum_gravity
