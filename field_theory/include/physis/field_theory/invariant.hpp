// Invariant field theories of rank r: their actions, Feynman expansions and exact series.
//
// A complex field φ of rank r has components φ_{a_1 ... a_r}, each index in a set of size N
// (rank 0 a scalar, 1 a vector, 2 a matrix, r >= 3 a tensor). A theory invariant under U(N)^r
// has the action (physis decision 0006, §2.3)
//
//   S[φ, φ̄] = φ̄·φ + Σ_B λ_B N^{−s_B} I_B(φ, φ̄),   Z = ∫ dφ dφ̄ e^{−S} / ∫ dφ dφ̄ e^{−φ̄·φ},
//
// with I_B the trace invariant of a bubble B (morphe topology/bubbles.hpp) and s_B its large-N
// scaling. With the Gaussian measure ⟨φ_a φ̄_b⟩ = δ_ab, Wick's theorem turns n_B copies of each
// bubble into a sum over the bijections σ₀ from their P black vertices to their P white ones,
// the propagators, each contraction a product of Kronecker deltas whose free index sums are
// the faces of colours (0, i) of the glued (r+1)-coloured graph:
//
//   ⟨Π_B I_B^{n_B}⟩ = Σ_{σ₀ ∈ S_P} N^{F₀(σ₀)},   F₀ = Σ_{i=1}^{r} F_{0i},
//
// so that Z = Σ_n Π_B (−λ_B N^{−s_B})^{n_B} / n_B! ⟨Π_B I_B^{n_B}⟩ exactly, a series whose
// coefficients are integer polynomials in N (Gurau 2016, Random Tensors, ch. 3; for matrices
// 't Hooft 1974 and Di Francesco, Ginsparg and Zinn-Justin 1995, §2). By the exponential
// formula the free energy log Z is the same sum over the σ₀ whose Feynman graph is connected
// (its bubbles linked by propagators): every σ₀ is a disjoint union of connected ones on
// labelled bubbles, and N^{F₀} is multiplicative.
//
// Isomorphism classes. The group H = Π_B S_{n_B} ≀ Aut(B) relabels the bubble copies and acts
// on each by its automorphisms; its orbits on the σ₀ are the unlabelled Feynman graphs Γ, and
// |Aut Γ| = |H| / |orbit|, so that Z = Σ_Γ Π_B (−λ_B)^{n_B} |Aut B|^{n_B} N^{F₀(Γ)} / |Aut Γ|.
// The orbits are found by union–find over the P! matchings (P <= 8), applying a generating
// set of H. The degree of Gurau (rank >= 3), or the genus of the ribbon graph (rank 2), is
// that of the glued coloured graph (morphe ColouredGraph::degree): with the melonic scaling
// the order λ^n carries N^{r − 2ω/(r−1)!} for connected graphs of quartic melonic theories.
//
// Default scalings: s_B = 0 in rank 0; p − 1 for a product of p traces in rank 1 (the vector
// model's 1/N); (r − 1)(p − 1) in rank r >= 2 (planar for matrices, melonic for tensors). For
// bubbles that are not melonic, such as the wheel K₃,₃, the optimal scaling is smaller
// (Bonzom, Delepouve and Rivasseau 2015); leading_power() finds the actual leading N-power.
//
// Cost: P! matchings, P <= 10 for the sums and P <= 8 for the classes.

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <morphe/arithmetic/rational.hpp>
#include <morphe/topology/bubbles.hpp>

namespace physis::field_theory {

using morphe::arithmetic::Integer;
using morphe::arithmetic::Rational;
using morphe::topology::ColouredGraph;

struct Interaction {
  std::string name;
  ColouredGraph bubble;  // r colours
  double coupling = 1;   // λ_B, for numerical evaluation; the series are in λ_B exactly
  long long scaling = 0; // s_B: λ_B enters as λ_B N^{−s_B}
};

struct InvariantTheory {
  std::string name;
  std::size_t rank = 0;
  std::vector<Interaction> interactions;

  void validate() const {
    if (interactions.empty()) throw std::invalid_argument("physis::field_theory::InvariantTheory: no interaction");
    for (const auto& I : interactions) {
      morphe::topology::validate_bubble(I.bubble);
      if (I.bubble.colours != rank) throw std::invalid_argument("physis::field_theory::InvariantTheory: a bubble of the wrong rank");
    }
  }
};

// The default large-N scaling of a bubble (see above).
inline long long default_scaling(const ColouredGraph& B) {
  const long long r = static_cast<long long>(B.colours), p = static_cast<long long>(B.p);
  if (r == 0) return 0;
  if (r == 1) return p - 1;
  return (r - 1) * (p - 1);
}

inline Interaction interaction(std::string name, ColouredGraph B) {
  const long long s = default_scaling(B);
  return Interaction{std::move(name), std::move(B), 1.0, s};
}

// ---- the classical models ----------------------------------------------------------------------------------------

// (φ̄φ)^k for a complex scalar.
inline InvariantTheory complex_scalar(std::size_t k = 2) {
  return {"complex scalar (φ̄φ)^" + std::to_string(k), 0, {interaction("(φ̄φ)^" + std::to_string(k), morphe::topology::scalar_bubble(k))}};
}
// (φ̄·φ)^k for a U(N) vector.
inline InvariantTheory vector_model(std::size_t k = 2) {
  return {"vector model (φ̄·φ)^" + std::to_string(k), 1, {interaction("(φ̄·φ)^" + std::to_string(k), morphe::topology::vector_bubble(k))}};
}
// tr (MM†)^k for a complex N × N matrix.
inline InvariantTheory matrix_model(std::size_t k = 2) {
  return {"matrix model tr (MM†)^" + std::to_string(k), 2, {interaction("tr (MM†)^" + std::to_string(k), morphe::topology::trace_bubble(k))}};
}
// The quartic melonic tensor model of rank r: one bubble of each colour, or of colour 1 only.
inline InvariantTheory quartic_melonic_model(std::size_t r, bool every_colour = true) {
  if (r < 2) throw std::invalid_argument("physis::field_theory::quartic_melonic_model: rank >= 2");
  InvariantTheory T{"quartic melonic tensor model, rank " + std::to_string(r), r, {}};
  for (std::size_t c = 0; c < (every_colour ? r : 1); c++)
    T.interactions.push_back(interaction("melon of colour " + std::to_string(c + 1), morphe::topology::quartic_melonic_bubble(r, c)));
  return T;
}
inline InvariantTheory wheel_model() { return {"sextic wheel K₃,₃, rank 3", 3, {interaction("K₃,₃", morphe::topology::wheel_bubble())}}; }
inline InvariantTheory necklace_model() { return {"sextic necklace, rank 3", 3, {interaction("necklace", morphe::topology::necklace_bubble())}}; }

// A theory by name: "scalar", "vector", "matrix", "melonic" (rank r, every colour), "melonic-1"
// (rank r, colour 1 only), "wheel", "necklace"; k the degree of the scalar, vector and matrix
// interactions.
inline InvariantTheory theory_by_name(const std::string& name, std::size_t rank = 3, std::size_t k = 2) {
  if (k < 1 || k > 6) throw std::invalid_argument("physis::field_theory::theory_by_name: 1 <= k <= 6");
  if (name == "scalar") return complex_scalar(k);
  if (name == "vector") return vector_model(k);
  if (name == "matrix") return matrix_model(k);
  if (name == "melonic") return quartic_melonic_model(rank, true);
  if (name == "melonic-1") return quartic_melonic_model(rank, false);
  if (name == "wheel") return wheel_model();
  if (name == "necklace") return necklace_model();
  throw std::invalid_argument("physis::field_theory::theory_by_name: no theory \"" + name + "\"");
}

// ---- the expansion -------------------------------------------------------------------------------------------------

using Polynomial = std::map<long long, Rational>;  // a polynomial in N: power -> coefficient

struct FeynmanClass {
  std::vector<std::size_t> sigma0;  // a representative matching (the least in lexicographic order)
  ColouredGraph graph;              // its (r+1)-coloured graph, colour 0 the propagators
  std::size_t labelled = 0;         // matchings in the class
  Integer automorphisms;            // |Aut Γ| = |H| / labelled
  long long faces = 0;              // F₀: the power of N
  long long scaled = 0;             // F₀ − Σ_B n_B s_B
  std::size_t components = 0;       // connected components of the Feynman graph
  bool connected = false;
  std::size_t degree = 0;           // Gurau's degree (rank >= 3) or the genus (rank 2); 0 below
};

struct Expansion {
  std::vector<std::size_t> counts;           // n_B for each interaction
  std::vector<ColouredGraph> bubbles;        // the bubble copies, in order
  std::vector<std::size_t> type;             // the interaction of each copy
  std::vector<std::size_t> black_bubble;     // the copy holding each black vertex
  std::vector<std::size_t> white_bubble;     // the copy holding each white vertex
  long long scaling = 0;                     // Σ_B n_B s_B
  Integer relabellings;                      // |H| = Π_B n_B! |Aut B|^{n_B}
  std::map<long long, Integer> all;          // Σ over all matchings of N^{F₀}, by power
  std::map<long long, Integer> connected;    // the same over connected matchings
  std::vector<FeynmanClass> classes;         // when computed (P <= 8)
  bool classes_computed = false;
};

namespace detail {

inline std::size_t cycles(const std::vector<std::size_t>& perm) {
  std::vector<char> seen(perm.size(), 0);
  std::size_t c = 0;
  for (std::size_t i = 0; i < perm.size(); i++) {
    if (seen[i]) continue;
    c++;
    for (std::size_t j = i; !seen[j]; j = perm[j]) seen[j] = 1;
  }
  return c;
}

// The rank of a permutation of {0..n-1} in lexicographic order (its Lehmer code).
inline std::uint32_t permutation_rank(const std::vector<std::size_t>& p) {
  const std::size_t n = p.size();
  std::uint32_t rank = 0;
  for (std::size_t i = 0; i < n; i++) {
    std::uint32_t smaller = 0;
    for (std::size_t j = i + 1; j < n; j++) smaller += p[j] < p[i];
    rank = rank * static_cast<std::uint32_t>(n - i) + smaller;
  }
  return rank;
}

// A generating set of a group of pairs (π, ρ) given by all its elements, by greedy closure.
inline std::vector<std::vector<std::size_t>> generators(const std::vector<std::vector<std::size_t>>& elements, std::size_t p) {
  auto compose = [p](const std::vector<std::size_t>& a, const std::vector<std::size_t>& b) {
    std::vector<std::size_t> c(2 * p);
    for (std::size_t i = 0; i < p; i++) c[i] = a[b[i]], c[p + i] = a[p + b[p + i]];
    return c;
  };
  std::set<std::vector<std::size_t>> span;
  std::vector<std::size_t> id(2 * p);
  for (std::size_t i = 0; i < p; i++) id[i] = i, id[p + i] = i;
  span.insert(id);
  std::vector<std::vector<std::size_t>> gens;
  for (const auto& e : elements) {
    if (span.count(e)) continue;
    gens.push_back(e);
    std::vector<std::vector<std::size_t>> frontier(span.begin(), span.end());
    while (!frontier.empty()) {
      std::vector<std::vector<std::size_t>> next;
      for (const auto& x : frontier)
        for (const auto& g : gens) {
          auto y = compose(g, x);
          if (span.insert(y).second) next.push_back(std::move(y));
        }
      frontier = std::move(next);
    }
  }
  return gens;
}

}  // namespace detail

// The Feynman expansion at n_B copies of each interaction B.
inline Expansion expand(const InvariantTheory& T, const std::vector<std::size_t>& counts, bool classes = true) {
  T.validate();
  if (counts.size() != T.interactions.size()) throw std::invalid_argument("physis::field_theory::expand: one count per interaction");
  Expansion E;
  E.counts = counts;
  E.relabellings = 1;
  std::vector<std::vector<std::vector<std::size_t>>> automorphisms(T.interactions.size());
  for (std::size_t t = 0; t < T.interactions.size(); t++) {
    const auto& B = T.interactions[t].bubble;
    automorphisms[t] = morphe::topology::bubble_automorphisms(B);
    E.relabellings *= morphe::arithmetic::factorial(counts[t]) *
                      morphe::arithmetic::pow(Integer(static_cast<long long>(automorphisms[t].size())), counts[t]);
    E.scaling += static_cast<long long>(counts[t]) * T.interactions[t].scaling;
    for (std::size_t k = 0; k < counts[t]; k++) E.bubbles.push_back(B), E.type.push_back(t);
  }
  std::size_t P = 0;
  std::vector<std::size_t> offset;
  for (std::size_t k = 0; k < E.bubbles.size(); k++) {
    offset.push_back(P);
    for (std::size_t b = 0; b < E.bubbles[k].p; b++) E.black_bubble.push_back(k), E.white_bubble.push_back(k);
    P += E.bubbles[k].p;
  }
  if (P == 0) {
    E.all[0] = 1;
    return E;
  }
  if (P > 10) throw std::invalid_argument("physis::field_theory::expand: more than 10 propagators (10! matchings)");
  const std::size_t r = T.rank;
  // The bubbles' colours as global permutations, inverted: inv[i][white] = black.
  std::vector<std::vector<std::size_t>> inv(r, std::vector<std::size_t>(P));
  for (std::size_t k = 0; k < E.bubbles.size(); k++)
    for (std::size_t i = 0; i < r; i++)
      for (std::size_t b = 0; b < E.bubbles[k].p; b++) inv[i][offset[k] + E.bubbles[k].sigma[i][b]] = offset[k] + b;
  auto faces = [&](const std::vector<std::size_t>& s0) {
    long long f = 0;
    std::vector<std::size_t> perm(P);
    for (std::size_t i = 0; i < r; i++) {
      for (std::size_t b = 0; b < P; b++) perm[b] = inv[i][s0[b]];
      f += static_cast<long long>(detail::cycles(perm));
    }
    return f;
  };
  const std::size_t n = E.bubbles.size();
  auto components = [&](const std::vector<std::size_t>& s0) {
    std::vector<std::size_t> parent(n);
    std::iota(parent.begin(), parent.end(), 0);
    auto root = [&](std::size_t x) {
      while (parent[x] != x) x = parent[x] = parent[parent[x]];
      return x;
    };
    std::size_t c = n;
    for (std::size_t b = 0; b < P; b++) {
      const std::size_t a = root(E.black_bubble[b]), w = root(E.white_bubble[s0[b]]);
      if (a != w) parent[a] = w, c--;
    }
    return c;
  };
  const bool want_classes = classes && P <= 8;
  std::vector<std::uint32_t> uf;
  std::vector<std::vector<std::size_t>> gens;  // relabellings (π on blacks, ρ on whites), length 2P
  if (want_classes) {
    std::size_t total = 1;
    for (std::size_t k = 2; k <= P; k++) total *= k;
    uf.resize(total);
    std::iota(uf.begin(), uf.end(), 0u);
    auto identity = [P] {
      std::vector<std::size_t> g(2 * P);
      for (std::size_t i = 0; i < P; i++) g[i] = i, g[P + i] = i;
      return g;
    };
    // Swaps of adjacent copies of one interaction.
    for (std::size_t k = 0; k + 1 < n; k++)
      if (E.type[k] == E.type[k + 1]) {
        auto g = identity();
        for (std::size_t b = 0; b < E.bubbles[k].p; b++) {
          std::swap(g[offset[k] + b], g[offset[k + 1] + b]);
          std::swap(g[P + offset[k] + b], g[P + offset[k + 1] + b]);
        }
        gens.push_back(g);
      }
    // Automorphisms of the first copy of each interaction.
    for (std::size_t t = 0; t < T.interactions.size(); t++) {
      std::size_t first = n;
      for (std::size_t k = 0; k < n && first == n; k++)
        if (E.type[k] == t) first = k;
      if (first == n) continue;
      const std::size_t p = E.bubbles[first].p;
      for (const auto& a : detail::generators(automorphisms[t], p)) {
        auto g = identity();
        for (std::size_t b = 0; b < p; b++) g[offset[first] + b] = offset[first] + a[b], g[P + offset[first] + b] = offset[first] + a[p + b];
        gens.push_back(g);
      }
    }
  }
  std::function<std::uint32_t(std::uint32_t)> find = [&](std::uint32_t x) {
    while (uf[x] != x) x = uf[x] = uf[uf[x]];
    return x;
  };
  std::vector<std::size_t> s0(P), image(P);
  std::iota(s0.begin(), s0.end(), 0);
  std::uint32_t index = 0;
  do {
    const long long f = faces(s0);
    E.all[f] += 1;
    if (components(s0) == 1) E.connected[f] += 1;
    if (want_classes) {
      // σ₀' = ρ σ₀ π⁻¹: σ₀'(π(b)) = ρ(σ₀(b)).
      for (const auto& g : gens) {
        for (std::size_t b = 0; b < P; b++) image[g[b]] = g[P + s0[b]];
        const std::uint32_t a = find(index), c = find(detail::permutation_rank(image));
        if (a != c) uf[std::max(a, c)] = std::min(a, c);
      }
    }
    index++;
  } while (std::next_permutation(s0.begin(), s0.end()));
  if (!want_classes) return E;
  // The classes: their sizes, and the least matching of each as representative.
  std::map<std::uint32_t, std::size_t> size;
  for (std::uint32_t i = 0; i < uf.size(); i++) size[find(i)]++;
  std::iota(s0.begin(), s0.end(), 0);
  index = 0;
  do {
    const std::uint32_t root = find(index);
    if (root == index) {
      FeynmanClass C;
      C.sigma0 = s0;
      C.graph = morphe::topology::glue(E.bubbles, s0);
      C.labelled = size[root];
      C.automorphisms = morphe::arithmetic::exact_div(E.relabellings, Integer(static_cast<unsigned long long>(C.labelled)));
      C.faces = faces(s0);
      C.scaled = C.faces - E.scaling;
      C.components = components(s0);
      C.connected = C.components == 1;
      if (r >= 2) C.degree = C.graph.degree();
      E.classes.push_back(std::move(C));
    }
    index++;
  } while (std::next_permutation(s0.begin(), s0.end()));
  E.classes_computed = true;
  return E;
}

// ---- series -----------------------------------------------------------------------------------------------------------

inline Polynomial operator*(const Polynomial& a, const Polynomial& b) {
  Polynomial c;
  for (const auto& [i, x] : a)
    for (const auto& [j, y] : b) c[i + j] += x * y;
  for (auto it = c.begin(); it != c.end();) it = it->second.is_zero() ? c.erase(it) : std::next(it);
  return c;
}
inline Polynomial& operator+=(Polynomial& a, const Polynomial& b) {
  for (const auto& [i, x] : b) a[i] += x;
  for (auto it = a.begin(); it != a.end();) it = it->second.is_zero() ? a.erase(it) : std::next(it);
  return a;
}
inline Polynomial scaled(const Polynomial& a, const Rational& c, long long shift = 0) {
  Polynomial b;
  if (c.is_zero()) return b;
  for (const auto& [i, x] : a) b[i + shift] = x * c;
  return b;
}
// The term of highest power: its power and coefficient (0, 0 for the zero polynomial).
inline std::pair<long long, Rational> leading(const Polynomial& a) {
  if (a.empty()) return {0, Rational(0)};
  return {a.rbegin()->first, a.rbegin()->second};
}
inline Rational evaluate(const Polynomial& a, const Rational& N) {
  Rational s = 0;
  for (const auto& [i, x] : a) s += x * morphe::arithmetic::pow(N, i);
  return s;
}

struct Series {
  std::vector<Polynomial> z;            // Z = Σ_n z_n (−λ)^n
  std::vector<Polynomial> free_energy;  // log Z = Σ_n w_n (−λ)^n, from the logarithm of z
  std::vector<Polynomial> connected;    // w_n from the connected matchings (equal to the above)
  long long scaling = 0;                // s_B: with λ N^{−s_B} the order n is shifted by −n s_B
};

// The series of a theory with one interaction, to order n_max: in the coupling λ itself (put
// λ N^{−s} by shifting the powers of order n by −n s).
inline Series series(const InvariantTheory& T, std::size_t order) {
  if (T.interactions.size() != 1) throw std::invalid_argument("physis::field_theory::series: one interaction");
  Series S;
  S.scaling = T.interactions[0].scaling;
  for (std::size_t n = 0; n <= order; n++) {
    const Expansion E = expand(T, {n}, false);
    const Rational inv(Integer(1), morphe::arithmetic::factorial(n));
    Polynomial z, c;
    for (const auto& [f, m] : E.all) z[f] = Rational(m) * inv;
    for (const auto& [f, m] : E.connected) c[f] = Rational(m) * inv;
    S.z.push_back(z);
    S.connected.push_back(n == 0 ? Polynomial{} : c);
  }
  // log z: n w_n = n z_n − Σ_{k<n} k w_k z_{n−k}.
  S.free_energy.assign(order + 1, Polynomial{});
  for (std::size_t n = 1; n <= order; n++) {
    Polynomial s = scaled(S.z[n], Rational(static_cast<long long>(n)));
    for (std::size_t k = 1; k < n; k++) s += scaled(S.free_energy[k] * S.z[n - k], Rational(-static_cast<long long>(k)));
    S.free_energy[n] = scaled(s, Rational(Integer(1), Integer(static_cast<long long>(n))));
  }
  return S;
}

// The highest N-power, after scaling, of the connected graphs at n copies of each interaction.
inline long long leading_power(const InvariantTheory& T, const std::vector<std::size_t>& counts) {
  const Expansion E = expand(T, counts, false);
  if (E.connected.empty()) throw std::invalid_argument("physis::field_theory::leading_power: no connected graph");
  return E.connected.rbegin()->first - E.scaling;
}

}  // namespace physis::field_theory
