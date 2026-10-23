// Random discrete spacetimes from the large-N limits: random planar quadrangulations, the
// planar Feynman graphs of the quartic matrix model and the discrete surfaces of
// two-dimensional quantum gravity; and random melonic graphs, the leading graphs of tensor
// models, whose continuum is the branched polymer.
//
// Quadrangulations (Cori and Vauquelin 1981; Schaeffer 1998; Chassaing and Schaeffer 2004,
// Probab. Theory Relat. Fields 128, 161–212; Le Gall 2013, Ann. Probab. 41, 2880–2960). A
// plane tree with n edges whose vertices carry integer labels, the root's 0 and each child's
// its parent's plus −1, 0 or 1, is read around its contour as 2n corners c_0, …, c_{2n−1}.
// With the labels shifted so that their minimum is 1, and an extra vertex ∂ labelled 0, each
// corner is joined by an arc to its successor: the first corner after it, cyclically, whose
// label is one less, or ∂ when its label is 1. The arcs form a planar quadrangulation with n
// faces, 2n edges and n + 2 vertices, each face containing one edge of the tree, and the label
// of a vertex is its graph distance from ∂. Uniform labelled trees give uniform rooted
// quadrangulations (the bijection is 2-to-1 onto rooted, pointed quadrangulations, all of
// which have n + 2 vertices). Their distances scale as n^{1/4}: the Hausdorff dimension of the
// Brownian map, the continuum of 2D quantum gravity, is 4.
//
// The faces are traced from the planar embedding of the arcs: the region outside the tree is
// a disc bounded by the contour (corners in counterclockwise order); arcs are non-crossing
// chords of it, and those to ∂ radii, which separate the chords (a chord never passes a corner
// of minimum label, labels changing by at most one along the contour). At a corner, the arcs
// are ordered counterclockwise by the position of their other end after the corner, the arc to
// ∂ first; around a tree vertex, its corners come in decreasing contour order; around ∂, the
// corners of minimum label in increasing order. A face is an orbit of φ = σ∘α (σ the rotation,
// α the reversal of a half-edge).
//
// Melonic graphs (Bonzom, Gurau, Riello and Rivasseau 2011, Nucl. Phys. B 853, 174–195). The
// rooted melonic (D+1)-coloured graphs with 2p vertices are in bijection with the
// (D+1)-ary trees with p nodes: the root is the elementary melon, and a child in slot c of a
// node is a melon inserted on the edge of colour c at that node's black vertex. Uniform trees
// come from the cycle lemma (Dvoretzky and Motzkin 1947): of the cyclic rotations of a word
// with p letters +D and Dp + 1 letters −1, exactly one is a Łukasiewicz word. Their number is
// the Fuss–Catalan number ((D+1)p)! / (p! (Dp + 1)!), and in the continuum they are branched
// polymers (Hausdorff dimension 2, spectral dimension 4/3; Gurau and Ryan 2014).

#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

#include <calx/random/stream.hpp>
#include <morphe/topology/coloured_graph.hpp>
#include <morphe/topology/cw_complex.hpp>

namespace physis::quantum_gravity {

struct Quadrangulation {
  morphe::topology::CWComplex complex;  // tree vertices 0..n in preorder, ∂ = n + 1; 2n edges; n faces
  std::vector<int> label;               // per vertex: its distance from ∂ (0 at ∂)
  std::size_t extra = 0;                // the vertex ∂
};

// The quadrangulation of a labelled plane tree given by its contour (a Dyck word of +1 and −1,
// length 2n) and the label increments of its non-root vertices in preorder (each −1, 0 or 1).
inline Quadrangulation quadrangulation(const std::vector<int>& dyck, const std::vector<int>& increments) {
  const std::size_t two_n = dyck.size(), n = two_n / 2;
  if (n == 0 || two_n % 2 || increments.size() != n) throw std::invalid_argument("physis::quantum_gravity::quadrangulation: a Dyck word of length 2n and n increments");
  // The tree: vertices in preorder, the contour's vertices u_0 .. u_{2n}.
  std::vector<std::size_t> parent(n + 1, SIZE_MAX), contour;
  std::vector<int> label(n + 2, 0);
  std::size_t at = 0, next = 1;
  int height = 0;
  for (int step : dyck) {
    contour.push_back(at);
    if (step > 0) {
      parent[next] = at;
      const int d = increments[next - 1];
      if (d < -1 || d > 1) throw std::invalid_argument("physis::quantum_gravity::quadrangulation: increments are -1, 0 or 1");
      label[next] = label[at] + d;
      at = next++;
      height++;
    } else {
      if (height == 0) throw std::invalid_argument("physis::quantum_gravity::quadrangulation: not a Dyck word");
      at = parent[at];
      height--;
    }
  }
  if (height != 0 || next != n + 1) throw std::invalid_argument("physis::quantum_gravity::quadrangulation: not a Dyck word");
  const int lo = *std::min_element(label.begin(), label.begin() + static_cast<long>(n + 1));
  for (std::size_t v = 0; v <= n; v++) label[v] = label[v] - lo + 1;
  const std::size_t extra = n + 1;
  label[extra] = 0;
  // Corner labels and successors (by the next occurrence of each label on the doubled contour).
  std::vector<int> cl(two_n);
  for (std::size_t i = 0; i < two_n; i++) cl[i] = label[contour[i]];
  const int top = *std::max_element(cl.begin(), cl.end());
  std::vector<std::size_t> next_at(static_cast<std::size_t>(top + 2), SIZE_MAX), succ(two_n, SIZE_MAX);
  for (std::size_t k = 2 * two_n; k-- > 0;) {
    const std::size_t i = k % two_n;
    if (k < two_n && cl[i] > 1) succ[i] = next_at[static_cast<std::size_t>(cl[i] - 1)] % two_n;
    next_at[static_cast<std::size_t>(cl[i])] = k;
  }
  // Half-edge 2e + end: arc e = corner e → succ, end 0 at its tail, end 1 at its head.
  Quadrangulation Q;
  Q.label = label;
  Q.extra = extra;
  auto& X = Q.complex;
  X.name = "random quadrangulation with " + std::to_string(n) + " faces";
  X.vertices = n + 2;
  for (std::size_t e = 0; e < two_n; e++) X.edges.push_back({contour[e], succ[e] == SIZE_MAX ? extra : contour[succ[e]]});
  // The arcs at each corner, ordered by the position of their other end after the corner.
  std::vector<std::vector<std::pair<std::size_t, std::size_t>>> at_corner(two_n);  // (key, half-edge)
  auto key = [&](std::size_t i, std::size_t p) { return (p + two_n - i - 1) % two_n + 1; };
  std::vector<std::size_t> minima;
  for (std::size_t e = 0; e < two_n; e++) {
    if (succ[e] == SIZE_MAX) {
      at_corner[e].push_back({0, 2 * e});
      minima.push_back(e);
    } else {
      at_corner[e].push_back({key(e, succ[e]), 2 * e});
      at_corner[succ[e]].push_back({key(succ[e], e), 2 * e + 1});
    }
  }
  // The rotation at each vertex: its corners in decreasing contour order; at ∂ the minima in
  // increasing order.
  std::vector<std::vector<std::size_t>> corners_of(n + 1);
  for (std::size_t i = 0; i < two_n; i++) corners_of[contour[i]].push_back(i);
  std::vector<std::size_t> sigma(2 * two_n);
  auto close = [&](const std::vector<std::size_t>& rot) {
    for (std::size_t k = 0; k < rot.size(); k++) sigma[rot[k]] = rot[(k + 1) % rot.size()];
  };
  for (std::size_t v = 0; v <= n; v++) {
    std::vector<std::size_t> rot;
    for (auto it = corners_of[v].rbegin(); it != corners_of[v].rend(); ++it) {
      auto items = at_corner[*it];
      std::sort(items.begin(), items.end());
      for (const auto& [k, h] : items) rot.push_back(h);
    }
    close(rot);
  }
  {
    std::vector<std::size_t> rot;
    for (auto e : minima) rot.push_back(2 * e + 1);
    close(rot);
  }
  // Faces: orbits of φ = σ ∘ α.
  std::vector<char> seen(2 * two_n, 0);
  for (std::size_t h0 = 0; h0 < 2 * two_n; h0++) {
    if (seen[h0]) continue;
    morphe::topology::Face f;
    f.base = h0 % 2 ? X.edges[h0 / 2][1] : X.edges[h0 / 2][0];
    for (std::size_t h = h0; !seen[h]; h = sigma[h ^ 1]) {
      seen[h] = 1;
      f.boundary.push_back({h / 2, h % 2 ? -1 : 1});
    }
    X.faces.push_back(f);
  }
  return Q;
}

// A uniform random quadrangulation with n faces: a uniform Dyck word (by the cycle lemma on n
// up-steps and n + 1 down-steps) and uniform increments.
inline Quadrangulation random_quadrangulation(std::size_t n, std::uint64_t seed) {
  if (n == 0) throw std::invalid_argument("physis::quantum_gravity::random_quadrangulation: n >= 1");
  calx::random::Stream z(seed, 0x71756164u, 0, 0);  // "quad"
  std::vector<int> w(2 * n + 1, -1);
  for (std::size_t i = 0; i < n; i++) w[i] = 1;
  for (std::size_t i = w.size(); i > 1; i--) std::swap(w[i - 1], w[z.below(i)]);
  // The rotation starting after the first minimum of the partial sums ends with −1 and stays
  // nonnegative before; dropping that last −1 leaves a Dyck word.
  int s = 0, lo = 1;
  std::size_t at = 0;
  for (std::size_t k = 0; k < w.size(); k++) {
    s += w[k];
    if (s < lo) lo = s, at = k;
  }
  std::vector<int> dyck;
  for (std::size_t k = 1; k < w.size(); k++) dyck.push_back(w[(at + k) % w.size()]);
  std::vector<int> inc(n);
  for (auto& d : inc) d = static_cast<int>(z.below(3)) - 1;
  return quadrangulation(dyck, inc);
}

// The rooted melonic (D+1)-coloured graph of a (D+1)-ary tree given as a Łukasiewicz word in
// preorder: 1 for a node (followed by its D + 1 slots), 0 for an empty slot.
inline morphe::topology::ColouredGraph melonic_graph(std::size_t D, const std::vector<int>& word) {
  using morphe::topology::ColouredGraph;
  if (word.empty() || word[0] != 1) throw std::invalid_argument("physis::quantum_gravity::melonic_graph: the word starts with a node");
  ColouredGraph G = morphe::topology::elementary_melon(D);
  // Stack of (black vertex of the node, next slot to fill).
  std::vector<std::pair<std::size_t, std::size_t>> stack{{0, 0}};
  for (std::size_t k = 1; k < word.size(); k++) {
    if (stack.empty()) throw std::invalid_argument("physis::quantum_gravity::melonic_graph: the word is longer than its tree");
    auto& [b, slot] = stack.back();
    const std::size_t c = slot++;
    const std::size_t parent = b;
    if (stack.back().second == D + 1) stack.pop_back();
    if (word[k] == 1) {
      const std::size_t q = G.p;  // the new melon's black vertex
      G = morphe::topology::insert_melon(G, c, parent);
      stack.push_back({q, 0});
    }
  }
  if (!stack.empty()) throw std::invalid_argument("physis::quantum_gravity::melonic_graph: the word ends before its tree");
  return G;
}

// A uniform random rooted melonic (D+1)-coloured graph with 2p vertices.
inline morphe::topology::ColouredGraph random_melonic_graph(std::size_t D, std::size_t p, std::uint64_t seed) {
  if (p == 0 || D < 1) throw std::invalid_argument("physis::quantum_gravity::random_melonic_graph: p >= 1 and D >= 1");
  calx::random::Stream z(seed, 0x6d656c6fu, 0, 0);  // "melo"
  const std::size_t m = (D + 1) * p + 1;
  std::vector<int> w(m, -1);
  for (std::size_t i = 0; i < p; i++) w[i] = static_cast<int>(D);
  for (std::size_t i = m; i > 1; i--) std::swap(w[i - 1], w[z.below(i)]);
  long s = 0, lo = 1;
  std::size_t at = 0;
  for (std::size_t k = 0; k < m; k++) {
    s += w[k];
    if (s < lo) lo = s, at = k;
  }
  std::vector<int> word;
  for (std::size_t k = 1; k <= m; k++) word.push_back(w[(at + k) % m] > 0 ? 1 : 0);
  return melonic_graph(D, word);
}

}  // namespace physis::quantum_gravity
