// The dimensions of discrete spaces: spectral, by diffusion, and Hausdorff, by the growth of
// balls. They are how a continuum is recognised in a discrete spacetime, and how candidates
// differ: causal dynamical triangulations have spectral dimension 4 at large scales and about 2
// at small ones (Ambjørn, Jurkiewicz and Loll 2005, Phys. Rev. Lett. 95, 171301); the branched
// polymers of melonic tensor models have d_s = 4/3 and d_H = 2 (Jonsson and Wheater 1998,
// Nucl. Phys. B 515, 549; Gurau and Ryan 2014, Ann. Henri Poincaré 15, 2085); random planar
// maps, 2D quantum gravity, have d_H = 4 and d_s = 2 (Ambjørn and Watabiki 1995; Le Gall
// 2013); a d-dimensional lattice has d for both.
//
// A discrete space is a graph: the 1-skeleton of a complex, the dual graph of its top cells
// (adjacent across a shared face of one dimension less), or a coloured graph itself.
//
// - **Spectral dimension.** The lazy random walk stays with probability 1/2 and otherwise
//   moves along an edge chosen uniformly among those at its vertex; laziness removes the
//   parity of bipartite graphs. Its return probability P(t), averaged over starting vertices,
//   falls on a d-dimensional space as t^{−d/2} for 1 ≪ t ≪ L² (diffusion; Durhuus, Jonsson and
//   Wheater 2006, J. Stat. Phys. 128, 1237), so
//     d_s(t) = −2 d ln P / d ln t,
//   estimated by centred differences over the window [t − w, t + w], w = max(1, ⌊window · t⌋),
//   which smooths the oscillations of small t. P is computed exactly: the distribution is
//   propagated from each starting vertex (p_{t+1}(w) = p_t(w)/2 + Σ_{v ∼ w} p_t(v)/(2 deg v)),
//   in O(t (V + E)) for each start; the starts are all vertices, or a sample drawn from a calx
//   stream. On a finite graph P(t) → deg-weighted 1/V and d_s → 0: the window of scales where
//   d_s is flat is the dimension's evidence.
// - **Hausdorff dimension.** The mean number of vertices within graph distance r of a centre,
//   ⟨|B(r)|⟩ ∼ r^{d_H}, by breadth-first search from each centre in O(V + E):
//     d_H(r) = d ln⟨|B(r)|⟩ / d ln r,  by centred differences over r ± 1.
//
// Both are estimates at finite size, to be read where they are flat; the checks compare the
// exact return probability on periodic lattices with its Fourier sum.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <vector>

#include <calx/random/stream.hpp>
#include <morphe/combinatorics/graph.hpp>
#include <morphe/topology/coloured_graph.hpp>
#include <morphe/topology/cw_complex.hpp>

namespace physis::spacetime {

// A graph by compressed rows: the neighbours of v are neighbour[offset[v] .. offset[v + 1]),
// a loop appearing twice (both its ends), parallel edges each time.
struct Graph {
  std::vector<std::size_t> offset{0}, neighbour;
  std::size_t vertices() const { return offset.size() - 1; }
  std::size_t degree(std::size_t v) const { return offset[v + 1] - offset[v]; }

  static Graph from_edges(std::size_t n, const std::vector<std::array<std::size_t, 2>>& edges) {
    Graph g;
    std::vector<std::size_t> deg(n, 0);
    for (const auto& e : edges) {
      if (e[0] >= n || e[1] >= n) throw std::out_of_range("physis::spacetime::Graph: an edge leaves the vertices");
      deg[e[0]]++, deg[e[1]]++;
    }
    g.offset.assign(n + 1, 0);
    for (std::size_t v = 0; v < n; v++) g.offset[v + 1] = g.offset[v] + deg[v];
    g.neighbour.resize(g.offset[n]);
    std::vector<std::size_t> fill(g.offset.begin(), g.offset.end() - 1);
    for (const auto& e : edges) {
      g.neighbour[fill[e[0]]++] = e[1];
      g.neighbour[fill[e[1]]++] = e[0];
    }
    return g;
  }
  static Graph from(const morphe::combinatorics::Multigraph& G) { return from_edges(G.vertices, G.edges); }
};

// The 1-skeleton of a complex.
inline Graph skeleton_graph(const morphe::topology::CWComplex& X) { return Graph::from_edges(X.vertices, X.edges); }

// The dual graph of a complex's top cells: for a 2-complex, its 2-cells, joined once for each
// pair of occurrences of an edge in different 2-cells' words (or in one word twice: a loop);
// for higher dimension, the top cells, joined across each (d−1)-cell for each pair of top cells
// whose boundary meets it.
inline Graph dual_graph(const morphe::topology::CWComplex& X) {
  const std::size_t d = X.dimension();
  std::vector<std::array<std::size_t, 2>> edges;
  auto join = [&](const std::vector<std::size_t>& holders) {
    for (std::size_t a = 0; a < holders.size(); a++)
      for (std::size_t b = a + 1; b < holders.size(); b++) edges.push_back({holders[a], holders[b]});
  };
  if (d < 2) return skeleton_graph(X);
  if (d == 2) {
    std::vector<std::vector<std::size_t>> holders(X.edges.size());
    for (std::size_t f = 0; f < X.faces.size(); f++)
      for (const auto& s : X.faces[f].boundary) holders[s.edge].push_back(f);
    for (const auto& h : holders) join(h);
    return Graph::from_edges(X.faces.size(), edges);
  }
  const auto& top = X.higher[d - 3];
  std::vector<std::vector<std::size_t>> holders(X.cells(d - 1));
  for (std::size_t c = 0; c < top.size(); c++)
    for (const auto& e : top[c]) {
      const long long k = e.value.to_int64();
      for (long long i = 0; i < (k < 0 ? -k : k); i++) holders[e.index].push_back(c);
    }
  for (const auto& h : holders) join(h);
  return Graph::from_edges(top.size(), edges);
}

// A coloured graph as a graph: black b is vertex b, white w is p + w.
inline Graph coloured_graph(const morphe::topology::ColouredGraph& G) {
  std::vector<std::array<std::size_t, 2>> edges;
  for (std::size_t c = 0; c < G.colours; c++)
    for (std::size_t b = 0; b < G.p; b++) edges.push_back({b, G.p + G.sigma[c][b]});
  return Graph::from_edges(2 * G.p, edges);
}

// Starting vertices: all when count >= V (or count == 0), else count distinct vertices drawn
// from a seeded stream.
inline std::vector<std::size_t> sample_vertices(std::size_t V, std::size_t count, std::uint64_t seed) {
  std::vector<std::size_t> all(V);
  std::iota(all.begin(), all.end(), 0);
  if (count == 0 || count >= V) return all;
  calx::random::Stream z(seed, 0x64696d73u, 0, 0);  // "dims"
  for (std::size_t i = 0; i < count; i++) std::swap(all[i], all[i + z.below(V - i)]);
  all.resize(count);
  return all;
}

// P(t), t = 0..steps, the mean return probability of the lazy walk from the starts.
inline std::vector<double> return_probability(const Graph& g, const std::vector<std::size_t>& starts, std::size_t steps) {
  const std::size_t V = g.vertices();
  for (std::size_t v = 0; v < V; v++)
    if (g.degree(v) == 0) throw std::invalid_argument("physis::spacetime::return_probability: an isolated vertex");
  std::vector<double> P(steps + 1, 0.0), p(V), q(V), share(V);
  for (std::size_t s : starts) {
    std::fill(p.begin(), p.end(), 0.0);
    p.at(s) = 1;
    P[0] += 1;
    for (std::size_t t = 1; t <= steps; t++) {
      for (std::size_t v = 0; v < V; v++) {
        q[v] = 0.5 * p[v];
        share[v] = 0.5 * p[v] / static_cast<double>(g.degree(v));
      }
      for (std::size_t v = 0; v < V; v++) {
        if (share[v] == 0) continue;
        for (std::size_t k = g.offset[v]; k < g.offset[v + 1]; k++) q[g.neighbour[k]] += share[v];
      }
      std::swap(p, q);
      P[t] += p[s];
    }
  }
  for (auto& x : P) x /= static_cast<double>(starts.size());
  return P;
}

struct SpectralOptions {
  std::size_t steps = 200;   // the largest diffusion time
  std::size_t starts = 16;   // starting vertices (0: all)
  std::uint64_t seed = 1;
  double window = 0.1;       // relative half-width of the window of the centred differences
};

struct SpectralDimension {
  std::vector<double> t, P;  // t = 0..steps, and P(t)
  std::vector<double> ds;    // d_s(t), NaN where the window does not fit
};

inline SpectralDimension spectral_dimension(const Graph& g, const SpectralOptions& o = {}) {
  SpectralDimension S;
  S.P = return_probability(g, sample_vertices(g.vertices(), o.starts, o.seed), o.steps);
  for (std::size_t t = 0; t <= o.steps; t++) S.t.push_back(static_cast<double>(t));
  S.ds.assign(o.steps + 1, std::nan(""));
  for (std::size_t t = 2; t <= o.steps; t++) {
    const std::size_t w = std::max<std::size_t>(1, static_cast<std::size_t>(o.window * static_cast<double>(t)));
    if (t < w + 1 || t + w > o.steps) continue;
    const double a = static_cast<double>(t - w), b = static_cast<double>(t + w);
    S.ds[t] = -2 * (std::log(S.P[t + w]) - std::log(S.P[t - w])) / (std::log(b) - std::log(a));
  }
  return S;
}

struct HausdorffOptions {
  std::size_t centres = 16;  // centres of the balls (0: all)
  std::uint64_t seed = 1;
  std::size_t radius = 0;    // the largest radius (0: the eccentricity of the centres)
};

struct HausdorffDimension {
  std::vector<double> r, volume;  // r = 0..R, and ⟨|B(r)|⟩
  std::vector<double> dh;         // d_H(r), NaN at the ends
};

inline HausdorffDimension hausdorff_dimension(const Graph& g, const HausdorffOptions& o = {}) {
  const std::size_t V = g.vertices();
  const auto centres = sample_vertices(V, o.centres, o.seed);
  std::vector<double> count;
  std::vector<std::size_t> dist(V);
  for (std::size_t c : centres) {
    std::fill(dist.begin(), dist.end(), SIZE_MAX);
    std::queue<std::size_t> q;
    dist[c] = 0;
    q.push(c);
    std::vector<double> shell;
    while (!q.empty()) {
      const std::size_t v = q.front();
      q.pop();
      if (shell.size() <= dist[v]) shell.resize(dist[v] + 1, 0.0);
      shell[dist[v]] += 1;
      if (o.radius && dist[v] >= o.radius) continue;
      for (std::size_t k = g.offset[v]; k < g.offset[v + 1]; k++)
        if (dist[g.neighbour[k]] == SIZE_MAX) dist[g.neighbour[k]] = dist[v] + 1, q.push(g.neighbour[k]);
    }
    if (count.size() < shell.size()) count.resize(shell.size(), 0.0);
    for (std::size_t r = 0; r < shell.size(); r++) count[r] += shell[r];
  }
  HausdorffDimension H;
  double cum = 0;
  for (std::size_t r = 0; r < count.size(); r++) {
    cum += count[r] / static_cast<double>(centres.size());
    H.r.push_back(static_cast<double>(r));
    H.volume.push_back(cum);
  }
  H.dh.assign(H.r.size(), std::nan(""));
  for (std::size_t r = 2; r + 1 < H.r.size(); r++)
    H.dh[r] = (std::log(H.volume[r + 1]) - std::log(H.volume[r - 1])) / (std::log(static_cast<double>(r + 1)) - std::log(static_cast<double>(r - 1)));
  return H;
}

}  // namespace physis::spacetime
