// Checks of physis/spacetime/dimensions.hpp: the return probability of the lazy walk on periodic
// lattices in one to four dimensions against the exact Fourier sum; the spectral dimension near d
// at intermediate scales on large tori and 1 on a long cycle; ball volumes on lattices against
// the exact count of lattice points within an L¹ distance, and the Hausdorff dimension they give;
// the square and cubic lattices self-dual (the dual graphs of their top cells give the same walk);
// a coloured graph's walk; and uniform random trees, whose dimensions approach those of branched
// polymers, d_s = 4/3 and d_H = 2.

#include <physis/spacetime/dimensions.hpp>
#include <morphe/topology/constructions.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <string>
#include <vector>

namespace sp = physis::spacetime;
namespace mt = morphe::topology;

namespace {

// P(t) of the lazy walk on the periodic lattice (ℤ/L)^d: (1/L^d) Σ_k λ_k^t,
// λ_k = 1/2 + (1/2d) Σ_i cos(2π k_i / L).
std::vector<double> fourier(std::size_t d, std::size_t L, std::size_t steps) {
  std::vector<double> P(steps + 1, 0.0);
  std::vector<std::size_t> k(d, 0);
  std::size_t total = 1;
  for (std::size_t i = 0; i < d; i++) total *= L;
  for (std::size_t idx = 0; idx < total; idx++) {
    std::size_t r = idx;
    double lambda = 0.5;
    for (std::size_t i = 0; i < d; i++) {
      lambda += std::cos(2 * std::numbers::pi * static_cast<double>(r % L) / static_cast<double>(L)) / (2.0 * static_cast<double>(d));
      r /= L;
    }
    double pw = 1;
    for (std::size_t t = 0; t <= steps; t++) P[t] += pw, pw *= lambda;
  }
  for (auto& x : P) x /= static_cast<double>(total);
  return P;
}

// The number of points of ℤ^d within L¹ distance r: Σ_k 2^k C(d, k) C(r, k).
double l1_ball(std::size_t d, std::size_t r) {
  auto C = [](double n, std::size_t k) {
    double c = 1;
    for (std::size_t i = 0; i < k; i++) c = c * (n - static_cast<double>(i)) / static_cast<double>(i + 1);
    return c;
  };
  double s = 0;
  for (std::size_t k = 0; k <= d; k++) s += std::pow(2.0, static_cast<double>(k)) * C(static_cast<double>(d), k) * C(static_cast<double>(r), k);
  return s;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- the return probability, exactly -------------------------------------------------------------------------------
  {
    const std::size_t sizes[5] = {0, 50, 16, 10, 6};
    for (std::size_t d = 1; d <= 4; d++) {
      const auto X = mt::torus_lattice(std::vector<std::size_t>(d, sizes[d]));
      const auto g = sp::skeleton_graph(X);
      const auto P = sp::return_probability(g, {0}, 60);
      const auto F = fourier(d, sizes[d], 60);
      double dev = 0;
      for (std::size_t t = 0; t <= 60; t++) dev = std::max(dev, std::fabs(P[t] / F[t] - 1));
      s.check("dimensions: lazy-walk return probability on the " + std::to_string(sizes[d]) + "^" + std::to_string(d) + " torus against the Fourier sum (relative, t <= 60)", dev, 1e-12);
    }
  }
  // ---- the spectral dimension near d ---------------------------------------------------------------------------------
  {
    // For 1 ≪ t ≪ L², P(t) = (d / (π t))^{d/2} (1 + O(1/t)) on ℤ^d (the lazy walk's variance is
    // 1/(2d) per coordinate and step), so d_s(t) = d + O(1/t). The checks take t well inside the
    // window and allow c d / t with c = 2, the size of the leading correction found by the
    // Fourier sums themselves.
    struct Case { std::size_t d, L, steps, t; };
    for (const Case c : {Case{1, 4000, 400, 200}, Case{2, 100, 300, 150}, Case{3, 30, 100, 50}, Case{4, 14, 40, 20}}) {
      const auto g = sp::skeleton_graph(mt::torus_lattice(std::vector<std::size_t>(c.d, c.L)));
      const auto S = sp::spectral_dimension(g, {c.steps, 1, 1, 0.1});
      const double dev = S.ds[c.t] - static_cast<double>(c.d);
      s.check("dimensions: spectral dimension of the " + std::to_string(c.L) + "^" + std::to_string(c.d) + " torus at t = " + std::to_string(c.t) + ", minus " +
                  std::to_string(c.d),
              dev, 2.0 * static_cast<double>(c.d) / static_cast<double>(c.t));
    }
    // A long cycle: one dimension at every intermediate scale.
    const auto ring = sp::skeleton_graph(mt::circle(10000));
    const auto S = sp::spectral_dimension(ring, {400, 1, 1, 0.1});
    double worst = 0;
    for (std::size_t t = 40; t <= 360; t += 20) worst = std::max(worst, std::fabs(S.ds[t] - 1));
    s.check("dimensions: a cycle of 10000 vertices has d_s = 1 for 40 <= t <= 360 (deviation, ~2/t)", worst, 0.05);
  }
  // ---- Hausdorff dimension -----------------------------------------------------------------------------------------------
  {
    int bad = 0;
    std::string detail;
    for (std::size_t d = 1; d <= 3; d++) {
      const std::size_t L = d == 1 ? 201 : d == 2 ? 61 : 21, R = (L - 1) / 2;
      const auto H = sp::hausdorff_dimension(sp::skeleton_graph(mt::torus_lattice(std::vector<std::size_t>(d, L))), {4, 3, R});
      for (std::size_t r = 0; r <= R; r++)
        if (H.volume[r] != l1_ball(d, r)) bad++, detail += std::to_string(d) + "/" + std::to_string(r) + " ";
    }
    s.require("dimensions: ball volumes on lattice tori equal the counts of lattice points within an L¹ distance (d = 1, 2, 3)", bad == 0, detail);
    // d_H from the exact counts approaches d as the lower terms of the polynomial fade: at r = 30
    // in two dimensions, 2r(2r + 1) / (2r² + 2r + 1) = 1.967 for the derivative itself.
    const auto H2 = sp::hausdorff_dimension(sp::skeleton_graph(mt::torus_lattice({101, 101})), {2, 5, 0});
    s.check("dimensions: Hausdorff dimension of the 101² torus at r = 30, minus 2 (lower terms ~ 1/r)", H2.dh[30] - 2, 0.06);
  }
  // ---- duals and coloured graphs -------------------------------------------------------------------------------------------
  {
    // The square and cubic lattices are self-dual: the walk on the dual graph of their top cells
    // is the walk on their vertices.
    const auto X = mt::torus_lattice({12, 12});
    const auto D = sp::dual_graph(X);
    bool deg = D.vertices() == X.faces.size();
    for (std::size_t v = 0; v < D.vertices(); v++) deg = deg && D.degree(v) == 4;
    const auto P = sp::return_probability(sp::skeleton_graph(X), {0}, 50), Q = sp::return_probability(D, {0}, 50);
    double dev = 0;
    for (std::size_t t = 0; t <= 50; t++) dev = std::max(dev, std::fabs(P[t] - Q[t]));
    s.require("dimensions: the dual graph of the 12² torus's plaquettes is 4-regular on its 144 faces", deg);
    s.check("dimensions: and its walk is the walk on the vertices (max |ΔP|)", dev, 1e-14);
    const auto Y = mt::torus_lattice({6, 6, 6});
    const auto D3 = sp::dual_graph(Y);
    bool deg3 = D3.vertices() == 216;
    for (std::size_t v = 0; v < D3.vertices(); v++) deg3 = deg3 && D3.degree(v) == 6;
    const auto P3 = sp::return_probability(sp::skeleton_graph(Y), {0}, 30), Q3 = sp::return_probability(D3, {0}, 30);
    double dev3 = 0;
    for (std::size_t t = 0; t <= 30; t++) dev3 = std::max(dev3, std::fabs(P3[t] - Q3[t]));
    s.check("dimensions: the cubes of the 6³ torus, joined across squares, walk as its vertices do (max |ΔP|; 6-regular)", deg3 ? dev3 : 1.0, 1e-14);
    // The elementary melon: two vertices joined by D + 1 edges; P(t) = 1/2 for t >= 1.
    const auto M = sp::coloured_graph(mt::elementary_melon(3));
    const auto PM = sp::return_probability(M, {0}, 10);
    double devm = 0;
    for (std::size_t t = 1; t <= 10; t++) devm = std::max(devm, std::fabs(PM[t] - 0.5));
    s.check("dimensions: the walk on the elementary 4-coloured melon returns with probability 1/2", devm, 1e-15);
  }
  // ---- branched polymers -----------------------------------------------------------------------------------------------------
  {
    // Uniform random labelled trees (Prüfer 1918), the critical Poisson Galton–Watson trees
    // conditioned on their size, the branched polymers of the melonic phase. Their dimensions
    // tend to d_s = 4/3 and d_H = 2 (Jonsson and Wheater 1998; Durhuus, Jonsson and Wheater
    // 2007); at 20000 vertices the scales are short and the estimates biased, so the checks
    // ask only that they lie between those of a line and of a plane, near the limits.
    const std::size_t n = 20000;
    calx::random::Stream z(11, 0x74726565u, 0, 0);  // "tree"
    std::vector<std::size_t> prufer(n - 2), deg(n, 1);
    for (auto& x : prufer) x = z.below(n), deg[x]++;
    std::vector<std::array<std::size_t, 2>> edges;
    std::size_t leaf = 0;
    while (deg[leaf] != 1) leaf++;
    std::size_t ptr = leaf;
    for (std::size_t x : prufer) {
      edges.push_back({leaf, x});
      if (--deg[x] == 1 && x < ptr) leaf = x;
      else {
        ptr++;
        while (deg[ptr] != 1) ptr++;
        leaf = ptr;
      }
    }
    edges.push_back({leaf, n - 1});
    const auto T = sp::Graph::from_edges(n, edges);
    const auto S = sp::spectral_dimension(T, {400, 24, 5, 0.2});
    const auto H = sp::hausdorff_dimension(T, {24, 5, 0});
    s.require("dimensions: the decoded Prüfer sequence is a tree: n − 1 edges, every ball reaching all n vertices",
              edges.size() == n - 1 && std::fabs(H.volume.back() - static_cast<double>(n)) < 1e-9 * static_cast<double>(n));
    double ds = 0, dh = 0;
    for (std::size_t t = 100; t <= 300; t += 50) ds += S.ds[t] / 5;
    for (std::size_t r = 8; r <= 16; r += 2) dh += H.dh[r] / 5;
    s.check("dimensions: a uniform random tree of 20000 vertices: d_s over 100 <= t <= 300, minus 4/3 (finite size)", ds - 4.0 / 3, 0.2);
    s.check("dimensions: and d_H over 8 <= r <= 16, minus 2 (finite size)", dh - 2, 0.4);
  }
  return s.done();
}
