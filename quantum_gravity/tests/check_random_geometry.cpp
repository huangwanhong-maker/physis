// Checks of physis/quantum_gravity random_geometry.hpp: every labelled tree with up to 4 edges
// gives a planar quadrangulation (a sphere with n + 2 vertices, 2n edges and n faces of
// degree 4, bipartite), whose labels are the distances from the extra vertex; random
// quadrangulations of 2000 faces likewise; the mean distance from the extra vertex grows as
// n^{1/4}; random melonic graphs are connected, of degree 0, and encode spheres; and the
// (D+1)-ary trees with p nodes give as many distinct rooted melonic graphs as the
// Fuss–Catalan numbers count.

#include <physis/quantum_gravity/random_geometry.hpp>
#include <morphe/combinatorics/graph.hpp>
#include <morphe/topology/invariants.hpp>
#include <morphe/topology/surfaces.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace pq = physis::quantum_gravity;
namespace mt = morphe::topology;

namespace {

// Whether a quadrangulation is what the bijection promises; a reason if not.
std::string defect(const pq::Quadrangulation& Q, std::size_t n) {
  const auto& X = Q.complex;
  try {
    X.validate();
  } catch (const std::exception& e) {
    return e.what();
  }
  if (X.vertices != n + 2 || X.edges.size() != 2 * n || X.faces.size() != n) return "counts";
  for (const auto& f : X.faces)
    if (f.boundary.size() != 4) return "a face of degree " + std::to_string(f.boundary.size());
  for (const auto& e : X.edges)
    if (std::abs(Q.label[e[0]] - Q.label[e[1]]) != 1) return "an edge between labels that do not differ by one";
  const auto S = mt::classify_surface(X);
  if (!S.surface || !S.closed || !S.orientable || S.genus != 0) return "not a sphere: " + (S.surface ? S.name : S.reason);
  const auto d = morphe::combinatorics::distances(X.one_skeleton(), Q.extra);
  for (std::size_t v = 0; v < X.vertices; v++)
    if (static_cast<int>(d[v]) != Q.label[v]) return "a label is not the distance from the extra vertex";
  return "";
}

void dyck_words(std::size_t n, const std::function<void(const std::vector<int>&)>& f) {
  std::vector<int> w;
  std::function<void(int, int)> go = [&](int up, int h) {
    if (w.size() == 2 * n) { f(w); return; }
    if (up < static_cast<int>(n)) { w.push_back(1); go(up + 1, h + 1); w.pop_back(); }
    if (h > 0) { w.push_back(-1); go(up, h - 1); w.pop_back(); }
  };
  go(0, 0);
}

// The canonical form of a connected coloured graph rooted at black vertex 0: vertices
// renumbered in breadth-first order through the colours in order.
std::vector<std::vector<std::size_t>> canonical(const mt::ColouredGraph& G) {
  const std::size_t p = G.p;
  std::vector<std::vector<std::size_t>> inv(G.colours, std::vector<std::size_t>(p));
  for (std::size_t c = 0; c < G.colours; c++)
    for (std::size_t b = 0; b < p; b++) inv[c][G.sigma[c][b]] = b;
  std::vector<std::size_t> nb(p, SIZE_MAX), nw(p, SIZE_MAX);
  std::size_t kb = 0, kw = 0;
  std::vector<std::pair<int, std::size_t>> queue{{0, 0}};
  nb[0] = kb++;
  for (std::size_t q = 0; q < queue.size(); q++) {
    const auto [white, v] = queue[q];
    for (std::size_t c = 0; c < G.colours; c++) {
      if (!white) {
        const std::size_t w = G.sigma[c][v];
        if (nw[w] == SIZE_MAX) nw[w] = kw++, queue.push_back({1, w});
      } else {
        const std::size_t b = inv[c][v];
        if (nb[b] == SIZE_MAX) nb[b] = kb++, queue.push_back({0, b});
      }
    }
  }
  std::vector<std::vector<std::size_t>> s(G.colours, std::vector<std::size_t>(p));
  for (std::size_t c = 0; c < G.colours; c++)
    for (std::size_t b = 0; b < p; b++) s[c][nb[b]] = nw[G.sigma[c][b]];
  return s;
}

void lukasiewicz(std::size_t D, std::size_t p, const std::function<void(const std::vector<int>&)>& f) {
  std::vector<int> w;
  // need: open slots still to fill; left: nodes still to place.
  std::function<void(std::size_t, std::size_t)> go = [&](std::size_t need, std::size_t left) {
    if (need == 0) { if (left == 0) f(w); return; }
    if (left > 0) { w.push_back(1); go(need - 1 + D + 1, left - 1); w.pop_back(); }
    if (need > 1 || left == 0) { w.push_back(0); go(need - 1, left); w.pop_back(); }
  };
  w.push_back(1);
  go(D + 1, p - 1);
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- quadrangulations ----------------------------------------------------------------------------------------------
  {
    std::size_t built = 0;
    std::string bad;
    for (std::size_t n = 1; n <= 4; n++)
      dyck_words(n, [&](const std::vector<int>& dyck) {
        std::vector<int> inc(n, -1);
        for (;;) {
          const auto Q = pq::quadrangulation(dyck, inc);
          const std::string d = defect(Q, n);
          if (!d.empty() && bad.empty()) bad = "n = " + std::to_string(n) + ": " + d;
          built++;
          std::size_t k = 0;
          while (k < n && inc[k] == 1) inc[k++] = -1;
          if (k == n) break;
          inc[k]++;
        }
      });
    // Σ_{n ≤ 4} Cat(n) 3^n = 3 + 18 + 135 + 1134 labelled trees.
    s.require("random geometry: the 1290 labelled trees with up to 4 edges give planar quadrangulations (sphere, n + 2 vertices, 2n edges, n faces of degree 4, bipartite) with labels the distances from the extra vertex",
              bad.empty() && built == 1290, bad.empty() ? std::to_string(built) + " trees" : bad);
    std::string rbad;
    for (std::uint64_t seed = 1; seed <= 5; seed++) {
      const auto Q = pq::random_quadrangulation(2000, seed);
      const std::string d = defect(Q, 2000);
      if (!d.empty()) rbad = d;
    }
    s.require("random geometry: five uniform random quadrangulations with 2000 faces are planar, bipartite, with labels the distances", rbad.empty(), rbad);
  }
  {
    // The mean distance from the extra vertex, averaged over samples, against n: it grows as
    // c n^α + d with α = 1/4, the inverse of the Hausdorff dimension 4 (Chassaing and Schaeffer
    // 2004, Theorem 4); the constant d is the offset of the labels (the minimum is 1, not 0)
    // and the first finite-size correction. α is fitted by least squares in (c, d) for each α on
    // a grid. With 150 samples at each of six sizes the statistical error of α is about 0.01;
    // the limit is 0.03.
    std::vector<double> ns, means;
    for (std::size_t n : {250u, 500u, 1000u, 2000u, 4000u, 8000u}) {
      double total = 0;
      const int samples = 150;
      for (int k = 0; k < samples; k++) {
        const auto Q = pq::random_quadrangulation(n, 1000 + 7919 * static_cast<std::uint64_t>(k) + n);
        double m = 0;
        for (std::size_t v = 0; v <= n; v++) m += Q.label[v];
        total += m / static_cast<double>(n + 1);
      }
      ns.push_back(static_cast<double>(n));
      means.push_back(total / samples);
    }
    double best_alpha = 0, best = INFINITY;
    for (double alpha = 0.10; alpha <= 0.40; alpha += 0.0005) {
      // Least squares of means ≈ c x + d, x = n^α.
      double sx = 0, sy = 0, sxx = 0, sxy = 0;
      const double m = static_cast<double>(ns.size());
      for (std::size_t i = 0; i < ns.size(); i++) {
        const double x = std::pow(ns[i], alpha);
        sx += x, sy += means[i], sxx += x * x, sxy += x * means[i];
      }
      const double c = (m * sxy - sx * sy) / (m * sxx - sx * sx), d = (sy - c * sx) / m;
      double res = 0;
      for (std::size_t i = 0; i < ns.size(); i++) res += std::pow(means[i] - c * std::pow(ns[i], alpha) - d, 2);
      if (res < best) best = res, best_alpha = alpha;
    }
    s.check("random geometry: the mean distance from the extra vertex grows as c n^α + d with α = 1/4, over n = 250..8000 (α - 1/4)", best_alpha - 0.25, 0.03);
  }
  // ---- melonic graphs ---------------------------------------------------------------------------------------------------
  {
    std::string bad;
    for (std::size_t D : {2u, 3u, 4u})
      for (std::uint64_t seed = 1; seed <= 4; seed++) {
        const auto G = pq::random_melonic_graph(D, D == 4 ? 12 : 30, seed);
        G.validate();
        if (!G.connected()) bad = "disconnected";
        else if (D >= 3 && G.degree() != 0) bad = "degree " + std::to_string(G.degree());
        else {
          const auto H = mt::homology(G.dual_complex(), false);
          for (std::size_t k = 0; k <= D; k++)
            if (H[k].rank != (k == 0 || k == D ? 1u : 0u) || !H[k].torsion.empty()) bad = "not a sphere in dimension " + std::to_string(D);
          if (D == 2 && mt::classify_surface(G.dual_complex()).genus != 0) bad = "not planar";
        }
      }
    s.require("random geometry: random melonic graphs (D = 2, 3, 4) are connected, of degree 0, and encode spheres", bad.empty(), bad);
  }
  {
    // (D+1)-ary trees with p nodes give distinct rooted melonic graphs, as many as Fuss–Catalan.
    bool ok = true;
    std::string detail;
    for (std::size_t p = 1; p <= 5; p++) {
      std::set<std::vector<std::vector<std::size_t>>> forms;
      std::size_t trees = 0;
      lukasiewicz(3, p, [&](const std::vector<int>& w) {
        const auto G = pq::melonic_graph(3, w);
        trees++;
        if (G.p != p || !G.connected() || G.degree() != 0) ok = false;
        forms.insert(canonical(G));
      });
      // The Fuss–Catalan number C(4p, p) / (3p + 1) (D = 3).
      double c = 1;
      for (std::size_t k = 1; k <= p; k++) c = c * static_cast<double>(3 * p + k) / static_cast<double>(k);
      const auto want = static_cast<std::size_t>(std::llround(c / static_cast<double>(3 * p + 1)));
      if (trees != want || forms.size() != want) ok = false;
      detail += std::to_string(forms.size()) + (p < 5 ? ", " : "");
    }
    s.require("random geometry: 4-ary trees with p = 1..5 nodes give distinct rooted melonic 4-coloured graphs, 1, 4, 22, 140, 969 (Fuss–Catalan)", ok, detail);
  }
  return s.done();
}
