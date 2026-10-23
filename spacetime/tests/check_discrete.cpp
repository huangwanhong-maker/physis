// Checks of physis/spacetime's discrete spacetimes.
//
// Regge geometry: volumes of regular simplices, and of random ones by Gram and by Cayley-
// Menger; the dihedral angles of the regular tetrahedron and 4-simplex, and the sine law of
// dihedral angles on a random simplex; the derivative of the volume against central
// differences; Schläfli's identity Σ_h A_h dθ_h = 0 on random 3- and 4-simplices; Gauss-
// Bonnet Σ_v ε_v = 2πχ on the octahedron, the 6-vertex projective plane and the 7-vertex
// torus with generic edge lengths; a flat triangulated cube, with zero deficits inside and
// the boundary term of Hartle and Sorkin 6πn on its twelve edges; the boundary of the
// 4-simplex (a 3-sphere), with equal deficits 2π - 3 arccos(1/3) and the Regge action's
// gradient by Schläfli against central differences.
//
// Causal sets: the Poisson number of sprinkled elements; the ordering fraction in two and
// four dimensions against f(d) = Γ(d+1)Γ(d/2)/(2Γ(3d/2)), and the Myrheim-Meyer dimension;
// links and interval abundances against the order (the transitive closure of the links is
// the order; Σ_k N_k = R); the longest chain in two dimensions against the mean longest
// increasing subsequence 2√N - 1.7711 N^{1/6}.

#include <physis/spacetime/causal_set.hpp>
#include <physis/spacetime/regge.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <numbers>
#include <vector>

namespace ps = physis::spacetime;

namespace {

const double pi = std::numbers::pi;

// Squared lengths of a k-simplex with every edge of squared length 1 + noise.
ps::SquaredLengths random_simplex(std::size_t k, calx::random::Stream& z, double noise = 0.2) {
  ps::SquaredLengths s(k + 1, std::vector<double>(k + 1, 0.0));
  for (std::size_t i = 0; i <= k; i++)
    for (std::size_t j = i + 1; j <= k; j++) s[i][j] = s[j][i] = 1 + noise * (2 * z.uniform() - 1);
  return s;
}

// Squared lengths for every edge of a complex, each 1 + noise.
ps::EdgeLengths random_lengths(const ps::ReggeComplex& c, calx::random::Stream& z, double noise) {
  ps::EdgeLengths L;
  for (const auto& [a, b] : c.edges()) L.set(a, b, 1 + noise * (2 * z.uniform() - 1));
  return L;
}

// Schläfli: Σ_h A_h ∂θ_h/∂s_e over the hinges of one simplex, for each edge e, relative to Σ_h A_h |∂θ_h/∂s_e|.
double schlafli(const ps::SquaredLengths& s) {
  const std::size_t k = s.size() - 1;
  double worst = 0;
  for (std::size_t a = 0; a <= k; a++)
    for (std::size_t b = a + 1; b <= k; b++) {
      double sum = 0, size = 0;
      const double h = 1e-6;
      for (std::size_t i = 0; i <= k; i++)
        for (std::size_t j = i + 1; j <= k; j++) {
          ps::SquaredLengths hinge;
          std::vector<std::size_t> v;
          for (std::size_t m = 0; m <= k; m++)
            if (m != i && m != j) v.push_back(m);
          hinge.assign(v.size(), std::vector<double>(v.size(), 0.0));
          for (std::size_t p = 0; p < v.size(); p++)
            for (std::size_t q = 0; q < v.size(); q++) hinge[p][q] = s[v[p]][v[q]];
          auto plus = s, minus = s;
          plus[a][b] = plus[b][a] = s[a][b] + h;
          minus[a][b] = minus[b][a] = s[a][b] - h;
          const double dtheta = (ps::dihedral_angle(plus, i, j) - ps::dihedral_angle(minus, i, j)) / (2 * h);
          const double A = ps::simplex_volume(hinge);
          sum += A * dtheta;
          size += A * std::fabs(dtheta);
        }
      worst = std::max(worst, std::fabs(sum) / size);
    }
  return worst;
}

// The Euler characteristic V - E + F of a list of triangles.
long euler(const std::vector<ps::Vertices>& tri) {
  std::map<std::pair<std::size_t, std::size_t>, int> edges;
  std::map<std::size_t, int> vertices;
  for (auto t : tri) {
    std::sort(t.begin(), t.end());
    for (auto v : t) vertices[v] = 1;
    edges[{t[0], t[1]}] = edges[{t[0], t[2]}] = edges[{t[1], t[2]}] = 1;
  }
  return static_cast<long>(vertices.size()) - static_cast<long>(edges.size()) + static_cast<long>(tri.size());
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream z(20261023, 7, 0, 0);

  // ---- simplices ---------------------------------------------------------------------------------------------------------
  {
    double dev = 0;
    for (std::size_t k = 1; k <= 5; k++) {
      ps::SquaredLengths unit(k + 1, std::vector<double>(k + 1, 1.0));
      for (std::size_t i = 0; i <= k; i++) unit[i][i] = 0;
      double f = 1;
      for (std::size_t i = 2; i <= k; i++) f *= static_cast<double>(i);
      const double want = std::sqrt(static_cast<double>(k + 1)) / (f * std::pow(2.0, k / 2.0));
      dev = std::max(dev, std::fabs(ps::simplex_volume(unit) / want - 1));
    }
    s.check("regge: volume of the regular k-simplex sqrt(k+1) / (k! 2^{k/2}), k = 1..5 (relative)", dev, 1e-14);
    double cm = 0;
    for (int t = 0; t < 50; t++) {
      const auto q = random_simplex(2 + t % 3, z);
      cm = std::max(cm, std::fabs(ps::cayley_menger_volume(q) / ps::simplex_volume(q) - 1));
    }
    s.check("regge: Cayley-Menger and Gram volumes agree on 50 random simplices (relative)", cm, 1e-12);
    ps::SquaredLengths tet(4, std::vector<double>(4, 1.0)), pent(5, std::vector<double>(5, 1.0));
    for (int i = 0; i < 4; i++) tet[i][i] = 0;
    for (int i = 0; i < 5; i++) pent[i][i] = 0;
    s.check("regge: dihedral angle of the regular tetrahedron arccos(1/3)", ps::dihedral_angle(tet, 0, 3) - std::acos(1.0 / 3), 1e-15);
    s.check("regge: of the regular 4-simplex arccos(1/4)", ps::dihedral_angle(pent, 1, 2) - std::acos(0.25), 1e-15);
    // The sine law: sin θ_ij = k V V_ij / ((k - 1) V_i V_j), V_i the facet opposite i, V_ij the hinge.
    double sine = 0;
    for (int t = 0; t < 30; t++) {
      const std::size_t k = 3 + t % 2;
      const auto q = random_simplex(k, z, 0.3);
      auto sub = [&](std::vector<std::size_t> drop) {
        std::vector<std::size_t> v;
        for (std::size_t m = 0; m <= k; m++)
          if (std::find(drop.begin(), drop.end(), m) == drop.end()) v.push_back(m);
        ps::SquaredLengths r(v.size(), std::vector<double>(v.size()));
        for (std::size_t p = 0; p < v.size(); p++)
          for (std::size_t w = 0; w < v.size(); w++) r[p][w] = q[v[p]][v[w]];
        return ps::simplex_volume(r);
      };
      const double want = k * ps::simplex_volume(q) * sub({0, 1}) / ((k - 1.0) * sub({0}) * sub({1}));
      sine = std::max(sine, std::fabs(std::sin(ps::dihedral_angle(q, 0, 1)) - want));
    }
    s.check("regge: the sine law of dihedral angles on 30 random 3- and 4-simplices", sine, 1e-13);
    double grad = 0;
    for (int t = 0; t < 10; t++) {
      const auto q = random_simplex(3 + t % 2, z);
      const auto g = ps::volume_gradient(q);
      for (std::size_t a = 0; a < q.size(); a++)
        for (std::size_t b = a + 1; b < q.size(); b++) {
          auto p = q, m = q;
          const double h = 1e-6;
          p[a][b] = p[b][a] = q[a][b] + h;
          m[a][b] = m[b][a] = q[a][b] - h;
          grad = std::max(grad, std::fabs((ps::simplex_volume(p) - ps::simplex_volume(m)) / (2 * h) - g[a][b]) / ps::simplex_volume(q));
        }
    }
    s.check("regge: ∂V/∂s_ab in closed form against central differences (relative to V)", grad, 1e-8);
    double sch = 0;
    for (int t = 0; t < 10; t++) sch = std::max(sch, schlafli(random_simplex(3 + t % 2, z)));
    s.check("regge: Schläfli's identity Σ_h A_h dθ_h = 0 on random 3- and 4-simplices (relative)", sch, 1e-8);
  }
  // ---- Gauss-Bonnet --------------------------------------------------------------------------------------------------------
  {
    const std::vector<ps::Vertices> octahedron{{0, 2, 4}, {0, 4, 3}, {0, 3, 5}, {0, 5, 2}, {1, 2, 4}, {1, 4, 3}, {1, 3, 5}, {1, 5, 2}};
    const std::vector<ps::Vertices> rp2{{0, 1, 3}, {0, 1, 5}, {0, 2, 4}, {0, 2, 5}, {0, 3, 4}, {1, 2, 3}, {1, 2, 4}, {1, 4, 5}, {2, 3, 5}, {3, 4, 5}};
    std::vector<ps::Vertices> torus_build;
    for (std::size_t i = 0; i < 7; i++) {
      torus_build.push_back({i, (i + 1) % 7, (i + 3) % 7});
      torus_build.push_back({i, (i + 2) % 7, (i + 3) % 7});
    }
    const std::vector<ps::Vertices>& torus = torus_build;
    double worst = 0;
    bool interior = true;
    for (const std::vector<ps::Vertices>* tri : {&octahedron, &rp2, &torus}) {
      const ps::ReggeComplex c(*tri);
      for (std::size_t h = 0; h < c.hinges().size(); h++) interior = interior && c.interior(h);
      for (int trial = 0; trial < 5; trial++) {
        const auto L = random_lengths(c, z, 0.25);
        const auto e = c.deficits(L);
        double sum = 0;
        for (double x : e) sum += x;
        worst = std::max(worst, std::fabs(sum - 2 * pi * static_cast<double>(euler(*tri))));
      }
    }
    s.check("regge: Gauss-Bonnet Σ ε_v = 2πχ on the octahedron, RP² (6 vertices) and the torus (7 vertices), random lengths", worst, 1e-12);
    s.require("regge: every vertex of a closed surface is an interior hinge", interior);
  }
  // ---- a flat triangulated cube ------------------------------------------------------------------------------------------
  {
    // The cube [0, n]³ cut into n³ unit cubes, each into six tetrahedra along its main diagonal
    // (Freudenthal's triangulation).
    const std::size_t n = 2, m = n + 1;
    auto id = [&](std::size_t x, std::size_t y, std::size_t w) { return x + m * (y + m * w); };
    std::vector<ps::Vertices> tets;
    const int perms[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    for (std::size_t x = 0; x < n; x++)
      for (std::size_t y = 0; y < n; y++)
        for (std::size_t w = 0; w < n; w++)
          for (const auto& p : perms) {
            std::array<std::size_t, 3> c{x, y, w};
            ps::Vertices t{id(c[0], c[1], c[2])};
            for (int k = 0; k < 3; k++) {
              c[p[k]]++;
              t.push_back(id(c[0], c[1], c[2]));
            }
            tets.push_back(t);
          }
    auto coords = [&](std::size_t v) { return std::array<double, 3>{double(v % m), double((v / m) % m), double(v / (m * m))}; };
    auto len = [&](std::size_t a, std::size_t b) {
      const auto p = coords(a), q = coords(b);
      return (p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2]);
    };
    const ps::ReggeComplex cube(tets);
    const auto e = cube.deficits(len);
    double inside = 0;
    std::size_t interior = 0;
    for (std::size_t h = 0; h < e.size(); h++)
      if (cube.interior(h)) { inside = std::max(inside, std::fabs(e[h])); interior++; }
    s.check("regge: a flat triangulated cube (48 tetrahedra) has zero deficit at every interior hinge", inside, 1e-14);
    s.require("regge: interior hinges found", interior > 0, std::to_string(interior) + " of " + std::to_string(e.size()));
    const auto A = cube.hinge_volumes(len);
    double scale = 0;
    for (std::size_t h = 0; h < A.size(); h++)
      if (cube.interior(h)) scale += 2 * pi * A[h];
    s.check("regge: its Regge action vanishes (|S| / Σ_interior 2π A_h)", cube.action(len) / scale, 1e-14);
    s.check("regge: with the boundary term, the cube's twelve edges each give (π/2) n (relative to 6πn)", cube.action(len, true) / (6 * pi * n) - 1, 1e-14);
  }
  // ---- the boundary of the 4-simplex ---------------------------------------------------------------------------------------
  {
    std::vector<ps::Vertices> sphere;
    for (std::size_t drop = 0; drop < 5; drop++) {
      ps::Vertices t;
      for (std::size_t v = 0; v < 5; v++)
        if (v != drop) t.push_back(v);
      sphere.push_back(t);
    }
    const ps::ReggeComplex c(sphere);
    auto unit = [](std::size_t, std::size_t) { return 1.0; };
    const auto e = c.deficits(unit);
    double dev = 0;
    for (double x : e) dev = std::max(dev, std::fabs(x - (2 * pi - 3 * std::acos(1.0 / 3))));
    s.check("regge: the boundary of the regular 4-simplex: every edge has deficit 2π - 3 arccos(1/3)", dev, 1e-14);
    const auto L = random_lengths(c, z, 0.15);
    const auto g = c.action_gradient(L);
    double worst = 0;
    for (const auto& [edge, value] : g) {
      const double h = 1e-6;
      ps::EdgeLengths p = L, q = L;
      p.set(edge.first, edge.second, L(edge.first, edge.second) + h);
      q.set(edge.first, edge.second, L(edge.first, edge.second) - h);
      worst = std::max(worst, std::fabs((c.action(p) - c.action(q)) / (2 * h) - value) / std::max(1.0, std::fabs(value)));
    }
    s.check("regge: on a 3-sphere with random lengths, ∂S/∂s_e by Schläfli against central differences of S", worst, 1e-8);
  }
  // ---- causal sets ---------------------------------------------------------------------------------------------------------
  {
    // The number of sprinkled elements: mean ρV over 40 sprinklings.
    const double rho = 500, tau = 2;
    const double V = ps::alexandrov_volume(2, tau);
    double sum = 0;
    const int m = 40;
    for (int k = 0; k < m; k++) {
      calx::random::Stream w(11, 1, static_cast<std::uint32_t>(k), 0);
      sum += static_cast<double>(ps::sprinkle_diamond<2>(tau, rho, w).size());
    }
    const double mean = sum / m, sd = std::sqrt(rho * V / m);
    s.check("causal sets: the number sprinkled is Poisson with mean ρV (deviation of the mean of 40, in standard errors)", (mean - rho * V) / sd, 4);
  }
  {
    s.check("causal sets: f(1) = 1, f(2) = 1/2, f(4) = 1/10", std::fabs(ps::expected_ordering_fraction(1) - 1) + std::fabs(ps::expected_ordering_fraction(2) - 0.5) + std::fabs(ps::expected_ordering_fraction(4) - 0.1), 1e-14);
    s.check("causal sets: the Myrheim-Meyer dimension inverts f", ps::myrheim_meyer_dimension(ps::expected_ordering_fraction(3.3)) - 3.3, 1e-10);
    // The ordering fraction of sprinklings into 2- and 4-dimensional intervals, 20 each of about
    // 1000 elements, against f(d), within four standard errors of the sample.
    auto ensemble = [&](auto dim) {
      constexpr std::size_t d = decltype(dim)::value;
      const double tau = 1.0, rho = 1000 / ps::alexandrov_volume(d, tau);
      std::vector<double> r;
      for (int k = 0; k < 20; k++) {
        calx::random::Stream w(12, static_cast<std::uint32_t>(d), static_cast<std::uint32_t>(k), 0);
        r.push_back(ps::CausalSet<d>(ps::sprinkle_diamond<d>(tau, rho, w)).ordering_fraction());
      }
      double mean = 0, var = 0;
      for (double x : r) mean += x / r.size();
      for (double x : r) var += (x - mean) * (x - mean) / (r.size() - 1);
      return std::pair{mean, std::sqrt(var / r.size())};
    };
    const auto [r2, e2] = ensemble(std::integral_constant<std::size_t, 2>{});
    const auto [r4, e4] = ensemble(std::integral_constant<std::size_t, 4>{});
    s.check("causal sets: ordering fraction in 2 dimensions against 1/2 (in standard errors of the mean of 20)", (r2 - 0.5) / e2, 4);
    s.check("causal sets: ordering fraction in 4 dimensions against 1/10 (in standard errors)", (r4 - 0.1) / e4, 4);
    s.check("causal sets: the Myrheim-Meyer dimension of the 4-dimensional sprinklings", ps::myrheim_meyer_dimension(r4) - 4, 0.1);
  }
  {
    // Links, intervals and the order on one small sprinkling.
    calx::random::Stream w(13, 0, 0, 0);
    const ps::CausalSet<3> c(ps::sprinkle_diamond<3>(1.0, 300 / ps::alexandrov_volume(3, 1.0), w));
    const auto links = c.links();
    const std::size_t n = c.size();
    std::vector<std::vector<bool>> reach(n, std::vector<bool>(n, false));
    for (const auto& [a, b] : links) reach[a][b] = true;
    for (std::size_t j = 0; j < n; j++)  // elements are in time order, so one pass in order closes the relation
      for (std::size_t i = j; i-- > 0;)
        if (reach[i][j])
          for (std::size_t k = 0; k < i; k++)
            if (reach[k][i]) reach[k][j] = true;
    bool same = true;
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++) same = same && (reach[i][j] == c.precedes(i, j));
    s.require("causal sets: the transitive closure of the links is the order (300 elements in 3 dimensions)", same, std::to_string(links.size()) + " links, " + std::to_string(c.relations()) + " relations");
    const auto N = c.interval_abundances(n);
    std::size_t total = 0;
    for (auto x : N) total += x;
    s.require("causal sets: N_0 is the number of links and Σ_k N_k = R", N[0] == links.size() && total == c.relations());
  }
  {
    // The longest chain in two dimensions, against 2√N - 1.7711 N^{1/6} (Baik, Deift and Johansson):
    // the mean over 12 sprinklings of (L - 2√N + 1.7711 N^{1/6}) / N^{1/6}, whose fluctuation is the
    // Tracy-Widom (β = 2) law of standard deviation 0.9018; allowed four standard errors and 0.3
    // for the corrections of lower order at N ≈ 3000.
    double z_sum = 0;
    const int m = 12;
    for (int k = 0; k < m; k++) {
      calx::random::Stream w(14, static_cast<std::uint32_t>(k), 0, 0);
      const ps::CausalSet<2> c(ps::sprinkle_diamond<2>(1.0, 6000, w));
      const double N = static_cast<double>(c.size()), n6 = std::pow(N, 1.0 / 6);
      z_sum += (static_cast<double>(c.longest_chain()) - 2 * std::sqrt(N) + 1.7711 * n6) / n6;
    }
    s.check("causal sets: the longest chain in 2 dimensions, mean of (L - 2√N + 1.7711 N^{1/6}) / N^{1/6}", z_sum / m, 4 * 0.9018 / std::sqrt(m) + 0.3);
  }
  return s.done();
}
