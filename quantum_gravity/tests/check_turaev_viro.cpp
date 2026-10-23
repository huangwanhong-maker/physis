// Checks of physis/quantum_gravity's Turaev–Viro state sum: independence of the triangulation
// (S³ as ∂Δ⁴, as the suspension of ∂Δ³, as the join of two triangles, and as the elementary
// melon of the Boulatov model), the values TV(S³) = w^{−2} and TV(S² × S¹) = 1 (on a coloured
// graph and on a simplicial product), and the lens spaces L(p, 1), built as quotients of joins
// of cycles by a free rotation, against |τ_r(L(p, 1))|² by the theorem of Turaev and Walker.

#include <physis/quantum_gravity/tensor_models.hpp>
#include <physis/quantum_gravity/turaev_viro.hpp>
#include <calx/verification/check.hpp>
#include <morphe/topology/constructions.hpp>
#include <morphe/topology/invariants.hpp>
#include <morphe/topology/manifolds.hpp>

#include <algorithm>
#include <chrono>
#include <map>
#include <optional>
#include <string>

namespace qg = physis::quantum_gravity;
namespace mt = morphe::topology;

namespace {

// The join of the cycles C_m and C_n (a triangulated S³), vertices a_i = i, b_j = m + j.
mt::SimplicialComplex join_of_cycles(std::size_t m, std::size_t n) {
  std::vector<mt::Simplex> f;
  for (std::size_t i = 0; i < m; i++)
    for (std::size_t j = 0; j < n; j++) f.push_back({i, (i + 1) % m, m + j, m + (j + 1) % n});
  return mt::SimplicialComplex::from_facets(f, m + n, "C" + std::to_string(m) + " * C" + std::to_string(n));
}

// The lens space L(p, q) as the quotient of C_{ps} * C_{pt} by the free rotation
// a_i ↦ a_{i+s}, b_j ↦ b_{j+qt}, the restriction of (z, w) ↦ (ζ z, ζ^q w) on S³ ⊂ ℂ².
qg::Triangulation3 lens(std::size_t p, std::size_t q, std::size_t s, std::size_t t) {
  const std::size_t m = p * s, n = p * t;
  const auto K = join_of_cycles(m, n);
  auto g = [&](std::size_t v, std::size_t k) { return v < m ? (v + k * s) % m : m + (v - m + k * q * t) % n; };
  auto canonical = [&](mt::Simplex x) {
    mt::Simplex best;
    for (std::size_t k = 0; k < p; k++) {
      mt::Simplex y;
      for (auto v : x) y.push_back(g(v, k));
      std::sort(y.begin(), y.end());
      if (best.empty() || y < best) best = y;
    }
    return best;
  };
  std::vector<std::map<mt::Simplex, std::size_t>> orbit(4);
  for (std::size_t d = 0; d <= 3; d++)
    for (const auto& x : K.simplices[d]) orbit[d].emplace(canonical(x), orbit[d].size());
  qg::Triangulation3 T;
  T.vertices = orbit[0].size();
  T.edges = orbit[1].size();
  auto edge = [&](std::size_t a, std::size_t b) { return orbit[1].at(canonical({std::min(a, b), std::max(a, b)})); };
  T.triangles.resize(orbit[2].size());
  for (const auto& [x, i] : orbit[2]) T.triangles[i] = {edge(x[0], x[1]), edge(x[0], x[2]), edge(x[1], x[2])};
  T.tetrahedra.resize(orbit[3].size());
  for (const auto& [x, i] : orbit[3]) T.tetrahedra[i] = {edge(x[0], x[1]), edge(x[0], x[2]), edge(x[0], x[3]), edge(x[1], x[2]), edge(x[1], x[3]), edge(x[2], x[3])};
  return T;
}

// The product ∂Δ³ × C_n, each prism [a b c] × [i, i+1] cut into three tetrahedra by the
// staircase (a_i b_i c_i c_{i+1}), (a_i b_i b_{i+1} c_{i+1}), (a_i a_{i+1} b_{i+1} c_{i+1}).
mt::SimplicialComplex sphere_times_circle(std::size_t n) {
  const auto S = mt::simplex_boundary(3);
  auto v = [&](std::size_t x, std::size_t i) { return x * n + i % n; };
  std::vector<mt::Simplex> f;
  for (const auto& t : S.simplices[2])
    for (std::size_t i = 0; i < n; i++) {
      const std::size_t a = t[0], b = t[1], c = t[2];
      f.push_back({v(a, i), v(b, i), v(c, i), v(c, i + 1)});
      f.push_back({v(a, i), v(b, i), v(b, i + 1), v(c, i + 1)});
      f.push_back({v(a, i), v(a, i + 1), v(b, i + 1), v(c, i + 1)});
    }
  return mt::SimplicialComplex::from_facets(f, 4 * n, "S2 x S1");
}

// Whether every 3-bubble of a 4-coloured graph is a 2-sphere, so that the graph encodes a
// closed 3-manifold: with one colour removed, the components have Σ (2 − 2g) = V − E + F.
bool bubbles_are_spheres(const mt::ColouredGraph& G) {
  for (std::size_t c = 0; c < 4; c++) {
    std::size_t components;
    G.components(15u & ~(1u << c), &components);
    long long F = 0;
    for (std::size_t i = 0; i < 4; i++)
      for (std::size_t j = i + 1; j < 4; j++)
        if (i != c && j != c) F += static_cast<long long>(G.faces(i, j).size());
    const long long V = 2 * static_cast<long long>(G.p), E = 3 * static_cast<long long>(G.p);
    if (V - E + F != 2 * static_cast<long long>(components)) return false;
  }
  return true;
}

std::string homology_text(const mt::CWComplex& X) {
  std::string s;
  const auto H = mt::homology(X, false);
  for (std::size_t k = 0; k < H.size(); k++) s += (k ? " ; " : "") + H[k].str();
  return s;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  const auto start = std::chrono::steady_clock::now();

  // ---- S³: four triangulations ------------------------------------------------------------------------------------------
  {
    const auto D4 = mt::simplex_boundary(4);
    const auto Su = mt::suspension(mt::simplex_boundary(3));
    const auto J = join_of_cycles(3, 3);
    const auto M = mt::elementary_melon(3);
    const auto MJ = mt::manifold_type(J);
    s.require("quantum_gravity: the join C3 * C3 is a closed combinatorial 3-manifold with the homology of S³",
              MJ.manifold && MJ.closed && homology_text(J.cw()) == "ℤ ; 0 ; 0 ; ℤ");
    for (int r = 3; r <= 6; r++) {
      const double want = 1 / physis::symmetry::total_quantum_dimension_squared(r);
      double worst = 0;
      for (const auto& T : {qg::triangulation(D4), qg::triangulation(Su), qg::triangulation(J), qg::triangulation(M)})
        worst = std::max(worst, std::fabs(qg::turaev_viro(T, r).value - want) / want);
      s.check("quantum_gravity: TV_" + std::to_string(r) + "(S³) = w^-2 on ∂Δ⁴, the suspension of ∂Δ³, C3 * C3 and the elementary melon (relative)", worst, 1e-11);
    }
  }
  // ---- a non-melonic sphere ---------------------------------------------------------------------------------------------
  {
    // The 4-coloured graph on two black vertices with colours {0, 1} parallel and {2, 3}
    // crossed: degree 1, not melonic, and still a 3-sphere.
    const mt::ColouredGraph P{4, 2, {{0, 1}, {0, 1}, {1, 0}, {1, 0}}};
    const auto X = P.dual_complex();
    const auto pi1 = mt::fundamental_group(X);
    s.require("quantum_gravity: the graph with colours {0, 1} parallel and {2, 3} crossed has degree 1, spherical bubbles and π1 = 1: a 3-sphere",
              P.degree() == 1 && bubbles_are_spheres(P) && homology_text(X) == "ℤ ; 0 ; 0 ; ℤ" && pi1.simplified.presentation.generators == 0);
    double worst = 0;
    for (int r = 3; r <= 7; r++) worst = std::max(worst, std::fabs(qg::turaev_viro(P, r).value * physis::symmetry::total_quantum_dimension_squared(r) - 1));
    s.check("quantum_gravity: TV_r = w^-2 on this triangulation of S³ too, r = 3..7 (relative)", worst, 1e-11);
  }
  // ---- S² × S¹ --------------------------------------------------------------------------------------------------------------
  {
    // The first 4-coloured graph, by enumeration, whose bubbles are spheres (so that it encodes
    // a closed 3-manifold) and whose fundamental group is ℤ: by Kneser's theorem and the
    // Poincaré conjecture, S² × S¹.
    std::optional<mt::ColouredGraph> found;
    for (std::size_t p = 2; p <= 4 && !found; p++)
      qg::for_each_coloured_graph(3, p, [&](const mt::ColouredGraph& G) {
        if (found || !G.connected() || !bubbles_are_spheres(G)) return;
        const auto pi1 = mt::fundamental_group(G.dual_complex());
        if (pi1.simplified.presentation.generators == 1 && pi1.simplified.presentation.relators.empty()) found = G;
      });
    s.require("quantum_gravity: a coloured graph of S² × S¹ is found among those on at most 8 vertices", found.has_value(),
              found ? homology_text(found->dual_complex()) + ", " + std::to_string(2 * found->p) + " vertices" : "none");
    if (found) {
      double worst = 0;
      for (int r = 3; r <= 7; r++) worst = std::max(worst, std::fabs(qg::turaev_viro(*found, r).value - 1));
      s.check("quantum_gravity: TV_r(S² × S¹) = 1 on its coloured graph, r = 3..7", worst, 1e-11);
    }
    const auto SS = sphere_times_circle(3);
    const auto MS = mt::manifold_type(SS);
    s.require("quantum_gravity: ∂Δ³ × C3, triangulated by staircases, is a closed 3-manifold with the homology of S² × S¹",
              MS.manifold && MS.closed && homology_text(SS.cw()) == "ℤ ; ℤ ; ℤ ; ℤ");
    s.check("quantum_gravity: TV_3(S² × S¹) = 1 on the 36 tetrahedra of ∂Δ³ × C3", qg::turaev_viro(SS, 3).value - 1, 1e-11);
  }
  // ---- lens spaces --------------------------------------------------------------------------------------------------------------
  {
    struct Case { std::size_t p, s, t; };
    for (const Case c : {Case{2, 2, 2}, Case{3, 1, 1}, Case{3, 2, 1}, Case{4, 1, 1}, Case{5, 1, 1}}) {
      const auto T = lens(c.p, 1, c.s, c.t);
      double worst = 0;
      std::string values;
      for (int r = 3; r <= 7; r++) {
        if (r == 7 && T.edges > 12) continue;
        const double got = qg::turaev_viro(T, r).value, want = qg::lens_space_tv(static_cast<int>(c.p), r);
        worst = std::max(worst, std::fabs(got - want));
        values += (values.empty() ? "" : ", ") + calx::verification::exponential(got);
      }
      s.check("quantum_gravity: TV_r(L(" + std::to_string(c.p) + ", 1)) = |τ_r|² on the quotient of C" + std::to_string(c.p * c.s) + " * C" +
                  std::to_string(c.p * c.t) + " (" + std::to_string(T.edges) + " edges), r = 3..7: " + values,
              worst, 1e-11);
    }
    // TV_3 of ℝP³ vanishes: the nontrivial class of H¹(ℝP³; ℤ/2) weighs −1.
    s.check("quantum_gravity: TV_3(ℝP³) = 0", qg::turaev_viro(lens(2, 1, 2, 2), 3).value, 1e-12);
  }
  const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  s.require("quantum_gravity: the Turaev–Viro checks run in under a minute", seconds < 60, std::to_string(seconds) + " s");
  return s.done();
}
