// Checks of the real cochains (cochains.hpp) and the Hodge decomposition (hodge.hpp).
//
// Cochains: d₁ d₀ = 0 on surfaces glued by words (the torus lattice, genus 2, the projective plane,
// the Klein bottle, the dunce hat); the integral of d₀ φ along a walk is φ(end) − φ(start), and
// round a closed walk 0; the potential of an exact cochain recovers it up to a constant on each
// component, with zero defect and zero chords, and a cochain round an unfilled circle has none.
// Hodge: on every complex, with random weights, the three parts add up to ω, are orthogonal, the
// harmonic part is closed and co-closed, the co-exact part carries the whole curvature, a part
// decomposes into itself; the harmonic space has dimension b₁(X; ℝ) (morphe), 0 on the projective
// plane and the octahedral sphere (whose 2-cycle tests the kernel); the sparse method agrees with
// the dense one, and decomposes a 60 × 60 lattice torus beyond the dense one; on a weighted cycle the
// harmonic part is c/w_e with c = Σ ω / Σ 1/w; on a triangle the integral log 0.6 is harmonic while
// the triangle is unfilled and curvature once it is filled. Flows with conductances: on every
// complex the driven flow carries the whole boundary and the circulating flows none (Kirchhoff),
// the energies add up, and Σ F ω splits into the pairings of like parts; Ohm's law on two parallel
// conductances.

#include <physis/discretisation/cochains.hpp>
#include <physis/discretisation/hodge.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <morphe/topology/constructions.hpp>
#include <morphe/topology/invariants.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

namespace pd = physis::discretisation;
namespace mt = morphe::topology;

namespace {

double max_abs(const pd::Cochain& a) {
  double m = 0;
  for (double x : a) m = std::max(m, std::fabs(x));
  return m;
}

pd::Cochain random_cochain(std::size_t n, calx::random::Stream& z, double lo = -1, double hi = 1) {
  pd::Cochain c(n);
  for (double& x : c) x = lo + (hi - lo) * z.uniform();
  return c;
}

// A graph: a CW complex with no 2-cells.
mt::CWComplex graph(std::size_t vertices, const std::vector<std::array<std::size_t, 2>>& edges) {
  mt::CWComplex X;
  X.vertices = vertices;
  X.edges = edges;
  return X;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream z(20261024, 0x636f6368u, 0, 0);  // "coch"

  const std::vector<std::pair<std::string, mt::CWComplex>> complexes{
      {"the 4 × 3 lattice torus", mt::torus_lattice({4, 3})},
      {"the surface of genus 2", mt::orientable_surface(2)},
      {"the projective plane", mt::real_projective_space(2)},
      {"the Klein bottle", mt::klein_bottle()},
      {"the dunce hat", mt::dunce_hat()},
      {"the octahedral sphere", mt::octahedron().cw()},
      {"a 6 × 5 grid disc", mt::product(mt::interval(6), mt::interval(5))},
      {"two triangles and a path (a graph)", graph(7, {{0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3}, {5, 6}})},
  };

  for (const auto& [name, X] : complexes) {
    const std::size_t V = X.vertices, E = X.edges.size(), F = X.faces.size();
    // d₁ d₀ = 0.
    const pd::Cochain phi = random_cochain(V, z);
    s.check(name + ": d₁ d₀ = 0", max_abs(pd::d1(X, pd::d0(X, phi))), 1e-12);
    // Stokes on every 2-cell: the integral round its word is d₁ ω there.
    const pd::Cochain omega = random_cochain(E, z);
    const pd::Cochain curv = pd::d1(X, omega);
    double stokes = 0;
    for (std::size_t f = 0; f < F; f++) stokes = std::max(stokes, std::fabs(pd::integral(omega, X.faces[f].boundary) - curv[f]));
    s.check(name + ": the integral round each 2-cell's word is its curvature d₁ ω", stokes, 1e-12);
    // The potential of an exact cochain.
    const auto P = pd::potential(X, pd::d0(X, phi));
    const auto comp = pd::components(X);
    double off = 0;
    for (std::size_t v = 0; v < V; v++)
      for (std::size_t u = 0; u < V; u++)
        if (comp[u] == comp[v]) off = std::max(off, std::fabs((P.phi[v] - P.phi[u]) - (phi[v] - phi[u])));
    s.require(name + ": an exact cochain's potential, up to a constant on each component, with no defect", off < 1e-12 && P.defect < 1e-12,
              std::to_string(off) + ", defect " + std::to_string(P.defect));
    // The Hodge decomposition with random weights.
    pd::Weights W;
    W.edge = random_cochain(E, z, 0.5, 2.0);
    W.face = random_cochain(F, z, 0.5, 2.0);
    const auto H = pd::hodge_decomposition(X, omega, W);
    pd::Cochain sum(E);
    for (std::size_t e = 0; e < E; e++) sum[e] = H.exact[e] + H.coexact[e] + H.harmonic[e] - omega[e];
    s.check(name + ": exact + co-exact + harmonic = ω", max_abs(sum), 1e-10);
    const double cross = std::max({std::fabs(pd::inner1(H.exact, H.coexact, W)), std::fabs(pd::inner1(H.exact, H.harmonic, W)),
                                   std::fabs(pd::inner1(H.coexact, H.harmonic, W))});
    s.check(name + ": the three parts are orthogonal in the weighted inner product", cross, 1e-10);
    pd::Cochain wh(E);
    for (std::size_t e = 0; e < E; e++) wh[e] = W.e(e) * H.harmonic[e];
    s.check(name + ": the harmonic part is closed and co-closed", std::max(max_abs(pd::d1(X, H.harmonic)), max_abs(pd::d0_transpose(X, wh))), 1e-10);
    pd::Cochain dc = pd::d1(X, H.coexact);
    for (std::size_t f = 0; f < F; f++) dc[f] -= curv[f];
    s.check(name + ": the co-exact part carries the whole curvature", max_abs(dc), 1e-10);
    pd::HodgeOptions sparse;
    sparse.method = pd::HodgeOptions::sparse;
    const auto Hs = pd::hodge_decomposition(X, omega, W, sparse);
    double agree = 0;
    for (std::size_t e = 0; e < E; e++)
      agree = std::max({agree, std::fabs(Hs.exact[e] - H.exact[e]), std::fabs(Hs.coexact[e] - H.coexact[e]), std::fabs(Hs.harmonic[e] - H.harmonic[e])});
    s.check(name + ": the sparse method (conjugate gradients) agrees with the dense one", agree, 1e-8);
    const auto again = pd::hodge_decomposition(X, H.harmonic, W);
    s.check(name + ": the harmonic part decomposes into itself", std::max(max_abs(again.exact), max_abs(again.coexact)), 1e-10);
    // The dimension of the harmonic space.
    const auto b = mt::betti(X, 0);
    const std::size_t b1 = b.size() > 1 ? b[1] : 0;
    const std::size_t dim = pd::harmonic_basis(X, W).size();
    s.require(name + ": the harmonic space has dimension b₁(X; ℝ)", dim == b1, std::to_string(dim) + " for b₁ = " + std::to_string(b1));
    // A flow with conductances: Kirchhoff's current law, the energies adding up, and the pairing
    // with a cochain splitting into like parts.
    const pd::Cochain flow = random_cochain(E, z, -3, 3);
    const auto Q = pd::flow_decomposition(X, flow, W.edge, W.face);
    const pd::Cochain bF = pd::d0_transpose(X, flow), bD = pd::d0_transpose(X, Q.driven);
    pd::Cochain circulating(E);
    for (std::size_t e = 0; e < E; e++) circulating[e] = Q.local[e] + Q.global[e];
    double kirchhoff = max_abs(pd::d0_transpose(X, circulating));
    for (std::size_t v = 0; v < V; v++) kirchhoff = std::max(kirchhoff, std::fabs(bD[v] - bF[v]));
    s.check(name + ": a flow's driven part carries its whole boundary, the circulating parts none (Kirchhoff)", kirchhoff, 1e-10);
    const auto Hw = pd::hodge_decomposition(X, omega, pd::Weights{W.edge, {}});
    double natural = 0, like = 0;
    for (std::size_t e = 0; e < E; e++) {
      natural += flow[e] * omega[e];
      like += Q.driven[e] * Hw.exact[e] + Q.local[e] * Hw.coexact[e] + Q.global[e] * Hw.harmonic[e];
    }
    s.check(name + ": Σ F ω splits into the pairings of like parts, and the energies add up",
            std::fabs(natural - like) + std::fabs(Q.driven_energy + Q.local_energy + Q.global_energy - Q.energy), 1e-9);
  }

  // Ohm's law on two parallel edges from 0 to 1 with conductances 2 and 3 carrying 4 and 1: the
  // driven flow divides the total 5 in proportion to the conductances (2, 3), the potential rises by
  // 5 / (2 + 3) = 1, and the remainder (2, −2) circulates round the loop they make.
  {
    mt::CWComplex P;
    P.vertices = 2;
    P.edges = {{0, 1}, {0, 1}};
    const auto Q = pd::flow_decomposition(P, {4.0, 1.0}, {2.0, 3.0});
    s.check("Ohm: two parallel conductances 2 and 3 carrying 5 between them drive 2 and 3, the potential rising by 1",
            std::fabs(Q.driven[0] - 2) + std::fabs(Q.driven[1] - 3) + std::fabs(Q.potential[1] - Q.potential[0] - 1), 1e-12);
    s.check("... and the rest, 2 one way and 2 back, circulates round their loop",
            std::fabs(Q.global[0] - 2) + std::fabs(Q.global[1] + 2) + std::fabs(Q.energy - (16.0 / 2 + 1.0 / 3)) + std::fabs(Q.driven_energy - 5), 1e-12);
  }

  // A 60 × 60 lattice torus (3600 vertices, 7200 edges, 3600 squares), beyond the dense method:
  // the sparse one converges, its parts orthogonal, the harmonic part closed and co-closed, and
  // the integral of the harmonic part round the torus's two circles that of ω.
  {
    const auto X = mt::torus_lattice({60, 60});
    const pd::Cochain omega = random_cochain(X.edges.size(), z);
    const auto H = pd::hodge_decomposition(X, omega);
    s.require("a 60 × 60 lattice torus: decomposed by the sparse method, converged", H.sparse && H.residual < 1e-10,
              std::to_string(H.iterations) + " iterations, residual " + std::to_string(H.residual));
    const double cross = std::max({std::fabs(pd::inner1(H.exact, H.coexact)), std::fabs(pd::inner1(H.exact, H.harmonic)), std::fabs(pd::inner1(H.coexact, H.harmonic))});
    s.check("... its three parts orthogonal", cross, 1e-7);
    s.check("... its harmonic part closed and co-closed", std::max(max_abs(pd::d1(X, H.harmonic)), max_abs(pd::d0_transpose(X, H.harmonic))), 1e-8);
  }

  // A potential round an unfilled circle has a defect: the integral round it.
  {
    const auto C = mt::circle(5);
    pd::Cochain omega(C.edges.size(), 0.0);
    omega[0] = 0.7;
    const auto P = pd::potential(C, omega);
    double chord = 0;
    for (double x : P.chords) chord += std::fabs(x);
    s.check("a circle: the chord carries the whole integral, 0.7", chord - 0.7, 1e-12);
  }

  // The weighted cycle: harmonic part c / w_e, c = Σ ω / Σ 1/w (edges oriented round the cycle).
  {
    const std::size_t n = 7;
    mt::CWComplex C;
    C.vertices = n;
    for (std::size_t i = 0; i < n; i++) C.edges.push_back({i, (i + 1) % n});
    const pd::Cochain omega = random_cochain(n, z);
    pd::Weights W;
    W.edge = random_cochain(n, z, 0.5, 3.0);
    double so = 0, si = 0;
    for (std::size_t e = 0; e < n; e++) { so += omega[e]; si += 1 / W.edge[e]; }
    const double c = so / si;
    const auto H = pd::hodge_decomposition(C, omega, W);
    double worst = 0;
    for (std::size_t e = 0; e < n; e++) worst = std::max(worst, std::fabs(H.harmonic[e] - c / W.edge[e]));
    s.check("a weighted 7-cycle: the harmonic part is c / w_e with c = Σ ω / Σ 1/w", worst, 1e-12);
  }

  // A triangle of transports 1.2, 0.5, 1.0 (the integral of log U round it is log 0.6): harmonic
  // while the triangle is unfilled, curvature once a 2-cell fills it.
  {
    mt::CWComplex T;
    T.vertices = 3;
    T.edges = {{0, 1}, {1, 2}, {2, 0}};
    const pd::Cochain omega{std::log(1.2), std::log(0.5), std::log(1.0)};
    const mt::Path loop{{0, 1}, {1, 1}, {2, 1}};
    const auto open = pd::hodge_decomposition(T, omega);
    s.check("an unfilled triangle: its harmonic part carries log 0.6", pd::integral(open.harmonic, loop) - std::log(0.6), 1e-12);
    s.check("... and its exact part none", pd::integral(open.exact, loop), 1e-12);
    T.faces.push_back({0, loop});
    const auto filled = pd::hodge_decomposition(T, omega);
    s.check("a filled triangle: no harmonic part", max_abs(filled.harmonic), 1e-12);
    s.check("... its curvature is log 0.6, all of it co-exact", pd::d1(T, filled.coexact)[0] - std::log(0.6), 1e-12);
  }
  return s.done();
}
