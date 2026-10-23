// Checks of physis/quantum_gravity spin_foam.hpp: the spin foams of ∂Δ⁴, of the suspended
// 3-sphere and of a melonic 4-coloured graph are dual 2-complexes with the right maps between
// faces and edges; the sum of their amplitudes over admissible spins is the Turaev–Viro
// invariant at r = 3, 4, 5; vertex amplitudes are the quantum 6j symbols of their spins, and
// tend to the classical ones as r grows; Metropolis sampling with weight |A| reproduces the
// exact distribution of a face's spin and the mean sign; quantum tetrahedra have the
// Cayley–Menger volume, and the Ponzano–Regge formula approaches the 6j symbol at large spins,
// its error halving as the spins double; the regular spin assignment on ∂Δ⁴ is a Regge
// geometry with the deficit angle 2π − 3 arccos(1/3) at every edge.

#include <physis/quantum_gravity/spin_foam.hpp>
#include <morphe/topology/constructions.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <map>
#include <numbers>
#include <string>
#include <vector>

namespace pq = physis::quantum_gravity;
namespace mt = morphe::topology;

namespace {

// The spin foam of a melonic 4-coloured graph: two melons inserted in the elementary one.
mt::ColouredGraph melon() {
  auto G = mt::elementary_melon(3);
  G = mt::insert_melon(G, 1, 0);
  G = mt::insert_melon(G, 3, 1);
  return G;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  const auto S3 = mt::simplex_boundary(4);
  const auto S3b = mt::suspension(mt::simplex_boundary(3));
  const auto fS3 = pq::spin_foam(S3), fS3b = pq::spin_foam(S3b), fM = pq::spin_foam(melon());

  // ---- the 2-complexes and their maps -------------------------------------------------------------------------------------
  {
    bool ok = true;
    std::string why;
    for (const auto* F : {&fS3, &fS3b, &fM}) {
      try {
        F->complex.validate();
      } catch (const std::exception& e) {
        ok = false, why = e.what();
      }
      if (F->faces() != F->triangulation.edges || F->vertices() != F->triangulation.tetrahedra.size()) ok = false, why = "counts";
      // Each vertex's six faces are the faces dual to its tetrahedron's edges, and each face of
      // the 2-complex passes through exactly the vertices whose tetrahedra contain its edge.
      for (std::size_t v = 0; v < F->vertices(); v++)
        for (int k = 0; k < 6; k++) {
          const std::size_t f = F->vertex_faces[v][static_cast<std::size_t>(k)];
          if (F->face_edge[f] != F->triangulation.tetrahedra[F->vertex_tetrahedron[v]][static_cast<std::size_t>(k)]) ok = false, why = "vertex faces";
          const auto& word = F->complex.faces[f].boundary;
          bool through = false;
          for (const auto& st : word) {
            const auto& ed = F->complex.edges[st.edge];
            through = through || ed[0] == v || ed[1] == v;
          }
          if (!through) ok = false, why = "a face does not pass through its vertex";
        }
    }
    s.require("spin foam: the foams of ∂Δ⁴, the suspended 3-sphere and a melonic graph are 2-complexes whose faces are dual to the edges and vertices to the tetrahedra", ok, why);
    s.require("spin foam: ∂Δ⁴ gives 5 vertices, 10 edges and 10 faces", fS3.complex.counts() == std::vector<std::size_t>{5, 10, 10});
  }
  // ---- amplitudes summed: the Turaev–Viro invariant -----------------------------------------------------------------------
  {
    double worst = 0;
    std::string detail;
    for (const auto* F : {&fS3, &fS3b, &fM})
      for (int r = 3; r <= 5; r++) {
        double sum = 0;
        pq::for_each_state(*F, r, [&](const std::vector<int>&, double a) { sum += a; });
        const double tv = pq::turaev_viro(F->triangulation, r).value;
        worst = std::max(worst, std::fabs(sum - tv) / std::fabs(tv));
        detail += std::to_string(sum).substr(0, 7) + " ";
      }
    s.check("spin foam: Σ_j A(j) is the Turaev-Viro invariant on three foams at r = 3, 4, 5 (relative)", worst, 1e-12);
  }
  // ---- vertex amplitudes --------------------------------------------------------------------------------------------------
  {
    double worst = 0;
    std::size_t states = 0;
    pq::for_each_state(fS3, 5, [&](const std::vector<int>& j, double) {
      if (++states % 37) return;
      const auto A = pq::amplitude(fS3, j, 5);
      for (std::size_t v = 0; v < fS3.vertices(); v++) {
        const auto x = fS3.vertex_spins(j, v);
        worst = std::max(worst, std::fabs(A.vertex[v] - physis::symmetry::quantum_6j(5, x[0], x[1], x[3], x[5], x[4], x[2])));
      }
    });
    s.check("spin foam: vertex amplitudes are the quantum 6j symbols {01 02 12; 23 13 03} of their spins", worst, 1e-15);
    // As r grows, the quantum weights tend to the classical ones (Ponzano-Regge), the
    // difference falling as (π/r)², since [n] = n (1 − (n² − 1) π² / (6 r²) + O(r⁻⁴)).
    const std::vector<int> j{2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
    const auto c = pq::amplitude(fS3, j, 0);
    auto error = [&](int r) {
      const auto q = pq::amplitude(fS3, j, r);
      return std::max(std::fabs(q.vertices / c.vertices - 1), std::fabs(q.faces / c.faces - 1));
    };
    s.check("spin foam: at r = 2000 the vertex and face weights are the classical ones (relative)", error(2000), 1e-4);
    s.check("spin foam: ... with an error falling as 1/r² (error(1000) / error(2000) - 4)", error(1000) / error(2000) - 4, 0.01);
  }
  // ---- sampling against exact enumeration ---------------------------------------------------------------------------------
  {
    const int r = 4;
    std::vector<double> exact(static_cast<std::size_t>(r - 1), 0.0);
    double sum_abs = 0, sum = 0;
    pq::for_each_state(fS3, r, [&](const std::vector<int>& j, double a) {
      exact[static_cast<std::size_t>(j[0])] += std::fabs(a);
      sum_abs += std::fabs(a);
      sum += a;
    });
    for (auto& e : exact) e /= sum_abs;
    pq::SampleOptions o;
    o.r = r;
    o.sweeps = 40000;
    o.thermalise = 500;
    o.seed = 3;
    const auto S = pq::sample(fS3, o);
    double dev = 0;
    for (int x = 0; x <= r - 2; x++) dev = std::max(dev, std::fabs(S.histogram[0][static_cast<std::size_t>(x)] - exact[static_cast<std::size_t>(x)]));
    // 40000 correlated sweeps of a chain accepting about half its moves: the standard error of a
    // frequency near 0.3 is below 0.005; the limit allows four of them.
    s.check("spin foam: Metropolis with weight |A| on ∂Δ⁴ at r = 4 reproduces the exact distribution of a face's spin (largest difference)", dev, 0.02);
    s.require("spin foam: the sampler moves (acceptance between 5% and 95%)", S.acceptance > 0.05 && S.acceptance < 0.95, std::to_string(S.acceptance));
  }
  {
    // At r = 6 the amplitudes of the suspended 3-sphere alternate (exactly, Σ A / Σ |A| = 0.949):
    // the sampler's mean sign and a face's distribution against enumeration of its 52769 states.
    const int r = 6;
    std::vector<double> exact(static_cast<std::size_t>(r - 1), 0.0);
    double sum_abs = 0, sum = 0;
    pq::for_each_state(fS3b, r, [&](const std::vector<int>& j, double a) {
      exact[static_cast<std::size_t>(j[3])] += std::fabs(a);
      sum_abs += std::fabs(a);
      sum += a;
    });
    for (auto& e : exact) e /= sum_abs;
    pq::SampleOptions o;
    o.r = r;
    o.sweeps = 40000;
    o.thermalise = 500;
    o.seed = 11;
    const auto S = pq::sample(fS3b, o);
    double dev = 0;
    for (int x = 0; x <= r - 2; x++) dev = std::max(dev, std::fabs(S.histogram[3][static_cast<std::size_t>(x)] - exact[static_cast<std::size_t>(x)]));
    s.check("spin foam: on the suspended 3-sphere at r = 6 the sampler's mean sign is Σ A / Σ |A| = " + std::to_string(sum / sum_abs).substr(0, 5), S.mean_sign - sum / sum_abs, 0.02);
    s.check("spin foam: ... and its distribution of a face's spin is the exact one (largest difference)", dev, 0.02);
  }
  // ---- quantum tetrahedra --------------------------------------------------------------------------------------------------
  {
    double worst = 0;
    for (const auto& j : std::vector<std::array<int, 6>>{{4, 4, 4, 4, 4, 4}, {2, 3, 3, 3, 3, 2}, {6, 5, 7, 5, 6, 4}, {10, 12, 10, 12, 10, 12}}) {
      const auto Q = pq::tetrahedron_geometry(j);
      if (!Q.euclidean) { worst = INFINITY; continue; }
      physis::spacetime::SquaredLengths sq(4, std::vector<double>(4, 0.0));
      for (int k = 0; k < 6; k++) {
        const auto [a, b] = pq::tetrahedron_slots[static_cast<std::size_t>(k)];
        sq[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] = sq[static_cast<std::size_t>(b)][static_cast<std::size_t>(a)] = Q.lengths[static_cast<std::size_t>(k)] * Q.lengths[static_cast<std::size_t>(k)];
      }
      worst = std::max(worst, std::fabs(Q.volume / physis::spacetime::cayley_menger_volume(sq) - 1));
    }
    const auto R = pq::tetrahedron_geometry({4, 4, 4, 4, 4, 4});
    const double l = 2.5;
    s.check("spin foam: quantum tetrahedra have the Cayley-Menger volume (relative)", worst, 1e-12);
    s.check("spin foam: the regular tetrahedron of spin 2 has volume ℓ³/(6√2), ℓ = 5/2 (relative)", R.volume / (l * l * l / (6 * std::sqrt(2.0))) - 1, 1e-12);
    // Opposite edges of length 9/2 and the other four 5/2: every triangle exists, the tetrahedron
    // does not (the two long edges cannot be that far apart), and its 6j symbol lies in the
    // classically forbidden region.
    s.require("spin foam: admissible spins whose lengths bound no tetrahedron are recognised as non-Euclidean", !pq::tetrahedron_geometry({8, 4, 4, 4, 4, 8}).euclidean);
  }
  {
    // At spins λ j the error of the Ponzano-Regge formula, relative to its envelope 1/√(12πV),
    // falls as 1/λ (Roberts 1999): the largest over eight consecutive scales, at λ = 20..27 and
    // 40..47.
    const std::vector<std::array<int, 6>> shapes{{2, 2, 2, 2, 2, 2}, {2, 3, 3, 3, 3, 2}};
    double small = 0, large = 0, ratio = 0;
    for (const auto& sh : shapes) {
      double e[2] = {0, 0};
      for (int pass = 0; pass < 2; pass++)
        for (int lam = pass ? 40 : 20; lam < (pass ? 48 : 28); lam++) {
          std::array<int, 6> j;
          for (int k = 0; k < 6; k++) j[static_cast<std::size_t>(k)] = sh[static_cast<std::size_t>(k)] * lam;
          const auto Q = pq::quantum_tetrahedron(j, 0);
          e[pass] = std::max(e[pass], std::fabs(Q.six_j - Q.ponzano_regge) / Q.envelope);
        }
      small = std::max(small, e[0]);
      large = std::max(large, e[1]);
      ratio = std::max(ratio, std::fabs(e[1] / e[0] - 0.5));
    }
    s.check("spin foam: the Ponzano-Regge formula at spins 20-27 (error relative to the envelope)", small, 0.05);
    s.check("spin foam: its error halves when the spins double (ratio - 1/2)", ratio, 0.15);
  }
  // ---- Regge geometry -------------------------------------------------------------------------------------------------------
  {
    const std::vector<int> j(10, 4);
    const auto R = pq::regge_geometry(fS3, j);
    const double eps = 2 * std::numbers::pi - 3 * std::acos(1.0 / 3.0);
    double worst = R.complete ? 0 : INFINITY;
    for (double d : R.deficit) worst = std::max(worst, std::fabs(d - eps));
    s.check("spin foam: the regular assignment j = 2 on ∂Δ⁴ is a Regge geometry with deficit 2π - 3 arccos(1/3) at every edge", worst, 1e-12);
    s.check("spin foam: ... and Regge action 10 ℓ ε, ℓ = 5/2 (relative)", R.action / (10 * 2.5 * eps) - 1, 1e-12);
  }
  return s.done();
}
