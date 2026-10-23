// Checks of physis/quantum_gravity beyond the Turaev–Viro sum: BF theory with finite groups
// (flat connections through π1 against enumeration, Z_BF of surfaces against Mednykh's
// formula Σ_ρ (dim ρ)^χ, the Boulatov amplitudes of S³ and S² × S¹ graphs, Dijkgraaf–Witten of
// lens spaces); the Feynman expansion of coloured tensor models (Gurau's scaling, the melonic
// leading order against the Fuss–Catalan numbers); and matrix models (gluings of polygons by
// genus against the Harer–Zagier recursion and the Catalan numbers, the genus of each gluing
// against morphe's classification of surfaces, and GUE moments by Monte Carlo).

#include <physis/quantum_gravity/bf.hpp>
#include <physis/quantum_gravity/matrix_models.hpp>
#include <physis/quantum_gravity/tensor_models.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <morphe/topology/constructions.hpp>
#include <morphe/topology/surfaces.hpp>

#include <cmath>
#include <numeric>
#include <string>

namespace qg = physis::quantum_gravity;
namespace mt = morphe::topology;
namespace ma = morphe::algebra;
using morphe::arithmetic::Integer;
using morphe::arithmetic::Rational;

int main() {
  calx::verification::Suite s;

  // ---- BF theory -----------------------------------------------------------------------------------------------------------
  {
    struct Case { mt::CWComplex X; ma::FiniteGroup G; };
    const std::vector<Case> cases = {
        {mt::polygon_complex("abAB"), ma::symmetric(3)},          {mt::torus_lattice({2, 2}), ma::cyclic(3)},
        {mt::polygon_complex("aa"), ma::quaternion()},            {mt::klein_bottle(), ma::symmetric(3)},
        {mt::dunce_hat(), ma::symmetric(3)},                      {mt::elementary_melon(3).two_complex(), ma::symmetric(3)},
        {mt::polygon_complex("abAB acAC"), ma::dihedral(4)},      {mt::real_projective_space(2), ma::cyclic(4)},
    };
    int bad = 0;
    std::string detail;
    for (const auto& c : cases) {
      const Integer a = qg::count_flat(c.X, c.G), b = qg::count_flat_brute(c.X, c.G);
      if (a != b) bad++, detail += c.X.name + ": " + a.str() + " vs " + b.str() + "; ";
    }
    s.require("quantum_gravity: flat connections through π1 agree with enumeration on 8 complexes and groups", bad == 0, detail);
    // Z_BF(Σ_g; G) = Σ_ρ (dim ρ)^χ: S3 has degrees 1, 1, 2; Q8 has 1, 1, 1, 1, 2.
    auto mednykh = [](const std::vector<int>& dims, long long chi) {
      Rational z = 0;
      for (int d : dims) z += morphe::arithmetic::pow(Rational(d), chi);
      return z;
    };
    bool ok = true;
    std::string got;
    for (std::size_t g = 0; g <= 3; g++) {
      const auto X = mt::orientable_surface(g);
      const long long chi = 2 - 2 * static_cast<long long>(g);
      const Rational zs = qg::bf_partition_function(X, ma::symmetric(3)), zq = qg::bf_partition_function(X, ma::quaternion());
      ok = ok && zs == mednykh({1, 1, 2}, chi) && zq == mednykh({1, 1, 1, 1, 2}, chi);
      got += zs.str() + " ";
    }
    s.require("quantum_gravity: Z_BF of closed surfaces of genus 0..3 is Σ_ρ (dim ρ)^χ for S3 and Q8 (Mednykh)", ok, got);
    // Boulatov amplitudes with a finite group: |G|^{v−1} |Hom(π1 M, G)|, v the 3-bubbles.
    const ma::FiniteGroup S3 = ma::symmetric(3);
    const auto melon = mt::elementary_melon(3);
    s.require("quantum_gravity: the Boulatov amplitude of the elementary melon (S³, 4 vertices) with S3 is 6³ = 216",
              qg::boulatov_amplitude(melon, S3) == Rational(216), qg::boulatov_amplitude(melon, S3).str());
    const mt::ColouredGraph pillow{4, 2, {{0, 1}, {0, 1}, {1, 0}, {1, 0}}};
    s.require("quantum_gravity: ... and of the non-melonic sphere graph, with 4 bubbles, the same 216", qg::boulatov_amplitude(pillow, S3) == Rational(216));
    s.require("quantum_gravity: the Dijkgraaf–Witten invariant of L(4, q) with ℤ/6 is gcd(4, 6) / 6 = 1/3",
              qg::dijkgraaf_witten(mt::lens_space(4), ma::cyclic(6)) == Rational(1, 3));
  }
  // ---- tensor models -----------------------------------------------------------------------------------------------------------
  {
    // Gurau's scaling: F − D(D−1)p/2 = D − 2ω/(D−1)! on every connected graph.
    int bad = 0;
    std::size_t graphs = 0;
    for (std::size_t D : {2, 3, 4})
      for (std::size_t p = 1; p <= (D == 4 ? 3u : 4u); p++)
        qg::for_each_coloured_graph(D, p, [&](const mt::ColouredGraph& G) {
          if (!G.connected()) return;
          graphs++;
          std::size_t fact = 1;
          for (std::size_t k = 2; k < D; k++) fact *= k;
          const long long want = static_cast<long long>(D) * static_cast<long long>(fact) - 2 * static_cast<long long>(G.degree());
          if (qg::amplitude_exponent(G) * static_cast<long long>(fact) != want) bad++;
        });
    s.require("quantum_gravity: every connected coloured graph scales as N^{D − 2ω/(D−1)!} (D = 2, 3, 4)", bad == 0,
              std::to_string(graphs) + " graphs, " + std::to_string(bad) + " failures");
    // The leading order of the free energy at (λλ̄)^p in D = 3: Fuss–Catalan(4, p) / p.
    bool leading = true;
    std::string coeffs;
    for (long long p = 1; p <= 4; p++) {
      const auto F = qg::free_energy(3, static_cast<std::size_t>(p));
      const Rational want(qg::fuss_catalan(4, p), Integer(p));
      leading = leading && F.count(0) && F.at(0) == want;
      coeffs += (F.count(0) ? F.at(0).str() : "?") + " ";
    }
    s.require("quantum_gravity: the melonic free energy of the D = 3 tensor model is 1, 2, 22/3, 35 at orders 1..4 (Fuss–Catalan / p)", leading, coeffs);
  }
  // ---- matrix models -----------------------------------------------------------------------------------------------------------
  {
    bool hz = true, catalan = true;
    std::string table;
    for (long long n = 1; n <= 7; n++) {
      const auto eps = qg::gluings_by_genus(static_cast<std::size_t>(n));
      for (long long g = 0; 2 * g <= n; g++) hz = hz && eps[static_cast<std::size_t>(g)] == qg::harer_zagier(g, n);
      catalan = catalan && eps[0] == morphe::arithmetic::exact_div(morphe::arithmetic::binomial(2 * n, n), Integer(n + 1));
      if (n == 4) for (const auto& e : eps) table += e.str() + " ";
    }
    s.require("quantum_gravity: gluings of the 2n-gon by genus follow the Harer–Zagier recursion, n <= 7 (n = 4: 14 70 21)", hz, table);
    s.require("quantum_gravity: the planar gluings are counted by the Catalan numbers", catalan);
    int mismatch = 0;
    for (std::size_t n = 1; n <= 4; n++)
      qg::for_each_pairing(n, [&](const qg::Pairing& pi) {
        const auto S = mt::classify_surface(qg::gluing_complex(pi));
        if (!S.surface || !S.closed || !S.orientable || S.genus != qg::gluing_genus(pi)) mismatch++;
      });
    s.require("quantum_gravity: each gluing's genus agrees with morphe's classification of the glued surface, n <= 4", mismatch == 0, std::to_string(mismatch) + " mismatches");
    // ⟨tr M^{2n}⟩ of the GUE of size 4, by Monte Carlo: within 5 standard errors.
    const std::size_t N = 4, samples = 40000;
    for (std::size_t n = 1; n <= 3; n++) {
      double sum = 0, sum2 = 0;
      calx::random::Stream zn(2026, 10, 23, static_cast<std::uint32_t>(n));
      for (std::size_t k = 0; k < samples; k++) {
        const double t = qg::trace_power(qg::gue_sample(N, zn), N, 2 * n);
        sum += t, sum2 += t * t;
      }
      const double mean = sum / samples, se = std::sqrt((sum2 / samples - mean * mean) / samples);
      const double want = qg::gue_moment(static_cast<long long>(n), static_cast<long long>(N)).to_double();
      s.check("quantum_gravity: ⟨tr M^" + std::to_string(2 * n) + "⟩ of the 4 × 4 GUE against Σ_g ε_g N^{n+1−2g} = " + std::to_string(static_cast<long long>(want)) +
                  " (in standard errors)",
              (mean - want) / se, 5);
    }
  }
  return s.done();
}
