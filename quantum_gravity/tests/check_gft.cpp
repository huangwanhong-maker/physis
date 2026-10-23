// Checks of physis/quantum_gravity/gft.hpp: flat connections on the Feynman 2-complexes of
// group field theories counted through π1 and by brute force (ℤ₂, ℤ₃, S₃); the amplitudes of
// melonic graphs, |G|^{p+2} in Boulatov's model and |G|^{3p+3} in Ooguri's; the abelian
// amplitudes from the homology of the 2-complex; without the gauge projection the amplitude
// |G|^{F₀} of the U(N) tensor model with N = |G|, so that the free energy of the group field
// theory without projection is the tensor model's series at N = |G|; with the trivial group,
// the count of graphs.

#include <physis/field_theory/invariant.hpp>
#include <physis/quantum_gravity/gft.hpp>
#include <calx/verification/check.hpp>

#include <string>
#include <vector>

namespace pq = physis::quantum_gravity;
namespace pf = physis::field_theory;
namespace ma = morphe::algebra;
namespace mt = morphe::topology;
using morphe::arithmetic::Integer;
using morphe::arithmetic::Rational;

int main() {
  calx::verification::Suite s;
  const ma::FiniteGroup Z2 = ma::cyclic(2), Z3 = ma::cyclic(3), S3 = ma::symmetric(3), one = ma::cyclic(1);

  // ---- flat connections, two ways ----------------------------------------------------------------------------------------
  {
    int bad = 0, complexes = 0;
    const auto T = pf::quartic_melonic_model(3, true);
    for (const std::vector<std::size_t>& counts : std::vector<std::vector<std::size_t>>{{1, 0, 0}, {1, 1, 0}, {2, 0, 0}, {0, 1, 1}}) {
      const auto E = pf::expand(T, counts);
      for (const auto& C : E.classes) {
        const auto X = pq::tensorial_two_complex(E, C.sigma0, 3);
        complexes++;
        for (const auto* G : {&Z2, &Z3, &S3})
          if (pq::count_flat(X, *G) != pq::count_flat_brute(X, *G)) bad++;
      }
    }
    for (std::size_t D : {3, 4})
      for (std::size_t p = 1; p <= 2; p++)
        pq::for_each_coloured_graph(D, p, [&](const pq::ColouredGraph& G) {
          if (!G.connected()) return;
          const auto X = G.two_complex();
          complexes++;
          if (pq::count_flat(X, Z2) != pq::count_flat_brute(X, Z2)) bad++;
          if (D == 3 && pq::count_flat(X, Z3) != pq::count_flat_brute(X, Z3)) bad++;
          if (D == 3 && p == 1 && pq::count_flat(X, S3) != pq::count_flat_brute(X, S3)) bad++;
        });
    s.check("gft: flat connections on " + std::to_string(complexes) + " Feynman 2-complexes, through π1 and by brute force (ℤ₂, ℤ₃, S₃; failures)", bad, 0);
  }
  // ---- melonic amplitudes ------------------------------------------------------------------------------------------------
  {
    int bad = 0, melonic = 0;
    for (std::size_t D : {3, 4})
      for (std::size_t p = 1; p <= (D == 3 ? 3u : 2u); p++)
        pq::for_each_coloured_graph(D, p, [&](const pq::ColouredGraph& G) {
          if (!G.connected() || G.degree() != 0) return;
          melonic++;
          for (const auto* H : {&Z2, &Z3, &S3}) {
            const long long e = D == 3 ? static_cast<long long>(p) + 2 : 3 * static_cast<long long>(p) + 3;
            if (pq::gft_amplitude(G.two_complex(), *H) != Rational(morphe::arithmetic::pow(Integer(static_cast<long long>(H->order())), e))) bad++;
          }
        });
    s.check("gft: melonic graphs have the amplitude |G|^{p+2} (Boulatov) and |G|^{3p+3} (Ooguri) (" + std::to_string(melonic) + " graphs; failures)", bad, 0);
  }
  // ---- abelian groups from homology --------------------------------------------------------------------------------------
  {
    int bad = 0;
    for (std::uint64_t seed = 1; seed <= 25; seed++) {
      const auto G = mt::random_coloured_graph(3, 3 + seed % 3, seed);
      if (!G.connected()) continue;
      const auto X = G.two_complex();
      const auto H = mt::homology(X, false);
      for (std::size_t N = 2; N <= 6; N++) {
        // |Hom(H1, ℤ_N)| = N^{b1} Π gcd(t, N), and A = N^{F − E + V − 1} |Hom|.
        Integer homs = morphe::arithmetic::pow(Integer(static_cast<long long>(N)), H[1].rank);
        for (const auto& t : H[1].torsion) homs *= morphe::arithmetic::gcd(t, Integer(static_cast<long long>(N)));
        const long long e = static_cast<long long>(X.faces.size()) - static_cast<long long>(X.edges.size()) + static_cast<long long>(X.vertices) - 1;
        const Rational want = Rational(homs) * morphe::arithmetic::pow(Rational(Integer(static_cast<long long>(N))), e);
        if (pq::gft_amplitude(X, ma::cyclic(N)) != want) bad++;
      }
    }
    s.check("gft: amplitudes with ℤ_N, N = 2..6, on random coloured graphs, from the homology of the 2-complex (failures)", bad, 0);
  }
  // ---- tensor models are group field theories without holonomy -----------------------------------------------------------
  {
    int bad = 0;
    const auto T = pf::quartic_melonic_model(3, false);
    const auto S = pf::series(T, 3);
    const auto M = pq::quartic_melonic_gft(3, false);
    for (std::size_t n = 1; n <= 3; n++)
      for (std::size_t N : {1, 2, 3}) {
        const Rational want = pf::evaluate(S.free_energy[n], Rational(static_cast<long long>(N)));
        if (pq::gft_free_energy(M, n, ma::cyclic(N), false) != want) bad++;
      }
    s.check("gft: without the gauge projection the free energy is the U(N) tensor model's at N = |G|, orders 1..3 (failures)", bad, 0);
    // With the trivial group the gauge projection is void: the same count at N = 1.
    bad = 0;
    for (std::size_t n = 1; n <= 3; n++)
      if (pq::gft_free_energy(M, n, one) != pf::evaluate(S.free_energy[n], Rational(1))) bad++;
    for (std::size_t p = 1; p <= 3; p++) {
      Rational count = 0;
      for (const auto& [omega, c] : pq::free_energy(3, p)) count += c;
      if (pq::gft_free_energy(pq::boulatov(), p, one) != count) bad++;
    }
    s.check("gft: with the trivial group the series counts graphs, tensorial and coloured (failures)", bad, 0);
    const Rational w = pq::gft_free_energy(pq::boulatov(), 2, Z2);
    s.require("gft: Boulatov's model with ℤ₂ at order (λλ̄)², the free energy", !w.is_zero(), w.str());
  }
  return s.done();
}
