// Checks of physis/quantum_gravity/stranded.hpp: Boulatov's tetrahedral vertex on real fields over
// finite groups. Its values on the constant field and on (−1)^a over ℤ₂ (64 and 0) and along every
// character of ℤ₂³; its invariance under G × G acting on the slot classes {1, 3}, {2} (for ℤ₂, ℤ₃
// and the non-abelian S₃) and the failure of invariance outside that subgroup; the gauge projection
// annihilating the functional odd under the simultaneous translation, for the abelian ℤ₂. Wick's
// theorem: Σ over pairings of |G|^F equal to the Gaussian moment ⟨V^n⟩ computed by brute force, for
// n = 1, 2 and G = 1, ℤ₂, ℤ₃, S₃, and n = 3 for G = 1, ℤ₂; the classes of pairings, each holding
// |Aut v|^n n!/|Aut Γ| of them; the vertex's automorphisms, the cyclic group of order 4; and the 24
// ways of joining two copies factor to factor: 4 give the two-vertex complex with six faces, Betti
// numbers (1, 0, 3) and the BF amplitude |G|^{F−E+V−1} of the 3-sphere, 4 give four faces, 16 three.

#include <physis/quantum_gravity/stranded.hpp>
#include <calx/verification/check.hpp>

#include <map>
#include <random>
#include <string>
#include <vector>

namespace pq = physis::quantum_gravity;
namespace ma = morphe::algebra;
namespace mt = morphe::topology;
using morphe::arithmetic::Integer;
using morphe::arithmetic::Rational;

int main() {
  calx::verification::Suite s;
  const ma::FiniteGroup one = ma::cyclic(1), Z2 = ma::cyclic(2), Z3 = ma::cyclic(3), S3 = ma::symmetric(3);
  const auto tet = pq::tetrahedral_vertex();

  // ---- values ----------------------------------------------------------------------------------------------------------------
  {
    std::vector<Rational> constant(8, Rational(1)), odd(8);
    for (std::size_t p = 0; p < 8; p++) odd[p] = pq::coordinates(p, 2, 3)[0] ? Rational(-1) : Rational(1);
    s.require("V_tet[1] = 64 and V_tet[(−1)^a] = 0 over ℤ₂", pq::interaction(tet, Z2, constant) == Rational(64) && pq::interaction(tet, Z2, odd) == Rational(0));
    // Along each character (−1)^{σ·x} of ℤ₂³: 64 for σ ∈ {000, 010, 101, 111}, 0 for the others.
    std::string got;
    for (std::size_t sigma = 0; sigma < 8; sigma++) {
      const auto sv = pq::coordinates(sigma, 2, 3);
      std::vector<Rational> chi(8);
      for (std::size_t p = 0; p < 8; p++) {
        const auto x = pq::coordinates(p, 2, 3);
        chi[p] = (sv[0] * x[0] + sv[1] * x[1] + sv[2] * x[2]) % 2 ? Rational(-1) : Rational(1);
      }
      got += pq::interaction(tet, Z2, chi) == Rational(64) ? "1" : pq::interaction(tet, Z2, chi) == Rational(0) ? "0" : "?";
    }
    s.require("V_tet along the characters of ℤ₂³ (σ = 000 … 111)", got == "10100101", got);
  }

  // ---- symmetries ------------------------------------------------------------------------------------------------------------
  {
    const auto classes = pq::slot_classes(tet);
    s.require("the slot classes of the tetrahedral vertex are {1, 3} and {2}",
              classes == std::vector<std::vector<std::size_t>>{{0, 2}, {1}});
    std::mt19937 rng(11);
    bool invariant = true, broken = false;
    for (const auto* G : {&Z2, &Z3, &S3}) {
      const std::size_t q = G->order();
      std::vector<Rational> psi(q * q * q);
      for (auto& v : psi) v = Rational(static_cast<long long>(rng() % 7) - 3);
      const Rational v0 = pq::interaction(tet, *G, psi);
      for (std::uint32_t g = 0; g < q; g++)
        for (std::uint32_t h = 0; h < q; h++) {
          invariant = invariant && pq::interaction(tet, *G, pq::translate(psi, *G, {g, h, g})) == v0;
          if (g != h) broken = broken || pq::interaction(tet, *G, pq::translate(psi, *G, {g, 0, h})) != v0;
        }
    }
    s.require("V_tet is invariant under (g, h, g) for ℤ₂, ℤ₃ and S₃", invariant);
    s.require("and not under translations that move slots 1 and 3 apart", broken);
    std::vector<Rational> psi(8);
    for (auto& v : psi) v = Rational(static_cast<long long>(rng() % 9) - 4);
    const auto eta = [](const std::vector<Rational>& f) {
      Rational e = 0;
      for (std::size_t p = 0; p < 8; p++) e = e + (pq::coordinates(p, 2, 3)[0] ? Rational(-1) : Rational(1)) * f[p];
      return e / Rational(8);
    };
    s.require("over ℤ₂ the gauge projection annihilates η[ψ] = 2⁻³ Σ (−1)^a ψ", eta(psi) != Rational(0) && eta(pq::gauge_project(psi, Z2, 3)) == Rational(0),
              "η[ψ] = " + eta(psi).str());
  }

  // ---- Wick's theorem --------------------------------------------------------------------------------------------------------
  {
    for (const auto* G : {&one, &Z2, &Z3, &S3})
      for (std::size_t n : {1, 2}) {
        if (G == &S3 && n == 2) continue;  // 6^12 assignments: too many for brute force here
        Integer wick = 0;
        pq::for_each_pairing(4 * n, [&](const pq::Pairing& p) {
          wick += morphe::arithmetic::pow(Integer(static_cast<long long>(G->order())), pq::stranded_graph(tet, n, p).faces());
        });
        const Integer brute = pq::gaussian_moment(tet, *G, n);
        s.require("⟨V_tet^" + std::to_string(n) + "⟩ over " + G->name + ": Wick pairings = brute force", wick == brute, wick.str());
      }
    for (const auto* G : {&one, &Z2}) {
      Integer wick = 0;
      pq::for_each_pairing(12, [&](const pq::Pairing& p) {
        wick += morphe::arithmetic::pow(Integer(static_cast<long long>(G->order())), pq::stranded_graph(tet, 3, p).faces());
      });
      s.require("⟨V_tet^3⟩ over " + G->name + ": Wick pairings = brute force", wick == pq::gaussian_moment(tet, *G, 3), wick.str());
    }
    const auto mel = pq::melonic_vertex();
    Integer wick = 0;
    pq::for_each_pairing(8, [&](const pq::Pairing& p) { wick += morphe::arithmetic::pow(Integer(3), pq::stranded_graph(mel, 2, p).faces()); });
    s.require("⟨V_melonic^2⟩ over ℤ₃: Wick pairings = brute force", wick == pq::gaussian_moment(mel, Z3, 2), wick.str());
  }

  // ---- classes and automorphisms ---------------------------------------------------------------------------------------------
  {
    const auto aut = pq::vertex_automorphisms(tet);
    s.require("Aut of the tetrahedral vertex: the 4 cyclic shifts of its factors", aut.size() == 4);
    for (std::size_t n : {1, 2, 3}) {
      const auto classes = pq::feynman_classes(tet, n);
      std::size_t total = 0, wreath = 1;
      for (std::size_t i = 0; i < n; i++) wreath *= 4 * (i + 1);
      bool ok = true;
      for (const auto& c : classes) {
        total += c.pairings;
        ok = ok && c.pairings * c.automorphisms == wreath;
      }
      std::size_t pairings = 1;
      for (std::size_t k = 4 * n - 1; k > 1; k -= 2) pairings *= k;
      s.require("n = " + std::to_string(n) + ": " + std::to_string(classes.size()) + " classes, each of |Aut v|^n n!/|Aut Γ| pairings",
                ok && total == pairings, std::to_string(total) + " pairings");
    }
  }

  // ---- two copies joined factor to factor -------------------------------------------------------------------------------------
  {
    std::map<std::size_t, std::size_t> by_faces;
    bool dipole_ok = false;
    std::vector<std::size_t> sigma{0, 1, 2, 3};
    do {
      pq::Pairing p(8);
      for (std::size_t i = 0; i < 4; i++) p[i] = 4 + sigma[i], p[4 + sigma[i]] = i;
      const auto S = pq::stranded_graph(tet, 2, p);
      by_faces[S.faces()]++;
      if (sigma == std::vector<std::size_t>{0, 1, 2, 3}) {
        const auto b = mt::betti(S.complex);
        bool lens = true;  // every face crosses two propagators, one each way: ∂f = ±(e_i − e_j)
        for (const auto& f : S.complex.faces) lens = lens && f.boundary.size() == 2 && f.boundary[0].sign == -f.boundary[1].sign;
        dipole_ok = S.complex.counts() == std::vector<std::size_t>{2, 4, 6} && b == std::vector<std::size_t>{1, 0, 3} && lens &&
                    pq::gft_amplitude(S.complex, Z2, true) == Rational(8) && pq::gft_amplitude(S.complex, Z2, false) == Rational(64) &&
                    pq::gft_amplitude(S.complex, S3, true) == Rational(216);
      }
    } while (std::next_permutation(sigma.begin(), sigma.end()));
    s.require("the identity gluing: (V, E, F) = (2, 4, 6), Betti (1, 0, 3), lens faces, BF amplitude 8 over ℤ₂ and 216 over S₃", dipole_ok);
    s.require("of the 24 gluings, 4 have six faces, 4 four, 16 three", by_faces == std::map<std::size_t, std::size_t>{{3, 16}, {4, 4}, {6, 4}});
  }
  return s.done();
}
