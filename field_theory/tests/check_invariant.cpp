// Checks of physis/field_theory/invariant.hpp: Gaussian moments in every rank (⟨(φ̄φ)^m⟩ = m!,
// ⟨(φ̄·φ)^m⟩ = N(N+1)…(N+m−1), ⟨tr (MM†)^k⟩ = Σ_σ N^{c(σ)+c(γσ⁻¹)} with Catalan leading
// coefficients); the free energy from the logarithm of Z equal to the sum over connected
// graphs; the isomorphism classes (their sizes and automorphisms reproducing Z); the genus
// expansion of the matrix model; Gurau's face formula and the melonic leading order of quartic
// tensor models, Cat_n / n.

#include <physis/field_theory/invariant.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <numeric>
#include <string>
#include <vector>

namespace pf = physis::field_theory;
namespace mt = morphe::topology;
using morphe::arithmetic::Integer;
using morphe::arithmetic::Rational;

namespace {

std::size_t cycles(const std::vector<std::size_t>& p) {
  std::vector<char> seen(p.size(), 0);
  std::size_t c = 0;
  for (std::size_t i = 0; i < p.size(); i++)
    if (!seen[i]) {
      c++;
      for (std::size_t j = i; !seen[j]; j = p[j]) seen[j] = 1;
    }
  return c;
}

Integer catalan(long long n) { return morphe::arithmetic::exact_div(morphe::arithmetic::binomial(2 * n, n), Integer(n + 1)); }

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- moments in each rank -------------------------------------------------------------------------------------------
  {
    bool ok = true;
    for (std::size_t m = 1; m <= 6; m++) {
      const auto E = pf::expand(pf::complex_scalar(m), {1}, false);
      ok = ok && E.all.size() == 1 && E.all.begin()->first == 0 && E.all.begin()->second == morphe::arithmetic::factorial(m);
    }
    s.require("invariant: rank 0, ⟨(φ̄φ)^m⟩ = m!, m = 1..6", ok);
  }
  {
    bool ok = true;
    for (std::size_t m = 1; m <= 7; m++) {
      // N(N+1)…(N+m−1), expanded.
      std::vector<Integer> rising{Integer(0), Integer(1)};
      for (std::size_t j = 1; j < m; j++) {
        std::vector<Integer> next(rising.size() + 1, Integer(0));
        for (std::size_t i = 0; i < rising.size(); i++) {
          next[i + 1] += rising[i];
          next[i] += rising[i] * Integer(static_cast<long long>(j));
        }
        rising = next;
      }
      const auto E = pf::expand(pf::vector_model(m), {1}, false);
      for (std::size_t i = 0; i < rising.size(); i++) {
        const auto it = E.all.find(static_cast<long long>(i));
        const Integer got = it == E.all.end() ? Integer(0) : it->second;
        ok = ok && got == rising[i];
      }
    }
    s.require("invariant: rank 1, ⟨(φ̄·φ)^m⟩ = N(N+1)…(N+m−1), m = 1..7", ok);
  }
  {
    bool ok = true, lead = true;
    for (std::size_t k = 1; k <= 8; k++) {
      // Σ_σ N^{c(σ) + c(γσ⁻¹)}, γ the k-cycle.
      std::map<long long, Integer> want;
      std::vector<std::size_t> sigma(k), inv(k), gs(k);
      std::iota(sigma.begin(), sigma.end(), 0);
      do {
        for (std::size_t i = 0; i < k; i++) inv[sigma[i]] = i;
        for (std::size_t i = 0; i < k; i++) gs[i] = (inv[i] + 1) % k;
        want[static_cast<long long>(cycles(sigma) + cycles(gs))] += 1;
      } while (std::next_permutation(sigma.begin(), sigma.end()));
      const auto E = pf::expand(pf::matrix_model(k), {1}, false);
      ok = ok && E.all == want;
      lead = lead && E.all.rbegin()->first == static_cast<long long>(k) + 1 && E.all.rbegin()->second == catalan(static_cast<long long>(k));
    }
    s.require("invariant: rank 2, ⟨tr (MM†)^k⟩ = Σ_σ N^{c(σ)+c(γσ⁻¹)}, k = 1..8", ok);
    s.require("invariant: rank 2, its leading term Catalan(k) N^{k+1}", lead);
  }
  // ---- the free energy, two ways --------------------------------------------------------------------------------------
  {
    int bad = 0;
    for (const auto& T : {pf::vector_model(2), pf::matrix_model(2), pf::quartic_melonic_model(3, false), pf::complex_scalar(2)}) {
      const auto S = pf::series(T, 5);
      for (std::size_t n = 1; n <= 5; n++)
        if (S.free_energy[n] != S.connected[n]) bad++;
    }
    s.check("invariant: log Z = the sum over connected graphs, to order 5, in four theories (orders that differ)", bad, 0);
  }
  // ---- classes and automorphisms --------------------------------------------------------------------------------------
  {
    int bad = 0;
    std::size_t classes = 0;
    for (const auto& T : {pf::vector_model(2), pf::matrix_model(2), pf::matrix_model(3), pf::quartic_melonic_model(3, false), pf::wheel_model(),
                          pf::complex_scalar(2)}) {
      const auto S = pf::series(T, std::min<std::size_t>(4, 8 / T.interactions[0].bubble.p));
      const long long autB = static_cast<long long>(mt::bubble_automorphisms(T.interactions[0].bubble).size());
      for (std::size_t n = 1; n <= 4 && n * T.interactions[0].bubble.p <= 8; n++) {
        const auto E = pf::expand(T, {n});
        std::size_t total = 0;
        pf::Polynomial z, w;
        for (const auto& C : E.classes) {
          total += C.labelled;
          if (C.automorphisms * Integer(static_cast<unsigned long long>(C.labelled)) != E.relabellings) bad++;
          const Rational weight(morphe::arithmetic::pow(Integer(autB), n), C.automorphisms);
          z[C.faces] += weight;
          if (C.connected) w[C.faces] += weight;
        }
        classes += E.classes.size();
        std::size_t P = n * T.interactions[0].bubble.p, fact = 1;
        for (std::size_t k = 2; k <= P; k++) fact *= k;
        if (total != fact) bad++;
        for (auto it = w.begin(); it != w.end();) it = it->second.is_zero() ? w.erase(it) : std::next(it);
        if (z != S.z[n] || w != S.free_energy[n]) bad++;
      }
    }
    s.check("invariant: classes partition the matchings, |Aut Γ| · |class| = |H|, Σ |Aut B|^n N^F / |Aut Γ| = z_n and w_n (" + std::to_string(classes) + " classes; failures)",
            bad, 0);
    // The vector model at first order: the two graphs N² and N, each with |Aut| = 2.
    const auto E = pf::expand(pf::vector_model(2), {1});
    s.require("invariant: (φ̄·φ)² at first order: two graphs, N² and N, |Aut| = 2 each",
              E.classes.size() == 2 && E.classes[0].automorphisms == Integer(2) && E.classes[1].automorphisms == Integer(2) &&
                  E.all == std::map<long long, Integer>{{1, Integer(1)}, {2, Integer(1)}});
  }
  // ---- the genus expansion of the matrix model -------------------------------------------------------------------------
  {
    int bad = 0;
    std::size_t planar = 0;
    for (std::size_t n = 1; n <= 4; n++) {
      const auto E = pf::expand(pf::matrix_model(2), {n});
      for (const auto& C : E.classes) {
        if (!C.connected) continue;
        if (C.scaled != 2 - 2 * static_cast<long long>(C.degree)) bad++;
        planar += C.degree == 0;
      }
    }
    s.check("invariant: tr (MM†)² with λ/N: every connected graph carries N^{2 − 2g}, g its genus (" + std::to_string(planar) + " planar classes; failures)", bad, 0);
  }
  // ---- tensors ---------------------------------------------------------------------------------------------------------
  {
    int bad = 0;
    for (std::size_t r : {3, 4})
      for (const std::vector<std::size_t>& counts : std::vector<std::vector<std::size_t>>{{1, 1, 0}, {2, 0, 0}, {1, 1, 1}, {2, 1, 0}}) {
        auto T = pf::quartic_melonic_model(r, true);
        std::vector<std::size_t> c(r, 0);
        for (std::size_t i = 0; i < 3; i++) c[i] = counts[i];
        const auto E = pf::expand(T, c);
        for (const auto& C : E.classes) {
          if (!C.connected) continue;
          // Gurau: F_total = r(r−1)P/2 + r − 2ω/(r−1)!, and the scaled power r − 2ω/(r−1)!.
          long long Ftot = 0;
          for (std::size_t i = 0; i <= r; i++)
            for (std::size_t j = i + 1; j <= r; j++) Ftot += static_cast<long long>(C.graph.faces(i, j).size());
          long long fact = 1;
          for (std::size_t k = 2; k + 1 <= r; k++) fact *= static_cast<long long>(k);
          const long long P = static_cast<long long>(C.graph.p), R = static_cast<long long>(r), w = static_cast<long long>(C.degree);
          if (Ftot * fact != (R * (R - 1) * P / 2 + R) * fact - 2 * w) bad++;
          if (C.scaled * fact != R * fact - 2 * w) bad++;
        }
      }
    s.check("invariant: quartic melonic models of rank 3 and 4: Gurau's face formula and N^{r − 2ω/(r−1)!} for every connected graph (failures)", bad, 0);
    // The melonic leading order of one quartic bubble: G = 1 − 2λG² (either black vertex of the
    // bubble may be the external one), so ∂_λ W = −N^r G² gives w_n = 2^{n−1} Cat_n / n at N^r.
    bool ok = true;
    std::string detail;
    const auto S = pf::series(pf::quartic_melonic_model(3, false), 5);
    for (std::size_t n = 1; n <= 5; n++) {
      const auto [power, coefficient] = pf::leading(S.free_energy[n]);
      const Rational want(morphe::arithmetic::pow(Integer(2), n - 1) * catalan(static_cast<long long>(n)), Integer(static_cast<long long>(n)));
      if (power - static_cast<long long>(n) * S.scaling != 3 || coefficient != want) ok = false;
      detail += coefficient.str() + " ";
    }
    s.require("invariant: one quartic melonic bubble in rank 3, the leading free energy 2^{n−1} Cat_n / n N³, n = 1..5", ok, detail);
    // One bubble of each of the D colours, each with the coupling λ: G = 1 − 2Dλ G², so the
    // leading coefficient of (−λ)^n is D (2D)^{n−1} Cat_n / n. From the expansion over every
    // distribution of the n bubbles among the colours, weighted 1 / Π n_c!.
    ok = true;
    detail.clear();
    for (std::size_t n = 1; n <= 3; n++) {
      Rational lead = 0;
      std::vector<std::size_t> c(3, 0);
      for (c[0] = 0; c[0] <= n; c[0]++)
        for (c[1] = 0; c[0] + c[1] <= n; c[1]++) {
          c[2] = n - c[0] - c[1];
          const auto E = pf::expand(pf::quartic_melonic_model(3, true), c, false);
          const auto it = E.connected.find(3 + E.scaling);
          if (it == E.connected.end()) continue;
          lead += Rational(it->second, morphe::arithmetic::factorial(c[0]) * morphe::arithmetic::factorial(c[1]) * morphe::arithmetic::factorial(c[2]));
        }
      const Rational want(Integer(3) * morphe::arithmetic::pow(Integer(6), n - 1) * catalan(static_cast<long long>(n)), Integer(static_cast<long long>(n)));
      ok = ok && lead == want;
      detail += lead.str() + " ";
    }
    s.require("invariant: the three quartic melonic bubbles of rank 3, leading free energy D (2D)^{n−1} Cat_n / n, n = 1..3", ok, detail);
    // Bubbles that are not melonic: with the melonic scaling a connected graph carries at most
    // N^{r + Σ_B 2ω_B/(r−2)!}, ω_B the bubble's own degree (Gurau's bound, the bubbles' faces
    // removed); at first order the wheel K₃,₃ (ω = 1) and the necklace.
    std::string powers;
    bool bound = true;
    for (const auto& T : {pf::wheel_model(), pf::necklace_model()}) {
      const long long wB = static_cast<long long>(T.interactions[0].bubble.degree());
      const auto E = pf::expand(T, {1});
      long long top = -1000;
      for (const auto& C : E.classes)
        if (C.connected) top = std::max(top, C.scaled), bound = bound && C.scaled <= 3 + 2 * wB;
      powers += T.interactions[0].name + ": ω_B = " + std::to_string(wB) + ", leading N^" + std::to_string(top) + "; ";
    }
    s.require("invariant: sextic bubbles at first order obey N^{r + 2ω_B/(r−2)!}", bound, powers);
  }
  return s.done();
}
