// Checks of physis/field_theory/large_n.hpp: the vector model's gap equation against the
// leading free energy of the exact expansion; the BIPZ eigenvalue density (normalised, Catalan
// moments at g = 0, the derivative of the planar free energy its fourth moment), the planar
// free energy's closed form against its series and the series against a direct count of
// planar quartic ribbon graphs, and the string susceptibility γ = −1/2 at g_c = −1/48; the
// melonic two-point function against its Fuss–Catalan series, its critical point, its exponent
// γ = 1/2, and the quartic melonic expansion of invariant.hpp.

#include <physis/field_theory/invariant.hpp>
#include <physis/field_theory/large_n.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <numeric>
#include <string>
#include <vector>

namespace pf = physis::field_theory;
using morphe::arithmetic::Integer;
using morphe::arithmetic::Rational;

namespace {

// Connected planar quartic ribbon graphs with k vertices: pairings of the 4k half-edges (each
// vertex's four in cyclic order), faces the cycles of γπ, planar when F = k + 2.
std::size_t planar_quartic(std::size_t k) {
  const std::size_t H = 4 * k;
  std::vector<std::size_t> pairing(H, SIZE_MAX);
  std::size_t count = 0;
  auto faces = [&] {
    std::vector<char> seen(H, 0);
    std::size_t f = 0;
    for (std::size_t h = 0; h < H; h++) {
      if (seen[h]) continue;
      f++;
      for (std::size_t x = h; !seen[x];) {
        seen[x] = 1;
        const std::size_t y = pairing[x];
        x = (y / 4) * 4 + (y % 4 + 1) % 4;  // γ after π
      }
    }
    return f;
  };
  auto connected = [&] {
    std::vector<std::size_t> parent(k);
    std::iota(parent.begin(), parent.end(), 0);
    auto root = [&](std::size_t x) {
      while (parent[x] != x) x = parent[x] = parent[parent[x]];
      return x;
    };
    for (std::size_t h = 0; h < H; h++) parent[root(h / 4)] = root(pairing[h] / 4);
    std::size_t c = 0;
    for (std::size_t v = 0; v < k; v++) c += root(v) == v;
    return c == 1;
  };
  std::function<void()> go = [&] {
    std::size_t first = 0;
    while (first < H && pairing[first] != SIZE_MAX) first++;
    if (first == H) {
      if (faces() == k + 2 && connected()) count++;
      return;
    }
    for (std::size_t other = first + 1; other < H; other++)
      if (pairing[other] == SIZE_MAX) {
        pairing[first] = other, pairing[other] = first;
        go();
        pairing[first] = pairing[other] = SIZE_MAX;
      }
  };
  go();
  return count;
}

double catalan(int n) { return std::exp(std::lgamma(2.0 * n + 1) - std::lgamma(n + 1.0) - std::lgamma(n + 2.0)); }

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- the vector model --------------------------------------------------------------------------------------------------
  {
    const auto S = pf::series(pf::vector_model(2), 5);
    const auto want = pf::vector_free_energy(5);
    bool ok = true;
    std::string detail;
    for (std::size_t n = 1; n <= 5; n++) {
      const auto [power, c] = pf::leading(S.free_energy[n]);
      ok = ok && power - static_cast<long long>(n) * S.scaling == 1 && c == want[n];
      detail += c.str() + " ";
    }
    s.require("large N: the vector model's leading free energy from the exact expansion is 2^{n−1} Cat_n / n, the gap equation's, n = 1..5", ok, detail);
    double worst = 0;
    for (double l : {-0.12, -0.05, 0.1, 0.7, 3.0}) {
      const double G = pf::vector_two_point(l);
      worst = std::max(worst, std::fabs(G * (1 + 2 * l * G) - 1));
    }
    s.check("large N: the gap equation G(1 + 2λG) = 1 holds at five couplings", worst, 1e-14);
    // Its series G = Σ Cat_m (−2λ)^m, summed at λ = 0.02.
    double sum = 0;
    for (int m = 0; m < 40; m++) sum += catalan(m) * std::pow(-0.04, m);
    s.check("large N: the gap equation's root is the series Σ Cat_m (−2λ)^m (relative)", std::fabs(pf::vector_two_point(0.02) / sum - 1), 1e-14);
  }
  // ---- the Hermitian quartic matrix model ---------------------------------------------------------------------------------
  {
    double worst = 0;
    for (double g : {-1.0 / 48, -0.015, -0.005, 0.0, 0.02, 0.5, 4.0})
      worst = std::max(worst, std::fabs(pf::bipz_moment(0, g) - 1));
    s.check("large N: the BIPZ density is normalised, ∫ ρ = 1, from g = g_c to g = 4", worst, 1e-10);
    double cat = 0;
    for (int j = 0; j <= 5; j++) cat = std::max(cat, std::fabs(pf::bipz_moment(2 * j, 0) - catalan(j)));
    s.check("large N: at g = 0 the moments of the semicircle are the Catalan numbers, to the 10th", cat, 1e-9);
    double der = 0;
    for (double g : {-0.015, 0.01, 0.3}) {
      const double h = 1e-4;
      auto E = [](double x) { return pf::bipz_free_energy(x); };
      const double dE = (8 * (E(g + h) - E(g - h)) - (E(g + 2 * h) - E(g - 2 * h))) / (12 * h);  // fourth order
      der = std::max(der, std::fabs(dE / pf::bipz_moment(4, g) - 1));
    }
    s.check("large N: dE₀/dg = ∫ λ⁴ ρ, the planar ⟨tr M⁴⟩/N (relative)", der, 1e-8);
    const auto e = pf::bipz_series(30);
    double series = 0;
    for (double g : {-0.01, 0.004}) {
      double sum = 0;
      for (std::size_t k = 1; k <= 30; k++) sum += e[k].to_double() * std::pow(g, static_cast<double>(k));
      series = std::max(series, std::fabs(sum / pf::bipz_free_energy(g) - 1));
    }
    s.check("large N: E₀'s closed form equals its series −Σ (−12g)^k (2k−1)!/(k!(k+2)!) (relative, 30 terms)", series, 1e-12);
    bool count = true;
    std::string detail;
    for (std::size_t k = 1; k <= 4; k++) {
      const std::size_t P = planar_quartic(k);
      // E₀ = −Σ (−g)^k / k! P_k.
      const Rational want = Rational(Integer(static_cast<long long>(k % 2 ? 1 : -1)) * Integer(static_cast<unsigned long long>(P)), morphe::arithmetic::factorial(k));
      count = count && want == e[k];
      detail += std::to_string(P) + " ";
    }
    s.require("large N: E₀'s coefficients count the connected planar quartic ribbon graphs, k = 1..4", count, detail);
    // The coefficients times g_c^k, exactly, so that they stay within doubles.
    std::vector<double> c;
    const auto big = pf::bipz_series(400);
    for (std::size_t k = 0; k < big.size(); k++)
      c.push_back(k == 0 ? 1.0 : (big[k] * morphe::arithmetic::pow(Rational(Integer(-1), Integer(48)), static_cast<long long>(k))).to_double());
    const double alpha = pf::exponent_from_series(c, 1.0);
    s.check("large N: the string susceptibility γ = 3 − α = −1/2 at g_c = −1/48, from 400 coefficients", (3 - alpha) - (-0.5), 1e-3);
  }
  // ---- melonic resummation -------------------------------------------------------------------------------------------------
  {
    double worst = 0, crit = 0;
    for (int m : {2, 4, 5}) {
      const auto c = pf::melonic_coefficients(m, 40);
      for (double frac : {-0.9, -0.3, 0.2, 0.5}) {
        const double x = frac * pf::melonic_critical(m);
        double sum = 0;
        for (std::size_t p = 0; p < c.size(); p++) sum += c[p].to_double() * std::pow(x, static_cast<double>(p));
        const double G = pf::melonic_two_point(x, m);
        worst = std::max(worst, std::max(std::fabs(1 + x * std::pow(G, m) - G), frac > -0.5 && frac < 0.6 ? std::fabs(sum / G - 1) : 0.0));
      }
      crit = std::max(crit, std::fabs(pf::melonic_two_point(pf::melonic_critical(m), m) - m / (m - 1.0)));
    }
    s.check("large N: G = 1 + x G^m solved, and equal to its Fuss–Catalan series, m = 2, 4, 5", worst, 1e-12);
    s.check("large N: at the critical point x_c = (m−1)^{m−1}/m^m, G = m/(m−1)", crit, 1e-6);
    double gamma = 0;
    for (int m : {2, 4}) {
      std::vector<double> c;
      for (const auto& x : pf::melonic_coefficients(m, 150)) c.push_back(x.to_double());
      gamma = std::max(gamma, std::fabs((2 - pf::exponent_from_series(c, pf::melonic_critical(m))) - 0.5));
    }
    s.check("large N: melonic graphs are branched polymers, γ = 2 − α = 1/2, m = 2, 4", gamma, 1e-3);
    // The exact expansion of one quartic melonic bubble, G = 1 − 2λG²: w_n = 2^{n−1} Cat_n / n,
    // the coefficients of −∫ G(−x/2)² from the melonic series with m = 2.
    const auto S = pf::series(pf::quartic_melonic_model(4, false), 4);
    const auto cm = pf::melonic_coefficients(2, 5);
    bool ok = true;
    for (std::size_t n = 1; n <= 4; n++) {
      const auto [power, c] = pf::leading(S.free_energy[n]);
      ok = ok && power - static_cast<long long>(n) * S.scaling == 4 && c == Rational(morphe::arithmetic::pow(Integer(2), n - 1) * cm[n], Integer(static_cast<long long>(n)));
    }
    s.require("large N: the rank-4 quartic melonic expansion leads with the melonic series of m = 2, n = 1..4", ok);
  }
  return s.done();
}
