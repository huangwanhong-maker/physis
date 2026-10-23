// Checks of physis/field_theory: the free propagator on lattice tori (dense inverse against
// the momentum sum) and on an irregular complex (the defining equation); the heat bath and
// Metropolis against the propagator within their errors; φ⁴ on one point against the
// zero-dimensional integral; Binder's cumulant of a Gaussian; the autocorrelation time of an
// AR(1) series against its closed form; Feynman graphs (Wick pairings grouped correctly, the
// sum of symmetry factors, the connected graphs against the exact log Z series) and the
// optimal truncation of the asymptotic series.

#include <physis/field_theory/feynman.hpp>
#include <physis/field_theory/scalar.hpp>
#include <physis/field_theory/statistics.hpp>
#include <physis/field_theory/zero_dimensional.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <morphe/topology/constructions.hpp>

#include <cmath>
#include <map>
#include <string>
#include <vector>

namespace ft = physis::field_theory;
namespace mt = morphe::topology;
using morphe::arithmetic::Integer;
using morphe::arithmetic::Rational;

int main() {
  calx::verification::Suite s;

  // ---- the free propagator ------------------------------------------------------------------------------------------
  {
    double worst = 0;
    for (const std::vector<std::size_t>& sides : {std::vector<std::size_t>{4, 4}, std::vector<std::size_t>{3, 4, 5}, std::vector<std::size_t>{2, 6}}) {
      const auto X = mt::torus_lattice(sides);
      const double m2 = 0.3;
      const auto G = ft::free_propagator(X, m2);
      for (std::size_t y = 0; y < X.vertices; y++) {
        const auto cx = ft::lattice_coordinates(0, sides), cy = ft::lattice_coordinates(y, sides);
        std::vector<long long> r(sides.size());
        for (std::size_t mu = 0; mu < sides.size(); mu++) r[mu] = static_cast<long long>(cx[mu]) - static_cast<long long>(cy[mu]);
        worst = std::max(worst, std::fabs(G(0, y) - ft::lattice_propagator(sides, m2, r)));
      }
    }
    s.check("field_theory: (Δ + m²)⁻¹ on 4×4, 3×4×5 and 2×6 lattice tori against the momentum sum (largest difference)", worst, 1e-13);
    const auto B = mt::barycentric_subdivision(mt::torus_7()).cw();
    const auto G = ft::free_propagator(B, 0.7);
    auto A = ft::laplacian(B);
    for (std::size_t i = 0; i < A.n; i++) A(i, i) += 0.7;
    const auto I = A * G;
    double res = 0, asym = 0;
    for (std::size_t i = 0; i < I.n; i++)
      for (std::size_t j = 0; j < I.n; j++) {
        res = std::max(res, std::fabs(I(i, j) - (i == j ? 1.0 : 0.0)));
        asym = std::max(asym, std::fabs(G(i, j) - G(j, i)));
      }
    s.check("field_theory: on the subdivided 7-vertex torus (42 vertices, degrees 4 to 12), (Δ + m²) G = I", res, 1e-12);
    s.check("field_theory: ... and G is symmetric", asym, 1e-13);
  }
  // ---- Monte Carlo of the free field ------------------------------------------------------------------------------------
  {
    const std::vector<std::size_t> sides{4, 4};
    const auto X = mt::torus_lattice(sides);
    const auto N = ft::neighbours(X);
    const ft::ScalarTheory T{0.3, 0};
    const auto G = ft::free_propagator(X, T.m2);
    calx::random::Stream z(11, 1, 0, 0);
    std::vector<double> phi(X.vertices, 0.0);
    for (int i = 0; i < 1000; i++) ft::heat_bath_sweep(T, N, phi, z);
    const std::size_t far = 2 * 4 + 2;  // the vertex at (2, 2)
    std::vector<double> c00, c01, c22, m2s, m4s;
    for (int i = 0; i < 40000; i++) {
      ft::heat_bath_sweep(T, N, phi, z);
      c00.push_back(phi[0] * phi[0]);
      c01.push_back(phi[0] * phi[1]);
      c22.push_back(phi[0] * phi[far]);
      const double M = ft::magnetisation(phi);
      m2s.push_back(M * M);
      m4s.push_back(M * M * M * M);
    }
    double worst = 0;
    for (auto [series, y] : {std::pair{&c00, std::size_t{0}}, std::pair{&c01, std::size_t{1}}, std::pair{&c22, far}}) {
      const auto e = ft::binned(*series, 50);
      worst = std::max(worst, std::fabs(e.mean - G(0, y)) / e.error);
    }
    s.check("field_theory: heat bath on the 4×4 torus: ⟨φ_0 φ_y⟩ at distance 0, 1, (2,2) against G, in binned standard errors", worst, 4.0, "σ");
    const double chi_exact = [&] {
      double t = 0;
      for (std::size_t i = 0; i < G.n; i++)
        for (std::size_t j = 0; j < G.n; j++) t += G(i, j);
      return t / static_cast<double>(G.n);
    }();
    const auto chi = ft::binned(m2s, 50);
    s.check("field_theory: ... the susceptibility ⟨M²⟩/V against Σ G / V = 1/m², in standard errors",
            std::fabs(chi.mean / 16 - chi_exact) / (chi.error / 16), 4.0, "σ");
    const auto U = ft::jackknife({m2s, m4s}, [](const std::vector<double>& m) { return ft::binder_cumulant(m[0], m[1]); }, 50);
    s.check("field_theory: ... Binder's cumulant of the Gaussian field is 0, in jackknife errors", std::fabs(U.mean) / U.error, 4.0, "σ");
    // Metropolis, the same field.
    std::vector<double> psi(X.vertices, 0.0), sq;
    std::size_t acc = 0;
    for (int i = 0; i < 1000; i++) ft::metropolis_sweep(T, N, psi, z, 1.5);
    for (int i = 0; i < 60000; i++) {
      acc += ft::metropolis_sweep(T, N, psi, z, 1.5);
      sq.push_back(ft::mean_square(psi));
    }
    const auto e = ft::estimate(sq);
    s.check("field_theory: Metropolis on the 4×4 torus: ⟨φ²⟩ against G(0,0), in standard errors from τ_int", std::fabs(e.mean - G(0, 0)) / e.error, 4.0, "σ");
    s.require("field_theory: ... with an acceptance between 0.2 and 0.8", acc > 0.2 * 60000 * 16 && acc < 0.8 * 60000 * 16,
              std::to_string(static_cast<double>(acc) / (60000.0 * 16)));
  }
  // ---- φ⁴ on one point --------------------------------------------------------------------------------------------------
  {
    mt::CWComplex point;
    point.vertices = 1;
    const auto N = ft::neighbours(point);
    const ft::ScalarTheory T{1.0, 2.0};
    calx::random::Stream z(12, 1, 0, 0);
    std::vector<double> phi{0.0}, p2, p4;
    for (int i = 0; i < 400000; i++) {
      ft::metropolis_sweep(T, N, phi, z, 2.0);
      p2.push_back(phi[0] * phi[0]);
      p4.push_back(phi[0] * phi[0] * phi[0] * phi[0]);
    }
    const auto e2 = ft::estimate(p2), e4 = ft::estimate(p4);
    const double w2 = ft::zero_dimensional_moment(2, 2.0), w4 = ft::zero_dimensional_moment(4, 2.0);
    s.check("field_theory: φ⁴ on one point (m² = 1, λ = 2): Metropolis ⟨φ²⟩ against the integral, in standard errors", std::fabs(e2.mean - w2) / e2.error, 4.0, "σ");
    s.check("field_theory: ... and ⟨φ⁴⟩", std::fabs(e4.mean - w4) / e4.error, 4.0, "σ");
  }
  // ---- statistics ------------------------------------------------------------------------------------------------------
  {
    // AR(1): x' = ρ x + √(1-ρ²) ε has τ_int = (1 + ρ) / (2 (1 - ρ)).
    calx::random::Stream z(13, 1, 0, 0);
    const double rho = 0.9;
    std::vector<double> x(400000);
    double v = 0;
    for (auto& y : x) y = v = rho * v + std::sqrt(1 - rho * rho) * z.normal();
    const double tau = ft::autocorrelation_time(x), want = (1 + rho) / (2 * (1 - rho));
    s.check("field_theory: the integrated autocorrelation time of AR(1), ρ = 0.9, against 9.5 (relative)", tau / want - 1, 0.1);
    // 200 bins of 2000 ≫ τ: the binned error is itself uncertain by 1/√(2·199) = 5 %.
    const auto a = ft::estimate(x), b = ft::binned(x, 200);
    s.check("field_theory: binned and windowed errors of the mean agree (relative)", b.error / a.error - 1, 0.2);
  }
  // ---- Feynman graphs ------------------------------------------------------------------------------------------------------
  {
    // Wick pairings grouped into labelled multigraphs, against the closed-form multiplicity.
    int bad = 0;
    for (auto [k, n] : {std::pair<std::size_t, std::size_t>{4, 2}, {4, 3}, {3, 4}, {6, 2}}) {
      std::map<ft::Adjacency, Integer> seen;
      ft::wick_pairings(k, n, [&](const std::vector<std::size_t>& p) { seen[ft::graph_of(p, k, n)] += Integer(1); });
      for (const auto& [a, c] : seen)
        if (ft::detail::pairings_of(a, k) != c) bad++;
      std::size_t labelled = 0;
      ft::labelled_vacuum_graphs(k, n, [&](const ft::Adjacency&) { labelled++; });
      if (labelled != seen.size()) bad++;
    }
    s.require("field_theory: Wick pairings of φ⁴ (orders 2, 3), φ³ (order 4), φ⁶ (order 2) produce each multigraph Π k! / (Π 2^a a! Π a!) times",
              bad == 0);
    const auto g1 = ft::vacuum_graphs(4, 1), g2 = ft::vacuum_graphs(4, 2);
    std::vector<std::string> w2;
    for (const auto& g : g2) w2.push_back(g.weight().str());
    std::sort(w2.begin(), w2.end());
    s.require("field_theory: φ⁴ vacuum graphs: the figure eight 1/8; at order 2, 1/16, 1/48 and 1/128",
              g1.size() == 1 && g1[0].automorphisms == Integer(8) && w2 == std::vector<std::string>{"1/128", "1/16", "1/48"});
    bool sums = true;
    for (auto [k, n] : {std::pair<std::size_t, std::size_t>{4, 0}, {4, 1}, {4, 2}, {4, 3}, {4, 4}, {4, 5}, {3, 2}, {3, 4}, {3, 6}, {6, 3}}) {
      Rational t = 0;
      for (const auto& g : ft::vacuum_graphs(k, n)) t += g.weight();
      const Rational want(morphe::arithmetic::double_factorial(static_cast<long long>(k * n) - 1),
                          morphe::arithmetic::pow(morphe::arithmetic::factorial(k), n) * morphe::arithmetic::factorial(n));
      sums = sums && t == want;
    }
    s.require("field_theory: Σ_Γ 1/|Aut Γ| = (kn - 1)!! / ((k!)^n n!) for φ³, φ⁴ and φ⁶ up to orders 6, 5 and 3", sums);
    const auto W = ft::log_series(ft::partition_series(5));
    bool connected = true;
    std::vector<std::size_t> counts;
    for (std::size_t n = 1; n <= 5; n++) {
      Rational t = 0;
      std::size_t c = 0;
      for (const auto& g : ft::vacuum_graphs(4, n))
        if (g.connected) t += g.weight(), c++;
      counts.push_back(c);
      if ((n % 2 ? -t : t) != W[n]) connected = false;
    }
    s.require("field_theory: the connected φ⁴ graphs sum to the coefficients of log Z (-1/8, 1/12, -11/96, ...) to order 5", connected,
              W[1].str() + ", " + W[2].str() + ", " + W[3].str() + ", " + W[4].str() + ", " + W[5].str());
    s.require("field_theory: connected φ⁴ vacuum graphs number 1, 2, 4, 10, 28 at orders 1 to 5",
              counts == std::vector<std::size_t>{1, 2, 4, 10, 28});
  }
  // ---- the asymptotic series --------------------------------------------------------------------------------------------------
  {
    s.check("field_theory: Z(0) = 1 by quadrature", ft::zero_dimensional_partition(0) - 1, 1e-14);
    const auto z = ft::partition_series(60);
    double worst = 0;
    for (double lambda : {0.05, 0.2, 0.5}) {
      const double Z = ft::zero_dimensional_partition(lambda);
      for (std::size_t N = 0; N < 50; N++) {
        const double err = std::fabs(Z - ft::partial_sum(z, N, lambda));
        const double next = std::fabs(z[N + 1].to_double()) * std::pow(lambda, static_cast<double>(N + 1));
        worst = std::max(worst, err / (next + 1e-14));
      }
    }
    s.check("field_theory: the truncated series of Z(λ), λ = 0.05, 0.2, 0.5, errs by less than its next term (the largest ratio)", worst, 1.0);
    // Optimal truncation near n = 3 / (2λ) leaves an error of order e^{-3/(2λ)}.
    const double lambda = 0.1, Z = ft::zero_dimensional_partition(lambda);
    double best = 1;
    for (std::size_t N = 0; N < 50; N++) best = std::min(best, std::fabs(Z - ft::partial_sum(z, N, lambda)));
    s.check("field_theory: optimal truncation at λ = 0.1 errs by about e^{-15} (the ratio to e^{-3/(2λ)})", best / std::exp(-1.5 / lambda), 1.0);
  }
  return s.done();
}
