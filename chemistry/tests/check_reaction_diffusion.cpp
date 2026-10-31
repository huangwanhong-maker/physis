// Checks of reaction–diffusion on graphs (reaction_diffusion.hpp) against closed forms and direct
// integration.
//
// - The path Laplacian: eigenvalues (2 − 2cos(πk/n))/spacing², rows summing to zero.
// - Two species against the Schnakenberg kinetics (Schnakenberg 1979; Murray 2003, ch. 2), u_t = a − u + u²v + Δu,
//   v_t = b − u²v + dΔv at a = 0.1, b = 0.9, whose Turing condition in a, b, d is solved here in closed
//   form, d_c = x², (b − a)x² − 2(a + b)²x − (a + b)³ = 0: the general formula's critical ratio agrees,
//   instability just above it and not just below; with the species exchanged, the ratio is 1/d_c; the
//   band's ends are the roots of h(k²); the fastest k² is the largest growth on a scan of the band.
// - turing_modes on the path Laplacian (54 vertices, spacing ½, d = 9): the growth of each mode equals
//   the two-species dispersion relation at its Λ; the unstable modes are exactly those with Λ in the
//   band, here mode 5 alone; it is the fastest.
// - Schnakenberg's kinetics as a reaction network on that chain, integrated by RK4 from a perturbation
//   of 10⁻³: the steady pattern's dominant cosine mode is 5, with three peaks, as predicted.
// - ReactionDiffusion's rate equals mass_action() at each vertex less a dense L D c; A + B ⇌ C on a
//   random graph keeps the totals of A + C and B + C to rounding over 2000 steps.
// - McKean's caricature, integrated from a step for D = ¼, 1, 3: at a = ½ the front comes to rest at the
//   closed form, its two interface sites at a_+ and a_− (so the closed form of the pinning interval is
//   read off an integration that does not know it); just inside the interval no site ever flips; just
//   outside, the state 1 (a < a_−) or 0 (a > a_+) invades.

#include <algorithm>
#include <cmath>
#include <numbers>
#include <string>
#include <vector>

#include <calx/linalg/eigen.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <physis/chemistry/reaction_diffusion.hpp>

using calx::linalg::Matrix;
using calx::verification::Suite;
namespace pc = physis::chemistry;
using morphe::arithmetic::Integer;
using morphe::arithmetic::Rational;

namespace {

// The Schnakenberg kinetics at (u₀, v₀) = (a + b, b/(a + b)²): its Jacobian.
Matrix schnakenberg_jacobian(double a, double b) {
  Matrix J(2);
  J(0, 0) = (b - a) / (a + b), J(0, 1) = (a + b) * (a + b);
  J(1, 0) = -2 * b / (a + b), J(1, 1) = -(a + b) * (a + b);
  return J;
}

std::vector<double> integrate_mckean(double D, double a, std::size_t n, double T, double h, std::vector<std::size_t>* highs = nullptr) {
  std::vector<double> u(n);
  for (std::size_t i = 0; i < n; i++) u[i] = i < n / 2 ? 1.0 : 0.0;
  const auto f = [&](std::span<const double> y) { return pc::mckean_rate(y, D, a); };
  const long steps = std::lround(T / h);
  for (long k = 0; k < steps; k++) {
    u = pc::runge_kutta_4(f, u, h);
    if (highs && k % 100 == 0) highs->push_back(static_cast<std::size_t>(std::count_if(u.begin(), u.end(), [&](double x) { return x >= a; })));
  }
  return u;
}

}  // namespace

int main() {
  Suite s;
  calx::random::Stream rng(20261030, 21, 0, 0);

  // ---- the path Laplacian ---------------------------------------------------------------------------------------------------------
  {
    const std::size_t n = 17;
    const double dx = 0.3;
    const Matrix L = pc::path_laplacian(n, dx);
    const auto ev = calx::linalg::symmetric_eigen(L).values;
    double worst = 0, rows = 0;
    for (std::size_t k = 0; k < n; k++) worst = std::max(worst, std::fabs(ev[k] - (2 - 2 * std::cos(std::numbers::pi * k / n)) / (dx * dx)));
    for (std::size_t i = 0; i < n; i++) {
      double r = 0;
      for (std::size_t j = 0; j < n; j++) r += L(i, j);
      rows = std::max(rows, std::fabs(r));
    }
    s.check("path Laplacian: eigenvalues (2 − 2cos(πk/n))/spacing² (relative to the largest), rows summing to zero", std::max(worst * dx * dx / 4, rows), 1e-13);
  }

  // ---- two species in closed form -------------------------------------------------------------------------------------------------
  const double a = 0.1, b = 0.9;
  const Matrix J = schnakenberg_jacobian(a, b);
  {
    const double x = ((a + b) * (a + b) + std::sqrt(std::pow(a + b, 4) + (b - a) * std::pow(a + b, 3))) / (b - a);
    const double dc = x * x;
    const auto T = pc::turing_two_species(J(0, 0), J(0, 1), J(1, 0), J(1, 1), 1, 20);
    s.check("Schnakenberg (a = 0.1, b = 0.9): the critical ratio against the closed form in a, b, d (relative)", std::fabs(T.critical_ratio / dc - 1), 1e-13);
    const auto above = pc::turing_two_species(J(0, 0), J(0, 1), J(1, 0), J(1, 1), 1, dc * (1 + 1e-6));
    const auto below = pc::turing_two_species(J(0, 0), J(0, 1), J(1, 0), J(1, 1), 1, dc * (1 - 1e-6));
    s.require("Schnakenberg: stable without diffusion; Turing-unstable just above d_c and not just below",
              T.stable_homogeneous && T.unstable && above.unstable && !below.unstable && below.stable_homogeneous);
    const auto mirror = pc::turing_two_species(J(1, 1), J(1, 0), J(0, 1), J(0, 0), 20, 1);
    s.check("the species exchanged (v activates): the critical ratio is 1/d_c, and still unstable", std::fabs(mirror.critical_ratio * dc - 1) + (mirror.unstable ? 0 : 1), 1e-13);

    const double P = 20, B = 20 * J(0, 0) + J(1, 1), det = J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);
    const auto h = [&](double k2) { return P * k2 * k2 - B * k2 + det; };
    const auto growth = [&](double k2) { return pc::turing_two_species_growth(J(0, 0), J(0, 1), J(1, 0), J(1, 1), 1, 20, k2); };
    s.check("d = 20: the band's ends are the roots of h(k²) = P k⁴ − B k² + Δ (relative to Δ)", std::max(std::fabs(h(T.k2_low)), std::fabs(h(T.k2_high))) / det, 1e-13);
    s.require("d = 20: growth inside the band, decay just outside", growth(0.5 * (T.k2_low + T.k2_high)) > 0 && growth(0.99 * T.k2_low) < 0 && growth(1.01 * T.k2_high) < 0);
    double best = -1, arg = 0;
    const int M = 200000;
    for (int i = 1; i < M; i++) {
      const double k2 = T.k2_low + (T.k2_high - T.k2_low) * i / M;
      if (growth(k2) > best) best = growth(k2), arg = k2;
    }
    s.check("d = 20: no k² of a scan of the band (2 × 10⁵ points) grows faster than the fastest k² in closed form", std::max(best - growth(T.k2_fastest), 0.0), 1e-15, "/s");
    s.require("d = 20: and the scan's best lies beside it (within 10⁻³ of the band)", std::fabs(arg - T.k2_fastest) < 1e-3 * (T.k2_high - T.k2_low));
    s.require("d = 20: the fastest k² lies below the onset wavenumber B/(2P)", T.k2_fastest < B / (2 * P) && T.k2_fastest > T.k2_low);
  }

  // ---- the modes of a chain ---------------------------------------------------------------------------------------------------------
  const std::size_t n = 54;
  const double spacing = 0.5, dv = 9;
  const Matrix L = pc::path_laplacian(n, spacing);
  {
    const std::vector<double> d{1, dv};
    const auto modes = pc::turing_modes(J, d, L);
    const auto T = pc::turing_two_species(J(0, 0), J(0, 1), J(1, 0), J(1, 1), 1, dv);
    double worst = 0;
    bool band = true;
    std::vector<std::size_t> unstable;
    for (std::size_t k = 0; k < n; k++) {
      const double Lk = modes.laplacian[k];
      worst = std::max(worst, std::fabs(modes.growth[k] - pc::turing_two_species_growth(J(0, 0), J(0, 1), J(1, 0), J(1, 1), 1, dv, Lk)));
      band = band && ((modes.growth[k] > 0) == (Lk > T.k2_low && Lk < T.k2_high));
      if (modes.growth[k] > 0) unstable.push_back(k);
    }
    s.check("turing_modes on the path Laplacian: each mode's growth is the dispersion relation at its Λ", worst, 1e-13);
    s.require("the unstable modes are exactly those with Λ in the band", band);
    s.require("here mode 5 alone, the fastest; Turing-unstable", unstable == std::vector<std::size_t>{5} && modes.fastest == 5 && modes.turing_unstable && modes.stable_homogeneous);
  }

  // ---- a pattern -------------------------------------------------------------------------------------------------------------------
  {
    pc::Network net;
    net.species = {"U", "V"};
    net.add("feed U", {}, {{"U", 1}}, Rational(Integer(1), Integer(10)));
    net.add("decay U", {{"U", 1}}, {}, Rational(1));
    net.add("autocatalysis", {{"U", 2}, {"V", 1}}, {{"U", 3}}, Rational(1));
    net.add("feed V", {}, {{"V", 1}}, Rational(Integer(9), Integer(10)));
    const pc::ReactionDiffusion rd(net, L, {1, dv});
    std::vector<double> c(2 * n);
    for (std::size_t v = 0; v < n; v++) c[2 * v] = a + b + 1e-3 * rng.normal(), c[2 * v + 1] = b / ((a + b) * (a + b)) + 1e-3 * rng.normal();
    for (int k = 0; k < 150000; k++) c = rd.step(c, 0.01);
    std::vector<double> u(n);
    double mean = 0;
    for (std::size_t v = 0; v < n; v++) u[v] = c[2 * v], mean += u[v] / n;
    std::size_t dominant = 0;
    double strongest = 0;
    for (std::size_t k = 1; k < n; k++) {
      double proj = 0;
      for (std::size_t v = 0; v < n; v++) proj += (u[v] - mean) * std::cos(std::numbers::pi * k * (v + 0.5) / n);
      if (std::fabs(proj) > strongest) strongest = std::fabs(proj), dominant = k;
    }
    int peaks = 0;
    for (std::size_t v = 0; v < n; v++) peaks += (v == 0 || u[v] > u[v - 1]) && (v + 1 == n || u[v] > u[v + 1]);
    double residual = 0;
    for (double r : rd.rate(c)) residual = std::max(residual, std::fabs(r));
    s.require("Schnakenberg on the chain: the steady pattern's dominant mode is the predicted 5, with three peaks",
              dominant == 5 && peaks == 3, "mode " + std::to_string(dominant) + ", " + std::to_string(peaks) + " peaks");
    s.check("Schnakenberg on the chain: the pattern is steady, |dc/dt| at t = 1500", residual, 1e-8);
    s.require("and of finite amplitude", *std::max_element(u.begin(), u.end()) - *std::min_element(u.begin(), u.end()) > 0.1);

    // rate() against mass_action() and a dense L D c.
    std::vector<double> x(2 * n);
    for (double& y : x) y = rng.uniform(0.2, 2.0);
    const auto r = rd.rate(x);
    double diff = 0;
    for (std::size_t v = 0; v < n; v++) {
      std::vector<double> cv{x[2 * v], x[2 * v + 1]}, dc(2);
      pc::mass_action(net, cv, dc);
      for (std::size_t sp = 0; sp < 2; sp++) {
        double lap = 0;
        for (std::size_t w = 0; w < n; w++) lap += L(v, w) * x[2 * w + sp];
        diff = std::max(diff, std::fabs(r[2 * v + sp] - (dc[sp] - (sp == 0 ? 1 : dv) * lap)));
      }
    }
    s.check("rate(): mass_action() at each vertex less L D c", diff, 1e-13);
  }

  // ---- conservation -------------------------------------------------------------------------------------------------------------------
  {
    pc::Network net;
    net.species = {"A", "B", "C"};
    net.add("bind", {{"A", 1}, {"B", 1}}, {{"C", 1}}, Rational(2));
    net.add("unbind", {{"C", 1}}, {{"A", 1}, {"B", 1}}, Rational(Integer(1), Integer(2)));
    const std::size_t V = 15;
    Matrix G(V);
    for (std::size_t i = 0; i < V; i++)
      for (std::size_t j = i + 1; j < V; j++)
        if (j == i + 1 || rng.uniform() < 0.2) {
          const double w = rng.uniform(0.2, 1.5);
          G(i, j) -= w, G(j, i) -= w, G(i, i) += w, G(j, j) += w;
        }
    const pc::ReactionDiffusion rd(net, G, {1.0, 0.5, 0.2});
    std::vector<double> c(3 * V);
    for (double& y : c) y = rng.uniform(0, 2);
    const auto totals = [&](const std::vector<double>& y) {
      double ac = 0, bc = 0;
      for (std::size_t v = 0; v < V; v++) ac += y[3 * v] + y[3 * v + 2], bc += y[3 * v + 1] + y[3 * v + 2];
      return std::vector<double>{ac, bc};
    };
    const auto t0 = totals(c);
    const double h = 0.5 * rd.diffusive_step();
    for (int k = 0; k < 2000; k++) c = rd.step(c, h);
    const auto t1 = totals(c);
    s.check("A + B ⇌ C on a random graph: the totals of A + C and B + C kept over 2000 RK4 steps (relative)",
            std::max(std::fabs(t1[0] / t0[0] - 1), std::fabs(t1[1] / t0[1] - 1)), 1e-13);
  }

  // ---- McKean's caricature: pinning ------------------------------------------------------------------------------------------------------
  {
    const std::size_t N = 80, front = N / 2 - 1;  // the last high site of the step is site 0 of the front
    double profile = 0, ends = 0;
    bool pinned = true, invades = true;
    std::string detail;
    for (double D : {0.25, 1.0, 3.0}) {
      const auto P = pc::mckean_pinning(D);
      const auto u = integrate_mckean(D, 0.5, N, 200, 0.05);
      ends = std::max({ends, std::fabs(u[front] - P.upper), std::fabs(u[front + 1] - P.lower)});
      for (std::size_t i = 20; i < 60; i++) profile = std::max(profile, std::fabs(u[i] - pc::mckean_front(static_cast<long long>(i) - static_cast<long long>(front), D)));
      for (double a_in : {P.lower + 0.01, P.upper - 0.01}) {
        std::vector<std::size_t> highs;
        integrate_mckean(D, a_in, N, 200, 0.05, &highs);
        pinned = pinned && std::all_of(highs.begin(), highs.end(), [&](std::size_t k) { return k == N / 2; });
      }
      std::vector<std::size_t> up, down;
      integrate_mckean(D, P.lower - 0.01, N, 200, 0.05, &up);
      integrate_mckean(D, P.upper + 0.01, N, 200, 0.05, &down);
      invades = invades && up.back() > N / 2 + 10 && down.back() + 10 < N / 2;
      detail += "D = " + std::to_string(D) + ": " + std::to_string(up.back()) + ", " + std::to_string(down.back()) + " high sites; ";
    }
    s.check("McKean, a = ½: the front at rest has u_0 = a_+ = ½(1 + 1/√(1 + 4D)) and u_1 = a_− (D = ¼, 1, 3)", ends, 1e-10);
    s.check("McKean, a = ½: the integrated front is the closed form μⁿ/(1 + μ), 1 − μ^{1−n}/(1 + μ)", profile, 1e-10);
    s.require("McKean: for a just inside (a_−, a_+) no site ever flips", pinned);
    s.require("McKean: just below a_− the state 1 invades, just above a_+ the state 0", invades, detail);
  }

  return s.done();
}
