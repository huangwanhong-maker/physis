// Checks of reaction–diffusion with kinetics of any form (field_kinetics.hpp; decision 0022).
//
// - Pure diffusion: a mode cos(2πi/nx) of a periodic grid decays at D (4/h²) sin²(π/nx) (relative 1e−6 against
//   e^{−λt}).
// - Mass: on a torus pure diffusion keeps the total to rounding.
// - Kinetics alone: a uniform state follows the local ODE x' = −x² exactly as 1/(1 + t).
// - A graph's Laplacian: the path of three vertices has rows (1, −1, 0), (−1, 2, −1), (0, −1, 1); a dish's grid has
//   no-flux walls, so diffusion keeps its total.
// - Fisher–KPP: a front of u' = D u'' + k u(1 − u) moves at 2√(D k) less Bramson's logarithmic correction
//   (within 1.5% over t ∈ [30, 60]), approaching 2√(D k) from below.
// - The split step: second order on a smooth problem; stable, and at each vertex's own equilibrium, for kinetics
//   that differ by vertex and are 36 times stiffer than the explicit method's bound.

#include <physis/chemistry/field_kinetics.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <numbers>

using namespace physis::chemistry;

int main() {
  calx::verification::Suite s;
  {
    const std::size_t nx = 32, ny = 8;
    const double h = 0.1, D = 0.3;
    FieldKinetics F(1, [](const double*, double* dc) { dc[0] = 0; }, SparseLaplacian::periodic_grid(nx, ny, h), {D});
    for (std::size_t j = 0; j < ny; j++)
      for (std::size_t i = 0; i < nx; i++) F.state()[j * nx + i] = std::cos(2 * std::numbers::pi * i / nx);
    double total0 = 0;
    for (double v : F.state()) total0 += v + 1;  // shifted, so that the total is not zero
    const double lambda = D * 4 / (h * h) * std::pow(std::sin(std::numbers::pi / nx), 2);
    F.run(1.0, 0.25 * F.diffusive_step());
    s.check("diffusion: a mode decays as e^{−λt} (relative)", F.at(0, 0) / std::exp(-lambda * F.time()) - 1, 1e-6);
    double total = 0;
    for (double v : F.state()) total += v + 1;
    s.check("diffusion: the total is kept on a torus (relative)", total / total0 - 1, 1e-14);
  }
  {
    FieldKinetics F(1, [](const double* c, double* dc) { dc[0] = -c[0] * c[0]; }, SparseLaplacian::periodic_grid(3, 3, 1.0), {1.0});
    for (double& v : F.state()) v = 1;
    F.run(2.0, 1e-3);
    s.check("kinetics: a uniform state follows x' = −x², 1/(1 + t)", F.at(4, 0) - 1.0 / 3, 1e-12);
  }
  {
    // A dish: no-flux walls. A corner has two neighbours and an interior site four; diffusion keeps the total.
    const auto L = SparseLaplacian::grid(5, 4, 1.0);
    s.require("dish: a corner has two neighbours, an interior site four", L.rows[0].size() == 3 && L.rows[6].size() == 5);
    FieldKinetics F(1, [](const double*, double* dc) { dc[0] = 0; }, L, {0.7});
    double total0 = 0;
    for (std::size_t v = 0; v < F.vertices(); v++) total0 += (F.state()[v] = std::sin(1.3 * v) + 2);
    F.run(5.0, 0.2 * F.diffusive_step());
    double total = 0;
    for (double c : F.state()) total += c;
    s.check("dish: diffusion keeps the total within no-flux walls (relative)", total / total0 - 1, 1e-14);
  }
  {
    const auto L = SparseLaplacian::graph(3, {{0, 1, 1.0}, {1, 2, 1.0}});
    auto entry = [&](std::size_t r, std::size_t c) {
      double v = 0;
      for (const auto& [k, w] : L.rows[r])
        if (k == c) v += w;
      return v;
    };
    s.require("graph: the path's Laplacian", entry(0, 0) == 1 && entry(0, 1) == -1 && entry(1, 1) == 2 && entry(1, 2) == -1 && entry(2, 2) == 1 && entry(0, 2) == 0);
  }
  {
    // A pulled front converges slowly: x(t) = 2√(Dk) t − (3/2)√(D/k) ln t + C (Bramson 1983, Mem. Amer. Math. Soc.
    // 44, 285; Ebert and van Saarloos 2000, Physica D 146, 1), so its mean speed over [t₁, t₂] is
    // 2√(Dk) − (3/2)√(D/k) ln(t₂/t₁)/(t₂ − t₁).
    const std::size_t n = 2000;
    const double h = 0.1, D = 1.0, k = 1.0;
    std::vector<std::tuple<std::size_t, std::size_t, double>> edges;
    for (std::size_t i = 0; i + 1 < n; i++) edges.push_back({i, i + 1, 1 / (h * h)});
    FieldKinetics F(1, [&](const double* c, double* dc) { dc[0] = k * c[0] * (1 - c[0]); }, SparseLaplacian::graph(n, edges), {D});
    for (std::size_t i = 0; i < 20; i++) F.state()[i] = 1;
    auto front = [&]() {
      for (std::size_t i = 0; i < n; i++)
        if (F.at(i, 0) < 0.5) return i * h;
      return n * h;
    };
    const double dt = 0.2 * F.diffusive_step();
    F.run(30.0, dt);
    const double x1 = front(), t1 = F.time();
    F.run(30.0, dt);
    const double t2 = F.time(), speed = (front() - x1) / (t2 - t1);
    const double bramson = 2 * std::sqrt(D * k) - 1.5 * std::sqrt(D / k) * std::log(t2 / t1) / (t2 - t1);
    s.check("Fisher–KPP: the front's speed, with Bramson's correction (relative)", speed / bramson - 1, 0.015);
    s.require("Fisher–KPP: and below 2√(D k), approaching it from below", speed < 2 * std::sqrt(D * k));
  }
  {
    // The split step: second order on a smooth problem (error ratio ≈ 4 when the step halves, against a fine run).
    auto make = [] {
      FieldKinetics F(2, [](const double* c, double* dc) { dc[0] = c[1] - c[0] * c[0]; dc[1] = -c[0]; }, SparseLaplacian::periodic_grid(16, 3, 0.5), {0.2, 0.1});
      for (std::size_t v = 0; v < F.vertices(); v++) F.state()[2 * v] = 0.5 + 0.3 * std::sin(2 * std::numbers::pi * (v % 16) / 16.0), F.state()[2 * v + 1] = 0.2;
      return F;
    };
    auto fine = make();
    fine.run(1.0, 1e-4);
    auto error = [&](double h) {
      auto F = make();
      F.run_split(1.0, h);
      double e = 0;
      for (std::size_t i = 0; i < F.state().size(); i++) e = std::max(e, std::fabs(F.state()[i] - fine.state()[i]));
      return e;
    };
    s.check("split step: second order (error ratio / 4 − 1)", error(0.05) / error(0.025) / 4 - 1, 0.2);
    // Stiff kinetics, vertex by vertex: c' = −10⁴ (c − v/10), at a step 36 times the explicit bound of RK4.
    FieldKinetics S(1, FieldKinetics::VertexKinetics([](std::size_t v, const double* c, double* dc) { dc[0] = -1e4 * (c[0] - 0.1 * v); }),
                    SparseLaplacian::periodic_grid(5, 3, 1.0), {0.01});
    S.run_split(0.1, 0.01);
    double worst = 0;
    bool finite = true;
    for (std::size_t v = 0; v < S.vertices(); v++) finite = finite && std::isfinite(S.at(v, 0));
    for (std::size_t v = 0; v < S.vertices(); v++) worst = std::max(worst, std::fabs(S.at(v, 0) - 0.1 * v));
    s.require("split step: stable for kinetics 36 times stiffer than RK4's bound", finite);
    s.check("split step: reaches each vertex's own equilibrium (diffusion's pull 1e−5 at most)", worst, 1e-4);
  }
  return s.done();
}
