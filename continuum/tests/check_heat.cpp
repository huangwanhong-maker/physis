// Checks of heat conduction (heat.hpp) against closed forms and its own invariants.
//
// - Two cells: the difference of their temperatures falls by R(hλ) a step, λ = G (1/C₁ + 1/C₂), their
//   energy constant.
// - Steady states that the two-point flux reproduces exactly: a linear profile between two held faces,
//   and a profile of two slopes across two materials, (k₁T₁ + k₂T₂)/(k₁ + k₂) at their interface, with
//   one flux through every section.
// - The budget of every step with held, exchanging, heated and insulated faces and sources, two
//   materials and a loose solve (10⁻⁶): storage + outflow = sources to the rounding of the stored energy,
//   and so for each material with the heat through their interface (formalism Propositions 5.3, 5.6).
// - Dissipation: ½ Σ C T² does not rise under Crank–Nicolson or implicit Euler at a hundred times the
//   explicit bound, nor under explicit Euler below 2/λ_max; above it, it grows (λ_max by power iteration;
//   the Gershgorin bound below it).
// - The decaying sine mode of the heat equation in the unit cube held at 0 on its faces,
//   sin πx sin πy sin πz e^{−3π²κt}: second-order convergence in the cells (n = 6, 12, 24), and in time
//   second order for Crank–Nicolson, first for implicit Euler.
// - The temperature as a field: a cochain read as its cells' densities, its versions committed.

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <numbers>
#include <string>
#include <vector>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <physis/continuum/heat.hpp>
#include <physis/discretisation/cochain_field.hpp>
#include <physis/fields/history.hpp>

using calx::verification::Suite;
namespace pc = physis::continuum;
namespace dz = physis::discretisation;
using pc::Boundary;
using pc::Vec3;

namespace {

const pc::Material unit{"unit", 1, 1, 1};

pc::HeatConduction::BoundaryOf held_on_x(double T1, double T2) {
  return [=](std::size_t, const Vec3&, const Vec3& out) {
    if (out[0] < -0.5) return Boundary::temperature(T1);
    if (out[0] > 0.5) return Boundary::temperature(T2);
    return Boundary::insulated();
  };
}

std::vector<double> random_temperatures(calx::random::Stream& rng, std::size_t n, double lo, double hi) {
  std::vector<double> T(n);
  for (double& t : T) t = rng.uniform(lo, hi);
  return T;
}

double total(const std::vector<double>& u) {
  double s = 0;
  for (double x : u) s += x;
  return s;
}

}  // namespace

int main() {
  Suite s;
  calx::random::Stream rng(20261029, 4, 0, 0);

  // ---- two cells ------------------------------------------------------------------------------------------------
  {
    const auto box = dz::uniform_box({0, 0, 0}, {2, 1, 1}, {2, 1, 1});
    const pc::Material m{"m", 1000, 1000, 1};
    const pc::HeatConduction heat(box.complex, {0, 0}, {m}, [](std::size_t, const Vec3&, const Vec3&) { return Boundary::insulated(); });
    const double C = 1e6, G = 1, lambda = G * (1 / C + 1 / C), h = 1e5;
    bool exact = true;
    double worst = 0, mean = 0;
    for (double theta : {0.0, 0.5, 1.0}) {
      auto st = heat.state(0, {310, 290});
      for (int k = 0; k < 10; k++) {
        const auto next = heat.step(st, h, theta, 1e-15).state;
        const double ratio = (next.temperature[0] - next.temperature[1]) / (st.temperature[0] - st.temperature[1]);
        worst = std::max(worst, std::fabs(ratio - calx::integrators::theta_factor(h * lambda, theta)));
        mean = std::max(mean, std::fabs(0.5 * (next.temperature[0] + next.temperature[1]) - 300));
        st = next;
      }
      exact = exact && heat.capacity()[0] == C && heat.conductance()[0] + heat.conductance()[1] >= 0;
    }
    s.check("two cells: their difference falls by R(hλ) a step, λ = G(1/C₁ + 1/C₂), θ = 0, ½, 1", worst, 1e-13);
    s.check("two cells: their mean stays", mean, 1e-12, "K");
  }

  // ---- steady states reproduced exactly -------------------------------------------------------------------------------
  {
    const auto box = dz::uniform_box({0, 0, 0}, {4, 1.5, 1}, {8, 3, 2});
    const pc::HeatConduction heat(box.complex, std::vector<int>(box.cells[3].size(), 0), {{"rock", 2700, 790, 2.9}}, held_on_x(300, 400));
    const auto T = heat.steady();
    const auto& c = heat.measures().centroid[3];
    double worst = 0;
    for (std::size_t i = 0; i < T.size(); i++) worst = std::max(worst, std::fabs(T[i] - (300 + 25 * c[i][0])));
    s.check("steady: a linear profile between faces held at 300 K and 400 K, at every cell", worst, 1e-10, "K");
    const auto J = heat.flux(T);
    const auto& F = heat.facets();
    double spread = 0;
    for (int plane = 1; plane < 8; plane++) {
      double q = 0;
      for (std::size_t f = 0; f < J.size(); f++)
        if (std::fabs(F.conormal[f][0]) > 0.5 && std::fabs(F.centroid[f][0] - 0.5 * plane) < 1e-12) q += J[f] * F.conormal[f][0];
      spread = std::max(spread, std::fabs(q - (-2.9 * 1.5 * 100 / 4)));
    }
    s.check("steady: the same heat through every section, k A ΔT/L = 108.75 W toward the cold face", spread, 1e-9, "W");

    const auto two = dz::uniform_box({0, 0, 0}, {4, 1, 1}, {8, 2, 2});
    std::vector<int> occupancy(two.cells[3].size());
    const auto M2 = dz::measure(two.complex);
    for (std::size_t i = 0; i < occupancy.size(); i++) occupancy[i] = M2.centroid[3][i][0] < 2 ? 0 : 1;
    const pc::HeatConduction layered(two.complex, occupancy, {{"k4", 1000, 1000, 4}, {"k1", 1000, 1000, 1}}, held_on_x(300, 400));
    const auto Tl = layered.steady();
    double profile = 0;
    for (std::size_t i = 0; i < Tl.size(); i++) {
      const double x = M2.centroid[3][i][0];
      profile = std::max(profile, std::fabs(Tl[i] - (x < 2 ? 300 + 10 * x : 320 + 40 * (x - 2))));
    }
    s.check("steady: across two materials, two slopes meeting at (k₁T₁ + k₂T₂)/(k₁ + k₂) = 320 K", profile, 1e-10, "K");
  }

  // ---- the budget of every step --------------------------------------------------------------------------------------
  {
    const auto box = dz::uniform_box({0, 0, 0}, {3, 2.5, 2}, {6, 5, 4});
    const auto M = dz::measure(box.complex);
    const std::size_t N = box.cells[3].size();
    std::vector<int> occupancy(N);
    for (std::size_t i = 0; i < N; i++) occupancy[i] = M.centroid[3][i][0] < 1.5 ? 0 : 1;
    std::vector<double> sources(N);
    for (double& q : sources) q = rng.uniform(0, 10);
    const pc::HeatConduction heat(box.complex, occupancy, {{"granite", 2700, 790, 2.9}, {"sandstone", 2200, 920, 1.8}},
                                  [](std::size_t, const Vec3&, const Vec3& out) {
                                    if (out[0] < -0.5) return Boundary::temperature(280);
                                    if (out[0] > 0.5) return Boundary::exchange(10, 300);
                                    if (out[2] < -0.5) return Boundary::flux(50);
                                    if (out[2] > 0.5) return Boundary::flux(-20);
                                    return Boundary::insulated();
                                  },
                                  sources);
    auto st = heat.state(0, random_temperatures(rng, N, 280, 320));
    double worst = 0, region = 0;
    const auto& F = heat.facets();
    for (int k = 0; k < 50; k++) {
      const auto step = heat.step(st, 600, 0.5, 1e-6);
      worst = std::max(worst, std::fabs(step.budget.residual()) / total(step.state.energy));
      // Each material's own budget, with the heat through their interface.
      std::array<double, 2> storage{}, outflow{}, sources_r{};
      for (std::size_t i = 0; i < N; i++) {
        storage[static_cast<std::size_t>(occupancy[i])] += step.state.energy[i] - st.energy[i];
        sources_r[static_cast<std::size_t>(occupancy[i])] += 600 * sources[i];
      }
      const auto& J = step.flux;
      for (std::size_t f = 0; f < J.size(); f++) {
        if (F.boundary(f)) {
          const auto [a, sign] = F.bounded(f);
          outflow[static_cast<std::size_t>(occupancy[a])] += 600 * sign * J[f];
        } else {
          const std::size_t a = F.cells[f][0], b = F.cells[f][1];
          if (occupancy[a] != occupancy[b]) {
            outflow[static_cast<std::size_t>(occupancy[a])] += 600 * J[f];
            outflow[static_cast<std::size_t>(occupancy[b])] -= 600 * J[f];
          }
        }
      }
      for (int r = 0; r < 2; r++) region = std::max(region, std::fabs(storage[static_cast<std::size_t>(r)] + outflow[static_cast<std::size_t>(r)] - sources_r[static_cast<std::size_t>(r)]) / total(step.state.energy));
      st = step.state;
    }
    s.check("budget: storage + outflow = sources at every step, held, exchanging, heated and insulated faces, a solve to 10⁻⁶ (of the stored energy)", worst, 1e-15);
    s.check("budget: and for each material, with the heat through their interface: what leaves one enters the other (of the stored energy)", region, 1e-15);

    const pc::HeatConduction closed(box.complex, occupancy, {{"granite", 2700, 790, 2.9}, {"sandstone", 2200, 920, 1.8}},
                                    [](std::size_t, const Vec3&, const Vec3&) { return Boundary::insulated(); });
    auto c0 = closed.state(0, random_temperatures(rng, N, 250, 350));
    const double U0 = total(c0.energy);
    auto c = c0;
    for (int k = 0; k < 100; k++) c = closed.step(c, 3600, 0.5, 1e-6).state;
    s.check("budget: an insulated body without sources keeps its energy over 100 steps, even with a solve to 10⁻⁶ (of the energy)", (total(c.energy) - U0) / U0, 1e-14);
  }

  // ---- dissipation and stability --------------------------------------------------------------------------------------
  {
    const auto box = dz::uniform_box({0, 0, 0}, {1, 0.8, 0.6}, {5, 4, 3});
    const std::size_t N = box.cells[3].size();
    const pc::HeatConduction heat(box.complex, std::vector<int>(N, 0), {{"copper", 8960, 385, 401}}, [](std::size_t, const Vec3&, const Vec3&) { return Boundary::insulated(); });
    // λ_max of C⁻¹L by power iteration.
    std::vector<double> x = random_temperatures(rng, N, -1, 1), y(N);
    double lambda = 0;
    for (int k = 0; k < 2000; k++) {
      calx::linalg::multiply(heat.matrix(), std::span<const double>(heat.matrix_values()), std::span<const double>(x), std::span<double>(y));
      double xLx = 0, xCx = 0, norm = 0;
      for (std::size_t i = 0; i < N; i++) xLx += x[i] * y[i], xCx += x[i] * heat.capacity()[i] * x[i];
      lambda = xLx / xCx;
      for (std::size_t i = 0; i < N; i++) y[i] /= heat.capacity()[i], norm = std::max(norm, std::fabs(y[i]));
      for (std::size_t i = 0; i < N; i++) x[i] = y[i] / norm;
    }
    const auto T0 = random_temperatures(rng, N, 290, 310);
    auto rises = [&](double h, double theta, int steps) {
      auto st = heat.state(0, T0);
      bool rose = false;
      double last = heat.quadratic(st.temperature);
      for (int k = 0; k < steps; k++) {
        st = heat.step(st, h, theta, 1e-14).state;
        const double e = heat.quadratic(st.temperature);
        if (e > last * (1 + 1e-14)) rose = true;
        last = e;
      }
      return std::pair{rose, last};
    };
    const double big = 100 * heat.explicit_step_bound();
    s.require("dissipation: ½ΣCT² never rises under Crank–Nicolson and implicit Euler at 100 times the explicit bound",
              !rises(big, 0.5, 200).first && !rises(big, 1.0, 200).first);
    s.require("dissipation: Gershgorin's explicit bound lies below 2/λ_max", heat.explicit_step_bound() <= 2 / lambda);
    const auto below = rises(0.95 * 2 / lambda, 0.0, 400), above = rises(1.05 * 2 / lambda, 0.0, 400);
    s.require("stability: explicit Euler at 0.95·2/λ_max never raises it; at 1.05·2/λ_max it grows beyond its start",
              !below.first && above.first && above.second > heat.quadratic(T0));
  }

  // ---- the decaying sine mode -----------------------------------------------------------------------------------------
  {
    const double lambda = 3 * std::numbers::pi * std::numbers::pi, t_end = 1 / lambda;
    auto run = [&](int n, int steps, double theta) {
      const auto box = dz::uniform_box({0, 0, 0}, {1, 1, 1}, {n, n, n});
      const pc::HeatConduction heat(box.complex, std::vector<int>(box.cells[3].size(), 0), {unit}, [](std::size_t, const Vec3&, const Vec3&) { return Boundary::temperature(0); });
      const auto& c = heat.measures().centroid[3];
      std::vector<double> T(c.size());
      for (std::size_t i = 0; i < c.size(); i++) T[i] = std::sin(std::numbers::pi * c[i][0]) * std::sin(std::numbers::pi * c[i][1]) * std::sin(std::numbers::pi * c[i][2]);
      auto st = heat.state(0, T);
      for (int k = 0; k < steps; k++) st = heat.step(st, t_end / steps, theta, 1e-13).state;
      return std::pair{st.temperature, c};
    };
    auto error = [&](int n) {
      const auto [T, c] = run(n, 64, 0.5);
      double e = 0;
      for (std::size_t i = 0; i < T.size(); i++)
        e = std::max(e, std::fabs(T[i] - std::exp(-lambda * t_end) * std::sin(std::numbers::pi * c[i][0]) * std::sin(std::numbers::pi * c[i][1]) * std::sin(std::numbers::pi * c[i][2])));
      return e;
    };
    const double e6 = error(6), e12 = error(12), e24 = error(24);
    s.check("sine mode: second order in the cells, log₂(e₆/e₁₂) − 2", std::log2(e6 / e12) - 2, 0.15);
    s.check("sine mode: and log₂(e₁₂/e₂₄) − 2", std::log2(e12 / e24) - 2, 0.1);
    s.check("sine mode: the error at n = 24 against the exact decay, of the amplitude e⁻¹ (a guard on its size, 1.4·10⁻³ when written)", e24 / std::exp(-1.0), 2e-3);
    // In time, against the same cells stepped finely.
    const auto reference = run(6, 1024, 0.5).first;
    auto in_time = [&](int steps, double theta) {
      const auto T = run(6, steps, theta).first;
      double e = 0;
      for (std::size_t i = 0; i < T.size(); i++) e = std::max(e, std::fabs(T[i] - reference[i]));
      return e;
    };
    s.check("time: Crank–Nicolson second order, log₂(e₈/e₁₆) − 2", std::log2(in_time(8, 0.5) / in_time(16, 0.5)) - 2, 0.05);
    s.check("time: implicit Euler first order, log₂(e₁₆/e₃₂) − 1", std::log2(in_time(16, 1.0) / in_time(32, 1.0)) - 1, 0.05);
  }

  // ---- the temperature as a field, and its versions -----------------------------------------------------------------
  {
    const auto box = std::make_shared<const dz::BoxComplex>(dz::uniform_box({0, 0, 0}, {2, 1, 1}, {4, 2, 2}));
    const pc::HeatConduction heat(box->complex, std::vector<int>(box->cells[3].size(), 0), {unit}, held_on_x(10, 20));
    const auto vol = heat.measures().of(3);
    auto as_field = [&](const pc::HeatConduction::State& st) {
      std::vector<double> c(st.temperature.size());
      for (std::size_t i = 0; i < c.size(); i++) c[i] = st.temperature[i] * vol[i];
      physis::fields::Descriptor d;
      d.name = "temperature";
      d.fibre.dimension = physis::units::dim::temperature;
      return dz::CochainField(d, box, 3, c);
    };
    physis::fields::Descriptor dh;
    dh.name = "temperature";
    physis::fields::History<dz::CochainField> history(dh);
    auto st = heat.state(0, std::vector<double>(box->cells[3].size(), 15));
    history.commit(0, as_field(st));
    st = heat.step(st, 0.1).state;
    history.commit(0.1, as_field(st));
    const Vec3 x{0.1, 0.3, 0.6};
    const auto L = box->locate(x);
    const auto now = history.at({0.1, x}), then = history.at({0.05, x});
    s.require("field: the temperature read in its cell at the committed instant; between commits the earlier version, with its note",
              std::fabs(now.value - st.temperature[L->cell]) < 1e-12 && std::fabs(then.value - 15) < 1e-12 && !then.note.empty() && history.versions() == 2);
  }
  return s.done();
}
