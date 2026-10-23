// Checks of physis/quantum_gravity/condensate.hpp: the conserved charge and energy along
// integrated modes; the closed form against Dormand–Prince, and for Q = 0 against
// ρ = A e^{mφ} + B e^{−mφ}; the effective Friedmann equation of one spin as an identity along
// the solution, and V''/V = 4m² + 2ΣV_j E_j / V by finite differences; the classical limit
// (V'/V)² → 12πG at large volume on both branches; the bounce at ρ_φ = ρ_c; a universe of two
// spins as the sum of its modes; no bounce without charge; the Landau condensate; and the
// classical limit's source, physis spacetime's flat FLRW with a ∝ t^{1/3}, whose Einstein
// tensor is that of a stiff fluid (a massless scalar) with (V'/V)² = 12πG in relational time.

#include <physis/quantum_gravity/condensate.hpp>
#include <physis/geometry/curvature.hpp>
#include <physis/spacetime/models.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>

namespace pq = physis::quantum_gravity;
namespace pg = physis::geometry;
namespace ps = physis::spacetime;
using std::numbers::pi;

int main() {
  calx::verification::Suite s;
  const double G = 1.0;
  const double m2 = 3 * pi * G, m = std::sqrt(m2);

  // ---- conservation and the closed form ------------------------------------------------------------------------------
  {
    const pq::CondensateMode md{1.5, 2.0, 0.7, 0.3};
    std::vector<double> grid;
    for (int i = 0; i <= 400; i++) grid.push_back(-1.2 + 2.4 * i / 400.0);
    const auto x = pq::integrate_mode(md, m2, grid);
    double dE = 0, drho = 0, dQ = 0;
    for (std::size_t i = 0; i < grid.size(); i++) {
      const double rho = x[i][0], dr = x[i][1];
      const double E = dr * dr + md.Q * md.Q / (rho * rho) - m2 * rho * rho;
      dE = std::max(dE, std::fabs(E - md.E) / (md.Q * md.Q / (rho * rho) + m2 * rho * rho));
      const auto ex = pq::mode_exact(md, m2, grid[i]);
      drho = std::max(drho, std::fabs(rho / ex.rho - 1));
      // θ' = Q/ρ²: the charge read back by a centred difference of the integrated phase.
      if (i > 0 && i + 1 < grid.size()) {
        const double dtheta = (x[i + 1][2] - x[i - 1][2]) / (grid[i + 1] - grid[i - 1]);
        dQ = std::max(dQ, std::fabs(dtheta * rho * rho / md.Q - 1));
      }
    }
    s.check("condensate: E = ρ'² + Q²/ρ² − m²ρ² conserved along Dormand–Prince (relative to its terms)", dE, 1e-10);
    s.check("condensate: Q = ρ² θ' conserved (centred differences, step 0.006)", dQ, 1e-3);
    s.check("condensate: Dormand–Prince against the closed form u = −E/2m² + C⁺e^{2mφ} + C⁻e^{−2mφ} (relative)", drho, 1e-10);
  }
  {
    // Q = 0: ρ = A e^{mφ} + B e^{−mφ}, with E = −4m²AB and the bounce at ln(B/A)/2m.
    const double A = 0.3, B = 2.0;
    const pq::CondensateMode md{0.5, 0.0, -4 * m2 * A * B, std::log(B / A) / (2 * m)};
    double dev = 0;
    std::vector<double> grid;
    for (double phi = -1; phi <= 1.0001; phi += 0.05) {
      grid.push_back(phi);
      const double want = A * std::exp(m * phi) + B * std::exp(-m * phi);
      dev = std::max(dev, std::fabs(pq::mode_exact(md, m2, phi).rho / want - 1));
    }
    const auto x = pq::integrate_mode(md, m2, grid);
    for (std::size_t i = 0; i < grid.size(); i++)
      dev = std::max(dev, std::fabs(x[i][0] / (A * std::exp(m * grid[i]) + B * std::exp(-m * grid[i])) - 1));
    s.check("condensate: without charge, ρ = A e^{mφ} + B e^{−mφ}, closed form and integrated (relative)", dev, 1e-10);
    const pq::CondensateMode bare{0.5, 0.0, 1.0, 0.0};
    s.require("condensate: without charge and with E > 0 the mode empties (u_min <= 0): no bounce", !(pq::mode_u_min(bare, m2) > 0));
  }
  // ---- one spin: the effective Friedmann equation and the bounce --------------------------------------------------------------
  {
    pq::Universe U;
    U.G = G;
    U.modes = {{2.0, 3.0, 0.0, 0.0}};
    const auto C = pq::evolve(U, -1.5, 1.5, 3001, true);
    double dev = 0;
    for (std::size_t i = 0; i < C.phi.size(); i++) dev = std::max(dev, std::fabs(C.hubble_squared[i] - C.friedmann[i]) / (4 * pi * G / 3));
    s.check("condensate: one spin, E = 0: (V'/3V)² = (4πG/3)(1 − ρ_φ/ρ_c) along the integrated solution (relative to 4πG/3)", dev, 1e-9);
    const double Vj = U.quantum_volume(2.0), rho_c = m2 / (2 * Vj * Vj);
    s.check("condensate: ρ_c = m²/(2V_j²) = 3πG/(2V_j²) (relative)", C.critical_density / rho_c - 1, 1e-15);
    s.check("condensate: the bounce is at ρ_φ = ρ_c (relative, parabola through the grid)", C.bounce_density / rho_c - 1, 1e-6);
    s.check("condensate: the bounce is at the mode's φ_b = 0", C.bounce_phi, 1e-6);
    s.require("condensate: a charged condensate bounces", C.bounces);
    // The classical limit on both branches: (V'/V)² → 12πG, V' > 0 after the bounce and < 0 before.
    const double rate = 9 * C.hubble_squared.back() / C.classical_rate - 1, rate0 = 9 * C.hubble_squared.front() / C.classical_rate - 1;
    s.check("condensate: (V'/V)² → 12πG at large volume, expanding branch (relative, φ = 1.5)", rate, 1e-6);
    s.check("condensate: and on the contracting branch (relative, φ = −1.5)", rate0, 1e-6);
    s.require("condensate: expanding after the bounce, contracting before", C.hubble.back() > 0 && C.hubble.front() < 0);
    s.check("condensate: V approaches the classical expanding solution V ∝ e^{√(12πG) φ} (relative, φ = 1.5)", C.volume.back() / C.classical_expanding.back() - 1, 1e-6);
  }
  {
    // E ≠ 0: the 4 V_j E / 9V term, and V''/V = 4m² + 2 V_j E / V by second differences.
    pq::Universe U;
    U.modes = {{1.0, 1.5, -2.0, 0.2}};
    const auto C = pq::evolve(U, -1.0, 1.0, 2001, true);
    double dev = 0, acc = 0;
    for (std::size_t i = 1; i + 1 < C.phi.size(); i++) {
      dev = std::max(dev, std::fabs(C.hubble_squared[i] - C.friedmann[i]) / (4 * pi * G / 3));
      const double h = C.phi[1] - C.phi[0];
      const double d2 = (C.volume[i + 1] - 2 * C.volume[i] + C.volume[i - 1]) / (h * h) / C.volume[i];
      acc = std::max(acc, std::fabs(d2 / C.acceleration[i] - 1));
    }
    s.check("condensate: one spin, E ≠ 0: (V'/3V)² = (4πG/3)(1 − ρ_φ/ρ_c) + 4V_jE/9V (relative to 4πG/3)", dev, 1e-9);
    s.check("condensate: V''/V = 4m² + 2V_jE/V (second differences, step 0.001)", acc, 1e-5);
  }
  // ---- two spins -------------------------------------------------------------------------------------------------------------
  {
    pq::Universe U, U1, U2;
    U.modes = {{0.5, 1.0, 0.4, -0.2}, {1.5, 2.5, -0.3, 0.4}};
    U1.modes = {U.modes[0]};
    U2.modes = {U.modes[1]};
    const auto C = pq::evolve(U, -1, 1, 401, true), C1 = pq::evolve(U1, -1, 1, 401), C2 = pq::evolve(U2, -1, 1, 401);
    double dev = 0, acc = 0;
    for (std::size_t i = 0; i < C.phi.size(); i++) {
      dev = std::max(dev, std::fabs(C.volume[i] / (C1.volume[i] + C2.volume[i]) - 1));
      if (i > 0 && i + 1 < C.phi.size()) {
        const double h = C.phi[1] - C.phi[0];
        const double d2 = (C.volume[i + 1] - 2 * C.volume[i] + C.volume[i - 1]) / (h * h) / C.volume[i];
        acc = std::max(acc, std::fabs(d2 / C.acceleration[i] - 1));
      }
    }
    s.check("condensate: two spins: V is the sum of the modes' volumes (relative)", dev, 1e-10);
    s.check("condensate: two spins: π_φ = Q₁ + Q₂", C.pi_phi - 3.5, 1e-15);
    s.check("condensate: two spins: V''/V = 4m² + 2ΣV_jE_j/V (second differences, step 0.005)", acc, 1e-4);
    // The corrections fall as e^{−2mφ} (the E terms) and e^{−4mφ}: far out they vanish.
    const auto far = pq::evolve(U, 2.4, 2.5, 3);
    s.check("condensate: two spins: (V'/V)² → 12πG at large volume (relative, φ = 2.5)", 9 * far.hubble_squared.back() / far.classical_rate - 1, 1e-6);
  }
  // ---- the Landau condensate ---------------------------------------------------------------------------------------------------
  {
    const auto L = pq::landau_condensate(-0.8, 0.5, 4), N = pq::landau_condensate(0.3, 0.5, 4);
    const double slope = 2 * (-0.8) * L.amplitude + 4 * 0.5 * std::pow(L.amplitude, 3);
    s.require("condensate: Landau: condensed for μ < 0, at a stationary point of lower potential; symmetric for μ > 0",
              L.condensed && std::fabs(slope) < 1e-14 && L.potential < 0 && !N.condensed);
  }
  // ---- the classical limit's source: flat FLRW with a ∝ t^{1/3} --------------------------------------------------------------
  {
    // A massless scalar is a stiff fluid, p = ρ: G^i_i = G_00. In relational time,
    // (V'/V)² = (3H)² / φ̇² = 9H² / 2ρ = 12πG, with 8πG ρ = G_00.
    const ps::FLRW<ps::PowerLaw> u{0.0, ps::PowerLaw{1.0, 1.0 / 3}};
    double stiff = 0, rel = 0;
    for (double t : {0.5, 1.0, 3.0}) {
      const auto k = pg::curvature(u, pg::Point<4>{t, 1.0, 1.0, 0.3});
      const auto Gm = pg::einstein_mixed(k);
      const double G00 = k.Einstein[0][0], H = 1.0 / (3 * t);
      stiff = std::max(stiff, std::fabs(Gm[1][1] / G00 - 1));
      const double rho = G00 / (8 * pi * G);
      rel = std::max(rel, std::fabs(9 * H * H / (2 * rho) / (12 * pi * G) - 1));
    }
    s.check("condensate: physis's FLRW with a ∝ t^{1/3} is a stiff fluid, G^i_i = G_00 (relative)", stiff, 1e-12);
    s.check("condensate: and in relational time (V'/V)² = 12πG, the condensate's classical limit (relative)", rel, 1e-12);
  }
  return s.done();
}
