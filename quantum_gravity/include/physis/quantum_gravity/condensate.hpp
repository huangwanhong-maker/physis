// The cosmology of group-field-theory condensates: a continuum spacetime from the condensed
// phase of the quanta of geometry.
//
// In group field theory the quanta of geometry are excitations of a field φ(g_1, ..., g_4) on
// SU(2)^4; a homogeneous universe is a condensate, a coherent state of many quanta, whose
// mean field σ obeys the classical equations of the GFT action, a Gross–Pitaevskii equation
// (Gielen, Oriti and Sindoni 2013, Phys. Rev. Lett. 111, 031301; 2014, JHEP 06, 013). Time is
// relational: a massless scalar field φ, added to the field's arguments, is the clock. For an
// isotropic condensate the field decomposes into one mode per spin j, σ_j(φ), and with the
// interactions negligible (the free regime, away from the largest densities) each mode obeys
// (Oriti, Sindoni and Wilson-Ewing 2016, Class. Quantum Grav. 33, 224001, §3)
//
//   A_j σ_j'' − B_j σ_j = 0,   so   σ_j'' = m_j² σ_j,   m_j² = B_j / A_j,
//
// a prime being d/dφ. Writing σ_j = ρ_j e^{iθ_j}, the modulus and phase obey
//
//   ρ_j'' − Q_j² / ρ_j³ − m_j² ρ_j = 0,   θ_j' = Q_j / ρ_j²,
//
// with two conserved quantities per mode, the "charge" Q_j and the "energy"
//
//   E_j = (ρ_j')² + Q_j² / ρ_j² − m_j² ρ_j².
//
// The modes decouple, and each has a closed form: u_j = ρ_j² obeys u_j'' = 2 E_j + 4 m_j² u_j
// (differentiate twice and use the two equations above), so
//
//   u_j(φ) = −E_j / (2 m_j²) + C⁺_j e^{2 m_j φ} + C⁻_j e^{−2 m_j φ},
//   C⁺_j C⁻_j = (E_j² + 4 m_j² Q_j²) / (16 m_j⁴),
//
// the product fixed by (u_j')² = 4 (E_j u_j − Q_j² + m_j² u_j²). With Q_j ≠ 0 (or E_j < 0) the
// minimum of u_j, (−E_j + √(E_j² + 4 m_j² Q_j²)) / (2 m_j²), is positive: the mode never
// empties, and the universe bounces. The bounce of mode j is placed at φ_j by
// C^±_j = c_j e^{∓2 m_j φ_j}, c_j = √(E_j² + 4 m_j² Q_j²) / (4 m_j²).
//
// Observables (OSW-E 2016, §4):
// - the volume V = Σ_j V_j ρ_j², with V_j the volume of one quantum of spin j. Here
//   V_j = v₀ j^{3/2} ℓ_P³, the scaling of the eigenvalues of the volume operator of loop
//   quantum gravity on a 4-valent node with equal spins (Rovelli and Smolin 1995; Brunnemann
//   and Thiemann 2006); v₀ (default 1) sets the unit;
// - the momentum of the clock π_φ = ħ Σ_j Q_j, and its energy density ρ_φ = π_φ² / (2 V²);
// - the effective Friedmann equations: with m_j² = m² for every j,
//     V'' / V = 4 m² + 2 Σ_j V_j E_j / V,
//     (V' / 3V)² = (2 Σ_j V_j ρ_j ρ_j' / 3V)²,
//   and for one spin, V = V_j ρ², the exact form
//     (V' / 3V)² = (4πG/3)(1 − ρ_φ / ρ_c) + 4 V_j E / (9 V),   ρ_c = m² / (2 V_j²) = 3πG / (2 V_j²),
//   when m² = 3πG (derived here from the equations above: (V'/3V)² = (4/9)(E/ρ² − Q²/ρ⁴ + m²)
//   and Q²/ρ⁴ = 2 V_j² ρ_φ). At large volume the corrections vanish and (V'/V)² → 4 m² = 12πG,
//   the classical Friedmann equation of a flat universe filled by a massless scalar written in
//   relational time: from H² = (8πG/3) ρ_φ and dφ/dt = π_φ / V, (dV/dφ / V)² = 9H² V² / π_φ²
//   = 12πG. Its solutions V ∝ e^{±√(12πG) φ} reach V = 0 at φ = ∓∞ in relational time (the
//   big bang, a ∝ t^{1/3} in proper time); the condensate's volume instead bounces where
//   V' = 0, at ρ_φ = ρ_c for one spin with E = 0 (OSW-E 2016, §4.2; Oriti, Sindoni and
//   Wilson-Ewing 2017, Class. Quantum Grav. 34, 04LT01).
//
// Units: Planck units, ħ = c = 1, with Newton's constant G a parameter (G = 1 by default);
// volumes in ℓ_P³, densities in Planck units. Interactions of the condensate (which modify the
// dynamics near the bounce and can produce recollapse or acceleration) are planned.
//
// Two computations are given: the closed form, and the integration of each mode by calx's
// Dormand–Prince pair outward from its bounce (stable, since both directions grow); the checks
// hold them against each other and against the conserved quantities and identities above.
//
// landau_condensate: the static mean field of a GFT potential μ|σ|² + λ|σ|^n (λ > 0, n > 2),
// condensed (σ ≠ 0 at the minimum) exactly when μ < 0, with |σ|^{n−2} = −2μ / (n λ) (Landau's
// theory of the transition; Pithis and Thürigen 2018, Phys. Rev. D 98, 126006, for its
// validity in GFT).

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <vector>

#include <calx/integrators/runge_kutta.hpp>

namespace physis::quantum_gravity {

// One isotropic mode of the condensate: spin j, charge Q, energy E, and the relational time
// of its bounce.
struct CondensateMode {
  double j = 0.5;           // spin (a half-integer >= 1/2)
  double Q = 1;             // conserved charge: the mode's share of the clock's momentum
  double E = 0;             // conserved energy
  double phi_bounce = 0;    // where u = ρ² is least
};

struct Universe {
  std::vector<CondensateMode> modes;
  double G = 1;          // Newton's constant, Planck units
  double m2 = -1;        // m² of every mode; negative: 3πG, the value that gives classical GR
  double v0 = 1;         // V_j = v0 j^{3/2}
  double mass2() const { return m2 > 0 ? m2 : 3 * std::numbers::pi * G; }
  double quantum_volume(double j) const { return v0 * std::pow(j, 1.5); }
  double pi_phi() const {
    double p = 0;
    for (const auto& m : modes) p += m.Q;
    return p;
  }
};

// The closed form of one mode: u = ρ², u', ρ, ρ' at φ.
struct ModeState {
  double u, du, rho, drho;
};

inline double mode_c(const CondensateMode& m, double mass2) { return std::sqrt(m.E * m.E + 4 * mass2 * m.Q * m.Q) / (4 * mass2); }

// u_j(φ) = −E/(2m²) + c e^{2m(φ−φ_j)} + c e^{−2m(φ−φ_j)}; ρ = √u and ρ' = u'/(2ρ).
inline ModeState mode_exact(const CondensateMode& m, double mass2, double phi) {
  const double mm = std::sqrt(mass2), c = mode_c(m, mass2);
  const double x = 2 * mm * (phi - m.phi_bounce);
  const double ep = std::exp(x), em = std::exp(-x);
  ModeState s;
  s.u = -m.E / (2 * mass2) + c * (ep + em);
  s.du = 2 * mm * c * (ep - em);
  s.rho = std::sqrt(std::max(s.u, 0.0));
  s.drho = s.rho > 0 ? s.du / (2 * s.rho) : std::numeric_limits<double>::infinity();
  return s;
}

// The least value of u: positive exactly when the mode bounces (Q ≠ 0, or E < 0).
inline double mode_u_min(const CondensateMode& m, double mass2) { return (-m.E + std::sqrt(m.E * m.E + 4 * mass2 * m.Q * m.Q)) / (2 * mass2); }

// The mode's equations as a first-order system in (ρ, ρ', θ).
struct ModeFlow {
  double Q, mass2;
  void operator()(const std::array<double, 3>& x, std::array<double, 3>& dx) const {
    dx[0] = x[1];
    dx[1] = Q * Q / (x[0] * x[0] * x[0]) + mass2 * x[0];
    dx[2] = Q / (x[0] * x[0]);
  }
};

// ρ, ρ', θ of a bouncing mode on the grid, by Dormand–Prince outward from its bounce
// (θ = 0 there). Grid points must be increasing.
inline std::vector<std::array<double, 3>> integrate_mode(const CondensateMode& m, double mass2, const std::vector<double>& grid,
                                                        calx::integrators::Tolerance tol = {1e-13, 1e-12}) {
  const double umin = mode_u_min(m, mass2);
  if (!(umin > 0)) throw std::invalid_argument("physis::quantum_gravity::integrate_mode: the mode does not bounce (Q = 0 and E >= 0): use mode_exact");
  const ModeFlow f{m.Q, mass2};
  std::vector<std::array<double, 3>> out(grid.size());
  const std::array<double, 3> start{std::sqrt(umin), 0.0, 0.0};
  const std::array<double, 3> scale{std::sqrt(umin), std::sqrt(umin * mass2), 1.0};
  // Forward from the bounce over the points after it, backward over those before. Backward,
  // the state (ρ, ρ', θ) keeps its derivatives with respect to φ and flows by −f in s = −φ.
  for (int dir : {+1, -1}) {
    std::array<double, 3> x = start;
    double at = m.phi_bounce, sub = 0;
    auto flow = [&](const std::array<double, 3>& y, std::array<double, 3>& dy) {
      f(y, dy);
      if (dir < 0)
        for (auto& v : dy) v = -v;
    };
    std::vector<std::size_t> order;
    for (std::size_t i = 0; i < grid.size(); i++)
      if ((dir > 0 && grid[i] >= m.phi_bounce) || (dir < 0 && grid[i] < m.phi_bounce)) order.push_back(i);
    if (dir < 0) std::reverse(order.begin(), order.end());
    for (std::size_t i : order) {
      const double h = std::fabs(grid[i] - at);
      if (h > 0) calx::integrators::dopri5<3>(flow, x, h, scale, tol, sub);
      at = grid[i];
      out[i] = x;
    }
  }
  return out;
}

// The universe on a grid of relational time.
struct CondensateCurves {
  std::vector<double> phi;
  std::vector<double> volume;            // V = Σ V_j ρ_j²
  std::vector<double> hubble;            // V' / 3V, the relational Hubble rate
  std::vector<double> hubble_squared;    // (V' / 3V)²
  std::vector<double> acceleration;      // V'' / V
  std::vector<double> density;           // ρ_φ = π_φ² / (2 V²)
  std::vector<double> classical_expanding, classical_contracting;  // Σ V_j C^±_j e^{±2mφ}: the classical branches
  std::vector<double> friedmann;         // (4πG/3)(1 − ρ_φ/ρ_c) + 4 Σ V_j E_j / (9 V): exact for one spin
  double pi_phi = 0;
  double critical_density = 0;           // ρ_c = m² / (2 V_j²) of the first mode
  bool bounces = false;
  double bounce_phi = 0, bounce_volume = 0, bounce_density = 0;  // where V is least on the grid (refined)
  double classical_rate = 0;             // 12πG, the classical (V'/V)²
};

// integrate: whether ρ_j come from Dormand–Prince (bouncing modes) or from the closed form.
inline CondensateCurves evolve(const Universe& U, double phi_min, double phi_max, std::size_t points, bool integrate = false) {
  if (U.modes.empty()) throw std::invalid_argument("physis::quantum_gravity::evolve: no modes");
  if (!(phi_max > phi_min) || points < 2) throw std::invalid_argument("physis::quantum_gravity::evolve: an increasing range of at least two points");
  const double m2 = U.mass2(), m = std::sqrt(m2);
  CondensateCurves C;
  C.pi_phi = U.pi_phi();
  C.classical_rate = 12 * std::numbers::pi * U.G;
  const double V1 = U.quantum_volume(U.modes[0].j);
  C.critical_density = m2 / (2 * V1 * V1);
  C.bounces = true;
  for (const auto& md : U.modes) C.bounces = C.bounces && mode_u_min(md, m2) > 0;
  for (std::size_t i = 0; i < points; i++) C.phi.push_back(phi_min + (phi_max - phi_min) * static_cast<double>(i) / static_cast<double>(points - 1));
  std::vector<std::vector<std::array<double, 3>>> integrated;
  if (integrate)
    for (const auto& md : U.modes) integrated.push_back(integrate_mode(md, m2, C.phi));
  double sumVE = 0;
  for (const auto& md : U.modes) sumVE += U.quantum_volume(md.j) * md.E;
  for (std::size_t i = 0; i < points; i++) {
    double V = 0, dV = 0, Vplus = 0, Vminus = 0;
    for (std::size_t k = 0; k < U.modes.size(); k++) {
      const auto& md = U.modes[k];
      const double Vj = U.quantum_volume(md.j);
      double rho, drho;
      if (integrate) rho = integrated[k][i][0], drho = integrated[k][i][1];
      else {
        const auto s = mode_exact(md, m2, C.phi[i]);
        rho = s.rho, drho = s.drho;
      }
      V += Vj * rho * rho;
      dV += 2 * Vj * rho * drho;
      const double c = mode_c(md, m2);
      Vplus += Vj * c * std::exp(2 * m * (C.phi[i] - md.phi_bounce));
      Vminus += Vj * c * std::exp(-2 * m * (C.phi[i] - md.phi_bounce));
    }
    const double H = dV / (3 * V);
    const double rho_phi = C.pi_phi * C.pi_phi / (2 * V * V);
    C.volume.push_back(V);
    C.hubble.push_back(H);
    C.hubble_squared.push_back(H * H);
    C.acceleration.push_back(4 * m2 + 2 * sumVE / V);
    C.density.push_back(rho_phi);
    C.classical_expanding.push_back(Vplus);
    C.classical_contracting.push_back(Vminus);
    C.friedmann.push_back(4 * std::numbers::pi * U.G / 3 * (1 - rho_phi / C.critical_density) + 4 * sumVE / (9 * V));
  }
  // The bounce: the least volume, refined by a parabola through the neighbours.
  const std::size_t k = static_cast<std::size_t>(std::min_element(C.volume.begin(), C.volume.end()) - C.volume.begin());
  C.bounce_phi = C.phi[k];
  C.bounce_volume = C.volume[k];
  if (k > 0 && k + 1 < points) {
    const double a = C.volume[k - 1], b = C.volume[k], c = C.volume[k + 1], h = C.phi[1] - C.phi[0];
    const double denom = a - 2 * b + c;
    if (denom > 0) {
      const double off = 0.5 * h * (a - c) / denom;
      C.bounce_phi = C.phi[k] + off;
      C.bounce_volume = b - (a - c) * (a - c) / (8 * denom);
    }
  }
  C.bounce_density = C.pi_phi * C.pi_phi / (2 * C.bounce_volume * C.bounce_volume);
  return C;
}

// The static mean field of μ|σ|² + λ|σ|^n.
struct LandauCondensate {
  bool condensed = false;
  double amplitude = 0;   // |σ| at the minimum
  double potential = 0;   // the potential there
};

inline LandauCondensate landau_condensate(double mu, double lambda, double n) {
  if (!(lambda > 0) || !(n > 2)) throw std::invalid_argument("physis::quantum_gravity::landau_condensate: λ > 0 and n > 2");
  LandauCondensate L;
  if (mu >= 0) return L;
  L.condensed = true;
  L.amplitude = std::pow(-2 * mu / (n * lambda), 1 / (n - 2));
  L.potential = mu * L.amplitude * L.amplitude + lambda * std::pow(L.amplitude, n);
  return L;
}

}  // namespace physis::quantum_gravity
