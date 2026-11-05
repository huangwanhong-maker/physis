// Checks of capillarity and slow flow in tubes (capillarity.hpp; decision 0022).
//
// - Young–Laplace: a water meniscus (γ = 0.0728 N/m, θ = 0) in a tube of radius 0.5 mm holds 291.2 Pa.
// - Hagen–Poiseuille: the conductance scales as r⁴ and as 1/L, 1/μ.
// - Washburn: L ∝ √t; the length at the time to fill a length is that length; water fills 10 cm of a 0.1 mm
//   capillary in 2 μ L² / (γ r) = 2.747 s.

#include <physis/continuum/capillarity.hpp>

#include <calx/verification/check.hpp>

#include <cmath>

using namespace physis::continuum;

int main() {
  calx::verification::Suite s;
  const double gamma = 0.0728, mu = 1.0e-3;
  s.check("Young–Laplace: 2γ/r for water in a 0.5 mm tube (Pa)", capillary_pressure(5e-4, gamma, 0) - 291.2, 1e-9, "Pa");
  const double g1 = poiseuille_conductance(1e-4, 0.1, mu), g2 = poiseuille_conductance(2e-4, 0.1, mu);
  s.check("Poiseuille: twice the radius, sixteen times the conductance", g2 / g1 - 16, 1e-12);
  s.check("Poiseuille: twice the length, half", poiseuille_conductance(1e-4, 0.2, mu) / g1 - 0.5, 1e-15);
  s.check("Washburn: four times the time, twice the length", washburn_length(4.0, 1e-4, gamma, 0.3, mu) / washburn_length(1.0, 1e-4, gamma, 0.3, mu) - 2, 1e-15);
  const double t = washburn_time(0.1, 1e-4, gamma, 0.3, mu);
  s.check("Washburn: the length at the time to fill it (m)", washburn_length(t, 1e-4, gamma, 0.3, mu) - 0.1, 1e-15, "m");
  s.check("Washburn: water fills 10 cm of a 0.1 mm capillary in 2.747 s", washburn_time(0.1, 1e-4, gamma, 0, mu) - 2 * mu * 0.01 / (gamma * 1e-4), 1e-15, "s");
  s.check("Washburn: that time, to four figures (s)", washburn_time(0.1, 1e-4, gamma, 0, mu) - 2.747, 5e-4, "s");
  return s.done();
}
