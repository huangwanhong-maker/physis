// Checks of the onset of convection and Lorenz's truncation (convection.hpp; decision 0022), against
// Chandrasekhar (1961, ch. II) and Lorenz (1963).
//
// - Free plates: the shooting's marginal Ra(a) equals (π² + a²)³/a² (relative 1e−8 at a = 1.5, 2.2, 3.5);
//   the least is 27π⁴/4 = 657.511 at a = π/√2 (relative 1e−7 and 1e−4).
// - Rigid plates: 1707.762 at a = 3.117 (within 0.01 and 0.001).
// - Rigid below, free above: 1100.65 at a = 2.682 (within 0.01 and 0.001).
// - Lorenz's parameters: at a = π/√2, b = 8/3 and r = Ra / (27π⁴/4); the time scale π²(3/2)κ/d²; the Nusselt
//   number on the convecting equilibria is 1 + 2(r − 1)/r.

#include <physis/continuum/convection.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <numbers>

using namespace physis::continuum;

int main() {
  calx::verification::Suite s;
  const double pi = std::numbers::pi;
  for (double a : {1.5, 2.2, 3.5}) {
    const double r = marginal_rayleigh(a, Plate::free, Plate::free);
    s.check("free plates: the marginal Ra(a) by shooting equals (π² + a²)³/a² (relative)", r / marginal_rayleigh_free(a) - 1, 1e-8);
  }
  const auto ff = critical_rayleigh(Plate::free, Plate::free);
  s.check("free plates: Ra_c = 27π⁴/4 (relative)", ff.rayleigh / (27 * std::pow(pi, 4) / 4) - 1, 1e-7);
  s.check("free plates: a_c = π/√2 (relative)", ff.wavenumber / (pi / std::numbers::sqrt2) - 1, 1e-4);
  const auto rr = critical_rayleigh(Plate::rigid, Plate::rigid);
  s.check("rigid plates: Ra_c = 1707.762 (Chandrasekhar 1961)", rr.rayleigh - 1707.762, 0.01);
  s.check("rigid plates: a_c = 3.117", rr.wavenumber - 3.117, 0.001);
  const auto rf = critical_rayleigh(Plate::rigid, Plate::free);
  s.check("rigid below, free above: Ra_c = 1100.65", rf.rayleigh - 1100.65, 0.01);
  s.check("rigid below, free above: a_c = 2.682", rf.wavenumber - 2.682, 0.001);
  const auto fr = critical_rayleigh(Plate::free, Plate::rigid);
  s.check("the layer turned over has the same onset", fr.rayleigh - rf.rayleigh, 1e-6);

  const double Rc = 27 * std::pow(pi, 4) / 4;
  const auto L = lorenz_parameters(28 * Rc, 10, pi / std::numbers::sqrt2, 1.43e-7, 0.01);
  s.check("Lorenz: b = 8/3 for the critical rolls of free plates", L.b - 8.0 / 3, 1e-14);
  s.check("Lorenz: r = Ra / Ra_c", L.r - 28, 1e-12);
  s.check("Lorenz: σ = Pr", L.sigma - 10, 0);
  s.check("Lorenz: dτ/dt = π² (3/2) κ / d² (relative)", L.time_scale / (pi * pi * 1.5 * 1.43e-7 / 1e-4) - 1, 1e-14);
  s.check("Lorenz: the Nusselt number on the convecting equilibria, 1 + 2(r − 1)/r", lorenz_nusselt(5 - 1, 5) - (1 + 2.0 * 4 / 5), 1e-15);
  return s.done();
}
