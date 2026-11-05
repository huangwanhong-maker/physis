// Capillarity and slow flow in tubes: the Young–Laplace pressure of a meniscus, Hagen and Poiseuille's
// conductance, and Washburn's imbibition (decision 0022).
//
// A meniscus in a tube of radius r, of surface tension γ and contact angle θ, holds a pressure jump
// 2 γ cos θ / r (Young 1805; Laplace 1806). Steady laminar flow through a tube of radius r and length L at
// viscosity μ carries Q = G Δp with G = π r⁴ / (8 μ L) (Hagen 1839; Poiseuille 1840). A liquid drawn into a
// horizontal tube by its meniscus alone, the column's inertia and the displaced air's viscosity neglected,
// advances as L² = γ r cos θ t / (2 μ) (Washburn 1921, Phys. Rev. 17, 273, eq. 6 without gravity), and so
// fills a tube of length L in 2 μ L² / (γ r cos θ). Units SI. Checks: tests/check_capillarity.cpp.
#pragma once

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace physis::continuum {

// Pa.
inline double capillary_pressure(double radius, double surface_tension, double contact_angle) {
  return 2 * surface_tension * std::cos(contact_angle) / radius;
}

// m³/(s Pa).
inline double poiseuille_conductance(double radius, double length, double viscosity) {
  if (!(radius > 0) || !(length > 0) || !(viscosity > 0)) throw std::invalid_argument("physis::continuum::poiseuille_conductance: positive radius, length and viscosity");
  return std::numbers::pi * radius * radius * radius * radius / (8 * viscosity * length);
}

// m: the length imbibed after t s.
inline double washburn_length(double t, double radius, double surface_tension, double contact_angle, double viscosity) {
  return std::sqrt(surface_tension * radius * std::cos(contact_angle) * t / (2 * viscosity));
}

// s: the time to imbibe a length L.
inline double washburn_time(double length, double radius, double surface_tension, double contact_angle, double viscosity) {
  return 2 * viscosity * length * length / (surface_tension * radius * std::cos(contact_angle));
}

}  // namespace physis::continuum
