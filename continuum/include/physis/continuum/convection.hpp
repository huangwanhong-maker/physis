// The onset of Rayleigh–Bénard convection, and Lorenz's truncation of it (decision 0022).
//
// A horizontal layer of a Boussinesq fluid, depth d, heated from below, is stable until its Rayleigh number
// Ra = g α ΔT d³ / (ν κ) passes a critical value. A normal mode W(z) e^{i a·x} (a the horizontal wavenumber
// times d, z ∈ [0, 1]) is marginal (neutral, with exchange of stabilities) when
//
//   (D² − a²)³ W = − Ra a² W,   D = d/dz          (Chandrasekhar 1961, Hydrodynamic and Hydromagnetic
//                                                   Stability, ch. II),
//
// with W = 0 and (D² − a²)² W = 0 (no vertical velocity, no temperature perturbation) at each plate, and
// DW = 0 at a rigid plate (no slip) or D²W = 0 at a free one (no stress). The equation has constant
// coefficients, so the state (W, DW, …, D⁵W) is carried from z = 0 to z = 1 exactly by the exponential of
// its 6 × 6 companion matrix (calx's expm). The three conditions at the bottom leave a three-dimensional
// family of starts; the three at the top make a 3 × 3 matrix whose determinant vanishes at the marginal
// Ra(a). The smallest positive root is found by a geometric scan and Brent's method, and its least value
// over a by Brent's minimiser.
//
// Values (Chandrasekhar 1961, ch. II): free–free, Ra(a) = (π² + a²)³ / a² exactly, least
// 27π⁴/4 = 657.511 at a = π/√2; rigid–rigid, 1707.762 at a = 3.117; rigid–free, 1100.65 at a = 2.682.
//
// Lorenz's truncation. Keeping the three modes of Saltzman (1962) between free plates, Lorenz (1963, J.
// Atmos. Sci. 20, 130) wrote X' = σ(Y − X), Y' = rX − XZ − Y, Z' = XY − bZ, with σ = ν/κ the
// Prandtl number, r = Ra / Ra_c(a), Ra_c(a) = π⁴(1 + a_L²)³/a_L², b = 4/(1 + a_L²), and the dimensionless
// time τ = π²(1 + a_L²) κ t / d², where a_L = a/π is his ratio of depth to half the rolls' wavelength. At
// a_L² = ½, the critical rolls of free plates, b = 8/3. X is the rolls' circulation, Y and Z the
// temperature's departures; in this scaling the Nusselt number, the heat carried relative to conduction, is
// 1 + 2Z/r, steady at 1 + 2(r − 1)/r on the convecting equilibria Z = r − 1 (1 < r below chaos). The truncation is a
// model of convection near onset between free plates; it is not the layer's flow at r = 28, which is
// three-dimensional turbulence. Checks: tests/check_convection.cpp.
#pragma once

#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>

#include <calx/linalg/dense.hpp>
#include <calx/linalg/expm.hpp>
#include <calx/optimisation/minimise.hpp>
#include <calx/solvers/roots.hpp>

namespace physis::continuum {

enum class Plate { rigid, free };

namespace convection_detail {

// The three boundary functionals at a plate, applied to the state (W, …, D⁵W).
inline std::array<double, 3> conditions(Plate p, double a, const double* y) {
  const double a2 = a * a;
  if (p == Plate::rigid) return {y[0], y[1], y[4] - 2 * a2 * y[2]};  // W, DW, (D² − a²)²W with W = 0
  return {y[0], y[2], y[4]};                                         // W, D²W, D⁴W
}

// A basis of the starts at the bottom plate that meet its three conditions.
inline std::array<std::array<double, 6>, 3> starts(Plate p, double a) {
  if (p == Plate::rigid) return {{{0, 0, 1, 0, 2 * a * a, 0}, {0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 1}}};
  return {{{0, 1, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 1}}};
}

inline double determinant3(const std::array<std::array<double, 3>, 3>& m) {
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
         m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

}  // namespace convection_detail

// The determinant whose zeros in Ra are the marginal states of wavenumber a, scaled to be of order one.
inline double marginal_determinant(double rayleigh, double a, Plate bottom, Plate top) {
  const double a2 = a * a;
  calx::linalg::Matrix M(6);
  for (int i = 0; i < 5; i++) M(i, i + 1) = 1;
  // D⁶W = 3a² D⁴W − 3a⁴ D²W + (a⁶ − Ra a²) W.
  M(5, 0) = a2 * a2 * a2 - rayleigh * a2, M(5, 2) = -3 * a2 * a2, M(5, 4) = 3 * a2;
  const calx::linalg::Matrix P = calx::linalg::expm(M);
  const auto s = convection_detail::starts(bottom, a);
  std::array<std::array<double, 3>, 3> B{};
  for (int c = 0; c < 3; c++) {
    double y[6];
    double scale = 0;
    for (int i = 0; i < 6; i++) {
      y[i] = 0;
      for (int j = 0; j < 6; j++) y[i] += P(i, j) * s[c][j];
      scale = std::max(scale, std::fabs(y[i]));
    }
    const auto f = convection_detail::conditions(top, a, y);
    for (int r = 0; r < 3; r++) B[r][c] = f[r] / scale;
  }
  return convection_detail::determinant3(B);
}

// The marginal Rayleigh number of wavenumber a (the smallest positive root), by a scan from 10 and Brent.
inline double marginal_rayleigh(double a, Plate bottom, Plate top) {
  if (!(a > 0)) throw std::invalid_argument("physis::continuum::marginal_rayleigh: a positive wavenumber");
  double lo = 10, flo = marginal_determinant(lo, a, bottom, top);
  for (int k = 0; k < 400; k++) {
    const double hi = lo * 1.04, fhi = marginal_determinant(hi, a, bottom, top);
    if ((flo < 0) != (fhi < 0)) {
      const auto root = calx::solvers::brent([&](double r) { return marginal_determinant(r, a, bottom, top); }, lo, hi, 1e-12);
      return root.x;
    }
    lo = hi, flo = fhi;
  }
  throw std::runtime_error("physis::continuum::marginal_rayleigh: no marginal state below Ra ≈ 1.6e8");
}

// Free–free plates in closed form.
inline double marginal_rayleigh_free(double a) { const double s = std::numbers::pi * std::numbers::pi + a * a; return s * s * s / (a * a); }

struct Onset {
  double rayleigh = 0, wavenumber = 0;
};

// The critical Rayleigh number and wavenumber: the least marginal Ra over a.
inline Onset critical_rayleigh(Plate bottom, Plate top) {
  const auto m = calx::optimisation::brent_minimise([&](double a) { return marginal_rayleigh(a, bottom, top); }, 1.0, 6.0, 1e-9);
  return {m.f, m.x};
}

// ---- Lorenz's truncation ----------------------------------------------------------------------------------------------

struct LorenzParameters {
  double sigma = 0, r = 0, b = 0;
  double time_scale = 0;  // dτ/dt, 1/s
};

// From the layer's Rayleigh and Prandtl numbers, the rolls' wavenumber a (times d; π/√2 for the critical rolls
// of free plates), and, for the time scale, the diffusivity κ (m²/s) and depth d (m).
inline LorenzParameters lorenz_parameters(double rayleigh, double prandtl, double a = std::numbers::pi / std::numbers::sqrt2, double kappa = 0, double depth = 0) {
  const double aL = a / std::numbers::pi, s = 1 + aL * aL, pi2 = std::numbers::pi * std::numbers::pi;
  LorenzParameters p;
  p.sigma = prandtl;
  p.r = rayleigh / (pi2 * pi2 * s * s * s / (aL * aL));
  p.b = 4 / s;
  p.time_scale = (kappa > 0 && depth > 0) ? pi2 * s * kappa / (depth * depth) : 0;
  return p;
}

// The Nusselt number of a state of Lorenz's model.
inline double lorenz_nusselt(double Z, double r) { return 1 + 2 * Z / r; }

}  // namespace physis::continuum
