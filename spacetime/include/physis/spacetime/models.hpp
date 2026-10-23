// Spacetime models: exact solutions of Einstein's equations as metric fields in charts.
//
// Each model is a metric field of physis `geometry` (an object with n = 4 and a templated
// operator() giving g_ab(x)), so its curvature and geodesics come from `geometry`. A model
// knows no body: its mass, charge, spin, curvature scale and scale factor are members set by
// the caller (physis invariant 2). Units (physis decision 0005, §2.4): coordinates are
// lengths in metres, time as x^0 = ct; a mass M enters as its gravitational length GM/c², a
// charge Q as sqrt(G / (4π ε0)) Q / c², a spin as a = J / (M c) with M the gravitational
// length; curvature radii ℓ are lengths. Signature (-, +, +, +). Angles θ (colatitude) and
// φ are in radians.
//
// The models, with their charts and closed forms (Misner, Thorne and Wheeler 1973; Wald
// 1984; Griffiths and Podolský 2009, Exact Space-Times in Einstein's General Relativity):
//
// - Minkowski space in Cartesian (ct, x, y, z) and spherical (ct, r, θ, φ) coordinates.
// - Schwarzschild in Schwarzschild coordinates (ct, r, θ, φ):
//     ds² = -(1 - 2M/r) c²dt² + dr² / (1 - 2M/r) + r² dΩ²                  (MTW 31.1)
//   and in ingoing Eddington-Finkelstein coordinates (v, r, θ, φ), v = ct + r + 2M ln|r/2M - 1|:
//     ds² = -(1 - 2M/r) dv² + 2 dv dr + r² dΩ²,                              (MTW Box 31.2)
//   regular across the horizon r = 2M. Photon sphere 3M, innermost stable circular orbit 6M;
//   Kretschmann scalar 48 M² / r⁶.
// - Reissner-Nordström: f = 1 - 2M/r + Q²/r², horizons r± = M ± sqrt(M² - Q²); Kretschmann
//   8 (6 M² r² - 12 M Q² r + 7 Q⁴) / r⁸ (Henry 2000, ApJ 535, 350).
// - Kerr in Boyer-Lindquist coordinates (ct, r, θ, φ), Σ = r² + a² cos²θ, Δ = r² - 2Mr + a²:
//     ds² = -(1 - 2Mr/Σ) c²dt² - (4 M a r sin²θ / Σ) c dt dφ + (Σ/Δ) dr² + Σ dθ²
//           + (r² + a² + 2 M a² r sin²θ / Σ) sin²θ dφ²                     (MTW 33.2)
//   horizons r± = M ± sqrt(M² - a²); the innermost stable circular orbits of Bardeen, Press and
//   Teukolsky 1972 (ApJ 178, 347, eq. 2.21); Kretschmann
//   48 M² (r² - a²c²) ((r² + a²c²)² - 16 r² a² c²) / (r² + a²c²)⁶, c = cos θ (Henry 2000).
// - FLRW, homogeneous and isotropic, with spatial curvature K (1/m²) and a dimensionless
//   scale factor a(ct) given as a template callable, in coordinates (ct, r, θ, φ):
//     ds² = -c²dt² + a² (dr² / (1 - K r²) + r² dΩ²)                        (MTW 27.26; Carroll 8.9)
//   Its Einstein tensor gives Friedmann's equations: G_00 = 3 (ȧ² + K) / a², G^i_i =
//   -(2 ä / a + (ȧ² + K) / a²), dots derivatives in ct (Carroll 2004, eqs. 8.44-8.46).
// - de Sitter space of radius ℓ (constant curvature, R = 12 / ℓ²): the static patch
//   (ct, r, θ, φ), ds² = -(1 - r²/ℓ²) c²dt² + dr² / (1 - r²/ℓ²) + r² dΩ², horizon at r = ℓ;
//   and the flat slicing (ct, x, y, z), ds² = -c²dt² + e^{2ct/ℓ} (dx² + dy² + dz²).
// - Anti-de Sitter space of radius ℓ (R = -12 / ℓ²) in global coordinates (ct, r, θ, φ):
//   ds² = -(1 + r²/ℓ²) c²dt² + dr² / (1 + r²/ℓ²) + r² dΩ².
//
// Each chart covers part of its manifold, and its components are singular where it ends
// (the horizon of Schwarzschild coordinates, the axis θ = 0 of every spherical chart).

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include "minkowski.hpp"

namespace physis::spacetime {

template <class T>
using Matrix4 = std::array<std::array<T, 4>, 4>;

// The gravitational constant and the speed of light (CODATA 2018), until physis `units`
// carries its versioned constants.
inline constexpr double gravitational_constant = 6.67430e-11;  // m³ / (kg s²)
inline constexpr double speed_of_light = 299792458.0;           // m / s, exact

// GM / c², the length by which a mass enters a spacetime model (m), for a mass in kg.
inline double gravitational_length(double mass_kg) { return gravitational_constant * mass_kg / (speed_of_light * speed_of_light); }

using MinkowskiCartesian = Minkowski<4>;

// Minkowski space in spherical coordinates (ct, r, θ, φ).
struct MinkowskiSpherical {
  static constexpr std::size_t n = 4;
  template <class T>
  Matrix4<T> operator()(const std::array<T, 4>& x) const {
    using std::sin;
    Matrix4<T> g{};
    const T& r = x[1];
    const T s = sin(x[2]);
    g[0][0] = T(-1.0);
    g[1][1] = T(1.0);
    g[2][2] = r * r;
    g[3][3] = r * r * s * s;
    return g;
  }
};

// A static, spherically symmetric metric -f dt² + dr² / f + r² dΩ², for the models below.
template <class T>
inline Matrix4<T> static_spherical(const std::array<T, 4>& x, const T& f) {
  using std::sin;
  Matrix4<T> g{};
  const T& r = x[1];
  const T s = sin(x[2]);
  g[0][0] = -f;
  g[1][1] = 1.0 / f;
  g[2][2] = r * r;
  g[3][3] = r * r * s * s;
  return g;
}

struct Schwarzschild {
  static constexpr std::size_t n = 4;
  double M = 0;  // gravitational length GM/c², m
  template <class T>
  Matrix4<T> operator()(const std::array<T, 4>& x) const { return static_spherical(x, T(1.0 - 2.0 * M / x[1])); }
  double horizon() const { return 2 * M; }
  double photon_sphere() const { return 3 * M; }
  double isco() const { return 6 * M; }
  // R_abcd R^abcd at radius r (1/m⁴).
  double kretschmann(double r) const { return 48 * M * M / std::pow(r, 6); }
};

// Schwarzschild in ingoing Eddington-Finkelstein coordinates (v, r, θ, φ).
struct SchwarzschildEddingtonFinkelstein {
  static constexpr std::size_t n = 4;
  double M = 0;  // m
  template <class T>
  Matrix4<T> operator()(const std::array<T, 4>& x) const {
    using std::sin;
    Matrix4<T> g{};
    const T& r = x[1];
    const T s = sin(x[2]);
    g[0][0] = -(1.0 - 2.0 * M / r);
    g[0][1] = g[1][0] = T(1.0);
    g[2][2] = r * r;
    g[3][3] = r * r * s * s;
    return g;
  }
  // The advanced time v of the Schwarzschild event (ct, r), outside the horizon.
  double advanced_time(double ct, double r) const { return ct + r + 2 * M * std::log(std::fabs(r / (2 * M) - 1)); }
};

struct ReissnerNordstrom {
  static constexpr std::size_t n = 4;
  double M = 0;  // m
  double Q = 0;  // charge as a length, sqrt(G / (4π ε0)) Q / c², m
  template <class T>
  Matrix4<T> operator()(const std::array<T, 4>& x) const {
    const T& r = x[1];
    return static_spherical(x, T(1.0 - 2.0 * M / r + Q * Q / (r * r)));
  }
  // The outer and inner horizons r± (none when Q > M).
  double outer_horizon() const { return M + std::sqrt(M * M - Q * Q); }
  double inner_horizon() const { return M - std::sqrt(M * M - Q * Q); }
  double kretschmann(double r) const { return 8 * (6 * M * M * r * r - 12 * M * Q * Q * r + 7 * Q * Q * Q * Q) / std::pow(r, 8); }
};

struct Kerr {
  static constexpr std::size_t n = 4;
  double M = 0;  // m
  double a = 0;  // spin J / (M c), m; |a| <= M for a black hole
  template <class T>
  Matrix4<T> operator()(const std::array<T, 4>& x) const {
    using std::cos;
    using std::sin;
    const T& r = x[1];
    const T s = sin(x[2]), c = cos(x[2]);
    const T Sigma = r * r + a * a * c * c;
    const T Delta = r * r - 2.0 * M * r + a * a;
    Matrix4<T> g{};
    g[0][0] = -(1.0 - 2.0 * M * r / Sigma);
    g[0][3] = g[3][0] = -2.0 * M * a * r * s * s / Sigma;
    g[1][1] = Sigma / Delta;
    g[2][2] = Sigma;
    g[3][3] = (r * r + a * a + 2.0 * M * a * a * r * s * s / Sigma) * s * s;
    return g;
  }
  double outer_horizon() const { return M + std::sqrt(M * M - a * a); }
  double inner_horizon() const { return M - std::sqrt(M * M - a * a); }
  // The innermost stable circular orbit, prograde (a > 0 co-rotating) or retrograde.
  double isco(bool prograde = true) const {
    const double chi = a / M;
    const double z1 = 1 + std::cbrt(1 - chi * chi) * (std::cbrt(1 + chi) + std::cbrt(1 - chi));
    const double z2 = std::sqrt(3 * chi * chi + z1 * z1);
    const double root = std::sqrt((3 - z1) * (3 + z1 + 2 * z2));
    return M * (3 + z2 + (prograde ? -root : root));
  }
  double kretschmann(double r, double theta) const {
    const double c = std::cos(theta), r2 = r * r, ac2 = a * a * c * c;
    return 48 * M * M * (r2 - ac2) * ((r2 + ac2) * (r2 + ac2) - 16 * r2 * ac2) / std::pow(r2 + ac2, 6);
  }
};

// Scale factors for FLRW, functions of x^0 = ct.
struct PowerLaw {
  double t0 = 1;      // the ct at which a = 1, m
  double p = 2.0 / 3;  // exponent: 2/3 for dust, 1/2 for radiation
  template <class T>
  T operator()(const T& ct) const { using std::pow; return pow(ct / t0, p); }
};
struct Exponential {
  double H = 1;  // Hubble rate H / c, 1/m
  template <class T>
  T operator()(const T& ct) const { using std::exp; return exp(H * ct); }
};

template <class ScaleFactor = PowerLaw>
struct FLRW {
  static constexpr std::size_t n = 4;
  double K = 0;    // spatial curvature, 1/m² (k / R0²)
  ScaleFactor a{};
  template <class T>
  Matrix4<T> operator()(const std::array<T, 4>& x) const {
    using std::sin;
    const T& r = x[1];
    const T A = a(x[0]);
    const T s = sin(x[2]);
    Matrix4<T> g{};
    g[0][0] = T(-1.0);
    g[1][1] = A * A / (1.0 - K * r * r);
    g[2][2] = A * A * r * r;
    g[3][3] = A * A * r * r * s * s;
    return g;
  }
};

struct DeSitterStatic {
  static constexpr std::size_t n = 4;
  double ell = 1;  // radius of curvature, m
  template <class T>
  Matrix4<T> operator()(const std::array<T, 4>& x) const { return static_spherical(x, T(1.0 - x[1] * x[1] / (ell * ell))); }
  double horizon() const { return ell; }
  double scalar_curvature() const { return 12 / (ell * ell); }
};

struct DeSitterFlat {
  static constexpr std::size_t n = 4;
  double ell = 1;  // m
  template <class T>
  Matrix4<T> operator()(const std::array<T, 4>& x) const {
    using std::exp;
    Matrix4<T> g{};
    const T e = exp(2.0 * x[0] / ell);
    g[0][0] = T(-1.0);
    g[1][1] = e;
    g[2][2] = e;
    g[3][3] = e;
    return g;
  }
  double scalar_curvature() const { return 12 / (ell * ell); }
};

struct AntiDeSitterGlobal {
  static constexpr std::size_t n = 4;
  double ell = 1;  // m
  template <class T>
  Matrix4<T> operator()(const std::array<T, 4>& x) const { return static_spherical(x, T(1.0 + x[1] * x[1] / (ell * ell))); }
  double scalar_curvature() const { return -12 / (ell * ell); }
};

}  // namespace physis::spacetime
