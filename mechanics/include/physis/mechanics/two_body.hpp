// The two-body problem: a point moving about a centre of attraction μ/r², its orbit described by
// classical elements, propagated in time, and the frames attached to it.
//
// Elements. An orbit is a conic of semi-latus rectum p and eccentricity e in a plane of inclination
// i and right ascension of the ascending node Ω, with argument of periapsis ω, the body at true
// anomaly ν; the semi-major axis is a = p/(1 − e²) (negative for a hyperbola, infinite for a
// parabola). From a state (r, v) they follow from the angular momentum h = r × v, the node vector
// n = ẑ × h and the eccentricity vector e = ((v² − μ/r) r − (r·v) v)/μ (Vallado 2013, Fundamentals
// of Astrodynamics and Applications, 4th ed., Algorithm 9, RV2COE); the state from the elements
// through the perifocal frame (Algorithm 10, COE2RV). Where an angle is undefined (a circular orbit
// has no periapsis, an equatorial one no node) the element taken in its place is named: the
// argument of latitude u = ω + ν (circular inclined), the true longitude of periapsis ϖ = Ω + ω
// (elliptical equatorial), the true longitude λ = Ω + ω + ν (circular equatorial), with the
// undefined angles set to zero; `Elements::kind` says which.
//
// Kepler's equation. The mean anomaly M = n(t − τ) is E − e sin E for an ellipse and e sinh H − H for
// a hyperbola; it is solved by Newton's method from Danby's start (Danby 1988, Fundamentals of
// Celestial Mechanics, §6.6), and a parabola's by Barker's equation in closed form.
//
// Propagation. A state is carried over a time Δt on any conic by universal variables: the
// universal anomaly χ solves √μ Δt = χ³ c₃(ψ) + (r₀·v₀/√μ) χ² c₂(ψ) + r₀ χ (1 − ψ c₃(ψ)), ψ = χ²/a,
// by Newton's method, and the Lagrange coefficients f, g, ḟ, ġ give the state (Vallado 2013,
// Algorithm 8, KEPLER; Battin 1999, §4.5), with the Stumpff functions c₂, c₃ by their series near
// ψ = 0.
//
// Frames of an orbit, each a right-handed orthonormal triad moving with the body (its origin the
// body, its axes as columns of the rotation to the inertial frame), with the rate at which it turns:
//   RTN (or RSW, Vallado §3.3): R along r, N along h, T = N × R;
//   LVLH (local vertical, local horizontal, the convention in which NASA flies the ISS): z toward
//     the centre (−R), y against the orbit normal (−N), x = y × z (along the velocity on a circle);
//   VNC (Vallado §3.3): V along v, N along h, C = V × N;
//   perifocal PQW: P toward periapsis (along e), W along h, Q = W × P; fixed while the orbit is.
// Their rates follow from differentiating the triad's defining vectors: r' = v, v' = a (the body's
// acceleration, the central one −μr/r³ unless another is given), h' = r × a, and for a unit vector
// û = w/|w|, û' = (w' − û(û·w'))/|w|.

#pragma once

#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string_view>

#include <physis/geometry/kinematics.hpp>
#include <physis/geometry/rotation.hpp>

namespace physis::mechanics {

using geometry::Mat3;
using geometry::StateTransform;
using geometry::Vec3;
using namespace geometry::operators;

inline double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
inline Vec3 cross(const Vec3& a, const Vec3& b) { return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}; }
inline double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

// ---- elements ---------------------------------------------------------------------------------------

enum class OrbitKind { elliptical_inclined, circular_inclined, elliptical_equatorial, circular_equatorial };

inline std::string_view name(OrbitKind k) {
  switch (k) {
    case OrbitKind::elliptical_inclined: return "elliptical inclined";
    case OrbitKind::circular_inclined: return "circular inclined (ω and ν undefined: u = ω + ν given in ν)";
    case OrbitKind::elliptical_equatorial: return "elliptical equatorial (Ω undefined: ϖ = Ω + ω given in ω)";
    default: return "circular equatorial (Ω and ω undefined: λ = Ω + ω + ν given in ν)";
  }
}

struct Elements {
  double p = 0;     // semi-latus rectum (m)
  double e = 0;     // eccentricity
  double i = 0;     // inclination (rad)
  double raan = 0;  // right ascension of the ascending node Ω (rad)
  double argp = 0;  // argument of periapsis ω (rad)
  double nu = 0;    // true anomaly ν (rad)
  OrbitKind kind = OrbitKind::elliptical_inclined;
  // The semi-major axis (m): negative for a hyperbola, infinite for a parabola.
  double a() const { return std::fabs(1 - e) < 1e-15 ? std::numeric_limits<double>::infinity() : p / (1 - e * e); }
};

inline double wrap_2pi(double x) {
  x = std::fmod(x, 2 * std::numbers::pi);
  return x < 0 ? x + 2 * std::numbers::pi : x;
}

// The elements of the state (r, v) about μ (RV2COE). Tolerances decide when an orbit counts as
// circular or equatorial.
inline Elements elements(const Vec3& r, const Vec3& v, double mu, double small = 1e-11) {
  const Vec3 h = cross(r, v);
  const double hn = norm(h), rn = norm(r), v2 = dot(v, v);
  if (!(hn > 0)) throw std::invalid_argument("physis: a rectilinear state has no orbital elements");
  const Vec3 n{-h[1], h[0], 0};  // ẑ × h
  const double nn = norm(n);
  const Vec3 ev = (1 / mu) * ((v2 - mu / rn) * r - dot(r, v) * v);
  Elements o;
  o.e = norm(ev);
  o.p = hn * hn / mu;
  o.i = std::acos(std::max(-1.0, std::min(1.0, h[2] / hn)));
  const bool circular = o.e < small, equatorial = nn < small * hn;
  if (!equatorial) {
    o.raan = std::atan2(n[1], n[0]);
    o.raan = wrap_2pi(o.raan);
  }
  if (!circular && !equatorial) {
    o.kind = OrbitKind::elliptical_inclined;
    o.argp = std::atan2(dot(cross(n, ev), h) / hn, dot(n, ev));
    o.nu = std::atan2(dot(cross(ev, r), h) / hn, dot(ev, r));
  } else if (circular && !equatorial) {
    o.kind = OrbitKind::circular_inclined;
    o.nu = std::atan2(dot(cross(n, r), h) / hn, dot(n, r));  // the argument of latitude u
  } else if (!circular && equatorial) {
    o.kind = OrbitKind::elliptical_equatorial;
    o.argp = std::atan2(ev[1], ev[0]) * (h[2] < 0 ? -1 : 1);  // the true longitude of periapsis ϖ
    o.nu = std::atan2(dot(cross(ev, r), h) / hn, dot(ev, r));
  } else {
    o.kind = OrbitKind::circular_equatorial;
    o.nu = std::atan2(r[1], r[0]) * (h[2] < 0 ? -1 : 1);  // the true longitude λ
  }
  o.argp = wrap_2pi(o.argp);
  o.nu = wrap_2pi(o.nu);
  return o;
}

// The rotation from the perifocal frame (P, Q, W) to the inertial one: R₃(−Ω) R₁(−i) R₃(−ω).
inline Mat3 perifocal_to_inertial(double raan, double i, double argp) {
  const double cO = std::cos(raan), sO = std::sin(raan), ci = std::cos(i), si = std::sin(i), cw = std::cos(argp), sw = std::sin(argp);
  return Mat3{{{cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si},
               {sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si},
               {sw * si, cw * si, ci}}};
}

// The state of the elements about μ (COE2RV): position (m) and velocity (m/s).
inline std::array<Vec3, 2> state(const Elements& o, double mu) {
  const double c = std::cos(o.nu), s = std::sin(o.nu), r = o.p / (1 + o.e * c), k = std::sqrt(mu / o.p);
  const Vec3 rp{r * c, r * s, 0}, vp{-k * s, k * (o.e + c), 0};
  double raan = o.raan, argp = o.argp;
  if (o.kind == OrbitKind::circular_inclined) argp = 0;
  if (o.kind == OrbitKind::elliptical_equatorial) raan = 0;
  if (o.kind == OrbitKind::circular_equatorial) raan = argp = 0;
  const Mat3 R = perifocal_to_inertial(raan, o.i, argp);
  return {R * rp, R * vp};
}

// ---- Kepler's equation ----------------------------------------------------------------------------------

// The eccentric anomaly E of the mean anomaly M on an ellipse (0 ≤ e < 1): E − e sin E = M.
inline double eccentric_anomaly(double M, double e) {
  if (e < 0 || e >= 1) throw std::invalid_argument("physis: an ellipse has 0 <= e < 1");
  const double Mw = std::remainder(M, 2 * std::numbers::pi);
  double E = Mw + 0.85 * e * (Mw >= 0 ? 1 : -1);  // Danby's start
  for (int k = 0; k < 50; k++) {
    const double f = E - e * std::sin(E) - Mw, fp = 1 - e * std::cos(E), fpp = e * std::sin(E);
    // Halley's step (third order): it converges from Danby's start for every e < 1.
    const double dE = -f / (fp - 0.5 * f * fpp / fp);
    E += dE;
    if (std::fabs(dE) < 1e-15 * std::max(1.0, std::fabs(E))) break;
  }
  return E + (M - Mw);
}

// The hyperbolic anomaly H of the mean anomaly M on a hyperbola (e > 1): e sinh H − H = M.
inline double hyperbolic_anomaly(double M, double e) {
  if (e <= 1) throw std::invalid_argument("physis: a hyperbola has e > 1");
  double H = std::asinh(M / e);  // a start that brackets the root from the side of convergence
  if (std::fabs(M) > 6 * e) H = std::copysign(std::log(2 * std::fabs(M) / e), M);
  for (int k = 0; k < 100; k++) {
    const double f = e * std::sinh(H) - H - M, fp = e * std::cosh(H) - 1;
    const double dH = -f / fp;
    H += dH;
    if (std::fabs(dH) < 1e-15 * std::max(1.0, std::fabs(H))) break;
  }
  return H;
}

// True anomaly from the eccentric (or hyperbolic) anomaly and back.
inline double true_from_eccentric(double E, double e) { return 2 * std::atan2(std::sqrt(1 + e) * std::sin(E / 2), std::sqrt(1 - e) * std::cos(E / 2)); }
inline double eccentric_from_true(double nu, double e) { return 2 * std::atan2(std::sqrt(1 - e) * std::sin(nu / 2), std::sqrt(1 + e) * std::cos(nu / 2)); }
inline double true_from_hyperbolic(double H, double e) { return 2 * std::atan(std::sqrt((e + 1) / (e - 1)) * std::tanh(H / 2)); }
inline double hyperbolic_from_true(double nu, double e) { return 2 * std::atanh(std::sqrt((e - 1) / (e + 1)) * std::tan(nu / 2)); }

// The mean anomaly of a true anomaly: E − e sin E, e sinh H − H, or Barker's D + D³/3 (D = tan ν/2).
inline double mean_anomaly(double nu, double e) {
  if (e < 1 - 1e-12) { const double E = eccentric_from_true(nu, e); return E - e * std::sin(E); }
  if (e > 1 + 1e-12) { const double H = hyperbolic_from_true(nu, e); return e * std::sinh(H) - H; }
  const double D = std::tan(nu / 2);
  return D + D * D * D / 3;
}

// The true anomaly of a mean anomaly on any conic (Barker's equation in closed form for a parabola).
inline double true_anomaly(double M, double e) {
  if (e < 1 - 1e-12) return true_from_eccentric(eccentric_anomaly(M, e), e);
  if (e > 1 + 1e-12) return true_from_hyperbolic(hyperbolic_anomaly(M, e), e);
  const double w = std::cbrt(1.5 * M + std::sqrt(1 + 2.25 * M * M));
  return 2 * std::atan(w - 1 / w);
}

// ---- propagation by universal variables -----------------------------------------------------------------

// The Stumpff functions c₂(ψ) = (1 − cos √ψ)/ψ and c₃(ψ) = (√ψ − sin √ψ)/√ψ³ (and their hyperbolic
// continuations), by their series near ψ = 0.
inline void stumpff(double psi, double& c2, double& c3) {
  if (psi > 1e-6) {
    const double s = std::sqrt(psi);
    c2 = (1 - std::cos(s)) / psi;
    c3 = (s - std::sin(s)) / (psi * s);
  } else if (psi < -1e-6) {
    const double s = std::sqrt(-psi);
    c2 = (1 - std::cosh(s)) / psi;
    c3 = (std::sinh(s) - s) / (-psi * s);
  } else {
    c2 = 0.5 - psi / 24 + psi * psi / 720;
    c3 = 1.0 / 6 - psi / 120 + psi * psi / 5040;
  }
}

// The state (r, v) carried over dt seconds about μ on its conic (KEPLER).
inline std::array<Vec3, 2> propagate(const Vec3& r0, const Vec3& v0, double dt, double mu) {
  if (dt == 0) return {r0, v0};
  const double sm = std::sqrt(mu), rn = norm(r0), rv = dot(r0, v0), v2 = dot(v0, v0);
  const double alpha = 2 / rn - v2 / mu;  // 1/a
  double chi;
  if (alpha > 1e-12 / rn) chi = sm * dt * alpha;  // ellipse
  else if (alpha < -1e-12 / rn) {                  // hyperbola
    const double a = 1 / alpha;
    chi = std::copysign(1.0, dt) * std::sqrt(-a) * std::log((-2 * mu * alpha * dt) / (rv + std::copysign(1.0, dt) * std::sqrt(-mu * a) * (1 - rn * alpha)));
  } else {  // parabola
    const Vec3 h = cross(r0, v0);
    const double p = dot(h, h) / mu, s = 0.5 * std::atan(1 / (3 * std::sqrt(mu / (p * p * p)) * dt));
    const double w = std::atan(std::cbrt(std::tan(s)));
    chi = std::sqrt(p) * 2 / std::tan(2 * w);
  }
  double c2 = 0, c3 = 0, psi = 0, r = rn;
  for (int k = 0; k < 100; k++) {
    psi = chi * chi * alpha;
    stumpff(psi, c2, c3);
    r = chi * chi * c2 + rv / sm * chi * (1 - psi * c3) + rn * (1 - psi * c2);
    const double t = chi * chi * chi * c3 + rv / sm * chi * chi * c2 + rn * chi * (1 - psi * c3);
    const double d = (sm * dt - t) / r;
    chi += d;
    if (std::fabs(d) < 1e-13 * std::max(1.0, std::fabs(chi))) break;
  }
  psi = chi * chi * alpha;
  stumpff(psi, c2, c3);
  r = chi * chi * c2 + rv / sm * chi * (1 - psi * c3) + rn * (1 - psi * c2);
  const double f = 1 - chi * chi / rn * c2, g = dt - chi * chi * chi / sm * c3;
  const double gd = 1 - chi * chi / r * c2, fd = sm / (r * rn) * chi * (psi * c3 - 1);
  return {f * r0 + g * v0, fd * r0 + gd * v0};
}

// ---- frames of an orbit -------------------------------------------------------------------------------------

enum class OrbitFrame { rtn, lvlh, vnc, perifocal };

namespace detail {
// A unit vector and its rate from a vector and its rate.
inline std::array<Vec3, 2> unit_and_rate(const Vec3& w, const Vec3& wd) {
  const double n = norm(w);
  const Vec3 u = (1 / n) * w;
  return {u, (1 / n) * (wd - dot(u, wd) * u)};
}
inline Mat3 columns(const Vec3& x, const Vec3& y, const Vec3& z) {
  return Mat3{{{x[0], y[0], z[0]}, {x[1], y[1], z[1]}, {x[2], y[2], z[2]}}};
}
}  // namespace detail

// The frame of the orbit at the state (r, v) with acceleration a, as the motion from the frame to the
// inertial one: x_inertial = R x_frame + r, with its rate.
inline StateTransform orbit_frame(OrbitFrame kind, const Vec3& r, const Vec3& v, const Vec3& a, double mu = 0) {
  const Vec3 h = cross(r, v), hd = cross(r, a);
  const auto [N, Nd] = detail::unit_and_rate(h, hd);
  StateTransform t;
  t.b = r;
  t.b_dot = v;
  switch (kind) {
    case OrbitFrame::rtn: {
      const auto [R, Rd] = detail::unit_and_rate(r, v);
      const Vec3 T = cross(N, R), Td = cross(Nd, R) + cross(N, Rd);
      t.r = detail::columns(R, T, N);
      t.r_dot = detail::columns(Rd, Td, Nd);
      break;
    }
    case OrbitFrame::lvlh: {
      const auto [R, Rd] = detail::unit_and_rate(r, v);
      const Vec3 z = -R, zd = -Rd, y = -N, yd = -Nd;
      const Vec3 x = cross(y, z), xd = cross(yd, z) + cross(y, zd);
      t.r = detail::columns(x, y, z);
      t.r_dot = detail::columns(xd, yd, zd);
      break;
    }
    case OrbitFrame::vnc: {
      const auto [V, Vd] = detail::unit_and_rate(v, a);
      const Vec3 C = cross(V, N), Cd = cross(Vd, N) + cross(V, Nd);
      t.r = detail::columns(V, N, C);
      t.r_dot = detail::columns(Vd, Nd, Cd);
      break;
    }
    case OrbitFrame::perifocal: {
      if (!(mu > 0)) throw std::invalid_argument("physis: the perifocal frame needs μ");
      const double rn = norm(r);
      const Vec3 e = (1 / mu) * ((dot(v, v) - mu / rn) * r - dot(r, v) * v);
      if (norm(e) < 1e-11) throw std::invalid_argument("physis: a circular orbit has no perifocal frame");
      // e = (v × h)/μ − r̂, so ė = (a × h + v × ḣ)/μ − r̂′: zero for the central force alone.
      const Vec3 rhat = (1 / rn) * r, rhat_dot = (1 / rn) * (v - dot(rhat, v) * rhat);
      const Vec3 ed = (1 / mu) * (cross(a, h) + cross(v, hd)) - rhat_dot;
      const auto [P, Pd] = detail::unit_and_rate(e, ed);
      const Vec3 Q = cross(N, P), Qd = cross(Nd, P) + cross(N, Pd);
      t.r = detail::columns(P, Q, N);
      t.r_dot = detail::columns(Pd, Qd, Nd);
      break;
    }
  }
  return t;
}

// The same, under the central force alone (a = −μ r / r³).
inline StateTransform orbit_frame(OrbitFrame kind, const Vec3& r, const Vec3& v, double mu) {
  const double rn = norm(r);
  return orbit_frame(kind, r, v, (-mu / (rn * rn * rn)) * r, mu);
}

}  // namespace physis::mechanics
