// Checks of the two-body problem (physis/mechanics): elements and states, Kepler's equation,
// propagation by universal variables, and the frames of an orbit.
//
// References: Vallado's RV2COE (the code of Fundamentals of Astrodynamics and Applications, through
// python-sgp4's copy; tests/fixtures/rv2coe-reference.txt, made by tests/make_fixtures.py) for the
// elements of a state, his Example 2-5 among them; the two-body equations of motion integrated by
// calx's Dormand–Prince pair at tight tolerance, as an independent method, for the propagation; the
// invariants (energy, angular momentum, the eccentricity vector); finite differences for the
// frames' rates.

#include <cmath>
#include <fstream>
#include <numbers>
#include <sstream>
#include <string>

#include <calx/integrators/runge_kutta.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <physis/mechanics/two_body.hpp>

using namespace physis::mechanics;
using calx::verification::Suite;

namespace {
constexpr double pi = std::numbers::pi;
constexpr double MU = 3.986004418e14;  // the Earth's GM (m³/s², IERS Conventions 2010, Table 1.1)
double max_abs(const Vec3& v) { return std::max({std::fabs(v[0]), std::fabs(v[1]), std::fabs(v[2])}); }
double angle_gap(double a, double b) { return std::fabs(std::remainder(a - b, 2 * pi)); }
}  // namespace

int main() {
  Suite s;

  // ---- against Vallado's RV2COE ----------------------------------------------------------------------------
  // His code's elements of his Example 2-5's state and of 119 random eccentric, inclined orbits
  // (tests/fixtures/rv2coe-reference.txt, through python-sgp4's copy, sgp4.ext.rv2coe).
  {
    std::ifstream in(PHYSIS_MECHANICS_FIXTURES "/rv2coe-reference.txt");
    std::string line;
    const double km = 1000, mu = 398600.4418e9;
    double worst_shape = 0, worst_angle = 0;
    int n = 0;
    Elements example{};
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream f(line);
      double x[13];
      for (double& v : x) f >> v;
      const Elements o = elements({x[0] * km, x[1] * km, x[2] * km}, {x[3] * km, x[4] * km, x[5] * km}, mu);
      if (n == 0) example = o;
      worst_shape = std::max({worst_shape, std::fabs(o.p / km - x[6]) / x[6], std::fabs(o.a() / km - x[7]) / std::fabs(x[7]), std::fabs(o.e - x[8])});
      worst_angle = std::max({worst_angle, std::fabs(o.i - x[9]), angle_gap(o.raan, x[10]), angle_gap(o.argp, x[11]), angle_gap(o.nu, x[12])});
      n++;
    }
    s.require("the reference file was read", n >= 100);
    s.check("elements against Vallado's RV2COE: p, a (relative) and e", worst_shape, 1e-12);
    s.check("elements against Vallado's RV2COE: i, Ω, ω, ν", worst_angle, 1e-12, "rad");
    // Example 2-5 as his code gives it: p = 11 067.798 km, e = 0.832853, i = 87.86913°.
    s.check("Vallado's Example 2-5: p = 11 067.798 km", example.p / km - 11067.798, 5e-4, "km");
    s.check("Vallado's Example 2-5: i = 87.86913°", example.i * 180 / pi - 87.86913, 5e-6, "deg");
  }

  // ---- elements and states, round trips over every kind of orbit -----------------------------------------
  {
    calx::random::Stream rng(20261028, 1, 0, 0);
    double worst_r = 0, worst_v = 0, worst_angles = 0;
    int kinds[4] = {0, 0, 0, 0};
    for (int k = 0; k < 400; k++) {
      Elements o;
      const int kind = k % 6;
      o.e = kind == 1 || kind == 3 ? 0 : kind == 5 ? 1.2 + 3 * rng.uniform() : 0.95 * rng.uniform();
      o.i = kind == 2 || kind == 3 ? 0 : kind == 4 ? pi : pi * rng.uniform();
      o.p = (6.6e6 + 4e7 * rng.uniform()) * (1 - o.e * o.e > 0 ? (1 - o.e * o.e) : 1);
      o.raan = 2 * pi * rng.uniform();
      o.argp = 2 * pi * rng.uniform();
      o.nu = o.e > 1 ? (rng.uniform() - 0.5) * 2 * std::acos(-1 / o.e) * 0.9 : 2 * pi * rng.uniform();
      const auto x = state(o, MU);
      const Elements b = elements(x[0], x[1], MU);
      kinds[int(b.kind)]++;
      const auto y = state(b, MU);
      worst_r = std::max(worst_r, max_abs(y[0] - x[0]) / std::max(1.0, max_abs(x[0])));
      worst_v = std::max(worst_v, max_abs(y[1] - x[1]) / std::max(1.0, max_abs(x[1])));
      if (b.kind == OrbitKind::elliptical_inclined && o.e > 1e-6 && o.i > 1e-6 && o.i < pi - 1e-6)
        worst_angles = std::max({worst_angles, angle_gap(b.raan, o.raan), angle_gap(b.argp, o.argp), angle_gap(b.nu, o.nu), std::fabs(b.i - o.i)});
    }
    s.check("elements → state → elements → state, 400 orbits (relative)", std::max(worst_r, worst_v), 1e-12);
    s.check("the angles of inclined, eccentric orbits are recovered", worst_angles, 1e-10, "rad");
    s.require("every kind of orbit was met: elliptical and circular, inclined and equatorial",
              kinds[0] > 0 && kinds[1] > 0 && kinds[2] > 0 && kinds[3] > 0);
  }

  // ---- Kepler's equation -----------------------------------------------------------------------------------
  {
    double worst = 0, worst_h = 0, worst_nu = 0;
    for (double e : {0.0, 0.1, 0.5, 0.9, 0.99, 0.999}) {
      for (int k = 0; k < 64; k++) {
        const double M = -pi + 2 * pi * k / 63.0;
        const double E = eccentric_anomaly(M, e);
        worst = std::max(worst, std::fabs(E - e * std::sin(E) - M));
        worst_nu = std::max(worst_nu, angle_gap(mean_anomaly(true_anomaly(M, e), e), M));
      }
    }
    for (double e : {1.01, 1.5, 3.0, 20.0})
      for (double M : {-50.0, -3.0, -0.1, 0.0, 0.2, 4.0, 100.0}) {
        const double H = hyperbolic_anomaly(M, e);
        worst_h = std::max(worst_h, std::fabs(e * std::sinh(H) - H - M) / std::max(1.0, std::fabs(M)));
      }
    s.check("Kepler's equation E − e sin E = M, e from 0 to 0.999", worst, 4e-15, "rad");
    s.check("the hyperbolic Kepler equation e sinh H − H = M (relative)", worst_h, 1e-14);
    s.check("mean → true → mean anomaly on ellipses", worst_nu, 1e-11, "rad");
    const double D = 0.7, M = D + D * D * D / 3;
    s.check("Barker's equation for a parabola, in closed form", true_anomaly(M, 1.0) - 2 * std::atan(D), 1e-14, "rad");
  }

  // ---- propagation against integration ----------------------------------------------------------------------
  {
    // Three conics, each carried over a day by universal variables and by Dormand–Prince at a tolerance
    // of 10⁻¹³: they agree to the integrator's error.
    struct Case { const char* name; Vec3 r, v; double dt; };
    const Case cases[] = {
        {"a low, nearly circular orbit", {6.778e6, 0, 0}, {0, 7.66e3, 0.3e3}, 86400},
        {"a Molniya orbit (e ≈ 0.74)", {6.9e6, 0, 0}, {0, 9.0e3, 5.7e3}, 86400},
        {"a hyperbolic flyby (v∞ ≈ 4.5 km/s)", {7.0e6, 0, 0}, {0, 0, 12.0e3}, 20000},
    };
    for (const Case& c : cases) {
      const auto u = propagate(c.r, c.v, c.dt, MU);
      calx::integrators::Vec<6> x{c.r[0], c.r[1], c.r[2], c.v[0], c.v[1], c.v[2]};
      auto f = [](const calx::integrators::Vec<6>& y, calx::integrators::Vec<6>& d) {
        const double r = std::sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]), k = -MU / (r * r * r);
        d = {y[3], y[4], y[5], k * y[0], k * y[1], k * y[2]};
      };
      double sub = 10;
      const std::array<double, 6> scale{1e3, 1e3, 1e3, 1, 1, 1};
      for (int step = 0; step < 200; step++) calx::integrators::dopri5<6>(f, x, c.dt / 200, scale, {1e-13, 1e-13}, sub);
      const double dr = max_abs(u[0] - Vec3{x[0], x[1], x[2]});
      s.check(std::string("propagation by universal variables against integration: ") + c.name, dr / norm(c.r), 2e-9);
      // The invariants: energy and angular momentum.
      const double e0 = dot(c.v, c.v) / 2 - MU / norm(c.r), e1 = dot(u[1], u[1]) / 2 - MU / norm(u[0]);
      s.check(std::string("energy kept: ") + c.name, std::fabs(e1 - e0) / std::fabs(e0), 1e-12);
      s.check(std::string("angular momentum kept: ") + c.name, max_abs(cross(u[0], u[1]) - cross(c.r, c.v)) / norm(cross(c.r, c.v)), 1e-12);
    }
    // There and back.
    const Vec3 r0{7.0e6, 1.0e6, -2.0e6}, v0{-1.0e3, 7.0e3, 1.5e3};
    const auto there = propagate(r0, v0, 12345.6, MU), back = propagate(there[0], there[1], -12345.6, MU);
    s.check("propagation there and back", max_abs(back[0] - r0), 1e-5, "m");
    // Kepler's law: after one period the state returns.
    const double a = 1 / (2 / norm(r0) - dot(v0, v0) / MU), T = 2 * pi * std::sqrt(a * a * a / MU);
    s.check("after one period the body is where it was", max_abs(propagate(r0, v0, T, MU)[0] - r0) / norm(r0), 1e-11);
  }

  // ---- frames of an orbit -----------------------------------------------------------------------------------
  {
    const Vec3 r{7.0e6, 1.0e6, -2.0e6}, v{-1.0e3, 7.0e3, 1.5e3};
    double ortho = 0, rate = 0;
    for (OrbitFrame k : {OrbitFrame::rtn, OrbitFrame::lvlh, OrbitFrame::vnc, OrbitFrame::perifocal}) {
      const StateTransform t = orbit_frame(k, r, v, MU);
      ortho = std::max({ortho, physis::geometry::orthogonality_defect(t.r), std::fabs(physis::geometry::determinant(t.r) - 1)});
      // The rate against a central difference along the orbit.
      const double h = 0.05;  // s: truncation h² ω³ / 6 ≈ 4 × 10⁻¹³ at ω ≈ 10⁻³ rad/s
      const auto a = propagate(r, v, -h, MU), b = propagate(r, v, h, MU);
      const Mat3 d = (1 / (2 * h)) * (orbit_frame(k, b[0], b[1], MU).r - orbit_frame(k, a[0], a[1], MU).r);
      rate = std::max(rate, physis::geometry::max_difference(d, t.r_dot));
    }
    s.check("the frames of an orbit are rotations", ortho, 1e-14);
    s.check("the frames' rates against central differences along the orbit", rate, 1e-11, "1/s");
    const StateTransform rtn = orbit_frame(OrbitFrame::rtn, r, v, MU), lvlh = orbit_frame(OrbitFrame::lvlh, r, v, MU);
    s.check("RTN's R is the radial direction", max_abs(rtn.direction({1, 0, 0}) - (1 / norm(r)) * r), 1e-15);
    s.check("LVLH's z points to the centre", max_abs(lvlh.direction({0, 0, 1}) + (1 / norm(r)) * r), 1e-15);
    const Vec3 rc{7.0e6, 0, 0}, vc{0, std::sqrt(MU / 7.0e6), 0};
    s.check("on a circle, LVLH's x is the velocity", max_abs(orbit_frame(OrbitFrame::lvlh, rc, vc, MU).direction({1, 0, 0}) - Vec3{0, 1, 0}), 1e-15);
    const Mat3 pqw = orbit_frame(OrbitFrame::perifocal, r, v, MU).r_dot;
    s.check("the perifocal frame stands still under the central force (to rounding)", physis::geometry::max_difference(pqw, physis::geometry::zero3()), 1e-18, "1/s");
    // The angular velocity of RTN under the central force is h / r² along the normal.
    const Vec3 w = orbit_frame(OrbitFrame::rtn, r, v, MU).angular_velocity(), h = cross(r, v);
    s.check("RTN turns at h / r² about the orbit's normal", max_abs(w - (1 / dot(r, r)) * h), 1e-18, "rad/s");
  }
  return s.done();
}
