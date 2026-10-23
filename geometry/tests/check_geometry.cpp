// Checks of physis/geometry: the Euclidean plane in polar coordinates is flat, with its
// Christoffel symbols in closed form; the round 2-sphere of radius a has scalar curvature
// 2/a² and R^θ_φθφ = sin²θ; the inverse, determinant and inertia of small forms; great
// circles as geodesics of the sphere, closing after one circumference and keeping g(u, u)
// and the angular momentum p_φ.

#include <physis/geometry/curvature.hpp>
#include <physis/geometry/geodesic.hpp>
#include <physis/geometry/metric.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>

namespace pg = physis::geometry;

int main() {
  calx::verification::Suite s;
  const double pi = std::numbers::pi;

  // ---- small forms -----------------------------------------------------------------------------------------------------
  {
    const pg::Tensor2<3> a{{{4, 1, -2}, {1, 3, 0.5}, {-2, 0.5, 5}}};
    const auto ai = pg::inverse(a);
    double dev = 0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        double p = 0;
        for (int k = 0; k < 3; k++) p += a[i][k] * ai[k][j];
        dev = std::max(dev, std::fabs(p - (i == j)));
      }
    s.check("geometry: inverse of a symmetric 3x3 form, |A A^-1 - I|", dev, 1e-15);
    s.check("geometry: determinant against the rule of Sarrus", pg::determinant(a) - (4 * (15 - 0.25) - 1 * (5 + 1) - 2 * (0.5 + 6)), 1e-13);
    const pg::Tensor2<4> eta{{{-1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}}};
    const pg::Tensor2<2> hyperbolic{{{0, 1}, {1, 0}}};  // the form 2xy, of signature (1, 1)
    s.require("geometry: inertia: Minkowski has one negative direction, 2xy one, a positive form none",
              pg::negative_directions(eta) == 1 && pg::negative_directions(hyperbolic) == 1 && pg::negative_directions(a) == 0);
  }
  // ---- the plane in polar coordinates -----------------------------------------------------------------------------------
  {
    const pg::PolarPlane plane;
    const pg::Point<2> x{1.7, 0.4};
    const auto G = pg::christoffel(plane, x);
    s.check("geometry: polar plane, Γ^r_φφ = -r", G[0][1][1] + 1.7, 1e-15);
    s.check("geometry: polar plane, Γ^φ_rφ = 1/r", G[1][0][1] - 1 / 1.7, 1e-15);
    const auto k = pg::curvature(plane, x);
    double worst = 0;
    for (const auto& a : k.Riemann)
      for (const auto& b : a)
        for (const auto& c : b)
          for (double v : c) worst = std::max(worst, std::fabs(v));
    s.check("geometry: the polar plane is flat, max |R^a_bcd|", worst, 1e-15);
  }
  // ---- the round sphere ------------------------------------------------------------------------------------------------
  {
    const pg::RoundSphere sphere{2.5};
    const double theta = 0.9;
    const auto k = pg::curvature(sphere, pg::Point<2>{theta, 1.1});
    s.check("geometry: the sphere of radius a has scalar curvature 2/a^2 (relative)", k.scalar / (2 / (2.5 * 2.5)) - 1, 1e-14);
    s.check("geometry: and R^θ_φθφ = sin^2 θ", k.Riemann[0][1][0][1] - std::sin(theta) * std::sin(theta), 1e-14);
    s.check("geometry: its Einstein tensor vanishes in two dimensions", std::fabs(k.Einstein[0][0]) + std::fabs(k.Einstein[1][1]) + std::fabs(k.Einstein[0][1]), 1e-14);
    s.check("geometry: the Kretschmann scalar of the sphere is R^2 = 4/a^4 (relative)", pg::kretschmann(k) / (4 / std::pow(2.5, 4)) - 1, 1e-14);
  }
  // ---- great circles ---------------------------------------------------------------------------------------------------
  {
    // A unit-speed geodesic of the sphere of radius a, from the equator heading north-east, is a
    // great circle: it returns after the arc length 2π a, keeps g(u, u) = 1 and p_φ.
    const double a = 1.3;
    const pg::RoundSphere sphere{a};
    const double heading = 0.6;  // angle from east
    pg::Phase<2> y = pg::phase<2>({pi / 2, 0.2}, {-std::sin(heading) / a, std::cos(heading) / a});
    const double p0 = pg::momentum(sphere, y, 1);
    double drift = 0, pdrift = 0;
    const std::size_t steps = 4000;
    const double h = 2 * pi * a / steps;
    pg::integrate_geodesic(sphere, y, h, steps, [&](double, const pg::Phase<2>& z) {
      drift = std::max(drift, std::fabs(pg::norm(sphere, z) - 1));
      pdrift = std::max(pdrift, std::fabs(pg::momentum(sphere, z, 1) - p0));
      return true;
    });
    const double phi = std::remainder(y[1] - 0.2, 2 * pi);
    s.check("geometry: a great circle closes after one circumference, |Δθ| + |Δφ| (RK4, 4000 steps)", std::fabs(y[0] - pi / 2) + std::fabs(phi), 1e-11);
    s.check("geometry: g(u, u) along it", drift, 1e-12);
    s.check("geometry: p_φ along it", pdrift, 1e-12);
  }
  return s.done();
}
