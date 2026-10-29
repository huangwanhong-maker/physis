// Checks of physis fields: sections over domains, their representations, capabilities, compositions,
// samplings, integral curves, versions and reconstructions, against closed forms.
//
// - A Gauss dipole V = a³ g₁⁰ z / r³ (B = −∇V = a³ g₁⁰ (3 z x / r⁵ − ẑ / r³)) as a formula, whose
//   derivative by dual numbers has zero divergence and curl, and as a spectral expansion of degree 1,
//   which must give the formula's field; its field lines lie on r = L cos² ψ (Chapman and Bartels
//   1940, §16.2); its field at the pole is twice that at the equator.
// - A rigid rotation v = ω × x, whose curl is 2ω and whose streamlines are circles of period 2π/|ω|.
// - Affine fields, which trilinear reconstruction reproduces exactly.

#include <array>
#include <cmath>
#include <numbers>
#include <string>
#include <vector>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <physis/fields/analytic.hpp>
#include <physis/fields/any.hpp>
#include <physis/fields/compose.hpp>
#include <physis/fields/curves.hpp>
#include <physis/fields/domains.hpp>
#include <physis/fields/history.hpp>
#include <physis/fields/potential.hpp>
#include <physis/fields/reductions.hpp>
#include <physis/fields/sampled.hpp>
#include <physis/fields/sampling.hpp>

using calx::verification::Suite;
namespace pf = physis::fields;
namespace pa = physis::analysis;
using pf::Mat3;
using pf::Vec3;
using namespace physis::geometry::operators;

namespace {

constexpr double A = 1737400.0, G10 = -3.0e-9;  // a reference sphere and a dipole coefficient (T)
constexpr double A3G = A * A * A * G10;

double rel(const Vec3& a, const Vec3& b) { return physis::geometry::norm(a - b) / std::max(physis::geometry::norm(b), 1e-300); }

pf::Descriptor described_as(std::string name, pf::FibreKind kind, physis::units::Dimension dim) {
  pf::Descriptor d;
  d.name = name;
  d.title = name;
  d.frame = "body";
  d.time = "none";
  d.fibre.kind = kind;
  d.fibre.dimension = dim;
  d.fibre.basis = "frame";
  return d;
}

// The dipole as a formula over any number type.
auto dipole_formula() {
  return [](double, const auto& x) {
    using std::sqrt;
    using S = std::decay_t<decltype(x[0])>;
    const S r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
    const S r = sqrt(r2);
    const S r3 = r2 * r, r5 = r3 * r2;
    return std::array<S, 3>{A3G * 3 * x[2] * x[0] / r5, A3G * 3 * x[2] * x[1] / r5, A3G * (3 * x[2] * x[2] / r5 - 1 / r3)};
  };
}

pa::HarmonicPotential dipole_expansion() {
  pa::HarmonicPotential p = pa::HarmonicPotential::zero(1, A);
  p.g[pa::HarmonicPotential::index(1, 0)] = G10;
  return p;
}

pa::HarmonicPotential random_expansion(int degree, std::uint32_t key) {
  calx::random::Stream rng(7, key, 0, 0);
  pa::HarmonicPotential p = pa::HarmonicPotential::zero(degree, A);
  for (int n = 1; n <= degree; n++)
    for (int m = 0; m <= n; m++) {
      const std::size_t k = pa::HarmonicPotential::index(n, m);
      p.g[k] = 1e-9 * rng.normal() / (n + 1);
      p.h[k] = m ? 1e-9 * rng.normal() / (n + 1) : 0.0;
    }
  return p;
}

Vec3 random_point(calx::random::Stream& rng, double rmin, double rmax) {
  const double z = rng.uniform(-1, 1), l = rng.uniform(-std::numbers::pi, std::numbers::pi), r = rng.uniform(rmin, rmax);
  const double s = std::sqrt(1 - z * z);
  return {r * s * std::cos(l), r * s * std::sin(l), r * z};
}

}  // namespace

int main() {
  Suite s;
  const auto T = physis::units::dim::magnetic_field;
  calx::random::Stream rng(20261029, 1, 0, 0);

  // ---- domains --------------------------------------------------------------------------------------------------
  {
    const pf::Shell shell{{}, A, 2 * A};
    const auto during = pf::During{0, 10};
    const auto both = pf::all_of(shell, during);
    const auto in = both.covers(pf::Query<Vec3>{5, {1.5 * A, 0, 0}}), low = both.covers(pf::Query<Vec3>{5, {0.5 * A, 0, 0}}),
               late = both.covers(pf::Query<Vec3>{11, {1.5 * A, 0, 0}});
    s.require("domains: a shell in a span of time covers its inside; below and after, outside with their reasons",
              in.ok() && low.quality == pf::Quality::outside && low.reason == shell.below && late.quality == pf::Quality::outside && late.reason == during.after);
  }

  // ---- analytic fields and their exact derivatives ---------------------------------------------------------------------
  const auto dipole = pf::analytic<Vec3>(described_as("dipole", pf::FibreKind::vector, T), dipole_formula(), pf::Shell{{}, 0.5 * A});
  {
    auto phi = pf::analytic<double>(described_as("phi", pf::FibreKind::scalar, {}), [](double, const auto& x) {
      using std::sin;
      return x[0] * x[0] * x[1] + sin(x[2]);
    });
    double worst = 0;
    for (int i = 0; i < 20; i++) {
      const Vec3 x{rng.uniform(-2, 2), rng.uniform(-2, 2), rng.uniform(-2, 2)};
      const auto g = phi.jacobian({0, x});
      worst = std::max(worst, rel(g.value, {2 * x[0] * x[1], x[0] * x[0], std::cos(x[2])}));
    }
    s.check("analytic: the gradient of x² y + sin z by dual numbers, against its closed form (20 points)", worst, 1e-15);
    double div = 0, cur = 0, fd = 0;
    for (int i = 0; i < 50; i++) {
      const Vec3 x = random_point(rng, 1.1 * A, 3 * A);
      const auto J = dipole.jacobian({0, x});
      const double scale = physis::geometry::norm(dipole.at({0, x}).value) / physis::geometry::norm(x);
      div = std::max(div, std::fabs(pf::divergence(J.value)) / scale);
      cur = std::max(cur, physis::geometry::norm(pf::curl(J.value)) / scale);
      const auto D = pf::central_differences(dipole, {0, x}, 1e-4 * physis::geometry::norm(x));
      double d = 0;
      for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) d = std::max(d, std::fabs(D.value[a][b] - J.value[a][b]) / scale);
      fd = std::max(fd, d);
    }
    s.check("analytic: the dipole's divergence by dual numbers, of |B|/r (50 points)", div, 1e-14);
    s.check("analytic: and its curl", cur, 1e-14);
    // Their truncation error is h²/6 times the third derivative, some 60 |B|/r³ for a dipole: 10⁻⁷ |B|/r at h = 10⁻⁴ r.
    s.check("analytic: central differences of step 10⁻⁴ r agree with the dual numbers, of |B|/r", fd, 2e-7);
    s.require("analytic: a formula of dual numbers is differentiable, a uniform field too, and both are pointwise",
              pf::capabilities_of<decltype(dipole)>().has(pf::Capability::jacobian) && pf::capabilities_of<decltype(dipole)>().has(pf::Capability::point));
    const pf::Uniform<Vec3> u(described_as("uniform", pf::FibreKind::vector, T), Vec3{1e-9, 2e-9, 3e-9});
    const auto uj = u.jacobian({0, {1, 2, 3}});
    s.require("uniform: one value everywhere, derivative zero", u.at({0, {5, 6, 7}}).value == Vec3{1e-9, 2e-9, 3e-9} && uj.value == Mat3{});
  }

  // ---- potential fields: the spectral representation -------------------------------------------------------------------
  const pf::PotentialField spectral_dipole(described_as("dipole expansion", pf::FibreKind::vector, T), dipole_expansion());
  const pf::PotentialField rough(described_as("degree 20", pf::FibreKind::vector, T), random_expansion(20, 2));
  {
    double worst = 0;
    for (int i = 0; i < 50; i++) {
      const Vec3 x = random_point(rng, A, 4 * A);
      worst = std::max(worst, rel(spectral_dipole.at({0, x}).value, dipole.at({0, x}).value));
    }
    s.check("potential: an expansion of degree 1 gives the dipole's formula (50 points)", worst, 1e-14);
    const pf::SphereGrid grid{{}, 1.01 * A, 18, 36};
    const auto fast = pf::sample(rough, grid, 0);
    double gap = 0;
    for (std::size_t k = 0; k < grid.size(); k++) gap = std::max(gap, rel(fast.values[k], rough.at({0, grid.point(k)}).value));
    s.check("potential: a sphere read row by row equals the expansion point by point (degree 20, 648 cells)", gap, 1e-13);
    const auto below = rough.at({0, {0.9 * A, 0, 0}});
    s.require("potential: inside the reference sphere, outside, with the reason",
              below.quality == pf::Quality::outside && below.reason == "below the expansion's reference sphere, inside its sources");
    double div = 0, cur = 0;
    for (int i = 0; i < 20; i++) {
      const Vec3 x = random_point(rng, 1.2 * A, 2 * A);
      const auto D = pf::central_differences(rough, {0, x}, 2e-4 * A);
      const double scale = physis::geometry::norm(rough.at({0, x}).value) / A;
      div = std::max(div, std::fabs(pf::divergence(D.value)) / scale);
      cur = std::max(cur, physis::geometry::norm(pf::curl(D.value)) / scale);
    }
    s.check("potential: a field of degree 20 has no divergence (by differences of step 2·10⁻⁴ a, of |B|/a)", div, 1e-5);
    s.check("potential: nor curl", cur, 1e-5);
    double minus_grad = 0;
    for (int i = 0; i < 10; i++) {
      const Vec3 x = random_point(rng, 1.2 * A, 2 * A);
      Vec3 g{};
      for (int a = 0; a < 3; a++) {
        Vec3 p = x, m = x;
        p[a] += 1e-4 * A;
        m[a] -= 1e-4 * A;
        g[a] = -(rough.potential({0, p}).value - rough.potential({0, m}).value) / (2e-4 * A);
      }
      minus_grad = std::max(minus_grad, rel(g, rough.at({0, x}).value));
    }
    s.check("potential: B = −∇V, the potential differenced against the field", minus_grad, 1e-6);
    double local = 0;
    for (int i = 0; i < 20; i++) {
      const Vec3 x = random_point(rng, A, 3 * A);
      const double r = physis::geometry::norm(x), psi = std::atan2(x[2], std::hypot(x[0], x[1])), lambda = std::atan2(x[1], x[0]);
      const auto f = pa::field_at(rough.expansion(), r, psi, lambda);
      const Vec3 enu = pf::LocalOfSphere{}({0, x}, rough.at({0, x}).value);
      local = std::max(local, rel(enu, {f.east, f.north, -f.down}));
    }
    s.check("potential: on the local east, north and up the field is analysis's (east, north, −down)", local, 1e-14);
  }

  // ---- reductions and compositions --------------------------------------------------------------------------------------
  {
    const double pole = pf::Magnitude{}({}, spectral_dipole.at({0, {0, 0, 2 * A}}).value), equator = pf::Magnitude{}({}, spectral_dipole.at({0, {2 * A, 0, 0}}).value);
    s.check("reductions: the dipole's field at the pole is twice that at the equator", pole / equator - 2, 1e-14);
    s.check("reductions: at the equator |B| = |g₁⁰| (a/r)³", equator / (std::fabs(G10) / 8) - 1, 1e-14);
    const double radial = pf::Radial{}({0, {0, 0, 2 * A}}, spectral_dipole.at({0, {0, 0, 2 * A}}).value);
    const double horizontal = pf::Horizontal{}({0, {2 * A, 0, 0}}, spectral_dipole.at({0, {2 * A, 0, 0}}).value);
    s.check("reductions: radial at the pole, B_r = 2 g₁⁰ (a/r)³ (inward for g₁⁰ < 0, as the Earth's)", radial / (2 * G10 / 8) - 1, 1e-14);
    s.check("reductions: horizontal at the equator, |g₁⁰| (a/r)³", horizontal / (std::fabs(G10) / 8) - 1, 1e-14);

    const auto magnitude = pf::reduce(rough, pf::Magnitude{});
    const pf::SphereGrid grid{{}, 1.5 * A, 12, 24};
    const auto fast = pf::sample(magnitude, grid, 0);
    double gap = 0;
    for (std::size_t k = 0; k < grid.size(); k++) gap = std::max(gap, std::fabs(fast.values[k] - magnitude.at({0, grid.point(k)}).value) / fast.values[k]);
    s.check("reduce: |B| on a sphere by the expansion's rows equals |B| point by point", gap, 1e-13);
    s.require("reduce: |B| keeps the expansion's sphere grid; its fibre is a scalar of B's dimension, named",
              pf::described(magnitude).capabilities.has(pf::Capability::sphere_grid) && magnitude.descriptor().fibre.kind == pf::FibreKind::scalar &&
                  magnitude.descriptor().fibre.dimension == T && magnitude.descriptor().name == "degree 20.magnitude");

    auto crust = described_as("crust", pf::FibreKind::vector, T);
    crust.provenance = pf::Provenance::assimilated;
    auto finer = described_as("finer", pf::FibreKind::vector, T);
    finer.provenance = pf::Provenance::synthetic;
    const auto total = pf::sum(described_as("total", pf::FibreKind::vector, T), pf::PotentialField(crust, random_expansion(20, 3)),
                               pf::PotentialField(finer, random_expansion(30, 4)));
    const auto d = total.descriptor();
    s.require("sum: a composite of the most derived kind among its parts, its parts named with theirs",
              d.provenance == pf::Provenance::synthetic && d.parts.size() == 2 && d.parts[0].provenance == pf::Provenance::assimilated &&
                  d.parts[1].name == "finer" && d.representation.kind == pf::RepresentationKind::composite);
    const pf::SphereGrid g2{{}, 1.2 * A, 6, 12};
    const auto st = pf::sample(total, g2, 0);
    double add = 0;
    for (std::size_t k = 0; k < g2.size(); k++) add = std::max(add, rel(st.values[k], total.at({0, g2.point(k)}).value));
    s.check("sum: sampled part by part, equal to the sum point by point", add, 1e-13);
    s.require("sum: it reads a sphere row by row when both its parts do, and claims no derivative they lack",
              pf::described(total).capabilities.has(pf::Capability::sphere_grid) && !pf::described(total).capabilities.has(pf::Capability::jacobian) &&
                  !pf::described(pf::sum(described_as("mixed", pf::FibreKind::vector, T), rough, pf::Uniform<Vec3>(described_as("u", pf::FibreKind::vector, T), Vec3{})))
                       .capabilities.has(pf::Capability::sphere_grid));
    bool threw = false;
    try {
      (void)pf::sum(described_as("bad", pf::FibreKind::vector, T), rough,
                    pf::Uniform<Vec3>(described_as("velocity", pf::FibreKind::vector, physis::units::dim::velocity), Vec3{}));
    } catch (const std::invalid_argument&) { threw = true; }
    s.require("sum: parts of different dimensions are refused", threw);
    s.require("sum: a point inside a part's sources has no value, with that part's reason",
              total.at({0, {0.5 * A, 0, 0}}).quality == pf::Quality::outside);

    // A profile on height, pulled back to space by h(x) = |x| − R.
    const auto profile = pf::analytic<double, double>(described_as("temperature", pf::FibreKind::scalar, physis::units::dim::temperature),
                                                      [](double, double h) { return 288.15 - 0.0065 * h; }, pf::Interval{-5000, 86000});
    const double R = 6371000;
    const auto air = pf::pullback<Vec3>(described_as("air", pf::FibreKind::scalar, {}), profile, [R](const pf::Query<Vec3>& q) {
      const double h = physis::geometry::norm(q.x) - R;
      if (h > 90000) return pf::refused<double>(pf::Quality::missing, "no height above the geoid here");
      return pf::valid(h);
    });
    const auto at1 = air.at({0, {R + 1000, 0, 0}}), high = air.at({0, {R + 87000, 0, 0}}), higher = air.at({0, {R + 95000, 0, 0}});
    s.require("pullback: a profile on height read in space: its value, its fibre; above its interval outside; where the map fails, the map's reason",
              std::fabs(at1.value - (288.15 - 6.5)) < 1e-12 && air.descriptor().fibre.dimension == physis::units::dim::temperature &&
                  high.quality == pf::Quality::outside && higher.quality == pf::Quality::missing);

    // A rigid rotation described in a turned and shifted frame: values turned, curl turned.
    const Vec3 omega{0.3, -0.2, 0.7};
    const auto spin = pf::analytic<Vec3>(described_as("spin", pf::FibreKind::vector, physis::units::dim::velocity), [omega](double, const auto& x) {
      using S = std::decay_t<decltype(x[0])>;
      return std::array<S, 3>{omega[1] * x[2] - omega[2] * x[1], omega[2] * x[0] - omega[0] * x[2], omega[0] * x[1] - omega[1] * x[0]};
    });
    const Mat3 R30 = physis::geometry::rotation_matrix(physis::geometry::Axis::x, 0.5) * physis::geometry::rotation_matrix(physis::geometry::Axis::z, 0.3);
    const Vec3 shift{10, -4, 2};
    const auto turned = pf::in_frame(spin, [&](double) { return pf::Motion{R30, shift}; }, "turned");
    double val = 0, crl = 0, dv = 0;
    for (int i = 0; i < 20; i++) {
      const Vec3 x{rng.uniform(-5, 5), rng.uniform(-5, 5), rng.uniform(-5, 5)};
      const Vec3 back = physis::geometry::transpose(R30) * (x - shift);
      val = std::max(val, rel(turned.at({0, x}).value, R30 * spin.at({0, back}).value));
      const auto J = turned.jacobian({0, x}).value;
      crl = std::max(crl, rel(pf::curl(J), 2.0 * (R30 * omega)));
      dv = std::max(dv, std::fabs(pf::divergence(J)));
    }
    s.check("in_frame: v′(x′) = R v(Rᵀ(x′ − b)) (20 points)", val, 1e-14);
    s.check("in_frame: the curl of a rigid rotation in the turned frame is 2Rω, by dual numbers turned J′ = R J Rᵀ", crl, 1e-14);
    s.check("in_frame: its divergence stays zero", dv, 1e-15);
    const auto boxed = pf::restrict(spin, pf::Box{{-1, -1, -1}, {1, 1, 1}}, "the unit box");
    s.require("restrict: outside the box, outside; inside, the field's value and derivative",
              boxed.at({0, {2, 0, 0}}).quality == pf::Quality::outside && boxed.jacobian({0, {0.5, 0, 0}}).ok() && boxed.descriptor().base.about == "the unit box");
  }

  // ---- samplings ----------------------------------------------------------------------------------------------------------
  {
    const pf::SphereGrid g{{}, 2.0, 4, 8};
    const Vec3 first = g.point(0), last = g.point(g.size() - 1);
    s.require("sphere grid: the first cell's centre at 67.5° N, 157.5° W, the last at 67.5° S, 157.5° E, on the sphere",
              std::fabs(g.latitude(0) - 67.5 * std::numbers::pi / 180) < 1e-15 && std::fabs(g.longitude(0) + 157.5 * std::numbers::pi / 180) < 1e-15 &&
                  std::fabs(physis::geometry::norm(first) - 2) < 1e-15 && last[2] < 0 && std::fabs(std::atan2(last[1], last[0]) - 157.5 * std::numbers::pi / 180) < 1e-14);
    pf::VolumeGrid v{{1, 2, 3}, physis::geometry::rotation_matrix(physis::geometry::Axis::z, 0.4), {10, 20, 30}, {5, 3, 7}};
    const Vec3 mid = v.offset(2, 1, 3);
    const Vec3 x{4, -1, 9};
    s.require("volume grid: with odd n the middle cell's offset is exactly zero; the first axis fastest; points back and forth",
              mid == Vec3{} && v.index(1, 0, 0) == 1 && v.index(0, 1, 0) == 5 && v.index(0, 0, 1) == 15 && v.point(v.index(2, 1, 3)) == Vec3{1, 2, 3} &&
                  rel(v.from_axes(v.to_axes(x)), x) < 1e-15);
    const pf::PlaneGrid p{{0, 0, 0}, {2, 0, 0}, {0, 4, 0}, 2, 2};
    s.require("plane grid: cells' centres at the quarters", p.point(0) == Vec3{0.5, 1, 0} && p.point(3) == Vec3{1.5, 3, 0});
  }

  // ---- integral curves ----------------------------------------------------------------------------------------------------
  {
    pf::CurveOptions o;
    o.step_fraction = 0.01;
    const double psi0 = 30 * std::numbers::pi / 180, r0 = A + 100e3;
    const Vec3 seed{r0 * std::cos(psi0), 0, r0 * std::sin(psi0)};
    const auto line = pf::integral_curve(spectral_dipole, seed, 0, o);
    const double L = r0 / (std::cos(psi0) * std::cos(psi0));
    double off = 0;
    for (const Vec3& x : line.points) {
      const double r = physis::geometry::norm(x), c = std::hypot(x[0], x[1]) / r;
      off = std::max(off, std::fabs(r / (c * c) - L) / L);
    }
    s.check("curves: a dipole's field line from 30° lies on r = L cos² ψ, by RK4 with steps of 0.01 r", off, 1e-9);
    const double end_back = physis::geometry::norm(line.points.front()) - A, end_ahead = physis::geometry::norm(line.points.back()) - A;
    s.require("curves: it ends at both ends at the edge of the expansion's base, within 2⁻¹² of a step above the sphere",
              line.backward.stop == pf::CurveStop::edge && line.forward.stop == pf::CurveStop::edge && end_back >= 0 && end_ahead >= 0 &&
                  end_back < 0.01 * A * std::pow(2.0, -11) && end_ahead < 0.01 * A * std::pow(2.0, -11) &&
                  line.forward.reason == "below the expansion's reference sphere, inside its sources");
    s.require("curves: the two ends lie in opposite hemispheres", line.points.front()[2] * line.points.back()[2] < 0);

    const Vec3 omega{0, 0, 2 * std::numbers::pi};  // one turn a second
    const auto rotation = pf::analytic<Vec3>(described_as("rotation", pf::FibreKind::vector, physis::units::dim::velocity), [omega](double, const auto& x) {
      using S = std::decay_t<decltype(x[0])>;
      return std::array<S, 3>{-omega[2] * x[1], omega[2] * x[0], S(0.0)};
    });
    pf::CurveOptions st;
    st.along_unit = false;
    st.step = 1e-3;
    st.direction = pf::CurveOptions::Direction::forward;
    st.max_steps = 1000;
    const auto turn = pf::integral_curve(rotation, {1, 0, 0.5}, 0, st);
    // RK4 maps a rotation by its stability function R(iθ) = 1 − θ²/2 + θ⁴/24 + i(θ − θ³/6), θ = ωh: after
    // N steps the radius is |R|^N and the angle N arg R (Hairer, Nørsett and Wanner 1993, §IV.2).
    const double theta = omega[2] * st.step;
    const double re = 1 - theta * theta / 2 + std::pow(theta, 4) / 24, im = theta - std::pow(theta, 3) / 6;
    const double N = 1000, radius_rk = std::pow(std::hypot(re, im), N), angle_rk = std::remainder(N * std::atan2(im, re), 2 * std::numbers::pi);
    const Vec3 end = turn.points.back();
    s.check("curves: a streamline of a rigid rotation after 1000 steps has RK4's radius |R(iωh)|^N", std::hypot(end[0], end[1]) - radius_rk, 1e-13);
    s.check("curves: and RK4's angle N arg R(iωh), a phase error of N(ωh)⁵/120 short of the turn", std::atan2(end[1], end[0]) - angle_rk, 1e-12, "rad");
    s.require("curves: on its plane", end[2] == 0.5);
    s.require("curves: ended by its count of steps", turn.forward.stop == pf::CurveStop::steps && turn.points.size() == 1001);
    const auto none = pf::integral_curve(spectral_dipole, {0.5 * A, 0, 0}, 0);
    s.require("curves: a seed without a value gives no curve, and says why", !none.seeded && none.points.empty() && none.seed_quality == pf::Quality::outside);
  }

  // ---- versions ---------------------------------------------------------------------------------------------------------------
  {
    using U = pf::Uniform<double>;
    pf::History<U> h(described_as("evolving", pf::FibreKind::scalar, physis::units::dim::temperature));
    h.commit(0, U(described_as("v0", pf::FibreKind::scalar, {}), 1.0));
    h.commit(10, U(described_as("v1", pf::FibreKind::scalar, {}), 2.0));
    const auto before = h.at({-1, {}}), between = h.at({5, {}}), at = h.at({10, {}});
    bool threw = false;
    try { h.commit(5, U(described_as("late", pf::FibreKind::scalar, {}), 3.0)); } catch (const std::invalid_argument&) { threw = true; }
    s.require("versions: before the first commit missing; between, the earlier version with its note; at a commit, that version; a commit back in time refused",
              before.quality == pf::Quality::missing && between.value == 1.0 && !between.note.empty() && at.value == 2.0 && at.note.empty() && threw &&
                  h.versions() == 2 && pf::capabilities_of<pf::History<U>>().has(pf::Capability::versions));
  }

  // ---- sampled fields and their reconstruction ------------------------------------------------------------------------------
  {
    const pf::VolumeGrid g{{0.5, -1, 2}, physis::geometry::rotation_matrix(physis::geometry::Axis::y, 0.3), {4, 6, 3}, {5, 4, 3}};
    auto affine = [](const Vec3& x) { return 2 + 0.5 * x[0] - 1.5 * x[1] + 0.25 * x[2]; };
    std::vector<double> v(g.size());
    for (std::size_t k = 0; k < g.size(); k++) v[k] = affine(g.point(k));
    const pf::Sampled<double> f(described_as("sampled", pf::FibreKind::scalar, {}), g, pf::Samples<double>::all_valid(v));
    double worst = 0;
    for (int i = 0; i < 50; i++) {
      // Between the outer centres: offsets within (n − 1)/n of half an edge.
      const Vec3 o{rng.uniform(-0.39, 0.39) * 4, rng.uniform(-0.37, 0.37) * 6, rng.uniform(-0.33, 0.33) * 3};
      const Vec3 x = g.from_axes(o);
      worst = std::max(worst, std::fabs(f.at({0, x}).value - affine(x)));
    }
    s.check("sampled: trilinear reconstruction reproduces an affine field between the cells' centres", worst, 1e-13);
    const auto edge = f.at({0, g.from_axes({1.95, 0, 0})}), beyond = f.at({0, g.from_axes({2.05, 0, 0})});
    s.require("sampled: in the outer half of a boundary cell, held, with its note; beyond the box, outside",
              edge.ok() && !edge.note.empty() && beyond.quality == pf::Quality::outside);
    std::vector<pf::Evaluation<double>> e(g.size());
    for (std::size_t k = 0; k < g.size(); k++) e[k] = pf::valid(v[k]);
    e[g.index(2, 2, 1)] = pf::refused<double>(pf::Quality::invalid, "a refused sample");
    const pf::Sampled<double> holed(described_as("holed", pf::FibreKind::scalar, {}), g, pf::samples_of(std::span<const pf::Evaluation<double>>(e)));
    const auto near_hole = holed.at({0, g.from_axes(g.offset(2, 2, 1))});
    s.require("sampled: a point whose surrounding samples include a refused one has none, with its reason",
              near_hole.quality == pf::Quality::invalid && near_hole.reason == "a refused sample");
  }

  // ---- fields behind one type ----------------------------------------------------------------------------------------------
  {
    const pf::VectorField erased = rough;
    const pf::VectorField formula = dipole;
    const auto caps = erased.descriptor().capabilities;
    s.require("any: an erased expansion keeps its capabilities (point, sphere grid; no exact derivative), a formula its derivative",
              caps.has(pf::Capability::point) && caps.has(pf::Capability::sphere_grid) && !caps.has(pf::Capability::jacobian) &&
                  formula.descriptor().capabilities.has(pf::Capability::jacobian) && !erased.differentiable() && formula.differentiable() &&
                  erased.jacobian({0, {2 * A, 0, 0}}).quality == pf::Quality::invalid);
    const pf::SphereGrid grid{{}, 1.3 * A, 9, 18};
    const auto a = pf::sample(erased, grid, 0), b = pf::sample(rough, grid, 0);
    double gap = 0;
    for (std::size_t k = 0; k < grid.size(); k++) gap = std::max(gap, rel(a.values[k], b.values[k]));
    s.check("any: an erased field samples a sphere as its field does", gap, 0);
    std::vector<pf::Query<Vec3>> q{{0, {2 * A, 0, 0}}, {0, {0, 0.5 * A, 0}}};
    std::vector<pf::Evaluation<Vec3>> out(2);
    erased.at_many(q, out);
    s.require("any: in a batch, each answer as the field gives it", out[0].ok() && out[1].quality == pf::Quality::outside);
    const auto names = pf::names(caps);
    s.require("capabilities: named", names.size() == 2 && names[0] == "point" && names[1] == "sphere grid");
  }
  return s.done();
}
