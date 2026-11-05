// Checks of waves on a surface (waves.hpp; decision 0022).
//
// - Dispersion: deep water, ω² = g k (relative 1e−9 at kh = 40, σ = 0); shallow water, ω = k √(g h)
//   (relative 2e−3 at kh = 0.01); the group velocity of deep gravity waves is half the phase velocity, of
//   capillary waves three halves (relative 1e−6), and agrees with the derivative of ω by differences.
// - The basin: a cosine mode is an exact eigenvector of the grid, and the leapfrog turns it at ω_h exactly
//   (to 1e−9 after 2000 steps); the grid's frequency tends to the continuum's at second order (the error
//   falls by 4 within 5% when the grid is refined twice); damping multiplies a mode each step by
//   √((1 − γh)/(1 + γh)); a step beyond the stable limit is refused; with no damping the discrete energy
//   is conserved (relative 1e−9) and a symmetric start stays symmetric.

#include <physis/continuum/waves.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <numbers>

using namespace physis::continuum;

int main() {
  calx::verification::Suite s;
  const double g = 9.80665;

  {
    const double k = 4.0;
    s.check("dispersion: deep water ω² = g k (relative)", std::pow(gravity_capillary_frequency(k, 10.0, g, 0.0), 2) / (g * k) - 1, 1e-9);
    const double h = 0.01 / k;
    s.check("dispersion: shallow water ω = k √(g h) (relative)", gravity_capillary_frequency(k, h, g, 0.0) / (k * std::sqrt(g * h)) - 1, 2e-3);
    const double cg = gravity_capillary_group_velocity(k, 10.0, g, 0.0), cp = gravity_capillary_frequency(k, 10.0, g, 0.0) / k;
    s.check("group velocity: deep gravity waves at half the phase velocity (relative)", cg / (0.5 * cp) - 1, 1e-6);
    const double kc = 5000.0, sigma = 0.0728, rho = 998.2;
    const double cgc = gravity_capillary_group_velocity(kc, 10.0, 0.0, sigma, rho), cpc = gravity_capillary_frequency(kc, 10.0, 0.0, sigma, rho) / kc;
    s.check("group velocity: capillary waves at three halves the phase velocity (relative)", cgc / (1.5 * cpc) - 1, 1e-6);
    const double k0 = 300.0, dk = 1e-3;
    const double num = (gravity_capillary_frequency(k0 + dk, 0.004) - gravity_capillary_frequency(k0 - dk, 0.004)) / (2 * dk);
    s.check("group velocity: equals dω/dk by differences, water 4 mm deep at k = 300 /m (relative)", gravity_capillary_group_velocity(k0, 0.004) / num - 1, 1e-6);
    s.check("damping: 2 ν k² for water at k = 300 /m (1/s)", viscous_damping_rate(300.0, 1.0e-6) - 0.18, 1e-12, "1/s");
  }

  {
    Basin::Options o;
    o.nx = 48, o.ny = 32, o.lx = 0.6, o.ly = 0.4, o.speed = 0.3, o.damping = 0;
    Basin b(o);
    const std::size_t m = 3, n = 2;
    b.set(b.mode(m, n));
    const double w = b.grid_frequency(m, n);
    const std::size_t steps = 2000;
    for (std::size_t k = 0; k < steps; k++) b.step();
    const auto e = b.mode(m, n);
    double err = 0;
    const double c = std::cos(w * b.time());
    for (std::size_t k = 0; k < e.size(); k++) err = std::max(err, std::fabs(b.surface()[k] - c * e[k]));
    // The start by Taylor expansion is exact to O(h⁴ ω⁴); the error of the phase stays at that size.
    s.check("basin: a grid mode turns at the leapfrog's frequency (largest error after 2000 steps)", err, 1e-9);

    auto rel_error = [&](std::size_t scale) {
      Basin::Options q = o;
      q.nx *= scale, q.ny *= scale;
      Basin c(q);
      // The grid's spatial error alone: the eigenvalue against the continuum's.
      const double lam = c.grid_eigenvalue(m, n), cont = std::pow(c.mode_frequency(m, n) / q.speed, 2);
      return std::fabs(lam / cont - 1);
    };
    s.check("basin: the grid's eigenvalue converges at second order (error ratio / 4 − 1)", rel_error(1) / rel_error(2) / 4 - 1, 0.05);

    Basin::Options d = o;
    d.damping = 0.7;
    Basin bd(d);
    bd.set(bd.mode(m, n));
    // Compare the amplitude every period of the mode, after the start: the envelope's ratio per step.
    double amp0 = 0;
    for (double v : bd.surface()) amp0 = std::max(amp0, std::fabs(v));
    const double gh = d.damping * bd.step_size();
    const double factor = std::sqrt((1 - gh) / (1 + gh));
    // The modal amplitude, from the projection onto the mode and the quadrature of its two last values.
    auto project = [&](const std::vector<double>& u) {
      double a = 0, nn = 0;
      for (std::size_t k = 0; k < u.size(); k++) a += u[k] * e[k], nn += e[k] * e[k];
      return a / nn;
    };
    std::vector<double> a;
    for (std::size_t k = 0; k < 600; k++) {
      bd.step();
      a.push_back(project(bd.surface()));
    }
    // a_n = A r^n cos(θ n + φ): r² = (a_{n+1}² − a_n a_{n+2}) / (a_n² − a_{n−1} a_{n+1}).
    const std::size_t p = 500;
    const double r2 = (a[p + 1] * a[p + 1] - a[p] * a[p + 2]) / (a[p] * a[p] - a[p - 1] * a[p + 1]);
    s.check("basin: damping multiplies a mode each step by √((1 − γh)/(1 + γh)) (relative)", std::sqrt(r2) / factor - 1, 1e-9);
    (void)amp0;

    Basin::Options bad = o;
    bad.step = 1.01 / (o.speed * std::sqrt(1 / std::pow(o.lx / o.nx, 2) + 1 / std::pow(o.ly / o.ny, 2)));
    bool refused = false;
    try {
      Basin x(bad);
    } catch (const std::invalid_argument&) {
      refused = true;
    }
    s.require("basin: a step beyond the stable limit is refused", refused);

    // The leapfrog's conserved energy: ½|(η^{n+1} − η^n)/h|² + ½ c² ⟨∇η^{n+1}, ∇η^n⟩ (a staggered form).
    Basin be(o);
    std::vector<double> start(o.nx * o.ny, 0.0);
    for (std::size_t j = 0; j < o.ny; j++)
      for (std::size_t i = 0; i < o.nx; i++) {
        const double x = (i + 0.5) * be.dx() - 0.3, y = (j + 0.5) * be.dy() - 0.2;
        start[j * o.nx + i] = std::exp(-(x * x + y * y) / 0.002);
      }
    be.set(start);
    auto staggered_energy = [&](const std::vector<double>& before, const std::vector<double>& now) {
      std::vector<double> lap(now.size());
      be.laplacian(now, lap);
      double kin = 0, pot = 0;
      for (std::size_t k = 0; k < now.size(); k++) {
        const double v = (now[k] - before[k]) / be.step_size();
        kin += 0.5 * v * v;
        pot -= 0.5 * o.speed * o.speed * before[k] * lap[k];
      }
      return kin + pot;
    };
    // E_n = ½|(η^{n} − η^{n−1})/h|² − ½ c² ⟨η^{n−1}, ∇²η^{n}⟩ is exactly conserved by leapfrog (∇² symmetric).
    std::vector<double> prev = be.surface();
    be.step();
    const double E0 = staggered_energy(prev, be.surface());
    for (int k = 0; k < 1500; k++) {
      prev = be.surface();
      be.step();
    }
    const double E1 = staggered_energy(prev, be.surface());
    s.check("basin: the leapfrog's discrete energy is conserved without damping (relative)", E1 / E0 - 1, 1e-9);
    double asym = 0;
    for (std::size_t j = 0; j < o.ny; j++)
      for (std::size_t i = 0; i < o.nx; i++) asym = std::max(asym, std::fabs(be.at(i, j) - be.at(o.nx - 1 - i, o.ny - 1 - j)));
    s.check("basin: a start symmetric about the centre stays symmetric", asym, 1e-12);
  }
  return s.done();
}
