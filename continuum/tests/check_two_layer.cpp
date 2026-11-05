// Checks of two-layer quasi-geostrophic flow (two_layer.hpp; decision 0022).
//
// - Linear stability: with β = 0, equal layers, no drag and U₁ = −U₂ = U, the growth rate of the linearised
//   system equals Phillips's k U √((k_d² − K²)/(k_d² + K²)) for K < k_d, k_d² = 2F (relative 1e−12 at several
//   wavevectors), and vanishes for K > k_d; drag lowers it; with no shear every mode is neutral.
// - The time stepper: a single small mode of the unstable band grows at the linear rate (relative 1% from the
//   energy's growth over a few e-foldings, the step 1/100 of the e-folding time).
// - Conservation: with U = 0, no drag and no damping, the energy is kept to 1e−8 relative over 400 steps of a
//   turbulent start, with β = 0 and with β = 2; the potential enstrophy too when β = 0.
// - A channel: a state kept odd in y stays odd through turbulent steps, and ψ vanishes on its walls (1e−12).
// - Inversion: a streamfunction set and read back is itself (to 1e−12, on modes within the truncation).

#include <physis/continuum/two_layer.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <numbers>

using physis::continuum::TwoLayerQG;

int main() {
  calx::verification::Suite s;
  const double pi = std::numbers::pi;

  TwoLayerQG::Parameters p;
  p.nx = p.ny = 32;
  p.F1 = p.F2 = 25;  // k_d² = 50
  p.U1 = 0.5, p.U2 = -0.5;
  const TwoLayerQG model(p);
  {
    double worst = 0;
    for (double k : {1.0, 2.0, 3.0, 5.0})
      for (double l : {0.0, 1.0, 4.0}) {
        const double K2 = k * k + l * l, kd2 = 50;
        if (K2 >= kd2) continue;
        const double phillips = k * 0.5 * std::sqrt((kd2 - K2) / (kd2 + K2));
        worst = std::max(worst, std::fabs(model.linear_growth(k, l).rate / phillips - 1));
      }
    s.check("linear: Phillips's growth rate k U √((k_d² − K²)/(k_d² + K²)) (largest relative error)", worst, 1e-12);
    s.check("linear: beyond the deformation scale the mode is neutral", model.linear_growth(6, 4).rate, 1e-12);
    TwoLayerQG::Parameters q = p;
    q.drag = 0.5;
    s.require("linear: drag lowers the growth", TwoLayerQG(q).linear_growth(3, 1).rate < model.linear_growth(3, 1).rate);
    TwoLayerQG::Parameters r = p;
    r.U1 = r.U2 = 0.3, r.beta = 1.5;
    s.check("linear: no shear, every mode neutral (Rossby waves)", TwoLayerQG(r).linear_growth(3, 1).rate, 1e-12);
  }

  {
    TwoLayerQG::Parameters q = p;
    const double k = 3, l = 1;  // wavevector (3, 1) on the 2π box
    const double sigma = model.linear_growth(k, l).rate;
    q.step = 0.01 / sigma;
    TwoLayerQG m(q);
    std::vector<double> psi1(q.nx * q.ny), psi2(psi1.size());
    for (std::size_t j = 0; j < q.ny; j++)
      for (std::size_t i = 0; i < q.nx; i++) {
        const double x = 2 * pi * i / q.nx, y = 2 * pi * j / q.ny;
        psi1[j * q.nx + i] = 1e-6 * std::cos(k * x + l * y);
        psi2[j * q.nx + i] = 1e-6 * std::cos(k * x + l * y + 1.0);
      }
    m.set_streamfunction(psi1, psi2);
    m.run(300);  // the start's decaying part dies away
    const double E1 = m.energy(), t1 = m.time();
    m.run(300);
    const double E2 = m.energy(), t2 = m.time();
    const double measured = std::log(E2 / E1) / (2 * (t2 - t1));
    s.check("time stepper: a small mode grows at the linear rate (relative)", measured / sigma - 1, 0.01);
  }

  for (double beta : {0.0, 2.0}) {
    TwoLayerQG::Parameters q = p;
    q.U1 = q.U2 = 0, q.beta = beta, q.step = 0.002;
    TwoLayerQG m(q);
    m.perturb(5.0, 7);
    const double E0 = m.energy(), Z0 = m.enstrophy();
    m.run(400);
    s.check(beta == 0 ? "conservation: energy kept, β = 0 (relative drift)" : "conservation: energy kept, β = 2 (relative drift)", m.energy() / E0 - 1, 1e-8);
    if (beta == 0) s.check("conservation: potential enstrophy kept, β = 0 (relative drift)", m.enstrophy() / Z0 - 1, 1e-8);
  }

  {
    // A channel: an odd state stays odd through turbulent steps, and the modes uniform across it are absent.
    TwoLayerQG::Parameters c = p;
    c.channel = true, c.step = 0.002;
    TwoLayerQG m(c);
    m.perturb(5.0, 13);
    m.run(200);
    const auto psi = m.streamfunction(1);
    double odd_error = 0, size = 0, uniform = 0;
    for (std::size_t j = 0; j < c.ny; j++)
      for (std::size_t i = 0; i < c.nx; i++) {
        const double a = psi[j * c.nx + i], b = psi[((c.ny - j) % c.ny) * c.nx + i];
        odd_error = std::max(odd_error, std::fabs(a + b)), size = std::max(size, std::fabs(a));
      }
    for (std::size_t i = 0; i < c.nx; i++) uniform = std::max(uniform, std::fabs(psi[i]));  // ψ on the wall y = 0
    s.check("channel: an odd state stays odd (relative)", odd_error / size, 1e-12);
    s.check("channel: ψ vanishes on the wall y = 0 (relative)", uniform / size, 1e-12);
  }

  {
    TwoLayerQG m(p);
    m.perturb(1.0, 11);
    const auto a = m.streamfunction(1), b = m.streamfunction(2);
    TwoLayerQG n(p);
    n.set_streamfunction(a, b);
    const auto c = n.streamfunction(1);
    double err = 0, size = 0;
    for (std::size_t k = 0; k < a.size(); k++) err = std::max(err, std::fabs(a[k] - c[k])), size = std::max(size, std::fabs(a[k]));
    s.check("inversion: a streamfunction set and read back is itself (relative)", err / size, 1e-12);
  }
  return s.done();
}
