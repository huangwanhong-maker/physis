// Checks of flows of any dimension, Lyapunov spectra and Floquet multipliers (flows.hpp, floquet.hpp;
// decision 0022).
//
// - rk4_step is fourth order: on x'' = −x, halving the step divides the error at t = 10 by 16 (within 5%).
// - Lorenz (σ = 10, r = 28, b = 8/3): the spectrum's sum is the divergence −41/3, to the integrator's error
//   (1e−5 at h = 0.005), whatever the path;
//   the middle exponent of a flow is 0 (within 0.005); the largest is 0.9056 (Sprott 2003, Chaos and Time-Series
//   Analysis, appendix A), within 0.02 over 2000 time units; the tangent dynamics by differences equal the exact
//   ones (to 1e−4 in every exponent); Kaplan and Yorke's dimension 2.06 (2 + 0.9056/14.572).
// - Lorenz 96 with n = 8, F = 8: the sum is −n.
// - A linear flow x' = A x with A upper triangular and stable (so that the state stays representable): the
//   exponents are A's diagonal, sorted.
// - The logistic map at r = 4: ln 2 (within 0.005 over 10⁶ iterations); Hénon (a = 1.4, b = 0.3): the sum is
//   ln 0.3 exactly (to 1e−9), the largest 0.4192 (within 0.01; Sprott 2003).
// - Floquet: a constant diagonal A has multipliers e^{a_i T}; Liouville's formula det M = exp ∫ tr A for a damped
//   Mathieu oscillator; the first tongue of x'' + (δ + ε cos t) x = 0 at ε = 0.1 lies between
//   δ = 1/4 ± ε/2 − ε²/8 (each boundary within 2e−4); a damped tongue's tip ε = 4γ/ω₀ (within 3% at γ/ω₀ = 0.01),
//   its multiplier negative: a subharmonic response.

#include <physis/dynamics/floquet.hpp>
#include <physis/dynamics/flows.hpp>

#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>

using namespace physis::dynamics;

int main() {
  calx::verification::Suite s;

  // ---- the integrator's order ----
  {
    auto osc = [](double, const State& x, State& dx) { dx[0] = x[1], dx[1] = -x[0]; };
    auto error = [&](double h) {
      State x{1, 0};
      integrate(osc, 0.0, 10.0, x, h);
      return std::hypot(x[0] - std::cos(10.0), x[1] + std::sin(10.0));
    };
    const double e1 = error(0.05), e2 = error(0.025);
    s.check("rk4: halving the step divides the error by 16 (ratio / 16 − 1)", e1 / e2 / 16 - 1, 0.05);
  }

  // ---- Lorenz ----
  const Lorenz lorenz;
  {
    LyapunovOptions o;
    o.step = 0.005, o.transient = 20, o.time = 2000, o.renormalise_every = 10;
    const auto L = lyapunov_spectrum(lorenz, State{1, 1, 20}, o);
    const double sum = L.exponents[0] + L.exponents[1] + L.exponents[2];
    // The classical Runge–Kutta step's tangent map has determinant e^{h tr J} to O((h|λ₃|)⁵/120) a step:
    // about 3e−6 per unit time at h = 0.005, whatever the path.
    s.check("Lorenz: the spectrum's sum is the divergence −(σ + 1 + b) = −41/3, to the integrator's error", sum + 41.0 / 3, 1e-5);
    s.check("Lorenz: the middle exponent of a flow is 0", L.exponents[1], 0.005);
    s.check("Lorenz: the largest exponent 0.9056 (Sprott 2003)", L.exponents[0] - 0.9056, 0.02);
    s.check("Lorenz: Kaplan–Yorke dimension 2 + λ₁/|λ₃| ≈ 2.062", kaplan_yorke(L.exponents) - 2.062, 0.005);
    o.time = 200;
    const auto A = lyapunov_spectrum(lorenz, State{1, 1, 20}, o);
    const auto B = lyapunov_spectrum(lorenz, [&](double, const State& x, const State& v) { return lorenz.jacobian_vector(x, v); }, State{1, 1, 20}, o);
    double d = 0;
    for (int j = 0; j < 3; j++) d = std::max(d, std::fabs(A.exponents[j] - B.exponents[j]));
    s.check("Lorenz: tangent dynamics by differences equal the exact ones (largest difference of an exponent)", d, 1e-4);
  }
  {
    const Lorenz96 l96;
    State x(8, 8.0);
    x[0] += 0.01;
    LyapunovOptions o;
    o.step = 0.005, o.transient = 10, o.time = 100;
    const auto L = lyapunov_spectrum(l96, x, o);
    double sum = 0;
    for (double l : L.exponents) sum += l;
    s.check("Lorenz 96 (n = 8): the spectrum's sum is −n", sum + 8, 1e-6);
    s.require("Lorenz 96 (n = 8, F = 8) is chaotic: λ₁ > 0", L.exponents[0] > 0.1);
  }
  {
    auto linear = [](double, const State& x, State& dx) {
      dx[0] = -0.1 * x[0] + 2 * x[1] - x[2];
      dx[1] = -0.7 * x[1] + 5 * x[2];
      dx[2] = -1.9 * x[2];
    };
    LyapunovOptions o;
    o.step = 0.01, o.time = 2000;
    const auto L = lyapunov_spectrum(linear, State{1, 1, 1}, o);
    // The frame's start costs O(ln κ / T), κ the eigenvectors' condition: about 1e−3 at T = 2000.
    s.check("linear flow: the exponents are the diagonal −0.1, −0.7, −1.9 (largest error)",
            std::max({std::fabs(L.exponents[0] + 0.1), std::fabs(L.exponents[1] + 0.7), std::fabs(L.exponents[2] + 1.9)}), 0.005);
  }

  // ---- maps ----
  {
    LyapunovOptions o;
    o.transient = 1000, o.time = 1e6, o.renormalise_every = 1;
    const auto L = lyapunov_spectrum_map(Logistic{}, State{0.1234}, o);
    s.check("logistic map at r = 4: the exponent is ln 2", L.exponents[0] - std::numbers::ln2, 0.005);
    o.time = 2e5;
    const auto H = lyapunov_spectrum_map(Henon{}, State{0.1, 0.1}, o);
    s.check("Hénon: the sum is ln b = ln 0.3 exactly", H.exponents[0] + H.exponents[1] - std::log(0.3), 1e-9);
    s.check("Hénon: the largest exponent 0.4192 (Sprott 2003)", H.exponents[0] - 0.4192, 0.01);
  }

  // ---- Floquet ----
  {
    const double T = 2.0;
    const auto F = floquet([](double, calx::linalg::Matrix& a) { a(0, 0) = 0.4, a(0, 1) = 0, a(1, 0) = 0, a(1, 1) = -1.1; }, 2, T);
    double m0 = std::max(F.multipliers[0].real(), F.multipliers[1].real()), m1 = std::min(F.multipliers[0].real(), F.multipliers[1].real());
    s.check("Floquet: constant A has multipliers e^{a T} (relative error of the larger)", m0 / std::exp(0.4 * T) - 1, 1e-10);
    s.check("Floquet: and of the smaller", m1 / std::exp(-1.1 * T) - 1, 1e-10);
    const double gamma = 0.05, Omega = 2.0;
    const auto D = parametric_oscillator(gamma, 1.0, 0.3, Omega);
    const double T2 = 2 * std::numbers::pi / Omega;
    const double det = D.monodromy(0, 0) * D.monodromy(1, 1) - D.monodromy(0, 1) * D.monodromy(1, 0);
    s.check("Floquet: Liouville, det M = e^{−2γT} (relative)", det / std::exp(-2 * gamma * T2) - 1, 1e-10);
  }
  {
    // x'' + (δ + ε cos t) x = 0: w0² = δ, w1² = −ε, Ω = 1.
    const double eps = 0.1;
    auto unstable = [&](double delta) { return parametric_oscillator(0, delta, -eps, 1.0, 4000).spectral_radius > 1 + 1e-9; };
    s.require("Mathieu: unstable at δ = 1/4 (inside the first tongue)", unstable(0.25));
    s.require("Mathieu: stable at δ = 0.18 and 0.32 (outside it)", !unstable(0.18) && !unstable(0.32));
    auto boundary = [&](double lo, double hi) {  // lo stable, hi unstable or the reverse; bisect
      const bool ulo = unstable(lo);
      for (int k = 0; k < 40; k++) {
        const double mid = 0.5 * (lo + hi);
        (unstable(mid) == ulo ? lo : hi) = mid;
      }
      return 0.5 * (lo + hi);
    };
    s.check("Mathieu: the tongue's upper boundary δ = 1/4 + ε/2 − ε²/8", boundary(0.25, 0.32) - (0.25 + eps / 2 - eps * eps / 8), 2e-4);
    s.check("Mathieu: the tongue's lower boundary δ = 1/4 − ε/2 − ε²/8", boundary(0.18, 0.25) - (0.25 - eps / 2 - eps * eps / 8), 2e-4);
  }
  {
    // x'' + 2γx' + ω₀²(1 − ε cos 2ω₀ t) x = 0: threshold ε = 4γ/ω₀.
    const double w0 = 1.0, gamma = 0.01;
    auto radius = [&](double e) { return parametric_oscillator(gamma, w0 * w0, e * w0 * w0, 2 * w0, 4000).spectral_radius; };
    double lo = 0, hi = 0.2;
    for (int k = 0; k < 50; k++) {
      const double mid = 0.5 * (lo + hi);
      (radius(mid) > 1 ? hi : lo) = mid;
    }
    s.check("damped Mathieu: the tongue's tip ε = 4γ/ω₀ (relative)", 0.5 * (lo + hi) / (4 * gamma / w0) - 1, 0.03);
    const auto F = parametric_oscillator(gamma, w0 * w0, 1.5 * 4 * gamma, 2 * w0, 4000);
    const auto& m = std::abs(F.multipliers[0]) > std::abs(F.multipliers[1]) ? F.multipliers[0] : F.multipliers[1];
    s.require("damped Mathieu: above the tip the growing multiplier is negative (subharmonic)", m.real() < 0 && std::fabs(m.imag()) < 1e-9 && std::abs(m) > 1);
  }
  return s.done();
}
