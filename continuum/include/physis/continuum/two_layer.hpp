// Two-layer quasi-geostrophic flow on a doubly periodic plane, by the pseudo-spectral method
// (decision 0022).
//
// The model (Phillips 1954, Tellus 6, 273; Pedlosky 1987, Geophysical Fluid Dynamics, 2nd ed., §6.16;
// Vallis 2017, Atmospheric and Oceanic Fluid Dynamics, 2nd ed., §5.6 and §9.6). Two layers of depths H₁
// (upper) and H₂ on an f-plane or β-plane, with streamfunctions ψᵢ and potential vorticities
//
//   q₁ = ∇²ψ₁ + F₁ (ψ₂ − ψ₁),   q₂ = ∇²ψ₂ + F₂ (ψ₁ − ψ₂),   Fᵢ = f₀² / (g′ Hᵢ),
//
// perturbing a uniform zonal flow Uᵢ, whose PV gradients are Q_y1 = β + F₁ (U₁ − U₂), Q_y2 = β − F₂ (U₁ − U₂):
//
//   ∂ₜqᵢ + J(ψᵢ, qᵢ) + Uᵢ ∂ₓqᵢ + Q_yi ∂ₓψᵢ = − δᵢ₂ r ∇²ψ₂ + D(qᵢ),
//
// J(a, b) = aₓ b_y − a_y bₓ, r a linear drag on the lower layer (an Ekman layer), D a hyperviscous damping
// of the smallest resolved scales. This is the form of pyqg's two-layer model (Abernathey et al., pyqg,
// doi:10.5281/zenodo.50569), whose conventions it follows.
//
// Inversion. For each wavevector (k, l), K² = k² + l² > 0, the 2 × 2 system [−(K² + F₁), F₁; F₂, −(K² + F₂)] ψ̂ = q̂
// has determinant K²(K² + F₁ + F₂); the mean (K = 0) is held at zero.
//
// Discretisation. Fourier modes on an nx × ny grid (powers of two); products in physical space with the
// 2/3 rule (modes beyond 2/3 of the largest wavenumber in either direction are set to zero), so that the
// semi-discrete system is a Galerkin truncation and keeps the energy
//
//   E = ½ ⟨δ₁ |∇ψ₁|² + δ₂ |∇ψ₂|² + δ₁ F₁ (ψ₁ − ψ₂)²⟩,   δᵢ = Hᵢ / (H₁ + H₂) = F_{3−i} / (F₁ + F₂),
//
// exactly when U = 0, r = 0 and D = 0 (with β or without), and the potential enstrophy ½⟨δ₁ q₁² + δ₂ q₂²⟩
// when also β = 0; the time step, the classical Runge–Kutta method, keeps them to its own error.
// D multiplies q̂ after each step by exp(−ρ h (K/K_c)⁸), ρ the declared damping rate at the truncation
// wavenumber K_c, so that it acts only on the smallest resolved scales.
//
// A channel. With `channel`, the state is kept odd in y, ψ(x, −y) = −ψ(x, y): by the method of images the box of
// width ly is then a channel of width ly/2 between free-slip walls at y = 0 and y = ly/2 where ψ = 0, whose modes
// are sin(jπy/(ly/2)). Oddness is preserved by the equations (J of two odd fields is odd; every other term has
// constant coefficients), and is restored after each step so that rounding cannot seed the even modes, which a
// channel forbids and which may grow faster (the modes uniform across it, l = 0).
//
// Linear stability. A mode e^{i(kx + ly)} grows at the real part of the eigenvalues of the 2 × 2 matrix of
// the linearised equations, linear_growth(k, l). For β = 0, equal layers (F₁ = F₂ = F), r = 0 and
// U₁ = −U₂ = U, its growth rate is k U √((k_d² − K²)/(k_d² + K²)) for K < k_d, k_d² = 2F, and zero beyond:
// Phillips's instability, a short-wave cut-off at the deformation scale.
//
// Units are the caller's, consistently. Checks: tests/check_two_layer.cpp.
#pragma once

#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <stdexcept>
#include <vector>

#include <calx/random/stream.hpp>
#include <calx/transforms/fft.hpp>

namespace physis::continuum {

class TwoLayerQG {
 public:
  using Complex = std::complex<double>;

  struct Parameters {
    std::size_t nx = 64, ny = 64;  // powers of two
    double lx = 2 * std::numbers::pi, ly = 2 * std::numbers::pi;
    double F1 = 1, F2 = 1;      // f₀² / (g′ Hᵢ)
    double beta = 0;            // planetary PV gradient
    double U1 = 0, U2 = 0;      // the background zonal flows
    double drag = 0;            // r, on the lower layer
    double step = 0.01;         // h
    double hyperviscous_rate = 0;  // ρ, at the truncation wavenumber
    bool channel = false;          // keep the state odd in y: a channel of width ly/2 between walls
  };

  struct Growth {
    double rate = 0;       // the largest real part
    double frequency = 0;  // its imaginary part (ω, with e^{i(kx + ly − ωt)})
  };

  explicit TwoLayerQG(const Parameters& p) : p_(p) {
    if (!calx::transforms::power_of_two(p.nx) || !calx::transforms::power_of_two(p.ny) || p.nx < 8 || p.ny < 8)
      throw std::invalid_argument("physis::continuum::TwoLayerQG: nx and ny powers of two, at least 8");
    if (!(p.lx > 0) || !(p.ly > 0) || !(p.F1 >= 0) || !(p.F2 >= 0) || !(p.step > 0) || p.drag < 0 || p.hyperviscous_rate < 0)
      throw std::invalid_argument("physis::continuum::TwoLayerQG: positive sides and step, non-negative F, drag and damping");
    const std::size_t n = p.nx * p.ny;
    k_.resize(n), l_.resize(n), K2_.resize(n), keep_.resize(n), filter_.resize(n);
    const double kmax = std::numbers::pi * static_cast<double>(p.nx) / p.lx, lmax = std::numbers::pi * static_cast<double>(p.ny) / p.ly;
    const double kc = std::min(2.0 / 3 * kmax, 2.0 / 3 * lmax);
    for (std::size_t j = 0; j < p.ny; j++)
      for (std::size_t i = 0; i < p.nx; i++) {
        const std::size_t m = j * p.nx + i;
        k_[m] = wavenumber(i, p.nx, p.lx), l_[m] = wavenumber(j, p.ny, p.ly);
        K2_[m] = k_[m] * k_[m] + l_[m] * l_[m];
        keep_[m] = std::fabs(k_[m]) <= 2.0 / 3 * kmax + 1e-12 && std::fabs(l_[m]) <= 2.0 / 3 * lmax + 1e-12 && i != p.nx / 2 && j != p.ny / 2;
        const double ratio = std::sqrt(K2_[m]) / kc;
        filter_[m] = std::exp(-p.hyperviscous_rate * p.step * std::pow(ratio, 8));
      }
    q1_.assign(n, 0.0), q2_.assign(n, 0.0);
    mirror_.resize(n);
    for (std::size_t j = 0; j < p.ny; j++)
      for (std::size_t i = 0; i < p.nx; i++) mirror_[j * p.nx + i] = ((p.ny - j) % p.ny) * p.nx + i;
  }

  const Parameters& parameters() const { return p_; }
  double time() const { return t_; }
  std::size_t nx() const { return p_.nx; }
  std::size_t ny() const { return p_.ny; }

  // ---- state ----------------------------------------------------------------------------------------------------------
  // Grids are row by row, index j nx + i at ((i) lx/nx, (j) ly/ny).
  void set_potential_vorticity(const std::vector<double>& q1, const std::vector<double>& q2) {
    q1_ = forward(q1), q2_ = forward(q2);
    truncate(q1_), truncate(q2_);
  }
  void set_streamfunction(const std::vector<double>& psi1, const std::vector<double>& psi2) {
    auto a = forward(psi1), b = forward(psi2);
    for (std::size_t m = 0; m < a.size(); m++) {
      q1_[m] = -K2_[m] * a[m] + p_.F1 * (b[m] - a[m]);
      q2_[m] = -K2_[m] * b[m] + p_.F2 * (a[m] - b[m]);
    }
    truncate(q1_), truncate(q2_);
  }
  // Small random PV in both layers, independent of the grid's order, from a keyed stream.
  void perturb(double amplitude, std::uint64_t seed) {
    calx::random::Stream rng(seed, 0x51C0u, 0x2Au, 0);
    std::vector<double> a(p_.nx * p_.ny), b(a.size());
    for (std::size_t m = 0; m < a.size(); m++) a[m] = amplitude * rng.normal(), b[m] = amplitude * rng.normal();
    set_potential_vorticity(a, b);
  }

  std::vector<double> potential_vorticity(int layer) const { return inverse(layer == 1 ? q1_ : q2_); }
  std::vector<double> streamfunction(int layer) const {
    std::vector<Complex> a, b;
    invert(q1_, q2_, a, b);
    return inverse(layer == 1 ? a : b);
  }
  // The zonal and meridional velocities of a layer's perturbation, u = −ψ_y, v = ψ_x.
  std::vector<double> velocity(int layer, int component) const {
    std::vector<Complex> a, b;
    invert(q1_, q2_, a, b);
    const auto& psi = layer == 1 ? a : b;
    std::vector<Complex> d(psi.size());
    for (std::size_t m = 0; m < d.size(); m++) d[m] = component == 0 ? -Complex(0, l_[m]) * psi[m] : Complex(0, k_[m]) * psi[m];
    return inverse(d);
  }

  // ---- invariants ---------------------------------------------------------------------------------------------------------
  double kinetic_energy(int layer) const {
    std::vector<Complex> a, b;
    invert(q1_, q2_, a, b);
    const auto& psi = layer == 1 ? a : b;
    const double d = layer == 1 ? delta1() : delta2();
    double s = 0;
    for (std::size_t m = 0; m < psi.size(); m++) s += K2_[m] * std::norm(psi[m]);
    return 0.5 * d * s / norm2();
  }
  double potential_energy() const {
    std::vector<Complex> a, b;
    invert(q1_, q2_, a, b);
    double s = 0;
    for (std::size_t m = 0; m < a.size(); m++) s += std::norm(a[m] - b[m]);
    return 0.5 * delta1() * p_.F1 * s / norm2();
  }
  double energy() const { return kinetic_energy(1) + kinetic_energy(2) + potential_energy(); }
  double enstrophy() const {
    double s1 = 0, s2 = 0;
    for (std::size_t m = 0; m < q1_.size(); m++) s1 += std::norm(q1_[m]), s2 += std::norm(q2_[m]);
    return 0.5 * (delta1() * s1 + delta2() * s2) / norm2();
  }

  // ---- time ---------------------------------------------------------------------------------------------------------------
  void step() {
    const std::size_t n = q1_.size();
    const double h = p_.step;
    std::vector<Complex> k1a(n), k1b(n), k2a(n), k2b(n), k3a(n), k3b(n), k4a(n), k4b(n), ya(n), yb(n);
    tendency(q1_, q2_, k1a, k1b);
    for (std::size_t m = 0; m < n; m++) ya[m] = q1_[m] + 0.5 * h * k1a[m], yb[m] = q2_[m] + 0.5 * h * k1b[m];
    tendency(ya, yb, k2a, k2b);
    for (std::size_t m = 0; m < n; m++) ya[m] = q1_[m] + 0.5 * h * k2a[m], yb[m] = q2_[m] + 0.5 * h * k2b[m];
    tendency(ya, yb, k3a, k3b);
    for (std::size_t m = 0; m < n; m++) ya[m] = q1_[m] + h * k3a[m], yb[m] = q2_[m] + h * k3b[m];
    tendency(ya, yb, k4a, k4b);
    for (std::size_t m = 0; m < n; m++) {
      q1_[m] += h / 6 * (k1a[m] + 2.0 * k2a[m] + 2.0 * k3a[m] + k4a[m]);
      q2_[m] += h / 6 * (k1b[m] + 2.0 * k2b[m] + 2.0 * k3b[m] + k4b[m]);
      if (p_.hyperviscous_rate > 0) q1_[m] *= filter_[m], q2_[m] *= filter_[m];
    }
    if (p_.channel) make_odd(q1_), make_odd(q2_);
    t_ += h;
  }
  void run(std::size_t steps) {
    for (std::size_t s = 0; s < steps; s++) step();
  }

  // ---- linear stability ---------------------------------------------------------------------------------------------------
  Growth linear_growth(double k, double l) const {
    const double K2 = k * k + l * l;
    if (K2 == 0) return {};
    const double det = K2 * (K2 + p_.F1 + p_.F2);
    // ψ̂ = P q̂.
    const double P11 = -(K2 + p_.F2) / det, P12 = -p_.F1 / det, P21 = -p_.F2 / det, P22 = -(K2 + p_.F1) / det;
    const double Qy1 = p_.beta + p_.F1 * (p_.U1 - p_.U2), Qy2 = p_.beta - p_.F2 * (p_.U1 - p_.U2);
    const Complex ik(0, k);
    const Complex b11 = -ik * (p_.U1 + Qy1 * P11), b12 = -ik * Qy1 * P12;
    const Complex b21 = (-ik * Qy2 + p_.drag * K2) * P21, b22 = -ik * p_.U2 + (-ik * Qy2 + p_.drag * K2) * P22;
    const Complex tr = b11 + b22, dt = b11 * b22 - b12 * b21;
    const Complex disc = std::sqrt(0.25 * tr * tr - dt);
    const Complex e1 = 0.5 * tr + disc, e2 = 0.5 * tr - disc;
    const Complex e = e1.real() >= e2.real() ? e1 : e2;
    return {e.real(), -e.imag()};
  }

  // The wavenumbers of the grid.
  double kx(std::size_t i) const { return wavenumber(i, p_.nx, p_.lx); }
  double ky(std::size_t j) const { return wavenumber(j, p_.ny, p_.ly); }

 private:
  static double wavenumber(std::size_t i, std::size_t n, double length) {
    const double m = i <= n / 2 ? static_cast<double>(i) : static_cast<double>(i) - static_cast<double>(n);
    return 2 * std::numbers::pi * m / length;
  }
  double delta1() const { return (p_.F1 + p_.F2) > 0 ? p_.F2 / (p_.F1 + p_.F2) : 0.5; }
  double delta2() const { return (p_.F1 + p_.F2) > 0 ? p_.F1 / (p_.F1 + p_.F2) : 0.5; }
  double norm2() const { const double n = static_cast<double>(p_.nx * p_.ny); return n * n; }

  void truncate(std::vector<Complex>& a) const {
    for (std::size_t m = 0; m < a.size(); m++)
      if (!keep_[m] || K2_[m] == 0) a[m] = 0;
    if (p_.channel) make_odd(a);
  }
  // The odd part in y: â(k, l) ← (â(k, l) − â(k, −l)) / 2.
  void make_odd(std::vector<Complex>& a) const {
    for (std::size_t m = 0; m < a.size(); m++)
      if (mirror_[m] > m) {
        const Complex odd = 0.5 * (a[m] - a[mirror_[m]]);
        a[m] = odd, a[mirror_[m]] = -odd;
      } else if (mirror_[m] == m) {
        a[m] = 0;
      }
  }

  void invert(const std::vector<Complex>& q1, const std::vector<Complex>& q2, std::vector<Complex>& psi1, std::vector<Complex>& psi2) const {
    psi1.assign(q1.size(), 0.0), psi2.assign(q1.size(), 0.0);
    for (std::size_t m = 0; m < q1.size(); m++) {
      if (K2_[m] == 0) continue;
      const double det = K2_[m] * (K2_[m] + p_.F1 + p_.F2);
      psi1[m] = (-(K2_[m] + p_.F2) * q1[m] - p_.F1 * q2[m]) / det;
      psi2[m] = (-p_.F2 * q1[m] - (K2_[m] + p_.F1) * q2[m]) / det;
    }
  }

  // The two-dimensional transform: rows, then columns.
  void transform(std::vector<Complex>& a, bool inv) const {
    const std::size_t nx = p_.nx, ny = p_.ny;
    std::vector<Complex> line(nx);
    for (std::size_t j = 0; j < ny; j++) {
      for (std::size_t i = 0; i < nx; i++) line[i] = a[j * nx + i];
      calx::transforms::fft(line, inv);
      for (std::size_t i = 0; i < nx; i++) a[j * nx + i] = line[i];
    }
    line.resize(ny);
    for (std::size_t i = 0; i < nx; i++) {
      for (std::size_t j = 0; j < ny; j++) line[j] = a[j * nx + i];
      calx::transforms::fft(line, inv);
      for (std::size_t j = 0; j < ny; j++) a[j * nx + i] = line[j];
    }
  }
  std::vector<Complex> forward(const std::vector<double>& u) const {
    if (u.size() != p_.nx * p_.ny) throw std::invalid_argument("physis::continuum::TwoLayerQG: a grid of nx × ny values");
    std::vector<Complex> a(u.begin(), u.end());
    transform(a, false);
    return a;
  }
  std::vector<double> inverse(std::vector<Complex> a) const {
    transform(a, true);
    std::vector<double> u(a.size());
    for (std::size_t m = 0; m < a.size(); m++) u[m] = a[m].real();
    return u;
  }

  // −J(ψ, q) of one layer, transformed and truncated.
  std::vector<Complex> advection(const std::vector<Complex>& psi, const std::vector<Complex>& q) const {
    const std::size_t n = psi.size();
    std::vector<Complex> u(n), v(n), qx(n), qy(n);
    for (std::size_t m = 0; m < n; m++) {
      u[m] = -Complex(0, l_[m]) * psi[m], v[m] = Complex(0, k_[m]) * psi[m];
      qx[m] = Complex(0, k_[m]) * q[m], qy[m] = Complex(0, l_[m]) * q[m];
    }
    transform(u, true), transform(v, true), transform(qx, true), transform(qy, true);
    std::vector<Complex> j(n);
    for (std::size_t m = 0; m < n; m++) j[m] = -(u[m].real() * qx[m].real() + v[m].real() * qy[m].real());
    transform(j, false);
    truncate(j);
    return j;
  }

  void tendency(const std::vector<Complex>& q1, const std::vector<Complex>& q2, std::vector<Complex>& d1, std::vector<Complex>& d2) const {
    std::vector<Complex> psi1, psi2;
    invert(q1, q2, psi1, psi2);
    d1 = advection(psi1, q1), d2 = advection(psi2, q2);
    const double Qy1 = p_.beta + p_.F1 * (p_.U1 - p_.U2), Qy2 = p_.beta - p_.F2 * (p_.U1 - p_.U2);
    for (std::size_t m = 0; m < q1.size(); m++) {
      if (!keep_[m] || K2_[m] == 0) { d1[m] = d2[m] = 0; continue; }
      const Complex ik(0, k_[m]);
      d1[m] += -ik * p_.U1 * q1[m] - ik * Qy1 * psi1[m];
      d2[m] += -ik * p_.U2 * q2[m] - ik * Qy2 * psi2[m] + p_.drag * K2_[m] * psi2[m];
    }
  }

  Parameters p_;
  double t_ = 0;
  std::vector<double> k_, l_, K2_, filter_;
  std::vector<char> keep_;
  std::vector<std::size_t> mirror_;  // the index of (k, −l)
  std::vector<Complex> q1_, q2_;
};

}  // namespace physis::continuum
