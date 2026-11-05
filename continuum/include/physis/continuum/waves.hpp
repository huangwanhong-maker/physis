// Linear waves on a liquid surface: the dispersion of gravity–capillary waves on a layer of finite depth,
// their viscous damping, and a basin, the damped wave equation of a surface in a rectangle with reflecting
// walls (decision 0022).
//
// Dispersion. Small waves e^{i(kx − ωt)} on an inviscid layer of depth h, density ρ and surface tension σ
// obey ω² = (g k + σ k³ / ρ) tanh(k h) (Lamb 1932, Hydrodynamics, 6th ed., §228 for gravity and §267 for
// capillarity). Long waves (kh ≪ 1) travel at √(g h) without dispersion.
//
// Damping. In a deep layer of small kinematic viscosity ν a wave's amplitude decays at the rate 2 ν k²
// (Lamb, §348–349): the bulk dissipation of irrotational motion. Boundary layers at the walls and floor,
// and surface films, add to it; a caller who knows them declares a larger rate.
//
// The basin. A surface η(x, y, t) over [0, Lx] × [0, Ly] obeys
//
//   η_tt = c² ∇²η − 2 γ η_t + s(x, y, t),
//
// with c the long-wave speed, γ the amplitude's decay rate and s a forcing (m/s²), and reflecting walls
// (∂η/∂n = 0). On a cell-centred grid of nx × ny cells the walls are a mirror (η_{−1} = η_0), so that
// the cosines cos(mπ(i + ½)/nx) cos(nπ(j + ½)/ny) are exact eigenvectors of the discrete Laplacian with
// eigenvalues −λ, λ = (2 − 2 cos(mπ/nx))/dx² + (2 − 2 cos(nπ/ny))/dy², which tend to the continuum's
// π²((m/Lx)² + (n/Ly)²) as the grid is refined (second order). The time step is leapfrog with the damping
// centred:
//
//   (1 + γh) η^{n+1} = 2 η^n − (1 − γh) η^{n−1} + h² (c² ∇²η^n + s^n),
//
// so that a mode's amplitude is multiplied each step by √((1 − γh)/(1 + γh)) and, undamped, turns at
// ω_h with cos(ω_h h) = 1 − h² c² λ / 2: second order in time, stable exactly when h c √(1/dx² + 1/dy²)
// ≤ 1 (the step is refused otherwise). Checks: tests/check_waves.cpp.
#pragma once

#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace physis::continuum {

// ω (rad/s) of a gravity–capillary wave of wavenumber k (1/m) on a layer of depth h (m).
inline double gravity_capillary_frequency(double k, double depth, double g = 9.80665, double surface_tension = 0.0728, double density = 998.2) {
  return std::sqrt((g * k + surface_tension * k * k * k / density) * std::tanh(k * depth));
}

// dω/dk, the group velocity, by the derivative of ω² / (2ω).
inline double gravity_capillary_group_velocity(double k, double depth, double g = 9.80665, double surface_tension = 0.0728, double density = 998.2) {
  const double th = std::tanh(k * depth), sech2 = 1 - th * th;
  const double a = g * k + surface_tension * k * k * k / density;
  const double da = g + 3 * surface_tension * k * k / density;
  const double w = std::sqrt(a * th);
  return (da * th + a * depth * sech2) / (2 * w);
}

// The amplitude's decay rate (1/s) of a wave of wavenumber k in a deep layer of kinematic viscosity ν (m²/s).
inline double viscous_damping_rate(double k, double kinematic_viscosity) { return 2 * kinematic_viscosity * k * k; }

class Basin {
 public:
  struct Options {
    std::size_t nx = 64, ny = 64;
    double lx = 1, ly = 1;  // m
    double speed = 1;       // c, m/s
    double damping = 0;     // γ, 1/s
    double step = 0;        // h, s; 0 for 0.9 of the stable limit
  };

  explicit Basin(const Options& o) : o_(o) {
    if (o.nx < 2 || o.ny < 2 || !(o.lx > 0) || !(o.ly > 0) || !(o.speed > 0) || o.damping < 0)
      throw std::invalid_argument("physis::continuum::Basin: at least 2 × 2 cells, positive sides and speed, damping ≥ 0");
    dx_ = o.lx / static_cast<double>(o.nx), dy_ = o.ly / static_cast<double>(o.ny);
    const double limit = 1 / (o.speed * std::sqrt(1 / (dx_ * dx_) + 1 / (dy_ * dy_)));
    h_ = o.step > 0 ? o.step : 0.9 * limit;
    if (h_ > limit * (1 + 1e-12)) throw std::invalid_argument("physis::continuum::Basin: the step exceeds the stable limit h c √(1/dx² + 1/dy²) ≤ 1");
    now_.assign(o.nx * o.ny, 0.0);
    before_ = now_;
  }

  std::size_t nx() const { return o_.nx; }
  std::size_t ny() const { return o_.ny; }
  double dx() const { return dx_; }
  double dy() const { return dy_; }
  double step_size() const { return h_; }
  double time() const { return t_; }
  const Options& options() const { return o_; }
  // η at the present step, row by row (index j nx + i, cell centre ((i + ½) dx, (j + ½) dy)).
  const std::vector<double>& surface() const { return now_; }
  double at(std::size_t i, std::size_t j) const { return now_[j * o_.nx + i]; }

  // Sets η and η_t at the present time; the previous step is taken from their Taylor expansion.
  void set(const std::vector<double>& eta, const std::vector<double>& eta_dot = {}) {
    if (eta.size() != now_.size() || (!eta_dot.empty() && eta_dot.size() != now_.size())) throw std::invalid_argument("physis::continuum::Basin::set: one value per cell");
    now_ = eta;
    std::vector<double> lap(now_.size());
    laplacian(now_, lap);
    for (std::size_t k = 0; k < now_.size(); k++) {
      const double v = eta_dot.empty() ? 0.0 : eta_dot[k];
      before_[k] = now_[k] - h_ * v + 0.5 * h_ * h_ * (o_.speed * o_.speed * lap[k] - 2 * o_.damping * v);
    }
  }

  // One step, with the forcing s (m/s², one value per cell) at the present time, or none.
  void step(const std::vector<double>* forcing = nullptr) {
    if (forcing && forcing->size() != now_.size()) throw std::invalid_argument("physis::continuum::Basin::step: one forcing per cell");
    std::vector<double> lap(now_.size()), next(now_.size());
    laplacian(now_, lap);
    const double c2 = o_.speed * o_.speed, gh = o_.damping * h_;
    for (std::size_t k = 0; k < now_.size(); k++) {
      const double s = forcing ? (*forcing)[k] : 0.0;
      next[k] = (2 * now_[k] - (1 - gh) * before_[k] + h_ * h_ * (c2 * lap[k] + s)) / (1 + gh);
    }
    before_.swap(now_);
    now_.swap(next);
    t_ += h_;
  }

  // The discrete Laplacian with mirror walls.
  void laplacian(const std::vector<double>& u, std::vector<double>& out) const {
    const std::size_t nx = o_.nx, ny = o_.ny;
    for (std::size_t j = 0; j < ny; j++)
      for (std::size_t i = 0; i < nx; i++) {
        const double c = u[j * nx + i];
        const double w = i > 0 ? u[j * nx + i - 1] : c, e = i + 1 < nx ? u[j * nx + i + 1] : c;
        const double s = j > 0 ? u[(j - 1) * nx + i] : c, n = j + 1 < ny ? u[(j + 1) * nx + i] : c;
        out[j * nx + i] = (w - 2 * c + e) / (dx_ * dx_) + (s - 2 * c + n) / (dy_ * dy_);
      }
  }

  // The continuum's frequency of mode (m, n): c π √((m/Lx)² + (n/Ly)²).
  double mode_frequency(std::size_t m, std::size_t n) const {
    const double a = static_cast<double>(m) / o_.lx, b = static_cast<double>(n) / o_.ly;
    return o_.speed * std::numbers::pi * std::sqrt(a * a + b * b);
  }
  // The grid's eigenvalue λ of mode (m, n), and the leapfrog's frequency of it, undamped.
  double grid_eigenvalue(std::size_t m, std::size_t n) const {
    const double pi = std::numbers::pi;
    return (2 - 2 * std::cos(pi * static_cast<double>(m) / static_cast<double>(o_.nx))) / (dx_ * dx_) +
           (2 - 2 * std::cos(pi * static_cast<double>(n) / static_cast<double>(o_.ny))) / (dy_ * dy_);
  }
  double grid_frequency(std::size_t m, std::size_t n) const {
    return std::acos(1 - 0.5 * h_ * h_ * o_.speed * o_.speed * grid_eigenvalue(m, n)) / h_;
  }
  // The cosine mode (m, n) on the grid, of unit amplitude.
  std::vector<double> mode(std::size_t m, std::size_t n) const {
    std::vector<double> u(o_.nx * o_.ny);
    const double pi = std::numbers::pi;
    for (std::size_t j = 0; j < o_.ny; j++)
      for (std::size_t i = 0; i < o_.nx; i++)
        u[j * o_.nx + i] = std::cos(pi * static_cast<double>(m) * (static_cast<double>(i) + 0.5) / static_cast<double>(o_.nx)) *
                           std::cos(pi * static_cast<double>(n) * (static_cast<double>(j) + 0.5) / static_cast<double>(o_.ny));
    return u;
  }

 private:
  Options o_;
  double dx_ = 0, dy_ = 0, h_ = 0, t_ = 0;
  std::vector<double> now_, before_;
};

}  // namespace physis::continuum
