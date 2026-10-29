// Heat conduction on a cellulation, as a balance law in conservation form (decision 0012;
// docs/theory/fields.md §10; formalism Definition 5.2, Propositions 5.3 and 5.6).
//
// The law. The energy density u, an n-form, balances its flux J, an (n − 1)-form, and its sources s:
// ∂ₜu + dJ = s, closed by Fourier's law J = −k ⋆ dT with u = ρ c T (Fourier 1822). On a cellulation the
// state is the energy uᵢ of each top cell (J), the flux J_f the heat through each facet along its
// conormal (W), and the closure the two-point flux of discrete exterior calculus with a diagonal Hodge
// star: across an interior facet between the cell a behind it and the cell b ahead,
//
//   J_f = G_f (T_a − T_b),   G_f = A_f / (d_a/k_a + d_b/k_b),
//
// A_f the facet's measure and d the distances from the cells' centroids to it along its conormal: the
// conductances of the two half-cells in series, so that across a change of material the flux is
// continuous and an affine temperature in each material is reproduced exactly. Tᵢ = uᵢ/Cᵢ with the
// heat capacity Cᵢ = ρ c Vᵢ of the cell's material (the occupancy label field, formalism Definition 4.1).
// A boundary facet is insulated, held at a temperature (G = k_a A_f/d_a), given a heat flux density into
// the body, or exchanging with an outside temperature through a coefficient h (A_f/(1/h + d_a/k_a)).
//
// The semi-discrete system is C Ṫ = −L T + b, L = d G dᵀ symmetric positive semidefinite (positive
// definite once any facet is held or exchanging), b the boundary's and the sources' forcing.
//
// The step. A θ-step (calx integrators) gives the temperatures T* of the new time; the energy is then
// updated by the fluxes, u′ = u − h Σ_f σ J_f(T_θ) + h s with T_θ = θ T* + (1 − θ) T, and T′ = u′/C. In
// exact arithmetic T′ = T*, the θ-method itself. In floating point the update is conservative whatever
// the linear solve's tolerance: every interior flux enters the cell behind and the cell ahead with
// opposite signs, so that the total energy changes by the boundary's outflow and the sources alone, to
// rounding (Proposition 5.3), and between materials what leaves one enters the other (Proposition 5.6).
// The solve's tolerance bounds only how far T′ is from the θ-method's.
//
// The checks also hold the scheme's quadratic dissipation, ½ Σ Cᵢ Tᵢ² non-increasing for an insulated
// body without sources at θ ≥ ½ for every step and at θ = 0 for h λ_max ≤ 2, and its convergence of
// second order to the decaying sine mode of the heat equation.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <calx/integrators/theta.hpp>
#include <calx/linalg/conjugate_gradient.hpp>
#include <calx/linalg/sparse.hpp>
#include <physis/discretisation/geometric.hpp>

namespace physis::continuum {

namespace dz = physis::discretisation;
using geometry::Vec3;

struct Material {
  std::string name;
  double density = 0;        // ρ, kg/m³
  double specific_heat = 0;  // c, J/(kg K)
  double conductivity = 0;   // k, W/(m K)
};

struct Boundary {
  enum class Kind : std::uint8_t { insulated, temperature, flux, exchange };
  Kind kind = Kind::insulated;
  double value = 0;        // temperature: the facet's (K); flux: the heat flux density into the body (W/m²); exchange: the outside temperature (K)
  double coefficient = 0;  // exchange: h (W/(m² K))

  static Boundary insulated() { return {}; }
  static Boundary temperature(double T) { return {Kind::temperature, T, 0}; }
  static Boundary flux(double into_body_W_m2) { return {Kind::flux, into_body_W_m2, 0}; }
  static Boundary exchange(double h_W_m2_K, double outside_K) { return {Kind::exchange, outside_K, h_W_m2_K}; }
};

// The energy of a step or a run (J): the change of what is stored, what flowed out through the boundary,
// what the sources gave. In exact arithmetic storage + outflow = sources.
struct Budget {
  double storage = 0;
  double outflow = 0;
  double sources = 0;
  double residual() const { return storage + outflow - sources; }
};

class HeatConduction {
 public:
  using BoundaryOf = std::function<Boundary(std::size_t facet, const Vec3& centroid, const Vec3& outward)>;

  struct State {
    double t = 0;
    std::vector<double> energy;       // uᵢ, J (from 0 K)
    std::vector<double> temperature;  // Tᵢ = uᵢ/Cᵢ, K
  };
  struct Step {
    State state;
    Budget budget;
    std::vector<double> flux;      // the heat through each facet along its conormal over the step (W), the θ-weighted
    std::size_t iterations = 0;    // of the solve
    double solve_residual = 0;     // its relative residual
  };

  HeatConduction(dz::GeometricComplex K, std::vector<int> occupancy, std::vector<Material> materials, const BoundaryOf& boundary,
                 std::vector<double> sources_W = {})
      : K_(std::move(K)), occupancy_(std::move(occupancy)), materials_(std::move(materials)), sources_(std::move(sources_W)) {
    M_ = dz::measure(K_);
    F_ = dz::facets(K_, M_);
    const std::size_t n = M_.dimension, N = K_.topology.cells(n);
    if (occupancy_.size() != N) throw std::invalid_argument("physis::continuum::HeatConduction: one material for each top cell");
    if (sources_.empty()) sources_.assign(N, 0.0);
    if (sources_.size() != N) throw std::invalid_argument("physis::continuum::HeatConduction: one source for each top cell");
    for (int m : occupancy_)
      if (m < 0 || static_cast<std::size_t>(m) >= materials_.size()) throw std::invalid_argument("physis::continuum::HeatConduction: a cell of no material");
    for (const auto& m : materials_)
      if (!(m.density > 0 && m.specific_heat > 0 && m.conductivity > 0)) throw std::invalid_argument("physis::continuum::HeatConduction: materials of positive ρ, c and k");
    const std::vector<double> volume = M_.of(n);
    capacity_.resize(N);
    for (std::size_t i = 0; i < N; i++) capacity_[i] = material(i).density * material(i).specific_heat * volume[i];
    // The conductance of each facet, and the boundary's forcing.
    conductance_.assign(F_.cells.size(), 0.0);
    boundary_.assign(F_.cells.size(), Boundary{});
    forcing_ = sources_;
    std::vector<std::uint32_t> r, c;
    std::vector<double> v;
    auto entry = [&](std::size_t i, std::size_t j, double x) {
      r.push_back(static_cast<std::uint32_t>(i));
      c.push_back(static_cast<std::uint32_t>(j));
      v.push_back(x);
    };
    for (std::size_t i = 0; i < N; i++) entry(i, i, 0.0);  // every diagonal entry stored
    for (std::size_t f = 0; f < F_.cells.size(); f++) {
      const double A = F_.measure[f];
      if (!F_.boundary(f)) {
        const std::size_t a = F_.cells[f][0], b = F_.cells[f][1];
        const double G = A / (half(a, f) / material(a).conductivity + half(b, f) / material(b).conductivity);
        conductance_[f] = G;
        entry(a, a, G), entry(b, b, G), entry(a, b, -G), entry(b, a, -G);
        continue;
      }
      const auto [a, sign] = F_.bounded(f);
      const Vec3 outward{sign * F_.conormal[f][0], sign * F_.conormal[f][1], sign * F_.conormal[f][2]};
      const Boundary B = boundary(f, F_.centroid[f], outward);
      boundary_[f] = B;
      const double k = material(a).conductivity, d = half(a, f);
      switch (B.kind) {
        case Boundary::Kind::insulated: break;
        case Boundary::Kind::temperature: conductance_[f] = k * A / d; break;
        case Boundary::Kind::exchange:
          if (!(B.coefficient > 0)) throw std::invalid_argument("physis::continuum::HeatConduction: an exchange coefficient > 0");
          conductance_[f] = A / (1 / B.coefficient + d / k);
          break;
        case Boundary::Kind::flux: forcing_[a] += B.value * A; break;
      }
      if (B.kind == Boundary::Kind::temperature || B.kind == Boundary::Kind::exchange) {
        entry(a, a, conductance_[f]);
        forcing_[a] += conductance_[f] * B.value;
      }
    }
    // One entry for each (row, column): the contributions summed.
    std::vector<std::uint64_t> order;
    const auto raw = calx::linalg::compress<std::uint32_t>(N, N, r, c, &order);
    std::vector<std::uint32_t> rr, cc;
    std::vector<double> vv;
    for (std::size_t i = 0; i < N; i++) {
      std::vector<std::pair<std::uint32_t, double>> row;
      for (std::uint64_t e = raw.row[i]; e < raw.row[i + 1]; e++) row.emplace_back(raw.col[e], v[order[e]]);
      std::stable_sort(row.begin(), row.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
      for (std::size_t e = 0; e < row.size(); e++) {
        if (e > 0 && row[e].first == row[e - 1].first) {
          vv.back() += row[e].second;
          continue;
        }
        rr.push_back(static_cast<std::uint32_t>(i));
        cc.push_back(row[e].first);
        vv.push_back(row[e].second);
      }
    }
    L_ = calx::linalg::compress<std::uint32_t>(N, N, rr, cc);
    values_ = std::move(vv);  // already in row order
  }

  std::size_t cells() const { return capacity_.size(); }
  const dz::GeometricComplex& complex() const { return K_; }
  const dz::Measures& measures() const { return M_; }
  const dz::Facets& facets() const { return F_; }
  const std::vector<double>& capacity() const { return capacity_; }
  const std::vector<double>& conductance() const { return conductance_; }
  const Material& material(std::size_t cell) const { return materials_[static_cast<std::size_t>(occupancy_[cell])]; }
  const calx::linalg::Csr<std::uint32_t>& matrix() const { return L_; }
  const std::vector<double>& matrix_values() const { return values_; }
  const std::vector<double>& forcing() const { return forcing_; }

  State state(double t, std::vector<double> temperature) const {
    if (temperature.size() != cells()) throw std::invalid_argument("physis::continuum::HeatConduction::state: one temperature for each cell");
    State s{t, std::vector<double>(cells()), std::move(temperature)};
    for (std::size_t i = 0; i < cells(); i++) s.energy[i] = capacity_[i] * s.temperature[i];
    return s;
  }

  // The heat through each facet along its conormal (W).
  std::vector<double> flux(std::span<const double> T) const {
    std::vector<double> J(F_.cells.size(), 0.0);
    for (std::size_t f = 0; f < J.size(); f++) {
      if (!F_.boundary(f)) {
        J[f] = conductance_[f] * (T[F_.cells[f][0]] - T[F_.cells[f][1]]);
        continue;
      }
      const auto [a, sign] = F_.bounded(f);
      J[f] = sign * outflow_through(f, a, T);
    }
    return J;
  }

  // The rate of each cell's energy, −Σ_f σ J_f + s (W).
  std::vector<double> rate(std::span<const double> T) const { return rate_of(flux(T)); }

  // The heat leaving through the boundary (W).
  double outflow(std::span<const double> T) const {
    double q = 0;
    for (std::size_t f = 0; f < F_.cells.size(); f++)
      if (F_.boundary(f)) q += outflow_through(f, F_.bounded(f).first, T);
    return q;
  }

  // One conservative θ-step of size h.
  Step step(const State& s, double h, double theta = 0.5, double tolerance = 1e-12) const {
    const auto th = calx::integrators::theta_step(L_, std::span<const double>(values_), std::span<const double>(capacity_), std::span<const double>(s.temperature),
                                                  std::span<const double>(forcing_), h, theta, tolerance);
    std::vector<double> Tth(cells());
    for (std::size_t i = 0; i < cells(); i++) Tth[i] = theta * th.y[i] + (1 - theta) * s.temperature[i];
    const std::vector<double> J = flux(Tth);
    const std::vector<double> du = rate_of(J);
    Step out;
    out.state.t = s.t + h;
    out.state.energy.resize(cells());
    out.state.temperature.resize(cells());
    for (std::size_t i = 0; i < cells(); i++) {
      out.state.energy[i] = s.energy[i] + h * du[i];
      out.state.temperature[i] = out.state.energy[i] / capacity_[i];
      out.budget.storage += out.state.energy[i] - s.energy[i];
      out.budget.sources += h * sources_[i];
    }
    // A boundary's outflow (a heat flux into the body counts as negative outflow).
    for (std::size_t f = 0; f < J.size(); f++)
      if (F_.boundary(f)) out.budget.outflow += h * F_.bounded(f).second * J[f];
    out.flux = J;
    out.iterations = th.iterations;
    out.solve_residual = th.relative_residual;
    return out;
  }

  // The steady temperatures, L T = b (positive definite when a facet is held or exchanging).
  std::vector<double> steady(double tolerance = 1e-14) const {
    const auto R = calx::linalg::conjugate_gradient(L_, std::span<const double>(values_), std::span<const double>(forcing_), tolerance);
    if (!R.converged) throw std::runtime_error("physis::continuum::HeatConduction::steady: the solve did not converge (no facet held or exchanging?)");
    return R.x;
  }

  // A step of explicit Euler that is stable whatever the temperatures, by Gershgorin's theorem: every
  // eigenvalue of C⁻¹L lies in a disc about Lᵢᵢ/Cᵢ of radius Σⱼ |Lᵢⱼ|/Cᵢ ≤ Lᵢᵢ/Cᵢ, so λ_max ≤ maxᵢ 2 Lᵢᵢ/Cᵢ and
  // h ≤ minᵢ Cᵢ/Lᵢᵢ keeps h λ_max ≤ 2.
  double explicit_step_bound() const {
    double h = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < cells(); i++)
      for (std::uint64_t e = L_.row[i]; e < L_.row[i + 1]; e++)
        if (L_.col[e] == i && values_[e] > 0) h = std::min(h, capacity_[i] / values_[e]);
    return h;
  }

  // ½ Σ Cᵢ Tᵢ² (J K): the quadratic functional the scheme dissipates.
  double quadratic(std::span<const double> T) const {
    double e = 0;
    for (std::size_t i = 0; i < cells(); i++) e += 0.5 * capacity_[i] * T[i] * T[i];
    return e;
  }

 private:
  // The distance from a cell's centroid to a facet along the facet's conormal.
  double half(std::size_t cell, std::size_t f) const {
    using namespace geometry::operators;
    const double d = std::fabs(geometry::dot(F_.centroid[f] - M_.centroid[M_.dimension][cell], F_.conormal[f]));
    if (!(d > 0)) throw std::invalid_argument("physis::continuum::HeatConduction: a cell's centroid on its facet");
    return d;
  }
  // The heat leaving cell a through boundary facet f (W).
  double outflow_through(std::size_t f, std::size_t a, std::span<const double> T) const {
    const Boundary& B = boundary_[f];
    switch (B.kind) {
      case Boundary::Kind::insulated: return 0;
      case Boundary::Kind::flux: return -B.value * F_.measure[f];
      default: return conductance_[f] * (T[a] - B.value);
    }
  }
  std::vector<double> rate_of(const std::vector<double>& J) const {
    std::vector<double> du = sources_;
    for (std::size_t f = 0; f < J.size(); f++) {
      if (F_.cells[f][0] != dz::no_cell) du[F_.cells[f][0]] -= J[f];  // the conormal points out of it
      if (F_.cells[f][1] != dz::no_cell) du[F_.cells[f][1]] += J[f];  // and into this one
    }
    return du;
  }

  dz::GeometricComplex K_;
  std::vector<int> occupancy_;
  std::vector<Material> materials_;
  std::vector<double> sources_;
  dz::Measures M_;
  dz::Facets F_;
  std::vector<double> capacity_, conductance_, forcing_;
  std::vector<Boundary> boundary_;
  calx::linalg::Csr<std::uint32_t> L_;
  std::vector<double> values_;
};

}  // namespace physis::continuum
