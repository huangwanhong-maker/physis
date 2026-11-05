// Reaction–diffusion with kinetics of any form, on a graph or a periodic grid (decision 0022).
//
// reaction_diffusion.hpp's ReactionDiffusion holds networks of mass action. Many media of chemistry are
// modelled by reduced kinetics that are not mass action (the Oregonator's rational terms; enzymes saturating by
// Michaelis and Menten): here the local law is any callable kinetics(c, dc) of the s species at a vertex, or
// kinetics(v, c, dc) that may differ from vertex to vertex (a pattern of light on a photosensitive medium, cells
// of different recipes), and
//
//   ∂c/∂t = R(c) − D L c,
//
// with L a graph Laplacian (positive semidefinite, as path_laplacian writes it) given by sparse rows and D the
// species' diffusivities. periodic_grid writes the five-point Laplacian of an nx × ny torus of spacing h, whose
// mode e^{i(kx + ly)} has eigenvalue (4/h²)(sin²(kh/2) + sin²(lh/2)); grid writes the same grid with no-flux walls
// (a dish), whose cosine modes have the same eigenvalues at k = πm/(nx h), l = πn/(ny h). The step is the classical Runge–Kutta
// method; diffusive_step is its stability bound for the diffusion alone (2.78 over D_max times the largest row
// sum), the kinetics' own stiffness being the caller's to respect. For stiff kinetics, step_split splits the
// diffusion from the kinetics and takes the kinetics implicitly (an L-stable SDIRK), so that only the diffusion
// bounds the step. Checks: tests/check_field_kinetics.cpp.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace physis::chemistry {

struct SparseLaplacian {
  std::size_t vertices = 0;
  std::vector<std::vector<std::pair<std::size_t, double>>> rows;  // (column, entry), the diagonal included

  static SparseLaplacian periodic_grid(std::size_t nx, std::size_t ny, double spacing) {
    if (nx < 3 || ny < 3 || !(spacing > 0)) throw std::invalid_argument("physis::chemistry::SparseLaplacian::periodic_grid: at least 3 × 3 and a positive spacing");
    SparseLaplacian L;
    L.vertices = nx * ny;
    L.rows.resize(L.vertices);
    const double w = 1 / (spacing * spacing);
    for (std::size_t j = 0; j < ny; j++)
      for (std::size_t i = 0; i < nx; i++) {
        const std::size_t v = j * nx + i;
        L.rows[v] = {{v, 4 * w}, {j * nx + (i + 1) % nx, -w}, {j * nx + (i + nx - 1) % nx, -w}, {((j + 1) % ny) * nx + i, -w}, {((j + ny - 1) % ny) * nx + i, -w}};
      }
    return L;
  }
  // The five-point Laplacian of an nx × ny grid of spacing h with no-flux walls (a dish): the grid's graph without
  // wrapping, so that a boundary vertex has fewer neighbours and nothing crosses the walls.
  static SparseLaplacian grid(std::size_t nx, std::size_t ny, double spacing) {
    if (nx < 2 || ny < 2 || !(spacing > 0)) throw std::invalid_argument("physis::chemistry::SparseLaplacian::grid: at least 2 × 2 and a positive spacing");
    std::vector<std::tuple<std::size_t, std::size_t, double>> edges;
    const double w = 1 / (spacing * spacing);
    for (std::size_t j = 0; j < ny; j++)
      for (std::size_t i = 0; i < nx; i++) {
        if (i + 1 < nx) edges.push_back({j * nx + i, j * nx + i + 1, w});
        if (j + 1 < ny) edges.push_back({j * nx + i, (j + 1) * nx + i, w});
      }
    return graph(nx * ny, edges);
  }
  // The Laplacian of a weighted graph: edges (a, b, w).
  static SparseLaplacian graph(std::size_t vertices, const std::vector<std::tuple<std::size_t, std::size_t, double>>& edges) {
    SparseLaplacian L;
    L.vertices = vertices;
    L.rows.resize(vertices);
    std::vector<double> degree(vertices, 0.0);
    for (const auto& [a, b, w] : edges) {
      if (a >= vertices || b >= vertices) throw std::invalid_argument("physis::chemistry::SparseLaplacian::graph: edges between vertices");
      L.rows[a].push_back({b, -w}), L.rows[b].push_back({a, -w});
      degree[a] += w, degree[b] += w;
    }
    for (std::size_t v = 0; v < vertices; v++) L.rows[v].push_back({v, degree[v]});
    return L;
  }
  double largest_row_sum() const {
    double m = 0;
    for (const auto& r : rows) {
      double s = 0;
      for (const auto& [c, v] : r) s += std::fabs(v);
      m = std::max(m, s);
    }
    return m;
  }
};

class FieldKinetics {
 public:
  using Kinetics = std::function<void(const double* c, double* dc)>;                           // the same at every vertex
  using VertexKinetics = std::function<void(std::size_t vertex, const double* c, double* dc)>;  // vertex by vertex

  FieldKinetics(std::size_t species, Kinetics kinetics, SparseLaplacian laplacian, std::vector<double> diffusion)
      : FieldKinetics(species, VertexKinetics([k = std::move(kinetics)](std::size_t, const double* c, double* dc) { k(c, dc); }), std::move(laplacian), std::move(diffusion)) {}
  FieldKinetics(std::size_t species, VertexKinetics kinetics, SparseLaplacian laplacian, std::vector<double> diffusion)
      : s_(species), R_(std::move(kinetics)), L_(std::move(laplacian)), D_(std::move(diffusion)) {
    if (species == 0 || D_.size() != species) throw std::invalid_argument("physis::chemistry::FieldKinetics: one diffusivity per species");
    c_.assign(s_ * L_.vertices, 0.0);
  }

  std::size_t species() const { return s_; }
  std::size_t vertices() const { return L_.vertices; }
  double time() const { return t_; }
  // Concentrations laid out vertex-major: c[v * species + k].
  const std::vector<double>& state() const { return c_; }
  std::vector<double>& state() { return c_; }
  double at(std::size_t vertex, std::size_t k) const { return c_[vertex * s_ + k]; }

  void rate(const std::vector<double>& c, std::vector<double>& dc) const {
    dc.assign(c.size(), 0.0);
    for (std::size_t v = 0; v < L_.vertices; v++) {
      R_(v, c.data() + v * s_, dc.data() + v * s_);
      for (const auto& [w, l] : L_.rows[v])
        for (std::size_t k = 0; k < s_; k++) dc[v * s_ + k] -= D_[k] * l * c[w * s_ + k];
    }
  }

  void step(double h) {
    const std::size_t n = c_.size();
    std::vector<double> k1, k2, k3, k4, y(n);
    rate(c_, k1);
    for (std::size_t i = 0; i < n; i++) y[i] = c_[i] + 0.5 * h * k1[i];
    rate(y, k2);
    for (std::size_t i = 0; i < n; i++) y[i] = c_[i] + 0.5 * h * k2[i];
    rate(y, k3);
    for (std::size_t i = 0; i < n; i++) y[i] = c_[i] + h * k3[i];
    rate(y, k4);
    for (std::size_t i = 0; i < n; i++) c_[i] += h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    t_ += h;
  }
  void run(double duration, double h) {
    const auto steps = static_cast<std::size_t>(std::llround(duration / h));
    for (std::size_t k = 0; k < steps; k++) step(h);
  }
  double diffusive_step() const {
    const double d = *std::max_element(D_.begin(), D_.end());
    return d > 0 ? 2.78 / (d * L_.largest_row_sum()) : INFINITY;
  }

  // A step for stiff kinetics: Strang's splitting, diffusion for h/2 by the classical Runge–Kutta method, the
  // kinetics for h at every vertex by Alexander's two-stage SDIRK (γ = 1 − √2/2; L-stable and of order two,
  // Alexander 1977, SIAM J. Numer. Anal. 14, 1006), each stage by a damped Newton's method on the vertex's s × s
  // system with the Jacobian by central differences, then diffusion for h/2 again. Second order overall; stable for
  // kinetics however stiff, the step bounded by the diffusion alone (diffusive_step()). Where a stage's Newton
  // iteration fails (kinetics with poles, such as the Oregonator's at u = −q, can throw an undamped step past them),
  // that vertex's kinetics are taken in two half steps instead, recursively, to a depth of 12. Species declared
  // concentrations (set_non_negative) keep Newton's iterates non-negative, so that it cannot converge to a root the
  // chemistry forbids, beyond a pole of the kinetics.
  void step_split(double h) {
    diffuse(0.5 * h);
    std::vector<double> y(s_);
    for (std::size_t v = 0; v < L_.vertices; v++) {
      double* x = c_.data() + v * s_;
      for (std::size_t k = 0; k < s_; k++) y[k] = x[k];
      react(v, y, h, 0);
      for (std::size_t k = 0; k < s_; k++) x[k] = y[k];
    }
    diffuse(0.5 * h);
    t_ += h;
  }
  // Declares the species concentrations, never negative (for step_split's Newton iterations).
  void set_non_negative(bool on = true) { non_negative_ = on; }
  void run_split(double duration, double h) {
    const auto steps = static_cast<std::size_t>(std::llround(duration / h));
    for (std::size_t k = 0; k < steps; k++) step_split(h);
  }

 private:
  // The kinetics alone at vertex v for a time h by the SDIRK, in place; halves on failure.
  void react(std::size_t v, std::vector<double>& y, double h, int depth) const {
    const double g = 1 - std::sqrt(0.5);
    std::vector<double> Y1(s_), Y2(s_), f1(s_), c(s_);
    bool ok = solve_stage(v, y, h * g, Y1);
    if (ok) {
      R_(v, Y1.data(), f1.data());
      for (std::size_t k = 0; k < s_; k++) c[k] = y[k] + h * (1 - g) * f1[k];
      ok = solve_stage(v, c, h * g, Y2);
    }
    if (ok) {
      y = Y2;
      return;
    }
    if (depth >= 12) throw std::runtime_error("physis::chemistry::FieldKinetics::step_split: Newton's method did not converge at a vertex, even in steps of h/4096");
    react(v, y, 0.5 * h, depth + 1);
    react(v, y, 0.5 * h, depth + 1);
  }
  // Diffusion alone for a time τ, by the classical Runge–Kutta method.
  void diffuse(double tau) {
    if (*std::max_element(D_.begin(), D_.end()) == 0) return;
    const std::size_t n = c_.size();
    auto rate_d = [&](const std::vector<double>& c, std::vector<double>& dc) {
      dc.assign(n, 0.0);
      for (std::size_t v = 0; v < L_.vertices; v++)
        for (const auto& [w, l] : L_.rows[v])
          for (std::size_t k = 0; k < s_; k++) dc[v * s_ + k] -= D_[k] * l * c[w * s_ + k];
    };
    std::vector<double> k1, k2, k3, k4, y(n);
    rate_d(c_, k1);
    for (std::size_t i = 0; i < n; i++) y[i] = c_[i] + 0.5 * tau * k1[i];
    rate_d(y, k2);
    for (std::size_t i = 0; i < n; i++) y[i] = c_[i] + 0.5 * tau * k2[i];
    rate_d(y, k3);
    for (std::size_t i = 0; i < n; i++) y[i] = c_[i] + tau * k3[i];
    rate_d(y, k4);
    for (std::size_t i = 0; i < n; i++) c_[i] += tau / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
  }
  // The residual Y − c − a R(Y) and its largest component.
  double residual(std::size_t v, const std::vector<double>& c, double a, const std::vector<double>& Y, std::vector<double>& r) const {
    R_(v, Y.data(), r.data());
    double norm = 0;
    for (std::size_t k = 0; k < s_; k++) {
      r[k] = Y[k] - c[k] - a * r[k];
      norm = std::max(norm, std::isfinite(r[k]) ? std::fabs(r[k]) : INFINITY);
    }
    return norm;
  }
  // Solves Y = c + a R(Y) at one vertex by a damped Newton's method from Y = c: a step that does not lower the
  // residual is halved, up to 20 times. Returns false if it does not converge in 40 iterations.
  bool solve_stage(std::size_t v, const std::vector<double>& c, double a, std::vector<double>& Y) const {
    Y = c;
    if (non_negative_)
      for (double& y : Y) y = std::max(y, 0.0);
    std::vector<double> r(s_), rp(s_), rm(s_), Yp(s_), Ym(s_), trial(s_), rt(s_);
    std::vector<double> J(s_ * s_);
    double norm = residual(v, c, a, Y, r);
    for (int it = 0; it < 40; it++) {
      double scale = 0;
      for (double y : Y) scale = std::max(scale, std::fabs(y));
      if (norm <= 1e-13 * (1 + scale)) return true;
      // I − a ∂R/∂Y by central differences.
      for (std::size_t j = 0; j < s_; j++) {
        const double d = 1e-7 * (1 + std::fabs(Y[j]));
        Yp = Y, Ym = Y;
        Yp[j] += d, Ym[j] -= d;
        R_(v, Yp.data(), rp.data());
        R_(v, Ym.data(), rm.data());
        for (std::size_t i = 0; i < s_; i++) J[i * s_ + j] = (i == j ? 1.0 : 0.0) - a * (rp[i] - rm[i]) / (2 * d);
      }
      // Gaussian elimination with partial pivoting on the small system J δ = −r.
      std::vector<double> b(s_);
      for (std::size_t k = 0; k < s_; k++) b[k] = -r[k];
      for (std::size_t col = 0; col < s_; col++) {
        std::size_t piv = col;
        for (std::size_t i = col + 1; i < s_; i++)
          if (std::fabs(J[i * s_ + col]) > std::fabs(J[piv * s_ + col])) piv = i;
        if (J[piv * s_ + col] == 0) return false;
        if (piv != col) {
          for (std::size_t j = 0; j < s_; j++) std::swap(J[col * s_ + j], J[piv * s_ + j]);
          std::swap(b[col], b[piv]);
        }
        for (std::size_t i = col + 1; i < s_; i++) {
          const double m = J[i * s_ + col] / J[col * s_ + col];
          for (std::size_t j = col; j < s_; j++) J[i * s_ + j] -= m * J[col * s_ + j];
          b[i] -= m * b[col];
        }
      }
      for (std::size_t i = s_; i-- > 0;) {
        double sum = b[i];
        for (std::size_t j = i + 1; j < s_; j++) sum -= J[i * s_ + j] * b[j];
        b[i] = sum / J[i * s_ + i];
      }
      // Damping: halve the step until the residual falls.
      double lambda = 1, trial_norm = INFINITY;
      for (int h = 0; h < 20; h++, lambda *= 0.5) {
        bool inside = true;
        for (std::size_t k = 0; k < s_; k++) trial[k] = Y[k] + lambda * b[k], inside = inside && (!non_negative_ || trial[k] >= 0);
        if (!inside) continue;
        trial_norm = residual(v, c, a, trial, rt);
        if (trial_norm < norm) break;
      }
      if (!(trial_norm < norm)) return false;
      if (non_negative_)
        for (double t : trial)
          if (t < 0) return false;
      Y = trial, r = rt, norm = trial_norm;
    }
    return false;
  }

  std::size_t s_;
  VertexKinetics R_;
  SparseLaplacian L_;
  std::vector<double> D_;
  std::vector<double> c_;
  double t_ = 0;
  bool non_negative_ = false;
};

// Stiff kinetics in a well-mixed volume: a single vertex without diffusion, stepped by step_split's L-stable SDIRK,
// so that a network whose fast reactions are many orders quicker than its slow ones (strand-displacement cascades)
// is integrated with steps set by the slow ones' accuracy alone. kinetics(c, dc) as for FieldKinetics.
inline std::vector<double> integrate_stiff(std::size_t species, FieldKinetics::Kinetics kinetics, std::vector<double> c, double duration, double h,
                                           bool non_negative = true) {
  FieldKinetics F(species, std::move(kinetics), SparseLaplacian::graph(1, {}), std::vector<double>(species, 0.0));
  F.state() = std::move(c);
  F.set_non_negative(non_negative);
  F.run_split(duration, h);
  return F.state();
}

}  // namespace physis::chemistry
