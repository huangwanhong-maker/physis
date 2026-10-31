// Reaction–diffusion on graphs: Turing's instability of a homogeneous state, the persistence of fronts
// on a chain (McKean's caricature and its pinning), and the integration of mass action with diffusion
// (decision 0016).
//
// Turing's instability (Turing 1952, Phil. Trans. R. Soc. B 237, 37; on networks Othmer and Scriven
// 1971, J. Theor. Biol. 32, 507, and Nakao and Mikhailov 2010, Nat. Phys. 6, 544). Species s = 1 … m
// with concentrations c_v,s at the vertices v of a graph obey dc_v/dt = f(c_v) − Σ_w L_vw D c_w, L the
// graph Laplacian (symmetric, positive semidefinite, zero row sums: L = degree − adjacency, or divided
// by a spacing squared) and D = diag(d). About a homogeneous steady state with Jacobian J, a
// perturbation along an eigenvector of L with eigenvalue Λ ≥ 0 evolves by J − Λ D, so it grows at the
// largest real part of the eigenvalues of J − Λ D: the dispersion relation, here evaluated at the
// graph's own Λ_α (calx's symmetric eigenproblem and its general eigenvalues). The state is Turing-
// unstable when it is stable without diffusion (Λ = 0) and some Λ_α > 0 grows.
//
// For two species in closed form (Murray 2003, Mathematical Biology II, 3rd ed., §2.3), with
// B = d_v J11 + d_u J22, P = d_u d_v, Δ = det J:
//   stable homogeneous ⇔ tr J < 0 and Δ > 0;
//   Turing-unstable ⇔ also B > 2 √(P Δ);
//   the unstable band of squared wavenumbers k² lies between the roots of h(k²) = P k⁴ − B k² + Δ;
//   at onset, the ratio d = d_v/d_u solves J11² d² + 2(2 J12 J21 − J11 J22) d + J22² = 0, whose roots
//   are d_± = (J11 J22 − 2 J12 J21 ± 2 √(−J12 J21 Δ)) / J11²: instability for d > d_+ when J11 > 0 (u
//   activates), for d < d_− when J22 > 0 (v activates), d_c = J22² / (−4 J12 J21) when J11 = 0;
//   the fastest-growing k² maximises λ(k²) = ½ (τ + √(τ² − 4h)), τ = tr J − (d_u + d_v) k². Setting
//   dλ/dk² = 0 in λ² − τλ + h = 0 gives λ = (B − 2P k²)/(d_u + d_v), and with it the quadratic
//   P(4P − s²) k⁴ + 2P(s tr J − 2B) k² + (B² − s B tr J + s² Δ) = 0, s = d_u + d_v (derived here; the
//   root in the band with B − 2P k² > 0). It lies below the onset wavenumber B/(2P) of the minimum of h.
//
// Fronts on a chain, and their pinning. McKean's caricature of a bistable medium (McKean 1970, Adv.
// Math. 4, 209), on a chain of sites coupled with strength D:
//
//   du_n/dt = D (u_{n+1} − 2u_n + u_{n−1}) − u_n + H(u_n − a),   a ∈ (0, 1),   H(x) = 1 for x ≥ 0, else 0,
//
// has the stable states 0 and 1. A front between them moves, one state invading the other, unless the
// sites' discreteness holds it: propagation failure (Keener 1987, SIAM J. Appl. Math. 47, 556; Fáth
// 1998, Physica D 116, 176). Derived here: a stationary decreasing front with u_n ≥ a for n ≤ 0 and
// u_n < a for n ≥ 1 solves, for n ≥ 1, the linear recurrence D(u_{n+1} − 2u_n + u_{n−1}) − u_n = 0,
// whose bounded solutions on n ≥ 0 are u_n = A μⁿ with μ the root in (0, 1) of D μ² − (2D + 1) μ + D = 0,
//
//   μ = (2D + 1 − √(1 + 4D)) / (2D) = 2D / (2D + 1 + √(1 + 4D));
//
// and v_n = 1 − u_n solves the same recurrence for n ≤ 0, v_n = C μ^{−n} on n ≤ 1. The two expressions
// agree at n = 0 and n = 1 exactly when A = 1/(1 + μ), C = μ/(1 + μ), so
//
//   u_n = μⁿ/(1 + μ) (n ≥ 0),   u_n = 1 − μ^{1−n}/(1 + μ) (n ≤ 0),   u_0 = 1/(1 + μ),   u_1 = μ/(1 + μ),
//
// and the front is consistent, H as assumed, exactly when u_1 < a ≤ u_0. Since 1 + μ = s(s − 1)/(2D)
// with s = √(1 + 4D) and s² − 1 = 4D,
//
//   a_−(D) = ½ (1 − 1/√(1 + 4D)),   a_+(D) = ½ (1 + 1/√(1 + 4D)):
//
// a front is pinned for a ∈ (a_−, a_+) (with H(0) = 1, also at a = a_+); for a < a_− the state 1
// invades, for a > a_+ the state 0. The interval shrinks to ½ as D → ∞, the continuum's single
// balanced value, and widens to (0, 1) as D → 0, uncoupled sites. From the step u_n = H(−n), the front
// for a inside never flips a site: the time derivative, antisymmetric about n = ½, stays of one sign on
// each side (a positive semigroup), so the sites rise or fall monotonically to the stationary front.
//
// Integration. ReactionDiffusion adds mass action of a Network (chemistry/kinetics.hpp; rates read as
// doubles) at each vertex to −Σ_w L_vw d_s c_w,s, advanced by the classical Runge–Kutta method with
// run-time sizes (runge_kutta_4, provisional here until calx has integrators of run-time size). RK4 is
// stable for the diffusion alone when h d_max λ_max(L) ≤ 2.78 (its stability interval on the negative
// real axis), λ_max(L) ≤ 2 max_v L_vv by Gershgorin's theorem: diffusive_step(). Linear invariants
// (the totals of conserved combinations, when L's columns sum to zero) are kept to rounding.
//
// Checks: tests/check_reaction_diffusion.cpp.

#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include <calx/linalg/dense.hpp>
#include <calx/linalg/eigen.hpp>
#include <calx/linalg/eigenvalues.hpp>
#include <physis/chemistry/kinetics.hpp>

namespace physis::chemistry {

// ---- Turing's instability ------------------------------------------------------------------------------------------------------

struct TuringModes {
  std::vector<double> laplacian;    // Λ_α ascending (calx symmetric_eigen of L)
  std::vector<double> growth;       // max Re eig(J − Λ_α D) (calx eigenvalues)
  std::size_t fastest = 0;          // α of the largest growth
  bool stable_homogeneous = false;  // max Re eig(J) < 0
  bool turing_unstable = false;     // stable at Λ = 0 and some Λ_α > 0 grows
};

// The largest real part of the eigenvalues of J − Λ D.
inline double turing_growth(const calx::linalg::Matrix& jacobian, std::span<const double> diffusion, double Lambda) {
  calx::linalg::Matrix M = jacobian;
  for (std::size_t s = 0; s < M.n; s++) M(s, s) -= Lambda * diffusion[s];
  double g = -std::numeric_limits<double>::infinity();
  for (const auto& z : calx::linalg::eigenvalues(M)) g = std::max(g, z.real());
  return g;
}

inline TuringModes turing_modes(const calx::linalg::Matrix& jacobian, std::span<const double> diffusion, const calx::linalg::Matrix& laplacian) {
  if (jacobian.n == 0 || diffusion.size() != jacobian.n) throw std::invalid_argument("physis::chemistry::turing_modes: one diffusivity for each species of the Jacobian");
  for (double d : diffusion)
    if (!(d >= 0 && std::isfinite(d))) throw std::invalid_argument("physis::chemistry::turing_modes: diffusivities ≥ 0");
  if (laplacian.n == 0) throw std::invalid_argument("physis::chemistry::turing_modes: a Laplacian of order ≥ 1");
  double scale = 0, asym = 0;
  for (std::size_t i = 0; i < laplacian.n; i++)
    for (std::size_t j = 0; j < laplacian.n; j++) scale = std::max(scale, std::fabs(laplacian(i, j))), asym = std::max(asym, std::fabs(laplacian(i, j) - laplacian(j, i)));
  if (asym > 1e-12 * scale) throw std::invalid_argument("physis::chemistry::turing_modes: a symmetric Laplacian");
  TuringModes out;
  out.laplacian = calx::linalg::symmetric_eigen(laplacian).values;
  out.growth.resize(out.laplacian.size());
  for (std::size_t a = 0; a < out.laplacian.size(); a++) {
    out.growth[a] = turing_growth(jacobian, diffusion, out.laplacian[a]);
    if (out.growth[a] > out.growth[out.fastest]) out.fastest = a;
  }
  out.stable_homogeneous = turing_growth(jacobian, diffusion, 0.0) < 0;
  out.turing_unstable = out.stable_homogeneous && out.growth[out.fastest] > 0;
  return out;
}

struct TwoSpeciesTuring {
  bool stable_homogeneous = false, unstable = false;
  double critical_ratio = 0;                         // d_c = d_v/d_u at onset (NaN when no ratio destabilises)
  double k2_low = 0, k2_high = 0, k2_fastest = 0;    // the unstable band and its fastest k² (0 unless unstable)
};

// λ(k²) = ½ (τ + √(τ² − 4h)): the larger real part of the two-species dispersion relation.
inline double turing_two_species_growth(double j11, double j12, double j21, double j22, double du, double dv, double k2) {
  const double tau = j11 + j22 - (du + dv) * k2;
  const double h = du * dv * k2 * k2 - (dv * j11 + du * j22) * k2 + (j11 * j22 - j12 * j21);
  const double disc = tau * tau - 4 * h;
  return disc >= 0 ? 0.5 * (tau + std::sqrt(disc)) : 0.5 * tau;
}

inline TwoSpeciesTuring turing_two_species(double j11, double j12, double j21, double j22, double du, double dv) {
  if (!(du > 0 && dv > 0 && std::isfinite(du) && std::isfinite(dv))) throw std::invalid_argument("physis::chemistry::turing_two_species: diffusivities > 0");
  TwoSpeciesTuring out;
  const double trace = j11 + j22, det = j11 * j22 - j12 * j21;
  const double B = dv * j11 + du * j22, P = du * dv;
  out.stable_homogeneous = trace < 0 && det > 0;
  out.unstable = out.stable_homogeneous && B > 0 && B * B > 4 * P * det;
  out.critical_ratio = std::numeric_limits<double>::quiet_NaN();
  if (out.stable_homogeneous && j12 * j21 < 0) {
    const double root = 2 * std::sqrt(-j12 * j21 * det), mid = j11 * j22 - 2 * j12 * j21;
    if (j11 > 0) out.critical_ratio = (mid + root) / (j11 * j11);
    else if (j11 < 0 && j22 > 0) out.critical_ratio = (mid - root) / (j11 * j11);
    else if (j11 == 0 && j22 > 0) out.critical_ratio = j22 * j22 / (-4 * j12 * j21);
  }
  if (!out.unstable) return out;
  const double disc = std::sqrt(B * B - 4 * P * det);
  out.k2_high = (B + disc) / (2 * P);
  out.k2_low = det / (P * out.k2_high);  // the product of the roots is Δ/P
  // The stationary points of λ(k²): P(4P − s²) x² + 2P(s tr J − 2B) x + (B² − s B tr J + s² Δ) = 0.
  const double s = du + dv, a2 = P * (4 * P - s * s), a1 = 2 * P * (s * trace - 2 * B), a0 = B * B - s * B * trace + s * s * det;
  std::vector<double> roots;
  if (a2 != 0) {
    const double D = a1 * a1 - 4 * a2 * a0;
    if (D >= 0) {
      const double q = -0.5 * (a1 + std::copysign(std::sqrt(D), a1));
      if (q != 0) roots.push_back(q / a2), roots.push_back(a0 / q);
    }
  } else if (a1 != 0) {
    roots.push_back(-a0 / a1);
  }
  const auto growth = [&](double x) { return turing_two_species_growth(j11, j12, j21, j22, du, dv, x); };
  double best = -std::numeric_limits<double>::infinity();
  for (double x : roots)
    if (x > out.k2_low && x < out.k2_high && B - 2 * P * x > 0 && growth(x) > best) best = growth(x), out.k2_fastest = x;
  if (!(best > -std::numeric_limits<double>::infinity())) {
    // A guard that the closed form does not need: golden-section search on the band.
    double lo = out.k2_low, hi = out.k2_high;
    const double g = 0.5 * (std::sqrt(5.0) - 1);
    for (int i = 0; i < 200 && hi - lo > 1e-15 * hi; i++) {
      const double x1 = hi - g * (hi - lo), x2 = lo + g * (hi - lo);
      if (growth(x1) < growth(x2)) lo = x1;
      else hi = x2;
    }
    out.k2_fastest = 0.5 * (lo + hi);
  }
  return out;
}

// The Laplacian of a chain of n vertices at the given spacing, free (Neumann) ends: L = (degree −
// adjacency)/spacing², eigenvalues (2 − 2 cos(πk/n))/spacing², k = 0 … n−1.
inline calx::linalg::Matrix path_laplacian(std::size_t n, double spacing) {
  if (n == 0 || !(spacing > 0 && std::isfinite(spacing))) throw std::invalid_argument("physis::chemistry::path_laplacian: n ≥ 1 vertices and a spacing > 0");
  calx::linalg::Matrix L(n);
  const double w = 1 / (spacing * spacing);
  for (std::size_t i = 0; i + 1 < n; i++) {
    L(i, i) += w, L(i + 1, i + 1) += w;
    L(i, i + 1) -= w, L(i + 1, i) -= w;
  }
  return L;
}

// ---- fronts on a chain: McKean's caricature ---------------------------------------------------------------------------------

struct Pinning {
  double lower = 0, upper = 0;  // a_−(D) = ½(1 − 1/√(1 + 4D)), a_+(D) = ½(1 + 1/√(1 + 4D))
  double decay = 0;             // μ ∈ [0, 1): the stationary front's ratio from one site to the next
};

inline Pinning mckean_pinning(double D) {
  if (!(D >= 0 && std::isfinite(D))) throw std::invalid_argument("physis::chemistry::mckean_pinning: a coupling D ≥ 0");
  const double s = std::sqrt(1 + 4 * D);
  return {0.5 * (1 - 1 / s), 0.5 * (1 + 1 / s), 2 * D / (2 * D + 1 + s)};
}

// The stationary front at site n: μⁿ/(1 + μ) for n ≥ 0, 1 − μ^{1−n}/(1 + μ) for n ≤ 0 (high for n ≤ 0).
inline double mckean_front(long long n, double D) {
  const double mu = mckean_pinning(D).decay;
  if (n >= 0) return std::pow(mu, static_cast<double>(n)) / (1 + mu);
  return 1 - std::pow(mu, static_cast<double>(1 - n)) / (1 + mu);
}

// The right-hand side of McKean's caricature on a chain with free (Neumann) ends.
inline std::vector<double> mckean_rate(std::span<const double> u, double D, double a) {
  const std::size_t n = u.size();
  std::vector<double> r(n);
  for (std::size_t i = 0; i < n; i++) {
    const double left = i > 0 ? u[i - 1] : u[i], right = i + 1 < n ? u[i + 1] : u[i];
    r[i] = D * (left - 2 * u[i] + right) - u[i] + (u[i] >= a ? 1.0 : 0.0);
  }
  return r;
}

// ---- integration --------------------------------------------------------------------------------------------------------------

// One step of the classical Runge–Kutta method for y′ = f(y), with run-time sizes; f(y) returns y′.
// Provisional here, until calx's integrators take run-time sizes.
template <class F>
std::vector<double> runge_kutta_4(F&& f, std::span<const double> y, double h) {
  const std::size_t n = y.size();
  const std::vector<double> k1 = f(y);
  std::vector<double> z(n);
  for (std::size_t i = 0; i < n; i++) z[i] = y[i] + 0.5 * h * k1[i];
  const std::vector<double> k2 = f(std::span<const double>(z));
  for (std::size_t i = 0; i < n; i++) z[i] = y[i] + 0.5 * h * k2[i];
  const std::vector<double> k3 = f(std::span<const double>(z));
  for (std::size_t i = 0; i < n; i++) z[i] = y[i] + h * k3[i];
  const std::vector<double> k4 = f(std::span<const double>(z));
  for (std::size_t i = 0; i < n; i++) z[i] = y[i] + h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
  return z;
}

class ReactionDiffusion {
 public:
  ReactionDiffusion(const Network& network, const calx::linalg::Matrix& laplacian, std::vector<double> diffusion)
      : species_(network.species.size()), vertices_(laplacian.n), diffusion_(std::move(diffusion)) {
    if (vertices_ == 0 || laplacian.a.size() != vertices_ * vertices_) throw std::invalid_argument("physis::chemistry::ReactionDiffusion: a square Laplacian of order ≥ 1");
    if (diffusion_.size() != species_) throw std::invalid_argument("physis::chemistry::ReactionDiffusion: one diffusivity for each species");
    for (double d : diffusion_)
      if (!(d >= 0 && std::isfinite(d))) throw std::invalid_argument("physis::chemistry::ReactionDiffusion: diffusivities ≥ 0");
    for (const auto& r : network.reactions) {
      Cached c{r.rate.to_double(), {}, {}};
      for (const auto& [i, a] : r.reactants) c.reactants.push_back({i, a});
      for (const auto& [i, b] : r.products) c.products.push_back({i, b});
      reactions_.push_back(std::move(c));
    }
    rows_.resize(vertices_);
    for (std::size_t v = 0; v < vertices_; v++) {
      for (std::size_t w = 0; w < vertices_; w++)
        if (laplacian(v, w) != 0) rows_[v].push_back({w, laplacian(v, w)});
      diagonal_ = std::max(diagonal_, laplacian(v, v));
    }
  }

  // dc/dt, c laid out vertex-major: c[v * species + s].
  std::vector<double> rate(std::span<const double> c) const {
    if (c.size() != vertices_ * species_) throw std::invalid_argument("physis::chemistry::ReactionDiffusion::rate: one concentration for each vertex and species");
    std::vector<double> dc(c.size(), 0.0);
    for (std::size_t v = 0; v < vertices_; v++) {
      const double* x = c.data() + v * species_;
      double* dx = dc.data() + v * species_;
      for (const auto& r : reactions_) {  // mass action, as mass_action() of kinetics.hpp
        double velocity = r.k;
        for (const auto& [i, a] : r.reactants)
          for (unsigned p = 0; p < a; p++) velocity *= x[i];
        for (const auto& [i, a] : r.reactants) dx[i] -= a * velocity;
        for (const auto& [i, b] : r.products) dx[i] += b * velocity;
      }
      for (const auto& [w, L] : rows_[v])
        for (std::size_t s = 0; s < species_; s++) dx[s] -= diffusion_[s] * L * c[w * species_ + s];
    }
    return dc;
  }

  // One step of RK4.
  std::vector<double> step(std::span<const double> c, double h) const {
    return runge_kutta_4([this](std::span<const double> y) { return rate(y); }, c, h);
  }

  std::size_t vertices() const { return vertices_; }
  std::size_t species() const { return species_; }
  // RK4's step bound for the diffusion alone, 2.78 / (d_max · 2 max_v L_vv) (infinite without diffusion).
  double diffusive_step() const {
    const double d = diffusion_.empty() ? 0.0 : *std::max_element(diffusion_.begin(), diffusion_.end());
    return d > 0 && diagonal_ > 0 ? 2.78 / (d * 2 * diagonal_) : std::numeric_limits<double>::infinity();
  }

 private:
  struct Cached {
    double k = 0;
    std::vector<std::pair<std::size_t, unsigned>> reactants, products;
  };
  std::size_t species_ = 0, vertices_ = 0;
  std::vector<double> diffusion_;
  std::vector<Cached> reactions_;
  std::vector<std::vector<std::pair<std::size_t, double>>> rows_;
  double diagonal_ = 0;
};

}  // namespace physis::chemistry
