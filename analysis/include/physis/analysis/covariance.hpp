// The covariance of a random potential field: the prior law of a field known in part, at any two points
// (decisions 0010 and 0013).
//
// An isotropic Gaussian potential field (random_field.hpp) has independent coefficients gₙᵐ, hₙᵐ of mean 0
// and variance σₙ² = E[Rₙ]/((n + 1)(2n + 1)). Its potential V = a Σₙ (a/r)ⁿ⁺¹ Σₘ (gₙᵐ cos mλ + hₙᵐ sin mλ) P̆ₙᵐ
// then has the covariance
//
//   K(x, y) = Cov(V(x), V(y)) = a² Σₙ σₙ² ρⁿ⁺¹ Pₙ(t),   ρ = a²/(|x||y|),   t = x̂ · ŷ,
//
// by the addition theorem of Schmidt's functions, Σₘ P̆ₙᵐ(sin ψ) P̆ₙᵐ(sin ψ′) cos m(λ − λ′) = Pₙ(cos γ)
// (Langel 1987, eq. 4.14; Winch et al. 2005, Geophys. J. Int. 160, 487–504, §2). It depends on the
// points only through their radii and the angle between them: the law is isotropic.
//
// The field B = −∇V has the covariance Cov(Bᵢ(x), Bⱼ(y)) = ∂²K/∂xᵢ∂yⱼ (the two signs cancel). With
// s = |x|, s′ = |y|, u = x̂, v = ŷ and K a function of (s, s′, t), the chain rule gives
//
//   ∂²K/∂xᵢ∂yⱼ = K_ss′ uᵢ vⱼ + K_st uᵢ βⱼ + K_s′t αᵢ vⱼ + K_tt αᵢ βⱼ + K_t [(δᵢⱼ − vᵢvⱼ) − uᵢ(uⱼ − t vⱼ)]/(s s′),
//
// αᵢ = ∂t/∂xᵢ = (vᵢ − t uᵢ)/s, βⱼ = ∂t/∂yⱼ = (uⱼ − t vⱼ)/s′, and the partial derivatives of K are six sums
// over the degrees,
//
//   K_t = B₀,  K_tt = C₀,  K_ss′ = A₂/(s s′),  K_st = −B₁/s,  K_s′t = −B₁/s′,
//   A_k = a² Σ σₙ² (n + 1)ᵏ ρⁿ⁺¹ Pₙ(t),  B_k = the same with Pₙ′(t),  C₀ = the same with Pₙ″(t), k = 0,
//
// computed in one pass by the three-term recurrence of Pₙ and the summed ones of its derivatives,
// P′ₙ₊₁ = P′ₙ₋₁ + (2n + 1) Pₙ and P″ₙ₊₁ = P″ₙ₋₁ + (2n + 1) P′ₙ, which hold at t = ±1 (Abramowitz and
// Stegun 1964, 8.5.3 and 22.7.10): the variance at a point (x = y, α = β = 0) needs no limit. O(N) per pair
// of points; the checks hold it to the sum over the coefficients, Σₖ σₖ² bₖ(x) bₖ(y)ᵀ, to rounding.
//
// A table. Points at one radius (a magnetometer's readings and a map at its height) need the sums only as
// functions of the angle θ between them. TabulatedCovariance computes them once at θ_k = k δ and
// interpolates by Lagrange's cubic through the four nearest, the sums being even in θ (the table is
// reflected below zero): its error is of order (N δ)⁴ of the sums, N the highest degree. Pairs at another
// radius, or farther apart than the table, are summed exactly.
//
// Cartesian components in the frame of the expansion (the body-fixed frame); B in the coefficients' unit.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace physis::analysis {

class PotentialCovariance {
 public:
  using Point = std::array<double, 3>;
  using Block = std::array<std::array<double, 3>, 3>;

  // From the expected spectrum E[Rₙ] at the reference sphere `radius`, for n_from ≤ n ≤ n_to.
  PotentialCovariance(const std::vector<double>& spectrum, int n_from, int n_to, double radius) : n_from_(n_from), n_to_(n_to), a_(radius) {
    if (n_from < 1 || n_to < n_from || n_to >= int(spectrum.size())) throw std::invalid_argument("physis::analysis::PotentialCovariance: degrees outside the spectrum");
    if (!(radius > 0)) throw std::invalid_argument("physis::analysis::PotentialCovariance: a positive radius");
    c_.assign(std::size_t(n_to) + 1, 0.0);
    for (int n = n_from; n <= n_to; n++) c_[n] = radius * radius * std::max(0.0, spectrum[n]) / ((n + 1.0) * (2.0 * n + 1.0));
  }

  int n_from() const { return n_from_; }
  int n_to() const { return n_to_; }
  double radius() const { return a_; }

  // The five sums over the degrees at radii s, s′ and t = cos γ.
  struct Sums {
    double A0 = 0, A2 = 0, B0 = 0, B1 = 0, C0 = 0;
  };

  // Cov(V(x), V(y)).
  double potential(const Point& x, const Point& y) const { return sums(x, y).A0; }

  // Cov(Bᵢ(x), Bⱼ(y)), B = −∇V, Cartesian.
  Block field(const Point& x, const Point& y) const {
    const Pair P = pair(x, y);
    return assemble(sums_at(P.s, P.sp, P.t), P);
  }

  // The geometry of a pair of points: radii, directions, t and the angle between them.
  struct Pair {
    double s = 0, sp = 0, t = 1, angle = 0;
    Point u{}, v{};
  };
  static Pair pair(const Point& x, const Point& y) {
    Pair P;
    P.s = norm(x), P.sp = norm(y);
    for (int i = 0; i < 3; i++) P.u[i] = x[i] / P.s, P.v[i] = y[i] / P.sp;
    P.t = std::clamp(P.u[0] * P.v[0] + P.u[1] * P.v[1] + P.u[2] * P.v[2], -1.0, 1.0);
    const Point c{P.u[1] * P.v[2] - P.u[2] * P.v[1], P.u[2] * P.v[0] - P.u[0] * P.v[2], P.u[0] * P.v[1] - P.u[1] * P.v[0]};
    P.angle = std::atan2(std::hypot(c[0], c[1], c[2]), P.t);  // keeps its digits when small
    return P;
  }

  // Cov(B(x), B(y)) from the sums, by the chain rule.
  static Block assemble(const Sums& S, const Pair& P) {
    const double s = P.s, sp = P.sp, t = P.t;
    const Point& u = P.u;
    const Point& v = P.v;
    const double Kt = S.B0, Ktt = S.C0, Kss = S.A2 / (s * sp), Kst = -S.B1 / s, Kspt = -S.B1 / sp;
    Block H{};
    for (int i = 0; i < 3; i++) {
      const double alpha = (v[i] - t * u[i]) / s;
      for (int j = 0; j < 3; j++) {
        const double beta = (u[j] - t * v[j]) / sp;
        const double ddt = ((i == j ? 1.0 : 0.0) - v[i] * v[j] - u[i] * (u[j] - t * v[j])) / (s * sp);
        H[i][j] = Kss * u[i] * v[j] + Kst * u[i] * beta + Kspt * alpha * v[j] + Ktt * alpha * beta + Kt * ddt;
      }
    }
    return H;
  }

  // E|B(x)|², the trace of field(x, x): at the reference sphere, Σₙ E[Rₙ].
  double total_variance(const Point& x) const {
    const Block H = field(x, x);
    return H[0][0] + H[1][1] + H[2][2];
  }

  // The sums at radii s, s′ and t = cos γ: O(N).
  Sums sums_at(double s, double sp, double t) const {
    const double rho = a_ * a_ / (s * sp);
    Sums S;
    // P, P′, P″ at n − 1 and n, from n = 0; ρⁿ⁺¹ alongside.
    double p0 = 1, p1 = t, d0 = 0, d1 = 1, e0 = 0, e1 = 0;
    double power = rho * rho;  // ρⁿ⁺¹ at n = 1
    for (int n = 1; n <= n_to_; n++) {
      if (n >= n_from_) {
        const double w = c_[n] * power, k = n + 1.0;
        S.A0 += w * p1;
        S.A2 += w * k * k * p1;
        S.B0 += w * d1;
        S.B1 += w * k * d1;
        S.C0 += w * e1;
      }
      const double p2 = ((2.0 * n + 1) * t * p1 - n * p0) / (n + 1.0);
      const double d2 = d0 + (2.0 * n + 1) * p1;
      const double e2 = e0 + (2.0 * n + 1) * d1;
      p0 = p1, p1 = p2, d0 = d1, d1 = d2, e0 = e1, e1 = e2;
      power *= rho;
    }
    return S;
  }

 private:
  static double norm(const Point& x) {
    const double r = std::hypot(x[0], x[1], x[2]);
    if (!(r > 0)) throw std::domain_error("physis::analysis::PotentialCovariance: a point at the centre");
    return r;
  }
  Sums sums(const Point& x, const Point& y) const {
    const double s = norm(x), sp = norm(y);
    const double t = std::clamp((x[0] * y[0] + x[1] * y[1] + x[2] * y[2]) / (s * sp), -1.0, 1.0);
    return sums_at(s, sp, t);
  }

  int n_from_, n_to_;
  double a_;
  std::vector<double> c_;  // a² σₙ²
};

// The covariance of points at one radius, tabulated in the angle between them.
class TabulatedCovariance {
 public:
  using Point = PotentialCovariance::Point;
  using Block = PotentialCovariance::Block;

  // The sums at `points` + 1 angles from 0 to max_angle, at the radius r.
  TabulatedCovariance(PotentialCovariance exact, double radius, double max_angle, int points)
      : exact_(std::move(exact)), r_(radius), step_(max_angle / points) {
    if (!(radius > 0 && max_angle > 0 && points >= 4)) throw std::invalid_argument("physis::analysis::TabulatedCovariance: a positive radius and angle, four points or more");
    table_.reserve(std::size_t(points) + 1);
    for (int k = 0; k <= points; k++) table_.push_back(exact_.sums_at(radius, radius, std::cos(k * step_)));
  }

  const PotentialCovariance& exact() const { return exact_; }
  double radius() const { return r_; }
  double max_angle() const { return step_ * double(table_.size() - 1); }

  // Cov(B(x), B(y)): from the table when both points are at its radius (to 10⁻⁹ of it) and the angle between
  // them within it, summed exactly otherwise.
  Block field(const Point& x, const Point& y) const {
    const auto P = PotentialCovariance::pair(x, y);
    const double tol = 1e-9 * r_;
    const double k = P.angle / step_;
    if (std::fabs(P.s - r_) > tol || std::fabs(P.sp - r_) > tol || k > double(table_.size()) - 3) return exact_.field(x, y);
    const int i = int(std::floor(k));
    const double f = k - i;
    // Lagrange's cubic through i − 1, i, i + 1, i + 2, reflected below zero (the sums are even in θ).
    const double w[4] = {-f * (f - 1) * (f - 2) / 6, (f + 1) * (f - 1) * (f - 2) / 2, -(f + 1) * f * (f - 2) / 2, (f + 1) * f * (f - 1) / 6};
    PotentialCovariance::Sums S;
    for (int j = 0; j < 4; j++) {
      const auto& T = table_[std::size_t(std::abs(i - 1 + j))];
      S.A0 += w[j] * T.A0, S.A2 += w[j] * T.A2, S.B0 += w[j] * T.B0, S.B1 += w[j] * T.B1, S.C0 += w[j] * T.C0;
    }
    return PotentialCovariance::assemble(S, P);
  }

 private:
  PotentialCovariance exact_;
  double r_, step_;
  std::vector<PotentialCovariance::Sums> table_;
};

}  // namespace physis::analysis
