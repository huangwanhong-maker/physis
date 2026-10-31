// Checks of the exact memory left by eliminating variables from a linear system (elimination.hpp).
//
// - On a random stable system of 2 kept and 3 eliminated variables, solved exactly (the matrix exponential
//   of the whole), the reduced equation ẋ = A x + B e^{Dt} y(0) + ∫₀ᵗ K(t − s) x(s) ds reproduces ẋ at
//   several times, the convolution by the trapezoid rule on a fine grid (error O(h²)).
// - K(0) = B C, and ∫₀^∞ K = −B D⁻¹ C, which makes A − B D⁻¹ C the Markovian limit.
// - Adiabatic elimination: with D made faster by 1/ε, the Markovian equation approaches the kept variables'
//   exact path, the error falling in proportion to ε.

#include <cmath>
#include <vector>

#include <calx/linalg/expm.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <physis/dynamics/elimination.hpp>

using calx::linalg::Matrix;
using calx::linalg::Rectangular;
using calx::verification::Suite;

namespace {

std::vector<double> mul(const Matrix& M, const std::vector<double>& v) {
  std::vector<double> out(M.n, 0.0);
  for (std::size_t i = 0; i < M.n; i++)
    for (std::size_t j = 0; j < M.n; j++) out[i] += M(i, j) * v[j];
  return out;
}

// D divided by `scale` (faster when scale < 1), B and C by its square root, so that B D⁻¹ C is unchanged.
physis::dynamics::LinearElimination system(calx::random::Stream& rng, double scale = 1) {
  physis::dynamics::LinearElimination L;
  L.A = Matrix(2), L.D = Matrix(3), L.B = Rectangular(2, 3), L.C = Rectangular(3, 2);
  for (std::size_t i = 0; i < 2; i++)
    for (std::size_t j = 0; j < 2; j++) L.A(i, j) = 0.3 * rng.normal() - (i == j ? 0.6 : 0);
  for (std::size_t i = 0; i < 3; i++)
    for (std::size_t j = 0; j < 3; j++) L.D(i, j) = (0.3 * rng.normal() - (i == j ? 1.5 : 0)) / scale;
  for (std::size_t i = 0; i < 2; i++)
    for (std::size_t a = 0; a < 3; a++) L.B(i, a) = 0.5 * rng.normal() / std::sqrt(scale), L.C(a, i) = 0.5 * rng.normal() / std::sqrt(scale);
  return L;
}

}  // namespace

int main() {
  Suite s;
  calx::random::Stream rng(41, 0, 0, 2);
  const auto L = system(rng);
  const Matrix M = L.full();
  const std::vector<double> z0{1.0, -0.5, 0.7, 0.2, -1.1};
  const std::vector<double> y0{0.7, 0.2, -1.1};
  const double h = 1e-3;
  const std::size_t steps = 2000;
  // The exact path on the grid.
  std::vector<std::vector<double>> z(steps + 1);
  const Matrix P = calx::linalg::expm(h * M);
  z[0] = z0;
  for (std::size_t k = 0; k < steps; k++) z[k + 1] = mul(P, z[k]);
  std::vector<Matrix> K(steps + 1);
  for (std::size_t k = 0; k <= steps; k++) K[k] = L.kernel(k * h);
  double worst = 0;
  for (std::size_t n : {500u, 1000u, 2000u}) {
    const double t = n * h;
    std::vector<double> xdot(2, 0.0), reduced(2, 0.0);
    for (std::size_t i = 0; i < 2; i++) {
      for (std::size_t j = 0; j < 2; j++) xdot[i] += L.A(i, j) * z[n][j];
      for (std::size_t a = 0; a < 3; a++) xdot[i] += L.B(i, a) * z[n][2 + a];
    }
    const auto noise = L.noise(t, y0);
    for (std::size_t i = 0; i < 2; i++) {
      reduced[i] = noise[i];
      for (std::size_t j = 0; j < 2; j++) reduced[i] += L.A(i, j) * z[n][j];
      for (std::size_t k = 0; k <= n; k++) {
        const double w = (k == 0 || k == n) ? 0.5 * h : h;
        for (std::size_t j = 0; j < 2; j++) reduced[i] += w * K[n - k](i, j) * z[k][j];
      }
    }
    for (std::size_t i = 0; i < 2; i++) worst = std::max(worst, std::fabs(reduced[i] - xdot[i]));
  }
  s.check("the reduced equation with its memory and noise reproduces ẋ", worst, 1e-6);

  double bc = 0;
  for (std::size_t i = 0; i < 2; i++)
    for (std::size_t j = 0; j < 2; j++) {
      double v = 0;
      for (std::size_t a = 0; a < 3; a++) v += L.B(i, a) * L.C(a, j);
      bc = std::max(bc, std::fabs(K[0](i, j) - v));
    }
  s.check("K(0) = B C", bc, 1e-14);
  {
    // ∫₀^∞ K by Simpson's rule to τ = 40, against A − (A − B D⁻¹ C).
    const auto mk = L.markovian();
    Matrix integral(2);
    const double dt = 0.01;
    for (std::size_t k = 0; k <= 4000; k++) {
      const double w = (k == 0 || k == 4000) ? 1 : (k % 2 ? 4 : 2);
      const auto Kt = L.kernel(k * dt);
      for (std::size_t i = 0; i < 2; i++)
        for (std::size_t j = 0; j < 2; j++) integral(i, j) += w * dt / 3 * Kt(i, j);
    }
    double e = 0;
    for (std::size_t i = 0; i < 2; i++)
      for (std::size_t j = 0; j < 2; j++) e = std::max(e, std::fabs(integral(i, j) - (L.A(i, j) - mk(i, j)) * -1));
    s.check("∫₀^∞ K = −B D⁻¹ C", e, 1e-9);
  }
  {
    // Adiabatic elimination: the error of the Markovian equation at t = 3, for D faster by 1/ε, y(0) at its
    // quasi-steady value −D⁻¹ C x(0) so that no initial layer remains.
    std::vector<double> err;
    for (double eps : {0.04, 0.02, 0.01}) {
      calx::random::Stream r2(41, 0, 0, 3);
      const auto Lf = system(r2, eps);
      const std::vector<double> x0{1.0, -0.5};
      const auto mk = Lf.markovian();
      // y(0) = −D⁻¹ C x0.
      Matrix rhs(3);
      for (std::size_t a = 0; a < 3; a++) rhs(a, 0) = -(Lf.C(a, 0) * x0[0] + Lf.C(a, 1) * x0[1]);
      const auto yq = calx::linalg::solve(Lf.D, rhs);
      const std::vector<double> zf{x0[0], x0[1], yq(0, 0), yq(1, 0), yq(2, 0)};
      const auto exact = mul(calx::linalg::expm(3.0 * Lf.full()), zf);
      const auto approx = mul(calx::linalg::expm(3.0 * mk), x0);
      err.push_back(std::hypot(exact[0] - approx[0], exact[1] - approx[1]));
    }
    s.check("adiabatic elimination: halving ε halves the error (ratio of errors − 2)", err[0] / err[1] - 2, 0.3);
    s.check("and again", err[1] / err[2] - 2, 0.3);
  }
  return s.done();
}
