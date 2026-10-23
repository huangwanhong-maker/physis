// Checks of physis/symmetry's groups: SU(2) as quaternions against its matrices, the double
// cover SU(2) → SO(3), the exponential and logarithm, characters and their Clebsch–Gordan
// series, Haar moments by Monte Carlo; SU(N) generators, exponential, reunitarisation and Haar
// moments; the Lorentz group, its cover by SL(2, ℂ), rapidity addition and the Wigner rotation
// against its closed form; Weyl integration: orthonormal characters and the character
// expansion of Wilson's weight against modified Bessel functions.

#include <physis/symmetry/lorentz.hpp>
#include <physis/symmetry/su2.hpp>
#include <physis/symmetry/sun.hpp>
#include <physis/symmetry/u1.hpp>
#include <physis/symmetry/weyl.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <string>

namespace ps = physis::symmetry;
using ps::Complex;
using ps::Matrix2c;

namespace {

double diff(const ps::SU2& x, const ps::SU2& y) {
  return std::max({std::fabs(x.a0 - y.a0), std::fabs(x.a1 - y.a1), std::fabs(x.a2 - y.a2), std::fabs(x.a3 - y.a3)});
}
double diff(const Matrix2c& x, const Matrix2c& y) {
  double d = 0;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++) d = std::max(d, std::abs(x[i][j] - y[i][j]));
  return d;
}
double diff(const ps::Matrix4& x, const ps::Matrix4& y) {
  double d = 0;
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) d = std::max(d, std::fabs(x[i][j] - y[i][j]));
  return d;
}
double diff(const ps::SUN& x, const ps::SUN& y) {
  double d = 0;
  for (std::size_t i = 0; i < x.a.size(); i++) d = std::max(d, std::abs(x.a[i] - y.a[i]));
  return d;
}
ps::Vector3 unit(calx::random::Stream& z) {
  for (;;) {
    ps::Vector3 v{z.normal(), z.normal(), z.normal()};
    const double n = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (n > 1e-6) return {v[0] / n, v[1] / n, v[2] / n};
  }
}
// The modified Bessel function I_n(x) by its series, all of whose terms are positive.
double bessel_i(int n, double x) {
  double term = std::pow(x / 2, n) / std::tgamma(n + 1.0), sum = term;
  for (int k = 1; k < 200; k++) {
    term *= (x / 2) * (x / 2) / (k * static_cast<double>(k + n));
    sum += term;
    if (term < 1e-18 * sum) break;
  }
  return sum;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream z(20261023, 1, 0, 0);
  constexpr double pi = std::numbers::pi;

  // ---- SU(2) -----------------------------------------------------------------------------------------------------------
  {
    double prod = 0, assoc = 0, inv = 0, round = 0, cover = 0, cover_hom = 0;
    for (int k = 0; k < 500; k++) {
      const ps::SU2 x = ps::haar_su2(z), y = ps::haar_su2(z), w = ps::haar_su2(z);
      // The quaternion product against the product of the matrices.
      Matrix2c mx = x.matrix(), my = y.matrix(), mxy{};
      for (int i = 0; i < 2; i++)
        for (int kk = 0; kk < 2; kk++)
          for (int j = 0; j < 2; j++) mxy[i][j] += mx[i][kk] * my[kk][j];
      prod = std::max(prod, diff((x * y).matrix(), mxy));
      assoc = std::max(assoc, diff((x * y) * w, x * (y * w)));
      inv = std::max(inv, diff(x * x.inverse(), ps::SU2::identity()));
      round = std::max(round, diff(ps::from_matrix(x.matrix()), x));
      // R(U)_ij = ½ tr(σ_i U σ_j U†) against Rodrigues' form.
      const auto R = ps::adjoint(x);
      const auto& sg = ps::pauli();
      Matrix2c ud{};
      for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++) ud[i][j] = std::conj(mx[j][i]);
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
          Matrix2c a{}, b{}, c{};
          for (int p = 0; p < 2; p++)
            for (int q = 0; q < 2; q++)
              for (int r = 0; r < 2; r++) a[p][q] += sg[i][p][r] * mx[r][q];
          for (int p = 0; p < 2; p++)
            for (int q = 0; q < 2; q++)
              for (int r = 0; r < 2; r++) b[p][q] += a[p][r] * sg[j][r][q];
          for (int p = 0; p < 2; p++)
            for (int q = 0; q < 2; q++)
              for (int r = 0; r < 2; r++) c[p][q] += b[p][r] * ud[r][q];
          cover = std::max(cover, std::fabs(0.5 * (c[0][0] + c[1][1]).real() - R[i][j]));
        }
      // R(xy) = R(x) R(y), and R(−x) = R(x).
      const auto Rx = ps::adjoint(x), Ry = ps::adjoint(y), Rxy = ps::adjoint(x * y), Rm = ps::adjoint(-1.0 * x);
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
          double v = 0;
          for (int q = 0; q < 3; q++) v += Rx[i][q] * Ry[q][j];
          cover_hom = std::max({cover_hom, std::fabs(v - Rxy[i][j]), std::fabs(Rm[i][j] - Rx[i][j])});
        }
    }
    s.check("symmetry: SU(2) quaternion product against the product of its 2x2 matrices, 500 Haar pairs", prod, 4e-15);
    s.check("symmetry: SU(2) associativity and inverses", std::max(assoc, inv), 4e-15);
    s.check("symmetry: SU(2) element read back from its matrix", round, 1e-16);
    s.check("symmetry: the adjoint map R_ij = tr(s_i U s_j U+)/2 against Rodrigues' form", cover, 4e-15);
    s.check("symmetry: SU(2) -> SO(3) a homomorphism with R(-U) = R(U)", cover_hom, 4e-15);
    // exp(θ e_z · T) is the right-handed rotation by θ about z.
    const auto Rz = ps::adjoint(ps::rotation_su2({0, 0, 1}, 0.7));
    s.check("symmetry: exp(0.7 e_z . T) covers the right-handed rotation by 0.7 about z",
            std::max({std::fabs(Rz[0][0] - std::cos(0.7)), std::fabs(Rz[1][0] - std::sin(0.7)), std::fabs(Rz[0][1] + std::sin(0.7)),
                      std::fabs(Rz[2][2] - 1)}),
            1e-15);
    double el = 0;
    for (int k = 0; k < 500; k++) {
      const ps::Vector3 n = unit(z);
      const double t = (2 * pi - 1e-3) * z.uniform();
      const ps::Vector3 x{t * n[0], t * n[1], t * n[2]};
      const ps::Vector3 y = ps::log_su2(ps::exp_su2(x));
      el = std::max({el, std::fabs(y[0] - x[0]), std::fabs(y[1] - x[1]), std::fabs(y[2] - x[2])});
    }
    s.check("symmetry: log(exp x) = x for |x| < 2 pi, 500 random x", el, 1e-11);

    // Characters: the closed form, the values at ±I, and the Clebsch–Gordan series χ_½ χ_j = χ_{j−½} + χ_{j+½}.
    double ch = 0, ends = 0, series = 0;
    for (int two_j = 0; two_j <= 12; two_j++) {
      for (int k = 1; k < 50; k++) {
        const double t = 2 * pi * k / 50.0 + 0.01;
        ch = std::max(ch, std::fabs(ps::character_su2(two_j, t) - std::sin((two_j + 1) * t / 2) / std::sin(t / 2)));
        if (two_j >= 1)
          series = std::max(series, std::fabs(ps::character_su2(1, t) * ps::character_su2(two_j, t) - ps::character_su2(two_j - 1, t) -
                                              ps::character_su2(two_j + 1, t)));
      }
      ends = std::max({ends, std::fabs(ps::character(two_j, ps::SU2::identity()) - (two_j + 1)),
                       std::fabs(ps::character(two_j, ps::SU2{-1, 0, 0, 0}) - (two_j % 2 ? -1 : 1) * (two_j + 1))});
    }
    s.check("symmetry: SU(2) characters by Chebyshev's recurrence against sin((2j+1)t/2)/sin(t/2), 2j <= 12", ch, 1e-12);
    s.check("symmetry: chi_j(I) = 2j + 1 and chi_j(-I) = (-1)^2j (2j + 1)", ends, 0);
    s.check("symmetry: the Clebsch-Gordan series chi_1/2 chi_j = chi_j-1/2 + chi_j+1/2", series, 1e-12);

    // Haar moments: E[a0²] = 1/4 (variance 1/16), E[χ_½] = E[χ_1] = 0 (variance 1); limits are five standard errors.
    const int M = 200000;
    double a2 = 0, c1 = 0, c2 = 0;
    for (int k = 0; k < M; k++) {
      const ps::SU2 u = ps::haar_su2(z);
      a2 += u.a0 * u.a0;
      c1 += ps::character(1, u);
      c2 += ps::character(2, u);
    }
    s.check("symmetry: SU(2) Haar measure, E[a0^2] - 1/4 (200000 draws; limit five standard errors)", a2 / M - 0.25, 5 * 0.25 / std::sqrt(M));
    s.check("symmetry: SU(2) Haar measure, E[chi_1/2] and E[chi_1] (limit five standard errors)", std::max(std::fabs(c1 / M), std::fabs(c2 / M)),
            5 / std::sqrt(M));
  }
  // ---- U(1) ------------------------------------------------------------------------------------------------------------
  {
    const ps::U1 a = ps::U1::from_angle(3.0), b = ps::U1::from_angle(1.0);
    s.check("symmetry: U(1) product wraps into (-pi, pi]", (a * b).theta - (4.0 - 2 * pi), 1e-15);
    const int M = 100000;
    double c = 0, c2 = 0;
    for (int k = 0; k < M; k++) {
      const ps::U1 u = ps::haar_u1(z);
      c += u.re_trace();
      c2 += u.re_trace() * u.re_trace();
    }
    s.check("symmetry: U(1) Haar measure, E[cos] = 0 and E[cos^2] = 1/2 (limit five standard errors)",
            std::max(std::fabs(c / M), std::fabs(c2 / M - 0.5)), 5 / std::sqrt(2.0 * M));
  }
  // ---- SU(N) -----------------------------------------------------------------------------------------------------------
  {
    double gen = 0, expd = 0, unit_def = 0;
    bool counts = true;
    for (std::size_t n = 2; n <= 5; n++) {
      const auto t = ps::generators(n);
      counts = counts && t.size() == n * n - 1;
      for (std::size_t a = 0; a < t.size(); a++) {
        gen = std::max({gen, diff(t[a], t[a].adjoint()), std::abs(t[a].trace())});
        for (std::size_t b = 0; b < t.size(); b++) gen = std::max(gen, std::abs((t[a] * t[b]).trace() - Complex(a == b ? 0.5 : 0.0)));
      }
      std::vector<double> x(n * n - 1);
      for (auto& v : x) v = z.normal();
      const ps::SUN U = ps::exp_sun(n, x);
      unit_def = std::max({unit_def, U.unitarity_defect(), std::abs(U.determinant() - Complex(1))});
    }
    s.require("symmetry: su(N) has N^2 - 1 generalised Gell-Mann generators, N = 2..5", counts);
    s.check("symmetry: the generators are Hermitian and traceless with tr(t_a t_b) = delta_ab / 2", gen, 1e-15);
    s.check("symmetry: exp(-i x^a t_a) is unitary with determinant 1, N = 2..5", unit_def, 1e-13);
    for (int k = 0; k < 50; k++) {
      const ps::Vector3 x{z.normal(), z.normal(), z.normal()};
      const ps::SUN U = ps::exp_sun(2, {x[0], x[1], x[2]});
      const Matrix2c m = ps::exp_su2(x).matrix();
      for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++) expd = std::max(expd, std::abs(U(i, j) - m[i][j]));
    }
    s.check("symmetry: exp_sun for N = 2 agrees with exp_su2", expd, 1e-13);
    double reu = 0, keep = 0;
    for (int k = 0; k < 100; k++) {
      const ps::SUN U = ps::haar_sun(3, z);
      ps::SUN V = U;
      for (auto& v : V.a) v += 1e-6 * Complex(z.normal(), z.normal());
      const ps::SUN W = V.reunitarised();
      reu = std::max({reu, W.unitarity_defect(), std::abs(W.determinant() - Complex(1))});
      keep = std::max(keep, diff(U.reunitarised(), U));
    }
    s.check("symmetry: a perturbed SU(3) element reunitarised is unitary with determinant 1", reu, 1e-14);
    s.check("symmetry: reunitarisation leaves an element of SU(3) as it is", keep, 1e-14);
    // Haar moments of SU(3): E|U_11|² = 1/3 (variance 2/12 − 1/9), E|tr U|² = 1 (variance 1), E tr U = 0.
    const int M = 40000;
    double u11 = 0, tr2 = 0;
    Complex tr = 0;
    for (int k = 0; k < M; k++) {
      const ps::SUN U = ps::haar_sun(3, z);
      u11 += std::norm(U(0, 0));
      tr2 += std::norm(U.trace());
      tr += U.trace();
    }
    s.check("symmetry: SU(3) Haar measure, E|U_11|^2 - 1/3 (40000 draws; limit five standard errors)", u11 / M - 1.0 / 3,
            5 * std::sqrt(2.0 / 12 - 1.0 / 9) / std::sqrt(M));
    s.check("symmetry: SU(3) Haar measure, E|tr U|^2 - 1 and |E tr U| (limit five standard errors)",
            std::max(std::fabs(tr2 / M - 1), std::abs(tr) / M), 5 / std::sqrt(M));
  }
  // ---- the Lorentz group -----------------------------------------------------------------------------------------------
  {
    double defect = 0, add = 0, inv = 0, cov = 0, hom = 0, kernel = 0, orth = 0, wig = 0, pure = 0;
    const ps::Matrix4 I = ps::multiply(ps::boost(0, {1, 0, 0}), ps::rotation({0, 0, 1}, 0));
    for (int k = 0; k < 300; k++) {
      const ps::Vector3 n1 = unit(z), n2 = unit(z), n3 = unit(z);
      const double f1 = 3 * z.uniform(), f2 = 3 * z.uniform(), th = 2 * pi * z.uniform();
      const ps::Matrix4 L = ps::multiply(ps::multiply(ps::boost(f1, n1), ps::rotation(n3, th)), ps::boost(f2, n2));
      defect = std::max(defect, ps::lorentz_defect(L) / (L[0][0] * L[0][0]));
      add = std::max(add, diff(ps::multiply(ps::boost(f1, n1), ps::boost(f2, n1)), ps::boost(f1 + f2, n1)) / std::cosh(f1 + f2));
      inv = std::max(inv, diff(ps::multiply(ps::lorentz_inverse(L), L), I) / (L[0][0] * L[0][0]));
      // The cover: boosts and rotations, a homomorphism on random elements of SL(2, ℂ), with kernel ±I.
      cov = std::max({cov, diff(ps::cover(ps::boost_sl2c(f1, n1)), ps::boost(f1, n1)) / std::cosh(f1),
                      diff(ps::cover(ps::rotation_su2(n3, th).matrix()), ps::rotation(n3, th))});
      Matrix2c A, B;
      for (auto* M : {&A, &B}) {
        for (auto& row : *M)
          for (auto& v : row) v = Complex(z.normal(), z.normal());
        const Complex d = std::sqrt((*M)[0][0] * (*M)[1][1] - (*M)[0][1] * (*M)[1][0]);
        for (auto& row : *M)
          for (auto& v : row) v /= d;
      }
      const ps::Matrix4 LA = ps::cover(A), LB = ps::cover(B), LAB = ps::cover(ps::multiply(A, B));
      double scale = 0;
      for (auto& row : LAB)
        for (double v : row) scale = std::max(scale, std::fabs(v));
      hom = std::max(hom, diff(ps::multiply(LA, LB), LAB) / scale);
      Matrix2c mA = A;
      for (auto& row : mA)
        for (auto& v : row) v = -v;
      kernel = std::max(kernel, diff(ps::cover(mA), LA));
      orth = std::max(orth, ps::lorentz_defect(LA) / (LA[0][0] * LA[0][0]) + (LA[0][0] >= 1 ? 0 : 1));
      // The Wigner rotation of two boosts against its closed form.
      const ps::Matrix4 C = ps::multiply(ps::boost(f2, n2), ps::boost(f1, n1));
      const auto d = ps::decompose(C);
      wig = std::max(wig, std::fabs(ps::rotation_angle(d.rotation) - ps::wigner_angle(f1, n1, f2, n2)));
      pure = std::max({pure, std::fabs(d.rotation[0][0] - 1), std::fabs(d.rotation[0][1]), std::fabs(d.rotation[1][0])});
    }
    s.check("symmetry: products of boosts and rotations preserve eta (relative to (L00)^2)", defect, 1e-14);
    s.check("symmetry: collinear rapidities add", add, 1e-14);
    s.check("symmetry: the inverse eta L^T eta", inv, 1e-14);
    s.check("symmetry: SL(2,C) -> SO+(1,3) covers boosts and rotations", cov, 1e-14);
    s.check("symmetry: the cover is a homomorphism on random elements of SL(2,C) (relative)", hom, 1e-13);
    s.check("symmetry: A and -A cover the same Lorentz transformation", kernel, 0);
    s.check("symmetry: the cover lands in SO+(1,3): eta preserved and L00 >= 1", orth, 1e-13);
    s.check("symmetry: the Wigner rotation of two boosts against its closed form (rad)", wig, 1e-10);
    s.check("symmetry: the Wigner rotation is a pure rotation", pure, 1e-11);
  }
  // ---- Weyl integration ------------------------------------------------------------------------------------------------
  {
    double orth = 0;
    for (int m = -3; m <= 3; m++)
      for (int n = -3; n <= 3; n++) {
        const auto c = ps::coefficient_u1(n, [&](double t) { return std::cos(m * t); });
        const auto c2 = ps::coefficient_u1(n, [&](double t) { return std::sin(m * t); });
        const Complex chi_m_coef = c + Complex(0, 1) * c2;  // ∫ e^{imθ} e^{−inθ}
        orth = std::max(orth, std::abs(chi_m_coef - Complex(m == n ? 1.0 : 0.0)));
      }
    s.check("symmetry: U(1) characters orthonormal under the Weyl integral, |m|, |n| <= 3", orth, 1e-13);
    double bes = 0;
    for (double beta : {0.5, 2.0, 5.0})
      for (int n = 0; n <= 6; n++) {
        const auto c = ps::coefficient_u1(n, [&](double t) { return std::exp(beta * std::cos(t)); });
        bes = std::max({bes, std::fabs(c.real() - bessel_i(n, beta)) / bessel_i(0, beta), std::fabs(c.imag())});
      }
    s.check("symmetry: U(1) character coefficients of exp(beta cos t) against I_n(beta), beta = 0.5, 2, 5 (relative to I_0)", bes, 1e-13);
    double orth2 = 0;
    for (int a = 0; a <= 6; a++)
      for (int b = 0; b <= 6; b++)
        orth2 = std::max(orth2, std::fabs(ps::coefficient_su2(b, [&](double t) { return ps::character_su2(a, t); }) - (a == b ? 1.0 : 0.0)));
    s.check("symmetry: SU(2) characters orthonormal under the Weyl integral, 2j <= 6", orth2, 1e-13);
    double bes2 = 0;
    for (double beta : {0.5, 2.0, 5.0}) {
      const double c0 = 2 * bessel_i(1, beta) / beta;
      for (int two_j = 0; two_j <= 6; two_j++) {
        const double c = ps::coefficient_su2(two_j, [&](double t) { return std::exp(beta * std::cos(t / 2)); });
        bes2 = std::max(bes2, std::fabs(c - 2 * (two_j + 1) * bessel_i(two_j + 1, beta) / beta) / c0);
      }
    }
    s.check("symmetry: SU(2) coefficients of exp((beta/2) tr U) against 2(2j+1) I_2j+1(beta) / beta (relative to c_0)", bes2, 1e-13);
    const double one = ps::integrate_su3([](double, double) { return 1.0; });
    auto tr2 = [](double a, double b) { return std::norm(std::polar(1.0, a) + std::polar(1.0, b) + std::polar(1.0, -a - b)); };
    const double m2 = ps::integrate_su3(tr2), m4 = ps::integrate_su3([&](double a, double b) { return tr2(a, b) * tr2(a, b); });
    s.check("symmetry: SU(3) Weyl integral: the volume 1, E|tr U|^2 = 1 and E|tr U|^4 = 2",
            std::max({std::fabs(one - 1), std::fabs(m2 - 1), std::fabs(m4 - 2)}), 1e-10);
  }
  return s.done();
}
