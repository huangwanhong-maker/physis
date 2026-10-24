// Checks of physis/symmetry's recoupling theory, exactly where the arithmetic allows: known
// values of Clebsch–Gordan coefficients and 6j symbols; the orthogonality of Clebsch–Gordan
// coefficients and of 6j symbols, the symmetries of 3j, 6j and 9j symbols, the Racah sum rule
// and the Biedenharn–Elliott identity, all as exact sums of surds; the 9j symbol against its
// defining sum in floating point and its reduction when a spin is zero; the Ponzano–Regge
// asymptotic formula at large spins; and the quantum 6j symbols at roots of unity: total
// quantum dimension, orthogonality, Biedenharn–Elliott, tetrahedral symmetry and the
// classical limit.

#include <physis/symmetry/quantum.hpp>
#include <physis/symmetry/recoupling.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <string>
#include <vector>

namespace ps = physis::symmetry;
using ps::Rational;
using ps::SignedSqrt;

namespace {

int phase(int two_x) { return (two_x / 2) % 2 == 0 ? 1 : -1; }  // (−1)^{x} for an integer x given as two_x

// The Ponzano–Regge value cos(Σ (j_e + ½) θ_e + π/4) / √(12πV) of {j1 j2 j3; j4 j5 j6} (two_j),
// with the tetrahedron whose edges are 01 = j1, 02 = j2, 12 = j3, 23 = j4, 13 = j5, 03 = j6 and
// θ_e its exterior dihedral angles; also the amplitude 1/√(12πV).
std::array<double, 2> ponzano_regge(const std::array<int, 6>& two_j) {
  std::array<double, 6> l;
  for (int e = 0; e < 6; e++) l[e] = two_j[e] / 2.0 + 0.5;
  const double d01 = l[0], d02 = l[1], d12 = l[2], d23 = l[3], d13 = l[4], d03 = l[5];
  using V = std::array<double, 3>;
  std::array<V, 4> P{};
  P[1] = {d01, 0, 0};
  const double x2 = (d01 * d01 + d02 * d02 - d12 * d12) / (2 * d01), y2 = std::sqrt(d02 * d02 - x2 * x2);
  P[2] = {x2, y2, 0};
  const double x3 = (d03 * d03 - d13 * d13 + d01 * d01) / (2 * d01);
  const double y3 = ((d03 * d03 + d02 * d02 - d23 * d23) / 2 - x3 * x2) / y2;
  const double z3 = std::sqrt(d03 * d03 - x3 * x3 - y3 * y3);
  P[3] = {x3, y3, z3};
  const double vol = d01 * y2 / 2 * z3 / 3;
  auto sub = [](const V& a, const V& b) { return V{a[0] - b[0], a[1] - b[1], a[2] - b[2]}; };
  auto dot = [](const V& a, const V& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; };
  auto exterior = [&](int i, int j, int k, int m) {
    V e = sub(P[j], P[i]);
    const double n = std::sqrt(dot(e, e));
    for (auto& c : e) c /= n;
    V u = sub(P[k], P[i]), w = sub(P[m], P[i]);
    const double pu = dot(u, e), pw = dot(w, e);
    for (int c = 0; c < 3; c++) { u[c] -= pu * e[c]; w[c] -= pw * e[c]; }
    return std::numbers::pi - std::acos(dot(u, w) / std::sqrt(dot(u, u) * dot(w, w)));
  };
  const std::array<double, 6> th{exterior(0, 1, 2, 3), exterior(0, 2, 1, 3), exterior(1, 2, 0, 3),
                                 exterior(2, 3, 0, 1), exterior(1, 3, 0, 2), exterior(0, 3, 1, 2)};
  double phi = std::numbers::pi / 4;
  for (int e = 0; e < 6; e++) phi += l[e] * th[e];
  const double amp = 1 / std::sqrt(12 * std::numbers::pi * vol);
  return {amp * std::cos(phi), amp};
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream z(20261023, 2, 0, 0);

  // ---- known values ----------------------------------------------------------------------------------------------------
  s.require("recoupling: <1/2 1/2 1/2 -1/2 | 1 0> = 1/sqrt 2 and <1/2 -1/2 1/2 1/2 | 0 0> = -1/sqrt 2",
            ps::clebsch_gordan(1, 1, 1, -1, 2, 0) == SignedSqrt{1, Rational(1, 2)} &&
                ps::clebsch_gordan(1, -1, 1, 1, 0, 0) == SignedSqrt{-1, Rational(1, 2)},
            ps::clebsch_gordan(1, -1, 1, 1, 0, 0).str());
  s.require("recoupling: <1 1 1 -1 | 0 0> = 1/sqrt 3", ps::clebsch_gordan(2, 2, 2, -2, 0, 0) == SignedSqrt{1, Rational(1, 3)},
            ps::clebsch_gordan(2, 2, 2, -2, 0, 0).str());
  s.require("recoupling: {1 1 1; 1 1 1} = 1/6", ps::wigner_6j(2, 2, 2, 2, 2, 2) == SignedSqrt{1, Rational(1, 36)}, ps::wigner_6j(2, 2, 2, 2, 2, 2).str());
  {
    // {a b c; 0 c b} = (−1)^{a+b+c} / √((2b+1)(2c+1)) (Edmonds eq. 6.3.2).
    bool ok = true;
    int n = 0;
    for (int a = 0; a <= 8; a++)
      for (int b = 0; b <= 8; b++)
        for (int c = 0; c <= 8; c++) {
          if (!ps::triangle(a, b, c)) continue;
          n++;
          const SignedSqrt want{phase(a + b + c), Rational(1, (b + 1) * (c + 1))};
          ok = ok && ps::wigner_6j(a, b, c, 0, c, b) == want;
        }
    s.require("recoupling: {a b c; 0 c b} = (-1)^(a+b+c) / sqrt((2b+1)(2c+1)), all " + std::to_string(n) + " admissible triads with 2j <= 8", ok);
  }
  // ---- Clebsch–Gordan orthogonality, exactly ---------------------------------------------------------------------------
  {
    bool ok = true;
    int sums = 0;
    for (int j1 = 0; j1 <= 5; j1++)
      for (int j2 = 0; j2 <= 5; j2++)
        for (int J = std::abs(j1 - j2); J <= j1 + j2; J += 2)
          for (int K = std::abs(j1 - j2); K <= j1 + j2; K += 2)
            for (int M = -std::min(J, K); M <= std::min(J, K); M += 2) {
              std::vector<SignedSqrt> t;
              for (int m1 = -j1; m1 <= j1; m1 += 2) t.push_back(ps::clebsch_gordan(j1, m1, j2, M - m1, J, M) * ps::clebsch_gordan(j1, m1, j2, M - m1, K, M));
              ok = ok && ps::sum(t) == (J == K ? SignedSqrt{1, Rational(1)} : SignedSqrt{});
              sums++;
            }
    s.require("recoupling: sum over m1 of <j1 m1 j2 m2|J M><j1 m1 j2 m2|K M> = delta_JK exactly, " + std::to_string(sums) + " sums, 2j <= 5", ok);
  }
  // ---- symmetries of the 3j symbol, exhaustively for 2j <= 4 ------------------------------------------------------------
  {
    bool ok = true;
    int n = 0;
    for (int a = 0; a <= 4; a++)
      for (int b = 0; b <= 4; b++)
        for (int c = 0; c <= 4; c++) {
          if (!ps::triangle(a, b, c)) continue;
          const int p = phase(a + b + c);
          for (int ma = -a; ma <= a; ma += 2)
            for (int mb = -b; mb <= b; mb += 2) {
              const int mc = -ma - mb;
              if (std::abs(mc) > c) continue;
              const SignedSqrt w = ps::wigner_3j(a, b, c, ma, mb, mc);
              if (w.is_zero()) continue;
              n++;
              ok = ok && ps::wigner_3j(b, c, a, mb, mc, ma) == w;                   // cyclic
              ok = ok && ps::wigner_3j(b, a, c, mb, ma, mc) == SignedSqrt{p * w.sign, w.square};    // odd permutation
              ok = ok && ps::wigner_3j(a, b, c, -ma, -mb, -mc) == SignedSqrt{p * w.sign, w.square};  // reflection
            }
        }
    s.require("recoupling: 3j symbols cyclic, and (-1)^(j1+j2+j3) under an odd permutation and under m -> -m, all " + std::to_string(n) +
                  " nonzero symbols with 2j <= 4",
              ok);
  }
  // ---- 6j symbols: orthogonality, symmetry, Racah, Biedenharn–Elliott --------------------------------------------------
  {
    bool ok = true;
    int sums = 0;
    for (int a = 0; a <= 4; a++)
      for (int b = 0; b <= 4; b++)
        for (int c = 0; c <= 4; c++)
          for (int d = 0; d <= 4; d++)
            for (int p = 0; p <= 8; p++)
              for (int q = 0; q <= 8; q++) {
                if (!ps::triangle(a, d, p) || !ps::triangle(c, b, p) || !ps::triangle(a, d, q) || !ps::triangle(c, b, q)) continue;
                std::vector<SignedSqrt> t;
                for (int x = 0; x <= 8; x++) t.push_back(Rational((x + 1) * (p + 1)) * (ps::wigner_6j(a, b, x, c, d, p) * ps::wigner_6j(a, b, x, c, d, q)));
                ok = ok && ps::sum(t) == (p == q ? SignedSqrt{1, Rational(1)} : SignedSqrt{});
                sums++;
              }
    s.require("recoupling: sum over x of (2x+1)(2p+1){a b x; c d p}{a b x; c d q} = delta_pq exactly, " + std::to_string(sums) + " sums", ok);

    bool sym = true, racah = true, be = true;
    int nontrivial_racah = 0, nontrivial_be = 0;
    for (int k = 0; k < 3000; k++) {
      std::array<int, 6> j;
      for (auto& v : j) v = static_cast<int>(z.below(7));
      const auto [a, b, c, d, e, f] = j;
      // Tetrahedral symmetry: permutations of the columns, and exchange of upper and lower in two columns.
      const SignedSqrt w = ps::wigner_6j(a, b, c, d, e, f);
      sym = sym && ps::wigner_6j(b, a, c, e, d, f) == w && ps::wigner_6j(a, c, b, d, f, e) == w && ps::wigner_6j(d, e, c, a, b, f) == w &&
            ps::wigner_6j(a, e, f, d, b, c) == w;
    }
    // Racah's sum rule, Σ_x (−1)^{p+q+x} (2x+1) {a b x; c d p}{a b x; d c q} = {a d p; b c q}, for every
    // a, b, c, d with 2j <= 4 and every p, q.
    for (int a = 0; a <= 4; a++)
      for (int b = 0; b <= 4; b++)
        for (int c = 0; c <= 4; c++)
          for (int d = 0; d <= 4; d++)
            for (int p = 0; p <= 8; p++)
              for (int q = 0; q <= 8; q++) {
                std::vector<SignedSqrt> t;
                bool parity = true;
                for (int x = 0; x <= 8; x++) {
                  const SignedSqrt term = ps::wigner_6j(a, b, x, c, d, p) * ps::wigner_6j(a, b, x, d, c, q);
                  if (term.is_zero()) continue;
                  if ((p + q + x) % 2) { parity = false; continue; }
                  t.push_back(Rational(phase(p + q + x) * (x + 1)) * term);
                }
                const SignedSqrt rhs = ps::wigner_6j(a, d, p, b, c, q);
                racah = racah && parity && ps::sum(t) == rhs;
                nontrivial_racah += !rhs.is_zero();
              }
    // Biedenharn–Elliott, Σ_x (−1)^{S+x} (2x+1) {a b x; c d p}{c d x; e f q}{e f x; b a r} = {p q r; e a d}{p q r; f b c},
    // S the sum of the nine spins: for 300 random a … f with 2j <= 5, every p, q, r.
    for (int k = 0; k < 300; k++) {
      std::array<int, 6> j;
      for (auto& v : j) v = static_cast<int>(z.below(6));
      const auto [a, b, c, d, e, f] = j;
      for (int p = 0; p <= 10; p++)
        for (int q = 0; q <= 10; q++)
          for (int r = 0; r <= 10; r++) {
            if (!ps::triangle(a, d, p) || !ps::triangle(b, c, p) || !ps::triangle(e, d, q) || !ps::triangle(f, c, q)) continue;
            const int S = a + b + c + d + e + f + p + q + r;
            std::vector<SignedSqrt> t;
            bool parity = true;
            for (int x = 0; x <= 10; x++) {
              const SignedSqrt term = ps::wigner_6j(a, b, x, c, d, p) * ps::wigner_6j(c, d, x, e, f, q) * ps::wigner_6j(e, f, x, b, a, r);
              if (term.is_zero()) continue;
              if ((S + x) % 2) { parity = false; continue; }
              t.push_back(Rational(phase(S + x) * (x + 1)) * term);
            }
            const SignedSqrt rhs = ps::wigner_6j(p, q, r, e, a, d) * ps::wigner_6j(p, q, r, f, b, c);
            be = be && parity && ps::sum(t) == rhs;
            nontrivial_be += !rhs.is_zero();
          }
    }
    s.require("recoupling: tetrahedral symmetry of 6j symbols, 3000 random arguments with 2j <= 6", sym);
    s.require("recoupling: Racah's sum rule exactly, every argument with 2j <= 4 (" + std::to_string(nontrivial_racah) + " with a nonzero side)",
              racah && nontrivial_racah > 100);
    s.require("recoupling: the Biedenharn-Elliott identity exactly, 300 random a ... f and every p, q, r (" + std::to_string(nontrivial_be) +
                  " with a nonzero side)",
              be && nontrivial_be > 100);
  }
  // ---- 9j symbols ------------------------------------------------------------------------------------------------------
  {
    double worst = 0;
    bool transpose = true, zero = true;
    int n = 0, nz = 0;
    // Arguments built from admissible rows and columns: j1, j2, j4, j5 drawn, the rest in their triangles.
    for (int k = 0; k < 20000 && n < 400; k++) {
      const int j1 = static_cast<int>(z.below(5)), j2 = static_cast<int>(z.below(5)), j4 = static_cast<int>(z.below(5)), j5 = static_cast<int>(z.below(5));
      auto pick = [&](int x, int y) { return std::abs(x - y) + 2 * static_cast<int>(z.below((x + y - std::abs(x - y)) / 2 + 1)); };
      const int j3 = pick(j1, j2), j6 = pick(j4, j5), j7 = pick(j1, j4), j8 = pick(j2, j5), j9 = pick(j3, j6);
      if (!ps::triangle(j7, j8, j9)) continue;
      const SignedSqrt w = ps::wigner_9j(j1, j2, j3, j4, j5, j6, j7, j8, j9);
      if (w.is_zero()) continue;
      n++;
      double sum = 0;
      for (int x = 0; x <= 16; x++)
        sum += (x % 2 ? -1 : 1) * (x + 1) * ps::wigner_6j(j1, j4, j7, j8, j9, x).value() * ps::wigner_6j(j2, j5, j8, j4, x, j6).value() *
               ps::wigner_6j(j3, j6, j9, x, j1, j2).value();
      worst = std::max(worst, std::fabs(sum - w.value()));
      transpose = transpose && ps::wigner_9j(j1, j4, j7, j2, j5, j8, j3, j6, j9) == w;
    }
    // {a b e; c d e; f f 0} = (−1)^{b+c+e+f} {a b e; d c f} / √((2e+1)(2f+1)) (Edmonds eq. 6.4.14).
    for (int a = 0; a <= 4; a++)
      for (int b = 0; b <= 4; b++)
        for (int c = 0; c <= 4; c++)
          for (int d = 0; d <= 4; d++)
            for (int e = 0; e <= 4; e++)
              for (int f = 0; f <= 4; f++) {
                const SignedSqrt six = ps::wigner_6j(a, b, e, d, c, f);
                const SignedSqrt nine = ps::wigner_9j(a, b, e, c, d, e, f, f, 0);
                if (six.is_zero()) { zero = zero && nine.is_zero(); continue; }
                nz++;
                const SignedSqrt want{phase(b + c + e + f) * six.sign, six.square / Rational((e + 1) * (f + 1))};
                zero = zero && nine == want;
              }
    s.check("recoupling: exact 9j symbols against their sum of three 6j symbols in floating point, " + std::to_string(n) + " symbols", worst, 1e-14);
    s.require("recoupling: 9j symbols invariant under transposition", transpose);
    s.require("recoupling: {a b e; c d e; f f 0} = (-1)^(b+c+e+f) {a b e; d c f} / sqrt((2e+1)(2f+1)), " + std::to_string(nz) + " nonzero cases", zero);
  }
  // ---- the Ponzano–Regge asymptotics -----------------------------------------------------------------------------------
  {
    // At spins λ j_e the error of the formula, relative to its amplitude, falls as 1/λ
    // (Roberts 1999, Geom. Topol. 3, 21–66; the next order is computed by Dupuis and Livine
    // 2009, Phys. Rev. D 80, 024035). Measured as the largest error over a window of eight
    // consecutive scales, so that a zero of the cosine does not decide it.
    const std::array<std::array<int, 6>, 2> shapes{{{2, 2, 2, 2, 2, 2}, {2, 3, 3, 2, 3, 3}}};
    std::array<double, 2> err_small{}, err_large{};
    for (int sh = 0; sh < 2; sh++) {
      for (int pass = 0; pass < 2; pass++) {
        const int lam0 = pass == 0 ? 20 : 40;
        double worst = 0;
        for (int lam = lam0; lam < lam0 + 8; lam++) {
          std::array<int, 6> j;
          for (int e = 0; e < 6; e++) j[e] = shapes[sh][e] * lam;
          const auto [pr, amp] = ponzano_regge(j);
          const double exact = ps::wigner_6j(j[0], j[1], j[2], j[3], j[4], j[5]).value();
          worst = std::max(worst, std::fabs(exact - pr) / amp);
        }
        (pass == 0 ? err_small : err_large)[sh] = worst;
      }
    }
    s.check("recoupling: Ponzano-Regge formula, regular tetrahedron, spins 20-27 (error relative to the amplitude)", err_small[0], 0.05);
    s.check("recoupling: Ponzano-Regge formula, irregular tetrahedron, spins 20-27 x (1, 1.5) (error relative to the amplitude)", err_small[1], 0.05);
    s.check("recoupling: the Ponzano-Regge error halves when the spins double: error(40-47) / error(20-27) - 1/2",
            std::max(std::fabs(err_large[0] / err_small[0] - 0.5), std::fabs(err_large[1] / err_small[1] - 0.5)), 0.15);
  }
  // ---- quantum 6j symbols ----------------------------------------------------------------------------------------------
  {
    double w2 = 0;
    for (int r = 3; r <= 30; r++) {
      double sum = 0;
      for (int two_j = 0; two_j <= r - 2; two_j++) sum += ps::quantum_dimension(two_j, r) * ps::quantum_dimension(two_j, r);
      w2 = std::max(w2, std::fabs(sum / ps::total_quantum_dimension_squared(r) - 1));
    }
    s.check("recoupling: total quantum dimension sum [2j+1]^2 = r / (2 sin^2(pi/r)), r = 3..30 (relative)", w2, 1e-13);

    double orth = 0;
    int sums = 0;
    for (int r : {5, 8, 13})
      for (int a = 0; a <= std::min(r - 2, 5); a++)
        for (int b = 0; b <= std::min(r - 2, 5); b++)
          for (int c = 0; c <= std::min(r - 2, 5); c++)
            for (int d = 0; d <= std::min(r - 2, 5); d++)
              for (int p = 0; p <= r - 2; p++)
                for (int q = 0; q <= r - 2; q++) {
                  if (!ps::admissible(a, d, p, r) || !ps::admissible(c, b, p, r) || !ps::admissible(a, d, q, r) || !ps::admissible(c, b, q, r)) continue;
                  double sum = 0;
                  for (int x = 0; x <= r - 2; x++)
                    sum += ps::quantum_dimension(x, r) * ps::quantum_dimension(p, r) * ps::quantum_6j(r, a, b, x, c, d, p) * ps::quantum_6j(r, a, b, x, c, d, q);
                  orth = std::max(orth, std::fabs(sum - (p == q ? 1 : 0)));
                  sums++;
                }
    s.check("recoupling: quantum 6j orthogonality, r = 5, 8, 13, " + std::to_string(sums) + " sums", orth, 1e-12);

    double be = 0, sym = 0;
    int nontrivial = 0;
    for (int r : {6, 7, 10}) {
      for (int k = 0; k < 100; k++) {
        std::array<int, 6> j;
        for (auto& v : j) v = static_cast<int>(z.below(std::min(r - 1, 6)));
        const auto [a, b, c, d, e, f] = j;
        const double w = ps::quantum_6j(r, a, b, c, d, e, f);
        sym = std::max({sym, std::fabs(ps::quantum_6j(r, b, a, c, e, d, f) - w), std::fabs(ps::quantum_6j(r, a, e, f, d, b, c) - w)});
        for (int p = 0; p <= r - 2; p++)
          for (int q = 0; q <= r - 2; q++)
            for (int rr = 0; rr <= r - 2; rr++) {
              if (!ps::admissible(a, d, p, r) || !ps::admissible(b, c, p, r) || !ps::admissible(e, d, q, r) || !ps::admissible(f, c, q, r)) continue;
              const int S = a + b + c + d + e + f + p + q + rr;
              double lhs = 0;
              for (int x = 0; x <= r - 2; x++) {
                const double term = ps::quantum_6j(r, a, b, x, c, d, p) * ps::quantum_6j(r, c, d, x, e, f, q) * ps::quantum_6j(r, e, f, x, b, a, rr);
                if (term == 0) continue;
                lhs += ((S + x) / 2 % 2 ? -1 : 1) * ps::quantum_dimension(x, r) * term;
              }
              const double rhs = ps::quantum_6j(r, p, q, rr, e, a, d) * ps::quantum_6j(r, p, q, rr, f, b, c);
              be = std::max(be, std::fabs(lhs - rhs));
              nontrivial += rhs != 0;
            }
      }
    }
    s.check("recoupling: quantum Biedenharn-Elliott identity, r = 6, 7, 10, 100 random a ... f each and every p, q, r (" + std::to_string(nontrivial) +
                " with a nonzero side)",
            nontrivial > 100 ? be : 1.0, 1e-12);
    s.check("recoupling: tetrahedral symmetry of quantum 6j symbols", sym, 1e-14);

    double lim = 0;
    for (int a = 0; a <= 6; a++)
      for (int b = 0; b <= 6; b++)
        for (int c = 0; c <= 6; c++)
          for (int d = 0; d <= 6; d++)
            for (int e = 0; e <= 6; e++)
              for (int f = 0; f <= 6; f++) lim = std::max(lim, std::fabs(ps::quantum_6j(20000, a, b, c, d, e, f) - ps::wigner_6j(a, b, c, d, e, f).value()));
    s.check("recoupling: quantum 6j symbols at r = 20000 against the classical ones, 2j <= 6 (error falls as 1/r^2)", lim, 1e-7);
  }

  // Invariants of tensor products: Catalan numbers for spins ½, the triangle conditions for three
  // spins, and for four the intermediate spins of a recoupling.
  {
    bool catalan = true;
    morphe::arithmetic::Integer c(1);  // C_k
    for (int n = 0; n <= 40; n += 2) {
      const int k = n / 2;
      if (k > 0) c = c * morphe::arithmetic::Integer(2 * (2 * k - 1)) / morphe::arithmetic::Integer(k + 1);
      catalan = catalan && ps::invariant_dimension(std::vector<int>(static_cast<std::size_t>(n), 1)) == c &&
                ps::invariant_dimension(std::vector<int>(static_cast<std::size_t>(n + 1), 1)).is_zero();
    }
    s.require("invariants: n spins ½ have the Catalan number C_{n/2} of them, n odd none (n <= 41)", catalan, "C_20 = " + c.str());
    bool three = true, four = true;
    for (int a = 0; a <= 8; a++)
      for (int b = 0; b <= 8; b++)
        for (int d = 0; d <= 8; d++) {
          three = three && ps::invariant_dimension({a, b, d}) == morphe::arithmetic::Integer(ps::triangle(a, b, d) ? 1 : 0);
          for (int e = 0; e <= 6; e++) {
            long long channels = 0;
            for (int x = 0; x <= 16; x++) channels += ps::triangle(a, b, x) && ps::triangle(x, d, e);
            four = four && ps::invariant_dimension({a, b, d, e}) == morphe::arithmetic::Integer(channels);
          }
        }
    s.require("invariants: three spins have one invariant exactly when they form a triad", three);
    s.require("invariants: four spins have one invariant for each intermediate spin of a recoupling", four);
  }
  return s.done();
}
