// The Levi-Civita connection and the curvature of a metric field, at a point.
//
// Given the components g_ab of a metric in a chart and their first and second partial
// derivatives, which calx's jets supply exactly to rounding (calx decision 0005), the
// tensors of curvature follow algebraically (Misner, Thorne and Wheeler 1973, Gravitation,
// §14.1 and Box 14.2; Carroll 2004, Spacetime and Geometry, §3.2-3.7). With the sign
// conventions of physis decision 0005, §2.4 (MTW's "+ + +"):
//
//   Christoffel symbols   Γ^a_bc = ½ g^ad (∂_b g_dc + ∂_c g_db - ∂_d g_bc)               (MTW 8.24b)
//   Riemann tensor        R^a_bcd = ∂_c Γ^a_db - ∂_d Γ^a_cb + Γ^a_ce Γ^e_db - Γ^a_de Γ^e_cb  (MTW 11.12)
//   Ricci tensor          R_bd = R^a_bad,   scalar R = g^bd R_bd                           (MTW 14.4-5)
//   Einstein tensor       G_ab = R_ab - ½ R g_ab                                          (MTW 14.7)
//   Weyl tensor (n >= 3)  C_abcd = R_abcd - 2/(n-2) (g_a[c R_d]b - g_b[c R_d]a)
//                                 + 2/((n-1)(n-2)) R g_a[c g_d]b                          (Wald 1984, 3.2.28)
//   Kretschmann scalar    K = R_abcd R^abcd
//
// with brackets antisymmetrising with weight ½. The derivatives of the Christoffel symbols
// need those of the inverse metric, ∂_e g^ad = -g^am (∂_e g_mn) g^nd. Every tensor is
// computed at one point of one chart; its accuracy is that of the metric's own arithmetic,
// a few units of rounding in the components, with the cancellation that a flat or vacuum
// metric's curvature (zero) suffers measured against the size of its terms.

#pragma once

#include <array>
#include <cstddef>

#include <calx/differentiation/jet.hpp>

#include "metric.hpp"

namespace physis::geometry {

template <std::size_t n>
struct Curvature {
  Tensor2<n> g{}, ginv{};  // g_ab, g^ab
  Tensor3<n> dg{};         // dg[c][a][b] = ∂_c g_ab
  Tensor3<n> Gamma{};      // Gamma[a][b][c] = Γ^a_bc
  Tensor4<n> Riemann{};    // Riemann[a][b][c][d] = R^a_bcd
  Tensor2<n> Ricci{};      // R_bd
  double scalar = 0;       // R
  Tensor2<n> Einstein{};   // G_ab
};

// The components of g at x.
template <MetricField M>
inline Tensor2<M::n> metric_at(const M& metric, const Point<M::n>& x) { return metric(x); }

// The Christoffel symbols Γ^a_bc at x, from first derivatives only (dual numbers).
template <MetricField M>
inline Tensor3<M::n> christoffel(const M& metric, const Point<M::n>& x) {
  constexpr std::size_t n = M::n;
  const auto z = calx::differentiation::variables<n, 1>(x);
  const auto G = metric(z);
  Tensor2<n> g;
  Tensor3<n> dg;  // dg[c][a][b]
  for (std::size_t a = 0; a < n; a++)
    for (std::size_t b = 0; b < n; b++) {
      g[a][b] = G[a][b].v;
      for (std::size_t c = 0; c < n; c++) dg[c][a][b] = G[a][b].d[c];
    }
  const Tensor2<n> gi = inverse(g);
  Tensor3<n> Gamma{};
  for (std::size_t a = 0; a < n; a++)
    for (std::size_t b = 0; b < n; b++)
      for (std::size_t c = b; c < n; c++) {
        double s = 0;
        for (std::size_t d = 0; d < n; d++) s += gi[a][d] * (dg[b][d][c] + dg[c][d][b] - dg[d][b][c]);
        Gamma[a][b][c] = Gamma[a][c][b] = 0.5 * s;
      }
  return Gamma;
}

// Every tensor of curvature at x, from the metric's jets of order 2.
template <MetricField M>
inline Curvature<M::n> curvature(const M& metric, const Point<M::n>& x) {
  constexpr std::size_t n = M::n;
  Curvature<n> k;
  const auto z = calx::differentiation::variables<n, 2>(x);
  const auto G = metric(z);
  Tensor4<n> ddg;  // ddg[e][c][a][b] = ∂_e ∂_c g_ab
  for (std::size_t a = 0; a < n; a++)
    for (std::size_t b = 0; b < n; b++) {
      k.g[a][b] = G[a][b].v;
      for (std::size_t c = 0; c < n; c++) {
        k.dg[c][a][b] = G[a][b].d[c];
        for (std::size_t e = 0; e < n; e++) ddg[e][c][a][b] = G[a][b].h[e][c];
      }
    }
  k.ginv = inverse(k.g);
  const auto& gi = k.ginv;
  const auto& dg = k.dg;
  // Christoffel symbols of the first kind, Γ_dbc, and their derivatives ∂_e Γ_dbc.
  Tensor3<n> first{};
  Tensor4<n> dfirst{};  // dfirst[e][d][b][c]
  for (std::size_t d = 0; d < n; d++)
    for (std::size_t b = 0; b < n; b++)
      for (std::size_t c = 0; c < n; c++) {
        first[d][b][c] = 0.5 * (dg[b][d][c] + dg[c][d][b] - dg[d][b][c]);
        for (std::size_t e = 0; e < n; e++) dfirst[e][d][b][c] = 0.5 * (ddg[e][b][d][c] + ddg[e][c][d][b] - ddg[e][d][b][c]);
      }
  // ∂_e g^ad = -g^am ∂_e g_mn g^nd.
  Tensor3<n> dgi{};  // dgi[e][a][d]
  for (std::size_t e = 0; e < n; e++)
    for (std::size_t a = 0; a < n; a++)
      for (std::size_t d = 0; d < n; d++) {
        double s = 0;
        for (std::size_t m = 0; m < n; m++)
          for (std::size_t q = 0; q < n; q++) s += gi[a][m] * dg[e][m][q] * gi[q][d];
        dgi[e][a][d] = -s;
      }
  // Γ^a_bc and ∂_e Γ^a_bc.
  Tensor4<n> dGamma{};  // dGamma[e][a][b][c]
  for (std::size_t a = 0; a < n; a++)
    for (std::size_t b = 0; b < n; b++)
      for (std::size_t c = 0; c < n; c++) {
        double s = 0;
        for (std::size_t d = 0; d < n; d++) s += gi[a][d] * first[d][b][c];
        k.Gamma[a][b][c] = s;
        for (std::size_t e = 0; e < n; e++) {
          double t = 0;
          for (std::size_t d = 0; d < n; d++) t += dgi[e][a][d] * first[d][b][c] + gi[a][d] * dfirst[e][d][b][c];
          dGamma[e][a][b][c] = t;
        }
      }
  // R^a_bcd = ∂_c Γ^a_db - ∂_d Γ^a_cb + Γ^a_ce Γ^e_db - Γ^a_de Γ^e_cb.
  const auto& Gm = k.Gamma;
  for (std::size_t a = 0; a < n; a++)
    for (std::size_t b = 0; b < n; b++)
      for (std::size_t c = 0; c < n; c++)
        for (std::size_t d = 0; d < n; d++) {
          double s = dGamma[c][a][d][b] - dGamma[d][a][c][b];
          for (std::size_t e = 0; e < n; e++) s += Gm[a][c][e] * Gm[e][d][b] - Gm[a][d][e] * Gm[e][c][b];
          k.Riemann[a][b][c][d] = s;
        }
  for (std::size_t b = 0; b < n; b++)
    for (std::size_t d = 0; d < n; d++) {
      double s = 0;
      for (std::size_t a = 0; a < n; a++) s += k.Riemann[a][b][a][d];
      k.Ricci[b][d] = s;
    }
  for (std::size_t b = 0; b < n; b++)
    for (std::size_t d = 0; d < n; d++) k.scalar += gi[b][d] * k.Ricci[b][d];
  for (std::size_t a = 0; a < n; a++)
    for (std::size_t b = 0; b < n; b++) k.Einstein[a][b] = k.Ricci[a][b] - 0.5 * k.scalar * k.g[a][b];
  return k;
}

// R_abcd = g_am R^m_bcd.
template <std::size_t n>
inline Tensor4<n> lowered_riemann(const Curvature<n>& k) {
  Tensor4<n> R{};
  for (std::size_t a = 0; a < n; a++)
    for (std::size_t b = 0; b < n; b++)
      for (std::size_t c = 0; c < n; c++)
        for (std::size_t d = 0; d < n; d++) {
          double s = 0;
          for (std::size_t m = 0; m < n; m++) s += k.g[a][m] * k.Riemann[m][b][c][d];
          R[a][b][c][d] = s;
        }
  return R;
}

// T^abcd from T_abcd, raising every index with g^ab.
template <std::size_t n>
inline Tensor4<n> raise_all(const Tensor4<n>& T, const Tensor2<n>& gi) {
  Tensor4<n> A = T, B{};
  // Raise one index at a time: O(n^5).
  for (int slot = 0; slot < 4; slot++) {
    for (std::size_t a = 0; a < n; a++)
      for (std::size_t b = 0; b < n; b++)
        for (std::size_t c = 0; c < n; c++)
          for (std::size_t d = 0; d < n; d++) {
            double s = 0;
            for (std::size_t m = 0; m < n; m++) {
              const std::size_t i[4] = {slot == 0 ? m : a, slot == 1 ? m : b, slot == 2 ? m : c, slot == 3 ? m : d};
              const std::size_t free = slot == 0 ? a : slot == 1 ? b : slot == 2 ? c : d;
              s += gi[free][m] * A[i[0]][i[1]][i[2]][i[3]];
            }
            B[a][b][c][d] = s;
          }
    A = B;
  }
  return A;
}

// The full contraction S_abcd T^abcd.
template <std::size_t n>
inline double full_contraction(const Tensor4<n>& S, const Tensor4<n>& T_upper) {
  double s = 0;
  for (std::size_t a = 0; a < n; a++)
    for (std::size_t b = 0; b < n; b++)
      for (std::size_t c = 0; c < n; c++)
        for (std::size_t d = 0; d < n; d++) s += S[a][b][c][d] * T_upper[a][b][c][d];
  return s;
}

// The Kretschmann scalar R_abcd R^abcd.
template <std::size_t n>
inline double kretschmann(const Curvature<n>& k) {
  const auto R = lowered_riemann(k);
  return full_contraction(R, raise_all(R, k.ginv));
}

// The Weyl tensor C_abcd, the trace-free part of the Riemann tensor (n >= 3).
template <std::size_t n>
inline Tensor4<n> weyl(const Curvature<n>& k) {
  static_assert(n >= 3, "physis::geometry::weyl: the Weyl tensor needs three dimensions or more");
  const auto R = lowered_riemann(k);
  const auto& g = k.g;
  const auto& Ric = k.Ricci;
  const double c1 = 2.0 / (n - 2), c2 = 2.0 / ((n - 1.0) * (n - 2.0));
  Tensor4<n> C{};
  for (std::size_t a = 0; a < n; a++)
    for (std::size_t b = 0; b < n; b++)
      for (std::size_t c = 0; c < n; c++)
        for (std::size_t d = 0; d < n; d++) {
          const double gR = 0.5 * (g[a][c] * Ric[d][b] - g[a][d] * Ric[c][b]) - 0.5 * (g[b][c] * Ric[d][a] - g[b][d] * Ric[c][a]);
          const double gg = 0.5 * (g[a][c] * g[d][b] - g[a][d] * g[c][b]);
          C[a][b][c][d] = R[a][b][c][d] - c1 * gR + c2 * k.scalar * gg;
        }
  return C;
}

// The mixed Einstein tensor G^a_b = g^ac G_cb.
template <std::size_t n>
inline Tensor2<n> einstein_mixed(const Curvature<n>& k) {
  Tensor2<n> m{};
  for (std::size_t a = 0; a < n; a++)
    for (std::size_t b = 0; b < n; b++)
      for (std::size_t c = 0; c < n; c++) m[a][b] += k.ginv[a][c] * k.Einstein[c][b];
  return m;
}

}  // namespace physis::geometry
