// Potential fields outside a sphere, from the coefficients of their spherical-harmonic expansion
// (decision 0010).
//
// A field B = −∇V that is harmonic outside a sphere of radius a (no sources there: ∇·B = 0 and
// ∇×B = 0) is the gradient of an internal potential
//
//   V(r, ψ, λ) = a Σₙ₌₁ᴺ (a/r)ⁿ⁺¹ Σₘ₌₀ⁿ (gₙᵐ cos mλ + hₙᵐ sin mλ) P̆ₙᵐ(sin ψ),
//
// with ψ the geocentric (planetocentric) latitude, λ the east longitude and P̆ₙᵐ the associated
// Legendre functions in Schmidt's semi-normalisation, without the Condon–Shortley phase (calx
// special_functions, legendre.hpp). This is the convention of geomagnetism, from Gauss to the IGRF,
// the WMM and the crustal fields of the Moon and Mars (Langel 1987, §4; Chulliat et al. 2025, §1.2);
// the gₙᵐ and hₙᵐ are in the field's own unit (tesla here), V in that unit times metres. Its
// components in the local frame of the sphere, north X, east Y and down Z, are
//
//   X = −Σₙ (a/r)ⁿ⁺² Σₘ (gₙᵐ cos mλ + hₙᵐ sin mλ) dP̆ₙᵐ/dψ,
//   Y =  Σₙ (a/r)ⁿ⁺² Σₘ m (gₙᵐ sin mλ − hₙᵐ cos mλ) P̆ₙᵐ / cos ψ,
//   Z = −Σₙ (n + 1) (a/r)ⁿ⁺² Σₘ (gₙᵐ cos mλ + hₙᵐ sin mλ) P̆ₙᵐ,
//
// (the WMM2025 report's eqs. 10–12 with its time dependence removed), finite at the poles through
// P̆ₙᵐ/cos ψ. The power of degree n at radius r is the Lowes–Mauersberger spectrum
//
//   Rₙ(r) = (n + 1) (a/r)²ⁿ⁺⁴ Σₘ [(gₙᵐ)² + (hₙᵐ)²],
//
// the mean over the sphere of |B|² carried by degree n (Lowes 1966): Σₙ Rₙ = ⟨|B|²⟩.
//
// A gravitational potential has the same form with g = ∇U and the fully normalised functions
// (P̄ₙᵐ = √(2n + 1) P̆ₙᵐ, coefficients C̄, S̄ scaled by GM/a); it is gravitation's to give them.
//
// Expansions are evaluated at a point (field_at) or on a grid of latitudes and longitudes
// (field_grid), the latter by rows: the Legendre functions once per latitude and the sum over the
// degrees once per order, so that a grid of L latitudes and M longitudes costs O(L N²) + O(L M N)
// rather than O(L M N²). Both hold to the degrees of the crustal models (449 for the Moon) at every
// latitude, the poles included, by calx's extended range.

#pragma once

#include <cmath>
#include <cstddef>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <calx/special_functions/legendre.hpp>

namespace physis::analysis {

// The coefficients of an internal potential, Schmidt semi-normalised, on the reference sphere a.
struct HarmonicPotential {
  int degree = 0;
  double radius = 0;          // a, m
  std::vector<double> g, h;   // at index(n, m); g[0] = h[0] = 0 (no monopole)

  static constexpr std::size_t index(int n, int m) { return std::size_t(n) * std::size_t(n + 1) / 2 + std::size_t(m); }
  static HarmonicPotential zero(int degree, double radius) {
    HarmonicPotential p;
    p.degree = degree;
    p.radius = radius;
    p.g.assign(index(degree, degree) + 1, 0.0);
    p.h.assign(p.g.size(), 0.0);
    return p;
  }
  // The expansion truncated to degree n ≤ degree.
  HarmonicPotential truncated(int n) const {
    if (n < 0 || n > degree) throw std::invalid_argument("physis::analysis: no such degree to truncate to");
    HarmonicPotential p = zero(n, radius);
    std::copy(g.begin(), g.begin() + std::ptrdiff_t(p.g.size()), p.g.begin());
    std::copy(h.begin(), h.begin() + std::ptrdiff_t(p.h.size()), p.h.begin());
    return p;
  }
};

// The sum of two expansions on one reference sphere, to the larger degree.
inline HarmonicPotential operator+(const HarmonicPotential& a, const HarmonicPotential& b) {
  if (a.radius != b.radius) throw std::invalid_argument("physis::analysis: expansions on different spheres");
  const HarmonicPotential& big = a.degree >= b.degree ? a : b;
  const HarmonicPotential& small = a.degree >= b.degree ? b : a;
  HarmonicPotential s = big;
  for (std::size_t k = 0; k < small.g.size(); k++) s.g[k] += small.g[k], s.h[k] += small.h[k];
  return s;
}

// SHTOOLS's text format with a header: "a, N" (a in metres, or in kilometres where the header
// says so, `header_km`), then lines "n m g h" for every degree and order; `scale` multiplies the
// coefficients into the field's unit (1 for a file in tesla, 1e-9 for one in nanotesla).
inline HarmonicPotential read_shtools(const std::string& path, double scale = 1, bool header_km = false) {
  std::ifstream in(path);
  if (!in) throw std::invalid_argument("physis::analysis: cannot read " + path);
  std::string line;
  if (!std::getline(in, line)) throw std::invalid_argument("physis::analysis: " + path + " is empty");
  for (char& c : line)
    if (c == ',') c = ' ';
  std::istringstream head(line);
  double a = 0;
  int n_max = 0;
  if (!(head >> a >> n_max) || !(a > 0) || n_max < 1) throw std::invalid_argument("physis::analysis: the header of " + path + " is not 'radius, degree'");
  HarmonicPotential p = HarmonicPotential::zero(n_max, header_km ? a * 1e3 : a);
  int n = 0, m = 0;
  double gv = 0, hv = 0;
  std::size_t seen = 0;
  while (in >> n >> m >> gv >> hv) {
    if (n < 0 || m < 0 || m > n || n > n_max) throw std::invalid_argument("physis::analysis: a row of " + path + " has no valid degree and order");
    p.g[HarmonicPotential::index(n, m)] = gv * scale;
    p.h[HarmonicPotential::index(n, m)] = hv * scale;
    seen++;
  }
  if (seen < HarmonicPotential::index(n_max, n_max)) throw std::invalid_argument("physis::analysis: " + path + " lacks coefficients");
  p.g[0] = p.h[0] = 0;
  return p;
}

// The potential and the field at a point, in the local frame of the sphere.
struct FieldValue {
  double potential = 0;            // V, in the field's unit times metres
  double north = 0, east = 0, down = 0;  // X, Y, Z
  double total() const { return std::sqrt(north * north + east * east + down * down); }
  double horizontal() const { return std::hypot(north, east); }
  double radial() const { return -down; }  // B_r
};

// At geocentric radius r (m), latitude ψ and east longitude λ (radians), r ≥ a for the expansion to hold.
inline FieldValue field_at(const HarmonicPotential& p, double r, double latitude, double longitude, int to_degree = -1) {
  const int N = to_degree < 0 ? p.degree : std::min(to_degree, p.degree);
  const calx::special_functions::SchmidtLegendre L = calx::special_functions::schmidt_legendre(N, latitude);
  const double q = p.radius / r;
  std::vector<double> cm(N + 1), sm(N + 1);
  for (int m = 0; m <= N; m++) cm[m] = std::cos(m * longitude), sm[m] = std::sin(m * longitude);
  FieldValue f;
  double power = q * q;  // (a/r)^(n+1) at n = 1
  for (int n = 1; n <= N; n++) {
    double v = 0, x = 0, y = 0;
    for (int m = 0; m <= n; m++) {
      const std::size_t k = HarmonicPotential::index(n, m);
      const double t = p.g[k] * cm[m] + p.h[k] * sm[m];
      v += t * L.p[k];
      x += t * L.dp[k];
      if (m) y += m * (p.g[k] * sm[m] - p.h[k] * cm[m]) * L.p_over_cos[k];
    }
    f.potential += p.radius * power * v;
    const double w = power * q;  // (a/r)^(n+2)
    f.north -= w * x;
    f.east += w * y;
    f.down -= (n + 1) * w * v;
    power *= q;
  }
  return f;
}

// The field on a grid of latitudes and longitudes (radians) at one radius: values[i * lons + j].
inline std::vector<FieldValue> field_grid(const HarmonicPotential& p, double r, const std::vector<double>& latitudes,
                                          const std::vector<double>& longitudes, int to_degree = -1) {
  const int N = to_degree < 0 ? p.degree : std::min(to_degree, p.degree);
  const double q = p.radius / r;
  std::vector<double> powers(N + 3);  // (a/r)ⁿ for n = 0 … N + 2
  powers[0] = 1;
  for (int n = 1; n <= N + 2; n++) powers[n] = powers[n - 1] * q;
  // cos mλ, sin mλ for each longitude, once, the orders of one longitude side by side (the sum over
  // the orders below runs along them).
  const std::size_t J = longitudes.size(), M = std::size_t(N) + 1;
  std::vector<double> cs(M * J), sn(M * J);
  for (std::size_t j = 0; j < J; j++)
    for (int m = 0; m <= N; m++) cs[j * M + std::size_t(m)] = std::cos(m * longitudes[j]), sn[j * M + std::size_t(m)] = std::sin(m * longitudes[j]);
  std::vector<FieldValue> out(latitudes.size() * J);
  // Per order m, the sums over the degrees of the coefficients of cos mλ and sin mλ in V, X, Y, Z.
  std::vector<double> vc(N + 1), vs(N + 1), xc(N + 1), xs(N + 1), yc(N + 1), ys(N + 1), zc(N + 1), zs(N + 1);
  for (std::size_t i = 0; i < latitudes.size(); i++) {
    const calx::special_functions::SchmidtLegendre L = calx::special_functions::schmidt_legendre(N, latitudes[i]);
    std::fill(vc.begin(), vc.end(), 0.0), std::fill(vs.begin(), vs.end(), 0.0), std::fill(xc.begin(), xc.end(), 0.0), std::fill(xs.begin(), xs.end(), 0.0);
    std::fill(yc.begin(), yc.end(), 0.0), std::fill(ys.begin(), ys.end(), 0.0), std::fill(zc.begin(), zc.end(), 0.0), std::fill(zs.begin(), zs.end(), 0.0);
    for (int n = 1; n <= N; n++) {
      const double pv = p.radius * powers[n + 1], w = powers[n + 2];
      for (int m = 0; m <= n; m++) {
        const std::size_t k = HarmonicPotential::index(n, m);
        const double gk = p.g[k], hk = p.h[k];
        vc[m] += pv * gk * L.p[k], vs[m] += pv * hk * L.p[k];
        xc[m] -= w * gk * L.dp[k], xs[m] -= w * hk * L.dp[k];
        if (m) yc[m] -= w * m * hk * L.p_over_cos[k], ys[m] += w * m * gk * L.p_over_cos[k];
        zc[m] -= (n + 1) * w * gk * L.p[k], zs[m] -= (n + 1) * w * hk * L.p[k];
      }
    }
    for (std::size_t j = 0; j < J; j++) {
      FieldValue f;
      const double* cj = &cs[j * M];
      const double* sj = &sn[j * M];
      for (int m = 0; m <= N; m++) {
        const double c = cj[m], s = sj[m];
        f.potential += vc[m] * c + vs[m] * s;
        f.north += xc[m] * c + xs[m] * s;
        f.east += yc[m] * c + ys[m] * s;
        f.down += zc[m] * c + zs[m] * s;
      }
      out[i * J + j] = f;
    }
  }
  return out;
}

// The Lowes–Mauersberger spectrum at radius r: R[n] for n = 0 … degree (R[0] = 0).
inline std::vector<double> lowes_spectrum(const HarmonicPotential& p, double r = 0) {
  if (r == 0) r = p.radius;
  const double q = p.radius / r;
  std::vector<double> R(p.degree + 1, 0.0);
  for (int n = 1; n <= p.degree; n++) {
    double s = 0;
    for (int m = 0; m <= n; m++) {
      const std::size_t k = HarmonicPotential::index(n, m);
      s += p.g[k] * p.g[k] + p.h[k] * p.h[k];
    }
    R[n] = (n + 1) * std::pow(q, 2 * n + 4) * s;
  }
  return R;
}

}  // namespace physis::analysis
