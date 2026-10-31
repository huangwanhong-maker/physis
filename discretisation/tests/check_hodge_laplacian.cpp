// Checks of the Hodge Laplacians of every degree (hodge_laplacian.hpp).
//
// - Hodge's theorem: with random weights the harmonic dimension of Δ_p is b_p — on a 3 × 2 × 2 cubical grid
//   (a ball: 1 0 0 0), on a 3 × 4 lattice torus (1 2 1), on the filled 3-simplex (1 0 0 0) and its 2-skeleton
//   (1 0 1), whose extra harmonic 2-cochain is what the dropped top cell filled.
// - Δ_p is self-adjoint and non-negative for ⟨·,·⟩_p, ⟨φ, Δφ⟩ = ‖d*φ‖² + ‖dφ‖², and its eigenvectors are
//   orthonormal for ⟨·,·⟩_p with Δu = λu.
// - Heat in four rooms (a 2 × 2 × 1 grid of cubes): with the stars 1/G on the faces and 1/C on the rooms,
//   Δ₃ Q = d G dᵀ (Q/C), the conduction of heat contents; the total flux Σ_c (Δ₃ Q)_c is carried by the outer
//   faces alone, Σ_f G_f T_c(f), the interior faces cancelling.

#include <cmath>
#include <map>
#include <vector>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <morphe/topology/constructions.hpp>
#include <morphe/topology/regular_complex.hpp>
#include <physis/discretisation/hodge_laplacian.hpp>

using calx::verification::Suite;
namespace pd = physis::discretisation;
namespace mt = morphe::topology;

namespace {

std::vector<std::vector<double>> random_weights(const morphe::algebra::ChainComplex& C, calx::random::Stream& rng) {
  std::vector<std::vector<double>> w(C.n.size());
  for (std::size_t q = 0; q < C.n.size(); q++)
    for (std::size_t i = 0; i < C.n[q]; i++) w[q].push_back(0.2 + 3 * rng.uniform());
  return w;
}

std::string harmonic(const morphe::algebra::ChainComplex& C, const std::vector<std::vector<double>>& w) {
  std::string s;
  for (std::size_t p = 0; p < C.n.size(); p++) s += (s.empty() ? "" : " ") + std::to_string(pd::hodge_laplacian(C, p, w).harmonic());
  return s;
}

mt::RegularComplex simplex3(bool with_top) {
  mt::RegularComplex X;
  std::map<unsigned, std::size_t> cell;
  for (std::size_t k = 0; k <= 3; k++)
    for (unsigned S = 1; S < 16; S++) {
      if (static_cast<std::size_t>(__builtin_popcount(S)) != k + 1 || (k == 3 && !with_top)) continue;
      std::vector<std::size_t> facets;
      if (k > 0)
        for (unsigned b = 0; b < 4; b++)
          if (S & (1u << b)) facets.push_back(cell.at(S ^ (1u << b)));
      cell[S] = X.add(k, "", facets);
    }
  mt::orient(X);
  return X;
}

double inner(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& w) {
  double s = 0;
  for (std::size_t i = 0; i < a.size(); i++) s += w[i] * a[i] * b[i];
  return s;
}
std::vector<double> mul(const calx::linalg::Matrix& M, const std::vector<double>& v) {
  std::vector<double> out(M.n, 0.0);
  for (std::size_t i = 0; i < M.n; i++)
    for (std::size_t j = 0; j < M.n; j++) out[i] += M(i, j) * v[j];
  return out;
}

}  // namespace

int main() {
  Suite s;
  calx::random::Stream rng(31, 0, 7, 1);
  const auto grid = mt::chain_complex(mt::cubical_grid(3, 2, 2).complex);
  s.require("Hodge's theorem: the grid of cubes, harmonic dimensions 1 0 0 0", harmonic(grid, random_weights(grid, rng)) == "1 0 0 0");
  const auto torus = mt::torus_lattice({3, 4}).chain_complex();
  s.require("the 3 × 4 lattice torus, 1 2 1", harmonic(torus, random_weights(torus, rng)) == "1 2 1", harmonic(torus, {}));
  const auto ball = mt::chain_complex(simplex3(true)), sphere = mt::chain_complex(simplex3(false));
  s.require("the filled 3-simplex 1 0 0 0, its 2-skeleton 1 0 1", harmonic(ball, random_weights(ball, rng)) == "1 0 0 0" && harmonic(sphere, random_weights(sphere, rng)) == "1 0 1");

  {
    const auto w = random_weights(grid, rng);
    double adjoint = 0, energy = 0, negative = 0, eigen = 0, ortho = 0;
    for (std::size_t p = 0; p < grid.n.size(); p++) {
      const auto H = pd::hodge_laplacian(grid, p, w);
      const std::size_t n = grid.n[p];
      std::vector<double> a(n), b(n);
      for (auto& x : a) x = rng.normal();
      for (auto& x : b) x = rng.normal();
      const auto La = mul(H.laplacian, a), Lb = mul(H.laplacian, b);
      adjoint = std::max(adjoint, std::fabs(inner(b, La, w[p]) - inner(Lb, a, w[p])));
      // ‖dφ‖² in ⟨·,·⟩_{p+1} and ‖d*φ‖² in ⟨·,·⟩_{p−1}.
      double up = 0, down = 0;
      if (p + 1 < grid.n.size()) {
        const auto d = pd::coboundary(grid, p);
        for (std::size_t t = 0; t < d.rows; t++) {
          double v = 0;
          for (std::size_t i = 0; i < n; i++) v += d(t, i) * a[i];
          up += w[p + 1][t] * v * v;
        }
      }
      if (p > 0) {
        const auto d = pd::coboundary(grid, p - 1);
        for (std::size_t q = 0; q < d.cols; q++) {
          double v = 0;
          for (std::size_t i = 0; i < n; i++) v += d(i, q) * w[p][i] * a[i];
          v /= w[p - 1][q];
          down += w[p - 1][q] * v * v;
        }
      }
      energy = std::max(energy, std::fabs(inner(a, La, w[p]) - up - down) / (1 + up + down));
      for (double v : H.values) negative = std::min(negative, v);
      for (std::size_t k = 0; k < n; k += std::max<std::size_t>(1, n / 7)) {
        std::vector<double> u(n);
        for (std::size_t i = 0; i < n; i++) u[i] = H.vectors(i, k);
        const auto Lu = mul(H.laplacian, u);
        for (std::size_t i = 0; i < n; i++) eigen = std::max(eigen, std::fabs(Lu[i] - H.values[k] * u[i]));
        for (std::size_t l = 0; l < n; l += std::max<std::size_t>(1, n / 5)) {
          std::vector<double> v(n);
          for (std::size_t i = 0; i < n; i++) v[i] = H.vectors(i, l);
          ortho = std::max(ortho, std::fabs(inner(u, v, w[p]) - (k == l ? 1.0 : 0.0)));
        }
      }
    }
    s.check("Δ_p is self-adjoint for ⟨·,·⟩_p", adjoint, 1e-10);
    s.check("⟨φ, Δφ⟩ = ‖d*φ‖² + ‖dφ‖² (relative)", energy, 1e-12);
    s.check("its eigenvalues are non-negative", negative, 1e-10);
    s.check("Δu = λu", eigen, 1e-9);
    s.check("the eigenvectors are orthonormal for ⟨·,·⟩_p", ortho, 1e-10);
  }
  {
    const auto G = mt::cubical_grid(2, 2, 1);
    const auto C = mt::chain_complex(G.complex);
    const auto& X = G.complex;
    std::vector<std::vector<double>> w(4);
    std::vector<double> cond, cap;
    for (std::size_t i = 0; i < X.size(); i++) {
      if (X.cells[i].dim == 2) cond.push_back(0.5 + rng.uniform());
      if (X.cells[i].dim == 3) cap.push_back(2 + rng.uniform());
    }
    w[2].resize(cond.size()), w[3].resize(cap.size());
    for (std::size_t f = 0; f < cond.size(); f++) w[2][f] = 1 / cond[f];
    for (std::size_t c = 0; c < cap.size(); c++) w[3][c] = 1 / cap[c];
    const auto H = pd::hodge_laplacian(C, 3, w);
    std::vector<double> Q(cap.size()), T(cap.size());
    for (std::size_t c = 0; c < cap.size(); c++) T[c] = 10 + 5 * rng.normal(), Q[c] = cap[c] * T[c];
    const auto LQ = mul(H.laplacian, Q);
    const auto d = pd::coboundary(C, 2);  // rooms × faces
    double worst = 0, outer = 0, total = 0;
    for (std::size_t c = 0; c < cap.size(); c++) {
      double v = 0;
      for (std::size_t f = 0; f < cond.size(); f++) {
        if (d(c, f) == 0) continue;
        double flux = 0;
        for (std::size_t c2 = 0; c2 < cap.size(); c2++) flux += d(c2, f) * T[c2];
        v += d(c, f) * cond[f] * flux;
      }
      worst = std::max(worst, std::fabs(LQ[c] - v));
      total += LQ[c];
    }
    for (std::size_t f = 0; f < cond.size(); f++) {
      std::size_t rooms = 0, room = 0;
      for (std::size_t c = 0; c < cap.size(); c++)
        if (d(c, f) != 0) rooms++, room = c;
      if (rooms == 1) outer += cond[f] * T[room];
    }
    s.check("heat in four rooms: Δ₃ Q = d G dᵀ (Q/C), conduction of heat contents", worst, 1e-10);
    s.check("the total flux is carried by the outer faces alone", total - outer, 1e-9);
  }
  return s.done();
}
