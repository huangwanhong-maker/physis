// Checks of the exact contraction of ℤ₂ gauge theories (contraction.hpp), and of couplings and
// frozen links in the gauge field (field.hpp, updates.hpp).
//
// - On a 3 × 3 lattice torus with a coupling for each plaquette, some negative, and three links
//   fixed: log Z, every ⟨h_f⟩ and ⟨u_e⟩, two-point functions and the law of the negative products
//   at a cut, against enumeration of all free links.
// - On closed surfaces with one coupling: the partition function and the mean plaquette against
//   Migdal and Witten's (exact.hpp), on a 4 × 4 torus and the surface of genus 2.
// - At strong coupling, log Z − β F tends to log of the number of flat connections, counted
//   through the fundamental group (exact.hpp).
// - With no link fixed, the spin networks crossing a cut are ℤ/2 1-cycles of the cut.
// - The derivatives carried by calx jets: the gradient of log Z is Σ ⟨h_f⟩ over each
//   parameter's plaquettes, and the Hessian matches differences of the gradient.
// - Metropolis with per-plaquette couplings and frozen links samples the same theory: the
//   frozen links never move, and the means agree with the contraction.

#include <physis/gauge/contraction.hpp>
#include <physis/gauge/exact.hpp>
#include <physis/gauge/field.hpp>
#include <physis/gauge/updates.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <morphe/algebra/finite_group.hpp>
#include <morphe/topology/constructions.hpp>

#include <array>
#include <cmath>
#include <map>
#include <numeric>
#include <string>
#include <vector>

namespace pg = physis::gauge;
namespace mt = morphe::topology;

namespace {

double max_difference(const std::vector<double>& a, const std::vector<double>& b) {
  double m = 0;
  for (std::size_t i = 0; i < a.size(); i++) m = std::max(m, std::fabs(a[i] - b[i]));
  return m;
}

// Everything by enumerating the free links.
struct Brute {
  double log_z = 0;
  std::vector<double> link, plaquette, law;
  std::vector<std::vector<double>> pair;
};

Brute brute(const pg::Z2Theory& T, const std::vector<std::vector<std::size_t>>& products = {}) {
  std::vector<std::size_t> free;
  std::vector<int> u(T.links, 1);
  for (std::size_t e = 0; e < T.links; e++) {
    if (T.fixed_value(e)) u[e] = T.fixed_value(e);
    else free.push_back(e);
  }
  Brute B;
  B.link.assign(T.links, 0);
  B.plaquette.assign(T.plaquettes.size(), 0);
  B.law.assign(products.size() + 1, 0);
  B.pair.assign(T.links, std::vector<double>(T.links, 0));
  double Z = 0;
  for (std::uint64_t mask = 0; mask < (std::uint64_t{1} << free.size()); mask++) {
    for (std::size_t i = 0; i < free.size(); i++) u[free[i]] = (mask >> i) & 1 ? -1 : 1;
    double logw = 0;
    std::vector<int> h(T.plaquettes.size(), 1);
    for (std::size_t f = 0; f < T.plaquettes.size(); f++) {
      for (std::size_t e : T.plaquettes[f]) h[f] *= u[e];
      logw += T.beta[f] * h[f];
    }
    const double w = std::exp(logw);
    Z += w;
    for (std::size_t e = 0; e < T.links; e++) B.link[e] += w * u[e];
    for (std::size_t f = 0; f < T.plaquettes.size(); f++) B.plaquette[f] += w * h[f];
    for (std::size_t x = 0; x < T.links; x++)
      for (std::size_t y = 0; y < T.links; y++) B.pair[x][y] += w * u[x] * u[y];
    std::size_t negative = 0;
    for (const auto& set : products) {
      int p = 1;
      for (std::size_t e : set) p *= u[e];
      negative += p < 0;
    }
    B.law[negative] += w;
  }
  B.log_z = std::log(Z);
  for (auto& x : B.link) x /= Z;
  for (auto& x : B.plaquette) x /= Z;
  for (auto& x : B.law) x /= Z;
  for (auto& row : B.pair)
    for (auto& x : row) x /= Z;
  return B;
}

std::vector<std::size_t> natural(std::size_t n) {
  std::vector<std::size_t> o(n);
  std::iota(o.begin(), o.end(), 0);
  return o;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // A 3 × 3 lattice torus with a coupling for each plaquette, some negative, and three links fixed.
  {
    const mt::CWComplex X = mt::torus_lattice({3, 3});
    std::vector<double> beta(X.faces.size());
    for (std::size_t f = 0; f < beta.size(); f++) beta[f] = -0.4 + 0.17 * static_cast<double>(f);
    std::vector<int> fixed(X.edges.size(), 0);
    fixed[0] = 1;
    fixed[5] = -1;
    fixed[11] = 1;
    const pg::Z2Theory T = pg::Z2Theory::of(X, beta, fixed);
    const auto order = natural(T.plaquettes.size());
    // Two-points between free links whose first plaquette comes no later than the other's last.
    std::vector<std::size_t> first(T.links, T.plaquettes.size()), last(T.links, 0);
    for (std::size_t k = 0; k < order.size(); k++)
      for (std::size_t e : T.plaquettes[order[k]]) {
        first[e] = std::min(first[e], k);
        last[e] = std::max(last[e], k);
      }
    std::vector<std::array<std::size_t, 2>> pairs;
    for (std::size_t x = 0; x < T.links; x++)
      for (std::size_t y = 0; y < T.links; y++)
        if (x != y && !fixed[x] && !fixed[y] && first[x] <= last[y] && first[x] < T.plaquettes.size()) pairs.push_back({x, y});
    pg::ContractionCut cut;
    cut.after = 4;
    cut.products = {{1, 2}, {3, 4, 6}, {0, 7}};
    for (const auto& set : cut.products) cut.keep.insert(cut.keep.end(), set.begin(), set.end());
    const pg::Contraction C = pg::contract(T, order, pairs, {cut});
    const Brute B = brute(T, cut.products);
    s.check("a 3 × 3 torus, couplings per plaquette, three links fixed: log Z against enumeration", C.log_z - B.log_z, 1e-10);
    s.check("the same: ⟨h_f⟩ of every plaquette", max_difference(C.plaquette, B.plaquette), 1e-10);
    s.check("the same: ⟨u_e⟩ of every link, the fixed ones their values", max_difference(C.link, B.link), 1e-10);
    double worst = 0;
    for (std::size_t k = 0; k < pairs.size(); k++) worst = std::max(worst, std::fabs(C.two_point[k] - B.pair[pairs[k][0]][pairs[k][1]]));
    s.check("the same: " + std::to_string(pairs.size()) + " two-point functions", worst, 1e-10);
    s.check("the same: the law of the negative products at a cut", max_difference(C.cuts[0].negative, B.law), 1e-10);
    double mean = 0;
    for (std::size_t k = 0; k < B.law.size(); k++) mean += static_cast<double>(k) * B.law[k];
    s.check("the same: their mean, by pairing alone", C.cuts[0].mean_negative - mean, 1e-10);
    s.check("log Z by the plain forward pass", pg::log_z(T, order) - B.log_z, 1e-10);

    // The derivatives, two parameters alternating over the plaquettes.
    std::vector<std::size_t> of(T.plaquettes.size());
    for (std::size_t f = 0; f < of.size(); f++) of[f] = f % 2;
    const std::array<double, 2> par{0.35, -0.2};
    pg::Z2Theory U = T;
    for (std::size_t f = 0; f < of.size(); f++) U.beta[f] = par[of[f]];
    const auto D = pg::log_z_derivatives<2>(U, order, of, par);
    const pg::Contraction CU = pg::contract(U, order);
    std::array<double, 2> sums{0, 0};
    for (std::size_t f = 0; f < of.size(); f++) sums[of[f]] += CU.plaquette[f];
    s.check("calx jets: log Z of the jet pass", D.log_z - CU.log_z, 1e-10);
    s.check("calx jets: ∂ log Z/∂β_p = Σ ⟨h_f⟩ over p's plaquettes", std::max(std::fabs(D.gradient[0] - sums[0]), std::fabs(D.gradient[1] - sums[1])), 1e-10);
    double hw = 0;
    const double h = 1e-5;
    for (std::size_t q = 0; q < 2; q++) {
      auto up = par, down = par;
      up[q] += h;
      down[q] -= h;
      const auto gu = pg::log_z_derivatives<2, 1>(U, order, of, up).gradient, gd = pg::log_z_derivatives<2, 1>(U, order, of, down).gradient;
      for (std::size_t p = 0; p < 2; p++) hw = std::max(hw, std::fabs(D.hessian[p][q] - (gu[p] - gd[p]) / (2 * h)));
    }
    s.check("calx jets: the Hessian of log Z against differences of the gradient", hw, 1e-7);

    // Metropolis with these couplings and frozen links.
    pg::GaugeField<pg::ZN> field(pg::ZN{2}, X);
    field.set_couplings(beta);
    for (std::size_t e = 0; e < X.edges.size(); e++)
      if (fixed[e]) field.freeze(e, fixed[e] < 0 ? 1u : 0u);
    calx::random::Stream z(20261024, 0, 0, 0);
    std::vector<double> mean_h(X.faces.size(), 0);
    const std::size_t sweeps = 40000;
    for (std::size_t k = 0; k < 2000; k++) pg::metropolis_sweep(field, 1.0, z);
    for (std::size_t k = 0; k < sweeps; k++) {
      pg::metropolis_sweep(field, 1.0, z);
      for (std::size_t f = 0; f < X.faces.size(); f++) mean_h[f] += field.plaquette(f);
    }
    for (auto& x : mean_h) x /= static_cast<double>(sweeps);
    bool still = true;
    for (std::size_t e = 0; e < X.edges.size(); e++)
      if (fixed[e]) still = still && field.link(e) == (fixed[e] < 0 ? 1u : 0u);
    s.require("Metropolis leaves the frozen links where they were frozen", still);
    s.check("Metropolis with per-plaquette couplings against the contraction, ⟨h_f⟩", max_difference(mean_h, C.plaquette), 0.03);

    // The heat bath of a finite group, with the same couplings and frozen links.
    std::fill(mean_h.begin(), mean_h.end(), 0.0);
    for (std::size_t k = 0; k < 2000; k++) pg::heat_bath_sweep(field, 1.0, z);
    for (std::size_t k = 0; k < sweeps; k++) {
      pg::heat_bath_sweep(field, 1.0, z);
      for (std::size_t f = 0; f < X.faces.size(); f++) mean_h[f] += field.plaquette(f);
    }
    for (auto& x : mean_h) x /= static_cast<double>(sweeps);
    still = true;
    for (std::size_t e = 0; e < X.edges.size(); e++)
      if (fixed[e]) still = still && field.link(e) == (fixed[e] < 0 ? 1u : 0u);
    s.require("the heat bath leaves the frozen links where they were frozen", still);
    s.check("the heat bath with per-plaquette couplings against the contraction, ⟨h_f⟩", max_difference(mean_h, C.plaquette), 0.03);
  }

  // A circle of four links in no 2-cell: its Wilson loop vanishes. The heat bath finds it; Metropolis
  // on ℤ₂, whose one proposal is always accepted there, flips every link at every sweep, and a loop
  // of even length never changes.
  {
    mt::CWComplex X;
    X.vertices = 4;
    X.edges = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    const mt::Path loop{{0, 1}, {1, 1}, {2, 1}, {3, 1}};
    pg::GaugeField<pg::ZN> heat(pg::ZN{2}, X), metropolis(pg::ZN{2}, X);
    calx::random::Stream z(20261024, 1, 0, 0);
    double wh = 0, wm = 0;
    const std::size_t sweeps = 20000;
    for (std::size_t k = 0; k < sweeps; k++) {
      pg::heat_bath_sweep(heat, 1.0, z);
      pg::metropolis_sweep(metropolis, 1.0, z);
      wh += heat.wilson_loop(loop);
      wm += metropolis.wilson_loop(loop);
    }
    s.check("a circle in no 2-cell: the heat bath's Wilson loop vanishes", wh / static_cast<double>(sweeps), 0.03);
    s.require("... where Metropolis on ℤ₂ never changes it (periodic)", std::fabs(wm) == static_cast<double>(sweeps));
  }

  // Closed surfaces, one coupling: Migdal and Witten.
  for (const auto& [name, X, chi] : std::vector<std::tuple<std::string, mt::CWComplex, long long>>{
           {"the 4 × 4 torus", mt::torus_lattice({4, 4}), 0}, {"the surface of genus 2", mt::orientable_surface(2), -2}}) {
    const double beta = 0.7;
    const pg::Z2Theory T = pg::Z2Theory::of(X, std::vector<double>(X.faces.size(), beta));
    const auto order = natural(T.plaquettes.size());
    const pg::Contraction C = pg::contract(T, order);
    const auto exact = pg::migdal_witten(pg::expansion_zn(2, beta), chi, X.faces.size());
    // Haar measures of mass 1: log Z less E log 2.
    s.check(name + ": log Z (Haar) against Migdal and Witten", C.log_z - static_cast<double>(X.edges.size()) * std::log(2.0) - exact.log_z, 1e-10);
    const double plaquette = std::accumulate(C.plaquette.begin(), C.plaquette.end(), 0.0) / static_cast<double>(C.plaquette.size());
    s.check(name + ": the mean plaquette against Migdal and Witten", plaquette - exact.mean_plaquette, 1e-10);
  }

  // Strong coupling: log Z − β F → log #flat.
  {
    const mt::CWComplex X = mt::torus_lattice({3, 4});
    const double beta = 40;
    const pg::Z2Theory T = pg::Z2Theory::of(X, std::vector<double>(X.faces.size(), beta));
    const double flat = std::stod(pg::count_flat_connections(X, morphe::algebra::cyclic(2)).str());
    s.check("a 3 × 4 torus at β = 40: log Z − β F is log of the flat connections counted through π₁", pg::log_z(T, natural(X.faces.size())) - beta * static_cast<double>(X.faces.size()) - std::log(flat), 1e-9,
            std::to_string(flat) + " flat connections");
  }

  // Gauge invariance at a cut: the spin networks are ℤ/2 1-cycles.
  {
    const mt::CWComplex X = mt::torus_lattice({4, 4});
    std::vector<double> beta(X.faces.size());
    for (std::size_t f = 0; f < beta.size(); f++) beta[f] = 0.2 + 0.05 * static_cast<double>(f % 5);
    const pg::Z2Theory T = pg::Z2Theory::of(X, beta);
    std::vector<pg::ContractionCut> cuts;
    for (std::size_t k = 1; k < X.faces.size(); k += 3) cuts.push_back({k, {}, {}});
    const pg::Contraction C = pg::contract(T, natural(X.faces.size()), {}, cuts, 1u << 20);
    bool cycles = true;
    std::size_t count = 0;
    for (const auto& law : C.cuts)
      for (const auto& [set, w] : law.networks) {
        std::map<std::size_t, int> degree;
        for (std::size_t e : set) {
          degree[X.edges[e][0]] ^= 1;
          degree[X.edges[e][1]] ^= 1;
        }
        for (const auto& [v, d] : degree) cycles = cycles && d == 0;
        count++;
      }
    s.require("with no link fixed, every spin network crossing a cut is a ℤ/2 1-cycle of the cut", cycles, std::to_string(count) + " networks at " + std::to_string(cuts.size()) + " cuts");
  }
  return s.done();
}
