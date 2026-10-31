// Checks of the master equation (master_equation.hpp) against closed forms and its theorems.
//
// - The generator's rows sum to zero.
// - Birth–death chains: the stationary law against the closed form of detailed balance,
//   π_i/π_0 = Π_{j<i} λ_j/μ_{j+1}, entry by entry to a relative 10⁻¹², over forty-six orders of magnitude
//   (GTH forms no differences); π Q = 0 on random chains; a reducible chain refused.
// - A driven ring of n states (k₊ one way, k₋ the other): one cycle, A = n ln(k₊/k₋), J = (k₊ − k₋)/n,
//   σ = (k₊ − k₋) ln(k₊/k₋), for n = 3 and 7.
// - Schnakenberg's identity σ = Σ_C J(C) A(C) on random graphs, every edge's term J ln(π k/π k) ≥ 0, and
//   one cycle for each chord (|E| − |V| + 1).
// - Detailed balance: rates from energies, k_ij = s_ij e^{(E_i − E_j)/2} with s symmetric, give the
//   Boltzmann law, σ = 0 and every affinity 0; one rate changed breaks Kolmogorov's criterion.
// - The H-theorem: D(p(t) ‖ π) non-increasing on random generators wherever it exceeds 10⁻¹⁰ nats, and
//   staying below once there (near 0, D is at the rounding of exp(Qt) against GTH's π, some 10⁻¹⁴);
//   falling towards 0.
// - evolve: against explicit Euler, which approaches it at order 1; and against the spectral
//   decomposition of a reversible generator, Q = Π^{−1/2} S Π^{1/2} with S symmetric (calx's
//   symmetric eigenproblem).
// - The relative entropy: zero on itself, a two-point closed form, infinite off the support.

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include <calx/linalg/eigen.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <physis/statistical_mechanics/master_equation.hpp>

using calx::linalg::Matrix;
using calx::verification::Suite;
namespace sm = physis::statistical_mechanics;

namespace {

// Random rates on a random connected graph, both directions positive on each edge, log-uniform.
Matrix random_rates(calx::random::Stream& rng, std::size_t n, double p) {
  Matrix k(n);
  const auto rate = [&] { return std::exp(rng.uniform(-2, 2)); };
  for (std::size_t i = 0; i + 1 < n; i++) k(i, i + 1) = rate(), k(i + 1, i) = rate();
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = i + 2; j < n; j++)
      if (rng.uniform() < p) k(i, j) = rate(), k(j, i) = rate();
  return k;
}

Matrix ring(std::size_t n, double kp, double km) {
  Matrix k(n);
  for (std::size_t i = 0; i < n; i++) k(i, (i + 1) % n) = kp, k((i + 1) % n, i) = km;
  return k;
}

}  // namespace

int main() {
  Suite s;
  calx::random::Stream rng(20261030, 18, 0, 0);

  // ---- the generator and the stationary law ---------------------------------------------------------------------------------
  {
    const Matrix k = random_rates(rng, 9, 0.4);
    const Matrix Q = sm::generator(k);
    double rows = 0;
    for (std::size_t i = 0; i < 9; i++) {
      double r = 0;
      for (std::size_t j = 0; j < 9; j++) r += Q(i, j);
      rows = std::max(rows, std::fabs(r));
    }
    s.check("the generator's rows sum to zero", rows, 1e-14);

    // A birth–death chain spanning sixty orders of magnitude.
    const std::size_t n = 30;
    Matrix bd(n);
    for (std::size_t i = 0; i + 1 < n; i++) bd(i, i + 1) = 0.5 + 0.1 * i, bd(i + 1, i) = 50.0 + i;
    const auto pi = sm::stationary(bd);
    std::vector<double> logw(n, 0.0);
    for (std::size_t i = 1; i < n; i++) logw[i] = logw[i - 1] + std::log(bd(i - 1, i) / bd(i, i - 1));
    double norm = 0;
    for (double w : logw) norm += std::exp(w);
    double worst = 0;
    for (std::size_t i = 0; i < n; i++) worst = std::max(worst, std::fabs(pi[i] / (std::exp(logw[i]) / norm) - 1));
    s.check("birth–death: π_i ∝ Π λ_j/μ_{j+1} entry by entry, π from 1 to 10⁻⁴⁶ (relative)", worst, 1e-12);

    double balance = 0;
    for (int trial = 0; trial < 20; trial++) {
      const Matrix kr = random_rates(rng, 4 + trial % 8, 0.5);
      const auto p = sm::stationary(kr);
      const Matrix Qr = sm::generator(kr);
      for (std::size_t j = 0; j < kr.n; j++) {
        double r = 0;
        for (std::size_t i = 0; i < kr.n; i++) r += p[i] * Qr(i, j);
        balance = std::max(balance, std::fabs(r));
      }
      balance = std::max(balance, std::fabs(std::accumulate(p.begin(), p.end(), 0.0) - 1));
    }
    s.check("random chains: π Q = 0 and Σ π = 1", balance, 1e-14);

    Matrix reducible(3);
    reducible(0, 1) = reducible(1, 0) = 1;
    reducible(2, 0) = 1;
    bool refused = false;
    try {
      sm::stationary(reducible);
    } catch (const std::invalid_argument&) {
      refused = true;
    }
    s.require("a reducible chain is refused", refused && !sm::irreducible(reducible));
  }

  // ---- driven rings -----------------------------------------------------------------------------------------------------------
  for (std::size_t n : {3, 7}) {
    const double kp = 2.0, km = 0.5;
    const auto I = sm::irreversibility(ring(n, kp, km));
    const std::string name = "a driven ring of " + std::to_string(n) + " states: ";
    s.require(name + "one cycle, of all its states", I.cycles.size() == 1 && I.cycles[0].states.size() == n);
    const auto& C = I.cycles[0];
    s.check(name + "|A| = n ln(k₊/k₋), |J| = (k₊ − k₋)/n, of one sign",
            std::max({std::fabs(std::fabs(C.affinity) - n * std::log(kp / km)), std::fabs(std::fabs(C.current) - (kp - km) / n), C.affinity * C.current < 0 ? 1.0 : 0.0}), 1e-14);
    s.check(name + "σ = (k₊ − k₋) ln(k₊/k₋) = J A", std::max(std::fabs(I.entropy_production - (kp - km) * std::log(kp / km)), std::fabs(I.entropy_production - C.current * C.affinity)), 1e-14);
  }

  // ---- Schnakenberg's identity on random graphs --------------------------------------------------------------------------------------
  {
    double identity = 0;
    bool positive = true, basis = true;
    for (int trial = 0; trial < 60; trial++) {
      const std::size_t n = 4 + trial % 9;
      const Matrix k = random_rates(rng, n, 0.45);
      const auto I = sm::irreversibility(k);
      double sum = 0;
      for (const auto& C : I.cycles) sum += C.current * C.affinity;
      identity = std::max(identity, std::fabs(I.entropy_production - sum) / std::max(I.entropy_production, 1e-300));
      std::size_t edges = 0;
      for (std::size_t i = 0; i < n; i++)
        for (std::size_t j = i + 1; j < n; j++)
          if (k(i, j) > 0) {
            edges++;
            const double J = I.current(i, j), X = std::log(I.pi[i] * k(i, j)) - std::log(I.pi[j] * k(j, i));
            positive = positive && J * X >= 0;
          }
      // A tree has no cycle, and every chain on a tree is in detailed balance.
      const bool tree = I.cycles.empty();
      basis = basis && I.cycles.size() == edges - n + 1 && (tree ? I.entropy_production < 1e-15 && I.detailed_balance : I.entropy_production > 0 && !I.detailed_balance);
    }
    s.check("Schnakenberg: σ = Σ_C J(C) A(C) over a fundamental cycle basis, 60 random graphs (relative)", identity, 1e-12);
    s.require("every edge's J ln(π_i k_ij / π_j k_ji) ≥ 0, so σ ≥ 0", positive);
    s.require("one cycle for each chord, |E| − |V| + 1; σ > 0 for random rates on a graph with cycles, 0 on a tree", basis);
  }

  // ---- detailed balance --------------------------------------------------------------------------------------------------------
  {
    const std::size_t n = 6;
    std::vector<double> E(n);
    for (double& e : E) e = rng.uniform(-3, 3);
    Matrix k(n);
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = i + 1; j < n; j++) {
        const double sij = std::exp(rng.uniform(-1, 1));
        k(i, j) = sij * std::exp((E[i] - E[j]) / 2);
        k(j, i) = sij * std::exp((E[j] - E[i]) / 2);
      }
    const auto I = sm::irreversibility(k);
    double Z = 0;
    for (double e : E) Z += std::exp(-e);
    double boltzmann = 0, affinity = 0;
    for (std::size_t i = 0; i < n; i++) boltzmann = std::max(boltzmann, std::fabs(I.pi[i] / (std::exp(-E[i]) / Z) - 1));
    for (const auto& C : I.cycles) affinity = std::max(affinity, std::fabs(C.affinity));
    s.check("rates from energies: the Boltzmann law (relative)", boltzmann, 1e-13);
    s.check("and σ = 0, every affinity 0", std::max(I.entropy_production, affinity), 1e-13);
    s.require("Kolmogorov's criterion holds; ten cycles in the complete graph of six", I.detailed_balance && I.cycles.size() == 10);
    k(1, 4) *= 1.5;
    const auto broken = sm::irreversibility(k);
    s.require("one rate changed: the criterion fails and σ > 0", !broken.detailed_balance && broken.entropy_production > 1e-6);
  }

  // ---- the H-theorem ---------------------------------------------------------------------------------------------------------------
  {
    double rise = 0;
    double last = 0, first = 0;
    for (int trial = 0; trial < 5; trial++) {
      const Matrix k = random_rates(rng, 6, 0.6);
      const auto pi = sm::stationary(k);
      std::vector<double> p0(6);
      double total = 0;
      for (double& x : p0) total += x = rng.uniform();
      for (double& x : p0) x /= total;
      double prev = sm::relative_entropy(p0, pi);
      if (trial == 0) first = prev;
      for (int step = 1; step <= 100; step++) {
        const double D = sm::relative_entropy(sm::evolve(k, p0, 0.05 * step), pi);
        // Where D is meaningful (> 10⁻¹⁰), it must not rise at all; below, it must stay below.
        rise = std::max(rise, prev > 1e-10 ? D - prev : D - 1e-10);
        prev = D;
      }
      if (trial == 0) last = prev;
    }
    s.check("the H-theorem: D(p(t) ‖ π) non-increasing on five random generators (its largest rise above 10⁻¹⁰ nats)", std::max(rise, 0.0), 1e-16, "nats");
    s.require("and falling towards 0", last < 1e-3 * first, std::to_string(first) + " → " + std::to_string(last));
  }

  // ---- evolve ----------------------------------------------------------------------------------------------------------------------
  {
    const Matrix k = random_rates(rng, 7, 0.5);
    const Matrix Q = sm::generator(k);
    std::vector<double> p0(7, 0.0);
    p0[2] = 1;
    const auto p = sm::evolve(k, p0, 1.0);
    const auto euler = [&](int steps) {
      std::vector<double> q = p0, dq(7);
      const double h = 1.0 / steps;
      for (int m = 0; m < steps; m++) {
        for (std::size_t j = 0; j < 7; j++) {
          dq[j] = 0;
          for (std::size_t i = 0; i < 7; i++) dq[j] += q[i] * Q(i, j);
        }
        for (std::size_t j = 0; j < 7; j++) q[j] += h * dq[j];
      }
      double e = 0;
      for (std::size_t j = 0; j < 7; j++) e = std::max(e, std::fabs(q[j] - p[j]));
      return e;
    };
    const double e1 = euler(4000), e2 = euler(8000);
    s.check("evolve: explicit Euler approaches p0 exp(Qt) at order 1", calx::verification::order(e1, e2) - 1, 0.02);
    s.check("evolve: and is within 10⁻³ of it at h = 1/8000", e2, 1e-3);

    // Reversible: rates from energies, exp(Qt) from the symmetric S = Π^{1/2} Q Π^{−1/2}.
    const std::size_t n = 6;
    std::vector<double> E(n);
    for (double& e : E) e = rng.uniform(-2, 2);
    Matrix kr(n);
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = i + 1; j < n; j++)
        if (rng.uniform() < 0.7 || j == i + 1) {
          const double sij = std::exp(rng.uniform(-1, 1));
          kr(i, j) = sij * std::exp((E[i] - E[j]) / 2), kr(j, i) = sij * std::exp((E[j] - E[i]) / 2);
        }
    const auto pi = sm::stationary(kr);
    const Matrix Qr = sm::generator(kr);
    Matrix S(n);
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++) S(i, j) = std::sqrt(pi[i]) * Qr(i, j) / std::sqrt(pi[j]);
    const auto eig = calx::linalg::symmetric_eigen(S);
    const double t = 0.7;
    std::vector<double> q0(n);
    double total = 0;
    for (double& x : q0) total += x = rng.uniform();
    for (double& x : q0) x /= total;
    // p(t)_j = Σ_i q0_i [Π^{−1/2} V e^{Λt} Vᵀ Π^{1/2}]_ij
    std::vector<double> spectral(n, 0.0);
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++) {
        double e = 0;
        for (std::size_t m = 0; m < n; m++) e += eig.vectors(i, m) * std::exp(eig.values[m] * t) * eig.vectors(j, m);
        spectral[j] += q0[i] * e * std::sqrt(pi[j]) / std::sqrt(pi[i]);
      }
    const auto pr = sm::evolve(kr, q0, t);
    double diff = 0, mass = 0;
    for (std::size_t j = 0; j < n; j++) diff = std::max(diff, std::fabs(pr[j] - spectral[j])), mass += pr[j];
    s.check("evolve: against the spectral decomposition of a reversible generator", diff, 1e-13);
    s.check("evolve: probability kept", mass - 1, 1e-14);
  }

  // ---- the relative entropy ----------------------------------------------------------------------------------------------------------
  {
    const std::vector<double> p{0.2, 0.8}, q{0.5, 0.5}, r{0.0, 1.0}, z{0.3, 0.0, 0.7}, w{0.3, 0.2, 0.5};
    const double closed = 0.2 * std::log(0.4) + 0.8 * std::log(1.6);
    s.check("D(p ‖ q) for two points against its closed form; D(p ‖ p) = 0", std::max(std::fabs(sm::relative_entropy(p, q) - closed), std::fabs(sm::relative_entropy(p, p))), 1e-16, "nats");
    s.require("infinite when p charges a point q does not, finite the other way",
              std::isinf(sm::relative_entropy(p, r)) && std::isfinite(sm::relative_entropy(r, p)) && std::isinf(sm::relative_entropy(w, z)) && std::isfinite(sm::relative_entropy(z, w)));
  }

  return s.done();
}
