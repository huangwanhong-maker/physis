// Checks of the entropy of mixing (mixing.hpp) against closed forms and the theorems of information.
//
// - A hand-computed table, [[2, 0], [1, 1]]: H(Q) = −(¾ ln ¾ + ¼ ln ¼), H(Q|X) = ½ ln 2, M their difference.
// - Proportional counts (every place in the whole's proportions): M = 0 and every local term 0.
// - Complete separation (each place one kind): H(Q|X) = 0, M = H(Q), normalised 1.
// - On random tables: H(Q) − H(Q|X) = M = Σ local, 0 ≤ M ≤ H(Q).
// - The chain rule over random partitions of the places: M = M_between + Σ_g p_g M_g exactly.
// - Invariance under permuting places and kinds.
// - The data-processing inequality: places exchanging their compositions by a random doubly stochastic
//   matrix (a mixture of permutations, Birkhoff) never raise M.

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <physis/statistical_mechanics/mixing.hpp>

using calx::linalg::Rectangular;
using calx::verification::Suite;
namespace sm = physis::statistical_mechanics;

namespace {

Rectangular random_counts(calx::random::Stream& rng, std::size_t places, std::size_t kinds) {
  Rectangular c(places, kinds);
  for (double& x : c.a) x = static_cast<double>(rng.below(20)) * (rng.uniform() < 0.8 ? 1 : 0);
  c(0, 0) += 1;
  return c;
}

std::vector<std::size_t> shuffled(calx::random::Stream& rng, std::size_t n) {
  std::vector<std::size_t> p(n);
  std::iota(p.begin(), p.end(), 0);
  for (std::size_t i = n; i-- > 1;) std::swap(p[i], p[rng.below(i + 1)]);
  return p;
}

}  // namespace

int main() {
  Suite s;
  calx::random::Stream rng(20261030, 20, 0, 0);

  // ---- closed forms -----------------------------------------------------------------------------------------------------------
  {
    const auto m = sm::mixing(Rectangular::from_rows({{2, 0}, {1, 1}}));
    const double H = -(0.75 * std::log(0.75) + 0.25 * std::log(0.25)), Hc = 0.5 * std::log(2.0);
    s.check("[[2, 0], [1, 1]]: H(Q), H(Q|X) and M = H(Q) − H(Q|X) in closed form",
            std::max({std::fabs(m.entropy - H), std::fabs(m.conditional - Hc), std::fabs(m.departure - (H - Hc))}), 1e-15, "nats");

    const auto p = sm::mixing(Rectangular::from_rows({{1, 2, 3}, {2, 4, 6}, {10, 20, 30}, {7, 14, 21}}));
    s.check("proportional counts: M = 0 and every local term 0", std::max(p.departure, *std::max_element(p.local.begin(), p.local.end())), 1e-15, "nats");

    const auto sep = sm::mixing(Rectangular::from_rows({{5, 0, 0}, {0, 3, 0}, {0, 0, 9}, {4, 0, 0}}));
    s.check("complete separation: H(Q|X) = 0, M = H(Q), normalised 1",
            std::max({std::fabs(sep.conditional), std::fabs(sep.departure - sep.entropy), std::fabs(sep.normalised - 1)}), 1e-15);
  }

  // ---- random tables: identities, the chain rule, invariance ---------------------------------------------------------------------------------
  {
    double identity = 0, chain = 0, invariance = 0;
    bool bounded = true;
    for (int trial = 0; trial < 200; trial++) {
      const std::size_t P = 3 + trial % 9, K = 2 + trial % 5;
      const auto c = random_counts(rng, P, K);
      const auto m = sm::mixing(c);
      const double local = std::accumulate(m.local.begin(), m.local.end(), 0.0);
      identity = std::max({identity, std::fabs(m.entropy - m.conditional - m.departure), std::fabs(local - m.departure)});
      bounded = bounded && m.departure >= 0 && m.departure <= m.entropy + 1e-15 && m.normalised >= 0 && m.normalised <= 1;

      std::vector<std::size_t> group(P);
      for (auto& g : group) g = rng.below(1 + P / 2);
      const auto d = sm::mixing_decomposition(c, group);
      chain = std::max(chain, std::fabs(d.total - d.between - d.within));

      const auto rp = shuffled(rng, P), kp = shuffled(rng, K);
      Rectangular perm(P, K);
      for (std::size_t x = 0; x < P; x++)
        for (std::size_t q = 0; q < K; q++) perm(rp[x], kp[q]) = c(x, q);
      invariance = std::max(invariance, std::fabs(sm::mixing(perm).departure - m.departure));
    }
    s.check("random tables: H(Q) − H(Q|X) = M = Σ_x p_x D(p(·|x) ‖ p(·))", identity, 1e-14, "nats");
    s.require("random tables: 0 ≤ M ≤ H(Q)", bounded);
    s.check("the chain rule: M = M_between + Σ_g p_g M_g over random partitions of the places", chain, 1e-14, "nats");
    s.check("M is invariant under permuting places and kinds", invariance, 1e-14, "nats");
  }

  // ---- the data-processing inequality ---------------------------------------------------------------------------------------------------------
  {
    double rise = 0;
    for (int trial = 0; trial < 300; trial++) {
      const std::size_t P = 4 + trial % 6, K = 3;
      const auto c = random_counts(rng, P, K);
      // W = Σ_r w_r Π_r, a mixture of permutation matrices: doubly stochastic.
      Rectangular W(P, P);
      double total = 0;
      std::vector<double> weight(4);
      for (double& w : weight) total += w = rng.uniform();
      for (double w : weight) {
        const auto pi = shuffled(rng, P);
        for (std::size_t x = 0; x < P; x++) W(pi[x], x) += w / total;
      }
      const Rectangular exchanged = W * c;  // place x′ receives Σ_x W(x′, x) of place x's composition
      rise = std::max(rise, sm::mixing(exchanged).departure - sm::mixing(c).departure);
    }
    s.check("data processing: a doubly stochastic exchange of compositions never raises M (largest rise)", std::max(rise, 0.0), 1e-15, "nats");
  }

  return s.done();
}
