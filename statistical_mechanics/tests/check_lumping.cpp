// Checks of lumping (lumping.hpp).
//
// - Two birth–death components, independent: lumping onto the first is exact, its lumped rates the first
//   component's own, the generators intertwining (defect 0) and so do the semigroups, exp(Qt)(f∘C) =
//   (exp(Q̄t)f)∘C, to 10⁻¹² for random f; made to depend on the second component, the first is no longer
//   lumpable and the test names the discrepancy, the difference of the two rates.
// - Mechanisms: each of two non-lumpable mechanisms whose failures cancel; their sum lumpable.
// - The aggregated process: on the identity partition the coarse density is the fine path density written
//   out by hand; under lumpability it equals the lumped chain's for every initial law with the same block
//   masses; without, two initial laws with the same block masses give different densities of one coarse
//   path — the coarse description has memory — and the dwell in a block is not exponential (squared
//   coefficient of variation 1 when lumpable, above 1 when not).
// - Approximate lumping: the aggregation by any weights of a lumpable chain is its lumped chain, discrepancy 0;
//   of one that is not, by its stationary law, keeps the blocks' stationary law; and on random chains the
//   blocks' law at T under the fine chain and under the aggregate differ by at most 1 − e^{−δT}.

#include <cmath>
#include <string>
#include <vector>

#include <calx/linalg/expm.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <physis/statistical_mechanics/lumping.hpp>

using calx::linalg::Matrix;
using calx::verification::Suite;
namespace sm = physis::statistical_mechanics;

namespace {

// States (i, j), i ∈ 0..n1−1, j ∈ 0..n2−1, index i * n2 + j. The first component moves i → i ± 1 at rates
// up(i, j), down(i, j); the second j → j ± 1 at u2, d2.
template <class Up, class Down>
Matrix product(std::size_t n1, std::size_t n2, Up up, Down down, double u2, double d2) {
  Matrix k(n1 * n2);
  for (std::size_t i = 0; i < n1; i++)
    for (std::size_t j = 0; j < n2; j++) {
      const std::size_t s = i * n2 + j;
      if (i + 1 < n1) k(s, s + n2) = up(i, j);
      if (i > 0) k(s, s - n2) = down(i, j);
      if (j + 1 < n2) k(s, s + 1) = u2;
      if (j > 0) k(s, s - 1) = d2;
    }
  return k;
}

double lumped_log_density(const Matrix& kb, const std::vector<double>& pb, const sm::CoarsePath& p) {
  double lp = std::log(pb[p.blocks[0]]);
  for (std::size_t s = 0; s < p.blocks.size(); s++) {
    double out = 0;
    for (std::size_t c = 0; c < kb.n; c++)
      if (c != p.blocks[s]) out += kb(p.blocks[s], c);
    lp -= out * p.dwell[s];
    if (s + 1 < p.blocks.size()) lp += std::log(kb(p.blocks[s], p.blocks[s + 1]));
  }
  return lp;
}

}  // namespace

int main() {
  Suite s;
  calx::random::Stream rng(17, 1, 2, 3);
  const std::size_t n1 = 4, n2 = 3;
  std::vector<std::size_t> first(n1 * n2);
  for (std::size_t x = 0; x < first.size(); x++) first[x] = x / n2;

  const auto k = product(n1, n2, [](std::size_t i, std::size_t) { return 1.0 + 0.5 * i; }, [](std::size_t i, std::size_t) { return 0.7 * i; }, 0.9, 1.3);
  const auto L = sm::strongly_lumpable(k, first);
  s.require("two independent components: the first lumps", L.lumpable);
  const auto kb = sm::lumped(k, first);
  double own = 0;
  for (std::size_t i = 0; i + 1 < n1; i++) own = std::max(own, std::fabs(kb(i, i + 1) - (1.0 + 0.5 * i)) + std::fabs(kb(i + 1, i) - 0.7 * (i + 1)));
  s.check("its lumped rates are the first component's own", own, 1e-14);
  s.check("the generators intertwine", sm::intertwining_defect(k, first, kb), 1e-14);
  {
    const auto E = calx::linalg::expm(0.7 * sm::generator(k)), Eb = calx::linalg::expm(0.7 * sm::generator(kb));
    double worst = 0;
    for (int trial = 0; trial < 5; trial++) {
      std::vector<double> f(n1);
      for (auto& x : f) x = rng.normal();
      for (std::size_t i = 0; i < k.n; i++) {
        double fine = 0, coarse = 0;
        for (std::size_t j = 0; j < k.n; j++) fine += E(i, j) * f[first[j]];
        for (std::size_t b = 0; b < n1; b++) coarse += Eb(first[i], b) * f[b];
        worst = std::max(worst, std::fabs(fine - coarse));
      }
    }
    s.check("and so do the semigroups: exp(Qt)(f∘C) = (exp(Q̄t) f)∘C", worst, 1e-12);
  }
  const auto kc = product(n1, n2, [](std::size_t i, std::size_t j) { return 1.0 + 0.5 * i + 0.4 * j; }, [](std::size_t i, std::size_t) { return 0.7 * i; }, 0.9, 1.3);
  const auto Lc = sm::strongly_lumpable(kc, first);
  s.require("made to depend on the second component, the first does not lump", !Lc.lumpable);
  s.check("the discrepancy is the difference of the rates, 0.4 × 2", Lc.discrepancy - 0.8, 1e-12);

  {
    Matrix m1(3), m2(3);
    m1(0, 2) = 1, m2(1, 2) = 1;
    const std::vector<std::size_t> part{0, 0, 1};
    s.require("neither mechanism lumps on its own", !sm::strongly_lumpable(m1, part).lumpable && !sm::strongly_lumpable(m2, part).lumpable && !sm::eventwise_lumpable({m1, m2}, part).lumpable);
    s.require("their sum does: the failures cancel (lumping each suffices, and is not necessary)", sm::strongly_lumpable(m1 + m2, part).lumpable);
  }

  // ---- the aggregated process ----
  {
    // The identity partition: the fine path density.
    std::vector<std::size_t> id(k.n);
    for (std::size_t i = 0; i < k.n; i++) id[i] = i;
    std::vector<double> p0(k.n, 0.0);
    p0[4] = 1;
    const sm::CoarsePath path{{4, 7, 6, 3}, {0.3, 0.8, 0.2, 0.5}};
    double by_hand = 0;
    for (std::size_t q = 0; q < path.blocks.size(); q++) {
      double out = 0;
      for (std::size_t j = 0; j < k.n; j++)
        if (j != path.blocks[q]) out += k(path.blocks[q], j);
      by_hand -= out * path.dwell[q];
      if (q + 1 < path.blocks.size()) by_hand += std::log(k(path.blocks[q], path.blocks[q + 1]));
    }
    s.check("on the identity partition, the fine path's log-density", sm::coarse_log_density(k, id, p0, path) - by_hand, 1e-12);
  }
  {
    const sm::CoarsePath path{{1, 2, 1, 0, 1}, {0.4, 0.25, 1.1, 0.3, 0.6}};
    std::vector<double> pa(k.n, 0.0), pc(k.n, 0.0);
    pa[1 * n2 + 0] = 1;                                    // block 1, second component at 0
    pc[1 * n2 + 2] = 0.5, pc[1 * n2 + 1] = 0.5;            // block 1, elsewhere
    std::vector<double> pb(n1, 0.0);
    pb[1] = 1;
    const double want = lumped_log_density(kb, pb, path);
    s.check("lumpable: the coarse density is the lumped chain's", sm::coarse_log_density(k, first, pa, path) - want, 1e-11);
    s.check("for every initial law with the same block masses", sm::coarse_log_density(k, first, pc, path) - want, 1e-11);
    const double a = sm::coarse_log_density(kc, first, pa, path), c = sm::coarse_log_density(kc, first, pc, path);
    s.require("not lumpable: the same coarse path, two initial laws of the same block masses, two densities (memory)", std::fabs(a - c) > 1e-3);
    const auto cv2 = [](const sm::PhaseType& P) { const double m1 = P.moment(1), m2 = P.moment(2); return m2 / (m1 * m1) - 1; };
    std::vector<double> alpha(n2, 1.0 / n2);
    s.check("lumpable: the dwell in a block is exponential (CV² = 1)", cv2(sm::dwell_law(k, first, 1, alpha)) - 1, 1e-10);
    s.require("not lumpable: it is not (CV² ≠ 1)", std::fabs(cv2(sm::dwell_law(kc, first, 1, alpha)) - 1) > 1e-3);
  }
  {
    std::vector<double> w(k.n);
    for (auto& x : w) x = 0.1 + rng.uniform();
    const auto ag = sm::aggregated(k, first, w);
    s.check("a lumpable chain aggregated by any weights is its lumped chain, discrepancy 0", sm::outward_discrepancy(k, first, ag) + sm::intertwining_defect(k, first, ag), 1e-13);
    const auto pi = sm::stationary(kc);
    const auto agc = sm::aggregated(kc, first, pi);
    const auto pib = sm::stationary(agc);
    std::vector<double> mass(n1, 0.0);
    for (std::size_t i = 0; i < kc.n; i++) mass[first[i]] += pi[i];
    double worst = 0;
    for (std::size_t b = 0; b < n1; b++) worst = std::max(worst, std::fabs(pib[b] - mass[b]));
    s.check("aggregated by its stationary law, a chain that does not lump keeps its blocks' stationary law", worst, 1e-12);
    bool bounded = true;
    double tightest = 0;
    for (int trial = 0; trial < 20; trial++) {
      Matrix r(6);
      for (std::size_t i = 0; i < 6; i++)
        for (std::size_t j = 0; j < 6; j++)
          if (i != j) r(i, j) = rng.uniform() < 0.6 ? 0.2 + 2 * rng.uniform() : 0.0;
      const std::vector<std::size_t> part{0, 0, 1, 1, 1, 2};
      std::vector<double> wt(6);
      for (auto& x : wt) x = 0.1 + rng.uniform();
      const auto co = sm::aggregated(r, part, wt);
      const double delta = sm::outward_discrepancy(r, part, co);
      std::vector<double> p0(6);
      double z = 0;
      for (auto& x : p0) z += x = rng.uniform();
      for (auto& x : p0) x /= z;
      std::vector<double> q0(3, 0.0);
      for (std::size_t i = 0; i < 6; i++) q0[part[i]] += p0[i];
      for (double T : {0.3, 1.0, 2.5}) {
        const auto pt = sm::evolve(r, p0, T), qt = sm::evolve(co, q0, T);
        std::vector<double> blocks(3, 0.0);
        for (std::size_t i = 0; i < 6; i++) blocks[part[i]] += pt[i];
        double tv = 0;
        for (std::size_t b = 0; b < 3; b++) tv += 0.5 * std::fabs(blocks[b] - qt[b]);
        const double bound = 1 - std::exp(-delta * T);
        bounded = bounded && tv <= bound + 1e-12;
        if (bound > 0) tightest = std::max(tightest, tv / bound);
      }
    }
    s.require("on 20 random chains the blocks' law differs by at most 1 − e^{−δT} at T = 0.3, 1, 2.5", bounded, "largest share of the bound used: " + std::to_string(tightest));
  }
  return s.done();
}
