// Checks of the probability of a record of events (event_chains.hpp).
//
// - On twenty random records of up to five cells of two or three states, with laws, transitions, observations and
//   pieces joining two cells, Z by the contraction equals Z by enumerating every first and last state of every cell.
// - A cell cut by an identity piece, its propagator split as P(τ₁) P(τ₂) = P(τ₁ + τ₂), gives the same Z: slicing time
//   is a refinement (the design's Theorem 4.5).
// - With stochastic kernels and propagators and no observation, Z = 1 (the design's Proposition 4.6).
// - Independent chains stay in factors of one cell: the largest factor has as many entries as one cell's states.
// - Faulty records refused: a cell ended before it begins, a kernel of the wrong size.

#include <cmath>
#include <functional>
#include <string>
#include <vector>

#include <calx/linalg/dense.hpp>
#include <calx/linalg/expm.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <physis/statistical_mechanics/event_chains.hpp>

using calx::verification::Suite;
namespace sm = physis::statistical_mechanics;

namespace {

calx::random::Stream rng(7031, 20, 0, 1);

calx::linalg::Matrix stochastic(std::size_t n) {
  calx::linalg::Matrix P(n);
  for (std::size_t i = 0; i < n; i++) {
    double s = 0;
    for (std::size_t j = 0; j < n; j++) s += (P(i, j) = rng.uniform(0.05, 1.0));
    for (std::size_t j = 0; j < n; j++) P(i, j) /= s;
  }
  return P;
}

// A generator's propagator at τ, by calx's exponential.
calx::linalg::Matrix propagator(const calx::linalg::Matrix& L, double tau) {
  auto A = L;
  for (double& v : A.a) v *= tau;
  return calx::linalg::expm(A);
}

calx::linalg::Matrix generator(std::size_t n) {
  calx::linalg::Matrix L(n);
  for (std::size_t i = 0; i < n; i++) {
    double s = 0;
    for (std::size_t j = 0; j < n; j++)
      if (j != i) s += (L(i, j) = rng.uniform(0.1, 1.5));
    L(i, i) = -s;
  }
  return L;
}

// Z by enumerating every first and last state of every cell.
double brute(const sm::EventChains& M) {
  const std::size_t C = M.cells.size();
  std::vector<std::size_t> first(C, 0), last(C, 0);
  double Z = 0;
  std::function<void(std::size_t)> go = [&](std::size_t k) {
    if (k == 2 * C) {
      double w = 1;
      for (std::size_t c = 0; c < C; c++) w *= M.cells[c].propagator(first[c], last[c]);
      for (const auto& q : M.pieces) {
        std::size_t a = 0;
        for (auto c : q.in) a = a * M.cells[c].states + last[c];
        for (auto c : q.out) a = a * M.cells[c].states + first[c];
        w *= q.kernel[a];
      }
      Z += w;
      return;
    }
    auto& v = k < C ? first[k] : last[k - C];
    const std::size_t n = M.cells[k < C ? k : k - C].states;
    for (v = 0; v < n; v++) go(k + 1);
  };
  go(0);
  return Z;
}

// A random record: each cell begins at a law or at a piece joining it to a cell that ends there, is observed at
// times, and ends at a sum.
sm::EventChains random_record(bool observe) {
  sm::EventChains M;
  const std::size_t cells = 2 + rng.next64() % 4;
  std::vector<std::size_t> live;
  auto add = [&]() {
    const std::size_t n = 2 + rng.next64() % 2;
    M.cells.push_back({n, propagator(generator(n), rng.uniform(0.2, 2.0))});
    return M.cells.size() - 1;
  };
  for (std::size_t k = 0; k < cells; k++) {
    if (!live.empty() && rng.uniform() < 0.4) {
      // A transition from a live cell to a new one (a rewrite carrying the state), or a piece joining two live cells
      // into two new ones (a rule over both).
      if (live.size() >= 2 && rng.uniform() < 0.5) {
        const std::size_t a = live[live.size() - 2], b = live.back();
        live.resize(live.size() - 2);
        const std::size_t c = add(), d = add();
        sm::ChainPiece q{{a, b}, {c, d}, {}};
        const std::size_t size = M.cells[a].states * M.cells[b].states * M.cells[c].states * M.cells[d].states;
        for (std::size_t i = 0; i < size; i++) q.kernel.push_back(rng.uniform(0.0, 1.0));
        M.pieces.push_back(q);
        live.push_back(c), live.push_back(d);
      } else {
        const std::size_t a = live.back();
        live.pop_back();
        const std::size_t c = add();
        sm::ChainPiece q{{a}, {c}, {}};
        for (std::size_t i = 0; i < M.cells[a].states; i++)
          for (std::size_t j = 0; j < M.cells[c].states; j++) q.kernel.push_back(rng.uniform(0.0, 1.0));
        M.pieces.push_back(q);
        live.push_back(c);
      }
    } else {
      const std::size_t c = add();
      sm::ChainPiece q{{}, {c}, {}};
      double s = 0;
      for (std::size_t i = 0; i < M.cells[c].states; i++) s += q.kernel.emplace_back(rng.uniform(0.1, 1.0));
      for (double& v : q.kernel) v /= s;
      M.pieces.push_back(q);
      live.push_back(c);
    }
    if (observe && !live.empty()) {
      // An observation of a live cell: it ends and a new one of its states begins, by diag ℓ.
      const std::size_t a = live.back();
      live.pop_back();
      const std::size_t n = M.cells[a].states;
      M.cells.push_back({n, propagator(generator(n), rng.uniform(0.2, 1.0))});
      const std::size_t c = M.cells.size() - 1;
      sm::ChainPiece q{{a}, {c}, std::vector<double>(n * n, 0.0)};
      for (std::size_t i = 0; i < n; i++) q.kernel[i * n + i] = rng.uniform(0.05, 1.0);
      M.pieces.push_back(q);
      live.push_back(c);
    }
  }
  for (auto c : live) M.pieces.push_back({{c}, {}, std::vector<double>(M.cells[c].states, 1.0)});
  return M;
}

bool refused(const sm::EventChains& M) {
  try {
    (void)sm::contract(M);
  } catch (const std::invalid_argument&) {
    return true;
  }
  return false;
}

}  // namespace

int main() {
  Suite s;
  // Against enumeration.
  {
    double worst = 0;
    for (int trial = 0; trial < 20; trial++) {
      const auto M = random_record(true);
      const double z = sm::contract(M).log_z, b = std::log(brute(M));
      worst = std::max(worst, std::fabs(z - b));
    }
    s.require("twenty random records: log Z by the contraction equals log Z by enumeration", worst < 1e-10, std::to_string(worst));
  }
  // A cell cut by an identity piece.
  {
    const auto L = generator(3);
    sm::EventChains A, B;
    A.cells = {{3, propagator(L, 2.0)}};
    A.pieces = {{{}, {0}, {0.2, 0.5, 0.3}}, {{0}, {}, {0.9, 0.1, 0.4}}};
    B.cells = {{3, propagator(L, 0.7)}, {3, propagator(L, 1.3)}};
    std::vector<double> I(9, 0.0);
    I[0] = I[4] = I[8] = 1;
    B.pieces = {{{}, {0}, {0.2, 0.5, 0.3}}, {{0}, {1}, I}, {{1}, {}, {0.9, 0.1, 0.4}}};
    s.require("a cell cut by an identity piece (P(0.7) P(1.3) = P(2)) gives the same Z",
              std::fabs(sm::contract(A).log_z - sm::contract(B).log_z) < 1e-13);
  }
  // No observation, stochastic kernels: Z = 1.
  {
    double worst = 0;
    for (int trial = 0; trial < 10; trial++) {
      sm::EventChains M;
      std::vector<std::size_t> live;
      for (int k = 0; k < 4; k++) {
        const std::size_t n = 2 + rng.next64() % 2;
        M.cells.push_back({n, propagator(generator(n), 1.0)});
        const std::size_t c = M.cells.size() - 1;
        if (!live.empty() && k % 2) {
          const std::size_t a = live.back();
          live.pop_back();
          sm::ChainPiece q{{a}, {c}, {}};
          for (std::size_t i = 0; i < M.cells[a].states; i++) {
            double t = 0;
            std::vector<double> row(n);
            for (auto& v : row) t += (v = rng.uniform(0.1, 1.0));
            for (auto v : row) q.kernel.push_back(v / t);
          }
          M.pieces.push_back(q);
        } else {
          std::vector<double> law(n, 1.0 / static_cast<double>(n));
          M.pieces.push_back({{}, {c}, law});
        }
        live.push_back(c);
      }
      for (auto c : live) M.pieces.push_back({{c}, {}, std::vector<double>(M.cells[c].states, 1.0)});
      worst = std::max(worst, std::fabs(sm::contract(M).log_z));
    }
    s.require("stochastic kernels and propagators, no observation: Z = 1", worst < 1e-13, std::to_string(worst));
  }
  // Independent chains stay apart.
  {
    sm::EventChains M;
    for (int k = 0; k < 40; k++) {
      M.cells.push_back({3, propagator(generator(3), 1.0)});
      M.pieces.push_back({{}, {M.cells.size() - 1}, {0.3, 0.3, 0.4}});
    }
    for (std::size_t c = 0; c < 40; c++) M.pieces.push_back({{c}, {}, {1.0, 0.5, 0.2}});
    const auto R = sm::contract(M);
    s.require("forty independent chains: the largest factor one cell's three states", R.largest == 3 && std::isfinite(R.log_z));
  }
  // Faulty records.
  {
    sm::EventChains early;
    early.cells = {{2, calx::linalg::Matrix::identity(2)}};
    early.pieces = {{{0}, {}, {1, 1}}, {{}, {0}, {0.5, 0.5}}};
    sm::EventChains size;
    size.cells = {{2, calx::linalg::Matrix::identity(2)}};
    size.pieces = {{{}, {0}, {0.5, 0.25, 0.25}}, {{0}, {}, {1, 1}}};
    s.require("a cell ended before it begins, and a kernel of the wrong size, are refused", refused(early) && refused(size));
  }
  return s.done();
}
