// Checks of transport on a graph (transport.hpp) against closed forms and the scheme's invariants.
//
// - Content: on a random graph with random conductances and an advection that is not divergence-free,
//   Σ V c is kept to rounding over 2000 steps; with sources it rises by h Σ s a step.
// - A uniform concentration is steady under a circulation (divergence-free advection) and not under a
//   random advection: the condition is the one stated.
// - Diffusion on a long chain from a unit amount at one vertex: the semi-discrete solution is the
//   discrete heat kernel e^{−2gt} I_n(2gt) (I_n the modified Bessel function, summed from its series);
//   Heun's steps reach it at order 2 in h, and keep the variance Σ n² c = 2gt exactly (each step adds
//   2gh, as Σ n² A²c = 2g Σ Ac = 0).
// - Upwind advection of a pulse along a chain: its mean moves at u/V exactly (Σ n Ac = (u/V) Σ c and
//   Σ n A²c = 0); at the explicit bound no concentration turns negative, at 1.5 times the bound some do.
// - A circulation with diffusion at the bound keeps every concentration within its initial range.
// - The steady state with sources summing to zero solves the linear system to rounding, has the given
//   content, and is where a long integration ends; matrix() reproduces rate().

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <physis/continuum/transport.hpp>

using calx::verification::Suite;
namespace pc = physis::continuum;

namespace {

// A random connected graph: a spanning path through a random order, then extra edges.
pc::GraphTransport random_graph(calx::random::Stream& rng, std::size_t n, std::size_t extra) {
  pc::GraphTransport T;
  T.volume.resize(n);
  for (double& v : T.volume) v = rng.uniform(0.5, 2.0);
  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), 0);
  for (std::size_t i = n; i-- > 1;) std::swap(order[i], order[rng.below(i + 1)]);
  for (std::size_t i = 0; i + 1 < n; i++) T.edges.push_back({order[i], order[i + 1]});
  while (T.edges.size() < n - 1 + extra) {
    const std::size_t a = rng.below(n), b = rng.below(n);
    if (a != b) T.edges.push_back({a, b});
  }
  for (std::size_t e = 0; e < T.edges.size(); e++) T.conductance.push_back(rng.uniform(0.1, 2.0));
  return T;
}

// e^{−x} I_n(x) by its series Σ_k (x/2)^{2k+n} / (k! (k+n)!), the first term in logarithms.
double scaled_bessel_i(int n, double x) {
  double term = std::exp(n * std::log(x / 2) - std::lgamma(n + 1.0) - x), sum = 0;
  for (int k = 0; k < 1000; k++) {
    sum += term;
    term *= (x / 2) * (x / 2) / ((k + 1.0) * (k + 1.0 + n));
    if (term < 1e-18 * sum) break;
  }
  return sum;
}

std::vector<double> chain_heat(double g, double h, double T, std::size_t n, std::size_t centre) {
  pc::GraphTransport chain;
  chain.volume.assign(n, 1.0);
  for (std::size_t i = 0; i + 1 < n; i++) chain.edges.push_back({i, i + 1});
  chain.conductance.assign(n - 1, g);
  std::vector<double> c(n, 0.0);
  c[centre] = 1;
  const long steps = std::lround(T / h);
  for (long k = 0; k < steps; k++) c = chain.step(c, h);
  return c;
}

}  // namespace

int main() {
  Suite s;
  calx::random::Stream rng(20261030, 17, 0, 0);

  // ---- content ---------------------------------------------------------------------------------------------------------------
  {
    auto T = random_graph(rng, 30, 40);
    for (std::size_t e = 0; e < T.edges.size(); e++) T.advection.push_back(rng.uniform(-1, 1));
    std::vector<double> c(30), src(30);
    for (double& x : c) x = rng.uniform(0, 3);
    for (double& x : src) x = rng.uniform(-0.2, 0.5);
    const double h = 0.5 * T.stable_step(), c0 = T.content(c);
    auto d = c;
    for (int k = 0; k < 2000; k++) c = T.step(c, h), d = T.step(d, h, src);
    s.check("a random graph, advection not divergence-free: Σ V c kept over 2000 steps (relative)", std::fabs(T.content(c) - c0) / c0, 1e-13);
    const double added = 2000 * h * std::accumulate(src.begin(), src.end(), 0.0);
    s.check("with sources: Σ V c rises by h Σ s a step (relative)", std::fabs(T.content(d) - c0 - added) / (c0 + added), 1e-13);
  }

  // ---- uniform states and divergence ------------------------------------------------------------------------------------------
  {
    auto T = random_graph(rng, 12, 10);
    // A circulation round the spanning path and back along a closing edge: divergence-free.
    pc::GraphTransport C = T;
    C.advection.assign(C.edges.size(), 0.0);
    for (std::size_t e = 0; e + 1 < 12; e++) C.advection[e] = 0.7;
    C.edges.push_back({C.edges[10][1], C.edges[0][0]});
    C.conductance.push_back(0.3);
    C.advection.push_back(0.7);
    const std::vector<double> uniform(12, 2.5);
    double worst = 0;
    for (double r : C.rate(uniform)) worst = std::max(worst, std::fabs(r));
    s.check("a uniform concentration is steady under a circulation", worst, 1e-15);
    pc::GraphTransport R = T;
    for (std::size_t e = 0; e < R.edges.size(); e++) R.advection.push_back(rng.uniform(-1, 1));
    double moved = 0;
    for (double r : R.rate(uniform)) moved = std::max(moved, std::fabs(r));
    s.require("and not under an advection with divergence", moved > 1e-3, "max |dc/dt| " + std::to_string(moved));

    // The discrete maximum principle at the explicit bound.
    std::vector<double> c(12);
    for (double& x : c) x = rng.uniform(1, 4);
    const double lo = *std::min_element(c.begin(), c.end()), hi = *std::max_element(c.begin(), c.end());
    double out = 0;
    const double h = C.stable_step();
    for (int k = 0; k < 500; k++) {
      c = C.step(c, h);
      for (double x : c) out = std::max({out, lo - x, x - hi});
    }
    s.check("a circulation with diffusion at the explicit bound: concentrations stay within their initial range", out, 1e-14);
  }

  // ---- diffusion on a chain: the discrete heat kernel --------------------------------------------------------------------------
  {
    const double g = 0.8, Tend = 20;
    const std::size_t n = 401, centre = 200;
    double e1 = 0, e2 = 0, var = 0;
    const auto c1 = chain_heat(g, 0.1, Tend, n, centre), c2 = chain_heat(g, 0.05, Tend, n, centre);
    for (std::size_t i = 0; i < n; i++) {
      const int m = std::abs(static_cast<int>(i) - static_cast<int>(centre));
      const double exact = scaled_bessel_i(m, 2 * g * Tend);
      e1 = std::max(e1, std::fabs(c1[i] - exact));
      e2 = std::max(e2, std::fabs(c2[i] - exact));
      var += static_cast<double>(m) * m * c2[i];
    }
    s.check("diffusion on a chain: Heun's steps reach the heat kernel e^{−2gt} I_n(2gt) at h = 0.05 s", e2, 3e-7);
    s.check("and at order 2 in h", calx::verification::order(e1, e2) - 2, 0.05);
    s.check("the variance Σ n² c = 2gt exactly", (var - 2 * g * Tend) / (2 * g * Tend), 1e-13);
  }

  // ---- upwind advection of a pulse ------------------------------------------------------------------------------------------------
  {
    const std::size_t n = 401;
    const double u = 0.6;
    pc::GraphTransport chain;
    chain.volume.assign(n, 1.0);
    for (std::size_t i = 0; i + 1 < n; i++) chain.edges.push_back({i, i + 1});
    chain.advection.assign(n - 1, u);
    std::vector<double> c0(n);
    for (std::size_t i = 0; i < n; i++) c0[i] = std::exp(-0.5 * std::pow((static_cast<double>(i) - 100) / 4, 2));
    const auto mean = [&](const std::vector<double>& c) {
      double m = 0, w = 0;
      for (std::size_t i = 0; i < n; i++) m += i * c[i], w += c[i];
      return m / w;
    };
    const double h = chain.stable_step();
    auto c = c0;
    double negative = 0;
    for (int k = 0; k < 200; k++) {
      c = chain.step(c, h);
      for (double x : c) negative = std::min(negative, x);
    }
    s.check("upwind advection: the pulse's mean moves at u/V, 200 steps at the bound", mean(c) - mean(c0) - u * 200 * h, 1e-10, "cells");
    s.require("at the explicit bound h = V/u no concentration turns negative", negative >= 0, "least " + std::to_string(negative));
    auto b = c0;
    double beyond = 0;
    for (int k = 0; k < 40; k++) {
      b = chain.step(b, 1.5 * h);
      for (double x : b) beyond = std::min(beyond, x);
    }
    s.require("at 1.5 times the bound some do", beyond < -1e-6, "least " + std::to_string(beyond));
  }

  // ---- the steady state -----------------------------------------------------------------------------------------------------------
  {
    auto T = random_graph(rng, 12, 14);
    for (std::size_t e = 0; e < T.edges.size(); e++) T.advection.push_back(rng.uniform(-0.8, 0.8));
    std::vector<double> src(12);
    double total = 0;
    for (double& x : src) total += x = rng.uniform(-1, 1);
    for (double& x : src) x -= total / 12;
    const auto cs = T.steady(src, 5.0);
    double residual = 0;
    for (double r : T.rate(cs, src)) residual = std::max(residual, std::fabs(r));
    s.check("steady: dc/dt = 0 at the solution of the linear system", residual, 1e-13);
    s.check("steady: with the content asked for", T.content(cs) - 5.0, 1e-13);
    double volume = 0;
    for (double v : T.volume) volume += v;
    std::vector<double> c(12, 5.0 / volume);
    const double h = 0.9 * T.stable_step();
    for (int k = 0; k < 20000; k++) c = T.step(c, h, src);
    double gap = 0;
    for (std::size_t i = 0; i < 12; i++) gap = std::max(gap, std::fabs(c[i] - cs[i]));
    s.check("a long integration ends at the steady state", gap, 1e-9);

    std::vector<double> x(12);
    for (double& v : x) v = rng.uniform(0, 2);
    const auto A = T.matrix();
    const auto r = T.rate(x, src);
    double diff = 0;
    for (std::size_t i = 0; i < 12; i++) {
      double Ax = src[i] / T.volume[i];
      for (std::size_t j = 0; j < 12; j++) Ax += A(i, j) * x[j];
      diff = std::max(diff, std::fabs(Ax - r[i]));
    }
    s.check("matrix(): A c + V⁻¹ s = rate(c, s)", diff, 1e-14);
  }

  return s.done();
}
