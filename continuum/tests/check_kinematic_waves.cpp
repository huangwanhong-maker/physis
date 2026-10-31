// Checks of kinematic waves (kinematic_waves.hpp) against closed forms, the integral form of the law,
// and the invariants of the scheme and its nodes.
//
// - The fluxes: ρ_c and q_max of Greenshields's and the triangular law in closed form; demand rising
//   and supply falling to meet at q_max.
// - Godunov's flux min(D(ρ_l), S(ρ_r)) equal to q of the Riemann solution at ξ = 0, for random states
//   of both laws (Lebacque's identity: two derivations of one number).
// - The Riemann solution: Greenshields's shock at v(1 − (ρ_l + ρ_r)/ρ_j); inside a fan, the density whose
//   wave speed is ξ; Lax's entropy condition at every shock; the triangular fan at ρ_c between −w and v;
//   and the integral form ∫_{−Ξ}^{Ξ} ρ(ξ) dξ = Ξ(ρ_l + ρ_r) + q(ρ_l) − q(ρ_r) for random states of both laws.
// - Convergence of Godunov's scheme on one link, its ends the Riemann problem's states (an entry with
//   demand D(ρ_l), an exit with supply S(ρ_r)): at a shock the L1 error, averaged over the second half
//   of the run (a discrete shock's profile cycles as it crosses the cells), falls at order 1; at a
//   centred rarefaction it falls as Δx log(ℓ/Δx), ℓ the fan's width (Harabetian 1988), so that the
//   order of E/log(ℓ/Δx) is 1 (calx's order).
// - A network with a merge and a diverge, time-varying demands and a choked exit: the budget closes to
//   rounding (10⁻¹³ of what passed), every density stays in [0, ρ_j], demand is refused at an entry once
//   the queue reaches it.
// - Junctions: under FIFO a diverge keeps its turning fractions exactly and meets the binding supply;
//   non-FIFO lets the free movement pass; a merge shares a short supply in proportion to demand; on
//   random 2 × 3 junctions under both rules no supply is exceeded, no demand exceeded, flow conserved.
// - A bottleneck (a link of lower jam density downstream): the queue's tail moves upstream at the
//   Rankine–Hugoniot speed (q_in − q_B)/(ρ_free − ρ_queue) of the closed form.
// - A scheduled change of a link's law: the queue discharges at the new capacity from the change on; a
//   step allowed by the links' own CFL bound is refused when a scheduled law is faster.

#include <algorithm>
#include <cmath>
#include <functional>
#include <numbers>
#include <string>
#include <vector>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <physis/continuum/kinematic_waves.hpp>

using calx::verification::Suite;
namespace pc = physis::continuum;
using pc::Flux;
using pc::KinematicNetwork;

namespace {

// Godunov's scheme on one link of length L centred on x = 0, from Riemann data, to time T; the L1
// error against the exact cell averages, at T or averaged over `samples` instants of the second half.
double riemann_error(const Flux& f, double rl, double rr, std::size_t N, double L, double T, int samples) {
  KinematicNetwork net(2, {pc::Link{0, 1, L, N, f}});
  net.entry(0, [&](double) { return f.demand(rl); });
  net.exit(1, [&](double) { return f.supply(rr); });
  const double dx = L / static_cast<double>(N);
  auto s = net.start();
  for (std::size_t k = 0; k < N; k++) s.density[0][k] = (-L / 2 + (k + 0.5) * dx) < 0 ? rl : rr;
  const auto steps = static_cast<std::size_t>(std::ceil(T / (0.9 * net.cfl())));
  const double h = T / static_cast<double>(steps);
  const auto error = [&](double t) {
    const int sub = 64;
    double e = 0;
    for (std::size_t k = 0; k < N; k++) {
      double avg = 0;
      for (int j = 0; j < sub; j++) avg += pc::riemann(f, rl, rr, (-L / 2 + (k + (j + 0.5) / sub) * dx) / t);
      e += std::fabs(s.density[0][k] - avg / sub) * dx;
    }
    return e;
  };
  double sum = 0;
  int count = 0;
  const std::size_t half = steps / 2, every = std::max<std::size_t>(1, (steps - half) / static_cast<std::size_t>(std::max(samples, 1)));
  for (std::size_t k = 1; k <= steps; k++) {
    s = net.step(s, h);
    if (samples > 1 && k >= half && (k - half) % every == 0) sum += error(s.t), count++;
  }
  return samples > 1 ? sum / count : error(s.t);
}

// The flows at a junction from given densities: incoming links' last cells and outgoing links' first cells.
struct Junction {
  std::vector<double> outflow, inflow, demand, supply;
};
Junction junction_flows(const std::vector<Flux>& in, const std::vector<Flux>& out, const std::vector<double>& rho_in, const std::vector<double>& rho_out,
                        const std::vector<std::vector<double>>& beta, KinematicNetwork::Rule rule) {
  // node 0 the junction; incoming link i from node 1 + i; outgoing link j to node 1 + in.size() + j.
  std::vector<pc::Link> links;
  for (std::size_t i = 0; i < in.size(); i++) links.push_back({1 + i, 0, 10, 1, in[i]});
  for (std::size_t j = 0; j < out.size(); j++) links.push_back({0, 1 + in.size() + j, 10, 1, out[j]});
  KinematicNetwork net(1 + in.size() + out.size(), links);
  net.junction(0, beta, rule);
  auto s = net.start();
  for (std::size_t i = 0; i < in.size(); i++) s.density[i][0] = rho_in[i];
  for (std::size_t j = 0; j < out.size(); j++) s.density[in.size() + j][0] = rho_out[j];
  const auto F = net.link_flows(s);
  Junction J;
  for (std::size_t i = 0; i < in.size(); i++) J.outflow.push_back(F[i][1]), J.demand.push_back(in[i].demand(rho_in[i]));
  for (std::size_t j = 0; j < out.size(); j++) J.inflow.push_back(F[in.size() + j][0]), J.supply.push_back(out[j].supply(rho_out[j]));
  return J;
}

template <class F>
bool throws(F&& f) {
  try {
    f();
  } catch (const std::invalid_argument&) {
    return true;
  }
  return false;
}

}  // namespace

int main() {
  Suite s;
  calx::random::Stream rng(20261030, 16, 0, 0);
  const Flux G = Flux::greenshields(1.4, 5.0);
  const Flux T = Flux::triangular(1.4, 0.7, 5.0);

  // ---- the fluxes ------------------------------------------------------------------------------------------------------------
  {
    s.check("Greenshields: ρ_c = ρ_j/2 and q_max = v ρ_j/4", std::max(std::fabs(G.critical() - 2.5), std::fabs(G.capacity() - 1.75)), 1e-15);
    s.check("triangular: ρ_c = w ρ_j/(v + w) and q_max = v w ρ_j/(v + w)", std::max(std::fabs(T.critical() - 0.7 * 5 / 2.1), std::fabs(T.capacity() - 1.4 * 0.7 * 5 / 2.1)), 1e-15);
    bool monotone = true;
    for (const Flux& f : {G, T}) {
      for (int k = 0; k < 500; k++) {
        const double a = f.jam_density * k / 500, b = f.jam_density * (k + 1) / 500;
        monotone = monotone && f.demand(b) >= f.demand(a) && f.supply(b) <= f.supply(a) && f.demand(a) <= f.capacity() && f.supply(a) <= f.capacity();
      }
      monotone = monotone && f.demand(f.critical()) == f.capacity() && f.supply(f.critical()) == f.capacity() && f.q(0) == 0 && std::fabs(f.q(f.jam_density)) < 1e-15;
    }
    s.require("demand rises and supply falls to meet at q_max; q(0) = q(ρ_j) = 0", monotone);
  }

  // ---- Godunov's flux is q of the Riemann solution at ξ = 0 ----------------------------------------------------------------
  {
    double worst = 0;
    for (const Flux& f : {G, T})
      for (int k = 0; k < 20000; k++) {
        const double rl = rng.uniform(0, f.jam_density), rr = rng.uniform(0, f.jam_density);
        worst = std::max(worst, std::fabs(pc::godunov_flux(f, rl, rr) - f.q(pc::riemann(f, rl, rr, 0.0))));
      }
    s.check("Godunov's flux min(D(ρ_l), S(ρ_r)) = q(ρ(ξ = 0)) of the Riemann solution, 40000 random pairs of both laws", worst, 1e-15, "/s");
  }

  // ---- the Riemann problem --------------------------------------------------------------------------------------------------
  {
    bool shocks = true, entropy = true;
    double fan = 0;
    for (int k = 0; k < 5000; k++) {
      const double a = rng.uniform(0, 5), b = rng.uniform(0, 5), rl = std::min(a, b), rr = std::max(a, b);
      if (rr - rl < 1e-6) continue;
      const double speed = 1.4 * (1 - (rl + rr) / 5);
      shocks = shocks && pc::riemann(G, rl, rr, speed - 1e-9) == rl && pc::riemann(G, rl, rr, speed + 1e-9) == rr;
      for (const Flux& f : {G, T}) {
        const double sigma = (f.q(rr) - f.q(rl)) / (rr - rl);
        entropy = entropy && f.dq(rl) >= sigma - 1e-12 && sigma >= f.dq(rr) - 1e-12;
      }
      // the fan of (rr, rl): inside it, the density whose wave speed is ξ
      const double lo = G.dq(rr), hi = G.dq(rl);
      for (int j = 1; j < 10; j++) {
        const double xi = lo + (hi - lo) * j / 10;
        fan = std::max(fan, std::fabs(G.dq(pc::riemann(G, rr, rl, xi)) - xi));
      }
    }
    s.require("Greenshields: a shock (ρ_l < ρ_r) at v(1 − (ρ_l + ρ_r)/ρ_j), 5000 random pairs", shocks);
    s.require("Lax's entropy condition q′(ρ_l) ≥ s ≥ q′(ρ_r) at every shock, both laws", entropy);
    s.check("Greenshields: inside a fan, q′(ρ(ξ)) = ξ", fan, 1e-14, "m/s");
    const double rc = T.critical();
    const bool tri = pc::riemann(T, 4.0, 0.5, -0.8) == 4.0 && pc::riemann(T, 4.0, 0.5, -0.6) == rc && pc::riemann(T, 4.0, 0.5, 1.3) == rc &&
                     pc::riemann(T, 4.0, 0.5, 1.5) == 0.5 && pc::riemann(T, 1.0, 0.5, 1.39) == 1.0 && pc::riemann(T, 1.0, 0.5, 1.41) == 0.5;
    s.require("triangular: a congested-to-free fan holds ρ_c between −w and v; on one branch, a jump at v", tri);

    // The integral form over [−Ξ, Ξ] beyond every wave: ∫ρ dξ = Ξ(ρ_l + ρ_r) + q(ρ_l) − q(ρ_r).
    double worst = 0;
    for (const Flux& f : {G, T})
      for (int k = 0; k < 40; k++) {
        const double rl = rng.uniform(0, f.jam_density), rr = rng.uniform(0, f.jam_density), X = 2 * f.max_speed();
        const int M = 400000;
        double integral = 0;
        for (int j = 0; j < M; j++) integral += pc::riemann(f, rl, rr, -X + (j + 0.5) * 2 * X / M) * 2 * X / M;
        const double expected = X * (rl + rr) + f.q(rl) - f.q(rr);
        worst = std::max(worst, std::fabs(integral - expected) / (X * f.jam_density));
      }
    s.check("the integral form ∫ρ(ξ)dξ = Ξ(ρ_l + ρ_r) + q(ρ_l) − q(ρ_r), 80 random problems of both laws (relative to Ξρ_j)", worst, 2e-5);
  }

  // ---- convergence of Godunov's scheme --------------------------------------------------------------------------------------
  {
    const double L = 200, Tend = 40;
    std::vector<double> shock, fan;
    for (std::size_t N : {100, 200, 400, 800}) {
      shock.push_back(riemann_error(G, 0.5, 3.0, N, L, Tend, 40));
      fan.push_back(riemann_error(G, 4.0, 0.5, N, L, Tend, 1));
    }
    double dev = 0;
    for (int i = 0; i < 3; i++) dev = std::max(dev, std::fabs(calx::verification::order(shock[i], shock[i + 1]) - 1));
    s.check("a shock (ρ 0.5 → 3 per m): the time-averaged L1 error falls at order 1 (N = 100 … 800)", dev, 0.1);
    const double ell = (G.dq(0.5) - G.dq(4.0)) * Tend;
    double dfan = 0, raw = 1;
    for (int i = 0; i < 3; i++) {
      const double dx1 = L / (100 << i), dx2 = dx1 / 2;
      dfan = std::max(dfan, std::fabs(calx::verification::order(fan[i] / std::log(ell / dx1), fan[i + 1] / std::log(ell / dx2)) - 1));
      raw = std::min(raw, calx::verification::order(fan[i], fan[i + 1]));
    }
    s.check("a centred rarefaction (ρ 4 → 0.5 per m): the L1 error falls as Δx log(ℓ/Δx), order 1 after the logarithm", dfan, 0.1);
    s.require("a centred rarefaction: its raw order lies between ½ and 1", raw > 0.5 && raw < 1, "least " + std::to_string(raw));
  }

  // ---- a network with a merge and a diverge: the budget ------------------------------------------------------------------------
  {
    // Nodes: entries 0, 1; merge 2; diverge 3; exits 4 (choked for a while), 5.
    const std::vector<pc::Link> links{{0, 2, 120, 24, G}, {1, 2, 80, 16, T}, {2, 3, 150, 30, G}, {3, 4, 60, 12, T}, {3, 5, 90, 18, G}};
    KinematicNetwork net(6, links);
    net.entry(0, [](double t) { return 1.2 + 0.8 * std::sin(t / 30); });
    net.entry(1, [](double t) { return t < 400 ? 1.4 : 0.3; });
    net.junction(2, {{1.0}, {1.0}});
    net.junction(3, {{0.35, 0.65}}, KinematicNetwork::Rule::fifo);
    net.exit(4, [](double t) { return t > 100 && t < 500 ? 0.05 : 10.0; });
    net.exit(5);
    KinematicNetwork::Budget budget;
    auto st = net.start();
    const double h = 0.9 * net.cfl();
    double low = 0, high = 0, residual = 0;
    for (int k = 0; k < 4000; k++) {
      st = net.step(st, h, &budget);
      residual = std::max(residual, std::fabs(budget.residual()) / (budget.initial + budget.entered));
      for (std::size_t l = 0; l < links.size(); l++)
        for (double r : st.density[l]) low = std::min(low, r), high = std::max(high, r - links[l].flux.jam_density);
    }
    s.check("a merge and a FIFO diverge, a choked exit: initial + entered − left − stored, relative to what passed", residual, 1e-13);
    s.check("every density stays in [0, ρ_j]", std::max(-low, high), 1e-12, "per m");
    s.require("the queue reaches an entry, whose demand is then refused", budget.refused > 1, "refused " + std::to_string(budget.refused));
  }

  // ---- junctions ------------------------------------------------------------------------------------------------------------
  {
    const double rc = G.critical();
    // A diverge, β = (0.3, 0.7), the first outgoing link nearly jammed.
    const auto fifo = junction_flows({G}, {G, G}, {rc}, {4.6, 0.0}, {{0.3, 0.7}}, KinematicNetwork::Rule::fifo);
    const double S1 = fifo.supply[0];
    s.check("FIFO diverge: the flows split exactly 0.3 : 0.7, and the blocked branch takes its supply",
            std::max(std::fabs(fifo.inflow[0] / fifo.outflow[0] - 0.3), std::fabs(fifo.inflow[0] - S1) / S1), 1e-15);
    const auto free = junction_flows({G}, {G, G}, {rc}, {4.6, 0.0}, {{0.3, 0.7}}, KinematicNetwork::Rule::non_fifo);
    s.check("non-FIFO diverge: the free movement passes (0.7 D), the blocked one takes its supply",
            std::max(std::fabs(free.inflow[1] - 0.7 * free.demand[0]), std::fabs(free.inflow[0] - S1)), 1e-15, "/s");
    // A merge into a nearly jammed link.
    const auto merge = junction_flows({G, T}, {G}, {2.0, 1.5}, {4.0}, {{1.0}, {1.0}}, KinematicNetwork::Rule::non_fifo);
    const double D1 = merge.demand[0], D2 = merge.demand[1], S = merge.supply[0];
    s.check("a merge shares a short supply in proportion to demand: D_i S/(D₁ + D₂)",
            std::max(std::fabs(merge.outflow[0] - D1 * S / (D1 + D2)), std::fabs(merge.outflow[1] - D2 * S / (D1 + D2))), 1e-15, "/s");

    bool honoured = true;
    double conserved = 0;
    for (int trial = 0; trial < 2000; trial++) {
      std::vector<std::vector<double>> beta(2, std::vector<double>(3));
      for (auto& row : beta) {
        double sum = 0;
        for (double& b : row) sum += b = rng.uniform() < 0.2 ? 0.0 : rng.uniform();
        if (sum == 0) row[0] = sum = 1;
        for (double& b : row) b /= sum;
      }
      const std::vector<double> rin{rng.uniform(0, 5), rng.uniform(0, 5)}, rout{rng.uniform(0, 5), rng.uniform(0, 5), rng.uniform(0, 5)};
      for (auto rule : {KinematicNetwork::Rule::non_fifo, KinematicNetwork::Rule::fifo}) {
        const auto J = junction_flows({G, T}, {T, G, G}, rin, rout, beta, rule);
        double in = 0, out = 0;
        for (std::size_t i = 0; i < 2; i++) honoured = honoured && J.outflow[i] <= J.demand[i] * (1 + 1e-14) && J.outflow[i] >= 0, out += J.outflow[i];
        for (std::size_t j = 0; j < 3; j++) honoured = honoured && J.inflow[j] <= J.supply[j] * (1 + 1e-14) + 1e-300 && J.inflow[j] >= 0, in += J.inflow[j];
        conserved = std::max(conserved, std::fabs(in - out));
      }
    }
    s.require("random 2 × 3 junctions under both rules: no supply exceeded, no demand exceeded", honoured);
    s.check("random 2 × 3 junctions: what leaves the incoming links enters the outgoing ones", conserved, 1e-14, "/s");

    KinematicNetwork diverge(3, {{0, 1, 10, 2, G}, {1, 2, 10, 2, G}, {1, 0, 10, 2, G}});
    s.require("an undeclared diverge is refused; so is a step beyond the CFL bound",
              throws([&] { diverge.step(diverge.start(), 0.1); }) && throws([&] {
                KinematicNetwork one(2, {{0, 1, 10, 2, G}});
                one.step(one.start(), 1.01 * one.cfl());
              }));
  }

  // ---- a bottleneck: spillback at the Rankine–Hugoniot speed ------------------------------------------------------------------
  {
    const Flux A = Flux::triangular(1.4, 0.7, 5.0), B = Flux::triangular(1.4, 0.7, 2.0);
    const double qin = 1.5, capB = B.capacity(), rho_free = qin / 1.4, rho_queue = 5.0 - capB / 0.7;
    const double expected = (qin - capB) / (rho_free - rho_queue);
    KinematicNetwork net(3, {{0, 1, 300, 300, A}, {1, 2, 100, 100, B}});
    net.entry(0, [&](double) { return qin; });
    net.exit(2);
    auto st = net.start();
    const double h = 0.9 * net.cfl();
    const double mid = 0.5 * (rho_free + rho_queue);
    const auto tail = [&] {  // where the density first crosses the midpoint, from upstream
      const auto& r = st.density[0];
      std::size_t i = 1;
      while (i < r.size() && r[i] <= mid) i++;
      return (i - 0.5) + (mid - r[i - 1]) / (r[i] - r[i - 1]);
    };
    double x1 = 0, x2 = 0;
    for (double target : {500.0, 1100.0}) {
      while (st.t < target - 1e-9) st = net.step(st, std::min(h, target - st.t));
      (target < 600 ? x1 : x2) = tail();
    }
    const double observed = (x2 - x1) / 600;
    s.check("a bottleneck: the queue's tail moves at (q_in − q_B)/(ρ_free − ρ_queue) = −0.2183 m/s, relative error", std::fabs(observed - expected) / std::fabs(expected), 1e-3);

    // The bottleneck opens at t = 1200: the queue discharges at the new capacity at once.
    KinematicNetwork open(3, {{0, 1, 300, 300, A}, {1, 2, 100, 100, B}});
    open.entry(0, [&](double) { return qin; });
    open.exit(2);
    open.schedule(1, [&](double t) { return t < 1200 ? B : A; });
    auto o = open.start();
    while (o.t < 1200 - 1e-9) o = open.step(o, std::min(h, 1200 - o.t));
    const double before = open.link_flows(KinematicNetwork::State{1199.0, o.density})[0][300];
    const double after = open.link_flows(KinematicNetwork::State{1200.0, o.density})[0][300];
    s.check("a scheduled opening: the junction passes q_B before it and the upstream link's capacity after",
            std::max(std::fabs(before - capB) / capB, std::fabs(after - A.capacity()) / A.capacity()), 1e-12);

    KinematicNetwork faster(2, {{0, 1, 100, 50, G}});
    faster.schedule(0, [&](double t) { return t < 10 ? G : Flux::greenshields(2.8, 5.0); });
    const double hb = faster.cfl();
    s.require("a scheduled faster law halves the CFL bound from its time on, and step() refuses the old step",
              std::fabs(faster.cfl(20) - hb / 2) < 1e-15 && !throws([&] { faster.step(faster.start(5), 0.9 * hb); }) && throws([&] { faster.step(faster.start(20), 0.9 * hb); }));
  }

  return s.done();
}
