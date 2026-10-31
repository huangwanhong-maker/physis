// Kinematic waves: a scalar conservation law with a concave flux, on one link and on networks of
// links joined at nodes (decision 0016).
//
// The law (Lighthill and Whitham 1955, Proc. R. Soc. A 229, 317; Richards 1956, Oper. Res. 4, 42). A
// density ρ(x, t) of walkers or vehicles per metre of a link, carried at the flow q(ρ) per second,
// obeys ∂ρ/∂t + ∂q(ρ)/∂x = 0. With q concave on [0, ρ_j] and q(0) = q(ρ_j) = 0, information travels
// at the kinematic wave speed q′(ρ), discontinuities at the Rankine–Hugoniot speed s = [q]/[ρ], and the
// entropy solution of a Riemann problem (ρ_l for x < 0, ρ_r for x > 0) is a shock when ρ_l < ρ_r and a
// rarefaction fan when ρ_l > ρ_r (LeVeque 2002, Finite Volume Methods for Hyperbolic Problems, ch. 11).
// Two fluxes are written:
//
//   Greenshields (1935, Proc. Highway Res. Board 14, 448):  q = v ρ (1 − ρ/ρ_j),  ρ_c = ρ_j/2,  q_max = v ρ_j/4;
//   triangular (Newell 1993, Transp. Res. B 27, 281):       q = min(v ρ, w (ρ_j − ρ)),  ρ_c = w ρ_j/(v + w),  q_max = v ρ_c.
//
// For a crowd, ρ is per square metre times the walkway's width, so that a link's flux carries its width.
//
// Godunov's scheme. A link is cut into equal cells; the flow through the interface between cells of
// densities ρ_l and ρ_r is Godunov's, the flux of the Riemann solution at x/t = 0, which for a concave
// q is min(D(ρ_l), S(ρ_r)) with the demand D(ρ) = q(min(ρ, ρ_c)) and the supply S(ρ) = q(max(ρ, ρ_c))
// (Daganzo 1994, Transp. Res. B 28, 269, the cell transmission model; Lebacque 1996, Proc. 13th ISTTT,
// 647). The cells' densities advance by ρ_k′ = ρ_k + (h/Δx)(F_k − F_{k+1}): every interface flow leaves
// one cell and enters the next, so the content Σ ρ Δx changes by what the nodes admit and release
// alone, to rounding. The scheme is monotone, and keeps every density in [0, ρ_j], under the
// Courant–Friedrichs–Lewy condition h max|q′| ≤ Δx (LeVeque 2002, ch. 12).
//
// Nodes. An entry admits min(its demand at t, the supply of its link's first cell); what it cannot
// admit is refused (counted, not queued). An exit takes from each of its incoming links D_i · min(1,
// S/Σ_k D_k), S its supply at t (unlimited by default). A junction with turning fractions β_ij (row i
// an incoming link, column j an outgoing one, Σ_j β_ij = 1) sends, with θ_j = min(1, S_j / Σ_k D_k β_kj):
//
//   non-FIFO (crowds, which overtake):        movement i → j carries D_i β_ij θ_j;
//   FIFO (Daganzo 1995, Transp. Res. B 29, 79): incoming i flows D_i · min_{j: β_ij > 0} θ_j, split by β_ij.
//
// Both honour every supply, conserve the flow through the node exactly, and never send more than a
// demand (Tampère et al. 2011, Transp. Res. B 45, 289, the requirements of node models); FIFO keeps
// the turning fractions of what moves, non-FIFO lets an unblocked movement pass a blocked one. A merge
// under either rule shares a short supply in proportion to the demands (Jin and Zhang 2003, Transp.
// Res. B 37, 521). A node left undeclared passes one incoming link into one outgoing link, merges
// several into one, and is closed (nothing passes) where it has no incoming or no outgoing link; a
// diverge must be declared with its turning fractions.
//
// Conventions. Densities are per metre (ρ_j per metre of link), flows per second, lengths in metres,
// speeds in metres per second. A state holds, for each link, its cells' densities from tail to head.
// A link's interfaces are numbered 0 (at its tail) to C (at its head); `link_flows` gives the flow
// through each. A node's incoming and outgoing links, and so the rows and columns of its β, are in the
// order of the links' indices. The boundary data (demands, supplies, scheduled fluxes) are read at the
// start of each step, held over it.
//
// Accuracy. First order: at a shock the L1 error is O(Δx); at a centred rarefaction O(Δx log(1/Δx))
// (Harabetian 1988, Commun. Math. Phys. 114, 527). Checks: tests/check_kinematic_waves.cpp.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace physis::continuum {

// ---- the flux law -------------------------------------------------------------------------------------------------------------

struct Flux {
  enum class Kind : std::uint8_t { greenshields, triangular };
  Kind kind = Kind::greenshields;
  double free_speed = 0;   // v_f, m/s
  double jam_density = 0;  // ρ_j, per metre of link (a crowd: per m² × width)
  double wave_speed = 0;   // w > 0, m/s, triangular only: q = min(v ρ, w (ρ_j − ρ))

  // q = v ρ (1 − ρ/ρ_j).
  static Flux greenshields(double v, double rho_j) {
    if (!(v > 0 && rho_j > 0 && std::isfinite(v) && std::isfinite(rho_j))) throw std::invalid_argument("physis::continuum::Flux::greenshields: v > 0 and ρ_j > 0");
    return {Kind::greenshields, v, rho_j, 0};
  }
  // q = min(v ρ, w (ρ_j − ρ)).
  static Flux triangular(double v, double w, double rho_j) {
    if (!(v > 0 && w > 0 && rho_j > 0 && std::isfinite(v) && std::isfinite(w) && std::isfinite(rho_j)))
      throw std::invalid_argument("physis::continuum::Flux::triangular: v > 0, w > 0 and ρ_j > 0");
    return {Kind::triangular, v, rho_j, w};
  }

  double q(double rho) const {
    if (kind == Kind::greenshields) return free_speed * rho * (1 - rho / jam_density);
    return std::min(free_speed * rho, wave_speed * (jam_density - rho));
  }
  // q′(ρ); at the triangular flux's kink ρ_c, the left derivative v.
  double dq(double rho) const {
    if (kind == Kind::greenshields) return free_speed * (1 - 2 * rho / jam_density);
    return rho <= critical() ? free_speed : -wave_speed;
  }
  // ρ_c, where q is greatest.
  double critical() const { return kind == Kind::greenshields ? 0.5 * jam_density : wave_speed * jam_density / (free_speed + wave_speed); }
  // q(ρ_c), per second.
  double capacity() const { return q(critical()); }
  // D(ρ) = q(min(ρ, ρ_c)), ρ read in [0, ρ_j].
  double demand(double rho) const { return q(std::min(std::max(rho, 0.0), critical())); }
  // S(ρ) = q(max(ρ, ρ_c)), ρ read in [0, ρ_j].
  double supply(double rho) const { return q(std::max(std::min(rho, jam_density), critical())); }
  // max |q′| over [0, ρ_j]: v (Greenshields), max(v, w) (triangular).
  double max_speed() const { return kind == Kind::greenshields ? free_speed : std::max(free_speed, wave_speed); }

  void validate(const char* where) const {
    const bool ok = free_speed > 0 && jam_density > 0 && std::isfinite(free_speed) && std::isfinite(jam_density) &&
                    (kind == Kind::greenshields || (wave_speed > 0 && std::isfinite(wave_speed)));
    if (!ok) throw std::invalid_argument(std::string(where) + ": a flux with v > 0, ρ_j > 0 (and w > 0 if triangular)");
  }
};

// Godunov's flux min(D(ρ_l), S(ρ_r)): the flux q(ρ(0)) of the Riemann solution at x/t = 0.
inline double godunov_flux(const Flux& f, double rho_left, double rho_right) { return std::min(f.demand(rho_left), f.supply(rho_right)); }

// The entropy solution of the Riemann problem (ρ_l for x < 0, ρ_r for x > 0) at x/t = ξ. A shock at
// s = (q(ρ_r) − q(ρ_l))/(ρ_r − ρ_l) when ρ_l < ρ_r (ρ_l for ξ < s, ρ_r from s on); a rarefaction when
// ρ_l > ρ_r: ρ_l for ξ ≤ q′(ρ_l), ρ_r for ξ ≥ q′(ρ_r), and between them the density whose wave speed
// is ξ, (ρ_j/2)(1 − ξ/v) for Greenshields and ρ_c for the triangular flux (its fan sits at the kink); on
// one linear branch of the triangular flux the fan degenerates to a discontinuity at speed v or −w.
inline double riemann(const Flux& f, double rho_left, double rho_right, double xi) {
  f.validate("physis::continuum::riemann");
  if (rho_left == rho_right) return rho_left;
  if (rho_left < rho_right) {
    const double s = (f.q(rho_right) - f.q(rho_left)) / (rho_right - rho_left);
    return xi < s ? rho_left : rho_right;
  }
  if (f.kind == Flux::Kind::greenshields) {
    if (xi <= f.dq(rho_left)) return rho_left;
    if (xi >= f.dq(rho_right)) return rho_right;
    return 0.5 * f.jam_density * (1 - xi / f.free_speed);
  }
  const double a = f.dq(rho_left), b = f.dq(rho_right);  // each v or −w, a ≤ b
  if (a == b) return xi < a ? rho_left : rho_right;
  if (xi <= a) return rho_left;
  if (xi >= b) return rho_right;
  return f.critical();
}

// ---- networks of links ---------------------------------------------------------------------------------------------------------

struct Link {
  std::size_t tail = 0, head = 0;  // nodes: the link runs from tail to head
  double length = 0;               // m
  std::size_t cells = 1;
  Flux flux;
};

class KinematicNetwork {
 public:
  enum class Rule : std::uint8_t { non_fifo, fifo };

  KinematicNetwork(std::size_t nodes, std::vector<Link> links) : links_(std::move(links)), nodes_(nodes), in_(nodes), out_(nodes), schedule_(links_.size()) {
    for (std::size_t l = 0; l < links_.size(); l++) {
      const Link& L = links_[l];
      if (L.tail >= nodes || L.head >= nodes) throw std::invalid_argument("physis::continuum::KinematicNetwork: a link's tail or head is not a node");
      if (!(L.length > 0 && std::isfinite(L.length)) || L.cells == 0) throw std::invalid_argument("physis::continuum::KinematicNetwork: a link of positive length and at least one cell");
      L.flux.validate("physis::continuum::KinematicNetwork");
      out_[L.tail].push_back(l);
      in_[L.head].push_back(l);
    }
  }

  // An entry: a node with no incoming link and one outgoing link, admitting demand(t) per second.
  void entry(std::size_t node, std::function<double(double t)> demand) {
    at(node, "entry");
    if (!in_[node].empty() || out_[node].size() != 1) throw std::invalid_argument("physis::continuum::KinematicNetwork::entry: an entry has no incoming link and one outgoing link");
    if (!demand) throw std::invalid_argument("physis::continuum::KinematicNetwork::entry: a demand");
    nodes_[node] = Node{Kind::entry, std::move(demand), {}, {}, Rule::non_fifo};
  }
  // An exit: a node with no outgoing link, taking at most supply(t) per second (empty: unlimited).
  void exit(std::size_t node, std::function<double(double t)> supply = {}) {
    at(node, "exit");
    if (!out_[node].empty() || in_[node].empty()) throw std::invalid_argument("physis::continuum::KinematicNetwork::exit: an exit has incoming links and no outgoing link");
    nodes_[node] = Node{Kind::exit, {}, std::move(supply), {}, Rule::non_fifo};
  }
  // A junction: beta[in][out] over the node's incoming and outgoing links in index order, rows summing to 1.
  void junction(std::size_t node, std::vector<std::vector<double>> beta, Rule rule = Rule::non_fifo) {
    at(node, "junction");
    if (beta.size() != in_[node].size()) throw std::invalid_argument("physis::continuum::KinematicNetwork::junction: one row of β for each incoming link");
    for (const auto& row : beta) {
      if (row.size() != out_[node].size()) throw std::invalid_argument("physis::continuum::KinematicNetwork::junction: one column of β for each outgoing link");
      double sum = 0;
      for (double b : row) {
        if (!(b >= 0 && b <= 1)) throw std::invalid_argument("physis::continuum::KinematicNetwork::junction: turning fractions in [0, 1]");
        sum += b;
      }
      if (std::fabs(sum - 1) > 1e-12) throw std::invalid_argument("physis::continuum::KinematicNetwork::junction: each incoming link's turning fractions sum to 1");
    }
    nodes_[node] = Node{Kind::junction, {}, {}, std::move(beta), rule};
  }
  // A link's flux law in time (a rule that opens a carriageway); if set, used instead of the link's flux.
  void schedule(std::size_t link, std::function<Flux(double t)> flux_at) {
    if (link >= links_.size()) throw std::invalid_argument("physis::continuum::KinematicNetwork::schedule: no such link");
    schedule_[link] = std::move(flux_at);
  }

  struct State {
    double t = 0;
    std::vector<std::vector<double>> density;  // per link, per cell from tail to head (per metre)
  };
  // The account of a run (numbers of walkers or vehicles). A budget whose fields are all zero takes, at
  // its first step, the content of that step's state as its initial content; each step then adds what
  // entered, left and was refused, and sets what is stored. initial + entered − left − stored = 0 to rounding.
  struct Budget {
    double entered = 0, left = 0, stored = 0, initial = 0;
    double refused = 0;  // demand an entry could not admit (lost, not queued)
    double residual() const { return initial + entered - left - stored; }
  };

  // The empty network at t0.
  State start(double t0 = 0) const {
    State s{t0, {}};
    s.density.reserve(links_.size());
    for (const Link& L : links_) s.density.emplace_back(L.cells, 0.0);
    return s;
  }

  // One step of Godunov's scheme; throws if h > cfl(s.t).
  State step(const State& s, double h, Budget* budget = nullptr) const {
    check(s, "physis::continuum::KinematicNetwork::step");
    if (!(h > 0 && std::isfinite(h))) throw std::invalid_argument("physis::continuum::KinematicNetwork::step: a step h > 0");
    if (h > cfl(s.t) * (1 + 1e-12)) throw std::invalid_argument("physis::continuum::KinematicNetwork::step: h exceeds the CFL bound min Δx / max|q′|");
    double entered = 0, left = 0, refused = 0;
    const auto F = flows(s, &entered, &left, &refused);
    State next{s.t + h, s.density};
    for (std::size_t l = 0; l < links_.size(); l++) {
      const double r = h / dx(l);
      for (std::size_t k = 0; k < links_[l].cells; k++) next.density[l][k] += r * (F[l][k] - F[l][k + 1]);
    }
    if (budget) {
      const bool fresh = budget->entered == 0 && budget->left == 0 && budget->stored == 0 && budget->initial == 0 && budget->refused == 0;
      if (fresh) budget->initial = content(s);
      budget->entered += h * entered;
      budget->left += h * left;
      budget->refused += h * refused;
      budget->stored = content(next);
    }
    return next;
  }

  // The CFL bound of the links' own fluxes, min over links of Δx / max|q′| (s); scheduled laws are not read.
  double cfl() const {
    double b = std::numeric_limits<double>::infinity();
    for (std::size_t l = 0; l < links_.size(); l++) b = std::min(b, dx(l) / links_[l].flux.max_speed());
    return b;
  }
  // The CFL bound of the laws in force at t (scheduled ones read at t): what step() requires of h.
  double cfl(double t) const {
    double b = std::numeric_limits<double>::infinity();
    for (std::size_t l = 0; l < links_.size(); l++) b = std::min(b, dx(l) / flux(l, t).max_speed());
    return b;
  }

  // The flow through each interface of each link (cells + 1 per link, from the tail), at the state s
  // and its time: the flows the next step from s uses, per second.
  std::vector<std::vector<double>> link_flows(const State& s) const {
    check(s, "physis::continuum::KinematicNetwork::link_flows");
    return flows(s, nullptr, nullptr, nullptr);
  }

  const std::vector<Link>& links() const { return links_; }
  std::size_t nodes() const { return nodes_.size(); }
  // A node's incoming and outgoing links, in index order (the rows and columns of its β).
  const std::vector<std::size_t>& incoming(std::size_t node) const { return in_.at(node); }
  const std::vector<std::size_t>& outgoing(std::size_t node) const { return out_.at(node); }
  // The flux law of a link in force at t.
  Flux flux(std::size_t link, double t) const {
    if (!schedule_.at(link)) return links_[link].flux;
    const Flux f = schedule_[link](t);
    f.validate("physis::continuum::KinematicNetwork: a scheduled flux");
    return f;
  }
  // A link's cell length Δx (m).
  double dx(std::size_t link) const { return links_.at(link).length / static_cast<double>(links_[link].cells); }
  // Σ ρ Δx over every cell: the number of walkers or vehicles on the network.
  double content(const State& s) const {
    double c = 0;
    for (std::size_t l = 0; l < links_.size(); l++)
      for (double r : s.density.at(l)) c += r * dx(l);
    return c;
  }

 private:
  enum class Kind : std::uint8_t { none, entry, exit, junction };
  struct Node {
    Kind kind = Kind::none;
    std::function<double(double)> demand, supply;
    std::vector<std::vector<double>> beta;
    Rule rule = Rule::non_fifo;
  };

  void at(std::size_t node, const char* what) const {
    if (node >= nodes_.size()) throw std::invalid_argument(std::string("physis::continuum::KinematicNetwork::") + what + ": no such node");
  }
  void check(const State& s, const char* where) const {
    bool ok = s.density.size() == links_.size();
    for (std::size_t l = 0; ok && l < links_.size(); l++) ok = s.density[l].size() == links_[l].cells;
    if (!ok) throw std::invalid_argument(std::string(where) + ": a state with one density for each cell of each link");
  }

  std::vector<std::vector<double>> flows(const State& s, double* entered, double* left, double* refused) const {
    const double t = s.t;
    std::vector<Flux> law(links_.size());
    for (std::size_t l = 0; l < links_.size(); l++) law[l] = flux(l, t);
    std::vector<std::vector<double>> F(links_.size());
    for (std::size_t l = 0; l < links_.size(); l++) {
      const auto& r = s.density[l];
      const std::size_t C = links_[l].cells;
      F[l].assign(C + 1, 0.0);
      for (std::size_t k = 1; k < C; k++) F[l][k] = godunov_flux(law[l], r[k - 1], r[k]);
    }
    for (std::size_t n = 0; n < nodes_.size(); n++) {
      const Node& N = nodes_[n];
      const auto& in = in_[n];
      const auto& out = out_[n];
      std::vector<double> D(in.size()), S(out.size());
      for (std::size_t i = 0; i < in.size(); i++) D[i] = law[in[i]].demand(s.density[in[i]].back());
      for (std::size_t j = 0; j < out.size(); j++) S[j] = law[out[j]].supply(s.density[out[j]].front());
      switch (N.kind) {
        case Kind::entry: {
          const double d = N.demand(t);
          if (!(d >= 0 && std::isfinite(d))) throw std::invalid_argument("physis::continuum::KinematicNetwork: an entry's demand is finite and ≥ 0");
          const double f = std::min(d, S[0]);
          F[out[0]][0] = f;
          if (entered) *entered += f;
          if (refused) *refused += d - f;
          break;
        }
        case Kind::exit: {
          double total = 0;
          for (double d : D) total += d;
          double S_exit = std::numeric_limits<double>::infinity();
          if (N.supply) {
            S_exit = N.supply(t);
            if (!(S_exit >= 0)) throw std::invalid_argument("physis::continuum::KinematicNetwork: an exit's supply is ≥ 0");
          }
          const double theta = total > S_exit ? S_exit / total : 1.0;
          for (std::size_t i = 0; i < in.size(); i++) {
            const double f = D[i] * theta;
            F[in[i]][links_[in[i]].cells] = f;
            if (left) *left += f;
          }
          break;
        }
        case Kind::junction:
        case Kind::none: {
          if (in.empty() || out.empty()) break;  // closed
          std::vector<std::vector<double>> ones;
          const std::vector<std::vector<double>>* beta = &N.beta;
          if (N.kind == Kind::none) {
            if (out.size() != 1) throw std::invalid_argument("physis::continuum::KinematicNetwork: a node with several outgoing links needs its turning fractions (junction)");
            ones.assign(in.size(), std::vector<double>{1.0});
            beta = &ones;
          }
          const auto& B = *beta;
          std::vector<double> theta(out.size(), 1.0);
          for (std::size_t j = 0; j < out.size(); j++) {
            double wanted = 0;
            for (std::size_t i = 0; i < in.size(); i++) wanted += D[i] * B[i][j];
            if (wanted > S[j]) theta[j] = S[j] / wanted;
          }
          std::vector<double> outflow(in.size(), 0.0), inflow(out.size(), 0.0);
          for (std::size_t i = 0; i < in.size(); i++) {
            double fifo = 1.0;
            if (N.rule == Rule::fifo)
              for (std::size_t j = 0; j < out.size(); j++)
                if (B[i][j] > 0) fifo = std::min(fifo, theta[j]);
            for (std::size_t j = 0; j < out.size(); j++) {
              if (B[i][j] == 0) continue;
              const double f = D[i] * B[i][j] * (N.rule == Rule::fifo ? fifo : theta[j]);
              outflow[i] += f;
              inflow[j] += f;
            }
          }
          for (std::size_t i = 0; i < in.size(); i++) F[in[i]][links_[in[i]].cells] = outflow[i];
          for (std::size_t j = 0; j < out.size(); j++) F[out[j]][0] = inflow[j];
          break;
        }
      }
    }
    return F;
  }

  std::vector<Link> links_;
  std::vector<Node> nodes_;
  std::vector<std::vector<std::size_t>> in_, out_;
  std::vector<std::function<Flux(double)>> schedule_;
};

}  // namespace physis::continuum
