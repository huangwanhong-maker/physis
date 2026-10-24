// A Monte Carlo run of a lattice gauge theory on a CW complex, with the group chosen at run
// time: the composition of field.hpp, updates.hpp, observables.hpp and exact.hpp that a
// program reading a specification (an engine, a notebook) calls once.
//
// The run starts cold (or hot), makes `thermalise` sweeps, then `sweeps` measured sweeps;
// a sweep is one heat-bath sweep followed by `overrelax` overrelaxation sweeps (U(1),
// SU(2)), or one Metropolis sweep with `hits` proposals per link. After each measured sweep
// it records the mean plaquette and the Wilson loops asked for. The result carries the
// final configuration (link class angles, raw links, plaquette values and angles), the
// history, the estimates with their errors (field_theory statistics), the U(1) topological
// charge on a closed orientable surface, and the exact mean plaquette where it is known: on
// a closed orientable surface (Migdal–Witten) for ℤ_N, U(1) and SU(2), and on a disc for
// these and SU(3).

#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <morphe/topology/surfaces.hpp>
#include <physis/field_theory/statistics.hpp>

#include "exact.hpp"
#include "observables.hpp"
#include "updates.hpp"

namespace physis::gauge {

struct GroupSpec {
  std::string kind = "U1";  // "Z", "U1", "SU2" or "SU"
  std::uint32_t n = 0;      // N for ℤ_N and SU(N)
};

struct RunOptions {
  double beta = 1;
  std::size_t thermalise = 100, sweeps = 500;
  std::uint64_t seed = 1;
  std::string update = "heatbath";  // "heatbath" or "metropolis"
  std::size_t overrelax = 1;        // overrelaxation sweeps after each heat-bath sweep (U(1), SU(2))
  double step = 0.5;                // the Metropolis step
  int hits = 4;                     // Metropolis proposals per link
  bool hot = false;                 // start from Haar-random links
  std::vector<Path> loops;          // Wilson loops to measure
};

struct RunResult {
  std::string group;
  std::vector<double> link_angle;          // class angle of each link, [0, 2π)
  std::vector<std::vector<double>> links;  // raw links: [θ] for U(1), [k] for ℤ_N, [a0, a1, a2, a3] for SU(2), [Re tr/N, Im tr/N] for SU(N)
  std::vector<double> plaquettes;          // Re tr U_f / N of the final configuration
  std::vector<double> plaquette_angle;     // class angle of each holonomy
  std::vector<double> history;             // the mean plaquette after each measured sweep
  field_theory::Estimate plaquette;        // its mean with error and autocorrelation time
  std::vector<field_theory::Estimate> loops;
  std::optional<double> charge;            // U(1) on a closed orientable surface
  std::optional<double> exact_plaquette;   // where known
  double acceptance = 1;
};

// The exact mean plaquette for the complex and group, where known.
inline std::optional<double> exact_mean_plaquette(const CWComplex& X, const GroupSpec& g, double beta) {
  const auto S = morphe::topology::classify_surface(X);
  if (!S.surface || !S.orientable) return std::nullopt;
  const bool closed = S.closed, is_disc = !S.closed && S.genus == 0 && S.boundary_components == 1;
  if (!closed && !is_disc) return std::nullopt;
  std::optional<CharacterExpansion> E;
  if (g.kind == "U1") E = expansion_u1(beta);
  else if (g.kind == "SU2" || (g.kind == "SU" && g.n == 2)) E = expansion_su2(beta);
  else if (g.kind == "Z") E = expansion_zn(g.n, beta);
  else if (g.kind == "SU" && g.n == 3 && is_disc) return free_plaquette_sun(3, beta);
  if (!E) return std::nullopt;
  return closed ? migdal_witten(*E, X.euler_characteristic(), X.faces.size()).mean_plaquette : disc(*E, X.faces.size()).mean_plaquette;
}

namespace detail {

inline std::vector<double> raw(const ZN&, std::uint32_t k) { return {static_cast<double>(k)}; }
inline std::vector<double> raw(const U1Group&, const symmetry::U1& u) { return {u.theta}; }
inline std::vector<double> raw(const SU2Group&, const symmetry::SU2& u) { return {u.a0, u.a1, u.a2, u.a3}; }
inline std::vector<double> raw(const SUNGroup& g, const symmetry::SUN& u) {
  const auto t = u.trace();
  return {t.real() / static_cast<double>(g.n), t.imag() / static_cast<double>(g.n)};
}

template <GaugeGroup G>
RunResult run(const CWComplex& X, G group, const RunOptions& o, const GroupSpec& spec) {
  GaugeField<G> U(std::move(group), X);
  calx::random::Stream z(o.seed, 0x67617567u, 0, 0);  // "gaug"
  if (o.hot) U.randomise(z);
  RunResult R;
  R.group = U.group().name();
  SweepStatistics stats;
  auto sweep = [&]() {
    if (o.update == "metropolis") {
      stats += metropolis_sweep(U, o.beta, z, o.step, o.hits);
      return;
    }
    if constexpr (StapleGroup<G>) {
      stats += heat_bath_sweep(U, o.beta, z, o.step);
      if constexpr (std::is_same_v<G, U1Group> || std::is_same_v<G, SU2Group>)
        for (std::size_t k = 0; k < o.overrelax; k++) overrelaxation_sweep(U);
    } else if constexpr (FiniteGaugeGroup<G>) {
      stats += heat_bath_sweep(U, o.beta, z);
    } else {
      stats += metropolis_sweep(U, o.beta, z, o.step, o.hits);
    }
  };
  for (std::size_t i = 0; i < o.thermalise; i++) sweep();
  stats = {};
  std::vector<std::vector<double>> loop_series(o.loops.size());
  for (std::size_t i = 0; i < o.sweeps; i++) {
    sweep();
    R.history.push_back(U.mean_plaquette());
    for (std::size_t l = 0; l < o.loops.size(); l++) loop_series[l].push_back(U.wilson_loop(o.loops[l]));
  }
  R.plaquette = field_theory::estimate(R.history);
  for (const auto& s : loop_series) R.loops.push_back(field_theory::estimate(s));
  R.link_angle = link_angles(U);
  for (const auto& u : U.links()) R.links.push_back(raw(U.group(), u));
  R.plaquettes = plaquettes(U);
  R.plaquette_angle = plaquette_angles(U);
  R.acceptance = stats.acceptance();
  if constexpr (std::is_same_v<G, U1Group>) {
    const auto S = morphe::topology::classify_surface(X);
    if (S.surface && S.closed && S.orientable) R.charge = std::round(topological_charge(U, S.orientation));
  }
  R.exact_plaquette = exact_mean_plaquette(X, spec, o.beta);
  return R;
}

}  // namespace detail

// The run, for a group named at run time.
inline RunResult run(const CWComplex& X, const GroupSpec& g, const RunOptions& o) {
  if (g.kind == "Z") {
    if (g.n < 2) throw std::invalid_argument("physis::gauge::run: Z_N needs N >= 2");
    return detail::run(X, ZN{g.n}, o, g);
  }
  if (g.kind == "U1") return detail::run(X, U1Group{}, o, g);
  if (g.kind == "SU2" || (g.kind == "SU" && g.n == 2)) return detail::run(X, SU2Group{}, o, g);
  if (g.kind == "SU") {
    if (g.n < 2 || g.n > 8) throw std::invalid_argument("physis::gauge::run: SU(N) needs 2 <= N <= 8");
    return detail::run(X, SUNGroup{g.n}, o, g);
  }
  throw std::invalid_argument("physis::gauge::run: unknown group \"" + g.kind + "\" (Z, U1, SU2, SU)");
}

}  // namespace physis::gauge
