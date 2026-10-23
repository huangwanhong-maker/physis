// Observables of lattice gauge fields.
//
// - The mean plaquette (1/F) Σ_f Re tr U_f / N, the lattice energy density.
// - Wilson loops, Re tr / N of the transport around a closed path (Wilson 1974); a loop that
//   winds around a non-contractible cycle (a generator of π_1 of the complex, as morphe's
//   fundamental_group() gives them) is a Polyakov loop when the cycle is a direction of
//   Euclidean time (Polyakov 1978), and in general a holonomy that no gauge transformation
//   removes.
// - The topological charge of a U(1) field on a closed orientable surface (Lüscher 1982,
//   Commun. Math. Phys. 85, 39; the geometric definition): with an orientation o of the
//   2-cells and arg U_f in (-π, π],
//     Q = (1/2π) Σ_f o_f arg U_f.
//   Every edge is traversed once in each direction by the oriented 2-cells, so the
//   unwrapped angles sum to zero and Q = -Σ_f o_f n_f, n_f the winding by which arg U_f was
//   wrapped: an integer, the first Chern number of the U(1) bundle the field defines when
//   no plaquette angle is near π.

#pragma once

#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <vector>

#include "field.hpp"

namespace physis::gauge {

template <GaugeGroup G>
double mean_plaquette(const GaugeField<G>& U) { return U.mean_plaquette(); }

template <GaugeGroup G>
double polyakov_loop(const GaugeField<G>& U, const Path& loop) { return U.wilson_loop(loop); }

// The plaquette values Re tr U_f / N and class angles of the holonomies, per 2-cell.
template <GaugeGroup G>
std::vector<double> plaquettes(const GaugeField<G>& U) {
  std::vector<double> p(U.complex().faces.size());
  for (std::size_t f = 0; f < p.size(); f++) p[f] = U.plaquette(f);
  return p;
}
template <GaugeGroup G>
std::vector<double> plaquette_angles(const GaugeField<G>& U) {
  std::vector<double> p(U.complex().faces.size());
  for (std::size_t f = 0; f < p.size(); f++) p[f] = U.group().class_angle(U.holonomy(f));
  return p;
}
template <GaugeGroup G>
std::vector<double> link_angles(const GaugeField<G>& U) {
  std::vector<double> a(U.links().size());
  for (std::size_t e = 0; e < a.size(); e++) a[e] = U.group().class_angle(U.link(e));
  return a;
}

// Q = (1/2π) Σ_f o_f arg U_f for a U(1) field, given an orientation of the 2-cells (morphe's
// classify_surface().orientation on a closed orientable surface).
inline double topological_charge(const GaugeField<U1Group>& U, const std::vector<int>& orientation) {
  if (orientation.size() != U.complex().faces.size()) throw std::invalid_argument("physis::gauge::topological_charge: one sign per 2-cell");
  double q = 0;
  for (std::size_t f = 0; f < orientation.size(); f++) q += orientation[f] * U.holonomy(f).theta;
  return q / (2 * std::numbers::pi);
}

}  // namespace physis::gauge
