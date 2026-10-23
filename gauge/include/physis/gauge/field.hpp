// Lattice gauge fields on a CW complex.
//
// A lattice gauge field (Wilson 1974, Phys. Rev. D 10, 2445) on a CW complex X with group G is
// a function U: 1-cells → G, the parallel transport along each edge from its tail to its
// head; along an edge traversed backwards the transport is U_e⁻¹. The transport along a
// path is the product of the transports of its steps in order, and the holonomy of a 2-cell
// is the transport along its attaching word from its base vertex (physis decision 0005,
// §2.4): the lattice field strength, a discrete curvature. Nothing here assumes a hypercubic
// lattice: every 2-complex of morphe carries a gauge theory, its plaquettes the 2-cells.
//
// A gauge transformation g: vertices → G acts by U_e ↦ g_tail U_e g_head⁻¹; the transport
// along a path from u to v becomes g_u T g_v⁻¹, so the holonomy of a 2-cell is conjugated
// by g_base and its trace is invariant: the Wilson action
//
//   S[U] = β Σ_f (1 - Re tr U_f / N)
//
// and the traces of closed transports (Wilson loops) are gauge invariant (Montvay and
// Münster 1994, eqs. 3.1–3.10; for a finite group Re tr/N is replaced by a class function).
//
// The action depends on a link U_e through the 2-cells whose words contain e. When e occurs
// once in a word, U_f = A U_e^{±1} B and Re tr U_f = Re tr(U_e Σ) with the staple Σ = B A for
// a forward step and (B A)⁻¹ for a backward one; summed over such 2-cells the local action is
// linear in U_e, -β/N Re tr(U_e Σ_e), the form heat bath and overrelaxation need. When e
// occurs twice in one word (the one-vertex torus a b a⁻¹ b⁻¹) it is not, and only
// Metropolis, which recomputes the holonomies, applies; linear() says which.

#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include <morphe/topology/cw_complex.hpp>

#include "group.hpp"

namespace physis::gauge {

using morphe::topology::CWComplex;
using morphe::topology::Path;
using morphe::topology::Step;

template <GaugeGroup G>
class GaugeField {
 public:
  using Element = typename G::Element;

  // The cold field: every link the identity.
  GaugeField(G group, const CWComplex& X) : group_(std::move(group)), X_(&X), links_(X.edges.size(), group_.identity()) {
    occurrences_.resize(X.edges.size());
    faces_of_.resize(X.edges.size());
    linear_.assign(X.edges.size(), 1);
    for (std::size_t f = 0; f < X.faces.size(); f++) {
      for (std::size_t i = 0; i < X.faces[f].boundary.size(); i++) {
        const std::size_t e = X.faces[f].boundary[i].edge;
        if (e >= X.edges.size()) throw std::invalid_argument("physis::gauge::GaugeField: a word names an edge that does not exist");
        occurrences_[e].push_back({f, i});
        if (faces_of_[e].empty() || faces_of_[e].back() != f) faces_of_[e].push_back(f);
        else linear_[e] = 0;  // e occurs twice in the word of f
      }
    }
  }

  const G& group() const { return group_; }
  const CWComplex& complex() const { return *X_; }
  std::vector<Element>& links() { return links_; }
  const std::vector<Element>& links() const { return links_; }
  const Element& link(std::size_t e) const { return links_[e]; }
  void set(std::size_t e, Element u) { links_[e] = std::move(u); }
  // The hot field: every link Haar-random.
  void randomise(Stream& z) {
    for (auto& u : links_) u = group_.random(z);
  }
  // Whether the action is linear in U_e (e occurs at most once in each word).
  bool linear(std::size_t e) const { return linear_[e]; }
  // The 2-cells whose words contain e (each once), and every occurrence (face, position).
  const std::vector<std::size_t>& faces_of(std::size_t e) const { return faces_of_[e]; }
  const std::vector<std::pair<std::size_t, std::size_t>>& occurrences(std::size_t e) const { return occurrences_[e]; }

  // The transport along a path: the product of its steps in order.
  Element transport(const Path& p) const {
    Element t = group_.identity();
    for (const auto& s : p) t = group_.mul(t, s.sign > 0 ? links_[s.edge] : group_.inv(links_[s.edge]));
    return t;
  }
  // The holonomy of 2-cell f, from its base along its word.
  Element holonomy(std::size_t f) const { return transport(X_->faces[f].boundary); }
  // Re tr U_f / N of 2-cell f.
  double plaquette(std::size_t f) const { return group_.re_trace(holonomy(f)); }
  double mean_plaquette() const {
    if (X_->faces.empty()) return 1;
    double s = 0;
    for (std::size_t f = 0; f < X_->faces.size(); f++) s += plaquette(f);
    return s / static_cast<double>(X_->faces.size());
  }
  // The Wilson action β Σ_f (1 - Re tr U_f / N).
  double action(double beta) const {
    double s = 0;
    for (std::size_t f = 0; f < X_->faces.size(); f++) s += 1 - plaquette(f);
    return beta * s;
  }
  // Re tr / N of the transport around a closed path: a Wilson (or Polyakov) loop.
  double wilson_loop(const Path& p) const { return group_.re_trace(transport(p)); }

  // U_e ↦ g_tail U_e g_head⁻¹.
  void transform(const std::vector<Element>& g) {
    if (g.size() != X_->vertices) throw std::invalid_argument("physis::gauge::GaugeField::transform: one element per vertex");
    for (std::size_t e = 0; e < links_.size(); e++)
      links_[e] = group_.mul(group_.mul(g[X_->edges[e][0]], links_[e]), group_.inv(g[X_->edges[e][1]]));
  }

  // Σ_{f ∋ e} Re tr U_f / N over the distinct 2-cells containing e: the part of -S/β that
  // depends on U_e.
  double local_plaquettes(std::size_t e) const {
    double s = 0;
    for (std::size_t f : faces_of_[e]) s += plaquette(f);
    return s;
  }

  // The staple Σ_e, with Σ_{f ∋ e} Re tr U_f / N = re_trace_with(U_e, Σ_e), for a linear edge.
  auto staple(std::size_t e) const
    requires StapleGroup<G>
  {
    auto sum = group_.zero_staple();
    for (const auto& [f, i] : occurrences_[e]) {
      const auto& word = X_->faces[f].boundary;
      // A: the steps before i; B: the steps after; Re tr(A U^s B) = Re tr(U^s B A).
      Element A = group_.identity(), B = group_.identity();
      for (std::size_t k = 0; k < i; k++) A = group_.mul(A, word[k].sign > 0 ? links_[word[k].edge] : group_.inv(links_[word[k].edge]));
      for (std::size_t k = i + 1; k < word.size(); k++) B = group_.mul(B, word[k].sign > 0 ? links_[word[k].edge] : group_.inv(links_[word[k].edge]));
      const Element BA = group_.mul(B, A);
      group_.add(sum, word[i].sign > 0 ? BA : group_.inv(BA));
    }
    return sum;
  }

 private:
  G group_;
  const CWComplex* X_;
  std::vector<Element> links_;
  std::vector<std::vector<std::pair<std::size_t, std::size_t>>> occurrences_;
  std::vector<std::vector<std::size_t>> faces_of_;
  std::vector<char> linear_;
};

}  // namespace physis::gauge
