// Group field theories: invariant field theories whose indices are group elements, and whose
// Feynman amplitudes are lattice gauge theories on the graphs' 2-complexes (physis decision
// 0006, §2.1, stages 1–3).
//
// A group field theory of rank d (Boulatov 1992, Mod. Phys. Lett. A 7, 1629; Ooguri 1992, Mod.
// Phys. Lett. A 7, 2799; Freidel 2005, Int. J. Theor. Phys. 44, 1769; Oriti 2009, in Approaches
// to Quantum Gravity, ch. 16) has a complex field φ(g_1, …, g_d) on G^d. Its kinetic term is the
// projection onto gauge-invariant fields, φ(g_1 h, …, g_d h) = φ(g_1, …, g_d), so that each
// propagator carries a group element h_e and identifies the strands g_i ↔ g_i h_e; its
// interactions are bubbles, as in invariant theories (field_theory/invariant.hpp), or, in the
// coloured simplicial models, the d + 1 fields of a d-simplex (tensor_models.hpp). Integrating
// the strands, every face (a closed strand, of colours (0, i) in a tensorial model, of any two
// colours in a coloured one) contributes δ(Π h_e^{±1}) of the holonomy along it, and the
// amplitude of a Feynman graph is the BF state sum on its 2-complex (bf.hpp):
//
//   A(Γ) = ∫ Π_e dh_e Π_f δ(h_f) = |G|^{F − E} #flat(Γ)   for a finite group G,
//
// a spin foam of flat connections: the 2-complex has a vertex for each bubble (tensorial) or
// each simplex (coloured), an edge for each propagator, a 2-cell for each face. Without the
// gauge projection every face carries δ(1) = |G| and A = |G|^F: the U(N)^d tensor model with
// N = |G|, so that the invariant theories are the group field theories of the trivial
// holonomy. With G of order 1 every amplitude is 1, and the series counts graphs.
//
// For a coloured graph that encodes a 3-manifold M (D = 3), the 2-complex is the dual
// 2-skeleton, π1 of it is π1(M), and A = |G|^{F − E + V − 1} |Hom(π1 M, G)|: for the melonic
// graphs, spheres, F = 3p + 3, E = 4p, V = 2p and A = |G|^{p+2}, the divergence of the
// Boulatov model growing with the order (Freidel and Louapre 2004, Nucl. Phys. B 662, 279).
//
// The free energy at order n is (1/n!) Σ over connected matchings of A, from the isomorphism
// classes of invariant.hpp for tensorial models (P <= 8) and from the coloured graphs with
// σ_0 = id, weighted 1/p!, for coloured ones (tensor_models.hpp).

#pragma once

#include <cstddef>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

#include <physis/field_theory/invariant.hpp>

#include "bf.hpp"
#include "tensor_models.hpp"

namespace physis::quantum_gravity {

struct GFTModel {
  enum class Kind { coloured, tensorial };
  std::string name;
  Kind kind = Kind::coloured;
  std::size_t dimension = 3;              // D: the dimension the Feynman graphs triangulate
  field_theory::InvariantTheory theory;   // a tensorial model's bubbles (rank D)
};

// Boulatov's model of three-dimensional gravity and Ooguri's of four-dimensional BF theory, in
// their coloured form (D + 1 fields).
inline GFTModel boulatov() { return {"Boulatov (coloured, D = 3)", GFTModel::Kind::coloured, 3, {}}; }
inline GFTModel ooguri() { return {"Ooguri (coloured, D = 4)", GFTModel::Kind::coloured, 4, {}}; }
// The quartic melonic group field theory of rank d, one bubble of each colour (or of colour 1).
inline GFTModel quartic_melonic_gft(std::size_t rank, bool every_colour = true) {
  return {"quartic melonic GFT, rank " + std::to_string(rank), GFTModel::Kind::tensorial, rank, field_theory::quartic_melonic_model(rank, every_colour)};
}
// "boulatov", "ooguri", "melonic" (rank d), "melonic-1" (rank d, colour 1 only).
inline GFTModel gft_by_name(const std::string& name, std::size_t rank = 3) {
  if (name == "boulatov") return boulatov();
  if (name == "ooguri") return ooguri();
  if (name == "melonic") return quartic_melonic_gft(rank, true);
  if (name == "melonic-1") return quartic_melonic_gft(rank, false);
  throw std::invalid_argument("physis::quantum_gravity::gft_by_name: no model \"" + name + "\"");
}

// The 2-complex of a tensorial model's Feynman graph: a vertex for each bubble copy, an edge for
// each propagator (from the copy of its black vertex to that of its white one), a 2-cell for each
// face of colours (0, i), its word the propagators met along it.
inline CWComplex tensorial_two_complex(const field_theory::Expansion& E, const std::vector<std::size_t>& sigma0, std::size_t rank) {
  const std::size_t P = sigma0.size();
  if (P != E.black_bubble.size()) throw std::invalid_argument("physis::quantum_gravity::tensorial_two_complex: σ₀ of the wrong size");
  CWComplex X;
  X.name = "Feynman 2-complex";
  X.vertices = E.bubbles.size();
  for (std::size_t b = 0; b < P; b++) {
    X.edges.push_back({E.black_bubble[b], E.white_bubble[sigma0[b]]});
    X.edge_names.push_back("e" + std::to_string(b + 1));
  }
  // The bubbles' colours as global maps white → black.
  std::vector<std::size_t> offset(E.bubbles.size(), 0);
  for (std::size_t k = 1; k < E.bubbles.size(); k++) offset[k] = offset[k - 1] + E.bubbles[k - 1].p;
  for (std::size_t i = 0; i < rank; i++) {
    std::vector<std::size_t> inv(P);
    for (std::size_t k = 0; k < E.bubbles.size(); k++)
      for (std::size_t b = 0; b < E.bubbles[k].p; b++) inv[offset[k] + E.bubbles[k].sigma[i][b]] = offset[k] + b;
    std::vector<char> seen(P, 0);
    for (std::size_t b0 = 0; b0 < P; b0++) {
      if (seen[b0]) continue;
      morphe::topology::Face f{E.black_bubble[b0], {}};
      for (std::size_t b = b0; !seen[b]; b = inv[sigma0[b]]) {
        seen[b] = 1;
        f.boundary.push_back({b, 1});  // along the propagator from black b to its white
      }
      X.faces.push_back(f);
    }
  }
  X.validate();
  return X;
}

// The amplitude of a Feynman 2-complex with a finite group: BF theory (with the gauge
// projection) or |G|^F (without it).
inline Rational gft_amplitude(const CWComplex& X, const FiniteGroup& G, bool gauge_projection = true) {
  if (gauge_projection) return bf_partition_function(X, G);
  return Rational(morphe::arithmetic::pow(Integer(static_cast<long long>(G.order())), X.faces.size()));
}

// The free energy's coefficient at order n: of (−λ)^n (every interaction with the coupling λ) for
// a tensorial model, of (λλ̄)^n for a coloured one; connected graphs only.
inline Rational gft_free_energy(const GFTModel& M, std::size_t n, const FiniteGroup& G, bool gauge_projection = true) {
  Rational total = 0;
  if (M.kind == GFTModel::Kind::coloured) {
    for_each_coloured_graph(M.dimension, n, [&](const ColouredGraph& Gamma) {
      if (Gamma.connected()) total += gft_amplitude(Gamma.two_complex(), G, gauge_projection);
    });
    return total / Rational(morphe::arithmetic::factorial(n));
  }
  // Every distribution of the n bubbles among the interactions, weighted 1 / Π n_B!.
  const std::size_t types = M.theory.interactions.size();
  std::vector<std::size_t> counts(types, 0);
  std::function<void(std::size_t, std::size_t)> go = [&](std::size_t t, std::size_t left) {
    if (t + 1 == types) {
      counts[t] = left;
      const auto E = field_theory::expand(M.theory, counts, true);
      if (!E.classes_computed) throw std::invalid_argument("physis::quantum_gravity::gft_free_energy: more than 8 propagators");
      Integer weight = 1;
      for (auto c : counts) weight *= morphe::arithmetic::factorial(c);
      for (const auto& C : E.classes)
        if (C.connected)
          total += gft_amplitude(tensorial_two_complex(E, C.sigma0, M.theory.rank), G, gauge_projection) *
                   Rational(Integer(static_cast<unsigned long long>(C.labelled)), weight);
      return;
    }
    for (std::size_t c = 0; c <= left; c++) {
      counts[t] = c;
      go(t + 1, left - c);
    }
  };
  go(0, n);
  return total;
}

}  // namespace physis::quantum_gravity
