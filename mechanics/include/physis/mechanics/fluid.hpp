// A planar tree of slender links in a fluid: the fluid's kinetic energy (added mass) and its drag, by strip
// theory (decision 0017; docs/theory/discrete-mechanics.md §8).
//
// The body. A slender part lies along a link's x from a point `from`, of length ℓ and height h (its extent
// across the plane of motion: a fish's depth, swimming in the horizontal plane). It is cut into strips of
// length ds; each strip moves with the velocity of its centre p, relative to the fluid w = ṗ − u, u the
// fluid's velocity there (a current), with components w_n across the part (along n = R(θ)(0, 1)) and w_t
// along it (along e = R(θ)(1, 0)).
//
// The reactive force: the fluid's kinetic energy. A plate of span h moving across itself sets the fluid
// moving, with the kinetic energy of the added mass ρπh²/4 per unit length (the two-dimensional plate:
// Newman 1977, Marine Hydrodynamics, §4.13; Lighthill 1960, J. Fluid Mech. 9, 305–317, for a slender
// swimmer); along itself, none. So
//
//   T_f = Σ ½ (ρπh²/4) ds w_n²,
//
// a term of the Lagrangian, written over any number type: its Euler–Lagrange equations give the added mass
// and the forces of its change, and the variational integrator moves it with the body.
//
// The resistive force: drag. Each strip feels the quadratic drag of the flow across it and the skin
// friction along it (Morison et al. 1950; Taylor 1952, Proc. R. Soc. A 214, 158–183, for strips),
//
//   f = −½ ρ C_n h ds |w_n| w_n n − ½ ρ C_f (2h) ds |w_t| w_t e,
//
// C_n of a plate across the flow (about 1 to 2 at a fish's Reynolds numbers) and C_f of friction on its two
// sides. |w| is regularised as √(w² + ε²), ε = 10⁻⁴ m/s, so that the force is smooth for jets. Its
// generalized force is Jᵀ f at each strip.
//
// Caveat. Strip theory treats each strip as if in a two-dimensional flow of its own: the flow about a body
// that the body has disturbed (its wake, the vortices a tail sheds) is not resolved, nor the momentum the
// trailing edge sheds apart from the strips' drag (Lighthill's large-amplitude theory). It is the model of
// resistive and reactive forces a swimming robot is first designed with, not a solution of the
// Navier–Stokes equations.

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <physis/mechanics/multibody.hpp>

namespace physis::mechanics {

struct SlenderPart {
  std::size_t link = 0;
  std::array<double, 2> from{};  // where it starts, in the link's frame (m)
  double length = 0;             // along the link's x (m)
  double height = 0;             // across the plane of motion (m)
};

struct Fluid {
  double density = 1000;          // ρ, kg/m³
  double normal_drag = 1.2;       // C_n
  double friction_drag = 0.02;    // C_f
  std::array<double, 2> current{};  // u, m/s, uniform
  int strips = 4;                 // a part's strips
  void validate() const {
    if (!(density > 0) || !(normal_drag >= 0) || !(friction_drag >= 0) || strips < 1) throw std::invalid_argument("physis::mechanics::Fluid: ρ > 0, drags ≥ 0, a strip at least");
  }
};

namespace fluid_detail {
template <class S>
S soft_abs(const S& w) {
  using std::sqrt;
  return sqrt(w * w + S(1e-8));
}
}  // namespace fluid_detail

// The fluid's kinetic energy, T_f = Σ ½ (ρπh²/4) ds w_n², over any number type.
template <std::size_t N, class S>
S fluid_kinetic(const PlanarTree<N>& tree, const std::vector<PlanarLinkState<S>>& st, const std::vector<SlenderPart>& parts, const Fluid& f) {
  using std::cos, std::sin;
  S T = S(0.0);
  for (const auto& b : parts) {
    const double ds = b.length / f.strips, ma = f.density * 3.141592653589793 * b.height * b.height / 4;
    const S c = cos(st[b.link].theta), s = sin(st[b.link].theta);
    for (int k = 0; k < f.strips; k++) {
      const std::array<double, 2> X{b.from[0] + (k + 0.5) * ds, b.from[1]};
      const auto pv = tree.point_velocity(st, b.link, X);
      const S wn = -(pv[0] - f.current[0]) * s + (pv[1] - f.current[1]) * c;
      T += wn * wn * (0.5 * ma * ds);
    }
  }
  return T;
}

// The drag's generalized force, over any number type.
template <std::size_t N, class S>
std::array<S, N> fluid_drag(const PlanarTree<N>& tree, const std::vector<PlanarLinkState<S>>& st, const std::vector<SlenderPart>& parts, const Fluid& f) {
  using std::cos, std::sin;
  std::array<S, N> Q;
  Q.fill(S(0.0));
  for (const auto& b : parts) {
    const double ds = b.length / f.strips;
    const double kn = 0.5 * f.density * f.normal_drag * b.height * ds, kt = 0.5 * f.density * f.friction_drag * 2 * b.height * ds;
    const S c = cos(st[b.link].theta), s = sin(st[b.link].theta);
    for (int k = 0; k < f.strips; k++) {
      const std::array<double, 2> X{b.from[0] + (k + 0.5) * ds, b.from[1]};
      const auto p = tree.point(st, b.link, X);
      const auto pv = tree.point_velocity(st, b.link, X);
      const S wx = pv[0] - f.current[0], wy = pv[1] - f.current[1];
      const S wn = -wx * s + wy * c, wt = wx * c + wy * s;
      const S fn = -kn * fluid_detail::soft_abs(wn) * wn, ft = -kt * fluid_detail::soft_abs(wt) * wt;
      const std::array<S, 2> force{fn * (-s) + ft * c, fn * c + ft * s};
      const auto g = tree.generalized_force(st, b.link, p, force);
      for (std::size_t i = 0; i < N; i++) Q[i] += g[i];
    }
  }
  return Q;
}

// The Lagrangian of a tree in a fluid: the tree's, and the fluid's kinetic energy.
template <std::size_t N>
struct TreeInFluid {
  const PlanarTree<N>* tree;
  const std::vector<SlenderPart>* parts;
  const Fluid* fluid;
  template <class S>
  S operator()(double, const std::array<S, N>& q, const std::array<S, N>& v) const {
    return tree->lagrangian(q, v) + fluid_kinetic<N>(*tree, tree->kinematics(q, v), *parts, *fluid);
  }
};

}  // namespace physis::mechanics
