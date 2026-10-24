// Gauge groups as the lattice gauge field needs them.
//
// A gauge group here is a value (it may carry an order N or a table) that knows its element
// type and the operations a lattice gauge theory uses (Montvay and Münster 1994, ch. 3;
// Creutz 1983, Quarks, Gluons and Lattices, ch. 7–9):
// - the product, inverse and identity;
// - the normalised real trace Re tr U / N of the defining representation, which enters the
//   Wilson action, or, for a finite group, a chosen real class function;
// - Haar-random elements, and symmetric proposals near an element for Metropolis;
// - a class angle in [0, 2π) to display an element by one number;
// - for the compact Lie groups, the "staple" type in which sums of elements live (ℂ for
//   U(1), the quaternions ℍ for SU(2), complex matrices for SU(N)), with Re tr(U Σ)/N.
//
// Adapters: ZN (ℤ_N ⊂ U(1), element k ↦ e^{2πik/N}), U1Group, SU2Group, SUNGroup (run-time
// N), and FiniteGroupAdapter over a morphe table (Re tr replaced by a class function, by
// default the normalised character of the regular representation, δ_{g,1}).

#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <string>
#include <vector>

#include <calx/random/stream.hpp>
#include <morphe/algebra/finite_group.hpp>
#include <physis/symmetry/su2.hpp>
#include <physis/symmetry/sun.hpp>
#include <physis/symmetry/u1.hpp>

namespace physis::gauge {

using calx::random::Stream;

template <class G>
concept GaugeGroup = requires(const G& g, const typename G::Element& a, const typename G::Element& b, Stream& z, double eps) {
  { g.identity() } -> std::convertible_to<typename G::Element>;
  { g.mul(a, b) } -> std::convertible_to<typename G::Element>;
  { g.inv(a) } -> std::convertible_to<typename G::Element>;
  { g.re_trace(a) } -> std::convertible_to<double>;
  { g.random(z) } -> std::convertible_to<typename G::Element>;
  { g.near(a, eps, z) } -> std::convertible_to<typename G::Element>;
  { g.class_angle(a) } -> std::convertible_to<double>;
  { g.name() } -> std::convertible_to<std::string>;
};

// A group whose Wilson action is linear in each link through a staple: heat bath and
// overrelaxation apply.
template <class G>
concept StapleGroup = GaugeGroup<G> && requires(const G& g, const typename G::Element& a, typename G::Staple& s) {
  { g.zero_staple() } -> std::convertible_to<typename G::Staple>;
  { g.add(s, a) };
  { g.re_trace_with(a, s) } -> std::convertible_to<double>;
};

// A finite group whose elements are 0, …, order() − 1: the heat bath enumerates them.
template <class G>
concept FiniteGaugeGroup = GaugeGroup<G> && std::same_as<typename G::Element, std::uint32_t> && requires(const G& g) {
  { g.order() } -> std::convertible_to<std::size_t>;
};

constexpr double two_pi = 2 * std::numbers::pi;

inline double angle_0_2pi(double t) {
  t = std::fmod(t, two_pi);
  return t < 0 ? t + two_pi : t;
}

// ℤ_N, the element k standing for e^{2πik/N}.
struct ZN {
  using Element = std::uint32_t;
  std::uint32_t n = 2;
  Element identity() const { return 0; }
  Element mul(Element a, Element b) const { return (a + b) % n; }
  Element inv(Element a) const { return (n - a) % n; }
  double re_trace(Element a) const { return std::cos(two_pi * a / n); }
  Element random(Stream& z) const { return static_cast<Element>(z.below(n)); }
  // A step of ±1 (or any element for N = 2): symmetric.
  Element near(Element a, double, Stream& z) const { return n == 2 ? (a + 1) % 2 : (z.uniform() < 0.5 ? (a + 1) % n : (a + n - 1) % n); }
  double class_angle(Element a) const { return two_pi * a / n; }
  std::string name() const { return "Z" + std::to_string(n); }
  std::size_t order() const { return n; }
};

struct U1Group {
  using Element = symmetry::U1;
  using Staple = std::complex<double>;
  Element identity() const { return Element::identity(); }
  Element mul(const Element& a, const Element& b) const { return a * b; }
  Element inv(const Element& a) const { return a.inverse(); }
  double re_trace(const Element& a) const { return a.re_trace(); }
  Element random(Stream& z) const { return symmetry::haar_u1(z); }
  Element near(const Element& a, double eps, Stream& z) const { return Element::from_angle(a.theta + eps * (2 * z.uniform() - 1)); }
  double class_angle(const Element& a) const { return angle_0_2pi(a.theta); }
  std::string name() const { return "U(1)"; }
  Staple zero_staple() const { return 0; }
  void add(Staple& s, const Element& a) const { s += a.value(); }
  double re_trace_with(const Element& a, const Staple& s) const { return (a.value() * s).real(); }
};

struct SU2Group {
  using Element = symmetry::SU2;
  using Staple = symmetry::SU2;  // a real multiple of an element of SU(2): a quaternion
  Element identity() const { return Element::identity(); }
  Element mul(const Element& a, const Element& b) const { return a * b; }
  Element inv(const Element& a) const { return a.inverse(); }
  double re_trace(const Element& a) const { return a.re_trace(); }
  Element random(Stream& z) const { return symmetry::haar_su2(z); }
  // R a with R = exp of a uniform vector in the ball of radius eps of su(2): symmetric.
  Element near(const Element& a, double eps, Stream& z) const {
    symmetry::Vector3 x;
    do {
      for (auto& c : x) c = 2 * z.uniform() - 1;
    } while (x[0] * x[0] + x[1] * x[1] + x[2] * x[2] > 1);
    for (auto& c : x) c *= eps;
    return (symmetry::exp_su2(x) * a).reprojected();
  }
  double class_angle(const Element& a) const { return angle_0_2pi(a.class_angle()); }
  std::string name() const { return "SU(2)"; }
  Staple zero_staple() const { return {0, 0, 0, 0}; }
  void add(Staple& s, const Element& a) const { s += a; }
  // Re tr(U Σ)/2 = (U Σ)_0 for quaternions.
  double re_trace_with(const Element& a, const Staple& s) const { return (a * s).a0; }
};

struct SUNGroup {
  using Element = symmetry::SUN;
  using Staple = symmetry::SUN;  // a complex matrix
  std::size_t n = 3;
  Element identity() const { return Element::identity(n); }
  Element mul(const Element& a, const Element& b) const { return a * b; }
  Element inv(const Element& a) const { return a.inverse(); }
  double re_trace(const Element& a) const { return a.re_trace(); }
  Element random(Stream& z) const { return symmetry::haar_sun(n, z); }
  // R a with R = exp(i ε Σ x_a t_a), x uniform in the cube: symmetric since x and -x are equally likely.
  Element near(const Element& a, double eps, Stream& z) const {
    std::vector<double> x(n * n - 1);
    for (auto& c : x) c = eps * (2 * z.uniform() - 1);
    return (symmetry::exp_sun(n, x) * a).reunitarised();
  }
  // 2 arccos of Re tr U / N, in [0, 2π]: for N = 2 the class angle.
  double class_angle(const Element& a) const { return 2 * std::acos(std::clamp(a.re_trace(), -1.0, 1.0)); }
  std::string name() const { return "SU(" + std::to_string(n) + ")"; }
  Staple zero_staple() const { return Element(n); }
  void add(Staple& s, const Element& a) const { s += a; }
  double re_trace_with(const Element& a, const Staple& s) const { return (a * s).re_trace(); }
};

// A finite group from its table; re_trace is a real class function, by default δ_{g,1}.
struct FiniteGroupAdapter {
  using Element = std::uint32_t;
  const morphe::algebra::FiniteGroup* group = nullptr;
  std::vector<double> weight;  // per element; empty: 1 at the identity, 0 elsewhere
  Element identity() const { return 0; }
  Element mul(Element a, Element b) const { return group->mul(a, b); }
  Element inv(Element a) const { return group->inv(a); }
  double re_trace(Element a) const { return weight.empty() ? (a == 0 ? 1.0 : 0.0) : weight[a]; }
  Element random(Stream& z) const { return static_cast<Element>(z.below(group->order())); }
  // Any other element, uniformly: symmetric.
  Element near(Element a, double, Stream& z) const {
    const auto n = group->order();
    if (n < 2) return a;
    return static_cast<Element>((a + 1 + z.below(n - 1)) % n);
  }
  double class_angle(Element a) const { return two_pi * a / static_cast<double>(group->order()); }
  std::string name() const { return group->name; }
  std::size_t order() const { return group->order(); }
};

}  // namespace physis::gauge
