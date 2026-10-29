// Fields given by formulas, and uniform fields (decision 0012; docs/theory/fields.md §3).
//
// An analytic field is a function fn(t, x) on its domain. Written once as a template over the number
// type of x (a generic lambda that calls sqrt, sin, exp unqualified after `using std::sqrt;`), it is
// evaluated on doubles for its values and on calx's dual numbers for its spatial derivative, exact to
// rounding (forward differentiation; Griewank and Walther 2008, ch. 3): a scalar's gradient, a vector's
// Jacobian, and from it the divergence and curl. A formula that does not take dual numbers is still a
// field; it has no derivative of its own, and says so by lacking the capability.
//
// A uniform field has one value on its domain (the interplanetary field over the Moon, a constant
// gravity of a small world).

#pragma once

#include <array>
#include <concepts>
#include <type_traits>
#include <utility>

#include <calx/differentiation/jet.hpp>
#include <physis/fields/domains.hpp>
#include <physis/fields/field.hpp>

namespace physis::fields {

namespace analytic_detail {

template <class P>
struct dual_point;
template <>
struct dual_point<double> {
  using type = calx::differentiation::Dual<1>;
  static type of(double x) { return calx::differentiation::variables<1, 1>(std::array<double, 1>{x})[0]; }
};
template <>
struct dual_point<Vec3> {
  using type = std::array<calx::differentiation::Dual<3>, 3>;
  static type of(const Vec3& x) { return calx::differentiation::variables<3, 1>(x); }
};

// The value and derivative of a field's formula at a point, read off its dual numbers.
template <class V, class P, class R>
Derivative<V, P> derivative_from(const R& r) {
  if constexpr (std::is_same_v<V, double> && std::is_same_v<P, double>) {
    return r.d[0];
  } else if constexpr (std::is_same_v<V, double> && std::is_same_v<P, Vec3>) {
    return Vec3{r.d[0], r.d[1], r.d[2]};
  } else {
    static_assert(std::is_same_v<V, Vec3> && std::is_same_v<P, Vec3>, "physis::fields: derivatives of scalars and of vectors in space");
    Mat3 j{};
    for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++) j[i][k] = r[i].d[k];
    return j;
  }
}

}  // namespace analytic_detail

// A formula fn(t, x) on a domain D. If fn also takes the dual numbers of x it is differentiable; a
// vector formula then returns std::array<S, 3> with S the number type of x's components, so that its
// three components have one type whatever constants it mixes in.
template <class V, class P, class Fn, class D = Everywhere>
class Analytic {
 public:
  using value_type = V;
  using point_type = P;

  Analytic(Descriptor d, Fn fn, D domain = {}) : d_(std::move(d)), fn_(std::move(fn)), domain_(std::move(domain)) {
    d_.representation.kind = RepresentationKind::analytic;
  }

  const Descriptor& descriptor() const { return d_; }
  const D& domain() const { return domain_; }

  Evaluation<V> at(const Query<P>& q) const {
    const Cover c = domain_.covers(q);
    if (!c.ok()) return refused<V>(c.quality, c.reason);
    return valid(static_cast<V>(fn_(q.t, q.x)));
  }

  Evaluation<Derivative<V, P>> jacobian(const Query<P>& q) const
    requires std::invocable<const Fn&, double, const typename analytic_detail::dual_point<P>::type&>
  {
    const Cover c = domain_.covers(q);
    if (!c.ok()) return refused<Derivative<V, P>>(c.quality, c.reason);
    const auto r = fn_(q.t, analytic_detail::dual_point<P>::of(q.x));
    return valid(analytic_detail::derivative_from<V, P>(r));
  }

 private:
  Descriptor d_;
  Fn fn_;
  D domain_;
};

template <class V, class P = Vec3, class Fn, class D = Everywhere>
Analytic<V, P, Fn, D> analytic(Descriptor d, Fn fn, D domain = {}) {
  return Analytic<V, P, Fn, D>(std::move(d), std::move(fn), std::move(domain));
}

// One value on a domain; its derivative is zero.
template <class V, class P = Vec3, class D = Everywhere>
class Uniform {
 public:
  using value_type = V;
  using point_type = P;

  Uniform(Descriptor d, V value, D domain = {}) : d_(std::move(d)), value_(value), domain_(std::move(domain)) {
    d_.representation.kind = RepresentationKind::uniform;
  }
  const Descriptor& descriptor() const { return d_; }
  Evaluation<V> at(const Query<P>& q) const {
    const Cover c = domain_.covers(q);
    if (!c.ok()) return refused<V>(c.quality, c.reason);
    return valid(value_);
  }
  Evaluation<Derivative<V, P>> jacobian(const Query<P>& q) const {
    const Cover c = domain_.covers(q);
    if (!c.ok()) return refused<Derivative<V, P>>(c.quality, c.reason);
    return valid(Derivative<V, P>{});
  }

 private:
  Descriptor d_;
  V value_;
  D domain_;
};

// The spatial derivative of a field that has no exact one, by central differences with step h along
// each axis: second-order in h, with a truncation error of order h² times the field's third
// derivative and a rounding error of order ε|φ|/h. Its answer carries the note that says so.
template <PointEvaluable F>
  requires std::is_same_v<typename F::point_type, Vec3>
Evaluation<Derivative<typename F::value_type, Vec3>> central_differences(const F& f, const Query<Vec3>& q, double h) {
  using V = typename F::value_type;
  using DV = Derivative<V, Vec3>;
  DV out{};
  for (int k = 0; k < 3; k++) {
    Query<Vec3> plus = q, minus = q;
    plus.x[k] += h;
    minus.x[k] -= h;
    const auto a = f.at(plus), b = f.at(minus);
    if (!a.ok()) return refused<DV>(a.quality, a.reason);
    if (!b.ok()) return refused<DV>(b.quality, b.reason);
    if constexpr (std::is_same_v<V, double>) {
      out[k] = (a.value - b.value) / (2 * h);
    } else {
      for (int i = 0; i < 3; i++) out[i][k] = (a.value[i] - b.value[i]) / (2 * h);
    }
  }
  return valid(out, "by central differences");
}

}  // namespace physis::fields
