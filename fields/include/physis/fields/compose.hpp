// Fields made from fields (decision 0012; docs/theory/fields.md §4).
//
//   reduce(f, ρ)        ρ ∘ f: a field of the reduction's fibre on f's base;
//   sum(a, b)           a + b on the intersection of their bases, a composite whose parts keep their
//                       provenance and whose own is the most derived of theirs;
//   pullback(f, m)      f ∘ m along a map of bases m: P → Q (the standard atmosphere, a field on height,
//                       carried to space by the height of each point above the geoid);
//   in_frame(f, M)      the same field described in another frame, M(t) = (R, b) the motion from f's
//                       frame to it, x′ = R x + b: points carried back, vectors and derivatives turned
//                       (a passive change of description: the world is unchanged, formalism §6.1);
//   restrict(f, D)      f on the part of its base that D covers.
//
// Each keeps what its parts can do where the composition allows it: a reduced field samples a sphere
// by its field's fast path; a field in another frame has its field's derivative, turned.

#pragma once

#include <algorithm>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <physis/fields/domains.hpp>
#include <physis/fields/field.hpp>
#include <physis/fields/reductions.hpp>
#include <physis/fields/sampling.hpp>

namespace physis::fields {

// The order of derivation of provenance, from the record of a measurement to a value without a law behind
// it: physical, archive, model, assimilated (a law joined to records), synthetic (statistics alone). A
// composite is of the most derived kind among its parts. The order is of kinds of evidence, not of
// accuracy: an assimilated estimate is usually nearer the truth than the model alone.
constexpr Provenance most_derived(Provenance a, Provenance b) { return static_cast<int>(a) >= static_cast<int>(b) ? a : b; }

// ---- reduce -------------------------------------------------------------------------------------------------

template <PointEvaluable F, Reduction R>
  requires std::is_same_v<typename F::value_type, typename R::input> && std::is_same_v<typename F::point_type, typename R::point_type>
class Reduced {
 public:
  using value_type = typename R::output;
  using point_type = typename F::point_type;
  // A reduced field samples a sphere as fast as its field does.
  static constexpr std::uint32_t extra_capabilities = capabilities_of<F>().bits & static_cast<std::uint32_t>(Capability::sphere_grid);

  Reduced(F f, R r) : f_(std::move(f)), r_(std::move(r)), d_(described(f_)) {
    d_.name += "." + std::string(r_.name());
    d_.title += ", " + std::string(r_.name());
    d_.fibre = r_.fibre(d_.fibre);
    d_.capabilities = {};
  }

  const Descriptor& descriptor() const { return d_; }
  const F& field() const { return f_; }
  const R& reduction() const { return r_; }

  Evaluation<value_type> at(const Query<point_type>& q) const {
    const auto a = f_.at(q);
    if (!a.ok()) return like<value_type>(a);
    return valid(static_cast<value_type>(r_(q, a.value)), a.note);
  }
  void at_many(std::span<const Query<point_type>> q, std::span<Evaluation<value_type>> out) const
    requires BatchEvaluable<F>
  {
    std::vector<Evaluation<typename F::value_type>> a(q.size());
    f_.at_many(q, std::span<Evaluation<typename F::value_type>>(a));
    for (std::size_t i = 0; i < q.size(); i++) out[i] = a[i].ok() ? valid(static_cast<value_type>(r_(q[i], a[i].value)), a[i].note) : like<value_type>(a[i]);
  }
  template <Sampling G>
  Samples<value_type> sample(const G& g, double t) const {
    const Samples<typename F::value_type> s = fields::sample(f_, g, t);
    Samples<value_type> out;
    out.values.resize(s.size());
    out.quality = s.quality;
    out.counts = s.counts;
    out.reasons = s.reasons;
    out.notes = s.notes;
    out.why = s.why;
    for (std::size_t i = 0; i < s.size(); i++)
      if (s.quality[i] == Quality::valid) out.values[i] = static_cast<value_type>(r_({t, g.point(i)}, s.values[i]));
    return out;
  }

 private:
  F f_;
  R r_;
  Descriptor d_;
};

template <PointEvaluable F, Reduction R>
Reduced<F, R> reduce(F f, R r) { return Reduced<F, R>(std::move(f), std::move(r)); }

// ---- sum ----------------------------------------------------------------------------------------------------

template <PointEvaluable A, PointEvaluable B>
  requires std::is_same_v<typename A::value_type, typename B::value_type> && std::is_same_v<typename A::point_type, typename B::point_type>
class Sum {
 public:
  using value_type = typename A::value_type;
  using point_type = typename A::point_type;

  // d names the composite; its fibre, parts, provenance and representation come from a and b, whose
  // fibres must have one dimension.
  Sum(Descriptor d, A a, B b) : a_(std::move(a)), b_(std::move(b)), d_(std::move(d)) {
    const Descriptor da = described(a_), db = described(b_);
    if (!(da.fibre.dimension == db.fibre.dimension)) throw std::invalid_argument("physis::fields::sum: parts of different dimensions");
    d_.fibre = da.fibre;
    d_.representation.kind = RepresentationKind::composite;
    d_.parts.clear();
    for (const Descriptor* x : {&da, &db}) {
      if (x->parts.empty()) d_.parts.push_back({x->name, x->provenance});
      else d_.parts.insert(d_.parts.end(), x->parts.begin(), x->parts.end());
    }
    d_.provenance = most_derived(da.provenance, db.provenance);
    // It reads a sphere as fast as its parts both do (sample() below reads each by its own method).
    d_.capabilities = {};
    if (da.capabilities.has(Capability::sphere_grid) && db.capabilities.has(Capability::sphere_grid)) d_.capabilities.add(Capability::sphere_grid);
  }

  const Descriptor& descriptor() const { return d_; }

  Evaluation<value_type> at(const Query<point_type>& q) const {
    const auto x = a_.at(q);
    if (!x.ok()) return x;
    const auto y = b_.at(q);
    if (!y.ok()) return y;
    return valid(add(x.value, y.value), x.note.empty() ? y.note : x.note);
  }
  template <Sampling G>
  Samples<value_type> sample(const G& g, double t) const
    requires std::is_same_v<point_type, Vec3>
  {
    const auto x = fields::sample(a_, g, t), y = fields::sample(b_, g, t);
    std::vector<Evaluation<value_type>> e(x.size());
    for (std::size_t i = 0; i < x.size(); i++) {
      if (x.quality[i] != Quality::valid) e[i] = refused<value_type>(x.quality[i], x.reason(i));
      else if (y.quality[i] != Quality::valid) e[i] = refused<value_type>(y.quality[i], y.reason(i));
      else e[i] = valid(add(x.values[i], y.values[i]), x.note(i).empty() ? y.note(i) : x.note(i));
    }
    return samples_of(std::span<const Evaluation<value_type>>(e));
  }

 private:
  static value_type add(const value_type& u, const value_type& v) {
    if constexpr (std::is_same_v<value_type, double>) return u + v;
    else return value_type{u[0] + v[0], u[1] + v[1], u[2] + v[2]};
  }
  A a_;
  B b_;
  Descriptor d_;
};

template <PointEvaluable A, PointEvaluable B>
Sum<A, B> sum(Descriptor d, A a, B b) { return Sum<A, B>(std::move(d), std::move(a), std::move(b)); }

// ---- pullback -----------------------------------------------------------------------------------------------

// m(q) → Evaluation<Q>: the point of f's base over q's point, or why there is none.
template <PointEvaluable F, class P, class M>
  requires std::is_invocable_r_v<Evaluation<typename F::point_type>, const M&, const Query<P>&>
class Pullback {
 public:
  using value_type = typename F::value_type;
  using point_type = P;

  // d names the new base and its frame; the values, and so the fibre, representation and provenance,
  // are f's.
  Pullback(Descriptor d, F f, M m) : f_(std::move(f)), m_(std::move(m)), d_(std::move(d)) {
    const Descriptor df = described(f_);
    d_.fibre = df.fibre;
    d_.representation = df.representation;
    d_.provenance = df.provenance;
    d_.parts = df.parts;
    if (d_.source.empty()) d_.source = df.source;
  }

  const Descriptor& descriptor() const { return d_; }

  Evaluation<value_type> at(const Query<P>& q) const {
    const Evaluation<typename F::point_type> y = m_(q);
    if (!y.ok()) return like<value_type>(y);
    auto v = f_.at({q.t, y.value});
    if (v.ok() && v.note.empty()) v.note = y.note;
    return v;
  }

 private:
  F f_;
  M m_;
  Descriptor d_;
};

template <class P, PointEvaluable F, class M>
Pullback<F, P, M> pullback(Descriptor d, F f, M m) { return Pullback<F, P, M>(std::move(d), std::move(f), std::move(m)); }

// ---- in another frame ---------------------------------------------------------------------------------------

// The motion from a field's frame to another at an instant: x′ = R x + b.
struct Motion {
  Mat3 R = geometry::identity3();
  Vec3 b{};
};

template <PointEvaluable F, class MotionOf>
  requires std::is_same_v<typename F::point_type, Vec3> && std::is_invocable_r_v<Motion, const MotionOf&, double>
class InFrame {
 public:
  using value_type = typename F::value_type;
  using point_type = Vec3;

  InFrame(F f, MotionOf motion, std::string frame) : f_(std::move(f)), m_(std::move(motion)), d_(described(f_)) {
    d_.frame = std::move(frame);
    d_.capabilities = {};
  }
  const Descriptor& descriptor() const { return d_; }

  Evaluation<value_type> at(const Query<Vec3>& q) const {
    const Motion M = m_(q.t);
    const auto v = f_.at({q.t, back(M, q.x)});
    if constexpr (std::is_same_v<value_type, Vec3>) {
      using namespace geometry::operators;
      if (v.ok()) return valid(Vec3(M.R * v.value), v.note);  // v′ = R v
    }
    return v;
  }
  Evaluation<Derivative<value_type, Vec3>> jacobian(const Query<Vec3>& q) const
    requires SpatiallyDifferentiable<F>
  {
    using namespace geometry::operators;
    const Motion M = m_(q.t);
    const auto j = f_.jacobian({q.t, back(M, q.x)});
    if (!j.ok()) return j;
    if constexpr (std::is_same_v<value_type, double>) return valid(Vec3(M.R * j.value), j.note);  // ∇′φ = R ∇φ
    else return valid(Mat3(M.R * j.value * geometry::transpose(M.R)), j.note);                     // J′ = R J Rᵀ
  }

 private:
  static Vec3 back(const Motion& M, const Vec3& x) {
    using namespace geometry::operators;
    return geometry::transpose(M.R) * (x - M.b);
  }
  F f_;
  MotionOf m_;
  Descriptor d_;
};

template <PointEvaluable F, class MotionOf>
InFrame<F, MotionOf> in_frame(F f, MotionOf motion, std::string frame) { return InFrame<F, MotionOf>(std::move(f), std::move(motion), std::move(frame)); }

// ---- restrict ------------------------------------------------------------------------------------------------

template <PointEvaluable F, class D>
  requires DomainOf<D, typename F::point_type>
class Restricted {
 public:
  using value_type = typename F::value_type;
  using point_type = typename F::point_type;

  Restricted(F f, D domain, std::string about = "") : f_(std::move(f)), domain_(std::move(domain)), d_(described(f_)) {
    if (!about.empty()) d_.base.about = std::move(about);
    d_.capabilities = {};
  }
  const Descriptor& descriptor() const { return d_; }

  Evaluation<value_type> at(const Query<point_type>& q) const {
    const Cover c = domain_.covers(q);
    if (!c.ok()) return refused<value_type>(c.quality, c.reason);
    return f_.at(q);
  }
  Evaluation<Derivative<value_type, point_type>> jacobian(const Query<point_type>& q) const
    requires SpatiallyDifferentiable<F>
  {
    const Cover c = domain_.covers(q);
    if (!c.ok()) return refused<Derivative<value_type, point_type>>(c.quality, c.reason);
    return f_.jacobian(q);
  }

 private:
  F f_;
  D domain_;
  Descriptor d_;
};

template <PointEvaluable F, class D>
Restricted<F, D> restrict(F f, D domain, std::string about = "") { return Restricted<F, D>(std::move(f), std::move(domain), std::move(about)); }

}  // namespace physis::fields
