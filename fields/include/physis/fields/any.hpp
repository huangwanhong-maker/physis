// Fields of any representation behind one type, for registries and for the boundary with JSON
// (decision 0012; docs/theory/fields.md §8).
//
// AnyField<V, P> holds a point-evaluable field of value type V on points P whatever its representation,
// and keeps its descriptor (with the capabilities read off its type), its batches, its exact derivative
// when it has one, and its own way of sampling a sphere. The call through the erased type is made once
// for a batch or a sampling, not once for each point, so that a registry of a world's fields costs no
// more than the fields themselves (numerical kernels keep the static types).

#pragma once

#include <memory>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

#include <physis/fields/field.hpp>
#include <physis/fields/sampling.hpp>

namespace physis::fields {

template <class V, class P = Vec3>
class AnyField {
 public:
  using value_type = V;
  using point_type = P;
  static constexpr bool type_erased = true;  // its capabilities are its field's, not its own type's

  AnyField() = default;

  template <PointEvaluable F>
    requires(!std::is_same_v<std::decay_t<F>, AnyField> && std::is_same_v<typename F::value_type, V> && std::is_same_v<typename F::point_type, P>)
  AnyField(F f) : m_(std::make_shared<const Holder<F>>(std::move(f))) {}

  explicit operator bool() const { return static_cast<bool>(m_); }
  const Descriptor& descriptor() const { return m_->descriptor(); }

  Evaluation<V> at(const Query<P>& q) const { return m_->at(q); }
  void at_many(std::span<const Query<P>> q, std::span<Evaluation<V>> out) const { m_->at_many(q, out); }
  bool differentiable() const { return descriptor().capabilities.has(Capability::jacobian); }
  Evaluation<Derivative<V, P>> jacobian(const Query<P>& q) const { return m_->jacobian(q); }

  Samples<V> sample(const SphereGrid& g, double t) const
    requires std::is_same_v<P, Vec3>
  {
    return m_->sample_sphere(g, t);
  }
  template <Sampling G>
  Samples<V> sample(const G& g, double t) const
    requires std::is_same_v<P, Vec3>
  {
    const auto q = queries_of(g, t);
    std::vector<Evaluation<V>> a(q.size());
    m_->at_many(std::span<const Query<P>>(q), std::span<Evaluation<V>>(a));
    return samples_of(std::span<const Evaluation<V>>(a));
  }

 private:
  struct Model {
    virtual ~Model() = default;
    virtual const Descriptor& descriptor() const = 0;
    virtual Evaluation<V> at(const Query<P>&) const = 0;
    virtual void at_many(std::span<const Query<P>>, std::span<Evaluation<V>>) const = 0;
    virtual Evaluation<Derivative<V, P>> jacobian(const Query<P>&) const = 0;
    virtual Samples<V> sample_sphere(const SphereGrid&, double) const = 0;
  };

  template <class F>
  struct Holder final : Model {
    F f;
    Descriptor d;
    explicit Holder(F x) : f(std::move(x)), d(described(f)) {}
    const Descriptor& descriptor() const override { return d; }
    Evaluation<V> at(const Query<P>& q) const override { return f.at(q); }
    void at_many(std::span<const Query<P>> q, std::span<Evaluation<V>> out) const override { evaluate_many(f, q, out); }
    Evaluation<Derivative<V, P>> jacobian(const Query<P>& q) const override {
      if constexpr (SpatiallyDifferentiable<F>) return f.jacobian(q);
      else return refused<Derivative<V, P>>(Quality::invalid, "the field has no exact spatial derivative");
    }
    Samples<V> sample_sphere(const SphereGrid& g, double t) const override {
      if constexpr (std::is_same_v<P, Vec3>) return fields::sample(f, g, t);
      else return {};
    }
  };

  std::shared_ptr<const Model> m_;
};

using ScalarField = AnyField<double, Vec3>;
using VectorField = AnyField<Vec3, Vec3>;

}  // namespace physis::fields
