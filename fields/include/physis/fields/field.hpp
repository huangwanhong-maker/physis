// Fields as sections: what a field is, independently of how it is stored, sampled, solved or drawn
// (decision 0012; docs/theory/fields.md).
//
// A field is a section φ ∈ Γ(E) of a bundle E → B over a base B, with values in a fibre that carries a
// physical dimension. Four things about a field are kept apart, because each can change while the
// others stay:
//
//   the base          where it is defined: a region of space, a surface, a curve, a worldline, a
//                     spacetime region, phase space, time alone, or an abstract complex that is itself
//                     the world (Base);
//   the fibre         what values it takes: a scalar, a vector or covector on named axes, a symmetric
//                     tensor, a k-form, a label in a finite set, a measure, a connection, a state, each
//                     with its unit (Fibre);
//   the representation how it is held for computation: an analytic function, a spectral expansion,
//                     values sampled on a grid with a stated reconstruction, a cochain on a cellulation,
//                     scattered points, a composite of parts, a uniform value (Representation);
//   the provenance    what kind of evidence its values are: physical, archive, model, assimilated or
//                     synthetic (Provenance), per field and per part.
//
// A field answers a query, an instant and a point of its base, with an Evaluation: a value, or one of
// three reasons why there is none, never a zero in its place. *outside*: the point is not in the base
// (the World Magnetic Model below 1 km under the ground, salinity in the air); *missing*: inside, but
// the data that would give the value are absent; *invalid*: inside, but the value is ill-defined there.
// A valid value may carry a note that qualifies it ("read at the sphere above").
//
// What a field can do is a set of capabilities, each a C++20 concept on its representation rather than
// a method every field must have: a cochain pairs with chains and is not evaluated at a point without a
// reconstruction; a Hodge star needs a metric; a spectral expansion samples a sphere faster than point by
// point. capabilities_of<F>() reads them off the type, and the descriptor reports them to whoever asks
// at run time.

#pragma once

#include <array>
#include <concepts>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include <physis/geometry/rotation.hpp>
#include <physis/units/quantity.hpp>

namespace physis::fields {

using geometry::Mat3;
using geometry::Vec3;

// ---- queries and their answers ------------------------------------------------------------------------

// An instant and a point of the base. t is the world's time coordinate (s); a field that does not
// change ignores it. P is the type of the base's points: Vec3 for space (m, on the axes of the field's
// frame), double for a line or for time alone.
template <class P>
struct Query {
  double t = 0;
  P x{};
};

enum class Quality : std::uint8_t { valid, outside, missing, invalid };

constexpr std::string_view name(Quality q) {
  switch (q) {
    case Quality::valid: return "valid";
    case Quality::outside: return "outside";
    case Quality::missing: return "missing";
    case Quality::invalid: return "invalid";
  }
  return "?";
}

// The answer to a query. reason and note are views of static text, so that answers are cheap to make
// in the millions; a caller that needs a reason with numbers in it composes it from the value.
template <class V>
struct Evaluation {
  V value{};
  Quality quality = Quality::valid;
  std::string_view reason{};  // why there is no value (quality ≠ valid)
  std::string_view note{};    // a qualification of a valid value

  constexpr bool ok() const noexcept { return quality == Quality::valid; }
};

template <class V>
constexpr Evaluation<V> valid(V value, std::string_view note = {}) {
  return {std::move(value), Quality::valid, {}, note};
}
template <class V>
constexpr Evaluation<V> refused(Quality quality, std::string_view reason) {
  return {V{}, quality, reason, {}};
}
// An answer of another value type that keeps the quality, reason and note of `of`.
template <class W, class V>
constexpr Evaluation<W> like(const Evaluation<V>& of, W value = W{}) {
  return {std::move(value), of.quality, of.reason, of.note};
}

// ---- the three axes, the representation and the provenance (formalism Definition 3.9) -------------------

enum class BaseKind : std::uint8_t { point, worldline, curve, surface, region, spacetime, phase_space, radiation, time, complex };
enum class FibreKind : std::uint8_t { scalar, vector, covector, symmetric_tensor, form, label, measure, connection, state };
enum class RepresentationKind : std::uint8_t { analytic, spectral, sampled, cochain, points, composite, uniform, observed };
enum class Provenance : std::uint8_t { physical, archive, model, assimilated, synthetic };

constexpr std::string_view name(BaseKind k) {
  constexpr std::array<std::string_view, 10> n{"point", "worldline", "curve", "surface", "region", "spacetime", "phase space", "radiation", "time", "complex"};
  return n[static_cast<std::size_t>(k)];
}
constexpr std::string_view name(FibreKind k) {
  constexpr std::array<std::string_view, 9> n{"scalar", "vector", "covector", "symmetric tensor", "form", "label", "measure", "connection", "state"};
  return n[static_cast<std::size_t>(k)];
}
constexpr std::string_view name(RepresentationKind k) {
  constexpr std::array<std::string_view, 8> n{"analytic", "spectral", "sampled", "cochain", "points", "composite", "uniform", "observed"};
  return n[static_cast<std::size_t>(k)];
}
constexpr std::string_view name(Provenance p) {
  constexpr std::array<std::string_view, 5> n{"physical", "archive", "model", "assimilated", "synthetic"};
  return n[static_cast<std::size_t>(p)];
}

struct Base {
  BaseKind kind = BaseKind::region;
  int dimension = 3;
  std::string about;  // its extent, in words: "the Earth from 1 km below the ellipsoid to 850 km above it"
};

struct Fibre {
  FibreKind kind = FibreKind::scalar;
  units::Dimension dimension{};
  int degree = 0;                   // k of a k-form
  std::string basis;                // the axes of a vector's components: "frame" (the field's), or a named local basis
  std::vector<std::string> labels;  // the finite set of a label field
};

struct Representation {
  RepresentationKind kind = RepresentationKind::analytic;
  std::string about;  // "spherical harmonics to degree 449 on a sphere of 1737.4 km"
};

// ---- capabilities -----------------------------------------------------------------------------------------

enum class Capability : std::uint32_t {
  point = 1u << 0,        // evaluated at a point
  batch = 1u << 1,        // evaluated at many points at once, faster than one by one
  jacobian = 1u << 2,     // its spatial derivative, exact to rounding (not by differences)
  chain = 1u << 3,        // paired with chains (a cochain: integrals over cells)
  coboundary = 1u << 4,   // has d
  hodge = 1u << 5,        // has a Hodge star (a metric and an orientation)
  sphere_grid = 1u << 6,  // samples a sphere faster than point by point (a spectral expansion)
  transfer = 1u << 7,     // carried to another cellulation, preserving a stated structure
  versions = 1u << 8,     // committed versions in time (an evolving field)
};

struct Capabilities {
  std::uint32_t bits = 0;
  constexpr bool has(Capability c) const noexcept { return (bits & static_cast<std::uint32_t>(c)) != 0; }
  constexpr void add(Capability c) noexcept { bits |= static_cast<std::uint32_t>(c); }
};

inline std::vector<std::string_view> names(Capabilities c) {
  constexpr std::array<std::pair<Capability, std::string_view>, 9> all{{{Capability::point, "point"},
                                                                         {Capability::batch, "batch"},
                                                                         {Capability::jacobian, "jacobian"},
                                                                         {Capability::chain, "chain"},
                                                                         {Capability::coboundary, "coboundary"},
                                                                         {Capability::hodge, "hodge"},
                                                                         {Capability::sphere_grid, "sphere grid"},
                                                                         {Capability::transfer, "transfer"},
                                                                         {Capability::versions, "versions"}}};
  std::vector<std::string_view> out;
  for (const auto& [cap, n] : all)
    if (c.has(cap)) out.push_back(n);
  return out;
}

// ---- the descriptor ---------------------------------------------------------------------------------------

struct Part {
  std::string name;
  Provenance provenance = Provenance::model;
};

struct Descriptor {
  std::string name;   // "earth.magnetic"
  std::string title;  // "The geomagnetic field (WMM2025)"
  std::string frame;  // the frame of positions and of vector components ("ITRS", "MOON_ME", "box")
  std::string time;   // the time coordinate of queries ("TT seconds from J2000.0"), or "none"
  Base base;
  Fibre fibre;
  Representation representation;
  Provenance provenance = Provenance::model;
  std::vector<Part> parts;  // a composite's parts, each with its own provenance
  std::string source;       // the model and its reference
  Capabilities capabilities;
};

// ---- the spatial derivative -------------------------------------------------------------------------------

// The type of the derivative of a V-valued field on points P: a scalar's gradient, a vector's Jacobian
// J[i][j] = ∂vᵢ/∂xⱼ.
template <class V, class P>
struct derivative_of;
template <>
struct derivative_of<double, double> { using type = double; };
template <>
struct derivative_of<double, Vec3> { using type = Vec3; };
template <>
struct derivative_of<Vec3, Vec3> { using type = Mat3; };
template <class V, class P>
using Derivative = typename derivative_of<V, P>::type;

inline double divergence(const Mat3& j) { return j[0][0] + j[1][1] + j[2][2]; }
inline Vec3 curl(const Mat3& j) { return {j[2][1] - j[1][2], j[0][2] - j[2][0], j[1][0] - j[0][1]}; }

// ---- the capabilities as concepts -------------------------------------------------------------------------

template <class F>
concept Field = requires(const F& f) {
  typename F::value_type;
  typename F::point_type;
  { f.descriptor() } -> std::convertible_to<const Descriptor&>;
};

template <class F>
concept PointEvaluable = Field<F> && requires(const F& f, const Query<typename F::point_type>& q) {
  { f.at(q) } -> std::same_as<Evaluation<typename F::value_type>>;
};

template <class F>
concept BatchEvaluable = PointEvaluable<F> && requires(const F& f, std::span<const Query<typename F::point_type>> q,
                                                       std::span<Evaluation<typename F::value_type>> out) {
  f.at_many(q, out);
};

template <class F>
concept SpatiallyDifferentiable = PointEvaluable<F> && requires(const F& f, const Query<typename F::point_type>& q) {
  { f.jacobian(q) } -> std::same_as<Evaluation<Derivative<typename F::value_type, typename F::point_type>>>;
};

template <class F>
concept Versioned = Field<F> && requires(const F& f) {
  { f.versions() } -> std::convertible_to<std::size_t>;
};

// The capabilities a type has, read off the concepts it satisfies; a type may add others it states
// itself (a cochain's pairing, a spectral expansion's sphere grid) through a static member
// `extra_capabilities`.
template <class F>
constexpr Capabilities capabilities_of() {
  Capabilities c;
  if constexpr (PointEvaluable<F>) c.add(Capability::point);
  if constexpr (BatchEvaluable<F>) c.add(Capability::batch);
  if constexpr (SpatiallyDifferentiable<F>) c.add(Capability::jacobian);
  if constexpr (Versioned<F>) c.add(Capability::versions);
  if constexpr (requires { F::extra_capabilities; }) c.bits |= F::extra_capabilities;
  return c;
}

// The field's descriptor, with its capabilities filled in from its type (a type-erased field already
// carries its field's).
template <Field F>
Descriptor described(const F& f) {
  Descriptor d = f.descriptor();
  if constexpr (!requires { F::type_erased; }) d.capabilities.bits |= capabilities_of<F>().bits;
  return d;
}

// ---- evaluation in bulk -----------------------------------------------------------------------------------

template <PointEvaluable F>
void evaluate_many(const F& f, std::span<const Query<typename F::point_type>> queries, std::span<Evaluation<typename F::value_type>> out) {
  if (queries.size() != out.size()) throw std::invalid_argument("physis::fields::evaluate_many: one answer for each query");
  if constexpr (BatchEvaluable<F>) {
    f.at_many(queries, out);
  } else {
    for (std::size_t i = 0; i < queries.size(); i++) out[i] = f.at(queries[i]);
  }
}

template <PointEvaluable F>
std::vector<Evaluation<typename F::value_type>> evaluate_many(const F& f, std::span<const Query<typename F::point_type>> queries) {
  std::vector<Evaluation<typename F::value_type>> out(queries.size());
  evaluate_many(f, queries, std::span<Evaluation<typename F::value_type>>(out));
  return out;
}

}  // namespace physis::fields
