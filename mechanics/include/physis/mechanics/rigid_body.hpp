// Rigid bodies: their mass properties, composed from parts, and their motion by a Lie group variational
// integrator (decision 0013; docs/theory/discrete-mechanics.md §4).
//
// A rigid body's configuration is the position x of its centre of mass and its attitude R ∈ SO(3), the
// active rotation that carries body coordinates to space coordinates (physis geometry's convention:
// v_space = R v_body). Its momenta are the linear momentum p (space) and the angular momentum Π = J Ω in
// body axes, with J the inertia tensor about the centre and Ω the body angular velocity, Ṙ = R Ω^.
//
// Mass properties. A body is a union of parts (solid boxes, point masses, or the cells of a material
// domain): masses add, the centre is the mass-weighted mean, and each part's inertia about its own centre
// is carried to the body's by the parallel axis theorem (Steiner), J = Σ J_i + m_i (|d_i|² I − d_i d_iᵀ),
// d_i = c_i − c. A solid box of sides a, b, c has J = (m/12) diag(b² + c², a² + c², a² + b²). The
// principal axes diagonalise J, by calx's symmetric eigensolver.
//
// Motion. The discrete Lagrangian of a rigid body (Lee, Leok and McClamroch 2007, Lie group variational
// integrators for the full body problem, Comput. Methods Appl. Mech. Engrg. 196, 2907–2924, §3)
//   L_d = (m/2h)|x_{k+1} − x_k|² + (1/h) tr[(I − F_k) J_d] − (h/2)(U_k + U_{k+1}),   F_k = R_kᵀ R_{k+1},
// with the nonstandard inertia J_d = ½ tr(J) I − J, has the discrete equations
//   h Π_k^ + (h²/2) M_k^ = F_k J_d − J_d F_kᵀ,   R_{k+1} = R_k F_k,
//   Π_{k+1} = F_kᵀ (Π_k + (h/2) M_k) + (h/2) M_{k+1},
//   p_{k+½} = p_k + (h/2) f_k,   x_{k+1} = x_k + (h/m) p_{k+½},   p_{k+1} = p_{k+½} + (h/2) f_{k+1},
// where f is the force (space) and M the moment about the centre (body axes). R stays on SO(3) exactly,
// since it is only ever multiplied by rotations; the map is symplectic; and when the wrench comes from a
// potential invariant under rotations of space, U(Q x, Q R) = U(x, R), the total angular momentum
// x × p + R Π is conserved exactly (the discrete Noether theorem; Lee et al. 2007, §3).
//
// The implicit equation for F_k is solved in Cayley's parameter f, F = (I − f^)⁻¹(I + f^) = I +
// c (f^ + f^ f^), c = 2/(1 + |f|²), the rotation by 2 atan|f| about f, whose quaternion is (1, f)
// normalised. For a symmetric A, a^ A + A a^ = ((tr A) I − A) a)^ and [a^, b^] = (a × b)^, so
//   vee(F J_d − J_d Fᵀ) = c (J f + f × J f),
// and the equation is G(f) = c (J f + f × J f) − h Π' = 0 with Π' = Π_k + (h/2) M_k: three smooth rational
// equations, solved by Newton's method on its exact Jacobian (calx jets) from f = (h/2) J⁻¹ Π'. The step
// must turn the body by less than π, which any useful step does.
//
// A wrench from a potential. potential_wrench takes U(x, R) written over any number type and returns
// f = −∂U/∂x and M = −∂U/∂η at η = 0 for the variation R(η) = R (I + η^), both by jets: the force and the
// moment of one function, as exact conservation needs. A wrench not from a potential (a motor, contact,
// drag) is any function of (t, x, R) returning a Wrench.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <calx/differentiation/jet.hpp>
#include <calx/linalg/dense.hpp>
#include <calx/linalg/eigen.hpp>
#include <calx/linalg/small.hpp>
#include <physis/geometry/rotation.hpp>
#include <physis/mechanics/variational.hpp>

namespace physis::mechanics {

using geometry::Mat3;
using geometry::Rotation;
using geometry::Vec3;

// ---- mass properties ----------------------------------------------------------------------------------------

struct MassProperties {
  double mass = 0;
  Vec3 centre{};   // the centre of mass, body coordinates
  Mat3 inertia{};  // about the centre, body axes: J = ∫ ρ (|r|² I − r rᵀ) dV
};

// A solid box [lower, upper] of uniform density (kg/m³).
inline MassProperties solid_box(const Vec3& lower, const Vec3& upper, double density) {
  const double a = upper[0] - lower[0], b = upper[1] - lower[1], c = upper[2] - lower[2];
  if (!(a > 0 && b > 0 && c > 0)) throw std::invalid_argument("physis::mechanics::solid_box: upper above lower on every axis");
  if (!(density >= 0)) throw std::invalid_argument("physis::mechanics::solid_box: a density of zero or more");
  MassProperties P;
  P.mass = density * a * b * c;
  P.centre = {0.5 * (lower[0] + upper[0]), 0.5 * (lower[1] + upper[1]), 0.5 * (lower[2] + upper[2])};
  const double k = P.mass / 12;
  P.inertia = {{{k * (b * b + c * c), 0, 0}, {0, k * (a * a + c * c), 0}, {0, 0, k * (a * a + b * b)}}};
  return P;
}

// A point mass (kg) at a point of the body.
inline MassProperties point_mass(const Vec3& at, double mass) {
  if (!(mass >= 0)) throw std::invalid_argument("physis::mechanics::point_mass: a mass of zero or more");
  MassProperties P;
  P.mass = mass;
  P.centre = at;
  return P;
}

// The body made of parts: masses add, the centre is their mean, and the inertias are carried to it.
inline MassProperties combine(const std::vector<MassProperties>& parts) {
  MassProperties B;
  for (const auto& p : parts) {
    B.mass += p.mass;
    for (int i = 0; i < 3; i++) B.centre[i] += p.mass * p.centre[i];
  }
  if (!(B.mass > 0)) throw std::invalid_argument("physis::mechanics::combine: a body of positive mass");
  for (int i = 0; i < 3; i++) B.centre[i] /= B.mass;
  for (const auto& p : parts) {
    const Vec3 d{p.centre[0] - B.centre[0], p.centre[1] - B.centre[1], p.centre[2] - B.centre[2]};
    const double d2 = geometry::dot(d, d);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) B.inertia[i][j] += p.inertia[i][j] + p.mass * ((i == j ? d2 : 0) - d[i] * d[j]);
  }
  return B;
}

// J = A diag(moments) Aᵀ, moments ascending, A a rotation whose columns are the principal axes.
struct PrincipalAxes {
  Vec3 moments{};
  Rotation axes;
};
inline PrincipalAxes principal_axes(const Mat3& J) {
  calx::linalg::Matrix A(3);
  for (std::size_t i = 0; i < 3; i++)
    for (std::size_t j = 0; j < 3; j++) A(i, j) = J[i][j];
  const auto e = calx::linalg::symmetric_eigen(A);
  Mat3 V{};
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) V[i][j] = e.vectors(std::size_t(i), std::size_t(j));
  if (geometry::determinant(V) < 0)
    for (int i = 0; i < 3; i++) V[i][2] = -V[i][2];
  return {{e.values[0], e.values[1], e.values[2]}, Rotation::from_matrix(V)};
}

// ---- the state and its wrench ---------------------------------------------------------------------------------

struct RigidBodyState {
  double t = 0;
  Vec3 x{};    // the centre of mass, space frame
  Vec3 p{};    // the linear momentum, space frame
  Rotation R;  // the attitude: body → space
  Vec3 Pi{};   // the angular momentum about the centre, body axes
};

struct Wrench {
  Vec3 force{};   // on the centre of mass, space frame
  Vec3 moment{};  // about the centre of mass, body axes
};

// The place in space of a point of the body (body coordinates, the same origin as MassProperties.centre).
inline Vec3 place(const RigidBodyState& s, const MassProperties& body, const Vec3& point) {
  using namespace geometry::operators;
  return s.x + s.R.apply(point - body.centre);
}

// The body angular velocity Ω = J⁻¹ Π.
inline Vec3 angular_velocity(const MassProperties& body, const Vec3& Pi) {
  return calx::linalg::solve<3>(body.inertia, Pi);
}

// The kinetic energy |p|²/2m + ½ Π · J⁻¹ Π.
inline double kinetic_energy(const MassProperties& body, const RigidBodyState& s) {
  return geometry::dot(s.p, s.p) / (2 * body.mass) + 0.5 * geometry::dot(s.Pi, angular_velocity(body, s.Pi));
}

// The angular momentum about the origin of space, x × p + R Π.
inline Vec3 angular_momentum(const RigidBodyState& s) {
  using namespace geometry::operators;
  return geometry::cross(s.x, s.p) + s.R.apply(s.Pi);
}

namespace rigid_detail {
template <class S>
using V3 = std::array<S, 3>;
template <class S>
using M3 = std::array<std::array<S, 3>, 3>;
}  // namespace rigid_detail

// The force −∂U/∂x and moment −∂U/∂η (R(η) = R(I + η^)) of a potential U(x, R) written over any number
// type: U(const std::array<S, 3>& x, const std::array<std::array<S, 3>, 3>& R) → S.
template <class Potential>
Wrench potential_wrench(const Potential& U, const Vec3& x, const Rotation& R) {
  using D = calx::differentiation::Dual<3>;
  const Mat3 r = R.matrix();
  Wrench w;
  {
    const auto xs = calx::differentiation::variables<3, 1>(x);
    rigid_detail::M3<D> rs;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) rs[i][j] = D(r[i][j]);
    const D u = U(xs, rs);
    for (int i = 0; i < 3; i++) w.force[i] = -u.d[i];
  }
  {
    rigid_detail::V3<D> xs{D(x[0]), D(x[1]), D(x[2])};
    const auto eta = calx::differentiation::variables<3, 1>(Vec3{0, 0, 0});
    // η^ = [[0, −η₃, η₂], [η₃, 0, −η₁], [−η₂, η₁, 0]]; R(η) = R + R η^.
    const rigid_detail::M3<D> hat{{{D(0), -eta[2], eta[1]}, {eta[2], D(0), -eta[0]}, {-eta[1], eta[0], D(0)}}};
    rigid_detail::M3<D> rs;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        D v(r[i][j]);
        for (int k = 0; k < 3; k++) v += r[i][k] * hat[k][j];
        rs[i][j] = v;
      }
    const D u = U(xs, rs);
    for (int i = 0; i < 3; i++) w.moment[i] = -u.d[i];
  }
  return w;
}

// The value of a potential U(x, R) at a configuration.
template <class Potential>
double potential_energy(const Potential& U, const Vec3& x, const Rotation& R) {
  const Mat3 r = R.matrix();
  return U(x, r);
}

// ---- the Lie group variational integrator ----------------------------------------------------------------------

class RigidBodyIntegrator {
 public:
  RigidBodyIntegrator(MassProperties body, double h, double tolerance = 1e-14, std::size_t max_iterations = 30)
      : body_(body), h_(h), tol_(tolerance), max_(max_iterations) {
    if (!(h > 0) || !std::isfinite(h)) throw std::invalid_argument("physis::mechanics::RigidBodyIntegrator: a positive step");
    if (!(body.mass > 0)) throw std::invalid_argument("physis::mechanics::RigidBodyIntegrator: a body of positive mass");
    const auto P = principal_axes(body.inertia);
    if (!(P.moments[0] > 0)) throw std::invalid_argument("physis::mechanics::RigidBodyIntegrator: an inertia tensor that is positive definite");
  }

  const MassProperties& body() const { return body_; }
  double step_size() const { return h_; }

  // One step under wrench(t, x, R) → Wrench.
  template <class WrenchFunction>
  Convergence step(RigidBodyState& s, WrenchFunction&& wrench) const {
    using namespace geometry::operators;
    const double h = h_;
    const Wrench w0 = wrench(s.t, s.x, s.R);
    const Vec3 p_half = s.p + (0.5 * h) * w0.force;
    const Vec3 Pi_half = s.Pi + (0.5 * h) * w0.moment;
    Convergence c;
    const Vec3 f = cayley_step(Pi_half, c);
    const Rotation F = Rotation::from_quaternion(1, f[0], f[1], f[2]);
    const Rotation R1 = s.R * F;
    s.R = Rotation::from_quaternion(R1.w(), R1.x(), R1.y(), R1.z());  // the product, its norm's rounding removed
    s.x = s.x + (h / body_.mass) * p_half;
    s.t += h;
    const Vec3 Pi_turned = F.inverse().apply(Pi_half);
    const Wrench w1 = wrench(s.t, s.x, s.R);
    s.p = p_half + (0.5 * h) * w1.force;
    s.Pi = Pi_turned + (0.5 * h) * w1.moment;
    return c;
  }

  // The free rotation of one step: Cayley's parameter f with c(f)(J f + f × J f) = h Π'.
  Vec3 cayley_step(const Vec3& Pi, Convergence& c) const {
    using D = calx::differentiation::Dual<3>;
    const Mat3& J = body_.inertia;
    const double h = h_;
    Vec3 f = angular_velocity(body_, Pi);
    for (auto& x : f) x *= 0.5 * h;
    double scale = 0;
    for (int i = 0; i < 3; i++) scale = std::max(scale, std::fabs(f[i]));
    for (;;) {
      if (c.iterations >= max_) throw std::runtime_error("physis::mechanics::RigidBodyIntegrator: the rotation of a step did not converge");
      const auto z = calx::differentiation::variables<3, 1>(f);
      std::array<D, 3> Jf;
      for (int i = 0; i < 3; i++) Jf[i] = z[0] * J[i][0] + z[1] * J[i][1] + z[2] * J[i][2];
      const D k = 2.0 / (1.0 + z[0] * z[0] + z[1] * z[1] + z[2] * z[2]);
      const std::array<D, 3> fxJf{z[1] * Jf[2] - z[2] * Jf[1], z[2] * Jf[0] - z[0] * Jf[2], z[0] * Jf[1] - z[1] * Jf[0]};
      calx::linalg::SmallMatrix<3> G{};
      calx::linalg::SmallVector<3> r{};
      for (int i = 0; i < 3; i++) {
        const D g = k * (Jf[i] + fxJf[i]);
        r[i] = h * Pi[i] - g.v;
        for (int j = 0; j < 3; j++) G[i][j] = g.d[j];
      }
      const auto df = calx::linalg::solve<3>(G, r);
      c.iterations++;
      double change = 0;
      for (int i = 0; i < 3; i++) f[i] += df[i], change = std::max(change, std::fabs(df[i])), scale = std::max(scale, std::fabs(f[i]));
      c.change = change;
      if (change <= tol_ * scale) return f;
    }
  }

 private:
  MassProperties body_;
  double h_, tol_;
  std::size_t max_;
};

}  // namespace physis::mechanics
